// Tests for the absent-slot sentinel fill (SortedVirtualSkmerList::fill_absent_sentinel).
//
// The fill clears each stored entry's absent (low/peripheral) flank bits from the
// build-time 0b11 padding down to 0 (the minimal completion). It was meant to make the
// list per-column monotone WITH holes included, so a binary search could navigate through
// entries lacking a k-mer at the column without find_closest_valid_skmer's linear scan.
//
// That goal is NOT achieved beyond trivially small inputs (see SubstrateHoleMonotonicity
// FailsAtScale and docs/sentinel_substrate.md). What these tests pin down is what is
// actually true:
//   * Construction-only safety (every k, m): query results are IDENTICAL with/without the
//     fill (the current query reads k-mers only at entries valid at the column, where the
//     bits are masked out, and skips holes); the list keeps its size; each entry's m_pair
//     is only lowered (present bits / sizes unchanged); self-query stays complete; round-trips.
//   * The valid-entry per-column monotonicity the current query relies on holds for all k
//     and at genome scale.
//   * The hole-INCLUSIVE monotonicity (the navigability the fill aimed for) fails at scale.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>
#include <random>
#include <filesystem>
#include <unistd.h>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {

template <typename kuint>
std::vector<km::Skmer<kuint>> enumerate(km::SkmerManipulator<kuint>& manip, std::string seq) {
    km::SeqSkmerator<kuint> rator{manip, seq};
    std::vector<km::Skmer<kuint>> out;
    for (const km::Skmer<kuint>& s : rator) out.push_back(s);
    return out;
}

template <typename kuint>
km::sortedlist::SortedVirtualSkmerList<kuint>
build(uint64_t k, uint64_t m, const std::vector<km::Skmer<kuint>>& enumeration, bool fill) {
    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);
    if (fill) list.fill_absent_sentinel();
    return list;
}

std::string random_seq(size_t len, uint32_t seed) {
    static const char alpha[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (size_t i = 0; i < len; i++) s[i] = alpha[d(rng)];
    return s;
}

std::string repeat_motif(const std::string& motif, size_t times) {
    std::string s;
    s.reserve(motif.size() * times);
    for (size_t i = 0; i < times; i++) s += motif;
    return s;
}

struct Config { uint64_t k, m; };

// Spread incl. degenerate tiny k-m for the universal safety tests.
const std::vector<Config> kAllConfigs = {
    {5, 2}, {8, 4}, {15, 7}, {21, 11}, {23, 11}, {31, 15},
};

// A spread of sequence shapes incl. low-complexity stressors for the sentinel.
std::vector<std::string> sample_sequences(uint64_t k) {
    const size_t L = 8 * k + 50;
    return {
        random_seq(L, 1234),
        random_seq(L, 99),
        random_seq(2 * L, 7),
        repeat_motif("ACGT", L / 4),
        repeat_motif("AACAG", L / 5),
        repeat_motif("GC", L / 2),                 // RC-palindrome-prone motif
        std::string(L, 'A'),                       // poly-A
        std::string(L, 'G'),                       // poly-G
        std::string(L / 2, 'A') + random_seq(L / 2, 55),
    };
}

// Count column-key monotonicity violations along the list. With include_holes=false only
// entries valid at the column are considered (the current-query invariant); with true,
// every entry is considered (the navigability the fill aimed for).
template <typename kuint>
long column_violations(const km::sortedlist::SortedVirtualSkmerList<kuint>& list,
                       km::SkmerManipulator<kuint>& manip, uint64_t k, uint64_t m,
                       bool include_holes) {
    using kpair = typename km::Skmer<kuint>::pair;
    const std::vector<kpair> kmask = manip.get_k_mask();
    const auto& v = list.get_list();
    long viol = 0;
    for (uint64_t c = 0; c <= k - m; c++) {
        kpair prev{};
        bool have = false;
        for (size_t i = 0; i < v.size(); i++) {
            if (!include_holes && !manip.has_valid_kmer(v[i], c)) continue;
            kpair key{v[i].m_pair & kmask[c]};
            if (have && key < prev) viol++;
            prev = key;
            have = true;
        }
    }
    return viol;
}

} // namespace

// ---- Universal: the fill must not change query answers (construction-only safety) ----
TEST(SentinelSubstrate, FillPreservesQueryResults) {
    using kuint = uint64_t;
    for (const Config& c : kAllConfigs) {
        km::SkmerManipulator<kuint> manip{c.k, c.m};
        for (const std::string& ref : sample_sequences(c.k)) {
            if (ref.size() < c.k) continue;
            std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, ref);
            if (enumeration.empty()) continue;

            auto unfilled = build<kuint>(c.k, c.m, enumeration, false);
            auto filled   = build<kuint>(c.k, c.m, enumeration, true);

            std::vector<std::string> queries = {ref, random_seq(5 * c.k, 4242),
                                                random_seq(5 * c.k, 888),
                                                repeat_motif("AC", 3 * c.k)};
            std::vector<km::Skmer<kuint>> qskmers;
            for (const std::string& q : queries) {
                if (q.size() < c.k) continue;
                std::vector<km::Skmer<kuint>> e = enumerate(manip, q);
                qskmers.insert(qskmers.end(), e.begin(), e.end());
            }
            if (qskmers.empty()) continue;

            ASSERT_EQ(unfilled.query_skmer_batch(qskmers), filled.query_skmer_batch(qskmers))
                << "fill changed query results for k=" << c.k << " m=" << c.m;
        }
    }
}

// ---- Universal: same size, present bits / sizes unchanged, m_pair only lowered ----
TEST(SentinelSubstrate, FillKeepsSizeAndOnlyLowers) {
    using kuint = uint64_t;
    for (const Config& c : kAllConfigs) {
        km::SkmerManipulator<kuint> manip{c.k, c.m};
        for (const std::string& ref : sample_sequences(c.k)) {
            if (ref.size() < c.k) continue;
            std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, ref);
            if (enumeration.empty()) continue;

            auto unfilled = build<kuint>(c.k, c.m, enumeration, false);
            auto filled   = build<kuint>(c.k, c.m, enumeration, true);
            const auto& a = unfilled.get_list();
            const auto& b = filled.get_list();
            ASSERT_EQ(a.size(), b.size());
            for (size_t i = 0; i < a.size(); i++) {
                EXPECT_EQ(a[i].m_pref_size, b[i].m_pref_size);
                EXPECT_EQ(a[i].m_suff_size, b[i].m_suff_size);
                EXPECT_FALSE(a[i].m_pair < b[i].m_pair)
                    << "fill raised entry " << i << " (k=" << c.k << " m=" << c.m << ")";
            }
        }
    }
}

// ---- Universal: self-query stays complete on the filled (shipped) substrate ----
TEST(SentinelSubstrate, SelfQueryAllPresentOnFilledList) {
    using kuint = uint64_t;
    for (const Config& c : kAllConfigs) {
        km::SkmerManipulator<kuint> manip{c.k, c.m};
        for (const std::string& ref : sample_sequences(c.k)) {
            if (ref.size() < c.k) continue;
            std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, ref);
            if (enumeration.empty()) continue;
            auto filled = build<kuint>(c.k, c.m, enumeration, true);
            for (const std::vector<uint8_t>& row : filled.query_skmer_batch(enumeration))
                for (uint8_t val : row)
                    EXPECT_NE(val, 0) << "self-query miss k=" << c.k << " m=" << c.m;
        }
    }
}

// ---- Universal: save/load round-trip on the filled list answers identically ----
TEST(SentinelSubstrate, RoundTripSaveLoadFilled) {
    using kuint = uint64_t;
    namespace fs = std::filesystem;
    const Config c{21, 11};
    km::SkmerManipulator<kuint> manip{c.k, c.m};
    std::string ref = random_seq(2000, 31337);
    std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, ref);
    ASSERT_FALSE(enumeration.empty());

    auto filled = build<kuint>(c.k, c.m, enumeration, true);
    fs::path tmp = fs::temp_directory_path() /
                   ("sentinel_rt_" + std::to_string(::getpid()) + ".sskm");
    km::sortedlist::VirtualSkmerSerializer<kuint>::save(filled, tmp.string());
    auto loaded = km::sortedlist::VirtualSkmerSerializer<kuint>::load(tmp.string());

    std::vector<km::Skmer<kuint>> qskmers = enumerate(manip, ref);
    std::vector<km::Skmer<kuint>> extra = enumerate(manip, random_seq(1500, 2024));
    qskmers.insert(qskmers.end(), extra.begin(), extra.end());
    EXPECT_EQ(filled.query_skmer_batch(qskmers), loaded.query_skmer_batch(qskmers));

    std::error_code ec;
    fs::remove(tmp, ec);
}

// ---- Current-query invariant (all k): per-column key non-decreasing over VALID entries ----
TEST(SentinelSubstrate, ValidEntriesPerColumnMonotone) {
    using kuint = uint64_t;
    for (const Config& c : kAllConfigs) {
        km::SkmerManipulator<kuint> manip{c.k, c.m};
        for (const std::string& ref : sample_sequences(c.k)) {
            if (ref.size() < c.k) continue;
            std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, ref);
            if (enumeration.empty()) continue;
            auto filled = build<kuint>(c.k, c.m, enumeration, true);
            EXPECT_EQ(column_violations(filled, manip, c.k, c.m, /*include_holes=*/false), 0)
                << "valid-entry column order broken k=" << c.k << " m=" << c.m;
        }
    }
}

// Same invariant at GENOME SCALE.
TEST(SentinelSubstrate, ValidEntriesPerColumnMonotoneAtScale) {
    using kuint = uint64_t;
    const Config c{21, 11};
    km::SkmerManipulator<kuint> manip{c.k, c.m};
    std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, random_seq(200000, 5));
    auto filled = build<kuint>(c.k, c.m, enumeration, true);
    EXPECT_EQ(column_violations(filled, manip, c.k, c.m, /*include_holes=*/false), 0)
        << "valid-entry column order broken at scale";
}

// IMPORTANT NEGATIVE RESULT (documents the limitation; see docs/sentinel_substrate.md):
// the sentinel fill does NOT make the list per-column monotone once HOLES are included,
// beyond trivially small inputs — so it does NOT enable a scan-free hole-aware binary
// search at genome scale. Holes (entries lacking a k-mer at a column) are placed by the
// columns where they ARE valid; at a column where they are absent, their high-order
// content is already out of order, and no absent-bit fill can fix high-order bits.
// (If this ever starts FAILING because the substrate became navigable, the idea was
// rescued — update the query and the docs accordingly.)
TEST(SentinelSubstrate, SubstrateHoleMonotonicityFailsAtScale) {
    using kuint = uint64_t;
    const Config c{21, 11};
    km::SkmerManipulator<kuint> manip{c.k, c.m};
    std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, random_seq(200000, 5));
    auto filled = build<kuint>(c.k, c.m, enumeration, true);
    EXPECT_GT(column_violations(filled, manip, c.k, c.m, /*include_holes=*/true), 0)
        << "holes unexpectedly monotone at scale — the sentinel idea may be rescuable";
}

// The hole-aware query (query_skmer_substrate) is exact WHEN the substrate is monotone,
// i.e. on tiny inputs below the break threshold. It is NOT correct at genome scale (see
// the SKLIB_BENCH micro-benchmark, which reports ~5% false negatives on a 2 Mb genome).
TEST(SentinelSubstrate, HoleAwareQueryExactOnTinyMonotoneInput) {
    using kuint = uint64_t;
    const Config c{21, 11};
    km::SkmerManipulator<kuint> manip{c.k, c.m};
    std::string ref = random_seq(1000, 9); // below the monotonicity break threshold
    std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, ref);
    auto filled = build<kuint>(c.k, c.m, enumeration, true);
    ASSERT_EQ(column_violations(filled, manip, c.k, c.m, /*include_holes=*/true), 0)
        << "precondition: tiny input must be hole-monotone";

    std::vector<km::Skmer<kuint>> q = enumeration;
    std::vector<km::Skmer<kuint>> e = enumerate(manip, random_seq(2000, 71));
    q.insert(q.end(), e.begin(), e.end());
    EXPECT_EQ(filled.query_skmer_batch(q), filled.query_skmer_batch_substrate(q))
        << "hole-aware query must match the existing query while the substrate is monotone";
}
