#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <random>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/SkmerPartialOrder.hpp>
#include <algorithms/VirtualSkmer.hpp>

using kuint = uint64_t;
using namespace km::experiment;

// ----------------------- compare_keys (4-way) ---------------------------

TEST(SkmerPartialOrder, compare_keys_less_greater)
{
    std::vector<uint8_t> a {0, 1, 2};
    std::vector<uint8_t> b {0, 2, 0};
    EXPECT_EQ(compare_keys(a, b), RealOrder::Less);    // differ at index 1: 1 < 2
    EXPECT_EQ(compare_keys(b, a), RealOrder::Greater);
}

TEST(SkmerPartialOrder, compare_keys_equal)
{
    std::vector<uint8_t> a {3, 0, 1};
    std::vector<uint8_t> b {3, 0, 1};
    EXPECT_EQ(compare_keys(a, b), RealOrder::Equal);
}

TEST(SkmerPartialOrder, compare_keys_prefix_is_incomparable)
{
    // A proper prefix must NOT be ordered — it is the absent-nucleotide ambiguity.
    std::vector<uint8_t> shortk {2, 1};
    std::vector<uint8_t> longk  {2, 1, 0, 3};
    EXPECT_EQ(compare_keys(shortk, longk), RealOrder::Incomparable);
    EXPECT_EQ(compare_keys(longk, shortk), RealOrder::Incomparable);
}

TEST(SkmerPartialOrder, compare_keys_divergence_beats_length)
{
    // Real divergence at index 1 decides even though one is shorter overall.
    std::vector<uint8_t> a {1, 0};
    std::vector<uint8_t> b {1, 2, 2, 2};
    EXPECT_EQ(compare_keys(a, b), RealOrder::Less); // 0 < 2 at index 1
}

// ----------------------- union-find ---------------------------

TEST(SkmerPartialOrder, union_find_basic)
{
    UnionFind uf;
    const uint32_t a = uf.make(), b = uf.make(), c = uf.make(), d = uf.make();
    EXPECT_NE(uf.find(a), uf.find(b));
    uf.unite(a, b);
    EXPECT_EQ(uf.find(a), uf.find(b));
    uf.unite(b, c);
    EXPECT_EQ(uf.find(a), uf.find(c)); // transitive closure
    EXPECT_NE(uf.find(a), uf.find(d)); // untouched stays separate
}

// ----------------------- partial reconciliation ---------------------------

// Enumerate the super-k-mers of one sequence (same path as FileSkmerator, per record).
static void enumerate_into(km::SkmerManipulator<kuint>& manip, std::string seq,
                           std::vector<km::Skmer<kuint>>& out)
{
    if (seq.length() < manip.k) return;
    km::SeqSkmerator<kuint> rator{manip, seq};
    for (const km::Skmer<kuint>& sk : rator) out.push_back(sk);
}

// Run both reconciliations on the same enumeration of `seqs` and return the two
// compacted lists plus the partial path's group_of.
static void run_both(uint64_t k, uint64_t m, const std::vector<std::string>& seqs,
                     std::vector<km::Skmer<kuint>>& prod_list,
                     std::vector<km::Skmer<kuint>>& part_list,
                     std::vector<uint32_t>& group_of)
{
    km::SkmerManipulator<kuint> manip{k, m};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const std::string& s : seqs) enumerate_into(manip, s, enumeration);

    km::sortedlist::SortedVirtualSkmerList<kuint> prod(k, m);
    prod.generate_sorted_list_from_enumeration(enumeration);
    prod_list = prod.get_list();

    km::sortedlist::SortedVirtualSkmerList<kuint> part(k, m);
    part.generate_partial_list_from_enumeration(enumeration, group_of);
    part_list = part.get_list();
}

static uint64_t count_real_groups(const std::vector<uint32_t>& group_of)
{
    std::unordered_map<uint32_t, uint64_t> counts;
    for (uint32_t g : group_of) counts[g]++;
    uint64_t n = 0;
    for (const auto& kv : counts) if (kv.second > 1) ++n;
    return n;
}

// The partial path reproduces production's compacted list exactly (order included):
// the colinear-chaining fusions and the emitted order are untouched.
TEST(SkmerPartialOrder, partial_reconciliation_matches_production)
{
    std::vector<km::Skmer<kuint>> a, b;
    std::vector<uint32_t> group_of;
    run_both(7, 3, {"ACGTACGTACGTACGT", "GGGGACGTACGTACGTGGGG",
                    "TTACGTACGTACTT", "ACGTACGTAC"}, a, b, group_of);

    ASSERT_FALSE(b.empty());
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i)
        EXPECT_TRUE(a[i] == b[i]) << "compacted skmer differs from production at " << i;
    ASSERT_EQ(group_of.size(), b.size());
}

// Same minimizer reached in different (non-fusable) contexts with different flank
// lengths is orderable only via the absent-nucleotide mask: the merge must refuse to
// order them and emit a group. This compact 3-sequence input is a known such case.
TEST(SkmerPartialOrder, partial_reconciliation_forms_group)
{
    std::vector<km::Skmer<kuint>> a, b;
    std::vector<uint32_t> group_of;
    run_both(7, 3, {"ATCACACCCAACCTTCAAATGCCG",
                    "GCCCTAACGCCCTAATCCTGCGCT",
                    "GGGGTTGCAGCGACCAGATG"}, a, b, group_of);

    // fusions still identical to production
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i)
        EXPECT_TRUE(a[i] == b[i]) << "compacted skmer differs from production at " << i;

    // and at least one mask-only orderable cluster is surfaced
    EXPECT_GE(count_real_groups(group_of), 1u);
}

// ----------------------- virtual nucleotide fill ---------------------------

static const std::vector<std::string> kFillSeqs {
    "ATCACACCCAACCTTCAAATGCCG",
    "GCCCTAACGCCCTAATCCTGCGCT",
    "GGGGTTGCAGCGACCAGATG"
};

// Filling writes only the absent slots: the bounds, the real (present) nucleotides, and the
// per-column order among the real k-mers are all preserved.
TEST(SkmerPartialOrder, fill_preserves_real_content_and_order)
{
    const uint64_t k = 7, m = 3;
    km::SkmerManipulator<kuint> manip{k, m};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const std::string& s : kFillSeqs) enumerate_into(manip, s, enumeration);

    km::sortedlist::SortedVirtualSkmerList<kuint> part(k, m);
    std::vector<uint32_t> group_of;
    part.generate_partial_list_from_enumeration(enumeration, group_of);

    const std::vector<km::Skmer<kuint>> before = part.get_list();
    std::vector<km::Skmer<kuint>> after = before;
    std::vector<uint32_t> fp, fs;
    km::experiment::fill_virtual_nucleotides<kuint>(manip, after, fp, fs);

    ASSERT_EQ(before.size(), after.size());
    for (size_t i = 0; i < after.size(); ++i) {
        EXPECT_EQ(before[i].m_pref_size, after[i].m_pref_size);
        EXPECT_EQ(before[i].m_suff_size, after[i].m_suff_size);
        km::Skmer<kuint> remasked = after[i];
        manip.mask_absent_nucleotides(remasked);
        EXPECT_TRUE(remasked == before[i]) << "real content changed at " << i;
    }
    size_t viol = 0;
    for (uint64_t c = 0; c <= k - m; ++c) {
        const km::Skmer<kuint>* prev = nullptr;
        for (const km::Skmer<kuint>& s : after) {
            if (!manip.has_valid_kmer(s, c)) continue;
            if (prev && manip.kmer_compare(*prev, s, c) > 0) ++viol;
            prev = &s;
        }
    }
    EXPECT_EQ(viol, 0u); // per-column order among real k-mers untouched
}

// The per-column virtual fills render lowercase.
TEST(SkmerPartialOrder, fill_renders_lowercase)
{
    const uint64_t k = 7, m = 3;
    km::SkmerManipulator<kuint> manip{k, m};
    std::vector<km::Skmer<kuint>> enumeration;
    std::mt19937 rng(7);
    static const char nt[] = {'A', 'C', 'G', 'T'};
    for (int s = 0; s < 10; ++s) {
        std::string seq;
        for (int j = 0; j < 120; ++j) seq += nt[rng() & 3u];
        enumerate_into(manip, seq, enumeration);
    }

    km::sortedlist::SortedVirtualSkmerList<kuint> part(k, m);
    std::vector<uint32_t> group_of;
    part.generate_partial_list_from_enumeration(enumeration, group_of);
    std::vector<km::Skmer<kuint>> skmers = part.get_list();
    std::vector<uint32_t> fp, fs;
    km::experiment::fill_virtual_nucleotides<kuint>(manip, skmers, fp, fs);

    bool found = false;
    for (size_t i = 0; i < skmers.size() && !found; ++i) {
        if (fp[i] == 0 && fs[i] == 0) continue;
        const std::string d = km::experiment::decode_with_dots(manip, skmers[i], fp[i], fs[i]);
        if (d.find_first_of("actg") != std::string::npos) found = true; // lowercase fill rendered
    }
    EXPECT_TRUE(found);
}

// On a larger random input: the per-column fill never mis-places a usable (valid OR fully
// filled) skmer (viol_usable == 0), while reducing the number of incomparable (truncated,
// unfilled) (skmer, column) pairs that would force a linear scan.
TEST(SkmerPartialOrder, fill_per_column_order_on_random)
{
    const uint64_t k = 7, m = 3;
    const int64_t flank = static_cast<int64_t>(k) - static_cast<int64_t>(m);
    km::SkmerManipulator<kuint> manip{k, m};
    std::mt19937 rng(12345);
    static const char nt[] = {'A', 'C', 'G', 'T'};
    std::vector<km::Skmer<kuint>> enumeration;
    for (int s = 0; s < 25; ++s) {
        std::string seq;
        for (int j = 0; j < 200; ++j) seq += nt[rng() & 3u];
        enumerate_into(manip, seq, enumeration);
    }
    ASSERT_FALSE(enumeration.empty());

    km::sortedlist::SortedVirtualSkmerList<kuint> part(k, m);
    std::vector<uint32_t> group_of;
    part.generate_partial_list_from_enumeration(enumeration, group_of);
    std::vector<km::Skmer<kuint>> skmers = part.get_list();

    // probe depth of each list position (node depth in the binary-search tree over [0, n-1])
    std::vector<int> depth(skmers.size(), 0);
    {
        struct R { int64_t lo, hi; int d; };
        std::vector<R> st{ { 0, static_cast<int64_t>(skmers.size()) - 1, 0 } };
        while (!st.empty()) {
            const R r = st.back(); st.pop_back();
            if (r.lo > r.hi) continue;
            const int64_t mid = (r.lo + r.hi) / 2;
            depth[mid] = r.d;
            st.push_back({ r.lo, mid - 1, r.d + 1 });
            st.push_back({ mid + 1, r.hi, r.d + 1 });
        }
    }
    std::vector<uint32_t> fp(skmers.size(), 0), fs(skmers.size(), 0);
    auto usable_at = [&](size_t i, uint64_t c) {
        if (manip.has_valid_kmer(skmers[i], c)) return true;
        const int64_t pref = skmers[i].m_pref_size, suff = skmers[i].m_suff_size;
        for (int64_t p = static_cast<int64_t>(c); p < flank - pref; ++p) if (!((fp[i] >> p) & 1u)) return false;
        for (int64_t p = flank - static_cast<int64_t>(c); p < flank - suff; ++p) if (!((fs[i] >> p) & 1u)) return false;
        return true;
    };
    auto measure = [&](const char* tag) {
        size_t viol_usable = 0, incomparable = 0, inc_d2 = 0, inc_d4 = 0;
        for (uint64_t c = 0; c <= k - m; ++c) {
            const km::Skmer<kuint>* prev = nullptr;
            for (size_t i = 0; i < skmers.size(); ++i) {
                if (!usable_at(i, c)) {
                    ++incomparable;
                    if (depth[i] <= 2) ++inc_d2;     // top 7 probe positions
                    if (depth[i] <= 4) ++inc_d4;     // top 31 probe positions
                    continue;
                }
                if (prev && manip.kmer_compare(*prev, skmers[i], c) > 0) ++viol_usable;
                prev = &skmers[i];
            }
        }
        std::cerr << tag << " viol_usable=" << viol_usable << " incomparable=" << incomparable
                  << "  (shallow probes incomparable: depth<=2=" << inc_d2 << ", depth<=4=" << inc_d4 << ")\n";
        return viol_usable;
    };
    measure("[production]");
    km::experiment::fill_virtual_nucleotides<kuint>(manip, skmers, fp, fs);
    EXPECT_EQ(measure("[after fill]"), 0u); // no usable skmer is ever mis-placed

    // Of the shallow-probe (depth<=2) incomparable left, how many COULD fit between their valid
    // anchors (recoverable) vs are genuinely out-of-order at that column (irreducible)?
    auto valid_at = [&](size_t i, uint64_t c) {
        return static_cast<int64_t>(skmers[i].m_pref_size) >= flank - static_cast<int64_t>(c)
            && static_cast<int64_t>(skmers[i].m_suff_size) >= static_cast<int64_t>(c);
    };
    size_t sh_inc = 0, sh_placeable = 0;
    for (uint64_t c = 0; c <= k - m; ++c)
        for (size_t i = 0; i < skmers.size(); ++i) {
            if (depth[i] > 2 || usable_at(i, c)) continue;
            ++sh_inc;
            int64_t A = -1, B = -1;
            for (int64_t j = static_cast<int64_t>(i) - 1; j >= 0; --j) if (valid_at(j, c)) { A = j; break; }
            for (size_t j = i + 1; j < skmers.size(); ++j) if (valid_at(j, c)) { B = static_cast<int64_t>(j); break; }
            if (A < 0 || B < 0) continue;
            const km::Skmer<kuint> save = skmers[i];                 // try filling i with A's values
            const int64_t pref = skmers[i].m_pref_size, suff = skmers[i].m_suff_size;
            for (int64_t p = static_cast<int64_t>(c); p < flank - pref; ++p)
                km::experiment::set_nibble(skmers[i], static_cast<uint64_t>(4 * p), km::experiment::nibble_at(skmers[A], static_cast<uint64_t>(4 * p)));
            for (int64_t p = flank - static_cast<int64_t>(c); p < flank - suff; ++p)
                km::experiment::set_nibble(skmers[i], static_cast<uint64_t>(4 * p + 2), km::experiment::nibble_at(skmers[A], static_cast<uint64_t>(4 * p + 2)));
            if (manip.kmer_compare(skmers[A], skmers[i], c) <= 0 && manip.kmer_compare(skmers[i], skmers[B], c) <= 0) ++sh_placeable;
            skmers[i] = save;
        }
    std::cerr << "shallow(depth<=2) incomparable=" << sh_inc << "  recoverable(fits anchors)=" << sh_placeable
              << "  irreducible=" << (sh_inc - sh_placeable) << "\n";
}
