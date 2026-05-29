// Micro-benchmark: existing query (query_skmer_batch, with find_closest_valid_skmer
// linear scans) vs the scan-free hole-aware query (query_skmer_batch_substrate) on the
// sentinel-filled substrate. Gated by SKLIB_BENCH (skipped otherwise) and meant to run
// from a Release build. Prints wall-clock + ns/super-k-mer to stderr; also asserts the
// two queries agree (correctness alongside the timing).

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {
std::string random_seq(size_t len, uint32_t seed) {
    static const char alpha[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (size_t i = 0; i < len; i++) s[i] = alpha[d(rng)];
    return s;
}
template <typename kuint>
std::vector<km::Skmer<kuint>> enumerate(km::SkmerManipulator<kuint>& manip, const std::string& seq) {
    std::string s = seq;
    km::SeqSkmerator<kuint> rator{manip, s};
    std::vector<km::Skmer<kuint>> out;
    for (const km::Skmer<kuint>& x : rator) out.push_back(x);
    return out;
}
double time_ms(const std::function<void()>& f) {
    auto t0 = std::chrono::high_resolution_clock::now();
    f();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}
} // namespace

TEST(SentinelBench, QueryThroughput) {
    if (std::getenv("SKLIB_BENCH") == nullptr)
        GTEST_SKIP() << "set SKLIB_BENCH=1 to run (Release build recommended)";

    using kuint = uint64_t;
    const uint64_t k = 21, m = 11;
    km::SkmerManipulator<kuint> manip{k, m};

    // Size sweep: at what genome size does hole-inclusive per-column monotonicity break?
    {
        using kpair = typename km::Skmer<kuint>::pair;
        const std::vector<kpair> kmask = manip.get_k_mask();
        for (size_t sz : {size_t{1000}, size_t{10000}, size_t{100000}, size_t{1000000}}) {
            std::string g = random_seq(sz, 4242);
            auto en = enumerate(manip, g);
            km::sortedlist::SortedVirtualSkmerList<kuint> L(k, m);
            L.generate_sorted_list_from_enumeration(en);
            L.fill_absent_sentinel();
            const auto& v = L.get_list();
            long va = 0, vv = 0;
            for (uint64_t col = 0; col <= k - m; col++) {
                kpair pa{}, pv{}; bool ha = false, hv = false;
                for (size_t i = 0; i < v.size(); i++) {
                    kpair key{v[i].m_pair & kmask[col]};
                    if (ha && key < pa) va++;
                    pa = key; ha = true;
                    if (manip.has_valid_kmer(v[i], col)) { if (hv && key < pv) vv++; pv = key; hv = true; }
                }
            }
            std::cerr << "[bench] sweep size=" << sz << " entries=" << v.size()
                      << " holes+valid_viol=" << va << " valid_only_viol=" << vv << std::endl;
        }
    }

    const size_t G = 2'000'000;

    std::string genome = random_seq(G, 20260528);
    std::vector<km::Skmer<kuint>> enumeration = enumerate(manip, genome);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);
    list.fill_absent_sentinel();
    std::cerr << "[bench] k=" << k << " m=" << m << " genome=" << G
              << " list_entries=" << list.size()
              << " query_skmers=" << enumeration.size() << std::endl;

    // Column monotonicity at scale (holes included): does the navigability invariant
    // that held on small test sequences still hold on a 2 Mb genome?
    {
        using kpair = typename km::Skmer<kuint>::pair;
        const std::vector<kpair> kmask = manip.get_k_mask();
        const auto& v = list.get_list();
        long viol_all = 0, viol_valid = 0;
        for (uint64_t col = 0; col <= k - m; col++) {
            kpair pa{}, pv{}; bool ha = false, hv = false;
            for (size_t i = 0; i < v.size(); i++) {
                kpair key{v[i].m_pair & kmask[col]};
                if (ha && key < pa) viol_all++;
                pa = key; ha = true;
                if (manip.has_valid_kmer(v[i], col)) {
                    if (hv && key < pv) viol_valid++;
                    pv = key; hv = true;
                }
            }
        }
        std::cerr << "[bench] column monotonicity at scale: holes+valid_viol=" << viol_all
                  << " valid_only_viol=" << viol_valid << std::endl;
    }

    std::vector<km::Skmer<kuint>> present = enumeration;
    std::vector<km::Skmer<kuint>> absent = enumerate(manip, random_seq(G, 777));

    for (auto* wl : {&present, &absent}) {
        const char* name = (wl == &present) ? "present" : "absent ";
        std::vector<std::vector<uint8_t>> r_old, r_new;
        double ms_old = time_ms([&]{ r_old = list.query_skmer_batch(*wl); });
        double ms_new = time_ms([&]{ r_new = list.query_skmer_batch_substrate(*wl); });
        // Count disagreements rather than assert, to quantify any scale effect.
        long mism = 0, total = 0;
        for (size_t i = 0; i < r_old.size(); i++) {
            total += static_cast<long>(r_old[i].size());
            for (size_t j = 0; j < r_old[i].size(); j++)
                if (r_old[i][j] != r_new[i][j]) mism++;
        }
        const double per_old = ms_old * 1e6 / static_cast<double>(wl->size());
        const double per_new = ms_new * 1e6 / static_cast<double>(wl->size());
        std::cerr << "[bench] " << name << " n=" << wl->size()
                  << " | existing " << ms_old << " ms (" << per_old << " ns/skmer)"
                  << " | substrate " << ms_new << " ms (" << per_new << " ns/skmer)"
                  << " | speedup x" << (ms_old / ms_new)
                  << " | mismatches=" << mism << "/" << total << std::endl;
    }

    // --- Rescue experiment: does fill + RE-SORT by (m_pair, pref, suff) yield a
    //     fully navigable substrate (and stay query-correct)? ---
    {
        using kpair = typename km::Skmer<kuint>::pair;
        std::vector<km::Skmer<kuint>> sorted = list.get_list(); // already sentinel-filled
        std::sort(sorted.begin(), sorted.end(),
                  [](const km::Skmer<kuint>& a, const km::Skmer<kuint>& b) {
                      if (a.m_pair == b.m_pair) {
                          if (a.m_pref_size == b.m_pref_size) return a.m_suff_size < b.m_suff_size;
                          return a.m_pref_size < b.m_pref_size;
                      }
                      return a.m_pair < b.m_pair;
                  });
        km::sortedlist::SortedVirtualSkmerList<kuint> rs(k, m);
        rs.add_list(sorted);

        const std::vector<kpair> kmask = manip.get_k_mask();
        long va = 0;
        for (uint64_t col = 0; col <= k - m; col++) {
            kpair p{}; bool h = false;
            for (size_t i = 0; i < sorted.size(); i++) {
                kpair key{sorted[i].m_pair & kmask[col]};
                if (h && key < p) va++;
                p = key; h = true;
            }
        }
        auto gt  = list.query_skmer_batch(present);            // ground truth (correct)
        auto rso = rs.query_skmer_batch(present);
        auto rsn = rs.query_skmer_batch_substrate(present);
        long m1 = 0, m2 = 0;
        for (size_t i = 0; i < gt.size(); i++)
            for (size_t j = 0; j < gt[i].size(); j++) {
                if (rso[i][j] != gt[i][j]) m1++;
                if (rsn[i][j] != gt[i][j]) m2++;
            }
        std::cerr << "[bench] fill+resort: holes+valid_viol=" << va
                  << " | existing-query mismatch_vs_gt=" << m1
                  << " | substrate-query mismatch_vs_gt=" << m2 << std::endl;
    }
}

// PRODUCTION benchmark: lib's fill_absent_interpolated + query_skmer_batch_bounded vs the
// existing query_skmer_batch, on a sentinel/interpolated-filled list. Asserts identical
// results (exactness) and reports throughput.
TEST(SentinelBench, BoundedQueryThroughput) {
    if (std::getenv("SKLIB_BENCH") == nullptr) GTEST_SKIP() << "set SKLIB_BENCH=1 (Release build)";
    using kuint = uint64_t;
    const uint64_t k = 21, m = 11;
    const size_t G = 2000000;
    // random + repetitive (long hole-runs: where the existing find_closest scans far)
    std::vector<std::pair<std::string, std::string>> genomes;
    genomes.emplace_back("random ", random_seq(G, 20260529));
    { std::string motif = random_seq(60, 1), s; s.reserve(G); while (s.size() < G) s += motif; genomes.emplace_back("repeat60", s); }
    { std::string motif = random_seq(300, 2), s; s.reserve(G); while (s.size() < G) s += motif; genomes.emplace_back("repeat300", s); }
    for (auto& [glabel, genome] : genomes) {
        km::SkmerManipulator<kuint> manip{k, m};
        auto en = enumerate(manip, genome);
        km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
        list.generate_sorted_list_from_enumeration(en);
        list.fill_absent_interpolated();
        std::cerr << "[boundedbench] " << glabel << " entries=" << list.size()
                  << " window_D=" << list.query_window() << std::endl;

        std::vector<km::Skmer<kuint>> present = en;
        std::vector<km::Skmer<kuint>> absent = enumerate(manip, random_seq(G, 13));
        for (auto* wl : {&present, &absent}) {
            const char* name = (wl == &present) ? "present" : "absent ";
            std::vector<std::vector<uint8_t>> r_old, r_new;
            double ms_old = time_ms([&]{ r_old = list.query_skmer_batch(*wl); });
            double ms_new = time_ms([&]{ r_new = list.query_skmer_batch_bounded(*wl); });
            ASSERT_EQ(r_old, r_new) << "bounded != existing on " << name << " (" << glabel << ")";
            std::cerr << "[boundedbench] " << name << " n=" << wl->size()
                      << " | existing " << ms_old << " ms (" << ms_old * 1e6 / wl->size() << " ns/skmer)"
                      << " | bounded " << ms_new << " ms (" << ms_new * 1e6 / wl->size() << " ns/skmer)"
                      << " | speedup x" << (ms_old / ms_new) << std::endl;
        }
    }
}
