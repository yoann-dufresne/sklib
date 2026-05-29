// Experiment for Yoann's interpolation idea: instead of filling absent slots with a
// constant (pure-min), fill them so the list stays NON-DECREASING in m_pair, with holes
// interpolated between their bracketing neighbours. Question: does an order-preserving
// fill make the per-column key monotone (holes included) — i.e. make the scan-free
// hole-aware query correct — at genome scale?
//
// We test the two EXTREME order-preserving fills, carry-min (lowest monotone completion)
// and carry-max (highest), plus their midpoint. They bracket every order-preserving fill
// (incl. linear interpolation): if all of them still leave per-column violations / query
// false-negatives, no order-preserving fill can make the substrate navigable; if one is
// clean, the idea works.
//
// Gated by SKLIB_BENCH; run from a Release build.

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <random>
#include <iostream>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {
using u64 = uint64_t;

std::string random_seq(size_t len, uint32_t seed) {
    static const char a[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (size_t i = 0; i < len; i++) s[i] = a[d(rng)];
    return s;
}
template <typename kuint>
std::vector<km::Skmer<kuint>> enumerate(km::SkmerManipulator<kuint>& m, const std::string& seq) {
    std::string s = seq;
    km::SeqSkmerator<kuint> r{m, s};
    std::vector<km::Skmer<kuint>> o;
    for (const km::Skmer<kuint>& x : r) o.push_back(x);
    return o;
}

inline bool ge128(u64 alo, u64 ahi, u64 blo, u64 bhi) { return ahi != bhi ? ahi > bhi : alo >= blo; }
inline bool le128(u64 alo, u64 ahi, u64 blo, u64 bhi) { return ahi != bhi ? ahi < bhi : alo <= blo; }
inline bool lt128(u64 alo, u64 ahi, u64 blo, u64 bhi) { return ahi != bhi ? ahi < bhi : alo <  blo; }

struct V128 { u64 lo, hi; };

// Smallest value with present bits fixed and absent bits ⊆ amask that is >= target.
// Start from H (all absent set), clear absent bits MSB->LSB while staying >= target.
V128 smallest_ge(u64 plo, u64 phi, u64 amlo, u64 amhi, u64 tlo, u64 thi) {
    u64 vlo = plo | amlo, vhi = phi | amhi;
    for (int b = 63; b >= 0; --b) { u64 bit = u64(1) << b; if (amhi & bit) { u64 nhi = vhi & ~bit; if (ge128(vlo, nhi, tlo, thi)) vhi = nhi; } }
    for (int b = 63; b >= 0; --b) { u64 bit = u64(1) << b; if (amlo & bit) { u64 nlo = vlo & ~bit; if (ge128(nlo, vhi, tlo, thi)) vlo = nlo; } }
    return {vlo, vhi};
}
// Largest value with present bits fixed and absent bits ⊆ amask that is <= target.
// Start from L (all absent clear), set absent bits MSB->LSB while staying <= target.
V128 largest_le(u64 plo, u64 phi, u64 amlo, u64 amhi, u64 tlo, u64 thi) {
    u64 vlo = plo & ~amlo, vhi = phi & ~amhi;
    for (int b = 63; b >= 0; --b) { u64 bit = u64(1) << b; if (amhi & bit) { u64 nhi = vhi | bit; if (le128(vlo, nhi, tlo, thi)) vhi = nhi; } }
    for (int b = 63; b >= 0; --b) { u64 bit = u64(1) << b; if (amlo & bit) { u64 nlo = vlo | bit; if (le128(nlo, vhi, tlo, thi)) vlo = nlo; } }
    return {vlo, vhi};
}

template <typename kuint>
long per_column_violations(const std::vector<km::Skmer<kuint>>& v,
                           km::SkmerManipulator<kuint>& manip, u64 k, u64 m) {
    using kpair = typename km::Skmer<kuint>::pair;
    const std::vector<kpair> kmask = manip.get_k_mask();
    long viol = 0;
    for (u64 c = 0; c <= k - m; c++) {
        kpair prev{}; bool have = false;
        for (size_t i = 0; i < v.size(); i++) {
            kpair key{v[i].m_pair & kmask[c]};
            if (have && key < prev) viol++;
            prev = key; have = true;
        }
    }
    return viol;
}

// Count mismatches between two batched query results.
long mismatches(const std::vector<std::vector<uint8_t>>& a,
                const std::vector<std::vector<uint8_t>>& b) {
    long mm = 0;
    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < a[i].size(); j++)
            if (a[i][j] != b[i][j]) mm++;
    return mm;
}
} // namespace

TEST(InterpExperiment, OrderPreservingFills) {
    if (std::getenv("SKLIB_BENCH") == nullptr)
        GTEST_SKIP() << "set SKLIB_BENCH=1 (Release build)";
    using kuint = uint64_t;
    const u64 k = 21, m = 11;
    for (size_t G : {size_t{200000}, size_t{1000000}}) {
        km::SkmerManipulator<kuint> manip{k, m};
        std::string genome = random_seq(G, 5);
        auto en = enumerate(manip, genome);
        km::sortedlist::SortedVirtualSkmerList<kuint> L(k, m);
        L.generate_sorted_list_from_enumeration(en);
        std::vector<km::Skmer<kuint>> base = L.get_list();   // unfilled (absent = 0b11)
        const size_t n = base.size();

        // Ground truth: existing (correct) query on the unfilled list.
        std::vector<km::Skmer<kuint>> present = en;
        std::vector<km::Skmer<kuint>> absent  = enumerate(manip, random_seq(G, 777));
        auto gt_present = L.query_skmer_batch(present);
        auto gt_absent  = L.query_skmer_batch(absent);

        // Precompute present bits and absent masks.
        std::vector<u64> plo(n), phi(n), amlo(n), amhi(n);
        for (size_t i = 0; i < n; i++) {
            auto am = manip.absent_slot_mask(base[i].m_pref_size, base[i].m_suff_size);
            amlo[i] = am.m_value[0]; amhi[i] = am.m_value[1];
            plo[i] = base[i].m_pair.m_value[0] & ~amlo[i];
            phi[i] = base[i].m_pair.m_value[1] & ~amhi[i];
        }

        // carry-min (forward): smallest monotone completion.
        std::vector<V128> lo(n);
        long stuck_min = 0;
        { u64 plo0 = 0, phi0 = 0;
          for (size_t i = 0; i < n; i++) {
              u64 tlo = (i == 0) ? 0 : lo[i-1].lo, thi = (i == 0) ? 0 : lo[i-1].hi;
              lo[i] = smallest_ge(plo[i], phi[i], amlo[i], amhi[i], tlo, thi);
              if (i && lt128(lo[i].lo, lo[i].hi, lo[i-1].lo, lo[i-1].hi)) stuck_min++;
              (void)plo0; (void)phi0;
          } }
        // carry-max (backward): largest monotone completion.
        std::vector<V128> hi(n);
        long stuck_max = 0;
        for (size_t ii = 0; ii < n; ii++) {
            size_t i = n - 1 - ii;
            u64 tlo = (ii == 0) ? ~u64(0) : hi[i+1].lo, thi = (ii == 0) ? ~u64(0) : hi[i+1].hi;
            hi[i] = largest_le(plo[i], phi[i], amlo[i], amhi[i], tlo, thi);
            if (ii && !le128(hi[i].lo, hi[i].hi, hi[i+1].lo, hi[i+1].hi)) stuck_max++;
        }
        // midpoint, kept monotone: smallest achievable >= max(prev, avg(lo[i], hi[i])).
        std::vector<V128> mid(n);
        for (size_t i = 0; i < n; i++) {
            // 128-bit average of lo[i] and hi[i]
            u64 sumlo = lo[i].lo + hi[i].lo;
            u64 carry = (sumlo < lo[i].lo) ? 1 : 0;
            u64 sumhi = lo[i].hi + hi[i].hi + carry;
            u64 avglo = (sumlo >> 1) | (sumhi << 63);
            u64 avghi = sumhi >> 1;
            u64 tlo = avglo, thi = avghi;
            if (i && ge128(mid[i-1].lo, mid[i-1].hi, tlo, thi)) { tlo = mid[i-1].lo; thi = mid[i-1].hi; }
            mid[i] = smallest_ge(plo[i], phi[i], amlo[i], amhi[i], tlo, thi);
        }

        auto eval = [&](const char* name, const std::vector<V128>& Vv, long stuck) {
            std::vector<km::Skmer<kuint>> filled = base;
            for (size_t i = 0; i < n; i++) { filled[i].m_pair.m_value[0] = Vv[i].lo; filled[i].m_pair.m_value[1] = Vv[i].hi; }
            long mono = 0;
            for (size_t i = 1; i < n; i++) if (lt128(Vv[i].lo, Vv[i].hi, Vv[i-1].lo, Vv[i-1].hi)) mono++;
            long colv = per_column_violations(filled, manip, k, m);
            km::sortedlist::SortedVirtualSkmerList<kuint> F(k, m);
            F.add_list(filled);
            long fn_present = mismatches(gt_present, F.query_skmer_batch_substrate(present));
            long fn_absent  = mismatches(gt_absent,  F.query_skmer_batch_substrate(absent));
            std::cerr << "[interp] G=" << G << " " << name
                      << " stuck=" << stuck << " mpair_inversions=" << mono
                      << " per_col_viol=" << colv
                      << " substrate_query_mismatch present=" << fn_present
                      << " absent=" << fn_absent << std::endl;
        };
        std::cerr << "[interp] G=" << G << " entries=" << n << std::endl;
        eval("carry-min", lo, stuck_min);
        eval("carry-max", hi, stuck_max);
        eval("midpoint ", mid, 0);
    }
}
