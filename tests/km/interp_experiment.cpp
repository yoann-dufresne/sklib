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

// Column-AWARE fill: for each column c, interpolate each hole's MISSING window slots
// between its valid-at-c neighbours (so column c's key lands between them); a slot
// constrained by several columns gets the average of their desired values. This targets
// per-column order directly, not the global m_pair order.
template <typename kuint>
std::vector<km::Skmer<kuint>> column_aware_fill(const std::vector<km::Skmer<kuint>>& base,
                                                km::SkmerManipulator<kuint>& manip,
                                                u64 k, u64 m) {
    const size_t n = base.size();
    const u64 flank = k - m;
    const size_t S = 2 * flank;                       // slot index: prefix s -> s, suffix s -> flank+s
    auto slotval = [&](const km::Skmer<kuint>& e, u64 off) -> unsigned {
        return off < 64 ? unsigned((e.m_pair.m_value[0] >> off) & 3ULL)
                        : unsigned((e.m_pair.m_value[1] >> (off - 64)) & 3ULL);
    };
    std::vector<float> sum(n * S, 0.0f), cnt(n * S, 0.0f);

    for (u64 c = 0; c <= flank; c++) {
        long il = -1;                                  // index of last valid-at-c entry
        // first pass forward; for runs between valid brackets, interpolate
        std::vector<long> valids;
        for (size_t idx = 0; idx < n; idx++)
            if (manip.has_valid_kmer(base[idx], c)) valids.push_back((long)idx);
        if (valids.empty()) continue;
        auto accum = [&](long h, long L, long R) {
            const u64 pref = base[h].m_pref_size, suff = base[h].m_suff_size;
            double f = (R > L) ? double(h - L) / double(R - L) : 0.0;
            // missing prefix slots in window [c, flank-1]: s in [c, flank-pref-1]
            for (u64 s = c; s + pref < flank; s++) {
                double lv = (L >= 0) ? slotval(base[L], 4 * s) : 0;
                double rv = (R >= 0) ? slotval(base[R], 4 * s) : lv;
                double d = (L >= 0) ? lv + f * (rv - lv) : rv;
                sum[h * S + s] += float(d); cnt[h * S + s] += 1;
            }
            // missing suffix slots in window: s in [flank-c, flank-suff-1]
            for (u64 s = flank - c; s + suff < flank; s++) {
                double lv = (L >= 0) ? slotval(base[L], 4 * s + 2) : 0;
                double rv = (R >= 0) ? slotval(base[R], 4 * s + 2) : lv;
                double d = (L >= 0) ? lv + f * (rv - lv) : rv;
                sum[h * S + flank + s] += float(d); cnt[h * S + flank + s] += 1;
            }
        };
        long prevValid = -1;
        for (long v : valids) {
            for (long h = prevValid + 1; h < v; h++) accum(h, prevValid, v);
            prevValid = v;
        }
        for (long h = prevValid + 1; h < (long)n; h++) accum(h, prevValid, -1);
        il = prevValid; (void)il;
    }

    std::vector<km::Skmer<kuint>> filled = base;
    for (size_t i = 0; i < n; i++) {
        const u64 pref = base[i].m_pref_size, suff = base[i].m_suff_size;
        for (u64 s = 0; s + pref < flank; s++) {       // absent prefix slot s
            unsigned v = cnt[i * S + s] > 0 ? unsigned(sum[i * S + s] / cnt[i * S + s] + 0.5f) : 0;
            v &= 3; u64 off = 4 * s;
            if (off < 64) { filled[i].m_pair.m_value[0] &= ~(3ULL << off); filled[i].m_pair.m_value[0] |= u64(v) << off; }
            else { u64 o = off - 64; filled[i].m_pair.m_value[1] &= ~(3ULL << o); filled[i].m_pair.m_value[1] |= u64(v) << o; }
        }
        for (u64 s = 0; s + suff < flank; s++) {       // absent suffix slot s
            unsigned v = cnt[i * S + flank + s] > 0 ? unsigned(sum[i * S + flank + s] / cnt[i * S + flank + s] + 0.5f) : 0;
            v &= 3; u64 off = 4 * s + 2;
            if (off < 64) { filled[i].m_pair.m_value[0] &= ~(3ULL << off); filled[i].m_pair.m_value[0] |= u64(v) << off; }
            else { u64 o = off - 64; filled[i].m_pair.m_value[1] &= ~(3ULL << o); filled[i].m_pair.m_value[1] |= u64(v) << o; }
        }
    }
    return filled;
}

// Max per-column displacement = how far each entry is from its sorted position in any
// column key (the window W a bounded query would need; 0 == fully per-column monotone).
template <typename kuint>
long max_col_displacement(const std::vector<km::Skmer<kuint>>& v,
                          km::SkmerManipulator<kuint>& manip, u64 k, u64 m) {
    using kpair = typename km::Skmer<kuint>::pair;
    const std::vector<kpair> kmask = manip.get_k_mask();
    const size_t n = v.size();
    std::vector<uint32_t> ord(n);
    long md = 0;
    for (u64 c = 0; c <= k - m; c++) {
        for (uint32_t z = 0; z < (uint32_t)n; z++) ord[z] = z;
        std::stable_sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
            return (v[a].m_pair & kmask[c]) < (v[b].m_pair & kmask[c]); });
        for (uint32_t r = 0; r < (uint32_t)n; r++) { long d = (long)ord[r] - (long)r; if (d < 0) d = -d; if (d > md) md = d; }
    }
    return md;
}

// Yoann's exact procedure (procedure_remplissage.txt): hierarchical interpolation anchored
// on COMPLETE entries + pinned boundaries (lst[0]->0b00, lst[n-1]->0b11). Process interleaved
// nucleotide positions low->high significance; in each window between consecutive anchors,
// resolve an entry once its highest-significance missing nucleotide reaches the current
// position, by index-interpolating between the current anchors (writing only its missing
// slots), then promote it to an anchor (and tighten the left endpoint).
template <typename kuint>
std::vector<km::Skmer<kuint>> interpolated_fill(const std::vector<km::Skmer<kuint>>& base,
                                                km::SkmerManipulator<kuint>& manip, u64 k, u64 m) {
    using u128 = unsigned __int128;
    const size_t n = base.size();
    const u64 flank = k - m;
    std::vector<km::Skmer<kuint>> lst = base;
    auto getv = [&](size_t i) -> u128 {
        return (u128(lst[i].m_pair.m_value[1]) << 64) | u128(lst[i].m_pair.m_value[0]);
    };
    auto write2 = [&](size_t x, u64 off, unsigned val) {
        if (off < 64) { lst[x].m_pair.m_value[0] = (lst[x].m_pair.m_value[0] & ~(u64(3) << off)) | (u64(val) << off); }
        else { u64 o = off - 64; lst[x].m_pair.m_value[1] = (lst[x].m_pair.m_value[1] & ~(u64(3) << o)) | (u64(val) << o); }
    };
    auto set_missing = [&](size_t x, u128 target) {
        const u64 pref = lst[x].m_pref_size, suff = lst[x].m_suff_size;
        for (u64 s = 0; s + pref < flank; s++) { u64 off = 4 * s;     write2(x, off, unsigned((target >> off) & 3)); }
        for (u64 s = 0; s + suff < flank; s++) { u64 off = 4 * s + 2; write2(x, off, unsigned((target >> off) & 3)); }
    };
    auto hi_missing = [&](size_t x) -> int {
        const u64 pref = lst[x].m_pref_size, suff = lst[x].m_suff_size;
        int hp = (pref < flank) ? int(4 * (flank - pref - 1))     : -1;
        int hs = (suff < flank) ? int(4 * (flank - suff - 1) + 2) : -1;
        return std::max(hp, hs);
    };
    std::vector<char> vb(n, 0);
    for (size_t i = 0; i < n; i++)
        if (lst[i].m_pref_size == (uint16_t)flank && lst[i].m_suff_size == (uint16_t)flank) vb[i] = 1;
    if (!vb[0])     { set_missing(0, u128(0)); vb[0] = 1; }
    if (!vb[n - 1]) { set_missing(n - 1, ~u128(0)); vb[n - 1] = 1; }

    std::vector<int> offsets;            // interleaved nucleotide positions, low->high
    for (u64 s = 0; s < flank; s++) { offsets.push_back(int(4 * s)); offsets.push_back(int(4 * s + 2)); }
    for (int c : offsets) {
        size_t left = 0; while (left < n && !vb[left]) left++;
        for (size_t b = left + 1; b < n; b++) {
            if (!vb[b]) continue;
            size_t li = left;
            for (size_t x = li + 1; x < b; x++) {
                if (hi_missing(x) == c) {
                    u128 vi = getv(li), vj = getv(b);
                    u128 target = (vj >= vi) ? vi + (vj - vi) * u128(x - li) / u128(b - li) : vi;
                    set_missing(x, target);
                    vb[x] = 1; li = x;
                }
            }
            left = b;
        }
    }
    return lst;
}

// forward decls (defined below)
template <typename kuint>
std::vector<km::Skmer<kuint>> midpoint_fill(const std::vector<km::Skmer<kuint>>&, km::SkmerManipulator<kuint>&);
template <typename kuint>
std::vector<std::vector<uint8_t>> bounded_query(const std::vector<km::Skmer<kuint>>&,
        km::SkmerManipulator<kuint>&, const std::vector<km::Skmer<kuint>>&, int);

// Extract MINIMAL examples where the interpolation fill leaves two consecutive entries out
// of order at some column c (key_c[i] > key_c[i+1]). Decoded slot-by-slot:
//   flanks printed high->low significance; UPPERCASE = present real nucleotide,
//   lowercase = filled hole, '_' = slot outside column c's window (masked out of key_c).
TEST(InterpExperiment, MinimalDisplacementExamples) {
    if (std::getenv("SKLIB_BENCH") == nullptr) GTEST_SKIP() << "set SKLIB_BENCH=1";
    using kuint = uint64_t;
    const u64 k = 8, m = 4, flank = k - m;
    km::SkmerManipulator<kuint> manip{k, m};
    auto en = enumerate(manip, random_seq(160, 3));
    km::sortedlist::SortedVirtualSkmerList<kuint> L(k, m);
    L.generate_sorted_list_from_enumeration(en);
    auto base = L.get_list();
    auto yo = interpolated_fill(base, manip, k, m);
    using kpair = typename km::Skmer<kuint>::pair;
    const auto kmask = manip.get_k_mask();
    static const char NUC[] = {'A', 'C', 'T', 'G'};
    auto nib = [&](const km::Skmer<kuint>& e, u64 off) {
        return unsigned((off < 64 ? (e.m_pair.m_value[0] >> off) : (e.m_pair.m_value[1] >> (off - 64))) & 3);
    };
    auto render = [&](const km::Skmer<kuint>& e, u64 c) {
        const u64 pref = e.m_pref_size, suff = e.m_suff_size;
        std::string s = "P:";
        for (long sl = (long)flank - 1; sl >= 0; sl--) {
            char ch = NUC[nib(e, 4 * sl)];
            if ((u64)sl < flank - pref) ch = char(ch - 'A' + 'a');     // filled hole
            if ((u64)sl < c) ch = '_';                                  // outside column c
            s += ch;
        }
        s += " S:";
        for (long sl = (long)flank - 1; sl >= 0; sl--) {
            char ch = NUC[nib(e, 4 * sl + 2)];
            if ((u64)sl < flank - suff) ch = char(ch - 'A' + 'a');
            if ((u64)sl < flank - c) ch = '_';
            s += ch;
        }
        return s;
    };
    int shown = 0;
    for (u64 c = 0; c <= flank && shown < 8; c++) {
        for (size_t i = 0; i + 1 < yo.size() && shown < 8; i++) {
            kpair ka{yo[i].m_pair & kmask[c]}, kb{yo[i + 1].m_pair & kmask[c]};
            if (kb < ka) {  // i is before i+1 in the list, but key_c(i) > key_c(i+1)
                const u64 miA = (yo[i].m_pair.m_value[0] >> (4 * flank)) & ((1ull << (2 * m)) - 1);
                const u64 miB = (yo[i + 1].m_pair.m_value[0] >> (4 * flank)) & ((1ull << (2 * m)) - 1);
                std::cerr << "\n[ex] VIOLATION col=" << c << " positions " << i << " < " << i + 1
                          << "  (mini A=" << miA << " B=" << miB << (miA == miB ? " same" : " DIFF") << ")\n";
                std::cerr << "  A[" << i << "] pref=" << yo[i].m_pref_size << " suff=" << yo[i].m_suff_size
                          << "  filled " << render(yo[i], c) << "  | unfilled " << render(base[i], c)
                          << "  key_c=0x" << std::hex << (uint64_t)(ka.m_value[0]) << std::dec << "\n";
                std::cerr << "  B[" << i + 1 << "] pref=" << yo[i + 1].m_pref_size << " suff=" << yo[i + 1].m_suff_size
                          << "  filled " << render(yo[i + 1], c) << "  | unfilled " << render(base[i + 1], c)
                          << "  key_c=0x" << std::hex << (uint64_t)(kb.m_value[0]) << std::dec << "\n";
                shown++;
            }
        }
    }
    std::cerr << "[ex] total shown=" << shown << " (list size " << yo.size() << ")\n";
    SUCCEED();
}

TEST(InterpExperiment, YoannFill) {
    if (std::getenv("SKLIB_BENCH") == nullptr) GTEST_SKIP() << "set SKLIB_BENCH=1 (Release build)";
    using kuint = uint64_t;
    const u64 k = 21, m = 11;
    for (size_t G : {size_t{200000}, size_t{1000000}, size_t{4000000}}) {
        km::SkmerManipulator<kuint> manip{k, m};
        auto en = enumerate(manip, random_seq(G, 5));
        km::sortedlist::SortedVirtualSkmerList<kuint> L(k, m);
        L.generate_sorted_list_from_enumeration(en);
        std::vector<km::Skmer<kuint>> base = L.get_list();
        auto present = en;
        auto absent = enumerate(manip, random_seq(G, 777));
        auto gt_present = L.query_skmer_batch(present);
        auto gt_absent = L.query_skmer_batch(absent);

        auto yo = interpolated_fill(base, manip, k, m);
        long disp = max_col_displacement(yo, manip, k, m);
        long colv = per_column_violations(yo, manip, k, m);
        long mid_disp = max_col_displacement(midpoint_fill(base, manip), manip, k, m);
        std::cerr << "[yoann] G=" << G << " entries=" << base.size()
                  << " yoann_max_disp=" << disp << " yoann_per_col_viol=" << colv
                  << " midpoint_max_disp=" << mid_disp << std::endl;
        for (int W : {0, 8, 16, 32}) {
            long fp = mismatches(gt_present, bounded_query(yo, manip, present, W));
            long fa = mismatches(gt_absent, bounded_query(yo, manip, absent, W));
            std::cerr << "[yoann] G=" << G << " W=" << W
                      << " present_mismatch=" << fp << " absent_mismatch=" << fa << std::endl;
        }
    }
}

// midpoint order-preserving fill (≈ Yoann's interpolation): per-entry value = achievable
// nearest the midpoint of the carry-min and carry-max monotone envelopes.
template <typename kuint>
std::vector<km::Skmer<kuint>> midpoint_fill(const std::vector<km::Skmer<kuint>>& base,
                                            km::SkmerManipulator<kuint>& manip) {
    const size_t n = base.size();
    std::vector<u64> plo(n), phi(n), amlo(n), amhi(n);
    for (size_t i = 0; i < n; i++) {
        auto am = manip.absent_slot_mask(base[i].m_pref_size, base[i].m_suff_size);
        amlo[i] = am.m_value[0]; amhi[i] = am.m_value[1];
        plo[i] = base[i].m_pair.m_value[0] & ~amlo[i];
        phi[i] = base[i].m_pair.m_value[1] & ~amhi[i];
    }
    std::vector<V128> lo(n), hi(n), mid(n);
    for (size_t i = 0; i < n; i++) { u64 tl = i ? lo[i-1].lo : 0, th = i ? lo[i-1].hi : 0; lo[i] = smallest_ge(plo[i], phi[i], amlo[i], amhi[i], tl, th); }
    for (size_t ii = 0; ii < n; ii++) { size_t i = n-1-ii; u64 tl = ii ? hi[i+1].lo : ~u64(0), th = ii ? hi[i+1].hi : ~u64(0); hi[i] = largest_le(plo[i], phi[i], amlo[i], amhi[i], tl, th); }
    for (size_t i = 0; i < n; i++) {
        u64 sl = lo[i].lo + hi[i].lo; u64 cy = (sl < lo[i].lo) ? 1 : 0; u64 sh = lo[i].hi + hi[i].hi + cy;
        u64 al = (sl >> 1) | (sh << 63), ah = sh >> 1;
        u64 tl = al, th = ah; if (i && ge128(mid[i-1].lo, mid[i-1].hi, tl, th)) { tl = mid[i-1].lo; th = mid[i-1].hi; }
        mid[i] = smallest_ge(plo[i], phi[i], amlo[i], amhi[i], tl, th);
    }
    std::vector<km::Skmer<kuint>> filled = base;
    for (size_t i = 0; i < n; i++) { filled[i].m_pair.m_value[0] = mid[i].lo; filled[i].m_pair.m_value[1] = mid[i].hi; }
    return filled;
}

// Hole-aware query with a bounded ±W fallback: binary-search by column key, then linear-
// scan a window around the convergence point for a valid match.
template <typename kuint>
std::vector<std::vector<uint8_t>> bounded_query(const std::vector<km::Skmer<kuint>>& list,
        km::SkmerManipulator<kuint>& manip, const std::vector<km::Skmer<kuint>>& qs, int W) {
    const int64_t N = (int64_t)list.size();
    std::vector<std::vector<uint8_t>> out(qs.size());
    for (size_t qi = 0; qi < qs.size(); qi++) {
        auto bnd = manip.get_valid_kmer_bounds(qs[qi]);
        if (bnd.second < bnd.first) continue;
        const uint64_t nk = bnd.second - bnd.first + 1;
        std::vector<uint8_t> res(nk, 0);
        for (uint64_t off = 0; off < nk; off++) {
            const uint64_t c = bnd.first + off;
            int64_t lo = 0, hi = N - 1, conv = 0; bool found = false;
            while (lo <= hi) {
                int64_t mid = lo + ((hi - lo) >> 1); conv = mid;
                int cmp = manip.kmer_compare(qs[qi], list[mid], c);
                if (cmp == 0) break; else if (cmp < 0) hi = mid - 1; else lo = mid + 1;
            }
            int64_t a = std::max<int64_t>(0, std::min(lo, conv) - W);
            int64_t b = std::min<int64_t>(N - 1, std::max(lo, conv) + W);
            for (int64_t i = a; i <= b; i++)
                if (manip.has_valid_kmer(list[i], c) && manip.kmer_compare(qs[qi], list[i], c) == 0) { found = true; break; }
            res[off] = found ? 1 : 0;
        }
        out[qi] = std::move(res);
    }
    return out;
}

TEST(InterpExperiment, BoundedFallback) {
    if (std::getenv("SKLIB_BENCH") == nullptr) GTEST_SKIP() << "set SKLIB_BENCH=1 (Release build)";
    using kuint = uint64_t;
    const u64 k = 21, m = 11;
    for (size_t G : {size_t{1000000}, size_t{4000000}}) {
        km::SkmerManipulator<kuint> manip{k, m};
        auto en = enumerate(manip, random_seq(G, 5));
        km::sortedlist::SortedVirtualSkmerList<kuint> L(k, m);
        L.generate_sorted_list_from_enumeration(en);
        auto filled = midpoint_fill(L.get_list(), manip);
        auto present = en;
        auto absent = enumerate(manip, random_seq(G, 777));
        auto gt_present = L.query_skmer_batch(present);
        auto gt_absent = L.query_skmer_batch(absent);
        std::cerr << "[bounded] G=" << G << " entries=" << filled.size() << std::endl;
        for (int W : {0, 8, 16, 32, 64}) {
            long fp = mismatches(gt_present, bounded_query(filled, manip, present, W));
            long fa = mismatches(gt_absent, bounded_query(filled, manip, absent, W));
            std::cerr << "[bounded] G=" << G << " W=" << W
                      << " present_mismatch=" << fp << " absent_mismatch=" << fa << std::endl;
        }
    }
}

TEST(InterpExperiment, ColumnAwareFill) {
    if (std::getenv("SKLIB_BENCH") == nullptr) GTEST_SKIP() << "set SKLIB_BENCH=1 (Release build)";
    using kuint = uint64_t;
    const u64 k = 21, m = 11;
    for (size_t G : {size_t{200000}, size_t{1000000}}) {
        km::SkmerManipulator<kuint> manip{k, m};
        auto en = enumerate(manip, random_seq(G, 5));
        km::sortedlist::SortedVirtualSkmerList<kuint> L(k, m);
        L.generate_sorted_list_from_enumeration(en);
        std::vector<km::Skmer<kuint>> base = L.get_list();
        auto present = en;
        auto absent = enumerate(manip, random_seq(G, 777));
        auto gt_present = L.query_skmer_batch(present);
        auto gt_absent = L.query_skmer_batch(absent);

        auto filled = column_aware_fill(base, manip, k, m);
        long colv = per_column_violations(filled, manip, k, m);
        km::sortedlist::SortedVirtualSkmerList<kuint> F(k, m);
        F.add_list(filled);
        long fnp = mismatches(gt_present, F.query_skmer_batch_substrate(present));
        long fna = mismatches(gt_absent, F.query_skmer_batch_substrate(absent));
        std::cerr << "[colaware] G=" << G << " entries=" << base.size()
                  << " per_col_viol=" << colv
                  << " substrate_query_mismatch present=" << fnp << " absent=" << fna << std::endl;
    }
}

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
            // max per-column displacement = window W a bounded ±W fallback would need to
            // be exact (how far an entry is from its sorted position in each column key).
            long maxdisp = 0;
            { using kpair = typename km::Skmer<kuint>::pair;
              const std::vector<kpair> kmask = manip.get_k_mask();
              std::vector<uint32_t> ord(n);
              for (u64 c2 = 0; c2 <= k - m; c2++) {
                  for (uint32_t z = 0; z < n; z++) ord[z] = z;
                  std::stable_sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
                      return (filled[a].m_pair & kmask[c2]) < (filled[b].m_pair & kmask[c2]); });
                  for (uint32_t r = 0; r < (uint32_t)n; r++) {
                      long d = (long)ord[r] - (long)r; if (d < 0) d = -d; if (d > maxdisp) maxdisp = d;
                  }
              }
            }
            km::sortedlist::SortedVirtualSkmerList<kuint> F(k, m);
            F.add_list(filled);
            long fn_present = mismatches(gt_present, F.query_skmer_batch_substrate(present));
            long fn_absent  = mismatches(gt_absent,  F.query_skmer_batch_substrate(absent));
            std::cerr << "[interp] G=" << G << " " << name
                      << " stuck=" << stuck << " mpair_inversions=" << mono
                      << " per_col_viol=" << colv << " max_col_displacement=" << maxdisp
                      << " substrate_query_mismatch present=" << fn_present
                      << " absent=" << fn_absent << std::endl;
        };
        std::cerr << "[interp] G=" << G << " entries=" << n << std::endl;
        eval("carry-min", lo, stuck_min);
        eval("carry-max", hi, stuck_max);
        eval("midpoint ", mid, 0);
    }
}
