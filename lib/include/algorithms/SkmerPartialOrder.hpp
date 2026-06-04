#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <ostream>
#include <unordered_map>

#include <io/Skmer.hpp>

#ifndef SKMERPARTIALORDER_H
#define SKMERPARTIALORDER_H

// EXPERIMENTAL — real-nucleotide ordering of super-k-mers, used to expose where the
// reconciliation orders skmers using *absent* (masked) nucleotides instead of real ones.
//
// The production list orders skmers by their interleaved m_pair *masks included*: absent
// flank nucleotides are filled with 0b11 (`mask_absent_nucleotides`), equal to G
// numerically, so two skmers that share their central nucleotides but have flanks of
// different lengths get ordered by that padding rather than by real nucleotides.
//
// `compare_real` orders strictly by the *real* nucleotides and returns INCOMPARABLE
// exactly when the first (most significant) difference falls on an absent slot — i.e.
// precisely when `Skmer::operator<` would have decided via the mask. The reconciliation
// (see VirtualSkmer's *_partial methods) uses this to record, during the merge, every
// ordering it would have settled with the mask; those links are unioned and surfaced as
// groups over the final compacted list.

namespace km
{
namespace experiment
{

enum class RealOrder { Less, Greater, Equal, Incomparable };

// Extract the 2-bit nucleotide stored at bit position `bitpos` of a skmer pair.
template<typename kuint>
inline uint8_t nibble_at(const km::Skmer<kuint>& sk, uint64_t bitpos) {
    return static_cast<uint8_t>(static_cast<kuint>(sk.m_pair >> bitpos) & static_cast<kuint>(0b11));
}

/** Real-nucleotide "significance key" of a skmer: its present nucleotides in the exact
 * order `Skmer::operator<` weighs them (most significant first), with absent flank
 * positions *omitted* (treated as end-of-string instead of 0b11 padding).
 *
 * Order = minimizer (top 2m bits, MSB nucleotide first) then flank slots from central
 * (high slot) to outer (slot 0); within a slot the suffix half (bits 4s+2) outranks the
 * prefix half (bits 4s). Built from the stored φ-permuted bits, so a lexicographic
 * compare on keys reproduces the production order on real nucleotides while making the
 * mask invisible. A slot is present per `mask_absent_nucleotides`: absent flank slots
 * are the low indices (prefix slot i absent iff i < flank - pref_size, likewise suffix).
 */
template<typename kuint>
std::vector<uint8_t> significance_key(const km::SkmerManipulator<kuint>& manip,
                                      const km::Skmer<kuint>& sk) {
    const int64_t flank    = static_cast<int64_t>(manip.k - manip.m);
    const int64_t sk_size  = static_cast<int64_t>(2 * manip.k - manip.m);
    const int64_t buf_pref = (sk_size + 1) / 2;
    const int64_t buf_suff = sk_size / 2;
    const int64_t pref     = static_cast<int64_t>(sk.m_pref_size);
    const int64_t suff     = static_cast<int64_t>(sk.m_suff_size);

    std::vector<uint8_t> key;
    key.reserve(static_cast<size_t>(manip.m) + sk.m_pref_size + sk.m_suff_size);

    for (int64_t s = buf_pref - 1; s >= 0; --s) {
        // suffix half is more significant than the prefix half of the same slot
        if (s <= buf_suff - 1 && s >= flank - suff)
            key.push_back(nibble_at(sk, static_cast<uint64_t>(4 * s + 2)));
        if (s >= flank - pref)
            key.push_back(nibble_at(sk, static_cast<uint64_t>(4 * s)));
    }
    return key;
}

/** Four-way comparison of two significance keys: lexicographic where a proper prefix is
 * INCOMPARABLE (not less) — the prefix relation is exactly the absent-nucleotide ambiguity. */
inline RealOrder compare_keys(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i)
        if (a[i] != b[i]) return a[i] < b[i] ? RealOrder::Less : RealOrder::Greater;
    if (a.size() == b.size()) return RealOrder::Equal;
    return RealOrder::Incomparable; // one is a proper real-nucleotide prefix of the other
}

// Compare two skmers by real nucleotides. INCOMPARABLE ⇔ operator< would have decided
// the order using an absent (masked) nucleotide.
template<typename kuint>
RealOrder compare_real(const km::SkmerManipulator<kuint>& manip,
                       const km::Skmer<kuint>& a, const km::Skmer<kuint>& b) {
    return compare_keys(significance_key(manip, a), significance_key(manip, b));
}

// Minimal union-find over freshly-minted ids, used to coalesce the incomparable pairs
// detected during the partial reconciliation into connected-component groups.
struct UnionFind {
    std::vector<uint32_t> parent;

    uint32_t make() {
        const uint32_t id = static_cast<uint32_t>(parent.size());
        parent.push_back(id);
        return id;
    }
    uint32_t find(uint32_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    }
    void unite(uint32_t a, uint32_t b) {
        a = find(a); b = find(b);
        if (a != b) parent[a] = b;
    }
};

// True when skmer i belongs to a group of size > 1 (an incomparable cluster).
inline std::vector<uint8_t> grouped_flags(const std::vector<uint32_t>& group_of) {
    std::unordered_map<uint32_t, uint64_t> counts;
    for (uint32_t g : group_of) counts[g]++;
    std::vector<uint8_t> grouped(group_of.size());
    for (size_t i = 0; i < group_of.size(); ++i) grouped[i] = counts[group_of[i]] > 1;
    return grouped;
}

// Overwrite the 2-bit nucleotide stored at bit `bitpos` of a skmer pair with `v` (0..3).
template<typename kuint>
inline void set_nibble(km::Skmer<kuint>& sk, uint64_t bitpos, int v) {
    using kpair = typename km::Skmer<kuint>::pair;
    sk.m_pair = sk.m_pair & ~(kpair(static_cast<kuint>(0b11)) << bitpos);
    sk.m_pair = sk.m_pair | (kpair(static_cast<kuint>(v & 0b11)) << bitpos);
}

/** Per-column "virtual nucleotide" fill.
 *
 * The list is sorted PER COLUMN: for each k-mer position c, the skmers that have a valid
 * k-mer at c are, in list order, non-decreasing by kmer_compare(.,.,c). A skmer truncated
 * at c (no valid k-mer there) forces the per-column dichotomic search into a linear scan.
 *
 * For each column c, a truncated skmer S sits between its two nearest valid-at-c neighbours
 * A (above) and B (below). Where A and B AGREE on a missing nucleotide of S's k-mer-at-c,
 * S's placement between them is unambiguous, so we propose that value; where they disagree
 * (or different columns propose conflicting values for the same physical slot), the slot is
 * left absent — S stays incomparable at the columns needing it. This NEVER reorders a real
 * k-mer (only absent slots are written; valid k-mers are untouched) and only ever places a
 * filled k-mer strictly between its neighbours, so the per-column order is preserved.
 *
 * pref/suff sizes are unchanged. `filled_pref[i]`/`filled_suff[i]` return the bitmask of
 * prefix/suffix flank slots actually filled for skmer i (slot p ↔ bit p), for display and
 * measurement. Slot encoding: prefix slot p at bit 4p, suffix slot p at bit 4p+2; a prefix
 * slot p is present iff p >= flank-pref, a suffix slot p iff p >= flank-suff.
 */
template<typename kuint>
void fill_virtual_nucleotides(const km::SkmerManipulator<kuint>& manip,
                              std::vector<km::Skmer<kuint>>& skmers,
                              std::vector<uint32_t>& filled_pref,
                              std::vector<uint32_t>& filled_suff) {
    const int64_t flank = static_cast<int64_t>(manip.k) - static_cast<int64_t>(manip.m);
    const size_t n = skmers.size();
    filled_pref.assign(n, 0);
    filled_suff.assign(n, 0);
    if (flank <= 0 || n == 0) return;

    // Probe priority. The query's per-column dichotomic search probes list index mid=(lo+hi)/2
    // first, then the quarters, etc. — the nodes of the binary-search tree over [0, n-1]. A
    // shallow node is hit by a large fraction of queries, so reverting (leaving incomparable)
    // a shallow skmer is far costlier than a deep one. weight[i] = 2^(maxdepth-depth(i)) makes
    // the correctness pass below sacrifice the deepest positions first.
    std::vector<int64_t> weight(n, 1);
    {
        std::vector<int32_t> depth(n, 0);
        int32_t maxd = 0;
        struct Range { int64_t lo, hi; int32_t d; };
        std::vector<Range> st{ { 0, static_cast<int64_t>(n) - 1, 0 } };
        while (!st.empty()) {
            const Range r = st.back(); st.pop_back();
            if (r.lo > r.hi) continue;
            const int64_t mid = (r.lo + r.hi) / 2;
            depth[mid] = r.d; if (r.d > maxd) maxd = r.d;
            st.push_back({ r.lo, mid - 1, r.d + 1 });
            st.push_back({ mid + 1, r.hi, r.d + 1 });
        }
        for (size_t i = 0; i < n; ++i)
            weight[i] = int64_t{1} << std::min<int32_t>(maxd - depth[i], 50);
    }

    // Per (skmer, slot) proposal: -1 unset, 0..3 a value. First column to want the slot wins
    // (we never drop a slot on conflict — that would leave shallow probe skmers unfilled); the
    // weighted correctness pass below reverts whatever still mis-places, deepest first.
    std::vector<int8_t> prop_pref(n * flank, -1), prop_suff(n * flank, -1);
    auto propose = [](std::vector<int8_t>& a, size_t idx, int v) {
        if (a[idx] == -1) a[idx] = static_cast<int8_t>(v);
    };

    std::vector<int64_t> above(n), below(n);
    for (int64_t c = 0; c <= flank; ++c) {
        // valid at c iff the k-mer-at-c fits: pref >= flank-c (enough prefix) and suff >= c.
        auto valid = [&](size_t i) {
            return static_cast<int64_t>(skmers[i].m_pref_size) >= flank - c
                && static_cast<int64_t>(skmers[i].m_suff_size) >= c;
        };
        int64_t last = -1;
        for (size_t i = 0; i < n; ++i) { above[i] = last; if (valid(i)) last = static_cast<int64_t>(i); }
        last = -1;
        for (size_t i = n; i-- > 0; ) { below[i] = last; if (valid(i)) last = static_cast<int64_t>(i); }

        for (size_t i = 0; i < n; ++i) {
            if (valid(i)) continue;                 // S has a real k-mer at c
            const int64_t A = above[i], B = below[i];
            const int64_t pref = skmers[i].m_pref_size, suff = skmers[i].m_suff_size;
            // Optimistic proposal: round S's missing outer toward the LOWER anchor A (or the
            // upper anchor B if there is no A). This makes S usable whenever its real centre
            // already ranks it between the anchors (then any outer is order-safe) or equals
            // theirs (then it matches A). Conflicting proposals from different columns drop
            // the slot; the LNDS correctness pass below reverts whatever still mis-places.
            auto constrain = [&](std::vector<int8_t>& prop, size_t idx, uint64_t bit) {
                int req = -1;
                if (A >= 0)      req = nibble_at(skmers[A], bit);
                else if (B >= 0) req = nibble_at(skmers[B], bit);
                if (req >= 0) propose(prop, idx, req);
            };
            for (int64_t p = c; p < flank - pref; ++p)            // missing prefix slots
                constrain(prop_pref, static_cast<size_t>(i) * flank + p, static_cast<uint64_t>(4 * p));
            for (int64_t p = flank - c; p < flank - suff; ++p)   // missing suffix slots
                constrain(prop_suff, static_cast<size_t>(i) * flank + p, static_cast<uint64_t>(4 * p + 2));
        }
    }

    for (size_t i = 0; i < n; ++i)
        for (int64_t p = 0; p < flank; ++p) {
            const int8_t v1 = prop_pref[static_cast<size_t>(i) * flank + p];
            if (v1 >= 0) { set_nibble(skmers[i], static_cast<uint64_t>(4 * p), v1);
                           filled_pref[i] |= (uint32_t{1} << p); }
            const int8_t v2 = prop_suff[static_cast<size_t>(i) * flank + p];
            if (v2 >= 0) { set_nibble(skmers[i], static_cast<uint64_t>(4 * p + 2), v2);
                           filled_suff[i] |= (uint32_t{1} << p); }
        }

    // Correctness pass. A filled skmer can still be mis-placed at a column where its *real*
    // central part already ranks it outside its neighbours (no fill can fix that), or where
    // a shared slot was forced to a value good for another column. For each column we keep a
    // LONGEST non-decreasing subsequence of its usable skmers (by kmer_compare) and revert
    // only the rest — and only the slots that column needs, so a skmer stays usable at the
    // columns that don't need them. Iterating (a revert can change other columns) converges:
    // every revert clears at least one filled slot.
    auto valid_at = [&](size_t i, int64_t c) {
        return static_cast<int64_t>(skmers[i].m_pref_size) >= flank - c
            && static_cast<int64_t>(skmers[i].m_suff_size) >= c;
    };
    auto usable = [&](size_t i, int64_t c) {
        if (valid_at(i, c)) return true;
        const int64_t pref = skmers[i].m_pref_size, suff = skmers[i].m_suff_size;
        for (int64_t p = c; p < flank - pref; ++p)        if (!((filled_pref[i] >> p) & 1u)) return false;
        for (int64_t p = flank - c; p < flank - suff; ++p) if (!((filled_suff[i] >> p) & 1u)) return false;
        return true;
    };
    auto revert_col = [&](size_t i, int64_t c) {                  // revert only column c's slots
        const int64_t pref = skmers[i].m_pref_size, suff = skmers[i].m_suff_size;
        for (int64_t p = c; p < flank - pref; ++p)
            if ((filled_pref[i] >> p) & 1u) {
                set_nibble(skmers[i], static_cast<uint64_t>(4 * p), 3);
                filled_pref[i] &= ~(uint32_t{1} << p);
            }
        for (int64_t p = flank - c; p < flank - suff; ++p)
            if ((filled_suff[i] >> p) & 1u) {
                set_nibble(skmers[i], static_cast<uint64_t>(4 * p + 2), 3);
                filled_suff[i] &= ~(uint32_t{1} << p);
            }
    };

    // Valid skmers are mandatory anchors (never reverted). They split each column's usable
    // sequence into segments; in a segment the filled skmers must lie within the surrounding
    // anchors' values, and among those we keep a LONGEST non-decreasing subsequence (so the
    // fewest fills are reverted) and revert the rest — column-c slots only, so a skmer stays
    // usable at columns that don't need them. Iterate until stable (each revert clears ≥1
    // filled slot, so it terminates).
    bool changed = true;
    std::vector<size_t> seg, cand;
    std::vector<int64_t> dp;
    std::vector<int> par;
    std::vector<char> kept;
    auto process_seg = [&](int64_t L, int64_t R, int64_t c) {
        cand.clear();
        for (size_t f : seg) {                              // drop fills outside the anchor bounds
            const bool below = (L >= 0 && manip.kmer_compare(skmers[L], skmers[f], c) > 0);
            const bool above = (R >= 0 && manip.kmer_compare(skmers[f], skmers[R], c) > 0);
            if (below || above) { revert_col(f, c); changed = true; }
            else cand.push_back(f);
        }
        const size_t s = cand.size();
        if (s < 2) return;
        // maximum-WEIGHT non-decreasing subsequence: keep the costliest-to-lose probe skmers,
        // revert the rest (so deep / rarely-probed positions are sacrificed first).
        dp.assign(s, 0);
        par.assign(s, -1);
        size_t best = 0;
        for (size_t j = 0; j < s; ++j) {
            dp[j] = weight[cand[j]];
            for (size_t k = 0; k < j; ++k)
                if (dp[k] + weight[cand[j]] > dp[j]
                    && manip.kmer_compare(skmers[cand[k]], skmers[cand[j]], c) <= 0) {
                    dp[j] = dp[k] + weight[cand[j]];
                    par[j] = static_cast<int>(k);
                }
            if (dp[j] > dp[best]) best = j;
        }
        kept.assign(s, 0);
        for (int cur = static_cast<int>(best); cur >= 0; cur = par[cur]) kept[cur] = 1;
        for (size_t j = 0; j < s; ++j)
            if (!kept[j]) { revert_col(cand[j], c); changed = true; }
    };
    auto run_revert = [&]() {
        changed = true;
        while (changed) {
            changed = false;
            for (int64_t c = 0; c <= flank; ++c) {
                int64_t prev_valid = -1;
                seg.clear();
                for (size_t i = 0; i < n; ++i) {
                    if (!usable(i, c)) continue;
                    if (valid_at(i, c)) { process_seg(prev_valid, static_cast<int64_t>(i), c);
                                          seg.clear(); prev_valid = static_cast<int64_t>(i); }
                    else seg.push_back(i);
                }
                process_seg(prev_valid, -1, c);
            }
        }
    };
    run_revert();

    // ---- Probe-skmer optimisation ----
    // The first dichotomic probes (list middle, quarters, …) are hit by (almost) every query,
    // so it pays to make those positions usable at as many columns as possible. For the
    // highest-weight (shallowest) positions, try re-filling the whole absent flank by copying a
    // neighbouring VALID skmer, and keep the copy that maximises the columns where this skmer
    // fits between its usable neighbours. The (weighted) correctness pass is then re-run to keep
    // everything consistent — it preserves these high-weight gains where it can.
    {
        std::vector<size_t> probes(n);
        for (size_t i = 0; i < n; ++i) probes[i] = i;
        std::sort(probes.begin(), probes.end(), [&](size_t a, size_t b) { return weight[a] > weight[b]; });
        const size_t LIMIT = std::min<size_t>(n, 64);

        auto set_ref = [&](size_t i, size_t V) {        // copy V's flank into i's absent slots
            const int64_t pref = skmers[i].m_pref_size, suff = skmers[i].m_suff_size;
            for (int64_t p = 0; p < flank - pref; ++p)
                set_nibble(skmers[i], static_cast<uint64_t>(4 * p), nibble_at(skmers[V], static_cast<uint64_t>(4 * p)));
            for (int64_t p = 0; p < flank - suff; ++p)
                set_nibble(skmers[i], static_cast<uint64_t>(4 * p + 2), nibble_at(skmers[V], static_cast<uint64_t>(4 * p + 2)));
            filled_pref[i] = (flank - pref > 0) ? ((uint32_t{1} << (flank - pref)) - 1) : 0;
            filled_suff[i] = (flank - suff > 0) ? ((uint32_t{1} << (flank - suff)) - 1) : 0;
        };
        auto fit_columns = [&](size_t i) {              // columns where i sits between its usable neighbours
            size_t cnt = 0;
            for (int64_t c = 0; c <= flank; ++c) {
                if (valid_at(i, c)) { ++cnt; continue; }
                int64_t L = -1, R = -1;
                for (int64_t j = static_cast<int64_t>(i) - 1; j >= 0; --j) if (usable(static_cast<size_t>(j), c)) { L = j; break; }
                for (size_t j = i + 1; j < n; ++j) if (usable(j, c)) { R = static_cast<int64_t>(j); break; }
                bool ok = true;
                if (L >= 0 && manip.kmer_compare(skmers[L], skmers[i], c) > 0) ok = false;
                if (R >= 0 && manip.kmer_compare(skmers[i], skmers[R], c) > 0) ok = false;
                if (ok) ++cnt;
            }
            return cnt;
        };

        std::vector<size_t> refs;
        for (size_t pi = 0; pi < LIMIT; ++pi) {
            const size_t i = probes[pi];
            if (static_cast<int64_t>(skmers[i].m_pref_size) == flank
             && static_cast<int64_t>(skmers[i].m_suff_size) == flank) continue; // nothing absent
            refs.clear();
            for (int64_t c = 0; c <= flank; ++c) {
                if (valid_at(i, c)) continue;
                for (int64_t j = static_cast<int64_t>(i) - 1; j >= 0; --j) if (valid_at(static_cast<size_t>(j), c)) { refs.push_back(static_cast<size_t>(j)); break; }
                for (size_t j = i + 1; j < n; ++j) if (valid_at(j, c)) { refs.push_back(j); break; }
            }
            km::Skmer<kuint> best = skmers[i];
            uint32_t bestP = filled_pref[i], bestS = filled_suff[i];
            size_t best_count = fit_columns(i);
            for (size_t V : refs) {
                set_ref(i, V);
                const size_t cnt = fit_columns(i);
                if (cnt > best_count) { best_count = cnt; best = skmers[i]; bestP = filled_pref[i]; bestS = filled_suff[i]; }
            }
            skmers[i] = best; filled_pref[i] = bestP; filled_suff[i] = bestS;
        }
    }
    run_revert();
}

// Human decode of a skmer. Mirrors VirtualSkmer's save_ascii convention (minimizer
// φ⁻¹-decoded, interleaved halves), padded full-width so shared centers align. Per slot:
//   present (real)                                  -> UPPERCASE
//   absent but filled (skmer `filled`, p>=flank-max) -> lowercase (interpolated value)
//   otherwise (grouped skmer, or slot never reached) -> '.'
// `max_pref`/`max_suff` are the list-wide maxima; below `flank-max` no skmer has the slot.
template<typename kuint>
std::string decode_with_dots(const km::SkmerManipulator<kuint>& manip,
                             const km::Skmer<kuint>& sk_in,
                             uint32_t filled_pref, uint32_t filled_suff) {
    static const char up[]  = {'A', 'C', 'T', 'G'};
    static const char low[] = {'a', 'c', 't', 'g'};
    km::Skmer<kuint> sk = sk_in;
    manip.unpermute_minimizer_slot(sk); // show the real minimizer, not the φ-mixed value

    const int64_t flank    = static_cast<int64_t>(manip.k - manip.m);
    const int64_t sk_size  = static_cast<int64_t>(2 * manip.k - manip.m);
    const int64_t buf_pref = (sk_size + 1) / 2;
    const int64_t buf_suff = sk_size / 2;
    const int64_t pref     = static_cast<int64_t>(sk.m_pref_size);
    const int64_t suff     = static_cast<int64_t>(sk.m_suff_size);

    std::string s;
    for (int64_t slot = 0; slot < buf_pref; ++slot) {
        if (slot >= flank - pref)                                s += up[nibble_at(sk, static_cast<uint64_t>(4 * slot))];
        else if (slot < flank && ((filled_pref >> slot) & 1u))   s += low[nibble_at(sk, static_cast<uint64_t>(4 * slot))];
        else                                                     s += '.';
    }
    s += ' ';
    for (int64_t slot = buf_suff - 1; slot >= 0; --slot) {
        if (slot >= flank - suff)                                s += up[nibble_at(sk, static_cast<uint64_t>(4 * slot + 2))];
        else if (slot < flank && ((filled_suff >> slot) & 1u))   s += low[nibble_at(sk, static_cast<uint64_t>(4 * slot + 2))];
        else                                                     s += '.';
    }
    return s;
}

/** Dump a compacted skmer list with its partial-order grouping. `group_of[i]` is the
 * (union-find canonical) group id of skmer i. Skmers sharing a group are incomparable
 * (the reconciliation could only have ordered them via the mask) and are printed as a
 * `GROUP (n=…)` block; comparable skmers print one per line. Items keep the list order
 * (first appearance of each group). */
template<typename kuint>
void dump_partial_ascii(const km::SkmerManipulator<kuint>& manip,
                        const std::vector<km::Skmer<kuint>>& skmers,
                        const std::vector<uint32_t>& group_of,
                        const std::vector<uint32_t>& filled_pref,
                        const std::vector<uint32_t>& filled_suff,
                        std::ostream& os) {
    std::vector<uint32_t> order;                              // group ids, first-appearance order
    std::unordered_map<uint32_t, std::vector<uint64_t>> members;
    for (uint64_t i = 0; i < skmers.size(); ++i) {
        const uint32_t g = group_of[i];
        if (members.find(g) == members.end()) order.push_back(g);
        members[g].push_back(i);
    }

    uint64_t n_groups = 0;
    for (uint32_t g : order) if (members[g].size() > 1) ++n_groups;

    os << "k=" << manip.k << " m=" << manip.m
       << " compacted_skmers=" << skmers.size()
       << " items=" << order.size()
       << " groups=" << n_groups << "\n";

    // `pos` is the skmer's index in the compacted list. Real nucleotides are UPPERCASE,
    // the per-column virtual fills lowercase, still-incomparable slots '.'.
    auto print_member = [&](uint64_t i, const char* indent) {
        os << indent << "pos " << i << "\t"
           << decode_with_dots(manip, skmers[i], filled_pref[i], filled_suff[i])
           << "  pref:" << skmers[i].m_pref_size
           << " suff:" << skmers[i].m_suff_size << "\n";
    };

    uint64_t idx = 0;
    for (uint32_t g : order) {
        const std::vector<uint64_t>& mem = members[g];
        if (mem.size() == 1) {
            os << idx << "\t";
            print_member(mem[0], "");
        } else {
            os << idx << "\tGROUP (n=" << mem.size() << ", positions";
            for (uint64_t i : mem) os << " " << i;
            os << "):\n";
            for (uint64_t i : mem) print_member(i, "\t  ");
        }
        ++idx;
    }
}

} // namespace experiment
} // namespace km

#endif // SKMERPARTIALORDER_H
