#include "algorithms/ColinearChaining.hpp"

// --- Colienar chaining algorithm and data structures---


namespace km
{
namespace chaining
{


/** Colinear chaining algorithm to select a compatible set of overlaps. The overlaps are compatible if their
 * coordinates are not crossing and if they do not have common coordinates.
 * WARNING: The algorithme will change the order of the input vector.
 *
 * @param begin iterator to the first element of the vector
 * @param end iterator to the last element of the vector
 *
 * @return a vector of overlaps containing compatible overlaps
 **/
std::vector<overlap> colinear_chaining(std::vector<overlap>::iterator begin, std::vector<overlap>::iterator end)
{
    size_t const n = static_cast<size_t>(end - begin);
    if (n == 0)
        return {};

    // Longest chain with STRICTLY increasing first AND second coordinate.
    // Process the overlaps in first-coordinate order (ties: larger second first, so
    // two overlaps that share a first coordinate are never chained to each other),
    // keeping a Fenwick tree of the best chain length over the second coordinate.
    // When overlap X is processed, the tree holds exactly the overlaps before X in
    // first order; among them, those with second < X.second also have first < X.first,
    // so they are the legal predecessors.
    //
    // We sort and index the caller's range in place (it is not reused after, and the
    // contract already documents that the input order changes) — no private copy, so
    // the per-overlap state is just the DP vectors below. get_candidate_overlaps' hash join
    // already emits candidates in (first asc, second desc) order, so the sort is usually a no-op;
    // an O(n) is_sorted guard skips the already-ordered O(n log n) sort (identical result).
    auto chain_less = [](overlap const& a, overlap const& b) {
        return a.first != b.first ? a.first < b.first : a.second > b.second;
    };
    if (!std::is_sorted(begin, end, chain_less))
        std::sort(begin, end, chain_less);

    // Coordinate-compress the second coordinates.
    std::vector<uint64_t> ys;
    ys.reserve(n);
    for (size_t i = 0; i < n; i++) ys.push_back(begin[i].second);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    auto rank = [&ys](uint64_t y) -> size_t {
        return static_cast<size_t>(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin());
    };
    size_t const S = ys.size();

    // Fenwick (1-indexed) prefix-maximum. A cell records the best chain ending at a
    // processed overlap, referenced by its index in the sorted range: longer wins;
    // ties go to the smaller (first, second) overlap, which reproduces the chain
    // tie-breaking the unit tests pin. Indices and lengths are 32-bit: both are
    // bounded by the candidate count, far below 2^31 for any in-RAM bucket, halving
    // the DP/Fenwick footprint that dominates repeat-rich buckets.
    struct Cell { uint32_t length; int32_t end_idx; };
    Cell const empty {0, -1};
    auto better = [begin](Cell const& a, Cell const& b) -> bool {
        if (a.length != b.length) return a.length > b.length;
        if (a.length == 0)        return false;
        overlap const& ae = begin[a.end_idx];
        overlap const& be = begin[b.end_idx];
        if (ae.first != be.first) return ae.first < be.first;
        return ae.second < be.second;
    };
    std::vector<Cell> bit(S + 1, empty);
    auto bit_update = [&](size_t pos, Cell const& v) {
        for (size_t i = pos; i <= S; i += i & (~i + 1))
            if (better(v, bit[i])) bit[i] = v;
    };
    auto bit_prefix_max = [&](size_t pos) -> Cell {
        Cell res {empty};
        for (size_t i = pos; i > 0; i -= i & (~i + 1))
            if (better(bit[i], res)) res = bit[i];
        return res;
    };

    // Per-overlap DP: predecessor index (or -1) and chain length, indexed by the
    // overlap's position in the sorted range.
    std::vector<int32_t> prev_idx(n, -1);
    std::vector<uint32_t> length(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        size_t const r = rank(begin[i].second);                   // ranks [0, r-1] have second < o.second
        Cell const best = (r == 0) ? empty : bit_prefix_max(r);   // Fenwick positions [1, r]
        length[i]   = best.length + 1;
        prev_idx[i] = (best.length == 0) ? int32_t{-1} : best.end_idx;
        bit_update(r + 1, Cell {length[i], static_cast<int32_t>(i)});
    }

    // Chain end: the longest chain, ties broken by larger (first, second).
    int32_t end_idx = -1;
    uint32_t max_len = 0;
    for (size_t i = 0; i < n; i++)
    {
        uint32_t const len = length[i];
        if (end_idx < 0 || len > max_len ||
            (len == max_len &&
             (begin[i].first > begin[static_cast<size_t>(end_idx)].first ||
              (begin[i].first == begin[static_cast<size_t>(end_idx)].first &&
               begin[i].second > begin[static_cast<size_t>(end_idx)].second))))
        {
            max_len = len;
            end_idx = static_cast<int32_t>(i);
        }
    }

    std::vector<overlap> chain(max_len);
    int32_t cur {end_idx};
    for (uint32_t i = max_len; i > 0; i--)
    {
        chain[i - 1] = begin[static_cast<size_t>(cur)];
        cur = prev_idx[static_cast<size_t>(cur)];
    }
    return chain;
}

}} // namespace sorting // namespace km
