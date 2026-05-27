#include "algorithms/ColinearChaining.hpp"

// --- Colienar chaining algorithm and data structures---


namespace km
{
namespace chaining
{

    static overlap const null_overlap = {UINT64_MAX, UINT64_MAX};


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
    std::vector<overlap> ov(begin, end);
    size_t const n = ov.size();
    if (n == 0)
        return {};

    // Longest chain with STRICTLY increasing first AND second coordinate.
    // Process the overlaps in first-coordinate order (ties: larger second first, so
    // two overlaps that share a first coordinate are never chained to each other),
    // keeping a Fenwick tree of the best chain length over the second coordinate.
    // When overlap X is processed, the tree holds exactly the overlaps before X in
    // first order; among them, those with second < X.second also have first < X.first,
    // so they are the legal predecessors.
    std::sort(ov.begin(), ov.end(), [](overlap const& a, overlap const& b) {
        return a.first != b.first ? a.first < b.first : a.second > b.second;
    });

    // Coordinate-compress the second coordinates.
    std::vector<uint64_t> ys;
    ys.reserve(n);
    for (overlap const& o : ov) ys.push_back(o.second);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    auto rank = [&ys](uint64_t y) -> size_t {
        return static_cast<size_t>(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin());
    };
    size_t const S = ys.size();

    // Fenwick (1-indexed) prefix-maximum. A cell records the best chain ending at a
    // processed overlap: longer wins; ties go to the smaller (first, second) overlap,
    // which reproduces the chain tie-breaking the unit tests pin.
    struct Cell { uint64_t length; overlap end; };
    Cell const empty {0, null_overlap};
    auto better = [](Cell const& a, Cell const& b) -> bool {
        if (a.length != b.length) return a.length > b.length;
        if (a.length == 0)        return false;
        if (a.end.first != b.end.first) return a.end.first < b.end.first;
        return a.end.second < b.end.second;
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

    std::unordered_map<overlap, overlap> previous;
    std::unordered_map<overlap, uint64_t> length;
    for (overlap const& o : ov)
    {
        size_t const r = rank(o.second);                          // ranks [0, r-1] have second < o.second
        Cell const best = (r == 0) ? empty : bit_prefix_max(r);   // Fenwick positions [1, r]
        uint64_t const len = best.length + 1;
        previous[o] = (best.length == 0) ? null_overlap : best.end;
        length[o]   = len;
        bit_update(r + 1, Cell {len, o});
    }

    // Chain end: the longest chain, ties broken by larger (first, second).
    overlap chain_end {null_overlap};
    uint64_t max_len {0};
    for (overlap const& o : ov)
    {
        uint64_t const len = length[o];
        if (len > max_len ||
            (len == max_len && (o.first > chain_end.first ||
                                (o.first == chain_end.first && o.second > chain_end.second))))
        {
            max_len = len;
            chain_end = o;
        }
    }

    std::vector<overlap> chain(max_len);
    overlap cur {chain_end};
    for (uint64_t i = max_len; i > 0; i--)
    {
        chain[i - 1] = cur;
        cur = previous[cur];
    }
    return chain;
}

}} // namespace sorting // namespace km
