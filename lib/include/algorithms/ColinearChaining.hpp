#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdint.h>

#ifndef COLINEARCHAINING_H
#define COLINEARCHAINING_H

namespace std {
    template <>
    struct hash<std::pair<unsigned long, unsigned long>> {
        std::size_t operator()(const std::pair<unsigned long, unsigned long>& p) const noexcept {
            return std::hash<unsigned long>{}(p.first) ^ (std::hash<unsigned long>{}(p.second) << 1);
        }
    };
}

// --- Colienar chaining algorithm and data structures---
namespace km
{
namespace chaining
{

using overlap = std::pair<uint64_t, uint64_t>;

/** Colinear chaining algorithm to select a compatible set of overlaps. The overlaps are compatible if their 
 * coordinates are not crossing and if they do not have common coordinates.
 * WARNING: The algorithme will change the order of the input vector.
 * 
 * @param begin iterator to the first element of the vector
 * @param end iterator to the last element of the vector
 * 
 * @return a vector of overlaps containing compatible overlaps
 **/
std::vector<overlap> colinear_chaining(std::vector<overlap>::iterator begin, std::vector<overlap>::iterator end);


/** Maximum chain of overlaps with STRICTLY increasing first AND second coordinate, via patience
 * sorting (longest strictly-increasing subsequence with reconstruction). Same optimization target as
 * colinear_chaining — a maximum chain, so the resulting super-k-mers are just as compact — but with
 * smaller constants: a single sort + one binary search per overlap, with no coordinate compression and
 * no Fenwick tree. The tie-breaking among equally-long chains may differ from colinear_chaining, so the
 * super-k-mer packing can differ (same k-mer SET); used on the set-operation re-compaction path, whose
 * correctness is the represented set, not the byte layout. WARNING: reorders the input range.
 */
inline std::vector<overlap> greedy_chaining(std::vector<overlap>::iterator begin,
                                            std::vector<overlap>::iterator end)
{
    const size_t n = static_cast<size_t>(end - begin);
    if (n == 0) return {};
    // (first asc, second desc): same-first overlaps then have descending second, so the strictly
    // increasing LIS on second can never chain two overlaps that share a first coordinate.
    std::sort(begin, end, [](const overlap& a, const overlap& b) {
        return a.first != b.first ? a.first < b.first : a.second > b.second;
    });
    std::vector<size_t> parent(n, SIZE_MAX);
    std::vector<size_t> tail;   // tail[L] = element index = current minimal tail of a length-(L+1) chain
    for (size_t i = 0; i < n; ++i) {
        const uint64_t s = begin[i].second;
        // first pile whose tail's second >= s (strictly increasing => lower_bound on second)
        size_t lo = 0, hi = tail.size();
        while (lo < hi) { const size_t mid = (lo + hi) / 2; if (begin[tail[mid]].second < s) lo = mid + 1; else hi = mid; }
        if (lo > 0) parent[i] = tail[lo - 1];
        if (lo == tail.size()) tail.push_back(i); else tail[lo] = i;
    }
    std::vector<overlap> chain;
    chain.reserve(tail.size());
    for (size_t cur = tail.back(); ; cur = parent[cur]) {
        chain.push_back(begin[cur]);
        if (parent[cur] == SIZE_MAX) break;
    }
    std::reverse(chain.begin(), chain.end());   // backtracked tail->head => restore ascending order
    return chain;
}


} // namespace chaining
} // namespace km

#endif
