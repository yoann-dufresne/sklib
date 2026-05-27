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


} // namespace chaining
} // namespace km

#endif
