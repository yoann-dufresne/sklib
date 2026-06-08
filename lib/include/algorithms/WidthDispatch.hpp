#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <io/wide_int.hpp>

#ifndef SKLIB_WIDTH_DISPATCH_HPP
#define SKLIB_WIDTH_DISPATCH_HPP

// Runtime selection of the packed-integer width for a sorted skmer list.
//
// A skmer is stored in a pair of `kuint` values, so its bit capacity is
// `16 * sizeof(kuint)` (two words). A full super-k-mer needs `2*(2k-m)` bits; a
// quotiented one (the top `b` φ-minimizer bits live in the bucket id) needs
// `2*(2k-m) - b`. We precompile the data structures for uint32_t / uint64_t /
// __uint128_t / kuint256 (pair capacities 64 / 128 / 256 / 512 bits) and pick,
// per requested k/m/b, the smallest type that fits. dispatch_width_bytes() turns
// the chosen byte width back into a concrete template instantiation at the call site.

namespace km
{
namespace sortedlist
{

// Bit capacity of a Skmer<kuint> pair (two kuint words).
constexpr uint64_t pair_capacity_bits(uint64_t kuint_bytes) { return 16 * kuint_bytes; }

// Smallest precompiled kuint width (in bytes: 4, 8, 16, or 32) whose pair holds
// `need_bits`. Throws if even kuint256's 512-bit pair is too small.
inline uint64_t select_width_bytes(uint64_t need_bits) {
    if (need_bits <= pair_capacity_bits(4))  return 4;   // pair<uint32_t>     (64 bits)
    if (need_bits <= pair_capacity_bits(8))  return 8;   // pair<uint64_t>     (128 bits)
    if (need_bits <= pair_capacity_bits(16)) return 16;  // pair<__uint128_t>  (256 bits)
    if (need_bits <= pair_capacity_bits(32)) return 32;  // pair<kuint256>     (512 bits)
    throw std::runtime_error(
        "k too large for the compiled integer widths: needs " + std::to_string(need_bits) +
        " bits but the widest backend (kuint256) holds 512. Reduce k, or raise m.");
}

// Number of high φ-minimizer bits a power-of-two prefix bucketing removes, i.e. the
// quotient bit count `b`. Mirrors make_prefix_bucketing exactly:
// effective_bits = min(floor(log2(buckets)), 2m).
inline uint64_t effective_bucket_bits(uint64_t m, uint64_t requested_buckets) {
    uint64_t bucket_bits = 0;
    while ((uint64_t{1} << (bucket_bits + 1)) <= std::max<uint64_t>(requested_buckets, 1))
        bucket_bits++;
    return std::min<uint64_t>(bucket_bits, 2 * m);
}

// Invoke `f.template operator()<kuint>()` with kuint = the type of the given byte width.
// `f` is a generic functor whose operator() is templated on the integer type; all branches
// must yield the same return type (typically void or int).
template<typename F>
decltype(auto) dispatch_width_bytes(uint64_t width_bytes, F&& f) {
    switch (width_bytes) {
        case 4:  return f.template operator()<uint32_t>();
        case 8:  return f.template operator()<uint64_t>();
        case 16: return f.template operator()<__uint128_t>();
        case 32: return f.template operator()<km::kuint256>();
        default:
            throw std::runtime_error("unsupported kuint width (bytes): " + std::to_string(width_bytes));
    }
}

} // namespace sortedlist
} // namespace km

#endif // SKLIB_WIDTH_DISPATCH_HPP
