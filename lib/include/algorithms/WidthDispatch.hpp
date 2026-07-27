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

// Generation (work) width for a given k/m, in bytes.
//
// Two independent constraints, and only the first used to be applied:
//
//  1. the interleaved super-k-mer must fit the PAIR: 2*(2k-m) <= 16*bytes;
//  2. the minimizer must fit a SINGLE WORD: 2m <= 8*bytes.
//
// (2) matters because the whole minimizer pipeline works on one word: `minimizer()` returns the
// pair's low word (`to_kuint()`), `phi`/`reverse_2m` mix a `kuint`, `permute_minimizer_slot` writes
// a `kuint`-wide ψ back into the 2m-bit slot, and `mmer_repeats` rolls the central m-mer in a
// `kuint`. Sizing on (1) alone lets 2m exceed the word whenever 3m > 2k, and then every one of
// those silently truncates. The damage is graded and was measured on ecoli:
//
//   * the top 2m - 8*bytes bits of every stored ψ slot are zero, so the minimizer-prefix bucketing
//     loses that many bits (k=40,m=33: 2 bits, 1024 of 4096 buckets actually used);
//   * once the WHOLE bucket prefix falls in that zero region (2m - b >= 8*bytes) the router
//     computed `mini >> shift` with shift >= width — undefined behaviour, and on x86 the masked
//     shift count returned the minimizer itself as a bucket id, indexing the per-bucket buffer
//     array out of bounds: SIGSEGV for e.g. `construct -k 42 -m 38`;
//   * and `mmer_repeats` rolling a 2m-bit m-mer in a narrower word misdetects the ambiguous
//     minimizer case, which changes the per-k-mer framing and breaks EXACTNESS — k=25,m=22 (a
//     44-bit minimizer in a uint32 word) answered 3 false negatives on ecoli's 4.5 M k-mers and
//     3 875 false positives out of 100 000 random k-mers.
//
// Taking the max of the two constraints removes all three at once. It cannot affect any k/m where
// the minimizer already fitted: 2m <= 8*bytes implies 4m <= 16*bytes, so select_width_bytes(4m)
// is then <= the width constraint (1) already picked, and the max is that same width.
inline uint64_t select_generation_width_bytes(uint64_t k, uint64_t m) {
    const uint64_t for_skmer {select_width_bytes(2 * (2 * k - m))};
    // pair capacity 16*bytes >= 4m  <=>  word capacity 8*bytes >= 2m
    const uint64_t for_minimizer {select_width_bytes(4 * m)};
    return std::max(for_skmer, for_minimizer);
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
