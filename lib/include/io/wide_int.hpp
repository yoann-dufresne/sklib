#ifndef SKLIB_WIDE_INT_HPP
#define SKLIB_WIDE_INT_HPP

// 256-bit packed-integer backend for the widest sorted skmer lists.
//
// A skmer is stored in a pair of `kuint` words, so its bit capacity is `16 * sizeof(kuint)`.
// The native backends top out at __uint128_t (256-bit pair); to index larger k-mers we need a
// 256-bit kuint, whose pair holds 512 bits -> supports `2*(2k-m) <= 512`. There is no native
// 256-bit integer, so we use the compiler's extended bit-precise integer `_BitInt(256)`. It
// supports every operator the kmer/skmer code applies (shifts, &|^~, +-*, comparisons) natively,
// so it drops straight into the templated `Skmer<kuint>` / `SkmerManipulator<kuint>` machinery.
//
// Requires a compiler with C++ _BitInt(256) support: Clang >= 16, or GCC >= 15 (GCC added _BitInt
// for C in 14 but for C++ only in 15; g++ <= 14 does not recognise the keyword in C++ at all). On
// Ubuntu 24.04, where GCC tops out at 14, build with Clang (`-DCMAKE_CXX_COMPILER=clang++-18`).
// Guard on the compiler version so an older toolchain fails with a clear message instead of a
// cryptic parse error. (We avoid keying on __BITINT_MAXWIDTH__: it is not reliably predefined in
// C++ mode, so testing it could misfire on a compiler that actually supports _BitInt.)

#if defined(__clang__)
#  if (__clang_major__ < 16)
#    error "sklib's 256-bit backend needs C++ _BitInt(256): Clang >= 16."
#  endif
#elif defined(__GNUC__)
#  if (__GNUC__ < 15)
#    error "sklib's 256-bit backend needs C++ _BitInt(256): GCC >= 15 or Clang >= 16. g++ <= 14 lacks C++ _BitInt; on Ubuntu 24.04 build with Clang, e.g. -DCMAKE_CXX_COMPILER=clang++-18."
#  endif
#endif

namespace km
{

// Widest precompiled packed-integer word. Its Skmer<kuint256>::pair spans 512 bits.
using kuint256 = unsigned _BitInt(256);

} // namespace km

#endif // SKLIB_WIDE_INT_HPP
