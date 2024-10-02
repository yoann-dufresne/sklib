#ifndef BIT_NOARCH_HPP
#define BIT_NOARCH_HPP

#include <cassert>
#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace dbglib::bit {

uint64_t pdep_u64(uint64_t a, uint64_t mask);

template <typename T>
constexpr int countr_zero(T v) noexcept {
    static_assert(std::is_unsigned<T>::value, "Value must be an unsigned integer type");
    return __builtin_ctz(v);
}

#if defined(__cplusplus) && (__cplusplus >= 202002L)
// FIXME downgrading project requirements to C++17 breaks skmer class.
// template <typename T>
// constexpr int std_popcount(T x) {return std::popcount(x);}
constexpr int std_popcount(uint8_t x) {return __builtin_popcount(x);}
constexpr int std_popcount(uint16_t x) {return __builtin_popcount(x);}
constexpr int std_popcount(uint32_t x) {return __builtin_popcount(x);}
constexpr int std_popcount(uint64_t x) {return __builtin_popcountll(x);}
#elif defined(__SSE4_2__)
inline int sse4_popcount(uint8_t x) {return _mm_popcnt_u32(x);}
inline int sse4_popcount(uint16_t x) {return _mm_popcnt_u32(x);}
inline int sse4_popcount(uint32_t x) {return _mm_popcnt_u32(x);}
inline int sse4_popcount(uint64_t x) {return _mm_popcnt_u64(x);}
#elif defined(__clang__) // use same extension as gcc
constexpr int clang_popcount(uint8_t x) {return __builtin_popcount(x);}
constexpr int clang_popcount(uint16_t x) {return __builtin_popcount(x);}
constexpr int clang_popcount(uint32_t x) {return __builtin_popcount(x);}
constexpr int clang_popcount(uint64_t x) {return __builtin_popcountll(x);}
#elif defined(__GNUC__) || defined(__GNUG__)
constexpr int gcc_popcount(uint8_t x) {return __builtin_popcount(x);}
constexpr int gcc_popcount(uint16_t x) {return __builtin_popcount(x);}
constexpr int gcc_popcount(uint32_t x) {return __builtin_popcount(x);}
constexpr int gcc_popcount(uint64_t x) {return __builtin_popcountll(x);}
#elif defined(_MSC_VER)
constexpr int mvsc_popcount(uint8_t x) {return __popcnt16(x);}
constexpr int mvsc_popcount(uint16_t x) {return __popcnt16(x);}
constexpr int mvsc_popcount(uint32_t x) {return __popcnt(x);}
constexpr int mvsc_popcount(uint64_t x) {return __popcnt64(x);}
#else
template <typename T>
constexpr int emulated_popcount(T x)
{
    x = x - ((x >> 1) & (T)~(T)0/3);
    x = (x & (T)~(T)0/15*3) + ((x >> 2) & (T)~(T)0/15*3);
    x = (x + (x >> 4)) & (T)~(T)0/255*15;
    return (T)(x * ((T)~(T)0/255)) >> (sizeof(T) - 1) * 8;
}
#endif

template <typename T>
constexpr int popcount(T x) noexcept
{
    #if __cplusplus >= 202002L
        return std_popcount(x);
    #elif __SSE4_2__
        return sse4_popcount(x);
    #elif defined(__clang__) // use same extension as gcc
        return clang_popcount(x);
    #elif defined(__GNUC__) || defined(__GNUG__)
        return gcc_popcount(x);
    #elif defined(_MSC_VER)
        return mvsc_popcount(x);
    #else
        return emulated_popcount(x);
    #endif
}

inline int select1(uint64_t x, std::size_t th)
{
#ifndef __BMI2__
    // Modified from: Bit Twiddling Hacks
    // https://graphics.stanford.edu/~seander/bithacks.html#SelectPosFromMSBRank
    unsigned int s;       // Output: Resulting position of bit with rank r [1-64]
    uint64_t a, b, c, d;  // Intermediate temporaries for bit count.
    unsigned int t;       // Bit count temporary.
    auto pc = popcount(x);
    assert(pc >= 0);
    assert(th < static_cast<std::size_t>(pc));
    th = pc - th;

    a = x - ((x >> 1) & ~0UL / 3);
    b = (a & ~0UL / 5) + ((a >> 2) & ~0UL / 5);
    c = (b + (b >> 4)) & ~0UL / 0x11;
    d = (c + (c >> 8)) & ~0UL / 0x101;
    t = (d >> 32) + (d >> 48);
    s = 64;
    s -= ((t - th) & 256) >> 3;
    th -= (t & ((t - th) >> 8));
    t = (d >> (s - 16)) & 0xff;
    s -= ((t - th) & 256) >> 4;
    th -= (t & ((t - th) >> 8));
    t = (c >> (s - 8)) & 0xf;
    s -= ((t - th) & 256) >> 5;
    th -= (t & ((t - th) >> 8));
    t = (b >> (s - 4)) & 0x7;
    s -= ((t - th) & 256) >> 6;
    th -= (t & ((t - th) >> 8));
    t = (a >> (s - 2)) & 0x3;
    s -= ((t - th) & 256) >> 7;
    th -= (t & ((t - th) >> 8));
    t = (x >> (s - 1)) & 0x1;
    s -= ((t - th) & 256) >> 8;
    return s - 1;
#else
    uint64_t i = 1ULL << th;
    asm("pdep %[x], %[mask], %[x]" : [x] "+r"(x) : [mask] "r"(i));
    asm("tzcnt %[bit], %[index]" : [index] "=r"(i) : [bit] "g"(x) : "cc");
    return i;
#endif
}

} // namespace dbglib::bit

#endif // BIT_NOARCH_HPP