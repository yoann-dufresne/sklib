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
template <typename T>
constexpr int std_popcount(T x) {return std::popcount(x);}
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

} // namespace dbglib::bit

#endif // BIT_NOARCH_HPP