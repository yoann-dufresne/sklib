#include <gtest/gtest.h>
#include <cstdint>
#include <vector>
#include <random>
#include <io/Skmer.hpp>

using km::Skmer;
using km::SkmerManipulator;

// φ must be a bijection on the 2m-bit minimizer space, invertible, with φ(0) != 0
// (otherwise poly-A = 0 would stay the global minimum and the reordering would be void).
TEST(MinimizerHash, BijectionFullEnumerationSmallM) {
    for (uint64_t m : std::vector<uint64_t>{1, 2, 3, 5, 7, 10, 11}) {
        const uint64_t k = m + 5;
        SkmerManipulator<uint64_t> manip(k, m);
        const uint64_t N = uint64_t{1} << (2 * m);
        std::vector<char> seen(N, 0);
        for (uint64_t x = 0; x < N; x++) {
            const uint64_t y = manip.phi(x);
            ASSERT_LT(y, N) << "phi output outside 2m-bit range (m=" << m << ")";
            ASSERT_EQ(seen[y], 0) << "phi not injective (m=" << m << ", x=" << x << ")";
            seen[y] = 1;
            ASSERT_EQ(manip.phi_inv(y), x) << "phi_inv o phi != id (m=" << m << ")";
        }
        EXPECT_NE(manip.phi(0), 0u) << "phi(0)==0 keeps poly-A minimal (m=" << m << ")";
    }
}

// For larger m full enumeration is too big: sample and check invertibility + range.
TEST(MinimizerHash, InvertibleSampledLargeM) {
    std::mt19937_64 rng(12345);
    for (uint64_t m : std::vector<uint64_t>{15, 20, 25, 30}) {
        const uint64_t k = (m + 2 <= 32) ? m + 2 : 32;
        SkmerManipulator<uint64_t> manip(k, m);
        const uint64_t mask =
            (2 * m >= 64) ? ~uint64_t{0} : ((uint64_t{1} << (2 * m)) - 1);
        for (int i = 0; i < 200000; i++) {
            const uint64_t x = rng() & mask;
            const uint64_t y = manip.phi(x);
            ASSERT_LE(y, mask);
            ASSERT_EQ(manip.phi_inv(y), x) << "phi_inv o phi != id (m=" << m << ")";
        }
        EXPECT_NE(manip.phi(0), 0u);
    }
}

// permute/unpermute_minimizer_slot must round-trip and leave the flank bits untouched.
TEST(MinimizerHash, PermuteSlotRoundTrip) {
    const uint64_t k = 12, m = 7;                 // k-m=5 -> slot at bit 20, 2m=14 bits
    SkmerManipulator<uint64_t> manip(k, m);
    const uint64_t shift = 4 * (k - m);           // 20
    const uint64_t mmask = (uint64_t{1} << (2 * m)) - 1; // minimizer 14 bits
    const uint64_t fmask = (uint64_t{1} << shift) - 1;   // flanks 20 bits
    std::mt19937_64 rng(999);
    for (int i = 0; i < 100000; i++) {
        const uint64_t mv = rng() & mmask;
        const uint64_t fv = rng() & fmask;
        Skmer<uint64_t> sk{typename Skmer<uint64_t>::pair((mv << shift) | fv)};
        manip.permute_minimizer_slot(sk);
        const uint64_t got = static_cast<uint64_t>(sk.m_pair);
        EXPECT_EQ((got >> shift) & mmask, manip.phi(mv));
        EXPECT_EQ(got & fmask, fv) << "flanks must be untouched";
        manip.unpermute_minimizer_slot(sk);
        EXPECT_EQ(static_cast<uint64_t>(sk.m_pair), (mv << shift) | fv);
    }
}
