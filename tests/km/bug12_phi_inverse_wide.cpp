// Regression tests for the inverse of the minimizer permutation phi at 2m > 192 (m >= 97).
//
// phi_inv undoes phi's two odd-constant multiplications with their inverses mod 2^bits(kuint),
// computed by mul_inverse's Newton iteration. It ran a fixed 6 steps: x = a is exact on 3 low bits
// and each step doubles that, so 3 * 2^6 = 192 exact bits, enough up to __uint128_t but not for the
// 256-bit kuint256 word. phi_inv(phi(x)) != x as soon as the 2m-bit minimizer exceeds 192 bits, so
// unpermute_minimizer_slot decoded wrong minimizers at m >= 97. The binary index, the query and the
// set operations never invert phi; the one production caller is the ASCII dump (`sskm construct
// --ascii`), which printed wrong minimizer nucleotides: at k=110, m=100 none of the dumped k-mers
// was in the input. mul_inverse now iterates until the exact bits cover the word.

#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <io/Skmer.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace {

const char NUC[4] = {'A', 'C', 'T', 'G'};

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = (c == 'A') ? 'T' : (c == 'T') ? 'A' : (c == 'C') ? 'G' : 'C';
    return r;
}

std::string canon(const std::string& s) { return std::min(s, revcomp(s)); }

template <typename kuint>
kuint random_word(std::mt19937_64& rng) {
    if constexpr (sizeof(kuint) <= 8) {
        return static_cast<kuint>(rng());
    } else {
        kuint x {0};
        for (size_t i {0}; i < sizeof(kuint); i += 8) x = static_cast<kuint>(x << 64) | static_cast<kuint>(rng());
        return x;
    }
}

template <typename kuint>
void expect_exact_inverses() {
    std::mt19937_64 rng(sizeof(kuint));
    std::vector<kuint> odd {static_cast<kuint>(0xbf58476d1ce4e5b9ULL), static_cast<kuint>(0x94d049bb133111ebULL)};
    for (int i {0}; i < 1000; ++i) odd.push_back(random_word<kuint>(rng) | kuint{1});
    for (const kuint a : odd)
        ASSERT_TRUE(static_cast<kuint>(a * km::SkmerManipulator<kuint>::mul_inverse(a)) == kuint{1})
            << "mul_inverse not exact at " << sizeof(kuint) * 8 << " bits";
}

// phi_inv(phi(x)) == x and the permute/unpermute slot round trip, on random 2m-bit minimizers.
template <typename kuint>
void expect_round_trips(uint64_t k, uint64_t m) {
    ASSERT_EQ(km::sortedlist::select_generation_width_bytes(k, m), sizeof(kuint)) << "k=" << k << " m=" << m;
    km::SkmerManipulator<kuint> manip {k, m};
    std::mt19937_64 rng(k * 1000 + m);
    const kuint mask {(2 * m >= sizeof(kuint) * 8) ? static_cast<kuint>(~kuint{0})
                                                   : static_cast<kuint>((kuint{1} << (2 * m)) - 1)};
    using kpair = typename km::Skmer<kuint>::pair;
    const uint64_t shift {4 * (k - m)};
    for (int i {0}; i < 2000; ++i) {
        const kuint x {static_cast<kuint>(random_word<kuint>(rng) & mask)};
        ASSERT_TRUE(manip.phi_inv(manip.phi(x)) == x) << "phi_inv o phi != id, k=" << k << " m=" << m;
        // a super-k-mer whose slot holds x, with random flanks
        km::Skmer<kuint> sk {(kpair {x} << shift) | (kpair {random_word<kuint>(rng), random_word<kuint>(rng)}
                                                     & (kpair {static_cast<kuint>(~kuint{0}), static_cast<kuint>(~kuint{0})}
                                                        >> (2 * sizeof(kuint) * 8 - shift)))};
        const km::Skmer<kuint> before {sk};
        manip.permute_minimizer_slot(sk);
        manip.unpermute_minimizer_slot(sk);
        ASSERT_TRUE(sk.m_pair == before.m_pair) << "slot round trip, k=" << k << " m=" << m;
    }
}

} // namespace

TEST(PhiInverseWide, MulInverseIsExactAtEveryWidth) {
    expect_exact_inverses<uint32_t>();
    expect_exact_inverses<uint64_t>();
    expect_exact_inverses<__uint128_t>();
    expect_exact_inverses<km::kuint256>();
}

TEST(PhiInverseWide, RoundTripsUpToTheFullWord) {
    expect_round_trips<uint32_t>(20, 16);
    expect_round_trips<uint64_t>(40, 32);
    expect_round_trips<__uint128_t>(80, 64);
    for (const auto& [k, m] : {std::pair<uint64_t, uint64_t>{70, 65}, {100, 96}, {100, 97}, {110, 100},
                               {127, 120}, {130, 127}, {150, 128}, {128, 128}})
        expect_round_trips<km::kuint256>(k, m);
}

// The ASCII dump (in-RAM construction) decodes back to exactly the input's k-mers.
TEST(PhiInverseWide, AsciiDumpDecodesTheInputKmers) {
    std::mt19937_64 rng(97);
    std::string seq(400, 'A');
    for (char& c : seq) c = NUC[rng() & 3u];
    const std::string fa {::testing::TempDir() + "bug12.fa"}, out {::testing::TempDir() + "bug12.txt"};
    { std::ofstream f(fa); f << ">s\n" << seq << "\n"; }
    for (const auto& [k, m] : {std::pair<uint64_t, uint64_t>{100, 96}, {100, 97}, {110, 100}, {127, 120},
                               {150, 128}, {128, 128}, {40, 20}}) {
        std::set<std::string> input;
        for (uint64_t i {0}; i + k <= seq.size(); ++i) input.insert(canon(seq.substr(i, k)));
        km::sortedlist::SortedListBuildParams p;
        p.k = k;
        p.m = m;
        p.input_path = fa;
        p.output_path = out;
        p.ascii = true;
        km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m), [&]<typename gen>() {
            km::sortedlist::build_sorted_list<gen, gen>(p, 0);
        });
        std::ifstream in(out);
        std::string line;
        std::getline(in, line);
        std::set<std::string> dumped;
        std::string nucl, pref, suff;
        while (in >> nucl >> pref >> suff)
            for (uint64_t i {0}; i + k <= nucl.size(); ++i) dumped.insert(canon(nucl.substr(i, k)));
        EXPECT_EQ(dumped, input) << "k=" << k << " m=" << m;
    }
}
