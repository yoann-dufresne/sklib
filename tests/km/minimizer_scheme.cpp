// Unit tests for the minimizer scheme (io/MinimizerScheme.hpp) and its hook into SkmerManipulator.
//
//  1. OcClassifier::class_of against a brute-force, string-level definition of the class.
//  2. Strand symmetry: class(x) == class(rc x) for every valid (m, s); and the wrong parity of s is
//     really asymmetric, so the check can detect the bug it guards against.
//  3. OcStream (the streaming class the Skmerator uses) equals class_of at EVERY position, on random
//     sequence and on poly-A / short tandem repeats, read with the Skmerator's lag.
//  4. The random rank is bit-identical to φ; the open-closed rank is (class << 2m) | φ and stays
//     injective on canonical m-mers.
//  5. SchemeParams: default_s, validate (range + parity), header decoding, widths.

#include <gtest/gtest.h>

#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <io/MinimizerScheme.hpp>
#include <io/Skmer.hpp>
#include <algorithms/WidthDispatch.hpp>

#include "scheme_test_utils.hpp"

namespace {

std::vector<uint8_t> random_nucls(size_t n, std::mt19937_64& rng) {
    std::vector<uint8_t> v(n);
    for (auto& x : v) x = static_cast<uint8_t>(rng() & 3u);
    return v;
}

std::vector<uint8_t> rc_nucls(const std::vector<uint8_t>& x) {
    std::vector<uint8_t> r(x.rbegin(), x.rend());
    for (auto& v : r) v = static_cast<uint8_t>(v ^ 2u);   // A=0<->T=2, C=1<->G=3
    return r;
}

// Brute-force class from the definition, for ANY (m, s, offset): canonical s-mer = the smaller of
// the s-mer and its reverse complement as linear codes (first nucleotide most significant), key =
// π_s(canonical) taken from a classifier with the same s, M = min over the m-s+1 s-mers, then the
// "attains the minimum" tests at `open_offset` (open) and at the two ends (closed).
uint8_t brute_class(const std::vector<uint8_t>& x, uint64_t s, uint64_t open_offset,
                    const km::OcClassifier& pi_source) {
    const uint64_t m = x.size();
    std::vector<uint64_t> keys;
    for (uint64_t i = 0; i + s <= m; ++i) {
        std::vector<uint8_t> u(x.begin() + i, x.begin() + i + s);
        const std::vector<uint8_t> r = rc_nucls(u);
        auto code = [](const std::vector<uint8_t>& v) {
            uint64_t c = 0;
            for (uint8_t n : v) c = (c << 2) | n;
            return c;
        };
        keys.push_back(pi_source.pi(std::min(code(u), code(r))));
    }
    const uint64_t M = *std::min_element(keys.begin(), keys.end());
    if (keys[open_offset] == M) return 0;
    if (keys.front() == M || keys.back() == M) return 1;
    return 2;
}

// Valid (m, s) pairs for a grid: every s of m's parity in [1, min(m-2, smax)].
std::vector<std::pair<uint64_t, uint64_t>> valid_pairs(uint64_t mmin, uint64_t mmax, uint64_t smax) {
    std::vector<std::pair<uint64_t, uint64_t>> out;
    for (uint64_t m = mmin; m <= mmax; ++m)
        for (uint64_t s = (m % 2 == 1 ? 1 : 2); s + 2 <= m && s <= smax; s += 2)
            out.push_back({m, s});
    return out;
}

std::vector<std::vector<uint8_t>> low_complexity_mmers(uint64_t m) {
    std::vector<std::vector<uint8_t>> v;
    const std::vector<std::vector<uint8_t>> units = {{0}, {1}, {2}, {3}, {0, 1}, {1, 3}, {0, 2},
                                                     {0, 0, 1}, {0, 1, 3}, {2, 2, 2, 1}, {0, 1, 2, 3}};
    for (const auto& u : units) {
        std::vector<uint8_t> x;
        while (x.size() < m) x.push_back(u[x.size() % u.size()]);
        v.push_back(x);
    }
    return v;
}

} // namespace

// ---- 1. class_of == brute force ----------------------------------------------------------------
TEST(MinimizerScheme, ClassOfMatchesBruteForce) {
    std::mt19937_64 rng(101);
    for (const auto& [m, s] : valid_pairs(3, 35, 11)) {
        const km::OcClassifier c(m, s);
        ASSERT_EQ(c.open_offset(), (m - s) / 2);
        std::vector<std::vector<uint8_t>> xs = low_complexity_mmers(m);
        for (int r = 0; r < 300; ++r) xs.push_back(random_nucls(m, rng));
        for (const auto& x : xs)
            ASSERT_EQ(c.class_of(x.data()), brute_class(x, s, (m - s) / 2, c)) << "m=" << m << " s=" << s;
    }
}

// The precomputed key table (s <= 8) and the per-call key agree on every s-mer, and key_at matches.
TEST(MinimizerScheme, KeyTableMatchesKeyFunction) {
    for (uint64_t s = 1; s <= km::OcClassifier::TABLE_MAX_S; ++s) {
        const km::OcClassifier c(s + 2, s);
        ASSERT_TRUE(c.has_table());
        std::set<uint64_t> canon_keys;
        for (uint64_t u = 0; u < (uint64_t{1} << (2 * s)); ++u) {
            const uint64_t rc = km::OcClassifier::rc_code(u, s);
            ASSERT_EQ(c.key_fast(u, rc), c.key_from_codes(u, rc)) << "s=" << s << " u=" << u;
            ASSERT_EQ(c.key_from_codes(u, rc), c.key_from_codes(rc, u));   // canonical
            canon_keys.insert(c.key_from_codes(u, rc));
            std::vector<uint8_t> nt(s);
            for (uint64_t j = 0; j < s; ++j) nt[j] = static_cast<uint8_t>((u >> (2 * (s - 1 - j))) & 3u);
            ASSERT_EQ(c.key_at(nt.data()), c.key_from_codes(u, rc));
        }
        // π_s is a bijection, so distinct canonical s-mers get distinct keys: 4^s/2 (+ palindromes).
        const uint64_t n = uint64_t{1} << (2 * s);
        const uint64_t palins = (s % 2 == 0) ? (uint64_t{1} << s) : 0;
        EXPECT_EQ(canon_keys.size(), (n + palins) / 2) << "s=" << s;
    }
    EXPECT_FALSE(km::OcClassifier(13, 9).has_table());
}

// ---- 2. strand symmetry ------------------------------------------------------------------------
TEST(MinimizerScheme, ClassIsStrandSymmetric) {
    std::mt19937_64 rng(202);
    for (const auto& [m, s] : valid_pairs(3, 64, km::SchemeParams::MAX_OC_S)) {
        const km::OcClassifier c(m, s);
        std::vector<std::vector<uint8_t>> xs = low_complexity_mmers(m);
        for (int r = 0; r < 30; ++r) xs.push_back(random_nucls(m, rng));
        for (const auto& x : xs) {
            const std::vector<uint8_t> rx = rc_nucls(x);
            ASSERT_EQ(c.class_of(x.data()), c.class_of(rx.data())) << "m=" << m << " s=" << s;
        }
    }
}

// The parity rule is load-bearing: with m - s odd the "open" offset floor((m-s)/2) has no mirror
// image, so the class differs between strands for a large fraction of m-mers (the bug the
// validation prevents). Same-parity s on the same m-mers: zero violations.
TEST(MinimizerScheme, WrongParityIsAsymmetric) {
    std::mt19937_64 rng(303);
    for (const auto& [m, s_bad, s_good] : std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>{
             {15, 4, 3}, {15, 2, 5}, {11, 4, 5}, {16, 5, 4}, {31, 4, 5}, {21, 6, 5}}) {
        const km::OcClassifier bad_pi(s_bad + 2, s_bad);    // only its π_s is used (depends on s only)
        const km::OcClassifier good(m, s_good);
        uint64_t bad_viol = 0, good_viol = 0;
        const int N = 2000;
        for (int r = 0; r < N; ++r) {
            const std::vector<uint8_t> x = random_nucls(m, rng);
            const std::vector<uint8_t> rx = rc_nucls(x);
            bad_viol  += brute_class(x, s_bad, (m - s_bad) / 2, bad_pi) != brute_class(rx, s_bad, (m - s_bad) / 2, bad_pi);
            good_viol += good.class_of(x.data()) != good.class_of(rx.data());
        }
        EXPECT_EQ(good_viol, 0u) << "m=" << m << " s=" << s_good;
        EXPECT_GT(bad_viol, static_cast<uint64_t>(N / 20)) << "m=" << m << " s=" << s_bad
            << ": wrong-parity class unexpectedly symmetric; the parity test would not catch the bug";
    }
}

// ---- 3. streaming == reference at every position ----------------------------------------------
TEST(MinimizerScheme, StreamMatchesReferenceEverywhere) {
    std::mt19937_64 rng(404);
    std::vector<std::vector<uint8_t>> seqs;
    for (int r = 0; r < 2; ++r) seqs.push_back(random_nucls(1500, rng));
    // poly-A, poly-C, short tandem repeats (period 1..6), and repeats embedded in random sequence.
    for (const std::vector<uint8_t>& unit : std::vector<std::vector<uint8_t>>{
             {0}, {1}, {0, 1}, {1, 3}, {0, 0, 1}, {0, 1, 3, 2}, {2, 0, 0, 1, 3}, {3, 3, 1, 0, 2, 1}}) {
        std::vector<uint8_t> v = random_nucls(100, rng);
        for (int i = 0; i < 200; ++i) v.push_back(unit[i % unit.size()]);
        const std::vector<uint8_t> tail = random_nucls(100, rng);
        v.insert(v.end(), tail.begin(), tail.end());
        for (int i = 0; i < 150; ++i) v.push_back(unit[i % unit.size()]);
        seqs.push_back(v);
    }
    // Every s parity/size regime: s = 1, table (s <= 8) and no-table (s > 8) keys, even and odd m.
    std::vector<std::pair<uint64_t, uint64_t>> grid;
    for (const uint64_t m : {3ull, 4ull, 5ull, 6ull, 7ull, 8ull, 9ull, 11ull, 15ull, 16ull, 19ull, 21ull, 31ull, 32ull})
        for (uint64_t s = (m % 2 == 1 ? 1 : 2); s + 2 <= m && s <= 11; s += 2) grid.push_back({m, s});
    for (const auto& [m, s] : grid) {
        const km::OcClassifier c(m, s);
        // Read the class back with every lag the Skmerator can use (0 .. k-m), here k-m = 9.
        const uint64_t lag_max = 9;
        km::OcStream st;
        st.init(&c, lag_max + 1);
        for (const auto& seq : seqs) {
            st.reset();
            for (uint64_t i = 0; i < seq.size(); ++i) {
                st.push(seq[i]);
                if (i + 1 < m) continue;
                const uint64_t j = i + 1 - m;
                ASSERT_EQ(st.class_at(j), c.class_of(seq.data() + j)) << "m=" << m << " s=" << s << " j=" << j;
                if (j >= lag_max)
                    ASSERT_EQ(st.class_at(j - lag_max), c.class_of(seq.data() + j - lag_max))
                        << "lagged read, m=" << m << " s=" << s << " j=" << j;
            }
            ASSERT_EQ(st.pushed(), seq.size());
        }
    }
}

// ---- 4. ranks -----------------------------------------------------------------------------------
TEST(MinimizerScheme, RandomRankIsPhiAndOcRankAddsClass) {
    using kuint = uint64_t;
    std::mt19937_64 rng(505);
    for (const auto& [k, m] : std::vector<std::pair<uint64_t, uint64_t>>{{31, 15}, {21, 11}, {15, 7}, {12, 6}}) {
        km::SkmerManipulator<kuint> rnd{k, m};
        km::SkmerManipulator<kuint> oc{k, m, 0, km::SchemeParams::open_closed(m)};
        ASSERT_FALSE(rnd.scheme().is_open_closed());
        ASSERT_TRUE(oc.scheme().is_open_closed());
        const km::OcClassifier& cls = *oc.oc_classifier();
        for (int r = 0; r < 500; ++r) {
            const std::vector<uint8_t> x8 = random_nucls(m, rng);
            const std::vector<kuint> S(x8.begin(), x8.end());
            const kuint r_rnd = rnd.mmer_canonical_rank(S, 0);
            const kuint r_oc  = oc.mmer_canonical_rank(S, 0);
            // random: exactly φ of the canonical interleaved code (the historical rank)
            km::Skmer<kuint> f = rnd.build_skmer_from_nucleotides(S, 0, m - 1, 0);
            rnd.canonicalize(f);
            ASSERT_EQ(r_rnd, rnd.phi(rnd.minimizer(f)));
            ASSERT_EQ(r_rnd, rnd.minimizer_rank(f));
            // open-closed: same φ in the low 2m bits, the class above
            ASSERT_EQ(r_oc & ((kuint{1} << (2 * m)) - 1), r_rnd);
            ASSERT_EQ(r_oc >> (2 * m), static_cast<kuint>(cls.class_of(x8.data())));
            ASSERT_EQ(r_oc, oc.minimizer_rank(f, cls.class_of(x8.data())));
            ASSERT_EQ(oc.minimizer_rank(f, 0), r_rnd);   // class 0 => bit-identical to random
        }
    }
}

// Equal rank <=> same canonical m-mer (the premise of the ambiguity/tie logic), by full enumeration.
TEST(MinimizerScheme, OcRankInjectiveOnCanonicalMmers) {
    using kuint = uint64_t;
    for (const uint64_t m : {5ull, 6ull, 7ull}) {
        const uint64_t k = 2 * m + 1;
        km::SkmerManipulator<kuint> oc{k, m, 0, km::SchemeParams::open_closed(m)};
        std::map<kuint, std::string> rank_to_canon;
        std::vector<uint64_t> class_count(3, 0);
        for (uint64_t code = 0; code < (uint64_t{1} << (2 * m)); ++code) {
            std::vector<kuint> S(m);
            std::string str;
            for (uint64_t j = 0; j < m; ++j) {
                S[j] = (code >> (2 * j)) & 3u;
                str.push_back("ACTG"[S[j]]);
            }
            const std::string canon = sktest::canonical(str);
            const kuint r = oc.mmer_canonical_rank(S, 0);
            ++class_count[r >> (2 * m)];
            auto [it, inserted] = rank_to_canon.emplace(r, canon);
            if (!inserted) ASSERT_EQ(it->second, canon) << "two different canonical m-mers share a rank, m=" << m;
        }
        std::set<std::string> canons;
        for (const auto& [r, c] : rank_to_canon) canons.insert(c);
        EXPECT_EQ(canons.size(), rank_to_canon.size()) << "one canonical m-mer with two ranks, m=" << m;
        // Class 2 needs an s-mer that is neither central nor an end: m - s + 1 >= 4 (not m=5, s=3).
        const uint64_t W = m - km::SchemeParams::default_s(m) + 1;
        for (int c = 0; c < (W >= 4 ? 3 : 2); ++c)
            EXPECT_GT(class_count[c], 0u) << "class " << c << " never occurs at m=" << m;
        if (W < 4) EXPECT_EQ(class_count[2], 0u);
    }
}

// ---- 5. SchemeParams ----------------------------------------------------------------------------
TEST(MinimizerScheme, DefaultSFollowsTheRule) {
    using SP = km::SchemeParams;
    for (uint64_t m = 5; m <= 17; m += 2) EXPECT_EQ(SP::default_s(m), 3u) << "m=" << m;
    for (uint64_t m = 19; m <= 63; m += 2) EXPECT_EQ(SP::default_s(m), 5u) << "m=" << m;
    for (uint64_t m = 8; m <= 64; m += 2) EXPECT_EQ(SP::default_s(m), 4u) << "m=" << m;
    EXPECT_EQ(SP::default_s(4), 2u);
    EXPECT_EQ(SP::default_s(6), 2u);
    EXPECT_EQ(SP::default_s(3), 1u);   // clamped: 3 > m-2
    EXPECT_EQ(SP::default_s(2), 0u);   // no valid s
    EXPECT_EQ(SP::default_s(1), 0u);
    for (uint64_t m = 3; m <= 128; ++m) {
        const SP p = SP::open_closed(m);
        EXPECT_TRUE(p.validation_error(m).empty()) << "default s invalid for m=" << m << ": " << p.validation_error(m);
        EXPECT_EQ((m - p.s) % 2, 0u);
    }
}

TEST(MinimizerScheme, ValidateRangeAndParity) {
    using SP = km::SchemeParams;
    EXPECT_NO_THROW(SP{}.validate(1));                          // random: always valid
    EXPECT_NO_THROW(SP::open_closed(15, 3).validate(15));
    EXPECT_NO_THROW(SP::open_closed(15, 13).validate(15));
    EXPECT_NO_THROW(SP::open_closed(16, 2).validate(16));
    EXPECT_NO_THROW(SP::open_closed(3, 1).validate(3));
    EXPECT_THROW(SP::open_closed(15, 4).validate(15), std::invalid_argument);   // parity
    EXPECT_THROW(SP::open_closed(16, 5).validate(16), std::invalid_argument);   // parity
    EXPECT_THROW(SP::open_closed(15, 15).validate(15), std::invalid_argument);  // s > m-2
    EXPECT_THROW(SP::open_closed(2, 2).validate(2), std::invalid_argument);     // m < 3
    EXPECT_THROW(SP::open_closed(71, 33).validate(71), std::invalid_argument);  // s > 32
    const std::string parity_msg = SP::open_closed(15, 4).validation_error(15);
    EXPECT_NE(parity_msg.find("parity"), std::string::npos) << parity_msg;
    EXPECT_NE(parity_msg.find("s=3"), std::string::npos) << "should suggest the default: " << parity_msg;
    // A manipulator refuses invalid parameters too (library callers bypass the CLI).
    EXPECT_THROW((km::SkmerManipulator<uint64_t>{31, 15, 0, SP::open_closed(15, 4)}), std::invalid_argument);
}

TEST(MinimizerScheme, NamesParseAndHeaderEncoding) {
    using SP = km::SchemeParams;
    EXPECT_EQ(SP::parse("random"), km::MinimizerScheme::random);
    EXPECT_EQ(SP::parse("open-closed"), km::MinimizerScheme::open_closed);
    EXPECT_FALSE(SP::parse("oc").has_value());
    EXPECT_EQ(SP{}.to_string(), "random");
    EXPECT_EQ(SP::open_closed(15).to_string(), "open-closed s=3");
    EXPECT_EQ(SP{}.class_bits(), 0u);
    EXPECT_EQ(SP::open_closed(15).class_bits(), 2u);
    EXPECT_EQ(SP{}.header_id(), 0u);
    EXPECT_EQ(SP{}.header_param(), 0u);
    EXPECT_EQ(SP::open_closed(31).header_id(), 1u);
    EXPECT_EQ(SP::open_closed(31).header_param(), 5u);
    EXPECT_EQ(SP::from_header(0, 0, 15, "f"), SP{});
    EXPECT_EQ(SP::from_header(1, 3, 15, "f"), SP::open_closed(15, 3));
    EXPECT_THROW(SP::from_header(2, 0, 15, "f"), std::runtime_error);   // unknown scheme
    EXPECT_THROW(SP::from_header(0, 3, 15, "f"), std::runtime_error);   // random with a parameter
    EXPECT_THROW(SP::from_header(1, 4, 15, "f"), std::runtime_error);   // wrong parity
    EXPECT_THROW(SP::from_header(1, 0, 15, "f"), std::runtime_error);   // s = 0
    EXPECT_NE(SP::open_closed(15, 3), SP::open_closed(15, 5));
    EXPECT_NE(SP::open_closed(15, 3), SP{});
}

// The open-closed rank needs 2m+2 bits in one generation word: the width bumps only where 2m sits
// within 2 bits of a word, and never for random (whose width must stay exactly as before).
TEST(MinimizerScheme, GenerationWidthHoldsTheRank) {
    using km::sortedlist::select_generation_width_bytes;
    using km::sortedlist::select_width_bytes;
    for (uint64_t k = 4; k <= 120; ++k)
        for (uint64_t m = 3; m < k; ++m) {
            if (2 * (2 * k - m) > 512 || 4 * m + 4 > 512) continue;
            const uint64_t r = select_generation_width_bytes(k, m);
            EXPECT_EQ(r, std::max(select_width_bytes(2 * (2 * k - m)), select_width_bytes(4 * m)))
                << "random width changed at k=" << k << " m=" << m;
            const uint64_t o = select_generation_width_bytes(k, m, km::SchemeParams::open_closed(m));
            EXPECT_GE(8 * o, 2 * m + 2) << "k=" << k << " m=" << m;
            EXPECT_GE(o, r);
            if (2 * m + 2 <= 8 * r) EXPECT_EQ(o, r) << "needless widening at k=" << k << " m=" << m;
        }
    EXPECT_EQ(select_generation_width_bytes(40, 32), 8u);
    EXPECT_EQ(select_generation_width_bytes(40, 32, km::SchemeParams::open_closed(32)), 16u);
    EXPECT_THROW(select_generation_width_bytes(128, 128, km::SchemeParams::open_closed(128)), std::runtime_error);
    EXPECT_NO_THROW(select_generation_width_bytes(127, 127, km::SchemeParams::open_closed(127)));
}
