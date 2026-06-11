#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include <filesystem>
#include <gtest/gtest.h>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <io/SkmerDigest.hpp>
#include <algorithms/WidthDispatch.hpp>

// Bit-exact regression guard for the super-k-mer PRODUCER (SeqSkmerator / FileSkmerator /
// SkmerManipulator). It pins the order-sensitive FNV-1a digest (io/SkmerDigest.hpp) of the whole
// yielded skmer stream over deterministic pseudo-random DNA, for a spread of (k,m) that exercise
// every precompiled integer width and the ambiguous/palindrome framing paths.
//
// WHY a brittle exact-value pin here, when the other Skmerator tests deliberately assert the
// order-/φ-agnostic *coverage invariant* instead: those guard semantic CORRECTNESS and must survive
// intended algorithm changes (φ ordering, per-k-mer splitting). THIS test guards a different
// contract -- that a *performance* rewrite of the producer keeps its output BYTE-IDENTICAL. Its job
// is precisely to break if the stream changes at all. The same digest is what the `sskm-produce`
// binary reports on whole genomes (chr21/celegans), so a green run here plus an unchanged genome
// digest certifies an optimization end-to-end.
//
// If the producer's output semantics are changed ON PURPOSE, regenerate the golden values below:
// each EXPECT also prints `GOLDEN k=.. m=.. digest=0x..`, so run the test, copy the printed values
// into kCases, and commit.

namespace {

// Deterministic pseudo-random DNA (xorshift64; pure integer arithmetic, identical on every
// little-endian target). Fixed seed -> the produced stream, and thus the digest, is reproducible.
std::string random_dna(std::size_t n, uint64_t seed) {
    static const char NUC[4] = {'A', 'C', 'G', 'T'};
    std::string s(n, 'A');
    uint64_t x = seed ? seed : 0x9E3779B97F4A7C15ull;
    for (std::size_t i = 0; i < n; ++i) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        s[i] = NUC[(x >> 11) & 0x3u];
    }
    return s;
}

// Digest of the skmer stream a SeqSkmerator yields over `seq`, at the SAME generation width the
// build (and sskm-produce) picks for these k/m: select_width_bytes(2*(2k-m)).
uint64_t digest_seq(uint64_t k, uint64_t m, const std::string& seq) {
    uint64_t out = km::SKMER_DIGEST_INIT;
    std::string s = seq; // SeqSkmerator binds a std::string&
    km::sortedlist::dispatch_width_bytes(
        km::sortedlist::select_width_bytes(2 * (2 * k - m)),
        [&]<typename kuint>() {
            km::SkmerManipulator<kuint> manip{k, m};
            km::SeqSkmerator<kuint> gen{manip, s};
            uint64_t h = km::SKMER_DIGEST_INIT;
            for (const km::Skmer<kuint>& sk : gen) km::skmer_digest_update(h, sk);
            out = h;
        });
    return out;
}

// Same digest, but read through FileSkmerator over a FASTA holding exactly `seq` as one record.
// For a single record longer than k this must equal digest_seq (the file path only adds record
// framing), which ties the binary's whole-genome digest to the in-memory golden.
uint64_t digest_file(uint64_t k, uint64_t m, const std::string& seq) {
    namespace fs = std::filesystem;
    const fs::path p = fs::temp_directory_path() /
        ("sklib_digest_" + std::to_string(::getpid()) + "_" + std::to_string(k) + "_" +
         std::to_string(m) + ".fa");
    {
        std::ofstream o(p, std::ios::binary | std::ios::trunc);
        o << ">seq\n" << seq << "\n";
    }
    uint64_t out = km::SKMER_DIGEST_INIT;
    km::sortedlist::dispatch_width_bytes(
        km::sortedlist::select_width_bytes(2 * (2 * k - m)),
        [&]<typename kuint>() {
            km::SkmerManipulator<kuint> manip{k, m};
            km::FileSkmerator<kuint> rator{manip, p.string()};
            uint64_t h = km::SKMER_DIGEST_INIT;
            for (const km::Skmer<kuint>& sk : rator) km::skmer_digest_update(h, sk);
            out = h;
        });
    std::error_code ec; fs::remove(p, ec);
    return out;
}

struct Case { uint64_t k, m; std::size_t len; uint64_t seed; uint64_t golden; };

// (k, m, length, seed, golden digest). Widths exercised: 4B (k21m11, k15m4, k11m2),
// 8B (k21m7, k31m15), 16B (k63m31). Small-m cases hit the palindrome/ambiguous framing path often.
// Golden values are bootstrapped from the reference (pre-optimization) producer; see file header.
const Case kCases[] = {
    {21, 11, 300000, 0xC0FFEEu,    0xc15f192ce818c0e6ULL},
    {21,  7, 200000, 0x1234567u,   0xf9ebdcc942af56f5ULL},
    {15,  4, 200000, 0xABCDEFu,    0xf4224abcdf6a75b0ULL},
    {11,  2, 100000, 0x55AA55AAu,  0x5bfbdf0cefa9cb9eULL},
    {31, 15, 150000, 0xDEADBEEFu,  0xd2fec771d94ab549ULL},
    {63, 31,  80000, 0x99887766u,  0x02a22d64bfc14d49ULL},
};

} // namespace

TEST(SkmeratorDigest, GoldenStreamPin) {
    for (const Case& c : kCases) {
        const std::string seq = random_dna(c.len, c.seed);
        const uint64_t d = digest_seq(c.k, c.m, seq);
        std::printf("GOLDEN k=%llu m=%llu digest=0x%016llx\n",
                    (unsigned long long)c.k, (unsigned long long)c.m, (unsigned long long)d);
        EXPECT_EQ(d, c.golden) << "producer stream changed for k=" << c.k << " m=" << c.m
                               << " (if intended, update kCases golden to 0x" << std::hex << d << ")";
    }
}

TEST(SkmeratorDigest, FileMatchesSeq) {
    for (const Case& c : kCases) {
        const std::string seq = random_dna(c.len, c.seed);
        EXPECT_EQ(digest_file(c.k, c.m, seq), digest_seq(c.k, c.m, seq))
            << "FileSkmerator vs SeqSkmerator digest mismatch for k=" << c.k << " m=" << c.m;
    }
}
