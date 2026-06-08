// Tests for minimizer-prefix quotienting (VSKMER_4) and runtime-selected record width.
//
// Core invariant: dropping the top b φ-minimizer bits from each record (they are implied by the
// bucket id) and storing in the narrowest integer that fits must NOT change any query result.
// We verify that the query output of a quotiented list is byte-identical to a full-width list
// built with the same k/m/buckets, across same-width (u32/u64/u256) and width-dropping
// (u128 -> u64, u256 -> u128) configurations, plus the bit-level behaviour of truncate_skmer and
// the width selectors.

#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <io/wide_int.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace {

std::string random_dna(uint64_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    static const char alphabet[] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(n);
    for (uint64_t i = 0; i < n; ++i) s.push_back(alphabet[rng() & 3u]);
    return s;
}

void write_fasta(const std::string& path, const std::string& seq) {
    std::ofstream f(path);
    f << ">seq\n" << seq << "\n";
}

// Build a list with explicit generation/storage widths and quotient bits (the dual-width path the
// CLI drives), returning its path.
template<typename gen, typename store>
std::string build_list(uint64_t k, uint64_t m, const std::string& seq,
                       uint64_t buckets, uint64_t b, const std::string& tag) {
    const std::string in_path  = ::testing::TempDir() + "qw_in_"  + tag + ".fa";
    const std::string out_path = ::testing::TempDir() + "qw_out_" + tag + ".sskm";
    write_fasta(in_path, seq);

    km::sortedlist::SortedListBuildParams params;
    params.k = k;
    params.m = m;
    params.input_path = in_path;
    params.output_path = out_path;
    params.ascii = false;
    params.buckets = buckets;
    params.has_output_file = true;
    km::sortedlist::build_sorted_list<gen, store>(params, b);
    return out_path;
}

// Query every super-k-mer of `qfa` against the list, returning the full presence-flag output as a
// string (the gen-aware sequential driver: parse at gen, route, down-convert to store, search).
template<typename gen, typename store>
std::string query_all(const std::string& list_path, const std::string& qfa) {
    auto reader = km::sortedlist::BucketedSkmerListReader<store>::open(list_path);
    std::stringstream ss;
    km::sortedlist::sequential_query<gen, store>(reader, qfa, ss);
    return ss.str();
}

uint64_t file_size(const std::string& p) {
    std::ifstream f(p, std::ios::binary | std::ios::ate);
    return static_cast<uint64_t>(f.tellg());
}

} // namespace

// ---- truncate_skmer bit behaviour ------------------------------------------------------------

// b == 0, store == gen: a width-preserving copy (pair and sizes unchanged).
TEST(TruncateSkmer, IdentityWhenNoQuotient) {
    const uint64_t k = 15, m = 5; // eff = 2*(2k-m) = 50 bits
    using kuint = uint64_t;
    km::Skmer<kuint>::pair p{0x0123456789ABCDEFull & ((1ull << 50) - 1), 0};
    km::Skmer<kuint> s{p, 3, 4};
    km::Skmer<kuint> t = km::truncate_skmer<kuint, kuint>(k, m, 0, s);
    EXPECT_EQ(t.m_pair.m_value[0], s.m_pair.m_value[0]);
    EXPECT_EQ(t.m_pair.m_value[1], s.m_pair.m_value[1]);
    EXPECT_EQ(t.m_pref_size, 3);
    EXPECT_EQ(t.m_suff_size, 4);
}

// b > 0, same width: keeps exactly the low 2*(2k-m)-b bits, clears the rest, keeps the sizes.
TEST(TruncateSkmer, DropsTopBitsSameWidth) {
    const uint64_t k = 15, m = 5, b = 6; // eff = 50 - 6 = 44 bits retained
    using kuint = uint64_t;
    const uint64_t raw = 0xFFFFFFFFFFFFFFFFull & ((1ull << 50) - 1);
    km::Skmer<kuint> s{km::Skmer<kuint>::pair{raw, 0}, 2, 1};
    km::Skmer<kuint> t = km::truncate_skmer<kuint, kuint>(k, m, b, s);
    EXPECT_EQ(t.m_pair.m_value[0], raw & ((1ull << 44) - 1));
    EXPECT_EQ(t.m_pref_size, 2);
    EXPECT_EQ(t.m_suff_size, 1);
}

// Width drop u128 -> u64: the low (2*(2k-m)-b) bits move into the narrower pair unchanged.
TEST(TruncateSkmer, WidthDrop128To64) {
    const uint64_t k = 40, m = 10, b = 12; // eff = 140 - 12 = 128 bits -> fits a u64 pair exactly
    using gen = __uint128_t;
    using store = uint64_t;
    // A value occupying bits across both 64-bit halves of the meaningful 128-bit range.
    __uint128_t lo = 0xDEADBEEF12345678ull;
    __uint128_t hi = 0x0FEDCBA987654321ull; // 60 bits -> within the retained 128 (lo 64 + hi 64)
    km::Skmer<gen>::pair p{static_cast<gen>(lo | (hi << 64))};
    km::Skmer<gen> s{p, 5, 6};
    km::Skmer<store> t = km::truncate_skmer<gen, store>(k, m, b, s);
    EXPECT_EQ(t.m_pair.m_value[0], static_cast<uint64_t>(lo));
    EXPECT_EQ(t.m_pair.m_value[1], static_cast<uint64_t>(hi));
    EXPECT_EQ(t.m_pref_size, 5);
    EXPECT_EQ(t.m_suff_size, 6);
}

// ---- width selectors --------------------------------------------------------------------------

TEST(WidthDispatch, SelectWidthBoundaries) {
    using km::sortedlist::select_width_bytes;
    EXPECT_EQ(select_width_bytes(1), 4u);
    EXPECT_EQ(select_width_bytes(64), 4u);
    EXPECT_EQ(select_width_bytes(65), 8u);
    EXPECT_EQ(select_width_bytes(128), 8u);
    EXPECT_EQ(select_width_bytes(129), 16u);
    EXPECT_EQ(select_width_bytes(256), 16u);
    EXPECT_EQ(select_width_bytes(257), 32u);   // beyond the __uint128_t pair -> kuint256
    EXPECT_EQ(select_width_bytes(512), 32u);
    EXPECT_THROW(select_width_bytes(513), std::runtime_error);
}

TEST(WidthDispatch, EffectiveBucketBits) {
    using km::sortedlist::effective_bucket_bits;
    EXPECT_EQ(effective_bucket_bits(10, 1), 0u);      // 1 bucket -> no prefix
    EXPECT_EQ(effective_bucket_bits(10, 4096), 12u);  // log2(4096), 12 < 2m=20
    EXPECT_EQ(effective_bucket_bits(3, 4096), 6u);    // clamped to 2m = 6
    EXPECT_EQ(effective_bucket_bits(8, 256), 8u);     // log2(256) == 2m
}

// ---- query equivalence: quotiented == full-width ----------------------------------------------

// Same width (u64): a quotiented list (b>0) returns exactly what the full-width list returns.
TEST(QuotientEquivalence, SameWidthU64) {
    constexpr uint64_t k = 25, m = 11; // 4k-2m = 78 -> u64; b = log2(256)=8 < 2m=22
    const std::string seq = random_dna(3000, 4242);
    const std::string qfa = ::testing::TempDir() + "qw_q_u64.fa";
    write_fasta(qfa, seq + random_dna(1500, 99)); // present + mostly-absent

    const std::string full = build_list<uint64_t, uint64_t>(k, m, seq, 256, 0, "u64_full");
    const std::string quot = build_list<uint64_t, uint64_t>(k, m, seq, 256, 8, "u64_quot");

    EXPECT_EQ((query_all<uint64_t, uint64_t>(full, qfa)),
              (query_all<uint64_t, uint64_t>(quot, qfa)));
}

// Small width (u32) quotient is also transparent.
TEST(QuotientEquivalence, SameWidthU32) {
    constexpr uint64_t k = 15, m = 5; // 4k-2m = 50 -> u32; b = log2(64)=6 < 2m=10
    const std::string seq = random_dna(2000, 7);
    const std::string qfa = ::testing::TempDir() + "qw_q_u32.fa";
    write_fasta(qfa, seq + random_dna(1000, 77));

    const std::string full = build_list<uint32_t, uint32_t>(k, m, seq, 64, 0, "u32_full");
    const std::string quot = build_list<uint32_t, uint32_t>(k, m, seq, 64, 6, "u32_quot");

    EXPECT_EQ((query_all<uint32_t, uint32_t>(full, qfa)),
              (query_all<uint32_t, uint32_t>(quot, qfa)));
}

// Width drop u128 -> u64: generate wide, store narrow. Output matches the full u128 list, and the
// quotiented file is strictly smaller on disk (24-byte vs 40-byte records).
TEST(QuotientEquivalence, WidthDrop128To64) {
    constexpr uint64_t k = 40, m = 10; // 4k-2m = 140 -> u128; b = log2(4096)=12; store 128 -> u64
    const std::string seq = random_dna(6000, 2024);
    const std::string qfa = ::testing::TempDir() + "qw_q_drop.fa";
    write_fasta(qfa, seq + random_dna(2500, 31));

    const std::string full = build_list<__uint128_t, __uint128_t>(k, m, seq, 4096, 0, "drop_full");
    const std::string quot = build_list<__uint128_t, uint64_t>(k, m, seq, 4096, 12, "drop_quot");

    EXPECT_EQ((query_all<__uint128_t, __uint128_t>(full, qfa)),
              (query_all<__uint128_t, uint64_t>(quot, qfa)));

    EXPECT_LT(file_size(quot), file_size(full)) << "narrower records should shrink the file";
}

// The quotiented (width-dropping) list must give identical parallel and sequential query output.
TEST(QuotientEquivalence, ParallelMatchesSequentialQuotiented) {
    constexpr uint64_t k = 40, m = 10;
    const std::string seq = random_dna(5000, 55);
    const std::string qfa = ::testing::TempDir() + "qw_par.fa";
    write_fasta(qfa, seq + random_dna(2000, 66));

    const std::string quot = build_list<__uint128_t, uint64_t>(k, m, seq, 4096, 12, "par_quot");

    const std::string expected = query_all<__uint128_t, uint64_t>(quot, qfa);
    for (unsigned threads : {2u, 4u, 8u}) {
        auto reader = km::sortedlist::BucketedSkmerListReader<uint64_t>::open(quot);
        std::stringstream got;
        km::sortedlist::parallel_query<__uint128_t, uint64_t>(reader, qfa, got, threads);
        ASSERT_EQ(expected, got.str()) << "parallel diverges at threads=" << threads;
    }
}

// Every k-mer of the source sequence is reported present in a quotiented list built from it.
TEST(QuotientEquivalence, AllSourceKmersPresent) {
    constexpr uint64_t k = 40, m = 10;
    const std::string seq = random_dna(4000, 9);
    const std::string qfa = ::testing::TempDir() + "qw_present.fa";
    write_fasta(qfa, seq);

    const std::string quot = build_list<__uint128_t, uint64_t>(k, m, seq, 4096, 12, "present_quot");
    const std::string out = query_all<__uint128_t, uint64_t>(quot, qfa);
    ASSERT_FALSE(out.empty());
    for (char c : out)
        ASSERT_TRUE(c == '1' || c == ',' || c == '\n')
            << "a source k-mer is reported absent ('0') in its own quotiented list";
}

// ---- 256-bit backend: large-k lists exceeding the __uint128_t pair (256-bit) -------------------

// Same width (kuint256): a quotiented large-k list (b>0) returns exactly what the full-width list
// returns. k=75 needs 2*(2k-m)=278 bits, beyond the 256-bit __uint128_t pair -> 32-byte records.
TEST(QuotientEquivalence, SameWidthU256) {
    constexpr uint64_t k = 75, m = 11; // 2*(2k-m)=278 -> kuint256; b=log2(4096)=12 < 2m=22
    const std::string seq = random_dna(8000, 1234);
    const std::string qfa = ::testing::TempDir() + "qw_q_u256.fa";
    write_fasta(qfa, seq + random_dna(3000, 88)); // present + mostly-absent

    const std::string full = build_list<km::kuint256, km::kuint256>(k, m, seq, 4096, 0, "u256_full");
    const std::string quot = build_list<km::kuint256, km::kuint256>(k, m, seq, 4096, 12, "u256_quot");

    EXPECT_EQ((query_all<km::kuint256, km::kuint256>(full, qfa)),
              (query_all<km::kuint256, km::kuint256>(quot, qfa)));
}

// Width drop kuint256 -> u128: generate at 512-bit, store at 256-bit. k=70 needs 258 bits (gen
// kuint256); after dropping b=12 the record fits 246 bits -> a __uint128_t store. Output matches
// the full kuint256 list and the narrower store shrinks the file. This is the case where the
// quotient-mask reference MUST be kuint256: the 258-bit full skmer overflows a __uint128_t pair.
TEST(QuotientEquivalence, WidthDrop256To128) {
    constexpr uint64_t k = 70, m = 11; // 2*(2k-m)=258 -> kuint256; b=12; store 246 -> __uint128_t
    const std::string seq = random_dna(8000, 2025);
    const std::string qfa = ::testing::TempDir() + "qw_q_drop256.fa";
    write_fasta(qfa, seq + random_dna(3000, 41));

    const std::string full = build_list<km::kuint256, km::kuint256>(k, m, seq, 4096, 0, "drop256_full");
    const std::string quot = build_list<km::kuint256, __uint128_t>(k, m, seq, 4096, 12, "drop256_quot");

    EXPECT_EQ((query_all<km::kuint256, km::kuint256>(full, qfa)),
              (query_all<km::kuint256, __uint128_t>(quot, qfa)));

    EXPECT_LT(file_size(quot), file_size(full)) << "narrower records should shrink the file";
}

// Every k-mer of the source sequence is reported present in a quotiented 256-bit list built from it.
TEST(QuotientEquivalence, AllSourceKmersPresentU256) {
    constexpr uint64_t k = 75, m = 11;
    const std::string seq = random_dna(6000, 17);
    const std::string qfa = ::testing::TempDir() + "qw_present_u256.fa";
    write_fasta(qfa, seq);

    const std::string quot = build_list<km::kuint256, km::kuint256>(k, m, seq, 4096, 12, "present_u256");
    const std::string out = query_all<km::kuint256, km::kuint256>(quot, qfa);
    ASSERT_FALSE(out.empty());
    for (char c : out)
        ASSERT_TRUE(c == '1' || c == ',' || c == '\n')
            << "a source k-mer is reported absent ('0') in its own quotiented 256-bit list";
}

// ---- backward compatibility: a hand-written VSKMER_3 file still reads as 64-bit, b=0 ----------

TEST(QuotientWidth, ReadsLegacyV3) {
    constexpr uint64_t k = 15, m = 5;
    const std::string seq = random_dna(800, 123);
    const std::string qfa = ::testing::TempDir() + "qw_v3_q.fa";
    write_fasta(qfa, seq);

    // Build a full-width u64 list, load its records, then re-emit them under the legacy V3 layout
    // (magic, k, m, count, n_buckets, directory, payload) — no store_width / quotient_bits fields.
    const std::string v4 = build_list<uint64_t, uint64_t>(k, m, seq, 1, 0, "v3src");
    auto loaded = km::sortedlist::VirtualSkmerSerializer<uint64_t>::load(v4);
    const auto& recs = loaded.get_list();

    const std::string v3_path = ::testing::TempDir() + "legacy_v3.bin";
    {
        std::ofstream out(v3_path, std::ios::binary);
        const uint64_t magic = km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3;
        const uint64_t kk = k, mm = m, count = recs.size(), n_buckets = 1;
        const km::sortedlist::BucketDirEntry dir{0, count};
        out.write(reinterpret_cast<const char*>(&magic), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&kk), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&mm), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&n_buckets), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&dir), sizeof(km::sortedlist::BucketDirEntry));
        out.write(reinterpret_cast<const char*>(recs.data()),
                  static_cast<std::streamsize>(count * sizeof(km::Skmer<uint64_t>)));
    }

    const km::sortedlist::ListHeaderInfo hdr = km::sortedlist::read_list_header(v3_path);
    EXPECT_EQ(hdr.store_width_bytes, 8u);
    EXPECT_EQ(hdr.quotient_bits, 0u);

    auto reader = km::sortedlist::BucketedSkmerListReader<uint64_t>::open(v3_path);
    EXPECT_EQ(reader.n_buckets(), 1u);
    EXPECT_EQ((query_all<uint64_t, uint64_t>(v3_path, qfa)),
              (query_all<uint64_t, uint64_t>(v4, qfa)));
}
