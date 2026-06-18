// Tests for the bucketed (VSKMER_3) sorted skmer list and its lazily-loaded reader.
//
// Core invariant: routing a query to its minimizer-prefix bucket and searching only that
// sub-list must return exactly what searching the whole globally-sorted list returns, for
// any bucket count. We verify that against the single-bucket list (loaded fully into RAM)
// as ground truth, plus the all-present property (every k-mer of the source sequence is in
// a list built from it), absence, and VSKMER_2 backward-read compatibility.

#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>

namespace {

using kuint = uint64_t;

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

// φ-permuted canonical super-k-mers of `seq`, exactly as the query path enumerates them.
std::vector<km::Skmer<kuint>> enumerate_skmers(uint64_t k, uint64_t m, std::string seq) {
    km::SkmerManipulator<kuint> manip{k, m};
    km::SeqSkmerator<kuint> rator{manip, seq};
    std::vector<km::Skmer<kuint>> out;
    for (const km::Skmer<kuint>& s : rator) out.push_back(s);
    return out;
}

// Build a bucketed list file from `seq` with `buckets` buckets (the real construction path),
// returning its path.
std::string build_list(uint64_t k, uint64_t m, const std::string& seq,
                       uint64_t buckets, const std::string& tag) {
    const std::string in_path  = ::testing::TempDir() + "bsl_in_"  + tag + ".fa";
    const std::string out_path = ::testing::TempDir() + "bsl_out_" + tag + ".sskm";
    write_fasta(in_path, seq);

    km::sortedlist::SortedListBuildParams params;
    params.k = k;
    params.m = m;
    params.input_path = in_path;
    params.output_path = out_path;
    params.ascii = false;
    params.buckets = buckets;
    params.has_output_file = true; // gates the bucketed (seekable) path
    km::sortedlist::build_sorted_list<kuint>(params);
    return out_path;
}

} // namespace

// The lazy bucketed reader must agree with the single global list for every bucket count,
// and every k-mer of the source sequence must be reported present.
TEST(BucketedSkmerList, ReaderMatchesGlobalAcrossBucketCounts) {
    constexpr uint64_t k{15}, m{5};
    const std::string seq = random_dna(800, 12345);
    const auto queries = enumerate_skmers(k, m, seq);
    ASSERT_FALSE(queries.empty());

    // Ground truth: a single-bucket list loaded entirely into RAM.
    auto truth = km::sortedlist::VirtualSkmerSerializer<kuint>::load(build_list(k, m, seq, 1, "truth"));

    for (uint64_t buckets : {uint64_t{1}, uint64_t{2}, uint64_t{8}, uint64_t{64}}) {
        const std::string path = build_list(k, m, seq, buckets, "b" + std::to_string(buckets));
        auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(path);

        for (const auto& q : queries) {
            const std::vector<uint8_t> expected = truth.query_skmer(q);
            const std::vector<uint8_t> got = reader.query_skmer(q);
            ASSERT_EQ(expected, got) << "bucket count " << buckets << " disagrees with the global list";
            for (uint8_t flag : got)
                ASSERT_EQ(flag, 1) << "a k-mer of the source sequence is reported absent (buckets=" << buckets << ")";
        }
    }
}

// Unrelated k-mers must route (possibly into empty buckets) and report not-found — and still
// match the global list exactly.
TEST(BucketedSkmerList, AbsentKmersMatchGlobalAndReportNotFound) {
    constexpr uint64_t k{15}, m{5};
    const std::string present = random_dna(800, 1);
    const std::string absent  = random_dna(400, 2); // disjoint 15-mers w.h.p.

    auto truth = km::sortedlist::VirtualSkmerSerializer<kuint>::load(build_list(k, m, present, 1, "abs_truth"));
    auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(build_list(k, m, present, 64, "abs"));

    uint64_t zeros = 0;
    for (const auto& q : enumerate_skmers(k, m, absent)) {
        const std::vector<uint8_t> expected = truth.query_skmer(q);
        const std::vector<uint8_t> got = reader.query_skmer(q);
        ASSERT_EQ(expected, got);
        for (uint8_t flag : got) zeros += (flag == 0);
    }
    ASSERT_GT(zeros, 0u) << "absent sequence unexpectedly fully present (test not exercising absence)";
}

// A legacy VSKMER_2 file (no bucket directory) must be read as a single bucket spanning the
// whole list, by both load() and the bucketed reader.
// Legacy VSKMER_2 files carry a φ-only minimizer slot; the ψ-permuted reader (VSKMER_5) MUST
// reject them (rebuild required) rather than silently mis-route/mis-decode. (Pre-V5 this read the
// V2 file as a single bucket.)
TEST(BucketedSkmerList, RejectsLegacyV2) {
    constexpr uint64_t k{15}, m{5};

    // Hand-write the legacy format: magic(VSKMER_2) + k + m + count(0) + no payload.
    const std::string v2_path = ::testing::TempDir() + "legacy_v2.bin";
    {
        std::ofstream out(v2_path, std::ios::binary);
        const uint64_t magic = km::sortedlist::util::ENDIANNESS_SANITY_INTEGER; // "VSKMER_2"
        const uint64_t kk = k, mm = m, count = 0;
        out.write(reinterpret_cast<const char*>(&magic), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&kk), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&mm), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
    }

    EXPECT_THROW(km::sortedlist::VirtualSkmerSerializer<kuint>::load(v2_path), std::runtime_error);
    EXPECT_THROW(km::sortedlist::BucketedSkmerListReader<kuint>::open(v2_path), std::runtime_error);
}

// parallel_query() must produce byte-identical output to the sequential reader.query(), in input
// order, for every thread count. A fresh reader per run also exercises concurrent bucket loading.
TEST(BucketedSkmerList, ParallelQueryMatchesSequential) {
    constexpr uint64_t k{15}, m{5};
    const std::string present = random_dna(2000, 11);
    const std::string absent  = random_dna(1500, 22); // mostly-absent 15-mers w.h.p.

    const std::string list_path = build_list(k, m, present, 64, "par");

    // Query file mixes present (1s) and absent (mostly 0s) sequences across two records.
    const std::string qfa = ::testing::TempDir() + "par_query.fa";
    {
        std::ofstream f(qfa);
        f << ">present\n" << present << "\n>absent\n" << absent << "\n";
    }

    std::stringstream expected;
    {
        auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(list_path);
        reader.query(qfa, expected);
    }
    ASSERT_FALSE(expected.str().empty()) << "ground-truth query produced no output";

    for (unsigned threads : {1u, 2u, 4u, 8u}) {
        auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(list_path);
        std::stringstream got;
        km::sortedlist::parallel_query<kuint>(reader, qfa, got, threads);
        ASSERT_EQ(expected.str(), got.str()) << "parallel output diverges at threads=" << threads;
    }
}

// Inputs that yield no super-k-mers (records shorter than k) must produce empty output without
// hanging or crashing the pipeline.
TEST(BucketedSkmerList, ParallelQueryHandlesEmptyInput) {
    constexpr uint64_t k{15}, m{5};
    const std::string list_path = build_list(k, m, random_dna(400, 3), 8, "par_empty_src");

    const std::string qfa = ::testing::TempDir() + "par_empty_query.fa";
    {
        std::ofstream f(qfa);
        f << ">tooshort\n" << random_dna(k - 1, 9) << "\n"; // shorter than k -> skipped
    }

    for (unsigned threads : {2u, 8u}) {
        auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(list_path);
        std::stringstream got;
        km::sortedlist::parallel_query<kuint>(reader, qfa, got, threads);
        ASSERT_TRUE(got.str().empty()) << "expected empty output at threads=" << threads;
    }
}

// A tiny batch size forces many batches, so the sink must reassemble many out-of-order batch
// completions; the in-order output must still match the sequential path exactly.
TEST(BucketedSkmerList, ParallelQueryReordersAcrossManyFlushes) {
    constexpr uint64_t k{15}, m{5};
    const std::string present = random_dna(4000, 31);
    const std::string list_path = build_list(k, m, present, 64, "par_flush");

    const std::string qfa = ::testing::TempDir() + "par_flush_query.fa";
    write_fasta(qfa, present);

    std::stringstream expected;
    {
        auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(list_path);
        reader.query(qfa, expected);
    }

    auto reader = km::sortedlist::BucketedSkmerListReader<kuint>::open(list_path);
    std::stringstream got;
    km::sortedlist::parallel_query<kuint>(reader, qfa, got, 8, /*batch_size=*/64);
    ASSERT_EQ(expected.str(), got.str());
}
