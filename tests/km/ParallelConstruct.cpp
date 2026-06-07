// Tests for the parallel per-bucket construction (phase 2 of build_bucketed; see
// ParallelConstruct.hpp). The whole contract of the optimization is: building with N worker
// threads must yield an index that is **byte-identical** to the sequential (single-thread) build,
// for any thread count, any bucketing, and any record width — only the build time changes. These
// tests lock that invariant, plus the degenerate single-bucket / empty-input edge cases that
// stress the producer/throttle/writer loop.

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/WidthDispatch.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>

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

// Build a bucketed index exactly like the CLI run_construct: choose the generation and storage
// integer widths from (k, m) and the quotient bits, then dispatch. `n_threads` selects the
// sequential (1) or parallel (>= 2) phase-2 path.
void build_index(uint64_t k, uint64_t m, const std::string& in_path,
                 const std::string& out_path, uint64_t buckets, unsigned n_threads) {
    km::sortedlist::SortedListBuildParams p;
    p.k = k;
    p.m = m;
    p.input_path = in_path;
    p.output_path = out_path;
    p.ascii = false;
    p.buckets = buckets;
    p.has_output_file = true;
    p.n_threads = n_threads;

    const uint64_t b  = km::sortedlist::quotient_bits_for(p);
    const uint64_t gw = km::sortedlist::select_width_bytes(2 * (2 * k - m));
    const uint64_t sw = km::sortedlist::select_width_bytes(2 * (2 * k - m) - b);
    km::sortedlist::dispatch_width_bytes(gw, [&]<typename gen>() {
        km::sortedlist::dispatch_width_bytes(sw, [&]<typename store>() {
            km::sortedlist::build_sorted_list<gen, store>(p, b);
        });
    });
}

std::string read_file_bytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// First differing byte offset, or npos if equal-prefix up to the shorter length and same length.
size_t first_diff(const std::string& a, const std::string& b) {
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i)
        if (a[i] != b[i]) return i;
    return (a.size() == b.size()) ? std::string::npos : n;
}

} // namespace

// The core invariant: the parallel phase-2 build is byte-identical to the sequential build, for
// every thread count, across record widths (uint64 and __uint128_t) and quotiented / default
// bucketings.
TEST(ParallelConstruct, ByteIdenticalToSequentialAcrossThreadsAndWidths) {
    struct Case { uint64_t k, m, buckets; const char* tag; };
    const std::vector<Case> cases = {
        {21, 11, 256,  "k21m11_q"},     // uint64 store, b>0 (quotiented)
        {31, 15, 4096, "k31m15"},       // uint64 store, default bucket count
        {15, 7,  64,   "k15m7"},        // small m, many records per bucket
        {45, 21, 256,  "k45m21_u128"},  // 2*(2k-m)=138 bits -> __uint128_t store, b>0
    };
    const std::string seq = random_dna(50000, 4242); // enough to populate many buckets

    for (const Case& c : cases) {
        const std::string in_path = ::testing::TempDir() + "pc_in_" + c.tag + ".fa";
        write_fasta(in_path, seq);

        const std::string seq_path = ::testing::TempDir() + "pc_seq_" + c.tag + ".sskm";
        build_index(c.k, c.m, in_path, seq_path, c.buckets, /*n_threads=*/1);
        const std::string truth = read_file_bytes(seq_path);
        ASSERT_FALSE(truth.empty()) << "sequential build produced no file for " << c.tag;

        for (unsigned threads : {2u, 4u, 8u, 16u}) {
            const std::string par_path = ::testing::TempDir() + "pc_par_" + c.tag
                                         + "_t" + std::to_string(threads) + ".sskm";
            build_index(c.k, c.m, in_path, par_path, c.buckets, threads);
            const std::string got = read_file_bytes(par_path);
            ASSERT_EQ(truth.size(), got.size())
                << "index size differs at " << c.tag << " threads=" << threads;
            const size_t d = first_diff(truth, got);
            ASSERT_EQ(d, std::string::npos)
                << "bytes differ at offset " << d << " for " << c.tag << " threads=" << threads;
        }
    }
}

// Degenerate shapes that stress the producer/throttle/writer loop: a single bucket (exactly one
// job) and an input shorter than k (zero super-k-mers => every bucket empty => zero jobs). Neither
// may hang or diverge from the sequential build.
TEST(ParallelConstruct, HandlesSingleBucketAndEmptyInput) {
    {   // buckets = 1: the parallel driver has exactly one job.
        const std::string in_path = ::testing::TempDir() + "pc_one_in.fa";
        write_fasta(in_path, random_dna(5000, 7));
        const std::string s1 = ::testing::TempDir() + "pc_one_seq.sskm";
        const std::string p1 = ::testing::TempDir() + "pc_one_par.sskm";
        build_index(21, 11, in_path, s1, /*buckets=*/1, /*n_threads=*/1);
        build_index(21, 11, in_path, p1, /*buckets=*/1, /*n_threads=*/8);
        ASSERT_EQ(read_file_bytes(s1), read_file_bytes(p1));
    }
    {   // input shorter than k: no super-k-mers, so all buckets are empty (zero jobs).
        const std::string in_path = ::testing::TempDir() + "pc_empty_in.fa";
        write_fasta(in_path, random_dna(20, 9)); // < k
        const std::string s1 = ::testing::TempDir() + "pc_empty_seq.sskm";
        const std::string p1 = ::testing::TempDir() + "pc_empty_par.sskm";
        build_index(31, 15, in_path, s1, /*buckets=*/256, /*n_threads=*/1);
        build_index(31, 15, in_path, p1, /*buckets=*/256, /*n_threads=*/8);
        ASSERT_EQ(read_file_bytes(s1), read_file_bytes(p1));
    }
}
