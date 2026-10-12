// End-to-end tests of the minimizer scheme option (construct --minimizer open-closed), driven
// through the same builders / query drivers / set operations the CLI uses.
//
// 10. Differential: the same FASTA indexed with random and with open-closed encodes the same k-mer
//     set, so the input's k-mers, their reverse complements, k-mers of a mutated copy and random
//     k-mers get IDENTICAL flags (and both equal the ground truth) — at every record width.
// 11. -t 1 vs -t 8 byte-identity under open-closed: construct (with an intra-sequence chunk seam),
//     query (with forced chunk seams) and set operations.
// 12. setop: open-closed results equal the random results as k-mer sets; mixed schemes are refused;
//     results carry the scheme.
// 13. Format: VSKMER_6 files (built by the pre-V7 code, tests/data/v6_payload_ref_*.sskm) still read
//     as random; a random V7 file's directory and payload are byte-identical to V6's; V7 open-closed
//     round-trips (bucketed, in-RAM and --max-ram paths); legacy / foreign / corrupted headers are
//     refused.

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include <io/MinimizerScheme.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/SetOperations.hpp>
#include <algorithms/WidthDispatch.hpp>

#include "scheme_test_utils.hpp"
#include "test_data.hpp"

namespace {

using sktest::canonical;
using sktest::random_dna;
using sktest::revcomp;

std::string tmp(const std::string& name) { return ::testing::TempDir() + "oce2e_" + name; }

// Input with random sequence plus repeat runs (homopolymers, STRs) that reach the ambiguous path.
std::vector<std::string> mixed_input(uint64_t k, uint64_t seed, uint64_t len) {
    std::string a = random_dna(len, seed);
    a += sktest::tile("A", 3 * k) + random_dna(200, seed + 1) + sktest::tile("ACG", 4 * k) +
         random_dna(200, seed + 2) + sktest::tile("AC", 3 * k) + random_dna(150, seed + 3);
    return {a, random_dna(len / 4, seed + 4), random_dna(k - 1, seed + 5) /* < k: skipped */};
}

std::set<std::string> canonical_kmers(const std::vector<std::string>& seqs, uint64_t k) {
    std::set<std::string> s;
    for (const std::string& x : seqs)
        for (size_t i = 0; i + k <= x.size(); ++i) s.insert(canonical(x.substr(i, k)));
    return s;
}

// A copy with one substitution every `every` bases (always a different base).
std::string mutate(std::string s, uint64_t every, uint64_t seed) {
    static const std::string nt = "ACGT";
    std::mt19937_64 rng(seed);
    for (size_t i = every / 2; i < s.size(); i += every) {
        const size_t cur = nt.find(s[i]);
        s[i] = nt[(cur + 1 + rng() % 3) % 4];
    }
    return s;
}

// One FASTA record per k-mer, so line i of the query output is the flag of k-mer i.
std::vector<std::string> single_kmer_queries(const std::vector<std::string>& input, uint64_t k,
                                             uint64_t n_random, uint64_t seed) {
    std::vector<std::string> q;
    for (const std::string& x : input)
        for (size_t i = 0; i + k <= x.size(); ++i) {
            q.push_back(x.substr(i, k));
            q.push_back(revcomp(x.substr(i, k)));
        }
    const std::string mut = mutate(input[0], 37, seed);
    for (size_t i = 0; i + k <= mut.size(); i += 3) q.push_back(mut.substr(i, k));
    std::mt19937_64 rng(seed);
    for (uint64_t r = 0; r < n_random; ++r) q.push_back(random_dna(k, rng()));
    return q;
}

struct Cfg { uint64_t k, m; const char* why; };
// Store widths at --buckets 256 (b = 8): 2*(2k-m) - 8 bits -> uint32 / uint64 / __uint128 / kuint256.
// (40,32): open-closed needs a 66-bit rank, so its GENERATION width is __uint128 vs uint64 for random.
const Cfg kWidths[] = {{15, 9, "uint32 store"}, {31, 15, "uint64 store"}, {45, 21, "u128 store"},
                       {90, 41, "kuint256 store"}, {40, 32, "open-closed bumps the generation width"}};

} // namespace

// ---- 10. differential ----------------------------------------------------------------------------
TEST(SchemeEndToEnd, RandomAndOpenClosedAnswerIdentically) {
    for (const Cfg& c : kWidths) {
        const std::string tag = std::to_string(c.k) + "_" + std::to_string(c.m);
        const std::vector<std::string> input = mixed_input(c.k, 100 + c.k, 12000);
        const std::string in_fa = tmp("diff_in_" + tag + ".fa");
        sktest::write_fasta(in_fa, input);
        const std::vector<std::string> q = single_kmer_queries(input, c.k, c.k == 31 ? 100000 : 20000, 7 + c.k);
        const std::string q_fa = tmp("diff_q_" + tag + ".fa");
        sktest::write_fasta(q_fa, q);

        const std::set<std::string> truth_set = canonical_kmers(input, c.k);
        std::vector<std::string> truth;
        uint64_t zeros = 0;
        for (const std::string& x : q) {
            truth.push_back(truth_set.count(canonical(x)) ? "1" : "0");
            zeros += truth.back() == "0";
        }
        ASSERT_GT(zeros, 0u);

        std::vector<std::vector<std::string>> answers;
        for (const km::SchemeParams& sp : {sktest::random_scheme(), sktest::oc_scheme(c.m)}) {
            const std::string idx = tmp("diff_" + tag + (sp.is_open_closed() ? "_oc" : "_rnd") + ".sskm");
            sktest::build_index(c.k, c.m, sp, in_fa, idx, /*buckets=*/256, /*threads=*/1);
            const km::sortedlist::ListHeaderInfo h = km::sortedlist::read_list_header(idx);
            EXPECT_EQ(h.scheme, sp);
            EXPECT_EQ(h.version, 7u);
            answers.push_back(sktest::lines_of(sktest::query_file(idx, q_fa, 1)));
            const std::string ctx = "k=" + std::to_string(c.k) + " m=" + std::to_string(c.m) + " " +
                                    sp.to_string() + " (" + c.why + ")";
            ASSERT_EQ(answers.back().size(), q.size()) << ctx;
            uint64_t wrong = 0;
            for (size_t i = 0; i < q.size(); ++i) wrong += answers.back()[i] != truth[i];
            EXPECT_EQ(wrong, 0u) << ctx << ": answers differ from the ground truth";
        }
        EXPECT_EQ(answers[0], answers[1]) << "random and open-closed disagree at k=" << c.k << " m=" << c.m;
    }
}

// ---- 11. thread-count byte identity ---------------------------------------------------------------
TEST(SchemeEndToEnd, OpenClosedConstructByteIdenticalAcrossThreads) {
    // One sequence longer than the 1 Mbp phase-1 chunk, so the parallel build splits it (seam), plus
    // short records.
    const std::vector<std::string> input = {random_dna((uint64_t{1} << 20) + 150000, 31),
                                            random_dna(5000, 32), sktest::tile("ACGTT", 900)};
    const std::string in_fa = tmp("thr_in.fa");
    sktest::write_fasta(in_fa, input);
    for (const Cfg& c : {Cfg{31, 15, ""}, Cfg{45, 21, ""}}) {
        const km::SchemeParams sp = sktest::oc_scheme(c.m);
        const std::string tag = std::to_string(c.k);
        const std::string s1 = tmp("thr_t1_" + tag + ".sskm"), s8 = tmp("thr_t8_" + tag + ".sskm");
        sktest::build_index(c.k, c.m, sp, in_fa, s1, 4096, 1);
        sktest::build_index(c.k, c.m, sp, in_fa, s8, 4096, 8);
        const std::string a = sktest::read_file_bytes(s1), b = sktest::read_file_bytes(s8);
        ASSERT_FALSE(a.empty());
        EXPECT_TRUE(a == b) << "open-closed index differs between -t 1 and -t 8 at k=" << c.k;
    }
}

TEST(SchemeEndToEnd, OpenClosedQueryByteIdenticalAcrossThreads) {
    const uint64_t k = 31, m = 15;
    const std::vector<std::string> input = mixed_input(k, 41, 60000);
    const std::string in_fa = tmp("qthr_in.fa"), idx = tmp("qthr.sskm"), q_fa = tmp("qthr_q.fa");
    sktest::write_fasta(in_fa, input);
    sktest::build_index(k, m, sktest::oc_scheme(m), in_fa, idx, 4096, 1);
    sktest::write_fasta(q_fa, {input[0], revcomp(input[0]), mutate(input[0], 53, 3), random_dna(20000, 44), input[1]});
    // Force intra-record chunking (seams every 3 kbp) in the parallel driver.
    ::setenv("SKLIB_QUERY_CHUNK_BP", "3000", 1);
    const std::string seq_out = sktest::query_file(idx, q_fa, 1);
    const std::string par_out = sktest::query_file(idx, q_fa, 8);
    ::unsetenv("SKLIB_QUERY_CHUNK_BP");
    ASSERT_FALSE(seq_out.empty());
    EXPECT_TRUE(seq_out == par_out) << "open-closed query output differs between -t 1 and -t 8";
}

// ---- 12. set operations -----------------------------------------------------------------------------
namespace {

struct SetopRun {
    km::sortedlist::SetSizes sizes;
    std::vector<std::string> out_paths;   // inter, union, diff_ab, diff_ba, xor
};

SetopRun run_multi(const std::string& a, const std::string& b, const std::string& tag, unsigned threads) {
    const km::sortedlist::ListHeaderInfo h = km::sortedlist::read_list_header(a);
    SetopRun r;
    for (const char* rel : {"inter", "union", "dab", "dba", "xor"}) r.out_paths.push_back(tmp("so_" + tag + "_" + rel + ".sskm"));
    km::sortedlist::dispatch_width_bytes(h.store_width_bytes, [&]<typename store>() {
        auto A = km::sortedlist::BucketedSkmerListReader<store>::open(a);
        auto B = km::sortedlist::BucketedSkmerListReader<store>::open(b);
        km::sortedlist::MultiSetOpRequest req;
        req.inter_out = r.out_paths[0];
        req.union_out = r.out_paths[1];
        req.diff_ab_out = r.out_paths[2];
        req.diff_ba_out = r.out_paths[3];
        req.sym_diff_out = r.out_paths[4];
        r.sizes = km::sortedlist::multi_setop<store>(A, B, req, threads).sizes;
    });
    return r;
}

} // namespace

TEST(SchemeEndToEnd, OpenClosedSetOpsMatchRandomAsKmerSets) {
    const uint64_t k = 31, m = 15;
    const std::string sa = random_dna(30000, 51) + sktest::tile("AAC", 200);
    const std::string sb = sa.substr(15000) + random_dna(15000, 52);
    const std::string fa_a = tmp("so_a.fa"), fa_b = tmp("so_b.fa");
    sktest::write_fasta(fa_a, {sa});
    sktest::write_fasta(fa_b, {sb});
    const std::set<std::string> A = canonical_kmers({sa}, k), B = canonical_kmers({sb}, k);
    uint64_t inter = 0;
    for (const std::string& x : A) inter += B.count(x);

    // Probe every k-mer of A and B (plus some absent ones) against each result list.
    std::vector<std::string> probes;
    for (const std::string* s : {&sa, &sb})
        for (size_t i = 0; i + k <= s->size(); ++i) probes.push_back(s->substr(i, k));
    for (int r = 0; r < 2000; ++r) probes.push_back(random_dna(k, 900 + r));
    const std::string probe_fa = tmp("so_probe.fa");
    sktest::write_fasta(probe_fa, probes);

    std::vector<std::vector<std::string>> flags_per_scheme;
    for (const km::SchemeParams& sp : {sktest::random_scheme(), sktest::oc_scheme(m)}) {
        const std::string tag = sp.is_open_closed() ? "oc" : "rnd";
        const std::string ia = tmp("so_a_" + tag + ".sskm"), ib = tmp("so_b_" + tag + ".sskm");
        sktest::build_index(k, m, sp, fa_a, ia, 256, 1);
        sktest::build_index(k, m, sp, fa_b, ib, 256, 1);
        const SetopRun r1 = run_multi(ia, ib, tag + "_t1", 1);
        const SetopRun r8 = run_multi(ia, ib, tag + "_t8", 8);
        EXPECT_EQ(r1.sizes.inter, inter) << tag;
        EXPECT_EQ(r1.sizes.only_a, A.size() - inter) << tag;
        EXPECT_EQ(r1.sizes.only_b, B.size() - inter) << tag;
        EXPECT_EQ(r8.sizes.inter, r1.sizes.inter);
        EXPECT_EQ(r8.sizes.only_a, r1.sizes.only_a);
        EXPECT_EQ(r8.sizes.only_b, r1.sizes.only_b);
        std::vector<std::string> flags;
        for (size_t i = 0; i < r1.out_paths.size(); ++i) {
            // (11) every set-op output is byte-identical across thread counts
            EXPECT_TRUE(sktest::read_file_bytes(r1.out_paths[i]) == sktest::read_file_bytes(r8.out_paths[i]))
                << tag << " output " << i << " differs between -t 1 and -t 8";
            const km::sortedlist::ListHeaderInfo h = km::sortedlist::read_list_header(r1.out_paths[i]);
            EXPECT_EQ(h.scheme, sp) << "set-op result must carry the inputs' scheme";
            const std::string out = sktest::query_file(r1.out_paths[i], probe_fa, 1);
            flags.push_back(out);
        }
        flags_per_scheme.push_back(flags);
    }
    for (size_t i = 0; i < flags_per_scheme[0].size(); ++i)
        EXPECT_TRUE(flags_per_scheme[0][i] == flags_per_scheme[1][i])
            << "set-op result " << i << ": open-closed and random encode different k-mer sets";
}

TEST(SchemeEndToEnd, SetOpRefusesMixedSchemes) {
    const uint64_t k = 31, m = 15;
    const std::string fa = tmp("mix.fa");
    sktest::write_fasta(fa, {random_dna(8000, 61)});
    const std::string ir = tmp("mix_rnd.sskm"), io = tmp("mix_oc.sskm"), io5 = tmp("mix_oc5.sskm");
    sktest::build_index(k, m, sktest::random_scheme(), fa, ir, 256, 1);
    sktest::build_index(k, m, sktest::oc_scheme(m, 3), fa, io, 256, 1);
    sktest::build_index(k, m, sktest::oc_scheme(m, 5), fa, io5, 256, 1);
    auto sizes = [&](const std::string& a, const std::string& b) {
        auto A = km::sortedlist::BucketedSkmerListReader<uint64_t>::open(a);
        auto B = km::sortedlist::BucketedSkmerListReader<uint64_t>::open(b);
        return km::sortedlist::set_sizes<uint64_t>(A, B, 1);
    };
    EXPECT_NO_THROW(sizes(io, io));
    EXPECT_THROW(sizes(ir, io), std::runtime_error);
    EXPECT_THROW(sizes(io, io5), std::runtime_error);   // same scheme, different s
    try {
        sizes(ir, io);
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("open-closed s=3"), std::string::npos) << e.what();
    }
}

// ---- 13. format --------------------------------------------------------------------------------------
namespace {

uint64_t u64_at(const std::string& bytes, size_t off) {
    uint64_t v;
    std::memcpy(&v, bytes.data() + off, sizeof v);
    return v;
}

void put_u64(std::string& bytes, size_t off, uint64_t v) { std::memcpy(bytes.data() + off, &v, sizeof v); }

void write_bytes(const std::string& path, const std::string& bytes) {
    std::ofstream o(path, std::ios::binary | std::ios::trunc);
    o.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

struct Fixture { uint64_t k, m, store; };
const Fixture kV6[] = {{15, 9, 4}, {31, 15, 8}, {45, 21, 16}, {90, 41, 32}};

} // namespace

// VSKMER_6 indexes written by the pre-V7 code are read as the random scheme and answer correctly.
TEST(SchemeFormat, V6StillReadableAsRandom) {
    const std::string fa = test_data("v6_payload_ref.fa");
    for (const Fixture& f : kV6) {
        const std::string v6 = test_data("v6_payload_ref_k" + std::to_string(f.k) + "m" + std::to_string(f.m) + ".sskm");
        const km::sortedlist::ListHeaderInfo h = km::sortedlist::read_list_header(v6);
        EXPECT_EQ(h.version, 6u);
        EXPECT_EQ(h.scheme, km::SchemeParams{});
        EXPECT_EQ(h.store_width_bytes, f.store);
        km::sortedlist::dispatch_width_bytes(h.store_width_bytes, [&]<typename store>() {
            auto r = km::sortedlist::BucketedSkmerListReader<store>::open(v6);
            EXPECT_EQ(r.format_version(), 6u);
            EXPECT_FALSE(r.scheme().is_open_closed());
        });
        const std::string out = sktest::query_file(v6, fa, 1);
        ASSERT_FALSE(out.empty());
        EXPECT_EQ(out.find('0'), std::string::npos) << "a k-mer of the source is missing from the V6 list k=" << f.k;
    }
}

// The random scheme's records are unchanged: a V7 build of the same input is the V6 file with only the
// magic changed and the two scheme words inserted before the directory.
TEST(SchemeFormat, RandomV7PayloadByteIdenticalToV6) {
    const std::string fa = test_data("v6_payload_ref.fa");
    for (const Fixture& f : kV6) {
        const std::string v6 = sktest::read_file_bytes(
            test_data("v6_payload_ref_k" + std::to_string(f.k) + "m" + std::to_string(f.m) + ".sskm"));
        const std::string p7 = tmp("v7_of_v6_" + std::to_string(f.k) + ".sskm");
        sktest::build_index(f.k, f.m, km::SchemeParams{}, fa, p7, /*buckets=*/16, /*threads=*/1);
        const std::string v7 = sktest::read_file_bytes(p7);
        ASSERT_GE(v6.size(), 56u);
        EXPECT_EQ(u64_at(v6, 0), km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V6);
        EXPECT_EQ(u64_at(v7, 0), km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V7);
        EXPECT_EQ(v7.substr(8, 48), v6.substr(8, 48)) << "k/m/count/n_buckets/width/b differ, k=" << f.k;
        EXPECT_EQ(u64_at(v7, 56), 0u);   // scheme_id = random
        EXPECT_EQ(u64_at(v7, 64), 0u);   // scheme_param
        EXPECT_EQ(v7.size(), v6.size() + 16);
        EXPECT_TRUE(v7.substr(72) == v6.substr(56)) << "directory or payload changed for the random scheme, k=" << f.k;
    }
}

// Open-closed round-trips through every writer/reader: bucketed (default), in-RAM (save/load, used for
// stdout / ASCII targets) and adaptive --max-ram.
TEST(SchemeFormat, OpenClosedV7RoundTrip) {
    const uint64_t k = 31, m = 15;
    const km::SchemeParams sp = sktest::oc_scheme(m, 5);
    const std::vector<std::string> input = mixed_input(k, 71, 20000);
    const std::string fa = tmp("rt_in.fa");
    sktest::write_fasta(fa, input);
    std::vector<std::string> probes;
    for (size_t i = 0; i + k <= input[0].size(); i += 7) probes.push_back(revcomp(input[0].substr(i, k)));
    const std::string probe_fa = tmp("rt_probe.fa");
    sktest::write_fasta(probe_fa, probes);

    const std::string bucketed = tmp("rt_bucketed.sskm"), inram = tmp("rt_inram.sskm"), maxram = tmp("rt_maxram.sskm");
    sktest::build_index(k, m, sp, fa, bucketed, 4096, 1);
    sktest::build_index(k, m, sp, fa, inram, 4096, 1, /*has_output_file=*/false);
    sktest::build_index(k, m, sp, fa, maxram, 4096, 1, true, /*max_ram_bytes=*/uint64_t{1} << 20);
    for (const std::string& p : {bucketed, inram, maxram}) {
        const km::sortedlist::ListHeaderInfo h = km::sortedlist::read_list_header(p);
        EXPECT_EQ(h.version, 7u) << p;
        EXPECT_EQ(h.scheme, sp) << p;
        const std::string out = sktest::query_file(p, probe_fa, 1);
        EXPECT_EQ(sktest::lines_of(out).size(), probes.size()) << p;
        EXPECT_EQ(out.find('0'), std::string::npos) << "an indexed k-mer is reported absent in " << p;
    }
    // The in-RAM list file is loadable as one global list; load() restores the scheme so the list's own
    // query(filename) frames queries with it.
    auto list = km::sortedlist::VirtualSkmerSerializer<uint64_t>::load(inram);
    EXPECT_EQ(list.scheme(), sp);
    std::stringstream ss;
    list.query(probe_fa, ss);
    EXPECT_EQ(ss.str().find('0'), std::string::npos);
    // ...and the raw header bytes are as documented: scheme_id = 1, scheme_param = s.
    const std::string bytes = sktest::read_file_bytes(bucketed);
    EXPECT_EQ(u64_at(bytes, 56), 1u);
    EXPECT_EQ(u64_at(bytes, 64), 5u);
}

TEST(SchemeFormat, RejectsLegacyForeignAndCorruptHeaders) {
    const std::string good = tmp("rej_good.sskm");
    const std::string fa = tmp("rej.fa");
    sktest::write_fasta(fa, {random_dna(3000, 81)});
    sktest::build_index(31, 15, sktest::oc_scheme(15), fa, good, 16, 1);
    const std::string bytes = sktest::read_file_bytes(good);
    auto expect_rejected = [&](const std::string& b, const std::string& what, const char* needle) {
        const std::string p = tmp("rej_" + what + ".sskm");
        write_bytes(p, b);
        try {
            km::sortedlist::read_list_header(p);
            ADD_FAILURE() << what << ": read_list_header accepted it";
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find(needle), std::string::npos) << what << ": " << e.what();
        }
        EXPECT_THROW(km::sortedlist::BucketedSkmerListReader<uint64_t>::open(p), std::runtime_error) << what;
        EXPECT_THROW(km::sortedlist::VirtualSkmerSerializer<uint64_t>::load(p), std::runtime_error) << what;
    };
    namespace u = km::sortedlist::util;
    for (const uint64_t legacy : {u::ENDIANNESS_SANITY_INTEGER, u::ENDIANNESS_SANITY_INTEGER_V3,
                                  u::ENDIANNESS_SANITY_INTEGER_V4, u::ENDIANNESS_SANITY_INTEGER_V5}) {
        std::string b = bytes;
        put_u64(b, 0, legacy);
        expect_rejected(b, "legacy" + std::to_string(legacy & 0xFF), "rebuild");
    }
    {
        std::string b = bytes;
        put_u64(b, 0, u::swap_endian(u::ENDIANNESS_SANITY_INTEGER_V7));
        expect_rejected(b, "swapped", "Endianness");
    }
    {
        std::string b = bytes;
        put_u64(b, 0, 0x1234567890ABCDEFULL);
        expect_rejected(b, "garbage", "Invalid file format");
    }
    {
        std::string b = bytes;
        put_u64(b, 56, 9);                    // unknown scheme id
        expect_rejected(b, "unknown_scheme", "Unknown minimizer scheme");
    }
    {
        std::string b = bytes;
        put_u64(b, 64, 4);                    // open-closed with s of the wrong parity for m = 15
        expect_rejected(b, "bad_s", "parity");
    }
}
