// Regression tests for the phantom k-mer emitted when a sequence starts with the m-mer of maximal
// minimizer rank.
//
// SeqSkmerator started every sequence with m_current_minimizer = ~0 as a "no minimizer yet"
// sentinel, and took the first candidate as the new minimizer only if its rank was strictly lower.
// But phi is a bijection on the 2m-bit minimizer space: when 2m fills the generation word (m = 64 at
// __uint128_t, k = 65..96; m = 128 at kuint256, k = 129..192), one canonical m-mer ranks exactly
// ~0. A sequence starting with it (either strand) took the equal-minimizer branch against the empty
// slot before the sequence, which handed the first candidate a full left flank over positions ahead
// of the sequence: one k-mer that is not in the input (k-m stale bases, then the m-mer) was emitted.
// Construction indexed it (set cardinalities off by one, the phantom in every materialized result
// and in the ASCII dump), and a query printed one flag more than the sequence has k-mers. At
// m = 16 (uint32) and m = 32 (uint64) the all-ones preimage is not the canonical orientation, so
// it is unreachable there. A has-minimizer flag now makes the first candidate the minimizer.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SetOperations.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace {

const char NUC[4] = {'A', 'C', 'T', 'G'};

// The canonical m-mers whose phi-rank is all-ones (computed with an exact phi inverse; each test
// re-checks the rank, so a change of phi shows up as a failed precondition, not a silent pass).
const std::string kMaxRank64 {"GTCGTTCCAAAACATATGAGCTACGCAACGCTATATGTCGCTCTCCTTCTTAGCCATTTATATC"};
const std::string kMaxRank128 {"GAGAAACTGGCCGGTATTTAAACGTAAGAGACTAGATGATATTTTGGAACGAGTTCACAAGCGACAGGTGCGAC"
                               "GATCGAACCGTCCCAAAGTGCGATAGTTACTTTTCGTGAGTCCGGAAAATCTGA"};

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = (c == 'A') ? 'T' : (c == 'T') ? 'A' : (c == 'C') ? 'G' : 'C';
    return r;
}

std::string canon(const std::string& s) { return std::min(s, revcomp(s)); }

std::string random_seq(std::mt19937_64& rng, uint64_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = NUC[rng() & 3u];
    return s;
}

template <typename kuint>
bool has_max_rank(uint64_t k, uint64_t m, const std::string& mmer) {
    const km::SkmerManipulator<kuint> manip {k, m};
    std::vector<kuint> S;
    for (char c : mmer) S.push_back(static_cast<kuint>((c >> 1) & 0b11));
    return manip.mmer_canonical_rank(S, 0) == static_cast<kuint>(~static_cast<kuint>(0));
}

// Number of k-mers carried by the yielded records of `seq` (no decoding needed).
template <typename kuint>
uint64_t emitted_kmers(km::SkmerManipulator<kuint>& manip, std::string seq) {
    km::SeqSkmerator<kuint> rator {manip, seq};
    uint64_t n {0};
    for (auto it = rator.begin(); !it.consumed(); ++it) {
        const km::Skmer<kuint> sk {*it};
        n += sk.m_pref_size + sk.m_suff_size - (manip.k - manip.m) + 1;
    }
    return n;
}

// Sequences starting (and, for contrast, ending) with the max-rank m-mer on either strand, on a
// fresh manipulator and after another sequence: every k-mer must be emitted exactly once.
template <typename kuint>
void expect_no_extra_kmer(uint64_t k, uint64_t m, const std::string& mmer) {
    ASSERT_EQ(km::sortedlist::select_generation_width_bytes(k, m), sizeof(kuint)) << "k=" << k << " m=" << m;
    ASSERT_TRUE(has_max_rank<kuint>(k, m, mmer)) << "precondition: phi changed, recompute the max-rank m-mer";
    std::mt19937_64 rng(k);
    for (int rep {0}; rep < 8; ++rep) {
        std::string seq {(rep & 1) ? revcomp(mmer) : mmer};
        seq += random_seq(rng, (k - m) + rng() % (2 * k));
        if (rep & 2) seq = revcomp(seq);   // the max-rank m-mer at the END: never affected
        km::SkmerManipulator<kuint> manip {k, m};
        if (rep & 4) (void)emitted_kmers(manip, random_seq(rng, 2 * k));
        EXPECT_EQ(emitted_kmers(manip, seq), seq.size() - k + 1) << "k=" << k << " m=" << m << " seq=" << seq;
    }
}

} // namespace

TEST(FirstMinimizerSentinel, NoExtraKmerAtM64) {
    for (const uint64_t k : {uint64_t{65}, uint64_t{80}, uint64_t{96}})
        expect_no_extra_kmer<__uint128_t>(k, 64, kMaxRank64);
}

TEST(FirstMinimizerSentinel, NoExtraKmerAtM128) {
    for (const uint64_t k : {uint64_t{129}, uint64_t{150}, uint64_t{192}})
        expect_no_extra_kmer<km::kuint256>(k, 128, kMaxRank128);
}

// The phantom k-mer, decoded: on the unfixed code the first record holds the k-m bases ahead of
// the sequence (zeros = 'A') followed by the max-rank m-mer.
TEST(FirstMinimizerSentinel, EveryEmittedKmerIsInTheInput) {
    const uint64_t k {80}, m {64};
    std::mt19937_64 rng(80);
    const std::string seq {kMaxRank64 + random_seq(rng, 92)};
    std::set<std::string> input;
    for (uint64_t i {0}; i + k <= seq.size(); ++i) input.insert(canon(seq.substr(i, k)));
    km::SkmerManipulator<__uint128_t> manip {k, m};
    std::string s {seq};
    km::SeqSkmerator<__uint128_t> rator {manip, s};
    for (auto it = rator.begin(); !it.consumed(); ++it) {
        km::Skmer<__uint128_t> sk {*it};
        manip.unpermute_minimizer_slot(sk);
        std::string str;
        for (__uint128_t v : manip.decode_to_nucleotides(sk)) str += NUC[static_cast<unsigned>(v)];
        for (uint64_t t {0}; t + k <= str.size(); ++t)
            EXPECT_TRUE(input.count(canon(str.substr(t, k)))) << "k-mer not in the input: " << str.substr(t, k);
    }
}

// End to end: the index holds exactly the input's k-mers and a self-query prints one '1' per k-mer.
TEST(FirstMinimizerSentinel, IndexAndQueryHoldOnlyInputKmers) {
    std::mt19937_64 rng(81);
    const std::string dir {::testing::TempDir()};
    for (const auto& [k, m, mmer] : {std::tuple<uint64_t, uint64_t, std::string>{80, 64, kMaxRank64},
                                     {150, 128, kMaxRank128}}) {
        const std::vector<std::string> recs {random_seq(rng, 2 * k), mmer + random_seq(rng, k + 7),
                                             revcomp(mmer) + random_seq(rng, k / 2)};
        const std::string fa {dir + "bug11.fa"}, idx {dir + "bug11.sskm"};
        { std::ofstream f(fa); for (size_t i {0}; i < recs.size(); ++i) f << ">r" << i << "\n" << recs[i] << "\n"; }
        std::set<std::string> input;
        uint64_t nkmers {0};
        for (const std::string& r : recs)
            for (uint64_t i {0}; i + k <= r.size(); ++i, ++nkmers) input.insert(canon(r.substr(i, k)));

        km::sortedlist::SortedListBuildParams p;
        p.k = k;
        p.m = m;
        p.input_path = fa;
        p.output_path = idx;
        p.buckets = 64;
        p.has_output_file = true;
        const uint64_t qb {km::sortedlist::quotient_bits_for(p)};
        const uint64_t gen_w {km::sortedlist::select_generation_width_bytes(k, m)};
        const uint64_t store_w {km::sortedlist::select_width_bytes(2 * (2 * k - m) - qb)};
        km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
            km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
                for (const unsigned threads : {1u, 4u}) {
                    p.n_threads = threads;
                    km::sortedlist::build_sorted_list<gen, store>(p, qb);
                    auto A {km::sortedlist::BucketedSkmerListReader<store>::open(idx)};
                    auto B {km::sortedlist::BucketedSkmerListReader<store>::open(idx)};
                    EXPECT_EQ(km::sortedlist::set_sizes<store>(A, B, threads).inter, input.size())
                        << "k=" << k << " m=" << m << " threads=" << threads;
                    std::stringstream out;
                    km::sortedlist::sequential_query<gen, store>(A, fa, out);
                    const std::string flags {out.str()};
                    EXPECT_EQ(static_cast<uint64_t>(std::count(flags.begin(), flags.end(), '1')), nkmers) << "k=" << k;
                    EXPECT_EQ(flags.find('0'), std::string::npos) << "k=" << k;
                }
            });
        });
    }
}
