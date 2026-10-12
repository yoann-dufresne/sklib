// Regression tests for m == k (window w = k-m+1 = 1), which used to crash `sskm construct`.
//
// The CLI accepted 1 <= m <= k, but SeqSkmerator's sliding window cannot express a window of one
// m-mer: its out-of-context test (m_ptr_current - m_ptr_min >= k - m) fired on every base, so
// recompute_minimizer ran on an empty window, read the stale rank of the slot about to be
// overwritten, and on an equal rank update_on_equal_mini compared a candidate with itself and set
// its left flank to -1, i.e. 65535 in the uint16_t size, which indexed the per-size mask tables out
// of bounds: SIGSEGV for every input (`sskm construct -k 5 -m 5`). Its "slot holds no k-mer" test
// (pref + suff < k - m) also cannot exist when k - m = 0. SeqSkmerator now yields each k-mer
// directly, as its own super-k-mer with empty flanks; the rest of the pipeline (one column) works
// unchanged. k = m = 1 stays rejected: the interleaved layout needs one prefix and one suffix slot
// (2k-m >= 2), so SkmerManipulator now throws std::invalid_argument for k < 2 or m outside [1, k].
//
// Pinned here: the producer emits each k-mer once, as itself; an index built at m == k answers
// exactly (sequential and parallel query, every store width, quotiented down to 0 retained bits);
// the set operations' cardinalities match the k-mer set algebra; invalid parameters are rejected.

#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
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

std::set<std::string> kmer_set(uint64_t k, const std::vector<std::string>& seqs) {
    std::set<std::string> out;
    for (const std::string& s : seqs)
        for (uint64_t i {0}; i + k <= s.size(); ++i) out.insert(canon(s.substr(i, k)));
    return out;
}

// k == m for every CLI generation width: 4 bytes (k <= 16), 8 (k <= 32), 16 (k <= 64), 32 (above),
// up to k = m = 128, the widest minimizer (2m = 256 bits).
const uint64_t kKs[] = {2, 3, 4, 5, 8, 11, 16, 17, 21, 32, 33, 47, 64, 65, 96, 97, 127, 128};

} // namespace

// The producer yields every k-mer exactly once, in input order, as itself (forward or reverse
// complement), with empty flanks and its last base as creation index.
TEST(MinimizerEqualsK, ProducerYieldsEachKmerAsItsOwnSuperKmer) {
    std::mt19937_64 rng(10);
    for (const uint64_t k : kKs) {
        const uint64_t m {k};
        std::vector<std::string> seqs {std::string(k, 'A'), std::string(3 * k + 1, 'C'), random_seq(rng, 4 * k + 3)};
        std::string tandem;
        while (tandem.size() < 3 * k) tandem += "ACGT";
        seqs.push_back(tandem);
        seqs.push_back(std::string(k - 1, 'G'));   // shorter than k: nothing
        km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m), [&]<typename kuint>() {
            km::SkmerManipulator<kuint> manip {k, m};
            for (std::string seq : seqs) {
                km::SeqSkmerator<kuint, true> rator {manip, seq};
                uint64_t i {0};
                for (auto it = rator.begin(); !it.consumed(); ++it, ++i) {
                    km::Skmer<kuint> sk {*it};
                    ASSERT_LT(i + k, seq.size() + 1) << "more records than k-mers, k=" << k << " seq=" << seq;
                    EXPECT_EQ(sk.m_pref_size, 0u);
                    EXPECT_EQ(sk.m_suff_size, 0u);
                    EXPECT_EQ(it.yielded_position(), static_cast<int64_t>(i + k - 1));
                    manip.unpermute_minimizer_slot(sk);
                    std::string got;
                    for (kuint v : manip.decode_to_nucleotides(sk)) got += NUC[static_cast<unsigned>(v)];
                    EXPECT_EQ(canon(got), canon(seq.substr(i, k))) << "k=" << k << " position " << i;
                }
                EXPECT_EQ(i, seq.size() >= k ? seq.size() - k + 1 : 0) << "k=" << k << " seq=" << seq;
            }
        });
    }
}

// End to end at m == k: build, query (sequential and parallel), set operations, against brute force.
TEST(MinimizerEqualsK, IndexQueryAndSetOperationsAreExact) {
    std::mt19937_64 rng(11);
    for (const uint64_t k : {uint64_t{2}, uint64_t{3}, uint64_t{5}, uint64_t{8}, uint64_t{16}, uint64_t{21},
                             uint64_t{31}, uint64_t{33}, uint64_t{64}, uint64_t{65}, uint64_t{96},
                             uint64_t{128}}) {
        const uint64_t m {k};
        std::vector<std::string> a {random_seq(rng, 300), std::string(2 * k, 'A'), random_seq(rng, 40)};
        std::vector<std::string> b {a[0].substr(50, 120), revcomp(a[1]), random_seq(rng, 200)};
        std::string pal;
        while (pal.size() < 2 * k + 3) pal += "AT";
        a.push_back(pal);
        const std::string dir {::testing::TempDir()};
        const std::string fa {dir + "bug10_a.fa"}, fb {dir + "bug10_b.fa"};
        { std::ofstream f(fa); for (size_t i {0}; i < a.size(); ++i) f << ">a" << i << "\n" << a[i] << "\n"; }
        { std::ofstream f(fb); for (size_t i {0}; i < b.size(); ++i) f << ">b" << i << "\n" << b[i] << "\n"; }
        const std::set<std::string> ka {kmer_set(k, a)}, kb {kmer_set(k, b)};

        // 1 bucket (no quotient), 16, and 4096 buckets: at small k the latter quotients ALL 2k bits
        // away (each bucket holds at most one k-mer, stored in 0 retained bits).
        for (const uint64_t buckets : {uint64_t{1}, uint64_t{16}, uint64_t{4096}}) {
            const std::string ia {dir + "bug10_a.sskm"}, ib {dir + "bug10_b.sskm"};
            km::sortedlist::SortedListBuildParams p;
            p.k = k;
            p.m = m;
            p.buckets = buckets;
            p.has_output_file = true;
            const uint64_t qb {km::sortedlist::quotient_bits_for(p)};
            const uint64_t gen_w {km::sortedlist::select_generation_width_bytes(k, m)};
            const uint64_t store_w {km::sortedlist::select_width_bytes(2 * (2 * k - m) - qb)};
            const std::string ctx {"k=m=" + std::to_string(k) + " buckets=" + std::to_string(buckets)};
            km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
                km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
                    for (const unsigned threads : {1u, 4u}) {
                        p.n_threads = threads;
                        p.input_path = fa; p.output_path = ia; km::sortedlist::build_sorted_list<gen, store>(p, qb);
                        p.input_path = fb; p.output_path = ib; km::sortedlist::build_sorted_list<gen, store>(p, qb);
                        auto A {km::sortedlist::BucketedSkmerListReader<store>::open(ia)};
                        auto B {km::sortedlist::BucketedSkmerListReader<store>::open(ib)};

                        // Query B's records against A: one flag per k-mer, '1' iff the k-mer is in A.
                        std::stringstream seq_out, par_out;
                        km::sortedlist::sequential_query<gen, store>(A, fb, seq_out);
                        km::sortedlist::parallel_query<gen, store>(A, fb, par_out, 4);
                        EXPECT_EQ(seq_out.str(), par_out.str()) << ctx;
                        std::string expected;
                        for (const std::string& s : b)
                            for (uint64_t i {0}; i + k <= s.size(); ++i)
                                expected += ka.count(canon(s.substr(i, k))) ? "1\n" : "0\n";
                        EXPECT_EQ(seq_out.str(), expected) << ctx << " threads=" << threads;

                        const km::sortedlist::SetSizes sz {km::sortedlist::set_sizes<store>(A, B, threads)};
                        uint64_t inter {0};
                        for (const std::string& x : ka) inter += kb.count(x);
                        EXPECT_EQ(sz.inter, inter) << ctx;
                        EXPECT_EQ(sz.only_a, ka.size() - inter) << ctx;
                        EXPECT_EQ(sz.only_b, kb.size() - inter) << ctx;

                        // A materialized union answers exactly too.
                        const std::string iu {dir + "bug10_union.sskm"};
                        EXPECT_EQ(km::sortedlist::set_union<store>(A, B, iu, false, threads), ka.size() + kb.size() - inter) << ctx;
                        auto U {km::sortedlist::BucketedSkmerListReader<store>::open(iu)};
                        std::stringstream u_out;
                        km::sortedlist::sequential_query<gen, store>(U, fa, u_out);
                        EXPECT_EQ(u_out.str().find('0'), std::string::npos) << ctx;
                    }
                });
            });
        }

        // In-RAM (ASCII) path: every line is one k-mer with empty flanks.
        km::sortedlist::SortedListBuildParams p;
        p.k = k;
        p.m = m;
        p.input_path = fa;
        p.output_path = dir + "bug10_a.txt";
        p.ascii = true;
        km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m), [&]<typename gen>() {
            km::sortedlist::build_sorted_list<gen, gen>(p, 0);
        });
        std::ifstream in(p.output_path);
        std::string line;
        std::getline(in, line);
        std::set<std::string> dumped;
        std::string nucl, pref, suff;
        while (in >> nucl >> pref >> suff) {
            EXPECT_EQ(pref + suff, "00");
            dumped.insert(canon(nucl));
        }
        EXPECT_EQ(dumped, ka) << "k=m=" << k;
    }
}

TEST(MinimizerEqualsK, InvalidParametersAreRejected) {
    EXPECT_THROW((km::SkmerManipulator<uint32_t>{1, 1}), std::invalid_argument);   // k = m = 1: 2k-m = 1 slot
    EXPECT_THROW((km::SkmerManipulator<uint32_t>{5, 0}), std::invalid_argument);
    EXPECT_THROW((km::SkmerManipulator<uint32_t>{5, 6}), std::invalid_argument);
    EXPECT_THROW((km::SkmerManipulator<uint64_t>{0, 0}), std::invalid_argument);
    EXPECT_NO_THROW((km::SkmerManipulator<uint32_t>{2, 2}));
    EXPECT_NO_THROW((km::SkmerManipulator<uint32_t>{2, 1}));
}
