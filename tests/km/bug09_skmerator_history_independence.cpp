// Regression tests for the history-dependent super-k-mer stream of SeqSkmerator.
//
// A SkmerManipulator outlives a sequence: FileSkmerator reuses it record after record, and every
// ParallelConstruct / ParallelQuery worker reuses its own across work items. SeqSkmerator's
// init_seq() did not clear it, so the window slots before a new sequence's start still held the
// previous sequence's tail. Those slots are masked at yield, but add_nucleotide picks each
// candidate's orientation on the unmasked pair: for an RC-palindromic minimizer (even m only) the
// stale flank decided the orientation, hence update_on_equal_mini's cuts, hence where the
// super-k-mers of the new sequence start and end. Frames stay canonical, so membership, the index
// and set operations were unaffected; but `sskm query` prints one line per super-k-mer, so the line
// grouping of a record depended on the record its worker had processed before. That broke the
// "byte-identical output regardless of thread count" invariant (query -t 1 vs -t N), and the
// sequential output of a record also depended on its predecessor in the file.
//
// Minimal case: k=4, m=2, "TATAT". Alone it yields TATA | ATAT; after "GGGG" on the same
// manipulator it yielded TATAT as one super-k-mer.
//
// The contract pinned here: the stream of a sequence is a function of the sequence alone.

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace {

const char NUC[4] = {'A', 'C', 'T', 'G'};

template <typename kuint>
using Raw = std::tuple<kuint, kuint, uint16_t, uint16_t>;

template <typename kuint>
Raw<kuint> raw(const km::Skmer<kuint>& s) {
    return {s.m_pair.m_value[1], s.m_pair.m_value[0], s.m_pref_size, s.m_suff_size};
}

// The yielded records of `seq`, on the manipulator `manip` (whatever it processed before).
template <typename kuint>
std::vector<Raw<kuint>> stream(km::SkmerManipulator<kuint>& manip, const std::string& seq) {
    std::string s {seq};
    km::SeqSkmerator<kuint> rator {manip, s};
    std::vector<Raw<kuint>> out;
    for (auto it = rator.begin(); !it.consumed(); ++it) out.push_back(raw(*it));
    return out;
}

// Same, but through a reused iterator re-primed with reset() (the FileSkmerator / worker pattern).
template <typename kuint>
std::vector<std::vector<Raw<kuint>>> streams_with_reset(uint64_t k, uint64_t m, const std::vector<std::string>& seqs) {
    km::SkmerManipulator<kuint> manip {k, m};
    std::string buf;
    km::SeqSkmerator<kuint> rator {manip, buf};
    auto it {rator.begin()};
    std::vector<std::vector<Raw<kuint>>> out;
    for (const std::string& s : seqs) {
        buf = s;
        it.reset();
        out.emplace_back();
        for (; !it.consumed(); ++it) out.back().push_back(raw(*it));
    }
    return out;
}

template <typename kuint>
std::vector<Raw<kuint>> fresh_stream(uint64_t k, uint64_t m, const std::string& seq) {
    km::SkmerManipulator<kuint> manip {k, m};
    return stream(manip, seq);
}

// (k, m) grid: every generation width, even m (RC-palindromic minimizers exist) and odd m, small m
// and m close to k.
struct KM { uint64_t k, m; };
const KM kGrid[] = {
    {4, 2}, {6, 2}, {9, 4}, {11, 5}, {12, 6}, {16, 8}, {21, 10},         // 4-byte generation width
    {21, 2}, {31, 14}, {31, 15}, {36, 2}, {40, 20},                     // 8 bytes
    {63, 8}, {63, 30}, {50, 41},                                        // 16 bytes
    {95, 40}, {127, 62}, {100, 94},                                     // 32 bytes (kuint256)
};

std::string tile(const std::string& unit, size_t len) {
    std::string s;
    while (s.size() < len) s += unit;
    return s.substr(0, len);
}

} // namespace

// The minimal case, on a fresh SeqSkmerator per sequence (the positional query / library path) and
// through FileSkmerator (the sequential `sskm query -i` / in-RAM construct path).
TEST(SkmeratorHistoryIndependence, MinimalReproducer) {
    using kuint = uint32_t;
    const uint64_t k {4}, m {2};
    const auto alone {fresh_stream<kuint>(k, m, "TATAT")};
    ASSERT_EQ(alone.size(), 2u);   // TATA | ATAT

    km::SkmerManipulator<kuint> manip {k, m};
    (void)stream(manip, "GGGG");
    EXPECT_EQ(stream(manip, "TATAT"), alone) << "the stream of TATAT depends on the sequence processed before it";

    const std::string fa {::testing::TempDir() + "bug09_minimal.fa"};
    { std::ofstream f(fa); f << ">a\nGGGG\n>b\nTATAT\n"; }
    km::SkmerManipulator<kuint> fmanip {k, m};
    km::FileSkmerator<kuint> file {fmanip, fa};
    std::vector<Raw<kuint>> all;
    for (const km::Skmer<kuint>& s : file) all.push_back(raw(s));
    std::vector<Raw<kuint>> expected {fresh_stream<kuint>(k, m, "GGGG")};
    expected.insert(expected.end(), alone.begin(), alone.end());
    EXPECT_EQ(all, expected) << "FileSkmerator: a record's stream depends on the previous record";
}

// Randomized: repeat-rich sequences (palindromic units, homopolymers, random flanks) streamed after
// a random history, through both a fresh SeqSkmerator on a reused manipulator and a reused iterator
// re-primed with reset(). Every stream must equal the fresh-manipulator stream.
TEST(SkmeratorHistoryIndependence, RandomizedHistory) {
    std::mt19937_64 rng(20261012);
    auto urand = [&](uint64_t lo, uint64_t hi) { return std::uniform_int_distribution<uint64_t>(lo, hi)(rng); };
    auto random_seq = [&](uint64_t n) {
        std::string s(n, 'A');
        for (char& c : s) c = NUC[urand(0, 3)];
        return s;
    };
    static const char* kUnits[] = {"AT", "TA", "GC", "CG", "ACGT", "AATT", "GGCC", "TGCA", "A", "C"};
    for (const KM& g : kGrid) {
        const uint64_t k {g.k}, m {g.m};
        const int reps {k > 40 ? 3 : 12};
        km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m), [&]<typename kuint>() {
            std::vector<std::string> seqs;
            for (int r {0}; r < reps; ++r) {
                seqs.push_back(random_seq(urand(k, 3 * k)));                                  // history
                std::string s {random_seq(urand(0, 3))};
                s += tile(kUnits[urand(0, 9)], urand(k, 4 * k));
                if (urand(0, 1)) s += random_seq(urand(1, k));
                seqs.push_back(s);
            }
            km::SkmerManipulator<kuint> shared {k, m};
            const auto reset_streams {streams_with_reset<kuint>(k, m, seqs)};
            for (size_t i {0}; i < seqs.size(); ++i) {
                const auto ref {fresh_stream<kuint>(k, m, seqs[i])};
                EXPECT_EQ(stream(shared, seqs[i]), ref) << "k=" << k << " m=" << m << " seq=" << seqs[i]
                                                        << (i ? " after " + seqs[i - 1] : std::string{});
                EXPECT_EQ(reset_streams[i], ref) << "reset() path, k=" << k << " m=" << m << " seq=" << seqs[i];
            }
        });
    }
}

// End to end: `sskm query` must print the same bytes whatever the thread count and the work-item
// size, and a record's lines must not depend on the records before it in the file.
TEST(SkmeratorHistoryIndependence, QueryOutputIndependentOfThreadsAndNeighbours) {
    std::mt19937_64 rng(7);
    auto urand = [&](uint64_t lo, uint64_t hi) { return std::uniform_int_distribution<uint64_t>(lo, hi)(rng); };
    std::vector<std::string> records;
    for (int i {0}; i < 400; ++i) {
        std::string r(urand(4, 40), 'A');
        for (char& c : r) c = NUC[urand(0, 3)];
        records.push_back(r);
        records.push_back(tile(urand(0, 1) ? "TA" : "GC", urand(4, 60)) + std::string(urand(0, 3), 'C'));
    }
    const std::string fa {::testing::TempDir() + "bug09_query.fa"};
    {
        std::ofstream f(fa);
        for (size_t i {0}; i < records.size(); ++i) f << ">r" << i << "\n" << records[i] << "\n";
    }
    for (const KM& g : {KM{4, 2}, KM{9, 4}, KM{21, 2}, KM{36, 2}, KM{63, 8}}) {
        const uint64_t k {g.k}, m {g.m};
        const std::string idx {::testing::TempDir() + "bug09_query_k" + std::to_string(k) + ".sskm"};
        km::sortedlist::SortedListBuildParams p;
        p.k = k;
        p.m = m;
        p.input_path = fa;
        p.output_path = idx;
        p.buckets = 16;
        p.has_output_file = true;
        const uint64_t b {km::sortedlist::quotient_bits_for(p)};
        const uint64_t gen_w {km::sortedlist::select_generation_width_bytes(k, m)};
        const uint64_t store_w {km::sortedlist::select_width_bytes(2 * (2 * k - m) - b)};
        km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
            km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
                km::sortedlist::build_sorted_list<gen, store>(p, b);
                auto reader {km::sortedlist::BucketedSkmerListReader<store>::open(idx)};
                std::stringstream seq_out;
                km::sortedlist::sequential_query<gen, store>(reader, fa, seq_out);

                // Each record queried on its own: concatenated, it must equal the whole-file output.
                std::string one_by_one;
                const std::string single {::testing::TempDir() + "bug09_single.fa"};
                for (const std::string& r : records) {
                    { std::ofstream f(single); f << ">r\n" << r << "\n"; }
                    std::stringstream o;
                    km::sortedlist::sequential_query<gen, store>(reader, single, o);
                    one_by_one += o.str();
                }
                EXPECT_EQ(seq_out.str(), one_by_one) << "a record's lines depend on its predecessors, k=" << k << " m=" << m;

                // Small work items spread the records over the workers in a run-dependent way.
                for (const char* chunk : {"1", "37", "300"}) {
                    setenv("SKLIB_QUERY_CHUNK_BP", chunk, 1);
                    for (unsigned threads : {2u, 8u}) {
                        std::stringstream par_out;
                        km::sortedlist::parallel_query<gen, store>(reader, fa, par_out, threads);
                        EXPECT_EQ(par_out.str(), seq_out.str()) << "-t " << threads << " with " << chunk
                                                                << " bp work items, k=" << k << " m=" << m;
                    }
                    unsetenv("SKLIB_QUERY_CHUNK_BP");
                }
            });
        });
    }
}
