// Regression tests for the duplicate k-mer emission of SeqSkmerator.
//
// When the minimizer went out of context, recompute_minimizer rescanned the window and walked every
// occurrence of the new minimal minimizer, splitting the k-mers shared by consecutive occurrences
// with update_on_equal_mini. Before each split it reset the PREVIOUS occurrence's left flank to its
// widest value. From the 3rd occurrence on, that undid the cut the previous split had just made, so
// two raw super-k-mers claimed the same k-mers and these k-mers were emitted twice. Trigger: >= 3
// occurrences of one canonical minimizer in the rescanned window (homopolymers, short tandem
// repeats, small m). Membership was unaffected (the duplicates carry the same canonical frame and
// construction deduplicates them), but `sskm query` prints one flag per EMITTED k-mer, so it printed
// more flags than the query has k-mers: 46 flags for the 25 k-mers of a 55 bp poly-A query.
//
// These tests pin the producer contract the bug broke:
//   (1) partition: every k-mer of a sequence is emitted exactly once. Records sharing a creation
//       index are the pieces canonical_pieces cut from ONE raw super-k-mer; raw super-k-mers come by
//       increasing creation index and tile the sequence left to right, the pieces of each one in
//       its canonical orientation (forward or reverse order);
//   (2) framing: each emitted k-mer is stored in the frame choose_kmer_minimizer derives from the
//       k-mer alone;
//   (3) strand invariance: the reverse complement of a sequence yields the same multiset of
//       canonical (k-mer, frame) records;
//   (4) query: an `sskm query` (sequential and parallel drivers) prints one flag per k-mer.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
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

const char NUC[4] = {'A', 'C', 'T', 'G'}; // sklib's 2-bit code: (ascii >> 1) & 0b11

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = (c == 'A') ? 'T' : (c == 'T') ? 'A' : (c == 'C') ? 'G' : 'C';
    return r;
}

std::string tile(const std::string& unit, size_t len) {
    std::string s;
    while (s.size() < len) s += unit;
    return s.substr(0, len);
}

// One yielded record, decoded: its nucleotides in stored orientation, its prefix size (the
// minimizer sits at nucl[pref]) and the creation index of its raw super-k-mer.
struct Emitted {
    std::string nucl;
    uint64_t pref;
    int64_t pos;
};

template <typename kuint>
std::vector<Emitted> emit(uint64_t k, uint64_t m, const std::string& seq) {
    km::SkmerManipulator<kuint> manip{k, m};
    std::string s = seq;
    km::SeqSkmerator<kuint, true> rator{manip, s};
    std::vector<Emitted> out;
    for (auto it = rator.begin(); !it.consumed(); ++it) {
        km::Skmer<kuint> sk = *it;
        manip.unpermute_minimizer_slot(sk);
        Emitted e {"", sk.m_pref_size, it.yielded_position()};
        for (kuint v : manip.decode_to_nucleotides(sk)) e.nucl += NUC[static_cast<unsigned>(v)];
        out.push_back(e);
    }
    return out;
}

uint64_t kmer_count(uint64_t k, const std::string& seq) { return seq.size() >= k ? seq.size() - k + 1 : 0; }

// Contract (1). Returns "" when the records tile the k-mers of `seq` exactly, else a reason.
std::string partition_error(uint64_t k, const std::string& seq, const std::vector<Emitted>& recs) {
    auto matches_at = [&](const Emitted& e, uint64_t q) {
        const uint64_t L {e.nucl.size()};
        return q + L <= seq.size() && (seq.compare(q, L, e.nucl) == 0 || seq.compare(q, L, revcomp(e.nucl)) == 0);
    };
    uint64_t next {0};    // first input k-mer not covered yet
    uint64_t emitted {0};
    int64_t last_pos {std::numeric_limits<int64_t>::min()};
    for (size_t i {0}; i < recs.size();) {
        size_t j {i};
        while (j < recs.size() && recs[j].pos == recs[i].pos) ++j;
        if (recs[i].pos <= last_pos)
            return "creation index " + std::to_string(recs[i].pos) + " does not increase";
        last_pos = recs[i].pos;
        for (size_t t {i}; t < j; ++t) {
            if (recs[t].nucl.size() < k) return "record " + recs[t].nucl + " is shorter than k";
            emitted += recs[t].nucl.size() - k + 1;
        }
        bool tiled {false};
        for (const bool reversed : {false, true}) {
            uint64_t q {next};
            bool ok {true};
            for (size_t c {0}; c < j - i && ok; ++c) {
                const Emitted& e {recs[reversed ? j - 1 - c : i + c]};
                ok = matches_at(e, q);
                q += e.nucl.size() - k + 1;
            }
            if (ok) { next = q; tiled = true; break; }
        }
        if (!tiled)
            return "raw super-k-mer created at " + std::to_string(recs[i].pos) + " (first piece " + recs[i].nucl +
                   ") does not continue the tiling at input k-mer " + std::to_string(next);
        i = j;
    }
    const uint64_t n {kmer_count(k, seq)};
    if (emitted != n)
        return "emitted " + std::to_string(emitted) + " k-mers for " + std::to_string(n) + " input k-mers";
    if (next != n) return "covered " + std::to_string(next) + " of " + std::to_string(n) + " input k-mers";
    return "";
}

template <typename kuint>
using FrameKey = std::tuple<kuint, kuint, uint16_t, uint16_t>;

// Contract (2): the single-k-mer frame of every emitted k-mer, as stored; counts the k-mers whose
// stored frame differs from choose_kmer_minimizer on the k-mer alone. choose_kmer_minimizer(S, t)
// reads only S[t..t+k-1] and the ranks of its own m-mers (a function of each m-mer alone), so the
// per-record rank table gives exactly the isolated k-mer's choice.
template <typename kuint>
std::vector<FrameKey<kuint>> frames(uint64_t k, uint64_t m, const std::vector<Emitted>& recs, uint64_t& mismatches) {
    const km::SkmerManipulator<kuint> manip{k, m};
    std::vector<FrameKey<kuint>> keys;
    for (const Emitted& e : recs) {
        std::vector<kuint> S;
        for (char c : e.nucl) S.push_back(static_cast<kuint>((c >> 1) & 0b11));
        if (S.size() < k) continue;
        std::vector<kuint> rank(S.size() - m + 1);
        for (uint64_t p {0}; p < rank.size(); ++p) rank[p] = manip.mmer_canonical_rank(S, p);
        for (uint64_t t {0}; t + k <= S.size(); ++t) {
            const km::Skmer<kuint> stored {manip.build_skmer_from_nucleotides(S, t, t + k - 1, e.pref)};
            uint64_t occ;
            km::orientation_t orient;
            manip.choose_kmer_minimizer(S, t, rank, occ, orient);
            km::Skmer<kuint> expected {manip.build_skmer_from_nucleotides(S, t, t + k - 1, occ)};
            if (orient == km::reverse_c) expected = manip.reverse_complement(expected);
            if (!(expected.m_pair == stored.m_pair) || expected.m_pref_size != stored.m_pref_size
                || expected.m_suff_size != stored.m_suff_size)
                ++mismatches;
            keys.emplace_back(stored.m_pair.m_value[1], stored.m_pair.m_value[0], stored.m_pref_size, stored.m_suff_size);
        }
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

// Contracts (1)-(3) on `seq` and on its reverse complement, at the generation width the CLI picks.
void expect_producer_contract(uint64_t k, uint64_t m, const std::string& seq) {
    const uint64_t width {km::sortedlist::select_generation_width_bytes(k, m)};
    km::sortedlist::dispatch_width_bytes(width, [&]<typename kuint>() {
        const std::string ctx {"k=" + std::to_string(k) + " m=" + std::to_string(m) + " width=" +
                               std::to_string(width) + " seq=" + seq};
        const std::string rc {revcomp(seq)};
        const std::vector<Emitted> fwd_recs {emit<kuint>(k, m, seq)};
        const std::vector<Emitted> rc_recs {emit<kuint>(k, m, rc)};
        EXPECT_EQ(partition_error(k, seq, fwd_recs), "") << ctx;
        EXPECT_EQ(partition_error(k, rc, rc_recs), "") << "reverse complement, " << ctx;
        uint64_t fwd_bad {0}, rc_bad {0};
        const auto fwd_frames {frames<kuint>(k, m, fwd_recs, fwd_bad)};
        const auto rc_frames {frames<kuint>(k, m, rc_recs, rc_bad)};
        EXPECT_EQ(fwd_bad, 0u) << "k-mers framed differently from choose_kmer_minimizer, " << ctx;
        EXPECT_EQ(rc_bad, 0u) << "k-mers framed differently from choose_kmer_minimizer (reverse complement), " << ctx;
        EXPECT_TRUE(fwd_frames == rc_frames) << "canonical (k-mer, frame) multisets differ between strands, " << ctx;
    });
}

// The (k, m) pairs used below: every CLI generation width (4/8/16/32 bytes), small m, m ~ k/2 and
// m close to k, even and odd k. m stays <= 96 because these tests predate the exact phi inverse
// at 2m > 192 (bug12_phi_inverse_wide.cpp), which decoding a record needs; the producer itself
// never inverts phi.
struct KM { uint64_t k, m; };
const KM kGrid[] = {
    {5, 2}, {9, 3}, {11, 5}, {12, 4}, {15, 7}, {16, 13}, {21, 11},        // 4-byte generation width
    {21, 2}, {31, 15}, {31, 7}, {32, 30}, {40, 20},                      // 8 bytes
    {63, 31}, {64, 5}, {50, 41},                                         // 16 bytes
    {95, 41}, {127, 61}, {100, 95},                                      // 32 bytes (kuint256)
};

} // namespace

// Minimal reproducers: the shortest one (k=4, m=1, 6 bp: CCGC was emitted twice; the trigger needs
// k-m >= 3 so the rescanned window can hold 3 occurrences), 13 bp at k=11/m=5 (3 k-mers, used to
// emit 4) and a 55 bp poly-A run at the default-like k=31/m=15 (25 k-mers, used to emit 46).
TEST(SkmeratorDuplicateEmission, MinimalReproducers) {
    for (const auto& [k, m, seq] : {std::tuple<uint64_t, uint64_t, std::string>{4, 1, "ACCGCA"},
                                    {11, 5, "AACCCCCCCCCCC"},
                                    {31, 15, "ACGTTGCA" + std::string(40, 'A') + "TTGCAGT"}}) {
        km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m), [&]<typename kuint>() {
            const std::vector<Emitted> recs {emit<kuint>(k, m, seq)};
            uint64_t emitted {0};
            for (const Emitted& e : recs) emitted += e.nucl.size() - k + 1;
            EXPECT_EQ(emitted, kmer_count(k, seq)) << "k=" << k << " m=" << m << " seq=" << seq;
        });
        expect_producer_contract(k, m, seq);
    }
}

// Homopolymers of every length around the window w = k-m+1 and around k, and di-/tri-/tetra-
// nucleotide tandem repeats, alone, flanked, and at either sequence end, over the whole grid. Poly-T
// and poly-G are covered as the reverse complements of poly-A and poly-C. The wide (k > 40) grid
// points, slow in Debug, skip the one-sided placements and the longest repeats.
TEST(SkmeratorDuplicateEmission, HomopolymersAndShortTandemRepeats) {
    const std::string left {"ACGTTGCA"}, right {"TTGCAGT"};
    for (const KM& g : kGrid) {
        const uint64_t k {g.k}, m {g.m}, w {k - m + 1};
        const bool wide {k > 40};
        for (const char c : {'A', 'C'})
            for (uint64_t len : {uint64_t{1}, uint64_t{3}, w - 1, w, w + 1, w + 2, k - 1, k, k + 1, 2 * k + 1}) {
                if (len == 0) continue;
                const std::string run(len, c);
                std::vector<std::string> seqs {left + run + right, run};
                if (!wide) { seqs.push_back(run + right); seqs.push_back(left + run); }
                for (const std::string& seq : seqs)
                    if (seq.size() >= k) expect_producer_contract(k, m, seq);
            }
        for (const std::string unit : {"AC", "AT", "GC", "CAG", "AAC", "ACGT", "AACG", "AATT"})
            for (uint64_t len : {k, 2 * k + 1, 3 * k}) {
                if (wide && len == 3 * k) continue;
                const std::string rep {tile(unit, len)};
                expect_producer_contract(k, m, left + rep + right);
                expect_producer_contract(k, m, rep);
            }
    }
}

// Randomized property test (fixed seeds): repeat-rich random sequences over the whole grid.
TEST(SkmeratorDuplicateEmission, RandomizedProperties) {
    std::mt19937_64 rng(20261011);
    auto urand = [&](uint64_t lo, uint64_t hi) { return std::uniform_int_distribution<uint64_t>(lo, hi)(rng); };
    auto random_seq = [&](uint64_t n, const std::string& alphabet) {
        std::string s(n, 'A');
        for (char& c : s) c = alphabet[urand(0, alphabet.size() - 1)];
        return s;
    };
    for (int rep {0}; rep < 6; ++rep) {
        for (const KM& g : kGrid) {
            const uint64_t k {g.k}, m {g.m};
            std::string seq {random_seq(urand(0, k), "ACGT")};
            // A random backbone with injected runs: homopolymers, short tandem repeats (period below
            // the window), RC-palindromes, and low-complexity stretches.
            for (uint64_t block {urand(1, 4)}; block > 0; --block) {
                switch (urand(0, 4)) {
                case 0: seq += std::string(urand(2, 2 * k), NUC[urand(0, 3)]); break;
                case 1: seq += tile(random_seq(urand(2, 4), "ACGT"), urand(k / 2, 3 * k)); break;
                case 2: seq += tile(random_seq(urand(1, std::max<uint64_t>(1, k - m)), "ACGT"), urand(k / 2, 3 * k)); break;
                case 3: { const std::string x {random_seq(urand(2, k), "ACGT")}; seq += x + revcomp(x); } break;
                default: seq += random_seq(urand(k, 3 * k), urand(0, 1) ? "AC" : "AT"); break;
                }
                seq += random_seq(urand(0, k / 2), "ACGT");
            }
            if (seq.size() < k) seq += random_seq(k - seq.size(), "ACGT");
            expect_producer_contract(k, m, seq);
        }
    }
}

// Contract (4) end to end: build an on-disk index from a repeat-rich FASTA and query that FASTA
// back with the CLI drivers. Each must print exactly one flag per input k-mer, all present.
TEST(SkmeratorDuplicateEmission, QueryPrintsOneFlagPerKmer) {
    std::vector<std::string> records {
        "ACGTTGCA" + std::string(40, 'A') + "TTGCAGT",
        "AACCCCCCCCCCC" + std::string(70, 'C') + "GATTACA",
        tile("CA", 90) + "GGT" + tile("AAT", 80),
        tile("ACGT", 100) + std::string(33, 'T'),
        "TTTGGGGCACT" + tile("GC", 64) + std::string(31, 'G') + tile("ATAT", 50),
    };
    const std::string fa {::testing::TempDir() + "bug08_query.fa"};
    {
        std::ofstream f(fa);
        for (size_t i {0}; i < records.size(); ++i) f << ">r" << i << "\n" << records[i] << "\n";
    }
    for (const KM& g : {KM{11, 5}, KM{31, 15}, KM{21, 4}, KM{63, 31}, KM{95, 41}}) {
        const uint64_t k {g.k}, m {g.m};
        uint64_t expected {0};
        for (const std::string& r : records) expected += kmer_count(k, r);

        const std::string idx {::testing::TempDir() + "bug08_query_k" + std::to_string(k) + ".sskm"};
        km::sortedlist::SortedListBuildParams p;
        p.k = k;
        p.m = m;
        p.input_path = fa;
        p.output_path = idx;
        p.ascii = false;
        p.buckets = 64;
        p.has_output_file = true;
        const uint64_t b {km::sortedlist::quotient_bits_for(p)};
        const uint64_t gen_w {km::sortedlist::select_generation_width_bytes(k, m)};
        const uint64_t store_w {km::sortedlist::select_width_bytes(2 * (2 * k - m) - b)};
        km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
            km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
                km::sortedlist::build_sorted_list<gen, store>(p, b);
                auto reader {km::sortedlist::BucketedSkmerListReader<store>::open(idx)};
                std::stringstream seq_out, par_out;
                km::sortedlist::sequential_query<gen, store>(reader, fa, seq_out);
                km::sortedlist::parallel_query<gen, store>(reader, fa, par_out, 4);
                for (const std::string& out : {seq_out.str(), par_out.str()}) {
                    const uint64_t ones = static_cast<uint64_t>(std::count(out.begin(), out.end(), '1'));
                    const uint64_t zeros = static_cast<uint64_t>(std::count(out.begin(), out.end(), '0'));
                    EXPECT_EQ(ones + zeros, expected) << "flags printed for " << expected << " k-mers, k=" << k << " m=" << m;
                    EXPECT_EQ(zeros, 0u) << "self-query reported absent k-mers, k=" << k << " m=" << m;
                }
                EXPECT_EQ(seq_out.str(), par_out.str()) << "k=" << k << " m=" << m;
            });
        });
    }
}
