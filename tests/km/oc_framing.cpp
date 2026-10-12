// Producer-level tests of the minimizer schemes (random and open-closed), on the raw super-k-mer
// stream the Skmerator yields (the stream both construction and query consume).
//
//  7. Exactly once, in order: decoding every yielded record back to the input, the records tile the
//     k-mer starts [0, n) contiguously — no k-mer missing, none emitted twice, no backward step.
//  8. Framing = the reference rule, and forwardness: every k-mer's minimizer occurrence and
//     orientation in the stream equal what the per-k-mer reference (choose_kmer_minimizer over the
//     memoized mmer_canonical_rank, i.e. the ambiguous path's class_of) selects — so the fast path's
//     STREAMING class agrees with the reference class at every position that matters — and the
//     selected occurrence never moves left.
//  9. Density smoke test: open-closed at (31,15) on 1 Mbp of random sequence gives ~0.0917 super-
//     k-mers per k-mer (random: ~0.111). A silently broken class falls back to ~0.111.
//
// The "Repeats" variants run on homopolymers and short tandem repeats, where a minimizer repeats
// >= 3 times in a window. They exercise the producer's duplicate-emission defect (notes/minimizers/
// 00-sklib-current-scheme.md §10, a separate fix in recompute_minimizer) and FAIL until that fix is
// in: a k-mer of the repeat is emitted in two records. The defect is scheme-independent.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <io/MinimizerScheme.hpp>
#include <algorithms/WidthDispatch.hpp>

#include "scheme_test_utils.hpp"

namespace {

using sktest::revcomp;

// Walk the stream of `seq` and check (7) + (8). Returns the number of records.
//
// Order: raw super-k-mers come out in sequence order, but the pieces an ambiguous one is split into
// (canonical_pieces) follow its CANONICAL orientation, i.e. reverse sequence order when that is the
// reverse complement. Pieces of one raw super-k-mer share its creation index (TrackPos), so the
// check is: each group of equal creation index tiles the next k-mers contiguously, in any order
// inside the group; groups follow each other in sequence order.
template<typename kuint>
uint64_t check_stream(uint64_t k, uint64_t m, const km::SchemeParams& sp, const std::string& seq,
                      const std::string& ctx) {
    km::SkmerManipulator<kuint> manip{k, m, 0, sp};
    std::string s = seq;
    km::SeqSkmerator<kuint, true> gen{manip, s};

    // Reference framing of every k-mer of the whole sequence (strand-invariant per-k-mer rule).
    std::vector<kuint> S(seq.size());
    for (size_t i = 0; i < seq.size(); ++i) S[i] = static_cast<kuint>((seq[i] >> 1) & 0b11);
    std::vector<kuint> rank_table(seq.size() - m + 1);
    for (uint64_t p = 0; p < rank_table.size(); ++p) rank_table[p] = manip.mmer_canonical_rank(S, p);

    struct Rec { int64_t created; std::string nts; uint64_t pref, suff; };
    std::vector<Rec> recs;
    for (auto it = gen.begin(); !it.consumed(); ++it) {
        km::Skmer<kuint> d = *it;
        manip.unpermute_minimizer_slot(d);                    // stored slot is ψ(min): decode the raw m-mer
        std::string Ds;
        for (kuint x : manip.decode_to_nucleotides(d)) Ds.push_back("ACTG"[static_cast<int>(x)]);
        recs.push_back({it.yielded_position(), Ds, d.m_pref_size, d.m_suff_size});
    }

    const uint64_t nk = seq.size() - k + 1;
    uint64_t next = 0;          // next k-mer start the stream must cover
    int64_t last_occ = -1;
    for (size_t g0 = 0; g0 < recs.size();) {
        size_t g1 = g0;
        while (g1 < recs.size() && recs[g1].created == recs[g0].created) ++g1;
        std::vector<bool> used(g1 - g0, false);
        for (size_t placed = 0; placed < g1 - g0; ++placed) {
            // the piece of this group that starts exactly at `next`
            size_t pick = g1;
            bool is_fwd = false, is_rc = false;
            for (size_t r = g0; r < g1 && pick == g1; ++r) {
                if (used[r - g0] || next + recs[r].nts.size() > seq.size()) continue;
                const std::string fwd = seq.substr(next, recs[r].nts.size());
                is_fwd = (recs[r].nts == fwd);
                is_rc = (recs[r].nts == revcomp(fwd));
                if (is_fwd || is_rc) pick = r;
            }
            if (pick == g1) {
                ADD_FAILURE() << ctx << ": no record of group #" << g0 << " (" << (g1 - g0) << " piece(s), first "
                              << recs[g0].nts << ") covers the next k-mer " << next << " (" << seq.substr(next, k)
                              << ") -- a k-mer emitted twice, missing, or out of order";
                return recs.size();
            }
            used[pick - g0] = true;
            const Rec& rec = recs[pick];
            const uint64_t L = rec.nts.size();
            EXPECT_GE(L, k) << ctx;
            const uint64_t occ_f = next + rec.pref, occ_r = next + rec.suff;
            for (uint64_t t = next; t + k <= next + L; ++t) {
                uint64_t occ;
                km::orientation_t o;
                manip.choose_kmer_minimizer(S, t, rank_table, occ, o);
                if (is_fwd && is_rc) {
                    EXPECT_TRUE(occ == occ_f || occ == occ_r) << ctx << " palindromic record, k-mer " << t;
                    continue;
                }
                const std::string kmer = seq.substr(t, k);
                const uint64_t got = is_fwd ? occ_f : occ_r;
                if (kmer == revcomp(kmer)) {
                    // RC-palindromic k-mer (even k): the frames at j and k-m-j are each other's reverse
                    // complement, i.e. the same canonical record; either occurrence is correct.
                    EXPECT_TRUE(got == occ || got == t + (k - m) - (occ - t)) << ctx << ": palindromic k-mer " << t;
                    continue;
                }
                EXPECT_EQ(occ, got) << ctx << ": k-mer " << t << " framed at a different "
                    "minimizer occurrence than the reference rule";
                EXPECT_EQ(o, is_fwd ? km::forward_c : km::reverse_c) << ctx << ": k-mer " << t << " orientation";
            }
            // Forwardness, in sequence order (the pieces are placed in sequence order here).
            if (!(is_fwd && is_rc)) {
                const int64_t occ = static_cast<int64_t>(is_fwd ? occ_f : occ_r);
                EXPECT_GE(occ, last_occ) << ctx << ": backward jump of the selected minimizer at k-mer " << next;
                last_occ = occ;
            }
            next += L - k + 1;
        }
        g0 = g1;
    }
    EXPECT_EQ(next, nk) << ctx << ": the stream stops before the last k-mer";
    return recs.size();
}

uint64_t check_stream_dispatch(uint64_t k, uint64_t m, const km::SchemeParams& sp, const std::string& seq) {
    const std::string ctx = "k=" + std::to_string(k) + " m=" + std::to_string(m) + " " + sp.to_string();
    uint64_t n = 0;
    km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m, sp),
        [&]<typename kuint>() { n = check_stream<kuint>(k, m, sp, seq, ctx); });
    return n;
}

std::vector<km::SchemeParams> both_schemes(uint64_t m) {
    std::vector<km::SchemeParams> v {km::SchemeParams{}};
    if (km::SchemeParams::default_s(m) != 0) v.push_back(km::SchemeParams::open_closed(m));
    return v;
}

struct KM { uint64_t k, m; };
// Every generation width: uint32 (21,11 / 15,7), uint64 (31,15), __uint128 (63,31 / 40,32 under OC),
// kuint256 (90,41). (40,32) is the one where open-closed bumps the generation width (66-bit rank).
const KM kGrid[] = {{21, 11}, {15, 7}, {31, 15}, {25, 9}, {32, 16}, {63, 31}, {40, 32}, {90, 41}};

} // namespace

// (7) + (8) on random sequence: every k-mer exactly once, in order, framed by the reference rule.
TEST(SchemeFraming, ExactlyOnceInOrderAndReferenceFraming_Random) {
    for (const KM& c : kGrid)
        for (const km::SchemeParams& sp : both_schemes(c.m))
            for (uint64_t seed : {11ull, 12ull}) {
                const std::string seq = sktest::random_dna(6000, seed * 1000 + c.k);
                check_stream_dispatch(c.k, c.m, sp, seq);
            }
}

// Same, with explicit expert s values (smallest, a no-table s > 8, even m) and several reads.
TEST(SchemeFraming, ExactlyOnceInOrderAndReferenceFraming_ExpertS) {
    const std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> cases = {
        {31, 15, 1}, {31, 15, 5}, {31, 15, 13}, {41, 21, 9}, {31, 16, 2}, {31, 16, 10}, {11, 5, 1}};
    for (const auto& [k, m, s] : cases)
        for (uint64_t seed : {21ull, 22ull, 23ull})
            check_stream_dispatch(k, m, km::SchemeParams::open_closed(m, s), sktest::random_dna(3000, seed + k + s));
}

// Repeat-rich input: the ambiguous path (repeated minimizers) and the streaming class on runs. The
// framing / strand checks hold today; the exactly-once part fails until the duplicate-emission fix
// in recompute_minimizer lands (see the file header).
TEST(SchemeFraming, ExactlyOnceInOrderAndReferenceFraming_Repeats) {
    const std::vector<std::string> units = {"A", "C", "AC", "AG", "ACG", "AAC", "ACGT", "AACG", "ATCGG"};
    for (const KM& c : {KM{21, 11}, KM{15, 7}, KM{31, 15}, KM{11, 5}, KM{12, 6}})
        for (const km::SchemeParams& sp : both_schemes(c.m))
            for (const std::string& u : units) {
                const std::string seq = sktest::random_dna(150, c.k + u.size()) + sktest::tile(u, 120) +
                                        sktest::random_dna(80, 7 * c.k + u.size()) + sktest::tile(u, 3 * c.k) +
                                        sktest::random_dna(60, 3 * c.k);
                check_stream_dispatch(c.k, c.m, sp, seq);
            }
}

// (9) Density smoke test on 1 Mbp of uniform random sequence.
TEST(SchemeFraming, OpenClosedDensitySmoke) {
    const uint64_t k = 31, m = 15;
    std::string seq = sktest::random_dna(1000000, 2026);
    const double nk = static_cast<double>(seq.size() - k + 1);
    auto count = [&](const km::SchemeParams& sp) {
        uint64_t n = 0;
        km::sortedlist::dispatch_width_bytes(km::sortedlist::select_generation_width_bytes(k, m, sp),
            [&]<typename kuint>() {
                km::SkmerManipulator<kuint> manip{k, m, 0, sp};
                km::SeqSkmerator<kuint> gen{manip, seq};
                for (const km::Skmer<kuint>& r : gen) { (void)r; ++n; }
            });
        return static_cast<double>(n) / nk;
    };
    const double d_rnd = count(km::SchemeParams{});
    const double d_oc  = count(km::SchemeParams::open_closed(m));   // s = 3
    std::printf("DENSITY k=31 m=15 random=%.5f open-closed(s=3)=%.5f (%.1f%%)\n", d_rnd, d_oc,
                100.0 * (d_oc / d_rnd - 1.0));
    EXPECT_NEAR(d_rnd, 0.1111, 0.1111 * 0.02);
    EXPECT_NEAR(d_oc, 0.0917, 0.0917 * 0.02);
}
