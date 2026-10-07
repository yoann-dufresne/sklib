// Regression test: the construction compaction extended a CLOSED virtual super-k-mer.
//
// One raw super-k-mer contributes k-mers to several columns under the same enumeration id. When a
// virtual super-k-mer stopped growing, the tail loops of merge_LList_column passed it through
// without closing it (expandable stayed true). A later column's overlap whose left k-mer had that
// same id then matched this stale record first (is_left compares last_id only), extended it past a
// gap of columns and left it at a list position that was no longer sorted for its new columns. The
// k-mer set stayed exact, but the per-column sorted invariant broke, so the set-op merge missed
// shared k-mers (6 of 844 M on Tara ocean reads, k=63, --buckets 262144).
//
// The 18 raw super-k-mers below are the delta-debugged reduction of that bucket (k=63, m=31, stored
// representation with the ψ slot). With greedy_chaining (the construction default) the compaction
// produced two records out of order at columns 25 and 27; colinear_chaining happened to avoid it.

#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

#include <io/Skmer.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {

using kuint = __uint128_t;
constexpr uint64_t K = 63, M = 31;

struct RawSkmer { uint64_t hi, lo_hi, lo_lo; uint16_t pref, suff; };

const RawSkmer kBucket[] = {
    {0x399bbf9800549480ULL, 0xc524904b331f8ccbULL, 0xbb7b333f7b3733fbULL, 15, 32},
    {0x399bbf98eca495c0ULL, 0x29721dee90a96e3aULL, 0xa5a95ad7619aa90aULL, 32, 32},
    {0x399bbf9bc409eb00ULL, 0x03228e0784526301ULL, 0x182f227a0829ceceULL, 32, 28},
    {0x399bbfa7fbb5d798ULL, 0x976e16773f7ff7b7ULL, 0x7bbff7b3bbb7f3ffULL, 6, 30},
    {0x399bbfdfe5104ba0ULL, 0xbdaa96a921aaa0b8ULL, 0x8fb9fd8e0c7bbbbfULL, 26, 32},
    {0x399bbfdff4312b40ULL, 0x33b9826e2f7a66aeULL, 0xcb3f7377bbf7bbb7ULL, 19, 32},
    {0x399bbff44f442c20ULL, 0xb48bf474afade37bULL, 0xe502682a7c560236ULL, 32, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65eba3bbbULL, 0xbf777733bbfb3777ULL, 12, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28a8ULL, 0xbf5745339be81667ULL, 31, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28a8ULL, 0xbf5745339be8166fULL, 31, 31},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28a8ULL, 0xbf5745339be81eefULL, 31, 29},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28a8ULL, 0xbf5745339be83777ULL, 28, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28a8ULL, 0xbf5745339be8deefULL, 31, 28},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28a8ULL, 0xbf5745379fecdeefULL, 31, 25},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb28abULL, 0xbf777733bbfb3777ULL, 15, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb3bbbULL, 0xbf777733bbbb3777ULL, 11, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d65ebb3bbbULL, 0xbf777733bbfb3777ULL, 11, 32},
    {0x399bbff876a8bd50ULL, 0xebb8f1d75fbb3bbbULL, 0xbf777733bbfb3777ULL, 9, 32},
};

std::vector<km::Skmer<kuint>> bucket_enumeration() {
    std::vector<km::Skmer<kuint>> e;
    for (const RawSkmer& r : kBucket) {
        const kuint lo {(static_cast<kuint>(r.lo_hi) << 64) | r.lo_lo};
        e.emplace_back(typename km::Skmer<kuint>::pair(lo, static_cast<kuint>(r.hi)), r.pref, r.suff);
    }
    km::sortedlist::sort_and_dedup(e);   // as BucketCompactor does before compacting
    return e;
}

}  // namespace

// At every column, the records holding a k-mer there must come in strictly increasing k-mer order:
// merge_columns and search_kmers_in_span both rely on it.
TEST(StaleExpandableVskmer, CompactionKeepsPerColumnOrder) {
    const std::vector<km::Skmer<kuint>> enumeration {bucket_enumeration()};
    km::SkmerManipulator<kuint> manip{K, M};
    for (const bool greedy : {true, false}) {
        km::sortedlist::SortedVirtualSkmerList<kuint> list(K, M);
        list.generate_sorted_list_from_enumeration(enumeration, greedy);
        const std::vector<km::Skmer<kuint>>& recs {list.get_list()};
        for (uint64_t c {0}; c <= K - M; ++c) {
            const km::Skmer<kuint>* prev {nullptr};
            for (const km::Skmer<kuint>& r : recs) {
                if (!manip.has_valid_kmer(r, c)) continue;
                if (prev) EXPECT_LT(manip.kmer_compare(*prev, r, c), 0)
                    << "column " << c << " out of order (greedy=" << greedy << ")";
                prev = &r;
            }
        }
    }
}

// Every k-mer of the input must be found again by a query on the compacted list.
TEST(StaleExpandableVskmer, EveryInputKmerIsFound) {
    const std::vector<km::Skmer<kuint>> enumeration {bucket_enumeration()};
    for (const bool greedy : {true, false}) {
        km::sortedlist::SortedVirtualSkmerList<kuint> list(K, M);
        list.generate_sorted_list_from_enumeration(enumeration, greedy);
        for (const km::Skmer<kuint>& s : enumeration)
            for (const uint8_t found : list.query_skmer(s))
                EXPECT_EQ(found, 1) << "k-mer of an input super-k-mer not found (greedy=" << greedy << ")";
    }
}
