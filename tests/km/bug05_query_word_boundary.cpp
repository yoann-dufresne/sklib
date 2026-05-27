// Regression test for issue #5:
//   `sskm query` returns false positives once the packed super-k-mer spans past
//   the first 64-bit word of Skmer::pair.  With m = 11:
//     k = 22  ->  2*(2k-m) = 66 bits  (fits in one word)      -> correct
//     k = 23  ->  2*(2k-m) = 70 bits  (spills into m_value[1]) -> false positives
//   Construction and true-positive lookups are correct; only *absent* k-mers are
//   wrongly reported present, so this is invisible to self-query / set-equality.
//
// This test currently FAILS (reports 2 false positives) and should pass once the
// high-word handling in the query comparison is fixed.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {

// Build a sorted list from `genome`, then count how many of `absent_kmers`
// (each of length k, none occurring in `genome`) the query wrongly reports found.
int query_false_positives(uint64_t k, uint64_t m,
                          std::string genome,
                          std::vector<std::string> absent_kmers) {
    using kuint = uint64_t;
    km::SkmerManipulator<kuint> manip{k, m};

    km::SeqSkmerator<kuint> gen{manip, genome};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& s : gen) enumeration.push_back(s);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    std::vector<km::Skmer<kuint>> queries;  // each k-long record -> 1 skmer / 1 k-mer
    for (std::string& kmer : absent_kmers) {
        km::SeqSkmerator<kuint> qg{manip, kmer};
        for (const km::Skmer<kuint>& s : qg) queries.push_back(s);
    }

    int false_positives = 0;
    for (const std::vector<uint8_t>& row : list.query_skmer_batch(queries))
        for (uint8_t found : row) false_positives += (found != 0);
    return false_positives;
}

}  // namespace

TEST(QueryWordBoundary, NoFalsePositivesAcrossSecondWord) {
    // 40 bp reference; the 16 query k-mers below are all length 23 and occur in
    // neither this sequence nor its reverse complement.
    std::string genome = "ACGTTGCATTACGGCATTGACCGTATGCATGCATGCTTAA";
    std::vector<std::string> absent_23 = {
        "GCTAAAGACAATTACATAACATA", "CACGTCAGCACGAAACTTGTTGG", "CCCAGTGTGAATCGCTTAAGGGT",
        "TAAGTAAGTGTGATGCATACGCC", "TTTACTTGCTGTGTCCACCCCAT", "CGGACTGGCATTTTTATTACACT",
        "CAGAAACAGAACTCGGGTAATTT", "TGACAGGTCACGCAGAGGCGCGC", "CCTCCTGAAGTGCGTGGACACTC",
        "GCTATGAATCTCTGATTTACCCA", "CTCTGCCAAACTCCAGCGCGGTC", "AGTTCCATCACCCTAAGTAACCG",
        "AATAATGCGTTCGCTCTATTGAC", "TACGACGCGCTCATTCCCTTGTC", "GGAGAGTTATGGAACAAGGACGC",
        "TGTCTGAGACTAGAAGACAGATA",
    };

    EXPECT_EQ(query_false_positives(23, 11, genome, absent_23), 0)
        << "absent 23-mers reported present (k=23, m=11 -> 70-bit super-k-mer "
           "crosses the 64-bit word boundary)";
}
