// Regression test for issue #6:
//   `sskm construct` silently drops a few low-complexity (poly-A, all-A
//   minimizer) k-mers from the sorted list on large repeat-rich input, even at
//   k = 21 / m = 11 (the query-safe regime).
//
// The bug needs the global arrangement of a large same-minimizer group, so it
// only reproduces at scale (~1 Mb of real repeat-rich sequence) and is therefore
// fixture-driven: point SKLIB_BUG2_FASTA at e.g. the first ~2 Mb of hg38 chr21's
// q-arm, sanitized to uppercase ACGT (split at non-ACGT runs).  Without the
// fixture the test is skipped.
//
// With that fixture this test currently FAILS: the 21-mer below occurs in the
// input but is missing from the constructed list.

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

TEST(ConstructHomopolymerDrop, ListContainsKnownPolyAKmer) {
    const char* env = std::getenv("SKLIB_BUG2_FASTA");
    if (env == nullptr) {
        GTEST_SKIP() << "set SKLIB_BUG2_FASTA to a >=1Mb repeat-rich ACGT FASTA "
                        "(e.g. first ~2 Mb of hg38 chr21 q-arm, sanitized)";
    }

    using kuint = uint64_t;
    const uint64_t k = 21, m = 11;  // 2*(2k-m) = 62 bits -> one word, query is reliable here
    std::string path(env);

    km::SkmerManipulator<kuint> manip{k, m};
    km::FileSkmerator<kuint> file_skmerator{manip, path};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& s : file_skmerator) enumeration.push_back(s);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    // Low-complexity 21-mer that occurs in the fixture but is dropped from the
    // list.  (For a different fixture, the exact dropped k-mer may differ; the
    // robust check is a full k-mer-set comparison, as scripts/large_scale_e2e.sh
    // does against KMC.)
    std::string kmer = "AAAAAAAAAAAACAAAAAAAA";
    km::SeqSkmerator<kuint> qg{manip, kmer};
    std::vector<km::Skmer<kuint>> query;
    for (const km::Skmer<kuint>& s : qg) query.push_back(s);

    std::vector<std::vector<uint8_t>> res = list.query_skmer_batch(query);
    ASSERT_EQ(res.size(), 1u);
    ASSERT_EQ(res[0].size(), 1u);
    EXPECT_EQ(res[0][0], 1)
        << "poly-A k-mer present in the input but missing from the constructed list";
}

// Self-contained reduction of #6 (no fixture): construct from a 24 bp sequence at
// k=5/m=2; the k-mer AAACA is in the input but is dropped from the list because
// colinear_chaining returns invalid overlaps for this column. Currently FAILS.
TEST(ConstructHomopolymerDrop, SmallSelfContainedRepro) {
    using kuint = uint64_t;
    const uint64_t k = 5, m = 2;
    km::SkmerManipulator<kuint> manip{k, m};

    std::string genome = "GGAGCCAAACAGGAGGAAAAGAGG";
    km::SeqSkmerator<kuint> gen{manip, genome};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& s : gen) enumeration.push_back(s);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    std::string kmer = "AAACA";  // present in genome (GGAGCC[AAACA]GGAGG...)
    km::SeqSkmerator<kuint> qg{manip, kmer};
    std::vector<km::Skmer<kuint>> query;
    for (const km::Skmer<kuint>& s : qg) query.push_back(s);

    bool found = false;
    for (const std::vector<uint8_t>& row : list.query_skmer_batch(query))
        for (uint8_t v : row) found = found || (v != 0);

    EXPECT_TRUE(found) << "AAACA is in the input but dropped from the list (issue #6)";
}
