// Regression test for issue #6:
//   `sskm construct` silently drops a few low-complexity (poly-A, all-A
//   minimizer) k-mers from the sorted list on large repeat-rich input, even at
//   k = 21 / m = 11 (the query-safe regime).

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

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
