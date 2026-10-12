// Regression test for issue #7:
//   A k-mer that IS present in the list must be queryable. When a k-mer's
//   minimizer occurs at several positions (repeat / homopolymer), construction
//   and a bare-k-mer query can frame the super-k-mer differently (different
//   m_pref_size/m_suff_size, i.e. different column), and query_skmer only
//   compares k-mers at matching columns -> false negative.
//
// Minimal case (k=5, m=2): construct from "GCGCGA"; the list stores GCGCG (the
// ASCII dump shows "GCGCG 2 1" and a KMC set-comparison shows no construction
// drop), but querying the bare k-mer GCGCG currently returns "not found".
//
// This test currently FAILS and should pass once the framing is made consistent.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

TEST(QueryContextFraming, PresentKmerIsQueryable) {
    using kuint = uint64_t;
    const uint64_t k = 5, m = 2;
    km::SkmerManipulator<kuint> manip{k, m};

    std::string genome = "GCGCGA";
    km::SeqSkmerator<kuint> gen{manip, genome};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& s : gen) enumeration.push_back(s);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    // GCGCG occurs in the genome and is stored in the list (ASCII: "GCGCG 2 1").
    std::string kmer = "GCGCG";
    km::SeqSkmerator<kuint> qg{manip, kmer};
    std::vector<km::Skmer<kuint>> query;
    for (const km::Skmer<kuint>& s : qg) query.push_back(s);

    bool found = false;
    for (const std::vector<uint8_t>& row : list.query_skmer_batch(query))
        for (uint8_t v : row) found = found || (v != 0);

    EXPECT_TRUE(found) << "GCGCG is present in the list but query reports it absent "
                          "(context-dependent super-k-mer framing)";
}

// The same context-free framing under the open-closed scheme (m >= 3): build from sequences whose
// minimizer repeats (and, for even m, is an RC-palindrome); every present k-mer, queried bare, is found.
TEST(QueryContextFraming, PresentKmerIsQueryable_OpenClosed) {
    using kuint = uint64_t;
    for (const auto& [k, m, genome] : std::vector<std::tuple<uint64_t, uint64_t, std::string>>{
             {7, 3, "GCGCGCGA"}, {9, 4, "GCGCGCGCGA"}, {11, 5, "AAAAAAAAAAAAACG"}, {9, 3, "ACACACACACACAT"}}) {
        km::SkmerManipulator<kuint> manip{k, m, 0, km::SchemeParams::open_closed(m)};
        std::string g = genome;
        km::SeqSkmerator<kuint> gen{manip, g};
        std::vector<km::Skmer<kuint>> enumeration;
        for (const km::Skmer<kuint>& s : gen) enumeration.push_back(s);

        km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
        list.generate_sorted_list_from_enumeration(enumeration);

        for (size_t i = 0; i + k <= genome.size(); ++i) {
            std::string kmer = genome.substr(i, k);
            km::SeqSkmerator<kuint> qg{manip, kmer};
            std::vector<km::Skmer<kuint>> query;
            for (const km::Skmer<kuint>& s : qg) query.push_back(s);
            bool found = false;
            for (const std::vector<uint8_t>& row : list.query_skmer_batch(query))
                for (uint8_t v : row) found = found || (v != 0);
            EXPECT_TRUE(found) << genome.substr(i, k) << " (k=" << k << ", m=" << m
                               << ", open-closed) is present but reported absent";
        }
    }
}
