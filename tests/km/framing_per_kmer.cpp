// TDD (red-first) regression tests for per-k-mer canonical framing — the
// completion of issue #7. A k-mer present in the list MUST be queryable,
// regardless of the surrounding context that framed its super-k-mer. This
// breaks today when a k-mer's (φ-)minimizer is a reverse-complement palindrome
// and/or repeats within the k-mer: the whole-super-k-mer canonicalize (the #7
// fix) cannot make every contained k-mer canonical, so the bare-k-mer query
// (canonical) misses it. These tests are expected RED until the construction
// stores each k-mer in its own canonical frame; then they guard against
// regression.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {

using kuint = uint64_t;

bool kmer_is_queryable(km::sortedlist::SortedVirtualSkmerList<kuint>& list,
                       km::SkmerManipulator<kuint>& manip, std::string kmer) {
    km::SeqSkmerator<kuint> qg{manip, kmer};
    std::vector<km::Skmer<kuint>> query;
    for (const km::Skmer<kuint>& s : qg) query.push_back(s);
    for (const std::vector<uint8_t>& row : list.query_skmer_batch(query))
        for (uint8_t v : row)
            if (v != 0) return true;
    return false;
}

// Coverage invariant: build the list from `seq`, then every length-k substring
// of `seq` (a k-mer that IS present) must be queryable as a bare k-mer.
void expect_all_kmers_queryable(uint64_t k, uint64_t m, std::string seq) {
    km::SkmerManipulator<kuint> manip{k, m};
    km::SeqSkmerator<kuint> gen{manip, seq};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& s : gen) enumeration.push_back(s);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    for (size_t i = 0; i + k <= seq.size(); i++) {
        std::string kmer = seq.substr(i, k);
        EXPECT_TRUE(kmer_is_queryable(list, manip, kmer))
            << "k=" << k << " m=" << m << " seq=" << seq
            << " kmer[" << i << "]=" << kmer << " is present but not queryable";
    }
}

}  // namespace

// 2-periodic repeats produce reverse-complement-palindrome minimizers (e.g. GC,
// AT) that repeat within each k-mer — the exact trigger of the residual.
TEST(FramingPerKmer, PalindromeRepeatMinimizers_k5m2) {
    expect_all_kmers_queryable(5, 2, "GCGCGA");
    expect_all_kmers_queryable(5, 2, "CGCGCGCG");
    expect_all_kmers_queryable(5, 2, "ATATATATA");
    expect_all_kmers_queryable(5, 2, "GCGCGCGCGC");
}

TEST(FramingPerKmer, PalindromeRepeatMinimizers_larger) {
    expect_all_kmers_queryable(7, 3, "GCGCGCGCGCGC");
    expect_all_kmers_queryable(9, 4, "ATATATATATATAT");
    expect_all_kmers_queryable(7, 3, "ACGCGCGCGTA");
}

// Mixed / less structured sequences should also round-trip with 0 false negatives.
TEST(FramingPerKmer, MixedSequences) {
    expect_all_kmers_queryable(5, 2, "ACGTACGTACGT");
    expect_all_kmers_queryable(7, 3, "ACGTTGCAACGTTGCA");
    expect_all_kmers_queryable(11, 5, "ACGTACGTACGTACGTACGT");
}
