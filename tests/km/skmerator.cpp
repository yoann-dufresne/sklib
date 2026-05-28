#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

using namespace std;

namespace {
// Coverage invariant (order/φ/split-agnostic): every k-mer of `seq`, once built into
// the list, must be queryable as a bare k-mer. This replaces the older exact-value /
// exact-count assertions, which encoded the raw-minimizer enumeration and broke once
// the order became φ-based and palindrome-minimizer super-k-mers are split per k-mer.
template <typename kuint>
void expect_covers_all_kmers(uint64_t k, uint64_t m, std::string seq) {
    km::SkmerManipulator<kuint> manip{k, m};
    km::SeqSkmerator<kuint> gen{manip, seq};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& s : gen) enumeration.push_back(s);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    for (size_t i = 0; i + k <= seq.size(); i++) {
        std::string kmer = seq.substr(i, k);
        km::SeqSkmerator<kuint> qg{manip, kmer};
        std::vector<km::Skmer<kuint>> q;
        for (const km::Skmer<kuint>& s : qg) q.push_back(s);
        bool found = false;
        for (const std::vector<uint8_t>& row : list.query_skmer_batch(q))
            for (uint8_t v : row) found = found || (v != 0);
        EXPECT_TRUE(found) << "k=" << k << " m=" << m << " seq=" << seq
                           << " kmer[" << i << "]=" << kmer << " present but not queryable";
    }
}
}  // namespace


// These sequences exercise specific Skmerator paths (equal consecutive minimizers,
// decreasing/increasing minimizer, out-of-context minimizer, homopolymers). The old
// assertions pinned the exact raw-minimizer enumeration (m_pair / pref / suff / counts),
// which is no longer meaningful now the order is φ-based and palindrome-minimizer
// super-k-mers are split per k-mer. We assert the order-agnostic coverage invariant
// instead: every k-mer of the sequence is queryable.
TEST(Skmerator, update_equal_mini_fwd_fwd) { expect_covers_all_kmers<uint16_t>(5, 2, "CCCCCC"); }
TEST(Skmerator, update_equal_mini_rev_rev) { expect_covers_all_kmers<uint16_t>(5, 2, "GGGGGG"); }
TEST(Skmerator, decreasing_minimizer)      { expect_covers_all_kmers<uint16_t>(5, 2, "CCCCAAAAA"); }
TEST(Skmerator, increasing_minimizer)      { expect_covers_all_kmers<uint16_t>(5, 2, "AACCCC"); }
TEST(Skmerator, outofcontext_minimizer)    { expect_covers_all_kmers<uint16_t>(8, 2, "AACAATAAGGGGGGG"); }

TEST(Skmerator, seq_test_5_2)
{
    std::string seq{"ATCGACTGTGTACACT"};
    for (uint64_t n{5} ; n<=seq.length() ; n++)
        expect_covers_all_kmers<uint16_t>(5, 2, seq.substr(0, n));
}


TEST(Skmerator, file_vs_seq)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};
    km::SkmerPrettyPrinter<kuint> pp {k, m};

    // --- Sequence ---
    std::string seq{"CCCCAAAAA"};
    km::SkmerManipulator<kuint> seq_manip {k, m};
    km::SeqSkmerator<kuint> seq_skmerator {seq_manip, seq};

    // Enumerates the superkmers from the sequence
    std::vector<km::Skmer<kuint> > seq_skmers {};
    for (km::Skmer<kuint> skmer : seq_skmerator)
        seq_skmers.push_back(skmer);

    // --- File ---
    std::string filename{"../tests/data/CCCCAAAAA.fa"};
    km::SkmerManipulator<kuint> file_manip {k, m};
    km::FileSkmerator<kuint> file_skmerator {file_manip, filename};


    // Enumerates the superkmers from the file
    std::vector<km::Skmer<kuint> > file_skmers {};
    for (km::Skmer<kuint> skmer : file_skmerator)
        file_skmers.push_back(skmer);

    // Comparison of size
    if (seq_skmers.size() != file_skmers.size())
    {
        std::cerr << "from sequence" << std::endl;
        for (const auto& skmer : seq_skmers)
        {
            pp << skmer;
            std::cerr << pp << " ";
        } std::cerr << std::endl;

        std::cerr << "from file" << std::endl;
        for (const auto& skmer : file_skmers)
        {
            pp << skmer;
            std::cerr << pp << " ";
        } std::cerr << std::endl;
    }
    ASSERT_EQ(seq_skmers.size(), file_skmers.size());

    // Compare skmers one by one
    for (uint64_t i{0} ; i<file_skmers.size() ; i++)
        ASSERT_TRUE(seq_manip.skmer_equals(seq_skmers[i], file_skmers[i]));
}

TEST(Skmerator, file_vs_seq1)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};
    km::SkmerPrettyPrinter<kuint> pp {k, m};

    // --- Sequence ---
    std::string seq{"ATCGACTGTGTACACT"};
    km::SkmerManipulator<kuint> seq_manip {k, m};
    km::SeqSkmerator<kuint> seq_skmerator {seq_manip, seq};

    // Enumerates the superkmers from the sequence
    std::vector<km::Skmer<kuint> > seq_skmers {};
    for (km::Skmer<kuint> const skmer : seq_skmerator)
        seq_skmers.emplace_back(skmer);

    // --- File ---
    std::string filename{"../tests/data/fasta00.fa"};
    km::SkmerManipulator<kuint> file_manip {k, m};
    km::FileSkmerator<kuint> file_skmerator {file_manip, filename};


    // Enumerates the superkmers from the file
    std::vector<km::Skmer<kuint> > file_skmers {};
    for (km::Skmer<kuint> const skmer : file_skmerator) {
        file_skmers.emplace_back(skmer);
    }

    // Comparison of size
    if (seq_skmers.size() != file_skmers.size())
    {
        std::cerr << "from sequence" << std::endl;
        for (const auto& skmer : seq_skmers)
        {
            pp << skmer;
            std::cerr << pp << " ";
        } std::cerr << std::endl;

        std::cerr << "from file" << std::endl;
        for (const auto& skmer : file_skmers)
        {
            pp << skmer;
            std::cerr << pp << " ";
        } std::cerr << std::endl;
    }
    ASSERT_EQ(seq_skmers.size(), file_skmers.size());

    // Compare skmers one by one
    for (uint64_t i{0} ; i<file_skmers.size() ; i++)
        ASSERT_TRUE(seq_manip.skmer_equals(seq_skmers[i], file_skmers[i]));
}


// ------------------------- BOUNDARY LENGTH TESTS -------------------------

TEST(Skmerator, empty_sequence)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{""};
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    uint64_t nb_skmer {0};
    for ([[maybe_unused]] km::Skmer<kuint> skmer : skmerator)
        nb_skmer += 1;

    EXPECT_EQ(nb_skmer, 0U);

    auto it = skmerator.begin();
    EXPECT_TRUE(it.consumed());
    EXPECT_TRUE(it == skmerator.end());
}


TEST(Skmerator, seq_shorter_than_k)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"AC"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    uint64_t nb_skmer {0};
    for ([[maybe_unused]] km::Skmer<kuint> skmer : skmerator)
        nb_skmer += 1;

    EXPECT_EQ(nb_skmer, 0U);
    EXPECT_TRUE(skmerator.begin().consumed());
}


TEST(Skmerator, seq_length_k_minus_one)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"ACGT"};  // length = k - 1
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    uint64_t nb_skmer {0};
    for ([[maybe_unused]] km::Skmer<kuint> skmer : skmerator)
        nb_skmer += 1;

    EXPECT_EQ(nb_skmer, 0U);
    EXPECT_TRUE(skmerator.begin().consumed());
}


TEST(Skmerator, seq_length_exactly_k)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"ACGTA"};  // length = k
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    uint64_t nb_skmer {0};
    for ([[maybe_unused]] km::Skmer<kuint> skmer : skmerator)
        nb_skmer += 1;

    EXPECT_EQ(nb_skmer, 1U);
}


// ------------------------- HOMOPOLYMER TESTS -------------------------

// Homopolymers: the minimizer (AA / its RC TT) is no longer raw 0 under φ, and the raw
// pref/suff shape is an implementation detail. Assert the coverage invariant instead.
TEST(Skmerator, all_A_homopolymer) { expect_covers_all_kmers<uint16_t>(5, 2, "AAAAAA"); }
TEST(Skmerator, all_T_homopolymer) { expect_covers_all_kmers<uint16_t>(5, 2, "TTTTTT"); }


TEST(Skmerator, all_A_equals_all_T_skmers)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip_a {k, m};
    std::string seq_a{"AAAAAA"};
    km::SeqSkmerator<kuint> skmerator_a {manip_a, seq_a};
    std::vector<km::Skmer<kuint> > skmers_a;
    for (km::Skmer<kuint> skmer : skmerator_a)
        skmers_a.emplace_back(skmer);

    km::SkmerManipulator<kuint> manip_t {k, m};
    std::string seq_t{"TTTTTT"};
    km::SeqSkmerator<kuint> skmerator_t {manip_t, seq_t};
    std::vector<km::Skmer<kuint> > skmers_t;
    for (km::Skmer<kuint> skmer : skmerator_t)
        skmers_t.emplace_back(skmer);

    // AAAAAA and TTTTTT are reverse-complements of each other, and both map to
    // the same canonical form at every sliding window. Their yielded skmer
    // streams must therefore be pairwise equal.
    ASSERT_EQ(skmers_a.size(), skmers_t.size());
    for (uint64_t i{0} ; i<skmers_a.size() ; i++)
        ASSERT_TRUE(manip_a.skmer_equals(skmers_a[i], skmers_t[i]));
}


// ------------------------- ITERATOR STATE TESTS -------------------------

TEST(Skmerator, iterator_idempotent_after_end)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"ATCGACTGTGTACACT"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    auto it = skmerator.begin();
    auto end = skmerator.end();

    while (not (it == end))
        ++it;

    EXPECT_TRUE(it.consumed());
    EXPECT_TRUE(it == end);

    // Extra increments on an already consumed iterator must be no-ops.
    for (int i{0} ; i<5 ; i++)
    {
        ++it;
        EXPECT_TRUE(it.consumed());
        EXPECT_TRUE(it == end);
    }
}


// ------------------------- FILE ITERATOR TESTS -------------------------

TEST(Skmerator, file_multiple_sequences)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    // fasta0.fa contains two records:
    //   >seq1  ATCGACTGTGTACACT
    //   >seq2  ACAACACTTACTAGGA
    std::string seq1{"ATCGACTGTGTACACT"};
    std::string seq2{"ACAACACTTACTAGGA"};

    km::SkmerManipulator<kuint> seq_manip {k, m};
    std::vector<km::Skmer<kuint> > expected;

    km::SeqSkmerator<kuint> seq1_skmerator {seq_manip, seq1};
    for (km::Skmer<kuint> skmer : seq1_skmerator)
        expected.emplace_back(skmer);

    km::SeqSkmerator<kuint> seq2_skmerator {seq_manip, seq2};
    for (km::Skmer<kuint> skmer : seq2_skmerator)
        expected.emplace_back(skmer);

    std::string filename{"../tests/data/fasta0.fa"};
    km::SkmerManipulator<kuint> file_manip {k, m};
    km::FileSkmerator<kuint> file_skmerator {file_manip, filename};

    std::vector<km::Skmer<kuint> > file_skmers;
    for (km::Skmer<kuint> skmer : file_skmerator)
        file_skmers.emplace_back(skmer);

    ASSERT_EQ(expected.size(), file_skmers.size());
    for (uint64_t i{0} ; i<expected.size() ; i++)
        ASSERT_TRUE(file_manip.skmer_equals(expected[i], file_skmers[i]));
}


TEST(Skmerator, file_skips_short_sequences)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    // short_skip.fa interleaves records whose length is below k with records
    // long enough to produce skmers. The file iterator must silently skip the
    // short ones and concatenate the skmer streams of the remaining records.
    std::string long1{"ATCGACTGTGTACACT"};
    std::string long2{"CCCCAAAAA"};

    km::SkmerManipulator<kuint> seq_manip {k, m};
    std::vector<km::Skmer<kuint> > expected;

    km::SeqSkmerator<kuint> long1_skmerator {seq_manip, long1};
    for (km::Skmer<kuint> skmer : long1_skmerator)
        expected.emplace_back(skmer);

    km::SeqSkmerator<kuint> long2_skmerator {seq_manip, long2};
    for (km::Skmer<kuint> skmer : long2_skmerator)
        expected.emplace_back(skmer);

    std::string filename{"../tests/data/short_skip.fa"};
    km::SkmerManipulator<kuint> file_manip {k, m};
    km::FileSkmerator<kuint> file_skmerator {file_manip, filename};

    std::vector<km::Skmer<kuint> > file_skmers;
    for (km::Skmer<kuint> skmer : file_skmerator)
        file_skmers.emplace_back(skmer);

    ASSERT_EQ(expected.size(), file_skmers.size());
    for (uint64_t i{0} ; i<expected.size() ; i++)
        ASSERT_TRUE(file_manip.skmer_equals(expected[i], file_skmers[i]));
}