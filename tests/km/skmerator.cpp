#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

using namespace std;


TEST(Skmerator, update_equal_mini_fwd_fwd)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"CCCCCC"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};
    km::SkmerPrettyPrinter<kuint> pp {k, m};


    //                         Prefix:         C   C   _   _             C   C   _   _
    //                         Suffix:       C   C   C   _             C   C   C   _
    const kuint expected_values[][2] { {0, 0b0101010101111111U}, {0, 0b0101010101111111U} };
    const uint64_t expected_prefixes[2] {1, 1};
    const uint64_t expected_suffixes[2] {2, 2};

    uint64_t nb_skmer {0};
    for ([[maybe_unused]]km::Skmer<kuint> skmer : skmerator)
    {
        pp << skmer;

        //                            Less significant             Most significant
        const kpair expected_pair{expected_values[nb_skmer][1], expected_values[nb_skmer][0]};
        ASSERT_EQ(expected_pair, skmer.m_pair);
        ASSERT_EQ(expected_prefixes[nb_skmer], skmer.m_pref_size);
        ASSERT_EQ(expected_suffixes[nb_skmer], skmer.m_suff_size);

        nb_skmer += 1;
    }

    EXPECT_EQ(nb_skmer, 2);
}


TEST(Skmerator, update_equal_mini_rev_rev)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"GGGGGG"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};
    km::SkmerPrettyPrinter<kuint> pp {k, m};


    //                         Prefix:         C   C   _   _             C   C   _   _
    //                         Suffix:       C   C   C   _             C   C   C   _
    const kuint expected_values[][2] { {0, 0b0101010101111111U}, {0, 0b0101010101111111U} };
    const uint64_t expected_prefixes[2] {1, 1};
    const uint64_t expected_suffixes[2] {2, 2};

    uint64_t nb_skmer {0};
    for ([[maybe_unused]]km::Skmer<kuint> skmer : skmerator)
    {
        pp << skmer;

        //                            Less significant             Most significant
        const kpair expected_pair{expected_values[nb_skmer][1], expected_values[nb_skmer][0]};
        ASSERT_EQ(expected_pair, skmer.m_pair);
        ASSERT_EQ(expected_prefixes[nb_skmer], skmer.m_pref_size);
        ASSERT_EQ(expected_suffixes[nb_skmer], skmer.m_suff_size);

        nb_skmer += 1;
    }

    EXPECT_EQ(nb_skmer, 2);
}

TEST(Skmerator, decreasing_minimizer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"CCCCAAAAA"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};
    km::SkmerPrettyPrinter<kuint> pp {k, m};


    //                         Prefix:         C   C   C   C             A   C   C   C             A   A   C   _
    //                         Suffix:       A   _   _   _             A   A   _   _             A   A   A   _
    const kuint expected_values[][2] { {0, 0b0001110111011101U}, {0, 0b0000000111011101U}, {0, 0b0000000000011111U}
    };

    uint64_t nb_skmer {0};
    for ([[maybe_unused]]km::Skmer<kuint> skmer : skmerator)
    {
        // pp << skmer;
        // cout << pp << endl;

        //                            Less significant             Most significant
        const kpair expected_pair{expected_values[nb_skmer][1], expected_values[nb_skmer][0]};
        ASSERT_EQ(expected_pair, skmer.m_pair);

        nb_skmer += 1;
    }
    EXPECT_EQ(nb_skmer, 3);
}


TEST(Skmerator, increasing_minimizer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"AACCCC"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};
    km::SkmerPrettyPrinter<kuint> pp {k, m};


    //                         Prefix:         A   _   _   _             A   _   _   _
    //                         Suffix:       A   C   C   C             C   C   C   C
    const kuint expected_values[][2] { {0, 0b0000011101110111U}, {0, 0b0100011101110111U}
    };

    uint64_t nb_skmer {0};
    for ([[maybe_unused]]km::Skmer<kuint> skmer : skmerator)
    {
        ASSERT_TRUE(nb_skmer < 2);

        // pp << skmer;
        // cout << pp << endl;

        //                            Less significant             Most significant
        const kpair expected_pair{expected_values[nb_skmer][1], expected_values[nb_skmer][0]};
        ASSERT_EQ(expected_pair, skmer.m_pair);

        nb_skmer += 1;
    }

    EXPECT_EQ(nb_skmer, 2);
}


TEST(Skmerator, outofcontext_minimizer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    const uint64_t k{8};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"AACAATAAGGGGGGG"};
    // cout << seq << endl;
    km::SeqSkmerator<kuint> skmerator {manip, seq};
    km::SkmerPrettyPrinter<kuint> pp {k, m};


    //                         Prefix:      A   _   _        _   _   _   _          A   C   A       _   _   _   _
    //                         Suffix:    A   C   A        A   T   A   A          A   T   A       A   G   G   G
    const kuint expected_values[][2] { {0b000001110011U, 0b0011101100110011U}, {0b000010010000U, 0b0011111111111111U}
        // Prefix:     A   T   A        _   _   _   _          C   C   C       _   _   _   _
        // Suffix:   A   G   G        G   G   G   G          C   C   C       C   T   _   _
        ,         {0b000011101100U, 0b1111111111111111U}, {0b010101010101, 0b0111101111111111U}
    };

    uint64_t nb_skmer {0};
    for ([[maybe_unused]]km::Skmer<kuint> skmer : skmerator)
    {
        ASSERT_TRUE(nb_skmer < 4);

        // pp << skmer;
        // cout << pp << endl;

        //                            Less significant             Most significant
        const kpair expected_pair{expected_values[nb_skmer][1], expected_values[nb_skmer][0]};
        ASSERT_EQ(expected_pair, skmer.m_pair);

        nb_skmer += 1;
    }

    ASSERT_EQ(nb_skmer, 4);
}


TEST(Skmerator, seq_test_5_2)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};
    km::SkmerPrettyPrinter<kuint> pp {k, m};

    // --- Sequence ---
    std::string seq{"ATCGACTGTGTACACT"};
    km::SkmerManipulator<kuint> seq_manip {k, m};
    uint64_t const expected_skmers[] = {1,2,2,3,3,4,4,4,4,5,5,5};
    
    for (uint seq_size{k} ; seq_size<=seq.length() ; seq_size++) {
        // cout << "seq_size: " << seq_size << " ---" << endl;
        std::string sub{seq.substr(0, seq_size)};
        // cout << sub << endl;
        km::SeqSkmerator<kuint> seq_skmerator {seq_manip, sub};

        // Enumerates the superkmers from the sequence
        std::vector<km::Skmer<kuint> > seq_skmers {};
        for (km::Skmer<kuint> const skmer : seq_skmerator) {
            // pp << skmer;
            // cout << pp << endl;
            seq_skmers.emplace_back(skmer);
        }

        ASSERT_EQ(seq_skmers.size(), expected_skmers[seq_size-k]);
    }
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

TEST(Skmerator, all_A_homopolymer)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"AAAAAA"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    uint64_t nb_skmer {0};
    for (km::Skmer<kuint> skmer : skmerator)
    {
        // Minimizer of an all-A window is AA = 0b0000 = 0 (forward canonical).
        ASSERT_EQ(manip.minimizer(skmer), static_cast<kuint>(0));
        // Same skmer shape as the CCCCCC / GGGGGG reference tests.
        ASSERT_EQ(skmer.m_pref_size, 1U);
        ASSERT_EQ(skmer.m_suff_size, 2U);

        nb_skmer += 1;
    }

    EXPECT_EQ(nb_skmer, 2U);
}


TEST(Skmerator, all_T_homopolymer)
{
    using kuint = uint16_t;

    const uint64_t k{5};
    const uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    std::string seq{"TTTTTT"};
    km::SeqSkmerator<kuint> skmerator {manip, seq};

    uint64_t nb_skmer {0};
    for (km::Skmer<kuint> skmer : skmerator)
    {
        // Reverse-complement of TT is AA => canonical minimizer is 0.
        ASSERT_EQ(manip.minimizer(skmer), static_cast<kuint>(0));
        ASSERT_EQ(skmer.m_pref_size, 1U);
        ASSERT_EQ(skmer.m_suff_size, 2U);

        nb_skmer += 1;
    }

    EXPECT_EQ(nb_skmer, 2U);
}


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