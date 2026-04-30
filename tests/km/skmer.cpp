#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <sstream>

#include <io/Skmer.hpp>

using namespace std;


TEST(SkmerManipulator, enumerate_kmer_full)
{
	using kuint = uint8_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{1};

	const string seq {"TCAAGCACTAA"};
    //                            prefix:        _   _       _   _  ,     _   _       _   _  ,     _   _       _   _  ,
    //                            suffix:          _       _   T    ,       _       T   C    ,       T       C   A    ,
    const kuint expected_fwd_values[][2] = {{0b00000000, 0b00001000},{0b00000000, 0b10000100},{0b00001000, 0b01000000},
    //                                           T   _       _   _  ,     C   T       _   _  ,     A   C       T   _  ,
    //                                             C       A   A    ,       A       A   G    ,       A       G   C    ,
                                            {0b00100100, 0b00000000},{0b00010010, 0b00001100},{0b00000001, 0b11100100},
    //                                           A   A       C   T  ,     G   A       A   C  ,     C   G       A   A  ,
    //                                             G       C   A    ,       C       A   C    ,       A       C   T    ,
                                            {0b00001100, 0b01010010},{0b00110100, 0b00000101},{0b00010011, 0b01001000},
    //                                           A   C       G   A  ,     C   A       C   G
    //                                             C       T   A    ,       T       A   A
                                            {0b00000101, 0b10110000},{0b00011000, 0b00010011}};

    //                            prefix:        _   _       _   A  ,     _   _       A   G  ,     _   A       G   T  ,
    //                            suffix:          _       _   _    ,       _       _   _    ,       _       _   _    ,
    const kuint expected_rev_values[][2] = {{0b00000000, 0b00000000},{0b00000000, 0b00000011},{0b00000000, 0b00110010},
    //                                           A   G       T   T  ,     G   T       T   C  ,     T   T       C   G  ,
    //                                             _       _   _    ,       A       _   _    ,       G       A   _    ,
                                            {0b00000011, 0b00100010},{0b00110010, 0b00100001},{0b00101110, 0b00010011},
    //                                           T   C       G   T  ,     C   G       T   G  ,     G   T       G   A  ,
    //                                             T       G   A    ,       T       T   G    ,       C       T   T    ,
                                            {0b00101001, 0b11110010},{0b00011011, 0b10101111},{0b00110110, 0b10111000},
    //                                           T   G       A   T  ,     G   A       T   T
    //                                             G       C   T    ,       T       G   C
                                            {0b00101111, 0b01001010},{0b00111000, 0b11100110}};

    km::SkmerManipulator<kuint> manip {k, m};

    uint64_t idx{0};
    for (const auto& letter : seq)
    {
    	kuint nucl {static_cast<kuint>((letter >> 1) & 0b11)};
        [[maybe_unused]]const auto& min_skmer {manip.add_nucleotide(nucl)};

        //                            Less significant             Most significant
        const kpair expected_fwd_pair{expected_fwd_values[idx][1], expected_fwd_values[idx][0]};
        ASSERT_EQ(manip.m_fwd.m_pair, expected_fwd_pair);

        const kpair expected_rev_pair{expected_rev_values[idx][1], expected_rev_values[idx][0]};
        ASSERT_EQ(manip.m_rev.m_pair, expected_rev_pair);

        // cout << "PAIR " << km::Skmer<kuint>::pair(expected_fwd_values[idx]) << endl;
        EXPECT_EQ(min_skmer.m_pair, std::min(
                    km::Skmer<kuint>::pair(expected_fwd_values[idx][1], expected_fwd_values[idx][0]),
                    km::Skmer<kuint>::pair(expected_rev_values[idx][1], expected_rev_values[idx][0])));

        idx += 1;
    }
}


// TEST(SkmerManipulator, enumerate_kmer_partial)
// {
//     using kuint = uint8_t;
//     using kpair = km::Skmer<kuint>::pair;
//     constexpr uint64_t k{4};
//     constexpr uint64_t m{2};

//     const string seq {"TCAAGCACTAA"};
//     //                            prefix:            _       _   _  ,         _       _   _  ,         _       _   _  ,
//     //                            suffix:          _       _   T    ,       _       T   C    ,       T       C   A    ,
//     const kuint expected_fwd_values[][2] = {{0b00000000, 0b00001000},{0b00000000, 0b10000100},{0b00001000, 0b01000000},
//     //                                               T       _   _  ,         C       T   _  ,         A       C   T  ,
//     //                                             C       A   A    ,       A       A   G    ,       A       G   C    ,
//                                             {0b00000110, 0b00000000},{0b00000001, 0b00101100},{0b00000000, 0b11010110},
//     //                                               A       A   C  ,         G       A   A  ,         C       G   A  ,
//     //                                             G       C   A    ,       C       A   C    ,       A       C   T    ,
//                                             {0b00001100, 0b01000001},{0b00000111, 0b00000100},{0b00000001, 0b01111000},
//     //                                               A       C   G  ,         C       A   C
//     //                                             C       T   A    ,       T       A   A
//                                             {0b00000100, 0b10010011},{0b00010001, 0b00000001}};

//     // //                            prefix:        _   _       _   A  ,     _   _       A   G  ,     _   A       G   T  ,
//     // //                            suffix:          _       _   _    ,       _       _   _    ,       _       _   _    ,
//     // const kuint expected_rev_values[][2] = {{0b00000000, 0b00000000},{0b00000000, 0b00000011},{0b00000000, 0b00110010},
//     // //                                           A   G       T   T  ,     G   T       T   C  ,     T   T       C   G  ,
//     // //                                             _       _   _    ,       A       _   _    ,       G       A   _    ,
//     //                                         {0b00000011, 0b00100010},{0b00110010, 0b00100001},{0b00101110, 0b00010011},
//     // //                                           T   C       G   T  ,     C   G       T   G  ,     G   T       G   A  ,
//     // //                                             T       G   A    ,       T       T   G    ,       C       T   T    ,
//     //                                         {0b00101001, 0b11110010},{0b00011011, 0b10101111},{0b00110110, 0b10111000},
//     // //                                           T   G       A   T  ,     G   A       T   T
//     // //                                             G       C   T    ,       T       G   C
//     //                                         {0b00101111, 0b01001010},{0b00111000, 0b11100110}};

//     km::SkmerManipulator<kuint> manip {k, m};

//     uint64_t idx{0};
//     for (const auto& letter : seq)
//     {
//         kuint nucl {static_cast<kuint>((letter >> 1) & 0b11)};
//         [[maybe_unused]]const auto& min_skmer {manip.add_nucleotide(nucl)};

//         //                            Less significant             Most significant
//         const kpair expected_fwd_pair{expected_fwd_values[idx][1], expected_fwd_values[idx][0]};
//         ASSERT_EQ(manip.m_fwd.m_pair, expected_fwd_pair);

//         // const kpair expected_rev_pair{expected_rev_values[idx][1], expected_rev_values[idx][0]};
//         // ASSERT_EQ(manip.m_rev.m_pair, expected_rev_pair);

//         // EXPECT_EQ(min_skmer.m_pair, std::min(km::Skmer<kuint>::pair(expected_fwd_values[idx]), km::Skmer<kuint>::pair(expected_rev_values[idx])));

//         idx += 1;
//     }
// }


TEST(SkmerManipulator, output)
{
    using kuint = uint16_t;
    // using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    const string seq {"TCAAGCA"};
    km::SkmerManipulator<kuint> manip {k, m};
    km::Skmer<kuint> min_skmer {};

    for (const auto& letter : seq)
    {
        kuint nucl {static_cast<kuint>((letter >> 1) & 0b11)};
        min_skmer = manip.add_nucleotide(nucl);
    }

    stringstream ss {};
    ss << manip;

    string output = ss.str();

    ASSERT_EQ(output, "[not interleaved: ATCA AGCA / TGCT TGAA]");
}

TEST(SkmerManipulator, minimizer_extraction)
{
    using kuint = uint16_t;
    // using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    const string seq {"TCAAGCATGTTAG"};
    km::SkmerManipulator<kuint> manip {k, m};

    for (uint64_t i {0}; i < k-m+1; i++)
    {
        kuint nucl {static_cast<kuint>((seq[i] >> 1) & 0b11)};
        manip.add_nucleotide(nucl);
    }

    kuint expected_minimizer[] {0b0011, 0b0001, 0b0000, 0b1001, 0b0111, 0b0001, 0b1000, 0b0001, 0b0100};

    for (uint64_t i {k-m+1}; i < seq.length(); i++)
    {
        kuint nucl {static_cast<kuint>((seq[i] >> 1) & 0b11)};
        manip.add_nucleotide(nucl);

        //cout << manip.minimizer() << endl;
        ASSERT_EQ(manip.minimizer(),expected_minimizer[i-k+m-1]);
    }
}


//     //                            prefix:            _       _   _  ,         _       _   _  ,         _       _   _  ,
//     //                            suffix:          _       _   T    ,       _       T   C    ,       T       C   A    ,

TEST(SkmerManipulator, overlap_masks_generation)
    {
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};

    std::array < kpair, 5 > const expected_kpairs { kpair(0b0011001100110011U,0), kpair(0b1111001100110000U,0),
                                              kpair(0b1111111100000000U,0), kpair(0b1111110011000000U,0),
                                              kpair(0b1100110011001100U,0) };
    // auto const masks {manip.generate_masks_sp()};
    auto const masks {manip.get_sp_mask()};

    for(uint64_t pos {0}; pos < expected_kpairs.size(); pos += 1){
        ASSERT_EQ(expected_kpairs[pos],masks[pos]);
    }

}

TEST(SkmerManipulator, generate_masks_k)
    {
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};

    std::array < kpair, 4 > const expected_kpairs { kpair(0b1111001100110011U,0), kpair(0b1111111100110000U,0),
                                              kpair(0b1111111111000000U,0), kpair(0b1111110011001100U,0) };
    auto const masks {manip.get_k_mask()};
    // auto const masks {manip.kmer_masks};
    for(uint64_t pos {0}; pos < expected_kpairs.size(); pos += 1){
        ASSERT_EQ(expected_kpairs[pos],masks[pos]);
    }

}

TEST(SkmerManipulator, generate_masks_nucleotide)
    {
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};

    std::array < kpair, 8 > const expected_kpairs { kpair(0b0000000000000011U,0), kpair(0b000000000110000U,0),
                                              kpair(0b0000001100000000U,0), kpair(0b0011000000000000U,0),
                                              kpair(0b1100000000000000U,0), kpair(0b0000110000000000U,0),
                                              kpair(0b0000000011000000U,0), kpair(0b0000000000001100U,0) };
    auto const masks {manip.get_n_mask()};

    for(uint64_t pos {0}; pos < expected_kpairs.size(); pos += 1){
        ASSERT_EQ(expected_kpairs[pos],masks[pos]);
    }

}

TEST(SkmerManipulator, extract_nucleotide_from_skmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};

    //                    Prefix:      A   C   G   T
    //                    Suffix:    A   C   T   C
    const kpair m_basic_pair {0b0000010110110110U, 0};
    const km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,4,4);

    std::array < kpair, 8 > const expected_kpairs { kpair(0b0000000000000010U,0), kpair(0b000000000110000U,0),
                                              kpair(0b0000000100000000U,0), kpair(0b0000000000000000U,0),
                                              kpair(0b0000000000000000U,0), kpair(0b0000010000000000U,0),
                                              kpair(0b0000000010000000U,0), kpair(0b0000000000000100U,0) };
    uint64_t position {0};
    for (kpair expected_nucleotide : expected_kpairs){
        ASSERT_EQ(expected_nucleotide,manip.extract_nucleotide(m_skmer, position));
        position++;
    }
}

TEST(SkmerManipulator, clean_nucleotide_at_position)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};

    //                    Prefix:      A   C   G   T
    //                    Suffix:    A   C   T   C
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,4,4);
    std::array < kpair, 8 > const expected_kpairs { kpair(0b0000010110110100U,0), kpair(0b0000010110000110U,0),
        kpair(0b0000010010110110U,0), kpair(0b0000010110110110U,0),
        kpair(0b0000010110110110U,0), kpair(0b0000000110110110U,0),
        kpair(0b0000010100110110U,0), kpair(0b0000010110110010U,0) };

    uint64_t position {0};
    for (kpair expected_nucleotide : expected_kpairs){
        manip.clean_nucleotide_position_skmer(m_skmer, position);
        ASSERT_EQ(expected_nucleotide,m_skmer.m_pair);
        position++;
        m_skmer = km::Skmer<kuint>(m_basic_pair,4,4);
    }
}


TEST(SkmerManipulator, HasValidSkmer1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,1,2);

    std::vector<bool> expected_results {false,false,true,false};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}
TEST(SkmerManipulator, HasValidSkmer2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,3,0);

    std::vector<bool> expected_results {true,false,false,false};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}
TEST(SkmerManipulator, HasValidSkmer3)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,2,1);

    std::vector<bool> expected_results {false,true,false,false};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}
TEST(SkmerManipulator, HasValidSkmer4)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,0,3);

    std::vector<bool> expected_results {false,false,false,true};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}

TEST(SkmerManipulator, HasValidSkmer5)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,3,3);

    std::vector<bool> expected_results {true,true,true,true};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}

TEST(SkmerManipulator, HasValidSkmer6)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,2,3);

    std::vector<bool> expected_results {false,true,true,true};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}

TEST(SkmerManipulator, HasValidSkmer7)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{3};

    km::SkmerManipulator<kuint> manip {k, m};
    const kpair m_basic_pair {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer = km::Skmer<kuint>(m_basic_pair,1,2);

    std::vector<bool> expected_results {false,true,true};
    for (size_t position {0}; position <= k-m; position++){
        ASSERT_EQ(expected_results[position],manip.has_valid_kmer(m_skmer,position));
    }
}

TEST(SkmerManipulator, KmerLessThanKmer1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:          A   T   T                C   C   C
    //                 Suffix:        A   _   _                A   _   _
    const kpair input_skmers[2]{{0b000011101110U, 0}, {0b000111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 1),
        km::Skmer<kuint>(input_skmers[1], 3, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}

TEST(SkmerManipulator, KmerLessThanKmer2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{3};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:       A   T             C   C
    //                 Suffix:     A   _             A   _
    const kpair input_skmers[2]{{0b00001110U, 0}, {0b00011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 2, 1),
        km::Skmer<kuint>(input_skmers[1], 2, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}

TEST(SkmerManipulator, KmerLessThanKmer3)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:       A   T   T   T             G   C   C   C
    //                 Suffix:     A   _   _   _             A   _   _   _
    const kpair input_skmers[2]{{0b0000111011101110U, 0}, {0b0011110111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 4, 1),
        km::Skmer<kuint>(input_skmers[1], 4, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}

TEST(SkmerManipulator, KmerLessThanKmer4)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:       A   T   T   _             G   C   C   _
    //                 Suffix:     A   T   _   _             A   A   _   _
    const kpair input_skmers[2]{{0b0000101011101111U, 0}, {0b0011000111011111U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 2),
        km::Skmer<kuint>(input_skmers[1], 3, 2)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}


TEST(SkmerManipulator, KmerLTKmer1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:          A   T   T                C   C   C
    //                 Suffix:        A   _   _                A   _   _
    const kpair input_skmers[2]{{0b000011101110U, 0}, {0b000111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 1),
        km::Skmer<kuint>(input_skmers[1], 3, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_lt_kmer(skmer_vector[0], position, skmer_vector[1], position));
}

TEST(SkmerManipulator, KmerLTKmer2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{3};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:       A   T             C   C
    //                 Suffix:     A   _             A   _
    const kpair input_skmers[2]{{0b00001110U, 0}, {0b00011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 2, 1),
        km::Skmer<kuint>(input_skmers[1], 2, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_lt_kmer(skmer_vector[0], position, skmer_vector[1], position));
}

TEST(SkmerManipulator, KmerLTKmer3)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:       A   T   T   T             G   C   C   C
    //                 Suffix:     A   _   _   _             A   _   _   _
    const kpair input_skmers[2]{{0b0000111011101110U, 0}, {0b0011110111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 4, 1),
        km::Skmer<kuint>(input_skmers[1], 4, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_lt_kmer(skmer_vector[0], position, skmer_vector[1], position));
}

TEST(SkmerManipulator, KmerLTKmer4)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:       A   T   T   _             G   C   C   _
    //                 Suffix:     A   T   _   _             A   A   _   _
    const kpair input_skmers[2]{{0b0000101011101111U, 0}, {0b0011000111011111U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 2),
        km::Skmer<kuint>(input_skmers[1], 3, 2)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_lt_kmer(skmer_vector[0], position, skmer_vector[1], position));
}


TEST(get_valid_kmer_bounds, getValidKmerBounds1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};

    const kpair input_skmer{0, 0};
    const std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmer, 2, 0),
        km::Skmer<kuint>(input_skmer, 2, 1),
        km::Skmer<kuint>(input_skmer, 2, 2),
        km::Skmer<kuint>(input_skmer, 1, 1),
        km::Skmer<kuint>(input_skmer, 1, 2),
        km::Skmer<kuint>(input_skmer, 0, 2),
    };

    const std::vector<std::pair<uint64_t, uint64_t>> expected_values{
        {0,0},
        {0,1},
        {0,2},
        {1,1},
        {1,2},
        {2,2},
    };

    for(size_t pos {0}; pos < expected_values.size(); pos++){
        ASSERT_EQ(expected_values[pos], manip.get_valid_kmer_bounds(skmer_vector[pos]));
    }
}

TEST(get_valid_kmer_bounds, getValidKmerBounds2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{3};

    km::SkmerManipulator<kuint> manip {k, m};

    const kpair input_skmer{0, 0};
    const std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmer, 2, 0),
        km::Skmer<kuint>(input_skmer, 2, 1),
        km::Skmer<kuint>(input_skmer, 2, 2),
        km::Skmer<kuint>(input_skmer, 1, 1),
        km::Skmer<kuint>(input_skmer, 1, 2),
        km::Skmer<kuint>(input_skmer, 0, 2),
    };

    const std::vector<std::pair<uint64_t, uint64_t>> expected_values{
        {0,0},
        {0,1},
        {0,2},
        {1,1},
        {1,2},
        {2,2},
    };

    for(size_t pos {0}; pos < expected_values.size(); pos++){
        ASSERT_EQ(expected_values[pos], manip.get_valid_kmer_bounds(skmer_vector[pos]));
    }
}


TEST(SkmerManipulator, KmerEqualsKmer1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:          A   T   T                A   T   T
    //                 Suffix:        A   _   _                A   _   _
    const kpair input_skmers[2]{{0b000011101110U, 0}, {0b000111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 1),
        km::Skmer<kuint>(input_skmers[1], 3, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}

TEST(SkmerManipulator, KmerEqualsKmer2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:          A   T   T                A   C   T
    //                 Suffix:        A   _   _                A   _   _
    const kpair input_skmers[2]{{0b000011101110U, 0}, {0b000111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 1),
        km::Skmer<kuint>(input_skmers[1], 3, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}


TEST(SkmerManipulator, KmerEqualsKmer3)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                 Prefix:          A   T   T                A   C   T
    //                 Suffix:        A   _   _                A   _   _
    const kpair input_skmers[2]{{0b000011101110U, 0}, {0b000111011101U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{
        km::Skmer<kuint>(input_skmers[0], 3, 1),
        km::Skmer<kuint>(input_skmers[1], 3, 1)
    };
    const uint64_t position{0};

    bool expected_result {true};

    ASSERT_EQ(expected_result, manip.kmer_compare(skmer_vector[0], skmer_vector[1], position) < 0);
}


// -------------------------------------------------------------------
//                     Skmer::pair bitwise operations
// -------------------------------------------------------------------

TEST(SkmerPair, ShiftRightWithinKuint)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    // Bits 16..47 all set; shift right 16 should land them at 0..31.
    kpair p{0xFFFF0000u, 0x0000FFFFu};
    p >>= 16;
    EXPECT_EQ(p.m_value[0], 0xFFFFFFFFu);
    EXPECT_EQ(p.m_value[1], 0x00000000u);
}

TEST(SkmerPair, ShiftRightAtKuintBoundary)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xDEADBEEFu, 0xCAFEBABEu};
    p >>= 32;
    EXPECT_EQ(p.m_value[0], 0xCAFEBABEu);
    EXPECT_EQ(p.m_value[1], 0x00000000u);
}

TEST(SkmerPair, ShiftRightBeyondKuintBoundary)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xDEADBEEFu, 0xCAFEBABEu};
    p >>= 40;
    EXPECT_EQ(p.m_value[0], 0x00CAFEBAu);
    EXPECT_EQ(p.m_value[1], 0x00000000u);
}

TEST(SkmerPair, ShiftRightByZeroIsIdentity)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xDEADBEEFu, 0xCAFEBABEu};
    p >>= 0;
    EXPECT_EQ(p.m_value[0], 0xDEADBEEFu);
    EXPECT_EQ(p.m_value[1], 0xCAFEBABEu);
}

TEST(SkmerPair, ShiftLeftWithinKuint)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    // Top nibble of m_value[0] should cross into m_value[1] after a 4-bit left shift.
    kpair p{0xA0000000u, 0x00000000u};
    p <<= 4;
    EXPECT_EQ(p.m_value[0], 0x00000000u);
    EXPECT_EQ(p.m_value[1], 0x0000000Au);
}

TEST(SkmerPair, ShiftLeftAtKuintBoundary)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xDEADBEEFu, 0xCAFEBABEu};
    p <<= 32;
    EXPECT_EQ(p.m_value[0], 0x00000000u);
    EXPECT_EQ(p.m_value[1], 0xDEADBEEFu);
}

TEST(SkmerPair, ShiftLeftBeyondKuintBoundary)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x000000DEu, 0x00000000u};
    p <<= 40;
    EXPECT_EQ(p.m_value[0], 0x00000000u);
    EXPECT_EQ(p.m_value[1], 0x0000DE00u);
}

TEST(SkmerPair, ShiftLeftByZeroIsIdentity)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xDEADBEEFu, 0xCAFEBABEu};
    p <<= 0;
    EXPECT_EQ(p.m_value[0], 0xDEADBEEFu);
    EXPECT_EQ(p.m_value[1], 0xCAFEBABEu);
}

TEST(SkmerPair, BitwiseNot)
{
    using kuint = uint8_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{static_cast<kuint>(0b10101010u), static_cast<kuint>(0b01010101u)};
    kpair q = ~p;
    EXPECT_EQ(q.m_value[0], static_cast<kuint>(0b01010101u));
    EXPECT_EQ(q.m_value[1], static_cast<kuint>(0b10101010u));
}

TEST(SkmerPair, AndWithUint64WipesHighHalf)
{
    // Pin documented behavior: operator& (uint64_t) forces m_value[1] to 0.
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xFFFFFFFFu, 0xFFFFFFFFu};
    kpair q = p & static_cast<uint64_t>(0xFFu);
    EXPECT_EQ(q.m_value[0], 0xFFu);
    EXPECT_EQ(q.m_value[1], 0x0u);
}

TEST(SkmerPair, AndWithPairKeepsBothHalves)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xFFFFFFFFu, 0xFFFFFFFFu};
    kpair mask{0x000000FFu, 0x0000FF00u};
    kpair r = p & mask;
    EXPECT_EQ(r.m_value[0], 0x000000FFu);
    EXPECT_EQ(r.m_value[1], 0x0000FF00u);
}

TEST(SkmerPair, OrAssignKuintOnlyTouchesLowHalf)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x0u, 0x0u};
    kuint v = 0xDEADBEEFu;
    p |= v;
    EXPECT_EQ(p.m_value[0], 0xDEADBEEFu);
    EXPECT_EQ(p.m_value[1], 0x0u);
}

TEST(SkmerPair, ImplicitConversionToUint64DropsHighHalf)
{
    using kuint = uint32_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0xDEADBEEFu, 0xCAFEBABEu};
    uint64_t v = static_cast<uint64_t>(p);
    EXPECT_EQ(v, static_cast<uint64_t>(0xDEADBEEFu));
}

TEST(SkmerPair, LessEqualSameHighPath)
{
    using kuint = uint8_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair a{static_cast<kuint>(0x10u), static_cast<kuint>(0x05u)};
    kpair b{static_cast<kuint>(0x10u), static_cast<kuint>(0x05u)};
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(b <= a);
}

TEST(SkmerPair, LessEqualDifferentHighPath)
{
    using kuint = uint8_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair a{static_cast<kuint>(0xFFu), static_cast<kuint>(0x03u)};
    kpair b{static_cast<kuint>(0x00u), static_cast<kuint>(0x05u)};
    EXPECT_TRUE(a <= b);
    EXPECT_FALSE(b <= a);
}

TEST(SkmerPair, HasherIsDeterministic)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    km::Skmer<kuint>::pair_hasher h{};
    kpair p{0x1234u, 0x5678u};
    kpair q{0x1234u, 0x5678u};
    EXPECT_EQ(h(p), h(q));
}


// -------------------------------------------------------------------
//                            Skmer class
// -------------------------------------------------------------------

TEST(Skmer, CopyConstructorPreservesSizes)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x1234u, 0};
    km::Skmer<kuint> original(p, 3, 5);
    km::Skmer<kuint> copy(original);
    EXPECT_EQ(copy.m_pair, p);
    EXPECT_EQ(copy.m_pref_size, 3);
    EXPECT_EQ(copy.m_suff_size, 5);
}

TEST(Skmer, MoveConstructorPreservesSizes)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x1234u, 0};
    km::Skmer<kuint> original(p, 3, 5);
    km::Skmer<kuint> moved(std::move(original));
    EXPECT_EQ(moved.m_pair, p);
    EXPECT_EQ(moved.m_pref_size, 3);
    EXPECT_EQ(moved.m_suff_size, 5);
}

TEST(Skmer, CopyAssignmentPreservesSizes)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x1234u, 0};
    km::Skmer<kuint> original(p, 3, 5);
    km::Skmer<kuint> assigned{};
    assigned = original;
    EXPECT_EQ(assigned.m_pair, p);
    EXPECT_EQ(assigned.m_pref_size, 3);
    EXPECT_EQ(assigned.m_suff_size, 5);
}

TEST(Skmer, EqualityRequiresBothSizesMatch)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x1234u, 0};
    km::Skmer<kuint> a(p, 3, 5);
    km::Skmer<kuint> b(p, 3, 5);
    km::Skmer<kuint> c(p, 4, 5);
    km::Skmer<kuint> d(p, 3, 4);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
    EXPECT_FALSE(a == d);
}

TEST(Skmer, LessThanIgnoresSizes)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    kpair p{0x1234u, 0};
    km::Skmer<kuint> a(p, 3, 5);
    km::Skmer<kuint> b(p, 0, 0);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(b < a);
}

TEST(Skmer, StreamOutput)
{
    using kuint = uint8_t;
    using kpair = km::Skmer<kuint>::pair;
    km::Skmer<kuint> s{kpair{static_cast<kuint>(0b00001111u), static_cast<kuint>(0b11110000u)}, 2, 3};
    stringstream ss;
    ss << s;
    EXPECT_EQ(ss.str(), "11110000 00001111 pref:2 suff:3");
}


// -------------------------------------------------------------------
//                       SkmerPrettyPrinter
// -------------------------------------------------------------------

TEST(SkmerPrettyPrinter, PrintsInterleavedSkmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};
    // Prefix slots 0..3 = A,C,T,G; suffix slots 0..3 = A,C,T,G.
    // Bits [1:0]=00 [3:2]=00 [5:4]=01 [7:6]=01 [9:8]=10 [11:10]=10 [13:12]=11 [15:14]=11
    const kpair value{static_cast<kuint>(0xFA50u), 0};
    km::Skmer<kuint> s{value, 3, 3};

    km::SkmerPrettyPrinter<kuint> pp{k, m};
    pp << s;
    stringstream ss;
    ss << pp;
    EXPECT_EQ(ss.str(), "[skmer not interleaved: ACTG GTCA]");
}


// -------------------------------------------------------------------
//                  SkmerManipulator untested methods
// -------------------------------------------------------------------

TEST(SkmerManipulator, AddEmptyNucleotideSaturatesBothStrands)
{
    // Pushing the 0b11 sentinel sk_size times fills every slot with 0b11 on both
    // strands; add_empty_nucleotide does NOT complement (the sentinel has no complement).
    using kuint = uint16_t;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};
    constexpr uint64_t sk_size = 2 * k - m;

    km::SkmerManipulator<kuint> manip{k, m};
    for (uint64_t i = 0; i < sk_size; i++)
        manip.add_empty_nucleotide();

    EXPECT_EQ(manip.m_fwd.m_pair.m_value[0], static_cast<kuint>(0xFFFFu));
    EXPECT_EQ(manip.m_rev.m_pair.m_value[0], static_cast<kuint>(0xFFFFu));
}

TEST(SkmerManipulator, AddEmptyNucleotideDiffersFromAddNucleotideOnReverse)
{
    // add_nucleotide(0b11) pushes 0b11 forward but its complement (0b01) on reverse.
    // add_empty_nucleotide pushes 0b11 on both strands.
    using kuint = uint16_t;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> a{k, m};
    km::SkmerManipulator<kuint> b{k, m};
    a.add_nucleotide(static_cast<kuint>(0b11));
    b.add_empty_nucleotide();

    EXPECT_EQ(a.m_fwd.m_pair, b.m_fwd.m_pair);
    EXPECT_NE(a.m_rev.m_pair, b.m_rev.m_pair);
}

TEST(SkmerManipulator, MinimizerOverloadAndNoArgConsistency)
{
    using kuint = uint16_t;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    const string seq{"TCAAGCA"};
    for (const auto& letter : seq)
    {
        kuint nucl{static_cast<kuint>((letter >> 1) & 0b11)};
        manip.add_nucleotide(nucl);
    }

    const kuint fwd_minimizer = static_cast<kuint>(manip.m_fwd.m_pair >> (4 * (k - m)));
    const kuint rev_minimizer = static_cast<kuint>(manip.m_rev.m_pair >> (4 * (k - m)));

    EXPECT_EQ(manip.minimizer(manip.m_fwd), fwd_minimizer);
    EXPECT_EQ(manip.minimizer(manip.m_rev), rev_minimizer);
    EXPECT_EQ(manip.minimizer(), std::min(fwd_minimizer, rev_minimizer));
}

TEST(SkmerManipulator, MaskAbsentNucleotidesFillsTrimmedSlots)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};

    // pref_size=1, suff_size=1 → absent prefix slots {0,1} and absent suffix slots {0,1}.
    // Prefix slot i sits at bits [4i+1:4i]; setting to 0b11 yields 0x03 at i=0 and 0x30 at i=1.
    // Suffix slot i sits at bits [4i+3:4i+2]; setting to 0b11 yields 0x0C at i=0 and 0xC0 at i=1.
    // Combined low byte: 0xFF.
    km::Skmer<kuint> s{kpair{0u, 0u}, 1, 1};
    manip.mask_absent_nucleotides(s);
    EXPECT_EQ(s.m_pair.m_value[0], static_cast<kuint>(0x00FFu));
    EXPECT_EQ(s.m_pair.m_value[1], static_cast<kuint>(0u));
}

TEST(SkmerManipulator, GetSkmerOfKmerPreservesKmerBitsAndSizes)
{
    using kuint = uint16_t;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    const string seq{"TCAAGCAT"};
    for (const auto& letter : seq)
    {
        kuint nucl{static_cast<kuint>((letter >> 1) & 0b11)};
        manip.add_nucleotide(nucl);
    }

    const auto k_masks = manip.get_k_mask();
    const km::Skmer<kuint> full = manip.m_fwd;

    for (uint64_t pos{0}; pos <= k - m; pos++)
    {
        const auto extracted = manip.get_skmer_of_kmer(full, pos);
        EXPECT_EQ(extracted.m_pair & k_masks[pos], full.m_pair & k_masks[pos])
            << "kmer_pos=" << pos;
        const uint16_t expected_pref = static_cast<uint16_t>(((2*k - m + 1)/2) - (m + 1)/2 - pos);
        const uint16_t expected_suff = static_cast<uint16_t>(pos);
        EXPECT_EQ(extracted.m_pref_size, expected_pref) << "kmer_pos=" << pos;
        EXPECT_EQ(extracted.m_suff_size, expected_suff) << "kmer_pos=" << pos;
    }
}

TEST(SkmerManipulator, ExtractPrefixSuffixCoversAllPositions)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    const auto sp_masks = manip.get_sp_mask();
    ASSERT_EQ(sp_masks.size(), k - m + 2u);

    km::Skmer<kuint> s{kpair{0xFFFFu, 0u}};
    for (uint64_t pos{0}; pos < sp_masks.size(); pos++)
    {
        const auto extracted = manip.extract_prefix_suffix(s, pos);
        EXPECT_EQ(extracted, s.m_pair & sp_masks[pos]) << "pos=" << pos;
    }
}

TEST(SkmerManipulator, ConcatenateSkmerPinsCurrentBehavior)
{
    // concatenate_skmer uses &= (not |=) and increments m_suff_size by 1.
    // This test pins the current observable behavior so any future change is deliberate.
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    km::Skmer<kuint> a{kpair{0xFFFFu, 0u}, 3, 2};
    km::Skmer<kuint> b{kpair{0xA5A5u, 0u}, 3, 3};

    manip.concatenate_skmer(a, b);

    EXPECT_EQ(a.m_pair, kpair(static_cast<kuint>(0xA5A5u), 0u));
    EXPECT_EQ(a.m_suff_size, 3);
    EXPECT_EQ(a.m_pref_size, 3);
}

TEST(SkmerManipulator, IsForwardMatchesReturnedSkmer)
{
    using kuint = uint16_t;
    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    const string seq{"ACGTACGTACGT"};
    for (const auto& letter : seq)
    {
        kuint nucl{static_cast<kuint>((letter >> 1) & 0b11)};
        const auto& returned = manip.add_nucleotide(nucl);
        if (manip.is_forward())
            EXPECT_EQ(returned.m_pair, manip.m_fwd.m_pair);
        else
            EXPECT_EQ(returned.m_pair, manip.m_rev.m_pair);
    }
}


TEST(SkmerManipulator, MinimumKMPlusOne)
{
    // Sanity check that the minimum skmer (k=m+1, so k-m=1) builds its mask
    // tables with the right sizes and enumerates nucleotides without crashing.
    using kuint = uint8_t;
    constexpr uint64_t k{2};
    constexpr uint64_t m{1};

    km::SkmerManipulator<kuint> manip{k, m};
    EXPECT_EQ(manip.get_sp_mask().size(), k - m + 2u);
    EXPECT_EQ(manip.get_k_mask().size(), k - m + 1u);
    EXPECT_EQ(manip.get_n_mask().size(), 2u * k - m);

    const string seq{"ACGT"};
    for (const auto& letter : seq)
    {
        kuint nucl{static_cast<kuint>((letter >> 1) & 0b11)};
        manip.add_nucleotide(nucl);
    }
}

TEST(SkmerManipulator, EnumerateAcrossKuintBoundary)
{
    using kuint = uint32_t;
    constexpr uint64_t k{10};
    constexpr uint64_t m{2};

    // 2 * sk_size bits of state > one kuint, so m_value[1] must be used.
    ASSERT_GT(2u * (2u*k - m), sizeof(kuint) * 8u);

    km::SkmerManipulator<kuint> manip{k, m};
    const string seq{"ACGTACGTACGTACGTACGT"};
    for (const auto& letter : seq)
    {
        kuint nucl{static_cast<kuint>((letter >> 1) & 0b11)};
        manip.add_nucleotide(nucl);
    }

    const kuint high_engaged = manip.m_fwd.m_pair.m_value[1] | manip.m_rev.m_pair.m_value[1];
    EXPECT_NE(high_engaged, static_cast<kuint>(0));

    const auto k_masks = manip.get_k_mask();
    const km::Skmer<kuint> full = manip.m_fwd;
    for (uint64_t pos{0}; pos <= k - m; pos++)
    {
        const auto extracted = manip.get_skmer_of_kmer(full, pos);
        EXPECT_EQ(extracted.m_pair & k_masks[pos], full.m_pair & k_masks[pos])
            << "kmer_pos=" << pos;
    }
}