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