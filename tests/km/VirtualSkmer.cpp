#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <array>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

using namespace std;
using kuint = uint16_t;
using kpair = km::Skmer<kuint>::pair;

std::array< std::array< km::Skmer<kuint>, 3>, 6> get_skmer_permutations (std::array< kpair, 3> const & kmer_triplet)
{

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};

    const std::array< std::string, 3 > skmer_strings {"AAAAA", "AACCC", "CCCCC"};

    std::array< std::array< km::Skmer<kuint>, 3>, 6> permutation_array {};

    permutation_array[0] = {km::Skmer<kuint>(kmer_triplet[0],0,3),
                            km::Skmer<kuint>(kmer_triplet[1],0,3),
                            km::Skmer<kuint>(kmer_triplet[2],0,3)};
    // for(auto el: permutation_array[0]){
    //     pp << el;
    //     std::cout << el << std::endl;
    // }

    for(int64_t i {1}; i < 6; i++){
        permutation_array[i] = permutation_array[i-1];
        std::next_permutation(permutation_array[i].begin(), permutation_array[i].end());
    }

    return permutation_array;
}


/** Testing the result of "has_valid_kmer" on 2 skmers.
 * TODO: move to the manipulator test file.
 */
TEST(VirtualSkmer, kmer_validation)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};

    //                  Prefix:         A   T   _   _             A   _   _   _   
    //                  Suffix:       A   C   C   C             C   C   C   C     
    const kpair input_skmers[2] { {0b0000011001110111U, 0}, {0b0100011101110111U, 0}} ;
    std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0],2,4), km::Skmer<kuint>(input_skmers[1],1,4)};

    const uint64_t kmer_positions {k - m + 1};

    const uint64_t expected_valid_kmers[kmer_positions][2]{{0,0},{0,0},{1,0},{1,1}};
    // 0 values map to false, else to true
    bool kmer_validity;
    for (uint64_t skmer_id {0}; skmer_id < 2; skmer_id++ ){
        for(uint64_t position{0}; position < kmer_positions; position++ ){
            kmer_validity = manip.has_valid_kmer(skmer_vector[skmer_id], position);
            ASSERT_EQ(kmer_validity, expected_valid_kmers[position][skmer_id]) << (string("has_valid_kmer returned an unexpected value for skmer ") + std::to_string(skmer_id) + " at position " + std::to_string(position));
        }
    }
}

/** Test the order of 2 kmers after a sort on one column */
TEST(VirtualSkmer, Single_kmer_sorting)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                 Prefix:          A   T   _   _             A   _   _   _   
    //                 Suffix:        A   C   C   C             C   C   C   C     
    const kpair input_skmers[2] { {0b0000011001110111U, 0}, {0b0100011101110111U, 0}} ;
    const uint64_t position {3};
    

    std::vector<km::Skmer<kuint> > skmer_vector{km::Skmer<kuint>(input_skmers[0],1,4), km::Skmer<kuint>(input_skmers[1],2,4)};
    std::vector<uint64_t> expected_order {0, 1};
    std::vector<uint64_t> ordered_kmers = km::sorting::sort_column(skmer_vector.begin(), skmer_vector.end(), position, manip);
    ASSERT_EQ(ordered_kmers,expected_order);

    std::vector<km::Skmer<kuint> > skmer_vector_rev{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};
    std::vector<uint64_t> expected_order_rev {1, 0};
    std::vector<uint64_t> ordered_kmers_rev = km::sorting::sort_column(skmer_vector_rev.begin(), skmer_vector_rev.end(), position, manip);
    ASSERT_EQ(ordered_kmers_rev, expected_order_rev);
}

/** Test sorting on a column for all possible permutations of 3 skmers */
TEST(VirtualSkmer, Three_kmer_sorting)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};

        //                      Prefix:                        A   _   _   _                 A   _   _   _  
    //                      Suffix:                      A   A   A   A                 A   C   C   C 
    std::array< kpair, 3 > const kmer_triplet { kpair(0b0000001100110011U,0), kpair(0b0000011101110111U,0),
        //                      Prefix:                        C   _   _   _
        //                      Suffix:                      C   C   C   C 
                                                    kpair(0b0101011101110111U,0) };

    std::vector<uint64_t> ordered_kmers {};

    uint64_t position {3};
    for(auto const & permuted: get_skmer_permutations(kmer_triplet)){
        ordered_kmers = km::sorting::sort_column(permuted.begin(), permuted.end(), position, manip);
        uint64_t loop_idx {0};
        
        for(auto & skmer_position: ordered_kmers){
            km::Skmer<kuint> const & curr_skmer { permuted[skmer_position] };
            kpair const & expected_value {kmer_triplet[loop_idx]};
            ASSERT_EQ(curr_skmer.m_pair,expected_value);

            loop_idx+=1;
        }
    }
}


/** Compare 2 kmers from to adjacent columns and try to overlap them.
 * No overlap is expected in this example.
 */
TEST(VirtualSkmer, get_candidate_overlaps__no_overlap)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:       A   C   _   _             A   _   _   _   
    //                    Suffix:     A   C   C   _             A   G   C   C     
    const kpair input_skmers[2] { {0b0000010101111111U, 0}, {0b0000111101110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0};
    std::vector<uint64_t> right_column_order {1};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {};

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        std::cerr << "Unexpected overlap: " << computed_overlaps[i].first << computed_overlaps[i].second << std::endl;
    }
    
    ASSERT_EQ(computed_overlaps.size(),0);

}


/** Compare 2 kmers from to adjacent columns and try to overlap them.
 * 1 overlap is expected in this example.
 */
TEST(VirtualSkmer, get_candidate_overlaps__1_overlap)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:       A   C   _   _             A   _   _   _   
    //                    Suffix:     A   C   C   _             A   C   C   C     
    const kpair input_skmers[2] { {0b0000010101111111U, 0}, {0b0000011101110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0};
    std::vector<uint64_t> right_column_order {1};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,1) };

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

/** Overlap test with 2 columns.
 * 1 kmer in the left column and 2 in the right column.
 */
TEST(VirtualSkmer, get_candidate_overlaps_1_2)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   _   _   _             A   _   _   _   
    //                    Suffix:    A   C   C   _             A   C   C   C             A   C   C   T
    const kpair input_skmers[3] { {0b0000010101111111U, 0}, {0b0000011101110111U, 0}, {0b0000011101111011U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],1,4),km::Skmer<kuint>(input_skmers[2],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0};
    std::vector<uint64_t> right_column_order {1,2};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,1),std::pair<uint64_t,uint64_t>(0,2)};

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

/** Overlap test with 2 columns.
 * 2 kmers in the left column and 1 in the right column.
 */
TEST(VirtualSkmer, get_candidate_overlaps_2_1)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   T   _   _             A   _   _   _   
    //                    Suffix:    A   C   C   _             A   C   C   _             A   C   C   C
    const kpair input_skmers[3] { {0b0000010101111111U, 0}, {0b0000011001111111U, 0}, {0b0000011101110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[2],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0,1};
    std::vector<uint64_t> right_column_order {2};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,2),std::pair<uint64_t,uint64_t>(1,2)};

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

/** Overlap test with 2 columns.
 * 2 kmers in the left column and 2 in the right column.
 * Expecting 2 parallel overlaps.
 */
TEST(VirtualSkmer, get_candidate_overlaps_2_2_parallel)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _
    //                    Suffix:    A   C   C   _             A   C   T   _             A   C   C   C
    const kpair input_skmers[4] { {0b0000010101111111U, 0}, {0b0000011010111111U, 0}, {0b0000011101110111U, 0},
    //                                 A   _   _   _   
    //                               A   C   T   C
                                  {0b0000011110110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],2,3), km::Skmer<kuint>(input_skmers[2],1,4), km::Skmer<kuint>(input_skmers[3],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0,1};
    std::vector<uint64_t> right_column_order {2,3};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,2),std::pair<uint64_t,uint64_t>(1,3)};

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

/** Overlap test with 2 columns.
 * 2 kmers in the left column and 2 in the right column.
 * Expecting 2 crossing overlaps.
 */
TEST(VirtualSkmer, get_candidate_overlaps_2_2_crossed)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _             
    //                    Suffix:    A   C   T   _             A   C   C   _             A   C   C   C                 
    const kpair input_skmers[4] { {0b0000010110111111U, 0}, {0b0000011001111111U, 0}, {0b0000011101110111U, 0}, 
    //                                 A   _   _   _   
    //                               A   C   T   C
                                  {0b0000011110110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],2,3), km::Skmer<kuint>(input_skmers[2],1,4), km::Skmer<kuint>(input_skmers[3],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0,1};
    std::vector<uint64_t> right_column_order {2,3};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,3),std::pair<uint64_t,uint64_t>(1,2)};

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

/** Overlap test with 2 columns.
 * 2 kmers in the left column and 2 in the right column.
 * Expecting 1 crossing overlaps.
 */
TEST(VirtualSkmer, get_candidate_overlaps_2_2_crossed_1)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _             
    //                    Suffix:    A   C   A   _             A   C   C   _             A   C   C   C                 
    const kpair input_skmers[4] { {0b0000010100111111U, 0}, {0b0000011001111111U, 0}, {0b0000011101110111U, 0}, 
    //                                 A   _   _   _   
    //                               A   C   T   C
                                  {0b0000011110110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],2,3), km::Skmer<kuint>(input_skmers[2],1,4), km::Skmer<kuint>(input_skmers[3],1,4)};
    const uint64_t left_column_position {2};
    std::vector<uint64_t> left_column_order {0,1};
    std::vector<uint64_t> right_column_order {2,3};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(1,2)};

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

/** Overlap test with 2 columns.
 * 2 kmers in the left column and 2 in the right column.
 * Expecting 2 crossing overlaps.
 */
TEST(VirtualSkmer, get_candidate_overlaps_2_2_crossed_beginning)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   C   C             A   T   T   T             A   T   T   _
    //                    Suffix:    A   _   _   _             A   _   _   _             A   C   _   _
    const kpair input_skmers[4] { {0b0000110111011101U, 0}, {0b0000111011101110U, 0}, {0b0000011011101111U, 0}, 
    //                                 A   C   C   _   
    //                               A   T   _   _
                                  {0b0000100111011111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],2,3), km::Skmer<kuint>(input_skmers[2],1,4), km::Skmer<kuint>(input_skmers[3],1,4)};
    const uint64_t left_column_position {0};
    std::vector<uint64_t> left_column_order {0,1};
    std::vector<uint64_t> right_column_order {2,3};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,3),std::pair<uint64_t,uint64_t>(1,2)};

    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}


/** Overlap test with 2 columns.
 * 2 kmers in the left column and 2 in the right column.
 * Expecting 2 crossing overlaps.
 */
TEST(VirtualSkmer, get_candidate_overlaps_2_2_crossed_end)
{   
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _             
    //                    Suffix:    A   C   T   C             A   C   C   C             A   C   C   C                 
    const kpair input_skmers[4] { {0b0000010110110111U, 0}, {0b0000011001110111U, 0}, {0b0000011101110111U, 0}, 
    //                                 A   _   _   _   
    //                               A   C   T   C
                                  {0b0000011110110111U, 0}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],2,3), km::Skmer<kuint>(input_skmers[2],1,4), km::Skmer<kuint>(input_skmers[3],1,4)};
    const uint64_t left_column_position {3};
    std::vector<uint64_t> left_column_order {0,1};
    std::vector<uint64_t> right_column_order {2,3};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {std::pair<uint64_t,uint64_t>(0,3),std::pair<uint64_t,uint64_t>(1,2)};
    
    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());

    for (uint64_t i {0}; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

TEST(VirtualSkmer, generate_virtual_skmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   G   T             
    //                    Suffix:    A   C   T   C                         
    const kpair input_skmer {0b0000010110110110U, 0};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmer,4,4)};
    uint64_t skmer_id {0};
    std::vector< km::sorting::Virtual_skmer<kuint> > extracted_skmers;

    for(int position = 0; position < 4; position++){
        extracted_skmers.push_back(km::sorting::generate_virtual_skmer(m_skmer_vector,manip,skmer_id,position));
    }

    //                    Prefix:      A   C   G   T             A   C   G   _            A   C   _   _             
    //                    Suffix:    A   _   _   _             A   C   _   _             A   C   T   _ 
    const kpair expected_kpairs[4] {{0b0000110111111110U, 0}, {0b0000010111111111U, 0}, {0b0000010110111111U, 0}, 
    //                                 A   _   _   _   
    //                               A   C   T   C
                                    {0b0000011110110111U, 0}};
    std::vector< km::sorting::Virtual_skmer<kuint> > expected_virtual_skmers {km::sorting::Virtual_skmer<kuint>(expected_kpairs[0],4,1,0),km::sorting::Virtual_skmer<kuint>(expected_kpairs[1],3,2,0),km::sorting::Virtual_skmer<kuint>(expected_kpairs[2],2,3,0),km::sorting::Virtual_skmer<kuint>(expected_kpairs[3],1,4,0) };

    for(size_t i {0}; i < expected_virtual_skmers.size(); i++){
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pair,extracted_skmers[i].skmer.m_pair);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pref_size,extracted_skmers[i].skmer.m_pref_size);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_suff_size,extracted_skmers[i].skmer.m_suff_size);
        // ASSERT_EQ(expected_virtual_skmers[i].last_id, extracted_skmers[i].last_id);
    }

}

TEST(VirtualSkmer, generate_virtual_skmer_2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{3};

    km::SkmerManipulator<kuint> manip {k, m};
    // km::SkmerPrettyPrinter<kuint> pp {k, m};

    //              Prefix:      A   C   G   T             
    //               Suffix:       C   T   C                         
    const kpair input_skmer {0b0000010110110110U, 0};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmer,4,3)};

    std::vector< km::sorting::Virtual_skmer<kuint> > extracted_skmers;

    for(int i = 0; i < 3; i++){
        extracted_skmers.push_back(km::sorting::generate_virtual_skmer(m_skmer_vector,manip,0,i));
    }

    //                    Prefix:      A   C   G   T             A   C   G   _            A   C   _   _             
    //                    Suffix:        C   _   _                 C   T   _                C   T   C
    const kpair expected_kpairs[3] {{0b0000010111111110U, 0}, {0b0000010110111111U, 0}, {0b0000010110110111U, 0}}; 
 
    std::vector< km::sorting::Virtual_skmer<kuint> > expected_virtual_skmers {km::sorting::Virtual_skmer<kuint>(expected_kpairs[0],4,1,0), km::sorting::Virtual_skmer<kuint>(expected_kpairs[1],3,2,0),km::sorting::Virtual_skmer<kuint>(expected_kpairs[2],2,3,0)};

    for(size_t i {0}; i < expected_virtual_skmers.size(); i++){
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pair,extracted_skmers[i].skmer.m_pair);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pref_size,extracted_skmers[i].skmer.m_pref_size);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_suff_size,extracted_skmers[i].skmer.m_suff_size);
    }
}

TEST(VirtualSkmer, add_kmer_to_virtual_skmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   G   T             
    //                    Suffix:    A   C   T   C                         
    const kpair m_vector_pair {0b0000010110110110U, 0};

    //                    Prefix:      A   C   G   _             
    //                    Suffix:    A   C   _   _      
    const kpair m_basic_pair {0b0000010111111111U, 0};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(m_vector_pair,4,4),km::Skmer<kuint>(m_basic_pair,3,2)};


    std::vector< km::sorting::Virtual_skmer<kuint> > extracted_skmers;
    uint64_t skmer_id_take_skmer {0};
    uint64_t skmer_id_base_skmer {1};
    uint64_t kmers_to_add[2] {2,3};

    km::sorting::Virtual_skmer<kuint> new_virtual_skmer(m_skmer_vector[1],skmer_id_base_skmer);
    for(uint64_t kmer_pos : kmers_to_add){
        km::sorting::add_kmer_virtual_skmer(new_virtual_skmer, m_skmer_vector, manip, skmer_id_take_skmer, kmer_pos);
        extracted_skmers.push_back(new_virtual_skmer);
    }

    //                      Prefix:      A   C   G   _             A   C   G   _             
    //                      Suffix:    A   C   T   _             A   C   T   C       
    const kpair expected_kpairs[4] {{0b0000010110111111U, 0}, {0b0000010110110111U, 0}};
 
    // km::sorting::Virtual_skmer<kuint>(expected_kpairs[0],4,2,0),
    std::vector< km::sorting::Virtual_skmer<kuint> > expected_virtual_skmers { km::sorting::Virtual_skmer<kuint>(expected_kpairs[0],3,3,0),km::sorting::Virtual_skmer<kuint>(expected_kpairs[1],3,4,0)};

    for(size_t i {0}; i < expected_virtual_skmers.size(); i++){
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pair,extracted_skmers[i].skmer.m_pair);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pref_size,extracted_skmers[i].skmer.m_pref_size);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_suff_size,extracted_skmers[i].skmer.m_suff_size);
        // ASSERT_EQ(expected_virtual_skmers[i].last_id, extracted_skmers[i].last_id);
    }
}

TEST(VirtualSkmer, empty_list_fill)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   _   _             A   T   _   _                   
    //                    Suffix:    A   C   T   C             A   C   C   C                          
    const kpair input_skmers[2] { {0b0000010110110111U, 0}, {0b0000011001110111U, 0}};
                                 //  0000010110111111          0000011001111111
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,3), km::Skmer<kuint>(input_skmers[1],2,3)};
    const std::vector<uint64_t> column{0,1};
    std::forward_list<km::sorting::Virtual_skmer<kuint>> computed_list;
    const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps;
    uint64_t column_position {2};
    km::sorting::merge_LList_column(m_skmer_vector, manip, computed_list, column, valid_overlaps, column_position);

    km::sorting::print_list(computed_list);

    //                       Prefix:      A   C   _   _             A   T   _   _                   
    //                       Suffix:    A   C   T   _             A   C   C   _                          
    const kpair expected_skmers[2] { {0b0000010110111111, 0}, {0b0000011001111111, 0}};
    const std::array<std::pair<uint64_t, uint64_t>, 2> prefix_suffix {{std::pair(2,3),std::pair(2,3)}};
    const std::vector<km::sorting::Virtual_skmer<kuint>> expected_list {km::sorting::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second , column[0]),km::sorting::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second , column[1])};

    auto vector_iterator = expected_list.begin();
    auto list_iterator = computed_list.begin();

    for(; vector_iterator != expected_list.end(); vector_iterator++){
        km::sorting::Virtual_skmer<kuint> expected_vskmer {*vector_iterator};
        km::sorting::Virtual_skmer<kuint> computed_vskmer {*list_iterator};

        ASSERT_EQ(expected_vskmer.skmer.m_pair, computed_vskmer.skmer.m_pair);
        ASSERT_EQ(expected_vskmer.skmer.m_pref_size, computed_vskmer.skmer.m_pref_size);
        ASSERT_EQ(expected_vskmer.skmer.m_suff_size, computed_vskmer.skmer.m_suff_size);
        ASSERT_EQ(expected_vskmer.last_id, computed_vskmer.last_id);
        
        list_iterator++;
    }
}

TEST(VirtualSkmer, merge_left_and_right_column_elements)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G     
    //                    Suffix:    A   G   G   G             A   T   G   T             A   C   T   T                
    const kpair input_skmers[3] { {0b0000110111011101U, 0}, {0b0000100111011010U, 0}, {0b0000010110011011U, 0}};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],4,4), km::Skmer<kuint>(input_skmers[1],4,4), km::Skmer<kuint>(input_skmers[2],4,4)};

    const std::vector<uint64_t> column{1,2};

    std::forward_list<km::sorting::Virtual_skmer<kuint>> computed_list {km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],4,1),0),0), km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1],4,1),0),1)};
    // km::Skmer<kuint> sk1 {manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],0))};

    // computed_list. km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],0)), 4, 1), 0), km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1],0), 4, 1), 1)};


    const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps {std::pair(0,1),std::pair(1,2)};
    uint64_t column_position {1};
    km::sorting::merge_LList_column(m_skmer_vector, manip, computed_list, column, valid_overlaps, column_position);

    km::sorting::print_list(computed_list);

    //                       Prefix:      A   C   C   C            A   C   C   T                   
    //                       Suffix:    A   T   _   _            A   C   _   _                          
    const kpair expected_skmers[2] { {0b0000100111011101, 0}, {0b0000010111011110, 0}};
    const std::array<std::pair<uint64_t, uint64_t>, 2> prefix_suffix {{std::pair(4,2), std::pair(4,2)}};
    const std::vector<km::sorting::Virtual_skmer<kuint>> expected_list {km::sorting::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second , column[0]),km::sorting::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second , column[1])};

    auto vector_iterator = expected_list.begin();
    auto list_iterator = computed_list.begin();

    for(; vector_iterator != expected_list.end(); vector_iterator++){
        km::sorting::Virtual_skmer<kuint> expected_vskmer {*vector_iterator};
        km::sorting::Virtual_skmer<kuint> computed_vskmer {*list_iterator};

        ASSERT_EQ(expected_vskmer.skmer.m_pair, computed_vskmer.skmer.m_pair);
        ASSERT_EQ(expected_vskmer.skmer.m_pref_size, computed_vskmer.skmer.m_pref_size);
        ASSERT_EQ(expected_vskmer.skmer.m_suff_size, computed_vskmer.skmer.m_suff_size);
        ASSERT_EQ(expected_vskmer.last_id, computed_vskmer.last_id);
        
        list_iterator++;
    }
}


TEST(VirtualSkmer, add_right_column_element)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G     
    //                    Suffix:    A   G   G   G             A   G   T   T             A   G   T   T       
    const kpair input_skmers[3] { {0b0000110111011101U, 0}, {0b0000110110011010U, 0}, {0b0000110110011011U, 0}};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],4,4), km::Skmer<kuint>(input_skmers[1],4,4), km::Skmer<kuint>(input_skmers[2],4,4)};

    const std::vector<uint64_t> column{1,2};

    std::forward_list<km::sorting::Virtual_skmer<kuint>> computed_list {km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],3,2),1),0), km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1],3,2),1),1)};

    const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps {std::pair(0,2)};
    uint64_t column_position {2};
    km::sorting::merge_LList_column(m_skmer_vector, manip, computed_list, column, valid_overlaps, column_position);


    km::sorting::print_list(computed_list);

    //                       Prefix:      A   C   _   _            A   C   C   _            A   C   C   _
    //                       Suffix:    A   G   T   _            A   G   T   _            A   G   _   _
    const kpair expected_skmers[3] { {0b0000110110111111, 0}, {0b0000110110011111, 0}, {0b0000110111011111, 0}};
    const std::array<std::pair<uint64_t, uint64_t>, 3> prefix_suffix {{std::pair(2,3), std::pair(3,3), std::pair(3,2)}};
    const std::vector<km::sorting::Virtual_skmer<kuint>> expected_list {km::sorting::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second , column[0]), km::sorting::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second , column[1]), km::sorting::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second , column[0])};

    auto vector_iterator = expected_list.begin();
    auto list_iterator = computed_list.begin();

    for(; vector_iterator != expected_list.end(); vector_iterator++){
        km::sorting::Virtual_skmer<kuint> expected_vskmer {*vector_iterator};
        km::sorting::Virtual_skmer<kuint> computed_vskmer {*list_iterator};

        ASSERT_EQ(expected_vskmer.skmer.m_pair, computed_vskmer.skmer.m_pair);
        ASSERT_EQ(expected_vskmer.skmer.m_pref_size, computed_vskmer.skmer.m_pref_size);
        ASSERT_EQ(expected_vskmer.skmer.m_suff_size, computed_vskmer.skmer.m_suff_size);
        ASSERT_EQ(expected_vskmer.last_id, computed_vskmer.last_id);
        
        list_iterator++;
    }
}

TEST(VirtualSkmer, add_left_column_element)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G     
    //                    Suffix:    A   G   G   G             A   G   T   T             A   G   T   T       
    const kpair input_skmers[3] { {0b0000110111011101U, 0}, {0b0000110110011010U, 0}, {0b0000110110011011U, 0}};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],4,4), km::Skmer<kuint>(input_skmers[1],4,4), km::Skmer<kuint>(input_skmers[2],4,4)};

    const std::vector<uint64_t> column{1,2};

    std::forward_list<km::sorting::Virtual_skmer<kuint>> computed_list {km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],2,3),2),0), km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1],2,3),2),1)};

    const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps {std::pair(1,1)};
    uint64_t column_position {3};
    km::sorting::merge_LList_column(m_skmer_vector, manip, computed_list, column, valid_overlaps, column_position);


    km::sorting::print_list(computed_list);

    //                       Prefix:      A   C   _   _            A   C   _   _            A   _   _   _
    //                       Suffix:    A   G   G   _            A   G   T   T            A   G   T   T
    const kpair expected_skmers[3] { {0b0000110111111111, 0}, {0b0000110110111011, 0}, {0b0000111110111011, 0}};
    const std::array<std::pair<uint64_t, uint64_t>, 3> prefix_suffix {{std::pair(2,3), std::pair(2,4), std::pair(1,4)}};
    const std::vector<km::sorting::Virtual_skmer<kuint>> expected_list {km::sorting::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second , 0), km::sorting::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second , 1), km::sorting::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second , 2)};

    auto vector_iterator = expected_list.begin();
    auto list_iterator = computed_list.begin();

    for(; vector_iterator != expected_list.end(); vector_iterator++){
        km::sorting::Virtual_skmer<kuint> expected_vskmer {*vector_iterator};
        km::sorting::Virtual_skmer<kuint> computed_vskmer {*list_iterator};

        ASSERT_EQ(expected_vskmer.skmer.m_pair, computed_vskmer.skmer.m_pair);
        ASSERT_EQ(expected_vskmer.skmer.m_pref_size, computed_vskmer.skmer.m_pref_size);
        ASSERT_EQ(expected_vskmer.skmer.m_suff_size, computed_vskmer.skmer.m_suff_size);
        ASSERT_EQ(expected_vskmer.last_id, computed_vskmer.last_id);
        
        list_iterator++;
    }
}

TEST(VirtualSkmer, no_element_pointed)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    km::SkmerPrettyPrinter<kuint> pp {k, m};
    //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G     
    //                    Suffix:    A   A   G   G             A   G   T   T             A   A   T   T       
    const kpair input_skmers[3] { {0b0000000111011101U, 0}, {0b0000110110011010U, 0}, {0b0000000110011011U, 0}};
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[0],4,4), km::Skmer<kuint>(input_skmers[1],4,4), km::Skmer<kuint>(input_skmers[2],4,4)};

    const std::vector<uint64_t> column{2,1};

    std::forward_list<km::sorting::Virtual_skmer<kuint>> computed_list {km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],2,3),2),0), km::sorting::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1],2,3),2),1)};

    const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps {};
    uint64_t column_position {3};
    km::sorting::merge_LList_column(m_skmer_vector, manip, computed_list, column, valid_overlaps, column_position);


    km::sorting::print_list(computed_list);

    //                       Prefix:      A   C   _   _            A   _   _   _            A   C   _   _
    //                       Suffix:    A   A   G   _            A   A   T   T            A   G   T   _
    const kpair expected_skmers[4] { {0b0000000111111111, 0}, {0b0000001110111011, 0}, {0b0000110110111111, 0}, {0b0000111110111011, 0}};
    const std::array<std::pair<uint64_t, uint64_t>, 4> prefix_suffix {{std::pair(2,3), std::pair(1,4), std::pair(2,3), std::pair(1,4)}};
    const std::vector<km::sorting::Virtual_skmer<kuint>> expected_list {km::sorting::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second , 0), km::sorting::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second , 2), km::sorting::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second , 1), km::sorting::Virtual_skmer<kuint>(expected_skmers[3], prefix_suffix[3].first, prefix_suffix[3].second , 1)};

    auto vector_iterator = expected_list.begin();
    auto list_iterator = computed_list.begin();

    for(; vector_iterator != expected_list.end(); vector_iterator++){
        km::sorting::Virtual_skmer<kuint> expected_vskmer {*vector_iterator};
        km::sorting::Virtual_skmer<kuint> computed_vskmer {*list_iterator};

        ASSERT_EQ(expected_vskmer.skmer.m_pair, computed_vskmer.skmer.m_pair);
        ASSERT_EQ(expected_vskmer.skmer.m_pref_size, computed_vskmer.skmer.m_pref_size);
        ASSERT_EQ(expected_vskmer.skmer.m_suff_size, computed_vskmer.skmer.m_suff_size);
        ASSERT_EQ(expected_vskmer.last_id, computed_vskmer.last_id);
        
        list_iterator++;
    }
}