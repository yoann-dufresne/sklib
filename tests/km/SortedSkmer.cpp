#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/SkmerSorting.hpp>

using namespace std;

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
const std::array< std::string, 3 > skmer_strings {"AAAAA", "AACCC", "CCCCC"};


std::array< std::array< km::Skmer<kuint>, 3>, 6> get_skmer_permutations (std::array< kpair, 3> const & kmer_triplet)
{

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
TEST(SkmerSorting, kmer_validation)
{
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
TEST(SkmerSorting, Single_kmer_sorting)
{
    
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
TEST(SkmerSorting, Three_kmer_sorting)
{
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
TEST(SkmerSorting, get_candidate_overlaps__no_overlap)
{   
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
TEST(SkmerSorting, get_candidate_overlaps__1_overlap)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_1_2)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_2_1)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_2_2_parallel)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_2_2_crossed)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_2_2_crossed_1)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_2_2_crossed_beginning)
{   
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
TEST(SkmerSorting, get_candidate_overlaps_2_2_crossed_end)
{   
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