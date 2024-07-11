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

    //                      Prefix:      A   _   _   _            A   _   _   _            C   _   _   _
    //                      Suffix:    A   A   A   A            A   C   C   C            C   C   C   C 
const std::array< kpair, 3 > kmer_triplet { kpair(0b0000001100110011U,0), kpair(0b0000011101110111U,0), kpair(0b0101011101110111U,0) };

const std::array< std::string, 3 > skmer_strings {"AAAAA", "AACCC", "CCCCC"};


std::array< std::array< km::Skmer<kuint>, 3>, 6> get_skmer_permutations (std::array< kpair, 3> const & kmer_triplet)
{

    std::array< std::array< km::Skmer<kuint>, 3>, 6> permutation_array {};

    permutation_array[0] = {km::Skmer<kuint>(kmer_triplet[0],0,3),
                            km::Skmer<kuint>(kmer_triplet[1],0,3),
                            km::Skmer<kuint>(kmer_triplet[2],0,3)};
    for(auto el: permutation_array[0]){
        pp << el;
        std::cout << el << std::endl;
    }

    for(int64_t i {1}; i < 6; i++){
        permutation_array[i] = permutation_array[i-1];
        std::next_permutation(permutation_array[i].begin(), permutation_array[i].end());
    }

    return permutation_array;
}



TEST(SkmerSorting, kmer_validation)
{
    //                         Prefix:         A   T   _   _             A   _   _   _   
    //                         Suffix:       A   C   C   C             C   C   C   C     
    const kpair input_skmers[2] { {0, 0b0000011001110111U}, {0, 0b0100011101110111U}} ;
    // km::Skmer<kuint> sk2(input_skmers[1],1,4);
    //const km::Skmer<kuint> skmers_array[2] = {sk1,sk2} ;
    std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0],2,4), km::Skmer<kuint>(input_skmers[1],1,4)};
    // std::cout << "Before change: " <<  m_skmer_vector[0].m_pref_size << std::endl;
    // m_skmer_vector[0].m_pref_size = 5;
    // std::cout << "After change: " << m_skmer_vector[0].m_pref_size << std::endl;
    // m_skmer_vector.push_back(sk1);


    const uint64_t kmer_positions {(2*k - m)};

    const uint64_t expected_valid_kmers[kmer_positions][2]{{0,0},{0,0},{1,0},{1,1},{1,1},{1,1},{1,1},{1,1}};
    // 0 values map to false, else to true
    bool kmer_validity;
    for (uint64_t skmer_id {0}; skmer_id < 2; skmer_id++ ){
        for(uint64_t position; position < kmer_positions; position++ ){
            kmer_validity = manip.has_valid_kmer(m_skmer_vector[skmer_id], position);
            ASSERT_EQ(kmer_validity, expected_valid_kmers[position][skmer_id]);
        }
    }
}

TEST(SkmerSorting, Single_kmer_sorting)
{
    
    //                         Prefix:         A   T   _   _             A   _   _   _   
    //                         Suffix:       A   C   C   C             C   C   C   C     
    const kpair input_skmers[2] { {0, 0b0000011001110111U}, {0, 0b0100011101110111U}} ;
    
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};

    //const uint64_t kmer_positions {(2*k - m)};
    const uint64_t position {3};
    std::vector<uint64_t> ordered_kmers = km::sorting::sort_column(m_skmer_vector.begin(), m_skmer_vector.end(), position, manip);
    std::vector<uint64_t> expected_order {1, 0};

    ASSERT_EQ(ordered_kmers,expected_order);
}

TEST(SkmerSorting, Three_kmer_sorting)
{
    std::vector<uint64_t> ordered_kmers {};

    uint64_t position {3};
    for(auto const & permuted: get_skmer_permutations(kmer_triplet)){

        ordered_kmers = km::sorting::sort_column(permuted.begin(), permuted.end(), position, manip);
        uint64_t loop_idx {0};

        for (auto & el: permuted){
            pp << el;
            std:: cout << pp;
        }
        
        for(auto & skmer_position: ordered_kmers){
            km::Skmer<kuint> const & curr_skmer { permuted[skmer_position] };
            kpair const & expected_value {kmer_triplet[loop_idx]};
            ASSERT_EQ(curr_skmer.m_pair,expected_value);

            loop_idx+=1;
        }
    }
}

TEST(SkmerSorting, get_canidate_overlaps_0)
{   
    //                         Prefix:         A   C   C   _             A   G   C   _   
    //                         Suffix:       A   C   _   _             A   C   _   _     
    const kpair input_skmers[2] { {0, 0b0000010111011111U}, {0, 0b0000011111011111U}} ;
    std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};
    const uint64_t left_column_position {1};
    std::vector<uint64_t> left_column_order {0};
    std::vector<uint64_t> right_column_order {1};
    auto computed_overlaps { km::sorting::get_candidate_overlaps(m_skmer_vector, manip, left_column_position, left_column_order, right_column_order) };
    std::vector<std::pair<uint64_t, uint64_t> > expected_overlaps {};
    ASSERT_EQ(computed_overlaps.size(),expected_overlaps.size());
    for (uint64_t i; i < computed_overlaps.size(); i+=1 ){
        ASSERT_EQ(computed_overlaps[i],expected_overlaps[i]);
    }
}

// TEST(SkmerSorting, get_canidate_overlaps_1)
// {   
//     //                         Prefix:         A   C   C   _           
//     //                         Suffix:       A   C   C   _   
//     const kpair input_skmers[2] { {0, 0b0000010101011111U}} ;
//     std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4)};
//     const uint64_t position {1};

// }

// TEST(SkmerSorting, get_canidate_overlaps_2)
// {   
//     //                         Prefix:         A   C   C   _             A   C   _   _   
//     //                         Suffix:       A   C   C   _             A   C   G   _     
//     const kpair input_skmers[2] { {0, 0b0000010101011111U}, {0, 0b0000010111111111U}} ;
//     std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};
//     const uint64_t position {1};

// }

// TEST(SkmerSorting, get_canidate_overlaps_3)
// {   
//     //                         Prefix:         A   C   C   _             A   C   G   _   
//     //                         Suffix:       A   C   C   _             A   C   _   _     
//     const kpair input_skmers[2] { {0, 0b0000010101011111U}, {0, 0b0000010111111111U}} ;
//     std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};
//     const uint64_t position {1};

// }

// TEST(SkmerSorting, get_canidate_overlaps_4)
// {   
//     //                         Prefix:         A   G   C   _             A   G   C   _   
//     //                         Suffix:       A   C   G   _             A   G   C   _     
//     const kpair input_skmers[2] { {0, 0b0000011111011111U}, {0, 0b0000010111111111U}} ;
//     std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};
//     const uint64_t position {1};

// }

// TEST(SkmerSorting, get_canidate_overlaps_5)
// {   
//     //                         Prefix:    A   C   T   _             A   C   A   _   
//     //                         Suffix:  A   C   C   _             A   C   A   _     
//     const kpair input_skmers[2] { {0, 0b0000010101101111U}, {0, 0b0000010100001111U}} ;
//     std::vector<km::Skmer<kuint> > m_skmer_vector{km::Skmer<kuint>(input_skmers[1],1,4), km::Skmer<kuint>(input_skmers[0],2,4)};
//     const uint64_t position {1};

// }

// TEST(SkmerSorting, Full_span_three_kmer_sorting)
// {
//     std::vector<uint64_t> ordered_kmers {};
//     std::array< km::Skmer<kuint> , 3 > skmers {km::Skmer<kuint>(),km::Skmer<kuint>(),km::Skmer<kuint>()};
//     std::array< km::SkmerManipulator<kuint>, 3> manipulators {km::SkmerManipulator<kuint>(5,2),km::SkmerManipulator<kuint>(5,2),km::SkmerManipulator<kuint>(5,2)};
//     const std::array< kpair, 3 > clean_triplet { kpair(0b0000000000000000U,0), kpair(0b0000010001000100U,0), kpair(0b0101010001000100U,0) };

//     for (int pos {0}; pos < k; pos+=1)
//     {
//         for(int64_t skmer_id {0}; skmer_id < manipulators.size(); skmer_id+=1){
//             auto letter {skmer_strings[skmer_id][pos]};
//             kuint nucl {static_cast<kuint>((letter >> 1) & 0b11)};
//             skmers[skmer_id] = manipulators[skmer_id].add_nucleotide(nucl);
//         }
//     }

//     for (uint64_t position {3}; position >=0; position-=1)
//     {
//         std::array<kpair, 3> pairs {skmers[0].m_pair, skmers[1].m_pair, skmers[2].m_pair};
//         for(auto const & permuted: get_skmer_permutations(pairs)){

//             ordered_kmers = km::sorting::sort_column(permuted.begin(), permuted.end(), position, manip);
//             uint64_t loop_idx {0};

//             for (auto & el: permuted){
//                 pp << el;
//                 std:: cout << pp;
//             }
            
//             for(auto & skmer_position: ordered_kmers){
//                 km::Skmer<kuint> const & curr_skmer { permuted[skmer_position] };
//                 kpair const & expected_value {clean_triplet[loop_idx] << (2 * ( 3 - position))};
//                 ASSERT_EQ(curr_skmer.m_pair,expected_value);

//                 loop_idx+=1;
//             }
//         }
//         for(int64_t skmer_id {0}; skmer_id < manipulators.size(); skmer_id+=1){

//             skmers[skmer_id] = manipulators[skmer_id].add_empty_nucleotide();
//         }
//     }

// }