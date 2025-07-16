// #include <iostream>
// #include <cstdlib>
// #include <gtest/gtest.h>
// #include <string>
// #include <array>

// #include <io/Skmer.hpp>
// #include <io/Skmerator.hpp>
// #include <io/io.hpp>
// // #include <algorithms/SkmerSorting.hpp>

// using namespace std;

// using kuint = uint16_t;
// using kpair = km::Skmer<kuint>::pair;


// /** TESTING WRITING AND READING OF A SIMPLE VECTOR.
//  * TODO: move to the manipulator test file.
//  */
// TEST(IO, Write_Read)
// {
//     constexpr uint64_t k{5};
//     constexpr uint64_t m{2};

//     km::SkmerManipulator<kuint> manip {k, m};

//     const std::string store_file {"test.bin"};
//     const size_t max_megs {10};
//     //                  Prefix:         A   T   _   _             A   _   _   _   
//     //                  Suffix:       A   C   C   C             C   C   C   C     
//     const kpair input_skmers[2] { {0b0000011001110111U, 0}, {0b0100011101110111U, 0}} ;
//     std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0],2,4), km::Skmer<kuint>(input_skmers[1],1,4)};

//     std::forward_list<km::sorting::Virtual_skmer<kuint>> v_skmer_list{km::sorting::Virtual_skmer<kuint>(skmer_vector[0],0), km::sorting::Virtual_skmer<kuint>(skmer_vector[1],1)};

//     int verification_writing = km::io::write_sorted_list(manip, v_skmer_list, max_megs, store_file);
//     ASSERT_EQ(verification_writing, 0) ;

//     // km::SkmerManipulator<kuint> read_manip {0,0};
//     // std::vector<km::sorting::Virtual_skmer<kuint>> read_vector;
//     uint64_t m_k, m_m {0};

//     std::vector<km::Skmer<kuint>> read_vector;
//     km::io::load_sorted_vector<kuint>(store_file, m_k, m_m, read_vector);
//     // ASSERT_EQ(verification_reading, 0) ;

//     // km::SkmerManipulator<kuint> read_manip = std::get<0>(std::move(read_tuple));
//     // std::vector<km::sorting::Virtual_skmer<kuint>> read_vector = std::get<1>(std::move(read_tuple));

//     ASSERT_EQ(manip.k, m_k);
//     ASSERT_EQ(manip.m, m_m);

//     // std::forward_list<km::Skmer<kuint>> read_v_skmer_list(read_vector.begin(), read_vector.end()) ;
//     auto it1 = read_vector.begin();
//     auto it2 = v_skmer_list.begin();

//     while(it1 != read_vector.end() && it2 != v_skmer_list.end() ){
//         ASSERT_EQ((*it1).m_pair, (*it2).skmer.m_pair);
//         it1++; it2++;
//     }

// }