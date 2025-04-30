#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <array>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <io/io.hpp>
#include <algorithms/SkmerSorting.hpp>

using namespace std;

using kuint = uint16_t;
using kpair = km::Skmer<kuint>::pair;

constexpr uint64_t k{5};
constexpr uint64_t m{2};

km::SkmerManipulator<kuint> manip {k, m};
km::SkmerPrettyPrinter<kuint> pp {k, m};


/** TESTING WRITING AND READING OF A SIMPLE VECTOR.
 * TODO: move to the manipulator test file.
 */
TEST(IO, Write_Read)
{
    //                  Prefix:         A   T   _   _             A   _   _   _   
    //                  Suffix:       A   C   C   C             C   C   C   C     
    const kpair input_skmers[2] { {0b0000011001110111U, 0}, {0b0100011101110111U, 0}} ;
    std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0],2,4), km::Skmer<kuint>(input_skmers[1],1,4)};

    std::forward_list<km::sorting::Virtual_skmer<kuint>> v_skmer_list{km::sorting::Virtual_skmer<kuint>(skmer_vector[0],0), km::sorting::Virtual_skmer<kuint>(skmer_vector[1],1)};

    int verification_writing = km::write_sorted_list(manip, v_skmer_list, "here.bin");
    ASSERT_EQ(verification_writing, 0) ;

    km::SkmerManipulator<kuint> read_manip;
    std::vector<km::sorting::Virtual_skmer<kuint>> read_vector;

    int verification_reading = km::load_sorted_vector(read_manip, read_vector, "here.bin");
    ASSERT_EQ(verification_reading, 0) ;

    ASSERT_EQ(manip.k, read_manip.k);
    ASSERT_EQ(manip.m, read_manip.m);

    std::forward_list<km::sorting::Virtual_skmer<kuint>> read_v_skmer_list(read_vector.begin(), read_vector.end()) ;
    auto it1 = read_v_skmer_list.begin();
    auto it2 = v_skmer_list.begin();

    while(it1 != read_v_skmer_list.end() || it2 != read_v_skmer_list.end() ){
        ASSERT_EQ(*it1.skmer, *it2.skmer);
        ASSERT_EQ(*it1.last_id, *it2.last_id);
        it1++; it2++;
    }

}