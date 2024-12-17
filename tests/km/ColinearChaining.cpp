#include <gtest/gtest.h>
#include <utility>
#include <iostream>

#include <algorithms/SkmerSorting.hpp>

using overlap = std::pair<uint64_t, uint64_t>;
using namespace km::sorting;
using namespace std;


/** Single pair
 */
TEST(ColinearChaining, single_overlap)
{   
    vector<overlap> overlaps {overlap(0,1)};
    
    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 1);
    ASSERT_EQ(chaining[0], overlap(0,1));
}


/** 2 parallel pairs
 */
TEST(ColinearChaining, parallel_overlap)
{   
    vector<overlap> overlaps {overlap(3,1), overlap(7, 12)};
    
    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    cout << "Chaining: " << endl;
    for(auto el: chaining){
        cout << el.first << " " << el.second << endl;
    }

    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], overlaps[0]);
    ASSERT_EQ(chaining[1], overlaps[1]);
}