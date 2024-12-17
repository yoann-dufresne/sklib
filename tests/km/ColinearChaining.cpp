#include <gtest/gtest.h>
#include <utility>

#include <algorithms/SkmerSorting.hpp>

using overlap = std::pair<uint64_t, uint64_t>;
using namespace km::sorting;


/** Single pair test
 */
TEST(ColinearChaining, single_overlap)
{   
    vector<overlap> overlaps {overlap(0,1)};
    
    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 1);
    ASSERT_EQ(chaining[0], overlap(0,1));
}