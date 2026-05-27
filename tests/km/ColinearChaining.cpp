#include <gtest/gtest.h>
#include <utility>
#include <iostream>
#include <cstdint>
#include <set>

#include <algorithms/ColinearChaining.hpp>

using overlap = std::pair<uint64_t, uint64_t>;
using namespace km::chaining;
using namespace std;


// --- Colinear chaining global tests ---

/** Single pair
 */
TEST(ColinearChaining, single_overlap)
{
    vector<overlap> overlaps {overlap(0, 1)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 1);
    ASSERT_EQ(chaining[0], overlap(0, 1));
}


/** 2 parallel pairs
 */
TEST(ColinearChaining, parallel_overlap1)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(7, 12)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);
}

TEST(ColinearChaining, parallel_overlap2)
{
    vector<overlap> overlaps {overlap(0, 0), overlap(3, 3)};
    vector<overlap> expected {overlap(0, 0), overlap(3, 3)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    // for(const overlap& el: chaining ){
    //     std::cout << "(" << el.first << "," << el.second << ")" << std::endl;
    // }

    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);
}

TEST(ColinearChaining, parallel_overlap3)
{
    vector<overlap> overlaps {overlap(1, 1), overlap(3, 3)};
    vector<overlap> expected {overlap(1, 1), overlap(3, 3)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    // for(const overlap& el: chaining ){
    //     std::cout << "(" << el.first << "," << el.second << ")" << std::endl;
    // }
    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);
}

TEST(ColinearChaining, parallel_overlap4)
{
    vector<overlap> overlaps {overlap(4, 1), overlap(5, 5)};
    vector<overlap> expected {overlap(4, 1), overlap(5, 5)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);
}

TEST(ColinearChaining, parallel_overlap5)
{
    vector<overlap> overlaps {overlap(4, 4), overlap(5, 11)};
    vector<overlap> expected {overlap(4, 4), overlap(5, 11)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);
}

TEST(ColinearChaining, parallel_overlap6)
{
    vector<overlap> overlaps {overlap(1, 1), overlap(3, 3), overlap(4, 4), overlap(5, 5)};
    vector<overlap> expected {overlap(1, 1), overlap(3, 3), overlap(4, 4), overlap(5, 5)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    // for(const overlap& el: chaining ){
    //     std::cout << "(" << el.first << "," << el.second << ")" << " -> ";
    // }
    // std::cout << std::endl;

    ASSERT_EQ(chaining.size(), 4);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);
    ASSERT_EQ(chaining[2], expected[2]);
    ASSERT_EQ(chaining[3], expected[3]);
}

TEST(ColinearChaining, parallel_overlap7)
{
    vector<overlap> overlaps {overlap(1, 1)};
    vector<overlap> expected {overlap(1, 1)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), 1);

    ASSERT_EQ(chaining[0], expected[0]);

}

TEST(ColinearChaining, parallel_overlap8)
{
    vector<overlap> overlaps {overlap(0, 0), overlap(1, 1)};
    vector<overlap> expected {overlap(0, 0), overlap(1, 1)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), 2);

    ASSERT_EQ(chaining[0], expected[0]);
    ASSERT_EQ(chaining[1], expected[1]);

}


/** 2 crossing pairs
 */
TEST(ColinearChaining, crossing_overlap)
{
    vector<overlap> overlaps {overlap(7, 1), overlap(3,12)};
    vector<overlap> expected {overlap(7, 1)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), 1);

    ASSERT_EQ(chaining[0], expected[0]);
}


/** V shape with left overlap colision
 */
TEST(ColinearChaining, leftV_overlap)
{
    vector<overlap> overlaps {overlap(7, 1), overlap(7,12)};
    vector<overlap> expected {overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    ASSERT_EQ(chaining[0], expected[0]);
}


/** V shape with right overlap colision
 */
TEST(ColinearChaining, rightV_overlap)
{
    vector<overlap> overlaps {overlap(3, 12), overlap(7,12)};
    vector<overlap> expected {overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    ASSERT_EQ(chaining[0], expected[0]);
}


/** left and right V shape
 */
TEST(ColinearChaining, leftV_rightV_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(3, 12), overlap(7,12)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** parallel left V shape
 */
TEST(ColinearChaining, parallel_leftV_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 7), overlap(5,12)};
    vector<overlap> expected {overlap(3, 1), overlap(5, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** parallel right V shape
 */
TEST(ColinearChaining, parallel_rightV_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 12), overlap(7,12)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}

/** left V parallel shape
 */
TEST(ColinearChaining, leftV_parallel_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(3, 7), overlap(5,12)};
    vector<overlap> expected {overlap(3, 1), overlap(5, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** right V parallel shape
 */
TEST(ColinearChaining, rightV_parallel_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 1), overlap(7,12)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** clique shape
 */
TEST(ColinearChaining, clique_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(7, 1), overlap(3, 12), overlap(7,12)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());

    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}

    // for (uint i{0}; i<chaining.size(); i++)
    // {
    //     std::cout << chaining[i].first << " " << chaining[i].second << std::endl;
    // }

/** 3 parallel pairs
 */
TEST(ColinearChaining, parallel_overlap_3)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 7), overlap(7, 12)};
    vector<overlap> expected {overlap(3, 1), overlap(5, 7), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 3);

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** 2 cross 1 parallel pairs
 */
TEST(ColinearChaining, cross2_parallel1_overlap)
{
    vector<overlap> overlaps {overlap(3, 7), overlap(5, 1), overlap(7, 12)};
    vector<overlap> expected {overlap(3, 7), overlap(7, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 2);

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** 2 cross 1 parallel pairs
 */
TEST(ColinearChaining, parallel1_cross2_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 12), overlap(7, 7)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 7)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 2);

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}

/** 1 overlap 2 cross
 */
TEST(ColinearChaining, doublecross_overlap)
{
    vector<overlap> overlaps {overlap(3, 12), overlap(5, 1), overlap(7, 7)};
    vector<overlap> expected {overlap(5, 1), overlap(7, 7)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 2);

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}

/** 1 overlap 2 cross
 */
TEST(ColinearChaining, doublecross_rev_overlap)
{
    vector<overlap> overlaps {overlap(3, 7), overlap(5, 12), overlap(7, 1)};
    vector<overlap> expected {overlap(3, 7), overlap(5, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 2);

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** triple cross
 */
TEST(ColinearChaining, triple_cross_overlap)
{
    vector<overlap> overlaps {overlap(3, 12), overlap(5, 7), overlap(7, 1)};
    vector<overlap> expected {overlap(7, 1)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** 4 parallel pairs
 */
TEST(ColinearChaining, parallel_overlap_4)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 7), overlap(7, 12), overlap(8, 15)};
    vector<overlap> expected {overlap(3, 1), overlap(5, 7), overlap(7, 12), overlap(8, 15)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), 4);

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** parallel cross parallel pairs
 */
TEST(ColinearChaining, parallel_cross_parallel_overlap)
{
    vector<overlap> overlaps {overlap(3, 1), overlap(5, 12), overlap(7, 7), overlap(8, 15)};
    vector<overlap> expected {overlap(3, 1), overlap(5, 12), overlap(8, 15)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** cross cross pairs
 */
TEST(ColinearChaining, cross_cross_overlap)
{
    vector<overlap> overlaps {overlap(3, 7), overlap(5, 1), overlap(7, 15), overlap(8, 12)};
    vector<overlap> expected {overlap(3, 7), overlap(8, 12)};

    auto chaining = colinear_chaining(overlaps.begin(), overlaps.end());
    ASSERT_EQ(chaining.size(), expected.size());

    for (uint i{0}; i<chaining.size(); i++)
    {
        ASSERT_EQ(chaining[i], expected[i]);
    }
}


/** Regression for issue #6. This candidate set (seen mid-construction on a chr21
 * region, reduced to GGAGCCAAACAGGAGGAAAAGAGG / k=5,m=2) drove colinear_chaining
 * to return overlaps that were never candidates -- e.g. (0,0) and the
 * null_overlap sentinel (2^64-1, 2^64-1) -- which corrupted merge_LList_column
 * and dropped a k-mer. The result must be a strictly-increasing chain made only
 * of input overlaps, of maximum length (here 4: (0,1)(2,4)(4,5)(5,6)).
 */
TEST(ColinearChaining, issue6_no_spurious_overlaps)
{
    vector<overlap> overlaps {overlap(1,0), overlap(0,1), overlap(2,4),
                              overlap(4,5), overlap(3,6), overlap(5,6)};
    std::set<overlap> input(overlaps.begin(), overlaps.end());

    auto chain = colinear_chaining(overlaps.begin(), overlaps.end());

    for (overlap const& o : chain)
        ASSERT_TRUE(input.count(o) == 1) << "spurious overlap (" << o.first << "," << o.second << ")";
    for (size_t i{1}; i<chain.size(); i++)
    {
        ASSERT_LT(chain[i-1].first,  chain[i].first);
        ASSERT_LT(chain[i-1].second, chain[i].second);
    }
    ASSERT_EQ(chain.size(), 4u);
}