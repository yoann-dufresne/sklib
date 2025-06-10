#include <gtest/gtest.h>
#include <utility>
#include <iostream>

#include <algorithms/ColinearChaining.hpp>

using overlap = std::pair<uint64_t, uint64_t>;
using namespace km::chaining;
using namespace std;


// --- Range Max Query tree tests ---


/** Construction of a basic tree with no score update
 */
TEST(RMQtree, construction)
{
    vector<overlap> overlaps {overlap(0,1), overlap(3, 1), overlap(7, 12)};
    RMQtree tree(overlaps.begin(), overlaps.end());

    auto nodes = tree.get_all_nodes();
    ASSERT_EQ(nodes.size(), 7);

    // Keys
    std::vector<overlap> expected_keys {overlap(0,UINT64_MAX), overlap(0,1), overlap(0,UINT64_MAX), overlap(0,1), overlap(3,1), overlap(7,12), overlap(UINT64_MAX, UINT64_MAX)};

    for(uint64_t i=0; i<nodes.size(); i++)
    {
        ASSERT_EQ(nodes[i].key, expected_keys[i]);
        ASSERT_EQ(nodes[i].score, 0);
    }
}


/** Construction of a basic tree with no score update
 */
TEST(RMQtree, basic_score_update)
{
    vector<overlap> overlaps {overlap(0,1), overlap(3, 1), overlap(7, 12)};
    // Expected keys. They should remain the same along all the update processes
    std::vector<overlap> expected_keys {overlap(0,UINT64_MAX), overlap(0,1), overlap(0,UINT64_MAX), overlap(0,1), overlap(3,1), overlap(7,12), overlap(UINT64_MAX, UINT64_MAX)};

    std::vector<uint_fast64_t> expected_scores[3] {
        {7, 7, 0, 7, 0, 0, 0},
        {7, 7, 0, 0, 7, 0, 0},
        {7, 0, 7, 0, 0, 7, 0}
    };

    uint64_t idx = 0;
    for (overlap const& o : overlaps)
    {
        RMQtree tree(overlaps.begin(), overlaps.end());

        tree.update(o, 7);

        std::vector<RMQnode> const& nodes = tree.get_all_nodes();
        for (uint64_t i=0; i<nodes.size(); i++)
        {
            ASSERT_EQ(nodes[i].key, expected_keys[i]);
            ASSERT_EQ(nodes[i].score, expected_scores[idx][i]);
        }

        idx += 1;
    }
}


/** Construction of a simple tree with a score update for each node 
 */
TEST(RMQtree, multiple_score_updates)
{
    vector<overlap> overlaps {overlap(0,1), overlap(3, 1), overlap(7, 12)};
    // Expected keys. They should remain the same along all the update processes
    std::vector<overlap> expected_keys {overlap(0,UINT64_MAX), overlap(0,1), overlap(0,UINT64_MAX), overlap(0,1), overlap(3,1), overlap(7,12), overlap(UINT64_MAX, UINT64_MAX)};

    RMQtree tree(overlaps.begin(), overlaps.end());
    tree.update(overlaps[0], 3);
    tree.update(overlaps[1], 7);
    tree.update(overlaps[2], 4);

    std::vector<uint_fast64_t> expected_scores {7, 7, 4, 3, 7, 4, 0};
    std::vector<RMQnode> const& nodes = tree.get_all_nodes();
    for (uint idx{0}; idx<nodes.size(); idx++)
    {
        ASSERT_EQ(nodes[idx].key, expected_keys[idx]);
        ASSERT_EQ(nodes[idx].score, expected_scores[idx]);
    }
}


/** RMQtree should not be changed by the queries. Some mutaion bugs led to this test */
TEST(RMQtree, rmq_right_stability)
{
    vector<overlap> overlaps {overlap(0,1), overlap(3, 1), overlap(7, 12)};
    // Expected keys. They should remain the same along all the update processes
    std::vector<overlap> expected_keys {overlap(0,UINT64_MAX), overlap(0,1), overlap(0,UINT64_MAX), overlap(0,1), overlap(3,1), overlap(7,12), overlap(UINT64_MAX, UINT64_MAX)};
    std::vector<uint64_t> expected_scores {7, 4, 7, 3, 4, 7, 0};

    RMQtree tree(overlaps.begin(), overlaps.end());
    tree.update(overlaps[0], 3);
    tree.update(overlaps[1], 4);
    tree.update(overlaps[2], 7);

    std::vector<RMQnode> const& nodes = tree.get_all_nodes();
    for (uint idx{0}; idx<nodes.size(); idx++)
    {
        ASSERT_EQ(nodes[idx].key, expected_keys[idx]);
        ASSERT_EQ(nodes[idx].score, expected_scores[idx]);
    }

    ASSERT_EQ(tree.rmq_right(11), 4);

    std::vector<RMQnode> const& after_nodes = tree.get_all_nodes();
    for (uint idx{0}; idx<after_nodes.size(); idx++)
    {
        ASSERT_EQ(after_nodes[idx].key, expected_keys[idx]);
        ASSERT_EQ(after_nodes[idx].score, expected_scores[idx]);
    }
}


/** Test the range query on multiple values before and after the right coordinates from the tree */
TEST(RMQtree, rmq_right)
{
    vector<overlap> overlaps {overlap(0,1), overlap(3, 1), overlap(7, 12)};
    // Expected keys. They should remain the same along all the update processes
    std::vector<overlap> expected_keys {overlap(0,UINT64_MAX), overlap(0,1), overlap(0,UINT64_MAX), overlap(0,1), overlap(3,1), overlap(7,12), overlap(UINT64_MAX, UINT64_MAX)};

    RMQtree tree(overlaps.begin(), overlaps.end());
    tree.update(overlaps[0], 3);
    tree.update(overlaps[1], 4);
    tree.update(overlaps[2], 7);

    ASSERT_EQ(tree.rmq_right(0), 0);
    ASSERT_EQ(tree.rmq_right(1), 4);
    ASSERT_EQ(tree.rmq_right(11), 4);
    ASSERT_EQ(tree.rmq_right(12), 7);
}


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
TEST(ColinearChaining, parallel_overlap)
{   
    vector<overlap> overlaps {overlap(3, 1), overlap(7, 12)};
    vector<overlap> expected {overlap(3, 1), overlap(7, 12)};
    
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