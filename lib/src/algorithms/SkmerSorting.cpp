

#include "algorithms/SkmerSorting.hpp"

// --- Colienar chaining algorithm and data structures---


namespace km
{
namespace sorting
{

    static overlap const null_overlap = {UINT64_MAX, UINT64_MAX};

/** Constructor of the RMQtree
 * 
 * @param begin iterator to the first element of the list
 * @param end iterator to the last element of the list
 */
RMQtree::RMQtree(std::vector<overlap>::iterator begin, std::vector<overlap>::iterator end)
{
    // Get the size of the list from the iterators
    uint64_t size = std::distance(begin, end);
    m_num_overlaps = size;

    // Get the number of leaves needed to have a complete binary tree
    m_num_leaves = 1;
    m_depth = 0;
    while (m_num_leaves < size)
    {
        m_depth++;
        m_num_leaves *= 2;
    }

    // Initialize the tree
    m_tree.resize(2 * m_num_leaves - 1);
    
    // Fill the leaves with the second value of the overlaps
    auto& current_overlap = begin;
    for (uint64_t i = 0; i < m_num_leaves; i++)
    {
        if (current_overlap != end)
        {
            m_tree[m_num_leaves - 1 + i].key = *current_overlap ;
            m_indexes[*current_overlap] = i;
            current_overlap++;
        }
        else
        {
            m_tree[m_num_leaves - 1 + i].key = null_overlap;
        }
    }

    // Fill the internal nodes with the segment right limit that they represent
    for (int64_t i = m_num_leaves - 2; i >= 0; i--)
    {
        m_tree[i].key.second = m_tree[2 * i + 2].key.second;
    }
}

/** Update the score of a given overlap of the tree
 * 
 * @param o overlap to update.
 * @param score new score to set.
 */
void RMQtree::update(overlap o, int64_t score)
{
    //                 ( internal nodes )
    uint64_t node_idx = m_num_leaves - 1 + m_indexes[o];
    m_tree[node_idx].score = score;

    // Update the parent nodes
    while (node_idx > 0)
    {
        node_idx = (node_idx - 1) / 2;
        m_tree[node_idx].score = std::max(m_tree[2 * node_idx + 1].score, m_tree[2 * node_idx + 2].score);
    }
}

uint64_t RMQtree::max_score() const
{
    return m_tree[0].score;
}

overlap RMQtree::get_max_overlap() const
{
    uint64_t idx = 0;
    uint64_t const score {m_tree[0].score};

    // Go down the tree, following the rightmost max score
    while (idx < m_num_leaves - 1)
    {
        auto left_idx = 2 * idx + 1;
        auto right_idx = 2 * idx + 2;

        RMQnode const& right_child = m_tree[right_idx];

        if (score == right_child.score)
            idx = right_idx;
        else
            idx = left_idx;
    }

    return m_tree[idx].key;
}

/** Get the maximum score in the range [0, second_coord]
 * 
 * @param second_coord Maximal second coord to consider
 * @return The max score in the range [0, second_coord]
 */
uint64_t RMQtree::rmq_right(uint64_t second_coord)
{
    // Go down to the leftmost non-compatible leaf
    uint64_t current_idx = 0;
    for (uint64_t lvl{0} ; lvl<m_depth ; lvl++)
    {
        auto left_idx = 2 * current_idx + 1;
        auto right_idx = 2 * current_idx + 2;

        RMQnode& left_child = m_tree[left_idx];

        if (second_coord < left_child.key.second)
            current_idx = left_idx;
        else
            current_idx = right_idx;
    }
    
    // Take the first compatible leaf
    RMQnode& current_node = m_tree[current_idx];
    if (second_coord < current_node.key.second)
    {
        if (current_idx == m_num_leaves - 1)
            return 0;
        current_idx -= 2;
        current_node = m_tree[current_idx];
    }

    // Go up to the root and get the max score
    uint64_t max_score = 0;
    while (current_idx > 0)
    {
        auto parent_idx = (current_idx - 1) / 2;
        bool is_right = (current_idx == 2 * parent_idx + 2);
        
        if (is_right)
        {
            RMQnode& left_child = m_tree[current_idx - 1];
            max_score = std::max(max_score, left_child.score);
        }

        current_idx = parent_idx;
    }

    return max_score;
}


// --- Max node Iterator ---

RMQtree::MaxValueIterator::MaxValueIterator(const RMQtree& tree, uint64_t score, uint64_t right_boundary)
    : tree(tree), m_score(score), m_right_boundary(right_boundary), m_leaf_index(0) {
    assert(score <= tree.m_tree[0].score);
    // Setup the first leaf
    m_leaf_index = 0;
    uint64_t const vect_idx = tree.m_num_leaves - 1;

    // Find the first leaf with the right score
    if (tree.m_tree[vect_idx].score != score)
        next_valid_max();

    // Verify the score
    assert(tree.m_tree[m_leaf_index + tree.m_num_leaves - 1].score == score);
}

// --- Définition des opérateurs ---

RMQtree::MaxValueIterator::reference RMQtree::MaxValueIterator::operator*() const {
    return tree.m_tree[m_leaf_index + tree.m_num_leaves - 1];
}

RMQtree::MaxValueIterator::pointer RMQtree::MaxValueIterator::operator->() const {
    return &(tree.m_tree[m_leaf_index + tree.m_num_leaves - 1]);
}

RMQtree::MaxValueIterator& RMQtree::MaxValueIterator::operator++() {
    next_valid_max();
    return *this;
}

RMQtree::MaxValueIterator RMQtree::MaxValueIterator::operator++(int) {
    MaxValueIterator temp = *this;
    ++(*this);
    return temp;
}

bool RMQtree::MaxValueIterator::operator==(const MaxValueIterator& other) const {
    return m_leaf_index == other.m_leaf_index;
}

// --- Find the next valid max scored overlap ---
void RMQtree::MaxValueIterator::next_valid_max() {
    // 1 - Go up the tree until we find a right sibling with the same score
    
    // 2 - Go down the tree to the leftmost leaf with a score >= to the current score

    // 3 - If score is not the right one, recursive call
    if (tree.m_tree[m_leaf_index + tree.m_num_leaves - 1].score != m_score)
        next_valid_max();
}


/** Colinear chaining algorithm to select a compatible set of overlaps. The overlaps are compatible if their 
 * coordinates are not crossing and if they do not have common coordinates.
 * WARNING: The algorithme will change the order of the input vector.
 * 
 * @param begin iterator to the first element of the vector
 * @param end iterator to the last element of the vector
 * 
 * @return a vector of overlaps containing compatible overlaps
 **/
std::vector<overlap> colinear_chaining(std::vector<overlap>::iterator begin, std::vector<overlap>::iterator end)
{
    // 1 - Sort the overlaps by the first coordinate.
    std::sort(begin, end, [](const overlap& a, const overlap& b) {
        if (a.first == b.first)
            return a.second < b.second;
        return a.first < b.first;
    });

    // 2 - Create a tree according to the order from 1 and initialize the scores with 0.
    RMQtree tree {begin, end};
    
    // 3 - Sort the overlaps by the second coordinate for the iteration.
    std::sort(begin, end, [](const overlap& a, const overlap& b) {
        if (a.second == b.second)
            return a.first < b.first;
        return a.second < b.second;
    });


    // 4 - For each overlap, update the score according to the best chaining.
    // 4 bis - Also register the previous compatible overlap
    unordered_map<overlap, overlap> previous_overlaps;
    for (auto it = begin; it != end; it++)
    {
        overlap const& current_overlap = *it;

        // Get previous max scores with compatible second coordinate
        uint64_t max_score = tree.rmq_right(it->second - 1);
        
        if (max_score == 0)
        {
            // No compatible overlap
            tree.update(current_overlap, 1);
            previous_overlaps[current_overlap] = null_overlap;
            continue;
        }

        overlap previous {null_overlap};
        for (auto max_it=tree.begin(max_score, it->second-1); max_it != tree.end(); max_it++)
        {
            RMQnode const& node = *max_it;
            // Is it compatible on the first coordinate?
            if (node.key.first < current_overlap.first)
            {
                // Update the score of the current overlap
                max_score = node.score + 1;
                previous = node.key;
                break;
            }
        }

        // Update the score
        tree.update(current_overlap, max_score);
        previous_overlaps[current_overlap] = previous;
    }


    uint64_t score = tree.max_score();
    std::vector<overlap> overlaps(score);
    // 5 - Get the last overlap of the chain
    overlaps[--score] = tree.get_max_overlap();

    // 6 - Get the chain of overlaps
    while (score > 0)
    {
        uint64_t const prev_score = score;
        overlaps[--score] = previous_overlaps[overlaps[prev_score]];
    }

    return overlaps;
}

}} // namespace sorting // namespace km
