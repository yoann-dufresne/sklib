#include <sstream>
#include <iostream>

#include "algorithms/ColinearChaining.hpp"

// --- Colienar chaining algorithm and data structures---


namespace km
{
namespace chaining
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

    // Get the number of leaves needed to have a complete binary tree (next power of 2))
    m_num_leaves = 1;
    m_depth = 0;
    while (m_num_leaves < size)
    {
        m_depth++;
        m_num_leaves *= 2;
    }

    // Initialize the tree (m_num_leaves - 1 internal nodes + m_num_leaves leaves)
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

    // Build the parallel tree sorted by key.second (ties broken by key.first).
    // Its internal nodes hold the max key.second of their subtree, which lets
    // rmq_right descend in O(log n) even when the primary tree's leaves have
    // non-monotonic key.second.
    std::vector<overlap> by_second;
    by_second.reserve(m_num_overlaps);
    for (uint64_t i = 0; i < m_num_overlaps; i++)
        by_second.push_back(m_tree[first_leaf_index() + i].key);
    std::sort(by_second.begin(), by_second.end(),
              [](overlap const& a, overlap const& b) {
                  if (a.second != b.second) return a.second < b.second;
                  return a.first < b.first;
              });

    m_second_tree.resize(2 * m_num_leaves - 1);
    for (uint64_t i = 0; i < m_num_leaves; i++)
    {
        if (i < m_num_overlaps)
        {
            m_second_tree[m_num_leaves - 1 + i].key = by_second[i];
            m_second_indexes[by_second[i]] = i;
        }
        else
        {
            m_second_tree[m_num_leaves - 1 + i].key = null_overlap;
        }
    }
    for (int64_t i = m_num_leaves - 2; i >= 0; i--)
    {
        m_second_tree[i].key.second = m_second_tree[2 * i + 2].key.second;
    }
}

/** Update the score of a given overlap of the tree
 * 
 * @param o overlap to update.
 * @param score new score to set.
 */
void RMQtree::update(overlap o, int64_t score)
{
    // Update the primary (key.first-sorted) tree.
    uint64_t node_idx = this->first_leaf_index() + m_indexes[o];
    m_tree[node_idx].score = score;
    while (node_idx > 0)
    {
        node_idx = (node_idx - 1) / 2;
        m_tree[node_idx].score = std::max(m_tree[2 * node_idx + 1].score, m_tree[2 * node_idx + 2].score);
    }

    // Mirror the update into the key.second-sorted tree.
    node_idx = this->first_leaf_index() + m_second_indexes[o];
    m_second_tree[node_idx].score = score;
    while (node_idx > 0)
    {
        node_idx = (node_idx - 1) / 2;
        m_second_tree[node_idx].score = std::max(m_second_tree[2 * node_idx + 1].score, m_second_tree[2 * node_idx + 2].score);
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
    while (idx < this->first_leaf_index())
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
uint64_t RMQtree::rmq_right(uint64_t second_coord) const
{
    if (m_num_overlaps == 0)
        return 0;

    // Descend the key.second-sorted tree to the boundary leaf, then walk back
    // up accumulating the max of every fully-included left subtree. O(log n).
    uint64_t current_idx = 0;
    for (uint64_t lvl{0}; lvl < m_depth; lvl++)
    {
        auto left_idx = 2 * current_idx + 1;
        auto right_idx = 2 * current_idx + 2;

        RMQnode const& left_child = m_second_tree[left_idx];

        if (second_coord < left_child.key.second)
            current_idx = left_idx;
        else
            current_idx = right_idx;
    }

    // Back off one leaf if we overshot the boundary.
    if (second_coord < m_second_tree[current_idx].key.second)
    {
        if (current_idx == m_num_leaves - 1)
            return 0;
        current_idx -= 1;
    }

    uint64_t max_score = m_second_tree[current_idx].score;
    while (current_idx > 0)
    {
        auto parent_idx = (current_idx - 1) / 2;
        bool is_right = (current_idx == 2 * parent_idx + 2);

        if (is_right)
        {
            RMQnode const& left_sibling = m_second_tree[current_idx - 1];
            max_score = std::max(max_score, left_sibling.score);
        }

        current_idx = parent_idx;
    }

    return max_score;
}


// --- Max node Iterator ---

RMQtree::MaxValueIterator::MaxValueIterator(const RMQtree& tree, uint64_t score, uint64_t right_boundary)
    : tree(tree), m_score(score), m_right_boundary(right_boundary), m_leaf_index(UINT64_MAX) {

    if (right_boundary == UINT64_MAX)
        return;

    assert(score <= tree.m_tree[0].score);

    // Try the first leaf; if it doesn't match on score or lies past the
    // right boundary, let next_valid_max() descend the tree to the next one.
    m_leaf_index = 0;
    RMQnode const& first_leaf = tree.m_tree[tree.first_leaf_index()];
    if (first_leaf.score != score || first_leaf.key.second > right_boundary)
        next_valid_max();
}

RMQtree::MaxValueIterator::MaxValueIterator(MaxValueIterator const& other)
    : tree(other.tree), m_score(other.m_score), m_right_boundary(other.m_right_boundary), m_leaf_index(other.m_leaf_index) {
        // std::cout << "Copy constructor" << std::endl;
    }

// --- Définition des opérateurs ---

RMQtree::MaxValueIterator::reference RMQtree::MaxValueIterator::operator*() const {
    // std::cout << "operator*" << std::endl;
    return tree.m_tree[tree.first_leaf_index() + m_leaf_index];
}

RMQtree::MaxValueIterator::pointer RMQtree::MaxValueIterator::operator->() const {
    // std::cout << "operator->" << std::endl;
    return &(tree.m_tree[tree.first_leaf_index() + m_leaf_index]);
}

RMQtree::MaxValueIterator& RMQtree::MaxValueIterator::operator++() {
    // bool is_end = ((*this) == tree.end());
    // std::cout << "operator++ " << m_right_boundary << " " << is_end  << std::endl;
    next_valid_max();
    return *this;
}

RMQtree::MaxValueIterator RMQtree::MaxValueIterator::operator++(int) {
    // std::cout << "before ++(int)" << std::endl;
    MaxValueIterator temp = *this;
    ++(*this);
    return temp;
}

bool RMQtree::MaxValueIterator::operator==(const MaxValueIterator& other) const {
    if (m_right_boundary == UINT64_MAX || other.m_right_boundary == UINT64_MAX)
    {
        // If one of the iterators is at the end, they are equal
        return m_right_boundary == other.m_right_boundary;
    }
    // std::cout << "== operator " << m_leaf_index << " " << other.m_leaf_index << std::endl;
    return m_leaf_index == other.m_leaf_index;
}

std::string RMQtree::toDot() const
{
    std::stringstream ss;

    ss << "digraph RMQtree {" << std::endl;
    for (uint64_t i = 0; i < m_tree.size(); i++)
    {
        RMQnode const& node = m_tree[i];
        ss << i << " [label=\"" << node.key.first << "," << node.key.second << " (" << node.score << ")\"];" << std::endl;
        if (i > 0)
        {
            uint64_t parent = (i - 1) / 2;
            ss << parent << " -> " << i << ";" << std::endl;
        }
    }
    ss << "}" << std::endl;

    return ss.str();
}

// --- Find the next valid max scored overlap ---
void RMQtree::MaxValueIterator::next_valid_max() {
    assert(m_right_boundary != UINT64_MAX);
    uint64_t current_idx = tree.first_leaf_index() + m_leaf_index;

    // 1 - Go up until we find a right sibling whose subtree might contain a
    //     leaf with score >= m_score.
    while (current_idx > 0)
    {
        auto parent_idx = (current_idx - 1) / 2;
        bool is_left = (current_idx == 2 * parent_idx + 1);

        if (is_left)
        {
            RMQnode const& right_child = tree.m_tree[current_idx + 1];
            if (right_child.score >= m_score)
            {
                current_idx = current_idx + 1;
                break;
            }
        }

        current_idx = parent_idx;
    }

    // 1 bis - Reaching the root means there is no further candidate.
    if (current_idx == 0)
    {
        m_right_boundary = UINT64_MAX;
        return;
    }

    // 2 - Descend to the leftmost leaf of that subtree with score >= m_score.
    while (current_idx < tree.m_num_leaves - 1)
    {
        auto left_idx = 2 * current_idx + 1;
        auto right_idx = 2 * current_idx + 2;

        RMQnode const& left_child = tree.m_tree[left_idx];

        if (left_child.score < m_score)
            current_idx = right_idx;
        else
            current_idx = left_idx;
    }

    // 3 - If the leaf fails the score or boundary test, advance past it and
    //     keep searching.  Do NOT conclude "end of iteration" from a single
    //     out-of-range leaf: leaves are sorted by key.first, so later leaves
    //     may still satisfy key.second <= right_boundary.
    RMQnode const& leaf = tree.m_tree[current_idx];
    if (leaf.score != m_score || leaf.key.second > m_right_boundary)
    {
        m_leaf_index = current_idx - tree.first_leaf_index() + 1;
        if (m_leaf_index >= tree.m_num_leaves)
        {
            m_right_boundary = UINT64_MAX;
            return;
        }
        next_valid_max();
        return;
    }

    m_leaf_index = current_idx - tree.first_leaf_index();
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
    
    // std::cout << "[colinear_chaining]::SORTED OVERLAPS" << std::endl;
    // for(auto el = begin; el != end; el++){
    //     std::cout << "(" << el->first << "," << el->second << ")" << " ";
    // }
    // std::cout << std::endl;

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
    std::unordered_map<overlap, overlap> previous_overlaps;
    for (auto it = begin; it != end; it++)
    {
        overlap const& current_overlap = *it;
        // std::cout << "\tPROCESSING OVERLAP: (" << current_overlap.first << "," << current_overlap.second << ")" << std::endl;
        
        // Get previous max scores with compatible second coordinate
        uint64_t max_score = tree.rmq_right(it->second - 1);
        // std::cout << "\t\tMax score: " << max_score << std::endl;
        
        if (max_score == 0)
        {
            // std::cout << "\t\tNo compatible overlap found. (score 0)" << std::endl;
            // No compatible overlap
            tree.update(current_overlap, 1);
            previous_overlaps[current_overlap] = null_overlap;
            continue;
        }
        
        overlap previous {null_overlap};
        // std::cout << "POUET" << std::endl;
        // std::cout << tree.toDot() << std::endl;

        for (auto max_it=tree.begin(max_score, it->second-1); max_it != tree.end(); max_it++)
        {
            // std::cout << "POUET PO" << std::endl;
            RMQnode const& node = *max_it;
            // std::cout << "\t\tChecking overlap: (" << node.key.first << "," << node.key.second << ") with score " << node.score << std::endl;

            // Is it compatible on the first coordinate?
            if (node.key.first < current_overlap.first)
            {
                // std::cout << "\t\tCompatible overlap: (" << node.key.first << "," << node.key.second << ")" << std::endl;
                // Update the score of the current overlap
                max_score = node.score + 1;
                previous = node.key;
                break;
            }
            else
            {
                // std::cout << "\t\tNot compatible overlap: (" << node.key.first << "," << node.key.second << ")" << std::endl;
                previous = previous_overlaps[node.key];
            }
        }
        // std::cout << "POUET POUET" << std::endl;

        // Update the score
        tree.update(current_overlap, max_score);

        previous_overlaps[current_overlap] = previous;

    }

    uint64_t score = tree.max_score();
    // std::cout << "[colinear_chaining]::MAX SCORE: " << score << std::endl;

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
