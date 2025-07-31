#include <vector>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <stdint.h>

#ifndef COLINEARCHAINING_H
#define COLINEARCHAINING_H

namespace std {
    template <>
    struct hash<std::pair<unsigned long, unsigned long>> {
        std::size_t operator()(const std::pair<unsigned long, unsigned long>& p) const noexcept {
            return std::hash<unsigned long>{}(p.first) ^ (std::hash<unsigned long>{}(p.second) << 1);
        }
    };
}

// --- Colienar chaining algorithm and data structures---
namespace km
{
namespace chaining
{

using overlap = std::pair<uint64_t, uint64_t>;

struct RMQnode
{
    // This pair is used to store an overlap in the leaves and a segment in the internal nodes
    std::pair<uint64_t, uint64_t> key;
    // Max score in the subtree (including the node)
    uint64_t score;

    RMQnode() : key({0, 0}), score(0) {};
};

/** Tree structure to support range maximum queries for ranges [0, x] where x is the second coordinate of an overlap.
 * WARNING: This RMQ tree structure is not for general purpose. It has been optimized for this specific use case and
 * may trigger exceptions on other contexts.
 */
class RMQtree
{
private:
    std::vector<RMQnode> m_tree;
    uint64_t m_num_overlaps;
    uint64_t m_num_leaves;
    uint64_t m_depth;

    std::unordered_map<overlap, uint64_t> m_indexes;

    /** Get the index of the first leaf node inside of the vector representing the tree.
     * @return the index of the first leaf node
     **/
    uint64_t first_leaf_index() const {
        // There is the number of leaves - 1 internal nodes, so the first leaf is at index m_num_leaves - 1
        return m_num_leaves - 1;
    }

public:

    /** Constructor of the RMQtree
     * 
     * @param begin iterator to the first element of the list
     * @param end iterator to the last element of the list
     */
    RMQtree(std::vector<overlap>::iterator begin, std::vector<overlap>::iterator end);

    /** Update the score of a given overlap of the tree
     * 
     * @param o overlap to update.
     * @param score new score to set.
     */
    void update(overlap o, int64_t score);

    /** Get the maximum score in the range [0, second_coord]
     * 
     * @param second_coord Maximal second coord to consider
     * @return The max score in the range [0, second_coord]
     */
    uint64_t rmq_right(uint64_t second_coord) const;

    /** Get the maximum score in the tree 
     * @return The max score in the tree
     */
    uint64_t max_score() const;

    /** Get the overlap with the maximum score in the tree 
     * @return The overlap with the maximum score
     */
    overlap get_max_overlap() const;

    /** Transform the tree into its dot representation
     * 
     * @return a string containing the dot representation of the tree
     */
    std::string toDot() const;

    // --- Iterator for the enumeration of max compatible overlaps ---

    class MaxValueIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = overlap;
        using difference_type = std::ptrdiff_t;
        using pointer = const RMQnode*;
        using reference = const RMQnode&;

        /** Constructor of the MaxValueIterator
         * 
         * @param tree The RMQtree to iterate on
         * @param score The score to consider for the enumeration
         * @param right_boundary Last right value from the leaves to consider
         **/
        MaxValueIterator(const RMQtree& tree, uint64_t score, uint64_t right_boundary);
        MaxValueIterator(MaxValueIterator const& other);

        reference operator*() const;
        pointer operator->() const;

        MaxValueIterator& operator++();   // Pré-incrémentation
        MaxValueIterator operator++(int); // Post-incrémentation

        bool operator==(const MaxValueIterator& other) const;
        bool operator!=(const MaxValueIterator& other) const {
            return !(*this == other);
        }

    private:
        RMQtree const& tree;
        uint64_t m_score;
        uint64_t m_right_boundary;
        // Index of the leaf (not the tree vector index)
        uint64_t m_leaf_index;

        void next_valid_max();
    };

    /** Get the first iterator for the enumeration of max compatible overlaps
     * 
     * @param score The score to consider for the enumeration
     * @param right_boundary First element of the leaves to exclude from the enumeration
     * @return An iterator to the first max compatible overlap
     **/
    MaxValueIterator begin(uint64_t score, uint64_t right_boundary) const {
        
        // return RMQtree::MaxValueIterator(*this, score, idx);
        return RMQtree::MaxValueIterator(*this, score, right_boundary);
    }
    MaxValueIterator end() const {
        // std::cout << "end()" << std::endl;
        return RMQtree::MaxValueIterator(*this, 0, UINT64_MAX);
    }

    /** Get the nodes of the tree. Function mostly used for testing.
     * 
     * @return a vector of RMQnode
     */
    std::vector<RMQnode> const& get_all_nodes() const { return m_tree; }

};

/** Colinear chaining algorithm to select a compatible set of overlaps. The overlaps are compatible if their 
 * coordinates are not crossing and if they do not have common coordinates.
 * WARNING: The algorithme will change the order of the input vector.
 * 
 * @param begin iterator to the first element of the vector
 * @param end iterator to the last element of the vector
 * 
 * @return a vector of overlaps containing compatible overlaps
 **/
std::vector<overlap> colinear_chaining(std::vector<overlap>::iterator begin, std::vector<overlap>::iterator end);


} // namespace chaining
} // namespace km

#endif
