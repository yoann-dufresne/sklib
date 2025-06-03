#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <forward_list>
#include <cassert>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

#ifndef SKMERSORTING_H
#define SKMERSORTING_H


namespace km
{
namespace sorting
{
template<typename kuint>
using kpair = km::Skmer<kuint>::pair;
using overlap = std::pair<uint64_t, uint64_t>;

template<typename kuint>
struct Virtual_skmer {
    Skmer<kuint> skmer;
    uint64_t last_id;

    Virtual_skmer() : skmer(), last_id(0) {}

    Virtual_skmer(km::Skmer<kuint>::pair kmer, uint16_t prefix, uint16_t suffix, uint64_t id_value)
        : skmer(kmer, prefix, suffix), last_id(id_value) {}
    Virtual_skmer(km::Skmer<kuint>& passed_skmer, uint64_t id_value)
        : skmer(passed_skmer), last_id(id_value) {}
    Virtual_skmer(km::Skmer<kuint>&& passed_skmer, uint64_t id_value)
        : skmer(passed_skmer), last_id(id_value) {}

    bool operator==(const Virtual_skmer& other) const {
        return (skmer == other.skmer && 
                last_id == other.last_id);
    }

};

template<typename kuint>
using LList = std::forward_list<Virtual_skmer<kuint>>;
template <class It, typename kuint>
class compare_kmer_skmer_pos;

/** Generates a new Virtual Skmer
 * @param skmer_enumeration start_position in the skmer generator
 * @param m_manip skmer manipulator
 * @param skmer_id end_positon in the skmer generator
 * @param kmer_pos position of the kmer in the skmer (column position)
 * @return a new Virtual_skmer 
 **/
template<typename kuint>
Virtual_skmer<kuint> generate_virtual_skmer(std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint>& m_manip, uint64_t skmer_id, uint64_t kmer_pos){
    km::sorting::Virtual_skmer<kuint> s_skmer(m_manip.get_skmer_of_kmer(skmer_enumeration[skmer_id], kmer_pos), skmer_id);
    return s_skmer;
}

/** Add k-mer to virtual super-k-mer
 * @param virtual_skmer the virtual super-k-mer that is being generated
 * @param skmer_enumeration enumeration of superkmers
 * @param m_manip super-k-mer manipulator
 * @param skmer_id the id of the skmer in the skmer_enumeration
 * @param kmer_pos position of the k-mer in the super-k-mer (column position)
 * @return a new virtual super-k-mer
 **/
template<typename kuint>
void add_kmer_virtual_skmer(Virtual_skmer<kuint> & virtual_skmer, std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint>& m_manip, uint64_t skmer_id, uint64_t kmer_pos){

    using kpair = typename Skmer<kuint>::pair;
    // Free the nucleotide slot in the skmer to accomodate the nucleotide associated with the kmer
    // As it would contain 1s in not used slots
    m_manip.clean_nucleotide_position_skmer(virtual_skmer.skmer, kmer_pos + m_manip.k - 1);
    // Extract the nucleotide from the "contigous" corresponding skmer
    kpair nucleotide {m_manip.extract_nucleotide(skmer_enumeration[skmer_id],kmer_pos + m_manip.k - 1)};
    // Add the nucleotide by OR logical operation
    std::cerr << "Pos:" << (kmer_pos + m_manip.k -1) << " , nucleotide: " << nucleotide << " , before: " << virtual_skmer.skmer.m_pair; 
    virtual_skmer.skmer.m_pair.m_value[0] |= nucleotide.m_value[0];
    virtual_skmer.skmer.m_pair.m_value[1] |= nucleotide.m_value[1];
    std::cerr << " ,after: " << virtual_skmer.skmer.m_pair << " , nucleotide[0]:" << nucleotide << std::endl;
    // HANDLING PREFIX / SUFFIX
    assert(virtual_skmer.skmer.m_suff_size < (2 * m_manip.k - m_manip.m + 1)/2 );
    virtual_skmer.skmer.m_suff_size += 1;

    return;
}

template<typename kuint>
void print_list(LList<kuint> const & my_list){
    std::cout << "list : {";
    for(char comma[3] = {'\0', ' ', '\0'}; Virtual_skmer<kuint> i : my_list){
        std::cout << comma << i.last_id << " : " << i.skmer.m_pair;
        comma[0] = ',';
    }
    std::cout << "}" << std::endl;
}

/** Sorts skmer ids based on the kmers they contain at a given positon.
 * @param start start_position in the skmer generator
 * @param end end_positon in the skmer generator
 * @param kmer_pos position of the kmer in the skmer (column position)
 * @param m_manip skmer manipulator
 * @return a vector of Virtual superkmer ids (if no kmer, no skmer id) 
 **/
template<class It, typename kuint>
std::vector<uint64_t> sort_column(It start, It end, uint64_t kmer_pos, SkmerManipulator<kuint>& m_manip)
{   
    // Accessing and comparing kmers in skmers (less than) is done by kmer_lt_kmer of skmermanipulator
    // 1st pass over the column: check which skmers are ok to be processed
    // Check if the first skmer has a kmer in this position
    std::vector<uint64_t> valid_skmer;
    uint64_t sk_id = 0;
    km::SkmerPrettyPrinter<kuint> pp {m_manip.k, m_manip.m};
    //Iterating over the range [start, end)
    for(It skmer = start; skmer != end; ++skmer)
    {
        // pp << *skmer;
        // std::cout << "checking kmer validity" << pp << std::endl;
        if (m_manip.has_valid_kmer(*skmer, kmer_pos)){
            valid_skmer.push_back(sk_id);
            // std::cout << "valid" << std::endl;
        }
        sk_id++;
    }

    // 2nd pass over the column: return ordered list 
    // For every "column" i.e. possible kmer in the skmer size
    // For every skmer that has a kmer in that column
    std::sort(valid_skmer.begin(), valid_skmer.end(),
            compare_kmer_skmer_pos<It, kuint>(kmer_pos, m_manip, start, end));

    // std::cout << "Virtual SKMER LIST - ( size: " << valid_skmer.size() << ") " << std::endl;
    // for (uint64_t i: valid_skmer) 
    //     std::cout << i << ' ';
    // std::cout << std::endl;
    
    return valid_skmer;
}

/** Returns candidate overlaps between two columns of Virtual skmer ids
 * @param skmer_enumeration the enumeration of skmer from the iterator
 * @param m_manip the skmer manipulator inizialized for the iterator
 * @param left_position the "column" position: i.e. the starting point of leftmost kmer considered for the overlap
 * @param left_column the list of skmers that have a valid kmer at the left position
 * @param right_column the list of skmers that have a valid kmer at the left position + 1 (contigous one)
 * @return a vector of pairs of candidate overlaps between the two columns
 **/
template<typename kuint>
std::vector<overlap> get_candidate_overlaps(std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint>& m_manip, uint64_t left_position, std::vector<uint64_t> const & left_column, std::vector<uint64_t> const & right_column)
{
    using kpair = typename Skmer<kuint>::pair;
    using kpairhash = typename Skmer<kuint>::pair_hasher;
    std::unordered_map< kpair, std::vector<uint64_t>, kpairhash > prefixes {};

    kpair suffix, prefix;
    std::vector<std::pair<uint64_t,uint64_t> > candidare_overlaps;
    typename std::unordered_map< kpair, std::vector<uint64_t>, kpairhash >::const_iterator matching_prefix;
    // First, there should be a function that extracts the k-1 prefix of the right column
    for (auto& skmer_id : right_column) {
        // std::cout << "pref" << std::endl;
        assert(skmer_id < skmer_enumeration.size());
        assert(skmer_id > 0);
        prefix = m_manip.extract_prefix_suffix(skmer_enumeration[skmer_id], left_position+1);
        prefixes[prefix].push_back(skmer_id);
    }

    // Second, there should be a function that extracts the k-1 suffix of the left column (same funct as before, just give param the place)
    for (auto& skmer_id : left_column) {
        // std::cout << "suff" << std::endl;
        suffix = m_manip.extract_prefix_suffix(skmer_enumeration[skmer_id], left_position+1);
        
        matching_prefix = prefixes.find (suffix);
        if (matching_prefix != prefixes.end()){
            for (auto& pref_sk_id: matching_prefix->second){
                candidare_overlaps.emplace_back(skmer_id,pref_sk_id);
            }
        }
    }
    return candidare_overlaps;
}

// --- reconciliation of kmers into skmers ---

// 1 - Intialize the linked list with the elements of the first column

/** Returns candidate overlaps between two columns of Virtual skmer ids
 * @param skmer_enumeration the enumeration of skmer from the iterator
 * @param manipulator the skmer manipulator inizialized for the iterator
 * @param column the list of skmers that have a valid kmer at the left position
 * @return a linked list containing the Virtual superkmer and the skmer id of the last kmer added 
 **/
// template<typename kuint>
// LList<kuint> generate_LList(std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint>& m_manip, std::vector<uint64_t> const & column)
// {
//     LList<kuint> Virtual_skmer_llist;
//     Virtual_skmer<kuint> this_skmer; 
//     for (std::vector<uint64_t>::reverse_iterator rit = column.rbegin(); rit != column.rend();  ++rit){
//         // get kmer of that column into a new skmer
//         this_skmer.last_id = *rit;
//         // this_skmer.skmer = manipulator.
//         Virtual_skmer_llist.push_front();
//     }
// }

// 2 - I take a pair of Virtual skmer columns (their position), the valid overlaps from the colinear chaining, and the skmers in input and output a linked-list of Virtual skmers

/** Returns the linked list resulting in the merging of the 2 columns
 * @param skmer_enumeration the enumeration of skmer from the iterator
 * @param manipulator the skmer manipulator inizialized for the iterator
 * @param list the linked list of skmers that have a valid kmer at the left position
 * @param column the list of skmers that have a valid kmer at the left position + 1 (contigous one)
 * @param valid_overlaps the list of overlaps between kmers of the two columns produced by the colinear chaining
 * @param column_pos the position of the column being introduced in the linked_list
 * @return a vector of pairs of candidate overlaps between the two columns
 **/
template<typename kuint>
void merge_LList_column(std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint> & m_manip, LList<kuint> & list, std::vector<uint64_t> const & column, std::vector<overlap> const & valid_overlaps, uint64_t const column_pos)
{   
    // assert(column_pos >= 0);
    assert(column_pos <= (m_manip.k - m_manip.m));

    // using kpair = typename Skmer<kuint>::pair;
    // using kpairhash = typename Skmer<kuint>::pair_hasher;

    // initiliazation of the iterators over the Linked List, Skmer_id column and valid overlaps vector
    // std::cerr << "VARIABLES INITIALIZATION" << std::endl;
    auto list_it_previous_element = list.before_begin(); // this specific iterator element is needed in some cases
    auto list_it = list.begin();
    auto column_it = column.begin();
    auto overlap_it = valid_overlaps.begin();

    // Until one of the two columns is used
    // std::cerr << "BEFORE WHILE LOOP" << std::endl;
    bool verification_list = (list_it != list.end()) ? true : false;
    bool verification_column = (column_it != column.end()) ? true : false;
    // std::cerr << "List_iterator: " << verification_list<< std::endl;
    // std::cerr << "Column_iterator: " << verification_column << std::endl;
    while (list_it != list.end() && column_it != column.end() and overlap_it != valid_overlaps.end() ){
        // std::cerr << "IN WHILE LOOP" << std::endl;

        // 1 - check if elements pointed in the left and right column are in the next valid overlap
        auto curr_value = *list_it;
        bool is_left_in_overlap = (curr_value.last_id == overlap_it->first); //? true : false;
        bool is_right_in_overlap =  (*column_it == (*overlap_it).second);

        // std::cerr << "LL ID: " << curr_value.last_id << " ; COL ID: " << *column_it << std::endl;
        // std::cerr << "OVERLAP: " << overlap_it->first << "," << overlap_it->second << std::endl;
        
        // CASE (A) IF BOTH ELEMENTS ARE POINTED, I MERGE THE VIRTUAL SKMER WITH THE KMER
        if ((is_left_in_overlap == is_right_in_overlap) && (is_left_in_overlap == true)){
            km::sorting::Virtual_skmer merging_virtual_skmer = *list_it;
            add_kmer_virtual_skmer(merging_virtual_skmer, skmer_enumeration,m_manip, *column_it, column_pos);
            merging_virtual_skmer.last_id = *column_it;
            *list_it = merging_virtual_skmer;
            // as we 'used' the element in the linked list, column and overlap list, we go to the next element in all 3;
            list_it_previous_element = list_it;
            ++list_it;
            ++column_it;
            ++overlap_it;
            // std::cerr << "CASE BOTH POINTED" << std::endl;
        }

        // CASE (B) THE ELEMENT IN THE LL IS POINTED, INSERT THE ELEMENT FROM THE COLUMN TO THE LL IN THE PLACE BEFORE
        else if (is_left_in_overlap)
        {
            list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
            list.insert_after(list_it, generate_virtual_skmer(skmer_enumeration,m_manip,*column_it, column_pos));
            ++list_it;
            ++column_it;
            // std::cerr << "CASE LL EL POINTED" << std::endl;
        }

        // CASE (C) THE ELEMENT IN THE COLUMN IS POINTED, I CAN INCREASE THE ITERATOR IN THE LL
        else if (is_right_in_overlap){ // if(is_right_in_overlap)
            list_it_previous_element = list_it;
            ++list_it;
            // std::cerr << "CASE COLUMN EL POINTED" << std::endl;
        }
        
        // CASE (D) BOTH ELEMENTS ARE NOT POINTED, INSERT BASED ON SKMER 
        // NOW COMMENTED AS <=> DEPENDS ON COLUMN POS
        else{
            // If the LLink skmer <= enumeration one, do like case (C)
            if (list_it->skmer <= m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it],column_pos)){
                list_it_previous_element = list_it;
                ++list_it;
            }
            // If the LLink skmer > enumeration one, do like case (B)
            else{
                list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
                list.insert_after(list_it, generate_virtual_skmer(skmer_enumeration,m_manip,*column_it, column_pos));
                ++list_it;
                ++column_it;
            }
        }
    }
    // std::cerr << "AFTER WHILE LOOP" << std::endl;
    // When exiting the while loop, one or both vectors are consumed. I add the final elements by consuming both separately so that if one is not consumed, it will be. This coincides with the special case of filling the LL the first time
    while (list_it != list.end() && column_it != column.end()){
        // std::cerr << "FILLING WHEN OVERLAP IS EMPTY" << std::endl;
        // list.insert_after(list_it, generate_virtual_skmer(skmer_enumeration, m_manip, *column_it, column_pos));
        // ++column_it;
        // ++list_it;
        if (list_it->skmer <= m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it],column_pos)){
            list_it_previous_element = list_it;
            ++list_it;
        }
        // If the LLink skmer > enumeration one, do like case (B)
        else{
            list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
            list.insert_after(list_it, generate_virtual_skmer(skmer_enumeration,m_manip,*column_it, column_pos));
            ++list_it;
            ++column_it;
        }
    }

    // First verify that the iterator over the valid overlaps has exausted eveything
    assert(overlap_it == valid_overlaps.end());

    // If there are still elements in the column, I add them into the linked list.
    list_it = list_it_previous_element; // I need to start to the element before the end.
    while (column_it != column.end()){
        // std::cerr << "NOW FILLING THE LINKED LIST FROM COLUMN" << std::endl;
        list.insert_after(list_it, generate_virtual_skmer(skmer_enumeration, m_manip, *column_it, column_pos));
        ++column_it;
        ++list_it;
    }
}



template <class It, typename kuint>
class compare_kmer_skmer_pos {
    uint64_t position;
    SkmerManipulator<kuint> & manipulator;
    const It start;
    const It end;

public:
    // the comparison function takes as argument 2 integers, a position and the vector of skmers. 
    // It compares the two skmers in the selected position and returns which one is before the other.
    compare_kmer_skmer_pos(uint64_t p, SkmerManipulator<kuint> & skmer_manipulator, const It start_skmer_en, const It end_skmer_en) // 
    : position(p), manipulator(skmer_manipulator), start(start_skmer_en), end(end_skmer_en) {}

    bool operator()(const uint64_t skmer_id_1,const uint64_t skmer_id_2) const {
        assert((start+skmer_id_1) < end);
        assert((start+skmer_id_2) < end);
        return manipulator.kmer_lt_kmer(*(start+skmer_id_1), position, *(start+skmer_id_2), position);
    }
};



// --- Colienar chaining algorithm and data structures---


struct RMQnode
{
    // This pair is used to store an overlap in the leaves and a segment in the internal nodes
    std::pair<uint64_t, uint64_t> key;
    // Max score in the subtree (including the node)
    uint64_t score;

    RMQnode() : key({0, 0}), score(0) {};
};

}// namespace sorting
}// namespace km

namespace std {
    template <>
    struct hash<std::pair<unsigned long, unsigned long>> {
        std::size_t operator()(const std::pair<unsigned long, unsigned long>& p) const noexcept {
            return std::hash<unsigned long>{}(p.first) ^ (std::hash<unsigned long>{}(p.second) << 1);
        }
    };
}

namespace km
{
namespace sorting
{

using overlap = std::pair<uint64_t, uint64_t>;

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

    MaxValueIterator begin(uint64_t score, uint64_t right_boundary) const {
        // std::cout << "begin()" << std::endl;
        return RMQtree::MaxValueIterator(*this, score, right_boundary);
    }
    MaxValueIterator end() const {
        // std::cout << "end()" << std::endl;
        return RMQtree::MaxValueIterator(*this, 0, this->m_num_leaves);
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


} // namespace sorting
} // namespace km

#endif