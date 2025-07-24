#include <vector>
#include <fstream>
#include <stdexcept>
#include <forward_list>
#include <gtest/gtest_prod.h>

#include <algorithms/ColinearChaining.hpp>
#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

#ifndef VIRTUALSKMER_H
#define VIRTUALSKMER_H


// namespace std {
//     template <>
//     struct hash<std::pair<unsigned long, unsigned long>> {
//         std::size_t operator()(const std::pair<unsigned long, unsigned long>& p) const noexcept {
//             return std::hash<unsigned long>{}(p.first) ^ (std::hash<unsigned long>{}(p.second) << 1);
//         }
//     };
// }

namespace km
{
namespace sortedlist
{

// template<typename kuint>
// using kpair = km::Skmer<kuint>::pair;
// using overlap = std::pair<uint64_t, uint64_t>;

// template <class It, typename kuint>
// class compare_kmer_skmer_pos {
//     uint64_t position;
//     SkmerManipulator<kuint> & manipulator;
//     const It start;
//     const It end;

// public:
//     // the comparison function takes as argument 2 integers, a position and the vector of skmers. 
//     // It compares the two skmers in the selected position and returns which one is before the other.
//     compare_kmer_skmer_pos(uint64_t p, SkmerManipulator<kuint> & skmer_manipulator, const It start_skmer_en, const It end_skmer_en) // 
//     : position(p), manipulator(skmer_manipulator), start(start_skmer_en), end(end_skmer_en) {}

//     bool operator()(const uint64_t skmer_id_1,const uint64_t skmer_id_2) const {
//         assert((start+skmer_id_1) < end);
//         assert((start+skmer_id_2) < end);
//         return manipulator.kmer_lt_kmer(*(start+skmer_id_1), position, *(start+skmer_id_2), position);
//     }
// };

// INTEGER USED TO CHECK THAT THE ENDIANESS IS CORRECT
constexpr uint64_t ENDINANESS_SANITY_INTEGER = 0x56534B4D45525F4DULL; // "VSKMER_M" in ASCII

// Helper function to check endiannesss
inline uint64_t swap_endian(uint64_t value) {
    return ((value & 0x00000000000000FFULL) << 56) |
           ((value & 0x000000000000FF00ULL) << 40) |
           ((value & 0x0000000000FF0000ULL) << 24) |
           ((value & 0x00000000FF000000ULL) << 8)  |
           ((value & 0x000000FF00000000ULL) >> 8)  |
           ((value & 0x0000FF0000000000ULL) >> 24) |
           ((value & 0x00FF000000000000ULL) >> 40) |
           ((value & 0xFF00000000000000ULL) >> 56);
}


// PAIR WRAPPER CLASS TO PROCESS PAIR OF COLUMNS
class SlidingWindow {
    std::vector<uint64_t> current;
    std::vector<uint64_t> next;

public:
    void slide(std::vector<uint64_t> new_next) {
        this->current = std::move(this->next);
        this->next = std::move(new_next);
    }    
    const std::vector<uint64_t>& left() const {return current;}
    const std::vector<uint64_t>& right() const {return next;}
};

// VIRTUAL SUPERKMER OBJECT
template<typename kuint>
struct Virtual_skmer {
    Skmer<kuint> skmer;
    uint64_t last_id;

    Virtual_skmer() : skmer(), last_id(0) {}

    Virtual_skmer(km::Skmer<kuint>& passed_skmer, SkmerManipulator<kuint>& m_manip, uint64_t kmer_pos, uint64_t id_value)
        : skmer(m_manip.get_skmer_of_kmer(passed_skmer,kmer_pos)), last_id(id_value) {}

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

    void inline add_kmer(std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint>& m_manip, uint64_t skmer_id, uint64_t kmer_pos){
        using kpair = typename Skmer<kuint>::pair;
        // Free the nucleotide slot in the skmer to accomodate the nucleotide associated with the kmer
        // As it would contain 1s in not used slots
        m_manip.clean_nucleotide_position_skmer(this->skmer, kmer_pos + m_manip.k - 1);
        // Extract the nucleotide from the "contigous" corresponding skmer
        kpair nucleotide {m_manip.extract_nucleotide(skmer_enumeration[skmer_id],kmer_pos + m_manip.k - 1)};
        // Add the nucleotide by OR logical operation
        this->skmer.m_pair.m_value[0] |= nucleotide.m_value[0];
        this->skmer.m_pair.m_value[1] |= nucleotide.m_value[1];
        // HANDLING PREFIX / SUFFIX
        assert(this->skmer.m_suff_size < (2 * m_manip.k - m_manip.m + 1)/2 );
        this->skmer.m_suff_size += 1;
        return;
    }

};

// VIRTUAL SUPERKMER LIST CLASS
// Forward declaring VirtualSkmerSerializer
template<typename kuint> class VirtualSkmerSerializer;

template<typename kuint>
class SortedVirtualSkmerList {
    using LList = std::vector<Virtual_skmer<kuint>>;
    using overlap = std::pair<uint64_t, uint64_t>;

    public:

    // Constructor with manipulator
    SortedVirtualSkmerList(uint64_t k, uint64_t m) : m_manip(k, m) {}

    void print_list() const{
        std::cout << "list : {";
        for(char comma[3] = {'\0', ' ', '\0'}; Skmer<kuint> i : m_skmer_list){
            std::cout << comma << i.m_pair;
            comma[0] = ',';
        }
        std::cout << "}" << std::endl;
    }

    void generate_sorted_list_from_enumeration(std::vector<Skmer<kuint> > const & skmer_enumeration) {
        //initialize the linked lists
        std::cout << "INIZIALIZING VALUES" << std::endl;
        
        // initialize columns ids, sliding window of column ids, vectors to store overlaps
        uint64_t right_column_position {0};
        uint64_t left_column_position {0};
        SlidingWindow window;
        std::vector<overlap> candidate_overlaps;
        std::vector<overlap> valid_overlaps;

        // 0 - sort the column ids based on kmers of the first column
        std::cout << "SLIDING FIRST WINDOW" << std::endl;
        window.slide(sort_column(skmer_enumeration.begin(), skmer_enumeration.end(), right_column_position));
        right_column_position++;

        LList m_vskmer_list;

        m_vskmer_list.reserve(window.right().size());
        for(const uint64_t el: window.right()){
            m_vskmer_list.emplace_back(m_manip.get_skmer_of_kmer(skmer_enumeration[el],0),el);
            std::cerr << "GET_SKMER_OF_KMER HAS PREFIX: " << m_manip.get_skmer_of_kmer(skmer_enumeration[el],0).m_pref_size << " AND SUFFIX: " << m_manip.get_skmer_of_kmer(skmer_enumeration[el],0).m_suff_size << std::endl;
        }

        // while there are columns, compute the next column, compute valid overlaps, merge them into VirtualSkmer
        while(right_column_position <= m_manip.k - m_manip.m ){
            std::cerr << "COL POSITION: " << right_column_position << " OUT OF " << m_manip.k - m_manip.m << std::endl;
            // 1 - sort the column ids based on kmers
            std::cout << "SLIDING WINDOW IN " << right_column_position << " ITERATION." << std::endl;
            window.slide(sort_column(skmer_enumeration.begin(), skmer_enumeration.end(), right_column_position));
            
            std::cout << "LEFT COLUMN:" << std::endl;
            for (const uint64_t el: window.left()){
                std::cout << "L: " << el << ";\t";
            }
            std::cout << std::endl;
            std::cout << "RIGHT COLUMN:" << std::endl;
            for (const uint64_t el: window.right()){
                std::cout << "R: " << el << ";\t";
            }
            std::cout << std::endl;

            // 2 - compute candidate overlaps for a pair of columns
            std::cout << "get_candidate_overlaps. Left column size: " << window.left().size() << "; Right column size: " << window.right().size() << std::endl;
            candidate_overlaps = get_candidate_overlaps(skmer_enumeration, left_column_position, window.left(), window.right());

            // 3 - get valid overlaps using colinear chaining
            std::cout << "colinear_chaining with colinear size of " << candidate_overlaps.size() <<  std::endl;
            for (auto overlap: candidate_overlaps){
                std::cout << "{" << overlap.first << "," << overlap.second << "}" << std::endl;
            }
            
            std::vector<overlap> valid_overlaps;
            if(candidate_overlaps.size() != 0){
                valid_overlaps = km::chaining::colinear_chaining(candidate_overlaps.begin(), candidate_overlaps.end());
            }
            else { valid_overlaps = candidate_overlaps;}
            
            std::cout << "VALID OVERLAPS: {";
            for (auto overlap: valid_overlaps){
                std::cout << "{" << overlap.first << "," << overlap.second << "},";
            }
            std::cout << "}. size: " << valid_overlaps.size() << std::endl;
            // 4 - reconcile kmers by merging columns
            std::cout << "merge_LList_column" << std::endl;
            merge_LList_column(skmer_enumeration, m_vskmer_list, window.right(), valid_overlaps, right_column_position);
            
            // go to next iteration
            left_column_position = right_column_position;
            right_column_position++;
        }
        m_skmer_list.reserve(m_vskmer_list.size());
        for(km::sortedlist::Virtual_skmer<kuint>& vskmer: m_vskmer_list){
            m_skmer_list.emplace_back(std::move(vskmer.skmer));
        }

    }

    std::vector<uint64_t> searchable_positions(
        uint64_t mean, 
        const std::vector<bool> & to_be_searched, 
        const std::vector<std::pair<int64_t,int64_t>> & binary_search_positions) const
    {
        std::vector<uint64_t> searchable;
        searchable.reserve(to_be_searched.size());

        for(uint64_t position {0}; position < to_be_searched.size(); position++){
            if (to_be_searched[position]){
                if (!(mean < binary_search_positions[position].first || mean > binary_search_positions[position].second)){
                    searchable.push_back(position);
                }
            }
        }
        return searchable;
    }

    template<typename T>
    void print_vector(const std::vector<T>& my_v) const {
        std::cout << " {";
        bool first = true;
        for(const T& el: my_v){
            if (!first) std::cout << ", ";
            std::cout << el;
            first = false;
        }
        std::cout << "}\n";
    }

    std::vector<bool> query_skmer(Skmer<kuint> query) const{
        if (m_skmer_list.size() == 0) return std::vector<bool>();
        // 1 CHECK BOUNDARIES SKMER TO EVALUATE WHICH KMERS INSIDE TO QUERY
        auto [query_start_position, query_end_position] = m_manip.get_valid_kmer_bounds(query);

        // PREPARE PARAMETERS FOR SEARCH
        const uint64_t tot_num_kmers_to_search {query_end_position - query_start_position + 1};
        if (tot_num_kmers_to_search <= 0) return std::vector<bool>();
        int mean {0};
        uint64_t current_priority_offset {0};
        uint64_t num_kmers_to_search {tot_num_kmers_to_search};
        std::vector<bool> result(tot_num_kmers_to_search,false); // to be returned
        std::vector<bool> to_search(tot_num_kmers_to_search,true); // to keep track which element I need to keep looking for
        std::vector<std::pair<int64_t,int64_t>> binary_search_boundaries(tot_num_kmers_to_search,{0,m_skmer_list.size()-1}); // to store if the kmer is less than or equal to the one I am looking for now.


        // 2 START BINARY SEARCH
        while(num_kmers_to_search > 0){

            // CASE 1 - I CANNOT FIND IN THE LIST THE KMER OF CURRENT_PRIORIY_OFFSET. I SET IT TO FALSE AND CONTINUE WITH ANOTHER IF POSSIBLE
            if (binary_search_boundaries[current_priority_offset].first > binary_search_boundaries[current_priority_offset].second){
                result[current_priority_offset] = false;
                to_search[current_priority_offset] = false;
                num_kmers_to_search--;
                if (num_kmers_to_search > 0){
                    for (uint64_t i {current_priority_offset}; i < to_search.size(); i++){
                        if (to_search[i]){
                            current_priority_offset = i;
                        }
                    }
                }
                else return result;
            }

            // UPDATE MEAN
            mean = (binary_search_boundaries[current_priority_offset].first + binary_search_boundaries[current_priority_offset].second) / 2;

            // COMPUTE POSITION TO UPDATE FOR BINARY SEARCH
            auto sp = searchable_positions(mean, to_search, binary_search_boundaries);

            // PRINT DEBUG
            std::cout << "MEAN: " << mean << "; NUM_ELEMENTS_TO_SEARCH: " << num_kmers_to_search << std::endl;
            std::cout << "SEARCHABLE POSITIONS: ";
            print_vector(sp);
            std::cout << "RESULT: ";
            print_vector(result);

            // IF NO POSITION FOR BINARY SEARCH I AM DONE, RETURN
            // if (sp.size() == 0){
            //     return result;
            // }

            for(const uint64_t valid_offset: sp){
                uint64_t position {query_start_position + valid_offset};
                if (m_manip.has_valid_kmer(query,position)){

                    // IF THE ELEMENT IS FOUND
                    if (m_manip.kmer_equals_to_kmer(query, m_skmer_list[mean], position)){
                        std::cout << "ELEMENT AT POS " << position << " FOUND." << std::endl;
                        result[valid_offset] = true;
                        to_search[valid_offset] = false;
                        num_kmers_to_search--;
                        // IF IT WAS THE ONE AT PRIORITY, GIVE PRIORITY TO NEXT
                        if(valid_offset == current_priority_offset){
                            for (uint64_t i {valid_offset}; i < to_search.size(); i++){
                                if (to_search[i]){
                                    current_priority_offset = i;
                                }
                            }
                        }
                    }
                    // UPDATING THE MAXIMUM
                    else if (m_manip.kmer_less_than_kmer(query, m_skmer_list[mean], position)){
                        //update values of binary_search_boundaries.second
                        std::cout << "UPDATING UPPER BOUND OF ELEMENT AT binary_search_boundaries OF " << valid_offset << " TO " << mean << std::endl;
                        binary_search_boundaries[valid_offset].second = mean - 1;
                    }
                    // UPDATING THE MINIMUM
                    else{
                        //update values of binary_search_boundaries.first
                        std::cout << "UPDATING LOWER BOUND OF ELEMENT AT binary_search_boundaries OF " << valid_offset << " TO " << mean << std::endl;
                        binary_search_boundaries[valid_offset].first = mean + 1;
                    }
                }
            }
        }
        std::cout << "END" << std::endl;
        // NOW FOR EACH ROUND OF BINARY SEARCH
        // TAKE A SKMER FROM THE SORTED LIST
        // FOR EVERY KMER I NEED TO QUERY IN SEARCH,
        // FOR ALL KMERS THAT HAVE SAME DIRECTION (NEED TO UNDERSTAND HOW TO STATE THIS)
        // CHECK IF IT HAS A VALID SKMER AND THEN TO <
        // UPDATE DIRECTION
        // UPDATE SEARCH
        // UPDATE RESULT
 
        // WHEN NO MORE SKMERS TO QUERY (search is all false)
        // RETURN 
        return result;
    }

    // std::vector<std::vector<uint64_t>> query_enumerated_skmers(std::vector<Skmer<kuint>>){

    // }
    
    void add_list(std::vector<Skmer<kuint>>&  list){
        m_skmer_list = list;
    }

    // Getter for the list 
    const std::vector<Skmer<kuint>>& get_list() const {
        return m_skmer_list;
    }

    // Size getter
    size_t size() const {
        return m_skmer_list.size();
    }

    ~SortedVirtualSkmerList() = default;

    // Move constructor
    SortedVirtualSkmerList(SortedVirtualSkmerList&& other) noexcept = default;

    // Move assignment
    SortedVirtualSkmerList& operator=(SortedVirtualSkmerList&& other) noexcept = default;


    private:
    friend class VirtualSkmerSerializer<kuint>;
    friend class SortedVirtualSkmerListTest;
    
    // KMER SORTING
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer1);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer2);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer3);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer4);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SingleKmerSorting);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SingleKmerSorting2);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SingleKmerSortingReversed);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, ThreeKmerSorting);

    // OVERLAP DETECTION
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlapNoOverlap);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap1Overlap);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap1Left2Right);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left1Right);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightParallel);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed1);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed2);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed3);

    // LISTS MERGING
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, MergeColumnFromEmpty1);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, MergeColumnFromEmpty2);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, MergeLeftRightColumnElements1);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, MergeLeftRightColumnsElements2);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, AddRightColumnElement);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, AddLeftColumnElement);
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, NoElementPointed);

    SkmerManipulator<kuint> m_manip;
    std::vector<Skmer<kuint>> m_skmer_list;           // Final storage, I do not need the uint in the virtual skmer

    /** Sorts skmer ids based on the kmers they contain at a given positon.
     * @param start start_position in the skmer generator
     * @param end end_positon in the skmer generator
     * @param kmer_pos position of the kmer in the skmer (column position)
     * @return a vector of Virtual superkmer ids (if no kmer, no skmer id) 
     **/
    template <class It>
    std::vector<uint64_t> sort_column(It start, It end, uint64_t kmer_pos){   
        // Accessing and comparing kmers in skmers (less than) is done by kmer_lt_kmer of skmermanipulator
        // 1st pass over the column: check which skmers are ok to be processed
        // Check if the first skmer has a kmer in this position
        std::vector<uint64_t> valid_skmer_ids;
        uint64_t sk_id = 0;

        //Iterating over the range [start, end)
        for(It skmer = start; skmer != end; ++skmer)
        {
            if (m_manip.has_valid_kmer(*skmer, kmer_pos)){
                valid_skmer_ids.push_back(sk_id);
            }
            sk_id++;
        }

        // 2nd pass over the column: return ordered list 
        // For every "column" i.e. possible kmer in the skmer size
        // For every skmer that has a kmer in that column
        std::sort(valid_skmer_ids.begin(), valid_skmer_ids.end(), 
            [this, kmer_pos, start](uint64_t id1, uint64_t id2){
                return m_manip.kmer_less_than_kmer(*(start + id1), *(start + id2), kmer_pos);
            });
        std::cout << "SORTED COLUMN: {" << std::endl;
        for (const uint64_t el: valid_skmer_ids){
            std::cout << "{" << el << "," << (start + el)->m_pair << "}\t";
        }
        std::cout << "}" << std::endl;
        return valid_skmer_ids;
    }

    /** Returns candidate overlaps between two columns of Virtual skmer ids
     * @param skmer_enumeration the enumeration of skmer from the iterator
     * @param left_position the "column" position: i.e. the starting point of leftmost kmer considered for the overlap
     * @param left_column the list of skmer_ids that have a valid kmer at the left position
     * @param right_column the list of skmer_ids that have a valid kmer at the left position + 1 (contigous one)
     * @return a vector of pairs of candidate overlaps between the two columns
     **/
    std::vector<overlap> get_candidate_overlaps(std::vector<Skmer<kuint> > const & skmer_enumeration, uint64_t left_position, std::vector<uint64_t> const & left_column, std::vector<uint64_t> const & right_column){
        using kpair = typename Skmer<kuint>::pair;
        using kpairhash = typename Skmer<kuint>::pair_hasher;
        std::unordered_map< kpair, std::vector<uint64_t>, kpairhash > prefixes {};

        kpair suffix, prefix;
        std::vector<std::pair<uint64_t,uint64_t> > candidate_overlaps;
        typename std::unordered_map< kpair, std::vector<uint64_t>, kpairhash >::const_iterator matching_prefix;
        // First, there should be a function that extracts the k-1 prefix of the right column
        for (auto& skmer_id : right_column) {
            // std::cout << "pref" << std::endl;
            assert(skmer_id < skmer_enumeration.size());
            assert(skmer_id >= 0);
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
                    candidate_overlaps.emplace_back(skmer_id,pref_sk_id);
                }
            }
        }
        return candidate_overlaps;
    }

    // Take a pair of Virtual skmer columns (their position), the valid overlaps from the colinear chaining, the skmers in input. It outputs a linked-list of Virtual skmers

    /** Returns the linked list resulting in the merging of the 2 columns
     * @param skmer_enumeration the enumeration of skmer from the iterator
     * @param list the linked list of skmers that have a valid kmer at the left position
     * @param column the list of skmers that have a valid kmer at the left position + 1 (contigous one)
     * @param valid_overlaps the list of overlaps between kmers of the two columns produced by the colinear chaining
     * @param column_pos the position of the column being introduced in the linked_list
     * @return a vector of pairs of candidate overlaps between the two columns
     **/
    void merge_LList_column(std::vector<Skmer<kuint>> const & skmer_enumeration, 
                       std::vector<Virtual_skmer<kuint>>& list,
                       std::vector<uint64_t> const & column, 
                       std::vector<overlap> const & valid_overlaps, 
                       uint64_t const column_pos)
{
    assert(column_pos <= (m_manip.k - m_manip.m));
    
    // Pre-reserve space to avoid reallocations
    list.reserve(list.size() + column.size());
    
    size_t list_idx = 0;
    size_t col_idx = 0;
    size_t overlap_idx = 0;
    
    // Main merge loop
    while (list_idx < list.size() && col_idx < column.size() && 
           overlap_idx < valid_overlaps.size()) {
        
        bool is_left = (list[list_idx].last_id == valid_overlaps[overlap_idx].first);
        bool is_right = (column[col_idx] == valid_overlaps[overlap_idx].second);
        
        if (is_left && is_right) {
            // CASE A: Both elements are in overlap - merge
            list[list_idx].add_kmer(skmer_enumeration, m_manip, column[col_idx], column_pos);
            list[list_idx].last_id = column[col_idx];
            list_idx++;
            col_idx++;
            overlap_idx++;
        }
        else if (is_left) {
            // CASE B: Only list element in overlap - insert from column
            assert(column[col_idx] < skmer_enumeration.size());
            auto skmer = m_manip.get_skmer_of_kmer(skmer_enumeration[column[col_idx]], column_pos);
            list.insert(list.begin() + list_idx,
                Virtual_skmer<kuint>(
                    skmer.m_pair,
                    skmer.m_pref_size,
                    skmer.m_suff_size,
                    column[col_idx]
                )
            );
            list_idx++; // Skip the newly inserted element
            col_idx++;
        }
        else if (is_right) {
            // CASE C: Only column element in overlap - advance list
            list_idx++;
        }
        else {
            // CASE D: Neither in overlap - compare and insert smaller
            assert(column[col_idx] < skmer_enumeration.size());
            auto col_skmer = m_manip.get_skmer_of_kmer(
                skmer_enumeration[column[col_idx]], column_pos);
            
            if (list[list_idx].skmer <= col_skmer) {
                list_idx++;
            } else {
                list.insert(list.begin() + list_idx,
                    Virtual_skmer<kuint>(
                        col_skmer.m_pair,
                        col_skmer.m_pref_size,
                        col_skmer.m_suff_size,
                        column[col_idx]
                    )
                );
                list_idx++;
                col_idx++;
            }
        }
    }
    
    // Handle remaining elements when overlaps are exhausted but both lists have elements
    while (list_idx < list.size() && col_idx < column.size()) {
        assert(column[col_idx] < skmer_enumeration.size());
        auto col_skmer = m_manip.get_skmer_of_kmer(
            skmer_enumeration[column[col_idx]], column_pos);
        
        if (list[list_idx].skmer <= col_skmer) {
            list_idx++;
        } else {
            list.insert(list.begin() + list_idx,
                Virtual_skmer<kuint>(
                    col_skmer.m_pair,
                    col_skmer.m_pref_size,
                    col_skmer.m_suff_size,
                    column[col_idx]
                )
            );
            list_idx++;
            col_idx++;
        }
    }
    
    // Append any remaining column elements
    // This handles the case where list is initially empty or we've reached the end
    while (col_idx < column.size()) {
        assert(column[col_idx] < skmer_enumeration.size());
        auto skmer = m_manip.get_skmer_of_kmer(skmer_enumeration[column[col_idx]], column_pos);
        list.emplace_back(
            skmer.m_pair,
            skmer.m_pref_size,
            skmer.m_suff_size,
            column[col_idx]
        );
        col_idx++;
    }
    
    // If list has remaining elements and column is exhausted, they're already in place - nothing to do
}

}; // end of my class


template<typename kuint>
class VirtualSkmerSerializer {
public:
    static void save(const SortedVirtualSkmerList<kuint>& list, const std::string& filename) {
        std::ofstream outFile(filename, std::ios::binary);
        
        if (outFile.fail()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return outFile.close();
        }

        // ENDIANESS CHECKING INTEGER
        outFile.write(reinterpret_cast<const char*>(&ENDINANESS_SANITY_INTEGER), sizeof(uint64_t));
        
        // Write k and m
        outFile.write(reinterpret_cast<const char*>(&list.m_manip.k), sizeof(uint64_t));
        outFile.write(reinterpret_cast<const char*>(&list.m_manip.m), sizeof(uint64_t));
        
        // Write count
        uint64_t count = list.size();
        outFile.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        
        if (!outFile) {
            std::cerr << "Error writing header to file: " << filename << std::endl;
            return;
        }

        // Calculate chunk parameters
        outFile.write(reinterpret_cast<const char*>(list.m_skmer_list.data()), 
                  count * sizeof(Skmer<kuint>));
        
        int went_good = outFile.good() ? true : false;
        if(!went_good){
            std::cerr << "Error in the writing of the skmer to disk to file: " << filename << std::endl;
        }
        return outFile.close();
    }
    
    /** Loads the sorted virtual skmer list from disk and creates a new manipulator
     * @returns a Sorted_Virtual_Skmer_List object
     * @param filename the path to the binary file to load
     * @throws std::runtime_error if file cannot be opened, is corrupted, or has endianness mismatch
     */
    static SortedVirtualSkmerList<kuint> load(const std::string& filename) {
        std::ifstream inFile(filename, std::ios::binary);
        
        if (inFile.fail()) {
            throw std::runtime_error("Error opening file for reading: " + filename);
        }

        // Read and check read_endianess_int
        uint64_t read_endianess_int;
        inFile.read(reinterpret_cast<char*>(&read_endianess_int), sizeof(uint64_t));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading magic number from file: " + filename);
        }
        
        // Check endianness
        if (read_endianess_int != ENDINANESS_SANITY_INTEGER) {
            uint64_t swapped_endianess_int = swap_endian(read_endianess_int);
            if (swapped_endianess_int == ENDINANESS_SANITY_INTEGER) {
                throw std::runtime_error("Endianness mismatch - file was written on a system with different endianness");
            } else {
                throw std::runtime_error("Invalid file format - ENDINANESS_SANITY_INTEGER mismatch");
            }
        }

        uint64_t file_k, file_m;
        inFile.read(reinterpret_cast<char*>(&file_k), sizeof(uint64_t));
        inFile.read(reinterpret_cast<char*>(&file_m), sizeof(uint64_t));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading k and m values from file: " + filename);
        }

        SortedVirtualSkmerList<kuint> m_virtual_skmer_list(file_k, file_m);

        uint64_t count;
        inFile.read(reinterpret_cast<char*>(&count), sizeof(uint64_t));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading count from file: " + filename);
        }

        // Read the skmer data
        m_virtual_skmer_list.m_skmer_list.reserve(count);
        m_virtual_skmer_list.m_skmer_list.resize(count);
        inFile.read(reinterpret_cast<char*>(m_virtual_skmer_list.m_skmer_list.data()), count * sizeof(Skmer<kuint>));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading virtual skmer data from file: " + filename);
        }
        
        inFile.close();
        return m_virtual_skmer_list;
    }
};


} // namespace sortedlist
} // namespace km

#endif