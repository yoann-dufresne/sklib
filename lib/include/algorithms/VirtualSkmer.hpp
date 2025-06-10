#include <vector>

#include "ColinearChaining.hpp"
#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

namespace km
{
namespace sortedlist
{


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

    void add_kmer(std::vector<Skmer<kuint> > const & skmer_enumeration, SkmerManipulator<kuint>& m_manip, uint64_t skmer_id, uint64_t kmer_pos){

        using kpair = typename Skmer<kuint>::pair;
        // Free the nucleotide slot in the skmer to accomodate the nucleotide associated with the kmer
        // As it would contain 1s in not used slots
        m_manip.clean_nucleotide_position_skmer(this->skmer, kmer_pos + m_manip.k - 1);
        // Extract the nucleotide from the "contigous" corresponding skmer
        kpair nucleotide {m_manip.extract_nucleotide(skmer_enumeration[skmer_id],kmer_pos + m_manip.k - 1)};
        // Add the nucleotide by OR logical operation
        // std::cerr << "Pos:" << (kmer_pos + m_manip.k -1) << " , nucleotide: " << nucleotide << " , before: " << virtual_skmer.skmer.m_pair; 
        this->skmer.m_pair.m_value[0] |= nucleotide.m_value[0];
        this->skmer.m_pair.m_value[1] |= nucleotide.m_value[1];
        // std::cerr << " ,after: " << virtual_skmer.skmer.m_pair << " , nucleotide[0]:" << nucleotide << std::endl;
        // HANDLING PREFIX / SUFFIX
        assert(this->skmer.m_suff_size < (2 * m_manip.k - m_manip.m + 1)/2 );
        this->skmer.m_suff_size += 1;
    
        return;
    }

};

// VIRTUAL SUPERKMER LIST CLASS

template<typename kuint>
class Sorted_Virtual_Skmer_List {

    using LList = std::forward_list<Virtual_skmer<kuint>>;
    using overlap = std::pair<uint64_t, uint64_t>;

    std::vector<Virtual_skmer<kuint>> m_list;
    SkmerManipulator<kuint>& m_manip;

    public:

    Sorted_Virtual_Skmer_List(SkmerManipulator<kuint>& manip) : m_manip(manip) {}

    void print_list(){
        std::cout << "list : {";
        for(char comma[3] = {'\0', ' ', '\0'}; Virtual_skmer<kuint> i : m_list){
            std::cout << comma << i.last_id << " : " << i.skmer.m_pair;
            comma[0] = ',';
        }
        std::cout << "}" << std::endl;
    }

    void generate_sorted_list_from_enumeration(std::vector<Skmer<kuint> > const & skmer_enumeration) {
        //initialize the linked list
        LList<kuint> merge_list{};
        
        // initialize columns ids, sliding window of column ids, vectors to store overlaps
        uint64_t right_column_position {0};
        uint64_t left_column_position {0};
        SlidingWindow window;
        std::vector<overlap> candidate_overlaps;
        std::vector<overlap> valid_overlaps;

        // 0 - sort the column ids based on kmers of the first column
        window.slide(sort_column(skmer_enumeration.begin(), skmer_enumeration.end(), right_column_position));
        right_column_position++;

        // while there are columns, compute the next column, compute valid overlaps, merge them into VirtualSkmer
        while(right_column_position < this->m_manip.sk_size ){
            // 1 - sort the column ids based on kmers
            window.slide(sort_column(skmer_enumeration.begin(), skmer_enumeration.end(), right_column_position));

            // 2 - compute candidate overlaps for a pair of columns
            candidate_overlaps = get_candidate_overlaps(skmer_enumeration, left_column_position, window.left(), window.right());

            // 3 - get valid overlaps using colinear chaining
            valid_overlaps = colinear_chaining(candidate_overlaps.begin(), candidate_overlaps.end());

            // 4 - reconcile kmers by merging columns
            merge_LList_column(skmer_enumeration, merge_list, window.right(), valid_overlaps, left_column_position);
            
            // go to next iteration
            left_column_position = right_column_position;
            right_column_position++;
        }

        // FINALLY store the sorted virtual superkmer list as vector
        m_list = std::vector<Virtual_skmer<kuint>>(merge_list.begin(), merge_list.end());
    }

    void load_from_disk(){}

    void dump_to_disk(){}



    private:

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

        // km::SkmerPrettyPrinter<kuint> pp {m_manip.k, m_manip.m};
        //Iterating over the range [start, end)
        for(It skmer = start; skmer != end; ++skmer)
        {
            // pp << *skmer;
            // std::cout << "checking kmer validity" << pp << std::endl;
            if (m_manip.has_valid_kmer(*skmer, kmer_pos)){
                valid_skmer_ids.push_back(sk_id);
                // std::cout << "valid" << std::endl;
            }
            sk_id++;
        }

        // 2nd pass over the column: return ordered list 
        // For every "column" i.e. possible kmer in the skmer size
        // For every skmer that has a kmer in that column
        // std::sort(valid_skmer_ids.begin(), valid_skmer_ids.end(),
                // compare_kmer_skmer_pos<It, kuint>(kmer_pos, this->m_manip, start, end));

        std::sort(valid_skmer_ids.begin(), valid_skmer_ids.end(), 
            [this, kmer_pos, start](uint64_t id1, uint64_t id2){
                return this->m_manip.kmer_lt_kmer(*(start + id1), kmer_pos, *(start + id2), kmer_pos);
            });

        // std::cout << "Virtual SKMER LIST - ( size: " << valid_skmer.size() << ") " << std::endl;
        // for (uint64_t i: valid_skmer) 
        //     std::cout << i << ' ';
        // std::cout << std::endl;
        
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
    void merge_LList_column(std::vector<Skmer<kuint> > const & skmer_enumeration, LList<kuint> & list, std::vector<uint64_t> const & column, std::vector<overlap> const & valid_overlaps, uint64_t const column_pos)
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
                Virtual_skmer<kuint> merging_virtual_skmer = *list_it;
                merging_virtual_skmer.add_kmer(skmer_enumeration, m_manip, *column_it, column_pos);
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
                list.emplace_after(list_it, m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
                *column_it);
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
                if (list_it->skmer <= this->m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it],column_pos)){
                    list_it_previous_element = list_it;
                    ++list_it;
                }
                // If the LLink skmer > enumeration one, do like case (B)
                else{
                    list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
                    list.emplace_after(list_it, m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
                *column_it);
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
            if (list_it->skmer <= this->m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it],column_pos)){
                list_it_previous_element = list_it;
                ++list_it;
            }
            // If the LLink skmer > enumeration one, do like case (B)
            else{
                list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
                list.emplace_after(list_it, m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
                *column_it);
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
            list.emplace_after(list_it, m_manip.get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
                *column_it);
            ++column_it;
            ++list_it;
        }
    }
}

}
}


// COMPARE KMER IN SUPERKMER CLASS
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