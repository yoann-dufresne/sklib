#include <vector>
#include <fstream>
#include <stdexcept>

#include "SkmerSorting.hpp"
#include "algorithms/ColinearChaining.hpp"
#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

#ifndef VIRTUALSKMER_H
#define VIRTUALSKMER_H

namespace km
{
namespace sortedlist
{
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
        this->skmer.m_pair.m_value[0] |= nucleotide.m_value[0];
        this->skmer.m_pair.m_value[1] |= nucleotide.m_value[1];
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

    public:

    // Constructor with manipulator - makes a copy
    Sorted_Virtual_Skmer_List(const SkmerManipulator<kuint>& manip) 
        : m_manip(std::make_unique<SkmerManipulator<kuint>>(manip)) {}
    
    // Constructor for loading from file
    Sorted_Virtual_Skmer_List(const std::string& filename) {
        load_from_disk(filename);
    }

    void print_list(){
        std::cout << "list : {";
        for(char comma[3] = {'\0', ' ', '\0'}; Virtual_skmer<kuint> i : m_skmer_list){
            std::cout << comma << i.last_id << " : " << i.skmer.m_pair;
            comma[0] = ',';
        }
        std::cout << "}" << std::endl;
    }

    void generate_sorted_list_from_enumeration(std::vector<Skmer<kuint> > const & skmer_enumeration) {
        //initialize the linked lists
        //<kuint> 
        LList merge_list{};
        
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
        while(right_column_position < m_manip->sk_size ){
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
        // FINALLY store the sorted virtual superkmer list as vector of Skmers 
        m_skmer_list.clear();
        m_skmer_list.reserve(std::distance(merge_list.begin(), merge_list.end()));
        for (const auto& vskmer : merge_list) {
            m_skmer_list.push_back(vskmer.skmer);
        }
    }

    int dump_to_disk(const std::string& filename, size_t chunk_size_mb = 100) {
        std::ofstream outFile(filename, std::ios::binary);
        
        if (outFile.fail()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return 1;
        }

        // ENDIANESS CHECKING INTEGER
        outFile.write(reinterpret_cast<const char*>(&ENDINANESS_SANITY_INTEGER), sizeof(uint64_t));
        
        // Write k and m
        outFile.write(reinterpret_cast<const char*>(&m_manip->k), sizeof(uint64_t));
        outFile.write(reinterpret_cast<const char*>(&m_manip->m), sizeof(uint64_t));
        
        // Write count
        uint64_t count = m_skmer_list.size();
        outFile.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        
        if (!outFile) {
            std::cerr << "Error writing header to file: " << filename << std::endl;
            return 1;
        }

        // Calculate chunk parameters
        constexpr uint64_t bytes_per_mb = 1000000;
        uint64_t bytes_per_vskmer = sizeof(Skmer<kuint>);
        uint64_t max_elements_per_chunk = (chunk_size_mb * bytes_per_mb) / bytes_per_vskmer;
        
        // Write data in chunks
        uint64_t elements_written = 0;
        while (elements_written < count) {
            uint64_t elements_to_write = std::min(max_elements_per_chunk, count - elements_written);
            
            outFile.write(reinterpret_cast<const char*>(m_skmer_list.data() + elements_written), 
                         elements_to_write * sizeof(Skmer<kuint>));
            
            if (!outFile) {
                std::cerr << "Error writing chunk to file: " << filename << std::endl;
                return 1;
            }
            
            elements_written += elements_to_write;
        }
        
        outFile.close();
        return 0;
    }

    // Getter for the list 
    const std::vector<Skmer<kuint>>& get_list() const {
        return m_skmer_list;
    }

    // Size getter
    size_t size() const {
        return m_skmer_list.size();
    }


    private:

    std::unique_ptr<SkmerManipulator<kuint>> m_manip;
    std::vector<Skmer<kuint>> m_skmer_list;           // Final storage, I do not need the uint in the virtual skmer

    /** Loads the sorted virtual skmer list from disk and creates a new manipulator
     * @param filename the path to the binary file to load
     * @throws std::runtime_error if file cannot be opened, is corrupted, or has endianness mismatch
     */
    void load_from_disk(const std::string& filename) {
        std::ifstream inFile(filename, std::ios::binary);
        
        if (inFile.fail()) {
            throw std::runtime_error("Error opening file for reading: " + filename);
        }

        // Read and check magic number
        uint64_t read_endianess_int;
        inFile.read(reinterpret_cast<char*>(&read_endianess_int), sizeof(uint64_t));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading magic number from file: " + filename);
        }
        
        // Check endianness
        if (read_endianess_int != ENDINANESS_SANITY_INTEGER) {
            uint64_t swapped_magic = swap_endian(read_endianess_int);
            if (swapped_magic == ENDINANESS_SANITY_INTEGER) {
                throw std::runtime_error("Endianness mismatch - file was written on a system with different endianness");
            } else {
                throw std::runtime_error("Invalid file format - ENDINANESS_SANITY_INTEGER mismatch");
            }
        }

        // Read k and m values
        uint64_t file_k, file_m;
        inFile.read(reinterpret_cast<char*>(&file_k), sizeof(uint64_t));
        inFile.read(reinterpret_cast<char*>(&file_m), sizeof(uint64_t));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading k and m values from file: " + filename);
        }
        
        // Create new manipulator with the k and m values from the file
        m_manip = std::make_unique<SkmerManipulator<kuint>>(file_k, file_m);
        
        // Log the loaded parameters
        std::cerr << "Loaded parameters from file: k=" << file_k << ", m=" << file_m << std::endl;

        // Read count
        uint64_t count;
        inFile.read(reinterpret_cast<char*>(&count), sizeof(uint64_t));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading count from file: " + filename);
        }

        // Read the skmer data
        m_skmer_list.resize(count);
        inFile.read(reinterpret_cast<char*>(m_skmer_list.data()), count * sizeof(Skmer<kuint>));
        
        if (inFile.fail()) {
            throw std::runtime_error("Error reading virtual skmer data from file: " + filename);
        }
        
        inFile.close();
    }

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
            if (m_manip->has_valid_kmer(*skmer, kmer_pos)){
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
                return this->m_manip->kmer_lt_kmer(*(start + id1), kmer_pos, *(start + id2), kmer_pos);
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
            prefix = m_manip->extract_prefix_suffix(skmer_enumeration[skmer_id], left_position+1);
            prefixes[prefix].push_back(skmer_id);
        }

        // Second, there should be a function that extracts the k-1 suffix of the left column (same funct as before, just give param the place)
        for (auto& skmer_id : left_column) {
            // std::cout << "suff" << std::endl;
            suffix = m_manip->extract_prefix_suffix(skmer_enumeration[skmer_id], left_position+1);
            
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
    //<kuint>
    void merge_LList_column(std::vector<Skmer<kuint> > const & skmer_enumeration, LList & list, std::vector<uint64_t> const & column, std::vector<overlap> const & valid_overlaps, uint64_t const column_pos)
    {   
        // assert(column_pos >= 0);
        assert(column_pos <= (m_manip->k - m_manip->m));

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
                merging_virtual_skmer.add_kmer(skmer_enumeration, *m_manip, *column_it, column_pos);
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
                list.emplace_after(list_it, this->m_manip->get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
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
                if (list_it->skmer <= this->m_manip->get_skmer_of_kmer(skmer_enumeration[*column_it],column_pos)){
                    list_it_previous_element = list_it;
                    ++list_it;
                }
                // If the LLink skmer > enumeration one, do like case (B)
                else{
                    list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
                    list.emplace_after(list_it, this->m_manip->get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
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
            if (list_it->skmer <= this->m_manip->get_skmer_of_kmer(skmer_enumeration[*column_it],column_pos)){
                list_it_previous_element = list_it;
                ++list_it;
            }
            // If the LLink skmer > enumeration one, do like case (B)
            else{
                list_it = list_it_previous_element; // I need to place it to the element before the one pointed.
                list.emplace_after(list_it, this->m_manip->get_skmer_of_kmer(skmer_enumeration[*column_it], column_pos), 
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
}; // end of my class

} // namespace sortedlist
} // namespace km

#endif