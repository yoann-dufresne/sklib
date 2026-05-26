#include <vector>
#include <array>
#include <fstream>
#include <stdexcept>
#include <forward_list>
#include <execution>
#include <future>
#include <chrono>


#include <algorithms/ColinearChaining.hpp>
#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

#ifndef VIRTUALSKMER_H
#define VIRTUALSKMER_H

// FRIEND_TEST is provided by gtest (<gtest/gtest_prod.h>). Non-test
// translation units that include this header don't pull in gtest, so
// stub the macro out to a no-op declaration that absorbs the trailing ';'.
#ifndef FRIEND_TEST
#define FRIEND_TEST(test_case_name, test_name) static_assert(true, "")
#endif

namespace km
{
namespace sortedlist
{

namespace util
{
// INTEGER USED TO CHECK THAT THE ENDIANESS IS CORRECT
constexpr uint64_t ENDIANNESS_SANITY_INTEGER = 0x56534B4D45525F4DULL; // "VSKMER_M" in ASCII
constexpr uint64_t MAX_POSSIBLE_KMERS = 64;
// Helper function to check endianness
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

template<typename T>
void print_vector(const std::vector<T>& my_v){
    std::cout << " {";
    bool first = true;
    for(const T& el: my_v){
        if (!first) std::cout << ", ";
        std::cout << el;
        first = false;
    }
    std::cout << "}\n";
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

// HELPER FUNCTIONS FOR SKMER QUERY
inline void setFalse(uint8_t& val) {val = 0; }
inline void setTrue(uint8_t& val) {val = 1; }

/** Updates the vector containing the currently searchable k-mer positions.
 * @param position_to_search position on the virtual s-kmer list that will be searched
 * @param to_be_searched vector of flags (uint8 as bool) of k-mers that have been found yet or not
 * @param binary_search_positions bound for binary search for each k-mer at that position
 * @param searchable_positions the vector of positions to be searched
 * @param tot_kmers_to_search number of total kmers being searched
 * @return a vector of Virtual superkmer ids (if no kmer, no skmer id)
 **/
inline uint64_t update_searchable_positions(
const int64_t position_to_search,
const uint8_t* to_be_searched,
const std::pair<int64_t,int64_t>* binary_search_positions,
uint64_t* searchable_positions,
const uint64_t tot_kmers_to_search)
{
    uint64_t count {0};
    for(uint64_t position {0}; position < tot_kmers_to_search; position++){
        if (to_be_searched[position] != 0 &&
            position_to_search >= binary_search_positions[position].first &&
            position_to_search <= binary_search_positions[position].second) {
            searchable_positions[count++] = position;
        }
    }
    return count;
}

/** Gets the next position of k-mer in the super-k-mer to be prioritized for the binary search. To harness cache-locality, k-mer positions to be searched in the same range as the previous positions are first looked to choose the new value.
 * @param old_position position on the queried s-k-mer that was searched
 * @param searchable_positions vector of positions of k-mer in the s-k-mer that could be searched in the same range as the old_position
 * @param searchable_count num valid elements in searchable_positions array
 * @param to_search vector that flags with 0 elements not to be searched anymore if found or not found.
 * @param tot_kmers_to_search number of total kmers being searched
 * @return a vector of Virtual superkmer ids (if no kmer, no skmer id)
 **/
inline uint64_t get_new_priority_value(
    const uint64_t old_position,
    const uint64_t* searchable_positions,
    const uint64_t searchable_count,
    const uint8_t* to_search,
    const uint64_t tot_kmers_to_search) {
    // FIRST CHECK IN SEARCHABLE POSITIONS
    for(uint64_t i {0}; i < searchable_count; i++){
        const uint64_t position = searchable_positions[i];
        if (position != old_position && to_search[position] != 0) return position;
    }
    // IF NO ELEMENT IN SEARCHABLE POSIIONS IS GOOD, LOOK FOR KMER STILL TO BE SEARCHED.
    for (uint64_t i {0}; i < tot_kmers_to_search; i++){
        if (to_search[i] != 0) return i;
    }
    // IF NO POSITION IS GOOD, RETURN OLD POSITION (IF I GET HERE, THE SEARCH SHOULD BE END)
    return old_position;
}

void print_query_results(const std::vector<std::vector<uint8_t>> & result_vector, std::ostream & os = std::cout){
    std::cerr << "PRINTING RESULT" << std::endl;
    if(result_vector.size() == 0){
        std::cout << "EMPTY RESULT VECTOR!" << std::endl;
        return;
    }
    for (const std::vector<uint8_t> & skmer_result: result_vector){
        if (skmer_result.size() == 0) continue;
        os << (bool)skmer_result[0];
        for (size_t i {1}; i < skmer_result.size(); i++){
            os << "," << (bool)skmer_result[i];
        }
        os << std::endl;
    }
}

}

// VIRTUAL SUPERKMER CLASS
template<typename kuint>
struct Virtual_skmer {
    Skmer<kuint> skmer;
    uint64_t last_id;
    bool expandable;

    Virtual_skmer() : skmer(), last_id(0), expandable(true) {}

    Virtual_skmer(km::Skmer<kuint>& passed_skmer, SkmerManipulator<kuint>& m_manip, uint64_t kmer_pos, uint64_t id_value)
        : skmer(m_manip.get_skmer_of_kmer(passed_skmer,kmer_pos)), last_id(id_value) {}

    Virtual_skmer(km::Skmer<kuint>::pair kmer, uint16_t prefix, uint16_t suffix, uint64_t id_value)
        : skmer(kmer, prefix, suffix), last_id(id_value), expandable(true) {}
    Virtual_skmer(km::Skmer<kuint>& passed_skmer, uint64_t id_value)
        : skmer(passed_skmer), last_id(id_value), expandable(true) {}
    Virtual_skmer(km::Skmer<kuint>&& passed_skmer, uint64_t id_value)
        : skmer(passed_skmer), last_id(id_value), expandable(true) {}

    bool operator==(const Virtual_skmer& other) const {
        return (skmer == other.skmer &&
                last_id == other.last_id);
    }

    friend std::ostream& operator<<(std::ostream& os, const Virtual_skmer<kuint>& vs)
    {
        os << vs.skmer << " last_id:" << vs.last_id;
        return os;
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
template<typename kuint> class VirtualSkmerSerializer;

template<typename kuint>
class SortedVirtualSkmerList {
    using LList = std::vector<Virtual_skmer<kuint>>;
    using overlap = std::pair<uint64_t, uint64_t>;

    public:

    // Constructor with manipulator
    SortedVirtualSkmerList(uint64_t k, uint64_t m) : m_manip(k, m) {}

    void print_list() const{
        km::sortedlist::util::print_vector(m_skmer_list);
    }

    void generate_sorted_list_from_enumeration(std::vector<Skmer<kuint> > const & skmer_enumeration) {
        // initialize columns ids, sliding window of column ids, vectors to store overlaps
        uint64_t right_column_position {0};
        uint64_t left_column_position {0};
        km::sortedlist::util::SlidingWindow window;
        std::vector<overlap> candidate_overlaps;
        std::vector<overlap> valid_overlaps;

        // 0 - sort the column ids based on kmers of the first column
        window.slide(sort_column(skmer_enumeration.begin(), skmer_enumeration.end(), right_column_position));
        // cout << "order col " << 0 << "\t";
        // for (const auto& x : window.right())
        //     cout << x << ",";
        // cout << endl;
        right_column_position++;

        LList vskmer_list;
        // SkmerPrettyPrinter<kuint> pp {5, 2};

        vskmer_list.reserve(window.right().size());
        for(const uint64_t el: window.right()){
            vskmer_list.emplace_back(m_manip.get_skmer_of_kmer(skmer_enumeration[el],0),el);
        }

        // cout << "VSkmer List " << 0 << endl;
        // for (const Virtual_skmer<kuint>& vskm : vskmer_list){
        //         pp << vskm.skmer;
        //         cout << pp << " last_id:" << vskm.last_id << endl;
        //         // cout << vskm << endl;
        //     }
        // cout << endl << endl;

        // while there are columns, compute the next column, compute valid overlaps, merge them into VirtualSkmer
        while(right_column_position <= m_manip.k - m_manip.m ){
            // 1 - sort the column ids based on kmers
            window.slide(sort_column(skmer_enumeration.begin(), skmer_enumeration.end(), right_column_position));
            // cout << "order col " << right_column_position << "\t";
            // for (const auto& x : window.right())
            //     cout << x << ",";
            // cout << endl;

            // 2 - compute candidate overlaps for a pair of columns
            candidate_overlaps = get_candidate_overlaps(skmer_enumeration, left_column_position, window.left(), window.right());
            // cout << "candidate overlaps for " << right_column_position << "\t";
            // for (const auto& o : candidate_overlaps)
            //     cout << "(" << o.first << ";" << o.second << ") , ";
            // cout << endl;

            // 3 - get valid overlaps using colinear chaining
            std::vector<overlap> valid_overlaps;
            if(candidate_overlaps.size() != 0){
                valid_overlaps = km::chaining::colinear_chaining(candidate_overlaps.begin(), candidate_overlaps.end());
            }
            else { valid_overlaps = candidate_overlaps;}
            // cout << "valid overlaps for " << right_column_position << "\t";
            // for (const auto& o : valid_overlaps)
            //     cout << "(" << o.first << ";" << o.second << ") , ";
            // cout << endl;


            // 4 - reconcile kmers by merging columns
            merge_LList_column(skmer_enumeration, vskmer_list, window.left(), window.right(), valid_overlaps, right_column_position);
            // cout << "VSkmer List " << right_column_position << endl;
            // for (const Virtual_skmer<kuint>& vskm : vskmer_list) {
            //     pp << vskm.skmer;
            //     cout << pp << " last_id:" << vskm.last_id << " expandable: " << vskm.expandable << endl;
            //     // cout << vskm << endl;
            // }
            // cout << endl;

            // go to next iteration
            left_column_position = right_column_position;
            right_column_position++;
            
            // cout << endl;
        }

        m_skmer_list.reserve(vskmer_list.size());

        for(km::sortedlist::Virtual_skmer<kuint>& vskmer: vskmer_list){
            m_skmer_list.emplace_back(std::move(vskmer.skmer));
        }

    }

    std::vector<uint8_t> query_skmer(const Skmer<kuint> query) const{
        // 1 CHECK BOUNDARIES SKMER TO EVALUATE WHICH KMERS INSIDE TO QUERY
        // std::cerr << "I AM IN" << std::endl;

        auto [query_start_position, query_end_position] = m_manip.get_valid_kmer_bounds(query);
        km::SkmerPrettyPrinter<kuint> pp {m_manip.k, m_manip.m};
        pp << query;
        // std::cerr << pp << " ";
        // std::cerr << "skmer properties[ PREFIX:" << query.m_pref_size << "; SUFFIX: " << query.m_suff_size << "; PAIR: " << query.m_pair << std::endl;
        // std::cerr << "query_start_position: " << query_start_position << "; query_end_position: " << query_end_position << std::endl;

        if(query_end_position < query_start_position) {
            // std::cerr << "RETURNING EMPTY" << std::endl;
            return std::vector<uint8_t>(0,0);
        }
        const uint64_t tot_num_kmers_to_search {query_end_position - query_start_position + 1};

        // RETURN ON EDGE CASES (SKMER DOES NOT CONTAIN A KMER OR LIST IS EMPTY)
        if (tot_num_kmers_to_search <= 0 || m_skmer_list.empty()){
            return std::vector<uint8_t>(std::max(0UL, tot_num_kmers_to_search), 0);
        }

        // PREPARE PARAMETERS FOR SEARCH
        uint64_t current_priority_offset {0};
        uint64_t searchable_position_count {tot_num_kmers_to_search};
        uint64_t num_kmers_to_search {tot_num_kmers_to_search};
        uint64_t current_searched_position_in_skmer;
        int64_t mean {0};

        // USING SMALL STACK ALLOCATED ARRAYS FOR FAST QUERY
        uint8_t result[km::sortedlist::util::MAX_POSSIBLE_KMERS] = {0};//'bool' to be returned
        uint8_t keep_searching[km::sortedlist::util::MAX_POSSIBLE_KMERS]; // keep track of which k-mers have been searched yet
        std::pair<int64_t,int64_t> binary_search_boundaries[km::sortedlist::util::MAX_POSSIBLE_KMERS];
        uint64_t positions_to_search[km::sortedlist::util::MAX_POSSIBLE_KMERS];
        // std::cerr << "INITIALZING ARRAYS" << std::endl;
        // std::cerr << "LIST SIZE: " << m_skmer_list.size() << "; tot_num_kmers_to_search: " << tot_num_kmers_to_search << std::endl;
        // ARRAYS INITIALIZATION
        const int64_t list_size = m_skmer_list.size();
        std::fill_n(keep_searching, tot_num_kmers_to_search, 1);
        std::fill_n(binary_search_boundaries, tot_num_kmers_to_search,
            std::make_pair(0LL, static_cast<int64_t>(list_size - 1)));
        for(size_t i {0}; i < tot_num_kmers_to_search; i++){
            positions_to_search[i] = i;
        }
        // std::cout << "STARTIG BYNARY SEARCH" << std::endl;
        // 2 START BINARY SEARCH
        while(num_kmers_to_search > 0){
            // VERIFY THAT THE CURRENT PRIORITY_POSITION IS STILL VALID AND UPDATE
            if (keep_searching[current_priority_offset] == 0) {
                current_priority_offset = km::sortedlist::util::get_new_priority_value(current_priority_offset, positions_to_search, searchable_position_count, keep_searching, tot_num_kmers_to_search);
            }

            // UPDATE MEAN
            const int64_t old_mean {mean};
            mean = (binary_search_boundaries[current_priority_offset].first + binary_search_boundaries[current_priority_offset].second) >> 1;
            if (mean == old_mean){
                int64_t new_pos = find_closest_valid_skmer(mean, binary_search_boundaries[current_priority_offset].first, binary_search_boundaries[current_priority_offset].second, query_start_position + current_priority_offset);
                // std::cout << "new_pos: " << new_pos << std::endl;
                if (new_pos < 0){
                    // there are no positions left, flag as not found and continue
                    km::sortedlist::util::setFalse(keep_searching[current_priority_offset]);
                    num_kmers_to_search--;
                    continue;
                }
                mean = (uint64_t)new_pos;
            }

            // COMPUTE POSITION TO UPDATE FOR BINARY SEARCH
            searchable_position_count = km::sortedlist::util::update_searchable_positions(mean, keep_searching, binary_search_boundaries, positions_to_search, tot_num_kmers_to_search);

            auto [queried_start_position, queried_end_position] = m_manip.get_valid_kmer_bounds(m_skmer_list[mean]);
            // std::cout << "searchable_position_count: " << searchable_position_count << "; num_kmers_to_search: " << num_kmers_to_search << std::endl;
            // std::cout << "mean: " << mean << "; old_mean: " << old_mean << std::endl;
            // std::cout << "current_priority_offset: " << current_priority_offset << std::endl;
            // std::cout << "QUERIED KMER START POSITION: " << queried_start_position << "; QUERIED KMER END POSITION: " << queried_end_position << std::endl;
            // for(size_t i {0}; i < tot_num_kmers_to_search; i++){
            //     std::cout << "[" << binary_search_boundaries[i].first << "," << binary_search_boundaries[i].second << "] - ";
            // }
            // std::cout << std::endl;

            for(uint64_t i {0}; i < searchable_position_count; i++){
                const uint64_t offset {positions_to_search[i]};
                current_searched_position_in_skmer = query_start_position + offset;
                // std::cout << "offset: " << offset << "; current_searched_position_in_skmer: " << current_searched_position_in_skmer << std::endl;
                if(current_searched_position_in_skmer >= queried_start_position && current_searched_position_in_skmer <= queried_end_position){

                    const int kmer_comparison {m_manip.kmer_compare(query, m_skmer_list[mean], current_searched_position_in_skmer)};
                    // std::cout << "kmer comparion: " << kmer_comparison << std::endl;
                    if (kmer_comparison == 0){
                        // FOUND. SET RESULT TO TRUE, KEEP_SEARCHING TO FALSE, UPDATE PRIORITY POSITION IF NECESSARY
                        km::sortedlist::util::setTrue(result[offset]);
                        km::sortedlist::util::setFalse(keep_searching[offset]);
                        num_kmers_to_search--;
                        if (current_searched_position_in_skmer == current_priority_offset){
                            current_priority_offset = km::sortedlist::util::get_new_priority_value(current_priority_offset, positions_to_search, searchable_position_count, keep_searching, tot_num_kmers_to_search);
                        }
                        continue;
                    }
                    else if (kmer_comparison < 0){
                        // UPDATE BINARY SEARCH UPPER BOUNDARY FOR THIS POSITION
                        binary_search_boundaries[offset].second = mean - 1;
                    }
                    else{
                        // UPDATE BINARY SEARCH LOWER BOUNDARY FOR THIS POSITION
                        binary_search_boundaries[offset].first = mean + 1;
                    }
                    if (binary_search_boundaries[offset].first > binary_search_boundaries[offset].second) {
                        km::sortedlist::util::setFalse(keep_searching[offset]);
                        num_kmers_to_search--;
                    }
                }
            }
        }

        return std::vector<uint8_t>(result, result + tot_num_kmers_to_search);
    }

    std::vector<std::vector<uint8_t>> query_skmer_batch(std::vector<Skmer<kuint>> query_skmers) const{
        const size_t size_query {query_skmers.size()};
        std::vector<std::vector<uint8_t>> result(size_query);
        auto start = std::chrono::high_resolution_clock::now();
        std::transform(std::execution::par, query_skmers.begin(), query_skmers.end(), result.begin(), [this](const Skmer<kuint>&  queried_skmer){ return this->query_skmer(queried_skmer); });
        auto end = std::chrono::high_resolution_clock::now();
        auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        // std::cerr << "[BATCHED QUERIES]:: Processed " << size_query << " elements in " << duration_s << " (" <<  duration_ns/size_query << " per element)" << std::endl;
        return result;
    }

    void query(const std::string filename, std::ostream& os = std::cout) {
        constexpr uint64_t MAX_INGESTED_SKMER {4096};
        //start enumeration from sequence
        km::FileSkmerator<kuint> file_skmerator {m_manip, filename};
        typename FileSkmerator<kuint>::Iterator it = file_skmerator.begin();

        //Enumerateing the superkmers from the file
        std::vector<km::Skmer<kuint>> skmer_bufferA;
        std::vector<km::Skmer<kuint>> skmer_bufferB;
        skmer_bufferA.reserve(MAX_INGESTED_SKMER);
        skmer_bufferB.reserve(MAX_INGESTED_SKMER);

        std::vector<km::Skmer<kuint>>* cur = &skmer_bufferA; // being filled
        std::vector<km::Skmer<kuint>>* work = &skmer_bufferB; //to be processed

        std::future<std::vector<std::vector<uint8_t>>> curr_task;

        km::Skmer<kuint> current_skmer;
        // std::cerr << "END OF INITIALIZATION" << std::endl;
        for (km::Skmer<kuint> const skmer : file_skmerator){
          cur->emplace_back(skmer);
          //when 1000 skmers have been loaded, dispatch query_skmer_batch
          if (cur->size() == MAX_INGESTED_SKMER){

            // if a process was already executing, expect it ends and dispatch results to outstream
            if (curr_task.valid()){
              km::sortedlist::util::print_query_results(curr_task.get(), os);
            }
            // swap the two buffers pointers
            std::swap(cur, work);

            // clear current for next calculation
            cur->clear();

            // dispatch batched query thread
            // std::cerr << "DISPATCHING" << std::endl;
            curr_task = std::async(std::launch::async,
                                    &SortedVirtualSkmerList<kuint>::query_skmer_batch, // member fn
                                    this,                         // object on which to call it
                                    std::move(*work));
          }
        }

        // taking care of the last elements in current, if present
        if (!cur->empty()) {
          if (curr_task.valid()){
            km::sortedlist::util::print_query_results(curr_task.get(), os);
          }
        //   std::cerr << "FINAL DISPATCHING" << std::endl;
          auto last = query_skmer_batch(*cur);
        //   std::cerr << "LAST RESULT VECTOR SIZE: " << last.size() << std::endl;
          km::sortedlist::util::print_query_results(last, os);
        }
        else if (curr_task.valid()) {
          km::sortedlist::util::print_query_results(curr_task.get(), os);
        }
    }

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

    uint64_t k() const { return m_manip.k; }
    uint64_t m() const { return m_manip.m; }

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
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, SingleKmerSortingAndUnique);
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

    // QUERY - BINARY SEARCH HELPERS
    FRIEND_TEST(SortedVirtualSkmerListPrivateTest, FindClosestValidSkmerUnderflow);

    SkmerManipulator<kuint> m_manip;
    std::vector<Skmer<kuint>> m_skmer_list;           // Final storage, I do not need the uint in the virtual skmer

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
                return m_manip.kmer_compare(*(start + id1), *(start + id2), kmer_pos) < 0;
            });

        // 3rd pass, remove k-mer duplicates
        auto last_value = std::unique(valid_skmer_ids.begin(), valid_skmer_ids.end(),
            [this, kmer_pos, start](uint64_t id1, uint64_t id2){
                return m_manip.kmer_compare(*(start + id1), *(start+id2), kmer_pos) == 0;
            });

        valid_skmer_ids.erase(last_value, valid_skmer_ids.end());

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
        uint64_t prefix_pos {0};
        for (auto& skmer_id : right_column) {
            // std::cout << "pref" << std::endl;
            assert(skmer_id < skmer_enumeration.size());
            assert(skmer_id >= 0);
            prefix = m_manip.extract_prefix_suffix(skmer_enumeration[skmer_id], left_position+1);
            prefixes[prefix].push_back(prefix_pos);
            prefix_pos++;
        }

        // Second, there should be a function that extracts the k-1 suffix of the left column (same funct as before, just give param the place)
        uint64_t suffix_pos {0};
        for (auto& skmer_id : left_column) {
            // std::cout << "suff" << std::endl;
            suffix = m_manip.extract_prefix_suffix(skmer_enumeration[skmer_id], left_position+1);

            matching_prefix = prefixes.find(suffix);
            if (matching_prefix != prefixes.end()){
                for (auto& pref_sk_id: matching_prefix->second){
                    candidate_overlaps.emplace_back(suffix_pos,pref_sk_id);
                }
            }
            suffix_pos++;
        }
        return candidate_overlaps;
    }

    int64_t find_closest_valid_skmer(const uint64_t position_in_list, const uint64_t minimum, const uint64_t maximum, const uint64_t kmer_position_in_skmer) const{
        // Scan down to `minimum` inclusive. `i` is unsigned, so guard the loop
        // with an explicit `i == minimum` break: a `i >= minimum` condition is
        // always true once minimum == 0 and `i--` would wrap past 0 to 2^64-1,
        // reading m_skmer_list out of bounds.
        for (uint64_t i {position_in_list}; ; i--){
            if (m_manip.has_valid_kmer(m_skmer_list[i], kmer_position_in_skmer)) return (int64_t)i;
            if (i == minimum) break;
        }
        for (uint64_t i {position_in_list}; i <= maximum; i++){
            if (m_manip.has_valid_kmer(m_skmer_list[i], kmer_position_in_skmer)) return (int64_t)i;
        }
        return -1UL;
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
                       std::vector<uint64_t> const & left_column, // BUG here ? Remaining of old version ?
                       std::vector<uint64_t> const & column,
                       std::vector<overlap> const & valid_overlaps,
                       uint64_t const column_pos)
{
    assert(column_pos <= (m_manip.k - m_manip.m));

    // Merge into a fresh buffer to avoid O(n) shifts from mid-vector inserts.
    LList merged;
    merged.reserve(list.size() + column.size());

    size_t list_idx = 0;
    size_t col_idx = 0;
    size_t overlap_idx = 0;

    // Main merge loop
    while (list_idx < list.size() && col_idx < column.size() &&
           overlap_idx < valid_overlaps.size()) {

        bool is_left = list[list_idx].expandable and (list[list_idx].last_id == left_column[valid_overlaps[overlap_idx].first] );
        bool is_right = (column[col_idx] == column[valid_overlaps[overlap_idx].second] );

        if (is_left && is_right) {
            // CASE A: Both elements are in overlap - merge
            list[list_idx].add_kmer(skmer_enumeration, m_manip, column[col_idx], column_pos);
            list[list_idx].last_id = column[col_idx];
            merged.push_back(std::move(list[list_idx]));
            list_idx++;
            col_idx++;
            overlap_idx++;
        }
        else if (is_left) {
            // CASE B: Only list element in overlap - emit column element,
            // keep list element pending against the same overlap entry.
            assert(column[col_idx] < skmer_enumeration.size());
            auto skmer = m_manip.get_skmer_of_kmer(skmer_enumeration[column[col_idx]], column_pos);
            merged.emplace_back(
                skmer.m_pair,
                skmer.m_pref_size,
                skmer.m_suff_size,
                column[col_idx]
            );
            col_idx++;
        }
        else if (is_right) {
            // CASE C: Only column element in overlap - advance list
            // Close the skmer
            list[list_idx].expandable = false;
            merged.push_back(std::move(list[list_idx]));
            list_idx++;
        }
        else {
            // CASE D: Neither in overlap - compare and emit smaller
            assert(column[col_idx] < skmer_enumeration.size());
            auto col_skmer = m_manip.get_skmer_of_kmer(
                skmer_enumeration[column[col_idx]], column_pos);

            if (list[list_idx].skmer <= col_skmer) {
                // Close the skmer
                list[list_idx].expandable = false;
                merged.push_back(std::move(list[list_idx]));
                list_idx++;
            } else {
                merged.emplace_back(
                    col_skmer.m_pair,
                    col_skmer.m_pref_size,
                    col_skmer.m_suff_size,
                    column[col_idx]
                );
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
            merged.push_back(std::move(list[list_idx]));
            list_idx++;
        } else {
            merged.emplace_back(
                col_skmer.m_pair,
                col_skmer.m_pref_size,
                col_skmer.m_suff_size,
                column[col_idx]
            );
            col_idx++;
        }
    }

    // Append any remaining column elements
    while (col_idx < column.size()) {
        assert(column[col_idx] < skmer_enumeration.size());
        auto skmer = m_manip.get_skmer_of_kmer(skmer_enumeration[column[col_idx]], column_pos);
        merged.emplace_back(
            skmer.m_pair,
            skmer.m_pref_size,
            skmer.m_suff_size,
            column[col_idx]
        );
        col_idx++;
    }

    // Append any remaining list elements
    while (list_idx < list.size()) {
        merged.push_back(std::move(list[list_idx]));
        list_idx++;
    }

    list = std::move(merged);
}

}; // end of my class


// VIRTUAL SUPERKMER LIST SERIALIZER
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
        outFile.write(reinterpret_cast<const char*>(&km::sortedlist::util::ENDIANNESS_SANITY_INTEGER), sizeof(uint64_t));

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

    static void save_ascii(const SortedVirtualSkmerList<kuint>& list, const std::string& filename) {
        std::ofstream outFile(filename);

        if (outFile.fail()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return outFile.close();
        }

        const uint64_t k = list.m_manip.k;
        const uint64_t m = list.m_manip.m;
        const uint64_t sk_size = 2 * k - m;
        const uint64_t flank = k - m;
        const uint64_t buf_pref = (sk_size + 1) / 2;
        const uint64_t buf_suff = sk_size / 2;
        uint64_t count = list.size();
        outFile << k << " " << m << " " << count << "\n";

        if (!outFile) {
            std::cerr << "Error writing header to file: " << filename << std::endl;
            return;
        }

        static const char nucleotides[] = {'A', 'C', 'T', 'G'};

        // The interleaved layout stores position p of the superkmer as follows:
        //   p in [0, buf_pref)   -> prefix slot p        at bits [4p+1 : 4p]
        //   p in [buf_pref, sk_size) -> suffix slot (sk_size-1-p) at bits [4s+3 : 4s+2]
        // The minimizer positions [flank, flank+m-1] are split across prefix and
        // suffix slots (not a contiguous bit block), so we iterate slots and the
        // minimizer nucleotides are interleaved along with the flanks.
        for (uint64_t i = 0; i < count; i++) {
            const auto& sk = list.m_skmer_list[i];

            // Prefix side: positions [flank - m_pref_size, buf_pref - 1]
            // (valid prefix flank, then prefix-side minimizer nucleotides)
            for (uint64_t slot = flank - sk.m_pref_size; slot < buf_pref; slot++)
                outFile << nucleotides[(sk.m_pair >> (4 * slot)) & 0b11UL];

            // Suffix side: positions (buf_pref, sk_size - 1 - (flank - m_suff_size)]
            // Iterate slot from (buf_suff - 1) down to (flank - m_suff_size).
            for (int64_t slot = static_cast<int64_t>(buf_suff) - 1;
                 slot >= static_cast<int64_t>(flank) - static_cast<int64_t>(sk.m_suff_size);
                 slot--)
                outFile << nucleotides[(sk.m_pair >> (4 * slot + 2)) & 0b11UL];

            outFile << " " << sk.m_pref_size << " " << sk.m_suff_size;
            outFile << "\n";
        }

        if (!outFile.good())
            std::cerr << "Error writing skmers to file: " << filename << std::endl;
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
        if (read_endianess_int != km::sortedlist::util::ENDIANNESS_SANITY_INTEGER) {
            uint64_t swapped_endianess_int = km::sortedlist::util::swap_endian(read_endianess_int);
            if (swapped_endianess_int == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER) {
                throw std::runtime_error("Endianness mismatch - file was written on a system with different endianness");
            } else {
                throw std::runtime_error("Invalid file format - ENDIANNESS_SANITY_INTEGER mismatch");
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
