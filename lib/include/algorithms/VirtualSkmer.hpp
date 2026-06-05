#include <vector>
#include <array>
#include <fstream>
#include <stdexcept>
#include <forward_list>
#include <algorithm>
#include <cstdint>
#include <atomic>
#include <mutex>
#include <memory>


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
// Doubles as a format-version marker: bumped from "VSKMER_M" to "VSKMER_2" when the
// on-disk minimizer slot became φ-permuted (minimizer-order hash). Old raw-order lists
// carry the previous magic, so load() rejects them (and vice-versa) instead of silently
// returning wrong query results.
constexpr uint64_t ENDIANNESS_SANITY_INTEGER = 0x56534B4D45525F32ULL; // "VSKMER_2" in ASCII
// Bumped to "VSKMER_3" when the on-disk layout gained a per-bucket directory: the sorted
// list is split into minimizer-prefix buckets, each a separately queryable sub-list. A V2
// file stays readable as a single bucket spanning the whole list (see load() / the reader).
constexpr uint64_t ENDIANNESS_SANITY_INTEGER_V3 = 0x56534B4D45525F33ULL; // "VSKMER_3" in ASCII
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

inline void print_query_results(const std::vector<std::vector<uint8_t>> & result_vector, std::ostream & os = std::cout){
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

// One per-bucket directory entry of the VSKMER_3 on-disk format. `mini_lower_bound` is the
// inclusive lower bound (in φ-permuted minimizer space) of the bucket's minimizer interval;
// entries are stored for every bucket (empty ones too) in strictly increasing order so a
// query routes with a single upper_bound. `count` is the bucket's super-k-mer count, from
// which the reader derives byte offsets by prefix sum.
struct BucketDirEntry {
    uint64_t mini_lower_bound;
    uint64_t count;
};
static_assert(sizeof(BucketDirEntry) == 2 * sizeof(uint64_t), "BucketDirEntry must be tightly packed");

// Walk outward from `position_in_list` (down to `minimum`, then up to `maximum`) to the
// closest super-k-mer that actually carries a k-mer at `kmer_position_in_skmer`. Span-based
// twin of the former SortedVirtualSkmerList member, so the search below works on any
// contiguous skmer range (the whole in-RAM list or a single lazily-loaded bucket).
template<typename kuint>
inline int64_t find_closest_valid_skmer_in_span(
    const SkmerManipulator<kuint>& manip, const Skmer<kuint>* list,
    const uint64_t position_in_list, const uint64_t minimum, const uint64_t maximum,
    const uint64_t kmer_position_in_skmer)
{
    // `i` is unsigned, so the down-scan guards on `i == minimum` (a `>= minimum` test would
    // wrap past 0 to 2^64-1 and read out of bounds).
    for (uint64_t i {position_in_list}; ; i--){
        if (manip.has_valid_kmer(list[i], kmer_position_in_skmer)) return (int64_t)i;
        if (i == minimum) break;
    }
    for (uint64_t i {position_in_list}; i <= maximum; i++){
        if (manip.has_valid_kmer(list[i], kmer_position_in_skmer)) return (int64_t)i;
    }
    return -1L;
}

// Per-column dichotomic search: for every valid k-mer position of `query`, binary-search the
// sorted skmer range [list, list+list_size) and report whether that k-mer is present. One
// parallel binary search per k-mer position shares probes for cache locality. Shared by the
// in-RAM SortedVirtualSkmerList::query_skmer and the disk-backed BucketedSkmerListReader, so
// the algorithm lives in exactly one place. Returns one flag per queried k-mer position.
// Fills `out` with one presence flag per valid k-mer position of `query`, reusing out's capacity so
// callers can avoid a per-query heap allocation. search_kmers_in_span() is the wrapper that returns
// a fresh vector (in-RAM list path and single-query callers).
template<typename kuint>
inline void search_kmers_in_span_into(
    const SkmerManipulator<kuint>& manip, const Skmer<kuint>* list, const int64_t list_size,
    const Skmer<kuint>& query, std::vector<uint8_t>& out)
{
    auto [query_start_position, query_end_position] = manip.get_valid_kmer_bounds(query);

    if(query_end_position < query_start_position) {
        out.clear();
        return;
    }
    const uint64_t tot_num_kmers_to_search {query_end_position - query_start_position + 1};

    // RETURN ON EDGE CASES (SKMER DOES NOT CONTAIN A KMER OR LIST/BUCKET IS EMPTY)
    if (tot_num_kmers_to_search <= 0 || list_size == 0){
        out.assign(std::max(0UL, tot_num_kmers_to_search), 0);
        return;
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
    // ARRAYS INITIALIZATION
    std::fill_n(keep_searching, tot_num_kmers_to_search, 1);
    std::fill_n(binary_search_boundaries, tot_num_kmers_to_search,
        std::make_pair(0LL, static_cast<int64_t>(list_size - 1)));
    for(size_t i {0}; i < tot_num_kmers_to_search; i++){
        positions_to_search[i] = i;
    }
    // START BINARY SEARCH
    while(num_kmers_to_search > 0){
        // VERIFY THAT THE CURRENT PRIORITY_POSITION IS STILL VALID AND UPDATE
        if (keep_searching[current_priority_offset] == 0) {
            current_priority_offset = km::sortedlist::util::get_new_priority_value(current_priority_offset, positions_to_search, searchable_position_count, keep_searching, tot_num_kmers_to_search);
        }

        // UPDATE MEAN
        const int64_t old_mean {mean};
        mean = (binary_search_boundaries[current_priority_offset].first + binary_search_boundaries[current_priority_offset].second) >> 1;
        if (mean == old_mean){
            int64_t new_pos = find_closest_valid_skmer_in_span(manip, list, mean, binary_search_boundaries[current_priority_offset].first, binary_search_boundaries[current_priority_offset].second, query_start_position + current_priority_offset);
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

        auto [queried_start_position, queried_end_position] = manip.get_valid_kmer_bounds(list[mean]);

        for(uint64_t i {0}; i < searchable_position_count; i++){
            const uint64_t offset {positions_to_search[i]};
            current_searched_position_in_skmer = query_start_position + offset;
            if(current_searched_position_in_skmer >= queried_start_position && current_searched_position_in_skmer <= queried_end_position){

                const int kmer_comparison {manip.kmer_compare(query, list[mean], current_searched_position_in_skmer)};
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

    out.assign(result, result + tot_num_kmers_to_search);
}

// Wrapper returning a fresh vector (in-RAM list path and single-query callers).
template<typename kuint>
inline std::vector<uint8_t> search_kmers_in_span(
    const SkmerManipulator<kuint>& manip, const Skmer<kuint>* list, const int64_t list_size,
    const Skmer<kuint>& query)
{
    std::vector<uint8_t> out;
    search_kmers_in_span_into<kuint>(manip, list, list_size, query, out);
    return out;
}

// VIRTUAL SUPERKMER CLASS
template<typename kuint>
struct Virtual_skmer {
    Skmer<kuint> skmer;
    // Index of the last enumerated super-k-mer merged into this virtual skmer.
    // It indexes the per-bucket (or whole-input, in-RAM path) enumeration, whose
    // size is far below 2^32, so 32 bits suffice and drop sizeof from 40 to 32
    // bytes (the 8-byte alignment floor) — shrinking the two ping-pong merge
    // buffers by 20%. Kept inline (not a parallel side-array) so the merge keeps
    // its locality.
    uint32_t last_id;
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

    // Thin wrapper over the shared span search on the whole in-RAM list.
    std::vector<uint8_t> query_skmer(const Skmer<kuint> query) const{
        return km::sortedlist::search_kmers_in_span<kuint>(
            m_manip, m_skmer_list.data(), static_cast<int64_t>(m_skmer_list.size()), query);
    }

    // Query a group of super-k-mers, one after another, returning per-skmer results.
    std::vector<std::vector<uint8_t>> query_skmer_batch(std::vector<Skmer<kuint>> query_skmers) const{
        std::vector<std::vector<uint8_t>> result(query_skmers.size());
        for (size_t i {0}; i < query_skmers.size(); ++i)
            result[i] = query_skmer(query_skmers[i]);
        return result;
    }

    void query(const std::string filename, std::ostream& os = std::cout) {
        constexpr uint64_t MAX_INGESTED_SKMER {4096};
        // Enumerate super-k-mers from the file and query them one after another,
        // flushing results in bounded groups to keep memory flat.
        km::FileSkmerator<kuint> file_skmerator {m_manip, filename};

        std::vector<km::Skmer<kuint>> buffer;
        buffer.reserve(MAX_INGESTED_SKMER);

        for (km::Skmer<kuint> const skmer : file_skmerator){
          buffer.emplace_back(skmer);
          if (buffer.size() == MAX_INGESTED_SKMER){
            km::sortedlist::util::print_query_results(query_skmer_batch(buffer), os);
            buffer.clear();
          }
        }

        // taking care of the last elements, if present
        if (!buffer.empty())
          km::sortedlist::util::print_query_results(query_skmer_batch(buffer), os);
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
    LList m_merge_scratch;                            // Reused merge target (ping-pong with the working list)

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
        using keyed_pos = std::pair<kpair, uint64_t>;

        // Index the right column's prefix (k-1)-mers in a sorted array of
        // (key, right_pos), replacing the per-column hash map (whose nodes
        // dominated peak RAM on repeat-rich buckets). Sorting by (key, pos) keeps
        // each key group's positions ascending.
        std::vector<keyed_pos> right_keys;
        right_keys.reserve(right_column.size());
        for (uint64_t pos {0}; pos < right_column.size(); pos++) {
            assert(right_column[pos] < skmer_enumeration.size());
            right_keys.emplace_back(
                m_manip.extract_prefix_suffix(skmer_enumeration[right_column[pos]], left_position + 1), pos);
        }
        std::sort(right_keys.begin(), right_keys.end(),
            [](keyed_pos const& a, keyed_pos const& b) {
                return a.first == b.first ? a.second < b.second : a.first < b.first;
            });

        // Walk the left column in order; each left suffix (k-1)-mer matches a
        // contiguous range of right_keys (binary search). Emitting (suffix_pos,
        // right_pos) this way reproduces the previous (left asc, right asc) order.
        std::vector<std::pair<uint64_t,uint64_t> > candidate_overlaps;
        for (uint64_t suffix_pos {0}; suffix_pos < left_column.size(); suffix_pos++) {
            assert(left_column[suffix_pos] < skmer_enumeration.size());
            const kpair suffix {m_manip.extract_prefix_suffix(skmer_enumeration[left_column[suffix_pos]], left_position + 1)};
            auto it = std::lower_bound(right_keys.begin(), right_keys.end(), suffix,
                [](keyed_pos const& e, kpair const& k) { return e.first < k; });
            for (; it != right_keys.end() && it->first == suffix; ++it)
                candidate_overlaps.emplace_back(suffix_pos, it->second);
        }
        return candidate_overlaps;
    }

    // Thin wrapper over the shared span helper, kept for the in-RAM list and its tests.
    int64_t find_closest_valid_skmer(const uint64_t position_in_list, const uint64_t minimum, const uint64_t maximum, const uint64_t kmer_position_in_skmer) const{
        return km::sortedlist::find_closest_valid_skmer_in_span<kuint>(
            m_manip, m_skmer_list.data(), position_in_list, minimum, maximum, kmer_position_in_skmer);
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

    // Merge into a reused scratch buffer (ping-pong): we read from `list` and
    // write to `m_merge_scratch`, then swap. Reusing the member buffer across
    // all k-m columns avoids reallocating a list-sized buffer each column,
    // which otherwise churns RSS and fragments the heap.
    LList& merged = m_merge_scratch;
    merged.clear();
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

    // Hand the merged result back to the caller and keep the now-stale `list`
    // storage as the scratch buffer for the next column.
    std::swap(list, merged);
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

        // VSKMER_3 header with a single bucket covering the whole list (lower bound 0).
        // Written inline with no seek, so non-seekable targets like /dev/stdout still work.
        const uint64_t count = list.size();
        const uint64_t n_buckets = 1;
        outFile.write(reinterpret_cast<const char*>(&km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3), sizeof(uint64_t));
        outFile.write(reinterpret_cast<const char*>(&list.m_manip.k), sizeof(uint64_t));
        outFile.write(reinterpret_cast<const char*>(&list.m_manip.m), sizeof(uint64_t));
        outFile.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        outFile.write(reinterpret_cast<const char*>(&n_buckets), sizeof(uint64_t));
        const BucketDirEntry only{0, count};
        outFile.write(reinterpret_cast<const char*>(&only), sizeof(BucketDirEntry));

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

    // ---- Incremental binary serialization (VSKMER_3, bucketed) ----
    // Streams the payload out one minimizer bucket at a time without ever holding the whole
    // list in RAM. Header layout: ENDIANNESS(8) + k(8) + m(8) + count(8) + n_buckets(8), then
    // the per-bucket directory (n_buckets × BucketDirEntry), then the raw Skmer payload in
    // bucket-id order. The total count and the directory are unknown when the header is
    // written, so reserve them and patch at the end. Requires a seekable (regular) file.
    static constexpr std::streamoff COUNT_OFFSET = 3 * sizeof(uint64_t);     // after endianness, k, m
    static constexpr std::streamoff DIRECTORY_OFFSET = 5 * sizeof(uint64_t); // after endianness, k, m, count, n_buckets

    // Byte size of the header (everything before the skmer payload) for a bucket count.
    static std::streamoff header_bytes(uint64_t n_buckets) {
        return DIRECTORY_OFFSET +
               static_cast<std::streamoff>(n_buckets) * static_cast<std::streamoff>(sizeof(BucketDirEntry));
    }

    static void write_header(std::ofstream& out, uint64_t k, uint64_t m, uint64_t count, uint64_t n_buckets) {
        out.write(reinterpret_cast<const char*>(&km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&k), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&m), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&n_buckets), sizeof(uint64_t));
        // Reserve the directory; patch_directory() fills it once per-bucket counts are known.
        const std::vector<BucketDirEntry> placeholder(n_buckets, BucketDirEntry{0, 0});
        if (n_buckets)
            out.write(reinterpret_cast<const char*>(placeholder.data()),
                      static_cast<std::streamsize>(n_buckets * sizeof(BucketDirEntry)));
    }

    static void append_payload(std::ofstream& out, const std::vector<Skmer<kuint>>& list) {
        if (list.empty()) return;
        out.write(reinterpret_cast<const char*>(list.data()),
                  static_cast<std::streamsize>(list.size() * sizeof(Skmer<kuint>)));
    }

    static void patch_count(std::ofstream& out, uint64_t total) {
        out.seekp(COUNT_OFFSET, std::ios::beg);
        out.write(reinterpret_cast<const char*>(&total), sizeof(uint64_t));
    }

    // Overwrite the reserved directory with the final per-bucket entries (one per bucket id,
    // empty buckets included, in increasing minimizer-lower-bound order).
    static void patch_directory(std::ofstream& out, const std::vector<BucketDirEntry>& dir) {
        out.seekp(DIRECTORY_OFFSET, std::ios::beg);
        if (!dir.empty())
            out.write(reinterpret_cast<const char*>(dir.data()),
                      static_cast<std::streamsize>(dir.size() * sizeof(BucketDirEntry)));
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
            // The stored minimizer slot is φ-permuted; undo it (φ⁻¹) before decoding so
            // the printed nucleotides are the real minimizer, not the mixed value. Bind a
            // const ref for the decode below (a non-const m_pair would make pair's implicit
            // operator uint64_t() viable and the `>>` overload ambiguous).
            Skmer<kuint> sk_decoded = list.m_skmer_list[i];
            list.m_manip.unpermute_minimizer_slot(sk_decoded);
            const Skmer<kuint>& sk = sk_decoded;

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

        // Accept both VSKMER_2 (single global list) and VSKMER_3 (bucketed). load() returns
        // the whole concatenated list either way (it ignores the bucket directory); the lazy
        // per-bucket reader lives in BucketedSkmerListReader.
        const bool is_v3 = (read_endianess_int == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3);
        const bool is_v2 = (read_endianess_int == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER);
        if (!is_v2 && !is_v3) {
            uint64_t swapped_endianess_int = km::sortedlist::util::swap_endian(read_endianess_int);
            if (swapped_endianess_int == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER ||
                swapped_endianess_int == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3) {
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

        if (is_v3) {
            // Skip n_buckets + the per-bucket directory; we read the full payload below.
            uint64_t n_buckets;
            inFile.read(reinterpret_cast<char*>(&n_buckets), sizeof(uint64_t));
            if (inFile.fail()) {
                throw std::runtime_error("Error reading bucket count from file: " + filename);
            }
            inFile.seekg(static_cast<std::streamoff>(n_buckets) * static_cast<std::streamoff>(sizeof(BucketDirEntry)),
                         std::ios::cur);
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


// Disk-backed, lazily-loaded reader over a bucketed (VSKMER_3) sorted skmer list. Each query
// is routed to its minimizer-prefix bucket and only that bucket's sub-list is read from disk
// (and cached), so query RAM is bounded by the touched buckets, not the whole index. A
// VSKMER_2 file is read as a single bucket spanning the whole list (queries load it all on
// first touch). Move-only (owns an ifstream); obtained via open().
template<typename kuint>
class BucketedSkmerListReader {
public:
    static BucketedSkmerListReader open(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if (in.fail())
            throw std::runtime_error("Error opening file for reading: " + filename);

        uint64_t magic;
        in.read(reinterpret_cast<char*>(&magic), sizeof(uint64_t));
        if (in.fail())
            throw std::runtime_error("Error reading magic number from file: " + filename);

        const bool is_v3 = (magic == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3);
        const bool is_v2 = (magic == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER);
        if (!is_v2 && !is_v3) {
            const uint64_t swapped = km::sortedlist::util::swap_endian(magic);
            if (swapped == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER ||
                swapped == km::sortedlist::util::ENDIANNESS_SANITY_INTEGER_V3)
                throw std::runtime_error("Endianness mismatch - file was written on a system with different endianness");
            throw std::runtime_error("Invalid file format - ENDIANNESS_SANITY_INTEGER mismatch");
        }

        uint64_t file_k, file_m, count;
        in.read(reinterpret_cast<char*>(&file_k), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&file_m), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&count), sizeof(uint64_t));
        if (in.fail())
            throw std::runtime_error("Error reading header from file: " + filename);

        std::vector<BucketDirEntry> dir;
        std::streamoff payload_start;
        if (is_v3) {
            uint64_t n_buckets;
            in.read(reinterpret_cast<char*>(&n_buckets), sizeof(uint64_t));
            if (in.fail())
                throw std::runtime_error("Error reading bucket count from file: " + filename);
            dir.resize(n_buckets);
            if (n_buckets)
                in.read(reinterpret_cast<char*>(dir.data()),
                        static_cast<std::streamsize>(n_buckets * sizeof(BucketDirEntry)));
            if (in.fail())
                throw std::runtime_error("Error reading bucket directory from file: " + filename);
            payload_start = VirtualSkmerSerializer<kuint>::header_bytes(n_buckets);
        } else {
            dir.push_back(BucketDirEntry{0, count});
            payload_start = static_cast<std::streamoff>(4 * sizeof(uint64_t));
        }

        return BucketedSkmerListReader(file_k, file_m, std::move(in), std::move(dir), payload_start);
    }

    uint64_t k() const { return m_manip.k; }
    uint64_t m() const { return m_manip.m; }
    uint64_t n_buckets() const { return m_n_buckets; }

    // Bucket id whose minimizer interval contains the query's φ-permuted minimizer. m_lower is
    // strictly increasing with m_lower[0] == 0, so upper_bound() - 1 lands on the right bucket.
    // Const and touches only immutable routing state, so it is safe to call concurrently (the
    // parallel query producer uses it to bucketize queries before handing them to consumers).
    uint64_t bucket_of(const Skmer<kuint>& query) const {
        if (m_n_buckets <= 1) return 0;
        const uint64_t mini = static_cast<uint64_t>(m_manip.minimizer(query));
        const auto it = std::upper_bound(m_lower.begin(), m_lower.end(), mini);
        return static_cast<uint64_t>(it - m_lower.begin()) - 1;
    }

    // Look up every k-mer of `query`, restricted to its bucket's sub-list (loaded on demand).
    std::vector<uint8_t> query_skmer(const Skmer<kuint>& query) {
        const std::vector<Skmer<kuint>>& span = bucket(bucket_of(query));
        return km::sortedlist::search_kmers_in_span<kuint>(
            m_manip, span.data(), static_cast<int64_t>(span.size()), query);
    }

    // Like query_skmer but writes flags into `out`, reusing its capacity (no per-query allocation).
    void query_skmer_into(const Skmer<kuint>& query, std::vector<uint8_t>& out) {
        const std::vector<Skmer<kuint>>& span = bucket(bucket_of(query));
        km::sortedlist::search_kmers_in_span_into<kuint>(
            m_manip, span.data(), static_cast<int64_t>(span.size()), query, out);
    }

    std::vector<std::vector<uint8_t>> query_skmer_batch(const std::vector<Skmer<kuint>>& query_skmers) {
        std::vector<std::vector<uint8_t>> result(query_skmers.size());
        for (size_t i {0}; i < query_skmers.size(); ++i)
            result[i] = query_skmer(query_skmers[i]);
        return result;
    }

    // Enumerate super-k-mers from a file and query them, flushing results in bounded groups.
    void query(const std::string& filename, std::ostream& os = std::cout) {
        constexpr uint64_t MAX_INGESTED_SKMER {4096};
        km::FileSkmerator<kuint> file_skmerator {m_manip, filename};

        std::vector<km::Skmer<kuint>> buffer;
        buffer.reserve(MAX_INGESTED_SKMER);

        for (km::Skmer<kuint> const skmer : file_skmerator){
            buffer.emplace_back(skmer);
            if (buffer.size() == MAX_INGESTED_SKMER){
                km::sortedlist::util::print_query_results(query_skmer_batch(buffer), os);
                buffer.clear();
            }
        }

        if (!buffer.empty())
            km::sortedlist::util::print_query_results(query_skmer_batch(buffer), os);
    }

private:
    BucketedSkmerListReader(uint64_t k, uint64_t m, std::ifstream&& in,
                            std::vector<BucketDirEntry> dir, std::streamoff payload_start)
        : m_manip(k, m), m_in(std::move(in)), m_n_buckets(dir.size())
    {
        m_lower.resize(m_n_buckets);
        m_count.resize(m_n_buckets);
        m_byte_offset.resize(m_n_buckets);
        std::streamoff off = payload_start;
        for (uint64_t b {0}; b < m_n_buckets; ++b) {
            m_lower[b] = dir[b].mini_lower_bound;
            m_count[b] = dir[b].count;
            m_byte_offset[b] = off;
            off += static_cast<std::streamoff>(dir[b].count) * static_cast<std::streamoff>(sizeof(Skmer<kuint>));
        }
        m_cache.assign(m_n_buckets, {});
        m_loaded = std::make_unique<std::atomic<uint8_t>[]>(m_n_buckets); // value-initialized to 0
        m_load_mtx = std::make_unique<std::mutex>();
    }

    // Lazily read bucket `b`'s sub-list from disk into the cache, then return it. Thread-safe:
    // cache hits are lock-free (acquire-load of the per-bucket flag); the first load of a bucket
    // takes m_load_mtx so concurrent consumers don't race on m_in / m_cache. m_cache is pre-sized
    // and never resized, and each entry is written exactly once before its flag is released, so the
    // returned reference stays valid for the reader's lifetime even as other buckets load.
    const std::vector<Skmer<kuint>>& bucket(uint64_t b) {
        if (m_loaded[b].load(std::memory_order_acquire))
            return m_cache[b];
        std::lock_guard<std::mutex> lock(*m_load_mtx);
        if (!m_loaded[b].load(std::memory_order_relaxed)) {
            std::vector<Skmer<kuint>>& dst = m_cache[b];
            dst.resize(m_count[b]);
            if (m_count[b]) {
                m_in.clear(); // drop any EOF/fail flag from a previous read
                m_in.seekg(m_byte_offset[b], std::ios::beg);
                m_in.read(reinterpret_cast<char*>(dst.data()),
                          static_cast<std::streamsize>(m_count[b] * sizeof(Skmer<kuint>)));
                if (m_in.fail())
                    throw std::runtime_error("Error reading bucket payload from sorted skmer list");
            }
            m_loaded[b].store(1, std::memory_order_release);
        }
        return m_cache[b];
    }

    SkmerManipulator<kuint> m_manip;
    std::ifstream m_in;
    uint64_t m_n_buckets {1};
    std::vector<uint64_t> m_lower;        // per-bucket minimizer lower bound (routing table)
    std::vector<uint64_t> m_count;        // per-bucket super-k-mer count
    std::vector<std::streamoff> m_byte_offset; // per-bucket payload offset in the file
    std::vector<std::vector<Skmer<kuint>>> m_cache; // lazily-filled per-bucket sub-lists
    // Per-bucket "is m_cache[b] populated" flags, atomic so query threads can check without taking
    // the lock. Held via unique_ptr (atomics aren't movable) to keep the reader move-constructible.
    std::unique_ptr<std::atomic<uint8_t>[]> m_loaded;
    std::unique_ptr<std::mutex> m_load_mtx; // serializes the first disk load of each bucket
};


} // namespace sortedlist
} // namespace km

#endif
