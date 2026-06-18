#include <iostream>
#include <sstream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <array>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include "test_data.hpp"

using namespace std;
using kuint = uint16_t;
using kpair = km::Skmer<kuint>::pair;

std::array<std::array<km::Skmer<kuint>, 3>, 6> get_skmer_permutations(std::array<kpair, 3> const &kmer_triplet)
{

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    km::SkmerPrettyPrinter<kuint> pp{k, m};

    const std::array<std::string, 3> skmer_strings{"AAAAA", "AACCC", "CCCCC"};

    std::array<std::array<km::Skmer<kuint>, 3>, 6> permutation_array{};

    permutation_array[0] = {km::Skmer<kuint>(kmer_triplet[0], 0, 3),
                            km::Skmer<kuint>(kmer_triplet[1], 0, 3),
                            km::Skmer<kuint>(kmer_triplet[2], 0, 3)};

    for (int64_t i{1}; i < 6; i++)
    {
        permutation_array[i] = permutation_array[i - 1];
        std::next_permutation(permutation_array[i].begin(), permutation_array[i].end());
    }

    return permutation_array;
}


/* TESTING VIRTUAL SUPER-K-MER STRUCTURE
*/
TEST(VirtualSkmer, GenerateVirtualSkmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip {k, m};
    //                    Prefix:      A   C   G   T
    //                    Suffix:    A   C   T   C
    const kpair input_skmer {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer(input_skmer,3,4);
    std::vector< km::sortedlist::Virtual_skmer<kuint> > extracted_skmers;

    for(int position = 0; position < 4; position++){
        extracted_skmers.emplace_back(m_skmer,manip, position, position);
    }

    //                    Prefix:      A   C   G   T             A   C   G   _            A   C   _   _
    //                    Suffix:    A   _   _   _             A   C   _   _             A   C   T   _
    const kpair expected_kpairs[4] {{0b0000110111111110U, 0}, {0b0000010111111111U, 0}, {0b0000010110111111U, 0},
    //                                 A   _   _   _
    //                               A   C   T   C
                                    {0b0000011110110111U, 0}};

    std::vector< km::sortedlist::Virtual_skmer<kuint> > expected_virtual_skmers {
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[0],3,0,0),
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[1],2,1,1),
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[2],1,2,2),
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[3],0,3,3)
    };

    for(size_t i {0}; i < expected_virtual_skmers.size(); i++){
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pair,extracted_skmers[i].skmer.m_pair);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pref_size,extracted_skmers[i].skmer.m_pref_size);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_suff_size,extracted_skmers[i].skmer.m_suff_size);
        ASSERT_EQ(expected_virtual_skmers[i].last_id, extracted_skmers[i].last_id);
    }

}

TEST(VirtualSkmer, GenerateVirtualSkmer2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{3};

    km::SkmerManipulator<kuint> manip {k, m};

    //              Prefix:      A   C   G   T
    //               Suffix:       C   T   C
    const kpair input_skmer {0b0000010110110110U, 0};
    km::Skmer<kuint> m_skmer(input_skmer,2,2);

    std::vector< km::sortedlist::Virtual_skmer<kuint> > extracted_skmers;

    for(int position = 0; position < 3; position++){
        // extracted_skmers.emplace_back(m_skmer,manip,0,position);
        extracted_skmers.emplace_back(m_skmer,manip, position, position);
    }

    //                    Prefix:      A   C   G   T             A   C   G   _            A   C   _   _
    //                    Suffix:        C   _   _                 C   T   _                C   T   C
    const kpair expected_kpairs[3] {{0b0000010111111110U, 0}, {0b0000010110111111U, 0}, {0b0000010110110111U, 0}};

    std::vector< km::sortedlist::Virtual_skmer<kuint> > expected_virtual_skmers {
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[0],2,0,0),
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[1],1,1,1),
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[2],0,2,2)};

    for(size_t i {0}; i < expected_virtual_skmers.size(); i++){
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pair,extracted_skmers[i].skmer.m_pair);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pref_size,extracted_skmers[i].skmer.m_pref_size);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_suff_size,extracted_skmers[i].skmer.m_suff_size);
    }
}

/** Testing the add_kmer method of the virtual super-k-mer
 */
TEST(VirtualSkmer, AddKmerToVirtualSkmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    km::SkmerPrettyPrinter<kuint> pp{k, m};
    //                    Prefix:      A   C   G   T
    //                    Suffix:    A   C   T   C
    const kpair m_vector_pair{0b0000010110110110U, 0};

    //                    Prefix:      A   C   G   _
    //                    Suffix:    A   C   _   _
    const kpair m_basic_pair{0b0000010111111111U, 0};
    std::vector<km::Skmer<kuint>> m_skmer_vector{
        km::Skmer<kuint>(m_vector_pair, 3, 3),
        km::Skmer<kuint>(m_basic_pair, 2, 1),
    };

    std::vector<km::sortedlist::Virtual_skmer<kuint>> extracted_skmers;
    uint64_t skmer_id_take_skmer{0};
    uint64_t skmer_id_base_skmer{1};
    uint64_t kmers_to_add[2]{2, 3};

    km::sortedlist::Virtual_skmer<kuint> new_virtual_skmer(m_skmer_vector[1], skmer_id_base_skmer);
    for (uint64_t kmer_pos : kmers_to_add)
    {
        new_virtual_skmer.add_kmer(m_skmer_vector, manip, skmer_id_take_skmer, kmer_pos);
        extracted_skmers.push_back(new_virtual_skmer);
    }

    //                      Prefix:      A   C   G   _             A   C   G   _
    //                      Suffix:    A   C   T   _             A   C   T   C
    const kpair expected_kpairs[4]{{0b0000010110111111U, 0}, {0b0000010110110111U, 0}};

    // km::sorting::Virtual_skmer<kuint>(expected_kpairs[0],4,2,0),
    std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_virtual_skmers{
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[0], 2, 2, 0),
        km::sortedlist::Virtual_skmer<kuint>(expected_kpairs[1], 2, 3, 0)};

    for (size_t i{0}; i < expected_virtual_skmers.size(); i++)
    {
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pair, extracted_skmers[i].skmer.m_pair);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_pref_size, extracted_skmers[i].skmer.m_pref_size);
        ASSERT_EQ(expected_virtual_skmers[i].skmer.m_suff_size, extracted_skmers[i].skmer.m_suff_size);
        // ASSERT_EQ(expected_virtual_skmers[i].last_id, extracted_skmers[i].last_id);
    }
}

/** Testing the result of "has_valid_kmer" on 2 skmers.
 * TODO: move to the manipulator test file.
 */
TEST(SortedVirtualSkmerListTest, KmerValidation)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::SkmerManipulator<kuint> manip{k, m};
    km::SkmerPrettyPrinter<kuint> pp{k, m};

    //                  Prefix:         A   T   _   _             A   _   _   _
    //                  Suffix:       A   C   C   C             C   C   C   C
    const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
    std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3)};

    const uint64_t kmer_positions{k - m + 1};

    const uint64_t expected_valid_kmers[kmer_positions][2]{{0, 0}, {0, 0}, {1, 0}, {1, 1}};
    // 0 values map to false, else to true
    bool kmer_validity;
    for (uint64_t skmer_id{0}; skmer_id < 2; skmer_id++)
    {
        for (uint64_t position{0}; position < kmer_positions; position++)
        {
            kmer_validity = manip.has_valid_kmer(skmer_vector[skmer_id], position);
            ASSERT_EQ(kmer_validity, expected_valid_kmers[position][skmer_id]) << (string("has_valid_kmer returned an unexpected value for skmer ") + std::to_string(skmer_id) + " at position " + std::to_string(position));
        }
    }
}

TEST(SortedVirtualSkmerListTest, InputOutput1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
    std::vector<km::Skmer<kuint>> skmer_enumeration{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3)};

    list.add_list(skmer_enumeration);

    const std::string tmp_path = ::testing::TempDir() + "input.bin";
    km::sortedlist::VirtualSkmerSerializer<kuint>::save(list, tmp_path);

    auto loaded_list = km::sortedlist::VirtualSkmerSerializer<kuint>::load(tmp_path);

    std::vector<km::Skmer<kuint>> m_loaded_list = loaded_list.get_list();

    for (size_t i{0}; i < skmer_enumeration.size(); i++)
    {
        ASSERT_EQ(skmer_enumeration[i].m_pref_size, m_loaded_list[i].m_pref_size);
        ASSERT_EQ(skmer_enumeration[i].m_suff_size, m_loaded_list[i].m_suff_size);
        ASSERT_EQ(skmer_enumeration[i].m_pair, m_loaded_list[i].m_pair);
    }
}

/** Test the order of 2 kmers after a sort on one column */
namespace km
{
    namespace sortedlist
    {
        class SortedVirtualSkmerListPrivateTest : public ::testing::Test
        {
        protected:
        };

        TEST_F(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer1)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   C             C   C   C   C
            const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
            const uint64_t position{0};

            std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3)};

            std::vector<uint64_t> expected_order{};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }
        TEST_F(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer2)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   C             C   C   C   C
            const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
            const uint64_t position{1};

            std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3)};

            std::vector<uint64_t> expected_order{};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer3)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   C             C   C   C   C
            const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
            const uint64_t position{2};

            std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3)};

            std::vector<uint64_t> expected_order{0};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, SortingColumnNoValidKmer4)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   _             C   C   C   C
            const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
            const uint64_t position{3};

            std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2), km::Skmer<kuint>(input_skmers[1], 0, 3)};

            std::vector<uint64_t> expected_order{1};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, SingleKmerSorting)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   C             C   C   C   C
            const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
            const uint64_t position{3};

            std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3)};
            std::vector<uint64_t> expected_order{0, 1};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            // std::vector<uint64_t> ordered_kmers = km::sorting::sort_column(skmer_vector.begin(), skmer_vector.end(), position, manip);
            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }
        TEST_F(SortedVirtualSkmerListPrivateTest, SingleKmerSorting2)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{4};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   T                C   C   C
            //                 Suffix:        A   _   _                A   _   _
            const kpair input_skmers[2]{{0b000011101110U, 0}, {0b000111011101U, 0}};
            const uint64_t position{0};

            std::vector<km::Skmer<kuint>> skmer_vector{
                km::Skmer<kuint>(input_skmers[1], 2, 0),
                km::Skmer<kuint>(input_skmers[0], 2, 0)
            };
            std::vector<uint64_t> expected_order{1, 0};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            // std::vector<uint64_t> ordered_kmers = km::sorting::sort_column(skmer_vector.begin(), skmer_vector.end(), position, manip);
            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, SingleKmerSortingAndUnique)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   C             C   C   C   C
            const kpair input_skmers[3]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0},{0b0000011001110111U, 0},};
            const uint64_t position{3};

            std::vector<km::Skmer<kuint>> skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 0, 3), km::Skmer<kuint>(input_skmers[0], 1, 3)};
            std::vector<uint64_t> expected_order{0, 1};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            // std::vector<uint64_t> ordered_kmers = km::sorting::sort_column(skmer_vector.begin(), skmer_vector.end(), position, manip);
            std::vector<uint64_t> ordered_kmers = m_list.sort_column(skmer_vector.begin(), skmer_vector.end(), position);
            ASSERT_EQ(ordered_kmers, expected_order);
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, SingleKmerSortingReversed)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                 Prefix:          A   T   _   _             A   _   _   _
            //                 Suffix:        A   C   C   C             C   C   C   C
            const kpair input_skmers[2]{{0b0000011001110111U, 0}, {0b0100011101110111U, 0}};
            const uint64_t position{3};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            std::vector<km::Skmer<kuint>> skmer_vector_rev{km::Skmer<kuint>(input_skmers[1], 0, 3), km::Skmer<kuint>(input_skmers[0], 1, 3)};
            std::vector<uint64_t> expected_order_rev{1, 0};
            std::vector<uint64_t> ordered_kmers_rev = m_list.sort_column(skmer_vector_rev.begin(), skmer_vector_rev.end(), position);
            ASSERT_EQ(ordered_kmers_rev, expected_order_rev);
        }

        /** Test sorting on a column for all possible permutations of 3 skmers */
        TEST_F(SortedVirtualSkmerListPrivateTest, ThreeKmerSorting)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                      Prefix:                        A   _   _   _                 A   _   _   _
            //                      Suffix:                      A   A   A   A                 A   C   C   C
            std::array<kpair, 3> const kmer_triplet{kpair(0b0000001100110011U, 0), kpair(0b0000011101110111U, 0),
                                                    //                      Prefix:                        C   _   _   _
                                                    //                      Suffix:                      C   C   C   C
                                                    kpair(0b0101011101110111U, 0)};

            std::vector<uint64_t> ordered_kmers{};
            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            uint64_t position{3};
            for (auto const &permuted : get_skmer_permutations(kmer_triplet))
            {
                ordered_kmers = m_list.sort_column(permuted.begin(), permuted.end(), position);
                uint64_t loop_idx{0};

                for (auto &skmer_position : ordered_kmers)
                {
                    km::Skmer<kuint> const &curr_skmer{permuted[skmer_position]};
                    kpair const &expected_value{kmer_triplet[loop_idx]};
                    ASSERT_EQ(curr_skmer.m_pair, expected_value);

                    loop_idx += 1;
                }
            }
        }

        /** Compare 2 kmers from to adjacent columns and try to overlap them.
         * No overlap is expected in this example.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlapNoOverlap)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};
            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            //                    Prefix:       A   C   _   _             A   _   _   _
            //                    Suffix:     A   C   C   _             A   G   C   C
            const kpair input_skmers[2]{{0b0000010101111111U, 0}, {0b0000111101110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2), km::Skmer<kuint>(input_skmers[1], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0};
            std::vector<uint64_t> right_column_order{1};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{};

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                std::cerr << "Unexpected overlap: " << computed_overlaps[i].first << computed_overlaps[i].second << std::endl;
            }

            ASSERT_EQ(computed_overlaps.size(), 0);
        }

        /** Compare 2 kmers from to adjacent columns and try to overlap them.
         * 1 overlap is expected in this example.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap1Overlap)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:       A   C   _   _             A   _   _   _
            //                    Suffix:     A   C   C   _             A   C   C   C
            const kpair input_skmers[2]{{0b0000010101111111U, 0}, {0b0000011101110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2), km::Skmer<kuint>(input_skmers[1], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0};
            std::vector<uint64_t> right_column_order{1};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 0)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 1 kmer in the left column and 2 in the right column.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap1Left2Right)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   _   _             A   _   _   _             A   _   _   _
            //                    Suffix:    A   C   C   _             A   C   C   C             A   C   C   T
            const kpair input_skmers[3]{{0b0000010101111111U, 0}, {0b0000011101110111U, 0}, {0b0000011101111011U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[1], 0, 3),
                                                            km::Skmer<kuint>(input_skmers[2], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0};
            std::vector<uint64_t> right_column_order{1, 2};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 0), std::pair<uint64_t, uint64_t>(0, 1)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 2 kmers in the left column and 1 in the right column.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left1Right)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   _   _             A   T   _   _             A   _   _   _
            //                    Suffix:    A   C   C   _             A   C   C   _             A   C   C   C
            const kpair input_skmers[3]{{0b0000010101111111U, 0}, {0b0000011001111111U, 0}, {0b0000011101110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[1], 0, 3),
                                                            km::Skmer<kuint>(input_skmers[2], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0, 1};
            std::vector<uint64_t> right_column_order{2};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 0), std::pair<uint64_t, uint64_t>(1, 0)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 2 kmers in the left column and 2 in the right column.
         * Expecting 2 parallel overlaps.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightParallel)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _
            //                    Suffix:    A   C   C   _             A   C   T   _             A   C   C   C
            const kpair input_skmers[4]{{0b0000010101111111U, 0}, {0b0000011010111111U, 0}, {0b0000011101110111U, 0},
                                        //                                 A   _   _   _
                                        //                               A   C   T   C
                                        {0b0000011110110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[1], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[2], 0, 3),
                                                            km::Skmer<kuint>(input_skmers[3], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0, 1};
            std::vector<uint64_t> right_column_order{2, 3};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 0), std::pair<uint64_t, uint64_t>(1, 1)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 2 kmers in the left column and 2 in the right column.
         * Expecting 2 crossing overlaps.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _
            //                    Suffix:    A   C   T   _             A   C   C   _             A   C   C   C
            const kpair input_skmers[4]{{0b0000010110111111U, 0}, {0b0000011001111111U, 0}, {0b0000011101110111U, 0},
                                        //                                 A   _   _   _
                                        //                               A   C   T   C
                                        {0b0000011110110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[1], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[2], 0, 3),
                                                            km::Skmer<kuint>(input_skmers[3], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0, 1};
            std::vector<uint64_t> right_column_order{2, 3};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 1), std::pair<uint64_t, uint64_t>(1, 0)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 2 kmers in the left column and 2 in the right column.
         * Expecting 1 crossing overlaps.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed1)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _
            //                    Suffix:    A   C   A   _             A   C   C   _             A   C   C   C
            const kpair input_skmers[4]{{0b0000010100111111U, 0}, {0b0000011001111111U, 0}, {0b0000011101110111U, 0},
                                        //                                 A   _   _   _
                                        //                               A   C   T   C
                                        {0b0000011110110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[1], 1, 2),
                                                            km::Skmer<kuint>(input_skmers[2], 0, 3),
                                                            km::Skmer<kuint>(input_skmers[3], 0, 3)};
            const uint64_t left_column_position{2};
            std::vector<uint64_t> left_column_order{0, 1};
            std::vector<uint64_t> right_column_order{2, 3};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(1, 0)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 2 kmers in the left column and 2 in the right column.
         * Expecting 2 crossing overlaps.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed2)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   C   C             A   T   T   T             A   T   T   _
            //                    Suffix:    A   _   _   _             A   _   _   _             A   C   _   _
            const kpair input_skmers[4]{{0b0000110111011101U, 0}, {0b0000111011101110U, 0}, {0b0000011011101111U, 0},
                                        //                                 A   C   C   _
                                        //                               A   T   _   _
                                        {0b0000100111011111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 2, 3), km::Skmer<kuint>(input_skmers[1], 2, 3), km::Skmer<kuint>(input_skmers[2], 1, 4), km::Skmer<kuint>(input_skmers[3], 1, 4)};
            const uint64_t left_column_position{0};
            std::vector<uint64_t> left_column_order{0, 1};
            std::vector<uint64_t> right_column_order{2, 3};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 1), std::pair<uint64_t, uint64_t>(1, 0)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        /** Overlap test with 2 columns.
         * 2 kmers in the left column and 2 in the right column.
         * Expecting 2 crossing overlaps.
         */
        TEST_F(SortedVirtualSkmerListPrivateTest, GetCandidateOverlap2Left2RightCrossed3)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            //                    Prefix:      A   C   _   _             A   T   _   _              A   _   _   _
            //                    Suffix:    A   C   T   C             A   C   C   C             A   C   C   C
            const kpair input_skmers[4]{{0b0000010110110111U, 0}, {0b0000011001110111U, 0}, {0b0000011101110111U, 0},
                                        //                                 A   _   _   _
                                        //                               A   C   T   C
                                        {0b0000011110110111U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 1, 3), km::Skmer<kuint>(input_skmers[2], 0, 3), km::Skmer<kuint>(input_skmers[3], 0, 3)};
            const uint64_t left_column_position{3};
            std::vector<uint64_t> left_column_order{0, 1};
            std::vector<uint64_t> right_column_order{2, 3};
            std::vector<std::pair<uint64_t, uint64_t>> computed_overlaps; m_list.get_candidate_overlaps(m_skmer_vector, left_column_position, left_column_order, right_column_order, computed_overlaps);
            std::vector<std::pair<uint64_t, uint64_t>> expected_overlaps{std::pair<uint64_t, uint64_t>(0, 1), std::pair<uint64_t, uint64_t>(1, 0)};

            ASSERT_EQ(computed_overlaps.size(), expected_overlaps.size());

            // candidate-overlap order is an implementation detail (both colinear_chaining and
            // greedy_chaining sort their input), so compare as sets.
            std::sort(computed_overlaps.begin(), computed_overlaps.end());
            std::sort(expected_overlaps.begin(), expected_overlaps.end());
            for (uint64_t i{0}; i < computed_overlaps.size(); i += 1)
            {
                ASSERT_EQ(computed_overlaps[i], expected_overlaps[i]);
            }
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, MergeColumnFromEmpty1)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            // km::SkmerManipulator<kuint> manip {k, m};
            //                    Prefix:      A   C   _   _             A   T   _   _
            //                    Suffix:    A   C   T   C             A   C   C   C
            const kpair input_skmers[2]{{0b0000010110110111U, 0}, {0b0000011001110111U, 0}};
            //  0000010110111111          0000011001111111
            std::vector<km::Skmer<kuint>> m_skmer_vector{km::Skmer<kuint>(input_skmers[0], 1, 3), km::Skmer<kuint>(input_skmers[1], 1, 3)};
            const std::vector<uint64_t> column{0, 1};
            const std::vector<uint64_t> left_c{};

            std::vector<km::sortedlist::Virtual_skmer<kuint>> computed_list;
            const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps;
            uint64_t column_position{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            m_list.merge_LList_column(m_skmer_vector, computed_list, left_c, column, valid_overlaps, column_position);

            //                       Prefix:      A   C   _   _             A   T   _   _
            //                       Suffix:    A   C   T   _             A   C   C   _
            const kpair expected_skmers[2]{{0b0000010110111111, 0}, {0b0000011001111111, 0}};
            const std::array<std::pair<uint64_t, uint64_t>, 2> prefix_suffix{{std::pair(1, 2), std::pair(1, 2)}};
            const std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_list{
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second, column[0]),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second, column[1])};

            auto expected_it = expected_list.begin();
            auto computed_it = computed_list.begin();

            for (; expected_it != expected_list.end(); expected_it++, computed_it++)
            {
                ASSERT_EQ(expected_it->skmer.m_pair, computed_it->skmer.m_pair);
                ASSERT_EQ(expected_it->skmer.m_pref_size, computed_it->skmer.m_pref_size);
                ASSERT_EQ(expected_it->skmer.m_suff_size, computed_it->skmer.m_suff_size);
                ASSERT_EQ(expected_it->last_id, computed_it->last_id);
            }
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, MergeColumnFromEmpty2)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            //                    Prefix:      A   G   C   _             A   G   C   T
            //                    Suffix:    A   C   _   _             A   _   _   _
            const kpair input_skmers[2]{{0b0000011111011111U, 0}, {0b0000111111011110U, 0}};

            std::vector<km::Skmer<kuint>> m_skmer_vector{
                km::Skmer<kuint>(input_skmers[0], 2, 1),
                km::Skmer<kuint>(input_skmers[1], 3, 0)
            };
            const std::vector<uint64_t> column{0};
            const std::vector<uint64_t> left_c{1};

            std::vector<km::sortedlist::Virtual_skmer<kuint>> computed_list {
                km::sortedlist::Virtual_skmer<kuint>(m_skmer_vector[1], 1)};
            const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps {std::pair(0,0)};
            uint64_t column_position{1};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);

            m_list.merge_LList_column(m_skmer_vector, computed_list, left_c, column, valid_overlaps, column_position);

            //                       Prefix:      A   G   C   T
            //                       Suffix:    A   C   _   _
            const kpair expected_skmers[1]{{0b0000011111011110U, 0}};
            const std::pair<uint64_t, uint64_t> prefix_suffix {std::pair(3, 1)};
            const std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_list{
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[0],
                                                    prefix_suffix.first,
                                                    prefix_suffix.second,
                                                    column[0])};

            ASSERT_EQ(expected_list.size(), computed_list.size());

            auto expected_it = expected_list.begin();
            auto computed_it = computed_list.begin();

            for (; expected_it != expected_list.end(); expected_it++, computed_it++)
            {
                std::cerr << "Computed has prfix: " << computed_it->skmer.m_pref_size << "; suffix: " << computed_it->skmer.m_suff_size << std::endl;
                ASSERT_EQ(expected_it->skmer.m_pair, computed_it->skmer.m_pair);
                ASSERT_EQ(expected_it->skmer.m_pref_size, computed_it->skmer.m_pref_size);
                ASSERT_EQ(expected_it->skmer.m_suff_size, computed_it->skmer.m_suff_size);
                ASSERT_EQ(expected_it->last_id, computed_it->last_id);
            }
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, MergeLeftRightColumnElements1)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::SkmerManipulator<kuint> manip(k, m);

            //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G
            //                    Suffix:    A   G   G   G             A   T   G   T             A   C   T   T
            const kpair input_skmers[3]{{0b0000110111011101U, 0}, {0b0000100111011010U, 0}, {0b0000010110011011U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{
                km::Skmer<kuint>(input_skmers[0], 3, 3),
                km::Skmer<kuint>(input_skmers[1], 3, 3),
                km::Skmer<kuint>(input_skmers[2], 3, 3)};

            const std::vector<uint64_t> column{2, 1, 0};
            const std::vector<uint64_t> left_c{0, 1};

            std::vector<km::sortedlist::Virtual_skmer<kuint>> computed_list{
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0], 3, 0), 0), 0),
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1], 3, 0), 0), 1),
            };

            const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps{std::pair(0, 0), std::pair(1, 1)};
            uint64_t column_position{1};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            m_list.merge_LList_column(m_skmer_vector, computed_list, left_c, column, valid_overlaps, column_position);

            //                       Prefix:     A   C   C   C           A   C   C   T            A   C   C   _
            //                       Suffix:   A   C   _   _           A   T   _   _            A   G   _   _
            const kpair expected_skmers[3]{{0b0000010111011101, 0}, {0b0000100111011110, 0}, {0b0000110111011111, 0}};
            const std::array<std::pair<uint64_t, uint64_t>, 3> prefix_suffix{{std::pair(3, 1), std::pair(3, 1), std::pair(2, 1)}};
            const std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_list{
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second, column[0]),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second, column[1]),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second, column[2])};

             auto expected_it = expected_list.begin();
            auto computed_it = computed_list.begin();

            ASSERT_EQ(expected_list.size(), computed_list.size());

            for (; expected_it != expected_list.end(); expected_it++, computed_it++)
            {
                ASSERT_EQ(expected_it->skmer.m_pair, computed_it->skmer.m_pair);
                ASSERT_EQ(expected_it->skmer.m_pref_size, computed_it->skmer.m_pref_size);
                ASSERT_EQ(expected_it->skmer.m_suff_size, computed_it->skmer.m_suff_size);
                ASSERT_EQ(expected_it->last_id, computed_it->last_id);
            }
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, AddRightColumnElement)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::SkmerManipulator<kuint> manip{k, m};
            //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G
            //                    Suffix:    A   G   G   G             A   G   T   T             A   G   T   T
            const kpair input_skmers[3]{{0b0000110111011101U, 0}, {0b0000110110011010U, 0}, {0b0000110110011011U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{
                km::Skmer<kuint>(input_skmers[0], 3, 3),
                km::Skmer<kuint>(input_skmers[1], 3, 3),
                km::Skmer<kuint>(input_skmers[2], 3, 3)};

            const std::vector<uint64_t> column{1, 2, 0};
            const std::vector<uint64_t> left_c{0, 1, 2};

            std::vector<km::sortedlist::Virtual_skmer<kuint>> computed_list{
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0], 2, 1), 1), 0),
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1], 2, 1), 1), 1)};

            const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps{std::pair(0, 1)};
            uint64_t column_position{2};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            m_list.merge_LList_column(m_skmer_vector, computed_list, left_c, column, valid_overlaps, column_position);

            //                       Prefix:      A   C   _   _            A   C   C   _            A   C   C   _
            //                       Suffix:    A   G   T   _            A   G   T   _            A   G   _   _
            const kpair expected_skmers[3]{{0b0000110110111111, 0}, {0b0000110110011111, 0}, {0b0000110111011111, 0}};
            const std::array<std::pair<uint64_t, uint64_t>, 3> prefix_suffix{{std::pair(1, 2), std::pair(2, 2), std::pair(2, 1)}};
            const std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_list{
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second, column[0]),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second, column[1]),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second, column[0])};

            auto expected_it = expected_list.begin();
            auto computed_it = computed_list.begin();

            for (; expected_it != expected_list.end(); expected_it++, computed_it++)
            {
                ASSERT_EQ(expected_it->skmer.m_pair, computed_it->skmer.m_pair);
                ASSERT_EQ(expected_it->skmer.m_pref_size, computed_it->skmer.m_pref_size);
                ASSERT_EQ(expected_it->skmer.m_suff_size, computed_it->skmer.m_suff_size);
                ASSERT_EQ(expected_it->last_id, computed_it->last_id);
            }
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, AddLeftColumnElement)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::SkmerManipulator<kuint> manip{k, m};
            km::SkmerPrettyPrinter<kuint> pp{k, m};
            //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G
            //                    Suffix:    A   G   G   G             A   G   T   T             A   G   T   T
            const kpair input_skmers[3]{{0b0000110111011101U, 0}, {0b0000110110011010U, 0}, {0b0000110110011011U, 0}};
            std::vector<km::Skmer<kuint>> m_skmer_vector{
                km::Skmer<kuint>(input_skmers[0], 3, 3),
                km::Skmer<kuint>(input_skmers[1], 3, 3),
                km::Skmer<kuint>(input_skmers[2], 3, 3)};

            const std::vector<uint64_t> column{1, 2, 0};
            const std::vector<uint64_t> left_c{1, 0};

            std::vector<km::sortedlist::Virtual_skmer<kuint>> computed_list{
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0], 1, 2), 2), 0),
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1], 1, 2), 2), 1)};

            const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps{std::pair(0, 0)};
            uint64_t column_position{3};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            m_list.merge_LList_column(m_skmer_vector, computed_list, left_c, column, valid_overlaps, column_position);

            //                       Prefix:      A   C   _   _            A   C   _   _            A   _   _   _
            //                       Suffix:    A   G   G   _            A   G   T   T            A   G   T   T
            const kpair expected_skmers[3]{{0b0000110111111111, 0}, {0b0000110110111011, 0}, {0b0000111110111011, 0}};
            const std::array<std::pair<uint64_t, uint64_t>, 3> prefix_suffix{{std::pair(1, 2),
                                                                              std::pair(1, 3),
                                                                              std::pair(0, 3)}};
            const std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_list{
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second, 0),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second, 1),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second, 2)};

            auto expected_it = expected_list.begin();
            auto computed_it = computed_list.begin();

            for (; expected_it != expected_list.end(); expected_it++, computed_it++)
            {
                ASSERT_EQ(expected_it->skmer.m_pair, computed_it->skmer.m_pair);
                ASSERT_EQ(expected_it->skmer.m_pref_size, computed_it->skmer.m_pref_size);
                ASSERT_EQ(expected_it->skmer.m_suff_size, computed_it->skmer.m_suff_size);
                ASSERT_EQ(expected_it->last_id, computed_it->last_id);
            }
        }

        TEST_F(SortedVirtualSkmerListPrivateTest, NoElementPointed)
        {
            using kuint = uint16_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            km::SkmerManipulator<kuint> manip {k, m};
            //                    Prefix:      A   C   C   C             A   C   C   T             A   C   C   G
            //                    Suffix:    A   A   G   G             A   G   T   T             A   A   T   T
            const kpair input_skmers[3] { {0b0000000111011101U, 0}, {0b0000110110011010U, 0}, {0b0000000110011011U, 0}};
            std::vector<km::Skmer<kuint> > m_skmer_vector{
                km::Skmer<kuint>(input_skmers[0],3,3),
                km::Skmer<kuint>(input_skmers[1],3,3),
                km::Skmer<kuint>(input_skmers[2],3,3)
            };

            const std::vector<uint64_t> column{2,1,0};
            const std::vector<uint64_t> left_c{2,1,0};

            std::vector<km::sortedlist::Virtual_skmer<kuint>> computed_list {
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[0],1,2),2),0),
                km::sortedlist::Virtual_skmer<kuint>(manip.get_skmer_of_kmer(km::Skmer<kuint>(input_skmers[1],1,2),2),1)
            };

            const std::vector<std::pair<uint64_t, uint64_t>> valid_overlaps {};
            uint64_t column_position {3};

            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            m_list.merge_LList_column(m_skmer_vector, computed_list, left_c, column, valid_overlaps, column_position);

            //                       Prefix:      A   C   _   _            A   _   _   _            A   C   _   _
            //                       Suffix:    A   A   G   _            A   A   T   T            A   G   T   _
            const kpair expected_skmers[4] { {0b0000000111111111, 0}, {0b0000001110111011, 0}, {0b0000110110111111, 0}, {0b0000111110111011, 0}};
            const std::array<std::pair<uint64_t, uint64_t>, 4> prefix_suffix {{std::pair(1,2), std::pair(0,3), std::pair(1,2), std::pair(0,3)}};
            const std::vector<km::sortedlist::Virtual_skmer<kuint>> expected_list {
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[0], prefix_suffix[0].first, prefix_suffix[0].second , 0),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[1], prefix_suffix[1].first, prefix_suffix[1].second , 2),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[2], prefix_suffix[2].first, prefix_suffix[2].second , 1),
                km::sortedlist::Virtual_skmer<kuint>(expected_skmers[3], prefix_suffix[3].first, prefix_suffix[3].second , 1)};

            auto expected_it = expected_list.begin();
            auto computed_it = computed_list.begin();

            for (; expected_it != expected_list.end(); expected_it++, computed_it++)
            {
                ASSERT_EQ(expected_it->skmer.m_pair, computed_it->skmer.m_pair);
                ASSERT_EQ(expected_it->skmer.m_pref_size, computed_it->skmer.m_pref_size);
                ASSERT_EQ(expected_it->skmer.m_suff_size, computed_it->skmer.m_suff_size);
                ASSERT_EQ(expected_it->last_id, computed_it->last_id);
            }
        }

        // REGRESSION TEST for find_closest_valid_skmer. When no skmer in the
        // [minimum, maximum] window has a valid kmer at the queried position and
        // minimum == 0, the descending loop underflowed its unsigned counter
        // past 0 (the `i >= minimum` guard is always true once minimum is 0),
        // reading m_skmer_list[2^64-1] — 8 bytes before the buffer — and
        // aborting under AddressSanitizer. It must stop at index 0 and return -1.
        TEST_F(SortedVirtualSkmerListPrivateTest, FindClosestValidSkmerUnderflow)
        {
            using kuint = uint64_t;
            using kpair = km::Skmer<kuint>::pair;

            constexpr uint64_t k{5};
            constexpr uint64_t m{2};

            // Three skmers with suffix size 0: none has a valid kmer at
            // position k-m = 3 (has_valid_kmer requires suff_size >= kmer_pos),
            // so the descending scan is forced all the way down to index 0.
            std::vector<km::Skmer<kuint>> skmers{
                km::Skmer<kuint>(kpair{0, 0}, 2, 0),
                km::Skmer<kuint>(kpair{0, 0}, 2, 0),
                km::Skmer<kuint>(kpair{0, 0}, 2, 0),
            };
            km::sortedlist::SortedVirtualSkmerList<kuint> m_list(k, m);
            m_list.add_list(skmers);

            // Search from the middle with lower bound 0: no hit anywhere.
            // Must return -1 rather than read out of bounds.
            const int64_t found = m_list.find_closest_valid_skmer(1, 0, 2, k - m);
            ASSERT_EQ(found, -1) << "no skmer has a valid kmer at position k-m";
        }

    }
}

// FULL LIST GENERATION
// generate_sorted_list_from_enumeration
TEST(SortedVirtualSkmerListTest, OneSkmerInList)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[1]{
    // P:    A   A   C    A
    // S:  A   C   C   _
        {0b0000010001011100U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 3, 2)
    };

    // EXPECTED
    const kpair expected_pairs[1]{
    // P:    A   A   C    A
    // S:  A   C   C   _
        {0b0000010001011100U, 0},
    };
    km::Skmer<kuint> expected_skmer(expected_pairs[0], 3, 2);

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), 1);

    ASSERT_EQ(expected_skmer.m_pref_size, m_output_list[0].m_pref_size);
    ASSERT_EQ(expected_skmer.m_suff_size, m_output_list[0].m_suff_size);
    ASSERT_EQ(expected_skmer.m_pair, m_output_list[0].m_pair);
}

TEST(SortedVirtualSkmerListTest, MergeTwoKmerInto1Right)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[2]{
    // P:    A   G   T    _
    // S:  C   T   T   _
        {0b0100101110101111U, 0},
    // P:    A   _   _    _
    // S:  C   T   T   G
        {0b0100101110111111U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 2, 2),
        km::Skmer<kuint>(input_skmers[1], 0, 3)
    };

    // EXPECTED
    const kpair expected_pairs[1]{
    // P:    A   G   T    _
    // S:  C   T   T   G
        {0b0100101110101111U, 0},
    };
    km::Skmer<kuint> expected_skmer(expected_pairs[0], 2, 3);

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), 1);

    ASSERT_EQ(expected_skmer.m_pref_size, m_output_list[0].m_pref_size);
    ASSERT_EQ(expected_skmer.m_suff_size, m_output_list[0].m_suff_size);
    ASSERT_EQ(expected_skmer.m_pair, m_output_list[0].m_pair);
}

TEST(SortedVirtualSkmerListTest, SortKmerSameColumn)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[3]{
    // P:    A   A   _    _
    // S:  C   C   A   _
        {0b0100010000111111U, 0},
    // P:    A   A   _    _
    // S:  A   C   C   _
        {0b0000010001111111U, 0},
    // P:    A   C   _    _
    // S:  A   C   C   _
        {0b0000010101111111U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 1, 2),
        km::Skmer<kuint>(input_skmers[1], 1, 2),
        km::Skmer<kuint>(input_skmers[2], 1, 2),
    };

    // EXPECTED
    const kpair expected_pairs[3]{
    // P:    A   A   _    _
    // S:  A   C   C   _
        {0b0000010001111111U, 0},
    // P:    A   C   _    _
    // S:  A   C   C   _
        {0b0000010101111111U, 0},
    // P:    A   A   _    _
    // S:  C   C   A   _
        {0b0100010000111111U, 0},
    };
    std::vector<km::Skmer<kuint>> expected_skmer{
        km::Skmer<kuint>(expected_pairs[0], 1, 2),
        km::Skmer<kuint>(expected_pairs[1], 1, 2),
        km::Skmer<kuint>(expected_pairs[2], 1, 2),
    };

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), expected_skmer.size());
    for (size_t i{0}; i < skmer_enumeration.size(); i++)
    {
        // std::cout << "Checking skmer #" << i << std::endl;
        ASSERT_EQ(expected_skmer[i].m_pref_size, m_output_list[i].m_pref_size);
        ASSERT_EQ(expected_skmer[i].m_suff_size, m_output_list[i].m_suff_size);
        ASSERT_EQ(expected_skmer[i].m_pair, m_output_list[i].m_pair);
    }
}

TEST(SortedVirtualSkmerListTest, SortAndCompactSuperKmers1)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[3]{
    // P:    A   A   _    _
    // S:  C   C   A   _
        {0b0100010000111111U, 0},
    // P:    A   A   C    _
    // S:  C   C   _   _
        {0b0100010011011111U, 0},
    // P:    T   T   C    _
    // S:  C   A   G   _
        {0b0110001011011111U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 1, 2),
        km::Skmer<kuint>(input_skmers[1], 2, 1),
        km::Skmer<kuint>(input_skmers[2], 2, 2),
    };

    //EXPECTED
    const kpair expected_pairs[2]{
    // P:    A   A   C    _
    // S:  C   C   A   _
        {0b0100010000011111U, 0},
    // P:    T   T   C    _
    // S:  C   A   G   _
        {0b0110001011011111U, 0},
    };
    std::vector<km::Skmer<kuint>> expected_skmer{
        km::Skmer<kuint>(expected_pairs[0], 2, 2),
        km::Skmer<kuint>(expected_pairs[1], 2, 2),
    };

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), expected_skmer.size());
    for (size_t i{0}; i < expected_skmer.size(); i++)
    {
        // std::cout << "Checking skmer #" << i << std::endl;
        ASSERT_EQ(expected_skmer[i].m_pref_size, m_output_list[i].m_pref_size);
        ASSERT_EQ(expected_skmer[i].m_suff_size, m_output_list[i].m_suff_size);
        ASSERT_EQ(expected_skmer[i].m_pair, m_output_list[i].m_pair);
    }
}


TEST(SortedVirtualSkmerListTest, SortAndCompactSuperKmers2)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{5};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[3]{
    // P:    A   C   T    _
    // S:  A   C   _   _
        {0b0000010111101111U, 0},
    // P:    A   C   T    G
    // S:  A   _   _   _
        {0b0000110111101111U, 0},
    // P:    A   C   _    _
    // S:  A   C   G   G
        {0b0000010111111111U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 2, 1),
        km::Skmer<kuint>(input_skmers[1], 3, 0),
        km::Skmer<kuint>(input_skmers[2], 1, 3),
    };

    //EXPECTED
    const kpair expected_pairs[1]{
    // P:    A   C   T    G
    // S:  A   C   G   G
        {0b0000010111101111U, 0},
    };
    std::vector<km::Skmer<kuint>> expected_skmer{
        km::Skmer<kuint>(expected_pairs[0], 3, 3)
    };

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), expected_skmer.size());

    for (size_t i{0}; i < expected_skmer.size(); i++)
    {
        // std::cout << "Checking skmer #" << i << std::endl;
        // std::cout << "Output has prefix: " << m_output_list[i].m_pref_size << "; suffix: " << m_output_list[i].m_suff_size << std::endl;
        ASSERT_EQ(expected_skmer[i].m_pair, m_output_list[i].m_pair);
        ASSERT_EQ(expected_skmer[i].m_pref_size, m_output_list[i].m_pref_size);
        ASSERT_EQ(expected_skmer[i].m_suff_size, m_output_list[i].m_suff_size);
    }
}

TEST(SortedVirtualSkmerListTest, SortAndCompactSuperKmers3)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[5]{
    // P:    A   C   _
    // S:  A   G   C
        {0b000011010111U, 0},
    // P:    C   G   A
    // S:  A   C   _
        {0b000101111100U, 0},
    // P:    A   C   C
    // S:  A   _   _
        {0b000011011101U, 0},
    // P:    C   C   _
    // S:  A   G   T
        {0b000111011011U, 0},
    // P:    C   _   _
    // S:  A   T   A
        {0b000110110011U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 1, 2),
        km::Skmer<kuint>(input_skmers[1], 2, 1),
        km::Skmer<kuint>(input_skmers[2], 2, 0),
        km::Skmer<kuint>(input_skmers[3], 1, 2),
        km::Skmer<kuint>(input_skmers[4], 0, 2)
    };

    //EXPECTED
    const kpair expected_pairs[4]{
    // P:    A   C   C
    // S:  A   G   C
        {0b000011010101U, 0},
    // P:    C   G   A
    // S:  A   C   _
        {0b000101111100U, 0},
    // P:    C   _   _
    // S:  A   T   A
        {0b000110110011U, 0},
    // P:    C   C   _
    // S:  A   G   T
        {0b000111011011U, 0},
    };
    std::vector<km::Skmer<kuint>> expected_skmer{
        km::Skmer<kuint>(expected_pairs[0], 2, 2),
        km::Skmer<kuint>(expected_pairs[1], 2, 1),
        km::Skmer<kuint>(expected_pairs[2], 0, 2),
        km::Skmer<kuint>(expected_pairs[3], 1, 2)
    };

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), expected_skmer.size());

    for (size_t i{0}; i < expected_skmer.size(); i++)
    {
        // std::cout << "Checking skmer #" << i << std::endl;
        // std::cout << "Output has prefix: " << m_output_list[i].m_pref_size << "; suffix: " << m_output_list[i].m_suff_size << std::endl;
        ASSERT_EQ(expected_skmer[i].m_pair, m_output_list[i].m_pair);
        ASSERT_EQ(expected_skmer[i].m_pref_size, m_output_list[i].m_pref_size);
        ASSERT_EQ(expected_skmer[i].m_suff_size, m_output_list[i].m_suff_size);
    }
}

TEST(SortedVirtualSkmerListTest, SortAndCompactSuperKmers4)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[2]{
    // P:    C   G   A
    // S:  A   _   _
        {0b000111111100U, 0},
    // P:    A   C   C
    // S:  A   _   _
        {0b000011011101U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 2, 0),
        km::Skmer<kuint>(input_skmers[1], 2, 0),
    };

    //EXPECTED
    const kpair expected_pairs[2]{
    // P:    A   C   C
    // S:  A   _   _
        {0b000011011101U, 0},
    // P:    C   G   A
    // S:  A   _   _
        {0b000111111100U, 0},
    };
    std::vector<km::Skmer<kuint>> expected_skmer{
        km::Skmer<kuint>(expected_pairs[0], 2, 0),
        km::Skmer<kuint>(expected_pairs[1], 2, 0),
    };

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), expected_skmer.size());

    for (size_t i{0}; i < expected_skmer.size(); i++)
    {
        // std::cout << "Checking skmer #" << i << std::endl;
        // std::cout << "Output has prefix: " << m_output_list[i].m_pref_size << "; suffix: " << m_output_list[i].m_suff_size << std::endl;
        ASSERT_EQ(expected_skmer[i].m_pair, m_output_list[i].m_pair);
        ASSERT_EQ(expected_skmer[i].m_pref_size, m_output_list[i].m_pref_size);
        ASSERT_EQ(expected_skmer[i].m_suff_size, m_output_list[i].m_suff_size);
    }
}

TEST(SortedVirtualSkmerListTest, EmptyInput)
{
    using kuint = uint16_t;
    // using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);


    std::vector<km::Skmer<kuint>> skmer_enumeration{};

    //EXPECTED EMPTY LIST OUT

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(),0);
}

//TEST(SortedVirtualSkmerListTest, AllKmersInvalid){}

// Giving as input the same k-mer twice, checking eforcement of uniqueness
TEST(SortedVirtualSkmerListTest, DuplicateKmers){
using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);

    //INPUT
    const kpair input_skmers[4]{
    // P:    C   G   A
    // S:  A   _   _
        {0b000111111100U, 0},
    // P:    A   C   C
    // S:  A   _   _
        {0b000011011101U, 0},
    // P:    A   C   C
    // S:  A   _   _
        {0b000011011101U, 0},
    // P:    C   G   A
    // S:  A   _   _
        {0b000111111100U, 0},
    };
    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_skmers[0], 2, 0),
        km::Skmer<kuint>(input_skmers[1], 2, 0),
        km::Skmer<kuint>(input_skmers[0], 2, 0),
        km::Skmer<kuint>(input_skmers[1], 2, 0),
    };

    //EXPECTED
    const kpair expected_pairs[2]{
    // P:    A   C   C
    // S:  A   _   _
        {0b000011011101U, 0},
    // P:    C   G   A
    // S:  A   _   _
        {0b000111111100U, 0},
    };
    std::vector<km::Skmer<kuint>> expected_skmer{
        km::Skmer<kuint>(expected_pairs[0], 2, 0),
        km::Skmer<kuint>(expected_pairs[1], 2, 0),
    };

    list.generate_sorted_list_from_enumeration(skmer_enumeration);

    std::vector<km::Skmer<kuint>> m_output_list = list.get_list();

    ASSERT_EQ(m_output_list.size(), expected_skmer.size());

    for (size_t i{0}; i < expected_skmer.size(); i++)
    {
        // std::cout << "Checking skmer #" << i << std::endl;
        // std::cout << "Output has prefix: " << m_output_list[i].m_pref_size << "; suffix: " << m_output_list[i].m_suff_size << std::endl;
        ASSERT_EQ(expected_skmer[i].m_pair, m_output_list[i].m_pair);
        ASSERT_EQ(expected_skmer[i].m_pref_size, m_output_list[i].m_pref_size);
        ASSERT_EQ(expected_skmer[i].m_suff_size, m_output_list[i].m_suff_size);
    }
}
// TEST(SortedVirtualSkmerListTest, MaximumKmerPosition)
// TEST(SortedVirtualSkmerListTest, MinimumKmerPosition)

// QUERY TESTS
// TEST(QueryTest, QueryEmptyList)
// TEST(QueryTest, QueryNotFound)
// TEST(QueryTest, QueryMultipleMatches)
// TEST(QueryTest, QueryBoundaryPositions)
// TEST(QueryTest, QueryBatchEmpty)
// TEST(QueryTest, QueryInvalidKmer)
TEST(QueryTest, SearchablePositions)
{
    uint64_t num_of_elements {6};
    int64_t mean {3};
    constexpr uint64_t max_array {32};
    uint8_t search_slots[max_array];
    std::fill_n(search_slots,num_of_elements,1);

    uint64_t no_search_slots[1] {1};
    for (const uint64_t el: no_search_slots){
        search_slots[el] = 0;
    }
    std::pair<int64_t,int64_t> binary_positions[max_array];
    binary_positions[0].first = 5;
    binary_positions[0].second = 10;
    binary_positions[1].first = 0;
    binary_positions[1].second = 5;
    binary_positions[2].first = 5;
    binary_positions[2].second = 10;
    binary_positions[3].first = 0;
    binary_positions[3].second = 5;
    binary_positions[4].first = 10;
    binary_positions[4].second = 20;
    binary_positions[5].first = 0;
    binary_positions[5].second = 5;



    uint64_t computed_searchable_positions[max_array];
    uint64_t num_searchable_positions {num_of_elements};
    num_searchable_positions = km::sortedlist::util::update_searchable_positions(mean, search_slots, binary_positions, computed_searchable_positions, num_of_elements);

    const std::vector<uint64_t> expected_result {3, 5};

    ASSERT_EQ(num_searchable_positions, expected_result.size());
    for(size_t i {0}; i < expected_result.size(); i++){
        ASSERT_EQ(computed_searchable_positions[i], expected_result[i]);
    }
}

TEST(QueryTest, QueryOneKmerinSkmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    const kpair input_pairs[4]{
    // P:    A   C   C
    // S:  A   G   C
        {0b000011010101U, 0},
    // P:    C   G   A
    // S:  A   C   _
        {0b000101111100U, 0},
    // P:    C   G   _
    // S:  A   T   A
        {0b000110110011U, 0},
    // P:    C   C   _
    // S:  A   G   T
        {0b000111011011U, 0},
    };

    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_pairs[0], 2, 2),
        km::Skmer<kuint>(input_pairs[1], 2, 1),
        km::Skmer<kuint>(input_pairs[2], 1, 2),
        km::Skmer<kuint>(input_pairs[3], 1, 2),
    };

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.add_list(skmer_enumeration);

    const kpair query_pairs[4]{
    // P:    C   G   A
    // S:  A   _   _
        {0b000111111100U, 0}, // IS IN
    // P:    C   G   _
    // S:  A   T   _
        {0b000110111111U, 0}, // IS IN
    // P:    C   _   _
    // S:  A   G   T
        {0b000111111011U, 0}, // IS IN
    // P:    A   A   _
    // S:  A   A   _
        {0b000000001111U, 0}, // IS NOT IN
    };

    std::vector<km::Skmer<kuint>> query_skmer{
        km::Skmer<kuint>(query_pairs[0], 2, 0),
        km::Skmer<kuint>(query_pairs[1], 1, 1),
        km::Skmer<kuint>(query_pairs[2], 0, 2),
        km::Skmer<kuint>(query_pairs[3], 1, 1),
    };

    std::vector<uint8_t> expected_result {1, 1, 1, 0};

    for(size_t i {0}; i < query_skmer.size(); i++){
        ASSERT_EQ(expected_result[i], list.query_skmer(query_skmer[i])[0]);
    }
}

TEST(QueryTest, QueryMultipleKmersInSkmer)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    const kpair input_pairs[4]{
    // P:    A   C   C
    // S:  A   G   C
        {0b000011010101U, 0},
    // P:    C   G   A
    // S:  A   C   _
        {0b000101111100U, 0},
    // P:    C   G   _
    // S:  A   T   A
        {0b000110110011U, 0},
    // P:    C   C   _
    // S:  A   G   T
        {0b000111011011U, 0},
    };

    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_pairs[0], 2, 2),
        km::Skmer<kuint>(input_pairs[1], 2, 1),
        km::Skmer<kuint>(input_pairs[2], 1, 2),
        km::Skmer<kuint>(input_pairs[3], 1, 2),
    };

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.add_list(skmer_enumeration);

    const kpair query_pairs[4]{
    // P:    C   G   A
    // S:  A   T   _
        {0b000110111100U, 0}, // IS IN
    // P:    C   T   _
    // S:  A   G   T
        {0b000111101011U, 0}, // FIRST NOT IN, SECOND NOT IN
    };

    std::vector<km::Skmer<kuint>> query_skmer{
        km::Skmer<kuint>(query_pairs[0], 2, 1),
        km::Skmer<kuint>(query_pairs[1], 1, 2),
    };

    std::vector<std::vector<uint8_t>> expected_result {{1, 1}, {0, 1}}; // 0 is false, 1 is true

    for(size_t i {0}; i < query_skmer.size(); i++){
        ASSERT_EQ(expected_result[i], list.query_skmer(query_skmer[i]));
    }
}

TEST(QueryTest, QueryMultipleKmersInSkmerBatch)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    const kpair input_pairs[4]{
    // P:    A   C   C
    // S:  A   G   C
        {0b000011010101U, 0},
    // P:    C   G   A
    // S:  A   C   _
        {0b000101111100U, 0},
    // P:    C   G   _
    // S:  A   T   A
        {0b000110110011U, 0},
    // P:    C   C   _
    // S:  A   G   T
        {0b000111011011U, 0},
    };

    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_pairs[0], 2, 2),
        km::Skmer<kuint>(input_pairs[1], 2, 1),
        km::Skmer<kuint>(input_pairs[2], 1, 2),
        km::Skmer<kuint>(input_pairs[3], 1, 2),
    };

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.add_list(skmer_enumeration);

    const kpair query_pairs[4]{
    // P:    C   G   A
    // S:  A   T   _
        {0b000110111100U, 0}, // IS IN
    // P:    C   T   _
    // S:  A   G   T
        {0b000111101011U, 0}, // FIRST NOT IN, SECOND NOT IN
    // P:    C   G   A
    // S:  A   T   _
        {0b000110111100U, 0}, // IS IN
    // P:    C   T   _
    // S:  A   G   T
        {0b000111101011U, 0}, // FIRST NOT IN, SECOND NOT IN
    };

    std::vector<km::Skmer<kuint>> query_skmer{
        km::Skmer<kuint>(query_pairs[0], 2, 1),
        km::Skmer<kuint>(query_pairs[1], 1, 2),
    };

    std::vector<std::vector<uint8_t>> expected_result {{1, 1}, {0, 1}};

    std::vector<std::vector<uint8_t>> computed_result = list.query_skmer_batch(query_skmer);
    ASSERT_EQ(expected_result.size(), computed_result.size());
    for(size_t i {0}; i < query_skmer.size(); i++){
        ASSERT_EQ(expected_result[i], computed_result[i]);
    }
}

// TOTAL QUERY TEST — captures the output and validates structure +
// self-consistency with query_skmer_batch on the same file-derived skmers.
TEST(QueryTest, QueryOutputWorking)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    const kpair input_pairs[4]{
    // P:    A   C   C
    // S:  A   G   C
        {0b000011010101U, 0},
    // P:    C   G   A
    // S:  A   C   _
        {0b000101111100U, 0},
    // P:    C   G   _
    // S:  A   T   A
        {0b000110110011U, 0},
    // P:    C   C   _
    // S:  A   G   T
        {0b000111011011U, 0},
    };

    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_pairs[0], 2, 2),
        km::Skmer<kuint>(input_pairs[1], 2, 1),
        km::Skmer<kuint>(input_pairs[2], 1, 2),
        km::Skmer<kuint>(input_pairs[3], 1, 2),
    };

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.add_list(skmer_enumeration);

    std::string filename{test_data("seq2.fa")};

    std::ostringstream captured;
    list.query(filename, captured);
    const std::string actual = captured.str();

    // Compute the reference output via the lower-level API on the same file.
    km::SkmerManipulator<kuint> manip{k, m};
    km::FileSkmerator<kuint> file_skmerator{manip, filename};
    std::vector<km::Skmer<kuint>> file_skmers;
    for (km::Skmer<kuint> const& s : file_skmerator) {
        file_skmers.push_back(s);
    }
    ASSERT_FALSE(file_skmers.empty()) << "test data file produced no skmers";

    auto batch_results = list.query_skmer_batch(file_skmers);
    std::ostringstream expected;
    km::sortedlist::util::print_query_results(batch_results, expected);

    ASSERT_FALSE(actual.empty()) << "query() wrote nothing to its ostream";
    ASSERT_EQ(actual, expected.str())
        << "query(filename, os) output diverges from query_skmer_batch + print_query_results";

    // Format sanity: every character must be '0', '1', ',', or '\n'.
    for (char c : actual) {
        ASSERT_TRUE(c == '0' || c == '1' || c == ',' || c == '\n')
            << "unexpected character in query output: 0x" << std::hex << (int)c;
    }
}

// ROUND-TRIP TEST — save then load the list and verify query_skmer returns
// identical results on both the in-memory and disk-roundtripped lists.
TEST(QueryTest, QuerySerializationRoundTrip)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    const kpair input_pairs[4]{
        {0b000011010101U, 0},
        {0b000101111100U, 0},
        {0b000110110011U, 0},
        {0b000111011011U, 0},
    };

    std::vector<km::Skmer<kuint>> skmer_enumeration{
        km::Skmer<kuint>(input_pairs[0], 2, 2),
        km::Skmer<kuint>(input_pairs[1], 2, 1),
        km::Skmer<kuint>(input_pairs[2], 1, 2),
        km::Skmer<kuint>(input_pairs[3], 1, 2),
    };

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.add_list(skmer_enumeration);

    const std::string tmp_path = ::testing::TempDir() + "roundtrip_query.bin";
    km::sortedlist::VirtualSkmerSerializer<kuint>::save(list, tmp_path);
    auto loaded = km::sortedlist::VirtualSkmerSerializer<kuint>::load(tmp_path);

    // Re-use the queries from QueryMultipleKmersInSkmer (mix of present + absent).
    const kpair query_pairs[2]{
        {0b000110111100U, 0}, // both present
        {0b000111101011U, 0}, // first absent, second present
    };
    std::vector<km::Skmer<kuint>> queries{
        km::Skmer<kuint>(query_pairs[0], 2, 1),
        km::Skmer<kuint>(query_pairs[1], 1, 2),
    };

    for (km::Skmer<kuint> const& q : queries) {
        auto in_memory = list.query_skmer(q);
        auto on_disk = loaded.query_skmer(q);
        ASSERT_EQ(in_memory, on_disk)
            << "round-tripped list disagrees with in-memory list";
        ASSERT_FALSE(in_memory.empty()) << "query returned no results";
    }
}

// CANONICAL TEST — a list built from a sequence S must yield the same
// matches when queried with skmers from S itself and from reverse-complement(S),
// because canonical encoding makes both produce the same skmers.
TEST(QueryTest, QueryCanonicalForwardReverse)
{
    using kuint = uint16_t;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    // Non-palindromic sequence with multiple distinct minimizers.
    std::string fwd_seq{"ACGTACGTACGTACGT"};

    km::SkmerManipulator<kuint> build_manip{k, m};
    std::string build_seq = fwd_seq;
    km::SeqSkmerator<kuint> build_rator{build_manip, build_seq};
    std::vector<km::Skmer<kuint>> build_skmers;
    for (km::Skmer<kuint> const& s : build_rator) {
        build_skmers.push_back(s);
    }
    ASSERT_FALSE(build_skmers.empty()) << "sequence produced no skmers";

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(build_skmers);

    // Forward query: every kmer is in the list by construction.
    km::SkmerManipulator<kuint> fwd_manip{k, m};
    std::string fwd_query_seq = fwd_seq;
    km::SeqSkmerator<kuint> fwd_rator{fwd_manip, fwd_query_seq};
    std::vector<km::Skmer<kuint>> fwd_skmers;
    for (km::Skmer<kuint> const& s : fwd_rator) {
        fwd_skmers.push_back(s);
    }
    auto fwd_results = list.query_skmer_batch(fwd_skmers);
    for (size_t i = 0; i < fwd_results.size(); i++) {
        for (size_t j = 0; j < fwd_results[i].size(); j++) {
            ASSERT_EQ(fwd_results[i][j], 1u)
                << "forward query missed kmer at skmer " << i << " position " << j;
        }
    }

    // Reverse-complement query: canonical encoding must yield the same matches.
    std::string rc_seq;
    rc_seq.reserve(fwd_seq.size());
    for (auto it = fwd_seq.rbegin(); it != fwd_seq.rend(); ++it) {
        switch (*it) {
            case 'A': rc_seq += 'T'; break;
            case 'C': rc_seq += 'G'; break;
            case 'G': rc_seq += 'C'; break;
            case 'T': rc_seq += 'A'; break;
            default: FAIL() << "unexpected nucleotide in input";
        }
    }
    km::SkmerManipulator<kuint> rc_manip{k, m};
    km::SeqSkmerator<kuint> rc_rator{rc_manip, rc_seq};
    std::vector<km::Skmer<kuint>> rc_skmers;
    for (km::Skmer<kuint> const& s : rc_rator) {
        rc_skmers.push_back(s);
    }
    auto rc_results = list.query_skmer_batch(rc_skmers);
    for (size_t i = 0; i < rc_results.size(); i++) {
        for (size_t j = 0; j < rc_results[i].size(); j++) {
            ASSERT_EQ(rc_results[i][j], 1u)
                << "reverse-complement query missed kmer at skmer " << i << " position " << j;
        }
    }
}

// ALL-ABSENT TEST — every queried kmer position must be eliminated by
// binary-search exhaustion (no match path), exercising the termination
// branch at lines 404–407 of VirtualSkmer.hpp.
TEST(QueryTest, QueryAllKmersAbsent)
{
    using kuint = uint16_t;
    using kpair = km::Skmer<kuint>::pair;

    constexpr uint64_t k{4};
    constexpr uint64_t m{2};

    // Single-skmer list containing only AAAA kmers (pair = 0 => all A).
    const kpair list_pair{0b000000000000U, 0};
    std::vector<km::Skmer<kuint>> list_skmers{
        km::Skmer<kuint>(list_pair, 2, 2),
    };
    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.add_list(list_skmers);

    // Query skmer whose valid kmer positions hold G nucleotides (0b11) —
    // disjoint from any AAAA kmer in the list.
    const kpair query_pair{0b000011111111U, 0};
    km::Skmer<kuint> query{query_pair, 1, 2};

    std::vector<uint8_t> result = list.query_skmer(query);

    ASSERT_EQ(result.size(), 2u) << "expected 2 valid kmer positions for pref=1 suff=2";
    ASSERT_EQ(result[0], 0u) << "kmer at position 0 should be absent";
    ASSERT_EQ(result[1], 0u) << "kmer at position 1 should be absent";
}

// SET OPERATION TESTS

// INTERSECTION TESTS
// TEST(SetOperationTest, IntersectionBasic)
// TEST(SetOperationTest, IntersectionEmpty)
// TEST(SetOperationTest, IntersectionDisjoint)


// UNION TESTS
/// TEST(SetOperationTest, UnionBasic)
// TEST(SetOperationTest, UnionEmpty)

// DIFF TESTS
// TEST(SetOperationTest, DifferenceBasic)
// TEST(SetOperationTest, DifferenceEmpty)
// TEST(SetOperationTest, DifferenceIdentical)