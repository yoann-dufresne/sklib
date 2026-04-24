#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <array>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/SetOperations.hpp>

using namespace std;
using kuint = uint16_t;
using kpair = km::Skmer<kuint>::pair;

constexpr uint64_t K = 5;
constexpr uint64_t M = 2;

km::sortedlist::SortedVirtualSkmerList<kuint>
make_list(const std::vector<km::Skmer<kuint>>& enumeration) {
    km::sortedlist::SortedVirtualSkmerList<kuint> list(K, M);
    list.generate_sorted_list_from_enumeration(enumeration);
    return list;
}

void expect_skmer_eq(const km::Skmer<kuint>& a,
                            const km::Skmer<kuint>& b) {
    EXPECT_EQ(a.m_pref_size, b.m_pref_size);
    EXPECT_EQ(a.m_suff_size, b.m_suff_size);
    EXPECT_EQ(a.m_pair,      b.m_pair);
}

inline km::Skmer<kuint> make_sk(uint64_t v, uint16_t pref, uint16_t suff) {
    return km::Skmer<kuint>(kpair{v, 0}, pref, suff);
}

km::Skmer<kuint> sk_one()
    //       Prefix:      A   C   G   T
    //       Suffix:    A   C   T   C
    { return make_sk(0b0000010110110110U, 3, 4); }
km::Skmer<kuint> sk_two()
     //        Prefix:    A   C   G   T
    //      Suffix:    A   C   T   C
    { return make_sk(0b0000010110110110U, 3, 3); }
km::Skmer<kuint> sk_three()
    //       Prefix:      A   C   G   _
    //       Suffix:    A   C   _   _
    { return make_sk(0b0000010111111111U, 2, 1); }


// ------------------------------------------------------------------
// 1. Empty-list handling
// ------------------------------------------------------------------
TEST(SetOps, Union_EmptyWithNonEmpty) {

    //                    Prefix:      A   C   G   T
    //                    Suffix:    A   C   T   C

    auto empty_list = make_list({});
    auto one_sk_list = make_list({sk_one(), sk_two(), sk_three()});

    auto res = km::sortedlist::setops::set_union(empty_list, one_sk_list);
    ASSERT_EQ(res.size(), one_sk_list.size());

    expect_skmer_eq(res.get_list()[0], one_sk_list.get_list()[0]);
}

TEST(SetOps, Intersection_EmptyWithNonEmpty) {
    auto empty_list = make_list({});
    auto one_sk_list = make_list({sk_one(), sk_two(), sk_three()});


    auto res = km::sortedlist::setops::set_intersection(empty_list, one_sk_list);
    EXPECT_EQ(res.size(), 0);
}

TEST(SetOps, Diff_EmptyMinusNonEmpty) {

    auto empty_list = make_list({});
    auto one_sk_list = make_list({sk_one(), sk_two(), sk_three()});

    auto res = km::sortedlist::setops::set_diff(empty_list, one_sk_list);
    EXPECT_EQ(res.size(), 0);
}

TEST(SetOps, Diff_NonEmptyMinusEmpty) {
    auto empty_list = make_list({});
    auto one_sk_list = make_list({sk_one(), sk_two(), sk_three()});

    auto res = km::sortedlist::setops::set_diff(one_sk_list, empty_list);
    ASSERT_EQ(res.size(), one_sk_list.size());
    expect_skmer_eq(res.get_list()[0], one_sk_list.get_list()[0]);

}

// ------------------------------------------------------------------
// 2. Disjoint minimizers  (A:{Y,Z}  B:{X})  for union, intersection, diff
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// 3. Identical lists  (exercises same-minimizer / full-dedup path)
// ------------------------------------------------------------------


// ------------------------------------------------------------------
// 4. Subset overlap  (A:{Y,Z}  B:{Y})
//    Y falls into the "same minimizer" slow path,
//    Z stays on the fast merge-scan path.
// ------------------------------------------------------------------
