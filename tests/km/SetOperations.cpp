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

inline km::Skmer<kuint> make_sk(uint16_t v, uint16_t pref, uint16_t suff) {
    return km::Skmer<kuint>(kpair{v, 0}, pref, suff);
}

km::Skmer<kuint> sk_one()
    //       Prefix:      A   C   G   T
    //       Suffix:    T   C   T   C
    { return make_sk(0b1000010110110110U, 3, 3); }
km::Skmer<kuint> sk_two()
     //        Prefix:    C   C   G   T
    //      Suffix:    G   C   T   C
    { return make_sk(0b1101010110110110U, 3, 3); }
km::Skmer<kuint> sk_three()
    //       Prefix:      T   C   G   _
    //       Suffix:    T   C   _   _
    { return make_sk(0b1010010111111111U, 2, 1); }

km::Skmer<kuint> sk_four()
    //       Prefix:      A   A   G   T
    //       Suffix:    A   A   T   C
    { return make_sk(0b0000000010110110U, 3, 3); }
km::Skmer<kuint> sk_five()
     //        Prefix:    A   C   G   T
    //      Suffix:    A   C   T   C
    { return make_sk(0b0000010110110110U, 3, 3); }
km::Skmer<kuint> sk_six()
    //       Prefix:      A   G   G   _
    //       Suffix:    A   G   _   _
    { return make_sk(0b0000111111111111U, 2, 1); }


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
// 2. Disjoint minimizers  (A:{Y,Z}  B:{X}), different minimizer (fast path)
// ------------------------------------------------------------------
TEST(SetOps, Union_Disjoint_diff_minimizer) {
    auto sk_list_one = make_list({sk_one(), sk_two()});
    auto sk_list_two = make_list({sk_three()});

    auto res = km::sortedlist::setops::set_union(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), sk_list_one.size() + sk_list_two.size() );
    expect_skmer_eq(res.get_list()[0], sk_one());
    expect_skmer_eq(res.get_list()[1], sk_three());
    expect_skmer_eq(res.get_list()[2], sk_two());
}

TEST(SetOps, Intersection_Disjoint_diff_minimizer) {
    auto sk_list_one = make_list({sk_one(), sk_two()});
    auto sk_list_two = make_list({sk_three()});

    auto res = km::sortedlist::setops::set_intersection(sk_list_one, sk_list_two);
    EXPECT_EQ(res.size(), 0);
}

TEST(SetOps, Diff_Disjoint_diff_minimizer) {
    auto sk_list_one = make_list({sk_one(), sk_two()});
    auto sk_list_two = make_list({sk_three()});

    auto res = km::sortedlist::setops::set_diff(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), sk_list_one.size());
    expect_skmer_eq(res.get_list()[0], sk_one());
    expect_skmer_eq(res.get_list()[1], sk_two());
}
// ------------------------------------------------------------------
// 3. Disjoint minimizers  (A:{Y,Z}  B:{X}), same minimizer (slow path)
// ------------------------------------------------------------------
TEST(SetOps, Union_Disjoint_same_minimizer) {
    auto sk_list_one = make_list({sk_four(), sk_six()});
    auto sk_list_two = make_list({sk_five()});

    auto res = km::sortedlist::setops::set_union(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), sk_list_one.size() + sk_list_two.size() );
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_five());
    expect_skmer_eq(res.get_list()[2], sk_six());
}

TEST(SetOps, Intersection_Disjoint_same_minimizer) {
    auto sk_list_one = make_list({sk_four(), sk_five()});
    auto sk_list_two = make_list({sk_four(), sk_six()});

    auto res = km::sortedlist::setops::set_intersection(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), 1);
    expect_skmer_eq(res.get_list()[0], sk_four());
}

TEST(SetOps, Diff_Disjoint_same_minimizer) {
    auto sk_list_one = make_list({sk_five(), sk_six()});
    auto sk_list_two = make_list({sk_four()});

    auto res = km::sortedlist::setops::set_diff(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), sk_list_one.size());
    expect_skmer_eq(res.get_list()[0], sk_five());
    expect_skmer_eq(res.get_list()[1], sk_six());
}

// ------------------------------------------------------------------
// 4. Identical lists  (exercises same-minimizer / full-dedup path)
// ------------------------------------------------------------------
TEST(SetOps, Union_Identical) {
    auto sk_list_one = make_list({sk_four(), sk_five()});
    auto sk_list_two = make_list({sk_five(), sk_four()});

    auto res = km::sortedlist::setops::set_union(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), sk_list_one.size());                 // deduplicated
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_five());
}

TEST(SetOps, Intersection_Identical) {
    auto sk_list_one = make_list({sk_four(), sk_one()});

    auto res = km::sortedlist::setops::set_intersection(sk_list_one, sk_list_one);
    ASSERT_EQ(res.size(), sk_list_one.size());
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_one());
}

TEST(SetOps, Diff_Identical) {
    auto sk_list_one = make_list({sk_four(), sk_one()});

    auto res = km::sortedlist::setops::set_diff(sk_list_one, sk_list_one);
    EXPECT_EQ(res.size(), 0);
}

// ------------------------------------------------------------------
// 5. Subset overlap  (A:{Y,Z}  B:{Y})
//    Y falls into the "same minimizer" slow path,
//    Z stays on the fast merge-scan path.
// ------------------------------------------------------------------
TEST(SetOps, Union_Subset) {
    auto sk_list_one = make_list({sk_two(), sk_five()});
    auto sk_list_two = make_list({sk_two()});

    auto res = km::sortedlist::setops::set_union(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), 2);
    expect_skmer_eq(res.get_list()[0], sk_five());
    expect_skmer_eq(res.get_list()[1], sk_two());
}

TEST(SetOps, Intersection_Subset) {
    auto sk_list_one = make_list({sk_two(), sk_five()});
    auto sk_list_two = make_list({sk_two()});

    auto res = km::sortedlist::setops::set_intersection(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), 1);
    expect_skmer_eq(res.get_list()[0], sk_two());
}

TEST(SetOps, Diff_Subset) {
    auto sk_list_one = make_list({sk_two(), sk_five()});
    auto sk_list_two = make_list({sk_two()});

    auto res = km::sortedlist::setops::set_diff(sk_list_one, sk_list_two);
    ASSERT_EQ(res.size(), 1);
    expect_skmer_eq(res.get_list()[0], sk_five());
}

// ------------------------------------------------------------------
// 6. File source round-trips
// ------------------------------------------------------------------
class SetOpsFile : public ::testing::Test {
protected:
    static std::string tmpname(const std::string& tag) {
        return std::string("vsk_test_") + tag + ".bin";
    }
    static void save(const std::string& path,
                     const km::sortedlist::SortedVirtualSkmerList<kuint>& list) {
        km::sortedlist::VirtualSkmerSerializer<kuint>::save(list, path);
    }
    static void cleanup(const std::string& path) { std::remove(path.c_str()); }
};


// mem + file
TEST_F(SetOpsFile, Union_MemWithFile) {
    auto a = make_list({sk_four(), sk_five(), sk_one(), sk_three()});
    auto b = make_list({sk_four(), sk_six(), sk_one(), sk_two()});
    std::string f = tmpname("b");
    save(f, b);

    auto res = km::sortedlist::setops::set_union(a, f);
    cleanup(f);

    ASSERT_EQ(res.size(), 6);
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_five());
    expect_skmer_eq(res.get_list()[2], sk_six());
    expect_skmer_eq(res.get_list()[3], sk_one());
    expect_skmer_eq(res.get_list()[4], sk_three());
    expect_skmer_eq(res.get_list()[5], sk_two());
}

TEST_F(SetOpsFile, Intersection_MemWithFile) {
    auto a = make_list({sk_four(), sk_five(), sk_one(), sk_three()});
    auto b = make_list({sk_four(), sk_six(), sk_one(), sk_two()});
    std::string f = tmpname("b");
    save(f, b);

    auto res = km::sortedlist::setops::set_intersection(a, f);
    cleanup(f);

    ASSERT_EQ(res.size(), 2);
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_one());
}

TEST_F(SetOpsFile, Diff_MemWithFile) {
    auto a = make_list({sk_four(), sk_five(), sk_one(), sk_three()});
    auto b = make_list({sk_four(), sk_six(), sk_one(), sk_two()});
    std::string f = tmpname("b");
    save(f, b);

    auto res = km::sortedlist::setops::set_diff(a, f);
    cleanup(f);

    ASSERT_EQ(res.size(), 2);
    expect_skmer_eq(res.get_list()[0], sk_five());
    expect_skmer_eq(res.get_list()[1], sk_three());
}

// file + mem
TEST_F(SetOpsFile, Union_FileWithMem) {
    auto a = make_list({sk_four(), sk_five(), sk_one(), sk_three()});
    auto b = make_list({sk_four(), sk_six(), sk_one(), sk_two()});
    std::string f = tmpname("a");
    save(f, a);

    auto res = km::sortedlist::setops::set_union(f, b);
    cleanup(f);

    ASSERT_EQ(res.size(), 6);
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_five());
    expect_skmer_eq(res.get_list()[2], sk_six());
    expect_skmer_eq(res.get_list()[3], sk_one());
    expect_skmer_eq(res.get_list()[4], sk_three());
    expect_skmer_eq(res.get_list()[5], sk_two());
}

TEST_F(SetOpsFile, Intersection_FileWithMem) {
    auto a = make_list({sk_four(), sk_five(), sk_one(), sk_three()});
    auto b = make_list({sk_four(), sk_six(), sk_one(), sk_two()});
    std::string f = tmpname("b");
    save(f, b);

    auto res = km::sortedlist::setops::set_intersection(f, a);
    cleanup(f);

    ASSERT_EQ(res.size(), 2);
    expect_skmer_eq(res.get_list()[0], sk_four());
    expect_skmer_eq(res.get_list()[1], sk_one());

}

TEST_F(SetOpsFile, Diff_FileWithMem) {
    auto a = make_list({sk_four(), sk_five(), sk_one(), sk_three()});
    auto b = make_list({sk_four(), sk_six(), sk_one(), sk_two()});
    std::string f = tmpname("b");
    save(f, b);

    auto res = km::sortedlist::setops::set_diff(f, a);
    cleanup(f);

    ASSERT_EQ(res.size(), 2);
    expect_skmer_eq(res.get_list()[0], sk_six());
    expect_skmer_eq(res.get_list()[1], sk_two());

}

