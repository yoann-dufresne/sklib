#include <gtest/gtest.h>
#include <iostream>

#include <dbg/datastruct/RSQF.hpp>
#include <dbg/datastruct/quotienting.hpp>

using namespace std;
using namespace dbglib;

static constexpr uint64_t q1 {6};
static constexpr uint64_t r1 {6};

// class QuotientFilter<q, r, LeftQuotienting>;

class QuotientFilter_6_6_Fixture : public QuotientFilter<q1, r1>, public testing::Test
{};

static constexpr uint64_t q2 {8};
static constexpr uint64_t r2 {3};

class QuotientFilter_8_3_Fixture : public QuotientFilter<q2, r2>, public testing::Test
{};


TEST_F(QuotientFilter_6_6_Fixture, init)
{
    ASSERT_EQ(std::size(m_rests), 1);
}

TEST_F(QuotientFilter_8_3_Fixture, init)
{
    ASSERT_EQ(std::size(m_rests), 4);
}

TEST(RSQF_insert, insert_free_space)
{
    QuotientFilter<7, 5> filter{};

    // Prepare the example to insert
    const uint64_t quotient {2};
    const uint64_t rest {0b11011};
    const uint64_t val {(quotient << 5) | rest};

    // Perform insertion
    const auto qr {LeftQuotienting::compute<7,5>(val)};
    filter.insert_in_free_space(qr);

    // Verify metadata
    ASSERT_EQ(filter.m_occupied.get(quotient), true);
    ASSERT_EQ(filter.m_runend.get(quotient), true);

    // Verify the rest value inserted
    const uint64_t block_idx {quotient / 64};
    ASSERT_EQ(filter.m_rests[block_idx].get(quotient % 64), rest);
}

TEST(RSQF_insert_and_shift, single_block)
{
    QuotientFilter<7, 5> filter{};

    const uint64_t quotient {2};
    const uint64_t rest1 {0b11011};
    const uint64_t val1 {(quotient << 5) | rest1};

    const auto qr {LeftQuotienting::compute<7, 5>(val1)};
    filter.insert_in_free_space(qr);

    // --- Insert at the beginning of a run ---
    const uint64_t rest2 {0b01000};
    const auto val2 {LeftQuotienting::compute<7, 5>((quotient << 5) | rest2)};
    // const uint64_t val2 {(quotient << 5) | rest2};
    filter.insert_and_shift(val2, quotient, quotient+1);

    // Verification of the run start
    ASSERT_EQ(filter.m_occupied.get(quotient), true);
    ASSERT_EQ(filter.m_occupied.get(quotient+1), false);
    // Verification of the run end
    ASSERT_EQ(filter.m_runend.get(quotient), false);
    ASSERT_EQ(filter.m_runend.get(quotient+1), true);
    // Are the rests at the right place ?
    const uint64_t block_idx {quotient / 64};
    ASSERT_EQ(filter.m_rests[block_idx].get(quotient % 64), rest2);
    ASSERT_EQ(filter.m_rests[block_idx].get((quotient + 1) % 64), rest1);

    // --- Insert at the beginning of a run and shift 2 values ---
    const uint64_t rest3 {0b00001};
    const auto val3 {LeftQuotienting::compute<7, 5>((quotient << 5) | rest3)};
    // const uint64_t val2 {(quotient << 5) | rest2};
    filter.insert_and_shift(val3, quotient, quotient+2);

    // Verification of the run start
    ASSERT_EQ(filter.m_occupied.get(quotient), true);
    // Verification of the run end
    ASSERT_EQ(filter.m_runend.get(quotient+1), false);
    ASSERT_EQ(filter.m_runend.get(quotient+2), true);
    // Are the rests at the right place ?
    ASSERT_EQ(filter.m_rests[block_idx].get(quotient % 64), rest3);
    ASSERT_EQ(filter.m_rests[block_idx].get((quotient + 1) % 64), rest2);
    ASSERT_EQ(filter.m_rests[block_idx].get((quotient + 2) % 64), rest1);

    // --- Insert at the end of a run ---
    const uint64_t rest4 {0b11111};
    const auto val4 {LeftQuotienting::compute<7, 5>((quotient << 5) | rest4)};
    // const uint64_t val2 {(quotient << 5) | rest2};
    filter.insert_and_shift(val4, quotient+3, quotient+3);

    // Verification of the run end
    ASSERT_EQ(filter.m_runend.get(quotient+2), false);
    ASSERT_EQ(filter.m_runend.get(quotient+3), true);
    // Are the rests at the right place ?
    ASSERT_EQ(filter.m_rests[block_idx].get((quotient+3) % 64), rest4);
}



TEST(RSQF_first_unused, in_block_no_offset)
{
    QuotientFilter<7, 5> filter{};
    // Filter filled at quotients 2/3/4                     qqqrrrrr
    filter.insert_in_free_space(LeftQuotienting::compute<7, 5>(0b01010001));
    filter.insert_in_free_space(LeftQuotienting::compute<7, 5>(0b01110001));
    filter.insert_in_free_space(LeftQuotienting::compute<7, 5>(0b10010001));

    ASSERT_EQ(filter.first_unused_slot(1), 1);
    ASSERT_EQ(filter.first_unused_slot(2), 5);
    ASSERT_EQ(filter.first_unused_slot(3), 5);
    ASSERT_EQ(filter.first_unused_slot(4), 5);
    ASSERT_EQ(filter.first_unused_slot(5), 5);
    ASSERT_EQ(filter.first_unused_slot(6), 6);
}

TEST(RSQF_first_unused, in_block_offset)
{
    QuotientFilter<7, 5> filter{};
    // Filter filled at quotients 2/3/4                     qqqrrrrr
    filter.insert_in_free_space(LeftQuotienting::compute<7, 5>(0b01010001));
    filter.insert_in_free_space(LeftQuotienting::compute<7, 5>(0b01110001));
    filter.insert_in_free_space(LeftQuotienting::compute<7, 5>(0b10010001));

    // Faking an offset
    filter.m_offsets[0] = 3;
    filter.m_runend.set(10);

    for (uint64_t i{0} ; i<=11 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), 11);
    for (uint64_t i{12} ; i<64 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), i);
}


TEST(RSQF_first_unused, out_of_block_offset)
{
    QuotientFilter<7, 5> filter{};

    // Faking a multiblock offset
    filter.m_offsets[0] = 70;
    filter.m_offsets[1] = 6;
    filter.m_runend.set(69);

    for (uint64_t i{0} ; i<=70 ; i++)
    {
        ASSERT_EQ(filter.first_unused_slot(i), 70);
    }
    for (uint64_t i{71} ; i<100 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), i);
}

TEST(RSQF_first_unused, out_of_block_no_offset)
{
    QuotientFilter<7, 5> filter{};

    // Faking a long run spreading on 2 blocks
    filter.m_occupied.set(10);
    filter.m_offsets[1] = 7;
    filter.m_runend.set(70);

    for (uint64_t i{0} ; i<10 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), i);
    for (uint64_t i{10} ; i<=71 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), 71);
    for (uint64_t i{72} ; i<100 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), i);
}

TEST(RSQF_first_unused, out_of_block_toric)
{
    QuotientFilter<7, 5> filter{};

    // Faking a long run spreading on 2 blocks
    filter.m_occupied.set(75);
    filter.m_offsets[0] = 11;
    filter.m_runend.set(10);

    for (uint64_t i{64} ; i<75 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), i);
    for (uint64_t i{75} ; i<128 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), 11);
    for (uint64_t i{0} ; i<=11 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), 11);
    for (uint64_t i{12} ; i<64 ; i++)
        ASSERT_EQ(filter.first_unused_slot(i), i);
}

TEST(RSQF_INSERT_POSITION, inside_block_insert)
{
    QuotientFilter<7, 5> filter{};
    filter.m_occupied.set(10);
    filter.m_runend.set(12);
    filter.m_rests[0].set(10, 3);
    filter.m_rests[0].set(11, 5);
    filter.m_rests[0].set(12, 7);

    QR<7, 5> element_before {9, 4};
    QR<7, 5> element_after {11, 4};

    ASSERT_EQ(filter.compute_insert_position(element_before), 9);
    ASSERT_EQ(filter.compute_insert_position(element_after), 13);

    QR<7, 5> element2 {10, 2};
    ASSERT_EQ(filter.compute_insert_position(element2), 10);
    QR<7, 5> element4 {10, 4};
    ASSERT_EQ(filter.compute_insert_position(element4), 11);
    QR<7, 5> element6 {10, 6};
    ASSERT_EQ(filter.compute_insert_position(element6), 12);
    QR<7, 5> element8 {10, 8};
    ASSERT_EQ(filter.compute_insert_position(element8), 13);
}

TEST(RSQF_INSERT_POSITION, shifted_next_block)
{
    QuotientFilter<7, 5> filter{};
    filter.m_occupied.set(61);
    filter.m_runend.set(63);
    filter.m_rests[0].set(61, 3);
    filter.m_rests[0].set(62, 5);
    filter.m_rests[0].set(63, 7);

    QR<7, 5> element_after {62, 4};
    ASSERT_EQ(filter.compute_insert_position(element_after), 64);

    QR<7, 5> element6 {61, 6};
    ASSERT_EQ(filter.compute_insert_position(element6), 63);
    QR<7, 5> element8 {61, 8};
    ASSERT_EQ(filter.compute_insert_position(element8), 64);
}


TEST(RSQF_full_insert, single_block)
{
    QuotientFilter<7, 5> filter{};

    for (uint64_t i{3} ; i<=7 ; i+=2)
    {
        const uint64_t element {(10UL << 5) | i};
        filter.insert(element);
    }

    ASSERT_EQ(filter.m_rests[0].get(9), 0);
    ASSERT_EQ(filter.m_rests[0].get(10), 3);
    ASSERT_EQ(filter.m_rests[0].get(11), 5);
    ASSERT_EQ(filter.m_rests[0].get(12), 7);
    ASSERT_EQ(filter.m_rests[0].get(13), 0);
}


TEST(RSQF_full_insert, block_overflow)
{
    QuotientFilter<7, 5> filter{};

    for (uint64_t i{3} ; i<=7 ; i+=2)
    {
        const uint64_t element {(62UL << 5) | i};
        filter.insert(element);
    }

    ASSERT_EQ(filter.m_rests[0].get(61), 0);
    ASSERT_EQ(filter.m_rests[0].get(62), 3);
    ASSERT_EQ(filter.m_rests[0].get(63), 5);
    ASSERT_EQ(filter.m_rests[1].get(0), 7);
    ASSERT_EQ(filter.m_rests[1].get(1), 0);

    ASSERT_EQ(filter.m_offsets[1], 1);
}


TEST(RSQF_full_insert, toric)
{
    QuotientFilter<7, 5> filter{};

    for (uint64_t i{3} ; i<=7 ; i+=2)
    {
        const uint64_t element {(127UL << 5) | i};
        filter.insert(element);
    }

    ASSERT_EQ(filter.m_rests[1].get(62), 0);
    ASSERT_EQ(filter.m_rests[1].get(63), 3);
    ASSERT_EQ(filter.m_rests[0].get(0), 5);
    ASSERT_EQ(filter.m_rests[0].get(1), 7);
    ASSERT_EQ(filter.m_rests[0].get(2), 0);
    
    ASSERT_EQ(filter.m_offsets[0], 2);
}


TEST(RSQF_get, single_slot)
{
    for (uint64_t insert_slot{0} ; insert_slot<128 ; insert_slot++)
    {
        QuotientFilter<7, 5> filter{};
        const uint64_t insert_value{(insert_slot << 5) | 0b10001};
        filter.insert(insert_value);

        for(uint64_t query_slot{0} ; query_slot<128 ; query_slot++)
        {
            const uint64_t query_value{(query_slot << 5) | 0b10001};
            ASSERT_EQ(filter.get(query_value), query_slot == insert_slot);
        }
    }
}


TEST(RSQF_get, triple_insertion)
{
    for (uint64_t insert_slot{62} ; insert_slot<128 ; insert_slot++)
    {
        QuotientFilter<7, 5> filter{};
        for (uint64_t i{1} ; i<=3 ; i++)
        {
            const uint64_t insert_value{(insert_slot << 5) | 0b10000 | i};
            filter.insert(insert_value);
        }

        for(uint64_t query_slot{62} ; query_slot<65 ; query_slot++)
        {
            const uint64_t query_value{(query_slot << 5) | 0b10010};
            ASSERT_EQ(filter.get(query_value), query_slot == insert_slot);
        }
    }
}

TEST(RSQF_resize, resize_split_runs)
{
    QuotientFilter<7, 5> filter{};

    uint64_t values[] {
        // Run in block 0
        (2UL << 5) | 0b00101, (2UL << 5) | 0b10101,
        // Run in block 1
        (73UL << 5) | 0b00001, (73UL << 5) | 0b10100
    };

    for (uint64_t i{0} ; i<4 ; i++)
        filter.insert(values[i]);

    QuotientFilter<8, 4> resized{filter};

    ASSERT_EQ(filter.size(), resized.size());

    for (uint64_t i{0} ; i<4 ; i++)
    {
        ASSERT_EQ(resized.get(values[i]), true);
    }
}

TEST(RSQF_resize, resize_same_run)
{
    QuotientFilter<7, 5> filter{};

    uint64_t values[] {
        // Run in block 0
        (2UL << 5) | 0b00101, (2UL << 5) | 0b01101,
        (2UL << 5) | 0b00001, (2UL << 5) | 0b00100
    };

    for (uint64_t i{0} ; i<4 ; i++)
        filter.insert(values[i]);

    QuotientFilter<8, 4> resized{filter};

    ASSERT_EQ(filter.size(), resized.size());

    for (uint64_t i{0} ; i<4 ; i++)
    {
        ASSERT_EQ(resized.get(values[i]), true);
    }
}

TEST(RSQF_resize, resize_overflow)
{
    QuotientFilter<7, 5> filter{};

    uint64_t values[] {
        (127UL << 5) | 0b00101, // slot 127          then 254
        (127UL << 5) | 0b11101, // slot 127 overflow then 255
        (127UL << 5) | 0b10001  // slot 127 overflow then 255 overflow
    };

    for (uint64_t i{0} ; i<3 ; i++)
        filter.insert(values[i]);

    QuotientFilter<8, 4> resized{filter};

    ASSERT_EQ(filter.size(), resized.size());

    for (uint64_t i{0} ; i<3 ; i++)
    {
        ASSERT_EQ(resized.get(values[i]), true);
    }
}

