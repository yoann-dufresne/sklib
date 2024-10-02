#include <gtest/gtest.h>
#include <dbg/datastruct/BitArray.hpp>
#include <random>

using namespace dbglib;

TEST(BitArray, init)
{
    BitArray<128> vector;
    ASSERT_EQ(vector.array_data().size(), 2);

    BitArray<129> vector_3uints;
    ASSERT_EQ(vector_3uints.array_data().size(), 3);
}

TEST(bitarray, one_bit_set)
{
    BitArray<128> bv;
    
    for (size_t i = 0; i < 128; ++i)
    {
        bv.set(i);
        if (i < 64) {
            ASSERT_EQ(bv.array_data()[0], 1ULL << i);
            ASSERT_EQ(bv.array_data()[1], 0);
        } else {
            ASSERT_EQ(bv.array_data()[0], 0);
            ASSERT_EQ(bv.array_data()[1], 1ULL << (i - 64));
        }
        bv.unset(i);
        ASSERT_EQ(bv.array_data()[0], 0);
        ASSERT_EQ(bv.array_data()[1], 0);
    }
}

TEST(bitarray, unset_set_full)
{
    BitArray<128> bv;
    uint64_t expected = 0;
    uint64_t max = ~expected;
    
    // --- Set the bitarray bit by bit and test it ---

    for (size_t i = 0; i < 128; ++i) {
        if (i == 64) expected = 0;
        expected <<= 1;
        expected |= 1;

        bv.set(i);
        if (i < 64) {
            ASSERT_EQ(bv.array_data()[0], expected);
            ASSERT_EQ(bv.array_data()[1], 0);
        } else {
            ASSERT_EQ(bv.array_data()[0], max);
            ASSERT_EQ(bv.array_data()[1], expected);
        }
    }

    for (size_t i = 0; i < 128; ++i) {
        if (i % 64 == 0) expected = max;
        expected <<= 1;

        bv.unset(i);
        if (i < 64) {
            ASSERT_EQ(bv.array_data()[0], expected);
            ASSERT_EQ(bv.array_data()[1], max);
        } else {
            ASSERT_EQ(bv.array_data()[0], 0);
            ASSERT_EQ(bv.array_data()[1], expected);
        }
    }
}

TEST(bitarray, get_1bit)
{
    BitArray<128> bv;
    uint64_t set_value;

    // --- Set the bitarray bit by bit and test it ---

    for (size_t i = 0; i < 128; ++i) {
        // Compute the value to simulate the bitarray
        if (i % 64 == 0) set_value = 1;
        else set_value <<= 1;

        // Set the bitarray value by replacing 64 bits chunks
        if (i < 64) {
            bv.modifiable_array_data()[0] = set_value;
            bv.modifiable_array_data()[1] = 0;
        } else {
            bv.modifiable_array_data()[0] = 0;
            bv.modifiable_array_data()[1] = set_value;
        }

        // Test all the 128 values
        for (uint64_t j = 0; j < 128; ++j) {
            ASSERT_EQ(bv.get(j), i==j);
        }
    }
}

TEST(bitarray, randomized_set_and_get)
{
    #define BV_SIZE 10000
    BitArray<BV_SIZE> bv;
    std::mt19937 gen(0);
    std::uniform_int_distribution<std::size_t> distrib(0, BV_SIZE - 1);
    std::vector<std::size_t> indexes;
    for (std::size_t i = 0; i < BV_SIZE / 2; ++i) {
        auto index = distrib(gen);
        indexes.push_back(index);
        bv.set(index);
    }
    std::sort(indexes.begin(), indexes.end());
    indexes.erase(unique(indexes.begin(), indexes.end() ), indexes.end());
    std::vector<std::size_t> tobechecked;
    for (std::size_t i = 0; i < bv.get_size(); ++i) {
        if (bv.get(i)) tobechecked.push_back(i);
    }
    ASSERT_EQ(indexes, tobechecked);
#undef BV_SIZE
}

TEST(bitarray_shift, single_uint)
{
    // 3 uints bitarray
    BitArray<192> bv;
    // Init
    for (uint64_t i{3} ; i<=6 ; i++) bv.set(i);

    // Large shift : Everythng should be shifted of 1 position
    bv.toric_right_shift(0, 10);
    ASSERT_EQ(bv.array_data()[0], 0b11110000UL);
    ASSERT_EQ(bv.array_data()[1], 0);
    ASSERT_EQ(bv.array_data()[2], 0);

    // Reinit
    bv.clear();
    for (uint64_t i{3} ; i<=6 ; i++) bv.set(i);

    // shift 1s : 1 bit should be lost by overflowing
    bv.toric_right_shift(3, 6);
    ASSERT_EQ(bv.array_data()[0], 0b01110000UL);
    ASSERT_EQ(bv.array_data()[1], 0);
    ASSERT_EQ(bv.array_data()[2], 0);

    // Reinit
    bv.clear();
    for (uint64_t i{3} ; i<=6 ; i++) bv.set(i);

    // Truncated shift : The less significant bit of the slice should be set to 0 after shift
    bv.toric_right_shift(5, 10);
    ASSERT_EQ(bv.array_data()[0], 0b11011000UL);
    ASSERT_EQ(bv.array_data()[1], 0);
    ASSERT_EQ(bv.array_data()[2], 0);
}


TEST(bitarray_shift, three_uint)
{
    // 3 uints bitarray
    BitArray<192> bv;
    // Init
    for (std::size_t i = 62; i <= 65; ++i) bv.set(i);
    for (std::size_t i = 126; i <= 129; ++i) bv.set(i);

    // Initial_verification
    const uint64_t expected_uints[] {(0b11UL << 62), (0b11UL << 62) | 0b11UL , 0b11UL};
    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.array_data()[i], expected_uints[i]);
    }

    // Large shift : Everythng should be shifted of 1 position
    bv.toric_right_shift(0, 190);
    const uint64_t expected_large_shift[] {(0b1UL << 63), (0b1UL << 63) | 0b111UL , 0b111UL};
    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.array_data()[i], expected_large_shift[i]);
    }

    // Reinit
    bv.clear();
    for (std::size_t i = 62; i <= 65; ++i) bv.set(i);
    for (std::size_t i = 126; i <= 129; ++i) bv.set(i);

    // Large shift : left half of the vector should be 1 position shifted
    bv.toric_right_shift(0, 100);
    const uint64_t expected_half_shift[] {(0b1UL << 63), (0b11UL << 62) | 0b111UL , 0b11UL};
    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.array_data()[i], expected_half_shift[i]);
    }

    // Reinit
    bv.clear();
    for (std::size_t i = 62; i <= 65; ++i) bv.set(i);
    for (std::size_t i = 126; i <= 129; ++i) bv.set(i);

    // Truncated shift : shift on the uint border
    bv.toric_right_shift(63, 100);
    const uint64_t expected_trunckated_shift[] {(0b1UL << 62), (0b11UL << 62) | 0b111UL , 0b11UL};
    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.array_data()[i], expected_trunckated_shift[i]);
    }
}


TEST(bitarray_shift, toric)
{
    // 2 uints bitarray
    BitArray<128> bv;
    // Init
    for (uint64_t i = 126; i <= 129; ++i) bv.set(i % 128);

    // Init checks
    ASSERT_EQ(bv.array_data()[0], 0b11UL);
    ASSERT_EQ(bv.array_data()[1], 0b11UL << 62);

    // Large shift : Everythng should be shifted of 1 position
    bv.toric_right_shift(120, 10);
    ASSERT_EQ(bv.array_data()[0], 0b111UL);
    ASSERT_EQ(bv.array_data()[1], 0b1UL << 63);

    // Reinit
    bv.clear();
    for (uint64_t i{126} ; i<=129 ; i++) bv.set(i%128);

    // Truncated shift : shift on the uint border
    // bv.toric_right_shift(127, 10);
    // ASSERT_EQ(bv.array_data()[0], 0b111UL);
    // ASSERT_EQ(bv.array_data()[1], 0b1UL << 62);
}


TEST(bitarray_rank, single_uint)
{
    // 3 uints bitarray
    BitArray<64> bv;
    // Init
    for (std::size_t i = 0; i < 64; i += 8) { 
        bv.set(i); 
        bv.set(i+7); 
    }

    ASSERT_EQ(bv.rank(0, 63), 16);
    ASSERT_EQ(bv.rank(7, 7), 1);
    ASSERT_EQ(bv.rank(31, 32), 2);
}

TEST(bitarray_rank, triple_uint)
{
    // 3 uints bitarray
    BitArray<192> bv;
    // Init
    for (uint64_t i{0} ; i<192 ; i+=8) { bv.set(i); bv.set(i+7); }

    // 2 uints ranks
    ASSERT_EQ(bv.rank(0, 127), 32);
    ASSERT_EQ(bv.rank(32, 95), 16);
    ASSERT_EQ(bv.rank(0, 64), 17);
    ASSERT_EQ(bv.rank(63, 127), 17);
    // 3 uints ranks
    ASSERT_EQ(bv.rank(0, 191), 48);
    ASSERT_EQ(bv.rank(32, 159), 32);
    ASSERT_EQ(bv.rank(63, 191), 33);
    ASSERT_EQ(bv.rank(0, 128), 33);
}

TEST(bitarray_rank, toric)
{
    // 3 uints bitarray
    BitArray<256> bv;
    // // Init
    for (uint64_t i{0} ; i<256 ; i+=8) { bv.set(i); bv.set(i+7); }

    ASSERT_EQ(bv.rank(128, 127), 64);
    ASSERT_EQ(bv.rank(64, 63), 64);
    ASSERT_EQ(bv.rank(255, 0), 2);
}


TEST(bitarray_rank, randomized)
{
#define BV_SIZE 3000
    BitArray<BV_SIZE> bv;
    std::mt19937 gen(0);
    std::uniform_int_distribution<std::size_t> distrib(0, BV_SIZE - 1);
    for (std::size_t i = 0; i < BV_SIZE / 2; ++i) {
        auto index = distrib(gen);
        bv.set(index);
    }
    std::size_t cur_rank = 0;
    bool ok = true;
    for (std::size_t i = 0; i < bv.get_size(); ++i) {
        if (bv.get(i)) ++cur_rank;
        if (ok and bv.rank(0, i) != cur_rank) {
            // std::cerr << "bv.rank(" << i << ") = " << bv.rank(0, i) << ", check = " << cur_rank << "\n";
            ok = false;
        }
    }
    ASSERT_TRUE(ok);
#undef BV_SIZE
}


TEST(bitarray_select, single_uint)
{
    // 3 uints bitarray
    BitArray<64> bv;
    // Init
    for (size_t i = 0; i < 64; i += 8) { 
        bv.set(i);
        bv.set(i+7);
    }

    for (std::size_t i = 0; i < 64 / 8; ++i) {
        ASSERT_EQ(bv.select(0, 2*i+1), 8*i  );
        ASSERT_EQ(bv.select(0, 2*i+2), 8*i+7);
    }
}

TEST(bitarray_select, toric)
{
    // 2 uints bitarray
    BitArray<128> bv;
    // Init
    for (std::size_t i = 0; i < 128; i += 8) { 
        bv.set(i);
        bv.set(i+7);
    }

    for (uint64_t i = 0; i < (128 / 8); ++i)
    {
        ASSERT_EQ(bv.select(64, 2*i+1), (64 + 8 * i    ) % 128 );
        ASSERT_EQ(bv.select(64, 2*i+2), (64 + 8 * i + 7) % 128 );
    }
}


TEST(bitarray_select, multiple_uints)
{
    // 4 uints bitarray
    BitArray<256> bv;
    // // Init
    for (std::size_t i = 0; i < 256; i += 8) { 
        bv.set(i);
        bv.set(i+7);
    }

    ASSERT_EQ(bv.select(0, 1), 0);
    ASSERT_EQ(bv.select(0, 64), 255);
    ASSERT_EQ(bv.select(128, 33), 0);
    ASSERT_EQ(bv.select(128, 34), 7);
    ASSERT_EQ(bv.select(128, 49), 64);
}


TEST(bitarray_select, randomized)
{
#define BV_SIZE 3000
    BitArray<BV_SIZE> bv;
    std::mt19937 gen(0);
    std::uniform_int_distribution<std::size_t> distrib(0, BV_SIZE - 1);
    std::vector<std::size_t> indexes;
    for (std::size_t i = 0; i < BV_SIZE / 2; ++i) {
        auto index = distrib(gen);
        indexes.push_back(index);
        bv.set(index);
    }
    std::sort(indexes.begin(), indexes.end());
    indexes.erase(unique(indexes.begin(), indexes.end() ), indexes.end());
    // std::cerr << "indexes.size = " << indexes.size() << "\n";
    for (std::size_t i = 0; i < indexes.size(); ++i) {
        // std::cerr << "i = " << i << ", bv.select(" << i << ") = " << bv.select(0, i+1) << ", indexes[" << i << "] = " << indexes.at(i) << std::endl;
        ASSERT_EQ(bv.select(0, i+1), indexes.at(i));
    }
#undef BV_SIZE
}


TEST(bitarray_first_one, single_uint)
{
    BitArray<64> bv;

    // Set all the 2^i bits
    for (std::size_t i = 0; i < 6; ++i)
        bv.set(1UL<<i);

    ASSERT_EQ(bv.first_one(0), 1);
    ASSERT_EQ(bv.first_one(2), 2);
    ASSERT_EQ(bv.first_one(3), 4);
    ASSERT_EQ(bv.first_one(5), 8);
    ASSERT_EQ(bv.first_one(9), 16);
    ASSERT_EQ(bv.first_one(17), 32);
    ASSERT_EQ(bv.first_one(33), 1);
    ASSERT_EQ(bv.first_one(63), 1);
}

TEST(bitarray_first_one, multiple_uints)
{
    BitArray<1024> bv;

    // Set all the 2^i bits
    for (std::size_t i = 0; i < 10; ++i)
        bv.set(1UL<<i);

    ASSERT_EQ(bv.first_one(0), 1);
    ASSERT_EQ(bv.first_one(2), 2);
    ASSERT_EQ(bv.first_one(3), 4);
    ASSERT_EQ(bv.first_one(5), 8);
    ASSERT_EQ(bv.first_one(9), 16);
    ASSERT_EQ(bv.first_one(17), 32);
    ASSERT_EQ(bv.first_one(33), 64);
    ASSERT_EQ(bv.first_one(65), 128);
    ASSERT_EQ(bv.first_one(129), 256);
    ASSERT_EQ(bv.first_one(257), 512);
    ASSERT_EQ(bv.first_one(513), 1);
    ASSERT_EQ(bv.first_one(1023), 1);
}