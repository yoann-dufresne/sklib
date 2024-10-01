#ifndef BITARRAY_H
#define BITARRAY_H

#include <cassert>
#include <bit>
#include <array>
#include <cstdint>
#include "bit.hpp"

namespace dbglib {

template <uint64_t size>
class BitArray
{	
	public:
		using underlying_type = uint64_t;
		using data_type = std::array<uint64_t, (size + 63UL) / 64>;

		BitArray () noexcept;
		data_type& modifiable_array_data() noexcept;
		const data_type& array_data() const noexcept;
		std::size_t first_one(const std::size_t from) const;
		void clear() noexcept;
		std::size_t get_size() const noexcept;
		void set(const uint64_t index);
		void unset(const uint64_t index);
		void set(const uint64_t index, bool val);
		bool get(const uint64_t index) const;
		std::size_t rank(const uint64_t start_idx, const uint64_t stop_idx) const;
		std::size_t select(const uint64_t start_idx, const uint64_t num_1s) const;
		void toric_right_shift(const uint64_t from, const uint64_t to);

	private:
		std::size_t select64(const uint64_t vector, const uint64_t num_1s) const;
		static constexpr std::size_t ut_bit_size = 8 * sizeof(underlying_type);
		static constexpr std::size_t block_size = size / ut_bit_size;
		data_type m_array; // Move this into private scope

		// having a static table containing all masks is actually slower than recomputing the mask on-the-fly since we incur in a cache miss
		// constexpr auto create_set_masks();
		// constexpr auto create_unset_masks();
		// static constexpr auto set_masks = create_set_masks();
		// static constexpr auto unset_masks = create_unset_masks();
};

template <uint64_t size>
BitArray<size>::BitArray() noexcept 
{
	clear();
}

template <uint64_t size>
typename BitArray<size>::data_type& 
BitArray<size>::modifiable_array_data() noexcept 
{
	return m_array;
}

template <uint64_t size>
const typename BitArray<size>::data_type& 
BitArray<size>::array_data() const noexcept 
{
	return m_array;
}

/** Get the position of the first bit set to 1, starting at position from. The search is toric.
 * @param from Bit position where to start the search.
 * @return The absolute position in the bitvector
 **/
template <uint64_t size>
std::size_t 
BitArray<size>::first_one(const std::size_t from) const 
{
	const uint64_t mask = (~static_cast<uint64_t>(0)) << (from % 64UL);
	std::size_t uint_idx = from / ut_bit_size;

	// First uint check
	uint64_t value {m_array[uint_idx] & mask};
	if (value == 0) {
		++uint_idx;
		if (uint_idx == block_size) {
			uint_idx = 0;
		}
		value = m_array[uint_idx];
	}

	// Following uints
	while (m_array[uint_idx] == 0) {
		++uint_idx;
		if (uint_idx == block_size) {
			uint_idx = 0;
		}
		value = m_array[uint_idx];
	}

	// Get the first one position
	return uint_idx * ut_bit_size + std::countr_zero(value);
}

template <uint64_t size>
void 
BitArray<size>::clear() noexcept 
{
	m_array.fill(0);
}

template <uint64_t size>
inline std::size_t 
BitArray<size>::get_size() const noexcept 
{
	return size;
}

template <uint64_t size>
void 
BitArray<size>::set(const uint64_t index) 
{
	m_array[index / ut_bit_size] |= static_cast<underlying_type>(1) << (index % ut_bit_size);
}

template <uint64_t size>
void 
BitArray<size>::unset(const uint64_t index) 
{
	m_array[index / ut_bit_size] &= ~(static_cast<underlying_type>(1) << (index % ut_bit_size));
}

template <uint64_t size>
void 
BitArray<size>::set(const uint64_t index, bool val) 
{
	if (val) set(index);
	else unset(index);
}

template <uint64_t size>
bool 
BitArray<size>::get(const uint64_t index) const 
{
	return (m_array.at(index / ut_bit_size) & (static_cast<underlying_type>(1) << (index % ut_bit_size))) != 0;
}

/** Shift the bitvector 1 bit to the right. from and to are included positions. Left bit is 0.
 * Right bit overflowing is discarded.
 * /!\ Do not work on toric queries that starts and ends in the same uint.
 * @param from First bit to shift.
 * @param to Last bit to shift (and to discard).
 **/
template <uint64_t size>
void 
BitArray<size>::toric_right_shift(const uint64_t from, const uint64_t to) 
{
	// std::cout << "shift from " << from << " to " << to << std::endl;
	const uint64_t first_uint {from / ut_bit_size};
	const uint64_t last_uint {to / ut_bit_size};
	uint64_t current_uint {first_uint};

	//  --- everything is in the same uint ---
	if ((first_uint == last_uint) and (from < to))
	{
		// - Prepare slice -
		// less significant part of the mask
		const uint64_t slice_right_mask {(~static_cast<uint64_t>(0UL)) >> ((ut_bit_size - 1) - (from % ut_bit_size))};
		// most significant part of the mask
		const uint64_t left_mask {(~static_cast<uint64_t>(0UL)) >> ((ut_bit_size - 1) - (to % ut_bit_size))};
		const uint64_t slice_full_mask {slice_right_mask ^ left_mask};
		// Get the slice to shift and shift it
		const uint64_t shifted_slice {(m_array[current_uint] << 1 ) & slice_full_mask};

		// - Prepare vector -
		// less significant part of the mask
		const uint64_t vector_right_mask {slice_right_mask >> 1};
		// most significant part of the mask
		const uint64_t vector_full_mask {~(left_mask ^ vector_right_mask)};
		// Remove the bits from the slice to shift
		m_array[current_uint] &= vector_full_mask;

		// - Merge everything -
		m_array[current_uint] |= shifted_slice;

		return;
	}

	// --- Multiple uint to modify ---
	uint64_t carry_bit {0};
	while (current_uint != last_uint) {
		// - Mask the beginning if this is the first uint to be shifted -
		if (current_uint == first_uint)
		{
			// Create mask to ignore less significant bits until from
			const uint64_t mask {(1UL << (from % ut_bit_size)) - 1UL};
			// extract and shift the slice
			const uint64_t shifted_slice {(m_array[current_uint] & ~mask) << 1};
			// Save the carry bit
			carry_bit = m_array[current_uint] >> (ut_bit_size - 1);
			// Remove the slice bits
			m_array[current_uint] &= mask;
			// Recompose the vector
			m_array[current_uint] |= shifted_slice;

			// increment
			if (current_uint + 1 == block_size) current_uint = 0;
			else ++current_uint;
			continue;
		}

		// - Shift a bitvector uint that is not in the extremities of the shifted area -
		assert(carry_bit <= 1);
		uint64_t new_value {carry_bit | m_array[current_uint] << 1};
		carry_bit = m_array[current_uint] >> (ut_bit_size - 1);
		m_array[current_uint] = new_value;

		// increment
		if (current_uint+1 == block_size) current_uint = 0;
		else ++current_uint;
	}

	// - shifts in the last uint -
	// Create mask to ignore most significant bits after to
	const uint64_t mask {(~static_cast<uint64_t>(0)) >> ((ut_bit_size - 1) - (to % ut_bit_size))};
	// extract and shift the slice
	const uint64_t shifted_slice {(m_array[current_uint] << 1) & mask};
	// Remove the slice bits
	m_array[current_uint] &= ~mask;
	// Recompose the vector
	m_array[current_uint] |= shifted_slice | carry_bit;
}

/** Count the number of bits set to 1 from start_idx to rank_idx. This rank operation is toric.
 * /!\ Known limitation: You cannot rank the whole bitvector in a toric way. If you are asking results on toric
 * query, the ranked bit must be in a different uint than the start bit.
 * @param start_idx The bit index from which the rank should start (Bit 0 in classic rank)
 * @param stop_idx Position in the bitvector where to stop the ranking (position included)
 * 
 * @return The number of 1s found between the 2 boundaries (included)
 **/
template <uint64_t size>
std::size_t
BitArray<size>::rank(const uint64_t start_idx, const uint64_t stop_idx) const 
{
	assert(start_idx < size);
	assert(stop_idx < size);

	const std::size_t start_block = start_idx / ut_bit_size;
	const std::size_t stop_block = stop_idx / ut_bit_size;
	const std::size_t start_local_idx = start_idx % ut_bit_size;
	const std::size_t stop_local_idx = stop_idx % ut_bit_size;

	std::size_t rank = bit::popcount(m_array[start_block] >> start_local_idx);
	if (start_idx == stop_idx) {
		rank -= bit::popcount(m_array[stop_block] >> stop_local_idx);
	} else if (start_idx < stop_idx) {
		for (std::size_t block_idx = start_block + 1; block_idx != stop_block; ++block_idx) {
			rank += bit::popcount(m_array[block_idx]);
		}
		rank += bit::popcount(m_array[stop_block] & ((static_cast<underlying_type>(1) << stop_local_idx) - 1));
	} else { // circular
		for (std::size_t block_idx = start_block + 1; block_idx != )
		for
	}
}
// {
// 	const std::size_t start_uint = start_idx / ut_bit_size;
// 	const std::size_t stop_uint = stop_idx / ut_bit_size;
// 	std::size_t current_uint = start_uint;

// 	// --- Both limits are in the same uint ---
// 	if (start_uint == stop_uint) {
// 		// Create a mask
// 		const uint64_t start_mask {(~static_cast<uint64_t>(0)) << (start_idx % ut_bit_size)};
// 		const uint64_t stop_mask {(~static_cast<uint64_t>(0)) >> ((ut_bit_size - 1) - (stop_idx % ut_bit_size))};
// 		const uint64_t mask {start_mask & stop_mask};

// 		// Return the number of 1s
// 		return std::popcount<uint64_t>(m_array[start_uint] & mask);
// 	}

// 	// --- Multiple uint coverage ---
// 	// Start uint
// 	const uint64_t start_mask {(~static_cast<uint64_t>(0)) << (start_idx % ut_bit_size)};
// 	uint64_t num_ones {static_cast<uint64_t>(std::popcount(m_array[start_uint] & start_mask))};

// 	if (current_uint + 1 == block_size) current_uint = 0;
// 	else ++current_uint;

// 	// Middle uints
// 	while (current_uint != stop_uint) {
// 		num_ones += std::popcount(m_array[current_uint]);

// 		// uint increment
// 		if (current_uint + 1 == block_size) current_uint = 0;
// 		else ++current_uint;
// 	}

// 	// Last uint
// 	const uint64_t stop_mask {(~static_cast<uint64_t>(0)) >> ((ut_bit_size - 1) - (stop_idx % ut_bit_size))};
// 	num_ones += std::popcount(m_array[stop_uint] & stop_mask);

// 	return num_ones;
// }

/** Return the position of the val-th bit set starting from index start_idx. The answer is computed using toricity
 * if needed.
 * /!\ Limitation: val == 0 leads to undefined behaviour.
 * @param start_idx The bit index from which the select should start (Bit 0 in classic select)
 * @param val The number of bits to count (Must be > 0)
 * 
 * @return The absolute position of the val-th set bit.
 **/
template <uint64_t size>
std::size_t 
BitArray<size>::select(const uint64_t start_idx, const uint64_t num_1s) const 
{
	const std::size_t start_uint {start_idx / ut_bit_size};
	std::size_t current_uint {start_uint};
	std::size_t remaining_1s {num_1s};

	// --- Multiple uint coverage ---
	// Start uint
	const uint64_t start_mask {(~static_cast<uint64_t>(0)) << (start_idx % ut_bit_size)};
	const uint64_t buffered_uint {m_array[start_uint] & start_mask};
	const uint64_t first_1s {static_cast<uint64_t>(std::popcount(buffered_uint))};

	// Update select
	if (first_1s >= remaining_1s) {
		return current_uint * ut_bit_size + select64(buffered_uint, num_1s);
	} else {
		remaining_1s -= first_1s;
	}

	// Update block
	if (current_uint + 1 == num_uint) current_uint = 0;
	else ++current_uint;

	// Middle uints
	while (current_uint != start_uint) {
		const uint64_t local_1s {static_cast<uint64_t>(std::popcount(m_array[current_uint]))};

		// Update select
		if (local_1s >= remaining_1s) {
			return current_uint * ut_bit_size + select64(m_array[current_uint], remaining_1s);
		} else {
			remaining_1s -= local_1s;
		}

		// uint increment
		if (current_uint + 1 == num_uint) current_uint = 0;
		else ++current_uint;
	}
	return size;
}

template <uint64_t size>
std::size_t 
BitArray<size>::select64(const uint64_t vector, const uint64_t num_1s) const
{
	uint64_t pdep_val = 
#if defined(FOR_X86_64_PLATFORM)
	_pdep_u64(1ULL << (num_1s - 1), vector);
#elif defined(FOR_ARM_PLATFORM)
	bit::pdep_u64(1ULL << (num_1s - 1), vector);
#else
	static_assert(0, "Something went wrong with the preprocessor macro declarations in CMake");
#endif
	auto ret = 
#if __cplusplus >= 202002L
	std::countr_zero<uint64_t>(pdep_val);
#else
	bit::countr_zero<uint64>(pdep_val);
#endif
	return ret;
}

// template <uint64_t size>
// constexpr auto 
// BitArray<size>::create_set_masks() {
// 	std::array<uint64_t, 64> array {};
// 	for (std::size_t i = 0; i < 64; ++i) {
// 		array[i] = 1UL << i;
// 	}
// 	return array;
// }

// template <uint64_t size>
// constexpr auto 
// BitArray<size>::create_unset_masks() {
// 	std::array<uint64_t, 64> array {};
// 	for (std::size_t i = 0; i < 64; ++i) {
// 		array[i] = ~set_masks[i];
// 	}
// 	return array;
// }

} // namespace dbglib

#endif // BITARRAY_HPP