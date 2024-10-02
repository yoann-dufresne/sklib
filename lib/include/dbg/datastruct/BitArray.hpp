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
		void set(const std::size_t index);
		void unset(const std::size_t index);
		void set(const std::size_t index, bool val);
		bool get(const std::size_t index) const;
		std::size_t rank(const std::size_t start_idx, const std::size_t stop_idx) const;
		std::size_t select(const std::size_t start_idx, const std::size_t rank) const;
		void toric_right_shift(const std::size_t from, const std::size_t to);

	private:
		underlying_type shift_slice(underlying_type v, std::size_t start, std::size_t stop);
		static constexpr std::size_t ut_bit_size = 8 * sizeof(underlying_type);
		static constexpr std::size_t block_size = size / ut_bit_size + (size % ut_bit_size == 0 ? 0 : 1);
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
BitArray<size>::first_one(const std::size_t start_idx) const 
{
	return select(start_idx, 1);
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
BitArray<size>::set(const std::size_t index) 
{
	m_array[index / ut_bit_size] |= static_cast<underlying_type>(1) << (index % ut_bit_size);
}

template <uint64_t size>
void 
BitArray<size>::unset(const std::size_t index) 
{
	m_array[index / ut_bit_size] &= ~(static_cast<underlying_type>(1) << (index % ut_bit_size));
}

template <uint64_t size>
void 
BitArray<size>::set(const std::size_t index, bool val) 
{
	if (val) set(index);
	else unset(index);
}

template <uint64_t size>
bool 
BitArray<size>::get(const std::size_t index) const 
{
	return (m_array.at(index / ut_bit_size) & (static_cast<underlying_type>(1) << (index % ut_bit_size))) != 0;
}

/** Shift the bitvector 1 bit to the right. 
 * start_idx and stop_idx are included positions. 
 * Inserted bit to the left is always 0.
 * Bit overflowing to the right is discarded.
 * /!\ Do not work on toric queries that starts and ends in the same uint.
 * @param start_idx First bit to shift.
 * @param stop_idx Last bit to shift (and to discard).
 **/
template <uint64_t size>
void 
BitArray<size>::toric_right_shift(const std::size_t start_idx, const std::size_t stop_idx) 
{
	assert(start_idx < size);
	assert(stop_idx < size);

	const std::size_t start_block = start_idx / ut_bit_size;
	const std::size_t stop_block = stop_idx / ut_bit_size;
	const std::size_t start_local_idx = start_idx % ut_bit_size;
	const std::size_t stop_local_idx = stop_idx % ut_bit_size;

	if (start_idx < stop_idx and start_block == stop_block) {
		m_array[start_block] = shift_slice(m_array[start_block], start_local_idx, stop_local_idx);
	} else {
		const underlying_type carry_mask = (static_cast<underlying_type>(1) << (ut_bit_size - 1));
		underlying_type carry;
		{ // shift first part
			carry = m_array[start_block] & carry_mask ? 1 : 0;
			m_array[start_block] = shift_slice(m_array[start_block], start_local_idx, ut_bit_size - 1);
		}
		{ // shift blocks between start_block and end_block
			auto block_idx = start_block + 1;
			if (block_idx >= block_size) block_idx = 0;
			while (block_idx != stop_block) {
				auto tmp_carry = m_array[block_idx] & carry_mask ? 1 : 0;
				m_array[block_idx] = (m_array[block_idx] << 1) | carry;
				carry = tmp_carry;
				++block_idx;
				if (block_idx >= block_size) block_idx = 0;
			}
		}
		{ // shift last part
			m_array[stop_block] = shift_slice(m_array[stop_block], 0, stop_local_idx);
			m_array[stop_block] |= carry;
		}
	}
}

/**
 * shift a slice of an integer to the left.
 * 0 is inserted to the right while the left-most bit of the slice is discarded.
 * Both start and stop indexes start from the "right" i.e.: [7 6 5 4 3 2 1 0]
 * Both indexes are INCLUSIVE.
 * @param start bit index from where to start
 * @param stop bit index of the discarded bit
 */
template <uint64_t size>
typename BitArray<size>::underlying_type 
BitArray<size>::shift_slice(underlying_type v, std::size_t start, std::size_t stop)
{
	assert(start < ut_bit_size);
	assert(stop < ut_bit_size);
	assert(start <= stop);

	underlying_type mask = ((static_cast<underlying_type>(1) << start) - 1);
	++stop;
	mask |= stop == ut_bit_size ? 0 : ~((static_cast<underlying_type>(1) << stop) - 1);
	auto slice = v & ~mask; // take bit slice
	v &= mask; // clear bit-slice from source
	slice = (slice << 1) & ~mask; // shift slice by one, clear right-most bit by re-applying the (inverted) mask
	v |= slice;
	return v;
}

// {
// 	// std::cout << "shift from " << from << " to " << to << std::endl;
// 	const uint64_t first_uint {from / ut_bit_size};
// 	const uint64_t last_uint {to / ut_bit_size};
// 	uint64_t current_uint {first_uint};

// 	//  --- everything is in the same uint ---
// 	if ((first_uint == last_uint) and (from < to)) {
// 		// - Prepare slice -
// 		// less significant part of the mask
// 		const uint64_t slice_right_mask {(~static_cast<uint64_t>(0UL)) >> ((ut_bit_size - 1) - (from % ut_bit_size))};
// 		// most significant part of the mask
// 		const uint64_t left_mask {(~static_cast<uint64_t>(0UL)) >> ((ut_bit_size - 1) - (to % ut_bit_size))};
// 		const uint64_t slice_full_mask {slice_right_mask ^ left_mask};
// 		// Get the slice to shift and shift it
// 		const uint64_t shifted_slice {(m_array[current_uint] << 1 ) & slice_full_mask};

// 		// - Prepare vector -
// 		// less significant part of the mask
// 		const uint64_t vector_right_mask {slice_right_mask >> 1};
// 		// most significant part of the mask
// 		const uint64_t vector_full_mask {~(left_mask ^ vector_right_mask)};
// 		// Remove the bits from the slice to shift
// 		m_array[current_uint] &= vector_full_mask;

// 		// - Merge everything -
// 		m_array[current_uint] |= shifted_slice;

// 		return;
// 	}

// 	// --- Multiple uint to modify ---
// 	uint64_t carry_bit  = 0;
// 	while (current_uint != last_uint) {
// 		// - Mask the beginning if this is the first uint to be shifted -
// 		if (current_uint == first_uint) {
// 			// Create mask to ignore less significant bits until from
// 			const uint64_t mask {(1UL << (from % ut_bit_size)) - 1UL};
// 			// extract and shift the slice
// 			const uint64_t shifted_slice {(m_array[current_uint] & ~mask) << 1};
// 			// Save the carry bit
// 			carry_bit = m_array[current_uint] >> (ut_bit_size - 1);
// 			// Remove the slice bits
// 			m_array[current_uint] &= mask;
// 			// Recompose the vector
// 			m_array[current_uint] |= shifted_slice;

// 			// increment
// 			if (current_uint + 1 == block_size) current_uint = 0;
// 			else ++current_uint;
// 			continue;
// 		}

// 		// - Shift a bitvector uint that is not in the extremities of the shifted area -
// 		assert(carry_bit <= 1);
// 		uint64_t new_value {carry_bit | m_array[current_uint] << 1};
// 		carry_bit = m_array[current_uint] >> (ut_bit_size - 1);
// 		m_array[current_uint] = new_value;

// 		// increment
// 		if (current_uint+1 == block_size) current_uint = 0;
// 		else ++current_uint;
// 	}

// 	// - shifts in the last uint -
// 	// Create mask to ignore most significant bits after to
// 	const uint64_t mask {(~static_cast<uint64_t>(0)) >> ((ut_bit_size - 1) - (to % ut_bit_size))};
// 	// extract and shift the slice
// 	const uint64_t shifted_slice {(m_array[current_uint] << 1) & mask};
// 	// Remove the slice bits
// 	m_array[current_uint] &= ~mask;
// 	// Recompose the vector
// 	m_array[current_uint] |= shifted_slice | carry_bit;
// }

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
BitArray<size>::rank(const std::size_t start_idx, const std::size_t stop_idx) const 
{
	assert(start_idx < size);
	assert(stop_idx < size);

	const std::size_t start_block = start_idx / ut_bit_size;
	const std::size_t stop_block = stop_idx / ut_bit_size;
	const std::size_t start_local_idx = start_idx % ut_bit_size;
	const std::size_t stop_local_idx = stop_idx % ut_bit_size;

	std::size_t rank = bit::popcount(m_array.at(start_block) >> start_local_idx);
	if (start_idx <= stop_idx) {
		if (start_block == stop_block) {
			auto shift = stop_local_idx + 1;
			rank -= shift == ut_bit_size ? 0 : bit::popcount(m_array.at(stop_block) >> shift);
		} else {
			for (std::size_t block_idx = start_block + 1; block_idx != stop_block; ++block_idx) {
				rank += bit::popcount(m_array.at(block_idx));
			}
			auto shift = stop_local_idx + 1;
			auto shifted = shift == ut_bit_size ? 0 : (static_cast<underlying_type>(1) << shift);
			rank += bit::popcount(m_array.at(stop_block) & (shifted - 1));
		}
	} else { // circular
		for (std::size_t block_idx = start_block + 1; block_idx != block_size; ++block_idx) rank += bit::popcount(m_array.at(block_idx));
		for (std::size_t block_idx = 0; block_idx != stop_block; ++block_idx) rank += bit::popcount(m_array.at(block_idx));
		auto shift = stop_local_idx + 1;
		auto shifted = shift == ut_bit_size ? 0 : (static_cast<underlying_type>(1) << shift);
		rank += bit::popcount(m_array.at(stop_block) & (shifted - 1));
	}
	return rank;
}

/** Return the position of the val-th bit set starting from index start_idx. The answer is computed using toricity
 * if needed.
 * /!\ Limitation: val == 0 leads to undefined behaviour.
 * @param start_idx The bit index from which the select should start (Bit 0 in classic select)
 * @param rank The number of bits to count (Must be > 0)
 * 
 * @return The absolute position of the rank-th set bit starting from start_idx. size if number of set bits < rank.
 **/
template <uint64_t size>
std::size_t 
BitArray<size>::select(const std::size_t start_idx, const std::size_t rank) const
{
	assert(start_idx < size);
	assert(rank > 0);

	const std::size_t start_block = start_idx / ut_bit_size;
	const std::size_t start_local_idx = start_idx % ut_bit_size;
	auto remaining = rank;
	
	auto block_itr = start_block;
	auto shift = start_local_idx;
	auto pcnt = bit::popcount(m_array[block_itr] >> shift);
	if (pcnt < remaining) {
		// the first iteration is unrolled for the bound check in the while loop
		remaining -= pcnt;
		++block_itr;
		if (block_itr == block_size) block_itr = 0;
		pcnt = bit::popcount(m_array[block_itr]);
		shift = 0;
		while(block_itr != start_block and pcnt < remaining) {
			remaining -= pcnt;
			++block_itr;
			if (block_itr == block_size) block_itr = 0;
			pcnt = bit::popcount(m_array[block_itr]);
		}
	}
	// we now know that the block we are point at has the bit we want
	return bit::select1(m_array[block_itr] >> shift, remaining - 1) + shift + block_itr * ut_bit_size;
}

} // namespace dbglib

#endif // BITARRAY_HPP