#ifndef PACKED_ARRAY_H
#define PACKED_ARRAY_H

#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>

/** Contains a packed array of 64 values of size slot_size. This is a building block for packed hash tables, counting bloom filters and quotient filter.
 **/
template<uint64_t slot_size>
struct PackedBlock
{
	std::array<uint64_t, slot_size> m_array{};
	
	PackedBlock() {};

	/** Get the value at index (converted in uint64_t).
	 **/
	uint64_t get(const size_t index) const
	{
		// cout << "get " << index << endl;
		// Get the position of the value in bits
		const auto first_bit_position { index * slot_size };
		const auto first_uint_position { first_bit_position / 64 };
		const auto bits_in_uint { std::min((64UL - (first_bit_position % 64)), slot_size) };
		const bool second_byte { bits_in_uint < slot_size };

		// cout << first_bit_position << " " << first_uint_position << endl;
		// cout << first_bit_position << "/" << slot_size << " = " << first_uint_position << endl;
		// cout << bits_in_uint << " " << second_byte << endl;

		// Get the first uint
		auto value { m_array[first_uint_position] };

		// Shift right if only one uint
		if (not second_byte)
		{
			const auto right_shift { 64UL - slot_size - (first_bit_position % 64) };
			value >>= right_shift;
		}
		// Shift left and add the second uint if needed
		else
		{
			const auto second_bits { slot_size - bits_in_uint };
			// Shift left to allow insertion
			value <<= second_bits;
			// Mask to remove right bits
			const auto left_mask { (1UL << second_bits) - 1 };
			const auto right_mask { ~left_mask };
			value &= right_mask;

			// Get the second uint
			auto right_value { m_array[first_uint_position+1] };
			// Shift right to align
			right_value >>= 64UL - second_bits;
			// Mask to remove the left bits
			right_value &= left_mask;

			// Compose the final value from the two parts
			value |= right_value;
		}

		// Mask to remove left bits
		constexpr uint64_t mask { (1UL << slot_size) - 1UL };
		value &= mask;

		// cout << "/get" << endl;
		return value;
	};


	/** Set the slot_size lower bits of the value to the index part of the array
	 * @param index Where to insert
	 * @param value Value to insert (slot_size lowest bits only)
	 **/
	void set(const size_t index, uint64_t value)
	{
		// Get the position of the value in bits
		const auto first_bit_position { index * slot_size };
		const auto first_uint_position { first_bit_position / 64 };
		const auto bits_in_uint { std::min((64UL - (first_bit_position % 64)), slot_size) };
		const bool second_byte { bits_in_uint < slot_size };

		if (not second_byte)
		{
			const auto empty_trailing_bits {64 - (first_bit_position % 64) - slot_size};
			// cout << "empty_trailing_bits " << 
			// Align the value
			value <<= empty_trailing_bits;
			// Mask already present bits
			uint64_t mask = first_bit_position == 0 ? ~0UL : (1UL << (empty_trailing_bits + slot_size)) - 1;
			mask ^= (1UL << empty_trailing_bits) - 1;
			mask = ~mask;
			m_array[first_uint_position] &= mask;
			// Insert the value
			m_array[first_uint_position] |= value;
		}
		else
		{
			// --- first uint ---
			// Right align the higher bits of the value
			const auto right_shift {slot_size - bits_in_uint};
			const uint64_t right_value {value >> right_shift};
			// Mask the right bits of the first uint
			const uint64_t right_mask {~((1UL << bits_in_uint) - 1)};
			m_array[first_uint_position] &= right_mask;
			// Insert the high bits
			m_array[first_uint_position] |= right_value;

			// --- second uint ---
			const auto bits_in_second_uint {slot_size - bits_in_uint};
			// Left align the lower bits of the value
			const auto left_shift {64 - bits_in_second_uint};
			const uint64_t left_value {value << left_shift};
			// Mask the left bits of the second uint
			const uint64_t left_mask { (1UL << (64 - bits_in_second_uint)) - 1 };
			m_array[first_uint_position+1] &= left_mask;
			// Insert the low bits
			m_array[first_uint_position+1] |= left_value;
		}
	};
};

// template<typename T>
// std::ostream& operator<<(std::ostream& os, PackedBlock<T>& kmer)
// {
// 	os << kmer.m_value;
//     return os;
// };

#endif
