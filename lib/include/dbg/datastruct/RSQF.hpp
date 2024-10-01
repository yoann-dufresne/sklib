#ifndef RSQF_H
#define RSQF_H

#include <array>
#include <cstdint>
#include <bit>

#include <dbg/datastruct/PackedArray.hpp>
#include <dbg/datastruct/BitArray.hpp>
#include <dbg/datastruct/quotienting.hpp>

namespace dbglib {

template<uint64_t q, uint64_t r, class Quotienting = LeftQuotienting>
class QuotientFilter
{
    public: // TODO switch to private
        // Number of blocks. Defined at compile time
        static constexpr uint64_t num_blocks {1UL << (q - 6)};
        static constexpr uint64_t num_slots {num_blocks * 64UL};
        // Storage area of rests
        std::array<PackedBlock<r>, QuotientFilter::num_blocks> m_rests {};
        // Occupied bitvector
        BitArray<1UL << q> m_occupied {};
        // Run ends bitvector
        BitArray<1UL << q> m_runend {};
        // Offset vector : An offset value corresponds to the offset to the next slot that can corresponds to the current block.
        std::array<uint64_t, QuotientFilter::num_blocks> m_offsets {};

        uint64_t m_num_elements;

public:
    QuotientFilter() : m_num_elements(0)
    {
        static_assert(64 >= q + r);
        static_assert(q >= 6, "Minimum structure size is 64 elements => q must be larger or equal to 6");
        static_assert(r < 64, "Rest size must be < 64 to be stored in uint64_t arrays");
    };


    /** Construct a filter from a half size one.
     **/
    QuotientFilter(const QuotientFilter<q-1, r+1, Quotienting>& other) : m_num_elements(0)
    {
        if (r == 0)
            throw std::runtime_error("Cannot create a QF with 0 bit size slots.");
        if (other.m_num_elements == 0)
            throw std::runtime_error("Cannot resize an empty filter.");

        static constexpr uint64_t other_num_slots {QuotientFilter<q-1, r+1, Quotienting>::num_slots};

        // 0 - Get the position of the first run of the QF
        const uint64_t first_occ {other.m_occupied.first_one(0)};
        const uint64_t first_end {other.m_runend.first_one(other.m_offsets[0])};

        // 1 - Setup the first run for the loop
        uint64_t current_occ {first_occ};
        uint64_t current_end {first_end};
        uint64_t first_after_last_run {other.m_offsets[0]};
        do {
            // Decide which is the first slot of the run
            uint64_t slot_idx {std::max(current_occ, first_after_last_run)};
            // Only if we are at the edge of the filter
            if (current_end < current_occ)
                // If we are in the middle of the run after a toric loop
                if (first_after_last_run < current_end)
                    slot_idx = first_after_last_run;

            // 2 - Enumerates the rests in between run start and run end
            for ( ; slot_idx!=((current_end+1)%other_num_slots) ; slot_idx++ )
            {
                // Apply toricity
                if (slot_idx == other_num_slots)
                    slot_idx = 0;

                //                                                             slot_idx % 64
                const QR<q-1,r+1> element {current_occ, other.m_rests[slot_idx/64UL].get(slot_idx & 0b111111)};
                const uint64_t value {(Quotienting::template recompose<q-1,r+1>(element))};
                // TODO: Change that to be more efficient
                this->insert(value);
            }

            // 3 - Go to the next run
            first_after_last_run = (current_end + 1) % other_num_slots;
            current_occ = other.m_occupied.first_one((current_occ + 1) % other_num_slots);
            current_end = other.m_runend.first_one((current_end + 1) % other_num_slots);
        }
        while (current_occ != first_occ);
    }


    uint64_t size() const
    {
        return m_num_elements;
    }

    /** True if the maximum capacity of the QF has been reached
     **/
    bool max_capacity() const
    {
        return static_cast<double>(m_num_elements) / static_cast<double>(num_slots) > .95;
    }

    // uint64_t * insert(const uint64_t * begin, const uint64_t * end)
    // {
    //     for (uint64_t* element = begin; element != end; ++element)
    //     {
    //         // If max capacity reached stop the insertion here and return the next pointer
    //         if (max_capacity())
    //             return element;
    //         // Insert the element
    //         this->insert(*element);
    //     }

    //     return end;
    // }


    /** Insert an element in the RSQF.
     * @param element The element to insert that is splitted into quotient and value
     **/
    void insert(const uint64_t value)
    {
        // 0 - Compute the quotienting
        const QR<q,r> element {Quotienting::template compute<q,r>(value)};

        // 1 - Get the slot where to insert the new element
        const uint64_t insert_position {compute_insert_position(element)};

        // 2 - Compute the position where to shift the colided element
        const uint64_t free_slot {first_unused_slot(insert_position)};

        // 3 - Insert and swap using the previously computed coordinates.
        insert_and_shift(element, insert_position, free_slot);
        
        // 4 - Update the offsets changed before the insertion position
        const uint64_t insertion_block {insert_position / 64UL};
        uint64_t current_block {element.quotient / 64UL};

        while (current_block != insertion_block)
        {
            // Jump to the next block
            current_block += 1;
            if (current_block == num_blocks)
                current_block = 0;

            // We changed block so offset has to be incremented
            m_offsets[current_block] += 1;
        }

        m_num_elements += 1;
    }


    /** Check the presence of the value inside of the filter.
     * @param value The value to look for
     * @return true if the value has been inserted/
     **/
    bool get(const uint64_t value) const
    {
        // 0 - Compute the quotienting
        const QR<q,r> element = (Quotienting::template compute<q,r>(value));

        // 1 - Get the slot where the value should be
        const uint64_t slot {compute_insert_position(element)};

        // 2 - Is it present ?
        return m_rests[slot / 64UL].get(slot % 64UL) == element.rest;
    }


    uint64_t first_unused_slot(const uint64_t start_quotient) const
    {
        // using namespace std;
        // std::cout << "first_unused_slot " << start_quotient << std::endl;
        uint64_t current_quotient {start_quotient};
        uint64_t current_block {start_quotient / 64UL};
        uint64_t current_relative_quotient {start_quotient % 64};
        uint64_t first_block_quotient {current_block * 64UL};
        uint64_t current_absolute_offset {m_offsets[current_block] + first_block_quotient};

        // cout << current_quotient << " " << current_block << " " << first_block_quotient << endl;
        
        // Compute the number of runs started in this block up to the current position
        uint64_t local_occ_rank {m_occupied.rank(first_block_quotient, current_quotient)};
        uint64_t current_runend_idx {0};
        uint64_t next_block {current_block};
        
        if (local_occ_rank > 0) {
            // Jump to the corresponding runend by performing a select on opened runs
            current_runend_idx = m_runend.select(first_block_quotient + m_offsets[current_block], local_occ_rank);
        }
        else if (m_offsets[current_block] > current_relative_quotient)
        {
            // No opened runs in this block but remaining ones in the offset
            current_runend_idx = (current_absolute_offset - 1) % num_slots;
        }
        else
            return start_quotient;

        next_block = current_runend_idx / 64UL;

        // If last jump was negative, un empty slot has been found
        while ((next_block != current_block) or (current_runend_idx >= current_quotient))
        {
            // Init new values
            current_quotient = current_runend_idx + 1UL;
            current_block = current_quotient / 64UL;
            current_relative_quotient = current_quotient % 64;
            first_block_quotient = current_block * 64UL;
            current_absolute_offset = m_offsets[current_block] + first_block_quotient;

            // Compute the number of runs started in this block up to the current position
            local_occ_rank = m_occupied.rank(first_block_quotient, current_quotient);

            if (local_occ_rank > 0)
                // Jump to the runend corresponding to the number of runs in this block
                current_runend_idx = m_runend.select((first_block_quotient + m_offsets[current_block]) % num_slots, local_occ_rank);
            else if (m_offsets[current_block] > current_relative_quotient)
                // No opened runs in this block but remaining ones in the offset
                current_runend_idx = (current_absolute_offset - 1) % num_slots;
            else
                return current_quotient;

            next_block = current_runend_idx / 64UL;
        }

        return current_quotient;
    }

    /** Search for the position where an element must be inserted.
     * @return The absolute position in the filter where the element has to be inserted.
     **/
    uint64_t compute_insert_position(const QR<q, r>& element) const
    {
        const uint64_t block_idx {element.quotient / 64UL};
        const uint64_t first_block_quotient {block_idx * 64UL};

        // Rank the number of run stating in this block before the element theoritical position
        const uint64_t block_runs {element.quotient == 0 ? 0 : m_occupied.rank(first_block_quotient, element.quotient-1)};
        uint64_t real_insertion_slot {element.quotient};

        // If there are runs to jump over
        if (block_runs > 0)
        {
            real_insertion_slot = (m_runend.select(first_block_quotient + m_offsets[block_idx], block_runs) + 1UL) % num_slots;
        }
        // If there is an offset to jump over
        else if (m_offsets[block_idx] > (element.quotient % 64UL))
        {
            real_insertion_slot = m_offsets[block_idx];
        }
        // else: current quotient is the right place to go

        // Find the position in the run
        if (m_occupied.get(element.quotient))
        {
            // Is it the run end ?
            while (element.rest > m_rests[real_insertion_slot / 64UL].get(real_insertion_slot % 64UL))
            {
                // End of the run. Must insert here
                if (m_runend.get(real_insertion_slot))
                {
                    return (real_insertion_slot + 1) % num_slots;
                }

                // Not the right place for insertion
                real_insertion_slot += 1;
                if (real_insertion_slot == num_slots)
                    real_insertion_slot = 0;
            }
        }

        return real_insertion_slot;
    }


    /** Insert the rest in the QF at the address represented by the quotient. Do not test the slot availablity and
     * insert as if it was free.
     * @param element The quotient/rest split of the value.
     * @param insertion_slot Real slot index where to insert the rest (if the value is shifted)
     **/
    void insert_in_free_space(const QR<q, r>& element, const uint64_t insertion_slot)
    {
        // Set metadata for the new run
        m_occupied.set(element.quotient);
        m_runend.set(insertion_slot);

        // Set the rest
        const uint64_t block_idx {static_cast<uint64_t>(insertion_slot) >> 6};
        const uint64_t relative_idx {element.quotient & 0b111111UL};
        m_rests[block_idx].set(relative_idx, element.rest);
    }

    void insert_in_free_space(const QR<q, r>& element)
    {
        insert_in_free_space(element, element.quotient);
    }


    /** Insert the rest at the insertion_idx position. Performs a cascading shift to move the rests that can already be
     * in place.
     * @param rest the rest to insert.
     * @param insertion_idx Index where the rest must be inserted
     * @param first_free_slot Last slot where the function write the casading slots shifted.
     **/
    void insert_and_shift(const QR<q, r>& element, const uint64_t insertion_idx, const uint64_t first_free_slot)
    {
        uint64_t to_insert {element.rest};

        uint64_t current_idx {insertion_idx};
        uint64_t current_block_index {insertion_idx >> 6};
        uint64_t current_64_index {insertion_idx & 0b111111};
        
        PackedBlock<r>& current_block = m_rests[current_block_index];

        // --- Loop over the slots to shift everything until we reach an empty slot ---
        while (first_free_slot != current_idx)
        {
            // Rewrite the current slot
            const uint64_t save {current_block.get(current_64_index)};
            current_block.set(current_64_index, to_insert);
            to_insert = save;

            // Jump to the next slot
            current_idx += 1;
            current_64_index = (current_64_index + 1) & 0b111111;

            // Block change
            if (current_64_index == 0)
            {
                current_block_index += 1;
                // Pushing a new rest in that block => updating the offset
                m_offsets[current_block_index] += 1;

                // In casee we hit the right of the datastructure => toricity
                if ((1UL << q) == current_idx)
                {
                    current_idx = 0;
                    current_block_index = 0;
                }

                // Getting the current block
                current_block = m_rests[current_block_index];
            }
        }

        // --- Write the last value ---
        current_block.set(current_64_index, to_insert);

        // --- shift the runends that has been moved ---
        if (insertion_idx != first_free_slot)
            m_runend.toric_right_shift(insertion_idx, first_free_slot);

        // --- Set the occupied and metadata bits ---
        // New run inserted
        if (not m_occupied.get(element.quotient))
        {
            m_occupied.set(element.quotient);
            m_runend.set(insertion_idx);
        }
        // Inserted at the end of an existing run
        else if (insertion_idx == first_free_slot)
        {
            m_runend.unset((insertion_idx + num_slots - 1) % num_slots);
            m_runend.set(insertion_idx);
        }
    }
};

} // namespace dbglib

#endif // RSQF_HPP
