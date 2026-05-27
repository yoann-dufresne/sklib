#include <iostream>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>

#include <io/Skmer.hpp>


#ifndef SKMERATOR_H
#define SKMERATOR_H



namespace km
{


template <typename kuint>
class FileSkmerator;



// ----------------------- SEQ ITERATOR ---------------------------

template <typename kuint>
class SeqSkmerator
{
protected:
    std::string empty_str{};

    SkmerManipulator<kuint>& m_manip;
    std::string& m_seq;
    Skmer<kuint> m_yielded_skmer;

public:

    SeqSkmerator(SkmerManipulator<kuint>& skmer_manipulator, std::string& sequence)
        : m_manip(skmer_manipulator), m_seq(sequence)
    {};

    SeqSkmerator(SkmerManipulator<kuint>& skmer_manipulator)
        : SeqSkmerator(skmer_manipulator, empty_str)
    {};

    SeqSkmerator& operator=(const SeqSkmerator& other)
    {
        m_seq = other.m_seq;
        m_yielded_skmer = other.m_yielded_skmer;

        return *this;
    }

    Skmer<kuint>& current_skmer()
    {
        return m_yielded_skmer;
    }

    struct Iterator
    {

    private:
        std::string empty_str{};
        // Sequence file related attributes
        SeqSkmerator<kuint>& m_rator;
        std::string& m_seq;
        int64_t m_remaining_nucleotides;
        bool m_consumed;
        bool m_isend;

        // Manipulator to compute interleaved along the sequence
        SkmerManipulator<kuint>& m_manip;
        kuint m_current_minimizer;

        // Buffered unmasked skmers
        uint64_t m_buffer_size;
        std::vector<Skmer<kuint> > m_skmer_buffer_array;
        std::vector<bool> m_skmer_orientation;
        uint64_t m_ptr_current;
        uint64_t m_ptr_min;
        uint64_t m_ptr_last_round;

    public:

        bool consumed() const
        {
            return m_isend;
        }

    protected:
        // Construct an iterator without control on the file stream
        Iterator(SeqSkmerator& skmerator, std::string& sequence)
            : m_rator(skmerator), m_seq(sequence), m_remaining_nucleotides(sequence.length())
            , m_consumed(false), m_isend(false)
            , m_manip(skmerator.m_manip), m_current_minimizer(~static_cast<kuint>(0))
            , m_buffer_size(2 * skmerator.m_manip.k - skmerator.m_manip.m)
            , m_skmer_buffer_array(std::vector<Skmer<kuint> >(m_buffer_size))
            , m_skmer_orientation(std::vector<bool>(m_buffer_size))
            , m_ptr_current(skmerator.m_manip.k - skmerator.m_manip.m)
            , m_ptr_min(0) // minimizer_position
            , m_ptr_last_round(0)
        {
            if (m_remaining_nucleotides < static_cast<int64_t>(m_manip.k))
            {
                m_consumed = true;
                m_isend = true;
                return;
            }

            this->init_seq();
            this->operator++();
        }

        Iterator(SeqSkmerator& skmerator)
            : Iterator(skmerator, empty_str)

        {};

        Iterator& operator=(const Iterator& other)
        {
            m_rator = other.m_rator;
            m_seq = other.m_seq;
            m_remaining_nucleotides = other.m_remaining_nucleotides;
            m_consumed = other.m_consumed;
            m_isend = other.m_isend;
            m_current_minimizer = other.m_current_minimizer;

            m_buffer_size = other.m_buffer_size;
            m_skmer_buffer_array = std::vector(other.m_skmer_buffer_array);
            m_skmer_orientation = std::vector<bool>(other.m_skmer_orientation);

            m_ptr_current = other.m_ptr_current;
            m_ptr_min = other.m_ptr_min;
            m_ptr_last_round = other.m_ptr_last_round;

            return *this;
        }

        void init_seq()
        {
            // Init the first kmer
            for (uint64_t idx{0} ; idx<m_manip.k-1 ; idx++)
            {
                // Nucl encoding. TODO: Move encoding to dedicated classes
                const kuint nucl {static_cast<kuint>((m_seq[idx] >> 1) & 0b11U)};
                m_manip.add_nucleotide(nucl);
            }

            m_current_minimizer = ~static_cast<kuint>(0);

            m_ptr_current = m_rator.m_manip.k - 2;
            m_ptr_min = m_rator.m_manip.k - 2;
            m_remaining_nucleotides -= m_manip.k - 1;
        }

        friend class SeqSkmerator;
        friend class FileSkmerator<kuint>;

    public:
        // Return kmer by value
        Skmer<kuint> operator*() const
        {
            // cout << "operator*" << endl;
            return m_rator.m_yielded_skmer;
        }

        Iterator& operator++()
        {
            // cout << "operator++" << endl;
            if (m_consumed) {
                m_isend = true;
                // debug_print_buffer();
                return *this;
            }

            const uint64_t k {m_manip.k};
            const uint64_t m {m_manip.m};
            // km::SkmerPrettyPrinter<kuint> pp {k, m};

            // End of the sequence => final yieldings
            if (m_remaining_nucleotides + k - m <= 0)
            {
                // cout << "final yielding" << endl;
                // debug_print_buffer();
                // Yield the stored but not returned skmers while sequence is already over
                do
                {
                    m_ptr_last_round += 1;
                    // cout << "m_ptr_last_round " << (m_ptr_last_round % m_buffer_size);
                    // cout << " m_ptr_current " << (m_ptr_current % m_buffer_size) << endl;

                    // Get yielding candidate
                    Skmer<kuint>& skmer {m_skmer_buffer_array[m_ptr_last_round % m_buffer_size]};

                    // Yield if needed
                    if (skmer.m_pref_size + skmer.m_suff_size >= k - m)
                    {
                        m_rator.m_yielded_skmer = skmer;
                        m_manip.mask_absent_nucleotides(m_rator.m_yielded_skmer);
                        m_manip.canonicalize(m_rator.m_yielded_skmer);

                        // Avoid jumping around the final condition of the do-while
                        if ((m_ptr_last_round % m_buffer_size) == (m_ptr_current % m_buffer_size))
                            m_consumed = true;

                        return *this;
                    }
                }
                while ((m_ptr_last_round % m_buffer_size) != (m_ptr_current % m_buffer_size));

                // debug_print_buffer();
                m_consumed = true;
                m_isend = true;
                return *this;
            }

            bool yield_needed{false};
            // km::SkmerPrettyPrinter<kuint> pp {k, m};
            while (m_remaining_nucleotides + k - m > 0)
            {
                // -- Save the skmer to eventually yield
                m_rator.m_yielded_skmer = m_skmer_buffer_array[(m_ptr_current + 1) % m_buffer_size];
                // pp << m_rator.m_yielded_skmer;
                // cout << pp << " " << ((m_ptr_current + 1) % m_buffer_size);
                // cout << " remaining " << m_remaining_nucleotides << endl;

                // -- On out of context minimizer
                if (m_ptr_current - m_ptr_min >= k - m)
                {
                    this->solve_out_of_context();
                }

                // Get the next candidate skmer
                Skmer<kuint>& candidate {compute_new_candidate_skmer()};
                auto const candidate_orient {m_skmer_orientation[m_ptr_current % m_buffer_size]};

                // Get the minimizer
                const kuint candidate_minimizer {m_manip.minimizer(candidate)};

                // -- A new minimizer has been discovered
                if (candidate_minimizer < m_current_minimizer) {
                    this->new_minimizer(candidate_minimizer);
                }

                // The candidate minimizer is the same than the previous minimizer
                else if (candidate_minimizer == m_current_minimizer)
                {
                    this->update_on_equal_mini(m_ptr_current);
                }

                // The candidate is not a minimizer regarding the ongoing one.
                else
                {
                    auto const dist_to_min {m_ptr_current - m_ptr_min};
                    update_skmer_left_size(candidate, candidate_orient, dist_to_min-1);
                }

                // -- Yield if needed
                if (m_rator.m_yielded_skmer.m_pref_size + m_rator.m_yielded_skmer.m_suff_size >= k - m)
                {
                    m_manip.mask_absent_nucleotides(m_rator.m_yielded_skmer);
                    m_manip.canonicalize(m_rator.m_yielded_skmer);
                    if ((m_remaining_nucleotides + k - m) == 0) {
                        yield_needed = true;
                        break;
                    }
                    return *this;
                }

            } // End of the while over the sequence

            // Update the last skmers that are not selected to remove their right part.
            for(uint64_t idx{m_ptr_min+1} ; idx<=m_ptr_current ; idx++)
            {
                Skmer<kuint>& unused_skmer = m_skmer_buffer_array[idx % m_buffer_size];
                unused_skmer.m_pref_size = unused_skmer.m_suff_size = 0;
            }

            m_ptr_last_round = m_ptr_current;
            // Recursive call to return the already computed skmer array
            if (yield_needed)
                return *this;
            else
                return this->operator++();
        }

        void solve_out_of_context()
        {
            const uint64_t k {m_manip.k};
            const uint64_t m {m_manip.m};

            // Set the suffix size of the ooc skmer
            Skmer<kuint>& ooc_skmer {m_skmer_buffer_array[m_ptr_min % m_buffer_size]};
            auto const ooc_orientation {m_skmer_orientation[m_ptr_min % m_buffer_size]};
            update_skmer_right_size(ooc_skmer, ooc_orientation, k - m);

            this->recompute_minimizer(m_ptr_min + 1, m_ptr_current);
        }

        void new_minimizer(const kuint minimizer)
        {
            const uint64_t k {m_manip.k};
            const uint64_t m {m_manip.m};

            // Limit the size on the right of the previous skmers
            for (uint64_t idx{0} ; idx < (k - m) ; idx++)
            {
                auto const pred_buff_idx {(m_ptr_current - (k - m) + idx) % m_buffer_size};
                Skmer<kuint>& predecessor {m_skmer_buffer_array[pred_buff_idx]};
                auto const pred_orient {m_skmer_orientation[pred_buff_idx]};
                // Mask the nucleotides from the predecessor that are part of this skmer.
                auto const right_size {get_skmer_right_size(predecessor, pred_orient)};
                update_skmer_right_size(predecessor, pred_orient, std::min(right_size
                                                  , static_cast<uint16_t>(k - m - idx - 1)));
            }

            // 4 - save the new current minimal skmer and minimizer
            m_ptr_min = m_ptr_current;
            m_current_minimizer = minimizer;
        }


        Skmer<kuint>& compute_new_candidate_skmer()
        {
            const uint64_t k {m_manip.k};
            const uint64_t m {m_manip.m};

            // -- Compute the new candidate skmer
            m_remaining_nucleotides -= 1;
            m_ptr_current += 1;

            const kuint nucl {
                    m_remaining_nucleotides >= 0 ?(static_cast<kuint>((m_seq[m_ptr_current] >> 1) & 0b11U)) : static_cast<kuint>(0b11U)};

            m_skmer_buffer_array[ m_ptr_current % m_buffer_size ] = (m_remaining_nucleotides >= 0) ? m_manip.add_nucleotide(nucl) : m_manip.add_empty_nucleotide();
            m_skmer_orientation[ m_ptr_current % m_buffer_size ] = m_manip.is_forward();

            Skmer<kuint> & skmer = m_skmer_buffer_array[ m_ptr_current % m_buffer_size ];
            orientation_t const orient = m_skmer_orientation[ m_ptr_current % m_buffer_size ];

            // Sequence left extremity limitation
            if (m_seq.length() - m_remaining_nucleotides >= 2 * k - m)
                update_skmer_left_size(skmer, orient, k-m);
            else
                update_skmer_left_size(skmer, orient, m_seq.length() - m_remaining_nucleotides - k);

            // Sequence right extremity limitation
            if (m_remaining_nucleotides >= 0)
                update_skmer_right_size(skmer, orient, k-m);
            else
                update_skmer_right_size(skmer, orient, k-m+m_remaining_nucleotides);

            return skmer;
        }

        void update_on_equal_mini(uint64_t ptr_equivalent)
        {
            const uint64_t k{m_manip.k};
            const uint64_t m{m_manip.m};

            const uint64_t pos_diff {ptr_equivalent - m_ptr_min};
            Skmer<kuint>& current_skmer {m_skmer_buffer_array[ptr_equivalent % m_buffer_size]};
            auto const curr_orient {m_skmer_orientation[ptr_equivalent % m_buffer_size]};
            Skmer<kuint>& prev_skmer {m_skmer_buffer_array[m_ptr_min % m_buffer_size]};
            auto const prev_orient {m_skmer_orientation[m_ptr_min % m_buffer_size]};


            // Define the kmers to compare
            uint64_t const start {std::max(pos_diff, k-m-get_skmer_left_size(prev_skmer, prev_orient))};
            uint64_t const stop {get_skmer_right_size(prev_skmer, prev_orient)};
            // Init outside of the first kmer (true if all the kmers are from the first skmer)
            uint64_t first_index_in_new_skmer {stop+1};

            // 1 - Scan all the kmer shared between the 2 skmers to decide where is the limit between the 2 skmers
            for (uint64_t i{start} ; i<=stop ; i++)
            {
                uint64_t const prev_skmer_km_pos {prev_orient == forward_c ? i : k - m - i};
                uint64_t const curr_skmer_km_pos {curr_orient == forward_c
                                                        ? i - pos_diff
                                                        : k - m - i + pos_diff
                                                 };

                // Get out of the loop after the changing point is reached
                if (m_manip.kmer_lt_kmer(current_skmer, curr_skmer_km_pos, prev_skmer, prev_skmer_km_pos))
                {
                    first_index_in_new_skmer = i;
                    break;
                }
            }

            // Invalidates skmers in between the 2 equivalent minimizers
            for (uint64_t idx{m_ptr_min+1} ; idx<ptr_equivalent ; idx++)
            {
                Skmer<kuint>& invalid {m_skmer_buffer_array[idx % m_buffer_size]};
                invalid.m_pref_size = invalid.m_suff_size = 0;
            }

            // Close the previous skmer
            update_skmer_right_size(prev_skmer, prev_orient, first_index_in_new_skmer - 1);
            update_skmer_left_size(current_skmer, curr_orient, k-m+pos_diff-first_index_in_new_skmer);

            m_ptr_min = ptr_equivalent;
        }

        /** Compute the dominant minimizer in the slice of the buffer. If multiple identical minimizers are present
         * it will initialize all the prefixes and suffixes of skmers needed to reach the last one.
         * @param buff_start First buffer coordinate to analyse (absolute value)
         * @param buff_stop Last buffer coordinate to analyse (absolute value)
         **/
        void recompute_minimizer(uint64_t buff_start, uint64_t buff_stop)
        {
            // 0 - Select first candidate
            Skmer<kuint>& first = m_skmer_buffer_array[buff_start % m_buffer_size];
            m_current_minimizer = m_manip.minimizer(first);
            m_ptr_min = buff_start;
            uint64_t first_equal_mini {m_ptr_min};
            uint64_t last_equal_mini {m_ptr_min};

            // 1 - Compute the minimizer in the window
            for (uint64_t idx{buff_start+1} ; idx<=buff_stop ; idx++)
            {
                Skmer<kuint>& candidate = m_skmer_buffer_array[idx % m_buffer_size];
                kuint candidate_minimizer = m_manip.minimizer(candidate);

                if (candidate_minimizer < m_current_minimizer)
                {
                    m_current_minimizer = candidate_minimizer;
                    first_equal_mini = idx;
                    last_equal_mini = idx;
                }
                else if (candidate_minimizer == m_current_minimizer)
                {
                    last_equal_mini = idx;
                }
            }

            for (uint64_t idx{buff_start} ; idx<=buff_stop ; idx++)
            {
                Skmer<kuint>& candidate = m_skmer_buffer_array[idx % m_buffer_size];

                // Before the last equivalent minimizer we can detroy the candidate skmers with larger minimizers
                if (idx < first_equal_mini)
                {
                    candidate.m_pref_size = candidate.m_suff_size = 0;
                }
                else if (idx <= last_equal_mini)
                {
                    kuint candidate_minimizer = m_manip.minimizer(candidate);

                    if (candidate_minimizer == m_current_minimizer)
                    {
                        Skmer<kuint> & mini_skmer {m_skmer_buffer_array[m_ptr_min % m_buffer_size]};
                        auto const mini_orient {m_skmer_orientation[m_ptr_min % m_buffer_size]};
                        // cout << mini_skmer << endl;

                        update_skmer_left_size(mini_skmer, mini_orient, m_ptr_min - buff_start);
                        if (idx != first_equal_mini)
                            this->update_on_equal_mini(idx);
                        else
                            m_ptr_min = idx;
                        // cout << mini_skmer << endl << endl;
                    }
                }
                // We now only have to update the left part of remaining skmers
                else
                {
                    auto const orient {m_skmer_orientation[idx % m_buffer_size]};
                    update_skmer_left_size(candidate, orient, idx - last_equal_mini - 1);
                }
            }
        }


        // Warning: This function suppose that we are comparing iterator over the same sequence.
        bool operator==(const Iterator& it) const
        {
            if (m_isend and it.m_isend)
                return true;
            else if (m_isend != it.m_isend)
                return false;
                
            return m_remaining_nucleotides == it.m_remaining_nucleotides;
        }

        void debug_print_buffer()
        {
            const uint64_t k{m_manip.k};
            const uint64_t m{m_manip.m};
            km::SkmerPrettyPrinter<kuint> pp {k, m};

            for (uint64_t i{0} ; i<m_buffer_size ; i++)
            {
                pp << m_skmer_buffer_array[i];
                std::cout << i << ":\t" << m_skmer_buffer_array[i] << " orientation:" << m_skmer_orientation[i] << "\t" << pp << std::endl;
            }
        }

        void update_skmer_right_size(Skmer<kuint>& skmer, orientation_t orientation, uint64_t size)
        {
            if (orientation == forward_c)
                skmer.m_suff_size = size;
            else
                skmer.m_pref_size = size;
        }

        void update_skmer_left_size(Skmer<kuint>& skmer, orientation_t orientation, uint64_t size)
        {
            if (orientation == forward_c)
                skmer.m_pref_size = size;
            else
                skmer.m_suff_size = size;
        }

        uint16_t get_skmer_left_size(Skmer<kuint> const & skmer, orientation_t orientation) const
        {
            if (orientation == forward_c)
                return skmer.m_pref_size;
            else
                return skmer.m_suff_size;
        }

        uint16_t get_skmer_right_size(Skmer<kuint> const & skmer, orientation_t orientation) const
        {
            if (orientation == forward_c)
                return skmer.m_suff_size;
            else
                return skmer.m_pref_size;
        }

    };

    Iterator begin() { return Iterator(*this, m_seq); }
    Iterator end() { return Iterator(*this); }
};


// ----------------------- FILE ITERATOR ---------------------------


template <typename kuint>
class FileSkmerator
{
protected:
    std::string m_filename;
    SkmerManipulator<kuint>& m_manip;

public:
    FileSkmerator(SkmerManipulator<kuint>& manipulator, const std::string& filename)
        : m_filename(std::string{filename}), m_manip(manipulator)
    {};

    FileSkmerator<kuint>& operator=(const FileSkmerator<kuint>& other)
    {
        m_filename = other.m_filename;
        m_manip = other.m_manip;

        return *this;
    }

    struct Iterator
    {

    private:
        // Sequence file related attributes
        FileSkmerator<kuint>& m_rator;
        std::unique_ptr<klibpp::SeqStreamIn> m_ptr;
        klibpp::KSeq m_record;
        // skmer sequence enumerator
        SeqSkmerator<kuint> m_seq_rator;
        typename SeqSkmerator<kuint>::Iterator m_skmer_iterator;

    protected:
        // Construct an iterator without control on the file stream
        Iterator(FileSkmerator<kuint>& skmerator, std::unique_ptr<klibpp::SeqStreamIn> stream_ptr)
            : m_rator(skmerator), m_ptr(std::move(stream_ptr))
            , m_seq_rator(skmerator.m_manip), m_skmer_iterator(m_seq_rator.begin())
        {
            if (m_ptr == nullptr)
                return;

            this->init_record();
        }

        // Construct a new file stream from the filename
        Iterator(FileSkmerator<kuint>& skmerator)
            : Iterator ( skmerator, std::make_unique<klibpp::SeqStreamIn>(skmerator.m_filename.c_str()) )
        {}

        void init_record()
        {
            do
            {
                if ((*m_ptr) >> m_record)
                {
                    // If sequence is too short, skip it
                    if (m_record.seq.length() < m_rator.m_manip.k) {
                        continue;
                    }

                    m_seq_rator = SeqSkmerator<kuint>(m_rator.m_manip, m_record.seq);
                    m_skmer_iterator = m_seq_rator.begin();

                }
                else
                {
                    // No more sequence to read from the file
                    m_ptr = nullptr;
                    return;
                }
            }
            while (m_skmer_iterator.consumed());
        }

        friend class FileSkmerator;

    public:
        // Return kmer by value
        Skmer<kuint> operator*()
        {
            return m_seq_rator.current_skmer();
        }

        Iterator& operator++()
        {
            m_skmer_iterator.operator++();

            // File already consumed
            if (m_ptr == nullptr)
                return *this;

            // Go to next sequence
            if (m_skmer_iterator.consumed())
            {
                this->init_record();

                // reached the end of the file while looking for the next sequence
                if (m_ptr == nullptr)
                    return *this;
            }

            return *this;
        }

        bool operator==(const Iterator& it) const
        {
            return m_rator.m_filename == it.m_rator.m_filename and m_ptr == it.m_ptr;
        }

    };


    Iterator begin() { return Iterator(*this); }
    Iterator end() { return Iterator(*this, nullptr); }
};


}
#endif