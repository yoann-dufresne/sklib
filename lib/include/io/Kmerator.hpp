#include <iostream>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>

#include <io/Kmer.hpp>


#ifndef KMERATOR_H
#define KMERATOR_H

using namespace std;


namespace km
{


template <typename kuint>
class FileKmerator
{
protected:
    const std::string& m_filename;
    KmerManipulator<kuint>& m_manip;
    Kmer<kuint> m_current_kmer;

public:
    FileKmerator(const std::string& filename, KmerManipulator<kuint>& manipulator)
        : m_filename(filename), m_manip(manipulator)
    {};

    struct Iterator
    {
    protected:
        // Construct an iterator without control on the file stream
        Iterator(FileKmerator& kmerator, std::unique_ptr<klibpp::SeqStreamIn> stream_ptr)
            : m_rator(kmerator), m_manip(kmerator.m_manip), m_ptr(std::move(stream_ptr))
            , m_remaining_kmers(0), m_seq_idx(0)
        {
            if (m_ptr == nullptr)
                return;

            this->init_record();

            const kuint nucl {(m_record.seq[m_seq_idx] >> 1) & 0b11U};
            m_rator.m_current_kmer = m_manip.add_nucleotide(nucl);

            m_remaining_kmers -= 1;
            m_seq_idx += 1;
        }
        // Construct a new file stream from the filename
        Iterator(FileKmerator& kmerator)
            : Iterator ( kmerator, std::make_unique<klibpp::SeqStreamIn>(kmerator.m_filename.c_str()) )
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

                    // Init the first kmer
                    m_manip.init_kmer();
                    for (m_seq_idx = 0 ; m_seq_idx<m_manip.k-1 ; m_seq_idx++)
                    {
                        // Nucl encoding. TODO: Move encoding to dedicated classes
                        const kuint nucl {(m_record.seq[m_seq_idx] >> 1) & 0b11U};
                        m_manip.add_nucleotide(nucl);
                    }

                    m_remaining_kmers = m_record.seq.length() - m_manip.k + 1;

                }
                else
                {
                    // No more sequence to read from the file
                    m_ptr = nullptr;
                    return;
                }
            }
            while (m_record.seq.length() < m_manip.k);
        }

        friend class FileKmerator;

    public:
        // Return kmer by value
        Kmer<kuint> operator*() const
        {
            return m_rator.m_current_kmer;
        }

        Iterator& operator++()
        {
            // File already consumed
            if (m_ptr == nullptr)
                return *this;

            // Go to next sequence
            if (m_remaining_kmers == 0)
            {
                this->init_record();
            }

            // Get the next kmer
            const kuint nucl {(m_record.seq[m_seq_idx] >> 1) & 0b11U};
            m_rator.m_current_kmer = m_manip.add_nucleotide(nucl);

            m_remaining_kmers -= 1;
            m_seq_idx += 1;

            return *this;
        }

        bool operator==(const Iterator& it) const
        {
            return m_rator.m_filename == it.m_rator.m_filename and m_ptr == it.m_ptr;
        }


    private:
        FileKmerator<kuint>& m_rator;
        KmerManipulator<kuint>& m_manip;
        std::unique_ptr<klibpp::SeqStreamIn> m_ptr;

        klibpp::KSeq m_record;
        uint64_t m_remaining_kmers;
        uint64_t m_seq_idx;
    };
    

    Iterator begin() { return Iterator(*this); }
    Iterator end() { return Iterator(*this, nullptr); }
};



};


#endif
