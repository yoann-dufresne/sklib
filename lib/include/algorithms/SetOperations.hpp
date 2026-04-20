#pragma once

#include <array>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <SortedVirtualSkmerList.hpp> // your existing header

namespace km
{
  namespace sortedlist
  {
    namespace setops
    {

      // ─────────────────────────────────────────────
      //  OPERATION TYPE
      // ─────────────────────────────────────────────

      enum class OpType
      {
        UNION,
        INTERSECTION,
        DIFF
      };

      // ─────────────────────────────────────────────
      //  IN-MEMORY SOURCE
      // ─────────────────────────────────────────────

      template <typename kuint>
      class InMemorySkmerSource
      {
      public:
        explicit InMemorySkmerSource(const std::vector<Skmer<kuint>> &list, uint64_t k, uint64_t m) : m_list(list), m_index(0), m_k(k), m_m(m) {}

        bool has_next() const
        {
          return m_index < m_list.size();
        }
        const Skmer<kuint> &current() const
        {
          return m_list[m_index];
        }
        void advance()
        {
          m_index++;
        }

        uint64_t k() const { return m_k; }
        uint64_t m() const { return m_m; }

      private:
        const std::vector<Skmer<kuint>> &m_list;
        size_t m_index;
        uint64_t m_k;
        uint64_t m_m;
      };

      // ─────────────────────────────────────────────
      //  FILE SOURCE
      // ─────────────────────────────────────────────

      template <typename kuint>
      class FileSkmerSource
      {
        static constexpr size_t BUFFER_SIZE = 512;

      public:
        explicit FileSkmerSource(const std::string &filepath)
            : m_buf_index(0), m_buf_valid(0), m_file_exhausted(false)
        {
          m_file.open(filepath, std::ios::binary);
          if (m_file.fail())
            throw std::runtime_error("Cannot open file: " + filepath);

          auto header = internal::read_and_validate_header<kuint>(m_file, filepath);
          m_k = header.k;
          m_m = header.m;
          m_total_count = header.count;
          // reuse logic from VirtualSkmerSerializer::load()
          // but DO NOT read skmers yet
          fill_buffer(); // read first BUFFER_SIZE skmers from file
        }
        ~FileSkmerSource()
        {
          if (m_file.is_open())
            m_file.close();
        }

        bool has_next() const
        {
          return m_buf_index < m_buf_valid;
        }

        const Skmer<kuint> &current() const
        {
          return m_buffer[m_buf_index];
        }

        void advance()
        {
          m_buf_index++;
          if (m_buf_index == m_buf_valid && !m_file_exhausted)
            fill_buffer();
        }

        uint64_t k() const { return m_k; }
        uint64_t m() const { return m_m; }

      private:
        void fill_buffer()
        {
          m_buf_index = 0;
          const uint64_t remaining = m_total_count - m_total_read;
          const uint64_t to_read = std::min(remaining, static_cast<uint64_t>(BUFFER_SIZE));

          m_file.read(reinterpret_cast<char *>(m_buffer.data()),
                      to_read * sizeof(Skmer<kuint>));

          m_buf_valid = m_file.gcount() / sizeof(Skmer<kuint>);
          m_total_read += m_buf_valid;
          m_file_exhausted = (m_total_read >= m_total_count);
        }

        std::ifstream m_file;
        std::array<Skmer<kuint>, BUFFER_SIZE> m_buffer;
        size_t m_buf_index;
        size_t m_buf_valid; // number of valid elements in buffer
        bool m_file_exhausted;
        uint64_t m_k;
        uint64_t m_m;
      };

      // ─────────────────────────────────────────────
      //  KMER CURSOR
      // ─────────────────────────────────────────────
      // Used to wrap a superkmer in the comparison between 2 in the set operations
      template <typename kuint>
      class KmerCursor
      {
      public:
        KmerCursor(const Skmer<kuint> &skmer, const SkmerManipulator<kuint> &manip) : m_skmer(skmer), m_manip(manip)
        {
          auto [start, end] = manip.get_valid_kmer_bounds(skmer);
          m_current_pos = start;
          m_end_pos = end;

          bool has_next_kmer() const
          {
            return m_current_pos <= m_end_pos;
          }
          uint64_t current_kmer_pos() const
          {
            return m_current_pos;
          }
          const Skmer<kuint> &skmer() const { return m_skmer; }
          void advance()
          {
            if (has_next_kmer())
            {
              m_current_pos++;
            }
          }

        private:
          const Skmer<kuint> &m_skmer;
          const SkmerManipulator<kuint> &m_manip;
          uint64_t m_current_pos;
          uint64_t m_end_pos;
        };

        // ─────────────────────────────────────────────
        //  MERGE ENGINE  (internal)
        // ─────────────────────────────────────────────

        namespace internal
        {

          // Handles the partial overlap (same minimizer, different pref/suff extent)
          // between two skmers at kmer granularity
          template <typename kuint>
          template <typename kuint>
          void merge_partial_skmers(
              const std::vector<Skmer<kuint>> &group_a,
              const std::vector<Skmer<kuint>> &group_b,
              const SkmerManipulator<kuint> &manip,
              std::vector<Skmer<kuint>> &output,
              OpType op)
          {
            // k-mer level carved skmer building + generate_sorted_list_from_enumeration

          }
          // Core merge-scan engine
          template <typename kuint, typename SourceA, typename SourceB>
          std::vector<Skmer<kuint>> merge_set_op(
              SourceA &source_a,
              SourceB &source_b,
              const SkmerManipulator<kuint> &manip,
              OpType op)
          {

            validate_sources(source_a, source_b);
            std::vector<Skmer<kuint>> output_vector;

            // pre-allocating vectors for slow-path to reduce reallocations
            std::vector<Skmer<kuint>> group_a;
            std::vector<Skmer<kuint>> group_b;

            while (source_a.has_next() && source_b.has_next())
            {
              const Skmer<kuint> &skmer_a = source_a.current();
              const Skmer<kuint> &skmer_b = source_b.current();
              kuint min_a = manip.minimizer(skmer_a);
              kuint min_b = manip.minimizer(skmer_b);

              if (min_a < min_b)
              {
                if (op == OpType::UNION || op == OpType::DIFF)
                  output_vector.push_back(skmer_a);
                source_a.advance();
              }
              else if (min_a > min_b)
              {
                if (op == OpType::UNION)
                  output_vector.push_back(skmer_b);
                source_b.advance();
              }
              else
              {
                group_a.clear();
                group_b.clear();

                while (source_a.has_next() && manip.minimizer(source_a.current()) == min_a)
                {
                  group_a.push_back(source_a.current());
                  source_a.advance();
                }
                while (source_b.has_next() && manip.minimizer(source_b.current()) == min_b)
                {
                  group_b.push_back(source_b.current());
                  source_b.advance();
                }

                merge_partial_skmers(group_a, group_b, manip, output_vector, op);
              }
            }

            if (op == OpType::UNION || op == OpType::DIFF)
              drain_source<kuint>(source_a, output_vector);
            if (op == OpType::UNION)
              drain_source<kuint>(source_b, output_vector);

            return output_vector;
          }

          // Validates that k and m match between two sources
          template <typename SourceA, typename SourceB>
          void validate_sources(const SourceA &a, const SourceB &b)
          {
            if (a.k() != b.k() || a.m() != b.m())
            {
              throw std::invalid_argument(
                  "Sources have incompatible k/m values: "
                  "k=" +
                  std::to_string(a.k()) + " vs k=" + std::to_string(b.k()) + ", "
                                                                             "m=" +
                  std::to_string(a.m()) + " vs m=" + std::to_string(b.m()));
            }
          }

          // Drains remaining elements of a source into the builder (for UNION / DIFF)
          template <typename kuint, typename Source>
          void drain_source(
              Source &source,
              std::vector<Skmer<kuint>> &output)
          { // direct append, no repacking needed
          }
        } // namespace internal

        // ─────────────────────────────────────────────
        //  PUBLIC API  — both in memory
        // ─────────────────────────────────────────────

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_union(
            const SortedVirtualSkmerList<kuint> &a,
            const SortedVirtualSkmerList<kuint> &b);

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_intersection(
            const SortedVirtualSkmerList<kuint> &a,
            const SortedVirtualSkmerList<kuint> &b);

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_diff(
            const SortedVirtualSkmerList<kuint> &a,
            const SortedVirtualSkmerList<kuint> &b);

        // ─────────────────────────────────────────────
        //  PUBLIC API  — one in memory, one from file
        // ─────────────────────────────────────────────

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_union(
            const SortedVirtualSkmerList<kuint> &a,
            const std::string &filepath_b);

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_intersection(
            const SortedVirtualSkmerList<kuint> &a,
            const std::string &filepath_b);

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_diff(
            const SortedVirtualSkmerList<kuint> &a,
            const std::string &filepath_b);

        // ─────────────────────────────────────────────
        //  PUBLIC API  — both from file
        // ─────────────────────────────────────────────

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_union(
            const std::string &filepath_a,
            const std::string &filepath_b,
            uint64_t k,
            uint64_t m);

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_intersection(
            const std::string &filepath_a,
            const std::string &filepath_b,
            uint64_t k,
            uint64_t m);

        template <typename kuint>
        SortedVirtualSkmerList<kuint> set_diff(
            const std::string &filepath_a,
            const std::string &filepath_b,
            uint64_t k,
            uint64_t m);

      } // namespace setops
    } // namespace sortedlist
  } // namespace km