#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

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
      //  IN-MEMORY SOURCE for public APIs
      // ─────────────────────────────────────────────

      template <typename kuint>
      class InMemorySkmerSource
      {
      public:
        explicit InMemorySkmerSource(const std::vector<Skmer<kuint>> &list,
                                     uint64_t k, uint64_t m)
            : m_list(list), m_index(0), m_k(k), m_m(m)
        {
        }

        bool has_next() const { return m_index < m_list.size(); }
        const Skmer<kuint> &current() const { return m_list[m_index]; }
        void advance() { ++m_index; }

        uint64_t k() const { return m_k; }
        uint64_t m() const { return m_m; }

      private:
        const std::vector<Skmer<kuint>> &m_list;
        size_t m_index;
        uint64_t m_k;
        uint64_t m_m;
      };

      // ─────────────────────────────────────────────
      //  FILE SOURCE for public APIs, uses buffer for lower memory footprint
      // ─────────────────────────────────────────────

      template <typename kuint>
      class FileSkmerSource
      {
        static constexpr size_t BUFFER_SIZE = 512;

      public:
        explicit FileSkmerSource(const std::string &filepath)
            : m_buf_index(0), m_buf_valid(0), m_file_exhausted(false),
              m_total_read(0), m_k(0), m_m(0), m_total_count(0)
        {
          m_file.open(filepath, std::ios::binary);
          if (m_file.fail())
            throw std::runtime_error("Cannot open file: " + filepath);

          auto header = internal::read_and_validate_header<kuint>(m_file, filepath);
          m_k = header.k;
          m_m = header.m;
          m_total_count = header.count;
          fill_buffer(); // read first chunk
        }

        ~FileSkmerSource()
        {
          if (m_file.is_open())
            m_file.close();
        }

        bool has_next() const { return m_buf_index < m_buf_valid; }
        const Skmer<kuint> &current() const { return m_buffer[m_buf_index]; }

        void advance()
        {
          ++m_buf_index;
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

          m_buf_valid = static_cast<size_t>(m_file.gcount()) / sizeof(Skmer<kuint>);
          m_total_read += m_buf_valid;
          if (m_total_read >= m_total_count)
            m_file_exhausted = true;
        }

        std::ifstream m_file;
        std::array<Skmer<kuint>, BUFFER_SIZE> m_buffer;
        size_t m_buf_index;
        size_t m_buf_valid;
        bool m_file_exhausted;
        uint64_t m_total_read;
        uint64_t m_k;
        uint64_t m_m;
        uint64_t m_total_count;
      };

      // ─────────────────────────────────────────────
      //  MERGE ENGINE
      // ─────────────────────────────────────────────

      namespace internal
      {
        // Handles the partial overlap (same minimizer, different pref/suff extent)
        // between two skmers at k-mer granularity
        template <typename kuint>
        void merge_partial_skmers(const std::vector<Skmer<kuint>> &group_a,
                                  const std::vector<Skmer<kuint>> &group_b,
                                  SkmerManipulator<kuint> &manip,
                                  std::vector<Skmer<kuint>> &output,
                                  OpType op)
        {
          if (op == OpType::UNION)
          {
            // By reusing generate_sorted_list_from_enumeration we ensure unicity. It does not use much the sorting property though
            std::vector<Skmer<kuint>> combined;
            auto before = group_a.size() + group_b.size();
            combined.reserve(group_a.size() + group_b.size());
            for (const auto &sk : group_a)
              combined.push_back(sk);
            for (const auto &sk : group_b)
              combined.push_back(sk);
            std::cout << "--- UNION-SORTED-LIST-GEN ---" << std::endl;
            SortedVirtualSkmerList<kuint> tmp(manip.k, manip.m);
            tmp.generate_sorted_list_from_enumeration(combined);
            auto after = tmp.size();
            if (after != group_a.size()){
              std::cerr << "Passed from " << before << " to " << after << " skmer." << std::endl;
              km::SkmerPrettyPrinter<kuint> pp {manip.k, manip.m};
              std::cout << "--- Fed in ---" << std::endl;
              for(km::Skmer<kuint> const skmer : combined){
                pp << skmer;
                std:: cout << pp << std::endl;
              }
              std::cout << "--- Got out ---" << std::endl;
              for(km::Skmer<kuint> const skmer : tmp.get_list()){
                pp << skmer;
                std:: cout << pp << std::endl;
              }
            }
            auto partial = tmp.steal_list();
            output.insert(output.end(),
                          std::make_move_iterator(partial.begin()),
                          std::make_move_iterator(partial.end()));
          }
          else
          {
            // INTERSECTION or DIFF: explode to singletons, filter, rebuild.
            std::vector<Skmer<kuint>> survivors;
            const uint64_t max_pos = manip.k - manip.m;
            // For every k-mer position I scan through the two lists and insert in a "survivors" group only the k-mers to be retained for new sorted virtual skmer list generation, based on diff or intersection.
            for (uint64_t p = 0; p <= max_pos; ++p)
            {
              size_t ia = 0, ib = 0;

              auto advance_a = [&]()
              {
                while (ia < group_a.size() && !manip.has_valid_kmer(group_a[ia], p))
                  ++ia;
              };
              auto advance_b = [&]()
              {
                while (ib < group_b.size() && !manip.has_valid_kmer(group_b[ib], p))
                  ++ib;
              };

              advance_a();
              advance_b();

              while (ia < group_a.size() && ib < group_b.size())
              {
                int cmp = manip.kmer_compare(group_a[ia], group_b[ib], p);

                if (cmp == 0)
                {
                  // k-mer present in both sides
                  if (op == OpType::INTERSECTION)
                    survivors.push_back(manip.get_skmer_of_kmer(group_a[ia], p));
                  // for DIFF: discard it
                  ++ia;
                  advance_a();
                  ++ib;
                  advance_b();
                }
                else if (cmp < 0)
                {
                  // A's k-mer is smaller → not in B
                  if (op == OpType::DIFF)
                    survivors.push_back(manip.get_skmer_of_kmer(group_a[ia], p));
                  ++ia;
                  advance_a();
                }
                else
                { // cmp > 0
                  // B's k-mer is smaller → not a match, move B forward
                  ++ib;
                  advance_b();
                }
              }

              // A has remaining valid k-mers at this position and B is exhausted
              if (op == OpType::DIFF)
              {
                while (ia < group_a.size())
                {
                  if (manip.has_valid_kmer(group_a[ia], p))
                    survivors.push_back(manip.get_skmer_of_kmer(group_a[ia], p));
                  ++ia;
                }
              }
            }

            // Re-compact the surviving singletons back into optimal skmers
            if (!survivors.empty())
            {
              SortedVirtualSkmerList<kuint> tmp(manip.k, manip.m);
              tmp.generate_sorted_list_from_enumeration(survivors);
              auto partial = tmp.steal_list();
              output.insert(output.end(),
                            std::make_move_iterator(partial.begin()),
                            std::make_move_iterator(partial.end()));
            }
          }
        }

        // Drains remaining elements of a source into the output vector
        template <typename kuint, typename Source>
        void drain_source(Source &source, std::vector<Skmer<kuint>> &output)
        {
          while (source.has_next())
          {
            output.push_back(source.current());
            source.advance();
          }
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
                std::to_string(a.k()) + " vs k=" + std::to_string(b.k()) +
                ", m=" + std::to_string(a.m()) + " vs m=" + std::to_string(b.m()));
          }
        }

        // Core merge-scan engine
        template <typename kuint, typename SourceA, typename SourceB>
        SortedVirtualSkmerList<kuint> merge_set_op(SourceA &source_a,
                                                   SourceB &source_b,
                                                   OpType op)
        {
          validate_sources(source_a, source_b);

          const uint64_t k = source_a.k();
          const uint64_t m = source_a.m();
          SkmerManipulator<kuint> manip(k, m);

          std::vector<Skmer<kuint>> output;
          std::vector<Skmer<kuint>> group_a;
          std::vector<Skmer<kuint>> group_b;

          while (source_a.has_next() && source_b.has_next())
          {
            // std::cout << "merge_loop" << std::endl;
            const Skmer<kuint> &sk_a = source_a.current();
            const Skmer<kuint> &sk_b = source_b.current();

            auto min_a = manip.minimizer(sk_a);
            auto min_b = manip.minimizer(sk_b);

            if (min_a < min_b)
            {
              // std::cout << "min_sk_a < min_sk_b" << std::endl;
              if (op != OpType::INTERSECTION)
                output.push_back(sk_a);
              source_a.advance();
            }
            else if (min_a > min_b)
            {
              // std::cout << "min_sk_a > min_sk_b" << std::endl;
              if (op == OpType::UNION)
                output.push_back(sk_b);
              source_b.advance();
            }
            else
            {
              // Same minimizer: collect entire runs from both sides
              // std::cout << "min_sk_a == min_sk_b" << std::endl;
              auto current_min = min_a;
              group_a.clear();
              group_b.clear();

              while (source_a.has_next() &&
                     manip.minimizer(source_a.current()) == current_min)
              {
                group_a.push_back(source_a.current());
                source_a.advance();
              }
              while (source_b.has_next() &&
                     manip.minimizer(source_b.current()) == current_min)
              {
                group_b.push_back(source_b.current());
                source_b.advance();
              }

              merge_partial_skmers(group_a, group_b, manip, output, op);
            }
          }

          // Drain whatever remains
          if (op == OpType::UNION || op == OpType::DIFF){
            // std::cout << "UNION/DIFF, draining A" << std::endl;
            drain_source(source_a, output);
          }
          if (op == OpType::UNION){
            // std::cout << "union, draining B" << std::endl;
            drain_source(source_b, output);
          }


          SortedVirtualSkmerList<kuint> result(k, m);
          result.add_list(std::move(output));
          return result;
        }

      } // namespace internal

      // ─────────────────────────────────────────────
      //  PUBLIC API  — both in memory
      // ─────────────────────────────────────────────

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_union(const SortedVirtualSkmerList<kuint> &a,
                                              const SortedVirtualSkmerList<kuint> &b)
      {
        InMemorySkmerSource<kuint> src_a(a.get_list(), a.get_k(), a.get_m());
        InMemorySkmerSource<kuint> src_b(b.get_list(), b.get_k(), b.get_m());
        return internal::merge_set_op<kuint>(src_a, src_b, OpType::UNION);
      }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_intersection(const SortedVirtualSkmerList<kuint> &a,
                                                     const SortedVirtualSkmerList<kuint> &b)
      {
        InMemorySkmerSource<kuint> src_a(a.get_list(), a.get_k(), a.get_m());
        InMemorySkmerSource<kuint> src_b(b.get_list(), b.get_k(), b.get_m());
        return internal::merge_set_op<kuint>(src_a, src_b, OpType::INTERSECTION);
      }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_diff(const SortedVirtualSkmerList<kuint> &a,
                                             const SortedVirtualSkmerList<kuint> &b)
      {
        InMemorySkmerSource<kuint> src_a(a.get_list(), a.get_k(), a.get_m());
        InMemorySkmerSource<kuint> src_b(b.get_list(), b.get_k(), b.get_m());
        return internal::merge_set_op<kuint>(src_a, src_b, OpType::DIFF);
      }

      // ─────────────────────────────────────────────
      //  PUBLIC API  — in-memory / file
      // ─────────────────────────────────────────────

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_union(const SortedVirtualSkmerList<kuint> &mem,
                                              const std::string &filepath)
      {
        FileSkmerSource<kuint> src(filepath);
        if (mem.get_k() != src.k() || mem.get_m() != src.m())
          throw std::invalid_argument("In-memory list k/m mismatch with file header");
        std::cout << "K is: " << mem.get_k() << ", M is: " <<  mem.get_m() << std::endl;
        InMemorySkmerSource<kuint> mem_src(mem.get_list(), mem.get_k(), mem.get_m());
        return internal::merge_set_op<kuint>(mem_src, src, OpType::UNION);
      }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_union(const std::string &filepath,
                                              const SortedVirtualSkmerList<kuint> &mem)
      { return set_union(mem, filepath); }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_intersection(const SortedVirtualSkmerList<kuint> &mem,
                                                     const std::string &filepath)
      {
        FileSkmerSource<kuint> src(filepath);
        if (mem.get_k() != src.k() || mem.get_m() != src.m())
          throw std::invalid_argument("In-memory list k/m mismatch with file header");

        InMemorySkmerSource<kuint> mem_src(mem.get_list(), mem.get_k(), mem.get_m());
        return internal::merge_set_op<kuint>(mem_src, src, OpType::INTERSECTION);
      }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_intersection(const std::string &filepath,
                                                     const SortedVirtualSkmerList<kuint> &mem)
      { return set_intersection(mem, filepath); }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_diff(const SortedVirtualSkmerList<kuint> &mem,
                                             const std::string &filepath)
      {
        FileSkmerSource<kuint> src(filepath);
        if (mem.get_k() != src.k() || mem.get_m() != src.m())
          throw std::invalid_argument("In-memory list k/m mismatch with file header");

        InMemorySkmerSource<kuint> mem_src(mem.get_list(), mem.get_k(), mem.get_m());
        return internal::merge_set_op<kuint>(mem_src, src, OpType::DIFF);
      }

      template <typename kuint>
      SortedVirtualSkmerList<kuint> set_diff(const std::string &filepath,
                                             const SortedVirtualSkmerList<kuint> &mem)
      {
        FileSkmerSource<kuint> src(filepath);
        if (mem.get_k() != src.k() || mem.get_m() != src.m())
          throw std::invalid_argument("In-memory list k/m mismatch with file header");

        InMemorySkmerSource<kuint> mem_src(mem.get_list(), mem.get_k(), mem.get_m());
        return internal::merge_set_op<kuint>(src, mem_src, OpType::DIFF); // file LHS, mem RHS
      }
    } // namespace setops
  } // namespace sortedlist
} // namespace km