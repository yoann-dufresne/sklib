#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <algorithms/VirtualSkmer.hpp>

// Set operations (intersection, union, difference) over two sorted super-k-mer lists.
//
// The atomic element of an operation is a *k-mer*, stored as (column, interleaved value) where the
// column is the minimizer's position inside the k-mer (0 .. k-m). Each k-mer has exactly one column,
// so the columns partition the k-mer space and every set operation decomposes per column:
//
//     op(A, B) = ⊎_bucket ⊎_column  merge2( A's k-mers at column c , B's k-mers at column c )
//
// Within a bucket the records valid at a column are sorted by SkmerManipulator::kmer_compare on that
// column — the very invariant the per-column dichotomic query (search_kmers_in_span) already relies
// on. So a two-cursor merge driven by the same kmer_compare is correct by that same invariant, and a
// single pass yields all three relations (A∩B, A\B, B\A). The variants that only need a count never
// materialize anything.
//
// Both lists must have been built with identical k, m and bucket layout (check_compatible): then
// bucket i of A and bucket i of B cover the same φ-minimizer range and store records at the same
// (quotiented) width, so records compare directly. This always holds for default-bucketed lists
// built with the same parameters.

namespace km {
namespace sortedlist {

// Cardinalities of the three disjoint relations, computed in one merge pass.
struct SetSizes {
    uint64_t inter {0};   // |A ∩ B|
    uint64_t only_a {0};  // |A \ B|
    uint64_t only_b {0};  // |B \ A|
    uint64_t uni() const { return inter + only_a + only_b; } // |A ∪ B|
};

// Advance `i` to the first index >= i whose record carries a valid k-mer at column `c` (or `n`).
template<typename store>
inline size_t set_next_valid(const km::SkmerManipulator<store>& cmp,
                             const km::Skmer<store>* L, size_t n, size_t i, uint64_t c) {
    while (i < n && !cmp.has_valid_kmer(L[i], c)) ++i;
    return i;
}

// Counting-sort a span's record indices into per-column groups: on return, column c's records are
// idx[off[c] .. off[c+1]) in ascending index order (== ascending k-mer order at column c, the sorted
// invariant). Replaces re-scanning the whole span from 0 for each column (the former set_next_valid
// hole-scan, (k-m+1)x over the span) with one O(span) bucketing pass + direct per-column walks.
template<typename store>
inline void build_column_csr(const km::SkmerManipulator<store>& cmp,
                             const km::Skmer<store>* L, size_t n, uint64_t ncols,
                             std::vector<uint32_t>& idx, std::vector<uint32_t>& off) {
    off.assign(ncols + 1, 0);
    for (size_t i {0}; i < n; ++i) {
        const auto [s, e] {cmp.get_valid_kmer_bounds(L[i])};
        if (e < s) continue;                       // no valid k-mer (unsigned: a wrapped s>e is skipped)
        for (uint64_t c {s}; c <= e && c < ncols; ++c) ++off[c + 1];
    }
    for (uint64_t c {0}; c < ncols; ++c) off[c + 1] += off[c];
    idx.resize(off[ncols]);
    std::vector<uint32_t> cur(off.begin(), off.end() - 1);
    for (size_t i {0}; i < n; ++i) {
        const auto [s, e] {cmp.get_valid_kmer_bounds(L[i])};
        if (e < s) continue;
        for (uint64_t c {s}; c <= e && c < ncols; ++c) idx[cur[c]++] = static_cast<uint32_t>(i);
    }
}

// Per-column two-cursor merge of two sorted spans. Works on any contiguous skmer range — an in-RAM
// SortedVirtualSkmerList (get_list().data()) or a lazily-loaded bucket — so the library and tests
// share one implementation. `sink` receives each emitted k-mer as (record, column) via only_a /
// only_b / both, in column-major order with ascending k-mer order within each column (identical to
// the previous set_next_valid implementation, so every sink sees the exact same call sequence).
template<typename store, typename Sink>
inline void merge_columns(const km::SkmerManipulator<store>& cmp,
                          const km::Skmer<store>* A, size_t nA,
                          const km::Skmer<store>* B, size_t nB,
                          Sink& sink) {
    const uint64_t k_minus_m {cmp.k - cmp.m};
    const uint64_t ncols {k_minus_m + 1};
    // Reused across buckets/calls (set-ops are sequential) to avoid re-allocating the index arrays.
    thread_local std::vector<uint32_t> idxA, offA, idxB, offB;
    build_column_csr<store>(cmp, A, nA, ncols, idxA, offA);
    build_column_csr<store>(cmp, B, nB, ncols, idxB, offB);
    for (uint64_t c {0}; c <= k_minus_m; ++c) {
        const uint32_t* ra {idxA.data() + offA[c]}; size_t na {offA[c + 1] - offA[c]};
        const uint32_t* rb {idxB.data() + offB[c]}; size_t nb {offB[c + 1] - offB[c]};
        size_t ia {0}, ib {0};
        while (ia < na && ib < nb) {
            const int r {cmp.kmer_compare(A[ra[ia]], B[rb[ib]], c)};
            if (r < 0)      { sink.only_a(A[ra[ia]], c); ++ia; }
            else if (r > 0) { sink.only_b(B[rb[ib]], c); ++ib; }
            else            { sink.both  (A[ra[ia]], c); ++ia; ++ib; }
        }
        for (; ia < na; ++ia) sink.only_a(A[ra[ia]], c);
        for (; ib < nb; ++ib) sink.only_b(B[rb[ib]], c);
    }
}

// Sink that only counts — backs the *_size variants (no allocation, nothing materialized).
template<typename store>
struct CountSink {
    uint64_t n_inter {0}, n_only_a {0}, n_only_b {0};
    void only_a(const km::Skmer<store>&, uint64_t) { ++n_only_a; }
    void only_b(const km::Skmer<store>&, uint64_t) { ++n_only_b; }
    void both  (const km::Skmer<store>&, uint64_t) { ++n_inter; }
};

// Sink that materializes the selected relation(s): each kept k-mer becomes a single-k-mer skmer
// (get_skmer_of_kmer) appended to `out`. The manipulator is held by non-const reference because
// get_skmer_of_kmer is non-const.
template<typename store>
struct CollectSink {
    km::SkmerManipulator<store>& manip;
    std::vector<km::Skmer<store>>& out;
    bool keep_inter, keep_only_a, keep_only_b;
    void only_a(const km::Skmer<store>& r, uint64_t c) { if (keep_only_a) out.push_back(manip.get_skmer_of_kmer(r, c)); }
    void only_b(const km::Skmer<store>& r, uint64_t c) { if (keep_only_b) out.push_back(manip.get_skmer_of_kmer(r, c)); }
    void both  (const km::Skmer<store>& r, uint64_t c) { if (keep_inter)  out.push_back(manip.get_skmer_of_kmer(r, c)); }
};

// Reject inputs the direct merge can't align. k and m must match for k-mers/columns to be comparable
// at all; n_buckets, quotient_bits and the per-bucket φ-minimizer lower bounds must match so bucket i
// of A and bucket i of B cover the same range and store records identically. All hold by construction
// for default-bucketed lists built with the same parameters.
template<typename store>
inline void check_compatible(const BucketedSkmerListReader<store>& A,
                             const BucketedSkmerListReader<store>& B) {
    if (A.k() != B.k() || A.m() != B.m())
        throw std::runtime_error(
            "set operation: both lists must share k and m (A: k=" + std::to_string(A.k()) +
            ",m=" + std::to_string(A.m()) + " ; B: k=" + std::to_string(B.k()) +
            ",m=" + std::to_string(B.m()) + ")");
    if (A.n_buckets() != B.n_buckets() || A.quotient_bits() != B.quotient_bits())
        throw std::runtime_error(
            "set operation: both lists must share the bucket layout (n_buckets/quotient_bits differ); "
            "rebuild both with identical --buckets and the same k/m");
    for (uint64_t i {0}; i < A.n_buckets(); ++i)
        if (A.bucket_lower(i) != B.bucket_lower(i))
            throw std::runtime_error(
                "set operation: bucket boundaries differ (adaptive --max-ram lists?); "
                "rebuild both with identical parameters");
}

// One merge pass over both lists, counting the three relations. Buckets are loaded on demand and
// released right after use, so peak RAM stays at ~one bucket pair rather than the whole lists.
template<typename store>
inline SetSizes set_sizes(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B) {
    check_compatible(A, B);
    const uint64_t nb {A.n_buckets()};
    km::SkmerManipulator<store> cmp{A.k(), A.m(), A.quotient_bits()};
    CountSink<store> sink;
    for (uint64_t i {0}; i < nb; ++i) {
        const std::vector<km::Skmer<store>>& sa {A.load_bucket(i)};
        const std::vector<km::Skmer<store>>& sb {B.load_bucket(i)};
        merge_columns<store>(cmp, sa.data(), sa.size(), sb.data(), sb.size(), sink);
        A.release_bucket(i);
        B.release_bucket(i);
    }
    return SetSizes{sink.n_inter, sink.n_only_a, sink.n_only_b};
}

// Materialize the selected relation(s) into a VSKMER_4 file reusing the inputs' bucket layout, and
// return the number of result k-mers (equal to the matching *_size). Per bucket: merge -> collect
// single-k-mer skmers -> re-compact into virtual super-k-mers via generate_sorted_list_from_enumeration
// (the same path construction uses; safe at b>0 since it touches only the quotient-regenerated masks)
// -> append. Records stay at the quotiented store width, so no extra truncation before writing.
// `no_compact`: skip the per-bucket super-k-mer re-compaction (the dominant cost). The merge already
// emits the kept k-mers; we sort that bucket's single-k-mer skmers by skmer order and write them as-is.
// The result stays a valid, queryable sorted list (each record holds exactly one k-mer, valid at one
// column; within a column skmer order == k-mer order, so the per-column query invariant holds), but the
// file is larger (one record per k-mer instead of ~4-5 k-mers per super-k-mer).
template<typename store>
inline uint64_t materialize_setop(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B,
                                  const std::string& out_path,
                                  bool keep_inter, bool keep_only_a, bool keep_only_b,
                                  bool no_compact = false) {
    check_compatible(A, B);
    const uint64_t k {A.k()}, m {A.m()}, b {A.quotient_bits()}, nb {A.n_buckets()};

    std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
    if (out.fail())
        throw std::runtime_error("Error opening output file for writing: " + out_path);
    VirtualSkmerSerializer<store>::write_header(out, k, m, /*count*/ 0, nb, sizeof(store), b);
    if (out.fail())
        throw std::runtime_error("Error writing header to file: " + out_path);

    // Directory carries every bucket's φ-minimizer lower bound (empty buckets included, so routing
    // covers the whole space); counts are patched in as each bucket is built.
    std::vector<BucketDirEntry> dir(nb);
    for (uint64_t i {0}; i < nb; ++i)
        dir[i] = BucketDirEntry{ A.bucket_lower(i), 0 };

    km::SkmerManipulator<store> manip{k, m, b};
    std::vector<km::Skmer<store>> collected;  // reused per-bucket scratch
    // One re-compaction worker reused across all buckets: its manipulator (and masks) are built once
    // here instead of per bucket, and its internal list/merge buffers grow to the largest bucket and
    // are then reused (generate_sorted_list_from_enumeration clears the list on entry).
    SortedVirtualSkmerList<store> sub(k, m, b);
    uint64_t total_records {0};   // virtual skmers written (what the reader offsets on)
    uint64_t total_kmers {0};     // result k-mers (== the matching _size)

    for (uint64_t i {0}; i < nb; ++i) {
        const std::vector<km::Skmer<store>>& sa {A.load_bucket(i)};
        const std::vector<km::Skmer<store>>& sb {B.load_bucket(i)};
        collected.clear();
        CollectSink<store> sink{manip, collected, keep_inter, keep_only_a, keep_only_b};
        merge_columns<store>(manip, sa.data(), sa.size(), sb.data(), sb.size(), sink);
        A.release_bucket(i);
        B.release_bucket(i);

        total_kmers += collected.size();
        if (collected.empty()) continue;

        if (no_compact) {
            // Sort the bucket's single-k-mer skmers by skmer order (Skmer::operator< on m_pair) and
            // write them directly — no super-k-mer assembly. Same k-mer set, more records.
            std::sort(collected.begin(), collected.end());
            VirtualSkmerSerializer<store>::append_payload(out, collected);
            if (out.fail())
                throw std::runtime_error("Error writing skmers to file: " + out_path);
            dir[i].count = collected.size();
            total_records += collected.size();
            continue;
        }

        sub.generate_sorted_list_from_enumeration(collected, /*greedy_chain*/ true);
        const std::vector<km::Skmer<store>>& sl {sub.get_list()};
        VirtualSkmerSerializer<store>::append_payload(out, sl);
        if (out.fail())
            throw std::runtime_error("Error writing skmers to file: " + out_path);
        dir[i].count = sub.size();
        total_records += sub.size();
    }

    VirtualSkmerSerializer<store>::patch_directory(out, dir);
    VirtualSkmerSerializer<store>::patch_count(out, total_records);
    if (out.fail())
        throw std::runtime_error("Error patching header/directory in file: " + out_path);
    out.close();
    return total_kmers;
}

// ---- public named operations (return the result k-mer count) ----
template<typename store>
inline uint64_t intersection(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, const std::string& out, bool no_compact = false) {
    return materialize_setop<store>(A, B, out, /*inter*/ true, /*only_a*/ false, /*only_b*/ false, no_compact);
}
template<typename store>
inline uint64_t set_union(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, const std::string& out, bool no_compact = false) {
    return materialize_setop<store>(A, B, out, true, true, true, no_compact);
}
// diff(A, B) = A \ B (asymmetric): keep only the k-mers present in A and absent from B.
template<typename store>
inline uint64_t difference(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, const std::string& out, bool no_compact = false) {
    return materialize_setop<store>(A, B, out, false, true, false, no_compact);
}

// ---- size-only variants (no materialization) ----
template<typename store>
inline uint64_t intersection_size(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B) {
    return set_sizes<store>(A, B).inter;
}
template<typename store>
inline uint64_t union_size(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B) {
    return set_sizes<store>(A, B).uni();
}
template<typename store>
inline uint64_t diff_size(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B) {
    return set_sizes<store>(A, B).only_a;
}

} // namespace sortedlist
} // namespace km
