#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <io/Skmer.hpp>
#include <algorithms/Parallel.hpp>
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

// Read bucket `b`'s payload from `R` into `dst` through a caller-owned file handle `fh`. This is the
// thread-safe alternative to R.load_bucket(b)/release_bucket(b) for the parallel set ops: each worker
// owns its own ifstream so positional reads on different buckets don't race, and dst is a reused
// per-worker buffer freed at the end (the shared bucket() cache, by contrast, never evicts and its
// release path is not thread-safe). Mirrors the seek+read in BucketedSkmerListReader::bucket().
template<typename store>
inline void read_bucket_into(const BucketedSkmerListReader<store>& R, std::ifstream& fh,
                             uint64_t b, std::vector<km::Skmer<store>>& dst) {
    const uint64_t cnt {R.bucket_count(b)};
    dst.resize(cnt);
    if (cnt) {
        fh.clear();                                  // drop any EOF/fail flag from a previous read
        fh.seekg(R.bucket_byte_offset(b), std::ios::beg);
        fh.read(reinterpret_cast<char*>(dst.data()),
                static_cast<std::streamsize>(cnt * sizeof(km::Skmer<store>)));
        if (fh.fail())
            throw std::runtime_error("Error reading bucket payload during set operation");
    }
}

// Serializes the per-bucket result payloads to `out` in strictly increasing bucket-id order while the
// workers finish in arbitrary order. The on-disk directory stores only per-bucket counts (the reader
// derives byte offsets by prefix sum), so payloads MUST be physically contiguous in bucket order.
// Out-of-order completions are buffered in `ready`; `cap` bounds how many completed-but-unwritten
// buckets may be held (back-pressure), keeping peak extra RAM at ~cap × one bucket payload. A worker
// whose bucket equals next_emit is always admitted, so the write frontier always advances (no
// deadlock). `dir` and `total_records` are filled here, under the lock, as buckets drain in order.
template<typename store>
struct OrderedPayloadWriter {
    std::ofstream& out;
    std::vector<BucketDirEntry>& dir;
    const std::string& out_path;
    const size_t cap;
    std::mutex mtx;
    std::condition_variable cv;
    std::map<uint64_t, std::vector<km::Skmer<store>>> ready;
    uint64_t next_emit {0};
    uint64_t total_records {0};   // virtual skmers written (what the reader offsets on)
    bool aborted {false};

    OrderedPayloadWriter(std::ofstream& o, std::vector<BucketDirEntry>& d,
                         const std::string& path, size_t capacity)
        : out(o), dir(d), out_path(path), cap(capacity) {}

    // Hand bucket `bucket_id`'s payload to the writer. Empty payloads are submitted too, so an empty
    // bucket still advances the write frontier and gets dir[id].count = 0. Throws (out of the calling
    // worker) on a write error; the caller catches, records it, and calls abort().
    void submit(uint64_t bucket_id, std::vector<km::Skmer<store>>&& payload) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]{ return aborted || ready.size() < cap || bucket_id == next_emit; });
        if (aborted) return;
        ready.emplace(bucket_id, std::move(payload));
        while (!ready.empty() && ready.begin()->first == next_emit) {
            std::vector<km::Skmer<store>>& p {ready.begin()->second};
            dir[next_emit].count = p.size();
            total_records += p.size();
            VirtualSkmerSerializer<store>::append_payload(out, p);
            if (out.fail())
                throw std::runtime_error("Error writing skmers to file: " + out_path);
            ready.erase(ready.begin());
            ++next_emit;
        }
        cv.notify_all();
    }

    // Wake every blocked worker so it returns without writing (used on the error path).
    void abort() {
        { std::lock_guard<std::mutex> lock(mtx); aborted = true; }
        cv.notify_all();
    }
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

// One merge pass over both lists, counting the three relations. Buckets are independent, so they are
// processed in parallel across `nthreads` workers (dynamic scheduling); each worker reads its buckets
// through its own file handles into reused buffers freed at the end, so peak RAM stays at ~nthreads
// bucket pairs rather than the whole lists. Per-worker counts are summed after the join — order does
// not matter for a count, so the result is identical to the sequential pass for any thread count.
template<typename store>
inline SetSizes set_sizes(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B,
                          unsigned nthreads = 8) {
    check_compatible(A, B);
    const uint64_t k {A.k()}, m {A.m()}, b {A.quotient_bits()}, nb {A.n_buckets()};
    const unsigned nthr {std::max(1u, nthreads)};

    std::vector<std::unique_ptr<km::SkmerManipulator<store>>> manips;   // owns raw masks: not movable
    std::vector<CountSink<store>> sinks(nthr);
    std::vector<std::vector<km::Skmer<store>>> bufA(nthr), bufB(nthr);
    std::vector<std::ifstream> fhA, fhB;
    manips.reserve(nthr); fhA.reserve(nthr); fhB.reserve(nthr);
    for (unsigned t {0}; t < nthr; ++t) {
        manips.push_back(std::make_unique<km::SkmerManipulator<store>>(k, m, b));
        fhA.emplace_back(A.path(), std::ios::binary);
        fhB.emplace_back(B.path(), std::ios::binary);
        if (fhA[t].fail() || fhB[t].fail())
            throw std::runtime_error("set operation: cannot open input list for parallel read");
    }

    std::atomic<bool> failed {false};
    std::exception_ptr eptr;
    std::mutex err_mtx;
    parallel_for_dynamic(nb, nthr, [&](uint64_t i, unsigned tid) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            read_bucket_into<store>(A, fhA[tid], i, bufA[tid]);
            read_bucket_into<store>(B, fhB[tid], i, bufB[tid]);
            merge_columns<store>(*manips[tid], bufA[tid].data(), bufA[tid].size(),
                                 bufB[tid].data(), bufB[tid].size(), sinks[tid]);
        } catch (...) {
            std::lock_guard<std::mutex> lk(err_mtx);
            if (!eptr) eptr = std::current_exception();
            failed.store(true, std::memory_order_relaxed);
        }
    });
    if (eptr) std::rethrow_exception(eptr);

    SetSizes total;
    for (unsigned t {0}; t < nthr; ++t) {
        total.inter  += sinks[t].n_inter;
        total.only_a += sinks[t].n_only_a;
        total.only_b += sinks[t].n_only_b;
    }
    return total;
}

// Materialize the selected relation(s) into a VSKMER_4 file reusing the inputs' bucket layout, and
// return the number of result k-mers (equal to the matching *_size). Buckets are independent and run
// in parallel across `nthreads` workers; the result is byte-identical for any thread count (each
// bucket's output is a pure function of its inputs, and the writer emits buckets in id order).
// Per bucket: merge -> collect
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
                                  bool no_compact = false, unsigned nthreads = 8) {
    check_compatible(A, B);
    const uint64_t k {A.k()}, m {A.m()}, b {A.quotient_bits()}, nb {A.n_buckets()};
    const unsigned nthr {std::max(1u, nthreads)};

    std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
    if (out.fail())
        throw std::runtime_error("Error opening output file for writing: " + out_path);
    VirtualSkmerSerializer<store>::write_header(out, k, m, /*count*/ 0, nb, sizeof(store), b);
    if (out.fail())
        throw std::runtime_error("Error writing header to file: " + out_path);

    // Directory carries every bucket's φ-minimizer lower bound (empty buckets included, so routing
    // covers the whole space); counts are patched in by the ordered writer as buckets drain in order.
    std::vector<BucketDirEntry> dir(nb);
    for (uint64_t i {0}; i < nb; ++i)
        dir[i] = BucketDirEntry{ A.bucket_lower(i), 0 };

    // Per-worker state, allocated once and reused across the buckets each worker happens to process.
    // The manipulator and the re-compaction worker own raw mask arrays and are not safely movable, so
    // they live behind unique_ptr (the vector never relocates them). The merge's CSR scratch is
    // thread_local inside merge_columns, so it is already per-worker.
    std::vector<std::unique_ptr<km::SkmerManipulator<store>>> manips;
    std::vector<std::unique_ptr<SortedVirtualSkmerList<store>>> subs;
    std::vector<std::vector<km::Skmer<store>>> collected(nthr), bufA(nthr), bufB(nthr);
    std::vector<uint64_t> kmers(nthr, 0);   // result k-mers per worker (== the matching _size, summed)
    std::vector<std::ifstream> fhA, fhB;
    manips.reserve(nthr); subs.reserve(nthr); fhA.reserve(nthr); fhB.reserve(nthr);
    for (unsigned t {0}; t < nthr; ++t) {
        manips.push_back(std::make_unique<km::SkmerManipulator<store>>(k, m, b));
        subs.push_back(std::make_unique<SortedVirtualSkmerList<store>>(k, m, b));
        fhA.emplace_back(A.path(), std::ios::binary);
        fhB.emplace_back(B.path(), std::ios::binary);
        if (fhA[t].fail() || fhB[t].fail())
            throw std::runtime_error("set operation: cannot open input list for parallel read");
    }

    // Workers finish out of order, but the count-only directory forces payloads to be written in
    // bucket order. The writer buffers out-of-order completions and drains the contiguous prefix; its
    // cap bounds the buffered payloads, so peak extra RAM is ~cap × one bucket payload. 2×nthreads is
    // ample headroom over the nthreads buckets in flight without holding many big payloads at once.
    OrderedPayloadWriter<store> writer(out, dir, out_path, static_cast<size_t>(nthr) * 2 + 1);

    std::atomic<bool> failed {false};
    std::exception_ptr eptr;
    std::mutex err_mtx;
    parallel_for_dynamic(nb, nthr, [&](uint64_t i, unsigned tid) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            read_bucket_into<store>(A, fhA[tid], i, bufA[tid]);
            read_bucket_into<store>(B, fhB[tid], i, bufB[tid]);
            std::vector<km::Skmer<store>>& col {collected[tid]};
            col.clear();
            CollectSink<store> sink{*manips[tid], col, keep_inter, keep_only_a, keep_only_b};
            merge_columns<store>(*manips[tid], bufA[tid].data(), bufA[tid].size(),
                                 bufB[tid].data(), bufB[tid].size(), sink);
            kmers[tid] += col.size();

            // The writer must own the payload (it may buffer it out of order until earlier buckets
            // drain, past this worker's next bucket), so move it out rather than copy: col and the
            // re-compaction worker's list are both refilled on the next bucket.
            std::vector<km::Skmer<store>> payload;
            if (!col.empty()) {
                if (no_compact) {
                    // Sort the single-k-mer skmers by skmer order (Skmer::operator< on m_pair) and
                    // write them as-is — no super-k-mer assembly. Same k-mer set, more records.
                    std::sort(col.begin(), col.end());
                    payload = std::move(col);       // col is re-cleared at the top of the next bucket
                } else {
                    subs[tid]->generate_sorted_list_from_enumeration(col, /*greedy_chain*/ true);
                    payload = subs[tid]->take_list();
                }
            }
            writer.submit(i, std::move(payload));   // empty payload still advances the write frontier
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (!eptr) eptr = std::current_exception();
            }
            failed.store(true, std::memory_order_relaxed);
            writer.abort();                         // unblock any worker waiting on the writer
        }
    });
    if (eptr) std::rethrow_exception(eptr);

    VirtualSkmerSerializer<store>::patch_directory(out, dir);
    VirtualSkmerSerializer<store>::patch_count(out, writer.total_records);
    if (out.fail())
        throw std::runtime_error("Error patching header/directory in file: " + out_path);
    out.close();

    uint64_t total_kmers {0};
    for (unsigned t {0}; t < nthr; ++t) total_kmers += kmers[t];
    return total_kmers;
}

// ---- public named operations (return the result k-mer count) ----
template<typename store>
inline uint64_t intersection(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, const std::string& out, bool no_compact = false, unsigned nthreads = 8) {
    return materialize_setop<store>(A, B, out, /*inter*/ true, /*only_a*/ false, /*only_b*/ false, no_compact, nthreads);
}
template<typename store>
inline uint64_t set_union(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, const std::string& out, bool no_compact = false, unsigned nthreads = 8) {
    return materialize_setop<store>(A, B, out, true, true, true, no_compact, nthreads);
}
// diff(A, B) = A \ B (asymmetric): keep only the k-mers present in A and absent from B.
template<typename store>
inline uint64_t difference(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, const std::string& out, bool no_compact = false, unsigned nthreads = 8) {
    return materialize_setop<store>(A, B, out, false, true, false, no_compact, nthreads);
}

// ---- size-only variants (no materialization) ----
template<typename store>
inline uint64_t intersection_size(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, unsigned nthreads = 8) {
    return set_sizes<store>(A, B, nthreads).inter;
}
template<typename store>
inline uint64_t union_size(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, unsigned nthreads = 8) {
    return set_sizes<store>(A, B, nthreads).uni();
}
template<typename store>
inline uint64_t diff_size(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B, unsigned nthreads = 8) {
    return set_sizes<store>(A, B, nthreads).only_a;
}

} // namespace sortedlist
} // namespace km
