#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
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

// --- optional env-gated phase timing for the mono-thread set-op profile -------------------------
// Attributes WALL time (including memory/I-O stalls, which perf's instruction sampling misses) to the
// four per-bucket phases: read / merge+collect / recompact / write. Enabled only when the env var
// SKLIB_UNION_PHASE_TIMING is set to a non-empty, non-"0" value AND nthreads==1 (the profiling mode);
// the gate is a single cached read so the timed path pays nothing when it is off. Meaningful only at
// -t1: the accumulators are plain doubles written from the inline (single-thread) bucket loop.
namespace detail {
inline bool setop_phase_timing_enabled() {
    static const bool on = [] {
        const char* e = std::getenv("SKLIB_UNION_PHASE_TIMING");
        return e && e[0] && !(e[0] == '0' && e[1] == '\0');
    }();
    return on;
}
struct SetopPhaseTimers {
    double read {0}, merge {0}, recompact {0}, write {0};
    using clock = std::chrono::steady_clock;
    static double secs(clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    }
    void report(const char* tag) const {
        const double tot {read + merge + recompact + write};
        const double d {tot > 0 ? tot : 1.0};
        std::cerr << "[setop-phase " << tag << "] total=" << tot << "s"
                  << "  read=" << read << "s (" << 100.0 * read / d << "%)"
                  << "  merge+collect=" << merge << "s (" << 100.0 * merge / d << "%)"
                  << "  recompact=" << recompact << "s (" << 100.0 * recompact / d << "%)"
                  << "  write=" << write << "s (" << 100.0 * write / d << "%)\n";
    }
};
} // namespace detail

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
    // Pass 1 — per-column counts via a difference array. A record valid at columns [s, hi] contributes
    // +1 to each; instead of the old per-column loop (O(hi-s+1), i.e. O(total k-mers) since A/B records
    // are full super-k-mers spanning many columns) apply it in O(1) as +1 at s and -1 at hi+1, then a
    // single prefix sum. The resulting CSR offsets `off` (and thus `idx`) are byte-identical.
    thread_local std::vector<int64_t> col_diff;
    col_diff.assign(ncols + 1, 0);
    for (size_t i {0}; i < n; ++i) {
        const auto [s, e] {cmp.get_valid_kmer_bounds(L[i])};
        if (e < s) continue;                       // no valid k-mer (unsigned: a wrapped s>e is skipped)
        const uint64_t hi {e < ncols ? e : ncols - 1};
        ++col_diff[s];
        --col_diff[hi + 1];
    }
    off.assign(ncols + 1, 0);
    int64_t running {0};
    for (uint64_t c {0}; c < ncols; ++c) { running += col_diff[c]; off[c + 1] = off[c] + static_cast<uint32_t>(running); }
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
// get_skmer_of_kmer is non-const. `col_count` (optional): per-column count of kept k-mers, filled as
// they are emitted (merge_columns visits columns in order, k-mers ascending within each), so its
// prefix sum gives the contiguous per-column blocks of `out` — handed to the re-compaction so it can
// skip re-discovering the column layout (see generate_sorted_list_from_enumeration's col_offsets).
template<typename store>
struct CollectSink {
    km::SkmerManipulator<store>& manip;
    std::vector<km::Skmer<store>>& out;
    bool keep_inter, keep_only_a, keep_only_b;
    std::vector<uint64_t>* col_count {nullptr};
    void only_a(const km::Skmer<store>& r, uint64_t c) { if (keep_only_a) { out.push_back(manip.get_skmer_of_kmer(r, c)); if (col_count) ++(*col_count)[c]; } }
    void only_b(const km::Skmer<store>& r, uint64_t c) { if (keep_only_b) { out.push_back(manip.get_skmer_of_kmer(r, c)); if (col_count) ++(*col_count)[c]; } }
    void both  (const km::Skmer<store>& r, uint64_t c) { if (keep_inter)  { out.push_back(manip.get_skmer_of_kmer(r, c)); if (col_count) ++(*col_count)[c]; } }
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
    // Per-worker scratch for the column-layout fast-path: col_count[tid] holds the per-column kept
    // count (filled by the sink), prefix-summed into col_off[tid] and handed to the re-compaction so
    // it skips re-scanning the enumeration for each column. k-m+1 columns.
    const uint64_t kmm {k - m};
    std::vector<std::vector<uint64_t>> col_count(nthr), col_off(nthr);
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
    detail::SetopPhaseTimers tmr;
    const bool timing {nthr == 1 && detail::setop_phase_timing_enabled()};
    using phase_clock = detail::SetopPhaseTimers::clock;
    parallel_for_dynamic(nb, nthr, [&](uint64_t i, unsigned tid) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            const phase_clock::time_point t0 {timing ? phase_clock::now() : phase_clock::time_point{}};
            read_bucket_into<store>(A, fhA[tid], i, bufA[tid]);
            read_bucket_into<store>(B, fhB[tid], i, bufB[tid]);
            std::vector<km::Skmer<store>>& col {collected[tid]};
            col.clear();
            std::vector<uint64_t>& cc {col_count[tid]};
            cc.assign(kmm + 1, 0);
            CollectSink<store> sink{*manips[tid], col, keep_inter, keep_only_a, keep_only_b, &cc};
            const phase_clock::time_point t1 {timing ? phase_clock::now() : t0};
            merge_columns<store>(*manips[tid], bufA[tid].data(), bufA[tid].size(),
                                 bufB[tid].data(), bufB[tid].size(), sink);
            kmers[tid] += col.size();
            const phase_clock::time_point t2 {timing ? phase_clock::now() : t0};

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
                    // Prefix-sum the per-column kept counts into contiguous block offsets so the
                    // re-compaction reads each column's ids directly instead of re-scanning col.
                    std::vector<uint64_t>& bo {col_off[tid]};
                    bo.resize(kmm + 2);
                    bo[0] = 0;
                    for (uint64_t c {0}; c <= kmm; ++c) bo[c + 1] = bo[c] + cc[c];
                    subs[tid]->generate_sorted_list_from_enumeration(col, /*greedy_chain*/ true, &bo);
                    payload = subs[tid]->take_list();
                }
            }
            const phase_clock::time_point t3 {timing ? phase_clock::now() : t0};
            writer.submit(i, std::move(payload));   // empty payload still advances the write frontier
            if (timing) {
                const phase_clock::time_point t4 {phase_clock::now()};
                tmr.read      += detail::SetopPhaseTimers::secs(t0, t1);
                tmr.merge     += detail::SetopPhaseTimers::secs(t1, t2);
                tmr.recompact += detail::SetopPhaseTimers::secs(t2, t3);
                tmr.write     += detail::SetopPhaseTimers::secs(t3, t4);
            }
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
    if (timing) tmr.report("materialize");

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

// ======================= combined (single-pass) set operations =======================
//
// One merge pass over both lists yields the three DISJOINT relations (both / only_a / only_b), and
// every set operation is a union of a subset of them (intersection={both}, union={both,only_a,only_b},
// A\B={only_a}, B\A={only_b}). So any requested subset of {A∩B, A∪B, A\B, B\A} — materialized and/or
// merely counted — can be produced in ONE pass instead of one pass per operation: the per-bucket read,
// the per-column CSR build and the merge are shared; only the per-output re-compaction is per-output.
// The four cardinalities always come for free (three counters), so the combined call also reports all
// sizes. This is the multi-operation analogue of materialize_setop, built from the same pieces
// (merge_columns, read_bucket_into, generate_sorted_list_from_enumeration, parallel_for_dynamic) and
// carrying the same byte-identical-for-any-thread-count guarantee, per output file.

// A request: each set output path is optional (unset ⇒ that relation is not materialized). With every
// path unset the call is a pure combined count (delegated to set_sizes; no file is touched).
struct MultiSetOpRequest {
    std::optional<std::string> inter_out;     // A ∩ B
    std::optional<std::string> union_out;     // A ∪ B
    std::optional<std::string> diff_ab_out;   // A \ B
    std::optional<std::string> diff_ba_out;   // B \ A
    bool no_compact {false};
    bool any_output() const { return inter_out || union_out || diff_ab_out || diff_ba_out; }
};

// All four cardinalities (always filled) plus the materialized result k-mer counts. Each *_kmers
// equals the matching cardinality (set-op results carry no duplicate k-mers); it is reported whether
// or not that relation was materialized.
struct MultiSetOpResult {
    SetSizes sizes;
    uint64_t inter_kmers {0}, union_kmers {0}, diff_ab_kmers {0}, diff_ba_kmers {0};
};

// Sink that fans each emitted k-mer out to every requested relation's buffer in one pass. A null
// target means that relation is not materialized (its count is still kept). The single-k-mer skmer is
// built once via get_skmer_of_kmer (a pure read of the manipulator's masks — verified side-effect
// free, takes its argument by value) and copied into each target, so requesting union alongside a
// partition costs one extract, not two.
template<typename store>
struct MultiCollectSink {
    km::SkmerManipulator<store>& manip;
    std::vector<km::Skmer<store>>* inter;     // nullable: A∩B
    std::vector<km::Skmer<store>>* uni;       // nullable: A∪B
    std::vector<km::Skmer<store>>* diff_ab;   // nullable: A\B
    std::vector<km::Skmer<store>>* diff_ba;   // nullable: B\A
    // Per-column kept counts, one per buffer (non-null iff the buffer is). Each buffer is filled in
    // column-major, k-mer-ascending order (merge_columns visits columns in order), so its prefix sum
    // gives the contiguous per-column blocks the re-compaction needs (the idea-#1 fast-path).
    std::vector<uint64_t>* cnt_inter {nullptr};
    std::vector<uint64_t>* cnt_uni {nullptr};
    std::vector<uint64_t>* cnt_ab {nullptr};
    std::vector<uint64_t>* cnt_ba {nullptr};
    uint64_t n_inter {0}, n_only_a {0}, n_only_b {0};

    void both(const km::Skmer<store>& r, uint64_t c) {
        ++n_inter;
        if (inter || uni) {
            const km::Skmer<store> s {manip.get_skmer_of_kmer(r, c)};
            if (inter) { inter->push_back(s); ++(*cnt_inter)[c]; }
            if (uni)   { uni->push_back(s);   ++(*cnt_uni)[c]; }
        }
    }
    void only_a(const km::Skmer<store>& r, uint64_t c) {
        ++n_only_a;
        if (diff_ab || uni) {
            const km::Skmer<store> s {manip.get_skmer_of_kmer(r, c)};
            if (diff_ab) { diff_ab->push_back(s); ++(*cnt_ab)[c]; }
            if (uni)     { uni->push_back(s);     ++(*cnt_uni)[c]; }
        }
    }
    void only_b(const km::Skmer<store>& r, uint64_t c) {
        ++n_only_b;
        if (diff_ba || uni) {
            const km::Skmer<store> s {manip.get_skmer_of_kmer(r, c)};
            if (diff_ba) { diff_ba->push_back(s); ++(*cnt_ba)[c]; }
            if (uni)     { uni->push_back(s);     ++(*cnt_uni)[c]; }
        }
    }
};

// One output file of a combined materialization. Value-semantic (pointers into the driver's stable
// `outs`/`dirs` vectors), so a std::vector of these is trivially resizable/assignable.
template<typename store>
struct MultiWriterChannel {
    std::ofstream* out {nullptr};
    std::vector<BucketDirEntry>* dir {nullptr};
    std::string path;
    uint64_t total_records {0};   // virtual skmers written to this file (what its reader offsets on)
};

// The combined analogue of OrderedPayloadWriter: emits each bucket's payloads to ALL output files in
// strictly increasing bucket-id order, under a SINGLE frontier (one mutex / cv / next_emit). A worker
// submits a *bundle* — one payload per channel for that bucket — and the writer, when the bucket is
// the next to emit, appends every channel's payload to its own file and patches every channel's
// directory entry. One shared frontier (rather than one writer per file) means one back-pressure
// dimension to reason about — exactly the proof the single-output writer already carries — and ~1/N
// the reorder-buffer RAM. Empty payloads/channels are submitted too, so every bucket advances the
// frontier and gets a count-0 directory entry in every file.
template<typename store>
struct MultiOrderedPayloadWriter {
    std::vector<MultiWriterChannel<store>>& channels;
    const size_t cap;
    std::mutex mtx;
    std::condition_variable cv;
    std::map<uint64_t, std::vector<std::vector<km::Skmer<store>>>> ready;   // bucket_id -> per-channel payloads
    uint64_t next_emit {0};
    bool aborted {false};

    MultiOrderedPayloadWriter(std::vector<MultiWriterChannel<store>>& ch, size_t capacity)
        : channels(ch), cap(capacity) {}

    void submit(uint64_t bucket_id, std::vector<std::vector<km::Skmer<store>>>&& bundle) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]{ return aborted || ready.size() < cap || bucket_id == next_emit; });
        if (aborted) return;
        ready.emplace(bucket_id, std::move(bundle));
        while (!ready.empty() && ready.begin()->first == next_emit) {
            std::vector<std::vector<km::Skmer<store>>>& b {ready.begin()->second};
            for (size_t ci {0}; ci < channels.size(); ++ci) {
                MultiWriterChannel<store>& ch {channels[ci]};
                std::vector<km::Skmer<store>>& p {b[ci]};
                (*ch.dir)[next_emit].count = p.size();
                ch.total_records += p.size();
                VirtualSkmerSerializer<store>::append_payload(*ch.out, p);
                if (ch.out->fail())
                    throw std::runtime_error("Error writing skmers to file: " + ch.path);
            }
            ready.erase(ready.begin());
            ++next_emit;
        }
        cv.notify_all();
    }

    void abort() {
        { std::lock_guard<std::mutex> lock(mtx); aborted = true; }
        cv.notify_all();
    }
};

// Compute any requested subset of {A∩B, A∪B, A\B, B\A} in a single pass and return all four sizes.
// With no output requested this is a pure count (set_sizes, nothing materialized). Otherwise each
// requested relation is written to its own VSKMER_4 file reusing the inputs' bucket layout; per
// output the file is byte-identical for any thread count. Buckets run in parallel across `nthreads`
// workers; each worker does ONE merge per bucket, fans the kept k-mers into per-output buffers, then
// re-compacts each (or, with no_compact, sorts each) reusing one SortedVirtualSkmerList sequentially,
// and submits the per-bucket bundle to the shared ordered writer.
template<typename store>
inline MultiSetOpResult multi_setop(BucketedSkmerListReader<store>& A, BucketedSkmerListReader<store>& B,
                                    const MultiSetOpRequest& req, unsigned nthreads = 8) {
    check_compatible(A, B);

    MultiSetOpResult result;

    // No materialization requested: a pure combined count. set_sizes is already the optimal
    // single-pass counter; reuse it and touch no files.
    if (!req.any_output()) {
        result.sizes = set_sizes<store>(A, B, nthreads);
        result.inter_kmers   = result.sizes.inter;
        result.union_kmers   = result.sizes.uni();
        result.diff_ab_kmers = result.sizes.only_a;
        result.diff_ba_kmers = result.sizes.only_b;
        return result;
    }

    const uint64_t k {A.k()}, m {A.m()}, b {A.quotient_bits()}, nb {A.n_buckets()};
    const unsigned nthr {std::max(1u, nthreads)};

    // Stable channel list (fixes bundle-slot ↔ file mapping). Per-relation channel index, -1 if that
    // relation is not requested; the sink routes by these.
    enum Rel { R_INTER, R_UNION, R_DIFF_AB, R_DIFF_BA };
    struct OutSpec { Rel rel; const std::string* path; };
    std::vector<OutSpec> specs;
    if (req.inter_out)   specs.push_back({R_INTER,   &*req.inter_out});
    if (req.union_out)   specs.push_back({R_UNION,   &*req.union_out});
    if (req.diff_ab_out) specs.push_back({R_DIFF_AB, &*req.diff_ab_out});
    if (req.diff_ba_out) specs.push_back({R_DIFF_BA, &*req.diff_ba_out});
    const size_t nout {specs.size()};

    // Reject two outputs writing to the same path (they would interleave and corrupt each other).
    for (size_t i {0}; i < nout; ++i)
        for (size_t j {i + 1}; j < nout; ++j)
            if (*specs[i].path == *specs[j].path)
                throw std::runtime_error("set operation: two results target the same output path: " + *specs[i].path);

    int ch_inter {-1}, ch_union {-1}, ch_ab {-1}, ch_ba {-1};
    for (size_t ci {0}; ci < nout; ++ci) {
        switch (specs[ci].rel) {
            case R_INTER:   ch_inter = static_cast<int>(ci); break;
            case R_UNION:   ch_union = static_cast<int>(ci); break;
            case R_DIFF_AB: ch_ab    = static_cast<int>(ci); break;
            case R_DIFF_BA: ch_ba    = static_cast<int>(ci); break;
        }
    }

    // Open one output stream per requested relation, write its header and pre-fill its directory with
    // every bucket's φ-minimizer lower bound (counts patched in as buckets drain in order). `outs` and
    // `dirs` are sized up front and never reallocated, so the pointers held by the channels stay valid.
    std::vector<std::ofstream> outs(nout);
    std::vector<std::vector<BucketDirEntry>> dirs(nout, std::vector<BucketDirEntry>(nb));
    for (size_t ci {0}; ci < nout; ++ci) {
        outs[ci].open(*specs[ci].path, std::ios::binary | std::ios::trunc);
        if (outs[ci].fail())
            throw std::runtime_error("Error opening output file for writing: " + *specs[ci].path);
        VirtualSkmerSerializer<store>::write_header(outs[ci], k, m, /*count*/ 0, nb, sizeof(store), b);
        if (outs[ci].fail())
            throw std::runtime_error("Error writing header to file: " + *specs[ci].path);
        for (uint64_t i {0}; i < nb; ++i)
            dirs[ci][i] = BucketDirEntry{ A.bucket_lower(i), 0 };
    }

    std::vector<MultiWriterChannel<store>> channels(nout);
    for (size_t ci {0}; ci < nout; ++ci)
        channels[ci] = MultiWriterChannel<store>{ &outs[ci], &dirs[ci], *specs[ci].path, 0 };
    MultiOrderedPayloadWriter<store> writer(channels, static_cast<size_t>(nthr) * 2 + 1);

    // Per-worker state, allocated once and reused across the buckets each worker happens to process.
    // Re-compaction uses ONE SortedVirtualSkmerList *per output channel* (not one shared across the
    // up-to-4 outputs): generate_sorted_list_from_enumeration reuses internal scratch, and recompacting
    // a different relation's enumeration in between perturbs the next result for wide/many-column
    // records. A dedicated recompactor per relation feeds it exactly the per-relation, bucket-ordered
    // sequence the single-op path feeds its recompactor, so each output is byte-identical to the
    // corresponding single-op materialization. The few extra recompactors per worker are lightweight.
    std::vector<std::unique_ptr<km::SkmerManipulator<store>>> manips;
    std::vector<std::vector<std::unique_ptr<SortedVirtualSkmerList<store>>>> subs(nthr);
    std::vector<std::vector<km::Skmer<store>>> bufA(nthr), bufB(nthr);
    // One collected buffer per output channel, per worker (reused across buckets).
    std::vector<std::vector<std::vector<km::Skmer<store>>>> collected(
        nthr, std::vector<std::vector<km::Skmer<store>>>(nout));
    // Per-channel per-column kept counts (cc) and their prefix-sum offsets (co), per worker — feed the
    // idea-#1 column-offset fast-path to each output's recompaction so it skips the per-column rescan.
    const uint64_t kmm {k - m};
    std::vector<std::vector<std::vector<uint64_t>>> collected_cnt(nthr, std::vector<std::vector<uint64_t>>(nout));
    std::vector<std::vector<std::vector<uint64_t>>> collected_off(nthr, std::vector<std::vector<uint64_t>>(nout));
    std::vector<uint64_t> acc_inter(nthr, 0), acc_a(nthr, 0), acc_b(nthr, 0);  // per-worker counts
    std::vector<std::ifstream> fhA, fhB;
    manips.reserve(nthr); fhA.reserve(nthr); fhB.reserve(nthr);
    for (unsigned t {0}; t < nthr; ++t) {
        manips.push_back(std::make_unique<km::SkmerManipulator<store>>(k, m, b));
        if (!req.no_compact) {                              // recompactors only needed when compacting
            subs[t].reserve(nout);
            for (size_t ci {0}; ci < nout; ++ci)
                subs[t].push_back(std::make_unique<SortedVirtualSkmerList<store>>(k, m, b));
        }
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

            std::vector<std::vector<km::Skmer<store>>>& col {collected[tid]};
            std::vector<std::vector<uint64_t>>& cc {collected_cnt[tid]};
            for (size_t ci {0}; ci < nout; ++ci) { col[ci].clear(); cc[ci].assign(kmm + 1, 0); }
            MultiCollectSink<store> sink{
                *manips[tid],
                ch_inter >= 0 ? &col[ch_inter] : nullptr,
                ch_union >= 0 ? &col[ch_union] : nullptr,
                ch_ab    >= 0 ? &col[ch_ab]    : nullptr,
                ch_ba    >= 0 ? &col[ch_ba]    : nullptr,
                ch_inter >= 0 ? &cc[ch_inter] : nullptr,
                ch_union >= 0 ? &cc[ch_union] : nullptr,
                ch_ab    >= 0 ? &cc[ch_ab]    : nullptr,
                ch_ba    >= 0 ? &cc[ch_ba]    : nullptr};
            merge_columns<store>(*manips[tid], bufA[tid].data(), bufA[tid].size(),
                                 bufB[tid].data(), bufB[tid].size(), sink);
            acc_inter[tid] += sink.n_inter;
            acc_a[tid]     += sink.n_only_a;
            acc_b[tid]     += sink.n_only_b;

            // Build this bucket's bundle: one payload per channel. The writer must own each payload
            // (it may buffer it out of order), so move it out — the collected buffers and the
            // re-compaction worker's list are both refilled on the next bucket.
            std::vector<std::vector<km::Skmer<store>>> bundle(nout);
            for (size_t ci {0}; ci < nout; ++ci) {
                std::vector<km::Skmer<store>>& cv {col[ci]};
                if (cv.empty()) continue;                     // empty payload still advances the frontier
                if (req.no_compact) {
                    std::sort(cv.begin(), cv.end());
                    bundle[ci] = std::move(cv);               // cv re-cleared at the top of the next bucket
                } else {
                    // Prefix-sum this channel's per-column counts into block offsets (idea #1).
                    std::vector<uint64_t>& bo {collected_off[tid][ci]};
                    bo.resize(kmm + 2);
                    bo[0] = 0;
                    for (uint64_t c {0}; c <= kmm; ++c) bo[c + 1] = bo[c] + cc[ci][c];
                    subs[tid][ci]->generate_sorted_list_from_enumeration(cv, /*greedy_chain*/ true, &bo);
                    bundle[ci] = subs[tid][ci]->take_list();  // dedicated per-output recompactor (see above)
                }
            }
            writer.submit(i, std::move(bundle));
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (!eptr) eptr = std::current_exception();
            }
            failed.store(true, std::memory_order_relaxed);
            writer.abort();
        }
    });
    if (eptr) std::rethrow_exception(eptr);

    for (size_t ci {0}; ci < nout; ++ci) {
        VirtualSkmerSerializer<store>::patch_directory(outs[ci], dirs[ci]);
        VirtualSkmerSerializer<store>::patch_count(outs[ci], channels[ci].total_records);
        if (outs[ci].fail())
            throw std::runtime_error("Error patching header/directory in file: " + *specs[ci].path);
        outs[ci].close();
    }

    SetSizes total;
    for (unsigned t {0}; t < nthr; ++t) {
        total.inter  += acc_inter[t];
        total.only_a += acc_a[t];
        total.only_b += acc_b[t];
    }
    result.sizes = total;
    result.inter_kmers   = total.inter;
    result.union_kmers   = total.uni();
    result.diff_ab_kmers = total.only_a;
    result.diff_ba_kmers = total.only_b;
    return result;
}

} // namespace sortedlist
} // namespace km
