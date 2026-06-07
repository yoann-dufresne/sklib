#ifndef SKLIB_PARALLEL_CONSTRUCT_HPP
#define SKLIB_PARALLEL_CONSTRUCT_HPP

// Parallel driver for phase 2 of the bucketed construction (see SortedSkmerListBuilder.hpp): the
// per-bucket compaction. Phase 2 loads each minimizer bucket, sorts+dedups it, runs the column
// algorithm (`generate_sorted_list_from_enumeration` — the construction hotspot, ~65-70% of build
// time) and down-converts each record to the storage width. Every bucket is fully independent, so
// the compaction runs on a pool of worker threads while the resulting payloads are written to the
// output file strictly in increasing bucket-id order. The on-disk bytes are therefore **identical**
// to the sequential path: only *when* a bucket is compacted changes — never *what* it contains nor
// *where* it lands in the file.
//
// Concurrency model (std-library only, no external deps):
//   - workers dynamically claim the next non-empty bucket via one atomic counter (so the ~0.2-0.3%-
//     of-work skew between buckets is load-balanced automatically — no static partition);
//   - each worker owns one reusable compactor (its SortedVirtualSkmerList scratch is reused across
//     the buckets that worker processes — the masks/manipulator are built once);
//   - a bounded reorder buffer holds at most `max_inflight` completed-but-unwritten bucket payloads;
//     a worker blocks before *starting* a bucket whose write-order index runs too far ahead of the
//     writer, which caps peak RAM to O(n_threads) bucket payloads + O(n_threads) compaction scratch;
//   - the calling (main) thread is the single writer: it emits payloads in id order, fills the
//     directory counts and accumulates the total.
//
// Deadlock-freedom: write-order indices (`seq`) are handed out 0,1,2,… in increasing bucket-id
// order; the worker holding the current `next_emit` bucket has `seq == next_emit`, so its throttle
// predicate `seq < next_emit + max_inflight` is always true — it is never blocked, the writer always
// gets the payload it is waiting for, advances, and frees a slot for the throttled workers.

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <io/Skmer.hpp>
#include <algorithms/VirtualSkmer.hpp> // BucketDirEntry

namespace km
{
namespace sortedlist
{

// Compact every non-empty bucket of `counts` in parallel and append the results to `out` in
// increasing bucket-id order (byte-identical to the sequential phase-2 loop).
//
//   make_compactor()  -> a fresh, thread-owned compactor; one is created per worker thread.
//   compactor.compact(id) -> std::vector<Skmer<store>>: loads bucket `id`, removes its temp file,
//       sorts+dedups, runs the column algorithm and truncates to the storage width. Returns an
//       empty vector for an empty bucket.
//
// Sets `dir[id].count` for every non-empty bucket (untouched for empty ones — distinct id per
// worker, so the writes are race-free) and returns the total record count. Rethrows the first
// worker exception (if any) after joining. The caller checks `out`'s error state.
template<typename store, typename MakeCompactor>
uint64_t parallel_build_phase2(std::ofstream& out,
                               std::vector<BucketDirEntry>& dir,
                               const std::vector<uint64_t>& counts,
                               uint64_t n_buckets,
                               unsigned n_threads,
                               MakeCompactor make_compactor)
{
    // Non-empty bucket ids in increasing order; the index into this list is the write-order `seq`.
    std::vector<uint64_t> jobs;
    jobs.reserve(n_buckets);
    for (uint64_t id = 0; id < n_buckets; ++id)
        if (counts[id] > 0) jobs.push_back(id);
    const uint64_t n_jobs = jobs.size();

    const unsigned nw = std::max(1u, n_threads);
    // Reorder-buffer depth: how many completed-but-unwritten buckets may pile up ahead of the
    // writer. The writer only needs ~nw slots to keep every worker fed; 2·nw gives a little slack
    // while keeping the buffered payloads (and thus peak RAM) bounded and tight even at high -t.
    const uint64_t max_inflight = std::max<uint64_t>(2ull * nw, 16);

    std::atomic<uint64_t> next_claim{0};
    std::atomic<uint64_t> total{0};

    std::mutex mtx;
    std::condition_variable cv_ready; // the payload for the current `next_emit` became available
    std::condition_variable cv_space; // `next_emit` advanced -> throttled workers may proceed
    std::unordered_map<uint64_t, std::string> ready; // seq -> serialized payload bytes (out-of-order)
    uint64_t next_emit = 0;
    std::exception_ptr first_err;

    auto worker = [&]() {
        try {
            auto comp = make_compactor();
            for (;;) {
                const uint64_t seq = next_claim.fetch_add(1, std::memory_order_relaxed);
                if (seq >= n_jobs) break;

                // Backpressure: stay within `max_inflight` of the writer. The holder of
                // seq == next_emit is never throttled, guaranteeing forward progress.
                {
                    std::unique_lock<std::mutex> lk(mtx);
                    cv_space.wait(lk, [&]{ return first_err || seq < next_emit + max_inflight; });
                    if (first_err) return;
                }

                const uint64_t id = jobs[seq];
                // `compact` returns a reference into the worker's own reused buffer, valid until the
                // next compact() call; we copy the bytes into `bytes` right here, before reusing it.
                const std::vector<Skmer<store>>& payload = comp.compact(id);
                dir[id].count = payload.size();                       // unique id per seq => race-free
                total.fetch_add(payload.size(), std::memory_order_relaxed);

                std::string bytes;
                if (!payload.empty())
                    bytes.assign(reinterpret_cast<const char*>(payload.data()),
                                 payload.size() * sizeof(Skmer<store>));

                std::lock_guard<std::mutex> lk(mtx);
                ready.emplace(seq, std::move(bytes));
                if (seq == next_emit) cv_ready.notify_one();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(mtx);
            if (!first_err) first_err = std::current_exception();
            cv_ready.notify_all();
            cv_space.notify_all();
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nw);
    for (unsigned t = 0; t < nw; ++t) workers.emplace_back(worker);

    // Single writer (this thread): emit payloads in seq (= bucket-id) order, holding the file I/O
    // outside the lock so workers never block on disk writes.
    {
        std::unique_lock<std::mutex> lk(mtx);
        std::vector<std::string> batch;
        while (next_emit < n_jobs && !first_err) {
            cv_ready.wait(lk, [&]{ return first_err || ready.find(next_emit) != ready.end(); });
            if (first_err) break;
            for (auto it = ready.find(next_emit); it != ready.end(); it = ready.find(next_emit)) {
                batch.push_back(std::move(it->second));
                ready.erase(it);
                ++next_emit;
            }
            lk.unlock();
            for (const std::string& s : batch)
                if (!s.empty()) out.write(s.data(), static_cast<std::streamsize>(s.size()));
            batch.clear();
            cv_space.notify_all(); // the writer advanced -> wake throttled workers
            lk.lock();
        }
    }

    for (std::thread& w : workers) w.join();
    if (first_err) std::rethrow_exception(first_err);
    return total.load();
}

} // namespace sortedlist
} // namespace km

#endif // SKLIB_PARALLEL_CONSTRUCT_HPP
