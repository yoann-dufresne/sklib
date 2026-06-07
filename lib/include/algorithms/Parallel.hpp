#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <thread>
#include <vector>

// Minimal dynamic-scheduled parallel-for, used by the set operations to process independent
// minimizer buckets concurrently. No thread pool, no third-party dependency — just std::thread
// over an atomic work counter, matching the rest of the library (ParallelQuery.hpp). Buckets are
// heavily skewed in size, so a static id-range split would stall on the big buckets; an atomic
// fetch_add hands each idle worker the next bucket, which keeps every thread busy.

namespace km {
namespace sortedlist {

// Run fn(i, tid) for every i in [0, n) across at most `nthreads` workers, scheduling indices
// dynamically (each worker grabs the next index with one atomic fetch_add). `tid` is the worker
// index in [0, nthreads) — use it to pick per-worker scratch so no two concurrent calls share
// mutable state. The calling thread participates as worker 0, so only nthreads-1 threads are
// spawned. With nthreads <= 1 (or n <= 1) the loop runs inline on the caller (no thread spawn),
// giving byte-for-byte the sequential behaviour.
//
// The helper is intentionally exception-naive: an exception escaping `fn` on a spawned thread
// would call std::terminate. Callers that can fail (I/O) must catch inside `fn` and surface the
// error themselves (the set-op workers do this via a shared exception_ptr).
template<typename Fn>
inline void parallel_for_dynamic(uint64_t n, unsigned nthreads, Fn&& fn) {
    if (nthreads <= 1 || n <= 1) {
        for (uint64_t i {0}; i < n; ++i) fn(i, 0u);
        return;
    }
    nthreads = static_cast<unsigned>(std::min<uint64_t>(nthreads, n));
    std::atomic<uint64_t> next {0};
    auto loop = [&](unsigned tid) {
        uint64_t i;
        while ((i = next.fetch_add(1, std::memory_order_relaxed)) < n) fn(i, tid);
    };
    std::vector<std::thread> workers;
    workers.reserve(nthreads - 1);
    for (unsigned t {1}; t < nthreads; ++t) workers.emplace_back(loop, t);
    loop(0u);                       // the calling thread is worker 0
    for (std::thread& w : workers) w.join();
}

} // namespace sortedlist
} // namespace km
