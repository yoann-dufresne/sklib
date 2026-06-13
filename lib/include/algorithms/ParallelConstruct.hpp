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
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <kseq++/seqio.hpp>
#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/SkmerBucketWriter.hpp>
#include <algorithms/VirtualSkmer.hpp> // BucketDirEntry

namespace km
{
namespace sortedlist
{

// ---------------------------------------------------------------------------
// Phase 1 (the producer) parallel driver — sharded multi-producer.
//
// build_bucketed's phase 1 (FASTA -> rolling minimizer -> super-k-mer -> per-minimizer disk bucket)
// is otherwise single-threaded and ~73-78% of construct wall at -t8 (CONSTRUCT_SCALING_DIAG.md). Each
// super-k-mer is independent, so this fans the per-sequence work across `nw` workers: the caller
// (reader) streams whole FASTA sequences onto a bounded queue, and each worker owns its own
// SkmerManipulator and SkmerBucketWriter, routing into a private shard directory
// (tmp_dir/shard_<tid>). Phase 2 loads + concatenates a bucket's shards before sort_and_dedup, which
// makes the result independent of routing order, so the final index is byte-identical to the
// sequential build and across thread counts (sha256-verified).
//
// A super-k-mer depends only on its own sequence's nucleotides (sequence ends are handled by the
// prefix/suffix-size tracking + mask_absent_nucleotides), so partitioning sequences across workers'
// manipulators yields the same (bucket, skmer) multiset the single-manip sequential path does; only
// the within-shard order differs, which sort_and_dedup normalizes. B(2) (skip consecutive-duplicate
// skmers) stays per-worker; B(1) (the phase-2 sort+unique) drops every remaining duplicate, so the
// output set is unchanged. Caller restricts this to the fixed power-of-two bucketing; --max-ram
// (adaptive, low-memory) keeps the sequential phase 1.
//
// Fills `out_shard_dirs` (one per worker) and aggregates per-bucket `counts`. Rethrows the first
// worker/reader exception after joining.
template<typename gen>
void parallel_build_phase1(const std::string& input_path,
                           uint64_t k, uint64_t m, uint64_t n_buckets,
                           const std::function<uint64_t(gen)>& bucket_of,
                           unsigned nw, size_t writer_budget,
                           const std::filesystem::path& tmp_dir,
                           std::vector<std::filesystem::path>& out_shard_dirs,
                           std::vector<uint64_t>& counts)
{
    namespace fs = std::filesystem;
    nw = std::max(2u, nw);

    out_shard_dirs.clear();
    for (unsigned t {0}; t < nw; ++t) {
        fs::path sd = tmp_dir / ("shard_" + std::to_string(t));
        fs::create_directories(sd);
        out_shard_dirs.push_back(sd);
    }

    std::mutex mtx;
    std::condition_variable cv_item;   // a sequence (or end-of-input) is available to consume
    std::condition_variable cv_space;  // a queue slot freed
    std::deque<std::string> queue;
    bool done_reading {false};
    std::exception_ptr first_err;
    // Bound the in-flight sequence bytes so peak RAM stays O(budget) — not the whole genome — even
    // on a huge multi-chromosome input (a count-only cap could hold several GB of chromosomes at
    // once). The reader always enqueues at least one sequence when the queue is empty, so a single
    // sequence larger than the budget still makes progress (bounding peak to ~budget + one sequence).
    const size_t max_queue {std::max<size_t>(2ull * nw, 8)};
    const size_t byte_budget {std::max<size_t>(size_t{64} << 20, writer_budget)};
    size_t queued_bytes {0};

    std::vector<std::vector<uint64_t>> worker_counts(nw);
    // Spread the phase-1 write budget across the workers; keep a sane floor for tiny budgets.
    const size_t per_worker_budget {std::max<size_t>(writer_budget / nw, size_t{1} << 20)};

    auto worker = [&](unsigned tid) {
        try {
            km::SkmerManipulator<gen> manip{k, m};                 // per-worker mutable encoder state
            km::sortedlist::SkmerBucketWriter<gen> writer{out_shard_dirs[tid], n_buckets, per_worker_budget};
            std::string seq;
            for (;;) {
                {
                    std::unique_lock<std::mutex> lk(mtx);
                    cv_item.wait(lk, [&]{ return !queue.empty() || done_reading || first_err; });
                    if (first_err) return;
                    if (queue.empty()) break;                      // done_reading && drained
                    seq = std::move(queue.front());
                    queue.pop_front();
                    queued_bytes -= seq.size();
                }
                cv_space.notify_one();

                // Fresh SeqSkmerator per sequence (its buffers are value-initialized == the
                // FileSkmerator reset() the sequential path uses); the manipulator is reused across
                // this worker's sequences exactly as the sequential path reuses one across the file.
                km::SeqSkmerator<gen> rator{manip, seq};
                km::Skmer<gen> prev{};
                bool have_prev {false};
                for (const km::Skmer<gen> sk : rator) {
                    if (have_prev && sk == prev) continue;          // B(2), per-worker
                    const gen mini = manip.minimizer(sk);
                    writer.append(bucket_of(mini), sk);
                    prev = sk;
                    have_prev = true;
                }
            }
            writer.close();
            std::vector<uint64_t>& wc = worker_counts[tid];
            wc.assign(n_buckets, 0);
            for (uint64_t id {0}; id < n_buckets; ++id) wc[id] = writer.bucket_count(id);
        } catch (...) {
            std::lock_guard<std::mutex> lk(mtx);
            if (!first_err) first_err = std::current_exception();
            cv_item.notify_all();
            cv_space.notify_all();
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nw);
    for (unsigned t {0}; t < nw; ++t) workers.emplace_back(worker, t);

    // Reader (this thread): stream sequences onto the bounded queue. kseq reuses rec.seq across
    // records, so each sequence is copied into the queue.
    try {
        klibpp::SeqStreamIn ksi(input_path.c_str());
        klibpp::KSeq rec;
        while (ksi >> rec) {
            if (rec.seq.length() < k) continue;                    // matches FileSkmerator::init_record
            std::unique_lock<std::mutex> lk(mtx);
            // Enqueue when under both caps, or unconditionally when the queue is empty (guarantees
            // progress for a single sequence larger than the byte budget).
            cv_space.wait(lk, [&]{ return queue.empty() ||
                                          (queue.size() < max_queue && queued_bytes < byte_budget) ||
                                          first_err; });
            if (first_err) break;
            queued_bytes += rec.seq.size();
            queue.emplace_back(rec.seq);
            cv_item.notify_one();
        }
    } catch (...) {
        std::lock_guard<std::mutex> lk(mtx);
        if (!first_err) first_err = std::current_exception();
    }
    {
        std::lock_guard<std::mutex> lk(mtx);
        done_reading = true;
    }
    cv_item.notify_all();

    for (std::thread& w : workers) w.join();
    if (first_err) std::rethrow_exception(first_err);

    counts.assign(n_buckets, 0);
    for (unsigned t {0}; t < nw; ++t)
        if (!worker_counts[t].empty())
            for (uint64_t id {0}; id < n_buckets; ++id) counts[id] += worker_counts[t][id];
}

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
