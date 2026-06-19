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
#include <memory>
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
// Phase 1 (the producer) parallel driver — sharded multi-producer with intra-sequence chunking.
//
// build_bucketed's phase 1 (FASTA -> rolling minimizer -> super-k-mer -> per-minimizer disk bucket)
// is otherwise single-threaded and ~73-78% of construct wall at -t8 (CONSTRUCT_SCALING_DIAG.md). Each
// super-k-mer is independent, so this fans the work across `nw` workers, each owning its own
// SkmerManipulator and SkmerBucketWriter shard (tmp_dir/shard_<tid>). Phase 2 loads + concatenates a
// bucket's shards before sort_and_dedup, which makes the result independent of routing order, so the
// final index is byte-identical to the sequential build and across thread counts (sha256-verified).
//
// Work unit = a chunk: a core range [a,b) of a sequence (~target_chunk bases). A sequence shorter
// than a chunk is one work unit (stage 1, per-sequence); a long sequence is split into many
// (stage 2), so a single dominant chromosome no longer caps the speedup. The reader streams one
// sequence at a time (shared_ptr, so the chunks of one sequence share its single copy and peak RAM
// stays ~the largest sequence) and emits chunk items onto a bounded queue.
//
// Byte-identical across the chunk seams: a worker processes the slice [a-margin, b+margin] (clamped
// to the sequence), so every super-k-mer whose creation index (the SeqSkmerator slot-write position,
// exposed via yielded_position() when TrackPos=true) lands in [a,b) is formed with full left/right
// context — a super-k-mer spans <= 2k-m, and margin = 4*(2k-m) also covers the iterator warm-up, so
// the warm-up's partial super-k-mers all map to orig < a and are filtered out. The creation index is
// the same original-sequence position whichever chunk computes it, so the [a,b) tiling assigns each
// super-k-mer to exactly one chunk (no dup, no drop) and reproduces the sequential super-k-mer set.
// B(2) (skip consecutive-duplicate skmers) is applied per chunk on the emitted stream; B(1) (the
// phase-2 sort+unique) drops the rest, so the output set is unchanged. Caller restricts this to the
// fixed power-of-two bucketing; --max-ram (adaptive, low-memory) keeps the sequential phase 1.
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

    // margin: full context for any super-k-mer with creation index in [a,b). A super-k-mer spans
    // <= 2k-m and the iterator warm-up is ~2k-m; 4*(2k-m) is comfortably safe and the redundant
    // per-chunk work (a few hundred bases) is negligible vs a ~1 Mbp chunk core.
    const int64_t margin {static_cast<int64_t>(4 * (2 * k - m))};
    const int64_t target_chunk {int64_t{1} << 20};   // ~1 Mbp cores: good load balance, ~0.02% overhead

    struct WorkItem { std::shared_ptr<const std::string> seq; int64_t a, b; };

    std::mutex mtx;
    std::condition_variable cv_item;   // a chunk (or end-of-input) is available to consume
    std::condition_variable cv_space;  // a queue slot freed
    std::deque<WorkItem> queue;
    bool done_reading {false};
    std::exception_ptr first_err;
    // Item-count cap; combined with reading one sequence at a time this bounds peak RAM to roughly
    // the largest sequence (a big sequence's chunks all share its single shared_ptr copy).
    const size_t max_items {std::max<size_t>(4ull * nw, 16)};

    std::vector<std::vector<uint64_t>> worker_counts(nw);
    // Spread the phase-1 write budget across the workers; keep a sane floor for tiny budgets.
    const size_t per_worker_budget {std::max<size_t>(writer_budget / nw, size_t{1} << 20)};

    auto worker = [&](unsigned tid) {
        try {
            km::SkmerManipulator<gen> manip{k, m};                 // per-worker mutable encoder state
            km::sortedlist::SkmerBucketWriter<gen> writer{out_shard_dirs[tid], n_buckets, per_worker_budget};
            // One SeqSkmerator + iterator per worker, bound to a reused `slice` buffer and re-primed
            // per chunk via reset() — exactly the FileSkmerator pattern. This keeps the iterator's
            // ring buffers (and the slice string) allocated once: constructing a fresh iterator per
            // chunk made the per-chunk allocations contend on the glibc malloc arenas at high -t
            // (perf: ~7% in malloc, growing with thread count), which dragged phase-1 scaling down.
            std::string slice;
            km::SeqSkmerator<gen, true> rator{manip, slice};
            auto it = rator.begin();
            for (;;) {
                WorkItem item;
                {
                    std::unique_lock<std::mutex> lk(mtx);
                    cv_item.wait(lk, [&]{ return !queue.empty() || done_reading || first_err; });
                    if (first_err) return;
                    if (queue.empty()) break;                      // done_reading && drained
                    item = std::move(queue.front());
                    queue.pop_front();
                }
                cv_space.notify_one();

                const int64_t L {static_cast<int64_t>(item.seq->size())};
                const int64_t slice_start {std::max<int64_t>(0, item.a - margin)};
                const int64_t slice_end   {std::min<int64_t>(L, item.b + margin)};
                slice.assign(*item.seq, static_cast<size_t>(slice_start),
                             static_cast<size_t>(slice_end - slice_start));
                it.reset();   // re-prime over the new slice contents, reusing the ring buffers

                // TrackPos=true: yielded_position() gives the super-k-mer's creation index within the
                // slice; + slice_start makes it the original-sequence position used for the [a,b) tiling.
                km::Skmer<gen> prev{};
                bool have_prev {false};
                for (; !it.consumed(); ++it) {
                    const int64_t orig {slice_start + it.yielded_position()};
                    // [a,b) tiling on the creation index. The last chunk (b==L) also claims the
                    // end-of-sequence super-k-mers, whose creation index lands in the empty-nucleotide
                    // flush *past* L (so an upper bound of L would drop them — a real k-mer loss).
                    if (orig < item.a || (item.b != L && orig >= item.b)) continue; // another chunk / warm-up
                    const km::Skmer<gen> sk {*it};
                    if (have_prev && sk == prev) continue;          // B(2), on the emitted stream
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

    // Reader (this thread): one sequence at a time, split into chunk items. kseq reuses rec.seq, so
    // each sequence is copied once into a shared_ptr the chunks reference.
    try {
        klibpp::SeqStreamIn ksi(input_path.c_str());
        klibpp::KSeq rec;
        while (ksi >> rec) {
            if (rec.seq.length() < k) continue;                    // matches FileSkmerator::init_record
            auto seq = std::make_shared<const std::string>(rec.seq);
            const int64_t L {static_cast<int64_t>(seq->size())};
            for (int64_t a {0}; a < L; a += target_chunk) {
                const int64_t b {std::min<int64_t>(L, a + target_chunk)};
                std::unique_lock<std::mutex> lk(mtx);
                cv_space.wait(lk, [&]{ return queue.size() < max_items || first_err; });
                if (first_err) break;
                queue.push_back(WorkItem{seq, a, b});
                cv_item.notify_one();
            }
            if (first_err) break;
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
// `serialize(job_id, payload, dir, bytes, scratch)` turns a compacted (coarse) bucket payload into the
// on-disk byte blob for that job and patches the affected dir counts: the identity for a non-decoupled
// build (one fine bucket == the job), or the per-fine-sub-bucket split+re-truncation when the build
// decouples a coarse compaction layout from a finer query layout. `bytes`/`scratch` are per-worker.
template<typename store, typename MakeCompactor, typename Serialize>
uint64_t parallel_build_phase2(std::ofstream& out,
                               std::vector<BucketDirEntry>& dir,
                               const std::vector<uint64_t>& counts,
                               uint64_t n_buckets,
                               unsigned n_threads,
                               MakeCompactor make_compactor,
                               Serialize serialize)
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
            std::string bytes;                     // reused: the serialized blob for one job
            std::vector<Skmer<store>> scratch;     // reused: repartition's fine-payload buffer
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
                total.fetch_add(payload.size(), std::memory_order_relaxed);
                // Serialize this coarse bucket into `bytes` and patch its fine dir counts. Coarse ids
                // are disjoint per seq and map to disjoint fine-id ranges, so the dir writes race-free.
                serialize(id, payload, dir, bytes, scratch);

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
