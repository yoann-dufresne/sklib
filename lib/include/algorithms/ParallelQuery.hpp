#ifndef SKLIB_PARALLEL_QUERY_HPP
#define SKLIB_PARALLEL_QUERY_HPP

// Multithreaded file-query for a bucketed sorted skmer list. A single producer thread reads the
// input and chunks the canonical super-k-mers into fixed-size, input-order batches; (n_threads - 1)
// consumer threads query those batches concurrently against a shared BucketedSkmerListReader and
// format each batch's results into one string; results are written to `os` in input order on the
// calling thread. Output is byte-identical to the sequential BucketedSkmerListReader::query().
//
// Concurrency model (std-library only, no external deps):
//   producer ──push──▶ bounded WorkQueue (FIFO, input-order batches) ──pop──▶ N consumers
//   consumers ─query_into(reused buf)+format(one string/batch)─▶ put(seq) ─▶ OrderedSink ──▶ os
//
// Batches are contiguous input ranges tagged with a sequence number, so the sink reorders by batch
// (a few thousand entries) instead of per query, and each consumer reuses one result buffer and one
// output string per batch — keeping per-query heap allocations out of the hot path. The reader is
// shared: query_into()/bucket_of_phi_min() are thread-safe (each bucket is loaded once under a lock,
// hits are lock-free), so consumers never re-read a bucket. The producer parses at the full `gen`
// width and down-converts each query to the stored `store` width before handing it off.

#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <deque>
#include <memory>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <ostream>
#include <utility>

#include <kseq++/seqio.hpp>
#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

#ifndef SKLIB_QUERY_BUCKET_SORT
// Consumer-side intra-window bucket-locality experiment (A/B). 0 = current behaviour: query each
// batch in input order (baseline; the #else path below is byte-for-byte the shipped loop). 1 =
// counting-sort each batch's items by bucket_id and query in bucket-grouped order into a reused
// per-slot result store, then format results back in ORIGINAL input order. Output is byte-identical
// either way (query_into is a pure function of (bucket_id, query); output order is reconstructed from
// the original position, never from execution order). Default 0 (off).
#define SKLIB_QUERY_BUCKET_SORT 0
#endif

namespace km
{
namespace sortedlist
{

namespace parallel_detail
{

// A unit of work: a contiguous run of input super-k-mers, tagged with its batch sequence number so
// the sink can re-emit batches in input order regardless of which consumer finishes first. Each
// item is (bucket id, down-converted query): the producer routes + truncates up front (it holds the
// wide `gen` skmer with the full minimizer), so consumers only search the stored `store` records.
template<typename store>
struct WorkBatch {
    uint64_t seq;
    std::vector<std::pair<uint64_t, Skmer<store>>> items;
};

// Append one presence vector to `out` exactly as the sequential print_query_results would: comma-
// separated booleans followed by '\n', or nothing for an empty result (a skmer with no valid k-mer
// position). Appending into a reused string keeps formatting allocation-free on the hot path.
inline void append_result(std::string& out, const std::vector<uint8_t>& result) {
    if (result.empty()) return;
    out += result[0] ? '1' : '0';
    for (size_t i {1}; i < result.size(); ++i) {
        out += ',';
        out += result[i] ? '1' : '0';
    }
    out += '\n';
}

// Bounded FIFO hand-off between the single producer and the consumers. push() blocks while full
// (backpressure caps in-flight RAM); pop() blocks until a batch is available and returns false once
// the producer is done and the queue is drained.
template<typename store>
class WorkQueue {
public:
    explicit WorkQueue(size_t capacity) : m_capacity(capacity) {}

    void push(WorkBatch<store>&& batch) {
        std::unique_lock<std::mutex> lock(m_mtx);
        m_not_full.wait(lock, [&]{ return m_queue.size() < m_capacity; });
        m_queue.push(std::move(batch));
        lock.unlock();
        m_not_empty.notify_one();
    }

    bool pop(WorkBatch<store>& out) {
        std::unique_lock<std::mutex> lock(m_mtx);
        m_not_empty.wait(lock, [&]{ return !m_queue.empty() || m_done; });
        if (m_queue.empty())
            return false; // producer finished and the queue is drained
        out = std::move(m_queue.front());
        m_queue.pop();
        lock.unlock();
        m_not_full.notify_one();
        return true;
    }

    void set_done() {
        {
            std::lock_guard<std::mutex> lock(m_mtx);
            m_done = true;
        }
        m_not_empty.notify_all();
    }

private:
    std::queue<WorkBatch<store>> m_queue;
    std::mutex m_mtx;
    std::condition_variable m_not_empty;
    std::condition_variable m_not_full;
    size_t m_capacity;
    bool m_done {false};
};

// Collects pre-formatted per-batch strings and writes them to `os` in batch (= input) order on the
// calling (main) thread. The mutex guards only the in-order bookkeeping; the actual file writes
// happen outside the lock, so consumers never block on output I/O.
class OrderedSink {
public:
    explicit OrderedSink(std::ostream& os) : m_os(os) {}

    // Called by consumers, in any order, exactly once per batch sequence number in [0, total).
    void put(uint64_t seq, std::string&& batch_text) {
        std::unique_lock<std::mutex> lock(m_mtx);
        m_ready.emplace(seq, std::move(batch_text));
        const bool wakes_writer = (seq == m_next_emit);
        lock.unlock();
        if (wakes_writer)
            m_cv.notify_one();
    }

    // Called by the producer once the total number of batches is known (after EOF).
    void set_total(uint64_t total) {
        {
            std::lock_guard<std::mutex> lock(m_mtx);
            m_total = total;
        }
        m_cv.notify_one();
    }

    // Drains batches in sequence order until every batch has been emitted. Runs on the main thread.
    void run() {
        std::vector<std::string> drained;
        std::unique_lock<std::mutex> lock(m_mtx);
        while (m_next_emit < m_total) {
            m_cv.wait(lock, [&]{
                return m_next_emit >= m_total || m_ready.find(m_next_emit) != m_ready.end();
            });
            // Detach the longest contiguous run of ready batches, then write it without the lock.
            for (auto it = m_ready.find(m_next_emit); it != m_ready.end(); it = m_ready.find(m_next_emit)) {
                drained.push_back(std::move(it->second));
                m_ready.erase(it);
                ++m_next_emit;
            }
            lock.unlock();
            for (const std::string& text : drained)
                if (!text.empty())
                    m_os.write(text.data(), static_cast<std::streamsize>(text.size()));
            drained.clear();
            lock.lock();
        }
    }

private:
    std::ostream& m_os;
    std::unordered_map<uint64_t, std::string> m_ready; // completed-but-not-yet-emitted batch texts
    uint64_t m_next_emit {0};
    uint64_t m_total {UINT64_MAX}; // unknown until the producer finishes reading the input
    std::mutex m_mtx;
    std::condition_variable m_cv;
};

} // namespace parallel_detail

// Query every super-k-mer of `filename` against `reader` using `n_threads` worker threads, writing
// presence vectors to `os` in input order. Equivalent in output to reader.query(filename, os).
//
// The parse is parallel too. It used to run on ONE producer thread that read the FASTA, drove the
// Skmerator, routed and truncated, while the other threads only searched. Once the search itself got
// cheap (per-bucket minimizer-prefix narrowing) that producer became the ceiling: it accounts for
// ~31% of the CPU work at k=31, which caps the whole query at ~3x however many threads are given
// (measured: -t4 2.7x, and -t6 SLOWER than -t4).
//
// Now a reader thread only slices the input into work items and every worker does parse + query for
// its own item. Two item shapes, both of which reproduce the sequential super-k-mer stream exactly:
//
//   * a run of WHOLE records (the common case — reads are far shorter than the chunk target). Each
//     record is enumerated in full, exactly as FileSkmerator does, so there is no seam to reason
//     about at all.
//   * one [a,b) sub-range of a record longer than the chunk target. The worker enumerates
//     [a-margin, b+margin] and keeps the super-k-mers whose CREATION index falls in [a,b) — the same
//     tiling ParallelConstruct uses, with the same margin of 4*(2k-m) (a super-k-mer spans <= 2k-m
//     and the iterator warm-up is ~2k-m). The creation index is the same original-sequence position
//     whichever chunk computes it, so each super-k-mer is claimed by exactly one chunk. The last
//     chunk (b == L) also claims the end-of-sequence flush, whose creation index runs past L.
//
// The iterator yields by increasing creation index and carries one creation index across the split
// pieces of an ambiguous super-k-mer, so within an item the order is the sequential order and the
// items partition the stream in input order. Emitting item texts by item index therefore reproduces
// the sequential output byte for byte — which is what the -t identity gate checks.
template<typename gen, typename store = gen>
void parallel_query(BucketedSkmerListReader<store>& reader, const std::string& filename,
                    std::ostream& os, unsigned n_threads, uint64_t batch_size = 4096) {
    using namespace parallel_detail;
    (void)batch_size;   // kept for API compatibility; the unit of work is now an input slice

    const unsigned nw = std::max(1u, n_threads);
    const uint64_t k = reader.k(), m = reader.m(), b = reader.quotient_bits();

    // Target bases per work item. Big enough that the per-item overhead (one queue hand-off, one
    // sink entry, one output string) is negligible, small enough to keep every worker fed and to
    // bound in-flight RAM.
    int64_t target_chunk {int64_t{1} << 20};
    if (const char* w = std::getenv("SKLIB_QUERY_CHUNK_BP")) {
        const long long v = std::atoll(w);
        if (v > 0) target_chunk = static_cast<int64_t>(v);
    }
    const int64_t margin {static_cast<int64_t>(4 * (2 * k - m))};

    struct WorkItem {
        std::vector<std::shared_ptr<const std::string>> whole;  // enumerated in full, in order
        std::shared_ptr<const std::string> big;                 // non-null => sub-range [a,b) of *big
        int64_t a {0}, b {0};
        uint64_t idx {0};
    };

    std::mutex mtx;
    std::condition_variable cv_item, cv_space;
    std::deque<WorkItem> queue;
    bool done_reading {false};
    std::exception_ptr first_err;
    const size_t max_items {std::max<size_t>(4ull * nw, 16)};

    OrderedSink sink(os);

    auto worker = [&]() {
        try {
            km::SkmerManipulator<gen> manip{k, m};
            std::string buf_whole, buf_big;
            km::SeqSkmerator<gen, false> rator_whole{manip, buf_whole};
            km::SeqSkmerator<gen, true>  rator_big{manip, buf_big};
            auto it_whole = rator_whole.begin();
            auto it_big = rator_big.begin();
            std::vector<uint8_t> res;     // reused presence buffer
            std::string text;             // reused per-item output

            for (;;) {
                WorkItem item;
                {
                    std::unique_lock<std::mutex> lk(mtx);
                    cv_item.wait(lk, [&]{ return !queue.empty() || done_reading || first_err; });
                    if (first_err) return;
                    if (queue.empty()) break;
                    item = std::move(queue.front());
                    queue.pop_front();
                }
                cv_space.notify_one();

                text.clear();
                auto run_one = [&](const km::Skmer<gen>& sk) {
                    const uint64_t bid {reader.route_minimizer(manip.minimizer(sk))};
                    const km::Skmer<store> trunc {km::truncate_skmer<gen, store>(k, m, b, sk)};
                    reader.query_into(bid, trunc, res);
                    append_result(text, res);
                };

                for (const std::shared_ptr<const std::string>& rec : item.whole) {
                    if (rec->length() < k) continue;      // matches FileSkmerator::init_record
                    buf_whole.assign(*rec);
                    it_whole.reset();
                    for (; !it_whole.consumed(); ++it_whole) run_one(*it_whole);
                }
                if (item.big) {
                    const int64_t L {static_cast<int64_t>(item.big->size())};
                    const int64_t s0 {std::max<int64_t>(0, item.a - margin)};
                    const int64_t s1 {std::min<int64_t>(L, item.b + margin)};
                    buf_big.assign(*item.big, static_cast<size_t>(s0), static_cast<size_t>(s1 - s0));
                    it_big.reset();
                    for (; !it_big.consumed(); ++it_big) {
                        const int64_t orig {s0 + it_big.yielded_position()};
                        if (orig < item.a || (item.b != L && orig >= item.b)) continue;
                        run_one(*it_big);
                    }
                }
                sink.put(item.idx, std::move(text));
                text = std::string();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(mtx);
            if (!first_err) first_err = std::current_exception();
            cv_item.notify_all();
            cv_space.notify_all();
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nw);
    for (unsigned t {0}; t < nw; ++t) workers.emplace_back(worker);

    // Reader (this thread): slice the input into items, in input order.
    uint64_t idx {0};
    auto emit = [&](WorkItem&& wi) {
        wi.idx = idx++;
        std::unique_lock<std::mutex> lk(mtx);
        cv_space.wait(lk, [&]{ return queue.size() < max_items || first_err; });
        if (first_err) return;
        queue.push_back(std::move(wi));
        cv_item.notify_one();
    };
    try {
        klibpp::SeqStreamIn ksi(filename.c_str());
        klibpp::KSeq rec;
        WorkItem acc;
        int64_t acc_bp {0};
        while (ksi >> rec) {
            if (rec.seq.length() < k) continue;             // matches FileSkmerator::init_record
            auto seq = std::make_shared<const std::string>(rec.seq);
            const int64_t L {static_cast<int64_t>(seq->size())};
            if (L <= target_chunk) {
                acc.whole.push_back(seq);
                acc_bp += L;
                if (acc_bp >= target_chunk) { emit(std::move(acc)); acc = WorkItem{}; acc_bp = 0; }
            } else {
                if (!acc.whole.empty()) { emit(std::move(acc)); acc = WorkItem{}; acc_bp = 0; }
                for (int64_t a {0}; a < L; a += target_chunk) {
                    WorkItem wi;
                    wi.big = seq;
                    wi.a = a;
                    wi.b = std::min<int64_t>(L, a + target_chunk);
                    emit(std::move(wi));
                    if (first_err) break;
                }
            }
            if (first_err) break;
        }
        if (!acc.whole.empty()) emit(std::move(acc));
    } catch (...) {
        std::lock_guard<std::mutex> lk(mtx);
        if (!first_err) first_err = std::current_exception();
    }
    {
        std::lock_guard<std::mutex> lk(mtx);
        done_reading = true;
    }
    cv_item.notify_all();
    sink.set_total(idx);

    sink.run();   // drain output in item order on this thread

    for (std::thread& w : workers) w.join();
    if (first_err) std::rethrow_exception(first_err);
}

// Sequential file query for the dual-width (gen >= store) path: parse at the full `gen` width,
// route on the full minimizer, down-convert to `store`, search, and stream results in input order.
// Byte-identical to BucketedSkmerListReader::query() when gen == store. Used by the CLI for -t 1
// and whenever the record width is narrower than the generation width.
template<typename gen, typename store = gen>
void sequential_query(BucketedSkmerListReader<store>& reader, const std::string& filename,
                      std::ostream& os) {
    using namespace parallel_detail;
    const uint64_t k = reader.k(), m = reader.m(), b = reader.quotient_bits();
    km::SkmerManipulator<gen> manip{k, m};
    km::FileSkmerator<gen> file_skmerator{manip, filename};

    std::vector<uint8_t> buf;
    std::string text;
    constexpr uint64_t FLUSH {4096};
    uint64_t since_flush {0};
    for (const km::Skmer<gen> skmer : file_skmerator) {
        const uint64_t bid {reader.route_minimizer(manip.minimizer(skmer))};
        const km::Skmer<store> trunc {km::truncate_skmer<gen, store>(k, m, b, skmer)};
        reader.query_into(bid, trunc, buf);
        append_result(text, buf);
        if (++since_flush >= FLUSH) {
            os.write(text.data(), static_cast<std::streamsize>(text.size()));
            text.clear();
            since_flush = 0;
        }
    }
    if (!text.empty())
        os.write(text.data(), static_cast<std::streamsize>(text.size()));
}

} // namespace sortedlist
} // namespace km

#endif // SKLIB_PARALLEL_QUERY_HPP
