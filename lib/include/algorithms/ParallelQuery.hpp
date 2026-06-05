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
// shared: query_skmer_into()/bucket_of() are thread-safe (each bucket is loaded once under a lock,
// hits are lock-free), so consumers never re-read a bucket.

#include <cstdint>
#include <string>
#include <vector>
#include <queue>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <ostream>
#include <utility>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace km
{
namespace sortedlist
{

namespace parallel_detail
{

// A unit of work: a contiguous run of input super-k-mers, tagged with its batch sequence number so
// the sink can re-emit batches in input order regardless of which consumer finishes first.
template<typename kuint>
struct WorkBatch {
    uint64_t seq;
    std::vector<Skmer<kuint>> skmers;
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
template<typename kuint>
class WorkQueue {
public:
    explicit WorkQueue(size_t capacity) : m_capacity(capacity) {}

    void push(WorkBatch<kuint>&& batch) {
        std::unique_lock<std::mutex> lock(m_mtx);
        m_not_full.wait(lock, [&]{ return m_queue.size() < m_capacity; });
        m_queue.push(std::move(batch));
        lock.unlock();
        m_not_empty.notify_one();
    }

    bool pop(WorkBatch<kuint>& out) {
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
    std::queue<WorkBatch<kuint>> m_queue;
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

// Query every super-k-mer of `filename` against `reader` using `n_threads` total worker threads
// (1 producer + the rest as parallel consumers), writing presence vectors to `os` in input order.
// Equivalent in output to reader.query(filename, os). `batch_size` is the number of super-k-mers per
// work batch (bounds per-batch RAM and sink granularity; mainly a test/tuning knob).
template<typename kuint>
void parallel_query(BucketedSkmerListReader<kuint>& reader, const std::string& filename,
                    std::ostream& os, unsigned n_threads, uint64_t batch_size = 4096) {
    using namespace parallel_detail;

    const unsigned n_consumers = (n_threads > 1) ? (n_threads - 1) : 1;
    if (batch_size == 0) batch_size = 1;
    const size_t queue_capacity = static_cast<size_t>(n_consumers) * 4 + 1;

    WorkQueue<kuint> queue(queue_capacity);
    OrderedSink sink(os);

    std::thread producer([&]{
        km::SkmerManipulator<kuint> manip{reader.k(), reader.m()};
        km::FileSkmerator<kuint> file_skmerator{manip, filename};

        std::vector<Skmer<kuint>> batch;
        batch.reserve(batch_size);
        uint64_t seq {0};

        for (const km::Skmer<kuint> skmer : file_skmerator) {
            batch.push_back(skmer);
            if (batch.size() >= batch_size) {
                queue.push(WorkBatch<kuint>{seq++, std::move(batch)});
                batch = std::vector<Skmer<kuint>>();
                batch.reserve(batch_size);
            }
        }
        if (!batch.empty())
            queue.push(WorkBatch<kuint>{seq++, std::move(batch)});

        queue.set_done();      // let consumers exit once the queue is drained
        sink.set_total(seq);   // let the sink stop once all `seq` batches are emitted
    });

    std::vector<std::thread> consumers;
    consumers.reserve(n_consumers);
    for (unsigned t {0}; t < n_consumers; ++t) {
        consumers.emplace_back([&]{
            WorkBatch<kuint> batch;
            std::vector<uint8_t> buf; // reused result buffer (no per-query allocation)
            std::string text;         // reused per-batch output
            while (queue.pop(batch)) {
                text.clear();
                for (const Skmer<kuint>& skmer : batch.skmers) {
                    reader.query_skmer_into(skmer, buf);
                    append_result(text, buf);
                }
                sink.put(batch.seq, std::move(text));
            }
        });
    }

    sink.run(); // drain output in order on this thread until every batch is emitted

    producer.join();
    for (std::thread& c : consumers)
        c.join();
}

} // namespace sortedlist
} // namespace km

#endif // SKLIB_PARALLEL_QUERY_HPP
