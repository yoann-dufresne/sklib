#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <cstdint>
#include <algorithm>

#include <io/Skmer.hpp>

#ifndef SKMERBUCKETWRITER_H
#define SKMERBUCKETWRITER_H

namespace km
{
namespace sortedlist
{

/** Streams super-k-mers into per-minimizer disk buckets during construction.
 *
 * Each bucket is a headerless file of raw Skmer records. Records are first
 * accumulated in a per-bucket in-memory buffer and flushed (append) once the
 * buffer fills, bounding phase-1 RAM to roughly `buffer_bytes_total`.
 *
 * Flushing opens the bucket file in append mode and closes it immediately, so
 * at most one file descriptor is held at a time regardless of bucket count
 * (avoids hitting `ulimit -n` for large bucket counts).
 *
 * Buckets are addressed by a contiguous id in [0, n_buckets); the phase-2
 * reader iterates them in increasing id order, which (for a monotone bucketing
 * of the minimizer) reproduces the global sorted order on concatenation.
 */
template<typename kuint>
class SkmerBucketWriter
{
public:
    SkmerBucketWriter(std::filesystem::path tmp_dir, uint64_t n_buckets,
                      size_t buffer_bytes_total = (size_t{32} << 20))
        : m_tmp_dir(std::move(tmp_dir)), m_n_buckets(n_buckets)
    {
        std::filesystem::create_directories(m_tmp_dir);

        // Per-bucket flush threshold, in records. Spread the global budget over
        // the buckets but keep a sane floor so tiny budgets still batch writes.
        const size_t per_bucket_bytes = std::max<size_t>(
            buffer_bytes_total / std::max<uint64_t>(n_buckets, 1), size_t{8} << 10);
        m_flush_records = std::max<size_t>(per_bucket_bytes / sizeof(Skmer<kuint>), 1);

        m_buffers.resize(n_buckets);
        m_counts.assign(n_buckets, 0);
    }

    void append(uint64_t bucket_id, const Skmer<kuint>& sk)
    {
        std::vector<Skmer<kuint>>& buf = m_buffers[bucket_id];
        buf.push_back(sk);
        if (buf.size() >= m_flush_records)
            flush_bucket(bucket_id);
    }

    // Flush every non-empty buffer to disk, then release the buffer memory so it
    // does not coexist with phase-2 per-bucket allocations. Bucket paths/counts
    // remain available afterwards. Idempotent.
    void close()
    {
        for (uint64_t id{0}; id < m_n_buckets; id++)
            if (!m_buffers[id].empty())
                flush_bucket(id);
        std::vector<std::vector<Skmer<kuint>>>().swap(m_buffers);
    }

    uint64_t n_buckets() const { return m_n_buckets; }

    // Deterministic bucket file name for an id (no directory). Phase 2 can
    // reconstruct paths from the tmp dir + id without keeping the writer alive.
    static std::string bucket_filename(uint64_t id)
    {
        char name[32];
        std::snprintf(name, sizeof(name), "bucket_%05llu.bin",
                      static_cast<unsigned long long>(id));
        return std::string(name);
    }

    std::filesystem::path bucket_path(uint64_t id) const
    {
        return m_tmp_dir / bucket_filename(id);
    }

    // Number of records routed to a bucket (0 => file never created).
    uint64_t bucket_count(uint64_t id) const { return m_counts[id]; }

    /** Read a whole bucket file back into a vector. Returns empty if absent. */
    static std::vector<Skmer<kuint>> load_bucket(const std::filesystem::path& path)
    {
        std::vector<Skmer<kuint>> out;
        std::error_code ec;
        const auto bytes = std::filesystem::file_size(path, ec);
        if (ec || bytes == 0) return out;

        if (bytes % sizeof(Skmer<kuint>) != 0)
            throw std::runtime_error("Corrupt bucket file (size not a multiple of sizeof(Skmer)): " + path.string());

        std::ifstream in(path, std::ios::binary);
        if (in.fail())
            throw std::runtime_error("Error opening bucket file for reading: " + path.string());

        out.resize(bytes / sizeof(Skmer<kuint>));
        in.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(bytes));
        if (in.fail())
            throw std::runtime_error("Error reading bucket file: " + path.string());
        return out;
    }

private:
    void flush_bucket(uint64_t id)
    {
        std::vector<Skmer<kuint>>& buf = m_buffers[id];
        std::ofstream out(bucket_path(id), std::ios::binary | std::ios::app);
        if (out.fail())
            throw std::runtime_error("Error opening bucket file for writing: " + bucket_path(id).string());
        out.write(reinterpret_cast<const char*>(buf.data()),
                  static_cast<std::streamsize>(buf.size() * sizeof(Skmer<kuint>)));
        if (out.fail())
            throw std::runtime_error("Error writing bucket file: " + bucket_path(id).string());
        m_counts[id] += buf.size();
        buf.clear(); // keep capacity for reuse
    }

    std::filesystem::path m_tmp_dir;
    uint64_t m_n_buckets;
    size_t m_flush_records;
    std::vector<std::vector<Skmer<kuint>>> m_buffers;
    std::vector<uint64_t> m_counts;
};

} // namespace sortedlist
} // namespace km

#endif
