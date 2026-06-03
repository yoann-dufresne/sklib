#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <functional>
#include <cstdint>
#include <stdexcept>
#include <unistd.h>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/SkmerBucketWriter.hpp>

#ifndef SORTEDSKMERLISTBUILDER_H
#define SORTEDSKMERLISTBUILDER_H

namespace km
{
namespace sortedlist
{

// Parameters driving the construction of a sorted super-k-mer list. This is the
// algorithmic contract, free of any CLI/option dependency: the caller resolves
// paths and budgets and hands over plain values.
struct SortedListBuildParams {
    uint64_t k = 0;
    uint64_t m = 0;
    std::string input_path;        // resolved (/dev/stdin if none)
    std::string output_path;       // resolved (/dev/stdout if none)
    bool ascii = false;
    uint64_t buckets = 4096;
    uint64_t max_ram_bytes = 0;    // 0 => no adaptive bucketing
    bool has_output_file = false;  // gates the bucketed (seekable regular file) path
    std::string tmp_dir;           // empty => next to the output file
};

// Sort by the total key (m_pair, pref_size, suff_size) then drop strictly-equal
// duplicates. Exact-duplicate super-k-mers (canonicalized => bit-identical)
// contribute identical k-mers in every column and are collapsed per-column by
// the construction anyway, so removing them up front is a pure optimization and
// does not change the final k-mer set (strategy B(1)).
template<typename kuint>
void sort_and_dedup(std::vector<km::Skmer<kuint>>& v) {
    std::sort(v.begin(), v.end(), [](const km::Skmer<kuint>& a, const km::Skmer<kuint>& b) {
        if (a.m_pair == b.m_pair) {
            if (a.m_pref_size == b.m_pref_size) return a.m_suff_size < b.m_suff_size;
            return a.m_pref_size < b.m_pref_size;
        }
        return a.m_pair < b.m_pair;
    });
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

// All-in-RAM construction (historical path). Used for ASCII output and for
// non-seekable / stdout targets where the bucketed path cannot patch the count.
template<typename kuint>
void build_in_ram(const SortedListBuildParams& params) {
    const uint64_t k = params.k;
    const uint64_t m = params.m;

    km::SkmerManipulator<kuint> manip{k, m};
    km::FileSkmerator<kuint> file_skmerator{manip, params.input_path};

    std::vector<km::Skmer<kuint>> skmer_enumeration;
    for (const km::Skmer<kuint>& skmer : file_skmerator)
        skmer_enumeration.push_back(skmer);

    km::sortedlist::SortedVirtualSkmerList<kuint> sorted_list(k, m);
    sorted_list.generate_sorted_list_from_enumeration(skmer_enumeration);

    if (params.ascii)
        km::sortedlist::VirtualSkmerSerializer<kuint>::save_ascii(sorted_list, params.output_path);
    else
        km::sortedlist::VirtualSkmerSerializer<kuint>::save(sorted_list, params.output_path);
}

// Removes a temporary directory tree on scope exit (success or exception).
struct TmpDirGuard {
    std::filesystem::path dir;
    bool armed{true};
    ~TmpDirGuard() {
        if (!armed) return;
        std::error_code ec;
        std::filesystem::remove_all(dir, ec);
    }
};

// Holds a minimizer -> bucket-id mapping plus the bucket count.
template<typename kuint>
struct Bucketing {
    uint64_t n_buckets;
    std::function<uint64_t(kuint)> bucket_of;
};

// Empirical phase-2 peak bytes per *raw* (pre-dedup) super-k-mer in a bucket.
// Drives --max-ram bucket sizing. Chosen conservatively (above the measured
// worst case) so --max-ram behaves as an approximate ceiling rather than a
// target; duplicate-heavy buckets cost even less. Very small budgets still hit
// an irreducible floor (one minimizer's worth + base process RSS).
inline constexpr uint64_t BYTES_PER_SKMER_PEAK = 512;
// Cap on adaptive bucket count, to keep the temp file count sane.
inline constexpr uint64_t MAX_ADAPTIVE_BUCKETS = 1u << 16;

// Strategy A v1: fixed bucketing on the top bits of the minimizer.
template<typename kuint>
Bucketing<kuint> make_prefix_bucketing(uint64_t m, uint64_t requested_buckets) {
    const uint64_t mini_bits = 2 * m;
    uint64_t bucket_bits = 0;
    while ((uint64_t{1} << (bucket_bits + 1)) <= std::max<uint64_t>(requested_buckets, 1))
        bucket_bits++;
    const uint64_t effective_bits = std::min<uint64_t>(bucket_bits, mini_bits);
    const uint64_t n_buckets = uint64_t{1} << effective_bits;
    const uint64_t shift = mini_bits - effective_bits;
    return Bucketing<kuint>{ n_buckets, [shift, n_buckets](kuint mini) -> uint64_t {
        if (n_buckets == 1) return 0;
        return static_cast<uint64_t>(mini >> shift);
    }};
}

// Strategy A v2: a histogram pass over the minimizers builds contiguous
// minimizer-prefix intervals each holding ~`target` raw super-k-mers, so the
// load is balanced and a hot prefix gets a narrow interval. Intervals stay
// contiguous in minimizer order (so concatenation reproduces the global sort)
// and never split a single minimizer.
template<typename kuint>
Bucketing<kuint> make_adaptive_bucketing(km::SkmerManipulator<kuint>& manip,
                                         const std::string& input_path,
                                         uint64_t m, uint64_t max_ram_bytes) {
    constexpr uint64_t H = 22;                       // 2^22 * 4 B = 16 MB histogram
    const uint64_t mini_bits = 2 * m;
    const uint64_t hist_bits = std::min<uint64_t>(mini_bits, H);
    const uint64_t low = mini_bits - hist_bits;      // minimizers per cell = 2^low
    const uint64_t cells = uint64_t{1} << hist_bits;

    std::vector<uint32_t> hist(cells, 0);
    uint64_t total = 0;
    {
        // Count the full stream (no B(2) dedup here): the histogram only needs
        // relative bin weights to place boundaries, and counting every skmer keeps
        // the boundaries — hence the whole --max-ram output — identical to the
        // pre-B(2) version. (B(2)'s consecutive-dup collapse happens in phase 1.)
        km::FileSkmerator<kuint> rator{manip, input_path};
        for (const km::Skmer<kuint>& sk : rator) {
            const kuint mini = manip.minimizer(sk);
            const uint64_t cell = (low >= 64) ? 0 : static_cast<uint64_t>(mini >> low);
            if (hist[cell] != UINT32_MAX) hist[cell]++;
            total++;
        }
    }

    if (total == 0)
        return Bucketing<kuint>{ 1, [](kuint) -> uint64_t { return 0; } };

    // Target raw skmers per bucket from the RAM budget, capped so the bucket
    // count stays bounded.
    uint64_t target = std::max<uint64_t>(max_ram_bytes / BYTES_PER_SKMER_PEAK, 1);
    if (total / target > MAX_ADAPTIVE_BUCKETS)
        target = (total + MAX_ADAPTIVE_BUCKETS - 1) / MAX_ADAPTIVE_BUCKETS;

    // Greedy contiguous partition: close an interval once it reaches `target`.
    std::vector<uint64_t> starts{0};
    uint64_t acc = 0;
    for (uint64_t cell = 0; cell < cells; cell++) {
        acc += hist[cell];
        if (acc >= target && cell + 1 < cells) {
            starts.push_back(cell + 1);
            acc = 0;
        }
    }

    const uint64_t n_buckets = starts.size();
    return Bucketing<kuint>{ n_buckets,
        [low, starts = std::move(starts)](kuint mini) -> uint64_t {
            const uint64_t cell = (low >= 64) ? 0 : static_cast<uint64_t>(mini >> low);
            return static_cast<uint64_t>(
                std::upper_bound(starts.begin(), starts.end(), cell) - starts.begin() - 1);
        }};
}

// Low-memory bucketed construction (strategy A v1 + B(1), optionally A v2).
// Phase 1 streams every super-k-mer into a per-minimizer disk bucket; phase 2
// loads one bucket at a time, sorts+dedups it, runs the unchanged column
// algorithm, and appends the sorted sub-list to the output. Buckets are keyed
// by a monotone prefix/interval of the minimizer, so concatenating them in
// increasing id order reproduces the global sorted list. Peak RAM is bounded by
// the largest bucket rather than by the whole genome.
template<typename kuint>
void build_bucketed(const SortedListBuildParams& params) {
    namespace fs = std::filesystem;

    const uint64_t k = params.k;
    const uint64_t m = params.m;
    const std::string& input_path  = params.input_path;
    const std::string& output_path = params.output_path;

    km::SkmerManipulator<kuint> manip{k, m};

    // ---- choose the minimizer -> bucket mapping ----
    // Adaptive bucketing (A v2) reads the input twice; the caller only sets
    // max_ram_bytes > 0 when that is possible (a real input file). Otherwise we
    // fall back to the fixed prefix bucketing (A v1).
    const bool use_adaptive = params.max_ram_bytes > 0;
    Bucketing<kuint> bucketing = use_adaptive
        ? make_adaptive_bucketing<kuint>(manip, input_path, m, params.max_ram_bytes)
        : make_prefix_bucketing<kuint>(m, params.buckets);
    const uint64_t n_buckets = bucketing.n_buckets;
    const auto& bucket_of = bucketing.bucket_of;

    // Scale the phase-1 write buffer down for tiny RAM budgets so it does not
    // dominate the requested budget.
    size_t writer_budget = size_t{32} << 20;
    if (use_adaptive)
        writer_budget = std::min<size_t>(writer_budget,
            std::max<size_t>(params.max_ram_bytes / 4, size_t{1} << 20));

    // ---- temporary directory (unique per process), cleaned on any exit ----
    fs::path tmp_base = !params.tmp_dir.empty() ? fs::path(params.tmp_dir)
                                                : fs::path(output_path).parent_path();
    if (tmp_base.empty()) tmp_base = fs::path(".");
    fs::path tmp_dir = tmp_base / (".sskm_buckets." + std::to_string(::getpid()));
    TmpDirGuard guard{tmp_dir};

    // ---- Phase 1: stream super-k-mers into buckets ----
    // Kept in its own scope so the writer's per-bucket buffers and the sequence
    // reader are fully released before phase 2 allocates, bounding the peak to
    // max(phase-1 buffers, largest bucket) rather than their sum.
    std::vector<uint64_t> counts(n_buckets, 0);
    {
        km::sortedlist::SkmerBucketWriter<kuint> writer{tmp_dir, n_buckets, writer_budget};
        km::FileSkmerator<kuint> file_skmerator{manip, input_path};
        // B(2): collapse runs of identical consecutive super-k-mers (tandem repeats
        // emit the same canonical skmer in a row) before they hit disk. Pure O(1)
        // filter; the per-bucket sort+unique (B(1)) drops the rest, so the output is
        // unchanged — this only shrinks the temporary bucket files and the phase-2
        // load.
        km::Skmer<kuint> prev{};
        bool have_prev = false;
        for (const km::Skmer<kuint>& sk : file_skmerator) {
            if (have_prev && sk == prev) continue;
            const kuint mini = manip.minimizer(sk);
            writer.append(bucket_of(mini), sk);
            prev = sk;
            have_prev = true;
        }
        writer.close();
        for (uint64_t id{0}; id < n_buckets; id++) counts[id] = writer.bucket_count(id);
    }

    // ---- Phase 2: build + append each bucket's sorted sub-list ----
    {
        std::ofstream out(output_path, std::ios::binary | std::ios::trunc);
        if (out.fail())
            throw std::runtime_error("Error opening output file for writing: " + output_path);
        km::sortedlist::VirtualSkmerSerializer<kuint>::write_header(out, k, m, 0);
        if (out.fail())
            throw std::runtime_error("Error writing header to file: " + output_path);

        uint64_t total = 0;
        for (uint64_t id{0}; id < n_buckets; id++) {
            if (counts[id] == 0) continue;
            const fs::path bpath = tmp_dir / km::sortedlist::SkmerBucketWriter<kuint>::bucket_filename(id);
            std::vector<km::Skmer<kuint>> vec =
                km::sortedlist::SkmerBucketWriter<kuint>::load_bucket(bpath);
            if (vec.empty()) continue;

            sort_and_dedup(vec);

            km::sortedlist::SortedVirtualSkmerList<kuint> sub(k, m);
            sub.generate_sorted_list_from_enumeration(vec);
            km::sortedlist::VirtualSkmerSerializer<kuint>::append_payload(out, sub.get_list());
            if (out.fail())
                throw std::runtime_error("Error writing skmers to file: " + output_path);
            total += sub.size();

            std::error_code ec;
            fs::remove(bpath, ec); // free disk as we go
        }

        km::sortedlist::VirtualSkmerSerializer<kuint>::patch_count(out, total);
        if (out.fail())
            throw std::runtime_error("Error patching count in file: " + output_path);
        out.close();
    }
}

// Dispatch to the right construction path. Bucketing needs a seekable regular
// file (to patch the header count) and the raw binary layout; ASCII output and
// stdout fall back to the all-in-RAM path.
template<typename kuint>
void build_sorted_list(const SortedListBuildParams& params) {
    const bool can_bucket = !params.ascii && params.has_output_file;
    if (can_bucket)
        build_bucketed<kuint>(params);
    else
        build_in_ram<kuint>(params);
}

} // namespace sortedlist
} // namespace km

#endif // SORTEDSKMERLISTBUILDER_H
