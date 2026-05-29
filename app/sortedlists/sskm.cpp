#include <CLI/CLI.hpp>
#include <string>
#include <optional>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <unistd.h>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/SkmerBucketWriter.hpp>

struct ConstructOptions {
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    int k = 0;
    int m = 0;
    bool ascii = false;
    uint64_t buckets = 4096;
    std::optional<std::string> max_ram;
    std::optional<std::string> tmp_dir;
};

struct QueryOptions {
    std::string list_file;
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    std::optional<std::string> sequence;
};

struct CLIResult {
    std::optional<ConstructOptions> construct;
    std::optional<QueryOptions> query;
};

CLIResult parse_cli(int argc, char** argv) {
    CLI::App app{
        "sskm — construct and query compact sorted super-k-mer lists.\n"
        "Packs DNA into 2-bit-per-nucleotide skmers, builds a disk-backed sorted list, "
        "and answers k-mer membership queries."
    };
    app.set_version_flag("--version,-V", std::string{SKLIB_VERSION});
    app.footer(
        "Examples:\n"
        "  sskm construct -k 21 -m 11 -f genome.fa -o genome.sskm\n"
        "  sskm query -l genome.sskm -i reads.fa -o hits.txt\n"
        "  sskm query -l genome.sskm ACGTACGTACGTACGTACGTA\n"
    );

    CLIResult result;

    // -----------------------
    // construct subcommand
    // -----------------------
    auto construct = app.add_subcommand("construct",
        "Construct a sorted super-k-mer list from a FASTA input.");

    ConstructOptions construct_opts;

    construct->add_option("-f,--file", construct_opts.input_file,
        "Input FASTA file (plain or gzip-compressed). Reads from stdin if omitted.");

    construct->add_option("-o,--output", construct_opts.output_file,
        "Output sorted skmer list. Writes to stdout if omitted.");

    construct->add_option("-k,--kmer-size", construct_opts.k,
        "k-mer length in nucleotides (1 <= k <= 32 with the default 64-bit backend).")
        ->required();

    construct->add_option("-m,--minimizer-size", construct_opts.m,
        "Minimizer length in nucleotides (1 <= m <= k). Smaller m yields longer skmers.")
        ->required();

    construct->add_flag("--ascii", construct_opts.ascii,
        "Write the output as human-readable ASCII instead of the default binary format.");

    construct->add_option("--buckets", construct_opts.buckets,
        "Number of on-disk minimizer buckets for the low-memory construction path "
        "(default 4096, rounded down to a power of two). More buckets => lower peak "
        "RAM. Use 1 for a single bucket. Binary output to a regular file only.");

    construct->add_option("--max-ram", construct_opts.max_ram,
        "Approximate target peak RAM for construction (e.g. 2G, 512M). When set, a "
        "histogram pass derives adaptive minimizer buckets balanced to this budget "
        "(overrides --buckets). Requires -f <file> (the input is read twice).");

    construct->add_option("--tmp-dir", construct_opts.tmp_dir,
        "Directory for temporary bucket files (default: next to the output file).");

    construct->footer(
        "Example:\n"
        "  sskm construct -k 21 -m 11 -f genome.fa -o genome.sskm\n"
    );

    construct->callback([&]() {
        result.construct = construct_opts;
    });

    // -----------------------
    // query subcommand
    // -----------------------
    auto query = app.add_subcommand("query",
        "Query k-mers against an existing sorted skmer list. "
        "Queries come either from a FASTA file (-i) or from a single sequence given as a positional argument.");

    QueryOptions query_opts;

    query->add_option("-l,--list", query_opts.list_file,
        "Sorted skmer list produced by `sskm construct` (binary format).")
        ->required();

    query->add_option("-i,--input", query_opts.input_file,
        "FASTA file whose k-mers are extracted and looked up in the list. "
        "Mutually exclusive with the positional `sequence` argument.");

    query->add_option("-o,--output", query_opts.output_file,
        "Output file for query hits. Writes to stdout if omitted.");

    query->add_option("sequence", query_opts.sequence,
        "Single DNA sequence to query. Mutually exclusive with -i/--input.");

    query->footer(
        "Examples:\n"
        "  sskm query -l genome.sskm -i reads.fa -o hits.txt\n"
        "  sskm query -l genome.sskm ACGTACGTACGTACGTACGTA\n"
    );

    query->callback([&]() {
        result.query = query_opts;
    });

    app.require_subcommand(true);
    // -----------------------
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        std::exit(app.exit(e));
    }

    // -----------------------
    // query subcommand
    // -----------------------
    if (result.query) {
        if(query_opts.sequence && query_opts.input_file) {
            throw CLI::ValidationError("Cannot use both sequence and -i");
        }
        if(!query_opts.sequence && !query_opts.input_file) {
            throw CLI::ValidationError("Provide a sequence or -i");
        }
    }

    return result;
}


// Sort by the total key (m_pair, pref_size, suff_size) then drop strictly-equal
// duplicates. Exact-duplicate super-k-mers (canonicalized => bit-identical)
// contribute identical k-mers in every column and are collapsed per-column by
// the construction anyway, so removing them up front is a pure optimization and
// does not change the final k-mer set (strategy B(1)).
template<typename kuint>
static void sort_and_dedup(std::vector<km::Skmer<kuint>>& v) {
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
static int run_construct_in_ram(const ConstructOptions& opts) {
    using kuint = uint64_t;
    const uint64_t k = static_cast<uint64_t>(opts.k);
    const uint64_t m = static_cast<uint64_t>(opts.m);

    const std::string input_path  = opts.input_file.value_or("/dev/stdin");
    const std::string output_path = opts.output_file.value_or("/dev/stdout");

    km::SkmerManipulator<kuint> manip{k, m};
    km::FileSkmerator<kuint> file_skmerator{manip, input_path};

    std::vector<km::Skmer<kuint>> skmer_enumeration;
    for (const km::Skmer<kuint>& skmer : file_skmerator)
        skmer_enumeration.push_back(skmer);

    km::sortedlist::SortedVirtualSkmerList<kuint> sorted_list(k, m);
    sorted_list.generate_sorted_list_from_enumeration(skmer_enumeration);

    // Order-preserving sentinel fill of the absent flank slots (substrate for hole-aware
    // queries). Skipped for --ascii, which inspects only the present nucleotides.
    if (!opts.ascii)
        sorted_list.fill_absent_interpolated();

    if (opts.ascii)
        km::sortedlist::VirtualSkmerSerializer<kuint>::save_ascii(sorted_list, output_path);
    else
        km::sortedlist::VirtualSkmerSerializer<kuint>::save(sorted_list, output_path);

    return 0;
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

// Parse a human size like "2G", "512M", "4096K", or a plain byte count.
// Returns 0 on parse failure.
static uint64_t parse_size(const std::string& s) {
    if (s.empty()) return 0;
    char* end = nullptr;
    const double value = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || value < 0) return 0;
    uint64_t mult = 1;
    switch (*end) {
        case 'k': case 'K': mult = 1ull << 10; break;
        case 'm': case 'M': mult = 1ull << 20; break;
        case 'g': case 'G': mult = 1ull << 30; break;
        case 't': case 'T': mult = 1ull << 40; break;
        case '\0': mult = 1; break;
        default: return 0;
    }
    return static_cast<uint64_t>(value * static_cast<double>(mult));
}

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
static constexpr uint64_t BYTES_PER_SKMER_PEAK = 512;
// Cap on adaptive bucket count, to keep the temp file count sane.
static constexpr uint64_t MAX_ADAPTIVE_BUCKETS = 1u << 16;

// Strategy A v1: fixed bucketing on the top bits of the minimizer.
template<typename kuint>
static Bucketing<kuint> make_prefix_bucketing(uint64_t m, uint64_t requested_buckets) {
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
static Bucketing<kuint> make_adaptive_bucketing(km::SkmerManipulator<kuint>& manip,
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
static int run_construct_bucketed(const ConstructOptions& opts) {
    using kuint = uint64_t;
    namespace fs = std::filesystem;

    const uint64_t k = static_cast<uint64_t>(opts.k);
    const uint64_t m = static_cast<uint64_t>(opts.m);
    const std::string input_path  = opts.input_file.value_or("/dev/stdin");
    const std::string output_path = *opts.output_file;

    km::SkmerManipulator<kuint> manip{k, m};

    // ---- choose the minimizer -> bucket mapping ----
    // --max-ram (A v2) needs to read the input twice, so it requires -f <file>;
    // otherwise fall back to the fixed prefix bucketing (A v1).
    bool use_adaptive = false;
    uint64_t max_ram_bytes = 0;
    if (opts.max_ram) {
        max_ram_bytes = parse_size(*opts.max_ram);
        if (max_ram_bytes == 0) {
            std::cerr << "Invalid --max-ram value: " << *opts.max_ram << std::endl;
            return 1;
        }
        if (!opts.input_file) {
            std::cerr << "--max-ram requires -f <file> (the input is read twice); "
                         "falling back to --buckets " << opts.buckets << std::endl;
        } else {
            use_adaptive = true;
        }
    }

    Bucketing<kuint> bucketing = use_adaptive
        ? make_adaptive_bucketing<kuint>(manip, input_path, m, max_ram_bytes)
        : make_prefix_bucketing<kuint>(m, opts.buckets);
    const uint64_t n_buckets = bucketing.n_buckets;
    const auto& bucket_of = bucketing.bucket_of;

    // Scale the phase-1 write buffer down for tiny RAM budgets so it does not
    // dominate the requested budget.
    size_t writer_budget = size_t{32} << 20;
    if (use_adaptive)
        writer_budget = std::min<size_t>(writer_budget,
            std::max<size_t>(max_ram_bytes / 4, size_t{1} << 20));

    // ---- temporary directory (unique per process), cleaned on any exit ----
    fs::path tmp_base = opts.tmp_dir ? fs::path(*opts.tmp_dir)
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
        if (out.fail()) {
            std::cerr << "Error opening output file for writing: " << output_path << std::endl;
            return 1;
        }
        km::sortedlist::VirtualSkmerSerializer<kuint>::write_header(out, k, m, 0);
        if (out.fail()) {
            std::cerr << "Error writing header to file: " << output_path << std::endl;
            return 1;
        }

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
            // Order-preserving sentinel fill before serialization (per bucket: each bucket
            // is a contiguous minimizer range, and the fill touches only sub-minimizer bits,
            // so concatenation stays globally sorted).
            sub.fill_absent_interpolated();
            km::sortedlist::VirtualSkmerSerializer<kuint>::append_payload(out, sub.get_list());
            if (out.fail()) {
                std::cerr << "Error writing skmers to file: " << output_path << std::endl;
                return 1;
            }
            total += sub.size();

            std::error_code ec;
            fs::remove(bpath, ec); // free disk as we go
        }

        km::sortedlist::VirtualSkmerSerializer<kuint>::patch_count(out, total);
        if (out.fail()) {
            std::cerr << "Error patching count in file: " << output_path << std::endl;
            return 1;
        }
        out.close();
    }

    return 0;
}

int run_construct(const ConstructOptions& opts) {
    // Bucketing needs a seekable regular file (to patch the header count) and
    // the raw binary layout. ASCII output and stdout fall back to the
    // all-in-RAM path.
    const bool can_bucket = !opts.ascii && opts.output_file.has_value();
    if (can_bucket)
        return run_construct_bucketed(opts);
    return run_construct_in_ram(opts);
}


int run_query(const QueryOptions& opts) {
    using kuint = uint64_t;

    auto list = km::sortedlist::VirtualSkmerSerializer<kuint>::load(opts.list_file);

    std::ofstream out_file;
    std::ostream* os = &std::cout;
    if (opts.output_file) {
        out_file.open(*opts.output_file);
        if (!out_file) {
            std::cerr << "Error opening output file: " << *opts.output_file << std::endl;
            return 1;
        }
        os = &out_file;
    }

    if (opts.input_file) {
        list.query(*opts.input_file, *os);
    } else {
        km::SkmerManipulator<kuint> manip{list.k(), list.m()};
        std::string seq = *opts.sequence;
        km::SeqSkmerator<kuint> seq_skmerator{manip, seq};

        std::vector<km::Skmer<kuint>> skmers;
        for (const km::Skmer<kuint>& s : seq_skmerator) {
            skmers.push_back(s);
        }

        auto results = list.query_skmer_batch(skmers);
        km::sortedlist::util::print_query_results(results, *os);
    }

    return 0;
}


int main(int argc, char* argv[]) {
    auto const parsed {parse_cli(argc, argv)};

    if (parsed.construct) {
        return run_construct(*parsed.construct);
    }
    if (parsed.query) {
        return run_query(*parsed.query);
    }

    return 0;
}
