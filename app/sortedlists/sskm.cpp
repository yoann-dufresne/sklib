#include <CLI/CLI.hpp>
#include <string>
#include <optional>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
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

// Low-memory bucketed construction (strategy A v1 + B(1)).
// Phase 1 streams every super-k-mer into a per-minimizer disk bucket; phase 2
// loads one bucket at a time, sorts+dedups it, runs the unchanged column
// algorithm, and appends the sorted sub-list to the output. Buckets are keyed
// by a monotone prefix of the minimizer, so concatenating them in increasing
// id order reproduces the global sorted list. Peak RAM is bounded by the
// largest bucket rather than by the whole genome.
static int run_construct_bucketed(const ConstructOptions& opts) {
    using kuint = uint64_t;
    namespace fs = std::filesystem;

    const uint64_t k = static_cast<uint64_t>(opts.k);
    const uint64_t m = static_cast<uint64_t>(opts.m);
    const std::string input_path  = opts.input_file.value_or("/dev/stdin");
    const std::string output_path = *opts.output_file;

    // ---- bucketing function: keep the top `effective_bits` of the minimizer ----
    const uint64_t mini_bits = 2 * m;
    uint64_t bucket_bits = 0;
    while ((uint64_t{1} << (bucket_bits + 1)) <= std::max<uint64_t>(opts.buckets, 1))
        bucket_bits++;
    const uint64_t effective_bits = std::min<uint64_t>(bucket_bits, mini_bits);
    const uint64_t n_buckets = uint64_t{1} << effective_bits;
    const uint64_t shift = mini_bits - effective_bits;
    auto bucket_of = [shift, n_buckets](kuint mini) -> uint64_t {
        if (n_buckets == 1) return 0;
        return static_cast<uint64_t>(mini >> shift);
    };

    // ---- temporary directory (unique per process), cleaned on any exit ----
    fs::path tmp_base = opts.tmp_dir ? fs::path(*opts.tmp_dir)
                                     : fs::path(output_path).parent_path();
    if (tmp_base.empty()) tmp_base = fs::path(".");
    fs::path tmp_dir = tmp_base / (".sskm_buckets." + std::to_string(::getpid()));
    TmpDirGuard guard{tmp_dir};

    km::SkmerManipulator<kuint> manip{k, m};

    // ---- Phase 1: stream super-k-mers into buckets ----
    // Kept in its own scope so the writer's per-bucket buffers and the sequence
    // reader are fully released before phase 2 allocates, bounding the peak to
    // max(phase-1 buffers, largest bucket) rather than their sum.
    std::vector<uint64_t> counts(n_buckets, 0);
    {
        km::sortedlist::SkmerBucketWriter<kuint> writer{tmp_dir, n_buckets};
        km::FileSkmerator<kuint> file_skmerator{manip, input_path};
        for (const km::Skmer<kuint>& sk : file_skmerator) {
            const kuint mini = manip.minimizer(sk);
            writer.append(bucket_of(mini), sk);
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
