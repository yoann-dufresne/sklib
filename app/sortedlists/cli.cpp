#include "cli.hpp"

#include <CLI/CLI.hpp>
#include <string>

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
        "k-mer length in nucleotides. The record integer width is selected automatically "
        "(uint32/uint64/__uint128) from k and m; the only limit is that a skmer fits the widest "
        "256-bit pair, i.e. 2*(2k-m) <= 256 (so k up to ~63 at small m, more as m grows).")
        ->required();

    construct->add_option("-m,--minimizer-size", construct_opts.m,
        "Minimizer length in nucleotides (1 <= m <= k). Smaller m yields longer skmers (more "
        "shared central nucleotides per group); typical values sit around m ~ k/2.")
        ->required();

    construct->add_flag("--ascii", construct_opts.ascii,
        "Write the output as human-readable ASCII instead of the default binary format.");

    construct->add_option("--buckets", construct_opts.buckets,
        "Number of minimizer buckets (default 4096, rounded down to a power of two). "
        "Splits the list into that many independently-sorted sub-lists, partitioned by the "
        "high-order bits of the minimizer; query routes to one sub-list. More buckets => "
        "lower peak construction RAM and a faster query. Use 1 for a single list. "
        "Binary output to a regular file only.");

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

    query->add_option("-t,--threads", query_opts.threads,
        "Total worker threads for file (-i) queries (default 8): 1 reads and bucketizes the input, "
        "the rest query buckets in parallel. Output stays in input order. 1 = sequential.")
        ->check(CLI::PositiveNumber);

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

    // -----------------------
    // setop subcommand
    // -----------------------
    auto setop = app.add_subcommand("setop",
        "Set operations between two sorted skmer lists A and B. Both must be built with the same k, m "
        "and --buckets. diff is asymmetric: A \\ B.");

    SetOpOptions setop_opts;

    setop->add_option("--op", setop_opts.op,
        "Operation: intersection | union | diff | intersection_size | union_size | diff_size. "
        "The *_size variants print only the cardinality and never materialize a list.")
        ->required()
        ->check(CLI::IsMember({"intersection", "union", "diff",
                               "intersection_size", "union_size", "diff_size"}));

    setop->add_option("-a,--list-a", setop_opts.list_a,
        "First sorted skmer list (A).")->required();

    setop->add_option("-b,--list-b", setop_opts.list_b,
        "Second sorted skmer list (B).")->required();

    setop->add_option("-o,--output", setop_opts.output_file,
        "Output sorted skmer list for the result (required for intersection/union/diff; "
        "ignored by the *_size variants, which print a count to stdout).");

    setop->add_flag("--no-compact", setop_opts.no_compact,
        "Skip re-compacting the result into super-k-mers: emit one record per result k-mer. "
        "Much faster (avoids the dominant cost) but a larger output file; still a valid, "
        "queryable sorted list. Ignored by the *_size variants.");

    setop->footer(
        "Examples:\n"
        "  sskm setop --op intersection -a a.sskm -b b.sskm -o inter.sskm\n"
        "  sskm setop --op diff_size -a a.sskm -b b.sskm\n"
    );

    setop->callback([&]() {
        result.setop = setop_opts;
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

    // -----------------------
    // setop subcommand
    // -----------------------
    if (result.setop) {
        const bool is_size = setop_opts.op.size() > 5 &&
            setop_opts.op.compare(setop_opts.op.size() - 5, 5, "_size") == 0;
        if (!is_size && !setop_opts.output_file) {
            throw CLI::ValidationError("setop " + setop_opts.op + " materializes a list; provide -o/--output");
        }
    }

    return result;
}
