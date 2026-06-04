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

    // -----------------------
    // experiment subcommand (EXPERIMENTAL)
    // -----------------------
    auto experiment = app.add_subcommand("experiment",
        "EXPERIMENTAL: real-nucleotide partial sort that groups incomparable super-k-mers "
        "(those orderable only via the absent-nucleotide mask) and dumps the grouped list as ASCII.");

    ExperimentOptions experiment_opts;

    experiment->add_option("-f,--file", experiment_opts.input_file,
        "Input FASTA file (plain or gzip-compressed). Reads from stdin if omitted.");

    experiment->add_option("-o,--output", experiment_opts.output_file,
        "Output ASCII file. Writes to stdout if omitted.");

    experiment->add_option("-k,--kmer-size", experiment_opts.k,
        "k-mer length in nucleotides (1 <= k <= 32 with the default 64-bit backend).")
        ->required();

    experiment->add_option("-m,--minimizer-size", experiment_opts.m,
        "Minimizer length in nucleotides (1 <= m <= k).")
        ->required();

    experiment->footer(
        "Example:\n"
        "  sskm experiment -k 7 -m 3 -f toy.fa -o groups.txt\n"
    );

    experiment->callback([&]() {
        result.experiment = experiment_opts;
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
