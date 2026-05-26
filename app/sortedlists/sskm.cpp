#include <CLI/CLI.hpp>
#include <string>
#include <optional>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

struct ConstructOptions {
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    int k = 0;
    int m = 0;
    bool ascii = false;
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


int run_construct(const ConstructOptions& opts) {
    using kuint = uint64_t;

    const uint64_t k = static_cast<uint64_t>(opts.k);
    const uint64_t m = static_cast<uint64_t>(opts.m);

    const std::string input_path  = opts.input_file.value_or("/dev/stdin");
    const std::string output_path = opts.output_file.value_or("/dev/stdout");

    km::SkmerManipulator<kuint> manip{k, m};
    km::FileSkmerator<kuint> file_skmerator{manip, input_path};

    km::SkmerPrettyPrinter<kuint> pp{k, m};

    std::vector<km::Skmer<kuint>> skmer_enumeration;
    for (const km::Skmer<kuint> & skmer : file_skmerator) {
        skmer_enumeration.push_back(skmer);
    }
    // exit(0);

    km::sortedlist::SortedVirtualSkmerList<kuint> sorted_list(k, m);
    sorted_list.generate_sorted_list_from_enumeration(skmer_enumeration);

    if (opts.ascii)
        km::sortedlist::VirtualSkmerSerializer<kuint>::save_ascii(sorted_list, output_path);
    else
        km::sortedlist::VirtualSkmerSerializer<kuint>::save(sorted_list, output_path);

    return 0;
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
