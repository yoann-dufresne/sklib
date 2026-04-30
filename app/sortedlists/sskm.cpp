#include <CLI/CLI.hpp>
#include <string>
#include <optional>
#include <cstdint>
#include <vector>

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
    CLI::App app{"sskm tool"};

    CLIResult result;

    // -----------------------
    // construct subcommand
    // -----------------------
    auto construct = app.add_subcommand("construct",
        "Construct a sorted super-kmer list");

    ConstructOptions construct_opts;

    construct->add_option("-f,--file", construct_opts.input_file,
        "Input fasta file (default: stdin)");

    construct->add_option("-o,--output", construct_opts.output_file,
        "Output file (default: stdout)");

    construct->add_option("-k", construct_opts.k, "k-mer size")
        ->required();

    construct->add_option("-m", construct_opts.m, "minimizer size")
        ->required();

    construct->add_flag("--ascii", construct_opts.ascii,
        "Write output in human-readable ASCII format instead of binary");

    construct->callback([&]() {
        result.construct = construct_opts;
    });

    // -----------------------
    // query subcommand
    // -----------------------
    auto query = app.add_subcommand("query",
        "Query kmers from a list");

    QueryOptions query_opts;

    query->add_option("-l,--list", query_opts.list_file,
        "Input ssk list")
        ->required();

    query->add_option("-i,--input", query_opts.input_file,
        "Input fasta file");

    query->add_option("-o,--output", query_opts.output_file,
        "Output file (default: stdout)");

    query->add_option("sequence", query_opts.sequence,
        "Sequence to query (if no -i)");

    // Constraints
    query->require_option(1, 2); // either sequence, or -i

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


int main(int argc, char* argv[]) {
    auto const parsed {parse_cli(argc, argv)};

    if (parsed.construct) {
        return run_construct(*parsed.construct);
    }

    return 0;
}
