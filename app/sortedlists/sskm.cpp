#include <CLI/CLI.hpp>
#include <string>
#include <optional>

struct ConstructOptions {
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    int k = 0;
    int m = 0;
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
    if(query_opts.sequence && query_opts.input_file) {
        throw CLI::ValidationError("Cannot use both sequence and -i");
    }
    if(!query_opts.sequence && !query_opts.input_file) {
        throw CLI::ValidationError("Provide a sequence or -i");
    }

    return result;
}


int main(int argc, char* argv[]) {
    auto const parsed {parse_cli(argc, argv)};

    return 0;
}
