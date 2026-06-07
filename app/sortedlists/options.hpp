#pragma once

#include <string>
#include <optional>
#include <cstdint>

// Data contract shared between the CLI parsing layer (cli.hpp) and the command
// implementations (commands.hpp). Holds the validated options for each
// subcommand, decoupled from how they were obtained.

struct ConstructOptions {
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    int k = 0;
    int m = 0;
    bool ascii = false;
    uint64_t buckets = 4096;
    std::optional<std::string> max_ram;
    std::optional<std::string> tmp_dir;
    unsigned int threads = 8;
};

struct QueryOptions {
    std::string list_file;
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    std::optional<std::string> sequence;
    unsigned int threads = 8;
};

struct SetOpOptions {
    // One of: intersection, union, diff, intersection_size, union_size, diff_size.
    // diff is asymmetric: A \ B (k-mers of A absent from B).
    std::string op;
    std::string list_a;                      // -a/--list-a
    std::string list_b;                      // -b/--list-b
    std::optional<std::string> output_file;  // -o/--output (required for the materializing ops)
    bool no_compact = false;                 // --no-compact: skip super-k-mer re-compaction of the result
    unsigned int threads = 8;                // -t/--threads: per-bucket parallel workers (output byte-identical)
};
