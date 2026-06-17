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
    // Single-op (legacy) mode: one of intersection, union, diff, xor, intersection_size, union_size,
    // diff_size, xor_size. diff is asymmetric: A \ B (k-mers of A absent from B); xor is the symmetric
    // difference A △ B (k-mers in exactly one list). Empty in combined mode.
    std::string op;
    std::string list_a;                      // -a/--list-a
    std::string list_b;                      // -b/--list-b
    std::optional<std::string> output_file;  // -o/--output (required for the single-op materializing ops)
    bool no_compact = false;                 // --no-compact: skip super-k-mer re-compaction of the result
    unsigned int threads = 8;                // -t/--threads: per-bucket parallel workers (output byte-identical)

    // Combined (single-pass) mode: produce any subset of these four results in ONE pass over A and B,
    // plus (with --sizes) all four cardinalities. Mutually exclusive with --op.
    std::optional<std::string> inter_out;    // --inter-out   (A ∩ B)
    std::optional<std::string> union_out;    // --union-out   (A ∪ B)
    std::optional<std::string> diff_ab_out;  // --diff-ab-out (A \ B)
    std::optional<std::string> diff_ba_out;  // --diff-ba-out (B \ A)
    std::optional<std::string> xor_out;      // --xor-out     (A △ B)
    bool sizes = false;                      // --sizes: print all cardinalities (single pass)
};
