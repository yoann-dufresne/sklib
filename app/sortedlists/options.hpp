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
};

struct QueryOptions {
    std::string list_file;
    std::optional<std::string> input_file;
    std::optional<std::string> output_file;
    std::optional<std::string> sequence;
};
