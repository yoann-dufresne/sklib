#pragma once

#include <optional>

#include "options.hpp"

// Result of parsing the command line: exactly one of the subcommand option
// blocks is populated (the chosen subcommand).
struct CLIResult {
    std::optional<ConstructOptions> construct;
    std::optional<QueryOptions> query;
    std::optional<ExperimentOptions> experiment;
};

// Parse argv into a CLIResult. On a parse/validation error this exits the
// process (via CLI11's exit handling), matching standard CLI behaviour.
CLIResult parse_cli(int argc, char** argv);
