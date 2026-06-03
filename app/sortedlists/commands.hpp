#pragma once

#include "options.hpp"

// Command implementations: the algorithmic work behind each subcommand. The CLI
// parsing layer hands over a validated options block and these drive the
// construction / query of sorted super-k-mer lists. Return 0 on success, a
// non-zero exit code on failure.

int run_construct(const ConstructOptions& opts);
int run_query(const QueryOptions& opts);
