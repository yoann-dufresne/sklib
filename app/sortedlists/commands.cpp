#include "commands.hpp"

#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>

// Parse a human size like "2G", "512M", "4096K", or a plain byte count.
// Returns 0 on parse failure. CLI-input concern, hence kept app-side.
static uint64_t parse_size(const std::string& s) {
    if (s.empty()) return 0;
    char* end = nullptr;
    const double value = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || value < 0) return 0;
    uint64_t mult = 1;
    switch (*end) {
        case 'k': case 'K': mult = 1ull << 10; break;
        case 'm': case 'M': mult = 1ull << 20; break;
        case 'g': case 'G': mult = 1ull << 30; break;
        case 't': case 'T': mult = 1ull << 40; break;
        case '\0': mult = 1; break;
        default: return 0;
    }
    return static_cast<uint64_t>(value * static_cast<double>(mult));
}

int run_construct(const ConstructOptions& opts) {
    using kuint = uint64_t;

    km::sortedlist::SortedListBuildParams params;
    params.k = static_cast<uint64_t>(opts.k);
    params.m = static_cast<uint64_t>(opts.m);
    params.input_path  = opts.input_file.value_or("/dev/stdin");
    params.output_path = opts.output_file.value_or("/dev/stdout");
    params.ascii = opts.ascii;
    params.buckets = opts.buckets;
    params.has_output_file = opts.output_file.has_value();
    if (opts.tmp_dir) params.tmp_dir = *opts.tmp_dir;

    // --max-ram (adaptive bucketing) reads the input twice, so it requires
    // -f <file>; otherwise warn and fall back to fixed --buckets (max_ram_bytes
    // stays 0, which the builder reads as "no adaptive bucketing").
    if (opts.max_ram) {
        const uint64_t max_ram_bytes = parse_size(*opts.max_ram);
        if (max_ram_bytes == 0) {
            std::cerr << "Invalid --max-ram value: " << *opts.max_ram << std::endl;
            return 1;
        }
        if (!opts.input_file) {
            std::cerr << "--max-ram requires -f <file> (the input is read twice); "
                         "falling back to --buckets " << opts.buckets << std::endl;
        } else {
            params.max_ram_bytes = max_ram_bytes;
        }
    }

    try {
        km::sortedlist::build_sorted_list<kuint>(params);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}

int run_query(const QueryOptions& opts) {
    using kuint = uint64_t;

    // Open the bucketed sorted list. The reader routes each query to its minimizer-prefix
    // bucket and loads only that bucket's sub-list from disk on demand.
    km::sortedlist::BucketedSkmerListReader<kuint> reader = [&]{
        try {
            return km::sortedlist::BucketedSkmerListReader<kuint>::open(opts.list_file);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
    }();

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
        // One reader+bucketize producer plus parallel consumers when a thread budget is given;
        // fall back to the sequential path for a single thread. Output stays in input order.
        if (opts.threads >= 2)
            km::sortedlist::parallel_query<kuint>(reader, *opts.input_file, *os, opts.threads);
        else
            reader.query(*opts.input_file, *os);
    } else {
        km::SkmerManipulator<kuint> manip{reader.k(), reader.m()};
        std::string seq = *opts.sequence;
        km::SeqSkmerator<kuint> seq_skmerator{manip, seq};

        std::vector<km::Skmer<kuint>> skmers;
        for (const km::Skmer<kuint>& s : seq_skmerator) {
            skmers.push_back(s);
        }

        auto results = reader.query_skmer_batch(skmers);
        km::sortedlist::util::print_query_results(results, *os);
    }

    return 0;
}
