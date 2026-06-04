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
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/SkmerPartialOrder.hpp>

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

int run_experiment(const ExperimentOptions& opts) {
    using kuint = uint64_t;

    const uint64_t k = static_cast<uint64_t>(opts.k);
    const uint64_t m = static_cast<uint64_t>(opts.m);
    const std::string input_path = opts.input_file.value_or("/dev/stdin");

    km::SkmerManipulator<kuint> manip{k, m};

    // Same enumeration as the in-RAM construction path: canonicalized, φ-permuted,
    // mask-filled super-k-mers straight from the sequence iterator.
    km::FileSkmerator<kuint> file_skmerator{manip, input_path};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& sk : file_skmerator)
        enumeration.push_back(sk);

    // Run the production reconciliation (all colinear-chaining fusions, identical compacted
    // list and order) while recording where the merge would have ordered skmers using
    // absent/masked nucleotides; those incomparable skmers are surfaced as groups.
    km::sortedlist::SortedVirtualSkmerList<kuint> sorted(k, m);
    std::vector<uint32_t> group_of;
    sorted.generate_partial_list_from_enumeration(enumeration, group_of);

    // Fill the absent flank nucleotides of the non-grouped super-k-mers with interpolated
    // (virtual) values written into the unused slot bits; pref/suff are left unchanged.
    // Work on a copy of the compacted list.
    std::vector<km::Skmer<kuint>> skmers = sorted.get_list();
    std::vector<uint32_t> filled_pref, filled_suff;
    km::experiment::fill_virtual_nucleotides<kuint>(manip, skmers, filled_pref, filled_suff);

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

    km::experiment::dump_partial_ascii<kuint>(manip, skmers, group_of, filled_pref, filled_suff, *os);
    return 0;
}
