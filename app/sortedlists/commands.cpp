#include "commands.hpp"

#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

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
        const uint64_t k = params.k;
        const uint64_t m = params.m;
        if (k == 0 || m == 0 || m > k) {
            std::cerr << "Require 1 <= m <= k" << std::endl;
            return 1;
        }
        // Generation needs the full 2*(2k-m) bits (the whole minimizer drives selection + bucket);
        // storage needs only 2*(2k-m)-b once the bucket id absorbs the top b φ-minimizer bits.
        // Pick the smallest precompiled integer width for each, then instantiate the build.
        const uint64_t b = km::sortedlist::quotient_bits_for(params);
        const uint64_t gen_w   = km::sortedlist::select_width_bytes(2 * (2 * k - m));
        const uint64_t store_w = km::sortedlist::select_width_bytes(2 * (2 * k - m) - b);
        km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
            km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
                km::sortedlist::build_sorted_list<gen, store>(params, b);
            });
        });
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}

int run_query(const QueryOptions& opts) {
    // Peek the header to learn k/m and the stored record width before instantiating any templated
    // reader. The generation width is derived from k/m (it must hold the full 2*(2k-m) bits); the
    // storage width and quotient bits come from the file.
    km::sortedlist::ListHeaderInfo hdr;
    try {
        hdr = km::sortedlist::read_list_header(opts.list_file);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

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

    // Generation width from k/m, but never narrower than the record width: queries are generated at
    // `gen` then down-converted to `store`, so we need gen >= store (a legacy 8-byte file with a
    // tiny k whose generation alone would fit 4 bytes is the one case where this max bites).
    const uint64_t gen_w = std::max<uint64_t>(
        km::sortedlist::select_width_bytes(2 * (2 * hdr.k - hdr.m)), hdr.store_width_bytes);
    const uint64_t store_w = hdr.store_width_bytes;

    try {
        km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
            km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
                // The reader routes each query to its minimizer-prefix bucket and loads only that
                // bucket's sub-list from disk on demand.
                auto reader = km::sortedlist::BucketedSkmerListReader<store>::open(opts.list_file);

                if (opts.input_file) {
                    // One reader+bucketize producer plus parallel consumers when a thread budget is
                    // given; sequential for a single thread. Output stays in input order either way.
                    if (opts.threads >= 2)
                        km::sortedlist::parallel_query<gen, store>(reader, *opts.input_file, *os, opts.threads);
                    else
                        km::sortedlist::sequential_query<gen, store>(reader, *opts.input_file, *os);
                } else {
                    // Single-sequence path: parse at gen width, route on the full minimizer, then
                    // down-convert each skmer to the stored width before searching.
                    const uint64_t k = reader.k(), m = reader.m(), b = reader.quotient_bits();
                    km::SkmerManipulator<gen> manip{k, m};
                    std::string seq = *opts.sequence;
                    km::SeqSkmerator<gen> seq_skmerator{manip, seq};

                    std::vector<std::vector<uint8_t>> results;
                    std::vector<uint8_t> buf;
                    for (const km::Skmer<gen>& s : seq_skmerator) {
                        const uint64_t bid {reader.bucket_of_phi_min(static_cast<uint64_t>(manip.minimizer(s)))};
                        const km::Skmer<store> trunc {km::truncate_skmer<gen, store>(k, m, b, s)};
                        reader.query_into(bid, trunc, buf);
                        results.push_back(buf);
                    }
                    km::sortedlist::util::print_query_results(results, *os);
                }
            });
        });
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
