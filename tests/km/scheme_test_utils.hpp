#ifndef SKLIB_SCHEME_TEST_UTILS_HPP
#define SKLIB_SCHEME_TEST_UTILS_HPP

// Shared helpers for the minimizer-scheme tests (minimizer_scheme.cpp, oc_framing.cpp,
// oc_end_to_end.cpp): deterministic DNA, reverse complement, FASTA I/O, and "build an index / query
// a FASTA exactly like the CLI does" for a given minimizer scheme.

#include <cstdint>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <io/MinimizerScheme.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace sktest {

inline std::string random_dna(uint64_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    static const char alphabet[] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(n);
    for (uint64_t i = 0; i < n; ++i) s.push_back(alphabet[rng() & 3u]);
    return s;
}

inline std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r)
        switch (c) {
            case 'A': c = 'T'; break;  case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;  case 'G': c = 'C'; break;
            default: break;
        }
    return r;
}

inline std::string canonical(const std::string& s) { return std::min(s, revcomp(s)); }

inline std::string tile(const std::string& unit, size_t len) {
    std::string s;
    while (s.size() < len) s += unit;
    return s.substr(0, len);
}

inline void write_fasta(const std::string& path, const std::vector<std::string>& seqs) {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    for (size_t i = 0; i < seqs.size(); ++i) f << ">s" << i << "\n" << seqs[i] << "\n";
}

inline std::string read_file_bytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

inline km::SchemeParams random_scheme() { return km::SchemeParams{}; }
inline km::SchemeParams oc_scheme(uint64_t m, uint64_t s = 0) { return km::SchemeParams::open_closed(m, s); }

// Build an index exactly like `sskm construct`: generation width from (k, m, scheme), storage width
// from the quotient bits. has_output_file=false takes the in-RAM path; max_ram_bytes>0 the adaptive one.
inline void build_index(uint64_t k, uint64_t m, const km::SchemeParams& scheme,
                        const std::string& in_path, const std::string& out_path,
                        uint64_t buckets = 4096, unsigned n_threads = 1,
                        bool has_output_file = true, uint64_t max_ram_bytes = 0) {
    km::sortedlist::SortedListBuildParams p;
    p.k = k;
    p.m = m;
    p.input_path = in_path;
    p.output_path = out_path;
    p.buckets = buckets;
    p.has_output_file = has_output_file;
    p.max_ram_bytes = max_ram_bytes;
    p.n_threads = n_threads;
    p.scheme = scheme;
    const uint64_t b  = km::sortedlist::quotient_bits_for(p);
    const uint64_t gw = km::sortedlist::select_generation_width_bytes(k, m, scheme);
    const uint64_t sw = km::sortedlist::select_width_bytes(2 * (2 * k - m) - b);
    km::sortedlist::dispatch_width_bytes(gw, [&]<typename gen>() {
        km::sortedlist::dispatch_width_bytes(sw, [&]<typename store>() {
            km::sortedlist::build_sorted_list<gen, store>(p, b);
        });
    });
}

// Query a FASTA against an index exactly like `sskm query -i` (widths and scheme from the header).
// n_threads == 1 -> sequential driver, else the parallel one. Returns the raw output text.
inline std::string query_file(const std::string& index_path, const std::string& query_fa, unsigned n_threads = 1) {
    const km::sortedlist::ListHeaderInfo h = km::sortedlist::read_list_header(index_path);
    const uint64_t gw = std::max<uint64_t>(
        km::sortedlist::select_generation_width_bytes(h.k, h.m, h.scheme), h.store_width_bytes);
    std::string out;
    km::sortedlist::dispatch_width_bytes(gw, [&]<typename gen>() {
        km::sortedlist::dispatch_width_bytes(h.store_width_bytes, [&]<typename store>() {
            auto reader = km::sortedlist::BucketedSkmerListReader<store>::open(index_path);
            std::stringstream ss;
            if (n_threads >= 2) km::sortedlist::parallel_query<gen, store>(reader, query_fa, ss, n_threads);
            else                km::sortedlist::sequential_query<gen, store>(reader, query_fa, ss);
            out = ss.str();
        });
    });
    return out;
}

inline std::vector<std::string> lines_of(const std::string& text) {
    std::vector<std::string> lines;
    std::stringstream ss(text);
    std::string line;
    while (std::getline(ss, line)) lines.push_back(line);
    return lines;
}

} // namespace sktest

#endif // SKLIB_SCHEME_TEST_UTILS_HPP
