// Regression tests for the two width bugs that the m ~ k/2 test grid never reached.
//
// Everything here is driven by k - m and by m/k, the two axes the rest of the suite holds fixed:
//
//  1. WORK width. The generation integer used to be sized for the interleaved super-k-mer only
//     (2*(2k-m) bits in the PAIR), ignoring that the whole minimizer pipeline works on a single
//     WORD: minimizer() returns the pair's low word, phi/reverse_2m mix a kuint,
//     permute_minimizer_slot writes a kuint-wide psi back, and mmer_repeats rolls the central
//     m-mer in a kuint. Whenever 3m > 2k the minimizer (2m bits) overflowed that word and all of
//     them silently truncated, which
//       - biased the stored psi (its top bits were always zero, so the minimizer-prefix bucketing
//         lost that many bits),
//       - made the router shift by >= the word width once the whole bucket prefix fell in the zero
//         region: undefined behaviour, an out-of-range bucket id, and SIGSEGV in construct,
//       - and misdetected the ambiguous-minimizer case, breaking EXACTNESS: measured on ecoli,
//         k=21 m=20 answered 8 884 false positives out of 200 000 random k-mers.
//     select_generation_width_bytes now takes the max of both constraints.
//
//  2. STORE width. The per-bucket minimizer-prefix table indexes the top bits of the STORED
//     minimizer (2m - b bits) read through minimizer(), i.e. again one word of the store type.
//     For small k - m the store shrinks faster than the minimizer (k=31 m=29 keeps a 46-bit
//     minimizer in a uint32 record), the shift ran past the word width and the table lookup
//     indexed far out of bounds: SIGSEGV in query. The table is now skipped in that case.
//
// The gate below is behavioural, not structural: for a grid that spans both regimes, every k-mer
// of the indexed sequence must be reported present, and k-mers that are NOT in it must not be.

#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/ParallelQuery.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace {

std::string rand_dna(uint64_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    static const char alphabet[] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(n);
    for (uint64_t i = 0; i < n; ++i) s.push_back(alphabet[rng() & 3u]);
    return s;
}

std::string canonical(const std::string& s) {
    std::string rc;
    rc.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it)
        rc.push_back(*it == 'A' ? 'T' : *it == 'T' ? 'A' : *it == 'C' ? 'G' : 'C');
    return std::min(s, rc);
}

// Build an index for (k, m) exactly as the CLI does — same width choices, same bucketing — then
// query `qseqs` through the same dual-width driver, returning one presence line per query record.
std::vector<std::string> build_and_query(uint64_t k, uint64_t m, const std::string& ref,
                                         const std::vector<std::string>& qseqs,
                                         const std::string& tag, uint64_t buckets) {
    const std::string dir = ::testing::TempDir();
    const std::string in_path  = dir + "mw_in_"  + tag + ".fa";
    const std::string q_path   = dir + "mw_q_"   + tag + ".fa";
    const std::string out_path = dir + "mw_out_" + tag + ".sskm";
    { std::ofstream f(in_path); f << ">ref\n" << ref << "\n"; }
    { std::ofstream f(q_path);
      for (size_t i = 0; i < qseqs.size(); ++i) f << ">q" << i << "\n" << qseqs[i] << "\n"; }

    km::sortedlist::SortedListBuildParams params;
    params.k = k;
    params.m = m;
    params.input_path = in_path;
    params.output_path = out_path;
    params.buckets = buckets;
    params.has_output_file = true;
    params.n_threads = 1;

    const uint64_t b = km::sortedlist::quotient_bits_for(params);
    const uint64_t gen_w   = km::sortedlist::select_generation_width_bytes(k, m);
    const uint64_t store_w = km::sortedlist::select_width_bytes(2 * (2 * k - m) - b);

    std::string out;
    km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
        km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
            km::sortedlist::build_sorted_list<gen, store>(params, b);
            auto reader = km::sortedlist::BucketedSkmerListReader<store>::open(out_path);
            std::stringstream ss;
            km::sortedlist::sequential_query<gen, store>(reader, q_path, ss);
            out = ss.str();
        });
    });

    std::vector<std::string> lines;
    std::stringstream ss(out);
    std::string line;
    while (std::getline(ss, line)) lines.push_back(line);
    return lines;
}

// (k, m) pairs spanning both failure modes plus healthy controls. Kept small so the suite stays fast.
struct KM { uint64_t k, m; const char* why; };
const KM kGrid[] = {
    {21, 11, "controle, m ~ k/2"},
    {31, 15, "controle, m ~ k/2"},
    {31, 29, "k-m=2 : minimizer stocke (46 bits) plus large que le mot du store (uint32)"},
    {33, 31, "k-m=2"},
    {25, 23, "k-m=2 et 3m > 2k"},
    {21, 17, "3m > 2k : la baseline tronquait le minimizer de travail"},
    {21, 20, "3m > 2k, troncature maximale (la baseline donnait 4% de faux positifs)"},
    {25, 22, "3m > 2k : la baseline plantait (construct)"},
    {42, 38, "3m > 2k a 64 bits : la baseline plantait (construct)"},
    {45, 43, "3m > 2k et k-m=2"},
};

} // namespace

// Every k-mer of the indexed sequence must be found, and absent k-mers must not be. This is the
// property both bugs violated: the work-width truncation produced false positives (and a crash in
// construct), the store-width one crashed the query.
TEST(MinimizerWidth, ExactAcrossKMinusMAndLargeM) {
    const std::string ref = rand_dna(60000, 7);
    // Two bucket counts on purpose. The default 4096 spreads this reference thin (a handful of
    // records per bucket), which leaves the per-bucket minimizer-prefix table BELOW its build
    // threshold — so it would never exercise the store-width bug. 8 buckets puts thousands of
    // records in each, the table is built, and the out-of-range lookup fires.
    for (const uint64_t buckets : {uint64_t{4096}, uint64_t{8}})
    for (const KM& c : kGrid) {
        // present: every k-mer of the reference, in a handful of records; absent: random k-mers
        // (a random k-mer of length >= 21 is essentially never in a 20 kb sequence).
        std::vector<std::string> present, absent;
        for (uint64_t i = 0; i + c.k <= ref.size(); i += 2003) present.push_back(ref.substr(i, c.k));
        std::mt19937 rng(1234 + c.k);
        std::set<std::string> ref_kmers;
        for (uint64_t i = 0; i + c.k <= ref.size(); ++i) ref_kmers.insert(canonical(ref.substr(i, c.k)));
        while (absent.size() < 400) {
            const std::string cand = rand_dna(c.k, rng());
            if (!ref_kmers.count(canonical(cand))) absent.push_back(cand);
        }

        std::vector<std::string> q = present;
        q.insert(q.end(), absent.begin(), absent.end());
        const std::string tag = std::to_string(c.k) + "_" + std::to_string(c.m) + "_" + std::to_string(buckets);
        const std::vector<std::string> lines = build_and_query(c.k, c.m, ref, q, tag, buckets);
        const std::string ctx = "k=" + std::to_string(c.k) + " m=" + std::to_string(c.m) +
                                " buckets=" + std::to_string(buckets) + " (" + c.why + ")";

        ASSERT_EQ(lines.size(), q.size()) << ctx << ": une ligne par requete attendue";
        for (size_t i = 0; i < present.size(); ++i)
            EXPECT_EQ(lines[i], "1") << ctx << ": faux negatif sur un k-mer du reference";
        for (size_t i = 0; i < absent.size(); ++i)
            EXPECT_EQ(lines[present.size() + i], "0") << ctx << ": faux positif sur un k-mer absent";
    }
}

// The width selector must satisfy BOTH constraints, and must not widen anything that already
// satisfied them (which is what keeps every healthy k/m byte-identical).
TEST(MinimizerWidth, GenerationWidthSatisfiesBothConstraints) {
    using km::sortedlist::select_generation_width_bytes;
    using km::sortedlist::select_width_bytes;
    for (uint64_t k = 4; k <= 100; ++k) {
        for (uint64_t m = 2; m < k; ++m) {
            if (2 * (2 * k - m) > 512 || 4 * m > 512) continue;   // outside the compiled widths
            const uint64_t g = select_generation_width_bytes(k, m);
            EXPECT_GE(16 * g, 2 * (2 * k - m)) << "k=" << k << " m=" << m << ": le skmer doit tenir dans la paire";
            EXPECT_GE(8 * g, 2 * m)            << "k=" << k << " m=" << m << ": le minimizer doit tenir dans un mot";
            const uint64_t skmer_only = select_width_bytes(2 * (2 * k - m));
            if (2 * m <= 8 * skmer_only)
                EXPECT_EQ(g, skmer_only) << "k=" << k << " m=" << m
                                         << ": ne doit pas elargir quand le minimizer tenait deja";
        }
    }
}
