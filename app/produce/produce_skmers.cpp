// sskm-produce -- isolated super-k-mer producer benchmark / verification harness.
//
// It runs ONLY the producer half of `sskm construct`'s Phase 1: it iterates a FileSkmerator over a
// FASTA (which drives SeqSkmerator::Iterator::operator++ and the SkmerManipulator per-nucleotide
// math) and consumes every yielded super-k-mer, with NO minimizer-bucket routing, sort, dedup or
// disk writer. CONSTRUCT_SCALING_DIAG.md showed that this producer is the serial floor of the build
// (~1/3 of the work, ~73-78% of wall at the default -t8) and that, inside it, FASTA parse is ~0% and
// I/O ~0.5% -- so this binary's wall time is essentially the SkmerManipulator compute we want to
// optimize, measured in isolation.
//
// To make the measurement trustworthy it also produces, for free, an order-sensitive FNV-1a digest
// of the whole skmer stream (see io/SkmerDigest.hpp): the digest is the anti-dead-code-elimination
// sink AND the bit-exact regression fingerprint. A monothread optimization of the producer is
// correct iff this digest is unchanged.
//
// The generation integer width is chosen exactly as the build does
// (select_width_bytes(2*(2k-m)) + dispatch_width_bytes), so the math here is the same instantiation
// that construct runs (e.g. k=21,m=11 -> 62 bits -> uint32_t pair).

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <string>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <io/SkmerDigest.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace
{

struct Args {
    uint64_t k = 0;
    uint64_t m = 0;
    std::string file;
    bool count_only = false; // --count: minimal sink (lower bound); default folds the digest
    int reps = 1;            // best (min) wall over reps is reported (warm-cache micro-bench)
};

[[noreturn]] void usage(const char* prog, int code)
{
    std::fprintf(stderr,
        "usage: %s -k <k> -m <m> -f <fasta> [--digest|--count] [--reps N]\n"
        "  Isolated super-k-mer producer: iterate FileSkmerator over <fasta> and consume each\n"
        "  yielded super-k-mer. No bucketing / sort / compaction.\n"
        "    --digest (default) order-sensitive FNV-1a fingerprint of the skmer stream (bit-exact\n"
        "                       regression check) + throughput.\n"
        "    --count            count only (minimal sink) -> lower bound of producer cost.\n"
        "    --reps N           run N passes; report the fastest (default 1).\n"
        "  Emits one machine-parseable line on stdout:\n"
        "    [producer] k=.. m=.. width_bytes=.. skmers=N digest=0x.. wall_s=.. Mskmer_s=..\n",
        prog);
    std::exit(code);
}

uint64_t parse_u64(const char* s)
{
    return std::strtoull(s, nullptr, 10);
}

} // namespace

int main(int argc, char** argv)
{
    Args args;

    for (int i = 1; i < argc; ++i) {
        const std::string a {argv[i]};
        auto value_of = [&](const char* name) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", name); std::exit(2); }
            return argv[++i];
        };
        if      (a == "-k")              args.k = parse_u64(value_of("-k"));
        else if (a == "-m")              args.m = parse_u64(value_of("-m"));
        else if (a == "-f")              args.file = value_of("-f");
        else if (a == "--digest")        args.count_only = false;
        else if (a == "--count")         args.count_only = true;
        else if (a == "--reps")          args.reps = static_cast<int>(parse_u64(value_of("--reps")));
        else if (a == "-h" || a == "--help") usage(argv[0], 0);
        // glued short forms (-k21, -m11, -f path), matching `sskm` ergonomics
        else if (a.rfind("-k", 0) == 0 && a.size() > 2) args.k = parse_u64(a.c_str() + 2);
        else if (a.rfind("-m", 0) == 0 && a.size() > 2) args.m = parse_u64(a.c_str() + 2);
        else if (a.rfind("-f", 0) == 0 && a.size() > 2) args.file = a.substr(2);
        else { std::fprintf(stderr, "unknown argument: %s\n", a.c_str()); usage(argv[0], 2); }
    }

    if (args.k == 0 || args.m == 0 || args.m > args.k || args.file.empty()) {
        std::fprintf(stderr, "error: require 1 <= m <= k and -f <fasta>\n");
        usage(argv[0], 2);
    }
    if (args.reps < 1) args.reps = 1;

    // Same width the build's generation phase uses: the full 2*(2k-m) interleaved bits.
    const uint64_t gen_w = km::sortedlist::select_width_bytes(2 * (2 * args.k - args.m));

    double  best_s  = -1.0;
    uint64_t count  = 0;
    uint64_t digest = 0;
    const bool count_only = args.count_only;

    for (int r = 0; r < args.reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        uint64_t c = 0;
        uint64_t h = km::SKMER_DIGEST_INIT;
        km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
            km::SkmerManipulator<gen> manip{args.k, args.m};
            km::FileSkmerator<gen>    rator{manip, args.file};
            uint64_t cc = 0;
            uint64_t hh = km::SKMER_DIGEST_INIT;
            if (count_only) {
                for (const km::Skmer<gen>& sk : rator) { (void)sk; ++cc; }
            } else {
                for (const km::Skmer<gen>& sk : rator) { km::skmer_digest_update(hh, sk); ++cc; }
            }
            c = cc;
            h = hh;
        });
        const auto t1 = std::chrono::steady_clock::now();
        const double s = std::chrono::duration<double>(t1 - t0).count();
        if (best_s < 0.0 || s < best_s) best_s = s;
        count  = c;
        digest = h;
    }

    const double mskmer_s = (best_s > 0.0) ? (static_cast<double>(count) / best_s / 1e6) : 0.0;
    std::printf("[producer] k=%llu m=%llu width_bytes=%llu skmers=%llu digest=0x%016llx "
                "wall_s=%.4f Mskmer_s=%.2f\n",
        static_cast<unsigned long long>(args.k), static_cast<unsigned long long>(args.m),
        static_cast<unsigned long long>(gen_w), static_cast<unsigned long long>(count),
        static_cast<unsigned long long>(digest), best_s, mskmer_s);
    return 0;
}
