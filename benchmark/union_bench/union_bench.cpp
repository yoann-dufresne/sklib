// union_bench — isolated mono-thread benchmark/verify harness for km::sortedlist::set_union.
//
// It includes the library header and calls the REAL set_union (no copy of the body), so every edit
// to SetOperations.hpp / VirtualSkmer.hpp / Skmer.hpp is measured as-is, free of the CLI / FASTA /
// construction noise. Inputs are two pre-built compatible lists (A, B); the union output goes to a
// caller-chosen path (use /dev/shm to remove write-I/O noise from timings).
//
//   union_bench --a A.sskm --b B.sskm --mode bench  [--warmup 2] [--reps 9] [--out /dev/shm/o]
//   union_bench --a A.sskm --b B.sskm --mode verify  --ref O_ref.sskm  [--out /dev/shm/o]
//
// bench  : warmup W untimed calls, then R timed calls of set_union(...,/*nthreads*/1); prints
//          median / min / stddev / MAD (seconds) + output k-mers + records + Mkmer/s (machine line).
// verify : one set_union, then checks content-equivalence against the frozen reference O_ref:
//          kmers == union_size(A,B), intersection_size(O_opt,O_ref) == that count (=> equal sets),
//          and records(O_opt) <= records(O_ref) (no compaction regression). Prints PASS/FAIL.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <algorithms/SetOperations.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/WidthDispatch.hpp>

namespace {

struct Args {
    std::string a, b, ref, out {"/dev/shm/union_bench_out.sskm"};
    std::string mode {"bench"};
    std::string op {"union"};   // union | intersection | diff | xor — all share materialize_setop
    unsigned warmup {2}, reps {9};
};

// The set operation under test + its cardinality oracle (a different sink path), selected by --op.
template<typename store>
uint64_t run_setop(const std::string& op, km::sortedlist::BucketedSkmerListReader<store>& A,
                   km::sortedlist::BucketedSkmerListReader<store>& B, const std::string& out) {
    if (op == "intersection") return km::sortedlist::intersection<store>(A, B, out, /*no_compact*/ false, 1);
    if (op == "diff")         return km::sortedlist::difference<store>(A, B, out, false, 1);
    if (op == "xor")          return km::sortedlist::symmetric_difference<store>(A, B, out, false, 1);
    return km::sortedlist::set_union<store>(A, B, out, false, 1);
}
template<typename store>
uint64_t setop_size(const std::string& op, km::sortedlist::BucketedSkmerListReader<store>& A,
                    km::sortedlist::BucketedSkmerListReader<store>& B) {
    if (op == "intersection") return km::sortedlist::intersection_size<store>(A, B, 1);
    if (op == "diff")         return km::sortedlist::diff_size<store>(A, B, 1);
    if (op == "xor")          return km::sortedlist::sym_diff_size<store>(A, B, 1);
    return km::sortedlist::union_size<store>(A, B, 1);
}

[[noreturn]] void usage(const char* prog, const std::string& msg) {
    std::cerr << msg << "\nusage: " << prog
              << " --a A.sskm --b B.sskm --mode bench|verify [--ref O_ref] [--out path]"
                 " [--warmup N] [--reps N]\n";
    std::exit(2);
}

Args parse(int argc, char** argv) {
    Args a;
    for (int i {1}; i < argc; ++i) {
        const std::string k {argv[i]};
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) usage(argv[0], std::string("missing value for ") + name);
            return argv[++i];
        };
        if      (k == "--a")      a.a = need("--a");
        else if (k == "--b")      a.b = need("--b");
        else if (k == "--ref")    a.ref = need("--ref");
        else if (k == "--out")    a.out = need("--out");
        else if (k == "--mode")   a.mode = need("--mode");
        else if (k == "--op")     a.op = need("--op");
        else if (k == "--warmup") a.warmup = static_cast<unsigned>(std::stoul(need("--warmup")));
        else if (k == "--reps")   a.reps = static_cast<unsigned>(std::stoul(need("--reps")));
        else usage(argv[0], "unknown argument: " + k);
    }
    if (a.a.empty() || a.b.empty()) usage(argv[0], "--a and --b are required");
    if (a.op != "union" && a.op != "intersection" && a.op != "diff" && a.op != "xor" && a.op != "multi" && a.op != "count") usage(argv[0], "--op must be union|intersection|diff|xor|multi|count");
    if (a.op == "count" && a.mode == "verify") usage(argv[0], "--op count supports --mode bench only (counts are cross-checked by KMC / value compare)");
    if (a.mode != "bench" && a.mode != "verify") usage(argv[0], "--mode must be bench or verify");
    if (a.op == "multi" && a.mode == "verify") usage(argv[0], "--op multi supports --mode bench only (verify multi via the sskm CLI + sha256)");
    if (a.mode == "verify" && a.ref.empty()) usage(argv[0], "--mode verify requires --ref");
    if (a.mode == "bench" && a.reps == 0) usage(argv[0], "--reps must be >= 1");
    return a;
}

struct Stats { double median, min, stddev, mad; };

Stats stats_of(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const size_t n {v.size()};
    Stats s{};
    s.min = v.front();
    s.median = (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    double mean {0};
    for (double x : v) mean += x;
    mean /= static_cast<double>(n);
    double acc {0};
    for (double x : v) acc += (x - mean) * (x - mean);
    s.stddev = (n > 1) ? std::sqrt(acc / static_cast<double>(n - 1)) : 0.0;
    std::vector<double> dev;
    dev.reserve(n);
    for (double x : v) dev.push_back(std::fabs(x - s.median));
    std::sort(dev.begin(), dev.end());
    s.mad = (n % 2) ? dev[n / 2] : 0.5 * (dev[n / 2 - 1] + dev[n / 2]);
    return s;
}

// Bench the combined single-pass operator (all four relations materialized at once) — the case
// multi_setop exists for. Times one multi_setop call writing 4 lists; reports the summed output.
template<typename store>
int run_multi_bench(const Args& a) {
    auto A = km::sortedlist::BucketedSkmerListReader<store>::open(a.a);
    auto B = km::sortedlist::BucketedSkmerListReader<store>::open(a.b);
    km::sortedlist::MultiSetOpRequest req;
    req.inter_out = a.out + ".inter"; req.union_out = a.out + ".union";
    req.diff_ab_out = a.out + ".ab";  req.diff_ba_out = a.out + ".ba";

    km::sortedlist::MultiSetOpResult r;
    for (unsigned w {0}; w < a.warmup; ++w) r = km::sortedlist::multi_setop<store>(A, B, req, 1);
    std::vector<double> times; times.reserve(a.reps);
    for (unsigned i {0}; i < a.reps; ++i) {
        const auto t0 {std::chrono::steady_clock::now()};
        r = km::sortedlist::multi_setop<store>(A, B, req, 1);
        const auto t1 {std::chrono::steady_clock::now()};
        times.push_back(std::chrono::duration<double>(t1 - t0).count());
    }
    const uint64_t kmers {r.inter_kmers + r.union_kmers + r.diff_ab_kmers + r.diff_ba_kmers};
    const uint64_t records {km::sortedlist::read_list_header(*req.inter_out).count
                          + km::sortedlist::read_list_header(*req.union_out).count
                          + km::sortedlist::read_list_header(*req.diff_ab_out).count
                          + km::sortedlist::read_list_header(*req.diff_ba_out).count};
    const Stats s {stats_of(times)};
    const double mkmer_s {s.median > 0 ? static_cast<double>(kmers) / s.median / 1e6 : 0.0};
    std::cerr << "[bench] MULTI store=" << sizeof(store) << "B  kmers(4 outputs)=" << kmers
              << "  records=" << records << "  reps=" << a.reps << "  warmup=" << a.warmup << "\n"
              << "[bench] median=" << s.median << "s  min=" << s.min << "s  stddev=" << s.stddev
              << "s  MAD=" << s.mad << "s  " << mkmer_s << " Mkmer/s\n";
    std::cout << "RESULT\tmedian_s=" << s.median << "\tmin_s=" << s.min << "\tstddev_s=" << s.stddev
              << "\tmad_s=" << s.mad << "\tkmers=" << kmers << "\trecords=" << records
              << "\tmkmer_s=" << mkmer_s << "\tstore=" << sizeof(store) << std::endl;
    return 0;
}

// Bench the combined-count path (set_sizes — what `--sizes` / every *_size variant calls). The merge
// is the whole cost here (no materialization), so this is the case the identical-record skip helps most.
template<typename store>
int run_count_bench(const Args& a) {
    auto A = km::sortedlist::BucketedSkmerListReader<store>::open(a.a);
    auto B = km::sortedlist::BucketedSkmerListReader<store>::open(a.b);
    km::sortedlist::SetSizes s{};
    for (unsigned w {0}; w < a.warmup; ++w) s = km::sortedlist::set_sizes<store>(A, B, 1);
    std::vector<double> times; times.reserve(a.reps);
    for (unsigned r {0}; r < a.reps; ++r) {
        const auto t0 {std::chrono::steady_clock::now()};
        s = km::sortedlist::set_sizes<store>(A, B, 1);
        const auto t1 {std::chrono::steady_clock::now()};
        times.push_back(std::chrono::duration<double>(t1 - t0).count());
    }
    const uint64_t kmers {s.uni()};
    const Stats st {stats_of(times)};
    const double mkmer_s {st.median > 0 ? static_cast<double>(kmers) / st.median / 1e6 : 0.0};
    std::cerr << "[bench] COUNT store=" << sizeof(store) << "B  inter=" << s.inter
              << " only_a=" << s.only_a << " only_b=" << s.only_b << " union=" << s.uni()
              << " xor=" << s.sym_diff() << "  reps=" << a.reps << "\n"
              << "[bench] median=" << st.median << "s  min=" << st.min << "s  stddev=" << st.stddev
              << "s  MAD=" << st.mad << "s  " << mkmer_s << " Mkmer/s\n";
    std::cout << "RESULT\tmedian_s=" << st.median << "\tmin_s=" << st.min << "\tstddev_s=" << st.stddev
              << "\tmad_s=" << st.mad << "\tkmers=" << kmers << "\trecords=0"
              << "\tinter=" << s.inter << "\tonly_a=" << s.only_a << "\tonly_b=" << s.only_b
              << "\tmkmer_s=" << mkmer_s << "\tstore=" << sizeof(store) << std::endl;
    return 0;
}

template<typename store>
int run_bench(const Args& a) {
    if (a.op == "multi") return run_multi_bench<store>(a);
    if (a.op == "count") return run_count_bench<store>(a);
    auto A = km::sortedlist::BucketedSkmerListReader<store>::open(a.a);
    auto B = km::sortedlist::BucketedSkmerListReader<store>::open(a.b);

    uint64_t kmers {0};
    for (unsigned w {0}; w < a.warmup; ++w)
        kmers = run_setop<store>(a.op, A, B, a.out);

    std::vector<double> times;
    times.reserve(a.reps);
    for (unsigned r {0}; r < a.reps; ++r) {
        const auto t0 {std::chrono::steady_clock::now()};
        kmers = run_setop<store>(a.op, A, B, a.out);
        const auto t1 {std::chrono::steady_clock::now()};
        times.push_back(std::chrono::duration<double>(t1 - t0).count());
    }
    const uint64_t records {km::sortedlist::read_list_header(a.out).count};
    const Stats s {stats_of(times)};
    const double mkmer_s {s.median > 0 ? static_cast<double>(kmers) / s.median / 1e6 : 0.0};

    // Human line (stderr) + one machine-parseable line (stdout) for the driver script.
    std::cerr << "[bench] store=" << sizeof(store) << "B  kmers=" << kmers
              << "  records=" << records << "  reps=" << a.reps << "  warmup=" << a.warmup << "\n"
              << "[bench] median=" << s.median << "s  min=" << s.min
              << "s  stddev=" << s.stddev << "s  MAD=" << s.mad
              << "s  " << mkmer_s << " Mkmer/s\n";
    std::cout << "RESULT\tmedian_s=" << s.median << "\tmin_s=" << s.min
              << "\tstddev_s=" << s.stddev << "\tmad_s=" << s.mad
              << "\tkmers=" << kmers << "\trecords=" << records
              << "\tmkmer_s=" << mkmer_s << "\tstore=" << sizeof(store) << std::endl;
    return 0;
}

template<typename store>
int run_verify(const Args& a) {
    auto A = km::sortedlist::BucketedSkmerListReader<store>::open(a.a);
    auto B = km::sortedlist::BucketedSkmerListReader<store>::open(a.b);

    // Reference cardinality, computed via the CountSink path (set_sizes) — a different sink than the
    // CollectSink materialization, so a count==count match is a semi-independent corroboration.
    const uint64_t expected {setop_size<store>(a.op, A, B)};

    const uint64_t kmers {run_setop<store>(a.op, A, B, a.out)};
    const uint64_t records_opt {km::sortedlist::read_list_header(a.out).count};
    const uint64_t records_ref {km::sortedlist::read_list_header(a.ref).count};

    // Content equality against the frozen reference: both lists carry `expected` distinct k-mers and
    // share the bucket layout, so |O_opt ∩ O_ref| == expected  <=>  identical k-mer sets.
    auto Oopt = km::sortedlist::BucketedSkmerListReader<store>::open(a.out);
    auto Oref = km::sortedlist::BucketedSkmerListReader<store>::open(a.ref);
    const uint64_t inter {km::sortedlist::intersection_size<store>(Oopt, Oref, 1)};

    // bits/kmer = compaction quality (payload bits / k-mers). It must not jump up — a fast but
    // poorly-compacted result (towards one record per k-mer) would inflate it. Guaranteed not to grow
    // by the records_opt <= records_ref gate (bits/kmer is proportional to records here).
    const double rec_bits {static_cast<double>(sizeof(km::Skmer<store>)) * 8.0};
    const double bpk_opt {expected ? records_opt * rec_bits / static_cast<double>(expected) : 0.0};
    const double bpk_ref {expected ? records_ref * rec_bits / static_cast<double>(expected) : 0.0};

    const bool ok_count   {kmers == expected};
    const bool ok_set     {inter == expected};
    const bool ok_compact {records_opt <= records_ref};
    const bool pass {ok_count && ok_set && ok_compact};

    std::cerr << "[verify] expected=" << expected << "  kmers=" << kmers
              << "  inter(O_opt,O_ref)=" << inter
              << "  records_opt=" << records_opt << "  records_ref=" << records_ref
              << "  bits/kmer opt=" << bpk_opt << " ref=" << bpk_ref << "\n";
    if (!ok_count)   std::cerr << "[verify] FAIL: union k-mer count != union_size(A,B)\n";
    if (!ok_set)     std::cerr << "[verify] FAIL: k-mer set differs from reference\n";
    if (!ok_compact) std::cerr << "[verify] FAIL: compaction regressed (more records than reference)\n";
    std::cout << (pass ? "VERIFY\tPASS" : "VERIFY\tFAIL")
              << "\texpected=" << expected << "\tkmers=" << kmers << "\tinter=" << inter
              << "\trecords_opt=" << records_opt << "\trecords_ref=" << records_ref
              << "\tbits_per_kmer=" << bpk_opt << std::endl;
    return pass ? 0 : 1;
}

} // namespace

int main(int argc, char** argv) {
    const Args a {parse(argc, argv)};
    try {
        const uint64_t store_w {km::sortedlist::read_list_header(a.a).store_width_bytes};
        const uint64_t store_w_b {km::sortedlist::read_list_header(a.b).store_width_bytes};
        if (store_w != store_w_b)
            throw std::runtime_error("A and B store records at different widths; rebuild with same k/m/--buckets");
        return km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() -> int {
            return a.mode == "bench" ? run_bench<store>(a) : run_verify<store>(a);
        });
    } catch (const std::exception& e) {
        std::cerr << "union_bench error: " << e.what() << std::endl;
        return 2;
    }
}
