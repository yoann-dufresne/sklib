#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/SetOperations.hpp>
#include <algorithms/SortedSkmerListBuilder.hpp>
#include <algorithms/VirtualSkmer.hpp>
#include <algorithms/WidthDispatch.hpp>

// Tests for the set operations (intersection / union / diff and the *_size variants) on two sorted
// super-k-mer lists. The oracle is independent of the merge: each list is reduced to the set of its
// k-mer keys (column, masked interleaved value) with std::set, and the three relations are computed
// by std::set membership. We then check that (a) the in-RAM merge_columns emits exactly those keys,
// (b) the *_size counts match, and (c) the file drivers agree across bucket counts/widths, with full
// content checked at buckets=1 (no quotient truncation, so keys compare directly).

namespace {

using km::sortedlist::BucketedSkmerListReader;
using km::sortedlist::SortedVirtualSkmerList;

// A k-mer key: (column, high word, low word) of the masked interleaved value. Comparable for std::set
// across every record width (uint32/uint64/__uint128). Only meaningful at b=0 (no quotienting), where
// the full value is retained — which is exactly where we use it (in-RAM lists and buckets=1 files).
template<typename kuint>
using Key = std::tuple<uint64_t, kuint, kuint>;

std::string random_dna(uint64_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    static const char alphabet[] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(n);
    for (uint64_t i = 0; i < n; ++i) s.push_back(alphabet[rng() & 3u]);
    return s;
}

void write_fasta(const std::string& path, const std::vector<std::string>& seqs) {
    std::ofstream f(path);
    for (size_t i = 0; i < seqs.size(); ++i) f << ">s" << i << "\n" << seqs[i] << "\n";
}

// Enumerate the (canonical, φ-permuted) super-k-mers of a set of sequences, exactly as construction
// does (one manipulator per sequence; k-mers never cross sequence boundaries).
template<typename kuint>
std::vector<km::Skmer<kuint>> enumerate(uint64_t k, uint64_t m, const std::vector<std::string>& seqs) {
    std::vector<km::Skmer<kuint>> out;
    for (const std::string& seq_in : seqs) {
        if (seq_in.size() < k) continue;
        std::string seq = seq_in;   // SeqSkmerator binds a non-const lvalue reference
        km::SkmerManipulator<kuint> manip{k, m};
        km::SeqSkmerator<kuint> rator{manip, seq};
        for (const km::Skmer<kuint>& s : rator) out.push_back(s);
    }
    return out;
}

// Set of k-mer keys carried by a span of records (one key per record × valid column).
template<typename kuint>
std::set<Key<kuint>> keys_of(km::SkmerManipulator<kuint>& manip, const km::Skmer<kuint>* L, size_t n) {
    const std::vector<typename km::Skmer<kuint>::pair> kmasks = manip.get_k_mask();
    const uint64_t km = manip.k - manip.m;
    std::set<Key<kuint>> ks;
    for (size_t i = 0; i < n; ++i)
        for (uint64_t c = 0; c <= km; ++c)
            if (manip.has_valid_kmer(L[i], c)) {
                const typename km::Skmer<kuint>::pair v = L[i].m_pair & kmasks[c];
                ks.insert(Key<kuint>{c, v.m_value[1], v.m_value[0]});
            }
    return ks;
}

// Sink that records each emitted k-mer as a key, per relation (deduping via std::set).
template<typename kuint>
struct CaptureSink {
    std::vector<typename km::Skmer<kuint>::pair> kmasks;
    std::set<Key<kuint>> inter, onlyA, onlyB;
    Key<kuint> mk(const km::Skmer<kuint>& r, uint64_t c) {
        const typename km::Skmer<kuint>::pair v = r.m_pair & kmasks[c];
        return Key<kuint>{c, v.m_value[1], v.m_value[0]};
    }
    void both  (const km::Skmer<kuint>& r, uint64_t c) { inter.insert(mk(r, c)); }
    void only_a(const km::Skmer<kuint>& r, uint64_t c) { onlyA.insert(mk(r, c)); }
    void only_b(const km::Skmer<kuint>& r, uint64_t c) { onlyB.insert(mk(r, c)); }
};

// Build a list to disk at a given bucket count; return its path.
std::string build_list_file(uint64_t k, uint64_t m, const std::vector<std::string>& seqs,
                            uint64_t buckets, const std::string& tag) {
    const std::string in  = ::testing::TempDir() + "setop_in_"  + tag + ".fa";
    const std::string out = ::testing::TempDir() + "setop_out_" + tag + ".sskm";
    write_fasta(in, seqs);
    km::sortedlist::SortedListBuildParams p;
    p.k = k; p.m = m; p.input_path = in; p.output_path = out;
    p.ascii = false; p.buckets = buckets; p.has_output_file = true;
    const uint64_t b = km::sortedlist::quotient_bits_for(p);
    const uint64_t gen_w   = km::sortedlist::select_width_bytes(2 * (2 * k - m));
    const uint64_t store_w = km::sortedlist::select_width_bytes(2 * (2 * k - m) - b);
    km::sortedlist::dispatch_width_bytes(gen_w, [&]<typename gen>() {
        km::sortedlist::dispatch_width_bytes(store_w, [&]<typename store>() {
            km::sortedlist::build_sorted_list<gen, store>(p, b);
        });
    });
    return out;
}

// Run set_sizes over two list files, dispatching on their (shared) store width.
km::sortedlist::SetSizes file_sizes(const std::string& a, const std::string& b) {
    const uint64_t w = km::sortedlist::read_list_header(a).store_width_bytes;
    km::sortedlist::SetSizes out;
    km::sortedlist::dispatch_width_bytes(w, [&]<typename store>() {
        auto A = BucketedSkmerListReader<store>::open(a);
        auto B = BucketedSkmerListReader<store>::open(b);
        out = km::sortedlist::set_sizes<store>(A, B);
    });
    return out;
}

// Materialize one operation ("intersection"|"union"|"diff"|"xor") to `out`; return result k-mer count.
uint64_t file_materialize(const std::string& op, const std::string& a, const std::string& b,
                          const std::string& out) {
    const uint64_t w = km::sortedlist::read_list_header(a).store_width_bytes;
    uint64_t n = 0;
    km::sortedlist::dispatch_width_bytes(w, [&]<typename store>() {
        auto A = BucketedSkmerListReader<store>::open(a);
        auto B = BucketedSkmerListReader<store>::open(b);
        if (op == "intersection") n = km::sortedlist::intersection<store>(A, B, out);
        else if (op == "union")   n = km::sortedlist::set_union<store>(A, B, out);
        else if (op == "diff")    n = km::sortedlist::difference<store>(A, B, out);
        else                      n = km::sortedlist::symmetric_difference<store>(A, B, out);
    });
    return n;
}

// Read a whole file into a string for byte-for-byte comparison of two output lists.
std::string slurp(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::string data;
    f.seekg(0, std::ios::end);
    const std::streamoff len = f.tellg();
    if (len > 0) { data.resize(static_cast<size_t>(len)); f.seekg(0); f.read(&data[0], len); }
    return data;
}

// Materialize one op with an explicit thread count / no_compact flag (the parallel set-op path).
uint64_t file_materialize_threads(const std::string& op, const std::string& a, const std::string& b,
                                  const std::string& out, unsigned nthreads, bool no_compact) {
    const uint64_t w = km::sortedlist::read_list_header(a).store_width_bytes;
    uint64_t n = 0;
    km::sortedlist::dispatch_width_bytes(w, [&]<typename store>() {
        auto A = BucketedSkmerListReader<store>::open(a);
        auto B = BucketedSkmerListReader<store>::open(b);
        if (op == "intersection") n = km::sortedlist::intersection<store>(A, B, out, no_compact, nthreads);
        else if (op == "union")   n = km::sortedlist::set_union<store>(A, B, out, no_compact, nthreads);
        else if (op == "diff")    n = km::sortedlist::difference<store>(A, B, out, no_compact, nthreads);
        else                      n = km::sortedlist::symmetric_difference<store>(A, B, out, no_compact, nthreads);
    });
    return n;
}

// set_sizes with an explicit thread count.
km::sortedlist::SetSizes file_sizes_threads(const std::string& a, const std::string& b, unsigned nthreads) {
    const uint64_t w = km::sortedlist::read_list_header(a).store_width_bytes;
    km::sortedlist::SetSizes out;
    km::sortedlist::dispatch_width_bytes(w, [&]<typename store>() {
        auto A = BucketedSkmerListReader<store>::open(a);
        auto B = BucketedSkmerListReader<store>::open(b);
        out = km::sortedlist::set_sizes<store>(A, B, nthreads);
    });
    return out;
}

// Run multi_setop over two list files (dispatching on their shared store width), materializing
// whichever of the four outputs are given (std::nullopt = not requested). Returns the full result.
km::sortedlist::MultiSetOpResult file_multi_setop(
        const std::string& a, const std::string& b,
        const std::optional<std::string>& inter, const std::optional<std::string>& uni,
        const std::optional<std::string>& diff_ab, const std::optional<std::string>& diff_ba,
        const std::optional<std::string>& sym_diff,
        unsigned nthreads, bool no_compact) {
    const uint64_t w = km::sortedlist::read_list_header(a).store_width_bytes;
    km::sortedlist::MultiSetOpResult res;
    km::sortedlist::dispatch_width_bytes(w, [&]<typename store>() {
        auto A = BucketedSkmerListReader<store>::open(a);
        auto B = BucketedSkmerListReader<store>::open(b);
        km::sortedlist::MultiSetOpRequest req;
        req.inter_out = inter; req.union_out = uni; req.diff_ab_out = diff_ab; req.diff_ba_out = diff_ba;
        req.sym_diff_out = sym_diff;
        req.no_compact = no_compact;
        res = km::sortedlist::multi_setop<store>(A, B, req, nthreads);
    });
    return res;
}

// The combined single-pass multi_setop must produce, for every bucket count, byte-for-byte the same
// files as the four individual single-op materializations (intersection / union / A\B / B\A), and the
// same four cardinalities and sizes. This is the headline guarantee that "one pass" == "four passes".
void expect_multi_matches_individual(uint64_t k, uint64_t m,
                                     const std::vector<std::string>& seqsA,
                                     const std::vector<std::string>& seqsB,
                                     const std::string& tag) {
    for (uint64_t buckets : {uint64_t{1}, uint64_t{4}, uint64_t{16}, uint64_t{64}}) {
        const std::string ta = tag + "_b" + std::to_string(buckets);
        const std::string pa = build_list_file(k, m, seqsA, buckets, ta + "_a");
        const std::string pb = build_list_file(k, m, seqsB, buckets, ta + "_b");
        const std::string td = ::testing::TempDir();

        // individual references (single-op, -t 1). B\A is diff with the inputs swapped.
        const std::string ri = td + "ref_i_" + ta, ru = td + "ref_u_" + ta;
        const std::string rdab = td + "ref_dab_" + ta, rdba = td + "ref_dba_" + ta, rsym = td + "ref_sym_" + ta;
        const uint64_t ni   = file_materialize_threads("intersection", pa, pb, ri, 1, false);
        const uint64_t nu   = file_materialize_threads("union",        pa, pb, ru, 1, false);
        const uint64_t ndab = file_materialize_threads("diff",         pa, pb, rdab, 1, false);
        const uint64_t ndba = file_materialize_threads("diff",         pb, pa, rdba, 1, false);
        const uint64_t nsym = file_materialize_threads("xor",          pa, pb, rsym, 1, false);

        // combined (single pass, all five outputs), -t 1
        const std::string ci = td + "cmb_i_" + ta, cu = td + "cmb_u_" + ta;
        const std::string cdab = td + "cmb_dab_" + ta, cdba = td + "cmb_dba_" + ta, csym = td + "cmb_sym_" + ta;
        const km::sortedlist::MultiSetOpResult r =
            file_multi_setop(pa, pb, ci, cu, cdab, cdba, csym, 1, false);

        EXPECT_EQ(slurp(ci),   slurp(ri))   << tag << " buckets=" << buckets << ": intersection bytes";
        EXPECT_EQ(slurp(cu),   slurp(ru))   << tag << " buckets=" << buckets << ": union bytes";
        EXPECT_EQ(slurp(cdab), slurp(rdab)) << tag << " buckets=" << buckets << ": A\\B bytes";
        EXPECT_EQ(slurp(cdba), slurp(rdba)) << tag << " buckets=" << buckets << ": B\\A bytes";
        EXPECT_EQ(slurp(csym), slurp(rsym)) << tag << " buckets=" << buckets << ": A△B bytes";

        EXPECT_EQ(r.inter_kmers,    ni)   << tag << " buckets=" << buckets;
        EXPECT_EQ(r.union_kmers,    nu)   << tag << " buckets=" << buckets;
        EXPECT_EQ(r.diff_ab_kmers,  ndab) << tag << " buckets=" << buckets;
        EXPECT_EQ(r.diff_ba_kmers,  ndba) << tag << " buckets=" << buckets;
        EXPECT_EQ(r.sym_diff_kmers, nsym) << tag << " buckets=" << buckets;

        const km::sortedlist::SetSizes s = file_sizes(pa, pb);
        EXPECT_EQ(r.sizes.inter,  s.inter)  << tag << " buckets=" << buckets;
        EXPECT_EQ(r.sizes.only_a, s.only_a) << tag << " buckets=" << buckets;
        EXPECT_EQ(r.sizes.only_b, s.only_b) << tag << " buckets=" << buckets;
    }
}

// The whole suite at a fixed record width `store` (chosen so the b=0 list fits it).
template<typename store>
void run_suite(uint64_t k, uint64_t m,
               const std::vector<std::string>& seqsA, const std::vector<std::string>& seqsB,
               const std::string& tag) {
    // ---- reference: keys from in-RAM (b=0) lists, set ops via std::set ----
    // Build in place (SortedVirtualSkmerList is not safely copy/movable — SkmerManipulator owns raw
    // mask arrays), so never construct via a by-value return.
    const std::vector<km::Skmer<store>> enA = enumerate<store>(k, m, seqsA);
    const std::vector<km::Skmer<store>> enB = enumerate<store>(k, m, seqsB);
    SortedVirtualSkmerList<store> A{k, m};
    A.generate_sorted_list_from_enumeration(enA);
    SortedVirtualSkmerList<store> B{k, m};
    B.generate_sorted_list_from_enumeration(enB);
    km::SkmerManipulator<store> manip{k, m};
    std::set<Key<store>> KA = keys_of<store>(manip, A.get_list().data(), A.size());
    std::set<Key<store>> KB = keys_of<store>(manip, B.get_list().data(), B.size());

    std::set<Key<store>> R_inter, R_onlyA, R_onlyB;
    for (const Key<store>& kk : KA) (KB.count(kk) ? R_inter : R_onlyA).insert(kk);
    for (const Key<store>& kk : KB) if (!KA.count(kk)) R_onlyB.insert(kk);
    const uint64_t exp_inter = R_inter.size();
    const uint64_t exp_onlyA = R_onlyA.size();
    const uint64_t exp_union = R_inter.size() + R_onlyA.size() + R_onlyB.size();
    const uint64_t exp_symdiff = R_onlyA.size() + R_onlyB.size();   // |A △ B|

    // ---- (a) in-RAM merge_columns emits exactly the reference keys ----
    {
        CaptureSink<store> sink{manip.get_k_mask(), {}, {}, {}};
        km::sortedlist::merge_columns<store>(manip, A.get_list().data(), A.size(),
                                             B.get_list().data(), B.size(), sink);
        EXPECT_EQ(sink.inter, R_inter) << tag << ": intersection keys";
        EXPECT_EQ(sink.onlyA, R_onlyA) << tag << ": A\\B keys";
        EXPECT_EQ(sink.onlyB, R_onlyB) << tag << ": B\\A keys";
    }
    // ---- (b) CountSink sizes match, and equal op-emission counts (so no duplicates) ----
    {
        km::sortedlist::CountSink<store> cs;
        km::sortedlist::merge_columns<store>(manip, A.get_list().data(), A.size(),
                                             B.get_list().data(), B.size(), cs);
        EXPECT_EQ(cs.n_inter, exp_inter) << tag;
        EXPECT_EQ(cs.n_only_a, exp_onlyA) << tag;
        EXPECT_EQ(cs.n_inter + cs.n_only_a + cs.n_only_b, exp_union) << tag;
        EXPECT_EQ(cs.n_only_a + cs.n_only_b, exp_symdiff) << tag;
    }

    // ---- (c) file drivers across bucket counts ----
    for (uint64_t buckets : {uint64_t{1}, uint64_t{4}, uint64_t{16}}) {
        const std::string ta = tag + "_a_b" + std::to_string(buckets);
        const std::string tb = tag + "_b_b" + std::to_string(buckets);
        const std::string pa = build_list_file(k, m, seqsA, buckets, ta);
        const std::string pb = build_list_file(k, m, seqsB, buckets, tb);

        // sizes
        const km::sortedlist::SetSizes s = file_sizes(pa, pb);
        EXPECT_EQ(s.inter, exp_inter) << tag << " buckets=" << buckets;
        EXPECT_EQ(s.only_a, exp_onlyA) << tag << " buckets=" << buckets;
        EXPECT_EQ(s.uni(), exp_union) << tag << " buckets=" << buckets;
        EXPECT_EQ(s.sym_diff(), exp_symdiff) << tag << " buckets=" << buckets;
        // op vs _size: materialized result count equals the corresponding _size
        const std::string pi = ::testing::TempDir() + "setop_res_i_" + ta + ".sskm";
        const std::string pu = ::testing::TempDir() + "setop_res_u_" + ta + ".sskm";
        const std::string pd = ::testing::TempDir() + "setop_res_d_" + ta + ".sskm";
        const std::string px = ::testing::TempDir() + "setop_res_x_" + ta + ".sskm";
        EXPECT_EQ(file_materialize("intersection", pa, pb, pi), exp_inter) << tag << " buckets=" << buckets;
        EXPECT_EQ(file_materialize("union", pa, pb, pu), exp_union) << tag << " buckets=" << buckets;
        EXPECT_EQ(file_materialize("diff", pa, pb, pd), exp_onlyA) << tag << " buckets=" << buckets;
        EXPECT_EQ(file_materialize("xor", pa, pb, px), exp_symdiff) << tag << " buckets=" << buckets;

        // ---- full content check at buckets=1 (b=0, full keys comparable to the reference) ----
        if (buckets == 1) {
            auto Ri = BucketedSkmerListReader<store>::open(pi);
            auto Ru = BucketedSkmerListReader<store>::open(pu);
            auto Rd = BucketedSkmerListReader<store>::open(pd);
            auto Rx = BucketedSkmerListReader<store>::open(px);
            km::SkmerManipulator<store> rmanip{k, m};
            const std::vector<km::Skmer<store>>& bi = Ri.load_bucket(0);
            const std::vector<km::Skmer<store>>& bu = Ru.load_bucket(0);
            const std::vector<km::Skmer<store>>& bd = Rd.load_bucket(0);
            const std::vector<km::Skmer<store>>& bx = Rx.load_bucket(0);
            std::set<Key<store>> got_inter = keys_of<store>(rmanip, bi.data(), bi.size());
            std::set<Key<store>> got_union = keys_of<store>(rmanip, bu.data(), bu.size());
            std::set<Key<store>> got_diff  = keys_of<store>(rmanip, bd.data(), bd.size());
            std::set<Key<store>> got_xor   = keys_of<store>(rmanip, bx.data(), bx.size());
            std::set<Key<store>> R_union = R_inter;
            R_union.insert(R_onlyA.begin(), R_onlyA.end());
            R_union.insert(R_onlyB.begin(), R_onlyB.end());
            std::set<Key<store>> R_symdiff = R_onlyA;
            R_symdiff.insert(R_onlyB.begin(), R_onlyB.end());
            EXPECT_EQ(got_inter, R_inter)   << tag << ": materialized intersection content";
            EXPECT_EQ(got_union, R_union)   << tag << ": materialized union content";
            EXPECT_EQ(got_diff,  R_onlyA)   << tag << ": materialized diff content";
            EXPECT_EQ(got_xor,   R_symdiff) << tag << ": materialized xor content";
        }
    }
}

// Pick the b=0 record width for (k,m) and run the suite at that width.
void run_dispatched(uint64_t k, uint64_t m,
                    const std::vector<std::string>& A, const std::vector<std::string>& B,
                    const std::string& tag) {
    const uint64_t w = km::sortedlist::select_width_bytes(2 * (2 * k - m));
    km::sortedlist::dispatch_width_bytes(w, [&]<typename store>() {
        run_suite<store>(k, m, A, B, tag);
    });
}

} // namespace

// ---- scenarios ----

TEST(SetOperations, EmptyVsEmpty) {
    run_dispatched(15, 5, {}, {}, "empty_empty");
}

TEST(SetOperations, EmptyVsNonEmpty) {
    const std::string x = random_dna(300, 1);
    run_dispatched(15, 5, {}, {x}, "empty_x");
    run_dispatched(15, 5, {x}, {}, "x_empty");
}

TEST(SetOperations, Identical) {
    const std::string x = random_dna(400, 2);
    run_dispatched(21, 11, {x}, {x}, "identical");
}

TEST(SetOperations, Disjoint) {
    // Long independent random sequences: sharing a k-mer by chance is negligible at k=25.
    run_dispatched(25, 11, {random_dna(500, 3)}, {random_dna(500, 9999)}, "disjoint");
}

TEST(SetOperations, PartialOverlap) {
    const std::string shared = random_dna(300, 10);
    const std::string onlyA  = random_dna(250, 11);
    const std::string onlyB  = random_dna(250, 12);
    run_dispatched(21, 7, {shared, onlyA}, {shared, onlyB}, "overlap");
}

TEST(SetOperations, Inclusion) {
    // B ⊂ A: B's sequences are a subset of A's, so A\B leftover, B\A empty, A∩B == B.
    const std::string s1 = random_dna(300, 20);
    const std::string s2 = random_dna(300, 21);
    run_dispatched(21, 11, {s1, s2}, {s1}, "inclusion");
}

TEST(SetOperations, DuplicatesInSource) {
    // The same sequence repeated must collapse to the same k-mer set.
    const std::string x = random_dna(200, 30);
    run_dispatched(15, 5, {x, x, x}, {x}, "dups");
}

TEST(SetOperations, SmallMinimizerAmbiguous) {
    // m ≤ 5 exercises the ambiguous-minimizer canonical framing (canonical_pieces).
    const std::string shared = random_dna(400, 40);
    run_dispatched(17, 3, {shared, random_dna(300, 41)}, {shared, random_dna(300, 42)}, "m3");
    run_dispatched(15, 5, {shared, random_dna(300, 43)}, {shared, random_dna(300, 44)}, "m5");
}

TEST(SetOperations, WideRecordsUint128) {
    // 2*(2k-m) = 138 > 128 selects the __uint128_t record width.
    const std::string shared = random_dna(600, 50);
    run_dispatched(40, 11, {shared, random_dna(400, 51)}, {shared, random_dna(400, 52)}, "wide");
}

// Lists built with mismatched parameters must be rejected, not silently mis-merged.
TEST(SetOperations, IncompatibleParametersThrow) {
    const std::string x = random_dna(300, 60);
    // different bucket counts (same k/m/width)
    const std::string a4  = build_list_file(15, 5, {x}, 4,  "inc_a4");
    const std::string b16 = build_list_file(15, 5, {x}, 16, "inc_b16");
    EXPECT_THROW(file_sizes(a4, b16), std::runtime_error);

    // different k (same width and bucket count)
    const std::string ak15 = build_list_file(15, 5, {x}, 4, "inc_k15");
    const std::string bk17 = build_list_file(17, 5, {x}, 4, "inc_k17");
    EXPECT_THROW(file_sizes(ak15, bk17), std::runtime_error);
}

// The parallel per-bucket path must produce byte-for-byte the same output file as the sequential
// path, for every operation, compaction mode, and bucket count — and the *_size counts must match
// too. This is the load-bearing guarantee that bucket parallelism did not change results. Several
// bucket counts (incl. 64 with up to 8 workers) exercise real out-of-order completion + the ordered
// writer's reorder buffer.
TEST(SetOperations, ParallelMatchesSequential) {
    const uint64_t k = 21, m = 11;
    const std::string shared = random_dna(4000, 700);
    const std::vector<std::string> seqsA = {shared, random_dna(4000, 701), random_dna(2000, 702)};
    const std::vector<std::string> seqsB = {shared, random_dna(4000, 703), random_dna(2000, 704)};

    for (uint64_t buckets : {uint64_t{4}, uint64_t{16}, uint64_t{64}}) {
        const std::string tag = "par_b" + std::to_string(buckets);
        const std::string pa = build_list_file(k, m, seqsA, buckets, tag + "_a");
        const std::string pb = build_list_file(k, m, seqsB, buckets, tag + "_b");

        for (const std::string& op : {std::string("intersection"), std::string("union"), std::string("diff")}) {
            for (bool nc : {false, true}) {
                const std::string base = ::testing::TempDir() + "setop_par_" + tag + "_" + op + (nc ? "_nc" : "");
                const std::string seq = base + "_t1.sskm";
                const uint64_t n1 = file_materialize_threads(op, pa, pb, seq, 1, nc);
                const std::string ref_bytes = slurp(seq);
                EXPECT_FALSE(ref_bytes.empty()) << op << " produced no output";
                for (unsigned t : {2u, 4u, 8u}) {
                    const std::string par = base + "_t" + std::to_string(t) + ".sskm";
                    const uint64_t nt = file_materialize_threads(op, pa, pb, par, t, nc);
                    EXPECT_EQ(n1, nt) << op << " nc=" << nc << " buckets=" << buckets << " t=" << t;
                    EXPECT_EQ(ref_bytes, slurp(par))
                        << op << " nc=" << nc << " buckets=" << buckets << " t=" << t << ": output not byte-identical";
                }
            }
        }

        // _size must be thread-count-independent too.
        const km::sortedlist::SetSizes s1 = file_sizes_threads(pa, pb, 1);
        for (unsigned t : {2u, 4u, 8u}) {
            const km::sortedlist::SetSizes st = file_sizes_threads(pa, pb, t);
            EXPECT_EQ(s1.inter, st.inter)   << "buckets=" << buckets << " t=" << t;
            EXPECT_EQ(s1.only_a, st.only_a) << "buckets=" << buckets << " t=" << t;
            EXPECT_EQ(s1.only_b, st.only_b) << "buckets=" << buckets << " t=" << t;
        }
    }
}

// ---- combined (single-pass) multi_setop ----

// The combined pass equals the four individual single-op materializations, byte-for-byte, across
// bucket counts and record widths, for a range of overlap regimes.
TEST(SetOperations, MultiMatchesIndividual_PartialOverlap) {
    const std::string shared = random_dna(600, 1000);
    expect_multi_matches_individual(21, 11,
        {shared, random_dna(500, 1001)}, {shared, random_dna(500, 1002)}, "multi_overlap");
}
TEST(SetOperations, MultiMatchesIndividual_Inclusion) {
    const std::string s1 = random_dna(400, 1010), s2 = random_dna(400, 1011);
    expect_multi_matches_individual(21, 11, {s1, s2}, {s1}, "multi_incl");      // B ⊂ A
}
TEST(SetOperations, MultiMatchesIndividual_Disjoint) {
    expect_multi_matches_individual(25, 11, {random_dna(500, 1020)}, {random_dna(500, 1021)}, "multi_disjoint");
}
TEST(SetOperations, MultiMatchesIndividual_EmptyAndIdentical) {
    const std::string x = random_dna(400, 1030);
    expect_multi_matches_individual(15, 5, {}, {}, "multi_empty_empty");
    expect_multi_matches_individual(15, 5, {}, {x}, "multi_empty_b");
    expect_multi_matches_individual(21, 11, {x}, {x}, "multi_identical");        // A == B
}
TEST(SetOperations, MultiMatchesIndividual_SmallMinimizerAmbiguous) {
    const std::string shared = random_dna(500, 1035);
    expect_multi_matches_individual(17, 3, {shared, random_dna(300, 1036)},
                                          {shared, random_dna(300, 1037)}, "multi_m3");
}
TEST(SetOperations, MultiMatchesIndividual_WideUint128) {
    const std::string shared = random_dna(700, 1040);   // 2*(2k-m)=138 bits -> __uint128_t records
    expect_multi_matches_individual(40, 11, {shared, random_dna(400, 1041)},
                                          {shared, random_dna(400, 1042)}, "multi_wide");
}

// Each combined output file is byte-identical for any thread count (the V2 guarantee), in both
// compaction modes. Several bucket counts exercise real out-of-order completion + the shared
// multi-output ordered writer's reorder buffer.
TEST(SetOperations, MultiParallelByteIdentical) {
    const uint64_t k = 21, m = 11;
    const std::string shared = random_dna(4000, 800);
    const std::vector<std::string> seqsA = {shared, random_dna(4000, 801), random_dna(2000, 802)};
    const std::vector<std::string> seqsB = {shared, random_dna(4000, 803), random_dna(2000, 804)};

    for (uint64_t buckets : {uint64_t{16}, uint64_t{64}}) {
        const std::string ta = "multipar_b" + std::to_string(buckets);
        const std::string pa = build_list_file(k, m, seqsA, buckets, ta + "_a");
        const std::string pb = build_list_file(k, m, seqsB, buckets, ta + "_b");
        for (bool nc : {false, true}) {
            const std::string base = ::testing::TempDir() + ta + (nc ? "_nc" : "");
            // reference at t = 1
            file_multi_setop(pa, pb, base + "_i1", base + "_u1", base + "_dab1", base + "_dba1", base + "_x1", 1, nc);
            const std::string ri = slurp(base + "_i1"),   ru   = slurp(base + "_u1");
            const std::string rdab = slurp(base + "_dab1"), rdba = slurp(base + "_dba1"), rx = slurp(base + "_x1");
            for (unsigned t : {2u, 4u, 8u}) {
                const std::string s = "_t" + std::to_string(t);
                file_multi_setop(pa, pb, base + "_i" + s, base + "_u" + s,
                                 base + "_dab" + s, base + "_dba" + s, base + "_x" + s, t, nc);
                EXPECT_EQ(ri,   slurp(base + "_i" + s))   << "buckets=" << buckets << " nc=" << nc << " t=" << t << " inter";
                EXPECT_EQ(ru,   slurp(base + "_u" + s))   << "buckets=" << buckets << " nc=" << nc << " t=" << t << " union";
                EXPECT_EQ(rdab, slurp(base + "_dab" + s)) << "buckets=" << buckets << " nc=" << nc << " t=" << t << " A\\B";
                EXPECT_EQ(rdba, slurp(base + "_dba" + s)) << "buckets=" << buckets << " nc=" << nc << " t=" << t << " B\\A";
                EXPECT_EQ(rx,   slurp(base + "_x" + s))   << "buckets=" << buckets << " nc=" << nc << " t=" << t << " A△B";
            }
        }
    }
}

// Requesting only a subset of outputs yields exactly the corresponding single-op files (proves the
// null-target routing in MultiCollectSink keeps the relations correctly separated).
TEST(SetOperations, MultiSubsetsMatchIndividual) {
    const uint64_t k = 21, m = 7, buckets = 16;
    const std::string shared = random_dna(800, 900);
    const std::string pa = build_list_file(k, m, {shared, random_dna(500, 901)}, buckets, "sub_a");
    const std::string pb = build_list_file(k, m, {shared, random_dna(500, 902)}, buckets, "sub_b");
    const std::string td = ::testing::TempDir();

    const std::string ref_i = td + "sub_ref_i", ref_u = td + "sub_ref_u", ref_dab = td + "sub_ref_dab",
                      ref_x = td + "sub_ref_x";
    file_materialize_threads("intersection", pa, pb, ref_i, 1, false);
    file_materialize_threads("union",        pa, pb, ref_u, 1, false);
    file_materialize_threads("diff",         pa, pb, ref_dab, 1, false);
    file_materialize_threads("xor",          pa, pb, ref_x, 1, false);

    // {inter} only
    const std::string only_i = td + "sub_only_i";
    file_multi_setop(pa, pb, only_i, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 1, false);
    EXPECT_EQ(slurp(only_i), slurp(ref_i)) << "inter-only must equal single-op intersection";

    // {union} only
    const std::string only_u = td + "sub_only_u";
    file_multi_setop(pa, pb, std::nullopt, only_u, std::nullopt, std::nullopt, std::nullopt, 1, false);
    EXPECT_EQ(slurp(only_u), slurp(ref_u)) << "union-only must equal single-op union";

    // {xor} only
    const std::string only_x = td + "sub_only_x";
    file_multi_setop(pa, pb, std::nullopt, std::nullopt, std::nullopt, std::nullopt, only_x, 1, false);
    EXPECT_EQ(slurp(only_x), slurp(ref_x)) << "xor-only must equal single-op symmetric difference";

    // {inter, diff_ab} together (no union, no diff_ba, no xor)
    const std::string c_i = td + "sub_pair_i", c_dab = td + "sub_pair_dab";
    const km::sortedlist::MultiSetOpResult r =
        file_multi_setop(pa, pb, c_i, std::nullopt, c_dab, std::nullopt, std::nullopt, 4, false);
    EXPECT_EQ(slurp(c_i),   slurp(ref_i))   << "pair: intersection";
    EXPECT_EQ(slurp(c_dab), slurp(ref_dab)) << "pair: A\\B";
    // sizes still come for free even when union / diff_ba / xor are not materialized
    const km::sortedlist::SetSizes s = file_sizes(pa, pb);
    EXPECT_EQ(r.union_kmers,    s.uni());
    EXPECT_EQ(r.diff_ba_kmers,  s.only_b);
    EXPECT_EQ(r.sym_diff_kmers, s.sym_diff());
}

// With no output requested the call is a pure combined count: sizes match the set_sizes oracle and
// no list file is created.
TEST(SetOperations, MultiSizesOnly) {
    const uint64_t k = 21, m = 11, buckets = 16;
    const std::string shared = random_dna(800, 910);
    const std::string pa = build_list_file(k, m, {shared, random_dna(400, 911)}, buckets, "szonly_a");
    const std::string pb = build_list_file(k, m, {shared, random_dna(400, 912)}, buckets, "szonly_b");
    const km::sortedlist::SetSizes s = file_sizes(pa, pb);
    const km::sortedlist::MultiSetOpResult r =
        file_multi_setop(pa, pb, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 4, false);
    EXPECT_EQ(r.sizes.inter,  s.inter);
    EXPECT_EQ(r.sizes.only_a, s.only_a);
    EXPECT_EQ(r.sizes.only_b, s.only_b);
    EXPECT_EQ(r.inter_kmers,    s.inter);
    EXPECT_EQ(r.union_kmers,    s.uni());
    EXPECT_EQ(r.diff_ab_kmers,  s.only_a);
    EXPECT_EQ(r.diff_ba_kmers,  s.only_b);
    EXPECT_EQ(r.sym_diff_kmers, s.sym_diff());
}
