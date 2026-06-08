# sklib — compact sorted k-mer library

A C++23 library for sorting and querying k-mers in a compact representation. sklib packs k-mers and super-k-mers (skmers) as 2-bit-per-nucleotide integers, builds a disk-backed, **bucketed** sorted skmer list (`VirtualSkmer`), answers membership queries — each routed to a single minimizer bucket — and computes **set operations** (intersection / union / difference, and their cardinalities) between two lists, all through the `sskm` CLI.

## Documentation

A page-by-page walkthrough of the internals lives in the [project wiki](https://github.com/yoann-dufresne/sklib/wiki): the 2-bit nucleotide encoding, the minimizer-centered interleaved super-k-mer layout, super-k-mer generation (k-mers are canonicalized at yield so a k-mer is always stored under a single orientation), sorted-list construction (with the Fenwick-tree colinear chaining), the bucketed on-disk format (the list is split into minimizer-prefix sub-lists, with width-selectable and minimizer-prefix-quotiented records), the k-mer search algorithm (each query is routed to one bucket; file queries run in parallel), the **set operations** between two lists (a **parallel per-bucket merge**, with a faster patience-sort chaining on the re-compaction path), and **complete cross-tool benchmarks** (set operations vs KMC / CBL / FMSI; construction & membership queries vs sshash / SBWT / CBL / BQF / FMSI).

## Dependencies

* zlib : https://zlib.net/ (apt install zlib1g-dev)
* bzip2 : https://sourceware.org/bzip2/ (apt install libbz2-dev)

Third-party dependencies fetched automatically: `CLI11`, `zlib`, `gtest`, `kseqpp`. The core library is header-only; `sskm query` (file path), `sskm setop` (per-bucket merge) and `sskm construct` (per-bucket compaction) parallelize with the C++ standard library only (`std::thread`; the CLI links `pthread`), so there is no TBB / parallel-STL runtime dependency.

## Compilation

```bash
mkdir -p build && cd build
cmake .. && make
```

> **Compiler / GCC-13 users.** `sklib` is C++23 (`CMAKE_CXX_STANDARD 23`). The `dev` and
> `main` branches track the latest toolchain and may require a compiler newer than GCC-13.
> A dedicated **`gcc-13` branch** — tagged **`v0.7.0-gcc13`** — holds a GCC-13-compatible
> build, kept in sync with `dev` feature-for-feature, so there are two working versions in
> parallel. Check out `gcc-13` if `dev`/`main` does not build with your compiler.

## Tests

* Run all tests
```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make -j && ./tests/sklib-tests
```

* Run a subset of tests
```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make && ./tests/sklib-tests --gtest_filter="regexp test names"
```

For example, to run tests on the (virtual) super-k-mer sorting:
```bash
./tests/sklib-tests --gtest_filter=*SkmerSorting*
```

## CLI — `sskm`

The `sskm` tool (built to `build/bin/sskm`) exposes the library's primary operations: `construct` builds a sorted super-k-mer list from a FASTA input, `query` looks k-mers up against an existing list, and `setop` computes set operations between two lists. Run `sskm --help`, `sskm construct --help`, `sskm query --help`, or `sskm setop --help` for the in-CLI reference; `sskm --version` prints the current version.

### Parameters `k` and `m`

* `k` — k-mer length in nucleotides. The record integer width is selected automatically (`uint32_t`/`uint64_t`/`__uint128_t`) from the packed skmer size: the only limit is that a skmer's `2·(2k−m)` bits fit the widest 256-bit `__uint128_t` pair, i.e. `2·(2k−m) ≤ 256` (so `k` up to ~63 at small `m`, more as `m` grows); `construct` reports a clear error above that.
* `m` — minimizer length in nucleotides, with `1 ≤ m ≤ k`. Smaller `m` produces longer skmers (more shared central nucleotides per group); typical values sit around `m ≈ k/2`.

### `sskm construct`

Builds a sorted skmer list from a FASTA file.

```bash
./bin/sskm construct -k 21 -m 11 -f input.fa -o out.sskm
```

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-k, --kmer-size <int>` | yes | — | k-mer length; the record width is auto-selected, the only cap is `2·(2k−m) ≤ 256` (see above). |
| `-m, --minimizer-size <int>` | yes | — | Minimizer length (1 ≤ m ≤ k). |
| `-f, --file <path>` | no | stdin | Input FASTA (plain or gzip). |
| `-o, --output <path>` | no | stdout | Output sorted skmer list. |
| `--ascii` | no | off | Write human-readable ASCII instead of the default binary format. |
| `--buckets <int>` | no | 4096 | Number of minimizer-prefix buckets (rounded down to a power of two). The list is split into that many independently-sorted sub-lists, partitioned by the high-order bits of the minimizer; a query is routed to a single sub-list. More buckets ⇒ lower peak construction RAM and a faster query, and drop more *quotient* bits per record (the bucket id implies the top `log2(buckets)` minimizer bits), which can shrink each record to a narrower integer width — a smaller index. `1` is a single list. |
| `--max-ram <size>` | no | — | Approximate target peak RAM (e.g. `2G`, `512M`); a histogram pass derives adaptive buckets balanced to this budget, overriding `--buckets`. Requires `-f <file>` (the input is read twice). |
| `--tmp-dir <path>` | no | next to output | Directory for the temporary bucket files. |

The default binary output embeds an endianness marker so files can be moved across machines. The ASCII format is intended for inspection and tests, not for downstream querying.

#### Bucketed construction and on-disk layout

By default, binary construction to a regular file (`-o`) is **disk-backed and low-memory**: super-k-mers are partitioned by the high-order bits of their minimizer into on-disk buckets, then each bucket is sorted, deduplicated, and built independently. The buckets are kept as **separate sorted sub-lists** in the output (the binary `VSKMER_4` format stores a per-bucket directory — see the wiki's *On-disk Format and Serialization* page), so a query can be routed to one sub-list instead of scanning the whole list. Peak construction RAM is bounded by the **largest bucket** rather than by the whole input, tunable with `--buckets` (or targeted with `--max-ram`); construction is also *faster* than an all-in-RAM build, since the column/chaining algorithm runs on small per-bucket inputs. Measured single-threaded in a Release build, the bucketed path keeps construction peak RSS to tens–hundreds of MB even on large genomes — far below an all-in-RAM build (e.g. human chr1 `k=21, m=11`: ~198 MB peak RSS, ~55 s single-threaded Release) — see [`benchmark/scripts/README.md`](benchmark/scripts/README.md). Temporary bucket files are removed automatically on success and on error.

`--ascii` and writing to stdout (no `-o`) fall back to an all-in-RAM build that emits a single bucket (a stream cannot seek back to patch the header's per-bucket directory and total count).

#### Minimizer ordering (compact index)

The order on minimizers uses a fixed, invertible hash (φ) rather than the raw lexicographic value. Hash-ordered minimizers have lower density than lexicographic ones, so the list holds **9–21% fewer super-k-mers** (a smaller index / fewer bits per k-mer) across genomes, and on repeat-rich genomes it also sharply lowers peak construction RAM (e.g. *C. elegans* `k=21, m=11`: ~341 MB → ~63 MB) by breaking the over-selection of low-complexity minimizers. Each k-mer is stored in a strand-invariant canonical frame derived from its own nucleotides; super-k-mers whose minimizer is *ambiguous* (an RC-palindrome, or repeated within the window) are re-framed per k-mer — picking the minimizer occurrence by minimal φ-rank, then most central, then smaller full interleaved value — and re-merged into compact super-k-mers. This makes queries exact at any `m` (so there is no minimum-`m` restriction) and also removes a rare construction drop of poly-A k-mers. The on-disk minimizer slot is hash-permuted and the binary file carries a format version. The current format is **`VSKMER_4`**: each record is **minimizer-prefix quotiented** (the top `b = log2(buckets)` φ-minimizer bits are constant within a bucket, so the bucket id implies them and they are dropped) and stored in a **runtime-selected integer width** — the smallest of `uint32_t`/`uint64_t`/`__uint128_t` whose pair holds the remaining `2·(2k−m) − b` bits. Legacy bucketed **`VSKMER_3`** and single-list **`VSKMER_2`** files still load (as non-quotiented 64-bit lists); pre-φ raw-order **`VSKMER_M`** files are rejected on load rather than queried incorrectly. Rebuild rejected lists — and, since the v0.4.2 framing fix re-frames a few rare *ambiguous* k-mers, rebuild any index built before v0.4.2.

### `sskm query`

Queries k-mers against an existing list, either from a FASTA file or from a single sequence given as a positional argument.

```bash
./bin/sskm query -l out.sskm -i queries.fa -o hits.txt
# or pass a sequence directly:
./bin/sskm query -l out.sskm ACGTACGTACGTACGTACGTA
```

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-l, --list <path>` | yes | — | Sorted skmer list produced by `sskm construct`. |
| `-i, --input <path>` | one of `-i` / `sequence` | — | FASTA file whose k-mers are extracted and looked up. |
| `sequence` (positional) | one of `-i` / `sequence` | — | Single DNA sequence to query. |
| `-o, --output <path>` | no | stdout | Output file for query hits. |
| `-t, --threads <int>` | no | 8 | Worker threads for the file (`-i`) path: 1 producer reads and buckets the input, the rest query in parallel; output stays in input order. `1` = sequential. The positional `sequence` path is always sequential. |

`-i` and the positional `sequence` argument are mutually exclusive — exactly one must be provided.

Internally, each query super-k-mer is **routed to the single bucket** its minimizer falls in (the high-order bits pick the sub-list); only that bucket is loaded from disk (lazily, then cached) and searched, so a query reads a fraction of the index rather than the whole file. File queries (`-i`) are multithreaded by default (`-t`, default 8): one thread reads and buckets the input into input-order batches while the rest query the shared (thread-safe) reader in parallel — the output is byte-identical to a sequential run and roughly 3× faster on real workloads. See the wiki's *Querying a Sorted List*, *The Kmer Search Algorithm*, and *Batched, streaming & parallel queries* pages.

### `sskm setop`

Set operations between two sorted skmer lists **A** and **B**: `intersection`, `union`, and `diff` (asymmetric, `A \ B`), plus `intersection_size`, `union_size`, and `diff_size` variants that report only the result **cardinality**. The atomic element is a *k-mer*: the operation decomposes per minimizer bucket and per minimizer-position column, and a single two-cursor merge over the two sorted lists yields all three relations (`A∩B`, `A\B`, `B\A`) in one pass.

```bash
# materialize A ∩ B into a new list
./bin/sskm setop --op intersection -a A.sskm -b B.sskm -o inter.sskm
# just the cardinality (nothing written) — e.g. the numerator/denominator of a Jaccard index
./bin/sskm setop --op union_size -a A.sskm -b B.sskm
# fast materialization when a larger output is acceptable (skip super-k-mer re-compaction)
./bin/sskm setop --op union --no-compact -a A.sskm -b B.sskm -o union.sskm
# parallel by bucket (default 8 threads); the output is byte-identical for any -t
./bin/sskm setop --op intersection -a A.sskm -b B.sskm -o inter.sskm -t 16
```

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--op <name>` | yes | — | One of `intersection`, `union`, `diff`, `intersection_size`, `union_size`, `diff_size`. `diff` is asymmetric: `A \ B` (k-mers of A absent from B). |
| `-a, --list-a <path>` | yes | — | First sorted skmer list. |
| `-b, --list-b <path>` | yes | — | Second sorted skmer list. |
| `-o, --output <path>` | for `intersection`/`union`/`diff` | — | Output list for the result; ignored by the `*_size` variants, which print a count to stdout. |
| `--no-compact` | no | off | Emit **one record per result k-mer** instead of re-compacting the result back into super-k-mers — much faster (it skips the dominant cost) at the price of a larger output file (~3–5×). The result is still a valid, queryable sorted list. Ignored by the `*_size` variants. |
| `-t, --threads <N>` | no | 8 | Worker threads for the per-bucket merge — the independent minimizer buckets are processed in parallel. The output is **byte-identical** regardless of `N`. |

**Both lists must be built with the same `k`, `m`, and `--buckets`**: the merge aligns bucket *i* of A with bucket *i* of B, so the bucket layouts must match (this is checked, and a mismatch is rejected — rebuild both with identical parameters). The independent buckets are merged **in parallel** (`-t`); each worker holds only its current bucket pair, so peak RAM stays low (a handful of buckets) regardless of list size — far below KMC/CBL — and grows modestly with the thread count.

Two cost regimes:

* **Cardinality (`*_size`)** counts the three relations in one merge pass **without writing anything** — the fast path for a Jaccard index, containment, or de-duplication, and the right choice when the result *set* is not needed.
* **Materialized (`intersection`/`union`/`diff`)** additionally re-compacts the kept k-mers back into super-k-mers and writes a new list. That re-compaction dominates the time; `--no-compact` skips it (one record per k-mer, ~3–5× larger output, ~3× faster) while keeping the output a valid queryable list.

Set operations scale ~×6–9 with `-t` (near-linear to 8 cores), which makes sklib the fastest set-op tool from a few cores up — including the high-overlap materialized intersection where KMC led single-core (chr1 ∩: 2.6 s vs KMC 2.7 s at 22 cores). Two benchmark reports accompany this feature: [`benchmark/results/reports/SETOPS_REPORT.md`](benchmark/results/reports/SETOPS_REPORT.md) compares sklib against KMC, CBL and FMSI (single- and multi-core), and [`benchmark/results/reports/SETOPS_BOTTLENECKS.md`](benchmark/results/reports/SETOPS_BOTTLENECKS.md) breaks down where the time goes and documents the implemented speedups. See also the wiki's *Set operations* and *Benchmark · Set operations* pages.
