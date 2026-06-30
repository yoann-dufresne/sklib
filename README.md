# sklib — compact sorted k-mer library

A C++23 library for sorting and querying k-mers in a compact representation. sklib packs k-mers and super-k-mers (skmers) as 2-bit-per-nucleotide integers, builds a disk-backed, **bucketed** sorted skmer list (`VirtualSkmer`), answers membership queries — each routed to a single minimizer bucket — and computes **set operations** (intersection / union / difference / symmetric difference, and their cardinalities) between two lists, all through the `sskm` CLI.

## Documentation

A page-by-page walkthrough of the internals lives in the [project wiki](https://github.com/yoann-dufresne/sklib/wiki): the 2-bit nucleotide encoding, the minimizer-centered interleaved super-k-mer layout, super-k-mer generation (k-mers are canonicalized at yield so a k-mer is always stored under a single orientation), sorted-list construction (with the Fenwick-tree colinear chaining), the bucketed on-disk format (the list is split into minimizer-prefix sub-lists, with width-selectable and minimizer-prefix-quotiented records), the k-mer search algorithm (each query is routed to one bucket; file queries run in parallel), the **set operations** between two lists (a **parallel per-bucket merge**, with a faster patience-sort chaining on the re-compaction path), and **complete cross-tool benchmarks** (set operations vs KMC / CBL / FMSI; construction & membership queries vs sshash / SBWT (C++ & Rust) / CBL / BQF / FMSI).

## Dependencies

* zlib : https://zlib.net/ (apt install zlib1g-dev)
* bzip2 : https://sourceware.org/bzip2/ (apt install libbz2-dev)

Third-party dependencies fetched automatically: `CLI11`, `zlib`, `gtest`, `kseqpp`. The core library is header-only; `sskm query` (file path), `sskm setop` (per-bucket merge) and `sskm construct` (per-bucket compaction) parallelize with the C++ standard library only (`std::thread`; the CLI links `pthread`), so there is no TBB / parallel-STL runtime dependency.

## Compilation

```bash
mkdir -p build && cd build
cmake .. && make
```

> **Compiler.** `sklib` is C++23 (`CMAKE_CXX_STANDARD 23`) and its widest k-mer backend
> `kuint256` is `unsigned _BitInt(256)`, so `dev`/`main` require **Clang ≥ 16 or GCC ≥ 15**
> (g++ ≤ 14 supports `_BitInt` only in C, not C++). CMake probes the compiler and **fails
> configuration with clear guidance** if it is too old. On Ubuntu 24.04 (GCC tops out at 14)
> build with Clang:
> ```bash
> apt install clang-18
> mkdir -p build && cd build
> cmake -DCMAKE_CXX_COMPILER=clang++-18 .. && make
> ```
> **GCC-13 users.** A dedicated **`gcc-13` branch** — tagged **`v0.7.0-gcc13`** — holds a
> GCC-13-compatible build, kept in sync with `dev` feature-for-feature, so there are two
> working versions in parallel. Check out `gcc-13` if `dev`/`main` does not build with your
> compiler.

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

* `k` — k-mer length in nucleotides. The record integer width is selected automatically (`uint32_t`/`uint64_t`/`__uint128_t`/`kuint256`) from the packed skmer size: the only limit is that a skmer's `2·(2k−m)` bits fit the widest 512-bit `kuint256` pair, i.e. `2·(2k−m) ≤ 512` (so `k` up to ~127 at small `m`, more as `m` grows); `construct` reports a clear error above that.
* `m` — minimizer length in nucleotides, with `1 ≤ m ≤ k`. Smaller `m` produces longer skmers (more shared central nucleotides per group); typical values sit around `m ≈ k/2`.

### `sskm construct`

Builds a sorted skmer list from a FASTA file.

```bash
./bin/sskm construct -k 21 -m 11 -f input.fa -o out.sskm
```

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-k, --kmer-size <int>` | yes | — | k-mer length; the record width is auto-selected, the only cap is `2·(2k−m) ≤ 512` (see above). |
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

The order on minimizers uses a fixed, invertible hash (φ) rather than the raw lexicographic value. Hash-ordered minimizers have lower density than lexicographic ones, so the list holds **9–21% fewer super-k-mers** (a smaller index / fewer bits per k-mer) across genomes, and on repeat-rich genomes it also sharply lowers peak construction RAM (e.g. *C. elegans* `k=21, m=11`: ~341 MB → ~63 MB) by breaking the over-selection of low-complexity minimizers. Each k-mer is stored in a strand-invariant canonical frame derived from its own nucleotides; super-k-mers whose minimizer is *ambiguous* (an RC-palindrome, or repeated within the window) are re-framed per k-mer — picking the minimizer occurrence by minimal φ-rank, then most central, then smaller full interleaved value — and re-merged into compact super-k-mers. This makes queries exact at any `m` (so there is no minimum-`m` restriction) and also removes a rare construction drop of poly-A k-mers. The on-disk minimizer slot is hash-permuted and the binary file carries a format version. The current format is **`VSKMER_4`**: each record is **minimizer-prefix quotiented** (the top `b = log2(buckets)` φ-minimizer bits are constant within a bucket, so the bucket id implies them and they are dropped) and stored in a **runtime-selected integer width** — the smallest of `uint32_t`/`uint64_t`/`__uint128_t`/`kuint256` (pair capacities 64/128/256/512 bits) whose pair holds the remaining `2·(2k−m) − b` bits. Legacy bucketed **`VSKMER_3`** and single-list **`VSKMER_2`** files still load (as non-quotiented 64-bit lists); pre-φ raw-order **`VSKMER_M`** files are rejected on load rather than queried incorrectly. Rebuild rejected lists — and, since the v0.4.2 framing fix re-frames a few rare *ambiguous* k-mers, rebuild any index built before v0.4.2.

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

Set operations between two sorted skmer lists **A** and **B**. The atomic element is a *k-mer*: the operation decomposes per minimizer bucket and per minimizer-position column, and a single two-cursor merge over the two sorted lists yields all three relations (`A∩B`, `A\B`, `B\A`) in one pass. There are two modes:

* **Single-op** (`--op <name>`): `intersection`, `union`, `diff` (asymmetric, `A \ B`), and `xor` (symmetric difference, `A △ B`), plus `intersection_size`, `union_size`, `diff_size`, and `xor_size` variants that report only the result **cardinality**.
* **Combined / single-pass** (any subset of `--inter-out`/`--union-out`/`--diff-ab-out`/`--diff-ba-out`/`--xor-out` and/or `--sizes`): materialize several relations — including **`B \ A`** — and/or report every cardinality in **one** merge pass. See below.

```bash
# materialize A ∩ B into a new list
./bin/sskm setop --op intersection -a A.sskm -b B.sskm -o inter.sskm
# materialize A △ B (symmetric difference: k-mers in exactly one list) into a new list
./bin/sskm setop --op xor -a A.sskm -b B.sskm -o xor.sskm
# just the cardinality (nothing written) — e.g. the numerator/denominator of a Jaccard index
./bin/sskm setop --op union_size -a A.sskm -b B.sskm
# fast materialization when a larger output is acceptable (skip super-k-mer re-compaction)
./bin/sskm setop --op union --no-compact -a A.sskm -b B.sskm -o union.sskm
# parallel by bucket (default 8 threads); the output is byte-identical for any -t
./bin/sskm setop --op intersection -a A.sskm -b B.sskm -o inter.sskm -t 16
# combined mode: several relations (incl. B \ A) in a single merge pass
./bin/sskm setop -a A.sskm -b B.sskm --inter-out i.sskm --union-out u.sskm --diff-ab-out dab.sskm
# combined mode: every cardinality (|A∩B|, |A∪B|, |A\B|, |B\A|, |A△B|, |A|, |B|) in one pass, nothing written
./bin/sskm setop -a A.sskm -b B.sskm --sizes
```

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--op <name>` | single-op mode | — | One of `intersection`, `union`, `diff`, `xor`, `intersection_size`, `union_size`, `diff_size`, `xor_size`. `diff` is asymmetric: `A \ B` (k-mers of A absent from B); `xor` is the symmetric difference `A △ B` (k-mers in exactly one list). Mutually exclusive with the combined-mode options. |
| `-a, --list-a <path>` | yes | — | First sorted skmer list. |
| `-b, --list-b <path>` | yes | — | Second sorted skmer list. |
| `-o, --output <path>` | for `intersection`/`union`/`diff`/`xor` | — | Single-op output list; ignored by the `*_size` variants, which print a count to stdout. |
| `--inter-out <path>` | combined mode | — | Write `A ∩ B` to this list. |
| `--union-out <path>` | combined mode | — | Write `A ∪ B` to this list. |
| `--diff-ab-out <path>` | combined mode | — | Write `A \ B` to this list. |
| `--diff-ba-out <path>` | combined mode | — | Write `B \ A` to this list. |
| `--xor-out <path>` | combined mode | — | Write `A △ B` (symmetric difference) to this list. |
| `--sizes` | combined mode | off | Print all five cardinalities (`|A∩B|`, `|A∪B|`, `|A\B|`, `|B\A|`, `|A△B|`) plus `|A|` and `|B|` to stdout, computed in the same pass. On its own it materializes nothing. |
| `--no-compact` | no | off | Emit **one record per result k-mer** instead of re-compacting the result back into super-k-mers — much faster (it skips the dominant cost) at the price of a larger output file (~3–5×). The result is still a valid, queryable sorted list. Ignored by the `*_size` variants; applies to every combined-mode output. |
| `-t, --threads <N>` | no | 8 | Worker threads for the per-bucket merge — the independent minimizer buckets are processed in parallel. Every output is **byte-identical** regardless of `N`. |

**Both lists must be built with the same `k`, `m`, and `--buckets`**: the merge aligns bucket *i* of A with bucket *i* of B, so the bucket layouts must match (this is checked, and a mismatch is rejected — rebuild both with identical parameters). The independent buckets are merged **in parallel** (`-t`); each worker holds only its current bucket pair, so peak RAM stays low (a handful of buckets) regardless of list size — far below KMC/CBL — and grows modestly with the thread count.

Two cost regimes:

* **Cardinality (`*_size`)** counts the three relations in one merge pass **without writing anything** — the fast path for a Jaccard index, containment, or de-duplication, and the right choice when the result *set* is not needed.
* **Materialized (`intersection`/`union`/`diff`/`xor`)** additionally re-compacts the kept k-mers back into super-k-mers and writes a new list. That re-compaction dominates the time; `--no-compact` skips it (one record per k-mer, ~3–5× larger output, ~3× faster) while keeping the output a valid queryable list.

#### Combined mode (single pass)

Instead of `--op`, pass any subset of `--inter-out`, `--union-out`, `--diff-ab-out`, `--diff-ba-out`, `--xor-out` (one result list each — note `B \ A` is only reachable here) and/or `--sizes`, and the single merge pass produces them all at once: the bucket read and the two-cursor merge are shared, so only the per-output re-compaction is paid per output, and all seven counts (`|A∩B|`, `|A∪B|`, `|A\B|`, `|B\A|`, `|A△B|`, `|A|`, `|B|`) come for free with `--sizes`. Combined mode is **mutually exclusive** with `--op`, and two outputs pointing at the same path are rejected. Each materialized file is **byte-identical** to the matching single-op materialization (for any `-t`). Combined counting costs about one merge pass (≈4× a single `*_size` run, since the merge is the whole cost), while combined materialization is only ≈1.1–1.3× a single materialization (re-compaction is per-output and cannot be shared). Cross-validated against KMC by `tests/setop_multi_verif.sh`.

Set operations scale ~×6–9 with `-t` (near-linear to 8 cores), which makes sklib the fastest set-op tool from a few cores up — including the high-overlap materialized intersection where KMC led single-core (chr1 ∩: 2.6 s vs KMC 2.7 s at 22 cores). Three benchmark reports accompany this feature: [`benchmark/results/journals/SETOPS_REPORT.md`](benchmark/results/journals/SETOPS_REPORT.md) compares sklib against KMC, CBL and FMSI (single- and multi-core), [`benchmark/results/journals/SETOPS_BOTTLENECKS.md`](benchmark/results/journals/SETOPS_BOTTLENECKS.md) breaks down where the time goes and documents the implemented speedups, and [`benchmark/results/journals/SETOPS_MULTI_REPORT.md`](benchmark/results/journals/SETOPS_MULTI_REPORT.md) covers combined single-pass mode vs the sequential single-op runs. See also the wiki's *Set operations* and *Benchmark · Set operations* pages.

## Benchmarks

A complete cross-tool benchmark suite lives in [`benchmark/`](benchmark/) — construction,
membership query (single & streamed), and set operations, comparing sklib against **KMC**,
sshash, **SBWT** (C++ & Rust), **CBL**, **BQF** and **FMSI**.

- **How to run it / what's measured:** [`benchmark/README.md`](benchmark/README.md) (experiments, metrics, capability matrix, datasets).
- **Latest full results** (the authoritative comparison): [`benchmark/results/runs/full_run_2026-06/RESULTS.md`](benchmark/results/runs/full_run_2026-06/RESULTS.md) — raw CSVs in [`…/data/`](benchmark/results/runs/full_run_2026-06/data/).
- **All runs** (index of every campaign): [`benchmark/results/runs/README.md`](benchmark/results/runs/README.md).

Benchmarking needs `kmc`/`kmc_tools` on `PATH` (membership-query oracle + set-op baseline);
optional competitors build via `bash benchmark/scripts/tools/setup.sh`.
