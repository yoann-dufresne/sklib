# sklib — compact sorted k-mer library

A C++23 library for sorting and querying k-mers in a compact representation. sklib packs k-mers and super-k-mers (skmers) as 2-bit-per-nucleotide integers, builds a disk-backed, **bucketed** sorted skmer list (`VirtualSkmer`), and answers membership queries — each routed to a single minimizer bucket — through the `sskm` CLI.

## Documentation

A page-by-page walkthrough of the internals lives in the [project wiki](https://github.com/yoann-dufresne/sklib/wiki): the 2-bit nucleotide encoding, the minimizer-centered interleaved super-k-mer layout, super-k-mer generation (k-mers are canonicalized at yield so a k-mer is always stored under a single orientation), sorted-list construction (with the Fenwick-tree colinear chaining), the bucketed on-disk format (the list is split into minimizer-prefix sub-lists), and the k-mer search algorithm (each query is routed to one bucket).

## Dependencies

* zlib : https://zlib.net/ (apt install zlib1g-dev)
* bzip2 : https://sourceware.org/bzip2/ (apt install libbz2-dev)

Third-party dependencies fetched automatically: `CLI11`, `zlib`, `gtest`, `kseqpp`. sklib is single-threaded (no parallel/TBB dependency).

## Compilation

```bash
mkdir -p build && cd build
cmake .. && make
```

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

The `sskm` tool (built to `build/bin/sskm`) exposes the library's two primary operations: `construct` builds a sorted super-k-mer list from a FASTA input, and `query` looks k-mers up against an existing list. Run `sskm --help`, `sskm construct --help`, or `sskm query --help` for the in-CLI reference; `sskm --version` prints the current version.

### Parameters `k` and `m`

* `k` — k-mer length in nucleotides. The record integer width is selected automatically (`uint32_t`/`uint64_t`/`__uint128_t`) from the packed skmer size, so roughly `1 ≤ k ≤ 63` (a skmer needs `2·(2k−m)` bits, capped at the 256-bit `__uint128_t` pair).
* `m` — minimizer length in nucleotides, with `1 ≤ m ≤ k`. Smaller `m` produces longer skmers (more shared central nucleotides per group); typical values sit around `m ≈ k/2`.

### `sskm construct`

Builds a sorted skmer list from a FASTA file.

```bash
./bin/sskm construct -k 21 -m 11 -f input.fa -o out.sskm
```

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-k, --kmer-size <int>` | yes | — | k-mer length (1 ≤ k ≤ 32). |
| `-m, --minimizer-size <int>` | yes | — | Minimizer length (1 ≤ m ≤ k). |
| `-f, --file <path>` | no | stdin | Input FASTA (plain or gzip). |
| `-o, --output <path>` | no | stdout | Output sorted skmer list. |
| `--ascii` | no | off | Write human-readable ASCII instead of the default binary format. |
| `--buckets <int>` | no | 4096 | Number of minimizer-prefix buckets (rounded down to a power of two). The list is split into that many independently-sorted sub-lists, partitioned by the high-order bits of the minimizer; a query is routed to a single sub-list. More buckets ⇒ lower peak construction RAM and a faster query; `1` is a single list. |
| `--max-ram <size>` | no | — | Approximate target peak RAM (e.g. `2G`, `512M`); a histogram pass derives adaptive buckets balanced to this budget, overriding `--buckets`. Requires `-f <file>` (the input is read twice). |
| `--tmp-dir <path>` | no | next to output | Directory for the temporary bucket files. |

The default binary output embeds an endianness marker so files can be moved across machines. The ASCII format is intended for inspection and tests, not for downstream querying.

#### Bucketed construction and on-disk layout

By default, binary construction to a regular file (`-o`) is **disk-backed and low-memory**: super-k-mers are partitioned by the high-order bits of their minimizer into on-disk buckets, then each bucket is sorted, deduplicated, and built independently. The buckets are kept as **separate sorted sub-lists** in the output (the binary `VSKMER_3` format stores a per-bucket directory — see the wiki's *On-disk Format and Serialization* page), so a query can be routed to one sub-list instead of scanning the whole list. Peak construction RAM is bounded by the **largest bucket** rather than by the whole input, tunable with `--buckets` (or targeted with `--max-ram`); construction is also *faster* than an all-in-RAM build, since the column/chaining algorithm runs on small per-bucket inputs. Measured single-threaded in a Release build, the bucketed path keeps construction peak RSS to tens–hundreds of MB even on large genomes — far below an all-in-RAM build (e.g. human chr1 `k=21, m=11`: ~198 MB peak RSS, ~55 s single-threaded Release) — see [`scripts/bench/README.md`](scripts/bench/README.md). Temporary bucket files are removed automatically on success and on error.

`--ascii` and writing to stdout (no `-o`) fall back to an all-in-RAM build that emits a single bucket (a stream cannot seek back to patch the header's per-bucket directory and total count).

#### Minimizer ordering (compact index)

The order on minimizers uses a fixed, invertible hash (φ) rather than the raw lexicographic value. Hash-ordered minimizers have lower density than lexicographic ones, so the list holds **9–21% fewer super-k-mers** (a smaller index / fewer bits per k-mer) across genomes, and on repeat-rich genomes it also sharply lowers peak construction RAM (e.g. *C. elegans* `k=21, m=11`: ~341 MB → ~63 MB) by breaking the over-selection of low-complexity minimizers. Each k-mer is stored in a strand-invariant canonical frame derived from its own nucleotides; super-k-mers whose minimizer is *ambiguous* (an RC-palindrome, or repeated within the window) are re-framed per k-mer — picking the minimizer occurrence by minimal φ-rank, then most central, then smaller full interleaved value — and re-merged into compact super-k-mers. This makes queries exact at any `m` (so there is no minimum-`m` restriction) and also removes a rare construction drop of poly-A k-mers. The on-disk minimizer slot is hash-permuted and the binary file carries a format version. The current format is **`VSKMER_3`** (it adds the per-bucket directory); older single-list φ-order files (**`VSKMER_2`**) still load — as a single bucket — while pre-φ raw-order files (**`VSKMER_M`**) are rejected on load rather than queried incorrectly. Rebuild rejected lists.

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

`-i` and the positional `sequence` argument are mutually exclusive — exactly one must be provided.

Internally, each query super-k-mer is **routed to the single bucket** its minimizer falls in (the high-order bits pick the sub-list); only that bucket is loaded from disk (lazily, then cached) and searched, so a query reads a fraction of the index rather than the whole file. See the wiki's *Querying a Sorted List* and *The Kmer Search Algorithm* pages.
