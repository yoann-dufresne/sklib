# sklib — compact sorted k-mer library

A C++23 library for sorting and querying k-mers in a compact representation. sklib packs k-mers and super-k-mers (skmers) as 2-bit-per-nucleotide integers, builds a disk-backed sorted skmer list (`VirtualSkmer`), and answers membership queries through the `sskm` CLI.

## Documentation

A page-by-page walkthrough of the internals lives in the [project wiki](https://github.com/yoann-dufresne/sklib/wiki): the 2-bit nucleotide encoding, the minimizer-centered interleaved super-k-mer layout, super-k-mer generation (k-mers are canonicalized at yield so a k-mer is always stored under a single orientation), sorted-list construction (with the Fenwick-tree colinear chaining), the on-disk format, and the k-mer search algorithm.

## Dependencies

* zlib : https://zlib.net/ (apt install zlib1g-dev)
* bzip2 : https://sourceware.org/bzip2/ (apt install libbz2-dev)

Third-party dependencies fetched automatically: `CLI11`, `zlib`, `gtest`, `kseqpp`. TBB (Intel Threading Building Blocks) is auto-detected if installed.

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

* `k` — k-mer length in nucleotides. With the default 64-bit backend the packed 2-bit encoding allows `1 ≤ k ≤ 32`.
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
| `--buckets <int>` | no | 4096 | Number of on-disk minimizer buckets for the low-memory path (rounded down to a power of two). More buckets ⇒ lower peak RAM; `1` is a single bucket. |
| `--max-ram <size>` | no | — | Approximate target peak RAM (e.g. `2G`, `512M`); a histogram pass derives adaptive buckets balanced to this budget, overriding `--buckets`. Requires `-f <file>` (the input is read twice). |
| `--tmp-dir <path>` | no | next to output | Directory for the temporary bucket files. |

The default binary output embeds an endianness marker so files can be moved across machines. The ASCII format is intended for inspection and tests, not for downstream querying.

#### Low-memory construction

By default, binary construction to a regular file (`-o`) is **disk-backed and low-memory**: super-k-mers are partitioned by minimizer into on-disk buckets, then each bucket is sorted, deduplicated, and built independently and appended to the output. Because the minimizer occupies the most significant bits of a super-k-mer, concatenating buckets in order reproduces the exact same globally-sorted list — the output format is unchanged and queries are unaffected. Peak RAM is therefore bounded by the **largest bucket** rather than by the whole input, and is tunable with `--buckets` (or targeted with `--max-ram`). On human chr21 (`k=31, m=13`) this cuts construction peak RSS from ~850 MB (all-in-RAM) to ~67 MB at the default `--buckets 4096`. Temporary bucket files are removed automatically on success and on error.

`--ascii` and writing to stdout (no `-o`) fall back to the historical all-in-RAM path (the disk-backed path patches the header record count with a seek, which a stream cannot do).

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
