# Full k-mer tooling benchmark — experiment & data documentation

Single-machine comparison of **sklib** against state-of-the-art k-mer set tools on
three workloads — **index construction**, **membership query**, and **set
operations** — over a ladder of genomes and k-mer sizes.

- **Machine:** `yoann-Precision-5490` — Intel® Core™ Ultra 7 165H, 22 logical
  cores, 62 GiB RAM, Linux 6.17.
- **Campaign window:** 2026-06-08 → 2026-06-16 (append-only CSVs; the canonical
  sklib build is **0.11.0**, see *Versioning* below).
- **Timing:** every measured command is wrapped in GNU `/usr/bin/time -v`; wall
  time is the median of `REPS=3` repetitions, peak RSS is the max over reps.
- **Raw data:** `data/{construct,query_single,query_stream,setop}.csv` (this folder).
- **Detailed results & analysis:** `RESULTS.md` (this folder).

## Tools compared

| Tool | Version (CSV `tool_version`) | Role | Query? | Set-ops? | k limit | Threads |
|------|------------------------------|------|--------|----------|---------|---------|
| **sklib** | `sklib-0.11.0` | structure under test (sorted skmer list) | ✅ | ✅ (+ size-only, joint) | k≤127 | construct & query & setop sweep |
| **kmc** | `kmc` | oracle + construct + set-ops | ❌ (oracle) | ✅ (materialize) | — | sweep |
| **cbl** | `cbl` | conway-bromage-lyndon | ✅ | ✅ (materialize) | **odd k ≤ 59** | 1 |
| **sshash** | `sshash-…` | minimizer hash (membership) | ✅ | ❌ | k≤31 | 1 |
| **sbwt** | `sbwt-plain-matrix` | C++ SBWT (algbio) | ✅ | ❌ | **k≤32** (MAX_KMER_LENGTH) | 1 |
| **sbwtrs** | `sbwt-rs-0.2.0` | **Rust SBWT** (jnalanko/sbwt-rs-cli) | ✅ | ✅ (materialize) | k≤255 (full grid) | sweep |
| **fmsi** | `fmsi-39c7` | FM-index masked superstring (KmerCamel) | ✅ | ✅ (materialize) | k≤127 | 1 |
| **bqf** | `bqf-s18(approx)` | quotient filter (**approximate**) | ✅ | ❌ | per-(k,z) | 1 |

`sbwt` (C++) and `sbwtrs` (Rust) are the **same algorithm, two implementations**,
kept as distinct tools so their rows never collide and the rewrite can be compared
head-to-head. Both build canonical indexes (`--add-reverse-complements` / `-r`),
matching sklib/kmc semantics.

## Datasets (sanitized FASTA)

| name | size | note |
|------|------|------|
| sarscov2 | 30 KB | tiny (startup-dominated) |
| ecoli | 4.6 MB | bacterial |
| yeast | 12 MB | eukaryote (small) |
| celegans | 98 MB | eukaryote (medium) |
| chr21 | 39 MB | human chromosome |
| chr1 | 224 MB | human chromosome (largest run) |
| chm13 | 3.0 GB | **deferred** — not in this run |

## Parameter grid

- **(k, m):** `15,7 · 21,11 · 31,15 · 41,19 · 51,25 · 63,31` (m = sklib minimizer length).
- **threads:** `1 2 4 8 16` (tools that honour `-t`; others run single-thread).
- **query presence** (`query_single`/`query_stream`): `0 25 50 75 100` % of queried
  k-mers present in the index.
- **set-op Jaccard target:** `0 0.1 0.3 0.5 0.7 0.9 1.0` (B is a mutated copy of A
  tuned to that Jaccard; the achieved Jaccard is recorded too).
- **set-op ops:** `inter, union, diffab (A\B), diffba (B\A)` and `joint` (all at once);
  **modes:** `materialize` (write the result) and `size` (cardinality only).
- **REPS = 3** (median).

## Experiments

1. **construct** — build the index from a genome; records size + bits/k-mer + time.
2. **query_single** — membership of isolated (shuffled) k-mers at each presence level.
3. **query_stream** — membership of consecutive (streamed) k-mers.
4. **setop** — two-list set operations; B is a Jaccard-tuned mutant of A. Authoritative
   cardinalities (the `result_kmers` column) come from one sklib `--sizes` pass and are
   shared across tools, so a tool's set-op result is *not* re-validated here — only its
   time/RAM is measured.

### Coverage caps & legitimate skips
- **fmsi set-ops:** capped at **yeast** (glacial at scale — chr1 ≈ days). construct+query full.
- **sbwtrs set-ops:** capped — full `sarscov2`, `ecoli` k15–41, `yeast` on a **light grid**
  (`JACCARD="0 0.5 1.0" THREADS="1 8"`); celegans/chr21/chr1 skipped (sbwtrs set-ops are
  ~as slow as fmsi). construct+query full on all 6 datasets.
- **sbwt (C++):** k ≤ 32 ⇒ only k = 15, 21, 31. No set-ops.
- **cbl:** odd k only, ≤ 59 ⇒ no k = 63. (All grid k are odd, so 15–51 covered.)
- **sshash:** k ≤ 31; also fails the *stream* query on `celegans k15` (skipped).
- **bqf:** approximate membership only; no set-ops; construct/query on a subset.
- **kmc:** the membership-query oracle, so it has no `query_*` rows (only construct/setop).

## Reproduce

```bash
# 1. build sklib (needs Clang ≥ 16 for the kuint256 backend)
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-18 -DWITH_TESTS=OFF
cmake --build build-bench -j --target sskm

# 2. build competitors (clones+builds into benchmark/data/tools_src/, writes tools.env)
bash benchmark/scripts/tools/setup.sh            # all; or e.g. `... setup.sh sbwtrs`

# 3. fetch genomes
bash benchmark/scripts/fetch_genomes.sh sarscov2 ecoli yeast celegans chr21 chr1

# 4. run an experiment (env-parametrized; defaults = full grid)
DATASETS="ecoli yeast" KM="21,11 31,15" TOOLS="sklib kmc" bash benchmark/scripts/construct.sh
#   drivers: construct.sh · query_single.sh · query_stream.sh · setop.sh
```

Results append to `benchmark/results/latest/<experiment>.csv` (git-ignored); a run is
**resumable** — a measurement already present (keyed by the identity columns below) is skipped.

---

# CSV schema

Four append-only CSVs, one per experiment. All share the identity prefix
`timestamp, host, cpu, tool, tool_version, dataset, k, m, threads`. A row is a single
measurement; **resume** dedups on (tool, tool_version, host, dataset, k, m, threads,
+ experiment-specific discriminators).

> ⚠️ **peak RSS units differ by file:** `construct`/`query_*` store **`peak_rss_kb`
> (kibibytes)**; `setop` stores **`peak_rss_mb` (mebibytes, integer)** — the setop
> driver converts before writing. Mind the column name.

### `construct.csv` (16 columns)
| # | column | meaning |
|---|--------|---------|
| 1 | timestamp | ISO-8601 local time the row was written |
| 2 | host | machine hostname |
| 3 | cpu | CPU model string |
| 4 | tool | `sklib`, `kmc`, `cbl`, `sshash`, `sbwt`, `sbwtrs`, `fmsi`, `bqf` |
| 5 | tool_version | version/variant tag (resume key) |
| 6 | dataset | genome name |
| 7 | k | k-mer length |
| 8 | m | sklib minimizer length (tools that ignore it still log the grid value) |
| 9 | threads | thread count for this measurement |
| 10 | time_s | median wall time (s) |
| 11 | **peak_rss_kb** | peak resident memory (**KB**) |
| 12 | index_bytes | on-disk index size (bytes); `NA` if not applicable |
| 13 | bits_per_kmer | `index_bytes*8 / n_kmers` |
| 14 | n_kmers | distinct k-mers in the index |
| 15 | n_superkmers | sklib super-k-mer count (`NA` for others) |
| 16 | throughput_Mkmer_s | `n_kmers / time_s / 1e6` |

### `query_single.csv` & `query_stream.csv` (16 columns, identical schema)
| # | column | meaning |
|---|--------|---------|
| 1–9 | timestamp … threads | identity prefix (as above) |
| 10 | presence | % of queried k-mers present in the index (0/25/50/75/100) |
| 11 | n_kmers | total k-mers queried |
| 12 | present_kmers | how many were present |
| 13 | absent_kmers | how many were absent |
| 14 | time_s | median wall time (s) |
| 15 | **peak_rss_kb** | peak resident memory (**KB**) |
| 16 | throughput_Mkmer_s | `n_kmers / time_s / 1e6` (higher = better) |

`query_single` = isolated/shuffled k-mers; `query_stream` = consecutive k-mers
(exercises streaming locality). `presence` discriminates the resume key.

### `setop.csv` (18 columns)
| # | column | meaning |
|---|--------|---------|
| 1–9 | timestamp … threads | identity prefix (as above) |
| 10 | op | `inter`, `union`, `diffab` (A\B), `diffba` (B\A), `joint` (all-at-once) |
| 11 | mode | `materialize` (write result list) or `size` (cardinality only; sklib only) |
| 12 | joint | `1` if this row is the combined single-pass op, else `0` |
| 13 | jaccard_target | requested Jaccard between A and B |
| 14 | jaccard_measured | achieved Jaccard (B is a tuned random mutant of A) |
| 15 | result_kmers | authoritative result cardinality (from sklib `--sizes`, shared across tools) |
| 16 | time_s | median wall time (s) |
| 17 | **peak_rss_mb** | peak resident memory (**MB** — note: not KB) |
| 18 | throughput_Mkmer_s | result throughput |

Resume key adds `op, mode, jaccard_target` to the identity prefix.

## Versioning note
The CSVs are append-only across the campaign, so `sklib` appears at **three**
`tool_version`s — `sklib-0.8.0`, `sklib-0.10.1`, `sklib-0.11.0` — from successive
re-measurements as the library was optimized. **`sklib-0.11.0` is the canonical
current build**; the older rows are kept for before/after history and should be
filtered out (`tool_version == "sklib-0.11.0"`) for the head-to-head comparison.
Competitor rows are single-version. A short-lived `sklib-0.12.0` (experimental
uint8/uint16 backends) was measured then reverted; its rows were purged.
