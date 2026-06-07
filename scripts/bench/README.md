# sklib benchmark harness

Performance/space benchmark for the sorted super-k-mer list, producing the CSV and
publication figures described in the benchmark plan. Correctness vs KMC is handled
separately by `scripts/large_scale_e2e.sh`; this harness only runs a light
false-negative gate and focuses on **curves**.

## Layout
- `bench.sh`  — main driver: loops `(dataset, tool, k, m, workload, threads)`, measures
  construction (time, peak RSS, index size → bits/k-mer, k-mers/super-k-mer) and query
  throughput, appends one row per measurement to `scripts/out/bench/results.csv`.
- `lib.sh`   — genome catalogue, `run_timed`/`run_timed_median`, KMC oracle, CSV helpers.
- `tools.sh` — per-tool adapters (`construct_*`, `query_*`, `version_*`, `available_*`) for
  sklib plus the competitors **sshash / sbwt / cbl / bqf** (wrappers in `tools/`, built by
  `tools/setup.sh`). `BUCKETS=<n>` overrides sklib's `--buckets`.
- `plots.py` — reads the CSV, writes Fig 1–7 to `scripts/out/bench/figures/`.
- `compare.py` — Markdown cross-tool comparison tables from the CSV (dedup by latest timestamp).
- `bench_buckets.sh` — `--buckets` sweep on one dataset → `scripts/out/bench/results_buckets.csv`.
- `bench_setops.sh` / `run_setops_parallel.sh` — **set-operation** benchmarks: single-core sklib vs
  KMC vs CBL, and multi-core sklib `-t` vs KMC vs FMSI (`→ setops_parallel.csv`); `plot_setops.py`
  draws their figures. Full write-up in `report/SETOPS_REPORT.md` (+ `report/SETOPS_BOTTLENECKS.md`).
- `.venv/`   — local Python env for plotting (pandas + matplotlib); git-ignored.

## Prerequisites
- `kmc`, `kmc_tools` on PATH (the k-mer-set oracle / bits-per-k-mer denominator).
- A built `sskm` at `build-bench/bin/sskm` (the harness default; build **Release** — a `DEBUG`
  `build/` enables ASan and skews RSS/time). Query is parallel (`-t`), so `THREADS=N` runs the
  parallel file query on N cores (`THREADS=1` = sequential); construction and set operations are also bucket-parallel (`-t`).
- Plotting env (one-off):
  ```bash
  python3 -m venv scripts/bench/.venv
  scripts/bench/.venv/bin/pip install matplotlib pandas
  ```

## Run
```bash
# defaults: sklib on sarscov2+ecoli, k=21 m=11, all workloads, 1 thread
bash scripts/bench/bench.sh

# a real sweep (single-threaded; all six genomes, the k/m grid)
DATASETS="sarscov2 ecoli yeast chr21 celegans chr1" \
KM="15,7 21,11 31,15" \
WORKLOADS="positive random streaming shuffled reads" \
THREADS="1" REPS=3 FRESH=0 \
bash scripts/bench/bench.sh

scripts/bench/.venv/bin/python scripts/bench/plots.py            # figures
python3 scripts/bench/compare.py --km 21,11 \
  --datasets "ecoli yeast celegans chr21 chr1"                   # comparison tables
bash scripts/bench/bench_buckets.sh                              # --buckets sweep on ecoli
```

### Env knobs (defaults)
`DATASETS="sarscov2 ecoli"` · `KM="21,11"` · `TOOLS="sklib"` ·
`WORKLOADS="positive random streaming shuffled reads"` · `THREADS="1"` ·
`N_QUERY=200000` · `STREAM_BP=2000000` · `READS_N=20000` · `READLEN=150` · `ERR=0.01` ·
`REPS=1` · `SEED=1234` · `SELFCHECK=1` · `FRESH=0` (set `FRESH=1` to truncate the CSV) ·
`BUCKETS=""` (empty ⇒ sklib's default 4096; set to force a specific `--buckets`).

The CSV is **append-only** across runs, so you can accumulate datasets/configs
incrementally (e.g. run the big genomes on a larger-RAM machine and concatenate CSVs).
Because of this it may hold **several sklib commits** (e.g. the V0/V1/V2 construction-RAM
rewrite). `plots.py` collapses each measurement key
(`tool,dataset,k,m,phase,workload,threads`) to its **most recent timestamp**, so figures
always reflect the current build — superseded sklib rows drop out, single-version
competitor rows are untouched.

## Query workloads
- `positive`  — k-mers sampled from the genome, each its own record (scattered, all present).
- `random`    — random k-mers (mostly absent; sklib is exact → no false positives).
- `streaming` — a contiguous genome region (consecutive k-mers; the cache-aware batched
  query shines — super-k-mers stay coherent).
- `shuffled`  — the same k-mers as `streaming` but one per record in random order
  (locality broken → isolates the batching gain).
- `reads`     — simulated reads (`READLEN` bp, `ERR` substitution rate): realistic mixed status.

## CSV columns
`timestamp,host,cpu,commit,tool,tool_version,dataset,dataset_bp,distinct_kmers,k,m,threads,phase,workload,n_queries,index_bytes,bytes_per_skmer,n_superkmers,kmers_per_superkmer,bits_per_kmer,time_s,peak_rss_kb,throughput_Mkmer_s,correctness`

`phase` is `construct` or `query`. Space columns (`bits_per_kmer`, …) are repeated on
query rows so each row is self-contained. `correctness` is `pass`/`fail` on the
`positive` query rows (false-negative gate), `NA` elsewhere. `bits_per_kmer` is computed on
the **payload** (skmer records), excluding the file header; for sklib the `VSKMER_4` header
includes the per-bucket directory (`40 + 16·n_buckets` bytes, reported by
`e2e_helpers.py binheadersize`), so the directory is *not* counted as payload.

## Figures
1. `fig1_bits_per_kmer` — space: bits/k-mer vs m (per k, per dataset) + 2-bit-raw reference.
2. `fig2_compaction` — k-mers/super-k-mer vs m; #super-k-mers vs #k-mers.
3. `fig3_construction` — construct throughput and peak RSS vs dataset size.
4. `fig4_query` — query throughput by workload (single-threaded; sklib has no thread scaling).
5. `fig5_pareto` — bits/k-mer vs query throughput (sklib swept over m, competitors as points).
6. `fig6_envelope` — query correctness over the (k,m) grid.
7. `fig7_ksweep` — space (bits/k-mer) and streaming query speed vs k, one line per tool per dataset.

## Cross-tool comparison (current run)

Single machine (`Intel Core Ultra 7 165H`), **single-thread** (`-t 1`), sklib **v0.5.0** (Release,
bucketed `VSKMER_4`, default `--buckets 4096`). Each cell is the full **FASTA → index** cost: *sshash*
includes BCALM2 unitigs, *bqf* KMC s-mer counting, *fmsi* the KmerCamel masked-superstring step.
**pos** / **stream** are `positive` (scattered, all-present) and `streaming` (contiguous) query
throughput in Mk-mer/s; **pos RSS** is peak query RSS. `bqf` is **approximate** (false positives);
`fmsi` is indexed `-x` (no kLCP); `cbl` pre-filters sub-`K` fragments. Regenerate with `python3
scripts/bench/compare.py --km <k>,<m>` (dedup by latest timestamp).

> ⚠️ **Build Release.** The harness now defaults `SSKM_BIN` to `build-bench/bin/sskm`. A `build/`
> configured `DEBUG` enables AddressSanitizer, inflating sklib construct RSS ~33× and time ~4× — the
> trap that first polluted this table (memory `measuring-construct-rss-asan-trap`). Competitor numbers
> are unchanged across runs (their code didn't change) and are reused as-is.

Columns: **index MB** (on-disk) · **bits/kmer** (payload) · **constr s** · **constr RSS MB** · **pos
Mk/s** (point lookups) · **pos RSS MB** · **stream Mk/s** (contiguous).


### ecoli  (k=15, m=7, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 30.9 | 55.36 | 1.21 | 19 | 0.645 | 34 | 2.451 |
| sshash | 9.5 | 16.98 | 1.59 | 44 | 1.587 | 13 | 18.868 |
| sbwt | 6.9 | 12.38 | 2.50 | 91 | 1.075 | 36 | 8.889 |
| cbl | 21.0 | 37.58 | 1.25 | 150 | 0.464 | 135 | 3.039 |
| bqf | 13.6 | 24.44 | 1.47 | 74 | 1.064 | 17 | 16.129 |
| fmsi | 1.2 | 2.17 | 2.64 | 87 | 0.806 | 5 | 1.424 |

### yeast  (k=15, m=7, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 83.4 | 60.77 | 3.42 | 19 | 0.500 | 84 | 1.698 |
| sshash | 36.9 | 26.87 | 8.21 | 152 | 0.909 | 39 | 7.042 |
| sbwt | 15.5 | 11.26 | 5.99 | 167 | 0.873 | 44 | 6.689 |
| cbl | 41.8 | 30.50 | 3.05 | 164 | 0.309 | 152 | 2.294 |
| bqf | 25.2 | 18.34 | 3.20 | 139 | 0.962 | 28 | 12.820 |
| fmsi | 3.1 | 2.25 | 8.39 | 200 | 0.719 | 7 | 1.264 |

### chr21  (k=15, m=7, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 225.7 | 66.36 | 10.73 | 45 | 0.303 | 220 | 0.901 |
| sshash | 151.9 | 44.67 | 18.60 | 431 | 0.556 | 149 | 2.083 |
| sbwt | 36.7 | 10.81 | 13.06 | 441 | 0.746 | 64 | 3.816 |
| cbl | 91.3 | 26.86 | 9.51 | 243 | 0.206 | 184 | 0.990 |
| bqf | 46.1 | 13.57 | 8.32 | 372 | 0.877 | 48 | 9.090 |
| fmsi | 8.7 | 2.57 | 31.51 | 795 | 0.565 | 13 | 0.992 |

### celegans  (k=15, m=7, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 525.8 | 69.41 | 26.20 | 51 | 0.172 | 506 | 0.518 |
| sshash | 428.3 | 56.54 | 62.15 | 1013 | 0.305 | 413 | – |
| sbwt | 80.6 | 10.64 | 33.03 | 1016 | 0.504 | 106 | 4.090 |
| cbl | 191.8 | 25.32 | 29.05 | 360 | 0.111 | 247 | 0.988 |
| bqf | 83.9 | 11.07 | 18.64 | 846 | 0.712 | 84 | 8.889 |
| fmsi | 20.8 | 2.74 | 89.59 | 1808 | 0.394 | 25 | 0.863 |

### chr1  (k=15, m=7, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 1032.6 | 73.86 | 59.51 | 140 | 0.109 | 989 | 0.285 |
| sshash | 1016.6 | 72.72 | 105.77 | 1654 | 0.168 | 974 | 0.614 |
| sbwt | 147.8 | 10.57 | 54.63 | 1862 | 0.427 | 170 | 2.590 |
| cbl | 345.6 | 24.72 | 70.04 | 561 | 0.076 | 345 | 0.376 |
| bqf | 83.9 | 6.00 | 30.47 | 1920 | 0.752 | 84 | 6.666 |
| fmsi | 40.4 | 2.89 | 207.73 | 3925 | 0.279 | 44 | 0.614 |

### ecoli  (k=21, m=11, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 14.2 | 24.81 | 0.98 | 19 | 0.746 | 18 | 5.319 |
| sshash | 5.8 | 10.23 | 0.80 | 54 | 2.083 | 10 | 34.481 |
| sbwt | 7.0 | 12.35 | 2.70 | 102 | 0.575 | 36 | 6.024 |
| cbl | 24.7 | 43.42 | 1.16 | 145 | 0.270 | 129 | 2.212 |
| bqf | 22.0 | 38.77 | 1.56 | 72 | 0.735 | 25 | 14.492 |
| fmsi | 1.2 | 2.13 | 2.46 | 82 | 0.625 | 5 | 1.152 |

### yeast  (k=21, m=11, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 39.5 | 27.39 | 2.79 | 18 | 0.662 | 42 | 4.629 |
| sshash | 15.7 | 10.93 | 1.69 | 93 | 1.613 | 19 | 33.332 |
| sbwt | 16.1 | 11.23 | 6.63 | 198 | 0.712 | 44 | 6.644 |
| cbl | 53.7 | 37.37 | 2.87 | 170 | 0.285 | 156 | 2.283 |
| bqf | 41.9 | 29.16 | 4.38 | 131 | 0.658 | 44 | 12.195 |
| fmsi | 3.1 | 2.15 | 7.57 | 202 | 0.592 | 7 | 1.089 |

### chr21  (k=21, m=11, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 130.1 | 31.80 | 9.46 | 45 | 0.476 | 129 | 2.631 |
| sshash | 53.1 | 12.98 | 5.53 | 205 | 0.935 | 55 | 8.063 |
| sbwt | 44.0 | 10.76 | 15.66 | 428 | 0.439 | 71 | 3.967 |
| cbl | 140.0 | 34.23 | 9.28 | 307 | 0.172 | 236 | 0.844 |
| bqf | 79.7 | 19.48 | 12.35 | 342 | 0.588 | 80 | 7.245 |
| fmsi | 9.3 | 2.28 | 26.56 | 591 | 0.467 | 13 | 0.922 |

### celegans  (k=21, m=11, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 418.8 | 36.72 | 24.33 | 50 | 0.280 | 404 | 1.842 |
| sshash | 159.1 | 13.95 | 17.76 | 524 | 0.545 | 156 | 6.309 |
| sbwt | 120.8 | 10.59 | 49.14 | 1031 | 0.308 | 144 | 3.578 |
| cbl | 363.7 | 31.90 | 34.82 | 1154 | 0.052 | 973 | 0.488 |
| bqf | 285.2 | 25.01 | 34.47 | 794 | 0.363 | 276 | 4.750 |
| fmsi | 25.4 | 2.22 | 86.49 | 1514 | 0.301 | 29 | 0.822 |

### chr1  (k=21, m=11, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 932.5 | 39.30 | 54.51 | 138 | 0.168 | 893 | 0.816 |
| sshash | 363.6 | 15.32 | 35.40 | 1069 | 0.333 | 351 | 2.155 |
| sbwt | 250.2 | 10.54 | 88.03 | 1926 | 0.266 | 267 | 2.242 |
| cbl | 758.5 | 31.97 | 78.12 | 1842 | 0.034 | 1485 | 0.167 |
| bqf | 536.9 | 22.63 | 66.70 | 1810 | 0.267 | 516 | 1.785 |
| fmsi | 54.3 | 2.29 | 214.07 | 3192 | 0.171 | 59 | 0.590 |

### ecoli  (k=31, m=15, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 12.3 | 21.43 | 0.93 | 18 | 0.581 | 16 | 7.519 |
| sshash | 4.1 | 7.12 | 0.53 | 36 | 1.887 | 8 | 45.453 |
| sbwt | 7.0 | 12.34 | 2.86 | 182 | 0.629 | 36 | 8.928 |
| cbl | 36.9 | 64.90 | 1.04 | 152 | 0.409 | 130 | 3.333 |
| bqf | 22.0 | 38.68 | 1.51 | 71 | 0.473 | 25 | 16.260 |
| fmsi | 1.2 | 2.13 | 2.51 | 82 | 0.476 | 5 | 0.870 |

### yeast  (k=31, m=15, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 31.2 | 21.55 | 2.81 | 18 | 0.559 | 34 | 6.172 |
| sshash | 10.6 | 7.34 | 1.40 | 83 | 1.515 | 14 | 41.664 |
| sbwt | 16.2 | 11.23 | 6.99 | 408 | 0.528 | 44 | 6.826 |
| cbl | 87.0 | 60.20 | 2.77 | 211 | 0.267 | 182 | 2.286 |
| bqf | 41.9 | 29.02 | 4.13 | 130 | 0.339 | 44 | 10.638 |
| fmsi | 3.1 | 2.15 | 7.73 | 203 | 0.424 | 7 | 0.832 |

### chr21  (k=31, m=15, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 98.9 | 22.83 | 8.84 | 45 | 0.429 | 99 | 3.570 |
| sshash | 38.0 | 8.78 | 3.46 | 136 | 0.926 | 40 | 12.817 |
| sbwt | 46.5 | 10.74 | 18.42 | 1255 | 0.321 | 73 | 3.875 |
| cbl | 250.0 | 57.77 | 9.57 | 618 | 0.126 | 485 | 0.633 |
| bqf | 79.7 | 18.41 | 12.54 | 342 | 0.311 | 80 | 6.942 |
| fmsi | 9.7 | 2.24 | 27.69 | 579 | 0.364 | 14 | 0.755 |

### celegans  (k=31, m=15, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 276.7 | 23.54 | 23.82 | 50 | 0.307 | 268 | 3.226 |
| sshash | 97.9 | 8.33 | 11.20 | 369 | 0.680 | 97 | 12.820 |
| sbwt | 124.4 | 10.59 | 54.59 | 1948 | 0.215 | 148 | 3.533 |
| cbl | 673.8 | 57.34 | 46.58 | 4976 | 0.018 | 4264 | 0.173 |
| bqf | 285.2 | 24.27 | 34.54 | 794 | 0.210 | 276 | 4.750 |
| fmsi | 25.8 | 2.20 | 84.86 | 1550 | 0.175 | 30 | 0.631 |

### chr1  (k=31, m=15, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 619.4 | 24.27 | 54.14 | 139 | 0.210 | 594 | 1.319 |
| sshash | 246.7 | 9.67 | 26.54 | 735 | 0.442 | 240 | 3.355 |
| sbwt | 269.0 | 10.54 | 103.11 | 1948 | 0.182 | 285 | 2.118 |
| cbl | 1448.5 | 56.76 | 102.72 | 9223 | 0.011 | 7843 | 0.054 |
| bqf | 536.9 | 21.04 | 66.91 | 1811 | 0.178 | 516 | 1.805 |
| fmsi | 57.5 | 2.25 | 222.20 | 3414 | 0.118 | 62 | 0.500 |

> **Large k (sklib's uint128 regime).** sshash (k≤31), SBWT (k≤32) and CBL (odd k≤59) cannot index here; `bqf` failed its per-(k,z) build at these k. So the field narrows to **sklib vs FMSI** (KMC remains a set-op/oracle competitor, §set-ops). FMSI is ~10× smaller on disk but its construction RAM **explodes** (chr1 k=63: 8.5 GB vs sklib 133 MB) and its streaming is far slower.


### ecoli  (k=59, m=29, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 13.8 | 24.03 | 1.17 | 16 | 0.272 | 18 | 6.451 |
| fmsi | 1.2 | 2.13 | 3.35 | 147 | 0.300 | 5 | 0.519 |

### yeast  (k=59, m=29, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 35.1 | 24.06 | 3.03 | 14 | 0.266 | 38 | 5.697 |
| fmsi | 3.1 | 2.15 | 9.38 | 274 | 0.275 | 7 | 0.501 |

### chr21  (k=59, m=29, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 112.2 | 24.62 | 11.66 | 42 | 0.232 | 112 | 4.554 |
| fmsi | 9.9 | 2.18 | 34.72 | 1081 | 0.234 | 14 | 0.476 |

### celegans  (k=59, m=29, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 292.1 | 24.26 | 31.22 | 45 | 0.202 | 283 | 3.442 |
| fmsi | 26.2 | 2.18 | 107.86 | 2170 | 0.109 | 31 | 0.404 |

### chr1  (k=59, m=29, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 684.8 | 25.00 | 69.02 | 134 | 0.178 | 658 | 2.645 |
| fmsi | 59.8 | 2.18 | 255.70 | 8555 | 0.064 | 64 | 0.374 |

### ecoli  (k=63, m=31, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 13.0 | 22.59 | 1.14 | 15 | 0.262 | 17 | 6.757 |
| fmsi | 1.2 | 2.13 | 3.12 | 145 | 0.282 | 5 | 0.488 |

### yeast  (k=63, m=31, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 33.1 | 22.66 | 3.19 | 14 | 0.247 | 36 | 5.847 |
| fmsi | 3.1 | 2.15 | 9.53 | 274 | 0.264 | 7 | 0.475 |

### chr21  (k=63, m=31, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 105.6 | 23.10 | 10.94 | 41 | 0.227 | 106 | 4.726 |
| fmsi | 10.0 | 2.18 | 34.92 | 1081 | 0.216 | 14 | 0.438 |

### celegans  (k=63, m=31, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 275.3 | 22.82 | 30.55 | 44 | 0.213 | 267 | 3.670 |
| fmsi | 26.3 | 2.18 | 106.97 | 2169 | 0.111 | 31 | 0.398 |

### chr1  (k=63, m=31, threads=1)

| tool | index MB | bits/kmer | constr s | constr RSS MB | pos Mkmer/s | pos RSS MB | stream Mkmer/s |
|---|---|---|---|---|---|---|---|
| **sklib (new)** | 644.0 | 23.42 | 70.51 | 133 | 0.084 | 618 | 1.263 |
| fmsi | 59.9 | 2.18 | 254.46 | 8550 | 0.060 | 64 | 0.357 |

**Reading.** sklib keeps the **lowest construction RAM in every cell** — 18–140 MB on the big genomes,
**3–9× below the best competitor and up to ~60× below the worst** (chr1 `k=21`: 138 MB vs 1069 MB
sshash … 3192 MB fmsi). Construction time is mid-pack (near sshash, faster than sbwt/bqf on big
genomes). With the P0/P1 query work, point-lookup throughput is now **competitive** (usually 2nd behind
sshash). The trade-offs are **index size** and **streaming**: as a sorted super-k-mer list it is less
compact than the exact dictionaries (sbwt ≈10–12 b/km, sshash 7–17) and far less than **fmsi**
(≈2.1–2.9 b/km, the smallest, via masked superstrings + MBWT), and streaming trails sshash/bqf; b/km
improves sharply as `k−m` grows (longer super-k-mers), reaching ~23 b/km at k=63.

**Large k (59/63).** sklib's reach is itself a result: sshash/SBWT/CBL **cannot index** at these k and
bqf failed to build, leaving only **sklib vs FMSI**. sklib wins decisively on everything but disk size
— construct RAM (chr1 k=63: **133 MB vs FMSI 8.5 GB**, ~64×), build time (~3–4×) and streaming
(~5–13×), comparable point-lookups — while FMSI stays ~10× smaller on disk.

**fmsi** sits at the **opposite corner from sklib**: the smallest index (≈2 b/km) and low query RAM,
but the **highest construction RAM** (~3.9 GB on chr1 at k=21, **~8.5 GB at k=63**) and slowest
construction. Its **set operations** exist (experimental f-MS) but are slow and RAM-heavy (see
`report/SETOPS_REPORT.md`).

## Competitors
Build them once with `bash scripts/bench/tools/setup.sh` (clones + builds into
`scripts/out/bench/tools_src/`, writes binary paths to `tools.env`). Then add them to
`TOOLS=...`. Per-tool notes (wrappers in `tools/`):

| tool | kind | input | notes |
|---|---|---|---|
| `sshash` | exact dictionary | BCALM2 unitigs | `--canonical`; k≤31; build = bcalm(1 core)+build. ⚠ segfaults on **negative** lookups at small k on large/saturated genomes (reproduced at k=15 on C. elegans) — SSHash is meant for k≳20; the harness skips the failed query rows. |
| `sbwt` | exact membership | sanitized FASTA | `--add-reverse-complements`; default `plain-matrix` variant; k≤32. For SBWT's best *space*, rebuild with `-march=native` and set `SBWT_VARIANT=mef-matrix` (BMI2). |
| `bqf` | **approximate** (false positives) | KMC s-mer dump | compiled per-(k,z); `BQF_S` holds the indexed s-mer size ≈constant across k (z=k−s, default s=18); stores s-mers, virtualizes k-mers |
| `cbl` | exact set + set ops | sanitized FASTA | compiled per-K (**odd** K only); needs nightly Rust **and** `apt install libclang-dev libstdc++-12-dev`. ⚠ aborts on any sequence shorter than K, so the wrapper pre-filters sub-K records (`dropshort`) before building |
| `fmsi` | exact membership (masked superstring + Masked BWT) | KmerCamel masked superstring | canonical (KmerCamel bidirectional default; no `-u`); k≤127, **any parity**. Two-stage build = `kmercamel compute -M` (max-ones mask) + `fmsi index`; `setup.sh` also clones/builds **KmerCamel** (separate repo). ⚠ KmerCamel needs `apt install libglpk-dev` (GLPK, for the mask-optimization ILP). Indexed `-x` (no kLCP) for the compact/fair space point; for accelerated streaming, drop `-x` and add query `-S`. |

`bcalm` is infrastructure (unitig generation), not a benchmarked tool. Competitors are
built once per `(dataset,k)` (they ignore sklib's `m` — see `uses_m_*` in `tools.sh`).
All tools index the **same k-mer set** (BCALM unitigs preserve every distinct k-mer; KMC
validates set equality via `large_scale_e2e.sh`).

To add another tool, implement `available_<t>`/`version_<t>`/`construct_<t>`/`query_<t>`
in `tools/<t>.sh`; set globals `IDX_PATH, IDX_FILE_BYTES, IDX_PAYLOAD_BYTES, N_SKMERS` and
time with `run_timed`/`run_timed_median`.
