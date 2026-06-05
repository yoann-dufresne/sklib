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
- `.venv/`   — local Python env for plotting (pandas + matplotlib); git-ignored.

## Prerequisites
- `kmc`, `kmc_tools` on PATH (the k-mer-set oracle / bits-per-k-mer denominator).
- A built `sskm` at `build/bin/sskm` (build **Release** for representative numbers). sklib is
  single-threaded, so `THREADS>1` only pins CPU affinity — it does not change throughput.
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
the **payload** (skmer records), excluding the file header; for sklib the `VSKMER_3` header
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

Single machine (`Intel Core Ultra 7 165H`), **single-thread**, sklib **v0.3.0** (bucketed
`VSKMER_3` format, default `--buckets 4096`). Each cell is the full **FASTA → index** cost:
*sshash* includes BCALM2 unitig generation, *bqf* includes KMC s-mer counting. **Query** =
`positive` workload (scattered, all-present), **1 thread**, peak RSS during the lookup. `bqf`
is **approximate** (false positives); `cbl` is exact but **omitted on chr1** (its ~0.02
Mkmer/s query would take hours, and it OOMs building chr1). Regenerate from `results.csv`
with `python3 scripts/bench/compare.py` (tables) or `plots.py` (figures).

Columns: **c.t** = construct time (s) · **c.RAM** = construct peak RSS (MB) · **idx** =
index size on disk (MB) · **b/km** = bits per k-mer (payload) · **q.t** = query time (s) ·
**q.RAM** = query peak RSS (MB).

### k=15, m=7
| dataset | tool | c.t | c.RAM | idx (MB) | b/km | q.t | q.RAM |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli (4.6 Mbp) | **sklib** | 1.1 | 29 | 46.4 | 83.0 | 0.20 | 49 |
|  | sshash | 1.6 | 44 | 9.5 | 17.0 | 0.06 | 13 |
|  | sbwt | 2.5 | 91 | 6.9 | 12.4 | 0.09 | 36 |
|  | bqf | 1.5 | 74 | 13.6 | 24.4 | 0.09 | 17 |
|  | cbl | 1.2 | 150 | 21.0 | 37.6 | 0.22 | 135 |
| yeast (12 Mbp) | **sklib** | 3.1 | 27 | 125.2 | 91.2 | 0.27 | 124 |
|  | sshash | 8.2 | 152 | 36.9 | 26.9 | 0.11 | 39 |
|  | sbwt | 6.0 | 167 | 15.5 | 11.3 | 0.11 | 44 |
|  | bqf | 3.2 | 139 | 25.2 | 18.3 | 0.10 | 28 |
|  | cbl | 3.0 | 164 | 41.8 | 30.5 | 0.32 | 152 |
| chr21 (40 Mbp) | **sklib** | 9.6 | 63 | 338.6 | 99.6 | 0.43 | 328 |
|  | sshash | 18.6 | 431 | 151.9 | 44.7 | 0.18 | 149 |
|  | sbwt | 13.1 | 441 | 36.7 | 10.8 | 0.13 | 64 |
|  | bqf | 8.3 | 372 | 46.1 | 13.6 | 0.11 | 48 |
|  | cbl | 9.5 | 243 | 91.3 | 26.9 | 0.49 | 184 |
| celegans (100 Mbp) | **sklib** | 25.5 | 72 | 788.8 | 104.1 | 0.77 | 756 |
|  | sshash | 62.1 | 1013 | 428.3 | 56.5 | 0.33 | 413 |
|  | sbwt | 33.0 | 1016 | 80.6 | 10.6 | 0.20 | 106 |
|  | bqf | 18.6 | 846 | 83.9 | 11.1 | 0.14 | 84 |
|  | cbl | 29.0 | 360 | 191.8 | 25.3 | 0.90 | 247 |
| chr1 (230 Mbp) | **sklib** | 54.9 | 208 | 1548.9 | 110.8 | 1.26 | 1480 |
|  | sshash | 105.8 | 1654 | 1016.6 | 72.7 | 0.59 | 974 |
|  | sbwt | 54.6 | 1862 | 147.8 | 10.6 | 0.23 | 170 |
|  | bqf | 30.5 | 1920 | 83.9 | 6.0 | 0.13 | 84 |

### k=21, m=11
| dataset | tool | c.t | c.RAM | idx (MB) | b/km | q.t | q.RAM |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli (4.6 Mbp) | **sklib** | 1.0 | 27 | 21.2 | 37.2 | 0.09 | 25 |
|  | sshash | 0.8 | 54 | 5.8 | 10.2 | 0.02 | 10 |
|  | sbwt | 2.7 | 102 | 7.0 | 12.3 | 0.09 | 36 |
|  | bqf | 1.6 | 72 | 22.0 | 38.8 | 0.07 | 25 |
|  | cbl | 1.2 | 145 | 24.7 | 43.4 | 0.19 | 129 |
| yeast (12 Mbp) | **sklib** | 2.6 | 25 | 59.2 | 41.1 | 0.20 | 61 |
|  | sshash | 1.7 | 93 | 15.7 | 10.9 | 0.06 | 19 |
|  | sbwt | 6.6 | 198 | 16.1 | 11.2 | 0.14 | 44 |
|  | bqf | 4.4 | 131 | 41.9 | 29.2 | 0.15 | 44 |
|  | cbl | 2.9 | 170 | 53.7 | 37.4 | 0.35 | 156 |
| chr21 (40 Mbp) | **sklib** | 8.2 | 61 | 195.2 | 47.7 | 0.28 | 191 |
|  | sshash | 5.5 | 205 | 53.1 | 13.0 | 0.11 | 55 |
|  | sbwt | 15.7 | 428 | 44.0 | 10.8 | 0.23 | 71 |
|  | bqf | 12.3 | 342 | 79.7 | 19.5 | 0.17 | 80 |
|  | cbl | 9.3 | 307 | 140.0 | 34.2 | 0.58 | 236 |
| celegans (100 Mbp) | **sklib** | 23.6 | 70 | 628.1 | 55.1 | 0.48 | 603 |
|  | sshash | 17.8 | 524 | 159.1 | 14.0 | 0.18 | 156 |
|  | sbwt | 49.1 | 1031 | 120.8 | 10.6 | 0.32 | 144 |
|  | bqf | 34.5 | 794 | 285.2 | 25.0 | 0.28 | 276 |
|  | cbl | 34.8 | 1154 | 363.7 | 31.9 | 1.94 | 973 |
| chr1 (230 Mbp) | **sklib** | 50.4 | 206 | 1399.0 | 59.0 | 0.80 | 1336 |
|  | sshash | 35.4 | 1069 | 363.6 | 15.3 | 0.30 | 351 |
|  | sbwt | 88.0 | 1926 | 250.2 | 10.5 | 0.38 | 267 |
|  | bqf | 66.7 | 1810 | 536.9 | 22.6 | 0.37 | 516 |

### k=31, m=15
| dataset | tool | c.t | c.RAM | idx (MB) | b/km | q.t | q.RAM |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli (4.6 Mbp) | **sklib** | 0.9 | 22 | 12.3 | 21.4 | 0.19 | 16 |
|  | sshash | 0.5 | 36 | 4.1 | 7.1 | 0.05 | 8 |
|  | sbwt | 2.9 | 182 | 7.0 | 12.3 | 0.16 | 36 |
|  | bqf | 1.5 | 71 | 22.0 | 38.7 | 0.21 | 25 |
|  | cbl | 1.0 | 152 | 36.9 | 64.9 | 0.24 | 130 |
| yeast (12 Mbp) | **sklib** | 2.4 | 19 | 31.2 | 21.6 | 0.21 | 34 |
|  | sshash | 1.4 | 83 | 10.6 | 7.3 | 0.07 | 14 |
|  | sbwt | 7.0 | 408 | 16.2 | 11.2 | 0.19 | 44 |
|  | bqf | 4.1 | 130 | 41.9 | 29.0 | 0.29 | 44 |
|  | cbl | 2.8 | 211 | 87.0 | 60.2 | 0.37 | 182 |
| chr21 (40 Mbp) | **sklib** | 8.0 | 55 | 98.9 | 22.8 | 0.26 | 99 |
|  | sshash | 3.5 | 136 | 38.0 | 8.8 | 0.11 | 40 |
|  | sbwt | 18.4 | 1255 | 46.5 | 10.7 | 0.31 | 73 |
|  | bqf | 12.5 | 342 | 79.7 | 18.4 | 0.32 | 80 |
|  | cbl | 9.6 | 618 | 250.0 | 57.8 | 0.79 | 485 |
| celegans (100 Mbp) | **sklib** | 21.3 | 64 | 276.7 | 23.5 | 0.33 | 268 |
|  | sshash | 11.2 | 369 | 97.9 | 8.3 | 0.15 | 97 |
|  | sbwt | 54.6 | 1948 | 124.4 | 10.6 | 0.46 | 148 |
|  | bqf | 34.5 | 794 | 285.2 | 24.3 | 0.48 | 276 |
|  | cbl | 46.6 | 4976 | 673.8 | 57.3 | 5.68 | 4264 |
| chr1 (230 Mbp) | **sklib** | 49.5 | 199 | 619.5 | 24.3 | 0.52 | 594 |
|  | sshash | 26.5 | 735 | 246.7 | 9.7 | 0.23 | 240 |
|  | sbwt | 103.1 | 1948 | 269.0 | 10.5 | 0.55 | 285 |
|  | bqf | 66.9 | 1811 | 536.9 | 21.0 | 0.56 | 516 |

**Reading.** sklib keeps the **lowest construction RAM in every cell** (≈3–9× below the best
competitor on the big genomes — e.g. chr1 `k=21`: 206 MB vs 1069–1926 MB), and construction
time is mid-pack (faster than sbwt/bqf, near sshash). The trade-offs are **index size** and
**query**: as a sorted super-k-mer list (24 B/super-k-mer) it is far less compact than the
exact dictionaries (sbwt ≈10–12 b/km, sshash 7–17 b/km) or the approximate filter (bqf),
though b/km improves sharply as `k−m` grows (longer super-k-mers → better compaction); and its
query — now a **per-bucket dichotomic search** (no longer a linear scan) — still trails
sshash/bqf in throughput and uses the most query RAM, because genome-wide `positive`/`random`
workloads touch every bucket and so load most of the index (localized queries load only the
buckets they hit).

## Competitors
Build them once with `bash scripts/bench/tools/setup.sh` (clones + builds into
`scripts/out/bench/tools_src/`, writes binary paths to `tools.env`). Then add them to
`TOOLS=...`. Per-tool notes (wrappers in `tools/`):

| tool | kind | input | notes |
|---|---|---|---|
| `sshash` | exact dictionary | BCALM2 unitigs | `--canonical`; k≤31; build = bcalm(1 core)+build. ⚠ segfaults on **negative** lookups at small k on large/saturated genomes (reproduced at k=15 on C. elegans) — SSHash is meant for k≳20; the harness skips the failed query rows. |
| `sbwt` | exact membership | sanitized FASTA | `--add-reverse-complements`; default `plain-matrix` variant; k≤32. For SBWT's best *space*, rebuild with `-march=native` and set `SBWT_VARIANT=mef-matrix` (BMI2). |
| `bqf` | **approximate** (false positives) | KMC s-mer dump | compiled per-(k,z); `BQF_S` holds the indexed s-mer size ≈constant across k (z=k−s, default s=18); stores s-mers, virtualizes k-mers |
| `cbl` | exact set + set ops | sanitized FASTA | compiled per-K (**odd** K only); needs nightly Rust **and** `apt install libclang-dev libstdc++-12-dev` |

`bcalm` is infrastructure (unitig generation), not a benchmarked tool. Competitors are
built once per `(dataset,k)` (they ignore sklib's `m` — see `uses_m_*` in `tools.sh`).
All tools index the **same k-mer set** (BCALM unitigs preserve every distinct k-mer; KMC
validates set equality via `large_scale_e2e.sh`).

To add another tool, implement `available_<t>`/`version_<t>`/`construct_<t>`/`query_<t>`
in `tools/<t>.sh`; set globals `IDX_PATH, IDX_FILE_BYTES, IDX_PAYLOAD_BYTES, N_SKMERS` and
time with `run_timed`/`run_timed_median`.
