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
- `tools.sh` — per-tool adapters (`construct_*`, `query_*`, `version_*`, `available_*`).
  Only **sklib** is implemented; competitors are stubs until `tools/` wrappers land.
- `plots.py` — reads the CSV, writes Fig 1–6 to `scripts/out/bench/figures/`.
- `.venv/`   — local Python env for plotting (pandas + matplotlib); git-ignored.

## Prerequisites
- `kmc`, `kmc_tools` on PATH (the k-mer-set oracle / bits-per-k-mer denominator).
- A built `sskm` at `build/bin/sskm`. **For query thread-scaling, build with TBB**
  (`sudo apt install libtbb-dev`, then reconfigure+rebuild Release) — otherwise
  `std::execution::par` runs serially and the threads axis is flat.
- Plotting env (one-off):
  ```bash
  python3 -m venv scripts/bench/.venv
  scripts/bench/.venv/bin/pip install matplotlib pandas
  ```

## Run
```bash
# defaults: sklib on sarscov2+ecoli, k=21 m=11, all workloads, 1 thread
bash scripts/bench/bench.sh

# a real sweep
DATASETS="ecoli yeast celegans" \
KM="21,7 21,9 21,11 21,13 21,15 31,11 31,15" \
WORKLOADS="positive random streaming shuffled reads" \
THREADS="1 4 22" REPS=3 FRESH=1 \
bash scripts/bench/bench.sh

scripts/bench/.venv/bin/python scripts/bench/plots.py
```

### Env knobs (defaults)
`DATASETS="sarscov2 ecoli"` · `KM="21,11"` · `TOOLS="sklib"` ·
`WORKLOADS="positive random streaming shuffled reads"` · `THREADS="1"` ·
`N_QUERY=200000` · `STREAM_BP=2000000` · `READS_N=20000` · `READLEN=150` · `ERR=0.01` ·
`REPS=1` · `SEED=1234` · `SELFCHECK=1` · `FRESH=0` (set `FRESH=1` to truncate the CSV).

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
`positive` query rows (false-negative gate), `NA` elsewhere.

## Figures
1. `fig1_bits_per_kmer` — space: bits/k-mer vs m (per k, per dataset) + 2-bit-raw reference.
2. `fig2_compaction` — k-mers/super-k-mer vs m; #super-k-mers vs #k-mers.
3. `fig3_construction` — construct throughput and peak RSS vs dataset size.
4. `fig4_query` — throughput by workload + thread scaling (needs TBB + multiple `THREADS`).
5. `fig5_pareto` — bits/k-mer vs query throughput (sklib swept over m, competitors as points).
6. `fig6_envelope` — query correctness over the (k,m) grid.
7. `fig7_ksweep` — space (bits/k-mer) and streaming query speed vs k, one line per tool per dataset.

## Cross-tool comparison (current run)

Single machine (`Intel Core Ultra 7 165H`), **single-thread construction**, sklib
`v0.1.0` (`7a1051d`, low-memory bucketed path, default `--buckets 4096`). Each cell is the
full **FASTA → index** cost: *sshash* includes BCALM2 unitig generation, *bqf* includes
KMC s-mer counting. **Query** = `positive` workload (scattered, all-present), **1 thread**,
peak RSS during the lookup. Query is measured on the grid only — chr21/chr1 are
**construct-only** here (`—`) because the query path is being reworked; `bqf` is
**approximate** (false positives); `cbl` **fails to build on chr1** (out of memory).
Regenerate from `results.csv` with `plots.py`.

Columns: **c.t** = construct time (s) · **c.RAM** = construct peak RSS (MB) · **idx** =
index size on disk (MB) · **b/km** = bits per k-mer · **q.t** = query time (s) · **q.RAM**
= query peak RSS (MB).

### k=15, m=7
| dataset | tool | c.t | c.RAM | idx (MB) | b/km | q.t | q.RAM |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli (4.6 Mbp) | **sklib** | 1.6 | **23** | 50.8 | 91.1 | 0.8 | 54 |
| | sshash | 1.6 | 44 | 9.5 | 17.0 | 0.1 | 13 |
| | sbwt | 2.5 | 91 | 6.9 | 12.4 | 0.1 | 36 |
| | bqf | 1.5 | 74 | 13.6 | 24.4 | 0.1 | 17 |
| | cbl | 1.2 | 150 | 21.0 | 37.6 | 0.2 | 135 |
| yeast (12 Mbp) | **sklib** | 4.8 | **20** | 137.2 | 100.0 | 1.4 | 137 |
| | sshash | 8.2 | 152 | 36.9 | 26.9 | 0.1 | 39 |
| | sbwt | 6.0 | 167 | 15.5 | 11.3 | 0.1 | 44 |
| | bqf | 3.2 | 139 | 25.2 | 18.3 | 0.1 | 28 |
| | cbl | 3.0 | 164 | 41.8 | 30.5 | 0.3 | 152 |
| chr21 (40 Mbp) | **sklib** | 10.8 | **51** | 370.2 | 108.9 | — | — |
| | sshash | 21.3 | 432 | 151.9 | 44.7 | — | — |
| | sbwt | 15.3 | 442 | 36.7 | 10.8 | — | — |
| | bqf | 8.7 | 372 | 46.1 | 13.6 | — | — |
| | cbl | 9.8 | 243 | 91.3 | 26.9 | — | — |
| celegans (100 Mbp) | **sklib** | 31.9 | **147** | 869.1 | 114.7 | 14.6 | 834 |
| | sshash | 62.1 | 1013 | 428.3 | 56.5 | 0.3 | 413 |
| | sbwt | 33.0 | 1016 | 80.6 | 10.6 | 0.2 | 106 |
| | bqf | 18.6 | 846 | 83.9 | 11.1 | 0.1 | 84 |
| | cbl | 29.0 | 360 | 191.8 | 25.3 | 0.9 | 247 |
| chr1 (230 Mbp) | **sklib** | 62.2 | **192** | 1702.8 | 121.8 | — | — |
| | sshash | 120.2 | 1653 | 1016.6 | 72.7 | — | — |
| | sbwt | 70.7 | 1862 | 147.8 | 10.6 | — | — |
| | bqf | 49.1 | 1920 | 83.9 | 6.0 | — | — |
| | cbl | — | — | — | — | — | — |

### k=21, m=11
| dataset | tool | c.t | c.RAM | idx (MB) | b/km | q.t | q.RAM |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli (4.6 Mbp) | **sklib** | 1.4 | **21** | 24.6 | 43.4 | 0.3 | 29 |
| | sshash | 0.8 | 54 | 5.8 | 10.2 | 0.0 | 10 |
| | sbwt | 2.7 | 102 | 7.0 | 12.3 | 0.1 | 36 |
| | bqf | 1.6 | 72 | 22.0 | 38.8 | 0.1 | 25 |
| | cbl | 1.2 | 145 | 24.7 | 43.4 | 0.2 | 129 |
| yeast (12 Mbp) | **sklib** | 4.0 | **16** | 69.6 | 48.4 | 1.8 | 72 |
| | sshash | 1.7 | 93 | 15.7 | 10.9 | 0.1 | 19 |
| | sbwt | 6.6 | 198 | 16.1 | 11.2 | 0.1 | 44 |
| | bqf | 4.4 | 131 | 41.9 | 29.2 | 0.2 | 44 |
| | cbl | 2.9 | 170 | 53.7 | 37.4 | 0.4 | 156 |
| chr21 (40 Mbp) | **sklib** | 9.2 | **56** | 226.2 | 55.3 | — | — |
| | sshash | 6.1 | 205 | 53.1 | 13.0 | — | — |
| | sbwt | 18.1 | 428 | 44.0 | 10.8 | — | — |
| | bqf | 14.0 | 342 | 79.7 | 19.5 | — | — |
| | cbl | 10.3 | 307 | 140.0 | 34.2 | — | — |
| celegans (100 Mbp) | **sklib** | 28.8 | **341** | 732.9 | 64.3 | 47.3 | 705 |
| | sshash | 17.8 | 524 | 159.1 | 14.0 | 0.2 | 156 |
| | sbwt | 49.1 | 1031 | 120.8 | 10.6 | 0.3 | 144 |
| | bqf | 34.5 | 794 | 285.2 | 25.0 | 0.3 | 276 |
| | cbl | 34.8 | 1154 | 363.7 | 31.9 | 1.9 | 973 |
| chr1 (230 Mbp) | **sklib** | 57.0 | **261** | 1584.8 | 66.8 | — | — |
| | sshash | 36.6 | 1069 | 363.6 | 15.3 | — | — |
| | sbwt | 98.6 | 1926 | 250.2 | 10.5 | — | — |
| | bqf | 100.9 | 1811 | 536.9 | 22.6 | — | — |
| | cbl | — | — | — | — | — | — |

### k=31, m=15
| dataset | tool | c.t | c.RAM | idx (MB) | b/km | q.t | q.RAM |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli (4.6 Mbp) | **sklib** | 1.4 | **19** | 14.4 | 25.3 | 0.7 | 19 |
| | sshash | 0.5 | 36 | 4.1 | 7.1 | 0.1 | 8 |
| | sbwt | 2.9 | 182 | 7.0 | 12.3 | 0.2 | 36 |
| | bqf | 1.5 | 71 | 22.0 | 38.7 | 0.2 | 25 |
| | cbl | 1.0 | 152 | 36.9 | 64.9 | 0.2 | 130 |
| yeast (12 Mbp) | **sklib** | 3.4 | **14** | 37.6 | 26.0 | 0.9 | 41 |
| | sshash | 1.4 | 83 | 10.6 | 7.3 | 0.1 | 14 |
| | sbwt | 7.0 | 408 | 16.2 | 11.2 | 0.2 | 44 |
| | bqf | 4.1 | 130 | 41.9 | 29.0 | 0.3 | 44 |
| | cbl | 2.8 | 211 | 87.0 | 60.2 | 0.4 | 182 |
| chr21 (40 Mbp) | **sklib** | 9.2 | **53** | 119.6 | 27.6 | — | — |
| | sshash | 3.8 | 136 | 38.0 | 8.8 | — | — |
| | sbwt | 20.9 | 1256 | 46.5 | 10.7 | — | — |
| | bqf | 13.8 | 342 | 79.7 | 18.4 | — | — |
| | cbl | 10.7 | 618 | 250.0 | 57.8 | — | — |
| celegans (100 Mbp) | **sklib** | 27.0 | **224** | 352.3 | 30.0 | 33.7 | 342 |
| | sshash | 11.2 | 369 | 97.9 | 8.3 | 0.1 | 97 |
| | sbwt | 54.6 | 1948 | 124.4 | 10.6 | 0.5 | 148 |
| | bqf | 34.5 | 794 | 285.2 | 24.3 | 0.5 | 276 |
| | cbl | 46.6 | 4976 | 673.8 | 57.3 | 5.7 | 4264 |
| chr1 (230 Mbp) | **sklib** | 51.3 | **204** | 750.5 | 29.4 | — | — |
| | sshash | 28.3 | 734 | 246.7 | 9.7 | — | — |
| | sbwt | 114.8 | 1948 | 269.0 | 10.5 | — | — |
| | bqf | 71.5 | 1811 | 536.9 | 21.0 | — | — |
| | cbl | — | — | — | — | — | — |

**Reading.** sklib has the **lowest construction RAM in every cell** (1.5–8.6× below the
best competitor, up to ~30× below the worst), is **competitive on construction time**, and
is the only tool that stays under ~260 MB up to chr1 (cbl OOMs). The trade-off is **index
size**: as a sorted super-k-mer list (16 B/super-k-mer) it is far less compact than the
exact dictionaries (sbwt ≈10–12 b/km, sshash 7–17 b/km) or the approximate filter (bqf);
sklib's b/km improves sharply as `k−m` grows (longer super-k-mers → better compaction).
The O(N)-per-k-mer query (`find_closest_valid_skmer` linear scan) is the current bottleneck
and the next target.

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
