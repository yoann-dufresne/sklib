# sklib benchmarks

Everything used to benchmark **sklib**'s sorted super-k-mer list (`sskm`) lives here:
the run scripts, the input-data catalogue, the measurement outputs, and the analysis
reports. The harness is organised around **four experiments**, each with a fixed set of
minimal metrics, each run across a uniform thread sweep and the whole `k/m` range, on
**all capable tools × all data**.

## The four experiments

| # | Script | Measures | Minimal metrics | Output |
|---|---|---|---|---|
| 1 | [`construct.sh`](scripts/construct.sh) | index construction | time · peak RSS · index size · bits/k-mer · throughput | `results/latest/construct.csv` |
| 2 | [`query_single.sh`](scripts/query_single.sh) | individual k-mer queries (one k-mer/record), swept over present-fraction | time · peak RSS · present-throughput · absent-throughput | `results/latest/query_single.csv` |
| 3 | [`query_stream.sh`](scripts/query_stream.sh) | stream (sequence) queries, swept over present-fraction | time · peak RSS · present-throughput · absent-throughput | `results/latest/query_stream.csv` |
| 4 | [`setop.sh`](scripts/setop.sh) | set ops {∪, ∩, A∖B, B∖A}, **unitary and joint**, materialize + cardinality, swept over Jaccard | time · peak RSS · throughput per Jaccard | `results/latest/setop.csv` |

For queries, **present-throughput** is the `presence=100` row and **absent-throughput** the
`presence=0` row (the present/absent k-mer counts are exact, by construction). For set ops,
**joint** is the combined single-pass (`sskm setop --*-out / --sizes`; KMC `kmc_tools simple`
with several ops) compared against the **sequential** sum of the four unitary passes.

## Sweep policy (every experiment, all env-overridable)

| Axis | Default | Notes |
|---|---|---|
| `THREADS` | `1 2 4 8 16` | mono-thread **and** multi-thread. Tools without `-t` (everything but sklib/KMC) are measured once at `t=1` and shown flat — still compared. |
| `KM` | `15,7 21,11 31,15 41,19 51,25 63,31` | whole `k` range (`m≈k/2`), up to sklib's `2*(2k-m)≤256` ceiling. At large `k` the field narrows automatically (see capability matrix). |
| `PRESENCE` | `0 25 50 75 100` | present-fraction (%) of the query set. |
| `JACCARD` | `0 0.1 0.3 0.5 0.7 0.9 1.0` | target overlap; realised via `mutate.py` (`rate = 1-(2J/(1+J))^(1/k)`), measured J recorded. |
| `DATASETS` | `sarscov2 ecoli yeast celegans chr21 chr1` | all datasets run on all capable tools. |
| `TOOLS` | `sklib sshash sbwt sbwtrs cbl bqf fmsi kmc` | gated per experiment by `can_{construct,query,setop}_<tool>`. |
| `REPS` | `3` | median wall-time per point (max RSS). |

The full matrix is `DATASETS × TOOLS × (k,m) × THREADS × (PRESENCE | JACCARD × ops)`; scope a
run by overriding any grid, e.g. `DATASETS=ecoli KM="21,11" THREADS="1 8" bash …/setop.sh`.

## Resumability

The global scripts are **interruptible and restartable without redoing work.** Each CSV is
append-only; on (re)start a script loads the rows already present and **skips any
measurement whose `(tool, tool_version, host, …)` key is already there** — i.e. *same tool
version + same machine ⇒ existing values are trusted as final*. A different build (version)
or a different host re-measures automatically. Built indexes, generated query sets and
mutated genomes are likewise cached (under the git-ignored `results/latest/` and
`data/genomes/mutants/`) and regenerated only when missing. So an interrupted overnight run
just picks up where it stopped.

## Metrics & CSV schemas

Measured with `/usr/bin/time -v` (wall time median over `REPS`, peak RSS), `bits/k-mer` on
the payload (`8·payload_bytes / distinct_kmers`, excludes the `VSKMER_4` header/directory),
throughput in `Mk-mer/s`. Set-op result cardinalities are authoritative from one sklib
`--sizes` pass (cross-validated `== KMC` by `tests/setop_multi_verif.sh`).

```
construct.csv     …,tool,tool_version,dataset,k,m,threads,time_s,peak_rss_kb,index_bytes,bits_per_kmer,n_kmers,n_superkmers,throughput_Mkmer_s
query_*.csv       …,tool,tool_version,dataset,k,m,threads,presence,n_kmers,present_kmers,absent_kmers,time_s,peak_rss_kb,throughput_Mkmer_s
setop.csv         …,tool,tool_version,dataset,k,m,threads,op,mode,joint,jaccard_target,jaccard_measured,result_kmers,time_s,peak_rss_kb,throughput_Mkmer_s
```
(each prefixed `timestamp,host,cpu`). `op ∈ inter|union|diffab|diffba|joint`; `mode ∈ materialize|size`.

## Competing tools — per-experiment capability

Built once into `data/tools_src/` by [`scripts/tools/setup.sh`](scripts/tools/setup.sh)
(writes `tools.env`); unavailable tools are skipped gracefully. `construct_<tool>` returning
non-zero (e.g. an out-of-range `k`) also skips just that `(tool,k)`.

| Tool | Construct | Query | Set ops | Joint set ops | Card-only | Multi-thread | k range |
|---|:--:|:--:|:--:|:--:|:--:|:--:|---|
| **sklib** (`sskm`) | ✓ | ✓ | ✓ | ✓ (`--*-out`/`--sizes`) | ✓ | ✓ (`-t`) | ≤127 (kuint256) |
| **KMC** | ✓ | — (oracle) | ✓ | ✓ (`kmc_tools simple`) | — | ✓ (`-t`) | ≤256 |
| sshash | ✓ | ✓ | — | — | — | — | ≤31 |
| sbwt *(C++)* | ✓ | ✓ | — | — | — | — | ≤32 (`MAX_KMER_LENGTH`) |
| sbwtrs *(Rust)* | ✓ | ✓ | ✓ (materialize) | — | — | ✓ (`-t`, construct & setop) | ≤255 |
| cbl | ✓ | ✓ | ✓ | — | — | — | **odd** ≤59 |
| bqf | ✓ | ✓ (approx) | — | — | — | — | per-(k,z); **≤31** on default sweep (see note) |
| fmsi | ✓ | ✓ | ✓ (experimental) | — | — | — | ≤127 |

`sbwt` (C++, algbio) and `sbwtrs` (Rust, jnalanko/sbwt-rs) are the **same algorithm, two
implementations**, kept as separate rows so they never collide and the rewrite can be compared
head-to-head; only the Rust one lifts the k≤32 cap and adds set ops + a thread sweep. So at
`k>32` the membership field is sklib + sbwtrs + fmsi (+ kmc for set ops); sklib, KMC and sbwtrs
vary with the thread sweep. KMC has no cardinality-only mode, so `mode=size` rows are sklib-only.

### Expected per-tool skips on the default `k`-sweep (these are *not* errors)

On `KM` up to `63,31`, each competitor refuses the `k` beyond its hard limit: the driver
logs `construct failed/unsupported (skipped)` (or a wrapper's `build failed`) in
`results/latest/run.log` and moves on — that one `(tool,k)` cell is simply absent from the
CSV. **sklib, KMC and fmsi cover the whole grid** (k=15…63), and none of these skips involve
a crash/OOM. The skips re-appear (and fail fast) on every *resume*, since failed cells are
not memoised.

| Tool | OK on default sweep | Skipped (logged) | Why |
|---|---|---|---|
| sshash | 15, 21, 31 | 41, 51, 63 | 64-bit k-mers ⇒ ≤31 (default build); logs `sshash: build failed` |
| sbwt | 15, 21, 31 | 41, 51, 63 | compiled `MAX_KMER_LENGTH=32` ⇒ ≤32; logs `sbwt: k>32 unsupported` |
| cbl | 15, 21, 31, 41, 51 | 63 | `build.rs` panics `K must be ≤ 59` (odd-`k` only, ≤59) |
| bqf | 15, 21, 31 | 41, 51, 63 | the per-(k,z) **binary compiles fine**; `bqf build` (the filter) fails because the wrapper holds `s=BQF_S=18` fixed for cross-`k` comparability, so `z=k−18` (≥23 at k≥41) exceeds fimpera's virtualization limit. Raise `BQF_S` to extend `k` — at the cost of comparability across `k`. bqf is *approximate* regardless. |

## Datasets

Catalogued genomes are declared in [`data/genomes.tsv`](data/genomes.tsv) and named in
`DATASETS=`; the harness downloads/sanitises on first use into
`data/genomes/<name>.sanitized.fa` (git-ignored). Pre-fetch any subset with
`bash scripts/fetch_genomes.sh [--list|--force] [names…]`. Custom datasets (pangenomes,
rice, metagenomes, real reads) and the per-Jaccard mutants are in
[`data/README.md`](data/README.md).

| Name | Organism | ~Size | Source |
|---|---|---|---|
| `sarscov2` / `ecoli` | SARS-CoV-2 / E. coli | 30 kb / 4.6 Mb | NCBI (`NC_045512.2` / `NC_000913.3`) |
| `yeast` / `celegans` | S. cerevisiae / C. elegans | 12 / 100 Mb | UCSC |
| `chr21` / `chr20` / `chr1` | human chromosome (hg38) | 40 / 64 / 230 Mb | UCSC |
| `chm13` | human T2T-CHM13v2 | 3.1 Gb | UCSC (hs1) |

## Layout

```
benchmark/
├── README.md                 # this index
├── scripts/                 # see scripts/README.md
│   ├── README.md             # detailed usage, env knobs, CSV schemas
│   ├── construct.sh  query_single.sh  query_stream.sh  setop.sh   # the 4 experiments
│   ├── lib.sh  tools.sh      # shared harness (grids, caches, resumability) + tool adapters
│   ├── genomes.sh  fetch_genomes.sh   # genome catalogue loader + CLI downloader
│   ├── e2e_helpers.py  mutate.py  plot.py                          # query-gen, mutants, figures
│   ├── tools/                # competitor adapters: sshash/sbwt/sbwtrs/cbl/bqf + setup.sh
│   │                         #   (bcalm.sh = preprocessor for sshash, not a benchmarked tool)
│   ├── producer/             # isolated super-k-mer producer throughput + digest gate
│   ├── profiling/            # construct/set-op perf: flamegraphs, scaling, bottleneck/ decomposition
│   ├── microbench/           # isolated set-op microbench (union_bench) + KMC-vs-sklib joint compare
│   └── verify/               # correctness vs the KMC oracle (large_scale_e2e, greedy chaining)
├── data/
│   ├── README.md             # dataset catalogue + download/prepare
│   ├── genomes.tsv           # genome catalogue (single source of truth; fetch_genomes.sh reads it)
│   ├── genomes/              # sanitise cache + mutants/ (git-ignored)
│   └── tools_src/            # built competitor binaries + tools.env (git-ignored)
└── results/                 # see results/README.md for the full map
    ├── runs/                 # full cross-tool campaigns + RUNS INDEX (committed)
    ├── journals/             # per-feature optimization write-ups (committed)
    ├── reference/            # frozen golden regression baselines + figures (committed)
    └── latest/               # fresh runs: CSVs, indexes/, queries/, figs/ (git-ignored)
```

**Results convention.** Scripts write to `results/latest/` (git-ignored). Promote a finished
campaign by copying its CSVs + a short `EXPERIMENT.md`/`RESULTS.md` into a new
`results/runs/<date>/` (and add a row to [`results/runs/README.md`](results/runs/README.md)), or
promote a curated regression snapshot to `results/reference/`, then commit.

## Quickstart

```bash
# 1. Build a Release sskm where the harness expects it (DEBUG/ASan skews RSS & time).
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release && cmake --build build-bench -j --target sskm

# 2. Fetch the genomes (or a subset). Optional: the harness also fetches on first use.
bash benchmark/scripts/fetch_genomes.sh                 # all; or e.g.: fetch_genomes.sh sarscov2 ecoli

# 3. Run the four experiments (needs kmc, kmc_tools on PATH). Override grids to scope.
bash benchmark/scripts/construct.sh
bash benchmark/scripts/query_single.sh
bash benchmark/scripts/query_stream.sh
bash benchmark/scripts/setop.sh
#   e.g. a quick slice:  DATASETS=ecoli KM="21,11" THREADS="1 8" bash benchmark/scripts/setop.sh

# 4. Figures (one-off venv for pandas + matplotlib).
python3 -m venv benchmark/scripts/.venv && benchmark/scripts/.venv/bin/pip install pandas matplotlib
benchmark/scripts/.venv/bin/python benchmark/scripts/plot.py    # -> results/latest/figs/

# 5. Competitors (optional; clones + builds into data/tools_src/).
bash benchmark/scripts/tools/setup.sh
```

Detailed usage and env knobs: **[`scripts/README.md`](scripts/README.md)**.

## Results & reports

- **Latest full run** (the authoritative cross-tool comparison): [`results/runs/full_run_2026-06/`](results/runs/full_run_2026-06/) — see its [`RESULTS.md`](results/runs/full_run_2026-06/RESULTS.md) and `data/*.csv`.
- **All runs** (index of every campaign, newest first): [`results/runs/README.md`](results/runs/README.md).
- **Optimization journals** (per-feature dev write-ups in [`results/journals/`](results/journals/)) — full index in [`results/journals/README.md`](results/journals/README.md); highlights:
  - [`SETOPS_MULTI_REPORT.md`](results/journals/SETOPS_MULTI_REPORT.md) — combined single-pass vs sequential set ops (sklib & KMC).
  - [`SETOPS_REPORT.md`](results/journals/SETOPS_REPORT.md) — set ops: sklib vs KMC, CBL, FMSI (single- and multi-core), scaling, overlap, memory.
  - [`SETOPS_BOTTLENECKS.md`](results/journals/SETOPS_BOTTLENECKS.md) — where set-op time goes and the implemented speedups.
  - [`CONSTRUCT_SPEEDUP.md`](results/journals/CONSTRUCT_SPEEDUP.md) — parallel per-bucket construction (`-t`).
