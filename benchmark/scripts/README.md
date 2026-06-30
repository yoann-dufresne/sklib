# Benchmark harness — scripts

Four experiment drivers, all sharing `lib.sh` (grids, caches, resumability, timing) and
`tools.sh` (per-tool adapters). See [`../README.md`](../README.md) for the conceptual
overview (experiments, metrics, capability matrix, datasets). This file is the operational
detail: how to run each script and what knobs exist.

## Prerequisites

- `kmc`, `kmc_tools`, `/usr/bin/time`, `python3`, `curl` on `PATH`. (On hosts without
  `/usr/bin/time` — e.g. RHEL — set `TIME_BIN` to a GNU `time`, such as `$CONDA_PREFIX/bin/time`;
  it must be GNU time, since peak RSS needs `-v`.)
- Genomes: catalogued in [`../data/genomes.tsv`](../data/genomes.tsv); fetched on demand by
  the drivers (via `genomes.sh`) or up-front with `bash fetch_genomes.sh [--list|--force] [names…]`.
- A **Release** `sskm` at `build-bench/bin/sskm` (override with `SSKM_BIN=…`). A DEBUG build
  enables AddressSanitizer and skews RSS/time, so do not point the harness at a debug binary.
  `cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release && cmake --build build-bench -j --target sskm`
- Competitors (optional): `bash tools/setup.sh` clones+builds sshash/sbwt/sbwtrs/cbl/bqf/fmsi into
  `../data/tools_src/` and writes `tools.env`. Unbuilt tools are skipped automatically.
- Figures: `python3 -m venv .venv && .venv/bin/pip install pandas matplotlib`.

## The four scripts

```bash
bash construct.sh        # exp 1 — build per (dataset,tool,k/m,threads)
bash query_single.sh     # exp 2 — individual k-mer queries over PRESENCE
bash query_stream.sh     # exp 3 — sequence queries over PRESENCE
bash setop.sh            # exp 4 — set ops over JACCARD (unitary + joint, materialize + size)
```

Each loops `dataset → tool (capability-gated) → k,m → threads`, builds/reuses a cached
index, and appends one row per measurement to its CSV in `../results/latest/`. Run them in
any order; query/set-op scripts build their indexes on demand (cached), so they are
self-contained.

## Env knobs (defaults)

Shared (in `lib.sh`): `DATASETS="sarscov2 ecoli yeast celegans chr21 chr1"` ·
`TOOLS="sklib sshash sbwt sbwtrs cbl bqf fmsi kmc"` · `KM="15,7 21,11 31,15 41,19 51,25 63,31"` ·
`THREADS="1 2 4 8 16"` · `PRESENCE="0 25 50 75 100"` · `JACCARD="0 0.1 0.3 0.5 0.7 0.9 1.0"` ·
`REPS=3` · `SEED=1234`. Query: `N_QUERY=200000` (single), `STREAM_RECS=2000`,
`STREAM_LEN=300` (stream). Paths: `SSKM_BIN`, `RESULTS` (default `../results/latest`),
`IDX_CACHE`, `QCACHE`, `GMUT`, `CSV` (per-script output override).

Override any to scope a run:
```bash
DATASETS=ecoli KM="21,11 31,15" THREADS="1 8" REPS=5 bash construct.sh
DATASETS=ecoli KM="21,11" PRESENCE="0 50 100" bash query_single.sh
DATASETS=ecoli KM="21,11" JACCARD="0.1 0.5 0.9" THREADS="1 8" bash setop.sh
```

## Resumability

CSVs are append-only. On start each script calls `load_done`, which keys the rows already
present by `(tool, tool_version, host, …)`; `is_done` then skips matching measurements. So
an interrupted run restarts without recomputing — and, on the **same tool version + same
host**, existing rows are taken as authoritative (no re-measure). A new build (version) or a
different machine yields a different key and is re-measured. Indexes (`IDX_CACHE`, keyed by
`tool/version/tag.k.m`), query FASTAs (`QCACHE`) and per-Jaccard mutants (`GMUT`) are cached
and regenerated only when missing. To force a full redo, delete the relevant CSV (and caches).

## Outputs (CSV schemas)

All rows are prefixed `timestamp,host,cpu`. `mrate` = Mk-mer/s, `bits_per_kmer` on payload.

| File | Key columns | Metric columns |
|---|---|---|
| `construct.csv` | `tool,tool_version,dataset,k,m,threads` | `time_s,peak_rss_kb,index_bytes,bits_per_kmer,n_kmers,n_superkmers,throughput_Mkmer_s` |
| `query_single.csv` / `query_stream.csv` | `…,threads,presence` | `n_kmers,present_kmers,absent_kmers,time_s,peak_rss_kb,throughput_Mkmer_s` |
| `setop.csv` | `…,threads,op,mode,jaccard_target` | `joint,jaccard_measured,result_kmers,time_s,peak_rss_kb,throughput_Mkmer_s` |

`op ∈ inter|union|diffab|diffba|joint` · `mode ∈ materialize|size`. Present-throughput is the
`presence=100` row, absent-throughput the `presence=0` row.

## Adding a tool

Implement `available_<t>` / `version_<t>` / `construct_<t> san k m wd` (set `IDX_PATH`,
`IDX_FILE_BYTES`, `IDX_PAYLOAD_BYTES`, `N_SKMERS`) / `query_<t> idx qfa k m cpus` in
`tools/<t>.sh`, then in `tools.sh` set its capability gates (`can_*`), `thread_list_<t>`,
and — for set ops — `setop_op_<t>` / `setop_joint_<t>` plus `setop_has_{size,joint}_<t>`.
Per-`k` limits live in the wrapper (e.g. `sbwt.sh` rejects `k>32`): just `return 1`.

## Beyond the four experiments — subfolders

The four drivers + their shared libs (`lib.sh`, `tools.sh`, `genomes.sh`), the genome fetcher
(`fetch_genomes.sh`), the query/mutant helpers (`e2e_helpers.py`, `mutate.py`) and the figure
generator (`plot.py`) sit at the top level. Everything else is grouped by role:

| Folder | What | Entry points |
|---|---|---|
| `tools/` | competitor adapters + one-shot setup | `setup.sh`, `{sshash,sbwt,sbwtrs,cbl,bqf,fmsi}.sh` (`bcalm.sh` = sshash preprocessor, not a competitor) |
| `producer/` | isolated super-k-mer producer throughput + bit-exact digest gate | `producer_{bench,median,setcheck,perf}.sh` |
| `profiling/` | construct/set-op perf attribution | `flamegraph_construct.sh`, `diag_perf.sh`, `construct_scaling.sh`, `diag_plot.py`; `bottleneck/` = set-op `_size`-vs-materialize decomposition + `perf` categorisation |
| `microbench/` | isolated single-thread set-op microbench + a focused KMC-vs-sklib joint compare | `union_bench.sh` / `union_ab.sh` (see [`../union_bench/README.md`](../union_bench/README.md)), `setop_joint_compare.sh` + `setop_joint_report.py` |
| `verify/` | correctness vs the KMC oracle (set-equality, determinism) | `large_scale_e2e.sh`, `greedy_chaining_verif.sh` |

- `plot.py [results_dir]` → per-experiment figures in `results/latest/figs/` (needs `.venv`).
- `verify/large_scale_e2e.sh` is correctness, not a benchmark; `tests/setop_multi_verif.sh`
  cross-validates combined set ops vs KMC.
