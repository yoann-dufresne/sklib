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
| `TOOLS` | `sklib sshash sbwt cbl bqf fmsi kmc` | gated per experiment by `can_{construct,query,setop}_<tool>`. |
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
| **sklib** (`sskm`) | ✓ | ✓ | ✓ | ✓ (`--*-out`/`--sizes`) | ✓ | ✓ (`-t`) | ≤63 (uint128) |
| **KMC** | ✓ | — (oracle) | ✓ | ✓ (`kmc_tools simple`) | — | ✓ (`-t`) | ≤256 |
| sshash | ✓ | ✓ | — | — | — | — | ≤31 |
| sbwt | ✓ | ✓ | — | — | — | — | ≤32 |
| cbl | ✓ | ✓ | ✓ | — | — | — | **odd** ≤59 |
| bqf | ✓ | ✓ (approx) | — | — | — | — | per-(k,z) |
| fmsi | ✓ | ✓ | ✓ (experimental) | — | — | — | ≤127 |

So at `k>32` the membership field is sklib + fmsi (+ kmc for set ops); only sklib/KMC vary
with the thread sweep. KMC has no cardinality-only mode, so `mode=size` rows are sklib-only.

## Datasets

Catalogued genomes are named in `DATASETS=`; the harness downloads/sanitises on first use
into `data/genomes/<name>.sanitized.fa` (git-ignored). Custom datasets (pangenomes, rice,
metagenomes, real reads) and the per-Jaccard mutants are in [`data/README.md`](data/README.md).

| Name | Organism | ~Size | Source |
|---|---|---|---|
| `sarscov2` / `ecoli` | SARS-CoV-2 / E. coli | 30 kb / 4.6 Mb | local fixture |
| `yeast` / `celegans` | S. cerevisiae / C. elegans | 12 / 100 Mb | UCSC |
| `chr21` / `chr20` / `chr1` | human chromosome (hg38) | 40 / 64 / 230 Mb | UCSC |
| `chm13` | human T2T-CHM13v2 | 3.1 Gb | UCSC (hs1) |

## Layout

```
benchmark/
├── README.md                 # this index
├── scripts/
│   ├── README.md             # detailed usage, env knobs, CSV schemas
│   ├── construct.sh  query_single.sh  query_stream.sh  setop.sh   # the 4 experiments
│   ├── lib.sh  tools.sh      # shared harness (grids, caches, resumability) + tool adapters
│   ├── e2e_helpers.py  mutate.py  plot.py                          # query-gen, mutants, figures
│   ├── large_scale_e2e.sh    # correctness vs KMC (not a benchmark)
│   ├── flamegraph_construct.sh
│   ├── tools/                # sshash/sbwt/cbl/bqf/fmsi/bcalm wrappers + setup.sh
│   └── bottleneck/           # set-op time decomposition + perf categorisation
├── data/
│   ├── README.md             # dataset catalogue + download/prepare
│   ├── genomes/              # sanitise cache + mutants/ (git-ignored)
│   └── tools_src/            # built competitor binaries + tools.env (git-ignored)
└── results/
    ├── reports/              # analysis write-ups (committed)
    ├── reference/            # curated CSV + figure snapshots (committed)
    └── latest/               # fresh runs: CSVs, indexes/, queries/, figs/ (git-ignored)
```

**Results convention.** Scripts write to `results/latest/` (git-ignored). Promote a run by
copying the relevant CSVs/figures to `results/reference/` and committing.

## Quickstart

```bash
# 1. Build a Release sskm where the harness expects it (DEBUG/ASan skews RSS & time).
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release && cmake --build build-bench -j --target sskm

# 2. Run the four experiments (needs kmc, kmc_tools on PATH). Override grids to scope.
bash benchmark/scripts/construct.sh
bash benchmark/scripts/query_single.sh
bash benchmark/scripts/query_stream.sh
bash benchmark/scripts/setop.sh
#   e.g. a quick slice:  DATASETS=ecoli KM="21,11" THREADS="1 8" bash benchmark/scripts/setop.sh

# 3. Figures (one-off venv for pandas + matplotlib).
python3 -m venv benchmark/scripts/.venv && benchmark/scripts/.venv/bin/pip install pandas matplotlib
benchmark/scripts/.venv/bin/python benchmark/scripts/plot.py    # -> results/latest/figs/

# 4. Competitors (optional; clones + builds into data/tools_src/).
bash benchmark/scripts/tools/setup.sh
```

Detailed usage and env knobs: **[`scripts/README.md`](scripts/README.md)**.

## Reports

- [`results/reports/SETOPS_MULTI_REPORT.md`](results/reports/SETOPS_MULTI_REPORT.md) — combined single-pass vs sequential set ops (sklib & KMC).
- [`results/reports/SETOPS_REPORT.md`](results/reports/SETOPS_REPORT.md) — set ops: sklib vs KMC, CBL, FMSI (single- and multi-core), scaling, overlap, memory.
- [`results/reports/SETOPS_BOTTLENECKS.md`](results/reports/SETOPS_BOTTLENECKS.md) — where set-op time goes and the implemented speedups.
- [`results/reports/CONSTRUCT_SPEEDUP.md`](results/reports/CONSTRUCT_SPEEDUP.md) — parallel per-bucket construction (`-t`).
