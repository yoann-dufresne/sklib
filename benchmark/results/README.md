# Benchmark results

This directory holds every committed benchmark artifact. There are three committed kinds, each
in its own subdirectory; two more directories are local scratch and **git-ignored**.

| Directory | Committed? | What it is |
|---|:--:|---|
| [`runs/`](runs/) | ✅ | **Full cross-tool campaigns** — the publishable, head-to-head comparisons of sklib vs the competitor tools. Each run is a self-contained dated folder (`EXPERIMENT.md` + `RESULTS.md` + `data/*.csv`). [`runs/README.md`](runs/README.md) is the **index of all runs**; start there. |
| [`journals/`](journals/) | ✅ | **Per-feature optimization write-ups** — internal dev journals tracking how one feature (set ops, construction, the producer, …) was profiled and sped up, with before/after numbers. See [`journals/README.md`](journals/README.md). |
| [`reference/`](reference/) | ✅ | **Frozen golden regression baselines** — curated CSV snapshots + figures used to catch regressions and to back the journals. *Not* a campaign. See [`reference/README.md`](reference/README.md). |
| `latest/` | ❌ git-ignored | The scratch the harness writes to: fresh CSVs, built indexes, query sets, logs, caches. Regenerable; safe to delete to reclaim disk. |
| `union_bench/` | ❌ git-ignored | Scratch for the set-op microbenchmark (`../../scripts/microbench/union_bench.sh`): built input lists + per-iteration TSVs. Regenerable. |

## Which numbers should I cite?

- For a **paper / external comparison** → the newest **reliable** row in [`runs/README.md`](runs/README.md)
  (currently [`runs/full_run_2026-06/`](runs/full_run_2026-06/)).
- For **"did this commit regress anything?"** → the matching CSV in [`reference/`](reference/).
- For **"how / why was feature X optimized?"** → the relevant journal in [`journals/`](journals/).

## Promoting a run

The harness writes to the git-ignored `latest/`. To publish a campaign, copy its
`{construct,query_single,query_stream,setop}.csv` into a new `runs/<YYYY-MM>[-tag]/data/`, add an
`EXPERIMENT.md` (machine, window, tools, grid) and `RESULTS.md` (tables/analysis), append a row to
[`runs/README.md`](runs/README.md), and commit. Cluster (Zeus) campaigns rsync back to the
git-ignored `benchmark_cluster/results/` first, then a curated snapshot is promoted here the same
way. A regression baseline is promoted to `reference/` instead.
