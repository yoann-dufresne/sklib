# Benchmark runs — index

Each row is one **full cross-tool campaign** (sklib vs competitors, all experiments), newest
first. A run is a self-contained folder: `EXPERIMENT.md` (machine, window, tools, parameter grid,
methodology), `RESULTS.md` (tables + analysis), and `data/*.csv` (raw measurements).

| Run | Window | Machine | sklib | Tools compared | Datasets | k grid | Status | Docs |
|---|---|---|---|---|---|---|:--:|---|
| [`full_run_2026-06/`](full_run_2026-06/) | 2026-06-08 → 06-16 | Precision-5490 · Core Ultra 7 165H · 22c / 62 GiB | 0.11.0 | sklib, KMC, CBL, sshash, SBWT (C++), sbwtrs (Rust), FMSI, BQF | ecoli, yeast, celegans, chr21, chr1 | 15·21·31·41·51·63 | ✅ **current** | [EXPERIMENT](full_run_2026-06/EXPERIMENT.md) · [RESULTS](full_run_2026-06/RESULTS.md) · [data/](full_run_2026-06/data/) |
| _cluster chm13 (Zeus)_ | _in progress_ | Lille HPC (Zeus), Slurm | 0.13.x | same set | + **chm13** (3 GB T2T) | 15…63 | 🚧 planned | — |

**Status legend.** ✅ current/reliable — cite these · 🚧 in progress · ⏳ superseded (kept for history).

## Notes

- `full_run_2026-06` is single-machine and **defers chm13** (3 GB needs the cluster). The
  in-progress Zeus campaign adds chm13 and re-measures sklib at the current release; when it lands
  it becomes a new `runs/<…>/` folder — fill in its row above and flip the statuses (the older run
  to ⏳ if superseded).
- The cluster overlay lives in the repo's sibling [`benchmark_cluster/`](../../../benchmark_cluster/);
  it reuses `benchmark/scripts/` unchanged and rsyncs results back to a git-ignored tree, from which
  a curated snapshot is promoted into a `runs/` folder (see [`../README.md`](../README.md)).
- Within `full_run_2026-06`, sklib appears at several `tool_version`s from successive
  re-measurements across the campaign window; filter `tool_version == "sklib-0.11.0"` for the
  canonical head-to-head (see that run's `EXPERIMENT.md` → *Versioning note*).
