# Optimization journals

Per-feature engineering write-ups: how one part of sklib was profiled and optimized, with
before/after numbers and the correctness gates used. These are **internal dev records**, not the
headline cross-tool comparison — for that, see [`../runs/`](../runs/). The raw data behind them
lives in [`../reference/`](../reference/) (plus the git-ignored microbench scratch).

### Set operations
- [`SETOPS_REPORT.md`](SETOPS_REPORT.md) — sklib vs KMC / CBL / FMSI on set ops: time, peak RAM, correctness; single- and multi-core, scaling, overlap, memory.
- [`SETOPS_BOTTLENECKS.md`](SETOPS_BOTTLENECKS.md) — `perf`/callgrind breakdown of set-op time (≈85–91 % is post-merge re-compaction) and the implemented speedups.
- [`SETOPS_MULTI_REPORT.md`](SETOPS_MULTI_REPORT.md) — combined single-pass set ops (`--*-out`/`--sizes`) vs the sequential single-op runs, sklib & KMC.
- [`SETOP_RECHAIN.md`](SETOP_RECHAIN.md) — what materializing a result costs: phase split (the re-chaining is ~60 % of the per-bucket wall at `-t1`), why `--no-compact` is *not* a shortcut (it loses in 395/400 configs), and the closure check `A ∪ A` = A's own record count.
- [`SETOPS_WIDE_WIDTH.md`](SETOPS_WIDE_WIDTH.md) — width-dispatched merge for the wide stores (`__uint128`/`kuint256`): cached masked k-mers, single 3-way compare.
- [`UNION_SPEEDUP.md`](UNION_SPEEDUP.md) — single-thread `set_union` optimization journal (driven by the `microbench/union_bench` harness).
- [`XOR_SPEEDUP.md`](XOR_SPEEDUP.md) — symmetric-difference (`xor`) tuning, incl. the identical-record drop fast path at high Jaccard.

### Construction & producer
- [`CONSTRUCT_SPEEDUP.md`](CONSTRUCT_SPEEDUP.md) — parallel per-bucket phase-2 compaction (`-t`); byte-identical index at any thread count.
- [`CONSTRUCT_SCALING_DIAG.md`](CONSTRUCT_SCALING_DIAG.md) — Amdahl analysis of construction scaling (the serial phase-1 producer caps speedup).
- [`PRODUCER_BOTTLENECKS.md`](PRODUCER_BOTTLENECKS.md) — flamegraph / `perf stat` of the isolated super-k-mer producer (instruction-throughput bound).
- [`PRODUCER_SPEEDUP.md`](PRODUCER_SPEEDUP.md) — six successive producer micro-optimizations (cumulative ~+58 %), digest-gated.
- [`GREEDY_DEFAULT.md`](GREEDY_DEFAULT.md) — switching the construction default to greedy chaining (identical record counts, ~3–7 % faster phase-2).

### Cross-cutting sweeps
- [`OPTIM_SWEEP_2026-07.md`](OPTIM_SWEEP_2026-07.md) — full from-scratch sweep over all three workloads (no prior journal read, everything re-measured): noise-floor protocol, time map per workload × width, 19 levers instructed (7 landed, 7 measured-and-rejected incl. PGO and record packing, 5 untested). construct +7…+49 %, query +18…+172 %, setop +1…+3,5 %; byte-identical output throughout. Also reports a pre-existing `construct` segfault for `2m − b ≥ 8·sizeof(gen)`.
