# Construction multithread-scaling diagnostic

**Question.** `sskm construct -t` plateaus well below linear. Is a single-threaded section / producer
thread the cause, and how much does it cap the speedup — on real data?

**TL;DR.** Yes. Construction is two phases; **only Phase 2 (per-bucket compaction) is parallel**.
**Phase 1 — the producer (FASTA parse → rolling minimizer → bucket routing) — is strictly
single-threaded**, and its wall time is *identical at every `-t`*. That serial producer is **~31–37 %
of the total work** (the true serial fraction `s`, agreed by three independent methods), which Amdahl-
caps the speedup at **~2.7–3.2×**. The realized maximum is **2.15–2.59×** at `-t22` (rising with
genome size); the shortfall below the ceiling is a **secondary** effect — Phase 2 itself loses
efficiency past ~8 threads (memory-bandwidth / hybrid-core bound: it keeps ~20 cores busy but does
~2× the CPU work at `-t22`). **At the default `-t8`, Phase 1 already accounts for ~73–78 % of wall
time**, so the producer is the dominant floor.

## Method

- **Binary:** `build-timing/bin/sskm`, clang-18 `-O3 -march=native` (Release). Profile build
  (`-O3 -ggdb3 -fno-inline`) for perf.
- **Instrumentation:** an env-gated `SKLIB_TIMING` phase split added to `build_bucketed`
  (`lib/include/algorithms/SortedSkmerListBuilder.hpp`) — a handful of `steady_clock` stamps at the
  Phase-1/Phase-2 boundaries, **no per-iteration cost**, **no behaviour change** (index byte-identical
  with the env on/off and across `-t`, sha256-verified).
- **Data:** real genomes, fetched + sanitized via `benchmark/scripts/genomes.sh` — ecoli (4.6 Mbp),
  yeast (12.2 Mbp), chr21 (40.1 Mbp), celegans (100.3 Mbp). k=21, m=11, default `--buckets 4096`.
- **Sweep:** `benchmark/scripts/construct_scaling.sh` — `-t ∈ {1,2,3,4,6,8,12,16,22}`, 3 reps (median),
  1 warmup/genome, tmp on NVMe. Per run: wall + CPU% + peak RSS (`/usr/bin/time -v`) and
  `phase1_s`/`phase2_s` (`SKLIB_TIMING`). 108 rows, 0 failures.
- **Machine:** Intel Core Ultra 7 165H — 22 logical cores (6 P-cores/12 threads + 8 E + 2 LP-E),
  62 GB RAM, NVMe. Inputs fit the page cache ⇒ Phase 1 is CPU-bound, not disk-bound.
- **Three independent estimates of the serial fraction `s`:** (1) timer — `p1/(p1+p2)` at `-t1`;
  (2) Amdahl least-squares fit of `speedup(t)=1/(s+(1-s)/t)`; (3) CPU-time — `p1 / cpu_seconds@t1`.

## Results — the producer is the floor

**Phase 1 wall time is flat across `-t`; only Phase 2 shrinks** (median seconds):

| genome | `phase1_s` (t1 → t22) | `phase2_s` (t1 → t8 → t22) |
|---|---|---|
| ecoli | 0.365 → 0.380 (flat) | 0.612 → 0.102 → 0.070 |
| yeast | 0.946 → 0.976 (flat) | 1.779 → 0.273 → 0.175 |
| chr21 | 3.231 → 3.231 (flat) | 5.734 → 0.935 → 0.530 |
| celegans | 7.492 → 7.693 (flat) | 16.557 → 2.780 → 1.564 |

Phase 1 does not move when you add threads — it is single-threaded. Phase 2 falls ~1/t (with
declining efficiency). So `-t` only attacks the smaller half of the work.

**Wall time and speedup vs `-t`:**

| genome | t1 | t2 | t4 | **t8** | t16 | t22 | max speedup |
|---|--:|--:|--:|--:|--:|--:|--:|
| ecoli | 0.985 s | 1.43× | 1.81× | **2.06×** | 2.13× | 2.15× | 2.15× |
| yeast | 2.740 s | 1.46× | 1.91× | **2.20×** | 2.34× | 2.37× | 2.37× |
| chr21 | 8.973 s | 1.46× | 1.85× | **2.12×** | 2.32× | 2.38× | 2.38× |
| celegans | 24.06 s | 1.49× | 1.94× | **2.31×** | 2.55× | 2.59× | 2.59× |

Past `-t8`, ≤12 % more speed for 2.75× the threads (and growing RAM) — `-t8` is the knee.

**At the default `-t8`, Phase 1 is the majority of wall time:** ecoli 77 %, yeast 78 %, chr21 78 %,
celegans 73 %.

### Serial fraction and ceiling — three methods agree

| genome | s (timer) | s (CPU@t1) | s (Amdahl fit) | ceiling 1/s | **realized max** |
|---|--:|--:|--:|--:|--:|
| ecoli | 0.374 | 0.374 | 0.425 | 2.68× | 2.15× |
| yeast | 0.347 | 0.349 | 0.384 | 2.88× | 2.37× |
| chr21 | 0.360 | 0.364 | 0.392 | 2.77× | 2.38× |
| celegans | 0.312 | 0.315 | 0.353 | 3.21× | 2.59× |

Timer and CPU-derived `s` match to within 0.01 — the serial fraction is real and ~⅓ of the work. The
Amdahl-fit `s` is slightly higher because it also absorbs Phase-2's imperfect scaling (below). Bigger
genomes have a **smaller** serial fraction (Phase 2's sort + column work grows a touch faster than the
linear parse), so they scale better — celegans reaches 2.59× vs ecoli 2.15×.

### CPU% — the serial signature, and a secondary Phase-2 ceiling

Total CPU% (from `/usr/bin/time -v`) tops out far below the 2200 % a 22-core machine could give:
**99 % at `-t1`** (one core — serial), rising only to **365–431 % at `-t22`**. Most wall time is
spent with a single core busy (Phase 1).

But Phase 2 also stops scaling cleanly past ~8 threads. Phase-2 efficiency `p2(1)/(t·p2(t))` falls
0.74–0.81 at `-t8` → **0.40–0.49 at `-t22`**, and total CPU-work in Phase 2 **inflates ~1.9–2.2×**
from `-t1` to `-t22` while ~19–21 cores stay busy. That signature (cores occupied, work expands) is
**memory-bandwidth saturation + the hybrid P/E topology** (beyond ~12 threads the work spills onto
slower E-cores and HT siblings). So two effects compound: a hard serial floor (Phase 1) *and*
diminishing Phase-2 returns at high `-t`, which is why the realized speedup sits below even the
Phase-1 Amdahl ceiling.

## Within-Phase-1 breakdown (perf)

Flamegraphs in `benchmark/results/latest/perf/` (chr21, celegans, `-t1`/`-t8`). Note: this is a hybrid
P/E CPU with split PMUs (`cpu_core`/`cpu_atom`), so cycle-based `-t1` captures under-sample; the clean
`-t1` attribution uses **`task-clock`** (time-based, topology-agnostic) on celegans.

- **perf agrees with the timers.** celegans `-t1`, task-clock: **Phase 1 = 28.6 %, Phase 2 = 70.6 %**
  of samples — matching the `SKLIB_TIMING` split (31 % / 69 %) within noise. The whole measurement
  chain is consistent.
- **The producer is compute-bound on minimizer math, not I/O or parsing.** Inside Phase 1: **FASTA
  parse ≈ 0 %**, **bucket routing/write = 0.5 %**. Effectively all of the producer's time is
  super-k-mer / minimizer computation — `SkmerManipulator` (`canonicalize`, `reverse_complement`,
  rolling-minimizer `new_minimizer`, `kmer_compare`, `mask_absent_nucleotides`) driven by
  `SeqSkmerator::Iterator::operator++`. So the serial floor is CPU work, and the page-cached FASTA
  read and the buffered bucket writes cost almost nothing.
- **Phase-2 hotspot = the column algorithm** (`generate_sorted_list_from_enumeration` → `sort_column`
  + `colinear_chaining` + `get_candidate_overlaps`), ~59 % of all cycles at `-t8` — the documented
  ~65–70 % build hotspot, now spread across the worker pool.
- **No synchronization bottleneck.** Lock/futex/condition-variable time is **0.0 %** at every `-t`.
  The single-writer mutex + bounded reorder buffer in `parallel_build_phase2` is *not* the limiter;
  the Phase-2 efficiency loss at high `-t` is the memory-bandwidth + hybrid-core effect (the
  ~1.9–2.2× CPU-work inflation), not contention.

## Diagnosis

The hypothesis is confirmed: **the single-threaded producer (Phase 1) is the scaling bottleneck.**
It is the loop in `build_bucketed` that streams the FASTA through `FileSkmerator` (parse + rolling
minimizer + super-k-mer assembly), computes the minimizer, and routes each record to a per-minimizer
on-disk bucket via `SkmerBucketWriter::append`. `-t` never touches it. A secondary ceiling
(Phase-2 memory-bandwidth/hybrid-core inefficiency) limits the high-`-t` tail.

## Recommendation

**Parallelize the producer (Phase 1).** It is the ~⅓ serial fraction that Amdahl-caps the build, and
perf shows it is pure CPU work (minimizer/super-k-mer computation) with ~0 % parse and ~0.5 % I/O —
i.e. **embarrassingly parallel**: each FASTA region's super-k-mers are independent.

Cleanest correctness-preserving design — **sharded multi-producer Phase 1**: split the FASTA across
worker threads (by sequence, or by byte offset re-synced to the next header), each thread computing
minimizers and routing into its **own** per-bucket writer (`bucket_<id>_<tid>.bin`); Phase 2 then
loads all shards of a bucket. The per-bucket `sort_and_dedup` (B(1)) already makes the output
**independent of routing order**, so the index stays byte-identical; only the B(2) consecutive-dup
*optimization* needs seam handling (it is an optimization, not a correctness requirement). Reuse the
existing `parallel_for_dynamic` (`Parallel.hpp`) for the chunk loop.

**Expected lift (order of magnitude).** Model `wall = phase1/eff₁(t) + phase2(t)`. If the producer
parallelizes at the efficiency Phase 2 already reaches at `-t8` (~0.74 — likely conservative, since
the producer is more compute- and less memory-bound than sorting), celegans `-t8` would go from
`7.66 + 2.78 = 10.43 s` to `≈1.27 + 2.78 = 4.05 s` — a **~2.6× build, i.e. ~5.9× vs `-t1`** (today:
2.31×). The knee at `-t8` roughly doubles. Beyond that, the **memory-bandwidth wall** (the Phase-2
CPU-work inflation) becomes the new ceiling, so the biggest, safest wins are at moderate `-t` (4–12),
not at `-t22`. A lighter complementary lever: a cheaper rolling minimizer (ntHash-style) shrinks the
serial floor directly.

The single-writer Phase 2 needs **no** change — perf found zero lock contention; its high-`-t`
inefficiency is hardware bandwidth, addressable only by doing less memory traffic (narrower records /
cache-friendlier compaction), not by more threads.

## Reproduce

```bash
# build (clang-18) with the SKLIB_TIMING phase split:
CC=clang-18 CXX=clang++-18 cmake -S . -B build-timing -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=OFF
cmake --build build-timing -j --target sskm

# fetch genomes + run the sweep + analyse:
bash benchmark/scripts/construct_scaling.sh           # -> benchmark/results/latest/construct_scaling.csv
python3 benchmark/scripts/diag_plot.py                # tables + figs/construct_{speedup,phase_split}.png

# phase split for a single build:
SKLIB_TIMING=1 build-timing/bin/sskm construct -k21 -m11 -f genome.fa -o /tmp/x.sskm -t8
#   -> [sklib-timing] ... phase1_s=… phase2_s=…

# perf flamegraphs (needs kernel.perf_event_paranoid <= 1):
bash benchmark/scripts/diag_perf.sh                   # -> benchmark/results/latest/perf/*.svg
```
