# Construction speedup — parallel per-bucket compaction (phase 2)

**TL;DR.** `sskm construct` is no longer single-threaded. The per-bucket compaction (phase 2 —
`generate_sorted_list_from_enumeration`, ~65-70 % of build time and the construction hotspot) now
runs on a worker pool, writing the result back in bucket-id order. The index is **byte-identical**
to the old sequential build for any thread count (sha256-verified) and peak RAM stays the same order
of magnitude (default `-t 8`: 1.0-1.9× the sequential RSS). End-to-end **~2-2.5×** on real genomes,
Amdahl-bounded by the still-sequential phase 1. On a fair same-machine multi-thread comparison,
sklib builds **3-8× faster than sbwt and sshash at 10-22× less RAM** — and even *single-threaded*
sklib beats both competitors running on 22 threads.

## What changed

- `lib/include/algorithms/ParallelConstruct.hpp` — `parallel_build_phase2()`: N workers dynamically
  claim buckets (one atomic counter ⇒ the ~0.2-0.3 %-of-work bucket skew self-balances), each owning
  a reusable `BucketCompactor` (`load_bucket` → `sort_and_dedup` →
  `generate_sorted_list_from_enumeration` → `truncate_skmer`); a bounded reorder buffer (`2·threads`)
  feeds a single writer thread that emits each bucket's payload **in increasing bucket-id order**.
- `SortedSkmerListBuilder.hpp` — `build_bucketed` phase 2 dispatches to the parallel driver for
  `n_threads ≥ 2`, keeps the (buffer-reusing) sequential loop for `-t 1`.
- `app/sortedlists/commands.cpp` — `mallopt(M_ARENA_MAX, 4)` on parallel builds (glibc): caps the
  per-thread malloc arenas that otherwise inflate peak RSS at high `-t` for no speed gain.
- CLI `construct -t/--threads` (default 8); harness knob `CONSTRUCT_THREADS`.
- Test `tests/km/ParallelConstruct.cpp` locks the byte-identical invariant.

Only *when* a bucket is compacted differs — never *what* it contains nor *where* it lands in the
file. The on-disk layout, compaction algorithm, bucketing and integer widths are unchanged.

## Correctness — identical index (the hard constraint)

`sha256sum` of the produced `.sskm`, on chr21 k=21 m=11, all equal to **one** hash:

```
baseline (dev)      da542f7d93b1c396…
opt -t 1            da542f7d93b1c396…
opt -t 8            da542f7d93b1c396…
opt -t 16 / -t 22   da542f7d93b1c396…
```

Also verified identical at k=41 m=21 (`__uint128_t` store) and across `uint64`/quotiented bucketings
in the gtest. ⇒ index byte-size, `bits/k-mer` and content are unchanged at every thread count.

## sklib — wall time vs threads (k=21, m=11)

Speedup vs the sequential build; RSS is the tightened binary's peak (arena-capped + `2·threads`
reorder buffer). Machine: 22 logical cores; ~10 % run-to-run wall variance.

| dataset | base (1 thr) | -t2 | -t4 | **-t8 (défaut)** | -t16 | -t22 | RSS base → -t8 → -t22 |
|---|--:|--:|--:|--:|--:|--:|--:|
| ecoli (4.7 Mbp) | 1.2 s | 1.3× | 1.8× | 1.5× | 2.1× | 1.7× | 18 → 18 → 22 MB |
| yeast (12 Mbp) | 3.0 s | 1.3× | 1.6× | **2.0×** | 2.05× | 2.05× | 18 → 18 → 34 MB |
| chr21 (40 Mbp) | 10.2 s | 1.35× | 1.7× | **1.97×** | 2.05× | 2.11× | 44 → 44 → 63 MB |
| celegans (102 Mbp) | 30.1 s | 1.6× | 2.1× | **2.45×** | 2.6× | 2.5× | 49 → 92 → 190 MB |
| chr1 (234 Mbp) | 55.7 s | 1.5× | 1.9× | **2.24×** | 2.42× | 2.44× | 138 → 183 → 410 MB |

- `-t 1` ≈ baseline (0.96-1.06×): no regression of the sequential path (on big genomes it is even
  slightly faster — celegans 28.4 vs 30.1 s, chr1 54.2 vs 55.7 s; the small-dataset jitter is noise).
- **`-t 8` is the sweet spot**: past it, time barely moves (Amdahl: phase 1 — FASTA parse + skmer
  generation + bucket routing — is sequential, ~40-50 % here) while RAM keeps climbing. `-t 16/-t 22`
  buy ≤ 8 % more speed for 2-4× the RAM.

### RAM tightening (peak RSS, before → after the two levers)

| dataset | base | -t8 | -t16 | -t22 |
|---|--:|--:|--:|--:|
| chr21 | 44 MB | 49 → **44** | 75 → **51** | 100 → **63** |
| celegans | 49 MB | 120 → **92** | 202 → **151** | 267 → **190** |
| chr1 | 138 MB | 227 → **183** | 403 → **319** | 519 → **410** |

`mallopt(M_ARENA_MAX,4)` was the bigger lever (glibc creates 8×cores arenas, each retaining its
high-water mark) and costs **zero** time (compaction is compute-bound, not malloc-bound); tightening
the reorder buffer `4·threads → 2·threads` trims the rest. Default `-t 8` RSS is now ~1.0-1.9× base.

## Fair multi-thread comparison vs competitors (honest cold FASTA → index)

Same machine, same inputs, same thread counts. **sklib** = `sskm construct -tN`; **sshash** =
`bcalm -nb-cores N` (FASTA→unitigs, mandatory) + `sshash build -t N`; **sbwt** = `sbwt build -t N`.
Time in s / peak RSS in MB.

| dataset | N | **sklib** | sshash (bcalm+build) | sbwt |
|---|--:|--:|--:|--:|
| chr21 | 1 | **8.5 / 44** | 71.6 / 2103 | 15.8 / 427 |
| chr21 | 8 | **4.4 / 44** | 24.3 / 590 | 12.2 / 1125 |
| chr21 | 22 | **3.8 / 63** | 24.3 / 658 | 11.5 / 1483 |
| celegans | 8 | **9.6 / 89** | 81.0 / 1066 | 31.4 / 1978 |
| celegans | 22 | **9.2 / 197** | 76.4 / 1283 | 30.1 / 2683 |
| chr1 | 8 | **24.7 / 179** | 173.4 / 2703 | 63.9 / 2456 |
| chr1 | 22 | **21.2 / 415** | 159.2 / 2703 | 62.2 / 3353 |

- At `-t 8` on chr1, sklib is **7.0× faster than sshash and 2.6× faster than sbwt**, at **15× / 14×
  less peak RAM**. On celegans: **8.4× / 3.3× faster**, **12× / 22× less RAM**.
- **Even single-threaded sklib beats both competitors at 22 threads** (chr21: 8.5 s vs sshash 24.3 s,
  sbwt 11.5 s). The gap *widens* with threads: sshash is bcalm-bound (plateaus ~24 s; bcalm = 67 s/2.1
  GB at 1 core, 21 s/562 MB at 8) and sbwt barely scales while its RAM balloons to 2.7-3.4 GB.
- **Methodology / honesty.** sshash needs unitig input, so the honest FASTA→index cost includes
  BCALM — as the repo's own wrapper intends. The previously-recorded "sshash ≈ 5.5 s" had cached or
  excluded that step; counting it (and using `bcalm -abundance-min 1` to match sklib's all-k-mers
  semantics) is what makes this apples-to-apples. If you *already* hold unitigs, sshash's build-only
  step is ~1-8 s — but that is not a FASTA→index cost. sklib's index is larger in bits/k-mer (a
  pre-existing structural trade-off, unchanged here); construction **time and RAM** are where it wins.

## Next lever

Phase 1 (FASTA parse → skmer generation → bucket routing, single producer) is now the floor.
Parallelizing it — a multi-producer parser, or offloading super-k-mer generation off the reader
thread — would lift the ~2-2.5× Amdahl ceiling. Deferred here: more invasive (ordering/determinism
of the bucket stream) and out of this iteration's conservative, byte-identical scope.

## Reproduce

```bash
# byte-identical check (any thread count -> one hash):
for t in 1 8 22; do build/bin/sskm construct -k21 -m11 -f genome.fa -o idx_t$t.sskm -t $t; done
sha256sum idx_t*.sskm

# sklib construction sweep via the harness (CONSTRUCT_THREADS feeds construct -t):
CONSTRUCT_THREADS=8 DATASETS="chr21 celegans chr1" KM="21,11" bash scripts/bench/bench.sh

# competitors, built once, then bcalm -nb-cores N + sshash build -t N / sbwt build -t N:
bash scripts/bench/tools/setup.sh bcalm sshash sbwt
```
