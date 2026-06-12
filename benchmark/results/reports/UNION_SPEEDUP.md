# `set_union` — monothread speedup (iteration journal)

**Goal.** Optimise, as far as possible and **single-threaded**, the union of two sorted super-k-mer
lists (`km::sortedlist::set_union` and its sub-functions) without ever regressing performance, and
while preserving the output's **content-equivalence**: a valid, maximally-compacted virtual
super-k-mer list (`VSKMER_4`, input bucket layout) whose k-mer set is exactly A∪B, queryable
identically. Byte-level packing may change; record count may not grow.

## Setup (held constant for the whole journal)

- **Machine.** Intel Core Ultra 7 165H (22 threads), Linux 6.17, pinned to one core (`taskset -c 5`).
- **Build profile.** `clang++-18` (18.1.3), `-DCMAKE_BUILD_TYPE=Release` (`-O3 -DNDEBUG`),
  `-DPORTABLE_BUILD=OFF` (`-march=native`), no LTO (the project's `ipo_supported` guard is unset, so
  `sskm` has none either — representative). Build dir `build-union`, `-DWITH_UNION_BENCH=ON`.
- **Harness.** `benchmark/union_bench/union_bench.cpp` calls the real `set_union<store>(A,B,out,
  /*no_compact*/false, /*nthreads*/1)` on two pre-built lists; output to `/dev/shm` (no write-I/O
  noise); inputs page-cache-warm. Driver: `benchmark/scripts/union_bench.sh`.
- **Data.** Dev on **chr21**, non-regression verdict on **celegans**. k=31 (m=15 → `uint64` store,
  8 B) and k=63 (m=31 → `__uint128` store, 16 B). A/B pairs at target **Jaccard ∈ {0.1, 0.5, 0.9}**
  (B = `mutate.py` copy of A). → **12 configs** (2 datasets × 2 k × 3 J); all must hold before commit.
- **Protocol.** warmup 2, reps 9; report **median / min / stddev / MAD** (seconds).
- **Decision threshold.** A gain is *confirmed* only if the **median improves by > 2× the observed
  stddev** on the dev configs **and** no config (the full 12) regresses (median_opt ≤ median_base +
  2σ). Within noise ⇒ **neutral**, not committed. Each idea is A/B'd back-to-back vs the current best
  committed state in the same machine state.
- **Correctness gates.** L0: `ctest` set-op suite green. L1 (every idea, 12 configs): harness
  `--verify` — union k-mer count == `union_size(A,B)`, `intersection_size(O_opt,O_ref) ==` that count
  (⇒ identical sets vs the frozen reference), records(O_opt) ≤ records(O_ref). L2 (pre-commit):
  `tests/setop_verif.sh` (KMC cross-check, independent of sklib's merge).

## Extraction validation

The harness output and `sskm setop -a A -b B --op union -o out -t 1` are **bit-identical** (sha256
`b9ba386e…` on chr21 k=31 J=0.5), proving the harness exercises the production path exactly.

## Baseline (commit `135283b`, monothread, `union_bench.sh all`)

| dataset | k | J | store | median s | stddev s | Mkmer/s | records |
|---|--:|--:|--:|--:|--:|--:|--:|
| chr21 | 31 | 0.1 | 8 | 6.761 | 0.057 | 9.89 | 8 311 037 |
| chr21 | 31 | 0.5 | 8 | 4.723 | 0.034 | 10.07 | 5 918 694 |
| chr21 | 31 | 0.9 | 8 | 3.403 | 0.009 | 10.78 | 4 410 010 |
| chr21 | 63 | 0.1 | 16 | 11.508 | 0.233 | 6.01 | 4 351 894 |
| chr21 | 63 | 0.5 | 16 | 7.396 | 0.174 | 6.74 | 3 145 316 |
| chr21 | 63 | 0.9 | 16 | 5.192 | 0.007 | 7.45 | 2 353 678 |
| celegans | 31 | 0.1 | 8 | 19.334 | 0.081 | 9.09 | 22 550 865 |
| celegans | 31 | 0.5 | 8 | 13.401 | 0.048 | 9.50 | 16 291 463 |
| celegans | 31 | 0.9 | 8 | 9.843 | 0.053 | 10.09 | 12 299 427 |
| celegans | 63 | 0.1 | 16 | 32.462 | 0.343 | 5.50 | 11 097 434 |
| celegans | 63 | 0.5 | 16 | 23.496 | 0.248 | 5.53 | 8 106 029 |
| celegans | 63 | 0.9 | 16 | 17.776 | 0.187 | 5.73 | 6 120 171 |

Throughput is a steady ~10 Mkmer/s at k=31 (`uint64`) but only ~6 at k=63 (`__uint128`): the wider
records make the re-compaction's `pair` arithmetic markedly more expensive per k-mer. Cross-run
absolute medians carry ~2-3 % drift (worse at k=63); **decisions use interleaved A/B**
(`union_ab.sh`), whose per-round interleaving cancels drift (bands: ~±0.3 % at k=31, ~±1.5 % at k=63).

## Profile (commit `135283b`)

Env-gated phase timer (`SKLIB_UNION_PHASE_TIMING=1`): the per-bucket **re-compaction**
(`generate_sorted_list_from_enumeration`) is **87-91 %** of the union; merge+collect 8-11 %; read+write
< 1.5 %. So all effort targets the re-compaction. `perf` self-time inside it (k31 / k63):

- `get_candidate_overlaps` (sort of `right_keys` + per-left `lower_bound`): ~48 % / ~42 %
- `sort_column` (re-scans the WHOLE enumeration via `has_valid_kmer` once **per column**): 11 % / **29 %**
  — `has_valid_kmer` alone is **22 % self** at k=63 (32 columns × full scan)
- `merge_LList_column` ~10 % / 7 %; `greedy_chaining` ~16 % / 10 %
- `Skmer::pair::operator<` (the 2-word compare, shared by every sort/search): 14 % / 19 % self

## Ideas

Profile-first; the obvious per-column sort is already short-circuited (`sort_column` `is_sorted`
fast-path) and `get_candidate_overlaps` already uses a sorted-array merge, so remaining wins are
subtler. Each entry: mechanism, numbers (interleaved A/B unless noted), status, commit.

### #1 — column-offset fast-path for `sort_column` — **COMMITTED** (`<hash>`)

**Mechanism.** The set-op re-compaction feeds `generate_sorted_list_from_enumeration` an enumeration
that is already **column-grouped and per-column sorted+distinct** (that is exactly what
`merge_columns` emits). The generic recompaction ignores this and re-discovers each column by scanning
the whole enumeration with `has_valid_kmer` once per column — `O(n·(k-m))`, the dominant cost,
disproportionately so at k=63 (twice the columns). The `CollectSink` now records the per-column kept
count; `materialize_setop` prefix-sums it into block offsets and passes them to the recompaction, which
takes column `c`'s ids as the contiguous range `[off[c], off[c+1])` — **`O(n)` total**, no
`has_valid_kmer` scan, no `is_sorted`, no `unique`. Produces the **same ids** the scan would, so the
output is **byte-identical** (`Skmer.hpp` CollectSink + `materialize_setop` in `SetOperations.hpp`;
`generate_sorted_list_from_enumeration`/`column_ids` in `VirtualSkmer.hpp`).

**Correctness.** Output **byte-identical** to baseline (sha256 match, chr21 k31 & k63 J0.5); harness
`--verify` PASS; 213/213 unit tests; KMC cross-check PASS at k=31 and k=63.

**Result (interleaved A/B, `/tmp/ub.timer` baseline → idea#1):**

| config | baseline s | idea#1 s | delta | band |
|---|--:|--:|--:|--:|
| chr21 k31 J0.5 | 4.762 | 4.190 | **−12.0 %** | ±2.25 % |
| chr21 k63 J0.5 | 7.352 | 5.323 | **−27.6 %** | ±0.37 % |

**Full 12-config sweep** (absolute medians, idea#1 vs baseline; ~2-3 % cross-run drift, but every gain
is far above it; all 12 `verify=PASS`, byte-identical):

| dataset | k | J | baseline s | idea#1 s | gain |
|---|--:|--:|--:|--:|--:|
| chr21 | 31 | 0.1 | 6.761 | 5.849 | −13.5 % |
| chr21 | 31 | 0.5 | 4.723 | 4.149 | −12.2 % |
| chr21 | 31 | 0.9 | 3.403 | 3.025 | −11.1 % |
| chr21 | 63 | 0.1 | 11.508 | 7.510 | −34.7 % |
| chr21 | 63 | 0.5 | 7.396 | 5.317 | −28.1 % |
| chr21 | 63 | 0.9 | 5.192 | 3.877 | −25.3 % |
| celegans | 31 | 0.1 | 19.334 | 16.354 | −15.4 % |
| celegans | 31 | 0.5 | 13.400 | 11.721 | −12.5 % |
| celegans | 31 | 0.9 | 9.843 | 8.695 | −11.7 % |
| celegans | 63 | 0.1 | 32.462 | 20.728 | −36.1 % |
| celegans | 63 | 0.5 | 23.496 | 14.854 | −36.8 % |
| celegans | 63 | 0.9 | 17.776 | 11.077 | −37.7 % |

k31 (`uint64`) −11 to −15 %; k63 (`__uint128`) **−25 to −38 %** — the fast-path removes the `O(n·(k-m))`
scan that hit k=63 (32 columns) hardest. Throughput at k=63 rose ~6 → ~9 Mkmer/s.

| # | idea (file) | mechanism | result (k31 / k63, J0.5) | status |
|--:|---|---|---|---|
| 1 | column-offset fast-path for `sort_column` (`SetOperations.hpp`, `VirtualSkmer.hpp`) | skip the per-column `has_valid_kmer` full-scan; ids = contiguous block from merge's per-column counts | **−12.0 % / −27.6 %** (byte-identical) | committed `<hash>` |
