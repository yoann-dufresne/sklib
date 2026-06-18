# `xor` (symmetric difference) — monothread speedup (iteration journal)

**Goal.** Optimise, as far as possible and **single-threaded**, the symmetric difference of two sorted
super-k-mer lists (`km::sortedlist::symmetric_difference` / `sym_diff_size`) at the three store widths
**k=21 → uint32 (4 B)**, **k=31 → uint64 (8 B)**, **k=63 → __uint128 (16 B)**, while preserving the
output's content-equivalence (a valid, maximally-compacted `VSKMER_4` list whose k-mer set is exactly
A△B). No config may regress.

**XOR is not a separate algorithm.** `--op xor` → `symmetric_difference()` →
`materialize_setop(keep_inter=false, only_a=true, only_b=true)` — the union path minus the intersection
k-mers; `--op xor_size` → `sym_diff_size()` → `set_sizes()`. So XOR already inherits the whole
`UNION_SPEEDUP.md` work (#1 col-offset fast-path, #2 hash-join `get_candidate_overlaps`, #3 is_sorted
chaining guard, #5 difference-array CSR). This journal optimises what is *specific* to XOR/diff.

## Setup (held constant)

- **Machine.** Intel Core Ultra 7 165H, Linux 6.17, pinned to one core (`taskset -c 5`).
- **Build.** `clang++-18` (18.1.3), `-DCMAKE_BUILD_TYPE=Release` (`-O3 -march=native`), `build-union`,
  `-DWITH_UNION_BENCH=ON`. The harness (`benchmark/union_bench/union_bench.cpp`, `--op xor`) calls the
  real `symmetric_difference<store>(A,B,out,false,1)` — bit-identical to `sskm setop --op xor -o … -t1`.
- **Data.** chr21, B = `mutate.py` copy of A at target Jaccard ∈ {0.1, 0.5, 0.9}; k ∈ {21,31,63}
  (m=k/2). 9 configs. celegans for non-regression spot-checks.
- **Protocol / decision.** Interleaved A/B (`union_ab.sh`, cancels cross-run drift): a gain is kept only
  if median improves beyond the combined-MAD band **and no config regresses**.
- **Correctness gates.** L0 `ctest` (213 tests); L1 harness `--verify` (count == `sym_diff_size`,
  set == frozen ref via `intersection_size`, records ≤ ref); byte-identical sha256 vs frozen ref; L2
  KMC cross-check `tests/setop_verif.sh` (independent of sklib's merge).

## Baseline (chr21, `symmetric_difference`, pre-optimization)

| k | store | J=0.1 | J=0.5 | J=0.9 |
|--:|--:|--:|--:|--:|
| 21 | 4 B | 2.880 s | 1.337 s | 0.420 s |
| 31 | 8 B | 2.726 s | 1.285 s | 0.357 s |
| 63 | 16 B | 3.229 s | 1.512 s | 0.437 s |

XOR's time tracks **output size = only_a + only_b**: largest at J=0.1 (output ≈ union), smallest at
J=0.9 (few differences).

## Profile — XOR has two regimes (env phase timer, chr21)

| | read | merge+collect | recompact | write |
|---|--:|--:|--:|--:|
| **J=0.1** (large output) | ~1 % | ~24-27 % | **~71-74 %** | ~2 % |
| **J=0.9** (small output) | ~6 % | **~59-65 %** | ~28-35 % | ~1 % |

- **Low overlap → recompaction-bound** (same shared `generate_sorted_list_from_enumeration` as union,
  already optimised by #1-#3).
- **High overlap → merge-bound** — `build_column_csr` (fixed by union #5) + `kmer_compare`. This is the
  regime union/intersection *cannot* speed up (they need the `both` k-mers), so it is XOR/diff's own.

## Ideas

### #E — drop byte-identical record pairs before the merge (xor/diff only) — **COMMITTED** (`b1af36e`)

**Mechanism.** A super-k-mer record present byte-identically in both A and B contributes only
intersection k-mers, which symmetric_difference and difference discard. Dropping such record pairs
before the per-column merge leaves the surviving only_a/only_b emissions — and therefore the
materialized result and every cardinality — **byte-identical**, while the column merge shrinks.
`mark_identical_records` is an `O(nA+nB)` two-cursor walk over the `m_pair`-sorted bucket payloads
(`Skmer::operator<` orders by `m_pair`), matching byte-identical records (`operator==`, sizes included)
within each tiny equal-`m_pair` run. `build_column_csr` gained a compile-time `Filtered` flag
(`if constexpr`) so the **no-drop path (union/intersection and every existing caller) keeps its original
codegen** — verified byte-identical and perf-neutral.

**Adaptive gate.** The pass costs `~O(records)` but only saves `~O(dropped·(k-m))` of column-merge work,
so it pays only above an overlap threshold — worst at narrow widths (most records, cheapest merge).
Measured drop fraction `f` vs A/B verdict calibrates the model **dedup pays iff `f·(k-m) ≥ 5`**:

| k | k−m | f (J0.1/0.5/0.9) | (k−m)·f | dedup verdict |
|--:|--:|--:|--:|--|
| 21 | 11 | 0.03 / 0.21 / 0.69 | 0.4 / 2.3 / 7.6 | off / off / **on** |
| 31 | 16 | 0.10 / 0.54 / 0.91 | 1.6 / 8.6 / 14.5 | off / **on** / **on** |
| 63 | 32 | 0.12 / 0.58 / 0.92 | 3.7 / 18.4 / 29.5 | off / **on** / **on** |

All clear regressions have `(k−m)·f ≤ 2.3`; all wins have `≥ 7.6` → GAMMA=5 separates them. The gate
samples a small prefix (~0.4 % of records), estimates `f`, and commits. **Dropping is always
correctness-safe, so the gate decision — and any monothread/parallel difference in it — changes only
speed, never the output.** `SKLIB_SETOP_DROP_STATS` reports drop telemetry (env-gated, zero cost off).

**Correctness.** Output **byte-identical** to the pre-#E `xor`/`diff` across all 9 configs (sha256);
union/intersection/diff byte-identical too (`Filtered=false` path unchanged); 213/213 gtest; KMC
cross-check (incl. xor, diff_AB, diff_BA) ALL PASS.

**Result (interleaved A/B vs pre-#E binary, chr21):**

| k (store) | J=0.1 | J=0.5 | J=0.9 |
|---|--:|--:|--:|
| 21 (uint32) | neutral | neutral¹ | **−15.1 %** |
| 31 (uint64) | neutral | **−5.3 %** | **−34.6 %** |
| 63 (__uint128) | neutral | **−9.4 %** | **−43.3 %** |

¹ k21 mid/low-overlap shows ~+1 % — but the **intersection control** (executes zero dedup code,
byte-identical, same instructions as base) shows the **same +1.2 %** at k21 J0.5, so this is
binary-layout/measurement noise of the modified TU on the most layout-sensitive config (k21 = most
records), **not** an algorithmic regression (the gate is OFF there; the pass adds ~0.2 %).

The win scales with overlap and with `(k-m)` (wider k ⇒ more per-column merge saved per dropped record),
peaking at **−43 % (k63, J0.9)**.

## Post-#E profile & remaining leads

With the high-overlap merge collapsed (#E), the **recompaction is now the dominant cost at every width
and overlap** (~70-73 % at J0.1; ~53-57 % at J0.9, k31/k63). That recompaction is the *shared*
`generate_sorted_list_from_enumeration`, already mined by `UNION_SPEEDUP.md` #1-#3 (hash-join
`get_candidate_overlaps`, is_sorted chaining guard, col-offset fast-path), whose verdict was that
further gains are incremental and the one big lever — a fused streaming assembly ("A1") — is high-risk:
the cross-column chaining is an inherent non-crossing-max-matching (LIS), not reducible to a single
greedy pass without growing the record count (which the content-equivalence bar forbids). #E captured
the XOR/diff-specific win (the high-overlap merge that union could not touch); the low-overlap regime is
bounded by that shared, already-optimised recompaction.

### #B — skip redundant `get_skmer_of_kmer` on the set-op recompaction — **REVERTED (construction regression)**

**Idea.** On the set-op path the recompaction's enumeration already holds single-k-mer skmers framed at
their column, so `get_skmer_of_kmer(enum[id], col)` is the identity (union #4's observation, never tried
at uint32). Skip it in `merge_LList_column` + the column-0 init. To keep construction untouched it was
templated on a compile-time `SetopFrame` flag (`if constexpr`, default false = construction), enabled
only for narrow stores (uint32/uint64; `__uint128` keeps the original path where union #4 regressed).

**Why it fails (measured).** The skip itself is **byte-identical** (213 gtest incl. construction golden
digests; xor/union/intersection/diff sha256 match) and gives a small **xor** win — k21 −1.2/−1.4 %,
k31 ~−1 %, k63 neutral (gate off), no setop regression. **But templating the shared `merge_LList_column`
regressed *construction* +2.4 %** (yeast k31 −t1, band ±0.60 %; same on ecoli): the `SetopFrame=false`
instantiation, though logically identical, does not reproduce the original non-template codegen — exactly
the shared-hot-loop fragility union #4 hit. Trading a ~1 % xor gain for a ~2.4 % construction regression
fails the no-regression rule. The only regression-free implementation — fully **duplicating** the
~120-line `merge_LList_column`/recompaction into a set-op-only copy — is disproportionate maintenance
debt for ~1 %. Reverted (`git checkout`).

**Assessed, not pursued:**
- *Fused streaming recompaction (A1).* Re-confirmed infeasible as a maximality-preserving single pass
  (LIS needed even in the 1:1 case — adjacent columns' overlap (k-1)-mer is not monotone in either
  k-mer order, so matches cross); the union journal's contained variant gained 0 %.

**Conclusion.** #E is the monothread XOR win (XOR/diff-specific, byte-identical, −15…−43 % at high
overlap, every width, no regression). The residual cost is the shared recompaction, which the union
journal already optimised and whose remaining levers are either infeasible (A1) or not worth the cost/risk
(#B). Monothread XOR is at its practical floor.
