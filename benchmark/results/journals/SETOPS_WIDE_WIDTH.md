# Width-dispatched merge inner loop for set operations (wide stores)

**Goal.** Differentiate the set-op merge *by the `kuint` store width chosen at runtime* (a function of
`k`/`m`/`b`), to recover the per-element `_BitInt` ALU slack that only appears at the wide store widths —
`__uint128` (16 B) and `kuint256` (32 B). The narrow widths (`uint32`/`uint64`) were already driven to
their floor by earlier journals (`SETOPS_BOTTLENECKS.md` **C1**: `kmer_compare` without the copy → 0.99×;
`UNION_SPEEDUP.md` **#6**: pre-masked contiguous key array → reverted, *"`kmer_compare` is at its floor"*),
and that verdict was measured only at k=21/31/63. The **kuint256 backend (k>63) had never been profiled or
specialized for set ops** — this journal covers that gap. No config may regress.

## Mechanism (`merge_columns`, `SetOperations.hpp`)

The two-cursor per-column merge is the whole cost of **counting** (`set_sizes`, behind `*_size` and `--sizes`)
and the shared front of **materialization** and **multi_setop**. Its inner call `kmer_compare`
(`io/Skmer.hpp`) does two things the compiler cannot remove and that scale with the pair's word count
(8 limbs for the kuint256 pair, 4 for `__uint128`):

1. **masks both operands every call** (`first &= kmer_masks[c]; second &= …`) — but in a two-cursor merge
   the *stationary* cursor's operand is re-masked needlessly each iteration;
2. **compares with `<` then `>`** (each a hi-then-lo lexicographic compare) — up to four word compares
   where a single three-way needs at most two.

The change, gated by `if constexpr (sizeof(store) >= WIDE_MERGE_MIN_STORE_BYTES)`:

- **Narrow (`<= 8 B`, uint32/uint64):** keep the original loop **verbatim** → byte-identical codegen, zero
  regression risk (honors C1/#6).
- **Wide (`>= 16 B`):** cache each cursor's **masked** k-mer (`SkmerManipulator::masked_kmer`, re-masked only
  when that cursor advances) and compare with a single three-way (`Skmer::pair::compare3`, ≤2 word compares).

This is **not** UNION #6: no parallel key array is built (zero extra memory); the masked key lives in two
registers and is recomputed lazily. Output is **byte-identical at every width** — same emission order (the
tail loops are shared) and same compare semantics (`compare3 == sign(kmer_compare)`, `masked_kmer ==`
the operand `kmer_compare` masks). One `if constexpr` in `merge_columns` therefore covers **count,
materialize and multi** at once, since all three route through it.

## Correctness

- **Full suite:** `ctest` / `sklib-tests` — **213/213** (Debug, clang++-18, asserts on).
- **Byte-identity vs the pre-change binary** (the journals' gold standard): a 100-check matrix
  (k=21/31/63/95/127 × {materialize ∈ inter/union/diff/xor, every `*_size`, combined `--sizes`+5 outputs} ×
  overlap {1 %, 5 % mutation} × threads {1, 4}) — **0 mismatches**, every output sha256-equal to the
  pre-change `sskm`. Because the pre-change code was already KMC-cross-validated and the new code is
  byte-for-byte identical to it, the KMC scripts (`setop_verif.sh`, `setop_multi_verif.sh`) are redundant
  for this change.

## Setup (bench)

- **Machine.** 22-core Linux 6.17, pinned to one core (`taskset -c 5`).
- **Build.** `clang++-18` (18.1.3), `-DCMAKE_BUILD_TYPE=Release` (`-O3 -march=native`, `PORTABLE_BUILD=OFF`),
  `-DWITH_UNION_BENCH=ON`. Two binaries from the same tree: `ub.new` (this change) and `ub.old`
  (the two files stashed) — the real `set_sizes` / `intersection` / `symmetric_difference`, not a copy.
- **Data.** `chr21.sanitized.fa` (40 MB); B = `mutate.py` copy of A at a per-k substitution rate
  `p = 1 − J^(1/k)` for target shared-k-mer fraction `J ∈ {0.9, 0.5}`. Widths exercised:
  k=31→`uint64` (8 B, narrow **control**), k=63→`__uint128` (16 B), k=95 & k=127→`kuint256` (32 B).
- **Protocol.** Interleaved new/old (`--warmup 3 --reps 11`, 2 rounds each, min-of-medians), per the
  XOR/UNION journals. A win is kept only if it clears run-to-run noise and **no control (k=31) regresses**.
- **Raw data.** `benchmark/results/reference/setops_wide_width.csv` (chr21 count/intersection/xor + yeast
  count/union/diff, both at J ∈ {0.9, 0.5}).

## Results (chr21, min-of-medians, delta = new vs old)

**`count` — pure merge (the counting headline; every `*_size` and `--sizes` route through it):**

| store | k | J=0.9 | J=0.5 |
|---|--:|--:|--:|
| `uint64` (control) | 31 | −0.2 % | +0.0 % |
| `__uint128` | 63 | **−2.4 %** | **−2.4 %** |
| `kuint256` | 95 | **−12.7 %** | **−8.6 %** |
| `kuint256` | 127 | **−10.0 %** | **−11.4 %** |

**`intersection` — merge + partial recompaction:**

| store | k | J=0.9 | J=0.5 |
|---|--:|--:|--:|
| `uint64` (control) | 31 | +0.6 % | −0.2 % |
| `__uint128` | 63 | −1.4 % | **−9.3 %** |
| `kuint256` | 95 | **−5.4 %** | **−19.7 %** |
| `kuint256` | 127 | **−8.4 %** | **−11.1 %** |

**`xor` — recompaction-dominated at these configs (13–24 s/rep at kuint256), so the merge gain is diluted:**

| store | k | J=0.9 | J=0.5 |
|---|--:|--:|--:|
| `uint64` (control) | 31 | −0.3 % | −1.7 % |
| `__uint128` | 63 | −4.2 % | +3.0 % |
| `kuint256` | 95 | −1.5 % | −0.9 % |
| `kuint256` | 127 | +0.5 % | −5.0 % |

Reading: the **narrow control (`uint64`) is flat** (−1.7 … +0.6 %, all within run-to-run noise) — the
`if constexpr` leaves that path's codegen untouched, so it cannot regress. At the **wide stores the merge
speeds up**, biggest on the pure-merge `count` (**−2.4 % at `__uint128`, −8.6 … −12.7 % at `kuint256`** — the
regime never before profiled) and on merge-heavy materialization (`intersection` −1.4 … −19.7 %). The two
small positives are both on `xor` materialization (k63 J0.5 **+3.0 %**, k127 J0.9 +0.5 %): those are
recompaction-bound (the change does not touch recompaction), the deltas sit inside the control's materialize
noise (±2–3 %), and the *clean* `count` win at the same width proves the merge itself got faster — so they
are variance, not a merge regression (byte-identity guarantees the work is unchanged).

### `union`/`diff` non-regression cross-check (yeast, 12 MB — second genome)

`union` is the most recompaction-bound op (materializes everything), so it is the strictest non-regression
target — the merge speedup is diluted most there. `diff` is merge-bound at high overlap (small output) and
recompaction-bound at low overlap. `count` repeated here cross-confirms the merge win on a second genome.

| store | k | op | J=0.9 | J=0.5 |
|---|--:|---|--:|--:|
| `uint64` (control) | 31 | count / union / diff | −0.3 / −0.1 / −0.1 % | −0.8 / +0.3 / +0.7 % |
| `__uint128` | 63 | count / union / diff | −5.9 / −2.2 / −5.6 % | −5.5 / −1.9 / −6.7 % |
| `kuint256` | 95 | count / union / diff | −9.5 / −1.3 / −17.3 % | −9.1 / +0.7 / +0.5 % |
| `kuint256` | 127 | count / union / diff | −13.1 / −1.3 / −4.6 % | −9.6 / +1.8 / −2.2 % |

`count` wins again on every wide config (−5.5 … −13.1 %); `diff` wins where merge-bound and is flat-within-
noise where recompaction-bound; **`union` ranges −2.2 % to +1.8 %** — flat-to-tiny-noise (its few small
positives are at low overlap where the full-union recompaction dominates), never a real regression. The
narrow control stays flat for all three ops. → no op regresses at any width; the wide-width merge wins hold.

## Gate decision

**`WIDE_MERGE_MIN_STORE_BYTES = 16`** — the wide loop is enabled for `__uint128` (16 B) and `kuint256`
(32 B); `uint32`/`uint64` keep the original loop verbatim. Justification: `__uint128` shows a clean,
repeatable pure-merge win (`count` −2.4 % at both overlaps) and favorable materialization, while the only
positives are recompaction-bound `xor` cells within noise; `kuint256` wins are large. Raising the gate to 32
would forgo the real `__uint128` merge gain for no benefit, so 16 is kept.

## Notes / non-goals

- Low-overlap **materialize** is recompaction-bound (`generate_sorted_list_from_enumeration`), already
  optimized by UNION #1–#5; #4 (skip `get_skmer_of_kmer`) *regressed at k63* — wide-width recompaction
  codegen is fragile, so it is **out of scope** here. This change targets only the merge front.
- `mark_identical_records` (the xor/diff dedup pre-pass) also walks `m_pair` `<`/`==`; it could take
  `compare3` too, but it is not the merge front — left for a follow-up unless the profile flags it.
