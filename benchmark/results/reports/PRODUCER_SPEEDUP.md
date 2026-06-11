# Super-k-mer producer — monothread speedup (iteration journal)

**Goal.** Shrink the serial producer floor (`CONSTRUCT_SCALING_DIAG.md`) with *single-threaded*,
**output-bit-identical** rewrites — the lever complementary to parallelizing Phase 1. Every change is
gated by the producer digest (`benchmark/results/reference/producer_digest.tsv`, chr21 + celegans) and
the full `ctest` suite; chr21 is the dev signal, celegans the verdict.

**Result.** Confirmed monothread wins, each **output-bit-identical** (digest unchanged on chr21 *and*
celegans, full `ctest` green) and measured in its own controlled back-to-back A/B (best-of-7, same
machine state, `sskm-produce`, each baseline = the previous committed state):

| # | change (file) | chr21 | celegans |
|--:|---|--:|--:|
| 1 | ring-buffer modulo `% (2k-m)` → power-of-two mask (`Skmerator.hpp`) | 2.31 → 2.64 Mskmer/s (**+14.3 %**) | 2.43 → 2.74 Mskmer/s (**+12.8 %**) |
| 2 | word-level `reverse_complement` — drop branchy `pair` shifts (`Skmer.hpp`) | 2.70 → 2.87 Mskmer/s (**+6.3 %**) | 2.81 → 2.97 Mskmer/s (**+5.7 %**) |
| 3 | precompute `mask_absent_nucleotides` flank masks — O(1) (`Skmer.hpp`) | 2.83 → 3.08 Mskmer/s (**+8.8 %**) | 2.95 → 3.18 Mskmer/s (**+7.8 %**) |
| 4 | word-parallel `reverse_complement` — nibble-swap + complement (`Skmer.hpp`) | 3.01 → 3.47 Mskmer/s (**+15.3 %**) | 3.15 → 3.56 Mskmer/s (**+13.0 %**) |
| 5 | cache minimizer φ-rank per buffer slot (`Skmerator.hpp`) | 3.37 → 3.60 Mskmer/s (**+6.8 %**) | 3.51 → 3.69 Mskmer/s (**+5.1 %**) |
| 6 | branch-free decode in `minimizer_is_ambiguous` (`Skmer.hpp`) | 3.44 → 3.53 Mskmer/s (**+2.6 %**) | 3.59 → 3.63 Mskmer/s (**+1.1 %**) |

**End-to-end, measured directly** (original producer @ commit 2441d4a vs current, same binary/harness,
best-of-11, same machine state):

| genome | original | current | speedup | wall |
|---|--:|--:|--:|--:|
| chr21    | 2.20 Mskmer/s | 3.51 Mskmer/s | **+59.5 %** | 3.063 s → 1.924 s (−37 %) |
| celegans | 2.28 Mskmer/s | 3.59 Mskmer/s | **+57.5 %** | 7.340 s → 4.667 s (−36 %) |

The digest is **identical between the original and the optimized producer** on both genomes
(chr21 `0xfd02872bc3f3ed34`, celegans `0x93d267f24daa8af8`) — an end-to-end proof that the whole
optimization series is output-bit-identical, so every `sskm construct` / `query` result built on the
producer is unchanged.

### Generalization across k / integer width

All the work above used k=21/m=11, which the build packs in **uint32** (width 4). The optimizations are
width-templated, so the same series was re-measured end-to-end (original @ 2441d4a vs current) at the
two other realistic widths — **k=31/m=15 → uint64 (width 8)** and **k=63/m=31 → __uint128 (width 16)**
— on both genomes (Mskmer/s, chr21 / celegans):

| k / m | width | original | current | **speedup** |
|---|---|--:|--:|--:|
| 21 / 11 | 4 (uint32)   | 2.20 / 2.28 | 3.51 / 3.59 | **+59.5 % / +57.5 %** |
| 31 / 15 | 8 (uint64)   | 1.60 / 1.68 | 2.54 / 2.57 | **+58.8 % / +53.0 %** |
| 63 / 31 | 16 (__uint128) | 0.52 / 0.47 | 0.90 / 0.85 | **+73.1 % / +80.9 %** |

Two takeaways: (1) the digest is **identical original-vs-current at every width** (the optimizations
stay output-bit-identical for uint32/uint64/__uint128, not just the k=21 case) — these reference
digests are pinned in `benchmark/results/reference/producer_digest.tsv`; (2) the speedup **holds and
even grows with width** — at width 16 it reaches +73–81 %, because the word-parallel `reverse_complement`
and the branch eliminations save more when each per-element op is heavier and the super-k-mer
(2k-m = 95 nt at k=63) is longer. So the wins are not a width-4 artifact; they generalize, with k=31
(the most realistic benchmark k) gaining ~+53–59 %.

Verify any width with e.g. `K=31 M=15 bash benchmark/scripts/producer_bench.sh` (digest-gated).

## What landed

**Ring-buffer modulo → power-of-two mask** (`lib/include/io/Skmerator.hpp`). The per-position super-
k-mer ring buffer was indexed `idx % m_buffer_size` with `m_buffer_size = 2k-m` (e.g. 31), a *runtime*
divisor the compiler cannot strength-reduce — it emitted a hardware `div` on the per-base hot path
(objdump: 202 `div` in the producer functions, several executed per base over 40–100 Mbp). The buffer
is now allocated at `m_buffer_capacity = std::bit_ceil(2k-m)` (next power of two) and indexed
`idx & m_buffer_mask`. Output-identical by construction: each absolute index still maps to one fixed
slot, and the live window (`≤ 2k-m ≤ capacity`) never aliases. After the change the producer functions
contain **0** `div`. This is the `perf stat` prediction realized: the loop was frontend/instruction-
throughput bound (IPC ~2.2, 0.36 % branch-miss, low backend), and removing a multi-cycle `div` per
base cut instructions on the critical path.

**Word-parallel `reverse_complement`** (`lib/include/io/Skmer.hpp`, rows 2 then 4). `canonicalize`
runs `reverse_complement` on every yielded super-k-mer (the per-yield path was ~37 % of the producer,
and RC alone was the single top symbol at ~18 %), and it rebuilt the RC slot-by-slot through the
**branch-heavy `pair` `>>`/`<<` operators**. The interleaving is built so that RC is a *within-lane*
prefix↔suffix swap plus a per-nucleotide complement — i.e. "swap the two 2-bit halves of every nibble
and XOR each with 0b10". That is now done **word-parallel** in a few ops per word
(`((w & 0x3..) << 2) | ((w & 0xC..) >> 2)` then `^ 0xA..`, masks precomputed per kuint), instead of a
~15-iteration per-lane loop. The one asymmetric case — the central self-mapping lane for odd 2k-m
(which is exactly k=21/m=11) — is fixed up explicitly, and a final `& m_mask` clears the tail the loop
left zero. Output-identical: digest unchanged and the full suite (strand-invariance, bug05/06/07,
framing) green. (An intermediate per-lane direct-word version, row 2, preceded this.)

**Precomputed `mask_absent_nucleotides`** (`lib/include/io/Skmer.hpp`). The "fill absent flank slots
with 0b11" step ran two per-call loops, and it is on the per-yield path *twice* (once in `operator++`,
once inside every `reverse_complement`). The flank fill depends only on the prefix/suffix size, so it
is now a pair of O(1) table lookups (`m_absent_pref_masks[pref] | m_absent_suff_masks[suff]`),
precomputed in the manipulator ctor next to `m_pref_masks`. Same bits set → output-identical.

**Cached minimizer φ-rank** (`lib/include/io/Skmerator.hpp`). `minimizer_rank` (a `phi()`
xorshift-multiply over the m-mer) is computed once per base in `operator++`; on every out-of-context
event `recompute_minimizer` re-scanned the window and recomputed `phi` for each slot again
(~13 % of the producer). A slot's minimizer bits don't change once its candidate is built (only the
prefix/suffix sizes do), so the per-base rank is now cached in a buffer-parallel array and
`recompute_minimizer` reads it instead of re-running `phi`. Output-identical (same ranks).

Also kept (correctness-neutral cleanup): **`m_skmer_orientation` `std::vector<bool>` →
`std::vector<uint8_t>`**. The bit-packed `vector<bool>` looked large in the `-fno-inline` profile, but
in the optimized Release build it was **perf-neutral** (measured). It is kept because it removes a
known anti-pattern the codebase already avoids elsewhere (`canonical_pieces` uses `vector<uint8_t>`
"not vector<bool>"), at zero cost and output-identical — not because it is a speedup.

## What was tried and rejected

- **`pair` / `Skmer` made trivially copyable (`= default` copy/move).** Semantically identical, digest
  green — but **regressed** ~13 % on both genomes (the defaulted copy also moves the `m_pad` tail bytes
  the hand-written copy skipped, and perturbed the move paths). **Reverted.** Lesson: the hand-written
  member-wise copies are deliberately leaner than the defaults here; do not "modernize" them.

## Method note (why the first profile misled, and how it was corrected)

The initial profile used the `-fno-inline` Profile build so `perf` could attribute per-function. That
**inflated inline-able ops** (`pair::operator[]`, the `pair` copy ctor, the shift operators) and the
`vector<bool>` bit-iterators — which is why the "obvious" `vector<bool>` fix turned out neutral in
Release. The decisive signal came from a **`-O3 -g` (inlining ON) build profiled by source line** plus
**`objdump` for the actual `div` instructions**: that pointed at the real, Release-present cost (the
modulo), which the bench then confirmed. Takeaway for the next iteration: validate every candidate on
*Release* throughput (chr21 then celegans), never on the Profile %.

## Remaining levers (not taken — diminishing returns / higher risk)

After the modulo fix, the Release source-line profile is dominated by the **`pair` shift operators**
(`Skmer.hpp operator<<=`/`operator>>=`, the cross-word transfer, ~20 % self) inside the per-base
`add_nucleotide`, and the **per-yield `canonicalize`/`reverse_complement`** (`finalize_and_yield`
~37 % children). Both are in `Skmer.hpp`, **shared with the query path** and exactness-critical:

- The per-base `pair` shift operators themselves (`operator<<=`/`operator>>=`) still carry the 4-way
  branch ladder; a single-word fast path needs a branch that stays correct for every backend width
  (at k=21/m=11 the 62-bit value genuinely spans both 32-bit words, so this is not a free win).
- `minimizer_is_ambiguous` runs per yield (O(L) decode + scan); at odd m it can only fire on a
  repeated minimizer, so the iterator could often skip it from state it already tracks.

These are smaller or riskier than the two wins above; the harness (digest + suite + bench) is in place
to attempt them safely. Iteration continues while gains stay clear on celegans.

## Iteration log — k=31 / k=63 phase (median-of-9)

The wins above were developed at k=21/m=11 (uint32). This phase targets the two realistic benchmark
points the construct path actually uses — **k=31/m=15 (uint64)** and **k=63/m=31 (__uint128)** — and
measures with `benchmark/scripts/producer_median.sh`: N independent process launches per genome,
**median** Mskmer/s with [min..max] spread (so a delta is judged against run-to-run noise), digest
gated on every run. Dev signal chr21, verdict celegans.

**Baseline @ b8f2780 (median-of-9, Mskmer/s):**

| combo | median | [min..max] |
|---|--:|--:|
| chr21 k=31    | 2.640 | [2.620..2.700] |
| celegans k=31 | 2.650 | [2.490..2.680] |
| chr21 k=63    | 0.920 | [0.900..0.930] |
| celegans k=63 | 0.860 | [0.860..0.870] |

Run-to-run noise is ~±2-3% (chr21), with occasional low outliers on celegans; treat <~3% as noise.

**Profile (build-relg `-O3 -g`, chr21, self%):** at **k=63 `minimizer_is_ambiguous` dominates (~49%)** —
it runs once per yield and m=31 makes its m-length loops (decode, center, 2×`rc_mmer`, fwd, roll)
expensive — with `operator++` ~18 %, `add_nucleotide` ~16 %. At k=31 the load is balanced:
`operator++` ~22 %, `add_nucleotide` ~21 %, `minimizer_is_ambiguous` ~19 %, `compute_new_candidate` ~13 %.

| # | idea | chr21 k31 | cel k31 | chr21 k63 | cel k63 | status |
|--:|---|--:|--:|--:|--:|---|
| B10 | `minimizer_is_ambiguous`: drop redundant 2nd `rc_mmer`+`min` — `{canon,rc_canon}`≡`{center,rc_center}` since `rc_mmer` is an involution, so compare against the pair directly (one fewer O(m) loop/yield) | 2.64→2.76 **+4.5 %** | 2.65→2.75 **+3.8 %** | 0.92→0.94 **+2.2 %** | 0.86→0.87 **+1.2 %** | **committed** — digest identical (chr21+celegans, both k), ctest green |
| A1 | `add_nucleotide`: scalar central-nucleotide transfer — the suffix↔prefix central nucleotide is one 2-bit value that never straddles a word, so precompute its word/offset once and move it with a per-word read/OR instead of 4 branchy full-`pair` shifts/base. Helps the wide (`__uint128`) path most | 2.76→2.80 +1.4 % (noise) | 2.75→2.75 ±0 | 0.94→0.98 **+4.3 %** | 0.87→0.90 **+3.4 %** | **committed** — digest identical, ctest green; clear k63 win, no k31 regression |
| C1 | `phi`: compute the mixer in `uint64_t` when 2m≤64 even if `kuint` is wider — only the low 2m bits survive the mask, so a 64-bit `mul` is bit-identical to the `__uint128` one and replaces a full 128-bit multiply per base. Pure k63 win (compile-time dead for k≤31, so that codegen is unchanged) | (same binary) | (same binary) | 0.98→1.02 **+4.1 %** | 0.90→0.95 **+5.6 %** | **committed** — digest identical, ctest green |

Best-known after C1 (median-of-9, Mskmer/s): chr21 k31 **2.80** · cel k31 **2.75** · chr21 k63 **1.02** · cel k63 **0.95**.
Cumulative vs baseline b8f2780: chr21 k31 **+6.1 %**, cel k31 **+3.8 %**, chr21 k63 **+10.9 %**, cel k63 **+10.5 %**.
(k31 unchanged by C1 — its `sizeof(kuint)>8` branch is compile-time dead at width 8; the ±noise wobble is the same binary.)

## Reproduce

```bash
CC=clang-18 CXX=clang++-18 cmake -S . -B build-timing -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=ON
cmake --build build-timing -j --target sskm-produce sklib-tests
ctest --test-dir build-timing --output-on-failure        # full suite incl. skmerator_digest
bash benchmark/scripts/producer_bench.sh                 # throughput + digest gate (chr21, celegans)
# A/B: git stash push -- lib/include/io/Skmerator.hpp ; rebuild ; bench ; git stash pop ; rebuild ; bench
```
