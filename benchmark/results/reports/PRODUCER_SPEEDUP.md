# Super-k-mer producer — monothread speedup (iteration journal)

**Goal.** Shrink the serial producer floor (`CONSTRUCT_SCALING_DIAG.md`) with *single-threaded*,
**output-bit-identical** rewrites — the lever complementary to parallelizing Phase 1. Every change is
gated by the producer digest (`benchmark/results/reference/producer_digest.tsv`, chr21 + celegans) and
the full `ctest` suite; chr21 is the dev signal, celegans the verdict.

**Result.** One confirmed win: replacing the ring-buffer's runtime modulo with a power-of-two mask
(`Skmerator.hpp`) gives **+12.8 % (celegans) / +14.3 % (chr21)** producer throughput, output
byte-identical. Controlled back-to-back A/B (best-of-7, same machine state, `sskm-produce`):

| genome | baseline wall | opt wall | baseline Mskmer/s | opt Mskmer/s | speedup | digest |
|---|--:|--:|--:|--:|--:|:--|
| chr21    | 2.9176 s | 2.5532 s | 2.31 | 2.64 | **+14.3 %** | identical (0xfd02872bc3f3ed34) |
| celegans | 6.9088 s | 6.1176 s | 2.43 | 2.74 | **+12.8 %** | identical (0x93d267f24daa8af8) |

The digest is unchanged on both genomes (and the in-RAM golden test + full suite stay green), so the
producer's output — hence every `sskm construct` / `query` result built on it — is exactly preserved.

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

- A single-word fast path for `2*(2k-m) ≤ 64` (high word always 0) would cut the per-base `pair`
  shift/copy work, but needs a branch that stays correct for every backend width.
- A word-level (bit-parallel) `reverse_complement` would cut per-yield cost vs the current slot loop.

These are larger, riskier rewrites with uncertain payoff; the harness (digest + suite + bench) is in
place to attempt them safely if a further push is wanted. Stopping here banks a clean, verified
~13 % monothread gain.

## Reproduce

```bash
CC=clang-18 CXX=clang++-18 cmake -S . -B build-timing -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=ON
cmake --build build-timing -j --target sskm-produce sklib-tests
ctest --test-dir build-timing --output-on-failure        # full suite incl. skmerator_digest
bash benchmark/scripts/producer_bench.sh                 # throughput + digest gate (chr21, celegans)
# A/B: git stash push -- lib/include/io/Skmerator.hpp ; rebuild ; bench ; git stash pop ; rebuild ; bench
```
