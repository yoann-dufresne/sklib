# Super-k-mer producer — within-producer bottleneck breakdown

**Question.** `CONSTRUCT_SCALING_DIAG.md` established that the *producer* (Phase 1: FASTA parse →
rolling minimizer → super-k-mer assembly) is the serial floor of `sskm construct` (~⅓ of the work,
~73–78 % of wall at `-t8`) and that, inside it, FASTA parse ≈ 0 % and I/O ≈ 0.5 % — i.e. it is pure
`SkmerManipulator` compute. **What exactly, inside the producer, is the slowest?** And can a
*single-threaded* rewrite shrink that floor (the complementary lever to parallelizing Phase 1)?

**TL;DR.** The producer's time is spent almost entirely in **packed-integer (`Skmer::pair`)
manipulation** and in **`std::vector<bool>` orientation bookkeeping** — not in parsing, hashing, or
the rare ambiguous-framing path. It is **instruction-throughput / frontend bound**, **not**
branch-misprediction bound (0.36 % miss) and **not** memory bound (≤12 % backend-bound). The two
highest-value, output-preserving monothread levers are therefore (1) **drop `std::vector<bool>`** for
the per-position orientation buffer (its bit-proxy iterators are ~13 % of producer samples), and (2)
**collapse the two-word `pair` work to one machine word** for the common `2*(2k-m) ≤ 64` case (the
`pair` shift/copy/`|=`/word-access family is >40 % of producer samples; for k=21/m=11 only 62 of the
64 pair bits are ever used, yet every op pays the generic two-word, multi-branch path).

## Method

- **Binary:** `build-profile/bin/sskm-produce` — the isolated producer (`app/produce/`,
  `FileSkmerator` only, no bucketing/sort/writer). Profile build (`-O3 -ggdb3 -fno-inline`, clang-18).
- **Profiling:** `benchmark/scripts/producer/producer_perf.sh` — `perf record --call-graph dwarf` (flamegraph +
  flat report) and `perf stat` (IPC, branch/cache misses, top-down) on chr21 (40.1 Mbp) and celegans
  (100.3 Mbp), k=21 m=11, single thread.
- **Throughput baseline:** `benchmark/scripts/producer/producer_bench.sh` (Release build), best-of-5, warm page
  cache. chr21 **2.44 Mskmer/s** (2.77 s, 6 750 903 skmers); celegans **2.53 Mskmer/s** (6.62 s,
  16 754 510 skmers). RSS ~53 MB (the input is streamed; no whole-genome buffering).
- **Caveat (read before trusting the absolute %s):** the Profile build is `-fno-inline` so that tiny
  functions get their own frames and `perf` can attribute them. This **inflates exactly the
  inline-able operations** (`pair::operator[]`, the `pair` copy ctor, the `pair` shift operators),
  which in the shipped Release build largely inline. So the per-symbol percentages below are a *map of
  where the work lives by family*, not the Release time split. The robust, build-independent signals
  are: (a) the **family composition** (pair manipulation ≫ everything; vector<bool> is a large
  separable block; reverse_complement/canonicalize run per yield); and (b) the **`perf stat`
  microarchitecture verdict**, which is structural. Each candidate below is validated by the *Release*
  throughput it actually delivers (chr21 signal, celegans verdict), never by the Profile %.
- celegans `perf record` landed on the `cpu_atom` (E-core) PMU and undersampled (63 samples) — its
  flat report is only a qualitative cross-check; **chr21 (cpu_core, 41 K samples) is authoritative.**

## Results — chr21 flat profile (cpu_core, 41 K samples, self-time families)

Aggregated by family (self %, summed across the family's symbols):

| Family | ~self % | Symbols |
|---|--:|---|
| **`pair` word access** | ~19 % | `std::array<unsigned int,2>::operator[]` (mut + const) |
| **`pair` shifts/copies/or** | ~40 %+ | `operator<<=` 9.3, copy-ctor 9.3, `operator>>=` 8.3, `operator\|=` 7.9, `operator<<` 3.1, `operator>>` 2.6 |
| **`std::vector<bool>` orientation** | ~13 % (children) | `vector<bool>::operator[]`, `_Bit_iterator::operator[]`, `operator+`, `_M_incr`, `begin` |
| **per-yield canonical** | — | `reverse_complement` 3.4 self / `canonicalize` 9.9 children, `mask_absent_nucleotides` 1.3 |
| **per-base assembly** | — | `add_nucleotide` 1.7 self, `compute_new_candidate_skmer` 2.9 self |
| **ambiguous framing (rare-ish)** | — | `canonical_pieces` 4.8 children, `choose_kmer_minimizer` 4.7 children |

Driver call tree: `operator++` (32.9 % children) → `compute_new_candidate_skmer` + `add_nucleotide`
(per base) and `finalize_and_yield` (13.8 % children, per yield) → `canonicalize` →
`reverse_complement`. The minimizer permutation `phi` does **not** surface as a hot symbol — the
rolling-minimizer order key is cheap relative to the pair shuffling.

## Microarchitecture (`perf stat`, structural — build-inlining-independent)

| metric | chr21 | celegans |
|---|--:|--:|
| IPC (cpu_core) | 2.19 | 2.20 |
| branch-misses | **0.36 %** | **0.36 %** |
| top-down frontend-bound | 54.5 % | 54.2 % |
| top-down backend-bound | 8.6 % | 12.3 % |
| top-down retiring | 40.4 % | 40.2 % |

The hot loop is **frontend-bound** (instruction supply), with **well-predicted branches** and **little
memory stall**. So the comment in `Skmer.hpp` that the pair shift is "branch-heavy" is right that it
has many branches, but the cost is the **branch/instruction *count*** (frontend pressure), **not
mispredictions**. The lever is therefore *fewer instructions per emitted super-k-mer*: shrink the
packed-integer work and remove the bit-proxy bookkeeping — exactly what (1) and (2) below do. (Note:
`-fno-inline` itself inflates frontend-bound; treat the 54 % as directional, the 0.36 % miss and low
backend-bound as solid.)

## Proposed monothread rewrites (each gated by the digest; chr21 = signal, celegans = verdict)

Ordered by value × safety. Each must keep the producer's output **bit-identical** (digest unchanged
on chr21 *and* celegans, full `ctest` green) — they are micro-architectural, not algorithmic.

1. **`m_skmer_orientation`: `std::vector<bool>` → `std::vector<uint8_t>`** (`Skmerator.hpp`).
   *Safest, do first.* The per-position orientation buffer is a bit-packed `vector<bool>`; every
   read/write goes through `_Bit_iterator` word-index + mask arithmetic (~13 % of samples, and it does
   **not** inline away in Release — the bit math is inherent). A byte buffer makes each access a plain
   load/store. Output-identical (same 0/1 values). The codebase already prefers this elsewhere
   (`canonical_pieces` uses `std::vector<uint8_t> orient` "not vector<bool> (proxy refs)") — this just
   applies the same fix to the producer's hot buffer.

2. **Single-word `pair` fast path for `2*(2k-m) ≤ 64`** (`Skmer.hpp`). The `pair` is
   `std::array<kuint,2>` with multi-branch two-word shift/mask operators; for k=21/m=11 the high word
   is always 0, so every shift/copy/`|=`/word-access does ~2× the necessary work plus branches.
   Options, least-to-most invasive: (a) make `pair` **trivially copyable** (`= default` the copy
   ctor/assignment) so copies lower to one 8-byte move; (b) branch-light the `>>=`/`<<=` hot paths;
   (c) a compile-time/inline single-word path when the high word is provably unused. Shared with the
   query path → re-run the **whole** suite, not just the digest.

3. **`reverse_complement`/`canonicalize` per yield** (`Skmer.hpp`) — runs once per emitted
   super-k-mer, rebuilding the RC slot-by-slot. A word-level reverse-complement (parallel bit ops)
   would cut the per-yield cost. Medium risk (shared, exactness-critical).

4. **Ring-buffer modulo `% m_buffer_size`** (`Skmerator.hpp`, buffer size `2k-m`, not a power of two)
   — several runtime 64-bit modulos per base. A power-of-two-sized buffer + mask removes the divisions.
   Output-identical (pure indexing). Lower expected payoff than 1–2 but trivial and safe.

5. **`add_nucleotide` maintains both strands fully** each base (`Skmer.hpp`) — defer/lighten the
   non-selected strand. Highest risk (it is the core encoder); revisit only if 1–4 leave a clear floor.

The ambiguous-framing path (`canonical_pieces`/`choose_kmer_minimizer`, ~5–9 % children) is **not** a
priority: at m=11 (odd → no RC-palindrome minimizers) it only triggers on repeated-minimizer k-mers,
and on real genomes that is a minority; it would matter more at small even m.

## Reproduce

```bash
CC=clang-18 CXX=clang++-18 cmake -S . -B build-profile -DCMAKE_BUILD_TYPE=Profile -DWITH_TESTS=OFF
cmake --build build-profile -j --target sskm-produce
bash benchmark/scripts/producer/producer_perf.sh     # -> benchmark/results/latest/perf-producer/*.{svg,flat.txt,stat.txt}
bash benchmark/scripts/producer/producer_bench.sh    # Release throughput + digest gate
```
