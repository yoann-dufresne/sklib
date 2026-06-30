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
  noise); inputs page-cache-warm. Driver: `benchmark/scripts/microbench/union_bench.sh`.
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

### #1 — column-offset fast-path for `sort_column` — **COMMITTED** (`7cf06f3`)

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

### #2 — hash join in `get_candidate_overlaps` — **COMMITTED** (`2ac3914`)

**Mechanism.** After #1, `get_candidate_overlaps` was the top hotspot (~54 % k31 / ~59 % k63): it
joined two adjacent columns on their shared overlap (k-1)-mer by **sorting** the right column's keys
(`O(R log R)`) then a **`lower_bound`** per left key (`O(L log R)`) — both hammering the wide
`pair::operator<` (15 % / 26 % self). Replaced with an **array-based chaining hash join**: build
`head[]`/`next[]` over the right keys (`O(R)`), probe each left key (`O(L)`) → `O(R+L)`, no sort, no
binary search. Two flat `int64` arrays reused across columns (`O(R)` memory), so unlike the original
per-key `unordered_map` it does **not** blow up peak RAM on repeat-rich buckets. A width-agnostic
`hash_kpair` folds each store word to 64 bits (uint64/__uint128/kuint256). Both `colinear_chaining`
and `greedy_chaining` **sort** their input, so the (different) candidate emit order selects the **same**
chain — output **byte-identical** for construction *and* set-ops (`VirtualSkmer.hpp`).

The `GetCandidateOverlap*` unit tests pinned the candidate *order*; that is now an implementation
detail, so they were made order-insensitive (compare as sets) — `tests/km/VirtualSkmer.cpp`.

**Correctness.** Output **byte-identical** (sha256 match k31 & k63 J0.5); bits/kmer unchanged
(23.89 / 24.24); 213/213 unit tests (construction golden digests included); KMC cross-check PASS.

**Result (interleaved A/B, idea#1 → idea#2):**

| config | idea#1 s | idea#2 s | delta | band |
|---|--:|--:|--:|--:|
| chr21 k31 J0.5 | 4.094 | 2.569 | **−37.3 %** | ±0.20 % |
| chr21 k63 J0.5 | 5.315 | 2.834 | **−46.7 %** | ±0.11 % |

Cumulative vs the original baseline (`135283b`): chr21 k31 J0.5 4.72 → 2.57 s (**−46 %**),
k63 J0.5 7.40 → 2.83 s (**−62 %**).

**Full 12-config sweep, cumulative vs baseline** (idea#1+#2; all 12 `verify=PASS`, byte-identical):

| dataset | k | baseline s (J.1/.5/.9) | idea#2 s (J.1/.5/.9) | gain |
|---|--:|--:|--:|--:|
| chr21 | 31 | 6.76 / 4.72 / 3.40 | 3.68 / 2.60 / 1.82 | **−45 to −46 %** |
| chr21 | 63 | 11.51 / 7.40 / 5.19 | 3.98 / 2.84 / 2.00 | **−61 to −65 %** |
| celegans | 31 | 19.33 / 13.40 / 9.84 | 10.06 / 7.28 / 5.14 | **−46 to −48 %** |
| celegans | 63 | 32.46 / 23.50 / 17.78 | 10.95 / 7.88 / 5.66 | **−66 to −68 %** |

k31 (`uint64`) ~−46 %; k63 (`__uint128`) **~−63 to −68 %** — the union is now ~2× (k31) to ~3× (k63)
faster than `135283b`.

### #3 — skip the redundant chaining sort (`is_sorted` guard) — **COMMITTED** (`2f379cc`)

**Mechanism.** After #2 the colinear-chaining of each column pair was the top recompaction cost; its
first step `std::sort`s the candidate overlaps by `(first asc, second desc)`. But #2's hash join
already emits them in exactly that order — the probe walks left indices ascending, and for each left
key the matching rights come out of the chain head-first, i.e. in **descending** index order (insertion
was ascending). So the sort is a no-op. Guard it with an `O(n)` `std::is_sorted` check in **both**
`greedy_chaining` (set-ops, `ColinearChaining.hpp`) and `colinear_chaining` (construction,
`ColinearChaining.cpp`): the already-ordered `O(n log n)` sort is skipped, the result is identical.
Byte-identical for set-ops *and* construction (golden digests unchanged).

**Correctness.** Output **byte-identical** (sha256 match k31 & k63 J0.5); bits/kmer unchanged; 213/213
unit tests (construction golden digests included); KMC cross-check PASS.

**Result (interleaved A/B, idea#2 → idea#3):**

| config | idea#2 s | idea#3 s | delta | band |
|---|--:|--:|--:|--:|
| chr21 k31 J0.5 | 2.627 | 2.249 | **−14.4 %** | ±0.45 % |
| chr21 k63 J0.5 | 2.851 | 2.430 | **−14.7 %** | ±0.40 % |

Cumulative vs `135283b`: chr21 k31 J0.5 4.72 → 2.25 s (**−52 %**), k63 J0.5 7.40 → 2.43 s (**−67 %**).

### #4 — skip redundant `get_skmer_of_kmer` in the recompaction — **REVERTED (regression)**

**Idea.** On the set-op path `col` already holds single-k-mer skmers valid at exactly their column, and
`get_skmer_of_kmer(S, column_pos)` is the identity on such an `S`, so `merge_LList_column` re-derives
what it already has. Replacing those calls with the element directly was **byte-identical** (sha256
match, 213 tests, golden digests) but **regressed k63**: interleaved A/B at chr21 J0.5 gave k31
**−3.5 %** but k63 **+1.0 %**. Tried three variants — return by value (k63 +1.02 %), by `const&` via a
reused scratch (+1.24 %), and a compile-time width gate `if constexpr (sizeof(kuint) <= 8)` so the wide
store is excluded (+1.30 %). The width-gated one *still* regressed k63, which is the tell: merely adding
the `frame` indirection + the `single_kmer_enum` parameter to the shared `merge_LList_column` hot loop
perturbs the `__uint128` codegen unfavourably, independent of the skip. The k31 win is unreachable
without a k63 regression in this shared loop → reverted per the no-regression rule (`git checkout`).

| # | idea (file) | mechanism | result (k31 / k63, J0.5) | status |
|--:|---|---|---|---|
| 1 | column-offset fast-path for `sort_column` (`SetOperations.hpp`, `VirtualSkmer.hpp`) | skip the per-column `has_valid_kmer` full-scan; ids = contiguous block from merge's per-column counts | **−12.0 % / −27.6 %** (byte-identical) | committed `7cf06f3` |
| 2 | hash join in `get_candidate_overlaps` (`VirtualSkmer.hpp`) | array-based chaining join, `O(R+L)` vs `O((R+L)·log R)`; drops the per-column sort + `lower_bound` | **−37.3 % / −46.7 %** (byte-identical) | committed `2ac3914` |
| 3 | skip redundant chaining sort (`ColinearChaining.hpp/.cpp`) | hash join already emits `(first asc, second desc)`; `is_sorted` guard skips the `O(n log n)` sort | **−14.4 % / −14.7 %** (byte-identical) | committed `2f379cc` |
| 4 | skip redundant `get_skmer_of_kmer` in `merge_LList_column` | identity on set-op single-k-mer skmers | **−3.5 % / +1.0 %** (k63 regressed, all 3 variants) | reverted-regression |

## State of play after #1–#3 (and the regenerated backlog)

**Cumulative:** chr21 k31 J0.5 4.72 → 2.25 s (**−52 %**), k63 J0.5 7.40 → 2.43 s (**−67 %**); celegans
similar (k31 −54 %, k63 −71 %). Throughput ~10 → ~21 Mkmer/s (k31), ~6 → ~20 Mkmer/s (k63). Output
**byte-identical** to `135283b` end-to-end (CLI `sskm setop --op union` sha256 match), so bits/kmer,
packing and the k-mer set are all unchanged. The whole `set_union` is still **87–91 % re-compaction**.

**Profile now (perf self-time, k31 / k63), all within the re-compaction:**

| function | self % | note |
|---|--:|---|
| `get_candidate_overlaps` (hash join) | 31 / 28 | already `O(R+L)`: extract + `hash_kpair` + chain walk |
| `merge_LList_column` | 21 / 23 | `add_kmer` (chain) + `get_skmer_of_kmer` (new vskmers); **k63-codegen-sensitive** |
| `merge_columns` (CSR + 2-cursor merge + CollectSink) | 18 / 22 | `kmer_compare` + `get_skmer_of_kmer` building `col` (necessary) |
| `greedy_chaining` | 16 / 16 | the LIS (patience sort) — sort already skipped by #3 |
| `build_column_csr` | 7 / 4 | 2-pass counting sort of A/B by column |

**Discarded (reasoned/measured), not committed:**
- *Drop the `get_candidate_overlaps` sort* (pre-#2): refuted — a probe showed `right_keys` is genuinely
  unsorted on big buckets (the overlap (k-1)-mer is not monotone in the k-mer order).
- *Skip `get_skmer_of_kmer` in `merge_LList`* (#4): byte-identical, k31 −3.5 % but k63 +1.0–1.3 % in 3
  variants — the shared wide-store hot loop is codegen-fragile. Reverted.

**Remaining leads — all sub-2 %, several with k63 codegen risk; profile-bounded, not yet worth a commit:**
- Precompute `~kmer_masks[c] & m_mask` (the get_skmer_of_kmer fill) per column: removes ~4 ops/call but
  it is 1 of several ops there (~<1 %, below the k63 band; touches the fragile path).
- `greedy_chaining` buffer reuse (`thread_local` parent/tail): saves only the ~0.4 % malloc.
- `build_column_csr` single-pass / cache `get_valid_kmer_bounds`: the re-read is cache-resident → ~neutral.
- Contiguous per-column **record** layout (scatter records, not indices) for cache-friendlier
  `merge_columns`: plausible but for k63 it 9× the scatter write traffic — likely a wash or worse.

The three algorithmic wins (eliminate the per-column scan, replace the join sort+search with a hash
join, skip the now-redundant chaining sort) captured the high-value monothread speedup; the
re-compaction is now memory/compute-bound on irreducible work. Further *union* gains look incremental.

## Extending to intersection & diff (the merge phase)

`intersection`/`difference` share `materialize_setop` with `union`, so #1–#3 already sped them up. But
the harness was union-only; it now takes `--op union|intersection|diff` (same for `union_bench.sh` /
`union_ab.sh` via `OP=`). The phase timer reveals a very different balance — union always materialises
*everything* so it is recompaction-bound, but intersection/diff materialise a **subset**, so the
**merge phase** (`merge_columns` + `build_column_csr`, which still scans *all* of A and B) dominates:

| op (chr21 k31) | J0.1 | J0.5 | J0.9 |
|---|--:|--:|--:|
| intersection | merge **70 %** | merge 39 % | recompact 74 % |
| diff | recompact 56 % | merge 56 % | merge **77 %** |
| union | recompact 74 % | recompact 72 % | recompact 76 % |

So intersection@low-J and diff@high-J are **merge-bound** — a hotspot untouched by #1–#3. `perf` put
`build_column_csr` at **44 %** self (k31 diff J0.9): its pass-1 counting loop increments `off[c+1]` for
every column a record spans, i.e. `O(total k-mers)` since A/B records are full super-k-mers.

### #5 — difference-array `build_column_csr` pass-1 count — **COMMITTED** (`120a2c6`)

**Mechanism.** A record valid at columns `[s,hi]` adds 1 to each per-column count — a **range update**,
applied in `O(1)` as `+1` at `s`, `−1` at `hi+1`, then one prefix sum, instead of the `O(hi−s+1)` loop.
Turns pass 1 from `O(total k-mers)` to `O(records)`; the scatter pass and the resulting `off`/`idx` are
unchanged, so **byte-identical** for every set op (`build_column_csr` in `SetOperations.hpp`). Shared,
so it speeds up **all** setops (and the `*_size` counts).

**Correctness.** Byte-identical (sha256) for union/intersection/diff at k31 & k63 (8 configs incl. the
merge-bound extremes); 213/213 unit tests; KMC cross-check PASS.

**Result (interleaved A/B vs #1–#3, chr21):**

| op / config | ref s | idea#5 s | delta | band |
|---|--:|--:|--:|--:|
| diff k31 J0.9 (merge 77 %) | 0.330 | 0.267 | **−19.1 %** | ±0.60 % |
| diff k63 J0.9 | 0.372 | 0.330 | **−11.1 %** | ±1.33 % |
| intersection k31 J0.1 (merge 70 %) | 0.889 | 0.802 | **−9.8 %** | ±0.61 % |
| intersection k63 J0.1 | 1.000 | 0.948 | **−5.2 %** | ±0.68 % |
| union k31 J0.5 | 2.215 | 2.080 | **−6.1 %** | ±0.36 % |
| union k63 J0.5 | 2.433 | 2.413 | −0.8 % | ±0.32 % |

Wins every op/width, no regressions — biggest on the merge-bound cases; even union gains (its merge is
~20–25 %). The merge phase's other half (`merge_columns`' `kmer_compare`: a masked wide-`pair` compare,
~60 % self at k63 intersection) is the next frontier but looks largely irreducible.

### #6 — pre-masked contiguous keys in `merge_columns` — **REVERTED (regression + RAM)**

**Idea.** `merge_columns`' compare reads each record through the `idx` indirection (random `A[idx]`)
and re-masks it (`& kmer_masks[c]`) every comparison. Prototype: have `build_column_csr` also emit a
**parallel pre-masked-key array** per column (the record's k-mer at `c`, byte-identical scatter order,
so per-column order is conserved → byte-identical output), and let the merge compare those contiguous,
already-masked keys directly — touching the record only for *kept* k-mers (the sink). Byte-identical
(7 ops/widths), 213 tests.

**Why it fails (measured).** The premise is wrong: `keys` holds one entry per **(record × column)**
incidence (≈ total k-mers), so it is *larger* than the `A` bucket it linearises (one entry per record).
At k63 (32 B/key) the "sequential" scan is over a bigger, cache-pressuring array than the small,
cache-resident bucket the random `A[idx]` already hits. Interleaved A/B vs #5: intersection **−8.5 %**
(k31) / −4.0 % (k63) but **diff +2.3 % / +13.2 %** and union +1.0 % / +1.2 % — regresses the op it
targeted. Peak RSS (celegans k63): diff 62→111 MB (**+79 %**), union 111→160 MB (+44 %) — the same
RAM blow-up the original code dropped its `unordered_map` to avoid. Net loss → reverted. Conclusion:
the merge's random access into a small cache-resident bucket is already near-optimal; `kmer_compare` is
at its floor.

| # | idea (file) | mechanism | result | status |
|--:|---|---|---|---|
| 5 | difference-array `build_column_csr` (`SetOperations.hpp`) | range-update pass-1 count, `O(records)` not `O(total k-mers)` | diff −19/−11 %, inter −10/−5 %, union −6/−0.8 % (byte-identical) | committed `120a2c6` |
| 6 | pre-masked contiguous keys in `merge_columns` | linearise + pre-mask the compare | inter −4..−9 % but **diff +2..+13 %**, union +1 %, **peak RAM +44..79 %** | reverted-regression |

## Combined operator `multi_setop` (the one-pass {∩, ∪, A\B, B\A})

`multi_setop` shares the merge + recompaction with the single-op path, so #2/#3/#5 already applied —
but **#1 did not**: it materialises up to four relations and re-compacts each with its own recompactor,
and those calls fed `generate_sorted_list_from_enumeration` **without** `col_offsets`, so every one of
the (up to 4) per-bucket recompactions still did the `O(n·(k-m))` `has_valid_kmer` re-scan. It was
therefore the *most* recompaction-starved set op.

### #7 — column-offset fast-path for `multi_setop` — **COMMITTED** (`b50346f`)

**Mechanism.** Each output buffer in `multi_setop` is column-grouped + per-column sorted+distinct (the
merge fans column-major into it), exactly like the single-op `col`. `MultiCollectSink` now keeps a
per-column kept count **per channel** (incremented as each k-mer is fanned out — `both`→inter&union,
`only_a`→A\B&union, `only_b`→B\A&union); `multi_setop` prefix-sums each into block offsets and passes
them to that channel's recompaction (idea #1, per output). Byte-identical — same ids, so each of the
four outputs is bit-for-bit unchanged (`MultiCollectSink` + `multi_setop` in `SetOperations.hpp`).

**Correctness.** All four outputs **byte-identical** (sha256) before/after across chr21/celegans ×
k31/k63 × J; and (as the single-op design already guaranteed) each equals the matching single-op
materialization; 213/213 tests; KMC `setop_multi_verif.sh` PASS.

**Result (interleaved A/B, idea#5 → #7, all four relations materialized):**

| config | idea#5 s | #7 s | delta | band |
|---|--:|--:|--:|--:|
| chr21 k31 J0.5 | 4.84 | 3.65 | **−24.5 %** | ±0.47 % |
| chr21 k63 J0.5 | 10.19 | 4.54 | **−55.4 %** | ±3.39 % |
| chr21 k31 J0.1 | 6.96 | 5.27 | **−24.3 %** | ±0.63 % |
| chr21 k63 J0.1 | 15.89 | 6.34 | **−60.1 %** | ±1.25 % |

k31 ~**−24 %**, k63 **−55 to −60 %** — the biggest single win in the journal, because `multi_setop` ran
the slow per-column scan **four times per bucket** and k63 (32 columns) suffered most. The harness now
takes `--op multi` (bench mode; 4 outputs, summed throughput).

| # | idea (file) | mechanism | result (k31 / k63) | status |
|--:|---|---|---|---|
| 7 | column-offset fast-path for `multi_setop` (`SetOperations.hpp`) | per-channel column counts → offsets to each output's recompaction | **−24 % / −55..−60 %** (byte-identical) | committed `b50346f` |
