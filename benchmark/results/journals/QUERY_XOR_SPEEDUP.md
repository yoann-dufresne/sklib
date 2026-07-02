# Streaming query — XOR multi-column short-circuit (iteration journal)

**Goal.** Speed up the streaming file query (`sskm query -i … -t N`) at the wide store widths
(**k=63 → __uint128 (16 B)**, k=127 → kuint256), where the per-column dichotomic search
(`search_kmers_in_span_into`, `VirtualSkmer.hpp`) is the hot loop. Idea (user's): use XOR to compare
**several k-mers at once within one probed super-k-mer**. Output must stay **byte-identical** (same
per-column presence flags), and no width may regress — the narrow (uint32/uint64) path must be untouched.

## The hot loop

`search_kmers_in_span_into` is a single dichotomic search over a bucket's sorted records, but at **each
probe** `list[mean]` an inner loop compares **all active columns** (k-mer positions of the query
super-k-mer, up to k−m+1) via `kmer_compare(query, list[mean], col)` — 2 masks + 2 word-compares per
column. `query` and `list[mean]` are **constant** across that loop; only the mask `kmer_masks[col]`
varies, and the minimizer + central nucleotides sit in the **high bits, shared by every column's mask**.

## Setup (held constant)

- **Machine.** Intel Core Ultra 7 165H, Linux 6.17, pinned to one core (`taskset -c 5`), `-t 1`.
- **Build.** `clang++-18` (18.1.3), `-DCMAKE_BUILD_TYPE=Release` (`-O3 -march=native`), `PORTABLE_BUILD=OFF`.
  A/B via a compile-time `-DSKLIB_QUERY_OPT=N` selector (0 = original loop, 1 = mask-cache, 2 = +XOR), so
  the baseline binary is the untouched loop.
- **Data.** Indexes at 4096 buckets (VSKMER_5): **chr21** (~17 KB/bucket, L1), **celegans** (~44 KB, L2),
  **chr1** (~102 KB, L2+/DRAM) at k63/m31 → __uint128; chr21 at k31/m15 → uint64 (narrow control). Stream
  queries: 40 000 records × 300 bp (`e2e_helpers.py simreads`), present-fraction p ∈ {0, 100} (p=0 = fully
  substituted reads, absent; p=100 = clean genome substrings). ~9.5 M k-mers/query (k63).
- **Protocol.** Two measurements: (a) **end-to-end wall** — interleaved **ABBA** order (cancels the
  turbo-warm 2nd-position bias that swamped a naive A,B loop, worst on the big chr1 index), page-cache
  pre-warmed, **min of 15**; (b) **search-loop time** — an env-gated internal `steady_clock` timer around
  `search_kmers_in_span_into` (isolates the optimized loop from process startup / index load / FASTA
  producer / output IO, ~40 % of wall), **min of 11**. Noise floor from base-vs-base A/B.
- **Correctness gates.** Byte-identical `-o` output vs the OPT=0 baseline (`cmp`, 18 configs); full `ctest`
  (DEBUG, OPT=2); `-t 1` vs `-t 8` determinism.

## Measure-first — is the loop even ALU-bound? (kill-gate)

A prior note flagged "query search memory-bound at big N", so before implementing, an env-gated instrument
(`SKLIB_QUERY_STATS`, since removed) measured the *algorithmic* opportunity per probe:

| k63 stream | C̄ (active cols/probe) | fire-rate (Dhi>Dlo) | compares removed by XOR |
|---|--:|--:|--:|
| p100 | 9.2 | 0.72 | **83.5 %** |
| p50  | 8.3 | 0.73 | **92.3 %** |
| p0   | 7.7 | 0.73 | **100 %** |

C̄ ≈ 8-9 (not the ~1-2 that would kill it), fire-rate ~0.72, and XOR would elide **83-100 % of the
per-column `kmer_compare`s**. Verdict: **GO** (the compares are real work worth removing; whether that
moves wall time is settled below).

## Design

**V1 — query-mask cache** (the `SetOperations.hpp` WIDE_MERGE lesson). `query & kmer_masks[col]` is
constant over the whole search, so precompute it once per query into `qk[col]`; the inner compare becomes
`qk[col].compare3(masked_kmer(list[mean], col))` — 1 mask + 1 three-way compare instead of 2 masks + `<`/`>`.

**V2 — XOR shared-high-bit short-circuit** (on top of V1). Per probe, `D = query.m_pair ^ list[mean].m_pair`
(one XOR). The active in-range columns are the contiguous ascending slice `[lo,hi)` of
`positions_to_search` (built ascending; the valid-window filter is an interval). With `Ma,Mb` the two
extreme columns' masks, `S = Ma & Mb` (bits shared by **all** active columns) and `V = Ma ^ Mb` (divergent
zone) by the monotone sliding-window mask geometry. If `(D & S) > (D & V)` the highest differing bit of D
is in the shared region → **every active column compares the same direction and none can be equal**, so a
single `(query&S).compare3(list&S)` resolves them all — replacing up to C̄ `kmer_compare`s with one.

- **The `Dhi>Dlo` guard is mandatory, not an optimization.** The interleaving (prefix in the low 2 bits,
  suffix in the high 2 bits of each nibble) can place a column's own bits (in V) *above* shared bits of S,
  so the naive `D & S != 0` would mis-order that column. Because S and V are disjoint bit-sets,
  `Dhi>Dlo ⟺ the top set bit of D∩(S∪V) lies in S`, which is exactly the common-direction condition.
- **Byte-identity.** Fire ⇒ `D∩S≠0` ⇒ `D∩M_c≠0` ∀ active c ⇒ no column is "found"; each gets the identical
  sign-based boundary update + collapse check as the per-column loop; columns outside `[lo,hi)` are
  untouched; `current_priority_offset` (only moved by the found branch) is unchanged. No fire → fall through
  to the V1 loop verbatim.
- **Wide-only gate.** All of V1/V2 is under `if constexpr (sizeof(kuint) >= 16)`; the narrow (uint32/uint64)
  instantiation is the original loop verbatim (mirrors WIDE_MERGE_MIN_STORE_BYTES).

## Correctness

- **Byte-identical** `-o` output vs OPT=0 across 18 configs (chr21/celegans × k63/k31 × p{0,50,100} ×
  -t{1,8}): V1 and V2 both `cmp`-equal. The `Dhi>Dlo` guard and the monotone-mask containment hold
  empirically at every width tested.
- **`ctest` 213/213** (DEBUG, default OPT=2).
- **Narrow untouched**: the search-loop timer reports **−0.0 %** at k31 — the `if constexpr` leaves the
  narrow codegen byte-for-byte identical (no layout regression).

## Results

Noise floor (base vs base, ABBA min-15): cele k63 +2.2 %, chr1 k63 +1.4 % — so wall gains below clear the band.

**Wall (end-to-end, base vs xor, min-of-15):**

| k63 index (regime) | p=100 | p=0 (absent) |
|---|--:|--:|
| chr21 (L1)         | **+5.0 %** | — |
| celegans (L2)      | **+7.2 %** | **+10.3 %** |
| chr1 (L2+/DRAM)    | **+7.5 %** | +6.0 % |
| chr21 **k31 narrow** | +1.7 % (≈ noise) | — |

**Search loop (internal timer, min-of-11):**

| k63 index | p=100 | p=0 |
|---|--:|--:|
| chr21    | +12.1 % | — |
| celegans | +15.2 % | +15.7 % |
| chr1     | +12.9 % | +14.1 % |
| **k31 narrow** | **−0.0 %** | — |

**V1 vs V2 split** (search loop, celegans k63 p100): mask-cache alone **+8.4 %**; XOR adds **+10.4 %** on
top of it → **+15.2 %** combined. Both contribute; the XOR short-circuit is worth its complexity.

The search loop is **+12-16 % faster** at every wide cache regime and grows with the absent fraction
(p=0 ⇒ 100 % fire). End-to-end that is **+5-10 %**, diluted because the search is only ~60 % of query wall
(the rest — FASTA producer, index load, output — is unchanged). The gain **survives the memory-bound
regime** (chr1, L2+/DRAM buckets): even there the C̄ compares against an already-loaded `list[mean]` are
enough ALU to matter, so removing them helps — the "memory-bound" kill risk did not materialize.

## Conclusion — **COMMITTED** (default `SKLIB_QUERY_OPT 2`)

XOR-comparing all columns of a probed super-k-mer at once (short-circuited on the shared high-bit region)
speeds the wide streaming-query search loop **+12-16 %** (**+5-10 % end-to-end**), byte-identical, with the
narrow path provably untouched. Enabled by default for wide records; set `-DSKLIB_QUERY_OPT=0` to disable.

**Remaining lead (not this journal).** The search win is capped because `update_searchable_positions` runs
`O(active columns)` **every probe** — with the compares gone it is now the per-probe hotspot. Making the
active-column set incremental (it only shrinks as columns resolve) would compound with this change; that is
the next lever for the streaming query.
