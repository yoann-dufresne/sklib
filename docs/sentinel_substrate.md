# Absent-slot sentinel — investigation & result

**Status: negative result.** The idea below was explored end-to-end. It is *safe* but does
**not** achieve its goal (scan-free queries) beyond trivially small inputs. This note
records what was tried, the measurements, the root cause, and the path forward.

## The idea (proposed direction)

A stored super-k-mer that covers few k-mers leaves its peripheral flank slots unused;
construction padded them with `0b11`. The membership query
(`SortedVirtualSkmerList::query_skmer`) binary-searches the sorted list; when it probes an
entry that has **no k-mer at the searched column** (a "hole"), `find_closest_valid_skmer`
does a **linear scan** to the nearest valid entry. The proposal: rewrite the unused
("absent") bits as an *order sentinel* so the binary search gets a usable direction at a
hole and the linear scan disappears — reusing wasted bits, no extra memory, format
unchanged.

## What was built

* **`fill_absent_sentinel()`** (`lib/include/algorithms/VirtualSkmer.hpp`): clears each
  entry's absent flank slots from `0b11` down to `0` (the minimal completion), applied at
  construction (per bucket in the disk-backed path; before save in the all-in-RAM binary
  path). Branch `feat/sentinel-construction`.
* **`query_skmer_substrate()`** (same header): a binary search that navigates through
  every entry by direction alone (no `find_closest_valid_skmer`), accepting a match only
  at a valid entry. Branch `feat/sentinel-query-e2e`.
* Tests `tests/km/sentinel_substrate.cpp`; micro-benchmark `tests/km/bench_substrate_query.cpp`
  (gated by `SKLIB_BENCH`).

## Why it is *safe* (construction-only change)

The current query compares k-mers only at entries **valid** at the column — where the
absent slots are masked out by `kmer_masks[c]` — and skips holes, so it never reads the
filled bits. The fill therefore leaves query results **byte-identical**. Verified at the
CLI level: `sskm` built with the fill produces self-query and random-absent-genome outputs
identical to the pristine `dev` build on a 2 Mb genome (`k=21, m=11`); the fill only lowers
each entry's `m_pair`, keeps the list size, and the existing test suite stays green.

## Why it does **not** work (the core finding)

For a binary search to navigate *through* holes, the per-column key
`key_c(x) = x.m_pair & kmer_masks[c]` must be **non-decreasing along the whole list,
holes included**. Two facts:

1. **Valid entries are already monotone — at every scale.** For entries that *have* a
   k-mer at column `c`, `key_c` is non-decreasing. This is the invariant the current query
   relies on; it holds for all measured `k, m` and all genome sizes
   (`valid_only_viol = 0` everywhere). Good.

2. **Holes are not, beyond trivially small inputs.** The sentinel fill makes `key_c`
   monotone with holes included only for tiny lists. It degrades quickly with size
   (`k=21, m=11`, random genome):

   | genome (bp) | list entries | hole+valid violations | valid-only violations |
   |------------:|-------------:|----------------------:|----------------------:|
   | 1 000       | 165          | 0                     | 0 |
   | 10 000      | 1 673        | 1                     | 0 |
   | 100 000     | 16 721       | 555                   | 0 |
   | 1 000 000   | 171 645      | 58 037                | 0 |
   | 2 000 000   | 352 114      | 213 283               | 0 |

   The hole-aware query is consequently **incorrect at scale**: on a 2 Mb genome it reports
   ~**5.3 % false negatives** on present k-mers (106 217 / 1 999 980), while being ~1.9–2.9×
   faster than the scan-based query (≈44 vs ≈87 ns/super-k-mer present; ≈37 vs ≈106 absent).
   Speed without correctness.

**Root cause.** `kmer_masks[c]` is **not** a high-order prefix of `m_pair`; in the
minimizer-centred interleaved layout it is a *scattered* mask (the minimizer plus an
interleaved window of flank slots). So (a) sorting by `m_pair` does **not** order the
per-column projections, and (b) a hole's position in the list is fixed by the columns where
it *is* valid — at a column where it is absent, its **high-order** content is already out
of place. No absent-bit fill can repair high-order bits, so no single fill makes all
columns monotone for holes. This is structural, not a tuning problem.

**Re-sorting does not rescue it.** Filling then `std::sort` by `(m_pair, pref, suff)` on a
2 Mb genome still leaves 59 326 hole violations (because of the scattered mask) **and
breaks the existing query** (16 589 mismatches vs ground truth — the construction's order
is essential and is *not* the `m_pair` order). So the earlier assumption that an
order-preserving fill exists is moot: the list was never in `m_pair` order to begin with.

## Conclusion

The order-sentinel fill cannot deliver scan-free queries. The fill is harmless but
provides no demonstrated query benefit; `feat/sentinel-construction` is best treated as a
revert candidate unless a use is found for canonical (zeroed) padding.

## Path forward (for the query-optimization phase)

The valid-entry-per-column monotonicity *does* hold, so the real task is **skipping holes
cheaply** instead of encoding direction in them:

* **Per-column "valid" index.** For each column `c`, a rank/select bit-vector marking
  entries valid at `c` (≈ `(k−m+1)` bits per entry total) lets the search jump to the
  nearest valid entry in O(1) and binary-search only valid entries. Costs memory; this is
  the most direct replacement for `find_closest_valid_skmer`.
* **Minimizer-range index** (a *different* lever): since all k-mers of a query super-k-mer
  share one minimizer and the list is grouped by minimizer, a `minimizer → [start,end)`
  table confines the whole search to one block — cache locality + smaller ranges. (This was
  the originally-considered idea before the sentinel direction.)
* Or accept the scan (it is bounded by the per-minimizer block in practice) and focus on
  cache locality / parallelism.

## Reproduce

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release -j
# honest unit tests (incl. the at-scale invariant and the documented limitation):
(cd build-release && ./tests/sklib-tests --gtest_filter='SentinelSubstrate.*')
# size sweep + 2 Mb timing/correctness + fill-then-resort experiment:
(cd build-release && SKLIB_BENCH=1 ./tests/sklib-tests --gtest_filter='SentinelBench.*')
```
