# Overnight summary — absent-slot sentinel (for Yoann)

**Bottom line: negative result.** Your "write an order-sentinel into the unused bits so the
dichotomy avoids the linear scan" idea was implemented and measured end-to-end. It is *safe*
but does **not** enable scan-free queries beyond trivially small inputs. Full technical
analysis + measurements: [`docs/sentinel_substrate.md`](docs/sentinel_substrate.md).

Both branches are derived from `dev`, pushed, and all tests are green.

## Branches / "versions"

| Branch | What it is | State | Recommendation |
|--------|-----------|-------|----------------|
| `feat/sentinel-construction` (**v1**) | construction-only: `fill_absent_sentinel()` clears absent flank slots to 0; honest tests + docs + investigation note | green; query results **byte-identical** to pristine `dev` (verified on a 2 Mb genome) | **revert candidate** — safe but no demonstrated benefit |
| `feat/sentinel-query-e2e` (**v3**) | v1 + `query_skmer_substrate()` (scan-free hole-aware binary search, PoC) + `SKLIB_BENCH` micro-benchmark + `docs/sentinel_substrate.md` | green (unit tests honest); the PoC query is ~2–3× faster **but incorrect at scale** | **don't merge the query**; keep the bench + doc as reference for the query phase |

(The optional "v2 precision-fill" branch from the plan was made moot: no single fill makes
the substrate navigable, so a smarter fill policy can't help — see below.)

## What was measured (k=21, m=11, random genome)

* Valid-entry per-column monotonicity (the invariant the current query relies on): **0
  violations at every scale** — solid.
* Hole-inclusive monotonicity (what the sentinel needed): breaks with size —
  0 (1 kb) → 1 (10 kb) → 555 (100 kb) → 58 037 (1 Mb) → 213 283 (2 Mb).
* PoC hole-aware query on 2 Mb: **~5.3 % false negatives** on present k-mers (0 on absent),
  while ~2–3× faster than the scan-based query. Speed without correctness.
* Fill **+ re-sort by `m_pair`** does not rescue it (still 59 326 violations) **and breaks
  the existing query** (16 589 mismatches) — the construction order is essential and is
  *not* the `m_pair` order.

## Why it fails (root cause)

`kmer_masks[c]` is **not** a high-order prefix of `m_pair`; in the minimizer-centred
interleaved layout it is a *scattered* mask. So sorting by `m_pair` does not order the
per-column projections, and a "hole" (entry lacking a k-mer at column `c`) is placed by the
columns where it *is* valid — at column `c` its **high-order** content is already out of
order, and no absent-bit fill can repair high-order bits. (My initial "interleaved magic:
sort by m_pair ⇒ every column sorted" claim was simply wrong; corrected in the commits.)

## Recommended next step (query-optimization phase)

The valid-entry monotonicity *does* hold, so the real task is **skipping holes cheaply**,
not encoding direction in them:
* **per-column rank/select index** (≈ `k−m+1` bits/entry) → jump to the nearest valid
  entry in O(1), binary-search valid entries only — the direct replacement for
  `find_closest_valid_skmer`; or
* **minimizer→[start,end) index** (a different lever) → confine each query super-k-mer's
  search to one minimizer block (cache locality + smaller ranges).

## Verify

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release -j
(cd build-release && ./tests/sklib-tests --gtest_filter='SentinelSubstrate.*')      # honest unit tests
(cd build-release && SKLIB_BENCH=1 ./tests/sklib-tests --gtest_filter='SentinelBench.*')  # sweep + 2 Mb timing/correctness (v3)
```

## Notes
* Format magic was NOT bumped: filled and unfilled lists stay mutually queryable by the
  current query (the fill is invisible to it).
* `build-release/` is untracked (not in `.gitignore`; only `build/`,`build-debug/`,
  `build-profile/` are). Local only.
* `dev` already exists on `origin` and equals `origin/main` (both at the φ-order 0.2.0
  milestone); branches are based on it.
