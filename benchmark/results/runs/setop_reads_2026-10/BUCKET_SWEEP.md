# Bucket-count sweep and the stale-expandable construction bug — 2026-10-06

Follow-up of the campaign (`EXPERIMENT.md`). At t=8, sklib's union throughput fell with the data
size (96 Mk-mer/s on chr1, 53 on ocean k=31), CPU-bound, while KMC stayed flat. The hypothesis was
that the default 4096 buckets make each bucket too large for the caches at that scale.

## Protocol

Same as the campaign: cold cache, t=8 pinned (`0,1,3,6,8,10,12,17`), t=1 unpinned, REPS=3.
- Each bucket count is a separate sklib "tool" (`sklib@b<N>`, built with `--buckets N`), so the
  variants alternate run by run.
- KMC is not re-run.
- Driver: `drivers/bucket_sweep.sh`. Data: `data/bucket_sweep_runs.csv` and
  `data/bucket_sweep_construct_runs.csv`.
- Runs: 2026-10-06, 06:01–07:25 and 18:20–22:21, with a pause at 20:17–20:45 for the diagnosis
  below.
- The 07:18 chr1 rows were measured while the laptop was in use. They are set aside in
  `data/bucket_sweep_runs_chr1_discarded.csv` and chr1 was re-measured at the end.

## Results (median wall seconds; % vs the 4096 default)

| Pair, op | 4096 | 16384 | 65536 | 262144 | other |
|---|---|---|---|---|---|
| chr1 k31 union t=8 (50 k k-mers/bucket at 4096) | **2.8** | 3.0 (+8 %) | 6.8* | — | 1024: 3.3, 256: 3.8 |
| gut k31 union t=8 (151 k) | 12.3 | 10.4 (−15 %) | **10.2 (−17 %)** | 10.8 | 1024: 13.9 |
| gut k31 joint t=8 | 18.6 | 15.4 | **15.1 (−19 %)** | 16.2 | 1024: 19.8 |
| HG002 k31 union t=8 (737 k) | 55.0 | **48.8 (−11 %)** | 56.4* | 61.3* | |
| HG002 k31 joint t=8 | 80.6 | **73.5 (−9 %)** | 76.0* | 111.7* | |
| ocean k31 union t=8 (984 k) | 113.5 | 96.2 | 87.3 (−23 %) | **85.2 (−25 %)** | |
| ocean k31 joint t=8 | 148.9 | 124.2 | **117.0 (−21 %)** | 124.6 | |
| ocean k31 union t=1 | 399.1 | 375.0 | 334.1 | **319.9 (−20 %)** | |
| ocean k63 union t=8 (825 k) | 97.3 | 77.1 | **69.7 (−28 %)** | wrong result† | |
| ocean k63 joint t=8 | 148.4 | 110.9 | **101.8 (−31 %)** | wrong result† | |
| `--sizes` (all pairs) | ≈ flat | ≈ flat | +10–25 % | +20–65 % | |

\* Noisy cells: %CPU fell to 90–580 % and some reps are outliers, likely because the laptop was in
use or because small cold-cache reads stalled.
† The 262144 cells at k=63 have a wrong result. See the bug below.

**Construction (t=8).** 65536 buckets is the fastest everywhere: −15 to −20 % against 4096, with
1.7 GB of RAM. 262144 is slower and needs 6.8 GB.

**Reading.**
- The optimum is around **10–60 k k-mers per bucket per input**.
- Buckets above about 200 k k-mers cost 15–30 % on materialized operations. That is the 4096
  default on every large set.
- Buckets below about 5 k k-mers cost again: a fixed per-bucket overhead plus small cold-cache
  reads (chr1 at 65536).
- Counting only (`--sizes`) does not gain. The whole effect is in the materialization: collecting
  the kept k-mers, re-chaining, writing.
- No single `--buckets` value fits every size: the best value goes from 4096 (chr1) to
  65536–262144 (ocean).

## Bug found: closed virtual super-k-mers left "expandable" (fixed in the working tree)

**Symptom.** At ocean k=63 with `--buckets 262144`, |A| and |B| matched KMC, but the merge missed 6
shared k-mers out of 844 M: |A∩B| was 6 too low, each difference 6 too high, and the union had 6
duplicates. The other bucket counts were exact.

**Localisation.**
- The bug is in the construction, not in the merge. One bucket of `ocean_DCM` violated the
  per-column sorted invariant (record order = k-mer order at every column), which both
  `merge_columns` and `search_kmers_in_span` rely on.
- Rebuilt with `SKLIB_NO_DECOUPLE=1`, the same bucket compacted at b=16 reproduced the 2 violations
  exactly; at b=18, 3 other violations appeared in another bucket. The cause therefore depends on
  the exact content of the compacted bucket, not on the decouple.
- The 5,751 raw super-k-mers of that bucket, extracted from the reads, reproduce the violations with
  greedy chaining (the construction default since `f6543af`) and not with colinear chaining.
- Delta-debugging reduced them to **18 super-k-mers**.

**Cause** (`merge_LList_column`, `lib/include/algorithms/VirtualSkmer.hpp`).
1. One raw super-k-mer gives k-mers to several columns under the same enumeration id.
2. When a virtual super-k-mer stopped growing, the two tail loops passed it through without setting
   `expandable = false`, unlike CASE C and CASE D.
3. A later column's overlap whose left k-mer had that id matched the stale record first, because
   `is_left` compares `last_id` only. The stale record was then extended past a gap of columns.
4. The k-mers were still right, since they belong to the same raw super-k-mer, but the record kept a
   list position that was wrong for its new columns.

**Effect.**
- Set ops miss shared k-mers.
- In the reduced case, 46 input k-mers are not found by a query (false negatives).
- It is very rare: 2 violations in 3.15 G adjacent pairs. A full scan found none in any ocean k=63
  list at 4096, 16384 or 65536 buckets.
- The set-op re-compaction is not exposed: its enumeration is single k-mers, one column each, so an
  id is never reused across columns.

**Fix.** Close every record passed through in the tail loops. 10 lines, uncommitted.
- Regression test: `tests/km/stale_expandable_vskmer.cpp`, built on the 18 super-k-mers. It
  fails before the fix (columns 25/27 out of order, 46 k-mers not found) and passes after.
- Full suite: 217/217 (Debug, clang-18).
- Ocean k=63 at 262144, rebuilt with the fix: 0 violations, and the set-op cardinalities equal the
  KMC-validated ones.
- Six campaign indexes rebuilt with the fix (default 4096) are byte-identical to the originals:
  chr1 k31, gut_A k31/k63, hg002_A k63, ocean_DCM k63, ocean_SRF k31. The campaign results stand.

Debug tools (git-ignored): `benchmark/results/latest/setop_reads/tmp/debug/`. They are
`bucket_diff`, `kmer_diff2`, `invariant`, `extract`, `ddmin` and `trace`.
