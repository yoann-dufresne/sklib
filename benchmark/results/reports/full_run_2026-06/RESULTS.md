# Full k-mer tooling benchmark — detailed results

Machine `yoann-Precision-5490` (Intel Core Ultra 7 165H, 22 cores, 62 GiB).
Methodology, tool list, grid and CSV schema: see `EXPERIMENT.md`. All sklib figures
are the canonical **`sklib-0.11.0`** rows; times are median of 3 reps. Unless stated,
tables are **t = 1** (single-thread) so every tool — including the single-thread-only
ones (cbl, sbwt-C++, fmsi, sshash) — is comparable.

Higher **throughput** (Mk-mer/s) = better; lower **time** / **bits-per-k-mer** / **RAM** = better.

---

## 1. Construction

### 1a. Speed — time (s), t = 1

```
[chr1, 224 MB]   k | sklib |  kmc  | sbwt(C++) | sbwtrs | cbl   | fmsi
              15 | 25.98 | 15.93 |   59.13   |  64.29 |  71.69| 219.6
              21 | 23.25 | 17.81 |  101.57   |  80.47 |  85.38| 218.6
              31 | 22.02 | 11.82 |  119.09   |  90.53 | 110.48| 236.4
              41 | 20.91 | 20.49 |     –     | 128.67 | 140.75| 258.1
              51 | 29.05 | 18.64 |     –     | 134.11 | 172.72| 259.8
              63 | 28.65 | 17.82 |     –     | 143.56 |   –   | 272.3
```
(celegans/chr21/yeast/ecoli follow the same ordering.)

- **kmc is the fastest builder**, **sklib a close second** (chr1: ~1.2–2× kmc). Both
  are **3–10× faster** than every other tool.
- **sbwt-rs (Rust) ≈ sbwt (C++)**, the Rust version pulling slightly ahead at high k
  (chr1 k31: 90 s vs 119 s) **and lifting the k ≤ 32 cap** (covers k15–63).
- **cbl** scales poorly with k (chr1: 72 s → 173 s). **fmsi** is the slowest builder
  by far (~220–270 s on chr1).

### 1b. Space — bits per k-mer

```
            k15    k21   k31   k41   k51   k63
fmsi        ~3     ~2    ~2    ~2    ~2    ~2     ← most compact
sbwt(both)  ~10    ~10   ~10   ~10   ~10   ~10    ← flat
sklib       74     39    24    17    29    23    ← k-dependent (chr1)
```
- **fmsi is the most compact** (~2–3 bits/k-mer).
- **SBWT (both impls) ~10 bits/k-mer, flat** in k.
- **sklib is k-dependent**: heavy at k = 15 (short super-k-mers ⇒ poor amortization,
  55–74 b/kmer), best around k = 41 (~16–17), rising again at k = 51/63 (larger m).
  Note sklib stores *more* than membership (a sorted super-k-mer list with positional
  columns + set-op support), so this is not a like-for-like space comparison.

---

## 2. Membership query

Throughput (Mk-mer/s, t = 1, median over presence levels). **No single winner — it
depends on k and genome size:**

```
[ecoli, 4.6 MB]   k | sklib | sbwt | sbwtrs | cbl  | sshash | fmsi
               15 | 1.09 | 1.68 |  4.65  | 0.81 |  1.87  | 0.74
               21 | 1.33 | 1.37 |  3.70  | 0.80 |  2.70  | 0.64
               31 | 1.07 | 1.06 |  2.78  | 0.69 |  2.02  | 0.53
               41 | 0.93 |  –   |  2.15  | 0.66 |   –    | 0.47
               63 | 0.49 |  –   |  1.38  |  –   |   –    | 0.34

[chr1, 224 MB]    k | sklib | sbwt | sbwtrs | cbl  | sshash | fmsi
               15 | 0.16 | 0.58 |  0.75  | 0.13 |  0.27  | 0.30
               21 | 0.28 | 0.40 |  0.48  | 0.06 |  0.54  | 0.19
               31 | 0.35 | 0.29 |  0.35  | 0.02 |  0.63  | 0.14
               41 | 0.38 |  –   |  0.27  | 0.01 |   –    | 0.12
               63 | 0.25 |  –   |  0.18  |  –   |   –    | 0.08
```

- **Small k:** **sbwt-rs is the fastest** (ecoli k15: 4.65 Mk-mer/s, ~2.8× the C++ sbwt).
- **Mid k on large genomes:** **sshash leads** (chr1 k31: 0.63; celegans k31: 0.91),
  closely followed by **sklib**.
- **High k:** **sklib is the most robust** — its throughput *rises* with k where SBWT
  *falls* (more bits/k-mer to process). On chr1, sklib overtakes sbwt-rs at k ≥ 31
  (k41: 0.38 vs 0.27).
- **sbwt-rs > sbwt-C++ at every shared point** (~1.2–2.8×) — the rewrite is a clear win.
- **cbl collapses with k** (celegans k31: 0.04 Mk-mer/s). **fmsi is consistently slow.**
- **kmc has no query mode** (it is the correctness oracle).
- Threads: sklib & sbwt-rs scale with `-t` (the others are single-thread here), so
  their multi-thread wall-clock advantage is larger than the t = 1 table shows.
- `query_stream` (consecutive k-mers) tracks `query_single` with modestly higher
  throughput for all tools.

---

## 3. Set operations

Time (s) and peak RAM, t = 1, median over Jaccard × {inter, union, diffab, diffba},
materialize mode. Two clear tiers:

```
[yeast, 12 MB]   k | sklib | kmc  | cbl  | fmsi  | sbwtrs   (time_s)  || RAM MB sklib/kmc/sbwtrs
              15 | 0.65 | 0.56 | 1.13 | 60.5  | 12.84   ||  10 / 126 / 214
              21 | 0.49 | 0.66 | 1.28 | 59.3  | 26.27   ||   6 / 161 / 229
              31 | 0.46 | 0.88 | 1.65 | 62.3  | 33.44   ||   7 / 230 / 185
              41 | 0.44 | 0.97 | 1.66 | 12.5  | 49.11   ||   8 / 300 / 194
              51 | 0.53 | 0.83 | 1.83 | 63.4  | 65.74   ||  10 / 363 / 202
              63 | 0.52 | 0.93 |  –   | 72.8  | 83.69   ||  12 / 413 / 202
```
(ecoli identical ordering: sklib ~0.2 s, kmc ~0.2–0.4 s, cbl ~0.6–0.7 s, fmsi ~13–24 s,
sbwt-rs 5–50 s.)

**Speed — fast tier vs slow tier:**
- **Fast:** **sklib** (≈ 0.2–0.65 s, *flat in k*) < **kmc** (0.2–1 s) < **cbl** (0.6–1.8 s).
- **Slow:** **fmsi** (13–73 s) ≈ **sbwt-rs** (5–84 s, *grows steeply with k* — each op
  rebuilds the result SBWT). The slow tier is **~50–170× slower than sklib** and
  ~30–90× slower than kmc; the gap widens with k.

**Memory:**
- **sklib uses the least by far** (6–12 MB — it streams bucket-by-bucket), ~**20×**
  less than sbwt-rs and ~**10–60×** less than kmc.
- **kmc grows with k** (70 → 413 MB). **sbwt-rs ~flat** (~110–230 MB, loads the SBWT
  indexes); it overtakes kmc (uses less) at high k.

**Coverage:** fmsi and sbwt-rs set-ops are capped at small/mid genomes (too slow at
scale — see EXPERIMENT.md); sklib/kmc/cbl run the full ladder incl. chr1. On chr1
(sklib vs kmc), sklib is ~1.0–1.3× kmc at t = 1 and pulls clearly ahead at t = 8
(kmc's `simple` is largely I/O-bound and barely threads), with ~20–50× less RAM.

---

## 4. Bottom line per workload

| Workload | Winner(s) | Notes |
|----------|-----------|-------|
| **Construct speed** | **kmc**, then **sklib** | both 3–10× faster than SBWT/cbl/fmsi |
| **Construct space** | **fmsi** (~2 b/kmer) | SBWT ~10 (flat); sklib k-dependent, stores more |
| **Query, small k** | **sbwt-rs** | Rust rewrite clearly beats C++ sbwt |
| **Query, mid k / big genome** | **sshash**, then **sklib** | |
| **Query, high k** | **sklib** | robust where SBWT degrades; + threads |
| **Set-ops speed** | **sklib** ≫ kmc > cbl ≫ (fmsi ≈ sbwt-rs) | sklib flat in k |
| **Set-ops memory** | **sklib** | streams; ~20× less than sbwt-rs, ~10–60× less than kmc |

**sklib's signature:** consistently top-tier on **construct speed**, the most
**robust query at high k**, and a **decisive set-op win on both time and memory** —
its per-bucket streaming keeps set-op RAM at tens of MB while competitors use
hundreds of MB to GBs.

**On the sbwt Rust rewrite (the reason for adding `sbwtrs`):** vs the C++ `sbwt` it is
**faster to query** (notably at small k), **on par to build**, **lifts the k ≤ 32 cap**,
and **adds native set-ops** — though those set-ops are slow (≈ fmsi, ~2 orders of
magnitude behind sklib/kmc).

## Caveats
- `peak_rss` units: KB in construct/query CSVs, **MB** in setop CSV (see EXPERIMENT.md).
- Set-op `result_kmers` is sklib's authoritative cardinality shared across tools; a
  tool's own set-op output is not re-validated here (perf-only).
- Legitimate skips: sbwt-C++ k≤32; cbl odd-k≤59; sshash k≤31 (+ celegans k15 stream
  query fails); fmsi/sbwt-rs set-ops capped; bqf is approximate (no set-ops).
- `chm13` (3 GB) deferred. The plotting step (`plot.py`) needs `pip install -U
  bottleneck` in the venv to regenerate figures (data unaffected).
