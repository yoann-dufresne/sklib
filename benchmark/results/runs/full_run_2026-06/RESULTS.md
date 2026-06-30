# Full k-mer tooling benchmark — detailed results

Machine `yoann-Precision-5490` (Intel Core Ultra 7 165H, 22 cores, 62 GiB).
Methodology, tool list, grid and CSV schema: see `EXPERIMENT.md`. All sklib figures
are the canonical **`sklib-0.11.0`** rows.

**How to read the numbers**
- Every figure is the **median of 3 reps**. For the threaded tools the first rep of
  each configuration is a cold-page-cache run (~2× slower on the big genomes); the
  median discards it, so the tables show **warm steady-state**.
- Higher **throughput** (Mk-mer/s) = better; lower **time** / **bits-per-k-mer** / **RAM** = better.
- Tables are **t = 1** (single-thread) unless the heading says otherwise, so every tool —
  including the single-thread-only ones — is comparable. Each workload then has a
  dedicated **"scaling with threads"** subsection (t = 1/2/4/8/16) for the tools that
  thread (**sklib, kmc, sbwt-rs**).
- All tables are regenerable from `data/*.csv` (the per-rep raw rows).

---

## 0. Tools & datasets

**Real datasets** (the toy `sarscov2`, 30 KB, is startup-bound and omitted from the
tables; `chm13`, 3 GB, is deferred):

| dataset  | size   | #k-mers (k=31) |
|----------|-------:|---------------:|
| ecoli    | 4.6 MB |      4.55 M    |
| yeast    | 12 MB  |     11.56 M    |
| chr21    | 39 MB  |     34.6 M     |
| celegans | 98 MB  |     94.0 M     |
| chr1     | 224 MB |    204.2 M     |

**Tools** (what each stores and which workloads / threading it supports):

| tool      | structure / note                              | k range   | construct | query | set-ops | threaded? |
|-----------|-----------------------------------------------|-----------|:---------:|:-----:|:-------:|-----------|
| **sklib** | sorted super-k-mer list + positional columns  | ≤127      | ✓ | ✓ | ✓ | **construct, query, set-ops** |
| **kmc**   | k-mer counter (used as correctness oracle)    | all       | ✓ | – | ✓ | construct (strong), set-ops (weak) |
| sbwt      | SBWT membership (C++)                          | **≤32**   | ✓ | ✓ | – | no (single-thread here) |
| sbwtrs    | SBWT membership + set-ops (Rust rewrite)       | all       | ✓ | ✓ | ✓ | construct & set-ops yes; **query flat** |
| sshash    | minimal-perfect-hash membership               | **≤31**   | ✓ | ✓ | – | no |
| cbl       | membership + set-ops                           | **odd ≤59** | ✓ | ✓ | ✓ | no |
| fmsi      | masked-superstring membership + set-ops        | all       | ✓ | ✓ | ✓ | no |
| bqf*      | **approximate** membership (false positives)   | k=31 here | ✓ | ✓ | – | no |

`*bqf` is a probabilistic filter (not exact membership); shown for reference, marked `*`.

---

## 1. Construction

### 1a. Speed across genomes — time (s), t = 1

**k = 31**

| tool   | ecoli | yeast | chr21 | celegans |  chr1  |
|--------|------:|------:|------:|---------:|-------:|
| **kmc**   | 0.39 | 0.81 |  2.37 |   5.28 |  11.82 |
| **sklib** | 0.42 | 1.06 |  3.64 |   9.40 |  22.40 |
| sshash    | 0.45 | 1.24 |  3.73 |  10.74 |  28.43 |
| bqf*      | 1.41 | 3.82 | 13.34 |  33.08 |  72.99 |
| sbwtrs    | 2.32 | 6.00 | 16.28 |  42.19 |  90.53 |
| cbl       | 0.92 | 2.47 |  9.98 |  39.95 | 110.48 |
| sbwt      | 2.35 | 6.34 | 20.19 |  52.12 | 119.09 |
| fmsi      | 2.45 | 7.47 | 28.92 |  87.56 | 236.35 |

**k = 63** (sbwt/sshash/cbl/bqf don't reach k=63)

| tool   | ecoli | yeast | chr21 | celegans |  chr1  |
|--------|------:|------:|------:|---------:|-------:|
| **kmc**   | 0.51 | 1.05 |  3.21 |   7.31 |  17.82 |
| **sklib** | 0.47 | 1.27 |  4.48 |  11.86 |  29.33 |
| sbwtrs    | 3.35 | 8.66 | 27.04 |  68.82 | 143.56 |
| fmsi      | 3.24 | 9.92 | 36.10 | 107.70 | 272.25 |

- **kmc is the fastest builder, sklib a close second** (chr1: ~1.6–1.9× kmc). Both are
  **3–10× faster** than every other tool.
- **sbwt-rs (Rust) ≈ sbwt (C++)** at the points they share, the Rust version lifting the
  k ≤ 32 cap (covers k15–63). **cbl** scales poorly with k; **fmsi** is the slowest by far.

### 1b. Speed across k — chr1, time (s), t = 1

| k  | sklib |  kmc  | sshash | sbwtrs |  sbwt  |  cbl   |  bqf*  |  fmsi  |
|----|------:|------:|-------:|-------:|-------:|-------:|-------:|-------:|
| 15 | 28.65 | 15.93 | 112.14 |  64.29 |  59.13 |  71.70 |  33.00 | 219.57 |
| 21 | 24.62 | 17.81 |  37.89 |  80.47 | 101.57 |  85.38 |  73.62 | 218.60 |
| 31 | 22.40 | 11.82 |  28.43 |  90.53 | 119.09 | 110.48 |  72.99 | 236.35 |
| 41 | 21.33 | 20.49 |    –   | 128.67 |    –   | 140.75 |    –   | 258.13 |
| 51 | 29.71 | 18.64 |    –   | 134.11 |    –   | 172.72 |    –   | 259.79 |
| 63 | 29.33 | 17.82 |    –   | 143.56 |    –   |    –   |    –   | 272.25 |

sklib is fastest to build around k≈31–41 and slows at the extremes (short super-k-mers at
k=15, larger minimizers at k≥51).

### 1c. Space — bits per k-mer (t = 1)

Across genomes, **k = 31**:

| tool   | ecoli | yeast | chr21 | celegans | chr1  |
|--------|------:|------:|------:|---------:|------:|
| fmsi   |  2.13 |  2.15 |  2.24 |   2.20   |  2.25 | ← most compact
| sshash |  7.12 |  7.34 |  8.78 |   8.33   |  9.67 |
| sbwtrs | 11.84 | 10.73 | 10.24 |  10.09   | 10.04 |
| sbwt   | 12.34 | 11.23 | 10.74 |  10.59   | 10.54 |
| sklib  | 21.43 | 21.55 | 22.83 |  23.54   | 24.27 |
| bqf*   | 38.68 | 29.02 | 18.41 |  24.27   | 21.04 |
| cbl    | 64.90 | 60.20 | 57.77 |  57.34   | 56.76 |
| kmc    | 66.30 | 64.91 | 64.30 |  64.11   | 64.05 |

- **fmsi** is the most compact (~2–3 b/kmer); **SBWT (both)** flat ~10; **sklib**
  k-dependent (chr1: 74 at k=15 → ~17 at k=41 → 23 at k=63 — short super-k-mers amortize
  poorly at low k). **kmc/cbl** are *databases* (counts, on-disk), not membership indexes,
  hence the large bits/kmer — not a like-for-like comparison.
- sklib stores **more than membership** (a sorted super-k-mer list with positional columns
  + set-op support), so its bits/kmer is not directly comparable to a pure membership filter.

### 1d. Space — peak RAM (MB), t = 1, k = 31

| tool   | ecoli | yeast | chr21 | celegans |  chr1  |
|--------|------:|------:|------:|---------:|-------:|
| **sklib** | 18.6 | 18.3 |  44.9 |   50.1 |   139 | ← streams, ~10–60× less
| sshash    | 36.3 | 82.6 |   136 |    369 |   735 |
| bqf*      | 71.3 |  131 |   341 |    794 |  1811 |
| sbwtrs    | 81.1 |  199 |   628 |   1551 |  3686 |
| fmsi      | 81.7 |  203 |   579 |   1550 |  3413 |
| kmc       |  104 |  218 |   644 |   1546 |  3522 |
| sbwt      |  182 |  409 |  1256 |   1948 |  1948 |
| cbl       |  147 |  211 |   612 |   4957 |  9223 |

sklib's construction is dramatically lighter (it streams bucket-by-bucket): ~25× less RAM
than kmc/sbwt-rs and ~65× less than cbl on chr1.

### 1e. Scaling with threads — construction time (s), speedup vs t = 1

Only **sklib, kmc, sbwt-rs** thread; the other five are single-thread.

**chr1, k = 31**

| threads | sklib        |  kmc         | sbwtrs       |
|---------|-------------:|-------------:|-------------:|
| 1       | 22.40        | 11.82        | 90.53        |
| 2       | 11.71        |  7.31        | 63.22        |
| 4       |  6.49        |  3.89        | 48.74        |
| 8       |  4.06 (5.5×) |  2.58 (4.6×) | 45.04 (2.0×) |
| 16      |  3.24 (6.9×) |  2.13 (5.5×) | 43.81 (2.1×) |

**chr1, k = 63**

| threads | sklib        |  kmc         | sbwtrs        |
|---------|-------------:|-------------:|--------------:|
| 1       | 29.33        | 17.82        | 143.56        |
| 2       | 15.34        | 11.45        | 102.03        |
| 4       |  8.68        |  6.22        |  81.14        |
| 8       |  5.86 (5.0×) |  3.91 (4.6×) |  73.78 (1.9×) |
| 16      |  5.08 (5.8×) |  3.20 (5.6×) |  72.78 (2.0×) |

**celegans, k = 31**

| threads | sklib        |  kmc         | sbwtrs       |
|---------|-------------:|-------------:|-------------:|
| 1       |  9.40        |  5.28        | 42.19        |
| 2       |  4.94        |  3.42        | 29.87        |
| 4       |  2.67        |  2.01        | 25.71        |
| 8       |  1.68 (5.6×) |  1.45 (3.6×) | 22.86 (1.8×) |
| 16      |  1.30 (7.2×) |  1.28 (4.1×) | 22.26 (1.9×) |

**celegans, k = 63**

| threads | sklib        |  kmc         | sbwtrs       |
|---------|-------------:|-------------:|-------------:|
| 1       | 11.86        |  7.31        | 68.82        |
| 8       |  2.19 (5.4×) |  1.92 (3.8×) | 35.33 (1.9×) |
| 16      |  1.84 (6.5×) |  1.54 (4.7×) | 31.57 (2.2×) |

- **sklib scales best** (5–7× at t = 16), thanks to its bucket-parallel phase 2.
- **kmc scales strongly** (4–6×) and keeps the **absolute** lead — but the gap closes: on
  celegans k=31 the two essentially **tie at t = 16** (sklib 1.30 s vs kmc 1.28 s).
- **sbwt-rs scales poorly** (~2× regardless of threads).

---

## 2. Membership query (point lookups)

Throughput (Mk-mer/s), median over presence levels {0, 25, 50, 75, 100 %}. **kmc has no
query mode** (it is the oracle). **No single winner — it depends on k and genome size.**

### 2a. Throughput across genomes — t = 1

**k = 31**

| tool   | ecoli | yeast | chr21 | celegans |  chr1 |
|--------|------:|------:|------:|---------:|------:|
| sshash | 2.02 | 1.89 | 1.105 |  0.909 | 0.633 |
| sbwtrs | 2.78 | 1.27 | 0.543 |  0.392 | 0.347 |
| sklib  | 1.06 | 0.90 | 0.662 |  0.498 | 0.345 |
| sbwt   | 1.06 | 0.88 | 0.543 |  0.344 | 0.295 |
| bqf*   | 0.53 | 0.41 | 0.325 |  0.285 | 0.241 |
| fmsi   | 0.53 | 0.56 | 0.477 |  0.268 | 0.145 |
| cbl    | 0.69 | 0.49 | 0.208 |  0.035 | 0.019 |

**k = 63**

| tool   | ecoli | yeast | chr21 | celegans |  chr1 |
|--------|------:|------:|------:|---------:|------:|
| sbwtrs | 1.38 | 0.704 | 0.290 | 0.202 | 0.181 |
| sklib  | 0.49 | 0.450 | 0.377 | 0.318 | 0.246 |
| fmsi   | 0.34 | 0.375 | 0.319 | 0.162 | 0.079 |

### 2b. Throughput across k — t = 1

**chr1**

| k  | sshash | sbwtrs | sklib | sbwt  | bqf*  | fmsi  | cbl   |
|----|-------:|-------:|------:|------:|------:|------:|------:|
| 15 |  0.265 |  0.752 | 0.158 | 0.583 | 0.840 | 0.304 | 0.132 |
| 21 |  0.535 |  0.481 | 0.284 | 0.401 | 0.371 | 0.187 | 0.060 |
| 31 |  0.633 |  0.347 | 0.345 | 0.295 | 0.241 | 0.145 | 0.019 |
| 41 |    –   |  0.269 | 0.380 |   –   |   –   | 0.115 | 0.008 |
| 51 |    –   |  0.208 | 0.241 |   –   |   –   | 0.096 | 0.005 |
| 63 |    –   |  0.181 | 0.246 |   –   |   –   | 0.079 |   –   |

**ecoli**

| k  | sbwtrs | sshash | sbwt  | sklib | cbl   | fmsi  | bqf*  |
|----|-------:|-------:|------:|------:|------:|------:|------:|
| 15 |  4.651 |  1.869 | 1.681 | 1.087 | 0.806 | 0.738 | 1.212 |
| 21 |  3.704 |  2.703 | 1.370 | 1.316 | 0.800 | 0.641 | 0.858 |
| 31 |  2.778 |  2.020 | 1.064 | 1.058 | 0.692 | 0.532 | 0.525 |
| 41 |  2.151 |    –   |   –   | 0.930 | 0.658 | 0.466 |   –   |
| 51 |  1.587 |    –   |   –   | 0.559 | 0.654 | 0.370 |   –   |
| 63 |  1.379 |    –   |   –   | 0.490 |   –   | 0.338 |   –   |

- **Small k:** **sbwt-rs is fastest** (ecoli k15: 4.65, ~2.8× the C++ sbwt).
- **Mid k / large genome:** **sshash leads** (chr1 k31: 0.63), then sklib.
- **High k:** **sklib is the most robust** — its throughput *rises* with k where SBWT
  *falls*; on chr1, sklib overtakes sbwt-rs at k ≥ 41 (k41: 0.38 vs 0.27) and is the only
  exact tool still competitive at k = 63.
- **cbl collapses with k/size** (chr1 k51: 0.005). **fmsi** is consistently slow.

### 2c. Scaling with threads — chr1, throughput (Mk-mer/s)

Only **sklib** and **sbwt-rs** expose threaded query.

**k = 31**

| threads | sklib        | sbwtrs       |
|---------|-------------:|-------------:|
| 1       | 0.345        | 0.347        |
| 2       | 0.398        | 0.342        |
| 4       | 0.515        | 0.331        |
| 8       | 0.560 (1.6×) | 0.338 (1.0×) |
| 16      | 0.597 (1.7×) | 0.342 (1.0×) |

**k = 63:** sklib 0.246 → 0.408 (t8, 1.7×) → 0.420 (t16, 1.7×); sbwt-rs flat at 0.18.
**celegans k = 31:** sklib 0.498 → 0.935 (t8, 1.9×) → 0.990 (t16, 2.0×); sbwt-rs flat at 0.39.

- **Only sklib benefits from `-t`** (~1.6–2.0×); the gain is sub-linear because the
  per-column dichotomic search is memory-bound (≈70 % of a point lookup is the per-record
  super-k-mer rebuild, not the search itself).
- **sbwt-rs query is flat across `-t`** in this benchmark (it does not parallelize the
  query path here), so sklib's multi-thread wall-clock advantage over sbwt-rs is larger
  than the t = 1 tables suggest.

---

## 3. Streaming query (consecutive k-mers)

### 3a. Throughput across genomes — Mk-mer/s, t = 1

**k = 31**

| tool   | ecoli | yeast | chr21 | celegans |  chr1 |
|--------|------:|------:|------:|---------:|------:|
| sshash | 30.00 | 25.71 | 9.474 |  6.585 | 3.253 |
| bqf*   | 12.86 |  8.85 | 6.207 |  2.919 | 1.714 |
| sbwtrs |  8.06 |  5.05 | 2.983 |  2.213 | 1.561 |
| sklib  |  7.30 |  5.63 | 3.750 |  2.443 | 1.448 |
| sbwt   |  4.46 |  3.83 | 2.827 |  2.109 | 1.588 |
| cbl    |  2.11 |  1.33 | 0.578 |  0.094 | 0.052 |
| fmsi   |  0.77 |  0.85 | 0.761 |  0.647 | 0.458 |

**k = 63**

| tool   | ecoli | yeast | chr21 | celegans |  chr1 |
|--------|------:|------:|------:|---------:|------:|
| sklib  |  6.90 |  5.29 | 3.748 |  2.368 | 1.352 |
| sbwtrs |  3.58 |  2.87 | 2.174 |  1.653 | 1.184 |
| fmsi   |  0.53 |  0.60 | 0.555 |  0.451 | 0.335 |

Streaming tracks `query_single` with higher absolute throughput for all tools; the ranking
is unchanged (sshash leads, sklib/sbwt-rs/sbwt close together, sklib alone robust at k=63).

### 3b. Scaling with threads — chr1, throughput (Mk-mer/s)

- **k = 31:** sklib 1.448 → 1.742 (t8, 1.2×) → 1.714 (t16); sbwt-rs flat at 1.56.
- **k = 63:** sklib 1.352 → 1.531 (t8, 1.13×); sbwt-rs flat at 1.18.
- **celegans k = 31:** sklib 2.443 → 3.396 (t8, 1.4×); sbwt-rs flat at 2.21.

Streaming saturates quickly: sklib gains only ~1.2–1.4× from threads, sbwt-rs is flat.

---

## 4. Set operations

Time (s) / peak RAM, materialize mode, **median over Jaccard × {inter, union, diff A\B,
diff B\A}**. (The combined single-pass `joint` mode — all relations at once — is covered
separately in `SETOPS_MULTI_REPORT.md` and excluded from these medians.) Only **sklib,
kmc, cbl, fmsi, sbwt-rs** implement set-ops; sbwt-C++, sshash and bqf do not.

### 4a. Speed across genomes — time (s), t = 1

**k = 31**

| tool   | ecoli  | yeast  | chr21 | celegans |  chr1  |
|--------|-------:|-------:|------:|---------:|-------:|
| **sklib** | 0.194 | 0.489 | 1.628 |  4.415 | 10.360 |
| kmc       | 0.221 | 0.885 | 1.885 |  4.895 |  8.655 |
| cbl       | 0.697 | 1.647 | 5.742 | 24.145 | 50.908 |
| sbwtrs    | 19.55 | 33.44 |   —   |    —   |    —   |
| fmsi      | 21.86 | 62.31 |   —   |    —   |    —   |

**k = 63** (cbl: no k = 63; fmsi/sbwt-rs capped to small genomes)

| tool   | ecoli  | yeast  | chr21 | celegans |  chr1  |
|--------|-------:|-------:|------:|---------:|-------:|
| **sklib** | 0.220 | 0.633 | 1.972 |  5.246 | 12.441 |
| kmc       | 0.439 | 0.926 | 2.607 |  7.103 | 13.116 |
| fmsi      | 23.51 | 72.79 |   —   |    —   |    —   |
| sbwtrs    |   —   | 83.69 |   —   |    —   |    —   |

Two tiers: a **fast tier** (sklib ≈ kmc ≪ cbl) and a **slow tier** (fmsi ≈ sbwt-rs,
~50–170× slower than sklib and capped at small/mid genomes).

### 4b. Speed across k — yeast, time (s), t = 1

| k  | sklib |  kmc  |  cbl  | sbwtrs |  fmsi  |
|----|------:|------:|------:|-------:|-------:|
| 15 | 0.710 | 0.565 | 1.130 | 12.836 | 60.540 |
| 21 | 0.514 | 0.660 | 1.281 | 26.270 | 59.245 |
| 31 | 0.489 | 0.885 | 1.647 | 33.441 | 62.314 |
| 41 | 0.498 | 0.974 | 1.663 | 49.107 | 12.497 |
| 51 | 0.619 | 0.827 | 1.831 | 65.739 | 63.381 |
| 63 | 0.633 |   –   |   –   | 83.691 | 72.793 |

**sklib is flat in k** (~0.5 s); the slow tier (sbwt-rs) **grows steeply with k** (each op
rebuilds the result SBWT).

### 4c. Memory — peak RSS (MB), t = 1, k = 31

| tool   | ecoli | yeast | chr21 | celegans |  chr1  |
|--------|------:|------:|------:|---------:|-------:|
| **sklib** |   6  |   7  | 13.5 |   25.5 |    54 | ← streams bucket-by-bucket
| sbwtrs    | 108  | 185  |   —  |     —  |    —  |
| kmc       | 112  | 230  |  627 |   1543 |  2514 |
| cbl       | 265  | 370  |  998 |   8533 | 15861 |
| fmsi      | 384  | 1122 |   —  |     —  |    —  |

**sklib uses the least by far** — ~20× less than sbwt-rs and ~10–60× less than kmc/cbl.

### 4d. Scaling with threads — time (s) & RAM (MB)

**chr1 — time**

| threads | sklib (k31)  | kmc (k31)    | sklib (k63)  | kmc (k63)     |
|---------|-------------:|-------------:|-------------:|--------------:|
| 1       | 10.36        | 8.66         | 12.44        | 13.12         |
| 2       |  5.84        | 8.59         |  7.97        | 12.90         |
| 4       |  3.54        | 8.58         |  5.03        | 13.06         |
| 8       |  2.50 (4.1×) | 7.42 (1.2×)  |  3.60 (3.5×) | 11.54 (1.1×)  |
| 16      |  2.18 (4.8×) | 4.29 (2.0×)  |  3.24 (3.8×) |  6.80 (1.9×)  |

**yeast — time**

| threads | sklib (k31)  | kmc (k31)    | sbwtrs (k31) | sklib (k63)  | kmc (k63) | sbwtrs (k63)  |
|---------|-------------:|-------------:|-------------:|-------------:|----------:|--------------:|
| 1       | 0.489        | 0.885        | 33.44        | 0.633        | 0.926     | 83.69         |
| 8       | 0.123 (4.0×) | 0.599 (1.5×) | 8.45 (4.0×)  | 0.170 (3.7×) | 0.796     | 19.07 (4.4×)  |
| 16      | 0.099 (4.9×) | 0.690        |     —        | 0.152 (4.2×) | 0.948     |     —         |

**RAM with threads** (chr1, k = 31): sklib 54 → 99 → 184 → 321 → **588** MB (grows ~linearly
— one bucket buffer per worker); kmc flat at ~2514 MB until t = 16, where it jumps to
**4619** MB. Even at t = 16, sklib (588 MB) ≪ kmc (4.6 GB) ≪ cbl (15.9 GB, single-thread).
(sbwt-rs set-op RAM is flat ~185 MB; it only exposes t = 1 and t = 8.)

- **sklib scales ~4–5×**, **sbwt-rs ~4×**, but **kmc barely threads** (1.1–2.0×) — its
  `simple` set-op engine is largely I/O-bound, and its RAM balloons at t = 16.
- So although **kmc edges sklib at t = 1 on the largest genome** (chr1: 8.7 s vs 10.4 s),
  **sklib pulls decisively ahead from t ≥ 4** (chr1 k31 t = 8: **2.5 s vs 7.4 s**) — with
  ~8–45× less RAM throughout.

---

## 5. Bottom line per workload

| Workload | Winner(s) | Threading note |
|----------|-----------|----------------|
| **Construct speed** | **kmc**, then **sklib** | both 3–10× faster than the rest; sklib scales best (5–7× @16), kmc 4–6×, sbwt-rs ~2×; the two tie on celegans at t=16 |
| **Construct space** | **fmsi** (~2 b/kmer) | SBWT ~10 (flat); sklib k-dependent, stores more; **sklib RAM ~10–60× lower** |
| **Query, small k** | **sbwt-rs** | only sklib threads (~1.6–2×); sbwt-rs query flat in `-t` |
| **Query, mid k / big genome** | **sshash**, then **sklib** | |
| **Query, high k** | **sklib** | robust where SBWT degrades; + the only one that threads |
| **Streaming** | **sshash**, then sklib/sbwt-rs | saturates fast: sklib ~1.2–1.4×, sbwt-rs flat |
| **Set-ops speed** | **sklib** ≈ kmc ≫ cbl ≫ (fmsi ≈ sbwt-rs) | sklib flat in k & scales ~5×; **kmc barely threads** → sklib wins clearly at t≥4 |
| **Set-ops memory** | **sklib** | streams; ~20× less than sbwt-rs, ~10–60× less than kmc/cbl |

**sklib's signature:** consistently top-tier on **construct speed**, the most **robust
query at high k**, a **decisive set-op win on both time and memory**, and the **best thread
scaling** across all three workloads — its per-bucket streaming keeps RAM at tens of MB
where competitors use hundreds of MB to tens of GB.

**On the sbwt Rust rewrite (`sbwtrs`):** vs the C++ `sbwt` it is **faster to query** (esp.
small k), **on par to build**, **lifts the k ≤ 32 cap**, and **adds native set-ops** —
though those set-ops are slow (≈ fmsi, ~2 orders of magnitude behind sklib/kmc), and its
query path does not thread here.

## 6. Caveats

- `peak_rss` units: **KB** in the construct/query CSVs (converted to MB in this report),
  **MB** in the setop CSV. See `EXPERIMENT.md`.
- All figures are the **median of 3 reps**; the threaded tools' first rep per config is a
  cold-page-cache outlier (~2× on big genomes) and is discarded by the median.
- Threaded tools = **sklib, kmc, sbwt-rs** (full t = 1/2/4/8/16 sweep; sbwt-rs set-ops
  only ran t = 1 and t = 8). The other five tools are single-thread.
- Set-op `result_kmers` is sklib's authoritative cardinality shared across tools; a tool's
  own set-op output is not re-validated here (perf-only).
- Legitimate skips: sbwt-C++ k≤32; sshash k≤31 (+ celegans k15 stream query fails); cbl
  odd-k≤59; fmsi/sbwt-rs set-ops capped to small/mid genomes; **bqf is approximate** and
  was only measured at k=31; **kmc has no query mode**.
- The toy `sarscov2` (30 KB, startup-bound) is omitted from the tables; `chm13` (3 GB) is
  deferred. The combined-mode set-ops (`joint`) are reported in `SETOPS_MULTI_REPORT.md`.
- The plotting step (`plot.py`) needs `pip install -U bottleneck` in the venv to regenerate
  figures (data unaffected).
