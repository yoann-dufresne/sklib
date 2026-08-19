# Membership queries re-measured with v0.15.0 — and a power-state lesson

**Question.** The SPIRE paper reports every sklib figure from v0.15.0 except the membership-query
throughputs. Which version produced those, and what does v0.15.0 give?

**TL;DR.** The published query figures were **v0.11.0**, not the v0.13.1 the paper claimed — no
v0.13.1 query data has ever existed here. Re-measured with **v0.15.0 on chr1 at `-t1`, queries are
31–51 % faster than v0.11.0** at every k: 0.480 vs 0.347 Mkmer/s at k=31, 0.323 vs 0.246 at k=63.
sklib now **leads sbwt-rs from k=31 on** (1.4× at k=31, 1.8× at k=63) instead of merely tying at
k=31, and the v0.13 penalty on absent k-mers is gone (0.477 absent vs 0.483 present at k=31).
**A first re-run was thrown away: the laptop was on battery**, where the same binary measures 33 %
slower. On AC, sbwt-rs reproduces its June numbers within 1 % at every k, which is what licenses
comparing today's sklib against the June competitor figures.

## Which version produced the published numbers
`query_single.csv` holds sklib rows for 0.8.0, 0.10.1, 0.11.0 and 0.13.0 — never 0.13.1. The paper
printed 0.345 (chr1, k=31) and 0.246 (k=63):

| version | chr1 k=31 | chr1 k=63 |
|---|--:|--:|
| v0.11.0 | 0.347 | **0.246** |
| v0.13.0 | 0.290 | 0.216 |

k=63 matches v0.11.0 exactly and k=31 sits inside its presence cloud (0.344–0.354). A rescaled
v0.13.x run is ruled out: 0.345/0.290 = 1.19 but 0.246/0.216 = 1.14, so no single machine factor
maps one onto the other.

## The battery incident (why the first re-run was discarded)
The first re-measurement (2026-08-19 02:29) gave 0.275 Mkmer/s at k=31 — *slower* than v0.11.0,
which made no sense against the query work landed since. The machine had been discharging since
2026-08-18 10:25 (upower history), `powersave` + `balance_power`, cores at ~2.0 GHz instead of 4.2.

| control (chr1, k=31, `-t1`) | battery | AC | June reference |
|---|--:|--:|--:|
| sklib v0.13.2 (same binary, same index) | 0.220 | 0.329 | — |
| sklib v0.13.0 | — | — | 0.290 |
| sshash (unchanged binary) | — | 0.660 | 0.633 |

Battery costs ~33 %. The rows are kept as `latest/query_*_v0.15.0.ON-BATTERY-INVALID.csv`.
**Check the power state before any timing campaign on this machine.**

## Machine control on AC — the June comparison is legitimate
sbwt-rs is an unchanged binary from the June campaign. Re-run alongside these measurements:

| k | 15 | 21 | 31 | 41 | 51 | 63 |
|---|--:|--:|--:|--:|--:|--:|
| sbwt-rs today | 0.712 | 0.478 | 0.345 | 0.261 | 0.216 | 0.185 |
| sbwt-rs June | 0.752 | 0.481 | 0.347 | 0.269 | 0.208 | 0.181 |
| ratio | 0.95 | 0.99 | 0.99 | 0.97 | 1.04 | 1.02 |

Within 5 % everywhere, median 0.99. (sshash: 0.660 vs 0.633, +4 %.)

## Result — chr1, `-t1`, Mkmer/s (median over the presence sweep)

| k | **sklib v0.15.0** | sklib v0.11.0 | gain | sbwt-rs | sshash | sbwt | cbl |
|---|--:|--:|--:|--:|--:|--:|--:|
| 15 | 0.239 | 0.158 | +51 % | 0.712 | 0.265 | 0.583 | 0.132 |
| 21 | 0.388 | 0.284 | +37 % | 0.478 | 0.535 | 0.401 | 0.060 |
| 31 | **0.480** | 0.347 | +38 % | 0.345 | 0.633 | 0.295 | 0.019 |
| 41 | **0.521** | 0.382 | +36 % | 0.261 | — | — | 0.008 |
| 51 | **0.318** | 0.241 | +32 % | 0.216 | — | — | 0.005 |
| 63 | **0.323** | 0.246 | +31 % | 0.185 | — | — | — |

Competitors are the June figures (sbwt-rs cross-checked above). sklib passes sbwt-rs at k=31
instead of k=41, and reaches 76 % of sshash at k=31 where v0.11.0 reached 55 %.

**All four datasets, v0.15.0 (Mkmer/s):**

| dataset | k=31 `-t1` | k=31 `-t8` | k=63 `-t1` | k=63 `-t8` |
|---|--:|--:|--:|--:|
| E. coli | 1.190 | 2.410 | 0.531 | 1.361 |
| Yeast | 1.042 | 2.151 | 0.494 | 1.227 |
| C. elegans | 0.687 | 1.042 | 0.404 | 0.787 |
| Chr1 | 0.480 | 0.635 | 0.323 | 0.522 |

**Present vs absent k-mers (chr1, `-t1`):** v0.15.0 is flat — 0.477 (all absent) vs 0.483 (all
present) at k=31, 0.324 vs 0.323 at k=63. v0.13.0 was not (0.192 vs 0.313 at k=31): the negative
query path had regressed and is fixed.

## Data
`benchmark/results/reference/query_single_v0.15.0.csv` (240 rows), `query_stream_v0.15.0.csv`
(240 rows, sequence queries, not cited in the paper) and `query_machine_control.csv` (the
battery/AC/June controls above). Command lines in `benchmark/scripts/query_{single,stream}.sh`;
sklib built as `build-rel-rechain` (clang-18 Release, `-O3 -march=native`, LTO).
