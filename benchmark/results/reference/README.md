# Reference — frozen golden baselines

Curated CSV/TSV snapshots and figures used as **regression anchors** and as the raw data behind
the [`../journals/`](../journals/). This is **not** a campaign run (for the cross-tool campaigns
see [`../runs/`](../runs/)); these numbers are promoted by hand from a known-good build and change
only deliberately. They are consumed by the producer/microbench harnesses (as digest & timing
baselines) and quoted in the journals.

## Set-operation CSVs
| File | What |
|---|---|
| `setops_realpairs.csv` | set ops on real genome pairs (ecoli variants, yeast, celegans, human chr). |
| `setops_overlap.csv` | set ops swept over Jaccard overlap (build_A / build_B / inter / inter_size). |
| `setops_scaling.csv` | scaling across k (15→51) on paired genomes. |
| `setops_ksweep.csv` | k-mer-size sweep (15→63). |
| `setops_threads.csv` | thread scaling (t ∈ 1,2,4,8,16,22). |
| `setops_parallel.csv` | threaded set ops, sklib `-t` vs KMC vs FMSI (backs `SETOPS_REPORT.md` §7). |
| `setops_multi_v1.csv` / `setops_multi_v2.csv` | combined single-pass set ops, single- / multi-thread. |
| `setops_wide_width.csv` | wide-store (`__uint128`/`kuint256`) merge timings (backs `SETOPS_WIDE_WIDTH.md`). |

## Bottleneck / profiling CSVs
| File | What |
|---|---|
| `bottleneck_decompose.csv` | set-op `_size` (merge-only) vs materialize split, per pair. |
| `bottleneck_construct.csv` | construction timing on E. coli variants. |
| `bottleneck_impl_results.csv` | before/after speedups for individual bottleneck fixes. |
| `bottleneck_reprofile_cycinst.csv` · `bottleneck_reprofile_postopt.csv` | cycle/instruction reprofiles. |
| `bottleneck_sweep_{k,m,buckets}.csv` | one-axis sweeps of the bottleneck measurement. |
| `bottleneck_perf_summary.txt` | formatted `perf stat` phase breakdown. |

## Producer
| File | What |
|---|---|
| `producer_digest.tsv` | **bit-exact FNV-1a digests** of the super-k-mer producer output (chr21, celegans at k=21/31/63). The `producer/producer_{bench,median}.sh` gate fails on any mismatch — it pins producer correctness across optimizations. |

## Figures (`figs/`)
`fig_time_vs_size.png`, `fig_ram_vs_size.png`, `fig_time_vs_overlap.png`, `fig_threads.png` —
the set-op figures embedded in [`../journals/SETOPS_REPORT.md`](../journals/SETOPS_REPORT.md).
