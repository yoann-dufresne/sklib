# Set operations on real read sets — sklib vs KMC, 2026-10

Set operations between **real pairs of k-mer sets**, 60× chr1 in bases, to show the SPIRE 2026
audience how the comparison holds at scale and on real data. The earlier campaigns pair a genome
with a mutated copy of itself (controlled Jaccard). Here the overlap comes from the data.

- **Started:** 2026-10-05 (download 17:22 CEST).
- **Machine:** `yoann-Precision-5490`, Intel Core Ultra 7 165H, 62 GiB, NVMe 1.9 TB. The run is
  on AC power and is paused (`setop_pairs.sh pause`) whenever the laptop is needed.
- **Tools:** sklib 0.15.0 at `0afc065` (`build-bench/bin/sskm`, Release, clang-18, native ISA) and
  KMC 3.2.4 (`kmc`, `kmc_tools`).
- **Scripts:** `benchmark/scripts/fetch_reads.sh` (data), `benchmark/scripts/setop_pairs.sh`
  (measurement) and `drivers/run_campaign.sh` (the order of the steps).

## Data

Six read sets of 13.5–15.8 Gbp each are catalogued in `benchmark/data/reads.tsv` and fetched from
ENA with an md5 check. They form three pairs:

| Pair | A | B | What the pair shows |
|---|---|---|---|
| `hg002_A:hg002_B` | HG002 SRR1766566 (flowcell 0018, library 2L1), 14.5 Gbp | HG002 SRR1766554 (flowcell 0029, library 2F1), 13.8 Gbp | One individual, two independent runs (GIAB HiSeq 2500 2×148, PRJNA200694). |
| `gut_A:gut_B` | HMP2 subject M2042, visits C2–C6 (5 runs), 13.5 Gbp | same subject, visits C19–C23 (5 runs), 15.4 Gbp | One gut, about a year apart (PRJNA398089, HiSeq 2000). |
| `ocean_SRF:ocean_DCM` | Tara station 023, surface, 15.8 Gbp | same station, DCM 55 m, 15.6 Gbp | High diversity, mostly singletons (PRJEB1787, 0.22–1.6 µm). |

Every FASTQ of a set (R1 and R2 of every run) goes into **one sanitized FASTA**
(`benchmark/data/genomes/<set>.sanitized.fa`): sequence lines are uppercased and split at every
non-ACGT run. Pairing is dropped, since a k-mer set ignores it. Both tools index that same file,
**every k-mer kept** (KMC `-ci1`, sklib has no count filter): these are raw read k-mer sets,
sequencing errors included.

## Protocol

- **Grid.**
  - Pairs × k/m ∈ {31/15, 63/31}, with k=31 for all pairs first.
  - Core: `union` and `joint` (∩, ∪, A∖B, B∖A materialized in one pass) at t=1 and t=8, plus
    sklib `--sizes` (cardinalities only).
  - E1: unitary `inter`, `diffab` and `diffba` at t=8.
  - REPS=3 for every cell. `data/setop_runs.csv` keeps **one row per rep**, and the medians are
    computed afterwards.
- **Cold cache.** The page cache is never warm:
  - before every timed run, `sync`, then the inputs are evicted (`posix_fadvise(DONTNEED)`);
  - after it, the outputs are deleted, then `sync`.

  At this scale a KMC database (8 B/k-mer at k=31, 16 at k=63) no longer fits in 62 GiB. A warm
  cache would serve sklib's smaller index from RAM and KMC's from disk, which nobody could
  reproduce. File-system I/O (GNU time, 512-byte blocks) is recorded in every row.
- **Construction.** Each index is built once, at t=8, from a cold input. Its time, peak RSS,
  %CPU, I/O and size are recorded in `data/construct_runs.csv`. sklib runs with
  `--buckets 4096` (the default) and `--tmp-dir` on the same NVMe; KMC with `-m12` (its default).
- **Pinning (laptop protocol).**
  - t=8 runs pinned to cpus `0,1,3,6,8,10,12,17` (the 6 physical P-cores plus one E-core per
    E-cluster), for both tools.
  - t=1 is not pinned, as in the deck. `kmc_tools -t1` keeps about 2.6 cores busy; the %CPU
    column shows it.
- **Interleaving.** Within a cell, sklib and KMC alternate run by run, and the order flips
  every rep.
- **Validity.** A run is discarded and redone if the pause flag was raised during it or if the
  laptop was on battery. The 2026-08-19 query campaign was lost to that. Sleep and lid-switch
  are inhibited for the whole campaign.
- **Correctness.**
  - One untimed sklib `--sizes` pass per (pair, k) gives |A|, |B|, |A∩B|, |A∪B|, |A∖B|, |B∖A| and
    |A△B| (`data/pair_sizes.csv`). |A| and |B| must equal KMC's (`kmc_tools info`), otherwise
    the pair is not measured.
  - Every KMC output is checked against these counts at every rep.
  - Every sklib output is checked once per operation, with a `--sizes` pass on the output.
  - The result goes in the `verified` column.
- **Drift anchor.** The deck's chr1 point (k=31, J=0.5 mutant, union, t=1/t=8, REPS=3) is
  measured with the same protocol at the start and at the end (`data/anchor_*_runs.csv`). It
  ties these numbers to the published ones, which were warm cache. A smoke run on 2026-10-05
  gave sklib 13.74 s and KMC 8.36 s at t=1, against 12.94 s (Aug) and 8.49 s (June) in the deck.

## Files

- `data/setop_runs.csv`: one row per (tool, pair, k, threads, op, rep).
- `data/construct_runs.csv`: one row per (tool, set, k).
- `data/pair_sizes.csv`: one row per (pair, k).
- `data/anchor_*`: the drift anchor.
- `logs/`: download, sanitize, measurement, pause log.
- `drivers/run_campaign.sh`: the campaign, resumable step by step.
