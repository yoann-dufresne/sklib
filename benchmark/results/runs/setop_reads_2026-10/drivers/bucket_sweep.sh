#!/usr/bin/env bash
# Follow-up of the campaign: does sklib's --buckets matter at scale, and where is the optimum?
#
# At t=8 sklib's union throughput fell from 96 Mk/s (chr1, 0.28 G k-mers) to 53 Mk/s (ocean k=31,
# 6.6 G) while staying CPU-bound (785 % CPU); KMC stayed at 38-49 Mk/s. Hypothesis: with the
# default 4096 buckets a bucket grows with the data (≈68 k k-mers per input on chr1, ≈1 M on ocean)
# and its working set leaves the caches (L2 2 MB/core cluster, L3 24 MB shared by 8 threads).
#
# Each bucket count is a separate sklib "tool" (sklib@b<N>, built with --buckets N), so the
# variants alternate run by run exactly like sklib/KMC did. Cold cache, t=8 pinned, REPS=3 —
# the campaign protocol. KMC is not re-run (its numbers do not depend on sklib's buckets).
# Resumable: re-run to continue. Pause: bash benchmark/scripts/setop_pairs.sh pause / resume.
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RUN="$(dirname "$HERE")"
REPO="$(cd -- "$RUN/../../../.." && pwd)"
export RUN_DIR="$RUN"
export SSKM_BIN="${SSKM_BIN:-$REPO/build-bench/bin/sskm}"
export CSV="$RUN/data/bucket_sweep_runs.csv" CCSV="$RUN/data/bucket_sweep_construct_runs.csv"
export SCSV="$RUN/data/pair_sizes.csv" REPS=3 KMC_CLEANUP=0
SP="$REPO/benchmark/scripts/setop_pairs.sh"
LOG="$RUN/logs/bucket_sweep.log"
say() { printf '[sweep %s] %s\n' "$(date '+%m-%d %H:%M:%S')" "$*"; }
sweep() {   # label pairs km tools ops
    say "$1: pairs=[$2] km=$3 tools=[$4] ops=[$5]"
    PAIRS="$2" KM="$3" TOOLS_SO="$4" OPS="$5" bash "$SP" >> "$LOG" 2>&1 || say "$1: setop_pairs exited non-zero (see $LOG)"
}
BIG="sklib sklib@b16384 sklib@b65536 sklib@b262144"     # 4096 (default) … 262144
SMALL="sklib@b256 sklib@b1024 sklib sklib@b16384 sklib@b65536"
T8="union@8 joint@8 sizes@8"

# stage 1 — t=8, largest expected effect first
sweep "ocean k31"  "ocean_SRF:ocean_DCM" "31,15" "$BIG"   "$T8"
sweep "hg002 k31"  "hg002_A:hg002_B"     "31,15" "$BIG"   "$T8"
sweep "gut k31"    "gut_A:gut_B"         "31,15" "sklib@b1024 $BIG" "$T8"
sweep "ocean k63"  "ocean_SRF:ocean_DCM" "63,31" "$BIG"   "$T8"
# stage 2 — does it also move t=1 (where KMC leads 2.9x on ocean)?
sweep "ocean k31 t1" "ocean_SRF:ocean_DCM" "31,15" "$BIG"   "union@1"
# chr1 last: seconds-long runs are the most sensitive to a busy laptop (the 2026-10-06 morning
# chr1 rows were discarded for that reason), and the end of the queue runs at night.
sweep "chr1 k31"   "chr1:chr1_mutJ05"    "31,15" "$SMALL" "$T8"
sweep "chr1 k31 t1"  "chr1:chr1_mutJ05"    "31,15" "$SMALL" "union@1"
say "bucket sweep done"
