#!/usr/bin/env bash
# Scaling sweep: same overlap regime (B = A mutated 1%) across size classes ecoli..chr1.
# Constant overlap isolates the size axis. Appends to a master CSV. Run from repo root.
set -uo pipefail
cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
G=scripts/out/e2e/genomes
OUTDIR=scripts/out/bench; mkdir -p "$OUTDIR"
OUT="$OUTDIR/setops_scaling.csv"; LOG="$OUTDIR/setops_scaling.log"
: > "$OUT"; : > "$LOG"
RATE=0.01; SEED=7; K=21; M=11
# label:genome:reps  (fewer reps for the big ones)
pairs="ecoli:ecoliK12:3 yeast:yeast:3 chr21:chr21:3 chr20:chr20:2 celegans:celegans:2 chr1:chr1:1"
hdr=1
for p in $pairs; do
    label=${p%%:*}; rest=${p#*:}; g=${rest%%:*}; reps=${rest#*:}
    A="$G/$g.sanitized.fa"
    B="$G/${g}.mut1pct.fa"
    [[ -f "$A" ]] || { echo "[scaling] MISSING $A" >&2; continue; }
    if [[ ! -f "$B" ]]; then
        echo "[scaling] mutating $g (1%)" >&2
        python3 scripts/bench/mutate.py "$A" "$RATE" "$SEED" > "$B" || { echo "[scaling] mutate $g FAILED" >&2; continue; }
    fi
    echo "[scaling] === $label (reps=$reps) ===" >&2
    BENCH_CSV_HEADER=$hdr bash scripts/bench/bench_setops.sh "${label}_mut1" "$A" "$B" "$K" "$M" "$reps" >> "$OUT" 2>> "$LOG"
    hdr=0
    echo "[scaling] $label DONE" >&2
done
echo "[scaling] ALL DONE -> $OUT" >&2
