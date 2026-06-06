#!/usr/bin/env bash
# Overlap sweep at a fixed size (chr21, ~40 Mbp): B = chr21 mutated at increasing rates, so the
# intersection ratio drops from 100% (identical) toward ~0. Plus a real disjoint pair (chr21 vs
# chr22). Shows how each tool's time/output scale with the result size / overlap.
set -uo pipefail
cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
G=scripts/out/e2e/genomes
OUT=scripts/out/bench/setops_overlap.csv; LOG=scripts/out/bench/setops_overlap.log
: > "$OUT"; : > "$LOG"
K=21; M=11; BASE="$G/chr21.sanitized.fa"; SEED=13
hdr=1
emit() { # label B reps
    BENCH_CSV_HEADER=$hdr bash scripts/bench/bench_setops.sh "$1" "$BASE" "$2" "$K" "$M" "$3" >> "$OUT" 2>> "$LOG"
    hdr=0; echo "[overlap] $1 DONE" >&2
}
# identical (overlap 100%)
emit "ov_identical" "$BASE" 2
# graded mutation
for rate in 0.005 0.01 0.05 0.20; do
    B="$G/chr21.mut${rate}.fa"
    [[ -f "$B" ]] || { echo "[overlap] mutating chr21 @ $rate" >&2; python3 scripts/bench/mutate.py "$BASE" "$rate" "$SEED" > "$B"; }
    emit "ov_mut${rate}" "$B" 2
done
# real disjoint-ish (different human chromosome, same size class)
emit "ov_disjoint_chr22" "$G/chr22.sanitized.fa" 2
echo "[overlap] ALL DONE -> $OUT" >&2
