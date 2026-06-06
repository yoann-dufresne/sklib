#!/usr/bin/env bash
# Real downloaded genome pairs (no synthetic mutation): bacterial strain pairs and human
# chromosome pairs. Varied, biologically real overlap. Appends to a master CSV.
set -uo pipefail
cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
G=scripts/out/e2e/genomes
OUT=scripts/out/bench/setops_realpairs.csv; LOG=scripts/out/bench/setops_realpairs.log
: > "$OUT"; : > "$LOG"
K=21; M=11
# label:genomeA:genomeB:reps
pairs="
ecoliK12_Sakai:ecoliK12:ecoliSakai:3
ecoliK12_UTI89:ecoliK12:ecoliUTI89:3
ecoliSakai_IAI39:ecoliSakai:ecoliIAI39:3
chr21_chr22:chr21:chr22:2
chr20_chr21:chr20:chr21:2
"
hdr=1
for p in $pairs; do
    label=${p%%:*}; r=${p#*:}; gA=${r%%:*}; r=${r#*:}; gB=${r%%:*}; reps=${r#*:}
    A="$G/$gA.sanitized.fa"; B="$G/$gB.sanitized.fa"
    [[ -f "$A" && -f "$B" ]] || { echo "[real] MISSING $A or $B" >&2; continue; }
    echo "[real] === $label ===" >&2
    BENCH_CSV_HEADER=$hdr bash scripts/bench/bench_setops.sh "$label" "$A" "$B" "$K" "$M" "$reps" >> "$OUT" 2>> "$LOG"
    hdr=0
    echo "[real] $label DONE" >&2
done
echo "[real] ALL DONE -> $OUT" >&2
