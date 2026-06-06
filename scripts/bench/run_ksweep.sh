#!/usr/bin/env bash
# k sweep on a fixed real pair (E. coli K12 vs Sakai): how time/output scale with k.
# CBL requires odd k in [1,59]; the wrapper builds a per-k CBL binary on demand.
set -uo pipefail
cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
G=scripts/out/e2e/genomes
OUT=scripts/out/bench/setops_ksweep.csv; LOG=scripts/out/bench/setops_ksweep.log
: > "$OUT"; : > "$LOG"
A="$G/ecoliK12.sanitized.fa"; B="$G/ecoliSakai.sanitized.fa"
# k:m  (k odd for CBL; m ~ k/2)
configs="15:7 21:11 31:15 41:19"
hdr=1
for c in $configs; do
    k=${c%%:*}; m=${c#*:}
    echo "[ksweep] === k=$k m=$m ===" >&2
    BENCH_CSV_HEADER=$hdr bash scripts/bench/bench_setops.sh "k${k}" "$A" "$B" "$k" "$m" 3 >> "$OUT" 2>> "$LOG"
    hdr=0
    echo "[ksweep] k=$k DONE" >&2
done
echo "[ksweep] ALL DONE -> $OUT" >&2
