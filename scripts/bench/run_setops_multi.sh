#!/usr/bin/env bash
# Drive bench_setops_multi.sh across the real-genome ladder (E. coli -> human chr1) at the thread
# counts in THREADS. V1 (monothread): THREADS="1" (pinned to one core). V2: THREADS="1 4 8 <nproc>"
# (unpinned, all cores). One CSV for all rows.
#
#   THREADS="1" ./run_setops_multi.sh                 # V1 mono
#   THREADS="1 4 8 $(nproc)" ./run_setops_multi.sh    # V2 scaling
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$HERE/../.." && pwd)"
source "$HERE/lib.sh"

G="$BENCH_GEN_DIR"
K="${K:-31}"; M="${M:-21}"
read -r -a THREADS_LIST <<< "${THREADS:-1}"
OUT="${OUT:-$BENCH_REPO_ROOT/scripts/out/bench/setops_multi.csv}"
mkdir -p "$(dirname "$OUT")"
echo "label,k,m,tool,mode,action,threads,sec,rss_mb,result_kmers" > "$OUT"

# label  Afile  Bfile  reps   (files relative to the genome dir)
PAIRS=(
  "ecoliK12_Sakai   ecoliK12.sanitized.fa   ecoliSakai.sanitized.fa  3"
  "ecoliK12_UTI89   ecoliK12.sanitized.fa   ecoliUTI89.sanitized.fa  3"
  "ecoliSakai_IAI39 ecoliSakai.sanitized.fa ecoliIAI39.sanitized.fa  3"
  "chr21_chr22      chr21.sanitized.fa      chr22.sanitized.fa       2"
  "chr20_chr21      chr20.sanitized.fa      chr21.sanitized.fa       2"
  "chr1_mut1pct     chr1.sanitized.fa       chr1.mut1pct.fa          1"
)

for t in "${THREADS_LIST[@]}"; do
  pin=""; [[ "$t" == "1" ]] && pin="taskset -c 0"   # mono runs pinned; multi-core unpinned
  for spec in "${PAIRS[@]}"; do
    read -r label Af Bf reps <<< "$spec"
    # Optional ONLY="lab1 lab2 ..." filter to run a subset (e.g. the scaling-relevant pairs).
    [[ -n "${ONLY:-}" && " $ONLY " != *" $label "* ]] && continue
    A="$G/$Af"; B="$G/$Bf"
    [[ -f "$A" && -f "$B" ]] || { warn "MISSING $A or $B — skip $label (t=$t)"; continue; }
    log ">>> $label  t=$t  reps=$reps"
    "$HERE/bench_setops_multi.sh" "$label" "$A" "$B" "$K" "$M" "$reps" "$t" "$pin" >> "$OUT" 2>>"${OUT%.csv}.log"
  done
done
log "DONE -> $OUT"
