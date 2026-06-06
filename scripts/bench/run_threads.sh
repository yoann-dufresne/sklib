#!/usr/bin/env bash
# Multi-core comparison: KMC is multithreaded by default; this measures how much it gains on the
# build and set-op steps at -t1 vs -t<all>, against sklib (single-thread v1) and CBL (single-thread
# by design). Not pinned (multi-core runs use all cores); machine must be otherwise idle.
#   CSV: pair,k,stage,op,tool,threads,sec,rss_mb,result_kmers
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$HERE/../.." && pwd)"
source "$HERE/lib.sh"
G="$BENCH_GEN_DIR"; SSKM="$BENCH_REPO_ROOT/build-bench/bin/sskm"
NT="$(nproc)"; K=21; M=11; RUN_REPS=3
OUT="$BENCH_REPO_ROOT/scripts/out/bench/setops_threads.csv"; LOG="${OUT%.csv}.log"
echo "pair,k,stage,op,tool,threads,sec,rss_mb,result_kmers" > "$OUT"; : > "$LOG"
emit(){ printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$1" "$K" "$2" "$3" "$4" "$5" "$6" "$(human_mb "$7")" "$8" >> "$OUT"; }

# pair: label A.fa B.fa
run_pair(){
    local label="$1" A="$2" B="$3"; local W="$BENCH_REPO_ROOT/scripts/out/bench/threads/${label}"
    rm -rf "$W"; mkdir -p "$W/kt"; local KT="$W/kt"
    log "=== threads pair $label ==="
    # sklib reference (single-thread)
    run_timed "$SSKM" construct -f "$A" -k $K -m $M -o "$W/A.sskm"; emit "$label" build A sklib 1 "$RUN_SEC" "$RUN_RSS_KB" NA
    run_timed "$SSKM" construct -f "$B" -k $K -m $M -o "$W/B.sskm"
    local rI rU rD
    rI=$("$SSKM" setop --op intersection_size -a "$W/A.sskm" -b "$W/B.sskm")
    rU=$("$SSKM" setop --op union_size        -a "$W/A.sskm" -b "$W/B.sskm")
    rD=$("$SSKM" setop --op diff_size         -a "$W/A.sskm" -b "$W/B.sskm")
    run_timed_median "$SSKM" setop --op intersection -a "$W/A.sskm" -b "$W/B.sskm" -o "$W/sI.sskm"; emit "$label" setop inter sklib 1 "$RUN_SEC" "$RUN_RSS_KB" "$rI"
    run_timed_median "$SSKM" setop --op union        -a "$W/A.sskm" -b "$W/B.sskm" -o "$W/sU.sskm"; emit "$label" setop union sklib 1 "$RUN_SEC" "$RUN_RSS_KB" "$rU"
    run_timed_median "$SSKM" setop --op diff         -a "$W/A.sskm" -b "$W/B.sskm" -o "$W/sD.sskm"; emit "$label" setop diff  sklib 1 "$RUN_SEC" "$RUN_RSS_KB" "$rD"
    run_timed_median "$SSKM" setop --op intersection_size -a "$W/A.sskm" -b "$W/B.sskm"; emit "$label" setop inter_size sklib 1 "$RUN_SEC" "$RUN_RSS_KB" "$rI"
    # KMC at t=1 and t=NT
    for t in 1 "$NT"; do
        run_timed kmc -k$K -ci1 -t"$t" -fm "$A" "$W/kA_$t" "$KT"; emit "$label" build A kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" NA
        run_timed kmc -k$K -ci1 -t"$t" -fm "$B" "$W/kB_$t" "$KT"
        run_timed_median kmc_tools -t"$t" simple "$W/kA_$t" "$W/kB_$t" intersect      "$W/kI_$t";  emit "$label" setop inter kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rI"
        run_timed_median kmc_tools -t"$t" simple "$W/kA_$t" "$W/kB_$t" union          "$W/kU_$t";  emit "$label" setop union kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rU"
        run_timed_median kmc_tools -t"$t" simple "$W/kA_$t" "$W/kB_$t" kmers_subtract "$W/kD_$t";  emit "$label" setop diff  kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rD"
    done
    rm -rf "$W/kt" "$W"/kA_* "$W"/kB_* "$W"/kI_* "$W"/kU_* "$W"/kD_* "$W"/*.sskm 2>/dev/null
}

run_pair ecoli   "$G/ecoliK12.sanitized.fa" "$G/ecoliSakai.sanitized.fa"
run_pair chr21x22 "$G/chr21.sanitized.fa"   "$G/chr22.sanitized.fa"
run_pair chr1    "$G/chr1.sanitized.fa"     "$G/chr1.mut1pct.fa"
log "threads sweep DONE -> $OUT"
