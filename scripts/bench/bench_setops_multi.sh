#!/usr/bin/env bash
# Benchmark COMBINED (single-pass) vs SEQUENTIAL (one pass per op) k-mer set operations.
#
#   bench_setops_multi.sh <label> <A.fa> <B.fa> <k> <m> [reps] [threads] [pin]
#
# The combined feature produces all four relations (A∩B, A∪B, A\B, B\A) in ONE merge pass; the
# baseline runs the four ops separately (four passes). KMC's `kmc_tools simple` ALSO supports several
# ops in one pass, so it is the head-to-head competitor (combined vs combined, and each tool's own
# sequential baseline). Indexes are built ONCE (separately measured), then the set ops run on the
# pre-built indexes. Emits CSV rows to stdout (header once if BENCH_CSV_HEADER=1), logs to stderr.
#
# CSV: label,k,m,tool,mode,action,threads,sec,rss_mb,result_kmers
#   mode   ∈ combined | sequential        (single pass vs one pass per op)
#   action ∈ build | materialize | sizes  (build = index A+B; materialize = the 4 result lists;
#                                           sizes = the 4 cardinalities, no list written)
#   result_kmers = |A∩B|+|A∪B|+|A\B|+|B\A| (a single denominator; authoritative from sklib --sizes,
#                  cross-validated == KMC by tests/setop_multi_verif.sh)
set -uo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BENCH_REPO_ROOT="$(cd -- "$HERE/../.." && pwd)"
export BENCH_REPO_ROOT
# shellcheck source=/dev/null
source "$HERE/lib.sh"

[[ $# -ge 5 ]] || die "usage: bench_setops_multi.sh <label> <A.fa> <B.fa> <k> <m> [reps] [threads] [pin]"
# NOTE: ${8-default} (single dash) defaults only when $8 is UNSET, not when it is an empty string —
# the driver passes an empty pin for multi-core runs, and ${8:-default} would have wrongly substituted
# the default (pinning every run to one core, hiding all thread scaling).
LABEL="$1"; A_FA="$2"; B_FA="$3"; K="$4"; M="$5"; RUN_REPS="${6:-3}"; THREADS="${7:-1}"; PIN="${8-taskset -c 0}"

SSKM="${SSKM:-$BENCH_REPO_ROOT/build-bench/bin/sskm}"
[[ -x "$SSKM" ]] || die "sskm (release) not found at $SSKM"
need_tools kmc kmc_tools /usr/bin/time

WD="$BENCH_REPO_ROOT/scripts/out/bench/setops_multi/${LABEL}_k${K}_t${THREADS}"
rm -rf "$WD"; mkdir -p "$WD"; KT="$WD/kmc_tmp"; mkdir -p "$KT"

[[ "${BENCH_CSV_HEADER:-0}" == "1" ]] && \
    echo "label,k,m,tool,mode,action,threads,sec,rss_mb,result_kmers"

row() { # tool mode action sec rss_kb result_kmers threads
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$LABEL" "$K" "$M" "$1" "$2" "$3" "${7:-$THREADS}" "$4" "$(human_mb "$5")" "$6"
}
# Sum of several run_timed_median wall times (sequential baseline); RSS = max across the calls.
SEQ_SEC=0; SEQ_RSS=0
seq_reset(){ SEQ_SEC=0; SEQ_RSS=0; }
seq_add(){ SEQ_SEC=$(awk -v a="$SEQ_SEC" -v b="$RUN_SEC" 'BEGIN{printf "%.3f",a+b}'); (( RUN_RSS_KB>SEQ_RSS )) && SEQ_RSS="$RUN_RSS_KB"; }

log "=== $LABEL  k=$K m=$M  reps=$RUN_REPS threads=$THREADS pin='$PIN' ==="

# ---- build indexes once (sklib + kmc), separately measured ----
run_timed "$SSKM" construct -f "$A_FA" -k "$K" -m "$M" -o "$WD/A.sskm" -t "$THREADS"; row sklib build buildA "$RUN_SEC" "$RUN_RSS_KB" NA
run_timed "$SSKM" construct -f "$B_FA" -k "$K" -m "$M" -o "$WD/B.sskm" -t "$THREADS"; row sklib build buildB "$RUN_SEC" "$RUN_RSS_KB" NA
run_timed kmc -k"$K" -ci1 -t"$THREADS" -fm "$A_FA" "$WD/kA" "$KT" >/dev/null 2>&1; row kmc build buildA "$RUN_SEC" "$RUN_RSS_KB" NA
run_timed kmc -k"$K" -ci1 -t"$THREADS" -fm "$B_FA" "$WD/kB" "$KT" >/dev/null 2>&1; row kmc build buildB "$RUN_SEC" "$RUN_RSS_KB" NA

# ---- authoritative result cardinalities (one combined --sizes pass; validated == KMC) ----
"$SSKM" setop -a "$WD/A.sskm" -b "$WD/B.sskm" --sizes >"$WD/sizes.txt" 2>/dev/null
sz(){ awk -v k="$1" '$1==k{print $2}' "$WD/sizes.txt"; }
RI=$(sz intersection); RU=$(sz union); RDab=$(sz diff_ab); RDba=$(sz diff_ba)
RES=$(( RI + RU + RDab + RDba ))
log "sizes: inter=$RI union=$RU A\\B=$RDab B\\A=$RDba  (sum=$RES)"

# ================= sklib =================
# combined materialize: one call, all four lists
run_timed_median $PIN "$SSKM" setop -a "$WD/A.sskm" -b "$WD/B.sskm" -t "$THREADS" \
    --inter-out "$WD/sI.sskm" --union-out "$WD/sU.sskm" --diff-ab-out "$WD/sDab.sskm" --diff-ba-out "$WD/sDba.sskm"
row sklib combined materialize "$RUN_SEC" "$RUN_RSS_KB" "$RES"
# sequential materialize: four separate calls (B\A swaps the inputs)
seq_reset
run_timed_median $PIN "$SSKM" setop --op intersection -a "$WD/A.sskm" -b "$WD/B.sskm" -o "$WD/qI.sskm"   -t "$THREADS"; seq_add
run_timed_median $PIN "$SSKM" setop --op union        -a "$WD/A.sskm" -b "$WD/B.sskm" -o "$WD/qU.sskm"   -t "$THREADS"; seq_add
run_timed_median $PIN "$SSKM" setop --op diff         -a "$WD/A.sskm" -b "$WD/B.sskm" -o "$WD/qDab.sskm" -t "$THREADS"; seq_add
run_timed_median $PIN "$SSKM" setop --op diff         -a "$WD/B.sskm" -b "$WD/A.sskm" -o "$WD/qDba.sskm" -t "$THREADS"; seq_add
row sklib sequential materialize "$SEQ_SEC" "$SEQ_RSS" "$RES"
# combined sizes: one --sizes pass (no list written)
run_timed_median $PIN "$SSKM" setop -a "$WD/A.sskm" -b "$WD/B.sskm" --sizes -t "$THREADS"
row sklib combined sizes "$RUN_SEC" "$RUN_RSS_KB" "$RES"
# sequential sizes: four *_size calls
seq_reset
run_timed_median $PIN "$SSKM" setop --op intersection_size -a "$WD/A.sskm" -b "$WD/B.sskm" -t "$THREADS"; seq_add
run_timed_median $PIN "$SSKM" setop --op union_size        -a "$WD/A.sskm" -b "$WD/B.sskm" -t "$THREADS"; seq_add
run_timed_median $PIN "$SSKM" setop --op diff_size         -a "$WD/A.sskm" -b "$WD/B.sskm" -t "$THREADS"; seq_add
run_timed_median $PIN "$SSKM" setop --op diff_size         -a "$WD/B.sskm" -b "$WD/A.sskm" -t "$THREADS"; seq_add
row sklib sequential sizes "$SEQ_SEC" "$SEQ_RSS" "$RES"

# ================= KMC =================
# combined materialize: all four ops in ONE kmc_tools simple pass
run_timed_median $PIN kmc_tools -t"$THREADS" simple "$WD/kA" "$WD/kB" \
    intersect "$WD/kI" union "$WD/kU" kmers_subtract "$WD/kDab" reverse_kmers_subtract "$WD/kDba"
row kmc combined materialize "$RUN_SEC" "$RUN_RSS_KB" "$RES"
# sequential materialize: four separate kmc_tools simple calls
seq_reset
run_timed_median $PIN kmc_tools -t"$THREADS" simple "$WD/kA" "$WD/kB" intersect      "$WD/kI2";   seq_add
run_timed_median $PIN kmc_tools -t"$THREADS" simple "$WD/kA" "$WD/kB" union          "$WD/kU2";   seq_add
run_timed_median $PIN kmc_tools -t"$THREADS" simple "$WD/kA" "$WD/kB" kmers_subtract "$WD/kDab2"; seq_add
run_timed_median $PIN kmc_tools -t"$THREADS" simple "$WD/kB" "$WD/kA" kmers_subtract "$WD/kDba2"; seq_add
row kmc sequential materialize "$SEQ_SEC" "$SEQ_RSS" "$RES"

rm -rf "$WD/kI" "$WD/kU" "$WD"/k*2* "$WD"/*.sskm 2>/dev/null || true
log "=== $LABEL done ==="
