#!/usr/bin/env bash
# Experiment 5 — COST OF THE RE-CHAINING on the set-operation path.
#
# A materialized set operation costs (a) the per-bucket merge that decides which k-mers the
# result holds, (b) the ordering of those k-mers into the output payload, and (c) the write.
# The three CLI modes isolate them:
#   size      (--op *_size)   : (a) only, nothing is written;
#   nocompact (--no-compact)  : (a) + a plain sort of the single-k-mer records + (c) on a big file;
#   compact   (default)       : (a) + the re-chaining into virtual super-k-mers + (c) on a small file.
# At -t1 the library also reports the per-phase wall split (SKLIB_UNION_PHASE_TIMING), which is
# captured here: read / merge+collect / order / write. Together they answer "what does making the
# result usable actually cost", and the two output files give the space each choice pays.
#
#   bash benchmark/scripts/setop_rechain.sh
#   DATASETS=ecoli KM="31,15" JACCARD="0.5" THREADS="8" bash …/setop_rechain.sh
#
# Indexes (A and the per-J B's) and mutated genomes are cached, and shared with setop.sh.
# RESUMABLE: rows already present for the same (version,host,…) key are trusted and skipped.
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/lib.sh"; source "$SCRIPT_DIR/tools.sh"
need_tools "$TIME_BIN" python3
[[ -x "$SSKM_BIN" ]] || die "sskm not found at $SSKM_BIN"

# This experiment is sklib-only: the modes it compares are internal to sklib.
DATASETS="${DATASETS:-ecoli yeast celegans chr1}"
KM="${KM:-31,15 63,31}"
THREADS="${THREADS:-1 8}"
OPS="${OPS:-inter union diffab diffba}"
MODES="${MODES:-size nocompact compact}"
# Phase timing is a handful of clock reads per bucket and the library only honours it at -t1.
export SKLIB_UNION_PHASE_TIMING=1

CSV="${CSV:-$RESULTS/setop_rechain.csv}"
csv_init "$CSV" "timestamp,host,cpu,tool_version,dataset,k,m,threads,op,mode,jaccard_target,jaccard_measured,result_kmers,time_s,peak_rss_mb,out_bytes,out_records,records_per_kmer,bits_per_kmer,ph_read,ph_merge,ph_order,ph_write,ph_total"
load_done "$CSV" tool_version host dataset k m threads op mode jaccard_target
SCR="${RECHAIN_SCRATCH:-$RESULTS/rechain_scratch}"; mkdir -p "$SCR"
trap 'rm -f "$SCR"/out.sskm' EXIT

median() { printf '%s\n' "$@" | sort -g | awk '{a[NR]=$1} END{ if(NR==0){print "NA"} else print (NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2 }'; }

# Run a command REPS times; keep the median wall time, the max peak RSS, and the median of each
# phase timer. Unlike lib.sh's run_timed_median this keeps stderr, where the phase line is printed.
# Sets: R_SEC R_RSS_KB R_READ R_MERGE R_ORDER R_WRITE R_TOTAL
run_reps() {
    local i t0 t1 status=0 maxrss=0 tlog line
    local times=() reads=() merges=() orders=() writes=() totals=()
    tlog="$(mktemp)"
    for (( i = 0; i < REPS; i++ )); do
        t0=$(date +%s.%N)
        "$TIME_BIN" -v "$@" >/dev/null 2>"$tlog"
        status=$?
        t1=$(date +%s.%N)
        (( status == 0 )) || { rm -f "$tlog"; return "$status"; }
        times+=( "$(awk "BEGIN{printf \"%.3f\", $t1-$t0}")" )
        local rss; rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog"); [[ -n "$rss" ]] || rss=0
        (( rss > maxrss )) && maxrss="$rss"
        # [setop-phase materialize] total=0.23s  read=0.008s (3.4%)  merge+collect=0.057s (24%) …
        line=$(grep -m1 '\[setop-phase' "$tlog" 2>/dev/null)
        if [[ -n "$line" ]]; then
            reads+=(  "$(sed -n 's/.* read=\([0-9.e+-]*\)s .*/\1/p'          <<<"$line")" )
            merges+=( "$(sed -n 's/.* merge+collect=\([0-9.e+-]*\)s .*/\1/p' <<<"$line")" )
            orders+=( "$(sed -n 's/.* recompact=\([0-9.e+-]*\)s .*/\1/p'     <<<"$line")" )
            writes+=( "$(sed -n 's/.* write=\([0-9.e+-]*\)s .*/\1/p'         <<<"$line")" )
            totals+=( "$(sed -n 's/.*total=\([0-9.e+-]*\)s .*/\1/p'          <<<"$line")" )
        fi
    done
    rm -f "$tlog"
    R_SEC=$(median "${times[@]}"); R_RSS_KB="$maxrss"
    if (( ${#reads[@]} )); then
        R_READ=$(median "${reads[@]}");  R_MERGE=$(median "${merges[@]}")
        R_ORDER=$(median "${orders[@]}"); R_WRITE=$(median "${writes[@]}")
        R_TOTAL=$(median "${totals[@]}")
    else
        R_READ=NA; R_MERGE=NA; R_ORDER=NA; R_WRITE=NA; R_TOTAL=NA
    fi
    return 0
}

emit() {   # op mode result_kmers
    local rpk bpk
    if [[ "$OUT_BYTES" == NA ]]; then rpk=NA; bpk=NA
    else
        rpk=$(awk -v r="$OUT_RECORDS" -v n="$3" 'BEGIN{ if(n+0>0) printf "%.5f", r/n; else printf "NA" }')
        bpk=$(bits_per_kmer "$OUT_BYTES" "$3")
    fi
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date -Is)" "$HOST" "$CPU" "$tver" "$dataset" "$k" "$m" "$th" \
        "$1" "$2" "$Jt" "$Jmeas" "$3" "$R_SEC" "$(human_mb "$R_RSS_KB")" \
        "$OUT_BYTES" "$OUT_RECORDS" "$rpk" "$bpk" \
        "$R_READ" "$R_MERGE" "$R_ORDER" "$R_WRITE" "$R_TOTAL" >> "$CSV"
}

# One timed set operation in one mode. Sets R_* and OUT_BYTES / OUT_RECORDS (NA for the
# cardinality-only mode, which writes nothing).
OUT_BYTES=NA; OUT_RECORDS=NA
run_mode() {                                       # A B op mode th out
    local A="$1" B="$2" op="$3" mode="$4" th="$5" out="$6" a b sop
    case "$op" in
      inter)  a="$A"; b="$B"; sop=intersection ;;
      union)  a="$A"; b="$B"; sop=union        ;;
      diffab) a="$A"; b="$B"; sop=diff         ;;
      diffba) a="$B"; b="$A"; sop=diff         ;;
      *) return 1 ;;
    esac
    OUT_BYTES=NA; OUT_RECORDS=NA
    case "$mode" in
      compact)   run_reps "$SSKM_BIN" setop --op "$sop" -a "$a" -b "$b" -o "$out" -t "$th" || return 1 ;;
      nocompact) run_reps "$SSKM_BIN" setop --op "$sop" -a "$a" -b "$b" -o "$out" -t "$th" --no-compact || return 1 ;;
      size)      run_reps "$SSKM_BIN" setop --op "${sop}_size" -a "$a" -b "$b" -t "$th" || return 1 ;;
      *) return 1 ;;
    esac
    if [[ "$mode" != size ]]; then
        [[ -s "$out" ]] || return 1
        # Payload bytes = file minus the real header (40 + 16*n_buckets): the space the records take.
        local hdr
        hdr=$(python3 "$BENCH_HELPER" binheadersize "$out" 2>/dev/null || echo 40)
        OUT_BYTES=$(( $(stat -c%s "$out") - hdr ))
        OUT_RECORDS=$(python3 "$BENCH_HELPER" bincount "$out" 2>/dev/null || echo NA)
    fi
    return 0
}

tver="$(csv_escape "$(version_sklib)")"
log "setop_rechain: datasets=[$DATASETS] km=[$KM] jaccard=[$JACCARD] threads=[$THREADS] modes=[$MODES] reps=$REPS -> $CSV"
for dataset in $DATASETS; do
    prepare_genome "$dataset" || { warn "$dataset: skip (genome prepare failed)"; continue; }
    san="$SAN"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"
        ensure_index sklib "$dataset" "$k" "$m" || { warn "$dataset k=$k: A index failed, skip"; continue; }
        skA="$IDX_PATH"
        for Jt in $JACCARD; do
            rate="$(jaccard_to_rate "$Jt" "$k")"
            ensure_mutant "$san" "$rate" "$SEED" || { warn "$dataset k=$k J=$Jt: mutate failed, skip"; continue; }
            mut="$MUT_FA"; tag="${dataset}.mut${rate}"
            ensure_index_fa sklib "$mut" "$tag" "$k" "$m" || { warn "$dataset k=$k J=$Jt: B index failed, skip"; continue; }
            skB="$IDX_PATH"
            # Authoritative result cardinalities + measured Jaccard, from ONE combined pass.
            "$SSKM_BIN" setop -a "$skA" -b "$skB" --sizes > "$SCR/sizes.txt" 2>/dev/null || { warn "sizes failed"; continue; }
            RI=$(awk '$1=="intersection"{print $2}' "$SCR/sizes.txt")
            RU=$(awk '$1=="union"{print $2}'        "$SCR/sizes.txt")
            RDab=$(awk '$1=="diff_ab"{print $2}'    "$SCR/sizes.txt")
            RDba=$(awk '$1=="diff_ba"{print $2}'    "$SCR/sizes.txt")
            Jmeas=$(awk -v i="${RI:-0}" -v u="${RU:-0}" 'BEGIN{printf "%.4f",(u>0?i/u:0)}')
            declare -A RES=( [inter]="$RI" [union]="$RU" [diffab]="$RDab" [diffba]="$RDba" )
            log "=== $dataset k=$k m=$m J=$Jt (rate=$rate, measured J=$Jmeas): inter=$RI union=$RU A\\B=$RDab B\\A=$RDba ==="
            for th in $THREADS; do
                (( th > NPROC )) && continue
                for op in $OPS; do
                    res="${RES[$op]:-NA}"
                    for mode in $MODES; do
                        key="$(mk_key "$tver" "$HOST" "$dataset" "$k" "$m" "$th" "$op" "$mode" "$Jt")"
                        is_done "$key" && continue
                        rm -f "$SCR/out.sskm"
                        if run_mode "$skA" "$skB" "$op" "$mode" "$th" "$SCR/out.sskm"; then
                            emit "$op" "$mode" "$res"
                            mark_done "$key"
                            log "  t=$th $op/$mode: ${R_SEC}s  out=${OUT_BYTES}B rec=${OUT_RECORDS}  phases(r/m/o/w)=${R_READ}/${R_MERGE}/${R_ORDER}/${R_WRITE}"
                        else
                            warn "  t=$th $op/$mode: FAILED (skipped)"
                        fi
                        rm -f "$SCR/out.sskm"
                    done
                done
            done
        done
    done
done
log "${C_G}setop_rechain done${C_0} -> $CSV ($(($(wc -l < "$CSV")-1)) rows)"
