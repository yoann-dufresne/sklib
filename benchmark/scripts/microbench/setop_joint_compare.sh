#!/usr/bin/env bash
# Focused KMC-vs-sklib comparison of the COMBINED ("joint") set operation only:
# the 4 operators of `kmc_tools simple` produced in one pass --
#   intersect, union, kmers_subtract (A\B), reverse_kmers_subtract (B\A).
# sklib equivalent: `sskm setop --inter-out --union-out --diff-ab-out --diff-ba-out`.
#
# Apples-to-apples == MATERIALIZE mode (KMC has no count-only mode). Reuses the whole
# benchmark harness (lib.sh/tools.sh) and only the joint-materialize path of setop.sh,
# so no unitary work is wasted. RESUMABLE (rows already present are trusted).
#
#   bash benchmark/scripts/microbench/setop_joint_compare.sh
#   DATASETS=ecoli KM="31,15" JACCARD="0.5" THREADS="1 8" bash …/setop_joint_compare.sh
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../../.." && pwd)"
# Capture any explicit user overrides BEFORE sourcing lib.sh (which fills these grids with
# its own wide defaults), so this script's focused defaults win when the user set nothing.
_U_DATASETS="${DATASETS:-}"; _U_KM="${KM:-}"; _U_JACCARD="${JACCARD:-}"; _U_THREADS="${THREADS:-}"; _U_TOOLS="${TOOLS:-}"
source "$SCRIPT_DIR/../lib.sh"; source "$SCRIPT_DIR/../tools.sh"
need_tools kmc kmc_tools "$TIME_BIN" python3
[[ -x "$SSKM_BIN" ]] || die "sskm not found at $SSKM_BIN"

# Focused defaults for this comparison (all env-overridable).
DATASETS="${_U_DATASETS:-chr21}"
KM="${_U_KM:-31,15 63,31 127,63}"
JACCARD="${_U_JACCARD:-0.1 0.5 0.9}"
THREADS="${_U_THREADS:-1 8}"
TOOLS="${_U_TOOLS:-sklib kmc}"

CSV="${CSV:-$RESULTS/setop_joint_compare.csv}"
csv_init "$CSV" "timestamp,host,cpu,tool,tool_version,dataset,k,m,threads,op,mode,joint,jaccard_target,jaccard_measured,result_kmers,time_s,peak_rss_mb,throughput_Mkmer_s"
load_done "$CSV" tool tool_version host dataset k m threads op mode jaccard_target
SCR="$SETOP_SCRATCH"; mkdir -p "$SCR"

emit_setop() {   # op mode joint res sec rss_kb
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date -Is)" "$HOST" "$CPU" "$tool" "$tver" "$dataset" "$k" "$m" "$th" \
        "$1" "$2" "$3" "$Jt" "$Jmeas" "$4" "$5" "$(human_mb "$6")" "$(mrate "$4" "$5")" >> "$CSV"
}

log "setop-joint-compare: datasets=[$DATASETS] tools=[$TOOLS] km=[$KM] jaccard=[$JACCARD] threads=[$THREADS] -> $CSV"
for dataset in $DATASETS; do
    prepare_genome "$dataset" || { warn "$dataset: skip (genome prepare failed)"; continue; }
    san="$SAN"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"; K_CUR="$k"
        ensure_index sklib "$dataset" "$k" "$m" || { warn "$dataset k=$k: sklib A index failed, skip"; continue; }
        skA="$IDX_PATH"
        for Jt in $JACCARD; do
            rate="$(jaccard_to_rate "$Jt" "$k")"
            ensure_mutant "$san" "$rate" "$SEED" || { warn "$dataset k=$k J=$Jt: mutate failed, skip"; continue; }
            mut="$MUT_FA"; tag="${dataset}.mut${rate}"
            ensure_index_fa sklib "$mut" "$tag" "$k" "$m" || { warn "$dataset k=$k J=$Jt: sklib B index failed, skip"; continue; }
            skB="$IDX_PATH"
            # authoritative cardinalities + measured Jaccard from ONE combined --sizes pass
            "$SSKM_BIN" setop -a "$skA" -b "$skB" --sizes > "$SCR/sizes.txt" 2>/dev/null || { warn "sizes failed"; continue; }
            RI=$(awk '$1=="intersection"{print $2}'  "$SCR/sizes.txt")
            RU=$(awk '$1=="union"{print $2}'         "$SCR/sizes.txt")
            RDab=$(awk '$1=="diff_ab"{print $2}'     "$SCR/sizes.txt")
            RDba=$(awk '$1=="diff_ba"{print $2}'     "$SCR/sizes.txt")
            Jmeas=$(awk -v i="${RI:-0}" -v u="${RU:-0}" 'BEGIN{printf "%.4f",(u>0?i/u:0)}')
            RES_JOINT=$(( ${RI:-0} + ${RU:-0} + ${RDab:-0} + ${RDba:-0} ))
            log "=== $dataset k=$k J=$Jt (rate=$rate, measured J=$Jmeas): inter=$RI union=$RU A\\B=$RDab B\\A=$RDba ==="

            for tool in $TOOLS; do
                can_setop_"$tool" || continue
                setop_has_joint_"$tool" || { warn "  $tool has no joint mode (skip)"; continue; }
                ensure_index "$tool" "$dataset" "$k" "$m" || { warn "  $tool A index failed (skip)"; continue; }; Aidx="$IDX_PATH"
                ensure_index_fa "$tool" "$mut" "$tag" "$k" "$m" || { warn "  $tool B index failed (skip)"; continue; }; Bidx="$IDX_PATH"
                tver="$(csv_escape "$(version_"$tool" 2>/dev/null)")"
                for th in $(thread_list_"$tool"); do
                    (( th > NPROC )) && continue
                    RUN_REPS="$REPS"
                    if ! is_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" joint materialize "$Jt")"; then
                        rm -rf "$SCR/j"; mkdir -p "$SCR/j"
                        if setop_joint_"$tool" "$Aidx" "$Bidx" materialize "$SCR/j" "$th"; then
                            emit_setop joint materialize 1 "$RES_JOINT" "$RUN_SEC" "$RUN_RSS_KB"
                            mark_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" joint materialize "$Jt")"
                            log "  ${C_G}$tool${C_0} t=$th J=$Jt: ${RUN_SEC}s rss=$(human_mb "$RUN_RSS_KB")MB"
                        else
                            warn "  $tool joint failed (t=$th J=$Jt)"
                        fi
                    fi
                done
            done
        done
    done
done
rm -rf "$SCR/j" 2>/dev/null || true
log "${C_G}setop-joint-compare done${C_0} -> $CSV ($(($(wc -l < "$CSV")-1)) rows)"
