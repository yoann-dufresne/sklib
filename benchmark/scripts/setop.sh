#!/usr/bin/env bash
# Experiment 4 — SET OPERATIONS.
#
# For each dataset A and a mutated copy B at a target Jaccard J ∈ JACCARD, benchmark the
# four operations {Union, Intersection, A\B, B\A} both UNITARY (one op per call) and JOINT
# (combined single pass), in MATERIALIZE and cardinality-only (SIZE) modes, across the
# thread sweep. sklib joint = `setop --*-out / --sizes`; KMC joint = `kmc_tools simple`
# with several ops; CBL/FMSI are single-op. Authoritative result cardinalities come from
# one sklib `--sizes` pass (cross-validated == KMC by tests/setop_multi_verif.sh).
#   bash benchmark/scripts/setop.sh
#   DATASETS=ecoli KM="21,11" JACCARD="0 0.5 1.0" THREADS="1 8" bash …/setop.sh
#
# Indexes (A and the per-J B's) and mutated genomes are cached. RESUMABLE: rows already
# present for the same (tool,version,host) are trusted and skipped on restart.
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/lib.sh"; source "$SCRIPT_DIR/tools.sh"
need_tools kmc kmc_tools "$TIME_BIN" python3
[[ -x "$SSKM_BIN" ]] || die "sskm not found at $SSKM_BIN"

CSV="${CSV:-$RESULTS/setop.csv}"
csv_init "$CSV" "timestamp,host,cpu,tool,tool_version,dataset,k,m,threads,op,mode,joint,jaccard_target,jaccard_measured,result_kmers,time_s,peak_rss_mb,throughput_Mkmer_s"
load_done "$CSV" tool tool_version host dataset k m threads op mode jaccard_target
SCR="$SETOP_SCRATCH"; mkdir -p "$SCR"

emit_setop() {   # op mode joint res sec rss_kb
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date -Is)" "$HOST" "$CPU" "$tool" "$tver" "$dataset" "$k" "$m" "$th" \
        "$1" "$2" "$3" "$Jt" "$Jmeas" "$4" "$5" "$(human_mb "$6")" "$(mrate "$4" "$5")" >> "$CSV"
}

log "setop: datasets=[$DATASETS] tools=[$TOOLS] km=[$KM] jaccard=[$JACCARD] threads=[$THREADS] -> $CSV"
for dataset in $DATASETS; do
    prepare_genome "$dataset" || { warn "$dataset: skip (genome prepare failed)"; continue; }
    san="$SAN"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"; K_CUR="$k"
        # sklib A index is needed for the authoritative sizes (and sklib is a tool too)
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
            declare -A RES=( [inter]="$RI" [union]="$RU" [diffab]="$RDab" [diffba]="$RDba" )
            log "=== $dataset k=$k J=$Jt (rate=$rate, measured J=$Jmeas): inter=$RI union=$RU A\\B=$RDab B\\A=$RDba ==="

            for tool in $TOOLS; do
                can_setop_"$tool" || continue
                ensure_index "$tool" "$dataset" "$k" "$m" || { warn "  $tool A index failed (skip)"; continue; }; Aidx="$IDX_PATH"
                ensure_index_fa "$tool" "$mut" "$tag" "$k" "$m" || { warn "  $tool B index failed (skip)"; continue; }; Bidx="$IDX_PATH"
                tver="$(csv_escape "$(version_"$tool" 2>/dev/null)")"
                for th in $(thread_list_"$tool"); do
                    (( th > NPROC )) && continue
                    RUN_REPS="$REPS"
                    # ---- unitary ops ----
                    for op in inter union diffab diffba; do
                        res="${RES[$op]:-NA}"
                        # materialize
                        if ! is_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" "$op" materialize "$Jt")"; then
                            rm -rf "$SCR/u"; mkdir -p "$SCR/u"
                            if setop_op_"$tool" "$Aidx" "$Bidx" "$op" materialize "$SCR/u/out" "$th"; then
                                emit_setop "$op" materialize 0 "$res" "$RUN_SEC" "$RUN_RSS_KB"
                                mark_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" "$op" materialize "$Jt")"
                            fi
                        fi
                        # size-only (sklib / cbl)
                        if setop_has_size_"$tool" && ! is_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" "$op" size "$Jt")"; then
                            if setop_op_"$tool" "$Aidx" "$Bidx" "$op" size "" "$th"; then
                                emit_setop "$op" size 0 "$res" "$RUN_SEC" "$RUN_RSS_KB"
                                mark_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" "$op" size "$Jt")"
                            fi
                        fi
                    done
                    # ---- joint (combined single pass) ----
                    if setop_has_joint_"$tool"; then
                        if ! is_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" joint materialize "$Jt")"; then
                            rm -rf "$SCR/j"; mkdir -p "$SCR/j"
                            if setop_joint_"$tool" "$Aidx" "$Bidx" materialize "$SCR/j" "$th"; then
                                emit_setop joint materialize 1 "$RES_JOINT" "$RUN_SEC" "$RUN_RSS_KB"
                                mark_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" joint materialize "$Jt")"
                            fi
                        fi
                        if setop_has_size_"$tool" && ! is_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" joint size "$Jt")"; then
                            if setop_joint_"$tool" "$Aidx" "$Bidx" size "$SCR/j" "$th"; then
                                emit_setop joint size 1 "$RES_JOINT" "$RUN_SEC" "$RUN_RSS_KB"
                                mark_done "$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" joint size "$Jt")"
                            fi
                        fi
                    fi
                done
            done
        done
    done
done
rm -rf "$SCR/u" "$SCR/j" 2>/dev/null || true
log "${C_G}setop done${C_0} -> $CSV ($(($(wc -l < "$CSV")-1)) rows)"
