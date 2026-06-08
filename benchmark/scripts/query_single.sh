#!/usr/bin/env bash
# Experiment 2 — INDIVIDUAL k-mer QUERIES (one k-mer per FASTA record, scattered).
#
# Sweeps the present-fraction p ∈ PRESENCE (0..100%): each query set mixes p% k-mers
# sampled from the genome (present) with the rest random (absent), shuffled. Per
# (dataset, tool, k/m, threads) records time, peak RSS, throughput and the present/absent
# split — present-throughput is the p=100 row, absent-throughput the p=0 row.
#   bash benchmark/scripts/query_single.sh
#   DATASETS=ecoli KM="21,11" PRESENCE="0 50 100" THREADS="1 4" bash …/query_single.sh
#
# The index is built once per (tool,dataset,k,m) and cached, so an interrupted run
# restarts cheaply. RESUMABLE: rows already present for the same (tool,version,host) skip.
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/lib.sh"; source "$SCRIPT_DIR/tools.sh"
need_tools kmc kmc_tools /usr/bin/time python3

CSV="${CSV:-$RESULTS/query_single.csv}"
csv_init "$CSV" "timestamp,host,cpu,tool,tool_version,dataset,k,m,threads,presence,n_kmers,present_kmers,absent_kmers,time_s,peak_rss_kb,throughput_Mkmer_s"
load_done "$CSV" tool tool_version host dataset k m threads presence

log "query_single: datasets=[$DATASETS] tools=[$TOOLS] km=[$KM] presence=[$PRESENCE] threads=[$THREADS] -> $CSV"
for dataset in $DATASETS; do
    prepare_genome "$dataset" || { warn "$dataset: skip (genome prepare failed)"; continue; }
    san="$SAN"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"
        for tool in $TOOLS; do
            can_query_"$tool" || continue
            if ! ensure_index "$tool" "$dataset" "$k" "$m"; then
                warn "  $dataset k=$k m=$m $tool: index unavailable/unsupported (skip)"; continue; fi
            idx="$IDX_PATH"; tver="$(csv_escape "$(version_"$tool" 2>/dev/null)")"
            for p in $PRESENCE; do
                qset_single "$san" "$k" "$p" "$N_QUERY" "$SEED"
                total=$(( QSET_PRES + QSET_ABS ))
                for th in $(thread_list_"$tool"); do
                    (( th > NPROC )) && continue
                    key="$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th" "$p")"
                    is_done "$key" && continue
                    cpus="0"; (( th > 1 )) && cpus="0-$((th-1))"
                    RUN_REPS="$REPS"
                    if query_"$tool" "$idx" "$QSET_FA" "$k" "$m" "$cpus"; then
                        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
                            "$(date -Is)" "$HOST" "$CPU" "$tool" "$tver" "$dataset" "$k" "$m" "$th" "$p" \
                            "$total" "$QSET_PRES" "$QSET_ABS" "$RUN_SEC" "$RUN_RSS_KB" "$(mrate "$total" "$RUN_SEC")" >> "$CSV"
                        mark_done "$key"
                        log "  $dataset k=$k $tool p=$p% t=$th: $(mrate "$total" "$RUN_SEC") Mkmer/s"
                    else
                        warn "  $dataset k=$k $tool p=$p t=$th: query failed (skip)"
                    fi
                done
            done
        done
    done
done
log "${C_G}query_single done${C_0} -> $CSV ($(($(wc -l < "$CSV")-1)) rows)"
