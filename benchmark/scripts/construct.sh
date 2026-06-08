#!/usr/bin/env bash
# Experiment 1 — INDEX CONSTRUCTION.
#
# For every (dataset, tool, k/m) and every thread count in the sweep, build the index
# and record: time, peak RSS, index size, bits/k-mer, throughput.
#   bash benchmark/scripts/construct.sh
#   DATASETS="ecoli" KM="21,11 31,15" THREADS="1 4" TOOLS="sklib kmc" bash …/construct.sh
#
# Tools that ignore threads (everything but sklib/kmc) are measured once (t=1) and shown
# flat against the others. Per-(tool,k) feasibility is enforced by the wrappers: an
# unsupported k makes construct_<tool> fail and the cell is skipped + warned.
#
# RESUMABLE: the CSV is append-only; on restart, rows already present for the same
# (tool, tool_version, host) are trusted and skipped — re-run after an interrupt without
# rebuilding what is already measured.
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/lib.sh"; source "$SCRIPT_DIR/tools.sh"
need_tools kmc kmc_tools /usr/bin/time python3

CSV="${CSV:-$RESULTS/construct.csv}"
csv_init "$CSV" "timestamp,host,cpu,tool,tool_version,dataset,k,m,threads,time_s,peak_rss_kb,index_bytes,bits_per_kmer,n_kmers,n_superkmers,throughput_Mkmer_s"
load_done "$CSV" tool tool_version host dataset k m threads
WORK="$RESULTS/work/construct"; mkdir -p "$WORK"

log "construct: datasets=[$DATASETS] tools=[$TOOLS] km=[$KM] threads=[$THREADS] -> $CSV"
for dataset in $DATASETS; do
    prepare_genome "$dataset" || { warn "$dataset: skip (genome prepare failed)"; continue; }
    san="$SAN"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"
        distinct="$(kmc_distinct_count "$san" "$k")"
        [[ "$distinct" =~ ^[0-9]+$ ]] || { warn "$dataset k=$k: distinct-count failed, skip"; continue; }
        for tool in $TOOLS; do
            can_construct_"$tool" || continue
            tver="$(csv_escape "$(version_"$tool" 2>/dev/null)")"
            for th in $(thread_list_"$tool"); do
                (( th > NPROC )) && continue
                key="$(mk_key "$tool" "$tver" "$HOST" "$dataset" "$k" "$m" "$th")"
                is_done "$key" && continue
                wd="$WORK/${tool}_${dataset}_${k}_${m}_t${th}"; rm -rf "$wd"; mkdir -p "$wd"
                RUN_REPS=1
                if CONSTRUCT_THREADS="$th" construct_"$tool" "$san" "$k" "$m" "$wd"; then
                    bpk="$(bits_per_kmer "${IDX_PAYLOAD_BYTES:-0}" "$distinct")"
                    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
                        "$(date -Is)" "$HOST" "$CPU" "$tool" "$tver" "$dataset" "$k" "$m" "$th" \
                        "$RUN_SEC" "$RUN_RSS_KB" "${IDX_FILE_BYTES:-NA}" "$bpk" "$distinct" \
                        "${N_SKMERS:-NA}" "$(mrate "$distinct" "$RUN_SEC")" >> "$CSV"
                    mark_done "$key"
                    log "  $dataset k=$k m=$m $tool t=$th: ${RUN_SEC}s ${IDX_FILE_BYTES:-?}B ${bpk}b/kmer"
                else
                    warn "  $dataset k=$k m=$m $tool t=$th: construct failed/unsupported (skipped)"
                fi
                rm -rf "$wd"
            done
        done
    done
done
log "${C_G}construct done${C_0} -> $CSV ($(($(wc -l < "$CSV")-1)) rows)"
