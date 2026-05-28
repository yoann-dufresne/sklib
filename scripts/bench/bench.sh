#!/usr/bin/env bash
# Large-scale benchmark harness for sklib's sorted super-k-mer list.
#
# For each (dataset, tool, k, m, query-workload, thread-count) it constructs the
# index, measures construction (time, peak RSS, index size -> bits/k-mer, k-mers
# per super-k-mer) and query throughput per workload, optionally checks query
# correctness against KMC, and appends one tidy row per measurement to a CSV that
# scripts/bench/plots.py turns into the publication figures.
#
# Correctness (set-equality vs KMC) is the job of scripts/large_scale_e2e.sh; this
# harness is about performance/space curves and only runs a light false-negative gate.
#
# Usage:
#   bash scripts/bench/bench.sh
#   DATASETS="ecoli yeast" KM="21,7 21,9 21,11 21,13 31,15" \
#     WORKLOADS="positive random streaming shuffled reads" THREADS="1 4 22" \
#     TOOLS="sklib" REPS=3 bash scripts/bench/bench.sh
#
# Env knobs (defaults):
#   DATASETS="sarscov2 ecoli"   KM="21,11"   TOOLS="sklib"
#   WORKLOADS="positive random streaming shuffled reads"   THREADS="1"
#   N_QUERY=200000  STREAM_BP=2000000  READS_N=20000  READLEN=150  ERR=0.01
#   REPS=1  SEED=1234  SELFCHECK=1  FRESH=0  (FRESH=1 truncates the CSV first)
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
cd "$BENCH_REPO_ROOT"
source "$SCRIPT_DIR/lib.sh"
source "$SCRIPT_DIR/tools.sh"

# ---- configuration ----
DATASETS="${DATASETS:-sarscov2 ecoli}"
KM="${KM:-21,11}"
TOOLS="${TOOLS:-sklib}"
WORKLOADS="${WORKLOADS:-positive random streaming shuffled reads}"
THREADS="${THREADS:-1}"
N_QUERY="${N_QUERY:-200000}"
STREAM_BP="${STREAM_BP:-2000000}"
READS_N="${READS_N:-20000}"
READLEN="${READLEN:-150}"
ERR="${ERR:-0.01}"
REPS="${REPS:-1}"
SEED="${SEED:-1234}"
SELFCHECK="${SELFCHECK:-1}"
FRESH="${FRESH:-0}"

OUT="$BENCH_REPO_ROOT/scripts/out/bench"
QDIR="$OUT/queries"
WORK="$OUT/work"
CSV="$OUT/results.csv"
mkdir -p "$QDIR" "$WORK"

need_tools kmc kmc_tools /usr/bin/time python3 curl gzip awk sort taskset

# ---- CSV setup + run metadata ----
CSV_HEADER="timestamp,host,cpu,commit,tool,tool_version,dataset,dataset_bp,distinct_kmers,k,m,threads,phase,workload,n_queries,index_bytes,bytes_per_skmer,n_superkmers,kmers_per_superkmer,bits_per_kmer,time_s,peak_rss_kb,throughput_Mkmer_s,correctness"
[[ "$FRESH" == "1" ]] && rm -f "$CSV"
[[ -f "$CSV" ]] || echo "$CSV_HEADER" > "$CSV"

HOST="$(hostname)"
CPU="$(csv_escape "$(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo)")"
COMMIT="$(git -C "$BENCH_REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo NA)"
NPROC="$(nproc)"

emit_row() {
    # positional: tool tool_version dataset dataset_bp distinct_kmers k m threads
    #             phase workload n_queries index_bytes bytes_per_skmer n_superkmers
    #             kmers_per_superkmer bits_per_kmer time_s peak_rss_kb throughput correctness
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date -Is)" "$HOST" "$CPU" "$COMMIT" "$@" >> "$CSV"
}

# ---- query-set construction (cached per dataset+k) -------------------------
# Builds positive/random/streaming/shuffled/reads FASTA + ".nq" k-mer counts.
build_query_sets() {
    local dataset="$1" k="$2" san="$3" base="$QDIR/$dataset.$k"
    local wl fa
    if [[ ! -f "$base.streaming.fa" ]]; then
        log "$dataset k=$k: building query workloads"
        python3 "$BENCH_HELPER" sample_positive "$san" "$k" "$N_QUERY" "$SEED" | seq2fa > "$base.positive.fa"
        python3 "$BENCH_HELPER" randkmers "$k" "$N_QUERY" "$SEED"               | seq2fa > "$base.random.fa"
        python3 "$BENCH_HELPER" prefix    "$san" "$STREAM_BP"                            > "$base.streaming.fa"
        python3 "$BENCH_HELPER" allkmers  "$base.streaming.fa" "$k" \
            | python3 "$BENCH_HELPER" shuffle /dev/stdin "$SEED" | seq2fa     > "$base.shuffled.fa"
        python3 "$BENCH_HELPER" simreads  "$san" "$READLEN" "$READS_N" "$ERR" "$SEED"    > "$base.reads.fa"
        for wl in positive random streaming shuffled reads; do
            python3 "$BENCH_HELPER" count_kmers "$base.$wl.fa" "$k" > "$base.$wl.nq"
        done
        python3 "$BENCH_HELPER" count_kmers "$san" "$k" > "$base.total.nq"   # construct throughput denominator
    fi
}

# ---- correctness gate (sklib only): positive set must return no 0 ----------
CORR="NA"
correctness_gate() {
    local idx="$1" k="$2" base="$3" hits="$WORK/gate_hits.txt"
    CORR="NA"
    [[ "$SELFCHECK" == "1" ]] || return 0
    if "$SSKM_BIN" query -l "$idx" -i "$base.positive.fa" -o "$hits" 2>/dev/null; then
        if [[ -s "$hits" ]] && ! grep -q '0' "$hits"; then CORR="pass"; else CORR="fail"; fi
    else
        CORR="fail"
    fi
    rm -f "$hits"
}

# ====================================================================
#  Main
# ====================================================================
log "datasets=[$DATASETS] km=[$KM] tools=[$TOOLS] workloads=[$WORKLOADS] threads=[$THREADS]"
log "reps=$REPS n_query=$N_QUERY stream_bp=$STREAM_BP reads=${READS_N}x${READLEN}@${ERR} -> $CSV"
warn "query threads use CPU affinity (taskset); needs TBB-linked sskm to scale (else serial)."

declare -A DONE_TK   # (tool|dataset|k) already measured -> skip extra m for m-independent tools

for dataset in $DATASETS; do
    if ! prepare_genome "$dataset"; then warn "$dataset: skipped (prepare failed)"; continue; fi
    san="$SAN"
    dbp=$(grep -v '^>' "$san" | tr -d '\n' | wc -c)
    log "=== $dataset: $dbp ACGT bp ==="

    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"
        distinct=$(kmc_distinct_count "$san" "$k")
        [[ "$distinct" == "NA" || -z "$distinct" ]] && { warn "$dataset k=$k: KMC distinct-count failed, skipping"; continue; }
        build_query_sets "$dataset" "$k" "$san"
        qbase="$QDIR/$dataset.$k"

        for tool in $TOOLS; do
            "available_$tool" || { warn "$dataset k=$k m=$m: tool '$tool' unavailable, skipping"; continue; }
            # m-independent tools (competitors) are built once per (dataset,k), not per m.
            tk="$tool|$dataset|$k"
            if ! "uses_m_$tool" && [[ -n "${DONE_TK[$tk]:-}" ]]; then continue; fi
            DONE_TK[$tk]=1
            tver="$(csv_escape "$("version_$tool")")"
            wdir="$WORK/${tool}_${dataset}_${k}_${m}"; rm -rf "$wdir"; mkdir -p "$wdir"

            # ---- construct (timed) ----
            log "$dataset k=$k m=$m [$tool] construct"
            RUN_REPS=1
            if ! "construct_$tool" "$san" "$k" "$m" "$wdir"; then
                warn "$dataset k=$k m=$m [$tool] construct failed/unsupported"; rm -rf "$wdir"; continue
            fi
            local_total=$(cat "$qbase.total.nq" 2>/dev/null || echo 0)
            bpk=$(bits_per_kmer "$IDX_PAYLOAD_BYTES" "$distinct")
            kps=$(div "$distinct" "${N_SKMERS:-0}")
            bps=$(div "$IDX_PAYLOAD_BYTES" "${N_SKMERS:-0}")
            ctput=$(mrate "$local_total" "$RUN_SEC")
            emit_row "$tool" "$tver" "$dataset" "$dbp" "$distinct" "$k" "$m" 1 \
                construct - - "$IDX_FILE_BYTES" "$bps" "${N_SKMERS:-NA}" "$kps" "$bpk" \
                "$RUN_SEC" "$RUN_RSS_KB" "$ctput" NA

            # ---- correctness gate (sklib) ----
            [[ "$tool" == "sklib" ]] && correctness_gate "$IDX_PATH" "$k" "$qbase"

            # ---- query workloads x threads ----
            for wl in $WORKLOADS; do
                qfa="$qbase.$wl.fa"; nq=$(cat "$qbase.$wl.nq" 2>/dev/null || echo 0)
                [[ -s "$qfa" ]] || { warn "  workload $wl missing, skip"; continue; }
                for th in $THREADS; do
                    (( th > NPROC )) && continue
                    local_cpus="0-$((th-1))"; [[ "$th" == "1" ]] && local_cpus="0"
                    RUN_REPS="$REPS"
                    if "query_$tool" "$IDX_PATH" "$qfa" "$k" "$m" "$local_cpus"; then
                        local_corr=NA; [[ "$wl" == "positive" ]] && local_corr="$CORR"
                        emit_row "$tool" "$tver" "$dataset" "$dbp" "$distinct" "$k" "$m" "$th" \
                            query "$wl" "$nq" "$IDX_FILE_BYTES" "$bps" "${N_SKMERS:-NA}" "$kps" "$bpk" \
                            "$RUN_SEC" "$RUN_RSS_KB" "$(mrate "$nq" "$RUN_SEC")" "$local_corr"
                        printf '  %s/%s th=%s: %s Mkmer/s (%.3fs, %s kmers) %s\n' \
                            "$wl" "$tool" "$th" "$(mrate "$nq" "$RUN_SEC")" "$RUN_SEC" "$nq" "$local_corr" >&2
                    else
                        warn "  query $wl/$tool th=$th failed"
                    fi
                done
            done
            rm -rf "$wdir"
        done
    done
done

log "${C_G}done${C_0} -> $CSV"
log "rows: $(($(wc -l < "$CSV") - 1))"
