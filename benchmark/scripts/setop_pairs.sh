#!/usr/bin/env bash
# Experiment 6 — SET OPERATIONS ON REAL PAIRS OF K-MER SETS (sklib vs KMC), at scale.
#
# setop.sh pairs a genome with a mutated copy of itself to sweep the Jaccard index. This script
# instead takes PAIRS of real datasets (read sets catalogued in benchmark/data/reads.tsv and
# fetched by fetch_reads.sh, or any <name>.sanitized.fa in the genome cache): the overlap is
# whatever the data gives. Both tools index the SAME sanitized FASTA, every k-mer kept (KMC -ci1).
#
# Protocol (the sets reach several billion k-mers, so the laptop protocol is made explicit):
#  - COLD cache by default: before every timed run the inputs are evicted from the page cache
#    (posix_fadvise DONTNEED, no root needed) and every output is synced and deleted after it.
#    A KMC database is 8 B/k-mer at k=31 and 16 at k=63 and no longer fits in RAM.
#  - sklib and KMC alternate run by run (the order flips every rep), so thermal drift hits both.
#  - Pinning per thread count: PIN_T<th> (cpu list for taskset; empty = not pinned).
#  - Every row records %CPU and the file-system I/O reported by GNU time.
#  - A run is DISCARDED (and redone) if the laptop was paused or on battery while it ran.
#  - Correctness: |A|,|B| and every result cardinality from one sklib --sizes pass must equal
#    KMC's (kmc_tools info on its databases/outputs); sklib's own outputs are checked once.
#
#   bash benchmark/scripts/setop_pairs.sh                        # run (resumable)
#   bash benchmark/scripts/setop_pairs.sh pause | resume | status
#   PAIRS="hg002_A:hg002_B" KM="31,15" OPS="union" THREADS="8" REPS=1 bash …/setop_pairs.sh
#
# OPS ⊆ {union inter diffab diffba joint sizes}, each optionally "@<threads>" ("inter@8",
# "union@1-8") to override THREADS for that op; joint = ∩,∪,A∖B,B∖A materialized in one pass
# (sklib --*-out, KMC kmc_tools simple with four outputs); sizes = sklib --sizes (count only).
# RESUMABLE: per-rep rows already in the CSV for the same (tool, version, host, …) are skipped.
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"

RUN_DIR="${RUN_DIR:-$BENCH_REPO_ROOT/benchmark/results/runs/setop_reads_2026-10}"
PAUSE_FILE="${PAUSE_FILE:-$RUN_DIR/PAUSE}"
PAUSE_LOG="${PAUSE_LOG:-$RUN_DIR/logs/pause.log}"

# ---- pause / resume / status: tiny subcommands, usable while a run is going on ----
case "${1:-run}" in
    pause)  mkdir -p "$(dirname "$PAUSE_LOG")"; touch "$PAUSE_FILE"; echo "$(date -Is) pause" >> "$PAUSE_LOG"
            echo "paused: the current measurement is discarded and the run waits for 'resume'"; exit 0 ;;
    resume) rm -f "$PAUSE_FILE"; echo "$(date -Is) resume" >> "$PAUSE_LOG"; echo "resumed"; exit 0 ;;
    status) [[ -e "$PAUSE_FILE" ]] && echo "PAUSED" || echo "running (no pause flag)"
            tail -n 5 "$RUN_DIR/logs/setop_pairs.log" 2>/dev/null; exit 0 ;;
    run) ;;
    *) echo "usage: $0 [run|pause|resume|status]" >&2; exit 2 ;;
esac

source "$SCRIPT_DIR/lib.sh"; source "$SCRIPT_DIR/tools.sh"
need_tools kmc kmc_tools "$TIME_BIN" python3 taskset
[[ -x "$SSKM_BIN" ]] || die "sskm not found at $SSKM_BIN"

PAIRS="${PAIRS:-hg002_A:hg002_B gut_A:gut_B ocean_SRF:ocean_DCM}"
KM="${KM:-31,15 63,31}"
THREADS="${THREADS:-1 8}"
OPS="${OPS:-union joint sizes}"
REPS="${REPS:-3}"
CACHE="${CACHE:-cold}"                         # cold | warm (warm = no eviction)
CONSTRUCT_THREADS="${CONSTRUCT_THREADS:-8}"
PIN_T1="${PIN_T1-}"                            # t=1 not pinned: kmc_tools -t1 runs ~2.6 cores' worth
PIN_T8="${PIN_T8-0,1,3,6,8,10,12,17}"          # the 6 physical P-cores + one E-core per E cluster
KMC_RAM_GB="${KMC_RAM_GB:-12}"                 # = KMC's default -m
KMC_CLEANUP="${KMC_CLEANUP:-1}"                # drop a pair's KMC databases once all its cells are done
MIN_FREE_GB="${MIN_FREE_GB:-300}"              # refuse to start a (pair,k) with less free disk than this
WORK="${WORK:-$BENCH_REPO_ROOT/benchmark/results/latest/setop_reads}"
GEN="$BENCH_GEN_DIR"
CSV="${CSV:-$RUN_DIR/data/setop_runs.csv}"
CCSV="${CCSV:-$RUN_DIR/data/construct_runs.csv}"
SCSV="${SCSV:-$RUN_DIR/data/pair_sizes.csv}"
mkdir -p "$WORK/idx" "$WORK/out" "$WORK/tmp" "$RUN_DIR/data" "$RUN_DIR/logs"

csv_init "$CSV"  "timestamp,host,cpu,tool,tool_version,pair,set_a,set_b,k,m,threads,pin,cache,op,rep,n_a,n_b,result_kmers,jaccard,time_s,peak_rss_kb,cpu_pct,fs_in_blocks,fs_out_blocks,out_bytes,verified"
csv_init "$CCSV" "timestamp,host,cpu,tool,tool_version,set,k,m,threads,pin,cache,time_s,peak_rss_kb,cpu_pct,fs_in_blocks,fs_out_blocks,index_bytes,n_superkmers,input_bytes"
csv_init "$SCSV" "timestamp,host,tool_version,pair,set_a,set_b,k,m,sk_a,sk_b,kmc_a,kmc_b,inter,union,diff_ab,diff_ba,xor,jaccard,check"
SKV="$(version_sklib)"; KMV="kmc-$(kmc 2>&1 | sed -n 's/.*ver\. \([0-9.]*\).*/\1/p' | head -1)"
# A tool is "kmc", "sklib" (default --buckets) or a sklib variant "sklib@b<N>" (built with --buckets N):
# variants are separate tools for the run loop, so they are interleaved run by run like sklib/KMC.
is_sklib()   { [[ "${1%%@*}" == sklib ]]; }
buckets_of() { [[ "$1" == *@b* ]] && echo "${1##*@b}"; return 0; }
tver() { is_sklib "$1" && echo "$SKV" || echo "$KMV"; }

# ---- small helpers -----------------------------------------------------------
ac_online() { cat /sys/class/power_supply/AC/online 2>/dev/null || echo 1; }
pin_of() { local v="PIN_T$1"; echo "${!v-}"; }
pin_csv() { [[ -n "$1" ]] && echo "${1//,/;}" || echo none; }   # cpu list without commas (CSV-safe)
evict() {   # drop clean page-cache pages of the given files (dirty ones: sync first)
    python3 - "$@" <<'EOF'
import os, sys
for p in sys.argv[1:]:
    try:
        fd = os.open(p, os.O_RDONLY); os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED); os.close(fd)
    except OSError as e:
        print(f"evict {p}: {e}", file=sys.stderr)
EOF
}
kmc_files() { printf '%s\n' "$1.kmc_pre" "$1.kmc_suf"; }
kmc_count() { kmc_tools info "$1" 2>/dev/null | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }
sk_count()  { "$SSKM_BIN" setop -a "$1" -b "$1" --sizes -t 8 2>/dev/null | awk '$1=="A"{print $2}'; }
pause_mtime() { stat -c %Y "$PAUSE_LOG" 2>/dev/null || echo 0; }
op_threads() { [[ "$1" == *@* ]] && echo "${1#*@}" | tr '-' ' ' || echo "$THREADS"; }   # "inter@8" / "union@1-8"
free_gb() { df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9'; }

# Block until the laptop is ours: no pause flag and on AC power.
wait_ready() {
    local said=0
    while [[ -e "$PAUSE_FILE" || "$(ac_online)" != 1 ]]; do
        (( said )) || { log "waiting: $([[ -e "$PAUSE_FILE" ]] && echo 'paused by user' || echo 'on battery')"; said=1; }
        sleep 20
    done
    (( said )) && log "resuming"
    return 0
}

# measure <cpus> <cmd…>: one run under GNU time (pinned if cpus non-empty). Sets M_SEC M_RSS
# M_CPU M_FSIN M_FSOUT; stdout goes to $WORK/last.stdout. Returns the command status, or 99 when
# the run must be discarded (paused / on battery while it ran).
measure() {
    local cpus="$1"; shift
    local pre=() tlog t0 t1 st pm0
    [[ -n "$cpus" ]] && pre=(taskset -c "$cpus")
    tlog="$WORK/last.time"; pm0=$(pause_mtime)
    t0=$(date +%s.%N)
    "${pre[@]}" "$TIME_BIN" -v "$@" >"$WORK/last.stdout" 2>"$tlog"; st=$?
    t1=$(date +%s.%N)
    M_SEC=$(awk "BEGIN{printf \"%.3f\", $t1-$t0}")
    M_RSS=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog")
    M_CPU=$(awk -F': ' '/Percent of CPU this job got/{gsub(/%/,"",$2); print $2}' "$tlog")
    M_FSIN=$(awk -F': ' '/File system inputs/{print $2}' "$tlog")
    M_FSOUT=$(awk -F': ' '/File system outputs/{print $2}' "$tlog")
    if (( st != 0 )); then warn "command failed ($st): $*"; tail -n 5 "$tlog" >&2; return "$st"; fi
    if [[ -e "$PAUSE_FILE" || "$(ac_online)" != 1 || "$(pause_mtime)" != "$pm0" ]]; then
        warn "run discarded (paused or on battery during the measurement): $*"; return 99
    fi
    return 0
}

prep_cold() {   # flush pending writes, then evict the inputs
    sync
    [[ "$CACHE" == cold ]] && evict "$@"
    return 0
}

# ---- index construction (measured once, t=CONSTRUCT_THREADS, cold input) -----
declare -A IDX
build_index() {   # tool set k m -> IDX[tool:set:k]
    local tool="$1" set="$2" k="$3" m="$4" san="$GEN/$2.sanitized.fa" d idx key
    [[ -s "$san" ]] || { warn "$set: $san missing (run fetch_reads.sh)"; return 1; }
    local nb; nb="$(buckets_of "$tool")"
    if is_sklib "$tool"; then d="$WORK/idx/sklib/$SKV/$set.k$k.m$m${nb:+.b$nb}"; idx="$d/index.sskm"
    else d="$WORK/idx/kmc/$KMV/$set.k$k"; idx="$d/db"; fi
    IDX["$tool:$set:$k"]="$idx"
    [[ -f "$d/.done" ]] && return 0
    mkdir -p "$d"; rm -rf "$d/kt"
    local pin; pin="$(pin_of "$CONSTRUCT_THREADS")"
    while :; do
        wait_ready; prep_cold "$san"
        log "construct $tool $set k=$k m=$m (t=$CONSTRUCT_THREADS)"
        if is_sklib "$tool"; then
            measure "$pin" "$SSKM_BIN" construct -k "$k" -m "$m" -f "$san" -o "$idx" -t "$CONSTRUCT_THREADS" --tmp-dir "$WORK/tmp" ${nb:+--buckets "$nb"}
        else
            mkdir -p "$d/kt"
            measure "$pin" kmc -k"$k" -ci1 -m"$KMC_RAM_GB" -t"$CONSTRUCT_THREADS" -fm "$san" "$idx" "$d/kt"
        fi
        local st=$?; rm -rf "$d/kt"
        (( st == 99 )) && continue
        (( st == 0 )) || return 1
        break
    done
    local bytes nsk="NA"
    if is_sklib "$tool"; then bytes=$(stat -c%s "$idx"); nsk=$(python3 "$BENCH_HELPER" bincount "$idx" 2>/dev/null || echo NA)
    else bytes=$(( $(stat -c%s "$idx.kmc_pre") + $(stat -c%s "$idx.kmc_suf") )); fi
    key="$(mk_key "$tool" "$(tver "$tool")" "$HOST" "$set" "$k" "$m")"
    if [[ -z "${_CDONE[$key]:-}" ]]; then
        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$(date -Is)" "$HOST" "$CPU" "$tool" "$(tver "$tool")" \
            "$set" "$k" "$m" "$CONSTRUCT_THREADS" "$(pin_csv "$pin")" "$CACHE" "$M_SEC" "$M_RSS" "$M_CPU" "$M_FSIN" "$M_FSOUT" \
            "$bytes" "$nsk" "$(stat -c%s "$san")" >> "$CCSV"
        _CDONE[$key]=1
    fi
    touch "$d/.done"
    log "  -> ${M_SEC}s, $(( bytes / 1000000 )) MB"
}

# ---- one timed set operation -----------------------------------------------------
# run_op tool op th -> runs it (pinned), leaves outputs under $WORK/out/<tool>.<op>/
op_cmd() {   # tool op th A B od -> fills CMD array
    local tool="$1" op="$2" th="$3" A="$4" B="$5" od="$6"
    if is_sklib "$tool"; then
        case "$op" in
            union)  CMD=("$SSKM_BIN" setop --op union        -a "$A" -b "$B" -o "$od/u.sskm"   -t "$th") ;;
            inter)  CMD=("$SSKM_BIN" setop --op intersection -a "$A" -b "$B" -o "$od/i.sskm"   -t "$th") ;;
            diffab) CMD=("$SSKM_BIN" setop --op diff         -a "$A" -b "$B" -o "$od/dab.sskm" -t "$th") ;;
            diffba) CMD=("$SSKM_BIN" setop --op diff         -a "$B" -b "$A" -o "$od/dba.sskm" -t "$th") ;;
            joint)  CMD=("$SSKM_BIN" setop -a "$A" -b "$B" -t "$th" --inter-out "$od/i.sskm" --union-out "$od/u.sskm" \
                         --diff-ab-out "$od/dab.sskm" --diff-ba-out "$od/dba.sskm") ;;
            sizes)  CMD=("$SSKM_BIN" setop -a "$A" -b "$B" --sizes -t "$th") ;;
        esac
    else
        case "$op" in
            union)  CMD=(kmc_tools -t"$th" simple "$A" "$B" union "$od/u") ;;
            inter)  CMD=(kmc_tools -t"$th" simple "$A" "$B" intersect "$od/i") ;;
            diffab) CMD=(kmc_tools -t"$th" simple "$A" "$B" kmers_subtract "$od/dab") ;;
            diffba) CMD=(kmc_tools -t"$th" simple "$A" "$B" reverse_kmers_subtract "$od/dba") ;;
            joint)  CMD=(kmc_tools -t"$th" simple "$A" "$B" intersect "$od/i" union "$od/u" \
                         kmers_subtract "$od/dab" reverse_kmers_subtract "$od/dba") ;;
            sizes)  return 1 ;;
        esac
    fi
    return 0
}
outputs_of() { case "$1" in union) echo u;; inter) echo i;; diffab) echo dab;; diffba) echo dba;; joint) echo "i u dab dba";; esac; }

# ---- main ----------------------------------------------------------------------------
load_done "$CSV" tool tool_version host pair k m threads op rep cache
declare -A _RDONE=(); for _k in "${!_DONE[@]}"; do _RDONE[$_k]=1; done
load_done "$CCSV" tool tool_version host set k m
declare -gA _CDONE=(); for _k in "${!_DONE[@]}"; do _CDONE[$_k]=1; done
unset _k
TOOLS_SO="${TOOLS_SO:-sklib kmc}"
log "setop_pairs: pairs=[$PAIRS] km=[$KM] threads=[$THREADS] ops=[$OPS] reps=$REPS cache=$CACHE sklib=$SKV kmc=$KMV"
log "  pins: t1='${PIN_T1:-none}' t8='${PIN_T8:-none}'; csv=$CSV; pause with: $0 pause"

for km in $KM; do
  k="${km%%,*}"; m="${km##*,}"
  for pair in $PAIRS; do
    a="${pair%%:*}"; b="${pair##*:}"; pname="$a:$b"
    # pending cells of this (pair,k)?
    pending=0
    for tok in $OPS; do op="${tok%%@*}"; for th in $(op_threads "$tok"); do for (( r = 1; r <= REPS; r++ )); do for tool in $TOOLS_SO; do
        [[ "$tool" == kmc && "$op" == sizes ]] && continue
        [[ -n "${_RDONE[$(mk_key "$tool" "$(tver "$tool")" "$HOST" "$pname" "$k" "$m" "$th" "$op" "$r" "$CACHE")]:-}" ]] || pending=1
    done; done; done; done
    (( pending )) || { log "== $pname k=$k: all cells done"; continue; }
    log "== $pname k=$k m=$m"
    (( $(free_gb) >= MIN_FREE_GB )) || die "only $(free_gb) GB free under $WORK (< MIN_FREE_GB=$MIN_FREE_GB): free space, then re-run"

    ok=1
    for set in "$a" "$b"; do for tool in $TOOLS_SO; do build_index "$tool" "$set" "$k" "$m" || ok=0; done; done
    (( ok )) || { warn "$pname k=$k: index build failed, skip"; continue; }
    SKT=""; for t in $TOOLS_SO; do is_sklib "$t" && { SKT="$t"; break; }; done
    SA="${IDX[$SKT:$a:$k]}"; SB="${IDX[$SKT:$b:$k]}"; KA="${IDX[kmc:$a:$k]:-}"; KB="${IDX[kmc:$b:$k]:-}"

    # authoritative cardinalities (one sklib --sizes pass, not timed) + cross-check with KMC
    szf="$WORK/sizes.$a.$b.k$k.tsv"
    if [[ ! -s "$szf" ]]; then
        "$SSKM_BIN" setop -a "$SA" -b "$SB" --sizes -t 8 > "$szf.tmp" 2>/dev/null && mv "$szf.tmp" "$szf" || { warn "sizes failed"; continue; }
        get() { awk -v n="$1" '$1==n{print $2}' "$szf"; }
        ka=NA; kb=NA; [[ -n "$KA" ]] && ka=$(kmc_count "$KA"); [[ -n "$KB" ]] && kb=$(kmc_count "$KB")
        check=ok; [[ "$ka" == "$(get A)" && "$kb" == "$(get B)" ]] || check=MISMATCH
        J=$(awk -v i="$(get intersection)" -v u="$(get union)" 'BEGIN{printf "%.4f", (u>0?i/u:0)}')
        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$(date -Is)" "$HOST" "$SKV" "$pname" "$a" "$b" "$k" "$m" \
            "$(get A)" "$(get B)" "$ka" "$kb" "$(get intersection)" "$(get union)" "$(get diff_ab)" "$(get diff_ba)" "$(get xor)" "$J" "$check" >> "$SCSV"
    fi
    declare -A RES=( [union]="$(awk '$1=="union"{print $2}' "$szf")" [inter]="$(awk '$1=="intersection"{print $2}' "$szf")"
                     [diffab]="$(awk '$1=="diff_ab"{print $2}' "$szf")" [diffba]="$(awk '$1=="diff_ba"{print $2}' "$szf")" )
    RES[i]="${RES[inter]}"; RES[u]="${RES[union]}"; RES[dab]="${RES[diffab]}"; RES[dba]="${RES[diffba]}"
    RES[joint]=$(( ${RES[inter]} + ${RES[union]} + ${RES[diffab]} + ${RES[diffba]} )); RES[sizes]="${RES[union]}"
    NA_=$(awk '$1=="A"{print $2}' "$szf"); NB_=$(awk '$1=="B"{print $2}' "$szf")
    J=$(awk -v i="${RES[inter]}" -v u="${RES[union]}" 'BEGIN{printf "%.4f", (u>0?i/u:0)}')
    if [[ -n "$KA" && ( "$(kmc_count "$KA")" != "$NA_" || "$(kmc_count "$KB")" != "$NB_" ) ]]; then
        warn "$pname k=$k: |A|/|B| differ between sklib and KMC — NOT measuring this pair"; continue
    fi
    log "   |A|=$NA_ |B|=$NB_ |A∩B|=${RES[inter]} |A∪B|=${RES[union]} J=$J"

    declare -A VERIFIED=()
    for tok in $OPS; do
      op="${tok%%@*}"
      for th in $(op_threads "$tok"); do
        pin="$(pin_of "$th")"
        for (( r = 1; r <= REPS; r++ )); do
          order="$TOOLS_SO"; (( r % 2 == 0 )) && order="$(printf '%s\n' $TOOLS_SO | tac | tr '\n' ' ')"
          for tool in $order; do
            [[ "$tool" == kmc && "$op" == sizes ]] && continue
            key="$(mk_key "$tool" "$(tver "$tool")" "$HOST" "$pname" "$k" "$m" "$th" "$op" "$r" "$CACHE")"
            [[ -n "${_RDONE[$key]:-}" ]] && continue
            if is_sklib "$tool"; then A="${IDX[$tool:$a:$k]}"; B="${IDX[$tool:$b:$k]}"; inputs=("$A" "$B")
            else A="$KA"; B="$KB"; mapfile -t inputs < <(kmc_files "$KA"; kmc_files "$KB"); fi
            od="$WORK/out/$tool.$op"; CMD=(); st=1
            while :; do
                rm -rf "$od"; mkdir -p "$od"
                op_cmd "$tool" "$op" "$th" "$A" "$B" "$od" || break
                wait_ready; prep_cold "${inputs[@]}"
                measure "$pin" "${CMD[@]}"; st=$?
                (( st == 99 )) && continue
                break
            done
            (( st == 0 )) || { warn "$tool $op t=$th rep $r failed"; rm -rf "$od"; continue; }
            # verify the materialized cardinalities (KMC: header, every rep; sklib: once per op)
            ver=ok; outb=0
            for o in $(outputs_of "$op"); do
                if [[ "$tool" == kmc ]]; then
                    c=$(kmc_count "$od/$o"); outb=$(( outb + $(stat -c%s "$od/$o.kmc_pre") + $(stat -c%s "$od/$o.kmc_suf") ))
                    [[ "$c" == "${RES[$o]}" ]] || { ver="MISMATCH:$o=$c"; }
                else
                    outb=$(( outb + $(stat -c%s "$od/$o.sskm") ))
                    if [[ -z "${VERIFIED[$tool:$op:$o]:-}" ]]; then
                        c=$(sk_count "$od/$o.sskm"); [[ "$c" == "${RES[$o]}" ]] && VERIFIED[$tool:$op:$o]=1 || ver="MISMATCH:$o=$c"
                    fi
                fi
            done
            [[ "$op" == sizes ]] && { [[ "$(awk '$1=="union"{print $2}' "$WORK/last.stdout")" == "${RES[union]}" ]] || ver="MISMATCH:sizes"; }
            [[ "$ver" == ok ]] || warn "$tool $op k=$k t=$th: $ver (expected per sizes pass)"
            printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$(date -Is)" "$HOST" "$CPU" "$tool" "$(tver "$tool")" \
                "$pname" "$a" "$b" "$k" "$m" "$th" "$(pin_csv "$pin")" "$CACHE" "$op" "$r" "$NA_" "$NB_" "${RES[$op]}" "$J" \
                "$M_SEC" "$M_RSS" "$M_CPU" "$M_FSIN" "$M_FSOUT" "$outb" "$ver" >> "$CSV"
            _RDONE[$key]=1
            log "   $op t=$th rep $r $tool: ${M_SEC}s  rss=$(( ${M_RSS:-0} / 1024 ))MB  cpu=${M_CPU}%  out=$(( outb / 1000000 ))MB  $ver"
            rm -rf "$od"; sync
          done
        done
      done
    done

    # all cells of this (pair,k) done? then the KMC databases can go (sklib indexes are kept)
    left=0
    for tok in $OPS; do op="${tok%%@*}"; for th in $(op_threads "$tok"); do for (( r = 1; r <= REPS; r++ )); do for tool in $TOOLS_SO; do
        [[ "$tool" == kmc && "$op" == sizes ]] && continue
        [[ -n "${_RDONE[$(mk_key "$tool" "$(tver "$tool")" "$HOST" "$pname" "$k" "$m" "$th" "$op" "$r" "$CACHE")]:-}" ]] || left=1
    done; done; done; done
    if (( ! left && KMC_CLEANUP )) && [[ -n "$KA" ]]; then
        log "   $pname k=$k complete: removing its KMC databases"
        rm -rf "$(dirname "$KA")" "$(dirname "$KB")"
    fi
  done
done
log "setop_pairs: done"
