#!/usr/bin/env bash
# Parallel cross-tool set-operation benchmark (sklib v0.6.0+ bucket-parallel set-ops).
#
# Since v0.6.0 `sskm setop` is parallel by bucket (-t). This sweeps sklib's threads against the two
# other set-op tools available here: KMC (kmc_tools simple, multithreaded) and FMSI (f-MS framework,
# serial). CBL is unavailable in this environment (skipped). Not pinned: multi-core runs use all
# cores, so the machine must be otherwise idle.
#
# Model (same as bench_setops.sh): build each index ONCE (separately measured), then run each set
# operation on the pre-built indexes. Result cardinality is sklib's `_size` (authoritative;
# cross-validated == KMC) and is recorded identically on every tool's op row.
#
# CSV: pair,k,stage,op,tool,threads,sec,rss_mb,result_kmers
#   stage ∈ build | setop ; op ∈ A|B (build) or inter|union|diff|*_size (setop).
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$HERE/../.." && pwd)"
source "$HERE/lib.sh"

K=21; M=11; NT="$(nproc)"
SKLIB_THREADS=(1 4 8 "$NT")          # sklib is the structure under test: full scaling curve
KMC_THREADS=(1 "$NT")                 # KMC: single vs all-core
FMSI_TIMEOUT="${FMSI_TIMEOUT:-900}"   # serial + RAM-heavy: cap each op at 15 min

SSKM="${SSKM:-$BENCH_REPO_ROOT/build-bench/bin/sskm}"
G="${SETOP_GEN_DIR:-$BENCH_GEN_DIR}"
FMSI_BIN="${FMSI_BIN:-}"; KMERCAMEL_BIN="${KMERCAMEL_BIN:-}"
[[ -x "$SSKM" ]] || die "sskm (release) not found at $SSKM"
need_tools kmc kmc_tools /usr/bin/time

OUT="${OUT:-$BENCH_REPO_ROOT/scripts/out/bench/setops_parallel.csv}"
LOG="${OUT%.csv}.log"; mkdir -p "$(dirname "$OUT")"
echo "pair,k,stage,op,tool,threads,sec,rss_mb,result_kmers" > "$OUT"; : > "$LOG"
emit(){ printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$1" "$K" "$2" "$3" "$4" "$5" "$6" "$(human_mb "$7")" "$8" >> "$OUT"; }

have_fmsi(){ [[ -n "$FMSI_BIN" && -x "$FMSI_BIN" && -n "$KMERCAMEL_BIN" && -x "$KMERCAMEL_BIN" ]]; }

# run_pair label A.fa B.fa reps do_fmsi
run_pair(){
    local label="$1" A="$2" B="$3" reps="$4" do_fmsi="$5"
    [[ -f "$A" && -f "$B" ]] || { warn "MISSING $A or $B — skip $label"; return; }
    local W="$BENCH_REPO_ROOT/scripts/out/bench/par/${label}"; rm -rf "$W"; mkdir -p "$W/kt"; local KT="$W/kt"
    RUN_REPS="$reps"
    log "=== pair $label  (reps=$reps, fmsi=$do_fmsi) ==="

    # ---- sklib: build once (construct is itself parallel; build at -t NT), then sweep set-op threads
    run_timed "$SSKM" construct -f "$A" -k $K -m $M -o "$W/A.sskm" -t "$NT"; emit "$label" build A sklib "$NT" "$RUN_SEC" "$RUN_RSS_KB" NA
    run_timed "$SSKM" construct -f "$B" -k $K -m $M -o "$W/B.sskm" -t "$NT"; emit "$label" build B sklib "$NT" "$RUN_SEC" "$RUN_RSS_KB" NA
    local rI rU rD
    rI=$("$SSKM" setop --op intersection_size -a "$W/A.sskm" -b "$W/B.sskm")
    rU=$("$SSKM" setop --op union_size        -a "$W/A.sskm" -b "$W/B.sskm")
    rD=$("$SSKM" setop --op diff_size         -a "$W/A.sskm" -b "$W/B.sskm")
    log "$label sizes: inter=$rI union=$rU diff=$rD"
    local t
    for t in "${SKLIB_THREADS[@]}"; do
        run_timed_median "$SSKM" setop --op intersection -a "$W/A.sskm" -b "$W/B.sskm" -o "$W/sI.sskm" -t "$t"; emit "$label" setop inter sklib "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rI"
        run_timed_median "$SSKM" setop --op union        -a "$W/A.sskm" -b "$W/B.sskm" -o "$W/sU.sskm" -t "$t"; emit "$label" setop union sklib "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rU"
        run_timed_median "$SSKM" setop --op diff         -a "$W/A.sskm" -b "$W/B.sskm" -o "$W/sD.sskm" -t "$t"; emit "$label" setop diff  sklib "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rD"
        run_timed_median "$SSKM" setop --op intersection_size -a "$W/A.sskm" -b "$W/B.sskm" -t "$t"; emit "$label" setop inter_size sklib "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rI"
        run_timed_median "$SSKM" setop --op union_size        -a "$W/A.sskm" -b "$W/B.sskm" -t "$t"; emit "$label" setop union_size sklib "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rU"
        run_timed_median "$SSKM" setop --op diff_size         -a "$W/A.sskm" -b "$W/B.sskm" -t "$t"; emit "$label" setop diff_size  sklib "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rD"
    done

    # ---- KMC: build + set-op at t=1 and t=NT
    for t in "${KMC_THREADS[@]}"; do
        run_timed kmc -k$K -ci1 -t"$t" -fm "$A" "$W/kA_$t" "$KT"; emit "$label" build A kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" NA
        run_timed kmc -k$K -ci1 -t"$t" -fm "$B" "$W/kB_$t" "$KT"; emit "$label" build B kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" NA
        run_timed_median kmc_tools -t"$t" simple "$W/kA_$t" "$W/kB_$t" intersect      "$W/kI_$t"; emit "$label" setop inter kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rI"
        run_timed_median kmc_tools -t"$t" simple "$W/kA_$t" "$W/kB_$t" union          "$W/kU_$t"; emit "$label" setop union kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rU"
        run_timed_median kmc_tools -t"$t" simple "$W/kA_$t" "$W/kB_$t" kmers_subtract "$W/kD_$t"; emit "$label" setop diff  kmc "$t" "$RUN_SEC" "$RUN_RSS_KB" "$rD"
    done

    # ---- FMSI (serial): kmercamel masked superstring + fmsi index (build), then inter/union/diff.
    # Build time/RSS = kmercamel + fmsi index (sum sec, max RSS), matching tools/fmsi.sh.
    if [[ "$do_fmsi" == "yes" ]] && have_fmsi; then
        local okfmsi=1 g src ks kr bs br
        for g in A B; do
            [[ "$g" == A ]] && src="$A" || src="$B"
            run_timed "$KMERCAMEL_BIN" compute -k $K -o /dev/null -M "$W/$g.ms.fa" "$src"; ks="$RUN_SEC"; kr="$RUN_RSS_KB"
            [[ -s "$W/$g.ms.fa" ]] || { warn "fmsi: kmercamel failed for $g on $label"; okfmsi=0; break; }
            run_timed "$FMSI_BIN" index -x -k $K "$W/$g.ms.fa" || { warn "fmsi: index failed for $g on $label"; okfmsi=0; break; }
            bs=$(awk -v a="$ks" -v b="$RUN_SEC" 'BEGIN{printf "%.3f", a+b}'); br=$(( kr > RUN_RSS_KB ? kr : RUN_RSS_KB ))
            emit "$label" build "$g" fmsi 1 "$bs" "$br" NA
        done
        if [[ "$okfmsi" == 1 ]]; then
            declare -A FM=( [inter]=inter [union]=union [diff]=diff ) FMR=( [inter]=$rI [union]=$rU [diff]=$rD )
            local op
            for op in inter union diff; do
                if run_timed timeout "$FMSI_TIMEOUT" "$FMSI_BIN" "${FM[$op]}" -p "$W/A.ms.fa" -p "$W/B.ms.fa" -r "$W/fm_$op" -k $K; then
                    emit "$label" setop "$op" fmsi 1 "$RUN_SEC" "$RUN_RSS_KB" "${FMR[$op]}"
                else
                    warn "fmsi $op on $label failed/timed out"; emit "$label" setop "$op" fmsi 1 NA NA "${FMR[$op]}"
                fi
            done
        fi
    elif [[ "$do_fmsi" == "yes" ]]; then
        warn "fmsi requested but FMSI_BIN/KMERCAMEL_BIN not set — skipping FMSI for $label"
    fi

    rm -rf "$W"   # free disk between pairs (indexes are large)
    log "=== $label DONE ==="
}

# label A B reps do_fmsi  (FMSI skipped on the largest pair: serial + multi-GB RAM)
run_pair ecoliK12_Sakai   "$G/ecoliK12.sanitized.fa"   "$G/ecoliSakai.sanitized.fa" 3 yes
run_pair ecoliK12_UTI89   "$G/ecoliK12.sanitized.fa"   "$G/ecoliUTI89.sanitized.fa" 3 yes
run_pair ecoliSakai_IAI39 "$G/ecoliSakai.sanitized.fa" "$G/ecoliIAI39.sanitized.fa" 3 yes
run_pair chr21_chr22      "$G/chr21.sanitized.fa"      "$G/chr22.sanitized.fa"      2 yes
run_pair chr20_chr21      "$G/chr20.sanitized.fa"      "$G/chr21.sanitized.fa"      2 no
run_pair chr1_mut1pct     "$G/chr1.sanitized.fa"       "$G/chr1.mut1pct.fa"         1 no
log "parallel set-op sweep DONE -> $OUT"
