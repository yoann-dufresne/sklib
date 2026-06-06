#!/usr/bin/env bash
# Benchmark k-mer SET OPERATIONS (intersection / union / difference) across tools:
#   sklib (sskm setop)  vs  KMC (kmc_tools simple)  vs  CBL (inter/merge/diff)
#
#   bench_setops.sh <label> <A.fa> <B.fa> <k> <m> [reps] [tools]
#
# Emits CSV rows to stdout (header once if BENCH_CSV_HEADER=1), logs to stderr.
# Model (fair, apples-to-apples): build each index ONCE (separately measured), then run each set
# operation on the PRE-BUILT indexes. Everything is pinned to a single core (taskset) and KMC runs
# -t1, because sklib's set ops are sequential in v1 — so this is a single-thread comparison.
# Result cardinality is sklib's `_size` (no materialization; cross-validated == KMC == CBL); it feeds
# the throughput denominator and is recorded identically on every tool's op row.
#
# CSV columns:
#   label,kA,kB,k,m,tool,action,sec,rss_mb,result_kmers,bytes,thr_Mkmer_s,bits_per_kmer
# action ∈ build_A|build_B | inter|union|diff (materialize) | inter_size|union_size|diff_size (sklib).
set -uo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BENCH_REPO_ROOT="$(cd -- "$HERE/../.." && pwd)"
export BENCH_REPO_ROOT
# shellcheck source=/dev/null
source "$HERE/lib.sh"

[[ $# -ge 5 ]] || die "usage: bench_setops.sh <label> <A.fa> <B.fa> <k> <m> [reps] [tools]"
LABEL="$1"; A_FA="$2"; B_FA="$3"; K="$4"; M="$5"; RUN_REPS="${6:-3}"; TOOLS_SEL="${7:-sklib kmc cbl}"
PIN="${BENCH_PIN:-taskset -c 0}"

SSKM="${SSKM:-$BENCH_REPO_ROOT/build-bench/bin/sskm}"
[[ -x "$SSKM" ]] || die "sskm (release) not found at $SSKM"
need_tools kmc kmc_tools /usr/bin/time python3
# CBL source (per-k binary built on demand)
ENVF="$BENCH_REPO_ROOT/scripts/out/bench/tools_src/tools.env"
[[ -f "$ENVF" ]] && source "$ENVF"

WD="$BENCH_REPO_ROOT/scripts/out/bench/setops/${LABEL}_k${K}"
rm -rf "$WD"; mkdir -p "$WD"
KT="$WD/kmc_tmp"; mkdir -p "$KT"

[[ "${BENCH_CSV_HEADER:-0}" == "1" ]] && \
    echo "label,kA,kB,k,m,tool,action,sec,rss_mb,result_kmers,bytes,thr_Mkmer_s,bits_per_kmer"

row() { # tool action sec rss_kb result_kmers bytes inkmers
    local tool="$1" act="$2" sec="$3" rss="$4" res="$5" bytes="$6" ink="$7"
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$LABEL" "$KA" "$KB" "$K" "$M" "$tool" "$act" "$sec" "$(human_mb "$rss")" \
        "$res" "$bytes" "$(mrate "$ink" "$sec")" "$(bits_per_kmer "$bytes" "$res")"
}

# ---- CBL per-k binary (odd k only) ----
cbl_bin=""
cbl_build_binary() {
    [[ -n "${CBL_SRC:-}" && -d "${CBL_SRC:-}" ]] || { warn "cbl: CBL_SRC unset"; return 1; }
    (( K % 2 == 1 )) || { warn "cbl: K must be odd (k=$K); skipping CBL"; return 1; }
    local cache="$CBL_SRC/bin/cbl_k$K"
    if [[ ! -x "$cache" ]]; then
        warn "cbl: building K=$K (one-off)"; mkdir -p "$CBL_SRC/bin"
        local gv cxxinc; gv="$(g++ -dumpversion | cut -d. -f1)"
        cxxinc="-I/usr/include/c++/$gv -I/usr/include/x86_64-linux-gnu/c++/$gv"
        ( cd "$CBL_SRC" && K="$K" BINDGEN_EXTRA_CLANG_ARGS="$cxxinc" \
            cargo +nightly build --release --example cbl ) >"$CBL_SRC/bin/build_k$K.log" 2>&1 \
            || { warn "cbl: build failed (see bin/build_k$K.log)"; return 1; }
        cp "$CBL_SRC/target/release/examples/cbl" "$cache" || return 1
    fi
    cbl_bin="$cache"
}

log "=== $LABEL  k=$K m=$M  reps=$RUN_REPS  tools='$TOOLS_SEL' ==="

# ---------------- sklib (always: provides authoritative cardinalities) ----------------
run_timed "$SSKM" construct -f "$A_FA" -k "$K" -m "$M" -o "$WD/A.sskm"; sk_bA=$RUN_SEC; sk_rA=$RUN_RSS_KB
run_timed "$SSKM" construct -f "$B_FA" -k "$K" -m "$M" -o "$WD/B.sskm"; sk_bB=$RUN_SEC; sk_rB=$RUN_RSS_KB
KA=$("$SSKM" setop --op intersection_size -a "$WD/A.sskm" -b "$WD/A.sskm")
KB=$("$SSKM" setop --op intersection_size -a "$WD/B.sskm" -b "$WD/B.sskm")
INK=$((KA + KB))
# authoritative result cardinalities (fast, no materialization)
RES_inter=$("$SSKM" setop --op intersection_size -a "$WD/A.sskm" -b "$WD/B.sskm")
RES_union=$("$SSKM" setop --op union_size        -a "$WD/A.sskm" -b "$WD/B.sskm")
RES_diff=$( "$SSKM" setop --op diff_size         -a "$WD/A.sskm" -b "$WD/B.sskm")
log "sizes: |A|=$KA |B|=$KB  inter=$RES_inter union=$RES_union diff=$RES_diff"

if [[ " $TOOLS_SEL " == *" sklib "* ]]; then
    row sklib build_A "$sk_bA" "$sk_rA" "$KA" "$(stat -c%s "$WD/A.sskm")" "$KA"
    row sklib build_B "$sk_bB" "$sk_rB" "$KB" "$(stat -c%s "$WD/B.sskm")" "$KB"
    declare -A SK_OP=( [inter]=intersection [union]=union [diff]=diff )
    for op in inter union diff; do
        eval "res=\$RES_$op"
        # materialize
        run_timed_median $PIN "$SSKM" setop --op "${SK_OP[$op]}" -a "$WD/A.sskm" -b "$WD/B.sskm" -o "$WD/sk_$op.sskm"
        row sklib "$op" "$RUN_SEC" "$RUN_RSS_KB" "$res" "$(stat -c%s "$WD/sk_$op.sskm")" "$INK"
        # size-only (no materialization) — sklib-specific fast path
        run_timed_median $PIN "$SSKM" setop --op "${SK_OP[$op]}_size" -a "$WD/A.sskm" -b "$WD/B.sskm"
        row sklib "${op}_size" "$RUN_SEC" "$RUN_RSS_KB" "$res" 0 "$INK"
    done
fi

# ---------------- KMC ----------------
if [[ " $TOOLS_SEL " == *" kmc "* ]]; then
    run_timed kmc -k"$K" -ci1 -t1 -fm "$A_FA" "$WD/kA" "$KT"; km_bA=$RUN_SEC; km_rA=$RUN_RSS_KB
    run_timed kmc -k"$K" -ci1 -t1 -fm "$B_FA" "$WD/kB" "$KT"; km_bB=$RUN_SEC; km_rB=$RUN_RSS_KB
    kbytes() { echo $(( $(stat -c%s "$1.kmc_pre") + $(stat -c%s "$1.kmc_suf") )); }
    row kmc build_A "$km_bA" "$km_rA" "$KA" "$(kbytes "$WD/kA")" "$KA"
    row kmc build_B "$km_bB" "$km_rB" "$KB" "$(kbytes "$WD/kB")" "$KB"
    declare -A KM_OP=( [inter]=intersect [union]=union [diff]=kmers_subtract )
    for op in inter union diff; do
        eval "res=\$RES_$op"
        run_timed_median $PIN kmc_tools -t1 simple "$WD/kA" "$WD/kB" "${KM_OP[$op]}" "$WD/km_$op"
        row kmc "$op" "$RUN_SEC" "$RUN_RSS_KB" "$res" "$(kbytes "$WD/km_$op")" "$INK"
    done
fi

# ---------------- CBL ----------------
if [[ " $TOOLS_SEL " == *" cbl "* ]]; then
    if cbl_build_binary; then
        python3 "$BENCH_HELPER" dropshort "$A_FA" "$K" > "$WD/A.cbl.fa"
        python3 "$BENCH_HELPER" dropshort "$B_FA" "$K" > "$WD/B.cbl.fa"
        run_timed $PIN "$cbl_bin" build "$WD/A.cbl.fa" -o "$WD/A.cbl" --canonical; cb_bA=$RUN_SEC; cb_rA=$RUN_RSS_KB
        run_timed $PIN "$cbl_bin" build "$WD/B.cbl.fa" -o "$WD/B.cbl" --canonical; cb_bB=$RUN_SEC; cb_rB=$RUN_RSS_KB
        row cbl build_A "$cb_bA" "$cb_rA" "$KA" "$(stat -c%s "$WD/A.cbl")" "$KA"
        row cbl build_B "$cb_bB" "$cb_rB" "$KB" "$(stat -c%s "$WD/B.cbl")" "$KB"
        declare -A CB_OP=( [inter]=inter [union]=merge [diff]=diff )
        for op in inter union diff; do
            eval "res=\$RES_$op"
            run_timed_median $PIN "$cbl_bin" "${CB_OP[$op]}" "$WD/A.cbl" "$WD/B.cbl" -o "$WD/cb_$op.cbl"
            row cbl "$op" "$RUN_SEC" "$RUN_RSS_KB" "$res" "$(stat -c%s "$WD/cb_$op.cbl")" "$INK"
        done
    else
        warn "cbl: skipped"
    fi
fi

log "=== $LABEL done ==="
