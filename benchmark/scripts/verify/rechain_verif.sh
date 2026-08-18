#!/usr/bin/env bash
# Correctness side of the re-chaining study (benchmark/scripts/setop_rechain.sh).
#
# The timing campaign compares a compacted set-operation result with the `--no-compact` one.
# That comparison is only meaningful if both are the SAME SET and both are usable, so this
# script checks, per dataset / k / operation:
#   1. same set        : |compact △ no-compact| = 0 and |compact ∩ no-compact| = |compact|
#   2. same queries    : `sskm query` over the source genome gives byte-identical answers
#   3. closed          : each result can feed a further set operation (result ∩ B), same
#                        cardinality from the compacted and the uncompacted operand
#   4. closure quality : re-chaining A ∪ A reproduces the record count of A's own index,
#                        i.e. a materialized result is as compact as a constructed one
#
#   bash benchmark/scripts/verify/rechain_verif.sh
#   DATASETS="ecoli yeast" KM="31,15" bash …/rechain_verif.sh
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BENCH_REPO_ROOT="$(cd -- "$SCRIPT_DIR/../../.." && pwd)"
source "$BENCH_REPO_ROOT/benchmark/scripts/lib.sh"
source "$BENCH_REPO_ROOT/benchmark/scripts/tools.sh"
[[ -x "$SSKM_BIN" ]] || die "sskm not found at $SSKM_BIN"

DATASETS="${DATASETS:-ecoli yeast celegans}"
KM="${KM:-31,15 63,31}"
OPS="${OPS:-union inter diffab}"
JV="${JV:-0.5}"                     # the pair used for checks 1-3 (J=1 is used for check 4)
WD="${RECHAIN_VERIF_WD:-$RESULTS/rechain_verif}"; mkdir -p "$WD"
PASS=0; FAIL=0
ok()   { printf '  \033[32mOK\033[0m   %s\n' "$*"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$*"; FAIL=$((FAIL+1)); }

card() {                                # <op> <A> <B> -> cardinality on stdout
    "$SSKM_BIN" setop --op "$1" -a "$2" -b "$3" -t 8 2>/dev/null | awk '{print $NF}' | tail -1
}
records() { python3 "$BENCH_HELPER" bincount "$1" 2>/dev/null || echo NA; }

for dataset in $DATASETS; do
    prepare_genome "$dataset" || { warn "$dataset: genome prepare failed"; continue; }
    san="$SAN"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"
        ensure_index sklib "$dataset" "$k" "$m" || { warn "$dataset k=$k: A index failed"; continue; }
        skA="$IDX_PATH"; recA="$(records "$skA")"
        rate="$(jaccard_to_rate "$JV" "$k")"
        ensure_mutant "$san" "$rate" "$SEED" || { warn "mutate failed"; continue; }
        ensure_index_fa sklib "$MUT_FA" "${dataset}.mut${rate}" "$k" "$m" || { warn "B index failed"; continue; }
        skB="$IDX_PATH"
        echo "== $dataset k=$k m=$m (J=$JV) =="
        for op in $OPS; do
            case "$op" in inter) sop=intersection;; union) sop=union;; diffab) sop=diff;; esac
            C="$WD/c.sskm"; N="$WD/n.sskm"
            "$SSKM_BIN" setop --op "$sop" -a "$skA" -b "$skB" -o "$C" -t 8 >/dev/null 2>&1 || { bad "$op: compact run failed"; continue; }
            "$SSKM_BIN" setop --op "$sop" -a "$skA" -b "$skB" -o "$N" -t 8 --no-compact >/dev/null 2>&1 || { bad "$op: no-compact run failed"; continue; }

            # 1. same set — one combined pass gives |C △ N|, |C| and |N|
            "$SSKM_BIN" setop -a "$C" -b "$N" --sizes > "$WD/sz.txt" 2>/dev/null
            x=$(awk '$1=="xor"{print $2}' "$WD/sz.txt")
            cc=$(awk '$1=="A"{print $2}' "$WD/sz.txt"); nn=$(awk '$1=="B"{print $2}' "$WD/sz.txt")
            if [[ "$x" == 0 && -n "$cc" && "$cc" == "$nn" ]]; then ok "$op: same set (|C △ N|=0, |C|=|N|=$cc)"
            else bad "$op: sets differ (xor=$x, |C|=$cc, |N|=$nn)"; fi

            # 2. same answers to the same queries
            "$SSKM_BIN" query -l "$C" -i "$san" -o "$WD/qc.txt" -t 8 >/dev/null 2>&1
            "$SSKM_BIN" query -l "$N" -i "$san" -o "$WD/qn.txt" -t 8 >/dev/null 2>&1
            if cmp -s "$WD/qc.txt" "$WD/qn.txt"; then ok "$op: identical query answers over $(basename "$san")"
            else bad "$op: query answers differ"; fi

            # 3. each result can feed a further operation
            fc=$(card intersection_size "$C" "$skB"); fn=$(card intersection_size "$N" "$skB")
            if [[ -n "$fc" && "$fc" == "$fn" ]]; then ok "$op: result feeds a further set operation (|result ∩ B|=$fc, both operands)"
            else bad "$op: further operation differs (compact=$fc, no-compact=$fn)"; fi
            rm -f "$C" "$N" "$WD/qc.txt" "$WD/qn.txt" "$WD/sz.txt"
        done

        # 4. closure quality: A ∪ A re-chained must be exactly A's own index
        rate1="$(jaccard_to_rate 1.0 "$k")"
        ensure_mutant "$san" "$rate1" "$SEED" >/dev/null 2>&1
        ensure_index_fa sklib "$MUT_FA" "${dataset}.mut${rate1}" "$k" "$m" >/dev/null 2>&1
        skA2="$IDX_PATH"
        "$SSKM_BIN" setop --op union -a "$skA" -b "$skA2" -o "$WD/uu.sskm" -t 8 >/dev/null 2>&1
        recU="$(records "$WD/uu.sskm")"
        if [[ "$recU" == "$recA" ]]; then ok "A ∪ A re-chains to A's own record count ($recA)"
        else bad "A ∪ A gives $recU records, A's index has $recA"; fi
        rm -f "$WD/uu.sskm"
    done
done
printf '\n%s passed, %s failed\n' "$PASS" "$FAIL"
[[ "$FAIL" == 0 ]]
