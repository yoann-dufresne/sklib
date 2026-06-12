#!/usr/bin/env bash
# Interleaved back-to-back A/B of two union_bench binaries on one (A,B) pair. Interleaving one timed
# rep from each binary per round cancels slow thermal/turbo/cache drift that corrupts "fresh run vs
# old table" comparisons (which is unreliable for the noisy k63 / __uint128 store). Reports each
# binary's median + a robust delta with a noise band, so a gain is only called when it clears the band.
#
#   benchmark/scripts/union_ab.sh <binA> <binB> <A.sskm> <B.sskm> [rounds] [label]
#   PIN=5 benchmark/scripts/union_ab.sh build-union/.../union_bench /tmp/ub.cand A.sskm B.sskm 15
#
# binA = reference (best known), binB = candidate. Positive delta% => candidate SLOWER (regression).
set -uo pipefail
binA="$1"; binB="$2"; A="$3"; B="$4"; ROUNDS="${5:-15}"; LABEL="${6:-ab}"
PIN="${PIN:-5}"
PINCMD=(); command -v taskset >/dev/null && PINCMD=(taskset -c "$PIN")

sample() {  # bin -> one timed median (reps=1, no warmup) in seconds on stdout
    "${PINCMD[@]}" "$1" --a "$A" --b "$B" --mode bench --warmup 0 --reps 1 --out "/dev/shm/ab_$$.sskm" 2>/dev/null \
        | awk -F'\t' '$1=="RESULT"{for(i=2;i<=NF;i++){split($i,a,"=");if(a[1]=="median_s")print a[2]}}'
}
median() { printf '%s\n' "$@" | sort -g | awk '{x[NR]=$1} END{print (NR%2)?x[(NR+1)/2]:(x[NR/2]+x[NR/2+1])/2}'; }
mad() {     # median absolute deviation around the median (arg1)
    local med="$1"; shift
    printf '%s\n' "$@" | awk -v m="$med" '{d=$1-m; if(d<0)d=-d; print d}' | sort -g \
        | awk '{x[NR]=$1} END{print (NR%2)?x[(NR+1)/2]:(x[NR/2]+x[NR/2+1])/2}'
}

# warm both (discarded)
sample "$binA" >/dev/null; sample "$binB" >/dev/null
declare -a SA=() SB=()
for ((r=0; r<ROUNDS; r++)); do
    SA+=("$(sample "$binA")")
    SB+=("$(sample "$binB")")
done
mA="$(median "${SA[@]}")"; mB="$(median "${SB[@]}")"
madA="$(mad "$mA" "${SA[@]}")"; madB="$(mad "$mB" "${SB[@]}")"
awk -v a="$mA" -v b="$mB" -v ma="$madA" -v mb="$madB" -v lbl="$LABEL" -v rounds="$ROUNDS" 'BEGIN{
    d=(b-a); pct=100*d/a;
    band=100*(ma+mb)/a;            # combined-MAD noise band, in % of A
    verdict = (pct < -band) ? "FASTER(candidate)" : (pct > band) ? "SLOWER(candidate)" : "neutral";
    printf "[ab %s] rounds=%d  ref(A)=%.4fs (MAD %.4f)  cand(B)=%.4fs (MAD %.4f)  delta=%+.2f%% band=±%.2f%% => %s\n",
           lbl, rounds, a, ma, b, mb, pct, band, verdict;
}'
