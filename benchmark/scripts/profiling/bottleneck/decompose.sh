#!/usr/bin/env bash
# Macro decomposition of sklib set-op time: _size (merge only) vs materialize.
# delta = materialize - size = collect(get_skmer_of_kmer) + recompact + write.
set -uo pipefail
WS=/tmp/sklib-bottleneck
SSKM=$WS/builds/rel/bin/sskm
L=$WS/data/lists
PIN="taskset -c 0"
OUT=$WS/results/decompose.csv
echo "pair,kA,kB,op,result_kmers,t_size,t_mat,delta,frac_post" > "$OUT"

# median wall seconds over $1 reps of the remaining args (warm cache: one warmup first)
med() { local reps="$1"; shift
  $PIN "$@" >/dev/null 2>&1   # warmup
  local t=()
  for ((i=0;i<reps;i++)); do
    local t0 t1; t0=$(date +%s.%N); $PIN "$@" >/dev/null 2>&1; t1=$(date +%s.%N)
    t+=("$(awk "BEGIN{printf \"%.4f\",$t1-$t0}")")
  done
  printf '%s\n' "${t[@]}" | sort -g | awk '{a[NR]=$1} END{print (NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}'
}

size_of() { "$SSKM" setop --op "$1" -a "$2" -b "$3" 2>/dev/null; }

run_pair() { # label A B reps
  local label="$1" A="$2" B="$3" reps="$4"
  [ -f "$A" ] && [ -f "$B" ] || { echo "[skip] $label (missing list)"; return; }
  local kA kB; kA=$(size_of intersection_size "$A" "$A"); kB=$(size_of intersection_size "$B" "$B")
  declare -A RES SOP MOP
  RES[inter]=$(size_of intersection_size "$A" "$B")
  RES[union]=$(size_of union_size        "$A" "$B")
  RES[diff]=$( size_of diff_size          "$A" "$B")
  SOP[inter]=intersection_size; SOP[union]=union_size; SOP[diff]=diff_size
  MOP[inter]=intersection;      MOP[union]=union;      MOP[diff]=diff
  local tmp="$WS/results/_out.sskm"
  for op in inter union diff; do
    local ts tm
    ts=$(med "$reps" "$SSKM" setop --op "${SOP[$op]}" -a "$A" -b "$B")
    tm=$(med "$reps" "$SSKM" setop --op "${MOP[$op]}" -a "$A" -b "$B" -o "$tmp")
    awk -v p="$label" -v ka="$kA" -v kb="$kB" -v op="$op" -v r="${RES[$op]}" \
        -v ts="$ts" -v tm="$tm" 'BEGIN{d=tm-ts; printf "%s,%s,%s,%s,%s,%.4f,%.4f,%.4f,%.3f\n",p,ka,kb,op,r,ts,tm,d,(tm>0?d/tm:0)}' | tee -a "$OUT"
  done
  rm -f "$tmp"
}

# --- scaling at ~constant overlap (self vs 1% mutant) ---
run_pair "scale_ecoli_4.5M"  "$L/ecoliK12.k21m11.sskm"  "$L/ecoliK12_mut01.k21m11.sskm" 7
run_pair "scale_yeast_11.5M" "$L/yeast.k21m11.sskm"     "$L/yeast_mut1.k21m11.sskm"     5
run_pair "scale_chr21_32.7M" "$L/chr21.k21m11.sskm"     "$L/chr21_mut1.k21m11.sskm"     4
[ -f "$L/celegans.k21m11.sskm" ] && run_pair "scale_celegans_91M" "$L/celegans.k21m11.sskm" "$L/celegans_mut1.k21m11.sskm" 3

# --- overlap sweep at constant size (ecoli) ---
run_pair "ov_self_100pct"  "$L/ecoliK12.k21m11.sskm" "$L/ecoliK12.k21m11.sskm"        7
run_pair "ov_mut005"       "$L/ecoliK12.k21m11.sskm" "$L/ecoliK12_mut005.k21m11.sskm" 7
run_pair "ov_mut01"        "$L/ecoliK12.k21m11.sskm" "$L/ecoliK12_mut01.k21m11.sskm"  7
run_pair "ov_mut05"        "$L/ecoliK12.k21m11.sskm" "$L/ecoliK12_mut05.k21m11.sskm"  7
run_pair "ov_real_Sakai"   "$L/ecoliK12.k21m11.sskm" "$L/ecoliSakai.k21m11.sskm"      7

# --- real low-overlap human pair ---
run_pair "real_chr21_chr22" "$L/chr21.k21m11.sskm" "$L/chr22.k21m11.sskm" 4
echo "DONE decompose"
