#!/usr/bin/env bash
# Driver for the isolated mono-thread set_union benchmark (union_bench harness).
#
# For each (dataset, k, J) it builds A (sskm construct), B (mutate A -> construct), and a FROZEN
# reference union O_ref (sskm setop --op union -t1, generated once from the current binary), then
# runs the harness in --verify (content-equivalence vs O_ref) and/or --bench (median/min/stddev/MAD)
# modes, pinned to one core. All inputs and O_ref are cached; O_ref is NEVER regenerated once it
# exists (the union k-mer set is invariant across optimizations), so it must be created from a
# pristine build BEFORE any optimization edit.
#
#   bash benchmark/scripts/union_bench.sh [prep|verify|bench|all]   (default: all)
#   DATASETS="chr21" KS="31" JS="0.5" bash benchmark/scripts/union_bench.sh bench
#
# Knobs (env): DATASETS, KS (k list; m=k/2 via the case below), JS (target Jaccard list),
#   SEED, WARMUP, REPS, PIN (taskset core), TAG (label written into the TSV), OUT_TSV.
set -uo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${BUILD:-$ROOT/build-union}"
SSKM="${SSKM:-$BUILD/bin/sskm}"
HARNESS="${HARNESS:-$BUILD/benchmark/union_bench/union_bench}"
GEN="$ROOT/benchmark/data/genomes"
MUT="$ROOT/benchmark/scripts/mutate.py"
S="${SCRATCH:-$ROOT/benchmark/results/union_bench/scratch}"
SHM="${SHM:-/dev/shm}"
mkdir -p "$S"

DATASETS="${DATASETS:-chr21 celegans}"
KS="${KS:-21 31 63}"             # 21->uint32(4B) 31->uint64(8B) 63->__uint128(16B), via m=k/2
JS="${JS:-0.1 0.5 0.9}"
SEED="${SEED:-1234}"
WARMUP="${WARMUP:-2}"
REPS="${REPS:-9}"
PIN="${PIN:-5}"
TAG="${TAG:-cur}"
OP="${OP:-union}"                 # union | intersection | diff | xor (all share materialize_setop)
OUT_TSV="${OUT_TSV:-$ROOT/benchmark/results/union_bench/last.tsv}"
MODE="${1:-all}"

for x in "$SSKM" "$HARNESS"; do
    [[ -x "$x" ]] || { echo "missing binary: $x (build with -DWITH_UNION_BENCH=ON)" >&2; exit 2; }
done
command -v taskset >/dev/null || PIN=""
PINCMD=(); [[ -n "$PIN" ]] && PINCMD=(taskset -c "$PIN")

m_of_k() { echo $(( $1 / 2 )); }                      # m ~ k/2  (31->15, 63->31)
rate_of() { python3 -c "k=$1;J=$2;s=2*J/(1+J);print(f'{1-s**(1/k):.6f}')"; }

log() { printf '[ub] %s\n' "$*" >&2; }

# Build the A/B indexes + frozen O_ref for one (dataset,k,J). Echoes "A B REF" paths.
prepare() {
    local ds="$1" k="$2" J="$3" m fa A mutfa Bf REF rate
    m="$(m_of_k "$k")"
    fa="$GEN/$ds.sanitized.fa"
    [[ -f "$fa" ]] || { log "missing genome $fa"; return 1; }
    A="$S/$ds.k$k.A.sskm"
    mutfa="$S/$ds.k$k.J$J.mut.fa"
    Bf="$S/$ds.k$k.J$J.B.sskm"
    REF="$S/$ds.k$k.J$J.$OP.ref.sskm"
    [[ -s "$A" ]]  || { log "construct A  $ds k=$k m=$m"; "$SSKM" construct -f "$fa" -k "$k" -m "$m" -o "$A" -t "$(nproc)" 2>/dev/null || return 1; }
    if [[ ! -s "$Bf" ]]; then
        rate="$(rate_of "$k" "$J")"
        log "mutate+construct B  $ds k=$k J=$J (rate=$rate)"
        python3 "$MUT" "$fa" "$rate" "$SEED" > "$mutfa" 2>/dev/null || return 1
        "$SSKM" construct -f "$mutfa" -k "$k" -m "$m" -o "$Bf" -t "$(nproc)" 2>/dev/null || return 1
    fi
    # Frozen reference: created once, never regenerated (set is optimization-invariant).
    [[ -s "$REF" ]] || { log "FREEZE O_ref (sskm $OP -t1)  $ds k=$k J=$J"; "$SSKM" setop -a "$A" -b "$Bf" --op "$OP" -o "$REF" -t 1 2>/dev/null || return 1; }
    echo "$A" "$Bf" "$REF"
}

printf 'tag\tdataset\tk\tm\tJ\tstore\tmedian_s\tmin_s\tstddev_s\tmad_s\tmkmer_s\trecords\tverify\n' > "$OUT_TSV"
printf '%-7s %-9s %-3s %-4s %-9s %-9s %-9s %-9s %-9s %-7s %s\n' \
    tag dataset k J median_s min_s stddev_s mad_s Mkmer/s store verify
fail=0
for ds in $DATASETS; do
  for k in $KS; do
    m="$(m_of_k "$k")"
    for J in $JS; do
        read -r A Bf REF < <(prepare "$ds" "$k" "$J") || { log "prepare failed $ds k=$k J=$J"; fail=1; continue; }
        verir="-"
        if [[ "$MODE" == verify || "$MODE" == all ]]; then
            vline="$("${PINCMD[@]}" "$HARNESS" --a "$A" --b "$Bf" --op "$OP" --mode verify --ref "$REF" --out "$SHM/ub_v.sskm" 2>/dev/null)"
            verir="$(awk -F'\t' '$1=="VERIFY"{print $2}' <<<"$vline")"
            [[ "$verir" == PASS ]] || { log "VERIFY $verir  $ds k=$k J=$J :: $vline"; fail=1; }
        fi
        med="-"; mn="-"; sd="-"; mad="-"; mk="-"; rec="-"; store="-"
        if [[ "$MODE" == bench || "$MODE" == all ]]; then
            rline="$("${PINCMD[@]}" "$HARNESS" --a "$A" --b "$Bf" --op "$OP" --mode bench --warmup "$WARMUP" --reps "$REPS" --out "$SHM/ub_b.sskm" 2>/dev/null)"
            eval "$(awk -F'\t' '$1=="RESULT"{for(i=2;i<=NF;i++){split($i,a,"=");printf "%s=%s\n",a[1],a[2]}}' <<<"$rline")"
            med="${median_s:--}"; mn="${min_s:--}"; sd="${stddev_s:--}"; mad="${mad_s:--}"
            mk="${mkmer_s:--}"; rec="${records:--}"; store="${store:--}"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$TAG" "$ds" "$k" "$m" "$J" "$store" "$med" "$mn" "$sd" "$mad" "$mk" "$rec" "$verir" >> "$OUT_TSV"
        printf '%-7s %-9s %-3s %-4s %-9s %-9s %-9s %-9s %-9s %-7s %s\n' \
            "$TAG" "$ds" "$k" "$J" "$med" "$mn" "$sd" "$mad" "$mk" "$store" "$verir"
    done
  done
done
log "wrote $OUT_TSV"
[[ "$fail" == 0 ]] || { log "SOME CONFIGS FAILED (verify or prepare)"; exit 1; }
log "done"
