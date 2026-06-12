#!/usr/bin/env bash
# Construction multithread-scaling diagnostic.
#
# For each (genome, thread count, rep): build the index under /usr/bin/time -v with SKLIB_TIMING=1
# and record wall time, the Phase-1 (serial producer) / Phase-2 (parallel compaction) split, CPU%
# and peak RSS. Phase times come from the `[sklib-timing]` stderr line emitted by build_bucketed;
# wall/CPU%/RSS from GNU time. Output is a tidy long-format CSV (one row per rep).
#
#   bash benchmark/scripts/construct_scaling.sh
#   GENOMES="ecoli chr21" THREADS="1 4 8" REPS=2 bash benchmark/scripts/construct_scaling.sh
#
# Inputs must already be fetched+sanitized (benchmark/scripts/genomes.sh prepare_genome).
set -uo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SSKM="${SSKM:-$ROOT/build-timing/bin/sskm}"
GDIR="${GDIR:-$ROOT/benchmark/data/genomes}"
OUT="${OUT:-$ROOT/benchmark/results/latest/construct_scaling.csv}"
TIME_BIN="${TIME_BIN:-/usr/bin/time}"
GENOMES="${GENOMES:-ecoli yeast chr21 celegans}"
THREADS="${THREADS:-1 2 3 4 6 8 12 16 22}"
K="${K:-21}"; M="${M:-11}"
REPS="${REPS:-3}"

command -v "$TIME_BIN" >/dev/null || { echo "need GNU time at $TIME_BIN" >&2; exit 1; }
[[ -x "$SSKM" ]] || { echo "no binary at $SSKM (build build-timing first)" >&2; exit 1; }

SCRATCH="$(mktemp -d /tmp/sklib-scaling.XXXXXX)"
trap 'rm -rf "$SCRATCH"' EXIT
mkdir -p "$(dirname "$OUT")"
echo "genome,n_bp,k,m,threads,rep,wall_s,phase1_s,phase2_s,cpu_pct,peak_rss_kb" > "$OUT"

bp_of() { awk '!/^>/{n+=length($0)} END{print n+0}' "$1"; }

for g in $GENOMES; do
  fa="$GDIR/$g.sanitized.fa"
  [[ -s "$fa" ]] || { echo "skip $g (missing $fa)" >&2; continue; }
  nbp="$(bp_of "$fa")"
  echo ">>> $g ($nbp bp), k=$K m=$M" >&2
  # One warmup (prime the page cache for this input); result discarded.
  SKLIB_TIMING=1 "$SSKM" construct -k "$K" -m "$M" -f "$fa" -o "$SCRATCH/w.sskm" -t 8 --tmp-dir "$SCRATCH" >/dev/null 2>&1
  rm -f "$SCRATCH/w.sskm"
  for th in $THREADS; do
    for rep in $(seq 1 "$REPS"); do
      log="$SCRATCH/run.log"
      t0="$(date +%s.%N)"
      SKLIB_TIMING=1 "$TIME_BIN" -v "$SSKM" construct -k "$K" -m "$M" -f "$fa" \
          -o "$SCRATCH/o.sskm" -t "$th" --tmp-dir "$SCRATCH" >/dev/null 2>"$log"
      rc=$?
      t1="$(date +%s.%N)"
      if [[ $rc -ne 0 ]]; then echo "  FAIL $g t=$th rep=$rep (rc=$rc)" >&2; sed -n '1,6p' "$log" >&2; continue; fi
      wall="$(awk "BEGIN{printf \"%.4f\", $t1-$t0}")"
      p1="$(awk -F'phase1_s=' '/sklib-timing/{split($2,a," "); print a[1]}' "$log")"
      p2="$(awk -F'phase2_s=' '/sklib-timing/{split($2,a," "); print a[1]}' "$log")"
      cpu="$(awk -F': ' '/Percent of CPU/{gsub(/%/,"",$2); print $2}' "$log")"
      rss="$(awk -F': ' '/Maximum resident set size/{print $2}' "$log")"
      echo "$g,$nbp,$K,$M,$th,$rep,$wall,${p1:-NA},${p2:-NA},${cpu:-NA},${rss:-NA}" >> "$OUT"
      echo "  t=$th rep=$rep wall=${wall}s p1=${p1}s p2=${p2}s cpu=${cpu}% rss=${rss}KB" >&2
      rm -f "$SCRATCH/o.sskm"
    done
  done
done
echo "DONE -> $OUT ($(($(wc -l < "$OUT")-1)) rows)" >&2
