#!/usr/bin/env bash
# Within-producer profiling for `sskm-produce` (the ISOLATED super-k-mer producer). Because there is
# no Phase-2 compaction here, the attribution falls straight onto the producer's own functions
# (add_nucleotide, reverse_complement/canonicalize, phi/minimizer_rank, minimizer_is_ambiguous,
# kmer_lt_kmer, mask_absent_nucleotides, compute_new_candidate_skmer) -- the breakdown the
# whole-construct diag_perf.sh could only aggregate. Emits a flamegraph, a flat perf report, and a
# `perf stat` (IPC + branch/cache misses) to test the "pair shift is branch-heavy" hypothesis noted
# in Skmer.hpp. Uses the Profile build (-O3 -ggdb3 -fno-inline). Monothread by construction.
#
# Requires kernel.perf_event_paranoid <= 1.
set -uo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SSKM="${SSKM:-$ROOT/build-profile/bin/sskm-produce}"
GDIR="$ROOT/benchmark/data/genomes"
OUT="${OUT:-$ROOT/benchmark/results/latest/perf-producer}"
FG="${FG:-$ROOT/thirdparty/FlameGraph}"
FREQ="${FREQ:-1999}"
GENOMES="${GENOMES:-chr21 celegans}"
K="${K:-21}"; M="${M:-11}"
mkdir -p "$OUT"
[[ -x "$SSKM" ]] || { echo "no binary at $SSKM (build sskm-produce in build-profile first)" >&2; exit 1; }
par="$(cat /proc/sys/kernel/perf_event_paranoid)"
(( par > 1 )) && { echo "perf_event_paranoid=$par (need <=1): sudo sysctl -w kernel.perf_event_paranoid=1" >&2; exit 1; }
if [[ ! -d "$FG" ]]; then
  echo "[perf] cloning FlameGraph" >&2
  git clone --depth 1 https://github.com/brendangregg/FlameGraph.git "$FG" || exit 1
fi
SCRATCH="$(mktemp -d /tmp/sklib-perf-prod.XXXXXX)"; trap 'rm -rf "$SCRATCH"' EXIT
for g in $GENOMES; do
  fa="$GDIR/$g.sanitized.fa"; [[ -s "$fa" ]] || { echo "skip $g (missing)" >&2; continue; }
  pd="$SCRATCH/$g.perf.data"
  echo ">>> perf record $g (freq=$FREQ, k=$K m=$M)" >&2
  perf record -F "$FREQ" --call-graph dwarf -o "$pd" -- \
    "$SSKM" -k "$K" -m "$M" -f "$fa" --digest >/dev/null 2>&1 \
    || { echo "  perf record failed for $g" >&2; continue; }
  perf script -i "$pd" 2>/dev/null | "$FG/stackcollapse-perf.pl" > "$OUT/$g.folded" 2>/dev/null
  "$FG/flamegraph.pl" --title "sskm-produce $g k$K m$M" "$OUT/$g.folded" > "$OUT/$g.svg" 2>/dev/null
  perf report -i "$pd" --stdio -g none --percent-limit 0.5 2>/dev/null > "$OUT/$g.flat.txt"
  # Counter-level view: IPC and branch/cache miss rates (the micro-arch signature of the hot loop).
  perf stat -o "$OUT/$g.stat.txt" -- \
    "$SSKM" -k "$K" -m "$M" -f "$fa" --digest >/dev/null 2>&1 \
    || echo "  (perf stat failed for $g)" >&2
  echo "  -> $OUT/$g.{folded,svg,flat.txt,stat.txt} (folded lines: $(wc -l < "$OUT/$g.folded"))" >&2
done
echo "PERF DONE -> $OUT" >&2
