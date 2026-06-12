#!/usr/bin/env bash
# perf flamegraph + numeric function attribution for `sskm construct`, to break down the Phase-1
# producer (FASTA parse / minimizer / bucket routing) vs the Phase-2 compaction, and to spot
# writer/mutex contention at high -t. Uses the Profile build (-O3 -ggdb3 -fno-inline). Run AFTER the
# scaling sweep (perf load would skew wall timings).
#
# Requires kernel.perf_event_paranoid <= 1.
set -uo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SSKM="${SSKM:-$ROOT/build-profile/bin/sskm}"
GDIR="$ROOT/benchmark/data/genomes"
OUT="${OUT:-$ROOT/benchmark/results/latest/perf}"
FG="${FG:-$ROOT/thirdparty/FlameGraph}"
FREQ="${FREQ:-1999}"
GENOMES="${GENOMES:-chr21 celegans}"
THREADS="${THREADS:-1 8}"
K=21; M=11
mkdir -p "$OUT"
[[ -x "$SSKM" ]] || { echo "no binary at $SSKM" >&2; exit 1; }
par="$(cat /proc/sys/kernel/perf_event_paranoid)"
(( par > 1 )) && { echo "perf_event_paranoid=$par (need <=1): sudo sysctl -w kernel.perf_event_paranoid=1" >&2; exit 1; }
if [[ ! -d "$FG" ]]; then
  echo "[perf] cloning FlameGraph" >&2
  git clone --depth 1 https://github.com/brendangregg/FlameGraph.git "$FG" || exit 1
fi
SCRATCH="$(mktemp -d /tmp/sklib-perf.XXXXXX)"; trap 'rm -rf "$SCRATCH"' EXIT
for g in $GENOMES; do
  fa="$GDIR/$g.sanitized.fa"; [[ -s "$fa" ]] || { echo "skip $g (missing)" >&2; continue; }
  for t in $THREADS; do
    tag="${g}_t${t}"; pd="$SCRATCH/$tag.perf.data"
    echo ">>> perf record $tag (freq=$FREQ)" >&2
    perf record -F "$FREQ" --call-graph dwarf -o "$pd" -- \
      "$SSKM" construct -k $K -m $M -f "$fa" -o "$SCRATCH/o.sskm" -t "$t" --tmp-dir "$SCRATCH" >/dev/null 2>&1 \
      || { echo "  perf record failed for $tag" >&2; continue; }
    rm -f "$SCRATCH/o.sskm"
    perf script -i "$pd" 2>/dev/null | "$FG/stackcollapse-perf.pl" > "$OUT/$tag.folded" 2>/dev/null
    "$FG/flamegraph.pl" --title "sskm construct $g -t$t" "$OUT/$tag.folded" > "$OUT/$tag.svg" 2>/dev/null
    perf report -i "$pd" --stdio -g none --percent-limit 0.3 2>/dev/null > "$OUT/$tag.flat.txt"
    echo "  -> $OUT/$tag.{folded,svg,flat.txt} (folded lines: $(wc -l < "$OUT/$tag.folded"))" >&2
  done
done
echo "PERF DONE -> $OUT" >&2
