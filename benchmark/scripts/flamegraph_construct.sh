#!/usr/bin/env bash
# Build sskm in Profile mode, run `sskm construct` on a FASTA under perf,
# and emit an SVG flame graph. Configurable via env vars (see below).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

K="${K:-21}"
M="${M:-11}"
INPUT="${INPUT:-benchmark/data/genomes/ecoli.sanitized.fa}"
FREQ="${FREQ:-997}"
BUILD_DIR="${BUILD_DIR:-build-profile}"
OUT_DIR="${OUT_DIR:-benchmark/results/latest}"
FLAMEGRAPH_DIR="${FLAMEGRAPH_DIR:-thirdparty/FlameGraph}"

if [[ ! -f "$INPUT" ]]; then
    echo "error: input FASTA not found: $INPUT" >&2
    echo "       fetch a catalogued genome first, e.g.: bash benchmark/scripts/fetch_genomes.sh ecoli" >&2
    exit 1
fi

if ! command -v perf >/dev/null 2>&1; then
    echo "error: perf not on PATH (try: sudo apt install linux-tools-generic)" >&2
    exit 1
fi

paranoid="$(cat /proc/sys/kernel/perf_event_paranoid)"
if (( paranoid > 1 )); then
    echo "error: kernel.perf_event_paranoid=$paranoid (need <=1 for user-space call graphs)" >&2
    echo "       sudo sysctl -w kernel.perf_event_paranoid=1" >&2
    exit 1
fi

if [[ ! -d "$FLAMEGRAPH_DIR" ]]; then
    echo "[flamegraph] cloning FlameGraph into $FLAMEGRAPH_DIR"
    git clone --depth 1 https://github.com/brendangregg/FlameGraph.git "$FLAMEGRAPH_DIR"
fi

echo "[flamegraph] configuring + building sskm (Profile) in $BUILD_DIR"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Profile -DWITH_TESTS=OFF >/dev/null
cmake --build "$BUILD_DIR" -j --target sskm

mkdir -p "$OUT_DIR"
PERF_DATA="$OUT_DIR/perf.data"
SSKM_OUT="$OUT_DIR/$(basename "${INPUT%.*}").sskm"
SVG_OUT="$OUT_DIR/flamegraph.svg"

echo "[flamegraph] perf record (k=$K m=$M, freq=$FREQ Hz, input=$INPUT)"
perf record -F "$FREQ" --call-graph dwarf -o "$PERF_DATA" -- \
    "$BUILD_DIR/bin/sskm" construct -k "$K" -m "$M" \
    -f "$INPUT" -o "$SSKM_OUT"

echo "[flamegraph] generating $SVG_OUT"
perf script -i "$PERF_DATA" \
  | "$FLAMEGRAPH_DIR/stackcollapse-perf.pl" \
  | "$FLAMEGRAPH_DIR/flamegraph.pl" \
        --title "sskm construct k=$K m=$M $(basename "$INPUT")" \
  > "$SVG_OUT"

echo
echo "flame graph: $REPO_ROOT/$SVG_OUT"
