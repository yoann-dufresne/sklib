#!/usr/bin/env bash
# Relaxed correctness gate for producer optimizations that CHANGE the super-k-mer segmentation.
#
# Once the producer is allowed to emit a DIFFERENT super-k-mer stream (different minimizer order /
# framing), the bit-exact digest no longer applies; the contract becomes: the FINAL set of canonical
# k-mers the index represents must still equal KMC's set for the genome. This script checks that for
# the current `sskm` build: it builds the index (ASCII), re-extracts the k-mers of the emitted
# super-k-mers with KMC, and diffs that set against KMC run directly on the genome. Segmentation-
# agnostic by construction (KMC re-derives the k-mers whatever the super-k-mers are).
#
#   bash benchmark/scripts/producer_setcheck.sh
#   GENOMES="chr21 celegans" K=31 M=15 bash benchmark/scripts/producer_setcheck.sh
set -uo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SSKM="${SSKM:-$ROOT/build-timing/bin/sskm}"
GDIR="${GDIR:-$ROOT/benchmark/data/genomes}"
GENOMES="${GENOMES:-chr21}"
K="${K:-21}"; M="${M:-11}"
command -v kmc >/dev/null      || { echo "need kmc in PATH" >&2; exit 1; }
command -v kmc_tools >/dev/null|| { echo "need kmc_tools in PATH" >&2; exit 1; }
[[ -x "$SSKM" ]] || { echo "no sskm at $SSKM (build it first)" >&2; exit 1; }
seq2fa() { awk '{print ">s" NR; print $0}'; }
TMP="$(mktemp -d /tmp/sklib-setcheck.XXXXXX)"; trap 'rm -rf "$TMP"' EXIT

# Sorted-unique canonical k-mer set of a FASTA via KMC. $1=in $2=fmt(-fa/-fm) $3=k $4=out
kmc_set() {
  local db="$TMP/db_$RANDOM"
  kmc -k"$3" -ci1 "$2" "$1" "$db" "$TMP" >/dev/null 2>&1 || return 1
  kmc_tools transform "$db" dump /dev/stdout 2>/dev/null | cut -f1 | LC_ALL=C sort -u > "$4"
  local rc=${PIPESTATUS[0]}; rm -f "$db".kmc_pre "$db".kmc_suf; return "$rc"
}

fail=0
for g in $GENOMES; do
  fa="$GDIR/$g.sanitized.fa"; [[ -s "$fa" ]] || { echo "skip $g (missing $fa)" >&2; continue; }
  ascii="$TMP/$g.ascii"
  if ! "$SSKM" construct -k "$K" -m "$M" -f "$fa" --ascii -o "$ascii" 2>/dev/null; then
    echo "  $g: construct failed" >&2; fail=1; continue
  fi
  # First field of each ASCII record is the decoded super-k-mer sequence; KMC re-extracts its k-mers.
  tail -n +2 "$ascii" | cut -d' ' -f1 | seq2fa > "$TMP/$g.skfa"
  if ! { kmc_set "$TMP/$g.skfa" -fa "$K" "$TMP/$g.setA" && kmc_set "$fa" -fm "$K" "$TMP/$g.setB"; }; then
    echo "  $g: KMC run failed" >&2; fail=1; continue
  fi
  if diff -q "$TMP/$g.setA" "$TMP/$g.setB" >/dev/null; then
    echo "  $g (k=$K m=$M): list k-mer set == KMC genome set ($(wc -l < "$TMP/$g.setB") k-mers) [OK]"
  else
    a=$(comm -23 "$TMP/$g.setA" "$TMP/$g.setB" | wc -l)
    b=$(comm -13 "$TMP/$g.setA" "$TMP/$g.setB" | wc -l)
    echo "  $g (k=$K m=$M): SET MISMATCH (list-only=$a, genome-only=$b) [FAIL]" >&2; fail=1
  fi
done
if [[ $fail -eq 0 ]]; then echo "SET CHECK OK"; else echo "SET CHECK FAILED" >&2; exit 1; fi
