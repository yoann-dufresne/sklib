#!/bin/bash
#
# Cross-validate sskm COMBINED set operations against KMC on real data.
#
#   tests/setop_multi_verif.sh <A.fa> <B.fa> [k] [m]
#
# Unlike setop_verif.sh (one sskm call per operation), this exercises the single-pass combined path:
# ONE `sskm setop` invocation produces all four result lists (A∩B, A∪B, A\B, B\A) plus --sizes. We
# then check, against KMC's k-mer set operations:
#   * cardinalities, from BOTH the --sizes stdout block AND the per-file "wrote N" stderr lines;
#   * content: every k-mer KMC puts in a result set is reported present when the corresponding sklib
#     result list is queried, and the result list's k-mer count equals KMC's (present-all + equal
#     count proves set equality).
#
# Requires kmc and kmc_tools on PATH. Override binaries with SSKM=... / SEQ2FA=... if needed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SSKM="${SSKM:-$REPO_ROOT/build/bin/sskm}"
SEQ2FA="${SEQ2FA:-$SCRIPT_DIR/sequences_2_fa.sh}"

kmc=$(command -v kmc)             || { echo "kmc not found on PATH"; exit 1; }
kmc_tools=$(command -v kmc_tools) || { echo "kmc_tools not found on PATH"; exit 1; }
[ -x "$SSKM" ]   || { echo "sskm not found/executable at $SSKM (build it, or set SSKM=...)"; exit 1; }
[ -f "$SEQ2FA" ] || { echo "sequences_2_fa.sh not found at $SEQ2FA (set SEQ2FA=...)"; exit 1; }

if [ $# -lt 2 ]; then echo "Usage: $0 <A.fa> <B.fa> [k] [m]"; exit 1; fi
A_FA="$1"; B_FA="$2"; K="${3:-31}"; M="${4:-21}"
[ "$M" -lt "$K" ] || { echo "Error: m ($M) must be < k ($K)"; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/kmc_tmp"
fail=0

echo "== build sklib lists (k=$K m=$M, default buckets) =="
"$SSKM" construct -f "$A_FA" -k "$K" -m "$M" -o "$WORK/A.sskm"
"$SSKM" construct -f "$B_FA" -k "$K" -m "$M" -o "$WORK/B.sskm"

echo "== KMC ground truth (all four relations) =="
"$kmc" -k"$K" -ci1 -fm "$A_FA" "$WORK/kmcA" "$WORK/kmc_tmp" >/dev/null 2>&1
"$kmc" -k"$K" -ci1 -fm "$B_FA" "$WORK/kmcB" "$WORK/kmc_tmp" >/dev/null 2>&1
# KMC also supports several ops in one `simple` pass — the head-to-head for the combined sskm path.
"$kmc_tools" simple "$WORK/kmcA" "$WORK/kmcB" \
    intersect "$WORK/kmcI" union "$WORK/kmcU" kmers_subtract "$WORK/kmcDab" >/dev/null 2>&1
"$kmc_tools" simple "$WORK/kmcB" "$WORK/kmcA" kmers_subtract "$WORK/kmcDba" >/dev/null 2>&1

dump() { "$kmc_tools" transform "$1" dump /dev/stdout | cut -f1 | LC_ALL=C sort; }
dump "$WORK/kmcI"   > "$WORK/I.txt"
dump "$WORK/kmcU"   > "$WORK/U.txt"
dump "$WORK/kmcDab" > "$WORK/Dab.txt"
dump "$WORK/kmcDba" > "$WORK/Dba.txt"

echo "== ONE combined sskm call: all four results + --sizes =="
# stdout = key<TAB>value sizes block; stderr = "<rel>: wrote N k-mers to <path>" per output.
"$SSKM" setop -a "$WORK/A.sskm" -b "$WORK/B.sskm" \
    --inter-out "$WORK/resI.sskm" --union-out "$WORK/resU.sskm" \
    --diff-ab-out "$WORK/resDab.sskm" --diff-ba-out "$WORK/resDba.sskm" \
    --sizes >"$WORK/sizes.txt" 2>"$WORK/wrote.txt"

size_of()  { awk -v k="$1" '$1==k{print $2}' "$WORK/sizes.txt"; }     # from --sizes stdout
wrote_of() { awk -v k="$1" '$1==k":"{print $3}' "$WORK/wrote.txt"; }  # from "rel: wrote N ..." stderr

check_size() {  # label  sklib_value  truth_file
    local truth; truth=$(wc -l < "$3")
    if [ "$2" = "$truth" ]; then printf "  OK   %-14s %s\n" "$1" "$2"
    else printf "  FAIL %-14s sklib=%s kmc=%s\n" "$1" "$2" "$truth"; fail=1; fi
}

echo "== cardinalities: combined --sizes vs KMC =="
check_size "intersection" "$(size_of intersection)" "$WORK/I.txt"
check_size "union"        "$(size_of union)"        "$WORK/U.txt"
check_size "diff_ab"      "$(size_of diff_ab)"      "$WORK/Dab.txt"
check_size "diff_ba"      "$(size_of diff_ba)"      "$WORK/Dba.txt"

echo "== cardinalities: combined per-file 'wrote N' vs KMC =="
check_size "wrote_inter"  "$(wrote_of intersection)" "$WORK/I.txt"
check_size "wrote_union"  "$(wrote_of union)"         "$WORK/U.txt"
check_size "wrote_dab"    "$(wrote_of diff_ab)"       "$WORK/Dab.txt"
check_size "wrote_dba"    "$(wrote_of diff_ba)"       "$WORK/Dba.txt"

# Content: materialized result equals KMC's set (|result| == |truth| and every truth k-mer present).
check_content() {  # label  result_list  truth_file
    local label="$1" out="$2" truth="$3"
    local n total absent
    n=$(wc -l < "$truth")
    if [ "$n" -eq 0 ]; then printf "  OK   %-14s (empty)\n" "$label"; return; fi
    "$SEQ2FA" < "$truth" > "$WORK/$label.fa"
    "$SSKM" query -l "$out" -i "$WORK/$label.fa" -t 1 >"$WORK/$label.q" 2>/dev/null
    total=$(wc -l < "$WORK/$label.q")
    absent=$(grep -c '0' "$WORK/$label.q" || true)
    if [ "$total" -eq "$n" ] && [ "$absent" -eq 0 ]; then
        printf "  OK   %-14s %s/%s present\n" "$label" "$n" "$n"
    else
        printf "  FAIL %-14s queried=%s truth=%s absent=%s\n" "$label" "$total" "$n" "$absent"; fail=1
    fi
}

echo "== content: KMC result k-mers are present in the combined sklib result lists =="
check_content "intersection" "$WORK/resI.sskm"   "$WORK/I.txt"
check_content "union"        "$WORK/resU.sskm"   "$WORK/U.txt"
check_content "diff_ab"      "$WORK/resDab.sskm" "$WORK/Dab.txt"
check_content "diff_ba"      "$WORK/resDba.sskm" "$WORK/Dba.txt"

if [ "$fail" -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; exit 1; fi
