#!/bin/bash
#
# Cross-validate sskm set operations against KMC on real data.
#
#   tests/setop_verif.sh <A.fa> <B.fa> [k] [m]
#
# Builds sklib lists for A and B (default --buckets, so the quotiented multi-bucket path is
# exercised end to end), then checks, against KMC's k-mer set operations:
#   * cardinalities: intersection_size / union_size / diff_size (both directions) and |A|, |B|;
#   * content: every k-mer KMC puts in a result set is reported present when the corresponding
#     sklib result list is queried, and the result list's k-mer count equals KMC's. Present-all
#     plus equal-count proves set equality (sklib result lists carry no duplicate k-mers).
#
# Requires kmc and kmc_tools on PATH. Override the binaries with SSKM=... / SEQ2FA=... if needed.

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

echo "== KMC ground truth =="
"$kmc" -k"$K" -ci1 -fm "$A_FA" "$WORK/kmcA" "$WORK/kmc_tmp" >/dev/null 2>&1
"$kmc" -k"$K" -ci1 -fm "$B_FA" "$WORK/kmcB" "$WORK/kmc_tmp" >/dev/null 2>&1
"$kmc_tools" simple "$WORK/kmcA" "$WORK/kmcB" intersect      "$WORK/kmcI"   >/dev/null 2>&1
"$kmc_tools" simple "$WORK/kmcA" "$WORK/kmcB" union          "$WORK/kmcU"   >/dev/null 2>&1
"$kmc_tools" simple "$WORK/kmcA" "$WORK/kmcB" kmers_subtract "$WORK/kmcDab" >/dev/null 2>&1
"$kmc_tools" simple "$WORK/kmcB" "$WORK/kmcA" kmers_subtract "$WORK/kmcDba" >/dev/null 2>&1

dump() { "$kmc_tools" transform "$1" dump /dev/stdout | cut -f1 | LC_ALL=C sort; }
dump "$WORK/kmcA"   > "$WORK/A.txt"
dump "$WORK/kmcB"   > "$WORK/B.txt"
dump "$WORK/kmcI"   > "$WORK/I.txt"
dump "$WORK/kmcU"   > "$WORK/U.txt"
dump "$WORK/kmcDab" > "$WORK/Dab.txt"
dump "$WORK/kmcDba" > "$WORK/Dba.txt"
# Symmetric difference A △ B = (A\B) ∪ (B\A). KMC's `simple` has no symdiff op, so merge the two
# subtract dumps — disjoint and already sorted, so `sort -m` is the union.
LC_ALL=C sort -m "$WORK/Dab.txt" "$WORK/Dba.txt" > "$WORK/Sym.txt"

check_size() {  # label  sklib_value  truth_file
    local truth; truth=$(wc -l < "$3")
    if [ "$2" -eq "$truth" ]; then printf "  OK   %-14s %s\n" "$1" "$2"
    else printf "  FAIL %-14s sklib=%s kmc=%s\n" "$1" "$2" "$truth"; fail=1; fi
}

echo "== cardinalities: sklib vs KMC =="
check_size "intersection" "$("$SSKM" setop --op intersection_size -a "$WORK/A.sskm" -b "$WORK/B.sskm")" "$WORK/I.txt"
check_size "union"        "$("$SSKM" setop --op union_size        -a "$WORK/A.sskm" -b "$WORK/B.sskm")" "$WORK/U.txt"
check_size "diff_AB"      "$("$SSKM" setop --op diff_size         -a "$WORK/A.sskm" -b "$WORK/B.sskm")" "$WORK/Dab.txt"
check_size "diff_BA"      "$("$SSKM" setop --op diff_size         -a "$WORK/B.sskm" -b "$WORK/A.sskm")" "$WORK/Dba.txt"
check_size "xor"          "$("$SSKM" setop --op xor_size          -a "$WORK/A.sskm" -b "$WORK/B.sskm")" "$WORK/Sym.txt"
check_size "|A|"          "$("$SSKM" setop --op intersection_size -a "$WORK/A.sskm" -b "$WORK/A.sskm")" "$WORK/A.txt"
check_size "|B|"          "$("$SSKM" setop --op intersection_size -a "$WORK/B.sskm" -b "$WORK/B.sskm")" "$WORK/B.txt"

# Materialize a result, then prove it equals KMC's set: |result| == |truth| and every truth k-mer is
# reported present by querying the result list.
check_content() {  # label  op  list_a  list_b  truth_file
    local label="$1" op="$2" la="$3" lb="$4" truth="$5"
    local out="$WORK/res_$label.sskm"
    local n got total absent msg
    n=$(wc -l < "$truth")
    msg=$("$SSKM" setop --op "$op" -a "$la" -b "$lb" -o "$out" 2>&1 >/dev/null || true)
    got=$(printf '%s' "$msg" | grep -oE '[0-9]+ k-mers' | grep -oE '^[0-9]+' || echo "?")
    [ "$got" = "$n" ] || { printf "  FAIL %-14s |result|=%s kmc=%s\n" "$label" "$got" "$n"; fail=1; }
    if [ "$n" -eq 0 ]; then printf "  OK   %-14s (empty)\n" "$label"; return; fi
    "$SEQ2FA" < "$truth" > "$WORK/$label.fa"
    "$SSKM" query -l "$out" -i "$WORK/$label.fa" -t 1 >"$WORK/$label.q" 2>/dev/null
    total=$(wc -l < "$WORK/$label.q")
    absent=$(grep -c '0' "$WORK/$label.q" || true)
    if [ "$total" -eq "$n" ] && [ "$absent" -eq 0 ] && [ "$got" = "$n" ]; then
        printf "  OK   %-14s %s/%s present, |result|=%s\n" "$label" "$n" "$n" "$got"
    else
        printf "  FAIL %-14s queried=%s truth=%s absent=%s\n" "$label" "$total" "$n" "$absent"; fail=1
    fi
}

echo "== content: KMC result k-mers are present in the sklib result list =="
check_content "intersection" "intersection" "$WORK/A.sskm" "$WORK/B.sskm" "$WORK/I.txt"
check_content "union"        "union"        "$WORK/A.sskm" "$WORK/B.sskm" "$WORK/U.txt"
check_content "diff_AB"      "diff"         "$WORK/A.sskm" "$WORK/B.sskm" "$WORK/Dab.txt"
check_content "diff_BA"      "diff"         "$WORK/B.sskm" "$WORK/A.sskm" "$WORK/Dba.txt"
check_content "xor"          "xor"          "$WORK/A.sskm" "$WORK/B.sskm" "$WORK/Sym.txt"

if [ "$fail" -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; exit 1; fi
