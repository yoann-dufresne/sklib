#!/usr/bin/env bash
# Median-of-N producer throughput + bit-exact digest gate (the PRODUCER_SPEEDUP loop's measurement
# workhorse). Unlike producer_bench.sh (best-of-REPS single process), this launches the producer as
# N independent processes and reports the MEDIAN Mskmer/s with min/max spread, so a change's delta can
# be judged against run-to-run noise (the task's "temps médian sur N runs" protocol). One untimed
# warmup pass per genome primes the page cache. Every run's stream digest is checked against the
# committed reference (benchmark/results/reference/producer_digest.tsv); any mismatch fails the run.
#
#   K=31 M=15 RUNS=9 bash benchmark/scripts/producer_median.sh
#   GENOMES=chr21 K=63 M=31 RUNS=7 SSKM=build-timing/bin/sskm-produce bash benchmark/scripts/producer_median.sh
set -uo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SSKM="${SSKM:-$ROOT/build-timing/bin/sskm-produce}"
GDIR="${GDIR:-$ROOT/benchmark/data/genomes}"
REF="${REF:-$ROOT/benchmark/results/reference/producer_digest.tsv}"
GENOMES="${GENOMES:-chr21 celegans}"
K="${K:-31}"; M="${M:-15}"
RUNS="${RUNS:-9}"     # recorded process launches per genome
REPS="${REPS:-1}"     # within-process passes; binary reports best-of-REPS per launch

[[ -x "$SSKM" ]] || { echo "no binary at $SSKM (build sskm-produce first)" >&2; exit 1; }

ref_digest_of() {
  awk -v g="$1" -v k="$2" -v m="$3" -F'\t' \
    '!/^#/ && $1==g && $2==k && $3==m {print $6; exit}' "$REF" 2>/dev/null
}

fail=0
printf '%-9s  k=%s m=%s  RUNS=%s (median Mskmer/s [min..max], digest gate)\n' "" "$K" "$M" "$RUNS" >&2
for g in $GENOMES; do
  fa="$GDIR/$g.sanitized.fa"
  [[ -s "$fa" ]] || { echo "skip $g (missing $fa)" >&2; continue; }
  ref="$(ref_digest_of "$g" "$K" "$M")"

  # warmup (prime page cache, not recorded)
  "$SSKM" -k "$K" -m "$M" -f "$fa" --digest --reps "$REPS" >/dev/null 2>&1

  rates=(); digest=""; dig_ok="OK"
  for ((i=0; i<RUNS; i++)); do
    line="$("$SSKM" -k "$K" -m "$M" -f "$fa" --digest --reps "$REPS" 2>/dev/null)"
    r="$(sed -n 's/.* Mskmer_s=\([0-9.]*\).*/\1/p' <<<"$line")"
    d="$(sed -n 's/.* digest=\(0x[0-9a-fA-F]*\).*/\1/p' <<<"$line")"
    [[ -z "$r" ]] && { echo "  FAIL $g (no output)" >&2; fail=1; break; }
    rates+=("$r"); digest="$d"
    if [[ -n "$ref" && "$d" != "$ref" ]]; then dig_ok="MISMATCH"; fail=1; fi
  done
  [[ ${#rates[@]} -eq 0 ]] && continue

  read -r med mn mx <<<"$(printf '%s\n' "${rates[@]}" | sort -n | awk '
    {a[NR]=$1} END{
      n=NR; mn=a[1]; mx=a[n];
      med=(n%2)?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
      printf "%.3f %.3f %.3f", med, mn, mx}')"
  [[ -z "$ref" ]] && dig_ok="NEW(ref=NA)"
  printf '  %-9s  %7.3f  [%7.3f .. %7.3f]  Mskmer/s   digest=%s [%s]\n' \
    "$g" "$med" "$mn" "$mx" "${digest:-?}" "$dig_ok" >&2
done
[[ $fail -ne 0 ]] && { echo "FAILED (digest mismatch or run error)" >&2; exit 1; }
exit 0
