#!/usr/bin/env bash
# Isolated super-k-mer PRODUCER throughput benchmark + bit-exact regression gate.
#
# Runs `sskm-produce` (FileSkmerator only -- no bucketing/sort/writer; see app/produce/) over each
# genome and records throughput (Mskmer/s), wall, peak RSS, and the order-sensitive FNV-1a digest of
# the whole skmer stream. The digest is checked against the committed baseline
# (benchmark/results/reference/producer_digest.tsv): a MISMATCH fails the run, so a producer
# optimization is measured and correctness-verified in one shot. chr21 is the dev signal, celegans
# the final verdict.
#
#   bash benchmark/scripts/producer/producer_bench.sh
#   GENOMES="chr21" REPS=7 SSKM=build-opt/bin/sskm-produce bash benchmark/scripts/producer/producer_bench.sh
#
# Inputs must already be fetched+sanitized (benchmark/scripts/genomes.sh). The binary's own
# best-of-REPS timing is used (warm page cache); GNU time adds peak RSS.
set -uo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
SSKM="${SSKM:-$ROOT/build-timing/bin/sskm-produce}"
GDIR="${GDIR:-$ROOT/benchmark/data/genomes}"
REF="${REF:-$ROOT/benchmark/results/reference/producer_digest.tsv}"
OUT="${OUT:-$ROOT/benchmark/results/latest/producer_bench.csv}"
TIME_BIN="${TIME_BIN:-/usr/bin/time}"
GENOMES="${GENOMES:-chr21 celegans}"
K="${K:-21}"; M="${M:-11}"
REPS="${REPS:-5}"

[[ -x "$SSKM" ]] || { echo "no binary at $SSKM (build sskm-produce first)" >&2; exit 1; }
command -v "$TIME_BIN" >/dev/null || { echo "need GNU time at $TIME_BIN" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
echo "genome,n_bp,k,m,reps,skmers,digest,ref_digest,digest_ok,best_wall_s,mskmer_s,peak_rss_kb" > "$OUT"

bp_of() { awk '!/^>/{n+=length($0)} END{print n+0}' "$1"; }
# Reference digest for (genome,k,m) from the committed baseline, or empty if none.
ref_digest_of() {
  awk -v g="$1" -v k="$2" -v m="$3" -F'\t' \
    '!/^#/ && $1==g && $2==k && $3==m {print $6; exit}' "$REF" 2>/dev/null
}

fail=0
for g in $GENOMES; do
  fa="$GDIR/$g.sanitized.fa"
  [[ -s "$fa" ]] || { echo "skip $g (missing $fa)" >&2; continue; }
  nbp="$(bp_of "$fa")"
  tlog="$(mktemp)"
  # The binary prints one `[producer] ...` line on stdout; GNU time writes -v to stderr.
  line="$("$TIME_BIN" -v "$SSKM" -k "$K" -m "$M" -f "$fa" --digest --reps "$REPS" 2>"$tlog")"
  rc=$?
  if [[ $rc -ne 0 || -z "$line" ]]; then echo "  FAIL $g (rc=$rc)" >&2; cat "$tlog" >&2; rm -f "$tlog"; fail=1; continue; fi

  skmers="$(sed -n 's/.* skmers=\([0-9]*\).*/\1/p' <<<"$line")"
  digest="$(sed -n 's/.* digest=\(0x[0-9a-fA-F]*\).*/\1/p' <<<"$line")"
  wall="$(sed -n 's/.* wall_s=\([0-9.]*\).*/\1/p' <<<"$line")"
  mrate="$(sed -n 's/.* Mskmer_s=\([0-9.]*\).*/\1/p' <<<"$line")"
  rss="$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog")"
  rm -f "$tlog"

  ref="$(ref_digest_of "$g" "$K" "$M")"
  if [[ -z "$ref" ]]; then ok="NEW"; ref="NA"
  elif [[ "$ref" == "$digest" ]]; then ok="OK"
  else ok="MISMATCH"; fail=1; fi

  echo "$g,$nbp,$K,$M,$REPS,$skmers,$digest,$ref,$ok,${wall:-NA},${mrate:-NA},${rss:-NA}" >> "$OUT"
  printf '  %-9s %sx  %8s skmers  %6s Mskmer/s  wall=%ss  rss=%sKB  digest=%s [%s]\n' \
    "$g" "$REPS" "${skmers:-?}" "${mrate:-?}" "${wall:-?}" "${rss:-?}" "${digest:-?}" "$ok" >&2
  [[ "$ok" == "MISMATCH" ]] && echo "    !! digest != reference $ref -- producer output CHANGED" >&2
done

echo "DONE -> $OUT" >&2
if [[ $fail -ne 0 ]]; then
  echo "FAILED: a run errored or a digest diverged from the baseline." >&2
  exit 1
fi
