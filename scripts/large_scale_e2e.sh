#!/usr/bin/env bash
# End-to-end, full-scale correctness + performance test for sklib's `sskm`.
#
# For each genome x (k,m) it runs five checks, all against an independent
# oracle (KMC) where applicable:
#   1. construction vs KMC      : the k-mer set of the sorted list == the k-mer
#                                 set of the genome (set equality).
#   2. self-query completeness  : querying the genome against its own list yields
#                                 only 1s (no false negatives).
#   3. mixed-query vs KMC        : a query set of {present, random} k-mers; every
#                                 present k-mer is 1, and every k-mer's status
#                                 matches KMC's ground truth (false-positive check).
#   4. serialization + determinism : construct twice -> byte-identical binary and
#                                 ASCII; binary/ASCII/header k-mer counts agree.
#   5. perf & resources         : wall time, peak RSS, list size, throughput.
#
# sklib silently encodes non-ACGT (e.g. N) as G while KMC drops N-containing
# k-mers, so every genome is first "sanitized": uppercased and split at non-ACGT
# runs into pure-ACGT records, fed identically to both tools.
#
# Usage:
#   bash scripts/large_scale_e2e.sh
#   GENOMES="ecoli yeast chr21 chr1" KM="21,11 31,13 32,17" bash scripts/large_scale_e2e.sh
#
# Env knobs (defaults): GENOMES="ecoli chr21"  KM="21,11 31,13"
#   N_PRESENT=5000  N_RANDOM=5000  SEED=1234  KEEP=0  (KEEP=1 keeps scratch dirs)
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

GENOMES="${GENOMES:-ecoli chr21}"
KM="${KM:-21,11 31,13}"
N_PRESENT="${N_PRESENT:-1500}"
N_RANDOM="${N_RANDOM:-1500}"
# Query is serial unless the binary links Intel TBB (std::execution::par falls back
# to serial otherwise), so it is the slow path. The self-query completeness check
# therefore runs on a bounded, structure-preserving prefix of the genome rather than
# the whole thing -- construction-vs-KMC already proves full k-mer-set coverage, so
# the self-query only needs to exercise the multi-k-mer query path. Raise for more.
SELFQUERY_BP="${SELFQUERY_BP:-10000}"
SEED="${SEED:-1234}"
KEEP="${KEEP:-0}"

REL_BIN="${REL_BIN:-$REPO_ROOT/build/bin/sskm}"
DBG_BIN="${DBG_BIN:-$REPO_ROOT/build-debug/bin/sskm}"
HELPER="$SCRIPT_DIR/e2e_helpers.py"
# tests/sequences_2_fa.sh is a bash `while read` loop -- far too slow for the ~30M
# skmer lines a human chromosome produces -- so we wrap lines into FASTA with awk.
seq2fa() { awk '{print ">s" NR; print $0}'; }

OUT="$REPO_ROOT/scripts/out/e2e"
GEN_DIR="$OUT/genomes"
KMC_TMP="$OUT/kmc_tmp"
mkdir -p "$GEN_DIR" "$KMC_TMP"

# ---- pretty logging + pass/fail accounting ----
if [[ -t 1 ]]; then C_G=$'\033[32m'; C_R=$'\033[31m'; C_B=$'\033[1m'; C_0=$'\033[0m'; else C_G=; C_R=; C_B=; C_0=; fi
PASS=0; FAIL=0
declare -a FAILED
CTX="-"
PERF_TSV="$OUT/perf_report.tsv"
: > "$PERF_TSV"

log()  { printf '%s[e2e]%s %s\n' "$C_B" "$C_0" "$*" >&2; }
ok()   { PASS=$((PASS+1)); printf '  %sPASS%s %s\n' "$C_G" "$C_0" "$1" >&2; }
bad()  { FAIL=$((FAIL+1)); printf '  %sFAIL%s %s\n' "$C_R" "$C_0" "$1" >&2; FAILED+=("[$CTX] $1"); }
die()  { printf '%sFATAL%s %s\n' "$C_R" "$C_0" "$*" >&2; exit 2; }

for tool in kmc kmc_tools /usr/bin/time python3 curl gzip cmp; do
    command -v "$tool" >/dev/null 2>&1 || die "required tool not found: $tool"
done
[[ -x "$REL_BIN" ]] || die "release binary missing: $REL_BIN (build it first)"
[[ -x "$DBG_BIN" ]] || die "debug binary missing: $DBG_BIN (build it first)"

# ---- genome catalogue ----
genome_url() {
    case "$1" in
        yeast) echo "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz" ;;
        chr21) echo "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz" ;;
        chr20) echo "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz" ;;
        chr1)  echo "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz" ;;
        *)     echo "" ;;
    esac
}
# Local genome fixtures (no download).
genome_local() {
    case "$1" in
        ecoli)    echo "$REPO_ROOT/data/ecoli.fa" ;;
        sarscov2) echo "$REPO_ROOT/data/sarscov2.fa" ;;
        *)        echo "" ;;
    esac
}
# Construct under the DEBUG build (asserts on) for small genomes, Release for big
# ones. Query is the slow, assert-light path, so it always uses the Release build.
construct_bin() {
    case "$1" in
        chr1|chr20|chr21) echo "$REL_BIN" ;;
        *)                echo "$DBG_BIN" ;;
    esac
}
query_bin() { echo "$REL_BIN"; }

# Resolve a genome name to a sanitized FASTA path (download + decompress + split).
# Sets global SAN.
prepare_genome() {
    local name="$1" raw san url
    san="$GEN_DIR/$name.sanitized.fa"
    SAN="$san"
    if [[ -f "$san" ]]; then log "$name: using cached sanitized FASTA"; return 0; fi

    local local_path
    local_path="$(genome_local "$name")"
    if [[ -n "$local_path" ]]; then
        raw="$local_path"
        [[ -f "$raw" ]] || { log "$name: local fixture missing: $raw"; return 1; }
    else
        url="$(genome_url "$name")"
        [[ -n "$url" ]] || { log "$name: unknown genome (no URL)"; return 1; }
        raw="$GEN_DIR/$name.fa"
        if [[ ! -f "$raw" ]]; then
            log "$name: downloading $url"
            curl -fSL --retry 3 "$url" -o "$GEN_DIR/$name.fa.gz" || { log "$name: download failed"; return 1; }
            gzip -dc "$GEN_DIR/$name.fa.gz" > "$raw" || { log "$name: gunzip failed"; return 1; }
        fi
    fi
    log "$name: sanitizing (uppercase + split at non-ACGT runs)"
    python3 "$HELPER" sanitize "$raw" > "$san.tmp" && mv "$san.tmp" "$san" || { log "$name: sanitize failed"; return 1; }
    return 0
}

# Build a sorted-unique list of canonical k-mers from a FASTA via KMC.
# $1=input fasta  $2=format flag (-fa/-fm)  $3=k  $4=output file
kmc_canon_set() {
    local in="$1" fmt="$2" kk="$3" out="$4"
    local db="$KMC_TMP/db_$$_$RANDOM"
    kmc -k"$kk" -ci1 "$fmt" "$in" "$db" "$KMC_TMP" >/dev/null 2>&1 || return 1
    kmc_tools transform "$db" dump /dev/stdout 2>/dev/null | cut -f1 | LC_ALL=C sort -u > "$out"
    local rc=${PIPESTATUS[0]}
    rm -f "$db".kmc_pre "$db".kmc_suf
    return "$rc"
}

# Run a command under /usr/bin/time -v; capture wall seconds + peak RSS (KB).
# Sets RUN_SEC, RUN_RSS_KB. Returns the command's exit status.
RUN_SEC=0; RUN_RSS_KB=0
run_timed() {
    local tlog t0 t1 status
    tlog="$(mktemp)"
    t0=$(date +%s.%N)
    /usr/bin/time -v "$@" >/dev/null 2>"$tlog"
    status=$?
    t1=$(date +%s.%N)
    RUN_SEC=$(awk "BEGIN{printf \"%.2f\", $t1-$t0}")
    RUN_RSS_KB=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog")
    [[ -n "$RUN_RSS_KB" ]] || RUN_RSS_KB=0
    rm -f "$tlog"
    return "$status"
}

human_mb() { awk "BEGIN{printf \"%.0f\", $1/1024}"; }   # KB -> MB
rate()     { awk "BEGIN{ if($2>0) printf \"%.2f\", $1/$2/1e6; else print \"-\" }"; }  # count/sec -> M/s

# ====================================================================
#  Per-(genome,k,m) run
# ====================================================================
run_one() {
    local genome="$1" k="$2" m="$3" san="$4" cbin="$5" qbin="$6" bp="$7"
    CTX="$genome k=$k m=$m"
    log "=== $CTX  (construct: $([[ "$cbin" == "$DBG_BIN" ]] && echo DEBUG || echo Release), query: Release) ==="
    local work="$OUT/work/${genome}_${k}_${m}"
    rm -rf "$work"; mkdir -p "$work"

    local sskm="$work/list.sskm" ascii="$work/list.ascii"

    # ---- step 0 + check 5: construct (timed) ----
    if ! run_timed "$cbin" construct -k "$k" -m "$m" -f "$san" -o "$sskm"; then
        bad "construct failed"; return
    fi
    local c_sec="$RUN_SEC" c_rss="$RUN_RSS_KB"
    "$cbin" construct -k "$k" -m "$m" -f "$san" --ascii -o "$ascii" 2>/dev/null || { bad "ascii construct failed"; return; }

    local nskmers; nskmers=$(python3 "$HELPER" bincount "$sskm" 2>/dev/null)
    [[ -n "$nskmers" ]] || { bad "could not read binary header"; return; }

    # ---- check 1: construction vs KMC (set equality) ----
    local setA="$work/setA_list.txt" setB="$work/setB_genome.txt" skfa="$work/skmers.fa"
    tail -n +2 "$ascii" | cut -d' ' -f1 | seq2fa > "$skfa"
    if kmc_canon_set "$skfa" -fa "$k" "$setA" && kmc_canon_set "$san" -fm "$k" "$setB"; then
        if diff -q "$setA" "$setB" >/dev/null; then
            ok "construction vs KMC: list k-mer set == genome k-mer set ($(wc -l < "$setB") k-mers)"
        else
            local onlyA onlyB
            onlyA=$(comm -23 "$setA" "$setB" | wc -l)
            onlyB=$(comm -13 "$setA" "$setB" | wc -l)
            bad "construction vs KMC: set mismatch (list-only=$onlyA, genome-only=$onlyB)"
        fi
    else
        bad "construction vs KMC: KMC run failed"
    fi

    # ---- check 2: self-query completeness on a bounded prefix sample (timed) ----
    local hits="$work/self_hits.txt" selfq="$work/selfquery.fa"
    python3 "$HELPER" prefix "$san" "$SELFQUERY_BP" > "$selfq"
    local self_bp; self_bp=$(grep -v '^>' "$selfq" | tr -d '\n' | wc -c)
    if run_timed "$qbin" query -l "$sskm" -i "$selfq" -o "$hits"; then
        local q_sec="$RUN_SEC" q_rss="$RUN_RSS_KB"
        if [[ -s "$hits" ]] && ! grep -q '0' "$hits"; then
            ok "self-query completeness: all $(wc -l < "$hits") skmer rows over ${self_bp}bp sample are 1 (no false negatives)"
        else
            local nzero; nzero=$(grep -o '0' "$hits" | wc -l)
            bad "self-query completeness: $nzero zero(s) found over ${self_bp}bp sample (false negatives)"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$genome" "$k" "$m" "$bp" "$nskmers" \
            "$c_sec" "$(human_mb "$c_rss")" "$(rate "$bp" "$c_sec")" \
            "$q_sec" "$(rate "$self_bp" "$q_sec")" >> "$PERF_TSV"
    else
        bad "self-query failed"
    fi

    # ---- check 3: mixed query (present + random) vs KMC truth ----
    local present="$work/present.txt" rand="$work/random.txt"
    local qk="$work/qkmers.txt" qfa="$work/Q.fa" qhits="$work/qhits.txt"
    local qcanon="$work/qcanon.txt" qcanon_su="$work/qcanon_su.txt"
    local present_canon="$work/present_canon.txt" truth="$work/truth.txt"
    if [[ -s "$setB" ]]; then
        local total np
        total=$(wc -l < "$setB")
        awk -v tot="$total" -v mm="$N_PRESENT" \
            'BEGIN{ if(mm>tot)mm=tot; if(mm<1)mm=1; step=tot/mm; for(i=0;i<mm;i++){t=int(i*step)+1; want[t]=1} }
             (NR in want){ print }' "$setB" > "$present"
        np=$(wc -l < "$present")
        python3 "$HELPER" randkmers "$k" "$N_RANDOM" "$SEED" > "$rand"
        cat "$present" "$rand" > "$qk"
        seq2fa < "$qk" > "$qfa"

        if "$qbin" query -l "$sskm" -i "$qfa" -o "$qhits" 2>/dev/null; then
            local nq nrows
            nq=$(wc -l < "$qk"); nrows=$(wc -l < "$qhits")
            if [[ "$nq" -ne "$nrows" ]]; then
                bad "mixed query: row count $nrows != query count $nq (1:1 mapping broken)"
            else
                # present subset must all be 1
                if head -n "$np" "$qhits" | grep -q '0'; then
                    bad "mixed query: a present k-mer returned 0 (false negative)"
                else
                    ok "mixed query: all $np present k-mers returned 1"
                fi
                # full status vs KMC ground truth
                python3 "$HELPER" canon "$qk" > "$qcanon"
                LC_ALL=C sort -u "$qcanon" > "$qcanon_su"
                LC_ALL=C comm -12 "$qcanon_su" "$setB" > "$present_canon"
                python3 "$HELPER" truth "$present_canon" "$qcanon" > "$truth"
                if diff -q "$truth" "$qhits" >/dev/null; then
                    local nfound; nfound=$(grep -c '^1$' "$truth")
                    ok "mixed query vs KMC: all $nq statuses match ($nfound present / $((nq-nfound)) absent)"
                else
                    local nmis; nmis=$(paste "$truth" "$qhits" | awk '$1!=$2{c++} END{print c+0}')
                    bad "mixed query vs KMC: $nmis/$nq status mismatches"
                fi
            fi
        else
            bad "mixed query failed"
        fi
    else
        bad "mixed query: genome k-mer set empty"
    fi

    # ---- check 4: serialization + determinism ----
    local sskm2="$work/list2.sskm" ascii2="$work/list2.ascii"
    "$cbin" construct -k "$k" -m "$m" -f "$san" -o "$sskm2" 2>/dev/null
    "$cbin" construct -k "$k" -m "$m" -f "$san" --ascii -o "$ascii2" 2>/dev/null
    local det_ok=1
    cmp -s "$sskm" "$sskm2" || { det_ok=0; }
    cmp -s "$ascii" "$ascii2" || { det_ok=0; }
    if [[ "$det_ok" -eq 1 ]]; then
        ok "determinism: construct x2 byte-identical (binary + ASCII)"
    else
        bad "determinism: re-construct differs (non-deterministic output)"
    fi
    # header / count consistency
    local acount fsize sz
    acount=$(head -1 "$ascii" | awk '{print $3}')
    fsize=$(stat -c%s "$sskm")
    local bhdr; bhdr=$(python3 "$HELPER" binheadersize "$sskm" 2>/dev/null || echo 32)
    if [[ "$nskmers" == "$acount" ]] && awk "BEGIN{exit !(($fsize-$bhdr)%$nskmers==0)}"; then
        sz=$(( (fsize-bhdr)/nskmers ))
        ok "serialization: binary count == ASCII count == $nskmers (sizeof Skmer=${sz}B, header=${bhdr}B)"
    else
        bad "serialization: count mismatch (binary=$nskmers ascii=$acount filesize=$fsize)"
    fi

    [[ "$KEEP" == "1" ]] || rm -rf "$work"
}

# ====================================================================
#  Main
# ====================================================================
log "genomes: $GENOMES   (k,m): $KM   present=$N_PRESENT random=$N_RANDOM seed=$SEED"
for genome in $GENOMES; do
    if ! prepare_genome "$genome"; then
        CTX="$genome"; bad "could not prepare genome (skipped)"; continue
    fi
    san="$SAN"
    bp=$(grep -v '^>' "$san" | tr -d '\n' | wc -c)
    cbin="$(construct_bin "$genome")"; qbin="$(query_bin "$genome")"
    log "$genome: $bp ACGT bp after sanitization"
    for km in $KM; do
        k="${km%%,*}"; m="${km##*,}"
        run_one "$genome" "$k" "$m" "$san" "$cbin" "$qbin" "$bp"
    done
done

# ---- summary ----
echo >&2
log "================= PERFORMANCE / RESOURCES ================="
if ldd "$REL_BIN" 2>/dev/null | grep -qi tbb; then
    log "query parallelism: TBB linked (std::execution::par is parallel)"
else
    log "query parallelism: ${C_R}TBB NOT linked${C_0} -> query runs SERIAL (install libtbb-dev + rebuild for parallel query)"
fi
log "selfquery_s/query_Mkmer/s below are over the ${SELFQUERY_BP}bp self-query sample, not the whole genome"
{
    printf 'genome\tk\tm\tbp\tskmers\tconstruct_s\tconstruct_MB\tconstruct_Mkmer/s\tselfquery_s\tquery_Mkmer/s\n'
    cat "$PERF_TSV"
} | column -t -s$'\t' >&2

echo >&2
log "===================== SUMMARY ====================="
log "checks passed: $PASS    failed: $FAIL"
if [[ "$FAIL" -gt 0 ]]; then
    for f in "${FAILED[@]}"; do printf '  %sFAIL%s %s\n' "$C_R" "$C_0" "$f" >&2; done
    exit 1
fi
log "${C_G}ALL CHECKS PASSED${C_0}"
exit 0
