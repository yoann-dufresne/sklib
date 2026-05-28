#!/usr/bin/env bash
# Shared library for the sklib benchmark harness (scripts/bench/bench.sh).
# Holds the genome catalogue, timing/measurement helpers, the KMC oracle, and
# CSV utilities. Sourced -- not executed directly.
#
# It deliberately reuses the same genome cache directory as the correctness
# harness (scripts/large_scale_e2e.sh) so multi-GB downloads are shared.

# ---- paths (BENCH_REPO_ROOT must be set by the caller) ----
: "${BENCH_REPO_ROOT:?source lib.sh with BENCH_REPO_ROOT set}"
BENCH_HELPER="$BENCH_REPO_ROOT/scripts/e2e_helpers.py"
BENCH_GEN_DIR="$BENCH_REPO_ROOT/scripts/out/e2e/genomes"   # shared with large_scale_e2e.sh
BENCH_KMC_TMP="$BENCH_REPO_ROOT/scripts/out/bench/kmc_tmp"
mkdir -p "$BENCH_GEN_DIR" "$BENCH_KMC_TMP"

# ---- pretty logging ----
if [[ -t 2 ]]; then C_G=$'\033[32m'; C_R=$'\033[31m'; C_Y=$'\033[33m'; C_B=$'\033[1m'; C_0=$'\033[0m'
else C_G=; C_R=; C_Y=; C_B=; C_0=; fi
log()  { printf '%s[bench]%s %s\n' "$C_B" "$C_0" "$*" >&2; }
warn() { printf '%s[bench]%s %s\n' "$C_Y" "$C_0" "$*" >&2; }
die()  { printf '%s[bench FATAL]%s %s\n' "$C_R" "$C_0" "$*" >&2; exit 2; }

need_tools() {
    local t
    for t in "$@"; do command -v "$t" >/dev/null 2>&1 || die "required tool not found: $t"; done
}

# tests/sequences_2_fa.sh is a slow bash loop; wrap raw sequence lines into FASTA with awk.
seq2fa() { awk '{print ">s" NR; print $0}'; }

# ---- genome catalogue -----------------------------------------------------
# Remote genomes resolvable by a single URL.
genome_url() {
    case "$1" in
        yeast)    echo "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz" ;;
        celegans) echo "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz" ;;
        chr21)    echo "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz" ;;
        chr20)    echo "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz" ;;
        chr1)     echo "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz" ;;
        chm13)    echo "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz" ;;
        *)        echo "" ;;
    esac
}
# Local genome fixtures (no download).
genome_local() {
    case "$1" in
        ecoli)    echo "$BENCH_REPO_ROOT/data/ecoli.fa" ;;
        sarscov2) echo "$BENCH_REPO_ROOT/data/sarscov2.fa" ;;
        *)        echo "" ;;
    esac
}

# Resolve a genome name to a sanitized (uppercased, split at non-ACGT) FASTA path.
# Sets the global SAN. Returns non-zero on failure.
prepare_genome() {
    local name="$1" raw san url local_path
    san="$BENCH_GEN_DIR/$name.sanitized.fa"
    SAN="$san"
    if [[ -f "$san" ]]; then log "$name: using cached sanitized FASTA"; return 0; fi

    local_path="$(genome_local "$name")"
    if [[ -n "$local_path" ]]; then
        raw="$local_path"
        [[ -f "$raw" ]] || { warn "$name: local fixture missing: $raw"; return 1; }
    else
        url="$(genome_url "$name")"
        [[ -n "$url" ]] || { warn "$name: unknown genome (no URL/fixture)"; return 1; }
        raw="$BENCH_GEN_DIR/$name.fa"
        if [[ ! -f "$raw" ]]; then
            log "$name: downloading $url"
            curl -fSL --retry 3 "$url" -o "$BENCH_GEN_DIR/$name.fa.gz" || { warn "$name: download failed"; return 1; }
            gzip -dc "$BENCH_GEN_DIR/$name.fa.gz" > "$raw" || { warn "$name: gunzip failed"; return 1; }
        fi
    fi
    log "$name: sanitizing (uppercase + split at non-ACGT runs)"
    python3 "$BENCH_HELPER" sanitize "$raw" > "$san.tmp" && mv "$san.tmp" "$san" || { warn "$name: sanitize failed"; return 1; }
    return 0
}

# ---- measurement ----------------------------------------------------------
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
    RUN_SEC=$(awk "BEGIN{printf \"%.3f\", $t1-$t0}")
    RUN_RSS_KB=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog")
    [[ -n "$RUN_RSS_KB" ]] || RUN_RSS_KB=0
    rm -f "$tlog"
    return "$status"
}

# Median wall time over N reps of a command (pins via the caller's wrapper).
# Sets RUN_SEC to the median and RUN_RSS_KB to the max RSS seen. Returns last status.
RUN_REPS=1
run_timed_median() {
    local reps="$RUN_REPS" i status=0 maxrss=0
    local times=()
    for ((i = 0; i < reps; i++)); do
        run_timed "$@"; status=$?
        times+=("$RUN_SEC")
        (( RUN_RSS_KB > maxrss )) && maxrss="$RUN_RSS_KB"
    done
    RUN_SEC=$(printf '%s\n' "${times[@]}" | sort -g | awk '{a[NR]=$1} END{print (NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}')
    RUN_RSS_KB="$maxrss"
    return "$status"
}

# ---- KMC oracle -----------------------------------------------------------
# Count distinct canonical k-mers in a FASTA (the bits/kmer denominator).
# Echoes the count.  $1=fasta $2=k
kmc_distinct_count() {
    local in="$1" kk="$2" db="$BENCH_KMC_TMP/cnt_$$_$RANDOM"
    if ! kmc -k"$kk" -ci1 -fm "$in" "$db" "$BENCH_KMC_TMP" >/dev/null 2>&1; then
        rm -f "$db".kmc_pre "$db".kmc_suf; echo "NA"; return 1
    fi
    # #unique k-mers counted is reported in the KMC database stats.
    local n
    n=$(kmc_tools transform "$db" dump /dev/stdout 2>/dev/null | wc -l)
    rm -f "$db".kmc_pre "$db".kmc_suf
    echo "$n"
}

# Dump the sorted-unique canonical k-mer set of a FASTA. $1=fasta $2=fmt(-fa/-fm) $3=k $4=out
kmc_canon_set() {
    local in="$1" fmt="$2" kk="$3" out="$4" db="$BENCH_KMC_TMP/set_$$_$RANDOM"
    kmc -k"$kk" -ci1 "$fmt" "$in" "$db" "$BENCH_KMC_TMP" >/dev/null 2>&1 || { rm -f "$db".kmc_*; return 1; }
    kmc_tools transform "$db" dump /dev/stdout 2>/dev/null | cut -f1 | LC_ALL=C sort -u > "$out"
    local rc=${PIPESTATUS[0]}
    rm -f "$db".kmc_pre "$db".kmc_suf
    return "$rc"
}

# ---- small numeric helpers ----
# bits per k-mer = 8 * payload_bytes / distinct_kmers
bits_per_kmer() { awk -v b="$1" -v n="$2" 'BEGIN{ if(n+0>0) printf "%.4f", 8*b/n; else printf "NA" }'; }
div()  { awk -v a="$1" -v b="$2" 'BEGIN{ if(b+0>0) printf "%.4f", a/b; else printf "NA" }'; }
# throughput in Mkmer/s = n / sec / 1e6
mrate() { awk -v n="$1" -v s="$2" 'BEGIN{ if(s+0>0) printf "%.3f", n/s/1e6; else printf "NA" }'; }
human_mb() { awk "BEGIN{printf \"%.0f\", $1/1024}"; }

# Replace CSV-hostile characters in a free-text field.
csv_escape() { printf '%s' "$1" | tr ',\n\t' ';  '; }
