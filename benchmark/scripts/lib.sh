#!/usr/bin/env bash
# Shared library for the sklib benchmark harness (construct.sh / query_*.sh / setop.sh).
# Holds the genome catalogue, timing/measurement helpers, the KMC oracle, and
# CSV utilities. Sourced -- not executed directly.
#
# It deliberately reuses the same genome cache directory as the correctness
# harness (benchmark/scripts/large_scale_e2e.sh) so multi-GB downloads are shared.

# ---- paths (BENCH_REPO_ROOT must be set by the caller) ----
: "${BENCH_REPO_ROOT:?source lib.sh with BENCH_REPO_ROOT set}"
BENCH_HELPER="$BENCH_REPO_ROOT/benchmark/scripts/e2e_helpers.py"
BENCH_GEN_DIR="$BENCH_REPO_ROOT/benchmark/data/genomes"   # shared with large_scale_e2e.sh
BENCH_KMC_TMP="$BENCH_REPO_ROOT/benchmark/results/latest/kmc_tmp"
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

# =====================================================================
#  Experiment harness shared layer
#  Used by construct.sh / query_single.sh / query_stream.sh / setop.sh.
#  Provides: the shared grids, the index/query/mutant caches, and the
#  RESUMABILITY contract -- an interrupted run restarts without redoing
#  work, and (same tool version + same host) rows are trusted as final.
# =====================================================================
BENCH_SCRIPTS="$(dirname "$BENCH_HELPER")"      # benchmark/scripts (mutate.py lives here)

# ---- run identity: a measurement is keyed by (tool, tool_version, host, ...) ----
HOST="${HOST:-$(hostname)}"
CPU="$(csv_escape "$(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo)")"
NPROC="$(nproc)"

# ---- shared grids (every one env-overridable to scope a run) ----
DATASETS="${DATASETS:-sarscov2 ecoli yeast celegans chr21 chr1}"
TOOLS="${TOOLS:-sklib sshash sbwt cbl bqf fmsi kmc}"    # gated per experiment by can_{construct,query,setop}_*
KM="${KM:-15,7 21,11 31,15 41,19 51,25 63,31}"          # whole-range (k,m) ladder, m~k/2
THREADS="${THREADS:-1 2 4 8 16}"                         # uniform thread sweep (1 = mono)
PRESENCE="${PRESENCE:-0 25 50 75 100}"                   # query present-fraction sweep (%)
JACCARD="${JACCARD:-0 0.1 0.3 0.5 0.7 0.9 1.0}"          # set-op overlap sweep (target J)
REPS="${REPS:-3}"; SEED="${SEED:-1234}"
N_QUERY="${N_QUERY:-200000}"          # individual k-mers per (presence) point
STREAM_RECS="${STREAM_RECS:-2000}"    # stream records per (presence) point
STREAM_LEN="${STREAM_LEN:-300}"       # bp per stream record

# ---- result + cache locations (under the git-ignored fresh-output dir) ----
RESULTS="${RESULTS:-$BENCH_REPO_ROOT/benchmark/results/latest}"
IDX_CACHE="${IDX_CACHE:-$RESULTS/indexes}"      # built indexes, keyed by tool+version+params
QCACHE="${QCACHE:-$RESULTS/queries}"            # generated query FASTAs
GMUT="${GMUT:-$BENCH_GEN_DIR/mutants}"          # mutated genomes for the Jaccard sweep
mkdir -p "$RESULTS" "$IDX_CACHE" "$QCACHE" "$GMUT"

# ---- resumability -----------------------------------------------------------
# CSVs are append-only. On (re)start, load_done() reads the rows already present and
# keys them by the identity columns; is_done() then skips any measurement whose
# (tool, tool_version, host, ...) key is already there. Because the key carries the
# tool version and the host, a different build or a different machine yields a
# different key and is re-measured -- exactly "same tool + same machine => keep".
declare -gA _DONE
mk_key() { local IFS='|'; printf '%s' "$*"; }          # stable join of identity fields
load_done() {                                          # load_done <csv> <key-colname...>
    local csv="$1"; shift; _DONE=()
    [[ -f "$csv" ]] || return 0
    local key
    while IFS= read -r key; do [[ -n "$key" ]] && _DONE["$key"]=1; done < <(
        awk -F, -v want="$*" '
            NR==1 { n=split(want,wc," "); for(i=1;i<=NF;i++) H[$i]=i; next }
            { k=""; for(j=1;j<=n;j++) k=k (j>1?"|":"") $(H[wc[j]]); print k }' "$csv")
    log "resume: ${#_DONE[@]} row(s) already in $(basename "$csv"); matching (tool,version,host,…) skipped"
}
is_done()   { [[ -n "${_DONE[$1]:-}" ]]; }
mark_done() { _DONE["$1"]=1; }
csv_init()  { [[ -f "$1" ]] || { mkdir -p "$(dirname "$1")"; printf '%s\n' "$2" > "$1"; }; }  # header once

# ---- index cache (resumable build) -----------------------------------------
# ensure_index_fa builds tool's index for a sanitized FASTA into a cache keyed by
# (tool, version, tag, k, m) and remembers IDX_PATH/sizes in a .idxmeta sidecar, so a
# restart reuses it instead of rebuilding. construct_<tool> (tools.sh) does the build.
ensure_index_fa() {                                    # <tool> <san_fa> <tag> <k> <m>
    local tool="$1" san="$2" tag="$3" k="$4" m="$5" ver cd meta
    ver="$(version_"$tool" 2>/dev/null | tr -c 'A-Za-z0-9._-' '_')"
    cd="$IDX_CACHE/$tool/${ver:-x}/$tag.k$k.m$m"; meta="$cd/.idxmeta"
    if [[ -f "$meta" ]]; then source "$meta"; [[ -s "$IDX_PATH" ]] && return 0; fi
    mkdir -p "$cd"
    construct_"$tool" "$san" "$k" "$m" "$cd" || return 1
    { printf "IDX_PATH=%q\n" "$IDX_PATH"; printf "IDX_FILE_BYTES=%q\n" "${IDX_FILE_BYTES:-0}"
      printf "IDX_PAYLOAD_BYTES=%q\n" "${IDX_PAYLOAD_BYTES:-0}"; printf "N_SKMERS=%q\n" "${N_SKMERS:-NA}"; } > "$meta"
    return 0
}
ensure_index() { prepare_genome "$2" || return 1; ensure_index_fa "$1" "$SAN" "$2" "$3" "$4"; }  # <tool> <dataset> <k> <m>

# ---- query-set generators (cached; regenerated only if missing) ------------
# qset_single: n one-k-mer-per-record queries at target present-fraction p% -- present
# k-mers sampled from the genome, absent ones random, shuffled together. Sets QSET_FA,
# QSET_PRES, QSET_ABS (exact counts, the throughput denominators).
qset_single() {                                        # <san> <k> <p> <n> <seed>
    local san="$1" k="$2" p="$3" n="$4" seed="$5" fa
    fa="$QCACHE/$(basename "${san%.fa}").k$k.single.p$p.n$n.s$seed.fa"
    QSET_FA="$fa"; QSET_PRES=$(( n * p / 100 )); QSET_ABS=$(( n - QSET_PRES ))
    [[ -s "$fa" ]] && return 0
    { (( QSET_PRES > 0 )) && python3 "$BENCH_HELPER" sample_positive "$san" "$k" "$QSET_PRES" "$seed"
      (( QSET_ABS  > 0 )) && python3 "$BENCH_HELPER" randkmers "$k" "$QSET_ABS" "$seed"; } \
      | python3 "$BENCH_HELPER" shuffle /dev/stdin "$seed" | seq2fa > "$fa.tmp" && mv "$fa.tmp" "$fa"
}
# qset_stream: nrec sequence records of len bp at target present-fraction p% -- present
# records are clean genome substrings (simreads err=0), absent ones fully substituted
# (err=1). Each record is homogeneous, so the k-mer present-fraction is exactly p.
qset_stream() {                                        # <san> <k> <p> <nrec> <len> <seed>
    local san="$1" k="$2" p="$3" nrec="$4" len="$5" seed="$6" fa npres nabs kpr
    fa="$QCACHE/$(basename "${san%.fa}").k$k.stream.p$p.r$nrec.l$len.s$seed.fa"
    QSET_FA="$fa"; npres=$(( nrec * p / 100 )); nabs=$(( nrec - npres )); kpr=$(( len - k + 1 ))
    QSET_PRES=$(( npres * kpr )); QSET_ABS=$(( nabs * kpr ))
    [[ -s "$fa" ]] && return 0
    { (( npres > 0 )) && python3 "$BENCH_HELPER" simreads "$san" "$len" "$npres" 0.0 "$seed"
      (( nabs  > 0 )) && python3 "$BENCH_HELPER" simreads "$san" "$len" "$nabs" 1.0 "$((seed+1))"; } \
      > "$fa.tmp" && mv "$fa.tmp" "$fa"
}

# ---- Jaccard sweep: mutated B genomes (cached) -----------------------------
# rate(J,k): a per-base substitution rate giving ~target Jaccard J between A and its
# mutant B. A k-mer survives iff none of its k bases mutate => shared fraction
# s=(1-r)^k, and J=s/(2-s) for equal-size sets, so r = 1 - (2J/(1+J))^(1/k).
jaccard_to_rate() { awk -v J="$1" -v k="$2" 'BEGIN{
    if(J>=1){print "0.000000";exit} s=2*J/(1+J); if(s<=0){print "1.000000";exit}
    printf "%.6f", 1-exp(log(s)/k) }'; }
ensure_mutant() {                                      # <san> <rate> <seed> -> MUT_FA
    local san="$1" rate="$2" seed="$3" fa
    fa="$GMUT/$(basename "${san%.fa}").mut$rate.s$seed.fa"; MUT_FA="$fa"
    [[ -s "$fa" ]] && return 0
    python3 "$BENCH_SCRIPTS/mutate.py" "$san" "$rate" "$seed" > "$fa.tmp" && mv "$fa.tmp" "$fa"
}
