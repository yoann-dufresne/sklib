#!/usr/bin/env bash
# Download (and sanitize) the real read sets catalogued in benchmark/data/reads.tsv.
#
# Each set pools one or more ENA runs; every FASTQ of every run is fetched from ENA, checked
# against the ENA md5, then all of them are converted into ONE sanitized FASTA
# (benchmark/data/genomes/<set>.sanitized.fa) that the harness reuses verbatim, exactly like a
# catalogued genome. "Sanitized" = uppercase + split at non-ACGT runs, one record per fragment,
# so sklib and KMC see the same k-mers. R1/R2 pairing is dropped: a k-mer set ignores it.
#
#   fetch_reads.sh                     # every catalogued set (skips what is already done)
#   fetch_reads.sh hg002_A hg002_B     # only these
#   fetch_reads.sh --download-only …   # fetch + md5 only, no FASTA
#   fetch_reads.sh --list              # print the catalogue and exit
#
# ENA throttles each connection (~1.3 MB/s measured 2026-10-05) but not the aggregate, so
# every file is fetched as CHUNKS byte ranges and PAR files run at once (PAR*CHUNKS
# connections, 16 by default). Interrupted downloads resume chunk by chunk.
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
CATALOGUE="${READS_TSV:-$REPO/benchmark/data/reads.tsv}"
GEN="${GENOMES_DIR:-$REPO/benchmark/data/genomes}"
RAW="${READS_RAW_DIR:-$GEN/reads_raw}"
PAR="${PAR:-4}"
CHUNKS="${CHUNKS:-4}"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

log() { printf '[fetch_reads %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
usage() { sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

catalogue_names() { awk '!/^#/ && NF {print $1}' "$CATALOGUE"; }
catalogue_runs()  { awk -v n="$1" '!/^#/ && $1==n {print $2}' "$CATALOGUE" | tr ',' ' '; }

# One file: CHUNKS ranged curls into part files (each resumable), concatenated, md5-checked.
# A finished file carries a .ok marker, so a re-run skips it.
fetch_file() {   # url md5 bytes dest
    local url="$1" md5="$2" bytes="$3" dest="$4" i lo hi pids=() step parts=()
    [[ -f "$dest.ok" ]] && return 0
    step=$(( (bytes + CHUNKS - 1) / CHUNKS ))
    for (( i = 0; i < CHUNKS; i++ )); do
        lo=$(( i * step )); hi=$(( lo + step - 1 )); (( hi >= bytes )) && hi=$(( bytes - 1 ))
        parts+=( "$dest.part$i" )
        # A part already at its full length is done; otherwise resume it from where it stopped.
        if [[ -f "$dest.part$i" ]] && (( $(stat -c%s "$dest.part$i") == hi - lo + 1 )); then continue; fi
        (
            # No curl --retry: a retry would re-append the already-written bytes to the part.
            # Each attempt instead resumes from the part's current length.
            for attempt in 1 2 3 4 5 6 7 8 9 10; do
                have=0; [[ -f "$dest.part$i" ]] && have=$(stat -c%s "$dest.part$i")
                (( have == hi - lo + 1 )) && exit 0
                curl -sS --fail --connect-timeout 30 --speed-limit 10000 --speed-time 120 \
                     -r "$(( lo + have ))-$hi" "$url" >> "$dest.part$i" && continue
                sleep $(( attempt * 10 ))
            done
            have=$(stat -c%s "$dest.part$i" 2>/dev/null || echo 0)
            (( have == hi - lo + 1 ))
        ) &
        pids+=( $! )
    done
    local ok=1 p; for p in "${pids[@]}"; do wait "$p" || ok=0; done
    (( ok )) || { log "FAILED (chunks): $url"; return 1; }
    cat "${parts[@]}" > "$dest.tmp" && rm -f "${parts[@]}"
    local got; got=$(md5sum "$dest.tmp" | cut -d' ' -f1)
    if [[ "$got" != "$md5" ]]; then
        log "md5 MISMATCH $dest ($got != $md5): removed, re-run to retry"; rm -f "$dest.tmp"; return 1
    fi
    mv "$dest.tmp" "$dest" && touch "$dest.ok"
    log "ok $(basename "$dest") ($(( bytes / 1000000 )) MB)"
}
export -f fetch_file log
export CHUNKS

# FASTQ(.gz)… -> one sanitized FASTA. Pure stream: zcat | mawk (≈ the e2e_helpers.py sanitize
# rule, without its 80-column wrapping, which no reader needs). Records are named ">r".
sanitize_set() {   # name files…
    local name="$1"; shift
    local out="$GEN/$name.sanitized.fa"
    [[ -s "$out" ]] && { log "$name: $out exists, skip sanitize"; return 0; }
    log "$name: sanitizing $# files -> $out"
    zcat -- "$@" | mawk 'NR % 4 == 2 {
            n = split(toupper($0), f, /[^ACGT]+/)
            for (i = 1; i <= n; i++) if (f[i] != "") printf ">r\n%s\n", f[i]
        }' > "$out.tmp" || { rm -f "$out.tmp"; log "$name: sanitize FAILED"; return 1; }
    mv "$out.tmp" "$out"
    log "$name: done, $(( $(stat -c%s "$out") / 1000000 )) MB"
}

download_only=0; names=()
while (( $# )); do
    case "$1" in
        --download-only) download_only=1 ;;
        --list) awk '!/^#/ && NF {printf "%-10s %6s Gbp  %s\n", $1, $3, $2}' "$CATALOGUE"; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1 (try --help)" ;;
        *) names+=( "$1" ) ;;
    esac
    shift
done
(( ${#names[@]} )) || mapfile -t names < <(catalogue_names)
mkdir -p "$RAW" "$GEN"

# 1. resolve every file of every requested set through the ENA API
jobs="$RAW/.jobs.tsv"; : > "$jobs"
declare -A SET_FILES=()
for name in "${names[@]}"; do
    runs="$(catalogue_runs "$name")"; [[ -n "$runs" ]] || die "$name: not in $CATALOGUE"
    for run in $runs; do
        row=$(curl -sS --retry 5 "$ENA_API?accession=$run&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes&format=tsv" | awk 'NR==2')
        [[ -n "$row" ]] || die "$run: ENA filereport returned nothing"
        IFS=$'\t' read -r _ ftp md5s sizes <<<"$row"
        IFS=';' read -r -a F <<<"$ftp"; IFS=';' read -r -a M <<<"$md5s"; IFS=';' read -r -a S <<<"$sizes"
        for i in "${!F[@]}"; do
            dest="$RAW/$(basename "${F[$i]}")"
            printf '%s\t%s\t%s\t%s\n' "https://${F[$i]}" "${M[$i]}" "${S[$i]}" "$dest" >> "$jobs"
            SET_FILES[$name]+="$dest "
        done
    done
done
total=$(awk -F'\t' '{s+=$3} END{printf "%.1f", s/1e9}' "$jobs")
log "${#names[@]} sets, $(wc -l < "$jobs") files, $total GB to fetch (minus what is already there); PAR=$PAR CHUNKS=$CHUNKS"

# 2. download (PAR files at once, CHUNKS connections each)
awk -F'\t' '{print $1" "$2" "$3" "$4}' "$jobs" \
    | xargs -P "$PAR" -L 1 bash -c 'fetch_file "$0" "$1" "$2" "$3"'
fail=0
while IFS=$'\t' read -r _ _ _ dest; do [[ -f "$dest.ok" ]] || { log "missing: $dest"; fail=1; }; done < "$jobs"
(( fail )) && die "some downloads failed; re-run the same command to resume"
(( download_only )) && { log "download-only: done"; exit 0; }

# 3. one sanitized FASTA per set
for name in "${names[@]}"; do
    # shellcheck disable=SC2086
    sanitize_set "$name" ${SET_FILES[$name]} || fail=1
done
exit "$fail"
