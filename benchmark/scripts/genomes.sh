#!/usr/bin/env bash
# Shared genome catalogue + fetch/sanitize logic for the sklib benchmark harness.
#
# Single source of truth: benchmark/data/genomes.tsv. This file is *sourced* by
# lib.sh and large_scale_e2e.sh (so they share one catalogue instead of each
# carrying a drifting copy) and run standalone by fetch_genomes.sh.
#
# Public surface:
#   genome_catalogue_names                 -> prints catalogued names, one per line
#   genome_row <name>                      -> sets GENOME_{NAME,KIND,SRC,SIZE,DESC}; rc!=0 if unknown
#   fetch_genome <name> [--force]          -> downloads+decompresses raw FASTA; sets RAW
#   prepare_genome <name> [--force]        -> fetch + sanitize; sets SAN (the cached path)
#
# Paths are derived from this file's location but honour an already-set
# BENCH_REPO_ROOT / BENCH_GEN_DIR / BENCH_HELPER (the experiment drivers export them).

_GENOMES_SH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
: "${BENCH_REPO_ROOT:=$(cd -- "$_GENOMES_SH_DIR/../.." && pwd)}"
GENOME_CATALOGUE="${GENOME_CATALOGUE:-$BENCH_REPO_ROOT/benchmark/data/genomes.tsv}"
GENOME_DIR="${BENCH_GEN_DIR:-$BENCH_REPO_ROOT/benchmark/data/genomes}"
GENOME_HELPER="${BENCH_HELPER:-$_GENOMES_SH_DIR/e2e_helpers.py}"
mkdir -p "$GENOME_DIR"

# Reuse the caller's logging if present (lib.sh / large_scale_e2e.sh have coloured
# log/warn/die); otherwise provide a plain fallback for standalone use.
declare -F log  >/dev/null || log()  { printf '[genomes] %s\n' "$*" >&2; }
declare -F warn >/dev/null || warn() { printf '[genomes] %s\n' "$*" >&2; }
declare -F die  >/dev/null || die()  { printf '[genomes FATAL] %s\n' "$*" >&2; exit 2; }

# Resolve <name> against the catalogue. Sets GENOME_NAME/KIND/SRC/SIZE/DESC globals.
# Returns 1 if the name is absent, 2 if the catalogue file is missing.
genome_row() {
    local name="$1" line
    [[ -f "$GENOME_CATALOGUE" ]] || { warn "catalogue not found: $GENOME_CATALOGUE"; return 2; }
    line="$(awk -v n="$name" '!/^[[:space:]]*#/ && NF && $1==n {print; exit}' "$GENOME_CATALOGUE")"
    [[ -n "$line" ]] || return 1
    # `read` puts the trailing remainder (the free-text description) into the last var.
    read -r GENOME_NAME GENOME_KIND GENOME_SRC GENOME_SIZE GENOME_DESC <<<"$line"
    return 0
}

genome_catalogue_names() {
    [[ -f "$GENOME_CATALOGUE" ]] || { warn "catalogue not found: $GENOME_CATALOGUE"; return 1; }
    awk '!/^[[:space:]]*#/ && NF {print $1}' "$GENOME_CATALOGUE"
}

# Download (and decompress) the raw FASTA for <name> into $GENOME_DIR/<name>.fa.
# Idempotent: a cached raw is reused unless --force. Sets RAW to the FASTA path.
fetch_genome() {
    local name="$1" force=0
    [[ "${2:-}" == "--force" ]] && force=1
    genome_row "$name" || { warn "$name: not in catalogue ($GENOME_CATALOGUE)"; return 1; }
    local raw="$GENOME_DIR/$name.fa"
    RAW="$raw"
    if (( ! force )) && [[ -s "$raw" ]]; then log "$name: raw FASTA cached"; return 0; fi
    command -v curl >/dev/null 2>&1 || { warn "curl not found"; return 1; }
    case "$GENOME_KIND" in
        url)
            command -v gzip >/dev/null 2>&1 || { warn "gzip not found"; return 1; }
            local gz="$GENOME_DIR/$name.fa.gz"
            log "$name: downloading $GENOME_SRC"
            curl -fSL --retry 3 "$GENOME_SRC" -o "$gz" || { warn "$name: download failed"; return 1; }
            gzip -dc "$gz" > "$raw.tmp" && mv "$raw.tmp" "$raw" \
                || { warn "$name: gunzip failed"; rm -f "$raw.tmp"; return 1; }
            rm -f "$gz"   # keep only the decompressed FASTA (chm13's .gz is ~900 MB)
            ;;
        ncbi)
            # NCBI E-utilities efetch: resolves RefSeq (NC_*) accessions directly to FASTA.
            local url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${GENOME_SRC}&rettype=fasta&retmode=text"
            log "$name: downloading NCBI $GENOME_SRC"
            curl -fSL --retry 3 "$url" -o "$raw.tmp" \
                || { warn "$name: NCBI download failed"; rm -f "$raw.tmp"; return 1; }
            [[ -s "$raw.tmp" ]] || { warn "$name: NCBI returned an empty FASTA for $GENOME_SRC"; rm -f "$raw.tmp"; return 1; }
            mv "$raw.tmp" "$raw"
            ;;
        *)
            warn "$name: unknown kind '$GENOME_KIND' (expected url|ncbi)"; return 1 ;;
    esac
    return 0
}

# Resolve <name> to a sanitized (uppercased, split at non-ACGT runs) FASTA, cached at
# $GENOME_DIR/<name>.sanitized.fa. Sets the global SAN. --force re-fetches + re-sanitizes.
prepare_genome() {
    local name="$1" force=0
    [[ "${2:-}" == "--force" ]] && force=1
    local san="$GENOME_DIR/$name.sanitized.fa"
    SAN="$san"
    if (( ! force )) && [[ -s "$san" ]]; then log "$name: using cached sanitized FASTA"; return 0; fi
    local fargs=(); (( force )) && fargs=(--force)
    fetch_genome "$name" "${fargs[@]}" || return 1
    log "$name: sanitizing (uppercase + split at non-ACGT runs)"
    python3 "$GENOME_HELPER" sanitize "$RAW" > "$san.tmp" && mv "$san.tmp" "$san" \
        || { warn "$name: sanitize failed"; rm -f "$san.tmp"; return 1; }
    return 0
}
