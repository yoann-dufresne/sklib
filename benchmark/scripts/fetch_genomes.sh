#!/usr/bin/env bash
# Download (and sanitize) the benchmark genomes catalogued in benchmark/data/genomes.tsv.
#
# Modular by design: fetch everything, a named subset, or force-refresh specific genomes
# (a corrupt download or an updated reference) without touching the rest. The catalogue
# itself lives in benchmark/data/genomes.tsv; the per-genome fetch logic in genomes.sh.
#
#   fetch_genomes.sh                  # all catalogued genomes (skips ones already cached)
#   fetch_genomes.sh ecoli chr21      # only these
#   fetch_genomes.sh --force chr21    # re-download chr21 (corruption / updated reference)
#   fetch_genomes.sh --raw-only ...   # download + decompress but skip the sanitize step
#   fetch_genomes.sh --list           # print the catalogue and exit
#
# Each genome is independent: one failure is reported and the rest still run.
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark/scripts/genomes.sh
source "$SCRIPT_DIR/genomes.sh"

usage() { sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

force=0 raw_only=0
names=()
while (( $# )); do
    case "$1" in
        --force)    force=1 ;;
        --raw-only) raw_only=1 ;;
        --list)
            printf '%-10s %-5s %-8s %s\n' NAME KIND SIZE SOURCE
            while read -r n; do
                genome_row "$n" || continue
                printf '%-10s %-5s %-8s %s\n' "$n" "$GENOME_KIND" "$GENOME_SIZE" "$GENOME_SRC"
            done < <(genome_catalogue_names)
            exit 0 ;;
        -h|--help)  usage; exit 0 ;;
        --)         shift; while (( $# )); do names+=("$1"); shift; done; break ;;
        -*)         die "unknown option: $1 (try --help)" ;;
        *)          names+=("$1") ;;
    esac
    shift
done

# No names given -> the whole catalogue.
if (( ${#names[@]} == 0 )); then
    mapfile -t names < <(genome_catalogue_names)
fi
(( ${#names[@]} > 0 )) || die "no genomes to fetch (empty/missing catalogue: $GENOME_CATALOGUE)"

fargs=(); (( force )) && fargs=(--force)
ok=0 fail=0; failed=()
for n in "${names[@]}"; do
    if (( raw_only )); then fetch_genome  "$n" "${fargs[@]}"; rc=$?
    else                    prepare_genome "$n" "${fargs[@]}"; rc=$?; fi
    if (( rc == 0 )); then ok=$((ok+1)); else fail=$((fail+1)); failed+=("$n"); fi
done

(( fail == 0 )) || warn "failed: ${failed[*]}"
log "done: $ok ok, $fail failed (cache: $GENOME_DIR)"
(( fail == 0 ))
