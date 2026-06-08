#!/usr/bin/env bash
# BCALM2 adapter — unitig construction. Not a benchmarked membership tool itself;
# it produces the stitched unitigs that are the *common input* for SSHash (and any
# other de-Bruijn-graph tool), guaranteeing all tools index the same k-mer set.
#
# VERIFY against `$BCALM_BIN --help` after building:
#   expected: bcalm -in <in.fa> -kmer-size <k> -abundance-min 1 -out <prefix>
#   output  : <prefix>.unitigs.fa   (one record per unitig; every distinct k-mer once)
# Sourced by tools.sh; BCALM_BIN comes from tools_src/tools.env.

bcalm_available() { [[ -n "${BCALM_BIN:-}" && -x "${BCALM_BIN:-}" ]]; }

# bcalm_unitigs <in.fa> <k> <out_prefix>
# Echoes the unitigs FASTA path; sets BCALM_SEC / BCALM_RSS (build cost, for fair
# construction timing of tools that consume unitigs). Cached if already built.
BCALM_SEC=0; BCALM_RSS=0
bcalm_unitigs() {
    local in="$1" k="$2" prefix="$3" uni="$3.unitigs.fa" prev="$PWD" st
    if [[ ! -s "$uni" ]]; then
        # bcalm scatters temp files in CWD -> run from the prefix dir (run_timed sets
        # globals in *this* shell, so cd here rather than in a subshell). Inputs are
        # absolute paths, so the relative -out basename is the only path that moves.
        cd "$(dirname "$prefix")" || return 1
        # Default 1 core: construction time stays comparable to sklib's serial build.
        # Set BCALM_CORES to parallelize the unitig step (e.g. for very large genomes).
        run_timed "$BCALM_BIN" -in "$in" -kmer-size "$k" \
            -abundance-min 1 -out "$(basename "$prefix")" -nb-cores "${BCALM_CORES:-1}"
        st=$?
        cd "$prev" || return 1
        BCALM_SEC="$RUN_SEC"; BCALM_RSS="$RUN_RSS_KB"
        [[ "$st" -eq 0 ]] || return 1
    fi
    [[ -s "$uni" ]] && { echo "$uni"; return 0; }
    return 1
}
