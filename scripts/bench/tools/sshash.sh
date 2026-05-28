#!/usr/bin/env bash
# SSHash adapter — exact, associative k-mer dictionary (github.com/jermp/sshash).
# Indexes the distinct k-mers of the de Bruijn graph; input = BCALM2 unitigs so the
# k-mer set matches sklib's exactly (KMC validates set equality elsewhere).
#
# VERIFY against `$SSHASH_BIN --help` / `build --help` / `query --help` after building;
# flag names have changed across SSHash releases. Current best-effort:
#   build : sshash build -i <unitigs.fa> -k <k> -m <minimizer> -o <index>
#   query : sshash query -i <index> -q <queries.fa>
# Construction time is bcalm (sequence->unitigs) + sshash build (unitigs->index), since
# that is the honest FASTA->index cost; peak RSS is the max of the two stages.
# Competitor => uses_m_sshash stays 1 (it ignores sklib's m; it has its own minimizer).

available_sshash() { [[ -n "${SSHASH_BIN:-}" && -x "${SSHASH_BIN:-}" ]] && bcalm_available; }

version_sshash() {
    local v; v="$("$SSHASH_BIN" --version 2>/dev/null | head -1)"
    echo "sshash-${v:-?}"
}

# SSHash's own minimizer length (independent of sklib's m). Typical: ~k-6, capped [1,20].
_sshash_m() { local mm=$(( $1 - 6 )); ((mm < 1)) && mm=1; ((mm > 20)) && mm=20; echo "$mm"; }

construct_sshash() {
    local san="$1" k="$2" m="$3" wd="$4" uni msh bsec brss
    uni="$(bcalm_unitigs "$san" "$k" "$wd/bcalm")" || { warn "sshash: bcalm unitigs failed"; return 1; }
    bsec="$BCALM_SEC"; brss="$BCALM_RSS"
    IDX_PATH="$wd/index.sshash"; msh="$(_sshash_m "$k")"
    # --canonical: match sklib/KMC canonical-k-mer membership (a k-mer == its rev-comp).
    run_timed "$SSHASH_BIN" build -i "$uni" -k "$k" -m "$msh" --canonical -o "$IDX_PATH" || { warn "sshash: build failed"; return 1; }
    [[ -s "$IDX_PATH" ]] || return 1
    RUN_SEC=$(awk -v a="$bsec" -v b="$RUN_SEC" 'BEGIN{printf "%.3f", a+b}')   # bcalm + sshash build
    (( brss > RUN_RSS_KB )) && RUN_RSS_KB="$brss"
    IDX_FILE_BYTES=$(stat -c%s "$IDX_PATH"); IDX_PAYLOAD_BYTES="$IDX_FILE_BYTES"; N_SKMERS=NA
    return 0
}

query_sshash() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    # --multiline: streaming/reads workloads are wrapped FASTA (>1 line/record), else
    # SSHash would only read the first line and undercount k-mers vs our n_queries.
    run_timed_median taskset -c "$cpus" "$SSHASH_BIN" query -i "$idx" -q "$qfa" --multiline
}
