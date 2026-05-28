#!/usr/bin/env bash
# SBWT adapter — exact, very compact k-mer membership via the Spectral Burrows-Wheeler
# Transform (github.com/algbio/SBWT). Builds directly from the sanitized FASTA (SBWT
# dedups k-mers); --add-reverse-complements gives canonical membership (a k-mer or its
# reverse complement counts as present), matching sklib/KMC.
#
# CLI (verified): build -i <fa> -o <idx> -k <k> [--variant V] [--add-reverse-complements]
#                 search -i <idx> -q <q.fa> -o <out>
# Uses the default plain-matrix variant (fast, robust). NOTE: the binary was built
# without BMI2, so the compact Elias-Fano variants (rrr-/mef-*) would be slow; for
# SBWT's best *space* point, rebuild with -march=native and set SBWT_VARIANT=mef-matrix.
# k must be <= 32 (MAX_KMER_LENGTH). Competitor => uses_m_sbwt stays 1.

SBWT_VARIANT="${SBWT_VARIANT:-plain-matrix}"

available_sbwt() { [[ -n "${SBWT_BIN:-}" && -x "${SBWT_BIN:-}" ]]; }
version_sbwt() { echo "sbwt-$SBWT_VARIANT"; }

construct_sbwt() {
    local san="$1" k="$2" m="$3" wd="$4"
    if (( k > 32 )); then warn "sbwt: k>32 unsupported (MAX_KMER_LENGTH=32)"; return 1; fi
    IDX_PATH="$wd/index.sbwt"
    run_timed "$SBWT_BIN" build -i "$san" -o "$IDX_PATH" -k "$k" \
        --variant "$SBWT_VARIANT" --add-reverse-complements --temp-dir "$wd" -t 1 \
        || { warn "sbwt: build failed"; return 1; }
    [[ -s "$IDX_PATH" ]] || return 1
    IDX_FILE_BYTES=$(stat -c%s "$IDX_PATH"); IDX_PAYLOAD_BYTES="$IDX_FILE_BYTES"; N_SKMERS=NA
    return 0
}

query_sbwt() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    run_timed_median taskset -c "$cpus" "$SBWT_BIN" search -i "$idx" -q "$qfa" -o /dev/null
}
