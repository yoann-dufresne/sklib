#!/usr/bin/env bash
# FMSI adapter — "FroM Superstring to Indexing" (github.com/OndrejSladky/fmsi): exact,
# highly space-efficient k-mer membership via masked superstrings + the Masked Burrows-
# Wheeler Transform (MBWT). Two-stage build like sshash: KmerCamel (github.com/OndrejSladky/
# kmercamel) computes a near-shortest masked superstring, then `fmsi index` builds the MBWT
# over it. Construction cost = kmercamel + fmsi index (both serial); peak RSS = max of the
# two stages -- the honest FASTA->index cost, as for sshash (bcalm+build) / bqf (kmc+build).
#
# Canonical membership: KmerCamel is bidirectional BY DEFAULT (a k-mer and its reverse
# complement are the same) => matches KMC/sklib/cbl/sbwt; do NOT pass -u. `-O` at query is
# required because the mask was built with -M (maximum number of ones).
#
# CLI (verified against the built binary's -h and KmerCamel's source; the ms path is the
# positional index prefix -- the `-p` form in KmerCamel's README is stale):
#   kmercamel compute -k <k> -o /dev/null -M ms.fa <in.fa>   # max-ones masked superstring -> ms.fa
#   fmsi index -x -k <k> ms.fa                                # writes ms.fa.fmsi.* index files
#   fmsi query -O -q <q.fa> -k <k> ms.fa                      # one bitstring per query line (stdout)
# Build dep: KmerCamel needs GLPK (apt install libglpk-dev) for the mask-optimization ILP.
# k up to 127, no odd/even constraint. Competitor => uses_m_fmsi stays 1.

available_fmsi() {
    [[ -n "${FMSI_BIN:-}"      && -x "${FMSI_BIN:-}" \
    && -n "${KMERCAMEL_BIN:-}" && -x "${KMERCAMEL_BIN:-}" ]]
}
version_fmsi() {   # FMSI prints "Version: <hash>" in its no-arg banner (to stderr)
    local v; v="$("$FMSI_BIN" 2>&1 | awk '/^Version:/{print $2; exit}')"
    echo "fmsi-${v:-?}"
}

construct_fmsi() {
    local san="$1" k="$2" m="$3" wd="$4"
    local ms="$wd/ms.fa" csec crss f bytes=0   # separate line: a single `local` can't reference wd it just set (set -u)
    # 1) masked superstring (canonical default; -M = max-ones mask, -o /dev/null drops the default-mask output)
    run_timed "$KMERCAMEL_BIN" compute -k "$k" -o /dev/null -M "$ms" "$san" \
        || { warn "fmsi: kmercamel compute failed"; return 1; }
    [[ -s "$ms" ]] || { warn "fmsi: empty masked superstring"; return 1; }
    csec="$RUN_SEC"; crss="$RUN_RSS_KB"
    # 2) MBWT index over the masked superstring (writes ms.fa.fmsi.* next to ms.fa).
    # -x: skip the optional kLCP array. The kLCP only accelerates *streamed* queries (via the
    # query-side -S flag) at a large space cost; FMSI's headline compact space is without it, so
    # building it -x keeps the bits/k-mer comparison fair (cf. the SBWT-variant note). For an
    # accelerated streaming config instead, drop -x here and add -S to query_fmsi (bigger index).
    run_timed "$FMSI_BIN" index -x -k "$k" "$ms" || { warn "fmsi: index failed"; return 1; }
    RUN_SEC=$(awk -v a="$csec" -v b="$RUN_SEC" 'BEGIN{printf "%.3f", a+b}')   # kmercamel + fmsi index
    (( crss > RUN_RSS_KB )) && RUN_RSS_KB="$crss"
    IDX_PATH="$ms"   # query takes the ms.fa stem; the real data lives in ms.fa.fmsi.*
    # Index size = the persistent .fmsi.* component files. The masked superstring ms.fa is the
    # intermediate representation (analogous to BCALM unitigs for sshash) and is excluded.
    shopt -s nullglob; for f in "$ms".fmsi.*; do bytes=$(( bytes + $(stat -c%s "$f") )); done; shopt -u nullglob
    (( bytes > 0 )) || { warn "fmsi: no .fmsi.* index files produced"; return 1; }
    IDX_FILE_BYTES="$bytes"; IDX_PAYLOAD_BYTES="$bytes"; N_SKMERS=NA
    return 0
}

query_fmsi() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    run_timed_median taskset -c "$cpus" "$FMSI_BIN" query -O -q "$qfa" -k "$k" "$idx"
}
