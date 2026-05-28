#!/usr/bin/env bash
# BQF adapter — Backpack Quotient Filter (github.com/vicLeva/bqf): compact *approximate*
# membership with abundances. APPROXIMATE: it can report false positives (a space/accuracy
# trade-off), unlike the exact tools — flag this when reading the Pareto figure.
#
# BQF inserts s-mers (s = k - z) and virtualizes k-mers (fimpera). It is specialized at
# COMPILE TIME by (k,z), so the wrapper builds a per-(k,z) binary on demand and caches it.
# Input to `bqf build` is a KMC dump of counted s-mers.
#
# CLI (verified): bqf build -q <quotient> -c <counter_bits> -k <k> -z <z> -i <smers.txt> -o <bqf>
#                 bqf query -b <bqf> -i <queries> -o <out>
# Construction cost = KMC s-mer counting + bqf build (both serial, like sklib/sshash).
# Competitor => uses_m_bqf stays 1.

# BQF indexes s-mers (s = k - z) and virtualizes k-mers. We hold the *indexed word
# size* s roughly constant across k (BQF_S, capped at k-1), so the space/FP trade-off
# is comparable across the k-sweep: z = k - s. (z=5 fixed degenerates at small k.)
BQF_S="${BQF_S:-18}"
BQF_COUNTER_BITS="${BQF_COUNTER_BITS:-5}"

available_bqf() { [[ -n "${BQF_SRC:-}" && -d "${BQF_SRC:-}/.git" ]] && command -v kmc >/dev/null 2>&1; }
version_bqf() { echo "bqf-s$BQF_S(approx)"; }

# Build (or reuse) a BQF binary compiled for this (k,z). Echoes its path.
_bqf_binary() {
    local k="$1" z="$2" bdir="$BQF_SRC/build_k${k}_z${z}" bin="$BQF_SRC/build_k${k}_z${z}/bin/bqf"
    if [[ ! -x "$bin" ]]; then
        warn "bqf: compiling binary for k=$k z=$z (one-off)"
        mkdir -p "$bdir"
        ( cd "$bdir" && cmake -DCMAKE_BUILD_TYPE=Release -DBQF_INDEX_K="$k" -DBQF_INDEX_Z="$z" .. \
            && make -j"$(nproc)" bqf ) >"$bdir/build.log" 2>&1
    fi
    [[ -x "$bin" ]] && echo "$bin"
}

construct_bqf() {
    local san="$1" k="$2" m="$3" wd="$4"
    # hold s ~= BQF_S (capped at k-1); z = k - s
    local s=$(( BQF_S <= k-1 ? BQF_S : k-1 )) z
    z=$(( k - s ))
    if (( s < 1 || z < 1 )); then warn "bqf: bad s=$s z=$z for k=$k"; return 1; fi
    local bin; bin="$(_bqf_binary "$k" "$z")" || { warn "bqf: per-(k,z) build failed (see build.log)"; return 1; }

    # 1) count s-mers with KMC (serial) and dump as the 'counted s-mers' BQF input
    local db="$wd/kmc_s" smers="$wd/smers.txt"
    run_timed kmc -k"$s" -ci1 -t1 -fm "$san" "$db" "$wd" || { warn "bqf: KMC s-mer count failed"; return 1; }
    local ksec="$RUN_SEC" krss="$RUN_RSS_KB"
    kmc_tools transform "$db" dump "$smers" >/dev/null 2>&1 || { warn "bqf: KMC dump failed"; return 1; }
    rm -f "$db".kmc_pre "$db".kmc_suf
    local n; n=$(wc -l < "$smers"); (( n < 1 )) && { warn "bqf: no s-mers"; return 1; }
    local q; q=$(awk -v n="$n" 'BEGIN{q=3; while((2^q)<=n) q++; print q}')   # smallest q with 2^q > n

    # 2) build the filter
    IDX_PATH="$wd/index.bqf"
    run_timed "$bin" build -q "$q" -c "$BQF_COUNTER_BITS" -k "$k" -z "$z" -i "$smers" -o "$IDX_PATH" \
        || { warn "bqf: build failed"; return 1; }
    [[ -s "$IDX_PATH" ]] || return 1
    RUN_SEC=$(awk -v a="$ksec" -v b="$RUN_SEC" 'BEGIN{printf "%.3f", a+b}')   # KMC + bqf build
    (( krss > RUN_RSS_KB )) && RUN_RSS_KB="$krss"
    BQF_BIN_CACHED="$bin"
    IDX_FILE_BYTES=$(stat -c%s "$IDX_PATH"); IDX_PAYLOAD_BYTES="$IDX_FILE_BYTES"; N_SKMERS=NA
    return 0
}

query_bqf() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    run_timed_median taskset -c "$cpus" "${BQF_BIN_CACHED:?bqf binary not set}" query -b "$idx" -i "$qfa" -o /dev/null
}
