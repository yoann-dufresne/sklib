#!/usr/bin/env bash
# CBL adapter — Conway-Bromage-Lyndon (github.com/imartayan/CBL): exact, fully dynamic
# k-mer set with high locality, and native union/intersection/difference. Closest in
# spirit to sklib (sorted exact set + planned set ops), so the most relevant exact rival.
#
# CBL specializes K at COMPILE TIME via the `K` env var (build.rs) and requires:
#   - Rust nightly (rustup install nightly)         [present]
#   - apt: libclang-dev libstdc++-12-dev            [USER must install; needs sudo]
# K MUST BE ODD and in [1,59]. The wrapper builds+caches a per-K binary on demand.
#
# CLI (from examples/cbl.rs): build <in> [-o idx] [--canonical] ; query <idx> <in>
#                             merge/inter/diff <a> <b> ...   (set ops, for future use)
# Competitor => uses_m_cbl stays 1.

available_cbl() {
    [[ -n "${CBL_SRC:-}" && -d "${CBL_SRC:-}" ]] && command -v cargo >/dev/null 2>&1 \
        && rustup toolchain list 2>/dev/null | grep -q nightly
}
version_cbl() { echo "cbl"; }

# Build (or reuse) a CBL binary compiled for this odd K. Echoes its cached path.
_cbl_binary() {
    local k="$1" cache="$CBL_SRC/bin/cbl_k$k"
    if [[ ! -x "$cache" ]]; then
        (( k % 2 == 1 )) || { warn "cbl: K must be odd (k=$k)"; return 1; }
        warn "cbl: building K=$k (one-off; needs libclang-dev + libstdc++-*-dev)"
        mkdir -p "$CBL_SRC/bin"
        # bindgen/libclang doesn't auto-find the GCC C++ stdlib headers (<algorithm>,
        # bits/c++config.h) -> point it at the default g++'s include dirs explicitly.
        local gv cxxinc
        gv="$(g++ -dumpversion 2>/dev/null | cut -d. -f1)"
        cxxinc="-I/usr/include/c++/$gv -I/usr/include/x86_64-linux-gnu/c++/$gv"
        ( cd "$CBL_SRC" && K="$k" \
            BINDGEN_EXTRA_CLANG_ARGS="${BINDGEN_EXTRA_CLANG_ARGS:-} $cxxinc" \
            cargo +nightly build --release --example cbl ) \
            >"$CBL_SRC/bin/build_k$k.log" 2>&1 \
            || { warn "cbl: cargo build failed (see bin/build_k$k.log)"; return 1; }
        cp "$CBL_SRC/target/release/examples/cbl" "$cache" || return 1
    fi
    echo "$cache"
}

construct_cbl() {
    local san="$1" k="$2" m="$3" wd="$4" bin
    (( k % 2 == 1 )) || { warn "cbl: K must be odd (k=$k), skipping"; return 1; }
    bin="$(_cbl_binary "$k")" || return 1
    IDX_PATH="$wd/index.cbl"
    # CBL aborts on any sequence shorter than K ("Sequence size (1) is smaller than K"),
    # e.g. 1-bp fragments left between N-runs by sanitize on human chromosomes. Drop sub-K
    # records first; they hold no k-mers, so the indexed k-mer set is unchanged.
    local cbl_in="$wd/cbl_input.fa"
    python3 "$BENCH_HELPER" dropshort "$san" "$k" > "$cbl_in" || { warn "cbl: dropshort failed"; return 1; }
    run_timed "$bin" build "$cbl_in" -o "$IDX_PATH" --canonical || { warn "cbl: build failed"; return 1; }
    [[ -s "$IDX_PATH" ]] || return 1
    CBL_BIN_CACHED="$bin"
    IDX_FILE_BYTES=$(stat -c%s "$IDX_PATH"); IDX_PAYLOAD_BYTES="$IDX_FILE_BYTES"; N_SKMERS=NA
    return 0
}

query_cbl() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    run_timed_median taskset -c "$cpus" "${CBL_BIN_CACHED:?cbl binary not set}" query "$idx" "$qfa"
}
