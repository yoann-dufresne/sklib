#!/usr/bin/env bash
# Per-tool adapters for scripts/bench/bench.sh. Each tool <T> defines:
#   available_<T>            -> 0 if the tool can run, else non-zero (bench skips it)
#   version_<T>              -> prints a version string
#   construct_<T> san k m wd -> builds the index; on success sets globals:
#                                 IDX_PATH, IDX_FILE_BYTES, IDX_PAYLOAD_BYTES, N_SKMERS
#                               (RUN_SEC / RUN_RSS_KB come from run_timed); returns status
#   query_<T> idx qfa k m cpus -> runs the query pinned to <cpus> over RUN_REPS reps;
#                                 sets RUN_SEC (median) / RUN_RSS_KB; returns status
#
# Only sklib is implemented; competitors are stubs disabled via available_* until
# their wrappers land (scripts/bench/tools/). Sourced -- not executed.

SSKM_BIN="${SSKM_BIN:-$BENCH_REPO_ROOT/build-bench/bin/sskm}"  # Release; build/ is a DEBUG build
SCRIPT_DIR_TOOLS="$BENCH_REPO_ROOT/scripts/bench/tools"

# uses_m_<tool>: 0 (true) if the tool's index depends on sklib's minimizer length m,
# so it must be (re)built for every m in the grid. Competitors don't use sklib's m, so
# they are built once per (dataset,k) and reused -- bench.sh dedups via this predicate.

# ---- sklib (the structure under test) -------------------------------------
available_sklib() { [[ -x "$SSKM_BIN" ]]; }
uses_m_sklib() { return 0; }

version_sklib() {
    local v
    v="$("$SSKM_BIN" --version 2>/dev/null | head -1)"; [[ -n "$v" ]] || v="unknown"
    echo "sklib-$v"   # sklib is fully sequential; no TBB tag (it is not used)
}

construct_sklib() {
    local san="$1" k="$2" m="$3" wd="$4"
    IDX_PATH="$wd/index.sskm"
    # Optional bucket count (BUCKETS env, empty => tool default 4096). Used by the bucket sweep.
    local bucket_args=()
    [[ -n "${BUCKETS:-}" ]] && bucket_args=(--buckets "$BUCKETS")
    # Optional construction thread count (CONSTRUCT_THREADS env, empty => binary default). The
    # per-bucket compaction is parallelized; the index is byte-identical for any thread count.
    local thread_args=()
    [[ -n "${CONSTRUCT_THREADS:-}" ]] && thread_args=(-t "$CONSTRUCT_THREADS")
    run_timed "$SSKM_BIN" construct -k "$k" -m "$m" -f "$san" -o "$IDX_PATH" "${bucket_args[@]}" "${thread_args[@]}" || return 1
    [[ -s "$IDX_PATH" ]] || return 1
    IDX_FILE_BYTES=$(stat -c%s "$IDX_PATH")
    # Payload = file minus the real header. V3 header is 40 + 16*n_buckets bytes (the per-bucket
    # directory), not the old fixed 32, so subtract the size the format actually reports.
    local hdr
    hdr=$(python3 "$BENCH_HELPER" binheadersize "$IDX_PATH" 2>/dev/null || echo 32)
    IDX_PAYLOAD_BYTES=$(( IDX_FILE_BYTES - hdr ))
    N_SKMERS=$(python3 "$BENCH_HELPER" bincount "$IDX_PATH" 2>/dev/null || echo NA)
    return 0
}

query_sklib() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    # Pass -t = number of pinned cores so the "threads" column is honest: cpus="0" -> -t 1
    # (true single-thread sequential reader), cpus="0-21" -> -t 22 (parallel file query).
    local ncpu; if [[ "$cpus" == *-* ]]; then ncpu=$(( ${cpus#*-} - ${cpus%%-*} + 1 )); else ncpu=1; fi
    run_timed_median taskset -c "$cpus" "$SSKM_BIN" query -t "$ncpu" -l "$idx" -i "$qfa" -o /dev/null
}

# ---- competitor stubs (filled by scripts/bench/tools/*.sh, task #6) --------
_stub_unavail() { return 1; }
_stub_version() { echo "n/a"; }
_stub_construct() { warn "construct: '$1' wrapper not implemented yet"; return 1; }
_stub_query()    { return 1; }

for _t in sshash bqf cbl sbwt bcalm fmsi; do
    eval "available_$_t() { _stub_unavail; }"
    eval "version_$_t()   { _stub_version; }"
    eval "construct_$_t() { _stub_construct $_t; }"
    eval "query_$_t()     { _stub_query; }"
    eval "uses_m_$_t()    { return 1; }"   # competitors ignore sklib's m
done
unset _t

# Real wrappers (built by scripts/bench/tools/setup.sh) override the stubs above.
_TOOLS_ENV="$BENCH_REPO_ROOT/scripts/out/bench/tools_src/tools.env"
[[ -f "$_TOOLS_ENV" ]] && source "$_TOOLS_ENV"
for _w in "$SCRIPT_DIR_TOOLS"/*.sh; do
    [[ -f "$_w" && "$(basename "$_w")" != "setup.sh" ]] && source "$_w"
done
unset _w
