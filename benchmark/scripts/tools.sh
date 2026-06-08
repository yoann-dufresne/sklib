#!/usr/bin/env bash
# Per-tool adapters for the benchmark harness (construct/query/setop). Each tool <T> defines:
#   available_<T>            -> 0 if the tool can run, else non-zero (bench skips it)
#   version_<T>              -> prints a version string
#   construct_<T> san k m wd -> builds the index; on success sets globals:
#                                 IDX_PATH, IDX_FILE_BYTES, IDX_PAYLOAD_BYTES, N_SKMERS
#                               (RUN_SEC / RUN_RSS_KB come from run_timed); returns status
#   query_<T> idx qfa k m cpus -> runs the query pinned to <cpus> over RUN_REPS reps;
#                                 sets RUN_SEC (median) / RUN_RSS_KB; returns status
#
# Only sklib is implemented; competitors are stubs disabled via available_* until
# their wrappers land (benchmark/scripts/tools/). Sourced -- not executed.

SSKM_BIN="${SSKM_BIN:-$BENCH_REPO_ROOT/build-bench/bin/sskm}"  # Release; build/ is a DEBUG build
SCRIPT_DIR_TOOLS="$BENCH_REPO_ROOT/benchmark/scripts/tools"

# uses_m_<tool>: 0 (true) if the tool's index depends on sklib's minimizer length m,
# so it must be (re)built for every m in the grid. Competitors don't use sklib's m, so
# they are built once per (dataset,k) and reused -- the harness dedups via this predicate.

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

# ---- competitor stubs (filled by benchmark/scripts/tools/*.sh, task #6) --------
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

# Real wrappers (built by benchmark/scripts/tools/setup.sh) override the stubs above.
_TOOLS_ENV="$BENCH_REPO_ROOT/benchmark/data/tools_src/tools.env"
[[ -f "$_TOOLS_ENV" ]] && source "$_TOOLS_ENV"
for _w in "$SCRIPT_DIR_TOOLS"/*.sh; do
    [[ -f "$_w" && "$(basename "$_w")" != "setup.sh" ]] && source "$_w"
done
unset _w

# =====================================================================
#  Experiment-layer adapters (4-experiment harness)
#  Index-handle convention: construct_<t> sets IDX_PATH to the tool's index
#  (sklib .sskm file, KMC db prefix, cbl .cbl file, fmsi .ms.fa) so the set-op
#  adapters can operate on cached handles. All adapters set RUN_SEC/RUN_RSS_KB.
# =====================================================================

# ---- KMC as a construction + set-op tool (membership query is N/A: it is the oracle) ----
available_kmc() { command -v kmc >/dev/null 2>&1 && command -v kmc_tools >/dev/null 2>&1; }
version_kmc()   { echo "kmc"; }
uses_m_kmc()    { return 1; }
construct_kmc() {                                  # san k m wd  -> IDX_PATH = db prefix
    local san="$1" k="$2" wd="$4" db tmp; db="$wd/kmcdb"; tmp="$wd/kt"; mkdir -p "$tmp"; IDX_PATH="$db"
    run_timed kmc -k"$k" -ci1 -t"${CONSTRUCT_THREADS:-1}" -fm "$san" "$db" "$tmp" || return 1
    [[ -s "$db.kmc_pre" ]] || return 1
    IDX_FILE_BYTES=$(( $(stat -c%s "$db.kmc_pre") + $(stat -c%s "$db.kmc_suf") ))
    IDX_PAYLOAD_BYTES="$IDX_FILE_BYTES"; N_SKMERS=NA; return 0
}

# ---- per-tool thread sweep: only sklib + kmc honour -t; the rest run single-thread ----
thread_list_sklib() { echo "$THREADS"; }
thread_list_kmc()   { echo "$THREADS"; }
for _t in sshash sbwt cbl bqf fmsi; do eval "thread_list_$_t() { echo 1; }"; done; unset _t

# ---- per-experiment capability gates (a non-zero construct_<t> additionally skips a (tool,k)) ----
can_construct_sklib(){ available_sklib; }; can_query_sklib(){ available_sklib; }; can_setop_sklib(){ available_sklib; }
can_construct_kmc(){ available_kmc; };      can_query_kmc(){ return 1; };          can_setop_kmc(){ available_kmc; }
for _t in sshash sbwt bqf; do
    eval "can_construct_$_t(){ available_$_t; }"; eval "can_query_$_t(){ available_$_t; }"; eval "can_setop_$_t(){ return 1; }"
done; unset _t
can_construct_cbl(){ available_cbl; };  can_query_cbl(){ available_cbl; };  can_setop_cbl(){ available_cbl; }
can_construct_fmsi(){ available_fmsi; }; can_query_fmsi(){ available_fmsi; }; can_setop_fmsi(){ available_fmsi; }
# size-only set-op (cardinality without materializing): sklib only; KMC/CBL/FMSI must materialize.
setop_has_size_sklib(){ return 0; }; for _t in kmc cbl fmsi; do eval "setop_has_size_$_t(){ return 1; }"; done; unset _t
# combined single-pass set ops: sklib + KMC (kmc_tools simple multi-op); CBL/FMSI single-op only.
setop_has_joint_sklib(){ return 0; }; setop_has_joint_kmc(){ return 0; }
for _t in cbl fmsi; do eval "setop_has_joint_$_t(){ return 1; }"; done; unset _t

# ---- set-op adapters: setop_op_<t> A B op mode out th   (op∈inter|union|diffab|diffba, mode∈materialize|size)
SETOP_SCRATCH="${SETOP_SCRATCH:-$RESULTS/setop_scratch}"
setop_op_sklib() {                                 # A B op mode out th
    local A="$1" B="$2" op="$3" mode="$4" out="$5" th="$6"
    case "$op:$mode" in
      inter:materialize)  run_timed_median "$SSKM_BIN" setop --op intersection -a "$A" -b "$B" -o "$out" -t "$th" ;;
      union:materialize)  run_timed_median "$SSKM_BIN" setop --op union        -a "$A" -b "$B" -o "$out" -t "$th" ;;
      diffab:materialize) run_timed_median "$SSKM_BIN" setop --op diff         -a "$A" -b "$B" -o "$out" -t "$th" ;;
      diffba:materialize) run_timed_median "$SSKM_BIN" setop --op diff         -a "$B" -b "$A" -o "$out" -t "$th" ;;
      inter:size)         run_timed_median "$SSKM_BIN" setop --op intersection_size -a "$A" -b "$B" -t "$th" ;;
      union:size)         run_timed_median "$SSKM_BIN" setop --op union_size        -a "$A" -b "$B" -t "$th" ;;
      diffab:size)        run_timed_median "$SSKM_BIN" setop --op diff_size         -a "$A" -b "$B" -t "$th" ;;
      diffba:size)        run_timed_median "$SSKM_BIN" setop --op diff_size         -a "$B" -b "$A" -t "$th" ;;
    esac
}
setop_joint_sklib() {                              # A B mode outdir th
    local A="$1" B="$2" mode="$3" od="$4" th="$5"
    if [[ "$mode" == materialize ]]; then mkdir -p "$od"
        run_timed_median "$SSKM_BIN" setop -a "$A" -b "$B" -t "$th" \
            --inter-out "$od/i.sskm" --union-out "$od/u.sskm" --diff-ab-out "$od/dab.sskm" --diff-ba-out "$od/dba.sskm"
    else
        run_timed_median "$SSKM_BIN" setop -a "$A" -b "$B" --sizes -t "$th"
    fi
}
# KMC: A,B are db prefixes; no size-only mode (must materialize), so size returns 1.
setop_op_kmc() {                                   # A B op mode out th
    local A="$1" B="$2" op="$3" mode="$4" out="$5" th="$6"
    [[ "$mode" == size ]] && return 1
    case "$op" in
      inter)  run_timed_median kmc_tools -t"$th" simple "$A" "$B" intersect              "$out" ;;
      union)  run_timed_median kmc_tools -t"$th" simple "$A" "$B" union                  "$out" ;;
      diffab) run_timed_median kmc_tools -t"$th" simple "$A" "$B" kmers_subtract         "$out" ;;
      diffba) run_timed_median kmc_tools -t"$th" simple "$A" "$B" reverse_kmers_subtract "$out" ;;
    esac
}
setop_joint_kmc() {                                # A B mode outdir th  (always materializes)
    local A="$1" B="$2" od="$4" th="$5"; mkdir -p "$od"
    run_timed_median kmc_tools -t"$th" simple "$A" "$B" \
        intersect "$od/i" union "$od/u" kmers_subtract "$od/dab" reverse_kmers_subtract "$od/dba"
}
# CBL: per-odd-K binary (cached at $CBL_SRC/bin/cbl_k<K>); A,B are .cbl handles; no size/joint.
_cbl_bin() {                                       # k -> echoes binary path (builds once)
    local k="$1" cache; [[ -n "${CBL_SRC:-}" && -d "${CBL_SRC:-}" ]] || return 1
    (( k % 2 == 1 )) || return 1; cache="$CBL_SRC/bin/cbl_k$k"
    if [[ ! -x "$cache" ]]; then
        mkdir -p "$CBL_SRC/bin"; local gv ci; gv="$(g++ -dumpversion | cut -d. -f1)"
        ci="-I/usr/include/c++/$gv -I/usr/include/x86_64-linux-gnu/c++/$gv"
        ( cd "$CBL_SRC" && K="$k" BINDGEN_EXTRA_CLANG_ARGS="$ci" cargo +nightly build --release --example cbl ) \
            >"$CBL_SRC/bin/build_k$k.log" 2>&1 || return 1
        cp "$CBL_SRC/target/release/examples/cbl" "$cache" || return 1
    fi
    echo "$cache"
}
setop_op_cbl() {                                   # A B op mode out th
    local A="$1" B="$2" op="$3" mode="$4" out="$5" bin; [[ "$mode" == size ]] && return 1
    bin="$(_cbl_bin "$K_CUR")" || return 1
    case "$op" in
      inter)  run_timed_median "$bin" inter "$A" "$B" -o "$out" ;;
      union)  run_timed_median "$bin" merge "$A" "$B" -o "$out" ;;
      diffab) run_timed_median "$bin" diff  "$A" "$B" -o "$out" ;;
      diffba) run_timed_median "$bin" diff  "$B" "$A" -o "$out" ;;
    esac
}
# FMSI: A,B are kmercamel masked-superstring .ms.fa handles; experimental ops, single-thread, no size/joint.
setop_op_fmsi() {                                  # A B op mode out th
    local A="$1" B="$2" op="$3" mode="$4" out="$5"; [[ "$mode" == size ]] && return 1
    [[ -n "${FMSI_BIN:-}" && -x "${FMSI_BIN:-}" ]] || return 1
    case "$op" in
      inter)  run_timed_median timeout "${FMSI_TIMEOUT:-900}" "$FMSI_BIN" inter -p "$A" -p "$B" -r "$out" -k "$K_CUR" ;;
      union)  run_timed_median timeout "${FMSI_TIMEOUT:-900}" "$FMSI_BIN" union -p "$A" -p "$B" -r "$out" -k "$K_CUR" ;;
      diffab) run_timed_median timeout "${FMSI_TIMEOUT:-900}" "$FMSI_BIN" diff  -p "$A" -p "$B" -r "$out" -k "$K_CUR" ;;
      diffba) run_timed_median timeout "${FMSI_TIMEOUT:-900}" "$FMSI_BIN" diff  -p "$B" -p "$A" -r "$out" -k "$K_CUR" ;;
    esac
}
