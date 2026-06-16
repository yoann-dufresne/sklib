#!/usr/bin/env bash
# sbwt-rs adapter — the Rust rewrite of SBWT (github.com/jnalanko/sbwt-rs-cli) which,
# unlike the C++ algbio/SBWT, exposes native k-mer SET OPERATIONS on prebuilt indexes
# (merge=union, intersect, difference=A\B) AND honours -t/--threads, so it can be swept
# across the thread grid like sklib/kmc. Added as a SEPARATE tool ("sbwtrs") so its rows
# never collide with the C++ "sbwt" rows; this also lets us compare the two on the plain
# construct/query operations to see if the Rust version is faster.
#
# CLI (verified against src/main.rs):
#   global:  sbwt -t <N> [-v] <subcommand> ...        # -t/--threads is GLOBAL, before the subcommand
#   build:   build -i <fa> -o <prefix> -k <k> [-r add-revcomp] [--temp-dir DIR]  # writes <prefix>.sbwt
#   query:   lookup -i <idx.sbwt> -q <q.fa> -o <out> [-b membership-only]
#   set-op:  intersect|merge|difference <sbwt1.sbwt> <sbwt2.sbwt> -o <out>        # positional inputs
#
# -r/--add-revcomp gives canonical membership (a k-mer or its RC counts as present),
# matching the C++ sbwt (--add-reverse-complements), sklib and KMC. No size-only or
# combined (joint) set-op mode. k is not compile-time capped at 32 (the Rust LongKmer
# supports k<=255), so unlike the C++ sbwt it can run the full k grid; a k it cannot
# handle just fails the build and the harness skips that (tool,k). Competitor => uses_m stays 1.

available_sbwtrs() { [[ -n "${SBWTRS_BIN:-}" && -x "${SBWTRS_BIN:-}" ]]; }
version_sbwtrs()   { echo "sbwt-rs-${SBWTRS_VER:-0.2.0}"; }
uses_m_sbwtrs()    { return 1; }
thread_list_sbwtrs() { echo "$THREADS"; }   # honours -t -> sweep the thread grid

can_construct_sbwtrs() { available_sbwtrs; }
can_query_sbwtrs()     { available_sbwtrs; }
can_setop_sbwtrs()     { available_sbwtrs; }
setop_has_size_sbwtrs()  { return 1; }   # no cardinality-only mode (must materialize)
setop_has_joint_sbwtrs() { return 1; }   # no combined single-pass mode

construct_sbwtrs() {
    local san="$1" k="$2" m="$3" wd="$4"
    IDX_PATH="$wd/index.sbwt"
    run_timed "$SBWTRS_BIN" -t "${CONSTRUCT_THREADS:-1}" build \
        -i "$san" -o "$wd/index" -k "$k" -r --temp-dir "$wd" \
        || { warn "sbwtrs: build failed (k=$k unsupported?)"; return 1; }
    [[ -s "$IDX_PATH" ]] || return 1
    IDX_FILE_BYTES=$(stat -c%s "$IDX_PATH"); IDX_PAYLOAD_BYTES="$IDX_FILE_BYTES"; N_SKMERS=NA
    return 0
}

query_sbwtrs() {
    local idx="$1" qfa="$2" k="$3" m="$4" cpus="$5"
    local ncpu; if [[ "$cpus" == *-* ]]; then ncpu=$(( ${cpus#*-} - ${cpus%%-*} + 1 )); else ncpu=1; fi
    run_timed_median taskset -c "$cpus" "$SBWTRS_BIN" -t "$ncpu" \
        lookup -i "$idx" -q "$qfa" -o /dev/null -b
}

setop_op_sbwtrs() {                                # A B op mode out th
    local A="$1" B="$2" op="$3" mode="$4" out="$5" th="$6"
    [[ "$mode" == size ]] && return 1              # materialize only
    case "$op" in
      inter)  run_timed_median "$SBWTRS_BIN" -t "$th" intersect  "$A" "$B" -o "$out" ;;
      union)  run_timed_median "$SBWTRS_BIN" -t "$th" merge      "$A" "$B" -o "$out" ;;
      diffab) run_timed_median "$SBWTRS_BIN" -t "$th" difference "$A" "$B" -o "$out" ;;
      diffba) run_timed_median "$SBWTRS_BIN" -t "$th" difference "$B" "$A" -o "$out" ;;
    esac
}
