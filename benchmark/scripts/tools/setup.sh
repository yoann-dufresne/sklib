#!/usr/bin/env bash
# Clone + build the competitor tools for the sklib benchmark into a git-ignored
# source tree, and write their resolved binary paths to tools.env (sourced by
# benchmark/scripts/tools.sh). Idempotent: skips clone/build when already present.
#
#   bash benchmark/scripts/tools/setup.sh [tool ...]      # default: all
#
# Tools: bcalm (unitigs), sshash, sbwt, bqf, cbl. CBL is K-specialized at compile
# time, so its per-k binaries are built on demand by the wrapper, not here.
set -uo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$HERE/../../.." && pwd)"
TOOLS_SRC="${TOOLS_SRC:-$REPO_ROOT/benchmark/data/tools_src}"
ENV_FILE="$TOOLS_SRC/tools.env"
LOG="$TOOLS_SRC/setup.log"
JOBS="${JOBS:-$(nproc)}"
mkdir -p "$TOOLS_SRC"
: > "$LOG"

say() { printf '[setup] %s\n' "$*" | tee -a "$LOG" >&2; }
run() { echo "+ $*" >>"$LOG"; "$@" >>"$LOG" 2>&1; }

clone() {  # clone <dir> <url>
    local dir="$1" url="$2"
    if [[ -d "$TOOLS_SRC/$dir/.git" ]]; then say "$dir: already cloned"; return 0; fi
    say "$dir: cloning $url"
    run git -C "$TOOLS_SRC" clone --recursive --depth 1 --shallow-submodules -j"$JOBS" "$url" "$dir"
}

record() { echo "$1=\"$2\"" >> "$ENV_FILE"; say "  -> $1=$2"; }

# ---------------------------------------------------------------- bcalm
setup_bcalm() {
    clone bcalm https://github.com/GATB/bcalm || return 1
    local bin="$TOOLS_SRC/bcalm/build/bcalm"
    if [[ ! -x "$bin" ]]; then
        say "bcalm: building"
        run mkdir -p "$TOOLS_SRC/bcalm/build"
        ( cd "$TOOLS_SRC/bcalm/build" && run cmake -DCMAKE_BUILD_TYPE=Release .. && run make -j"$JOBS" )
    fi
    [[ -x "$bin" ]] && record BCALM_BIN "$bin" || say "bcalm: BUILD FAILED (see log)"
}

# ---------------------------------------------------------------- sshash
setup_sshash() {
    clone sshash https://github.com/jermp/sshash || return 1
    local bdir="$TOOLS_SRC/sshash/build" bin
    if ! ls "$bdir"/sshash "$bdir"/bin/sshash >/dev/null 2>&1; then
        say "sshash: building"
        run mkdir -p "$bdir"
        ( cd "$bdir" && run cmake -DCMAKE_BUILD_TYPE=Release .. && run make -j"$JOBS" )
    fi
    bin="$(find "$bdir" -maxdepth 2 -type f -name sshash -perm -u+x 2>/dev/null | head -1)"
    [[ -n "$bin" ]] && record SSHASH_BIN "$bin" || say "sshash: BUILD FAILED (see log)"
}

# ---------------------------------------------------------------- sbwt
setup_sbwt() {
    clone SBWT https://github.com/algbio/SBWT || return 1
    local bdir="$TOOLS_SRC/SBWT/build" bin
    if ! find "$bdir" -name sbwt -perm -u+x >/dev/null 2>&1; then
        say "sbwt: building (MAX_KMER_LENGTH=32)"
        rm -rf "$bdir"; mkdir -p "$bdir"   # clean dir: stale CMake state can fail the build
        ( cd "$bdir" && run cmake -DCMAKE_BUILD_TYPE=Release -DMAX_KMER_LENGTH=32 .. && run make -j"$JOBS" )
    fi
    bin="$(find "$bdir" -maxdepth 3 -type f -name sbwt -perm -u+x 2>/dev/null | head -1)"
    [[ -n "$bin" ]] && record SBWT_BIN "$bin" || say "sbwt: BUILD FAILED (see log)"
}

# ---------------------------------------------------------------- bqf (clone only)
setup_bqf() {
    clone bqf https://github.com/vicLeva/bqf || return 1
    # BQF is specialized at compile time by (k,z) (-DBQF_INDEX_K/-DBQF_INDEX_Z), so the
    # wrapper compiles a per-(k,z) binary on demand; here we just record the source tree.
    record BQF_SRC "$TOOLS_SRC/bqf"
}

# ---------------------------------------------------------------- cbl (clone only)
setup_cbl() {
    clone CBL https://github.com/imartayan/CBL || return 1
    record CBL_SRC "$TOOLS_SRC/CBL"   # K-specialized; wrapper builds per-k binaries
}

# ---------------------------------------------------------------- fmsi (+ kmercamel)
setup_fmsi() {
    # FMSI indexes a masked superstring computed by KmerCamel (separate repo, NOT a submodule).
    clone fmsi      https://github.com/OndrejSladky/fmsi      || return 1
    clone kmercamel https://github.com/OndrejSladky/kmercamel || return 1
    local fbin="$TOOLS_SRC/fmsi/fmsi" kbin="$TOOLS_SRC/kmercamel/kmercamel"
    [[ -x "$fbin" ]] || { say "fmsi: building";      ( cd "$TOOLS_SRC/fmsi"      && run make -j"$JOBS" ); }
    [[ -x "$kbin" ]] || { say "kmercamel: building"; ( cd "$TOOLS_SRC/kmercamel" && run make -j"$JOBS" ); }
    [[ -x "$fbin" ]] && record FMSI_BIN      "$fbin" || say "fmsi: BUILD FAILED (see log)"
    # KmerCamel needs GLPK (glpk.h / -lglpk) for the mask-optimization ILP; on Debian/Ubuntu
    # this is 'apt install libglpk-dev' (needs sudo). Its Makefile's /opt/homebrew paths are
    # harmless on Linux (gcc ignores missing -I/-L dirs).
    [[ -x "$kbin" ]] && record KMERCAMEL_BIN "$kbin" \
        || say "kmercamel: BUILD FAILED -- install GLPK ('apt install libglpk-dev'), then re-run; see log"
}

TOOLS=("$@"); [[ ${#TOOLS[@]} -eq 0 ]] && TOOLS=(bcalm sshash sbwt bqf cbl fmsi)
: > "$ENV_FILE"
for t in "${TOOLS[@]}"; do
    say "=== $t ==="
    "setup_$t" || say "$t: setup error"
done
say "done. binary paths in $ENV_FILE :"
cat "$ENV_FILE" >&2
