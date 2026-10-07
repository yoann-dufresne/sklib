#!/usr/bin/env bash
# Set operations on real read-set pairs, sklib vs KMC — the whole campaign, in order.
# Resumable: re-run this script and every step skips what is already measured.
#
#   setsid nohup systemd-inhibit --what=sleep:idle:handle-lid-switch --who=sklib-bench \
#       --why="setop campaign" --mode=block bash drivers/run_campaign.sh > logs/campaign.log 2>&1 &
#
# Pause (the laptop is needed): bash benchmark/scripts/setop_pairs.sh pause    … then: resume
# A paused measurement is thrown away and redone; nothing measured while paused is kept.
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RUN="$(dirname "$HERE")"
REPO="$(cd -- "$RUN/../../../.." && pwd)"
export RUN_DIR="$RUN"
export SSKM_BIN="${SSKM_BIN:-$REPO/build-bench/bin/sskm}"     # Release, clang-18, built from 0afc065
SP="$REPO/benchmark/scripts/setop_pairs.sh"
say() { printf '[campaign %s] %s\n' "$(date '+%m-%d %H:%M:%S')" "$*"; }

# 0. data: let a download already in flight finish (two fetchers would append to the same
#    parts), then resume/verify every download and build the six sanitized FASTAs.
if [[ -f "$RUN/logs/fetch.pid" ]] && kill -0 "$(cat "$RUN/logs/fetch.pid")" 2>/dev/null; then
    say "waiting for the running download (pid $(cat "$RUN/logs/fetch.pid"))"
    while kill -0 "$(cat "$RUN/logs/fetch.pid")" 2>/dev/null; do sleep 60; done
fi
say "fetch + sanitize"
bash "$REPO/benchmark/scripts/fetch_reads.sh" >> "$RUN/logs/fetch.log" 2>&1 || { say "fetch failed, see logs/fetch.log"; exit 1; }

# 1. drift anchor (start): the deck's chr1 point (k=31, J=0.5 mutant, union, t=1/t=8) under
#    this campaign's protocol, so the new numbers can be laid next to the published ones.
say "anchor (start)"
CSV="$RUN/data/anchor_start_runs.csv" CCSV="$RUN/data/anchor_construct_runs.csv" SCSV="$RUN/data/anchor_pair_sizes.csv" \
PAIRS="chr1:chr1_mutJ05" KM="31,15" OPS="union" THREADS="1 8" REPS=3 \
    bash "$SP" >> "$RUN/logs/setop_pairs.log" 2>&1

# 2. the campaign: core (union + joint at t=1/t=8, sklib --sizes) then E1 (∩, A∖B, B∖A at t=8),
#    k=31 for the three pairs first, then k=63.
for km in "31,15" "63,31"; do
    say "pairs, k/m=$km"
    PAIRS="hg002_A:hg002_B gut_A:gut_B ocean_SRF:ocean_DCM" KM="$km" \
    OPS="union joint sizes inter@8 diffab@8 diffba@8" THREADS="1 8" REPS=3 \
        bash "$SP" >> "$RUN/logs/setop_pairs.log" 2>&1 || say "setop_pairs exited non-zero for k/m=$km (see log)"
done

# 3. drift anchor (end)
say "anchor (end)"
CSV="$RUN/data/anchor_end_runs.csv" CCSV="$RUN/data/anchor_construct_runs.csv" SCSV="$RUN/data/anchor_pair_sizes.csv" \
PAIRS="chr1:chr1_mutJ05" KM="31,15" OPS="union" THREADS="1 8" REPS=3 \
    bash "$SP" >> "$RUN/logs/setop_pairs.log" 2>&1
say "campaign done"
