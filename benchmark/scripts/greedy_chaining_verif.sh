#!/bin/bash
# Greedy-vs-colinear compaction-size test on 20 distinct genomes >300MB.
# For each genome: build with colinear_chaining (default) and greedy_chaining
# (SKLIB_CONSTRUCT_GREEDY=1), compare the VSKMER_4 record count. Deletes each
# genome+index after testing so peak disk stays ~one genome.
set -u
ROOT=/home/yoann/bioinfo/sklib
SSKM=$ROOT/build-rel-clang/bin/sskm
WD=$ROOT/benchmark/data/greedy_test
mkdir -p "$WD"
LOG="$WD/results.tsv"
DL="$WD/dl.log"
: > "$DL"
printf '#%-12s %8s %4s %4s %14s %14s %s\n' genome size_MB k m colinear greedy match > "$LOG"
cnt(){ od -An -t u8 -j 24 -N 8 "$1" | tr -d ' '; }   # record count @ header offset 24
K=31; M=15
n=0
test_one(){ # name fasta
  local name="$1" fa="$2"
  local szMB=$(( $(stat -c%s "$fa") / 1000000 ))
  if ! $SSKM construct -k$K -m$M -f "$fa" -o "$WD/col.sskm" -t16 2>/dev/null; then
    printf '%-13s %8s %4s %4s %14s %14s %s\n' "$name" "$szMB" "$K" "$M" BUILD_FAIL - - >> "$LOG"; return 1; fi
  if ! SKLIB_CONSTRUCT_GREEDY=1 $SSKM construct -k$K -m$M -f "$fa" -o "$WD/grd.sskm" -t16 2>/dev/null; then
    printf '%-13s %8s %4s %4s %14s %14s %s\n' "$name" "$szMB" "$K" "$M" - BUILD_FAIL - >> "$LOG"; return 1; fi
  local c g mt
  c=$(cnt "$WD/col.sskm"); g=$(cnt "$WD/grd.sskm")
  if [ "$c" = "$g" ]; then mt=YES; else mt="NO(diff=$((g-c)))"; fi
  printf '%-13s %8s %4s %4s %14s %14s %s\n' "$name" "$szMB" "$K" "$M" "$c" "$g" "$mt" >> "$LOG"
  rm -f "$WD/col.sskm" "$WD/grd.sskm"
  n=$((n+1))
}
# 1) chm13 (already present, 3GB)
[ -f "$ROOT/benchmark/data/genomes/chm13.sanitized.fa" ] && test_one chm13 "$ROOT/benchmark/data/genomes/chm13.sanitized.fa"
# 2) download diverse UCSC genomes >300MB (skip download/too-small failures), stop at 20
DBS="mm39 rn7 canFam6 bosTau9 susScr11 galGal6 danRer11 xenTro10 rheMac10 panTro6 felCat9 equCab3 oviAri4 monDom5 ornAna2 gasAcu1 fr3 oryLat2 taeGut2 anoCar2 gorGor6 ponAbe3 calJac4 loxAfr3 nomLeu3 micMur2"
for db in $DBS; do
  [ "$n" -ge 20 ] && break
  url="https://hgdownload.soe.ucsc.edu/goldenPath/$db/bigZips/$db.fa.gz"
  echo "[$(date +%H:%M:%S)] downloading $db ..." >> "$DL"
  if ! timeout 1800 curl -sf "$url" -o "$WD/$db.fa.gz"; then echo "  DL_FAIL $db" >> "$DL"; rm -f "$WD/$db.fa.gz"; continue; fi
  if ! gunzip -f "$WD/$db.fa.gz"; then echo "  GUNZIP_FAIL $db" >> "$DL"; rm -f "$WD/$db.fa.gz" "$WD/$db.fa"; continue; fi
  sz=$(stat -c%s "$WD/$db.fa" 2>/dev/null || echo 0)
  if [ "$sz" -lt 300000000 ]; then echo "  TOO_SMALL $db ($((sz/1000000))MB)" >> "$DL"; rm -f "$WD/$db.fa"; continue; fi
  echo "[$(date +%H:%M:%S)] testing $db ($((sz/1000000))MB) ..." >> "$DL"
  test_one "$db" "$WD/$db.fa"
  rm -f "$WD/$db.fa"
done
printf 'DONE: %s genomes tested\n' "$n" >> "$LOG"
