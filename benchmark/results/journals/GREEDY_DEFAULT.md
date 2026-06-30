# Construction default: greedy_chaining (was colinear_chaining)

**Question.** Construction's per-column reconciliation used the Fenwick `colinear_chaining` (pinned for
byte-reproducibility), while the set-op path uses the faster patience-sort `greedy_chaining`. Both pick
a *maximum* chain, so the represented k-mer set is identical — but does greedy yield the **same number
of compactions** (record count) as colinear across diverse large genomes? If yes, make greedy the
default (it is faster). If no, keep colinear.

**TL;DR.** **Yes — identical on all 20 genomes >300 MB tested** (and 6 smaller ones), with the **same
k-mer set** (setop `diff = 0` both ways) and **~3–7 % faster phase 2**. Greedy is now the construction
default (`v0.11.0`). The byte *packing* differs, so indexes built from v0.11 on are **not byte-identical**
to older ones (same content, queryable-identical); `SKLIB_CONSTRUCT_COLINEAR=1` reproduces the old
colinear bytes exactly.

## Why the count is invariant (not just empirically)
Each column pair's candidate overlaps are computed from the **raw bucket enumeration** (fixed by the
k-mer set), not from the partially-merged list. Both algorithms select a **maximum** chain (same length)
per column ⇒ same total merges ⇒ `#records = #kmers − #merges` is the same. Only *which* overlaps are
chosen (tie-breaking among co-optimal chains) differs ⇒ different packing, identical count. Confirmed
end-to-end by setop set-difference = 0.

## Method
`benchmark/scripts/verify/greedy_chaining_verif.sh`: for each genome, build with colinear (default at the time)
and greedy (`SKLIB_CONSTRUCT_GREEDY` during the experiment), compare the VSKMER_4 record count, delete
genome+index (peak disk ≈ one genome). k=31/m=15, `-t16`, clang-18 Release+LTO. Genomes fetched from
UCSC goldenPath (>300 MB unpacked).

## Result — 20/20 genomes >300 MB, record count colinear == greedy

| genome | MB | records (both) | genome | MB | records (both) |
|---|--:|--:|---|--:|--:|
| chm13 (human) | 3169 | 385 553 171 | felCat9 (cat) | 2572 | 336 745 643 |
| mm39 (mouse) | 2782 | 325 922 218 | equCab3 (horse) | 2557 | 336 271 159 |
| rn7 (rat) | 2700 | 332 414 150 | oviAri4 (sheep) | 2667 | 331 446 386 |
| canFam6 (dog) | 2359 | 319 481 550 | monDom5 (opossum) | 3677 | 476 150 976 |
| bosTau9 (cow) | 2770 | 331 814 197 | ornAna2 (platypus) | 2041 | 223 727 686 |
| susScr11 (pig) | 2551 | 330 336 150 | gasAcu1 (stickleback) | 472 | 50 275 498 |
| galGal6 (chicken) | 1086 | 135 268 168 | fr3 (fugu) | 399 | 39 212 794 |
| danRer11 (zebrafish) | 1712 | 165 705 936 | oryLat2 (medaka) | 886 | 82 862 584 |
| xenTro10 (xenopus) | 1480 | 153 629 523 | taeGut2 (zebra finch) | 1257 | 148 169 720 |
| rheMac10 (macaque) | 3030 | 384 413 877 | panTro6 (chimp) | 3111 | 385 502 079 |

**0 mismatches.** Also identical on ecoli/yeast/celegans/chr21/chr1 (incl. chr1 at k63/m31).

**Speed (phase 2, greedy vs colinear, interleaved A/B):** celegans k31 ~−6 %, chr1 k31 ~−7 %,
chr1 k63 ~−3.5 % — greedy faster in nearly every rep.

## Notes
- `ctest` 213/213 green with greedy default (no test pins the compaction packing; the byte-exact
  `skmerator_digest` pin is the *producer* stream, which is upstream of the chaining and unchanged).
- Greedy construction stays **deterministic across `-t`** (chr21 `-t1` == `-t8`).
- Escape hatch `SKLIB_CONSTRUCT_COLINEAR=1` reproduces the pre-0.11 colinear bytes (sha256-verified).
