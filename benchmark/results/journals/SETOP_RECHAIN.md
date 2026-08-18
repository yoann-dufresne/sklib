# What materializing a set operation costs — and what `--no-compact` really buys

**Question.** A materialized set operation ends by re-chaining the kept k-mers into virtual
super-k-mers, with the construction routine. How much of the operation is that re-chaining, and
is `--no-compact` (one record per result k-mer) the shortcut its help text claims — "much faster
but a larger output file"?

**TL;DR.** The re-chaining **is** the dominant phase — **49 % (a difference) to 74 % (a union)** of
the per-bucket wall time at `-t1`, against 22–42 % for the merge that actually computes the result. But **skipping it does not
pay**: `--no-compact` still has to *sort* the kept k-mers, which costs as much as chaining them at
k=31 and **more** at k=63 (chaining is 0.63–0.77× the plain sort there), and it then writes
**8× (k=31) to 15× (k=63)** more bytes. End to end the uncompacted mode is **1.14–2.34× slower**
(median 1.52 over the whole grid); it was faster in **5 of 400** configurations, all of them under
80 ms and by at most 11 %. The `--no-compact` help text is therefore wrong and should be reworded.
Materializing costs **3.3× a cardinality-only pass** (median; 6.9× for union, 2.5× for a difference),
and a materialized result is **exactly as compact as a constructed index**: `A ∪ A` re-chains to A's
own record count on all 8 dataset/k configurations.

## Method
`benchmark/scripts/setop_rechain.sh` — the same grid as `setop.sh` (E. coli, yeast, C. elegans,
chr1 × k=31/m=15, k=63/m=31 × t=1,8 × {∩, ∪, A\B, B\A} × J ∈ {0, .1, .3, .5, .7, .9, 1}, B a mutated
copy of A at the substitution rate giving the target Jaccard), each point the median of `REPS=3`,
in three modes:

| mode | CLI | does |
|---|---|---|
| `size` | `--op *_size` | read + merge only, writes nothing |
| `nocompact` | `--op … --no-compact` | read + merge + `std::sort` of the kept k-mers + write |
| `compact` | `--op …` (default) | read + merge + greedy re-chaining + write |

sklib **v0.15.0**, clang-18 Release (`-O3 -march=native`, LTO), `build-rel-rechain`, on
yoann-Precision-5490 (Core Ultra 7 165H, 22 threads, 62 GiB, Linux 6.17), output on the NVMe.
1344 rows, kept as `benchmark/results/reference/setops_rechain.csv` (working copy under the
git-ignored `benchmark/results/rechain/`); the 144 rows whose result is empty
(the two differences at J=1) are dropped from every aggregate below.
At `-t1` the library's own `SKLIB_UNION_PHASE_TIMING=1` split (read / merge+collect / order / write)
is captured per run, which is where the phase shares come from.

## Result 1 — where the time goes (`-t1`, % of the per-bucket wall, median over ops × J)

| dataset | k | mode | total (s) | read | merge+collect | **order** | write |
|---|--:|---|--:|--:|--:|--:|--:|
| E. coli | 31 | compact | 0.161 | 5.3 % | 35.3 % | **56.9 %** | 2.5 % |
| E. coli | 31 | nocompact | 0.190 | 4.2 % | 34.8 % | 49.5 % | 11.8 % |
| E. coli | 63 | compact | 0.181 | 4.7 % | 30.8 % | **62.0 %** | 2.5 % |
| E. coli | 63 | nocompact | 0.274 | 3.0 % | 26.0 % | 55.4 % | 15.3 % |
| Yeast | 31 | compact | 0.390 | 2.9 % | 36.8 % | **58.1 %** | 2.3 % |
| Yeast | 31 | nocompact | 0.509 | 2.3 % | 34.8 % | 51.8 % | 11.7 % |
| Yeast | 63 | compact | 0.435 | 2.7 % | 33.0 % | **61.5 %** | 2.3 % |
| Yeast | 63 | nocompact | 0.806 | 1.6 % | 30.7 % | 53.8 % | 14.0 % |
| C. elegans | 31 | compact | 3.369 | 1.8 % | 37.0 % | **59.2 %** | 1.9 % |
| C. elegans | 31 | nocompact | 4.291 | 1.5 % | 38.1 % | 51.7 % | 9.1 % |
| C. elegans | 63 | compact | 3.690 | 1.8 % | 35.4 % | **60.7 %** | 2.0 % |
| C. elegans | 63 | nocompact | 6.279 | 1.1 % | 32.1 % | 53.4 % | 12.3 % |
| Chr1 | 31 | compact | 8.024 | 1.6 % | 35.3 % | **61.1 %** | 1.7 % |
| Chr1 | 31 | nocompact | 9.154 | 1.5 % | 35.3 % | 54.1 % | 8.7 % |
| Chr1 | 63 | compact | 8.878 | 1.5 % | 37.6 % | **59.0 %** | 1.7 % |
| Chr1 | 63 | nocompact | 15.045 | 0.9 % | 28.0 % | 59.4 % | 12.1 % |

Over the 200 compacted `-t1` points the order phase is 27–81 % of the wall, median **60.6 %**.

Broken down by operation instead (medians over the four datasets × J, which agree within two points
of each other), the share tracks how much of the operands the result keeps — a union keeps
everything and leaves the merge only 22 %:

| operation | k | read | merge | **order** | write | sort/chain | write nc/c |
|---|--:|--:|--:|--:|--:|--:|--:|
| A ∩ B | 31 | 2.7 % | 32.2 % | **62.7 %** | 2.1 % | 1.01× | 5.4× |
| A ∩ B | 63 | 2.6 % | 29.8 % | **64.2 %** | 2.2 % | 1.51× | 10.1× |
| A ∪ B | 31 | 1.6 % | 22.2 % | **74.0 %** | 2.4 % | 1.01× | 6.0× |
| A ∪ B | 63 | 1.6 % | 22.3 % | **73.3 %** | 2.3 % | 1.54× | 11.9× |
| A \ B | 31 | 3.8 % | 41.7 % | **49.3 %** | 2.0 % | 0.98× | 4.7× |
| A \ B | 63 | 4.0 % | 36.5 % | **56.0 %** | 2.1 % | 1.44× | 8.5× |
| B \ A | 31 | 3.7 % | 40.9 % | **51.0 %** | 2.0 % | 0.98× | 4.8× |
| B \ A | 63 | 3.9 % | 35.9 % | **56.8 %** | 2.1 % | 1.43× | 9.1× |

## Result 2 — chaining vs the sort you cannot avoid (`-t1`, median over ops × J)

The uncompacted mode must order the same k-mers; the only difference is *how*.

| dataset | k | chaining (s) | plain sort (s) | chain/sort | write nocompact/compact |
|---|--:|--:|--:|--:|--:|
| E. coli | 31 | 0.104 | 0.099 | 1.06 | 5.0× |
| E. coli | 63 | 0.117 | 0.153 | **0.77** | 8.7× |
| Yeast | 31 | 0.252 | 0.285 | 0.95 | 5.8× |
| Yeast | 63 | 0.288 | 0.422 | **0.69** | 9.6× |
| C. elegans | 31 | 2.362 | 2.470 | 0.98 | 5.2× |
| C. elegans | 63 | 2.395 | 3.864 | **0.64** | 10.5× |
| Chr1 | 31 | 5.462 | 5.631 | 1.02 | 5.6× |
| Chr1 | 63 | 5.854 | 9.269 | **0.63** | 11.8× |

The chaining works per column pair on already-ordered runs, so it beats one global sort of the whole
bucket as soon as the columns get numerous (k=63, m=31 ⇒ 33 columns): at k=63 chaining the k-mers is
**cheaper than sorting them**, and it produces 15× less output.

## Result 3 — end to end (s, median over ops × J)

| dataset | k | t | count-only | uncompacted | compacted | nc/c | compact/count |
|---|--:|--:|--:|--:|--:|--:|--:|
| E. coli | 31 | 1 | 0.048 | 0.211 | 0.172 | 1.17 | 3.03 |
| E. coli | 31 | 8 | 0.015 | 0.080 | 0.048 | 1.62 | 2.75 |
| E. coli | 63 | 1 | 0.044 | 0.307 | 0.189 | 1.54 | 3.29 |
| E. coli | 63 | 8 | 0.016 | 0.121 | 0.053 | 2.28 | 2.45 |
| Yeast | 31 | 1 | 0.112 | 0.566 | 0.405 | 1.26 | 3.26 |
| Yeast | 31 | 8 | 0.025 | 0.177 | 0.108 | 1.54 | 3.85 |
| Yeast | 63 | 1 | 0.097 | 0.890 | 0.444 | 1.76 | 3.54 |
| Yeast | 63 | 8 | 0.026 | 0.266 | 0.113 | 2.34 | 2.77 |
| C. elegans | 31 | 1 | 0.966 | 4.612 | 3.427 | 1.17 | 3.39 |
| C. elegans | 31 | 8 | 0.172 | 1.249 | 0.843 | 1.35 | 4.53 |
| C. elegans | 63 | 1 | 0.767 | 6.867 | 3.771 | 1.79 | 3.19 |
| C. elegans | 63 | 8 | 0.170 | 2.096 | 0.908 | 2.14 | 3.03 |
| Chr1 | 31 | 1 | 2.121 | 9.860 | 8.163 | 1.14 | 3.38 |
| Chr1 | 31 | 8 | 0.414 | 2.694 | 1.908 | 1.27 | 4.63 |
| Chr1 | 63 | 1 | 1.849 | 16.058 | 8.994 | 1.68 | 3.28 |
| Chr1 | 63 | 8 | 0.438 | 4.856 | 2.299 | 1.96 | 3.17 |

Per operation (median over the whole grid): uncompacted/compacted = 1.73 (∪), 1.53 (∩), 1.43 (A\B),
1.47 (B\A); compacted/count-only = 6.87 (∪), 4.45 (∩), 2.48 (A\B), 2.54 (B\A).

On chr1 alone, per operation (median over J):

| operation | k | t1 count | t1 uncomp. | t1 comp. | t8 count | t8 uncomp. | t8 comp. | bits/k-mer |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| A ∩ B | 31 | 1.890 | 9.244 | 7.781 | 0.292 | 2.271 | 1.744 | 19.2 |
| A ∩ B | 63 | 1.614 | 14.587 | 8.260 | 0.284 | 4.363 | 1.970 | 17.7 |
| A ∪ B | 31 | 2.109 | 14.546 | 12.943 | 0.385 | 4.075 | 3.212 | 18.7 |
| A ∪ B | 63 | 1.823 | 26.059 | 13.086 | 0.408 | 7.573 | 3.595 | 17.1 |
| A \ B | 31 | 2.429 | 6.439 | 5.726 | 0.482 | 1.619 | 1.310 | 20.6 |
| A \ B | 63 | 2.143 | 8.803 | 5.611 | 0.541 | 3.089 | 1.407 | 19.5 |
| B \ A | 31 | 2.427 | 6.882 | 6.136 | 0.489 | 1.749 | 1.383 | 20.3 |
| B \ A | 63 | 2.173 | 8.778 | 5.924 | 0.539 | 2.859 | 1.469 | 19.2 |

## Result 4 — output size (median over ops × J)

| dataset | k | records/k-mer | bits/k-mer | records/k-mer (nc) | bits/k-mer (nc) | ratio |
|---|--:|--:|--:|--:|--:|--:|
| E. coli | 31 | 0.122 | 17.69 | 1.000 | 144.14 | 8.15× |
| E. coli | 63 | 0.065 | 17.68 | 1.000 | 272.14 | 15.39× |
| Yeast | 31 | 0.122 | 17.65 | 1.000 | 144.05 | 8.16× |
| Yeast | 63 | 0.065 | 17.60 | 1.000 | 272.05 | 15.46× |
| C. elegans | 31 | 0.129 | 18.61 | 1.000 | 144.01 | 7.74× |
| C. elegans | 63 | 0.065 | 17.56 | 1.000 | 272.01 | 15.49× |
| Chr1 | 31 | 0.133 | 19.20 | 1.000 | 144.00 | 7.50× |
| Chr1 | 63 | 0.065 | 17.73 | 1.000 | 272.00 | 15.34× |

Uncompacted, one record per k-mer at the same 18-byte (k=31) / 34-byte (k=63) slot: 144 and 272
bits per k-mer exactly. The compacted result sits at the construction's own 17.6–19.2 bits/k-mer.

## Result 5 — closure

`A ∪ A` (the J=1 pair) re-chains to **exactly** the record count of A's own index, on all eight
dataset/k configurations (E. coli 508 430 / 268 739, yeast 1 298 130 / 687 476, C. elegans
11 527 083 / 5 733 614, chr1 25 806 166 / 13 415 299). A materialized result is therefore not a
degraded list: repeated operations do not accumulate compaction loss.

`benchmark/scripts/verify/rechain_verif.sh` additionally checks, per dataset/k/op, that the compacted
and uncompacted results are the same set (`|C △ N| = 0`, `|C| = |N|`), that `sskm query` answers
identically from both, and that both can feed a further set operation: **80 checks, 0 failures**
over the four datasets × k=31/63 × {∪, ∩, A\B} (`benchmark/results/rechain/verif.log`).

## Notes
- Phase timing is `-t1` only (the accumulators are plain doubles written from the inline loop), so
  Results 1 and 2 are single-thread; Results 3–5 cover t=1 and t=8.
- The write phase is measured against a warm page cache on NVMe; on slower storage the uncompacted
  mode's 5–12× bigger write would hurt more, never less.
- **Action taken:** `--no-compact`'s CLI help ("Much faster (avoids the dominant cost) but a larger
  output file") is inaccurate — the phase it avoids is replaced by a sort of the same k-mers, and the
  larger file costs more than the chaining saves. Reworded to say it trades ~8×/15× space for no
  reliable time gain, and is only there to keep the merge output raw.
