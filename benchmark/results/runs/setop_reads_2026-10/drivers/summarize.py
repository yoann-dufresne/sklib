#!/usr/bin/env python3
"""Markdown tables for RESULTS.md from the per-rep CSVs of this run (medians over reps).

    python3 drivers/summarize.py > /tmp/tables.md
"""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
PAIRS = [("hg002_A:hg002_B", "HG002"), ("gut_A:gut_B", "gut"), ("ocean_SRF:ocean_DCM", "ocean")]


def load(name):
    with open(DATA / name) as f:
        return list(csv.DictReader(f))


def med(xs):
    return st.median(xs) if xs else float("nan")


def main():
    runs = defaultdict(list)
    rss = defaultdict(int)
    for r in load("setop_runs.csv"):
        key = (r["pair"], r["k"], r["op"], r["threads"], r["tool"])
        runs[key].append(float(r["time_s"]))
        rss[key] = max(rss[key], int(r["peak_rss_kb"]))
    sizes = {(r["pair"], r["k"]): r for r in load("pair_sizes.csv")}
    cons = defaultdict(dict)
    for r in load("construct_runs.csv"):
        cons[(r["set"], r["k"])][r["tool"]] = r

    print("### Pairs\n")
    print("| pair | k | \\|A\\| | \\|B\\| | \\|A∩B\\| | \\|A∪B\\| | J |")
    print("|---|---|---:|---:|---:|---:|---:|")
    for pair, label in PAIRS:
        for k in ("31", "63"):
            s = sizes[(pair, k)]
            g = lambda x: f"{int(s[x]) / 1e9:.2f} G"
            print(f"| {label} | {k} | {g('sk_a')} | {g('sk_b')} | {g('inter')} | {g('union')} | {float(s['jaccard']):.2f} |")

    print("\n### Set operations (median s, 3 reps; ratio = KMC / sklib, > 1 means sklib faster)\n")
    print("| pair | k | op | t | sklib | KMC | ratio | sklib RSS | KMC RSS |")
    print("|---|---|---|---|---:|---:|---:|---:|---:|")
    for pair, label in PAIRS:
        for k in ("31", "63"):
            for op, th in (("union", "8"), ("joint", "8"), ("inter", "8"), ("diffab", "8"), ("diffba", "8"),
                           ("union", "1"), ("joint", "1")):
                s, c = runs.get((pair, k, op, th, "sklib")), runs.get((pair, k, op, th, "kmc"))
                if not s or not c:
                    continue
                ms, mc = med(s), med(c)
                print(f"| {label} | {k} | {op} | {th} | {ms:.1f} | {mc:.1f} | {mc / ms:.2f}× | "
                      f"{rss[(pair, k, op, th, 'sklib')] / 1024:.0f} MB | {rss[(pair, k, op, th, 'kmc')] / 1024:.0f} MB |")
            for th in ("8", "1"):
                s = runs.get((pair, k, "sizes", th, "sklib"))
                if s:
                    print(f"| {label} | {k} | sizes (count only) | {th} | {med(s):.1f} | — | — | "
                          f"{rss[(pair, k, 'sizes', th, 'sklib')] / 1024:.0f} MB | — |")

    print("\n### Construction (t=8, one run, cold input)\n")
    print("| set | k | sklib s | KMC s | sklib index | KMC index | sklib bits/k-mer | KMC bits/k-mer | sklib RSS | KMC RSS |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    nk = {}
    for (pair, k), s in sizes.items():
        nk[(s["set_a"], k)] = int(s["sk_a"])
        nk[(s["set_b"], k)] = int(s["sk_b"])
    for (st_, k), d in sorted(cons.items()):
        if "sklib" not in d or "kmc" not in d or (st_, k) not in nk:
            continue
        a, b, n = d["sklib"], d["kmc"], nk[(st_, k)]
        print(f"| {st_} | {k} | {float(a['time_s']):.0f} | {float(b['time_s']):.0f} | "
              f"{int(a['index_bytes']) / 1e9:.1f} GB | {int(b['index_bytes']) / 1e9:.1f} GB | "
              f"{8 * int(a['index_bytes']) / n:.1f} | {8 * int(b['index_bytes']) / n:.1f} | "
              f"{int(a['peak_rss_kb']) / 1024 ** 2:.1f} GB | {int(b['peak_rss_kb']) / 1024 ** 2:.1f} GB |")

    print("\n### Drift anchor (chr1 k=31, J=0.5 mutant, union; median s)\n")
    print("| when | sklib t=1 | KMC t=1 | sklib t=8 | KMC t=8 |")
    print("|---|---:|---:|---:|---:|")
    for name, when in (("anchor_start_runs.csv", "start (2026-10-05)"), ("anchor_end_runs.csv", "end (2026-10-06)")):
        a = defaultdict(list)
        for r in load(name):
            a[(r["tool"], r["threads"])].append(float(r["time_s"]))
        print(f"| {when} | {med(a[('sklib', '1')]):.2f} | {med(a[('kmc', '1')]):.2f} | "
              f"{med(a[('sklib', '8')]):.2f} | {med(a[('kmc', '8')]):.2f} |")


if __name__ == "__main__":
    main()
