#!/usr/bin/env python3
"""Cross-tool comparison tables from scripts/out/bench/results.csv.

Dedups each (tool,dataset,k,m,phase,workload,threads) to its most recent timestamp (same rule
as plots.py), so sklib shows the latest build and competitors their single measurement. Prints
a Markdown table per (dataset, k, m): space (index size, bits/k-mer), construction (time, peak
RSS), and query (throughput + peak RSS) for the positive and streaming workloads.

Usage: compare.py --km 21,11 --datasets "ecoli yeast chr21" [--threads 1]
"""
import csv, sys, argparse, os

CSV = os.path.join(os.path.dirname(__file__), "..", "out", "bench", "results.csv")

def load(path):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (r["tool"], r["dataset"], r["k"], r["m"], r["phase"], r["workload"], r["threads"])
            cur = rows.get(key)
            if cur is None or r["timestamp"] > cur["timestamp"]:
                rows[key] = r
    return rows

def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def mb(x):
    v = num(x)
    return f"{v/1e6:.1f}" if v is not None else "–"

def rss_mb(x):
    v = num(x)
    return f"{v/1024:.0f}" if v is not None else "–"

def f2(x, nd=2):
    v = num(x)
    return f"{v:.{nd}f}" if v is not None else "–"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--km", default="21,11")
    ap.add_argument("--datasets", default="ecoli yeast celegans chr21 chr1")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--tools", default="sklib sshash sbwt cbl bqf")
    ap.add_argument("--since", default="", help="only show a dataset if its sklib construct row timestamp >= this")
    args = ap.parse_args()
    k, m = args.km.split(",")
    th = args.threads
    tools = args.tools.split()

    rows = load(CSV)
    for ds in args.datasets.split():
        # only show datasets whose latest sklib construct row is from this run (>= --since)
        sk = rows.get(("sklib", ds, k, m, "construct", "-", th))
        if sk is None or sk["timestamp"] < args.since:
            continue
        print(f"\n### {ds}  (k={k}, m={m}, threads={th})\n")
        hdr = ("| tool | index MB | bits/kmer | constr s | constr RSS MB "
               "| pos Mkmer/s | pos RSS MB | stream Mkmer/s |")
        print(hdr)
        print("|" + "---|" * 8)
        for t in tools:
            c = rows.get((t, ds, k, m, "construct", "-", th))
            qp = rows.get((t, ds, k, m, "query", "positive", th))
            qs = rows.get((t, ds, k, m, "query", "streaming", th))
            if c is None and qp is None:
                continue
            label = "**sklib (new)**" if t == "sklib" else t
            idx = mb(c["index_bytes"]) if c else "–"
            bpk = f2(c["bits_per_kmer"]) if c else "–"
            cs  = f2(c["time_s"]) if c else "–"
            crss= rss_mb(c["peak_rss_kb"]) if c else "–"
            pth = f2(qp["throughput_Mkmer_s"], 3) if qp else "–"
            prss= rss_mb(qp["peak_rss_kb"]) if qp else "–"
            sth = f2(qs["throughput_Mkmer_s"], 3) if qs else "–"
            print(f"| {label} | {idx} | {bpk} | {cs} | {crss} | {pth} | {prss} | {sth} |")

if __name__ == "__main__":
    main()
