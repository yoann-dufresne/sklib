#!/usr/bin/env python3
"""Aggregate benchmark/results/rechain/setop_rechain.csv into the tables of the re-chaining study.

Three views, all medians over the four operations and the seven Jaccard targets:
  1. phase split at -t1 (read / merge+collect / order / write) for both materializing modes;
  2. end-to-end wall time of the three modes (cardinality-only, uncompacted, compacted);
  3. output size (records per k-mer, bits per k-mer) of the two materializing modes.

    python3 benchmark/scripts/rechain_report.py [csv] [--tex]
"""
import csv
import sys
from collections import defaultdict
from statistics import median

DATASET_ORDER = ["ecoli", "yeast", "celegans", "chr1"]
DATASET_TEX = {"ecoli": r"\emph{E.~coli}", "yeast": "Yeast",
               "celegans": r"\emph{C.~elegans}", "chr1": "Chr1"}
OP_ORDER = ["inter", "union", "diffab", "diffba"]
MODES = ["size", "nocompact", "compact"]


def load(path, keep_empty=False):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                r["k"] = int(r["k"]); r["threads"] = int(r["threads"])
                r["time_s"] = float(r["time_s"])
            except ValueError:
                continue
            for f in ("out_bytes", "out_records", "result_kmers"):
                r[f] = int(r[f]) if r[f] not in ("", "NA") else None
            for f in ("records_per_kmer", "bits_per_kmer",
                      "ph_read", "ph_merge", "ph_order", "ph_write", "ph_total"):
                r[f] = float(r[f]) if r[f] not in ("", "NA") else None
            if not keep_empty and r["result_kmers"] == 0:
                continue      # J=1 differences are empty: nothing is ordered or written
            rows.append(r)
    return rows


def med(vals):
    vals = [v for v in vals if v is not None]
    return median(vals) if vals else None


def fmt(v, nd=2):
    return "NA" if v is None else f"{v:.{nd}f}"


def group(rows, keyf):
    g = defaultdict(list)
    for r in rows:
        g[keyf(r)].append(r)
    return g


def datasets_present(rows):
    seen = {r["dataset"] for r in rows}
    return [d for d in DATASET_ORDER if d in seen]


def ks_present(rows):
    return sorted({r["k"] for r in rows})


def threads_present(rows):
    return sorted({r["threads"] for r in rows})


# ---- view 1: phase split at -t1 -------------------------------------------
def phase_table(rows):
    out = []
    sel = [r for r in rows if r["threads"] == 1 and r["ph_total"] is not None]
    g = group(sel, lambda r: (r["dataset"], r["k"], r["mode"]))
    for ds in datasets_present(sel):
        for k in ks_present(sel):
            for mode in ("compact", "nocompact"):
                rs = g.get((ds, k, mode), [])
                if not rs:
                    continue
                shares = {}
                for ph in ("ph_read", "ph_merge", "ph_order", "ph_write"):
                    shares[ph] = med([100.0 * r[ph] / r["ph_total"] for r in rs if r["ph_total"]])
                out.append({"dataset": ds, "k": k, "mode": mode, "n": len(rs),
                            "total": med([r["ph_total"] for r in rs]),
                            "order_s": med([r["ph_order"] for r in rs]),
                            **shares})
    return out


# ---- view 1b: chaining vs the plain sort of the same k-mers ----------------
def order_table(rows):
    out = []
    sel = [r for r in rows if r["threads"] == 1 and r["ph_order"] is not None]
    g = group(sel, lambda r: (r["dataset"], r["k"], r["mode"]))
    for ds in datasets_present(sel):
        for k in ks_present(sel):
            pairs = defaultdict(dict)
            for mode in ("compact", "nocompact"):
                for r in g.get((ds, k, mode), []):
                    pairs[(r["op"], r["jaccard_target"])][mode] = r
            ratios = [v["compact"]["ph_order"] / v["nocompact"]["ph_order"]
                      for v in pairs.values()
                      if "compact" in v and "nocompact" in v and v["nocompact"]["ph_order"]]
            wratios = [v["nocompact"]["ph_write"] / v["compact"]["ph_write"]
                       for v in pairs.values()
                       if "compact" in v and "nocompact" in v and v["compact"]["ph_write"]]
            if ratios:
                out.append({"dataset": ds, "k": k, "n": len(ratios),
                            "chain_over_sort": med(ratios),
                            "chain_s": med([v["compact"]["ph_order"] for v in pairs.values() if "compact" in v]),
                            "sort_s": med([v["nocompact"]["ph_order"] for v in pairs.values() if "nocompact" in v]),
                            "write_ratio": med(wratios)})
    return out


# ---- view 2: end-to-end wall time -----------------------------------------
def time_table(rows):
    out = []
    g = group(rows, lambda r: (r["dataset"], r["k"], r["threads"], r["mode"]))
    for ds in datasets_present(rows):
        for k in ks_present(rows):
            for th in threads_present(rows):
                e = {"dataset": ds, "k": k, "threads": th}
                for mode in MODES:
                    rs = g.get((ds, k, th, mode), [])
                    e[mode] = med([r["time_s"] for r in rs])
                    e[mode + "_n"] = len(rs)
                # Paired per-(op,J) ratios, so the ratio is not a ratio of medians.
                pairs = defaultdict(dict)
                for mode in MODES:
                    for r in g.get((ds, k, th, mode), []):
                        pairs[(r["op"], r["jaccard_target"])][mode] = r["time_s"]
                e["nocompact_over_compact"] = med([v["nocompact"] / v["compact"]
                                                   for v in pairs.values()
                                                   if "nocompact" in v and "compact" in v and v["compact"]])
                e["compact_over_size"] = med([v["compact"] / v["size"]
                                              for v in pairs.values()
                                              if "size" in v and "compact" in v and v["size"]])
                if e["compact"] is not None:
                    out.append(e)
    return out


# ---- view 3: output size ---------------------------------------------------
def size_table(rows):
    out = []
    sel = [r for r in rows if r["out_bytes"] is not None and r["threads"] == 1]
    g = group(sel, lambda r: (r["dataset"], r["k"], r["mode"]))
    for ds in datasets_present(sel):
        for k in ks_present(sel):
            e = {"dataset": ds, "k": k}
            for mode in ("compact", "nocompact"):
                rs = g.get((ds, k, mode), [])
                e[mode + "_rpk"] = med([r["records_per_kmer"] for r in rs])
                e[mode + "_bpk"] = med([r["bits_per_kmer"] for r in rs])
            pairs = defaultdict(dict)
            for mode in ("compact", "nocompact"):
                for r in g.get((ds, k, mode), []):
                    pairs[(r["op"], r["jaccard_target"])][mode] = r["bits_per_kmer"]
            e["ratio"] = med([v["nocompact"] / v["compact"] for v in pairs.values()
                              if "compact" in v and "nocompact" in v and v["compact"]])
            if e["compact_bpk"] is not None:
                out.append(e)
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0] if args else "benchmark/results/rechain/setop_rechain.csv"
    rows = load(path)
    print(f"# {len(rows)} rows from {path}\n")

    print("== 1. phase split, -t1 (median over ops x Jaccard targets, % of in-loop wall) ==")
    print(f"{'dataset':10} {'k':>3} {'mode':10} {'n':>4} {'total_s':>8} "
          f"{'read%':>7} {'merge%':>7} {'order%':>7} {'write%':>7}")
    for e in phase_table(rows):
        print(f"{e['dataset']:10} {e['k']:>3} {e['mode']:10} {e['n']:>4} {fmt(e['total'], 3):>8} "
              f"{fmt(e['ph_read'], 1):>7} {fmt(e['ph_merge'], 1):>7} "
              f"{fmt(e['ph_order'], 1):>7} {fmt(e['ph_write'], 1):>7}")

    print("\n== 1b. ordering the kept k-mers: chaining vs plain sort, -t1 (median over ops x J) ==")
    print(f"{'dataset':10} {'k':>3} {'n':>4} {'chain_s':>9} {'sort_s':>9} {'chain/sort':>11} {'write nc/c':>11}")
    for e in order_table(rows):
        print(f"{e['dataset']:10} {e['k']:>3} {e['n']:>4} {fmt(e['chain_s'], 3):>9} {fmt(e['sort_s'], 3):>9} "
              f"{fmt(e['chain_over_sort']):>11} {fmt(e['write_ratio'], 1):>11}")

    print("\n== 2. end-to-end wall time, s (median over ops x Jaccard targets) ==")
    print(f"{'dataset':10} {'k':>3} {'t':>3} {'size':>8} {'nocomp':>8} {'compact':>8} "
          f"{'nc/c':>6} {'c/size':>7} {'n':>4}")
    for e in time_table(rows):
        print(f"{e['dataset']:10} {e['k']:>3} {e['threads']:>3} {fmt(e['size'], 3):>8} "
              f"{fmt(e['nocompact'], 3):>8} {fmt(e['compact'], 3):>8} "
              f"{fmt(e['nocompact_over_compact']):>6} {fmt(e['compact_over_size']):>7} {e['compact_n']:>4}")

    print("\n== 3. output size (median over ops x Jaccard targets) ==")
    print(f"{'dataset':10} {'k':>3} {'c_rec/kmer':>11} {'c_bits/kmer':>12} "
          f"{'nc_rec/kmer':>12} {'nc_bits/kmer':>13} {'ratio':>6}")
    for e in size_table(rows):
        print(f"{e['dataset']:10} {e['k']:>3} {fmt(e['compact_rpk'], 3):>11} {fmt(e['compact_bpk']):>12} "
              f"{fmt(e['nocompact_rpk'], 3):>12} {fmt(e['nocompact_bpk']):>13} {fmt(e['ratio']):>6}")

    # Spread of the interesting ratios across the grid, for the prose.
    print("\n== 4. spread over the whole grid ==")
    g = group(rows, lambda r: (r["dataset"], r["k"], r["threads"], r["op"], r["jaccard_target"]))
    nc, cs, ordsh = [], [], []
    for key, rs in g.items():
        by = {r["mode"]: r for r in rs}
        if "nocompact" in by and "compact" in by and by["compact"]["time_s"]:
            nc.append(by["nocompact"]["time_s"] / by["compact"]["time_s"])
        if "size" in by and "compact" in by and by["size"]["time_s"]:
            cs.append(by["compact"]["time_s"] / by["size"]["time_s"])
        if "compact" in by and by["compact"]["ph_total"]:
            ordsh.append(100.0 * by["compact"]["ph_order"] / by["compact"]["ph_total"])
    for name, vals in (("nocompact/compact time", nc), ("compact/size time", cs),
                       ("order share %, t=1 compact", ordsh)):
        if vals:
            vals.sort()
            print(f"{name:28} min={vals[0]:.2f}  med={median(vals):.2f}  max={vals[-1]:.2f}  n={len(vals)}")

    print("\n== 5. closure: records of A u A vs records of A's own index ==")
    for ds, k, got, want in closure_check(rows):
        flag = "OK" if got == want else "MISMATCH"
        print(f"{ds:10} k={k:<3} union(A,A)={got}  index(A)={want}  {flag}")


def closure_check(rows, idx_root=None):
    """union at J=1 rebuilds A: its record count must equal A's own index."""
    import glob
    import os
    out = []
    sel = [r for r in rows if r["op"] == "union" and r["jaccard_target"] in ("1.0", "1")
           and r["mode"] == "compact" and r["out_records"]]
    idx_root = idx_root or "benchmark/results/rechain/indexes/sklib"
    for r in sel:
        pat = os.path.join(idx_root, "*", f"{r['dataset']}.k{r['k']}.m{r['m']}", ".idxmeta")
        meta = glob.glob(pat)
        if not meta:
            continue
        n = None
        for line in open(meta[0]):
            if line.startswith("N_SKMERS="):
                n = int(line.split("=", 1)[1].strip().strip("'\""))
        out.append((r["dataset"], r["k"], r["out_records"], n))
    return sorted(set(out))


if __name__ == "__main__":
    main()
