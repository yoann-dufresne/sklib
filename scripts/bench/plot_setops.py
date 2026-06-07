#!/usr/bin/env python3
"""Generate the figures for the set-operations benchmark report from the sweep CSVs.

    python3 scripts/bench/plot_setops.py

Reads scripts/out/bench/setops_{scaling,overlap,ksweep,threads}.csv and writes PNGs to
scripts/out/bench/figs/. Missing CSVs are skipped. Uses the Agg backend (no display).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BENCH = os.path.join(os.path.dirname(__file__), "..", "out", "bench")
FIGS = os.path.join(BENCH, "figs")
os.makedirs(FIGS, exist_ok=True)

COL = {"sklib": "#1f77b4", "kmc": "#d62728", "cbl": "#2ca02c", "fmsi": "#9467bd"}
MARK = {"sklib": "o", "kmc": "s", "cbl": "^", "fmsi": "D"}


def load(name):
    p = os.path.join(BENCH, f"setops_{name}.csv")
    if not os.path.isfile(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def cell(rows, label, tool, action, field="sec"):
    for r in rows:
        if r["label"] == label and r["tool"] == tool and r["action"] == action:
            try:
                return float(r[field])
            except (ValueError, KeyError):
                return None
    return None


def labels_in_order(rows, by="kA"):
    seen = {}
    for r in rows:
        seen.setdefault(r["label"], float(r[by]))
    return sorted(seen, key=lambda l: seen[l])


# ---- Fig 1+2: time and RAM vs genome size (scaling, intersection) ----
def fig_scaling():
    rows = load("scaling")
    if not rows:
        return
    labs = labels_in_order(rows, "kA")
    x = [cell(rows, l, "sklib", "build_A", "result_kmers") / 1e6 for l in labs]  # |A| Mkmer
    for metric, field, ylab, fname, title in [
        ("time", "sec", "set-op wall time (s)", "fig_time_vs_size.png",
         "Intersection time vs dataset size (k=21, single core, B=A mutated 1%)"),
        ("ram", "rss_mb", "peak RSS (MB)", "fig_ram_vs_size.png",
         "Peak memory during intersection vs dataset size"),
    ]:
        plt.figure(figsize=(7, 5))
        for tool in ("sklib", "kmc", "cbl"):
            y = [cell(rows, l, tool, "inter", field) for l in labs]
            plt.plot(x, y, marker=MARK[tool], color=COL[tool], label=tool, lw=2)
        if metric == "time":
            y = [cell(rows, l, "sklib", "inter_size", field) for l in labs]
            plt.plot(x, y, marker="x", color=COL["sklib"], ls="--", label="sklib _size (no materialize)", lw=2)
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("|A| distinct k-mers (millions)"); plt.ylabel(ylab)
        plt.title(title); plt.grid(True, which="both", ls=":", alpha=0.5); plt.legend()
        plt.tight_layout(); plt.savefig(os.path.join(FIGS, fname), dpi=120); plt.close()


# ---- Fig 3: intersection time vs overlap (Jaccard) ----
def fig_overlap():
    rows = load("overlap")
    if not rows:
        return
    labs = labels_in_order(rows, "kA")  # all same base size; reorder by jaccard below
    pts = []
    for l in labs:
        inter = cell(rows, l, "sklib", "inter", "result_kmers")
        uni = cell(rows, l, "sklib", "union", "result_kmers")
        if inter and uni:
            pts.append((100.0 * inter / uni, l))
    pts.sort()
    jac = [p[0] for p in pts]; labs = [p[1] for p in pts]
    plt.figure(figsize=(7, 5))
    for tool in ("sklib", "kmc", "cbl"):
        y = [cell(rows, l, tool, "inter", "sec") for l in labs]
        plt.plot(jac, y, marker=MARK[tool], color=COL[tool], label=tool, lw=2)
    y = [cell(rows, l, "sklib", "inter_size", "sec") for l in labs]
    plt.plot(jac, y, marker="x", color=COL["sklib"], ls="--", label="sklib _size", lw=2)
    plt.xlabel("Jaccard similarity (%)  =  |A∩B| / |A∪B|")
    plt.ylabel("intersection wall time (s)")
    plt.title("Intersection time vs overlap (chr21 base, 32.7M k-mers, k=21)")
    plt.grid(True, ls=":", alpha=0.5); plt.legend()
    plt.tight_layout(); plt.savefig(os.path.join(FIGS, "fig_time_vs_overlap.png"), dpi=120); plt.close()


# ---- Fig 4: threads scaling — sklib (bucket-parallel) vs KMC vs FMSI, one panel per pair ----
def fig_threads():
    rows = load("parallel")
    if not rows:
        return
    pairs = []
    for r in rows:
        if r["pair"] not in pairs:
            pairs.append(r["pair"])

    def series(pair, tool, op="inter"):
        pts = []
        for r in rows:
            if r["pair"] == pair and r["tool"] == tool and r["op"] == op and r["stage"] == "setop":
                try:
                    pts.append((int(r["threads"]), float(r["sec"])))
                except ValueError:
                    pass  # NA (timed out) or non-numeric
        pts.sort()
        return [p[0] for p in pts], [p[1] for p in pts]

    n = len(pairs); ncol = min(3, n); nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow), squeeze=False)
    for idx, pair in enumerate(pairs):
        ax = axes[idx // ncol][idx % ncol]
        for tool in ("sklib", "kmc"):
            xs, ys = series(pair, tool)
            if xs:
                ax.plot(xs, ys, marker=MARK[tool], color=COL[tool], label=tool, lw=2)
        xs, ys = series(pair, "sklib", "inter_size")
        if xs:
            ax.plot(xs, ys, marker="x", color=COL["sklib"], ls="--", label="sklib _size", lw=1.5)
        xs, ys = series(pair, "fmsi")            # serial reference (single point at t=1)
        if xs:
            ax.scatter([1], [ys[0]], color=COL["fmsi"], marker=MARK["fmsi"], s=70, zorder=5, label="fmsi (serial)")
        ax.set_title(pair); ax.set_xlabel("threads"); ax.set_ylabel("intersection time (s)")
        ax.set_xscale("log", base=2); ax.set_yscale("log")
        ax.set_xticks([1, 4, 8, 22]); ax.set_xticklabels(["1", "4", "8", "22"])
        ax.grid(True, which="both", ls=":", alpha=0.5); ax.legend(fontsize=8)
    for idx in range(n, nrow * ncol):
        axes[idx // ncol][idx % ncol].axis("off")
    fig.suptitle("Set-op intersection: thread scaling — sklib (bucket-parallel, v0.5.1) vs KMC vs FMSI (k=21)")
    plt.tight_layout(); plt.savefig(os.path.join(FIGS, "fig_threads.png"), dpi=120); plt.close()


def main():
    fig_scaling(); fig_overlap(); fig_threads()
    print("figures written to", FIGS)
    for f in sorted(os.listdir(FIGS)):
        print(" ", f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
