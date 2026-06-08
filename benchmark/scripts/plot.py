#!/usr/bin/env python3
"""Figures for the 4-experiment benchmark harness.

Reads the four CSVs written by construct.sh / query_single.sh / query_stream.sh /
setop.sh and draws one figure family per experiment under <results>/figs/:

  fig_construct.png       time / peak RSS / bits-per-kmer vs threads, per (dataset,k)
  fig_query_single.png    throughput vs presence-rate (line per thread), per (dataset,k)
  fig_query_stream.png    same, for stream queries
  fig_setop_jaccard.png   set-op time vs Jaccard (line per op), per (dataset,k)
  fig_setop_joint.png     joint single-pass vs sequential (sum of unitary), per tool

Missing CSVs are skipped. Single-thread-only tools (one threads value) draw flat.
    python3 benchmark/scripts/plot.py [results_dir]      # default: results/latest
Dedups by the latest timestamp per identity, so re-runs/accumulated CSVs stay clean.
"""
import os
import sys
import math

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "results", "latest")
FIGS = os.path.join(RES, "figs")
os.makedirs(FIGS, exist_ok=True)


def load(name, keys):
    """Load <RES>/<name>.csv, keep the latest row per identity tuple `keys`."""
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        print("skip (missing):", path)
        return None
    df = pd.read_csv(path)
    if df.empty:
        return None
    if "timestamp" in df:
        df = df.sort_values("timestamp").drop_duplicates(keys, keep="last")
    return df


def facets(df):
    """(dataset,k) panels, sorted; returns the list and a grid size."""
    combos = sorted(df[["dataset", "k"]].drop_duplicates().itertuples(index=False, name=None),
                    key=lambda x: (x[0], x[1]))
    n = len(combos)
    cols = min(3, n) or 1
    rows = math.ceil(n / cols)
    return combos, rows, cols


def grid(df, title):
    combos, rows, cols = facets(df)
    fig, ax = plt.subplots(rows, cols, figsize=(5 * cols, 3.4 * rows), squeeze=False)
    fig.suptitle(title, fontsize=13)
    axes = [ax[i // cols][i % cols] for i in range(rows * cols)]
    for a in axes[len(combos):]:
        a.axis("off")
    return fig, combos, axes


def plot_construct():
    df = load("construct.csv", ["tool", "tool_version", "host", "dataset", "k", "m", "threads"])
    if df is None:
        return
    for metric, ylab, logy in [("time_s", "build time (s)", True),
                               ("peak_rss_kb", "peak RSS (MB)", True),
                               ("bits_per_kmer", "bits / k-mer", False)]:
        fig, combos, axes = grid(df, "Construction — %s" % ylab)
        for a, (ds, k) in zip(axes, combos):
            sub = df[(df.dataset == ds) & (df.k == k)]
            for tool in sorted(sub.tool.unique()):
                t = sub[sub.tool == tool].sort_values("threads")
                y = t[metric] / 1024.0 if metric == "peak_rss_kb" else t[metric]
                a.plot(t.threads, y, marker="o", label=tool)
            a.set_title("%s k=%d" % (ds, k)); a.set_xlabel("threads"); a.set_ylabel(ylab)
            if logy:
                a.set_yscale("log")
            a.grid(True, alpha=.3); a.legend(fontsize=7)
        fig.tight_layout(rect=[0, 0, 1, .96])
        out = os.path.join(FIGS, "fig_construct_%s.png" % metric)
        fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def plot_query(name):
    df = load(name + ".csv", ["tool", "tool_version", "host", "dataset", "k", "m", "threads", "presence"])
    if df is None:
        return
    fig, combos, axes = grid(df, "%s — throughput vs present-fraction" % name)
    for a, (ds, k) in zip(axes, combos):
        sub = df[(df.dataset == ds) & (df.k == k)]
        for tool in sorted(sub.tool.unique()):
            for th in sorted(sub.threads.unique()):
                t = sub[(sub.tool == tool) & (sub.threads == th)].sort_values("presence")
                if len(t):
                    a.plot(t.presence, t.throughput_Mkmer_s, marker="o",
                           label="%s t=%d" % (tool, th))
        a.set_title("%s k=%d" % (ds, k)); a.set_xlabel("present k-mers (%)")
        a.set_ylabel("Mk-mer/s"); a.grid(True, alpha=.3); a.legend(fontsize=6)
    fig.tight_layout(rect=[0, 0, 1, .96])
    out = os.path.join(FIGS, "fig_%s.png" % name)
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def plot_setop():
    df = load("setop.csv", ["tool", "tool_version", "host", "dataset", "k", "m",
                            "threads", "op", "mode", "jaccard_target"])
    if df is None:
        return
    mat = df[df["mode"] == "materialize"].copy()
    # (1) time vs measured Jaccard, line per op, sklib at max threads
    fig, combos, axes = grid(mat, "Set-ops (materialize) — time vs Jaccard, sklib @ max threads")
    for a, (ds, k) in zip(axes, combos):
        sub = mat[(mat.dataset == ds) & (mat.k == k) & (mat.tool == "sklib")]
        if len(sub):
            th = sub.threads.max(); sub = sub[sub.threads == th]
            for op in ["inter", "union", "diffab", "diffba", "joint"]:
                t = sub[sub.op == op].sort_values("jaccard_measured")
                if len(t):
                    a.plot(t.jaccard_measured, t.time_s, marker="o", label=op)
        a.set_title("%s k=%d" % (ds, k)); a.set_xlabel("Jaccard"); a.set_ylabel("time (s)")
        a.grid(True, alpha=.3); a.legend(fontsize=7)
    fig.tight_layout(rect=[0, 0, 1, .96])
    out = os.path.join(FIGS, "fig_setop_jaccard.png"); fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)

    # (2) joint single-pass vs sequential (sum of the 4 unitary ops), per tool, at one mid Jaccard
    js = mat[mat.op == "joint"]
    if len(js):
        ds, k = js.iloc[0][["dataset", "k"]]
        sub = mat[(mat.dataset == ds) & (mat.k == k)]
        jt = sorted(sub.jaccard_target.unique())[len(sub.jaccard_target.unique()) // 2]
        sub = sub[sub.jaccard_target == jt]
        tools = sorted(sub.tool.unique()); xs = range(len(tools)); w = .35
        joint = [sub[(sub.tool == t) & (sub.op == "joint")].time_s.min() for t in tools]
        seq = [sub[(sub.tool == t) & (sub.op.isin(["inter", "union", "diffab", "diffba"]))]
               .groupby("threads").time_s.sum().min() for t in tools]
        fig, a = plt.subplots(figsize=(1.6 * len(tools) + 3, 4))
        a.bar([x - w / 2 for x in xs], seq, w, label="sequential (4 passes)")
        a.bar([x + w / 2 for x in xs], joint, w, label="joint (1 pass)")
        a.set_xticks(list(xs)); a.set_xticklabels(tools)
        a.set_ylabel("time (s)"); a.set_title("%s k=%d J~%s — joint vs sequential" % (ds, k, jt))
        a.grid(True, axis="y", alpha=.3); a.legend()
        fig.tight_layout(); out = os.path.join(FIGS, "fig_setop_joint.png")
        fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    plot_construct()
    plot_query("query_single")
    plot_query("query_stream")
    plot_setop()
    print("figures ->", FIGS)
