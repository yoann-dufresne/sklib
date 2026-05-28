#!/usr/bin/env python3
"""Turn scripts/bench/results.csv into the publication figures (Fig 1-7).

Run with the bench venv:
    scripts/bench/.venv/bin/python scripts/bench/plots.py [results.csv] [out_dir]

Each figure is defensive: it is skipped (with a note) when the CSV lacks the data
it needs, so the script works whether you have sklib alone or sklib + competitors,
one (k,m) point or a full sweep, one thread count or a scaling series.
"""
import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

NUMERIC = ["dataset_bp", "distinct_kmers", "k", "m", "threads", "n_queries",
           "index_bytes", "bytes_per_skmer", "n_superkmers", "kmers_per_superkmer",
           "bits_per_kmer", "time_s", "peak_rss_kb", "throughput_Mkmer_s"]


def load(path):
    df = pd.read_csv(path)
    for c in NUMERIC:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # The CSV is append-only and may hold several sklib commits (e.g. the
    # construction-RAM rewrite V0/V1/V2). The figures show the *current* state,
    # so collapse each measurement to its most recent timestamp. Competitors,
    # measured once, are untouched; superseded sklib rows drop out.
    if "timestamp" in df.columns:
        key = ["tool", "dataset", "k", "m", "phase", "workload", "threads"]
        df = (df.assign(_ts=pd.to_datetime(df.timestamp, errors="coerce", utc=True))
                .sort_values("_ts")
                .drop_duplicates(subset=key, keep="last")
                .drop(columns="_ts"))
    return df


def save(fig, out_dir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def _markers():
    return ["o", "s", "^", "D", "v", "P", "X", "*"]


def fig1_space(df, out_dir):
    """bits/k-mer vs m, one line per k, faceted by dataset; competitors as ref points."""
    c = df[df.phase == "construct"].dropna(subset=["bits_per_kmer", "m"])
    if c.empty:
        print("  [fig1] no construct rows -> skip"); return
    datasets = sorted(c.dataset.unique())
    fig, axes = plt.subplots(1, len(datasets), figsize=(5 * len(datasets), 4.2), squeeze=False)
    for ax, ds in zip(axes[0], datasets):
        sub = c[c.dataset == ds]
        sk = sub[sub.tool == "sklib"]
        for k in sorted(sk.k.dropna().unique()):
            kk = sk[sk.k == k].sort_values("m")
            ax.plot(kk.m, kk.bits_per_kmer, marker="o", label=f"sklib k={int(k)}")
            ax.axhline(2 * k, ls=":", lw=0.8, color="grey")
            ax.annotate(f"2-bit raw k={int(k)}", (kk.m.min(), 2 * k), fontsize=7, color="grey", va="bottom")
        # competitors: a marker per (tool,k) at their (constant) bits/kmer
        for tool in sorted(sub[sub.tool != "sklib"].tool.unique()):
            t = sub[sub.tool == tool]
            ax.scatter(t.m, t.bits_per_kmer, marker="*", s=90, label=tool, zorder=5)
        ax.set_title(ds); ax.set_xlabel("minimizer length m"); ax.set_ylabel("bits / k-mer")
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle("Fig 1 — Space: bits per k-mer vs m")
    save(fig, out_dir, "fig1_bits_per_kmer")


def fig2_compaction(df, out_dir):
    """k-mers per super-k-mer vs m (per k, per dataset) + n_superkmers vs n_kmers."""
    c = df[(df.phase == "construct") & (df.tool == "sklib")].dropna(subset=["kmers_per_superkmer", "m"])
    if c.empty:
        print("  [fig2] no sklib construct rows -> skip"); return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for ds in sorted(c.dataset.unique()):
        sub = c[c.dataset == ds]
        for k in sorted(sub.k.dropna().unique()):
            kk = sub[sub.k == k].sort_values("m")
            ax1.plot(kk.m, kk.kmers_per_superkmer, marker="o", label=f"{ds} k={int(k)}")
    ax1.set_xlabel("minimizer length m"); ax1.set_ylabel("mean k-mers / super-k-mer")
    ax1.set_title("compaction vs m"); ax1.grid(alpha=0.3); ax1.legend(fontsize=7)

    sc = c.dropna(subset=["n_superkmers", "distinct_kmers"])
    for ds in sorted(sc.dataset.unique()):
        sub = sc[sc.dataset == ds]
        ax2.scatter(sub.distinct_kmers, sub.n_superkmers, label=ds)
    if not sc.empty:
        lim = [sc.distinct_kmers.min(), sc.distinct_kmers.max()]
        ax2.plot(lim, lim, ls=":", color="grey", label="y=x (1 kmer/skmer)")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("# distinct k-mers"); ax2.set_ylabel("# super-k-mers")
    ax2.set_title("super-k-mers vs k-mers"); ax2.grid(alpha=0.3); ax2.legend(fontsize=7)
    fig.suptitle("Fig 2 — Compaction")
    save(fig, out_dir, "fig2_compaction")


def fig3_construction(df, out_dir):
    """Construction throughput and peak RSS vs dataset size (scaling)."""
    c = df[df.phase == "construct"].dropna(subset=["dataset_bp"])
    if c.empty:
        print("  [fig3] no construct rows -> skip"); return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    marks = _markers()
    for i, tool in enumerate(sorted(c.tool.unique())):
        t = c[c.tool == tool].sort_values("dataset_bp")
        tp = t.dropna(subset=["throughput_Mkmer_s"])
        ax1.plot(tp.dataset_bp, tp.throughput_Mkmer_s, marker=marks[i % len(marks)], ls="-", label=tool)
        rs = t.dropna(subset=["peak_rss_kb"])
        ax2.plot(rs.dataset_bp, rs.peak_rss_kb / 1024.0, marker=marks[i % len(marks)], ls="-", label=tool)
    for ax in (ax1, ax2):
        ax.set_xscale("log"); ax.set_xlabel("dataset size (bp)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax1.set_ylabel("construct throughput (Mkmer/s)"); ax1.set_title("construction speed vs scale")
    ax2.set_yscale("log"); ax2.set_ylabel("peak RSS (MB)"); ax2.set_title("construction memory vs scale")
    fig.suptitle("Fig 3 — Construction scaling")
    save(fig, out_dir, "fig3_construction")


def fig4_query(df, out_dir):
    """Query throughput per workload (bars) + thread scaling (if available)."""
    q = df[df.phase == "query"].dropna(subset=["throughput_Mkmer_s"])
    if q.empty:
        print("  [fig4] no query rows -> skip"); return
    # representative dataset = most distinct k-mers
    ds = q.loc[q.distinct_kmers.idxmax(), "dataset"] if q.distinct_kmers.notna().any() else q.dataset.iloc[0]
    qd = q[q.dataset == ds]
    tmax = int(qd.threads.max())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    bars = qd[qd.threads == tmax].groupby("workload").throughput_Mkmer_s.max()
    ax1.bar(bars.index, bars.values, color="steelblue")
    ax1.set_ylabel("throughput (Mkmer/s)")
    ax1.set_title(f"{ds}: query throughput by workload (threads={tmax})")
    ax1.tick_params(axis="x", rotation=20); ax1.grid(alpha=0.3, axis="y")

    threads = sorted(qd.threads.dropna().unique())
    if len(threads) > 1:
        for wl in sorted(qd.workload.unique()):
            w = qd[qd.workload == wl].groupby("threads").throughput_Mkmer_s.max().sort_index()
            ax2.plot(w.index, w.values, marker="o", label=wl)
        ax2.set_xlabel("threads (CPU affinity)"); ax2.set_ylabel("throughput (Mkmer/s)")
        ax2.set_title(f"{ds}: thread scaling"); ax2.grid(alpha=0.3); ax2.legend(fontsize=8)
    else:
        ax2.text(0.5, 0.5, "single thread count\n(rebuild with TBB +\nset THREADS for scaling)",
                 ha="center", va="center", transform=ax2.transAxes, color="grey")
        ax2.set_axis_off()
    fig.suptitle("Fig 4 — Query performance")
    save(fig, out_dir, "fig4_query")


def fig5_pareto(df, out_dir, workload="streaming"):
    """Pareto: bits/k-mer (x) vs query throughput (y). sklib swept over m; competitors as points."""
    c = df[df.phase == "construct"]
    q = df[(df.phase == "query") & (df.workload == workload)]
    if c.empty or q.empty:
        print(f"  [fig5] need construct + '{workload}' query rows -> skip"); return
    tmax = int(q.threads.max())
    q = q[q.threads == tmax]
    # join space (from construct) with throughput (from query) on the index key
    key = ["tool", "dataset", "k", "m"]
    merged = pd.merge(c[key + ["bits_per_kmer"]].dropna(),
                      q[key + ["throughput_Mkmer_s"]].dropna(), on=key)
    if merged.empty:
        print("  [fig5] no matching construct/query pairs -> skip"); return
    fig, ax = plt.subplots(figsize=(7, 5))
    marks = _markers()
    datasets = sorted(merged.dataset.unique())
    cmap = plt.get_cmap("tab10")
    for di, ds in enumerate(datasets):
        sub = merged[merged.dataset == ds]
        sk = sub[sub.tool == "sklib"].sort_values("m")
        if not sk.empty:
            ax.plot(sk.bits_per_kmer, sk.throughput_Mkmer_s, "-o", color=cmap(di), label=f"sklib {ds}")
            for _, r in sk.iterrows():
                ax.annotate(f"m={int(r.m)}", (r.bits_per_kmer, r.throughput_Mkmer_s), fontsize=6)
        for ti, tool in enumerate(sorted(sub[sub.tool != "sklib"].tool.unique())):
            t = sub[sub.tool == tool]
            ax.scatter(t.bits_per_kmer, t.throughput_Mkmer_s, marker=marks[1 + ti % 6],
                       color=cmap(di), s=80, edgecolor="k", label=f"{tool} {ds}")
    ax.set_xlabel("bits / k-mer  (smaller = more compact)")
    ax.set_ylabel(f"query throughput Mkmer/s  ({workload}, threads={tmax})")
    ax.set_title("Fig 5 — Space / speed Pareto"); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    save(fig, out_dir, "fig5_pareto")


def fig6_envelope(df, out_dir):
    """Correctness over the (k,m) grid for sklib (valid-envelope / limitations view)."""
    g = df[(df.phase == "query") & (df.workload == "positive") & (df.tool == "sklib")]
    g = g.dropna(subset=["k", "m"])
    g = g[g.correctness.isin(["pass", "fail"])]
    if g.empty:
        print("  [fig6] no positive-workload correctness rows -> skip"); return
    datasets = sorted(g.dataset.unique())
    fig, axes = plt.subplots(1, len(datasets), figsize=(4 * len(datasets), 3.6), squeeze=False)
    for ax, ds in zip(axes[0], datasets):
        sub = g[g.dataset == ds]
        piv = sub.assign(ok=(sub.correctness == "pass").astype(int)) \
                 .pivot_table(index="k", columns="m", values="ok", aggfunc="min")
        im = ax.imshow(piv.values, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels([int(x) for x in piv.columns])
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels([int(x) for x in piv.index])
        ax.set_xlabel("m"); ax.set_ylabel("k"); ax.set_title(ds)
    fig.suptitle("Fig 6 — Query correctness over (k,m)  [green=pass]")
    save(fig, out_dir, "fig6_envelope")


def fig7_ksweep(df, out_dir):
    """k-sweep: bits/k-mer and streaming throughput vs k, one line per tool, per dataset."""
    c = df[df.phase == "construct"].dropna(subset=["bits_per_kmer", "k"])
    if c.empty or c.k.nunique() < 2:
        print("  [fig7] need >=2 k values -> skip"); return
    # best (min bits) construct per (tool,dataset,k); collapses any m-sweep to its best point
    cc = c.loc[c.groupby(["tool", "dataset", "k"]).bits_per_kmer.idxmin()]
    q = df[(df.phase == "query") & (df.workload == "streaming")].dropna(subset=["throughput_Mkmer_s", "k"])
    if not q.empty:
        q = q[q.threads == q.threads.max()]
        qq = q.groupby(["tool", "dataset", "k"]).throughput_Mkmer_s.max().reset_index()
    datasets = sorted(cc.dataset.unique()); tools = sorted(cc.tool.unique())
    marks = _markers()
    fig, axes = plt.subplots(2, len(datasets), figsize=(4.5 * len(datasets), 7), squeeze=False)
    for j, ds in enumerate(datasets):
        for ti, tool in enumerate(tools):
            s = cc[(cc.dataset == ds) & (cc.tool == tool)].sort_values("k")
            if not s.empty:
                axes[0][j].plot(s.k, s.bits_per_kmer, marker=marks[ti % len(marks)], label=tool)
            if not q.empty:
                sq = qq[(qq.dataset == ds) & (qq.tool == tool)].sort_values("k")
                if not sq.empty:
                    axes[1][j].plot(sq.k, sq.throughput_Mkmer_s, marker=marks[ti % len(marks)], label=tool)
        axes[0][j].set_title(ds); axes[0][j].set_xlabel("k"); axes[0][j].set_ylabel("bits / k-mer")
        axes[1][j].set_xlabel("k"); axes[1][j].set_ylabel("streaming Mkmer/s")
        for ax in (axes[0][j], axes[1][j]):
            ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle("Fig 7 — k-sweep: space (top) and streaming query speed (bottom) vs k")
    save(fig, out_dir, "fig7_ksweep")


def main():
    csv = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(__file__), "..", "out", "bench", "results.csv")
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(csv), "figures")
    os.makedirs(out_dir, exist_ok=True)
    if not os.path.exists(csv):
        sys.exit(f"results CSV not found: {csv}")
    df = load(csv)
    print(f"loaded {len(df)} rows from {csv}; writing figures to {out_dir}")
    for fn in (fig1_space, fig2_compaction, fig3_construction, fig4_query, fig5_pareto, fig6_envelope, fig7_ksweep):
        try:
            fn(df, out_dir)
        except Exception as e:  # one bad figure shouldn't sink the rest
            print(f"  [{fn.__name__}] error: {e}")


if __name__ == "__main__":
    main()
