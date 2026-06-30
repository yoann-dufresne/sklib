#!/usr/bin/env python3
"""Analyse construct_scaling.csv: median per (genome,threads), speedup, Amdahl fit, phase split,
Phase-2 efficiency, and a CPU%-derived serial-fraction cross-check. Prints a text report and
writes two PNGs (speedup-vs-threads + Amdahl fit; stacked phase1/phase2 wall vs threads).

  python3 benchmark/scripts/profiling/diag_plot.py [csv] [outdir]
"""
import sys, csv, math
from collections import defaultdict

CSV = sys.argv[1] if len(sys.argv) > 1 else "benchmark/results/latest/construct_scaling.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "benchmark/results/latest/figs"

def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0: return float("nan")
    return xs[n//2] if n % 2 else 0.5*(xs[n//2-1]+xs[n//2])

# rows[genome][threads] = list of dict(wall,p1,p2,cpu,rss)
rows = defaultdict(lambda: defaultdict(list))
nbp = {}
order = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        g = r["genome"]; t = int(r["threads"])
        if g not in nbp:
            nbp[g] = int(r["n_bp"]); order.append(g)
        def fnum(x):
            try: return float(x)
            except: return float("nan")
        rows[g][t].append(dict(wall=fnum(r["wall_s"]), p1=fnum(r["phase1_s"]),
                               p2=fnum(r["phase2_s"]), cpu=fnum(r["cpu_pct"]),
                               rss=fnum(r["peak_rss_kb"])))

def amdahl_fit(ts, sp):
    """Least-squares serial fraction s for speedup(t)=1/(s+(1-s)/t). Grid + refine, no scipy."""
    def sse(s):
        return sum((1.0/(s+(1-s)/t) - y)**2 for t, y in zip(ts, sp))
    lo, hi = 1e-4, 0.999
    for _ in range(60):
        cands = [lo + (hi-lo)*i/20 for i in range(21)]
        best = min(cands, key=sse)
        step = (hi-lo)/20
        lo, hi = max(1e-4, best-step), min(0.999, best+step)
    return 0.5*(lo+hi)

summary = []
print(f"# Construction scaling diagnostic\n\nsource: {CSV}\n")
for g in order:
    ts = sorted(rows[g])
    med = {t: dict(wall=median([x["wall"] for x in rows[g][t]]),
                   p1=median([x["p1"] for x in rows[g][t]]),
                   p2=median([x["p2"] for x in rows[g][t]]),
                   cpu=median([x["cpu"] for x in rows[g][t]]),
                   rss=median([x["rss"] for x in rows[g][t]])) for t in ts}
    w1 = med[ts[0]]["wall"]  # threads==min (==1 expected)
    p2_1 = med[1]["p2"] if 1 in med else med[ts[0]]["p2"]
    p1_1 = med[1]["p1"] if 1 in med else med[ts[0]]["p1"]
    sp = [w1/med[t]["wall"] for t in ts]
    s_fit = amdahl_fit(ts, sp)
    s_timer = p1_1/(p1_1+p2_1) if (p1_1+p2_1) > 0 else float("nan")
    Mbp = nbp[g]/1e6
    print(f"## {g}  ({Mbp:.1f} Mbp)")
    print(f"{'thr':>4} {'wall_s':>8} {'speedup':>8} {'phase1_s':>9} {'phase2_s':>9} "
          f"{'p2_eff':>7} {'cpu%':>6} {'cpu_s':>7} {'rss_MB':>7}")
    for t in ts:
        m = med[t]
        eff = (p2_1/(t*m["p2"])) if m["p2"] > 0 else float("nan")
        cpu_s = m["cpu"]/100.0*m["wall"] if m["cpu"]==m["cpu"] else float("nan")
        print(f"{t:>4} {m['wall']:>8.3f} {w1/m['wall']:>8.2f} {m['p1']:>9.3f} {m['p2']:>9.3f} "
              f"{eff:>7.2f} {m['cpu']:>6.0f} {cpu_s:>7.2f} {m['rss']/1024:>7.0f}")
    # CPU cross-check. At t=1 total CPU-seconds == wall == W (= p1 + p2@1), so a CPU-based serial
    # fraction p1/W independently corroborates the timer split. The high-t plateau then exposes
    # Phase-2 parallel overhead: if Phase 2 were overhead-free, total CPU-s would stay == W; any
    # excess is wasted work (memory-bandwidth stalls / hybrid E-cores / sync), shown as inflation.
    thi = ts[-1]
    W_cpu = med[1]["cpu"]/100.0*med[1]["wall"]
    s_cpu = p1_1/W_cpu if W_cpu > 0 else float("nan")
    ph2_cpu_hi = med[thi]["cpu"]/100.0*med[thi]["wall"] - p1_1   # CPU-s spent in Phase 2 at top -t
    ph2_cpu_1  = W_cpu - p1_1
    infl = ph2_cpu_hi/ph2_cpu_1 if ph2_cpu_1 > 0 else float("nan")
    avg_cores_hi = ph2_cpu_hi/med[thi]["p2"] if med[thi]["p2"] > 0 else float("nan")
    print(f"\n  serial fraction s:  timer(p1/(p1+p2@1))={s_timer:.3f}   "
          f"Amdahl-fit={s_fit:.3f}   CPU@t1-derived={s_cpu:.3f}")
    print(f"  Amdahl ceiling 1/s: timer={1/s_timer:.2f}x   fit={1/s_fit:.2f}x   "
          f"| measured max speedup={max(sp):.2f}x @ t={ts[sp.index(max(sp))]}")
    print(f"  Phase-2 @t{thi}: eff={p2_1/(thi*med[thi]['p2']):.2f}, ~{avg_cores_hi:.1f} cores busy, "
          f"CPU-work inflation {infl:.2f}x vs t1\n")
    summary.append((g, Mbp, s_timer, s_fit, 1/s_fit, max(sp)))

print("## Summary across genomes")
print(f"{'genome':>10} {'Mbp':>7} {'s_timer':>8} {'s_fit':>7} {'ceiling':>8} {'max_obs':>8}")
for g, Mbp, st, sf, ceil, mx in summary:
    print(f"{g:>10} {Mbp:>7.1f} {st:>8.3f} {sf:>7.3f} {ceil:>7.2f}x {mx:>7.2f}x")

# ---- plots ----
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import os
    os.makedirs(OUT, exist_ok=True)
    # speedup vs threads + Amdahl
    fig, ax = plt.subplots(figsize=(7,5))
    for g in order:
        ts = sorted(rows[g]); med_wall = {t: median([x["wall"] for x in rows[g][t]]) for t in ts}
        w1 = med_wall[ts[0]]; sp = [w1/med_wall[t] for t in ts]
        s = amdahl_fit(ts, sp)
        line, = ax.plot(ts, sp, "o-", label=f"{g} ({nbp[g]/1e6:.0f} Mbp), 1/s={1/s:.1f}x")
        tt = [ts[0]+ (ts[-1]-ts[0])*i/50 for i in range(51)]
        ax.plot(tt, [1/(s+(1-s)/t) for t in tt], "--", color=line.get_color(), alpha=0.5)
    ax.plot(sorted({t for g in order for t in rows[g]}),
            sorted({t for g in order for t in rows[g]}), ":", color="gray", label="linear ideal")
    ax.set_xlabel("threads (-t)"); ax.set_ylabel("speedup vs -t1")
    ax.set_title("sskm construct — multithread scaling (dashed = Amdahl fit)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{OUT}/construct_speedup.png", dpi=120)
    # stacked phase1/phase2 vs threads (one genome per subplot)
    n = len(order)
    fig, axs = plt.subplots(1, n, figsize=(3.2*n, 4), sharey=False)
    if n == 1: axs = [axs]
    for ax, g in zip(axs, order):
        ts = sorted(rows[g])
        p1 = [median([x["p1"] for x in rows[g][t]]) for t in ts]
        p2 = [median([x["p2"] for x in rows[g][t]]) for t in ts]
        ax.bar(range(len(ts)), p1, label="phase1 (serial producer)", color="#d62728")
        ax.bar(range(len(ts)), p2, bottom=p1, label="phase2 (parallel)", color="#1f77b4")
        ax.set_xticks(range(len(ts))); ax.set_xticklabels(ts)
        ax.set_title(f"{g} ({nbp[g]/1e6:.0f} Mbp)"); ax.set_xlabel("-t")
    axs[0].set_ylabel("wall time (s)"); axs[0].legend(fontsize=8)
    fig.suptitle("Phase-1 (serial) stays flat; only Phase-2 (parallel) shrinks with -t")
    fig.tight_layout(); fig.savefig(f"{OUT}/construct_phase_split.png", dpi=120)
    print(f"\nplots: {OUT}/construct_speedup.png  {OUT}/construct_phase_split.png")
except Exception as e:
    print(f"\n[plot skipped: {e}]")
