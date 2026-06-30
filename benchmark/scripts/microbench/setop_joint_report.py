#!/usr/bin/env python3
"""Summarise setop_joint_compare.csv: sklib vs KMC on the combined (joint) 4-op
materialization. Prints a markdown table per (k, threads, Jaccard) with time and
peak-RSS ratios (sklib/kmc; <1 means sklib faster / leaner) and a verdict."""
import csv, sys
from collections import defaultdict

path = sys.argv[1] if len(sys.argv) > 1 else "benchmark/results/latest/setop_joint_compare.csv"
rows = list(csv.DictReader(open(path)))
# key: (k, threads, jaccard_target) -> tool -> row
cell = defaultdict(dict)
for r in rows:
    if r["op"] != "joint" or r["mode"] != "materialize":
        continue
    cell[(int(r["k"]), int(r["threads"]), r["jaccard_target"])][r["tool"]] = r

def f(x):
    try: return float(x)
    except: return float("nan")

print("# KMC vs sklib — opération combinée (joint 4-op), matérialisation, chr21\n")
print("Ratio = sklib / kmc.  **<1 ⇒ sklib plus rapide / plus sobre** ; >1 ⇒ KMC meilleur.\n")
print("| k | t | J | sklib s | kmc s | **ratio t** | sklib RSS Mo | kmc RSS Mo | ratio RSS | result_kmers |")
print("|--:|--:|--:|--:|--:|:--:|--:|--:|:--:|--:|")
faster = same_work = total = 0
for key in sorted(cell):
    k, t, J = key
    s, m = cell[key].get("sklib"), cell[key].get("kmc")
    if not (s and m):
        continue
    total += 1
    rt = f(s["time_s"]) / f(m["time_s"]) if f(m["time_s"]) else float("nan")
    rr = f(s["peak_rss_mb"]) / f(m["peak_rss_mb"]) if f(m["peak_rss_mb"]) else float("nan")
    if rt < 1: faster += 1
    if s["result_kmers"] == m["result_kmers"]: same_work += 1
    flag = "✅" if rt < 1 else "❌"
    print(f"| {k} | {t} | {J} | {s['time_s']} | {m['time_s']} | {flag} {rt:.2f}× | "
          f"{s['peak_rss_mb']} | {m['peak_rss_mb']} | {rr:.2f}× | {s['result_kmers']} |")

print(f"\n**Bilan** : sklib plus rapide que KMC dans **{faster}/{total}** cellules. "
      f"Cardinalités identiques (même travail) dans {same_work}/{total} cellules.")
