#!/usr/bin/env python3
"""Convert a callgrind output file into FlameGraph "folded stack" lines.

Reads a callgrind dump (events: Ir) and reconstructs an inclusive call tree by
proportional attribution, then emits `frame1;frame2;... <self_cost>` lines that
FlameGraph/flamegraph.pl turns into an SVG. Handles callgrind name compression
(`fn=(ID) name` then `fn=(ID)`).

Attribution: each function's self cost is split among the call paths that reach
it, proportionally to the inclusive cost flowing through each path. Recursion is
broken by attributing the back-edge as self cost at the recursion point so totals
are conserved and the graph stays finite.

Usage: callgrind_to_folded.py cg.out [--phases] > folded.txt
       --phases  print an inclusive per-phase table to stderr instead of folding
"""
import sys
import re
from collections import defaultdict


def parse(path):
    self_cost = defaultdict(int)
    edges = defaultdict(lambda: defaultdict(int))
    fn_names = {}   # id -> name
    cur_fn = None
    pending_callee = None

    def resolve(spec, table):
        # spec is the text after 'fn='/'cfn='; forms: "(ID) name", "(ID)", "name"
        m = re.match(r"^\((\d+)\)\s*(.*)$", spec)
        if m:
            idx, name = m.group(1), m.group(2)
            if name:
                table[idx] = name
                return name
            return table.get(idx, "(" + idx + ")")
        return spec

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("fn="):
                cur_fn = resolve(line[3:], fn_names)
                pending_callee = None
            elif line.startswith("cfn="):
                pending_callee = resolve(line[4:], fn_names)
            elif line.startswith("calls="):
                pass
            elif line[0] in "0123456789+-*":
                # cost line: "<position spec> <Ir>"; position can be absolute,
                # relative (+N/-N), or '*'. Events: Ir => the cost is the last token.
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    ir = int(parts[-1])
                except ValueError:
                    continue
                if pending_callee is not None:
                    edges[cur_fn][pending_callee] += ir
                    pending_callee = None
                elif cur_fn is not None:
                    self_cost[cur_fn] += ir
    fns = set(self_cost) | set(edges)
    for c in edges:
        fns.update(edges[c])
    inclusive = {fn: self_cost.get(fn, 0) + sum(edges[fn].values()) for fn in fns}
    return self_cost, edges, inclusive, fns


def short(fn):
    fn = re.sub(r"\s*\[/.*?\]", "", fn)
    fn = fn.replace("km::sortedlist::", "").replace("km::", "")
    fn = fn.replace("SortedVirtualSkmerList<unsigned long>::", "")
    fn = fn.replace("SkmerManipulator<unsigned long>::", "")
    return fn.replace(";", ":")


def main():
    path = sys.argv[1]
    phases = "--phases" in sys.argv[2:]
    self_cost, edges, inclusive, fns = parse(path)

    if phases:
        total = sum(self_cost.values())
        def find_incl(sub):
            return max((inclusive[k] for k in inclusive if sub in k), default=0)
        rows = [
            ("index load",                 "Serializer<unsigned long>::load"),
            ("enumeration (File++)",       "FileSkmerator<unsigned long>::Iterator::operator++"),
            ("  canonicalize",             "::canonicalize"),
            ("  reverse_complement",       "::reverse_complement"),
            ("  add_nucleotide",           "::add_nucleotide"),
            ("search query_skmer",         "query_skmer(Skmer"),
            ("  find_closest_valid_skmer", "find_closest_valid_skmer"),
            ("  has_valid_kmer",           "has_valid_kmer"),
            ("  kmer_compare",             "kmer_compare"),
            ("output print_query_results", "print_query_results"),
        ]
        sys.stderr.write(f"total self Ir = {total:,}\n")
        for name, sub in rows:
            v = find_incl(sub)
            sys.stderr.write(f"  {name:<30} {v:>14,}  {100*v/total:5.1f}%\n")
        return

    root_arg = None
    for a in sys.argv[2:]:
        if a.startswith("--root="):
            root_arg = a[len("--root="):]

    if root_arg:
        roots = sorted((fn for fn in fns if root_arg in fn),
                       key=lambda x: -inclusive.get(x, 0))[:1]
    else:
        callees = set()
        for c in edges:
            callees.update(edges[c])
        roots = [fn for fn in fns if fn not in callees and inclusive.get(fn, 0) > 0]
    folded = defaultdict(int)
    sys.setrecursionlimit(1000000)

    def walk(fn, stack, flow):
        inc = inclusive.get(fn, 0)
        if inc <= 0 or flow <= 0:
            return
        frac = flow / inc
        key = ";".join(stack)
        s = int(self_cost.get(fn, 0) * frac)
        if s > 0:
            folded[key] += s
        for callee, ecost in edges[fn].items():
            if callee in stack:
                folded[key] += int(ecost * frac)
                continue
            walk(callee, stack + [short(callee)], ecost * frac)

    for r in sorted(roots, key=lambda x: -inclusive[x]):
        walk(r, [short(r)], inclusive[r])

    for k, v in folded.items():
        if v > 0:
            print(f"{k} {v}")


if __name__ == "__main__":
    main()
