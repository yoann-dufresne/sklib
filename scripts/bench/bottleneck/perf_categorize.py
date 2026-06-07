#!/usr/bin/env python3
# Aggregate `perf report --stdio --no-children -g none` self-time % into phase buckets.
import re, subprocess, sys

data = sys.argv[1]
raw = subprocess.run(["perf","report","-i",data,"--stdio","--no-children","-g","none"],
                     capture_output=True, text=True).stdout

# (category, regex on demangled symbol). First match wins; order matters.
RULES = [
    ("recompact: colinear_chaining", r"colinear_chaining"),
    ("recompact: get_candidate_overlaps", r"get_candidate_overlaps"),
    ("recompact: sort_column",      r"sort_column"),
    ("recompact: merge_LList_column",r"merge_LList_column"),
    ("recompact: gen_sorted_list",  r"generate_sorted_list_from_enumeration"),
    ("recompact: Virtual_skmer/vec",r"Virtual_skmer|emplace_back|_M_realloc"),
    ("std::sort (introsort/insert)",r"__introsort_loop|__insertion_sort|__final_insertion|std::sort|__heap"),
    ("std::lower_bound/binsearch",  r"lower_bound|__lower_bound"),
    ("skmer: get_skmer_of_kmer",    r"get_skmer_of_kmer"),
    ("skmer: extract_*",            r"extract_kmer|extract_prefix_suffix|extract_nucleotide"),
    ("skmer: add/clean/masks",      r"add_nucleotide|clean_nucleotide|generate_masks|init_skmer|SkmerManipulator"),
    ("skmer: kmer_compare",         r"kmer_compare"),
    ("skmer: has_valid_kmer/bounds",r"has_valid_kmer|get_valid_kmer_bounds"),
    ("skmer: pair shift/cmp ops",   r"Skmer<[^>]*>::pair::operator|::pair::"),
    ("merge driver (materialize/merge_columns/set_next_valid/CountSink/CollectSink)",
                                    r"materialize_setop|merge_columns|set_next_valid|set_sizes|CountSink|CollectSink|intersection|union|difference|diff_size|union_size"),
    ("I/O read (bucket load)",      r"bucket|load_bucket|ifstream|_M_extract|filebuf|::read|sgetn|memcpy.*read"),
    ("I/O write (append/serialize)",r"append_payload|patch_|write_header|ofstream|filebuf.*write|VirtualSkmerSerializer"),
    ("alloc: malloc/free/new",      r"\bmalloc\b|_int_malloc|_int_free|\bfree\b|operator new|operator delete|cfree|tcache|arena|brk|mmap|munmap|sbrk"),
    ("libc memmove/memcpy/memset",  r"\bmemmove\b|\bmemcpy\b|__memmove|__memcpy|\bmemset\b|__memset"),
    ("kernel/[unknown]",            r"^\[k\]|\[unknown\]|0x[0-9a-f]+$|^0xffff"),
]

def demangle(name):
    return name  # perf already demangles with c++filt available

agg = {}
total = 0.0
unmatched = []
for line in raw.splitlines():
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    m = re.match(r'^([\d.]+)%\s+\S+\s+(\S+)\s+\[.\]\s+(.*)$', line)
    if not m:
        # kernel symbols: [k]
        m = re.match(r'^([\d.]+)%\s+\S+\s+(\S+)\s+\[k\]\s+(.*)$', line)
        if not m:
            continue
        pct = float(m.group(1)); sym = "[k] "+m.group(3)
    else:
        pct = float(m.group(1)); sym = m.group(3)
    total += pct
    placed = False
    for cat, rx in RULES:
        if re.search(rx, sym):
            agg[cat] = agg.get(cat, 0.0) + pct
            placed = True
            break
    if not placed:
        agg.setdefault("OTHER", 0.0)
        agg["OTHER"] += pct
        unmatched.append((pct, sym))

print(f"# {data}   (sum of self% = {total:.1f})")
# group recompact vs merge vs skmer-helpers vs io vs alloc
GROUPS = {
    "RE-COMPACTION (rebuild super-k-mers)": ["recompact: colinear_chaining","recompact: get_candidate_overlaps","recompact: sort_column","recompact: merge_LList_column","recompact: gen_sorted_list","recompact: Virtual_skmer/vec","std::sort (introsort/insert)","std::lower_bound/binsearch"],
    "SKMER HELPERS": ["skmer: get_skmer_of_kmer","skmer: extract_*","skmer: add/clean/masks","skmer: kmer_compare","skmer: has_valid_kmer/bounds","skmer: pair shift/cmp ops"],
    "MERGE DRIVER": ["merge driver (materialize/merge_columns/set_next_valid/CountSink/CollectSink)"],
    "I/O": ["I/O read (bucket load)","I/O write (append/serialize)"],
    "ALLOC/MEM": ["alloc: malloc/free/new","libc memmove/memcpy/memset"],
    "KERNEL/OTHER": ["kernel/[unknown]","OTHER"],
}
for g, cats in GROUPS.items():
    gtot = sum(agg.get(c,0.0) for c in cats)
    if gtot < 0.01: continue
    print(f"\n== {g}: {gtot:.1f}%")
    for c in cats:
        if agg.get(c,0.0) >= 0.05:
            print(f"   {agg[c]:6.2f}%  {c}")
if unmatched:
    print("\n# top unmatched:")
    for pct, sym in sorted(unmatched, reverse=True)[:8]:
        print(f"   {pct:6.2f}%  {sym[:120]}")
