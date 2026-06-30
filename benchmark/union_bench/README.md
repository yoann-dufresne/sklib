# union_bench — isolated set-op microbenchmark

A small standalone harness that benchmarks and verifies sklib's in-memory set operations
**directly**, with none of the CLI / FASTA-parsing / construction noise. It `#include`s the real
`SetOperations.hpp` / `VirtualSkmer.hpp` / `Skmer.hpp` and calls the actual `set_union` /
`materialize_setop`, so every edit to those headers is measured as-is. Single-threaded by design
(it isolates the per-element merge cost, not parallel scaling).

> **The name is historical.** It started as a *union*-only harness; it now covers
> **union / intersection / diff / xor** via `--op` (all share `materialize_setop`). The directory,
> the `WITH_UNION_BENCH` CMake flag and the driver scripts keep the old name for continuity.

## Build & run

Off by default — enable with the CMake option, then drive it through the script:

```bash
cmake -S . -B build-union -DWITH_UNION_BENCH=ON && cmake --build build-union -j
# prep inputs + verify (content-equivalence vs a frozen reference) + bench (median/min/stddev/MAD):
bash benchmark/scripts/microbench/union_bench.sh all
# A/B two binaries on one pair, interleaved to cancel thermal drift:
benchmark/scripts/microbench/union_ab.sh <binA> <binB> A.sskm B.sskm [rounds]
```

The driver builds the A/B input lists (via `mutate.py` at a target Jaccard) and a **frozen**
reference output, caches everything under the git-ignored `benchmark/results/union_bench/`, and
writes per-iteration TSVs there. It is the data source behind
[`../results/journals/UNION_SPEEDUP.md`](../results/journals/UNION_SPEEDUP.md) and
[`../results/journals/XOR_SPEEDUP.md`](../results/journals/XOR_SPEEDUP.md).

## Files

- `union_bench.cpp` — the harness (`--mode bench` / `--mode verify`, `--op union|intersection|diff|xor`).
- `CMakeLists.txt` — builds the `union_bench` target only when `-DWITH_UNION_BENCH=ON`.
