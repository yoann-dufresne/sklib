# Supplementary material — SPIRE 2026 paper

Supplementary material for **"Binary Search and Set Operations on Compacted
k-mer Lists"** (Y. Dufresne and F. Andreace), to appear in the proceedings of
[SPIRE 2026](https://www.springer.com/series/558) (LNCS, Springer).

The proceedings version carries no appendix, so the material it cites as
"supplementary material, Sect. A" through "Sect. J" lives here. The section
lettering below matches those citations exactly.

📄 **[supplementary.pdf](supplementary.pdf)**

| Sect. | Content |
|---|---|
| A | Compared tools (representation, *k* range, query and set-op support) |
| B | Set operations: thread scaling on chr1, and cost across operations and overlap |
| C | Worked example of the merge-like scan computing A ∩ B |
| D | Cost of materializing a set operation: phase breakdown, compacted vs uncompacted output, closure check |
| E | Implementation optimizations: bucketing and quotienting |
| F | Canonical k-mers and minimizer position |
| G | Record layout and space accounting |
| H | Construction space: bits per k-mer, all tools |
| I | Experimental setup, mutation model, tool versions and build flags |
| J | Raw data and reproduction: where every reported number comes from |

Section J maps each result to its CSV and script under [`benchmark/`](../../benchmark/);
the campaign write-ups are in [`benchmark/results/journals/`](../../benchmark/results/journals/).

## Rebuilding the PDF

`supplementary.bbl` is committed, so only Springer's LNCS class file is needed:
download `llncs.cls` from the [LNCS authors' page](https://www.springer.com/gp/computer-science/lncs/conference-proceedings-guidelines),
drop it next to `supplementary.tex`, then

```bash
latexmk -pdf supplementary.tex
```

## Citing

Please cite the preprint — see [Citation](../../README.md#citation) in the top-level README.
