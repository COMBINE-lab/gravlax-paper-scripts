# Storing the quotient, not the labels: UMI adjacency as a sufficient statistic

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (66,601,887 reads, 1,225 cells, 70% saturation)
**Code:** `src/crates/aie/src/graph.rs` (`aie umi-graph`), `build --no-umi-collapse`
**Results:** `results/umigraph-{full,ingest,full-policies}.json`, `results/archive-full-graphbase.json`

## Result: the archive gets **2.05× smaller and 4.05× more accurate at the same time**

This is not a trade-off. Replacing stored UMI *values* with the UMI *adjacency graph* halves the
archive and quadruples grouping fidelity, because the graph defers merge decisions to replay
instead of freezing one at ingest.

| Archive variant | Size | vs CRAM 3.1 | vs Malva | vs FASTQ | Grouping error |
|---|---:|---:|---:|---:|---:|
| E1 + raw UMIs (Gate G-E) | 119.5 MB | 11.76× | 3.29× | 42.2× | 0.939% |
| **E1 no-UMI + adjacency graph** | **58.4 MB** | **24.06×** | **6.73×** | **86.3×** | **0.232%** |

Graph-variant breakdown: 54.93 MB base (no UMI stream) + **3.46 MB graph**.

## The argument

Replay never reads a UMI value. Values decide exactly two things: whether two molecules in a
(cell, gene) share a UMI, and whether one sits within a substitution of a more abundant one. Both are
properties of a *relation over molecules*, not per-molecule payload — and the relation is sparse,
because UMI space (4¹² ≈ 16.7M) dwarfs per-cell molecule count (~10⁴).

Measured, on the annotation-free ingest at full scale:

| Window | same-value links | 1-substitution edges |
|---|---:|---:|
| 10 kb | 349,787 | 1,076,634 |
| 100 kb | 426,116 | 1,078,540 |
| 1 Mb | 468,094 | 1,080,072 |
| 3 Mb | 476,687 | 1,081,023 |
| chromosome-wide | 529,444 | 1,098,826 |

Across **23.1M base molecules** the entire graph is ~1.56M relations. Encoded and zstd-compressed it
is **3.46 MB = 1.20 bits/molecule, against 27.00 bits/molecule for raw values — 19.4×.** On the
oracle BAM the figures are 3.69 MB and 1.34 bits/molecule (18.2×).

**The structure is intrinsically local**, which is what makes a bounded window safe: widening from
10 kb to chromosome-wide adds only 2% more edges. A 3 Mb window therefore forecloses essentially
nothing — it comfortably exceeds any human gene span, so replay can serve any annotation whose
genes fit inside it.

## Sufficiency, and a prediction of mine that was wrong

I predicted graph replay would be **exact**. It is not — it is 0.232%, not 0%. The first run was
worse still (1.05%) and I recorded that before diagnosing it. The diagnosis held: the residual was
policy mismatch, not missing information. STARsolo's oracle runs `--soloUMIfiltering
MultiGeneUMI_CR`, which keeps a multi-gene UMI only for its best-supported gene; our replay was
counting such UMIs in every gene they touched. There are **82,260 multi-gene (cell, UMI) pairs of
11,207,287** — 0.73%, small in count but concentrated where it matters.

| 1MM policy | MultiGeneUMI_CR | median per-cell L1 | p90 |
|---|---|---:|---:|
| strict (Gate G-D greedy) | no | 1.9134% | 2.3053% |
| tie-merging (CellRanger-like) | no | 1.0490% | 1.4052% |
| strict | yes | 1.0885% | 1.3133% |
| **tie-merging** | **yes** | **0.2318%** | **0.3421%** |
| *Gate G-D frozen grouping (reference)* | — | *0.9390%* | *1.1388%* |

Two policy details each worth ~1 percentage point, and they compose: **tie-merging alone 1.05%,
the filter alone 1.09%, both together 0.23%.** Gate G-D's residual 0.939% was therefore never
evidence about information sufficiency — it was our 1MM greedy declining to merge equally-abundant
neighbours, plus the missing multi-gene filter. That reading of G-D should be corrected in the paper.

The remaining 0.232% is not yet explained. Candidates, in order of suspicion: CellRanger's 1MM_CR
consults **base qualities** when choosing an absorber and we ignore them; and our base molecules use
2 kb single-linkage loci, so a gene whose reads split across a wider gap contributes two nodes where
STARsolo sees one. Both are testable and neither implies missing information.

## Why this is a method contribution, not a coding trick

The archive stores a **quotient structure** rather than the labels that induce it. That reframing
does three things at once:

1. **Compression** — 19.4× on what was 56% of the archive.
2. **Fidelity** — merge decisions stay open until an annotation is supplied, so the archive is
   strictly more replayable than one with a frozen grouping. This is why size and accuracy improve
   together rather than trading off.
3. **A theorem-shaped claim** with an explicit annotation-class parameter, in the spirit of project
   plan §22: *the UMI-adjacency graph restricted to a window W, together with per-molecule read
   counts, is a sufficient statistic for reproducing 1MM-style UMI collapse under every annotation
   whose genes span ≤ W.* The window ladder above is the empirical support: the relation is local,
   so bounded W costs almost nothing.

## Consequences for the project

- **Gate G-E should be restated** at 24.06× vs CRAM 3.1 (was 11.76×), still with UMIs' full
  re-collapse capability retained rather than frozen.
- **Gate G-D's error should be restated** at 0.232% (was 0.939%), and its residual re-attributed
  from "irreducible annotation-dependence" to "unmatched UMI policy". The error budget improves:
  grouping error drops well below the alignment ceiling (1.13%), making the **alignment ceiling
  unambiguously the binding constraint** against the 2.43% annotation signal.
- The base molecule count rises (19.9M → 21.2M) because nothing is collapsed at ingest. That cost is
  already included in the 58.4 MB figure.

## Caveats

- One dataset, one depth. The graph's sparsity is a birthday-problem consequence of molecules per
  cell versus UMI space, so a deeper library or a longer UMI changes it predictably — but it should
  be measured, not extrapolated, given how badly dev-scale extrapolation failed for Gate G-D.
- The 0.232% is measured against STARsolo only. The alevin-fry oracle family uses different UMI
  resolution (including transcript-level ambiguity) and will need its own sufficiency test.
- The graph encoding here is deliberately plain (varint molecule-index pairs). It is a price tag,
  not a proposed format; a real design would fold it into the block structure the query work needs.
