# Downstream clustering, marker, and velocity fidelity verdict

**Date:** 2026-08-31  
**Preregistered thresholds:** commit `d3e9b3a`  
**Frozen protocol:** commit `07d37b4`  
**Primary result:** `results/post-v1-downstream-fidelity.json`  
**Post-hoc localization:** `results/post-v1-downstream-cluster-diagnostics.json`

## Verdict

The preregistered human Gene downstream gate **FAILS**: D0 passes, but D1 and D2′ miss the ARI,
NMI, and median-neighborhood thresholds. The mouse stress test also fails its partition and
pseudobulk thresholds. These misses are retained as failures; no threshold was changed after
outcomes were observed.

The failure is narrower than a general loss of biological signal. Across D0/D1/D2′, mean top-25
marker overlap is 0.962–0.987, median marker-effect Pearson is 0.9988–0.9997, and pseudobulk
Pearson is 0.99996–0.99999. D5 has perfect mean top-25 overlap and 0.9992 marker-effect Pearson,
although its pseudobulk Pearson/Spearman (0.9967/0.9980) miss the frozen 0.999 thresholds.

The separately defined velocity consequence gate **PASSES all four datasets**. Median projected
steady-state velocity-vector cosine is 0.9931–0.9986; every cell has positive oracle/replay
cosine; median top-3 transition Jaccard is 1.0 and mean Jaccard is 0.873–0.914. This does not erase
the component-level MARGINAL verdict for mouse velocity: it says the measured S/U deviations do
not materially change this fixed downstream consequence probe.

| dataset | Gene ARI | NMI | median 15-NN Jaccard | mean marker top-25 Jaccard | pseudobulk Pearson | velocity cosine | Gene verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| D0 | 0.952 | 0.959 | 0.875 | 0.981 | 0.999986 | 0.9986 | PASS |
| D1 | 0.834 | 0.901 | 0.765 | 0.987 | 0.999963 | 0.9969 | FAIL |
| D2′ | 0.619 | 0.745 | 0.765 | 0.962 | 0.999993 | 0.9967 | FAIL |
| D5 mouse 5′ | 0.609 | 0.680 | 0.875 | 1.000 | 0.996728 | 0.9931 | FAIL |

## Sensitivity and deterministic replication

The structured 20% cross-cluster replacement control is detected on every confirmatory dataset:
its ARI is 0.485–0.568 and its median neighbor Jaccard is 0.579, at least 0.05 below replay as
frozen. Thus the metric panel is not simply incapable of detecting downstream disruption.

The complete four-dataset run was repeated from the clean protocol commit. Every non-timing JSON
field is identical. Primary/replicate external wall times are 28.26/27.55 s and peak RSS is
3,014,024/3,086,512 KiB. Primary result SHA-256 is
`06e017cda86a585aa9b7463a4f1a914a60af54bb0b3e0edfabf33311bd24196f`.

## Failure localization

Post-hoc diagnostics preserve the frozen FAIL and ask where the partition changes occur.
Fixed-oracle-centroid label agreement is 99.84% (D0), 99.09% (D1), 95.09% (D2′), and 98.85%
(D5). The independent Leiden differences are dominated by split/merge boundaries within coherent
compartments: an NK cluster in D0, an HLA/S100 myeloid cluster in D1, related neuronal clusters in
D2′, and related B-cell clusters in D5. Mismatched cells do not show materially larger per-cell
total or moved-mass deviations than stable cells. This supports a precise conclusion: marker,
pseudobulk, reference-mapping, and the tested velocity consequences are stable, but independently
recomputed high-resolution Leiden partitions are not invariant to the small count perturbation.

## Paper consequence

Do not write that small UMI-mass deviation guarantees identical clustering. Report the failed
partition gate alongside the marker and velocity passes. The result strengthens the paper by
showing where conditional sufficiency ends: the stored evidence is sufficient for the tested
count and consequence policies, while a discontinuous downstream algorithm can amplify small
fresh-alignment differences near partition boundaries.
