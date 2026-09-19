# Stable macro-partition retention after the frozen Leiden stress test

**Verdict date:** 2026-08-31  
**Protocol lock:** `a83fe33`  
**Oracle selection lock:** `7d36bf8`  
**Final evaluator implementation:** `de8f603`  
**Selection:** `results/post-v1-downstream-multiresolution-selection.json`  
**Evaluation:** `results/post-v1-downstream-multiresolution-evaluation.json`

## Chronology and scope

The frozen downstream gate used 30 PCs, Leiden resolution 0.8, and one seed. After it failed, a
post-hoc sweep showed that its absolute ARIs were often near the within-arm stochastic floor and
that stable lower-resolution partitions existed. This finding motivated, but preceded, the
locked protocol. The analysis below is therefore a locked **post-exploratory robustness
analysis**, not an untouched preregistration. The frozen result remains unchanged and is useful
as a high-resolution clustering-instability reference, not as the primary measure of what
Gravlax retains.

## Oracle-only lock

Commit `a83fe33` fixed the candidate dimensions and resolutions, six deterministic 80% cell
subsamples, optimizer seeds, cluster-size guards, split-half marker requirements, thresholds, and
tie-break before the selector ran. The selector opened no replay matrices. It independently
refit scaling and PCA in every cell subsample and locked:

| dataset | PCs | resolution | clusters | optimizer ARI | subsample ARI | adjacent-dimension | adjacent-resolution | marker Jaccard |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| D0 | 40 | 0.60 | 7 | 0.985 | 0.940 | 0.978 | 0.998 | 0.724 |
| D1 | 30 | 0.40 | 9 | 0.935 | 0.929 | 0.947 | 0.994 | 0.786 |
| D2′ | 40 | 0.25 | 9 | 0.941 | 0.872 | 0.860 | 0.920 | 0.724 |
| D5 mouse 5′ | 20 | 0.35 | 4 | 0.989 | 0.991 | 0.997 | 0.997 | 0.786 |

Every selected partition has 4–15 groups, at least 20 cells per group, no group exceeding 80% of
cells, and split-half marker-effect Pearson correlation above 0.90. Selection took 3:36 wall time
and 1.55 GB peak RSS. The committed result records source digests and all candidate diagnostics.

## Locked replay result

Replay evaluation opened only after the selections and generic evaluator were committed. Twenty
seeds were run in each arm. Cross-arm agreement is judged both by the medoid partitions and
relative to within-arm variability.

| dataset | within ARI oracle/replay | all-pair cross ARI | excess loss | medoid ARI / NMI | 15-NN Jaccard | marker top-25 | control ARI delta | partition |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| D0 | 0.995 / 0.995 | 0.992 | 0.003 | 0.992 / 0.988 | 0.875 | 1.000 | 0.373 | PASS |
| D1 | 0.951 / 0.959 | 0.947 | 0.008 | 0.978 / 0.970 | 0.765 | 0.991 | 0.325 | PASS |
| D2′ | 0.946 / 0.925 | 0.912 | 0.024 | 0.965 / 0.929 | 0.765 | 0.989 | 0.354 | PASS |
| D5 mouse 5′ | 0.944 / 0.944 | 0.946 | -0.001 | 0.997 / 0.990 | 0.875 | 0.970 | 0.326 | PASS |

All four pass the locked partition-retention, marker-retention, and control-sensitivity
components. Replay exceeds the corruption control by 0.186–0.296 in median neighbor Jaccard as
well as by 0.325–0.373 in ensemble cross ARI. The local 15-NN overlap remains 0.765 on D1 and D2′,
so replay does change some local neighbors; those changes do not alter the stable macro-partition.

The composite downstream verdict remains FAIL because it deliberately retains unrelated frozen
pseudobulk thresholds: D1 Spearman is 0.998862 versus the 0.999 bar, and mouse Pearson/Spearman
remain 0.996728/0.997959. These are unchanged count-level findings, not clustering failures.

## Paper interpretation

The supported claim is:

> Gravlax replay retains every nontrivial partition that passed an oracle-only stability and
> marker-coherence lock across the three human datasets and the mouse 5′ stress test. Stable
> macro-partition medoid ARI is 0.965–0.997. Independently recomputed higher-resolution Leiden
> partitions can be unstable even within one arm and should not be used as an absolute fidelity
> oracle.

Do not claim that arbitrary resolution choices, all fine subtypes, or biological annotations are
preserved. No cell-type truth is available in this gate. The resolution-0.8 result belongs in the
reproducibility record and as a brief instability sensitivity observation, not as the headline
measure of Gravlax retention.
