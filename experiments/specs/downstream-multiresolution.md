# Locked multiresolution downstream-stability analysis

**Lock date:** 2026-08-31  
**Status at lock:** oracle selection and replay evaluation not yet run  
**Manifest:** `experiments/manifests/downstream-multiresolution.yaml`

## Scope and chronology

The original downstream gate fixed 30 PCs, Leiden resolution 0.8, two optimizer iterations, and
one seed. After that gate failed, post-hoc seed, dimension, and resolution sweeps showed that the
absolute ARI at resolution 0.8 was often below the within-graph stochastic ceiling and that stable
lower-resolution partitions might exist. Those exploratory results informed this protocol.
Consequently this is a **locked post-exploratory robustness analysis**, not an untouched
preregistration. The original result and verdict remain immutable and are retained as a reference
for the instability of independently recomputed high-resolution clustering.

The primary question here is narrower and scientifically better posed: does replay preserve a
nontrivial macro-partition that is demonstrably stable in the oracle data itself?

## Oracle-only selection

The selector must not read replay matrices. For every dataset, retain genes detected in at least
10 oracle cells, normalize each cell to 10,000 counts, apply `log1p`, and choose at most 2,000
Seurat-flavor HVGs from the complete oracle matrix. Evaluate the manifest's fixed grid of PCA
dimensions and Leiden resolutions using 15-neighbor graphs and two Leiden iterations.

For each grid point:

1. Fit the PCA on the complete oracle matrix and run the five locked full-data seeds. The medoid
   partition is the seed partition having greatest mean ARI to the others.
2. Select six deterministic 80% cell subsamples by barcode hash. Refit centering, scaling, and PCA
   independently within every subsample and run the three locked seeds. Compare subsample medoid
   partitions only on cells present in both subsamples.
3. Record optimizer-seed stability, cell-subsample stability, cluster count and size, and agreement
   of the full-data medoid with its best immediately adjacent dimension and resolution.
4. Split cells deterministically into two barcode-hash halves. Holding the full-data oracle medoid
   groups fixed, compute group-versus-rest effects independently in the halves. Record median
   top-25 positive-marker Jaccard and median effect Pearson correlation across eligible groups.

A candidate is eligible only if it meets every numerical threshold in the manifest. The cluster
count, minimum-size, and maximum-fraction guards prevent a trivial coarse solution; split-half
markers require the selected groups to have reproducible molecular distinctions. The exact
selection and tie-breaking rule is frozen in the manifest. A dataset with no eligible candidate
has no stable partition and receives no replay partition verdict.

## Locked replay evaluation

The selector writes only the chosen dimension and resolution plus complete oracle diagnostics.
The evaluator validates the selector and manifest digests, then opens replay matrices. It fits the
oracle scaling and 40-component PCA once, projects replay into that space, slices the locked number
of PCs, independently constructs 15-neighbor graphs, and runs the 20 locked seeds in each arm.

The primary partition statistics are:

- cross-arm ARI between the oracle and replay medoid partitions;
- median ARI over every cross-arm seed pair;
- excess cross-arm loss, defined as the mean of the two within-arm median ARIs minus the median
  cross-arm ARI.

The partition passes only when medoid cross-arm ARI is at least 0.90 and excess loss is at most
0.05. This combines an absolute guard with a comparison to the stochastic ceiling. Absolute
neighbor Jaccard remains descriptive because it measures local graph changes rather than whether
a stable macro-partition changed.

Marker and pseudobulk comparisons retain the prior fixed-oracle-group thresholds. The registered
20% corruption control is reconstructed using the locked oracle clustering; replay must exceed it
by 0.05 in both cross-arm ARI and median neighbor Jaccard. D0, D1, and D2′ must all pass for the
human macro-partition claim. Mouse 5′ is reported separately.

## Interpretation

A pass establishes retention of stable molecular macrostructure relative to fresh STARsolo. It
does not establish biological correctness, discover every subtype, or imply that arbitrary
high-resolution clustering is invariant. Resolution 0.8 remains useful as a documented
counterexample: a downstream algorithm can amplify small count changes in a regime where the
algorithm is itself unstable.
