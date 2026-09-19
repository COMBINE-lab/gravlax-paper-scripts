# Downstream clustering, marker, and velocity fidelity gate

Status: complete 2026-08-31. Human Gene partition gate FAILS on D1/D2′ while marker and
pseudobulk summaries remain highly concordant; the velocity consequence gate PASSES all three
human datasets and the mouse stress test. Thresholds below are unchanged from preregistration
commit `d3e9b3a`.

## Outcome

D0 passes every Gene threshold. D1 and D2′ fail the independent Leiden/neighbor portion (ARI
0.834/0.619; median 15-NN Jaccard 0.765/0.765), so the preregistered human Gene verdict is FAIL.
Their mean top-25 marker overlaps remain 0.987/0.962, median marker-effect correlations exceed
0.9988, and pseudobulk Pearson correlations exceed 0.99996. D5 similarly fails partition and
pseudobulk thresholds despite 0.875 median neighbor Jaccard and perfect mean top-25 marker overlap.
The structured positive control is detected in every confirmatory dataset.

The separately frozen steady-state velocity consequence probe passes all four datasets: median
cell-vector cosine is 0.993–0.999 and median top-3 transition Jaccard is 1.0 (mean 0.873–0.914).
Post-hoc failure localization, which does not change the verdict, finds split/merge instability
inside coherent myeloid, neuronal, and B-cell compartments; fixed oracle-centroid label agreement
is 95.1–99.8%, and mismatched cells do not have materially larger count errors. See
`results/post-v1-downstream-fidelity.json` and
`results/post-v1-downstream-cluster-diagnostics.json`.

## Question and scope

The archive/fresh-pipeline comparisons move 0.22–0.75% of Gene UMI mass and larger fractions in
some Velocyto components. This gate asks whether those differences change common biological
summaries. It does not use archive-versus-BAM identity—which is already exact—as a proxy for this
question.

The confirmatory datasets are the frozen roadmap's human D0, D1, and D2′ experiments under
GENCODE v49. The public mouse 5′ v2 D5 experiment is a cross-species/chemistry stress test. All
comparisons use the fresh STARsolo Gene-filtered barcode set, identical gene order, fixed seeds,
and no cell-type labels. Consequently this establishes stability of downstream conclusions, not
biological correctness of an annotation or clustering.

## Fixed Gene pipeline

For each dataset:

1. Subset fresh STARsolo and Gravlax raw matrices to the oracle's called cells and common genes.
2. Retain genes detected in at least 10 oracle cells. Normalize every cell to 10,000 counts and
   apply `log1p`.
3. Select up to 2,000 highly variable genes from the oracle only using Scanpy's Seurat flavor.
   Learn centering, nonzero-variance scaling, and 30 principal components from the oracle only;
   project replay and the positive control into this same coordinate system.
4. Build independent 15-nearest-neighbor graphs in the shared PCA space. Run Leiden at resolution
   0.8 and seed 1729 independently on the oracle, replay, and control graphs.
5. For marker comparisons, hold the oracle Leiden groups fixed in all arms. Rank by absolute
   group-versus-rest log-normalized mean difference and compare the top 25 genes per group. This
   separates marker stability from cluster-label instability.

Report adjusted Rand index (ARI), normalized mutual information (NMI), per-cell neighbor Jaccard,
top-25 marker Jaccard, marker-effect Pearson correlation, and Pearson/Spearman pseudobulk
correlation. Cluster labels are compared without manual relabeling.

## Sensitivity control

After oracle clustering, select exactly 20% of cells by the lowest SHA-256 hash of
`dataset|barcode`. Replace each selected expression profile with the profile of a deterministic
donor from the next nonempty oracle cluster, cycling donors with replacement. Process this arm
identically. The control is intentionally not a biological alternative; it tests whether the
chosen stability metrics detect structured cell-profile corruption.

The gate is sensitive only if replay exceeds the control by at least 0.05 in ARI and median
neighbor Jaccard on every confirmatory dataset. A sensitivity failure prevents a downstream
stability claim even if replay appears concordant.

## Gene gates

A dataset passes Gene downstream fidelity only if all of the following hold:

- ARI ≥ 0.95 and NMI ≥ 0.95;
- median per-cell neighbor Jaccard ≥ 0.80;
- mean top-25 marker Jaccard across oracle groups ≥ 0.80;
- median per-group marker-effect Pearson correlation ≥ 0.98;
- pseudobulk Pearson and Spearman correlations ≥ 0.999;
- the sensitivity-control deltas above pass.

D0/D1/D2′ must all pass for the human downstream-stability claim. D5 receives the same thresholds
but is reported separately: a pass strengthens the mouse 5′ Gene claim, while failure bounds that
claim without invalidating the human result.

## Fixed velocity comparison

Velocity is evaluated separately and does not inherit the Gene verdict. On the same cells and
oracle HVGs, retain genes detected with spliced counts in at least 20 cells and unspliced counts in
at least 10. Normalize spliced and unspliced layers jointly to 10,000 molecules per cell. Using the
fixed oracle Gene neighbor graph, average each cell with its 15 neighbors. Fit one steady-state
coefficient per gene from oracle totals, `gamma_g = sum(U_g) / sum(S_g)`, and construct residual
vectors `U - gamma*S` for both arms. Project residuals through the fixed oracle Gene scaling and
PCA loadings.

Report cell-wise cosine agreement and, within the same fixed oracle-neighbor candidates, the
Jaccard overlap of the three highest positive-cosine destinations. A dataset passes velocity
consequence fidelity if median cosine ≥ 0.90, at least 95% of cells have positive cosine, and
median top-3 transition Jaccard ≥ 0.80. Results are reported per dataset; a miss is MARGINAL when
median cosine remains ≥0.75 and STOP otherwise. This is a transparent steady-state consequence
probe, not a claim that a complete dynamical scVelo model is exact.

## Failure interpretation

- Gene failure with source-exact replay means a small fresh-alignment deviation is biologically
  structured; identify the responsible genes/cells and narrow the paper claim.
- Velocity failure does not negate Gene replay. It defines the boundary of the current evidence
  reduction and alignment protocol.
- A control-sensitivity failure means the metric panel is too weak and must not be rescued by
  changing thresholds after seeing outcomes.
- Results must retain exact package versions, seeds, input hashes, per-dataset metrics, and the
  complete positive-control outcomes in compact JSON.
