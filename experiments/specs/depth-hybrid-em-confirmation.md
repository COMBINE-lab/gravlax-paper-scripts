# Locked untouched confirmation: monotone depth-hybrid EM

**Locked:** 2026-08-31 before fresh D4 GENCODE-v49 quantification, D4 group construction, or any
D4 EM measurement. The implementation is Gravlax
`69a18c2fc227691793907a2a35b3b11a19b8dffc`; its selected development parameters are frozen.
D4 is the primary confirmation and D3 remains held back unless D4 passes.

## Why this gate differs from the D0 development gate

The D0 hybrid improved mean negative log loss by 0.513% over the posterior-mean Dirichlet proxy
and real groups improved the hybrid by 0.654% over shuffled groups. Its registered
shuffle-explained ratio nevertheless failed because that ratio divided the *model-form* gain by a
different, *biological-group-information* gain. That failure remains the formal D0 verdict and is
not relabeled. The untouched test instead predeclares the two scientifically distinct contrasts:

1. **Model form:** hybrid with real groups versus proxy with the same real groups.
2. **Group information:** hybrid with real groups versus hybrid with shuffled groups.

The complete two-by-two hybrid/proxy by real/shuffled factorial is retained for interaction
diagnosis. Fixed convex with real groups is rerun only for the locked evidence-depth guardrail.

## Frozen model and grouping

The transition is `s(d)=d^2/(d^2+64^2)`. It mixes fixed convex
`(lambda_cell,lambda_group,lambda_sample,kappa)=(0.2,0.6,0.2,80)` and the posterior-mean proxy
`(kappa_cell,kappa_group)=(64,80)` from one shared iterative state. No model parameter is fit on D4
or D3.

Fresh STARsolo with GENCODE v49 determines only the called-cell barcode set. Gravlax replays the
archive against v49 for expression values; STARsolo expression values and labels never enter
grouping. Before grouping, exclude the union of genes in all three masked candidate sets. From
the remaining replay matrix, require detection in at least 10 cells, normalize each cell to
10,000 and apply `log1p`, select up to 2,000 replay-derived HVGs, standardize, and fit up to 50
PCs.

Sweep dimensions 20, 30, 40, and 50 and resolutions 0.20 through 1.00 by 0.05. For every point,
run Leiden seeds 0--19 and choose the partition with maximum mean ARI to the other 19 partitions
as its medoid. A point is eligible only if its medoid has 3--12 groups, every group has at least
20 cells, and the median of all 190 pairwise seed ARIs is at least 0.90. Select the eligible point
with maximum median pairwise ARI, then lower dimension, resolution closest to 0.50, and lower
resolution. A deterministic size-preserving permutation supplies the shuffled control. If no
point is eligible, stop without relaxing the rule.

## Decisions

D4 passes only if (i) hybrid-real beats proxy-real on all three mask seeds with at least 0.2% mean
relative NLL improvement and (ii) hybrid-real beats hybrid-shuffled on all three seeds with at
least 0.2% mean relative NLL improvement. Mean Brier must not increase for either contrast;
hybrid top-1 may fall at most 0.10 percentage point versus proxy. Each of the three shallow depth
strata may be at most 0.25% worse than fixed convex, the `16+` stratum may be at most 0.25% worse
than proxy, peak RSS must remain below 4 GiB, and all comparisons are seed-paired.

Only a full D4 pass opens D3. D3 is a directional replication: each primary contrast must have a
positive mean NLL improvement and win at least two of three seeds, with the same Brier, top-1,
depth, and memory guardrails. A replicated pass supports the depth hybrid as the preferred
*masked-recovery evaluator*, not as proof of quantification accuracy; production emission remains
pooled. A full hierarchical Dirichlet implementation is licensed only if real group information
replicates and the hybrid retains structured group-, depth-, or candidate-class-dependent
miscalibration. Otherwise the hierarchy branch closes here.
