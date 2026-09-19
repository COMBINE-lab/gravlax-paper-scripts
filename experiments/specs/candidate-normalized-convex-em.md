# Locked gate: candidate-normalized convex partial pooling

**Locked:** 2026-08-31 after implementation, but before a parameter grid or any D4/D3 EM
measurement. Gravlax commit `f6fd759ae118049ae0a11ef9f50b0d32d5daee07` is the implementation
under test. A single release-mode D0 smoke run at the CLI defaults was observed before this lock;
its exact values are disclosed in the manifest. D0 is development-only, so that observation is
not confirmation evidence.

## Why this model

The registered additive hierarchy failed because its pseudo-counts were normalized over all genes
before being restricted to a target candidate set. A pseudo-count of 80 consequently had median
effective mass `3.60e-5` in the tested classes. Increasing the mass repaired much of the
calibration, but a tuned no-group arm explained most of that improvement and D1/D2-prime had
already been inspected. This gate tests the simpler model suggested by that diagnosis.

For target class `k` in cell `c`, let `C_k` be its candidate genes and let `pi_cg`, `Pi_hg`, and
`Pi_g` be the current cell, group, and sample abundances. Define

```text
p0(g | C_k) = (Pi_g + epsilon) / sum[j in C_k](Pi_j + epsilon)
pc(g | C_k) = pi_cg / sum[j in C_k] pi_cj
n_hg,-c     = max(Pi_hg - pi_cg, 0)
ph(g | C_k) = (n_hg,-c + kappa p0(g | C_k))
              / (sum[j in C_k] n_hj,-c + kappa)
q(g | c,C_k) = lambda_c pc + lambda_h ph
             + (1 - lambda_c - lambda_h) p0.
```

Empty cell or leave-one-cell-out group components fall back to `p0`. Cells absent from the group
map transfer `lambda_h` to the sample component. All terms are distributions on the same
candidate set and the weights lie on a simplex. At `(lambda_c,lambda_h)=(0,0)`, the model equals
pooled candidate probabilities. The packed implementation makes two passes over each class and
stores no per-target probability vectors. `--convex-only` prevents a grid from recomputing the six
reference arms. Production `--emit` remains pooled.

## Development and leakage control

D0 alone selects parameters. The grid crosses cell weights `{0,.05,.10,.20}` with group weights
equal to `{0,.25,.50,.75,1}` of the remaining simplex mass. For a positive group weight,
`kappa` is in `{0,5,20,80}`; group-zero configurations run once at `kappa=20`. Select the lowest
mean negative log loss across mask seeds 7, 17, and 29, breaking exact ties by lower group weight,
then lower cell weight, then lower `kappa`. The comparator is the best group-zero configuration
from this same grid, not untuned pooled EM.

Groups retain the earlier leakage controls: construct them only from the pre-EM Gravlax Gene
matrix, exclude the union of genes occurring in masked candidate sets, use the fixed called-cell
population, and select the medoid over 20 Leiden seeds. A deterministic size-preserving label
permutation is the group-information control. Group construction and any dimension/resolution
choice must be frozen without examining convex-EM outcomes.

D1 and D2-prime are descriptive only because both were inspected during the additive scale audit.
D4 (500 PBMC 3' LT Chromium X) is the untouched primary confirmation. Before unblinding it,
generate a fresh GENCODE-v49 STARsolo quantification to fix only its called-cell barcodes, replay
its archive against v49, and freeze a stable macro-partition without inspecting EM scores. Sweep
dimensions `{20,30,40,50}` and resolutions 0.20 through 1.00 by 0.05 over 20 Leiden seeds. A
configuration is eligible only with 3--12 clusters, at least 20 cells per cluster, and median
pairwise seed ARI at least 0.90. Select maximum median ARI, breaking ties by lower dimension,
resolution closest to 0.50, then lower resolution; use that configuration's ARI-medoid labels. If
none is eligible, D4 is STOP rather than relaxing the rule. D3 (PBMC 10k v3) remains untouched as
the larger replication and uses the same rule only if D4 passes. Fresh STARsolo expression values
or cluster labels never enter the EM model.

## Metrics and decisions

Preserve per-seed top-1, truth probability, negative log loss, multiclass Brier, calibration bins,
wall time, peak RSS, group/candidate digests, binary digest, code commit, command lines, and raw
logs. Masked recovery is a controlled evaluation of ambiguous assignment, not ground-truth
quantification accuracy.

D4 passes only if real groups beat the tuned no-group arm and shuffled labels on every seed, the
mean relative NLL improvement over no-group is at least 0.5%, shuffled labels explain no more than
half of the group increment, mean Brier does not increase, top-1 falls by at most 0.20 percentage
point, and peak RSS remains below 4 GiB. A pass triggers D3 replication; production adoption still
requires D3 to preserve the direction and real labels to beat shuffled labels on at least two of
three seeds. A repeatable proper-score gain that misses a magnitude or top-1 guardrail is marginal.
No real-over-shuffle/no-group increment is STOP.

## Reserved hierarchical Dirichlet escalation

The full model remains a designed fallback rather than the default next step:

```text
theta_0 ~ Dirichlet(eta)
theta_h ~ Dirichlet(kappa_0 theta_0)
theta_c ~ Dirichlet(kappa_h theta_h)
P(class k | theta_c) = sum[g in C_k] theta_cg.
```

Concentrations would be estimated by empirical Bayes or variational inference over packed
equivalence classes. Escalate only if real groups reproducibly beat shuffled and no-group controls
but a fixed mixture leaves group-depth-, candidate-class-, or uncertainty-dependent calibration
errors. It then requires new untouched data or nested cross-fitting. Do not escalate merely because
the generative model is more formal; if D4 shows no real group signal, stop the hierarchy ladder.
