# Locked D0 screen: posterior-mean Dirichlet proxy

**Locked:** 2026-08-31 after implementation and a disclosed default-parameter smoke run, but
before the parameter grid. Gravlax commit
`16c012db53cd4fa4dad30dc2d6f8321de2c73f9f` is the implementation under test. This is a quick,
development-only D0 gate; it is not confirmation and does not unblind D4.

## Question and model

This gate asks whether evidence-adaptive shrinkage has enough incremental value to justify a full
hierarchical Dirichlet implementation. It evaluates a packed posterior-mean proxy, not the full
variational model. For a target candidate set `C`, let `p0` be the candidate-normalized sample
distribution, `n_h,-c` the leave-one-cell-out group mass, and `n_c` the target-cell mass. Then

```text
p_h(g | C) = [n_h,-c(g) + kappa_0 p0(g | C)]
             / [N_h,-c(C) + kappa_0]
p_c(g | C) = [n_c(g) + kappa_h p_h(g | C)]
             / [N_c(C) + kappa_h].
```

The model scores with `p_c`. It therefore adapts continuously: high-depth cells rely more on
their own fitted evidence, low-depth cells shrink toward their group, and low-depth groups shrink
toward the sample. Empty levels fall back to their parent. The target cell is removed from its
group contribution. Two packed passes per class avoid per-target probability vectors.

## Locked selection and controls

D0 crosses `kappa_h={1,4,16,64,256}` and `kappa_0={1,5,20,80,320}` for mask seeds 7, 17, and 29.
Select minimum mean negative log loss, breaking exact ties by lower `kappa_h`, then lower
`kappa_0`. Rerun the already selected fixed convex model
`(cell,group,sample)=(0.2,0.6,0.2), kappa=80` on real and deterministically shuffled groups with
the same binary so both comparators contain the newly added evidence-depth strata. Run the
selected proxy with shuffled groups. Group construction, masking, and leakage controls are
unchanged from the locked convex gate.

Evidence depth is the initial fitted unique-candidate mass available in the target cell, before
the evaluated EM mode. Report `0--1`, `(1,4]`, `(4,16]`, and `>16` strata. This measure is
mode-independent, so paired depth comparisons are valid.

## Decision

Escalation requires every criterion below:

1. Mean NLL improves at least 0.25% over fixed convex.
2. The selected real-group proxy beats its shuffled-group counterpart on every seed.
3. The shuffled-label proxy explains no more than half of the proxy's absolute NLL gain over its
   matched fixed-convex control.
4. Mean Brier does not increase and top-1 drops by no more than 0.10 percentage point.
5. Among strata with at least 1,000 observations, at least two improve in NLL, at least one of the
   first three strata improves by 0.25%, and none regresses by more than 0.5%.

Here the shuffle-explained fraction is
`(NLL_fixed,shuffled - NLL_proxy,shuffled) /
 (NLL_fixed,real - NLL_proxy,real)`. It is intentionally not clipped; a nonpositive real gain is
an automatic stop.

Passing this D0 screen only earns a separately locked implementation and validation gate for the
full hierarchical Dirichlet model. A positive real-over-shuffle signal that misses a magnitude,
depth, or guardrail criterion is marginal and does not justify implementation yet. No aggregate
gain, no consistent real-over-shuffle signal, or a shuffle-explained fraction above 0.5 is STOP.
Production remains pooled in all cases.
