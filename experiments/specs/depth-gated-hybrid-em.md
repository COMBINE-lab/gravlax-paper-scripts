# Locked D0 gate: monotone evidence-depth hybrid

**Locked:** 2026-08-31 after implementation and tests, before any D0 hybrid run. Gravlax commit
`9ccc3f65bd665dbbc30dc84de36bb2254d2cdbe7` is the implementation under test. D0 selected both
constituents and exposed their depth tradeoff, so this remains development rather than independent
confirmation. D4 is not touched by this gate.

## Model

Freeze the previously selected fixed-convex parameters at
`(lambda_cell,lambda_group,lambda_sample,kappa)=(0.2,0.6,0.2,80)` and the selected posterior-mean
proxy concentrations at `(kappa_cell,kappa_group)=(64,80)`. Let `d` be the initial fitted unique
mass over the target candidate set, calculated before any evaluated mode changes its state. For a
single transition scale `D>0`, use

```text
s(d;D) = d^2 / (d^2 + D^2)
q_hybrid = (1-s) q_fixed-convex + s q_Dirichlet-proxy.
```

The gate is monotone, equals one-half at `d=D`, and is candidate-independent within a target, so
the hybrid remains normalized. Both component probabilities are recomputed from the hybrid's one
shared current EM state; it is not an after-the-fact mixture of separately fitted outputs. The
square is frozen rather than selected because it yields a smooth but useful separation around the
observed 16-unit boundary without adding a second tuning parameter.

## Development grid and controls

On D0 only, cross `D={2,4,8,16,32,64,128,256,512,1024}` with mask seeds 7, 17, and 29. Select
minimum mean negative log loss, breaking an exact tie toward the larger scale (closer to the
already safer fixed-convex arm). Rerun fixed convex and the selected Dirichlet proxy on real and
size-preserving shuffled groups with the same binary, then run the selected hybrid on shuffled
groups. Preserve the existing candidate-gene exclusion, called-cell population, group map, masks,
and four evidence-depth strata.

The earlier depth audit implies a hard piecewise upper reference: fixed convex in `0-1`, `(1,4]`,
and `(4,16]`, and the proxy in `16+`. This is explicitly a post-hoc oracle reference, not a model,
selection arm, or confirmation claim. Its sole purpose is to quantify how much of the known
constituent tradeoff a one-parameter monotone transition recovers.

## Decision

Pass requires every item below:

1. Mean NLL improves at least 0.2% over the selected proxy, the better aggregate constituent.
2. The hybrid recovers at least half of the hard-stratum oracle's NLL gain over the proxy.
3. Each shallow stratum is no more than 0.25% worse than fixed convex, and `16+` is no more than
   0.25% worse than the proxy.
4. Real groups beat shuffled groups on all three seeds. The shuffled-label gain over the proxy
   explains no more than half of the real-group hybrid gain over the proxy.
5. Mean multiclass Brier does not increase versus the proxy, and top-1 falls by no more than 0.10
   percentage point versus the better constituent.

All comparisons use matching seeds and the newly rerun constituents. The shuffle-explained
fraction is `(NLL_proxy,shuffled - NLL_hybrid,shuffled) /
(NLL_proxy,real - NLL_hybrid,real)` and is not clipped. PASS earns a separately locked untouched
D4 confirmation; it does not license a full variational Dirichlet model. A positive aggregate
gain that misses a depth or guardrail criterion is MARGINAL. No aggregate gain or failure of the
real-over-shuffled control is STOP. Production remains pooled.
