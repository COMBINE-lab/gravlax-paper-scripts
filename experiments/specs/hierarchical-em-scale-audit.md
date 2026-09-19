# Post-hoc hierarchical-EM scale audit

**Status:** exploratory diagnosis, not a new confirmatory gate. The registered experiment and its
STOP verdict remain in `experiments/specs/hierarchical-em.md` and
`results/post-v1-hierarchical-em.json`.

## Why this audit was necessary

The registered additive model uses whole-transcriptome frequencies:

```text
w(c,g) = local(c,g)
       + alpha_group  * group(h(c),g) / total_group(h(c))
       + alpha_global * global(g)     / total_global
       + epsilon.
```

Its D0 optimum occurred at the largest registered values, 80/20. Those numbers are total prior
masses over all genes, not the effective prior mass within one target's candidate set. The audit
therefore measures the actual prior contribution beside a unit local count and extends the scale
geometrically. It also includes a group-free `0/25,600` arm so that calibration from stronger
cell/global shrinkage is not misattributed to biological groups.

## Frozen reconstruction

- Runner: `scripts/149_hierarchical_em_scale_audit.sh`
- Summarizer: `scripts/150_summarize_hierarchical_em_scale_audit.py`
- Compact result: `results/post-v1-hierarchical-em-scale-audit.json`
- Raw outputs: `runs/post-v1/hierarchical-em-scale-audit-r1` (excluded from Git)
- Seeds: 7, 17, 29
- D0 grid: total masses 6,400, 25,600, and 102,400 with group fractions 0, 0.2, 0.5, and 0.8,
  plus an 80%-group point at total mass 409,600.
- Exploratory D0 selection: minimum mean negative log loss, with lower group then global mass as
  deterministic tie-breakers.
- Comparators: pooled, calibrated cell/global (`0/25,600`), real hierarchy, and size-preserving
  shuffled hierarchy.

The candidate-gene exclusion, replay-only group labels, called-cell denominator, packed target
representation, ten iterations, and one-group implementation invariant are inherited unchanged
from the registered experiment.

## Inferential limit

D1 and D2-prime were already inspected when the scale audit was designed. Their enlarged-scale
results can diagnose whether under-scaling explains the registered failure, but cannot rescue the
original gate or serve as untouched confirmation. A positive hierarchy claim requires a newly
registered, candidate-normalized model and genuinely untouched data or independently generated
ambiguity truth.
