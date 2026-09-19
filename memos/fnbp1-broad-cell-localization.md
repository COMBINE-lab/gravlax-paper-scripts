# D5 FNBP1 broad-cell localization

## Outcome

The fixed FNBP1 cassette event was localized in the existing D5 normal-brain sample with the
unchanged `aie` executable.  The marker rule and all thresholds were frozen before event evidence
was queried by group.  FNBP1 was excluded from both the marker modules and the expression
normalization denominator.

The exact grouped query reproduced the bulk result: 72 inclusion-only, two exclusion-only, zero
both, and 97.3% conservative inclusion among 74 informative molecules.  Event-bearing cells occur
in every broad compartment.  Oligodendrocyte/OPC has the largest count (30) and the highest
descriptive detection rate (9.36 per 1,000 full-H5 nuclei), closely followed by astrocytes (8.88
per 1,000).

| Frozen broad group | Full H5 nuclei | Archive-selected | Event cells / 1,000 full nuclei | Include | Exclude | Both | Usage (Wilson 95%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Neuronal | 923 | 635 | 4.33 | 4 | 0 | 0 | 1.00 (0.51--1.00) |
| Oligodendrocyte/OPC | 3,206 | 2,345 | 9.36 | 30 | 0 | 0 | 1.00 (0.89--1.00) |
| Astrocyte | 788 | 571 | 8.88 | 7 | 0 | 0 | 1.00 (0.65--1.00) |
| Microglia/immune | 1,222 | 386 | 3.27 | 2 | 2 | 0 | 0.50 (0.15--0.85) |
| Vascular | 286 | 167 | 3.50 | 1 | 0 | 0 | 1.00 (0.21--1.00) |
| Unresolved | 3,533 | 2,280 | 7.64 | 27 | 0 | 0 | 1.00 (0.88--1.00) |
| Ambiguous/doublet | 303 | 200 | 3.30 | 1 | 0 | 0 | 1.00 (0.21--1.00) |

Both exclusion-only molecules occur in the broad microglia/immune group, which also contains two
inclusion-only molecules.  This is an interesting lead, not a differential-splicing result: its
usage estimate is based on four molecules and its 95% Wilson interval spans 0.15--0.85.

## Denominators and limitations

The targeted archive dictionary contains 6,584 of the 10,261 called nuclei because nuclei with no
local read are absent.  Rates above deliberately divide event-bearing cells by the complete H5
group size rather than by that selected dictionary.  The marker classifier is intentionally
conservative: 3,533 nuclei are unresolved and 303 are ambiguous/doublet.  These are broad
marker-rule labels rather than an independently validated cell-type annotation.

This is one donor and only 74 informative molecules.  It establishes descriptive localization and
generates a specific donor-resolved follow-up hypothesis; it cannot establish cell-type-specific
splicing, a sample-level contrast, or mechanism.  A suitable next gate would query only this fixed
event in independent normal-brain donors, first testing whether the microglial exclusion signal is
recurrent and then estimating compartment effects across donors.

## Reproduction

- Lock: `experiments/manifests/fnbp1-broad-cell-localization.json`
- Specification: `experiments/specs/fnbp1-broad-cell-localization.md`
- Label builder: `scripts/186_build_fnbp1_broad_cell_groups.py`
- Existing-executable runner: `scripts/187_run_fnbp1_broad_cell_localization.sh`
- Archive-scope filter: `scripts/188_filter_fnbp1_archive_groups.py`
- Compact summarizer: `scripts/189_summarize_fnbp1_broad_cell_localization.py`
- Results: `results/post-v1-fnbp1-broad-cell-localization.{json,tsv}`

The final end-to-end grouped event query took 0.08 seconds and peaked at 135,112 KiB RSS.  The executable
SHA-256 was `e0b7da33af840f0def87729af43d3f0daef5e7b59af025ef85b45a9871c94ffd`.
