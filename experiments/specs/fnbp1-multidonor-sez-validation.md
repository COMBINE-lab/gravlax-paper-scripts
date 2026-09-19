# Prospective multi-donor FNBP1 validation in adult human SEZ

## Registration and ordering

The machine-readable gate and event-blind feasibility implementation were committed as
`95b79ad9190063d3a6702f8bdfab7d3cf9cd1cb6` at `2026-09-01T11:01:55-04:00` and pushed before
any GSE234790 raw FASTQ was acquired and before any FNBP1 sequence or event evidence was queried.
The fixed event is
`cassette:chr9:129908999-129923843:129924026-129924959:129908999-129924959`.

The registration already sets the inferential microglia threshold to six eligible donors. This
paragraph makes the mathematical reason explicit without changing the registered rule: for a
two-sided exact sign test, four all-negative signs give `p=0.125` and five give `p=0.0625`; neither
can attain `p<=0.05`. Therefore four or five eligible donors are a **descriptive STOP**, not a FAIL.
PASS/FAIL inference begins only with at least six eligible donors. This clarification and the
event-blind feasibility result were committed before raw acquisition or event inspection.

## Cohort

GSE234790 / PRJNA983239 contains one 10x Chromium Single Cell 3-prime v3.1 experiment from each of
eight independent adult human subependymal-zone donors: four aged 16--22 years and four aged
44--53 years. The public processed matrix contains 36,626 published-QC nuclei. ENA resolves the
study to eight experiments and 40 paired runs totaling 977,312,622 read pairs and 44,692,024,247
compressed FASTQ bytes; every FASTQ has an MD5 identity.

Donor is the unit of replication. Age is metadata only and is not an outcome. The event, groups,
thresholds, and STOP rules cannot be changed after event evidence is opened.

## Event-blind feasibility

Script `190_build_fnbp1_multidonor_feasibility.py` reads only public metadata and expression. It
applies the D5 broad-cell rule independently within each donor. FNBP1 is absent from the marker
modules and normalization denominator. To make non-use auditable, the script recognizes the
single FNBP1 count row only to skip it wholesale: its values are never tokenized, counted, or
retained. The script reads no FASTQ, alignment, archive, junction, or event input.

The feasibility gate passes only if donor identities remain distinct, at least four donors have
at least 1,000 nuclei, every locked marker is available, at least four donors receive at least 50
microglia/immune labels, and all paired FASTQ MD5 identities are complete. Failure means raw
acquisition and event evaluation stop.

## Fixed outcome analysis

Reads will be processed separately by donor. Molecules are deduplicated by donor, called-cell
barcode, and UMI. For every donor and frozen broad group, report inclusion-only, exclusion-only,
both, informative cells, and conservative inclusion usage. Never treat cells or molecules as
replicates.

The brain-wide outcome is evaluable with at least four donors having at least ten informative
molecules. It passes when at least 75% of eligible donors have usage above 0.50 and median donor
usage is at least 0.75.

The microglia localization outcome requires at least six donors having both at least four
microglia/immune informative molecules and at least ten informative molecules among other resolved
groups. Below six, report exact counts only. If evaluable, the estimand is each donor's
microglia/immune usage minus usage among all other resolved groups. It passes when at least 75% of
signs are negative and the two-sided exact sign-test p-value is at most 0.05.

No result establishes mechanism, an age effect, generalization beyond SEZ, or cell-level
independence.
