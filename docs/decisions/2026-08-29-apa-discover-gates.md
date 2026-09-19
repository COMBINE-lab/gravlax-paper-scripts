# Frozen gates — APA and DISCOVER capability demonstrations

**Date:** 2026-08-29 · **Status:** FROZEN before implementation, per standing practice.

Both capabilities are demonstrated against **labeled truth we already possess**: the withheld
panel `annotations/withheld/w1.manifest.tsv` (900 real expressed features held out of v49), whose
oracle counts are known. No synthetic data, no new sequencing.

## APA — 3′-end usage from the archive

The archive stores each chain's span extremes; a molecule's 3′-most coordinate (strand-aware:
max end on +, min start on −) is its polyadenylation-side evidence. The count matrix cannot
represent this at all — stated as the negative control, not tested.

| Gate | Test | PASS | STOP |
|---|---|---|---|
| **APA-1 fidelity** | 3′-site assignments computed from the archive ≡ computed from the BAM rows (same `flatten` bridge) | byte-identical site tables | any mismatch |
| **APA-2 labeled recovery** | For the **300 utr3-withheld genes** (terminal exons trimmed 500 bp in w1): the apa query on each gene's v49 span must place site mass inside the trimmed window — 3′ usage the w1 annotation cannot see | ≥90% of the 300 genes show ≥1 site with ≥3 UMIs inside their trimmed window; aggregate in-window UMI mass within 30% of the 19,662 UMIs the w1 oracle measurably lost | <70% of genes, or mass off by >2× |
| **APA-3 sanity** | Genome-wide: fraction of 3′-site UMI mass within 200 bp of an annotated v49 transcript 3′ end | ≥60% (10x 3′ chemistry is end-biased but has internal priming) | <30% — sites would not be 3′ ends at all |

## DISCOVER — unannotated transcription

A molecule unclaimed by the supplied annotation (no overlap with any transcript span, loose and
therefore conservative about novelty) is candidate evidence; single-linkage clusters of unclaimed
molecules are candidate loci.

| Gate | Test | PASS | STOP |
|---|---|---|---|
| **DISC-1 labeled recall** | `discover --gtf w1`: the 300 withheld whole genes are real expressed loci absent from w1 | ≥95% of withheld genes with oracle ≥50 UMIs recovered as a candidate locus overlapping the true gene span | <80% |
| **DISC-2 quantification** | Candidate-locus UMI counts vs the v49 oracle counts of the matched withheld genes | median relative error ≤15% (locus boundaries are data-driven, not annotation-driven, so exact equality is not expected) | >50% |
| **DISC-3 precision structure** | Classify all candidates against v49: (a) withheld-gene hits, (b) hits on other v49 genes absent/changed in w1, (c) overlapping nothing in v49 (potentially genuinely novel) | reported, with (a)+(b) as the labeled-precision proxy | — |

## Parameters (frozen)

APA site clustering gap: 24 bp. Discover: molecule clustering gap 1,000 bp, min 10 UMIs per
candidate (class-deduplicated), claim test = overlap with any transcript span of the supplied GTF,
same strand OR either strand for claiming (claim generously, discover conservatively: a molecule is
claimed if it overlaps a span on EITHER strand — antisense transcription is real but reporting it
as "novel" against the sense gene would flatter recall).

All full-scale; dev-scale results are for debugging only (standing rule, five instances now).
