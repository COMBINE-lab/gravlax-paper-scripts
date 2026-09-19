# Prospective SEZ transcript-end atlas v1

## Question and protocol fit

The primary biological question is whether the transition from published
Astro/NSC populations to mature neuronal populations in the adult human
subependymal zone is accompanied by recurrent distal RNA 3′-end usage. This is
the protocol-native alternative to a genome-wide internal-splicing screen in
10x 3′ v3.1 data. The analysis measures alignment-supported terminal usage; it
does not call an alignment endpoint a cleavage site without additional
sequence or catalogue support.

GSE234790 contains eight independent donors. Published clusters 0, 7, and 16
define `astro_nsc`; clusters 2, 3, 10, 12, and 19 define `mature_neuron`. The
mapping is frozen from the source publication before inspecting any
whole-genome endpoint result. Donors A--D are youth and E--H are middle-aged.
Age is secondary and descriptive because four-versus-four replication is too
small for a new genome-wide interaction screen.

All 36,626 published-QC barcodes remain in the archive whitelist. Published
clusters outside the named analysis groups are labeled `other`; they are not
dropped merely because they are outside the primary contrast.

## Common terminal-site catalogue

One annotation-free archive represents each donor. GENCODE v49 supplies gene
spans, strand, terminal-exon structure, and labels, but never supplies the
observed site coordinates. Molecules are assigned by transcript-oriented
aligned 3′ base to same-strand gene windows extended by 2 kb. Ambiguous
overlapping windows are reported and excluded from the primary analysis rather
than multiply assigned.

Exact endpoint coordinates are pooled only for catalogue construction and are
clustered once across all donors with a 24-bp maximum adjacent gap. Every
sample × group × site count is then emitted, including explicit zeros. A site
enters the recurrent catalogue with at least 10 UMIs in at least three donors.
The primary high-confidence catalogue additionally requires no internal-
priming flag and either a same-strand PolyASite site within 50 bp or a canonical
AAUAAA-family signal 5--60 bp upstream plus recurrence in at least four donors.
Unmatched sites remain in descriptive output.

`endpoint` means the transcript-oriented terminal aligned base. `site` means a
cluster of those endpoints. `catalogue-supported site` and `tail-supported
site` are distinct labels. `cleavage site` is forbidden unless tail sequence or
an orthogonal 3′-end assay supports it.

## Primary estimand and inference

For every eligible gene, sites are ordered in transcript direction. The distal
usage index for a sample/group is the UMI-weighted site rank scaled to [0,1]. A
donor is eligible for that gene when both primary groups have at least 20
high-confidence site UMIs. A gene is tested with at least six eligible donors,
at least two high-confidence sites, and at least 20 distal-site UMIs across the
eligible rows.

The per-gene effect is the mean paired donor difference
`mature_neuron - astro_nsc`. The confirmatory test is a two-sided paired
Student t test over donor-level distal usage indices; an exact two-sided
sign-flip p-value over the same paired differences is emitted as calibration,
not substituted after results are seen. BH correction covers all tested genes.
A reported gene-level discovery requires q ≤ 0.05, absolute mean effect ≥ 0.10,
the same sign in at least six eligible donors, and no leave-one-donor-out mean
sign reversal.

The preregistered global directional hypothesis is that each donor's median
eligible-gene shift is positive. A protocol-level directional PASS requires at
least seven of eight donor medians above zero and an across-donor median of at
least 0.02. This one-sided direction was fixed from prior neural-differentiation
APA literature; it is reported separately from per-gene two-sided tests.

Age effects are descriptive. They may rank terminal programs for future
validation but cannot create a paper claim, enter the primary correction
family, or rescue a failed lineage analysis.

## Accuracy and negative controls

- Synthetic archives must reproduce exact site counts, explicit zeros,
  transcript-direction ordering, internal-priming flags, motif orientation,
  and paired-test effect signs.
- On a frozen real-data subset, the optimized reducer must be byte-identical to
  an independent slow reference reducer.
- Primary sites must have a higher PolyASite match rate than internal-priming
  flagged sites; otherwise site accuracy is `STOP` and no biological test is
  promoted.
- A within-donor shuffled group map must yield no more than 5% q ≤ 0.05 genes
  and must not pass the global directional rule.
- Results are repeated at site gaps 12, 24, and 48 bp. A reported gene cannot
  reverse effect sign at both neighboring gaps.
- A sparse raw-read audit will measure non-templated poly(A) soft clips before
  any archive-format change. Tail evidence is an optional accuracy tier, not a
  prerequisite invented after observing the biological result.

Alternative-last-exon paths are a secondary capability gate. A path couples a
same-molecule terminal junction chain to a terminal-site cluster. It is
promoted only if at least 100 genes have two recurrent paths supported in at
least three donors; otherwise the result remains APA-only without penalty.

## Performance and failure interpretation

The engine scans each archive once, aggregates exact endpoint coordinates
before site construction, and never stores a cell × site matrix. The fixed
host target for eight donors is ≤90 s wall time and ≤8 GiB peak RSS, with at
least a 2× wall-time or RSS improvement over eight independent legacy
`apa-test` reductions. A hard union-site limit fails without truncation.

Biological underpowering is not a software failure. A failure of PolyASite/IP
separation is an accuracy failure. Count disagreement, strand errors, unstable
site construction, pseudoreplication, missing zero rows, or silent truncation
is a software failure. Failure of the neural directional hypothesis is a
valid negative biological result and cannot be repaired by selecting another
cell contrast after inspection.
