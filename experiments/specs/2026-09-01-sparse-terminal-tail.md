# Prospective pilot: sparse terminal poly(A)-tail evidence

Status: **frozen before looking at cohort enrichment results**
Base commit: `e378516a0fb6602e63d626d0c49fb11e58f463c8`
Protocol: 10x Chromium 3' gene expression, `SoloStrand=Forward`

## Question

Can a deliberately sparse record of untemplated terminal-tail evidence strengthen cleavage-site
queries without turning an annotation-independent archive into a sequence archive?

This pilot may inspect the input BAM sequence only at the transcript-3'-terminal soft clip. It may
retain the cleavage anchor, strand, corrected cell barcode, UMI, clipped length, tail-base count,
and terminal homopolymer run. It must not retain bases, qualities, aligned sequence, or general
soft clips. The default `.aie` format and ingest path do not change.

## Frozen evidence semantics

SAM/BAM `SEQ` is alignment-oriented. Under `SoloStrand=Forward`, an alignment's strand is the
transcript strand. A forward alignment therefore contributes only its trailing soft clip and is
tested for A; a reverse alignment contributes only its leading soft clip and is tested for T. The
cleavage anchor is the 0-based exclusive aligned end on `+` and the 0-based inclusive aligned start
on `-`, matching Gravlax's transcript-end convention.

A read is tail-positive exactly when all of the following hold:

- the relevant terminal operation is a soft clip of at least 6 nt (hard clips do not count);
- at least 80% of clipped bases are the oriented tail base; and
- the transcript-distal homopolymer run is at least 4 nt.

Only primary mapped records with an ACGT UMI and a barcode that is exact or uniquely one mismatch
from the supplied whitelist are eligible. Evidence is deduplicated by
`(reference, anchor, strand, corrected barcode, UMI)`, retaining the strongest signal for a key.

The optional sidecar is a versioned, zstd-compressed stream of fixed-width 16-byte records. Each
record contains a 16-bit reference id, 32-bit anchor, packed 32-bit barcode, packed 32-bit UMI, and
a 16-bit saturated signal word (strand, clip length, tail count, terminal run). Reference names and
the frozen thresholds are in the header. This is an experimental rebuildable sidecar, not an E1
archive section.

## Frozen evaluation

Data are the eight GSE234790 adult human SEZ donor BAMs already used for the PolyASite-mixture
analysis. PolyASite 2.0 is the external catalogue; GRCh38 is used only to flag genomic
internal-priming candidates by the existing Gravlax rule (at least 12 transcript-oriented A bases
in the next 20 nt or an A run of at least 8 in the next 140 nt).

For every terminally soft-clipped read meeting the 6-nt length floor, compare an external site
(anchor within 25 bp) with its transcript-downstream 500-bp shifted control. A control is excluded
when it lies within 100 bp of any external site. Candidate discrimination is evaluated on events
exclusive to one class using (i) clip length and (ii) the continuous oriented-tail score, ordered
lexicographically by tail fraction then terminal run. Report AUC, the tail-positive odds ratio,
and exact counts. Candidate summaries are molecule-deduplicated.

NTRK2 is a secondary interpretation check, not a promotion gate: count tail-confirmed molecules
within 25 bp of the nine fitted NTRK2 catalogue sites and compare the proximal shorter-coding and
distal MANE/full-length-like terminal groups using the already frozen cell groups. A terminal-tail
witness does not phase the terminal exon to the coding exons and must not be called TrkB.T1 or
TrkB.FL molecule-level evidence.

## Prospective gates

The pilot is promoted only if all primary gates pass:

1. orientation and packing round-trip tests pass, and reverse-complemented synthetic records give
   the same transcript-oriented call;
2. tail-positive molecules are at least 5-fold enriched at external sites versus shifted controls,
   with at least 100 external-site witnesses across the cohort;
3. the continuous tail score improves external-site discrimination by at least 0.05 AUC over clip
   length alone;
4. tail-positive support is depleted, not enriched, at genome-flagged internal-priming candidates;
5. the compressed sidecar is at most 1% of the corresponding `.aie` bytes and no retained field
   contains nucleotide sequence or quality.

If any primary gate fails, the result remains an ingest audit and the archive format is unchanged.
NTRK2 can strengthen the story only if it has at least 20 witnesses and its cell-group direction is
consistent in at least six donors; failure here does not override the primary decision.
