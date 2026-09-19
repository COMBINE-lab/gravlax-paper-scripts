# Prospective promotion: sparse terminal-tail archive evidence

Status: **frozen before production-layout measurements**
Pilot base: `db817cf967c4c6d8df442bb1a73efa81e1152909`
Protocol: 10x Chromium 3' gene expression, `SoloStrand=Forward`

## Promotion question

Can the successful `AIETAIL1` pilot be represented inside an optional `.aie` extension by attaching
each retained terminal-tail witness to the exact molecule record that ingest emits, without
retaining cell barcode, UMI, nucleotide sequence, quality, or general soft-clip sequence?

The proposed compact layout is a per-chunk sparse set of local molecule ordinals followed by one
16-bit oriented tail signal per selected ordinal. It is enabled only by an explicit ingest flag.
With the flag absent, ingest must execute the existing extraction and writing path and produce the
same archive bytes as the frozen binary.

## Frozen semantic invariant

The read-level evidence definition, orientation, thresholds, cleavage coordinate, and strongest-
signal ordering are exactly those in `2026-09-01-sparse-terminal-tail.md`. Promotion must not change
the evidence unit merely to fit the representation.

A retained pilot witness has key
`(reference, cleavage anchor, strand, corrected barcode, UMI)`. An ordinal attachment is exact only
if all of the following are true:

1. the witness maps to one and only one emitted molecule record;
2. its cleavage anchor is recoverable unambiguously from that record's retained placement geometry;
3. every retained witness for the record has the same recoverable anchor and the same strongest
   16-bit signal; and
4. archive-only candidate matching reproduces the pilot witness key, rather than substituting the
   molecule's genomic sort anchor or a different representative read's endpoint.

The 16-bit signal remains `(strand, saturated clip length, saturated oriented-tail count,
saturated terminal run)`. It may not encode a candidate, gene, annotation, cell, UMI, sequence, or
quality. Candidate catalogues remain query-time inputs.

## Structural hard stop

Before implementing the section, construct forward- and reverse-strand BAM fixtures and audit the
real donor-A extraction. If a tail-positive read can be absorbed into an emitted molecule while
its cleavage endpoint is not retained, or if one emitted molecule can carry distinct retained
cleavage anchors, the proposed ordinal-plus-16-bit layout is **not lossless**. Stop promotion and
report the minimal counterexample. Do not choose a representative endpoint, union candidate hits,
or keep only the strongest anchor.

In that event, specify the smallest lossless alternative. The anticipated fallback is a sparse
ordinal set plus a per-selected-molecule event list of signed cleavage-anchor deltas and 16-bit
signals. That is a separate layout decision and requires approval because it changes the proposed
representation and permits more than one event per molecule.

## Target optional layout, conditional on exactness

If and only if the structural gate passes:

- archive meta declares `terminal_tail_layout = "chunk-ordinal-signal-v1"`, its frozen thresholds,
  and its `SoloStrand` semantics;
- a selected chunk has one named `tail.c<N>` section; absent sections mean no evidence;
- each section independently guards its magic/version, chunk molecule count, selected count,
  strictly increasing in-range local ordinals, exact signal count, and absence of trailing bytes;
- the ordinary ten-stream molecule chunk is unchanged, and readers unaware of the optional
  section continue to skip it;
- query resolves ordinal to the decoded molecule, then cell id to the existing packed-barcode
  dictionary. No barcode or UMI is duplicated in the extension.

## Prospective acceptance gates

1. **Default invariance.** On the exact forward/reverse fixture and donor A, ingest without the
   flag is byte-identical to the frozen implementation (whole-file BLAKE3).
2. **Exact attachment.** Fixture and donor-A archive-only decoded witnesses are set-identical to
   the BAM-derived pilot oracle under the frozen witness key. Counts alone are insufficient.
3. **Orientation.** A forward trailing-A and reverse leading-T record describing the same
   transcript-oriented event decode to the expected strand-specific cleavage positions. Wrong-edge
   and wrong-base controls do not attach.
4. **Guards.** Tests reject unsupported layout/version, out-of-range or repeated ordinals,
   signal-count mismatch, truncated payload, unexpected section without matching meta, missing
   declared sections, and trailing bytes.
5. **Archive-only reduction.** Candidate/site-set and group reduction runs using only `.aie`, a
   site-set TSV, and a barcode-to-group TSV. It matches the sidecar oracle without original BAM or
   redundant barcode/UMI fields.
6. **Privacy/scope.** No optional section or metadata field contains bases, qualities, general
   sequence, read names, raw barcodes, or UMI values.
7. **Size.** Optional compressed bytes, including directory/meta growth, are at most 0.25% of the
   donor-A archive and smaller than the corresponding `AIETAIL1` sidecar.
8. **Runtime.** Median of three donor-A runs has at most 10% ingest wall-time overhead. Report
   extraction/write components, archive-only candidate/group query wall time, and peak RSS.
9. **Pilot preservation.** Donor-A external/control counts, tail-positive odds ratio, clip AUC,
   continuous-tail AUC, internal-priming split, and NTRK2 witness result reproduce the frozen pilot
   semantics. Promotion cannot turn the failed NTRK2 secondary gate into a positive claim.

The production layout is promoted only if every applicable gate passes. A structural-hard-stop
result supersedes the size/runtime measurements because benchmarking a lossy substitute would not
answer the promotion question.
