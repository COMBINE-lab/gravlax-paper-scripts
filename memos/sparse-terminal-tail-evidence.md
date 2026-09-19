# Sparse terminal-tail evidence: passing biological pilot, failed compact promotion

## Executive decision

The sparse terminal-tail **biological pilot passes** its prospectively frozen primary gates across
all eight GSE234790 SEZ donors. The signal is strong enough to support an optional,
annotation-independent cleavage-candidate quality-control capability: oriented untemplated tails
are enriched at external PolyASite candidates, improve candidate discrimination beyond soft-clip
length, and are depleted at genome-flagged internal-priming sites for 0.1261% of the matching
archive bytes.

Two deliberately separate conclusions remain negative:

1. The NTRK2 secondary interpretation check has only three witnesses, all proximal and all in one
   donor. It cannot validate the existing NTRK2 terminal-group contrast or support a molecule-level
   TrkB.T1/TrkB.FL claim.
2. The proposed production encoding—one selected archive molecule ordinal followed by one 16-bit
   signal—is provably lossy under the current molecule reducer. Its frozen structural hard stop
   fired before size or timing measurements. No experimental archive code should be merged.

These outcomes are not contradictory. The read-level evidence is useful; a particular compressed
attachment model cannot represent it exactly.

## What was frozen before measurement

The pilot specification was committed as
`ff0cecc75c31aaa9d1c67f47de243cca0ab15500` on Gravlax base
`e378516a0fb6602e63d626d0c49fb11e58f463c8`. Under 10x 3′
`SoloStrand=Forward`, only the transcript-3′ terminal soft clip is inspected: trailing A on a
forward alignment and leading T on a reverse alignment. A positive read requires a clip of at
least 6 nt, at least 80% oriented tail base, and a transcript-distal homopolymer run of at least
4 nt. Cleavage anchors are 0-based exclusive aligned ends on `+` and 0-based inclusive aligned
starts on `-`.

Evidence is deduplicated by
`(reference, cleavage anchor, strand, corrected barcode, UMI)`, retaining the strongest signal.
The experimental `AIETAIL1` sidecar stores fixed 16-byte records containing only reference,
anchor, strand, packed corrected barcode/UMI, and bounded clip/tail/run statistics. It stores no
bases, qualities, read names, aligned sequence, or general soft clips, and it does not change the
default `.aie` path.

The primary gates required correct orientation/packing, at least 5-fold external-site enrichment
with at least 100 witnesses, at least +0.05 AUC over clip length, depletion at internal-priming
candidates, and no more than 1% storage. The NTRK2 check was explicitly secondary.

## Eight-donor result

Across 715,967,790 mapped primary reads, 510,389,785 reads had an eligible called-cell barcode and
UMI. Of 2,932,436 eligible terminally soft-clipped reads, 144,673 passed the tail rule and yielded
136,763 exact-key witnesses. There were 53,346 witnesses near external PolyASite candidates.

After exact-key deduplication, external candidates had 53,345 positives among 186,470 witnesses
(28.61%), versus 763 among 20,386 shifted controls (3.74%), for an odds ratio of 10.299. The
continuous oriented-tail score increased AUC from 0.526843 for clip length alone to 0.654823, a
gain of 0.127980. Genome-flagged internal-priming candidates had 2,085 positives among 23,177
witnesses (9.00%), versus 51,260 among 163,293 non-flagged witnesses (31.39%).

The result is donor-robust under the frozen gates. Every donor exceeds both 5-fold enrichment and
+0.05 AUC improvement, and every donor shows internal-priming depletion. The weakest enrichment
is 6.962-fold and the weakest AUC gain is 0.087829. The largest per-donor storage fraction is
0.137969%.

The eight sidecars total 1,671,772 bytes against 1,326,036,691 archive bytes (0.126073%). Their
uncompressed record payload is 2,188,208 bytes. Standalone, single-threaded audits took 25.49 to
127.56 seconds per donor, 400.09 seconds in aggregate; production extraction would be fused into
initial ingest and would not require a second BAM pass.

The extraction point matters. The annotation-free STAR BAM has full CIGAR, strand, CR/UR, and
alignment-oriented sequence, so the statistic is annotation-independent. The sequence-free,
function-matched reduced CRAM cannot reconstruct it. If retained, this bounded evidence must be
computed during the original ingest (or from another sequence-bearing alignment), before that
reduced representation is produced.

## NTRK2 secondary check

Within 25 bp of the nine frozen NTRK2 sites, only donor G contributes witnesses: one Astro/NSC and
two mature-neuron molecules, all in the proximal shorter-coding terminal group. No donor has a
direct distal MANE/full-length-like witness. The observed total of three is far below the frozen
20-witness and six-donor requirements.

This does not invalidate the fragment-mixture NTRK2 result because cleavage-spanning fragments are
a small, distinct subset of fragmented 3′ libraries. It does prohibit using terminal-tail evidence
as NTRK2 confirmation. More generally, a tail-confirmed terminal site does not phase a terminal
exon to coding exons and must not be described as molecule-level TrkB isoform evidence.

## Why one 16-bit signal per molecule ordinal is lossy

The production gate was frozen in
`50b50a2360d11caf85eb3b1c56c16033f6863546`. The proposed representation selected local molecule
ordinals and stored one 16-bit tail signal per ordinal. It was conditional on exact recoverability
of the frozen read-level witness key.

The actual molecule reducer groups reads by corrected cell, chromosome, strand, locus, UMI, and
junction chain, then retains only two span representatives and a read weight. A tail-positive
middle read can therefore be summarized by an emitted `MolRec` whose retained geometries do not
contain that read's cleavage endpoint. One molecule can also summarize several distinct tail
anchors. This is missing information and one-to-many multiplicity, not an implementation bug.

The executable counterexample in
`d8b7482e1995433d9b017a098ab65b1c700cc2f2` passes tagged records through BAM I/O, barcode
correction, shape interning, locus/chain reduction, renumbering, and final molecule sorting:

- on `+`, spans `90..190`, `100..140`, and `110..130` collapse to one weight-3 molecule; the
  middle read has trailing `10S` A and cleavage anchor 140, while retained endpoints are 190 and
  130;
- on `-`, starts 290, 300, and 310 collapse likewise; the middle read has leading `10S` T and
  anchor 300, while retained starts are 290 and 310.

Using `MolRec::anchor()`, selecting the closest representative, retaining only the strongest
signal, or unioning candidate hits all change the candidate-level evidence. The structural gate
therefore hard-stops the ordinal-plus-one-signal layout. No production section was written, and
size/runtime gates were correctly not measured for a lossy substitute.

## Minimal exact redesign

The smallest exact representation is a sparse **per-selected-molecule event list**:

1. Carry every positive contributing read's `(cleavage anchor, u16 signal)` through molecule
   reduction.
2. Deduplicate globally by `(cell, UMI class, chromosome, strand, cleavage anchor)`, retaining the
   frozen strongest signal. Attach a duplicate witness to one deterministic canonical record.
3. Sort molecule records and their event lists together.
4. For each selected archive chunk, encode strictly increasing delta-coded local molecule
   ordinals, an event count per selected ordinal, and each event as a signed varint anchor delta
   from `MolRec::anchor()` followed by the existing 16-bit signal.
5. Optionally add event-coordinate extents/postings so site queries can skip chunks whose molecule
   sort anchors are remote from the requested cleavage sites.

This design remains annotation-independent and sequence-free, but it permits multiple exact
events per molecule. Its size and ingest overhead must be measured prospectively; the pilot's
0.1261% sidecar ratio must not be presented as the cost of this unimplemented archive layout.

## Claim discipline and next decision

Safe claim: a bounded terminal-tail audit provides a strong, donor-consistent signal for
external-candidate quality control and sparse cleavage-site discovery in 3′ tagged data.

Unsafe claims: that the optional evidence validates NTRK2, phases coding isoforms, is already part
of `.aie`, or can be represented exactly by one signal per current molecule record.

The current decision is to preserve the successful `AIETAIL1` evidence pilot as an optional,
rebuildable result and leave default Gravlax unchanged. Production work should resume only after
explicitly approving and prospectively gating the event-list representation. Exact metrics and
counterexample structure are in `results/post-v1-sparse-terminal-tail.json`; excluded inputs and
outputs are pinned in `manifests/sparse-terminal-tail-artifacts.tsv`.
