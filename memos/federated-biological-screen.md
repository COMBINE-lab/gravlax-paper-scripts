# Federated biological event screen

## Verdict

The event federation is already capable of a stronger biological story than the locked AKR1A1
vignette.  A post-gate exploratory whole-genome screen across four PBMC archives recovered the
exact FYB1 exon-skipping program as its dominant locus and CD47 among the next highest-effect
loci.  Only after coordinate ranking did we check the literature: both were independently
validated by RT-PCR in T cells and monocytes (Ner-Gaon et al., 2023,
doi:10.3389/fimmu.2023.1116392).  The FYB1 direction is identical in all four libraries: cassette
inclusion is 5.7--8.8% in marker-defined T cells and 98.3--99.1% in monocytes.  Direct aggregation
of STAR's separate SJ matrices gives skip-junction fractions of 96.8--97.5% in T cells and
0.6--3.9% in monocytes.  The latter is a reducer audit against a second representation of the
same alignments, not independent biological validation; the published RT-PCR is the orthogonal
validation.

The four-PBMC whole-genome scan took 43.78 s as 24 sequential chromosome queries (2.32 GiB
maximum RSS); the six-shard called-cell scan took 73.71 s (2.57 GiB).  These are sums of
per-chromosome process times, not a parallelized best case.

The six-shard called-cell contrast also recovers established tissue/cell-compartment programs.
MYL6 cassette inclusion is 69.9--77.6% in the two brain-derived samples and 3.3--4.3% in four
PBMCs.  MYL6 is a headline single-cell splicing locus with RNA-FISH validation (Olivieri et al.,
2021, doi:10.7554/eLife.70692).  RPS24 and PPP1R12A, also highlighted by the SpliZ studies, occur
in the strict candidate set.  This comparison is descriptive: the two brain-derived libraries
are one glioblastoma cell sample and one normal brain-nucleus sample, so tissue, malignancy,
cell composition, and nucleus/cell protocol are inseparable.

## New candidate

The highest-effect cross-tissue locus is FNBP1.  Both brain-derived archives include an annotated
183-nt coding exon (cassette usage 94.1--100%); four PBMC archives mostly skip it (4.0--13.1%).
The exact STAR SJ matrices reproduce the switch: the long skip-junction fraction is 0--7.4% in
the brain-derived samples and 88.3--96.2% in PBMC.  The exon is in-frame (61 amino acids) on the
minus-strand FNBP1 transcripts.  A targeted search found extensive FNBP1 isoform annotation and a
mouse glial Fnbp1 splicing report, but not this human normal-brain/PBMC program.  It is therefore
a candidate, not a discovery claim.

Before opening another dataset, the manifest locks an untouched confirmation in the public 10x
2026 normal-human-brain GEM-X Epi Multiome GEX library.  Passing requires at least 50 informative
called-nucleus UMIs, archive cassette usage at least 0.75, STAR SJ inclusion proxy at least 0.50,
and direction agreement.  Failure is retained and reported; the dataset cannot be silently
replaced.

## Prospective confirmation result

The frozen test passes without changing the event, dataset, or any threshold.  The untouched
public library contains 10,261 called nuclei.  An archive made from called-nucleus reads at the
locked window reports 72 inclusion-only and 2 exclusion-only molecules: 74 informative UMIs and
97.3% cassette inclusion.  Reconstructing the original barcode/UMI and cDNA reads and realigning
them against the genome with annotation-free, two-pass STAR produces exactly the same archive
totals.  STAR's independently aggregated SJ matrix contains 61 inclusion-right versus 2 skip
UMIs (96.8%); direct aggregation of corrected CB/UB pairs from the public Cell Ranger BAM gives
62 versus 2 (96.9%).  Every locked criterion passes.

The STAR audit is deliberately targeted: reads were selected from the public BAM at the already
locked FNBP1 window and then realigned without a GTF.  It validates direction and event-reducer
semantics but cannot detect reads that Cell Ranger placed elsewhere, so it is not presented as a
whole-transcriptome mapping comparison.  The defensible biological statement is that the
predeclared FNBP1 normal-brain cassette program prospectively replicates in a third brain sample;
novelty, cell-type specificity, and mechanism remain open.

## What belongs in the paper now

The FYB1/CD47 rediscovery is ready as an external-positive-control capability result, explicitly
post hoc and without a convenience p-value.  MYL6 provides a second, cross-tissue positive
control but should remain descriptive.  FNBP1 can now enter the main text as a prospectively
replicated candidate, with the targeted-alignment and biological-interpretation limits stated.

Large per-chromosome JSON is intentionally excluded under `runs/post-v1/`.  The repository keeps
the runner, exact thresholds, compact summary, and ranked TSV tables.
