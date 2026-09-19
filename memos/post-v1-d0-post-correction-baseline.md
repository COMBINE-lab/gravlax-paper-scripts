# D0 post-correction BAM/CRAM baseline

Date: 2026-08-30  
Result: `results/post-v1-d0-post-correction-baseline.json`  
Protocol: `scripts/114_d0_post_correction_baseline.sh`  
Gravlax: `7e8e9abf26a67201a5f57999b723ecd15faaed95`

## Verdict

The fairest storage comparator so far materially narrows, but does not erase, the `.aie` size
advantage. D0's 111,087,381-byte archive is 2.324 times smaller than a 258,151,959-byte CRAM 3.1
archive-mode encoding of the post-correction molecule abstraction (equivalently, `.aie` uses
56.97% fewer bytes). The explicit BAM is 688,313,625 bytes. This replaces the 6.83×
ingest-equivalent-CRAM ratio for claims specifically about post-correction Gene-replay evidence.

The result is stronger scientifically and less dramatic numerically. The ingest-equivalent CRAM
remains a useful, clearly labeled information-richer comparator: it is 758,514,044 bytes because
it retains raw barcode/UMI evidence and every alignment instance. Moving to the post-correction
representation removes 65.97% of those CRAM bytes before `.aie`-specific compression.

## Representation and controls

The archive has opaque global UMI-class ids rather than nucleotide UMI strings. The exporter
therefore uses sequence-free alignment records plus local tags for cell, class, molecule, weight,
kind, group, alternative, and anchor; unmapped suffix records carry explicit 1MM edges. It does
not fabricate `UB` values. The contract is in Gravlax `docs-notes/molecule-bam.md`.

The accepted v3 run passed all of the following:

- BAM, CRAM, and decoded-BAM `samtools quickcheck`;
- byte-identical normalized SAM record streams from BAM and CRAM;
- byte-identical `matrix.mtx`, `features.tsv`, and `barcodes.tsv` against the frozen D0 replay
  reference from both the molecule BAM and the CRAM-decoded BAM;
- the established 10,393,179-UMI, 3,600,217-entry D0 replay totals.

Two fail-closed development runs improved the representation before acceptance. The first exposed
real multimapper patterns with duplicate alternatives at anchor geometry; the exporter now picks
one deterministic anchor without deleting the duplicate. The second showed CRAM normalizing the
undefined MAPQ of unmapped edge records from 255 to 0; the exporter now writes canonical MAPQ 0.
Neither failure was accepted or used for the reported result.

## Claim boundary

Call this artifact **post-correction molecule CRAM**, not generic CRAM and not a universally
interoperable molecule format. Generic SAM tooling can inspect and compress it, but Gravlax's
local-tag contract is required to reconstruct molecule and UMI-edge semantics. It is function
matched for the evidence needed by the tested Gene replay. It does not carry `.aie`'s query
indexes, and the prototype records the reference digest beside the artifact rather than embedding
the archive genome signature.

Accordingly, the defensible D0 storage statement is:

> For the tested post-correction Gene-replay abstraction, `.aie` was 2.32× smaller than a CRAM 3.1
> archive-mode encoding of the same explicit molecule placements and UMI graph.

Generalization to D1 and D2′ is still mandatory before this becomes a multi-dataset headline.

## Performance note

Single-run warm diagnostics were 17.03 s for molecule-BAM replay, 32.22 s for BAM-to-CRAM
compression, and 12.16 s for CRAM decoding; direct CRAM replay is not implemented. These are not
publication timings. The accepted performance protocol still requires randomized warm replicates
and a recorded non-destructive cold-cache measurement. The much faster native `.aie` replay also
benefits from its purpose-built representation and must not be attributed to compression alone.
