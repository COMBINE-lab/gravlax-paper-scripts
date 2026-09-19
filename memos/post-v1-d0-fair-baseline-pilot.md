# Post-v1 D0 fair-baseline pilot

Date: 2026-08-30  
Verdict: correctness gate passed; repeated performance gate remains open.

## What was removed

The reduced BAM retains alignment reference/start, flags, CIGAR, QNAME, and `CR`, `CY`, `UR`,
and `NH`. It sets MAPQ to zero, removes read sequence and quality, and removes every other
auxiliary tag. These are exactly the fields consumed by the current Gravlax BAM extraction path,
but equivalence was established by output comparison rather than assumed from the field audit.

## Correctness result

- A development-set preflight reduced 288,389,557 bytes to 125,298,035 bytes and produced
  byte-identical full- and reduced-BAM archives.
- On D0, independent full- and reduced-BAM ingests produced the same 111,079,672-byte archive
  with SHA-256 `a664edd16cecc86347341bff177f184f67ad53007d3c194ed8b41f0f4d49c82c`.
- Stamping the reference-genome signatures into the reduced-input archive reproduced the frozen
  canonical `d0.aie` byte-for-byte. The initial 7,709-byte discrepancy was entirely the absent
  genome-signature metadata, not molecular content.
- The reduced BAM and its CRAM decode to the same complete headerless SAM stream, SHA-256
  `fbcf85ca51d6f24318a6921fc0ff47aba2b0e7e9a220c669b811dd7674b66933`.
- Archive replay and direct replay from both BAMs reproduced the frozen reference matrix,
  features, and barcodes byte-for-byte.

## Storage result

| artifact | bytes | relative to `.aie` |
|---|---:|---:|
| all-tag annotation-free BAM | 4,208,722,967 | 37.89× |
| all-tag archive-mode CRAM 3.1 | 1,404,884,801 | 12.65× |
| ingest-equivalent reduced BAM | 1,843,272,819 | 16.59× |
| ingest-equivalent reduced CRAM 3.1 | 758,514,044 | 6.83× |
| genome-stamped `.aie` | 111,087,381 | 1.00× |

Removing unused content shrinks BAM 2.28× and CRAM 1.85×. Gravlax remains 6.83× smaller than
the reduced CRAM. This is already a fairer comparison, but the reduced alignment files remain
information-richer: they retain raw barcode correction inputs, names, and every alignment
instance. Call this **ingest-equivalent**, not fully function matched.

## Preliminary runtime signal

One warm/uncontrolled pass took 60.97 s from the all-tag BAM, 58.68 s from the reduced BAM, and
2.52 s from `.aie`. The reduced BAM cut combined user+system CPU time from 140.06 s to 96.66 s,
but wall time by only 4%. This points to record parsing and evidence construction, rather than
raw byte volume alone, as the dominant BAM replay cost. These observations are not manuscript
numbers until the frozen randomized warm/cold replication protocol is complete.

## Decision

The D0 field set is accepted for D1 and D2' correctness pilots. Before scaling performance runs:

1. add counts-only STARsolo to the same D0 protocol;
2. implement a CRAM-fed replay arm or explicitly time CRAM-to-BAM/SAM conversion;
3. freeze a post-correction baseline if it can be represented without restoring unused content;
4. execute randomized five-warm/one-cold repeats with physical-I/O instrumentation.

