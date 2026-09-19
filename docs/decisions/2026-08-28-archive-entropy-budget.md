# Archive entropy budget, and why Gate G-D and Gate G-E are coupled

**Date:** 2026-08-28 · **Status:** design note, to be tested by G-E

## The per-molecule budget

An E1 molecule record carries: cell barcode, UMI, chromosome, strand, aligned blocks, splice
junctions, and supporting read count. Estimating the entropy of each field for D0 (pbmc_1k_v3,
~30M molecules once intronic/intergenic molecules are retained):

| Field | Coding | Est. bits |
|---|---|---|
| UMI (12 bp) | effectively random | **24** |
| cell barcode | sort by cell, delta → amortized | ~8–17 |
| position | sort by (chrom, start), delta-encode; mean gap ~100 bp | ~8–10 |
| molecule shape (block lengths + junction offsets) | interned; 10x 3' reads are overwhelmingly one 91 bp block | ~1–2 |
| strand | 1 bit | 1 |
| n_reads | varint, small values | ~2 |
| | | **≈ 40 bits ≈ 5 bytes** |

At ~30M molecules that projects to **~150 MB**, against Malva's measured 375 MB on the identical
dataset and a CRAM 3.1 that should land in the 1.5–2.5 GB range. That would clear the frozen
G-E bar of ≥10× vs CRAM. It is a projection, not a measurement, and G-E exists to test it.

## The coupling

**The UMI is 24 of those ~40 bits — well over half the archive.**

But the UMI sequence has exactly one job: defining which reads are the same molecule. Once
molecules are grouped, replay never needs the UMI *string* again — it needs the molecule's
evidence and its identity, not its label.

So the UMI is only worth storing if a future annotation might force us to **redo the collapse**.
Per decision D5 it might: STARsolo collapses per (cell, gene), so a different annotation regroups
reads differently, and reproducing that faithfully needs the UMIs back.

This makes the two gates one question:

- **If G-D passes** — annotation-independent grouping tracks the per-gene collapse closely enough —
  then the collapse can be frozen at ingest and **the UMI sequence dropped entirely**, cutting the
  record from ~5 bytes to ~2, a further ~1.6–2×.
- **If G-D fails**, UMIs must be retained to redo the collapse per annotation, the record roughly
  doubles, and the storage claim weakens correspondingly.

G-E must therefore measure **both** variants:

1. `E1+umi` — UMI retained; collapse is replayable per annotation.
2. `E1-umi` — molecule identity only; collapse frozen at ingest.

Reporting both also makes the paper's storage claim honest rather than a single cherry-picked
number, which the prior-art scan flagged as a weakness of the original "5–100×" framing.

## Consequence for the encoder

The writer is parameterised on whether the UMI field is emitted, and the two variants are built
from the same in-memory molecule set so their sizes are directly comparable. No other field
differs between them.
