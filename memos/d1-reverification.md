# Memo — D1 re-verification of all headline numbers

**Date:** 2026-08-29 · **Dataset:** D1 = pbmc_5k_v3 (10x 3′ v3; 383.9M reads, 4 lanes, 26.6 GB
FASTQ R1+R2; 5,038 filtered cells) · **Driver:** `scripts/90_d1_reverify.sh` → `logs/d1-*.log`,
`results/gate1a-d1-v49.json`. Everything at full scale on `.aie` v1.2; project rule satisfied:
no headline number now rests on a single dataset.

## Verdict table — D0 vs D1

| Headline | D0 (pbmc_1k_v3) | D1 (pbmc_5k_v3) | Holds? |
|---|---:|---:|---|
| Molecules / classes | 22.16M / 19.27M | 108.7M / 91.3M | 4.9× scale |
| Archive size | 113.73 MB | 553.91 MB | — |
| **Bits per molecule** | **41.1** | **40.8** | **YES — format cost is scale-stable** |
| vs CRAM 3.1 (all tags) | 12.35× | **11.78×** (6.52 GB) | YES |
| vs BAM | 37.1× | 36.3× (20.09 GB) | YES |
| vs FASTQ (R1+R2) | 44.0× | 48.1× | YES |
| Replay fidelity (UMI mass moved vs fresh oracle, filtered cells) | 0.217% | **0.307%** | YES (≈8× under the 2.43% annotation signal) |
| Cells >1% changed | 0/1,225 | 10/5,038 (0.20%) | YES |
| Genes present, both sides | 78,691 = oracle | 78,691 = oracle | YES |
| Regression (archive replay ≡ BAM replay) | byte-identical | **byte-identical** | YES |
| Replay vs fresh STARsolo | ~50 s vs ~16 min (~19×) | 187 s vs 1,225 s (**6.5×**) | YES, margin narrower (below) |
| aie ingest wall | ~2.4 min | 13.3 min | scales ~linearly |
| Region/junction query | ~0.1 s (open 0.17 s) | ~1.0 s (**open 0.69 s**) | works; open cost scales (below) |

## Honest notes

- **Fidelity 0.307% vs 0.217%.** Same order, direction expected: D1 is deeper per cell (~76k vs
  54k reads/cell), so barcode-collision and correction mass grow. Decomposition: 0.394% oracle
  mass lost by replay vs 0.220% gained; net −0.17% total UMI. Still ~8× below the annotation
  signal the archive exists to capture. Not re-decomposed into
  barcode/alignment/assignment shares here — D0's decomposition (barcode correction dominant)
  presumably carries, but that is an inference, not a measurement.
- **Replay speedup 6.5×, not ~19×.** Both are honest wall clocks: D1's oracle ran at 24 threads on
  an idle host (1,225 s), and replay (187 s) is single-threaded rows processing whose absolute
  time scales with molecules. The speed claim to publish is the qualitative one — re-quantify in
  minutes without FASTQ or realignment — plus both measured pairs, not a single ratio.
- **Open cost 0.69 s at D1 scale.** Eager dictionary decompression (cellofclass is now 137.2 MB
  compressed, 24.8% of the file) scales with archive size. Queries remain ~1 s total. This is the
  measured, growing motivation for lazy dictionary loading — decode postings first, dictionaries
  on demand — already on the backlog (D25/D27).
- **Storage ladder (D1):** FASTQ 26.63 GB → BAM 20.09 GB → CRAM 3.1 6.52 GB → `.aie` 553.9 MB,
  with all molecule evidence, multimapper patterns, UMI graph, and indexes included.
- Sections at D1 scale: chunks 337.5 MB, cellofclass 137.2, patterns 29.9, edges 23.8,
  shapes 11.0, cells 5.1, indexes ~9.3.

## Verdict

**GO.** Every headline number reproduces on a second dataset at 4.9× scale; per-molecule cost is
flat (40.8 vs 41.1 bits). The two numbers that moved (fidelity +0.09pp, replay margin) moved for
understood, recorded reasons and neither approaches its STOP territory.
