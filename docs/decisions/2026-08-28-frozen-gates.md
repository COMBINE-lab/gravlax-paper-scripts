# Frozen gate thresholds — recorded BEFORE any measurement

**Date:** 2026-08-28 · **Host:** analysis-host · **Status:** FROZEN

Per the practice that worked on the predecessor likelihood-coresets project (its
`memos/gate0-decision.md`), all pass/stop thresholds are written down before the measurements that
test them, so that a failure cannot be renegotiated after the fact.

There is no acceptable "continue because the idea is interesting" outcome after a fired STOP
trigger. A MARGINAL result requires an explicit written capability argument in the gate memo,
signed off before proceeding.

---

## G-B · Annotation divergence — is the capability worth anything?

Oracle: STARsolo on D0 (pbmc_1k_v3) under A1=v49, A2-far=v32, A2-near=v48.

| Criterion | PASS | STOP |
|---|---|---|
| UMI mass changing gene assignment, v32→v49 | ≥ 3% | < 1% |
| Expressed genes changing count > 10%, v32→v49 | ≥ 5% | < 2% |

STOP fires only if **both** stop conditions hold. Also report, without thresholds: cells with >1%
total count change; the taxonomy of causes; and the same metrics for v48→v49 (the harder near pair).

*Prior signal (measured 2026-08-28, before running): GENCODE v32 → v49 grows 60,609 → 78,691 genes
(+29.8%) and 227,462 → 507,365 transcripts (+123%). v48 → v49 is flat in genes (78,686 → 78,691)
but +31% in transcripts. Divergence is expected to clear the bar comfortably; the measurement is
still required because transcript growth need not move a 3′-biased count matrix.*

## G-C · Alignment invariance — highest-risk cheap gate

One base index (no GTF). Three mapping configs: (i) no sjdb + `--twopassMode Basic`,
(ii) `--sjdbGTFfile` A1, (iii) `--sjdbGTFfile` A2-far. Compare per-read annotation-independent
evidence tuples: aligned blocks, splice junctions, strand, primary placement.

| Criterion | PASS | STOP |
|---|---|---|
| Gene-assigned reads with identical tuples, (i) vs (ii) | ≥ 99.0% | < 97% |
| Gene-assigned reads with identical tuples, (i) vs (iii) | ≥ 98.5% | < 97% |

A STOP here means annotation-free ingest has an irreducible fidelity ceiling that no amount of
stored evidence can repair, and the project stops on day 2.

## G-D · Annotation-independent UMI grouping — an error bound, not a capability

On the *same* A1 alignments: G1 = STARsolo per-(cell,gene) collapse (`CB`/`UB`/`GX`);
G2 = collapse on cell + UMI + genomic locus + strand.

| Criterion | PASS | STOP |
|---|---|---|
| Median per-cell L1 relative error, G1 vs G2 | < 1% | — |
| G1↔G2 discrepancy vs the A1↔A2 annotation delta from G-B | strictly smaller | ≥ annotation delta |

The second row is the criterion that matters: our grouping error must sit **below the signal we
claim to capture**. Prior art (UMI-tools 2017, sc-SPLASH BKC 2026, UMI-nea 2025) means
gene-agnostic collapsing is not itself a contribution; only the bound is.

Also report: fraction of (CB,UMI) keys spanning >1 gene, and >1 genomic locus.

## G-E · Storage projection — the gate that killed the predecessor project

Real E1 archive for D0 with a real encoder (interned dictionaries, delta coding, zstd).

Baselines, all on the identical dataset:
- compressed FASTQ — **5.0 GB** (measured, on disk)
- **CRAM 3.1** of the annotation-free BAM, `--capture-all-tags` so CB/UB survive — primary baseline
- Genozip
- Malva index — **375 MB** (measured on this host, `${MALVA_ROOT}/indices/pbmc_1k_v3`)

| Criterion | PASS | MARGINAL | STOP |
|---|---|---|---|
| Archive size vs CRAM 3.1 | ≥ 10× smaller | 5–10× | < 5× |
| Archive size vs Malva index | ≤ 375 MB | — | — |

*Honest prior estimate: 15,519,436 assigned UMIs (measured from the existing simpleaf run), but the
archive must also retain intronic/intergenic/antisense molecules that gene-level quantification
discards — precisely because a future annotation may claim them. At ~30M records × 12 B ≈ 360 MB,
which is roughly Malva's 375 MB. The bar is tight and that is the expected outcome, not a surprise.*

## Gate -1A/B/C · Replay fidelity (source doc thresholds, unchanged)

- **-1A** A1 self-replay: > 99.5% agreement on high-confidence uniquely assignable molecules;
  very high gene-level UMI count agreement; every discrepancy classifiable.
- **-1B** A2 cross-replay: strong accuracy specifically on the annotation-sensitive subset
  S₁₂ = {m : Q(D,A1) ≠ Q(D,A2)}. Overall agreement is not acceptable evidence on its own.
- **-1C** A3 unseen features: high precision/recall of supporting molecules per introduced change class.

Correlation (Pearson/Spearman) alone is never an acceptable metric at any gate.

## Gate 0 · Full audit

Entered only on a G-E **PASS**, not on MARGINAL.
