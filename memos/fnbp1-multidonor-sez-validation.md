# Prospective multi-donor FNBP1 validation in adult human SEZ

## Registered outcome

The fixed FNBP1 cassette event and unchanged broad-cell marker rule were registered and pushed
before any GSE234790 sequence or FNBP1 event evidence was inspected. Donor, not cell or molecule,
is the replicate. The primary brain-replication gate is **STOP_UNDERPOWERED**: only donor G has the
registered minimum of 10 informative molecules, so the required four eligible donors are not
available. This is not a failed replication test, and the eligibility threshold must not be
relaxed after seeing the data.

The microglia-localization gate is **STOP_DESCRIPTIVE**. No donor has any informative molecule in
the frozen microglia/immune group, so none meets the registered requirement of at least four
microglial and ten other-resolved informative molecules. No donor-level sign test was performed.

| Donor | Published-QC nuclei | Read pairs | Include | Exclude | Informative | Usage | Primary eligible | Microglial informative |
|---|---:|---:|---:|---:|---:|---:|:---:|---:|
| A | 4,275 | 66,841,000 | 4 | 0 | 4 | 1.00 | no | 0 |
| B | 1,892 | 64,522,235 | 0 | 0 | 0 | -- | no | 0 |
| C | 2,888 | 66,508,966 | 1 | 0 | 1 | 1.00 | no | 0 |
| D | 7,259 | 297,115,616 | 9 | 0 | 9 | 1.00 | no | 0 |
| E | 1,159 | 52,244,860 | 6 | 0 | 6 | 1.00 | no | 0 |
| F | 2,966 | 64,280,614 | 4 | 0 | 4 | 1.00 | no | 0 |
| G | 14,213 | 310,977,420 | 13 | 0 | 13 | 1.00 | yes | 0 |
| H | 1,974 | 54,821,911 | 2 | 0 | 2 | 1.00 | no | 0 |

Across donors, 39 informative molecules were observed and all were inclusion-only. Seven of eight
donors contributed at least one inclusion molecule. These are useful descriptive observations,
but pooling the 39 molecules would pseudoreplicate the donor-level experiment and cannot convert
the registered STOP into a PASS. Likewise, the absence of microglial event molecules does not
establish absence of a compartment effect.

## Technical audit

All 977,312,622 paired reads from the 40 ENA runs were aligned without an annotation. Full-genome
STAR BAM output was streamed through the fixed FNBP1 window, so no genome-wide BAM was
materialized. Gravlax and STAR's independent SJ-matrix reduction of those same alignments agree
exactly in every donor: 39 inclusion and zero skip UMIs in total. This validates the event reducer
on the prospective panel; it is not an independent-mapper validation.

The eight donor alignments were scheduled concurrently. Per-donor alignment wall times sum to
4,069.42 seconds; the slowest donor took 1,334.03 seconds, and the largest observed per-worker peak
RSS was 37,152,996 KiB (35.43 GiB). The streamed target BAMs total 1,719,237 bytes and the called-cell
target archives total 175,312 bytes; these targeted-window sizes are execution artifacts, not a
whole-dataset storage comparison.

The first parallel wrapper reached the end of all donor work in 22:24.87, then exited 1 because
the original result summarizer misparsed GNU `time`'s elapsed-time label. The donor outputs were
not rerun or discarded. Commit `6a784cb` fixes only that parser. A checkpointed resume exited 0 in
6.49 seconds at 37,820 KiB peak RSS, and a separate direct summarization exited 0 in 6.74 seconds at
38,148 KiB. The donor-group count table was byte-identical across the two successful summaries.
The failed wrapper, successful resume, and direct-summary timings are retained separately rather
than presenting the failed wrapper as a successful end-to-end timing.

## Interpretation and next decision

This cohort cannot answer the registered donor-level biological question at its observed
junction depth. It nevertheless supplies a clean prospective technical audit and descriptive
directional evidence consistent with high inclusion. If the biological gate is pursued again,
the next cohort should be selected using event-blind junction-depth proxies at matched genes and
3-prime distances, with at least four donors predicted to exceed ten bulk informative molecules.
A microglial contrast additionally needs at least six donors predicted to exceed both frozen
compartment minima. A deeper, longer-read, or full-length normal-brain design is more promising
than adding similarly shallow nuclei from this cohort.

## Reproduction

- Prospective lock: `experiments/manifests/fnbp1-multidonor-sez-validation.json`
- Protocol: `experiments/specs/fnbp1-multidonor-sez-validation.md`
- Event-blind feasibility and frozen groups: script `190` and
  `results/post-v1-fnbp1-multidonor-sez-feasibility/`
- Checksum-verified ENA acquisition: script `191` and
  `results/post-v1-fnbp1-multidonor-sez-acquisition.json`
- Streamed donor runner and fixed window: script `192` and
  `experiments/manifests/fnbp1-window.grch38.bed`
- Registered summarizer: script `193`
- Compact results: `results/post-v1-fnbp1-multidonor-sez-validation.{json,tsv}`

Raw FASTQs, STAR outputs, target BAMs, and archives remain excluded from Git.
