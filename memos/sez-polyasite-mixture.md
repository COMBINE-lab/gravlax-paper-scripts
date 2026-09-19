# Protocol-aware PolyASite mixture: result and interpretation

## Verdict

The raw endpoint experiment correctly diagnosed a measurement-model problem: 10x 3′ cDNA is
fragmented, so a read's terminal aligned base is not the RNA cleavage coordinate. Broad endpoint
clusters therefore cannot support a cleavage-site story. The replacement analysis is a
protocol-aware mixture model whose candidate identities come from PolyASite and whose observation
kernel is learned from the archive. It passes every frozen method-amendment control.

The production implementation is Gravlax commit
`e378516a0fb6602e63d626d0c49fb11e58f463c8`; that full identity is also embedded in the compact
result.

This result is exploratory. The original preregistered global-lengthening gate failed: all eight
donor median effects are positive, but the across-donor median is 0.01281 rather than the required
0.02. No later analysis rescues or redefines that gate.

## Accuracy controls

- GENCODE v49 supplies only uniquely assigned, same-strand terminal regions; PolyASite 2.0 supplies
  candidate site identities. Candidates with genomic internal-priming signatures are removed.
- A label-blind empirical distribution of fragment distance upstream of a candidate is learned
  from molecules compatible with exactly one candidate. Each donor is fitted using the kernel
  learned from the other seven donors.
- The held-out kernel beats a uniform distance model in all eight donors, by 0.418--0.598
  bits per unique-candidate UMI. The learned kernel peaks 150--160 bp upstream, as expected for
  fragmented 3′ libraries.
- The primary 24-bp merge analysis contains 45,705 recurrent candidates, 35.7 million assigned
  expected UMIs, 2,182 tested genes, and 95 reported genes. The within-donor shuffled-label arm
  reports zero genes.
- Of the 95 primary genes, 93 (97.9%) retain their effect sign at both the 12-bp and 48-bp
  candidate-merging sensitivities. Those arms report 92 and 96 genes, respectively.
- A deterministic two-site simulation defines the resolution boundary. At 1,000 UMIs, mean
  absolute mixture error falls from 0.054 at 50-bp separation to 0.013 at 1,000 bp; at 200 UMIs it
  falls from 0.098 to 0.029. Close, low-depth sites should therefore be interpreted cautiously.

## Speed

The frozen cross-fitted implementation required 37.16 s and 3,389,336 KiB maximum RSS. Caching
the sparse endpoint-to-candidate likelihood graph, locating compatible candidates by binary
search, and accumulating the pooled kernel from donor kernels reduced the identical primary run
to 22.02 s (1.69× faster) at 3,660,468 KiB maximum RSS. `sites.tsv`, `genes.tsv`, and
`fragment-kernel.tsv` are byte-identical across implementations.

Eight independent legacy `apa-test` scans require 114.00 s and 2,432,640 KiB. The production
cohort command is therefore 5.18× faster while adding a common recurrent catalogue, cross-fitted
protocol model, paired-donor inference, and complete zero-explicit sample/group/site output. The
additional 1.17 GiB peak RSS remains below the frozen 8-GiB gate.

## Concrete biological story: NTRK2 receptor architecture

`NTRK2` is the strongest result in the primary table (effect +0.2819, eight of eight donor signs,
paired t = 32.48, q = 1.48e-5, exact sign-flip p = 0.0078125). A structure-only validation—written
after the method was fixed and therefore explicitly exploratory—clusters protein-coding GENCODE
transcript ends and assigns fitted sites within 2 kb:

- the proximal GENCODE-Primary group ends at chr9:84,811,121--84,815,706 and has a median coding
  sequence of 1,431 bases;
- the distal group containing the MANE Select transcript ends at
  chr9:85,021,896--85,027,054 and has a median coding sequence of 2,466 bases.

Across the nine fitted sites adjacent to those two groups, the distal MANE/full-length-like share
is 20.0% in Astro/NSCs (564 distal versus 2,260 proximal expected UMIs) and 78.1% in mature
neurons (3,146 versus 880). Every donor shifts in the same direction; the median paired change is
+56.8 percentage points and its exact sign-flip p-value is 0.0078125.

This is a compelling receptor-architecture result because human NTRK2 truncated variants use an
alternative terminal exon and lack the kinase-containing 3′ exons, while the full-length receptor
is the canonical neuronal BDNF receptor. Independent work reports predominant truncated TrkB.T1
in astrocytes and predominant TrkB.FL in neurons. Gravlax recovers the matching pattern in adult
human SEZ nuclei in seconds, after the candidate catalogue and annotation are supplied at query
time.

The defensible phrase is “proximal shorter-coding versus distal MANE/full-length-like terminal
group.” PolyASite-mixture counts do not phase a terminal site to every upstream coding exon, so
the present data do not by themselves prove a complete TrkB.T1 or TrkB.FL molecule. Junction or
long-read validation would be needed for literal isoform labels.

## Reproduction

- `scripts/197_build_sez_full_archives.sh`: annotation-free full-archive construction.
- `scripts/200_run_sez_polyasite_mixture_crossfit.sh`: frozen primary, 12/48-bp sensitivities,
  shuffled-label null, and legacy timing.
- `scripts/199_simulate_polyasite_mixture_accuracy.py`: deterministic mixture-recovery surface.
- `scripts/201_validate_ntrk2_terminal_architecture.py`: GENCODE-derived NTRK2 terminal groups.
- `scripts/203_benchmark_polyasite_mixture_hotpath.sh`: optimized rerun and byte-identity guard.
- `scripts/202_summarize_sez_polyasite_mixture.py`: fail-closed compact result.

Compact machine-readable outputs are `results/post-v1-sez-polyasite-mixture.json`,
`results/post-v1-sez-polyasite-mixture.tsv`, and
`results/post-v1-sez-polyasite-mixture-ntrk2.{json,tsv}`. Raw archives, FASTQs, BAMs, and full
site tables remain ignored.
