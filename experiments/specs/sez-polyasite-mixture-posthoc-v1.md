# Protocol-aware SEZ PolyASite mixture audit

The registered endpoint analysis exposed a method failure rather than a
biological failure: in fragmented 10x 3′ libraries, aligned boundaries occupy
hundreds of bases upstream of a polyadenylation site. Single-link endpoint
clusters therefore chain across broad terminal regions and cannot be promoted
as cleavage-site calls.

The amended analysis retains the frozen eight-donor Astro/NSC to mature-neuron
contrast, replicate definition, eligibility thresholds, and inference. It
changes the measurement model. PolyASite supplies orthogonal candidate site
identities; GENCODE supplies only unique terminal regions. Genomic A-rich
candidates are removed. A label-blind empirical distribution of
transcript-upstream fragment distance is learned from molecules compatible
with exactly one candidate. Each donor is then deconvolved by EM with the
kernel learned from the other seven donors. Thus group labels neither define
sites nor train the observation model.

This is explicitly a post-hoc, method-amended analysis. The original global
directional threshold remains binding and cannot be rescued. The biological
result may be promoted as exploratory, donor-replicated evidence only if the
cross-fitted 12/24/48-bp catalogue sensitivities, shuffled-label null,
held-out predictive checks, simulation recovery, and performance gates pass.
