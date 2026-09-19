# Promotion test for hybrid EM

**Locked:** 2026-08-31, before any new D0/D4 grid run, any direct hybrid-versus-pooled
measurement, or any cell-clustered paired interval. D0 and D4 are development datasets; D4's
earlier untouched test is complete and fully disclosed. D3 remains untouched.

## Question

The original depth hybrid had lower mean negative log loss than the posterior-mean Dirichlet
proxy on D0 and D4, and on D4 it also had slightly better top-1 accuracy. Its D4 Brier score was
worse by 1.32e-5, however, and it had never been compared directly with pooled EM using paired
target-level uncertainty. The present test asks the decision-relevant question: is there one
predeclared monotone transition that is a better estimator than the current pooled default and
the proxy, rather than merely winning one aggregate metric?

## Development family and selection

The fixed convex and proxy constituents remain unchanged. For initial fitted unique evidence
depth `d`, evaluate

`s(d) = d^p / (d^p + D^p)` and
`q_hybrid = (1-s(d)) q_convex + s(d) q_proxy`.

The grid crosses `D` in `{2,4,8,16,32,64,128,256}` with `p` in
`{0.5,1,2,4,8}` on the already defined D0 and D4 masks (seeds 7, 17, and 29) and frozen real
group maps. A configuration is eligible only if, separately in both datasets, its target-count
weighted mean NLL is below both pooled and proxy, its mean Brier is no greater than both, and its
top-1 accuracy is no more than 0.10 percentage point below the better reference. Among eligible
configurations, select the minimum target-count-weighted NLL over all six dataset/seed runs;
break an exact tie by lower power and then lower scale.

This is a joint D0/D4 development selection, not a new confirmation. If no configuration is
eligible, stop the simple depth-only family without relaxing the constraints after inspection.
A later model may add a predeclared disagreement or confidence variable, but it requires its own
locked development rule.

## Paired audit and untouched confirmation

Rerun the selected configuration together with all reference modes so every masked target has
candidate and reference losses from the same cell and seed. Report mean paired NLL and Brier
differences and 95% normal intervals using cells as independent clusters. Promotion development
requires the upper bound to be below zero for NLL and at most zero for Brier against both pooled
and proxy in D0 and D4, with the same top-1 and 4-GiB peak-RSS guardrails.

Only after the model, parameters, and D3 gate are committed and pushed may D3 be grouped and
evaluated once. The primary D3 comparison is hybrid versus the existing pooled default: both
paired NLL and Brier upper confidence limits must be at most zero (strictly below zero for NLL).
The hybrid must also beat the proxy in mean NLL without worse mean Brier, preserve top-1 within
0.10 percentage point of either reference, and real groups must beat their frozen shuffled
control on NLL in at least two seeds. Passing promotes the hybrid as the default evaluator and
licenses an opt-in recovered-count emission mode. Failure retains pooled as the default and
labels the hybrid experimental.

These masked-target proper scores evaluate recovery of deliberately hidden molecular evidence;
they do not establish biological truth or fresh-aligner equivalence.

## Transparent protocol amendment

After the first original-hybrid D0 baseline seed had completed—but before grid selection or any
selected-candidate audit—the across-seed interval calculation was made explicit. Seed masks reuse
cells and partially overlap targets, so the point estimate is target-count weighted and its
standard error is the weighted sum of the three cell-clustered standard errors, equivalent to a
conservative perfect-positive-correlation assumption. At amendment time, the only new paired
result observed was D0 seed 7 hybrid versus pooled: NLL difference -0.00179157 (95% CI
[-0.00264674,-0.000936401]) and Brier difference -0.00110260
([-0.00172448,-0.000480716]). This amendment does not change the grid eligibility rule.

### Robustness amendment after the first selected audit

The original objective selected `D=32,p=2`. It passed 12 of 13 development-audit criteria; its
D4 mean Brier difference versus proxy was favorable (-3.94355e-5), but the conservative interval
included zero ([-1.00325e-4,2.14544e-5]). Before paired results were computed for any other grid
candidate, one final development selection was locked: among the 24 configurations satisfying
the original mean-score eligibility constraints, maximize the smaller relative mean Brier gain
versus proxy across D0 and D4, breaking a tie by lower joint NLL, lower power, and lower scale.
Audit that candidate with the identical paired criteria. If it fails, stop the simple depth-only
family. This is transparently post-hoc model development on D0/D4; only D3 can confirm it.

The robustness rule selected `D=8,p=8`. In its paired audit, all 13 criteria passed. Against
pooled, its conservative NLL/Brier intervals were [-0.00283,-0.000862] and
[-0.00194,-0.000467] on D0, and [-0.00318,-0.00157] and [-0.00177,-0.00109] on D4. Against the
proxy, both proper-score intervals also excluded zero in the favorable direction on both data
sets. The model and the original D3 acceptance rule are therefore frozen. D3 grouping uses the
unchanged candidate-excluded replay-only protocol in `depth-hybrid-em-confirmation.yaml`; no D3
EM score, group partition, or fresh-v49 quantification has yet been observed.

## D3 result

The formal D3 verdict is **FAIL**, and is not relabeled: the all-mode and hybrid-only runs reached
7,807,888 KiB peak RSS, exceeding the frozen 4-GiB cap. Every scientific criterion passed. The
hybrid beat pooled by -0.00247187 NLL (conservative 95% CI [-0.00359083,-0.00135291]) and
-0.00168846 Brier ([-0.00252929,-0.000847637]); it also beat the proxy by -0.000208191 NLL
([-0.000244221,-0.000172162]) and -0.0000905761 Brier
([-0.000109928,-0.0000712242]). Top-1 guardrails passed, and real groups beat shuffled groups in
all three seeds with 1.074% relative NLL gain. Thus D3 is a scientific replication coupled to an
operational memory failure. Any memory remediation is a separately registered engineering test,
not a second untouched statistical confirmation.

## Separately registered memory remediation

The predeclared disk-support-shard gate subsequently **passed** at Gravlax commit
`93def86f4226ecbc3a8a1fbdfa6ca16a99127fd1`. The large D3 shuffled-control seed-7 run retained
byte-identical metrics and stdout while reducing peak RSS from 7,653,508 to 3,846,008 KiB
(49.75%) at 35.86 s wall time, below the locked 4-GiB and 60-s limits. All 106 Rust tests passed;
the D0 in-memory regression was also byte-identical and used 794,232 KiB. The implementation
spills exact batch-compacted support records into 64 temporary cell-shard streams, finalizes one
shard at a time, reduces the spill-path decode window to eight access units, releases archive
dictionaries before finalization, and applies a command-local two-arena glibc allocator bound.
Temporary support files are removed on success or error.

This engineering PASS removes the observed operational blocker but does not retroactively change
the frozen D3 verdict. Taken together, the untouched scientific replication and independently
registered memory result support promoting `D=8,p=8` as the preferred masked-recovery evaluator.
They do not yet validate replacing pooled EM for emitted biological count matrices; that remains
a separate downstream-output gate.
