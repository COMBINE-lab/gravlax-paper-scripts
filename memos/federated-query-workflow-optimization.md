# Federated query workflow optimization

Date: 2026-09-01

This change optimizes the reproducibility driver and summarizer without changing the frozen v1
result. The existing raw screen is under
`runs/post-v1/federated-biological-screen-r1` in the analysis workspace.

## Driver behavior

`scripts/180_run_federated_biology_screen.sh` now interprets `THREADS` as a fixed aggregate
thread budget and `WORKERS` as the maximum number of concurrent chromosome jobs. The defaults are
24 total threads and four workers, giving each Rayon process six threads. If workers exceed the
thread budget, the effective worker count is capped at the budget. Chromosome results, timing, and
stderr are written to same-directory temporary files and renamed only after a successful command;
the JSON is followed by an atomic completion marker containing the query/concurrency signature.
Incomplete commands and legacy outputs lacking the new signature therefore remain resumable and
are not mistaken for a completed optimized query.

The runner supplies the engine's row-denominator pushdown gate as
`--min-row-informative 10` for every T and M row in the PBMC arm and 20 for every called-cell row
in the six-shard arm. It records the sum of individual command wall times and maximum single-job
RSS as `summed_command_wall_seconds` and `maximum_peak_rss_kib`, plus the end-to-end driver's wall
interval (group preparation through summarization) and thread allocation in the compact JSON
result. Legacy results retain their original `sequential_chromosome_wall_seconds` field exactly.

## Count semantics under pushdown

Legacy, non-pushdown inputs retain the exact v1 `events_screened` meaning and output. With
pushdown active, the summarizer deliberately does not relabel the smaller emitted set as the old
screen denominator. It reports:

- `catalogue_recurrent_events`: the sum of `planning.candidate_events`;
- `denominator_qualified_events`: emitted events after row pushdown;
- `legacy_summed_minimum_events: null`, explicitly marked not evaluated because early pushdown
  avoids computing that historical intermediate set; and
- the active `min_row_informative` threshold.

The frozen legacy denominators remain 20,028 PBMC events and 46,767 six-shard events. They are not
synthesized for a pushdown run.

## Streaming summarizer benchmark

Both implementations were run once on the same 394 MiB of legacy chromosome JSON, including the
same STAR-SJ audit, on the same host. `/usr/bin/time -v` reported:

| implementation | wall time | maximum RSS |
|---|---:|---:|
| original full-document retention | 13.33 s | 826,756 KiB |
| incremental event visitor | 12.29 s | 32,336 KiB |

The streaming implementation reduced peak RSS by 794,420 KiB (96.1%, 25.6-fold) and was 1.04 s
faster in this single paired run. More importantly, its memory is bounded by one event, the parser
buffer, duplicate identifiers, passing candidate summaries, and selected-event summaries rather
than both complete event dictionaries.

Reproduce the new measurement from the workspace root with:

```bash
/usr/bin/time -v -o /tmp/federated-stream.time \
  python3 gravlax-paper-scripts/scripts/181_summarize_federated_biology.py \
  --root "$PWD" \
  --screen runs/post-v1/federated-biological-screen-r1 \
  --out-json /tmp/summary.json \
  --out-pbmc /tmp/pbmc4-consistent.tsv \
  --out-tissue /tmp/brain-associated.tsv
```

The generated summary and both TSVs were byte-identical to their tracked frozen counterparts:

- summary: `84ade0154abf9f1a1bd7591593d73b3521f607b22bfc23cc5f05c2b1e1eab05c`
- PBMC candidates: `4caf5db65601894b9a897a6cf6b7e70e4458d0c38d1594873d6a856ef5e724b9`
- tissue candidates: `622a090ce71a4621dfcb1a665f2005e02507d8ea93b927b6066e2a492f643002`

Run the focused parser/count tests with:

```bash
python3 -m unittest tests/test_federated_biology_summary.py
```

## Exact engine and full-workflow gates

The final implementation was checked after removing an initially over-conservative local-catalogue
presence condition: row pushdown now implements exactly the frozen downstream predicate on observed
evidence.  This matters because an event can satisfy the cross-sample recurrence rule and have a
well-supported row in a sample where `present` is false under the stricter local catalogue rule.
FNBP1 is the concrete regression case.  Its six per-sample denominators are 25, 145, 85, 36, 192,
and 32, although two samples mark the locally discovered event definition absent.

On the six-sample chr19 arm, the original executable/result, packed reducer without pushdown, and
the exact row-gated reducer gave:

| implementation | emitted events | wall time | maximum RSS | JSON bytes |
|---|---:|---:|---:|---:|
| original reducer | 2,393 | 5.67 s | 2,696,248 KiB | 15,165,992 |
| packed group reducer | 2,393 | 4.88 s | 1,685,504 KiB | 15,165,992 |
| packed + exact row pushdown | 683 | 4.73 s | 1,738,916 KiB | 4,388,100 |

The packed default JSON is byte-identical to the original.  The 683 pushed-down event objects are
exactly equal, in the same order, to independently applying the locked six-row denominator filter
to the original 2,393 events.  The pushdown measurement followed the default measurement and should
not be interpreted as a cold-cache runtime comparison; its main demonstrated benefit here is
avoiding later reductions and shrinking serialized output 3.46-fold.  The packed reducer itself
lowers this single-run maximum RSS by 37.5% and wall time by 13.9% relative to the recorded original.

The full 24-chromosome workflow then ran with four workers and a fixed 24-thread aggregate budget.
It completed in 45.34 s by the internal driver clock (45.62 s externally).  Both final candidate
TSVs are byte-identical to the frozen outputs: 784 PBMC candidates (44 with minimum absolute delta
at least 0.10; 14 at least 0.15) and 66 tissue candidates.  The row gates retain 4,703 of 22,953
recurrent PBMC events and 10,760 of 54,328 recurrent tissue events.  Raw chromosome JSON falls from
413,115,915 to 102,166,302 bytes (75.3%, 4.04-fold smaller).  Summed per-chromosome query time falls
from 43.78 to 41.54 s for PBMC and from 73.71 to 67.48 s for tissue; maximum single-command RSS falls
from 2,436,292 to 1,554,168 KiB and from 2,696,248 to 1,914,832 KiB, respectively.  Because jobs run
concurrently, these maximum-command figures are not an aggregate simultaneous-memory measurement.

Compact machine-readable metrics and identities are in
`results/post-v1-federated-query-workflow-optimization.json`; the implementation is Gravlax
commit `c6b68cce2945eb4ca5b1ebc55400f7f066c2863a`.
