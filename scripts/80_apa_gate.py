#!/usr/bin/env python
"""Gate APA-2: for the 300 utr3-withheld genes, does the archive's apa query place 3'-site mass
inside the trimmed 500 bp window — usage the w1 annotation cannot see?

Per gene: window = the 3'-most 500 bp of its v49 span (strand-aware; every transcript's terminal
exon was trimmed 500 bp, so the gene-level trimmed region is contained here). PASS (frozen):
>=90% of genes show >=1 site with >=3 UMIs in-window; aggregate in-window mass within 30% of the
19,662 UMIs the w1 oracle measurably lost.
"""
import os
import subprocess
import sys
from pathlib import Path

P = Path(os.environ["GRAVLAX_PROJECT_ROOT"])
AIE = P / "src/target/release/aie"
ARCHIVE = P / "runs/archive/d0.aie"

# Gene spans and strands from the v49 GTF for the manifest's utr3 genes.
targets = set()
for line in (P / "annotations/withheld/w1.manifest.tsv").read_text().splitlines():
    if line.startswith("utr3\t"):
        targets.add(line.split("\t")[1])

spans = {}
for line in open(P / "annotations/gencode.v49.annotation.gtf"):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if f[2] != "gene":
        continue
    i = f[8].find('gene_id "')
    g = f[8][i + 9 : f[8].index('"', i + 9)]
    if g in targets:
        spans[g] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
missing = targets - set(spans)
if missing:
    sys.exit(f"{len(missing)} manifest genes not found in v49")

hit_genes = 0
in_window_mass = 0
for g, (chrom, start, end, strand) in sorted(spans.items()):
    win = (end - 500, end) if strand == "+" else (start, start + 500)
    out = subprocess.run(
        [str(AIE), "query", str(ARCHIVE), "apa", f"{chrom}:{start}-{end}",
         "--strand", strand, "--tsv"],
        capture_output=True, text=True, check=True,
    ).stdout
    got = False
    for row in out.splitlines():
        _, s, e, _, u, _ = row.split("\t")
        s, e, u = int(s), int(e), int(u)
        # Site overlaps the trimmed window.
        if e >= win[0] and s < win[1]:
            in_window_mass += u
            if u >= 3:
                got = True
    if got:
        hit_genes += 1

print(f"utr3 genes with >=1 in-window site (>=3 UMIs): {hit_genes}/300 ({100*hit_genes/300:.1f}%)  [PASS >=90%]")
print(f"aggregate in-window 3'-site UMI mass: {in_window_mass}  (w1 oracle lost 19,662)  "
      f"ratio {in_window_mass/19662:.2f}  [PASS within 30%]")
