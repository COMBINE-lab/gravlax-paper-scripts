#!/usr/bin/env python
"""Build a withheld-annotation panel: A_minus = A1 with real, expressed features held out.

This replaces the source plan's synthetic A3 edits (docs/implementation_plan.md §4) with the
held-out-truth design of the oarfish genome-projection / TranSigner evaluations: every withheld
feature is real biology with real expression, its identity is recorded, and the withholding
fraction is a dial. Replaying the FULL annotation from an archive whose oracle ran on A_minus —
and vice versa — measures exactly how well the archive recovers features it has never seen,
stratified by change class.

Classes produced (deterministic under --seed):
  gene      — the entire gene is absent from A_minus (the "new gene in the future" case)
  isoform   — one transcript of a multi-transcript gene is absent (the "new isoform" case)
  utr3      — every transcript of the gene has its 3'-terminal exon trimmed by --utr-trim bp
              (the "3' end moved" case; the highest-signal class for 10x 3' data)

Genes are sampled only from a supplied expression table (oracle matrix gene sums), stratified by
expression tercile within the expressed set, so withheld features are guaranteed to carry signal.

Usage:
  70_withhold_annotation.py A1.gtf --expr gene_sums.tsv --out-prefix annotations/withheld/w1 \
      [--n-gene 300] [--n-isoform 300] [--n-utr3 300] [--utr-trim 500] [--seed 20260828]

gene_sums.tsv: two columns, unversioned gene id \t total UMI count (produce with
scripts/71_gene_sums.py from a STARsolo raw matrix).
"""
import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path


def gene_of(attrs: str):
    i = attrs.find('gene_id "')
    if i < 0:
        return None
    j = attrs.index('"', i + 9)
    return attrs[i + 9 : j]


def tx_of(attrs: str):
    i = attrs.find('transcript_id "')
    if i < 0:
        return None
    j = attrs.index('"', i + 15)
    return attrs[i + 15 : j]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf", type=Path)
    ap.add_argument("--expr", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--n-gene", type=int, default=300)
    ap.add_argument("--n-isoform", type=int, default=300)
    ap.add_argument("--n-utr3", type=int, default=300)
    ap.add_argument("--utr-trim", type=int, default=500)
    ap.add_argument("--min-expr", type=int, default=20,
                    help="minimum oracle UMI count for a gene to be eligible — withholding an "
                         "unexpressed feature tests nothing")
    ap.add_argument("--seed", type=int, default=20260828)
    args = ap.parse_args()

    expr = {}
    for line in args.expr.read_text().splitlines():
        g, n = line.split("\t")
        expr[g] = float(n)

    lines = args.gtf.read_text().splitlines()
    tx_by_gene = defaultdict(set)
    exon_count = defaultdict(int)
    for l in lines:
        if l.startswith("#"):
            continue
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g, t = gene_of(f[8]), tx_of(f[8])
        if g and t:
            tx_by_gene[g].add(t)
            exon_count[t] += 1

    def unv(g):
        return g.split(".", 1)[0]

    eligible = [g for g in tx_by_gene if expr.get(unv(g), 0) >= args.min_expr]
    if not eligible:
        sys.exit("no eligible genes — check --expr ids against the GTF")
    # Stratify by expression tercile so the withheld set spans the dynamic range instead of being
    # dominated by whatever the sampler happens to hit.
    eligible.sort(key=lambda g: expr[unv(g)])
    n = len(eligible)
    tiers = [eligible[: n // 3], eligible[n // 3 : 2 * n // 3], eligible[2 * n // 3 :]]

    rng = random.Random(args.seed)
    picked = set()

    def sample(count, pred):
        out = []
        for tier in tiers:
            pool = [g for g in tier if g not in picked and pred(g)]
            rng.shuffle(pool)
            take = pool[: count // 3 + 1]
            out.extend(take)
            picked.update(take)
        return out[:count]

    drop_genes = set(sample(args.n_gene, lambda g: True))
    iso_genes = sample(args.n_isoform, lambda g: len(tx_by_gene[g]) >= 2)
    utr_genes = set(sample(args.n_utr3, lambda g: True))

    # For isoform withholding, drop the transcript with the most exons that is not the gene's only
    # multi-exon transcript — a real structural isoform, deterministically chosen.
    drop_tx = {}
    for g in iso_genes:
        txs = sorted(tx_by_gene[g], key=lambda t: (-exon_count[t], t))
        drop_tx[txs[0]] = g
    iso_genes = set(iso_genes)

    out_gtf = args.out_prefix.with_suffix(".gtf")
    out_manifest = args.out_prefix.with_suffix(".manifest.tsv")
    out_gtf.parent.mkdir(parents=True, exist_ok=True)

    # utr3: find, per transcript of each utr_gene, the 3'-terminal exon (strand-aware) and trim it.
    # Terminal exons shorter than the trim are left intact rather than deleted, so the transcript
    # keeps its exon count and the change stays a pure end-shift.
    term_exon = {}  # (tx) -> (line_index, start, end, strand)
    if utr_genes:
        for idx, l in enumerate(lines):
            if l.startswith("#"):
                continue
            f = l.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            g = gene_of(f[8])
            if g not in utr_genes:
                continue
            t = tx_of(f[8])
            s, e, strand = int(f[3]), int(f[4]), f[6]
            cur = term_exon.get(t)
            is_term = cur is None or (strand == "+" and e > cur[2]) or (strand == "-" and s < cur[1])
            if is_term:
                term_exon[t] = (idx, s, e, strand)

    trimmed_lines = {}
    for t, (idx, s, e, strand) in term_exon.items():
        if e - s + 1 <= args.utr_trim + 20:
            continue
        f = lines[idx].split("\t")
        if strand == "+":
            f[4] = str(e - args.utr_trim)
        else:
            f[3] = str(s + args.utr_trim)
        trimmed_lines[idx] = "\t".join(f)

    kept, m = [], []
    for idx, l in enumerate(lines):
        if l.startswith("#"):
            kept.append(l)
            continue
        f = l.split("\t")
        if len(f) < 9:
            kept.append(l)
            continue
        g, t = gene_of(f[8]), tx_of(f[8])
        if g in drop_genes:
            continue
        if t is not None and t in drop_tx:
            continue
        kept.append(trimmed_lines.get(idx, l))

    for g in sorted(drop_genes):
        m.append(f"gene\t{g}\t-\t{expr[unv(g)]:.0f}")
    for t, g in sorted(drop_tx.items()):
        m.append(f"isoform\t{g}\t{t}\t{expr[unv(g)]:.0f}")
    for g in sorted(utr_genes):
        m.append(f"utr3\t{g}\t-\t{expr[unv(g)]:.0f}")

    out_gtf.write_text("\n".join(kept) + "\n")
    out_manifest.write_text("class\tgene\ttranscript\toracle_umis\n" + "\n".join(m) + "\n")
    print(f"withheld: {len(drop_genes)} genes, {len(drop_tx)} isoforms, "
          f"{len(utr_genes)} 3'-trims ({sum(1 for _ in trimmed_lines)} exon lines trimmed)")
    print(f"wrote {out_gtf} and {out_manifest}")


if __name__ == "__main__":
    main()
