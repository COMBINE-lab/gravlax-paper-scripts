#!/usr/bin/env python
"""Gate G-B: how much does a count matrix change between GENCODE releases?

The archive's whole value proposition is replaying annotation changes. If changing the annotation
barely moves the matrix, there is nothing worth replaying and the project stops. No published study
measures this for scRNA-seq across GENCODE versions, so this is also a result in its own right.

Genes are matched by *unversioned* Ensembl id (ENSG00000123456.7 -> ENSG00000123456): GENCODE bumps
the version suffix on unchanged genes constantly, and matching on the full id would report churn
that is purely cosmetic.

Usage: 30_gateb_compare.py <solo_dir_A> <solo_dir_B> <label_A> <label_B> [--out results/gateb.json]
where each solo_dir is a STARsolo `Solo.out/Gene/raw` directory.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp


def strip_version(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def load(solo_raw: Path):
    """Load a STARsolo raw matrix as CSR (genes x barcodes), with its gene and barcode labels."""
    mtx = scipy.io.mmread(str(solo_raw / "matrix.mtx")).tocsr()
    genes = [strip_version(l.split("\t")[0]) for l in
             (solo_raw / "features.tsv").read_text().splitlines() if l]
    barcodes = [l.strip() for l in (solo_raw / "barcodes.tsv").read_text().splitlines() if l]
    if mtx.shape != (len(genes), len(barcodes)):
        sys.exit(f"shape {mtx.shape} does not match {len(genes)} genes x {len(barcodes)} barcodes")
    return mtx, genes, barcodes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir_a", type=Path)
    ap.add_argument("dir_b", type=Path)
    ap.add_argument("label_a")
    ap.add_argument("label_b")
    ap.add_argument("--cells", type=Path,
                    help="Solo.out/Gene/filtered/barcodes.tsv; restrict to real cells. The "
                         "unfiltered set is ~240k mostly-empty droplets of ambient RNA.")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    A, genes_a, bc_a = load(args.dir_a)
    B, genes_b, bc_b = load(args.dir_b)

    if bc_a != bc_b:
        sys.exit("barcode lists differ; both runs must use the same whitelist")

    if args.cells:
        keep = {l.strip() for l in args.cells.read_text().splitlines() if l.strip()}
        cols = [i for i, b in enumerate(bc_a) if b in keep]
        if not cols:
            sys.exit("no barcodes from --cells found in the matrices")
        A, B = A[:, cols], B[:, cols]
        bc_a = [bc_a[i] for i in cols]
        print(f"restricted to {len(cols)} cells from {args.cells}")

    idx_a = {g: i for i, g in enumerate(genes_a)}
    idx_b = {g: i for i, g in enumerate(genes_b)}
    shared = sorted(set(idx_a) & set(idx_b))
    only_a = sorted(set(idx_a) - set(idx_b))
    only_b = sorted(set(idx_b) - set(idx_a))

    Sa = A[[idx_a[g] for g in shared], :]
    Sb = B[[idx_b[g] for g in shared], :]

    tot_a, tot_b = A.sum(), B.sum()
    # Count mass living in genes that exist in only one annotation. This mass is unreachable from
    # the other annotation's count matrix by any means -- it is the clearest statement of what a
    # matrix-only archive loses.
    mass_only_a = A[[idx_a[g] for g in only_a], :].sum() if only_a else 0
    mass_only_b = B[[idx_b[g] for g in only_b], :].sum() if only_b else 0

    diff = abs(Sa - Sb)
    l1 = diff.sum()
    # Total UMI mass that moves: L1/2 counts each relocated UMI once (out of one gene, into another).
    moved_frac = (l1 / 2) / max(tot_a, 1)

    # The shared-gene L1 above ignores mass sitting in genes that exist in only one annotation --
    # which is precisely the mass a count matrix can never recover. The honest figure is the L1 over
    # the *union* of genes, treating an absent gene as zero.
    l1_union = l1 + mass_only_a + mass_only_b
    moved_union_frac = (l1_union / 2) / max(tot_a, 1)
    # L1 conflates reassignment with net gain, so decompose it: v49 simply has more genes to assign
    # to, and that is a different phenomenon from a UMI moving between two genes both annotations know.
    d = (Sb - Sa)
    gained = d[d > 0].sum() + mass_only_b
    lost = -d[d < 0].sum() + mass_only_a
    reassigned = min(gained, lost)
    net = abs(gained - lost)

    ga = np.asarray(Sa.sum(axis=1)).ravel()
    gb = np.asarray(Sb.sum(axis=1)).ravel()
    expressed = (ga + gb) > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(ga - gb) / np.maximum(ga, 1)
    changed_10 = int((expressed & (rel > 0.10)).sum())
    n_expressed = int(expressed.sum())

    ca = np.asarray(Sa.sum(axis=0)).ravel()
    cb_ = np.asarray(Sb.sum(axis=0)).ravel()
    live = (ca + cb_) > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        crel = np.abs(ca - cb_) / np.maximum(ca, 1)
    cells_1pct = int((live & (crel > 0.01)).sum())

    res = {
        "label_a": args.label_a,
        "label_b": args.label_b,
        "genes_a": len(genes_a), "genes_b": len(genes_b),
        "genes_shared": len(shared), "genes_only_a": len(only_a), "genes_only_b": len(only_b),
        "total_umi_a": float(tot_a), "total_umi_b": float(tot_b),
        "umi_mass_in_genes_only_a": float(mass_only_a),
        "umi_mass_in_genes_only_b": float(mass_only_b),
        "l1_shared": float(l1),
        "l1_union": float(l1_union),
        "umi_mass_moved_frac": float(moved_frac),
        "umi_mass_moved_union_frac": float(moved_union_frac),
        "mass_gained": float(gained),
        "mass_lost": float(lost),
        "mass_reassigned": float(reassigned),
        "mass_net_gain": float(net),
        "expressed_genes": n_expressed,
        "genes_changed_gt10pct": changed_10,
        "genes_changed_gt10pct_frac": changed_10 / max(n_expressed, 1),
        "cells_changed_gt1pct": cells_1pct,
        "cells_live": int(live.sum()),
    }

    print(f"=== Gate G-B: {args.label_a} vs {args.label_b} ===")
    print(f"genes: {len(genes_a)} vs {len(genes_b)}  (shared {len(shared)}, "
          f"only-{args.label_a} {len(only_a)}, only-{args.label_b} {len(only_b)})")
    print(f"total UMI: {tot_a:.0f} vs {tot_b:.0f}")
    print(f"UMI mass in genes present in only one annotation: "
          f"{mass_only_a:.0f} ({100*mass_only_a/max(tot_a,1):.3f}%) / "
          f"{mass_only_b:.0f} ({100*mass_only_b/max(tot_b,1):.3f}%)")
    print(f"UMI mass moved (shared genes):   {100*moved_frac:.3f}%")
    print(f"UMI mass moved (gene UNION):     {100*moved_union_frac:.3f}%   [PASS >=3%, STOP <1%]")
    print(f"  of which gained by {args.label_b}: {gained:.0f} ({100*gained/max(tot_a,1):.3f}%)")
    print(f"  of which lost by   {args.label_a}: {lost:.0f} ({100*lost/max(tot_a,1):.3f}%)")
    print(f"  reassigned (min):              {reassigned:.0f} ({100*reassigned/max(tot_a,1):.3f}%)")
    print(f"  net gain (|gained-lost|):      {net:.0f} ({100*net/max(tot_a,1):.3f}%)")
    print(f"expressed genes changing >10%:   {changed_10}/{n_expressed} "
          f"({100*changed_10/max(n_expressed,1):.3f}%)   [PASS >=5%, STOP <2%]")
    print(f"cells changing >1% total count:  {cells_1pct}/{int(live.sum())}")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(res, indent=2))
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
