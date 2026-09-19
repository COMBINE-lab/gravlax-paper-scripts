#!/usr/bin/env python
"""Gate DISC-3: classify results/discover-w1.tsv candidates against v49 gene spans.

Buckets: (a) overlaps a withheld gene, (b) overlaps another v49 gene, (c) no v49 overlap.
Bucket (c) is further decomposed by strand-aware proximity (<= 5 kb downstream of a v49 3' end)
because 10x 3' readthrough is the expected dominant source of "novel" loci.

Run 2026-08-29 (full D0): 9,388 candidates / 664,615 UMIs -> a 17.0% loci / 33.9% mass,
b 3.2% / 7.1%, c 79.7% / 59.0%; within c, 57.9% of loci (56.9% of mass) are <=5 kb downstream
of a v49 3' end; 3,151 truly intergenic loci (169,116 UMIs). Memo: memos/gate-apa-discover.md.
"""
import os
import bisect

P = os.environ["GRAVLAX_PROJECT_ROOT"]

wh = set(l.split("\t")[1] for l in open(P + "/annotations/withheld/w1.manifest.tsv") if l.startswith("gene\t"))
v49 = {}
for line in open(P + "/annotations/gencode.v49.annotation.gtf"):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    i = f[8].find('gene_id "')
    g = f[8][i + 9 : f[8].index('"', i + 9)]
    v49.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4]), f[6], g in wh))
for ch in v49:
    v49[ch].sort()

def classify(ch, s, e):
    lst = v49.get(ch, [])
    i = bisect.bisect_left(lst, (e, 0, "", False))
    w = o = False
    for j in range(max(0, i - 2000), min(len(lst), i + 1)):
        if lst[j][0] < e and lst[j][1] > s:
            if lst[j][3]:
                w = True
            else:
                o = True
    return w, o

def near_3p(ch, s, e, d=5000):
    lst = v49.get(ch, [])
    i = bisect.bisect_left(lst, (e + d, 0, "", False))
    for j in range(max(0, i - 2000), min(len(lst), i + 1)):
        gs, ge, st, _ = lst[j]
        if st == "+" and 0 <= s - ge < d:
            return True
        if st == "-" and 0 <= gs - e < d:
            return True
    return False

a = b = c = c3 = 0
ua = ub = uc = uc3 = 0
for l in open(P + "/results/discover-w1.tsv"):
    ch, s, e, st, u, cc = l.rstrip("\n").split("\t")
    s, e, u = int(s), int(e), int(u)
    w, o = classify(ch, s, e)
    if w:
        a += 1; ua += u
    elif o:
        b += 1; ub += u
    elif near_3p(ch, s, e):
        c3 += 1; uc3 += u
    else:
        c += 1; uc += u

n = a + b + c + c3
tu = ua + ub + uc + uc3
print(f"candidates {n}: withheld {a} ({100*a/n:.1f}%), other-v49 {b} ({100*b/n:.1f}%), "
      f"3'-readthrough {c3} ({100*c3/n:.1f}%), intergenic {c} ({100*c/n:.1f}%)")
print(f"UMI mass {tu}: withheld {ua} ({100*ua/tu:.1f}%), other-v49 {ub} ({100*ub/tu:.1f}%), "
      f"3'-readthrough {uc3} ({100*uc3/tu:.1f}%), intergenic {uc} ({100*uc/tu:.1f}%)")
