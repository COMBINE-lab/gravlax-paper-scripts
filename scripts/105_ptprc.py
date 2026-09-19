#!/usr/bin/env python
"""PTPRC (CD45) vignette: per-population 3'-site usage from the index, both PBMC datasets."""
import os
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

P=os.environ["GRAVLAX_PROJECT_ROOT"]
AIE=P+"/src/target/release/aie"
SPAN="chr1:198629899-198759346"
plt.rcParams.update({"font.size":9,"font.family":"sans-serif","axes.spines.top":False,
    "axes.spines.right":False,"axes.spines.left":False,"svg.fonttype":"none","axes.axisbelow":True})
fig,axes=plt.subplots(1,2,figsize=(6.4,2.4),sharey=True)
for ax,(tag,arch,gfile) in zip(axes,[("D0",P+"/runs/archive/d0.aie",P+"/runs/replay/apa-groups-d0.tsv"),
                                     ("D1",P+"/runs/archive/d1.aie",P+"/runs/replay/apa-groups-d1.tsv")]):
    r=subprocess.run([AIE,"query",arch,"apa",SPAN,"--strand","+","--tsv","--groups",gfile],
                     capture_output=True,text=True)
    hdr=None; sites=[]
    for line in r.stdout.splitlines():
        f=line.rstrip("\n").split("\t")
        if line.startswith("#"): hdr=f; continue
        sites.append((int(f[1]),int(f[hdr.index("T")]),int(f[hdr.index("M")])))
    sites.sort()
    # keep sites with >=1% of either population's mass
    tT=sum(s[1] for s in sites); tM=sum(s[2] for s in sites)
    keep=[(p,t/tT,m/tM) for p,t,m in sites if t/max(tT,1)>=0.01 or m/max(tM,1)>=0.01]
    pos=np.array([(p-198629899)/1e3 for p,_,_ in keep])
    ax.bar(pos-0.8,[t*100 for _,t,_ in keep],1.5,label="T cells",color="#1a6faf",align="edge")
    ax.bar(pos+0.7,[m*100 for _,_,m in keep],1.5,label="monocytes",color="#c98f1a",align="edge")
    ax.set_xlim(-4,133)
    ax.set_xlabel("3′-site position in PTPRC (kb from gene start)",fontsize=8)
    ax.tick_params(axis="x",labelsize=7.5)
    ax.set_title(tag,fontsize=9)
    ax.yaxis.grid(True,color="#e3e3e3",lw=0.7)
    ax.tick_params(length=0)
axes[0].set_ylabel("share of population's\nPTPRC 3′ mass (%)")
axes[1].legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(P+"/paper/figs/fig_ptprc.svg")
print("PTPRC vignette written")
