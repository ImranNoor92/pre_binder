#!/usr/bin/env python3
"""Bucket C + binding-energy figure (Rosetta InterfaceAnalyzer on as-designed complexes).
2x2, A4-friendly compact, contrasting palette by monomer group. reports/bucketC_interface.{pdf,png,svg}
"""
from __future__ import annotations
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/data/binder_software/pre-binder"
rows = list(csv.DictReader(open(f"{ROOT}/outputs/bucketC_metrics.csv")))
rows.sort(key=lambda r: float(r["dG_separated"]))   # most favorable first
tags = [r["tag"].replace("design_", "d") for r in rows]
x = np.arange(len(rows))

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
    "font.size":8,"axes.titlesize":10,"axes.labelsize":9,"xtick.labelsize":7,"ytick.labelsize":7,
    "legend.fontsize":7,"text.color":"black","axes.labelcolor":"black","axes.titlecolor":"black",
    "xtick.color":"black","ytick.color":"black","axes.edgecolor":"black",
    "pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
GCOL = {"6/6":"#0072B2","5/6":"#009E73","4/6":"#E69F00","3/6":"#CC79A7"}
cols = [GCOL.get(r["group"], "#777") for r in rows]
CUT = "#b3001b"

fig, ax = plt.subplots(2, 2, figsize=(7.3, 6.6))
fig.suptitle("Bucket C — interface geometry & binding energy (Rosetta, as-designed poses, n = 20)",
             fontsize=10.5, fontweight="bold", y=0.995)

def bars(a, key, title, ylabel, cut=None, cutlabel=None, conv=float):
    v = np.array([conv(r[key]) for r in rows])
    a.bar(x, v, color=cols, edgecolor="black", linewidth=0.3)
    if cut is not None:
        a.axhline(cut, color=CUT, ls="--", lw=1.3, label=cutlabel); a.legend(frameon=False, loc="best")
    a.set_xticks(x); a.set_xticklabels(tags, rotation=90, fontsize=7)
    a.set_ylabel(ylabel); a.set_title(title); a.tick_params(width=0.8, length=2.5)

bars(ax[0,0], "dG_separated", "(A) Binding energy  dG_separated", "REU  (negative = favorable)",
     cut=0.0, cutlabel="favorable < 0")
bars(ax[0,1], "sc_value", "(B) Shape complementarity", "sc_value", cut=0.65, cutlabel="good ≥ 0.65")
bars(ax[1,0], "dSASA", "(C) Buried interface area", "dSASA (Å²)")
bars(ax[1,1], "interface_hbonds", "(D) Interface H-bonds", "count", conv=lambda s:int(float(s)))

# group legend on panel C
import matplotlib.patches as mp
ax[1,0].legend(handles=[mp.Patch(color=c, label=g) for g,c in GCOL.items()],
               frameon=False, loc="upper left", title="monomers", fontsize=7)

fig.text(0.5, 0.004, "Rosetta InterfaceAnalyzer (pack_separated) on binder threaded into the RFd docking pose. "
         "All dG positive (repulsive); contact area forms but interface is energetically unfavorable. "
         "Absolute dG inflated by un-relaxed threading; sign and ranking are the message.",
         ha="center", fontsize=7, color="black")
fig.tight_layout(rect=[0,0.012,1,0.975])
for ext in ("pdf","png","svg"): fig.savefig(f"{ROOT}/reports/bucketC_interface.{ext}", dpi=300)
print("wrote reports/bucketC_interface.{pdf,png,svg}")
