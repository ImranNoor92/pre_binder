#!/usr/bin/env python3
"""Compact re-plot of the target profile (reads outputs/target_profile.csv; no recompute).
Smaller canvas, publication contrasting palette. Keeps font sizes. Output: reports/target_profile_v2.{pdf,png,svg}
"""
from __future__ import annotations
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/data/binder_software/pre-binder"
HOTSPOT = 111
LOOPS = [(47, 51), (105, 120), (276, 296), (332, 335)]

# SASA + hydropathy come from the precomputed CSV; the model's "B-factor" column is actually
# pLDDT (mean 94.7, neg. corr with SASA) -> use (100 - pLDDT) as a predicted-flexibility proxy.
x, sasa, kd = [], [], []
for r in csv.DictReader(open(f"{ROOT}/outputs/target_profile.csv")):
    x.append(int(r["resid"])); sasa.append(float(r["rel_sasa"])); kd.append(float(r["kd_hydropathy_w7"]))
x = np.array(x); sasa = np.array(sasa); kd = np.array(kd)

# read pLDDT per residue (mean over 6 chains) directly from the model
import collections
plddt_by_res = collections.defaultdict(list)
for l in open(f"{ROOT}/inputs/151lp3t3_hexamer_6chain.pdb"):
    if l.startswith("ATOM") and l[12:16].strip() == "CA":
        plddt_by_res[int(l[22:26])].append(float(l[60:66]))
plddt = np.array([np.mean(plddt_by_res[r]) for r in x])
flex = 100.0 - plddt                      # low confidence ~ flexible/disordered

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
    "font.size":8,"axes.titlesize":10,"axes.labelsize":9,"xtick.labelsize":7,"ytick.labelsize":7,
    "legend.fontsize":7,"text.color":"black","axes.labelcolor":"black","axes.titlecolor":"black",
    "xtick.color":"black","ytick.color":"black","axes.edgecolor":"black",
    "pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})

LOOPC="#f3df9f"; HOT="#b3001b"          # gold shading, strong red hotspot
C_RMSF="#1d3557"; C_SASA="#0b7a6f"      # navy, teal
C_PHOB="#e07b39"; C_PHIL="#3b6ea5"      # orange (hydrophobic) vs blue (hydrophilic)

fig, ax = plt.subplots(3, 1, figsize=(6.9, 6.4), sharex=True)
fig.suptitle("AAV T=3 hexamer — per-residue target profile", fontsize=11, fontweight="bold", y=0.995)

def shade(a, lab=False):
    for j,(s0,e0) in enumerate(LOOPS):
        a.axvspan(s0, e0, color=LOOPC, alpha=0.85, lw=0,
                  label=("exposed loop" if (lab and j==0) else None))
    a.axvline(HOTSPOT, color=HOT, ls="--", lw=1.3, label=("Asp111 hotspot" if lab else None))
    a.tick_params(width=0.8, length=2.5)
    a.margins(x=0)

shade(ax[0], True)
ax[0].fill_between(x, 0, flex, color=C_RMSF, alpha=0.18, lw=0)
ax[0].plot(x, flex, color=C_RMSF, lw=1.0)
ax[0].set_ylabel("100 − pLDDT"); ax[0].set_ylim(0, max(flex.max()*1.1, 10))
ax[0].set_title("(A) Predicted flexibility (model confidence; higher = more flexible)")
ax[0].legend(frameon=False, loc="upper right", ncol=2, handlelength=1.4)

shade(ax[1])
ax[1].fill_between(x, 0, sasa, color=C_SASA, alpha=0.18, lw=0)
ax[1].plot(x, sasa, color=C_SASA, lw=1.0)
ax[1].axhline(0.25, color="#444444", ls=":", lw=0.9, label="exposed > 0.25")
ax[1].set_ylabel("Relative SASA"); ax[1].set_ylim(0,1)
ax[1].set_title("(B) Solvent exposure (assembled hexamer)")
ax[1].legend(frameon=False, loc="upper right", handlelength=1.4)

shade(ax[2])
ax[2].fill_between(x, 0, kd, where=(kd>=0), color=C_PHOB, alpha=0.85, lw=0, label="hydrophobic")
ax[2].fill_between(x, 0, kd, where=(kd<0), color=C_PHIL, alpha=0.85, lw=0, label="hydrophilic")
ax[2].axhline(0, color="black", lw=0.6)
ax[2].set_ylabel("Kyte–Doolittle"); ax[2].set_xlabel("VP3 residue number")
ax[2].set_title("(C) Hydrophobicity (window = 7)")
ax[2].legend(frameon=False, loc="upper right", ncol=2, handlelength=1.4)

fig.tight_layout(rect=[0,0,1,0.975])
for ext in ("pdf","png","svg"): fig.savefig(f"{ROOT}/reports/target_profile_v2.{ext}", dpi=300)
print("wrote reports/target_profile_v2.{pdf,png,svg}")
