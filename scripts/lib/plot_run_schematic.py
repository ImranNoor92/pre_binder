#!/usr/bin/env python3
"""Methods/workflow schematic for the central single-binder run: conditions, parameters,
assumptions only (no results). A4 portrait, Liberation Sans (Arial-metric), pure-black text,
muted greyscale. Outputs reports/run_schematic.{pdf,png,svg}.
"""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "text.color": "black", "pdf.fonttype": 42, "svg.fonttype": "none"})

FILL, EDGE, ACC = "#eef2f6", "#33414f", "#cfd8e0"
fig, ax = plt.subplots(figsize=(8.27, 11.69))
ax.set_xlim(0, 100); ax.set_ylim(0, 130); ax.axis("off")

ax.text(50, 127, "Central single-binder design — workflow and parameters",
        ha="center", va="top", fontsize=11, fontweight="bold")
ax.text(50, 123.3, "AAV T=3 capsid hexamer · binder engaging all six VP3 protomers from the exterior",
        ha="center", va="top", fontsize=8.5, color="#333333")

def box(y, h, title, lines, fill=FILL, w=70, x=15, tsize=9.5, lsize=8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2",
                 linewidth=1.0, edgecolor=EDGE, facecolor=fill))
    ax.text(x + 2.5, y + h - 2.2, title, ha="left", va="top", fontsize=tsize, fontweight="bold")
    ax.text(x + 2.5, y + h - 6.2, "\n".join(lines), ha="left", va="top", fontsize=lsize, color="black")

def arrow(y0, y1, x=50):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>", mutation_scale=14,
                 linewidth=1.2, color=EDGE))

# Target header
box(110, 11, "Target  ·  full VP3 T=3 capsid hexamer (151lp3t3_hexamer.pdb)", [
    "6 chains A–F, 504 aa each (3024 aa); flattened from 6-MODEL to one frame",
    "Hotspots: Asp111 of all six chains (exposed central convergence, ~22 Å cluster)",
    "Reduced to ≤50 Å of central axis (~1158 aa) for backbone generation; full hexamer kept for validation",
], fill=ACC, w=86, x=7, tsize=9.5, lsize=7.8)
arrow(110, 104)

box(92, 12, "1 · RFdiffusion  (backbone generation)", [
    "Checkpoint: Complex_base   ·   binder length 80–120 aa",
    "noise_scale_ca = 0, noise_scale_frame = 0",
    "200 backbones · 2-GPU split · multi-crop receptor contig (R50 target)",
])
arrow(92, 86)

box(74, 12, "2 · ProteinMPNN  (sequence design)", [
    "8 sequences per backbone  →  1600 sequences",
    "Binder chain G designed; target chains A–F held fixed",
    "Sampling temperature 0.1 · seed 37",
])
arrow(74, 68)

box(56, 12, "3 · AlphaFold2 binder-only  (foldability filter)", [
    "af2_initial_guess, monomer mode · single-sequence, no MSA",
    "model_1_ptm · 3 recycles · binder predicted alone (no target)",
    "Screens fold confidence (pLDDT) and self-consistency (Cα RMSD to design)",
])
arrow(56, 50)

box(38, 12, "4 · ColabDesign / AF2-multimer  (interface validation)", [
    "multimer_v3 weights · 3 recycles · initial-guess enabled",
    "Six-chain target supplied as a held template (assembly cannot drift)",
    "Metrics: i-pTM, i-pAE, ipSAE, pDockQ2 (Buckets B/C)",
])

# Assumptions panel
box(4, 28, "Assumptions", [
    "(i)   Target conformation is fixed — the binder must not perturb it.",
    "(ii)  Binding occurs from the capsid exterior (top), not through interior loops.",
    "(iii) Specificity arises by engaging multiple protomers at a site present only in the assembled hexamer.",
    "(iv)  The central Asp111 convergence is an accessible all-six-monomer epitope.",
    "(v)   Binder foldability is assessed before (and independently of) the binding evaluation.",
    "(vi)  Target reduction for backbone generation changes compute only, not geometry;",
    "        the full hexamer is the structure of record for validation.",
], fill="#f4f1ea", w=86, x=7, tsize=9.5, lsize=7.8)

fig.savefig("/data/binder_software/pre-binder/reports/run_schematic.pdf", dpi=300, bbox_inches="tight")
fig.savefig("/data/binder_software/pre-binder/reports/run_schematic.png", dpi=300, bbox_inches="tight")
fig.savefig("/data/binder_software/pre-binder/reports/run_schematic.svg", bbox_inches="tight")
print("wrote reports/run_schematic.{pdf,png,svg}")
