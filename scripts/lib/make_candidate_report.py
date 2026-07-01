#!/usr/bin/env python3
"""Consolidated report of the AF2-multimer-evaluated candidates (top-5-per-group, n=20).
Merges binder-only foldability (Bucket A: pLDDT, RMSD) with interface confidence
(Bucket B: i-pTM, i-pAE in Å, ipSAE, pDockQ2). Ranked by i-pTM (best interface first).
Outputs outputs/top_candidates_report.csv + reports/top_candidates_table.{pdf,png,svg}.
All candidates are below binding thresholds (reported for completeness).
"""
from __future__ import annotations
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/data/binder_software/pre-binder"

# Bucket A foldability per tag
fold = {}
for r in csv.DictReader(open(f"{ROOT}/outputs/14_af2_binderonly/ranked_binderonly.csv")):
    fold[r["tag"]] = (r["plddt_binder"], r["binder_rmsd"])
# binder length per tag
blen = {}
for line in open(f"{ROOT}/outputs/15_colabdesign/results_bucketB.jsonl"):
    d = json.loads(line); blen[d["tag"]] = d.get("binder_len", "")

rows = []
for r in csv.DictReader(open(f"{ROOT}/outputs/bucketB_metrics.csv")):
    tag = r["tag"]; pl, rm = fold.get(tag, ("", ""))
    rows.append({
        "design": tag, "group": r["group"], "monomers": f"{r['n_monomers']}/6",
        "binder_len": blen.get(tag, ""),
        "plddt": round(float(pl), 1) if pl else "",
        "rmsd_A": round(float(rm), 2) if rm else "",
        "i_ptm": round(float(r["i_ptm"]), 3),
        "i_pae_A": round(float(r["i_pae"]) * 31.0, 1),   # ColabDesign i_pae is /31 -> back to Å
        "ipSAE": round(float(r["ipSAE"]), 3),
        "pDockQ2": round(float(r["pDockQ2"]), 3),
    })
rows.sort(key=lambda x: -x["i_ptm"])
for i, r in enumerate(rows, 1): r["rank"] = i

cols = ["rank", "design", "monomers", "binder_len", "plddt", "rmsd_A", "i_ptm", "i_pae_A", "ipSAE", "pDockQ2"]
with open(f"{ROOT}/outputs/top_candidates_report.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n"); w.writeheader()
    for r in rows: w.writerow({k: r[k] for k in cols})
print(f"wrote outputs/top_candidates_report.csv ({len(rows)} candidates)")

# ---- table figure ----
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "text.color": "black", "pdf.fonttype": 42, "svg.fonttype": "none"})
HEAD = ["Rank", "Design", "Mono.", "Len", "pLDDT", "RMSD\n(Å)", "i-pTM", "i-pAE\n(Å)", "ipSAE", "pDockQ2"]
cell = [[r["rank"], r["design"].replace("design_", "d"), r["monomers"], r["binder_len"],
         r["plddt"], r["rmsd_A"], r["i_ptm"], r["i_pae_A"], r["ipSAE"], r["pDockQ2"]] for r in rows]

fig, ax = plt.subplots(figsize=(7.6, 7.4)); ax.axis("off")
ax.set_title("AF2-multimer-evaluated candidates (top 5 per monomer-group, n = 20)",
             fontsize=11, fontweight="bold", pad=14)
t = ax.table(cellText=cell, colLabels=HEAD, loc="center", cellLoc="center")
t.auto_set_font_size(False); t.set_fontsize(8); t.scale(1, 1.35)
GCOL = {"6/6": "#d7e6f2", "5/6": "#d8efe5", "4/6": "#fcefd6", "3/6": "#f3e1ec"}
ncol = len(HEAD)
for (rr, cc), cellobj in t.get_celld().items():
    cellobj.set_edgecolor("#888888"); cellobj.set_linewidth(0.5)
    if rr == 0:
        cellobj.set_facecolor("#33414f"); cellobj.set_text_props(color="white", fontweight="bold", fontsize=8)
        cellobj.set_height(0.075)
    else:
        cellobj.set_facecolor(GCOL.get(rows[rr-1]["monomers"], "#ffffff"))
ax.text(0.5, -0.02,
        "Foldability (pLDDT, RMSD) from AF2 binder-only; interface metrics from ColabDesign/AF2-multimer "
        "(hexamer templated). Pass bars: i-pTM ≥ 0.5, i-pAE ≤ ~10 Å, ipSAE ≥ 0.5, pDockQ2 ≥ 0.23.\n"
        "All 20 candidates fall below every binding threshold; row shading = monomer group. Listed for completeness.",
        transform=ax.transAxes, ha="center", va="top", fontsize=7, color="black")
fig.tight_layout()
for ext in ("pdf", "png", "svg"): fig.savefig(f"{ROOT}/reports/top_candidates_table.{ext}", dpi=300, bbox_inches="tight")
print("wrote reports/top_candidates_table.{pdf,png,svg}")
