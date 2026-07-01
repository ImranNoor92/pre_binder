#!/usr/bin/env python3
"""Bucket A (binder-intrinsic, target-independent) QC for ALL central-run binders.

Sources: outputs/14_af2_binderonly/scores.sc (AF2 binder-only: pLDDT, intra-binder PAE, RMSD)
         outputs/13_mpnn_central/seqs/*.fa (the MPNN sequences -> developability via ProtParam)
Output : outputs/bucketA_metrics.csv  +  reports/bucketA_binder_intrinsic.{pdf,png}

Figure: 3x3 panels, A4 portrait, Liberation Sans (Arial-metric), pure-black text, muted greyscale.
"""
from __future__ import annotations
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = "/data/binder_software/pre-binder"
SC = f"{ROOT}/outputs/14_af2_binderonly/scores.sc"
SEQDIR = f"{ROOT}/outputs/13_mpnn_central/seqs"

# ---------- gather ----------
def parse_scores():
    rows, hdr = {}, None
    for line in open(SC):
        if not line.startswith("SCORE:"):
            continue
        t = line.split()[1:]
        if t[0] == "binder_aligned_rmsd":
            hdr = t; continue
        if hdr and len(t) == len(hdr):
            d = dict(zip(hdr, t)); tag = re.sub(r"_af2pred$", "", d["description"])
            rows[tag] = d
    return rows

def seqs_by_tag():
    out = {}
    for fa in glob.glob(f"{SEQDIR}/design_*.fa"):
        name = os.path.basename(fa)[:-3]
        recs, cur = [], None
        for line in open(fa):
            if line.startswith(">"): cur = []; recs.append(cur)
            elif cur is not None: cur.append(line.strip())
        seqs = ["".join(r) for r in recs][1:]   # skip native poly-G
        for i, s in enumerate(seqs):
            out[f"{name}_s{i}"] = s
    return out

def aliphatic_index(seq):
    n = len(seq)
    f = lambda a: 100.0 * seq.count(a) / n
    return f("A") + 2.9 * f("V") + 3.9 * (f("I") + f("L"))

def main():
    sc = parse_scores(); sq = seqs_by_tag()
    rows = []
    for tag, d in sc.items():
        s = sq.get(tag)
        if not s or set(s) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        pa = ProteinAnalysis(s)
        try: charge = pa.charge_at_pH(7.0)
        except Exception: charge = float("nan")
        rows.append({
            "tag": tag,
            "plddt": float(d["plddt_binder"]),
            "pae_binder": float(d["pae_binder"]),
            "rmsd": float(d["binder_aligned_rmsd"]),
            "length": len(s),
            "instability": pa.instability_index(),
            "aliphatic": aliphatic_index(s),
            "gravy": pa.gravy(),
            "net_charge_pH7": charge,
            "pI": pa.isoelectric_point(),
        })
    with open(f"{ROOT}/outputs/bucketA_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"computed Bucket A for {len(rows)} binders -> outputs/bucketA_metrics.csv")

    # ---------- style ----------
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "text.color": "black", "axes.labelcolor": "black", "axes.titlecolor": "black",
        "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8,
    })
    BAR = "#9aa7b3"; EDGE = "black"; CUT = "#8b0000"   # muted slate bars, dark-red cutoff lines
    A = {k: np.array([r[k] for r in rows], float) for k in rows[0] if k != "tag"}

    fig, ax = plt.subplots(3, 3, figsize=(8.27, 11.69))
    fig.suptitle("Bucket A — binder-intrinsic quality (central-run binders, n = %d)" % len(rows),
                 fontsize=11, fontweight="bold", y=0.985)

    def hist(a, data, xlabel, title, cut=None, cutlabel=None, bins=35):
        a.hist(data, bins=bins, color=BAR, edgecolor=EDGE, linewidth=0.4)
        if cut is not None:
            a.axvline(cut, color=CUT, ls="--", lw=1.2, label=cutlabel)
            a.legend(frameon=False, loc="best")
        a.set_xlabel(xlabel); a.set_ylabel("Number of binders"); a.set_title(title)
        a.tick_params(width=0.8, length=3)

    # A pLDDT
    hist(ax[0,0], A["plddt"], "Binder pLDDT", "(A) Fold confidence",
         cut=80, cutlabel="good fold ≥ 80")
    # B RMSD
    hist(ax[0,1], np.clip(A["rmsd"],0,10), "RMSD to design (Å)", "(B) Self-consistency",
         cut=2.0, cutlabel="pass < 2 Å")
    # C scatter pLDDT vs RMSD
    a = ax[0,2]
    a.scatter(A["plddt"], np.clip(A["rmsd"],0,10), s=5, c="#4d4d4d", alpha=0.35,
              edgecolor="none", rasterized=True)
    a.axhline(2.0, color=CUT, ls="--", lw=1.0); a.axvline(80, color=CUT, ls="--", lw=1.0)
    npass = int(((A["plddt"]>80)&(A["rmsd"]<2)).sum())
    a.set_xlabel("Binder pLDDT"); a.set_ylabel("RMSD to design (Å)")
    a.set_title("(C) Foldability landscape"); a.tick_params(width=0.8, length=3)
    a.text(0.04,0.95,f"pass: {npass}/{len(rows)}", transform=a.transAxes, va="top", fontsize=7)
    # D intra-binder PAE
    hist(ax[1,0], A["pae_binder"], "Intra-binder PAE (Å)", "(D) Topology confidence")
    # E instability
    hist(ax[1,1], A["instability"], "Instability index", "(E) Predicted stability",
         cut=40, cutlabel="stable < 40")
    # F aliphatic
    hist(ax[1,2], A["aliphatic"], "Aliphatic index", "(F) Predicted thermostability")
    # G GRAVY
    hist(ax[2,0], A["gravy"], "GRAVY (hydropathy)", "(G) Solubility / aggregation",
         cut=0.0, cutlabel="soluble < 0")
    # H net charge
    hist(ax[2,1], A["net_charge_pH7"], "Net charge at pH 7", "(H) Charge / solubility",
         cut=0.0, cutlabel="neutral")
    # I length
    hist(ax[2,2], A["length"], "Sequence length (residues)", "(I) Binder size")

    fig.text(0.5, 0.005,
             "A/B/C/D from AlphaFold2 binder-only prediction; E–I from sequence (Biopython ProtParam). "
             "Higher better: pLDDT, aliphatic. Lower better: RMSD, PAE, instability.",
             ha="center", fontsize=7, color="black")
    fig.tight_layout(rect=[0, 0.012, 1, 0.975])
    os.makedirs(f"{ROOT}/reports", exist_ok=True)
    fig.savefig(f"{ROOT}/reports/bucketA_binder_intrinsic.pdf", dpi=300)
    fig.savefig(f"{ROOT}/reports/bucketA_binder_intrinsic.png", dpi=300)
    fig.savefig(f"{ROOT}/reports/bucketA_binder_intrinsic.svg")
    print("wrote reports/bucketA_binder_intrinsic.{pdf,png,svg}")

    # quick text summary
    def pct(m): return 100.0*np.mean(m)
    print(f"  pLDDT>80: {pct(A['plddt']>80):.0f}% | RMSD<2: {pct(A['rmsd']<2):.0f}% | "
          f"both: {pct((A['plddt']>80)&(A['rmsd']<2)):.0f}%")
    print(f"  instability<40 (stable): {pct(A['instability']<40):.0f}% | GRAVY<0 (soluble): {pct(A['gravy']<0):.0f}%")

if __name__ == "__main__":
    main()
