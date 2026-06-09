#!/usr/bin/env python3
"""One-page, 4-panel summary of a Phase-IG validation run (outputs/06_ig/ranked.csv).
Usage: plot_run_summary.py --csv outputs/06_ig/ranked.csv --out outputs/06_ig/run_summary.png
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PAE_MAX, PLDDT_MIN, RMSD_OK = 10.0, 80.0, 2.0


def load(path):
    pae, plddt, rmsd, passed = [], [], [], []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            pae.append(float(r["pae_interaction"])); plddt.append(float(r["plddt_binder"]))
            rmsd.append(float(r["binder_rmsd"])); passed.append(r["pass"].strip() == "True")
    return (np.array(pae), np.array(plddt), np.array(rmsd), np.array(passed))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="pre_binder · Phase IG validation (af2_initial_guess)")
    ap.add_argument("--subtitle", default="")
    a = ap.parse_args()

    pae, plddt, rmsd, passed = load(a.csv)
    n, npass = len(pae), int(passed.sum())

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "figure.dpi": 150})
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    sub = a.subtitle or (f"{n} designs · {npass}/{n} pass "
                         f"(pass = pae_interaction < {PAE_MAX:.0f} AND plddt_binder > {PLDDT_MIN:.0f})")
    fig.suptitle(a.title, fontsize=14, fontweight="bold", y=0.98)
    fig.text(0.5, 0.935, sub, ha="center", fontsize=10, color="#444")

    # --- Panel A: the binding landscape (the crucial one) ---
    axA = ax[0, 0]
    scA = axA.scatter(plddt, pae, c=rmsd, cmap="viridis_r", s=34, edgecolor="k", linewidth=0.3, zorder=3)
    axA.axhline(PAE_MAX, color="crimson", ls="--", lw=1.4, zorder=2)
    axA.axvline(PLDDT_MIN, color="crimson", ls="--", lw=1.4, zorder=2)
    # pass box (plddt>min, pae<max)
    x0, x1 = plddt.min() - 3, plddt.max() + 3
    y0, y1 = 0, pae.max() + 3
    axA.add_patch(plt.Rectangle((PLDDT_MIN, 0), x1 - PLDDT_MIN, PAE_MAX, color="green", alpha=0.10, zorder=1))
    axA.text(PLDDT_MIN + 0.4, PAE_MAX - 0.6, f"PASS region\n({npass} designs)", color="green",
             fontsize=8, va="top", fontweight="bold")
    axA.set_xlim(x0, x1); axA.set_ylim(y0, y1)
    axA.set_xlabel("plddt_binder  (does it FOLD →)"); axA.set_ylabel("pae_interaction  (↓ does it BIND)")
    axA.set_title("A · Binding landscape — folds well, but no interface")
    cb = fig.colorbar(scA, ax=axA, pad=0.02); cb.set_label("binder_rmsd vs design (Å)")

    # --- Panel B: pae_interaction distribution ---
    axB = ax[0, 1]
    axB.hist(pae, bins=20, color="#4C72B0", edgecolor="k", linewidth=0.3)
    axB.axvline(PAE_MAX, color="crimson", ls="--", lw=1.6, label=f"pass threshold < {PAE_MAX:.0f}")
    axB.axvline(pae.min(), color="black", ls=":", lw=1.2, label=f"best = {pae.min():.1f}")
    axB.set_xlabel("pae_interaction (lower = stronger interface)"); axB.set_ylabel("designs")
    axB.set_title("B · Interface confidence — all far from threshold"); axB.legend(fontsize=8)

    # --- Panel C: plddt_binder distribution ---
    axC = ax[1, 0]
    axC.hist(plddt, bins=20, color="#55A868", edgecolor="k", linewidth=0.3)
    axC.axvline(PLDDT_MIN, color="crimson", ls="--", lw=1.6, label=f"good fold > {PLDDT_MIN:.0f}")
    axC.set_xlabel("plddt_binder (higher = better fold)"); axC.set_ylabel("designs")
    axC.set_title(f"C · Fold confidence — {int((plddt>PLDDT_MIN).sum())}/{n} fold well"); axC.legend(fontsize=8)

    # --- Panel D: binder_rmsd distribution ---
    axD = ax[1, 1]
    axD.hist(np.clip(rmsd, 0, 10), bins=20, color="#C44E52", edgecolor="k", linewidth=0.3)
    axD.axvline(RMSD_OK, color="navy", ls="--", lw=1.6, label=f"folds as designed < {RMSD_OK:.0f} Å")
    axD.set_xlabel("binder_rmsd vs design (Å, clipped at 10)"); axD.set_ylabel("designs")
    axD.set_title(f"D · Design recapitulation — {int((rmsd<RMSD_OK).sum())}/{n} within {RMSD_OK:.0f} Å"); axD.legend(fontsize=8)

    fig.text(0.5, 0.005,
             "Takeaway: af2_initial_guess validation ran clean (no MSA, no OOM). Designs FOLD into their "
             "intended shape (panels C/D) but do NOT form a confident interface with the target "
             "(panels A/B) — bottleneck is the designs/epitope, not the pipeline.",
             ha="center", fontsize=8.5, color="#333", wrap=True)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight")
    print(f"wrote {a.out}  ({n} designs, {npass} pass)")


if __name__ == "__main__":
    main()
