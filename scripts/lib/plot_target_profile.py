#!/usr/bin/env python3
"""Target characterization (3 stacked panels vs residue) to justify the hotspot critique:
  (A) B-factor-derived RMSF  — residue mobility (from deposited B-factors; no MD available)
  (B) relative SASA          — exposure in the ASSEMBLED hexamer (Shrake-Rupley, mean of 6 chains)
  (C) Kyte-Doolittle hydropathy (window-smoothed) — hydrophobicity
Exposed loop regions shaded; Asp111 hotspot marked. A4 portrait, Liberation Sans, pure-black.
Output: reports/target_profile.{pdf,png,svg}  +  outputs/target_profile.csv
"""
from __future__ import annotations
import csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = "/data/binder_software/pre-binder"
PDB = f"{ROOT}/inputs/151lp3t3_hexamer_6chain.pdb"
CHAINS = "ABCDEF"
HOTSPOT = 111
LOOPS = [(47, 51), (105, 120), (276, 296), (332, 335)]   # exposed top loops to highlight

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
      'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
MAXASA = {'A':129.,'R':274.,'N':195.,'D':193.,'C':167.,'E':223.,'Q':225.,'G':104.,'H':224.,'I':197.,
 'L':201.,'K':236.,'M':224.,'F':240.,'P':159.,'S':155.,'T':172.,'W':285.,'Y':263.,'V':174.}

def main():
    s = PDBParser(QUIET=True).get_structure("hex", PDB)[0]
    print("computing SASA on assembled hexamer (Shrake-Rupley)...")
    ShrakeRupley().compute(s, level="R")

    # per-residue, averaged over the 6 chains
    resids = sorted({r.id[1] for r in s["A"] if r.id[0] == " "})
    bf, rsasa, seq = {}, {}, {}
    for rid in resids:
        bvals, svals, aa = [], [], None
        for c in CHAINS:
            try: res = s[c][(" ", rid, " ")]
            except KeyError: continue
            rn = res.resname
            if rn not in THREE2ONE: continue
            aa = THREE2ONE[rn]
            ca = res["CA"] if "CA" in res else None
            if ca is not None: bvals.append(ca.get_bfactor())
            svals.append(res.sasa / MAXASA[aa])
        if aa is None: continue
        bf[rid] = np.mean(bvals) if bvals else np.nan
        rsasa[rid] = np.clip(np.mean(svals), 0, 1) if svals else np.nan
        seq[rid] = aa

    x = np.array(sorted(seq))
    bfa = np.array([bf[i] for i in x])
    rmsf = np.sqrt(3.0 * bfa / (8.0 * math.pi**2))          # B = 8pi^2/3 <u^2>
    sasa = np.array([rsasa[i] for i in x])
    kd_raw = np.array([KD[seq[i]] for i in x])
    w = 7; kern = np.ones(w)/w
    kd = np.convolve(kd_raw, kern, mode="same")

    with open(f"{ROOT}/outputs/target_profile.csv", "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["resid","aa","rmsf_from_bfactor","rel_sasa","kd_hydropathy_w7"])
        for i, r in enumerate(x): wr.writerow([r, seq[r], round(float(rmsf[i]),3), round(float(sasa[i]),3), round(float(kd[i]),3)])
    print(f"wrote outputs/target_profile.csv ({len(x)} residues)")

    # ---- figure ----
    plt.rcParams.update({
        "font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
        "font.size":8,"axes.titlesize":10,"axes.labelsize":9,"xtick.labelsize":7,"ytick.labelsize":7,
        "legend.fontsize":7,"text.color":"black","axes.labelcolor":"black","axes.titlecolor":"black",
        "xtick.color":"black","ytick.color":"black","axes.edgecolor":"black",
        "pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
    LINE="#33414f"; LOOPC="#cdd7df"; HOT="#8b0000"
    fig, ax = plt.subplots(3, 1, figsize=(8.27, 11.69), sharex=True)
    fig.suptitle("AAV T=3 hexamer — per-residue target profile (one VP3 protomer)",
                 fontsize=11, fontweight="bold", y=0.985)

    def shade(a, label_loops=False):
        for j,(s0,e0) in enumerate(LOOPS):
            a.axvspan(s0, e0, color=LOOPC, alpha=0.7, lw=0,
                      label=("exposed loop" if (label_loops and j==0) else None))
        a.axvline(HOTSPOT, color=HOT, ls="--", lw=1.2,
                  label=("hotspot Asp111" if label_loops else None))
        a.tick_params(width=0.8, length=3)

    ax[0].plot(x, rmsf, color=LINE, lw=0.8); shade(ax[0], True)
    ax[0].set_ylabel("RMSF (Å)"); ax[0].set_title("(A) Residue mobility — from deposited B-factors")
    ax[0].legend(frameon=False, loc="upper right", ncol=2)

    ax[1].plot(x, sasa, color=LINE, lw=0.8); shade(ax[1])
    ax[1].axhline(0.25, color="#555555", ls=":", lw=0.9, label="exposed > 0.25")
    ax[1].set_ylabel("Relative SASA"); ax[1].set_ylim(0,1)
    ax[1].set_title("(B) Solvent exposure in assembled hexamer"); ax[1].legend(frameon=False, loc="upper right")

    ax[2].fill_between(x, 0, kd, where=(kd>=0), color="#b0b7be", lw=0, label="hydrophobic")
    ax[2].fill_between(x, 0, kd, where=(kd<0), color="#dfe4e8", lw=0, label="hydrophilic")
    ax[2].plot(x, kd, color=LINE, lw=0.7); shade(ax[2])
    ax[2].axhline(0, color="black", lw=0.6)
    ax[2].set_ylabel("Kyte–Doolittle (w=7)"); ax[2].set_xlabel("VP3 residue number")
    ax[2].set_title("(C) Hydrophobicity"); ax[2].legend(frameon=False, loc="upper right", ncol=2)
    ax[2].set_xlim(x.min(), x.max())

    fig.text(0.5,0.005,"Shaded = exposed top loops; dashed red = Asp111 hotspot. RMSF from B (RMSF=√(3B/8π²)); "
             "relative SASA = Shrake-Rupley / max ASA (Tien 2013), mean of 6 chains.",
             ha="center", fontsize=7, color="black")
    fig.tight_layout(rect=[0,0.012,1,0.975])
    for ext in ("pdf","png","svg"): fig.savefig(f"{ROOT}/reports/target_profile.{ext}", dpi=300)
    print("wrote reports/target_profile.{pdf,png,svg}")
    # quick textless check at hotspot
    i = list(x).index(HOTSPOT)
    print(f"  at Asp111: rel_SASA={sasa[i]:.2f}, KD(w7)={kd[i]:.2f}, RMSF={rmsf[i]:.2f}")

if __name__ == "__main__":
    main()
