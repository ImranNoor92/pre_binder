#!/usr/bin/env python3
"""Compute Bucket B interface metrics from saved ColabDesign complex predictions and plot.

For each outputs/15_colabdesign/bucketB/<tag>/ (pae.npy [A, raw], plddt.npy [0-1], complex.pdb
with chain A=target, chain B=binder):
  - i_pTM, i_pAE       : from results_bucketB.jsonl (ColabDesign log)
  - ipSAE (d0res)      : Dunbrack interface score, dilution-robust  (max better, ~>0.5 good)
  - pDockQ2 (Zhu 2023) : interface quality (max better, >0.23 ~ acceptable interface)

Figure: 2x2, A4 portrait, Liberation Sans (Arial-metric), pure-black, muted greyscale.
"""
from __future__ import annotations
import csv, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/data/binder_software/pre-binder"
BDIR = f"{ROOT}/outputs/15_colabdesign/bucketB"
JSONL = f"{ROOT}/outputs/15_colabdesign/results_bucketB.jsonl"
PAE_CUT = 10.0


def d0(L):
    L = max(int(L), 1)
    return max(1.24 * (L - 15) ** (1.0 / 3.0) - 1.8, 1.0) if L > 27 else 1.0

def ptm(x, d): return 1.0 / (1.0 + (x / d) ** 2)

def ipsae_asym(pae, idx1, idx2):
    best = 0.0
    for i in idx1:
        row = pae[i, idx2]
        m = row < PAE_CUT
        n = int(m.sum())
        if n == 0: continue
        v = ptm(row[m], d0(n)).mean()
        if v > best: best = v
    return best

def parse_complex(pdb):
    """Return per-CA chain list and CB coords (CA for Gly), in file order."""
    chain, cb = [], []
    cur = {}
    order = []
    for l in open(pdb):
        if not l.startswith("ATOM"): continue
        rid = (l[21], l[22:26]); an = l[12:16].strip()
        if rid not in cur:
            cur[rid] = {}; order.append(rid)
        cur[rid][an] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        cur[rid]["_ch"] = l[21]
    for rid in order:
        a = cur[rid]
        chain.append(a["_ch"])
        cb.append(a.get("CB", a.get("CA", (np.nan,)*3)))
    return np.array(chain), np.array(cb, float)


def pdockq2(pae, plddt100, cb, chain, binder="B", target="A"):
    bi = np.where(chain == binder)[0]; ti = np.where(chain == target)[0]
    D = np.sqrt(((cb[bi][:, None] - cb[ti][None]) ** 2).sum(-1))
    bj, tj = np.where(D < 8.0)
    if len(bj) == 0: return 0.0
    bset = bi[np.unique(bj)]; tset = ti[np.unique(tj)]
    paes = np.array([pae[bi[a], ti[b]] for a, b in zip(bj, tj)])
    paes = np.concatenate([paes, [pae[ti[b], bi[a]] for a, b in zip(bj, tj)]])  # both directions
    paeterm = (1.0 / (1.0 + (paes / 10.0) ** 2)).mean()
    plddt_int = plddt100[np.concatenate([bset, tset])].mean()
    x = paeterm * plddt_int
    L, x0, k, b = 1.31034849, 84.7326239, 0.0747157696, 0.00501886443
    return L / (1.0 + np.exp(-k * (x - x0))) + b


def main():
    scal = {json.loads(l)["tag"]: json.loads(l) for l in open(JSONL)} if os.path.exists(JSONL) else {}
    rows = []
    for d in sorted(glob.glob(f"{BDIR}/design_*")):
        tag = os.path.basename(d)
        if not all(os.path.exists(f"{d}/{f}") for f in ("pae.npy", "plddt.npy", "complex.pdb")):
            continue
        pae = np.load(f"{d}/pae.npy"); plddt = np.load(f"{d}/plddt.npy")
        plddt100 = plddt * 100 if np.nanmax(plddt) <= 1.5 else plddt
        chain, cb = parse_complex(f"{d}/complex.pdb")
        bidx = np.where(chain == "B")[0]; tidx = np.where(chain == "A")[0]
        ips = max(ipsae_asym(pae, bidx, tidx), ipsae_asym(pae, tidx, bidx))
        pdq = pdockq2(pae, plddt100, cb, chain)
        s = scal.get(tag, {})
        rows.append({"tag": tag, "group": s.get("group", ""), "n_monomers": s.get("n_monomers", ""),
                     "i_ptm": s.get("i_ptm", float("nan")), "i_pae": s.get("i_pae", float("nan")),
                     "ipSAE": round(float(ips), 4), "pDockQ2": round(float(pdq), 4)})
    rows.sort(key=lambda r: -r["ipSAE"])
    with open(f"{ROOT}/outputs/bucketB_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"computed Bucket B for {len(rows)} designs -> outputs/bucketB_metrics.csv")
    for r in rows[:5]:
        print(f"  {r['tag']}: i_ptm={r['i_ptm']} ipSAE={r['ipSAE']} pDockQ2={r['pDockQ2']}")

    # ---- figure (2x2, A4 portrait) ----
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "text.color": "black", "axes.labelcolor": "black",
        "axes.titlecolor": "black", "xtick.color": "black", "ytick.color": "black",
        "axes.edgecolor": "black", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
    EDGE, CUT = "black", "#b3001b"
    # colorblind-safe qualitative palette by monomer group (Okabe-Ito)
    GCOL = {"6of6": "#0072B2", "5of6": "#009E73", "4of6": "#E69F00", "3of6": "#CC79A7"}
    cols = [GCOL.get(r["group"], "#777777") for r in rows]
    tags = [r["tag"].replace("design_", "d") for r in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(2, 2, figsize=(7.3, 6.6))
    fig.suptitle("Bucket B — interface confidence (top-5/group finalists, n = %d)" % len(rows),
                 fontsize=11, fontweight="bold", y=0.995)

    def barpanel(a, key, title, ylabel, cut, cutlabel, ymax=None):
        vals = np.array([r[key] for r in rows], float)
        a.bar(x, vals, color=cols, edgecolor=EDGE, linewidth=0.3)
        a.axhline(cut, color=CUT, ls="--", lw=1.3, label=cutlabel)
        a.set_xticks(x); a.set_xticklabels(tags, rotation=90, fontsize=7)
        a.set_ylabel(ylabel); a.set_title(title)
        if ymax: a.set_ylim(0, ymax)
        a.legend(frameon=False, loc="upper right"); a.tick_params(width=0.8, length=2.5)

    barpanel(ax[0,0], "i_ptm", "(A) i-pTM (dilution-prone)", "i-pTM", 0.5, "pass ≥ 0.5", 0.6)
    barpanel(ax[0,1], "ipSAE", "(B) ipSAE (dilution-robust)", "ipSAE", 0.5, "good ≥ 0.5", 0.6)
    barpanel(ax[1,0], "pDockQ2", "(C) pDockQ2 (interface quality)", "pDockQ2", 0.23, "acceptable ≥ 0.23", 0.3)
    # D: agreement scatter i_ptm vs ipSAE, colored by group
    a = ax[1,1]
    for g, gc in GCOL.items():
        gr = [r for r in rows if r["group"] == g]
        if gr: a.scatter([r["i_ptm"] for r in gr], [r["ipSAE"] for r in gr], s=26, c=gc,
                          edgecolor="black", linewidth=0.3, label=g)
    a.plot([0,0.6],[0,0.6], color="#999999", ls=":", lw=0.8)
    a.axhline(0.5, color=CUT, ls="--", lw=1.0); a.axvline(0.5, color=CUT, ls="--", lw=1.0)
    a.set_xlim(0,0.6); a.set_ylim(0,0.6)
    a.set_xlabel("i-pTM"); a.set_ylabel("ipSAE")
    a.set_title("(D) Agreement (by group)"); a.legend(frameon=False, loc="upper left", title="monomers")
    a.tick_params(width=0.8, length=2.5)

    fig.text(0.5, 0.004, "ColabDesign/AF2-multimer, hexamer templated (held). ipSAE/pDockQ2 interface-specific "
             "(not diluted). Higher = better; all designs fall far below every threshold.",
             ha="center", fontsize=7, color="black")
    fig.tight_layout(rect=[0, 0.01, 1, 0.975])
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{ROOT}/reports/bucketB_interface.{ext}", dpi=300)
    print("wrote reports/bucketB_interface.{pdf,png,svg}")


if __name__ == "__main__":
    main()
