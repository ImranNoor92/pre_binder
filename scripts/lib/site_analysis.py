#!/usr/bin/env python3
"""Step 1 — target site / electrostatics / shape mapping for hexamer-specific binder design.

For one VP3 protomer (chain A), evaluated in the assembled hexamer, compute per residue:
  1. exposure        : relative SASA (Shrake-Rupley / max ASA, Tien 2013), in the hexamer
  2. assembly_burial : rSASA(isolated chain A) - rSASA(hexamer); high = inter-protomer seam
  3. electrostatic   : screened-Coulomb potential phi at the side-chain charge centre
  4. concavity       : local heavy-atom density within 10 A (high among exposed = groove)
  5. n_chains        : distinct chains with atoms within 12 A (>=2 = multi-protomer seam)
Candidate spots = exposed AND (groove OR seam) AND charge-coherent; clustered into sites.

Outputs: outputs/site_analysis.csv, outputs/site_candidates.csv, reports/site_map.{pdf,png,svg},
and chainA PDBs with B-factor = phi and = concavity for PyMOL surface colouring.
Reproduce: Biopython + scipy only. (Gold standard for #3: pdb2pqr -> APBS Poisson-Boltzmann.)
"""
from __future__ import annotations
import csv, numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = "/data/binder_software/pre-binder"
PDB = f"{ROOT}/inputs/151lp3t3_hexamer_6chain.pdb"
CHAINS = "ABCDEF"; PROT = "A"
HOTSPOT = 111; OLD_HOTSPOTS = (105, 120)
KCOUL, EPS, LAMBDA_D = 332.0, 78.0, 10.0     # kcal*A/mol/e^2, water dielectric, Debye length (A)
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
MAXASA = {'A':129.,'R':274.,'N':195.,'D':193.,'C':167.,'E':223.,'Q':225.,'G':104.,'H':224.,'I':197.,
 'L':201.,'K':236.,'M':224.,'F':240.,'P':159.,'S':155.,'T':172.,'W':285.,'Y':263.,'V':174.}

def charged_groups(model):
    """(coord, q) for ionisable groups across the whole hexamer (physiological)."""
    g = []
    for ch in model:
        for res in ch:
            rn = res.resname
            try:
                if rn == "ASP": g.append((res["CG"].coord, -1.0))
                elif rn == "GLU": g.append((res["CD"].coord, -1.0))
                elif rn == "LYS": g.append((res["NZ"].coord, +1.0))
                elif rn == "ARG": g.append((res["CZ"].coord, +1.0))
                elif rn == "HIS": g.append((res["NE2"].coord, +0.1))
            except KeyError:
                pass
    return np.array([c for c, _ in g]), np.array([q for _, q in g])

def main():
    s = PDBParser(QUIET=True).get_structure("h", PDB)[0]
    print("SASA (hexamer)..."); ShrakeRupley().compute(s, level="R")
    rsasa_hex = {r.id[1]: r.sasa / MAXASA[THREE2ONE[r.resname]]
                 for r in s[PROT] if r.resname in THREE2ONE}
    # isolated protomer SASA
    iso = PDBParser(QUIET=True).get_structure("i", PDB)[0]
    for ch in list(iso):
        if ch.id != PROT: iso.detach_child(ch.id)
    ShrakeRupley().compute(iso, level="R")
    rsasa_iso = {r.id[1]: r.sasa / MAXASA[THREE2ONE[r.resname]]
                 for r in iso[PROT] if r.resname in THREE2ONE}

    # all heavy atoms + per-chain trees for concavity & multi-chain
    allxyz, allchain = [], []
    for ch in s:
        for res in ch:
            for a in res:
                if a.element != "H": allxyz.append(a.coord); allchain.append(ch.id)
    allxyz = np.array(allxyz); allchain = np.array(allchain); tree = cKDTree(allxyz)
    cg_xyz, cg_q = charged_groups(s)

    rows = []
    for res in s[PROT]:
        if res.resname not in THREE2ONE: continue
        rid = res.id[1]; aa = THREE2ONE[res.resname]
        ca = res["CA"].coord if "CA" in res else next(iter(res)).coord
        # side-chain charge centre ~ use CB (or CA for Gly) as the probe point
        probe = res["CB"].coord if "CB" in res else ca
        # 3 electrostatics
        d = np.linalg.norm(cg_xyz - probe, axis=1); d = np.clip(d, 1.5, None)
        phi = float(np.sum(KCOUL * cg_q * np.exp(-d / LAMBDA_D) / (EPS * d)))
        # 4 concavity = heavy-atom density within 10 A
        dens = len(tree.query_ball_point(ca, 10.0))
        # 5 multi-chain within 12 A
        idx = tree.query_ball_point(ca, 12.0)
        nch = len(set(allchain[idx]) )
        rows.append({"resid": rid, "aa": aa,
                     "rsasa": round(rsasa_hex.get(rid, np.nan), 3),
                     "assembly_burial": round(rsasa_iso.get(rid, 0) - rsasa_hex.get(rid, 0), 3),
                     "phi": round(phi, 3), "concavity": dens, "n_chains": nch})

    R = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != "aa"}
    exposed = R["rsasa"] > 0.25
    dens_thr = np.nanpercentile(R["concavity"][exposed], 66) if exposed.any() else 1e9
    groove = R["concavity"] >= dens_thr
    seam = R["n_chains"] >= 2
    charged = np.abs(R["phi"]) > np.nanpercentile(np.abs(R["phi"]), 60)
    candidate = exposed & (groove | seam) & charged
    for i, r in enumerate(rows): r["candidate"] = bool(candidate[i])

    with open(f"{ROOT}/outputs/site_analysis.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    # cluster candidate residues into sites by sequence contiguity (gap<=3)
    cand_ids = [rows[i]["resid"] for i in range(len(rows)) if candidate[i]]
    sites, cur = [], []
    for rid in cand_ids:
        if cur and rid - cur[-1] > 3: sites.append(cur); cur = []
        cur.append(rid)
    if cur: sites.append(cur)
    site_rows = []
    for st in sites:
        m = np.isin(R["resid"], st)
        site_rows.append({"residues": f"{st[0]}-{st[-1]}", "n_res": len(st),
                          "mean_phi": round(float(R["phi"][m].mean()), 2),
                          "charge": "acidic(-)" if R["phi"][m].mean() < 0 else "basic(+)",
                          "mean_rsasa": round(float(R["rsasa"][m].mean()), 2),
                          "mean_nchains": round(float(R["n_chains"][m].mean()), 1)})
    site_rows.sort(key=lambda x: (-x["n_res"], -abs(x["mean_phi"])))
    with open(f"{ROOT}/outputs/site_candidates.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(site_rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(site_rows)
    print(f"\n{len(rows)} residues; {int(candidate.sum())} candidate residues in {len(sites)} sites")
    print("Top candidate sites:")
    for r in site_rows[:8]:
        print(f"  {r['residues']:>10}  {r['n_res']:>2} res  phi={r['mean_phi']:>6} ({r['charge']})  "
              f"rSASA={r['mean_rsasa']}  chains={r['mean_nchains']}")
    print(f"\nHotspot check:")
    for rid in [HOTSPOT] + list(range(*OLD_HOTSPOTS, 5)):
        i = np.where(R["resid"] == rid)[0]
        if len(i): i = i[0]; print(f"  res {rid} ({rows[i]['aa']}): rSASA={R['rsasa'][i]:.2f} "
            f"phi={R['phi'][i]:.2f} concavity={R['concavity'][i]} chains={R['n_chains'][i]} cand={rows[i]['candidate']}")

    # write B-factor PDBs (phi, concavity) on chain A for surface viewing
    phi_by = {rows[i]["resid"]: R["phi"][i] for i in range(len(rows))}
    con_by = {rows[i]["resid"]: R["concavity"][i] for i in range(len(rows))}
    for tag, mp in [("phi", phi_by), ("concavity", con_by)]:
        out = []
        for l in open(PDB):
            if l.startswith("ATOM") and l[21] == PROT:
                rid = int(l[22:26]); v = mp.get(rid, 0.0)
                out.append(f"{l[:60]}{v:6.2f}{l[66:].rstrip()}")
        open(f"{ROOT}/outputs/chainA_{tag}.pdb", "w").write("\n".join(out) + "\n")

    # ---- figure ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
        "font.size":8,"axes.titlesize":9.5,"axes.labelsize":9,"xtick.labelsize":7,"ytick.labelsize":7,
        "legend.fontsize":7,"text.color":"black","axes.labelcolor":"black","axes.titlecolor":"black",
        "xtick.color":"black","ytick.color":"black","axes.edgecolor":"black",
        "pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
    x = R["resid"]
    fig, ax = plt.subplots(5, 1, figsize=(7.0, 8.8), sharex=True)
    fig.suptitle("Target site mapping — VP3 protomer in assembled hexamer", fontsize=11, fontweight="bold", y=0.997)
    def deco(a):
        for st in sites: a.axvspan(st[0], st[-1], color="#f3df9f", alpha=0.55, lw=0)
        a.axvline(HOTSPOT, color="#b3001b", ls="--", lw=1.1)
        a.axvspan(*OLD_HOTSPOTS, color="#cccccc", alpha=0.35, lw=0)
        a.margins(x=0); a.tick_params(width=0.8, length=2.5)
    panels = [("rsasa","Relative SASA","#0b7a6f",(0,1)),
              ("assembly_burial","Assembly burial (ΔrSASA)","#7b5ea7",None),
              ("phi","Electrostatic φ (− acidic / + basic)","#3b6ea5",None),
              ("concavity","Concavity (atom density ≤10 Å)","#c4622d",None),
              ("n_chains","# chains ≤12 Å","#444444",None)]
    for a,(k,lab,c,yl) in zip(ax, panels):
        if k=="phi":
            a.fill_between(x,0,R[k],where=(R[k]>=0),color="#3b6ea5",alpha=.8,lw=0)
            a.fill_between(x,0,R[k],where=(R[k]<0),color="#c0504d",alpha=.8,lw=0); a.axhline(0,color="k",lw=.5)
        else: a.plot(x,R[k],color=c,lw=0.9)
        deco(a); a.set_ylabel(lab)
        if yl: a.set_ylim(*yl)
    ax[-1].set_xlabel("VP3 residue number")
    fig.text(0.5,0.004,"Gold = candidate spots (exposed + groove/seam + charged); grey = old inner hotspots 105–120; "
             "red dashed = Asp111. φ from screened Coulomb; concavity = local atom density.",ha="center",fontsize=7)
    fig.tight_layout(rect=[0,0.012,1,0.985])
    for e in ("pdf","png","svg"): fig.savefig(f"{ROOT}/reports/site_map.{e}", dpi=300)
    print("wrote reports/site_map.{pdf,png,svg}; outputs/site_analysis.csv; outputs/site_candidates.csv; chainA_{phi,concavity}.pdb")

if __name__ == "__main__":
    main()
