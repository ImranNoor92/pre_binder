#!/usr/bin/env python3
"""T=1 vs T=3 interface comparison + combined site-selection map.

Question: which surface spots are at an inter-monomer junction in the T=3 hexamer but NOT in the
natural T=1 capsid (1LP3) -> genuinely assembly(T=3)-unique seams.

T=3 interface : in 151lp3t3 hexamer, chain A residues with heavy-atom contact (<5 A) to chains B-F.
T=1 interface : rebuild the 1LP3 (AAV2) capsid from its 60 BIOMT operators; chain-A residues with
                heavy-atom contact (<5 A) to any neighbouring copy.
Numbering      : 1LP3 (80-598) -> hexamer (1-504) by sequence alignment.
Combined with  : exposure/phi/concavity/n_chains (site_analysis.csv) + real MD RMSF (averaged/residue).
Output: outputs/t1t3_site_table.csv  +  reports/site_selection_combined.{pdf,png,svg}
"""
from __future__ import annotations
import csv, numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner

ROOT = "/data/binder_software/pre-binder"
HEX = f"{ROOT}/inputs/151lp3t3_hexamer_6chain.pdb"
T1 = f"{ROOT}/inputs/1LP3_T1_aav2.pdb"
RMSF = f"{ROOT}/inputs/rmsf_BB_COM_atomcenter.xvg"
CONTACT = 5.0
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def heavy(res):
    return [a.coord for a in res if a.element != "H"]

def chain_atoms_resids(model, cid):
    xyz, rid = [], []
    for res in model[cid]:
        if res.resname not in THREE2ONE: continue
        for a in res:
            if a.element != "H": xyz.append(a.coord); rid.append(res.id[1])
    return np.array(xyz), np.array(rid)

def seq_of(model, cid):
    items = [(r.id[1], THREE2ONE[r.resname]) for r in model[cid] if r.resname in THREE2ONE]
    return "".join(a for _, a in items), [n for n, _ in items]

def parse_biomt(fn):
    ops = {}
    for l in open(fn):
        if l[11:23].strip().startswith("BIOMT"):
            s = l.split()
            row = int(s[2][-1]); op = int(s[3]); vals = list(map(float, s[4:8]))
            ops.setdefault(op, np.zeros((3, 4)))[row - 1] = vals
    return [ops[k] for k in sorted(ops)]

def main():
    p = PDBParser(QUIET=True)
    hexm = p.get_structure("h", HEX)[0]
    t1m = p.get_structure("t", T1)[0]

    # --- T=3 interface contacts per chain-A residue ---
    aA_xyz, aA_rid = chain_atoms_resids(hexm, "A")
    other = np.vstack([chain_atoms_resids(hexm, c)[0] for c in "BCDEF"])
    tree_o = cKDTree(other)
    t3 = {}
    for x, r in zip(aA_xyz, aA_rid):
        if tree_o.query_ball_point(x, CONTACT): t3[r] = t3.get(r, 0) + 1

    # --- T=1: build capsid neighbours from BIOMT, contacts per 1LP3 chain-A residue ---
    ops = parse_biomt(T1)
    t1A_xyz, t1A_rid = chain_atoms_resids(t1m, "A")
    ref_tree = cKDTree(t1A_xyz)
    t1 = {}
    for op in ops:
        R, t = op[:, :3], op[:, 3]
        img = (t1A_xyz @ R.T) + t
        if np.allclose(img, t1A_xyz, atol=1e-3): continue   # identity = self
        # any image atom within CONTACT of the reference monomer?
        nb = cKDTree(img)
        pairs = ref_tree.query_ball_tree(nb, CONTACT)
        for i, hit in enumerate(pairs):
            if hit: t1[t1A_rid[i]] = t1.get(t1A_rid[i], 0) + 1
    print(f"T=3: {len(t3)} interface residues; T=1: {len({k for k,v in t1.items() if v})} interface residues")

    # --- map 1LP3 numbering -> hexamer numbering by sequence alignment ---
    sH, nH = seq_of(hexm, "A"); sT, nT = seq_of(t1m, "A")
    aln = PairwiseAligner(); aln.mode = "global"; aln.open_gap_score = -5; aln.extend_gap_score = -0.5
    a = aln.align(sH, sT)[0]
    map_t1_to_hex = {}
    for (h0, h1), (t0, t1b) in zip(a.aligned[0], a.aligned[1]):
        for k in range(h1 - h0):
            map_t1_to_hex[nT[t0 + k]] = nH[h0 + k]
    t1_in_hex = {}
    for r1, v in t1.items():
        if r1 in map_t1_to_hex: t1_in_hex[map_t1_to_hex[r1]] = t1_in_hex.get(map_t1_to_hex[r1], 0) + v

    # --- real MD RMSF averaged per residue (504 x 180), nm -> Angstrom ---
    rmsf_acc = {}
    for l in open(RMSF):
        l = l.strip()
        if not l or l[0] in "#@": continue
        r, v = l.split()[:2]; r = int(r)
        rmsf_acc.setdefault(r, []).append(float(v) * 10.0)
    rmsf = {r: float(np.mean(v)) for r, v in rmsf_acc.items()}

    # --- merge with site_analysis.csv (rsasa, phi, concavity, n_chains) ---
    site = {int(r["resid"]): r for r in csv.DictReader(open(f"{ROOT}/outputs/site_analysis.csv"))}
    rows = []
    for rid in sorted(site):
        s = site[rid]
        t3c = t3.get(rid, 0); t1c = t1_in_hex.get(rid, 0)
        rows.append({"resid": rid, "aa": s["aa"], "rsasa": float(s["rsasa"]),
            "phi": float(s["phi"]), "concavity": int(s["concavity"]), "n_chains": int(s["n_chains"]),
            "rmsf": round(rmsf.get(rid, np.nan), 2), "t3_iface": t3c, "t1_iface": t1c,
            "t3_unique": bool(t3c > 0 and t1c == 0)})
    # candidate T=3-unique seam spots: exposed + concave + charged + T3-unique
    dens_thr = np.percentile([r["concavity"] for r in rows if r["rsasa"] > 0.25], 60)
    phi_thr = np.percentile([abs(r["phi"]) for r in rows], 60)
    for r in rows:
        r["candidate"] = bool(r["rsasa"] > 0.25 and r["t3_unique"] and r["concavity"] >= dens_thr
                              and abs(r["phi"]) >= phi_thr)
    with open(f"{ROOT}/outputs/t1t3_site_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    cand = [r["resid"] for r in rows if r["candidate"]]
    # cluster candidates into sites
    sites, cur = [], []
    for rid in cand:
        if cur and rid - cur[-1] > 3: sites.append(cur); cur = []
        cur.append(rid)
    if cur: sites.append(cur)
    print(f"\n{int(sum(r['t3_unique'] for r in rows))} T=3-unique interface residues; "
          f"{len(cand)} top candidates in {len(sites)} seams:")
    for st in sorted(sites, key=len, reverse=True)[:8]:
        m = [r for r in rows if r["resid"] in st]
        ph = np.mean([r["phi"] for r in m])
        print(f"  {st[0]}-{st[-1]} ({len(st)} res)  phi={ph:+.2f} ({'acidic' if ph<0 else 'basic'})  "
              f"rSASA={np.mean([r['rsasa'] for r in m]):.2f}")

    # ---- combined figure ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
        "font.size":8,"axes.titlesize":9.5,"axes.labelsize":8.5,"xtick.labelsize":7,"ytick.labelsize":7,
        "legend.fontsize":7,"text.color":"black","axes.labelcolor":"black","axes.titlecolor":"black",
        "xtick.color":"black","ytick.color":"black","axes.edgecolor":"black",
        "pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
    x = np.array([r["resid"] for r in rows])
    g = lambda k: np.array([r[k] for r in rows], float)
    fig, ax = plt.subplots(6, 1, figsize=(7.2, 9.4), sharex=True)
    fig.suptitle("Site selection for a T=3-specific binder (VP3 protomer)", fontsize=11, fontweight="bold", y=0.998)
    def deco(a):
        for st in sites: a.axvspan(st[0], st[-1], color="#fff2a8", alpha=0.9, lw=0)
    ax[0].plot(x, g("rsasa"), color="#0b7a6f", lw=.9); ax[0].set_ylabel("Rel. SASA\n(exposure)"); ax[0].set_ylim(0,1)
    ax[1].fill_between(x,0,g("phi"),where=(g("phi")>=0),color="#3b6ea5",alpha=.85,lw=0)
    ax[1].fill_between(x,0,g("phi"),where=(g("phi")<0),color="#c0504d",alpha=.85,lw=0); ax[1].axhline(0,color="k",lw=.5)
    ax[1].set_ylabel("Electrostatic φ\n(− acid / + base)")
    ax[2].plot(x, g("concavity"), color="#c4622d", lw=.9); ax[2].set_ylabel("Concavity\n(atom density)")
    ax[3].plot(x, g("rmsf"), color="#7b5ea7", lw=.9); ax[3].set_ylabel("MD RMSF (Å)\n(flexibility)")
    ax[4].plot(x, g("t3_iface"), color="#1d3557", lw=.9, label="T=3 hexamer")
    ax[4].plot(x, g("t1_iface"), color="#e07b39", lw=.9, label="T=1 capsid")
    ax[4].set_ylabel("Interface\ncontacts"); ax[4].legend(frameon=False, loc="upper right", ncol=2)
    ax[5].fill_between(x,0,g("t3_unique"),color="#2a7a2a",lw=0,step="mid")
    ax[5].set_ylabel("T=3-UNIQUE\nseam"); ax[5].set_ylim(0,1.2); ax[5].set_yticks([0,1])
    for a in ax: deco(a); a.margins(x=0); a.tick_params(width=0.8,length=2.5)
    ax[-1].set_xlabel("VP3 residue number")
    fig.text(0.5,0.004,"Yellow = recommended spots: exposed + concave + charged + interface present in T=3 hexamer "
             "but ABSENT in T=1 capsid (assembly-unique). Panel 5: T=3 (navy) vs T=1 (orange) interface contacts.",
             ha="center", fontsize=7)
    fig.tight_layout(rect=[0,0.012,1,0.985])
    for e in ("pdf","png","svg"): fig.savefig(f"{ROOT}/reports/site_selection_combined.{e}", dpi=300)
    print("\nwrote reports/site_selection_combined.{pdf,png,svg}; outputs/t1t3_site_table.csv")

if __name__ == "__main__":
    main()
