#!/usr/bin/env python3
"""Stage 2 (model) — build the C3 ring of binder arms on the full hexamer + clash check.

The hexamer (151lp3t3_hexamer_6chain.pdb) is a trimer of dimers with EXACT C3: the operator that
maps chain A->C also maps F->B (RMSD 0.00); A->E / F->D is its square. The dimers are (A,F),(C,B),(E,D).

For a validated arm docked on the A+F dimer (RFd recenters its own frame, so we first superpose the
design's A+F onto the hexamer's A+F to seat the arm on the hexamer), we apply the C3 operator and its
square to make 3 arms = the ring engaging all 6 patches. Then we clash-check:
  (1) ring self-clash  — do the 3 arms collide (esp. near the 3-fold centre)?
  (2) arm vs NON-cognate hexamer subunits — does an arm bump neighbours it wasn't designed against?
A clean model = few/no clashes in both. Cognate contacts (arm0·A,F etc.) are the designed interface.

Usage: 27_ring_hexamer.py --binder-pdb 08_rosetta/relaxed/design_34_s10.pdb --binder-chain B \
         --hexamer inputs/151lp3t3_hexamer_6chain.pdb --out 09_ring/design_34_s10_ring.pdb
"""
from __future__ import annotations
import argparse
import numpy as np

def load(path):
    """Return dict chain -> list of (resnum, atomname, np.xyz, raw_line) and CA-by-residue arrays."""
    atoms = {}
    for l in open(path):
        if not l.startswith("ATOM"): continue
        ch = l[21]; xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        atoms.setdefault(ch, []).append((int(l[22:26]), l[12:16].strip(), xyz, l))
    return atoms

def ca_array(chain_atoms, common=None):
    d = {rn: xyz for rn, an, xyz, _ in chain_atoms if an == "CA"}
    rns = sorted(d) if common is None else common
    return np.array([d[r] for r in rns]), rns

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Qc - R @ Pc
    rmsd = np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean())
    return R, t, rmsd

def apply(R, t, X): return (R @ X.T).T + t

def clash_count(A, B, thr):
    """count heavy-atom pairs <thr between coord arrays A (M,3) and B (N,3); return count, mindist."""
    cnt = 0; mind = 1e9
    for i in range(0, len(A), 200):
        chunk = A[i:i+200]
        d = np.sqrt(((chunk[:, None, :] - B[None, :, :]) ** 2).sum(2))
        cnt += int((d < thr).sum()); mind = min(mind, float(d.min()))
    return cnt, mind

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binder-pdb", required=True)
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--hexamer", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--clash", type=float, default=2.5)
    a = ap.parse_args()

    hexa = load(a.hexamer)
    design = load(a.binder_pdb)

    # common residues across hexamer chains + design A/F
    hchains = ["A", "B", "C", "D", "E", "F"]
    common = sorted(set.intersection(*[{rn for rn, an, x, l in hexa[c] if an == "CA"} for c in hchains]))

    # C3 operator from the hexamer: A -> C (verified to also send F -> B)
    Aca, _ = ca_array(hexa["A"], common); Cca, _ = ca_array(hexa["C"], common)
    R3, t3, rc3 = kabsch(Aca, Cca)
    Fca_h, _ = ca_array(hexa["F"], common); Bca_h, _ = ca_array(hexa["B"], common)
    fb_rmsd = np.sqrt(((apply(R3, t3, Fca_h) - Bca_h) ** 2).sum(1).mean())
    print(f"C3 operator A->C: rmsd {rc3:.3f}; same op sends F->B rmsd {fb_rmsd:.3f} (expect ~0)")

    # placement: design (A+F) CA -> hexamer (A+F) CA
    dA, rA = ca_array(design["A"]); dF, rF = ca_array(design["F"])
    hA = np.array([x for rn, an, x, l in hexa["A"] if an == "CA" and rn in rA])
    hF = np.array([x for rn, an, x, l in hexa["F"] if an == "CA" and rn in rF])
    Pd = np.vstack([dA, dF]); Qh = np.vstack([hA, hF])
    Rp, tp, rp = kabsch(Pd, Qh)
    print(f"placement (design A+F -> hexamer A+F): rmsd {rp:.3f} (expect ~0, target was rigid)")

    # seat the binder on the hexamer, then C3-replicate
    binder = design[a.binder_chain]
    b0 = np.array([apply(Rp, tp, x) for rn, an, x, l in binder])
    b1 = apply(R3, t3, b0)
    b2 = apply(R3, t3, b1)
    copies = {"P": b0, "Q": b1, "R": b2}
    cognate = {"P": ("A", "F"), "Q": ("C", "B"), "R": ("E", "D")}

    # ---- write hexamer + ring ----
    out = []
    ser = 0
    for c in hchains:
        for rn, an, x, l in hexa[c]:
            ser += 1
            out.append(f"{l[:6]}{ser:5d}{l[11:54]}{l[54:]}".rstrip("\n"))
        out.append("TER")
    binfo = [(rn, l) for rn, an, x, l in binder]
    for cid, coords in copies.items():
        for (rn, an, x, l), nx in zip(binder, coords):
            ser += 1
            newl = f"{l[:6]}{ser:5d}{l[11:21]}{cid}{l[22:30]}{nx[0]:8.3f}{nx[1]:8.3f}{nx[2]:8.3f}{l[54:].rstrip(chr(10))}"
            out.append(newl)
        out.append("TER")
    out.append("END")
    import os
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write("\n".join(out) + "\n")
    print(f"wrote {a.out}  (hexamer A-F + ring arms P,Q,R)")

    # ---- clash analysis ----
    hex_all = {c: np.array([x for rn, an, x, l in hexa[c]]) for c in hchains}
    print(f"\nclash threshold = {a.clash} A (heavy-atom)")
    print("--- ring self-clash (arm vs arm) ---")
    for x, y in [("P", "Q"), ("Q", "R"), ("P", "R")]:
        cnt, mind = clash_count(copies[x], copies[y], a.clash)
        print(f"  arm {x} vs arm {y}: {cnt} clashing atoms, min dist {mind:.2f} A")
    print("--- arm vs NON-cognate hexamer subunits (designed-against are excluded) ---")
    for cid in copies:
        noncog = [c for c in hchains if c not in cognate[cid]]
        tot = 0; worst = 1e9
        for c in noncog:
            cnt, mind = clash_count(copies[cid], hex_all[c], a.clash)
            tot += cnt; worst = min(worst, mind)
        cogcnt = sum(clash_count(copies[cid], hex_all[c], a.clash)[0] for c in cognate[cid])
        print(f"  arm {cid} (cognate {cognate[cid]}): {tot} clashes vs non-cognate {noncog}, "
              f"min {worst:.2f} A  | (cognate-interface clashes: {cogcnt})")

if __name__ == "__main__":
    main()
