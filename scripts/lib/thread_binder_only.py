#!/usr/bin/env python3
"""Thread a ProteinMPNN sequence onto the binder backbone and emit it ALONE (chain A, 1 chain).

For the binder-only AF2 foldability / self-consistency filter: af2_initial_guess on a single
chain runs in monomer mode and reports plddt + RMSD of the prediction vs this threaded backbone
(= the RFdiffusion design). No target included.

Usage: thread_binder_only.py --backbone design_0.pdb --binder-chain G --seq AVAAL... --out design_0_s1.pdb
"""
from __future__ import annotations
import argparse
from pathlib import Path

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", required=True)
    ap.add_argument("--binder-chain", default="G")
    ap.add_argument("--seq", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    atoms = [l for l in Path(a.backbone).read_text().splitlines()
             if l.startswith("ATOM") and l[21] == a.binder_chain]
    if not atoms:
        raise SystemExit(f"no binder atoms (chain {a.binder_chain}) in {a.backbone}")
    order, seen = [], set()
    for l in atoms:
        if l[22:26] not in seen:
            seen.add(l[22:26]); order.append(l[22:26])
    res_idx = {r: i for i, r in enumerate(order)}
    seq = a.seq.strip()
    if len(seq) != len(order):
        raise SystemExit(f"seq length {len(seq)} != binder residues {len(order)}")
    out, ser = [], 0
    for l in atoms:
        ri = res_idx[l[22:26]]; ser += 1
        out.append(f"{l[:6]}{ser:5d}{l[11:17]}{ONE_TO_THREE[seq[ri]]:>3s} A{ri+1:4d}{l[26:]}")
    out += ["TER", "END"]
    Path(a.out).write_text("\n".join(out) + "\n")
    print(f"wrote {a.out} ({len(order)} res, binder-only)")


if __name__ == "__main__":
    main()
