#!/usr/bin/env python3
"""
Build a 2-chain PDB for af2_initial_guess: binder subunit (threaded with a designed sequence)
as chain A FIRST, then one target chain as chain B. af2_initial_guess requires exactly
binder+target, binder first, and reads the sequence from residue names.

Threading = rewrite the binder backbone's residue names to the ProteinMPNN sequence; pyrosetta
rebuilds sidechains from those names. Only backbone atoms are needed for the initial guess.

Usage:
  thread_and_pair.py --backbone design_0.pdb --binder-chain G --seq AVAAL... \
                     --target-pdb hexamer.pdb --target-chain A --out design_0_s1.pdb
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


def chain_atoms(pdb: Path, chain: str):
    return [l for l in pdb.read_text().splitlines()
            if l.startswith("ATOM") and l[21] == chain]


def residue_order(atoms):
    seen, order = set(), []
    for l in atoms:
        r = l[22:26]
        if r not in seen:
            seen.add(r); order.append(r)
    return order


def emit(atoms, new_chain, start_serial, seq=None):
    """Renumber residues 1..N on new_chain; if seq given, rewrite residue names to thread it."""
    order = residue_order(atoms)
    res_idx = {r: i for i, r in enumerate(order)}
    if seq is not None and len(seq) != len(order):
        raise SystemExit(f"seq length {len(seq)} != binder residues {len(order)}")
    out, ser = [], start_serial
    for l in atoms:
        ri = res_idx[l[22:26]]
        resname = ONE_TO_THREE[seq[ri]] if seq is not None else l[17:20].strip()
        ser += 1
        # cols: 1-6 rec | 7-11 serial | 12-17 (space+atom+altloc) | 18-20 resName | 21 space | 22 chain | 23-26 resSeq | 27+ rest
        out.append(f"{l[:6]}{ser:5d}{l[11:17]}{resname:>3s} {new_chain}{ri+1:4d}{l[26:]}")
    return out, ser


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", required=True)
    ap.add_argument("--binder-chain", default="G")
    ap.add_argument("--seq", required=True, help="designed binder subunit sequence (1-letter)")
    ap.add_argument("--target-pdb", required=True)
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    binder = chain_atoms(Path(a.backbone), a.binder_chain)
    target = chain_atoms(Path(a.target_pdb), a.target_chain)
    if not binder:
        raise SystemExit(f"no binder atoms (chain {a.binder_chain}) in {a.backbone}")
    if not target:
        raise SystemExit(f"no target atoms (chain {a.target_chain}) in {a.target_pdb}")

    lines, ser = emit(binder, "A", 0, seq=a.seq.strip())     # binder FIRST, threaded
    lines.append("TER")
    tlines, _ = emit(target, "B", ser)                        # target second, native
    lines += tlines + ["TER", "END"]
    Path(a.out).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
