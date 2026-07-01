#!/usr/bin/env python3
"""Build a 2-chain complex for af2_initial_guess from a central-binder RFd design.

af2_initial_guess takes exactly 2 chains. The RFd output (design_N.pdb) holds the binder
(chain G) + the 6 target chains (A-F), all in one consistent frame. We emit:
  chain A = binder, threaded with the MPNN sequence (FIRST chain -> binderlen)
  chain B = all 6 target chains merged, unique continuous residue numbers (spatial gaps
            become AF2 chain-breaks). pae_interaction = binder vs whole hexamer surface.

NOTE: like any af2_initial_guess multi-protomer case, check target_aligned_rmsd afterwards —
if the (trimmed) hexamer fragments drift, the interface score is not meaningful.

Usage: build_binder_complex.py --design design_148.pdb --binder-chain G --seq AVAAL... --out design_148_s4_cplx.pdb
"""
from __future__ import annotations
import argparse
from pathlib import Path

ONE_TO_THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
"G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
"S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
TARGET_CHAINS = "ABCDEF"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--binder-chain", default="G")
    ap.add_argument("--seq", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    lines = Path(a.design).read_text().splitlines()
    binder = [l for l in lines if l.startswith("ATOM") and l[21] == a.binder_chain]
    if not binder:
        raise SystemExit(f"no binder chain {a.binder_chain} in {a.design}")

    out = []; ser = 0
    # binder -> chain A, threaded
    order, seen = [], set()
    for l in binder:
        if l[22:26] not in seen:
            seen.add(l[22:26]); order.append(l[22:26])
    ridx = {r: i for i, r in enumerate(order)}
    seq = a.seq.strip()
    if len(seq) != len(order):
        raise SystemExit(f"seq length {len(seq)} != binder residues {len(order)}")
    for l in binder:
        i = ridx[l[22:26]]; ser += 1
        out.append(f"{l[:6]}{ser:5d}{l[11:17]}{ONE_TO_THREE[seq[i]]:>3s} A{i+1:4d}{l[26:]}")
    out.append("TER")
    # 6 target chains -> merged chain B, unique continuous residue numbers.
    # Drop contiguous crops < MIN_SEG residues: tiny 1-2-res fragments from the R50 trim
    # become isolated dipeptides that break PyRosetta (RamaPrePro / cyclic dipeptide). They
    # are peripheral bits far from the central site and contribute nothing to the interface.
    MIN_SEG = 4
    rnum = len(order)
    for c in TARGET_CHAINS:
        tch = [l for l in lines if l.startswith("ATOM") and l[21] == c]
        # group lines by contiguous residue-number runs
        runs, cur, prevr = [], [], None
        for l in tch:
            r = int(l[22:26])
            if prevr is None or r == prevr:
                cur.append(l)
            elif r == prevr + 1:
                cur.append(l)
            else:
                runs.append(cur); cur = [l]
            prevr = r
        if cur: runs.append(cur)
        for run in runs:
            nres = len({l[22:26] for l in run})
            if nres < MIN_SEG:
                continue
            prev = None
            for l in run:
                if l[22:26] != prev:
                    rnum += 1; prev = l[22:26]
                ser += 1
                out.append(f"{l[:6]}{ser:5d}{l[11:21]}B{rnum:4d}{l[26:]}")
    out += ["TER", "END"]
    Path(a.out).write_text("\n".join(out) + "\n")
    print(f"wrote {a.out} (binder {len(order)} aa chain A + merged target chain B)")


if __name__ == "__main__":
    main()
