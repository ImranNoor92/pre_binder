#!/usr/bin/env python3
"""
Compute the 4 validation metrics from an AlphaFold-multimer output directory.
Runs under .venv-af2 (needs numpy, biopython, freesasa).

Metrics (per the report):
  binder_plddt   mean pLDDT over the binder chain residues          (filter: > 70 on 0-100)
  iptm           interface pTM from the top-ranked model            (filter: > 0.65)
  min_subunit_sasa   smallest per-subunit interface dSASA (A^2)     (filter: all 3 > 200)
  rmsd_to_design CA-RMSD of predicted binder vs trimerized design   (filter: < 3.0 A)

Usage:
  af2_metrics.py --af2-dir DIR --layout layout.json [--design-pdb trimer.pdb] [--out-json m.json]
"""
from __future__ import annotations
import argparse
import json
import pickle
import sys
import tempfile
from pathlib import Path

import numpy as np


def load_top_result(af2_dir: Path):
    ranking = json.load((af2_dir / "ranking_debug.json").open())
    top = ranking["order"][0]
    with (af2_dir / f"result_{top}.pkl").open("rb") as fh:
        result = pickle.load(fh)
    # AF2 writes the top model as ranked_0.pdb
    pdb = af2_dir / "ranked_0.pdb"
    return top, result, pdb


def read_ca(pdb: Path, chain: str):
    """Return {resnum: (x,y,z)} for CA atoms of `chain`."""
    out = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
            out[int(line[22:26])] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out


def detect_binder_chain(pdb: Path, binder_len: int) -> str:
    """AF2 relabels output chains (e.g. inputs A-F,G come out B-G,H). The binder is the chain
    whose CA count matches the fused-trimer length; fall back to the last chain seen."""
    ca, order = {}, []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c = line[21]
            if c not in ca:
                order.append(c)
            ca[c] = ca.get(c, 0) + 1
    for c, n in ca.items():
        if n == binder_len:
            return c
    return order[-1] if order else "A"


def subunit_resranges(layout):
    """1-based residue ranges of the 3 subunits within the binder chain (linkers excluded)."""
    L = layout["binder"]["len"]
    sl = layout["subunit_len"]
    ln = len(layout["linker"])
    starts = [1, sl + ln + 1, 2 * (sl + ln) + 1]
    return [(s, s + sl - 1) for s in starts]  # inclusive


def subunit_dsasa(pdb: Path, binder_chain: str, ranges):
    """Per-subunit interface dSASA = SASA(subunit isolated) - SASA(subunit in complex)."""
    import freesasa
    lines = pdb.read_text().splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]

    def calc(sel_lines):
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=True) as tf:
            tf.write("\n".join(sel_lines) + "\nEND\n")
            tf.flush()
            st = freesasa.Structure(tf.name)
            res = freesasa.calc(st)
            return res.totalArea()

    results = []
    for (lo, hi) in ranges:
        sub_lines = [l for l in atom_lines
                     if l[21] == binder_chain and lo <= int(l[22:26]) <= hi]
        complex_minus = atom_lines  # full complex
        iso = calc(sub_lines)
        # subunit area within complex: total(complex) - total(complex without subunit)
        without = [l for l in atom_lines
                   if not (l[21] == binder_chain and lo <= int(l[22:26]) <= hi)]
        in_complex = calc(complex_minus) - calc(without)
        results.append(round(iso - in_complex, 1))
    return results


def rmsd_to_design(pred_pdb: Path, design_pdb: Path, pred_chain: str, design_chain: str, ranges):
    """CA-RMSD of predicted binder vs trimerized design over subunit residues (Kabsch)."""
    from Bio.SVDSuperimposer import SVDSuperimposer
    pc = read_ca(pred_pdb, pred_chain)
    dc = read_ca(design_pdb, design_chain)
    common = [r for (lo, hi) in ranges for r in range(lo, hi + 1) if r in pc and r in dc]
    if len(common) < 3:
        return None
    P = np.array([pc[r] for r in common])
    D = np.array([dc[r] for r in common])
    sup = SVDSuperimposer()
    sup.set(D, P)
    sup.run()
    return round(float(sup.get_rms()), 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af2-dir", required=True)
    ap.add_argument("--layout", required=True)
    ap.add_argument("--design-pdb", default=None)
    ap.add_argument("--binder-chain", default="G")
    ap.add_argument("--plddt-min", type=float, default=70.0)
    ap.add_argument("--iptm-min", type=float, default=0.65)
    ap.add_argument("--sasa-min", type=float, default=200.0)
    ap.add_argument("--rmsd-max", type=float, default=3.0)
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args()

    af2_dir = Path(args.af2_dir)
    layout = json.load(open(args.layout))
    top, result, pred_pdb = load_top_result(af2_dir)

    b = layout["binder"]
    plddt = np.asarray(result["plddt"])
    binder_plddt = float(plddt[b["start"]: b["start"] + b["len"]].mean())
    iptm = float(result.get("iptm", float("nan")))
    ptm = float(result.get("ptm", float("nan")))

    ranges = subunit_resranges(layout)
    # AF2 relabels chains, so detect the binder chain in the prediction; the design PDB uses --binder-chain.
    pred_binder = detect_binder_chain(pred_pdb, b["len"]) if pred_pdb.exists() else args.binder_chain
    sasa = subunit_dsasa(pred_pdb, pred_binder, ranges) if pred_pdb.exists() else None
    rmsd = (rmsd_to_design(pred_pdb, Path(args.design_pdb), pred_binder, args.binder_chain, ranges)
            if args.design_pdb and pred_pdb.exists() else None)

    checks = {
        "plddt": binder_plddt > args.plddt_min,
        "iptm": iptm > args.iptm_min,
        "sasa": (sasa is not None and min(sasa) > args.sasa_min),
        "rmsd": (rmsd is not None and rmsd < args.rmsd_max),
    }
    m = {
        "af2_dir": str(af2_dir), "top_model": top,
        "binder_plddt": round(binder_plddt, 2), "iptm": round(iptm, 4), "ptm": round(ptm, 4),
        "subunit_sasa": sasa, "min_subunit_sasa": (min(sasa) if sasa else None),
        "rmsd_to_design": rmsd,
        "checks": checks, "pass": all(checks.values()),
    }
    print(json.dumps(m))
    if args.out_json:
        Path(args.out_json).write_text(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()
