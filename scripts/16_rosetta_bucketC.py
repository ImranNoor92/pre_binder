#!/usr/bin/env python3
"""Bucket C + binding energy — Rosetta InterfaceAnalyzer on the designed complexes (binder threaded
onto the RFdiffusion docking pose against the target). One tool gives both interface geometry and
the interface binding energy (dG_separated). Scored with pack_separated (interface repack); these
are the as-designed interfaces, not the failed AF2 poses.

Outputs outputs/bucketC_metrics.csv (incremental). Run with af2ig python (pyrosetta).
"""
from __future__ import annotations
import csv, json, os, subprocess, sys
import numpy as np
import pyrosetta as pr

ROOT = "/data/binder_software/pre-binder"
REPORT = f"{ROOT}/outputs/top_candidates_report.csv"
SEQDIR = f"{ROOT}/outputs/13_mpnn_central/seqs"
INDIR = f"{ROOT}/outputs/16_rosetta/in"
OUT = f"{ROOT}/outputs/bucketC_metrics.csv"
JSONL = f"{ROOT}/outputs/bucketC.jsonl"
os.makedirs(INDIR, exist_ok=True)

pr.init("-ignore_unrecognized_res -load_PDB_components false -mute all -holes:dalphaball /dev/null", silent=True)
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
HYDRO = set("ACFILMPVWY")
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def seq_for(tag):
    bb, s = tag.rsplit("_s", 1)
    recs, cur = [], None
    for line in open(f"{SEQDIR}/{bb}.fa"):
        if line.startswith(">"): cur = []; recs.append(cur)
        elif cur is not None: cur.append(line.strip())
    return ["".join(r) for r in recs][1:][int(s)]


def score(pdb):
    pose = pr.pose_from_pdb(pdb)
    iam = InterfaceAnalyzerMover(); iam.set_interface("A_B")     # binder A vs target B
    iam.set_scorefunction(pr.get_fa_scorefxn())
    iam.set_compute_packstat(True); iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True); iam.set_compute_interface_sc(True); iam.set_pack_separated(True)
    iam.apply(pose)
    d = iam.get_all_data()
    # interface binder residues (chain A within 5A CB of chain B) for hydrophobicity
    return {
        "dG_separated": round(iam.get_interface_dG(), 1),
        "dSASA": round(iam.get_interface_delta_sasa(), 0),
        "dG_dSASA_ratio": round(d.dG_dSASA_ratio * 100, 3),
        "sc_value": round(d.sc_value, 3),
        "interface_hbonds": int(d.interface_hbonds),
        "interface_nres": int(sum(d.interface_nres)),
        "packstat": round(iam.get_interface_packstat(), 3),
    }


def main():
    cands = list(csv.DictReader(open(REPORT)))
    done = {json.loads(l)["tag"] for l in open(JSONL)} if os.path.exists(JSONL) else set()
    inc = open(JSONL, "a", buffering=1)
    for c in cands:
        tag = c["design"]
        if tag in done: continue
        bb = tag.rsplit("_s", 1)[0]
        seq = seq_for(tag)
        pdb = f"{INDIR}/{tag}.pdb"
        if not os.path.exists(pdb):
            subprocess.run([sys.executable, f"{ROOT}/scripts/lib/build_binder_complex.py",
                            "--design", f"{ROOT}/outputs/11_rfd_central/{bb}.pdb",
                            "--binder-chain", "G", "--seq", seq, "--out", pdb], check=True,
                           stdout=subprocess.DEVNULL)
        try:
            r = {"tag": tag, "group": c["monomers"], **score(pdb)}
        except Exception as e:
            r = {"tag": tag, "group": c["monomers"], "error": str(e)[:120]}
        print(json.dumps(r), flush=True); inc.write(json.dumps(r) + "\n")
    inc.close()
    # write CSV ranked by dG_separated (most favorable first)
    rows = [json.loads(l) for l in open(JSONL)]
    rows = [r for r in rows if "dG_separated" in r]
    rows.sort(key=lambda x: x["dG_separated"])
    keys = ["tag","group","dG_separated","dSASA","dG_dSASA_ratio","sc_value","interface_hbonds","interface_nres","packstat"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {OUT} ({len(rows)} designs)")


if __name__ == "__main__":
    main()
