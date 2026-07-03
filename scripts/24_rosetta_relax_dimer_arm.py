#!/usr/bin/env python3
"""Stage 1 Step 4b — FastRelax + re-score the top dimer-arm hits.

Pack-only InterfaceAnalyzer (Step 4) gave all-positive dG_separated because unrelaxed RFd/MPNN
poses are clash-dominated. Here we take the best N designs, apply FastRelax with coordinate
constraints to the STARTING coords (relieves clashes locally without letting the fold drift apart),
then re-run InterfaceAnalyzer (AF_B) to get a fair dG_separated. This is the standard binder-eval
protocol. Reuses the complexes already built in 04_rosetta/in/.

Run with af2ig python (pyrosetta):
  SHARD=0 NSHARD=8 /home/.../af2ig/bin/python scripts/24_rosetta_relax_dimer_arm.py     # a shard
  /home/.../af2ig/bin/python scripts/24_rosetta_relax_dimer_arm.py merge                # combine
Env: RELAX_TOPN (default 20) = how many top-dG designs to relax; REPEATS (default 1) FastRelax rounds.
Output: outputs/C3_symmetric_Binder_2026_07_02/04_rosetta/relaxed_metrics.csv
"""
from __future__ import annotations
import csv, glob, json, os, sys

ROOT = "/data/binder_software/pre-binder"
CAMP = f"{ROOT}/outputs/C3_symmetric_Binder_2026_07_02"
PACK_CSV = f"{CAMP}/04_rosetta/interface_metrics.csv"
INDIR = f"{CAMP}/04_rosetta/in"
RELAXDIR = f"{CAMP}/04_rosetta/relaxed"
OUT = f"{CAMP}/04_rosetta/relaxed_metrics.csv"
os.makedirs(RELAXDIR, exist_ok=True)

TOPN    = int(os.environ.get("RELAX_TOPN", "20"))
REPEATS = int(os.environ.get("REPEATS", "1"))
SHARD   = int(os.environ.get("SHARD", "0"))
NSHARD  = int(os.environ.get("NSHARD", "1"))
JSONL   = f"{CAMP}/04_rosetta/relaxed_s{SHARD}.jsonl"

MERGE = len(sys.argv) > 1 and sys.argv[1] == "merge"
if not MERGE:
    import pyrosetta as pr
    pr.init("-ignore_unrecognized_res -load_PDB_components false -mute all -holes:dalphaball /dev/null",
            silent=True)
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, InitializeFromCommandline)
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover


def geom_checks(pose):
    info = pose.pdb_info()
    binder_carboxyl, k42_nz, anchor_atoms, binder_heavy = [], [], [], []
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i); ch = info.chain(i); rnum = info.number(i); rn = res.name3()
        if ch == "B":
            for a in range(1, res.nheavyatoms() + 1):
                binder_heavy.append(res.xyz(a))
            if rn in ("ASP", "GLU"):
                for an in ("OD1", "OD2", "OE1", "OE2"):
                    if res.has(an): binder_carboxyl.append(res.xyz(an))
        if ch in ("A", "F") and rnum == 42 and rn == "LYS" and res.has("NZ"):
            k42_nz.append(res.xyz("NZ"))
        if ch in ("A", "F") and rnum in (40, 155):
            for a in range(5, res.nheavyatoms() + 1):
                anchor_atoms.append(res.xyz(a))
    def mind(pa, pb):
        m = 99.0
        for p in pa:
            for q in pb:
                d = (p - q).norm()
                if d < m: m = d
        return m
    sb = mind(binder_carboxyl, k42_nz) if (binder_carboxyl and k42_nz) else 99.0
    n_anchor = 0
    for p in binder_heavy:
        if any((p - q).norm() < 4.5 for q in anchor_atoms):
            n_anchor += 1
    return round(sb, 2), n_anchor


def interface_movemap(pose, shell=8.0):
    """Free the binder (chain B, bb+chi) and the target sidechains within `shell` A of it.
    Keep the bulk of the 504-res A/F target rigid — it's the fixed capsid epitope, and relaxing
    all ~1000 target residues is both physically wrong and ~10x slower. Returns a MoveMap."""
    info = pose.pdb_info()
    binder_xyz = []
    for i in range(1, pose.total_residue() + 1):
        if info.chain(i) == "B":
            res = pose.residue(i)
            for a in range(1, res.nheavyatoms() + 1):
                binder_xyz.append(res.xyz(a))
    from pyrosetta.rosetta.core.kinematics import MoveMap
    mm = MoveMap()
    for i in range(1, pose.total_residue() + 1):
        ch = info.chain(i)
        mm.set_bb(i, False)                     # ALL backbones rigid: binder = as-designed (AF2-
                                               # validated) pose, target = fixed capsid epitope.
                                               # Free-bb relax from a clashing start explodes.
        if ch == "B":
            mm.set_chi(i, True)                # binder sidechains: repack + minimize
            continue
        res = pose.residue(i)
        near = any(any((res.xyz(a) - q).norm() < shell for q in binder_xyz)
                   for a in range(1, res.nheavyatoms() + 1))
        mm.set_chi(i, near)                     # only interface-shell target sidechains move
    return mm


def relax_and_score(pdb, relaxed_out):
    pose = pr.pose_from_pdb(pdb)
    sfxn = pr.get_fa_scorefxn()
    # Sidechain-only optimization: repack all rotamers (no design), then minimize chi with the
    # backbone rigid. Relieves the rotamer clashes that pack-only (discrete) can't, without any
    # backbone motion (so dSASA/fold stay as-designed). REPEATS rounds of pack+min.
    tf = TaskFactory(); tf.push_back(InitializeFromCommandline()); tf.push_back(RestrictToRepacking())
    prm = PackRotamersMover(sfxn); prm.task_factory(tf)
    mm = interface_movemap(pose)                 # bb all False; chi = binder + interface-shell target
    minm = MinMover(mm, sfxn, "lbfgs_armijo_nonmonotone", 1e-4, True)
    for _ in range(max(1, REPEATS)):
        prm.apply(pose); minm.apply(pose)
    pose.dump_pdb(relaxed_out)
    iam = InterfaceAnalyzerMover(); iam.set_interface("AF_B")
    iam.set_scorefunction(pr.get_fa_scorefxn())
    iam.set_compute_packstat(True); iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True); iam.set_compute_interface_sc(True); iam.set_pack_separated(True)
    iam.apply(pose)
    d = iam.get_all_data()
    sb, n_anchor = geom_checks(pose)
    return {
        "dG_relaxed": round(iam.get_interface_dG(), 1),
        "dSASA": round(iam.get_interface_delta_sasa(), 0),
        "sc_value": round(d.sc_value, 3),
        "interface_hbonds": int(d.interface_hbonds),
        "delta_unsat_hbonds": int(d.delta_unsat_hbonds),
        "packstat": round(iam.get_interface_packstat(), 3),
        "salt_bridge_K42_dist": sb,
        "salt_bridge_K42": bool(sb < 4.0),
        "anchor_contacts": n_anchor,
    }


def merge():
    pack = {r["tag"]: r for r in csv.DictReader(open(PACK_CSV))}
    out, seen = [], set()
    for jf in sorted(glob.glob(f"{CAMP}/04_rosetta/relaxed_s*.jsonl")):
        for l in open(jf):
            r = json.loads(l)
            if "dG_relaxed" in r and r["tag"] not in seen:
                seen.add(r["tag"])
                r["dG_packonly"] = pack.get(r["tag"], {}).get("dG_separated", "")
                out.append(r)
    out.sort(key=lambda x: x["dG_relaxed"])
    keys = ["tag","backbone","dG_relaxed","dG_packonly","dSASA","sc_value","interface_hbonds",
            "delta_unsat_hbonds","packstat","salt_bridge_K42","salt_bridge_K42_dist","anchor_contacts"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    neg = [r for r in out if r["dG_relaxed"] < 0]
    negsb = [r for r in neg if r["salt_bridge_K42"] and r["delta_unsat_hbonds"] <= 5]
    print(f"wrote {OUT} ({len(out)} relaxed)")
    print(f"  {len(neg)} now have NEGATIVE dG_relaxed; "
          f"{len(negsb)} of those also have K42 salt bridge AND buns<=5")
    for r in out[:12]:
        print(f"  {r['tag']}: dG {r['dG_packonly']} -> {r['dG_relaxed']}  sc={r['sc_value']} "
              f"buns={r['delta_unsat_hbonds']} K42sb={r['salt_bridge_K42']} anchors={r['anchor_contacts']}")


def main():
    rows = list(csv.DictReader(open(PACK_CSV)))[:TOPN]     # already sorted by dG ascending
    mine = [r for i, r in enumerate(rows) if i % NSHARD == SHARD]
    print(f"[shard {SHARD}/{NSHARD}] relaxing {len(mine)} of top {TOPN} (REPEATS={REPEATS})")
    done = {json.loads(l)["tag"] for l in open(JSONL)} if os.path.exists(JSONL) else set()
    inc = open(JSONL, "a", buffering=1)
    for n, c in enumerate(mine, 1):
        tag = c["tag"]
        if tag in done: continue
        pdb = f"{INDIR}/{tag}.pdb"
        try:
            r = {"tag": tag, "backbone": c["backbone"],
                 **relax_and_score(pdb, f"{RELAXDIR}/{tag}.pdb")}
        except Exception as e:
            r = {"tag": tag, "backbone": c["backbone"], "error": str(e)[:140]}
        print(f"[shard {SHARD} {n}/{len(mine)}] " + json.dumps(r), flush=True)
        inc.write(json.dumps(r) + "\n")
    inc.close()


if __name__ == "__main__":
    merge() if MERGE else main()
