#!/usr/bin/env python3
"""Stage 1 Step 4 — Rosetta InterfaceAnalyzer on the dimer-arm designs (the binding test).

For each foldable candidate (design_N_sK) from the AF2 fold filter: thread the MPNN sequence onto
the binder backbone (chain B) IN THE COMPLEX, keeping the A+F target dimer verbatim, then score the
binder-vs-dimer interface (AF_B) with InterfaceAnalyzerMover (pack_separated repack). This is the
as-designed docking pose from RFdiffusion, not an AF2-multimer re-prediction (per instruction).

Reports the decision metrics from PLAN.md:
  - dG_separated       (want NEGATIVE = favorable binding)
  - sc_value           (shape complementarity)
  - interface_hbonds, packstat, dSASA
  - delta_unsat_hbonds (buried unsatisfied polars — must be LOW)
  - salt_bridge_K42    (binder Asp/Glu carboxylate within 4.0 A of Lys42 NZ on A or F)
  - anchor_contacts    (binder heavy atoms within 4.5 A of Leu40 / Met155 sidechains = grip)

Run with the af2ig python (has pyrosetta):
  /home/a-mxn833/mambaforge/envs/af2ig/bin/python scripts/23_rosetta_dimer_arm.py
Env: MAXN=<n> caps the number scored (0 = all foldable); TIER=strong restricts to plddt>90 & rmsd<1.5.
Output: outputs/C3_symmetric_Binder_2026_07_02/04_rosetta/interface_metrics.csv
"""
from __future__ import annotations
import csv, glob, json, os, sys

ROOT = "/data/binder_software/pre-binder"
CAMP = f"{ROOT}/outputs/C3_symmetric_Binder_2026_07_02"
RANKED = f"{CAMP}/03_af2/ranked_binderonly.csv"
SEQDIR = f"{CAMP}/02_mpnn/seqs"
RFDDIR = f"{CAMP}/01_rfd"
OUTDIR = f"{CAMP}/04_rosetta"
INDIR  = f"{OUTDIR}/in"
OUT    = f"{OUTDIR}/interface_metrics.csv"
os.makedirs(INDIR, exist_ok=True)

MAXN   = int(os.environ.get("MAXN", "0"))         # 0 = all foldable
TIER   = os.environ.get("TIER", "")               # "strong" = plddt>90 & rmsd<1.5
SHARD  = int(os.environ.get("SHARD", "0"))        # this worker's index
NSHARD = int(os.environ.get("NSHARD", "1"))       # total workers (candidates split by index)
JSONL  = f"{OUTDIR}/metrics_s{SHARD}.jsonl"

MERGE = len(sys.argv) > 1 and sys.argv[1] == "merge"
if not MERGE:
    import pyrosetta as pr
    pr.init("-ignore_unrecognized_res -load_PDB_components false -mute all -holes:dalphaball /dev/null",
            silent=True)
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

ONE_TO_THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
"G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
"S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}


def seq_for(tag):
    """tag = design_N_sK -> the K-th designed MPNN sequence (record 0 is the native poly-G)."""
    bb, s = tag.rsplit("_s", 1)
    recs, cur = [], None
    for line in open(f"{SEQDIR}/{bb}.fa"):
        if line.startswith(">"): cur = []; recs.append(cur)
        elif cur is not None: cur.append(line.strip())
    return ["".join(r) for r in recs][1:][int(s)]


def build_complex(design_pdb, seq, out):
    """Thread MPNN seq onto chain B; keep target chains A + F verbatim. Emit A, F, B."""
    lines = open(design_pdb).read().splitlines()
    binder = [l for l in lines if l.startswith("ATOM") and l[21] == "B"]
    order, seen = [], set()
    for l in binder:
        if l[22:26] not in seen:
            seen.add(l[22:26]); order.append(l[22:26])
    if len(seq) != len(order):
        raise ValueError(f"seq {len(seq)} != chain B {len(order)}")
    ridx = {r: i for i, r in enumerate(order)}
    out_lines, ser = [], 0
    # target A + F unchanged
    for c in ("A", "F"):
        for l in lines:
            if l.startswith("ATOM") and l[21] == c:
                ser += 1
                out_lines.append(f"{l[:6]}{ser:5d}{l[11:]}")
        out_lines.append("TER")
    # binder chain B, threaded
    for l in binder:
        i = ridx[l[22:26]]; ser += 1
        out_lines.append(f"{l[:6]}{ser:5d}{l[11:17]}{ONE_TO_THREE[seq[i]]:>3s} B{i+1:4d}{l[26:]}")
    out_lines += ["TER", "END"]
    open(out, "w").write("\n".join(out_lines) + "\n")


def geom_checks(pose):
    """Salt bridge to Lys42 (either dimer chain) and hydrophobic grip on Leu40/Met155."""
    info = pose.pdb_info()
    # map: collect binder (chain B) carboxyl O atoms, and target anchor sidechain atoms
    binder_carboxyl, k42_nz, anchor_atoms = [], [], []
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i); ch = info.chain(i); rnum = info.number(i); rn = res.name3()
        if ch == "B" and rn in ("ASP", "GLU"):
            for an in ("OD1", "OD2", "OE1", "OE2"):
                if res.has(an): binder_carboxyl.append(res.xyz(an))
        if ch in ("A", "F") and rnum == 42 and rn == "LYS" and res.has("NZ"):
            k42_nz.append(res.xyz("NZ"))
        if ch in ("A", "F") and rnum in (40, 155):   # Leu40 / Met155 sidechain heavy atoms
            for a in range(5, res.nheavyatoms() + 1):  # skip backbone N,CA,C,O,CB start ~5
                anchor_atoms.append(res.xyz(a))
    def mind(pts_a, pts_b):
        m = 99.0
        for p in pts_a:
            for q in pts_b:
                d = (p - q).norm()
                if d < m: m = d
        return m
    sb = mind(binder_carboxyl, k42_nz) if (binder_carboxyl and k42_nz) else 99.0
    # count binder heavy atoms within 4.5 A of an anchor sidechain atom
    n_anchor = 0
    binder_heavy = []
    for i in range(1, pose.total_residue() + 1):
        if info.chain(i) != "B": continue
        res = pose.residue(i)
        for a in range(1, res.nheavyatoms() + 1):
            binder_heavy.append(res.xyz(a))
    for p in binder_heavy:
        for q in anchor_atoms:
            if (p - q).norm() < 4.5:
                n_anchor += 1; break
    return round(sb, 2), n_anchor


def score(pdb):
    pose = pr.pose_from_pdb(pdb)
    iam = InterfaceAnalyzerMover()
    iam.set_interface("AF_B")                       # target dimer (A+F) vs binder (B)
    iam.set_scorefunction(pr.get_fa_scorefxn())
    iam.set_compute_packstat(True); iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True); iam.set_compute_interface_sc(True); iam.set_pack_separated(True)
    iam.apply(pose)
    d = iam.get_all_data()
    sb, n_anchor = geom_checks(pose)
    return {
        "dG_separated": round(iam.get_interface_dG(), 1),
        "dSASA": round(iam.get_interface_delta_sasa(), 0),
        "sc_value": round(d.sc_value, 3),
        "interface_hbonds": int(d.interface_hbonds),
        "delta_unsat_hbonds": int(d.delta_unsat_hbonds),
        "packstat": round(iam.get_interface_packstat(), 3),
        "interface_nres": int(sum(d.interface_nres)),
        "salt_bridge_K42_dist": sb,
        "salt_bridge_K42": bool(sb < 4.0),
        "anchor_contacts": n_anchor,
    }


def select_candidates():
    rows = list(csv.DictReader(open(RANKED)))
    cands = [r for r in rows if r["foldable"] == "True"]
    if TIER == "strong":
        cands = [r for r in cands
                 if float(r["plddt_binder"]) > 90 and float(r["binder_rmsd"]) < 1.5]
    if MAXN: cands = cands[:MAXN]
    return cands


def merge():
    """Combine all shard jsonls -> interface_metrics.csv, ranked by dG (most negative first)."""
    out, seen = [], set()
    for jf in sorted(glob.glob(f"{OUTDIR}/metrics_s*.jsonl")):
        for l in open(jf):
            r = json.loads(l)
            if "dG_separated" in r and r["tag"] not in seen:
                seen.add(r["tag"]); out.append(r)
    out.sort(key=lambda x: x["dG_separated"])
    keys = ["tag","backbone","plddt","dG_separated","dSASA","sc_value","interface_hbonds",
            "delta_unsat_hbonds","packstat","interface_nres","salt_bridge_K42",
            "salt_bridge_K42_dist","anchor_contacts"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    npos = sum(r["dG_separated"] < 0 for r in out)
    nsb = sum(r["dG_separated"] < 0 and r.get("salt_bridge_K42") for r in out)
    print(f"wrote {OUT} ({len(out)} scored)")
    print(f"  {npos} with NEGATIVE dG_separated; {nsb} of those also satisfy the Lys42 salt bridge")
    print("  best 10 by dG_separated:")
    for r in out[:10]:
        print(f"  {r['tag']}: dG={r['dG_separated']} sc={r['sc_value']} "
              f"buns={r['delta_unsat_hbonds']} K42sb={r['salt_bridge_K42']} anchors={r['anchor_contacts']}")


def main():
    cands = select_candidates()
    mine = [c for i, c in enumerate(cands) if i % NSHARD == SHARD]
    print(f"[shard {SHARD}/{NSHARD}] scoring {len(mine)} of {len(cands)} foldable candidates")
    done = {json.loads(l)["tag"] for l in open(JSONL)} if os.path.exists(JSONL) else set()
    inc = open(JSONL, "a", buffering=1)
    for n, c in enumerate(mine, 1):
        tag = c["tag"]
        if tag in done: continue
        bb = tag.rsplit("_s", 1)[0]
        pdb = f"{INDIR}/{tag}.pdb"
        try:
            if not os.path.exists(pdb):
                build_complex(f"{RFDDIR}/{bb}.pdb", seq_for(tag), pdb)
            r = {"tag": tag, "backbone": bb, "plddt": c["plddt_binder"], **score(pdb)}
        except Exception as e:
            r = {"tag": tag, "backbone": bb, "error": str(e)[:140]}
        print(f"[shard {SHARD} {n}/{len(mine)}] " + json.dumps(r), flush=True)
        inc.write(json.dumps(r) + "\n")
    inc.close()


if __name__ == "__main__":
    if MERGE:
        merge()
    else:
        main()
