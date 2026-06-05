#!/usr/bin/env python3
"""
Orchestrate the AlphaFold-multimer validation phases. Runs under .venv-af2.
Shells out to the single-venv helpers so each tool runs in its own environment:
  ProteinMPNN  -> lib/mpnn_subunit.sh   (.venv-rfd-gpu / torch)
  AlphaFold    -> lib/run_af2.sh        (.venv-af2 / jax, trial_1 install)
  FASTA build  -> lib/seqtools.py
  metrics      -> lib/af2_metrics.py

Modes:
  gate   Phase 2b — 1 quick MPNN seq per backbone, AF2, 4 filters -> validated backbone list.
  final  Phase 4  — AF2 each Phase-3 MPNN seq, 4 filters, dimer-only ACID TEST, rank.
"""
from __future__ import annotations
import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

LIB = Path(__file__).resolve().parent
PROJECT = LIB.parent.parent
OUT = PROJECT / "outputs"
TARGET_PDB = PROJECT / "inputs" / "1lp3_hexamer_trimmed_fixed.pdb"
PY = sys.executable  # .venv-af2 python

HEX = "ABCDEF"
DIMER = "AE"          # one dimer pair (A<->E) for the specificity acid test
IPTM_DROP_MIN = 0.15  # required interface-pTM drop hexamer->dimer for hexamer-specificity


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, **kw)


def mpnn_seqs(backbone: Path, outdir: Path, n: int, gpu: int):
    sh([str(LIB / "mpnn_subunit.sh"), str(backbone), str(outdir), str(n), str(gpu)])
    fa = outdir / "seqs" / f"{backbone.stem}.fa"
    out = sh([PY, str(LIB / "seqtools.py"), "mpnn-seqs", str(fa)], capture_output=True)
    return [s for s in out.stdout.splitlines() if s.strip()]


def af2_predict(subunit_seq: str, target_chains: str, tag: str,
                af2_root: Path, design_pdb: Path | None, gpu: int) -> dict:
    """Build FASTA, run AF2 on (target_chains + fused trimer), return metrics dict."""
    af2_root.mkdir(parents=True, exist_ok=True)
    fdir = af2_root / "_fasta"; fdir.mkdir(exist_ok=True)
    fasta = fdir / f"{tag}.fasta"
    layout = fdir / f"{tag}.layout.json"
    sh([PY, str(LIB / "seqtools.py"), "build-fasta",
        "--subunit-seq", subunit_seq, "--target-pdb", str(TARGET_PDB),
        "--target-chains", target_chains, "--out", str(fasta), "--layout", str(layout)])
    sh([str(LIB / "run_af2.sh"), str(fasta), str(af2_root), str(gpu)])
    af2_dir = af2_root / tag
    cmd = [PY, str(LIB / "af2_metrics.py"), "--af2-dir", str(af2_dir), "--layout", str(layout)]
    if design_pdb:
        cmd += ["--design-pdb", str(design_pdb)]
    res = sh(cmd, capture_output=True)
    return json.loads(res.stdout.strip().splitlines()[-1])


def cmd_gate(args):
    designs = sorted((OUT / "01_rfdiffusion_pilot").glob("design_*.pdb"))
    mpnn_dir = OUT / "02b_gate_mpnn"
    af2_root = OUT / "02b_af2"
    valid_dir = OUT / "02b_af2_validated"; valid_dir.mkdir(parents=True, exist_ok=True)
    rows, validated = [], []
    for bb in designs:
        name = bb.stem
        trimer = OUT / "02a_trimerized" / f"{name}_trimer.pdb"
        print(f"[gate] {name}")
        seqs = mpnn_seqs(bb, mpnn_dir, 1, args.gpu)
        if not seqs:
            print(f"  ! no MPNN seq for {name}"); continue
        m = af2_predict(seqs[0], HEX, f"{name}_hex", af2_root, trimer, args.gpu)
        m["design"] = name
        rows.append(m)
        print(f"  plddt={m['binder_plddt']} iptm={m['iptm']} "
              f"min_sasa={m['min_subunit_sasa']} rmsd={m['rmsd_to_design']} pass={m['pass']}")
        if m["pass"]:
            validated.append(name)
            shutil.copy(trimer, valid_dir / f"{name}_trimer.pdb")
    write_csv(OUT / "02b_af2_metrics.csv", rows)
    (valid_dir / "validated.txt").write_text("\n".join(validated) + ("\n" if validated else ""))
    print(f"\n[gate] {len(validated)}/{len(designs)} backbones validated -> {valid_dir}")


def cmd_final(args):
    valid_list = OUT / "02b_af2_validated" / "validated.txt"
    if valid_list.exists() and valid_list.read_text().strip():
        names = [l.strip() for l in valid_list.read_text().splitlines() if l.strip()]
    else:
        names = [p.stem for p in sorted((OUT / "01_rfdiffusion_pilot").glob("design_*.pdb"))]
        print("[final] no gate list; falling back to all designs")
    mpnn_dir = OUT / "03_mpnn_sequences"
    af2_root = OUT / "04_af2"
    final_dir = OUT / "04_final_ranked"; final_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in names:
        bb = OUT / "01_rfdiffusion_pilot" / f"{name}.pdb"
        trimer = OUT / "02a_trimerized" / f"{name}_trimer.pdb"
        fa = mpnn_dir / "seqs" / f"{name}.fa"
        seqs = ([s for s in sh([PY, str(LIB/"seqtools.py"), "mpnn-seqs", str(fa)],
                                capture_output=True).stdout.splitlines() if s.strip()]
                if fa.exists() else mpnn_seqs(bb, mpnn_dir, args.num_seq, args.gpu))
        for i, seq in enumerate(seqs):
            tag = f"{name}_s{i}"
            print(f"[final] {tag}")
            m = af2_predict(seq, HEX, f"{tag}_hex", af2_root, trimer, args.gpu)
            row = {"design": name, "seq_idx": i, **{f"hex_{k}": v for k, v in m.items()
                   if k in ("binder_plddt", "iptm", "ptm", "min_subunit_sasa", "rmsd_to_design", "pass")}}
            row["hex_sum_sasa"] = sum(m["subunit_sasa"]) if m["subunit_sasa"] else 0
            if m["pass"]:
                d = af2_predict(seq, DIMER, f"{tag}_dimer", af2_root, None, args.gpu)
                drop = round(m["iptm"] - d["iptm"], 4)
                row.update(dimer_iptm=d["iptm"], iptm_drop=drop,
                           specific=bool(drop >= IPTM_DROP_MIN))
                rmsd = m["rmsd_to_design"] or 99.0
                row["score"] = round(m["iptm"] * (m["binder_plddt"] / 100.0)
                                     * (1.0 / max(rmsd, 0.5)) * (row["hex_sum_sasa"] / 1000.0), 4)
            else:
                row.update(dimer_iptm=None, iptm_drop=None, specific=False, score=0.0)
            rows.append(row)
            print(f"  hex_pass={m['pass']} iptm_drop={row.get('iptm_drop')} "
                  f"specific={row['specific']} score={row['score']}")
    # rank the hexamer-specific passers
    keep = sorted([r for r in rows if r["hex_pass"] and r["specific"]],
                  key=lambda r: r["score"], reverse=True)
    for rank, r in enumerate(keep, 1):
        src = af2_root / f"{r['design']}_s{r['seq_idx']}_hex" / "ranked_0.pdb"
        if src.exists():
            shutil.copy(src, final_dir / f"rank{rank:02d}_{r['design']}_s{r['seq_idx']}.pdb")
    write_csv(OUT / "04_final_metrics.csv", rows)
    print(f"\n[final] {len(keep)} hexamer-specific designs ranked -> {final_dir}")


def write_csv(path: Path, rows: list):
    if not rows:
        path.write_text(""); return
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "design", k))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (json.dumps(v) if isinstance(v, (list, dict)) else v)
                        for k, v in r.items()})
    print(f"  wrote {path} ({len(rows)} rows)")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    g = sub.add_parser("gate"); g.add_argument("--gpu", type=int, default=0); g.set_defaults(func=cmd_gate)
    f = sub.add_parser("final"); f.add_argument("--gpu", type=int, default=0)
    f.add_argument("--num-seq", type=int, default=8); f.set_defaults(func=cmd_final)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
