#!/usr/bin/env python3
"""Build a C3-trimer + target complex for the af2_initial_guess specificity test.

Takes a validated per-subunit binder (chain A of an IG-input PDB, which carries the
MPNN sequence in the hexamer x,y frame) and replicates it by exact C3 symmetry about
the Z axis (the hexamer 3-fold; centroid is at x=0,y=0). The three copies become ONE
chain (chain A) so af2_initial_guess treats the whole trimer as the binder
(predict.py: binderlen = first chain). No linker residues are added — the large
inter-subunit gaps become AlphaFold chain-breaks, which is the faithful model for the
binding test. The fused single-chain synthesis construct (with GS linkers) is a
separate, downstream concern.

Target chains are read from the hexamer PDB and translated into the binder frame
(the IG-input frame = hexamer translated by (0,0,-Z_SHIFT)), then relabelled to come
AFTER the binder chain so they are scored as "target".

Output PDB layout:  chain A = trimer (3*Nsub residues) ; chains B,C,... = target chains.

Usage:
  build_trimer_complex.py --subunit-pdb outputs/06_ig/inputs/design_16_s0.pdb \
      --subunit-chain A --target-pdb inputs/1lp3_hexamer_trimmed_fixed.pdb \
      --target-chains ADF --z-shift 168.41 --out outputs/07_specificity/design_16_s0_hex.pdb
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path

C3_ANGLES = (0.0, 120.0, 240.0)


def read_atoms(path: Path, chain: str):
    """Return list of raw ATOM lines for a given chain (ATOM records only)."""
    out = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            out.append(line)
    return out


def rotate_z(x: float, y: float, theta_deg: float):
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return c * x - s * y, s * x + c * y


def set_xyz(line: str, x: float, y: float, z: float) -> str:
    return f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"


def set_chain_res_serial(line: str, chain: str, resnum: int, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:21]}{chain}{resnum:4d}{line[26:]}"


def get_xyz(line: str):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def get_resnum(line: str) -> int:
    return int(line[22:26])


def build(subunit_lines, target_blocks, z_shift):
    """subunit_lines: ATOM lines of one binder subunit (in binder frame).
    target_blocks: list of (orig_chain, [lines]) from hexamer, to translate by -z_shift.
    Returns list of output PDB lines."""
    out = []
    serial = 0

    # --- binder: 3 C3 copies, one continuous chain A, no linker (gaps -> AF2 chainbreaks) ---
    resnums = sorted({get_resnum(l) for l in subunit_lines})
    res_min = resnums[0]
    n_res = len(resnums)
    out_res = 0
    for theta in C3_ANGLES:
        for line in subunit_lines:
            x, y, z = get_xyz(line)
            nx, ny = rotate_z(x, y, theta)
            serial += 1
            new_res = out_res + (get_resnum(line) - res_min) + 1
            nl = set_xyz(line, nx, ny, z)
            nl = set_chain_res_serial(nl, "A", new_res, serial)
            out.append(nl)
        out_res += n_res
    out.append(f"TER   {serial+1:5d}      {subunit_lines[-1][17:20]} A{out_res:4d}")

    # --- target: ALL chains merged into ONE chain B with unique continuous residue
    # numbers. af2_initial_guess (predict.py) supports only 2 chains and splits binder
    # vs target by chain LETTER; the spatially-separate protomers (zero cross-contacts)
    # become AF2 chain-breaks within chain B. Residue numbers must be globally unique or
    # predict.py renumbers/aborts, so we continue numbering past the binder.
    tgt_res = out_res
    last_res = out_res
    for orig_chain, lines in target_blocks:
        prev = None
        for line in lines:
            r = get_resnum(line)
            if r != prev:
                tgt_res += 1
                prev = r
            x, y, z = get_xyz(line)
            serial += 1
            nl = set_xyz(line, x, y, z - z_shift)
            nl = set_chain_res_serial(nl, "B", tgt_res, serial)
            out.append(nl)
        last_res = tgt_res
    out.append(f"TER   {serial+1:5d}      {target_blocks[-1][1][-1][17:20]} B{last_res:4d}")

    out.append("END")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subunit-pdb", required=True)
    ap.add_argument("--subunit-chain", default="A")
    ap.add_argument("--target-pdb", required=True)
    ap.add_argument("--target-chains", required=True, help="e.g. ADF (hex) or AE (dimer)")
    ap.add_argument("--z-shift", type=float, default=168.41,
                    help="z translation to bring hexamer into the binder (IG-input) frame")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    subunit = read_atoms(Path(a.subunit_pdb), a.subunit_chain)
    if not subunit:
        raise SystemExit(f"No chain {a.subunit_chain} atoms in {a.subunit_pdb}")

    target_blocks = []
    for ch in a.target_chains:
        lines = read_atoms(Path(a.target_pdb), ch)
        if not lines:
            raise SystemExit(f"No chain {ch} atoms in {a.target_pdb}")
        target_blocks.append((ch, lines))

    out_lines = build(subunit, target_blocks, a.z_shift)
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text("\n".join(out_lines) + "\n")

    nsub = len({get_resnum(l) for l in subunit})
    print(f"wrote {outp}  (binder 3x{nsub}={3*nsub} res chain A; "
          f"target chains {a.target_chains} merged -> chain B)")


if __name__ == "__main__":
    main()
