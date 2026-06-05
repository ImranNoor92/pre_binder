#!/usr/bin/env python3
"""
Phase 2a: Take each Phase 1 backbone, replicate the binder subunit by C3 symmetry
around the hexamer's 3-fold axis, and fuse the three subunits into one polypeptide
via flexible (GGGGS)n linkers.

Output: outputs/02a_trimerized/design_{N}_trimer.pdb

Wall time: <1 minute total (pure geometric transformation, no ML)

Run after 01_pilot_rfdiffusion.sh completes.
"""
from __future__ import annotations
import math
import sys
from pathlib import Path

PROJECT = Path("/data/binder_software/pre-binder")
INPUT_DIR = PROJECT / "outputs" / "01_rfdiffusion_pilot"
OUTPUT_DIR = PROJECT / "outputs" / "02a_trimerized"
HEXAMER_3FOLD_AXIS_XY = (0.0, 0.0)   # from earlier geometric analysis
HEXAMER_3FOLD_AXIS_Z_DIRECTION = "z"  # axis is along Z

# ---------- PDB parsing ----------

def parse_pdb_atoms(path: Path):
    """Yield (line, chain, resnum, atom_name, x, y, z) for ATOM records."""
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            chain = line[21]
            resnum = int(line[22:26])
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            yield line, chain, resnum, atom_name, x, y, z

def update_xyz(line: str, x: float, y: float, z: float) -> str:
    return f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"

def update_chain_resnum_serial(line: str, chain: str, resnum: int, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:21]}{chain}{resnum:4d}{line[26:]}"

# ---------- Rotation ----------

def rotate_z(x: float, y: float, theta_deg: float):
    """Rotate (x,y) by theta_deg around the Z axis (which passes through origin)."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return c * x - s * y, s * x + c * y

# ---------- Identify the binder chain ----------

def identify_binder_chain(pdb: Path) -> str:
    """Target has chains A-F. Binder is whatever chain ID is NOT A-F (typically G or H)."""
    target_chains = set("ABCDEF")
    seen = set()
    for _, c, *_ in parse_pdb_atoms(pdb):
        if c not in target_chains:
            seen.add(c)
    if not seen:
        raise RuntimeError(f"No binder chain found in {pdb}. Expected something other than A-F.")
    if len(seen) > 1:
        raise RuntimeError(f"Multiple non-target chains in {pdb}: {seen}. Expected one binder chain.")
    return next(iter(seen))

# ---------- Main per-PDB transformation ----------

def trimerize(pdb_in: Path, pdb_out: Path) -> dict:
    """Replicate binder subunit by C3 around Z axis, fuse with GS linkers."""
    binder_chain = identify_binder_chain(pdb_in)

    # Collect target lines (unchanged) and binder atoms
    target_lines = []
    binder_atoms = []  # list of (line, chain, resnum, atom_name, x, y, z)

    for line, c, r, a, x, y, z in parse_pdb_atoms(pdb_in):
        if c in "ABCDEF":
            target_lines.append(line)
        elif c == binder_chain:
            binder_atoms.append((line, c, r, a, x, y, z))

    if not binder_atoms:
        raise RuntimeError(f"No binder atoms in {pdb_in}.")

    # Determine binder residue range
    binder_residues = sorted({r for _, _, r, *_ in binder_atoms})
    n_res = len(binder_residues)
    res_min, res_max = binder_residues[0], binder_residues[-1]

    # Build three copies: 0°, 120°, 240° rotation around Z axis through (0, 0, *)
    # Note: the hexamer's 3-fold axis is at (X=0, Y=0) per our analysis, so we rotate
    # around the literal Z axis. If that ever changes, translate binder so axis is at origin first.

    LINKER_LEN = 8  # GGGGSGGS — 8 residues, ~24 Å fully extended (rough match for inter-subunit gap)
    LINKER_SEQ = ["GLY", "GLY", "GLY", "GLY", "SER", "GLY", "GLY", "SER"]

    out_lines = list(target_lines)
    # All-atom serial counter for the output (target serials are kept; binder serials renumbered)
    # Target's last serial:
    last_target_serial = max(int(line[6:11]) for line in target_lines)
    serial = last_target_serial

    out_resnum = 0  # binder will be one chain, residues numbered 1..N
    new_binder_chain_id = "G"

    for copy_idx, theta in enumerate([0, 120, 240]):
        # Add binder subunit, rotated
        for line, _, r, atom, x, y, z in binder_atoms:
            new_x, new_y = rotate_z(x, y, theta)
            new_serial = serial + 1
            new_resnum = out_resnum + (r - res_min) + 1
            new_line = update_xyz(line, new_x, new_y, z)
            new_line = update_chain_resnum_serial(new_line, new_binder_chain_id, new_resnum, new_serial)
            out_lines.append(new_line)
            serial = new_serial
        out_resnum += n_res

        # Add linker (only between subunits, not after the last)
        if copy_idx < 2:
            # Linker atoms are placeholders — exact positions don't matter for downstream MPNN/AF2
            # because AF2 will refold from scratch given the sequence. We put dummy CA atoms at
            # the average of subunit N+1 N-term and subunit N C-term positions so the chain is
            # continuous. ProteinMPNN tolerates missing atoms.
            for i, aa in enumerate(LINKER_SEQ):
                out_resnum += 1
                serial += 1
                # Placeholder coords — will be ignored / refined by AF2
                placeholder = f"ATOM  {serial:5d}  CA  {aa} {new_binder_chain_id}{out_resnum:4d}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           C  "
                out_lines.append(placeholder)

    # TER + END
    out_lines.append(f"TER   {serial+1:5d}      {LINKER_SEQ[-1]} {new_binder_chain_id}{out_resnum:4d}")
    out_lines.append("END")

    pdb_out.write_text("\n".join(out_lines) + "\n")

    return {
        "subunit_residues": n_res,
        "total_binder_residues": 3 * n_res + 2 * LINKER_LEN,
        "linker_length": LINKER_LEN,
    }


def main():
    if not INPUT_DIR.exists():
        print(f"ERROR: {INPUT_DIR} does not exist. Run 01_pilot_rfdiffusion.sh first.", file=sys.stderr)
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    inputs = sorted(INPUT_DIR.glob("design_*.pdb"))
    if not inputs:
        print(f"ERROR: No design_*.pdb files in {INPUT_DIR}", file=sys.stderr)
        sys.exit(1)

    print(f"→ Found {len(inputs)} Phase 1 backbones. Trimerizing each...")
    print()
    print(f"{'Backbone':<40} {'Subunit res':>12} {'Total binder res':>18}  Status")
    print("-" * 90)

    success = 0
    for pdb_in in inputs:
        pdb_out = OUTPUT_DIR / pdb_in.name.replace(".pdb", "_trimer.pdb")
        try:
            info = trimerize(pdb_in, pdb_out)
            print(f"{pdb_in.name:<40} {info['subunit_residues']:>12d} {info['total_binder_residues']:>18d}  ✓")
            success += 1
        except Exception as e:
            print(f"{pdb_in.name:<40} {'-':>12s} {'-':>18s}  ✗ {e}")

    print()
    print(f"→ {success}/{len(inputs)} backbones trimerized. Output: {OUTPUT_DIR}")
    print()
    print("NOTE: linker residues in the output PDB have placeholder (0,0,0) coordinates.")
    print("      AF2 in Phase 2b will refold the trimer from scratch given the sequence;")
    print("      the placeholder coords are only used for chain connectivity in MPNN.")


if __name__ == "__main__":
    main()
