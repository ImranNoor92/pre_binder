#!/usr/bin/env python3
"""
Pure-stdlib sequence utilities for the pre-binder pipeline (runs under any python).

Subcommands:
  target-seqs   PDB chains A-F  -> per-chain sequences (printed / JSON)
  mpnn-seqs     ProteinMPNN .fa -> the designed binder-subunit sequences (skips native record)
  build-fasta   subunit seq + target PDB + target chains -> AF2-multimer FASTA + layout JSON

The binder is a single fused polypeptide: subunit - LINKER - subunit - LINKER - subunit.
Because the trimer is built by exact C3 replication (02_trimerize_replicate.py), the three
subunits are identical, so one designed subunit sequence is replicated to all three copies.
This is equivalent to ProteinMPNN tied-positions for a replication-built homotrimer.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def chain_sequences(pdb: Path) -> "dict[str, str]":
    """Return {chain: one-letter sequence} in residue order, CA-defined, from ATOM records."""
    seqs: "dict[str, list[tuple[int, str]]]" = {}
    seen: "dict[str, set]" = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain = line[21]
        resnum = int(line[22:26])
        resname = line[17:20].strip()
        key = (chain, resnum)
        seen.setdefault(chain, set())
        if resnum in seen[chain]:
            continue
        seen[chain].add(resnum)
        seqs.setdefault(chain, []).append((resnum, THREE_TO_ONE.get(resname, "X")))
    return {c: "".join(aa for _, aa in sorted(v)) for c, v in seqs.items()}


def cmd_target_seqs(args):
    seqs = chain_sequences(Path(args.pdb))
    chains = list(args.chains) if args.chains else sorted(seqs)
    out = {c: seqs[c] for c in chains if c in seqs}
    if args.json:
        print(json.dumps(out))
    else:
        for c, s in out.items():
            print(f">{c}\n{s}")


def cmd_mpnn_seqs(args):
    """Print designed subunit sequences from a ProteinMPNN .fa (one per line), skipping native."""
    records = []
    cur = None
    for line in Path(args.fasta).read_text().splitlines():
        if line.startswith(">"):
            cur = []
            records.append(("header", line, cur))
        elif cur is not None:
            cur.append(line.strip())
    seqs = ["".join(body) for _, _, body in records]
    # record 0 is the native (poly-Gly) input sequence -> skip it
    designed = seqs[1:]
    if args.max:
        designed = designed[: args.max]
    for s in designed:
        print(s)


def cmd_build_fasta(args):
    subunit = args.subunit_seq.strip()
    linker = args.linker
    binder_seq = (subunit + linker) * 2 + subunit  # 3 subunits, 2 linkers
    target = chain_sequences(Path(args.target_pdb))
    target_chains = list(args.target_chains)

    fasta_lines = []
    layout = {"target_chains": [], "linker": linker, "subunit_len": len(subunit)}
    offset = 0
    for c in target_chains:
        s = target[c]
        fasta_lines.append(f">{c}")
        fasta_lines.append(s)
        layout["target_chains"].append({"chain": c, "len": len(s), "start": offset})
        offset += len(s)
    fasta_lines.append(f">{args.binder_id}")
    fasta_lines.append(binder_seq)
    # per-subunit windows (0-based, within the full concatenated sequence)
    L, Ln = len(subunit), len(linker)
    subunit_starts = [offset, offset + L + Ln, offset + 2 * (L + Ln)]
    layout["binder"] = {
        "chain": args.binder_id, "len": len(binder_seq), "start": offset,
        "subunit_windows": [[s, s + L] for s in subunit_starts],
    }
    Path(args.out).write_text("\n".join(fasta_lines) + "\n")
    Path(args.layout).write_text(json.dumps(layout, indent=2))
    sys.stderr.write(
        f"Wrote {args.out}: {len(target_chains)} target chains + binder "
        f"({len(binder_seq)} aa = 3x{L} + 2x{Ln}); layout -> {args.layout}\n"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("target-seqs")
    t.add_argument("pdb")
    t.add_argument("--chains", default="ABCDEF")
    t.add_argument("--json", action="store_true")
    t.set_defaults(func=cmd_target_seqs)

    m = sub.add_parser("mpnn-seqs")
    m.add_argument("fasta")
    m.add_argument("--max", type=int, default=0)
    m.set_defaults(func=cmd_mpnn_seqs)

    b = sub.add_parser("build-fasta")
    b.add_argument("--subunit-seq", required=True)
    b.add_argument("--target-pdb", required=True)
    b.add_argument("--target-chains", required=True, help="e.g. ABCDEF or AE")
    b.add_argument("--linker", default="GGGGSGGS")
    b.add_argument("--binder-id", default="G")
    b.add_argument("--out", required=True)
    b.add_argument("--layout", required=True)
    b.set_defaults(func=cmd_build_fasta)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
