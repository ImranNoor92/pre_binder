#!/usr/bin/env python3
"""Parse an af2_initial_guess scores.sc, apply pass filters, rank, write CSV.
Pass = pae_interaction < pae_max AND plddt_binder > plddt_min (dl_binder_design convention).
Ranked best-first by pae_interaction (lower = stronger interface)."""
from __future__ import annotations
import argparse, csv, re
from pathlib import Path


def parse_sc(path: Path):
    rows, header = [], None
    for line in path.read_text().splitlines():
        if not line.startswith("SCORE:"):
            continue
        toks = line.split()[1:]
        if header is None and "description" in toks:
            header = toks; continue
        if header and len(toks) == len(header):
            rows.append(dict(zip(header, toks)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sc", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--pae-max", type=float, default=10.0)
    ap.add_argument("--plddt-min", type=float, default=80.0)
    a = ap.parse_args()

    rows = parse_sc(Path(a.sc))
    out = []
    for r in rows:
        tag = re.sub(r"_af2pred$", "", r.get("description", ""))
        m = re.match(r"(design_\d+)_s(\d+)", tag)
        design, seq = (m.group(1), int(m.group(2))) if m else (tag, -1)
        pae = float(r["pae_interaction"]); plddt = float(r["plddt_binder"])
        out.append({
            "design": design, "seq_idx": seq, "tag": tag,
            "pae_interaction": round(pae, 2), "plddt_binder": round(plddt, 2),
            "binder_rmsd": round(float(r.get("binder_aligned_rmsd", "nan")), 2),
            "pass": bool(pae < a.pae_max and plddt > a.plddt_min),
        })
    out.sort(key=lambda x: (not x["pass"], x["pae_interaction"]))  # passers first, then best pae
    with open(a.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rank", "design", "seq_idx", "tag",
                                           "pae_interaction", "plddt_binder", "binder_rmsd", "pass"])
        w.writeheader()
        for i, r in enumerate(out, 1):
            w.writerow({"rank": i, **r})
    npass = sum(r["pass"] for r in out)
    print(f"{npass}/{len(out)} designs pass (pae_interaction<{a.pae_max}, plddt_binder>{a.plddt_min}) -> {a.out_csv}")
    for r in out[:5]:
        print(f"  rank {out.index(r)+1}: {r['tag']}  pae_i={r['pae_interaction']} plddt_b={r['plddt_binder']} pass={r['pass']}")


if __name__ == "__main__":
    main()
