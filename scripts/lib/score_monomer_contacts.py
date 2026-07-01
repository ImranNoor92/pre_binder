#!/usr/bin/env python3
"""Score how many of the 6 hexamer monomers each RFdiffusion backbone's binder contacts.

The 'binds all 6 monomers' goal is a property of the backbone geometry, not the sequence,
so we measure it directly on the RFd output (binder = chain G, target = chains A-F). For each
design we report: #monomers contacted, per-chain CA contacts, binder centrality (r_xy of the
binder centroid), and whether it caps from the exterior (binder beyond the Asp111 cluster).

Usage: score_monomer_contacts.py --dir outputs/11_rfd_central --out outputs/11_rfd_central/monomer_contacts.csv
"""
from __future__ import annotations
import argparse, csv, glob, math, os
import numpy as np

BINDER_CHAIN = "G"
TARGET_CHAINS = "ABCDEF"
CA_CONTACT = 8.0     # CA-CA; tight enough to mean a real contact patch (10A was uninformative:
                     # a centrally-placed binder is within 10A of all 6 convergent loops trivially)
MIN_PATCH = 5        # >= this many CA<8 contacts to a chain == that chain is "substantially engaged"


def load(path):
    ca = {}
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            ca.setdefault(l[21], {})[int(l[22:26])] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return ca


def score(path):
    ca = load(path)
    if BINDER_CHAIN not in ca:
        return None
    G = np.array(list(ca[BINDER_CHAIN].values()))
    gc = G.mean(0)
    per = {}
    nmon = 0   # chains substantially engaged (>= MIN_PATCH CA<8 contacts)
    for c in TARGET_CHAINS:
        if c not in ca:
            per[c] = 0; continue
        T = np.array(list(ca[c].values()))
        n = int((np.sqrt(((G[:, None] - T[None]) ** 2).sum(-1)) < CA_CONTACT).sum())
        per[c] = n
        if n >= MIN_PATCH:
            nmon += 1
    weakest = min(per[c] for c in TARGET_CHAINS)   # contacts to the least-engaged chain
    # exterior check: binder beyond the Asp111 cluster along the body->cluster axis
    body = np.vstack([np.array(list(ca[c].values())) for c in TARGET_CHAINS if c in ca])
    bc = body.mean(0)
    cl = np.mean([ca[c][111] for c in TARGET_CHAINS if c in ca and 111 in ca[c]], axis=0)
    up = (cl - bc); up = up / (np.linalg.norm(up) + 1e-9)
    exterior = float(np.dot(gc - bc, up)) > float(np.dot(cl - bc, up))
    return {
        "design": os.path.basename(path).replace(".pdb", ""),
        "binder_len": len(G),
        "n_monomers": nmon,
        "weakest_chain": weakest,
        "binder_rxy": round(float(math.hypot(gc[0] - bc[0], gc[1] - bc[1])), 1),
        "exterior_cap": exterior,
        **{f"contacts_{c}": per[c] for c in TARGET_CHAINS},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for p in sorted(glob.glob(os.path.join(a.dir, "design_*.pdb")),
                    key=lambda x: int(x.split("design_")[-1].split(".")[0])):
        r = score(p)
        if r:
            rows.append(r)
    if not rows:
        print("no designs found"); return
    # rank: most chains engaged, then best-balanced (weakest chain), exterior caps, then central
    rows.sort(key=lambda r: (-r["n_monomers"], -r["weakest_chain"], not r["exterior_cap"], r["binder_rxy"]))
    keys = ["design", "binder_len", "n_monomers", "weakest_chain", "exterior_cap", "binder_rxy"] + \
           [f"contacts_{c}" for c in TARGET_CHAINS]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    dist = {}
    for r in rows:
        dist[r["n_monomers"]] = dist.get(r["n_monomers"], 0) + 1
    print(f"scored {len(rows)} designs -> {a.out}")
    print("monomers-contacted distribution:", dict(sorted(dist.items(), reverse=True)))
    print("top 5:")
    for r in rows[:5]:
        print(f"  {r['design']}: {r['n_monomers']}/6 monomers, exterior={r['exterior_cap']}, rxy={r['binder_rxy']}")


if __name__ == "__main__":
    main()
