#!/usr/bin/env python3
"""Run ColabDesign/AF2-multimer binding validation over the top-5-per-group shortlist.

For each candidate (design_N_si): target hexamer chains A-F templated (held), binder chain G
backbone as initial guess, MPNN sequence set; report i_pTM / i_pAE / pLDDT / pTM. Ranks by i_pTM.
Run with the BindCraft python (colabdesign env).
"""
from __future__ import annotations
import csv, glob, os, json
from colabdesign import mk_afdesign_model, clear_mem

ROOT = "/data/binder_software/pre-binder"
SHORTLIST = f"{ROOT}/outputs/14_af2_binderonly/top_af2_binder_only/top5_per_group.csv"
BBDIR = f"{ROOT}/outputs/11_rfd_central"
SEQDIR = f"{ROOT}/outputs/13_mpnn_central/seqs"
PARAMS = "/data/binder_software/BindCraft/params"
OUT = f"{ROOT}/outputs/15_colabdesign"
os.makedirs(OUT, exist_ok=True)


def seq_for(tag):
    bb, s = tag.rsplit("_s", 1)
    idx = int(s)
    recs, cur = [], None
    for line in open(f"{SEQDIR}/{bb}.fa"):
        if line.startswith(">"):
            cur = []; recs.append(cur)
        elif cur is not None:
            cur.append(line.strip())
    seqs = ["".join(r) for r in recs][1:]   # skip native
    return seqs[idx]


def main():
    cands = list(csv.DictReader(open(SHORTLIST)))
    rows = []
    # incremental log so a kill never loses completed predictions
    inc = open(f"{OUT}/results_incremental.jsonl", "a", buffering=1)
    done = set()
    try:
        for line in open(f"{OUT}/results_incremental.jsonl"):
            done.add(json.loads(line)["tag"])
    except FileNotFoundError:
        pass
    for c in cands:
        if c["tag"] in done:
            continue
        tag = c["tag"]; bb = c["backbone"]
        pdb = f"{BBDIR}/{bb}.pdb"
        seq = seq_for(tag)
        try:
            clear_mem()
            m = mk_afdesign_model(protocol="binder", use_multimer=True, use_initial_guess=True,
                                  num_recycles=3, data_dir=PARAMS)
            m.prep_inputs(pdb_filename=pdb, target_chain="A,B,C,D,E,F",
                          binder_chain="G", binder_len=len(seq))
            m.predict(seq=seq, models=[0], num_recycles=3, verbose=False)
            log = m.aux["log"]
            r = {"tag": tag, "group": c["group"], "n_monomers": c["n_monomers"],
                 "i_ptm": round(float(log["i_ptm"]), 4), "i_pae": round(float(log["i_pae"]), 4),
                 "plddt": round(float(log["plddt"]), 4), "ptm": round(float(log["ptm"]), 4)}
        except Exception as e:
            r = {"tag": tag, "group": c["group"], "n_monomers": c["n_monomers"],
                 "i_ptm": "", "i_pae": "", "plddt": "", "ptm": "", "error": str(e)[:120]}
        print(json.dumps(r), flush=True); inc.write(json.dumps(r) + "\n"); rows.append(r)
    inc.close()
    # merge any previously-completed rows from the incremental log
    seen = {r["tag"] for r in rows}
    for line in open(f"{OUT}/results_incremental.jsonl"):
        d = json.loads(line)
        if d["tag"] not in seen:
            rows.append(d); seen.add(d["tag"])
    rows.sort(key=lambda x: -(x["i_ptm"] if isinstance(x["i_ptm"], float) else -1))
    keys = ["rank", "tag", "group", "n_monomers", "i_ptm", "i_pae", "plddt", "ptm"]
    with open(f"{OUT}/colabdesign_ranked.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for i, r in enumerate(rows, 1):
            w.writerow({"rank": i, **r})
    print(f"\nwrote {OUT}/colabdesign_ranked.csv")


if __name__ == "__main__":
    main()
