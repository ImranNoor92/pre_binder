#!/usr/bin/env python3
"""Bucket B (interface confidence) — re-predict the top-5-per-group shortlist with ColabDesign/
AF2-multimer (target hexamer templated) and SAVE the PAE matrix + complex structure so we can
compute the dilution-robust interface metrics (ipSAE, pDockQ2) in addition to i_pTM / i_pAE.

Saves per design under outputs/15_colabdesign/bucketB/<tag>/: pae.npy, plddt.npy, complex.pdb
and appends scalars to results_bucketB.jsonl (incremental/resumable). Run with BindCraft python.
"""
from __future__ import annotations
import csv, glob, os, json
import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

ROOT = "/data/binder_software/pre-binder"
SHORT = f"{ROOT}/outputs/14_af2_binderonly/top_af2_binder_only/top5_per_group.csv"
BBDIR = f"{ROOT}/outputs/11_rfd_central"
SEQDIR = f"{ROOT}/outputs/13_mpnn_central/seqs"
PARAMS = "/data/binder_software/BindCraft/params"
OUT = f"{ROOT}/outputs/15_colabdesign/bucketB"
os.makedirs(OUT, exist_ok=True)
JSONL = f"{ROOT}/outputs/15_colabdesign/results_bucketB.jsonl"


def seq_for(tag):
    bb, s = tag.rsplit("_s", 1)
    recs, cur = [], None
    for line in open(f"{SEQDIR}/{bb}.fa"):
        if line.startswith(">"): cur = []; recs.append(cur)
        elif cur is not None: cur.append(line.strip())
    return ["".join(r) for r in recs][1:][int(s)]


def main():
    done = set()
    if os.path.exists(JSONL):
        done = {json.loads(l)["tag"] for l in open(JSONL)}
    inc = open(JSONL, "a", buffering=1)
    for c in csv.DictReader(open(SHORT)):
        tag = c["tag"]
        if tag in done:
            continue
        bb = c["backbone"]; seq = seq_for(tag); L = len(seq)
        d = f"{OUT}/{tag}"; os.makedirs(d, exist_ok=True)
        try:
            clear_mem()
            m = mk_afdesign_model(protocol="binder", use_multimer=True, use_initial_guess=True,
                                  num_recycles=3, data_dir=PARAMS)
            m.prep_inputs(pdb_filename=f"{BBDIR}/{bb}.pdb", target_chain="A,B,C,D,E,F",
                          binder_chain="G", binder_len=L)
            m.predict(seq=seq, models=[0], num_recycles=3, verbose=False)
            aux = m.aux; log = aux["log"]
            np.save(f"{d}/pae.npy", np.asarray(aux["pae"], dtype=np.float32))
            np.save(f"{d}/plddt.npy", np.asarray(aux["plddt"], dtype=np.float32))
            m.save_pdb(f"{d}/complex.pdb")
            r = {"tag": tag, "group": c["group"], "n_monomers": c["n_monomers"],
                 "binder_len": L, "i_ptm": round(float(log["i_ptm"]), 4),
                 "i_pae": round(float(log["i_pae"]), 4), "ptm": round(float(log["ptm"]), 4),
                 "plddt": round(float(log["plddt"]), 4)}
        except Exception as e:
            r = {"tag": tag, "group": c["group"], "error": str(e)[:140]}
        print(json.dumps(r), flush=True); inc.write(json.dumps(r) + "\n")
    inc.close()
    print("done")


if __name__ == "__main__":
    main()
