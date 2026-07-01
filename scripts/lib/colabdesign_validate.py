#!/usr/bin/env python3
"""Validate a central-binder design against the full hexamer with ColabDesign / AF2-multimer.

Why this and not af2_initial_guess: the multimer model is trained on assemblies and (crucially)
ColabDesign keeps the target as a TEMPLATE, so the multi-protomer hexamer is held in place and
cannot drift (the failure mode that made af2_initial_guess unusable here). Reports interface
metrics: i_pTM (higher better) and i_pAE (lower better), plus binder pLDDT and pTM.

Runs in the BindCraft conda env (colabdesign + jax + multimer weights).
Usage (BindCraft python):
  colabdesign_validate.py --pdb design_148.pdb --target-chains A,B,C,D,E,F --binder-chain G \
      --seq <binder_seq> --params /data/binder_software/BindCraft/params [--recycles 3]
Prints one JSON line of metrics.
"""
from __future__ import annotations
import argparse, json
from colabdesign import mk_afdesign_model, clear_mem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--target-chains", required=True)   # e.g. A,B,C,D,E,F
    ap.add_argument("--binder-chain", default="G")
    ap.add_argument("--seq", required=True)
    ap.add_argument("--params", default="/data/binder_software/BindCraft/params")
    ap.add_argument("--recycles", type=int, default=3)
    ap.add_argument("--model", type=int, default=0)      # multimer_v3 model index
    a = ap.parse_args()

    clear_mem()
    m = mk_afdesign_model(protocol="binder", use_multimer=True, use_initial_guess=True,
                          num_recycles=a.recycles, data_dir=a.params)
    # target chains templated (held rigid); binder chain provides the initial-guess backbone
    m.prep_inputs(pdb_filename=a.pdb, target_chain=a.target_chains,
                  binder_chain=a.binder_chain, binder_len=len(a.seq.strip()))
    m.predict(seq=a.seq.strip(), models=[a.model], num_recycles=a.recycles, verbose=False)
    log = m.aux["log"]
    out = {k: round(float(log[k]), 4) for k in ("plddt", "ptm", "i_ptm", "pae", "i_pae") if k in log}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
