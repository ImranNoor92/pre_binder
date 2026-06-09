#!/usr/bin/env bash
# Shared paths and tool call-sites for the pre-binder pipeline.
# Sourced by the phase scripts. Single source of truth for "where each tool lives".

# --- Project ---
export PROJECT=/data/binder_software/pre-binder
export TARGET_PDB="$PROJECT/inputs/1lp3_hexamer_trimmed_fixed.pdb"

# --- RFdiffusion (Phase 1) ---
export RFD_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python
export RFD_SCRIPT=/data/rfdiffusion/scripts/run_inference.py
export RFD_CKPT=/data/rfdiffusion/models/Complex_base_ckpt.pt

# --- ProteinMPNN (Phase 3) — runs under the RFdiffusion env (needs torch) ---
export MPNN_DIR=/data/rfdiffusion/external/ProteinMPNN
export MPNN_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python
export MPNN_RUN="$MPNN_DIR/protein_mpnn_run.py"
export MPNN_HELPERS="$MPNN_DIR/helper_scripts"

# --- AlphaFold-multimer (Phases 2b, 4) — DeepMind code on PYTHONPATH via wrapper ---
export AF2_PY=/data/rfdiffusion/.venv-af2/bin/python
export AF2_WRAPPER=/data/rfdiffusion/trials/trial_1/run_alphafold_wrapper.py   # injects /data/alphafold_code/alphafold
export AF2_DB=/data/alphafold_db
export AF2_MAX_TEMPLATE_DATE="${AF2_MAX_TEMPLATE_DATE:-2026-05-08}"

# --- af2_initial_guess validator (dl_binder_design) — single-seq, templated, NO MSA ---
# Set up once with scripts/setup_af2ig.sh. Validates one binder subunit vs one target chain.
export IG_PY=/home/a-mxn833/mambaforge/envs/af2ig/bin/python
export IG_DIR="$PROJECT/external/dl_binder_design/af2_initial_guess"
export IG_PREDICT="$IG_DIR/predict.py"
# dl_binder_design pass thresholds (per the paper): strong interface pae_interaction<10, plddt_binder>80
export IG_PAE_MAX="${IG_PAE_MAX:-10}"
export IG_PLDDT_MIN="${IG_PLDDT_MIN:-80}"

# Binder is always the single non-target chain produced by RFdiffusion.
export BINDER_CHAIN=G
# The target chain each subunit is designed against (Phase 1 hotspots are on chain A).
export TARGET_CHAIN="${TARGET_CHAIN:-A}"
# Inter-subunit linker used by 02_trimerize_replicate.py (GGGGS GGS = 8 aa).
export LINKER_SEQ=GGGGSGGS

mkdir -p "$PROJECT/outputs" "$PROJECT/logs"
