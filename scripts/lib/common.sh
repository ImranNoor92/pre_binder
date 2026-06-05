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

# Binder is always the single non-target chain produced by RFdiffusion.
export BINDER_CHAIN=G
# Inter-subunit linker used by 02_trimerize_replicate.py (GGGGS GGS = 8 aa).
export LINKER_SEQ=GGGGSGGS

mkdir -p "$PROJECT/outputs" "$PROJECT/logs"
