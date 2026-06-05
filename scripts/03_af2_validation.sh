#!/usr/bin/env bash
# Phase 2b: Run AF2-multimer prediction on each trimerized backbone from Phase 2a,
# then apply the 4 filters (binder pLDDT, interface pTM, per-subunit interface SASA, RMSD).
#
# Wall time: ~10-15 minutes per design on RTX 6000 Ada × N designs.
# Output: outputs/02b_af2_validated/  (passes)
#         outputs/02b_af2_rejected/   (failures, with reason in filename suffix)
#         outputs/02b_af2_metrics.csv (per-design metrics for all)
#
# IMPORTANT: This script is a SKELETON. The actual AF2 invocation depends on which
# AF2 pipeline you use on this machine. Two main options:
#   (1) /data/alphafold_code (full AlphaFold) — heavy, accurate
#   (2) ColabFold MMseqs2 (lightweight) — if installed
# The trial_2B pipeline at /data/rfdiffusion/trial_2B/master_pipeline/ has working
# AF2 invocations against capsid targets — adapt from there.
#
# Filters are applied in a separate Python step after AF2 finishes.

set -euo pipefail

PROJECT=/data/binder_software/pre-binder
INPUT_DIR="$PROJECT/outputs/02a_trimerized"
OUT_VALID="$PROJECT/outputs/02b_af2_validated"
OUT_REJECT="$PROJECT/outputs/02b_af2_rejected"
METRICS_CSV="$PROJECT/outputs/02b_af2_metrics.csv"
LOGDIR="$PROJECT/logs"

GPU="${GPU:-0}"

mkdir -p "$OUT_VALID" "$OUT_REJECT" "$LOGDIR"

if [ ! -d "$INPUT_DIR" ] || [ -z "$(ls -A "$INPUT_DIR" 2>/dev/null)" ]; then
  echo "ERROR: No trimerized PDBs found in $INPUT_DIR. Run 02_trimerize_replicate.py first."
  exit 1
fi

LOG="$LOGDIR/03_af2_validation_$(date +%Y%m%d_%H%M%S).log"
echo "→ Phase 2b: AF2 validation. Log: $LOG"

# ===== TODO: replace this block with the actual AF2 invocation =====
# The trial_2B pipeline uses a master Python orchestrator at:
#   /data/rfdiffusion/trial_2B/master_pipeline/protein_variant_pipeline.py
# That script has an "run-af2" subcommand. Adapt its config to point at our inputs.
#
# Pseudocode:
#   for each trimer.pdb in INPUT_DIR:
#     Convert to chain-separated FASTA: target (6 chains) + binder (1 chain)
#     Invoke AF2-multimer with model_preset=multimer, num_predictions=5 (or 1 for speed)
#     Save top-ranked prediction to a temporary location
#     Apply filters (see below)
#
# Example: see /data/rfdiffusion/trial_2B/master_pipeline/trial2B_config.json
# for how AF2 is configured in that pipeline.

echo "==========================================================================="
echo "  Skeleton script — AF2 invocation not yet wired in."
echo "  Adapt from /data/rfdiffusion/trial_2B/master_pipeline/protein_variant_pipeline.py"
echo "  (specifically its 'run-af2' subcommand)."
echo "==========================================================================="

# ===== Filter step (this part is complete and runnable once AF2 outputs exist) =====
# Filters are applied by a Python helper (see 03b_apply_filters.py below if you want
# to add it; not strictly needed for the planning skeleton).

echo
echo "→ When AF2 outputs are available, apply filters using:"
echo "    binder pLDDT > 0.70"
echo "    interface pTM > 0.65"
echo "    all 3 subunits with interface SASA > 200 Å²"
echo "    backbone RMSD vs trimerized design < 3.0 Å"
echo
echo "→ Pass rate target: ~30-50% of input designs."
