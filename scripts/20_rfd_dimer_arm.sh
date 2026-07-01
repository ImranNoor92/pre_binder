#!/usr/bin/env bash
# STAGE 1 — one binder "arm" against the A+F dimer patch (residues 38-42 + 155).
# Standard RFdiffusion binder design (NO symmetry) — the reliable mode. Purpose: prove the
# 38/40/42/155 patch is bindable at all before attempting the C3 ring (Stage 2).
# Patch chemistry: Leu40 + Met155 (hydrophobic anchors), Tyr41 (aromatic), Lys42 (salt-bridge
# handle), Asn38/Gln (polar) -> a mixed, designable interface (unlike the failed central Asp site).
#
# Output: outputs/20_dimer_arm/design_{N}.pdb (+ .trb)
# Usage : [GPU=0] [NUM_DESIGNS=20] [BINDER_MIN=90] [BINDER_MAX=150] [DESIGN_STARTNUM=0] bash scripts/20_rfd_dimer_arm.sh
set -euo pipefail
PROJECT=/data/binder_software/pre-binder
TARGET="$PROJECT/inputs/151lp3t3_dimerAF.pdb"          # A+F dimer (the 28 A hotspot-patch pair)
OUTDIR="$PROJECT/outputs/20_dimer_arm"
LOGDIR="$PROJECT/logs"
RFD_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python
RFD_SCRIPT=/data/rfdiffusion/scripts/run_inference.py
RFD_CKPT="${RFD_CKPT:-/data/rfdiffusion/models/Complex_base_ckpt.pt}"

GPU="${GPU:-0}"; NUM_DESIGNS="${NUM_DESIGNS:-20}"; DESIGN_STARTNUM="${DESIGN_STARTNUM:-0}"
BINDER_MIN="${BINDER_MIN:-90}"; BINDER_MAX="${BINDER_MAX:-150}"
# hotspots = hydrophobic anchors (Leu40, Met155) + salt-bridge handle (Lys42) on BOTH dimer chains
HOTSPOTS="${HOTSPOTS:-[A40,A42,A155,F40,F42,F155]}"

mkdir -p "$OUTDIR" "$LOGDIR"
LOG="$LOGDIR/20_rfd_dimer_arm_$(date +%Y%m%d_%H%M%S).log"
echo "-> RFdiffusion Stage-1 arm vs A+F dimer. Log: $LOG"
echo "   GPU=$GPU designs=$NUM_DESIGNS binder=$BINDER_MIN-$BINDER_MAX hotspots=$HOTSPOTS"

CUDA_VISIBLE_DEVICES=$GPU \
"$RFD_PY" "$RFD_SCRIPT" \
  inference.ckpt_override_path="$RFD_CKPT" \
  inference.output_prefix="$OUTDIR/design" \
  inference.input_pdb="$TARGET" \
  inference.num_designs="$NUM_DESIGNS" \
  inference.design_startnum="$DESIGN_STARTNUM" \
  "contigmap.contigs=[A1-504/0 F1-504/0 ${BINDER_MIN}-${BINDER_MAX}]" \
  "ppi.hotspot_res=$HOTSPOTS" \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
  2>&1 | tee "$LOG"

N=$(find "$OUTDIR" -maxdepth 1 -name 'design_*.pdb' 2>/dev/null | wc -l)
echo "-> Done. $N backbones in $OUTDIR"
