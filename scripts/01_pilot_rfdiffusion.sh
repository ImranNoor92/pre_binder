#!/usr/bin/env bash
# Phase 1: RFdiffusion pilot — generate 10 single-subunit binder backbones against chain A of the hexamer.
# Wall time: ~1-2.5 hours on one RTX 6000 Ada.
# Output: outputs/01_rfdiffusion_pilot/design_{0..9}.pdb (+ .trb metadata)

set -euo pipefail

PROJECT=/data/binder_software/pre-binder
TARGET=/data/binder_software/pre-binder/inputs/1lp3_hexamer_trimmed_fixed.pdb
OUTDIR="$PROJECT/outputs/01_rfdiffusion_pilot"
LOGDIR="$PROJECT/logs"
RFD_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python    # RFdiffusion env (torch/e3nn/dgl) — NOT .venv-af2
RFD_SCRIPT=/data/rfdiffusion/scripts/run_inference.py
# rfdiffusion is a non-editable pip install, so its built-in relative model lookup
# (site-packages/models) is empty. Point it at the real weights. Complex_base is the
# checkpoint RFdiffusion auto-selects for PPI/complex design with hotspots.
RFD_CKPT="${RFD_CKPT:-/data/rfdiffusion/models/Complex_base_ckpt.pt}"

GPU="${GPU:-1}"            # override with GPU=0 ./01_pilot_rfdiffusion.sh
NUM_DESIGNS="${NUM_DESIGNS:-10}"
BINDER_MIN="${BINDER_MIN:-60}"
BINDER_MAX="${BINDER_MAX:-90}"

mkdir -p "$OUTDIR" "$LOGDIR"

# Idempotency: skip if pilot is already complete (10 PDBs present)
# (find, not ls glob: ls exits non-zero on no match, which trips set -e/pipefail on the first run)
EXISTING=$(find "$OUTDIR" -maxdepth 1 -name 'design_*.pdb' 2>/dev/null | wc -l)
if [ "$EXISTING" -ge "$NUM_DESIGNS" ]; then
  echo "→ Already have $EXISTING designs in $OUTDIR. Skipping. Delete them or change OUTDIR to re-run."
  exit 0
fi

LOG="$LOGDIR/01_pilot_rfdiffusion_$(date +%Y%m%d_%H%M%S).log"
echo "→ Launching RFdiffusion pilot. Log: $LOG"
echo "  GPU=$GPU  designs=$NUM_DESIGNS  binder length=$BINDER_MIN-$BINDER_MAX"

CUDA_VISIBLE_DEVICES=$GPU \
"$RFD_PY" "$RFD_SCRIPT" \
  inference.ckpt_override_path="$RFD_CKPT" \
  inference.output_prefix="$OUTDIR/design" \
  inference.input_pdb="$TARGET" \
  inference.num_designs="$NUM_DESIGNS" \
  "contigmap.contigs=[A70-150/0 B70-150/0 C70-150/0 D70-150/0 E70-150/0 F70-150/0 ${BINDER_MIN}-${BINDER_MAX}]" \
  "ppi.hotspot_res=[A105,A107,A109,A111,A114,A115]" \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0 \
  2>&1 | tee "$LOG"

# Post-flight check
N=$(find "$OUTDIR" -maxdepth 1 -name 'design_*.pdb' 2>/dev/null | wc -l)
echo
echo "→ Done. $N PDBs in $OUTDIR"
if [ "$N" -lt "$NUM_DESIGNS" ]; then
  echo "⚠ Expected $NUM_DESIGNS designs, got $N. Check $LOG for errors."
  exit 1
fi
