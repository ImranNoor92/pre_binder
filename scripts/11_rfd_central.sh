#!/usr/bin/env bash
# Phase 1 (new campaign): RFdiffusion — ONE binder that caps the central top of the hexamer,
# contacting the convergent Asp111 cluster of all six monomers (no trimerization).
#
# Strategy change from the retired C3-subunit approach (see context/new_hexamer_model.pdf):
#  - target = NEW full-VP3 hexamer, UNTRIMMED (all 6 chains, residues 1-504). The only
#    preprocessing is flattening the source 6-MODEL file into one model with 6 unique chains
#    (151lp3t3_hexamer_6chain.pdb) so RFdiffusion reads it as an assembly. No residues removed.
#  - hotspots = residue 111 (exposed, convergent at the central axis) of ALL six chains
#  - the binder is expected to sit FLAT on top (perpendicular-ish to the radial chains), capping
#    the central pocket — engaging multiple monomers at once => hexamer-specific by construction.
#  - we deliberately AVOID the old buried/inner hotspots (105-108,113-120) that made the previous
#    binder run parallel into a single monomer's loops.
#  NOTE: the full target is 3024 residues (504x6) — large for RFdiffusion. Watch GPU memory/time.
#
# Output: outputs/11_rfd_central/design_{N}.pdb (+ .trb)
# Usage : [GPU=1] [NUM_DESIGNS=8] [BINDER_MIN=80] [BINDER_MAX=120] [DESIGN_STARTNUM=0] bash scripts/11_rfd_central.sh
set -euo pipefail

PROJECT=/data/binder_software/pre-binder
# Full hexamer reduced to the binding region (all 6 chains, everything within 50A of the central
# site, ~1158 res) — RFdiffusion can't fit the full 3024-res target on a 48GB GPU. The full,
# untrimmed hexamer (151lp3t3_hexamer_6chain.pdb) is retained for downstream validation.
TARGET="$PROJECT/inputs/151lp3t3_hexamer_R50_6chain.pdb"
CONTIG_TARGET=$(cat "$PROJECT/inputs/151lp3t3_hexamer_R50.contig")
OUTDIR="$PROJECT/outputs/11_rfd_central"
LOGDIR="$PROJECT/logs"
RFD_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python
RFD_SCRIPT=/data/rfdiffusion/scripts/run_inference.py
RFD_CKPT="${RFD_CKPT:-/data/rfdiffusion/models/Complex_base_ckpt.pt}"

GPU="${GPU:-1}"
NUM_DESIGNS="${NUM_DESIGNS:-8}"
DESIGN_STARTNUM="${DESIGN_STARTNUM:-0}"
BINDER_MIN="${BINDER_MIN:-80}"
BINDER_MAX="${BINDER_MAX:-120}"
# Central convergent hotspot: Asp111 of all six chains (the only exposed residue reaching the axis-top).
HOTSPOTS="${HOTSPOTS:-[A111,B111,C111,D111,E111,F111]}"

mkdir -p "$OUTDIR" "$LOGDIR"

EXISTING=$(find "$OUTDIR" -maxdepth 1 -name 'design_*.pdb' 2>/dev/null | wc -l)
if [ "$EXISTING" -ge "$((DESIGN_STARTNUM+NUM_DESIGNS))" ]; then
  echo "→ Already have $EXISTING designs in $OUTDIR. Skipping. Delete them or change OUTDIR to re-run."
  exit 0
fi

LOG="$LOGDIR/11_rfd_central_$(date +%Y%m%d_%H%M%S).log"
echo "→ RFdiffusion central-cap run. Log: $LOG"
echo "  GPU=$GPU designs=$NUM_DESIGNS start=$DESIGN_STARTNUM binder=$BINDER_MIN-$BINDER_MAX hotspots=$HOTSPOTS"

CUDA_VISIBLE_DEVICES=$GPU \
"$RFD_PY" "$RFD_SCRIPT" \
  inference.ckpt_override_path="$RFD_CKPT" \
  inference.output_prefix="$OUTDIR/design" \
  inference.input_pdb="$TARGET" \
  inference.num_designs="$NUM_DESIGNS" \
  inference.design_startnum="$DESIGN_STARTNUM" \
  "contigmap.contigs=[${CONTIG_TARGET} ${BINDER_MIN}-${BINDER_MAX}]" \
  "ppi.hotspot_res=$HOTSPOTS" \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0 \
  2>&1 | tee "$LOG"

N=$(find "$OUTDIR" -maxdepth 1 -name 'design_*.pdb' 2>/dev/null | wc -l)
echo
echo "→ Done. $N PDBs in $OUTDIR"
