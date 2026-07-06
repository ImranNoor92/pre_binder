#!/usr/bin/env bash
# Stage 2 — RFdiffusion C3-symmetric motif scaffolding: connect the 3 validated arms into one
# rigid ring. Each arm (96 aa) is held as a FIXED motif at its C3 position on the hexamer patch;
# RFd diffuses a symmetric connecting scaffold (BEFORE/AFTER counts below) that reaches toward the
# central hole and packs the 3 subunits together (olig_contacts potentials → compact rigid core).
# No target in this step → designs must be clash-filtered against the hexamer afterward (script 27
# logic) and checked that the arms stay on the patches.
#
# Input : outputs/.../10_ring_scaffold/arm_motif_c3.pdb (3 arms as chain A, C3 about Z, centered)
# Output: outputs/.../10_ring_scaffold/<run>/ring_*.pdb
# Usage : [GPU=0] [NUM=4] [BEFORE=15] [AFTER=30] [RUN=test] bash scripts/28_rfd_ring_scaffold.sh
set -euo pipefail
PROJECT=/data/binder_software/pre-binder
CAMP="$PROJECT/outputs/C3_symmetric_Binder_2026_07_02"
DIR="$CAMP/10_ring_scaffold"
RFD_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python
RFD_SCRIPT=/data/rfdiffusion/scripts/run_inference.py
CKPT=/data/rfdiffusion/models/Base_epoch8_ckpt.pt      # symmetric-motif scaffolding checkpoint

GPU="${GPU:-0}"; NUM="${NUM:-4}"; BEFORE="${BEFORE:-15}"; AFTER="${AFTER:-30}"; RUN="${RUN:-test}"
OUT="$DIR/$RUN"; mkdir -p "$OUT"
# per-subunit contig: BEFORE diffused / arm motif (96) / AFTER diffused / chainbreak — x3 (C3)
CONTIG="[${BEFORE}/A1-96/${AFTER}/0 ${BEFORE}/A101-196/${AFTER}/0 ${BEFORE}/A201-296/${AFTER}/0]"
echo "→ C3 ring scaffold: $NUM designs, scaffold ${BEFORE}+${AFTER}/subunit, arm motif fixed."
echo "  contig: $CONTIG"

CUDA_VISIBLE_DEVICES=$GPU "$RFD_PY" "$RFD_SCRIPT" \
  inference.symmetry="C3" \
  inference.num_designs="$NUM" \
  inference.output_prefix="$OUT/ring" \
  inference.input_pdb="$DIR/arm_motif_c3.pdb" \
  inference.ckpt_override_path="$CKPT" \
  "contigmap.contigs=$CONTIG" \
  'potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.06"]' \
  potentials.olig_intra_all=True potentials.olig_inter_all=True \
  potentials.guide_scale=2 potentials.guide_decay="quadratic" \
  hydra.run.dir="$DIR/_hydra/\${now:%Y-%m-%d_%H-%M-%S}" \
  2>&1 | tail -20

N=$(find "$OUT" -maxdepth 1 -name 'ring_*.pdb' 2>/dev/null | wc -l)
echo "→ Done. $N ring designs in $OUT"
