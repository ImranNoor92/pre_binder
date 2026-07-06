#!/usr/bin/env bash
# Stage 2 — RFd C3 symmetric motif scaffolding WITH capsid context (target-aware workaround).
# Per the RFd README: symmetric motif scaffolding diffuses scaffold AROUND the fixed (symmetrized)
# motif and propagates the asymmetric unit by canonical C3 (Z axis). We fold the capsid INTO the
# fixed motif so RFd sees it and won't grow scaffold into the capsid; C3 closes the ring.
#
# Asymmetric unit (untested multi-chain): capsid context (D, 214 res, fixed) + binder (BEFORE
# scaffold / arm A1-96 fixed / AFTER scaffold). x3 by C3. Input already centered on Z at origin.
#
# Input : outputs/.../10_ring_scaffold/sym_ctx_input.pdb  (arms A/B/C, capsid ctx D/E/F)
# Output: outputs/.../10_ring_scaffold/<run>/ring_*.pdb
# Usage : [GPU=0] [NUM=2] [BEFORE=15] [AFTER=30] [RUN=ctx_test] bash scripts/29_rfd_ring_ctx.sh
set -euo pipefail
PROJECT=/data/binder_software/pre-binder
CAMP="$PROJECT/outputs/C3_symmetric_Binder_2026_07_02"
DIR="$CAMP/10_ring_scaffold"
RFD_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python
RFD_SCRIPT=/data/rfdiffusion/scripts/run_inference.py
CKPT=/data/rfdiffusion/models/Base_epoch8_ckpt.pt

GPU="${GPU:-0}"; NUM="${NUM:-2}"; BEFORE="${BEFORE:-15}"; AFTER="${AFTER:-30}"; RUN="${RUN:-ctx_test}"
STARTNUM="${STARTNUM:-0}"; CTXN="${CTXN:-214}"; INPUT="${INPUT:-$DIR/sym_ctx_input.pdb}"
OUT="$DIR/$RUN"; mkdir -p "$OUT"
# per C3 unit: capsid context (D1-CTXN, fixed) + chainbreak + [BEFORE scaffold / arm / AFTER scaffold] + chainbreak
CONTIG="[D1-${CTXN}/0 ${BEFORE}/A1-96/${AFTER}/0 E1-${CTXN}/0 ${BEFORE}/B1-96/${AFTER}/0 F1-${CTXN}/0 ${BEFORE}/C1-96/${AFTER}/0]"
echo "→ C3 ring + capsid context: $NUM designs, scaffold ${BEFORE}+${AFTER}/subunit."
echo "  contig: $CONTIG"

# USE_OLIG=1 keeps the olig_contacts potential (pulls subunits toward their centre — only correct
# when the ring centre is solvent). For the capsid ring the centre is INTERIOR, so default OFF.
USE_OLIG="${USE_OLIG:-0}"
POT=()
if [ "$USE_OLIG" = "1" ]; then
  POT=('potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.06"]'
       potentials.olig_intra_all=True potentials.olig_inter_all=True
       potentials.guide_scale=2 potentials.guide_decay="quadratic")
fi
CUDA_VISIBLE_DEVICES=$GPU "$RFD_PY" "$RFD_SCRIPT" \
  inference.symmetry="C3" \
  inference.model_only_neighbors=True \
  inference.num_designs="$NUM" \
  inference.design_startnum="$STARTNUM" \
  inference.output_prefix="$OUT/ring" \
  inference.input_pdb="$INPUT" \
  inference.ckpt_override_path="$CKPT" \
  "contigmap.contigs=$CONTIG" \
  "${POT[@]}" \
  hydra.run.dir="$DIR/_hydra/\${now:%Y-%m-%d_%H-%M-%S}" \
  2>&1 | tail -30

N=$(find "$OUT" -maxdepth 1 -name 'ring_*.pdb' 2>/dev/null | wc -l)
echo "→ Done. $N designs in $OUT"
