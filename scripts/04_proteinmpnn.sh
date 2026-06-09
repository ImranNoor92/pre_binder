#!/usr/bin/env bash
# Phase 3 — ProteinMPNN sequence design on the Phase-2b-validated backbones.
# Designs the binder subunit (chain G) only, with target chains A-F fixed. One subunit is
# designed and reused in all three C3 copies (equivalent to tied positions for a
# replication-built homotrimer — see lib/seqtools.py and the README "Where ... is called").
#
# Input : outputs/02b_af2_validated/validated.txt  (falls back to all Phase-1 designs)
# Output: outputs/03_mpnn_sequences/seqs/<design>.fa  (8 sequences/backbone)
# Usage : [GPU=0] [NUM_SEQ_PER_BACKBONE=8] bash scripts/04_proteinmpnn.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"

GPU="${GPU:-0}"; NSEQ="${NUM_SEQ_PER_BACKBONE:-8}"
VALID="$PROJECT/outputs/02b_af2_validated/validated.txt"
OUTDIR="$PROJECT/outputs/03_mpnn_sequences"

if [ -s "$VALID" ]; then
  mapfile -t NAMES < "$VALID"
  echo "→ Phase 3: MPNN on $(wc -l < "$VALID") validated backbone(s), $NSEQ seqs each."
else
  mapfile -t NAMES < <(find "$PROJECT/outputs/01_rfdiffusion_pilot" -maxdepth 1 -name 'design_*.pdb' -printf '%f\n' | sed 's/\.pdb$//' | sort)
  echo "⚠ No gate list at $VALID — running MPNN on all ${#NAMES[@]} Phase-1 backbones."
fi

for name in "${NAMES[@]}"; do
  [ -n "$name" ] || continue
  bash "$HERE/lib/mpnn_subunit.sh" "$PROJECT/outputs/01_rfdiffusion_pilot/$name.pdb" "$OUTDIR" "$NSEQ" "$GPU"
done
echo "→ Phase 3 complete. FASTAs in $OUTDIR/seqs/"
