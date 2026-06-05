#!/usr/bin/env bash
# Phase 4: Re-run AF2-multimer on each MPNN-designed sequence and re-apply the 4 filters.
# Final survivors go to outputs/04_final_ranked/, sorted by combined score.
#
# Wall time: ~10-15 min per sequence × N sequences.
# Output: outputs/04_final_ranked/

set -euo pipefail

PROJECT=/data/binder_software/pre-binder
INPUT_DIR="$PROJECT/outputs/03_mpnn_sequences"
OUT_FINAL="$PROJECT/outputs/04_final_ranked"
LOGDIR="$PROJECT/logs"

GPU="${GPU:-0}"

mkdir -p "$OUT_FINAL" "$LOGDIR"

if [ ! -d "$INPUT_DIR" ]; then
  echo "ERROR: $INPUT_DIR does not exist. Run Phase 3 first."
  exit 1
fi

echo "→ Phase 4: AF2 re-validation on MPNN-designed sequences."
echo "  Input: $INPUT_DIR"
echo "  Output: $OUT_FINAL"

# ===== TODO: same as 03_af2_validation.sh =====
# Use the same AF2 invocation as Phase 2b, adapted to read MPNN sequences from the
# FASTA files in $INPUT_DIR/<backbone>/seqs/ rather than from the trimerized PDB.
#
# For each backbone/sequence:
#   1. Combine target hexamer sequences (chains A-F, locked) with the MPNN binder sequence
#   2. Run AF2-multimer
#   3. Apply the same 4 filters (binder pLDDT > 0.70, i_pTM > 0.65, etc.)
#   4. Compute combined score: i_pTM * pLDDT * (1/max(RMSD, 0.5)) * (sum_iSASA/1000)
#   5. Write top-passing designs to $OUT_FINAL with rank prefix

echo "==========================================================================="
echo "  Skeleton script — AF2 invocation not yet wired in."
echo "  Reuses the AF2 setup from 03_af2_validation.sh."
echo "==========================================================================="

echo
echo "→ ACID TEST (the most important validation):"
echo "  After Phase 4 survivors are written, run an additional dimer-only check:"
echo "  For each final design, re-predict AF2 against ONLY chains A+E (one dimer pair)"
echo "  and confirm interface pTM drops by ≥ 0.15 compared to the full hexamer prediction."
echo "  If a design binds the dimer just as well, it's NOT hexamer-specific — discard it."
