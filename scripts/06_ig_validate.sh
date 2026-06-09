#!/usr/bin/env bash
# Phase IG — validate each designed subunit against its target chain with af2_initial_guess
# (single-sequence, templated, NO MSA). Replaces the vanilla full-MSA AF2 gate/revalidation
# (which got OOM-killed, was slow, and didn't reproduce de-novo designs).
#
# Flow: for every RFdiffusion backbone x every ProteinMPNN sequence (Phase 3), thread the
# sequence onto the binder backbone, pair it with target chain A, batch-predict with
# af2_initial_guess (one model load for all), then filter+rank by interface quality.
#
# Output: outputs/06_ig/ranked.csv  +  outputs/06_ig/out/<tag>_af2pred.pdb (ranked designs)
# Usage : [GPU=0] [IG_PAE_MAX=10] [IG_PLDDT_MIN=80] bash scripts/06_ig_validate.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"

GPU="${GPU:-0}"
MPNN_DIR="$PROJECT/outputs/03_mpnn_sequences/seqs"
IGOUT="$PROJECT/outputs/06_ig"
INDIR="$IGOUT/inputs"; PREDOUT="$IGOUT/out"
mkdir -p "$INDIR" "$PREDOUT"

[ -d "$MPNN_DIR" ] || { echo "ERROR: no MPNN sequences at $MPNN_DIR. Run scripts/04_proteinmpnn.sh first."; exit 1; }

echo "→ Building threaded binder+target inputs..."
n=0
for fa in "$MPNN_DIR"/design_*.fa; do
  [ -e "$fa" ] || continue
  design=$(basename "$fa" .fa)
  bb="$PROJECT/outputs/01_rfdiffusion_pilot/$design.pdb"
  [ -f "$bb" ] || { echo "  ! missing backbone $bb"; continue; }
  i=0
  while IFS= read -r seq; do
    [ -n "$seq" ] || continue
    out="$INDIR/${design}_s${i}.pdb"
    # target chain comes from the SAME backbone PDB so the binder keeps its designed pose
    # relative to the target (af2_initial_guess uses those coords as the initial guess).
    [ -f "$out" ] || "$IG_PY" "$HERE/lib/thread_and_pair.py" \
      --backbone "$bb" --binder-chain "$BINDER_CHAIN" --seq "$seq" \
      --target-pdb "$bb" --target-chain "$TARGET_CHAIN" --out "$out"
    i=$((i+1)); n=$((n+1))
  done < <("$IG_PY" "$HERE/lib/seqtools.py" mpnn-seqs "$fa")
done
echo "  built $n threaded complexes in $INDIR"

echo "→ Running af2_initial_guess over the batch (one model load)..."
SC="$IGOUT/scores.sc"
cd "$IG_DIR"
CUDA_VISIBLE_DEVICES=$GPU "$IG_PY" predict.py \
  -pdbdir "$INDIR" -outpdbdir "$PREDOUT" -scorefilename "$SC"

echo "→ Ranking..."
"$IG_PY" "$HERE/lib/ig_rank.py" --sc "$SC" --out-csv "$IGOUT/ranked.csv" \
  --pae-max "$IG_PAE_MAX" --plddt-min "$IG_PLDDT_MIN"
echo "→ Done. Ranked designs: $IGOUT/ranked.csv ; structures: $PREDOUT/"
