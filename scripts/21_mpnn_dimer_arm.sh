#!/usr/bin/env bash
# Stage 1 Step 2 — ProteinMPNN over the dimer-arm backbones. Designs the binder (chain B),
# target chains A+F held fixed. Batched (one model load). 8 seq/backbone.
# Output: outputs/C3_symmetric_Binder_2026_07_02/02_mpnn/seqs/design_*.fa
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; source "$HERE/lib/common.sh"
GPU="${GPU:-0}"; NSEQ="${NSEQ:-8}"; TEMP="${TEMP:-0.1}"; BCHAIN="${BCHAIN:-B}"
CAMPAIGN="$PROJECT/outputs/C3_symmetric_Binder_2026_07_02"   # one campaign = one folder
INDIR="$CAMPAIGN/01_rfd"; OUTDIR="$CAMPAIGN/02_mpnn"          # backbones in, sequences out
mkdir -p "$OUTDIR"
N=$(find "$INDIR" -maxdepth 1 -name 'design_*.pdb' | wc -l)
echo "-> MPNN (batched) over $N backbones, $NSEQ seq each, binder=chain $BCHAIN, target A+F fixed."
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
$MPNN_PY "$MPNN_HELPERS/parse_multiple_chains.py" --input_path="$INDIR" --output_path="$work/parsed.jsonl" >/dev/null
$MPNN_PY "$MPNN_HELPERS/assign_fixed_chains.py" --input_path="$work/parsed.jsonl" \
  --output_path="$work/assigned.jsonl" --chain_list "$BCHAIN" >/dev/null
CUDA_VISIBLE_DEVICES=$GPU $MPNN_PY "$MPNN_RUN" \
  --jsonl_path "$work/parsed.jsonl" --chain_id_jsonl "$work/assigned.jsonl" \
  --out_folder "$OUTDIR" --num_seq_per_target "$NSEQ" --sampling_temp "$TEMP" --seed 37 --batch_size 1
NF=$(find "$OUTDIR/seqs" -name 'design_*.fa' 2>/dev/null | wc -l)
echo "-> Done. $NF FASTAs in $OUTDIR/seqs/ (expected $N)"
