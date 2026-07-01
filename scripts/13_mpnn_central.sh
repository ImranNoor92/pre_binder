#!/usr/bin/env bash
# Phase 3 (central campaign) — ProteinMPNN sequence design over ALL central-binder backbones.
# Designs the binder (chain G) only, fixing target chains A-F. Batched: all backbones parsed
# into one jsonl and run in a single model load (fast), one FASTA per backbone out.
#
# Input : outputs/11_rfd_central/design_*.pdb
# Output: outputs/13_mpnn_central/seqs/design_*.fa  (NSEQ designs each)
# Usage : [GPU=0] [NSEQ=8] [TEMP=0.1] bash scripts/13_mpnn_central.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"

GPU="${GPU:-0}"; NSEQ="${NSEQ:-8}"; TEMP="${TEMP:-0.1}"
INDIR="$PROJECT/outputs/11_rfd_central"
OUTDIR="$PROJECT/outputs/13_mpnn_central"
mkdir -p "$OUTDIR"

N=$(find "$INDIR" -maxdepth 1 -name 'design_*.pdb' | wc -l)
echo "→ MPNN (batched) over $N backbones, $NSEQ seqs each, binder=chain $BINDER_CHAIN, target A-F fixed."

work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
$MPNN_PY "$MPNN_HELPERS/parse_multiple_chains.py" \
  --input_path="$INDIR" --output_path="$work/parsed.jsonl" >/dev/null
$MPNN_PY "$MPNN_HELPERS/assign_fixed_chains.py" \
  --input_path="$work/parsed.jsonl" --output_path="$work/assigned.jsonl" \
  --chain_list "$BINDER_CHAIN" >/dev/null
CUDA_VISIBLE_DEVICES=$GPU $MPNN_PY "$MPNN_RUN" \
  --jsonl_path "$work/parsed.jsonl" \
  --chain_id_jsonl "$work/assigned.jsonl" \
  --out_folder "$OUTDIR" \
  --num_seq_per_target "$NSEQ" \
  --sampling_temp "$TEMP" --seed 37 --batch_size 1

NF=$(find "$OUTDIR/seqs" -name 'design_*.fa' 2>/dev/null | wc -l)
echo "→ Done. $NF FASTAs in $OUTDIR/seqs/ (expected $N)"
