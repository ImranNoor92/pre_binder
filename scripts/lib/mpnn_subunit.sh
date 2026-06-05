#!/usr/bin/env bash
# Design the binder subunit (chain G) of one RFdiffusion backbone with ProteinMPNN,
# fixing the target chains A-F. One subunit only — the C3 trimer reuses this sequence
# in all three copies (equivalent to tied positions for a replication-built homotrimer).
#
# Usage: mpnn_subunit.sh <backbone.pdb> <out_dir> [num_seq=8] [gpu=0] [temp=0.1]
# Output: <out_dir>/seqs/<name>.fa   (record 0 = native poly-Gly, records 1..N = designs)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/common.sh"

PDB="$1"; OUTDIR="$2"; NSEQ="${3:-8}"; GPU="${4:-0}"; TEMP="${5:-0.1}"
name="$(basename "$PDB" .pdb)"
mkdir -p "$OUTDIR"

if [ -f "$OUTDIR/seqs/$name.fa" ]; then
  echo "  → MPNN already done for $name (skipping)"; exit 0
fi

work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
mkdir -p "$work/pdbs"; cp "$PDB" "$work/pdbs/"

$MPNN_PY "$MPNN_HELPERS/parse_multiple_chains.py" \
  --input_path="$work/pdbs" --output_path="$work/parsed.jsonl" >/dev/null
$MPNN_PY "$MPNN_HELPERS/assign_fixed_chains.py" \
  --input_path="$work/parsed.jsonl" --output_path="$work/assigned.jsonl" \
  --chain_list "$BINDER_CHAIN" >/dev/null
CUDA_VISIBLE_DEVICES=$GPU $MPNN_PY "$MPNN_RUN" \
  --jsonl_path "$work/parsed.jsonl" \
  --chain_id_jsonl "$work/assigned.jsonl" \
  --out_folder "$OUTDIR" \
  --num_seq_per_target "$NSEQ" \
  --sampling_temp "$TEMP" --seed 37 --batch_size 1 >/dev/null
echo "  → MPNN $name: $NSEQ sequences -> $OUTDIR/seqs/$name.fa"
