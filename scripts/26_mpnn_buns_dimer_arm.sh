#!/usr/bin/env bash
# Stage 1 ROUND 2 — buns-focused ProteinMPNN on the backbones that reached negative Rosetta dG.
# Round-1 MPNN was run plain (no bias) → interfaces were greasy (high buried-unsat polars). Here we
# bias the sequence design AWAY from the usual buried-unsatisfied-polar culprits (Asn/Gln/Ser/Thr)
# while leaving Asp/Glu/Lys/Arg fully available so the Lys42 salt bridge can still form. Deeper
# sampling (50 seq/backbone vs 8) on the ~14 backbones we know dock well. Target A+F held fixed.
#
# Input : backbone list (arg 1, default 06_mpnn_buns_backbones.txt) + 01_rfd/design_*.pdb
# Output: <campaign>/06_mpnn_buns/seqs/design_*.fa
# Usage : [GPU=0] [NSEQ=50] [BIAS=-1.0] bash scripts/26_mpnn_buns_dimer_arm.sh [backbones.txt]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; source "$HERE/lib/common.sh"
GPU="${GPU:-0}"; NSEQ="${NSEQ:-50}"; TEMP="${TEMP:-0.1}"; BCHAIN="${BCHAIN:-B}"; BIAS="${BIAS:--1.0}"
CAMP="$PROJECT/outputs/C3_symmetric_Binder_2026_07_02"
RFDDIR="$CAMP/01_rfd"; OUTDIR="$CAMP/06_mpnn_buns"
BBLIST="${1:-$CAMP/06_mpnn_buns_backbones.txt}"
mkdir -p "$OUTDIR"
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
STAGE="$work/pdbs"; mkdir -p "$STAGE"

# stage only the chosen backbones (parse_multiple_chains reads a whole dir)
n=0; while read -r bb; do [ -z "$bb" ] && continue; cp "$RFDDIR/$bb.pdb" "$STAGE/"; n=$((n+1)); done < "$BBLIST"
echo "→ buns-MPNN over $n backbones, $NSEQ seq each, binder=chain $BCHAIN, target A+F fixed."
echo "  bias: N,Q,S,T = $BIAS (suppress buried-unsat polar culprits); D,E,K,R untouched (keep K42 salt bridge)."

# global AA bias: down-weight the amide/hydroxyl polars most prone to being buried & unsatisfied
"$MPNN_PY" "$MPNN_HELPERS/make_bias_AA.py" \
  --output_path "$work/bias.jsonl" \
  --AA_list "N Q S T" --bias_list "$BIAS $BIAS $BIAS $BIAS"

"$MPNN_PY" "$MPNN_HELPERS/parse_multiple_chains.py" --input_path="$STAGE" --output_path="$work/parsed.jsonl" >/dev/null
"$MPNN_PY" "$MPNN_HELPERS/assign_fixed_chains.py" --input_path="$work/parsed.jsonl" \
  --output_path="$work/assigned.jsonl" --chain_list "$BCHAIN" >/dev/null

CUDA_VISIBLE_DEVICES=$GPU "$MPNN_PY" "$MPNN_RUN" \
  --jsonl_path "$work/parsed.jsonl" --chain_id_jsonl "$work/assigned.jsonl" \
  --bias_AA_jsonl "$work/bias.jsonl" \
  --out_folder "$OUTDIR" --num_seq_per_target "$NSEQ" --sampling_temp "$TEMP" --seed 37 --batch_size 1

NF=$(find "$OUTDIR/seqs" -name 'design_*.fa' 2>/dev/null | wc -l)
echo "→ Done. $NF FASTAs in $OUTDIR/seqs/ ($((n*NSEQ)) sequences total)."
