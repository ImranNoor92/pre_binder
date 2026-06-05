#!/usr/bin/env bash
# Phase 3: Run ProteinMPNN on each Phase 2b-validated backbone to design sequences.
# Tie the 3 binder subunits so they remain identical (symmetric homotrimer).
# Lock the target chain sequences.
#
# Wall time: ~30 sec per design × N validated designs.
# Output: outputs/03_mpnn_sequences/  (FASTA per backbone + per-sequence repacked PDBs)

set -euo pipefail

PROJECT=/data/binder_software/pre-binder
INPUT_DIR="$PROJECT/outputs/02b_af2_validated"
OUTDIR="$PROJECT/outputs/03_mpnn_sequences"
LOGDIR="$PROJECT/logs"
MPNN_DIR=/data/rfdiffusion/external/ProteinMPNN
MPNN_PY=/data/rfdiffusion/.venv-rfd-gpu/bin/python   # ProteinMPNN runs from the RFdiffusion env (torch) — NOT .venv-af2

GPU="${GPU:-0}"
NUM_SEQ_PER_BACKBONE="${NUM_SEQ_PER_BACKBONE:-8}"
SAMPLING_TEMP="${SAMPLING_TEMP:-0.1}"

mkdir -p "$OUTDIR" "$LOGDIR"

if [ ! -d "$INPUT_DIR" ]; then
  echo "ERROR: $INPUT_DIR does not exist. Run Phase 2b first."
  exit 1
fi

BACKBONES=("$INPUT_DIR"/*.pdb)
if [ ${#BACKBONES[@]} -eq 0 ] || [ ! -f "${BACKBONES[0]}" ]; then
  echo "ERROR: No PDBs in $INPUT_DIR. Phase 2b filter rejected everything; rerun with looser filters or revisit pilot."
  exit 1
fi

echo "→ Phase 3: ProteinMPNN on ${#BACKBONES[@]} validated backbones."
echo "  Sequences per backbone: $NUM_SEQ_PER_BACKBONE"
echo "  Sampling temp: $SAMPLING_TEMP"
echo "  Tied positions: subunit 1 / 2 / 3 of binder chain (homotrimer constraint)"
echo "  Fixed chains: A, B, C, D, E, F (target sequence locked)"

# Per-backbone loop. For each, build a "tied positions" JSON specifying that
# the three binder subunits are equivalent.

for pdb in "${BACKBONES[@]}"; do
  basename=$(basename "$pdb" .pdb)
  echo
  echo "  → $basename"

  # Generate tied_positions JSON for this backbone
  # The exact format depends on ProteinMPNN's version. Helper script is in MPNN_DIR/helper_scripts/.
  # The script make_tied_positions_dict.py builds this from a comma-separated spec.

  # For our case: binder chain G has residues 1..N (3 subunits + 2 linkers).
  # If each subunit is M residues and linker is 8: subunit 1 = 1..M, linker = M+1..M+8,
  # subunit 2 = M+9..2M+8, linker = 2M+9..2M+16, subunit 3 = 2M+17..3M+16.
  # Tied: position i in subunit 1 ↔ position i in subunit 2 ↔ position i in subunit 3.

  # ===== TODO: build tied_positions.jsonl per backbone =====
  # Adapt the pattern from trial_2B's MPNN step:
  #   /data/rfdiffusion/trial_2B/master_pipeline/protein_variant_pipeline.py
  # → look for the "run-mpnn" or equivalent subcommand and how it constructs the
  #   tied_positions specification.

  # Skeleton invocation:
  CUDA_VISIBLE_DEVICES=$GPU \
  "$MPNN_PY" "$MPNN_DIR/protein_mpnn_run.py" \
    --pdb_path "$pdb" \
    --pdb_path_chains "G" \
    --chain_id_jsonl /tmp/${basename}_chain.jsonl \
    --fixed_positions_jsonl /tmp/${basename}_fixed.jsonl \
    --tied_positions_jsonl /tmp/${basename}_tied.jsonl \
    --num_seq_per_target $NUM_SEQ_PER_BACKBONE \
    --sampling_temp $SAMPLING_TEMP \
    --seed 37 \
    --batch_size 1 \
    --out_folder "$OUTDIR/$basename/" \
    || { echo "  ✗ MPNN failed for $basename"; continue; }

  echo "  ✓ $basename done. Sequences in $OUTDIR/$basename/seqs/"
done

echo
echo "→ Phase 3 complete. Inspect FASTA files in $OUTDIR/<backbone>/seqs/"
