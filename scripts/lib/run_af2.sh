#!/usr/bin/env bash
# Run AlphaFold-multimer on one FASTA, reusing the proven trial_1 install.
# Usage: run_af2.sh <input.fasta> <output_dir> [gpu]
# AF2 code: /data/alphafold_code/alphafold (injected by run_alphafold_wrapper.py)
# Weights/DBs: /data/alphafold_db   |   env: /data/rfdiffusion/.venv-af2
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/common.sh"

FASTA="$1"; OUTDIR="$2"; GPU="${3:-0}"
mkdir -p "$OUTDIR"

export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="$PROJECT/logs/.mplconfig"; mkdir -p "$MPLCONFIGDIR"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export CUDA_VISIBLE_DEVICES="$GPU"

# Idempotency: skip if this FASTA already produced a ranked structure.
name="$(grep -m1 '^>' "$FASTA" >/dev/null 2>&1 && basename "$FASTA" .fasta || basename "$FASTA" .fasta)"
if [ -f "$OUTDIR/$name/ranking_debug.json" ]; then
  echo "  → AF2 already done for $name (skipping)"; exit 0
fi

"$AF2_PY" -u "$AF2_WRAPPER" \
  --fasta_paths="$FASTA" \
  --output_dir="$OUTDIR" \
  --data_dir="$AF2_DB" \
  --max_template_date="$AF2_MAX_TEMPLATE_DATE" \
  --db_preset=full_dbs \
  --model_preset=multimer \
  --num_multimer_predictions_per_model=1 \
  --jackhmmer_binary_path=/usr/bin/jackhmmer \
  --hhblits_binary_path=/usr/bin/hhblits \
  --hhsearch_binary_path=/usr/bin/hhsearch \
  --hmmsearch_binary_path=/usr/bin/hmmsearch \
  --hmmbuild_binary_path=/usr/bin/hmmbuild \
  --kalign_binary_path=/usr/bin/kalign \
  --uniref90_database_path="$AF2_DB/uniref90/uniref90.fasta" \
  --mgnify_database_path="$AF2_DB/mgnify/mgy_clusters_2022_05.fa" \
  --bfd_database_path="$AF2_DB/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt" \
  --uniref30_database_path="$AF2_DB/uniref30/UniRef30_2021_03" \
  --uniprot_database_path="$AF2_DB/uniprot/uniprot.fasta" \
  --pdb_seqres_database_path="$AF2_DB/pdb_seqres/pdb_seqres.txt" \
  --template_mmcif_dir="$AF2_DB/pdb_mmcif/mmcif_files" \
  --obsolete_pdbs_path="$AF2_DB/pdb_mmcif/obsolete.dat" \
  --use_precomputed_msas=True \
  --use_gpu_relax=False \
  --models_to_relax=none
