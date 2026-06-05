#!/usr/bin/env bash
# Phase 4 — final AF2 re-validation + specificity acid test + ranking.
# For each Phase-3 MPNN sequence: fold trimer + hexamer with AF2-multimer, apply the same
# 4 filters; for every hexamer-pass, re-fold against ONE dimer pair (chains A+E) and require
# interface pTM to drop by >= 0.15 (the operational definition of hexamer-specificity).
# Survivors are ranked by  iptm * (pLDDT/100) * (1/max(RMSD,0.5)) * (sum_subunit_SASA/1000).
#
# Output: outputs/04_final_ranked/rankNN_<design>_sK.pdb  +  outputs/04_final_metrics.csv
# Usage : [GPU=0] [NUM_SEQ=8] bash scripts/05_af2_revalidation.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"
"$AF2_PY" "$HERE/lib/af2_phase.py" final --gpu "${GPU:-0}" --num-seq "${NUM_SEQ:-8}"
