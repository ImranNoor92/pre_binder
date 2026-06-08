#!/usr/bin/env bash
# Run the validation half end-to-end: Phase 2b (gate) -> Phase 3 (MPNN) -> Phase 4 (AF2+acid+rank).
# Designed to run detached inside tmux so it survives logout (see README "Running long jobs").
# Appends a structured entry to RUNLOG.md at start and finish.
#
# Usage:  GPU=0 bash scripts/run_all.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"
GPU="${GPU:-0}"; NUM_SEQ="${NUM_SEQ:-8}"
ts() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { printf '%s | %s\n' "$(ts)" "$1" | tee -a "$PROJECT/RUNLOG.md"; }

log "RUN START  run_all (gate→mpnn→final)  GPU=$GPU  NUM_SEQ=$NUM_SEQ  host=$(hostname)  pid=$$"

GPU="$GPU" bash "$HERE/03_af2_validation.sh";  rc=$?; log "phase2b gate  exit=$rc"
[ $rc -eq 0 ] || { log "RUN ABORTED at gate"; exit $rc; }

GPU="$GPU" bash "$HERE/04_proteinmpnn.sh";     rc=$?; log "phase3 mpnn   exit=$rc"
[ $rc -eq 0 ] || { log "RUN ABORTED at mpnn"; exit $rc; }

GPU="$GPU" NUM_SEQ="$NUM_SEQ" bash "$HERE/05_af2_revalidation.sh"; rc=$?; log "phase4 final  exit=$rc"

npass=$(awk -F, 'NR>1 && $0 ~ /True/{n++} END{print n+0}' "$PROJECT/outputs/04_final_metrics.csv" 2>/dev/null || echo 0)
log "RUN COMPLETE  ranked dir=outputs/04_final_ranked  (see outputs/04_final_metrics.csv)"
