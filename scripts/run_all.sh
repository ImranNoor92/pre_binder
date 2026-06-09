#!/usr/bin/env bash
# Run the validation half end-to-end with the af2_initial_guess validator:
#   Phase 3 (ProteinMPNN, all backbones) -> Phase IG (thread + initial-guess + rank).
# This is the standard RFdiffusion->MPNN->AF2 order. The old vanilla full-MSA AF2 gate
# (scripts/03, 05) is retired — it was OOM-killed, slow, and didn't reproduce de-novo designs.
#
# Launch durably (survives logout) via systemd-run --user — see README "Running long jobs".
# Appends a structured entry to RUNLOG.md at start and after each phase.
#
# Usage:  GPU=0 bash scripts/run_all.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"
GPU="${GPU:-0}"; NUM_SEQ="${NUM_SEQ:-8}"
ts() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { printf '%s | %s\n' "$(ts)" "$1" | tee -a "$PROJECT/RUNLOG.md"; }

log "RUN START  run_all (mpnn→initial-guess)  GPU=$GPU  NUM_SEQ=$NUM_SEQ  host=$(hostname)  pid=$$"

GPU="$GPU" NUM_SEQ_PER_BACKBONE="$NUM_SEQ" bash "$HERE/04_proteinmpnn.sh"; rc=$?; log "phase3 mpnn        exit=$rc"
[ $rc -eq 0 ] || { log "RUN ABORTED at mpnn"; exit $rc; }

GPU="$GPU" bash "$HERE/06_ig_validate.sh"; rc=$?; log "phaseIG initial-guess exit=$rc"
[ $rc -eq 0 ] || { log "RUN ABORTED at initial-guess"; exit $rc; }

npass=$(awk -F, 'NR>1 && $NF=="True"{n++} END{print n+0}' "$PROJECT/outputs/06_ig/ranked.csv" 2>/dev/null || echo 0)
log "RUN COMPLETE  $npass design(s) pass  →  outputs/06_ig/ranked.csv"
