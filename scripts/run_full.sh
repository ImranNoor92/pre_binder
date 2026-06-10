#!/usr/bin/env bash
# Full pipeline from scratch: RFdiffusion (split across both GPUs) -> ProteinMPNN -> af2_initial_guess.
# Launch durably via systemd-run --user (see README "Running long jobs").
#
# Usage:  TOTAL=200 bash scripts/run_full.sh      (optionally: HOTSPOTS='[A107,A116,A117,A120]')
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"
TOTAL="${TOTAL:-200}"; HALF=$(( TOTAL / 2 )); REST=$(( TOTAL - HALF ))
export HOTSPOTS="${HOTSPOTS:-[A107,A116,A117,A120]}"
ts() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { printf '%s | %s\n' "$(ts)" "$1" | tee -a "$PROJECT/RUNLOG.md"; }

log "RUN START run_full  RFd=$TOTAL (2-GPU split)  hotspots=$HOTSPOTS  host=$(hostname)  pid=$$"

# Phase 1 — split: GPU0 designs [0..HALF), GPU1 designs [HALF..TOTAL)
GPU=0 NUM_DESIGNS="$HALF" DESIGN_STARTNUM=0      bash "$HERE/01_pilot_rfdiffusion.sh" & p0=$!
GPU=1 NUM_DESIGNS="$REST" DESIGN_STARTNUM="$HALF" bash "$HERE/01_pilot_rfdiffusion.sh" & p1=$!
wait $p0; r0=$?; wait $p1; r1=$?
ndes=$(find "$PROJECT/outputs/01_rfdiffusion_pilot" -maxdepth 1 -name 'design_*.pdb' | wc -l)
log "phase1 rfdiffusion exit=$r0/$r1  backbones=$ndes"
{ [ "$r0" -eq 0 ] && [ "$r1" -eq 0 ]; } || { log "RUN ABORTED at rfdiffusion"; exit 1; }

# Phase 3 (MPNN, all backbones) -> Phase IG (thread + initial-guess + rank)
GPU=0 NUM_SEQ="${NUM_SEQ:-8}" bash "$HERE/run_all.sh"
log "RUN FULL DONE"
