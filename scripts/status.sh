#!/usr/bin/env bash
# One-shot, unambiguous status of the pipeline run. Read-only.
#   bash scripts/status.sh           # print once
#   watch -n 30 bash scripts/status.sh   # live refresh every 30s
#
# Verdict logic (don't trust the log timestamp alone — MSA is silent for 10+ min):
#   ALIVE  = systemd unit active AND a worker (jackhmmer/hhblits/run_alphafold/mpnn) is burning CPU
#   IDLE?  = unit active but no busy worker right now (brief between-step gap, or wedged)
#   DEAD   = unit not active
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE/.."
UNIT=prebinder
bold() { printf '\033[1m%s\033[0m\n' "$1"; }

state=$(systemctl --user is-active "$UNIT" 2>/dev/null || echo "none")
since=$(systemctl --user show "$UNIT" -p ActiveEnterTimestamp --value 2>/dev/null)

# busiest worker process in the run (MSA = jackhmmer/hhblits; inference = python wrapper)
read -r wpid wcpu wcmd < <(ps -eo pid,pcpu,comm --no-headers 2>/dev/null \
  | grep -E "jackhmmer|hhblits|hhsearch|run_alphafo|protein_mpnn|python" \
  | sort -k2 -nr | head -1)
busy=0; [ -n "${wcpu:-}" ] && awk "BEGIN{exit !($wcpu>20)}" && busy=1

bold "── pipeline status ──"
echo  "unit        : $state   (since ${since:-?})"
if [ "$state" = "active" ] && [ "$busy" = 1 ]; then
  bold "VERDICT     : ✅ ALIVE — working"
  echo "worker      : $wcmd  pid=$wpid  cpu=${wcpu}%"
elif [ "$state" = "active" ]; then
  bold "VERDICT     : ⏳ active but no busy worker right now (between steps, or check again in 60s)"
  echo "worker      : ${wcmd:-none}  cpu=${wcpu:-0}%"
else
  bold "VERDICT     : ❌ NOT RUNNING (dead/finished) — see 'systemctl --user status $UNIT'"
fi

# progress (current pipeline: RFdiffusion -> ProteinMPNN -> af2_initial_guess)
bb=$(find outputs/01_rfdiffusion_pilot -maxdepth 1 -name 'design_*.pdb' 2>/dev/null | wc -l)
mpnn=$(ls outputs/03_mpnn_sequences/seqs/*.fa 2>/dev/null | wc -l)
igin=$(ls outputs/06_ig/inputs/*.pdb 2>/dev/null | wc -l)
igout=$(ls outputs/06_ig/out/*_af2pred.pdb 2>/dev/null | wc -l)
echo  "Phase 1 RFd : $bb backbones generated"
echo  "Phase 3 MPNN: $mpnn backbones sequenced"
echo  "Phase IG    : $igout/$igin initial-guess predictions"
if [ -f outputs/06_ig/ranked.csv ]; then
  np=$(awk -F, 'NR>1 && $NF=="True"{n++} END{print n+0}' outputs/06_ig/ranked.csv)
  echo "RESULT      : $np pass  →  outputs/06_ig/ranked.csv"
fi

LOG=$(ls -t logs/run_full_*.log logs/run_all_*.log 2>/dev/null | head -1)
if [ -n "${LOG:-}" ]; then
  age=$(( $(date +%s) - $(stat -c %Y "$LOG") ))
  echo "log         : $LOG  (last write ${age}s ago)"
  echo "last line   : $(tail -1 "$LOG" | sed 's/.*] //' | cut -c1-70)"
fi
