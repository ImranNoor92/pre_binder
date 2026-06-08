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

# progress
gate=$(ls outputs/02b_af2/*/ranking_debug.json 2>/dev/null | wc -l)
fin=$(ls outputs/04_af2/*_hex/ranking_debug.json 2>/dev/null | wc -l)
cache=$([ -d outputs/_msa_cache/target_A ] && echo "cached" || echo "building")
echo  "phase 2b    : $gate/10 backbones predicted   | target MSA: $cache"
echo  "phase 4     : $fin hex predictions done"
[ -f outputs/02b_af2_validated/validated.txt ] && echo "validated   : $(wc -l < outputs/02b_af2_validated/validated.txt) backbone(s)"
[ -d outputs/04_final_ranked ] && echo "final ranked: $(ls outputs/04_final_ranked/*.pdb 2>/dev/null | wc -l) design(s)"

# what AF2 is doing right now + log freshness
cur=$(ls -d outputs/02b_af2/design_*_hex outputs/04_af2/*_hex 2>/dev/null | tail -1)
[ -n "${cur:-}" ] && echo "current     : $(basename "$cur")  msas: $(ls "$cur"/msas/*/ 2>/dev/null | grep -c hits)"
LOG=$(ls -t logs/run_all_*.log 2>/dev/null | head -1)
if [ -n "${LOG:-}" ]; then
  age=$(( $(date +%s) - $(stat -c %Y "$LOG") ))
  echo "log         : $LOG  (last write ${age}s ago)"
  echo "last line   : $(tail -1 "$LOG" | sed 's/.*] //' | cut -c1-70)"
fi
