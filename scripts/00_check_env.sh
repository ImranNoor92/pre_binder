#!/usr/bin/env bash
# Verify all tools, paths, and inputs exist before running anything.
# Safe to run any number of times. Run this FIRST before any other script.

set -u
fail=0

say() { printf '%-50s %s\n' "$1" "$2"; }
check() {
  if eval "$1" >/dev/null 2>&1; then say "$2" "✓ OK"; else say "$2" "✗ MISSING"; fail=$((fail+1)); fi
}

echo "=== RFdiffusion install ==="
check "[ -f /data/rfdiffusion/scripts/run_inference.py ]" "run_inference.py present"
check "[ -d /data/rfdiffusion/external/ProteinMPNN ]"    "ProteinMPNN present"
check "[ -f /data/rfdiffusion/external/ProteinMPNN/protein_mpnn_run.py ]" "  protein_mpnn_run.py present"

echo
echo "=== Conda / venv environments ==="
# RFdiffusion + ProteinMPNN run from .venv-rfd-gpu (torch/e3nn/dgl); AF2 runs from .venv-af2 (jax).
# We verify the binaries AND that the key modules actually import — a present binary with
# missing packages is the failure mode that previously slipped through this check.
RFD_VENV=/data/rfdiffusion/.venv-rfd-gpu
AF2_VENV=/data/rfdiffusion/.venv-af2
check "[ -x $RFD_VENV/bin/python ]" "RFdiffusion venv python ($RFD_VENV/bin/python)"
if [ -x "$RFD_VENV/bin/python" ]; then
  say "  Python version" "$($RFD_VENV/bin/python --version 2>&1)"
fi
check "$RFD_VENV/bin/python -c 'import torch,e3nn,dgl,hydra'" "  RFdiffusion deps import (torch,e3nn,dgl,hydra)"
check "$RFD_VENV/bin/python -c 'import torch; assert torch.cuda.is_available()'" "  torch CUDA available"
check "[ -x $AF2_VENV/bin/python ]" "AF2 venv python ($AF2_VENV/bin/python)"
check "$AF2_VENV/bin/python -c 'import jax,jaxlib'" "  AF2 deps import (jax,jaxlib)"

echo
echo "=== GPU availability ==="
check "command -v nvidia-smi" "nvidia-smi available"
if command -v nvidia-smi >/dev/null 2>&1; then
  echo "  Current GPU state:"
  nvidia-smi --query-gpu=index,name,utilization.gpu,memory.free --format=csv,noheader | sed 's/^/    /'
fi

echo
echo "=== Inputs ==="
TARGET_PDB=/data/binder_software/pre-binder/inputs/1lp3_hexamer_trimmed_fixed.pdb
check "[ -f $TARGET_PDB ]" "target PDB ($TARGET_PDB)"
if [ -f "$TARGET_PDB" ]; then
  CHAINS=$(awk '/^ATOM/{c[substr($0,22,1)]=1} END{for(k in c) printf "%s ",k}' "$TARGET_PDB")
  say "  Chains in target" "$CHAINS"
  RES=$(awk '/^ATOM/ && substr($0,22,1)=="A"{r[substr($0,23,4)+0]=1} END{n=0; min=99999; max=-99999; for(k in r){n++; kn=k+0; if(kn<min)min=kn; if(kn>max)max=kn} printf "%d residues, range %d-%d", n, min, max}' "$TARGET_PDB")
  say "  Chain A" "$RES"
fi

echo
echo "=== AlphaFold (for Phase 2b/4 validation) ==="
check "[ -d /data/alphafold_code ]" "alphafold_code present"
check "[ -d /data/alphafold_db ]"   "alphafold_db present"

echo
echo "=== ProteinMPNN model weights ==="
check "[ -d /data/rfdiffusion/external/ProteinMPNN/ca_model_weights ]" "MPNN CA model weights"

echo
echo "=== RFdiffusion model weights ==="
check "[ -d /data/rfdiffusion/models ]" "RFdiffusion models directory"
if [ -d /data/rfdiffusion/models ]; then
  N=$(ls /data/rfdiffusion/models/*.pt 2>/dev/null | wc -l)
  say "  .pt files" "$N"
fi

echo
echo "=== Summary ==="
if [ $fail -eq 0 ]; then
  echo "✓ All checks passed. Safe to proceed with 01_pilot_rfdiffusion.sh"
  exit 0
else
  echo "✗ $fail check(s) failed. Resolve before proceeding."
  exit 1
fi
