#!/usr/bin/env bash
# Reproducibly set up the af2_initial_guess validator (dl_binder_design) on this machine.
# external/ is git-ignored (vendored 282-file repo + a multi-GB conda env), so this script
# records exactly how the working setup was built. Run once.
#
# Result: a dedicated conda env `af2ig` that can run af2_initial_guess (single-sequence,
# templated, NO MSA — so no systemd-oomd kills, ~50 s/design). BindCraft env is untouched.
set -euo pipefail

CONDA=/home/a-mxn833/mambaforge/bin/conda
ENVS=/home/a-mxn833/mambaforge/envs
EXT=/data/binder_software/pre-binder/external
IG="$ENVS/af2ig/bin/python"

mkdir -p "$EXT" && cd "$EXT"

# 1. Vendor dl_binder_design (af2_initial_guess + bundled AlphaFold + silent_tools)
[ -d dl_binder_design ] || git clone --depth 1 https://github.com/nrbennet/dl_binder_design.git

# 2. Dedicated env: clone BindCraft (gives pyrosetta + jax0.4.28 + haiku + dm-tree + ml_collections)
$CONDA env list | grep -qw af2ig || $CONDA create --clone BindCraft -n af2ig -y

# 3. The clone pulls mismatched jax-cuda12 plugins (break jaxlib). Remove them and restore the
#    known-working jax/jaxlib from BindCraft (monolithic jaxlib==0.4.28+cuda12.cudnn89).
"$ENVS/af2ig/bin/pip" uninstall -y jax-cuda12-plugin jax-cuda12-pjrt 2>/dev/null || true
BCSP="$ENVS/BindCraft/lib/python3.10/site-packages"
IGSP="$ENVS/af2ig/lib/python3.10/site-packages"
rm -rf "$IGSP"/jaxlib "$IGSP"/jaxlib-*.dist-info "$IGSP"/jax "$IGSP"/jax-*.dist-info
cp -a "$BCSP"/jaxlib "$BCSP"/jaxlib-*.dist-info "$BCSP"/jax "$BCSP"/jax-*.dist-info "$IGSP"/

# 4. AF2's tf-based feature code needs tensorflow (CPU-only avoids a CUDA clash with jaxlib)
"$ENVS/af2ig/bin/pip" install -q tensorflow-cpu

# 5. Patch vendored AlphaFold for modern Biopython (SCOPData moved to PDBData)
MMCIF="$EXT/dl_binder_design/af2_initial_guess/alphafold/data/mmcif_parsing.py"
grep -q "PDBData as SCOPData" "$MMCIF" || sed -i \
  '/^from Bio.Data import SCOPData/c\try:\n    from Bio.Data import SCOPData\nexcept ImportError:\n    from Bio.Data import PDBData as SCOPData' "$MMCIF"

# 6. AF2 model_1_ptm weights (initial-guess uses only this model)
mkdir -p "$EXT/dl_binder_design/af2_initial_guess/model_weights/params"
ln -sf /data/alphafold_db/params/params_model_1_ptm.npz \
  "$EXT/dl_binder_design/af2_initial_guess/model_weights/params/params_model_1_ptm.npz"

# 7. Verify
"$IG" - <<'PY'
import sys; sys.path.insert(0, "/data/binder_software/pre-binder/external/dl_binder_design/af2_initial_guess")
import jax, tensorflow, pyrosetta
from alphafold.model import model
print("af2ig OK:", jax.__version__, jax.devices())
PY
echo "✓ af2ig ready."
