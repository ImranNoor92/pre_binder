#!/usr/bin/env bash
# Phase 2b — backbone gate. For each RFdiffusion backbone: design 1 quick ProteinMPNN
# sequence for the subunit, fold the C3-fused trimer + hexamer with AlphaFold-multimer,
# and apply the 4 filters (binder pLDDT>70, interface pTM>0.65, per-subunit SASA>200 A^2,
# RMSD-to-design<3 A). Passing backbones -> outputs/02b_af2_validated/.
#
# AF2 needs a sequence to fold, so this stage produces one (the report's "validate the
# trimerized backbone" step); Phase 3 then designs the fuller sequence set on survivors.
#
# Wall time: dominated by AF2 (~minutes-hours/design; target MSA computed once, reused).
# Usage: [GPU=0] bash scripts/03_af2_validation.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"
"$AF2_PY" "$HERE/lib/af2_phase.py" gate --gpu "${GPU:-0}"
