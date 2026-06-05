# Target PDB

**Source:** `/data/binder_software/pre-binder/inputs/1lp3_hexamer_trimmed_fixed.pdb`

Same hexamer used by trial 5 and trial 6 (single-MODEL, 6 chains A-F, residues 70-150). We point to it rather than copying because:
- Avoids a stale copy if the target ever needs revising
- Trial 6's PDB has been validated end-to-end through BindCraft already

**To use:** Scripts reference this path via the `$TARGET` variable. No copying needed.

**To replace the target** (e.g., if you discover a different epitope range works better):
1. Update `$TARGET` in `scripts/01_pilot_rfdiffusion.sh` (and others)
2. Re-verify residue ranges in `scripts/02_trimerize_replicate.py` if you change the chain residue ranges
3. Re-check the C3 axis location — if the new PDB isn't centered at (0, 0, ~180.8), update the rotation reference in `02_trimerize_replicate.py`
