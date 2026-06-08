# RUNLOG — pre_binder pipeline

Append-only, chronological (oldest first). One block per run. Keep entries terse.
**Schema:** `date` · `phase` · `status` · `inputs` · `params` · `Δ from prev` · `result` · `artifacts`.
Status ∈ {OK, FAILED, RUNNING}. The `run_all.sh` launcher auto-appends timestamped lines to **Event log** at the bottom.

---

### 2026-06-05 · Phase 1 — RFdiffusion · **OK**
- inputs: `inputs/1lp3_hexamer_trimmed_fixed.pdb` (hexamer A–F)
- params: `NUM_DESIGNS=10`, binder 60–90 aa, hotspots A105,107,109,111,114,115, `ckpt=Complex_base`, GPU 1
- Δ: first run; fixed `ls`-glob idempotency bug + `inference.ckpt_override_path`
- result: 10/10 backbones
- artifacts: `outputs/01_rfdiffusion_pilot/design_{0..9}.pdb`

### 2026-06-05 · Phase 2a — C3 trimerize · **OK**
- inputs: 10 Phase-1 backbones
- params: 120°/240° about axis (0,0,180.8); linker GGGGSGGS; binder chain G
- Δ: none
- result: 10/10 trimers (subunits 62–90 aa; trimers 202–286 aa)
- artifacts: `outputs/02a_trimerized/design_*_trimer.pdb`

### 2026-06-05 · Phase 2b — AF2 gate (attempt 1) · **FAILED (environment)**
- inputs: 10 trimerized backbones
- params: AF2 multimer `full_dbs`, GPU 0
- Δ: first AF2 wiring
- result: **0/10.** Process died on session logout (`Linger=no`); log froze mid-MSA on design_0 after ~7 min. No compute lost to error — just not persisted.
- fix applied → see attempt 2

### 2026-06-08 · Phase 2b→4 — full validation (run_all) · **RUNNING**
- inputs: 10 trimerized backbones
- params: GPU 0; AF2 multimer `full_dbs`; MPNN temp 0.1 seed 37, `NUM_SEQ=8`; filters pLDDT>70, iPTM>0.65, per-subunit ΔSASA>200 Å², RMSD<3 Å; acid test iPTM drop ≥0.15 (hexamer→A+E dimer)
- Δ from prev: **(1)** `loginctl enable-linger` + run inside **tmux** so it survives logout; **(2)** **target-MSA reuse** — the 6 identical target chains' MSA is computed once and reused across all designs (only the binder MSA recomputes); **(3)** chained gate→mpnn→final via `scripts/run_all.sh`
- result: _(running)_
- artifacts: `logs/run_all_*.log`, `outputs/02b_*`, `outputs/03_mpnn_sequences/`, `outputs/04_final_ranked/`, `outputs/04_final_metrics.csv`

---

## Event log (auto-appended by run_all.sh)
2026-06-08 12:02:03 EDT | RUN START  run_all (gate→mpnn→final)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=3640884
