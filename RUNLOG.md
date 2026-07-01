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

### 2026-06-08 12:02 · Phase 2b→4 (run_all, tmux) · **FAILED (process reaped)**
- result: died ~12:27, 25 min into design_0 (at the uniprot pairing MSA). tmux session vanished.
- cause: the tmux server was spawned inside the agent's sandboxed shell and got reaped with that process tree. `enable-linger` guards against *logout*, not against the launching shell's teardown.
- salvage: design_0 target MSA (uniref90/mgnify/bfd/pdb) completed on disk → reused by next attempt.

### 2026-06-08 12:37 · Phase 2b→4 — full validation (run_all) · **RUNNING**
- inputs: 10 trimerized backbones
- params: GPU 0; AF2 multimer `full_dbs`; MPNN temp 0.1 seed 37, `NUM_SEQ=8`; filters pLDDT>70, iPTM>0.65, per-subunit ΔSASA>200 Å², RMSD<3 Å; acid test iPTM drop ≥0.15 (hexamer→A+E dimer)
- Δ from prev: launched via **`systemd-run --user --unit=prebinder`** (owned by the user systemd manager, cgroup `user@…/prebinder.service`) — survives both shell teardown and logout. Target-MSA reuse + chained gate→mpnn→final as before.
- result: **FAILED** — `systemd-oomd` OOM-killed it 2026-06-09 02:21 (266 procs, status 9/KILL) after designs 0–1, on design_2's MSA. Cause: full-BFD `hhblits` on the de-novo binder spikes memory pressure (box has only 2 GB swap). Also design_0 had failed all filters (pLDDT 34, RMSD 23 Å) — vanilla full-MSA AF2 doesn't reproduce de-novo designs.

### 2026-06-09 10:19 · Phase 3 → IG — validation via **af2_initial_guess** · **RUNNING**
- inputs: 10 RFdiffusion backbones (subunit vs one target chain)
- params: GPU 0; ProteinMPNN temp 0.1 seed 37, `NUM_SEQ=8`; validator = dl_binder_design `af2_initial_guess` (single-seq, templated, **no MSA**), model_1_ptm, 3 recycles; pass = `pae_interaction<10` AND `plddt_binder>80`
- Δ from prev: **retired vanilla full-MSA AF2** (OOM + slow + non-reproducing). New env `af2ig` (see `scripts/setup_af2ig.sh`). Order is now standard RFdiffusion→MPNN→AF2. ~2–3 s/design, no MSA so **no OOM**. Per-subunit interface validation; hexamer-specificity rests on the C3 construction (optional multichain re-check later).
- result: **COMPLETE, 0/80 pass** (80 designs in ~14 min, no crashes/OOM). Binders **fold** well (plddt_binder up to 83; binder_rmsd mostly 1–2 Å) but **don't bind**: best pae_interaction = 25.6, uniformly ~26 across all 80 (need < 10). Bottleneck = the *designs* (epitope / backbones / only 10), not the tooling. Best: `design_9_s4` (pae_i 25.6, plddt 82, rmsd 1.5).
- next: extend the target hotspot patch; scale backbone count. The C3-trimer hexamer-specificity test is still unbuilt (per-subunit validation only).
- artifacts: `outputs/06_ig/ranked.csv`, `outputs/06_ig/out/*_af2pred.pdb`, summary figure `outputs/06_ig/run_summary.png`

### 2026-06-10 09:31 · Phase 1→IG — 200 backbones, hydrophobic hotspots · **OK — 18/1600 pass** 🎯
- inputs: 200 RFdiffusion backbones (2-GPU split, 100/GPU), **hotspots `[A107,A116,A117,A120]`** (hydrophobic V/I/A/L, replaces polar 105–115)
- params: ProteinMPNN 8 seq/backbone (1600 designs); af2_initial_guess pass = pae_interaction<10 AND plddt_binder>80
- Δ from prev: **new hydrophobic-centred site** (paper: site needs ≥3 hydrophobics; old 105–115 had 1) + **20× scale** (10→200 backbones). Runtime ~7 h, durable (`systemd-run`), no OOM.
- result: **18/1600 PASS across 8 backbones** (design_16, 123, 69, 188, 47, 63, 172, 193). Best `design_16_s0`: pae_interaction **6.79**, plddt_binder 85.8, rmsd 0.6 Å. vs old polar 105–115 run = **0/80, best pae 25.6** → the hotspot change worked. (40/1600 reach pae<10; 18 also clear plddt>80.)
- artifacts: `outputs/06_ig/ranked.csv`, `reports/2026-06-10_phaseIG_200_hydrophobic.png`; old run archived at `outputs/_archive_2026-06-09_hs105-115/`
- note: RESULT line in the auto event-log said "0 pass" — a `\r\n` CSV counting bug (now fixed, `lineterminator="\n"`); true count is 18. Still per-subunit validation; C3-trimer specificity test next.

### 2026-06-11 · C3-trimer hexamer-specificity test · **ABANDONED (method limitation)**
- inputs: top-4 IG backbones (design_16/123/69/188); `scripts/lib/build_trimer_complex.py` (C3-replicate validated subunit about Z, place target chains)
- attempt: validate trimer-vs-hexamer (A,D,F) vs trimer-vs-dimer (A+E) interface with af2_initial_guess
- result: **af2_initial_guess cannot do the joint prediction** — it is the AF2 *monomer* net (single-seq) holding *one rigid target*; given 3 disconnected binder subunits + ≥3 non-contacting target protomers it rearranges everything (target self-RMSD 64 Å; pae_interaction 27.9 / binder_rmsd 24 / plddt 62 vs 6.79 / 0.6 / 85.8 per-subunit). Root cause: the engaged protomers A/D/F **don't contact each other** (16 Å apart) so nothing constrains their geometry. Architectural, not fixable.
- finding: the C3 specificity is geometric/by-construction (3 subunits engage on hexamer, 1 on isolated dimer) — the per-subunit interface was already af2ig-validated; a joint ML re-check is unnecessary and unreliable. → pivoted to a new model (next).

### 2026-06-14 · NEW CAMPAIGN · Phase 1 — single central binder vs full hexamer · **RUNNING**
- strategy: drop C3 subunit+trimerization. Design **one binder that caps the hexamer top and engages all six monomers** (see `context/new_hexamer_model.pdf`).
- target: full VP3 T=3 capsid hexamer `inputs/151lp3t3_hexamer.pdb` (6-MODEL → flattened to 6 unique chains A–F, 504 aa each → `inputs/151lp3t3_hexamer_6chain.pdb`). The full 3024-res target **OOMs RFdiffusion (>48 GB)**, so for backbone gen only it is reduced to everything within 50 Å of the central site (~1158 res, 193/chain, `inputs/151lp3t3_hexamer_R50_6chain.pdb`, multi-crop contig); the full hexamer is retained for validation.
- params: `ckpt=Complex_base`; **hotspots `[A111,B111,C111,D111,E111,F111]`** = exposed Asp111, the only residue that converges at the central axis-top (six copies in a ~22 Å cluster) — replaces the old buried/inner 105–120 that caused parallel single-monomer binding; binder **80–120 aa**; noise_ca/frame=0; **200 backbones**, 2-GPU split (100/GPU); ~14 min/design (~1 day); durable `systemd-run --user` units `rfd_central_gpu0/1`.
- pilot (design_0, 89 aa): caps centrally (binder centroid r_xy 2.7 Å), **exterior-facing**, contacts **4/6 monomers** at backbone level (poised over all six Asp111, 4.6–6.8 Å), target held rigid (motif RMSD 0.34 Å).
- plan (staged): 200 RFd → geometric monomers-contacted score (`scripts/lib/score_monomer_contacts.py`) → **1600 MPNN** (8/bb, binder chain G, target fixed) → **AF2 binder-only** foldability filter (`scripts/lib/thread_binder_only.py` + predict.py monomer; rank by plddt + self-consistency RMSD) → top 5 → **AF2 complex** (binder + 6 target chains merged to one chain) for binding + all-6 confirmation.
- downstream tools pre-tested on design_0 (MPNN ✓, binder-only AF2 ✓, scorer ✓).
- artifacts: `scripts/11_rfd_central.sh`, `scripts/lib/{score_monomer_contacts,thread_binder_only}.py`, `outputs/11_rfd_central/`; prior phase-IG candidates archived at `outputs/_archive_2026-06-10_phaseIG_hydrophobic/` (18 accepted).

---

## Event log (auto-appended by run_all.sh)
2026-06-08 12:02:03 EDT | RUN START  run_all (gate→mpnn→final)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=3640884
2026-06-08 12:36:58 EDT | RUN START  run_all (gate→mpnn→final)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=3650791
2026-06-08 13:46:08 EDT | phase2b gate  exit=1
2026-06-08 13:46:08 EDT | RUN ABORTED at gate
2026-06-08 18:46:45 EDT | RUN START  run_all (gate→mpnn→final)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=3744628
2026-06-09 10:17:22 EDT | RUN START  run_all (mpnn→initial-guess)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=3977158
2026-06-09 10:19:09 EDT | RUN START  run_all (mpnn→initial-guess)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=3977741
2026-06-09 10:25:10 EDT | phase3 mpnn        exit=0
2026-06-09 10:39:48 EDT | phaseIG initial-guess exit=0
2026-06-09 10:39:48 EDT | RUN COMPLETE  0 design(s) pass  →  outputs/06_ig/ranked.csv
2026-06-10 09:31:24 EDT | RUN START run_full  RFd=200 (2-GPU split)  hotspots=[A107,A116,A117,A120]  host=caspbioa01.as.acorn.miami.edu  pid=121526
2026-06-10 13:57:04 EDT | phase1 rfdiffusion exit=0/0  backbones=200
2026-06-10 13:57:04 EDT | RUN START  run_all (mpnn→initial-guess)  GPU=0  NUM_SEQ=8  host=caspbioa01.as.acorn.miami.edu  pid=717545
2026-06-10 14:19:36 EDT | phase3 mpnn        exit=0
2026-06-10 16:44:13 EDT | phaseIG initial-guess exit=0
2026-06-10 16:44:13 EDT | RUN COMPLETE  0 design(s) pass  →  outputs/06_ig/ranked.csv
2026-06-10 16:44:13 EDT | RUN FULL DONE
