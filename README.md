# pre_binder — C3-symmetric trimer binder for the AAV T=3 capsid hexamer

A staged pipeline (RFdiffusion → C3 replication → AlphaFold-multimer → ProteinMPNN → AF2 re-validation) for designing a binder that is **hexamer-specific by construction**: it engages the six-chain "rosette" that exists only in the T=3 icosahedral capsid, and *cannot* bind an isolated VP3 dimer (which would also occur in T=1 particles).

> **Background.** This approach was adopted after BindCraft single-chain hallucination (Trials 5 and 6) failed to produce a hexamer-specific binder for this target — a single chain cannot bridge the cross-dimer-pair patches of a C3-symmetric target without colliding with the intervening chains. The full retrospective and design rationale are in [`docs/`](docs/) and the source report (`AAV_T3_binder_design_report`).

---

## The idea in one paragraph

The hexamer has exact three-fold rotational symmetry about an axis through `(0, 0, 180.8) Å`. Instead of one binder that must reach three patches at once, we design **one small subunit** against a single chain (no bridging requirement), then **replicate it by 120°/240° rotation** about that verified C3 axis and fuse the three copies into one polypeptide with flexible linkers. On the full hexamer all three subunits engage simultaneously (strong, multivalent binding); on an isolated dimer at most one subunit finds a partner, so affinity collapses — **specificity by construction**. The symmetry step is a deterministic geometric rotation, not machine learning, which sidesteps RFdiffusion's known issues combining symmetric-PPI mode with hotspots.

---

## Pipeline

| Phase | Script | What it does | Engine | Status |
|------:|--------|--------------|--------|--------|
| 1  | `scripts/01_pilot_rfdiffusion.sh` | Design ~10 single-subunit binder backbones against chain A (B–F as fixed steric context, hotspots `A105,107,109,111,114,115`) | RFdiffusion (`.venv-rfd-gpu`) | ✅ **runnable** |
| 2a | `scripts/02_trimerize_replicate.py` | Replicate each subunit by C3 (0°/120°/240°) about the hexamer axis; fuse with `(GGGGS)`-type linkers | pure Python | ✅ **runnable** |
| 2b | `scripts/03_af2_validation.sh` | Backbone gate: 1 quick MPNN seq → AF2-multimer (hexamer+trimer) → 4 filters (binder pLDDT > 70, interface pTM > 0.65, per-subunit interface SASA > 200 Å², RMSD vs design < 3 Å) | ProteinMPNN + AlphaFold | ✅ **wired** |
| 3  | `scripts/04_proteinmpnn.sh` | Sequence design of the subunit (8 seqs); one subunit designed and reused in all 3 copies = tied positions for a replication-built homotrimer; target chains fixed | ProteinMPNN (`.venv-rfd-gpu`) | ✅ **wired** |
| 4  | `scripts/05_af2_revalidation.sh` | AF2 each MPNN sequence; re-apply 4 filters; rank survivors by combined score | AlphaFold (`.venv-af2`) | ✅ **wired** |
| Acid test | (in `05`) | Re-predict each survivor against one dimer pair (A+E only); require interface pTM to drop ≥ 0.15 — the operational test of hexamer-specificity | AlphaFold | ✅ **wired** |

**Status:** all phases wired and runnable. Shared helpers live in [`scripts/lib/`](scripts/lib/); see "Where … are called from" below. Phases 1 and 2a have been run (10 backbones, 10 trimers); the AF2 phases are compute-heavy (see the efficiency note).

---

## Environment

This repo does **not** ship its own environment — it reuses the existing RFdiffusion install on this machine. Three venvs under `/data/rfdiffusion/`:

| venv | Provides | Used by |
|------|----------|---------|
| `/data/rfdiffusion/.venv-rfd-gpu` | RFdiffusion + ProteinMPNN (torch 1.12+cu116, e3nn, dgl, se3) | Phases 1, 3 |
| `/data/rfdiffusion/.venv-af2`     | AlphaFold (jax 0.4.26); AF2 code at `/data/alphafold_code`, DB at `/data/alphafold_db` | Phases 2b, 4 |

Model weights live at `/data/rfdiffusion/models/` (the pip-installed `rfdiffusion` package's built-in relative lookup is empty, so Phase 1 passes the checkpoint explicitly via `inference.ckpt_override_path`).

Verify everything before running:

```bash
bash scripts/00_check_env.sh   # checks binaries, real module imports, CUDA, weights, inputs
```

---

## Where RFdiffusion / ProteinMPNN / AF2 are called from

Nothing is vendored into this repo — every engine is invoked from the shared
`/data/rfdiffusion` and `/data/alphafold_*` installs. All call-sites are defined once in
[`scripts/lib/common.sh`](scripts/lib/common.sh); the phase scripts source that file.

| Engine | Env (interpreter) | Code / entry point | Weights / DB | Invoked by |
|--------|-------------------|--------------------|--------------|------------|
| **RFdiffusion** | `/data/rfdiffusion/.venv-rfd-gpu/bin/python` | `/data/rfdiffusion/scripts/run_inference.py` | `/data/rfdiffusion/models/Complex_base_ckpt.pt` (passed via `inference.ckpt_override_path`) | `scripts/01_pilot_rfdiffusion.sh` |
| **ProteinMPNN** | `/data/rfdiffusion/.venv-rfd-gpu/bin/python` | `/data/rfdiffusion/external/ProteinMPNN/protein_mpnn_run.py` (+ `helper_scripts/`) | bundled `vanilla_model_weights/` | `scripts/lib/mpnn_subunit.sh` (used by Phases 2b & 3) |
| **AlphaFold-multimer** | `/data/rfdiffusion/.venv-af2/bin/python` | `/data/alphafold_code/alphafold/run_alphafold.py` via `…/trials/trial_1/run_alphafold_wrapper.py` (injects it on `PYTHONPATH`) | `/data/alphafold_db` (`model_preset=multimer`, `full_dbs`) | `scripts/lib/run_af2.sh` (used by Phases 2b & 4) |

Why two Python envs: `.venv-rfd-gpu` carries torch/e3nn/dgl (RFdiffusion + ProteinMPNN);
`.venv-af2` carries jax/haiku/openmm (AlphaFold). They are mutually exclusive, so each phase
shells out to the helper that runs under the correct interpreter — orchestration
([`scripts/lib/af2_phase.py`](scripts/lib/af2_phase.py)) runs under `.venv-af2` and calls the
MPNN helper as a subprocess.

> **Target-MSA reuse (implemented):** the six target chains are identical, so their (expensive,
> `full_dbs`) MSA is computed once on the first design, cached at `outputs/_msa_cache/target_A/`,
> and reused for every later design (`af2_phase.py:seed_target_msa`/`save_target_msa`). Only the
> de-novo binder chain's MSA is recomputed per design. If the cache is stale, AF2 just recomputes
> (`--use_precomputed_msas`) — it never fails. `reduced_dbs` is *not* available here (`small_bfd`
> is not downloaded). The remaining per-design cost is the binder MSA; installing
> `af2_initial_guess` (single-sequence + templating) would remove that too.

---

## Running long jobs (so they survive logout)

The AF2 phases run for hours. A plain background job — even under `tmux` — **gets killed when
the launching shell's process tree is torn down** (both failures are recorded in `RUNLOG.md`,
2026-06-05 logout and 2026-06-08 tmux-reaped). The durable fix is to hand the job to the
**systemd user manager** so it's owned by `user@.service`, not by your shell:

```bash
cd /data/binder_software/pre-binder
loginctl enable-linger "$USER"          # once: lets user processes outlive logout

systemd-run --user --unit=prebinder --collect \
  --working-directory="$PWD" --setenv=GPU=0 --setenv=NUM_SEQ=8 \
  bash -c 'exec bash scripts/run_all.sh >> logs/run_all_$(date +%Y%m%d_%H%M%S).log 2>&1'
```

`run_all.sh` chains Phase 2b → 3 → 4 and appends a timestamped block to **`RUNLOG.md`** at start
and after each phase. To run a single phase instead, call its `scripts/0X_*.sh` directly.

### Checking the run (don't trust the log timestamp)

During MSA the log can sit silent for 10–15 min and the GPU shows 0% — that looks dead but
isn't. Use `scripts/status.sh`, which checks the authoritative signals (systemd unit state +
whether a `jackhmmer`/`hhblits`/`run_alphafold` worker is burning CPU) and prints a flat verdict:

```bash
bash scripts/status.sh                 # one-shot: ✅ ALIVE / ⏳ idle / ❌ NOT RUNNING + progress
watch -n 30 bash scripts/status.sh     # live, refreshing

# bare-bones, no script:
systemctl --user is-active prebinder   # "active" = alive; anything else = dead/finished
systemctl --user status  prebinder     # full detail + recent output
systemctl --user stop    prebinder     # kill it
```

---

## Quick start

```bash
cd /data/binder_software/pre-binder

bash   scripts/00_check_env.sh                                 # ~30 s — must be all ✓
bash   scripts/01_pilot_rfdiffusion.sh  2>&1 | tee logs/01.log # ~2–4 h on 1 GPU → 10 backbones
python scripts/02_trimerize_replicate.py                       # ~1 min  → trimerized PDBs
GPU=0  bash scripts/03_af2_validation.sh   2>&1 | tee logs/03.log # Phase 2b gate (MPNN+AF2 filter)
GPU=0  bash scripts/04_proteinmpnn.sh      2>&1 | tee logs/04.log # Phase 3   (8 seqs/survivor)
GPU=0  bash scripts/05_af2_revalidation.sh 2>&1 | tee logs/05.log # Phase 4   (AF2 + acid test + rank)
```

Overrides: Phase 1 takes `GPU` (default 1), `NUM_DESIGNS`, `BINDER_MIN`/`BINDER_MAX` (default 60–90); the AF2/MPNN phases take `GPU` (default 0) and `NUM_SEQ_PER_BACKBONE`/`NUM_SEQ` (default 8). Every script is idempotent — it skips work whose output already exists, so a re-run resumes.

---

## Target geometry

Source structure: `inputs/1lp3_hexamer_trimmed_fixed.pdb` (copied from the AAV T=3 analysis).

| | |
|---|---|
| Topology | 6 chains (A–F), residues 70–150 each, single model |
| Symmetry | C3, axis through ≈ `(0, 0, 180.8) Å` |
| Dimer pairs | A↔E, B↔D, C↔F |
| Hotspots | residues 105–115; the six solvent-exposed ones (`105,107,109,111,114,115`) are used |

---

## Repository layout

```
pre_binder/
├── README.md
├── docs/
│   ├── 01_workflow.md             # phase-by-phase technical detail
│   ├── 02_rfdiffusion_contigs.md  # contig syntax for this target
│   └── 03_validation_filters.md   # AF2 filter rationale
├── RUNLOG.md                      # append-only, structured per-run log
├── scripts/
│   ├── 00_check_env.sh            # environment + inputs verification
│   ├── 01_pilot_rfdiffusion.sh    # Phase 1   RFdiffusion
│   ├── 02_trimerize_replicate.py  # Phase 2a  C3 trimerize
│   ├── 03_af2_validation.sh       # Phase 2b  gate (MPNN+AF2 filter)
│   ├── 04_proteinmpnn.sh          # Phase 3   ProteinMPNN
│   ├── 05_af2_revalidation.sh     # Phase 4   AF2 + acid test + rank
│   ├── run_all.sh                 # chains 2b→3→4 (launch via systemd-run)
│   ├── status.sh                  # ALIVE/DEAD + progress check
│   └── lib/                       # common.sh, seqtools.py, af2_metrics.py,
│                                  #   run_af2.sh, mpnn_subunit.sh, af2_phase.py
├── inputs/
│   ├── 1lp3_hexamer_trimmed_fixed.pdb
│   ├── hotspot_residues.txt
│   └── target_pdb_link.md
├── outputs/                       # generated; git-ignored
└── logs/                          # generated; git-ignored
```

`outputs/` and `logs/` are git-ignored (regenerable, potentially large AF2 artifacts). Everything needed to reproduce a run is committed.
