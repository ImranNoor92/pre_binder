# pre_binder — C3-symmetric trimer binder for the AAV T=3 capsid hexamer

A staged pipeline (RFdiffusion → C3 replication → AlphaFold-multimer → ProteinMPNN → AF2 re-validation) for designing a binder that is **hexamer-specific by construction**: it engages the six-chain "rosette" that exists only in the T=3 icosahedral capsid, and *cannot* bind an isolated VP3 dimer (which would also occur in T=1 particles).

> **Background.** This approach was adopted after BindCraft single-chain hallucination (Trials 5 and 6) failed to produce a hexamer-specific binder for this target — a single chain cannot bridge the cross-dimer-pair patches of a C3-symmetric target without colliding with the intervening chains. The full retrospective and design rationale are in [`docs/`](docs/) and the source report (`AAV_T3_binder_design_report`).

---

## The idea in one paragraph

The hexamer has exact three-fold rotational symmetry about an axis through `(0, 0, 180.8) Å`. Instead of one binder that must reach three patches at once, we design **one small subunit** against a single chain (no bridging requirement), then **replicate it by 120°/240° rotation** about that verified C3 axis and fuse the three copies into one polypeptide with flexible linkers. On the full hexamer all three subunits engage simultaneously (strong, multivalent binding); on an isolated dimer at most one subunit finds a partner, so affinity collapses — **specificity by construction**. The symmetry step is a deterministic geometric rotation, not machine learning, which sidesteps RFdiffusion's known issues combining symmetric-PPI mode with hotspots.

---

## Pipeline

The active pipeline is the standard de-novo binder flow **RFdiffusion → ProteinMPNN → AF2 (initial guess)**:

| Phase | Script | What it does | Engine | Status |
|------:|--------|--------------|--------|--------|
| 1  | `scripts/01_pilot_rfdiffusion.sh` | Design single-subunit binder backbones against chain A (B–F as fixed steric context, hotspots `A105,107,109,111,114,115`) | RFdiffusion (`.venv-rfd-gpu`) | ✅ run (10 backbones) |
| 2a | `scripts/02_trimerize_replicate.py` | Replicate each subunit by C3 (0/120/240°) about the hexamer axis; fuse with `(GGGGS)` linkers | pure Python | ✅ run (10 trimers) |
| 3  | `scripts/04_proteinmpnn.sh` | Sequence-design the subunit (8 seqs); one subunit designed and reused in all 3 copies (= tied positions for a replication-built homotrimer); target fixed | ProteinMPNN (`.venv-rfd-gpu`) | ✅ **active** |
| IG | `scripts/06_ig_validate.sh` | Thread each MPNN seq onto the backbone, pair with one target chain, validate with **af2_initial_guess** (single-seq, templated, no MSA); filter `pae_interaction<10` & `plddt_binder>80`; rank | af2_initial_guess (`af2ig` env) | ✅ **active** |

Run end-to-end with `scripts/run_all.sh` (Phase 3 → IG). Helpers live in [`scripts/lib/`](scripts/lib/).

> **Retired:** `scripts/03_af2_validation.sh` / `05_af2_revalidation.sh` (full-MSA AF2-multimer gate + dimer-pTM acid test). Kept in-tree for reference but **not used** — see the Environment section for why (OOM + slow + non-reproducing). Hexamer-specificity now rests on the C3 construction, with an optional multi-chain finalist re-check available later.

---

## Environment

This repo does **not** ship its own environment — it reuses the existing RFdiffusion install on this machine. Three venvs under `/data/rfdiffusion/`:

| env | Provides | Used by |
|------|----------|---------|
| `/data/rfdiffusion/.venv-rfd-gpu` | RFdiffusion + ProteinMPNN (torch 1.12+cu116, e3nn, dgl, se3) | Phases 1, 3 |
| `af2ig` conda env (`~/mambaforge/envs/af2ig`) | af2_initial_guess validator: jax 0.4.28+cuda, pyrosetta, tensorflow-cpu, bundled AlphaFold | Phase IG |
| `/data/rfdiffusion/.venv-af2` | (legacy) full-MSA AlphaFold-multimer — retired as the validator, see below | — |

Model weights: RFdiffusion at `/data/rfdiffusion/models/` (Phase 1 passes `inference.ckpt_override_path`); AF2 `model_1_ptm` at `/data/alphafold_db/params/`. The `af2ig` env is built reproducibly by **`scripts/setup_af2ig.sh`** (clones the BindCraft conda env, repairs jax, adds tensorflow-cpu, patches the vendored AlphaFold, links weights). `external/` (the vendored `dl_binder_design` + that env) is git-ignored.

> **Why af2_initial_guess, not full-MSA AF2 (learned the hard way):** vanilla AlphaFold-multimer with full MSA was (1) **OOM-killed** by `systemd-oomd` (full-BFD `hhblits` on the de-novo binder, only 2 GB swap), (2) ~1 h/design, and (3) didn't reproduce de-novo designs (RMSD ~23 Å). `af2_initial_guess` is single-sequence + templated on the design: **no MSA → no OOM**, ~2–3 s/design, and it validates the actual designed interface. It requires binder+target as exactly **2 chains** (binder = chain A), so we validate **each subunit against one target chain**; hexamer-specificity rests on the C3 construction (an optional multi-chain re-check can be run on finalists later). Full RUNLOG: `RUNLOG.md` 2026-06-08/09.

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
| **ProteinMPNN** | `/data/rfdiffusion/.venv-rfd-gpu/bin/python` | `/data/rfdiffusion/external/ProteinMPNN/protein_mpnn_run.py` (+ `helper_scripts/`) | bundled `vanilla_model_weights/` | `scripts/lib/mpnn_subunit.sh` (Phase 3) |
| **af2_initial_guess** | `~/mambaforge/envs/af2ig/bin/python` | `external/dl_binder_design/af2_initial_guess/predict.py` (bundled AlphaFold, single-seq + initial guess, no MSA) | `model_1_ptm` (→ `/data/alphafold_db/params`) | `scripts/06_ig_validate.sh` (Phase IG) |
| **AlphaFold-multimer** *(legacy)* | `/data/rfdiffusion/.venv-af2/bin/python` | `/data/alphafold_code/alphafold/run_alphafold.py` via the `trial_1` wrapper | `/data/alphafold_db` (`full_dbs`) | `scripts/lib/run_af2.sh` — retired as validator (kept for an optional multi-chain finalist check) |

Each tool runs under its own interpreter (mutually exclusive deps), so the phase scripts shell
out to the helper that runs under the correct env.

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
bash   scripts/setup_af2ig.sh                                  # once — builds the af2ig validator env
# validation half (run durably via systemd-run — see "Running long jobs"):
GPU=0  bash scripts/04_proteinmpnn.sh      2>&1 | tee logs/04.log # Phase 3   (8 seqs/backbone)
GPU=0  bash scripts/06_ig_validate.sh      2>&1 | tee logs/06.log # Phase IG  (thread + initial-guess + rank)
# or both at once:  GPU=0 bash scripts/run_all.sh
```

Overrides: Phase 1 takes `GPU` (default 1), `NUM_DESIGNS`, `BINDER_MIN`/`BINDER_MAX` (default 60–90); Phase 3/IG take `GPU` (default 0), `NUM_SEQ_PER_BACKBONE`/`NUM_SEQ` (default 8), and `IG_PAE_MAX`/`IG_PLDDT_MIN` (default 10 / 80). Every script is idempotent — it skips work whose output already exists, so a re-run resumes.

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
