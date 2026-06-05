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
| 2b | `scripts/03_af2_validation.sh` | AF2-multimer predict hexamer+trimer; apply 4 filters (binder pLDDT > 0.70, interface pTM > 0.65, per-subunit interface SASA > 200 Å², RMSD vs design < 3 Å) | AlphaFold (`.venv-af2`) | ⚠️ **skeleton** — AF2 call + filters not yet wired |
| 3  | `scripts/04_proteinmpnn.sh` | Sequence design with the three subunits **tied** (preserves homotrimer symmetry); target chains fixed | ProteinMPNN (`.venv-rfd-gpu`) | ⚠️ env wired; tied/fixed/chain JSONLs not yet generated |
| 4  | `scripts/05_af2_revalidation.sh` | Re-predict each MPNN sequence; re-apply the 4 filters; rank survivors | AlphaFold (`.venv-af2`) | ⚠️ **skeleton** — reuses Phase 2b AF2 |
| Acid test | (in `05`) | Re-predict each survivor against a single dimer pair (A+E only); a true hexamer-specific binder must lose ≥ 0.15 interface pTM | AlphaFold | ⚠️ described, not yet implemented |

**Runnable today:** Phases 1 and 2a. Phases 2b/3/4 need their AF2 invocation, the MPNN tied-positions JSONL builder, and the dimer-only acid test implemented before an end-to-end run. The reference orchestration to adapt from is `/data/rfdiffusion/trial_2B/master_pipeline/`.

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

## Quick start

```bash
cd /data/binder_software/pre-binder

bash   scripts/00_check_env.sh                                 # ~30 s — must be all ✓
bash   scripts/01_pilot_rfdiffusion.sh  2>&1 | tee logs/01.log # ~2–4 h on 1 GPU → 10 backbones
python scripts/02_trimerize_replicate.py                       # ~1 min  → trimerized PDBs
# Phases 2b–4 require the AF2/MPNN wiring noted above before they will run end-to-end:
bash   scripts/03_af2_validation.sh     2>&1 | tee logs/03.log
bash   scripts/04_proteinmpnn.sh        2>&1 | tee logs/04.log
bash   scripts/05_af2_revalidation.sh   2>&1 | tee logs/05.log
```

Useful overrides for Phase 1: `GPU=0` (default 1), `NUM_DESIGNS=N`, `BINDER_MIN`/`BINDER_MAX` (default 60–90). Each script is idempotent — it skips work whose output already exists.

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
├── scripts/
│   ├── 00_check_env.sh            # environment + inputs verification
│   ├── 01_pilot_rfdiffusion.sh    # Phase 1
│   ├── 02_trimerize_replicate.py  # Phase 2a
│   ├── 03_af2_validation.sh       # Phase 2b (skeleton)
│   ├── 04_proteinmpnn.sh          # Phase 3
│   └── 05_af2_revalidation.sh     # Phase 4 (skeleton)
├── inputs/
│   ├── 1lp3_hexamer_trimmed_fixed.pdb
│   ├── hotspot_residues.txt
│   └── target_pdb_link.md
├── outputs/                       # generated; git-ignored
└── logs/                          # generated; git-ignored
```

`outputs/` and `logs/` are git-ignored (regenerable, potentially large AF2 artifacts). Everything needed to reproduce a run is committed.
