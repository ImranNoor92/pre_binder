# RFdiffusion contig syntax for this target

The "contig string" is RFdiffusion's way of describing what to build. This document explains the exact contig we'll use for Phase 1 (single-subunit binder against chain A of the hexamer).

---

## Anatomy of the contig

For a single-chain PPI binder design against a multi-chain target, the contig has this structure:

```
contigmap.contigs="[<chain1_range>/0 <chain2_range>/0 ... <binder_length_range>]"
```

Where:
- `<chainX_range>` — a target chain range, e.g. `A70-150` means "chain A residues 70 to 150, as fixed positions (RFdiffusion will not change them)"
- `/0` — a chainbreak (no peptide bond) between this chain and the next. Without this, RFdiffusion would interpret the next contig element as the same chain.
- `<binder_length_range>` — e.g. `60-90` means "generate a new chain of 60-90 residues" (RFdiffusion samples a length from this range each iteration)

The last element (no chain letter, just a range) becomes the new binder chain.

---

## Our specific contig

For the hexamer target with all 6 chains kept as context:

```
contigmap.contigs="[A70-150/0 B70-150/0 C70-150/0 D70-150/0 E70-150/0 F70-150/0 60-90]"
```

Decoded:
- `A70-150/0` — chain A residues 70-150 (fixed). The full range of chain A in our PDB.
- `B70-150/0`, `C70-150/0`, `D70-150/0`, `E70-150/0` — same for chains B through E.
- `F70-150/0` — chain F residues 70-150 (fixed). The `/0` here is technically not necessary on the last target chain, but harmless.
- `60-90` — generate a new chain of 60-90 residues. This is the binder. RFdiffusion will pick a length from this range each design.

The binder chain ID is auto-assigned by RFdiffusion (usually 'G' since A-F are taken).

---

## Hotspots

Hotspots tell the diffusion model "bias the binder to make contacts with these specific residues":

```
ppi.hotspot_res="[A105,A107,A109,A111,A114,A115]"
```

We're using **6 residues on chain A only** — the ones from our SASA analysis that have rSASA > 30% (well-exposed). The fully or partially buried residues (108, 110, 113) are omitted because the binder can't physically reach them without burying through the dimer interface.

**Why not all 6 chains?** Because:
1. RFdiffusion's PPI hotspot mode is per-chain biased; with 66 hotspot residues spread across 6 chains, the bias diffuses and the model has trouble committing to one binding face.
2. We *want* the binder to commit to chain A here — Phase 2 handles the C3 replication to engage the other chains.

---

## Noise scale flags

```
denoiser.noise_scale_ca=0
denoiser.noise_scale_frame=0
```

These set noise to zero during denoising. The RFdiffusion paper found that for PPI design, zero noise gives higher quality designs at the cost of less diversity. Since we'll have many backbones from Phase 1 to filter through anyway, we prioritize quality.

If you find designs are too similar to each other (low diversity), bump these up to `0.5` for more exploration.

---

## Number of designs and timing

```
inference.num_designs=10        # pilot
inference.num_designs=100       # full run, after pilot validates the setup
```

Per RFdiffusion's typical timing:
- Each design takes ~5-15 min on the RTX 6000 Ada (the bigger the binder length, the longer)
- 10 designs ≈ 1-2.5 hours
- 100 designs ≈ 10-25 hours

These run sequentially within one RFdiffusion process. To parallelize, launch two processes on two GPUs with non-overlapping seed ranges.

---

## Full Phase 1 command (preview, not for execution yet)

```bash
cd /data/binder_software/pre-binder

CUDA_VISIBLE_DEVICES=0 \
/data/rfdiffusion/.venv-af2/bin/python /data/rfdiffusion/scripts/run_inference.py \
  inference.output_prefix=outputs/01_rfdiffusion_pilot/design \
  inference.input_pdb=/data/binder_software/pre-binder/inputs/1lp3_hexamer_trimmed_fixed.pdb \
  inference.num_designs=10 \
  'contigmap.contigs=[A70-150/0 B70-150/0 C70-150/0 D70-150/0 E70-150/0 F70-150/0 60-90]' \
  'ppi.hotspot_res=[A105,A107,A109,A111,A114,A115]' \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
```

**Caveat about the venv:** The active rfdiffusion venv at `/data/rfdiffusion/.venv-af2` is named "af2" but is the one your earlier rfdiffusion job (trial_2B) used. Verify by running `00_check_env.sh` before launch. If RFdiffusion has its own venv elsewhere (e.g. `.venv-rfd`), use that instead.

---

## How to read the output PDBs

Each Phase 1 output looks like:

```
outputs/01_rfdiffusion_pilot/design_0.pdb
outputs/01_rfdiffusion_pilot/design_0.trb        # trajectory data (binary)
outputs/01_rfdiffusion_pilot/design_1.pdb
...
```

Open `design_0.pdb` in PyMOL or ChimeraX. You should see:
- Chains A-F: the original hexamer
- Chain G (or similar): the new binder, docked onto chain A's 105-115 face

Backbone atoms only (C, N, CA, O). All residue names are `GLY` (placeholder). The actual sequence comes from MPNN in Phase 3.

The `.trb` files contain extra metadata (per-residue confidence, sampled length, RNG seed, etc.) — usable for filtering but not needed for the basic pipeline.
