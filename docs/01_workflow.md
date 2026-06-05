# Workflow — phase-by-phase technical detail

This document explains *what each phase does, why it does it that way, and what the output looks like*. Read after the top-level README.

---

## Phase 1: RFdiffusion single-subunit binder

**What it does:** Generate ~10 backbone PDBs where each one is a small (50–80 residue) binder docked onto chain A of the hexamer, with hotspot bias toward residues 105-115.

**Why this and not symmetric directly:** Per the RFdiffusion README, symmetric mode + hotspots is buggy. We get a more reliable result by designing one subunit first and then composing them symmetrically post-hoc.

**Why chain A specifically:** All 6 chains have equivalent 105-115 patches (just rotated by C3 around the axis). Chain A is one of the "upper" symmetry-class chains (A, D, F). Using A means we'll later replicate by 120° to get F+D copies positioned at the C3-related sites. (Or we can target one of the "lower" class — B, C, E — by replicating from B; same result by symmetry.)

**Input:**
- Target PDB: `/data/binder_software/pre-binder/inputs/1lp3_hexamer_trimmed_fixed.pdb` (the full hexamer — RFdiffusion needs to see the steric context of B, C, D, E, F when designing the binder against A)
- Hotspots: `A105, A107, A109, A111, A114, A115` (the exposed residues from the SASA analysis, skipping the buried 108 and 110)
- Binder length: range 60-90 residues (small subunits — each only needs to engage ~11 hotspot residues)

**Output:** 10 PDB files in `outputs/01_rfdiffusion_pilot/`. Each contains:
- All 6 chains of the target (unchanged)
- 1 new chain (binder subunit) docked on chain A's 105-115 face
- Backbone-only (sequences are placeholder glycines — MPNN comes later)

**Quick sanity check after this phase:** All 10 backbones should have the binder reasonably positioned near A105-115 with no severe clashes. If the success rate is low (say, <50% of designs are clean), the epitope itself may be problematic and we revisit hotspot choice before scaling up.

---

## Phase 2: C3 replication + AF2 validation

**Phase 2a — geometric replication (Python, ~1 min):**

For each Phase 1 backbone:
1. Identify the binder subunit's atoms.
2. Compute the hexamer's 3-fold axis. This is the line through (0, 0, z_avg) along the Z direction (verified from the geometry: 3 dimer-pair centroids are at exactly 120° intervals around this axis).
3. Generate two additional copies of the binder subunit by rotating 120° and 240° around the 3-fold axis.
4. Fuse the three subunits into a single polypeptide via flexible linkers:
   - Linker: `(GGGGS)n` — 5-15 residues depending on N-to-C terminus distances after rotation. Compute the distance from the C-terminus of subunit_1 to the N-terminus of subunit_2 after rotation, then use 1 linker residue per ~3.5 Å as a rough rule.
   - If the geometry doesn't allow fusion (termini point the wrong way), skip this backbone — its topology isn't suitable for a fused trimer. You'd need to express it as 3 separate chains instead.
5. Write a new PDB: target hexamer + fused trimeric binder.

**Output:** PDBs in `outputs/02a_trimerized/`. Each represents one candidate fused-trimer design.

**Phase 2b — AF2-multimer validation (~10-15 min per design):**

For each trimerized backbone, run AF2-multimer prediction with:
- Chains 1-6: target hexamer
- Chain 7: the fused trimeric binder

Then evaluate:

| Filter | Threshold | Why |
|--------|-----------|-----|
| Binder pLDDT | > 0.70 | The binder must be well-folded |
| Interface pTM (binder ↔ target) | > 0.65 | Confident interface prediction |
| All 3 subunits in contact | per-subunit interface SASA > 200 Å² | This is what makes it hexamer-specific — if only 1-2 subunits contact, it could bind a dimer too |
| RMSD vs. designed trimer | < 3.0 Å | AF2's prediction matches the symmetric design |

**Output:** PDBs that pass filters → `outputs/02b_af2_validated/`. Designs failing any filter → `outputs/02b_af2_rejected/` (kept for debugging).

**Expected yield:** ~30-50% of Phase 2a designs pass Phase 2b. If yield is <20%, the C3 replication geometry may be incompatible with how AF2 wants to fold the trimer — would need to reduce linker length or try a different subunit length distribution.

---

## Phase 3: ProteinMPNN sequence design

**What it does:** Replace the placeholder sequences (post-RFdiffusion) with realistic amino-acid sequences optimized for stability and the backbone shape.

**Key MPNN flags for this case:**
- `--fixed_chain_list "A B C D E F"` — lock the target sequence (we're not redesigning the hexamer)
- `--design_only_positions` — restrict design to the binder chain
- `--tied_positions` — equivalent positions in the three binder subunits must have the same amino acid (preserves the symmetry — each subunit is identical at the sequence level too)
- `--num_seq_per_target 8` — sample 8 sequences per backbone (gives MPNN room to find a good one)
- `--sampling_temp 0.1` — low temperature for confidence over diversity

**Output:** FASTA files with 8 sequences per backbone in `outputs/03_mpnn_sequences/`, plus repacked PDBs (each backbone with its top sequence threaded back on).

**Note on tied positions:** This is the most important MPNN constraint for this design. Without tying, the three subunits could diverge into different sequences, which (a) breaks the C3 symmetry assumption and (b) defeats the whole point — the binder would no longer be a homotrimer. The trial_2B pipeline at `/data/rfdiffusion/trial_2B/master_pipeline/` has a working example of tied positions for the analogous "A/C and B/D windows" case.

---

## Phase 4: AF2 re-validation on MPNN sequences

**Why again?** MPNN scores sequences against the backbone, but doesn't predict whether AF2 will fold them back into that backbone. ~20-40% of MPNN sequences fail to fold to the design when AF2 predicts them.

**Same AF2 setup as Phase 2b**, applied per MPNN sequence:
- Target hexamer + binder (with MPNN-assigned sequence)
- Same 4 filters

**Final output:** `outputs/04_final_ranked/` containing PDBs that pass all filters in both Phase 2 and Phase 4, sorted by combined score (interface pTM × binder pLDDT × inverse RMSD).

**Expected funnel** (if everything goes well):
- 10 Phase 1 backbones (pilot)
- → ~5-7 trimerizable in Phase 2a
- → ~3-5 pass AF2 validation in Phase 2b
- → ~30 MPNN sequences in Phase 3 (~8 per validated backbone)
- → ~10-15 pass final AF2 in Phase 4

For 50 final designs, scale Phase 1 to ~50 backbones (Phase 1 has the lowest computational cost per design).

---

## What success looks like

A pre-binder result we'd be confident shipping for wet-lab validation:
- ≥3 distinct backbone topologies (so we're not betting on one design)
- All pass the 4 AF2 filters
- All have inter-subunit interface area > 500 Å² (real interfaces, not glancing contacts)
- All show < 2.5 Å backbone RMSD between the RFdiffusion design and the final AF2 prediction
- ≥1 of them shows a clear "hexamer vs dimer" specificity in a final AF2 dimer-only re-prediction (binder pLDDT or i_pTM should drop significantly when target is reduced to a single dimer pair)

The last criterion is the **acid test** for hexamer-specificity and should be added as a separate validation step after Phase 4.
