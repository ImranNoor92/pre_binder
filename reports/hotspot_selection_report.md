# Hotspot / Site Selection for a T=3-specific AAV binder — working report
*Date: 2026-06-22. Covers the pivot away from the central (Asp111) design toward a seam-based, charge-complementary site.*

---

## 1. Background — why we moved off the central design

The previous campaign designed a single binder to **cap the centre of the AAV T=3 capsid hexamer**, steered to the only fully-exposed residue that converges from all six VP3 protomers: **Asp111** (a ~22 Å central cluster). Backbones were made with RFdiffusion (200), sequence-designed with ProteinMPNN (1,600 sequences), foldability-screened with AlphaFold2, and the best 20 were validated against the hexamer.

The result was a clean, three-way **negative**:

- **Binder quality (intrinsic):** the binders were *good proteins* — ~21 % folded confidently and to-design (pLDDT > 80, Cα RMSD < 2 Å), and ~99 % were predicted stable. (One flag: ~86 % skewed hydrophobic overall — a mild solubility/aggregation risk.)
- **Binding confidence (AlphaFold2-multimer, target held as template):** every candidate failed — i-pTM 0.21–0.30 (pass ≈ 0.5), **ipSAE = 0**, **pDockQ2 ≈ 0.009** (acceptable ≈ 0.23). The interface-specific scores (ipSAE, pDockQ2) confirmed this was *not* a large-target dilution artifact.
- **Interface energy/geometry (Rosetta InterfaceAnalyzer):** every design had a **positive (repulsive) binding energy** (dG_separated +53 to +591 REU), mediocre shape complementarity (0.36–0.60), and almost no interface hydrogen bonds (0–4), despite burying real surface area (830–1,355 Å²).

**Diagnosis:** the binders fold and make contact, but the **site is the problem**. The central all-six epitope is **flat, polar, and charged** (an Asp cluster). Binding affinity is normally driven by the hydrophobic effect plus shape complementarity, so a flat charged patch offers nothing to grip — the documented hard/failure case for de novo binder design.

This matched the project's two prior failures, which turn out to be **two different problems**:
1. **Original subunit + C3-replication:** hotspots were on *buried interior* residues (105–120), so the binder docked *through* a monomer's inner loops — useless, "AF score trash."
2. **Central design:** hotspots were on the *exposed but flat, polar* tip (Asp111) — nothing to grip.

---

## 2. Purpose

Find a **better binding site** for a **T=3-hexamer-specific** binder — one that avoids both prior failure modes — and characterise it rigorously enough to justify the choice in a methods section. Concretely, a site that is simultaneously:

- **exposed** (reachable, not an interior loop — avoids failure #1),
- **concave** (a groove to grip, not a flat tip — avoids failure #2),
- **charge-patterned** (so we can design a *charge-complementary* binder face rather than chase hydrophobicity that isn't there),
- **assembly-unique** (present in the T=3 hexamer but **not** in the natural T=1 capsid — the source of specificity).

---

## 3. New idea / brainstorming

Several threads were weighed before settling on the plan:

**(a) Supervisor input — embrace the chemistry.** Rather than hunt for absent hydrophobicity, *design for shape + charge complementarity* (e.g., a basic binder face against an acidic patch), with the hard constraint that **no buried polar/charged group is left unsatisfied** (every buried polar atom must hydrogen-bond or salt-bridge). Target geometric features — **depressions, grooves, seams, and regions that differ from the trimer/pentamer/T=1**.

**(b) Symmetric ring binder (6WVS-style / C3).** The idea of a ring-shaped, symmetric binder that matches the capsid's symmetry. Findings that shaped the decision:
- The target hexamer is **C3, not C6** (verified: a 120° rotation maps chain A onto chain C at 0.0 Å RMSD; 60° gives 12.3 Å). It is a **trimer of dimers** (A/C/E and B/D/F).
- **6WVS** is a ~200-residue de novo **TIM-barrel** (8-stranded β-barrel ring), not a C8 toroid — useful only as a "large, ring-like topology" reference.
- **RFdiffusion cannot cleanly design a symmetric binder against a target**: its symmetric mode makes free-standing oligomers, and its docs warn that symmetry + hotspots + PPI "interact weirdly / break." This is the same limitation that drove the original C3-replication route.

**(c) C3-replication (the original route) — and its trap.** Designing one subunit and replicating it by C3 previously produced copies that were **too far apart**, needing **flexible linkers** to fuse into one chain, which was **hard to make rigid** — and it failed.

**(d) The resolution — decouple "where" from "what shape," and target a seam.** Two independent questions:
- **Where** to bind → a **seam** (a groove at the junction of two adjacent monomers). A seam only exists in the assembled hexamer, so **a single binder on a T=3-unique seam is hexamer-specific by itself** — no symmetry and no trimerization required (avoiding the trap entirely).
- **What shape** the binder is → start with a **single** binder (RFdiffusion's reliable standard mode); add symmetry/multivalency only later as an optional avidity booster.

This sequencing also de-risks correctly: we have never successfully bound *any* site, so the priority is to prove a bindable, T=3-unique seam exists before stacking on the most failure-prone machinery (symmetry). Finding such a seam (Step 1 below) is a prerequisite for *every* downstream architecture.

---

## 4. Methods

### 4.1 Per-residue site descriptors (one VP3 protomer, in the assembled hexamer)
Computed on chain A of `151lp3t3_hexamer_6chain.pdb`:

- **Exposure** — relative SASA = Shrake–Rupley SASA ÷ residue maximum ASA (Tien et al. 2013), computed in the full hexamer (Biopython `ShrakeRupley`).
- **Assembly burial** — ΔrSASA = rSASA(isolated protomer) − rSASA(hexamer); high values mark inter-protomer **seams**.
- **Electrostatics (φ)** — screened-Coulomb potential at each residue's side-chain charge centre: φ = Σ k·qⱼ·exp(−r/λ_D)/(ε·r) over all ionisable groups (Asp/Glu −1, Lys/Arg +1, His +0.1, termini ±1), with ε = 78 and Debye length λ_D = 10 Å. *(Gold-standard upgrade, not used here: PDB2PQR → APBS Poisson–Boltzmann.)*
- **Concavity** — local heavy-atom density within 10 Å (scipy `cKDTree`); among exposed residues, high density = groove, low = convex tip.
- **Hexamer contact breadth** — number of distinct chains with atoms within 12 Å (≥2 = multi-protomer seam).
- **Flexibility** — **real MD RMSF**, averaged per residue over the 180-protomer T=3 capsid trajectory (GROMACS; nm → Å). *(This replaced an earlier pLDDT-based proxy.)*

### 4.2 T=1 vs T=3 assembly-uniqueness
**Structures:** T=3 = the VP3 hexamer; T=1 = the natural AAV2 capsid (**PDB 1LP3**), with the full capsid reconstructed by applying its **60 icosahedral symmetry operators** (`BIOMT`, REMARK 350) to the deposited chain A. Contacts = heavy-atom distance < 5 Å. 1LP3 numbering (80–598) was mapped to the hexamer numbering (1–504) by global sequence alignment (Biopython `PairwiseAligner`).

Two definitions of "T=3-unique" were applied (the second is preferred):

1. **Binary (presence/absence).** A residue is unique if it is at a T=3 interface (contacts chains B–F) **and not at any T=1 interface** (contacts no neighbouring capsid copy). *Limitation:* a residue at an interface in **both** assemblies is called "shared" even if it contacts a **different partner** in each — so it under-reports.

2. **Partner-pairing (contact identity).** For each residue, record the **set of partner residue numbers** it contacts in T=3 vs anywhere in the full T=1 capsid. A residue is unique if `T3_partners − T1_partners` is non-empty — i.e., it contacts a partner in the hexamer that it **never** contacts in T=1. This distinguishes genuinely new seams from shared interfaces that merely reuse the same residues.

### 4.3 Candidate site call
A residue is a candidate if it is **T=3-unique AND exposed (rSASA > 0.25) AND concave (density ≥ 60th percentile among exposed) AND charged (|φ| ≥ 60th percentile)**. Adjacent candidates (gap ≤ 3) were merged into seams.

---

## 5. Results

### 5.1 Site map quantitatively explains the prior failures
*(Figure: `reports/site_map.pdf` — exposure, electrostatics, concavity, hexamer contact breadth.)*

| residue | exposure | φ | concavity | #chains | interpretation |
|---|---|---|---|---|---|
| **Lys105** (old inner hotspot) | **0.00 (buried)** | +0.2 | high (204) | 2 | binder had to dive inside → failure #1 |
| **Asp111** (central tip) | 0.69 | **−3.3 (very acidic)** | **85 (low → convex tip)** | 5 | exposed + charged + multi-chain, but a **flat bump** → failure #2 |

The mapping objectively shows the missing ingredient was **concavity** (a groove), which we now screen for.

### 5.2 T=1 vs T=3 — the partner-pairing test reveals real unique seams
*(Figure: `reports/site_selection_combined.pdf` — 6 panels: exposure, electrostatics, concavity, MD RMSF, T=3-vs-T=1 partner counts, and the partner-pairing T=3-unique track.)*

| metric | binary | **partner-pairing (preferred)** |
|---|---|---|
| T=3-unique residues | 39 | **105** |
| exposed unique | 2 | **22** |
| candidate seams | 1 | **5 (7 residues)** |

Context: the natural T=1 capsid buries most of the VP3 surface at interfaces (365/504 residues), so the binary test is conservative; the partner-pairing test recovers genuinely new seams it misses.

### 5.3 Candidate T=3-unique seams (shortlist)
| seam | charge | exposure | MD RMSF | note |
|---|---|---|---|---|
| **482–484** | **basic (+)** | 0.54 | **0.9 Å (rigid)** | top pick — exposed, rigid, away from the failed centre; design an **acidic** face |
| **108–109** | acidic (−) | 0.48 | 1.1 Å | central; unique but adjacent to the failed Asp111 → design a **basic** face |
| 338 | acidic (−) | 0.39 | 1.0 Å | secondary |
| 144 | basic (+) | 0.26 | 0.8 Å | secondary, rigid |
| 115 | acidic (−) | 0.36 | 1.1 Å | secondary |

**Recommendation:** lead with **seam 482–484** — exposed, rigid (low MD RMSF), genuinely T=3-unique, and basic (so a charge-complementary acidic binder face), and *not* the central polar region that already failed.

### 5.4 Data and visualisation outputs
- Per-residue table (all metrics): `outputs/t1t3_site_table.csv` and `outputs/site_analysis.csv`.
- ChimeraX surface colouring (B-factor): `outputs/chainA_t3unique.pdb` (# T=3-unique partner contacts), `chainA_phi.pdb` (charge), `chainA_concavity.pdb`.
- Figures: `reports/site_map.{pdf,png,svg}`, `reports/site_selection_combined.{pdf,png,svg}`.

---

## 6. Tools, data, and methods used

**Structures / data**
- AAV T=3 hexamer model `151lp3t3_hexamer_6chain.pdb` (Cioffi & Luque 2026).
- Natural T=1 capsid: **PDB 1LP3** (AAV2), capsid rebuilt from its 60 `BIOMT` icosahedral operators.
- **MD RMSF** from a GROMACS simulation of the 180-protomer T=3 capsid (`rmsf_BB_COM_atomcenter.xvg`).

**Design / validation engines (used earlier; context for the pivot)**
- **RFdiffusion** (Complex_base) — backbone generation.
- **ProteinMPNN** — sequence design.
- **AlphaFold2** — `af2_initial_guess` (model_1_ptm) for binder-only foldability; **ColabDesign/AlphaFold2-multimer** (multimer_v3, target held as template) for interface confidence.
- **PyRosetta `InterfaceAnalyzer`** — dG_separated, shape complementarity, buried SASA, interface H-bonds, buried-unsatisfied polars.
- Interface metrics: **ipSAE** (Dunbrack et al. 2025) and **pDockQ2** (Zhu et al. 2023) — interface-specific scores robust to large-target dilution.

**Analysis libraries**
- **Biopython** — `ShrakeRupley` SASA, `ProtParam` developability descriptors, `PairwiseAligner` for T=1↔T=3 numbering.
- **scipy** — `cKDTree` (contacts/neighbour density), spatial geometry.
- **matplotlib** — figures (A4/compact, Arial-metric, vector PDF/SVG).
- *Documented but not installed (upgrade path):* PDB2PQR + APBS (Poisson–Boltzmann electrostatics); fpocket/MSMS (formal pocket/groove detection).

**Reproducibility (scripts)**
- Per-residue descriptors: `scripts/lib/site_analysis.py`
- T=1 vs T=3, binary: `scripts/lib/t1_t3_compare.py`
- T=1 vs T=3, partner-pairing (preferred): `scripts/lib/t1_t3_partner.py`

---

## 7. Next step (not yet done)
Pick the seam (lead candidate **482–484**), then design **one** binder against it with a **charge-complementary** face and **no buried unsatisfied polars** (standard RFdiffusion PPI mode → ProteinMPNN with interface charge bias → AlphaFold-multimer + Rosetta filters that explicitly reward salt bridges / penalise buried-unsat). Symmetry/multivalency remains an optional later booster, achieved geometrically rather than via RFdiffusion's symmetric-PPI mode.
