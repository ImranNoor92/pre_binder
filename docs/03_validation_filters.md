# AF2 validation filter rationale

The Phase 2b and Phase 4 AF2-multimer validation steps each apply four filters. This document explains why each threshold is where it is.

---

## Filter 1: Binder pLDDT > 0.70

**What it measures:** AF2's per-residue confidence in its prediction of the binder's structure. Averaged across all binder residues.

**Why 0.70:** Below 0.70 means AF2 doesn't trust the predicted fold — the sequence may not fold into the designed shape. Above 0.80 is "high confidence" but rare for de novo designs; 0.70-0.80 is the sweet spot for design pipelines.

**What rejection looks like:** A loopy, disordered-looking binder where multiple regions are predicted at <0.50.

---

## Filter 2: Interface pTM > 0.65

**What it measures:** AF2's confidence that the *interface* between binder and target is correctly placed (separate from the binder's own fold). Computed only over residue pairs across the binder-target interface.

**Why 0.65:** This is more lenient than the binder pLDDT threshold because interface prediction is intrinsically harder than monomer fold prediction. Most published binder design papers use 0.60-0.70 here.

**What rejection looks like:** Binder is well-folded but AF2 doesn't think it's actually engaging the target — the binding mode is uncertain or the binder is "hovering" without making confident contacts.

---

## Filter 3: All 3 subunits in contact (per-subunit interface SASA > 200 Å²)

**What it measures:** For each of the three binder subunits, the amount of solvent-accessible surface area that gets buried upon binding. Computed by SASA difference: SASA(binder alone) - SASA(binder bound to target).

**Why 200 Å²:** A small but nontrivial interface. Below 200 Å² typically means a glancing contact (a single salt bridge or one helix tip touching). 500-1500 Å² is a "normal" protein-protein interface.

**Why all 3 must hit threshold:** This is the **hexamer-specificity filter**. If only 1 or 2 subunits make real contacts (with the third just dangling), the design effectively becomes a dimer-binder with extra trailing structure. The whole point of the C3-symmetric design is that the binder *requires* all three target chains to be present. This filter enforces that requirement.

**What rejection looks like:** Two subunits dock cleanly onto A and C, but the third subunit (intended for E) is rotated away or has tiny contact area — meaning the binder could still function as a dimer-binder.

---

## Filter 4: RMSD vs. designed structure < 3.0 Å

**What it measures:** Backbone-RMSD (Cα atoms) between the AF2 prediction and the RFdiffusion + C3-replicated input structure. Computed after global alignment on the target.

**Why 3.0 Å:** AF2 nearly always shifts the design slightly during prediction. < 2.0 Å is excellent agreement; 2.0-3.0 Å is acceptable; > 3.0 Å suggests AF2 is predicting a fundamentally different conformation than what RFdiffusion designed (and we should not trust that the design's properties carry over).

**What rejection looks like:** RFdiffusion designed the binder as an alpha-bundle bridging A and C; AF2 predicts it as an extended loop wrapping around the target — same sequence, totally different conformation. The design's planned properties don't apply to the AF2-predicted conformation.

---

## Putting them together

A design that passes all 4 means:
- ✓ The binder folds confidently (Filter 1)
- ✓ It binds the target confidently (Filter 2)
- ✓ It engages all three target chains, not just one (Filter 3)
- ✓ Its predicted conformation matches its designed conformation (Filter 4)

That's a defensible "this design should work as intended" claim.

---

## Suggested per-design score for final ranking

After all filters pass, rank by a combined score that rewards stronger interfaces and tighter agreement:

```
score = interface_pTM * binder_pLDDT * (1 / max(RMSD, 0.5)) * (sum_interface_SASA / 1000)
```

Top-10 by this score is what would go to wet-lab.

---

## The acid test: dimer-only re-prediction

Above all the filters, the **most important** validation for hexamer-specificity is to re-run AF2 on each surviving design against **only one dimer pair** of the target (e.g., chains A + E only) and confirm:

- Interface pTM drops by ≥ 0.15 (i.e., AF2 is much less confident it binds the dimer)
- At least one subunit has near-zero interface SASA (the subunit that would've engaged C, D, E, or F has no partner)

If the design *still* shows interface pTM > 0.65 against the dimer alone, it's not actually hexamer-specific in practice — it's a dimer-binder that happened to fold as a trimer. Discard those.

This dimer-only validation is added as an explicit script step at the end of Phase 4.
