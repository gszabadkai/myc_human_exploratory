---
date: 2026-10-01
status: >
  DECLARATION - pre-data. Nothing below has been computed. Exposure-side only;
  no outcome variable is read, and none is joined.
repo: myc_human_exploratory. `myc_human_validation` and `myc_mouse` read-only.
relates-to:
  - E11, E14 (the conditioning asymmetry; compartment specificity)
  - E34 (OXPHOS-ordered in 16 of 18 cells; MYC four to seven times smaller and
    identical at both respiratory levels)
  - E35 (SATURATED; the coupling flat across genomic burden tertiles)
  - myc_mouse, the gate model: wild-type gland slope -0.012, p 0.95; the
    coupling appears only in Myc+ animals
decides:
  - The two questions, the composition controls, and the three readings,
    including the one that damages the model rather than failing to support it.
next-action: commit alone, then data note, then script unrun, then result.
---

# E36 — is the coupling present in normal breast?

## 0. Why this is the load-bearing analysis

Everything run so far has been confirmatory. **This one can come back against
the model**, which is what makes it worth running.

The mouse finds the respiratory state and the apoptotic configuration unrelated
in wild-type gland (slope −0.012, p 0.95) and coupled only where MYC is
expressed. Human tumours show the coupling irrespective of MYC level (E34) and
across the full range of genomic burden (E35). The standing reconciliation is
that **the human comparison lies entirely above the condition the mouse
manipulates**, because the mouse contrasts a pre-tumour gland with transformed
tissue whilst the human contrasts tumours with tumours.

That reconciliation rests on a premise nobody has measured. E35 removed one
alternative — that the coupling is confined to high-burden disease — and left
the premise where it was, because **a cohort that never straddles the condition
produces E35's picture whether the condition exists or not.**

Normal tissue is the only available comparison below the tumour range.

## 1. Two questions, one scoring

TCGA-BRCA carries ~113 solid-tissue-normal samples.

**Q1. Is the OXPHOS-to-configuration coupling present in normal breast?**
The human analogue of the mouse wild-type gland.

**Q2. Where do tumour MYC levels sit relative to normal breast?**
This tests the premise directly. It has never been measured: every MYC quantity
in this arm is scored **within** a tumour cohort, and GSVA is cohort-relative by
construction, so "MYC-low" has carried no information about absolute level.

**Both require a single joint scoring of tumours and normals**, which
cohort-relativity requires anyway. Score once, over the combined set, and say so.

## 2. Composition is the design problem, not a caveat

Normal-adjacent breast is largely stroma and adipose; a tumour is epithelium.
**Any difference in a mitochondrial-to-apoptotic relationship could be tissue
composition rather than biology.** Three mitigations, all declared now:

**(a) Epithelial content, not ABSOLUTE purity.** ABSOLUTE is defined for tumours
and does not transfer to normal tissue. Use an epithelial-content estimate
computable in both (a deconvolution estimate, or an epithelial marker score)
and say which was used and why. Report the comparison with and without it.

**(b) A positive control pathway, from E14's thirty comparators.** Choose a
mitochondrial programme whose cytosolic half has no reason to differ between
normal and tumour, declared **before** the fit. **If it differs as much as
apoptosis does, the analysis is measuring composition** and the apoptotic result
is uninterpretable. This is a stop condition, not a footnote.

**(c) Report the normal-sample n at every step.** At ~113 the normal-side
interval will be wide, and a wide interval around zero is not evidence of
absence.

## 3. Readings, fixed now

| reading | condition | meaning |
|---|---|---|
| **ABSENT IN NORMAL** | the coupling does not clear zero in normals and does in tumours, **and** the positive control does not differ | supports the conditionality: the coupling is a property of transformed tissue. The predicted outcome |
| **WEAKER IN NORMAL** | present in normals but materially smaller, control clean | partial support. The coupling is graded rather than switched |
| **PRESENT IN NORMAL** | present at comparable strength | **the mouse wild-type result and the human normal result disagree.** The conditionality claim is damaged, not merely unproven |
| **UNINTERPRETABLE** | the positive control differs between normal and tumour | composition. No reading taken |

**PRESENT IN NORMAL is a real possibility and it is recorded here so that it
cannot be reframed afterwards.** If it occurs, the manuscript says the two
species disagree and that the mouse carries the conditionality by intervention
whilst the human data do not corroborate it.

## 4. For Q2, the MYC comparison

Report the MYC activity distributions of normals and of each tumour quadrant on
the common scale, as medians with IQR and as the percentage of tumours in each
quadrant exceeding the normal median. All three estimators (`M_a`, `M_b`, and
`FELSHER__PROLIFSTRIP`), since the proliferation-stripped version is the one that
matters for a claim about MYC rather than growth.

**No verdict on Q2.** It is descriptive, and it either does or does not support
the premise. If a substantial fraction of Q1 and Q2 tumours sit at or below the
normal median, the threshold reconciliation fails and that must be written.

## 5. What this cannot do

- Normal-adjacent tissue is not normal tissue: it is field-adjacent, from breasts
  that grew a tumour.
- n ≈ 113, unreplicated. SCAN-B carries no normals.
- Cross-sectional, bulk, transcript-level. The configuration is not lethality.
- It cannot show the barrier operated. It can show whether the configuration is
  distributed as a barrier would predict.

## 6. What a reading licenses

**ABSENT or WEAKER IN NORMAL** licenses one sentence: that the coupling observed
across human tumours is not a property of breast tissue generally. That is the
human analogue of the mouse gland result and it is the first direct support for
the conditionality from human data.

**PRESENT IN NORMAL** licenses a limitation, stated in the Discussion rather than
buried: the human data do not corroborate the conditionality, and the mouse
carries it alone, by intervention.

**No reading licenses a treatment-selection claim**, and none concerns killing.
