---
date: 2026-09-30
status: >
  DECLARATION - pre-data. Nothing below has been computed. Exposure-side only;
  no outcome variable is read, and none is joined.
repo: myc_human_exploratory. `myc_human_validation` and `myc_mouse` read-only.
relates-to:
  - E11 (the conditioning asymmetry; the machinery is ordered by OXPHOS, not MYC)
  - E10, E14 (the twelve transcripts; compartment specificity)
  - E33 (NEITHER SEPARATES; the configuration leans MB1, as OXPHOS does)
  - myc_human_validation/docs/2026-08-30_STATE_frozen_and_H4_buffer_declaration.md
decides:
  - The four contrasts, their predicted directions, and the reading rule, all
    fixed before any group mean is computed.
  - What each verdict would mean for the therapeutic stratification, written
    now so that reading is not constructed after the fact.
next-action: commit alone, then data note, then script unrun, then result.
---

# E34 — the apoptotic configuration across the MYC x OXPHOS quadrants

## 0. Why this exists, stated plainly

E11 established the central human finding through a conditioning ladder on two
continuous axes: removing MYC left the OXPHOS ordering intact (0.453 to 0.485
and 0.489 to 0.525), removing OXPHOS abolished the MYC ordering (0.187 to
−0.043 and 0.137 to −0.058). E33 reached the same conclusion through an
orthogonal variable.

**This analysis is that finding in a form a reader can see**, and it isolates
one cell the continuous treatment cannot: **MYC-low with OXPHOS-high**. That
cell is what the therapeutic argument in the Discussion turns on and `STATE`
cannot supply it, because `STATE` level 1 collapses across OXPHOS.

It is therefore **largely confirmatory and presentational**, and the result note
must say so. What is new is the Q2 cell and the group-level form.

**It replaces the bicluster route deliberately.** That route required the reader
to accept MCbiclust, a fork assignment produced by another study, a UF/LF/none
grouping with roughly 600 unassigned samples, and a PAM50 that agrees with ours
on 71.7% of 920. The quadrant requires none of it.

## 1. The variable

Built from the frozen constructor's own MYC and OXPHOS calls: split `STATE`
level 1 by the OXPHOS call it already computes, merge levels 3 and 4.

```
Q1 MYC-low  OXPHOS-low      Q2 MYC-low  OXPHOS-high
Q3 MYC-high OXPHOS-low      Q4 MYC-high OXPHOS-high
```

**Assert before anything else** that Q3 equals `STATE` level 2 exactly and Q4
equals levels 3 + 4 exactly, on both instruments. E33 verified this in TCGA;
verify again wherever it is rebuilt.

**The quadrant needs no BUFFER**, which makes it more portable than `STATE`: it
carries none of the tertile-fallback problem (kappa 0.221 against copy number)
into any cohort.

**Both cohorts, both instruments.** TCGA and SCAN-B, GSVA and mitoPPS. In TCGA
the two instruments agreed on 85% of quadrant assignments (E33); report the
equivalent in SCAN-B. Never score SCAN-B without `scanb_pheno.rds$symbol_map`.

## 2. Readouts

| | |
|---|---|
| **primary** | the configuration composite, as E11 defines it |
| reported always | trigger arm (`BBC3`, `BID`, `BIK`, `BAD`) and guardian arm (`BCL2L1`, `MCL1`), separately |
| reported, not the verdict | the twelve transcripts individually. **No per-gene FDR** |

## 3. The four contrasts, and what each is for

| contrast | what it asks | predicted |
|---|---|---|
| **Q2 against Q4** | **THE ROW.** Same OXPHOS status, different MYC | **no difference**, or much smaller than the OXPHOS contrasts |
| Q1 against Q3 | control. Same OXPHOS status, different MYC | no difference |
| Q1 against Q2 | the OXPHOS effect where MYC is low | **positive** |
| Q3 against Q4 | the OXPHOS effect where MYC is high | **positive** |

## 4. Reading rule, fixed now

No equivalence testing; the rule is a comparison of magnitudes.

| | condition |
|---|---|
| **OXPHOS-ORDERED** | both OXPHOS contrasts positive with CI excluding zero, in **both cohorts on both instruments**, **and** both MYC contrasts smaller in absolute value than both OXPHOS contrasts |
| **MYC-DEPENDENT** | Q4 exceeds Q2 with CI excluding zero and of comparable magnitude to the OXPHOS contrasts |
| **MIXED** | anything else, including cohort or instrument disagreement |

**OXPHOS-ORDERED is the predicted outcome** given E11 and E33. A confirmatory
analysis that confirms is not a null result, but it is also not a new finding,
and the note must not present it as one.

## 5. Adjustment, and the confounder that is worse here than elsewhere

Q2 and Q4 differ sharply in subtype: in TCGA, Q2 is 110 LumA of 180 whilst Q4 is
86 Basal and 80 LumB of 312. **So the primary contrast is heavily
subtype-confounded**, and unusually the confounding runs *against* the predicted
null — subtype would tend to create a Q2-versus-Q4 difference where the model
predicts none.

- **PAM50 adjustment is mandatory**, and within-subtype estimates are reported
  alongside as the result rather than as a robustness check.
- Purity and leukocyte fraction as in E11.
- **Proliferation companion, mandatory**, not optional: the primary repeated
  with `PROLIF_DISJOINT` added. Reported beside the primary, never chosen
  between.

## 6. What each verdict means for the therapeutic stratification

Fixed now, so that neither reading is constructed after the numbers arrive.

**If OXPHOS-ORDERED.** The configuration is present in MYC-low, OXPHOS-high
tumours to the same degree as in MYC-high ones. The human transcriptome
therefore gives **no support** to the claim that the standing pro-apoptotic
state is MYC-restricted. That claim rests on the mouse, where MYC-dependence is
demonstrated by intervention, and the Discussion must say so explicitly. This is
a **limitation to state, not a contradiction**: transcript abundance is not a
functional readout, which the sibling study's null on functional BH3 assays
already established.

**If MYC-DEPENDENT.** The human data corroborate the MYC restriction and the
stratification gains human support. Given E11 and E33 this would be surprising
and would need the subtype-stratified and proliferation-adjusted versions to
hold before it is believed.

**Neither verdict licenses a treatment-selection claim.** No cohort available
carries an OXPHOS-directed or BH3-directed intervention.

## 7. What this cannot do

Cross-sectional, descriptive, bulk transcriptomes, two cohorts, in tumours that
have already passed the barrier. It describes where the configuration is found.
It cannot show that the barrier operated, and it cannot speak to function.
