---
date: 2026-10-22
status: >
  PROPOSAL, nothing built. No script written, none run, nothing committed.
  Every plotted value below was read from a saved object or table
  (results/e34_quadrant_configuration.rds, outputs/tables/E40_pcr_meta.csv,
  E43/E44 notes and E44_contrast_metabric.csv). Exponentiations were done
  in-session from those values. ONE number was computed in-session and is in
  no note: rho(OX, PROLIF) in METABRIC (section 4.2); it is not citable until
  recorded.
authority: >
  2026-10-04_handoff.md section 2 remains the authority on what a panel may
  assert. v5 section 3 (2026-10-08_section4_narrative_v5_submission_draft.md)
  is the specification this amends. Nothing here relaxes either.
posture: EXPLORATORY. Display only. No panel proposed here requires a new fit.
---

# Fig. 4g-4i: the options, the narrative they serve, and a proposed final specification

## 0. Summary

**All three panels can be built from saved objects with no new fit**, provided
two parts of v5's specification change. As written, neither can be built:

1. **4g's row blocks.** v5 asks for "MitoCarta members first and cytosolic
   effectors second". But **11 of the 12 configuration genes are in Human
   MitoCarta 3.0**; only `BMF` is not. Read from
   `data/mitocarta_human/Human.MitoCarta3.0.xls`, sheet A. The "cytosolic half"
   in paragraph 3 is the **24-gene cytosolic half of E14's 44-gene machinery**
   (`FAS`, `TNF`, `FADD`, `XIAP`, `APAF1`, ...). Apart from `BMF`, none of those
   genes is among the twelve. **Proposal: block the rows by function
   (pro-apoptotic, anti-apoptotic)**, and move the compartment argument to ED
   Fig. 7.
2. **4i's GSE25066 DRFS row.** **Declaration 15.4 excludes it** ("If it is ever
   computed it is reported whatever it shows"). E43 never saved its PROLIF
   term, and E44 records it as "absent, by declaration". **Proposal: replace it
   with the two declared comparators of the METABRIC primary**, which the
   narrative already cites ("in the primary stratum and in both declared
   comparators").

Three further changes:

- **One axis convention for both forests.** v5's 4i axis ("log hazard ratio
  ... 1 marked") contradicts itself.
- **4h keeps `m0` on the main panel.** v5's fallback of moving it to ED Fig. 8
  should not be used.
- **4g's strip carries both `m1` and `m3`.** E34's declaration requires the
  proliferation companion "beside the primary, never chosen between".

**Once drawn, the panels would contradict six sentences in v5** (section 4).
**Five decisions are the author's** (section 5).

---

## 1. What exists

| source | what it offers | standing |
|---|---|---|
| v5 section 3 | specification for 4g, 4h and 4i, plus a two-panel fallback (keep 4g; merge 4h and 4i, response at about 2x the width) | the specification. Nothing is built |
| `2026-10-04_handoff.md` 2.1 | the two-exposure panel (OX and PROLIF, both endpoints, same models), with the caption superseded by E44's tested contrast | absorbed into v5's 4h and 4i |
| `docs/figures/2026-10-01_E34_configuration_by_quadrant.png` | point-and-interval adjusted means for the composite, trigger arm and guardian arm. x is MYC low/high, colour is OXPHOS; both cohorts, GSVA, `m1` | "built as a candidate main-figure panel" (`2026-09-30_e34_data.md:248`). **Its caption carries a verdict and interpretation, which v5's legend rule forbids**, and it has no gene rows and no `m3` |
| E33, E35, E36, E38 figures | quadrant map, burden coupling, normal comparison, diagnostics | not 4g material. E36 licenses nothing, E38 is methodological only, and E35 is at most ED Fig. 7 |
| E39-E44 | **no figure exists.** The numbers are in `results/e40_ladder.rds`, `e43_proliferation_coefficient.rds` and `e44_contrast_test.rds` | 4h and 4i start the "write the key figure twice" convention |

Earlier narrative drafts (v2, v4 and the 2026-09-24 figure plan) exist only in
chat. They were not consulted. v5 supersedes v4 by its own frontmatter.

**One design rule carries over from the E34 figure.** It joins no quadrant
pairs with a line, because "NO PRODUCT TERM is fitted". The proposal below
keeps that rule: no connector between Q1 and Q3 or between Q2 and Q4.

---

## 2. Fig. 4g: the configuration by quadrant

**Source.** `results/e34_quadrant_configuration.rds$groups`, `subtype == "all"`.
It already holds adjusted means for all 12 genes, the composite and both arms,
for both cohorts, both instruments, and `m1` and `m3`. **No new fit is needed.**
**Snapshot:** `data/from_validation/` (E34). The joint snapshot is not used.

**Layout, kept from v5.**

- A heatmap with a strip beneath it, sharing the column axis.
- Eight columns: TCGA | SCAN-B, each ordered Q1, Q3 ‖ Q2, Q4.
- A thin rule separates the respiratory levels and a wider gap separates the
  cohorts.
- With this order, the MYC pairs sit side by side within each respiratory
  level. That matters, because **the MYC contrasts are the informative ones**:
  the OXPHOS contrasts are near-guaranteed by construction (E34 1.2).

**Changes.**

1. **Rows are blocked by function, not by compartment.**
   - Pro-apoptotic: `BAD`, `BIK`, `BBC3`, `BID`, then `BCL2L11`, `BMF`,
     `PMAIP1`.
   - Anti-apoptotic: `BCL2L1`, `MCL1`, then `BCL2`, `BCL2L2`, `BCL2A1`.
   - The order is fixed, not data-driven: composite members come first in each
     block.
   - A narrow annotation column marks the six composite genes with their sign
     (+ or −). This makes v5's "the strip is a different object from the rows"
     visible on the panel.
2. **Cell value.** The `m1` adjusted mean in within-cohort SD units (Blom),
   which is what v5 calls "standardised difference from the cohort mean". Use
   one diverging scale, white at 0, with symmetric limits of about ±0.6. The
   largest |value| is 0.50 for genes and 0.58 for the composite.
3. **Strip.** The composite per column: the `m1` estimate as a filled point with
   a 95% CI, and `m3` (+ `PROLIF_DISJOINT`) as an open point dodged beside it.
   Do not join the quadrants with lines.
4. **Instrument.** **GSVA**, chosen because the outcome cohorts in 4h and 4i
   have no mitoPPS, so one instrument then runs through 4g-4i. **The choice is
   not made on the result.**
   - **It does flatter the narrative.** E34 2.2 orders GSVA below mitoPPS on
     the MYC/OXPHOS ratio in every comparable cell.
   - **Two consequences follow.** The same layout on mitoPPS goes in ED Fig. 7.
     The text must say the MYC contrast is larger on mitoPPS.

**What the panel will show.** These are the `m1` GSVA values the cells and strip
will carry, read from the object:

| | TCGA Q1 | Q3 | Q2 | Q4 | SCAN-B Q1 | Q3 | Q2 | Q4 |
|---|---|---|---|---|---|---|---|---|
| composite `m1` | -0.58 | -0.44 | +0.41 | +0.40 | -0.36 | -0.24 | +0.23 | +0.37 |
| composite `m3` | -0.65 | -0.37 | +0.35 | +0.47 | -0.44 | -0.19 | +0.19 | +0.45 |
| `BAD` | -0.48 | -0.50 | +0.42 | +0.29 | -0.33 | -0.34 | +0.38 | +0.31 |
| `BCL2L1` | -0.34 | -0.33 | +0.41 | +0.26 | -0.23 | -0.18 | +0.24 | +0.21 |
| `BID` | -0.44 | -0.12 | +0.16 | +0.35 | -0.16 | +0.02 | -0.05 | +0.15 |
| `PMAIP1` | +0.18 | +0.04 | +0.08 | -0.19 | +0.09 | 0.00 | +0.09 | -0.11 |
| `BCL2L11` | +0.37 | +0.20 | -0.39 | -0.16 | -0.02 | +0.12 | -0.02 | -0.04 |
| n (`m1`) | 293 | 170 | 167 | 308 | 1,025 | 545 | 546 | 1,027 |

The colour change falls cleanly at the respiratory rule for `BAD`, `BCL2L1`,
`BIK`, `MCL1` and the composite. **It does not fall there for `BID` in SCAN-B,
for `PMAIP1`, or for `BCL2L11` in SCAN-B.** Those genes are MYC-ordered or flat
at group level, which matches E34 section 7. The heatmap is honest about this,
and that is one reason section 4.1 asks paragraph 3 to cite ED Fig. 10, not
4g, for the per-gene directions.

**Legend (descriptive only) must state:**

- the cohorts and n: TCGA 938 and SCAN-B 3,143, the model's complete cases;
- that the quadrants are median splits of MYC activity (`M_a`) and OXPHOS
  subunits (GSVA);
- the covariates for each cohort and rung:
  - TCGA `m1` is PAM50, purity and leukocyte fraction;
  - **SCAN-B `m1` is PAM50 only, because SCAN-B has no purity estimate**;
  - the open marker adds `PROLIF_DISJOINT`;
- **that the six composite genes and their signs were selected by their
  association with OXPHOS in these cohorts (Methods).** This describes how the
  strip was constructed. It is not an interpretation, and without it the
  reader cannot see that the OXPHOS separation is built in;
- the units, and that the bars are 95% CIs.

---

## 3. Fig. 4h and Fig. 4i: the decomposition

### 3.1 One axis convention for both

**Proposal: a ratio scale on a log axis in both panels**, with 1 marked:

- 4h shows the odds ratio for pCR per within-cohort SD;
- 4i shows the hazard ratio for breast-cancer death per within-cohort SD.

The axes are not shared, as v5 requires. The reasons for this choice:

- v5's "1 marked" only makes sense on a ratio axis;
- the text already quotes HRs for METABRIC;
- E44 already reports the pCR contrast as a ratio of odds ratios.

The alternative, log units with 0 marked in both panels, is equally valid.
**What is not valid is v5's mixture.**

**Add a direction cue under each axis.**

- 4h: "← lower pCR | higher pCR →".
- 4i: "← longer survival | shorter survival →".

**This is what makes the inversion legible.** On both panels PROLIF sits right of
1 and OX left of it. Without the cues, the two panels look like the same
pattern, although a higher OR is benefit and a higher HR is harm. An axis label
is description, so the legend rule allows it.

**Markers, the same in both panels.**

- **OX `m1`**: filled circle.
- **OX `m0`**: open circle, in 4h only.
- **PROLIF `m1`**: filled triangle.
- Use two Okabe-Ito colours for the two terms.
- Print no coefficient values on the panel. Print the contrast only.

### 3.2 Fig. 4h: pCR, three cohorts

**Source.** E40 (OX `m0`, `m1` and the meta-analysis), E43 (PROLIF), E44
(the contrast). No new fit is needed.

**Rows**, as in v5:

1. GSE194040 (988 patients, 319 pCR)
2. GSE164458 (482, 236)
3. GSE25066 (465, 93)
4. a rule, then **pooled**, in bold.

**Keep `m0` on the main panel** as an open marker, joined to OX `m1` by a thin
connector. That connector joins two rungs of one term, not two quadrants, so the
E34 rule does not apply. Do not use v5's ED fallback, for two reasons:

- the licence requires `m0` beside `m1` (`2026-10-04_handoff.md` 2.1);
- paragraph 4 cites 4h for "respiration alone bore no consistent relationship".

**Pooled summaries follow declaration section 9's rule:** random effects
(DerSimonian-Laird) at I2 ≥ 50, fixed effect otherwise. Print the I2 for each
term on the pooled row.

| pooled | log scale | ratio scale | summary used | I2 |
|---|---|---|---|---|
| OX `m0` | +0.0535 [-0.0498, 0.1567] | OR 1.05 [0.95, 1.17] | fixed | 48.8 |
| OX `m1` | -0.1616 [-0.2813, -0.0418] | OR 0.85 [0.76, 0.96] | fixed (equal to random) | 0.0 |
| PROLIF `m1` | +0.5422 [0.3067, 0.7778] | OR 1.72 [1.36, 2.18] | **random** | 66.9 |

**The contrast, `PROLIF − OX`.** Show it in a right-hand column on every row,
not only in the bracket on the pooled row:

| row | ratio of ORs [95% CI] |
|---|---|
| GSE194040 | 1.72 [1.28, 2.31] |
| GSE164458 | 2.07 [1.42, 3.04] |
| GSE25066 | 2.40 [1.44, 4.00] |
| **pooled** (fixed, I2 0) | **1.93 [1.56, 2.39]**, also on the bracket |

The per-row column is what puts paragraph 4's "the cohorts disagree about the
size of the proliferation coefficient and agree about the size of the gap" on
the panel.

**Legend:**

- the cohorts, n and events;
- that the scores are per within-cohort SD;
- the marker key;
- the contents of `m0` and `m1`;
- the meta-analysis rule;
- that the contrast is the within-model difference, estimated with the
  covariance taken from each fit's variance matrix;
- that the bars are 95% CIs.

**Avoid "control"** anywhere in the legend: PROLIF was declared as a covariate.

### 3.3 Fig. 4i: METABRIC breast-cancer-specific survival

**Source.** E43's and E44's METABRIC fits. No new fit is needed.

**Rows, changed from v5:**

| row | n / events | OX HR | PROLIF HR | ratio of HRs |
|---|---|---|---|---|
| ER-positive (IHC), **primary** | 1,505 / 431 | 0.896 [0.803, 0.9996] | 1.229 [1.074, 1.406] | **1.37 [1.13, 1.67]**, also on the bracket |
| ER-positive (`ER.Expr`), declared sensitivity | 1,512 / 432 | 0.929 [0.833, 1.036] | 1.181 [1.030, 1.355] | 1.27 [1.04, 1.56] |
| all patients (pooled ladder), declared comparator | 1,981 / 623 | 0.953 [0.872, 1.041] | 1.134 [1.009, 1.274] | 1.19 [1.005, 1.41] |
| *rule* | | | | |
| ER-negative (IHC) | 435 / 183 | 1.067 [0.909, 1.252] | 0.963 [0.765, 1.212] | 0.90 [0.65, 1.25] |

The ER-negative `ER.Expr` row goes in ED Fig. 8d.

**Why these rows.**

- **They are the three cells the narrative's sentence names.**
- **They solve the "equal weight overstates OX" problem
  (`2026-10-04_handoff.md` 2.1) structurally rather than by styling:**
  - OX's interval excludes 1 in the primary row only, by 3.7e-04;
  - it includes 1 in both comparators;
  - the contrast excludes 1 in all three.
- That is declaration 15.3 drawn as a figure: read the contrast, not the
  coefficient.
- De-emphasising the OX marker graphically would need a legend sentence to
  explain it, and the legend rule forbids one.

**Do not print the OX interval on the panel.** At three decimals the primary
row's upper bound rounds to 1.000.

**GSE25066 DRFS is dropped.**

- Declaration 15.4 excluded it from the contrast. The cohort has 103 events, 63
  of them in the ER-negative stratum that section 10 assigned to it.
- No PROLIF term for it exists on disk.
- v5's paragraph 4 never mentions it.
- **If the author wants it:** first a declaration section 16, then an E45
  refit on the **pre-specified 465** set, never `s508`, reported whatever it
  shows. Its OX-only ladder already belongs in ED Fig. 8b.

**Legend:**

- the cohort, endpoint, n and events per row;
- both ER definitions;
- the model contents, PAM50 spine plus `PROLIF_DISJOINT`;
- the contrast and covariance statement as in 4h;
- that the bars are 95% CIs.

### 3.4 Layout

- **Top row:** 4g at full width.
- **Bottom row:** 4h and 4i side by side, at about 3:2. 4i now has four rows and
  a contrast column, so v5's 2:1 is too narrow for it.
- Give one marker key for 4h and 4i.
- **If the budget forces two human panels, v5's fallback stands**: keep 4g, and
  merge 4h and 4i into one panel with two sub-plots.

---

## 4. Sentences in v5 that the drawn panels would contradict

The panels will be read against the text. Each item below is visible on the
panel as specified.

### 4.1 Paragraphs 2 and 3, against 4g

1. **"the two OXPHOS-high groups carrying the configuration equally"**
   - In SCAN-B at `m1`, Q2 is +0.23 and Q4 +0.37: `Q4 − Q2` is +0.134
     [0.037, 0.231].
   - At `m3` the gap widens in both cohorts.
   - E34 9.1 does not license "zero", and its section 4 says the manuscript
     "must say 'much smaller' and not 'none'".
   - "Equally" reads as "none".
2. **"adjusted for intrinsic subtype, tumour purity, leukocyte fraction and
   proliferation"**
   - SCAN-B has no purity or leukocyte estimate.
   - Proliferation is the `m3` companion, not the primary.
3. **"MYC contributing four to seven times less"**
   - Using E34 2.2's own ratio, neither rung gives that range in both cohorts:
     GSVA `m1` gives about 0 and 4.5-fold less, GSVA `m3` 7.4-fold and
     2.4-fold.
   - On mitoPPS the ratio is smaller still.
   - Quote the per-cohort ratio at the rung drawn.
4. **"the same amount at either respiratory level"**
   - This holds in 6 of 8 pooled cells.
   - **The exception is TCGA/GSVA, which is the cell 4g draws**: `m1` −0.002
     against +0.134, `m3` +0.115 against +0.278.
5. **Paragraph 3 cites Fig. 4g for ten per-gene directions.**
   - Those directions replicate in sign as proliferation-partial correlations
     (synthesis 3.1).
   - Several are small in SCAN-B: `BID` +0.091, `BCL2L11` −0.025, `PMAIP1`
     −0.064.
   - The quadrant heatmap shows them MYC-ordered or flat (section 2).
   - **Cite ED Fig. 10 for the per-gene directions.** The compartment sentences
     and "conditioning each axis on the other" (E11) should cite ED Fig. 7.

### 4.2 Paragraph 4, against 4h and 4i

6. **"correlating at Spearman 0.36 to 0.49" in all four cohorts**
   - The three neoadjuvant values are recorded (0.3575, 0.3812, 0.4934).
   - **METABRIC's is in no note or table.**
   - Computed read-only this session from `e41_metabric_ladder.rds$frame`
     (n 1,981), it is **0.381**, so the claim holds.
   - It must be recorded before it is cited. The build script can print it.

Not contradicted by the panels, but worth a second look:

- **"held in every subtype with enough events to fit".** TNBC is heterogeneous
  (I2 62.1, cohort signs disagreeing), not a held effect (handoff 2.3). "No
  subtype ran the other way" is what E42 licenses.
- **Units.** If the ratio axis is adopted, add the ORs beside the log-odds
  figures already quoted: PROLIF 1.72, contrast 1.93.

---

## 5. Decisions for the author

1. **4g row blocks.** By function (proposed), or another partition. MitoCarta
   versus cytosolic cannot be drawn from these twelve genes.
2. **4g instrument.** GSVA in the main panel and mitoPPS in ED Fig. 7
   (proposed, on cross-panel consistency), or both instruments in the strip.
3. **Axis convention for 4h and 4i.** Ratio on a log axis with 1 marked
   (proposed), or log units with 0 marked.
4. **4i rows.** The METABRIC primary, its two comparators and ER-negative IHC
   (proposed), or a declared E45 to add GSE25066 DRFS.
5. **Outside 4g-4i, but blocking the narrative. The validation study's place in
   the paper.**
   - CLAUDE.md says `myc_human_validation` "does not go in the paper".
   - v5 uses it three times:
     - paragraph 3's MYC/BCL2L1/MCL1 co-amplification (validation G2);
     - paragraph 4's "contrary to the prediction the pre-registered analysis had
       fixed";
     - ED Fig. 9, the H4 forest.
   - The validation repo's own write-up offers exactly that: one ED figure and
     one paragraph (`handoff-2026-10-21_validation_audit.md` section 8).
   - One of the two documents needs to change.
   - If ED Fig. 9 stays, the audit's contradictions 1 and 2 travel with it:
     - H4 is "two-instrument" in the draft but "GSVA alone" in the amendment;
     - the draft gives 1,978 patients but 1,940 were fitted.

---

## 6. Build, once the decisions are taken

- **Write one display script, `E45_fig4_human_panels.R`** (the next free
  number), under Option A: Claude Code writes it and the author sources it.
- **Inputs.** It reads only the E34, E40, E43 and E44 objects, and the E41
  frame for the METABRIC correlation.
- **Computation.** It estimates nothing. It only exponentiates and lays out.
- **No declaration is needed**, on the E15 precedent for display scripts, unless
  decision 4 adds the GSE25066 row. That row is a fit and needs one.
- **Assertions, in the manner of E44's component check:**
  - every plotted value equals its stored value;
  - 4i contains no `s508` cell and no GSE25066 row;
  - 4h uses the pre-specified 465 set;
  - the pooled summary type follows the I2 ≥ 50 rule.
- **Outputs:**
  - the panels, to `outputs/figures/` and to tracked `docs/figures/`;
  - the plotting frame, per R rule 4;
  - a digest table, written before the author's run.
