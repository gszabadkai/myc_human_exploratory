---
date: 2026-11-06
status: >
  PANEL INVENTORY and STATISTICS AUDIT for the legends still to be written:
  Fig. 2, Fig. 3, Fig. 4, Extended Data Fig. 4 and Extended Data Fig. 6.
  Every panel is listed with what it draws, its n, its error bars and its test,
  read from the author's assembled PDFs and from NM_all_figs.prism. Eleven
  flags, F1 to F11, need the author's decision; several change what a legend
  may say. The author confirms this before the wording is final.
relates-to:
  - docs/legends/README.md (the rules and the tools)
  - docs/legends/2026-11-06_legends_draft.md (the legends themselves)
  - myc_mouse 2ebf646 figures/natmetab/LEGENDS.md (the author's legend rules, read-only)
---

# Panel inventory and statistics audit

## 0. What was read

- **The assembled figures**, from the author's Drive folder: `NM_Fig2.pdf` and
  `NM_Extended_Data_Fig4.pdf` (2026-10-25), `NM_Fig3.pdf` (2026-10-29),
  `NM_Fig4.pdf` and `NM_Extended_Data_Fig6.pdf` (2026-11-06).
- **`NM_all_figs.prism`** (2026-11-06 14:58): 43 data sheets, 40 analyses,
  43 graphs. Its sheet titles carry the panel letters, so each panel is matched
  to its sheet by name.
- **The superseded legends**, in `11.8.26_MK_Myc_Paper_figures_tab.pdf`: the old
  Figures 3 to 5 and Supplementary Figures 3 to 5.

**Every *P* printed on a panel was matched to the Prism analysis that produced
it**, except where a flag says otherwise. The matches are exact to the printed
precision.

## 1. The error bars, measured

**These were not assumed.** Each panel's drawn error bars were measured from the
PDF's own vector geometry and converted into data units against that panel's
axis, then compared with the s.d., the s.e.m. and the 95% confidence interval
computed from the Prism values.

| panel measured | drawn half-length | s.d. | s.e.m. | verdict |
|---|---|---|---|---|
| Fig. 2b, Bcl-xL at 36 h | 4.92 | **4.95** | 2.86 | s.d. |
| Fig. 2b, EV at 36 h | 4.41 | **4.44** | 2.56 | s.d. |
| Fig. 2d, EV | 1.89 | **1.91** | 1.10 | s.d. |
| Fig. 2d, Pgc-1α | 8.53 | **8.59** | 4.96 | s.d. |
| Fig. 2e, Pgc-1α at 36 h | 4.67 | 2.92 | 1.46 | **95% CI (4.65)** |
| Fig. 2e, Pgc-1α at 64 h | 5.08 | 3.25 | 1.62 | **95% CI (5.16)** |
| ED Fig. 4a, all three bars | 0.0882, 0.436, 0.1494 | **0.0875, 0.435, 0.148** | 0.0505, 0.251, 0.0853 | s.d. |
| Fig. 4a, all three bars | 1.97, 0.944, 6.98 (×10⁸) | **1.97, 0.934, 6.96** | 0.569, 0.259, 2.10 | s.d. |
| Fig. 4b, four groups at day 26 | 59.1, 79.4, 89.6, 73.3 | 199, 263, 314, 244 | **59.9, 79.3, 90.6, 73.6** | s.e.m. |

**So the convention is:**
- **bar charts and the apoptosis time course: mean ± s.d.;**
- **tumour-growth and lung-colonisation curves: mean ± s.e.m.;**
- **Fig. 2e, the confluency curve: mean ± 95% confidence interval** (flag F9).

**Panels not individually measured** (Fig. 2g to 2j, Fig. 3, Fig. 4c to 4e,
ED Fig. 4b to 4g) are written to the convention of their own kind. The author
confirms, or names any that differ.

## 2. Fig. 2

**Title, as the Doc has it:** OXPHOS drives apoptosis in Myc expressing MECs via
NADH oxidation.

| panel | what it draws | *n* | centre, error | test, as run in Prism |
|---|---|---|---|---|
| **a** | immunoblot, iMMEC + *Tre-Myc* + EV or Bcl-xL, Dox 0 to 48 h: Myc, Bcl-xL, cParp, Atp5f1a, Uqcrc2, Mt-co1, Sdhb, Ndufb8, Actb | `[...]` | — | — |
| **b** | % apoptotic cells against Dox time, EV and Bcl-xL | 3 | mean ± s.d. | two-way ANOVA, Šídák; EV vs Bcl-xL within each time |
| **c** | immunoblot, MYAZ + EV or Pgc-1α: Pgc-1α, the five OXPHOS subunits, Bcl2l11 (BimEL, BimL), Bbc3 (PUMA), Bax, Actb | `[...]` | — | — |
| **d** | % apoptotic cells, MYAZ, EV and Pgc-1α | 3 | mean ± s.d. | one-way ANOVA, Dunnett against EV — **F1** |
| **e** | % confluency over 64 h, MYAZ, EV and Pgc-1α | 4 | mean ± 95% CI — **F9** | two-way repeated-measures ANOVA, two-stage FDR — **F3** |
| **f** | immunoblot, iMMEC + *Tre-Myc* + EV or Pgc-1α, Dox − or + | `[...]` | — | — |
| **g** | % apoptotic cells, iMMEC + *Tre-Myc*, EV and Pgc-1α × Dox | 3 | mean ± s.d. | two-way ANOVA, Tukey — **F1, F7** |
| **h** | the same, EV and Ndi1 × Dox | 4 | mean ± s.d. | two-way ANOVA, Tukey — **F7** |
| **i** | % apoptotic cells, MYAZ + *Tre-Myc*, EV, Ndi1, mt-lbNOX, c-lbNOX × Dox | 3 | mean ± s.d. | two-way ANOVA, Tukey — **F7** |
| **j** | % apoptotic cells against ABT-737, iMMEC and MYAZ | 2 | mean ± s.d. | **none in the file — F4** |

## 3. Extended Data Fig. 4

| panel | what it draws | *n* | centre, error | test |
|---|---|---|---|---|
| **a** | *Cox4i1* mRNA, fold change, MYAZ: EV, Pgc-1α, Ndi1 | 3 | mean ± s.d. | one-way ANOVA, Dunnett against EV |
| **b** | *Ndufs3* mRNA, the same | 3 | mean ± s.d. | one-way ANOVA, Dunnett against EV |
| **c** | OCR, basal, maximal and ATP-linked, EV and Pgc-1α | 3 | mean ± s.d. | two-way repeated-measures ANOVA, Benjamini–Hochberg FDR; the panel prints *q* — **F8** |
| **d** | *Bbc3* mRNA, fold change: EV, Pgc-1α, Ndi1 | 3 | mean ± s.d. | one-way ANOVA, Dunnett against EV — **F1** |
| **e** | three graphs: basal OCR, EV and Ndi1; OCR normalised to basal, ± rotenone, EV and Ndi1; ATP/ADP ratio, EV and Ndi1 | 3, 3, **2** | mean ± s.d. | one-way ANOVA with Dunnett; two-way ANOVA with Tukey for the rotenone pair — **F2, F5, F6** |
| **f** | ΔΨm, normalised TMRM fluorescence, EV and Ndi1 | **2** | mean ± s.d. | one-way ANOVA, Dunnett — **F2, F5, F6** |
| **g** | mitochondrial NAD(P)H/NAD(P)⁺, EV and Ndi1 | **2** | mean ± s.d. | one-way ANOVA, Dunnett — **F2, F5, F6** |

## 4. Fig. 3

**Title, as the Doc has it:** High Myc and OXPHOS levels during tumour
progression are associated with a disrupted apoptotic machinery.

| panel | what it draws | *n* | centre, error | test |
|---|---|---|---|---|
| **a** | % apoptotic cells, iMMEC + *Tre-Myc* + EV or Bcl-xL: GFP, Myc 48 h, Myc continuous | 3 | mean ± s.d. | two-way repeated-measures ANOVA, Dunnett within each vector |
| **b** | immunoblot of the same cells | `[...]` | — | — |
| **c** | Myc expression relative to GFP, the same groups | 3 | mean ± s.d. | two-way repeated-measures ANOVA, Tukey |
| **d** | immunoblot, iMMEC + *Tre-Myc* + EV or Pgc-1α, Dox − or + | `[...]` | — | — |
| **e** | % apoptotic cells, EV and Pgc-1α × Dox, after adaptation | 3 | mean ± s.d. | two-way ANOVA, Tukey — **F7** |
| **f** | three graphs: *p19arf*, *Cdkn1a* and *Myc* mRNA, EV and the Pgc-1α surviving fraction | 3 | mean ± s.d. | one-way ANOVA, Dunnett against EV — **F1** |
| **g** | tumour volume over time, EV and Pgc-1α; an immunoblot of those tumours; *Ppargc1a* mRNA (CPM) | 10–11 mice; `[...]`; 7 tumours | mean ± s.e.m.; —; median | mixed-effects model — **F3**; —; one-way ANOVA, Dunnett — **F1, F2** |
| **h** | tumour volume, EV + EV, Bcl-xL + EV, Bcl-xL + Pgc-1α; an immunoblot; PUMA protein; *Cdkn2a* mRNA (CPM) | 10 mice; `[...]`; 8; 7 | mean ± s.e.m.; —; mean ± s.d.; mean ± s.d. | two-way repeated-measures ANOVA — **F3**; —; one-way ANOVA, Dunnett; one-way ANOVA, Dunnett |
| **i** | OXPHOS priority (mitoPPS) against a matched control, and *Ppargc1a* (CPM), over the Pgc-1α, Bcl-xL and *p19arf*-knockout arms | 3 and 7 | **F10** | **F10** |
| **j** | tumour volume, NT + EV, NT + Pgc-1α, *p19arf* KO + EV, *p19arf* KO + Pgc-1α | 10–11 mice | mean ± s.e.m. | mixed-effects model — **F3** |

## 5. Fig. 4

**Title, as the Doc has it:** OXPHOS coupling to apoptotic regulation extends
from mouse to human breast tumours beyond its association with proliferation.

**Its letters are now settled by the assembled figure** (flag F11): the mouse
panels are **a** to **e** and the human panels **f** to **j**.

| panel | what it draws | *n* | centre, error | test |
|---|---|---|---|---|
| **a** | lung colonisation after intravenous injection, MYAZ: EV + EV, Bcl-xL + EV, Bcl-xL + Pgc-1α | 12, 13, 11 mice | mean ± s.d. | one-way ANOVA, Tukey |
| **b** | orthotopic tumour volume: EV EV, EV Ndi1, Bcl-xL EV, Bcl-xL Ndi1 | 11–12 mice | mean ± s.e.m. | mixed-effects model — **F3** |
| **c** | spontaneous lung metastasis, % area per lung, the same four groups, with representative lungs | 11, 10, 8, 10 mice | mean ± s.d. | one-way ANOVA, Tukey |
| **d** | lung colonisation over time, EV and Ndi1 | 9–11 and 6–9 mice, falling over time | mean ± s.e.m. | **none in the file — F4** |
| **e** | lung colonisation, relative burden: EV, mt-lbNOX, c-lbNOX | 10, 9, 10 mice | mean ± s.d. | one-way ANOVA, Tukey |
| **f** | fifteen BCL-2 family genes and the signature, standardised β on OXPHOS, MYC activity and proliferation, TCGA-BRCA and SCAN-B | 938 and 3,143 tumours | — | linear models, `E50` from `E49` |
| **g** | *BBC3*:*BCL2L1* against OXPHOS priority, both cohorts | 938 and 3,143 | — | `E52`; *P*, two-sided *t* test of the drawn slope |
| **h** | odds ratios for *BCL2L1* gain and *MCL1* amplification with *MYC* amplification | 1,043 tumours | 95% CI | `E54`, from the validation study's saved result |
| **i** | odds ratios for pathological complete response, three neoadjuvant cohorts and the pooled row | 988, 482, 465 | 95% CI | `E57` |
| **j** | hazard ratios for breast-cancer death, METABRIC | 1,505, 1,512, 1,981, 435 | 95% CI | `E57` |

**The human panels' numbers are already verified in this repo**, each against
digests committed before its run: `E50` (4f), `E52` (4g), `E54` (4h) and `E57`
(4i, 4j). Their legend drafts are in those scripts' notes and are condensed for
the legend.

## 6. Extended Data Fig. 6

| panel | what it draws | drawn by |
|---|---|---|
| **a** | correlations among the scores, both cohorts, with the neoadjuvant cohorts beneath | `E53` |
| **b** | the OXPHOS association of *BCL2L1* and *BBC3*, and their difference | `E53` |
| **c** | unadjusted Spearman correlations of the fifteen genes and the signature with all seven measures | `E55` |
| **d** | the fifteen genes across the four MYC × OXPHOS groups, both splits, with the signature beneath | `E56` |

## 7. The flags

**Each needs the author's decision. F1 to F7 change what a legend may say.**

### F1. Several tests include groups that are not drawn

The adjusted *P* on these panels comes from a family that contains groups the
panel does not show. The printed values are the correct members of that family,
but the correction is over more comparisons than the panel displays.

| panel | drawn | also in the test |
|---|---|---|
| Fig. 2d | EV, Pgc-1α | Nrf1 |
| Fig. 2g | EV ± Dox, Pgc-1α ± Dox | Nrf1 ± Dox |
| Fig. 3f, all three | EV, Pgc-1α | Ndi1, Nrf1 |
| Fig. 3g, *Ppargc1a* | EV, Pgc-1α | Ndi1, Nrf1 |
| ED Fig. 4d | EV, Pgc-1α, Ndi1 | Nrf1 |

**The decision:** either the legend says the test also included the other arms,
or the test is re-run on the drawn groups and the *P* values are updated. **The
draft names the test as it was run and does not mention the undrawn arms**,
which is the smaller claim but leaves the family unstated.

### F2. Duplicated columns inflate the comparison family

| sheet | what is duplicated |
|---|---|
| ED Fig. 4e, basal OCR | Ndi1 appears twice, identical |
| ED Fig. 4e, ATP/ADP | Ndi1 twice, identical |
| ED Fig. 4e, rotenone | three Basal/+Rot pairs, two of them identical |
| ED Fig. 4f, TMRM | Ndi1 twice, identical |
| ED Fig. 4g, redox ratio | Ndi1 twice, identical |
| Fig. 3g, *Ppargc1a* CPM | the Ndi1 column equals the Pgc-1α column |

Each duplicate is counted as a separate comparison. The effect on the printed
value is small, but the stated family is wrong, and for ED Fig. 4e to 4g the
drawn comparison is a single pair, for which an unpaired *t* test is the natural
test rather than a one-way ANOVA with Dunnett over a duplicate.

### F3. Five analyses have no stored results

`isCalculated` is false and the results file is empty for: Fig. 2e's multiple
*t* tests, and the growth models of Fig. 3g, Fig. 3h, Fig. 3j and Fig. 4b. **So
the *P* values printed on those panels could not be checked against the file.**
Their test and design are readable from the saved parameters and are used in the
draft. Opening the file in Prism so it recalculates would let them be checked.

### F4. Two panels print *P* values with no analysis in the file

**Fig. 2j** prints four and **Fig. 4d** two. Neither has an analysis sheet. The
draft leaves the test as `[F4]`.

### F5. *n* = 2 in three Extended Data panels

ED Fig. 4e (ATP/ADP), 4f (TMRM) and 4g (redox ratio) have two values per group.
**The superseded Supplementary Fig. 5 legend gave *n* = 3 for TMRM and the redox
ratio.** An s.d. and a Dunnett test on two points are both drawn from them.

### F6. Three of those panels normalise the control to a constant

In ED Fig. 4e (rotenone, and ATP/ADP) and 4f (TMRM) every control value is
exactly 100 or 1, so the control's s.d. is zero and the test compares a constant
with a variable. The draft says "normalised to EV", which is what is drawn.

### F7. Column data analysed as a two-way ANOVA

Fig. 2g, 2h, 2i and Fig. 3e are single column families analysed as two-way
ANOVAs, so Prism's row factor is the replicate index. **Its *P* is meaningless**
(0.29, 0.65, 0.66, 0.08) and is not printed on the panels. The drawn comparisons
are the column Tukey ones and are right. The draft names "two-way ANOVA with
Tukey's multiple comparisons test", which is what was run.

### F8. The unit of *n*, and ED Fig. 4c's *q*

Methods says *n* ≥ 3 **biological replicates**. The draft uses "biological
replicates" for the cell work and "mice" for the *in vivo* work. **Two to
settle:**
- Fig. 2e's *n* = 4 — wells of one experiment, or four experiments?
- ED Fig. 4c is corrected by false discovery rate, so its printed values are
  *q*, not *P*. The panel prints them as "P". The draft calls them *q* and names
  the method.

### F9. Fig. 2e's error bars are 95% confidence intervals

Measured, twice, against the Prism values. Every other cell panel is mean ± s.d.
**Confirm this is intended**, since it is the only panel of its kind.

### F10. Fig. 3i is not in the Prism file

Its two parts are the retained-dose series from the mouse repo's script 52. Its
intervals are a model's, not Prism's. **The author supplies what they are**, or
points at the saved object and I read it from `myc_mouse` read-only.

### F11. The letters are settled, and three citations follow from them

The assembled Fig. 4 makes the human panels **f, g, h, i, j**. This answers the
question the text ruling left open (the handoff's section 0, items 1, 2 and 6):
- **R6** cites the PUMA relative to BCL-xL panel: now **Fig. 4g**;
- **the co-amplification sentence** cites **Fig. 4h**;
- **the outcome forests** are **Fig. 4i** and **Fig. 4j**;
- **R13's citation is still open.** Its claim is the joint model, whose
  coefficients **Fig. 4f** draws; ED Fig. 6c keeps only the unadjusted side.

**Extended Data Fig. 6 is a to d**, so R1 and R15 cite **6a**, and R5 cites
**6b**.

## 8. What the old legends gave that the new ones must not lose

The superseded legends are the author's own, in `11.8.26_MK_Myc_Paper_figures_tab.pdf`.
Everything in them that is a fact about a panel is carried into the draft except
these, which need the author:

- **"500 ng/mL of dox"** — in Methods, so the draft says "Dox" and leaves the
  concentration to Methods;
- **"14 days"** for the continuously induced survivors, and **"48 hrs"** for the
  short induction — both carried;
- **"n=3, representative western blot shown"** — carried as `[...]`, because the
  blots' *n* is the author's to confirm;
- **"randomly chosen tumours per group"** for the *Ppargc1a* measurement —
  carried;
- **rotenone 1 µM**, **TMRM 25 nM** — in Methods;
- **the old Supplementary Fig. 5I**, cyto- against mito-LbNOX on metastasis, is
  now **Fig. 4e**;
- **the old Figure 4A–C**, Pgc-1α alone *in vivo*, is now **Fig. 3g**.
