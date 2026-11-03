---
date: 2026-11-02
amended: 2026-11-03, the author's amendments and a section 8 addendum (see approval)
status: >
  DECLARATION, pre-read, APPROVED WITH AMENDMENTS by the author. No E48
  quantity has been computed, printed or plotted: no per-gene proliferation
  coefficient, no within-model contrast, no Part C fit. Section 2 lists
  everything that WAS seen, corrected from the record. Section 2.9 records
  what those numbers already imply, and section 11 states the expectation
  from it.
repo: myc_human_exploratory. myc_human_validation (d3ac60e) and myc_mouse are read-only.
governed-by: docs/2026-11-02_separation_line_brief.md (cbf2e4c), the author's decision record. Its section 4 governs every word written about E48's results.
relates-to:
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (dcb58a7; addendum 3b5e221); the format model
  - docs/2026-10-23_e47_result.md (E47: MIXED on S6, which stands)
  - docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md sections 14 and 15 (E43 and E44: the precedents for Part A)
  - docs/2026-10-22_fig4_human_panels_proposal.md sections 7 and 8 (E45's and E46's specifications, which 4g reuses)
  - docs/2026-11-02_human_section_reconciliation_proposal.md (its AD HOC numbers are listed in 2.3)
decides:
  - The readouts, the quantities, the reading rule and the outcome-to-text table, Parts B and C, the displays (4g in three layouts; ED-a to ED-d), and the saved objects and assertions, all before any E48 quantity is read.
approval: >
  2026-11-03, by the author, with nine amendments:
  1. Section 11 restated from what the seen coefficients imply: the identity,
     the implied sign per cell, the expected labels (MCL1 left open between
     DISTINCT and NOT SEPARATED), the agreement check, a Methods sentence.
  2. Section 8 row 3 accepted. Any gene not OPPOSITE gets one descriptive
     sentence; the separation and closing sentences name the same genes.
  3. Section 7.4's tally accepted; its sentence is optional.
  4. Part B unavailable, no refit. E39's gate values join ED-a. The optional
     pCR MYC sentence leaves the text.
  5. C2 dropped, with its reason (10.2).
  6. C4 added (10.4).
  7. C1 kept; TP53 mRNA's membership of the proliferation score is noted.
  8. C3 and its replicate gate kept.
  9. The corrected seen-before list accepted. STYLE.md will be in the repo
     before stage 3 and then governs prose.
  10. At the OK: section 8's closing paragraph, nothing is discarded. To
      match it, section 5 extracts the twelve genes on both MYC scores, ED-d
      shows both, and S-1 displays C1 and C4.
  11. After 1ccf805 and before the author's run: section 8's addendum. Each
      read gene is named with its own label, and "the BCL-2 family" is never
      used for them.
next-action: >
  The author sources scripts/E48_separation_and_fig4g.R (1ccf805). The run is
  verified against docs/2026-11-03_e48_prerun_digests.md, then stage 3.
---

# E48: OXPHOS, proliferation and MYC activity on BBC3, BCL2L1 and MCL1, each with the others held fixed; and Fig. 4g

## Terms used from here on

| term | meaning |
|---|---|
| **OXPHOS level** | the `OXPHOS subunits` score on GSVA (VST input). "Level" in panels and text |
| **OXPHOS priority** | the same arm on mitoPPS (linear input). "Priority" in panels and text |
| **MYC activity score** (repo: **estimator**) | a score of MYC activity. **Never a "signature".** Four are fitted in E47: `M_a` = `FELSHER__MITOSTRIP` (61 genes, the primary); the **least-entangled score** = `MYC_UP.V1_UP__MITOSTRIP` (178); HALLMARK = `HALLMARK_MYC_TARGETS_V1__MITOSTRIP` (177); the CollecTRI regulon = `M_b__MITOSTRIP` |
| **proliferation** | `PROLIF_DISJOINT` (318 genes), the term `PROLIF` in E47's fits |
| **m1, m3** | E47's rungs: `m1` = OXPHOS + MYC + spine; `m3` = `m1` + proliferation |
| **the difference** | proliferation's coefficient minus OXPHOS's, **within one `m3` fit**, E44's orientation |
| **block** | one cohort-by-OXPHOS-measure combination; there are four |
| **separable** | in E39's sense only: below 0.60 fine; 0.60-0.80 not cleanly separable; 0.80 or above not identifiable |

## 0. Why this exists

**The author's line (brief section 0):** in human breast tumours OXPHOS, MYC
activity and proliferation rise together; modelled together, OXPHOS and
proliferation run in opposite directions, on the BCL-2 family and on treatment
response; and the apparent association of MYC activity with the BCL-2 family
depends on how much proliferation its score carries.

**The quantity that line needs has never been declared.** Proliferation's own
coefficient sits in every one of E47's 272 `m3` fits and has been printed for
the signature only (E47 result 3.1, RETRIEVED). **Per gene it has never been
looked at**, and no within-model contrast has ever been computed on an E47 fit.
On the clinical side the same two quantities exist and are declared: E43
retrieved proliferation's coefficient and E44 tested the difference.

**E48 puts the BCL-2-family half on the same footing** (Part A), records what
the MYC-on-response check can and cannot do (Part B), declares three
sensitivities (Part C), and builds Fig. 4g so that it shows the line on its own.

## 1. Posture

**EXPLORATORY and POST HOC.** It is written after E47's result and in response
to it. Nothing here is pre-registered. It fixes its rule before its own
quantities are read, which is all a post-hoc analysis can do.

- **Part A fits nothing.** It recovers a declared covariate's coefficient and
  within-model contrasts from E47's verified saved fits (coefficients and
  `vcov`), on the E43 and E44 precedents.
- **Part B is unavailable** (section 9) and fits nothing.
- **Part C is three sensitivities, C1, C3 and C4**, declared here,
  **reported and never read** (section 10). Each involves new fits. **C2 is
  declared and not run** (10.2).
- **E47 stands.** It is not refitted, its MIXED verdict is not amended, and E48
  does not read S6.
- **Proliferation's coefficient is "a post-hoc characterisation of a declared
  covariate"** (brief 4.7), as E43's was.
- **Main effects only. No product term is fitted.**
- **`E48` is confirmed as the next free script number**: no `E48` script, no
  `E48` branch locally or on origin (checked 2026-11-02).

## 2. What was seen before this was written

**The brief's 5.2 list, corrected from the record.** The record searched: this
repo's notes; E47's script, whose console output prints only the OXPHOS and
MYC rows; and the four session transcripts on disk, which cover 2026-10-08
onward. **For earlier sessions the notes are the record.**

### 2.1 E47's result note in full

**This includes proliferation's coefficient on S6 at `m3`, which was printed
once, for 8 cells** (`M_a` and the least-entangled score, both cohorts, both
measures): -0.221 to -0.379, every CI excluding zero (E47 result 3.1).
**No other readout's proliferation coefficient was printed**; the extraction
filtered `readout == "S6"`.

### 2.2 E47's OXPHOS and MYC coefficients: more than the brief lists

- **Every readout on `M_a` and on the least-entangled score, at `m1` and
  `m3`, both coefficients.** That is all twelve genes, S6 and both arms. They
  were printed on 2026-11-02 for the reconciliation proposal.
- **E47 Displays A, B and C**, viewed. Display B shows MYC's coefficient for
  every readout on all four scores, at `m1` and `m3`.
- **Counts over the 44 non-void specifications**, both coefficients, for
  `BAD`, `BBC3`, `BCL2`, `BCL2L1`, `BID`, `MCL1` and S6.
- **MYC's coefficient on `BBC3` and `BCL2L1`**, broken down by all four scores
  and all three rungs, including `m0`.

### 2.3 The reconciliation proposal's AD HOC correlations, all computed on `e47$frame`

- **Raw Spearman:** the full block of OXPHOS level, OXPHOS priority, `M_a`,
  the least-entangled score and proliferation, against those five and S6,
  `BBC3` and `BCL2L1`. It includes:
  - proliferation with `BBC3` -0.30 / -0.15, with `BCL2L1` -0.13 / -0.03, and
    with S6 -0.09 / -0.04;
  - level with priority 0.88 / 0.82;
  - `M_a` with the least-entangled score 0.84 / 0.79.
- **Partial on proliferation:**
  - `M_a` with `BBC3` 0.174 / 0.105, and with `BCL2L1` -0.028 / 0.029;
  - level with `BBC3` 0.370 / 0.170, and with `BCL2L1` 0.427 / 0.348.
- **`MCL1` was in none of them.**

### 2.4 E47's gates and diagnostics

- **`e47$gates$g3`, all 16 rows**: OXPHOS-MYC, MYC-proliferation and
  OXPHOS-proliferation, for every score, cohort and measure.
- **D1 and D2.**

### 2.5 E46

- **The figure**: twelve genes by quadrant, both measures, and the signature
  strip at `m1` and `m3`.
- **The numeric quadrant means for `BBC3` and `BCL2L1`** in all four blocks
  (reconciliation proposal 3.1).

### 2.6 E40, E43, E44 and E45

Their result notes, and E45's figure. **This includes E43's per-cohort
proliferation coefficients on pCR and on METABRIC.**

### 2.7 Earlier work

- **E22, on `BBC3`.** Adding `PROLIF_DISJOINT` moved MYC's coefficient by
  +0.332 and OXPHOS's by -0.015, on E22's own rulers and analysis set. Its
  twelve genes' MYC coefficients at `+PROLIF` are listed in E47 declaration
  2(b).
- **The provenance table (E47 declaration 4a.2).** It gives each gene's
  correlation with OXPHOS and with MYC, partial on proliferation.
- **E10, E11 and E15 record that `TP53` and `BIRC5` are members of
  `PROLIF_DISJOINT`.** Verified for both proliferation sets on 2026-11-02.

### 2.8 The structure checks for this declaration (section 6)

These read names, dimensions, term labels and covariate counts only. **No
coefficient and no contrast was read.**

### Confirmed unseen

- **Any per-gene proliferation coefficient in any E47 fit.**
- **Any within-model contrast on any E47 fit.**
- **The raw correlation of `MCL1` with proliferation.**
- **E40's `m2` MYC coefficient.** It exists nowhere on disk (section 9).
- **Anything in this repo stratified or adjusted by TP53 status.** The atlas
  strata were all, ER+, ER-, LumA, LumB, HER2, Basal and Normal. TP53 status
  was a covariate in the frozen validation study's Block C models, whose
  coefficients are not read here.

### 2.9 What the seen numbers already imply: arithmetic, not a fit

**The identity.** E47's `m1` and `m3` fits are nested on the same tumours, so
by least-squares algebra:

```
beta_MYC(m3) - beta_MYC(m1) = - beta_PROLIF(m3) * delta
```

- **`delta` is MYC's coefficient when proliferation is regressed on `m1`'s
  terms** (OXPHOS, MYC and the spine).
- **In words, adding proliferation moves MYC's coefficient by minus
  proliferation's coefficient times `delta`.**
- **`delta > 0`:** MYC activity predicts proliferation positively (rho
  0.615-0.696 for `M_a`, 0.413-0.522 for the least-entangled score).
- **So the direction of MYC's `m1`-to-`m3` shift fixes the sign of
  proliferation's coefficient:** a rise means negative, a fall positive.
- **`delta` was not computed, so no size is implied.** E48 computes it and
  checks the identity in every cell (G-A4).

The shifts, re-read from `e47$results`, and the implied sign of
proliferation's `m3` coefficient. Blocks are in the order TCGA level, TCGA
priority, SCAN-B level, SCAN-B priority.

| gene | MYC score | MYC's shift, `m1` to `m3` | implied sign of proliferation's coefficient |
|---|---|---|---|
| `BBC3` | `M_a` | +0.182 / +0.191 / +0.108 / +0.114 | neg / neg / neg / neg |
| `BBC3` | least-entangled | +0.097 / +0.094 / +0.036 / +0.041 | neg / neg / neg / neg |
| `BCL2L1` | `M_a` | +0.051 / +0.071 / +0.091 / +0.101 | neg / neg / neg / neg |
| `BCL2L1` | least-entangled | +0.043 / +0.045 / +0.041 / +0.045 | neg / neg / neg / neg |
| `MCL1` | `M_a` | **+0.069 / +0.053** / -0.039 / -0.048 | **neg / neg** / pos / pos |
| `MCL1` | least-entangled | **+0.005 / +0.003** / -0.021 / -0.024 | **neg / neg** / pos / pos. TCGA's shifts are near zero |

**What follows is section 11.** In brief:
- **the signs for `BBC3` and `BCL2L1` were foreseeable**, and this declaration
  cannot claim otherwise;
- **`MCL1` is not expected OPPOSITE;**
- **the agreement check is expected to hold in all 12 cells.**

**A sign difference between the cohorts for `MCL1` is a pattern the brief's
composition alternative names** (4.4). TCGA adjusts for purity and leukocyte
fraction; SCAN-B cannot. **C4 (10.4) tests whether it follows the
adjustment, and the result note must weigh it whatever E48 shows.**

## 3. Readouts, and why these three

**READ: `BBC3`, `BCL2L1` and `MCL1`.** Why these three, fixed before reading:
- **`BBC3` and `BCL2L1` encode PUMA and BCL-xL**, the proteins examined in the
  mouse (Figs 2 and 3);
- **`BCL2L1` and `MCL1` are the primary BCL-2-family measures** fixed by the
  27 August validation plan, in its section 9;
- **none of them is chosen by this analysis.**

**REPORTED, not read: the other nine genes** (`BAD`, `BIK`, `BID`,
`BCL2L11`, `BMF`, `PMAIP1`, `BCL2`, `BCL2L2`, `BCL2A1`). **No per-gene FDR.**

**NOT USED: S6, S5, S7 and both arms.** Their signs follow OXPHOS by
construction.
- E48 reads no fit of theirs.
- They are drawn in no E48 panel.
- **This is asserted** (section 14).

**G2 carries from E47.** None of the twelve sits in any MYC score, in
`OXPHOS subunits` or in `PROLIF_DISJOINT`.

## 4. Model and inputs

**E47's fits, unchanged**, read from `results/e47_mutual_adjustment.rds`
(verified 22 of 22 against digests committed at `9251355`):

- **The model:** `readout ~ OX + MYC + spine (+ PROLIF at m3)`, OLS, with every
  variable Blom-scored per cohort.
- **The spine:**
  - TCGA: PAM50 + purity + leukocyte fraction;
  - SCAN-B: PAM50, with no purity estimate (trap 2).
- **The analysis sets are E34's**: TCGA 938, SCAN-B 3,143.
- **OXPHOS:** level and priority.
- **The MYC scores:**
  - `M_a` is the primary;
  - the least-entangled score is used for the agreement check (7.2);
  - HALLMARK and the CollecTRI regulon are reported only. **HALLMARK's `m3`
    cells are void** under E47's G3 (68 fits).
- **The snapshot** is `data/from_validation/`, via E34. The joint snapshot is
  not read.

**Part A fits nothing.**

## 5. Quantities

**Every quantity is computed only for the gene fits.** Proliferation's
coefficient and the difference are extracted for:
- **the twelve genes on `M_a` and on the least-entangled score** at `m3`, in
  all four blocks. The twelve on the least-entangled score follow from section
  8's closing paragraph;
- `BBC3`, `BCL2L1` and `MCL1` on the CollecTRI regulon at `m3`, in all four
  blocks.

HALLMARK's `m3` cells are void, so they are extracted and flagged void, never
read.

1. **Proliferation's coefficient at `m3`** (term `PROLIF`). It has the
   proliferation, OXPHOS and MYC terms and the spine all in the model.
   - It is read from E47's coefficient matrix: estimate, SE, and t-based 95%
     CI on the fit's residual df.
   - It is recomputed from `coef` and `vcov` as a gate (6.2).
2. **The within-model difference, proliferation minus OXPHOS, at `m3`.**
   - It is a general linear contrast: `L beta`, with `L` = +1 on `PROLIF` and
     -1 on `OX`.
   - `Var = L V L'` = `Var(PROLIF) + Var(OX) - 2 Cov(PROLIF, OX)`, from the
     fit's saved `vcov`.
   - The 95% CI is t-based, on the fit's residual df.
   - Reported beside it, as E44 15.5 did: the correlation between the two
     coefficient estimates, and the naive SE. **The naive SE is for
     comparison only and is never the statistic.**
   - **Never the difference of two separately summarised estimates.**
3. **Re-read, not new: the OXPHOS and MYC coefficients at `m1` and `m3`.**
   They are asserted equal to `e47$results`.
4. **Descriptive: Spearman correlations among both OXPHOS measures, the four
   MYC scores and proliferation**, per cohort, on the analysis sets.
   - **`g3` holds** OXPHOS-MYC (per measure and score), MYC-proliferation and
     OXPHOS-proliferation.
   - **`g3` does not hold** level-priority or MYC-MYC. These are computed on
     `e47$frame`.
   - **Gate:** the `g3` pairs, recomputed the same way, must reproduce `g3`
     exactly.
5. **E47's VIF for every read cell**, from `e47$vif`, as `GVIF^(1/(2*Df))`
   for multi-df terms.

## 6. Gates and structure checks

### 6.1 Structure, checked 2026-11-02 (names, dimensions and term labels only)

**`results/e47_mutual_adjustment.rds`:**

| element | what it holds |
|---|---|
| `coefficient_matrix` | 5,440 rows for 816 fits (272 per rung); columns `fit_id`, cohort, instrument, estimator, rung, readout, term, estimate, se, t, p, `ci_lo`, `ci_hi`, n, **df** |
| `vcov` | a list of 816 square matrices. **For all 544 `m1` and `m3` fits, the dimnames equal the fit's coefficient terms** |

**Term sets:**

| fit | terms |
|---|---|
| TCGA `m3` | `(Intercept)`, `OX`, `MYC`, PAM50 (Basal, HER2, LumB, Normal; LumA is the reference), `purity`, `leuko`, `PROLIF` |
| SCAN-B `m3` | the same without `purity` and `leuko` |
| `m1` | the same without `PROLIF` |

**`PROLIF` is in all 272 `m3` fits**, in both the coefficients and the `vcov`,
**and in none of the `m1` fits.**

**Void cells:** 68, all HALLMARK `m3`.

**`PROLIF` is `PROLIF_DISJOINT`** (E47 `PROLIF_SET`). It shares 0 genes with
`M_a`, 2 with the least-entangled score and 41 with HALLMARK.

**So Part A is computable without a refit.**

### 6.2 Gates in the script, before any quantity is read

If G-A1, G-A2 or G-A4 fails, nothing is read.

- **G-A1, identity.** The re-read OXPHOS and MYC coefficients equal
  `e47$results` exactly, for every gene fit read.
- **G-A2, recomputation.** For every cell extracted:
  - each stored `PROLIF` and `OX` estimate equals the fit's `coef`;
  - each SE equals `sqrt(diag(vcov))`;
  - each stored CI equals one recomputed with `qt(0.975, df)`.

  The tolerance is 1e-10.
- **G-A3, void.** Void cells come from `e47$results$void`. They are never
  read, never counted, and never drawn as read.
- **G-A4, nesting.** For every gene fit read, MYC's and OXPHOS's shifts from
  `m1` to `m3` must each equal minus proliferation's `m3` coefficient times
  that term's coefficient in the regression of proliferation on `m1`'s terms,
  in the same tumours, to 1e-10.
  - **A failure means the two rungs are not nested on one sample.**
  - **`delta`'s sign is reported per cell**, beside section 11, which assumes
    it is positive.
- **G3 is carried.** `M_a`'s correlation with proliferation is 0.615 to
  0.696, "not cleanly separable".
  - It is printed beside every read cell, with E47's VIF.
  - **The Methods says so.**

## 7. Reading rule, FIXED BEFORE READING

### 7.1 Primary cells and labels

**The primary cells:** 3 genes x 2 cohorts x 2 OXPHOS measures = **12 cells**,
on `M_a` at `m3`. **Each gene is labelled over its 4 cells.**

| label | condition |
|---|---|
| **OPPOSITE** | **In all 4 cells**: proliferation's CI excludes zero with the sign opposite to OXPHOS's `m3` coefficient in that cell, **and** the difference's CI excludes zero |
| **DISTINCT** | not OPPOSITE, but **the difference's CI excludes zero in all 4 cells, with one sign** |
| **NOT SEPARATED** | anything else |

- **"Opposite to OXPHOS's" uses OXPHOS's own `m3` coefficient in that cell.**
  If its CI includes zero there, "opposite" is undefined and the cell fails
  OPPOSITE. In all 44 specifications seen, OXPHOS's CI excludes zero for these
  three genes.
- **The nine other genes get no label.**

### 7.2 The agreement check

- **The test.** On the least-entangled score at `m3`, proliferation's
  coefficient must carry the same sign as on `M_a` in all 12 cells. **Sign
  only.**
- **If any cell disagrees,** ESTIMATOR-DEPENDENT is printed beside that gene's
  label. **No manuscript sentence may then state the label without naming the
  dependence.**

### 7.3 Order of reading

The label is assigned first. The agreement check is then applied to it.
**Neither CollecTRI nor HALLMARK enters the rule.**

### 7.4 The nine other genes: a descriptive tally, fixed now

The brief's draft says that OXPHOS and proliferation "ran in opposite
directions for n of them". **So that n is not chosen after reading, it is
defined now** (accepted by the author, 2026-11-03):
- **A gene counts toward n** if, in all 4 `M_a` `m3` cells, both
  coefficients' CIs exclude zero with opposite signs.
- **No contrast enters the tally, and no label.**
- **The tally is descriptive only.**
- **The sentence is optional.** ED-d may instead be cited from the OXPHOS
  sentence.

## 8. What each outcome licenses, fixed now

**Rows 1 to 3 are exclusive on the labels. Row 4 applies in addition to
whichever of them holds.**

| result | text |
|---|---|
| **1. All three OPPOSITE**, agreement holds | Abstract: **"opposite associations"**. Brief section 8's sentence "Proliferation ran the other way ..." **as drafted** |
| **2. None NOT SEPARATED, at least one DISTINCT** | **"distinct associations"**. Name the genes on which proliferation ran the other way |
| **3. Any gene NOT SEPARATED** | **That gene leaves the separation sentence.** If it is `BCL2L1`, use the reconciliation proposal's section 4 wording for `BCL2L1`. **The clinical paragraph is unaffected.** **The abstract names the genes on which proliferation ran the other way, and does not say "the BCL-2 family"** (accepted by the author, 2026-11-03) |
| **4. Any ESTIMATOR-DEPENDENT** | State the label with the dependence named, or do not state it |

**Two further rules, whatever the outcome** (the author, 2026-11-03):
- **Any gene that is not OPPOSITE gets one descriptive sentence**, giving the
  direction of its proliferation association in each cohort, **with no
  inference.**
- **The separation sentence and the closing sentence name the same genes.**

**No outcome licenses:**
- a causal, mediation, interaction or "primed" claim;
- **any link from the BCL-2 family to outcome** (brief 4.3). The two
  separations run in parallel; the data do not chain them.

**Part C never changes a label.** **Brief section 4 governs every word**,
including its item 9's list of words that stay out.

**Nothing is discarded** (the author, 2026-11-03). Nothing computed under this
declaration is discarded. Every coefficient and difference, for all twelve
genes, both MYC scores and C1, C3 and C4, is saved, tabulated and displayed,
and may be described in the text whatever its label. The labels decide only
the words used for the three read genes. Patterns this declaration did not
anticipate are reported in an exploratory section of the result note, stated
as hypotheses, and interpreted there against the mouse data.

**ADDENDUM, 2026-11-03 (the author), after `1ccf805` and before the author's
run: each read gene is named, with its own label.**

> In every row, the abstract and the separation and closing sentences name the
> read genes, each described by its own label: 'opposite' for a gene labelled
> OPPOSITE, 'distinct' for one labelled DISTINCT, and one descriptive sentence,
> with no inference, for a gene NOT SEPARATED. The phrase 'the BCL-2 family' is
> not used for the three read genes in any row.

- **When it was written.** After the script and its pre-run digests were
  committed (`1ccf805`), and before the author's run.
  - The author had seen no E48 output.
  - Claude Code had looked at the dry-run panels for legibility only. They
    show the plotted values. No result was reported.
- **What it changes.** The wording of rows 1 to 3, and nothing else.
  - **It supersedes "the BCL-2 family" for the three read genes** wherever
    brief section 8's drafts use it: the abstract and the closing paragraph.
  - **It changes no label, no part of section 7's rule, no computed quantity
    and no line of the script.**
- **The script's saved `spec$outcome_table` predates this addendum.** It is
  left unchanged, so the committed digests still hold. **Where the two
  differ, this declaration governs.**

**CORRECTION, 2026-11-03, at the author's request: the commit followed the
author's run.** "Before the author's run" is true of the author's text, not of
its commit. All times are UTC:

| time | event |
|---|---|
| 01:13:26 | the author's clarification arrived |
| about 01:13:35 | the author's run began. The dry runs took about 19 s from start to first figure, so this is after the clarification |
| 01:13:54 to 01:14:04 | the run wrote its outputs |
| 01:14:19 | this addendum was committed (`258a496`) |

**So git cannot show that the addendum came first.** This is the same lapse as
E45's digests (proposal 7.6). **It changes nothing the run applied**: the
addendum concerns wording only, and the script was not changed.

## 9. Part B: MYC activity on pCR. UNAVAILABLE

**The structure check (brief 5.9.1), 2026-11-02:**

| object | what it holds |
|---|---|
| **E40** (`results/e40_ladder.rds`) | `$pcr` and `$meta` hold **the OXPHOS row only** at every rung. The `m2` rows carry OXPHOS's estimate, and their `terms` string records MYC as fitted. There are **no fitted model objects, no coefficient matrix and no `vcov`.** `$discordant` is a 54-row subset, not the frame |
| **E43** (`results/e43_proliferation_coefficient.rds`) | full coefficient matrices, VIFs and the fitted objects, **for `m1` only**: eight fits of `OX + subtype_term (+ treatment) + PROLIF`, **with no MYC term** |
| **E41** | no MYC term |

**So, under brief 5.9.3, Part B is recorded as "unavailable".** E48 fits
nothing for it. **A refit needs its own declaration.**

**Recorded now for that declaration:**
- **E40's MYC score is `M-a`, the stripped Felsher estimator** ("M-a is the
  MYC estimator throughout here", E39 declaration section 4). It comes from
  the H4 frames built by the validation study's script 13.
- **E39's gate gives `rho(MYC, PROLIF_DISJOINT)`** as:
  - **0.7797** in GSE25066 (n 508);
  - **0.7851** in GSE194040 (988);
  - **0.8236** in GSE164458 (482).
- **Under brief 5.9.2's rule GSE164458 would be void** (rho 0.80 or above),
  leaving two cohorts.

**Recorded as unavailable. No refit** (the author, 2026-11-03).

**E39's three gate values go into ED-a** (12.3), beside TCGA's and SCAN-B's.
They are the record that the reference MYC score is barely separable from
proliferation in the neoadjuvant cohorts.

**The text consequence:**
- **the optional pCR MYC sentence (brief section 8) leaves the text;**
- **E40's MYC-adjusted OR, 0.87, stays in the text**, as the brief says.

## 10. Part C: sensitivities. New fits, declared here, REPORTED and NEVER READ

### 10.1 C1, TP53 status (TCGA only). AVAILABLE

**The variable.** `TP53_status`, from
`data/from_validation/tcga_brca_covariates.rds$covariates` (snapshot
`d3ac60e`).
- **Its definition.** The validation study built it in its `scripts/03` from
  MC3 non-silent calls, under its FILTER rule (`scripts/02`). It is `mutant`,
  `wildtype`, or NA where MC3 has no coverage.
- **The join.** It is joined on the patient barcode, as E34 joins PAM50,
  purity and leukocyte fraction. PAM50 agrees in 938 of 938 tumours across the
  join.
- **In the 938 tumours: 325 mutant, 613 wild-type, 0 missing**; all are
  MC3-covered.

**The model.** The 6 TCGA primary cells (3 genes x 2 measures, `M_a`, `m3`)
are refitted with `TP53_status` added, wild-type as reference:
`readout ~ OX + MYC + PAM50 + purity + leuko + PROLIF + TP53_status`.

**Reported:**
- **proliferation's coefficient and the difference, beside the primary;**
- TP53 status's own coefficient, which is saved with the full matrix,
  **printed and not read.**

**Gate.** The same code without the TP53 term must reproduce E47's stored
`m3` coefficients for those six fits, exactly.

**Recorded caveats:**
- **`TP53` mRNA is one of `PROLIF_DISJOINT`'s 318 genes.** Mutation status is a
  genomic variable and a different one.
- **MC3 calls come from the tumour DNA aliquot**, which may differ from the RNA
  vial (E33's aliquot note).
- **C1 asks whether proliferation's coefficient moves when TP53 status is held
  fixed.** It is not a mediation analysis.

### 10.2 C2, a second proliferation score. NOT RUN: a declared omission

**Dropped by the author, 2026-11-03.**

**What it would have used.** `PROLIF_STD` exists in both cohorts and contains
none of the twelve genes. But it is `PROLIF_DISJOINT` plus nine genes
(`CTPS1`, `DNMT1`, `HMGA1`, `NCL`, `NOP56`, `PRMT5`, `RANBP1`, `TFDP1` and
`UCK2`), and **all nine are in `M_a`** (`FELSHER__MITOSTRIP`, 61 genes).

**Why it is dropped: C2 would test the disjoint construction, not
proliferation.**
- With `M_a` in the model, it would reintroduce the shared measurement that
  D7 removed.
- So any change in proliferation's coefficient could come from the nine
  shared genes rather than from a second measure of proliferation.

**No second proliferation score is fitted in E48.** One free of `M_a`'s genes
would need its own declaration.

### 10.3 C3, "more steeply" at `m3`. AVAILABLE

**The code.** E47 section 11's bootstrap (script lines 683-724) is inline code,
not a function.
- E48 copies it with one change: the design matrix is built from `m3`'s
  covariates, which add `PROLIF`.
- **The seeds are `PROJECT_SEED` (20260831, E00) + k**, for k = 1 to 4 in
  E47's cell order. So each `m3` replicate resamples the same tumours as E47's
  `m1` replicate k.

**Gates:**
1. **The copied code, run at `m1`, must reproduce `e47$boot_replicates` and
   `e47$bcl2l1_bbc3` exactly.** This proves the copy is faithful.
2. **The `m3` design matrix must reproduce E47's stored `m3` OXPHOS
   coefficients** for both genes, within E47's tolerance.

**The estimate.** `delta = beta_OX(BCL2L1) - beta_OX(BBC3)` at `m3`, from
2,000 paired tumour-level resamples, with a percentile 95% CI.

**The reading, descriptive** (E47 section 9's rule): "more steeply" holds at
`m3` if `delta > 0` with its CI above zero in all four cells.
- **If it holds,** the text may say **"with the other two held fixed"**.
- **If not,** it says **"at fixed MYC activity"** and cites `m1`. TCGA's `m1`
  margins are thin, with lower bounds of 0.014 and 0.015.

### 10.4 C4, TCGA without the composition adjustment. Declared here, REPORTED and NEVER READ

**Added by the author, 2026-11-03.**

**Purpose.** Whether `MCL1`'s TCGA/SCAN-B split (2.9) follows the composition
adjustment rather than the cohort. TCGA's spine adjusts for purity and
leukocyte fraction; SCAN-B's cannot (trap 2).

**The model.** TCGA is refitted with SCAN-B's spine, PAM50 only, **on E34's
938 tumours**:

```
readout ~ OX + MYC + PAM50 + PROLIF
```

That is 3 genes x 2 OXPHOS measures x 2 MYC scores (`M_a` and the
least-entangled score), at `m3`: **12 fits.**

**Reported, beside the primary** (E47's `m3` with the full spine):
- OXPHOS's and proliferation's coefficients;
- the difference.

**It is never read, and it changes no label.**

**Gate.** The same code with purity and leukocyte fraction restored must
reproduce E47's stored `m3` coefficients for those 12 fits, exactly.

**Expectation (the author's).** **`MCL1`'s TCGA proliferation coefficient
moves towards SCAN-B's sign, which is positive, when the adjustment is
removed.**

**What it cannot do.** Removing two covariates does not make TCGA into SCAN-B:
their platforms and their tumours differ.
- **A move towards SCAN-B's sign** is consistent with the composition reading
  and does not establish it.
- **No move** points to the cohort rather than this adjustment.
- **Neither establishes a mechanism.**

## 11. Expectation, stated so that it can be wrong

**Restated by the author's amendment of 2026-11-03, from what the seen
coefficients imply (2.9). It replaces the brief's 5.7.**

### 11.1 How it is derived

E47's `m1` and `m3` fits are nested on the same tumours. **Adding
proliferation moves MYC's coefficient by minus proliferation's coefficient
times MYC's coefficient in the regression of proliferation on the other
terms** (`delta`):

```
beta_MYC(m3) - beta_MYC(m1) = - beta_PROLIF(m3) * delta
```

- **MYC activity predicts proliferation positively, so `delta > 0`.** The
  direction of MYC's `m1`-to-`m3` shift therefore fixes the sign of
  proliferation's coefficient.
- **The shifts were seen (2.2). `delta` was not computed.** So the signs are
  implied and the sizes are not.
- **One relation between sizes does follow.** `delta` depends only on the
  exposures and the tumours, so it is common to every gene in a cell. Within
  a cell, the implied sizes therefore stand in the ratio of the genes' shifts.
- **E48 checks the identity and `delta`'s sign in every cell** (G-A4).

### 11.2 The expectation, gene by gene

**`BBC3` and `BCL2L1`.**
- **Proliferation's coefficient is negative in all 8 cells on both MYC
  scores**, inferable from the seen shifts.
- **OXPHOS's is positive in every cell**, so the signs are opposite.
- **E48 adds the intervals and the within-model differences.**
- **Expected label: OPPOSITE for both.**
- **`BCL2L1` is the gene most likely to miss it on an interval.** On `M_a`
  its shifts are smaller than `BBC3`'s in all four cells (0.051-0.101 against
  0.108-0.191), so its implied coefficient is smaller in size in each.

**`MCL1`, on the reference score (`M_a`).**
- **Implied negative in TCGA, the same sign as OXPHOS's, and positive in
  SCAN-B. So not OPPOSITE.**
- **The label turns on TCGA's size, which was not inferable:**
  - **DISTINCT** if proliferation's TCGA coefficient is smaller in size than
    OXPHOS's, so that the difference keeps one sign in all four cells, with
    each CI excluding zero;
  - **NOT SEPARATED** otherwise.
- **Left open.**

**`MCL1`, on the least-entangled score.** The draft's "near zero in TCGA,
positive in SCAN-B" described the implied sign loosely. Precisely:
- **MYC's shift** is +0.005 and +0.003 in TCGA (level, priority), and -0.021
  and -0.024 in SCAN-B;
- **so proliferation's implied sign is negative, negative, positive,
  positive**, as on `M_a`.
- **TCGA's two implied coefficients are near zero in size.**

### 11.3 The agreement check

**Expected to hold in all 12 cells.** The implied signs on the least-entangled
score match `M_a`'s in every cell:
- **`BBC3` and `BCL2L1`:** negative throughout;
- **`MCL1`:** the same split, negative in TCGA and positive in SCAN-B.

**The weak point is `MCL1`'s two TCGA cells.** Their implied coefficients on
the least-entangled score are near zero, so agreement there rests on a sign
and not on a measurable association.

### 11.4 Methods text (draft; STYLE.md applies at stage 3)

> This analysis was post hoc, and the directions of the proliferation
> coefficients could already be inferred from previously reported
> coefficients for MYC activity. Their confidence intervals and the
> within-model differences were computed under a reading rule fixed before
> any of them was examined.

**This amendment changes neither the reading rule (section 7), the three genes
(section 3) nor the layouts (section 12).** The expectation is recorded so the
result can be read against it, as E43 14.5 and E47 section 11 were.

## 12. Displays, declared before any quantity is read

**The author fixes 4g's final design after seeing the three layouts with the
text** (brief 2.6). E48 renders all three.

### 12.1 The grammar, taken from E45 exactly

| element | encoding |
|---|---|
| **OXPHOS** | Okabe-Ito bluish green `#009E73`, a circle. **Open: without the proliferation term (`m1`). Filled: with it (`m3`)**, joined by a thin connector. This is 4h's meaning: before and after proliferation enters the model |
| **Proliferation** | vermillion `#D55E00`, a filled triangle (`m3` only; `m1` has no such term) |
| **MYC activity** (layout 2 only) | a square, open (`m1`) to filled (`m3`). **Its colour is fixed by rule now:** the first of `#0072B2`, `#CC79A7`, `#E69F00` that passes the palette checker against the other two |

**Rules for every panel:**
- **Connectors join only two rungs of one term.** Never different terms,
  cohorts or genes.
- **Every panel E48 renders, 4g and Extended Data alike, follows E45's fill
  convention.** E47's Displays A and B use the opposite one.

**The axis:**
- a vertical line at 0;
- the title "Association with expression (SD per SD of each score)";
- the cue "← lower expression | higher expression →" beneath it;
- **limits symmetric about 0, computed from the data to contain every CI,
  rounded outward to 0.1.** This is a rule, not a value.

**The columns:**
- **left:** "ρ, OXPHOS–proliferation", for each row's cohort and measure, from
  `g3`;
- **right:** "Proliferation − OXPHOS (95% CI)", the `m3` difference on `M_a`,
  as text.

**Row labels:**
- the genes "BBC3 (PUMA)", "BCL2L1 (BCL-xL)" and "MCL1";
- the sub-rows "TCGA · level", "TCGA · priority", "SCAN-B · level" and
  "SCAN-B · priority".

**No statistical notation on any panel**: no β, m1, m3, M_a or Δ.

**Size: 178 mm wide, with no text below 5.5 pt. Legends are descriptive
only.**

### 12.2 The three layouts. Their content is fixed now

**Layout 1, two scores (the baseline).**
- 3 genes x 4 blocks = 12 rows.
- OXPHOS open to filled, and proliferation filled.
- The ρ column on the left and the difference column on the right.

**Layout 2, three scores.**
- Layout 1, plus MYC activity (`M_a`) open to filled on every row.
- A second left column, "ρ, MYC–proliferation", per cohort.

**Layout 3, the whole argument in one panel.**
- **OXPHOS level only.** The clinical cohorts have no mitoPPS.
- **Three stacked facets**, each with its own x-axis and direction cue:

  | facet | rows | axis | right column |
  |---|---|---|---|
  | the BCL-2 family | 3 genes x 2 cohorts | SD per SD, 0 marked | the difference |
  | pCR | the three cohorts and the pooled estimate | odds ratio on a log axis, 1 marked | the ratio of ORs |
  | METABRIC | the ER-positive (IHC) primary row | hazard ratio on a log axis, 1 marked | the ratio of HRs |

- **The log-axis limits** follow the containing-and-rounding rule on the log
  scale.
- **A ρ column runs on every row.**
- **The clinical values are read from E45's `plot_frame`, and their ρ from
  E45's `exposure_rho`.**
- **The legend states the asymmetry.** The BCL-2 facet holds MYC activity
  fixed; the clinical facets come from models without MYC. **E40's
  MYC-adjusted OR, 0.87, goes in the text.**
- **If layout 3 is chosen,** the author decides whether 4h and 4i stay in the
  main figure.

### 12.3 Extended Data, rendered by E48

| | content |
|---|---|
| **ED-a** | Spearman correlations among the scores (5.4), per cohort and measure, with E39's bands. Each MYC score is named with its proliferation entanglement. **Beside TCGA's and SCAN-B's values: E39's gate values for the reference MYC score against proliferation in the three neoadjuvant cohorts** (0.7797 GSE25066, 0.7851 GSE194040, 0.8236 GSE164458), read from `e39_prepare_and_gate.rds$gate`. They are the record that it is barely separable from proliferation there (section 9) |
| **ED-b** | E47 Display C (`BCL2L1` against `BBC3`), relabelled in plain words, **values identical to `e47$displays$C_coefficients` and `$C_difference`**, with C3's `m3` result beside it |
| **ED-c** | MYC activity's coefficient on all four scores, ordered by entanglement, at `m1` and `m3`, for `BBC3`, `BCL2L1` and `MCL1`, with proliferation's `m3` coefficient beside it. **Void cells are drawn faded and labelled**, as in E47 Display B |
| **ED-d** | all twelve genes, three coefficients at `m3` (OXPHOS, proliferation and MYC activity), **on `M_a` and on the least-entangled score**, with the difference and its CI beside them. **No signature or arm rows** |
| **S-1** | the sensitivities C1 and C4 beside the primary, TCGA only: OXPHOS's and proliferation's coefficients and the difference, for the three read genes. **A display for the result note and the review pack; whether it enters the paper is the author's call** |

- **E47's Displays A and B remain E47's verified record.** For the paper, ED-c
  and ED-d replace them.
- **E46's re-render is open item 1** (brief section 10), not E48.
- **Co-amplification is brief section 9**, not E48.

## 13. What this cannot do

- **It is cross-sectional, on bulk transcripts, and post hoc.** Nothing is
  "primed" (N3).
- **Main effects only.** It does not test the mouse's MYC requirement, which
  acts at the death step (brief 4.6).
- **"Separable" is in E39's sense only.** `M_a` is the axis least separable
  from proliferation, and that is a finding (brief 4.1).
- **"The proliferation-associated component of MYC activity scores"**, never
  "the proliferation effects of MYC" (brief 4.2).
- **No link from the BCL-2 family to outcome, and no mediation** (brief 4.3).
  E40's `m3` term is close to uncorrelated with OXPHOS by construction. So its
  failure to move OXPHOS's coefficient says nothing about mediation.
- **The result note must name the p53 alternative for `BBC3`**, the
  composition alternative for `MCL1`, and `BCL2L1` as the weakest expected
  (brief 4.4).
- **Two cohorts, not pre-registered.**

## 14. Saved, per R rules 4 and 5, and the assertions

### 14.1 Saved

| element | contents |
|---|---|
| **the analysis frame** | one row per tumour: `sample_id`, cohort, `in_set`; both OXPHOS scores and all four MYC scores; `PROLIF` (Blom); `TP53_status` (TCGA); PAM50, purity and leukocyte fraction; the twelve genes |
| **the re-read coefficients** | OXPHOS and MYC at `m1` and `m3`, and proliferation at `m3`, for every gene fit read |
| **the contrast table** | the difference, SE, CI, df, the coefficient correlation and the naive SE, per cell extracted |
| **the label table** | the primary cells, per-gene labels, the agreement check, the G3 band and VIF per read cell |
| **the nine-gene tally** | 7.4 |
| **`part_b`** | "unavailable", with section 9's structure facts |
| **C1, C3, C4** | tables. For every new fit: the **full coefficient matrix**, the **full `vcov`** and the **full VIF vector** |
| **`part_c2`** | the omission and its reason (10.2) |
| **`nesting`** | G-A4: `delta` per cell and term, and the identity's maximum deviation |
| **`fit_digests`** | `coef()` and `vcov()` per new fit. **The fitted objects themselves are excluded from digests** |
| **C3's bootstrap replicates** | |
| **the correlation table** | ED-a, with E39's three gate values |
| **plotting frames** | 4g layouts 1 to 3, ED-a to ED-d, and S-1 |
| **legend drafts** | descriptive only |
| **`assertions`** | each recorded with what it compared and the maximum deviation |
| **`figure_files`** | names only |
| **`spec`** | the reading rule (section 7) and the outcome table (section 8) verbatim |
| **`built`** | |

**Digests are committed BEFORE the author's run.** `built` and the fitted
objects are excluded.

### 14.2 Assertions

**From the brief (5.11):**
1. **The re-read OXPHOS and MYC coefficients equal `e47$results` exactly**
   (G-A1).
2. **Each difference, recomputed independently from `coef` and `vcov`, equals
   the stored one.**
3. **Every drawn point equals its stored value**, read back from
   `ggplot_build()` (after the log transform in layout 3).
4. **No void cell is drawn as read.**
5. **Layout 3's clinical values equal E45's `plot_frame` exactly**, and its
   clinical ρ equal E45's `exposure_rho`.
6. **No interval is clipped by an axis.**
7. **The relabelled Display C (ED-b) has values identical to E47's display
   frames.**

**From this declaration:**

| assertion | what it checks |
|---|---|
| G-A2 | the stored proliferation and OXPHOS estimates, SEs and CIs reproduce from `coef`, `vcov` and df |
| `g3` | the `g3` pairs computed on `e47$frame` reproduce `g3` exactly |
| no signature reads | **no S6, S5, S7 or arm fit is read, and none is drawn** |
| C1 | the join is 938 of 938; PAM50 agrees after the label prefix; without the TP53 term, E47's six `m3` fits reproduce exactly |
| G-A4 | the `m1`-to-`m3` shifts of MYC and OXPHOS equal minus proliferation's coefficient times `delta`, to 1e-10, in every gene fit read |
| C4 | with purity and leukocyte fraction restored, the 12 TCGA fits reproduce E47's `m3` coefficients exactly |
| ED-a | E39's three gate values equal `e39_prepare_and_gate.rds$gate` |
| C3 | the `m1` copy reproduces `e47$boot_replicates` and `$bcl2l1_bbc3`; the `m3` design reproduces E47's `m3` OXPHOS coefficients |
| axes | the limits are symmetric, contain every CI, and are rounded outward to 0.1 |
| text and size | width 178 mm; no text below 5.5 pt; **no β, m1, m3, M_a or Δ in any label of any paper panel** |
| palette | the palette checker passes, with layout 2's colour chosen by 12.1's rule |
| encoding | open = without proliferation and filled = with it, in every E48 panel; connectors join only two rungs of one term |
