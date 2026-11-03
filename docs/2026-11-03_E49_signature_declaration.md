---
date: 2026-11-03
status: >
  DECLARATION, pre-read, APPROVED WITH EDITS by the author (STOP 1,
  2026-11-03). No new E49 quantity has been computed, printed or plotted. That
  covers the signature's within-model differences, S5's and S7's
  proliferation coefficients, every effector quantity and every display.
  Section 2 lists what WAS seen, corrected from the record. Section 7.1
  records the structure checks, which read names, dimensions, term labels and
  counts only.
repo: myc_human_exploratory. myc_human_validation (d3ac60e) and myc_mouse are read-only.
governed-by:
  - docs/2026-11-03_signature_and_fig4g_brief.md (eb45a9b), the author's decision record, as amended by the author on 2026-11-03. Section 3 records the amendment.
  - docs/2026-11-02_separation_line_brief.md section 4. It governs every word written about results, with the new brief's section 5 added.
relates-to:
  - docs/2026-11-02_E48_separation_declaration.md (2941438; addendum 258a496). It is the format model, and the source of section 8's "nothing is discarded".
  - docs/2026-11-03_e48_result.md. Its Part 2 holds the twelve genes' values that 4g re-reads.
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (dcb58a7). Its addendum 4a.7 is the signature's exact rule.
  - docs/2026-10-23_e47_result.md (MIXED on S6, which stands)
  - docs/2026-10-22_six_gene_composite_provenance.md, sections 4 and 5
  - data/annotations/bh3_binding_code.csv and its README, new, committed with this declaration
decides:
  - Before any new E49 quantity is read, it fixes the readouts, the quantities, the gates, the displays, the placement rule, the binding-code annotation, and the saved objects and assertions. The displays are 4g in two versions and the Extended Data.
  - There is no reading rule and no label.
approval: >
  2026-11-03, by the author, with four edits:
  1. Claude Code's five additions are accepted: the signature's difference
     treated as circular (8.2); ED-j; both versions of every coefficient
     heatmap; the display details; the README for data/annotations/.
  2. Section 14 carries the author's two sentences for brief section 8, from
     the chat of 2026-11-03, with their basis.
  3. The section mark is not on the effectors: 8.4 adds a sentence, 9.1's
     legend note is reworded, and 12 records the decision and its reason.
  4. The design sketch is v2 (2.5). It is committed separately, as a picture
     of intent, never a source of values.
  This declaration, data/annotations/bh3_binding_code.csv and
  data/annotations/README.md are committed together, and nothing else.
next-action: >
  Stage 2: scripts/E49_signature_and_fig4g.R, two dry runs with every output
  redirected, and the digests, committed before the author's run.
---

# E49: the OXPHOS-associated BCL-2 signature beside OXPHOS, MYC activity and proliferation; the three effectors; and Fig. 4g

## Terms used from here on

| term | meaning |
|---|---|
| **OXPHOS level** | the `OXPHOS subunits` score on GSVA (VST input). "Level" in panels and text |
| **OXPHOS priority** | the same arm on mitoPPS (linear input). "Priority" in panels and text |
| **MYC activity score** (repo: **estimator**) | one of E47's four scores. **Never a "signature".** See the four rows below |
| **the reference score** | `M_a` = `FELSHER__MITOSTRIP`, 61 genes. It was the primary in E47 and E48 |
| **the least-entangled score** | `MYC_UP.V1_UP__MITOSTRIP`, 178 genes. E47 designated it by proliferation-gene share (1.5%) |
| **the HALLMARK score** | `HALLMARK_MYC_TARGETS_V1__MITOSTRIP` |
| **the CollecTRI regulon** | `M_b__MITOSTRIP` |
| **proliferation** | `PROLIF_DISJOINT` (318 genes), the term `PROLIF` in E47's fits |
| **m1, m3** | E47's rungs: `m1` = OXPHOS + MYC + spine; `m3` = `m1` + proliferation |
| **the signature** | **S6**, by E47 4a.7's exact rule: `BAD`+, `BIK`+, `BBC3`+, `BID`+, `BCL2L1`+, `MCL1`−. **Its variants:** S5, without `BID`; S7, with `BMF`−. **Membership is fixed** |
| **signature sign** | the direction in which a gene enters the signature, from its association with OXPHOS. It is not the gene's function |
| **the effectors** | `BAX`, `BAK1` and `BOK`. **They are outside every signature** |
| **4g's model** | for each row, the `m3` fit on the least-entangled score (amendment item 1, section 3) |
| **the reference model** | the `m3` fit on the reference score, as declared in E48 |
| **the difference** | proliferation's coefficient minus OXPHOS's, **within one `m3` fit** (E44's orientation) |
| **block** | one cohort-by-OXPHOS-measure combination. There are four |
| **the binding code** | "binds BCL-xL" and "binds MCL1", read from `data/annotations/bh3_binding_code.csv` |
| **separable** | in E39's sense only: below 0.60, fine; 0.60 to 0.80, not cleanly separable; 0.80 or above, not identifiable |

## 0. Why this exists

**The author's decision (brief section 1):**
- the six-gene signature returns as Fig. 4g's summary;
- 4g becomes a heatmap in Fig. 1g's layout;
- the text reads PUMA against BCL-xL.

**A declaration is needed before any of that can carry a sentence.**
- E47 result 3.1 says that the signature's proliferation coefficient was RETRIEVED and not declared, "and no manuscript sentence may rest on it without its own declaration". **E49 is that declaration.**
- It also fits the three effectors, so that 4g carries Fig. 1g's four groups.
- It renders the displays.

## 1. Posture

**EXPLORATORY and POST HOC.** It is written after E47 and E48, and in response to them.

- **Part A fits nothing.**
  - It recovers the signature's quantities from E47's verified saved fits (coefficients and `vcov`). The precedents are E43, E44 and E48.
  - It re-reads the twelve genes' values from E48's verified saved object.
- **Part B fits only the three effectors.** It uses E47's exact model, behind a gate that reproduces E47's `BBC3` fits (section 7.2).
- **There is no reading rule and no label** (brief 6.1). Everything is reported and may be described.
  - **The words are governed by the 2026-11-02 brief's section 4**, with the new brief's section 5 added.
  - The two claims retracted in brief 4.8 are used nowhere.
- **E47 stands, MIXED, and E48 stands**, with its declared reading and every saved result. **Neither is refitted**, except the effector fits and their reproduction gate.
- **The signature's membership is fixed.** The effectors are never added to it.
- **Main effects only. No product term is fitted.**
- **`E49` is confirmed as the next free script number.** On 2026-11-03 there was no `E49` script, branch or note.

## 2. What was seen before this was written

**The brief's 6.2 list, corrected from the record.** The record searched:
- this repo's notes;
- E47's and E48's result notes;
- section 2 of E48's declaration, which searched the session transcripts from 2026-10-08 onward;
- this session's own record, 2026-11-03.

### 2.1 The signature

- **Proliferation's coefficient on S6 at `m3`.**
  - It was printed once per cell, on 2026-10-23, for 8 cells: the reference and least-entangled scores, both cohorts, both measures.
  - E47 result 3.1 records the range, -0.22 to -0.38, with every CI excluding zero.
  - **It was not printed for the regulon or HALLMARK.**
- **S6's OXPHOS and MYC coefficients.**
  - Printed on 2026-11-02 for the reference and least-entangled scores, at `m1` and `m3`.
  - E47 result section 2 gives the reference score at all rungs. Its section 3.2 gives the least-entangled score: +0.215 to +0.424 at `m3`, with the `m1` values.
  - **Its MYC coefficient on all four scores at `m1` and `m3` was drawn in E47 Display B, which was viewed.**
  - E47 3.2 counts it: positive in 8 of 8 regulon cells, 7 of them with the CI above zero, and in 8 of 8 least-entangled cells.
- **S5 and S7.**
  - Their MYC coefficients on all four scores were drawn in Display B.
  - Their reference-score values and their ordering are in E47 section 4.
  - Their OXPHOS coefficients on the reference score were drawn in Display A.
- **S6's raw Spearman correlations** with OXPHOS level and priority, the reference and least-entangled scores, and proliferation were computed AD HOC on 2026-11-02. The reconciliation proposal, section 1.1, prints three pairs:

  | with | TCGA / SCAN-B |
  |---|---|
  | OXPHOS level | 0.55 / 0.40 |
  | reference score | 0.04 / 0.10 |
  | proliferation | -0.09 / -0.04 |

### 2.2 The twelve genes

- **E48 Part 2.** At `m3`, on the reference and least-entangled scores and in all four blocks: OXPHOS, proliferation and MYC activity, and the difference. OXPHOS and MYC activity at `m1`.
  - **These are 4g's twelve-gene values.**
  - **The author's design sketch is drawn from them** (section 2.5).
- **E48 ED-c and result 2.3:**
  - the three read genes' MYC coefficients at `m1` and `m3` on all four scores;
  - their proliferation coefficients at `m3`.
- **E47 Displays A and B; E10's partial correlations** (the E47 4a.2 table); **E22's coefficients.**
- **Raw Spearman correlations.**
  - Of `BBC3` and `BCL2L1` with the five scores in 2.1 (AD HOC, 2026-11-02).
  - Of all twelve with MYC activity on the full Felsher score, and with OXPHOS, on the full cohorts (Phase 1 D3 and the Phase 2 notes, 2026-09-01). Those scores and sets are not E49's.
- **AD HOC, 2026-11-03.** Purity's and leukocyte fraction's coefficients for `BBC3`, `BCL2L1` and `MCL1` in E47's TCGA `m3` reference-score fits (E48 result 1.7).

### 2.3 The effectors. The brief's "every effector quantity" unseen is NOT right

**Seen in 2026-09, on other scores and sets:**
- **raw Spearman correlations of `BAX` and `BAK1` with MYC activity**, on the full Felsher score (Phase 1 D3);
- **their correlations with OXPHOS**, in the Phase 2 S1 strata and its module table;
- **E10 scored the machinery genes**, `BOK` among them.
  - Its notes flag `HRK`, not `BOK`, as below the 25th expression percentile.
- **E15 drew their OXPHOS-against-MYC partial correlations.**
  - `BOK` is one of six genes whose difference disagrees in sign between cohorts.

**Unseen:**
- every effector coefficient from E47's model, at any rung and on any score;
- the effectors' within-model differences;
- their raw correlations on E34's analysis sets with the four scores;
- their expression percentiles on these matrices.

**The structure check (7.1) printed no effector value.**

### 2.4 What the seen numbers already imply

1. **The placement rule's outcome is known before E49 runs** (amendment item 4; section 8.3).
   - S6's MYC coefficient is positive in all 8 regulon cells and all 8 least-entangled cells (E47 3.2).
   - So at `m3` the two scores agree in sign in all four blocks, and **the MYC column stays in 4g.**
   - **The rule therefore cannot be called blind.** E49 re-reads the values and applies the rule as written.
2. **The Fig. 1i-style slope is known** (6.8).
   - By construction it equals `beta_OX(BBC3) - beta_OX(BCL2L1)` from E47's `m3` priority fits on the reference score.
   - That is minus C3's `m3` priority quantity, which E48's result note reports with its bootstrap interval (its 1.6).
   - **Only the decile means are new.**
3. **The signature's difference.**
   - Its point value can be roughly inferred from the seen OXPHOS and proliferation coefficients. Its interval cannot.
   - **It inherits the circularity of the OXPHOS term** (8.2).
4. **4g's twelve-gene cells are all seen values** (E48 Part 2).
   - The S6 strip's points on 4g's model were also seen, as part of E47's printed cells.
   - **What 4g adds that has not been seen:** the three effector rows.

### 2.5 The author's design sketch, v2

`docs/2026-11-03_fig4g_fig1g_layout_SKETCH_v2.png`, committed separately (the author, 2026-11-03).
- **It is a picture of intent, and never a source of values.**
- **It shows the amendment's three blocks**, with one model per row.
- **The effector rows are blank.**
- Its values, like the first sketch's, come from E48 Part 2.1, which section 2.2 lists as seen.
- **The first sketch, `docs/2026-11-03_fig4g_fig1g_layout_SKETCH.png`**, predates the amendment: it shows four column blocks. It was never tracked and is removed once v2 is in place.

## 3. The amendment, recorded here (the author, 2026-11-03)

**It arrived in the session prompt that opened E49.** It is applied in this declaration and is NOT written into the brief, which stays as committed (`eb45a9b`). Verbatim:

> Amendment to docs/2026-11-03_signature_and_fig4g_brief.md (the author, 2026-11-03). Apply it to the E49 declaration and record it there:
>
> 1. Fig. 4g has three exposure blocks, OXPHOS, MYC activity and proliferation, each split by cohort: six columns. All three come from one model per row, the m3 fits on the least-entangled MYC score (MYC_UP.V1_UP__MITOSTRIP), each coefficient with the other two held fixed. Reason: the reference score correlates with proliferation at 0.61-0.70, in E39's not-cleanly-separable band, so it cannot show MYC activity separated from proliferation; the least-entangled score was designated by proliferation-gene share in E47, before E48 was read.
> 2. The signature strip uses the same single model: OXPHOS (grey, defines the signature), MYC activity, proliferation.
> 3. The same heatmap from the reference-score model goes to Extended Data, beside ED-c's MYC ladder.
> 4. Placement rule, fixed now: report S6's m3 MYC coefficient on the CollecTRI regulon (re-read from E47). If its sign matches the least-entangled score's in all four blocks, keep the MYC column in 4g. Otherwise, 4g shows OXPHOS and proliferation only, from the reference-score model as declared in E48, and MYC activity appears only in Extended Data.
> 5. The 4g legend names the MYC score used and points to Extended Data for the others. The result note discloses that the CollecTRI regulon, though as weakly correlated with proliferation, does not share the least-entangled score's direction for BCL2L1 and MCL1 (E48 H2).
> 6. Brief section 8: add the score-choice sentence after the modelling sentence, and replace the MYC sentence, as drafted in the author's chat of 2026-11-03.

**What it changes in the brief:**

| brief | before | after |
|---|---|---|
| 0.2 and 7, columns | OXPHOS, MYC activity on two scores, proliferation | **three blocks, six columns**: OXPHOS, MYC activity (least-entangled), proliferation |
| 7, source of each column | OXPHOS and proliferation from the reference fits; each MYC column from its own score's fit | **all three from 4g's model**, one fit per row |
| 7, the strip | S6 in each of eight columns | **S6 from 4g's model**: OXPHOS (grey), MYC activity, proliferation |
| 6.8 | (none) | **adds the reference-model heatmap in Extended Data**, beside ED-c |
| (none) | (none) | **adds the placement rule** (8.3) |
| 7, legend | (none) | **names the MYC score used** and points to Extended Data for the others |
| 8 | the drafts | **a score-choice sentence and a new MYC sentence**, from the author's chat |

**The two section 8 sentences are not in the repo.** Stage 3 needs them from the author.

## 4. Readouts

- **S6, S5 and S7, exactly as E47 built them.** Read from E47's coefficient matrix and frame. **Their membership does not change.**
- **The twelve genes.** Read from E48's saved coefficients.
- **`BAX`, `BAK1` and `BOK`**, fitted new (5.2). **They stay outside every signature, whatever their association.** The Methods says so.
- **No arm is read or drawn.** The trigger-arm claim is retracted (the 2026-11-02 brief, 1.4).

## 5. Model and inputs

### 5.1 Part A, re-read

**E47's saved object**, `results/e47_mutual_adjustment.rds`, verified 22 of 22 at `9251355`. Its elements:
- `coefficient_matrix`: 816 fits, every term, with `n` and `df`;
- `vcov`: 816 matrices;
- `fit_digests`;
- `frame`;
- `gates$g3`.

**E48's saved object**, `results/e48_separation_and_fig4g.rds`, verified 33 of 33 at `1ccf805`:
- `coefficients`: the twelve genes on the reference and least-entangled scores;
- `c3` and `c3_boot_replicates_m3`;
- `frame`.

**The model, unchanged:**
- `readout ~ OX + MYC + spine (+ PROLIF at m3)`, OLS with an intercept;
- every variable Blom-scored per cohort;
- the spine is PAM50, purity and leukocyte fraction in TCGA, and PAM50 in SCAN-B;
- **E34's analysis sets: TCGA 938, SCAN-B 3,143.**

### 5.2 Part B, the effector fits

- **Expression:** `data/from_validation/tcga_brca_linear.rds` and `scanb_linear.rds`.
  - Linear DESeq2-normalised counts, genes by samples, taken to `log2(x + 1)`.
  - SCAN-B resolves through `scanb_pheno.rds$symbol_map`. **The three need no map entry** (7.1).
- **Scoring:** E34's recipe, verbatim as E47 copies it.
  - Blom over the whole cohort: 1,095 tumours in TCGA and 3,207 in SCAN-B.
  - The fits then subset to the analysis sets.
- **Exposures, spine and analysis sets** are read from E47's frame and never rescored.
- **The fits:** 3 effectors x 2 cohorts x 2 measures x 4 MYC scores x {`m1`, `m3`} = **96 fits.** `m0` is not fitted.
- **Everything else is E47's:**
  - `lm`, with t-based 95% CIs;
  - the full coefficient matrix, `vcov` and VIF vector of every fit;
  - the 96 fitted objects;
  - a `fit_digests` table.

## 6. Quantities

**All are at `m3` unless stated, in all four blocks.**

1. **S6, S5 and S7: proliferation's coefficient and the difference**, on the reference and least-entangled scores.
   - The difference is computed by E48's method: `L V L'`, with a t-based 95% CI on the residual df.
   - Proliferation's coefficient is also read on the regulon and HALLMARK, for the ladder (ED-j). HALLMARK's `m3` cells are void.
2. **S6, S5 and S7: MYC activity's coefficient at `m1` and `m3`, on all four scores**, re-read from E47. HALLMARK's `m3` cells are void.
3. **S6, S5 and S7: OXPHOS's coefficient**, in every fit read. **It is displayed grey and never read** (8.2).
4. **Raw Spearman correlations**, per cohort, on the analysis sets. They are descriptive.
   - **Of:** S6, S5, S7 and all fifteen genes.
   - **With seven exposures:** OXPHOS level, OXPHOS priority, each of the four MYC scores, and proliferation.
5. **The effectors:**
   - OXPHOS's, MYC activity's and proliferation's coefficients at `m1` and `m3`;
   - the difference at `m3`;
   - on all four scores, in all four blocks. HALLMARK's `m3` cells are void.
6. **The effectors' expression levels**: each one's median percentile among all genes, per cohort.
   - They are reported beside the coefficients. This is the check E10 made for the twelve at the 25th percentile.
   - **It is reported, not a gate.**
7. **The placement rule's inputs** (8.3): S6's `m3` MYC coefficient on the regulon and on the least-entangled score, in all four blocks, and the rule's outcome.
8. **The Fig. 1i-style display**, per cohort (brief 6.4.6):
   - **the relation:** `Blom(BBC3) - Blom(BCL2L1)` against OXPHOS priority, with reference MYC activity, proliferation and the spine held fixed. Both sides are residualised on those terms;
   - **what is drawn:** the residuals' decile means and the fitted slope;
   - **the slope's interval:** C3's bootstrap interval, negated to this orientation and read from E48;
   - **also saved and displayed:** the partial correlation.

## 7. Gates and structure checks

### 7.1 Structure, checked 2026-11-03 (names, dimensions, term labels and counts only)

| object | what it holds |
|---|---|
| **E47** | **All 144 signature fits**: S6, S5 and S7, at `m0`, `m1` and `m3`, on 2 cohorts x 2 measures x 4 scores. Each has a full coefficient vector with an intercept, `n` and `df`, and a square `vcov` named by its terms in the same order. **`PROLIF` is in all 48 `m3` fits and in none of the 96 others.** HALLMARK's `m3` signature cells are flagged void. `fit_digests` covers all 816 fits. `fits_S6` holds the 48 S6 objects. The frame holds S6, S5, S7, the arms, the twelve genes and every exposure, for 4,302 tumours, 4,081 of them in the analysis sets |
| **the expression matrices** | **`BAX`, `BAK1` and `BOK` are present in both**, under E34's recipe. TCGA is 18,115 genes x 1,095 tumours and SCAN-B 18,153 x 3,207; both are linear DESeq2-normalised and hold every cohort sample. **Each effector resolves to exactly one row, with no symbol-map entry needed.** Each has no NA and is non-constant over the cohort, a TRUE/FALSE check only. As a control, the twelve resolve the same way. **None of the three sits in any of the four MYC scores, in `OXPHOS subunits` or in `PROLIF_DISJOINT`** |
| **E48** | `coefficients` holds the twelve genes on the reference and least-entangled scores, 48 rows each: 12 genes x 4 blocks. Every column 4g needs is present, with no NA: OXPHOS, MYC and proliferation at `m3`, with CIs; OXPHOS and MYC at `m1`; and the difference. `c3` holds `m1` and `m3` on both measures, and `c3_boot_replicates_m3` holds four replicate sets. **E48's frame equals E47's on all twelve genes and on the five exposures compared**, after aligning by sample. Those five are both OXPHOS scores, the reference and least-entangled scores, and proliferation. It holds the analysis sets only |

### 7.2 Gates in the script, before any quantity is used

| gate | what it checks |
|---|---|
| **G-A1** | the re-read S6, S5 and S7 coefficients and `vcov`s equal E47's: digest equality against `fit_digests`, every fit read |
| **G-A2** | each difference is recomputed independently from `Var(P) + Var(O) - 2 Cov(P, O)`, and equals the `L V L'` contrast |
| **G-A3** | the twelve genes' values used in every display equal E48's saved `coefficients` exactly |
| **G-B1** | the gene-scoring code, run on the twelve, reproduces E47's frame exactly: maximum absolute difference below 1e-10 in both cohorts. **The effectors are scored only by the function that passed G-B1** |
| **G-B2** | the fitting code, refitting `BBC3` in all 32 combinations (2 cohorts x 2 measures x 4 scores x {`m1`, `m3`}), reproduces E47's coefficient and `vcov` digests exactly. **If G-B2 fails, no effector is fitted** |
| **G-B3** | the effectors are in none of the six fitted predictor sets (G2's rule) |
| **G-B4** | G3: the effector cells void wherever `rho(MYC, PROLIF)` is 0.80 or above, which is HALLMARK's `m3` |
| **G-C1** | the raw S6 correlations equal the reconciliation proposal's to 2 dp: level 0.55 / 0.40, reference score 0.04 / 0.10, proliferation -0.09 / -0.04 |
| **G-C2** | the Fig. 1i-style slope equals `beta_OX(BBC3) - beta_OX(BCL2L1)` from E47's `m3` priority fits, to 1e-10. The interval drawn equals C3's, negated, exactly |
| **G-C3** | the placement rule is evaluated from E47's coefficient matrix, and its inputs and outcome are printed and saved |
| **membership** | S6, S5 and S7 equal E47's `SIGS`, and no effector is in any of them |

**The display assertions, as E48's, plus:**
- **values and drawing:**
  - every drawn value equals its stored value, read back from `ggplot_build()`;
  - the two 4g versions draw identical values;
  - **version (b) has a dot exactly where the CI excludes zero**;
- **the heatmap scale:**
  - E46's diverging scale, `#2166AC` through white to `#B2182B`;
  - limits symmetric, containing every coloured value and rounded outward to 0.1;
- **the grey cells**, the signature's OXPHOS cells, are drawn in a grey that is not on the scale, and none of them sets the limits;
- **the left-hand columns:**
  - the binding code is drawn only on BH3-only rows and only where the table holds `TRUE` or `FALSE`. A blank row is drawn blank;
  - the signature-sign column equals E47's `SIGS` for S6;
  - § appears after `BAD` and `BID` only;
- **size and type:**
  - the width is 178 mm, and no text is below 5.5 pt;
  - no β, m1, m3, M_a or Δ appears in any label;
  - every colour comes from the validated palette.

## 8. No reading rule: what is reported, and what governs the words

### 8.1 Everything is reported

- **There is no label** (brief 6.1).
- **Nothing is discarded** (E48 section 8). Every computed quantity is saved, tabulated and displayed, and may be described.

### 8.2 Circularity

- **The signature's OXPHOS coefficient is never read**, for S6, S5 or S7 alone (brief section 2 and 5.1). It is drawn grey and labelled "defines the signature".
- **Added here: the signature's difference inherits that circularity.**
  - The difference is proliferation's coefficient minus OXPHOS's, so it contains the OXPHOS coefficient.
  - It is computed, saved and displayed, flagged as such.
  - **No sentence may use it as evidence that OXPHOS and proliferation separate on the signature.**
  - What may be said about the signature and proliferation rests on proliferation's coefficient alone, worded raw and adjusted (brief section 2).

### 8.3 The placement rule (amendment item 4), fixed now

| condition | 4g | MYC activity |
|---|---|---|
| **S6's `m3` MYC coefficient on the regulon has the same sign as on the least-entangled score in all four blocks** | **three blocks from 4g's model**, as section 9.1 draws it | in 4g |
| **otherwise** | **OXPHOS and proliferation only, from the reference model** (as declared in E48), with the strip likewise | Extended Data only. Both full three-block heatmaps go to Extended Data |

- **Its outcome is already implied by E47 3.2** (2.4, item 1).
- **The script implements both branches and asserts which one it drew.**

### 8.4 The words

- **The governing sections:**
  - **the 2026-11-02 brief's section 4 governs**, with the new brief's section 5 added;
  - **the two retractions in brief 4.8 are never used**, in any document.
- **No inference beyond the data:**
  - **binding classes are never inferred**;
  - **no sentence infers a net apoptotic load** from the `BAD`, `BID`, `BIM` or `MCL1` transcripts (brief 5.3);
  - **`BAD` is reported with its explanation** in the same sentence or the next (brief 5.4);
  - **No sentence reads the effectors' transcripts as their activation.** (The author, 2026-11-03.)
- **What the text and legend must carry:**
  - **4g's legend names the MYC score used and points to Extended Data for the others** (amendment item 5);
  - **the result note discloses** that the regulon, though as weakly correlated with proliferation as the least-entangled score, does not share its direction for `BCL2L1` and `MCL1` (E48 H2);
  - **any pattern visible in 4g that the text's summary does not cover gets one descriptive clause** (brief 5.6).

## 9. Displays, declared before any quantity is read

**The author chooses between 4g's two versions at STOP 3.** Every heatmap of model coefficients is rendered in both versions:
- 4g;
- ED-e and ED-f;
- ED-h's adjusted block.

So the author's choice carries to every one of them.

### 9.1 Fig. 4g (main)

**Rows**, in Fig. 1g's groups and order, with vertical group labels as in Fig. 1g:

| group | genes |
|---|---|
| Effectors | `BAX`, `BAK1`, `BOK` |
| Activators | `BID`, `BBC3` (PUMA), `BCL2L11` (BIM) |
| Sensitisers | `BAD`, `PMAIP1` (NOXA), `BMF`, `BIK` |
| Anti-apoptotic | `MCL1`, `BCL2L2` (BCL-w), `BCL2A1` (A1), `BCL2`, `BCL2L1` (BCL-xL) |

**Left of the heatmap:**
- **"binds BCL-xL" and "binds MCL1"**: filled for `TRUE` and open for `FALSE`.
  - They appear on BH3-only rows only.
  - They are blank where the table is blank, so `BID`, `BIK` and `BMF` are blank until the author fills them.
- **"signature sign"**: + or −, for S6's six genes only.
- **§ after `BAD` and `BID`.**

**Heatmap columns** (amendment item 1):
- **OXPHOS | MYC activity | proliferation**, each split into TCGA and SCAN-B: **six columns**;
- OXPHOS level, at `m3`;
- **all three from 4g's model, the least-entangled score, each coefficient with the other two held fixed.**

**Colour:**
- E46's diverging scale (`#2166AC`, white, `#B2182B`);
- limits symmetric about zero, containing every coloured value, rounded outward to 0.1;
- **not Fig. 1g's brown to green**, because the estimand differs.

**The signature strip**, beneath and aligned to the six columns:
- **what it shows:** S6 from 4g's model, its estimate and 95% CI as text in each cell;
- **the OXPHOS cells:**
  - grey, off the colour scale;
  - their value in grey text;
  - labelled "defines the signature".

**The two versions:**
- **(a)** unmarked;
- **(b)** a small dot on each gene cell whose 95% CI excludes zero.

The strip carries its CIs as text in both.

**Size:** 178 mm wide, with no text below 5.5 pt.

**The legend is descriptive only.** It carries:
- the § note: "§ BH3-only proteins whose activity is switched after translation (BAD: off by phosphorylation; BID: on by cleavage to tBID)" (the author, 2026-11-03);
- the binding code's source, and that blank means not yet verified;
- **the MYC score used, with a pointer to Extended Data for the others**;
- OXPHOS level, with priority in Extended Data;
- the spine;
- n per cohort.

### 9.2 Extended Data, rendered by E49

**The letters continue E48's a to d. All of them are placeholders.**

| | content |
|---|---|
| **ED-e** | **4g on OXPHOS priority**, 4g's model (brief 6.8) |
| **ED-f** | **4g from the reference model**, OXPHOS level, with its strip. It sits beside ED-c's MYC ladder (amendment item 3) |
| **ED-g** | **S5, S6 and S7 strips side by side**: on 4g's model and the reference model, on both measures |
| **ED-h** | **raw correlations beside the adjusted coefficients**, for S6, S5, S7 and the fifteen genes. It shows all seven exposures' raw Spearman correlations per cohort, and beside them 4g's coefficients |
| **ED-i** | **intervals for all fifteen genes**: E48's ED-d extended with the effectors, in Fig. 1g's groups. OXPHOS, proliferation and MYC activity at `m3`, with the difference, on the reference and least-entangled scores, in all four blocks. **The twelve are asserted equal to E48's** |
| **ED-j** | **added here: the MYC ladder for S6, S5, S7 and the effectors**, in ED-c's format. MYC activity at `m1` and `m3` on all four scores, with proliferation at `m3`; void cells faded. **It displays quantities 6.2 and 6.5, which no listed display shows, and it is what a sentence on the signature's MYC ladder would cite** |
| **ED-k** | **the Fig. 1i-style human panel**, per cohort, on OXPHOS priority (6.8) |

- **Kept from E48, unchanged:**
  - ED-a, the score correlations;
  - ED-b, `BCL2L1` against `BBC3`;
  - ED-c, the MYC ladder;
  - S-1, the sensitivities.
- **Kept from E46:** the quadrant map, as built (brief decision 6).
- **Retired from the paper:** E48's layouts 1 to 3, which stay as E48's record.

## 10. The binding-code annotation

**`data/annotations/bh3_binding_code.csv`**, created with this declaration, with a provenance `README.md`.
- **Columns:** `gene`, `binds_BCL2L1`, `binds_MCL1`, `source` and `verified_by`.
- **The verified rows**, from the abstract of Chen et al. 2005 (PMID 15694340), each quoted in its row:

  | gene | binds BCL-xL | binds MCL1 |
  |---|---|---|
  | `BBC3` | `TRUE` | `TRUE` |
  | `BCL2L11` | `TRUE` | `TRUE` |
  | `BAD` | `TRUE` | `FALSE` |
  | `PMAIP1` | `FALSE` | `TRUE` |

  Claude Code re-checked all four against the PubMed abstract on 2026-11-03.
- **`BID`, `BIK` and `BMF` are blank.** The author supplies their classes from the paper's full binding data. **Until then nothing is drawn for them, and a class is never inferred.**
- **The script reads the table and asserts what it draws** (7.2).

## 11. Expectation, stated so that it can be wrong

**What can still be wrong:**
- **S5's proliferation coefficient:** negative in every cell. It shares five of S6's genes, and only `BID`, whose own coefficient is small, leaves.
- **S7's proliferation coefficient:** negative in TCGA. Smaller in size in SCAN-B, where `BMF`, entering with −, has a proliferation coefficient that changes sign between cohorts (E48 Part 2).
- **The effectors' OXPHOS coefficients:** **`BAX` and `BAK1` rise with OXPHOS** at fixed MYC activity and proliferation, in both cohorts. The basis is the raw correlations seen in 2026-09 (2.3).
- **No expectation is stated** for `BOK`, or for any effector's MYC or proliferation coefficient.

**What is already known, and recorded, not expected (2.4):**
- the Fig. 1i-style slope is negative in both cohorts;
- the placement rule keeps the MYC column.

## 12. What this cannot do

- **It is cross-sectional, on bulk transcripts, and post hoc.** Nothing is "primed" (N3).
- **Main effects only.** It does not test the mouse's MYC requirement, which acts at the death step.
- **The signature's OXPHOS association is circular**, and its difference inherits that (8.2).
- **Some transcripts do not measure activity:**
  - **`BAD` and `BID`'s transcripts do not measure their activity** (the §);
  - **nor do the effectors' transcripts measure their activation.** `BAX` and `BAK` are activated after translation.
    - **The decision (the author, 2026-11-03): the § is NOT put on the effectors.**
    - **The reason.** The § flags BH3-only genes whose transcripts a reader might take as a pro-apoptotic load. The effectors' transcripts are never read at all: no sentence reads them as activation (8.4).
- **The binding code comes from peptide affinities** (Chen et al. 2005), not from these tumours.
- **Two cohorts, not pre-registered.**

## 13. Saved, per R rules 4 and 5, and the assertions

| element | contents |
|---|---|
| **the analysis frame** | E47's frame, plus the three effectors' Blom values. One row per tumour, with `in_set` |
| **the re-read tables** | S6, S5 and S7: every coefficient read, with its fit identity. The twelve genes as used |
| **the contrast table** | the difference for S6, S5, S7 and the effectors: SE, CI, df, the coefficient correlation and the naive SE |
| **the effector fits** | for all 96: the full coefficient matrices, `vcov`s and VIF vectors; the fitted objects; `fit_digests` |
| **the raw correlations** | 18 readouts x 7 exposures x 2 cohorts |
| **the effector expression percentiles** | 6.6 |
| **the placement rule** | its inputs, its outcome, and the branch drawn |
| **the Fig. 1i-style frame** | residuals, deciles, the slope, the partial correlation and C3's interval |
| **the binding code, as read** | |
| **plotting frames** | 4g in both versions, and ED-e to ED-k |
| **legend drafts** | descriptive, with the § note and the score named |
| **`assertions`** | each with what it compared and the maximum deviation |
| **`figure_files`, `spec` and `built`** | `spec` holds the posture, the placement rule and the amendment, verbatim |

**Digests are committed BEFORE the author's run.** `built` and the fitted objects are excluded.

## 14. The text

**What stage 3 finalises against the verified numbers:**
- **brief section 8's drafts, with amendment item 6 applied** (14.1);
- **one descriptive sentence on the effectors**;
- **the Methods notes of brief section 8.**

### 14.1 Amendment item 6: the author's two sentences (chat of 2026-11-03)

**After the modelling sentence of the signature paragraph:**

> Two of four MYC activity scores, including the widely used Hallmark set, correlated with proliferation too closely to be separated from it, so we used one containing few proliferation genes (Extended Data Fig. Xa).

- **It does not claim "the fewest" proliferation genes.** The regulon's share is not on record.

**Replacing the MYC sentence:**

> On the score with few proliferation genes, and on a second, equally separable score, MYC activity went with the signature in the same direction as OXPHOS, whereas the other two scores showed this less (Fig. 4g and Extended Data Fig. Xj).

- **Its basis is E47 3.2.** The signature's MYC coefficient is positive in all 8 regulon cells and all 8 least-entangled cells.
- **It concerns the signature only.**
- **The regulon's different direction for `BCL2L1` and `MCL1` stays disclosed**, as amendment item 5 requires.

**What governs it:**
- STYLE.md;
- the 2026-11-02 brief's section 4, with the new brief's section 5.
