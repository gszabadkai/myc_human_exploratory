---
date: 2026-11-02
status: >
  AUTHOR'S DECISION RECORD AND BRIEF. Written in the claude.ai Project chat
  with the author, from the repo's notes and the Doc (read only; no repo
  access from the chat). It sets the narrative line for the human half of
  Fig. 4 and specifies E48. E48 is NOT declared, written or run. Nothing here
  is a new result. Commit it as written, alone.
decisions: 2026-11-02, by the author (section 2)
supersedes: >
  - the narrative LINE and the figure plan (section 6) of
    docs/2026-11-02_human_section_reconciliation_proposal.md. Its corrections
    are ADOPTED (section 1.4).
  - docs/2026-10-23_handoff.md section 2 for Fig. 4g, and its section 7 item 1
    (decision 5), within the scope of section 2 item 4 below.
relates-to:
  - docs/2026-10-23_e47_result.md (E47: MIXED on the signature; "more steeply" licensed at m1)
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (dcb58a7; addendum 3b5e221)
  - docs/2026-11-02_human_section_reconciliation_proposal.md
  - docs/2026-10-22_fig4_human_panels_proposal.md (sections 7-8: the E45 and E46 specifications)
  - the E43 and E44 notes (the proliferation coefficient recovered; the within-model contrast)
  - docs/2026-10-22_six_gene_composite_provenance.md (co-amplification numbers, 4.2 and 4.6)
  - the Doc "11.8.26_MK_Myc_Paper", Nat Metab version (Drive file ID 1M-kEC1d_VJmbTG-SyYYfykZ3F1_aEiDtu0xtwFaPtJo)
next-action: >
  Claude Code checks structure, drafts the E48 declaration from sections 5-7,
  and stops for the author's approval (see the session prompt).
---

# The separation line for Fig. 4's human panels, and the brief for E48

## 0. The line in three sentences

1. In human breast tumours, OXPHOS, MYC activity and proliferation rise together.
2. Modelled together, OXPHOS and proliferation run in opposite directions, on the BCL-2 family and on treatment response.
3. MYC activity scores carry proliferation to different degrees, so the apparent association of MYC activity with the BCL-2 family depends on how much proliferation its score carries.

**E48's job:** put the one missing quantity, proliferation's own coefficient on the BCL-2 family, on a declared footing, and build Fig. 4g so that it shows the line on its own.

## 1. Why this line

### 1.1 The sentence it replaces

"The apoptotic signature tracks OXPHOS but not MYC" (the abstract, the Fig. 4 human paragraphs and the closing paragraph).

- **"Not MYC" is contradicted.** E47 (the mutual-adjustment analysis of 23 October) found MYC at fixed OXPHOS above zero in 4 of 8 primary cells, and in 8 of 8 on the least-entangled MYC score. Its result note (9.1) withdraws "lost after adjustment for OXPHOS".
- **"Tracks OXPHOS" is constructed** at signature level, because the six genes and their signs were selected on OXPHOS in these same cohorts.
- **It is retired. Do not reintroduce it in any form.**

### 1.2 It is the paper's own line

The Doc (Nat Metab version, read 2026-11-02) already argues it in five places. The BCL-2 paragraph is the one place out of step.

- **Abstract, sentence 2:** "OXPHOS supports tumour proliferation by providing electron acceptors for biosynthesis, but whether it also constrains transformation remains unclear."
- **Abstract, sentence 9:** "The distinction between specific OXPHOS effects independent of its proliferative benefits **also** emerged clinically." The "also" promises a first instance in the BCL-2 result. That result currently adjusts for MYC instead.
- **Fig. 1a,c:** at constant MYC dose, the apoptotic response was lost in the adult gland whilst proliferation was maintained.
- **Fig. 1f:** the WT fall in OXPHOS was about five times the fall in proliferation, tested as OXPHOS minus proliferation (the bracket).
- **Closing paragraph:** "the role of OXPHOS in mammary tumourigenesis cannot be reduced to its support of proliferation."

**Estimand caution.** Fig. 1f and the human contrasts share a logical form, not an estimand. They may be paralleled in words, by direction only, and never drawn on a shared axis.

### 1.3 What the verified fits already show

| | without the proliferation term | with it | source |
|---|---|---|---|
| OXPHOS on the signature (S6, M_a). Order: TCGA GSVA / TCGA mitoPPS / SCAN-B GSVA / SCAN-B mitoPPS | +0.614 / +0.498 / +0.383 / +0.304 | +0.599 / +0.493 / +0.392 / +0.331 | E47 result §2: shift −0.015 to +0.027 |
| MYC activity on the signature, same cells | −0.076 / +0.046 / +0.033 / +0.127 | +0.026 / +0.175 / +0.165 / +0.269 | E47 result §2 and D2: shift +0.10 to +0.14 |
| OXPHOS on pCR, pooled | OR 1.05 [0.95, 1.17] (m0) | OR 0.85 [0.76, 0.96] (m1) | E40 |

**The reading is the arithmetic of omitted variables, not biology.** Proliferation correlates with both scores and runs against OXPHOS on both readouts. It hides in whichever correlated score is modelled without it:
- on the BCL-2 family it sat inside MYC activity, because E47's models always contain MYC;
- on pCR it sat inside OXPHOS, because E40's m0 has no MYC term.

**This is where "not MYC" came from.** Raw, the signature correlates with MYC activity at only 0.04 / 0.10 (TCGA / SCAN-B; reconciliation proposal 1.1, AD HOC). MYC's proliferation-associated component and its remainder cancel.

**The one quantity the line needs is not declared.** Proliferation's own coefficient on the signature is −0.22 to −0.38, CI excluding zero, in all 8 retrieved cells (E47 result 3.1). It is RETRIEVED, not declared. **Per gene, it has never been looked at.**

### 1.4 The two proposals of 2026-11-02

**Claude Code's reconciliation proposal (`docs/2026-11-02_human_section_reconciliation_proposal.md`).**

ADOPTED, its corrections:
- its section 5 removal table;
- its 3.2: the adult-gland *Bbc3*/*Bcl2l1* ratio points at the wrong path, so bridge to Fig. 3h,i;
- the ten-gene list holds for six genes;
- "restricted to apoptosis among four control pathways" comes out;
- the trigger-arm hypothesis stays out of the Discussion (its 8.4);
- the signature leaves the main text;
- both OXPHOS measures are reported every time.

SUPERSEDED:
- **its line** ("OXPHOS tracks PUMA and, more steeply, BCL-xL; MYC unresolved") **and its figure plan** (Display C as 4g). They leave the BCL-2 paragraph asking OXPHOS against MYC, whilst the abstract, Fig. 1 and the clinical paragraph ask OXPHOS against proliferation;
- **its pooled 44-specification count for MYC's coefficient.** It mixes rungs and scores that differ for a known, monotone reason. **Summarise MYC by score, and never pool across the entanglement gradient.**

**A chat proposal earlier the same day** claimed that at fixed OXPHOS, MYC activity loads the BH3-only arm and not BCL2L1. **It is RETRACTED.** It read M_a only. On the least-entangled score, BCL2L1 also rises with MYC at fixed OXPHOS (3 of 4 blocks) and MCL1 falls (4 of 4). Those values come from the reconciliation proposal 8.4, read from `e47$results`, and are recorded in no note. Do not use the claim anywhere.

## 2. Decisions taken by the author, 2026-11-02

1. **The separation line (section 0) is the narrative line** for the human half of Fig. 4.
2. **E48 goes ahead** (section 5).
3. **Fig. 4g becomes the separation panel** (section 6). E47 Display C (BCL2L1 against BBC3) and E46's quadrant map go to Extended Data.
4. **Decision 5 is taken, for co-amplification only.** The MYC/BCL2L1/MCL1 co-amplification result (validation study G2) enters the paper as one main-text sentence and one Extended Data panel, not as a main panel. Nothing else from the validation study is admitted by this decision. CLAUDE.md's "does not go in the paper" is amended to this scope (section 9).
5. **The Fig. 4 title will change.** The wording waits until the figure is final.
6. **The author fixes 4g's final design after seeing the drafts together with the text** (section 6.5).

## 3. The evidence, cell by cell

Each exposure is read with the others held fixed.

| exposure | BCL-2 family | pCR (three neoadjuvant cohorts) | METABRIC, breast-cancer-specific survival |
|---|---|---|---|
| **OXPHOS** | BBC3 and BCL2L1 up, MCL1 down, CI excluding zero in all 44 specifications (reconciliation proposal 1.3). BCL2L1 more steeply than BBC3: LICENSED at m1 (E47 §7). **Declared** | OR 0.85 [0.76, 0.96], I² 0, all three cohorts agree in sign. **Declared** (E40) | HR 0.896 [0.803, 0.9996]: no claim under the declared rule (E41, E43) |
| **Proliferation** | S6: −0.22 to −0.38 in 8 of 8 cells. **RETRIEVED.** Per gene: **never read; E48 Part A** | OR 1.72 [1.36, 2.18]. A post-hoc characterisation of a declared covariate (E43) | HR 1.229 [1.074, 1.406] (E43) |
| **MYC activity** | Depends on the score's proliferation content. On the least-entangled score, positive on S6 in 8 of 8 (E47 §3.2). **Descriptive** | E40's m2 adds MYC (OXPHOS OR 0.87). MYC's own coefficient was **never read; E48 Part B, if saved** | no MYC term fitted |
| **Proliferation against OXPHOS, within one model** | **never computed; E48 Part A** | ratio of ORs 1.93 [1.56, 2.39] (E44) | ratio of HRs 1.37 [1.13, 1.67] (E44) |

## 4. Rules for interpreting and writing

These rules apply to E48's declaration, its result note, every legend and every manuscript sentence.

1. **"Separable" in E39's sense only.** E39's bands: below 0.60 fine; 0.60 to 0.80 not cleanly separable; 0.80 or above not identifiable.
   - OXPHOS against proliferation (0.23 to 0.42) is separable everywhere.
   - M_a against proliferation (0.615 to 0.696) is not cleanly separable. Write "MYC activity was the axis least separable from proliferation". That is the finding, not a weakness.
2. **Write "the proliferation-associated component of MYC activity scores"**, never "the proliferation effects of MYC", in Results. The data cannot distinguish proliferation driven by MYC from proliferation with other causes, or from target genes the scores share. The causal version is a hypothesis for the closing paragraph, and Fig. 1a,c is its mouse support.
3. **No "thus", "because" or "explains" leading from the BCL-2 family to outcome.**
   - If E48 returns the expected signs, proliferation goes with lower BBC3 and higher MCL1. Read naively, that is less pro-apoptotic, yet proliferation predicts better response. The BCL-2 configuration cannot then explain proliferation's clinical association; the sensitivity of cycling cells to chemotherapy is the obvious alternative.
   - **The two separations run in parallel. The data do not chain them, and no mediation claim is made.**
   - E40's m3 term is the mean of z_BCL2L1 and z_MCL1, two genes that run in opposite directions with OXPHOS. It is close to uncorrelated with OXPHOS by construction, so its failure to move OXPHOS's coefficient says nothing about mediation.
   - The Doc's "or BCL2L1 expression" has no fitted model behind it: E40 has no BCL2L1-only rung.
4. **The result note must name these alternatives:**
   - **BBC3 and p53.** PUMA is a canonical p53 target, and proliferative breast tumours are enriched for TP53 mutation. A negative proliferation coefficient on BBC3 may reflect p53 status (Part C1).
   - **MCL1 and composition.** MCL1 is lymphocyte-expressed (provenance note section 7, limitation 8). Leukocyte fraction is adjusted in TCGA only, and coarsely, so compare TCGA (adjusted) with SCAN-B (not).
   - **BCL2L1** is the weakest expected: its raw correlation with proliferation is −0.13 / −0.03 (AD HOC).
5. **These are transcripts, not function.** Nothing is "primed" (N3). The data are cross-sectional, and the analysis is post hoc.
6. **Main effects only; no interaction was fitted.** The human data do not test the mouse's MYC requirement, which acts at the death step. The Doc already says so: "induction of cell death by this trigger was determined by the Myc status rather than by PUMA upregulation alone".
7. **Proliferation's coefficient is "a post-hoc characterisation of a declared covariate"**, as in the v5 Methods (statement 3) and the E43 precedent.
8. **Never set a mouse coefficient beside a human one.** No axis is shared between estimand types (SD per SD, odds ratio, hazard ratio).
9. **These words stay out:** "independently of MYC", "MYC expression" (trap 4), "comparable", "invariant", "lost after adjustment", "signature" for any MYC score, "primed", "priming".
10. **E47 stands.** Its MIXED verdict on S6 is not amended, and E48 does not read S6.

## 5. The E48 specification (to become the declaration)

### 5.1 Title and posture

**"E48: OXPHOS, proliferation and MYC activity on BBC3, BCL2L1 and MCL1, each with the others held fixed; and Fig. 4g."**

EXPLORATORY and POST HOC.
- **Part A recovers a declared covariate's coefficient, and within-model contrasts, from E47's verified saved fits** (coefficients and vcov), on the E43 and E44 precedents. **Part A fits nothing.**
- **Parts B and C are conditional** (5.9, 5.10).
- E48 also absorbs the reconciliation proposal's re-render of Displays A and B without the signature rows (section 7). The E46 re-render waits on open item 1.
- Confirm that E48 is still the next free script number.

### 5.2 Seen before writing: correct this list from your own record

**Seen:**
- (a) E47's result note in full, including proliferation's coefficient on S6 at m3 (−0.22 to −0.38, 8 retrieved cells, on M_a and on the least-entangled score).
- (b) E47's per-gene OXPHOS and MYC coefficients (result §5; Display A).
- (c) From the reconciliation proposal:
  - the AD HOC raw correlations of BBC3 and BCL2L1 with proliferation (−0.30 / −0.15 and −0.13 / −0.03);
  - the 44-specification counts;
  - the least-entangled score's MYC coefficients for BCL2L1, MCL1 and the guardian arm.
- (d) E40, E43, E44 and E45.

**Believed unseen:**
- per-gene proliferation coefficients in any E47 fit;
- any within-model contrast on E47 fits;
- the raw correlation of MCL1 with proliferation;
- E40's m2 MYC coefficient;
- anything stratified or adjusted by TP53.

### 5.3 Readouts

**READ: BBC3, BCL2L1 and MCL1.** Why these three, fixed before reading:
- BBC3 and BCL2L1 encode PUMA and BCL-xL, the proteins examined in the mouse (Figs 2 and 3);
- BCL2L1 and MCL1 are the primary BCL-2-family measures fixed by the 27 August validation plan, in its section 9;
- none of them is chosen by this analysis.

**REPORTED, not read:** the other nine genes. No per-gene FDR.

**NOT USED:** S6, S5, S7 and both arms, whose signs follow OXPHOS by construction.

### 5.4 Model and inputs

E47's fits, unchanged, read from `results/e47_mutual_adjustment.rds`:
- `readout ~ OX + MYC + spine (+ PROLIF_DISJOINT at m3)`, OLS, every variable Blom-scored per cohort;
- E34's analysis sets, TCGA 938 and SCAN-B 3,143;
- OXPHOS on GSVA ("level") and on mitoPPS ("priority");
- the MYC scores:
  - **M_a** (`FELSHER__MITOSTRIP`), the primary;
  - **`MYC_UP.V1_UP__MITOSTRIP`**, the least entangled, for the agreement check;
  - CollecTRI and HALLMARK, reported only. HALLMARK's m3 cells are void under G3.

### 5.5 Quantities

1. **Proliferation's coefficient at m3**, with proliferation, OXPHOS, MYC activity and the spine all in the model.
2. **The within-model difference, proliferation minus OXPHOS, at m3.**
   - Var = Var(prolif) + Var(OX) − 2 Cov(prolif, OX), from the saved vcov.
   - t-based 95% CI on the fit's residual df.
   - The orientation is E44's.
3. **Re-read, not new:** the OXPHOS and MYC coefficients at m1 and m3. Assert that they equal `e47$results`.
4. **Descriptive:** Spearman correlations among both OXPHOS scores, the four MYC scores and `PROLIF_DISJOINT`, per cohort. Read them from `e47$gates$g3` if it holds them; otherwise compute them on `e47$frame`.

### 5.6 Primary cells and the reading rule, FIXED BEFORE READING

**Primary cells:** 3 genes × 2 cohorts × 2 OXPHOS measures = 12 cells, on M_a at m3. Each gene is labelled over its 4 cells:

| label | condition |
|---|---|
| **OPPOSITE** | in all 4 cells, proliferation's CI excludes zero with the sign opposite to OXPHOS's, and the difference's CI excludes zero |
| **DISTINCT** | not OPPOSITE, but the difference's CI excludes zero in all 4 cells, with one sign |
| **NOT SEPARATED** | anything else |

- **The agreement check.** On `MYC_UP.V1_UP` at m3, proliferation's coefficient must carry the same sign as on M_a in all 12 cells. Any disagreement is printed beside that gene's label as ESTIMATOR-DEPENDENT, and no manuscript sentence may then state the label without naming the dependence.
- **G3 is carried.** M_a's correlation with proliferation (0.615 to 0.696) puts the M_a cells in "not cleanly separable", and the Methods says so. Print E47's VIF for every read cell.

### 5.7 The expectation, stated so that it can be wrong

**All three genes OPPOSITE**: proliferation negative for BBC3 and BCL2L1, positive for MCL1, with agreement in all 12 cells. Basis: 5.2(a) and (c). BCL2L1 is the weakest.

### 5.8 What each outcome licenses, fixed now

| result | text |
|---|---|
| all three OPPOSITE, agreement holds | abstract: "opposite associations"; section 8's sentence "Proliferation ran the other way ..." as drafted |
| none NOT SEPARATED, at least one DISTINCT | "distinct associations"; name the genes on which proliferation ran the other way |
| any gene NOT SEPARATED | that gene leaves the separation sentence. If it is BCL2L1, use the reconciliation proposal's §4 wording for BCL2L1. The clinical paragraph is unaffected |
| any ESTIMATOR-DEPENDENT | state the label with the dependence named, or do not state it |

**No outcome licenses** a causal, mediation, interaction or "primed" claim, or any link from the BCL-2 family to outcome.

### 5.9 Part B: MYC activity on pCR (conditional)

1. **Check structure only, before the declaration is committed:**
   - Do E40's or E43's saved objects hold the m2 fits (OX + PROLIF + MYC + subtype), with MYC's coefficient and vcov?
   - Which MYC score did E40 use?
   - What is E39's gate correlation of that score with proliferation in each pCR cohort?
2. **If the fits are saved:**
   - recover MYC's coefficient (log OR per SD) per cohort, void wherever ρ ≥ 0.80;
   - pool it by the E40 declaration's section 9 rule: fixed effect at I² < 50, DerSimonian–Laird random effects otherwise;
   - REPORT it; do not read it.
   - **The interpretation is fixed now.** A pooled CI including OR 1 means "MYC activity was not associated with response at fixed OXPHOS and proliferation". A CI excluding 1 means the direction and size are described, with no further claim.
3. **If the fits are not saved,** record "unavailable". No refit is made in E48; a refit needs its own declaration.
4. **METABRIC has no MYC term,** and nothing is fitted.

### 5.10 Part C: sensitivities (new fits, declared here, reported and never read)

- **C1, TP53 (TCGA only).** If TP53 mutation calls exist in the repo's TCGA snapshot, refit the 6 TCGA primary cells with TP53 status added, and report proliferation's coefficient and the difference beside the primary. If they are absent, record that and **download nothing.**
- **C2, a second proliferation score.** If `PROLIF_STD` (falsifier F2, E47 result §10) exists and passes G2 for the three genes, refit the 12 primary cells with it in place of `PROLIF_DISJOINT`, and report.
- **C3, "more steeply" at m3.** E47 §9 licensed BCL2L1 against BBC3 at m1 only, but the paragraph's frame is m3.
  - Repeat the comparison at m3 by E47's method: 2,000 paired tumour-level bootstrap resamples, seeded as E47 seeded, with both gene models refitted on each resample.
  - Report it beside m1, applying E47 §9's rule descriptively.
  - If it holds in all four cells, the text may say "with the other two held fixed". Otherwise it says "at fixed MYC activity" and cites m1. TCGA's m1 margins are thin (lower bounds 0.014 and 0.015).

### 5.11 Saved objects, assertions and digests

As in E47 (R rules 4 and 5):
- plotting frames, the contrast table, the label table and all re-read coefficients;
- the assertions, each recorded with what it compared;
- digests committed **before** the author's run, with fit objects excluded.

**The assertions must include:**
- the re-read OXPHOS and MYC coefficients equal `e47$results` exactly;
- each difference recomputed from coef and vcov equals the stored one;
- every drawn point equals its stored value, read back from `ggplot_build()`;
- no void cell is drawn as read;
- layout 3's clinical values equal E45's plotting frame exactly;
- no interval is clipped by an axis;
- the relabelled Display C has values identical to E47's display frame.

## 6. The Fig. 4g design brief

### 6.1 What it must do

**Simple, and strong.**
- Layouts 1 and 2 carry the BCL-2 half of the argument on their own, and complete it with 4h (pCR) and 4i (METABRIC) through a shared grammar.
- Layout 3 carries the whole argument in one panel.

**What a reader should take from it in five seconds:** OXPHOS and proliferation rise together (left column), yet with the others held fixed they sit on opposite sides of zero (markers), by a margin printed on the right.

### 6.2 The grammar, taken from E45 exactly

- **OXPHOS:** Okabe-Ito bluish green `#009E73`, a circle.
  - **Filled:** with the proliferation term (m3).
  - **Open:** without it (m1), joined to the filled marker by a thin connector. This is the same meaning as 4h's open marker: before and after proliferation enters the model.
  - Connect only two rungs of one term. Never join different terms, cohorts or genes.
- **Proliferation:** vermillion `#D55E00`, a filled triangle.
- **MYC activity (layout 2 only):** a third Okabe-Ito colour that passes the palette checker against the other two; a square, open (m1) to filled (m3).
- **NOTE: E47's Displays A and B use the opposite fill convention** (filled = m1). Every panel E48 renders, 4g and the Extended Data displays alike, follows E45: open = without proliferation, filled = with it.
- **The axis:**
  - a vertical line at 0;
  - the title "Association with expression (SD per SD of each score)";
  - a direction cue beneath it: "← lower expression | higher expression →";
  - limits symmetric about 0, computed from the data to contain every CI, rounded to 0.1. This is a rule, not a value.
- **The right column:** "Proliferation − OXPHOS (95% CI)".
- **The left column:** "ρ, OXPHOS–proliferation", for each row's cohort and measure.
- **Row labels in plain words:**
  - "BBC3 (PUMA)", "BCL2L1 (BCL-xL)", "MCL1";
  - sub-rows "TCGA · level", "TCGA · priority", "SCAN-B · level", "SCAN-B · priority", where level is GSVA and priority is MitoPPS, both defined in the legend.
- **No statistical notation on the panel:** no β, m1, m3, M_a or Δ.
- **Size and type:** 178 mm wide (E45's width), with no text below 5.5 pt. The legend is descriptive only.

### 6.3 Three layouts

Their content is fixed now. The author chooses after the run, as E47 §8a allowed for placement.

**Layout 1: two scores (the recommended baseline).**
- 3 genes × 4 blocks = 12 rows.
- OXPHOS open to filled, and proliferation filled.
- The ρ column on the left and the difference column on the right.

**Layout 2: three scores.**
- Layout 1, plus MYC activity (M_a) open to filled on every row.
- A second left column, "ρ, MYC–proliferation".
- It shows where proliferation was hiding.

**Layout 3: the whole argument in one panel.**
- OXPHOS level (GSVA) only, because the clinical cohorts have no MitoPPS.
- Three stacked facets, each with its own x-axis and direction cue:
  - the BCL-2 family: 3 genes × 2 cohorts, SD per SD, 0 marked;
  - pCR: the three cohorts and the pooled estimate, odds ratio on a log axis, 1 marked;
  - METABRIC: the ER-positive (IHC) primary row, hazard ratio on a log axis, 1 marked.
- The right column in each facet's own scale: the difference, the ratio of ORs and the ratio of HRs.
- A ρ column on every row.
- **The clinical values are read from E45's plotting frame**, and ρ from E45's `exposure_rho`.
- **A known asymmetry, which the legend states:** the BCL-2 facet holds MYC activity fixed, whereas the clinical facets come from models without MYC. E40's MYC-adjusted OR, 0.87, goes in the text.
- If layout 3 is chosen, the author decides whether 4h and 4i stay in the main figure or move to Extended Data.

### 6.4 Wireframes

These are schematic. **They contain no data, and marker positions are arbitrary.**

**Layout 1** (layout 2 adds an open-to-filled square per row, and a second ρ column):

```
                        ρ, OXPHOS–      ← lower expression | higher expression →    Proliferation − OXPHOS
                        proliferation                      |                            (95% CI)
BBC3 (PUMA)
  TCGA · level             r.rr          ▲                 |        ○──●              d.dd [l.ll, u.uu]
  TCGA · priority          r.rr          ▲                 |        ○──●              d.dd [l.ll, u.uu]
  SCAN-B · level           r.rr          ▲                 |        ○─●               d.dd [l.ll, u.uu]
  SCAN-B · priority        r.rr          ▲                 |        ○─●               d.dd [l.ll, u.uu]
BCL2L1 (BCL-xL)
  (four rows, as above)
MCL1
  (four rows, as above)
                                     −x.x                  0                  +x.x
                                 Association with expression (SD per SD of each score)

  ○──● OXPHOS, without → with the proliferation term          ▲ Proliferation
```

**Layout 3:**

```
                        ρ, OXPHOS–                                                  right column
                        proliferation
BCL-2 family (OXPHOS level)          ← lower expression | higher expression →      Proliferation − OXPHOS
  BBC3 · TCGA              r.rr                            |                         d.dd [l.ll, u.uu]
  BBC3 · SCAN-B            r.rr                            |
  BCL2L1 · TCGA            r.rr                            |
  BCL2L1 · SCAN-B          r.rr                            |
  MCL1 · TCGA              r.rr                            |
  MCL1 · SCAN-B            r.rr                            |
                                          SD per SD, 0 marked
Chemotherapy response                ← lower pCR | higher pCR →                    Proliferation / OXPHOS
  GSE194040                r.rr                            |                         ratio of ORs
  GSE164458                r.rr                            |
  GSE25066                 r.rr                            |
  Pooled                                                   |
                                          odds ratio (log axis), 1 marked
Breast-cancer-specific survival      ← longer survival | shorter survival →        Proliferation / OXPHOS
  METABRIC, ER+ (IHC)      r.rr                            |                         ratio of HRs
                                          hazard ratio (log axis), 1 marked
```

### 6.5 The review pack (stage 3)

One page per layout:
- the panel at final size;
- beneath it, section 8's text finalised against the verified numbers, with the panel citations that the layout implies;
- the legend draft;
- STYLE metrics for every paragraph.

Make it a PDF with Quarto or R Markdown if either is available; otherwise a markdown file with the PNGs. Save it under `docs/` with the date.

## 7. Extended Data displays rendered by E48

| | content |
|---|---|
| **ED-a** | Spearman correlations among the scores, per cohort and measure, with E39's bands; each MYC score named, with its entanglement |
| **ED-b** | E47 Display C (BCL2L1 against BBC3), relabelled in plain words with its values unchanged, and C3's m3 result beside it |
| **ED-c** | MYC activity's coefficient on all four scores, ordered by proliferation entanglement, at m1 and m3, for BBC3, BCL2L1 and MCL1, with proliferation's coefficient beside it |
| **ED-d** | all twelve genes, three coefficients (OXPHOS, proliferation, and MYC on M_a) at m3; no signature or arm rows |

E46's quadrant map is open item 1. Co-amplification is section 9, not E48.

## 8. Text drafts, to be finalised against E48

**Written against STYLE.md.**
- **Brackets** mark clauses that wait on E48 or on section 9.
- **The "Fig. 4g" citations** shift with the layout chosen; the review pack shows each version.
- **Every Extended Data panel needs a main-text citation,** because legends describe panels only.
- **The sentences to remove** are the reconciliation proposal's section 5 table, unchanged.

**Abstract, sentence 8** (replaces the "4,302" sentence):

> Extending this relationship to 4,081 human breast tumours[8,9], OXPHOS and proliferation had [opposite] associations with the BCL-2 family when modelled together with MYC activity, OXPHOS tracking BBC3 and, more steeply, BCL2L1, encoding BCL-xL.

Optional sentence 8b: *The association of MYC activity with these genes depended on how much proliferation its score carried.*

Sentence 9 stays as written; its "also" now connects. Sentence 10 is unchanged. **Sentence 11:**

> Together, these findings establish an apoptotic constraint on MYC-driven transformation and reveal that, in human breast tumours, OXPHOS is linked to apoptotic regulation and treatment response separably from proliferation.

**The last sentence of the mouse Fig. 4 paragraph** (replaces "established human cancers exhibiting elevated Myc expression ..."):

> Consequently, human tumours sustaining high OXPHOS were predicted to carry raised BCL-xL alongside raised PUMA, distinguishable from the consequences of the proliferation that OXPHOS also supports.

**The separation paragraph** (replaces the Doc's first human paragraph):

> In tumours, however, OXPHOS rises with both MYC activity and proliferation, so their associations can only be interpreted once the three are separated. In TCGA-BRCA (n = 938) and SCAN-B (n = 3,143), OXPHOS was scored as subunit expression and as its priority within the mitochondrial programme, as in the mouse gland (Fig. 1e). OXPHOS correlated moderately with proliferation and with MYC activity (Spearman ρ = 0.23–0.43), whereas MYC activity scores tracked proliferation more closely (ρ = 0.41–0.84; Extended Data Fig. Xa). We therefore modelled BBC3, BCL2L1 and MCL1 expression on all three together, adjusting for intrinsic subtype in both cohorts and for tumour purity and leukocyte fraction in TCGA. With the other two held fixed, higher OXPHOS was associated with higher BBC3 and BCL2L1 and with lower MCL1, in both cohorts and on both measures (Fig. 4g). [Proliferation ran the other way at fixed OXPHOS and MYC activity, being associated with lower BBC3 and BCL2L1 and higher MCL1 (Fig. 4g).] [Across the nine other BCL-2 family genes, OXPHOS and proliferation ran in opposite directions for n of them (Extended Data Fig. Xd).] The association of MYC activity accordingly depended on how much proliferation its score carried (Extended Data Fig. Xc). [On the MYC activity score least entangled with proliferation, that association ran in the same direction as OXPHOS.] The correlation of OXPHOS with proliferation therefore concealed opposite associations of the two with the BCL-2 family.

Notes on this paragraph:
- the ρ ranges come from the reconciliation proposal's reading of `g3`, and E48's ED-a must reproduce them, with every MYC score labelled;
- the second bracketed MYC sentence rests on values recorded in no note until E48 records them;
- the nine-gene sentence is descriptive, with no reading rule applied.

**The OXPHOS-side paragraph** (replaces the Doc's second human paragraph):

> Of the genes that rose with OXPHOS, BCL2L1 did so more steeply than BBC3 in both cohorts and on both measures [with the other two held fixed / at fixed MYC activity] (Extended Data Fig. Xb). This resembled the mouse tumours in which high OXPHOS persisted only alongside Bcl-xL (Fig. 3h,i). [At the level of copy number, MYC-amplified tumours carried an excess of BCL2L1 amplification across subtypes (Extended Data Fig. Xe). Its size depended on which measure of genomic instability was adjusted for, giving an odds ratio of 1.77 (95% CI 1.26–2.47) or 1.36 (0.98–1.88).] A similar state has been described in quiescent, chemotherapy-persistent triple-negative breast cancer cells, in which high OXPHOS accompanies dependence on BCL-xL[26].

Choose the bracket in the first sentence by C3's result. The co-amplification sentences depend on section 9's verification. Ref. 26 is the author's citation and has not been re-read here.

**The clinical paragraph, edits only.**

The opening second sentence becomes:

> Having separated OXPHOS from proliferation on the apoptotic machinery, we next asked whether the two are also separable in their association with chemotherapy response.

The composite/BCL2L1 sentence becomes the reconciliation proposal's wording:

> Additional adjustment for MYC activity, and then for the mean expression of BCL2L1 and MCL1, did not materially alter the OXPHOS–response association (pooled odds ratios 0.87 and 0.88).

If Part B yields a value, add one sentence, worded by 5.9's fixed interpretation:

> [MYC activity itself was not associated with response at fixed OXPHOS and proliferation (pooled odds ratio x.xx, 95% CI x.xx–x.xx).]

**Keep** "Their transcriptional coupling therefore concealed a clinically relevant distinction." It deliberately echoes the separation paragraph's last sentence.

**The closing paragraph, its human sentence:**

> In human breast tumours, OXPHOS and proliferation, although correlated, had opposite associations with both the BCL-2 family and chemotherapy response, whilst the association of MYC activity depended on how much proliferation its score carried.

**Metrics, measured on the two human paragraphs with every bracket included:** 15 sentences, 335 words, mean 22.3, median 23, range 15 to 29, 27% under 20 words. No sentence is under the 12-word floor, and there are no em-dashes and no serial commas. Re-measure after finalising.

**STYLE.md, if it is not in the repo:**
- Sentences: mean 22–25 words, median 21–24, range 12–38. **No sentence under 12 words.** Aim for 35–45% under 20 words; result-dense text has run near 30%.
- No em-dashes. British spelling, and "whilst". No serial comma.
- Colons only to introduce a list. (i)/(ii)/(iii) for enumerations inside a sentence.
- No one-sentence paragraphs, and no short emphatic sentences.
- Connectives: However, Moreover, Thus, Accordingly, Indeed, Notably, Altogether, In addition.
- Legends are descriptive only.

## 9. Co-amplification (decision 5): scope and constraints. NOT part of E48

- **Source:** validation study G2, in `myc_human_validation` (read-only, frozen at `d3ac60e`). These numbers come via the provenance note, 4.2 and 4.6.
  - **BCL2L1:** OR 1.77 (1.26–2.47) under `ANEUPLOIDY_SCORE`, and 1.36 (0.984–1.88) under `FRACTION_GENOME_ALTERED`. Both covariates were pre-specified, so **report both, never a bare OR.**
  - **MCL1:** OR 1.69 (1.17–2.44) at high-level amplification (+2), but **it inverts at any gain (≥ +1, OR 0.503).** It is basal-specific: Basal OR 3.12, Breslow–Day by TNBC p 0.022.
  - The Doc's "MCL1 co-amplification was restricted to the basal subtype" omits the threshold inversion. Either state the inversion, or drop the MCL1 clause, and ref. 27 with it if nothing else cites it.
- **Verify before any panel is built:**
  - none of these numbers has been verified from this repo (reconciliation proposal, section 6);
  - locate the G2 objects in the validation repo and check them;
  - read `handoff-2026-10-21_validation_audit.md` for anything bearing on G2. Its contradictions 1 and 2 concern H4, which this decision does not admit.
- **BID's pre-declared negative-control role** concerns the validation plan's priming endpoints, not G2. Because only G2 is admitted, the paper need not mention it.
- **The build is a separate display task,** specified after E48: one Extended Data panel and one main-text sentence (section 8).
- **CLAUDE.md:** amend "the validation study does not go in the paper" to "G2 co-amplification only: one ED panel and one main-text sentence (author, 2026-11-02)".

## 10. Open items

1. **E46 in Extended Data**, the author's call. Either:
   - as built, with the signature strip, in which case the E47 declaration's 4a.7 Methods paragraph goes into Methods; or
   - re-rendered without the strip, so that the signature leaves the paper (a display-only task, which can also switch E46 to E45's fill convention).
2. **Signature transparency:** whether Methods records that a signature built from these genes was analysed (E47), with its construction sentence.
3. **REFEREE_QUESTIONS, after E48:**
   - R8: no interaction was tested;
   - R9: the answer is Fig. 1a,c, where proliferation was maintained whilst apoptosis was lost;
   - a new entry on p53 and PUMA.
4. **The Fig. 4 title**, when the figure is final.
5. **Carried over:** LbNOX; the companion-paper figures; the legend pass on the mouse paragraphs; A7.

## 11. Do not

- read any E48 quantity before the declaration is committed: inspect structure only;
- refit E47, or amend its MIXED verdict;
- use S6, S5, S7 or the arms as readouts, or draw them in any paper panel;
- reintroduce the retracted trigger-arm claim (1.4);
- choose 4g's layout: render all three;
- change content or layout because of a value;
- download data (TP53), or edit the Google Doc;
- run numbered scripts except as redirected dry runs (Option A).
