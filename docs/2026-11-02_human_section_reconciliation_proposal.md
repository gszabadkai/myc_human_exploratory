---
date: 2026-11-02
status: >
  PROPOSAL for the author, written at the author's request. It answers: what
  do we say about MYC, OXPHOS and the BCL-2-family genes in human tumours,
  and how do we present it? NO ANALYSIS was run for it and no script was
  written. Every number in the drop-in text (section 4) is a saved quantity
  in a verified object. A few numbers in sections 1 and 2 were computed
  read-only from E47's saved frame, for the author's understanding. They are
  marked AD HOC, their code is in section 9, and they are not for the
  manuscript.
reads:
  - the Google Doc "11.8.26_MK_Myc_Paper" (Nat Metab version), as read on 2026-11-02. Read, never edited
  - results/e47_mutual_adjustment.rds (verified 22 of 22)
  - results/e46_fig4_quadrant_panels.rds (verified 10 of 10)
updates: >
  docs/2026-10-23_handoff.md sections 5 and 7. Section 8 below maps each of
  their items to what this note proposes.
relates-to:
  - docs/2026-10-23_e47_result.md
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (4a.2, 4a.7)
  - docs/2026-10-22_six_gene_composite_provenance.md
  - docs/2026-10-22_fig4_human_panels_proposal.md
next-action: >
  The author amends the Doc's human paragraphs from section 4, then decides
  the figure plan (section 6) and the four items in section 8. Figure-panel
  proposals follow in a later session.
---

# The human section: what to say about MYC and OXPHOS, the text, and the figures

## 0. The verdict

**Retire "the apoptotic signature tracks OXPHOS but not MYC"** from the
abstract, the Fig. 4 human paragraphs and the conclusions. Two sentences
replace it:

1. **OXPHOS tracks `BBC3` (PUMA) and, more steeply, `BCL2L1` (BCL-xL), after
   adjustment for MYC activity.** This holds in both cohorts and on both OXPHOS
   measures, with or without proliferation adjustment, on all four MYC
   activity estimators. The CI clears zero in all 44 model specifications for
   both genes.
2. **MYC activity's own association with these genes is unresolved, not
   absent.** Its direction depends on which MYC estimator is used and on
   whether proliferation is accounted for.

Sentence 1 carries the abstract, the Results, Fig. 4g and the conclusion.
Sentence 2 needs one Results sentence pointing to Extended Data. **Neither
contradicts the mouse** (section 3).

## 1. The numbers behind it

### 1.1 Raw Spearman correlations, TCGA / SCAN-B

Analysis sets: TCGA 938, SCAN-B 3,143. OXPHOS is the level score (GSVA); MYC
activity is the reference estimator, `M_a`.

| | OXPHOS (level) | MYC activity | proliferation |
|---|---|---|---|
| OXPHOS (level) | — | 0.39 / 0.43 | 0.24 / 0.40 |
| MYC activity | | — | 0.70 / 0.61 |
| `BCL2L1` | 0.38 / 0.31 | -0.11 / 0.00 | -0.13 / -0.03 |
| `BBC3` | 0.27 / 0.10 | -0.09 / -0.01 | -0.30 / -0.15 |
| signature (S6) | 0.55 / 0.40, *by construction* | 0.04 / 0.10 | -0.09 / -0.04 |

- **The first two rows are saved** in `e47$gates$g3`.
- **The gene rows are AD HOC** (section 9): for the author, not for the paper.

### 1.2 How to read the table

- **Raw, the genes do look "OXPHOS, not MYC".** That is where the sentence
  came from.
- **MYC's near-zero hides two opposing parts.** The reference MYC estimator
  largely tracks proliferation (rho 0.61-0.70), and proliferation goes with
  *lower* `BBC3`.
  - With proliferation partialled out, MYC-`BBC3` becomes **+0.16 / +0.11**
    (the provenance table, E47 declaration 4a.2).
  - On the estimator least tied to proliferation (`MYC_UP.V1_UP`), with each
    held fixed, the SCAN-B signature tracks MYC at least as strongly as
    OXPHOS: **+0.33 against +0.27** on level, **+0.37 against +0.23** on
    MitoPPS (E47, `m1`).
- **OXPHOS is only modestly tied to proliferation, and its association barely
  moves** (1.3).

**This explanation is for the author.** The proliferation coefficient behind
it is RETRIEVED, not declared (E47 result 3.1), so no manuscript sentence may
rest on it.

### 1.3 The models (E47, verified)

The model is `readout ~ OX + MYC + spine (+ PROLIF)`, OLS, on rank-normalised
scores. The 44 specifications are:
- 4 MYC estimators;
- 3 rungs: `m0` none, `m1` spine, `m3` spine plus proliferation;
- 2 OXPHOS measures;
- 2 cohorts.

That makes 48, minus the 4 void HALLMARK `m3` cells.

**OXPHOS with MYC held fixed: the robust half.**
- **`BBC3`, `BCL2L1` and `BAD` are positive, and `MCL1` negative, with the CI
  excluding zero in all 44.** `BID` is positive in 40, `BCL2` negative in 40.
- **On `M_a`:**
  - `BCL2L1` is +0.26 to +0.44 at `m1` and +0.28 to +0.44 at `m3`;
  - `BBC3` is +0.10 to +0.36 at `m1` and +0.11 to +0.34 at `m3`.

  **Adding proliferation moves no value by more than 0.03.**

**MYC with OXPHOS held fixed: the unresolved half.**

| gene | range over the 44 | CI above zero | CI below zero |
|---|---|---|---|
| `BBC3` | -0.45 to +0.30 | 21 | 16 |
| `BCL2L1` | -0.37 to +0.28 | 9 | 22 |
| `BID` | +0.06 to +0.55 | **44** | 0 |

The sign follows the estimator's proliferation entanglement and the
proliferation adjustment (E47 result 3.1, 3.2). `BID` is the one gene
positive everywhere.

## 2. The signature's provenance, and the two OXPHOS measures

### 2.1 The signature

**S6 is `BAD`, `BCL2L1`, `BIK`, `BBC3` and `BID` (+1), and `MCL1` (-1).**
- **How the genes were chosen.** From the author's twelve BCL-2-family genes
  (E10), by their own correlation with OXPHOS in these same two cohorts:
  GSVA, partial on proliferation, mean of TCGA and SCAN-B.
- **No selection rule was written down at the time.** The exact rule that
  reproduces the choice is a reconstruction (E47 4a.7): MitoCarta members
  with `|mean rho| >= 0.20`, each signed by its direction. `BID` qualifies by
  0.010.

**Three consequences:**
1. "The signature tracks OXPHOS" is guaranteed.
2. "But not MYC" compared a guaranteed number with an unconstrained one.
3. The informative test, MYC with OXPHOS held fixed, was **MIXED**: above
   zero in 4 of 8 cells, below in 1.

**The new text needs no signature.** Its claim rests on the two genes the
mouse named:
- **PUMA**, the trigger that Pgc-1alpha induces (Fig. 2d-g);
- **BCL-xL**, the protection that sustains OXPHOS-high tumours (Fig. 3h,i;
  4a).

### 2.2 The two OXPHOS measures

- **Level (GSVA):** how high the nuclear-encoded OXPHOS subunits are,
  relative to the other tumours in the cohort.
- **Priority (MitoPPS):** how strongly OXPHOS is prioritised relative to the
  tumour's other mitochondrial pathways. It is blind to level by design
  (CLAUDE.md trap 6).

How they relate:
- **They agree closely**: rho 0.88 / 0.82 (AD HOC).
- **Priority overlaps less with MYC activity**: rho 0.30 / 0.31 against
  0.39 / 0.43 for level (`g3`).

**Report both, every time.**
- Level is the only measure the clinical cohorts have; they are arrays.
- MitoPPS is the measure used on the mouse gland (Fig. 1i).

Holding on both is what links the mouse figure to the clinical one. **They
are two views of one quantity, not two replications.**

## 3. The mouse-human gap comes from the wording

| mouse | human | |
|---|---|---|
| Pgc-1alpha raised PUMA in iMMECs independently of Myc (Fig. 2f,g) | OXPHOS tracks `BBC3` at fixed MYC activity | consistent |
| High OXPHOS was sustained only with Bcl-xL, PUMA raised too (Fig. 3h,i; 4a) | `BCL2L1` tracks OXPHOS more steeply than `BBC3` | consistent, and the strongest bridge |
| Death required Myc (Fig. 2f-i) | cannot be tested | one Discussion sentence |

**The Doc already places the Myc requirement at death, not at PUMA:**
"induction of cell death by this trigger was determined by the Myc status
rather than by PUMA upregulation alone." Tumour transcripts measure the
regulators, not death, and established tumours are the survivors.

### 3.1 The one real difference, which helps the paper if stated plainly

**The human coupling is not confined to MYC-high tumours.** Within each MYC
half, OXPHOS-high tumours exceed OXPHOS-low ones by about the same amount, in
SD:

| gene, MYC half | TCGA level | TCGA priority | SCAN-B level | SCAN-B priority |
|---|---|---|---|---|
| `BCL2L1`, MYC-low | +0.75 | +0.53 | +0.46 | +0.39 |
| `BCL2L1`, MYC-high | +0.59 | +0.51 | +0.39 | +0.30 |
| `BBC3`, MYC-low | +0.53 | +0.39 | +0.19 | +0.20 |
| `BBC3`, MYC-high | +0.54 | +0.42 | +0.14 | +0.11 |

These are differences of E46's verified adjusted quadrant means, `Q2 - Q1` and
`Q4 - Q3`. They are descriptive and untested.

- **This is the Doc's "broader association".** Ref. 26 supports it:
  OXPHOS-high, BCL-xL-dependent persisters with suppressed MYC.
- **E33 already forbids implying confinement** to MYC-coupled tumours.
- **A reviewer may set it against Fig. 1i's genotype x OXPHOS term.** The
  answer: Fig. 1i compares transgenic with wild-type glands across
  development, which no tumour cohort can reproduce. E47 fitted main effects
  only; a formal interaction would need its own declaration.

### 3.2 "Paralleling the lower *Bbc3*/*Bcl2l1* ratio in the adult mouse gland" links the human result to the wrong path

- **In the adult gland the ratio fell *with* OXPHOS deprioritisation**: the
  first path (Fig. 1h,i).
- **The human result is about OXPHOS-*high* tumours**: the second path.
- **The human estimates compare strength of association in SD units.** They
  are not a ratio.

**Cite Fig. 3h,i instead.**

## 4. Text to paste into the Doc

**Abstract:**

> Extending this relationship beyond the experimental Myc models, OXPHOS
> tracked *BBC3* and, more steeply, *BCL2L1*, encoding BCL-xL, across 4,081
> human breast tumours, even after adjustment for subtype and MYC activity.

**Results, paragraph 1** (replaces "Thus, we next asked..." through "...four
control mitochondrial pathways"):

> We next asked whether the coupling of OXPHOS to PUMA and BCL-xL observed in
> the Myc models extends to established human breast tumours. In TCGA-BRCA
> (n = 938) and SCAN-B (n = 3,143), we scored OXPHOS in two ways: as the
> expression level of nuclear-encoded OXPHOS subunits relative to other
> tumours (GSVA), and as the priority of OXPHOS within each tumour's
> mitochondrial programme (MitoPPS), the measure applied to the mouse gland
> (Fig. 1i). OXPHOS correlated moderately with MYC activity (Spearman
> ρ = 0.30–0.43) and with proliferation (ρ = 0.23–0.42). We therefore
> modelled *BBC3* and *BCL2L1* expression on OXPHOS and MYC activity together,
> adjusting for PAM50 subtype and, in TCGA, tumour purity and leukocyte
> fraction. In both cohorts and with both measures, tumours with higher OXPHOS
> expressed more *BBC3* and, more steeply, more *BCL2L1* (0.10–0.36 and
> 0.26–0.44 standard deviations per standard deviation of OXPHOS; Fig. 4g).
> Both associations were essentially unchanged by additional adjustment for
> proliferation (Extended Data Fig. Xa). The association of MYC activity with
> these genes at fixed OXPHOS was not consistent: its direction depended on
> the MYC activity score used and on adjustment for proliferation, with which
> MYC activity scores were correlated (ρ = 0.41–0.84; Extended Data Fig. Xb).
> Across the twelve core BCL-2-family genes, *BAD* and *BID* were also higher
> in OXPHOS-high tumours, and *MCL1* and *BCL2* lower, in both cohorts and with
> both measures; the other six were inconsistent (Extended Data Fig. Xa). The
> anti-apoptotic profile associated with OXPHOS was therefore selective:
> *BCL2L1* rose while *MCL1* and *BCL2* fell.

**Results, paragraph 2:**

> *BCL2L1* rose more steeply with OXPHOS than *BBC3* in both cohorts and with
> both measures (difference, 0.08–0.19 standard deviations; Fig. 4g).
> OXPHOS-high human tumours thus resembled the mouse tumours in which Bcl-xL
> co-expression allowed high OXPHOS and raised PUMA to be sustained
> (Fig. 3h,i). This configuration was present in tumours with low as well as
> high MYC activity (Extended Data Fig. Xc), as in quiescent,
> doxorubicin-persistent triple-negative breast cancer cells, in which an
> OXPHOS-high state accompanies BCL-xL dependence despite suppressed MYC
> activity[26].

The co-amplification sentences go only if decision 5 admits them (section 6,
4h).

**Clinical paragraph, two edits:**

> Because the association of OXPHOS with *BBC3* and *BCL2L1* held after
> adjustment for MYC activity and was not explained by proliferation, we next
> asked whether its association with chemotherapy response is also distinct
> from that of proliferation.

> Additional adjustment for MYC activity, and then for the mean expression of
> *BCL2L1* and *MCL1*, did not materially alter the OXPHOS–response
> association (pooled odds ratios 0.87 and 0.88)...

**Conclusion paragraph:**

> In human breast tumours, higher OXPHOS was associated with higher *BBC3* and,
> more steeply, *BCL2L1* expression after adjustment for MYC activity, and with
> chemotherapy response in the direction opposite to proliferation, extending
> the OXPHOS–apoptosis relationship beyond the Myc-driven experimental setting
> and revealing clinical information distinct from proliferation.

**Discussion, new:**

> The human data cannot test whether OXPHOS-induced death requires MYC, as it
> does in the mouse: established tumours are the survivors of that constraint,
> and bulk transcriptomes report the regulators of apoptosis, not apoptosis
> itself. They do show the configuration the mouse experiments predict for
> tumours that maintain high OXPHOS: higher *BCL2L1* alongside higher *BBC3*
> expression.

**Fig. 4g legend** (if 4g is E47 Display C, section 6):

> **g**, Association of OXPHOS with *BCL2L1* (BCL-xL) and *BBC3* (PUMA)
> expression at fixed MYC activity in TCGA-BRCA (n = 938) and SCAN-B
> (n = 3,143), with OXPHOS scored as subunit expression level (GSVA) or
> priority (MitoPPS). Left, coefficients from linear models of each gene on
> OXPHOS and MYC activity, adjusted for PAM50 subtype and, in TCGA, tumour
> purity and leukocyte fraction; variables rank-normalised within cohort, so
> coefficients are in SD per SD. Right, *BCL2L1* minus *BBC3*, 95% intervals
> from 2,000 paired bootstrap resamples. Bars, 95% CI.

**Words avoided on purpose:**
- **"independently of MYC"**, which reads as a contrast with MYC and is not
  licensed (E47 result 9.1);
- **"at any level of MYC activity"**, which implies no interaction, and none
  was fitted;
- **"MYC expression"**, which is trap 4.

## 5. What comes out of the main text

| sentence | why |
|---|---|
| "Unexpectedly... rather than MYC activity" | MYC's association is unresolved, not absent |
| "the two OXPHOS-high groups showed comparable..." | true in 1 of 4 blocks. On MitoPPS, MYC-high exceeds MYC-low among OXPHOS-high tumours by 0.19 and 0.23 SD on the signature (E47 D1) |
| "Reciprocal adjustment... lost after adjustment for OXPHOS" | contradicted by E47; withdraw |
| the mitochondrial-ribosome control | the loss it was meant to qualify never happened |
| "invariant after... purity, and immune infiltration" | SCAN-B has neither; the new text states the adjustment instead |
| "restricted to apoptosis among four control pathways" | that result (E14) is about a different 44-gene set, and E38 found it was not apoptosis-specific in normal tissue |
| the ten-gene list "reproduced in both cohorts" | holds for six genes. `BIK`, `BMF`, `BCL2L11` and `PMAIP1` each fail in at least one cohort or measure (E47 Display A) |
| "paralleling the... adult mouse gland" | wrong path (3.2) |
| "anti-apoptotic composite score or *BCL2L1* expression" | E40 has no `BCL2L1`-only model; its `m3` adds the mean of `BCL2L1` and `MCL1` |
| "MYC expression", "4,302" | should be MYC activity, and 4,081 tumours analysed |
| the signature | the claim now rests on two genes the mouse named |

## 6. Figure plan

| panel | content | source | state |
|---|---|---|---|
| **4g**, main | the two genes at fixed MYC, plus their difference | E47 Display C | built and verified. Relabel the jargon ("beta_OX\|MYC", "m1, M_a") |
| 4h | co-amplification | validation study | blocked on decision 5. If admitted: Extended Data, with BOTH pre-specified ORs, 1.77 (1.26-2.47) and 1.36 (0.98-1.88). The "MCL1 restricted to basal" clause is not verified from this repo |
| 4i, 4j | pCR and METABRIC forests | E45 | unchanged |
| ED Xa | all twelve genes, both coefficients, with and without proliferation | E47 Display A | built; the author's decision of 2026-10-23 |
| ED Xb | MYC's coefficient on four MYC estimators | E47 Display B | built |
| ED Xc | quadrant heatmap | E46 | built. Add the median-split legend sentence (handoff 2026-10-23, 5.I) |
| Supp. table | all 44 specifications | `outputs/tables/E47_coefficients_OX_MYC.csv` | built |

**This reverses the 2026-10-23 arrangement, which put E46 in the main figure
as 4g. It is the author's call.** Why Display C:
- **It shows exactly the Results sentence.**
- **It matches the format of 4i and 4j**, so the human half of Fig. 4 becomes
  three effect-size panels.
- **It needs no signature and no caveat.**

**Why E46 moves to Extended Data.** It invites a MYC comparison it cannot
support: on GSVA its MYC step is mostly leftover OXPHOS (E47 D1). In Extended
Data it carries the MYC-low point (3.1) well.

**If 4h goes, 4i and 4j become 4h and 4i.** No filename changes: files are
named by content.

## 7. How it fits the paper, and how clear it becomes

- **One human claim**, worded the same way in the abstract, the Results, the
  figure and the conclusion. Each part maps to a mouse panel (2f,g; 3h,i;
  4a): the abstract's "second path", seen in patients.
- **The human section stops trying to decide what MYC does.** That is the
  mouse's job, and it does it experimentally. The human data extend the
  OXPHOS coupling and the protected state to patients.
- **The clinical section follows** with a one-clause bridge and is otherwise
  unchanged.
- **What is lost** is the "unexpected" hook. **What is gained** is that every
  sentence survives a reader checking it against the Extended Data, and
  nothing rests on a circular score.

**Soft spots, for the legend or Methods only:**
- TCGA's difference margin is thin, with lower bounds of 0.014 and 0.015;
- `BBC3`'s association is weak in SCAN-B (0.10 to 0.14);
- purity and immune adjustment exist in TCGA only;
- everything is post hoc and cross-sectional, and these are transcripts.

## 8. Decisions for the author

1. **4g:** E47 Display C (recommended) or E46.
2. **The signature.** Either:
   - **drop it from the paper entirely (recommended).** That needs a
     display-only `E48` to re-render Displays A and B, and E46, without the
     score rows;
   - **or keep it in Extended Data** with the E47 4a.7 Methods paragraph,
     construction sentence included.
3. **Decision 5**, the validation study's place in the paper, which decides
   4h.
4. **Keep the trigger-arm hypothesis out of the Discussion.**
   - **What the clean estimator shows.** On the estimator least tied to
     proliferation, `BCL2L1` also rises with MYC at fixed OXPHOS:
     - in 3 of 4 blocks, at both `m1` and `m3` (+0.08 to +0.28 where the CI
       clears zero);
     - the guardian arm in 4 of 4;
     - `MCL1` falls in 4 of 4.
   - **So the hypothesis's "not `BCL2L1`" half holds on `M_a` only.** E47
     result section 5 reports the arms on `M_a` only, so the note does not
     show this. CLAUDE.md's summary line on the arms is likewise
     `M_a`-specific.
   - **Recommended:** a dated pointer in E47 result section 10. **Not yet
     written: the author's call.**

**Against the 2026-10-23 handoff:**
- **Section 5:**
  - **settled by section 4's text:** rows A, B, C, D, F, G and H;
  - **row I** is met by the new 4g legend, plus E46's sentence in Extended
    Data;
  - **row E** stays on decision 5.
- **Section 7:**
  - **item 2** is answered (item 4 above);
  - **items 3 and 4** are proposed (section 6);
  - **item 5**, the push, was done at `d17e1a1`;
  - **items 1, 6 and 7** are unchanged.

## 9. Where every number comes from

| number | source |
|---|---|
| rho OXPHOS-MYC, OXPHOS-proliferation, MYC-proliferation (all estimators) | `e47$gates$g3`, saved |
| OXPHOS and MYC coefficients, the 44 specifications | `e47$results`, saved |
| `BCL2L1` minus `BBC3` and its CI | `e47$bcl2l1_bbc3`, saved |
| MYC step among OXPHOS-high tumours on the signature (0.19, 0.23) | `e47$d1_leakage`, saved |
| the within-half OXPHOS steps (3.1) | `e46$heatmap_cells`, saved; differences of adjusted means |
| MYC-`BBC3` partial on proliferation, +0.16 / +0.11 | E47 declaration 4a.2 (the provenance table) |
| pooled pCR odds ratios 0.85 / 0.87 / 0.88 | E40, `exp()` of the pooled log-ORs -0.1616 / -0.1345 / -0.1283 |
| co-amplification ORs | validation study G2, via `2026-10-22_six_gene_composite_provenance.md` 4.2 |
| **AD HOC**: the gene rows of 1.1, rho(level, priority) 0.88 / 0.82 | computed read-only from `e47$frame` on 2026-11-02, code below |

**The AD HOC numbers fit nothing and change no saved object.** They reproduce
the documented ones where they overlap. Partial on proliferation:
- MYC-`BBC3` is +0.174 / +0.105, against 4a.2's +0.163 / +0.110;
- OXPHOS-`BCL2L1` is +0.427 / +0.348, against +0.426 / +0.350.

The small gaps are the analysis sets. **For the manuscript, any of them must
come from a scripted display first.**

```r
suppressMessages(library(dplyr))
x  <- readRDS("results/e47_mutual_adjustment.rds")
fr <- x$frame %>% filter(in_set)                 # TCGA 938, SCAN-B 3,143
v  <- c("OX_gsva", "OX_mitopps", "MYC_M_a", "MYC_clean", "PROLIF",
        "S6", "BBC3", "BCL2L1")
for (co in c("TCGA", "SCAN-B")) {
  d <- fr %>% filter(cohort == co) %>% select(all_of(v))
  print(round(cor(d, method = "spearman", use = "pairwise.complete.obs"), 2))
}

# partial Spearman on proliferation: correlate the residuals of the ranks
pc <- function(a, b, z) {
  rz <- rank(z)
  cor(resid(lm(rank(a) ~ rz)), resid(lm(rank(b) ~ rz)))
}

# the 44 specifications, per gene
x$results %>%
  filter(readout %in% c("BBC3", "BCL2L1"), !void) %>%
  group_by(readout) %>%
  summarise(n = n(), ox_above = sum(ox_lo > 0),
            myc_above = sum(myc_lo > 0), myc_below = sum(myc_hi < 0))

# the within-half OXPHOS steps (3.1)
readRDS("results/e46_fig4_quadrant_panels.rds")$heatmap_cells %>%
  filter(gene %in% c("BBC3", "BCL2L1")) %>%
  select(instrument, cohort, gene, q, est) %>%
  tidyr::pivot_wider(names_from = q, values_from = est) %>%
  mutate(step_MYC_low = Q2 - Q1, step_MYC_high = Q4 - Q3)
```

**Mouse figure references are as of the 2026-11-02 read of the Doc:**
- Fig. 1h,i;
- Fig. 2d-i;
- Fig. 3h,i;
- Fig. 4a;
- the line-77 quotation.

**If the mouse panels are re-lettered, update section 3 and the drop-in
text.**
