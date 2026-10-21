---
date: 2026-10-08
purpose: >
  Section 4 narrative, version 5 — the submission draft. Supersedes
  2026-10-08_section4_narrative_v4_final_draft.md.
status: >
  FINAL DRAFT for initial submission. Word limits relaxed; completeness of the
  argument prioritised over compression.
changes-from-v4:
  - Paragraph 4's distant-outcome passage rewritten against E43 and E44. The
    METABRIC coefficient is no longer reported as a result, because
    E39/E40/E41 section 7.2 forbids a distant-outcome claim of any kind.
  - The proliferation coefficient (E43) enters as the comparator and as an
    unplanned positive control.
  - The respiration-proliferation contrast (E44) enters as a tested quantity,
    replacing the interval-overlap statement v4 relied on.
  - Fig. 4h and 4i now carry both terms, since the contrast is the result.
excludes: >
  X01 and everything in it, which is manuscript-excluded by its own posture.
recovered: >
  Reconstructed 2026-10-21 from the chat "Section 4 handoff and next steps"
  (2026-10-08). Goes in docs/ as 2026-10-08_section4_narrative_v5.md.
---

# Section 4 — narrative v5, submission draft

## 1. Standing rules applied throughout

**Figure legends are descriptive only.** A legend states what the panel shows —
axes, markers, groups, n, the adjustment set, what the error bars represent —
and nothing else. No result, no effect size offered as a finding, no
qualification, no interpretation. Everything of that kind lives in the main
text, which is why the paragraphs below are longer than a Letter's usual share.

**Three prohibitions from the arm's own result notes govern paragraph 4.**

1. **No distant-outcome claim of any kind for respiration.** Nothing met the
   reading rule in any pre-specified set, in either outcome cohort. The METABRIC
   ER-positive cell is one boundary cell of thirty, excluding zero by 3.7 × 10⁻⁴,
   unreproduced by its own declared sensitivity, with the adjacent rungs
   including the null.
2. **A significant contrast does not make either coefficient a finding.** E44
   licenses a statement that the two coefficients differ, not that either is
   established.
3. **The two endpoints are different quantities.** No sentence may say the pCR
   finding was replicated, confirmed or extended in METABRIC.

---

## 2. The narrative

### Paragraph 1 — the mouse in vivo consequence

*Unchanged. See `2026-09-24_section4_narrative_and_figure_plan.md`.*

**One flag carried forward.** The LbNOX sentences assert that the growth
advantage requires ATP-coupled respiration, and the panel supporting it cannot
distinguish an absent advantage from an advantage offset by an apoptotic cost,
because LbNOX was never tested with Bcl-xL. The two source documents also
disagree on whether mtLbNOX affects apoptosis. Unresolved, and it affects the
Discussion's mechanism sentence.

### Paragraph 2 — the configuration is general, and does not follow MYC

> We next asked whether the relationship between respiration and the apoptotic
> machinery is confined to the window the mouse defines or holds across
> established tumours. In TCGA-BRCA (n = 1,095) and SCAN-B (n = 3,207) we
> divided tumours at the median of MYC activity and of OXPHOS status into four
> groups. The machinery separated by respiratory state and not by MYC, with the
> two OXPHOS-high groups carrying the configuration equally and the two
> OXPHOS-low groups equally lacking it (Fig. 4g). Sixteen of eighteen readable
> comparisons ordered by respiratory state, with MYC contributing four to seven
> times less and the same amount at either respiratory level. These estimates
> are adjusted for intrinsic subtype, tumour purity, leukocyte fraction and
> proliferation. Conditioning each axis on the other gave the same answer, since
> removing MYC left the OXPHOS ordering intact whilst removing OXPHOS abolished
> the MYC ordering. The configuration that MYC exposes in the mouse gland is
> therefore general among tumours that reach a high respiratory state.

### Paragraph 3 — what the configuration is, and the lesion that accompanies it

> The configuration itself replicated across both cohorts, with BAD, BCL2L1,
> BIK, BBC3 and BID higher in OXPHOS-high tumours and MCL1, BMF, BCL2L11, BCL2
> and PMAIP1 lower (Fig. 4g). A tumour retaining a high respiratory state cannot
> dispense with the mitochondrial half of this machinery, and the cytosolic half
> is the only part it can lose without cost. The compartment difference is an
> ordering rather than a difference in magnitude. Neither cytosolic median was
> individually distinguishable from its own null, and the argument rests on
> their being the extreme of thirty comparator pathways in both cohorts.
> Consistent with the buffering the mouse requires, tumours carrying MYC
> amplification showed an excess of both BCL2L1 and MCL1 amplification
> (Extended Data Fig. 9).

### Paragraph 4 — separating the two consequences, and the close

> Respiration carries a second consequence alongside the apoptotic one, since it
> supports the biosynthesis on which proliferation depends. Whether the two
> consequences are separable, and what each predicts clinically, requires
> cohorts in which treatment outcome is recorded. We therefore examined 1,935
> patients across three neoadjuvant cohorts, in which pathological complete
> response to chemotherapy was scored, together with a long-term outcome cohort.
> Respiration and proliferation proved separable in all four, correlating at
> Spearman 0.36 to 0.49, so a model can distinguish their contributions.
> Respiration alone bore no consistent relationship to pathological complete
> response, the three cohort estimates disagreeing in sign and heterogeneity
> reaching 49 per cent. Adjusting for proliferation brought them into agreement,
> giving a pooled log odds of −0.162 per standard deviation of respiration, an
> odds ratio of 0.85 (Fig. 4h). Heterogeneity fell to zero, all three cohorts
> agreed in sign, and the relationship held in every subtype with enough events
> to fit. Agreement in sign across all three cohorts was fixed in advance as the
> condition for reading the pooled estimate as a signal. The relationship is
> invisible without adjustment, because the two consequences oppose one another
> and cancel. Proliferation carried in the same models ran the other way,
> predicting better response at a pooled log odds of +0.542. The difference
> between the two coefficients is itself estimable, and on this endpoint it
> pools at 0.66 log odds per standard deviation with a confidence interval of
> 0.44 to 0.87. The cohorts disagree about the size of the proliferation
> coefficient and agree about the size of the gap, which pools without
> heterogeneity. Neither the guardian composite nor BCL2L1 alone modified the
> respiratory relationship, contrary to the prediction the pre-registered
> analysis had fixed in advance. We next examined long-term survival rather than
> treatment response, in the METABRIC cohort of 1,981 patients. This endpoint
> records outcome under mixed and largely historical care, so neither result
> replicates the other. Proliferation again behaved canonically, predicting
> shorter breast cancer survival in receptor positive disease at a hazard ratio
> of 1.23. The respiratory coefficient ran the other way at 0.90, and the
> difference between the two again excludes zero, in the primary stratum and in
> both declared comparators. That respiratory estimate does not meet the reading
> rule declared for this analysis, since it lies at the boundary, its declared
> sensitivity does not reproduce it, and the adjacent rungs include the null. We
> therefore report no association between respiration and distant outcome, and
> read the contrast with proliferation rather than the coefficient itself.
> Respiration therefore inverts, on both endpoints, a pair of associations that
> these models recover canonically in proliferation. Crucially, a high
> respiratory state does not establish dependence on respiration, and the two
> are conflated in current proposals to target oxidative phosphorylation.
> Inhibition would withdraw the apoptotic liability together with the
> biosynthetic benefit it targets, whilst the cytotoxic direction on a MYC
> background is the opposite one.

**Metrics.** 22 sentences, 466 words, mean 21.2, median 21, range 15 to 33,
nothing under 12 or over 38, no serial comma. Share under 20 words is 32 per
cent against the corpus band of 35 to 45, which is expected where nearly every
sentence carries a result.

**Three sentences are optional if the paragraph runs long.** The heterogeneity
sentence ("The cohorts disagree about the size…"), the pre-registered negative
("Neither the guardian composite…"), and the sign-agreement sentence
("Agreement in sign across all three cohorts…") can each move to Methods or an
Extended Data note without damaging the argument. Cutting all three brings the
paragraph to roughly 400 words.

### Discussion — one short paragraph

> The two mitochondrial biclusters previously identified in human breast tumours
> separate the two axes examined here, one indexing the MYC component and the
> other the respiratory component on which the apoptotic configuration sits. The
> upper fork of the respiratory bicluster has been reported to carry worse
> outcome, and that axis correlates more strongly with proliferation than with
> the respiratory score used here, which may reconcile the two directions.
> Respiration therefore supports two opposing consequences in an established
> tumour, and clinical associations that do not separate them report their sum.
> A high respiratory state marks a tumour that has already neutralised the
> accompanying death signal, so inhibiting respiration would remove that signal
> rather than exploit it. The cytotoxic direction our mouse data identify is the
> opposite of inhibition, namely raising mitochondrial NADH oxidation, which was
> lethal only on a MYC background. Guardian antagonism remains the alternative
> target.

**Two notes.** The therapeutic sentence is here because it has been removed from
the abstract, and it uses "remove" rather than "relieve" — relieving a burden
and relieving a repression point in opposite directions. The mechanism sentence
no longer says "without proton pumping", because that formulation leans on the
mtLbNOX result the two source documents dispute.

**Guardrails, unchanged.** MYC-selectivity is the claim, not that uncoupling
kills cells. Do not attribute the published prognostic result to buffering.

---

## 3. Main figure — the human arm

Panel letters assume the mouse arm occupies 4a to 4f.

### Fig. 4g — the configuration by quadrant

**Layout.** A heatmap with a point-and-interval strip beneath it, sharing one
column axis.

- **Columns**, eight, grouped by cohort then quadrant: TCGA then SCAN-B, each
  Q1 MYC-low/OXPHOS-low, Q3 MYC-high/OXPHOS-low, Q2 MYC-low/OXPHOS-high,
  Q4 MYC-high/OXPHOS-high. **Order them so the two OXPHOS-high columns are
  adjacent**, with a thin rule between the respiratory levels. The colour change
  should fall at that rule and not at the MYC boundary.
- **Rows**, twelve, the configuration genes in two blocks separated by a gap,
  MitoCarta members first and cytosolic effectors second. Label the blocks.
- **Colour**, standardised difference from the cohort mean, diverging, white at
  zero, one shared colour bar.
- **Strip beneath**, the six-gene signed composite per column, point with 95 per
  cent interval.

**Legend, descriptive only.** Cohorts and n. The four quadrants and the median
splits producing them. That rows are the twelve configuration genes and the
strip the six-gene signed composite, which is a different object from the rows.
The covariates in the model. What the colour scale and error bars represent.

### Fig. 4h — the decomposition on chemotherapy response

**CHANGED FROM v4: the panel now carries both terms.** The contrast is the
result, so the opposition must read off the panel rather than from the text.

**Layout.** A paired forest, two markers per row.

- **Rows**, four: GSE194040, GSE164458, GSE25066, pooled. The pooled row is
  separated by a rule and set in bold.
- **Within each row, two terms**: respiration and proliferation, from the same
  model, as two markers distinguished by shape or fill, each with its interval.
- **The unadjusted respiration estimate** appears as a third, open marker on each
  row, joined to the adjusted one by a thin connector, so the m0-to-m1 movement
  stays visible. If three markers per row is too dense, move the unadjusted
  estimates to ED Fig. 8 and say so in the legend.
- **A bracket** spanning the two adjusted markers on the pooled row, labelled
  with the contrast and its interval.
- **x-axis**, log odds ratio per within-cohort standard deviation, zero marked.
- **Right-hand annotation**, n and pCR events per cohort.

**Legend, descriptive only.** The three cohorts, n and events. That markers are
log odds ratios per within-cohort standard deviation of each score. Which
markers denote which term and which adjustment. The model contents. The
meta-analysis method. That the bracket shows the within-model difference between
the two terms, estimated with the covariance from the model's variance matrix.
What the bars represent.

### Fig. 4i — the same decomposition on distant outcome

**Layout.** The same grammar as 4h, given less width, so the reader re-learns
nothing.

- **Rows**, three: METABRIC breast-cancer-specific survival in receptor positive
  disease (the primary), the same in receptor negative disease, and GSE25066
  distant relapse-free survival.
- **Within each row, two terms**, as in 4h, with the same shapes.
- **A bracket** on the receptor positive row, labelled with the contrast and its
  interval.
- **x-axis**, log hazard ratio per within-cohort standard deviation, 1 marked.
  **Not shared with 4h** — odds ratios and hazard ratios are different
  quantities and a shared axis would imply otherwise.
- **Right-hand annotation**, n and events per row.

**Legend, descriptive only.** Cohorts, endpoints, n and events. The receptor
status definition used for the strata. The model contents. That the bracket
shows the within-model difference, estimated with the covariance. What the bars
represent.

### If the display budget forces two human panels

Keep **4g** and merge **4h and 4i into one panel with two sub-plots**, the
response forest at roughly twice the width of the outcome forest. Do not drop 4h
to keep 4i.

---

## 4. Extended Data — four figures

Legends descriptive only, as above.

**ED Fig. 7 — the quadrant result in full.** The conditioning ladder with
matched nulls, both instruments, both cohorts. The eighteen-cell table with the
two exceptions. The thirty comparator pathways behind the compartment ordering.
The renormalisation check. The burden strata. The bicluster partials.

**ED Fig. 8 — the decomposition in full.** (a) The identifiability gate: the
three pairwise correlations in all four cohorts, with the declared bands marked.
(b) The complete ladder, every rung, both endpoints, every cohort and stratum,
with intervals and the variance inflation factors. (c) The subtype scope check,
per-subtype forest with the pooled estimate as a reference line. (d) The
contrast in every fit, with the naive and covariance-corrected standard errors
side by side. (e) Proportional hazards diagnostics and the declared
sensitivities.

**ED Fig. 9 — the pre-registered analysis, and the lesion.** The H4 forest
across all three cohorts with the meta-analysis. Alongside it, BCL2L1 and MCL1
co-amplification with MYC.

**ED Fig. 10 — instruments, coverage and robustness.** The twelve configuration
genes with intervals in both cohorts. The METABRIC coverage audit and the
restriction check. Purity, leukocyte and subtype sensitivities. The MYC
estimator concordance.

---

## 5. Methods — three statements that must appear

1. **The contrast was estimated as a within-model linear combination**, with the
   covariance between the two coefficients taken from each fit's variance
   matrix. The naive calculation ignoring that covariance is anti-conservative
   here, by 14 to 23 per cent on the standard error in every fit.
2. **The pooled contrast is the meta-analysis of the per-cohort contrasts**, not
   the difference of two separately pooled estimates, which has no valid
   standard error because the two estimates come from overlapping patients.
3. **Proliferation was carried as a declared covariate, not as a control.** Its
   behaviour is reported as a post-hoc characterisation of a declared term.

---

## 6. Open items before submission

**The normal-breast comparison.** It licenses no claim, so it supports no
sentence and no panel. Either omit it from the figures and record in Methods
that the analysis was performed, that a positive control was declared in
advance, that the control fired and that no reading was taken; or give it one
Extended Data panel showing the dissociation beside the control that reproduced
it. I lean to the first. Author's decision.

**The LbNOX arm.** Paragraph 1's ATP-coupling claim rests on a panel that cannot
distinguish an absent advantage from an offset one, and the two source documents
disagree on whether mtLbNOX affects apoptosis. The minimum repair is to report
the observation and state that LbNOX was not tested with Bcl-xL. Cleaved
caspase-3 and Ki67 staining of the existing LbNOX tumours would settle it
without a new in vivo arm.

**The companion paper's figures.** The Discussion characterises the published
bicluster survival result, and that characterisation is author-sourced and has
never been verified against the PDF from inside any repository. Check panel
labels, strata, endpoints and hazard ratios before the sentence is final.

**BCL2L1 exceeding BBC3.** If the abstract or Results state that the guardian
tracks respiration more steeply than the trigger, file the paired contrast in
the notes first. It is currently recorded as unreported.

**The legend rule across the whole manuscript.** Applied here to Section 4 only.
The mouse paragraphs and their legends need the same pass.
