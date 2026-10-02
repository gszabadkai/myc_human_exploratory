---
date: 2026-10-02
status: DECLARATION - read rules fixed before any new fit. Commit this before
  the script is written.
posture: >
  EXPLORATORY and POST-HOC throughout. This is a descriptive decomposition of an
  exposure, not a fifth hypothesis. Plan section 2 forbids post-hoc hypotheses;
  nothing here carries a pass, a fail or a falsification criterion, and no result
  from it may be written as though it did.
relates-to:
  - 2026-08-31_block_F1_result_H4_not_supported.md (section 5, the seen numbers)
  - 2026-08-30_STATE_frozen_and_H4_buffer_declaration.md
  - 2026-08-28_D5_cohort_selection.md (power, and the "not powered to exclude"
    reporting rule)
  - 2026-08-28_D7_proliferation_covariate.md (PROLIF_DISJOINT, its definition and
    why it is specific to M-a)
  - 2026-08-29_escape_reading_declaration.md (section 6, the Block C main-effect
    discard and the rule governing retrieval)
  - 2026-09-02_e11_prolif_adjusted.md (the exposure-side proliferation result)
  - 2026-08-29_G3_result_forkscale_availability.md (section 5, the METABRIC file)
decides:
  - E39 runs on GSE25066 first, both endpoints, same patients.
  - METABRIC is declared IN, on design grounds, before any outcome is seen.
  - m3 (BUFFER_c) is retained as a rung despite its weak construct.
  - ER stratification is retained as descriptive, with the ER-positive stratum
    declared non-estimable in GSE25066 in advance.
next-action: write the script. Nothing is fitted until this note is committed.
---

# E39 - decomposing the respiratory axis against two endpoints

## 1. The question

Respiration in breast tumours is coupled to proliferation, and proliferation
carries opposite consequences for the two outcomes: proliferative tumours respond
better to neoadjuvant chemotherapy and relapse sooner. The clinical associations
of oxidative phosphorylation reported so far have not separated these
contributions.

This note declares a descriptive decomposition of the respiratory axis into what
is carried by proliferation and what is not, read against both endpoints, in
cohorts already on disk.

It is not a test of anything. It exists because the current draft contains a
sentence - "respiration did not itself predict response" - whose supporting
quantity was never fitted with a covariate set, an interval or a pooled estimate.

## 2. Disclosure - what has already been seen

The OXPHOS coefficient **at fixed MYC and BUFFER_c** was printed once and
survives only as prose in the Block F1 note section 5:

```
GSE194040  -0.074   p 0.32
GSE164458  -0.066   p 0.54
GSE25066   -0.120   p 0.37
```

No CIs, no meta-analysis, no subtype or treatment adjustment. It came from an
`if (FALSE)` sandbox at `scripts/13_outcome_models.R:745` which does not execute
in a normal run, so those three numbers have no reproducible provenance.

**Those three point estimates are therefore NOT blind.** Everything else declared
below is: every CI, every pooled estimate, every heterogeneity statistic, every
proliferation-adjusted estimate, and the whole DRFS side.

No reading rule in section 9 is conditioned on those three numbers.

### 2.1 Retrieval of a main effect is pre-authorised

The Block C main-effect discard was declared in advance twice - in
`2026-08-29_script09_build_spec.md:311-313` ("extracting the interaction
coefficient only") and in `block_c_models.rds$spec$coprimary` ("one
pre-registered coefficient per instrument") - and then restated in
`2026-08-29_escape_reading_declaration.md:197-201`, which anticipates retrieval:

> Those main effects are relevant to a constraint reading, and retrieving them is
> a one-line change. They must not be inspected before this note is committed,
> and when retrieved they inform E2's interpretation without redefining it.

So estimating a main effect here is a pre-authorised move in this arm, not a
reversal of a decision. What was never separately declared is any rationale
specific to excluding the OXPHOS main effect; it falls out of a one-coefficient
pre-registration.

## 3. Identifiability gate - run and read before any model

Compute, **per cohort and in all three neoadjuvant cohorts**, Spearman between
OXPHOS and PROLIF_DISJOINT, between OXPHOS and MYC, and between MYC and
PROLIF_DISJOINT. Report all three with n, in GSE25066, GSE194040 and GSE164458.
The band in the table below is read per cohort; a cohort may be identifiable
where another is not, and the consequence then applies to that cohort alone.

| rho(OXPHOS, PROLIF_DISJOINT) | consequence |
|---|---|
| >= 0.80 | the ladder is not identifiable in that cohort. Report m0 only there and interpret no rung beyond it. |
| 0.60 to 0.80 | report the ladder with variance inflation factors beside every coefficient, and state in the result note that the rungs are not cleanly separable. |
| < 0.60 | report the ladder as specified. |

This gate is descriptive and stops nothing else in the arm. It is a
**stop-and-check**: the gate output is read before any outcome model is fitted.

### 3.1 Gate result, recorded 2026-10-02

E39 ran the gate (`results/e39_prepare_and_gate.rds`; declaration at `536c29f`).

| cohort | n | rho(OX, PROLIF_DISJOINT) | band |
|---|---|---|---|
| GSE25066 | 508 | **0.358** | < 0.60, ladder as specified |
| GSE194040 | 988 | **0.381** | < 0.60, ladder as specified |
| GSE164458 | 482 | **0.493** | < 0.60, ladder as specified |

**The ladder is identifiable in all three cohorts.** Respiration and
proliferation share **13%, 15% and 24% of RANK variance** respectively (Spearman
rho squared; these are rank statistics, not Pearson r squared, and the phrase is
written that way deliberately). They are substantially separate axes in these
tumours. **This had not previously been measured in this arm and is reportable in
its own right.**

### 3.2 What the gate found that it was not looking for

`rho(MYC, PROLIF_DISJOINT)` is **0.780, 0.785 and 0.824**. GSE164458 is in the
`>= 0.80` band; the other two are in the `0.60-0.80` band.

**The gate table above bands only `rho(OXPHOS, PROLIF_DISJOINT)`.** The bands are
applied to this pair for description; no declared consequence attaches to them
here, and the consequence that does attach is specified in 5.3.

This is **D7's problem in a form D7 did not address.** D7 removed the nine genes
shared between the stripped Felsher estimator and the proliferation covariate,
which fixed overlapping *measurement*. It could not fix the underlying
*correlation*, because MYC drives proliferation and G1 measured that enrichment
at **8.6-fold**. `PROLIF_DISJOINT` is **disjoint by construction and
near-collinear by biology.**

`rho(OX, MYC)` is **0.358, 0.358 and 0.506**, so OX is only moderately correlated
with either, and **the OX coefficient should remain estimable at every rung.**
The consequence is for **attribution, not estimation**, and is specified in 5.3.

## 4. Exposure, scoring and scale

- **OXPHOS** is `OXPHOS subunits`, GSVA, **reused from the H4 frames rather than
  rescored**, so the exposure is byte-identical to the one H4 used. Any rescoring
  would make this a different analysis.
- **mitoPPS does not exist in these cohorts** (amendment A1). This is
  **SINGLE-INSTRUMENT throughout** and must not be written as though it had
  cleared the arm's two-instrument bar. Blocks C, B and G each required agreement
  between two instruments; this cannot and does not.
- **Proliferation** is `PROLIF_DISJOINT`, the 318-gene score already standard
  across this arm, scored into **ALL THREE** neoadjuvant cohorts from the
  on-disk expression matrices with the same GSVA parameters used elsewhere
  (`kcdf = "Gaussian"`, log-scale input). Recovery, reported in the result note
  as the exact figures:

  | cohort | PROLIF_DISJOINT recovered |
  |---|---|
  | GSE25066 | **293 of 318 (92.1%)** |
  | GSE194040 | **305 of 318 (95.9%)** |
  | GSE164458 | **311 of 318 (97.8%)** |

  `PROLIF_STD` recovers 300 of 327 (91.7%) in GSE25066.

  **Why all three, resolving the ambiguity against section 3.** `m1` is a rung
  of the three-cohort pCR ladder. Scoring proliferation in one cohort only would
  give a **three-cohort m0 against a one-cohort m1**, and the section 9
  meta-analysis would then be pooling at one rung and not at the next. **That is
  not a ladder.** Section 3's "per cohort" is the governing wording and section
  4's earlier "into GSE25066" is superseded by this bullet.
- **PROLIF_DISJOINT is defined against the stripped Felsher estimator (M-a)** per
  D7 and must never be used with M-b. M-a is the MYC estimator throughout here.
  Verified by construction rather than by label: the **9 genes** excluded from
  `PROLIF_STD` to form it - `CTPS1, DNMT1, HMGA1, NCL, NOP56, PRMT5, RANBP1,
  TFDP1, UCK2` - are all in M-a, and `PROLIF_DISJOINT` intersects M-a in **0**
  genes.

### 4.1 Coverage caveat, declared here

GSE25066 carries **57 of 89 OXPHOS subunits (64.0%)** - the coverage figure that
excluded it from H4's specificity arm. The exposure is measured on two thirds of
the gene set, on a GPL96 array, while every TCGA and SCAN-B number in the
manuscript uses the full set.

**No coefficient from GSE25066 may be compared numerically with a TCGA or SCAN-B
coefficient**, and no claim may rest on the GSE25066 exposure being the same
instrument as the one used elsewhere in the arm.

### 4.2 Scale

GSVA scores are cohort-relative. Every exposure is standardised to **within-cohort
SD** before fitting, so coefficients read per within-cohort SD and never as raw
values across cohorts. This is a condition for the meta-analysis in section 9 to
mean anything.

### 4.3 OX arrives already standardised, and is NOT re-standardised

`script 13` applies `.z()` when it builds the H4 frames, so `OX` arrives with
**sd = 1.000000** in all three cohorts (mean ~1e-16). **E39 does not
re-standardise it.** Only `PROLIF_DISJOINT` is standardised to within-cohort SD.

This satisfies both requirements at once: section 4's demand that the exposure be
**byte-identical** to the one H4 used, and section 4.2's demand that every
exposure read in **within-cohort SD units**. Re-standardising would be
arithmetically a no-op and would break the first for no gain in the second.

### 4.4 The GSVA call, fixed exactly

E39 matches **`scripts/13_outcome_models.R`**, not `E02`. The two differ, and
script 13 is the one that produced the `OX` and `MYC` in the H4 frames:

```r
M <- .impute(neo$cohorts[[cn]]$expr, cn)   # A5: drop genes >5% NA, impute at
                                            # the gene median. Called even where
                                            # it is a no-op, so provenance is
                                            # identical to script 13.
sets[[".PIN_A"]] <- rownames(M)[c(TRUE, FALSE)]
sets[[".PIN_B"]] <- rownames(M)[c(FALSE, TRUE)]
par <- GSVA::gsvaParam(exprData = M, geneSets = sets, kcdf = "Gaussian",
                       minSize = 3L, maxSize = Inf)
s   <- GSVA::gsva(par, verbose = FALSE)
```

log2 input; universe pinned with `.PIN_A` / `.PIN_B`; **a silently dropped set is
a hard stop**. `.impute()` is a verified no-op in GSE25066 (0 NA) but is called
anyway - it may not be one in the other two cohorts.

**It was not.** E39's run found that in **GSE194040 one gene exceeded the 5% NA
threshold and was dropped, and 2,773 values were imputed at the gene median.**
Had PROLIF been scored on the raw `$expr`, it would have come from a different
matrix than the `OX` beside it in that cohort. **Record this in the result
note.**

## 5. The ladder

Same spine at every rung, each rung adding one term:

```
m0:  endpoint ~ OX + subtype + treatment
m1:  + PROLIF_DISJOINT
m2:  + MYC
m3:  + BUFFER_c
```

Report the **OX coefficient with its CI at every rung, in every cohort**, plus the
pooled estimate and heterogeneity (Q, I2, both fixed and random) wherever more
than one cohort contributes.

What each step answers:

- **m0** is the quantity the draft sentence "respiration did not itself predict
  response" refers to. It has never been fitted.
- **m0 -> m1** is how much of the respiratory axis is proliferation.
- **m1 -> m2** is whether what remains is MYC. E34 predicts almost no movement,
  since the configuration tracks respiratory state and not MYC.
- **m2** is the estimand of the three seen numbers in section 2.
- **m2 -> m3** is whether the buffer accounts for it.

### 5.1 m3 is retained, with its construct weakness recorded

`BUFFER_c` is the mean of `z_MCL1` and `z_BCL2L1`. The Block F1 note section 3.1
records that the two limbs correlate at **0.019, 0.120 and -0.016** across the
three cohorts, so the composite averages two variables carrying different
information rather than two readings of one latent capacity.

m3 is retained anyway, because the movement from m2 to m3 is informative about
the exposure regardless of whether the composite names a real construct. **No
result from m3 may be described as being about "buffering capacity."** It is
about adjustment for two anti-apoptotic transcripts.

### 5.2 The spine is cohort-specific, by audit rather than by assumption

**Before fitting, each cohort's spine terms are audited for constancy and any
constant term is omitted explicitly, never left for R to drop silently.** Known
cases:

- `treatment` is **constant in GSE25066** (taxane-anthracycline, all 508).
- `subtype` is **constant in GSE164458** (BrighTNess is TNBC-only; 482 of 482).
  Asserted and reported rather than assumed.
- `treatment` carries **13 levels in GSE194040** and **3 in GSE164458**, both
  fitted.

`subtype` in GSE25066 is **collapsed to two levels**, `HRpos_HER2neg` against
`TNBC`, and the **5 HER2-positive patients are dropped, not pooled**
(`HRneg_HER2pos` n = 4, `HRpos_HER2pos` n = 1). A factor level of n = 1 is not
fittable at this event count.

After the drop, **`subtype2` and `er_primary` partition GSE25066 identically**,
278 / 187 on both. **This is an identity by construction, not a coincidence**:
`scripts/12_fetch_neoadjuvant_cohorts.R` builds `subtype` from `er_status_ihc`
and `her2_status` **only** - `pr_status_ihc` never enters - so once the
HER2-positive patients are dropped, `subtype2` is a relabelling of
`er_status_ihc`.

It is **asserted in the script anyway**, because the assertion guards the two
things that could still break it: a change in the upstream subtype definition,
and the NA path (`is.na(er_status_ihc) | is.na(her2_status)` maps to NA, so a
patient could in principle be callable on one and not the other). Where the
assertion holds, **the two are one variable**: `subtype2` is redundant in the
pooled model and constant inside each ER stratum, and is omitted there
explicitly.

### 5.3 MYC and proliferation are near-collinear: what m2 and m3 can be read as

At `rho(MYC, PROLIF_DISJOINT)` of **0.78 to 0.82**, MYC is largely a second
measurement of proliferation in these cohorts.

**The m1 -> m2 step therefore does not decompose by variable.** The question it
answers is restated here, before the fit:

- **Declared OUT**: *"whether what remains after proliferation is MYC"*. **That
  attribution is not available at this collinearity and must not be written.**
- **Declared IN**: whether the OX coefficient is **robust to the addition of a
  covariate strongly correlated with the one already in the model.**

The same applies to **m2 -> m3**, with `BUFFER_c` added to an already
near-collinear pair.

**Variance inflation factors for the OX term are reported at every rung, in every
cohort, alongside the coefficient.** Where the spine carries factors the
generalised form is used, and the measure reported is named in the output. **A
VIF above 5 for OX at any rung is reported in the result note as a caveat on that
rung's interpretability; it does not stop the rung being reported.**

**Nothing here licenses dropping MYC or PROLIF from the ladder.** Both stay, and
the limitation is recorded rather than engineered away.

## 6. Endpoints and analysis sets

| endpoint | cohorts | n | model |
|---|---|---|---|
| pCR | GSE194040 / GSE164458 / GSE25066 | 988 / 482 / 470 | logistic |
| DRFS | GSE25066 only | 470 primary | Cox |

### 6.1 The DRFS set is the pCR model set, after the HER2 drop

**Primary analysis frame is 465 patients carrying 103 events.** This is the 470
complete cases less the 5 HER2-positive patients dropped under 5.2, which costs
**1 event**. Both endpoints run on these same 465 patients, so the sign pair in
9.1 is interpretable.

| set | n | events |
|---|---|---|
| all GSE25066 | 508 | 111 |
| pCR complete cases (the 470) | 470 | 104 |
| **primary, after the HER2 drop** | **465** | **103** |

Declared sensitivities, unchanged:

- **490**, the 508 less the 18 without a callable subtype.
- **508**, reachable only by a model without `subtype`.

**DRFS is complete for all 508**, so the primary frame is a strict subset and the
**43 excluded patients carry 8 events.**

### 6.2 Follow-up

Median 2.74 y, IQR 1.74-4.11, max 7.44. By stratum: **3.03 y ER-positive against
2.22 y ER-negative** (max 7.36 / 7.44). This asymmetry is the subject of
section 8.

### 6.3 Proportional hazards

Tested by Schoenfeld residuals and **reported whether or not it holds**. If it
fails for OX, report restricted mean survival time at 3 y **alongside** the hazard
ratio, not instead of it.

## 7. ER stratification - descriptive, and retained

`er_status_ihc` is **primary**: 279 positive and 191 negative in the 470, carrying
40 and 64 events. `subtype` in the H4 frames is built from it, so this keeps the
stratification consistent with the rest of the arm.

The array call `esr1_status` is the **declared sensitivity**, giving 269 / 201
with 37 / 67 events in the 470.

> **CORRECTED 2026-10-02, after E40's first run.** This section previously said
> `esr1_status` *"disagrees on roughly ten patients"*. **Ten is the NET MARGIN
> SHIFT, not the discordance.** In the 465-patient primary frame the two calls
> disagree on **54 patients carrying 11 events** - **22** IHC-negative and
> array-positive, **32** IHC-positive and array-negative - which nets to the 10
> the margins show. **A sensitivity analysis that swaps the ER definition
> reclassifies 54 patients, not 10**, and the sensitivity must be read with that
> figure.

Full counts, recorded before the fit:

| stratum | 470: n | events | **465 (primary): n** | **events** |
|---|---|---|---|---|
| ER-positive (IHC) | 279 | 40 | **278** | **40** |
| ER-negative (IHC) | 191 | 64 | **187** | **63** |
| **total** | **470** | **104** | **465** | **103** |

`esr1_status`, the declared sensitivity, gives **269 / 201 with 37 / 67 events**
in the 470.

The 4 indeterminate and 2 missing all fall outside the 470, so the model set is
cleanly P/N.

### 7.1 `subtype2` is omitted in BOTH stratifications, for two different reasons

**Decided 2026-10-02, after E40's first run stopped on an assertion that
exposed the asymmetry.**

Section 5.2 establishes that `subtype2` is a relabelling of `er_status_ihc`, so
inside an `er_primary` stratum it is **constant** and omitting it is forced.
**That does not carry over to `esr1_status`**: the two calls disagree on 54
patients, so `subtype2` genuinely varies inside the esr1 strata -

| `esr1_status` | n | events | HRpos_HER2neg | TNBC |
|---|---|---|---|---|
| P | 268 | 37 | 246 | **22** |
| N | 197 | 66 | **32** | 165 |

**`subtype2` is omitted there as well, by decision rather than by constancy**, so
that the sensitivity is the same model as the primary stratification and the two
differ only in how patients are assigned to strata. Fitting it in one and not
the other would confound a change of stratifier with a change of adjustment.

**The cost is recorded rather than hidden**: in the esr1 strata the estimates are
**unadjusted for subtype**, and the 22 and 32 reclassified patients carry their
subtype imbalance into them. **The script asserts constancy only for
`er_primary`** - asserting it for `esr1_status` would be asserting something
false - and reports the composition above beside the esr1 output so the omission
is visible where it is read.

**ER-stratified estimates are reported as descriptive, not as a test.** A formal
`ER x OX` interaction is NOT declared and will not be fitted: at 104 events it is
not estimable, and the published prior does not fix a direction (Tao 2017
meta-analysis; Keam 2011 reporting the discordance inside TNBC against Wu 2022
reporting it in ER-positive disease and not in TNBC).

The stratification is retained because the ER-negative stratum is where this
cohort has its information - 63 of 104 events in 187 of 470 patients - and
reporting it pooled would hide that.

## 8. Detectable effect, declared before the fit

For a standardised continuous exposure at 80% power, alpha 0.05:

| stratum | events | detectable HR per SD | with adjustment inflating variance twofold |
|---|---|---|---|
| **all, 465** | **103** | 1.32 | 1.47 |
| ER-negative | **63** | **1.43** | 1.64 |
| ER-positive | 40 | 1.56 | 1.87 |

**The conclusion is unchanged by the HER2 drop.**

**GSE25066 answers the ER-negative half and cannot answer the ER-positive half.**
The ER-positive stratum detects nothing below HR 1.6 to 1.9, and its follow-up is
shortest (2.22 y against 3.03 y is the wrong way round, but the stratum's events
are few either way) exactly where ER-positive relapse is latest.

**A null in the ER-positive stratum of GSE25066 carries no information and must be
reported in those words.** This is the D5 section 6.3 rule applied in advance:
not powered to exclude a modest effect.

## 9. Reading rules

These fix how numbers are read, not whether anything passes.

- A coefficient is read as a signal only where its **CI excludes zero** AND, for
  the three-cohort pCR analysis, **the sign agrees across all three cohorts**.
  Three cohorts disagreeing in sign is read as heterogeneity, as in Block F1
  section 2, where the pooled signal came entirely from the smallest cohort while
  the primary cohort sat at exactly zero.
- At **I2 >= 50%** the random-effects estimate is the summary and the
  fixed-effect estimate is reported beside it, never instead of it.
- **No p value is described as a near miss.** Block F1's 0.052 is the precedent
  for why.
- Every gap is reported with both component positions (standing rule).
- Negative results are explicitly recorded (standing rule).

### 9.1 The sign pair - descriptive only

Read across the two endpoints as follows, and only descriptively:

| pCR | DRFS | consistent with |
|---|---|---|
| better | worse | the proliferation route |
| no effect | worse | the metastasis route the mouse data predict |
| better | better | an unbuffered apoptotic liability |

**None of these three is a test of anything.** The table exists so that the
reading is fixed before the numbers arrive, not so that a pattern can be claimed.

## 10. METABRIC - declared IN

METABRIC is **declared in, not conditional**, on design grounds fixed before any
outcome variable is touched: the ER-positive stratum is not estimable in GSE25066
(section 8), and twenty years of follow-up is the only available instrument for
late ER-positive relapse. SCAN-B carries no survival endpoint of any kind
(`scanb_pheno.rds$pheno`, 3,207 x 17, no time or event column), and TCGA survival
is forbidden by plan section 3.

Deciding this now rather than after GSE25066 runs makes it a declared division of
labour and not a rescue.

**Division of labour, fixed here.** GSE25066 carries the sign-pair argument in
ER-negative disease, where both endpoints exist in the same patients. METABRIC
carries the ER-positive half.

### 10.1 Specification

- Expression and clinical from the **cBioPortal datahub GitHub-LFS mirror**, with
  a provenance README in the house pattern. The S3 datahub returns 403.
- **Breast-cancer-specific survival is primary**; OS is a sensitivity, because
  twenty years of follow-up in an older cohort loads OS with non-cancer death.
- Same ladder, same within-cohort standardisation, same reading rules.
- Proliferation scored the same way, from the same 318-gene set.

### 10.2 Forkscale rung - SECONDARY, and it is NOT F3

MB1 forkscale enters as an additional rung, asking whether the respiratory axis
adds anything over the published fork axis. F3-pre left this **INTERMEDIATE** at
Spearman 0.529 (GSVA) and 0.418 (mitoPPS) - above the independence line and below
the redundancy line on both instruments.

**This is NOT F3.** F3 was specified in `myc_human_validation`, which is frozen at
`d3ac60e`. This is a separate exploratory analysis and must be named as one in
every sentence that mentions it.

### 10.3 Forkscale traps, verified 2026-10-02

The fork data is `data/menegollo_biclusters/METABRIC_starting_data_final_corr_groups.Rdata`,
untracked and gitignored by name. `all.clinical.df` is 1,981 x 58 and carries
`sample` plus the forkscale columns and fork calls. Identifiers are of the form
`MB-0000`, 1,981 unique of 1,981, matching the cBioPortal METABRIC patient-ID
form.

Three traps, carried forward:

1. **One `Inf` per `.log` column**, at different rows: MB1 row 290, MB2 row 85,
   MB3 row 848. Drop the affected row from any analysis using that variant and
   record which.
2. **`MB3.pc1` is NOT negated in this object.** The negated copy is stored
   separately as `MB3.pc1.rev`. `cor(MB3.pc1, MB3.forkscale) = -0.279` and
   `max|MB3.forkscale - MB3.pc1/MB3.index| = 87.31`. **This corrects the G3
   note**, which records METABRIC as negating `MB3.pc1` before forming forkscale;
   following G3 literally gives a reversed axis with no error thrown. MB1 and MB2
   both reproduce exactly (max difference 0), matching G3's TCGA finding.
3. **Forkscale is severely skewed by construction**, so Spearman and not Pearson.

Not verified, because it is not on the critical path: whether
`MB3.forkscale == MB3.pc1.rev / MB3.index`. MB3 is not used by the secondary
rung. Verify before any use of MB3.

## 11. What this cannot license

- **Nothing about H1 through H4.** All four have been tested and none is
  supported. This note does not reopen any of them.
- **No causal claim.** The mouse manipulations are the controlled experiment;
  these cohorts are confounded observation. The asymmetry is deliberate and is
  stated rather than hidden.
- **No claim of two-instrument agreement** (section 4).
- **No comparison of any coefficient across cohorts as a raw value** (section 4.2).
- **No reinstatement** of any claim on the retired list in the human-arm handoff.
- **No "primed" or "priming" language** (standing rule N3).
- **No per-gene FDR** (standing rule).

## 12. Dependencies

**`E39` and `E40` introduce `survival` as the first survival dependency in this
arm.** Established by search before it was added: `coxph|survfit|Surv(` returns
**zero hits across every branch of all three repos** - `myc_human_exploratory`,
`myc_human_validation` and `myc_mouse` - so no script anywhere has fitted a
survival model or joined a time-to-event endpoint before this one.

`survival` is added to `.pkg_analysis` in `scripts/E00_setup_packages.R`, the
tier that is **checked for presence and attached by the scripts that need it**,
in the same commit as this amendment. It is not added to `.pkg_core`: `E39`
itself fits nothing and does not need it attached, and only `E40` will.

`survival` ships with R and is not a new install; the entry exists so that a
missing or broken installation fails at `E00` with a named package rather than
at the first `coxph()` call.
