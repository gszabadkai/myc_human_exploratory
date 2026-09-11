# E26 - the MYC separation survives composition adjustment, and the control says the adjustment is real

**Run and written 2026-09-11. EXPLORATORY, POST-HOC, DESCRIPTIVE. Nothing here
is pre-registered and nothing here is a hypothesis test.**

Script: `scripts/E26_purity_and_composition.R`. **Sourced by the author
2026-09-11**; `results/purity_and_composition.rds`,
`outputs/tables/E26_adjusted_contrasts.csv` and
`outputs/mitopps_rank/E26_composition_adjustment.pdf` are on disk. The author's
run reproduces **all 26 saved elements identically** against the dry run,
seeded 5,000-draw bootstraps included, so **the dry-run caveat is discharged**.

**TCGA only, and that is structural.** SCAN-B has no purity estimate and nothing
is imputed (CLAUDE.md trap 2), so this arm cannot replicate. SCAN-B's E25
numbers remain unadjusted and uncheckable. That is a limitation of the data, not
a choice.

---

## 0. The answer in six lines

| | |
|---|---|
| **Does the positive control survive?** | **Yes, and this is the gate.** Proliferation keeps **96%** of its MYC difference under purity and leukocyte adjustment, same sign, interval clear of zero, while those adjusters explain **8.8%** of its variance. |
| **Are the plan's adjusters an adequate composition control?** | **No, and the pre-declared rule caught it.** They leave **93%** of the adipocyte composition difference standing, so the adipose-adjusted models became primary by a rule fixed before any outcome. |
| **Where is the exposure actually?** | **Adipose, not purity.** The MYC groups differ by **0.058** in purity, **-0.002** in leukocyte fraction and **-0.743** standard units in the adipocyte tracer. |
| **Does the separation survive?** | **Yes, entirely.** 44.4 percentile points unadjusted, **43.0** under purity and leukocyte, **47.2** under all three adjusters. |
| **Is the respiratory chain hit harder than the panel?** | **No.** It retains 83% against a panel median of 86%, the **41st percentile** of the panel's attenuation distribution. Not an outlier. |
| **Does the cross-species agreement survive?** | **Yes, under every model.** rho with the mouse Myc effect goes 0.39 and 0.47 unadjusted to **0.33 and 0.42** fully adjusted, every interval clear of zero. |

**The one-sentence reading: the MYC separation in the mitochondrial ordering is
not a composition artefact, and the positive control is what licenses saying
so.**

---

## 1. Why a control had to come first

**Trap 15 is the reason this script exists.** In the mouse fat pad, `Mki67`
went from tau +0.313 to -0.096 after composition adjustment while the adjusters
explained only **15%** of its variance. An adjustment that destroys an
association it barely explains is over-correcting, and the only way to know is
to carry a quantity whose answer is not in doubt **through** the adjustment and
read it first.

**The positive control is proliferation.** MYC drives proliferation. That is not
in doubt and it is not a composition artefact. If purity adjustment destroys it,
the adjustment is removing biology and nothing below it is readable.

| model | MYC-high minus MYC-low | 95% CI | retained | adjusters' R-squared |
|---|---|---|---|---|
| unadjusted | 0.581 | 0.542 to 0.621 | 1.00 | - |
| purity + leukocyte | 0.557 | 0.518 to 0.597 | **0.96** | **0.088** |
| adipose | 0.509 | 0.467 to 0.551 | 0.88 | 0.246 |
| all three | 0.496 | 0.454 to 0.538 | 0.85 | 0.284 |

**GATE 1 PASSED.** This is the mouse diagnostic run in reverse. There, a
comparable share of explained variance destroyed the association. Here it leaves
it almost untouched. **The adjustment is not over-correcting, and that is what
makes every number below readable.**

---

## 2. The exposure check, which reshaped the script

An adjuster cannot change much if the groups do not differ on it. So the
exposure was sized before anything was adjusted.

| covariate | MYC-high | MYC-low | difference | 95% CI |
|---|---|---|---|---|
| purity | 0.599 | 0.541 | **+0.058** | 0.031 to 0.085 |
| leukocyte fraction | 0.226 | 0.228 | **-0.002** | -0.023 to 0.018 |
| adipocyte tracer | -0.349 | 0.394 | **-0.743** | -0.858 to -0.627 |

**The covariates the plan named barely differ between the MYC groups. The
adipocyte tracer differs enormously.** That is the fat-pad concern from the
mouse appearing in the human data, and it is where the risk actually sits.

**GATE 2 fired, as it was written to.** The rule, fixed before any number, was
that an adjustment leaving 50% or more of the tracer's MYC difference is not an
adequate composition control. Purity and leukocyte fraction leave **93%**. They
correlate with the tracer only at -0.266 and +0.200. **So the primary adjusted
reading became the model carrying all three adjusters, by rule rather than by
choice after the fact.**

### 2.1 What the tracer is, and what it cannot do

Sixteen canonical adipocyte markers - `ADIPOQ`, `LEP`, `PLIN1`, `PLIN4`,
`FABP4`, `CFD`, `CIDEC`, `LIPE`, `PPARG`, `CD36`, `AQP7`, `GPD1`, `TRARG1`,
`THRSP`, `LPL`, `PNPLA2` - all present in the TCGA matrix, scored as the mean
z-score per sample.

**It is declared in the script and is not a pinned snapshot.** CLAUDE.md says
consume the snapshots rather than rebuild. That rule governs claim-bearing gene
sets; this is a covariate, never scored, never on the panel, never in a result.
The departure is deliberate and is recorded rather than hidden.

**It cannot separate "less adipose tissue in the section" from "tumour cells
expressing fewer adipocyte genes".** In bulk breast tissue these markers are
overwhelmingly adipocyte-derived, so it is mostly composition - but *mostly* is
the honest word, and an adjustment on it may remove some real biology along with
the fat. That direction of error makes the surviving separation a conservative
reading rather than an inflated one.

---

## 3. The separation survives, and the panel says why that is not luck

| model | MYC-low | MYC-high | separation | retained |
|---|---|---|---|---|
| unadjusted | 27.8 | 72.2 | **44.4** | 1.00 |
| purity + leukocyte | 28.5 | 71.5 | 43.0 | 0.97 |
| adipose | 25.0 | 75.0 | 50.0 | 1.13 |
| **all three** | 26.4 | 73.6 | **47.2** | **1.06** |

**Adjusting for adipose does not shrink the separation. It widens it slightly.**

### 3.1 Attenuation is read against the panel, never alone

Every adjustment shrinks every coefficient somewhat, so a number shrinking is
not evidence. What would be evidence is the primary shrinking **more than the
panel does**. That is the same device as E21's size-matched null.

| model | primary retains | panel median | primary's panel percentile | outlier |
|---|---|---|---|---|
| purity + leukocyte | 0.918 | 0.969 | 23 | no |
| adipose | 0.865 | 0.887 | 42 | no |
| all three | 0.828 | 0.855 | **41** | no |

**The respiratory chain attenuates like a typical member of the panel.** The
adjustment is non-specific, which is what an adjustment removing composition
rather than programme looks like.

**A limitation of the attenuation ratio, recorded.** It is a ratio of
coefficients, so it is unstable where the unadjusted coefficient is near zero:
8 of 142 pathways have an absolute retained value above 2, and the full range
runs from -19.5 to +7.8. **R1 reads only the median and a percentile rank, both
of which are immune to those tails**, and the ratio is never averaged.

---

## 4. Two independent checks that do not rely on adjustment at all

**The high-purity subset.** At purity 0.65 or above, n = 235 (85 MYC-low, 150
MYC-high), the separation is **41.5** percentile points against 44.4 in the full
analysis set. Restricting composition by design gives the same answer as
modelling it.

**The DNA-level split, which composition cannot manufacture.** A copy-number
call is not an expression measurement, so it cannot be produced by stromal or
adipose content in the way an activity signature might be.

| model | not amplified (n=789) | MYC-amplified (n=218) |
|---|---|---|
| unadjusted | 13.7 | 86.3 |
| purity + leukocyte | 13.0 | 87.0 |
| adipose | 12.3 | 87.7 |
| all three | **11.6** | **88.4** |

**The copy-number split widens under every adjustment.** Adjustment makes this
reading stronger, not weaker.

---

## 5. The cross-species agreement survives every model

This is the quantity E25's conclusion rests on, so it is the one that had to
hold.

| model | vs `Myc_effect_6W` | vs `Myc_effect_12W` |
|---|---|---|
| unadjusted | 0.39 [0.22, 0.54] | 0.47 [0.30, 0.61] |
| purity + leukocyte | 0.39 [0.23, 0.54] | 0.47 [0.31, 0.60] |
| adipose | 0.33 [0.17, 0.48] | 0.41 [0.26, 0.55] |
| **all three** | **0.33 [0.17, 0.49]** | **0.42 [0.26, 0.56]** |

**Every interval is clear of zero under every model.** The two timeline
contrasts stay negative throughout, -0.39 to -0.52, so the opposition between
the MYC axis and the developmental axis that E25 reported is also unaffected by
composition adjustment.

---

## 6. One readout that composition does hit

**`mtDNA-encoded OXPHOS subunits` is the exception and it is worth recording.**
Its MYC contrast retains 1.08 under purity and leukocyte adjustment but only
**0.243** and **0.271** under the adipose models, the 12th and 15th percentile
of the panel. It is the only readout whose behaviour changes character when
adipose enters.

That is consistent with what the pathway is: adipose tissue is mitochondrion-
rich and mtDNA-encoded transcripts dominate by abundance (trap 8), so a large
part of this pathway's apparent MYC association is adipose content. **The
nuclear-encoded primary set does not behave this way**, which is the contrast
that matters, and trap 8's rule that the two are never pooled is what made the
distinction visible at all.

---

## 7. What this does not settle

**SCAN-B cannot be checked.** It carries no purity estimate and nothing is
imputed. Its E25 numbers stay unadjusted. The only reassurance is indirect: the
two cohorts agreed at +0.88 in E25, and TCGA - the one that can be checked -
survives. **That is not the same as SCAN-B surviving.**

**The tracer is one instrument.** A deconvolution-based stromal and adipose
estimate would be an independent second opinion and is not run here.

**Genome doublings are not a covariate** - the author ruled that out - and
ploidy is available and deliberately unused. Only the composition adjusters the
question is about enter the models.

**The coverage effect is not an adjustment effect.** The adjusted models run on
the **675** samples with complete covariates, not E25's 730, and that alone
moves the unadjusted MYC-low percentile from **16.5 to 27.8**. Every adjusted
row is scored against E26's own unadjusted baseline and never against E25's
numbers. GATE 0 proves the machinery separately: intercept-only residuals
reproduce E25's percentiles at **max |delta| = 0**, because mitoPPS is centred
at exactly 1 per pathway and the two orderings are the same object.

---

## 8. Constraints observed

- **TCGA only**, and the asymmetry is stated wherever a number appears.
- **The control is read first** and is a gate, per trap 15.
- **Attenuation is read against the panel**, never alone.
- **Rank only** across the species boundary; mouse mitoPPS not recomputed;
  artefact md5 re-checked at run time.
- **GSVA and mitoPPS never share an input object**: the control is GSVA on VST,
  the panel is mitoPPS on linear DESeq2-normalised counts.
- **`mt-` / `MT-` held separately** (trap 8), which is what made section 6
  visible.
- **No ortholog call anywhere.** Tripwire clean.
- **N3.** Nothing here is called primed.
- **The figure carries no significance marks.**

## 9. Provenance

- Script: `scripts/E26_purity_and_composition.R`, seed 1, 5,000 bootstrap draws.
- Covariates: `data/from_validation/tcga_brca_covariates.rds` (purity,
  leukocyte fraction, `MYC_amp`).
- Control: `PROLIF_DISJOINT` from `tcga_brca_mito_scores.rds$gsva_cov`.
- Tracer: 16 markers from `data/from_validation/tcga_brca_linear.rds`.
- Panel and groups: inherited from `results/mitopps_group_rankings.rds`.
- Saved: `results/purity_and_composition.rds` (26 elements),
  `outputs/tables/E26_adjusted_contrasts.csv`,
  `outputs/mitopps_rank/E26_composition_adjustment.pdf`.
