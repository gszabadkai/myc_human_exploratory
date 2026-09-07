# E22 - the FOXO3-repression arm does not transfer to human tumours

**2026-09-07. EXPLORATORY. Nothing here is pre-registered.**

Script: `scripts/E22_bbc3_direct_myc_effect.R`. **Sourced by the author
2026-09-07**, and `results/bbc3_direct_myc_effect.rds` plus its two tables are
on disk. The numbers below were written from a dry run beforehand; the real run
reproduces **all 25 non-timestamp saved objects identically**, so every number
here is now read from a saved object and **the dry-run caveat is discharged**.

**N3 throughout.** Every number here is a transcript association. The word
"primed" is never written of a transcript.

---

## WHAT THIS IS NOT

> This is a MAIN-EFFECTS model. It is a different estimand from the
> pre-registered `MYC:OXPHOS` interaction on `PRIME`, which is null and
> permanently closed at `d3ac60e`.
>
> - No `M:O` product term appears anywhere.
> - `PRIME` is not computed.
> - Nothing was read, referenced or imported from `myc_human_validation` - the
>   purity-provenance question was deliberately settled within this repo instead.
> - **Nothing produced here bears on the registered interaction result.**

---

## 0. The four things this note carries

| | |
|---|---|
| **1. The reportable negative** | **The FOXO3-repression arm does not transfer to human tumours.** Reading (ii) of the pre-specification. `beta1` does not merely go to zero on adjustment - it **reverses sign**, and within Luminal it is positive in 0 of 6 cells. |
| **2. The positive result** | **`beta2 > 0` in 50 of 50 cells**, every interval excluding zero, +0.108 to +0.550, and **essentially unmoved by adjustment** (`d beta2 = -0.015` on the primary cell). OXPHOS tracks `BBC3` at fixed MYC, robustly. |
| **3. The marginal correction** | **`BBC3` is weakly NEGATIVELY correlated with MYC, 10 of 10 cells - not uncorrelated.** A manuscript wording change, not a footnote. |
| **4. A standalone finding** | **`rho(MYC transcript dose, ox_rel) < 0`**, interval excluding zero, in four of four pooled and Luminal cells, while every activity estimator is positive. **This makes D2 consequential for the H4 STATE definition.** |

---

## 1. The pre-specification, and which reading fired

Fixed in the script header **before** the adjusted fits were run:

```
beta1 = partial MYC    -> BBC3 at fixed OXPHOS    PREDICTED NEGATIVE
beta2 = partial OXPHOS -> BBC3 at fixed MYC       PREDICTED POSITIVE
```

- **(i) SUPPORTED** - `beta1` negative, interval excluding zero, on all three
  claim estimators after `PROLIF_DISJOINT` + purity, pooled **and** within
  Luminal.
- **(ii) NOT TRANSFERRED** - `beta1 -> 0` on adjustment where it was not clearly
  present unadjusted on the same samples. A reportable negative.
- **(iii) AMBIGUOUS** - `beta1 -> 0` having been present unadjusted.
  Proliferation may be confounder or mediator and the model cannot separate
  them. Under (iii) the FOXO3 mediation leg is the tiebreaker, gating on the
  **unadjusted** `beta1`.

**The classifier fired (ii) on both rulers. The stage 3 gate is CLOSED, and the
decision to accept that was taken deliberately rather than mechanically** -
section 2 is the reasoning, and it is recorded because the classifier's label
understates the case.

The ladder runs on **one fixed subset per cohort** - TCGA n = 1,020
purity-complete, SCAN-B n = 3,207 - so (ii) and (iii) are separated by a
computed fact rather than a judgement. Comparing an adjusted value at n = 1,020
against an unadjusted one at n = 1,095 would confound the covariate with the
sample set (E16 section 6's discipline).

---

## 2. Why (iii) cannot apply - four arguments, none of them a technicality

The classifier reached (ii) because `all(unadjusted negative)` was **5 of 6**,
not 6 - the failing cell being **TCGA dose, unadjusted, n = 1,020: `-0.018`
`[-0.077, +0.041]`**. Taken alone that is a borderline cell and a thin reason to
choose between two readings. It is not the reason.

**(a) The mediator argument has no purchase on estimators that never carried the
effect.** (iii) is an argument about what adjustment *destroys*. Two readouts
show **no negative `beta1` unadjusted at all**: TCGA dose (above), and
`sig_clean` - `MYC_UP.V1_UP`, the least proliferation-entangled signature in the
panel at 1.5 pct - which is `-0.007 [-0.069, +0.055]` in TCGA and **`+0.114`
`[+0.077, +0.150]`, positive with an interval excluding zero**, in SCAN-B.
Adjustment cannot have destroyed what was not there.

**(b) Mediation attenuates. It does not cross zero and land significantly
opposite.** On the signature estimator `beta1` goes `-0.200 -> +0.136` in TCGA
and `-0.060 -> +0.106` in SCAN-B, **intervals excluding zero on the positive
side**. That is not the signature of removing a mediated path.

**(c) The covariate acts on the MYC readout, not on a shared causal node.**
`PROLIF_DISJOINT` moves `beta1` by **+0.332** and `beta2` by **-0.015** in the
same cell. A mediator on the `MYC -> BBC3` path would not leave the OXPHOS arm
untouched while moving the MYC arm by twenty times as much.

**(d) `PROLIF_DISJOINT` is the wrong variable to be the FOXO3 mediator.** The
proposed path is `MYC -> PI3K/AKT -> FOXO3 nuclear exclusion -> BBC3`.
`PROLIF_DISJOINT` is `PROLIF_STD` (HALLMARK `E2F_TARGETS` + `G2M_CHECKPOINT`,
327 genes) minus the 9 proliferation genes of `FELSHER_61`. Measured here:
**318 of 318 of its genes are inside the E2F+G2M union, 5 of 318 are in
`HALLMARK_PI3K_AKT_MTOR_SIGNALING`, and `FOXO3` is not in it.** It is a
cell-cycle *output* composite, not a PI3K/AKT signalling composite. It cannot
stand in for the mediator the hypothesis names.

**And (e), which is not an argument about mediation at all: Luminal is 0 of 6.**
In the stratum the effect had to hold in, all three claim estimators are
**positive** in TCGA with intervals excluding zero - dose `+0.157`, regulon
`+0.207`, signature `+0.218` - and no cell in either cohort is negative.

**Conclusion.** (ii) stands, and it stands on substance. The stage 3 FOXO3
mediation leg is not run.

---

## 3. Result 1 - the reportable negative, and its scope

### The ladder on `BBC3`, `ox_rel`

| cohort | estimator | unadjusted | `+PROLIF` | `+PROLIF+purity` |
|---|---|---|---|---|
| TCGA | dose | -0.018 [-0.077,+0.041] | +0.071 | **+0.087** [+0.031,+0.143] |
| TCGA | **signature** | **-0.200** [-0.260,-0.140] | +0.132 | **+0.136** [+0.058,+0.214] |
| TCGA | regulon | **-0.211** [-0.273,-0.149] | -0.025 | +0.078 [-0.003,+0.159] |
| SCAN-B | dose | **-0.043** [-0.077,-0.009] | **-0.036** [-0.070,-0.002] | - |
| SCAN-B | **signature** | **-0.060** [-0.096,-0.023] | **+0.106** [+0.062,+0.150] | - |
| SCAN-B | regulon | **-0.168** [-0.203,-0.133] | **-0.101** [-0.139,-0.064] | - |

`ox_lvl` gives the same pattern. **Purity does almost nothing**: `|d beta1| <=
0.02` everywhere except the regulon (`+0.103` on `ox_rel`, `+0.077` on
`ox_lvl`). Leukocyte fraction, carried as a **named sensitivity** on its own
n = 1,007 subset and not as the reported model, moves `<= 0.003` except the
regulon (`+0.055`).

**The SCAN-B bound, and the reasoning rather than just the number.** SCAN-B has
no purity estimate at all (0 of 3,207) and nothing was substituted. The only
evidence available about what that costs is the TCGA purity rung, and it is used
as a **bound**, not a correction: purity's largest effect on `beta1` anywhere in
TCGA is **0.103**, and for the two signature estimators **<= 0.007**. Since the
SCAN-B conclusion does not turn on a shift of that size - its signature estimate
moves from `-0.060` to `+0.106`, a swing of 0.166 - **the SCAN-B PROLIF-only
fits are interpretable**, and that is why they are reported rather than withheld.

### Within stratum

TCGA Luminal (n = 682): dose `+0.157 [+0.082,+0.232]`, regulon `+0.207
[+0.095,+0.318]`, signature `+0.218 [+0.122,+0.313]`. SCAN-B Luminal
(n = 2,436): signature `+0.164` positive; regulon and dose span zero.
**0 of 6 Luminal cells negative.**

Basal is all-positive with interval widths **0.18 to 0.65** (TCGA n = 162,
SCAN-B n = 317). **A Basal null is weak evidence and is written as such, never
as failure to hold**; the widths are quoted for exactly that reason.

### Specificity - and it relocates the effect

Top rung, `ox_rel`, 3 claim estimators x 2 cohorts = 6 cells per gene:

| gene | side | cells `beta1 < 0` | interval excludes 0 | range |
|---|---|---|---|---|
| **`BCL2L1`** | anti | **6** | **6** | [-0.315, -0.049] |
| `BIK` | pro | 5 | 5 | [-0.246, +0.025] |
| `BCL2L11` | pro | 6 | 5 | [-0.183, -0.013] |
| **`BBC3`** | pro | **2** | **2** | [-0.101, +0.136] |
| **`MCL1`** | anti | **0** | **0** | [+0.010, +0.299] |
| `BID` | pro | 0 | 0 | [+0.093, +0.530] |

**A negative direct MYC effect at fixed OXPHOS is not a generic property of this
roster** - `MCL1`, `BID` and `BCL2A1` are 0 of 6. **But `BBC3` is not the gene
that carries it. `BCL2L1` is**, 6 of 6 on both rulers. The effect the FOXO3 arm
predicts for the trigger is found on the **anti-apoptotic guardian** instead.

### The scope of the negative

This is a **cross-sectional, post-selection, transcript-level** result and none
of those three words is decoration.

- **Cross-sectional.** Every tumour is one timepoint. A direct MYC effect that
  operates transiently, or during transformation rather than in an established
  tumour, would not appear here.
- **Post-selection.** These are tumours that survived. A `MYC -> FOXO3 -> BBC3`
  arm that kills the cells carrying it is selected against, so its absence in a
  prevalent tumour cohort is expected under the hypothesis, not only under its
  negation. **This is the single most important limitation and it is not
  resolvable in this design.**
- **Transcript-level.** FOXO3 repression by MYC is proposed to act through
  **nuclear exclusion**, a post-translational event. A transcript-level model
  cannot see it, and `FOXO3` mRNA is not a measure of FOXO3 activity. The
  mediation leg was designed to use a regulon score precisely because of this,
  and it is not run.

**So the negative is about the human tumour transcriptome, not about the
biology of the mouse gland.** It says the arm does not transfer to this
measurement in this material - it does not say the arm is absent in the mouse,
where it was observed under a different design.

---

## 4. Result 2 - `beta2`, and it is the positive result

**`beta2 > 0` in 50 of 50 cells** - 5 estimators x 2 rulers x every rung x 2
cohorts - with **every interval excluding zero**, ranging +0.108 to +0.550.

Top rung, `ox_rel`:

| cohort | dose | signature | regulon | sig_clean | sig_entangled |
|---|---|---|---|---|---|
| TCGA | +0.336 | +0.304 | +0.305 | +0.263 | +0.350 |
| SCAN-B | +0.214 | +0.201 | +0.231 | +0.166 | +0.259 |

**The robustness statement, and it is the one to quote:** under adjustment for
proliferation, `beta2` moves by **`-0.015`** on the primary cell (TCGA,
signature, `ox_rel`), against `+0.332` for `beta1` in the same fit. Across all
cells `d beta2` from proliferation runs `-0.087` to `+0.089` and from purity
`-0.029` to `+0.006`. **The OXPHOS arm is not a proliferation artefact and it is
not a purity artefact**, and unlike `beta1` it does not depend on which MYC
estimator is used.

**The manuscript sentence this licenses**, and only this one:

> *In human breast tumours, `BBC3` transcript abundance rises with the OXPHOS
> axis at fixed MYC activity, in two cohorts, on both respiratory rulers and on
> all five MYC estimators, and the association is unchanged by adjustment for
> proliferation and tumour purity.*

**What it does NOT license.** It is not a statement about apoptotic priming, or
about protein, or about the `MYC x OXPHOS` interaction. It is a partial
association among transcripts. **N3 applies: do not write "primed".**

---

## 5. Result 3 - the marginal correction

`rho(BBC3, MYC)`, unadjusted, Spearman beside the Blom-Pearson the identity uses:

| cohort | dose | signature | regulon | sig_clean | sig_entangled |
|---|---|---|---|---|---|
| TCGA | -0.048 | -0.117 | -0.105 | **+0.064** | -0.174 |
| SCAN-B | -0.057 | -0.004 | -0.158 | **+0.130** | -0.136 |

**Negative in 10 of 10 estimator-by-cohort cells on Spearman**, and Blom-Pearson
agrees (largest gap 0.043). The two positive-reading cells are `sig_clean`,
which is a *sensitivity* signature and not one of the three that carry the
claim - and it is named as a tension in section 7, not hidden.

> **MANUSCRIPT WORDING CHANGE.** Any sentence saying `BBC3` is *"uncorrelated
> with MYC"* in human tumours is **wrong as written** and must read *"weakly
> negatively correlated with MYC"*. The magnitude is small - `|rho| <= 0.17` -
> so "weakly" is doing real work and should not be dropped either.

### The reconciliation, reported once as accounting

For a two-predictor OLS on equal-variance variables, `beta1 + beta2 * r(M,O) ==
r(Y,M)` **identically**. It holds here to `8.3e-16` in all 20 cells. **It cannot
fail and it is not evidence.** What it delivers is how negative `beta1` has to
be: `beta2 * r(M,O)` is +0.055 to +0.156 on the non-vacuous rows, so `beta1` had
to be -0.060 to -0.265 to land on the observed marginal.

Rows with `|r(M,O)| < 0.05` are flagged **VACUOUS** - all three dose rows - and
must never be read as agreement: with `r(M,O) ~ 0` there is nothing to reconcile
and `beta1` is simply the marginal.

---

## 6. Result 4 - `rho(MYC dose, ox_rel) < 0`, and why it matters elsewhere

Reported in its own right, 5 estimators x 2 rulers x 2 cohorts x 3 strata with
intervals, in `$rho_myc_ox`.

| cell | `rho(dose, ox_rel)` |
|---|---|
| TCGA pooled | **-0.103** [-0.161, -0.044] |
| TCGA Luminal | **-0.163** [-0.235, -0.089] |
| SCAN-B pooled | **-0.048** [-0.083, -0.014] |
| SCAN-B Luminal | **-0.043** [-0.083, -0.003] |

On `ox_lvl` it spans zero almost everywhere (-0.016 / -0.003 pooled). **Every
activity estimator is positive throughout**, +0.11 to +0.55.

**So "MYC-high tumours are OXPHOS-high" is a statement about MYC ACTIVITY that
does not survive being made about MYC DOSE, and on the relative ruler dose
points the other way.** That is CLAUDE.md trap 4 with a sign attached.

**The consequence, and it is why this is logged rather than buried: D2 is
consequential for the H4 STATE definition.** Which estimator defines the state
decides the sign of that state's relationship to the respiratory axis. Any
argument that conditions on a "MYC-high, OXPHOS-high" state - the collider
argument included - inherits an estimator choice it has not declared.
**Logged; not acted on here.**

---

## 7. `sig_clean` - a named tension, not a footnote

`MYC_UP.V1_UP__MITOSTRIP` is the least proliferation-entangled signature in the
panel (1.5 pct, against 14.8 for FELSHER and 23.5 for HALLMARK V1). It is a
reported sensitivity and it **splits from the claim estimators**:

| | unadjusted | top rung | Luminal |
|---|---|---|---|
| `sig_clean` TCGA | -0.007 | **+0.249** | **+0.300** |
| `sig_clean` SCAN-B | **+0.114** | **+0.220** | **+0.260** |
| `sig_entangled` TCGA | -0.367 | -0.077 | -0.045 |
| `sig_entangled` SCAN-B | -0.253 | **-0.200** | **-0.219** |

**`beta1`'s sign tracks proliferation entanglement across the five estimators.**
The cleanest MYC readout says the direct effect is **positive**; the most
entangled says it is **negative**; the three claim estimators sit between them
and change sign on adjustment.

Two things follow and both should be written.

1. **It strengthens the negative.** If the effect were a MYC effect it should be
   *clearest* on the least-contaminated estimator. It is absent there
   unadjusted and positive there adjusted.
2. **It is a referee's first question** and it is answered here rather than
   found later. It is also a fifth independent instance of CLAUDE.md trap 3 -
   the estimators disagree by a factor of thirty in entanglement, and here that
   disagreement changes a sign.

---

## 8. What this note licenses, and what it does not

**Licensed:**

- The `beta2` sentence in section 4, verbatim.
- *"A negative direct MYC effect on `BBC3` at fixed OXPHOS is not supported in
  human breast tumours; the mouse gland's FOXO3-repression arm does not transfer
  to this measurement."* - with section 3's three scope words attached.
- The `BBC3`-MYC wording change in section 5.
- *"Within this BCL2-family roster the gene carrying a concordant negative
  direct MYC effect is `BCL2L1`, not `BBC3`."* - as an **observation from a
  specificity control**, explicitly not as a tested hypothesis.

**NOT licensed:**

- Anything about the `MYC x OXPHOS` **interaction**. Different estimand,
  pre-registered elsewhere, null, closed. **This note does not bear on it.**
- Any claim that the FOXO3 arm is absent **in the mouse**, or absent in human
  biology as opposed to in the human tumour transcriptome.
- Any use of `beta1` from a single estimator. Its sign tracks entanglement
  (section 7) and no single value is interpretable alone.
- **`BCL2L1` as a finding.** It emerged from a specificity control, in the same
  data, after the fact. It is a lead, and section 9 says what would make it more
  than that.
- Anything using the word "primed" of a transcript. **N3.**

---

## 9. What would change these answers

1. **A time course, or a pre-malignant series.** The post-selection limitation
   is structural: a lethal arm is selected against in prevalent tumours. Nothing
   cross-sectional can fix it.
2. **A PI3K/AKT covariate, or a FOXO3 regulon score.** `PROLIF_DISJOINT` is not
   the mediator the hypothesis names (section 2d). The FOXO3 regulon is in the
   snapshot - 196 CollecTRI targets, 179 after the mitochondrial strip - and was
   **not** run because the gate closed. If a PI3K/AKT axis moved `beta1` the way
   `PROLIF_DISJOINT` does, that would be a different and much more relevant
   result.
3. **`BBC3` protein, or FOXO3 localisation.** The proposed mechanism is
   post-translational and this measurement cannot see it.
4. **A third cohort**, since both cohorts here agree and both are bulk breast.
5. **`BCL2L1` as a pre-specified endpoint** rather than a specificity readout -
   see `E23`, which tests the mouse's own `MCL1:BCL2L1` guardian balance and is
   licensed by the mouse gate-model document rather than by this table.

---

## 10. What was NOT done, listed so it stays a decision

No interaction term, no `PRIME`, no Johnson-Neyman. No stage 3 FOXO3 mediation -
the gate closed and it is not run. No per-gene FDR across the twelve, as
specified. `genome_doublings` was available and is **not** a covariate. SCAN-B
purity was **not** imputed. `M_b__FULL` and `M_b__PROLIFSTRIP` were **not**
fitted - both contain `BBC3` and are disqualified as the exposure for this
outcome. Nothing was read from `myc_human_validation`, including for the purity
provenance question. The mouse was not read by the script. No ortholog function
called, in either direction.

---

## 11. Two things about the measurement, recorded because they nearly bit

**Ties break Blom scoring's equal-variance property.** SCAN-B carries ties in
`BCL2A1`, `BIK` and `PMAIP1` - **not** `BBC3` - which spread the variances by
`4.4e-4`. Scores are therefore Blom-transformed **and then scaled to unit
variance**; without that step the section 5 identity would hold only to `~1e-4`.
`$ties` carries the fractions. The identity is claimed for the **pooled**
unadjusted fit only, because stage 2 subsets pooled scores per E19's
no-re-scoring rule and subsetting breaks equal variance again.

**Eleven of the twelve outcome genes sit inside `ox_rel`'s denominator** - all
but `BMF` - with weight `-1/1047`, because they are MitoCarta members. Measured
rather than argued away with a leave-one-out rebuild per gene: worst shift across
120 cells is `|d beta1| = 0.0018` and `|d beta2| = 0.0039`. Standard `ox_rel`
stays primary for comparability with E10, E16, E20 and E21.

**Purity provenance.** No ingest or download script in this repo creates the
column; `E01:202` reads it from `data/from_validation/`, whose README names the
column and no algorithm. **Methods must say "purity, source undocumented in this
repo."** It sits beside `ploidy` and `genome_doublings`, which is the triple
ABSOLUTE emits - recorded as **inference from column composition**, never as
documentation.

---

## 12. Where the numbers live

| | |
|---|---|
| script | `scripts/E22_bbc3_direct_myc_effect.R` |
| object | `results/bbc3_direct_myc_effect.rds` - `$ladder_pooled`, `$deltas`, `$ladder_strata`, `$spec_summary`, `$reading`, `$verdict`, `$rho_myc_ox`, `$marginals`, `$identity_chk`, `$loo`, `$purity_cost`, `$leuko_sens`, `$ties`, `$contamination`, `$settings`, `$rules` |
| tables | `outputs/tables/E22_fit_a.csv`, `outputs/tables/E22_rho_myc_oxphos.csv` |
| figures | none |
| the pre-specification | the script header, fixed before the fits |
| the mouse document, read read-only and cited as conclusions only | `myc_mouse` `docs/2026-09-02_myc_oxphos_priming_gate_model.md` |
