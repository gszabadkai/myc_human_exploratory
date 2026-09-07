# E23 - the guardian balance is a second negative, not a confirmation

**2026-09-07. EXPLORATORY. Nothing here is pre-registered.**

Script: `scripts/E23_guardian_balance.R`. **Sourced by the author 2026-09-07**,
after `E22`, and `results/guardian_balance.rds` plus its table are on disk. The
numbers below were written from a dry run that read `E22`'s dry-run object; the
real run - reading `E22`'s real object - reproduces **all 18 non-timestamp saved
objects identically**, so every number here is now read from a saved object and
**the dry-run caveat is discharged**.

**N3 throughout.** These are transcript associations. "Primed" is never written
of a transcript.

---

## WHAT THIS IS NOT

> A MAIN-EFFECTS model, and a different estimand from the pre-registered
> `MYC:OXPHOS` interaction on `PRIME`, which is null and permanently closed at
> `d3ac60e`. No product term, no `PRIME`, nothing read from
> `myc_human_validation`, **and nothing here bears on the registered
> interaction result.**

---

## 0. The answer in five lines

| | |
|---|---|
| **The verdict** | **SECOND NEGATIVE.** Both pre-specified failure readings fired. This is written as a negative and **not** softened into "attenuated". |
| **Why it looks like a confirmation at first** | `beta1 > 0` with intervals excluding zero in **5 of 6** pooled cells on `ox_rel` and **6 of 6** on `ox_lvl`. A headline reader would call this confirmed. |
| **F1, entanglement tracking** | **FIRED** on `ox_lvl` (monotone in both cohorts, spread 0.354). Its positive control - E22's `BBC3`, known to be entanglement-tracking - fires on both rulers, so the rule is working. |
| **F2, Luminal reversal** | **FIRED on both rulers.** `beta1 <= 0` within Luminal in **4 of 6** cells on `ox_rel` and 2 of 6 on `ox_lvl`. |
| **The decisive cell** | **`sig_clean`, the least proliferation-entangled estimator, gives `beta1` NEGATIVE**: `-0.280 [-0.315, -0.244]` pooled in SCAN-B, `-0.377 [-0.414, -0.340]` in SCAN-B Luminal and `-0.271 [-0.343, -0.199]` in TCGA Luminal. **The prediction was `beta1 > 0`.** |

**The same artefact as `BBC3` in E22: `beta1`'s sign tracks how much
proliferation the MYC estimator carries.** Two endpoints, opposite predicted
directions, and the same failure.

---

## 1. What licensed this, what did not, and why it was never an independent test

**LICENSED BY** `myc_mouse` `docs/2026-09-02_myc_oxphos_priming_gate_model.md`,
read read-only and cited as conclusions only. There
`B = log2(Mcl1) - log2(Bcl2l1)` is the guardian balance, promoted to co-primary
in section 4.3, and the model's most robust endpoint - top of section 3.3's
menu in both the full and the within-timepoint regime, with `Mcl1` at the 98.9th
percentile of its expression-matched null and `Bcl2l1` at the 2.1st.

**NOT LICENSED BY** E22's specificity table. That is where the endpoint was
*noticed*, and noticing is not licensing. The distinction is not pedantry - it
is what stops the same data being used twice.

**AND IT WAS NEVER AN INDEPENDENT TEST IN THIS DATA, which the script header
says before any number.** E22 had already fitted both halves in the same
samples: `BCL2L1` negative in 12 of 12 cells with every interval excluding zero,
`MCL1` positive in 12 of 12 point estimates. A difference of those two is
positive with near-certainty **before anything is fitted**. Confirming
`beta1 > 0` here was therefore worth almost nothing, and section 5 measures
exactly how little.

**The estimand also differs from the mouse's.** The mouse's `bMX` is an
*interaction* coefficient; `beta1` here is a *main effect*. What transferred was
the predicted **sign**, never the estimand. No interaction is fitted anywhere.

**So the script's whole value was that the battery could fail. It did.**

---

## 2. The result that would have been reported

The ratio, every rung, claim estimators, `ox_rel`:

| cohort | estimator | unadjusted | `+PROLIF` | `+PROLIF+purity` |
|---|---|---|---|---|
| TCGA | dose | +0.213 [+0.157,+0.269] | +0.194 | **+0.168** [+0.112,+0.225] |
| TCGA | signature | +0.210 [+0.152,+0.269] | +0.234 | **+0.226** [+0.146,+0.306] |
| TCGA | regulon | +0.381 [+0.324,+0.439] | +0.413 | **+0.387** [+0.307,+0.466] |
| SCAN-B | dose | +0.054 [+0.021,+0.087] | **+0.048** [+0.016,+0.081] | - |
| SCAN-B | signature | +0.125 [+0.090,+0.159] | **+0.036** [-0.006,+0.079] | - |
| SCAN-B | regulon | +0.302 [+0.269,+0.335] | **+0.279** [+0.244,+0.314] | - |

Positive throughout, stable under adjustment, both cohorts. On `ox_lvl` all six
top-rung intervals exclude zero. The leave-one-out - both halves removed from
`ox_rel`'s denominator, since both are MitoCarta members - moves nothing
(`<= 0.002`).

**Read alone, that is a confirmation.** The battery is what stops it being
written as one, and this is the clearest argument this repo has produced for
pre-specifying a battery rather than a single coefficient.

---

## 3. F1 - entanglement tracking, with its positive control

`beta1` across the three signature estimators ordered by proliferation
entanglement (1.5, 14.8, 23.5 pct), top rung:

| | `sig_clean` 1.5 | `signature` 14.8 | `sig_entangled` 23.5 |
|---|---|---|---|
| TCGA `ox_rel` | **-0.051** | +0.226 | +0.186 |
| SCAN-B `ox_rel` | **-0.280** | +0.036 | +0.183 |
| TCGA `ox_lvl` | +0.044 | +0.332 | +0.397 |
| SCAN-B `ox_lvl` | **-0.213** | +0.112 | +0.357 |

| target | ruler | cohorts monotone | min spread | **fires** |
|---|---|---|---|---|
| **E22 `BBC3` (control)** | `ox_rel` | 2 of 2 | 0.326 | **YES** |
| **E22 `BBC3` (control)** | `ox_lvl` | 2 of 2 | 0.450 | **YES** |
| `MCL1_minus_BCL2L1` | `ox_rel` | 1 of 2 | 0.277 | no |
| **`MCL1_minus_BCL2L1`** | `ox_lvl` | **2 of 2** | **0.354** | **YES** |

**The positive control fires on both rulers.** That matters: F1 is a heuristic,
not a test - with three points, monotone-by-chance is 1 in 3 per cohort - so it
was pre-specified with a control precisely so that a non-firing result could not
be read as reassurance. It detects the pattern it was built for.

**F1 fires on `ox_lvl` and not on `ox_rel`**, where TCGA is non-monotone
(`-0.051`, `+0.226`, `+0.186`). But the failure does not rest on the monotone
count. **On `ox_rel` the cleanest estimator is `-0.051` in TCGA and `-0.280` in
SCAN-B - the wrong sign, decisively so in SCAN-B - while the most entangled is
`+0.186` and `+0.183`.** The ordering is broken only by the middle estimator
overshooting the top one in TCGA.

---

## 4. F2 - Luminal reversal, and the decisive cell

`ox_rel`, top rung, within Luminal:

| cohort | dose | signature | regulon | `sig_clean` | `sig_entangled` |
|---|---|---|---|---|---|
| TCGA (n = 682) | -0.060 | -0.033 | +0.099 | **-0.271** [-0.343,-0.199] | -0.071 |
| SCAN-B (n = 2,436) | -0.003 | **-0.077** [-0.123,-0.032] | +0.131 | **-0.377** [-0.414,-0.340] | +0.114 |

**`beta1 <= 0` in 4 of 6 claim-estimator cells on `ox_rel` and 2 of 6 on
`ox_lvl`. F2 fires on both rulers.** Only the regulon stays positive in Luminal,
and only on one ruler with an interval that covers zero in TCGA.

**And the decisive cell is `sig_clean` in Luminal**: `-0.271` and `-0.377`, both
intervals comfortably excluding zero, in the large homogeneous stratum where
stromal and immune composition vary least. The least-contaminated MYC readout,
in the cleanest stratum, gives the **opposite** of the predicted sign in both
cohorts.

**Basal is all over the place with wide intervals** - TCGA n = 162, SCAN-B
n = 317, widths **0.195 to 0.571** - and the two cohorts disagree (TCGA dose
`+0.287`, regulon `-0.025`; SCAN-B dose `+0.083`, regulon `+0.211`). **A Basal
result at these widths is weak evidence and is written as such**, in either
direction. It neither rescues nor condemns the endpoint.

---

## 5. Does the ratio add anything to its two halves?

Measured, because the header claimed it would not. Top rung, all rulers and
estimators, 30 cells:

| | |
|---|---|
| `sign(beta1_ratio) == sign(beta1_MCL1 - beta1_BCL2L1)` | **30 of 30** |
| median `\|ratio - halves difference\|` | **0.120** |
| direction of the gap | the ratio is systematically **smaller in magnitude** than the halves difference |

**The ratio's SIGN is fully determined by the two coefficients E22 already
reported.** Its magnitude is attenuated by roughly a third - because
Blom-scoring a difference is not the difference of Blom scores - but no cell
anywhere changes direction. So as evidence about MYC this endpoint is a
re-expression of E22, exactly as the header said it would be, and its
contribution is the **framing** one: the human arm is now on the mouse's own
endpoint rather than a human-chosen surrogate.

**Which member carries it**, top rung, `ox_rel` - and it is `BCL2L1` as much as
`MCL1`:

| cohort | estimator | `MCL1` | `BCL2L1` |
|---|---|---|---|
| TCGA | dose | +0.112 [+0.051,+0.173] | -0.153 [-0.209,-0.097] |
| TCGA | signature | +0.190 [+0.104,+0.276] | -0.176 [-0.255,-0.097] |
| TCGA | regulon | +0.299 [+0.212,+0.386] | -0.315 [-0.394,-0.235] |
| SCAN-B | dose | +0.032 [-0.002,+0.065] | -0.051 [-0.084,-0.019] |
| SCAN-B | signature | +0.010 [-0.035,+0.054] | -0.049 [-0.092,-0.006] |
| SCAN-B | regulon | +0.226 [+0.189,+0.263] | -0.221 [-0.257,-0.185] |

**100 of E22's single-gene cells reproduce here exactly** (max `|delta| = 0`),
which is the control that these are E22's fits and not a lookalike.

---

## 6. `beta2` - unpredicted, and an unanticipated sign agreement

**`beta2` was deliberately not predicted**, and it is reported and never scored.
The mouse's own statement for this endpoint is a *crossover*: `bX = -0.541` at
MYC-null with the sign reversed by MYC, so pooled over a mixed-MYC tumour cohort
its direction is not fixed in advance.

Observed: **`beta2` is negative throughout**, `-0.314` to `-0.463` across every
cohort, rung and claim estimator. So in human tumours the OXPHOS axis is
associated with a **lower** `MCL1:BCL2L1` balance at fixed MYC.

**That agrees in sign with the mouse's `bX` at MYC-null.** It is recorded
because it is striking and because someone will notice it - **but it was not
predicted, it was not scored, and it must not be presented as a confirmation.**
A quantity nobody committed to in advance, in a cohort pooled over the whole MYC
range, agreeing in sign with a mouse coefficient defined at MYC-null, is an
observation and a candidate hypothesis. **If it is to become anything, it must
be pre-specified somewhere else first.**

---

## 7. What this licenses, and what it does not

**Licensed:**

- *"The mouse gate model's guardian balance `MCL1:BCL2L1` does not transfer to
  human breast tumours as a MYC main effect: the association is positive pooled
  but reverses within Luminal and its sign tracks the proliferation entanglement
  of the MYC estimator."*
- *"Two endpoints with opposite predicted directions - `BBC3` and
  `MCL1:BCL2L1` - fail in the same way, and the shared failure is that `beta1`'s
  sign is a function of which MYC estimator is used."* **This is the more
  general finding and it is worth more than either endpoint.**

**NOT licensed:**

- **Any confirmatory reading of the pooled `beta1 > 0`.** It was near-guaranteed
  by E22 before it was fitted (section 1), and the battery it was pre-specified
  against fails.
- Any statement about the `MYC x OXPHOS` **interaction**. Different estimand,
  closed elsewhere.
- **The word "attenuated".** The pre-specification says a firing failure reading
  is a second negative and is written as one.
- The `beta2` / mouse `bX` sign agreement as evidence (section 6).
- Anything about `BCL2L1` alone. **The pre-specified object was the pair**, and
  it stays the pair.
- Anything using "primed" of a transcript. **N3.**

---

## 8. What would change this answer

1. **A MYC estimator panel that is not confounded with proliferation.** This is
   now the binding constraint on the whole arm - two endpoints have failed on
   it. `beta1`'s sign is a function of estimator entanglement, and until that is
   broken no main-effect claim about MYC on any BCL2-family transcript is
   interpretable.
2. **A proliferation covariate that is not a cell-cycle *output* composite.**
   `PROLIF_DISJOINT` is 318 of 318 inside HALLMARK E2F+G2M (E22 section 2d).
   Adjusting a MYC effect for E2F output is close to adjusting MYC for itself.
3. **Protein, or a functional readout.** Both endpoints here are transcript
   ratios and the mouse's are too, but the mouse's licensing came from a
   perturbation design that this cross-sectional cohort cannot imitate.
4. **A stratum-native analysis.** Every reversal here is a Luminal reversal, and
   Luminal is 64 pct of TCGA and 76 pct of SCAN-B. What is different about the
   pooled fit is composition, and nothing here decomposes it.

---

## 9. What was NOT done, listed so it stays a decision

No interaction term, no `PRIME`, no Johnson-Neyman, no mediation. No per-gene
FDR. `beta2` was not scored against any prediction. **`BCL2L1` was not followed
up alone** - the pre-specified object is the pair and it stays the pair. SCAN-B
purity was not imputed. Nothing was read from `myc_human_validation`. The mouse
document was read **read-only** and cited as conclusions, never as values. No
ortholog function called, in either direction.

---

## 10. Where the numbers live

| | |
|---|---|
| script | `scripts/E23_guardian_balance.R` |
| object | `results/guardian_balance.rds` - `$pooled`, `$strata_fits`, `$additivity`, `$f1`, `$f2`, `$succ`, `$verdicts`, `$repro`, `$ties`, `$contamination`, `$settings`, `$rules` |
| table | `outputs/tables/E23_guardian_balance.csv` |
| figures | none |
| the control | `$repro`: 100 of E22's single-gene cells, max \|delta\| = 0; and `$f1`'s `BBC3` control row, which must fire |
| the pre-specification | the script header, fixed before the fits |
| what it follows | `docs/2026-09-07_e22_bbc3_direct_myc_effect.md` |
| the mouse document, read-only, conclusions only | `myc_mouse` `docs/2026-09-02_myc_oxphos_priming_gate_model.md` |
