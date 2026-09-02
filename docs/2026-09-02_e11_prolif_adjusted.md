---
date: 2026-09-02
script: scripts/E11_prolif_adjusted_machinery.R
status: WRITTEN AND VERIFIED, NOT YET RUN BY THE AUTHOR
posture: EXPLORATORY - nothing here is pre-registered
supersedes: the interpretation of R1 in docs/2026-09-02_e10_machinery_and_priming.md
---

# E11 - after proliferation, does OXPHOS still organise the apoptotic
# machinery, and does MYC?

The question, as put on 2026-09-02: *"after correcting to proliferation, OXPHOS
correlates with the apoptotic machinery, while MYC does not."*

**The comparison is real and it survives every correction. The word
"correlates" is where it needs restating, in two places.** Numbers from a dry
run against both cohorts; nothing is written into the repo until the script is
sourced.

---

## P1. The comparison holds, and proliferation is not what makes it

Spread of the 44 per-gene correlations - the claim as one number:

| cohort | axis | raw SD | adjusted SD | frac \|rho\| > 0.2, adjusted |
|---|---|---|---|---|
| TCGA | OXPHOS | 0.301 | **0.313** | 0.66 |
| TCGA | MYC | 0.227 | **0.191** | 0.34 |
| SCAN-B | OXPHOS | 0.243 | **0.245** | 0.50 |
| SCAN-B | MYC | 0.190 | **0.157** | 0.25 |

Adjusting for proliferation (`PROLIF_DISJOINT`, 318 genes, partial Spearman,
applied identically to both axes) leaves the OXPHOS spread untouched and
narrows the MYC spread by about 16%. The localisation ordering behaves the same
way: OXPHOS 0.460 -> 0.492 (TCGA) and 0.539 -> 0.535 (SCAN-B); MYC 0.161 ->
0.216 and 0.267 -> 0.189.

It survives the second correction as well - `FELSHER__PROLIFSTRIP` and
`FELSHER__BOTHSTRIP`, which remove the E2F/G2M genes from the estimator rather
than partialling them out - and in TCGA it survives purity and leukocyte
fraction on top of proliferation (n = 1,007; OXPHOS split 0.406, MYC 0.205).
**E1 warned the two corrections disagree in sign; here they agree**, which is
the first time in this study they have.

**So proliferation is not the story.** The MYC column was already the weaker one
before any adjustment, and adjustment sharpens rather than creates the contrast.

## P2. THE CONTROL PASSES - adjusted MYC has not been emptied

The obvious boring explanation: a MYC signature is partly a proliferation
signature, so partialling proliferation out could leave a score that correlates
with nothing. Mean \|rho\| above an expression-matched null (2,000 draws, 20
ventiles), as z:

| set | MYC raw | MYC adjusted | OXPHOS raw | OXPHOS adjusted |
|---|---|---|---|---|
| mitoribosome (83) | 8.4 / 6.9 | **6.5 / 8.7** | 11.9 / 12.0 | 10.0 / 10.7 |
| OXPHOS subunits (89) | 3.4 / 2.8 | 0.7 / 4.1 | 15.6 / 18.3 * | 15.9 / 18.3 * |
| apoptotic machinery (44) | -0.8 / -1.0 | **0.1 / 0.7** | 0.0 / 1.1 | 1.3 / 1.5 |

(SCAN-B / TCGA. `*` is the OXPHOS arm against the OXPHOS score - its own genes.)

The mitoribosome control stays at z = 6.5-8.7 on adjusted MYC. **The adjustment
did not silence the estimator**, so the machinery's flatness on that axis is
about the machinery.

## P3. FIRST RESTATEMENT - MYC never correlated with the machinery

Raw z for the machinery on MYC is **-0.8 / -1.0** - at or slightly below an
expression-matched null before any correction is applied, and it is the same
under both stripped estimators (-0.7 to -0.9).

So the finding is not *"proliferation adjustment removes MYC's correlation with
the machinery"*. It is **"MYC never had one above background, and the
correction is not what removes it."** That is a cleaner claim and a stronger
one - it does not depend on the covariate being right.

## P4. SECOND RESTATEMENT, AND IT IS THE ONE THAT BITES - the OXPHOS side is
## mitochondrial composition, not an apoptotic programme

Two nulls E10 did not run.

**The machinery is not special on OXPHOS either.** Its mean \|rho\| with OXPHOS
is 0.274 (TCGA, adjusted) against an expression-matched null that sits at
essentially the same place: z = 1.3 / 1.5. An OXPHOS score correlates with a
large fraction of the transcriptome, and 44 genes of this expression profile
reach that value by construction.

**And the localisation split is what composition gives.** Holding the
composition fixed - 20 MitoCarta and 24 non-MitoCarta genes, expression-matched
within each half:

| axis | pool | observed | null mean +/- SD | z |
|---|---|---|---|---|
| OXPHOS, adjusted, TCGA | all MitoCarta | 0.453 | 0.397 +/- 0.130 | 0.43 |
| OXPHOS, adjusted, TCGA | no OXPHOS/mitoribo | 0.453 | 0.289 +/- 0.142 | **1.16** |
| OXPHOS, adjusted, SCAN-B | all MitoCarta | 0.489 | 0.452 +/- 0.117 | 0.32 |
| OXPHOS, adjusted, SCAN-B | no OXPHOS/mitoribo | 0.489 | 0.362 +/- 0.135 | **0.94** |
| MYC, adjusted, TCGA | no OXPHOS/mitoribo | 0.187 | 0.248 +/- 0.143 | -0.43 |
| MYC, adjusted, SCAN-B | no OXPHOS/mitoribo | 0.137 | 0.230 +/- 0.141 | -0.66 |

The strict pool removes every OXPHOS and mitoribosome gene from the null,
because those correlate with an OXPHOS score by definition while the
machinery's own 20 mitochondrial genes are BCL2-family members and caspases.
Even so, the observed split sits about **one SD** above it - in the same
direction in all four cohort x adjustment cells, which is worth noting, but not
separable from composition.

**So E10 R1's 0.453 is a fact about mitochondrial localisation and expression,
not about apoptosis.** The E10 note now carries this correction at the top.

## P5. What is left, stated so it can be falsified

**The apoptotic machinery is arranged along the OXPHOS axis and is not arranged
along MYC activity, and neither pattern is explained by proliferation.** The
OXPHOS arrangement is the arrangement mitochondrial localisation produces; the
MYC non-arrangement is *below* what composition predicts, in both cohorts.

Both halves replicate: the 44 per-gene values correlate 0.88-0.92 between TCGA
and SCAN-B under every adjustment.

That is a defensible statement and it is narrower than the one it started as.
The interesting residue is the MYC half: a set of mitochondrial genes that
tracks the OXPHOS axis about as expected but tracks MYC *less* than its
composition predicts. If that is real it is a dissociation, not an absence.

### What would falsify it

- **P1** dies if the contrast disappears inside the luminal compartment. E10
  fig9/fig10 already show the MYC column is largely between-subtype for the
  BCL2-family genes; the same split has not been run on the 44.
- **P4** would be overturned by a null built on a curated non-apoptotic
  mitochondrial gene set with matched sub-compartment composition. The MOM /
  IMS / MIM distribution of the machinery's 20 is not matched by the current
  draw, only MitoCarta membership is.
- **P5**'s residue dies if the ~1 SD excess does not grow when the null is made
  stricter still, or if it disappears under a different MYC estimator. The panel
  is scored; this is a filter away.

## What was not done

- **No strata.** All-samples only, both cohorts. Given E10's compartment result
  this is the first thing to add.
- **`PROLIF_STD` agrees with `PROLIF_DISJOINT` throughout** (differences in the
  third decimal), so the covariate choice does not matter here. It did in E1.
- **Purity is TCGA-only** (n = 1,007). SCAN-B has no purity estimate and it is
  never imputed, so that row cannot replicate - trap 2.
- `TP53` and `BIRC5` are in the proliferation covariate and partly adjusted for
  themselves. Every summary is reported with and without them; dropping them
  changes the SD by at most 0.017.
