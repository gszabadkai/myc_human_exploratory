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

### P4a. A third null, matching the sub-compartment - and the residue grows

Section 3.2, added 2026-09-02. E10 found a ladder by depth into the organelle
(not in MitoCarta -0.100, MOM +0.085, IMS +0.115, MIM +0.443), so 20
mitochondrial genes drawn without regard to compartment are not the machinery's
20, which are **13 MOM, 5 IMS and 2 MIM**. This null reproduces that
distribution exactly, expression-matching within each compartment, still from
the strict pool.

| axis, adjusted | cohort | observed | null mean +/- SD | z | percentile |
|---|---|---|---|---|---|
| OXPHOS | TCGA | 0.453 | 0.276 +/- 0.134 | **1.32** | 90th |
| OXPHOS | SCAN-B | 0.489 | 0.312 +/- 0.132 | **1.34** | 91st |
| MYC | TCGA | 0.187 | 0.218 +/- 0.141 | -0.22 | 39th |
| MYC | SCAN-B | 0.137 | 0.146 +/- 0.147 | -0.07 | 46th |

**The null FELL when the compartment was matched** - 0.29-0.37 under the strict
pool, 0.28-0.31 here - because the machinery is MOM-heavy and MOM has the
lowest median of the three compartments, while MitoCarta at large is
matrix-heavy. So z rises to about **1.3, at the 90th percentile of draws, in
both cohorts and under both adjustments**.

That is the hardest null available without a new annotation source, and the
answer is the same in kind: **consistent, in the right direction, and not
separable.** P4 stands - but the residue does not disappear when composition is
matched more tightly, it grows a little, which is the opposite of what a pure
composition artefact should do. **This is the single number most worth
attacking next**, and the way to attack it is the sub-compartment ladder over
more than 2 MIM genes.

## P4b. Inside a single subtype - and pooling was UNDERSTATING the MYC split

Section 4.2, added 2026-09-02, because P1 was measured on all samples and E10
fig9 had just shown that for the BCL2-family ratios the pooled MYC value sits
outside the range of both its compartments for 27 of 39 ratios in TCGA.

**The spread contrast holds inside both compartments**, which is what P1 needed:

| cohort | stratum | n | SD MYC | SD OXPHOS | ratio |
|---|---|---|---|---|---|
| TCGA | Luminal | 696 | 0.188 | 0.324 | 1.7x |
| TCGA | Basal | 171 | 0.238 | 0.357 | 1.5x |
| SCAN-B | Luminal | 2,436 | 0.160 | 0.254 | 1.6x |
| SCAN-B | Basal | 317 | 0.162 | 0.211 | 1.3x |

**The localisation split does not behave like the spread, and the direction is
the surprise.** Pooling *understates* the MYC split:

| cohort | axis | pooled | Luminal | Basal | pooled outside both? |
|---|---|---|---|---|---|
| TCGA | MYC | 0.187 | 0.298 | 0.323 | **yes, below both** |
| SCAN-B | MYC | 0.137 | 0.295 | 0.194 | **yes, below both** |
| TCGA | OXPHOS | 0.453 | 0.464 | 0.349 | no |
| SCAN-B | OXPHOS | 0.489 | 0.442 | 0.367 | yes, above both |

So the MYC/OXPHOS gap on the split narrows from about 0.30 pooled to about 0.17
inside Luminal in TCGA. **Within a subtype the machinery is more
MYC-organised than the pooled number says**, roughly doubling from 0.14-0.19 to
0.29-0.32 in three of four cells. This is the mirror image of D3/S1: there
pooling inflated a MYC correlation, here it suppresses one.

Against sub-compartment-matched nulls computed within each stratum, everything
is still inside: OXPHOS z 0.51-1.42, MYC z 0.12-0.70.

**This softens P3.** MYC's flatness against the machinery is partly a pooling
artefact of the opposite sign to the one this study has been guarding against.
The z against composition is unchanged in its conclusion, but the raw numbers
are not what the all-samples reading suggested.

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

- **P1** SURVIVED that test on 2026-09-02 - see P4b. The spread contrast holds
  inside both compartments at 1.3x to 1.7x. It would still die if it
  disappeared inside LumA alone, which has not been run.
- **P4** got that null on 2026-09-02 - see P4a - and the residue grew rather
  than vanished, to z ~1.3 at the 90th percentile in all four cells. It would
  now be overturned by a curated non-apoptotic mitochondrial comparator set
  rather than a random draw, or by widening the MIM stratum past 2 genes.
- **P5**'s residue dies if the ~1 SD excess does not grow when the null is made
  stricter still, or if it disappears under a different MYC estimator. The panel
  is scored; this is a filter away.

## What was not done

- ~~No strata.~~ Added 2026-09-02 as section 4.2; see P4b. LumA and LumB
  separately have still not been run, and P4b's Basal numbers rest on 171 TCGA
  samples.
- **`PROLIF_STD` agrees with `PROLIF_DISJOINT` throughout** (differences in the
  third decimal), so the covariate choice does not matter here. It did in E1.
- **Purity is TCGA-only** (n = 1,007). SCAN-B has no purity estimate and it is
  never imputed, so that row cannot replicate - trap 2.
- `TP53` and `BIRC5` are in the proliferation covariate and partly adjusted for
  themselves. Every summary is reported with and without them; dropping them
  changes the SD by at most 0.017.
