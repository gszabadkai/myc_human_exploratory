---
date: 2026-09-02
script: scripts/E09_correlation_measures.R
status: RUN 2026-09-02 10:42 - results/correlation_measures.rds + 2 figures
posture: EXPLORATORY - nothing here is pre-registered
---

# E09 - is Spearman hiding anything?

The author's third phase-2 question: every correlation in this study is a
Spearman. Is there a specific advantage to that, or would another measure
reveal something else?

**Answer: it is the right default, the choice is immaterial on GSVA, and where
the measures disagree Spearman is demonstrably the correct one.** 220 pairs, 5
MYC estimators x 22 partners x 2 cohorts.

---

## C1. The four measures agree almost completely

| pair | Spearman of the two, over 220 pairs | max abs difference |
|---|---|---|
| spearman vs bicor | **0.999** | 0.033 |
| spearman vs kendall | 0.998 | 0.204 |
| spearman vs pearson | 0.996 | 0.093 |
| pearson vs bicor | 0.995 | 0.078 |

**The Kendall row is a scale difference, not a disagreement.** Over the 108
pairs with `|rho| > 0.2` the ratio `kendall / spearman` is **0.681 +/- 0.033**,
against the bivariate-normal expectation of about 0.68. Kendall reproduces
Spearman's ordering exactly and reports it on a compressed scale. It is a
consistency check that passed and carries no independent information.

## C2. Every material disagreement is a heavy tail, and Spearman wins it

The reading rule was written into the script before it ran: *pearson departs
from spearman while bicor sits WITH spearman -> heavy-tail artefact; bicor sits
WITH pearson -> magnitude carries real information ranks discard.*

**In all 12 of the largest Spearman-Pearson gaps, bicor sits with Spearman.**
Not most - all of them.

| cohort | estimator | partner | spearman | pearson | bicor | gap |
|---|---|---|---|---|---|---|
| TCGA | `HALLMARK_MYC_TARGETS_V1__FULL` | mitopps::OXPHOS subunits | 0.531 | 0.438 | **0.503** | -0.093 |
| SCAN-B | `M_b__MITOSTRIP` | mitopps::Mitochondrial ribosome | 0.262 | 0.347 | **0.287** | +0.085 |
| TCGA | `HALLMARK_MYC_TARGETS_V1__FULL` | mitopps::Fatty acid oxidation | -0.421 | -0.344 | **-0.422** | +0.077 |
| TCGA | `FELSHER__MITOSTRIP` | mitopps::OXPHOS subunits | 0.306 | 0.241 | **0.299** | -0.065 |

Across all 220 pairs bicor sides with Spearman 144 times against Pearson 76,
but the 76 sit in the regime where the gap is under 0.01 and the question does
not arise.

## C3. The disagreement is an INSTRUMENT effect, and it names mitoPPS

Mean absolute departure from Spearman, by instrument:

| instrument | mean \|pearson - spearman\| | max | mean \|bicor - spearman\| |
|---|---|---|---|
| **mitopps** | **0.029** | **0.093** | 0.009 |
| content | 0.020 | 0.069 | 0.008 |
| zmean | 0.017 | 0.076 | 0.009 |
| gsva | **0.009** | 0.038 | 0.006 |

Bicor stays within 0.006-0.009 of Spearman on every instrument; only Pearson
moves, and it moves most on mitoPPS.

**This turns CLAUDE.md trap 6 from an assertion into a measurement.** mitoPPS is
a linear composition statistic with a long right tail, and Pearson reads that
tail. On GSVA the choice of measure is worth less than 0.01 and is genuinely
immaterial; on mitoPPS it is worth up to 0.093 and Spearman is the one to trust.

## C4. Nothing in the atlas is non-monotone

Largest gain in R-squared from a natural spline over a straight line, on ranks,
across all 220 pairs: **0.051**. The headline `FELSHER__MITOSTRIP` vs
`gsva::OXPHOS subunits` pair returns **0.012**.

The decile profile of that pair is a straight rise with no plateau at either
end:

| MYC decile | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| TCGA | -0.369 | -0.216 | -0.113 | -0.119 | 0.021 | 0.002 | 0.138 | 0.052 | 0.080 | 0.215 |
| SCAN-B | -0.362 | -0.221 | -0.147 | -0.073 | -0.029 | 0.013 | 0.058 | 0.144 | 0.118 | 0.175 |

**There is no threshold and no saturation.** MYC's biology is dose-dependent
and this was a real blind spot to check rather than a theoretical one; the check
came back empty.

The largest spline gains that do exist are all `M_b__MITOSTRIP` against mitoPPS
and zmean arms (0.028-0.051). E3 already established `M_b` is the activated and
weaker half of its own regulon, so a slightly odd shape there is consistent
with what is known about that estimator and is not large enough to act on.

---

## What this licenses, and what it does not

- **The study's Spearman default is vindicated and needs no hedging.** Report it
  without a sensitivity caveat.
- **Do not report a Pearson against a mitoPPS score anywhere.** C3 is the
  measured reason, and E10's Pearson panels deliberately exclude mitoPPS for
  it.
- **Do not claim the measures "agree" without naming the instrument.** They
  agree to 0.009 on GSVA and to 0.029 on mitoPPS, and only the first of those
  is negligible.
- **This does not license reporting Pearson as a confirmation.** 220 pairs
  correlating at 0.996 are not 220 independent confirmations of anything; they
  are one statement about the shape of the data.

## What would falsify it

- C4 dies if a non-monotone relationship exists in a stratum but not pooled. The
  spline probe was run on all samples only, and a threshold present only in
  Basal would be invisible here.
- C2/C3 die if the mitoPPS tail is a property of this pathway universe rather
  than of composition statistics generally - testable by rebuilding mitoPPS on a
  different arm set.
