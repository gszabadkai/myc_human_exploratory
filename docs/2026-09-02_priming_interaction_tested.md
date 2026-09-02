---
date: 2026-09-02
status: a tested hypothesis and its falsifiers. NOT a result note.
purpose: stop the MYC x OXPHOS interaction question being re-asked without its
         controls, and correct the localisation over-claim
---

# Two corrections and one hypothesis that did not survive

## C1. "Outer mitochondrial membrane" was an over-claim. It is CO-REGULATION

MitoCarta is a **proteome catalogue**: membership means a protein has been
detected at mitochondria, not that it sits there constitutively. Half the
BCL2 family translocates - BAX is cytosolic until activated, BID must be cleaved
to tBID first, BAD is held by 14-3-3 when phosphorylated. A static localisation
reading of these genes is not defensible.

**And this is a transcriptomic study.** A transcript has no idea where its
protein ends up. So what MitoCarta membership actually marks, at this level of
measurement, is **membership of the nuclear-encoded mitochondrial regulon** -
the set of transcripts co-regulated with mitochondrial biogenesis.

THE REFRAME IS STRONGER THAN WHAT IT REPLACES. It explains why the composition
null came out flat: of course the machinery is not special, because any
transcript in that regulon tracks an OXPHOS score. The finding is about
transcriptional co-regulation, which is a claim this data type can support,
rather than about protein localisation, which it cannot.

Wherever a figure or note says "mitochondrial genes", read "transcripts of the
nuclear-encoded mitochondrial regulon". The sub-compartment ladder in E10
(MOM/IMS/MIM) should be read the same way and rests on 13/5/2 genes.

## C2. The ratios keep MYC-specific signal that single genes lose, and that is
## a property of ratios rather than of priming

| | MYC, prolif removed | + OXPHOS removed | change |
|---|---|---|---|
| the 44 genes, localisation split | 0.187 | **-0.043** | destroyed |
| the 35 ratios, mean \|rho\| | 0.142 | **0.113** | -20% |

The asymmetry is real and the explanation is mechanical. A ratio subtracts two
transcripts. Whatever the two share cancels - and what they mostly share is the
OXPHOS-borrowed component, because that component is a scaled copy of the OXPHOS
column and hits pro and anti genes alike. What does NOT cancel is MYC's own
component, which E12 showed is the larger of the two and is unsorted.

**So a ratio is an enrichment device for the non-shared part of an association.**
That is why MYC survives in the ratios and not in the genes, and it is not
evidence of a priming-specific mechanism.

## H1. "MYC x OXPHOS interaction sets the priming threshold" - TESTED, NOT
## SUPPORTED

The hypothesis: OXPHOS sets the level, and the MYC x OXPHOS interaction sets
where the balance tips. Worth testing, and tested four ways.

**The encouraging first look.** Fitting `ratio ~ prolif + MYC + OXPHOS +
MYC:OXPHOS` over the 35 ratios: 9 of 35 interactions clear \|t\| > 2 in TCGA and
11 of 35 in SCAN-B, the 35 coefficients correlate **r = 0.49** across cohorts,
and five replicate with the same sign and \|t\| > 2 in both - `BID/BCL2`,
`BMF/BCL2`, `BBC3/BCL2`, `BAD/BCL2` (all positive) and `PMAIP1/MCL1` (negative).
Four of the five share BCL2 as denominator.

**Then the falsifiers.**

| test | result |
|---|---|
| **F1 curvature.** A product term absorbs nonlinearity in either main effect. Refit with natural splines on both. | **FAILS.** In TCGA the linear and spline interaction coefficients correlate **0.125** - which ratios carry it changes almost completely. Three of the five named above collapse: `BID/BCL2` 0.081 to -0.003, `BAD/BCL2` 0.059 to 0.012, `BMF/BCL2` 0.065 to 0.019 |
| **F2 estimator.** Repeat with `MYC_UP.V1_UP`, 1.5% proliferation entanglement. | **FAILS.** Cross-cohort replication **reverses**, r = -0.33 on the linear term and -0.42 on ranks |
| **F3 tertiles.** The same hypothesis without a model: split tumours by OXPHOS tertile, ask whether MYC's association with the ratio strengthens. | **FAILS.** Mean rho by tertile in TCGA is 0.010 / 0.003 / **0.031** - non-monotone, dipping in the middle. SCAN-B moves the other way (mean high-minus-low -0.008 against TCGA +0.021) |
| **F4 prior.** CLAUDE.md trap 1: the sibling PRE-REGISTERED study tested MYC x OXPHOS on apoptotic priming and found it null. | Different endpoint - functional priming there, transcript ratios here - so this neither confirms nor contradicts it. **Nothing here overturns that null and nothing here may be written as if it did.** |

**Verdict: the apparent interaction is mostly curvature in the main effects.**
The TCGA mid-tertile dip is exactly what a product term misreads as an
interaction, and it is why F1 breaks it.

### What would test it properly

A transcript ratio cannot see priming. Priming is a post-translational,
protein-interaction property - which BH3 profiling measures and RNA does not.
The tertile design is the right shape, but it needs that endpoint, or a
perturbation that moves OXPHOS and asks whether the MYC-to-priming slope
changes. **That is the mouse arm's question, not this one's.**

## What DID survive, and it is narrow

Two priming genes have a MYC association that survives conditioning on OXPHOS,
in both cohorts, under all three MYC estimators (`FELSHER__MITOSTRIP`,
`MYC_UP.V1_UP__FULL`, `HALLMARK_MYC_TARGETS_V1__FULL`):

| gene | MYC \| prolif + OXPHOS, across 2 cohorts x 3 estimators |
|---|---|
| **BID** | +0.12 to +0.28, positive in all six |
| **PMAIP1 (NOXA)** | -0.05 to -0.21, negative in all six |

Neither is a member of any of the three estimators, so this is not self-overlap.

**`BCL2` does not survive.** Its MYC association looks strong on the reference
estimator (-0.22 TCGA) but falls to -0.06 and -0.01 under the other two after
conditioning. The reading that "MYC represses BCL2, and every pro/BCL2 ratio
inherits it" is estimator-dependent and should not be written.

And note that BID and NOXA are both pro-apoptotic BH3-only proteins moving in
OPPOSITE directions with MYC. Whatever this is, it is not a coherent shift in
the priming balance - it is two genes.
