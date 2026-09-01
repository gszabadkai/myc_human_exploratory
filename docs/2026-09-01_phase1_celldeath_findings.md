# Phase 1 findings - the cell-death axis

2026-09-01. Source: `results/celldeath_axis.rds`, built by `E05`.
The plane it sits on: `docs/2026-08-31_phase1_atlas_findings.md`.

> **EVERYTHING BELOW IS HYPOTHESIS-GENERATING.** Nothing is pre-registered. The
> grid is 73,422 cells. Claims are made only where they hold across both
> cohorts, and each carries what would kill it.

**Read D0 first. It reframes every other number on this page.**

---

## D0. Most of the death axis on this plane is mitochondrial content

Rank the 18 scored death sets by how strongly they correlate with MYC, and
separately by what fraction of their genes are MitoCarta genes. The two orderings
are nearly the same:

```
rho( set's correlation , set's mitochondrial fraction ), over the 18 sets

                    vs MYC        vs OXPHOS subunits
TCGA                 0.603              0.806
SCAN-B               0.692              0.783

excluding the three _MITO strata (15 sets):
TCGA                 0.470              0.688
SCAN-B               0.584              0.645
```

The four sets that clear |z| > 2 against an expression-matched null **in both
cohorts** are exactly the four most mitochondrial ones:

```
                                  z TCGA   z SCAN-B   mitochondrial fraction
CDC_PROSURVIVAL_APOPTOSIS_MITO     3.30      3.71            0.951
CDC_PRODEATH_APOPTOSIS_MITO        2.83      3.01            0.949
TANG_MPT_DRIVEN_NECROSIS           2.52      2.35            0.767
TANG_CUPROPTOSIS                   2.71      2.95            0.556
```

- The `_MITO` strata are MitoCarta-restricted **by construction**. Both sit at
  +0.80 against OXPHOS subunits and their pro-death minus pro-survival contrast
  is **+0.002 / -0.034** - zero. They measure mitochondrial content, not death.
- `TANG_MPT_DRIVEN_NECROSIS` is 23/30 MitoCarta: `ATP5F1A`, `ATP5F1B`, `ATP5PF`,
  `VDAC1/2/3`, `SLC25A4/5/6`, `MCU`, `PPIF`, `SPG7`. Three of those are
  ATP-synthase subunits, i.e. **members of the OXPHOS arm it is being correlated
  against.**
- `TANG_CUPROPTOSIS` is 15/27 MitoCarta and its core is the lipoylated TCA
  machinery: `FDX1`, `LIAS`, `LIPT1`, `DLD`, `DLAT`, `PDHA1`, `PDHB`. Its
  correlation with MYC (+0.514 / +0.427, the highest of any death set) is a
  statement about **mitochondrial metabolism**, not about copper-dependent death.

**So the headline "death programmes track OXPHOS" is, for most of these sets,
a restatement of "mitochondrial gene sets track the mitochondrial programme".**
Every number below is reported against that background, and the one finding that
escapes it is D1.

The expression-matched null itself replicates well - z agrees between cohorts at
Spearman **0.94** across the 18 sets - so it is a reliable filter, not noise.

---

## D1. The pro-death / pro-survival skew - and it is an OXPHOS finding, not a MYC one

The unrestricted apoptosis strata (502 pro-death, 584 pro-survival genes, both
only ~7% mitochondrial, so D0 does not apply) go opposite ways:

```
                            prodeath   prosurvival   contrast
TCGA    raw     vs MYC        -0.029       0.339      -0.368
TCGA    raw     vs OXPHOS     -0.158       0.255      -0.413
SCAN-B  raw     vs MYC        -0.065       0.371      -0.437
SCAN-B  raw     vs OXPHOS     -0.268       0.283      -0.551
```

**MYC-high and OXPHOS-high tumours carry a pro-survival-skewed apoptotic
transcriptome.** Replicated in both cohorts on both axes, and the contrast is
the right quantity to read because both strata come from the same curation and
differ only in annotated direction of effect.

**Then the adjustments split them apart, and this is the important part:**

```
                            contrast vs MYC     contrast vs OXPHOS
TCGA    raw                    -0.368              -0.413
TCGA    proliferation-adj      -0.050              -0.320
TCGA    purity + leukocyte     -0.361              -0.409
SCAN-B  raw                    -0.437              -0.551
SCAN-B  proliferation-adj      -0.125              -0.372
```

Against **MYC** the contrast collapses under proliferation adjustment - to -0.05
and -0.13, i.e. essentially gone. Against **OXPHOS** it survives, losing about a
fifth in TCGA and a third in SCAN-B. Purity and leukocyte fraction cost it
almost nothing.

> **The defensible claim is: tumours with a high OXPHOS programme carry a
> pro-survival-skewed apoptotic transcriptome, and this is not explained by
> proliferation, purity or immune infiltrate. The apparent MYC version of the
> same statement is largely proliferation.**

That is a different claim from the one the plane's F1 supports, and the two
should not be run together. It also sits oddly beside F1: MYC and OXPHOS
correlate at +0.39/+0.43, yet the death skew attaches to one and not the other.

**Caveat that must travel with D1.** The gene-level null test is far more modest
than the GSVA contrast. `CDC_PROSURVIVAL_APOPTOSIS` sits at z = +1.58 / +2.59 and
`CDC_PRODEATH_APOPTOSIS` at z = -1.70 / -1.89. Directionally consistent in both
cohorts, but small. A GSVA score is a relative-rank enrichment and amplifies;
**-0.37 is not an effect size.** The sign and its replication are the result.

**What would falsify D1**
1. A better proliferation control removing the OXPHOS contrast too - it already
   removes the MYC one.
2. The skew being an annotation artefact: if `pro-survival` genes are on average
   more highly expressed or more mitochondrial than `pro-death` genes, the
   contrast would follow from that. The mitochondrial fractions are close (0.070
   vs 0.076) but **expression has not been matched between the two strata**, and
   that is the single cheapest remaining test.
3. Failure at protein level, where BCL2-family behaviour is largely
   post-translational.

---

## D2. Trap 10 vindicated: TANG_FERROPTOSIS does not survive its null

`TANG_FERROPTOSIS` correlates +0.398 (TCGA) and +0.292 (SCAN-B) with MYC as a
GSVA score. Against 2,000 size- and expression-matched draws its genes give
**z = +1.33 and +0.98** - indistinguishable from 418 random genes of the same
expression.

**That correlation should not be reported.** It is what CLAUDE.md trap 10 was
written for, and it is the clearest demonstration in the study that a large gene
set's GSVA correlation can be an artefact of size and expression.

`TANG_APOPTOSIS` (608) is marginal at +1.81 / +2.32.
`TANG_AUTOPHAGY_DEPENDENT_CELL_DEATH` (876) is the only set **below** its null in
both cohorts (-2.24 / -1.53) - its genes correlate with MYC *less* than matched
random genes.

---

## D3. The BCL2 family splits by function, and the split is not mitochondrial content

Gene by gene against MYC activity (FELSHER_61), TCGA / SCAN-B:

```
up with MYC          BID    +0.425 / +0.346      pro-death, activator
                     BAK1   +0.412 / +0.354      pro-death, effector
                     BAX    +0.282 / +0.335      pro-death, effector
                     HRK    +0.276 / +0.148      pro-death, BH3-only
                     BCL2A1 +0.203 / +0.192      pro-survival

flat                 BIK    -0.044 / +0.012
                     BAD    -0.050 / -0.025
                     BCL2L11 (BIM)  -0.100 / -0.059     pro-death, BH3-only
                     BMF    -0.100 / -0.109
                     PMAIP1 (NOXA)  -0.104 / -0.091
                     BBC3   (PUMA)  -0.117 / -0.004
                     BCL2L1 (BCL-xL) -0.124 / +0.008    pro-survival
                     MCL1   +0.014 / -0.015              pro-survival

down with MYC        BCL2L2 -0.211 / -0.144
                     BCL2   -0.369 / -0.288     pro-survival
```

Two things are worth saying and one is worth not saying.

**The split is by function, not by compartment.** 13 of the 15 are MitoCarta
genes, so D0's confound cannot explain why `BAK1` goes up and `BCL2` goes down -
both are mitochondrial. Effectors and the activator `BID` rise; the guardian
`BCL2` falls.

**The BH3-only sensitizers are flat.** `BCL2L11` (BIM), `BBC3` (PUMA) and
`PMAIP1` (NOXA) sit within +/-0.12 of zero in both cohorts, as do `MCL1` and
`BCL2L1`. Whatever MYC does to apoptotic sensitivity in these tumours, **it is
not visible in the transcript levels of the sensitizers**.

> **QUALIFIED 2026-09-01 by `E08`. The MYC column above is inflated by pooling.**
> Split by stratum, every one is weaker than the pooled value in both cohorts:
> `BCL2` runs -0.15 (ERpos), -0.18 (ERneg), -0.11 (Luminal) and **-0.009 (LumA)**
> against the -0.369 quoted here. ER-negative and basal tumours are both MYC-high
> and BCL2-low, so pooling manufactures most of it. `BID` and `BAK1` are inflated
> the same way, more mildly.
>
> **The OXPHOS axis is stratum-stable and is where this claim belongs**: `BAX`
> 0.45-0.58 across every stratum, `BCL2L11` -0.26 to -0.35, `MCL1` -0.21 to
> -0.31. See S1 in `docs/2026-09-01_phase2_estimator_findings.md`.

**What must not be said:** that this explains, supports or contradicts the
validation study. That study tested a specific pre-registered hypothesis about
the `MYC x OXPHOS` interaction on apoptotic priming and found it null. This is an
unregistered marginal correlation of transcript levels in the same cohorts. It
is a different quantity, and CLAUDE.md trap 1 applies with full force.

Note also that the CollecTRI regulon estimator `M_b` gives a much sharper version
of the same pattern (`BAK1` +0.538, `BID` +0.556, `BCL2A1` +0.658, `BCL2`
-0.550) despite being the *weakest* estimator on the OXPHOS plane. Worth
understanding before either number is used.

---

## D4. CICD cannot be scored on the pro-survival side, and the data say so

`CDC_PROSURVIVAL_CICD` is four genes. E05 refuses to score it. Here is why, from
the genes themselves, against MYC (TCGA / SCAN-B):

```
BCL2    -0.369 / -0.288
NFE2L2  -0.118 / -0.043
GPX4    -0.012 /  0.000        (but +0.561 / +0.396 against OXPHOS)
GCH1    +0.135 / +0.194
```

They do not move together on either axis, and `GPX4` is near-zero with MYC while
being one of the strongest OXPHOS correlates in the whole overlay. A four-gene
mean of these would be a number with no referent.

The 13-gene pro-death side is scored (+0.291 / +0.344 with MYC) but is itself
heterogeneous - `BAK1` +0.412, `PPIF` +0.382, `TFRC` +0.339 at one end,
`BNIP3L` -0.394, `BCL2` -0.369, `BECN1` -0.243 at the other - and its
expression-matched null is cohort-discordant (z +1.01 TCGA, +2.11 SCAN-B).
**CICD remains the axis of most interest and the weakest measured (trap 9).**

`MIF` is absent from the TCGA matrix entirely and is missing from TCGA's CICD
numbers.

---

## D5. The product terms carry nothing

All 12 `MYC x OXPHOS` product terms against all 18 death sets, both cohorts, sit
between -0.16 and +0.08, with no structure.

This is **not** evidence about an interaction, and must never be reported as
such - a product of two z-scores is large in the both-low corner as well as the
both-high one, which is why the plan's product term is a poor instrument for the
question it superficially resembles. The validation study tested the interaction
properly and found it null. The correct statement here is simply that the
product term is uninformative and should be dropped from phase 2.

---

## What phase 2 should take from this

1. **D1 is the only death finding that survives D0**, and it is an OXPHOS
   finding. Test it properly: match expression between the pro-death and
   pro-survival strata, then re-run.
2. **Every death set should carry its mitochondrial fraction** as a reported
   covariate. Without it, "death tracks OXPHOS" is unreadable.
3. **The expression-matched null should be standard**, not a check. It agreed
   across cohorts at 0.94 and it killed one headline number.
4. The `BAK1`/`BID`/`BAX` up, `BCL2` down pattern is the most mechanistically
   suggestive thing here and is *not* a content artefact. It deserves a protein-
   level or functional follow-up rather than more correlation.
