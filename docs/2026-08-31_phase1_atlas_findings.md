# Phase 1 findings - the correlation atlas

2026-08-31. Source: `results/correlation_atlas.rds`, built by `E03` from `E02`.

> **EVERYTHING BELOW IS HYPOTHESIS-GENERATING.** Nothing in this repo is
> pre-registered. The atlas is 68,255 cells and no single one of them is a
> result. What follows is claimed only where it holds across **both cohorts**,
> **multiple instruments** and **the estimator panel** - and each claim carries
> the observation that would kill it. Those falsifiers are written here *before*
> `E04` and `E05` are built, so the exploratory phase can hand a real hypothesis
> to a confirmatory one.

The atlas reproduces the validation study's three published anchors to within
2-3e-4 when they are read back out of the finished table, which is what
establishes that it is indexed correctly and not merely computed correctly.

---

## F1. The MYC-OXPHOS correlation is real, replicated, and not merely proliferation

**What holds.** Against `OXPHOS subunits`, every one of the 18 MYC signatures is
positive in both cohorts on all four instruments. Raw GSVA spans 0.049 to 0.670
(SCAN-B) and 0.181 to 0.652 (TCGA). Across the 18 arms the two cohorts agree at
Spearman 0.92 (GSVA) and 0.94 (zmean).

It is **not** proliferation wearing MYC's name:

| | SCAN-B | TCGA |
|---|---|---|
| signatures with prolif-adjusted rho >= 0.2 | 15 / 18 | 17 / 18 |
| median prolif-adjusted rho | 0.323 | 0.452 |
| `MYC_UP.V1_UP` (1.5% entangled), raw | 0.502 | 0.552 |
| `MYC_UP.V1_UP`, proliferation-adjusted | 0.401 | 0.515 |
| `MYC_UP.V1_UP`, purity + leukocyte-adjusted | - | 0.538 |

The least proliferation-entangled signature in the compendium gives one of the
**highest** correlations, and loses only 0.04-0.10 when proliferation is
projected out. In TCGA the purity and leukocyte adjustment costs almost nothing
(FELSHER_61 0.388 -> 0.376).

**The honest complication, which must be reported with it.**
`slope_vs_entanglement` - the Spearman, across the 18 signatures, between a
signature's proliferation content and the rho it produces - is **positive in all
20 cohort x instrument x adjustment combinations**, +0.18 to +0.54. On average,
more entangled signatures do correlate more.

But the relationship is **not monotone and the extremes contradict it**:
`rho_least` (1.5% entangled) exceeds `rho_most` (47.6%) in **7 of 8** raw
combinations. The positive slope is carried by two weak signatures in the
*low*-entanglement half - `ELLWOOD_MYC_TARGETS_UP` (13 genes, 0.049 / 0.181) and
`ALFANO_MYC_TARGETS` (0.160 / 0.320) - not by a gradient. And it is not set
size: rho-vs-size is 0.05-0.18 where rho-vs-entanglement is 0.30-0.54.

So: **between-signature spread (0.05 to 0.70) is far larger than the
entanglement trend, and entanglement does not explain the correlation away.**
CLAUDE.md trap 3 stands - never quote one signature - but the panel as a whole
survives.

**What would falsify F1**
1. A proliferation estimator better than `PROLIF_DISJOINT` (e.g. a cell-cycle
   phase deconvolution, or Ki67 where measured) that removes the correlation for
   `MYC_UP.V1_UP` as well as for the entangled sets.
2. `MYC_UP.V1_UP` turning out to be entangled with something else that tracks
   OXPHOS - it is a perturbation-derived set and its 198 genes have not been
   audited against, say, mitochondrial mass or ribosome biogenesis.
3. The correlation vanishing in a cohort scored on a different platform
   (METABRIC microarray would be the cheap test, and is out of scope for
   phase 1 by decision, not by oversight).

---

## F2. MYC mRNA carries no OXPHOS signal at all. The gap is the result.

`log2(MYC)` against `OXPHOS subunits`: **-0.011** (SCAN-B), **-0.032** (TCGA).
Against MYC *activity* on the same arm in the same samples: **+0.427** and
**+0.389**.

Across all 18 arms, `log2(MYC)` never exceeds |rho| = 0.12 in SCAN-B or 0.21 in
TCGA, and its largest value is on `Folate and 1-C`, not on any OXPHOS arm. The
null is not specific to OXPHOS; the transcript simply does not track the
mitochondrial programme.

This replicates CLAUDE.md trap 4 in a second cohort. **Anyone plotting MYC
expression against OXPHOS sees nothing, and would conclude the opposite of what
the activity signatures say.**

**What would falsify F2**: MYC protein or a MYC-amplification-stratified
analysis showing that the mRNA is informative once its dynamic range is
accounted for. Partial support against this already exists in the atlas -
`M_c_call`, the copy-number call, gives +0.137 on `OXPHOS subunits` where the
mRNA gives -0.032, so the *locus* carries some signal that the *transcript level*
does not.

---

## F3. The two genomes are regulated differently - the sharpest finding, and it replicates to two decimal places

Nuclear-encoded and mtDNA-encoded OXPHOS do not behave alike, and the split is
almost identical in the two cohorts:

| FELSHER_61, raw, stratum all | content | gsva | mitopps | zmean |
|---|---|---|---|---|
| `OXPHOS subunits` SCAN-B / TCGA | 0.433 / 0.434 | 0.427 / 0.388 | 0.305 / 0.306 | 0.425 / 0.371 |
| `mtDNA-encoded OXPHOS` SCAN-B / TCGA | 0.043 / 0.056 | 0.061 / 0.068 | **-0.081 / -0.088** | 0.026 / 0.028 |

On the composition instrument the sign **flips**. Across 19 estimators, mitoPPS
returns a negative mtDNA correlation in 16/19 (SCAN-B) and 13/19 (TCGA); the
other three instruments return essentially none. Read as mitoPPS is built to be
read (trap 6): **MYC-high tumours de-prioritise the mtDNA-encoded programme
relative to the rest of the mitochondrion**, while raising the nuclear-encoded
one.

**And the 13 genes do not move together.** Against FELSHER_61, raw:

```
MT-CO2  +0.240 / +0.243      <- positive under 20/20 and 21/21 estimators
MT-ND2  +0.080 / +0.124
MT-ND4  +0.059 / +0.101
...
MT-ND1  -0.041 / -0.057
MT-CO1  -0.007 / -0.073
MT-ND6  -0.045 / -0.154
MT-ND5  -0.160 / -0.143      <- negative under 19/20 and 20/21 estimators
```

`MT-CO2` agrees between cohorts to **0.003**. `MT-ND5` to 0.017.

**Two confounders tested and cleared.**

- *Purity and infiltrate.* In TCGA the split **survives** purity + leukocyte
  adjustment: `MT-CO2` 0.243 -> 0.230, `MT-ND5` -0.143 -> -0.162, `MT-ND6`
  -0.154 -> -0.160. Not stromal or adipose contamination.
- *3' bias and polyadenylation.* There is **no positional gradient**. Ordered
  along the heavy strand, the maximum is at position 4 (`MT-CO2`) and position 1
  (`MT-ND1`) is already negative. `MT-ND6` - the only light-strand, only
  non-polyadenylated protein gene, and therefore the obvious library-chemistry
  suspect - is *not* the extreme in SCAN-B (-0.045).
- *Not tested, and it matters:* proliferation adjustment removes `MT-ND6`'s
  negative almost entirely (-0.154 -> -0.008 in TCGA) but leaves `MT-ND5`
  (-0.143 -> -0.103) and `MT-CO2` (+0.243 -> +0.166). So the ND6 signal and the
  ND5/CO2 signal are not the same phenomenon.

**Why this is interesting rather than merely odd.** `MT-CO2` and `MT-ND5` are
both on the **heavy-strand polycistron** - one primary transcript. A stable,
cross-cohort, opposite-signed correlation between two genes of the same
transcript cannot be explained by mtDNA copy number or by heavy-strand
transcription rate. It has to be post-transcriptional, or an artefact.

**What would falsify F3**
1. **Copy number.** If MYC-high tumours simply carry fewer mitochondrial
   genomes, all 13 should move together. They do not - `MT-CO2` is +0.24 while
   `MT-ND5` is -0.15 - so uniform copy number is already excluded. The remaining
   test is a per-sample mtDNA-content proxy as a covariate: if adjusting for it
   collapses the *spread* among the 13, the finding is a copy-number artefact
   with gene-specific quantification noise on top.
2. **Quantification.** `MT-ND5` and `MT-ND6` overlap on opposite strands, so
   multi-mapping and strandedness could produce exactly this. Test: re-quantify
   the 13 from a stranded pipeline, or check whether the `MT-ND5`/`MT-ND6`
   correlation with each other is anomalous relative to the other 11.
3. **The mitoPPS sign flip is a composition artefact.** mitoPPS is a ratio
   against the rest of the mitochondrial programme. If the nuclear programme
   simply rises more, the mtDNA arm must fall in composition terms with no
   biology at all. **This is the most likely benign explanation and it is not
   yet excluded.** Test: hold the mitoPPS denominator fixed - score the mtDNA
   arm against a mitochondrial universe with the OXPHOS subunits removed - and
   see whether the negative survives.
4. Failure to reproduce in a third cohort.

Falsifier 3 is the one to run first. It is cheap and it is the difference
between a finding and an arithmetic identity.

---

## F4. There are three instruments here, not four

`gsva` and `zmean` agree at median Spearman **0.994** across the 18 arms
(minimum 0.953 over 21 estimators). They are not independent measurements.

Every other pair sits at 0.87-0.90 median but with minima near 0.50, so
disagreement is estimator-specific rather than uniform. `mitopps` is the
instrument that replicates *least* across cohorts (239 of 360 arm-cells
"replicated" against `content`'s 287) and it is also the only one that changes
the sign of the mtDNA arm - i.e. exactly where it disagrees is where it is
designed to say something the others cannot.

CLAUDE.md trap 5 quoted 0.24-0.94 for GSVA-vs-mitoPPS agreement of the **scores**.
This is a different quantity - agreement of the **answers**, arm ordering against
the same estimator - and it is higher. Both belong on the figure and they must
not be conflated.

---

## F5. Present in every subtype, weakest in ER-negative

`FELSHER_61` / `OXPHOS subunits` / GSVA / raw, by stratum:

```
           SCAN-B   TCGA
all         0.427  0.388
ERpos       0.447  0.408
ERneg       0.242  0.337
LumA        0.368  0.418
LumB        0.384  0.272
HER2        0.359  0.290   (TCGA n = 78, CI 0.067 to 0.485)
Basal       0.381  0.443
Normal      0.391  0.273   (TCGA n = 36, CI -0.068 to 0.556)
```

Median across the whole 18-signature panel: ERpos 0.479 / 0.512, ERneg 0.315 /
0.416. **The correlation is not carried by one subtype**, which is the main
thing worth knowing, and ER-negative is consistently the weakest stratum in both
cohorts. TCGA's HER2 and Normal strata are too small to interpret and are shown
with their intervals for that reason.

---

## F6. The mitochondrial ribosome beats OXPHOS, in both cohorts

`Mitochondrial ribosome` 0.617 / 0.590 against `OXPHOS subunits` 0.427 / 0.388,
and it is the top arm of 18 in both cohorts. Consistent with MYC's canonical
ribosome-biogenesis role extending to the mitochondrial ribosome - and it sits
awkwardly beside F3, because the mitoribosome exists to translate the very
mtDNA-encoded transcripts whose programme is de-prioritised. **That tension is
the most interesting thing in the atlas and phase 2 should go at it.**

## F7. Fatty acid oxidation is the only consistently negative arm

-0.140 / -0.121, and the only arm whose entanglement slope is negative (-0.31
TCGA, -0.14 SCAN-B): the *more* proliferation-entangled the signature, the more
negative the FAO correlation. Read with trap 2 - adipose is FAO-high and breast
is the worst TCGA tissue for infiltrate - but note it replicates in SCAN-B,
where there is no purity estimate to blame it on.

---

## What is NOT claimed

- **No interaction is claimed anywhere.** The validation study found the
  `MYC x OXPHOS` interaction on apoptotic priming null. Everything above is
  correlation. CLAUDE.md trap 1: these are different questions and neither
  confirms nor contradicts the other.
- **No mitoPPS value is compared across cohorts** - only patterns (trap 6).
- **No p-value appears above.** The atlas carries one per cell so a cell can be
  described; with 68,255 cells it cannot select one.
- Nothing about the death axis. That is `E05`.
