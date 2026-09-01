# Phase 2 findings - the anatomy of the MYC estimator panel

2026-09-01. Source: `results/estimator_anatomy.rds`, built by `E06` on the
relabelled panel. Answers two of the handoff's section 5 questions.

> EVERYTHING BELOW IS HYPOTHESIS-GENERATING. Nothing here is pre-registered.
> F1's own resolution is in `docs/2026-08-31_phase1_atlas_findings.md`; this
> note is about the estimators rather than about the correlation.

---

## E1. Adjusting for proliferation and removing proliferation genes are different operations, and they disagree in sign

`E06` prints the two side by side for the first time. TCGA, OXPHOS subunits,
GSVA:

```
                          FULL raw   FULL prolif-ADJUSTED   PROLIFSTRIP raw
FELSHER                     0.404          0.344                0.416
dang_myc_core_genes         0.340          0.256                0.373
MYC_UP.V1_UP                0.552          0.515                0.552
MUHAR_MYC_SIGNATURE         0.469          0.422                0.462
PID_MYC_ACTIV_PATHWAY       0.307          0.205                0.272
HALLMARK_MYC_TARGETS_V1     0.603          0.718                0.671
DANG_MYC_TARGETS_UP         0.652          0.715                0.688
```

Look at `FELSHER`. It is **the one signature `PROLIF_DISJOINT` is genuinely
disjoint from** - the covariate was built by deleting exactly its 9
proliferation genes. Adjusting still *lowers* its correlation by 0.060 while
deleting its own proliferation genes *raises* it by 0.012. Same for
`dang_myc_core_genes` (-0.084 vs +0.033) and `MYC_UP.V1_UP`, whose 3
proliferation genes make the strip a no-op while the adjustment costs 0.037.

In SCAN-B the adjustment is uniformly more destructive: it lowers rho for 17 of
18 signatures, by as much as 0.24 (`YU_MYC_TARGETS_UP` 0.426 -> 0.183) and 0.21
(`PID_MYC_ACTIV_PATHWAY` 0.427 -> 0.209), where stripping lowers it for 13.

**Why they differ.** Adjusting removes shared *variance*; stripping removes
shared *genes*. MYC activity genuinely drives proliferation, so projecting
proliferation out of the plane projects out part of MYC with it. Deleting the
proliferation genes from the signature does not - the remaining genes still
carry whatever MYC-driven proliferation there is in the samples.

**Consequence.** F1's original "proliferation-adjusted" column was
**understating** the correlation, not correcting it, and doing so unevenly
across a panel the covariate was disjoint from in only one member. The
`__PROLIFSTRIP` and `__BOTHSTRIP` variants are the right instrument for that
question and the adjustment should be reported beside them, never instead.

---

## E2. Internal coherence measures proliferation content, not signature quality

The obvious reading of the handoff's question - are `ELLWOOD` and `ALFANO` weak?
- was to check whether their genes co-vary. They barely do. But neither do the
genes of several signatures that work perfectly well:

```
                          coherence   frac_prolif   rho_TCGA / rho_SCANB
PID_MYC_PATHWAY             -0.002        8.3%        0.451 / 0.449
ELLWOOD_MYC_TARGETS_UP       0.004        7.7%        0.181 / 0.049
ALFANO_MYC_TARGETS           0.007        5.0%        0.320 / 0.160
BILD_MYC_ONCOGENIC_SIG       0.021        3.9%        0.398 / 0.352
MYC_UP.V1_UP                 0.041        1.5%        0.552 / 0.502
...
HALLMARK_MYC_TARGETS_V2      0.256       19.3%        0.580 / 0.584
YU_MYC_TARGETS_UP            0.409       47.6%        0.317 / 0.426
```

**Coherence tracks proliferation entanglement at Spearman 0.641 (TCGA) / 0.618
(SCAN-B)** and tracks the OXPHOS correlation at only 0.253.

That is not a nuisance, it is the explanation. Proliferation genes are tightly
co-expressed, so a signature loaded with them looks internally coherent. A
signature of genuine, functionally diverse MYC targets does not. **Low
coherence in a MYC signature means low proliferation content, and is closer to
a virtue than a defect.** `PID_MYC_PATHWAY` has the lowest coherence in the
panel and a rho of 0.45 that replicates to three decimal places.

**Agreement with the rest of the panel is the discriminator that works**
(rho 0.637 with the OXPHOS correlation, against coherence's 0.253):

- **`ELLWOOD_MYC_TARGETS_UP` fails on every count.** Agreement 0.202, the lowest
  in the panel by a wide margin - the next is 0.485 - on 13 genes that do not
  co-vary, with rho 0.181 / 0.049 and the third-worst cross-cohort gap.
  **It should be dropped from the panel**, and F1's entanglement slope
  recomputed without it.
- **`ALFANO_MYC_TARGETS` is unreliable rather than weak.** Agreement 0.492 is
  mid-panel, but it has the worst cross-cohort gap of any signature (0.320 vs
  0.160), it loses the most of any signature to BOTHSTRIP in SCAN-B (0.219, to
  -0.059), and 238 genes with coherence 0.007 is a set with no internal
  structure at all. Report it, do not lean on it.
- **`MYC_UP.V1_UP` and `BILD` are vindicated.** Low coherence, high agreement
  (0.749, 0.666), and rho that replicates. Their low coherence is their low
  proliferation content, which is exactly why F1 leaned on `MYC_UP.V1_UP`.

---

## E3. M_b is the ACTIVATED half of its own regulon, and that is the weak half

The handoff asked why `M_b` behaves differently. The strip is part of it -
`M_b__FULL` 0.320 / 0.223 against `M_b__MITOSTRIP` 0.240 / 0.117, so the
mitochondrial strip costs the regulon two to three times what it costs any GSVA
signature. But the bigger reason is structural.

Splitting the regulon by its own sign and scoring each half **unsigned**, as a
mean of per-gene z:

```
                          TCGA    SCAN-B
activated half (736 genes)
  vs M_b                  0.955   0.934      <- M_b IS this half
  vs OXPHOS subunits      0.246   0.100
  vs mitochondrial ribo   0.308   0.114
repressed half (72 genes)
  vs M_b                  0.068   0.086      <- M_b barely sees this half
  vs OXPHOS subunits     -0.518  -0.578
  vs mitochondrial ribo  -0.606  -0.682
activated MINUS repressed
  vs OXPHOS subunits      0.595   0.604
  vs mitochondrial ribo   0.716   0.711
```

**The MYC-OXPHOS signal in this regulon is carried far more by MYC's repressed
targets going down (-0.52 / -0.58) than by its activated targets going up
(+0.25 / +0.10).** And `M_b` is, empirically, the activated half: it correlates
0.955 / 0.934 with it and 0.068 / 0.086 with the repressed half.

ULM is signed, so the repressed targets *should* add to the score. They are
outnumbered **10.2 to 1** (736 against 72), so the fitted slope is dominated by
the activated set - and 86 of those 736 are edges flagged both stimulatory and
inhibitory that the upstream sign rule silently assigns +1.

**A balanced activated-minus-repressed contrast beats `M_b` by a factor of 2.5
to 5** on this plane (0.595 / 0.604 against 0.240 / 0.117), and beats every GSVA
signature except `MENSSEN` on the mitoribosome (0.716 / 0.711).

**Two things must travel with that.** The contrast is heavily
proliferation-loaded - it correlates 0.660 / 0.669 with `PROLIF_DISJOINT`, more
than it does with OXPHOS - so it is not a clean MYC estimator, it is a strong
one. And the 72-gene repressed set is small, so its own score is noisier than
its correlations suggest.

**What would falsify E3:** a regulon with a balanced sign split, or CollecTRI's
`COLLECTRI_MYC_STIM` (739 genes, activation-only, already in the snapshot and
never scored here) behaving like the activated half as predicted. That is the
cheapest next test and it needs one more variant in `E02`.

---

## E4. Four FULL signatures contain a BCL2-family gene, and E05's overlay must say so

`ALFANO`(MCL1), `BILD`(BID), `DANG_MYC_TARGETS_UP`(BAX) and
`PID_MYC_ACTIV_PATHWAY`(BAX, PMAIP1). A correlation between one of those
signatures and the gene it contains is not an independent observation.

`M_b__MITOSTRIP` contains none of the 15 - all 15 are MitoCarta genes and the
mitochondrial strip removed every one that was in the regulon - so **D3's
BCL2-family numbers are independent**. `M_b__FULL` retains them and must never
be used for a BCL2-family claim. `E06` stops if that ever changes.

---

## What phase 2 should do next

1. **Drop `ELLWOOD` from the panel** and recompute F1's entanglement slope.
   Flag `ALFANO` as unreliable rather than dropping it.
2. **Score `COLLECTRI_MYC_STIM`** as a fifth base in `E02`. It is the
   activation-only regulon, it is already in the snapshot, and E3 predicts it
   will behave like the activated half - about +0.25 / +0.10 on OXPHOS.
3. **Report the balanced regulon contrast beside `M_b`**, with its
   proliferation loading stated.
4. **Never report a proliferation-adjusted rho without the `__PROLIFSTRIP`
   value beside it.** E1 shows they disagree in sign.

---

# Phase 2 findings - the anatomy of the mtDNA divergence

From `E07` on `results/mtdna_anatomy.rds`. Answers handoff candidate 1:
`MT-CO1` and `MT-CO2` are adjacent and go opposite ways, and the earlier note
could only say their mutual correlation was unremarkable - an absence of
evidence, not a test.

## M1. The bicistronic control settles it

`MT-ATP8`/`MT-ATP6` and `MT-ND4L`/`MT-ND4` are each **one mature mRNA**. Two
genes on one transcript cannot be transcribed apart, so the gap between them in
MYC deviation is the scale of purely non-transcriptional variation - the tightest
available control, measured in the same samples on the same axis.

```
pair            gap in MYC deviation      one mRNA?
                TCGA      SCAN-B
ATP8 vs ATP6    0.047     0.148           yes
ND4L vs ND4     0.158     0.135           yes
CO1 vs CO2      0.735     0.526           no
```

**The CO1-CO2 gap is 4.7x (TCGA) and 3.6x (SCAN-B) the largest bicistronic
gap.** Genes that physically cannot be separated differ by at most 0.16. `CO1`
and `CO2` differ by 0.5 to 0.7.

`CO1` and `CO2` are on the same heavy-strand polycistron but are excised into
separate mature mRNAs at the flanking tRNAs. So their divergence involves
something that happens **after the polycistron is cut** - differential stability,
processing or degradation of separate mature transcripts - and it is far outside
what shared-transcript arithmetic can produce.

## M2. It replicates, and it is not positional

- **Per-gene MYC deviation replicates across cohorts at Spearman 0.888**
  (Pearson 0.917) over the 12 trusted genes.
- Deviation against position along the heavy strand: **-0.105 / -0.140**. No
  3' gradient, no decay ramp.
- Deviation against expression level: **-0.021 / +0.140**. Not abundance.

## M3. The reproducible axis is PC2, and PC1 is a cohort artefact

Conditioning on mtDNA content leaves 20.0% / 19.0% of the rank variance. Within
that residual:

```
        variance explained    loadings replicate    aligns with MYC deviation
        TCGA   SCAN-B         across cohorts        TCGA      SCAN-B
PC1     23.4%  20.3%          -0.168                -0.350    +0.168
PC2     18.0%  16.6%          +0.629                +0.762    +0.615
```

**The larger axis does not replicate and its alignment with the MYC deviation
flips sign between cohorts. The second one does both.** PC1 correlates with
every nuclear arm at about -0.41 in TCGA and -0.09 in SCAN-B - a 4-fold
asymmetry that marks it as cohort-specific, most likely purity or composition.

`PC2` tracks proliferation (0.298 / 0.425), MYC activity (0.270 / 0.254), the
mitoribosome (0.130 / 0.283) and fatty-acid oxidation negatively (-0.120 /
-0.113) - consistent in sign on almost everything, modest in size.

**The size caveat matters.** PC2 is 17-18% of a residual that is itself 20% of
the gene-level variance, so it accounts for roughly **3.5%** of the variation in
the 13 genes. It is a reproducible axis, not a large one.

## M4. What this does and does not license

**Licensed:** there is a reproducible, non-positional, post-transcriptional
divergence within the mitochondrial genome that tracks MYC activity, with
`MT-CO2` at one pole and `MT-CO1`, `MT-ND1` and `MT-ND5` at the other, and it is
larger than anything shared-transcript variation can explain.

**Not licensed:** any statement about mechanism. Nothing here distinguishes
differential mRNA stability from differential processing from differential
degradation, and the axis carrying it is small.

**`MT-ND6` remains excluded.** The raw correlation matrices re-confirm it:
0.46-0.86 with the other genes in TCGA against 0.12-0.36 in SCAN-B, and it is
the highest residual pair with `MT-ND5` (0.75) in TCGA - the two overlapping
genes, coupled in one cohort and detached in the other.

## What would falsify M1-M3

1. A stranded, total-RNA dataset. It resolves `MT-ND6` outright and tests
   whether `CO1`/`CO2` survives without polyA selection. **This is the one test
   that matters most and phase 1 cannot do it.**
2. A third cohort: the deviation vector should reproduce at ~0.9.
3. Direct mtDNA copy number from WGS rather than the expression proxy.
4. Ribosome profiling or protein, to separate "not transcribed" from "not
   translated" - the fork none of this can resolve.

---

# Phase 2 findings - strata, and the death axis at gene resolution

From `E08` on `results/strata_and_death_genes.rds`. Answers handoff candidate 4
and the author's questions 1 and 2.

## S1. BCL2's correlation with MYC is mostly BETWEEN subtypes. D3 must be qualified.

`D3` reported `BCL2` at **-0.369 / -0.288** against MYC activity. Split by
stratum, **every single stratum in both cohorts is weaker than the pooled
value**:

```
BCL2 vs MYC       TCGA     SCAN-B
all              -0.369   -0.288
ERpos            -0.150   -0.132
ERneg            -0.180   -0.175
Luminal          -0.110   -0.083
LumA             -0.009   -0.003     <- essentially zero
LumB             -0.142   -0.146
Basal            -0.068   -0.165
HER2             -0.173   -0.252
Normal           -0.170   -0.106
```

ER-negative and basal tumours are both MYC-high and BCL2-low, so pooling
manufactures most of the correlation. **Within LumA - the largest single
subtype in both cohorts - it is zero.**

**The same inflation, milder, affects the pro-death side.** `BID` 0.42 pooled
against 0.20-0.26 in the luminal strata; `BAK1` 0.41 against 0.15-0.27.

**But the OXPHOS axis is stratum-stable.** `BAX` runs 0.45-0.58 across every
stratum against a pooled 0.48; `BCL2L11` -0.26 to -0.35 against -0.30; `MCL1`
-0.21 to -0.31 against -0.25.

> **D3's BCL2-family claim should be made against OXPHOS, not against MYC.**
> Against OXPHOS the numbers hold within every subtype. Against MYC they are
> substantially a between-subtype artefact of pooling, and `BCL2` in particular
> must not be quoted at -0.37.

This is the same shape as D1, which also turned out to be an OXPHOS finding
rather than a MYC one.

## S2. The D1 contrast is not self-overlap, and that is now measured

Only **three** genes of the 1,086 are in the `OXPHOS subunits` arm they are
correlated against: `NDUFA2` and `NDUFS3` on the pro-survival side, `CYCS` on
the pro-death side. They are also the three largest contributors by |rho|
(0.65-0.73), so they dominate the "top contributors" list and that list should
not be read as biology.

Deleting them barely moves the contrast, and deleting **every** MitoCarta gene
makes it slightly stronger:

```
OXPHOS-axis contrast     all genes   minus the 3 OXPHOS   minus all 79 MitoCarta
TCGA                      -0.0442      -0.0427              -0.0472
SCAN-B                    -0.0572      -0.0557              -0.0620
```

**D1 is not a self-overlap artefact.** D0's confound does not reach it.

## S3. But the contrast is carried by the tails more than the middle

`frac_positive` shifts across the whole distribution (pro-death 0.38-0.42
against pro-survival 0.45-0.50) and the top decile holds only 0.24-0.26 of the
|rho| mass, both of which suggest a broad shift. **The trim says otherwise.**
Deleting the top N genes by |rho| from each side:

```
dropped   SCAN-B contrast   TCGA contrast
0          -0.0572           -0.0442
10         -0.0473           -0.0353
25         -0.0427           -0.0228
50         -0.0278           -0.0194
100        -0.0093           -0.0123
```

Removing 100 of ~500 genes per side - 20% - removes **84% (SCAN-B) and 72%
(TCGA)** of the contrast. A uniform shift would survive that, because the shift
would be present in the remaining 400 genes too.

> **The honest statement: a broad shift in sign (`frac_positive`) sits on top of
> a contrast whose MAGNITUDE lives in the strongly-correlated tails of both
> strata.** It is neither "a few genes" nor "every gene equally". D1's direction
> is robust; its size is a tail property and should not be quoted as a
> per-gene effect.

## S4. The ER-negative gap is real, is not range restriction, and is arm-specific

**Not range restriction.** `Basal` has the *smallest* MYC spread of any stratum
(IQR 0.317 TCGA / 0.294 SCAN-B) and among the *highest* median rho (0.483 /
0.384). `ERneg`'s spread is mid-range (0.476 / 0.523) and its median rho is the
lowest (0.406 / 0.297). Restricted range would predict the opposite ordering.

**And ERneg is ~70% Basal** (144 of 217 in TCGA, 280 of 472 in SCAN-B), so the
two rows are not independent readings - yet they point opposite ways, which is
itself the finding: the gap is an ER-negative property, not a basal one.

**The gap is not uniform across the mitochondrion.** ERpos minus ERneg, by arm:

```
largest gap                SCAN-B   TCGA        smallest gap             SCAN-B   TCGA
ROS and glutathione         0.315   0.189       Mitochondrial ribosome    0.045   0.013
Fatty acid oxidation        0.301   0.197       OXPHOS assembly factors   0.081   0.011
mtDNA-encoded OXPHOS        0.268   0.114       Folate and 1-C            0.003  -0.135
TCA cycle                   0.257   0.145       CII subunits              0.081  -0.069
```

**The mitochondrial ribosome has essentially no gap** (0.573 vs 0.618 in SCAN-B;
0.582 vs 0.594 in TCGA) - MYC's coupling to it is identical in ER-positive and
ER-negative disease. The metabolic arms carry the whole gap, and in ER-negative
tumours fatty-acid oxidation goes strongly negative (-0.286 / -0.159) where in
ER-positive it is flat.

That maps onto F1's resolution: the mitoribosome is the arm MYC couples to most
tightly and most universally.

## S5. The Luminal stratum is larger than its parts, and that is arithmetic

`Luminal` (LumA + LumB) gives 0.458 in SCAN-B against LumA 0.375 and LumB 0.398
- **higher than either component**. In TCGA 0.491 against 0.493 and 0.384.

Combining two subtypes that differ in mean MYC and mean OXPHOS adds
between-group variance to both axes, which inflates the correlation. The
Luminal number is therefore not a within-compartment statement and must be read
beside LumA and LumB, never instead of them. Same caution as S1.

## What phase 2 should do next

1. **Amend D3 in the cell-death note**: the BCL2-family claim is an OXPHOS
   claim. `BCL2` at -0.37 against MYC is pooling.
2. **Report D1's contrast with its trim curve**, not as a single number.
3. The **ER-negative fatty-acid-oxidation reversal** (-0.29 / -0.16 against ~0
   in ER-positive) is the largest arm-level stratum effect in the study and
   nothing has been done with it.

