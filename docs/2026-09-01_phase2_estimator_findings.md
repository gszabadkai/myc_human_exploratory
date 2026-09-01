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
