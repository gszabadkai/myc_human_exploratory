---
date: 2026-10-22
status: >
  DECLARATION, pre-fit. No model in this document has been fitted, and no
  number below comes from a fit made for it. Section 2 lists everything that
  WAS seen before it was written.
repo: myc_human_exploratory. myc_human_validation and myc_mouse are read-only.
relates-to:
  - E34 (the quadrant means and contrasts this resolves)
  - E46 (Fig. 4g, which draws them on both instruments)
  - E22 (the same main-effects model, gene by gene, on two other OXPHOS rulers)
  - E11 (the conditioning ladder, which concerns gene ORDERING, not level)
decides:
  - The model, exposures, readouts, signature variants, MYC estimator ladder,
    gates, displays and reading rule, all fixed before any fit.
approval: >
  2026-10-22, by the author, in four steps:
  1. APPROVED sections 0-13 as first drafted, including section 7's reading
     rule.
  2. INCLUDED section 9, BCL2L1 against BBC3.
  3. ASKED for sections 4a (gene choices, the sign, Methods text) and 8a (the
     individual genes shown).
  4. DECIDED to add both rule-based variants (4a.6), and to rename "the
     configuration composite" as "the signature", with "signature sign" for
     its signs.
next-action: >
  Commit ALONE. Then the script, unrun; redirected dry runs; digests committed
  BEFORE the author's run; the author's run; the result note.
---

# E47 — MYC and OXPHOS on the BCL-2-family signature, each with the other held fixed

## Terms used from here on

| term | meaning |
|---|---|
| **the signature (S6)** | the six-gene BCL-2-family signature. **Called "the configuration composite" in E33, E34 and every earlier note**, and `config_comp` in code. The renaming is the author's, 2026-10-22. Code identifiers are unchanged |
| **signature sign** | the +1 or -1 with which a gene enters a signature (4a.4) |
| **S5, S7** | the two rule-based variants (4a.6) |
| **trigger arm, guardian arm** | the signature's two parts (section 4) |
| **MYC estimator** | a score of MYC activity, such as `M_a`. **Never called a "signature" here**, so that the word means one thing |

**In the manuscript, write "BCL-2-family signature" at first use.** "Signature"
alone collides with the MYC activity signatures.

## 0. Why this exists

**Fig. 4g (E46) shows a MYC increment at both OXPHOS levels**, most visibly on
mitoPPS. The narrative of record says the opposite: in the abstract, OXPHOS
tracks the signature's genes "even after adjustment for ... MYC expression"; in
paragraph 1, "the association with MYC activity was lost after adjustment for
OXPHOS". Before 4g is kept, redrawn or dropped, the author needs the answer to
two questions:

1. **Does OXPHOS track the signature with MYC in the model?**
2. **Does MYC track it with OXPHOS in the model?**

**Three existing analyses come close. None of them answers the question.**

- **E34's quadrant contrasts do not hold the other axis fixed.** A median split
  leaves the other axis free to vary within each half. The axes correlate
  (`rho(M_a, OXPHOS subunits)` is 0.388 / 0.427 on GSVA and 0.306 / 0.305 on
  mitoPPS, TCGA / SCAN-B, atlas). So within the OXPHOS-high half, MYC-high
  tumours are also more OXPHOS-high. Section 2(d) estimates how much MYC
  contrast that leftover OXPHOS alone would produce.
- **E11's conditioning asymmetry is about how the 44 machinery genes are
  ORDERED across MitoCarta membership** (0.187 to -0.043 for MYC), not about any
  tumour's level. E11 itself records that MYC's own per-gene component "is real
  ... larger than the borrowed part" (reconciliation note 4.1). **The Doc's
  "lost after adjustment" sentence is a translation of the ordering result into
  a level claim, which E11 does not make.**
- **E22 fitted `gene ~ MYC + OXPHOS` for all twelve genes**, but not the
  signature or the arms, and with three differences from 4g:
  - its OXPHOS rulers are `ox_rel` and `ox_lvl`, not 4g's GSVA and mitoPPS;
  - it has no PAM50 term;
  - its analysis set is not E34's.

## 1. Posture

**EXPLORATORY and POST-HOC.** It was written after E34's and E46's group means
had been seen, and in response to them. Nothing here is pre-registered. It
fixes its rule before its own fits, which is all a post-hoc analysis can do.

**MAIN EFFECTS ONLY. No product term is fitted** (section 3).

## 2. What was seen before this was written

**(a) E34's group means and contrasts, both instruments, `m1` and `m3`.**

The MYC contrasts on the signature at `m1`, as `Q4 - Q2` / `Q3 - Q1`:

| | GSVA | mitoPPS |
|---|---|---|
| TCGA | -0.002 / +0.134 | +0.185 / +0.163 |
| SCAN-B | +0.134 / +0.128 | +0.232 / +0.222 |

At `m3` the MYC contrasts are all larger, by about +0.12. The OXPHOS contrasts
are +0.40 to +0.98.

**(b) E22's gene-level coefficients**: twelve genes, five estimators, two
rulers, three rungs.

At `+PROLIF` on `M_a`, MYC at fixed OXPHOS, `beta_MYC|OX`, for the signature's
six genes:

| gene | signature sign | `beta_MYC\|OX` |
|---|---|---|
| `BID` | + | +0.29 to +0.39; positive on all four activity estimators, both cohorts, both rulers |
| `BBC3` | + | +0.05 to +0.13; sign depends on the estimator |
| `BAD` | + | -0.02 to +0.11; sign depends on the estimator |
| `BIK` | + | -0.19 to +0.02 |
| `BCL2L1` | + | -0.05 to -0.28 on `M_a` |
| `MCL1` | − | +0.01 to +0.28 |

**So within the signature, MYC's gene-level coefficients run in opposite
directions.** OXPHOS at fixed MYC, `beta_OX|MYC`, is positive for all six once
the signature sign is applied. For `BMF`, which is not in the signature,
`beta_MYC|OX` is -0.05 to -0.23.

**(c) The atlas's `rho(M_a, OXPHOS subunits)`.**

| | GSVA | mitoPPS |
|---|---|---|
| raw, TCGA / SCAN-B | 0.388 / 0.427 | 0.306 / 0.305 |
| proliferation-adjusted, TCGA / SCAN-B | 0.320 / 0.256 | 0.199 / **0.070** |

**(d) Arithmetic on the numbers above, not a fit.** Assume standard bivariate
normal scores with correlation `rho`. Within the OXPHOS-high half, MYC-high
tumours then average a higher OXPHOS score than MYC-low tumours by:

| `rho` | 0.20 | 0.30 | 0.388 | 0.427 |
|---|---|---|---|---|
| gap (SD) | 0.12 | 0.18 | 0.23 | 0.25 |

This was checked by simulation, 0.2305 at 0.388.

Multiply the gap by the per-SD OXPHOS slope implied by E34's OXPHOS contrasts.
That gives the MYC contrast leftover OXPHOS alone would produce, at about:
- **0.13 for TCGA/GSVA**, against observed -0.002 / +0.134;
- 0.10 for SCAN-B/GSVA, against +0.13;
- 0.08 for TCGA/mitoPPS, against +0.17 to +0.19;
- **0.05 for SCAN-B/mitoPPS, against +0.22 to +0.23.**

**This is an approximation**: it ignores the `m1` covariates and assumes
linearity. It motivated this declaration, and it is NOT a result.

**(e) The reconciliation note's section 3.1 table** (reproduced in 4a.2), from
which the signature's genes and signs were read. The variants of 4a.6 were
defined after it, and after (b), had been seen.

Nothing else was seen. **In particular, no continuous model of the signature,
its variants or either arm has ever been fitted.**

## 3. The model

For each cohort (TCGA, SCAN-B), OXPHOS instrument (GSVA, mitoPPS), MYC
estimator (section 5) and rung:

```
readout ~ OX + MYC + spine          OLS; t-based 95% CI, as E34
```

- **readout**: as section 4.
  - The signature (S6), its arms and the twelve genes are E34's values, from
    `e34_quadrant_configuration.rds$analysis_frames`, **unchanged**. These are
    Blom normal scores at unit variance, per cohort.
  - **S6 and the arms are also rebuilt here from expression, and must
    reproduce E34's values exactly** (gate G4).
  - **S5 and S7 are built here by the same function** (4a.6).
- **OX**: the `OXPHOS subunits` score on that instrument, **the same continuous
  score E34 split** (gate G1).
  - GSVA on VST; mitoPPS on linear counts.
  - **Read, never rescored.**
- **MYC**: the estimator's score as E02 computed it (section 5), read, never
  rescored.
- **Both exposures are Blom-scored per cohort**, exactly as E34 scored its
  readouts and `PROLIF`: once per cohort, after which the fits subset and
  never rescore. **So every readout and every exposure is on one convention.**
- **Rungs:**
  - `m0`: no covariates;
  - `m1`: E34's spine, which is PAM50 + purity + leukocyte fraction in TCGA,
    and PAM50 only in SCAN-B (no purity estimate; trap 2);
  - `m3`: `m1` + `PROLIF_DISJOINT`.

  `m0` is reported and never read. **`m1` and `m3` are read, side by side, and
  never chosen between.**
- **Analysis sets: E34's `m1` complete cases, TCGA 938 and SCAN-B 3,143.** They
  are identical across instruments, estimators, rungs and readouts, so every
  comparison is on the same tumours.
- **Snapshot: `data/from_validation/` and E02's scores, as E34.** Genes are on
  log2(linear DESeq2-normalised + 1). SCAN-B symbols resolve through
  `scanb_pheno.rds$symbol_map` (trap 7). The joint snapshot is not read.

**No product term.** The question is about main effects. Interaction has been
asked twice already:
- the pre-registered sibling study found the MYC x OXPHOS interaction on
  priming null, and E22 records it as "permanently closed";
- this arm tested it on the priming ratios and did not support it
  (reconciliation note N1: "the apparent interaction is curvature in the main
  effects").

**Adding a product term would need its own declaration.**

## 4. Readouts

| role | readout |
|---|---|
| **primary, the verdict** | **the signature, S6**: `BBC3`, `BID`, `BIK`, `BAD`, `BCL2L1` and `MCL1`, with their signature signs (4a.4) |
| reported always, with the rule applied descriptively | **the two variants, S5 and S7** (4a.6) |
| reported always, separately | **the trigger arm** (`BBC3`, `BID`, `BIK`, `BAD`) and **the guardian arm** (`BCL2L1` against `MCL1`) |
| reported, never the verdict | the twelve genes individually. **No per-gene FDR** |
| secondary, included by the author | `BCL2L1` against `BBC3` (section 9) |

**Which genes, why, how the signature is built, and what its signs mean:
section 4a.**

**The arms are what separate "a pro-apoptotic transcriptome" from "its
compensation by BCL-xL".**
- The trigger arm is the mean of the four BH3 genes' z-scores.
- The guardian arm is `(z_BCL2L1 - z_MCL1) / 2`, that is, BCL-xL relative to
  MCL1.
- The signature contains both, so it cannot separate them.

## 4a. The genes: which, why, and the signature sign

*Added at the author's request, 2026-10-22. It records provenance and defines
terms. Its last parts are the two variants and draft Methods text.*

### 4a.1 The twelve

They are **the author's lists**: E10, 2026-09-02,
`scripts/E10_machinery_measures_and_priming.R:117-118`. They are the core
BCL-2-family regulators of apoptotic priming.

| group | genes (protein) |
|---|---|
| **BH3-only, pro-apoptotic (7)** | `BCL2L11` (BIM), `BID`, `BBC3` (PUMA), `BAD`, `BIK`, `BMF`, `PMAIP1` (NOXA) |
| **anti-apoptotic (5)** | `BCL2`, `BCL2L1` (BCL-xL), `MCL1`, `BCL2L2` (BCL-w), `BCL2A1` (BFL-1) |

- **`BCL2L2` was first listed on both sides.** The author confirmed the
  pro-apoptotic entry was an error (E10, R5b).
- **The effectors (`BAX`, `BAK`, `BOK`) and `HRK` are not in the lists.**
- **All twelve clear the 25th expression percentile in both cohorts** (E10), so
  no pattern below is a low-expression artefact.
- **Earlier notes call the signature's four BH3 genes "sensitisers".** That was
  loose. `BID` is usually classed as an activator, and `BBC3` (PUMA) as an
  activator or a sensitiser depending on the study. **Write "BH3-only
  genes".**

### 4a.2 The six genes of the signature, and how they were chosen

**They come from one table**: reconciliation note
(`2026-09-03_human_arm_for_mouse_reconciliation.md`) section 3.1. It gives
each gene's Spearman correlation with the OXPHOS subunit score (GSVA), **partial
on `PROLIF_DISJOINT`**, in each cohort, ordered by the mean of the two:

| gene | function | OXPHOS: TCGA / SCAN-B | **mean** | MYC: TCGA / SCAN-B | mean | S6 | S5 | S7 |
|---|---|---|---|---|---|---|---|---|
| `BAD` | pro | +0.525 / +0.482 | **+0.503** | +0.162 / +0.114 | +0.138 | +1 | +1 | +1 |
| `BCL2L1` | anti | +0.426 / +0.350 | **+0.388** | -0.030 / +0.033 | +0.002 | +1 | +1 | +1 |
| `BIK` | pro | +0.338 / +0.213 | **+0.276** | +0.008 / +0.051 | +0.029 | +1 | +1 | +1 |
| `BBC3` | pro | +0.365 / +0.170 | **+0.267** | +0.163 / +0.110 | +0.136 | +1 | +1 | +1 |
| `BID` | pro | +0.329 / +0.091 | **+0.210** | **+0.357 / +0.243** | **+0.300** | +1 | **–** | +1 |
| `BCL2A1` | anti | -0.004 / -0.107 | -0.056 | +0.095 / +0.014 | +0.055 | – | – | – |
| `BCL2L2` | anti | -0.145 / -0.008 | -0.077 | -0.073 / -0.080 | -0.076 | – | – | – |
| `PMAIP1` | pro | -0.147 / -0.064 | -0.105 | -0.194 / -0.223 | -0.209 | – | – | – |
| `BCL2` | anti | -0.160 / -0.074 | -0.117 | -0.222 / -0.160 | -0.191 | – | – | – |
| `BCL2L11` | pro | -0.314 / -0.025 | -0.170 | -0.162 / -0.039 | -0.101 | – | – | – |
| `BMF` | pro | -0.200 / -0.223 | **-0.212** | -0.150 / -0.079 | -0.115 | – | – | **−1** |
| `MCL1` | anti | -0.241 / -0.292 | **-0.266** | +0.081 / -0.037 | +0.022 | −1 | −1 | −1 |

**The history:**
1. **The reconciliation note named the result** "a specific, replicating
   transcript configuration - sensitisers and `BCL2L1` up, `MCL1` down, in
   OXPHOS-high tumours" (its section 1).
2. **E33's declaration (2026-09-29) fixed it as two arms**: trigger, `BBC3`,
   `BID`, `BIK`, `BAD`; guardian, `BCL2L1`, `MCL1`.
3. **E33's data note (author's decision, 2026-09-30) fixed the signs and the
   scoring**, before any arm was correlated with anything.
4. **E34 adopted the signed score as its primary readout** (author's decision,
   2026-09-30, before any data was read).

**NO NUMERIC SELECTION RULE WAS EVER WRITTEN DOWN.** Described against the
table, the six are:
- **the five genes with the highest mean correlation with OXPHOS**, all at
  +0.21 or above, which are higher in OXPHOS-high tumours;
- **`MCL1`, the most negatively correlated gene** (-0.27). It is the gene
  whose movement opposite to `BCL2L1` makes the guardian "switch" the
  signature was named for.

**Two facts a reader or reviewer will ask about, recorded here:**

- **`BID` is the only member more associated with MYC than with OXPHOS**
  (+0.300 against +0.210 on the means). It is **not** one of the five genes
  this table calls "OXPHOS-led", meaning `|rho(OXPHOS) - rho(MYC)| >= 0.20` in
  at least one cohort: `BCL2L1`, `BAD`, `BIK`, `BBC3` and `MCL1`. It entered
  as the fourth BH3 gene that rises with OXPHOS. **E22 shows `BID` carrying
  most of MYC's gene-level association within the signature (section 2b), so
  its inclusion bears directly on this analysis's primary question.**
- **`BMF` is excluded, although its OXPHOS association is consistent in both
  cohorts (-0.200 / -0.223) and nearly as strong as `MCL1`'s.** **No written
  reason exists.** By the logic of the signature sign it would enter with -1.

### 4a.3 The six on Fig. 4g that are not in the signature

| gene | why it is not in the signature, as the table shows |
|---|---|
| `BCL2A1` | its OXPHOS association is near zero and inconsistent between cohorts (-0.004 / -0.107) |
| `BCL2L2` | weak and inconsistent (-0.145 / -0.008) |
| `PMAIP1` | weak, and it tracks MYC more than OXPHOS (-0.209 against -0.105) |
| `BCL2` | weak, and it tracks MYC more than OXPHOS (-0.191 against -0.117) |
| `BCL2L11` | strong in TCGA only (-0.314 against -0.025) |
| `BMF` | **consistent and moderately strong; excluded without a written reason** (4a.2). **It enters S7** (4a.6) |

**All six are shown individually on Fig. 4g and fitted individually here.**

### 4a.4 The signature sign

**Definition.** A gene's **signature sign** `s_g` is the direction in which its
expression moves in OXPHOS-high tumours: the sign of its mean correlation with
OXPHOS in the 4a.2 table.
- `+1` if it is higher there: `BBC3`, `BID`, `BIK`, `BAD`, `BCL2L1`.
- `-1` if it is lower: `MCL1`, and `BMF` where it enters S7.

**A signature score, for each tumour, built in each cohort separately:**

```
score = Blom( (1/k) * sum over the signature's k genes of  s_g * z_g )
```

`z_g` is the gene's log2 normalised expression, standardised across the whole
cohort. Blom is a rank normalisation within the cohort. S6 has `k = 6`; the
variants have 5 and 7.

**Four properties:**

1. **It is a direction, not a function.** `+1` does not mean pro-apoptotic and
   `-1` does not mean anti-apoptotic. The signature's two anti-apoptotic genes
   take opposite signs: `BCL2L1` is `+1` and `MCL1` is `-1`. **A high score
   means "this tumour shows the OXPHOS-associated pattern": four BH3 genes up,
   and BCL-xL up relative to MCL1. It does NOT mean "more pro-apoptotic".**
2. **It is fixed, not estimated.** S6's signs were set once (E33's data note,
   2026-09-30). The variants' signs follow from the same table. **None is
   re-derived from any analysis's data.**
3. **It is why the OXPHOS association is expected by construction** (E34 1.2).
   The selection used the GSVA OXPHOS score, so this applies most directly on
   GSVA. mitoPPS was not selected on directly, but the two OXPHOS scores
   correlate.
4. **On Fig. 4g, each gene row shows that gene's own level, unsigned.** The
   signature-sign column says how each gene enters the signature. So `MCL1`'s
   row is red where `MCL1` is high, which lowers the score.

### 4a.5 Methods text, draft for the manuscript

> **BCL-2-family genes and signature.** We examined twelve BCL-2-family genes
> that set the threshold for mitochondrial apoptosis: seven BH3-only proteins
> (BCL2L11/BIM, BID, BBC3/PUMA, BAD, BIK, BMF and PMAIP1/NOXA) and five
> anti-apoptotic proteins (BCL2, BCL2L1/BCL-xL, MCL1, BCL2L2/BCL-w and
> BCL2A1/BFL-1). For each gene we computed the Spearman correlation of its
> expression with the OXPHOS subunit score (GSVA), partial on a proliferation
> score, in TCGA-BRCA and SCAN-B. Six genes defined an OXPHOS-associated
> BCL-2-family signature: the five most positively correlated (BAD, BCL2L1,
> BIK, BBC3 and BID; mean rho 0.21 to 0.50), and MCL1, the most negatively
> correlated (mean rho -0.27), whose movement opposite to BCL2L1 marks a switch
> between the two anti-apoptotic proteins. Each gene was assigned a signature
> sign, +1 for the five and -1 for MCL1, recording the direction in which it
> moves with OXPHOS; the sign records direction, not pro- or anti-apoptotic
> function. A tumour's signature score was the mean of the six signed z-scores
> (log2-normalised expression, standardised across the cohort),
> rank-normalised within each cohort. Because these genes and signs were chosen
> by their association with OXPHOS in the same cohorts, the signature's
> association with OXPHOS is expected by construction; the informative quantity
> is its association with MYC activity at fixed OXPHOS. To test whether
> conclusions depended on the choice of genes, all analyses were repeated with
> two rule-based variants: the five genes whose correlation with OXPHOS
> exceeded their correlation with MYC by at least 0.20 in at least one cohort
> (BID omitted), and the seven genes with a mean absolute correlation with
> OXPHOS of at least 0.20 (BMF added, signature sign -1). All twelve genes are
> also reported individually.

**This paragraph must NOT be shortened to drop the construction sentence.** It
is what makes the OXPHOS association honest.

### 4a.6 The two rule-based variants

*Added at the author's decision, 2026-10-22.* **No numeric rule chose S6**
(4a.2). Each variant applies one explicit rule to the 4a.2 table, so the
question "why these genes?" has a stated answer, and the analysis shows
whether the answer to section 0 depends on it.

| name | rule, applied to the 4a.2 table | genes, with signature signs | differs from S6 by |
|---|---|---|---|
| **S6, the signature (primary)** | as chosen in 2026-09 (4a.2) | `BAD`+, `BIK`+, `BBC3`+, `BID`+, `BCL2L1`+, `MCL1`− | — |
| **S5, OXPHOS-led** | `\|rho(OXPHOS) - rho(MYC)\| >= 0.20` in at least one cohort: **the reconciliation note's own criterion** | `BAD`+, `BIK`+, `BBC3`+, `BCL2L1`+, `MCL1`− | **`BID` removed** |
| **S7, threshold** | `\|mean rho(OXPHOS)\| >= 0.20` | `BAD`+, `BIK`+, `BBC3`+, `BID`+, `BCL2L1`+, `BMF`−, `MCL1`− | **`BMF` added** |

- **Both rules select exactly the genes listed.** Checked against all twelve
  rows of the table.
- **`BBC3` enters S5 at the margin**: 0.202 in TCGA, 0.060 in SCAN-B.
- **Both are built exactly as S6** (4a.4), by the function that must first
  reproduce S6 exactly (G4).
- **BOTH ARE POST HOC.** They were defined after the table and E22's gene-level
  coefficients had been seen. They test the selection; they are not new
  signatures to prefer. **S6 remains the primary and its verdict is the
  verdict.**

## 5. The MYC estimator ladder, which N9 makes mandatory

N9 (reconciliation note): **on BCL2-family transcripts, the sign of MYC's main
effect is a function of the estimator's proliferation entanglement.** "No
single-estimator MYC main effect ... is interpretable without its entanglement
spread." So E22's four activity estimators are fitted:

| role | estimator | proliferation entanglement |
|---|---|---|
| **primary** | `FELSHER__MITOSTRIP` (`M_a`), the axis of 4g | 14.8% |
| ladder | `MYC_UP.V1_UP__MITOSTRIP` | 1.5%, the cleanest |
| ladder | `HALLMARK_MYC_TARGETS_V1__MITOSTRIP` | 23.5% |
| ladder | `M_b__MITOSTRIP`, the CollecTRI regulon | E22's "regulon" |

**None contains any of the twelve genes** (E22 `$contamination`; asserted
again, G2). E22's own estimator labels ("dose", "signature", "regulon") are not
used here, because "signature" now names the BCL-2-family readout.

**Not fitted:**
- **`M_b__FULL` and `M_b__PROLIFSTRIP`.** Each contains 5 of the 12.
- **log2 `MYC` mRNA.** It is not activity (trap 4).

## 6. Gates, before any readout is modelled

- **G1, identity of the exposures.** E34's quadrant labels must be a threshold
  split of this script's continuous scores, on each cohort and instrument:
  - every MYC-high tumour's `M_a` exceeds every MYC-low tumour's;
  - likewise for OX on that instrument;
  - the analysis sets are E34's, 938 and 3,143.

  **If G1 fails, nothing is fitted.**
- **G2, self-overlap.** None of the twelve may sit in:
  - any estimator in section 5;
  - `OXPHOS subunits`;
  - `PROLIF_DISJOINT`.

  This is asserted.

  **One known indirect overlap is recorded rather than rebuilt.** mitoPPS
  scores OXPHOS relative to the rest of the mitochondrial programme, where 11
  of the 12 sit, at roughly 1/1,000 weight each. E22's leave-one-out rebuild of
  the analogous `ox_rel` moved every coefficient by at most **0.004**.
- **G3, identifiability, with E39's bands.** Spearman correlations among OX,
  MYC and PROLIF are computed per cohort, instrument and estimator.
  - **Bands:** below 0.60, fine; 0.60-0.80, not cleanly separable; 0.80 or
    above, **not identifiable**.
  - **If `rho(MYC, PROLIF)` is 0.80 or above, that estimator's `m3` cells in
    that cohort are VOID.** They are reported and never read.
  - The full VIF vector is saved for every fit (R rule 5). A multi-df term is
    read as `GVIF^(1/(2*Df))` against `sqrt(5)`.
- **G4, reconstruction of the readouts.** S6 and both arms, rebuilt here from
  expression by E34's recipe, must equal E34's saved values. The maximum
  absolute difference must be below 1e-10 in both cohorts. **S5 and S7 are
  built only by the function that passed G4. If G4 fails, nothing is
  fitted.**

## 7. Reading rule, fixed now: `beta_MYC|OX` on the signature

**Primary cells: `M_a` x 2 cohorts x 2 instruments x {`m1`, `m3`} = 8.**
Cells made void by G3 are excluded, and every count below is over **non-void**
cells. Each ladder estimator likewise has 8 cells. **The rule is applied in
this order, and the first verdict that matches is the verdict.**

| verdict | condition |
|---|---|
| **MYC ADDS** | On `M_a`, `beta_MYC > 0` with its CI excluding zero in **every** primary cell, **and** each of the other three estimators gives `beta_MYC > 0` in at least three-quarters of its cells (sign only) |
| **ESTIMATOR-DEPENDENT** | On `M_a`, every primary cell has its CI above zero, but at least one other estimator gives `beta_MYC < 0` in more than a quarter of its cells. **This is the N9 pattern** |
| **INSTRUMENT-SPLIT** | On `M_a`, every cell of one instrument (2 cohorts x 2 rungs) has its CI above zero, and every cell of the other has its CI including zero or lying below it |
| **OXPHOS ONLY** | On `M_a`, `beta_MYC`'s CI includes zero, or `beta_MYC < 0`, in at least three-quarters of the primary cells |
| **MIXED** | anything else |

MYC ADDS and ESTIMATOR-DEPENDENT are complementary once `M_a` passes, so the
order decides nothing between them.

**THE VERDICT IS S6's ALONE.** The same rule is applied to S5 and S7
**descriptively**, and their labels are printed beside S6's.
- **If either variant's label differs from S6's**, the result note says so in
  its first section.
- **No manuscript sentence may then state S6's verdict without naming that
  dependence**, for example: "MYC's contribution is carried by `BID`."

**`beta_OX|MYC` is reported in every cell and NEVER read.** It is
near-guaranteed: the signature's genes and signs were selected on OXPHOS in
these same cohorts (E34 1.2).

**`beta_MYC / beta_OX` is printed as a description only.** Because `beta_OX`
is inflated by construction, the ratio understates MYC's share. It is never
quoted without that sentence.

**The arms and the twelve genes get the same table, with the rule NOT
applied.**

## 8. Diagnostics, reported, no verdict

- **D1, leakage.** Within E34's quadrants on each instrument, report:
  - the mean OX difference between MYC-high and MYC-low tumours, at each
    OXPHOS level;
  - that difference times `beta_OX|MYC` at `m1`, which is the MYC contrast
    leftover OXPHOS would produce;
  - E34's observed MYC contrasts, beside it.

  **This reconciles 4g with this analysis**, and replaces section 2(d)'s
  approximation with a measurement.
- **D2.** The movement of `beta_MYC` and `beta_OX` from `m0` to `m1` to `m3`.

## 8a. Displays, declared before any fit

*Added at the author's request, 2026-10-22: "to convince the reader and
reviewer we need to show the individual genes, too".* The content is fixed
here so that it cannot be chosen after the numbers are seen. **Whether a
display goes in the main figure or Extended Data is the author's call after
the verdict.**

**Display A: the genes, both coefficients.**
- **Rows:**
  - at the top, S6, S5, S7, the trigger arm and the guardian arm;
  - then the twelve genes, blocked by function exactly as on Fig. 4g, with
    the signature-sign column.
- **Two panels sharing the rows:**
  - left, **OXPHOS with MYC held fixed** (`beta_OX|MYC`);
  - right, **MYC with OXPHOS held fixed** (`beta_MYC|OX`).
- **Per row, four points**, one for each cohort-by-instrument block.
  - The MYC estimator is `M_a`.
  - `m1` is drawn filled and `m3` open, beside it. **No line joins any
    points.**
- **Bars are 95% CIs.** A vertical line marks 0. The two panels share one x
  axis, in SD per SD.
- **The legend says** that `beta_OX|MYC` is expected by construction for the
  signature genes and the scores built from them (4a.4), and that no per-gene
  FDR is applied.

**Display B: the estimator ladder for `beta_MYC|OX`** (Extended Data).
- Rows as Display A.
- Four points per row, one per MYC estimator, ordered by proliferation
  entanglement, for each cohort-by-instrument block, at `m1` and `m3`.
- **This is the figure that shows whether any MYC statement depends on how
  MYC activity is measured** (N9).

**Display C: section 9's contrast.** `BCL2L1` against `BBC3`: the two
`beta_OX|MYC` and their bootstrap difference, per cohort-by-instrument block.

## 9. Secondary, INCLUDED by the author on 2026-10-22: `BCL2L1` against `BBC3`

The abstract says OXPHOS tracks `BCL2L1` "more steeply" than `BBC3`. **That is
recorded in no note** (evening handoff, item D). **The author's reason for
including it: the statement needs proof before it stays.**

**The estimate.**
- `delta = beta_OX|MYC(BCL2L1) - beta_OX|MYC(BBC3)`, at `m1` on `M_a`, per
  cohort and instrument.
- Its 95% CI comes from **2,000 paired tumour-level bootstrap resamples**,
  seeded with `PROJECT_SEED`. Both gene models are refitted on each resample.
- **A naive difference of two separate CIs is not used.** The two models share
  their tumours.

**Reading.** "Steeper" is licensed only if `delta > 0` with its CI excluding
zero in **all four** cohort-by-instrument cells.

**Already seen.** E22's point estimates favour `BCL2L1` in all four of its
cells:
- +0.44 / +0.50 against +0.24 / +0.38 on `ox_lvl`;
- +0.36 / +0.42 against +0.20 / +0.31 on `ox_rel`.

**Why this comparison is cleaner than the signature.** Both genes were named by
the mouse, not selected on human OXPHOS, so it carries less of E34 1.2's
circularity. It does not escape it entirely: `BCL2L1`'s human OXPHOS
association was first seen in E10.

## 10. What each verdict means for Fig. 4g and the narrative, fixed now

**MYC ADDS.** The signature rises with both axes, OXPHOS more steeply. Then:
- 4g is kept as drawn. The text says MYC adds a smaller, consistent increment
  at fixed OXPHOS.
- The abstract's "even after adjustment for ... MYC" and the Doc's
  "comparable", "independent of MYC" and "lost after adjustment" sentences are
  withdrawn.
- **The mouse comparison may say "both contribute". It may NOT say "both are
  required"**, which is an interaction claim that was not tested here and is
  null where it was.

**ESTIMATOR-DEPENDENT.** MYC's contribution is not interpretable, because it is
the N9 pattern for a third time.
- The text says the signature tracks OXPHOS, and that whether MYC adds depends
  on how MYC activity is measured.
- 4g is kept with that in its legend.

**INSTRUMENT-SPLIT.** MYC adds when OXPHOS is measured one way and not the
other. The two instruments ask different questions:
- **mitoPPS asks priority** relative to the rest of the mitochondrial
  programme, blind to mass (trap 6);
- **GSVA asks level.**

If MYC adds on mitoPPS only, the reading is: **at fixed OXPHOS priority, MYC
adds; at fixed OXPHOS level, it does not.** Then:
- the Doc's "OXPHOS priority promotes a pro-apoptotic transcriptome" names its
  instrument;
- 4g's two panels show exactly that split;
- **it must be written as instrument-dependent**, never as one instrument's
  result alone.

**OXPHOS ONLY.** At fixed OXPHOS, MYC adds nothing measurable to the signature.
Then:
- **4g's visible MYC increment is leftover OXPHOS from the median split**, and
  D1 measures how much.
- The "OXPHOS rather than MYC" reading stands for the signature. 4g's legend
  says the quadrants do not hold the other axis fixed, or the panel shows the
  continuous estimates beside the quadrants.

**MIXED.** Reported as mixed, cell by cell.

**IF S5 OR S7 DISAGREES WITH S6, whatever S6's verdict:** the MYC statement
depends on which genes define the signature. **The text must name the gene
that carries it.**

**None of these licenses a functional, causal, prognostic or
treatment-selection claim.**

## 11. Expectation, stated so that it can be wrong

From E22's signs (section 2b):
- S6's `beta_MYC|OX` on `M_a` is a net of opposing genes, with `BID` positive
  and the guardian limb negative;
- the guardian limb's MYC sign is the one N9 found to flip with entanglement.

From section 2(d): part of 4g's GSVA increment is leftover OXPHOS, and less of
its mitoPPS increment is.

**The expected verdict is therefore INSTRUMENT-SPLIT, with mitoPPS adding and
GSVA not, or ESTIMATOR-DEPENDENT. It is not MYC ADDS.**

**For the variants:**
- **S5 is expected to sit BELOW S6 on `beta_MYC|OX`**, because it drops `BID`,
  which carries MYC.
- **S7 is expected to sit ABOVE S6**, because `BMF`'s `beta_MYC|OX` is negative
  at gene level and it enters with -1.

**So the variants should bracket S6.** If all three share a label, the gene
choice does not drive the answer.

This is recorded so that the result can be read against it, as E43 14.5 was.

## 12. What this cannot do

- **It is cross-sectional, on bulk transcripts.** Nothing is "primed" (N3).
- **"MYC adds at fixed OXPHOS" does not mean MYC acts independently of
  OXPHOS.** A common cause and mediation give identical partial associations
  (N6).
- **It cannot say "both are required."** That is an interaction claim, and the
  analysis does not test interaction.
- **It covers two cohorts and is not pre-registered.**

## 13. Saved, per R rules 4 and 5

- **The analysis frame**, one row per tumour. It holds:
  - `sample_id`, cohort, both OX scores and all four MYC estimators, each raw
    and Blom-scored;
  - `PROLIF`, PAM50, purity and leukocyte fraction;
  - S6, S5, S7, both arms and the twelve genes;
  - both quadrant labels.
- **The full coefficient matrix of every fit**, every term. That is 17
  readouts, giving 816 fits.
- **The full variance-covariance matrix of every fit.** This is what any later
  within-model contrast needs; E44 could only be computed because E43 kept it.
- **The full VIF vector of every fit.**
- **The fitted model objects for the 48 S6 fits only**: 2 cohorts x 2
  instruments x 4 estimators x 3 rungs. Keeping all 816 would be about 400 MB.
  Every other fit keeps its coefficients and full `vcov()`, so nothing is
  discarded that a later contrast could need.
- **A `fit_digests` table** of `coef()` and `vcov()` for every fit. **The fitted
  objects are excluded from digest comparison**, since they are not
  digest-stable.
- **Section 9's bootstrap replicates**, so the interval can be recomputed
  without refitting.
