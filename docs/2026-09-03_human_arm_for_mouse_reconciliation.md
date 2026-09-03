---
date: 2026-09-03 (rewritten, priming-first)
purpose: ONE document. What the human arm knows about APOPTOTIC PRIMING
         REGULATION BY MYC AND OXPHOS, every data source it rests on, and a
         model stated so a mouse model can be applied to it, confronted with
         it, and eventually merged into it.
supersedes: the 2026-09-03 morning version of this file, which was organised
         around the compartment split. That result is intact and is now
         section 4; priming is the centre of gravity.
sources: 2026-09-02_e10_machinery_and_priming.md,
         2026-09-02_e11_prolif_adjusted.md,
         2026-09-02_e14_curated_comparators.md,
         2026-09-02_priming_interaction_tested.md,
         2026-09-02_paper_opening_human.md,
         2026-09-02_handoff_evening.md,
         2026-09-03_e15_two_axis_gene_view.md
posture: EXPLORATORY. Nothing in the human arm is pre-registered.
---

# The human arm, for reconciliation with the mouse

## Apoptotic priming regulation by MYC and OXPHOS

**I have not read the mouse repo in preparing this.** It is not attached and
CLAUDE.md keeps it read-only. This is one side of the comparison, written so the
other side can be laid against it without either arm re-deriving what the first
one meant.

**The intended use** is that an existing mouse model, built from experiments, is
applied to human data. So this document is heavy on **what the human variables
actually are** - how each axis is defined, on what scale, with what covariate,
in how many tumours - because a model can only be applied to variables whose
construction is known. Section 2 is that inventory. Section 3 is the priming
analysis. Section 5 is a human-side model to confront the mouse one with.

---

## 0. How to read this, and the status of every number

| | |
|---|---|
| **Every number below is read from a saved object**, not from a dry run. Where a number has no script behind it, it says so |
| **E10** | re-sourced 2026-09-04. `results/machinery_and_priming.rds` is current and carries `additive_fit`, `diff_tab` with all three axis-difference forms, the `*` / `**` / `^` / border cell marks, figure 8, and section 5.0a-ii's `$gain_rank_agree`, `$gain_by_discord`, `$gain_one_cohort` and `$priming$discord` |
| **E11, E12, E13, E14** | all run by the author, unchanged since 2026-09-02 |
| **E15** | re-sourced 2026-09-04. All five figures on disk are the current ones, `fig5` included |
| **The one number here with no script behind it** | section 3.8, the cognate-pairing test. Computed in session from `results/machinery_and_priming.rds$priming`; recipe given; listed in section 10 |

**Two cohorts throughout.** Two arms of evidence only ever count when they agree;
a single cohort is a description. Where one number is given it is TCGA; where
two are given separated by `/` they are **TCGA / SCAN-B** in that order.

**Everything is proliferation-adjusted unless stated.** Partial Spearman on
`PROLIF_DISJOINT`. See 2.4 - this is not a detail, it is the reason the two axes
are comparable at all.

---

## 1. The finding, priming-first

> In two independent human breast cancer cohorts, **tumour OXPHOS status sets
> the level of individual BCL2-family transcripts, gene by gene, and does not
> set a balance between them.** The twelve transcripts span the OXPHOS axis from
> `BAD` +0.50 to `MCL1` -0.27 with pro- and anti-apoptotic members interleaved
> across the whole range, and the 35 pro/anti ratios built from them are
> **92-95% explained by an additive model of numerator and denominator identity
> with no interaction term** - so at most 5-8% of any ratio can be pair-specific.
> MYC activity does not participate: its association with the same twelve is
> much narrower, and its ordering of the wider machinery disappears entirely
> when OXPHOS is removed while OXPHOS's survives when MYC is removed. The
> `MYC x OXPHOS` interaction on these ratios was tested four ways and broke
> under three of them. What the human data do support is a **specific,
> replicating transcript configuration** - sensitisers and `BCL2L1` up, `MCL1`
> down, in OXPHOS-high tumours - and that configuration, not a balance, is what
> a functional assay should be pointed at.

Priming itself is **not measured here** and cannot be. It is a
post-translational, protein-interaction property. Section 5 says exactly which
step of the model the human arm can and cannot reach.

---

## 2. Data sources - everything a mouse model would need to be applied here

### 2.1 The two cohorts

| | TCGA-BRCA | SCAN-B |
|---|---|---|
| accession | TCGA | **GSE202203** |
| tumours | **1,095** | **3,207** |
| genes in matrix | 18,115 | 18,153 |
| PAM50 | yes | yes |
| ER call | `er_call` Positive/Negative | `er_call` 1/0 |
| **tumour purity** | yes (n = 1,007 with both) | **absent, never imputed** |
| **leukocyte fraction** | yes | **absent** |
| Luminal (LumA + LumB) | 696 | 2,436 |
| Basal | 171 | 317 |
| symbol build | current | **2014 UCSC - must be harmonised** |

Both matrices come from the snapshot `data/from_validation/`, taken from
`myc_human_validation` at pinned commit
**`d3ac60ec06c8e07df0626cbc063a52186ee167c5`** on 2026-08-31. Nothing is read
across repos at runtime.

Two scales, and they are **opposite requirements that must never share an
object**: `*_vst.rds` is log (VST) and feeds GSVA; `*_linear.rds` is DESeq2
-normalised linear and feeds mitoPPS and every per-gene value.

**SCAN-B symbols.** 19 of the 89 `OXPHOS subunits` genes carry pre-2018
ATP-synthase names. Unharmonised, the exposure covers 0.775 of the set instead
of 0.989 and Complex V loses its F1 head and c-ring. `scanb_pheno.rds$symbol_map`
is mandatory. If a mouse model is projected onto SCAN-B by symbol, **this is the
first thing that will silently break.**

### 2.2 The MYC axis

| | |
|---|---|
| **reference estimator** | **`FELSHER__MITOSTRIP`**, 61 genes |
| base set | Felsher integrative MYC signature, 67 genes |
| what `__MITOSTRIP` means | the 6 MitoCarta members removed, so the MYC axis shares **0 genes** with the OXPHOS arm, **0** with the mitoribosome and **0** with the 44 machinery genes |
| instrument | GSVA on VST, `kcdf = "Gaussian"`, scored in one batch per cohort |
| proliferation entanglement | 14.8% of the base set is HALLMARK E2F/G2M |
| source | `data/genesets_from_library_human/felsher_integrative_signature.csv`, `mammary_geneset_library` v1.0, commit `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b`, 2026-06-20 |

**Two secondary estimators are used wherever an estimator-dependence check is
made**: `MYC_UP.V1_UP__FULL` (198 genes, **1.5%** entanglement - the least
entangled available) and `HALLMARK_MYC_TARGETS_V1__FULL` (200, 23.5%). The full
compendium is 16 signatures x 4 strip variants = 72 sets in
`results/set_definitions.rds$myc_sets`, spanning **1.5% to 47.6%** proliferation
entanglement. A MYC result that does not survive the low-entanglement estimator
is a proliferation result.

**`M_b`, the CollecTRI regulon estimator, is tabled and never plotted** for the
machinery: `M_b__MITOSTRIP` contains 7 of the 44 genes.

**MYC mRNA is not MYC activity.** `rho(log2 MYC, OXPHOS subunits) = -0.032`
against `+0.388` for the activity signature. This is N7 and it matters for any
mouse model expressed in terms of *Myc* transcript.

### 2.3 The OXPHOS axis

| | |
|---|---|
| **reference arm** | **`OXPHOS subunits`, 89 genes** - the narrow set, subunits only |
| instrument | GSVA on VST, same batch as the MYC axis |
| NOT used as the axis | MitoCarta's `OXPHOS` umbrella, which adds assembly factors |
| mtDNA-encoded genes | held in a **separate synthetic arm** of 13 `MT-` genes and never pooled with the nuclear subunits |
| source of the gene list | MitoCarta 3.0 MitoPathways sheet |
| other instruments available | mitoPPS (linear), `content` (log2 summed), `zmean` - all four in `results/scanb_scores.rds` and the TCGA snapshot |

`rho(MYC axis, OXPHOS axis)` after proliferation adjustment = **0.320 / 0.256**.
The two axes are correlated, which is why every claim below is a *conditioning*
claim and not a comparison of two marginal correlations.

**One self-overlap to carry forward:** `CYCS` is 1 of the 89. Its OXPHOS value is
partly self-correlation and it is marked on every panel.

### 2.4 The covariate, and what "adjusted" means

| | |
|---|---|
| primary covariate | **`PROLIF_DISJOINT`, 318 genes** |
| why this one | 0 genes shared with `FELSHER__MITOSTRIP`, 0 with the 89 OXPHOS-arm genes, 0 with the 83 mitoribosome genes. **The same adjustment can therefore be applied to both axes**, which is what makes them comparable |
| secondary | `PROLIF_STD`, 327 genes - shares 9 with the MYC reference, so it is carried beside and never alone |
| method | partial Spearman = Pearson on residuals of ranks |
| the two corrections disagree | partialling on a score and stripping genes from the estimator **disagree in sign**. Both are computed everywhere; neither is reportable without the other |
| self-adjustment | `TP53` and `BIRC5` are IN the covariate. Marked on every figure; every summary given with and without them |
| TCGA-only extra | purity + leukocyte fraction, n = 1,007 |

### 2.5 The priming gene set and the ratio grid

**The author's lists**, 7 pro-apoptotic x 5 anti-apoptotic = **35 ratios**, no
self-pair, no empty cell.

| side | genes |
|---|---|
| pro-apoptotic (numerator) | `BCL2L11` `BMF` `PMAIP1` `BBC3` `BID` `BAD` `BIK` |
| anti-apoptotic (denominator) | `BCL2` `BCL2L1` `MCL1` `BCL2L2` `BCL2A1` |

- **A ratio is `log2(pro) - log2(anti)`** on `log2(linear DESeq2-normalised + 1)`,
  i.e. a difference of logs, correlated by Spearman, partialled on
  `PROLIF_DISJOINT`.
- **`BCL2L2` was corrected 2026-09-02.** It was supplied on both lists; Bcl-w is
  a multidomain guardian, not a BH3-only sensitiser, and the author confirmed
  the pro-side entry was an error. It is anti-apoptotic only.
- **All 12 clear the 25th expression percentile in both cohorts**, so none of
  this is a low-expression artefact. (`HRK` and `FASLG` in the wider 44 do not,
  and are starred everywhere.)
- **Compartment: 11 of the 12 are MitoCarta MOM; only `BMF` is absent.** So
  **30 of the 35 ratios are MOM/MOM** and cancel mitochondrial abundance by
  construction. This is the basis of the E13 control (3.4).
- If a seventh numerator is ever wanted, `HRK` and `BOK` are already scored -
  `HRK` would need its low-expression flag carried through.

### 2.6 The wider death axis

- **44 canonical machinery genes** = the members of the curated death set
  carrying a `family_pathway` label, i.e. the genes that *are* apoptotic
  machinery rather than curated alongside it. 20 in MitoCarta, 24 not.
  Modules: mitochondrial/intrinsic 21, death receptor/extrinsic 11, IAP/NF-kB 7,
  effector caspase 3, other 2.
- **1,232-gene curated set** with `effect` (pro-death 512 / pro-survival 587),
  `pathway` (apoptosis / CICD) and `family_pathway` columns:
  `data/genesets_celldeath_human/cell_death_genes_consolidated.csv`, from
  `myc_mouse` at commit `6a9c7dd513800a2a433934314a87d161ce98caa2`, **human-native
  columns, never remapped**.
- **15 Tang modalities**, Tang et al. 2024, **doi:10.1016/j.csbj.2024.08.012**.

### 2.7 Pinned external catalogues

| catalogue | pin | used for |
|---|---|---|
| **Human MitoCarta 3.0** (Rath et al. 2021, *NAR*) | MD5 `3c0bd24e362238216e142bc708e41286`, 1,136 symbols, 154 MitoPathways | regulon membership, sub-compartment (MOM/IMS/MIM), the 30-pathway ladder, the OXPHOS arm |
| **MSigDB** via `msigdbr` | **2026.1.Hs** | Reactome C2 CP:REACTOME comparator programmes; Hallmark E2F/G2M for entanglement |
| **CollecTRI** (OmniPath) | SHA-256 `d72660703f7cc...1d287d` | `M_b` estimator; **and the untouched next analysis, section 10** |
| **Menegollo / Bentham et al. 2024**, *Cancer Res* CAN-23-3172 | commit `8fdbb3437ae5537055d5d5429411bdb3b333c04a` | the study's starting point; the forkscale/MCbiclust axis is phase 2 and out of scope |
| **`myc_human_validation`** | frozen at `d3ac60e` | the pre-registered sibling study - N2 |

### 2.8 Scripts, saved objects, tables and figures

| script | what it is for | saved object |
|---|---|---|
| `E10_machinery_measures_and_priming.R` | **the priming script.** 12 components, 35 ratios, gain, strata, additive fit, cell marks | `results/machinery_and_priming.rds` |
| `E11_prolif_adjusted_machinery.R` | the 44 under adjustment; matched nulls; the conditioning ladder; bootstrap | `results/prolif_adjusted_machinery.rds` |
| `E12_borrowing_explainer.R` | decomposes MYC's per-gene correlation into borrowed and own | `results/borrowing_explainer.rds` |
| `E13_priming_and_content.R` | **is the OXPHOS association mitochondrial mass?** MOM/MOM ratios + random-pair null | `results/priming_and_content.rds` |
| `E14_curated_comparators.R` | named non-apoptotic membrane-spanning comparators; infiltrate falsifier | `results/curated_comparators.rds` |
| `E15_two_axis_gene_view.R` | display only; the 44 on both axes, per gene | `results/two_axis_gene_view.rds` |

**Objects a mouse-model application would read first**, by name:

| object | rows | what it holds |
|---|---|---|
| `machinery_and_priming.rds$component_cor` | 144 | **the 12 genes x 2 axes x 2 cohorts, with 95% Fisher-z CIs.** The primary human priming table |
| `$priming` | 420 | the 35 ratios x axes x cohorts: `rho`, `ci_lo/hi`, `rho_pro`, `rho_anti`, `coexpr`, `best_component`, `gain`, `discord` |
| `$priming_strata` | 1,260 | the same, split `all` / `Luminal` / `Basal` |
| `$additive_fit` | 4 | the R2 of `rho ~ pro + anti` per axis per cohort |
| `$gain_summary` | 70 | per ratio, `gain` in each cohort and `both_positive` |
| `$gain_rank_agree` | 2 | per axis, the both/one/neither split and the **cross-cohort rank agreement of `gain`** (3.3a) |
| `$gain_by_discord` | 8 | gain split by whether the two components lean opposite ways (3.3b) |
| `$gain_one_cohort` | 19 | the ratios that gain in exactly one cohort - **read with 3.3a, never as winners** |
| `$diff_tab` | 70 | the three ways to collapse the two axes into one number, per ratio |
| `$gene_lean` | 12 | per gene, `rho(OXPHOS) - rho(MYC)` per cohort and whether it is OXPHOS-led |
| `$between_test` | 140 | pooled value against its own Luminal and Basal values |
| `priming_and_content.rds$ratio_rho` | 70 | ratios with **`ox_condM`** and **`myc_condO`** - the conditioned versions |
| `priming_and_content.rds$mom_null` | 2 | the 500-draw random MOM/MOM null |
| `borrowing_explainer.rds$decomp` | 2 | SD of the OXPHOS column, the MYC column, the borrowed part and MYC's own |

**Tables on disk:** `outputs/tables/E10_priming_ratios.csv`,
`E10_canonical_machinery.csv`, `E11_gene_rho_by_adjustment.csv`,
`E14_comparator_splits.csv`, `E14_mitocarta_pathway_ladder.csv`,
`E15_gene_rho_two_axes.csv`.

**Priming figures:** `E10_fig3` the 35-cell ratio heatmap; `E10_fig4` does a
ratio beat its parts; `E10_fig5` the 12 components with intervals; `E10_fig6`
the heatmap by compartment; `E10_fig7` the 12 by compartment; `E10_fig8` the two
axes collapsed to one number, three ways; `E13` the MOM/MOM control.

---

## 3. The priming analyses, in the order they answer each other

### 3.1 The 12 transcripts on their own - the primary human table

Partial Spearman on `PROLIF_DISJOINT`. Ordered by the cross-cohort mean on
OXPHOS. **This is the table to map a mouse model onto.**

| gene | side | OXPHOS TCGA | OXPHOS SCAN-B | mean | MYC TCGA | MYC SCAN-B | mean |
|---|---|---|---|---|---|---|---|
| **`BAD`** | pro | **+0.525** | **+0.482** | **+0.503** | +0.162 | +0.114 | +0.138 |
| **`BCL2L1`** | **anti** | **+0.426** | **+0.350** | **+0.388** | -0.030 | +0.033 | +0.002 |
| `BIK` | pro | +0.338 | +0.213 | +0.276 | +0.008 | +0.051 | +0.029 |
| `BBC3` | pro | +0.365 | +0.170 | +0.267 | +0.163 | +0.110 | +0.136 |
| `BID` | pro | +0.329 | +0.091 | +0.210 | **+0.357** | **+0.243** | **+0.300** |
| `BCL2A1` | anti | -0.004 | -0.107 | -0.056 | +0.095 | +0.014 | +0.055 |
| `BCL2L2` | anti | -0.145 | -0.008 | -0.077 | -0.073 | -0.080 | -0.076 |
| `PMAIP1` | pro | -0.147 | -0.064 | -0.105 | **-0.194** | **-0.223** | **-0.209** |
| `BCL2` | anti | -0.160 | -0.074 | -0.117 | -0.222 | -0.160 | -0.191 |
| `BCL2L11` | pro | -0.314 | -0.025 | -0.170 | -0.162 | -0.039 | -0.101 |
| `BMF` | pro | -0.200 | -0.223 | -0.212 | -0.150 | -0.079 | -0.115 |
| **`MCL1`** | **anti** | **-0.241** | **-0.292** | **-0.266** | +0.081 | -0.037 | +0.022 |

**Three readings, and they are the foundation of everything after.**

1. **Neither side moves as a block.** The most positive and the most negative
   anti-apoptotic genes are `BCL2L1` (+0.39) and `MCL1` (-0.27); the most
   positive and most negative pro-apoptotic are `BAD` (+0.50) and `BMF` (-0.21).
   The functional dichotomy the ratios assume is **not present at gene level.**
2. **The OXPHOS range is 0.77 wide; the MYC range is 0.51 and is narrower in
   both cohorts.** The five genes that are *individually* OXPHOS-led -
   `|rho(OXPHOS) - rho(MYC)| >= 0.20` in at least one cohort - are `BCL2L1`
   (+0.455 / +0.317), `BAD` (+0.363 / +0.368), `BIK` (+0.330 / +0.163), `BBC3`
   (+0.202 / +0.061) and **`MCL1` (-0.321 / -0.254, leaning hard to OXPHOS while
   running the other way)**. The other seven are not.
3. **The two genes MYC holds on to** are `BID` (+0.300) and `PMAIP1` (-0.209),
   and they are the only two whose MYC association survives conditioning on
   OXPHOS in both cohorts under all three estimators. **They move in opposite
   directions**, and both are pro-apoptotic BH3-only proteins. That is two
   genes, not a shift in balance.

Source: `machinery_and_priming.rds$component_cor`; figure `E10_fig5`.

### 3.2 The 35-ratio grid, and the ceiling on what a ratio can carry

The question a ratio is supposed to answer is whether the *balance* moves. It
does not carry that information here.

**Fit `rho ~ numerator identity + denominator identity`, no interaction term,
over all 35 ratios:**

| axis | cohort | R2 additive | numerator only | denominator only |
|---|---|---|---|---|
| MYC | TCGA | **0.946** | 0.642 | 0.305 |
| MYC | SCAN-B | **0.942** | 0.740 | 0.202 |
| OXPHOS | TCGA | **0.939** | 0.729 | 0.210 |
| OXPHOS | SCAN-B | **0.921** | 0.577 | 0.345 |

**92.1 to 94.6% of every ratio in the grid is predicted by which two genes it is
made of, with no interaction term at all.** At most **5-8%** of a ratio can be
pair-specific, and that is the hard ceiling on anything "priming balance" could
mean in this data. The numerator carries about twice the denominator.

This regenerates from `E10` section 5.0b (`$additive_fit`) as of 2026-09-03. It
is the single most important number in the priming section and **the sentence to
keep if only one survives editing.**

### 3.3 Does a ratio beat its own two genes?

`gain` = `|rho of the ratio|` minus `|rho of its stronger component|`. The
falsifier was written before the numbers: **a ratio is only worth reporting as a
ratio where `gain > 0` in both cohorts.**

| axis | n | `gain > 0` in both | the ratios |
|---|---|---|---|
| **OXPHOS** | 35 | **0** | none |
| **MYC** | 35 | **5** | `BBC3/BCL2` (+0.052/+0.062), `BAD/BCL2` (+0.040/+0.050), `BAD/BCL2L1` (+0.025/+0.003), `BBC3/BCL2L2` (+0.002/+0.019), `BBC3/BCL2L1` (+0.006/+0.015) |

**On the OXPHOS axis not one ratio of the 35 adds information to its stronger
component in both cohorts.** The strongest-looking cells are the components
talking: the best OXPHOS ratio is `BAD/MCL1` at +0.507 and `BAD` alone is +0.525.

The five that gain on MYC are all under **+0.07**, and they reshuffled
completely when the proliferation covariate was added - none of them is one of
the three the unadjusted run named. **A list that changes identity under a
covariate is a description of noise; the safest thing is not to name them.**

#### 3.3a What that zero does and does not mean

The both-cohorts rule is what decides what gets reported and it does not move.
But `0` and *nothing there* are different claims, and three things behind the
zero are themselves reproducible.

**First, the distribution, not just its sign.** Splitting the 35 ratios by how
many cohorts they gain in:

| axis | gain in both | gain in one | gain in neither | mean gain (T / S) | best single cell |
|---|---|---|---|---|---|
| **OXPHOS** | 0 | **10** | 25 | -0.116 / -0.110 | **+0.090** |
| **MYC** | 5 | 9 | 21 | -0.050 / -0.031 | +0.062 |

Two readings run in opposite directions and both are true. The OXPHOS
distribution sits **two to three times as far below zero** as the MYC one -
ratios cost more on the OXPHOS axis, which is the negative result. And yet **the single
largest gain anywhere in the grid is on OXPHOS, not MYC** (+0.090 against
+0.062). The zero in the table is a statement about reproducibility, not about
size.

**Second, `gain` is a stable quantity that happens to be negative.** The
*ordering* of the 35 ratios by gain agrees between the cohorts even though the
sign does not: Spearman **+0.408** on OXPHOS, **+0.560** on MYC. If gain were
noise around zero the ordering would not survive at all. It is not that the
measurement is unstable; it is that the whole distribution is shifted down and
the top of a reproducible ordering brushes zero one cohort at a time.

**Third, the near miss is very near.** `BAD/MCL1`, the strongest OXPHOS priming
ratio in the grid (+0.507 / +0.467), fails the test by **-0.018 / -0.015** -
a margin of about 3% of the correlation it is made of. The
verdict "the ratio is its components" is right, but it is right by a whisker for
the cell that matters most, and that is worth saying out loud in the paper
rather than only reporting the 0.

**The ten that gain in exactly one cohort** (proliferation-adjusted; `rho` and
`best component` are from the cohort that gains):

| ratio | gains in | gain there | gain in the other | ratio rho | best component | pro / anti rho |
|---|---|---|---|---|---|---|
| `BCL2L11/BCL2L1` | TCGA | **+0.090** | -0.176 | -0.515 | 0.426 | -0.314 / **+0.426** |
| `BID/MCL1` | TCGA | **+0.077** | -0.033 | +0.406 | 0.329 | +0.329 / **-0.241** |
| `BID/BCL2A1` | SCAN-B | +0.045 | -0.179 | +0.152 | 0.107 | +0.091 / -0.107 |
| `BMF/BCL2L1` | SCAN-B | +0.042 | -0.009 | -0.392 | 0.350 | -0.223 / +0.350 |
| `BBC3/MCL1` | TCGA | +0.040 | -0.017 | +0.405 | 0.365 | +0.365 / -0.241 |
| `BIK/MCL1` | TCGA | +0.027 | -0.007 | +0.365 | 0.338 | +0.338 / -0.241 |
| `BIK/BCL2L2` | TCGA | +0.012 | -0.105 | +0.350 | 0.338 | +0.338 / -0.145 |
| `BBC3/BCL2L2` | TCGA | +0.009 | -0.080 | +0.374 | 0.365 | +0.365 / -0.145 |
| `BBC3/BCL2` | SCAN-B | +0.008 | -0.020 | +0.179 | 0.170 | +0.170 / -0.074 |
| `BBC3/BCL2A1` | SCAN-B | +0.008 | -0.168 | +0.178 | 0.170 | +0.170 / -0.107 |

**None of these is a result and none should be quoted as one.** Six of the ten
lose in the other cohort by more than they win by in this one, and the top of
the list is the worst offender: `BCL2L11/BCL2L1` is the largest gain in the grid
and the largest loss on the list.

But the list is not arbitrary either, which is the point of writing it down:

#### 3.3b Which cells can gain at all, and why the OXPHOS ones fail

**27 of the 29 positive gains in the grid belong to a ratio whose numerator and
denominator lean OPPOSITE ways on the axis - and on the OXPHOS axis it is 10 of
10.** Splitting the 35 ratios by whether `sign(rho_pro) != sign(rho_anti)`:

| axis | cohort | components | n | n with gain > 0 | mean gain | best |
|---|---|---|---|---|---|---|
| OXPHOS | TCGA | same direction | 16 | **0** | -0.193 | -0.044 |
| OXPHOS | TCGA | **opposite** | 19 | 6 | -0.052 | +0.090 |
| OXPHOS | SCAN-B | same direction | 16 | **0** | -0.159 | -0.016 |
| OXPHOS | SCAN-B | **opposite** | 19 | 4 | -0.069 | +0.045 |
| MYC | TCGA | same direction | 17 | **0** | -0.093 | -0.010 |
| MYC | TCGA | **opposite** | 18 | 9 | -0.009 | +0.052 |
| MYC | SCAN-B | same direction | 17 | 2 | -0.057 | +0.015 |
| MYC | SCAN-B | **opposite** | 18 | 8 | -0.006 | +0.062 |

The two exceptions are both `MYC` in SCAN-B (`BBC3/BCL2L1` +0.015,
`BAD/BCL2L1` +0.004) and both have `BCL2L1` as the denominator, whose MYC rho
there is **+0.033** - close enough to zero that calling the pair concordant is a
coin flip on a sign, and small enough that the gains are trivial.

Half of this is arithmetic and must be labelled as such: dividing two things
that move together subtracts the shared signal, dividing two things that move
apart adds them. Nobody should be surprised.

**The biology is in which pairs are discordant.** On the OXPHOS axis the
discordant block is essentially the *pro-death-up / anti-death-down* quadrant -
`{BAD, BBC3, BIK, BID} x {MCL1, BCL2, BCL2L2, BCL2A1}`, plus the three cells
where the two down-going sensitisers meet `BCL2L1`. **The OXPHOS axis splits the
BCL2 family rather than lifting it**, which is 3.1 and section 5 restated at the
level of pairs, and it is why 19 of 35 OXPHOS cells are even eligible to gain.

So the failure of the ratios on the OXPHOS axis is **not** "OXPHOS does not
reach priming". It is narrower and more useful: *within a discordant pair, the
axis's effect is already fully carried by whichever of the two genes it moves
more.* That is the additive fit of 3.5 (92.1-94.6%) stated one pair at a time,
and it is why the components, not the ratios, are what section 5's model is
built on.

**Why the two best OXPHOS cells do not replicate is checkable, not mysterious.**
`BID/MCL1` and `BBC3/MCL1` are textbook priming ratios - sensitiser up,
guardian down - and both gain in TCGA. Their denominator is stable across
cohorts (`MCL1` -0.241 / -0.292); their **numerators are not**. `BID` falls from
**+0.329 to +0.091** and `BBC3` from **+0.365 to +0.170** between TCGA and
SCAN-B. The pair loses its discordance strength on the numerator side, so the
ratio has nothing left to add. **That is a falsifiable statement**: a third
cohort in which `BID` and `BBC3` correlate with OXPHOS at TCGA strength should
show `BID/MCL1` and `BBC3/MCL1` gaining; one in which they behave like SCAN-B
should not. Nothing about the ratio itself needs to be invoked.

Figure `E10_fig4` draws this as a scatter against the diagonal. Figures
`E10_fig3` and `E10_fig6` carry the same test onto the heatmap cells:

| mark | meaning |
|---|---|
| `*` | in this cohort, `gain > 0` **and** `\|rho\| >= 0.30` |
| `**` | that cell's two genes are **each** OXPHOS-led on their own (3.1 reading 2) |
| `^` on an axis label | that gene passes the `**` gene test |
| heavy border | the `*` test passes in **both** cohorts |

**Not one MYC cell is marked**, in either heatmap, at any stratum, in either
cohort. All 22 marks on `E10_fig6` and all 7 on `E10_fig3` are OXPHOS. The
0.30 and 0.20 lines are hand-drawn and are not tests.

### 3.4 Is the OXPHOS association just mitochondrial mass? (E13)

The most obvious deflationary reading of everything above is *"the OXPHOS score
is a mitochondrial content index and the BCL2 family is along for the ride."*
**The priming ratios can test that and the 44-gene analysis could not**, because
11 of the 12 are MOM: a MOM/MOM ratio cancels anything that scales with how many
mitochondria a tumour has.

| | TCGA | SCAN-B |
|---|---|---|
| mean `\|rho(OXPHOS)\|` of the 30 MOM/MOM ratios | **0.246** | **0.147** |
| 500 random MOM/MOM pairs, null mean | 0.200 +/- 0.039 | 0.153 +/- 0.034 |
| z | **+1.18** | **-0.16** |
| percentile of draws below | 0.878 | 0.480 |

**Two conclusions, and the second is the one people forget.**

1. **The ratios are not flat, so the OXPHOS axis is not a mitochondrial mass
   index.** It resolves *within* a single compartment.
2. **But random MOM/MOM pairs are not flat either, and the BCL2-family excess
   does not replicate** (z +1.18 in TCGA, -0.16 in SCAN-B). **The
   within-compartment resolving power belongs to the axis, not to the death
   machinery.** Whatever separates `BAD` from `MCL1` on the OXPHOS axis is a
   general property of that axis and has not been shown to be about apoptosis.

Source: `priming_and_content.rds$mom_null`; figure `E13`.

### 3.5 Conditioning, at the ratio level

The two axes correlate at 0.320 / 0.256, so a marginal comparison cannot say
which one carries a ratio. Conditioning can.

| | TCGA | SCAN-B |
|---|---|---|
| mean `\|rho(ratio, OXPHOS)\|` | 0.229 | 0.142 |
| ...with MYC removed | 0.208 | 0.135 |
| **kept** | **91%** | **95%** |
| mean `\|rho(ratio, MYC)\|` | 0.142 | 0.105 |
| ...with OXPHOS removed | 0.113 | 0.087 |
| **kept** | **79%** | **83%** |

And at the top of the range the OXPHOS ratios barely move at all:
`BCL2L11/BCL2L1` -0.515 -> **-0.511**; `BAD/MCL1` +0.507 -> **+0.508**;
`BMF/BCL2L1` -0.416 -> -0.396. MYC's equivalents fall: `BID/MCL1` +0.222 ->
+0.106; `BIK/MCL1` -0.021 -> **-0.156** (it changes sign).

**But note the asymmetry is much weaker at ratio level than at gene level**, and
that is not an inconsistency - it is a property of ratios, and it has an
explanation:

> **A ratio subtracts two transcripts. What the two share cancels - and what
> they mostly share is the OXPHOS-borrowed component, because that component is
> a scaled copy of the OXPHOS column and hits pro- and anti-apoptotic genes
> alike. What does not cancel is MYC's own component, which is the larger of the
> two and is unsorted. So a ratio is an enrichment device for the non-shared
> part of an association.**

That is `priming_interaction_tested.md` C2 and it is why MYC survives in the
ratios and not in the genes. **It is a property of ratios, not evidence of a
priming-specific mechanism**, and any mouse-side ratio analysis will inherit it.

Cross-cohort replication of the 35 ratio values: OXPHOS **0.826** raw / 0.800
conditioned; MYC **0.930** / 0.857 (Spearman over the 35).

Sources: `priming_and_content.rds$ratio_rho`; `borrowing_explainer.rds$decomp`.

### 3.6 `MYC x OXPHOS` on the ratios - tested four ways, not supported

The hypothesis was worth testing: *OXPHOS sets the level, and the interaction
sets where the balance tips.*

**The encouraging first look.** Fitting `ratio ~ prolif + MYC + OXPHOS +
MYC:OXPHOS` over the 35: 9 of 35 interactions clear `|t| > 2` in TCGA and 11 of
35 in SCAN-B, the 35 coefficients correlate **r = 0.49** across cohorts, and five
replicate with the same sign - `BID/BCL2`, `BMF/BCL2`, `BBC3/BCL2`, `BAD/BCL2`
(positive) and `PMAIP1/MCL1` (negative). Four of the five share `BCL2`.

**Then the falsifiers.**

| test | result |
|---|---|
| **F1 curvature.** A product term absorbs nonlinearity in either main effect. Refit with natural splines on both | **FAILS.** Linear and spline interaction coefficients correlate **0.125** in TCGA - which ratios carry it changes almost completely. Three of the five collapse: `BID/BCL2` 0.081 -> -0.003, `BAD/BCL2` 0.059 -> 0.012, `BMF/BCL2` 0.065 -> 0.019 |
| **F2 estimator.** Repeat with `MYC_UP.V1_UP` (1.5% entanglement) | **FAILS.** Cross-cohort replication **reverses**: r = **-0.33** on the linear term, -0.42 on ranks |
| **F3 tertiles.** No model: split tumours by OXPHOS tertile, ask whether MYC's association with the ratio strengthens | **FAILS.** TCGA mean rho by tertile 0.010 / 0.003 / **0.031** - non-monotone, dipping in the middle. SCAN-B moves the other way (high-minus-low -0.008 against TCGA +0.021) |
| **F4 prior.** The sibling **pre-registered** study tested `MYC x OXPHOS` on *functional* apoptotic priming and found it null | Different endpoint, so neither confirms the other. **Nothing here overturns that null and nothing may be written as if it did** |

**Verdict: the apparent interaction is curvature in the main effects.** The TCGA
mid-tertile dip is exactly what a product term misreads as an interaction, and
it is why F1 breaks it.

Source: `priming_interaction_tested.md` H1.

### 3.7 Subtype - and the only three ratios that survive everything

Pooled values are not within-subtype values. **15 of 35 OXPHOS ratios in TCGA
and 19 of 35 in SCAN-B sit outside the range of their own Luminal and Basal
values** - i.e. they are reading a difference *between* subtypes.

Cell marks by stratum, on the OXPHOS axis (`E10_fig6`):

| stratum | marked (`*`) TCGA | SCAN-B | marked in **both** |
|---|---|---|---|
| pooled `all` | 7 | 4 | **0** |
| Luminal | 4 | 3 | **0** |
| **Basal** | 6 | 11 | **3** |

**Only three cells in the entire 420-cell grid pass in both cohorts, and all
three are Basal, all three OXPHOS:**

| ratio | TCGA (n = 171) | SCAN-B (n = 317) | numerator | denominator |
|---|---|---|---|---|
| **`BBC3/BCL2`** | **+0.494** [0.36, 0.61] | **+0.316** [0.21, 0.41] | `BBC3` +0.412 / +0.200 | `BCL2` -0.367 / -0.257 |
| **`BID/BCL2`** | **+0.489** [0.36, 0.60] | **+0.312** [0.21, 0.41] | `BID` +0.453 / +0.189 | `BCL2` -0.367 / -0.257 |
| **`BBC3/MCL1`** | **+0.446** [0.31, 0.56] | **+0.301** [0.19, 0.40] | `BBC3` +0.412 / +0.200 | `MCL1` -0.211 / -0.222 |

`BBC3/MCL1` also carries `**` in both cohorts - both its genes are individually
OXPHOS-led - and is **the only cell in the figure carrying all three marks.**

**Read this against where the marks are densest.** Basal is 171 TCGA and 317
SCAN-B tumours, so it is also where the intervals are widest (about +/- 0.15 in
TCGA). Three cells out of 70 in the stratum with the least power is a
hypothesis, not a result. **But it is the sharpest hypothesis the priming side
produces, and if a mouse model is basal-like this is where the comparison should
be made.**

Note in the other direction: the **compartment split of the wider 44 on the MYC
axis is higher inside subtypes than pooled** (TCGA 0.187 pooled -> 0.298 Luminal
-> 0.323 Basal; SCAN-B 0.137 -> 0.295 -> 0.194), while OXPHOS's is flat to
slightly lower (0.453 -> 0.464 -> 0.349; 0.489 -> 0.442 -> 0.367). Every one of
those MYC values is still inside its compartment-matched null (z 0.12-0.70).

Sources: `$priming_strata`, `$between_test`, `$strata_marked`; `E10_fig6`.

### 3.8 A sub-model tested and rejected: binding-specificity pairing

**The obvious mechanistic refinement**, and it would have been a strong result:
*OXPHOS status moves each guardian together with its own cognate sensitiser, so
a ratio against the cognate partner cancels and a ratio against a non-cognate
partner adds.* Canonical selectivity (Certo et al. 2006; Chen et al. 2005):
`BAD` and `BMF` -> BCL2 / BCL-xL / BCL-w; `PMAIP1` -> MCL1 / A1; `BIK` -> BCL-xL
/ BCL-w. `BCL2L11`, `BBC3` and `BID` are promiscuous and excluded.

| | n | mean `gain` TCGA | mean `gain` SCAN-B |
|---|---|---|---|
| cognate pairs | 10 | **-0.116** | **-0.174** |
| non-cognate pairs | 10 | -0.109 | -0.102 |

**Not supported.** The difference is negligible in TCGA and points the *wrong*
way in SCAN-B. **The transcript configuration is not organised by binding
specificity**, which removes one of the more attractive bridges from a
transcript pattern to a functional one.

**This number has no script behind it.** Computed in session on 2026-09-03 from
`machinery_and_priming.rds$priming` filtered to `axis == "OXPHOS"` and
`pro %in% c("BAD","BMF","PMAIP1","BIK")`, grouping `gain` by whether `anti` is
in the sensitiser's cognate list. Listed in section 10 to be folded into `E10`.

### 3.9 What survives which filter

| filter | survivors |
|---|---|
| all 35 ratios | 35 |
| ...`gain > 0` in both cohorts, OXPHOS | **0** (but 10 in one cohort, and the best single cell in the grid is one of them - 3.3a) |
| ...`gain > 0` in both cohorts, MYC | 5, all under +0.07, all reshuffled by the covariate |
| ...components lean opposite ways, OXPHOS | 19 of 35 - **the only cells that ever gain**, 10 of 10 (3.3b) |
| ...`gain > 0` AND `\|rho\| >= 0.30`, pooled | 7 cells, all OXPHOS, **0 in both cohorts** |
| ...same, any stratum | 22 cells, all OXPHOS |
| ...same, in **both cohorts** | **3, all Basal OXPHOS** |
| ...and both genes individually OXPHOS-led | **1: `BBC3/MCL1` in Basal** |

Pooled, a priming ratio is **either replicable or large, never both**: the five
that beat their parts in both cohorts are all MYC and top out at `|rho|` 0.27,
while the seven that reach 0.30 are all OXPHOS and none beats its parts twice.

---

## 4. The context: the machinery-wide result, condensed

The priming genes are 12 of a wider 44, and the wider result is what licenses
reading the priming table as an OXPHOS phenomenon at all. Full treatment:
`e11_prolif_adjusted.md`, `e14_curated_comparators.md`, `e15_two_axis_gene_view.md`.

### 4.0 The claim ladder, for reference

The candidate claim was *"the apoptotic machinery is driven, or at least more
strongly correlated, by OXPHOS than directly by MYC."* Five rungs; section 8
refers to them by number.

| # | statement | status | what licenses it |
|---|---|---|---|
| 1 | The 44 genes correlate more strongly with OXPHOS than with MYC activity | **SUPPORTED** | bootstrap SD ratio **1.640 [1.473, 1.784]** and **1.558 [1.461, 1.650]**; holds inside Luminal and Basal separately |
| 2 | ...and it is not proliferation | **SUPPORTED** | unchanged by `PROLIF_DISJOINT`, by the `__PROLIFSTRIP` / `__BOTHSTRIP` estimators, and in TCGA by purity + leukocyte fraction on top |
| 3 | ...and it is OXPHOS **rather than** MYC | **SUPPORTED - the strong result** | the conditioning asymmetry in 4.1. **The only rung that required a test rather than a description, and the only one whose falsifier was written before the answer was seen** |
| 4 | ...and it is specific to apoptosis | **SPLIT IN TWO** | 4a NOT supported for the regulon half (4.3, first paragraph); 4b SUPPORTED for the cytosolic half (4.3, the table) |
| 5 | "driven" | **NOT TESTED, and untestable here** | cross-sectional bulk tumour RNA. Every causal word must come from the mouse arm |

**Rung 3 is the one to carry into the mouse comparison**, and section 5's model
is what rung 3 plus the priming analyses add up to.

### 4.1 The dissociation, and the conditioning asymmetry

The 44 per-gene correlations spread **1.6x wider** on OXPHOS than on MYC (SD
0.313 / 0.245 against 0.191 / 0.157; tumour-level bootstrap SD ratio
**1.640 [1.473, 1.784]** and **1.558 [1.461, 1.650]**). But the axes correlate,
so the test is conditioning:

| | TCGA | SCAN-B |
|---|---|---|
| OXPHOS split on MitoCarta membership | 0.453 | 0.489 |
| ...**with MYC removed** | **0.485** | **0.525** |
| MYC split | 0.187 | 0.137 |
| ...**with OXPHOS removed** | **-0.043** | **-0.058** |

**MYC's ordering was inherited; OXPHOS's was not.** All five bootstrap contrasts
exclude their null in both cohorts.

**The control that makes this mean anything** (`E11_fig2`): after the same
adjustment, MYC still tracks the **mitoribosome** at z = 8.66 / 6.54 against its
matched null, while the machinery sits at z = 0.67 / 0.13. The MYC axis was not
emptied - it simply does not order these genes.

**And MYC's own component is real.** Decomposing MYC's per-gene correlation:
the borrowed part (`rho(MYC,OXPHOS) x rho(OXPHOS,gene)`) has SD **0.100 / 0.063**;
what remains after conditioning has SD **0.149 / 0.135**, far above the ~0.03
sampling noise. **It is larger than the borrowed part and it is unsorted** - its
gap across the compartment split is -0.019 / +0.004. *Something* orders MYC's
associations with these genes and this study has not asked what.

### 4.2 What is ordered is regulon membership

MitoCarta membership predicts a gene's position on the OXPHOS axis at
**0.45-0.49**; its annotated direction of effect at only **0.22-0.28**.

**Read "mitochondrial" as co-regulation, not localisation.** MitoCarta is a
proteome catalogue and half the BCL2 family translocates. At transcript level
membership marks **membership of the nuclear-encoded mitochondrial regulon**.
This matters for the mouse step: the claim this data type supports is
transcriptional co-regulation, not protein localisation.

### 4.3 The cytosolic anomaly - the one apoptosis-specific thing

The mitochondrial half of the 44 is **not** special: at its expression- and
sub-compartment-matched null (z +0.48 / +0.62), 63rd percentile of 30 MitoCarta
leaf pathways, and *below* mitophagy's mitochondrial half.

The cytosolic half is. Mean rho with OXPHOS against its own null:

| programme | cytosolic half TCGA | SCAN-B |
|---|---|---|
| **apoptotic machinery** | **-0.106** (z **-1.19**) | **-0.094** (z **-1.39**) |
| mitophagy | +0.143 (z +2.47) | +0.096 (z +1.72) |
| Fe-S cluster assembly | +0.118 (z +1.71) | +0.159 (z +2.80) |
| isozyme pairs (declared ceiling) | +0.186 (z +2.71) | +0.193 (z +3.36) |

**0th percentile of 30 pathways in both cohorts.** Not infiltrate: purity +
leukocyte fraction moves it -0.109 -> -0.091 (n = 1,007); deleting the
death-receptor module leaves 14 genes at -0.111 / -0.087; the most negative
module is the **intrinsic** pathway's cytosolic members (`APAF1`, `BMF`, `HRK`)
at -0.332 / -0.257.

### 4.4 Where the BCL2 family sits inside this - and why 4.1-4.3 cannot order it

**This is the hinge between section 4 and section 3, and it is easy to miss.**

11 of the 12 priming genes are in MitoCarta and all 11 are MOM. So the
compartment rule that explains the wider 44 **cannot separate the priming
genes** - they are all on the same side of it. Yet they span +0.50 to -0.27.

So the priming ordering is a *different* phenomenon from the compartment split,
sitting inside one of its halves, and E13 (3.4) is what characterises it: real,
within-compartment, and **a general property of the OXPHOS axis rather than a
property of the death machinery**.

---

## 5. A MODEL of apoptotic priming regulation by MYC and OXPHOS

Stated as an architecture with separately falsifiable arrows, so the mouse model
can accept, reject or extend each one. **Nothing here is established; it is the
most economical structure consistent with every observation in sections 3 and 4,
including the negatives.**

### 5.1 The architecture

```
                            MYC ACTIVITY
                     (FELSHER__MITOSTRIP, 61 genes)
                                 |
        (a) rho = 0.32 / 0.26    |    (b) SD 0.149 / 0.135, UNSORTED
           +---------------------+----------------------+
           |                                            |
           v                                            |
  MITOCHONDRIAL BIOGENESIS PROGRAMME                     |
  (read out as OXPHOS subunits GSVA, 89 genes)           |
           |                                            |
     +-----+--------------------+                       |
     |                          |                       |
(c)  v                     (d)  v                  (b)  v
REGULON TRANSCRIPTS        CYTOSOLIC APOPTOTIC    a real, comparably
20 of the 44               EFFECTORS             large association with
11 of the 12 priming       24 of the 44          the same transcripts,
genes                      1 of the 12 (BMF)     by a route this study
     UP                    APAF1 BMF HRK,        has NOT identified
median +0.25 / +0.20       death-receptor /
GENERIC - any regulon      NF-kB arm
member does this                DOWN
     |                     median -0.16 / -0.15
     |                     SPECIFIC - no other
     |                     programme does this
     |                          |
     +-----------+--------------+
                 |
      (e) WITHIN the regulon half, the OXPHOS axis resolves gene
          by gene - a GENERAL property of the axis (E13), not of
          the death machinery
                 |
                 v
   BCL2-FAMILY TRANSCRIPT CONFIGURATION
   up:   BAD +0.50  BCL2L1 +0.39  BIK +0.28  BBC3 +0.27  BID +0.21
   down: MCL1 -0.27  BMF -0.21  BCL2L11 -0.17  BCL2 -0.12  PMAIP1 -0.11
   flat: BCL2A1 -0.06  BCL2L2 -0.08
   NOT a balance: pro and anti interleave; the ratio grid is 92-95% additive
                 |
                 v  (f) THE STEP THIS STUDY CANNOT TAKE
        APOPTOTIC PRIMING (protein, BH3 profiling)
```

### 5.2 Each arrow, its evidence, and its status

| arrow | statement | human evidence | status |
|---|---|---|---|
| **(a)** | MYC activity and mitochondrial biogenesis co-vary | `rho = 0.320 / 0.256` after proliferation adjustment | **correlation only.** Direction not identifiable from cross-sectional bulk RNA (N6) |
| **(b)** | MYC has its own association with the apoptotic transcripts, **off the mitochondrial axis** | own component SD **0.149 / 0.135** > borrowed **0.100 / 0.063**, both >> 0.03 sampling noise; gap across the compartment split **-0.019 / +0.004** | **real and unexplained.** The largest hole in the human model |
| **(c)** | The biogenesis programme raises regulon transcripts, including 11 of the 12 priming genes | median +0.249 / +0.197 | **real but GENERIC.** At its matched null (z **+0.48 / +0.62**), 63rd percentile of 30 pathways, below mitophagy's mitochondrial half (+0.210 / +0.184) |
| **(d)** | ...and lowers cytosolic apoptotic effectors | mean **-0.106 / -0.094** against its own null (median -0.164 / -0.145), **0th percentile of 30 pathways**, unique among four curated membrane-spanning programmes, survives purity + leukocyte and death-receptor deletion | **the apoptosis-specific claim** |
| **(e)** | Inside the regulon half the OXPHOS axis resolves gene by gene, and this is a property of the axis | MOM/MOM ratios not flat (0.246 / 0.147) but random MOM/MOM pairs equally not flat (0.200 / 0.153), excess z **+1.18 / -0.16 - does not replicate** | **the ordering is real; its apoptotic specificity is NOT shown** |
| **(f)** | The configuration reaches protein and sets priming | **none.** Transcript abundance only | **NOT MEASURED. The mouse arm's step** |

### 5.3 What the model says about priming - five propositions

**M1. OXPHOS status sets the LEVEL of individual BCL2-family transcripts, gene
by gene, and does not set a BALANCE.**
Evidence: pro and anti interleave across the whole range (3.1); the ratio grid
is 92-95% additive with no interaction term (3.2); **0 of 35** ratios beat both
their components in both cohorts on OXPHOS (3.3).

**M2. The configuration it sets is specific, directional and replicates.**
In OXPHOS-high tumours: `BAD`, `BCL2L1`, `BIK`, `BBC3`, `BID` higher; `MCL1`,
`BMF`, `BCL2L11`, `BCL2`, `PMAIP1` lower. **The axis of the configuration is
`BCL2L1` up with `MCL1` down** - the two anti-apoptotic genes that move most,
and in opposite directions. This is the one thing here a functional assay can be
pointed at.

**M3. MYC does not participate in the configuration.**
Its range over the 12 is narrower in both cohorts; its ordering of the wider 44
goes to **-0.04 / -0.06** when OXPHOS is removed while OXPHOS's *rises* when MYC
is removed (4.1). The two genes MYC holds on to move in **opposite** directions
(3.1 reading 3), so it is not a coherent balance shift. But note **(b)**: MYC's
own association is real and simply is not organised along this axis.

**M4. There is no `MYC x OXPHOS` interaction on the transcript configuration.**
Tested four ways; three falsifiers broke it; the apparent interaction is
curvature in the main effects (3.6). Independently, the sibling **pre-registered**
study found the *functional* interaction on BH3 priming null (N2).

**M5. The MYC signal that appears in ratios is a ratio artefact.**
A ratio cancels the shared, OXPHOS-borrowed component and enriches the
non-shared, MYC-own one (3.5). **Any ratio-based priming analysis in any species
will show this**, and it is not evidence of a priming-specific mechanism.

### 5.4 What the model forbids

If any of these appears in a merged human-mouse write-up, the human arm is being
misquoted.

| forbidden | why |
|---|---|
| "OXPHOS-high tumours are more primed" | priming is post-translational; N3 |
| any `MYC x OXPHOS` interaction language | N1 exploratory, N2 pre-registered |
| "the pro-apoptotic genes rise with OXPHOS" | they do not move as a block; `BMF` -0.21 against `BAD` +0.50 |
| "MYC represses BCL2, shifting the balance to death" | estimator-dependent; N5 |
| reading (c) as apoptosis-specific | it is what any regulon member does |
| "the priming ratio increases with OXPHOS" | name the transcript that moves |

### 5.5 Verifiable predictions

Ordered by how much a mouse experiment could add.

**V1 - THE DEPENDENCY PREDICTION. The model's sharpest output.**
If the configuration in M2 reaches protein, **OXPHOS-high tumours should be more
BCL-xL-dependent and less MCL1-dependent than OXPHOS-low tumours.** Testable by
BH3 profiling with selective peptides (HRK for BCL-xL, MS1 for MCL1) or by
sensitivity to selective mimetics (A-1331852 vs S63845). **This is the one
prediction that converts a transcript pattern into a therapeutic statement**, and
the human arm cannot test it.

**V2 - THE NULL THAT MUST HOLD.**
BH3 profiling stratified by OXPHOS should show **no `MYC x OXPHOS` interaction**
on overall priming. Two independent human results point the same way (N1
exploratory transcript, N2 pre-registered functional). **If the mouse shows a
robust interaction, that is a real disagreement** and either the human arm is
wrong or the systems differ in a way worth naming.

**V3 - THE COMPARTMENT PREDICTION.**
In any system where an OXPHOS axis can be defined: regulon members up, cytosolic
effectors down. **The down leg is the novel one** - `Apaf1`, `Bmf`, `Hrk` and the
death-receptor / NF-kB arm below the axis while regulon members sit above.

**V4 - THE PERTURBATION PREDICTION, which is what would make "sets" causal.**
Move OXPHOS and read the 12 transcripts. **Lowering** it - ETC inhibition, or a
genetic biogenesis perturbation - should take `BCL2L1`, `BAD`, `BIK` and `BBC3`
down with it and bring `MCL1` up; raising it should do the reverse. **If they do not move, the human correlation is a property of
tumour heterogeneity (subtype, stroma, infiltrate) and not a regulatory
relationship**, and M1/M2 both fall.

**V5 - AN ALREADY-REJECTED REFINEMENT, stated so it is not re-proposed.**
The configuration is **not** organised by BH3 binding specificity: cognate
sensitiser/guardian pairs cancel no more than non-cognate ones (3.8). A mouse
model that predicts cognate pairing should expect the human data to disagree.

**V6 - THE SUBTYPE PREDICTION.**
The three priming ratios that survive every filter do so **only in Basal**
(`BBC3/BCL2`, `BID/BCL2`, `BBC3/MCL1`). **If the mouse model is basal-like, that
is where the comparison should be made; if it is luminal, the human arm predicts
these ratios will be flat.**

**V7 - THE ARROW (b) PREDICTION, which is an invitation rather than a test.**
MYC has a real, larger-than-borrowed association with these transcripts that is
not organised by the mitochondrial axis. **Anything the mouse arm can say about
what MYC does to the apoptotic machinery off the OXPHOS axis is the single
highest-value addition to this model.**

### 5.6 Where the model is weakest, ranked

1. **Arrow (b) is a named unknown.** MYC's own component has SD 0.149 / 0.135 and
   this study has not asked what orders it.
2. **Arrow (e) is not shown to be apoptosis-specific.** The within-compartment
   ordering of the BCL2 family is real; E13's excess over random MOM/MOM pairs
   does not replicate.
3. **Arrow (f) is unmeasured**, so every word between "transcript configuration"
   and "priming" is inference.
4. **V6 rests on the least powerful stratum.** Basal is 171 / 317 tumours.
5. **SCAN-B has no purity estimate**, so the infiltrate control for arrow (d) is
   TCGA-only.
6. **One comparator set is author-curated in-script** (the isozyme pairs) and is
   a declared ceiling, not a pinned catalogue.

---

## 6. The standing negatives

**These must survive into the mouse step intact.** The commonest way a two-arm
study goes wrong is a mouse result being read as confirming something the human
arm never claimed.

| # | the negative | strength |
|---|---|---|
| **N1** | **The `MYC x OXPHOS` interaction on the priming ratios is not supported.** Tested four ways, failed three (3.6). The apparent interaction is curvature in the main effects | tested, three falsifiers |
| **N2** | **The sibling PRE-REGISTERED study found the functional `MYC x OXPHOS` interaction on apoptotic priming null.** Different endpoint from N1 - BH3 priming there, transcript ratios here - so neither confirms the other. **Nothing in this exploratory arm overturns that null and nothing may be written as if it did** | pre-registered, frozen at `d3ac60e` |
| **N3** | **Priming is not measurable in transcript abundance.** It is post-translational and protein-interaction-level. Write "carries a higher `BAD`/`MCL1` transcript ratio", never "is more primed" | definitional |
| **N4** | **The priming ratios carry no pair-specific information.** 92.1-94.6% additive; at most 5-8% pair-specific; 0 of 35 beat both components on OXPHOS | measured (3.2, 3.3) |
| **N5** | **"MYC represses BCL2" is estimator-dependent.** -0.22 on the reference estimator, -0.06 and -0.01 on the other two after conditioning. Do not write it | measured |
| **N6** | **Mediation versus confounding is not identifiable.** `MYC -> OXPHOS -> genes` and a common cause give identical partial correlations. What IS ruled out is OXPHOS acting *through* MYC | structural |
| **N7** | **MYC mRNA is not MYC activity.** `rho(log2 MYC, OXPHOS subunits) = -0.032` against `+0.388` for the activity signature | measured |
| **N8** | **The within-compartment ordering is a property of the OXPHOS axis, not of the death machinery.** Random MOM/MOM pairs are as non-flat as the BCL2-family ones; the excess does not replicate | measured (3.4) |

**What DID survive on the priming side** is narrow and should be carried as two
genes, not as a programme: **`BID`** (+0.12 to +0.28) and **`PMAIP1`/NOXA**
(-0.05 to -0.21) keep a MYC association after conditioning on OXPHOS, in both
cohorts under all three estimators. Neither is a member of any estimator, so it
is not self-overlap. Both are pro-apoptotic BH3-only proteins **moving in
opposite directions**.

---

## 7. What goes in the paper

### 7.1 The machinery-wide sentences

> In two independent human breast cancer cohorts (TCGA-BRCA n=1,095, SCAN-B
> n=3,207), transcript levels of the canonical apoptotic machinery are ordered
> along tumour OXPHOS status and not along MYC activity. The two axes are
> themselves correlated, but conditioning OXPHOS on MYC leaves the ordering
> intact while conditioning MYC on OXPHOS abolishes it, and neither is explained
> by proliferation, tumour purity or immune infiltrate. What predicts a gene's
> position is membership of the nuclear-encoded mitochondrial regulon - regulon
> members rise with OXPHOS, cytosolic ones fall - and not whether it promotes or
> prevents death. The regulon members are ordered no more steeply than those of
> any comparably composed mitochondrial programme, but the cytosolic members are
> the only ones among four curated membrane-spanning programmes to run against
> OXPHOS, which localises what is specific to apoptosis on the cytosolic side.

### 7.2 The priming subsection - descriptive, four sentences

1. *The family does not move as a block and functional class does not predict
   position.* The twelve transcripts span OXPHOS from `BAD` **+0.503** to `MCL1`
   **-0.266**; the two most positive are `BAD` (pro) and `BCL2L1` (anti,
   **+0.388**), the two most negative `MCL1` (anti) and `BMF` (pro, **-0.212**).
   Narrower on MYC, `BID` **+0.300** to `PMAIP1` **-0.209**.
2. *The ratio matrix is additive.* **92.1-94.6%** of all 35 values from numerator
   and denominator identity alone, numerator carrying about twice the
   denominator. **Present ratios as a compact display of component correlations,
   not as a measurement of priming.** Keep this sentence if only one survives.
3. *The OXPHOS association is not mitochondrial mass* - MOM/MOM ratios are not
   flat - *but neither is it apoptosis-specific*: random MOM/MOM pairs are
   equally non-flat and the excess does not replicate.
4. *The two axes differ in the same direction as the wider machinery*, with
   `BID` and `PMAIP1` as the two named exceptions moving oppositely.

### 7.3 Figures

| slot | figure | note |
|---|---|---|
| **Main A** | `E11_fig9` A - the 44 on the MYC x OXPHOS plane, square identically-scaled axes | the cloud is taller than wide and the colour separates vertically |
| **Main B** | `E11_fig9` B - the conditioning ladder with the matched-null band and bootstrap intervals | **not optional.** It turns A from a description into a test |
| **Main / supp** | `E14_fig6` - the one-panel specificity figure | four programmes, only apoptosis below zero |
| **Priming, main or supp** | `E10_fig5` - the 12 transcripts with 95% intervals, by side | this is the priming result. Everything else in the priming section is derived from it |
| Priming supp | `E10_fig4` - does a ratio beat its stronger component | the falsifier, drawn |
| Priming supp | `E10_fig3` / `E10_fig6` - the ratio heatmaps with `*` `**` `^` and border marks | **not one MYC cell is marked anywhere**; the three double-cohort cells are all Basal OXPHOS |
| Priming supp | `E10_fig8` - the two axes collapsed to one number, three ways | shows why it must be `\|rho(OXPHOS)\| - \|rho(MYC)\|` and not the signed difference; the two disagree in sign on 22 of 70 cells |
| Priming supp | `E13` - the MOM/MOM control against its random-pair null | N8 |
| Supplementary | `E11_fig2` mitoribosome control; `E14_fig5` infiltrate falsifier | **each must be cited in its principal figure's legend** |
| Supplementary | `E15_fig1` per-gene dumbbell; `E15_fig4` two-bar summary (**carry its SD-not-SE sentence**); `E15_fig5` the difference alone | for a reader who wants to know which genes carry it |

Delete on-figure titles for submission. Tables as listed in 2.8.

### 7.4 Statistics

The permutation nulls **are** the test - 2,000 expression-matched draws, and for
the split a draw matched on sub-mitochondrial compartment as well; 500 draws for
the MOM/MOM pair null. Report foreground z and percentile; **do not add p-values
on top.** Intervals come from **1,000 tumour-level bootstrap resamples** (tumours,
not genes - the 44 are co-expressed and a gene-level bootstrap would return an
interval far too narrow).

**Two things not to do.** No per-gene p-values or FDR across the 44 or the 35 -
that is the grid-of-cells trap and it invites gene-picking. And do not report the
composition null as "not significant" and move on; report it as the bound it is.

### 7.5 Language rules

| do not write | write instead |
|---|---|
| "the priming ratio increases with OXPHOS" | name the transcript that moves |
| "OXPHOS-high tumours are more primed" | "carry a higher `BAD`/`MCL1` transcript ratio" |
| "MYC represses BCL2, shifting the balance to death" | nothing - N5 |
| any interaction language | nothing - N1, N2 |
| "drives", "engages", "activates" | "is ordered along", "tracks", "co-varies with" |
| "mitochondrial genes" (as localisation) | "transcripts of the nuclear-encoded mitochondrial regulon" |

---

## 8. THE HANDOVER

### 8.1 What the human arm cannot decide, by construction

| question | why not | what would settle it |
|---|---|---|
| Does OXPHOS **cause** the configuration? | cross-sectional bulk RNA; mediation and confounding give identical partial correlations (N6) | **V4** - a perturbation that moves OXPHOS and reads the 12 |
| Does the configuration reach **protein**? | transcript abundance only | protein-level readout of the 12 |
| Does it set **priming**? | priming is post-translational (N3) | **V1** - BH3 profiling stratified by OXPHOS |
| Does `MYC x OXPHOS` set a threshold? | transcript version failed three falsifiers (N1); functional version pre-registered null (N2) | **V2** - BH3 profiling under a perturbation that moves OXPHOS |
| What is arrow **(b)**? | not asked here | **V7** - anything at all |
| Is this breast-specific? | two human breast cohorts only | a different tissue, or normal tissue |

### 8.2 The claims most worth taking to the mouse

**H1 - the dissociation.** *The apoptotic machinery is ordered along OXPHOS and
not along MYC, and conditioning is what separates them.* The only claim with a
pre-written falsifier; replicates across cohorts (per-gene values correlate
**0.87-0.92** between them under every adjustment).

**H2 - the cytosolic anomaly.** *Cytosolic apoptotic transcripts run against
OXPHOS while regulon members run with it, and no other membrane-spanning
programme does this.* The novel claim, and the one with the cleanest prediction
(**V3**).

**H3 - the priming configuration.** *OXPHOS sets BCL2-family transcript levels
gene by gene, not as a balance; the axis of the configuration is `BCL2L1` up
with `MCL1` down.* This is what section 5 adds, and **V1** is how it is tested.

### 8.3 What counts as agreement - and what does not

**Counts as agreement**

- The compartment split reproduces **in sign and direction**: regulon members up
  along an OXPHOS axis, cytosolic members down.
- The ordering survives conditioning on a mouse MYC estimator, and the mouse MYC
  ordering does not survive conditioning on OXPHOS. **The conditioning asymmetry
  is the claim, not the marginal correlations.**
- The M2 configuration reproduces: sensitisers and *Bcl2l1* up, *Mcl1* down.
- A curated non-apoptotic membrane-spanning comparator, built species-natively,
  fails to reproduce the cytosolic negative.

**Does NOT count as agreement, however tempting**

- *"MYC correlates with OXPHOS in mouse too."* The human point is that the
  machinery's ordering is OXPHOS's, not MYC's. A bare axis correlation confirms
  nothing.
- *"Apoptotic genes correlate with OXPHOS in mouse."* Arrow (c) says the regulon
  half does that by composition. Without a matched null and a named comparator
  this is not evidence.
- *"The `BAD`/`MCL1` ratio tracks OXPHOS in mouse."* So does `BAD` alone, better.
  A ratio result needs its `gain` (3.3) before it means anything.
- *"`MYC x OXPHOS` is significant in mouse."* N1 and N2. It needs F1, F2 and F3
  before it means anything.
- Any priming claim from RNA. N3.

**Would count as disagreement, and would be important**

- Mouse shows the machinery ordering is **MYC's** after conditioning. Rung 3
  would need restating and the human result would look like a property of human
  tumour heterogeneity.
- Mouse shows the cytosolic members track OXPHOS **positively**. H2 becomes a
  human-tumour phenomenon.
- **Mouse shows a robust `MYC x OXPHOS` interaction on functional priming.** This
  is the highest-stakes disagreement available, because it contradicts a
  pre-registered null as well as an exploratory one.
- Mouse shows the ratio *balance* moves while its components do not. That would
  say the additive ceiling (N4) is a property of bulk human tumours rather than
  of the biology.

### 8.4 Traps specific to the cross-species step

1. **No ortholog function is called in this repo, in either direction.** The
   reconciliation is a comparison of **conclusions**, not of gene lists. Any
   gene-level mapping happens elsewhere. The tripwire here stays clean.
2. **Human and mouse MitoCarta are different inventories.** Regulon membership
   is not automatically the same object across species; rebuild it
   species-natively and expect the count in each half to differ.
3. **Never pool GSVA or mitoPPS values across cohorts** - and a species is a
   cohort. Compare correlations, patterns and rankings.
4. **The composition null must be rebuilt species-natively too**, or the
   comparison is circular.
5. **The human confounds are human.** Purity and leukocyte fraction mattered
   here; a mouse model has different ones, and the absence of these two is not
   the absence of confounding.
6. **A mouse `Bcl2l2` / `Bcl2a1` mapping is not one-to-one.** Mouse has several
   *Bcl2a1* paralogues. The human grid is 7 x 5; a mouse grid will not be.
7. **The sibling study is closed.** `myc_human_validation` at `d3ac60e` is
   pre-registered, found nothing supported, does not go in the paper, and is not
   reopened. **N2 outranks anything exploratory that appears to contradict it.**
8. **Different measurement types are not a disagreement.** The human arm reads
   bulk tumour RNA. If the mouse arm reads protein, function or perturbation, a
   null there does not refute a correlation here and a positive there does not
   confirm one. Say which claim each addresses before comparing.

---

## 9. Superseded statements - read this, not that

| source | superseded statement | read instead |
|---|---|---|
| this file, 2026-09-03 morning version | organised around the compartment split, with priming as section 5.3 | sections 3 and 5 here. The compartment result is unchanged and is section 4 |
| `paper_opening_human.md`, rung 4 | "the ordering is no steeper than any matched gene set ... the bulk transcriptome reports mitochondrial content" | 4.3. True of the regulon half, false of the cytosolic half |
| `paper_opening_human.md`, open item 1 | "the ~1.3 SD residue ... would be settled by a curated non-apoptotic comparator" | **answered.** `e14_curated_comparators.md` |
| `handoff_evening.md`, the finding paragraph | "...and it is not specific to apoptosis" | section 1 and 4.3 |
| `e10_machinery_and_priming.md` R1 | the 0.453 split as a localisation-organised death programme | 4.2. It is regulon membership |
| `e10_machinery_and_priming.md` R1 | "protein localisation" language | co-regulation. `priming_interaction_tested.md` C1 |
| `e10_machinery_and_priming.md` R4 | the Pearson spot check | withdrawn; superseded by E09 |
| `e10_machinery_and_priming.md` R5 | "6 ratios gain on OXPHOS, 3 on MYC" | **adjusted: 0 on OXPHOS, 5 on MYC** (3.3) |
| any note before 2026-09-02 evening | "outer mitochondrial membrane", "acts at the mitochondrion" | "nuclear-encoded mitochondrial regulon" |

---

## 10. Still open on the human side

1. **A cytosolic stress programme that is not apoptotic**, as the next
   comparator for arrow (d). If a proteotoxic or integrated-stress-response
   module also runs against OXPHOS, then H2's property is "cytosolic and
   stress-responsive" rather than "cytosolic and apoptotic". **The single most
   valuable next human analysis**, and it should be done before H2 is leaned on
   hard in either arm.

   - **Its other half: what arrow (b) actually is.** `data/collectri_human/` is
     pinned and has never been used for this. Two forms - which CollecTRI
     regulons contain the 44, and whether per-sample regulon activity tracks
     either half or MYC's own component. **Falsifier for the "two regulons"
     language:** if no regulon separates the halves beyond the mitochondrial
     ones, the phrase goes and "different correlate" replaces it.
     `e15_two_axis_gene_view.md` V5.

2. **Fold the cognate-pairing test (3.8) into `E10`** so it regenerates from
   code. It is currently the only number in this document without a script.
3. **`LumA` alone** - the homogeneous stratum where stromal and immune
   composition vary least. Named in E14 as the falsifier for H2, and it is also
   where V6 should be checked from the other side.
4. **The `BID` / `BBC3` numerator test named in 3.3b is the cheapest real
   falsifier left in the priming section**, and it needs a third cohort rather
   than more work on these two. METABRIC is out of scope, so this is a note for
   whoever reopens that decision.
5. Score `COLLECTRI_MYC_STIM` (739 genes, in the snapshot, never scored) as a
   fifth base; drop `ELLWOOD` and recompute the entanglement slope. Both need a
   pipeline re-run.
6. The ER-negative fatty-acid-oxidation reversal, untouched.
7. A stranded total-RNA dataset for `MT-ND6` and `CO1`/`CO2`.

**Out of scope until reopened:** MCbiclust / forkscale, survival, treatment,
METABRIC, DepMap, causal modelling on human data.
