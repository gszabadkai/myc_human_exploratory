---
date: 2026-09-03
purpose: ONE document. What the human arm concludes, what goes in the paper,
         and what it hands to the mouse arm to reconcile against.
supersedes: nothing - it summarises and points. Where a source note and this
         document disagree, THIS ONE is current and section 8 says why.
sources: 2026-09-02_e10_machinery_and_priming.md,
         2026-09-02_e11_prolif_adjusted.md,
         2026-09-02_e14_curated_comparators.md,
         2026-09-02_paper_opening_human.md,
         2026-09-02_priming_interaction_tested.md,
         2026-09-02_handoff_evening.md
posture: EXPLORATORY. Nothing in the human arm is pre-registered.
---

# The human arm, stated for reconciliation with the mouse

**I have not read the mouse repo in preparing this.** It is not attached and
CLAUDE.md keeps it read-only. This document is one side of the comparison,
written so the other side can be laid against it without either arm having to
re-derive what the first one meant.

---

## 0. Status of every number below

| | |
|---|---|
| **E09 - E14** | all run by the author. Every number below is read from a saved object, not from a dry run |
| **E10** | re-sourced 2026-09-03 00:43 and now carries `covariate: PROLIF_DISJOINT`. Its verification check passed exactly as predicted: ratios beating both components **0 of 35 on OXPHOS, 5 of 35 on MYC** |
| **E14** | on disk is the current build - `plot_sets` of four, `mito_pair_md_agree = FALSE`. Both verdicts reproduce: `APOPTOSIS-SPECIFIC in both cohorts`, `SURVIVES` on the infiltrate falsifier |
| **E15** | written 2026-09-03 and **not yet run by the author**. A display script: it recomputes nothing, and asserts its numbers are bit-equal to E14's. See `e15_two_axis_gene_view.md` |
| **one small re-run outstanding** | `E10` gained an `additive_fit` block on 2026-09-03 so that section 5.3's headline number regenerates from code instead of from a note. **~15 s.** The values do not change - the block reproduces 0.921-0.946 exactly against the object already on disk |

Two cohorts throughout: TCGA-BRCA n=1,095, SCAN-B/GSE202203 n=3,207. Two arms of
evidence only ever count when they agree; a single cohort is a description.

---

## 1. The finding, in one paragraph

> In two independent human breast cancer cohorts, transcripts of the canonical
> apoptotic machinery are ordered along tumour OXPHOS status and not along MYC
> activity. The two axes are themselves correlated (0.26-0.32 after
> proliferation) and rank the 44 genes almost identically (0.61-0.73), so
> comparing them separately could not settle which carries the ordering -
> conditioning does. Removing MYC leaves the OXPHOS ordering intact; removing
> OXPHOS abolishes MYC's. Neither is proliferation, tumour purity or immune
> infiltrate. What predicts a gene's position is whether it belongs to the
> nuclear-encoded mitochondrial regulon, not whether it promotes or prevents
> death. **The mitochondrial members behave as any mitochondrial programme's
> members do. What no other programme does is send its cytosolic members the
> other way.**

That last sentence is the 2026-09-02 evening change and it is section 3.2.

---

## 2. The claim ladder

Five rungs. The candidate claim was *"the apoptotic machinery is driven, or at
least more strongly correlated, by OXPHOS than directly by MYC."*

| # | statement | status | what licenses it |
|---|---|---|---|
| 1 | The 44 genes correlate more strongly with OXPHOS than with MYC activity | **SUPPORTED** | SD of the 44 per-gene rho 0.313 vs 0.191 (TCGA), 0.245 vs 0.157 (SCAN-B). Tumour-level bootstrap SD ratio **1.64 [1.47, 1.78]** and **1.56 [1.46, 1.65]**. Holds inside Luminal and Basal separately (1.3x to 1.7x) |
| 2 | ...and it is not proliferation | **SUPPORTED** | Unchanged by partialling `PROLIF_DISJOINT` (318 genes, disjoint from the MYC estimator, the OXPHOS arm and the mitoribosome), by the `__PROLIFSTRIP` / `__BOTHSTRIP` estimators, and in TCGA by purity + leukocyte fraction on top (n = 1,007) |
| 3 | ...and it is OXPHOS **rather than** MYC | **SUPPORTED - the strong result** | Conditioning. OXPHOS split 0.453 -> **0.485** (TCGA) and 0.489 -> **0.525** (SCAN-B) with MYC removed; MYC split 0.187 -> **-0.043** and 0.137 -> **-0.058** with OXPHOS removed. MYC's ordering was inherited, OXPHOS's was not |
| 4 | ...and it is specific to apoptosis | **SPLIT IN TWO - see 3.1 and 3.2** | 4a NOT supported for the mitochondrial half; 4b SUPPORTED for the cytosolic half |
| 5 | "driven" | **NOT TESTED, and untestable here** | Cross-sectional bulk tumour RNA. Every causal word must come from the mouse arm |

**Rung 3 is the one to carry into the mouse comparison.** It is the only rung
that required a test rather than a description, and it is the only one whose
falsifier was written before the answer was seen.

---

## 3. What the ordering actually is

Not pro-death versus pro-survival. **Membership of the nuclear-encoded
mitochondrial regulon.** MitoCarta membership predicts a gene's position on the
OXPHOS axis at 0.45-0.49; its annotated direction of effect at only 0.22-0.28.
MitoCarta is independent of the death curation and agrees with it on 43 of 44
genes.

**Read "mitochondrial" as co-regulation, not localisation.** MitoCarta is a
proteome catalogue and half the BCL2 family translocates - BAX is cytosolic
until activated, BID must be cleaved, BAD is held by 14-3-3. A transcript has no
idea where its protein ends up. This reframe matters for the mouse step: the
claim this data type can support is **transcriptional co-regulation**, not
protein localisation.

`APAF1` is the clean illustration and it cuts the right way: intrinsic-pathway
by Reactome, absent from MitoCarta because the apoptosome is cytosolic, and the
most negative of the 44 (-0.452). Pathway membership gets it wrong, regulon
membership gets it right.

### 3.1 The mitochondrial half is NOT special (rung 4a)

| test | result |
|---|---|
| against its own expression- and sub-compartment-matched null | +0.167 / +0.143 against +0.139 / +0.117, **z +0.48 / +0.60** |
| among 30 MitoCarta leaf pathways | **63rd percentile in both cohorts** |
| against mitophagy's mitochondrial half | mitophagy is **higher**: +0.210 / +0.184 |

Everything E11 measured about the mitochondrial 20 stands, and stands as a
negative. An OXPHOS score correlates with a large fraction of the
transcriptome; 20 mitochondrial genes of this expression profile reach that
value by construction.

### 3.2 The cytosolic half IS specific (rung 4b) - the 2026-09-02 result

The residue E11 left was closed by replacing the random null with **named
non-apoptotic programmes that span the outer membrane**. Protein import cannot
serve - 62 of its 62 genes are in MitoCarta, so it has no cytosolic half at all,
and the same is true of cristae formation (31/31) and calcium transport (21/22).
Mitophagy can: 16 outer-membrane genes against 22 cytosolic effectors, zero
genes shared with the 44.

Split on OXPHOS, proliferation-adjusted:

| programme | TCGA | SCAN-B |
|---|---|---|
| **apoptotic machinery (44)** | **0.453** | **0.489** |
| Fe-S cluster assembly | 0.397 | -0.013 |
| mitophagy | 0.112 | 0.219 |
| isozyme pairs (declared ceiling) | -0.092 | -0.022 |
| *random matched null* | *0.276 +/- 0.131* | *0.319 +/- 0.126* |

**And the reason is entirely the cytosolic half.** Mean rho with OXPHOS,
observed against its own null:

| programme | cytosolic half, TCGA | SCAN-B |
|---|---|---|
| **apoptotic machinery** | **-0.106** (null -0.036, z -1.2) | **-0.094** (z -1.3) |
| mitophagy | +0.143 (z +2.4) | +0.096 (z +1.8) |
| Fe-S cluster assembly | +0.118 (z +1.7) | +0.159 (z +2.9) |
| isozyme pairs | +0.186 (z +2.7) | +0.193 (z +3.4) |

On the MitoCarta pathway ladder the machinery's cytosolic half is at the **0th
percentile of 30 pathways in both cohorts**; its mitochondrial half is at the
63rd.

**It is not immune infiltrate**, and that had to be tested because the 24
cytosolic genes are `FAS`, `FASLG`, `TNF`, `TNFRSF1A/10A/10B`, `TRADD`, `FADD`,
`CASP10`, `CFLAR`, `BIRC2/3`, `XIAP`, `NFKB1`, `RELA` - inflammatory as much as
apoptotic, and trap 2 puts rho(OXPHOS, leukocyte fraction) at -0.158:

- purity + leukocyte fraction on top of proliferation: -0.109 -> **-0.091**
- delete the death-receptor module entirely: 14 genes at **-0.111 / -0.087**,
  9 of 14 negative in both cohorts
- most negative module is the **intrinsic** pathway's cytosolic members
  (`APAF1`, `BMF`, `HRK`) at **-0.332 / -0.257**

**Two things that did not go the way they were predicted to, recorded because
they matter more than the ones that did.** (i) The isozyme-pair ceiling - same
reaction either side of the membrane - was predicted to give a large split and
gave **none** (-0.09 / -0.02, both halves positive). Being in the mitochondrial
regulon does not by itself make a split when the cytosolic counterpart does the
same job. (ii) **MYC orders mitophagy and does not order apoptosis**: mitophagy
splits at 0.233 / 0.306 on the MYC axis, above its null, while the machinery
sits at 0.187 / 0.137, at or below its own.

---

## 4. The standing negatives

**These must survive into the mouse step intact.** The commonest way a
two-arm study goes wrong is a mouse result being read as confirming something
the human arm never claimed.

| # | the negative | strength |
|---|---|---|
| N1 | **The MYC x OXPHOS interaction on the priming ratios is not supported.** Tested four ways and failed three: splines (linear vs spline coefficients correlate **0.125**), a low-entanglement estimator (cross-cohort replication **reverses**, r = -0.33), OXPHOS tertiles (non-monotone, TCGA +0.021 vs SCAN-B -0.008). The apparent interaction is curvature in the main effects | tested, three falsifiers |
| N2 | **The sibling PRE-REGISTERED study found the functional MYC x OXPHOS interaction on apoptotic priming null.** Different endpoint from N1 - BH3 priming there, transcript ratios here - so neither confirms the other. **Nothing in the human exploratory arm overturns that null and nothing may be written as if it did** | pre-registered, frozen at `d3ac60e` |
| N3 | **Priming is not measurable in transcript abundance.** It is post-translational and protein-interaction-level. Write "carries a higher `BAD`/`MCL1` transcript ratio", never "is more primed" | definitional |
| N4 | **The priming ratios carry no pair-specific information.** An additive model of numerator and denominator identity - no interaction term - explains **92.1 to 94.6%** of all 35 ratios in both cohorts on both axes, the numerator carrying about twice the denominator (R2 0.58-0.74 against 0.20-0.35). So at most 5-8% of a ratio can be pair-specific, which is the ceiling on anything "priming" could mean here. Adjusted, **0 of 35** beat both components on OXPHOS and 5 of 35 on MYC | measured |
| N5 | **"MYC represses BCL2" is estimator-dependent.** -0.22 on the reference estimator, -0.06 and -0.01 on the other two after conditioning. Do not write it | measured |
| N6 | **Mediation versus confounding is not identifiable.** MYC -> OXPHOS -> genes and a common cause give identical partial correlations. What IS ruled out is OXPHOS acting through MYC | structural |
| N7 | **MYC mRNA is not MYC activity.** rho(log2 MYC, OXPHOS subunits) = -0.032 against +0.388 for the activity signature. Anyone plotting MYC expression sees nothing | measured |

What DID survive on the priming side is narrow and should be carried as two
genes, not as a programme: **`BID`** (+0.12 to +0.28) and **`PMAIP1`/NOXA**
(-0.05 to -0.21) keep a MYC association after conditioning on OXPHOS, in both
cohorts under all three estimators. Both are pro-apoptotic BH3-only proteins
**moving in opposite directions**, so it is not a coherent shift in balance.

---

## 5. What goes in the paper

### 5.1 The four sentences

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

**The fourth sentence replaces the earlier one** ("...indicating that the bulk
transcriptome reports mitochondrial content rather than a selectively engaged
death programme"). That earlier version conceded too much: it was written when
the only comparator was a random draw.

### 5.2 Figures

| slot | figure | note |
|---|---|---|
| **Main, panel A** | `E11_fig9` A - the 44 on the MYC x OXPHOS plane, square identically-scaled axes, coloured by regulon membership, adjusted column, both cohorts | the cloud is taller than wide and the colour separates vertically |
| **Main, panel B** | `E11_fig9` B - the conditioning ladder, four rows, two cohorts, with the compartment-matched null band and bootstrap intervals | **not optional.** It is what turns A from a description into a test |
| **Main or supplementary** | `E14_fig6` - the one-panel specificity figure | x = what a programme's regulon members do, y = what its cytosolic members do. Four programmes, only apoptosis below zero |
| Supplementary | `E15_fig1` - the 44 as one bar per gene, MYC point and OXPHOS point, sorted by the difference, both cohorts | the per-gene version of panel A, for a reader who wants to know WHICH genes carry it. Cite `E15_fig3` with it if the difference is quoted in the text |
| Supplementary | `E11_fig2` - the mitoribosome control | **must be cited in panel B's legend.** It rules out "the adjustment emptied the MYC score" |
| Supplementary | `E14_fig5` - the infiltrate falsifier | **must be cited in `E14_fig6`'s legend.** The panel's whole weight rests on that negative not being infiltrate |
| Supplementary | `E14_fig1` head-to-head, `E14_fig3` pathway ladder | the three-step version of the same argument |
| Supplementary | `E13` panel A or `E10_fig5` - the 12 BCL2-family transcripts ranked | "pro- and anti-apoptotic members are interleaved across the whole range" |
| Supplementary | `E10_fig3` ratio matrix, `E10_fig6` luminal/basal split | descriptive only |

Delete on-figure titles for submission. Tables:
`E14_comparator_splits.csv`, `E14_mitocarta_pathway_ladder.csv`,
`E11_gene_rho_by_adjustment.csv`, `E15_gene_rho_two_axes.csv`.

### 5.3 The priming subsection - descriptive, three sentences

**Requires E10 to be re-sourced first (section 0).** In order:

1. *The family does not move as a block and functional class does not predict
   position.* Confirmed against the adjusted object: the twelve transcripts
   span OXPHOS from `BAD` **+0.503** to `MCL1` **-0.266**; the two most
   positive are `BAD` (pro) and `BCL2L1` (anti, +0.388), the two most negative
   `MCL1` (anti) and `BMF` (pro, -0.212). Narrower on MYC, `BID` **+0.300** to
   `PMAIP1` **-0.209**.
2. *The ratio matrix is additive.* **92.1-94.6%** of all 35 values from
   numerator and denominator identity alone, numerator carrying about twice
   the denominator. **Present ratios as a compact display of component
   correlations, not as a measurement of priming.** This is the sentence to
   keep if only one survives editing, and as of 2026-09-03 it regenerates from
   `E10` section 5.0b rather than from a note.
3. *The two axes differ in the same direction as the wider machinery*, with
   `BID` and `PMAIP1` as the two named exceptions moving oppositely.

**If the five MYC gainers are ever named, use the ADJUSTED five.** Adjustment
changed their identity completely, not just their count: they are now
`BBC3/BCL2` (+0.052 / +0.062), `BAD/BCL2` (+0.040 / +0.050), `BAD/BCL2L1`,
`BBC3/BCL2L2` and `BBC3/BCL2L1` - all under +0.07, and none of them is one of
the three the unadjusted run named. A list that reshuffles entirely under a
covariate is a description of noise, and the safest thing is not to name them
at all.

### 5.4 Statistics

The permutation nulls **are** the test - 2,000 expression-matched draws, and
for the split a draw matched on sub-mitochondrial compartment as well.
Foreground z and percentile; do not add p-values on top. Intervals come from
**1,000 tumour-level bootstrap resamples** (tumours, not genes - the 44 are
co-expressed and a gene-level bootstrap would return an interval far too
narrow). All five contrasts exclude their null in both cohorts.

**Two things not to do.** No per-gene p-values or FDR across the 44 - that is
the grid-of-cells trap and it invites gene-picking. And do not report the
composition null as "not significant" and move on; report it as the bound it is.

### 5.5 Language rules

| do not write | write instead |
|---|---|
| "the priming ratio increases with OXPHOS" | name the transcript that moves |
| "MYC represses BCL2, shifting the balance to death" | nothing - N5 |
| "OXPHOS-high tumours are more primed" | "carry a higher `BAD`/`MCL1` transcript ratio" |
| any interaction language | nothing - N1, N2 |
| "drives", "engages", "activates" | "is ordered along", "tracks", "co-varies with" |
| "mitochondrial genes" (as localisation) | "transcripts of the nuclear-encoded mitochondrial regulon" |

---

## 6. Known softness, stated so the mouse arm is not asked to confirm noise

1. **The comparator field is noisy.** `Fe-S cluster assembly` is 0.397 in TCGA
   and -0.013 in SCAN-B. The machinery is the only set that reproduces tightly
   (0.453 / 0.489).
2. **The mitophagy comparator is softer than the figures show.** Its PINK1/PRKN
   subset agrees on the anchor statistic `split` but **disagrees in sign on
   `med_diff` in TCGA** (+0.184 against -0.030). `split` is the criterion
   because it has been the anchor since E11, not because it is the one that
   passed - and I had seen both before writing the criterion down.
3. **Three comparators cannot make a p-value.** The deliverable is a ranking
   that reproduces across two cohorts.
4. **The sub-compartment ladder rests on 2 inner-membrane genes** in the 44.
   E14's pathway ladder is the wide version and it agrees (MIM 0.18-0.25 > MOM
   0.12-0.13 > matrix 0.10-0.12), but the 44-gene version alone is a direction,
   not a gradient.
5. **`LumA` alone has never been run**; `Basal` rests on 171 TCGA samples.
6. **One set is author-curated in-script** - the isozyme pairs - and is the only
   set here not from a pinned catalogue. It is a declared ceiling.
7. **`CYCS` is one of the 89 OXPHOS-arm genes**, so its OXPHOS value is partly
   self-correlation; marked on every panel.

---

## 7. THE HANDOVER - what the mouse arm is being asked

### 7.1 What the human arm cannot decide, by construction

| question | why not | what would settle it |
|---|---|---|
| Does OXPHOS **cause** the ordering? | cross-sectional bulk RNA; mediation and confounding give identical partial correlations (N6) | a perturbation that moves OXPHOS and reads the machinery |
| Is the cytosolic anti-correlation **functional**? | transcript abundance, no protein, no activity | protein-level or functional readout |
| Does MYC x OXPHOS set a priming threshold? | priming is post-translational (N3); the transcript version failed three falsifiers (N1) | BH3 profiling under a perturbation that moves OXPHOS |
| Is any of this specific to breast, or to tumours? | two human breast cohorts only | a different tissue, or normal tissue |

### 7.2 The two claims most worth taking to the mouse

**H1 - the dissociation.** *The apoptotic machinery is ordered along OXPHOS and
not along MYC, and conditioning is what separates them.* This is rung 3, the
only rung with a pre-written falsifier, and it replicates across two cohorts
(the 44 per-gene values correlate 0.88-0.92 between them under every
adjustment).

**H2 - the cytosolic anomaly.** *Cytosolic apoptotic transcripts run against
OXPHOS while regulon members run with it, and no other membrane-spanning
programme does this.* This is the novel claim, it is the one a reviewer will
push on, and it is the one with the cleanest prediction: in any system where an
OXPHOS axis can be defined, `Apaf1`, `Bmf`, `Hrk` and the death-receptor / NF-kB
arm should sit **below** it while the regulon members sit above.

### 7.3 What would count as agreement - and what would not

**Counts as agreement**

- The compartment split reproduces **in sign and in direction**: regulon members
  up along an OXPHOS axis, cytosolic members down.
- The ordering survives conditioning on a mouse MYC estimator, and the mouse MYC
  ordering does not survive conditioning on OXPHOS. **The conditioning
  asymmetry is the claim, not the marginal correlations.**
- A curated non-apoptotic membrane-spanning comparator, built species-natively,
  fails to reproduce the cytosolic negative.

**Does NOT count as agreement, however tempting**

- *"MYC correlates with OXPHOS in mouse too."* The human arm's point is that the
  machinery's ordering is OXPHOS's and not MYC's. A bare axis correlation
  confirms nothing on either side.
- *"Apoptotic genes correlate with OXPHOS in mouse."* Rung 4a says the regulon
  half does that by composition. Without a matched null and a named comparator
  this is not evidence.
- *"MYC x OXPHOS is significant in mouse."* N1 and N2. A significant interaction
  in one model is what the human arm already found and then broke with three
  falsifiers; the mouse version needs the same three before it means anything.
- Any priming claim from RNA. N3.

**Would count as disagreement, and would be important**

- Mouse shows the machinery ordering is **MYC's** after conditioning. That would
  say the human result is a property of human tumour heterogeneity - subtype,
  stroma, infiltrate - rather than of the regulon, and rung 3 would need
  restating.
- Mouse shows the cytosolic members track OXPHOS **positively**. H2 would be a
  human-tumour phenomenon, most likely inflammatory-microenvironment-driven
  despite the purity and leukocyte adjustment.

### 7.4 Traps specific to the cross-species step

1. **No ortholog function is called in this repo, in either direction.** The
   reconciliation is a comparison of **conclusions**, not of gene lists. Any
   gene-level mapping happens in the mouse repo or a third place. The tripwire
   here stays clean.
2. **Human and mouse MitoCarta are different inventories.** "Regulon membership"
   is not automatically the same object across species; it must be rebuilt
   species-natively, and the count of members in each half will differ.
3. **Never pool GSVA or mitoPPS values across cohorts** - and a species is a
   cohort. Compare correlations, patterns and rankings.
4. **The composition null must be rebuilt species-natively too.** Importing the
   human null would make the comparison circular.
5. **The human confounds are human.** Purity and leukocyte fraction were the
   ones that mattered here; a mouse model has different ones and the absence of
   these two is not the absence of confounding.
6. **The sibling study is closed.** `myc_human_validation` at `d3ac60e` is
   pre-registered, found nothing supported, does not go in the paper, and is not
   reopened. Its result is trustworthy in a way nothing in this repo is,
   *because it was declared first* - so N2 outranks anything exploratory that
   appears to contradict it.
7. **Different measurement types are not a disagreement.** The human arm reads
   bulk tumour RNA. If the mouse arm reads protein, function or perturbation,
   a null there does not refute a correlation here and a positive there does not
   confirm one. Say which claim each is addressing before comparing them.

---

## 8. Superseded statements - read this, not that

| source | superseded statement | read instead |
|---|---|---|
| `paper_opening_human.md`, rung 4 and the fourth sentence | "the ordering is no steeper than any matched gene set ... the bulk transcriptome reports mitochondrial content rather than a selectively engaged death programme" | section 3.1 **and** 3.2. True of the regulon half, false of the cytosolic half |
| `paper_opening_human.md`, open item 1 | "the ~1.3 SD residue ... would be settled by a curated non-apoptotic comparator" | **answered.** `e14_curated_comparators.md` |
| `handoff_evening.md`, the finding paragraph | "...and it is not specific to apoptosis" | section 1 of this document |
| `e10_machinery_and_priming.md` R1 | the 0.453 split as a localisation-organised death programme | the correction at the head of that note, then 3.1 here. It is regulon membership, and the split itself is section 3.2 |
| `e10_machinery_and_priming.md` R1, "protein localisation" | localisation language throughout | co-regulation. `priming_interaction_tested.md` C1 |
| `e10_machinery_and_priming.md` R4 | the Pearson spot check | withdrawn; superseded by E09 |
| `e10_machinery_and_priming.md` R5 | "6 ratios gain on OXPHOS, 3 on MYC" | adjusted: **0 on OXPHOS, 5 on MYC** |
| any note before 2026-09-02 evening | "outer mitochondrial membrane", "acts at the mitochondrion" | "nuclear-encoded mitochondrial regulon" |

---

## 9. Still open on the human side

1. **A cytosolic stress programme that is not apoptotic**, as the next
   comparator. If a proteotoxic or integrated-stress-response module also runs
   against OXPHOS, then H2's property is "cytosolic and stress-responsive"
   rather than "cytosolic and apoptotic". **This is the single most valuable
   next human analysis** and it should be done before H2 is leaned on hard in
   either arm.

   - **Its other half: what the shared correlate of the cytosolic genes
     actually is.** `data/collectri_human/` is pinned, snapshotted and has
     never been used for this. Two forms - which CollecTRI regulons contain the
     44, and whether per-sample regulon activity tracks either half. The
     comparator asks whether the property is "cytosolic and apoptotic" or
     "cytosolic and stress-responsive"; the regulon test asks what the
     correlate is. **Falsifier for the "two regulons" language:** if no regulon
     separates the halves beyond the mitochondrial ones, the phrase goes and
     "different correlate" replaces it. `e15_two_axis_gene_view.md` V5.

2. **`LumA` alone**, the homogeneous stratum where stromal and immune
   composition vary least. Named in E14 as the falsifier for H2.
3. Score `COLLECTRI_MYC_STIM` (739 genes, in the snapshot, never scored) as a
   fifth base; drop `ELLWOOD` and recompute the entanglement slope. Both need a
   pipeline re-run.
4. The ER-negative fatty-acid-oxidation reversal, untouched.
5. A stranded total-RNA dataset for `MT-ND6` and `CO1`/`CO2`.

**Out of scope until reopened:** MCbiclust / forkscale, survival, treatment,
METABRIC, DepMap, causal modelling on human data.
