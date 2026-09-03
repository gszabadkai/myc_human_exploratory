---
date: 2026-09-02
script: scripts/E10_machinery_measures_and_priming.R
status: WRITTEN AND VERIFIED, NOT YET RUN BY THE AUTHOR
posture: EXPLORATORY - nothing here is pre-registered
---

# E10 - the canonical machinery re-annotated, on both axes, under both measures,
# and the BCL2-family priming ratios

Answers handoff open item **5.4** plus the author's three follow-ups. Numbers
below come from a dry run of the finished script against the real cohorts; the
figures and tables land in the repo only when the author sources it.

**Read every number here as one cell of a large grid.** 44 genes x 6 axes x 2
measures x 2 cohorts of gene cells and 39 x 6 x 2 x 2 ratio cells. Structure and
cross-cohort agreement are the result; no interval here is one.

---

> **CORRECTION, 2026-09-02, from E11 section 3.1.** R1 below is arithmetically
> right and its interpretation was wrong. The 0.453 split is what ANY 20
> MitoCarta plus 24 non-MitoCarta genes of this expression profile produce - a
> composition-matched null gives 0.40-0.46 and the observed value sits inside
> it (z = 0.06 to 0.44). Under a stricter null pool that removes every OXPHOS
> and mitoribosome gene the null falls to 0.29-0.37 and the observed sits about
> one SD above (z = 0.83 to 1.16), consistently in all four cohort x adjustment
> cells but not separably. **S6 replicates under an independent annotation, and
> it is not specific to apoptosis.** Read R1 as "the machinery behaves like
> mitochondrial genes of its expression profile", not as "the machinery has a
> localisation-organised death programme". E10 did not run that null; E11 does.

## R1. S6 SURVIVES an independent localisation, and it survives exactly

The claim under test (S6): among the 44 genes carrying a `family_pathway` label,
the split that predicts the sign of the OXPHOS correlation is **where the gene
acts** (0.453), not its **annotated direction of effect** (0.225). Both labels
came out of the same CSV, so both could be one curator's habit.

Replacing the localisation with **MitoCarta 3.0** - mass spectrometry, GFP
imaging and targeting-signal prediction, which has never heard of apoptosis:

| predictor | source | Spearman with the OXPHOS rho |
|---|---|---|
| acts at the mitochondrion | the curation under test | **0.453** |
| **in MitoCarta 3.0** | **independent** | **0.453** |
| in a Reactome INTRINSIC pathway | independent | 0.110 |
| in a Reactome EXTRINSIC pathway | independent | -0.428 |
| annotated pro-death | the curation under test | 0.225 |

MitoCarta and the curation agree on **43 of 44** genes; the single disagreement
is `BCL2A1`, which MitoCarta places at the outer membrane and the curation does
not. The independent predictor is not merely close to the original, it is
numerically identical, and both beat direction of effect by a factor of two.

**S6 stands as a replication, and it is now a claim about protein localisation
rather than about one curation.** Whether protein localisation says anything
about APOPTOSIS specifically is a different question, and the correction above
is the answer: it does not.

### R1a. The compartment ladder

MitoCarta's sub-mitochondrial call orders the 44 monotonically:

| compartment | n | median rho with OXPHOS |
|---|---|---|
| not in MitoCarta | 24 | -0.100 |
| outer membrane (MOM) | 13 | +0.085 |
| intermembrane space (IMS) | 5 | +0.115 |
| inner membrane (MIM) | 2 | **+0.443** |

The ladder runs the right way and the deeper the protein sits, the more it
tracks OXPHOS - but MIM is `AIFM1` and `HTRA2` alone. **Two genes are a
direction to test, not a gradient.** Named as the weakest link in R1.

---

## R2. APAF1, the named anomaly, is resolved - by the distinction between
## what a gene does and where it is

`APAF1` sits at -0.452, the most negative of the 44, inside a
`mitochondrial / intrinsic` module whose median is +0.115.

- **Reactome keeps it an anomaly.** It is in the apoptosome and cytochrome-c
  sub-pathways, so Reactome calls it intrinsic, and it is still the outlier.
- **MitoCarta removes the anomaly entirely.** APAF1 is **not in MitoCarta** -
  the apoptosome assembles in the cytosol. On the localisation split it is
  correctly on the non-mitochondrial side, where the median is -0.100 and
  -0.452 is unremarkable.

So the two annotations disagree about APAF1 *because they answer different
questions*, and the one that predicts the correlation is localisation. This is
the same conclusion as R1 arrived at from the other end.

### R2a. Reactome's apoptosis annotation has a systematic hole

**8 of the 44 are in no Reactome apoptosis pathway at all**: `AIFM1`, `BCL2A1`,
`BOK`, `ENDOG`, `HRK`, `HTRA2`, `MCL1`, `NOL3`. Five of those eight are
MitoCarta-positive and they carry the **highest** median of any module (+0.214).
Reactome's reaction-level curation is built around the caspase cascade and
misses the mitochondrial matrix effectors - which are precisely the genes
carrying the effect. This is why Reactome INTRINSIC membership scores only
0.110: its "intrinsic" set includes cytosolic apoptosome components and excludes
matrix proteins.

**Reactome EXTRINSIC membership at -0.428 is the useful half of the
re-annotation** and is close in magnitude to the localisation split.

---

## R3. The split is an OXPHOS phenomenon. On the MYC axis it nearly vanishes

New, from Q-b. The same 44 genes against `FELSHER__MITOSTRIP`:

| axis | measure | vs MitoCarta | vs pro-death | median mito | median non-mito |
|---|---|---|---|---|---|
| OXPHOS | spearman | **0.453** | 0.225 | +0.147 | -0.100 |
| OXPHOS | pearson | 0.442 | 0.233 | +0.156 | -0.101 |
| MYC | spearman | **0.140** | 0.149 | +0.048 | +0.006 |
| MYC | pearson | 0.101 | 0.145 | +0.056 | +0.021 |

On MYC the localisation split collapses to 0.14 and stops beating direction of
effect. **Whatever organises the apoptotic machinery on this plane is organised
by the OXPHOS axis, not by MYC activity** - which is consistent with D3/S1,
where the BCL2 column against MYC turned out to be between-subtype pooling while
the OXPHOS column was stratum-stable.

Two cells on the MYC panel must not be read: `BIRC5` (+0.601) and `TP53` are
HALLMARK E2F/G2M genes, and a MYC signature is entangled with proliferation by
construction (trap 3). They are crossed on figs 2 and 3.

---

## R4. WITHDRAWN 2026-09-02 - superseded by E09

R4 was a 44-gene spot check finding that Pearson changed almost nothing (largest
gap 0.102, everything else under 0.07, the S6 split 0.442 against 0.453). E09
then answered the same question properly over 220 pairs and all four
instruments: the measures correlate at 0.996, bicor sides with Spearman in all
twelve of the largest disagreements, and nothing is non-monotone. See
`docs/2026-09-02_e09_correlation_measures.md`.

**The Pearson computation and its four figures were removed from E10 on
2026-09-02**, on the author's instruction and because E09 makes them redundant.
The script is Spearman throughout and emits 7 figures rather than 10. The
finding "the measure does not matter" is E09's, not this script's - and E09
qualifies it: the measure DOES matter on mitoPPS, where the mean departure is
0.029 and reaches 0.093.

---

## R5. The priming ratios mostly do NOT beat their own component genes

**35 ratios** (7 pro x 5 anti). The falsifier was written before the numbers:
`gain` = |rho of the ratio| - |rho of its stronger component|, and a ratio is
only worth reporting as a ratio where gain is positive **in both cohorts**.

| axis | n ratios | gain > 0 in both | median gain TCGA | median gain SCAN-B |
|---|---|---|---|---|
| OXPHOS | 35 | 6 | -0.125 | -0.062 |
| MYC | 35 | 3 | -0.062 | -0.020 |

**The median ratio is worse than one of its two genes on both axes.** The
strongest-looking cells are the components talking: the best OXPHOS ratio is
`BAD/MCL1` at +0.411, and `BAD` alone is +0.402.

The few that do gain in both cohorts, and by how little:

- **OXPHOS**: `BBC3/BCL2` (+0.094 / +0.059), `BIK/BCL2` (+0.055 / +0.073),
  `BID/MCL1` (+0.069 / +0.047), `BAD/MCL1` (+0.011 / +0.008),
  `BMF/BCL2L1` (+0.008 / +0.067), `PMAIP1/BCL2A1` (+0.002 / +0.005)
- **MYC**: `PMAIP1/BCL2A1` (+0.049 / +0.029), `BMF/BCL2A1` (+0.033 / +0.044),
  `PMAIP1/MCL1` (+0.011 / +0.001)

Every one of those pairs a pro-apoptotic gene with an anti-apoptotic gene that
moves the *opposite* way on that axis, which is exactly when a ratio can add
rather than cancel. **The three or four with gains above ~0.05 are the only ones
worth carrying forward, and they are hypothesis-generating.**

### R5a. What the components do on their own (Spearman, mean over cohorts)

| gene | side | MYC | OXPHOS |
|---|---|---|---|
| BCL2L1 | anti | -0.058 | **+0.342** |
| BCL2A1 | anti | +0.198 | +0.033 |
| BCL2L2 | anti | -0.177 | -0.125 |
| BCL2 | anti | **-0.328** | -0.196 |
| MCL1 | anti | 0.000 | -0.253 |
| BAD | pro | -0.037 | **+0.402** |
| BID | pro | **+0.385** | +0.277 |
| BIK | pro | -0.016 | +0.244 |
| BBC3 | pro | -0.061 | +0.178 |
| PMAIP1 | pro | -0.097 | -0.068 |
| BCL2L11 | pro | -0.079 | -0.171 |
| BMF | pro | -0.105 | -0.214 |

The pro-apoptotic list does not move as a block and neither does the
anti-apoptotic one: `BCL2L1` (anti) sits at +0.342 on OXPHOS while `BCL2` (anti)
sits at -0.196, and `BAD` (pro) at +0.402 while `BMF` (pro) at -0.214. **The
functional dichotomy the ratios assume is not present in the data at gene
level**, which is the same conclusion S6 reached for the wider machinery and the
reason most of the ratios cancel.

All 12 genes clear the 25th expression percentile in both cohorts, so none of
this is a low-expression artefact.

### R5b. `BCL2L2` - RESOLVED 2026-09-02

It was originally supplied on both lists. The first run flagged it as
canonically **anti**-apoptotic - Bcl-w is a multidomain guardian, not a BH3-only
sensitiser - and the author confirmed the pro-side entry was an error. **It is
anti-apoptotic only.** The grid is 7 x 5 = 35 with no self-pair and no empty
cell, and `BCL2L2/BCL2A1` leaves the MYC gainer list, taking it from 4 to 3.

If a seventh pro-apoptotic gene is wanted, `HRK` and `BOK` are both already
scored in section 4. `HRK` would need its low-expression flag carried through -
it is below the 25th percentile in both cohorts.

### R5c. The two tests, put on the heatmaps - ADDED 2026-09-03

Author's request: stop making a reader hold figure 4 in their head while looking
at figure 3 or figure 6. Both heatmaps now carry the test as a cell mark.

| mark | meaning |
|---|---|
| `*` | in **this cohort**, `gain > 0` **and** `\|rho\| >= 0.30` - the ratio beats the stronger of its own two genes, and the effect is not small |
| **heavy black border** | the same ratio does that in **both cohorts**, at the same stratum and on the same axis |

The 0.30 floor is **hand-drawn, is not a test, and moving it moves the marks.**
It is `|rho|` rather than `rho` so a ratio running strongly the other way is
marked too - five of the marked cells are negative, and hiding them would make
the mark mean "large and positive" while the caption said "large".

**What the marks show, and it is a cleaner statement than R5's table:**

1. **Not one MYC cell is marked**, in either heatmap, at any stratum, in either
   cohort. All 22 marked cells on figure 6 and all 7 on figure 3 are OXPHOS.
2. **The pooled heatmap has no bordered cell at all, and the two halves of that
   failure pull opposite ways.** The 5 ratios that beat their parts in both
   cohorts are all on MYC and their `|rho|` tops out at **0.27**; the 7 cells
   that reach 0.30 are all on OXPHOS and not one of them beats its parts in both
   cohorts. *Pooled, a priming ratio is either replicable or large, never both.*
3. **Three cells clear both bars in both cohorts, and all three are Basal
   OXPHOS**: `BBC3/BCL2`, `BID/BCL2`, `BBC3/MCL1`. Nothing in Luminal, nothing
   pooled.
4. **Read (3) against where the marks are densest.** Basal is 171 TCGA and 317
   SCAN-B samples, so it is also where the intervals are widest - about +/- 0.15
   in TCGA. A mark in Basal is a weaker claim than the same mark would be in
   Luminal, and three cells out of 70 is a hypothesis, not a result.

Saved as `priming_marked`, `strata_marked`, `mark_summary`, `mark_both_list` and
`pooled_two_way`; the constant is `settings$mark_min_abs_rho`.

**Also fixed in the same pass:** figure 3's subtitle had been clipped, not
wrapped, since the figure was written - ggplot clips a long subtitle the same
way it clips a caption, and every run before 2026-09-03 silently lost its last
three words.

---

## R6. Housekeeping the script does that is worth knowing about

- **`CYCS` is 1 of the 89 genes in the `OXPHOS subunits` arm.** Its OXPHOS value
  (+0.425) is partly a correlation with itself. E08 fig6 did not mark this; all
  four E10 machinery panels do.
- **`FELSHER__MITOSTRIP` contains none of the 44 and none of the 12**, so every
  MYC panel here is clean of self-overlap. `M_b__MITOSTRIP` contains seven
  (`BIRC5 CFLAR FAS FASLG TNF TNFRSF10B TP53`) and is therefore tabled and never
  plotted.
- **`FASLG` and `HRK` sit below the 25th expression percentile** and are starred
  on every machinery panel.
- The script **asserts** that its Spearman values reproduce E08's per-gene
  numbers exactly (max difference 0e+00 over 44 genes) before reading anything.
  E08 computed them on the raw linear matrix through the atlas engine; E10
  computes them on log2(x + 1) through `stats::cor`.

---

## What would falsify each of these

- **R1** DIED, in the sense that mattered, on 2026-09-02: the split is inside a
  composition-matched null (E11 section 3.1). What survives is the replication
  across annotation sources. It would die further if the MIM/IMS ladder does
  not survive being widened past two genes.
- **R3** dies if the split reappears on MYC under an estimator other than
  `FELSHER__MITOSTRIP`. The panel is scored; this is a filter away.
- **R5** is already stated as a negative and needs no falsifier. The positive
  residue - three or four ratios with gain above 0.05 - dies if it does not
  survive stratification, which E08's strata machinery can test directly.

## What was NOT done

- **No strata.** Every number here is all-samples. Given D3/S1, the MYC column
  in particular deserves the LumA / Luminal / Basal breakdown before anything in
  R3 is believed.
- **No purity or proliferation adjustment.** Raw only, both cohorts.
- **No expression-matched null.** These are named genes chosen a priori, not
  sets, so the E05 null does not apply - but the 39 ratios are a panel and a
  permutation over pro/anti labels would sharpen R5.
