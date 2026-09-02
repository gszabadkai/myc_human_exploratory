---
date: 2026-09-02
script: scripts/E14_curated_comparators.R
status: WRITTEN AND VERIFIED BY DRY RUN, NOT YET RUN BY THE AUTHOR
posture: EXPLORATORY - nothing here is pre-registered
answers: open item 1 of docs/2026-09-02_handoff_evening.md
revises: the last clause of the one-paragraph finding in that handoff
---

# E14 - the residue, closed with a named comparator instead of more draws

Open item 1 asked whether "find another control" meant a second random null.
It did not. The random null and a curated comparator ask different questions,
and only the second is the one the paper's sentence makes.

| | |
|---|---|
| the random null asks | would ANY 20 mitochondrial plus 24 cytosolic genes of this expression profile and compartment split this way? |
| the sentence asks | does any other PROGRAMME spanning the outer membrane split this way? |

A random draw's members share nothing but expression and compartment. A real
pathway's members are **co-regulated with each other**, and that is precisely
the property that could make a split. The null used until now did not have it.

## C0. Protein import cannot be the comparator, and the reason is structural

It was the first candidate and it fails on a fact, not on a judgement.
`REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT` is 62 genes in the TCGA matrix and
**62 of 62 are in MitoCarta** - TOM, TIM, SAM, PAM, the presequence proteases,
`TOMM34` with them. There is no cytosolic half, so the split statistic cannot
be computed at all. The same is true of cristae formation (31 of 31) and
mitochondrial calcium transport (21 of 22). All three are **level** comparators
and are used as such in section 6, never as split comparators.

**Mitophagy is the one that works.** `REACTOME_MITOPHAGY`, 38 genes: 16
MitoCarta (15 outer membrane, 1 matrix) against 22 cytosolic - ATG5, ATG12,
ATG9A, LC3, OPTN, SQSTM1, ULK1, TBK1, CK2, SRC, the ubiquitin conjugation
machinery - and **zero genes in common with the 44**. It has the machinery's
shape, and its cytosolic members translocate to the outer membrane on
activation, which is the same hazard `BAX`, `BID` and `BAD` carry.

## C1. THE ANSWER. The machinery out-splits every programme tested

OXPHOS axis, proliferation-adjusted, the `split` statistic. Random null in the
last column for scale.

| set | mito/cyt | TCGA | SCAN-B | matched null (TCGA) |
|---|---|---|---|---|
| **apoptotic machinery (44)** | 20/24 | **0.453** | **0.489** | 0.276 +/- 0.131 |
| Fe-S cluster assembly | 13/10 | 0.397 | -0.013 | 0.331 +/- 0.185 |
| mitophagy | 16/22 | 0.112 | 0.219 | 0.168 +/- 0.149 |
| mitophagy, PINK1/PRKN only | 14/17 | 0.065 | 0.196 | 0.182 +/- 0.163 |
| isozyme pairs (declared ceiling) | 33/16 | -0.092 | -0.022 | 0.304 +/- 0.127 |

The falsifier written into the script header before the run required a gap over
mitophagy of more than 0.12 in **both** cohorts. It is **0.341 (TCGA) and 0.270
(SCAN-B)**. The verdict the script prints is `APOPTOSIS-SPECIFIC in both
cohorts`.

**And note what the null column does.** The random matched null sits at
0.28-0.32 - *higher than three of the four real programmes*. That is why more
draws could never have settled this: a random cytosolic gene is an average
cytosolic gene, while a real pathway's cytosolic arm is a co-regulated,
growth-associated gene that tracks OXPHOS positively, which SHRINKS the split.
The synthetic null was not lenient, it was pointing the wrong way.

## C2. THE MECHANISM, AND IT IS NOT WHERE E11 LOOKED

Section 5 splits the statistic into its two halves against their own nulls.
Mean rho with OXPHOS, adjusted, observed (null) z:

| set | mitochondrial half | cytosolic half |
|---|---|---|
| **apoptotic machinery** | +0.167 (+0.139) **z +0.48** | **-0.106** (-0.036) **z -1.22** |
| mitophagy | +0.210 (+0.105) z +1.57 | +0.143 (-0.008) z +2.41 |
| mitophagy, PINK1/PRKN | +0.185 (+0.111) z +1.08 | +0.155 (+0.002) z +2.16 |
| Fe-S cluster assembly | +0.313 (+0.168) z +1.96 | +0.118 (-0.035) z +1.74 |
| isozyme pairs | +0.146 (+0.186) z -0.90 | +0.186 (-0.013) z +2.70 |

(TCGA; SCAN-B agrees in sign on every cell.)

**The machinery's mitochondrial half is unremarkable.** It sits at its null,
and mitophagy's mitochondrial half is *higher*. Everything E11 measured about
the mitochondrial 20 stands, and stands as a negative.

**The machinery's cytosolic half is the anomaly, and it is the only one.** Every
other programme's cytosolic arm tracks OXPHOS positively and above its own
null; the machinery's runs against it. The gap between the two cytosolic arms
is 0.20-0.25 in rho units, in both cohorts.

The ladder says the same thing from a different direction. Against 30 MitoCarta
leaf pathways:

| | TCGA | SCAN-B |
|---|---|---|
| machinery, mitochondrial half | 63rd percentile | 63rd percentile |
| machinery, cytosolic half | **0th** | **0th** |

## C3. THE CONFOUND, TESTED IN THE SAME SCRIPT - it is not infiltrate

The 24 cytosolic genes are `FAS`, `FASLG`, `TNF`, `TNFRSF1A`, `TNFRSF10A/B`,
`TRADD`, `FADD`, `CASP10`, `CFLAR`, `BIRC2/3`, `XIAP`, `NFKB1`, `RELA` - as
much an inflammatory module as an apoptotic one. Trap 2 puts `rho(OXPHOS
subunits, leukocyte fraction)` at -0.158, so infiltrate alone would give a
negative cytosolic half, a large split and an "apoptosis-specific" verdict with
no apoptosis in it. **The comparator cannot control for this - mitophagy's
cytosolic arm is housekeeping and carries no immune signal. Only the covariate
can.** Section 5.1 was declared in the header before the run for that reason.

| | TCGA |
|---|---|
| cytosolic half, proliferation-adjusted | -0.109 |
| **the same, plus purity and leukocyte fraction (n = 1,007)** | **-0.091** |

It barely moves. And removing the death-receptor module entirely still leaves
14 genes at -0.111 (TCGA) and -0.087 (SCAN-B), **9 of 14 negative in both**.

By module, the most negative group is not the immune one:

| module | n | TCGA | SCAN-B |
|---|---|---|---|
| mitochondrial / intrinsic, cytosolic (`APAF1`, `BMF`, `HRK`) | 3 | **-0.332** | **-0.257** |
| IAP / NF-kB | 7 | -0.128 | -0.077 |
| death receptor / extrinsic | 10 | -0.098 | -0.104 |
| effector caspase | 2 | +0.057 | -0.017 |

SCAN-B has no purity estimate and it is never imputed, so the covariate test is
TCGA-only and the module breakdown is what both cohorts carry.

## C4. Two things worth writing down because they were not predicted

**The isozyme ceiling failed, and the miss is informative.** The script header
predicted a large split there - same reaction either side of the membrane, so
compartment is the only thing that varies, and several members are TCA enzymes
or mitochondrial tRNA synthetases. It came out at **-0.09 / -0.02**: no split
at all. Both halves are positive (mito +0.15/+0.18, cytosolic +0.19/+0.19). So
**being in the nuclear-encoded mitochondrial regulon does not by itself produce
a split** when the cytosolic counterpart does the same job. The prediction is
left in the script header as written; this is the record of it being wrong.

**MYC orders mitophagy and does not order apoptosis.** On the MYC axis
mitophagy splits at 0.233 (TCGA) / 0.306 (SCAN-B), above its null (z +0.99 /
+1.86, `z_md` +0.68 / +2.62), while the machinery sits at 0.187 / 0.137, at or
below its own (z -0.24 / -0.16). This is an independent answer to the worry E11
P2 raised: the adjusted MYC estimator does order a mitochondrial programme's
genes by compartment. Just not the apoptotic one.

## C5. What this changes in the paper, and what it does not

The handoff's one-paragraph finding ends: *"and the ordering is no steeper than
an expression- and compartment-matched gene set gives, so it is real, it is
OXPHOS's, and it is not specific to apoptosis."* **The last clause has to go**,
and the sentence that replaces it is narrower and better:

> The machinery's mitochondrial members sit where any mitochondrial
> programme's members sit. What no other programme does is send its cytosolic
> members the other way: mitophagy, iron-sulfur cluster assembly and
> compartment-matched isozyme pairs all have cytosolic arms that track OXPHOS
> positively, and the apoptotic machinery's runs against it, in both cohorts,
> after proliferation, purity and leukocyte fraction.

Unchanged: E11 P4/P4a's arithmetic (the 1.3 SD against a random null is still
1.35), the conditioning result (OXPHOS rather than MYC), P3's restatement, the
mitoribosome control, and every caveat about mediation not being identifiable.

## C5a. The supplementary panel

`E14_fig6_specificity_one_panel` is the whole argument in one square panel:
each programme is a point whose **x** is what its mitochondrial members do with
OXPHOS and whose **y** is what its cytosolic members do. Both axes are the same
quantity on the same scale, so "they agree horizontally and disagree
vertically" is a property of the picture and not of how it was drawn, and both
zeros are in view.

Four programmes, not five - the PINK1/PRKN subset sits almost on top of
mitophagy, because it *is* mitophagy minus the receptor arm, so it is named in
the caption with its values instead of plotted. Label positions are fixed
coordinates rather than repelled: a seed-driven layout moves between runs and a
figure that goes in a paper must not.

`fig1` (the head-to-head), `fig2` (the two halves) and `fig5` (the infiltrate
falsifier) make the same argument in three steps and stay the main-text route.
**`fig6` must not be shown without `fig5` cited in its legend** - the panel's
whole weight rests on the negative not being infiltrate.

## C6. What would falsify the new claim

Written before the next analysis, as the rule requires.

1. **A non-death cytosolic stress programme with a negative arm.** If a
   cytosolic pathway that is stress-responsive but not apoptotic - a
   proteotoxic or integrated-stress-response module - also runs against OXPHOS,
   then the property is "cytosolic and stress-responsive", not "cytosolic and
   apoptotic". This is the single most valuable next comparator.
2. **It dies if it is a handful of genes.** It is not, on this run: 9 of 14
   non-death-receptor cytosolic genes are negative in both cohorts. Re-check
   after any change to the 44.
3. **It dies if it does not survive inside one subtype.** LumA alone has never
   been run (open item 3) and is the homogeneous stratum where stromal and
   immune composition vary least. Do this before writing anything.
4. **Four comparators cannot make a p-value.** `Fe-S cluster assembly` is
   0.397 in TCGA and -0.013 in SCAN-B - the comparator field is noisy, and the
   machinery is the only set here that reproduces tightly (0.453 / 0.489). The
   claim rests on the ranking reproducing, not on a test.

## C7. Traps this script added

1. **`.gene_rows` names its matrix by the ORIGINAL symbol.** Taking
   `match(rownames(mat), rownames(L))` after it and dropping the NAs makes the
   index vector shorter than the compartment vector beside it; `split()` then
   RECYCLES the labels with a warning and no error, and the compartments come
   out attached to the wrong genes. Caught in the first dry run. `.row_index()`
   in E14 resolves and indexes in one step and every caller asserts the
   lengths. **E11 sections 3.1, 3.2, 4.2 and 4.3 use the same pattern** - it is
   harmless there only because all 44 genes resolve identically in both
   matrices, which is luck, not design.
2. **Trap 7 reaches sets invented after E02.** The pinned symbol map covers the
   arm sets, the MitoCarta pathways and the estimators - not new lists. Eleven
   cytosolic tRNA synthetases are `AARS`, `DARS`, `LARS` in SCAN-B's 2014 build,
   so the isozyme set arrived as 60 genes in TCGA and 49 in SCAN-B **with the
   whole loss on one half**. E14 now cuts every comparator to the genes
   resolvable in both matrices before measuring anything, and prints what it
   cut.
3. **Declaring six criteria and implementing four.** The header listed six
   eligibility criteria; the filter checked three. Criterion 6 is now in the
   table. It changes no verdict - the four eligible sets carry zero OXPHOS-arm
   and zero mitoribosome genes - but a header that does not match its code is
   how a reader is misled.
