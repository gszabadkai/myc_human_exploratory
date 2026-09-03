# E15 - the two axes on one row per gene

**2026-09-03. EXPLORATORY. Nothing here is pre-registered.**

Script: `scripts/E15_two_axis_gene_view.R`. Written this session, **not yet run
by the author** - the numbers below are from a dry run with output paths
redirected to the scratchpad, and the repo was left untouched. Sourcing it takes
a few seconds: it reads three saved objects and computes no correlation.

Entry point for the whole human arm remains
`docs/2026-09-03_human_arm_for_mouse_reconciliation.md`. This note is narrower:
it exists because the author asked for the OXPHOS and MYC per-gene correlations
**side by side on one figure**, which no previous figure does.

---

## Why the figure did not already exist

| figure | what it shows | what it cannot show |
|---|---|---|
| `E08_fig6`, `E10_fig1` | the 44 against **OXPHOS**, one row per gene | nothing about MYC |
| `E10_fig2` | the 44 against **MYC**, rows in the same order | requires laying two pages side by side |
| `E11_fig1` | the 44 as a **cloud**, MYC on x and OXPHOS on y | which gene is which; only the shape of the cloud |

The question "does a gene move between the axes, and which ones" was
answerable from the saved tables and not from any picture. E15 is that picture.

## Where every number comes from

Nothing is recomputed. `E11$gene_tab` supplies every rho (partial Spearman on
`PROLIF_DISJOINT`, MYC axis `FELSHER__MITOSTRIP`); `E08$expr_rank` supplies the
expression percentile that stars the noisy genes; MitoCarta sheet 2 supplies the
sub-compartment label. **Section 1 asserts that the localisation split
recomputed from E11's per-gene numbers is bit-equal to the value E14 published
for the same 44 genes, in all eight cohort x axis x adjustment cells.** If a
future edit to E11 changes those numbers, E15 stops rather than drawing them.

Outputs: `results/two_axis_gene_view.rds`,
`outputs/tables/E15_gene_rho_two_axes.csv`, five figures.

---

## V1. What the figures show

`E15_fig1` - one bar per gene, hollow point = MYC, filled point = OXPHOS, rows
sorted by the cross-cohort mean of the difference, both cohorts, same order.
`E15_fig2` - the same 44 as a numbered heatmap, raw and adjusted, halves
stacked. `E15_fig3` - the arithmetic check in V2. `E15_fig4` - the two halves
as bars, mean +/- SD, with every gene drawn on top (V1b). `E15_fig5` - the
difference alone, gene by gene, both cohorts in one panel (V1c).

Adjusted, TCGA / SCAN-B. The split is Spearman of the quantity with MitoCarta
membership across the 44 genes - the study's anchor statistic:

| quantity | split | median, mito half | median, cytosolic half |
|---|---|---|---|
| **MYC** | 0.187 / 0.137 | +0.117 / +0.047 | +0.011 / -0.005 |
| **OXPHOS** | 0.453 / 0.489 | +0.249 / +0.197 | -0.164 / -0.145 |
| **OXPHOS - MYC** | 0.471 / 0.485 | +0.058 / +0.063 | -0.098 / -0.144 |

**The two middle rows are the same statement twice.** Subtracting MYC from
OXPHOS leaves the ordering where it was - 0.453 to 0.471 in TCGA, 0.489 to 0.485
in SCAN-B - because MYC carries almost none of it. That is a cleaner way to say
"OXPHOS orders these genes and MYC does not" than comparing two separate splits,
because it is computed inside each gene rather than across two summaries.

The difference **replicates**: Spearman 0.83 between cohorts, 38 of 44 agreeing
in sign. The six that do not - `BCL2L11`, `BOK`, `CASP10`, `CASP8`, `DIABLO`,
`XIAP` - are drawn with a dotted bar and are all small-gap genes.

## V1b. The same thing as two bars - `E15_fig4`

Requested 2026-09-03: mean +/- SD of the proliferation-adjusted per-gene
difference, mitochondrial against cytosolic. The two components are drawn beside
it on the same y axis, because a reader shown only the difference cannot tell
whether it comes from OXPHOS moving or from MYC moving, and that is the point.

Adjusted, TCGA / SCAN-B. `sep` is the gap between the two group means:

| | mitochondrial (n=20) | cytosolic (n=24) | sep |
|---|---|---|---|
| **MYC** | +0.078 +/- 0.196 / +0.055 +/- 0.152 | +0.002 +/- 0.183 / +0.003 +/- 0.161 | **+0.076 / +0.053** |
| **OXPHOS** | +0.167 +/- 0.286 / +0.143 +/- 0.226 | -0.106 +/- 0.283 / -0.094 +/- 0.209 | **+0.273 / +0.237** |
| **OXPHOS - MYC** | +0.089 +/- 0.188 / +0.088 +/- 0.171 | -0.108 +/- 0.186 / -0.096 +/- 0.161 | **+0.197 / +0.184** |

Three readings, in order of how much weight they take:

1. **MYC's two bars are on the same side of zero and 0.05-0.08 apart. OXPHOS's
   straddle zero and are 0.24-0.27 apart.** MYC separates the halves three to
   four times less, and does not reverse their sign.
2. **The error bar is +/- one SD of the genes in the group, not an uncertainty
   on the mean.** At 0.16-0.19 it is nearly as wide as the distance between the
   means, and the bars overlap: **these halves separate on average and not gene
   by gene.** Every gene is drawn on top so that is visible rather than
   inferred. The SE of each mean is about a fifth of the SD (0.033-0.042 for the
   difference).
3. **A mean and a rank statistic disagree slightly here, and both are right.**
   Subtracting MYC costs the mean-scale separation +0.08 / +0.05 (0.273 to
   0.197) while it costs the rank-based split nothing (0.453 to 0.471; 0.489 to
   0.485). Means subtract exactly and ranks do not - the two weight genes
   differently. **Neither should be quoted as the other**, and both are in
   `bar_diff` and `split_tab`.

Robustness: deleting all eight flagged genes - `CYCS`, `TP53`, `BIRC5` and the
five below the 25th expression percentile - leaves the means at +0.089 / +0.069
(mitochondrial) and -0.107 / -0.091 (cytosolic). The difference is not carried
by them. **No test is reported**, here or anywhere in E15: two groups of a
curated 44, in a study whose atlas is a grid of thousands of cells. The reading
is the direction, its size relative to the spread, and that both cohorts show
it.

## V1c. The difference alone - `E15_fig5`

Requested 2026-09-03, then revised: a Cleveland plot of the differences, gene
names on the x axis coloured by MitoCarta membership, **both cohorts in one
panel**, no droplines to zero, narrower. It is figure 1 with the two component
points removed; what is left is the quantity the argument is about, and dropping
the components buys the room to put both cohorts together, so the pair of points
per gene reads as the replication directly rather than by looking between two
panels.

Two construction decisions worth recording.

**One gene order, shared** - the cross-cohort mean, as in figures 1, 2 and 4 -
rather than each cohort sorting on its own values. Per-cohort sorting gives a
tidier monotone descent and makes the two incomparable, which is the wrong trade
in a study whose unit of evidence is cross-cohort agreement. Neither cohort is
therefore monotone; **where the descent breaks is where the cohorts differ, and
that departure is the replication being displayed.**

**The gene names are drawn as data, not as axis text.** Colouring 44 axis labels
needs a vectorised `element_text(colour = ...)`, which ggplot2 warns is
unsupported and whose recycling order is not guaranteed - a silently re-ordered
vector would colour every name wrongly and still render. Drawing them from the
same rows that position the points cannot come apart from them. This is the same
rule figure 1 follows with its margin squares.

What the eye gets from it that the summaries do not: the red names crowd the
left of the axis and the grey ones the right, and the exceptions are legible by
name rather than buried in a mean. `NOL3` and `BIRC5` are non-MitoCarta and sit
among the largest positive differences; `MCL1`, `BCL2A1`, `BCL2L11` and `BID`
are MitoCarta members at the negative end.

## V2. The one arithmetic caveat, and it is real

`rho(g, OXPHOS) - rho(g, MYC)` is **not a difference in a common unit.** The 44
correlations spread wider on OXPHOS (SD 0.313 / 0.245) than on MYC (0.191 /
0.157), so a gene earns part of its bar by sitting on the wider axis.

Dividing each axis by the SD of its own 44 values before differencing takes the
split from **0.471 / 0.485 to 0.349 / 0.356**. So roughly a quarter of the
apparent difference is the spread, and the rest is not - and what remains still
sits far above the MYC row (0.187 / 0.137). Both versions are saved; `gap_z` is
in the CSV.

**This check was run while the figures were being designed, so its answer was
known before it was written down.** It is a caveat quantified, not a prediction
tested, and it is not a result.

## V3. The structure inside each half

Adjusted medians, TCGA / SCAN-B.

**The mitochondrial half, by MitoCarta sub-compartment.** The direction of E10's
depth ladder reproduces in both cohorts - deeper is higher on OXPHOS - but it
rests on two genes at the deep end and must not be quoted as a gradient:

| | n | OXPHOS | MYC |
|---|---|---|---|
| MOM | 13 | +0.242 / +0.091 | +0.095 / +0.033 |
| IMS | 5 | +0.247 / +0.212 | +0.002 / +0.018 |
| MIM | 2 | +0.444 / +0.384 | +0.245 / +0.156 |

MOM's OXPHOS median differs between cohorts by 0.15, which is more than the
spacing between the rungs. Read as: no reliable gradient at this resolution.

**The cytosolic half, by curation module.** The negative arm is **not** the
death receptors alone:

| | n | OXPHOS | MYC |
|---|---|---|---|
| IAP / NF-kB | 7 | -0.292 / -0.147 | -0.006 / -0.015 |
| death receptor / extrinsic | 10 | -0.167 / -0.145 | +0.039 / +0.032 |
| intrinsic, cytosolic members (`APAF1`, `BMF`, `HRK`) | 3 | -0.200 / -0.223 | -0.150 / -0.079 |
| effector caspase | 2 | +0.057 / -0.017 | -0.034 / -0.066 |

This is the gene-level version of E14's module-deletion test: removing the
death-receptor module left 14 genes at -0.111 / -0.087 with 9 of 14 negative in
both cohorts. The negative is spread across the cytosolic modules.

**Annotated direction of effect still does not do the work.** Inside the
mitochondrial half pro-death runs +0.291 / +0.207 against pro-survival -0.082 /
-0.090 on OXPHOS, which looks like a split until the pro-survival group is
counted: it is four genes, and one of them (`BCL2L1`, +0.388) has the largest
gap in the whole table. Inside the cytosolic half both directions are negative.

## V4. The reading, and what "a different regulon" can and cannot mean

The author's question is whether the two halves sit on different regulons and
what that would mean. The figure supports a narrower and more useful version of
that, and it inverts the intuitive one.

1. **The red half is on the obvious regulon, and that is not a finding.**
   MitoCarta membership at transcript level marks membership of the
   nuclear-encoded mitochondrial regulon, and an OXPHOS score is the best
   available bulk readout of that regulon. So a positive rho is what ANY member
   of it gives. E11 and E14 measured exactly that: the machinery's mitochondrial
   half sits **at** its compartment-matched null (z +0.48 / +0.60, 63rd
   percentile of 30 MitoCarta pathways), and mitophagy's mitochondrial half is
   higher. Asking which regulon those genes belong to has a known and
   uninteresting answer.

2. **The grey half is where the result lives.** No comparator programme's
   cytosolic half is negative: apoptosis is at -0.106 / -0.094 against +0.10 to
   +0.19 for mitophagy, Fe-S assembly, mitochondrial biogenesis and the isozyme
   ceiling - **0th percentile in both cohorts**. So the right question is not
   "what regulon are the mitochondrial apoptotic genes in" but **what makes the
   cytosolic apoptotic genes run AGAINST OXPHOS when other cytosolic gene sets
   do not.**

3. **On MYC there is no split to explain.** Both halves sit within 0.12 of zero
   on the adjusted MYC axis in both cohorts. Whatever orders these genes, it is
   not MYC activity, and the paper should not describe the result as a
   difference between two MYC-driven programmes.

Candidate explanations for point 2, with their current status:

| | status |
|---|---|
| immune and stromal composition - these are inflammatory-adjacent transcripts and infiltrate is anti-correlated with OXPHOS-high epithelium | **tested, survives.** Purity + leukocyte fraction moves the cytosolic mean from -0.109 to -0.091 (n=1,007, TCGA only). Not excluded - SCAN-B has no purity estimate |
| a shared cytosolic stress / NF-kB programme that is itself anti-correlated with OXPHOS | **untested.** This is the open item, below |
| expression level or set size | **addressed.** Expression-matched nulls in E11 and E14; and every comparator's cytosolic half is positive |
| the genes are apoptotic | **the claim.** It is what remains after the three above |

## V5. What would test point 2, written before it is run

`data/collectri_human/` is already snapshotted and pinned and has never been
used for this. Two forms, both cheap and both human-native:

1. **Membership.** For each of the 44, which CollecTRI regulons contain it. If
   the two halves differ only in the mitochondrial-biogenesis regulons, then
   "different regulon" is a restatement of MitoCarta membership and the phrase
   should be dropped in favour of "different correlate".
2. **Activity.** Per-sample regulon activity (`decoupleR::run_ulm`, as E02
   already does for `M_b`), correlated with the mean of each half. If a
   stress-responsive or NF-kB regulon tracks the cytosolic half's mean and runs
   negative against OXPHOS, that is the mechanism and it is measurable.

**Falsifier for the "two regulons" language:** if no regulon separates the two
halves beyond the mitochondrial ones, the language goes. The figure would still
stand - the ordering is a fact - but the explanation would not.

This is a sharper form of open item 1 in the synthesis doc (a non-apoptotic
cytosolic stress comparator), and it should be done before either arm leans on
H2 hard. **The two are complementary: the comparator asks whether the property
is "cytosolic and apoptotic" or "cytosolic and stress-responsive"; the regulon
test asks what the shared correlate actually is.**

## V6. Traps

1. **No cell of this grid is a finding.** 44 genes x 2 cohorts x 2 adjustments.
   The 44 were not selected on the statistic shown, which is better than E08's
   driver lists, but the figure is a description of an ordering.
2. **Five genes are starred** - `DIABLO`, `ENDOG`, `FASLG`, `HRK`, `TNF` sit
   below the 25th expression percentile in at least one cohort. Three of them
   carry large gaps. Their rho is largely quantisation noise against a score.
3. **`CYCS` is marked `[ox]`** - it is 1 of the 89 OXPHOS-arm genes and is
   partly correlated with itself. `TP53` and `BIRC5` are marked `[p]` - they are
   in the 318-gene proliferation covariate and are partly adjusted for
   themselves.
4. **The figure is not evidence that MYC is unrelated to these genes.** E11
   figure 2 is the control: the adjusted MYC axis still tracks the mitoribosome,
   so this is a flat MYC column and not an emptied MYC score. An emptied score
   would look identical here.
5. **Proliferation adjustment does not simply shrink the MYC column.** In TCGA
   the mitochondrial half's MYC median goes from +0.011 raw to +0.117 adjusted -
   it moves up. E11 figure 3 is that rearrangement drawn as line crossings.
6. **`mito_class` is an ordered factor on purpose.** Left as a character string
   it sorts alphabetically, puts the cytosolic block on top of figure 2, and
   runs the opposite way to figure 1 while the subtitle claims a shared row
   order.

## V7. What this changes in the paper

Nothing yet, and deliberately: every number here was already established. It
earns **supplementary figure slots** - `E15_fig1` as the per-gene version of
`E11_fig9` panel A, for a reader who wants to know which genes carry it, and
`E15_fig4` as the two-bar summary of the same thing for a reader who does not,
`E15_fig5` as the per-gene difference on its own - plus
`E15_gene_rho_two_axes.csv` as the underlying table. `E15_fig3` belongs in
the supplement beside them if the difference is quoted anywhere in the text.

`E15_fig4` is the most compact statement of the result in the whole study and
is the candidate if a main-figure panel is ever needed at this altitude. It must
carry its SD-not-SE sentence wherever it goes; without it the overlap between
the two distributions is invisible and the figure overclaims.

`E15_fig2` is a table drawn as a figure. It is the right object for a reviewer
and probably the wrong one for a page.
