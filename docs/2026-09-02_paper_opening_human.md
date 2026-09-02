---
date: 2026-09-02
purpose: how the human data can open, and what each candidate sentence is
         actually licensed by
status: framing note, not a result. Results are in the E09/E10/E11 notes.
---

# The human opening statement: what is and is not licensed

The candidate claim: *"the apoptotic machinery is driven, or at least more
strongly correlated, by OXPHOS than directly by MYC."*

**Two of its three parts survive. The third does not, and one part is stronger
than the claim as written.**

---

## The claim ladder

| # | statement | status |
|---|---|---|
| 1 | The 44 canonical apoptosis genes correlate more strongly with OXPHOS than with MYC activity | **SUPPORTED.** SD of the 44 per-gene rho 0.313 vs 0.191 (TCGA) and 0.245 vs 0.157 (SCAN-B); 66% vs 34% and 50% vs 25% exceed \|rho\| 0.2 |
| 2 | ...and this is not proliferation | **SUPPORTED.** Unchanged by partialling `PROLIF_DISJOINT`, by the `__PROLIFSTRIP`/`__BOTHSTRIP` estimators, and in TCGA by purity + leukocyte fraction on top |
| 3 | ...and it is OXPHOS **rather than** MYC | **SUPPORTED, and this is the strong result.** The two axes correlate (0.26-0.32 adjusted) and rank the 44 genes at 0.61-0.73, so the separate comparison could not have settled this. Conditioning OXPHOS on MYC leaves the localisation split intact (0.453 -> 0.485 TCGA, 0.489 -> 0.525 SCAN-B). Conditioning MYC on OXPHOS takes it to **zero** (0.187 -> -0.043, 0.137 -> -0.058). MYC's ordering was inherited; OXPHOS's was not |
| 4 | ...and OXPHOS therefore engages the death programme | **NOT SUPPORTED.** The machinery's mean \|rho\| with OXPHOS is not above an expression-matched null (z 1.3-1.5), its spread is not either (z_sd ~ 0), and the localisation split is ~1.3 SD above a null matched for expression AND sub-mitochondrial compartment - the 90th percentile, consistent across four cells, not separable |
| 5 | "driven" | **NOT TESTED.** Cross-sectional bulk transcriptomes. And CLAUDE.md trap 1: the sibling pre-registered study found the MYC x OXPHOS *interaction* on apoptotic priming null. Never write a causal verb here |

## What the ordering actually is

Not pro-death versus pro-survival. **Where the protein sits.** A gene's position
on the OXPHOS axis is predicted by MitoCarta membership at 0.45-0.49
(independent of the death curation, agreeing with it on 43 of 44 genes) and by
its annotated direction of effect at only 0.22-0.28. Mitochondrial members are
positive (median +0.20 to +0.25), non-mitochondrial ones negative (-0.16 to
-0.19), and the pro-death/pro-survival distinction cuts across both.

`APAF1` is the clean illustration: intrinsic-pathway by Reactome, absent from
MitoCarta because the apoptosome is cytosolic, and the most negative of the 44.
Pathway membership gets it wrong; localisation gets it right.

## The sentence I would write

> In two independent human breast cancer cohorts (TCGA-BRCA n=1,095, SCAN-B
> n=3,207), transcript levels of the canonical apoptotic machinery are ordered
> along tumour OXPHOS status and not along MYC activity. The two axes are
> themselves correlated, but conditioning OXPHOS on MYC leaves the ordering
> intact while conditioning MYC on OXPHOS abolishes it, and neither is explained
> by proliferation, tumour purity or immune infiltrate. What predicts a gene's
> position is where its protein acts - mitochondrial members rise with OXPHOS,
> cytosolic ones fall - and not whether it promotes or prevents death. The
> ordering is, however, no steeper than that of any gene set matched for
> expression and sub-mitochondrial compartment, indicating that the bulk
> transcriptome reports mitochondrial content rather than a selectively engaged
> death programme.

Four sentences. The fourth is what makes the first three publishable rather than
over-read, and it sets up whatever the mouse arm can show functionally: **the
human bulk transcriptome cannot see a death programme; it sees organelle
content.** That is a real limit on what transcriptomic evidence can be asked to
do, and it is worth stating as a finding rather than conceding as a weakness.

## The figure

**`E11_fig9_paper_figure1` is it, composed and ready.** Two panels, and the
second is not optional.

- **(A)** the 44 genes on the MYC x OXPHOS plane, square axes, identically
  scaled, coloured by MitoCarta membership, **proliferation-adjusted column
  only, both cohorts**. One image carries the descriptive claim: the cloud is
  taller than it is wide and the colour separates vertically, not horizontally.
  The SDs and their bootstrap interval are printed in the panel.
- **(B)** the conditioning ladder. Four rows, two cohorts: OXPHOS survives
  conditioning on MYC, MYC collapses to zero conditioned on OXPHOS, grey band is
  the compartment-matched null. This is what turns (A) from a description into a
  test, and the band is what stops a reviewer reading more into it than is
  there. The conditioned split gap and its interval are printed in the panel.

Provenance of the two panels, for anyone wanting the exploratory versions:
(A) is the adjusted column of `E11_fig1` - the raw column is a methods point and
belongs in supplementary; (B) is `E11_fig8` with the intervals added.
**`E11_fig2` (the mitoribosome control) belongs in supplementary and must be
cited in (B)'s legend** - it is what rules out "the adjustment emptied the MYC
score". The on-figure title should be deleted for submission; most journals set
it from the legend.

## The two versions of "mitochondrial genes correlate more"

The quantity panel B plots was, until 2026-09-02, the Spearman of the 44
per-gene correlations against the 0/1 MitoCarta label - a rank-biserial
correlation. It is a defensible statistic and it is a terrible axis label: it
describes a recipe rather than a quantity anyone thinks in.

The figures now plot **delta-rho** instead - the median correlation of the
mitochondrial genes minus the median of the cytosolic ones. Same units as the
correlations themselves, so it reads straight off panel A, and the axis becomes
`delta-rho (median rho of mitochondrial minus cytosolic genes)`, set in
plotmath so the symbol renders. Naming both groups in the label is deliberate:
`delta-rho` alone would be read as a change BETWEEN the rows of the panel, which
are four different models.

| | rank split | z | median difference | z |
|---|---|---|---|---|
| OXPHOS | 0.453 | 1.35 | 0.414 | 2.15 |
| OXPHOS, MYC removed | 0.485 | 1.63 | 0.376 | 2.14 |
| MYC | 0.187 | -0.25 | 0.106 | 0.02 |
| MYC, OXPHOS removed | -0.043 | -1.13 | -0.019 | -1.11 |

(TCGA; SCAN-B is the same pattern with z 2.49 and 2.35 on the OXPHOS rows.)

**THE CHOICE WAS MADE AFTER BOTH HAD BEEN SEEN, AND THE ONE NOW PLOTTED HAS THE
LARGER z.** That is the shape of statistic-shopping, so the claim is NOT moved
onto it: **the composition bound stays anchored on the rank split**, where the
OXPHOS rows sit ~1.4-1.6 SD above the null and the honest word remains "not
separable". Both are carried in every saved table, and if the median difference
is ever quoted, this paragraph has to be quoted with it.

The qualitative conclusion is identical under both, which is the reason the
switch is safe to make at all.

## Does it need a statistical test?

**The permutation nulls already ARE the test, and they are the right one.** 2,000
expression-matched draws, and for the localisation split a draw matched on
sub-mitochondrial compartment as well. Foreground the z and the percentile;
do not add p-values on top of them.

**The interval that was missing is now in `E11` section 4.4.** 1,000
tumour-level bootstrap resamples per cohort - resampling TUMOURS and not genes,
because the 44 are co-expressed and a gene-level bootstrap would treat 44
correlated observations as 44 independent ones and return an interval far too
narrow.

| contrast | tested against | TCGA | SCAN-B |
|---|---|---|---|
| SD ratio, OXPHOS / MYC | 1 | **1.64 [1.47, 1.78]** | **1.56 [1.46, 1.65]** |
| rank split gap | 0 | 0.27 [0.21, 0.33] | 0.35 [0.30, 0.40] |
| rank split gap, mutually conditioned | 0 | 0.53 [0.44, 0.58] | 0.58 [0.52, 0.64] |
| **delta-rho gap** | 0 | 0.31 [0.25, 0.37] | 0.29 [0.24, 0.33] |
| **delta-rho gap, mutually conditioned** | 0 | **0.40 [0.34, 0.47]** | **0.31 [0.28, 0.35]** |

All ten exclude their null value, and the conditioned gap is larger than the
unconditioned one on both statistics - conditioning does not weaken the
contrast, it sharpens it. **The delta-rho rows are the ones printed on the
figure**, because an interval must be in the units of the axis it sits beside;
the rank rows are what the composition bound is anchored on. Two independent
cohorts agreeing remains the stronger evidence; these are what a reviewer will
ask for.

**Two things NOT to do.**

- Do not put per-gene p-values or an FDR across the 44 anywhere. That is the
  grid-of-cells trap and it invites exactly the gene-picking the study has
  avoided. The per-gene values carry Fisher-z intervals already; that is enough
  to describe a gene and not enough to select one.
- Do not report the composition null as "not significant" and move on. Report it
  as the bound it is - the ordering is real and is OXPHOS's, and it is not
  specific to apoptosis. Both halves are the finding.

## Open, if the claim is to be pushed further

1. The ~1.3 SD residue above the compartment-matched null is the only thing
   standing between "OXPHOS orders these genes" and "OXPHOS orders these genes
   *because they are apoptotic*". It would be settled by a curated non-apoptotic
   mitochondrial comparator set rather than a random draw.
2. `LumA` alone has not been run; `Basal` rests on 171 TCGA samples.
3. Nothing here is functional. The ordering is a correlation in bulk tumour
   RNA and every causal word must come from elsewhere.
