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

## The priming subsection: how to write it descriptively

Decided 2026-09-02: the BCL2-family results stay **descriptive**. The
interaction hypothesis was tested and did not survive
(`docs/2026-09-02_priming_interaction_tested.md`), and priming is a
post-translational property that transcript abundance cannot measure. What the
data support is a description of how these transcripts sit on the two axes.

### RESOLVED 2026-09-02: E10 IS NOW PROLIFERATION-ADJUSTED

Every value E10 plots is a partial Spearman on `PROLIF_DISJOINT`, the same
covariate E11 and E13 use, so all three sit on the same footing. The unadjusted
pass is still computed and saved as `gene_cor_raw`, because section 4's
reproduction check is against E08, which has no covariate to reproduce - it
still returns `max |difference| = 0e+00`.

**The adjustment changed two numbers that matter and both got cleaner**, so the
text below already uses the adjusted values: the count of ratios beating their
components went from 6 to **0** on the OXPHOS axis and from 3 to 5 on MYC.

### The three things worth saying, in order

**1. The family does not move as a block, and functional class does not predict
position.** This is the observation everything else rests on, and it is the one
that most needs stating before any ratio is shown.

> Across both cohorts the twelve BCL2-family transcripts span the OXPHOS axis
> from `BAD` (+0.50) to `MCL1` (-0.27) after adjustment for proliferation, and
> their position is not predicted by whether they promote or prevent apoptosis:
> the two most positive are `BAD` (pro-apoptotic) and `BCL2L1`
> (anti-apoptotic), and the two most negative are `MCL1` (anti-apoptotic) and
> `BMF` (pro-apoptotic). The same transcripts occupy a markedly narrower range
> on the MYC axis (+0.30 to -0.21).

**2. The ratio matrix is almost perfectly additive in its two components.**
This replaces an earlier draft that called it "organised by its denominators" -
that was checked after adjustment and is not right. The quantitative version is
better anyway, because it puts a number on how little the pairing adds.

> A pairwise ratio inherits the correlations of both its members and adds almost
> nothing to them. An additive model of numerator and denominator identity
> accounts for 92 to 95% of the variance across the 35 ratios, in both cohorts
> and on both axes, with the pro-apoptotic numerator contributing about twice as
> much as the anti-apoptotic denominator (R-squared 0.58 to 0.74 against 0.20 to
> 0.35). Once proliferation is accounted for, **none** of the 35 ratios on the
> OXPHOS axis exceeds both of its component transcripts in both cohorts, and 5
> of 35 do so on the MYC axis. We therefore present the ratios as a compact
> display of the component correlations, and not as a measurement of apoptotic
> priming - a post-translational property that transcript abundance cannot
> address.

That 92-95% is the sentence to keep if only one survives editing. It says the
matrix carries no pair-specific information, which is the whole reason the
subsection is descriptive.

**3. The two axes differ in the same direction as the wider machinery.** This is
what ties the subsection back to the paper's spine.

> The contrast between the two axes reproduces that seen across the wider
> apoptotic machinery: the transcripts spread further along OXPHOS than along
> MYC, and conditioning each axis on the other removes the MYC association while
> leaving the OXPHOS one intact. Two transcripts are exceptions and retain a
> MYC association after OXPHOS is accounted for, in both cohorts and under all
> three MYC estimators: `BID` (+0.12 to +0.28) and `PMAIP1`/NOXA (-0.05 to
> -0.21). Both are pro-apoptotic BH3-only proteins and they move in opposite
> directions, so this is not a coherent shift in the balance.

### What NOT to write

| | |
|---|---|
| "the priming ratio increases with OXPHOS" | The ratio is its two components added. Say which transcript moves |
| "MYC represses BCL2, shifting the balance towards death" | `BCL2` is -0.22 on the reference estimator and -0.06 and -0.01 on the other two after conditioning. Estimator-dependent; do not write it |
| "OXPHOS-high tumours are more primed" | Priming is not measurable here. Write "carry a higher `BAD`/`MCL1` transcript ratio" and let the reader draw the inference |
| any interaction language | Tested, failed three falsifiers, and the sibling pre-registered study found the functional version null |

### The figure legend, in one sentence each

- **individual transcripts** (`E13` panel A, or `E10_fig5`): "The twelve
  BCL2-family transcripts ranked by their correlation with OXPHOS, showing that
  pro- and anti-apoptotic members are interleaved across the whole range."
- **the ratio matrix** (`E10_fig3`): "All 35 pro-over-anti transcript ratios
  against each axis, adjusted for proliferation. An additive model of the two
  component transcripts accounts for 92-95% of these values; see Methods."
- **the compartment split** (`E10_fig6`): supplementary. It is the same matrix
  inside luminal and basal tumours and its main use is to show that the pooled
  MYC column is a between-subtype effect.

## Open, if the claim is to be pushed further

1. The ~1.3 SD residue above the compartment-matched null is the only thing
   standing between "OXPHOS orders these genes" and "OXPHOS orders these genes
   *because they are apoptotic*". It would be settled by a curated non-apoptotic
   mitochondrial comparator set rather than a random draw.
2. `LumA` alone has not been run; `Basal` rests on 171 TCGA samples.
3. Nothing here is functional. The ordering is a correlation in bulk tumour
   RNA and every causal word must come from elsewhere.
