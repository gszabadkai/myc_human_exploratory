# Phase 1 findings - the correlation atlas

2026-08-31. Source: `results/correlation_atlas.rds`, built by `E03` from `E02`.

> **EVERYTHING BELOW IS HYPOTHESIS-GENERATING.** Nothing in this repo is
> pre-registered. The atlas is 68,255 cells and no single one of them is a
> result. What follows is claimed only where it holds across **both cohorts**,
> **multiple instruments** and **the estimator panel** - and each claim carries
> the observation that would kill it. Those falsifiers are written here *before*
> `E04` and `E05` are built, so the exploratory phase can hand a real hypothesis
> to a confirmatory one.

The atlas reproduces the validation study's three published anchors to within
2-3e-4 when they are read back out of the finished table, which is what
establishes that it is indexed correctly and not merely computed correctly.

---

## F1. The MYC-OXPHOS correlation is real, replicated, and not merely proliferation

**What holds.** Against `OXPHOS subunits`, every one of the 18 MYC signatures is
positive in both cohorts on all four instruments. Raw GSVA spans 0.049 to 0.670
(SCAN-B) and 0.181 to 0.652 (TCGA). Across the 18 arms the two cohorts agree at
Spearman 0.92 (GSVA) and 0.94 (zmean).

It is **not** proliferation wearing MYC's name:

| | SCAN-B | TCGA |
|---|---|---|
| signatures with prolif-adjusted rho >= 0.2 | 15 / 18 | 17 / 18 |
| median prolif-adjusted rho | 0.323 | 0.452 |
| `MYC_UP.V1_UP` (1.5% entangled), raw | 0.502 | 0.552 |
| `MYC_UP.V1_UP`, proliferation-adjusted | 0.401 | 0.515 |
| `MYC_UP.V1_UP`, purity + leukocyte-adjusted | - | 0.538 |

The least proliferation-entangled signature in the compendium gives one of the
**highest** correlations, and loses only 0.04-0.10 when proliferation is
projected out. In TCGA the purity and leukocyte adjustment costs almost nothing
(FELSHER_61 0.388 -> 0.376).

**The honest complication, which must be reported with it.**
`slope_vs_entanglement` - the Spearman, across the 18 signatures, between a
signature's proliferation content and the rho it produces - is **positive in all
20 cohort x instrument x adjustment combinations**, +0.18 to +0.54. On average,
more entangled signatures do correlate more.

But the relationship is **not monotone and the extremes contradict it**:
`rho_least` (1.5% entangled) exceeds `rho_most` (47.6%) in **7 of 8** raw
combinations. The positive slope is carried by two weak signatures in the
*low*-entanglement half - `ELLWOOD_MYC_TARGETS_UP` (13 genes, 0.049 / 0.181) and
`ALFANO_MYC_TARGETS` (0.160 / 0.320) - not by a gradient. And it is not set
size: rho-vs-size is 0.05-0.18 where rho-vs-entanglement is 0.30-0.54.

So: **between-signature spread (0.05 to 0.70) is far larger than the
entanglement trend, and entanglement does not explain the correlation away.**
CLAUDE.md trap 3 stands - never quote one signature - but the panel as a whole
survives.

**What would falsify F1**
1. A proliferation estimator better than `PROLIF_DISJOINT` (e.g. a cell-cycle
   phase deconvolution, or Ki67 where measured) that removes the correlation for
   `MYC_UP.V1_UP` as well as for the entangled sets.
2. `MYC_UP.V1_UP` turning out to be entangled with something else that tracks
   OXPHOS - it is a perturbation-derived set and its 198 genes have not been
   audited against, say, mitochondrial mass or ribosome biogenesis.
3. The correlation vanishing in a cohort scored on a different platform
   (METABRIC microarray would be the cheap test, and is out of scope for
   phase 1 by decision, not by oversight).

---

## F2. MYC mRNA carries no OXPHOS signal at all. The gap is the result.

`log2(MYC)` against `OXPHOS subunits`: **-0.011** (SCAN-B), **-0.032** (TCGA).
Against MYC *activity* on the same arm in the same samples: **+0.427** and
**+0.389**.

Across all 18 arms, `log2(MYC)` never exceeds |rho| = 0.12 in SCAN-B or 0.21 in
TCGA, and its largest value is on `Folate and 1-C`, not on any OXPHOS arm. The
null is not specific to OXPHOS; the transcript simply does not track the
mitochondrial programme.

This replicates CLAUDE.md trap 4 in a second cohort. **Anyone plotting MYC
expression against OXPHOS sees nothing, and would conclude the opposite of what
the activity signatures say.**

**What would falsify F2**: MYC protein or a MYC-amplification-stratified
analysis showing that the mRNA is informative once its dynamic range is
accounted for. Partial support against this already exists in the atlas -
`M_c_call`, the copy-number call, gives +0.137 on `OXPHOS subunits` where the
mRNA gives -0.032, so the *locus* carries some signal that the *transcript level*
does not.

---

## F3. The two genomes are regulated differently - the sharpest finding, and it replicates to two decimal places

Nuclear-encoded and mtDNA-encoded OXPHOS do not behave alike, and the split is
almost identical in the two cohorts:

| FELSHER_61, raw, stratum all | content | gsva | mitopps | zmean |
|---|---|---|---|---|
| `OXPHOS subunits` SCAN-B / TCGA | 0.433 / 0.434 | 0.427 / 0.388 | 0.305 / 0.306 | 0.425 / 0.371 |
| `mtDNA-encoded OXPHOS` SCAN-B / TCGA | 0.043 / 0.056 | 0.061 / 0.068 | **-0.081 / -0.088** | 0.026 / 0.028 |

On the composition instrument the sign **flips**. Across 19 estimators, mitoPPS
returns a negative mtDNA correlation in 16/19 (SCAN-B) and 13/19 (TCGA); the
other three instruments return essentially none. Read as mitoPPS is built to be
read (trap 6): **MYC-high tumours de-prioritise the mtDNA-encoded programme
relative to the rest of the mitochondrion**, while raising the nuclear-encoded
one.

**And the 13 genes do not move together.** Against FELSHER_61, raw:

```
MT-CO2  +0.240 / +0.243      <- positive under 20/20 and 21/21 estimators
MT-ND2  +0.080 / +0.124
MT-ND4  +0.059 / +0.101
...
MT-ND1  -0.041 / -0.057
MT-CO1  -0.007 / -0.073
MT-ND6  -0.045 / -0.154
MT-ND5  -0.160 / -0.143      <- negative under 19/20 and 20/21 estimators
```

`MT-CO2` agrees between cohorts to **0.003**. `MT-ND5` to 0.017.

**Two confounders tested and cleared.**

- *Purity and infiltrate.* In TCGA the split **survives** purity + leukocyte
  adjustment: `MT-CO2` 0.243 -> 0.230, `MT-ND5` -0.143 -> -0.162, `MT-ND6`
  -0.154 -> -0.160. Not stromal or adipose contamination.
- *3' bias and polyadenylation.* There is **no positional gradient**. Ordered
  along the heavy strand, the maximum is at position 4 (`MT-CO2`) and position 1
  (`MT-ND1`) is already negative. `MT-ND6` - the only light-strand, only
  non-polyadenylated protein gene, and therefore the obvious library-chemistry
  suspect - is *not* the extreme in SCAN-B (-0.045).
- *Not tested, and it matters:* proliferation adjustment removes `MT-ND6`'s
  negative almost entirely (-0.154 -> -0.008 in TCGA) but leaves `MT-ND5`
  (-0.143 -> -0.103) and `MT-CO2` (+0.243 -> +0.166). So the ND6 signal and the
  ND5/CO2 signal are not the same phenomenon.

**Why this is interesting rather than merely odd.** `MT-CO2` and `MT-ND5` are
both on the **heavy-strand polycistron** - one primary transcript. A stable,
cross-cohort, opposite-signed correlation between two genes of the same
transcript cannot be explained by mtDNA copy number or by heavy-strand
transcription rate. It has to be post-transcriptional, or an artefact.

**What would falsify F3**
1. **Copy number.** If MYC-high tumours simply carry fewer mitochondrial
   genomes, all 13 should move together. They do not - `MT-CO2` is +0.24 while
   `MT-ND5` is -0.15 - so uniform copy number is already excluded. The remaining
   test is a per-sample mtDNA-content proxy as a covariate: if adjusting for it
   collapses the *spread* among the 13, the finding is a copy-number artefact
   with gene-specific quantification noise on top.
2. **Quantification.** `MT-ND5` and `MT-ND6` overlap on opposite strands, so
   multi-mapping and strandedness could produce exactly this. Test: re-quantify
   the 13 from a stranded pipeline, or check whether the `MT-ND5`/`MT-ND6`
   correlation with each other is anomalous relative to the other 11.
3. **The mitoPPS sign flip is a composition artefact.** mitoPPS is a ratio
   against the rest of the mitochondrial programme. If the nuclear programme
   simply rises more, the mtDNA arm must fall in composition terms with no
   biology at all. **This is the most likely benign explanation and it is not
   yet excluded.** Test: hold the mitoPPS denominator fixed - score the mtDNA
   arm against a mitochondrial universe with the OXPHOS subunits removed - and
   see whether the negative survives.
4. Failure to reproduce in a third cohort.

Falsifier 3 is the one to run first. It is cheap and it is the difference
between a finding and an arithmetic identity.

---

## F3 VERDICT - `E03b` run, 2026-08-31

The original claim above is left standing as written. This is what the
falsifiers returned, recorded beside it rather than folded into it. **Two parts
of F3 come out stronger than they went in, one comes out weaker, and one gene
has to be dropped.**

| Falsifier | Criterion | Result |
|---|---|---|
| 3a mitoPPS denominator, no OXPHOS gene in it | mtDNA arm rho still < 0 in both cohorts | SURVIVES |
| 3b instrument-free ratio vs a non-OXPHOS mitochondrion | log2 ratio rho still < 0 in both cohorts | SURVIVES |
| 1 mtDNA copy number as a covariate | MT-CO2 still > 0 and MT-ND5 still < 0 in both | SURVIVES |
| 2 quantification | judgement, not pass/fail | **MT-ND6 fails it** |

### 3a WEAKENED, and the mitoPPS framing should be demoted

Taking OXPHOS out of the denominator does move the number, by about a fifth:

```
mtDNA arm, FELSHER_61      full     no_oxname   no_oxgene
TCGA                      -0.088     -0.078      -0.071
SCAN-B                    -0.081     -0.073      -0.065
```

and across all 20 estimators the median goes from -0.024 to **-0.004** in TCGA
and -0.039 to -0.016 in SCAN-B, with the count of negative estimators falling
14 -> 11 and 16 -> 14. So part of the mitoPPS negative genuinely *was* the
composition identity, and what remains at the arm level is small enough that
**the mitoPPS sign flip should no longer be quoted as the headline.** It is a
weak effect on a composition instrument, and section 3b says the same thing far
better.

### 3b SURVIVES, and reframes the finding into something sharper

Dropping mitoPPS for a log2 ratio anyone can recompute:

```
FELSHER_61 against          TCGA    SCAN-B    estimators negative (of 20)
mtDNA / all genes          +0.004   +0.023      4 / 5     median +0.078 / +0.058
mtDNA / rest of mitochondrion -0.122  -0.120   19 / 18    median -0.065 / -0.080
mtDNA / mitochondrion, no OXPHOS -0.110 -0.111 16 / 18    median -0.047 / -0.060
mtDNA / nuclear OXPHOS subunits  -0.189 -0.163 20 / 19    median -0.143 / -0.154
mtDNA / mitochondrial ribosome   -0.205 -0.209 20 / 19    median -0.145 / -0.165
```

**Against the whole transcriptome the mtDNA programme is flat to slightly
positive. Against the mitochondrion it is negative. Against the mitoribosome it
is most negative of all, in 20 of 20 and 19 of 20 estimators.**

That changes what should be said. The right statement is **not** that MYC
suppresses the mtDNA-encoded programme - the `content` instrument (+0.056 /
+0.043) and the whole-transcriptome ratio both say it is flat. The right
statement is:

> **MYC-high tumours build more nuclear-encoded OXPHOS and more mitochondrial
> ribosome, and do not raise the mtDNA-encoded transcripts to match.** The
> mtDNA programme is not pushed down; everything around it goes up and it stays
> where it is.

That also resolves the F6 tension the other way round from how it was posed.
The mitoribosome is not merely the top arm - it is the arm the mtDNA programme
falls furthest behind. More machinery for translating mtDNA transcripts, no
more transcripts.

Note what this contrast is and is not: `mtDNA / mitoribosome` is negative
largely *because* the mitoribosome is strongly MYC-positive. That is the claim
stated as a contrast, not independent evidence for it. The independent part is
`vs_all_genes` being ~0, which is what licenses "flat" rather than "suppressed".

### 1 SURVIVES, and inverts: copy number was HIDING the split

Conditioning on mtDNA content does not collapse the per-gene spread. It nearly
doubles it:

```
                         TCGA raw -> content-adj   SCAN-B raw -> content-adj
MT-CO2                    +0.243  ->  +0.440        +0.240  ->  +0.413
MT-ND2                    +0.124  ->  +0.178        +0.080  ->  +0.099
MT-ND1                    -0.057  ->  -0.198        -0.041  ->  -0.131
MT-CO1                    -0.073  ->  -0.295        -0.007  ->  -0.113
MT-ND5                    -0.143  ->  -0.266        -0.160  ->  -0.337

median spread over the 13   0.408 ->   0.717         0.398  ->   0.744
```

The shared component - which is what a copy-number difference would be - was
**masking** the within-mtDNA structure, not creating it. This is the opposite of
what the copy-number explanation predicts, and it is the strongest single result
in E03b. What is left after conditioning is within-mtDNA regulation: `MT-CO2`
up, `MT-CO1`, `MT-ND1`, `MT-ND5` down, in both cohorts.

### 2 FAILS FOR MT-ND6, WHICH IS THEREFORE DROPPED FROM THE CLAIM

The 13 x 13 structure singles out exactly one gene, and it is not the one the
original note leaned on hardest:

```
                    TCGA     SCAN-B
median pairwise rho  0.742    0.751
MT-ND5 with others   0.651    0.680     <- normal in both
MT-CO2 with others   0.731    0.691     <- normal in both
MT-ND6 with others   0.576    0.219     <- NOT normal in SCAN-B
MT-ND5 - MT-ND6      0.862    0.354
```

In SCAN-B `MT-ND6` correlates 0.12 to 0.36 with every other mtDNA gene while
every other pair sits at 0.6 to 0.94. It is effectively **decoupled from the rest
of the mitochondrial genome in one cohort and not the other**. That is what a
technical artefact looks like, and `MT-ND6` is the obvious candidate: the only
light-strand, only non-polyadenylated protein gene, in two cohorts with
different library chemistry. Its MYC correlation is correspondingly the worst
replicator of the 13 (-0.154 vs -0.045 raw; -0.251 vs -0.056 adjusted).

**`MT-ND6` is dropped from F3.** It is reported as a measured artefact, not as
evidence. `MT-ND5` is unaffected - its correlation structure is normal in both
cohorts and its MYC correlation replicates and strengthens under conditioning.

### What F3 now says, and what would still falsify it

1. The mtDNA-encoded programme is **flat** with MYC activity in absolute terms,
   and **falls behind** the nuclear-encoded mitochondrial programme, most
   sharply behind the mitochondrial ribosome. Replicated in two cohorts across
   19-20 of 20 MYC estimators.
2. Within the 13 genes there is real regulation that is **not** copy number and
   is exposed rather than hidden by conditioning on it: `MT-CO2` up, `MT-CO1`,
   `MT-ND1`, `MT-ND5` down. `MT-CO2` and `MT-CO1` are adjacent on the same
   transcript and go opposite ways, which remains the most striking and the
   least explained part of it. Their mutual correlation (0.77 / 0.78) is
   unremarkable, so it is not obviously a mismapping artefact - but that is an
   absence of evidence, not a clean test.
3. The mitoPPS sign flip is demoted to a footnote.

**Still open, and now the sharpest tests:**
- Re-quantify the 13 from a stranded, total-RNA (non-polyA) dataset. That kills
  the `MT-ND6` question outright and tests whether `MT-CO1` / `MT-CO2` survive.
- A third cohort, ideally on a different platform.
- Direct mtDNA copy number (WGS-derived) rather than the expression proxy used
  here.
- Protein-level or ribosome-profiling data would separate "not transcribed" from
  "not translated", which is the mechanistic fork this cannot resolve.

---

## F4. There are three instruments here, not four

`gsva` and `zmean` agree at median Spearman **0.994** across the 18 arms
(minimum 0.953 over 21 estimators). They are not independent measurements.

Every other pair sits at 0.87-0.90 median but with minima near 0.50, so
disagreement is estimator-specific rather than uniform. `mitopps` is the
instrument that replicates *least* across cohorts (239 of 360 arm-cells
"replicated" against `content`'s 287) and it is also the only one that changes
the sign of the mtDNA arm - i.e. exactly where it disagrees is where it is
designed to say something the others cannot.

CLAUDE.md trap 5 quoted 0.24-0.94 for GSVA-vs-mitoPPS agreement of the **scores**.
This is a different quantity - agreement of the **answers**, arm ordering against
the same estimator - and it is higher. Both belong on the figure and they must
not be conflated.

---

## F5. Present in every subtype, weakest in ER-negative

`FELSHER_61` / `OXPHOS subunits` / GSVA / raw, by stratum:

```
           SCAN-B   TCGA
all         0.427  0.388
ERpos       0.447  0.408
ERneg       0.242  0.337
LumA        0.368  0.418
LumB        0.384  0.272
HER2        0.359  0.290   (TCGA n = 78, CI 0.067 to 0.485)
Basal       0.381  0.443
Normal      0.391  0.273   (TCGA n = 36, CI -0.068 to 0.556)
```

Median across the whole 18-signature panel: ERpos 0.479 / 0.512, ERneg 0.315 /
0.416. **The correlation is not carried by one subtype**, which is the main
thing worth knowing, and ER-negative is consistently the weakest stratum in both
cohorts. TCGA's HER2 and Normal strata are too small to interpret and are shown
with their intervals for that reason.

---

## F6. The mitochondrial ribosome beats OXPHOS, in both cohorts

`Mitochondrial ribosome` 0.617 / 0.590 against `OXPHOS subunits` 0.427 / 0.388,
and it is the top arm of 18 in both cohorts. Consistent with MYC's canonical
ribosome-biogenesis role extending to the mitochondrial ribosome - and it sits
awkwardly beside F3, because the mitoribosome exists to translate the very
mtDNA-encoded transcripts whose programme is de-prioritised. **That tension is
the most interesting thing in the atlas and phase 2 should go at it.**

## F7. Fatty acid oxidation is the only consistently negative arm

-0.140 / -0.121, and the only arm whose entanglement slope is negative (-0.31
TCGA, -0.14 SCAN-B): the *more* proliferation-entangled the signature, the more
negative the FAO correlation. Read with trap 2 - adipose is FAO-high and breast
is the worst TCGA tissue for infiltrate - but note it replicates in SCAN-B,
where there is no purity estimate to blame it on.

---

## What is NOT claimed

- **No interaction is claimed anywhere.** The validation study found the
  `MYC x OXPHOS` interaction on apoptotic priming null. Everything above is
  correlation. CLAUDE.md trap 1: these are different questions and neither
  confirms nor contradicts the other.
- **No mitoPPS value is compared across cohorts** - only patterns (trap 6).
- **No p-value appears above.** The atlas carries one per cell so a cell can be
  described; with 68,255 cells it cannot select one.
- Nothing about the death axis. That is `E05`.
