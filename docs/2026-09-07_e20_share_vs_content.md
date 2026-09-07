# E20 - does MYC track compartment SHARE or compartment CONTENT?

**2026-09-07. EXPLORATORY. Nothing here is pre-registered.**

Script: `scripts/E20_share_vs_content.R`. **Sourced by the author 2026-09-07**,
and `results/share_vs_content.rds` plus
`outputs/tables/E20_share_vs_content.csv` are on disk. The numbers below were
written from a dry run beforehand; the real run reproduces **all 14
non-timestamp saved objects identically** - the seeded bootstrap included - so
every number here is now read from a saved object and **the dry-run caveat is
discharged**.

**THIS IS A RE-READ, NOT A NEW ANALYSIS.** E16's separability table was produced
to answer one question - which respiratory ruler is least MYC-entangled - and it
is read here for a different one. Nothing is re-scored: mitoPPS and GSVA are
read from their saved objects, and the two composites E16 built are rebuilt to
E16's own recipe and then **checked against E16's saved numbers**. **48 of 48
cells reproduce, max |delta| = 0 exactly** (section 3). One quantity is new -
`ox_nuc_mtrib` - and it is the only new computation in the script.

**N3 throughout.** Every number here is a transcript correlation. The word
"primed" appears nowhere as a description of a tumour, a cell or an animal.
**The mouse is not read, compared or cited as a value by this script.**

---

## 0. The answer in six lines

| | |
|---|---|
| **Does MYC track content or share?** | **ESTIMATOR-DEPENDENT, and that is the result.** Content wins on all three signature estimators in both cohorts; the **CollecTRI regulon dissents in both cohorts**. 6 of 8 adjudicating cells say content; the 2 regulon cells say **neither** - share against `ox_rel`, content against `ox_nuc_mtrib`. |
| **Is the dissent noise?** | **No. It replicates**, and it is E16 section 3.1's pattern reappearing on a new question - because it *is* E16's table, asked differently. |
| **Against the mitoribosome specifically?** | **One-directional, 8 of 8.** `rho(ox_lvl) - rho(ox_nuc_mtrib)` is positive in every adjudicating cell, both cohorts, raw and adjusted - **including the regulon**. |
| **Why the two share rulers disagree** | The three signature estimators track the **mitoribosome** more than OXPHOS (`mtrib_lvl - ox_lvl` positive, 6 of 6). The regulon is the one estimator that does **not** (+0.004 / +0.003, both intervals covering 0). |
| **"Overwhelmed, not undone"?** | **NOT AVAILABLE** on the rule fixed before the panel was computed. Two of its three conditions fail. Section 5 says exactly which, and what would settle it. |
| **What would settle it** | **`COLLECTRI_MYC_STIM`, 739 genes, in the snapshot and never scored.** The entire verdict hangs on one regulon estimator. E16 already named this the cheapest thing that could overturn its own section 3.1; it is now the cheapest thing that could settle E20 too. |

**The prompt's mouse premise - that development lowers OXPHOS *share* between 6W
and 12W while compartment *content* rises - is taken as given and is NOT checked
here.** The mouse repo was not read by this script. The human side of the
argument therefore stands or falls on its own.

---

## 1. What the panel is, and why the design is clean

**The three rulers are the same numerator with three denominators.** That is the
whole design, and it is why the contrast isolates the thing it claims to.

| ruler | definition | reads as |
|---|---|---|
| `ox_lvl` | `comp(OXPHOS subunits)` | **CONTENT** |
| `ox_rel` | `comp(OXPHOS subunits) - comp(rest of MitoCarta)` | **SHARE**, against the whole compartment |
| **`ox_nuc_mtrib`** | **`comp(OXPHOS subunits) - comp(mitoribosome)`** | **SHARE**, against the one denominator the reprioritisation claim is about. **NEW** |
| `mtrib_lvl` | `comp(mitoribosome)` | context: the new denominator alone |
| `ox_gsva`, `ox_ppd` | read, not rebuilt | context: the incumbent and mitoPPS |

`comp()` is E16's `.comp`, copied verbatim: the mean per sample of each gene's
z-score across samples, on `log2(linear DESeq2-normalised + 1)`. That is the
mouse's `comp_e`.

**One departure from the mouse, named as a decision.** The mouse's helper is
`relify(s) = comp(s) - comp(mito_all \ s)` - a set against the *rest of the
compartment*. `ox_nuc_mtrib` is **not** `relify("MITOCARTA_MITOCHONDRIAL_
RIBOSOME")`, which is a different object. What is matched to the mouse is the
**composite recipe**; the denominator is chosen to state the reprioritisation
axis directly, OXPHOS against the mitoribosome, rather than each against the
whole. Both readings are legitimate and the script carries `ox_rel` beside it so
the difference between them is a number rather than an argument.

### What the new ruler is built from

| | nominal | TCGA used | SCAN-B used | coverage |
|---|---|---|---|---|
| numerator `OXPHOS subunits` | 89, **0 mtDNA-encoded** | 88 | 88 | 0.989 / 0.989 |
| denominator `Mitochondrial ribosome` | 83, **0 mtDNA-encoded** | 82 | 83 | 0.988 / 1.000 |

Asserted rather than described: neither half carries an `MT-` gene (trap 8 is
not engaged), the two halves are **disjoint**, and the mitoribosome sits inside
the same pinned MitoCarta inventory that is `ox_rel`'s denominator - so
`ox_nuc_mtrib` is a *narrowing* of `ox_rel`'s denominator, not a different
universe. SCAN-B is resolved through `scanb_scores.rds$symbol_map` throughout;
unharmonised the numerator covers 0.775 instead of 0.989 (trap 7).

### What the new ruler shares with what it is measured against

| estimator | n | x numerator | x rest-of-MitoCarta | **x mitoribosome** |
|---|---|---|---|---|
| `FELSHER__MITOSTRIP` | 61 | 0 | 0 | **0** |
| `MYC_UP.V1_UP__MITOSTRIP` | 178 | 0 | 0 | **0** |
| `HALLMARK_MYC_TARGETS_V1__MITOSTRIP` | 177 | 0 | 0 | **0** |
| `M_b__MITOSTRIP` (CollecTRI) | 811 | 0 | 0 | **0** |
| `PROLIF_DISJOINT` | 318 | 0 | 10 | **0** |

**The mitoribosome column is zero everywhere, including the covariate.** The new
ruler shares no gene with anything it is measured against, on either half.

### Is it actually a new ruler?

Spearman, **within cohort only** - a value comparison across cohorts is
meaningless here and across species is forbidden.

| pair | TCGA | SCAN-B |
|---|---|---|
| `ox_lvl` vs `ox_rel` | 0.941 | 0.924 |
| **`ox_lvl` vs `ox_nuc_mtrib`** | **0.402** | **0.408** |
| **`ox_rel` vs `ox_nuc_mtrib`** | **0.542** | **0.567** |
| `ox_lvl` vs `mtrib_lvl` | 0.905 | 0.865 |
| `ox_nuc_mtrib` vs `mtrib_lvl` | 0.024 | -0.058 |
| `ox_lvl` vs `ox_gsva` | 0.986 | 0.975 |

**Yes.** `ox_nuc_mtrib` agrees with `ox_rel` at 0.54 / 0.57, far below the 0.92 /
0.94 that makes `ox_lvl` and `ox_rel` near-neighbours, and it is all but
orthogonal to its own denominator. It is a third column, not a relabelling.

**And one caveat that follows from the same table.** `ox_lvl` and `mtrib_lvl`
correlate at 0.905 / 0.865, so `ox_nuc_mtrib` is a difference between two highly
correlated composites. Its variance is what *separates* two co-regulated
programmes, and nothing here establishes that what separates them is biological
rather than technical - transcript length, GC content and expression level all
differ systematically between MRPs and OXPHOS subunits. It correlates 0.40 with
content in both cohorts, so it is not noise; that it is *biology* is an
assumption, not a result.

---

## 2. The panel

Spearman. **Both raw and proliferation-adjusted are reported**: raw is the
literal re-read of E16's section 3 table, and the adjusted panel - partial
Spearman on `PROLIF_DISJOINT`, 318 genes, as E16 section 4 adjusts - is what the
verdict is scored on. The `prolif` column is raw only; partialling
`PROLIF_DISJOINT` out of a correlation *with* `PROLIF_DISJOINT` is degenerate
and is left `NA` rather than printed as a near-zero someone might read.

### RAW - and rows 1, 2, 5 and 6 are E16's table, bit-equal

**TCGA**

| ruler | myc_mRNA | myc_msigdb | myc_regulon | myc_felsher | myc_lowent | prolif |
|---|---|---|---|---|---|---|
| `ox_lvl` | -0.016 | 0.525 | 0.279 | 0.367 | 0.471 | 0.208 |
| `ox_rel` | -0.103 | 0.442 | 0.333 | 0.259 | 0.329 | 0.159 |
| **`ox_nuc_mtrib`** | **-0.333** | **-0.219** | **-0.120** | **-0.395** | **-0.284** | **-0.333** |
| `mtrib_lvl` | 0.132 | **0.676** | 0.354 | **0.583** | **0.657** | 0.371 |
| `ox_gsva` | -0.032 | 0.558 | 0.240 | 0.388 | 0.479 | 0.237 |
| `ox_ppd` | -0.014 | 0.497 | 0.292 | 0.306 | 0.396 | 0.240 |

**SCAN-B**

| ruler | myc_mRNA | myc_msigdb | myc_regulon | myc_felsher | myc_lowent | prolif |
|---|---|---|---|---|---|---|
| `ox_lvl` | -0.003 | 0.552 | 0.190 | 0.422 | 0.452 | 0.383 |
| `ox_rel` | -0.048 | 0.470 | 0.254 | 0.328 | 0.338 | 0.334 |
| **`ox_nuc_mtrib`** | **-0.096** | **-0.156** | **-0.078** | **-0.282** | **-0.192** | **-0.196** |
| `mtrib_lvl` | 0.045 | **0.693** | 0.250 | **0.623** | **0.612** | 0.521 |
| `ox_gsva` | -0.011 | 0.577 | 0.117 | 0.427 | 0.430 | 0.394 |
| `ox_ppd` | -0.025 | 0.484 | 0.175 | 0.305 | 0.338 | 0.416 |

Raw, `ox_nuc_mtrib` is **negative on every estimator in both cohorts**, 10 of 10
including `prolif`. Part of that is proliferation: the mitoribosome is more
proliferation-entangled than OXPHOS (0.371 / 0.521 against 0.208 / 0.383), so a
ruler that subtracts it partly measures "not proliferating". That is what the
adjusted panel is for, and it moves the cells a long way.

### PROLIFERATION-ADJUSTED - the panel the verdict is scored on

**TCGA**

| ruler | myc_mRNA | myc_msigdb | myc_regulon | myc_felsher | myc_lowent |
|---|---|---|---|---|---|
| `ox_lvl` | -0.062 | 0.646 | 0.207 | 0.318 | 0.435 |
| `ox_rel` | -0.141 | 0.562 | **0.296** | 0.210 | 0.293 |
| `ox_nuc_mtrib` | -0.286 | **+0.109** | 0.057 | -0.239 | -0.137 |
| `mtrib_lvl` | 0.061 | 0.711 | 0.210 | 0.489 | 0.586 |

**SCAN-B**

| ruler | myc_mRNA | myc_msigdb | myc_regulon | myc_felsher | myc_lowent |
|---|---|---|---|---|---|
| `ox_lvl` | -0.005 | 0.458 | 0.033 | 0.257 | 0.349 |
| `ox_rel` | -0.053 | 0.369 | **0.131** | 0.166 | 0.233 |
| `ox_nuc_mtrib` | -0.097 | +0.014 | 0.006 | -0.209 | -0.124 |
| `mtrib_lvl` | 0.050 | 0.549 | 0.036 | 0.450 | 0.510 |

The three panel rulers with intervals, adjusted, on the four adjudicating
estimators:

| cohort | estimator | `ox_lvl` | `ox_rel` | `ox_nuc_mtrib` |
|---|---|---|---|---|
| TCGA | `myc_msigdb` | +0.646 [+0.606, +0.682] | +0.562 [+0.517, +0.605] | **+0.109 [+0.050, +0.167]** |
| TCGA | `myc_regulon` | +0.207 [+0.149, +0.263] | +0.296 [+0.240, +0.350] | +0.057 [-0.002, +0.116] |
| TCGA | `myc_felsher` | +0.318 [+0.263, +0.372] | +0.210 [+0.152, +0.266] | -0.239 [-0.295, -0.182] |
| TCGA | `myc_lowent` | +0.435 [+0.383, +0.484] | +0.293 [+0.236, +0.347] | -0.137 [-0.195, -0.078] |
| SCAN-B | `myc_msigdb` | +0.458 [+0.429, +0.486] | +0.369 [+0.338, +0.399] | +0.014 [-0.020, +0.049] |
| SCAN-B | `myc_regulon` | +0.033 [-0.002, +0.067] | +0.131 [+0.097, +0.165] | +0.006 [-0.029, +0.040] |
| SCAN-B | `myc_felsher` | +0.257 [+0.224, +0.290] | +0.166 [+0.132, +0.200] | -0.209 [-0.242, -0.175] |
| SCAN-B | `myc_lowent` | +0.349 [+0.317, +0.380] | +0.233 [+0.199, +0.266] | -0.124 [-0.158, -0.090] |

### Which estimators adjudicate, and the rule that decided it

Fixed before the panel was computed: an estimator adjudicates only if
`|raw rho with ox_lvl| >= 0.15` **in both cohorts**. The difference between two
near-zero correlations has a direction and no meaning, and E16 section 3.1 made
exactly this argument in prose about `myc_mRNA` and declined to lean on it. Here
it is a number, fixed in advance rather than applied afterwards.

| estimator | min |rho| across cohorts | adjudicates |
|---|---|---|
| `myc_mRNA` | 0.003 | **no** |
| `myc_msigdb` | 0.525 | yes |
| `myc_regulon` | 0.190 | yes |
| `myc_felsher` | 0.367 | yes |
| `myc_lowent` | 0.452 | yes |

`myc_mRNA` is reported in every table and scored in none - trap 4, and E16's own
handling of it.

---

## 3. The contrast

`rho(est, ox_lvl) - rho(est, ox_rel)` and `- rho(est, ox_nuc_mtrib)`.

**The interval is a paired sample-level bootstrap, B = 2,000, seed 20260831.**
The three rulers share their numerator and correlate at 0.40-0.94, so the two
marginal Fisher-z intervals in section 2 say nothing useful about their
difference - they would be far too wide. The bootstrap resamples samples and
recomputes both rhos on the same resample, which is the paired quantity.
Percentile interval, **descriptive, not a test**. "Material" means the interval
excludes 0 **and** `|delta| >= 0.05`, the bar borrowed from E16's own
`R1_MAX_DELTA` rather than invented here.

### `ox_lvl - ox_rel` - content minus whole-compartment share, adjusted

| cohort | estimator | delta | bootstrap 95% | material |
|---|---|---|---|---|
| TCGA | `myc_msigdb` | +0.083 | [+0.065, +0.102] | yes |
| TCGA | **`myc_regulon`** | **-0.089** | **[-0.109, -0.069]** | **yes** |
| TCGA | `myc_felsher` | +0.109 | [+0.088, +0.132] | yes |
| TCGA | `myc_lowent` | +0.142 | [+0.122, +0.165] | yes |
| SCAN-B | `myc_msigdb` | +0.089 | [+0.075, +0.102] | yes |
| SCAN-B | **`myc_regulon`** | **-0.098** | **[-0.112, -0.084]** | **yes** |
| SCAN-B | `myc_felsher` | +0.091 | [+0.077, +0.104] | yes |
| SCAN-B | `myc_lowent` | +0.116 | [+0.102, +0.130] | yes |

**On this contrast, six cells say content and two say share - and the two are
the same estimator in both cohorts.** Every one of the eight is material, so this is not
a tie decided by noise - it is a genuine, replicating disagreement between
estimators. Raw gives the same eight signs.

### `ox_lvl - ox_nuc_mtrib` - content minus mitoribosome-relative share, adjusted

| cohort | estimator | delta | bootstrap 95% | material |
|---|---|---|---|---|
| TCGA | `myc_msigdb` | +0.537 | [+0.485, +0.586] | yes |
| TCGA | `myc_regulon` | +0.150 | [+0.092, +0.211] | yes |
| TCGA | `myc_felsher` | +0.558 | [+0.506, +0.611] | yes |
| TCGA | `myc_lowent` | +0.572 | [+0.517, +0.625] | yes |
| SCAN-B | `myc_msigdb` | +0.444 | [+0.415, +0.474] | yes |
| SCAN-B | `myc_regulon` | +0.027 | [-0.007, +0.063] | **no** |
| SCAN-B | `myc_felsher` | +0.466 | [+0.436, +0.496] | yes |
| SCAN-B | `myc_lowent` | +0.473 | [+0.442, +0.505] | yes |

**Positive in 8 of 8, and material in 7.** The one non-material cell is the
regulon in SCAN-B, and it is non-material because it is near zero, not because
it points the other way. Raw gives the same eight signs and all eight material.

**So the regulon's dissent is specific to the whole-compartment denominator.**
Against the mitoribosome, every estimator agrees that MYC tracks content.

### `mtrib_lvl - ox_lvl` - context, and the reason the two share rulers differ

Not a verdict input. `ox_nuc_mtrib` is a difference of two composites, so its
sign alone cannot say whether MYC is failing to raise OXPHOS or is raising the
mitoribosome harder. This is that decomposition. Adjusted:

| cohort | estimator | delta | bootstrap 95% | material |
|---|---|---|---|---|
| TCGA | `myc_msigdb` | +0.065 | [+0.047, +0.084] | yes |
| TCGA | **`myc_regulon`** | **+0.003** | **[-0.023, +0.027]** | **no** |
| TCGA | `myc_felsher` | +0.171 | [+0.149, +0.194] | yes |
| TCGA | `myc_lowent` | +0.151 | [+0.130, +0.175] | yes |
| SCAN-B | `myc_msigdb` | +0.091 | [+0.074, +0.109] | yes |
| SCAN-B | **`myc_regulon`** | **+0.004** | **[-0.015, +0.023]** | **no** |
| SCAN-B | `myc_felsher` | +0.193 | [+0.173, +0.212] | yes |
| SCAN-B | `myc_lowent` | +0.160 | [+0.141, +0.179] | yes |

**This is the mechanism of the dissent, and it replicates exactly.** The three
signature estimators put MYC nearer the **mitoribosome** than OXPHOS, in both
cohorts, materially. The CollecTRI regulon puts it **equidistant** - +0.003 and
+0.004, both intervals covering zero, in both cohorts.

That is CLAUDE.md trap 3 with a named cause rather than a warning. The
thirty-fold spread in proliferation entanglement across MYC signatures has a
companion: **the signature estimators and the regulon estimator disagree about
which mitochondrial programme MYC is nearest to**, and that disagreement is what
produces the split verdict two tables up.

---

## 4. The verdict, on rules fixed before any number existed

The rule, transcribed from section 0 of the script:

> CONTENT if `rho(ox_lvl)` exceeds **both** share rulers in every adjudicating
> estimator-by-cohort cell; SHARE if the reverse in every such cell;
> **ESTIMATOR-DEPENDENT otherwise**. MATERIAL is separate: `|delta| >= 0.05` and
> a bootstrap interval excluding 0, in every such cell.

| item | value |
|---|---|
| **verdict** | **ESTIMATOR-DEPENDENT** |
| adjudicating estimators | `myc_msigdb`, `myc_regulon`, `myc_felsher`, `myc_lowent` |
| dissenting cells | `myc_regulon` / SCAN-B; `myc_regulon` / TCGA |
| material in every adjudicating cell | **no** - smallest `|delta|` 0.027 against a bar of 0.05 |
| `rho(ox_nuc_mtrib) <= 0` or interval covering 0, everywhere | **no** |
| **"overwhelmed, not undone" available** | **no** |

**The prompt's instruction is followed literally: the estimators disagree, the
disagreeing estimator is named, and this stops here.** No tie-break is invented
after the fact. That is what section 0 of the script exists to prevent.

### One thing this result is not

**It is not independent corroboration of E16.** `rho(MYC, ox_lvl) -
rho(MYC, ox_rel)` is the same table E16 section 3.1 reports as an entanglement
difference, read with a sign instead of an absolute value. E16 found that
`ox_rel` buys separability from every signature estimator and loses it against
the regulon, in both cohorts; E20 finds content beating share on every signature
estimator and losing on the regulon, in both cohorts. **These are one finding,
not two.** The genuinely new content of this note is the mitoribosome column -
`ox_nuc_mtrib` and `mtrib_lvl` - which E16 does not carry.

---

## 5. The formulation, and exactly where it failed

The rule was fixed before the panel was computed. The script states it as two
clauses - "(a) CONTENT and MATERIAL, and (b) `rho(ox_nuc_mtrib) <= 0` or its CI
covering 0, in every adjudicating cell of both cohorts" - and it is unpacked
here into the three things that can independently fail.

- **(a) FAILS.** The regulon dissents in both cohorts, materially.
- **(b) FAILS**, on one cell: the regulon in SCAN-B, `ox_lvl - ox_nuc_mtrib` =
  +0.027 [-0.007, +0.063].
- **(c), the script's clause (b), FAILS** on one cell: `rho(myc_msigdb, ox_nuc_mtrib)` in TCGA is
  **+0.109 [+0.050, +0.167]** - positive, interval excluding zero. On the
  Hallmark signature in TCGA, MYC *does* raise the OXPHOS share against the
  mitoribosome, which is the opposite of what the formulation needs.

**Two details a reader needs, neither of which rescues anything.**

The (c) failure **does not replicate**: the SCAN-B counterpart is +0.014 [-0.020,
+0.049], covering zero. A single non-replicating cell breaks the rule as
written, and the rule is applied as written - but the honest description is "one
cohort, one estimator", not "MYC raises OXPHOS share".

And the (a) failure is **larger than the (b) and (c) failures put together**.
The regulon's `-0.089 / -0.098` on `ox_lvl - ox_rel` is a real, replicating,
material disagreement about direction. That is the reason the formulation is
unavailable; the other two are close calls on single cells.

**What the numbers would have to look like for the formulation to be
available**, written down now so the next analysis can be scored against it
rather than argued into:

1. `ox_lvl - ox_rel` positive and material in **all 8** adjudicating cells -
   which requires the regulon column to change sign, in both cohorts.
2. `rho(est, ox_nuc_mtrib)` at or below zero, or interval-covering-zero, in all
   8 - which requires the TCGA `myc_msigdb` cell to fall from +0.109 to
   interval-covering-zero.

---

## 6. What would change this answer

1. **`COLLECTRI_MYC_STIM`, 739 genes, in `data/collectri_human/` and never
   scored.** The whole verdict turns on **one** regulon estimator, `M_b`.
   Scoring a second one is the cheapest thing that could settle whether the
   dissent is a property of *regulon estimators* or a property of *this
   regulon*. **E16 section 9 already named this the cheapest thing that could
   overturn its own section 3.1; E20 gives the same item a second, independent
   reason.** Handoff section 5 items 1 and 5.
2. **A fourth signature estimator.** Three signature estimators agreeing is
   three points, and they share genes with each other. `myc_msigdb` is 23.5%
   proliferation-entangled and `myc_lowent` is 1.5%, so the panel does span the
   entanglement axis - but it does not span the *curation* axis.
3. **A denominator between the two.** `ox_rel`'s denominator is 1,047 genes and
   `ox_nuc_mtrib`'s is 83. The regulon dissents on the first and not the second,
   and nothing here says where between them it changes its mind. The mouse's
   `relify()` builds exactly the intermediate objects.
4. **A technical control on `ox_nuc_mtrib`.** Section 1's caveat: it is a
   difference of two composites correlating at 0.865 / 0.905, and transcript length,
   GC content and expression level all differ systematically between MRPs and
   OXPHOS subunits. A size- and expression-matched pseudo-set drawn from
   MitoCarta would say how much of its variance survives those.
5. **SCAN-B purity.** Not available (trap 2), never imputed. E16 check 4 showed
   the twelve-gene block is infiltrate-robust in TCGA, but no ruler in *this*
   note has been run through that control.

---

## 7. What was NOT done, listed so it stays a decision

No re-scoring of anything - GSVA and mitoPPS were read. No pooling of scores
across cohorts, and no comparison of values across them. No interaction model,
no MYC stratification, no subtype split, no new endpoint or ratio, no
Johnson-Neyman. No per-gene work: the twelve BCL2-family genes do not appear in
this note at all. The mouse was not read, compared or cited as a value. Nothing
written to `myc_human_validation` or to `myc_mouse`. No ortholog function called,
in either direction - the tripwire returns nothing.

---

## 8. Where the numbers live

| | |
|---|---|
| script | `scripts/E20_share_vs_content.R` |
| object | `results/share_vs_content.rds` - `$panel`, `$contrast`, `$repro`, `$ruler_agreement`, `$adjudicate`, `$verdict_panel`, `$verdicts`, `$verdict`, `$dissent`, `$formulation_available`, `$set_sizes`, `$overlap_audit`, `$settings`, `$rules` |
| table | `outputs/tables/E20_share_vs_content.csv` |
| figures | none - this is a re-read, and every number in it is a table |
| the control | `$repro`: 48 of E16's raw separability cells, max \|delta\| = 0 exactly - **confirmed in the author's run** |
| read for a different question | `docs/2026-09-04_e16_respiratory_rulers.md` sections 3 and 3.1 |
