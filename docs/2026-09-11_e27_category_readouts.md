# E27 - the seven level-1 categories: the internal null holds on the rule as written, and does not survive the stricter statistic

**2026-09-11. EXPLORATORY, POST-HOC, DESCRIPTIVE**, as the whole of Phase 2 is.
Nothing here is pre-registered and nothing is a hypothesis test. R1 to R5 are
**reading rules fixed in advance so that the outcome could not be chosen
afterwards. They are not predictions and none of them was predicted.**

Script: `scripts/E27_category_readouts.R`. **Sourced by the author 2026-09-11**;
`results/category_readouts.rds`, `outputs/tables/E27_category_readouts.csv` and
two figures under `outputs/mitopps_rank/` are on disk. The author's run
reproduces **all 21 saved elements identically** against the dry run, so every
number here is read from the saved object.

**Nothing was re-scored.** The level-1 categories are themselves rows of the
142-pathway universe and were already scored. `E25` is untouched. Design
rationale: `docs/2026-09-11_level1_category_comparison_rationale.md`.

---

## 0. The answer in seven lines

| | |
|---|---|
| **Does OXPHOS fall across the mouse wild-type window?** | **Yes.** -0.0926 [-0.1841, -0.0033] in mitoPPS units, interval clear of zero. `OXPHOS subunits` alone: -0.1095 [-0.2173, -0.0030]. |
| **Does the internal null hold?** | **On the rule as written, yes. On a stricter statistic, no.** Central dogma moves -0.0108 [-0.0654, +0.0545], which lies inside the equivalence bound of 0.0926. But the **paired** difference of magnitudes is 0.0819 [-0.0194, +0.1584] and **covers zero**. |
| **Does Metabolism rise?** | **Yes, barely.** +0.0313 [0.0003, 0.0609]. On its own this confirms nothing. |
| **Is the reallocation clean?** | **No, and R3's second branch is live.** Three of the six non-OXPHOS categories fall alongside OXPHOS. This is not one-down-one-up. |
| **Is the mouse-human ordering difference real at category level?** | **Yes.** OXPHOS above central dogma in mouse 6W Myc+ (82.7 against 70.8) and below it in **both** human cohorts (62.3 against 79.9; 67.3 against 82.0). |
| **Do the OXPHOS halves separate?** | **In three of four windows, and the flatness is real not structural.** Assembly factors' mean absolute deviation from baseline is **0.0095** in human and **0.0208** in mouse, against subunits at **0.0503** and **0.0657**. |
| **Is the partition disjoint?** | **No.** OXPHOS and Metabolism share **23** genes, 147 genes sit in more than one category, and 13 catalogue genes sit in none. |

**The one-sentence reading: on the native MitoCarta partition the mouse
developmental change is visible and the internal null does the work the design
asked of it under the pre-specified rule, but not under a statistic that
propagates uncertainty in both quantities, and the reallocation is not clean.**

---

## 1. The control, first

**88 recomputed percentiles reproduce `E25` exactly, max |delta| = 0** - eleven
readouts across eight group orderings. If this had failed every number below
would have been a rebuild rather than a re-read. It is the control on the whole
script.

---

## 2. The partition, and an asymmetry wider than E24 recorded

Read from source on both sides: the human MitoCarta 3.0 workbook for the human
hierarchy, the artefact's own `pathway_levels` for the mouse.

| readout | role | human | mouse | delta |
|---|---|---|---|---|
| OXPHOS | level-1, design | 156 | 155 | +1 |
| Mitochondrial central dogma | level-1, design | 231 | 230 | +1 |
| Metabolism | level-1, design | 459 | 462 | **-3** |
| Mitochondrial dynamics and surveillance | level-1, context | 104 | 102 | **+2** |
| Protein import, sorting and homeostasis | level-1, context | 86 | 88 | **-2** |
| Signaling | level-1, context | 48 | 47 | +1 |
| Small molecule transport | level-1, context | 85 | 85 | 0 |
| OXPHOS subunits | child of OXPHOS | 89 | 89 | 0 |
| OXPHOS assembly factors | child of OXPHOS | 68 | 67 | +1 |
| mtDNA-encoded OXPHOS subunits | synthetic, outside the umbrella | 13 | 13 | 0 |
| Mitochondrial ribosome | grandchild of central dogma | 83 | 83 | 0 |

**Seven of eleven readouts differ in size across species, not two.** E24
recorded the umbrella at 156 against 155 and assembly factors at 68 against 67
and judged them harmless because they sat outside the primary set. **They no
longer do, and the asymmetry is wider than that pair.** Metabolism differs by
three genes and two context categories by two.

**A note on the saved object.** Its `$notes` field quotes only the pair E24 had
recorded. That sentence is true but narrower than the `$partition` table beside
it, **and the table is the authority.** The script was deliberately not edited
after the author's run, so that script and saved object stay in exact
correspondence.

**Which gene differs cannot be reported.** Naming it requires aligning mouse and
human symbols, which is an ortholog mapping whatever it is called, and the
tripwire governs. Within-species membership and cross-species counts only.

**The mitoribosome is a grandchild**, not a child: `Mitochondrial central dogma >
Translation > Mitochondrial ribosome`.

---

## 3. Overlap and containment

**The level-1 partition is not disjoint**, and the overlap sits on the pair the
reallocation claim is about.

| pair | shared genes | share of the first | share of the second |
|---|---|---|---|
| Metabolism and Small molecule transport | 41 | 8.9% | 48.2% |
| **OXPHOS and Metabolism** | **23** | **14.7%** | **5.0%** |
| Signaling and Small molecule transport | 20 | 41.7% | 23.5% |
| Dynamics and Protein import | 16 | 15.4% | 18.6% |

Union of the seven: 1,022 against a size sum of 1,169, so **147 genes sit in
more than one category**. Of 1,035 catalogue genes, **13 sit in no level-1
category at all**.

**Reported, not gated**, as the rule required. **The direction matters and is
not neutral:** shared genes push two categories toward the *same* value, so
OXPHOS falling while Metabolism rises happens **despite** the overlap rather
than because of it. The overlap is a conservative confound for a dissociation
and an inflating one only for concordance.

**The qualifier that argument needs:** it holds if the 23 shared genes behave
like OXPHOS-typical genes. If they are atypical it weakens. **A disjointified
recomputation would settle it and is feasible on the human side only** - E24's
PART A2 rebuilt the human universe bit-equal, but the mouse artefact carries
pathway scores and not gene-level values, so the mouse half would mean touching
that repo. **Not run. Available if the objection is pressed.**

**Containment is exact.** 89 subunits plus 68 assembly factors minus one gene in
both equals the umbrella's 156, with **zero** umbrella genes in neither half and
**zero** mtDNA-encoded genes inside it. **The umbrella therefore carries no
independent information.** It is an aggregate of a mobile half and an immobile
one and is never a third measurement. Level-1 categories contain their
descendants by design; that is containment, not overlap, and is expected.

---

## 4. R1, R2, R3, R5 - on values, within species

The mouse wild-type window, all seven categories, in mitoPPS units:

| category | change | 95% CI | clear of zero |
|---|---|---|---|
| **OXPHOS** | **-0.0926** | -0.1841 to -0.0033 | **yes** |
| Protein import, sorting and homeostasis | -0.0555 | -0.1178 to +0.0092 | no |
| Signaling | -0.0318 | -0.0903 to +0.0285 | no |
| **Mitochondrial central dogma** | **-0.0108** | -0.0654 to +0.0545 | no |
| Mitochondrial dynamics and surveillance | +0.0077 | -0.0650 to +0.0836 | no |
| Small molecule transport | +0.0118 | -0.0314 to +0.0527 | no |
| **Metabolism** | **+0.0313** | +0.0003 to +0.0609 | **yes** |

**R1 holds.** OXPHOS falls with an interval clear of zero. `OXPHOS subunits`
alone falls further, -0.1095 [-0.2173, -0.0030].

### 4.1 R2, and the part that must not be reported alone

**On the pre-specified rule: STRONG.** Central dogma's interval lies entirely
inside the equivalence bound of 0.0926, the magnitude of the OXPHOS change.

**On a statistic that propagates uncertainty in both quantities: not
established.** The paired difference of magnitudes, computed on the same
bootstrap resamples, is **0.0819 [-0.0194, +0.1584]** and **covers zero**.

The two disagree because the equivalence bound treats the OXPHOS change as a
fixed number while the paired statistic does not. **The pre-specified rule
stands as the verdict, because it was fixed before any number and changing it
afterwards is exactly what the rule exists to prevent. Both numbers are
reported, and the honest summary is that the null panel works on the rule as
written and does not separate the two categories at n = 6 under the stricter
statistic.** This is never written as a demonstrated null.

### 4.2 The alternative null, read and not designed for

**`OXPHOS assembly factors` moves -0.0016 [-0.0511, +0.0556]** across the same
window, against central dogma's -0.0108. A sibling inside the OXPHOS parent
controls for the respiratory branch specifically rather than for mitochondrial
pathways in general, and it moves less than the pre-specified null does. Its
paired difference behaves the same way, 0.0910 [-0.0168, +0.1666], covering
zero.

**This was read, not designed for**, at the author's instruction, and it is an
observation rather than a replacement for the pre-specified null.

### 4.3 R3 and R5

**Metabolism rises**, +0.0313, interval clear of zero but only just. **On its own
this confirms nothing** - a rise is partly forced by the relativity of the score,
which is what R2 exists to control.

**R5: the reallocation is not clean, and R3's second branch is live.** Of the six
non-OXPHOS categories three rise and three fall, and **Protein import and
Signaling fall alongside OXPHOS**. Only OXPHOS and Metabolism have intervals
clear of zero. So this is not one-down-one-up; it is a broad settling in which
OXPHOS falls furthest and Metabolism rises most.

### 4.4 Where the two estimands disagree

**On `Mitochondrial dynamics and surveillance`**: +0.0077 on values, exactly 0.00
on percentiles. Its percentile is identical at 6W and 12W wild-type. This is the
only disagreement among the seven, it is small on both scales, and the interval
covers zero on values, so nothing rests on it. An exact percentile tie is coded
**flat**, not down; coding it as a fall would invent a direction the data does
not have.

---

## 5. R4 - the ordering difference is present at category level

| group | OXPHOS | central dogma | OXPHOS above? |
|---|---|---|---|
| mouse 6W Myc+ | 82.7 | 70.8 | **yes** |
| TCGA MYC-high | 62.3 | 79.9 | no |
| SCAN-B MYC-high | 67.3 | 82.0 | no |

**What `E25`'s figure showed for one pathway is also true of the parent
categories, in both cohorts.** The same thing appears in the contrast orderings:

| contrast | kind | OXPHOS | central dogma | Metabolism |
|---|---|---|---|---|
| `Myc_effect_6W` | genotype, within timepoint | 69.4 | 68.0 | 39.8 |
| `Myc_effect_12W` | genotype, within timepoint | 72.9 | 64.4 | 37.0 |
| TCGA MYC-high minus low | MYC activity | 68.7 | **77.8** | 43.3 |
| SCAN-B MYC-high minus low | MYC activity | 73.6 | **83.5** | 46.1 |
| `Temporal_Myc-` | **TIMELINE, batch-confounded** | 13.7 | 42.6 | 69.4 |
| `Temporal_Myc+` | **TIMELINE, batch-confounded** | 11.6 | 29.2 | 67.3 |

**In the mouse genotype contrasts OXPHOS sits above central dogma; in both human
MYC contrasts it sits below.** Neither outcome was predicted and neither is the
better result. It is one descriptive observation.

The table also shows the two axes pointing opposite ways at category level: the
MYC contrasts put OXPHOS and central dogma high and Metabolism low, and the
timeline contrasts do the reverse.

---

## 6. R4b - the halves separate, and the flatness is real

| window | species | subunits | assembly factors | gap | 95% CI | separate |
|---|---|---|---|---|---|---|
| mouse WT 6W to 12W | mouse | -0.1095 | -0.0016 | -0.1079 | -0.2401 to +0.0221 | no |
| mouse Myc+ 6W to 12W | mouse | -0.1533 | -0.0290 | -0.1243 | -0.2550 to -0.0157 | **yes** |
| TCGA MYC-high minus low | human | +0.1137 | +0.0171 | +0.0966 | +0.0720 to +0.1212 | **yes** |
| SCAN-B MYC-high minus low | human | +0.0873 | +0.0209 | +0.0665 | +0.0569 to +0.0761 | **yes** |

**Three of four windows separate**, named individually rather than generalised.
The exception is the mouse wild-type window, where the gap is the second largest
of the four but the interval covers zero at n = 6.

**The flatness is real, not percentile room.** That was the caveat R4b was
written on the value scale to answer, and the baseline check settles it: mean
absolute deviation from the cohort baseline of 1 is **0.0095** for assembly
factors in human against **0.0503** for subunits, and **0.0208** against
**0.0657** in mouse. **Assembly factors are five times less mobile in human and
three times less in mouse**, on a scale with no ceiling or floor.

**So the respiratory chain and the machinery that builds it are prioritised
differently** - a mass-versus-function distinction one level below the one the
title makes, on the native ontology, in both species.

---

## 7. What this panel is for, and what it is not for

**The sentence it is intended to carry:**

> On the native MitoCarta level-1 partition, the mouse developmental window
> lowers the priority of OXPHOS while central dogma does not move measurably and
> metabolism rises - and within OXPHOS the move is carried by the respiratory
> subunits, not by the assembly machinery.

**The sentences it must not be used to carry:**

1. **It does not establish the developmental change.** The
   expression-matched-null analysis in `myc_mouse` remains the primary evidence.
   **This is a presentation of an existing claim on a cleaner ontology, not a
   new result.** Both stay reportable and **if they disagree that disagreement
   is the finding**, not something to be resolved by choosing the friendlier one.
2. **It does not demonstrate that central dogma is unchanged.** The
   pre-specified rule returns STRONG; the stricter paired statistic does not.
   The honest form is that central dogma's change is small and not separable
   from OXPHOS's at n = 6 under a statistic that propagates both uncertainties.
3. **It does not show a clean reallocation to metabolism.** Three of six
   non-OXPHOS categories fall alongside OXPHOS.
4. **It does not cross the species boundary on values.** The eight-group figure
   is on percentiles for that reason, and the value figure gives mouse and human
   independent axes.
5. **It says nothing about purity.** `E26` addressed composition for `OXPHOS
   subunits` in TCGA only. **Whether these category readouts survive adjustment
   is open**, and SCAN-B structurally cannot be checked at all.

---

## 8. The batch confounder, and what this design does to it

The mouse timeline is **batch-confounded**, batch equals timepoint. The author's
ruling stands: temporal numbers are reported in full, the confounder is stated
wherever they appear, and **nothing is withdrawn on its account.**

**What this design adds, stated as a narrowing and not an exclusion:** a batch
effect would have to act **differentially across three mitochondrial categories
in the same libraries** - lowering OXPHOS, raising Metabolism, leaving central
dogma untouched - and within OXPHOS it would have to move the subunits five
times as far as the assembly factors. That is a far stronger requirement than a
global batch effect. **It does not eliminate the confound; it narrows what the
confound would have to be.**

---

## 9. Two things worth knowing before the figure is used

**The mtDNA-encoded pathway is too volatile to read across the window.** Its
wild-type change is +0.525 with an interval of -0.530 to +1.993. Thirteen genes
on a ratio-based score gives an interval wider than the quantity. It appears in
the decomposition panel as context and no statement rests on it.

**Disclosure on what was visible in advance.** The eleven readouts' **percentile**
positions were computable from `E25`'s saved object - they are rows of its
`$orderings` - though `E25`'s note reported only four of them. They were
surfaced during planning, **before** this run rather than after it. The reading
rules were fixed by the author before any number was seen, and that is what
protects them.

---

## 10. Constraints observed

- **Nothing re-scored.** The categories were already scored; this re-reads and
  re-aggregates. The mouse artefact is read-only and its md5 was re-checked.
- **`E25` not modified.**
- **Rank only across the species boundary.** No mitoPPS value appears in a
  cross-species table, axis or panel.
- **TCGA and SCAN-B separate**, never averaged, never pooled.
- **No ortholog call.** Tripwire clean; each species' hierarchy read from its
  own source.
- **`mt-` / `MT-` held separately** (trap 8).
- **N3.** Transcript associations throughout. "Primed" appears nowhere.
- **The figures carry no significance marks.**

## 11. Provenance

- Script: `scripts/E27_category_readouts.R`, seed 1, 5,000 bootstrap draws.
- Inputs: `results/mitopps_group_rankings.rds` (E25),
  `data/from_myc_mouse/mitopps_scores.rds` (md5
  `b8a125af4bc0d5909ad03f6e126e2890`),
  `data/from_validation/tcga_brca_mito_scores.rds`, `results/scanb_scores.rds`,
  `results/new_set_scores.rds`, `data/mitocarta_human/Human.MitoCarta3.0.xls`.
- Saved: `results/category_readouts.rds` (21 elements),
  `outputs/tables/E27_category_readouts.csv`,
  `outputs/mitopps_rank/E27_category_positions.pdf`,
  `outputs/mitopps_rank/E27_equivalence_values.pdf`.
