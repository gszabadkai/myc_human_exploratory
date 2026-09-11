# E25 - the pathway ordering, group by group: MYC, not age, is what the human cohorts match

**Run and written 2026-09-11. EXPLORATORY, POST-HOC, DESCRIPTIVE. Nothing here
is pre-registered and nothing here is a hypothesis test.**

Script: `scripts/E25_group_ranking_comparison.R`. **Sourced by the author
2026-09-11**; `results/mitopps_group_rankings.rds`,
`outputs/tables/E25_group_orderings.csv` and two figures under
`outputs/mitopps_rank/` are on disk. The numbers below were first written from a
dry run in the scratchpad. The author's run reproduces **all 24 saved elements
identically**, the seeded 5,000-draw bootstraps included, so every number here
is read from the saved object and **the dry-run caveat is discharged**.

**RANK, NEVER VALUE.** No mitoPPS value is compared across the species boundary.
Every cross-species quantity is an ordering, a rank correlation between
orderings, or a percentile position inside one. Mouse and human were never
scored together; TCGA and SCAN-B are never averaged.

**This note supersedes E24's readings.** See
`docs/2026-09-10_e24_mitopps_rank_comparison.md`, which now carries a banner.

---

## 0. The answer in six lines

| | |
|---|---|
| **What was actually compared?** | The **ordering of the mitochondrial pathway panel** in each of the four mouse groups and in MYC-high against MYC-low tumours in each human cohort. Eight orderings of the same 142 shared pathways. |
| **Does each human group match a mouse group?** | **Yes, and the assignment is separated from its runner-up in all four cases.** MYC-low matches **mouse 12W WT**; MYC-high matches **mouse 6W Myc+**. Identical in both cohorts. |
| **Where does the respiratory chain sit?** | MYC-low **16.6** (TCGA) and **15.1** (SCAN-B) against mouse 12W WT at **7.4**. MYC-high **68.0** and **74.3** against mouse 6W Myc+ at **89.1**. |
| **Do the MYC contrasts order the panel alike?** | **Yes, both cohorts, against the mouse Myc effect at both timepoints.** rho 0.39 to 0.54, every interval clear of zero. |
| **Is the rank size-structural, as in E24?** | **No. The problem is gone, not managed.** No arm exceeds \|rho\| = 0.21 against log set size, against 0.54 and 0.57 for E24's statistic. |
| **Was it MYC all along?** | **That is what this says.** Pooled over MYC status, as E24 had them, the human cohorts average into the middle of the panel and read as nothing in particular. Split on MYC, they separate into two groups that match two different mouse groups. |

**The one-sentence reading: human breast tumours do not sit at a single position
on the mouse's developmental axis - they sit at two, and which one depends on
MYC, not on anything about age.**

---

## 1. Why the estimand changed, and why that is a withdrawal rather than a tweak

E24 ranked **one pathway inside each sample** and summarised the percentile by
group. The mouse figures
(`fig02_reallocation_ranked.R`, `figS6_reallocation_temporal_ranked.R` in
`myc_mouse`) do something else: they rank **pathways** by a **group contrast** in
mitoPPS. Those are different statistics, and only one of them is size-biased.

| rho(statistic, log set size) | mouse | TCGA | SCAN-B |
|---|---|---|---|
| E24, within-sample percentile | +0.165 | **+0.536** | **+0.571** |
| E25, group ordering | +0.04 to -0.14 | -0.201 | -0.182 |
| E25, group contrast | +0.03 to -0.15 | +0.192 | +0.196 |

**E24 compared two differently biased statistics across the species boundary.**
That is a stronger objection than its own size gate raised, because a size gate
can be remedied with a size-matched panel while an asymmetry cannot.

**The cause is the rank transform, not mitoPPS.** The script asserts it rather
than assuming it: mitoPPS is centred at **exactly 1 for every pathway** in both
species (max \|delta\| 3.3e-16 mouse, and asserted in both human cohorts), and
both implementations take the **mean** per set - `myc_mouse` script 08 line 386
and this repo's `.path_scores`. So the first moment is size-free on both sides
and **there is no asymmetry in how set size enters the score.** What is not
size-free is the cross-sample **variance**, which shrinks with set size
(rho(sd, log size) = -0.509 in the mouse). Pairwise ratios are right-skewed, so
small high-variance pathways sit below the crowd more often than above, and a
**within-sample** rank transform converts that variance difference into a
mean-rank difference. It scales with cohort heterogeneity, which is why human
tumours showed three times the mouse effect. Ranking pathways by a group mean or
a group contrast does not do it.

**So the author's second point is answered, and the answer is not the one the
question implied.** There is no missing size normalisation in the mouse set. The
size problem was created by E24's own choice of statistic and it disappears when
the statistic matches the mouse figures'.

---

## 2. The three gates, all prior to any reading

**GATE A - the mouse side reproduces.** Recomputing all four contrasts from
`mitopps_scores` reproduces the artefact's own `$mitopps_pairwise` at
**max \|delta\| = 4.44e-16** for every contrast. If this had failed, the snapshot
was being read wrongly and nothing below it would mean anything.

**GATE B - the panels and the hierarchy agree.** 142 shared pathways, matching
E24's saved intersection exactly. The **full MitoCarta hierarchy path** - not
merely the depth - agrees on all **141** shared catalogue pathways, **7 level-1,
38 level-2, 96 level-3 in both species**, each read from its own native source:
the human MitoCarta 3.0 workbook for the human side, the artefact's own
`pathway_levels` for the mouse side. **No ortholog function is called.**

**GATE C - the size gate, on both estimands and every panel.** Clear in every
arm: the largest absolute correlation anywhere on the primary panel is 0.201.

### 2.1 The antichain, and a correction to my own first draft

The author's rule was that size should apply only to single pathways. That is
implemented as the **antichain**: the **114** pathways of which no member
contains another, built independently from each species' hierarchy and
**identical from both**.

**It is not a single level, and a level-only panel gets it wrong.** The first
draft used the level-3 slice and silently dropped the primary set, which sits at
**level 2 with no descendants in the panel** - a legitimate antichain member that
no single level contains. The script now reports readout membership per panel so
an absence can never again be an empty table.

**The antichain barely moves anything**, which closes the hierarchy question
rather than managing it: the primary's position shifts from 60.9 to 59.2 in
mouse 6W WT, 7.4 to 8.3 in 12W WT, 16.5 to 17.1 in TCGA MYC-low and 74.3 to 76.8
in SCAN-B MYC-high.

**Re-ranking, never re-scoring.** Every panel is a subset of the same scores. The
scoring universe stays 142 human and 144 mouse. Mouse mitoPPS is not recomputed
in this repo, and for symmetry neither is human.

---

## 3. The result: two human groups, two mouse groups

### 3.1 The mapping

| human group | nearest mouse ordering | rho | 95% CI | runner-up | its rho | separated |
|---|---|---|---|---|---|---|
| TCGA MYC-low | **mouse 12W WT** | 0.61 | 0.47 to 0.73 | mouse 6W WT | 0.08 | yes |
| SCAN-B MYC-low | **mouse 12W WT** | 0.64 | 0.52 to 0.75 | mouse 6W WT | 0.16 | yes |
| TCGA MYC-high | **mouse 6W Myc+** | 0.46 | 0.30 to 0.60 | mouse 12W Myc+ | 0.00 | yes |
| SCAN-B MYC-high | **mouse 6W Myc+** | 0.51 | 0.36 to 0.64 | mouse 12W Myc+ | 0.09 | yes |

**Two of the four mouse groups have no human counterpart at all.** Every
cross-species pair involving **mouse 6W WT** or **mouse 12W Myc+** has an
interval covering zero. The pre-window wild-type gland and the late transgenic
tumour match nothing on the human side.

### 3.2 The respiratory chain's position

| arm | shared (142) | 95% CI | antichain (114) |
|---|---|---|---|
| mouse 6W WT | 60.9 | 20.1 to 93.3 | 59.2 |
| mouse 12W WT | **7.4** | 7.4 to 27.8 | 8.3 |
| mouse 6W Myc+ | **89.1** | 70.8 to 91.2 | 88.2 |
| mouse 12W Myc+ | 32.7 | 14.4 to 62.3 | 32.9 |
| TCGA MYC-low | **16.5** | 9.5 to 28.5 | 17.1 |
| TCGA MYC-high | **68.0** | 56.7 to 77.1 | 69.7 |
| SCAN-B MYC-low | **15.1** | 12.3 to 18.7 | 16.2 |
| SCAN-B MYC-high | **74.3** | 65.8 to 78.5 | 76.8 |

The mouse intervals are wide because n = 6 per arm. The human intervals are
narrow because n is 365 and 1,069 per group, which makes the human positions
precise and the mouse reference positions still imprecise.

### 3.3 The structure that makes the mapping credible

The human MYC axis reproduces the mouse's **most extreme internal contrast**.

| pair | rho |
|---|---|
| mouse 12W WT against mouse 6W Myc+ | **-0.95** |
| SCAN-B MYC-low against SCAN-B MYC-high | **-0.95** |
| TCGA MYC-low against TCGA MYC-high | -0.92 |
| TCGA MYC-low against SCAN-B MYC-low | +0.88 |
| TCGA MYC-high against SCAN-B MYC-high | +0.88 |

The two human cohorts agree with each other at +0.88 on both sides of the split,
and each human pole is a near-mirror of the other exactly as the two matched
mouse groups are. **R2 is satisfied: the cohorts agree, and the agreement is not
a coincidence of one estimator.**

---

## 4. The contrast comparison - the mouse figures' own statistic

| human contrast | mouse contrast | kind | rho | 95% CI | agrees |
|---|---|---|---|---|---|
| TCGA MYC-high minus low | `Myc_effect_6W` | genotype, within timepoint | 0.39 | 0.23 to 0.54 | yes |
| TCGA MYC-high minus low | `Myc_effect_12W` | genotype, within timepoint | 0.47 | 0.31 to 0.61 | yes |
| SCAN-B MYC-high minus low | `Myc_effect_6W` | genotype, within timepoint | 0.46 | 0.31 to 0.60 | yes |
| SCAN-B MYC-high minus low | `Myc_effect_12W` | genotype, within timepoint | 0.54 | 0.40 to 0.65 | yes |

**Q1 is answered the same way in both cohorts: the human MYC contrast orders the
mitochondrial panel like the mouse Myc effect, at both mouse timepoints.**

Where the primary set falls inside each contrast's own ordering, primary panel:

| contrast | kind | percentile |
|---|---|---|
| `Myc_effect_6W` | genotype, within timepoint | 75.7 |
| `Myc_effect_12W` | genotype, within timepoint | 78.5 |
| TCGA MYC-high minus low | MYC activity, within cohort | 71.5 |
| SCAN-B MYC-high minus low | MYC activity, within cohort | 80.6 |
| `Temporal_Myc-` | TIMELINE, batch-confounded | 9.5 |
| `Temporal_Myc+` | TIMELINE, batch-confounded | 8.8 |

**MYC promotes the respiratory chain in all four genotype and activity
contrasts, and the timeline demotes it in both.** Those are opposite directions,
which is the whole finding in one table.

---

## 5. The timeline, included by decision, with its confounder stated

**The author ruled that the timeline comparisons stay in**, on the grounds that
reprioritisation is a timeline effect measurable only that way and that DESeq2
normalisation is what handles the batch. They are computed and reported in full
and nothing is withdrawn on the confounder's account. The confounder is recorded
wherever a temporal number appears, here and in the saved object and on both
figures.

**What the confounder is.** In `myc_mouse` the time axis is batch-confounded,
batch equals timepoint. That repo's own `figS6` header calls the individual 6W to
12W changes descriptive rather than claims, and its evidence audit notes the
batch axis is genotype-independent and so cancels inside each timepoint's
contrast. The genotype contrast within a timepoint is therefore the cleaner axis.
That is recorded as a fact about the two axes, not as a reason to prefer one.

**What the timeline says here.** Both temporal contrasts run **opposite** to the
human MYC contrast: -0.42 and -0.52 against TCGA, -0.38 and -0.51 against
SCAN-B, all four intervals clear of zero. The wild-type ordering does change
across the timeline while staying positively correlated with itself
(6W WT against 12W WT, rho +0.33).

**How to read the two together.** The human MYC contrast matches the genotype
axis positively and the timeline axis negatively, and both readings point the
same way: what human tumours share with the mouse is the **MYC** move, not the
**developmental** move. The timeline's inclusion strengthens that rather than
complicating it, because the timeline is where the human cohorts fail to match.

---

## 6. What moves the answer, and what does not

**The estimator does, and it is trap 3 and trap 4 exactly.** MYC-high minus
MYC-low shift in the primary's percentile, TCGA and SCAN-B:

| estimator | TCGA shift | SCAN-B shift |
|---|---|---|
| `MYC_UP.V1_UP` (1.5% entangled) | +79.6 | +67.6 |
| `FELSHER` (14.8%, the primary) | +51.4 | +59.2 |
| `HALLMARK_MYC_TARGETS_V1` (23.5%) | +85.9 | +79.6 |
| CollecTRI regulon | +57.7 | +28.2 |
| **MYC mRNA** | **-5.6** | **-33.8** |

**All four activity estimators agree in sign and are large. MYC mRNA gives the
opposite sign in both cohorts.** Trap 4 says MYC mRNA is not MYC activity and the
difference is total; here the gap is a sign flip, and the regulon's weaker
SCAN-B value is consistent with E21's finding that it is size-dominated.

**The DNA-level split is sharper than any activity split.** TCGA
`MYC_amp`: **13.0** in the 816 non-amplified against **85.6** in the 224
amplified. This is the closest human analogue of a transgene and it is the one
SCAN-B structurally cannot provide.

**The antichain does not move it** (section 2.1). **Dropping the mtDNA pathway
does not either** - it is a named sensitivity panel and the primary's ordering is
unchanged in kind.

**The mtDNA-encoded pathway runs opposite to the nuclear subunits throughout**,
which is trap 8 behaving as documented: mouse 12W WT **99.6** and 12W Myc+
**98.9** against 6W WT 6.7 and 6W Myc+ 0.4; human MYC-low 91.2 (TCGA) and 70.1
(SCAN-B) against MYC-high 13.0 and 12.3. **The mitoribosome tracks the nuclear
subunits, not the mtDNA pathway** (MYC-high 81.3 and 89.1 against MYC-low 12.3
and 8.8).

---

## 7. What this does not address, and what would falsify it

> **ANSWERED 2026-09-11 by `docs/2026-09-11_e26_purity_and_composition.md`.**
> The purity question below was run and **the separation survives**. The gate
> that licenses saying so: proliferation keeps 96 per cent of its MYC difference
> under adjustment while the adjusters explain 8.8 per cent of its variance, so
> the adjustment is not over-correcting. The separation goes 44.4 unadjusted to
> 43.0 and 47.2 adjusted, the respiratory chain attenuates to the 41st
> percentile of the panel's own attenuation distribution, and the cross-species
> agreement holds under every model. **E26 also found that the exposure is
> adipose rather than purity** - the MYC groups differ by 0.058 in purity and
> 0.743 standard units in an adipocyte marker score - and the pre-declared
> tracer rule promoted the adipose-adjusted models accordingly. **SCAN-B still
> cannot be checked** and that paragraph below stands unchanged.

**Purity is not addressed and every human number above is unadjusted.** That is
the next script's subject and **trap 15 governs it**: the mouse fat-pad result
put `Mki67` at tau +0.313 before adjustment and -0.096 after, while the adjusters
explained only 15% of its variance, so a **positive control is declared and read
first**, before any adjusted OXPHOS number is looked at. TCGA carries purity,
leukocyte fraction and ploidy for 1,020 of 1,095 samples; **SCAN-B carries none
and nothing is imputed**, so the purity arm is structurally TCGA-only.

**This is where the finding is most exposed.** MYC-high tumours are more
proliferative and plausibly more tumour-cell-pure, and adipose is OXPHOS-high, so
part of the MYC-high to MYC-low separation could be composition rather than
programme. The falsification is concrete: if the separation collapses under
purity adjustment in TCGA **while the pre-declared positive control survives**,
the mapping in section 3 is a composition artefact. If the control also collapses,
the adjustment is over-correcting and neither reading stands.

**Two other falsifiers.** A human cohort scored on this panel whose MYC-high
ordering fails to match mouse 6W Myc+, which would make section 3.1 cohort-
specific. And a mouse experiment with more animals in which the 12W WT ordering
stops being the near-mirror of 6W Myc+, which is the structure that makes the
human split interpretable at all.

---

## 8. Constraints observed

- **Rank only.** No mitoPPS value crosses the species boundary in the object, the
  table or the figures.
- **No pooling.** Mouse and human never scored together; TCGA and SCAN-B never
  averaged; species is a cohort (trap 6).
- **The mouse artefact is read-only input.** md5 `b8a125af4bc0d5909ad03f6e126e2890`
  re-checked against `data/from_myc_mouse/README.md` at run time; mouse mitoPPS
  not recomputed; the mouse repo not read at run time and its working tree and
  HEAD (`129b0c0`) exactly as found.
- **No ortholog call anywhere.** Tripwire clean. Each species' hierarchy read
  from its own native source, and their agreement is GATE B.
- **SCAN-B scored through `scanb_pheno.rds$symbol_map`** (trap 7), inherited
  with the saved universe.
- **`mt-` / `MT-` genes held separately in both species** (trap 8), and their
  pathway is a named sensitivity here.
- **No MYC statement rests on one signature** (trap 3): five estimators reported.
- **N3.** Nothing here is called primed. Every number is a transcript ordering.
- **The figures carry no significance marks** and nothing is read off them
  statistically, per the author's instruction that they are for insight.

## 9. Provenance and reproduction

- Script: `scripts/E25_group_ranking_comparison.R`, seed 1, 5,000 bootstrap draws.
- Mouse input: `data/from_myc_mouse/mitopps_scores.rds`, produced by `myc_mouse`
  `scripts/08_mitoPPS_analysis.R` at commit `0b9b28e`.
- Human inputs: the saved 142-pathway mitoPPS universes, TCGA's from
  `data/from_validation/tcga_brca_mito_scores.rds` and SCAN-B's from
  `results/scanb_scores.rds`; MYC estimators from `results/new_set_scores.rds`
  and `results/scanb_scores.rds`; the amplification call from
  `data/from_validation/tcga_brca_myc_scores.rds`.
- Hierarchy: `data/mitocarta_human/Human.MitoCarta3.0.xls` sheet C for the human
  side, the artefact's own `pathway_levels` for the mouse side.
- Saved: `results/mitopps_group_rankings.rds` (24 elements),
  `outputs/tables/E25_group_orderings.csv`,
  `outputs/mitopps_rank/E25_group_orderings.pdf`,
  `outputs/mitopps_rank/E25_contrast_rank_vs_rank.pdf`.
