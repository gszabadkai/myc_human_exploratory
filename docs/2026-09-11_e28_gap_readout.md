# E28 - the gap readout does not survive its own tests, and the note says so

**2026-09-11. EXPLORATORY, POST-HOC, DESCRIPTIVE**, as all of Phase 2 is.
Nothing here is pre-registered and nothing is a hypothesis test. R1 to R5 were
fixed before any number so the outcome could not be chosen afterwards, and
**none of them was a prediction.**

**THE GAP IS POST HOC.** It was derived after `E27` ran and was not among its
reading rules. It is labelled post hoc in the saved object, the verdict, both
figures and here.

Script: `scripts/E28_gap_readout.R`. **Sourced by the author 2026-09-11**; the
author's run reproduces **all 30 saved elements identically** against the dry
run, so every number below is read from the saved object. Design rationale:
`docs/2026-09-11_gap_readout_rationale.md`.

---

## 0. The answer in seven lines

| | |
|---|---|
| **Do the controls pass?** | **Both, exactly.** 88 positions reproduce `E27` at max delta 0, and the rebuilt `ox_rel` returns `E16`'s saved rho at delta 0 in both cohorts. |
| **Does differencing protect against composition?** | **NO.** The gap moves **2.11** percentile points under composition adjustment against **1.41** and **0.70** for the two positions it is built from. **The rationale's justification for the readout does not hold.** |
| **Is the mouse developmental gap resolvable?** | **No.** -42.3 [-80.3, +7.0] at n = 6. It cannot stand alone. |
| **Does development separate from MYC in the mouse?** | **Not established.** The point-estimate ratio of 3.3 to 8.6 has a separation interval of [-33.8, +66.9]. |
| **Is the human MYC difference a robust readout?** | **No. This is the substantive finding.** The sign depends on the estimator, and the **copy-number split disagrees in sign with the primary**. |
| **Is the gap a property of the compartment?** | **In the mouse by point estimate only** - all six comparators agree in sign but only **two of six** have intervals clear of zero. In the human contrasts only **four of six** agree in sign. |
| **Does the aggregate agree with `ox_rel`?** | **Yes, in direction, in both cohorts.** Corroboration on a second ruler and nothing more. |

**The one-sentence reading: the gap is a vivid way to draw a pattern that
already exists in `E25` and `E27`, and it is not a quantity that carries
evidence of its own.**

---

## 1. Both controls pass, and one of them was worth building

**Against `E27`:** 88 positions, eleven readouts across eight groups, max
|delta| = **0**.

**The rebuilt `ox_rel`:** rho against `ox_ppd` returns **0.9129712** in TCGA and
**0.8851938** in SCAN-B, **delta 0** against `E16`'s saved values, which are
read from `E16`'s object at run time rather than transcribed.

**Under which rule the rebuild sits, stated because it matters.** "Nothing is
re-scored" protects mitoPPS - the mouse artefact and the 142-pathway universe.
`ox_rel` is a gene-set composite rebuilt from `E16`'s own recipe on the human
linear matrices, which is a different object, and the control is what makes the
claim checkable rather than asserted. **Human side only**: the mouse artefact
carries pathway scores and not gene-level values, so a mouse `ox_rel` would mean
touching that repo.

---

## 2. The protection check fails, and it was the reason for the readout

The rationale justified the gap by saying two categories in the same ordering
move together under composition, so differencing cancels it. **That argument was
tested rather than accepted, and it does not hold.**

| quantity | maximum move under composition adjustment |
|---|---|
| OXPHOS position | 1.41 |
| central dogma position | 0.70 |
| **the gap** | **2.11** |

**Differencing does not cancel composition; it adds.** The two positions move in
opposite directions under adjustment, so their difference moves further than
either. The reason the original argument fails is structural: **a percentile is
already a within-ordering rank**, so a monotone panel-wide shift changes nothing
and differencing adds nothing there. The real threat is differential movement,
and differencing does not touch it.

**And the absolute scale, which that verdict alone does not carry.** The gap
moves 2.11 points against a gap difference of 15.5 on the same samples - about
**14%**. So differencing does not help, but composition is a small perturbation
either way. **Both statements belong in any use of this readout.**

This arm is **TCGA only**, because SCAN-B has no purity estimate and nothing is
imputed (trap 2). It runs on the 675 samples with complete adjusters rather than
the 730 in the two MYC groups, which is why its unadjusted gap difference is
15.5 and not section 4's 20.4. That is the coverage effect `E26` recorded, not
an adjustment effect.

---

## 3. The mouse side: neither R1 nor R2 is resolvable

| contrast | kind | change | 95% CI | clear of zero |
|---|---|---|---|---|
| WT window 6W to 12W | TIMELINE, batch-confounded | -42.3 | -80.3 to +7.0 | no |
| Myc+ window 6W to 12W | TIMELINE, batch-confounded | -34.5 | -73.2 to +11.3 | no |
| MYC at 6W | genotype, within timepoint | -12.7 | -57.0 to +36.6 | no |
| MYC at 12W | genotype, within timepoint | -4.9 | -43.7 to +35.2 | no |
| **TCGA MYC-high minus low** | MYC activity | **-20.4** | **-32.4 to -2.1** | **yes** |
| **SCAN-B MYC-high minus low** | MYC activity | **-14.1** | **-23.2 to -6.3** | **yes** |

**R1: covers zero.** This does not refute the point estimate. It means the mouse
gap contrast **cannot stand alone** and must be read beside the
expression-matched-null analysis in `myc_mouse`, **which remains the primary
evidence for the developmental change.**

**R2: the separation is not established.** The point estimates give a ratio of
3.3 at 6W and 8.6 at 12W, but the paired separation intervals are [-33.8, +66.9]
and [-26.1, +69.0]. **The ratio is recorded as a point estimate and must not be
quoted as a result.**

Of the eight per-group gaps, only the two human MYC-high groups have intervals
clear of zero. **Every mouse group gap covers zero at n = 6.**

The two temporal contrasts are batch-confounded, batch equals timepoint. That
ruling stands: they are reported in full and nothing is withdrawn on their
account.

---

## 4. R3 - the human MYC difference is instrument-dependent, and this is the finding

| instrument | proliferation entanglement | TCGA | SCAN-B |
|---|---|---|---|
| `MYC_UP.V1_UP` | 1.5% | **+14.8** | -2.1 |
| `FELSHER`, the primary | 14.8% | **-20.4** | **-14.1** |
| `HALLMARK_MYC_TARGETS_V1` | 23.5% | **+23.2** | **+21.8** |
| CollecTRI regulon | - | +11.3 | -6.3 |
| MYC mRNA | - | -57.0 | -54.2 |
| **TCGA copy number** | **not an expression programme** | **+25.4** | not available |

**Only the primary and MYC mRNA agree with the primary's sign in both cohorts.
`HALLMARK` gives the opposite sign in both. And the copy-number split, the one
instrument not defined on the expression matrix at all, disagrees in sign with
the primary.**

**The failure mode is not cleanly proliferation, and the note must not claim it
is.** Among the three instruments whose entanglement is measured, the gap
difference is **not monotone** in it: +14.8 at 1.5%, -20.4 at 14.8%, +23.2 at
23.5%. A proliferation artefact would not look like that.

**The strata say the same thing from the other side.** MYC was re-split inside
each proliferation tertile so the groups stay balanced.

| cohort | low proliferation | mid | high |
|---|---|---|---|
| TCGA | **+30.3** | **+16.9** | -12.7 |
| SCAN-B | -9.2 | -12.7 | **-45.8** |

**In TCGA the gap difference is positive at low and mid proliferation and turns
negative only at high. In SCAN-B it is negative throughout and grows with
proliferation. The two cohorts disagree.** If it were simply proliferation,
holding proliferation constant would remove it everywhere; instead it reverses
in one cohort.

**R3's second branch applies: the human gap difference is not a robust readout
of MYC.** What it is instead - instrument-dependent in a way this script does not
explain - is left open rather than named.

---

## 5. R4 - the comparator panel, and a distinction the point estimates hide

**Central dogma is a middling comparator, rank 3 of 6** by distance from OXPHOS.
Protein import is closest and Signaling is also closer. **The ordering is
identical on the value scale and the percentile scale**, which was checked
rather than assumed, so the rationale's conclusion survives its own mixing of
the two.

**In the mouse wild-type window all six comparators agree in sign - and only two
of six have intervals clear of zero.**

| comparator | change | 95% CI | clear of zero |
|---|---|---|---|
| Metabolism | -59.2 | -101.4 to -11.3 | **yes** |
| Small molecule transport | -51.4 | -88.7 to -0.7 | **yes** |
| Dynamics and surveillance | -46.5 | -95.1 to +5.6 | no |
| **central dogma** | **-42.3** | -80.3 to +7.0 | no |
| Signaling | -30.3 | -88.0 to +23.9 | no |
| Protein import | -22.5 | -40.1 to +9.9 | no |

So "a property of the compartment" is true **of the point estimates** in the
mouse and not of the intervals. **In the two human contrasts the reverse holds:
all six are clear of zero but only four of six agree in sign**, with Protein
import and central dogma negative and the other four positive. The human claim
is therefore narrower and depends on which comparator is chosen.

**The six gaps are reported individually and are never averaged into a
headline.**

---

## 6. The aggregate, and R5

**The aggregate is `ox_rel` rebuilt in mitoPPS space**, and the numerator matters.
With subunits, the analogue is close; with the umbrella it is diluted by **8.5
percentile points on average**, which is exactly what `E27` predicted since
assembly factors barely move while subunits do. **Subunits are the primary
numerator.**

**R5: the aggregate and `ox_rel` agree in direction in both cohorts** - the
aggregate rises by 38.4 and 44.1 percentile points where `ox_rel` rises by 0.212
and 0.253 in its own units. Only direction is compared, because the units differ.

**That is corroboration on a second ruler and nothing more.** `E16`'s R1 declared
the two rulers **not interchangeable**, and agreement on one contrast does not
overturn it.

---

## 7. C3 is uninformative, exactly as declared in advance

| quantity | unadjusted | adjusted |
|---|---|---|
| POSITIVE CONTROL: TCGA copy-number split | +25.4 | **+56.3** |
| the gap, TCGA | -20.4 | -2.1 |
| the gap, SCAN-B | -14.1 | **-26.8** |

**Under proliferation adjustment quantities grew rather than attenuated, the
positive control included.** An adjustment that doubles its own control is not
behaving as an attenuation, and nothing is read from it. **The retained ratio is
deliberately not leaned on**: above 1 means the quantity grew, and the unadjusted
values are small enough that the ratio is unstable.

This is the ambiguity the header declared before any number, not a surprise.
**C1 and C2 are primary because they remove nothing; C3 is secondary and
ambiguous by construction.** The pre-declared reading - C3 collapsing while C1
and C2 survive is over-adjustment rather than refutation - does not arise here,
because C1 and C2 did not survive either.

**The positive control's limitation, stated wherever it is used:** the
copy-number split is independent of the expression matrix **in measurement, not
in biology**. Amplified tumours proliferate more. The mitoribosome was not used
as a control because it sits inside central dogma and is one arm of the primary
gap, which would be circular.

---

## 8. What this carries, and what it does not

**The sentence it can carry, and it is a modest one:**

> Across the mouse developmental window the respiratory chain loses position
> relative to every other level-1 category, by point estimate, with two of six
> comparisons resolvable at n = 6.

**Everything below this line it cannot carry.**

1. **It carries nothing on the human side.** The sign of the human gap
   difference depends on which MYC estimator defines the groups, and the one
   DNA-level instrument disagrees with the primary. **No human gap number from
   this script should appear in a manuscript.**
2. **It is not independent evidence for the developmental change.** The
   expression-matched-null analysis in `myc_mouse` remains the primary evidence.
   This is a second presentation whose interval covers zero.
3. **It is not protected against composition.** Its own justification failed.
   Composition is a small perturbation, about 14%, but differencing does not
   reduce it.
4. **It is not a second line of evidence via the aggregate.** That is `ox_rel`
   on a different ruler, and the two rulers were already declared not
   interchangeable.
5. **The ratio of developmental to MYC effect must not be quoted.** The
   separation is not established.
6. **It says nothing about whether human tumours arrived at any configuration by
   a developmental route.** That was never in scope and no reweighting of these
   numbers reaches it.

**What survives from Phase 2 is unaffected.** `E25`'s group-ordering result and
`E26`'s composition check stand on their own estimands and are not built on the
gap. Nothing in `docs/2026-09-11_phase2_conclusion.md` changes.

---

## 9. Two corrections to the rationale, recorded rather than absorbed

**A rounding slip.** The group table gives `mouse 6W WT` as +24.7 where the
computed value is 24.6479, a delta of 0.052 that sits just outside the
one-decimal boundary. Every other cell reproduces. The assertion caught it
rather than tolerating it.

**Mixed scales.** Section 3's distance table is computed on mitoPPS values while
the gap is on percentiles. Both are now reported and **the comparator ordering
is identical on the two scales**, so the conclusion survives - but it survives by
luck rather than by construction, and the check is now in the script.

**And two definitions that must not be conflated.** The readout is a
**difference of group gaps**. The gap computed inside a contrast's own ordering
is a different statistic: **-42.3 against -28.9** for the wild-type window, and
it even changes sign for the two genotype contrasts. Both are carried in the
object so the second is never quoted by mistake.

---

## 10. Constraints observed

- **The gap is post hoc** and labelled so everywhere it appears.
- **Nothing is re-scored.** `ox_rel` is a composite rebuilt from its own recipe,
  with the control reported first, human side only.
- **No mitoPPS value crosses the species boundary.** The gap is in percentile
  units, which is what lets mouse and human sit in one table.
- **TCGA and SCAN-B separate**, never averaged.
- **`E25`, `E26` and `E27` are not modified.**
- **No ortholog call.** Tripwire clean.
- The mouse timeline is **batch-confounded**; nothing withdrawn on that account.
- **Purity is not addressed** beyond the protection check, which is TCGA-only.
- **N3.** Transcript associations. "Primed" appears nowhere.

## 11. Provenance

- Script: `scripts/E28_gap_readout.R`, seed 1, 5,000 paired bootstrap draws.
- Reads: `results/mitopps_group_rankings.rds`, `results/purity_and_composition.rds`,
  `results/category_readouts.rds`, `results/respiratory_rulers.rds`,
  `data/from_myc_mouse/mitopps_scores.rds` (md5
  `b8a125af4bc0d5909ad03f6e126e2890`), the validation snapshots and
  `results/set_definitions.rds`.
- Saved: `results/gap_readout.rds` (30 elements),
  `outputs/tables/E28_gap_readout.csv`,
  `outputs/mitopps_rank/E28_gap_by_group.pdf`,
  `outputs/mitopps_rank/E28_comparator_panel.pdf`.
