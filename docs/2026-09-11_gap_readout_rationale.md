# The OXPHOS gap readout — where we stand, and what E28 tests

**2026-09-11.** Background for the Claude Code session that builds `E28`. Read
this before the prompt.

---

## 1. Where Phase 2 stands

`E24` compared a within-sample percentile and was superseded — wrong estimand,
size-biased asymmetrically across species. `E25` compared **group orderings** and
found human tumours split on MYC: MYC-low tumours order the mitochondrial
compartment like a mature wild-type mouse gland, MYC-high like a young
MYC-driven one. `E26` showed that split is not composition: the positive control
kept 96% of its MYC difference while the adjusters explained 8.8% of its
variance, and the copy-number split widened under adjustment.

`E27` read eleven readouts instead of one, on the native MitoCarta level-1
partition. Three results from it matter here:

- **Assembly factors barely move anywhere.** Mean absolute deviation from
  baseline: human subunits 0.050 against assembly 0.0095; mouse 0.066 against
  0.021. The respiratory chain moves; the machinery that builds it does not.
- **The OXPHOS umbrella is exactly its two halves** (89 + 68 − 1 = 156), so it
  carries no independent information and is retired as a readout.
- **The ordering of OXPHOS and central dogma flips** between mouse 6W Myc+ and
  human MYC-high, in both cohorts.

---

## 2. The quantity this is about

That last result generalises into a single number: **OXPHOS position minus
central dogma position, within the same ordering.**

| group | OXPHOS | central dogma | **gap** |
|---|---|---|---|
| mouse 6W WT | 60.2 | 35.6 | **+24.7** |
| mouse 6W Myc+ | 82.7 | 70.8 | **+12.0** |
| TCGA MYC-low | 25.0 | 22.2 | +2.8 |
| SCAN-B MYC-low | 19.4 | 20.1 | −0.7 |
| SCAN-B MYC-high | 67.3 | 82.0 | **−14.8** |
| TCGA MYC-high | 62.3 | 79.9 | **−17.6** |
| mouse 12W WT | 13.7 | 31.3 | **−17.6** |
| mouse 12W Myc+ | 34.2 | 56.7 | **−22.5** |

**In the mouse it falls across the developmental window in both genotypes, by
three to eight times more than MYC moves it at either timepoint:**

| contrast | change in gap |
|---|---|
| WT window, 6W → 12W | **−42.3** |
| Myc+ window, 6W → 12W | **−34.5** |
| MYC effect at 6W | −12.7 |
| MYC effect at 12W | −4.9 |
| TCGA MYC-high − low | −20.4 |
| SCAN-B MYC-high − low | −14.1 |

**Why it is worth having.** It states Movement 1's claim in one number on the
native ontology, without expression-matched nulls. And it is better protected
than a single position: both categories sit in the same ordering, so composition
that shifts the whole panel shifts them together and largely cancels — which is
the objection that sank the cross-species *level* comparison in `E24`.

**Its limitation, and it is the reason for test 2 below.** The gap also falls
under MYC in the mouse. So a low gap in human MYC-high tumours is predicted by
the MYC axis on its own, and the human data alone cannot separate "these tumours
sit where development put the mouse gland" from "these tumours are MYC-high and
MYC lowers the gap". The magnitudes lean the right way — the human MYC effect,
−14 to −20, is larger than any mouse MYC effect and closer to the developmental
scale — but that is a magnitude comparison in percentile space across species,
which is soft.

**The gap was derived after `E27` ran.** It was not among R1–R5. Reportable, but
labelled post hoc throughout.

---

## 3. The three tests

### Test 1 — intervals

Everything above is eight point estimates with no uncertainty, and the mouse ones
come from n = 6. Bootstrap the gap per group and per contrast, resampling
**samples within group** and recomputing both positions and their difference, so
the pairing is preserved.

**Expect the mouse intervals to be wide, and pre-declare what that means.**
`E27`'s paired statistic already showed this: the |ΔOXPHOS| − |Δcentral dogma|
difference was 0.082 [−0.019, 0.158] and covered zero, where the fixed-bound
version passed. The paired number is the honest one.

An interval covering zero at n = 6 does not refute the point estimate; it means
the mouse gap contrast **cannot stand alone** and must be read beside the
expression-matched-null analysis in `myc_mouse`, which is where the developmental
change was established and remains the primary evidence. If the interval is
clear of zero, the gap is a second, independent presentation of it.

### Test 2 — the proliferation confound

MYC dose tracks proliferation monotonically, and proliferation drives
mitochondrial biogenesis, so a MYC-high group could show a low gap for reasons
that have nothing to do with the respiratory chain.

**The honest difficulty: proliferation and MYC are entangled biologically, not
merely by measurement.** Adjusting the gap on proliferation risks removing MYC
itself, which makes a collapse under adjustment ambiguous rather than
informative. The three approaches are therefore not equivalent:

**P1 — the estimator ladder (primary).** `E25` carries five MYC estimators, and
their MYC-high minus MYC-low shift ranges from +28 to +86 percentile points, with
**MYC mRNA giving −5.6 and −33.8**, a sign flip. Recompute the gap difference on
each. If it survives on the least proliferation-entangled readouts, it is not
simply proliferation. This is the cleanest of the three because it removes
nothing.

**P2 — stratification (primary).** Within a narrow proliferation band, does
MYC-high still show a lower gap? This holds proliferation roughly constant
without subtracting anything from MYC. Report stratum sizes; they will be small.

**P3 — adjustment (secondary, and ambiguous by construction).** `E26`'s
machinery exists and its GATE 0 proved intercept-only residuals reproduce `E25`
exactly. But trap 15 applies with a twist: **the usual positive control is
proliferation, and here proliferation is the adjuster.** A different control is
needed, and every candidate MYC-driven quantity is itself proliferation-
entangled. Whatever is chosen, its non-independence must be stated, and a
collapse under P3 with P1 and P2 surviving is read as over-adjustment, not as
refutation.

### Test 3 — the other categories

Central dogma is one comparator of six. Two things need doing.

**Report all six gaps separately.** From `E27`'s mouse WT window values, distance
from OXPHOS's −0.0926:

| comparator | its own change | distance from OXPHOS |
|---|---|---|
| Protein import, sorting and homeostasis | −0.0555 | **0.037 — closest** |
| Signaling | −0.0318 | 0.061 |
| **Mitochondrial central dogma** | −0.0108 | **0.082** |
| Mitochondrial dynamics and surveillance | +0.0077 | 0.100 |
| Small molecule transport | +0.0118 | 0.104 |
| Metabolism | +0.0313 | **0.124 — furthest** |

**Note the correction:** Protein import is closest to OXPHOS as expected, but
**Signaling is also closer than central dogma**. So central dogma sits third of
six, not second. That makes it a **middling** comparator rather than a
conservative one, and the six gaps will span a wide range. They are variable and
should be reported individually, never averaged into a single headline.

**And the aggregate is not an independent result.** OXPHOS against the mean of
the other six is, conceptually, `ox_rel` — the respiratory chain against the rest
of MitoCarta — rebuilt in mitoPPS space. If the two agree, that is corroboration
on a second ruler and worth stating. **It is not a second line of evidence**, and
must not be presented as one.

---

## 4. What E28 can and cannot settle

**Can:** whether the gap's developmental fall in the mouse is resolvable at
n = 6; whether the human MYC difference survives the proliferation ladder and
stratification; whether the pattern is a property of the compartment or of the
central dogma pairing specifically.

**Cannot:** whether human tumours arrived at the post-window configuration by a
developmental route. That is a cross-sectional comparison across species and
states, and no reweighting of these numbers resolves it.

**Standing:** exploratory, post-hoc, descriptive. The reading rules exist to stop
an outcome being chosen after the fact, not because any outcome is expected.
N3 throughout.
