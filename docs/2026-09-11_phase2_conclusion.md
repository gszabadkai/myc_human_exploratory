# Phase 2 conclusion - human breast tumours match the mouse MYC move, not the mouse developmental move

**2026-09-11. EXPLORATORY, POST-HOC, DESCRIPTIVE. Nothing in Phase 2 is
pre-registered, nothing is a hypothesis test, and every statement below is
hypothesis-generating.**

This closes the cross-species mitoPPS rank comparison. Three scripts, `E24`,
`E25` and `E26`, with `E24` superseded in part by `E25` and `E26` testing
whether `E25`'s answer is composition.

**Read this, then `docs/2026-09-11_e25_group_ranking_comparison.md` for the
result and `docs/2026-09-11_e26_purity_and_composition.md` for the check.**

---

## 1. The question, and the answer

**The question.** Development lowers the priority of the respiratory chain
within the mitochondrial compartment in mouse mammary epithelium. Do human
breast tumours sit at that lowered position, or back at the pre-window one?

**The answer is that the question had a false premise, and correcting it gives a
cleaner result than the original question could have.** Human breast tumours do
not sit at a single position on the mouse's developmental axis. They sit at two,
and which one depends on MYC.

| human group | nearest mouse group | rho | 95% CI | runner-up rho | separated |
|---|---|---|---|---|---|
| TCGA MYC-low | mouse 12W WT | 0.61 | 0.47 to 0.73 | 0.08 | yes |
| SCAN-B MYC-low | mouse 12W WT | 0.64 | 0.52 to 0.75 | 0.16 | yes |
| TCGA MYC-high | mouse 6W Myc+ | 0.46 | 0.30 to 0.60 | 0.00 | yes |
| SCAN-B MYC-high | mouse 6W Myc+ | 0.51 | 0.36 to 0.64 | 0.09 | yes |

**Two of the four mouse groups match nothing human.** Every cross-species pair
involving mouse 6W wild-type or mouse 12W Myc+ has an interval covering zero.

**The two axes point opposite ways.** The human MYC contrast agrees with the
mouse genotype effect at both timepoints (rho 0.39 to 0.54, intervals clear of
zero) and disagrees with both timeline contrasts (-0.38 to -0.52). The
respiratory chain sits at the 76th and 79th percentile of the two mouse genotype
contrasts and the 10th and 9th of the two temporal ones.

**So: what human tumours share with the mouse is the MYC move, not the
developmental move.**

---

## 2. Why the first attempt got a different answer

`E24` pooled the human cohorts over MYC status and compared a within-sample
percentile. It concluded that TCGA sat between the two mouse positions and near
neither, that SCAN-B sat near the 6W position, and that the cohorts disagreed.
**Two things were wrong with it, in order of weight.**

**The estimand was the wrong one.** The mouse figures rank *pathways* by a
*group contrast*. `E24` ranked *one pathway* inside each *sample*. Those
statistics behave differently, and not symmetrically between the species:

| rho(statistic, log set size) | mouse | TCGA | SCAN-B |
|---|---|---|---|
| `E24`, within-sample percentile | +0.165 | +0.536 | +0.571 |
| `E25`, group contrast | +0.03 to -0.15 | +0.192 | +0.196 |

**`E24` compared two differently size-biased statistics across the species
boundary.** The cause is the rank transform, not mitoPPS: mitoPPS is centred at
exactly 1 per pathway in both species and both implementations take the mean per
set, so the first moment is size-free on both sides. What is not size-free is
the cross-sample variance, and a within-sample rank transform converts that into
a mean-rank difference that scales with cohort heterogeneity.

**The cohorts were pooled over MYC status.** Averaging two real and opposite
positions produced a middling one that read as nothing in particular.

**`E24`'s panel reconciliation stands and is reused unchanged**: 142 shared
pathways, the identical eight dropping under the gene-count filter in both
species, `OXPHOS subunits` at 89 in both, `mt-`/`MT-` in one synthetic pathway
of 13 in both.

---

## 3. What was checked, and what each check cost

**The panels are the same object.** The full MitoCarta hierarchy path agrees on
all 141 shared catalogue pathways, 7 level-1, 38 level-2 and 96 level-3 in both
species, each read from its own native source. The derived 114-member antichain
is identical from both. **No ortholog function is called anywhere.**

**The mouse side reproduces.** Recomputing all four contrasts from the snapshot
reproduces the artefact's own pairwise table at max |delta| 4.4e-16.

**Set size does not drive it.** On the corrected estimand nothing exceeds 0.21
against log set size, where `E24`'s statistic was 0.54 and 0.57.

**Hierarchy does not drive it.** The non-overlapping antichain moves the primary
position by one to two percentile points.

**Composition does not drive it.** `E26`'s full result is in its own note. The
gate that licenses reading it: proliferation keeps **96%** of its MYC difference
under purity and leukocyte adjustment while those adjusters explain **8.8%** of
its variance, which is the mouse fat-pad diagnostic coming out the other way.
The separation goes 44.4 unadjusted to 43.0 and 47.2 adjusted; the respiratory
chain attenuates to the **41st percentile** of the panel's own attenuation
distribution, not an outlier; the cross-species agreement survives every model;
the high-purity subset gives 41.5 points at n=235; and the copy-number split
**widens** under adjustment, 13.7 against 86.3 becoming 11.6 against 88.4.

**The estimator does matter, and that is trap 3 and trap 4 working.** Four MYC
activity estimators give the MYC-high minus MYC-low shift as +28 to +86
percentile points. **MYC mRNA gives -5.6 and -33.8** - a sign flip. No statement
here rests on one signature, and none rests on MYC mRNA.

---

## 4. What is not settled

**SCAN-B cannot be composition-checked.** It carries no purity estimate and
nothing is imputed. Its numbers stay unadjusted. The reassurance is indirect:
the two cohorts agreed at +0.88 in `E25` and the checkable one survived. That is
not the same as SCAN-B surviving.

**The mouse timeline is batch-confounded**, batch equals timepoint. The author
ruled that the timeline comparisons stay in, on the grounds that
reprioritisation is a timeline effect measurable only that way and that DESeq2
normalisation is what handles the batch. Every temporal quantity is reported in
full with the confounder stated and nothing withdrawn on its account. **The
genotype contrast within a timepoint remains the cleaner axis, and the
conclusion above rests on that one.**

**The mouse arms are n = 6.** Every mouse interval is wide. The human intervals
are narrow because n is 365 and 1,069 per group, so the human positions are
precise and the mouse reference positions are not.

**The adipocyte tracer is one instrument** and cannot separate "less adipose
tissue in the section" from "tumour cells expressing fewer adipocyte genes". A
deconvolution-based estimate would be an independent second opinion and is not
run.

**This is a transcript ordering and nothing more.** No flux, no respiration, no
priming. N3 throughout.

---

## 5. What would falsify it

1. **A human cohort whose MYC-high ordering fails to match mouse 6W Myc+.**
   That would make the mapping cohort-specific rather than general.
2. **A mouse experiment with more animals in which 12W WT stops being the
   near-mirror of 6W Myc+** (rho -0.95 here). That structure is what makes the
   human split interpretable at all.
3. **A deconvolution-based composition estimate that collapses the separation
   while proliferation survives.** `E26` rules this out for purity, leukocyte
   fraction and an adipocyte marker score; it does not rule it out for every
   possible composition instrument.
4. **A demonstration that the group-contrast estimand carries a bias that the
   size and hierarchy checks missed.** The estimand was chosen to match the
   mouse figures, and its size-freedom is measured rather than assumed, but it
   has not been stress-tested the way `E24`'s was.

---

## 6. The one sentence, and what it is worth

**MYC-high human breast tumours order the mitochondrial compartment the way the
MYC-driven mouse gland does, MYC-low tumours order it the way the older
wild-type gland does, and the developmental axis that motivated the comparison
is the one axis human tumours do not match.**

**It is one descriptive observation, generated post hoc, in an openly
exploratory repo.** It is not a test, it did not have a pre-registered
prediction, and both possible outcomes were written down before the numbers were
seen because neither was predicted. What it is good for is handing a confirmatory
study a hypothesis that is sharp enough to fail: **the mitochondrial
reprioritisation shared between mouse and human is MYC's, not development's.**

---

## 7. The record

| | |
|---|---|
| `E24` | `scripts/E24_mitopps_rank_comparison.R`, `docs/2026-09-10_e24_mitopps_rank_comparison.md` - **superseded in part**, banner names the sections |
| `E25` | `scripts/E25_group_ranking_comparison.R`, `docs/2026-09-11_e25_group_ranking_comparison.md` - **the result** |
| `E26` | `scripts/E26_purity_and_composition.R`, `docs/2026-09-11_e26_purity_and_composition.md` - **the composition check** |
| inputs | `data/from_myc_mouse/` (read-only artefact, md5 pinned), `data/from_validation/`, `results/scanb_scores.rds`, `results/new_set_scores.rds`, `data/mitocarta_human/` |
| standing | mouse repo never written to, working tree and HEAD `129b0c0` as found; ortholog tripwire clean; rank only across the species boundary; TCGA and SCAN-B never averaged |
