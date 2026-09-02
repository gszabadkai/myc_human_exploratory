---
date: 2026-09-02 (evening)
status: disposable - delete once used
purpose: ONE script to re-run. The science has moved a long way; the paper
         framing is written and lives in three notes.
supersedes: docs/2026-09-02_handoff.md (deleted)
---

# Handoff - one script to re-run, and the paper has a spine now

`main` level with `origin`, tree clean. **`E00`-`E13` are all written and all
have been run.** Only `E10` is out of date.

---

## Resume prompt

```
Read CLAUDE.md, then docs/2026-09-02_handoff_evening.md.

This repo is EXPLORATORY. Nothing is pre-registered. The sibling repo
myc_human_validation is a COMPLETED pre-registered study, frozen at d3ac60e,
not reopened, and its no-post-hoc rule does NOT apply here.

Everything is run except E10, which needs re-sourcing because it was made
proliferation-adjusted after my last run. I source scripts in Positron; you
write and edit them, you do not run them.
```

---

## 1. RUN THIS

| script | why | time |
|---|---|---|
| **`E10_machinery_measures_and_priming.R`** | **REQUIRED.** Your run is unadjusted; it is now partial Spearman on `PROLIF_DISJOINT` so its panels match E11 and E13. Numbers change | ~15 s |
| **`E14_curated_comparators.R`** | **NEW, 2026-09-02 evening.** Closes open item 1 and revises the paper's last clause - see `docs/2026-09-02_e14_curated_comparators.md`. Dry-run verified end to end | ~4 min |
| `E13_priming_and_content.R` | optional. Only two panel subtitles changed (the localisation wording); no number moves | ~1 min |

`E09`, `E11`, `E12` are current - do not re-run.

**Check as E10 scrolls past:** `max |difference| = 0e+00` (E10 still reproduces
E08 through a separately computed unadjusted pass), `31 anchor pathways ... all
present`, and `n_gain_in_both` = **0 on OXPHOS, 5 on MYC**. If the last one
still says 6 and 3, the adjustment did not take.

**Check as E14 scrolls past:** the eligibility table drops `protein import`,
`cristae formation` and `mitochondrial Ca transport` for having no cytosolic
half; `VERDICT AGAINST THE FALSIFIER` reads `APOPTOSIS-SPECIFIC in both
cohorts`; and `VERDICT ON THE CONFOUND` reads `SURVIVES`. If the second says
`INFILTRATE`, stop and read section C3 of the E14 note before believing
anything downstream of it.

---

## 2. THE PAPER. Start with the synthesis, then the three notes

**`docs/2026-09-03_human_arm_for_mouse_reconciliation.md` is the entry point.**
It summarises all of the below into one claim ladder, states what goes in the
paper, lists the standing negatives, and sets out what the mouse arm is being
asked to reconcile against. Section 8 of it is a superseded-statements table -
where it and any note below disagree, it is current.

### The notes it summarises, read in this order

| note | what it is |
|---|---|
| `2026-09-02_paper_opening_human.md` | **the main one.** The claim ladder - which of the five component claims each result licenses - the four-sentence opening, the figure, the statistics, and drafted text for the priming subsection |
| `2026-09-02_priming_interaction_tested.md` | the MYC x OXPHOS interaction, tested and failed. Read before anyone re-asks it |
| `2026-09-02_e11_prolif_adjusted.md` | P1-P5 with the numbers |

### The finding, in one paragraph

**Read `docs/2026-09-02_e14_curated_comparators.md` C5 first - it replaces the
last clause of this paragraph.**

In two cohorts the canonical apoptotic machinery is ordered along OXPHOS and not
along MYC. The two axes are correlated (0.26-0.32 adjusted) and rank the 44
genes almost identically (0.61-0.73), so the separate comparison could not have
settled which carries it - **conditioning does**: removing MYC leaves the OXPHOS
ordering intact (0.453 -> 0.485), removing OXPHOS abolishes MYC's
(0.187 -> -0.043). Neither is proliferation, purity or infiltrate. What predicts
a gene's position is whether it belongs to the nuclear-encoded mitochondrial
regulon, not whether it promotes or prevents death. **And the ordering is no
steeper than an expression- and compartment-matched gene set gives**, so it is
real, it is OXPHOS's, and it is not specific to apoptosis.

### The figure

`E11_fig9_paper_figure1` - A the plane, B the conditioning ladder, both with
bootstrap intervals printed. **B is not optional**, and `E11_fig2` (the
mitoribosome control) must be cited in its legend. Delete the on-figure title
for submission.

---

## 3. What each script now contributes

| | |
|---|---|
| `E09` | Spearman vindicated. The measures agree at 0.996 over 220 pairs and bicor sides with Spearman in all 12 largest gaps. **The one place measure matters is mitoPPS** (mean gap 0.029, up to 0.093) - never report a Pearson against a mitoPPS score. Nothing is non-monotone |
| `E10` | the 44 re-annotated from Reactome and MitoCarta; both axes; the 35 priming ratios and their luminal/basal split |
| `E11` | the conditioning test (4.3), the bootstrap intervals (4.4), three progressively stricter composition nulls, the compartment split |
| `E12` | the borrowing decomposition, made legible for collaborators |
| `E13` | the content control the 44-gene analysis could not run |
| `check_sandboxes.R` | maintenance, not pipeline. Run after renaming any saved column |

### Findings that changed today

- **E10 R1** (localisation split 0.453, replicated under MitoCarta) - **still
  true, and reframed twice.** First, it is inside a composition-matched null, so
  not apoptosis-specific. Second, MitoCarta is a proteome catalogue and half the
  BCL2 family translocates, so at transcript level it marks **mitochondrial
  regulon membership**, not localisation. The reframe is stronger: it explains
  why the null came out flat.
- **E10 R4** withdrawn - superseded by E09.
- **The priming ratios** are additive in their components: an additive model of
  numerator and denominator explains **92-95%** of all 35 values. After
  adjustment **none** of the 35 beats both its components on OXPHOS.
- **`BID` and `PMAIP1`/NOXA** are the only priming genes whose MYC association
  survives conditioning on OXPHOS in both cohorts under all three estimators
  (+0.12 to +0.28 and -0.05 to -0.21). **`BCL2` does not** - it is
  estimator-dependent, so "MYC represses BCL2" must not be written.

Phase 1 and earlier phase 2 findings (F1-F7, D0-D5, E1-E4, M1-M4, S1-S6) stand
as recorded in their own dated notes. F1 resolved, F3 reframed, D1 reinterpreted,
D3 qualified, S6 superseded by P4/P4a.

---

## 4. Open

1. ~~**The ~1.3-1.6 SD residue** above the compartment-matched null.~~
   **ANSWERED by `E14`, 2026-09-02.** The comparator is mitophagy, not protein
   import - 62 of the 62 import genes are in MitoCarta, so it has no cytosolic
   half and the split statistic does not exist for it. The machinery
   out-splits every programme tested, and the reason is its CYTOSOLIC half,
   the only one in the comparison that runs against OXPHOS. Not infiltrate.
   See `docs/2026-09-02_e14_curated_comparators.md`, and note that its **C5
   revises the last clause of the finding paragraph above**. The successor
   open item is C6.1: a cytosolic stress programme that is not apoptotic, as
   the next comparator.
2. **Mediation versus confounding is not identifiable here.** MYC -> OXPHOS ->
   genes and a common cause give identical partial correlations. What IS ruled
   out is OXPHOS acting through MYC (conditioning on MYC leaves OXPHOS intact).
   The perturbation that separates them is the mouse arm's.
3. `LumA` alone has not been run; `Basal` rests on 171 TCGA samples.
   **Promoted** - E14 C6.3 makes it the falsifier for the new claim.
4. **Score `COLLECTRI_MYC_STIM`** (739 genes, in the snapshot, never scored) as
   a fifth base in `E02` - needs a pipeline re-run.
5. **Drop `ELLWOOD`** from the panel, recompute F1's entanglement slope - same.
6. The ER-negative fatty-acid-oxidation reversal, untouched.
7. A stranded total-RNA dataset for `MT-ND6` and `CO1`/`CO2`.

**Out of scope until reopened:** MCbiclust/forkscale, survival, treatment,
METABRIC, DepMap, causal modelling.

---

## 5. Traps

Carried forward: stale saved artefacts; `slice(-seq_len(0))` returns no rows;
a message that disagrees with its code; prose that prejudges its own test;
`frac_prolif` is 0 in stripped variants; `thin` must be per cohort.

Added since:

1. **`tibble()` evaluates in sequence** - a new column shadows an outer object
   of the same name. A column `mitocarta` beside a data frame `mitocarta` read
   the logical vector.
2. **`pivot_longer(names_sep = "_")` splits on EVERY underscore.** `adj.
   PROLIF_DISJOINT` produced an `NA` facet. Use `names_pattern`.
3. **ggplot CLIPS captions, it does not wrap.** Hard-wrap at ~85 characters for
   9 inches, ~125 for 13. `paste0` adds no space, so one rendered line can span
   two source strings.
4. **Subsetting a named vector by an absent name returns an element named
   `<NA>`.** `split()` then drops the group silently and `cor()` returns NA with
   a warning and no error - eight cells read NA and the table still printed.
5. **`if (FALSE)` sandboxes are the one thing nothing ever runs.** Use
   `scripts/check_sandboxes.R` after renaming any saved column.
6. **An interval must be in the units of the axis it is printed beside.** A
   rank-split CI was printed on a median-difference axis: 0.53 beside points
   whose difference is 0.40.
7. **Two panels in one paper cannot be on different footings.** E10 was raw
   while E11 and E13 were adjusted; caught only when drafting the text.
8. **Choosing a statistic after seeing its z is statistic-shopping.** The
   median difference is plotted for legibility and the rank split still anchors
   the claim, because the median difference happens to have the larger z. Both
   are carried in every table.
9. **Check a drawn claim against the drawing.** "The matrix is organised by its
   denominators" was asserted from eyeballing a heatmap and was wrong; the
   decomposition gave 92-95% additive instead.

---

## 6. Rules

- **Option A.** Claude Code writes and edits the numbered scripts and does not
  run them. Infrastructure it may execute.
- **Exploratory.** Nothing is pre-registered. Report structure, gradients and
  cross-cohort reproducibility. Write the falsifier BEFORE the next analysis -
  `E11` sections 3.1/3.2 and `E13` section 3 are what that looks like.
- **The MYC naming contract.** `__FULL` / `__MITOSTRIP` / `__PROLIFSTRIP` /
  `__BOTHSTRIP`, no bare names, use the `E00` constants.
- **Never write to `myc_human_validation` or `myc_mouse`.**
- **Human only.** Ortholog tripwire narrowed to calls.
- **Scale discipline.** GSVA log VST, mitoPPS linear DESeq2, never one object.
  Rank correlations are immune; E09 showed Pearson is not, on mitoPPS.
- **Never pool GSVA or mitoPPS values across cohorts.**
- R: no `print(n = X)` after `head()`; `dplyr::count()`; ASCII-only strings;
  every numbered script ends with an `if (FALSE)` sandbox.
- **Gene lists selected on a statistic are descriptions, not findings.**
