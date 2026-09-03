---
date: 2026-09-03 (evening)
status: disposable - delete once used
purpose: ONE short script to re-run, and where the science stands after a day
         spent on the priming side.
supersedes: docs/2026-09-02_handoff_evening.md - which is KEPT, because the
         synthesis document cites it, but a banner has been put on it. Do not
         act on its "RUN THIS" table; both scripts in it have been run.
---

> **SUPERSEDED 2026-09-04 by `docs/2026-09-04_handoff.md`.** Everything still
> live here is carried there, including the same two-script "RUN THIS" table.
> **Do not act on section 6** - it says 15 commits are unpushed, which was true
> when written and is not now; `origin/main` is at `17c09ca`. Kept only because
> `docs/2026-09-02_handoff_evening.md` points at it.

# Handoff - two short re-runs, and the model is written

`main`, tree clean. `E00`-`E15` are all written. **`E15` and `E10` are both a
re-source behind**, and in both cases no existing number changes.

**Updated 2026-09-03, later:** `E10` gained script section **5.0a-ii** in answer
to "the OXPHOS row of the gain table is 0 in both cohorts - is a gain in one
cohort worth a mention?" It is; sections 3.3a and 3.3b of the synthesis document
are the write-up, verified against a redirected dry run that wrote nothing into
the repo. Re-source `E10` to put the three new objects on disk.

---

## Resume prompt

```
Read CLAUDE.md, then docs/2026-09-03_handoff_evening.md, then
docs/2026-09-03_human_arm_for_mouse_reconciliation.md.

This repo is EXPLORATORY. Nothing is pre-registered. The sibling repo
myc_human_validation is a COMPLETED pre-registered study, frozen at d3ac60e,
not reopened, and its no-post-hoc rule does NOT apply here.

Everything is run except E15 and E10, both of which need re-sourcing and
neither of which changes an existing number. I source scripts in Positron; you
write and edit them, you do not run them.

The next step is confronting the human model in section 5 of the synthesis
document with a mouse model built from experiments. Details of the human data
matter for that; section 2 of the same document is the inventory.
```

---

## 1. RUN THIS

| script | why | time |
|---|---|---|
| **`E15_two_axis_gene_view.R`** | figure 5 was redrawn after your 14:51 run: both cohorts in **one** panel, no droplines to zero, narrower, gene names coloured by MitoCarta membership. **No number changes** - it is a display script and recomputes nothing | ~10 s |
| **`E10_machinery_measures_and_priming.R`** | new section **5.0a-ii**: what the OXPHOS `gain > 0 in both` = 0 does and does not mean. Adds `$gain_rank_agree`, `$gain_by_discord`, `$gain_one_cohort` and a `discord` column on `$priming` / `$priming_strata`. **Every existing object is bit-identical** - checked object by object against the 17:00 run | ~20 s |

**Check as it scrolls past:**
`anchor OK: the split recomputed here is bit-equal to E14's for all 8 cohort x
axis x adjustment cells`. If that line does not appear, E11 or E14 has moved
underneath E15 and nothing downstream is comparable.

**`scripts/check_sandboxes.R` reports `E10 ... FAILED` until you re-source it**,
and that is the check working, not breaking: E10's sandbox now names
`x$gain_one_cohort`, which the saved object will not carry until the re-run. The
three new sandbox lines were verified against the dry-run object. After the
re-source it goes back to `OK`; if it does not, something else moved.

From `E10`, three lines worth reading as they scroll: the `gain` rank agreement
between cohorts (**+0.408** OXPHOS, +0.560 MYC), the discordance table, and
`positive gains outside the discordant block: 2 of 29`.

**Everything else is current. Do not re-run E09, E11-E14.** `E11`, `E12`, `E13`,
`E14` are unchanged since 2026-09-02.

---

## 2. WHERE TO START READING

**`docs/2026-09-03_human_arm_for_mouse_reconciliation.md` is the entry point and
it was rewritten today.** It is now organised around **apoptotic priming
regulation by MYC and OXPHOS**, because that is what the mouse comparison turns
on. Four things in it are worth knowing before opening anything else:

| section | what it is |
|---|---|
| **2** | **the data-source inventory.** Cohorts, how each axis is defined and on what scale, the covariate and why it is disjoint from both axes, the priming lists and the ratio definition, pinned catalogue SHAs, and which saved object and column to read first. Written so a mouse model can be applied to human variables whose construction is known |
| **3** | **the priming analyses**, in the order they answer one another. 3.1 is the 12-transcript table and is the table to map a mouse model onto |
| **4** | the compartment split, condensed. Unchanged science, demoted to context. **4.4 is the hinge** - 11 of the 12 priming genes are on the same side of the compartment rule, so it cannot order them |
| **5** | **the model.** Six labelled arrows, five propositions M1-M5, what the model forbids, seven predictions V1-V7 |

Per-script notes, unchanged unless marked:
`e10_machinery_and_priming.md` (**R5c and R5d added today**),
`e11_prolif_adjusted.md`, `e14_curated_comparators.md`,
`e15_two_axis_gene_view.md` (**new today**),
`priming_interaction_tested.md`, `paper_opening_human.md`.

---

## 3. WHAT CHANGED TODAY

| | |
|---|---|
| **`E15` is new** | the 44 machinery genes on BOTH axes at once, which no previous figure showed. Five figures: a per-gene dumbbell sorted by the OXPHOS-minus-MYC difference, a numbered heatmap, the scale-artefact check, a two-bar mean +/- SD summary, and the difference alone with both cohorts in one panel |
| **`E10` gained the cell marks** | figures 3 and 6 now carry figure 4's test directly: `*` where a ratio beats its own two genes and reaches `\|rho\| >= 0.30`, `**` where both its genes are individually OXPHOS-led, `^` on the axis labels of those genes, a heavy border where the `*` test passes in **both** cohorts. **Not one MYC cell is marked anywhere.** Three cells pass in both cohorts and all three are Basal OXPHOS |
| **`E10` gained figure 8** | the two axes collapsed into one number per ratio, drawn the **three** possible ways, because the choice is not cosmetic - the signed difference and the which-axis difference disagree in sign on 22 of 70 cells |
| **`E10` gained `additive_fit`** | the paper's 92-95% headline now regenerates from code instead of from a note |
| **the synthesis document was rewritten** | priming-first, with the data inventory and the model |
| **one figure-3 bug fixed** | its subtitle had been *clipped* rather than wrapped since the figure was written, silently losing its last three words on every run |

**One number in the whole document has no script behind it**: the cognate-pairing
test, synthesis section 3.8. It is a **negative** - the transcript configuration
is not organised by BH3 binding specificity - so it lowers a claim rather than
raising one, but it should still be folded into `E10`.

---

## 4. OPEN, in the order I would do them

1. **The mouse confrontation.** This is what the synthesis document was rewritten
   for. Section 5 is the human model; sections 8.2 and 8.3 say what would count
   as agreement and what would not. **8.3's "does NOT count as agreement" list is
   the one to read first** - four of the five tempting-but-empty comparisons are
   named there.
2. **The CollecTRI regulon analysis.** `data/collectri_human/` is pinned,
   snapshotted and has never been used for this. Two questions in one: what
   arrow (b) is - MYC's own, larger-than-borrowed, completely unsorted
   association with these transcripts - and what the cytosolic half's shared
   correlate is. **The single most valuable next human analysis**, and it has a
   falsifier already written (synthesis section 10.1).
3. **Fold the cognate-pairing test into `E10`** so section 3.8 regenerates.
   Small: filter `$priming` to `axis == "OXPHOS"` and the four selective
   sensitisers, group `gain` by cognate membership.
4. **A cytosolic stress programme that is not apoptotic**, as the next comparator
   for arrow (d). If a proteotoxic or integrated-stress-response module also runs
   against OXPHOS, the property is "cytosolic and stress-responsive" rather than
   "cytosolic and apoptotic".
5. **`LumA` alone** - the homogeneous stratum, and where prediction V6 should be
   checked from the other side.
6. Score `COLLECTRI_MYC_STIM` (739 genes, in the snapshot, never scored); drop
   `ELLWOOD`; recompute the entanglement slope. Needs a pipeline re-run.
7. The ER-negative fatty-acid-oxidation reversal, untouched.

**Out of scope until reopened:** MCbiclust / forkscale, survival, treatment,
METABRIC, DepMap, causal modelling on human data.

---

## 5. THINGS THAT WILL BITE

1. **No ortholog function anywhere, in either direction.** The reconciliation is
   a comparison of **conclusions**, not of gene lists. The tripwire is
   `grep -rnE "(mouse_to_human|human_to_mouse|ortholog[s]?)[[:space:]]*\(" scripts/`
   and it must stay empty.
2. **Never score SCAN-B without `scanb_pheno.rds$symbol_map`.** If a mouse model
   is projected onto SCAN-B by symbol, this is the first thing that breaks
   silently - coverage drops from 0.989 to 0.775 and Complex V loses its F1 head.
3. **Never pool GSVA or mitoPPS values across cohorts** - and a species is a
   cohort. Compare correlations, patterns and rankings.
4. **N2 outranks anything exploratory.** The sibling study is pre-registered,
   found the functional `MYC x OXPHOS` interaction on apoptotic priming null, and
   is frozen. Nothing here overturns it and nothing may be written as if it did.
5. **Priming is not measurable in transcript abundance.** Write "carries a higher
   `BAD`/`MCL1` transcript ratio", never "is more primed".
6. **Option A holds.** Claude Code writes and edits the numbered scripts; the
   author sources them in Positron. Infrastructure - git, snapshots, provenance
   READMEs, docs - Claude Code may do directly.

---

## 6. Git

15 commits ahead of `origin/main`, oldest first:

```
e029552  Add the one-panel supplementary figure for the specificity claim
b17fcdf  Take the PINK1/PRKN variant off every figure, and say what that costs
53c9947  Write the human-arm synthesis, aimed at the mouse reconciliation
bd42e4c  Put the paper's additive-ratio number in a script, and confirm the E10 re-run
283333f  Put OXPHOS and MYC on the same row, one row per gene
36e5ee5  Add the two-bar summary: mean +/- SD of the per-gene difference, by half
789dddb  Add the difference on its own, gene by gene, a panel per cohort
448ecf8  Put both cohorts in one panel on fig5, and colour the gene names
6fcf4cd  Put figure 4's test onto the two priming heatmaps
7cf45c0  Add the ** mark: starred cells whose two genes are both OXPHOS-led
6cca4c6  Collapse the two axes into one number per ratio, both possible ways
5db8f0c  Add the ratio as a third panel - the better sentence, the worse picture
6695d07  Rewrite the reconciliation document around priming regulation
```

**Nothing has been pushed.** `results/` and `outputs/` are gitignored and
regenerable; `data/from_validation/` is regenerable by re-copying from
`myc_human_validation` at `d3ac60e`.
