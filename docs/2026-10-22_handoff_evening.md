---
date: 2026-10-22
status: >
  HANDOFF at a clean stop. FIGURE WORK IN PROGRESS, now paused while the
  author rethinks the narrative.
  - Fig. 4i and 4j (the Doc's lettering) are BUILT by E45, run by the author
    and VERIFIED.
  - Fig. 4g is WRITTEN as E46, on both instruments, and NOT yet run by the
    author. Its digests were committed BEFORE that run, at e12e94f.
  - No phase is open and nothing is in flight.
written-from: >
  main at e12e94f, one commit after 9cef8ae. Saved objects were read with
  Rscript readRDS only. Every numbered script was dry-run with its outputs
  redirected to the session scratchpad; nothing was written in the repo
  except what was committed.
supersedes: >
  2026-10-22_handoff.md, as the entry point for Fig. 4 work. That file's
  sections 2-6 are now partly out of date (see section 7 below). Its
  section 4 is still a valid list of traps, with the correction in 7.1.
relates-to:
  - docs/2026-10-22_fig4_human_panels_proposal.md (sections 7 and 8 hold
    today's decisions, corrections and digests)
  - docs/2026-10-04_handoff.md section 2 (still the authority on what a
    panel may assert)
  - docs/2026-10-01_e34_result.md section 9 (what 4g may assert)
  - docs/handoff-2026-10-21.md (still the orientation entry point for the
    repo)
next-action: >
  1. The author rethinks the narrative in the Google Doc. Section 5 lists
     every sentence the built panels contradict, or that has no fitted
     quantity behind it.
  2. The author sources scripts/E46_fig4_quadrant_heatmap.R. Claude Code
     verifies the run against the committed digests (section 6).
  3. Decision 5 (the validation study in the paper) is taken. It now blocks
     a MAIN panel: the Doc's 4h.
---

# Handoff 2026-10-22, evening: two of three human panel sets built; the narrative is next

## 1. Where to start

**Read section 5 first.** It is the input to the narrative rewrite. Then:

1. **The Google Doc** "11.8.26_MK_Myc_Paper" (Nat Metab version). It is the
   narrative of record, by the author's decision today.
   - It is **not in the repo.** Read it through the Google Drive connector,
     file ID `1M-kEC1d_VJmbTG-SyYYfykZ3F1_aEiDtu0xtwFaPtJo`.
   - The human section is the "Main" text, Fig. 4 paragraphs.
   - The "Human" draft further down is older and is not the narrative.
2. **Proposal sections 7 and 8** record today's decisions, two corrections to
   this morning's notes, and the E45 and E46 digests.
3. **`2026-10-04_handoff.md` section 2** is unchanged in authority.

## 2. State

| item | state |
|---|---|
| `main` | `e12e94f`, ahead of origin; **not pushed** |
| next free script number | **`E47`** |
| **Fig. 4i** pCR forest and **Fig. 4j** METABRIC forest | **BUILT by `E45`.** Run by the author and **VERIFIED**: 9 of 9 digests, 27 of 27 assertions, PNGs byte-identical to the dry run. Its digests were committed *after* the run (proposal 7.6) |
| **Fig. 4g** configuration by quadrant, GSVA and mitoPPS | **WRITTEN as `E46`, NOT RUN by the author.** 19 assertions, 10 element digests, **committed before the run** at `e12e94f` |
| **Fig. 4h** in the Doc, the MYC / BCL2L1 / MCL1 co-amplification | **NOT BUILT, and it cannot be built from this repo.** It is validation-study G2, so it is blocked on decision 5 |
| ED Figs 7-10 | not specified beyond v5 section 4 |
| phase | none open. `E31` is still the only analysis off `main` |

**Tracked figures:**
- `docs/figures/2026-10-22_E45_outcome_forests.{png,pdf}` is committed.
- `docs/figures/2026-10-22_E46_quadrant_configuration.{png,pdf}` will appear
  when the author runs E46. **Commit it then.**

## 3. The narrative of record, and its lettering

The Doc supersedes v5 (`docs/2026-10-08_section4_narrative_v5_submission_draft.md`).

| Doc | v5 | content | built by |
|---|---|---|---|
| 4g | 4g | configuration by quadrant, **now two panels side by side**: OXPHOS by GSVA and by mitoPPS | `E46` (unrun) |
| 4h | ED 9 | MYC / BCL2L1 / MCL1 co-amplification | nothing; decision 5 |
| 4i | 4h | OXPHOS and proliferation on pCR | `E45` |
| 4j | 4i | the same on METABRIC breast-cancer-specific survival | `E45` |

**If 4g's two panels take two letters, every later letter shifts.** The scripts
name their files by content for exactly this reason. **Any reletter is the
author's call.**

## 4. The five decisions

| # | decision | status |
|---|---|---|
| 1 | 4g row blocks | **TAKEN**: by function, pro- then anti-apoptotic, composite members first |
| 2 | 4g instrument | **TAKEN, differently from the proposal**: both instruments in the main figure, side by side |
| 3 | axis for the forests | **TAKEN**: ratios on a log axis, 1 marked |
| 4 | 4j rows | **TAKEN** as proposed, following the Doc: the METABRIC primary, its two comparators, ER-negative. No GSE25066 DRFS row |
| 5 | the validation study's place in the paper | **OPEN. Since the Doc made co-amplification a main panel (4h), it now blocks a figure, not only text.** CLAUDE.md says the validation study "does not go in the paper". If it goes in, the audit's contradictions 1 and 2 travel with it (proposal section 5, item 5) |

**One design change the author did not ask for, and why: there is no bracket
on 4i** (proposal 7.2).
- A bracket between the pooled markers spans **2.02** on the axis.
- Its label would read **1.93**, the meta-analysis of the per-cohort
  contrasts.
- The contrast therefore runs in a column on every row instead.

## 5. FOR THE NARRATIVE REWRITE: what the Doc says against what was fitted

**The Doc was read, never edited.** Citations are to the Doc's Fig. 4
paragraphs (paragraphs 1-4) and its abstract. **Items A-F would be visible to a
reader comparing the text with the panels.** Items G-J are true now and need
no action beyond care.

### A. "The two OXPHOS-high groups showed comparable apoptotic signatures" (paragraph 1)

**It holds in one of the four cohort-by-instrument blocks.** E34's
MYC-at-OXPHOS-high contrast, `Q4 - Q2`, on the composite at `m1`:

| | GSVA | mitoPPS |
|---|---|---|
| TCGA | **-0.002** [-0.169, 0.164] | +0.185 [0.012, 0.358] |
| SCAN-B | +0.134 [0.037, 0.231] | +0.232 [0.134, 0.331] |

On mitoPPS the means are visibly apart: Q2 +0.21 against Q4 +0.40 in TCGA,
and +0.09 against +0.32 in SCAN-B. **Proliferation adjustment (`m3`) enlarges
every MYC contrast, by about +0.12.**

**E34 licenses "much smaller", never "none" or "comparable"** (E34 result 9.1
and 9.2). Its licensed sentence quotes GSVA means only. **With both
instruments on the page, any sentence quoting the means must quote both.**

### B. "Remained invariant after controlling for cell proliferation rate, tumour purity, and immune cell infiltration" (paragraph 1)

- **SCAN-B has no purity or leukocyte estimate** (trap 2). Its `m1` is
  subtype only.
- **Proliferation is not invariant for the MYC contrasts**, which grow by
  about 0.12 under `m3`. The OXPHOS contrasts move little, but they are
  near-guaranteed by construction (E34 1.2).

**Say which quantity is invariant, and in which cohort.**

### C. "The association with MYC activity was lost after adjustment for OXPHOS" (paragraph 1)

This is E11's reciprocal partials. **It was not re-checked today.** E34 shows
the MYC contrast is small but present in 6 of 8 pooled cells. **Check E11's
numbers before keeping "lost".**

### D. "BCL2L1 showed a steeper association with OXPHOS than BBC3" (abstract and paragraph 2)

**This is not recorded in any note.** v5's pre-submission list requires it to
be filed before the abstract or Results state it. It needs a dated note with
its numbers, both cohorts and both instruments, before it stays.

### E. The co-amplification panel, Fig. 4h (paragraph 2)

This is validation-study G2, so it is decision 5.

### F. "Additional adjustment for the anti-apoptotic composite score or BCL2L1 expression did not materially alter the OXPHOS-response association" (paragraph 3)

- **The "anti-apoptotic composite" is E40's `m3` term `BUFFER_c`.** That is
  the mean of `z_MCL1` and `z_BCL2L1`, added on top of MYC (`m2`). Pooled
  `m3` is -0.1283 [-0.252, -0.005], against `m1` -0.1616. **The sentence is
  backed for that rung.**
- **Declaration 5.1 forbids describing `m3` as about "buffering capacity".**
  It is "adjustment for two anti-apoptotic transcripts", whose limbs
  correlate at -0.02 to 0.12.
- **"Or BCL2L1 expression" has NO fitted quantity behind it.** E40 has no
  `BCL2L1`-only rung or sensitivity; `z_BCL2L1` exists only inside
  `BUFFER_c`. **Drop the clause, or declare and fit the model first.**
- The sentence also omits that `m2` and `m3` both carry MYC.

### G. "After adjustment for subtypes and MYC expression" (abstract)

**It should say MYC activity.** MYC mRNA has essentially no correlation with
OXPHOS (rho -0.032; trap 4), so "expression" makes the claim about the wrong
quantity.

### H. "Across 4,302 human breast tumours" (abstract)

This is the cohort total, 1,095 + 3,207. **The adjusted quadrant models use
938 + 3,143 = 4,081.** If the sentence is about the adjusted result, use 4,081
or say "complete cases".

### I. Paragraph 1's per-gene directions

These cite an ED figure, which is correct: ED Fig. 10, the partial
correlations.

**In 4g's heatmap, `BID`, `PMAIP1` and `BCL2L11` are MYC-ordered or flat in
SCAN-B** (proposal section 2). A reader holding the main panel against
paragraph 2's list will see that, so the list should cite the ED figure
explicitly.

### J. Consistent with what was fitted, no action needed

- "OXPHOS alone showed no consistent association with pCR". `m0` fails the
  rule in both directions.
- The pooled OR of 0.85, and "the difference between their associations was
  statistically supported" (E44).
- The METABRIC sentence. It names the exposure, gives OXPHOS no survival
  claim, and rests on the direct comparison.
- 1,935 neoadjuvant patients, and METABRIC n 1,981.
- "Positively but moderately correlated": 0.358 to 0.493, all four cohorts
  recorded.

**Not checked against the Doc today:**
- paragraph 1's "restricted to apoptosis among four control mitochondrial
  pathways";
- the mouse paragraphs;
- the Discussion.

## 6. When the author runs E46: verification

Source `scripts/E46_fig4_quadrant_heatmap.R` in Positron. Then compare every
element except `built` against proposal 8.3:

```r
x   <- readRDS("results/e46_fig4_quadrant_panels.rds")
rec <- c(heatmap_cells = "98bad28c4bcf5fdd9362969cdf6fa32d",
         strip         = "a65aacc6636ededef218ee137062dde4",
         contrasts_e34 = "950c6e26c84d262bab136e538b88c82a",
         agreement_e34 = "7609a32fc8b58d09bb3d66ba61fc5616",
         sign_check    = "d0cc6384f4301602f1a4e8655d422fd1",
         analysis_n    = "ed776bc0b02f02ceb57ed48dcb1ae647",
         assertions    = "419a73c7e4fb2c4e899d789a3eaa2561",
         legend_draft  = "09e563930a25476928525b78e64ccf70",
         figure_files  = "24f2b96d6e876b1cfd5bbcfa200eec49",
         spec          = "4f47e6400aff77ff37472bd2b9ac6ed1")
got <- vapply(names(rec), function(k) digest::digest(x[[k]]), "")
table(got == rec)        # expect TRUE 10
```

**Then:**
1. Commit the tracked figure.
2. Record the verification in proposal section 8 or a result note.
3. Mark 4g BUILT in CLAUDE.md.

**If the narrative rewrite changes 4g's content, change the script before
running it**, then re-dry-run and re-commit its digests. Changes of this kind
include a different rung, different rows, or dropping an instrument.

## 7. Corrections made today to this morning's notes

1. **METABRIC's `rho(OX, PROLIF)` = 0.381 WAS recorded**, at
   `2026-10-03_e41_handoff.md:106` and in E41's `$exposure_corr`.
   `2026-10-22_handoff.md` 4.4 and proposal 4.2 item 6 said it was not.
2. **The METABRIC `m1` model includes treatment group**, as well as PAM50
   and `PROLIF`. Proposal 3.3's legend omitted it.
3. **Decision 2 went the other way from the proposal**: both instruments are
   main, not GSVA main and mitoPPS in ED.
4. `2026-10-22_handoff.md` section 2 calls Fig. 4g-4i "not built", and its
   decisions are now taken. **This file supersedes it for Fig. 4 work.**

## 8. Still open

- **Decision 5.**
- **ED Figs 7-10**, unspecified. mitoPPS is no longer ED Fig. 7's role for
  4g, so ED Fig. 7 needs re-scoping.
- **Push to origin**, when the author wants it.
- Carried from `handoff-2026-10-21.md` 8.1, unchanged:
  - the four unexecuted falsifier lists;
  - LbNOX;
  - the companion-paper figures;
  - the legend pass on the mouse paragraphs.

  "BCL2L1 exceeding BBC3" is now item D above.

## 9. Commit trail

| commit | what |
|---|---|
| `9cef8ae` | `E45`: Fig. 4i and 4j, run and verified; proposal section 7; CLAUDE.md |
| `e12e94f` | `E46`: Fig. 4g on both instruments, written and not run; proposal section 8 with the digests, committed before the run |
| this one | this handoff, and CLAUDE.md pointed at it |
