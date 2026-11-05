---
date: 2026-11-05
status: >
  The author's ruling that ED-k becomes a main panel, redrawn with both cohorts
  overlaid and a P value per cohort, and the four choices that specify it.
  Recorded verbatim and committed alone, before E51 is written. It is E51's
  specification.
relates-to:
  - docs/2026-11-03_stop3_rulings_consolidated.md (cb9b7ce; 1.4, the Extended Data list, which this changes for ED-k)
  - docs/2026-11-03_e49_result.md (1.7, the relation ED-k draws)
  - docs/Fig. 4 narrative v6 human arm and clinical link.md (c3b3482)
  - myc_mouse 2ebf646, figures/natmetab/fig1i_coupling.R and FIGURE_RULES.md (Fig. 1i, the panel this one is set beside; read only)
---

# ED-k becomes a main panel (the author, 2026-11-05)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message of 2026-11-05.
- **Section 3** gives the author's answers to four questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown.

**What it changes:**
- **ED-k leaves Extended Data and becomes a main panel.** It is redrawn by a new display script, `E51`, with TCGA-BRCA and SCAN-B overlaid in one panel and no title.
- **Each cohort's *P* value is printed inside the panel.** It is the two-sided t-test of the drawn slope, a quantity not computed before, so E51 declares it.
- **The model stays ED-k's**, on the reference MYC activity score, and the colours are a dark and a mid slate.
- **ED-b stays in Extended Data.**
- **Two text changes follow**, for a later ruling on the text:
  - sentence R6 cites the new panel in place of Extended Data Fig. Xk;
  - the main panel's letter is to be set.

**Nothing in sections 2 and 3 is edited here.**

## 2. The author's message, verbatim

> and ED-k, `2026-11-03_E49_puma_bclxl_against_priority.pdf` has the same message, so I would use that as a main panel. However, I would need a version which overlays the two sources TCGA and SCAN-B in one panel, with different colours. That figure reports the slopes with CI and partial correlation. would it be more simple just to report a p-value? Place the p values inside the panel, no need for title

## 3. The author's answers, verbatim

**Question 1:**

> For the new main panel (ED-k redrawn: TCGA-BRCA and SCAN-B overlaid, PUMA relative to BCL-xL against OXPHOS priority), which P value should it print for each cohort?

**The author's choice:**

> Model t-test (Recommended): Two-sided t-test of the drawn slope in its own model; the same test as the partial correlation. Exact, two significant figures, italic P, as Fig. 1i's rules ask. Computed from E49's saved residuals, no refit.

**Question 2:**

> ED-k's lines adjust for MYC activity on the reference (Felsher) score, whereas Fig. 4g uses MYC_UP.V1_UP. Which model should the main panel use?

**The author's choice:**

> Reference score (Recommended): ED-k's verified model, as declared: nothing is refitted, the slopes and intervals are the stored ones, and the legend names the MYC score.

**Question 3:**

> Which two colours should distinguish TCGA-BRCA and SCAN-B? The Okabe-Ito hues already mean mouse genotypes and cell-line constructs in Figs 1-3.

**The author's choice:**

> Dark and mid slate (Recommended): #243642 and #8FA3B0, the ends of Fig. 1's unsigned palette. Outside every mouse and construct hue, and told apart by lightness under any colour-vision deficiency.

**Question 4:**

> ED-b (BCL2L1 against BBC3) shows the same difference on OXPHOS level and priority; the new panel shows only priority. What happens to ED-b?

**The author's choice:**

> Keep ED-b (Recommended): Sentence R5 cites it for 'BCL2L1 rising more steeply with OXPHOS than BBC3', which rests on OXPHOS level, and the new panel does not show level.
