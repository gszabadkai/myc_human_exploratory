---
date: 2026-11-06
status: >
  The author's ruling on the last two main panels of Fig. 4, the outcome
  forests, and the four choices that specify them. Recorded verbatim and
  committed alone, before anything is written. It is the specification of
  E57, a display script that redraws E45's two forests as two separate
  panels.
relates-to:
  - docs/2026-10-22_fig4_human_panels_proposal.md (E45's specification, its decisions in section 7 and its recorded digests)
  - scripts/E45_fig4_human_panels.R and docs/figures/2026-10-22_E45_outcome_forests.{pdf,png} (the panels E57 redraws)
  - docs/2026-11-05_e54_result.md (Fig. 4x, whose axis length, type, line weights, circles and dark slate E57 follows)
  - docs/2026-11-05_E51_main_panel_ruling.md (question 3: the dark and mid slates, the ends of Fig. 1's unsigned palette)
---

# Fig. 4's outcome forests as two panels (the author, 2026-11-06)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message.
- **Section 3** gives the author's answers to four questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown.

**The "last two main figure panels"** are E45's forests: pathological complete response (pCR) in three neoadjuvant cohorts with their pooled row, and breast-cancer-specific survival in METABRIC. The Google Doc lettered them 4i and 4j.

**What it changes, by a new display script, `E57`:**
- **Two separate panels,** one per forest. E45 drew them side by side in one figure.
- **The colours** (question 1):
  - **OXPHOS** takes Fig. 4x's dark slate (#243642): filled when adjusted for proliferation, open when not;
  - **proliferation** takes mid slate (#8FA3B0), the other end of Fig. 1's unsigned palette.

  They replace E45's Okabe-Ito green and orange, which Figs 1 to 3 use for the PGC-1α construct and the MYC-positive 6-week mice.
- **All symbols are circles.** Proliferation's triangle becomes a circle, told apart by its colour.
- **No bold characters.**
- **The keys:**
  - "OXPHOS" for the model without the proliferation term;
  - "OXPHOS, adjusted for proliferation" for the model with it;
  - "Proliferation, adjusted for OXPHOS" for proliferation in that model (question 4).
- **The size** (question 2): a 26-mm axis, with Fig. 4x's type sizes, line weights and 1.05-mm circles.
- **The direction cues stay** (question 3).
- **Otherwise the content is E45's:**
  - the rows and their labels;
  - the ratios and intervals;
  - the contrast column;
  - the axes' scales.

**Nothing in sections 2 and 3 is edited here.**

## 2. The author's message, verbatim

> OK, for the last two main figure panels we need to use `2026-10-22_E45_outcome_forests.pdf`. I need two separate panels. Find the best matching colour code, no bold characters, otherwise the content is OK. THe legend keys again OXPHOS; OXPHOS, adjusted for proliferation. THe symbols whould be all circle, but proliferation should be different colour.

## 3. The author's answers, verbatim

**Question 1:**

> Colours for the two outcome forests. E45's green (#009E73) and orange (#D55E00) are the paper's Okabe-Ito hues for the PGC-1α construct and the MYC-positive 6-week mice in Figs 1–3. In every option OXPHOS takes Fig. 4x's dark slate (#243642), filled when adjusted for proliferation and open when not. Which colour should proliferation take? (Fig. 1g's brown is ruled out: it is too close to dark slate under colour-vision deficiency, OKLab dE 6–9.)

**The author's choice:**

> Mid slate (Recommended): #8FA3B0, the other end of Fig. 1's unsigned palette, the second colour of the PUMA/BCL-xL panel. Separated from dark slate by lightness under every colour-vision deficiency (OKLab dE 38), but only 2.6:1 against white, so its bars read light.

**Question 2:**

> Size of the two forest panels. E52's panel and Fig. 4x are drawn at the Prism graphs' size: a 26-mm axis, 5-pt text and 6-pt titles. E45's forests have axes of about 47 and 36 mm and 6.5-pt text, and are 91 and 84 mm wide. Which size should the two separate panels take?

**The author's choice:**

> 26-mm axis (Recommended): Fig. 4x's axis length, type sizes, line weights and 1.05-mm circles. Each panel about 65–70 mm wide, so the two fit side by side with room to spare.

**Question 3:**

> Direction cues on the outcome forests. E45 draws arrows under each axis: 'lower pCR' / 'higher pCR' and 'longer survival' / 'shorter survival'. A higher odds ratio is benefit and a higher hazard ratio is harm, which is why they were added. You had similar arrows removed from Extended Data Fig. 6b. Keep them here?

**The author's choice:**

> Keep them (Recommended): Content as E45 has it. At a 26-mm axis they are set at 5 pt, and 'longer survival' runs a little past the axis's left end.

**Question 4:**

> Proliferation's key on the outcome forests. Your OXPHOS keys are 'OXPHOS' (no proliferation term) and 'OXPHOS, adjusted for proliferation'. Proliferation's odds and hazard ratios come from the adjusted model, with OXPHOS in it. What should proliferation's key say?

**The author's choice:**

> Adjusted for OXPHOS (Recommended): 'Proliferation, adjusted for OXPHOS', parallel with the OXPHOS keys. A bare 'Proliferation' would read as unadjusted, as 'OXPHOS' does.
