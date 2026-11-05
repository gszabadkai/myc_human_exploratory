---
date: 2026-11-05
status: >
  The author's ruling on Extended Data Fig. 6c and 6d, and the four choices
  that specify them. Recorded verbatim and committed alone, before anything
  is declared or written. It is the specification of E55 (Fig. 6c) and E56
  (Fig. 6d, which first fits the three effectors under its own declaration).
relates-to:
  - docs/2026-11-03_e49_result.md (E49; its raw correlations, the left side of ED-h)
  - docs/2026-11-03_e50_result.md (E50; Fig. 4g, the heatmap whose genes, groups, size and colours these panels follow)
  - docs/2026-10-22_fig4_human_panels_proposal.md and the E46 record (E46, the quadrant panel 6d redraws)
  - docs/2026-10-01_e34_result.md (E34, the quadrant model 6d's new fits repeat)
---

# Extended Data Fig. 6c and 6d (the author, 2026-11-05)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message.
- **Section 3** gives the author's answers to four questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown.

**"E52" in the message is read as Fig. 4g (E50's heatmap),** as question 1 stated. E52 is the PUMA/BCL-xL scatter and has no gene list. The author did not object.

**What it changes:**
- **Fig. 6c**, by a new display script, `E55`:
  - it shows the left side of E49's raw-beside-adjusted heatmap: unadjusted Spearman correlations with all seven measures, in both cohorts (question 1);
  - it uses Fig. 4g's 15 genes, order and four groups, with one signature row (S6) at the bottom;
  - it takes 4g's heatmap size and colour code, with colour only and no printed values.
- **Fig. 6d**, by `E56`, which first fits and then draws:
  - **the effectors BAX, BAK1 and BOK are fitted** with E34's exact quadrant model on the same tumours. This is a new quantity, post hoc and labelled so, and its declaration is committed before the fit (question 2);
  - **the panel then redraws E46's quadrant panel** with 4g's 15 genes and four groups, in 4g's colour code. Both OXPHOS splits (GSVA and mitoPPS) sit side by side at full width (question 3);
  - **the colour key reads "z-score"** in place of "Adjusted mean (within-cohort SD)";
  - **the dot plot's key reads "Signature" and "proliferation adjusted"** (question 4).

**Nothing in sections 2 and 3 is edited here.**

## 2. The author's message, verbatim

> ok. next are two extended data panels. ED Fig. 6c.: for this we need to split the `2026-11-03_E49_raw_beside_adjusted_unmarked.pdf` panel. WE use only the left side, spearman for each MYC activity measure. THe genes and order should be the same as in E52, show only one signature at the bottom. THe size of the heatmap the same, the colour code the same. ED Fig. 6d.: we use the `2026-10-22_E46_quadrant_configuration.pdf`. Also here the same genes, same grouping as in E52, same colour code. instead of the legend key `adjusted mean (within cohort SD)` just write `z-score`. For the dot-plot below legend key, just say signature and adjusted to proliferation.

## 3. The author's answers, verbatim

**Question 1:**

> For ED Fig. 6c, taken from the left side of E49's raw-beside-adjusted heatmap (Spearman correlations, unadjusted), which columns should it keep? I take 'E52' in your message to mean Fig. 4g (E50's heatmap): its 15 genes, order, four groups, one signature row at the bottom, heatmap size and colours.

**The author's choice:**

> All seven columns: Also OXPHOS level, OXPHOS priority and proliferation: 14 columns in 4g's heatmap size, cells too narrow for printed values (colour only).

**Question 2:**

> ED Fig. 6d redraws E46's quadrant panel with Fig. 4g's genes and groups. But E46 (and E34 before it) computed quadrant means for only 12 genes: the effectors BAX, BAK1 and BOK were added to the paper later and were never fitted by quadrant. What should 6d show?

**The author's choice:**

> Fit the effectors (Recommended): Declare first, then fit BAX, BAK1 and BOK with E34's exact quadrant model on the same tumours, so 6d carries 4g's 15 genes in its four groups. A new, post hoc quantity, labelled as such.

**Question 3:**

> E46 shows the quadrants twice, with OXPHOS split by GSVA (level) and by mitoPPS (priority), side by side at full width (about 180 mm), each with its signature dot plot below. What should ED Fig. 6d show?

**The author's choice:**

> Both splits, full width (Recommended): As E46: level and priority side by side, both cohorts, about 180 mm wide, with the signature dot plots below.

**Question 4:**

> For the key under 6d's signature dot plot (now 'Signature, m1' filled and 'Signature, m1 + proliferation (m3)' open), you asked for 'signature' and 'adjusted to proliferation'. Which wording?

**The author's choice:**

> Signature / proliferation adjusted (Recommended): Matches 6b's key, where you asked for 'proliferation adjusted'.
