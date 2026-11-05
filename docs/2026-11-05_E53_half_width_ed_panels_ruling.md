---
date: 2026-11-05
status: >
  The author's ruling that Extended Data Fig. 6a and 6b, E48's ED-a and ED-b,
  are redrawn at half width with edited labels, and the four choices that
  complete it. Recorded verbatim and committed alone, before E53 is written.
  It is E53's specification.
relates-to:
  - docs/2026-11-03_e48_result.md (E48, run and verified; ED-a and ED-b are its panels)
  - docs/2026-11-03_e48_prerun_digests.md (1ccf805; the digests E53 checks)
  - docs/2026-11-03_stop3_rulings_consolidated.md (cb9b7ce; 1.4, the Extended Data list, which keeps ED-a and ED-b)
  - myc_mouse 2ebf646, figures/natmetab/FIGURE_RULES.md, "Page (Nature)" (Extended Data at most 180 mm wide; read only)
---

# Extended Data Fig. 6a and 6b at half width (the author, 2026-11-05)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message.
- **Section 3** gives the author's answers to four questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown, with the preview of the chosen option where there was one.

**The author's numbering:** ED-a (`E48_score_correlations`) is Extended Data Fig. 6a, and ED-b (`E48_bcl2l1_against_bbc3`) is Fig. 6b. Files stay named by content.

**What it changes.** A new display script, `E53`, redraws both panels from E48's verified saved object. Each panel is at most 88 mm wide, so that two fit the 180-mm Extended Data page.
- **6a:**
  - the TCGA-BRCA and SCAN-B triangles stay side by side, with 4.6-mm cells;
  - the three neoadjuvant values move to a row of tiles under the colour key (question 1).
- **6b:**
  - the direction arrows and their words ("lower expression", "higher expression", "BBC3 steeper", "BCL2L1 steeper") are removed;
  - the row labels become "TCGA: OXPHOS level", "TCGA: OXPHOS priority", "SCAN-B: OXPHOS level" and "SCAN-B: OXPHOS priority" (question 3);
  - in the key, BCL2L1's two entries become "BCL2L1" and "proliferation adjusted", and BBC3's and the difference's follow the same pattern, as "BBC3" and "Difference" (question 4);
  - the x axes start at −0.1: the coefficients run to 0.6, and the difference to 0.3 (question 2).

**Nothing in sections 2 and 3 is edited here.**

## 2. The author's message, verbatim

> now the first two supplementaries: ED Fig 6a: `E48_score_correlations` 6b: `E48_bcl2l1_against_bbc3`. Both should be reduced in size. they are currently full width, should be reduced to half width if possible. For 6b, the arrows with explanantion (lower.., higer.., steeper) from the panel should be removed the set labesl should read: `TCGA: OXPHOS priority, TCGA: OZPHOS levels, etc. The lgegend key: BCL2L1, without proliferation term should be simply BCL2L1, and the BCL2L1 with proliferation term, should read proliferation adjusted.

## 3. The author's answers, verbatim

**Question 1:**

> ED Fig. 6a (E48_score_correlations) cannot keep its three blocks side by side at half width: at the 5-pt floor they need about 100 mm, and half width is 88 mm. Which layout should the new script (E53) draw?

**The author's choice:**

> Side by side, smaller (Recommended): TCGA-BRCA and SCAN-B triangles side by side as now, with 4.6-mm cells. The three neoadjuvant values move to a row of tiles under the colour key. About 85 x 65 mm.

**The preview shown with it:**

```
         TCGA-BRCA (n = 938)  SCAN-B (n = 3,143)
OX prior.  #                  #
MYC_UP     # #                # #
Felsher    # # #              # # #
Hallmark   # # # #            # # # #
CollecTRI  # # # # #          # # # # #
Prolif.    # # # # # #        # # # # # #
          (labels at 40 deg)  (labels at 40 deg)

 Spearman rho, absolute  [ ] <0.60  [:] 0.60-0.80  [#] >=0.80
 Neoadjuvant cohorts, reference MYC score vs proliferation
   GSE25066 [:] 0.78  GSE194040 [:] 0.79  GSE164458 [#] 0.82

                                        ~85 x 65 mm
```

**Question 2:**

> With the direction arrows removed from ED Fig. 6b (E48_bcl2l1_against_bbc3), its symmetric x axes (-0.6 to +0.6, and -0.3 to +0.3) leave the left half of each panel empty, because every estimate is positive. At half width, which x axes should it use?

**The author's choice:**

> Start at -0.1 (Recommended): Coefficients from -0.1 to 0.6, the difference from -0.1 to 0.3. The zero line stays, with room on its left, so each interval is seen clear of it, and the data spread over about twice the width.

**The preview shown with it:**

```
TCGA: OXPHOS level      |      o-+-o
                        |    <>-+-<>
                     ---+----+----+----+
                     -0.1 0  0.2  0.4  0.6
```

**Question 3:**

> For ED Fig. 6b's row labels you wrote 'TCGA: OXPHOS priority, TCGA: OXPHOS levels'. Fig. 6a and Fig. 4g say 'OXPHOS level'. Which should the four labels use? The order stays as now: level, then priority, TCGA then SCAN-B.

**The author's choice:**

> OXPHOS level (Recommended): TCGA: OXPHOS level / TCGA: OXPHOS priority / SCAN-B: OXPHOS level / SCAN-B: OXPHOS priority. Singular, as in 6a and 4g.

**Question 4:**

> In ED Fig. 6b's key, BCL2L1's two entries become 'BCL2L1' and 'proliferation adjusted'. The script would give BBC3 and the difference (black symbols) the same pattern. What should the difference's entry say?

**The author's choice:**

> Difference (Recommended): E48's word for the black symbols. The right panel's axis title already says 'BCL2L1 minus BBC3'.

**The preview shown with it:**

```
(o) BCL2L1       (*) proliferation adjusted
<>  BBC3         <*> proliferation adjusted
(o) Difference   (*) proliferation adjusted

open = without the proliferation term
filled = with it
```
