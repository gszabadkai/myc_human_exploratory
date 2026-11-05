---
date: 2026-11-05
status: >
  The author's ruling that E51's panel is redrawn at the size of Fig. 4's
  Prism graphs, with a 26-mm axis and axis titles without "adjusted", and the
  two choices that complete it. Recorded verbatim and committed alone, before
  E52 is written. It is E52's specification.
relates-to:
  - docs/2026-11-05_E51_main_panel_ruling.md (86f38aa; E51's specification, which stands except where this changes it)
  - docs/2026-11-05_e51_result.md (E51, run and verified)
  - myc_mouse 2ebf646, figures/natmetab/FIGURE_RULES.md, "The Prism format" (read only)
---

# E51's panel at the Prism graphs' size (the author, 2026-11-05)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message.
- **Section 3** gives the author's answers to two questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown, including the preview of the chosen layout.

**What it changes from E51** (whose specification otherwise stands):
- **The data area becomes 26 x 26 mm,** the axis length of Fig. 4's Prism graphs. The page is the panel's natural size.
- **"adjusted" leaves both axis titles.** The legend carries the adjustment.
- **The two *P* labels, with the cohort names,** move to a small key above the data area. They are still part of the panel.
- **The axes and ticks stay at 0.53 pt.**

**It supersedes the earlier note** that the author would remove "adjusted" in Illustrator. E52 draws the titles without it.

**Nothing in sections 2 and 3 is edited here.**

## 2. The author's message, verbatim

> actually, the pnael is too large, it should fit in the overall prism panel sizes, so the axis length should be 26 mm, and in this case you can also edit the axis title, removing `adjusted`

## 3. The author's answers, verbatim

**Question 1:**

> At a 26-mm axis the one-line P labels ("TCGA-BRCA  P = 0.017", 18.8 mm each) would sit on the top SCAN-B points. Where should the P values go?

**The author's choice:**

> Above the graph: The two labels, with cohort names, sit as a small key above the data area, still part of the panel. The y axis stays at ±0.3 and nothing overlaps; the panel is about 6 mm taller.

**The preview shown with it:**

```
 ━ TCGA-BRCA  P = 0.017
 ━ SCAN-B  P = 7.8×10⁻⁹
  +0.3 ┤
       │  ○ ○
       │ •  ○ •
     0 ┼────•○──────
       │       ○ • ○
  −0.3 ┤          • ○
       └──┬──┬──┬──┬
       −2  −1  0  +1  +2
```

**Question 2:**

> Fig. 1's figure rules match an R panel's axis and tick lines to its Prism neighbours: they measure 0.9 pt in Fig. 1 and 0.72 pt in Fig. 3. Fig. 4's Prism graphs have not been measured. What weight should this panel's axes and ticks take?

**The author's choice:**

> Keep 0.53 pt: The R panels' default, as E51 drew it. The rules note it reads visibly thinner beside Prism graphs.
