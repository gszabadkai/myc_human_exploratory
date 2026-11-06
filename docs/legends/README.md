# Figure legends for the Nat Metab Letter

**What this folder is.** The source of the figure legends that are not yet
written, and the two tools that read the author's own files to get their
numbers. Fig. 1 and Extended Data Fig. 1 to 3 are already written, by the
author, in the Doc. This folder covers **Fig. 2, Fig. 3, Fig. 4, Extended Data
Fig. 4 and Extended Data Fig. 6**. Extended Data Fig. 5 is skipped: the author
says it is not finished (2026-11-06).

**Nothing here is written to the Doc.** The author pastes.

## The rules these legends follow

**The author's rules for legends** are in the mouse repo, read-only at
`2ebf646:figures/natmetab/LEGENDS.md` (the author, 2026-10-16):

> describe what is on each panel and give what is needed to read it; no extra
> detail and no interpretation, except the figure title and, where useful, a
> short sentence naming what a panel shows. Nature Metabolism also asks for the
> centre value, error bars, *n*, the test and the *P* values, and for keys in
> words, not symbols.

**The author's naming decisions** (2026-11-06):
- the labels on the figures are followed;
- mouse proteins in title case (Myc, Bcl-xL, Pgc-1α, Ndi1); mouse RNA in title
  case italics (*Bbc3*, *Ppargc1a*);
- human proteins in capitals; human RNA and DNA in capital italics;
- **no titles for Extended Data figures**;
- blot and image details are placeholders, `[...]`, for the author to fill.

**From `2ebf646:figures/natmetab/FIGURE_RULES.md`**, which the panels were drawn
to: exact *P* to two significant figures with an italic *P*; a test named
wherever a glyph or a bracket marks significance; black text only.

**The Doc's Methods** says the exact test and *n* are reported in the legend,
so every data panel here carries both.

## Where the numbers come from

| what | source |
|---|---|
| panels, their letters and what is drawn | the assembled PDFs in the author's Drive folder `1qLbI6nQlLczmUBrmWyMvMmNPJoG6Fano` |
| *n*, groups, tests, corrections, adjusted *P* | `NM_all_figs.prism` in that folder (modified 2026-11-06 14:58) |
| what the error bars are | measured from the PDFs' own vector geometry against the Prism values |
| the human panels (Fig. 4f to 4j, ED Fig. 6) | this repo's verified result notes and pre-run notes for E50 and E52 to E57 |
| method essentials | the Doc's Methods, CURRENT VERSION tab |
| the superseded legends | `11.8.26_MK_Myc_Paper_figures_tab.pdf`, the old Figures 3 to 5 and Supplementary Figures 3 to 5 |

**The two tools** read those files and compute nothing of their own:

```
python3 tools/prism_inventory.py NM_all_figs.prism [title prefix ...]
python3 tools/pdf_error_bars.py NM_Fig4.pdf <xmin> <xmax> <ymin> <ymax>
```

`prism_inventory.py` prints, per data sheet, the groups, the *n*, the mean,
s.d. and s.e.m., and then each analysis on that sheet with its test, design,
post-test and every adjusted *P* or *q* with the two groups it compares.
`pdf_error_bars.py` prints every drawn two-point segment in a window of a page,
in points, so an error bar can be read against its own axis.

**Neither file is in this repo.** Both are the author's, in Drive. The tools
take a path.

## Read in this order

1. **`2026-11-06_panel_inventory_and_stats.md`** — every panel, what it draws,
   its *n*, its test, and **eleven flags that need the author's decision**.
   Several flags change what a legend may say, so they come first.
2. **`2026-11-06_legends_draft.md`** — the legends, with `[...]` where the
   author fills in and `[Fn]` where a flag is unresolved.
