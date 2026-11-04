# Source of the E50 review pack

This folder generates `docs/2026-11-03_e50_review_pack.qmd`, which renders to
the PDF beside it. It is not a numbered script: it computes nothing from data.

**The text is the author's STOP 3 rulings'** (`docs/2026-11-03_stop3_rulings_consolidated.md`,
`cb9b7ce`), as their addendum amends it (`docs/2026-11-03_stop3_rulings_addendum.md`,
`3f9ffbe`). It is section 3's sentences, with E47 declaration 4a.7, item 5, as the
Methods definition, and section 4's edits to the 4g legend. It was first E49's
(2026-11-03, `2a939bf`); the rulings replaced it.

**The files:**
- **`make_pack.py`** writes the `.qmd`. It reads:
  - **the text**, from its single source, E49's
    `../2026-11-03_e49_review_pack_src/pack_text.py`. That file holds:
    - the sentences, each with its id: A1 and A2, ML, R1 to R16, C1 to C3, Z1 and
      Z2, the Methods definition D1 to D12, and N1 to N3. No sentence is
      Claude Code's;
    - the rulings' edits to the 4g legend, items 1 to 5 applied and item 6
      waiting on the binding classes;
    - `check_against_rulings()`, which compares every sentence and edit with the
      rulings, the addendum, E47's declaration (4a.7, item 5) and the binding ruling. `make_pack.py`
      runs it first and stops on any difference;
    - `outside_range()`, which lists the sentences outside STYLE.md's 12 to 38
      words. The pack reports them; it does not stop on them;
    - the STYLE.md word counter.
  - **this folder's `legends.json`**;
  - **E49's and E48's `legends.json`**, for ED-i to ED-k and ED-a and ED-b;
  - **`data/annotations/bh3_binding_code.csv`**, only to check that
    `BINDING_CLASSES_ENTERED` in `pack_text.py` matches the table. It stops if
    they disagree.
- **`legends.json`** holds E50's legend drafts for 4g, ED-e, ED-f, ED-g and
  ED-h, verbatim from `results/e50_fig4g_in_fig1g_style.rds` (`$legend_draft`).
  It was checked identical, entry by entry, on 2026-11-03. The 4g legend shown
  in the pack is this draft with the rulings' edits applied.

**What the pack shows.** The nine Extended Data panels in the paper (the rulings'
1.4). Where each is cited is read from the text and the legend, and checked
against 1.4. ED-c and S-1 stay in the repo, out of the paper.

**To regenerate**, from the repo root:

```
python3 -I -B docs/2026-11-03_e49_review_pack_src/pack_text.py
python3 -I -B docs/2026-11-03_e50_review_pack_src/make_pack.py docs/2026-11-03_e50_review_pack_src docs/2026-11-03_e50_review_pack.qmd
/Applications/Positron.app/Contents/Resources/app/quarto/bin/quarto render docs/2026-11-03_e50_review_pack.qmd --to typst
```

- The first line checks the text against the rulings, the addendum, E47's
  declaration and the binding ruling, and prints the metrics.
- `-B` keeps Python from writing a `__pycache__` into `docs/`. The render cache
  is gitignored.

**The binding classes are entered** (`docs/2026-11-04_binding_classes_ruling.md`,
`66c205a`). `BINDING_CLASSES_ENTERED` is True, so legend edit 6 applies, and the
binding ruling adds edit 7, which cites Ku et al. 2010 and Wang et al. 2023 beside
Chen et al. 2005. `make_pack.py` stops if the table's `BID`, `BIK` and `BMF` rows
differ from the ruling, or if the flag and the table disagree.
