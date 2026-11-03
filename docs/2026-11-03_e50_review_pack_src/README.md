# Source of the E50 review pack

This folder generates `docs/2026-11-03_e50_review_pack.qmd`, which renders to
the PDF beside it. It is not a numbered script: it computes nothing from data.

**The files:**
- **`make_pack.py`** writes the `.qmd`. It reads:
  - **the text**, from its single source, E49's
    `../2026-11-03_e49_review_pack_src/pack_text.py`. That file holds:
    - the result note's Part 4, sentence by sentence;
    - the proposed legend additions;
    - the STYLE.md word counter.
  - **this folder's `legends.json`**;
  - **E49's and E48's `legends.json`**, for ED-i to ED-k and ED-a to ED-c.
- **`legends.json`** holds E50's legend drafts for 4g, ED-e, ED-f, ED-g and
  ED-h, verbatim from `results/e50_fig4g_in_fig1g_style.rds` (`$legend_draft`).
  It was checked identical, entry by entry, on 2026-11-03.

**The text is E49's, not yet amended.** The author's STOP 3 rulings on it are
not in the repo. When they are supplied:
1. Edit the sentences, or `LEGEND_ADDITIONS`, in E49's `pack_text.py`, and
   re-measure: `python3 -I -B docs/2026-11-03_e49_review_pack_src/pack_text.py`.
2. Regenerate and render from the repo root:

   ```
   python3 -I -B docs/2026-11-03_e50_review_pack_src/make_pack.py docs/2026-11-03_e50_review_pack_src docs/2026-11-03_e50_review_pack.qmd
   /Applications/Positron.app/Contents/Resources/app/quarto/bin/quarto render docs/2026-11-03_e50_review_pack.qmd --to typst
   ```

   `-B` keeps Python from writing a `__pycache__` into `docs/`. The render
   cache is gitignored.
3. Carry the same edits into Part 4, sections 4.1 and 4.3 of
   `docs/2026-11-03_e49_result.md`.
