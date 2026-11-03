# Source of the E49 review pack

This folder generates `docs/2026-11-03_e49_review_pack.qmd`, which renders to
the PDF beside it. It is not a numbered script: it computes nothing from data.

**2026-11-03, after STOP 3: `pack_text.py` now holds the author's STOP 3
rulings' text** (`docs/2026-11-03_stop3_rulings_consolidated.md`, `cb9b7ce`).
- It is the single source for E49's result note, Part 4, and for the E50 review
  pack (`docs/2026-11-03_e50_review_pack_src/`). Run it directly to check the
  text against the rulings file and print the metrics.
- **`make_pack.py` here stays pinned to E49's text**, and stops with a message if
  run against the new `pack_text.py`.
- **To rebuild E49's pack as committed,** run this folder from `c117f8f`:

  ```
  git worktree add ../e49_pack_c117f8f c117f8f
  cd ../e49_pack_c117f8f
  python3 -I -B docs/2026-11-03_e49_review_pack_src/make_pack.py docs/2026-11-03_e49_review_pack_src docs/2026-11-03_e49_review_pack.qmd
  ```

  From `c117f8f`'s files this reproduces the committed `.qmd` byte for byte
  (checked 2026-11-03). Remove the worktree afterwards with `git worktree remove`.

What follows describes the folder as it was at `c117f8f`.

**The files:**
- **`pack_text.py`** holds the result note's Part 4 text:
  - each sentence, with its id (S1, S2b, B3 ...);
  - `LEGEND_ADDITIONS`, the proposed additions to the 4g legend (result note
    4.3), which are not in the script's draft;
  - the STYLE.md word counter, `metrics()`.

  Run it directly to print the metrics.
- **`make_pack.py`** writes the `.qmd`. It reads:
  - `pack_text.py`;
  - `legends.json`;
  - E48's `../2026-11-03_e48_review_pack_src/legends.json`, for ED-a to ED-c.
- **`legends.json`** holds E49's legend drafts, verbatim from
  `results/e49_signature_and_fig4g.rds` (`$legend_draft`). It was checked
  identical, entry by entry, on 2026-11-03.

**The check.** Regenerating the `.qmd` from these files reproduces the committed
one byte for byte (2026-11-03).

**After the author's decisions:**
1. Edit the sentences, or `LEGEND_ADDITIONS`, in `pack_text.py`.
2. Regenerate and render from the repo root:

   ```
   python3 -I -B docs/2026-11-03_e49_review_pack_src/pack_text.py
   python3 -I -B docs/2026-11-03_e49_review_pack_src/make_pack.py docs/2026-11-03_e49_review_pack_src docs/2026-11-03_e49_review_pack.qmd
   /Applications/Positron.app/Contents/Resources/app/quarto/bin/quarto render docs/2026-11-03_e49_review_pack.qmd --to typst
   ```

   - The first line prints the metrics; the next two write the `.qmd` and render it.
   - `-B` keeps Python from writing a `__pycache__` into `docs/`.
   - The render cache is gitignored.
3. Carry the same edits into Part 4, section 4.1 and section 4.3 of
   `docs/2026-11-03_e49_result.md`. That text was printed from the same `PARAS`
   list and `LEGEND_ADDITIONS`.

**The word rule.** A word is any whitespace token that contains a letter or
digit. Parenthetical figure citations and reference numbers are not counted.
