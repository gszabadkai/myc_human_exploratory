---
date: 2026-11-06
status: >
  HANDOFF for the figure legends. The legends for Fig. 2, Fig. 3, Fig. 4,
  Extended Data Fig. 4 and Extended Data Fig. 6 are drafted, the author's
  answers to flags F1 to F11 are applied, and the drafts are checked against
  the corrected Prism file and the assembled PDFs of 2026-11-06 evening.
  Nothing is pending from Claude Code. Two findings need the author: Fig. 3i's
  error bars are drawn the wrong way up (F13), and Fig. 4b's printed P values
  are not in the Prism file (F12). Extended Data Fig. 5 is skipped, unfinished.
  Two commits are UNPUSHED.
entry-point-for: the figure legends
relates-to: >
  docs/2026-11-03_e50_handoff.md is still the entry point for the panels and
  for the text ruling. This file does not replace it; it covers the legends
  only.
---

# Handoff, 2026-11-06: the legends are drafted, and two panels need the author

## 0. Resuming: start here

**Nothing is pending from Claude Code.** Everything below is the author's.

1. **Fig. 3i's error bars are drawn the wrong way up** (section 3). A Prism fix,
   not a repo one. **It changes which arms exclude zero.**
2. **Fig. 4b's three printed *P* values are not in the corrected Prism file**
   (section 3). Probably the separate mixed-effects file; say so and the flag
   is struck.
3. **The gaps in the drafts** (section 4): seven blot *n*, the ABT-737 exposure
   time, Fig. 3e's passaging detail, Fig. 4c's stain and scale bar, and the
   `[14 days]` carried from the old legend.
4. **The legends themselves**, to edit and paste. They are never pasted from
   here: the Doc is the author's and is read, never edited.
5. **The text ruling is still open** (section 5), and one part of it has moved:
   the author has edited the Results, so `pack_text.py` is behind the Doc.
6. **Two commits are unpushed** (section 6). `origin/main` is at `26a55f6`.

## 1. What this session did

- **Wrote the legends** for Fig. 2, Fig. 3, Fig. 4, Extended Data Fig. 4 and
  Extended Data Fig. 6, to the author's rules
  (`myc_mouse 2ebf646:figures/natmetab/LEGENDS.md`) and naming decisions of
  2026-11-06. Extended Data Fig. 5 is skipped; the author says it is not
  finished. Fig. 1 and Extended Data Fig. 1 to 3 are the author's own, already
  in the Doc.
- **Built the inventory and the statistics audit** behind them, and raised
  eleven flags. **The author answered all eleven**, inline in section 7 of the
  inventory, then recalculated the statistics in Prism and reassembled the
  figures.
- **Rewrote the legends to those answers**, against
  `NM_all_figs_corrected.prism` (2026-11-06 17:36) and the figure PDFs of 18:01
  to 18:03, and recorded panel by panel what changed, in section 9 of the
  inventory.
- **Matched every printed *P* to its own panel** by position on the page, and
  then to the analysis behind it. Every value on Fig. 2, Fig. 3 and Extended
  Data Fig. 4 is now accounted for. That check caught two errors in the
  inventory's own first reading, and the Fig. 3i finding.

## 2. Read in this order

1. **`docs/legends/README.md`** — the rules, where every number comes from, and
   the tools.
2. **`docs/legends/2026-11-06_panel_inventory_and_stats.md`** — every panel,
   what it draws, its *n*, its error bars and its test.
   - **Sections 1 to 6** are the morning's reading.
   - **Section 7** holds the eleven flags with **the author's replies inline**.
   - **Section 9 is what the corrected file and the new PDFs say**, and wins
     wherever it differs. F12 and F13 are at its end.
3. **`docs/legends/2026-11-06_legends_draft.md`** — the legends. `[...]` is a
   gap for the author; `[F12]` and `[F13]` are the two open flags.
4. **`docs/legends/tools/`** — three readers and a checker. They compute
   nothing of their own:
   - `prism_inventory.py` — per sheet: groups, *n*, mean, s.d., s.e.m.; per
     analysis: test, design, post-test, method, every adjusted *P* or *q* with
     the two groups it compares, and every *t* test with its *t*, df and *P*;
   - `pdf_error_bars.py` — every drawn two-point segment in a window of a page,
     in points, so a bar can be read against its own axis;
   - `pdf_pvalue_map.py` — every *P* and *q* the page prints, with the panel it
     sits in, from the panel letters' own positions;
   - `check_legends.py` — 68 checks over the drafts: every panel letter present
     and in order, *n*, centre, error bars and a test on every data panel,
     italic *P*, no interpretive verb, no title on an Extended Data figure, and
     every flag cited defined in the inventory. **68 of 68 pass.**

**The author's files are not in this repo.** They are in Drive, folder
`1qLbI6nQlLczmUBrmWyMvMmNPJoG6Fano`, and the tools take a path.

## 3. The two findings that need the author

### F13. Fig. 3i's asymmetric error bars are drawn the wrong way up

**Measured on all four arms**, from the panel's own vector geometry, with the
axis calibrated from its own ticks (0.0 at 99.28 pt, 97.95 pt per unit):

| arm | drawn | the mouse table's interval |
|---|---|---|
| NT + Pgc-1α | [−0.011, +0.075] | [+0.030, +0.116] |
| Bcl-xL | [−0.050, +0.049] | [−0.009, +0.128] |
| Bcl-xL + Pgc-1α | [+0.518, +0.664] | [+0.475, +0.622] |
| *p19arf* KO + Pgc-1α | [+0.036, +0.206] | **[−0.064, +0.105]** |

Every bar matches the table's two offsets **exchanged**, to within 0.002 of a
data unit, and none matches the table itself. The cause is the paste order:
Prism's `y_plus_minus` subcolumns are **mean, +error, −error**, and
`Fig3i_prism_mitopps_contrast.csv` writes `err_minus` first.

**It changes a conclusion.** For *p19arf* KO + Pgc-1α the interval includes zero
and is drawn excluding it; for NT + Pgc-1α it excludes zero and is drawn
including it. `myc_mouse 2ebf646:figures/natmetab/fig3i_prism_tables.R` warned
about the neighbouring error — symmetrising on `err_plus` — in the same arm, for
the same reason.

**The fix is in Prism:** exchange the two error subcolumns of
`3i Orthotopic OXPHOS mitoPPS`, the only sheet in the file entered as mean with
asymmetric error. **Extended Data Fig. 5b comes from the same script's CSV in
the same column order and needs the same check when that figure is finished.**

### F12. Fig. 4b's printed values are not in the corrected file

Fig. 4b prints *P* < 0.0001, *P* = 0.0002 and *P* < 0.0001. The corrected file's
mixed-effects analysis of that sheet has 36 comparisons (Tukey, every pair at
every time) and **its smallest adjusted *P* is 0.00094**; its fixed effects are
time 5.2 × 10⁻³⁶ and group 0.0030. **The likely answer is the one F4 already
gave** for Fig. 2j and Fig. 4d: a separate mixed-effects file. If so the legend
is already right and the flag is struck.

## 4. What the drafts still need from the author

- **Seven blot *n***, as `[n = ...]`: Fig. 2a, 2c, 2f; Fig. 3b, 3d, and the
  middle immunoblots of 3g and 3h. Two of them also carry
  `representative blot shown`, from the old legend.
- **Fig. 2j: the ABT-737 exposure time**, as `[... h]`.
- **Fig. 3a: `[14 days]`**, carried from the superseded legend for the
  continuously induced survivors. Confirm or correct.
- **Fig. 3e: the passaging detail**, as `[...]` — how far the cells were
  passaged before they stopped responding to Pgc-1α.
- **Fig. 4c: the stain and the scale bar** for the representative lungs.

## 5. The text ruling, still open

**`docs/2026-11-03_e50_handoff.md` section 0 holds the list.** Two changes since
it was written:

- **R13's citation is settled** (the author, 2026-11-06, option 3): R13 cites
  **Fig. 4f**, and **Extended Data Fig. 6c is cited from a different sentence**,
  which the author will write. Fig. 4's letters are settled by the assembled
  figure: mouse **a** to **e**, human **f** (heatmap), **g** (PUMA against
  BCL-xL), **h** (co-amplification), **i** and **j** (the forests); Extended
  Data Fig. 6 is **a** to **d**.
- **The author has since edited the Results**, so the sentences and the citation
  order in `docs/2026-11-03_e49_review_pack_src/pack_text.py` are behind the
  Doc. **A ruling must catch them up before that file is used for anything.**
  Until then it is a record of what was ruled, not of what the paper says.

The rest of the list is unchanged: Fig. 4h's sentence, the legends and letters
for E52 to E57's panels, E53's three open choices, v6's six proposed changes,
whether "no coloured text" reaches the white values on the darkest cells, and
whether E57's grey cue arrows should be black.

## 6. The record

| commit | what |
|---|---|
| `df34f81` | the inventory, the statistics audit, the draft legends, three tools |
| `1b3b38a` | the legends rewritten to the author's answers; section 9; F12 and F13; `pdf_pvalue_map.py` |

**Both are unpushed.** `origin/main` is at `26a55f6`. Push only when the author
asks.

## 7. What carries forward

- **The legend rules are the author's**, and a legend says what is on the panel
  and what is needed to read it: centre, error bars, *n*, the test, the *P*, and
  keys in words. No interpretation. Methods carry the rest.
- **Read a printed value's panel from its position**, not from the reading order
  of the PDF's text. Both of the inventory's own errors came from reading order.
- **Measure an error bar rather than assume it.** This session: Fig. 2e is a 95%
  confidence interval where every other cell panel is mean ± s.d., and Fig. 3i's
  bars are inverted. Neither is visible in the numbers alone.
- **A Prism sheet can carry more arms than its panel draws**, and a duplicated
  column counts as a comparison. Both were flagged, and the author re-ran the
  tests that needed it.
- **`myc_mouse` is read-only from here**, through
  `git -C /Users/gs/code/myc_mouse show <ref>:<path>`. It is where Fig. 3i's
  intervals are defined, and it carries an instruction for that legend.
