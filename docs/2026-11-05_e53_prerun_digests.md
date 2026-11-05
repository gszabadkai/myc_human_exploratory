---
date: 2026-11-05
status: >
  PRE-RUN NOTE for E53, committed with the script BEFORE the author's run.
  Ten dry runs, with every output redirected to a scratch folder; runs 8, 9
  and 10 are identical. It gives the digests and checksums the author's run
  must reproduce.
relates-to:
  - docs/2026-11-05_E53_half_width_ed_panels_ruling.md (7325ff1, the specification)
  - scripts/E53_half_width_ed_panels.R
  - docs/2026-11-03_e48_prerun_digests.md (1ccf805, the E48 digests E53 checks)
  - docs/2026-11-03_e48_result.md (E48, whose ED-a and ED-b these panels redraw)
  - myc_mouse 2ebf646, figures/natmetab/_style.R and FIGURE_RULES.md (read only)
---

# E53 pre-run note: Extended Data Fig. 6a and 6b at half width

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E53_half_width_ed_panels.R`. It takes a few seconds.
3. **It needs only** `results/e48_separation_and_fig4g.rds`. It checks that object against E48's committed digests and stops if any differs.

**It writes:**
- `results/e53_half_width_ed_panels.rds`;
- `outputs/tables/E53_assertions.csv`;
- each panel twice, as a PDF with live Helvetica and a PNG preview:
  - **6a:** `outputs/figures/E53_score_correlations_half_width.{pdf,png}` and `docs/figures/2026-11-05_E53_score_correlations_half_width.{pdf,png}`;
  - **6b:** `outputs/figures/E53_bcl2l1_against_bbc3_half_width.{pdf,png}` and `docs/figures/2026-11-05_E53_bcl2l1_against_bbc3_half_width.{pdf,png}`.

  The `docs/figures` copies are tracked.

**E48's files are not touched.** They stay as E48's verified record.

## 1. What changes from E48

| | E48 (ED-a, ED-b) | E53 (Fig. 6a, 6b) |
|---|---|---|
| width | 178 mm each | **6a 85.4 mm, 6b 86.4 mm**, both under the ruling's 88 |
| page | 178 × 73.7 mm and 178 × 62.2 mm | **242 × 184 pt (85.4 × 64.9 mm)** and **245 × 167 pt (86.4 × 58.9 mm)**, each the natural size |
| type | 5.5 to 6.5 pt; titles bold | **5 pt; titles 6 pt, plain** (Fig. 1's figure rules) |
| PDF | `cairo_pdf` | **base `pdf()`, live Helvetica**, for Illustrator |
| 6a's layout | three blocks side by side | **the two triangles side by side, 4.6-mm cells; the neoadjuvant values in a row of tiles under the colour key** (question 1) |
| 6b's arrows | "lower/higher expression", "BBC3/BCL2L1 steeper" | **removed** |
| 6b's rows | "TCGA · level" and so on | **"TCGA: OXPHOS level"**, "TCGA: OXPHOS priority", "SCAN-B: OXPHOS level", "SCAN-B: OXPHOS priority" (question 3) |
| 6b's key | "BCL2L1, without the proliferation term" and five more | **"BCL2L1", "BBC3", "Difference", each beside "proliferation adjusted"** (question 4) |
| 6b's x axes | −0.6 to 0.6, and −0.3 to 0.3 | **−0.1 to 0.6, and −0.1 to 0.3**, on one scale (question 2) |

**Unchanged:**
- the data: every value, estimate and interval is E48's, from its saved display frames;
- 6a's bands, greys and inks, its labels and titles, and its column labels at 40°;
- 6b's marks: colours, shapes, sizes, the rows' offsets, the intervals, the zero line and the connectors;
- 6b's axis titles, in E48's words.

## 2. The dry runs

**Dry run 1 stopped at a check: 6b was 90.3 mm wide.**
- **Why:** patchwork adds its own outer margin, 5.5 pt (1.9 mm) a side, from its default theme. E52's page carries the same margin.
- **The fix:** 6b sets that margin to 0. The two plots' own margins frame the page. It is now 86.4 mm wide.

**Dry runs 2 to 5 stopped at checks of mine.** Each was a fault in a check or in the drawing, and none in the data.
- **Run 2:** the exact-strings check compared a named vector with an unnamed one. `gsub()` keeps names. Both sides are now unnamed.
- **Run 3:** the PDF writes text matrices to two decimals, so a rotated 5-pt label reads back as 4.997 pt. The size checks now allow 0.01 pt.
- **Run 4:** the row-label check also picked up five column labels, which share their words. It now reads only unrotated text.
- **Run 5:** each row label sat one descent below its target.
  - **Why:** for a plain string, grid's `vjust = 0` sets the *baseline* at `y`. My finding from E52, that it sets the lowest point there, holds only for plotmath expressions.
  - **The fix:** labels are drawn at their baselines. The rotated labels' extents now count their descenders.

**Dry run 6 passed all 42 checks, but its preview showed a fault: the key's "Difference" row had no symbols.**
- **Why:** since ggplot2 3.5, a layer draws key symbols only for the levels in its own data, unless `show.legend = TRUE`. The key is drawn from the left panel, which holds no differences.
- **The fix:** the left panel's points set `show.legend = TRUE`, and a 43rd check confirms the key's six entries.

**Dry run 7 passed all 43 checks, and its previews were right.**

**Dry runs 8 and 9 are identical:** the six saved elements' digests, both PNGs in both folders, both PDFs with their timestamps zeroed, and the table.

**Dry run 10 sourced the script in one R session, as Positron does,** then scanned the workspace.
- **The scan:** 278 objects, none holding invalid UTF-8 (E51's lesson).
- **It matched dry runs 8 and 9** in all 15 digests and checksums.

**Both PDFs were previewed with Apple's engine** (`qlmanage -t`), as Fig. 1's figure rules ask.

## 3. What the checks read back from the PDFs

**6a:**
- the PDF shows exactly 73 strings: the two titles, the 6 row and 12 column labels, the 45 values, the key and the neoadjuvant row;
- the four titles are 6 pt and everything else 5 pt;
- the twelve column labels are at 40°, and no other text is rotated;
- 45 cells measure 4.6 mm, and the key's three squares 3 mm;
- each row label's capitals are centred on its row.

**6b:**
- the PDF shows exactly 20 strings: the four row labels, the ticks, the two axis titles and the six key entries. So none of the arrows' words, and none of E48's key texts;
- the axis titles are 6 pt and everything else 5 pt;
- the rows run from TCGA's OXPHOS level at the top to SCAN-B's priority;
- the two x axes measure 37.8 and 21.6 mm, 54 mm per unit, so both panels share one scale;
- each zero line spans its panel, 0.1 of a unit from its axis's start;
- each axis title's two lines are 2.4 mm apart, and the key lies below them.

**Both:** the page is at most 88 mm wide, the text is live Helvetica and not embedded, and the two copies of each file are identical.

## 4. The choices the ruling leaves open

The script's header and comments give each one with its reason. Each is the author's to reverse.

**6a:**
- **The colour key and the neoadjuvant row start at the page's left margin,** where the longest row label starts.
- **The neoadjuvant row:** each cohort's name sits to the left of its tile, and the value is inside the tile, as in the triangles.
- **The titles are plain 6 pt.** E48's were bold 6.5 pt.
- **The white line between cells is 0.35.** E48's was 0.5, at twice the size.
- **Every label is set on its baseline,** with its capitals centred on its row, tile or square.

**6b:**
- **Both panels share one scale,** 54 mm per unit, so the difference reads at the coefficients' scale.
- **Ticks are labelled at 0, 0.2 and so on, with unlabelled ticks at the odd tenths.** So −0.1, where each axis starts, is a tick without a label, as Fig. 1's axes end. The ruling's preview sketched "−0.1" labelled. Labelling it is a one-line change.
- **The key has three rows of two,** aligned with the left panel's data area. The preview's two closing lines ("open = …", "filled = …") explained the symbols to the author and are not drawn. The legend draft says what they mean.
- **The axis titles break after "MYC activity" and after "BBC3".**
- **E48's marks are kept at their full-width sizes:** points 1.5, intervals 0.4, connectors 0.3 at 60%, and the zero line grey45 at 0.3. The axes take Fig. 1's 0.53-pt line; E48's were 0.85 pt.

## 5. The legend drafts

E48's, with one word changed in 6a's: the neoadjuvant values are now "Below", not "Right". The script checks that nothing else differs.

> **6a:** Spearman correlations among the two OXPHOS scores, four MYC activity scores and proliferation in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143), shaded by band (below 0.60; 0.60 to 0.80; 0.80 or above). Each MYC activity score is labelled with the percentage of its genes in proliferation gene sets (Hallmark E2F targets and G2M checkpoint). Below, the reference MYC activity score against proliferation in three neoadjuvant cohorts.

> **6b:** OXPHOS coefficients for BCL2L1 and BBC3 with MYC activity held fixed, without (open) and with (filled) the proliferation term, in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143), and the difference between them, with 95% percentile intervals from 2,000 paired bootstrap resamples of tumours. Bars on the coefficients, model-based 95% confidence intervals.

## 6. The digests the author's run must reproduce

`digest::digest()` of each saved element, from dry runs 8, 9 and 10. `run_info` is excluded.

| element | digest |
|---|---|
| `frames` | `1148ab9713dfb71fa30dce2e544e1a20` |
| `layout` | `7d45fc213533a3818beb1457f6d0301d` |
| `legend_drafts` | `dea4c3726b9a76c00f0bca3e26a02ea7` |
| `assertions` | `3307e23f3963430ece145ec1b062ddd0` |
| `figure_files` | `e58898ff3ed27b715581e54bea10237a` |
| `spec` | `413ecd99d8d50211dab02aca33a83723` |

**No check's label is built from an output path**, so the dry run's labels are the real run's (handoff lesson 3).

## 7. The checksums

| file | md5 |
|---|---|
| 6a's PNG, in both folders | `270fb8f9222e518debaccf530ff38a80` |
| 6a's PDF, in both folders, timestamps zeroed | `2d877145a51b6015a540e0d579ebfc6f` |
| 6b's PNG, in both folders | `80c63bbb0c78c6a38deebc6ce05b787c` |
| 6b's PDF, in both folders, timestamps zeroed | `f15e91a8a00b83c2506464e70bfe1175` |
| `E53_assertions.csv` | `151b53518b290c49db9f6e5539bc9fd2` |

**The PDFs' masked checksums** are E51's: the md5 of the file with NUL bytes removed and the digits after `/CreationDate (D:` and `/ModDate (D:` replaced by a single 0.

## 8. After the author's run

1. **Compare** the six digests (section 6) and the checksums (section 7).
2. **Commit** the four tracked figures in `docs/figures/` with a short result note.
3. **Update** CLAUDE.md and the handoff: E53's panels replace E48's ED-a and ED-b for the paper, as Extended Data Fig. 6a and 6b.
