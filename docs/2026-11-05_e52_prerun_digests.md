---
date: 2026-11-05
status: >
  PRE-RUN NOTE for E52, committed with the script BEFORE the author's run.
  Six dry runs, with every output redirected to a scratch folder; runs 4, 5
  and 6 are identical. It gives the digests and checksums the author's run
  must reproduce.
relates-to:
  - docs/2026-11-05_E52_prism_size_ruling.md (9c2ca44, the specification)
  - docs/2026-11-05_E51_main_panel_ruling.md (86f38aa, which stands except where E52's ruling changes it)
  - scripts/E52_puma_bclxl_panel_26mm.R
  - docs/2026-11-05_e51_prerun_digests.md (1b0be30, the E51 digests E52 checks)
  - myc_mouse 2ebf646, figures/natmetab/_style.R, fig1h_bcl2_fade.R and FIGURE_RULES.md (read only)
---

# E52 pre-run note: E51's panel at the Prism graphs' size

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E52_puma_bclxl_panel_26mm.R`. It takes a few seconds.
3. **It needs only** `results/e49_signature_and_fig4g.rds`. It checks that object against E49's committed digests and stops if any differs.
   - It does not read E51's saved object.
   - It recomputes E51's quantities and checks them against E51's committed digests.

**It writes:**
- `results/e52_puma_bclxl_panel_26mm.rds`;
- `outputs/tables/E52_slopes_and_p.csv` and `E52_assertions.csv`;
- the panel twice, as a PDF with live Helvetica and a PNG preview:
  - `outputs/figures/E52_puma_bclxl_against_priority_26mm.{pdf,png}`;
  - `docs/figures/2026-11-05_E52_puma_bclxl_against_priority_26mm.{pdf,png}`, which is tracked.

**E51's files are not touched.** They stay as E51's verified record.

## 1. What changes from E51

| | E51 | E52 |
|---|---|---|
| data area | 46.75 mm square | **26 mm square** (the ruling) |
| page | 190 × 165 pt, Fig. 4g's height | **117 × 128 pt (41.3 × 45.2 mm)**, the panel's natural size |
| x title | "OXPHOS priority (mitoPPS, SD), adjusted" | **"OXPHOS priority" over "(mitoPPS, SD)"** |
| y title | "*BBC3:BCL2L1* (SD)," over "adjusted" | **"*BBC3:BCL2L1* (SD)"**, on one line |
| *P* labels | inside the data area, top right | **a key above the data area** (the ruling, question 1) |

**Unchanged:**
- the data, the *P* values, the colours, the points and the fitted lines;
- the axes: x from −2 to +2, y from −0.3 to +0.3;
- the type: 5 pt for tick labels and *P*, 6 pt for axis titles;
- the axes and ticks at 0.53 pt (the ruling, question 2);
- the legend draft (section 5).

## 2. The dry runs

**Dry run 1 stopped at a check of mine**, on the x title's line spacing.
- **The check:** that the two lines are 2.4 mm apart, baseline to baseline. Fig. 1's figure rules ask for "2.4 mm baseline to baseline (6 pt at 1.15)".
- **What happened:** R spaces lines 1.2 × the type size apart at lineheight 1. So lineheight 1.15 gave 8.28 pt (2.92 mm).
- **The fix:** the lineheight is now derived from the gap: 2.4 mm ÷ (1.2 × 6 pt) = 0.945. Fig. 1h sets its two-line title 2.4 mm apart (`fig1h_bcl2_fade.R`, `LEAD <- 2.4`).

**Dry run 2 passed every check, but SCAN-B's key line sat high.**
- **Why:** each row was centred on its bounding box. The exponent in SCAN-B's row makes that box taller, so its key line sat about 0.4 mm above the middle of "SCAN-B". E51's labels had the same offset, inside the data area.
- **The fix:**
  - the rows are now set on their baselines, 2.8 mm apart;
  - each key line is level with the middle of its row's capitals;
  - two new checks read both back from the PDF.

**Dry run 3 passed all 38 checks.** Then three checks' records were tidied. Each had stored a margin, not a deviation, under `max_abs_dev`, one of them negative; the margin now goes in `detail`. No drawing changed.

**Dry runs 4 and 5 are identical:**
- the eight saved elements' digests;
- the PNG, byte for byte;
- the PDF, with its two timestamps zeroed;
- both tables.

**Dry run 6 sourced the script in one R session, as Positron does,** then scanned the workspace.
- **The scan:** 210 objects, none holding invalid UTF-8 (E51's lesson).
- **It matched dry runs 4 and 5** in every digest and checksum.

**The PDF was previewed with Apple's engine** (`qlmanage -t`), as Fig. 1's figure rules ask. It shows:
- the key above the data area, TCGA-BRCA's row on top, each key line at the middle of its row's capitals;
- the x title on two lines, and no "adjusted";
- the italic *P*, the hyphens in "TCGA-BRCA" and "SCAN-B" (soft hyphens in the PDF), the minus signs on the axes, and "× 10⁻⁹" with its exponent at full size.

## 3. What the checks read back from the PDF

**E52 keeps E51's checks of the data and the *P* values**, and checks that the slopes, points, lines and legend draft are E51's verified ones, by digest.

**These checks read the PDF the author will open in Illustrator:**
- **the data area is 26 mm square.** Its two axis lines and two zero lines measure 73.70 pt and start on its edges, and no straight line is longer. The largest deviation is 0.0008 pt;
- **the key's rows lie above the data area,** TCGA-BRCA's on top. The lower baseline is 1.2 mm above the data area, and the two are 2.8 mm apart;
- **each key line starts at the data area's left edge,** level with the middle of its row's capitals;
- **the x title's two lines lie below the data area,** 2.4 mm apart;
- **no text contains "adjusted";**
- **the page is the panel's natural size,** rounded up to whole points.

## 4. The choices the ruling leaves open

The script's header and comments give each one with its reason. Each is the author's to reverse.
- **The key is left-aligned with the data area,** starting at the y axis, as in the ruling's preview. It is 22.1 mm wide, within the 26-mm axis.
- **Its rows:**
  - TCGA-BRCA on top, as in E51;
  - baselines 2.8 mm apart, the lower one 1.2 mm above the data area;
  - a 2.5-mm key line in the cohort's colour, 0.8 mm before black 5-pt text.
- **Above the key,** the theme's 1-mm top margin is kept. The key makes the panel 5.3 mm taller.
- **The x title breaks after "priority":** "OXPHOS priority" (15.9 mm) over "(mitoPPS, SD)" (13.7 mm), centred under the axis.
- **The y title stays on one line**, 18.6 mm along the 26-mm axis.
- **The page is 117 × 128 pt (41.27 × 45.16 mm).** Base `pdf()` writes the page box in whole points, so the natural size, 41.26 × 45.16 mm, is rounded up.

## 5. The legend draft: E51's, unchanged

It already states the adjustment that leaves the axis titles: "each adjusted for MYC activity on the reference score, proliferation and intrinsic subtype, and for tumour purity and leukocyte fraction in TCGA-BRCA". E52 checks it against E51's digest.

> BBC3 (PUMA) relative to BCL2L1 (BCL-xL) against OXPHOS priority (mitoPPS), each adjusted for MYC activity on the reference score, proliferation and intrinsic subtype, and for tumour purity and leukocyte fraction in TCGA-BRCA; TCGA-BRCA (n = 938), dark, and SCAN-B (n = 3,143), light. Each cohort is standardised within itself, and the two are not pooled. Points, means of the ten deciles of adjusted OXPHOS priority; lines, the fitted slopes, -0.09 (95% CI -0.15 to -0.02) in TCGA-BRCA and -0.12 (-0.16 to -0.07) in SCAN-B, with intervals from 2,000 paired bootstrap resamples of tumours. P, two-sided t-test of the slope.

## 6. The digests the author's run must reproduce

`digest::digest()` of each saved element, from dry runs 4, 5 and 6. `run_info` is excluded.

| element | digest | the same as E51's |
|---|---|---|
| `slopes` | `dec549b6e83e26be49aaaf9d49ce34a5` | yes |
| `deciles` | `85136d7800b63e5dd8d69b064c90acd8` | yes |
| `lines` | `ae9cdfec16f9f8682d3eea90a37f70fe` | yes |
| `layout` | `d52c33a31639f3933bd492cffba2c3d6` | no |
| `legend_draft` | `8ccf7b2f800aad901bc5d4cc199f9288` | yes |
| `assertions` | `f71e9868d0def84e5eee6506289f993a` | no |
| `figure_files` | `52fb33386987f5c144e54b9eb9af384f` | no |
| `spec` | `7b5056aab2607c89318d352ea690139e` | no |

**No check's label is built from an output path**, so the dry run's labels are the real run's (handoff lesson 3).

## 7. The checksums

| file | md5 |
|---|---|
| the PNG, in both folders | `1e031002ab3d5af17f8dffd84771c115` |
| the PDF, in both folders, timestamps zeroed | `d78da647dec8172ecf16893fe8e8cfec` |
| `E52_slopes_and_p.csv` | `177544668c89c1e65ec137b0c62e4ab2`, byte-identical to E51's table |
| `E52_assertions.csv` | `9e4c62cdd94742614b03e40ee0bc685d` |

**The PDF's masked checksum** is E51's: the md5 of the file with NUL bytes removed and the digits after `/CreationDate (D:` and `/ModDate (D:` replaced by a single 0.

## 8. After the author's run

1. **Compare** the eight digests (section 6) and the checksums (section 7).
2. **Commit** the tracked figure pair in `docs/figures/` with a short result note.
3. **Update** CLAUDE.md, the handoff and E51's result note: E52's panel replaces E51's for the paper.
4. **The text still follows, for a later ruling:** R6's citation of the new panel, the panel's letter, and the legend.
