---
date: 2026-11-05
status: >
  PRE-RUN NOTE for E51, committed with the script BEFORE the author's run.
  Three dry runs, with every output redirected to a scratch folder; the second
  and third are identical. It gives the digests and checksums the author's run
  must reproduce.
relates-to:
  - docs/2026-11-05_E51_main_panel_ruling.md (86f38aa, the specification)
  - scripts/E51_puma_bclxl_main_panel.R
  - docs/2026-11-03_e49_prerun_digests.md (a63ef6c, the E49 digests E51 checks)
  - myc_mouse 2ebf646, figures/natmetab/fig1i_coupling.R, _style.R and FIGURE_RULES.md (read only)
---

# E51 pre-run note: ED-k as a main panel, both cohorts overlaid, with each P

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E51_puma_bclxl_main_panel.R`. It takes a few seconds.
3. **It needs only** `results/e49_signature_and_fig4g.rds`. It checks that object against E49's committed digests and stops if any differs.

**It writes:**
- `results/e51_puma_bclxl_main_panel.rds`;
- `outputs/tables/E51_slopes_and_p.csv` and `E51_assertions.csv`;
- the panel twice, as a PDF with live Helvetica and a PNG preview:
  - `outputs/figures/E51_puma_bclxl_against_priority_overlaid.{pdf,png}`;
  - `docs/figures/2026-11-05_E51_puma_bclxl_against_priority_overlaid.{pdf,png}`, which is tracked.

## 1. The dry runs

**Dry run 1 stopped at a check of mine**, after every gate on the data had passed.
- **The check:** that the PNG is the page's size at 600 dpi. It expected rounding, but R's PNG writer truncates, so 1,375.98 rows become 1,375. The check now allows less than one pixel.
- **Two fixes followed from reading its PDF:**
  - **Base `pdf()` writes the page box in whole points, truncating** (189.92 pt became 189). The page is now set in whole points, 0.001 pt over, so nothing is cut. Its height is Fig. 4g's PDF page, 165 pt.
  - **The x-axis title is written as a kerned `TJ` array,** not a plain string. The PDF read-back now decodes both forms.

**Dry runs 2 and 3 pass all 26 assertions, and are identical:**
- the eight saved elements' digests;
- the PNG, byte for byte;
- the PDF, with its two timestamps zeroed;
- both tables.

Between them, only a warning was silenced: `useBytes = TRUE` in one search of the PDF's bytes.

**The PDF was previewed with Apple's engine** (`qlmanage -t`), as Fig. 1's figure rules ask. It shows:
- the italic *P*;
- the hyphens in "TCGA-BRCA" and "SCAN-B", which are soft hyphens in the PDF;
- the minus signs on the axes;
- "× 10⁻⁹" with its exponent at full size.

## 2. What it computes

**Nothing is fitted.** The one new quantity is each cohort's *P*, declared in the ruling (question 1).
- **What it is:** the two-sided t-test of the drawn slope in its own model, ED-k's (`m3` on OXPHOS priority, the reference MYC score).
- **How it is computed:** from E49's saved residuals and the model's residual degrees of freedom.
- **It is checked two ways** (assertions 9 to 11):
  - the slope of residual on residual equals the stored coefficient exactly (Frisch-Waugh-Lovell);
  - the slope's t equals the partial correlation's t, to 5.3 × 10⁻¹⁵.

| cohort | n | residual df | slope [bootstrap 95% CI], as stored | partial r | t | *P* |
|---|---|---|---|---|---|---|
| TCGA-BRCA | 938 | 928 | −0.086 [−0.153, −0.022] | −0.079 | −2.40 | 0.017 |
| SCAN-B | 3,143 | 3,135 | −0.117 [−0.158, −0.074] | −0.103 | −5.79 | 7.8 × 10⁻⁹ |

**The panel prints *P* = 0.017 and *P* = 7.8 × 10⁻⁹.** The slopes and their intervals go in the legend.

## 3. The choices the ruling leaves open

The script's header gives each one with its reason.
- **The page is 190 × 165 pt (67.0 × 58.2 mm).**
  - **The height** is Fig. 4g's PDF page.
  - **The data area** is a square 46.75 mm side, as large as that height allows, measured as Fig. 1i measures it.
- **The style is Fig. 1i's,** from `myc_mouse` `2ebf646`. The theme, axes and *P* formatter are copied unchanged from `_style.R`.
  - **Points:** shape 21, size 1.4, outlined grey20.
  - **Lines:** 0.55, solid in both cohorts.
  - **Zero lines:** grey85.
  - **Type:** 5 pt for tick labels and *P*, 6 pt for axis titles.
- **Points are the ten decile means per cohort that E49 stored.** The lines pass through the origin, because the residuals have mean 0, and run over each cohort's decile range.
- **The axes are not on one scale:** x runs from −2 to +2 and y from −0.3 to +0.3.
- **The *P* labels are at the top right of the data area, in black.** Each begins with a key line in its cohort's colour.
  - Assertion 16 checks they overlap no point and no line.
  - There is no title, no other key and no *R*².
- **The colours:** TCGA-BRCA `#243642`, drawn over SCAN-B `#8FA3B0`.
- **The axis titles:** "OXPHOS priority (mitoPPS, SD), adjusted" and, on two lines, "*BBC3:BCL2L1* (SD)," and "adjusted". They are live text, for the author to edit in Illustrator.

**The legend draft, for the author:**

> BBC3 (PUMA) relative to BCL2L1 (BCL-xL) against OXPHOS priority (mitoPPS), each adjusted for MYC activity on the reference score, proliferation and intrinsic subtype, and for tumour purity and leukocyte fraction in TCGA-BRCA; TCGA-BRCA (n = 938), dark, and SCAN-B (n = 3,143), light. Each cohort is standardised within itself, and the two are not pooled. Points, means of the ten deciles of adjusted OXPHOS priority; lines, the fitted slopes, -0.09 (95% CI -0.15 to -0.02) in TCGA-BRCA and -0.12 (-0.16 to -0.07) in SCAN-B, with intervals from 2,000 paired bootstrap resamples of tumours. P, two-sided t-test of the slope.

## 4. The digests the author's run must reproduce

`digest::digest()` of each saved element, from dry runs 2 and 3. `run_info` is excluded.

| element | digest |
|---|---|
| `slopes` | `dec549b6e83e26be49aaaf9d49ce34a5` |
| `deciles` | `85136d7800b63e5dd8d69b064c90acd8` |
| `lines` | `ae9cdfec16f9f8682d3eea90a37f70fe` |
| `layout` | `18dcfa64172b163d15ef39dc0997f22d` |
| `legend_draft` | `8ccf7b2f800aad901bc5d4cc199f9288` |
| `assertions` | `a1ea712a0b7ddab7112106bb19e0b94c` |
| `figure_files` | `ce7502a4974d004f5ca3c976e9b05fd0` |
| `spec` | `d5ed00268258e5c8a578c331597525fe` |

**No assertion's label is built from an output path**, so the dry run's labels are the real run's (E50's lesson).

## 5. The checksums

| file | md5 |
|---|---|
| the PNG, in both folders | `49757a2995ed14c02bc9ff51d23435e3` |
| the PDF, in both folders, timestamps zeroed | `9be53025ea6a5bf8c4723b641d01da68` |
| `E51_slopes_and_p.csv` | `177544668c89c1e65ec137b0c62e4ab2` |
| `E51_assertions.csv` | `fdcfb3b29b8269e99216154143854823` |

**The PDF's masked checksum** is the md5 of the file with NUL bytes removed and the digits after `/CreationDate (D:` and `/ModDate (D:` replaced by a single 0. This is the script's own `.pdf_md5()`.

## 6. After the author's run

1. **Compare** the eight digests (section 4) and the checksums (section 5).
2. **Commit** the tracked figure pair in `docs/figures/` with a short result note.
3. **Then the text** follows, for a later ruling: R6's citation of the new panel, the panel's letter, and the legend.
