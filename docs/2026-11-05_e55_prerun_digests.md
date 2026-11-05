---
date: 2026-11-05
status: >
  PRE-RUN NOTE for E55, committed with the script BEFORE the author's run.
  Five dry runs, with every output redirected to a scratch folder; runs 3, 4
  and 5 are identical. It gives the digests and checksums the author's run
  must reproduce.
relates-to:
  - docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md (dd5004a, the specification)
  - scripts/E55_raw_correlations_heatmap.R
  - docs/2026-11-03_e49_prerun_digests.md (a63ef6c, the E49 digests E55 checks)
  - scripts/E50_fig4g_in_fig1g_style.R (Fig. 4g's grammar, whose constants and helpers E55 copies and checks)
---

# E55 pre-run note: Extended Data Fig. 6c, the unadjusted correlations

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E55_raw_correlations_heatmap.R`. It takes a few seconds.
3. **It needs** `results/e49_signature_and_fig4g.rds` and `results/e48_separation_and_fig4g.rds`. It checks both against their committed digests and stops if any differs.

**It writes:**
- `results/e55_raw_correlations_heatmap.rds`;
- `outputs/tables/E55_raw_correlations.csv` and `E55_assertions.csv`;
- the panel twice, as a PDF with live Helvetica and a PNG preview:
  - `outputs/figures/E55_raw_correlations_heatmap.{pdf,png}`;
  - `docs/figures/2026-11-05_E55_raw_correlations_heatmap.{pdf,png}`, which is tracked.

## 1. What it draws

**E49's unadjusted Spearman correlations,** the left side of ED-h, in Fig. 4g's grammar.
- **Rows:** Fig. 4g's fifteen genes in its four groups and order, then one row for the signature (S6, "Signature") at the bottom. E49's S5 and S7 are left out.
- **Columns:** all seven measures (the ruling, question 1), in two cohort blocks, TCGA then SCAN-B, as E49 drew them:
  - OXPHOS level and OXPHOS priority;
  - the four MYC activity scores, by Fig. 6a's names, with the share of their genes in proliferation gene sets;
  - proliferation.
- **The heatmap takes Fig. 4g's size,** 58.1 × 51.8 mm. Its rows keep 4g's height and gaps, and its fourteen columns share 4g's width, at 4.10 mm each. **No values are printed.**
- **The colours are Fig. 4g's:** Fig. 1g's palette at ±0.6, with 4g's colour bar.
- **The page is 87.3 × 80.9 mm.** Its width is within 4g's 87.6 mm; the vertical column labels make it taller.
- **The data:** E49's `raw_correlations`, checked by digest and against the values E49 drew in ED-h, all 224. Nothing is computed.

## 2. The dry runs

**Dry run 1 stopped at a check:** the PDF's page box was read with a leading empty field. E53 and E54's reader drops it; this one did not. Fixed.

**It also showed the page at 90.4 mm wide,** wider than 4g's.
- **Why:** the colour bar's title, "Spearman rho" (Fig. 6a's words), is 11.2 mm wide on one line, and 4g's "SD per SD" is 8.4 mm.
- **The fix:** the title is now set on two lines, "Spearman" over "rho", at 8.2 mm. A check now holds the page to 88 mm.

**Dry run 2 passed all 22 checks.** Dry runs 3 and 4 are identical: the seven saved elements' digests, the PNG in both folders, the PDF with its timestamps zeroed, and both tables.

**Dry run 5 sourced the script in one R session, as Positron does:** 179 objects, none holding invalid UTF-8. It matched dry runs 3 and 4 in all 13 digests and checksums.

**The PDF was previewed with Apple's engine** (`qlmanage -t`).

## 3. What the checks read back from the PDF and the PNG

- **The cells:** 224, each 4.10 × 2.98 mm, in its colour and place. The colour bar is a gradient.
- **The labels:** every one is live Helvetica (no image, no embedded or Type 3 font), at its size, angle and place.
  - The gene symbols are oblique, and nothing else is.
  - The fourteen column labels stand at 90°.
- **The page** is the panel's size, in whole points.
- **The PNG** is the panel at 600 dpi, and every cell shows its colour exactly.

The layout checks also hold the heatmap to 4g's width and height, the rows to 4g's gaps, and the labels clear of one another, of the cells and of the bar. The 54 lines copied from E50 (the palette and the text helpers) are verbatim.

## 4. The choices the ruling leaves open

The script's header and comments give each one with its reason. Each is the author's to reverse.
- **The columns are grouped by cohort,** TCGA then SCAN-B, as E49 drew them. 4g groups by measure, but with seven measures that would need fourteen cohort labels above columns 4 mm wide.
- **The column labels stand vertically,** as in E49, reading upwards from just above the cells.
- **The MYC scores carry Fig. 6a's names,** with their proliferation shares, so the two panels match. E49's labels had no percentages.
- **The colour bar's title** is "Spearman" over "rho" (section 2).
- **Two correlations exceed ±0.6 and take the end colour:** BCL2A1 with the CollecTRI MYC regulon, +0.64 in TCGA-BRCA and +0.61 in SCAN-B. Raising the scale would change the colour code, which the ruling keeps as 4g's. The legend draft says so.
- **PDF and PNG only.** E50 also wrote an SVG for 4g; none is written here.

## 5. The legend draft (Claude Code's, for the author)

> Unadjusted Spearman correlations of fifteen BCL-2 family genes and the signature with OXPHOS (level and priority), four MYC activity scores and proliferation, in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143). Each MYC activity score is labelled with the percentage of its genes in proliferation gene sets. Rows and colours as in Fig. 4g; colours saturate at +/-0.6, which 2 correlations exceed (BCL2A1 with the CollecTRI MYC regulon, TCGA, +0.64; BCL2A1 with the CollecTRI MYC regulon, SCAN-B, +0.61).

## 6. The digests the author's run must reproduce

`digest::digest()` of each saved element, from dry runs 3, 4 and 5. `run_info` is excluded.

| element | digest |
|---|---|
| `values` | `79b613997fc59763950a468df47b4775` |
| `beyond_scale` | `3f070b62c469c06c9cb3f0208c6882fc` |
| `layout` | `347df987297cdc20d14f51b04a0dbae0` |
| `legend_draft` | `e8719dd9d46dfaac01583108fa45039e` |
| `assertions` | `50a2031b5f4fc1c7030e5c8b25fc2e04` |
| `figure_files` | `df385777c612c0441e170f61ba0a1957` |
| `spec` | `23d63f602f46f55e7f99148334355fac` |

**No check's label is built from an output path**, so the dry run's labels are the real run's.

## 7. The checksums

| file | md5 |
|---|---|
| the PNG, in both folders | `26c617ab045e2e8dc93f278a0101f14e` |
| the PDF, in both folders, timestamps zeroed | `0e44eb30758b71583330a8d1a6ddbfff` |
| `E55_raw_correlations.csv` | `2379fd28f60252b732db645462b918b7` |
| `E55_assertions.csv` | `f5c8ae47c34218ac105fb33808f24da5` |

## 8. After the author's run

1. **Compare** the seven digests (section 6) and the checksums (section 7).
2. **Commit** the tracked figure pair with a short result note.
