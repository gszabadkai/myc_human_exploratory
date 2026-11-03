---
date: 2026-11-03
status: >
  E50 WRITTEN and NOT YET RUN BY THE AUTHOR. These digests are recorded and
  committed BEFORE that run. E50 is display-only: it fits nothing and reads
  E49's verified object. Verification follows the author's run.
relates-to:
  - docs/2026-11-03_E50_display_amendment.md (6f5bf4d, the specification, committed alone)
  - scripts/E50_fig4g_in_fig1g_style.R
  - docs/2026-11-03_e49_prerun_digests.md (a63ef6c, the digests E50's input gate checks)
---

# E50: digests recorded before the author's run

## 0. Running it

**Restart R, then source `scripts/E50_fig4g_in_fig1g_style.R`.**
- It needs `results/e49_signature_and_fig4g.rds`, the verified E49 object, and nothing else.
- A run takes about 15 seconds.
- **A fresh session matters** (section 2.2). An earlier `cairo_pdf` in the same session can move the text in a later cairo PNG, and the PNGs are compared byte for byte.

## 1. What the dry runs established

**Every output of every run went to the session scratchpad.** Nothing was written in the repo.

| run | what happened |
|---|---|
| first attempt | stopped while reading its own SVG back (2.1) |
| dry run 1 | **112 of 112 assertions passed.** The panels were looked at; the PNG preview's colours were wrong (2.2) |
| dry runs 2 and 3 | **114 of 114.** In each run, ED-g's two PNG copies differed (2.2) |
| dry runs 4 and 5, final version | **116 of 116 in each; no warnings.** 12 of 12 saved elements agree. 24 of 24 compared files are identical: the PNGs, SVGs and tables byte for byte, and the three base-`pdf()` PDFs after their timestamps are masked |

**The elements excluded from comparison:** `run_info`, which holds the build time, the bitmap type and the R version.

### 1.1 The input gate

- **E49's object is the verified one.** Six of its elements are read: `all_readouts_m3`, `placement_rule`, `raw_correlations`, `displays`, `legend_draft` and `assertions`. Each matches the digest committed at `a63ef6c`, before E49's run.
- **E49's 192 assertions all passed.**
- **4g is E49's branch:** three blocks, all from the least-entangled model.
- **The rows are Fig. 1g's groups and order,** as E49 drew them (E49 declaration 9.1).
- **Every value E50 draws is one E49 stored.** For 4g, ED-e and ED-f the check is against E49's own display tables, the signature's grey cells included. For ED-g and ED-h it is against `all_readouts_m3` and `raw_correlations`.

### 1.2 What the assertions check

**The amendment's three assertions (item 5), on the layout and on the files themselves:**
- **Every printed value equals its stored value, rounded.**
  - On the layout, every label is checked: signed, one decimal, "0" when it rounds to zero, and within 0.05 of its value.
  - Each value is then read back from the SVG and from the PDF.
- **The palette's three anchors equal item 2's colours.** These are Fig. 1g's source, `div_cols` in `myc_mouse` `2ebf646` `figures/natmetab/_style.R` (section 3).
  - The scale returns the anchors at -L, 0 and +L.
  - Every cell's colour equals the colour that Fig. 1g's own scale call, `ggplot2::scale_fill_gradientn`, gives its value.
  - The colour bar's gradient carries the palette.
- **The cell size, gaps and fonts match item 3:**
  - **cells:** 96 per panel, each 26.71 x 8.46 pt, with no border;
  - **gaps:** none within a group; 0.8 mm (2.27 pt) between groups and between blocks; 1.6 mm before the signature;
  - **type:** Helvetica 5 pt, with 6 pt for the score headers only, and oblique for the gene symbols only.

  Each is checked on the layout and read back from the SVG and the PDF.

**And:**
- **The files.**
  - **PDF:** the text is live Helvetica, with no image and no Type 3 font.
  - **SVG:** it has no `<path>`, so no letter is outlined.
  - **PNG:** each cell's pixels equal its colour exactly.
  - Each SVG and PNG is identical in `outputs/figures` and `docs/figures`.
- **The legibility.**
  - No two labels overlap, and no label touches a cell or the colour bar.
  - Every value fits its cell and reads at a contrast of 4.5 or more; the lowest is 4.63.
  - Everything lies inside the panel.
- **ED-g and ED-h, re-coloured.**
  - Every cell sits where E49 drew it, with E49's text, in its actual colour.
  - Nothing is drawn in E49's grey `#A6A6A6`.
- **The legends.**
  - Each replaced sentence occurred exactly once in E49's drafts.
  - The author's sentence stands in 4g's, ED-g's and ED-h's legends.
  - No legend mentions grey, a dot or a marked version.
  - Nothing else changed: each edited legend is reconstructed into E49's exactly.

## 2. Changes before the final dry runs

### 2.1 A stop before a run completed

1. **Reading the SVG back.** `xml2::xml_parent()` returns each parent once, not once per text, so the reader built a short column. It now asks each text node for its own parent.

### 2.2 After the panels were looked at

**The panels were looked at for layout, legibility and colour only.**
- The PDF was also rendered by macOS's own PDF renderer, and the SVG by QuickLook. Both draw live text in Helvetica, in the same places.
- Viewing showed the values. Every one is a quantity E49 stored, verified and reported. **Nothing in content or layout was changed because of a value.**

1. **The PNG previews now go through cairo.**
   - **The problem:** R's default bitmap device here, quartz, colour-manages the PNG and darkens the palette. A cell meant to be `#624E3F` was written as `#4F3D30`.
   - **The fix:** cairo writes every colour exactly. An assertion now reads one pixel per cell back from each PNG.
   - **This affects E49 too.** Its PNG previews went through quartz, so they show the same darkening. Its PDFs are unaffected.
   - **This may bear on the author's measured anchors.** The values measured from the Fig. 1 PDF (#462915, #F9F9F9, #318866) are darker than Fig. 1g's source values (#4A3525, #FAFAFA, #2A8A6D). That would fit a colour-managed reading, but it is a guess and is not tested here.
2. **All ED-g and ED-h PNGs are written before any cairo PDF.**
   - **The problem:** in one R session, a cairo PNG drawn after a cairo PNG and a `cairo_pdf` sets ggplot's two-line cell text higher than the first PNG did. So ED-g's copy in `docs/figures` differed from its copy in `outputs/figures`.
   - **Reproduced in isolation:** PNG, `cairo_pdf`, PNG changes the second PNG. PNG then PNG does not, nor does a `cairo_pdf` first, nor a base `pdf()`.
   - **The panels drawn by grid in Helvetica were not affected** (4g, ED-e and ED-f).
   - **Two assertions now require each PNG and SVG to be identical in both folders.**

## 3. The palette, and the choices E50 makes

**Item 2 puts Fig. 1g's source first, and the source is in the mouse repo.** So E50 uses its palette, not the measured values:

| anchor | Fig. 1g's source (used) | measured from the Fig. 1 PDF (not used) | largest channel difference |
|---|---|---|---|
| negative limit | `#4A3525` | `#462915` | 16 |
| zero | `#FAFAFA` | `#F9F9F9` | 1 |
| positive limit | `#2A8A6D` | `#318866` | 7 |

**Fig. 1g's scale interpolates in CIELAB, not in RGB.** It is ggplot2's `scale_fill_gradientn`, and E50 follows it. Item 2's "linear in RGB" belongs to the measured branch.

**Choices the amendment leaves open, each the author's to reverse:**
- **One colour scale for 4g, ED-e and ED-f: ±0.6**, by item 2's rule over all three panels' values.
  - So a colour means the same value in each panel.
  - On its own values ED-e would be ±0.5.
- **ED-g's own scale is ±0.9.** Its "proliferation minus OXPHOS" cells, no longer grey, reach -0.88.
- **ED-h's scales are ±0.7** for the raw correlations and **±0.6** for the adjusted coefficients.
- **The gap before the signature row** is twice the group gap, 1.6 mm.
- **Value ink is black, or white where white has the higher contrast.**
  - Fig. 1g needs only black, because its limits are ±2.
  - At ±0.6 the brown end is too dark for black, so 9 of the 288 cells take white.
- **Type follows item 3's two sizes.**
  - 6 pt for the score headers only.
  - 5 pt for everything else. That includes the group names and the bar's title, which Fig. 1g sets at 6 pt.
- **The score header sits 2.4 pt above the cohort labels.** That is Fig. 1g's strip margin; Fig. 1g has no second header tier.
- **No cell borders** (item 3). Fig. 1g's tiles have a 0.64 pt white border.
- **ED-h is drawn once, unmarked**, because there are no dots.
- **ED-g and ED-h keep E49's layout, labels, text and PDF device** (`cairo_pdf`). They take the same ink rule.

**From Fig. 1g's code, measured on 2026-11-03.** Its panel was rebuilt with dummy values and written to an uncompressed PDF, and these were read from it:
- the cell size;
- the 0.8 mm gap;
- the label offsets: gene symbols 1.2 pt to the left, column labels 1.2 pt above, group names 2.4 pt to the right;
- the bar's geometry: 2 x 18 mm, ticks 1.70 pt long in grey30, tick labels 2.99 pt from the bar, and the key 1 mm past the wider of the first two group names.

**The panels' size:** 87.6 x 58.2 mm. All three are identical in size and positions.

## 4. Digests, recorded before the author's run

| element | digest |
|---|---|
| `values` | `78a99ee849f5c523cef42916e7aa8707` |
| `layouts` | `eb9e716586176f0b0f8c738f0e078dd4` |
| `palette` | `92fd15d514d0e4618dbe72f120297897` |
| `axes` | `bb248ae6fc4e3056523b8eabc0191398` |
| `ed_g` | `6d56c02405fd322a260fed766b4a7d25` |
| `ed_h_raw` | `b18eae8e8d50531d58b7e253072f32bb` |
| `ed_h_adjusted` | `ad41c2a89b756397c3011cd35876f40e` |
| `output_checks` | `fdbf873f2dd24283b84c0b58b3c612f5` |
| `legend_draft` | `eb6455f5be282070224f8fc223001666` |
| `assertions` | `8bd997a016e917ad0858396e0565cea9` |
| `figure_files` | `4817b8704dbf194426028e2b1b6d0b31` |
| `spec` | `fb0e8c586ddcafddb1ee969cfada9aed` |

## 5. Verifying the author's run

```r
x  <- readRDS("results/e50_fig4g_in_fig1g_style.rds")
el <- setdiff(names(x), "run_info")
sapply(el, function(k) digest::digest(x[[k]]))   # compare with section 4
```

**The files are compared with section 6:**
- **the PNGs, SVGs and the two tables,** byte for byte;
- **the three base-`pdf()` PDFs,** byte for byte after the digits of `/CreationDate` and `/ModDate` are replaced by zeros:

  ```r
  pdf_md5 <- function(f) {
    x <- readBin(f, "raw", file.info(f)$size)
    for (key in c("/CreationDate (", "/ModDate (")) {
      for (s in grepRaw(key, x, fixed = TRUE, all = TRUE)) {
        i <- s + nchar(key)
        while (x[i] != charToRaw(")")) { x[i] <- charToRaw("0"); i <- i + 1L }
      }
    }
    digest::digest(x, algo = "md5", serialize = FALSE)
  }
  ```

- **ED-g's and ED-h's cairo PDFs are not compared.** Cairo compresses their metadata, so their timestamps cannot be masked. E49's PDFs were not compared either.

## 6. File checksums from dry run 5

**Each file is identical in `outputs/figures/` and, with the `2026-11-03_` prefix, in `docs/figures/`.**

| file | md5 |
|---|---|
| `E50_bcl2_family_heatmap_level.pdf` (masked) | `8596db5a587ac6d6245e453f13176050` |
| `E50_bcl2_family_heatmap_level.svg` | `40d17cf762e8209ee35a1686a1aa6c93` |
| `E50_bcl2_family_heatmap_level.png` | `f76fcc12a93fec02630728b60a0d619f` |
| `E50_bcl2_family_heatmap_priority.pdf` (masked) | `376abe107e1c0fb0d0f6a4d8f54e7ae2` |
| `E50_bcl2_family_heatmap_priority.svg` | `801d3a85d2a6c9660b503c8512aa6843` |
| `E50_bcl2_family_heatmap_priority.png` | `155d63be36e02b0a3a51985d3714210c` |
| `E50_bcl2_family_heatmap_reference_score.pdf` (masked) | `1f7f92d53aebad96207691457457e17c` |
| `E50_bcl2_family_heatmap_reference_score.svg` | `885864d1be8c9f055061435f8bc4809e` |
| `E50_bcl2_family_heatmap_reference_score.png` | `4cb4ef63a842163584e3ed4c45f05f43` |
| `E50_signature_strips.png` | `def837722081339891abdd7cf8d6d63b` |
| `E50_raw_beside_adjusted.png` | `2d306de1bf99ed4d155154d1acb18112` |

**The tables, in `outputs/tables/`:**

| file | md5 |
|---|---|
| `E50_heatmap_cells.csv` | `a600ab2e7b07d20e17726ced8b53fcd8` |
| `E50_assertions.csv` | `d17f4fd17f5b6882856b378f08e181ff` |

## 7. The tracked figures

The author's run writes them to `docs/figures/`. **Commit them once the run is verified.**

| file stem (`2026-11-03_` prefix) | files | display |
|---|---|---|
| `E50_bcl2_family_heatmap_level` | `.pdf`, `.svg`, `.png` | Fig. 4g |
| `E50_bcl2_family_heatmap_priority` | `.pdf`, `.svg`, `.png` | ED-e |
| `E50_bcl2_family_heatmap_reference_score` | `.pdf`, `.svg`, `.png` | ED-f |
| `E50_signature_strips` | `.png`, `.pdf` | ED-g |
| `E50_raw_beside_adjusted` | `.png`, `.pdf` | ED-h |

**For Illustrator:**
- **The PDFs** reference Helvetica as a standard font, without embedding it. Illustrator supplies the Mac's Helvetica and keeps the text editable, as for Fig. 1g's panels.
- **The SVGs** carry named groups: cells, values, row labels, group names, column headers and the colour bar.
- **Both are at final size.**
