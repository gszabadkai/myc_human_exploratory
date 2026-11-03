---
date: 2026-11-03
status: >
  E49 WRITTEN and NOT YET RUN BY THE AUTHOR. These digests are recorded and
  committed BEFORE that run. No result is reported here; the result note
  follows the author's verified run (stage 3).
relates-to:
  - docs/2026-11-03_E49_signature_declaration.md (7bd7865, committed with data/annotations/)
  - docs/2026-11-03_signature_and_fig4g_brief.md (eb45a9b)
  - scripts/E49_signature_and_fig4g.R
---

# E49: digests recorded before the author's run

## 1. What the dry runs established

**Every output of every run went to the session scratchpad.** Nothing was written in the repo.

| run | what happened |
|---|---|
| three early attempts | stopped at assertions before completing. See 2.1 |
| the first complete run | **192 of 192 assertions passed, with no warnings.** Its panels were looked at, and seven layout problems were fixed (2.2) |
| dry runs 1 and 2, final version | **24 of 24 saved elements agree**; **192 of 192 assertions pass** in each; **no warnings**; **33 of 33 PNGs and tables byte-identical**: 12 PNGs in each figure directory and 9 CSVs |

**The elements excluded from comparison:**
- `built`;
- `fits_effectors`, 96 `lm` objects, which are not digest-stable (R rule 5). `fit_digests` is their stable handle.

### 1.1 The gates, all as the declaration requires

**Part A, the re-reads:**

| gate | what it checks | result |
|---|---|---|
| structure (7.1) | the 144 signature fits | every `vcov` is square and named by its terms; `PROLIF` is in the 48 `m3` fits and in none of the others |
| G-A1 | every signature fit read, at `m1` and `m3` | coefficients and `vcov` equal E47's by digest; SEs and CIs recompute exactly |
| G-A2 | every signature and effector difference | the independent recomputation equals the `L V L'` contrast |
| G-A3 | the twelve genes' values from E48 | equal E47's saved fits exactly |

**Part B, the effector fits:**

| gate | what it checks | result |
|---|---|---|
| G-B1 | the gene-scoring code | reproduces all twelve genes in E47's frame, in both cohorts, below 1e-10. The effectors are scored by that function |
| G-B2 | the fitting code, refitting `BBC3` in all 32 combinations | **reproduces E47's coefficients and `vcov` by digest, exactly** |
| G-B3 | the effectors against every fitted predictor set | in none of them |
| G-B4 | the void cells | only HALLMARK's `m3` cells are void, as in E47 |

**The checks against earlier numbers:**

| gate | what it checks | result |
|---|---|---|
| G-C1 | the raw S6 correlations | equal the reconciliation proposal's to 2 dp |
| G-C2 | the Fig. 1i-style slope | equals E47's `m3` priority OXPHOS coefficient for `BBC3` minus `BCL2L1`, to 1e-10. C3's interval is negated exactly |
| G-C3 | the placement rule | evaluated from E47. **It drew the branch declaration 2.4 already recorded: the MYC column stays in 4g** |

**Membership and colour:**
- **membership:** S6, S5 and S7 equal E47's `SIGS`, and no effector is in any of them;
- **the grey:** the circular cells' grey sits at least dE 19.1 from every colour on the scale.

**Most of the 192 assertions are display checks.** They confirm:
- every drawn cell and point equals its stored value, read back from `ggplot_build()`. For a heatmap that means each tile's colour is its value on the scale;
- the two versions of 4g, ED-e, ED-f and ED-h draw identical values, and version (b) has a dot exactly where the CI excludes zero;
- the circular cells are grey and off the scale;
- binding marks appear only on the four verified BH3-only rows, and the blank rows are drawn blank;
- the signature-sign column is S6's, and § appears on `BAD` and `BID` only;
- every limit is symmetric, contains every value and is rounded outward to 0.1;
- no text is below 5.5 pt, and no label carries beta, m1, m3, M_a or Delta;
- every panel is 178 mm wide;
- no legend draft uses a banned word, none calls a MYC score a signature, and the § note is the author's verbatim.

## 2. Changes before the final dry runs

### 2.1 Stops before a run completed

1. **Text below 5.5 pt.** `theme_void` and `theme_classic` leave unused elements at 0.8 of the base size, which is 5.2 pt. Every theme now sets every text element explicitly.
2. **Finding a layer.** ggplot2 4 does not keep a layer's data identical to the frame passed to it. Each layer's data now carries an explicit `.layer` tag, and the checks find layers by tag.
3. **ED-h's width.** Its two parts summed to 170.4 mm. The right part now takes the rest of the 178 mm, and the assertion requires an exact sum.

### 2.2 Legibility, after the first complete run

**The panels were looked at for legibility, clipping, palette and text size only.**
- **Viewing them necessarily showed the plotted values**, including quantities new to E49:
  - S5's and S7's proliferation coefficients;
  - the signatures' differences;
  - the effectors' coefficients, as colours;
  - the raw correlations;
  - the Fig. 1i-style decile means and slopes.
- **None is reported here, and nothing in content or layout was changed because of a value.**
- **The console logs were read only for assertion, warning and error lines.**

1. **4g, ED-e and ED-f: the left-hand headers are spread out.** "binds MCL1" and "signature sign" almost touched. The heatmap now starts 2.5 mm further right.
2. **4g, ED-e and ED-f: the empty band below the strip is smaller.** It was about 15 mm deep above the colour bar.
3. **ED-g: the content sits 2 mm lower, and the panel is 4 mm taller.** The "OXPHOS level" title was clipped at the top edge.
4. **ED-h: the rotated column headers have room.**
   - They were clipped at the top and collided with the cohort labels.
   - The significance dots moved into the cell corner, where they no longer touch the digits.
   - The empty band above the legends is smaller.
5. **ED-i: the direction cues sit 0.3 row-units lower.** They were tight under the last row's markers.
6. **ED-k: the slope annotation is now the subtitle.**
   - In one panel the fitted line crossed it.
   - The annotation now sits off the plotting area, so its place does not depend on the data.

## 3. Digests, recorded before the author's run

| element | digest |
|---|---|
| `signature_coefficients` | `4abd71bf1215901636a8b08adcaf5daa` |
| `genes12` | `4bdf30d2cb47a5b53b4320d1eabef304` |
| `effector_coefficients` | `ba4b1e6c4975ed773292f8520e534d5b` |
| `all_readouts_m3` | `7023883b873b17020e800838db591684` |
| `effector_coefficient_matrix` | `63189ae2502320a307dd92fa6f2bbd21` |
| `effector_vcov` | `d280a60734c51272c6b015be87b3bf88` |
| `effector_vif` | `eb58b00dd2def842bb1ec08cc7eb44ff` |
| `fit_digests` | `4af192a7d279f65922ad072f18d8c511` |
| `gate_a1` | `c13ab634d5ab4e03a1aac2710d753651` |
| `gate_a3` | `cda2d55edfba6565f7adce3c86c29aee` |
| `gate_b1` | `0f912f8ac07b376d8e82dbd2b3ef1b4a` |
| `gate_b2` | `4cc18961e033e313ba4e518da9ff4c84` |
| `placement_rule` | `82d62a1e160e2ff889e0aa854872690f` |
| `raw_correlations` | `f9f32a13550b07075b1690cb7d0a6a36` |
| `gate_c1` | `863f9bd8009cc560381e598636a7169b` |
| `expression_percentiles` | `ba77cc9d243eb80f9d9cfd026bd09a1f` |
| `fig1i` | `f9a71873b59f3d12a9dad7641c4c0eda` |
| `binding_code` | `fe47804b8722e2b0f61f5f90a7a4c64c` |
| `frame` | `0c86a29dbd4e6f6c97ab4ca80e3468b3` |
| `displays` | `419f5c39819c1f951129a794948b6294` |
| `legend_draft` | `6b7639c40dcde9ae894183f76b367829` |
| `assertions` | `2710612444bdbe6592173817b73c2d86` |
| `figure_files` | `26bc126367a079eef5d4fe0a0f779104` |
| `spec` | `ec618f2d9510e6dd1bb67248ca6afff3` |

## 4. Verifying the author's run

```r
x  <- readRDS("results/e49_signature_and_fig4g.rds")
el <- setdiff(names(x), c("built", "fits_effectors"))
sapply(el, function(k) digest::digest(x[[k]]))   # compare with section 3
```

**The PNGs and the nine tables are compared byte for byte** against the md5s in section 5:
- the PNGs in `outputs/figures/` and `docs/figures/`;
- the tables, `outputs/tables/E49_*.csv`.

The PDFs carry their creation time and are not compared.

## 5. File checksums from dry run 2

**The PNGs.** Each is identical in `outputs/figures/` and, with the `2026-11-03_` prefix, in `docs/figures/`.

| file | md5 |
|---|---|
| `E49_fig4g_heatmap_unmarked.png` | `c36c710030bf53930478218c66c9e793` |
| `E49_fig4g_heatmap_marked.png` | `2c1cb69999d89bb322ff5a6ea308ff68` |
| `E49_heatmap_priority_unmarked.png` | `af9dd4b4d1914feaf4d3556939af8a4c` |
| `E49_heatmap_priority_marked.png` | `d164dc7ded1fd3586542ee20714ac750` |
| `E49_heatmap_reference_score_unmarked.png` | `8ca4857d077a24699781ecd691dbdd26` |
| `E49_heatmap_reference_score_marked.png` | `eb9078abb0fe7db4c7c3ec13357035cc` |
| `E49_signature_strips.png` | `b975db756d0c72229c2e96efb37c13c8` |
| `E49_raw_beside_adjusted_unmarked.png` | `a1292b6a9828f31ee32cc0a804f66a4d` |
| `E49_raw_beside_adjusted_marked.png` | `e77de1ba874b784988b57ad049053bb4` |
| `E49_fifteen_genes_intervals.png` | `7a9758b3a50afd87712d0b564b58205a` |
| `E49_myc_ladder_signatures_effectors.png` | `cde12c948b8483229258d34a4dbeab09` |
| `E49_puma_bclxl_against_priority.png` | `2ae69de0d62007c0cee7504e3b376420` |

**The tables:**

| file | md5 |
|---|---|
| `E49_signature_coefficients.csv` | `710b19b1136d634b356e69b162adb890` |
| `E49_effector_coefficients.csv` | `2e8bc94d6320b4bed967e8ef1f8dda61` |
| `E49_all_readouts_m3.csv` | `bd89ab7e810a973811f4c1458a724f71` |
| `E49_raw_correlations.csv` | `6f2142ce3919c710f84f4e64b2118238` |
| `E49_expression_percentiles.csv` | `56f9870a362bdb190cfd6edf2b739af2` |
| `E49_placement_rule.csv` | `b7816d68317b6137c9b855f4760ea3f2` |
| `E49_fig1i_style.csv` | `ed97865b67af62e61b3789a4399a2160` |
| `E49_fig1i_style_deciles.csv` | `a694cf62bdd4970182538b730d58ee79` |
| `E49_assertions.csv` | `d5fa5c59edb3aaa5a239308a9964ffbb` |

## 6. The tracked figures

The author's run writes them to `docs/figures/`, each as `.png` and `.pdf`. **Commit them once the run is verified.**

| file stem (`2026-11-03_` prefix) | display |
|---|---|
| `E49_fig4g_heatmap_unmarked` | Fig. 4g, version (a) |
| `E49_fig4g_heatmap_marked` | Fig. 4g, version (b) |
| `E49_heatmap_priority_unmarked`, `_marked` | ED-e |
| `E49_heatmap_reference_score_unmarked`, `_marked` | ED-f |
| `E49_signature_strips` | ED-g |
| `E49_raw_beside_adjusted_unmarked`, `_marked` | ED-h |
| `E49_fifteen_genes_intervals` | ED-i |
| `E49_myc_ladder_signatures_effectors` | ED-j |
| `E49_puma_bclxl_against_priority` | ED-k |

## 7. Run time and memory

- **A dry run takes about one minute.**
- **Part B loads each cohort's linear matrix once**: 18,115 x 1,095 for TCGA, and 18,153 x 3,207 for SCAN-B. It frees each one after use.
