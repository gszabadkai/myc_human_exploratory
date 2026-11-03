---
date: 2026-11-03
status: >
  E48 WRITTEN and NOT YET RUN BY THE AUTHOR. These digests are recorded and
  committed BEFORE that run. No result is reported here; the result note
  follows the author's verified run (stage 3).
relates-to:
  - docs/2026-11-02_E48_separation_declaration.md (2941438, committed alone)
  - docs/2026-11-02_separation_line_brief.md (cbf2e4c)
  - scripts/E48_separation_and_fig4g.R
---

# E48: digests recorded before the author's run

## 1. What the dry runs established

**There were three dry runs, every output redirected to the session
scratchpad.** Nothing was written in the repo.

| run | version | result |
|---|---|---|
| 1 | first | 175 of 175 assertions passed, with two warnings and several legibility problems in the panels (section 2) |
| 2 and 3 | final | **33 of 33 saved elements agree**; **175 of 175 assertions pass**; **no warnings**; **26 of 26 PNGs and tables byte-identical** (8 PNGs in each figure directory, 10 CSVs) |

Run 1's two warnings came from the G-A2 digest check: a numeric vector was
coerced by `all()`. The check itself passed, and it is now a typed table.

**The elements excluded from comparison:**
- `built`;
- `fits_c1` (6 `lm` objects) and `fits_c4` (12), which are not digest-stable
  (R rule 5). `fit_digests` is their stable handle.

### 1.1 The gates, all as the declaration requires

| gate | what it compared | result |
|---|---|---|
| structure (6.1) | all 544 `m1` and `m3` fits | square `vcov`s named by their terms; `PROLIF` in all 272 `m3` fits and none of the `m1` fits |
| G-A1 | 1,440 re-read OXPHOS and MYC values against `e47$results` | maximum deviation 0 |
| G-A2 | 240 fits against E47's `fit_digests` | every coefficient vector and `vcov` matches; SEs equal `sqrt(diag(vcov))` to 2.8e-17; CIs recomputed with `qt(0.975, df)` to 0 |
| G-A3 | | the only void cells read are HALLMARK's `m3` cells |
| G-A4 | 240 omitted-variable identities | hold to 1.6e-16. **`delta` is positive in 16 of 16 cells**, as section 11 assumes |
| the difference | recomputed independently from the three variance elements | equals the contrast `L V L'` to 6.9e-18 |
| no signature reads | | **Part A reads gene fits only. No S6, S5, S7 or arm fit is read or drawn** |
| Part B | | confirmed unavailable: E40 saves no fitted model and no `vcov`; E43's eight fits are `m1` only, with no MYC term |
| ED-a | E39's three gate values | equal the declaration's to 4 dp |
| C1 | | the join is 938 of 938 on the patient barcode; PAM50 agrees in 938 of 938; TP53 status is 325 mutant and 613 wild-type, none missing. Without the TP53 term, the six fits reproduce E47's coefficients and `vcov`s exactly |
| C4 | | with purity and leukocyte fraction restored, the 12 fits reproduce E47's exactly |
| C3 | | the copied code at `m1` reproduces E47's 4 x 2,000 bootstrap replicates and its four section-9 rows exactly. The `m3` design matrix reproduces E47's stored `m3` OXPHOS coefficients |
| `g3` | the 48 pairs recomputed on E47's frame | reproduce E47's `g3` exactly |
| layout 3 | | the clinical values equal E45's plotting frame exactly, and the contrasts are E45's labels unchanged |
| ED-b | | the relabelled Display C values are identical to E47's display frames |
| legend drafts | | none uses a word brief 4.9 keeps out |

**The 128 display assertions check:**
- every drawn point and interval equals its stored value, read back from
  `ggplot_build()`, including after the log transform;
- open is without the proliferation term and filled is with it;
- no void cell is drawn as read;
- every data colour comes from the validated palette;
- every axis is symmetric, contains every interval and is rounded outward
  to 0.1;
- no text is below 5.5 pt (the smallest is 5.50 pt);
- no label carries beta, m1, m3, M_a or Delta;
- the width is 178 mm.

**Palette (12.1).** The dataviz checker passed `#009E73, #D55E00, #0072B2`
on all pairs, light surface, with no warning. `#0072B2` was the first
candidate under the declared rule.

## 2. Legibility changes after run 1

The rendered panels were looked at for legibility, clipping, palette and text
size only. **Viewing them necessarily showed the plotted values. Nothing in
content or layout was changed because of a value, and no result is reported
here.**

1. **The two-line column headers now share one weight, bold.** plotmath
   ignores `fontface`, so line 1 had rendered plain beside a bold line 2.
2. **The legends of layout 2 and ED-b wrap to two rows.** At 178 mm they were
   clipped at both ends.
3. **Layout 3's facet titles sit on two lines.** Two were clipped by their
   column.
4. **The direction cues sit 0.45 row-units lower in every panel.** They
   touched the last row's markers. The figure heights were raised to match.
5. **ED-c's facets are spaced 5 mm apart.** Their tick labels had collided.
6. **The difference panels of ED-d and S-1 carry no direction cue.** A
   difference of two coefficients has no expression direction, so the cue was
   wrong, and in ED-d it was also clipped. The axis title "Proliferation minus
   OXPHOS" carries the sign. ED-d's narrow difference panels also have three
   ticks and a two-line axis title.
7. **ED-d's three read genes are now bold italic.** `bold(italic())` had
   rendered them plain.

## 3. Digests, recorded before the author's run

| element | digest |
|---|---|
| `labels` | `dbd7ab8fbd45f6dbb335cf7c3abbc7b8` |
| `primary_cells` | `a098f6b1845421bee42e1af07f4b5909` |
| `agreement` | `cd09069f09067a2e1a88d69cae21db24` |
| `text_row` | `9cb853ba20cc6e762749e895ab17a8d9` |
| `nine_gene_tally` | `2f467c4537b73e98af4482209910b9d6` |
| `coefficients` | `7229ecddb0bc3b04d36b0fb91e5d8e91` |
| `nesting` | `53a793f55138c565852f13f8e50c44f4` |
| `gate_a2` | `162eb0078aa56e2e6af119f59ecfc4f1` |
| `part_b` | `671d7ed1394c2c217c422926073c3a99` |
| `c1` | `22d9e600a83bbd84a88b74095d444426` |
| `c1_coefficient_matrix` | `ae294da83b8099479d79ba7014ba456b` |
| `c1_vcov` | `0c9d33c15adff822aec7852a75ea543d` |
| `c1_vif` | `afc23b07e486862be61ad8415909e6c3` |
| `c1_gate` | `8da712e7a8159dbc2589293e7ec34321` |
| `part_c2` | `8c711da19168b9d8ca385781a7f01e64` |
| `c3` | `c079b13c18eee9a09875dc7e6571612e` |
| `c3_holds_m3` | `bb73ad91bcb7e948250d465016f7bea8` |
| `c3_boot_replicates_m3` | `abe91bd185702784c73aa0a754e202df` |
| `c4` | `57864a4d12aaf58d6b5c4ca3bd69af8f` |
| `c4_coefficient_matrix` | `2b537a068e1564036d14b0603d77c618` |
| `c4_vcov` | `7a4106fe16a11f4bfe44d15aed30f82a` |
| `c4_vif` | `3a0de0a1ff5cd615466d619b5019d05a` |
| `c4_gate` | `ea0a1517e5386fa64766394cd0e58446` |
| `sensitivities` | `51523cf0c65a27867ba2bcbb7b31e58d` |
| `fit_digests` | `22025b3b83dfcb66619d9da7635bc378` |
| `correlations` | `360159b45ec5bac05c0e95d48d354622` |
| `e39_gate` | `50143900aef967cbfcad7008963a26ee` |
| `frame` | `ab23f60dd60e1f829944351f8040ae90` |
| `displays` | `c8397b9106fbe770476a7b10bb513f14` |
| `legend_draft` | `b9dca928c086b687785da7436aca297e` |
| `assertions` | `f6a6ee8b7d1edcbfad1a277df3b6ef56` |
| `figure_files` | `6f40e9b9454c254ed2e36dfe79858b64` |
| `spec` | `9039744d4243578c0d074af120f40c61` |

## 4. Verifying the author's run

```r
x  <- readRDS("results/e48_separation_and_fig4g.rds")
el <- setdiff(names(x), c("built", "fits_c1", "fits_c4"))
sapply(el, function(k) digest::digest(x[[k]]))   # compare with section 3
```

**The PNGs in `outputs/figures/` and `docs/figures/`, and the ten tables in
`outputs/tables/E48_*.csv`, are compared byte for byte with the dry run.**
PDFs carry their creation time and are not compared.

## 5. The tracked figures

These are written by the author's run into `docs/figures/`, each as `.png` and
`.pdf`. **Commit them after the run is verified.**

| file stem (`2026-11-03_` prefix) | display |
|---|---|
| `E48_separation_layout1_two_scores` | Fig. 4g, layout 1 |
| `E48_separation_layout2_three_scores` | Fig. 4g, layout 2 |
| `E48_separation_layout3_with_outcomes` | Fig. 4g, layout 3 |
| `E48_score_correlations` | ED-a |
| `E48_bcl2l1_against_bbc3` | ED-b |
| `E48_myc_score_ladder` | ED-c |
| `E48_twelve_genes_three_coefficients` | ED-d |
| `E48_sensitivities_tp53_composition` | S-1 |

## 6. The result note's structure, fixed by the author on 2026-11-03, before the run

**Stage 3 writes `docs/<date>_e48_result.md` in four parts:**

1. **The declared reading, briefly:**
   - the three labels;
   - the agreement check;
   - C1, C3 and C4.
2. **Exploration.** All twelve genes, three coefficients each. Any pattern,
   expected or not, is described as a hypothesis.
3. **Against the mouse.** By direction only, never with coefficients side by
   side, and including disagreements:
   - Fig. 1a,c: proliferation maintained whilst the apoptotic response was
     lost;
   - Fig. 1f: OXPHOS against proliferation in the developing gland;
   - Fig. 2f,g: Pgc-1alpha raised PUMA independently of Myc;
   - Fig. 3h,i: high OXPHOS sustained only with Bcl-xL;
   - the orthotopic tumours, where Bbc3 fell as OXPHOS rose, the opposite of
     the human direction (provenance note 6.3);
   - C1 against the mouse's p53 observations: Foxo3 in the gland, and the
     unchanged p53 targets in the orthotopic arms;
   - the death step, which the human data cannot test.
4. **The text.** The simplest version the numbers allow, drafted from brief
   section 8, for the review pack. **Only this and the figure go into the
   paper.**

**The review pack (brief 6.5)** has one page per layout:
- the panel at final size;
- the finalised text with that layout's citations;
- the descriptive legend draft;
- STYLE metrics per paragraph.

**STYLE.md, which the author will add before stage 3, governs its prose.**
