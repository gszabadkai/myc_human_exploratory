---
date: 2026-11-06
status: >
  PRE-RUN NOTE for E57, committed with the script BEFORE the author's run.
  AMENDED for the SECOND run. The first run matched this note's first version
  (dba2962) and is verified (6d988fa). The author then ruled that all text be
  black and the pCR panel's connectors thinner and dashed (7ed8316). Sections 6
  and 7 now give the second version's digests and checksums, from its dry
  runs 3, 4 and 5, which are identical.
relates-to:
  - docs/2026-11-06_E57_outcome_forest_panels_ruling.md (cb3c7d3, the specification)
  - docs/2026-11-06_E57_black_text_dashed_connector_ruling.md (7ed8316, the amendment)
  - docs/2026-11-06_e57_result.md (6d988fa, the first run, verified)
  - scripts/E57_outcome_forest_panels.R
  - docs/2026-10-22_fig4_human_panels_proposal.md (section 7.6: the digests of E45's object, which E57 checks)
  - scripts/E45_fig4_human_panels.R (the forests redrawn), scripts/E54_myc_coamplification_panel.R (Fig. 4x's grammar), scripts/E52_puma_bclxl_panel_26mm.R (the two slates)
---

# E57 pre-run note: Fig. 4's outcome forests, as two panels

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E57_outcome_forest_panels.R`. It takes a few seconds.
3. **It needs** `results/e45_fig4_human_panels.rds`. It checks five of its elements against the digests recorded before E45's run, and stops if any differs.

**It writes:**
- `results/e57_outcome_forest_panels.rds`;
- `outputs/tables/E57_drawn_values.csv` and `E57_assertions.csv`;
- each panel twice, as a PDF with live Helvetica and a PNG preview:
  - `outputs/figures/E57_pcr_forest.{pdf,png}` and `E57_metabric_survival_forest.{pdf,png}`;
  - `docs/figures/2026-11-06_E57_pcr_forest.{pdf,png}` and `2026-11-06_E57_metabric_survival_forest.{pdf,png}`, which are tracked.

## 1. What it draws

**E45's two forests, each as its own panel** (the ruling's message):
- **pathological complete response** in three neoadjuvant cohorts, with their pooled row;
- **breast-cancer-specific survival** in METABRIC.

**The changes the ruling asks for:**
- **Colours** (question 1):
  - OXPHOS in Fig. 4x's dark slate (#243642), open without the proliferation term and filled with it;
  - proliferation in mid slate (#8FA3B0).
- **All circles**, with proliferation told apart by colour.
- **No bold**: the PDF read-back checks the font of every string.
- **All text black** (the amendment, `7ed8316`). In the first version three kinds of text were grey, as in E45: the rows' second lines, I² and the cues' words. The PDF read-back checks the colour of every string.
- **The keys:** "OXPHOS", "OXPHOS, adjusted for proliferation" and "Proliferation, adjusted for OXPHOS" (question 4). METABRIC's panel has no model without proliferation, so its key has the last two.
- **Fig. 4x's size** (question 2):
  - a 26-mm axis;
  - 5-pt text with 6-pt titles;
  - its line weights and 1.05-mm circles;
  - the key above the data area.
- **The direction cues stay** (question 3).

**The rest is E45's:**
- the rows and their labels, and the column headers;
- the log axes, with E45's limits and breaks;
- the line from OXPHOS without to with proliferation, now 0.26 pt and dashed (the amendment);
- the rule before the fourth row;
- the contrast column.

**E45's row geometry is kept**, scaled to 5.6 mm per unit.

**The pages:**
- **pCR:** 63.5 × 54.7 mm (180 × 155 pt).
- **METABRIC:** 68.1 × 51.9 mm (193 × 147 pt).
- **Both data areas** are 26 × 33.0 mm. Aligned by their bottom edges, the two panels have their axes and rows in line; a check reads this from the PDFs.

**It fits nothing.**
- **Every value drawn is E45's**, checked by digest: 20 markers with their intervals, and 8 contrasts.
- **Every ratio and bound** is checked to be exp() of E45's stored log value, exactly.
- **Every contrast label** is E45's stored contrast, exponentiated.
- **No OXPHOS or proliferation ratio is printed**, as in E45.

## 2. The dry runs

**Before the first run, re-reading the script fixed three things:**
- a leftover duplicate line;
- line widths compared as exact floating-point numbers;
- a variable named `c`.

**Dry run 1 stopped in the layout,** after the checks of sections 1 and 2 had passed. Inside `tibble()`, the name `key` meant the column, not the key's table. The table is now `kt`.

**Dry run 2 stopped at a check: "no two strings come within 0.2 mm".**
- **The cause:** the raised "2" of I², set at the full 5 pt, reaches 2.12 mm above its baseline, which is exactly the line step. It touched "Pooled" above it.
- **The fix:** the pooled row's two lines are now 0.4 mm further apart, still centred on the row (section 4).

**Dry run 3 stopped at a check of the two panels against each other.**
- **The cause:** the check compared label baselines, which the pooled row's spacing now moves.
- **The fix:** it now compares the row, cue and rule heights. A new check reads the two axes' heights from the PDFs.

**Dry runs 4 and 5 stopped at the tick-label and I² checks.** "2" appears twice in the pCR panel, as a tick label and in I². Each check now looks on its own side of the axis.

**Dry run 6 passed every check of each panel,** but the two axes sat 0.03 pt apart.
- **The cause:** rounding each page up to whole points left 0.37 and 0.31 pt, spread above and below the panel.
- **The fix:** each page's top and right margins now take that fraction of a point, so each page is its panel's natural size.

**Dry run 7 passed all 75 checks,** and both panels were previewed. Dry runs 8 and 9 were identical to it.

**Then the check labels were changed to name the pCR panel "pCR"** (they had read "PCR"). That changed the `assertions` digest.

**Dry runs 10 and 11 are identical:**
- the seven saved elements' digests;
- each PNG in both folders;
- each PDF with its timestamps zeroed;
- both tables.

**Dry run 12 sourced the script in one R session, as Positron does:** 269 objects, none holding invalid UTF-8. It matched dry runs 10 and 11 in every digest and checksum.

**Also checked:**
- the PDFs, previewed with Apple's engine;
- the ortholog tripwire, which is clean.

### The second version (the amendment, `7ed8316`)

**What changed in the script:**
- every string is black. The cues' arrows are lines and stay grey30;
- the connectors are 0.26 pt and dashed;
- two new read-backs: the colour of every string the PDF shows, and the dash of every line and circle.

**Its dry run 1 stopped at the new check that every string is black.**
- **The cause:** the reader took each string's colour where its text block began, but R sets the colour inside the block, just before the string.
- **The fix:** it now reads the colour at the show operator (`Tj`, `TJ`).

**Its dry runs 2 to 5 passed all 80 checks.**
- **Runs 3 and 4 are identical** in every digest and checksum.
- **Run 5 sourced the script in one R session, as Positron does,** and matched them: 276 objects, none holding invalid UTF-8.
- **The PDFs were previewed again** with Apple's engine.

**The values drawn, the contrasts and the legend drafts are byte-identical to the first version's:** `frame`, `contrast` and `legend_draft` keep their digests.

## 3. What the checks hold (80)

**What is copied** (section 1 of the script):
- **From E54:** 42 objects, which are Fig. 4x's sizes, type, marks and helpers.
  - Each was evaluated from E54's own script and compared with E57's, value for value and function body for body.
- **From E45:** 12 values, which are its rows, row geometry, keys' order and axes. All its titles, headers, axis titles and cues are E45's, character for character, as are its greys.
- **From E52:** the two slates.
- **From the ruling:** the four answers and the keys' words, read from its file.
- **From the amendment:** its two requests, read from its file.

**E45's values** (section 2):
- five elements of E45's object, by digest, and its own 27 assertions;
- every ratio, bound and contrast, as in section 1;
- every interval inside E45's axis limits.

**The layout, in mm** (section 3):
- no two strings within 0.2 mm of each other;
- no string within 0.2 mm of a marker, an interval or a cue's arrow;
- the columns apart;
- the tick labels at least 0.5 mm apart;
- everything inside the page's margins.

**The PDF:**
- **The page:** its size in whole points.
- **The strings:**
  - exactly the strings expected, and no others;
  - all live Helvetica, none embedded and **none bold**;
  - **all black** (the amendment);
  - only the I of I² oblique;
  - 6 pt for the titles, 5 pt for the rest, none rotated.
- **The axes and lines:**
  - the 26-mm axis, and the line at 1 in its place;
  - every label, header, contrast, key line, title and cue at its place, to 0.02 pt;
  - each tick label centred under its tick;
  - the axis title's lines 2.4 mm apart.
- **The circles:** every one 1.05 mm across, at its ratio and row, in its key's fill and stroke.
- **The intervals:** each one's ends at its bounds, in its key's colour.
- **The rest:**
  - the rule, the connectors and the cues' arrows, each at its place, colour and width;
  - the connectors dashed, 1.02 pt on and off, and every other line and every circle solid (the amendment);
  - the key 1.2 mm above the data area, its lines 2.8 mm apart.

**The PNG:** its size at 600 dpi, and every circle's centre showing its fill exactly.

**Across the two panels:**
- the same row heights;
- the two axes at the same height above their pages' bottom edges.

**No check's label is built from an output path**, so the dry run's labels are the real run's.

## 4. The choices the ruling leaves open

The script's header gives each one with its reason. Each is the author's to reverse.
- **Each page is its panel's natural size.** The two pages differ (63.5 and 68.1 mm wide), and aligned by their bottom edges, their axes and rows are in line.
- **The title is kept, plain, at 6 pt, at the page's left margin above the key.** E45 centred it over the forest, and had the key below.
- **The tick labels drop trailing zeros**, as Fig. 4x's do: "0.5", not "0.50". At 26 mm, "1.50" and "2.00" would touch.
- **The line from OXPHOS without to with proliferation is dark slate, 0.26 pt and dashed** (the amendment).
  - **The pattern** is R's "dashed": dashes and gaps of four line widths, 1.02 pt.
  - **Why dark slate:** E45's 60% opacity would make a grey close to proliferation's mid slate.
- **The "2" of I² is full size and raised,** as E52 sets its exponents, and nothing is under 5 pt. So the pooled row's two lines are 0.4 mm further apart than the other rows'.
- **All text is black** (the amendment). E45's greys stay on lines only: the cues' arrows in grey30 and the rule in grey75. Making the arrows black is a one-line change.
- **The cues are scaled with the type:** arrowheads of 0.8 mm (E45: 1.1 mm) at 0.53 pt (E45: 0.64 pt), with each text 1 mm above its arrow.
  - **On the survival panel,** "longer survival" runs about 1.3 mm past the axis's left end, as question 3 said it would.
- **PDF and PNG only.** E45 wrote no SVG either.

## 5. The legend drafts (E45's, edited for the new keys)

These are E45's drafts with five edits. Each edit is applied only if its old text occurs exactly once.
- **The symbols:** "triangle" is replaced, and the keys are named.
- **The m0 and m1 labels** are dropped from the I² sentence and the contrast sentence.
- **The METABRIC model** loses its "(m1)".

**pCR:**

> Odds ratios for pathological complete response (pCR) per within-cohort standard deviation of the OXPHOS subunit score (GSVA) and of the proliferation score, in three neoadjuvant chemotherapy cohorts (GSE194040, n = 988, 319 pCR; GSE164458, n = 482, 236 pCR; GSE25066, n = 465, 93 pCR). Open circles, OXPHOS, from a logistic model with receptor subtype and treatment arm where each varies within the cohort; dark circles, OXPHOS, adjusted for proliferation, from the same model with the proliferation score added; light circles, proliferation, adjusted for OXPHOS, from that model. Pooled, inverse-variance meta-analysis, random effects (DerSimonian-Laird) where I2 >= 50% and fixed effect otherwise; I2 is given for the three, in that order (49%, fixed effect; 0%, fixed effect; 67%, random effects). Right, the ratio of the proliferation to the OXPHOS odds ratio within each adjusted model, its variance taken from that model's covariance matrix; the pooled value is the fixed effect meta-analysis of the three within-model ratios (I2 0%). Bars, 95% confidence intervals.

**METABRIC:**

> Hazard ratios for breast-cancer-specific death per within-cohort standard deviation of the OXPHOS subunit score (GSVA) and of the proliferation score in METABRIC, from Cox models with PAM50 subtype, treatment group, OXPHOS and proliferation. Rows: ER-positive by immunohistochemistry (n = 1,505, 431 events); ER-positive by expression (ER.Expr) (n = 1,512, 432 events); all patients in one model (n = 1,981, 623 events); ER-negative by immunohistochemistry (n = 435, 183 events). Dark circles, OXPHOS, adjusted for proliferation; light circles, proliferation, adjusted for OXPHOS. Right, the ratio of the proliferation to the OXPHOS hazard ratio within each model, its variance taken from that model's covariance matrix. The endpoint, treatment and model differ from those in the pCR panel, and the axes are not shared. Bars, 95% confidence intervals.

**What travels with these panels is unchanged** (CLAUDE.md, the E39 to E44 line):
- a significant contrast makes neither coefficient a finding;
- nothing here supports a distant-outcome claim for respiration;
- METABRIC does not replicate the pCR analysis;
- both panels use one instrument only, GSVA.

## 6. The digests the author's run must reproduce

`digest::digest()` of each saved element, from the second version's dry runs 3, 4 and 5. `run_info` is excluded.
- **The first version's digests** are in this note's first version (`dba2962`).
- **`frame`, `contrast` and `legend_draft` are unchanged** from the first version.

| element | digest |
|---|---|
| `frame` | `75f37d50c274fd4b070a52bd475202c6` |
| `contrast` | `4b472baf850f2b56ba9beeb73e25cc8c` |
| `layout` | `c4965986e2973151110d9f07f752dba5` |
| `legend_draft` | `4a6d2851ab3cb826a45bfa807bd8d9e6` |
| `assertions` | `ee8e9aeb0b76341e784ccb04d70652ab` |
| `figure_files` | `f6673f04e3d2b9ab3580683f2225a941` |
| `spec` | `c87d18a5110f4b241d3d970bdb4c3904` |

## 7. The checksums

| file | md5 |
|---|---|
| `E57_pcr_forest.png`, in both folders | `38f120d6b7e0fcc4934e4cf368f6e951` |
| `E57_pcr_forest.pdf`, in both folders, timestamps zeroed | `447570ae413b393bee60710aa2e10ad0` |
| `E57_metabric_survival_forest.png`, in both folders | `49722c4ab764f418024c0d5757e17545` |
| `E57_metabric_survival_forest.pdf`, in both folders, timestamps zeroed | `727cc61b93024a9861ff2b42d77f1ac2` |
| `E57_drawn_values.csv` | `46140675ee0b127ac5529c5851a734d3` |
| `E57_assertions.csv` | `c07b2a2dbc4ccb1229a0110d70e09101` |

## 8. After the author's second run

1. **Compare** the seven digests (section 6) and the checksums (section 7).
2. **Commit** the four tracked figures, and add the run to the result note.
3. **Record** that E57's panels replace E45's figure for the paper's last two Fig. 4 panels.
4. **What still follows, for the text ruling:**
   - the panels' letters;
   - the legends (section 5);
   - the sentences that cite them.
