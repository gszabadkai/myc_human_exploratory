---
date: 2026-11-06
status: >
  PRE-RUN NOTE for E56, committed with the script BEFORE the author's run.
  Six dry runs, with every output redirected to a scratch folder; runs 4, 5
  and 6 are identical. It gives the digests and checksums the author's run
  must reproduce.
relates-to:
  - docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md (dd5004a, the specification)
  - docs/2026-11-06_E56_effector_quadrant_means_declaration.md (922f2ad, the declared fit)
  - scripts/E56_quadrant_configuration_fifteen_genes.R
  - scripts/E34_quadrant_configuration.R and scripts/E46_fig4_quadrant_heatmap.R (the model and the panel)
  - scripts/E50_fig4g_in_fig1g_style.R (Fig. 4g's grammar)
---

# E56 pre-run note: Extended Data Fig. 6d, the quadrants with fifteen genes

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E56_quadrant_configuration_fifteen_genes.R`. It takes about a minute, because it reads the two linear expression snapshots.
3. **It needs:**
   - `results/e34_quadrant_configuration.rds` and `results/e46_fig4_quadrant_panels.rds`, which it checks against their records;
   - the frozen snapshot in `data/from_validation/`.

**It writes:**
- `results/e56_quadrant_configuration_fifteen_genes.rds`;
- `outputs/tables/E56_effector_quadrant_means.csv`, `E56_heatmap_values.csv` and `E56_assertions.csv`;
- the panel twice, as a PDF with live Helvetica and a PNG preview:
  - `outputs/figures/E56_quadrant_configuration_fifteen_genes.{pdf,png}`;
  - `docs/figures/2026-11-06_E56_quadrant_configuration_fifteen_genes.{pdf,png}`, which is tracked.

**On the declaration's date:** its file is dated 2026-11-06, but it was committed at 23:55 GMT on 2026-11-05 (`922f2ad`). Git records the true time. It was committed before this script was written and before any effector was fitted.

## 1. What it does, in order

1. **Checks E34's and E46's saved objects,** each by the checksum its own record uses.
   - E34's result note records `serialize(version = 3)` md5s, not `digest::digest()`. By that method its object matches.
   - E46's matches the final version in its record.
2. **Gate G1:** the expression path reproduces E34's twelve genes in its saved frames. **Maximum deviation 0.**
3. **Gate G2:** E34's model path, copied verbatim, reproduces E34's stored m1 group means for the twelve genes: 768 values. **Maximum deviation 0.**
4. **Gate G3:** the same path reproduces E46's signature strip, m1 and m3: 128 values. **Maximum deviation 0.** E46's heatmap cells are also checked to be E34's stored means.
5. **Fits the declared new quantity:** BAX, BAK1 and BOK on E34's m1, 48 group means, on the same tumours as E34's fits (n = 938 and 3,143).
   - **Range of the means across the 16 cells:** BAX −0.58 to +0.48, BAK1 −0.40 to +0.33, BOK −0.23 to +0.20. None exceeds ±0.6.
   - **What is saved (R rules 4 and 5):**
     - the coefficient vectors and variance matrices, and a digest table;
     - the refitted `lm` objects, shown to have identical coefficients, saved but not digested;
     - the analysis frame with the effectors' columns added.
6. **Draws Fig. 6d.**

## 2. What it draws

- **Rows:** Fig. 4g's fifteen genes in its four groups and order, with the group names on the right, as in 4g.
- **Columns:** both OXPHOS splits, GSVA then mitoPPS, at full width (the ruling, question 3). Each cohort's four quadrants are in E46's order: OXPHOS-low (MYC-low, MYC-high), then OXPHOS-high.
- **Cells:** the m1 group means, in Fig. 1g's palette at ±0.6, colour only, as E46 drew them. Twelve genes come from E34 via E46; the three effectors come from this script's declared fit.
- **The colour key reads "z-score"** (the ruling). Its colour bar is 4g's, at the top right.
- **Below:** E46's signature strip, the means with 95% confidence intervals.
  - Filled is m1; open is m3.
  - The key reads "Signature" and "proliferation adjusted" (question 4).
  - The circles are 1.05 mm across, as in Fig. 4x.
- **The page is 179.9 × 93.2 mm.** The columns are 8.51 mm wide and the rows 4g's height. The strip's axis runs from −0.8 to +0.8, with E46's ticks.

## 3. The dry runs

**Dry run 1 stopped before running.** The file-writing tool had turned two Unicode escapes in comments into literal characters. They are now written as code-point names, and the script is ASCII.

**Dry run 2 stopped after the gates and the fit, in the layout.** E46's strip carries its own `x` column, which collided with the layout's on a join. The strip's `x` is now dropped.

**Dry run 3 stopped in the layout:** a local `Y_TITLE` (the colour bar's) hid the strip's axis title. The latter is now `STRIP_TITLE`. Then all 25 checks passed. The split titles were moved 0.6 mm further from the cohort headers.

**Dry runs 4 and 5 are identical:**
- the eleven saved elements' digests;
- the PNG in both folders;
- the PDF with its timestamps zeroed;
- all three tables.

**Dry run 6 sourced the script in one R session, as Positron does:** 228 objects, none holding invalid UTF-8. It matched dry runs 4 and 5 in all 18 digests and checksums.

**Also checked:**
- the 54 lines copied from E50 and the 40 from E34 are verbatim;
- the ortholog tripwire is clean.

**The PDF was previewed with Apple's engine.**

## 4. What the checks read back from the PDF and the PNG

- **The cells:** 240, each in its size, colour and place. The colour bar is a gradient.
- **The labels:** every one is live Helvetica (no image, no embedded or Type 3 font), at its size and angle. The gene symbols are oblique, and nothing else is.
- **The strip:** its 32 intervals have their lengths, and its 34 circles (32 means and the key's 2) are 1.05 mm across.
- **The page** is the panel's size, within 180 mm.
- **The PNG** is the panel at 600 dpi, and every cell shows its colour exactly.

## 5. The choices the ruling leaves open

The script's header and comments give each one with its reason. Each is the author's to reverse.
- **The group names are on the right,** as 4g has them. E46 had two groups, rotated, on the left.
- **No signature-sign column.** 4g draws none; the author adds annotations in Illustrator.
- **No values printed in the cells,** as E46 drew them.
- **E46's split titles are kept:** "OXPHOS split by GSVA" and "OXPHOS split by mitoPPS". The rest of ED Fig. 6 says "OXPHOS level" and "OXPHOS priority".
- **The strip keeps E46's design,** at Fig. 4x's circle size.

## 6. The legend draft (Claude Code's, for the author)

> Fifteen BCL-2 family genes across the four MYC x OXPHOS groups of tumours, in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143). Left, the OXPHOS score by GSVA; right, by mitoPPS; MYC activity is the same in both. Columns: within each cohort, OXPHOS-low then OXPHOS-high tumours, MYC-low then MYC-high within each. Cells, z-scores: adjusted group means of Blom-transformed expression in within-cohort standard deviations, from a linear model with group, PAM50 subtype, tumour purity and leukocyte fraction in TCGA-BRCA, and group and PAM50 subtype in SCAN-B. The three effectors (BAX, BAK1, BOK) were fitted for this figure on the same model, post hoc. Rows and colours as in Fig. 4g. Strip, the signature: filled, the same model; open, proliferation adjusted; bars, 95% confidence intervals.

## 7. The digests the author's run must reproduce

`digest::digest()` of each saved element, from dry runs 4, 5 and 6. `run_info` and the `lm` objects (`fits`) are excluded: a fitted model is not digest-stable (R rule 5). Its coefficients and variance are digested in `coefficients` and `fit_digests`.

| element | digest |
|---|---|
| `effector_means` | `a337ae7e265f9e7f42c4fad5de3032d1` |
| `gates` | `8239d3eca27394c243f7a583dfc5fa09` |
| `values` | `a8639b6ba2599f57e074b48ce4ba76b9` |
| `frame` | `80a3e45e21d6aecdd6b5333552dca01b` |
| `coefficients` | `ba0d1413fac72548e607cca41bc8235a` |
| `fit_digests` | `ea2a8422515a590912bef1a99a38e6e3` |
| `layout` | `4969bf0289240b0525f4664b6a51163c` |
| `legend_draft` | `86e8d2bc1145bed84aa298d964e7cde6` |
| `assertions` | `9bb26c561040a4fbcc2669ebd9876048` |
| `figure_files` | `79e3800063a9974d42ef673c76da6cc7` |
| `spec` | `184f9e7ce0130e43a98defa725fc5f93` |

## 8. The checksums

| file | md5 |
|---|---|
| the PNG, in both folders | `7483bb2c8f7977d71ba33f420b817543` |
| the PDF, in both folders, timestamps zeroed | `5bd9dca83838f8e5a1a167788f846be6` |
| `E56_effector_quadrant_means.csv` | `f34a937b557549cc9f125a46378d43b6` |
| `E56_heatmap_values.csv` | `1230940ce55f2a8432f78a0f61226785` |
| `E56_assertions.csv` | `8b54faf74a33476d049dfad974095dc7` |

## 9. After the author's run

1. **Compare** the eleven digests (section 7) and the checksums (section 8).
2. **Commit** the tracked figure pair with a result note.
   - The note records the effectors' means as drawn, not read.
   - It says whether the declaration's section 6 expectation held: BAX's and BAK1's OXPHOS-high quadrants above their OXPHOS-low ones.
