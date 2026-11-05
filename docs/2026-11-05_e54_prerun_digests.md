---
date: 2026-11-05
status: >
  PRE-RUN NOTE for E54, committed with the script BEFORE the author's run.
  AMENDED for the SECOND run: the author's point-size ruling (6550011) made
  the circles smaller after the first run, which matched this note's first
  version (d290bc3) and is verified (ec9b590). Sections 8 and 9 now give the
  second version's digests and checksums, from dry runs 9, 10 and 11, which
  are identical.
relates-to:
  - docs/2026-11-05_E54_coamplification_panel_ruling.md (8fde82a, the specification; it amends decision 5)
  - docs/2026-11-05_E54_point_size_ruling.md (6550011, the amendment: circles 1.05 mm across)
  - docs/2026-11-05_e54_result.md (ec9b590, the first run, verified)
  - data/from_validation_g2/README.md (c1dceed, the snapshot E54 draws)
  - scripts/E54_myc_coamplification_panel.R
  - myc_human_validation d3ac60e, docs/2026-08-28_G2_result.md (the committed note E54 checks against; read only)
  - docs/2026-11-05_E52_prism_size_ruling.md (9c2ca44, the 26-mm panel this one is sized like)
---

# E54 pre-run note: Fig. 4x, co-amplification with MYC

## 0. Running it

1. **Restart R.**
2. **Source** `scripts/E54_myc_coamplification_panel.R`. It takes a few seconds.
3. **It needs only** the tracked snapshot `data/from_validation_g2/g2_cnv_cooccurrence.rds`. It checks the file's md5 and stops if it differs.

**It writes:**
- `results/e54_myc_coamplification_panel.rds`;
- `outputs/tables/E54_odds_ratios.csv` and `E54_assertions.csv`;
- the panel twice, as a PDF with live Helvetica and a PNG preview:
  - `outputs/figures/E54_myc_coamplification_odds_ratios.{pdf,png}`;
  - `docs/figures/2026-11-05_E54_myc_coamplification_odds_ratios.{pdf,png}`, which is tracked.

## 1. What it draws

**Twelve adjusted odds ratios from the validation study's G2**, in TCGA-BRCA. Each compares MYC-amplified tumours (GISTIC +2) with the others, for BCL2L1 gain (+1 or more) and MCL1 amplification (+2).

| gene | rows | covariate | odds ratio [95% CI] | n | co-altered |
|---|---|---|---|---|---|
| BCL2L1 gain | all tumours | aneuploidy score | 1.77 [1.26, 2.47] | 1,043 | 132 |
| | all tumours | fraction genome altered | 1.36 [0.98, 1.88] | 1,043 | 132 |
| | Basal | aneuploidy score | 2.69 [1.33, 5.62] | 169 | 52 |
| | Her2 | aneuploidy score | 0.99 [0.31, 3.20] | 78 | 9 |
| | LumA | aneuploidy score | 1.70 [0.88, 3.26] | 490 | 26 |
| | LumB | aneuploidy score | 0.90 [0.44, 1.86] | 197 | 34 |
| MCL1 amplification | all tumours | aneuploidy score | 1.69 [1.17, 2.44] | 1,043 | 53 |
| | all tumours | fraction genome altered | 1.60 [1.10, 2.31] | 1,043 | 53 |
| | Basal | aneuploidy score | 3.12 [1.53, 6.56] | 169 | 33 |
| | Her2 | aneuploidy score | 1.22 [0.24, 4.85] | 78 | 3 |
| | LumA | aneuploidy score | 0.68 [0.25, 1.56] | 490 | 6 |
| | LumB | aneuploidy score | 1.12 [0.46, 2.58] | 197 | 10 |

**What the rows are in G2:**
- **the rows for all tumours** are its pre-registered primary tests, under both of its pre-specified covariates;
- **the subtype rows** were declared there as secondary and exploratory.

**The page is 113 × 131 pt (39.9 × 46.2 mm),** with a data area of 26 × 30.8 mm. E52's panel is 41.3 × 45.2 mm.

## 2. Where its numbers come from, and how they are checked

- **The source is the snapshot of G2's saved result** (`c1dceed`), with its md5 checked at the start. Its README records its provenance.
- **Every number drawn that the committed G2 note prints is checked against the note,** at the note's precision: 24 numbers. They are the twelve odds ratios, and the intervals of the four rows for all tumours and the two Basal rows. The note gives no intervals for Her2, LumA or LumB, so those come from the snapshot alone.
- **The snapshot itself was checked against the note on 2026-11-05:** 46 of 46 numbers match (its README).

## 3. The dry runs

**Dry run 1 passed all 26 checks, but its preview showed a fault:** the hyphen in "co-occurrence" printed as a minus sign. That is base `pdf()`'s encoding, as in E51 to E53. The PDF now writes a soft hyphen, and the PNG a hyphen.

**Dry runs 2 to 5 passed.** Dry runs 3 and 4 were identical, and dry run 5, sourced in one session, left no invalid UTF-8 in the workspace.

**Then one check was strengthened.** Its label said the gene names were the only italic text, but it tested only that an oblique font was present.
- The PDF read-back now records each string's font.
- The check now tests that exactly "BCL2L1" and "MCL1" are in Helvetica's oblique, and every other word in Helvetica.
- The drawing did not change: the PNG and the PDF are identical to dry runs 3 to 5. Only the check table changed.

**Dry runs 6 and 7 are identical:** the seven saved elements' digests, the PNG in both folders, the PDF with its timestamps zeroed, and both tables.

**Dry run 8 sourced the script in one R session, as Positron does,** then scanned the workspace: 194 objects, none holding invalid UTF-8. It matched dry runs 6 and 7 in all 13 digests and checksums.

**The PDF was previewed with Apple's engine** (`qlmanage -t`).

**Then the first run, and the amendment.**
- **The author's first run matched this note's first version** in every element and file (`ec9b590`).
- **Then the author ruled that the circles be smaller** (`6550011`): 1.05 mm across, outline included, in place of 1.6 mm. Point size 1.5 with a 0.45 outline became 1.0 with a 0.3 outline. The key's symbols are now placed by the circle's true width.
- **Two checks were added:**
  - the size is the ruling's, as its own file records;
  - every circle measured in the PDF, path and outline, is 1.05 mm across.

  Tried on the first version's PDF, the reader finds its 14 circles at 1.58 mm, so the check can tell the sizes apart.
- **Dry runs 9 and 10 are identical, with all 28 checks passing.** The circles measure 1.051 to 1.055 mm.
- **Dry run 11 sourced the script in one session:** 198 objects, none holding invalid UTF-8. It matched dry runs 9 and 10 in all 13 digests and checksums.
- **Unchanged from the first version:** every number drawn, the frame, the note checks, the legend draft and the odds-ratio table.

## 4. What the checks read back from the PDF

- the PDF shows exactly its 21 words: the headers, the row labels, the ticks, the axis title and the key;
- the gene headers and the axis title are 6 pt, everything else 5 pt, and nothing is rotated;
- *BCL2L1* and *MCL1* are in Helvetica's oblique, and every other word in Helvetica;
- the x axis measures 26 mm, and the line at 1 spans the data area at its middle;
- every interval's two ends sit where its odds ratios fall on the 26-mm log2 axis, to 0.02 pt;
- the axis title's two lines are 2.4 mm apart;
- the key lies 1.2 mm above the data area, its rows 2.8 mm apart;
- BCL2L1's rows come before MCL1's;
- the 14 circles, 12 marks and the key's 2, are 1.05 mm across, outline included.

## 5. The choices the ruling leaves open

The script's header and comments give each one with its reason. Each is the author's to reverse.
- **The axis is log2, from 0.125 to 8,** so that 1 sits in the middle. It is labelled at 0.25, 1 and 4, with unlabelled ticks at each doubling between.
- **No values are printed beside the rows.** The ruling's preview listed them at the right. Printing them would widen the panel by about 10 mm.
- **The Normal stratum is left out** (n = 35), as the preview said.
- **The all-tumours result takes two rows,** aneuploidy score (filled) above fraction genome altered (open), with one "All tumours" label between them.
- **A two-line key above the data area** names the two covariates, as E52's key names the cohorts: "adjusted for aneuploidy score" and "adjusted for fraction genome altered".
- **Each gene's header sits in its own empty row,** at 6 pt with the gene name in italics. It runs into the data area's left side and ends 6.9 mm before the line at 1.
- **The axis title is Claude Code's wording:** "Odds ratio for co-occurrence" over "with MYC amplification".
- **The marks are TCGA-BRCA's dark slate** (`#243642`), the colour E52's key gives TCGA-BRCA. The line at 1 is E48's grey.
- **The circles are 1.05 mm across** (the point-size ruling). Intervals and axes follow E48 and E52.

## 6. What the text ruling must weigh (the ruling, question 4)

**The panel shows what the author's sentence describes, and some things the sentence does not say.**
- **BCL2L1 in all tumours:** 1.77 [1.26, 2.47] adjusted for aneuploidy score, but 1.36 [0.98, 1.88] for fraction genome altered. The second interval reaches 1, and the panel shows it.
- **BCL2L1 by subtype:**
  - only Basal's interval excludes 1;
  - Her2 is 0.99 and LumB 0.90;
  - "across subtypes" rests on a homogeneity test that did not reject (Breslow-Day p = 0.33), which the panel does not show.
- **MCL1 by subtype:**
  - only Basal's interval excludes 1;
  - Breslow-Day p is 0.051 across PAM50 subtypes, and 0.022 by TNBC status.
- **MCL1 at any gain inverts** (OR 0.50 [0.35, 0.72]). The panel shows amplification (+2) only. The 2026-11-02 brief (section 9) asked that the sentence state the inversion or drop the MCL1 clause.
- **The subtype rows are secondary and exploratory** in G2. The legend draft says so.

## 7. The legend draft (Claude Code's, for the text ruling)

> BCL2L1 gain and MCL1 amplification in MYC-amplified tumours of TCGA-BRCA (n = 1,043): odds ratios against tumours without MYC amplification, from logistic models adjusted for aneuploidy score (filled) or, in all tumours, for fraction genome altered (open), with 95% confidence intervals, on a log scale. Below each gene, PAM50 subtypes, adjusted for aneuploidy score (Basal, n = 169; Her2, n = 78; LumA, n = 490; LumB, n = 197). Copy number from GISTIC calls: MYC and MCL1 amplification, +2; BCL2L1 gain, +1 or more. The tests in all tumours were pre-registered in the validation study; the subtype estimates were declared there as secondary and exploratory.

## 8. The digests the author's run must reproduce

**For the second run.** `digest::digest()` of each saved element, from dry runs 9, 10 and 11. `run_info` is excluded.

| element | digest | the same as the first version's |
|---|---|---|
| `frame` | `5eca0135d88aab31ab1357a02235669a` | yes |
| `note_checks` | `927d07b3eea27df51d152d7a6a607070` | yes |
| `layout` | `9d2ae89d96aac3a3448461ff6ff66f33` | no: the key's symbols moved |
| `legend_draft` | `11d3e89301706c8aca6526b759687567` | yes |
| `assertions` | `8ac01a9bcfdd42d0cef23b38d4c4f685` | no: 28 checks, not 26 |
| `figure_files` | `b1200d88e7686a1a75669560ee6e38d2` | no |
| `spec` | `91cab01f93747d94991778c018edd0ea` | no: the point size |

**No check's label is built from an output path**, so the dry run's labels are the real run's (handoff lesson 3).

## 9. The checksums

| file | md5 |
|---|---|
| the PNG, in both folders | `ddbd6fb1d3436714e881a9b954a43ec9` |
| the PDF, in both folders, timestamps zeroed | `862b9ef86ad97a5d2876c7b3c9a21623` |
| `E54_odds_ratios.csv` | `e7ce3260de106ae539387f1190ee0c1c`, unchanged |
| `E54_assertions.csv` | `5d4ee7c78516b6339b4c1b72e30a6822` |

**The PDF's masked checksum** is E51's: the md5 of the file with NUL bytes removed and the digits after `/CreationDate (D:` and `/ModDate (D:` replaced by a single 0.

## 10. After the author's second run

1. **Compare** the seven digests (section 8) and the checksums (section 9).
2. **Commit** the tracked figure pair, and add the second run to the result note.
3. **Then the text ruling:**
   - the sentence, against section 6;
   - the legend;
   - the panel's letter.
