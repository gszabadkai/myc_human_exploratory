# CLAUDE.md — MYC / OXPHOS / cell death in human breast cancer (exploratory)

Rules and context for Claude Code sessions in this repo. Read this first, then
`docs/2026-08-31_aim.md`.

## What this repo is

An **openly exploratory** study of how MYC activity and OXPHOS relate in human
breast tumours, and how the cell-death programme sits on that plane. Two
cohorts, TCGA-BRCA (n=1,095) and SCAN-B / GSE202203 (n=3,207). Starting point:
Menegollo, Bentham et al., *Cancer Res* 2024 (CAN-23-3172).

## THIS IS EXPLORATORY. THAT IS THE POINT, AND IT IS ALSO THE RISK.

The most important rule, and the mirror image of the sibling repo's.

- **Nothing here is pre-registered.** Every finding is hypothesis-generating and
  must be labelled as such, in notes, in figures and in conversation.
- **Multiple comparisons are the default state, not an exception.** The atlas is
  21 MYC estimators x 4 instruments x ~30 arms x 8 strata. Any single cell of it
  is uninteresting. Report structure, gradients and reproducibility across
  cohorts — never a p-value plucked from the grid.
- **"Significant" is not a result here.** Consistency across cohorts,
  instruments and estimators is. If something holds in TCGA and SCAN-B, on GSVA
  and mitoPPS, across low- and high-entanglement MYC signatures — that is worth
  something. One cell of the grid is not.
- **When something looks real, say what would falsify it** — and write that down
  *before* the next analysis, so the exploratory phase can hand a real
  hypothesis to a confirmatory one.

### The sibling repo, and how to read it

`myc_human_validation` (`/Users/gs/code/myc_human_validation`, frozen at
`d3ac60e`) is a **completed, pre-registered** study of one narrow hypothesis. It
found nothing supported. It is not reopened, and it does not go in the paper
**except for G2 co-amplification: one MAIN panel (Fig. 4x) and one main-text
sentence** (the author's decision 5, 2026-11-02, amended on 2026-11-05 from an
Extended Data panel: `docs/2026-11-05_E54_coamplification_panel_ruling.md`,
`8fde82a`). **G2's saved result is snapshotted** in
`data/from_validation_g2/` (`c1dceed`), and all 46 numbers its committed note
prints match the snapshot. **G2's pooled tests were pre-registered; its subtype
strata and Breslow-Day tests were declared secondary and exploratory.**

- Its `CLAUDE.md` forbids post-hoc hypotheses. That rule is correct *there* and
  **actively wrong here**. Do not import it.
- Its *results* are trustworthy in a way nothing here is, because they were
  declared first. **Do not import a conclusion from it without checking whether
  it was pre-registered.** Its dated notes say which.
- Read it read-only:
  `git -C /Users/gs/code/myc_human_validation show d3ac60e:<path>`.
  Never write to it. It is not attached to sessions here.

## THIS IS A HUMAN REPO

- Human MitoCarta 3.0, human gene symbols, human-native gene sets only.
- The mouse repo (`/Users/gs/code/myc_mouse`) is **not** attached.
  Read it read-only via `git -C ... show <ref>:<path>` if ever needed. It moved
  there from `/Users/gs/G/data/MK_myc_2022/myc_mouse`, which is now renamed
  `myc_mouse_OLD_do_not_use`; **never point anything at an `_OLD_` copy.**
  Dated notes before 2026-09-28 still give the old path, as written.
- **No ortholog function is called anywhere in this repo**, in either direction.
  The check is for *calls*, not for the word — comments asserting the rule are
  the reason it must be narrowed to a `(`, or the tripwire always fires:

  ```
  grep -rnE "(mouse_to_human|human_to_mouse|ortholog[s]?)[[:space:]]*\(" scripts/
  ```

  must return nothing. It catches `mouse_to_human(`, `ortholog (` and
  `convert_orthologs(`; it ignores prose. The cell-death and MYC sets are
  human-native — see their READMEs for why that is established rather than
  assumed.

## Current phase

**NEWEST, 2026-11-06, small hours: `E57`, Fig. 4's last two main panels.
Its first run is VERIFIED** (7 of 7 digests, 75 of 75 assertions;
`docs/2026-11-06_e57_result.md`, `6d988fa`).

**The author then amended it**
(`docs/2026-11-06_E57_black_text_dashed_connector_ruling.md`, `7ed8316`):
- **no coloured text:** every string is black;
- **the pCR panel's connectors** are thinner (0.26 pt) and dashed.

**The second version is RUN and VERIFIED** (7 of 7 digests, 80 of 80
assertions, every file identical to the dry run; the result note, section
3). **Its panels are Fig. 4's last two**, in place of E45's figure:
`docs/figures/2026-11-06_E57_{pcr_forest,metabric_survival_forest}.{pdf,png}`.
The specification is `docs/2026-11-06_E57_outcome_forest_panels_ruling.md`
(`cb3c7d3`).

**Open for the author:**
- does "no coloured text" also reach the white values on the darkest cells?
  They are in Fig. 4g, with ED-e and ED-f (E50), and in ED Fig. 6a (E53).
  Neither is changed;
- the cues' arrows on E57's panels are lines, and stay grey30. Should they
  be black?
- **What it draws:** E45's two outcome forests, each as its own panel:
  - pCR in three neoadjuvant cohorts;
  - METABRIC breast-cancer-specific survival.
- **At Fig. 4x's size:** a 26-mm axis, 5-pt text, 1.05-mm circles, and the
  key above the data area.
- **The colours:**
  - OXPHOS in dark slate, open without the proliferation term and filled
    with it;
  - proliferation in mid slate.
- **The marks and text:**
  - all circles;
  - no bold, and all text black (the amendment);
  - the pCR connectors dashed at 0.26 pt (the amendment).
- **The keys:** "OXPHOS", "OXPHOS, adjusted for proliferation" and
  "Proliferation, adjusted for OXPHOS". The direction cues stay.
- **It fits nothing.** Its values are E45's, checked by digest.
- **Its 80 checks** include, in the PDF and the PNG:
  - every string's font (none bold) and colour (all black);
  - every line's dash;
  - every circle's colour.
- **What travels with the forests is unchanged:**
  - a significant contrast makes neither coefficient a finding;
  - nothing here supports a distant-outcome claim for respiration;
  - METABRIC does not replicate pCR;
  - GSVA only.

**Earlier, 2026-11-06, just after midnight: `E55` and `E56` are RUN and
VERIFIED; their panels are Extended Data Fig. 6c and 6d** (the author's
ruling, `docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md`, `dd5004a`). The
author's reply: "done, both are OK." The text ruling follows; the
handoff's section 0 lists what it must settle.
- **`E55` (Fig. 6c): 7 of 7 digests, 22 of 22 assertions, every file
  identical to the dry run** (`docs/2026-11-06_e55_result.md`). It draws
  E49's unadjusted Spearman correlations (all seven measures, both
  cohorts) for Fig. 4g's fifteen genes and the signature, in 4g's grammar,
  heatmap size and colours, colour only. It fits nothing.
  - **R13's citation waits on the author.** R13 cites ED-h for a
    joint-model claim. 6c keeps only ED-h's unadjusted side, and Fig. 4g
    draws the coefficients (result note, section 4).
- **`E56` (Fig. 6d): 11 of 11 digests, 25 of 25 assertions, every file
  identical to the dry run** (`docs/2026-11-06_e56_result.md`). It redraws
  E46's quadrant panel with 4g's fifteen genes, groups and colours, both
  splits at full width (179.9 mm); key "z-score", strip key "Signature"
  and "proliferation adjusted". It settles the 2026-11-02 brief's open
  item 1, E46 in Extended Data.
  - **Its one new quantity was declared first**
    (`docs/2026-11-06_E56_effector_quadrant_means_declaration.md`,
    `922f2ad`): the effectors' (BAX, BAK1, BOK) m1 quadrant means on E34's
    model, post hoc, **drawn, not read**. Three gates reproduced E34's
    twelve genes, E34's stored means and E46's strip with deviation 0
    before the fit.
  - **The declaration's expectation held, and licenses nothing:** BAX's and
    BAK1's OXPHOS-high quadrants sit above their OXPHOS-low ones in all four
    blocks. It came from E49 on the same tumours, so it is a consistency
    check, not a replication.
- **The open choices of both pre-run notes are read as accepted** by the
  author's reply, which came after they were listed for confirmation. Each
  is still the author's to reverse.

**Earlier, 2026-11-05, night: `E54` is RUN and VERIFIED, twice** (first
version 7 of 7 digests, 26 of 26 assertions; then, after the author made the
circles smaller, `6550011`, the second version 7 of 7 and 28 of 28;
`docs/2026-11-05_e54_result.md`). **The second version is Fig. 4x.** Its
tracked figures are
`docs/figures/2026-11-05_E54_myc_coamplification_odds_ratios.{pdf,png}`.
- **Decision 5 is amended** (the author's ruling,
  `docs/2026-11-05_E54_coamplification_panel_ruling.md`, `8fde82a`): the
  co-amplification result (validation study G2) becomes a MAIN panel,
  Fig. 4x, with a 26-mm axis like E52's.
- **What it draws:** adjusted odds ratios for BCL2L1 gain and MCL1
  amplification with MYC amplification in TCGA-BRCA (n = 1,043):
  - all tumours, under both pre-specified covariates (aneuploidy score,
    filled; fraction genome altered, open);
  - then Basal, Her2, LumA and LumB, declared secondary and exploratory in
    G2.
- **Its numbers** come from G2's saved result, snapshotted with its
  provenance in `data/from_validation_g2/` (`c1dceed`, md5-checked). All 46
  numbers the committed G2 note prints match the snapshot, and E54 checks
  the 24 it draws.
- **The author's sentence waits on a text ruling** (the ruling's question 4).
  The pre-run note's section 6 lists what it must weigh:
  - BCL2L1 crosses 1 under fraction genome altered;
  - only Basal's interval excludes 1 for either gene;
  - "across subtypes" rests on a homogeneity test that did not reject;
  - MCL1 inverts at any gain.
- **Choices the ruling left open** are in the pre-run note, section 5. One is
  that no values are printed beside the rows, though the preview showed them.

**Earlier on 2026-11-05, late evening: `E53` is RUN and VERIFIED** (6 of 6
digests, 43 of 43 assertions, every file identical to the dry run;
`docs/2026-11-05_e53_result.md`). **Its panels replace E48's ED-a and ED-b
for the paper, as Extended Data Fig. 6a and 6b.** Its digests were committed
before the run (`docs/2026-11-05_e53_prerun_digests.md`, `f3745f0`).
- **What it draws:** E48's ED-a and ED-b, which the author numbers Extended
  Data Fig. 6a and 6b, at half width (at most 88 mm), by the author's ruling
  (`docs/2026-11-05_E53_half_width_ed_panels_ruling.md`, `7325ff1`):
  - **6a:** the two triangles side by side with 4.6-mm cells; the
    neoadjuvant values in a row of tiles under the colour key;
  - **6b:** no direction arrows; rows "TCGA: OXPHOS level" and so on; the
    key "BCL2L1", "BBC3", "Difference", each beside "proliferation
    adjusted"; x axes from -0.1.
- **It fits nothing.** Its data are E48's saved display frames, checked by
  digest. Its legend drafts are E48's, with 6a's "Right" now "Below".
- **Three choices from the pre-run note's section 4 are still to confirm
  with the author:** -0.1 as an unlabelled tick; no "open/filled" lines in
  the key; plain titles.

**Earlier on 2026-11-05, evening: `E52` is RUN and VERIFIED** (8 of 8
digests, 38 of 38 assertions, every file identical to the dry run;
`docs/2026-11-05_e52_result.md`). **Its panel replaces E51's for the
paper.** The tracked figures are
`docs/figures/2026-11-05_E52_puma_bclxl_against_priority_26mm.{pdf,png}`.
Its digests were committed before the run
(`docs/2026-11-05_e52_prerun_digests.md`, `bd0245b`).
- **What it draws:** E51's panel at the size of Fig. 4's Prism graphs, by
  the author's ruling (`docs/2026-11-05_E52_prism_size_ruling.md`,
  `9c2ca44`):
  - a 26-mm square data area, on a page of the panel's natural size;
  - axis titles without "adjusted", which the legend carries;
  - the two *P* labels, with the cohort names, in a key above the data area;
  - the axes and ticks still at 0.53 pt.
- **It fits nothing.** Its slopes, points, lines, *P* values and legend
  draft are E51's, checked by digest.
- **The Illustrator edit of E51's titles is no longer needed.** The text
  still follows, for a later ruling: R6's citation, the panel's letter and
  the legend.

**Earlier on 2026-11-05: `E51` is RUN and VERIFIED** (8 of 8 digests,
26 of 26 assertions, every file identical to the dry run;
`docs/2026-11-05_e51_result.md`). ED-k becomes a main panel (the author's
ruling, `docs/2026-11-05_E51_main_panel_ruling.md`, `86f38aa`). Its text
follows, for a later ruling: R6's citation, the panel's letter and the
legend.
- **What it draws:** PUMA relative to BCL-xL against OXPHOS priority, with
  TCGA-BRCA and SCAN-B overlaid in one panel, in Fig. 1i's style. Each
  cohort's *P* is printed inside the panel, and there is no title.
- **Its one new quantity is that *P*:** the two-sided t-test of the drawn
  slope, from E49's saved residuals, with no refit. TCGA-BRCA *P* = 0.017,
  SCAN-B *P* = 7.8 x 10^-9.
- **Kept as declared:** ED-k's model and the stored slopes. ED-b stays in
  Extended Data.
- **Its digests were committed before the run**
  (`docs/2026-11-05_e51_prerun_digests.md`, `1b0be30`), and the run matched
  them.
- **The crash after each run is fixed** (`4924f29`, confirmed by the
  author). The cause was the PDF's raw bytes left in the workspace: they are
  not valid UTF-8, and Positron's Variables pane crashed on them. **In any
  script, read a binary file inside a function and return only ASCII**
  (handoff, lesson 5).

**Earlier on 2026-11-05: the author's v6 narrative is committed, and its
changes wait on the author's ruling.**
- **v6** is `docs/Fig. 4 narrative v6 human arm and clinical link.md`
  (`c3b3482`). It is the author's draft, and it proposes six changes to the
  ruled text: five sharpening changes and a re-cut of the four Methods
  sentences outside STYLE.md's range. Until the author rules, the ruled text
  in `pack_text.py` stands. A ruling is recorded verbatim, alone, then applied
  (handoff section 0).
- **A narrative-shift explainer for co-authors** is published privately
  ("Fig. 4 Narrative Shift"). Its source is in
  `docs/2026-11-05_narrative_shift_src/` (`9d43357`). It reads only verified
  saved objects.

**2026-11-04: the author's STOP 3 rulings, their addendum and the
binding ruling are RECORDED and APPLIED. Nothing is pending from Claude
Code.** The entry
point is **`docs/2026-11-03_e50_handoff.md`**. Behind it:
1. **`docs/2026-11-03_stop3_rulings_consolidated.md`** (`cb9b7ce`), the
   rulings, verbatim. They replace every chat version of items 1 to 6 and
   amendments A to F. **`docs/2026-11-03_stop3_rulings_addendum.md`**
   (`3f9ffbe`), verbatim, amends them;
2. **`docs/2026-11-03_e50_review_pack.pdf`**, 4g at final size, with the
   rulings' text, the edited 4g legend and the nine Extended Data panels;
3. **`docs/2026-11-03_e49_result.md`, Part 4**, the text for the paper. It
   supersedes the text committed there at `4fd661b`.

**What the rulings decide:**
- **E50's displays are accepted:**
  - Fig. 1g's source palette;
  - the shared ±0.6 scale;
  - the 1.6 mm gap;
  - white ink on 9 cells;
  - ED-g's ±0.9 scale (addendum 4).
- **THE REMAINING DISPLAY CHANGES ARE THE AUTHOR'S, IN ILLUSTRATOR, OUTSIDE
  THE REPO'S SCRIPTS. No new display script.**
  - the group names and the colour bar's title at 6 pt;
  - 4g's annotations: the binding code, the signature signs, the § marks and
    any protein names;
  - ED-g's "Proliferation minus OXPHOS" column, removed. It stays in E49's
    tables.
- **Extended Data in the paper: ED-a, ED-b and ED-e to ED-k**, each cited.
  - ED-c and S-1 stay in the repo, out of the paper.
  - ED-d stays superseded by ED-i.
  - The handoff's section 1 maps each panel to its file and its citation.
- **H2 IS ACCEPTED.**
  - **The signature's associations with OXPHOS and with proliferation are what
    its selection intended,** never independent evidence.
  - **The independent evidence:**
    - `BBC3` and `BCL2L1`;
    - MYC activity;
    - the effectors.
  - **The H2 test is NOT declared.** It is the answer if a referee asks.
- **The text** is the rulings' section 3, as the addendum amends it. **The 4g
  legend** is E50's draft with the rulings' section 4 edits.
- **The addendum answers the four observations left after `da548bd`:**
  - R2 gives the search's cohorts, n = 1,095 and 3,207. The model sets stay in
    the 4g legend;
  - the brief's rule 5.4 is reworded. `BAD` is explained "in the sentence that
    sets its rise against BCL-xL's, or the next";
  - N4 to N6 are dropped. **The Methods definition is E47 4a.7, item 5,
    verbatim, with two changes:**
    - the construction sentence now includes proliferation;
    - the variants sentence cites ED-g;
  - ED-g keeps its ±0.9 scale.

**How it was applied** (`da548bd`, then `8d7cd31`):
- **The text has one source, `docs/2026-11-03_e49_review_pack_src/pack_text.py`.**
  It checks every sentence and legend edit against the rulings, the addendum
  and E47's declaration: 32 comparisons, all exact. **No sentence is Claude
  Code's.**
- **Four sentences of the Methods definition are outside STYLE.md's 12 to 38
  words,** for the author to revise (handoff section 0):
  - D3 (58) and D10 (43), to split;
  - D4 (8) and D12 (7), under the hard floor.
- **THE BINDING CLASSES ARE ENTERED** (the author's binding ruling,
  `docs/2026-11-04_binding_classes_ruling.md`, `66c205a`; applied in
  `adda767`).
  - `BID` binds both. `BIK` and `BMF` bind BCL-xL and are coded as not
    binding MCL1. Their binding is weak, not absent, and the table notes it.
  - The sources are Ku et al. 2010 and Wang et al. 2023. The 4g legend and N2
    now cite them beside Chen 2005.
  - Legend edit 6 applies. E49, which read the table with these rows blank, is
    not re-run.

**Earlier the same evening: `E50` was RUN and VERIFIED.**

**What E50 is:** a display-only script, `scripts/E50_fig4g_in_fig1g_style.R`.
- **Its specification** is the author's further amendment to the STOP 3
  rulings, verbatim in **`docs/2026-11-03_E50_display_amendment.md`**
  (`6f5bf4d`).
  - The amendment replaces three rulings: the version choice, the grey cells
    and the E50 specification.
  - **So `E50` is not the test of H2.** If that test is declared, it takes
    `E51` or the next free number.
- **It verified** against digests committed before the run (`eae8ef5`):
  - 116 of 116 assertions;
  - 11 of 12 saved elements by digest. The twelfth, the assertions table,
    differs only in one label that names its output folder;
  - every comparable figure file matches its checksum.
- **What it draws, from E49's numbers and nothing else:**
  - **no grey anywhere, and no dots;**
  - **4g, ED-e and ED-f in Fig. 1g's style,** with no annotation columns, as
    a PDF with live Helvetica, an SVG and a PNG, at 87.6 x 58.2 mm;
  - **ED-g and ED-h re-coloured.**
- **The binding code, the signature signs, the § marks and any protein names
  are the author's, added in Illustrator.**
- **Its palette is Fig. 1g's source** (`myc_mouse` `2ebf646`,
  `_style.R` `div_cols`), interpolated in CIELAB as Fig. 1g's scale is. It is
  NOT the values the author measured from the PDF.
- **The choices the amendment left open** are in the pre-run note
  (`docs/2026-11-03_e50_prerun_digests.md`), section 3. Each is the author's
  to reverse.
- **Two device findings** (result note, section 4):
  - **the default quartz PNG device darkens colours,** so E48's and E49's PNG
    previews are darkened. Their PDFs are not;
  - **a `cairo_pdf` between two cairo PNGs moves ggplot's multi-line text.**

**E50 then waited on the STOP 3 rulings.** They arrived consolidated and are
applied (see the 2026-11-04 entry above).

**Earlier the same day, in the afternoon: `E49` is RUN and VERIFIED, and the
session stopped at STOP 3 for the author.** Where this block conflicts with
E50's amendment or the STOP 3 rulings, they win. The conflicts are the grey
cells, the two versions, the drawn binding-code columns, the palette, the
text and H2. The entry point was
**`docs/2026-11-03_e49_handoff.md`**. Behind it, in this order:
1. **`docs/2026-11-03_e49_result.md`**, the result note.
   - Part 4 is the text for the paper.
   - Its 4.1 lists every change from the brief's drafts, and each "author's
     call".
2. **`docs/2026-11-03_e49_review_pack.pdf`**, 4g in two versions, with the text,
   the legend draft and the Extended Data.
3. **`docs/2026-11-03_signature_and_fig4g_brief.md`**, the author's brief. Its
   amendment is recorded in the E49 declaration's section 3, not in the brief.

**`E49` is EXPLORATORY and POST HOC, with no reading rule and no label.** It
verified at:
- 24 of 24 digests, committed before the run at `a63ef6c`;
- 192 of 192 assertions;
- 33 of 33 PNGs and tables byte-identical to the dry run.

Part A re-read E47's saved fits. Part B fitted only the three effectors
(96 fits), behind a gate that reproduced E47's `BBC3` fits by digest.

**What now stands:**
- **THE SIGNATURE (S6) IS RESTORED as Fig. 4g's summary** (brief decision 1),
  with its membership fixed: `BAD`, `BIK`, `BBC3`, `BID` and `BCL2L1` enter
  +, `MCL1` enters -.
  - **Its OXPHOS coefficient is circular.** It is drawn grey and never read,
    and so is its difference (declaration 8.2).
  - **Its MYC activity coefficient on 4g's model is positive in all four
    blocks**, +0.21 to +0.42. The regulon agrees in sign, so the placement
    rule KEPT the MYC column.
  - **Its proliferation coefficient is negative in all 24 cells**, -0.22 to
    -0.38. It is now DECLARED (E49), but read it with the caveat at the end
    of this list.
  - Unadjusted, the signature barely correlates with proliferation (-0.09 /
    -0.04).
- **THE NEW 4g is a heatmap in Fig. 1g's layout.**
  - **Rows:** fifteen genes in Fig. 1g's four groups (effectors, activators,
    sensitisers, anti-apoptotic).
  - **Left:** the BH3 binding code (`data/annotations/bh3_binding_code.csv`;
    `BID`, `BIK` and `BMF` are blank until the author supplies them) and the
    signature signs. § marks `BAD` and `BID`.
  - **Columns:** three exposure blocks x two cohorts: OXPHOS level, MYC
    activity on the least-entangled score, proliferation.
  - **The model:** one `m3` model per row, on the least-entangled score.
  - **Beneath:** the S6 strip.
  - **Two versions, the author's choice at STOP 3:** (a) unmarked; (b) with a
    dot where the CI excludes zero.
  - **In Extended Data:** OXPHOS priority (ED-e) and the reference model
    (ED-f, beside E48's ED-c).
- **The effectors `BAX`, `BAK1` and `BOK` are outside the signature, and never
  added to it.**
  - They rise with OXPHOS and MYC activity in both cohorts.
  - Proliferation goes with lower `BAX` and `BOK`.
  - **`BAK1`'s proliferation sign splits by cohort.**
  - **No sentence reads an effector's transcript as its activation.**
- **THE TWO RETRACTIONS (brief 4.8) are used nowhere, in any document:**
  1. "Its [BAD's] rise increases reliance on BCL-xL rather than overriding
     it";
  2. the BAD-based argument for BCL-xL over MCL1, that BCL-xL can absorb
     the rising load and MCL1 cannot.
- **The text reads PUMA against BCL-xL.**
  - **BAD's rise is answered by its phosphorylation** (Zha 1996), in the same
    sentence or the next. **The addendum (`3f9ffbe`) rewords the rule:** "in
    the sentence that sets its rise against BCL-xL's, or the next".
  - **MCL1 is answered without BAD**, from the binding code (Chen 2005).
  - **PUMA relative to BCL-xL FALLS as OXPHOS priority rises** (ED-k). That
    is the reverse of the gland's Fig. 1i, and the text says so.
- **THE H2 CAVEAT, AD HOC and not declared** (result note, Part 2 H2).
  - **What it is:** part of S6's proliferation coefficient can be inherited
    from its selection rule, which partialled on proliferation.
  - **How much:** at the extreme, up to about half of the coefficient on 4g's
    model in TCGA, and a third in SCAN-B.
  - **What remains:** negative in every cell, but with no interval.
  - **MYC activity's coefficient is unaffected.**
  - **The test is the remainder's interval and a same-rule control
    signature.** It would need its own number, now `E51` or later.
  - **The rulings accept H2 and do NOT declare the test** (rulings 2.4).

**STOP 3 was then held.** The author ruled, then amended three rulings; the
amendment is E50's specification. **Everything STOP 3 asked is now settled**
by the consolidated rulings (`cb9b7ce`; see the 2026-11-04 entry above), and the binding
classes by the binding ruling (`66c205a`). What STOP 3 asked of the author:
- 4g's version, now replaced: there are no dots;
- the text (result note 4.1's "author's call" items, and the legend additions
  in 4.3);
- whether to declare the H2 test;
- the binding classes for `BID`, `BIK` and `BMF`. The author fills them in
  `data/annotations/bh3_binding_code.csv` and draws the code in Illustrator.

**Earlier the same day, before E49 ran: the author changed course on Fig. 4g**,
in `docs/2026-11-03_signature_and_fig4g_brief.md`. What follows is that
brief's record.

**What the brief decides:**
- **The six-gene signature (S6) returns as 4g's summary.**
  - Its OXPHOS coefficient is circular. It is drawn grey and never read.
  - Its proliferation and MYC coefficients carry the message.
- **4g becomes a heatmap in Fig. 1g's layout**:
  - the four groups of Fig. 1g as rows;
  - the BH3 binding code and the signature signs on the left;
  - a strip for the signature beneath.
- **The text reads PUMA against BCL-xL.**
  - BAD's rise is answered by its phosphorylation.
  - MCL1 is answered without BAD.
  - **Two claims from the chat are retracted** (brief 4.8). Neither is used
    anywhere.
- **It specifies `E49`**, now declared (`7bd7865`), run and verified:
  - Part A recovers the signature's quantities from E47's saved fits;
  - Part B fits the three effectors, BAX, BAK1 and BOK;
  - then the displays.

**An amendment of the same day changes the brief, and is recorded in the E49
declaration, not in the brief.** It changes three things:
- 4g's columns: three exposure blocks, from one model per row on the
  least-entangled MYC score;
- a placement rule for the MYC column;
- two sentences of section 8.

**So do not build 4g from the brief's section 7 alone.**

**What the brief supersedes:**
- **the 2026-11-03 handoff's sections 0 and 3.** E48's layout choice and its
  text choices are no longer pending;
- **E48 result Part 4 as the paper's text**;
- **three parts of the 2026-11-02 brief:**
  - 1.4's "the signature leaves the main text";
  - section 6, the 4g design;
  - section 8, the text drafts.

  **That brief's section 4 still governs every word**, with the new brief's
  section 5 added.

**E49's stages 1 and 2 are done**: the declaration and the binding-code table
(`7bd7865`), then the script and its pre-run digests (`a63ef6c`).

**Earlier the same day: E48 is RUN and VERIFIED, and the session stopped at
STOP 3 for the author to fix Fig. 4g's final design.** The new brief has
since replaced that STOP 3. What follows is E48's record.

**The entry point for E48's work was `docs/2026-11-03_handoff.md`.** Behind
it, in this order:
1. **`docs/2026-11-02_separation_line_brief.md`**, the author's decision
   record.
   - **Where it conflicts with anything below, the brief wins.**
   - **Its section 4 governs every word written about results.**
2. **`docs/2026-11-03_e48_result.md`**, the E48 result.
   - Its Part 2 holds all twelve genes' coefficients, which E49's 4g
     re-reads.
   - Its Part 4 was the text for the paper. The new brief's section 8
     replaces it.
3. **`docs/2026-11-03_e48_review_pack.pdf`**, the three layouts of 4g.
   - They are retired from the paper by the new brief.
   - They remain E48's verified record.

**For the whole earlier chain**, `docs/handoff-2026-10-21.md` is still the
consolidated record. Its predecessor `docs/2026-10-04_handoff.md` keeps
authority over what the section-4 figure panels may assert (its section 2).

**The line, as E48 left it (brief section 0):**
- in human breast tumours, OXPHOS, MYC activity and proliferation rise
  together;
- modelled together, OXPHOS and proliferation run in opposite directions, on
  `BBC3` and `BCL2L1` and on treatment response;
- MYC activity's association depends on how much proliferation its score
  carries.

**`E48` is RUN and VERIFIED**: 33 of 33 digests, committed before the run at
`1ccf805`, and 175 of 175 assertions. It is EXPLORATORY and POST HOC.
- **The declared reading, on `M_a` at `m3`:**
  - `BBC3` OPPOSITE;
  - `BCL2L1` OPPOSITE;
  - `MCL1` NOT SEPARATED.

  On the least-entangled score, agreement is 12 of 12. **Outcome row 3
  applies.**
- **The words that follow** (declaration section 8 and its addendum):
  - "opposite" for `BBC3` and `BCL2L1`, each named;
  - one descriptive sentence for `MCL1`, with no inference: proliferation is
    lower in TCGA and higher in SCAN-B;
  - **never "the BCL-2 family" for the three read genes.**
- **C3 holds at `m3`**, so "more steeply ... with the other two held fixed" is
  licensed.
- **C1 (TP53) and C4 (composition) change no label.** **C4's expectation was
  WRONG**: removing TCGA's purity and leukocyte adjustment moves `MCL1` away
  from SCAN-B's sign.
- **Nothing links the three genes to outcome**, and there is no causal,
  mediation, interaction or "primed" claim.
- **Part 2 of the result note holds four hypotheses, H1 to H4, each with a
  falsifier.** None of them is a finding.

**The trigger-arm claim is RETRACTED (brief 1.4) and stays out of the paper.**
- **The claim:** at fixed OXPHOS, MYC activity goes with the BH3-only
  (trigger) arm and not with `BCL2L1`.
- **It holds on `M_a` only.** On the least-entangled score, `BCL2L1` also
  rises with MYC activity, and `MCL1` falls. That is so at `m1` (E47) and
  with proliferation held fixed (E48).
- **E47's "The arms" bullet below is its record, not a licence.**

**Decision 5 is taken, for G2 co-amplification only** (the author,
2026-11-02).
- It admits one main-text sentence and, **since its amendment of 2026-11-05,
  one MAIN panel (Fig. 4x) in place of an Extended Data panel** (`8fde82a`).
  Nothing else from the validation study.
- **Its numbers are verified against the committed G2 note** (46 of 46), from
  the snapshot in `data/from_validation_g2/` (`c1dceed`). E54 draws the panel.

**E48's STOP 3 is superseded by the new brief.** It had the author choose
4g's layout and rule on the text choices in result note 4.1. The open items
it listed still stand:
- the author edits the Doc, which is never edited from here;
- the co-amplification display task, done by E54 (Fig. 4x; see its entry above);
- brief section 10's open items:
  - E46 in Extended Data, now settled: it is Extended Data Fig. 6d, drawn by
    `E56`;
  - signature transparency;
  - REFEREE_QUESTIONS: R8, R9 and a p53/PUMA entry;
  - the Fig. 4 title;
  - the carried-over items.

**`docs/STYLE.md` governs manuscript prose.** It is the author's, at v0.1.
Read it before drafting any sentence for the paper.

**What the brief superseded**, kept here so that the older entries below read
correctly:
- the reconciliation proposal's line and figure plan, whose corrections it
  adopts;
- the figure map's "Fig. 4g (`E46`) is MAIN":
  - **4g is now the separation panel, from E48**;
  - E46's quadrant map and E47 Display C go to Extended Data;
- the 2026-10-23 handoff's section 2 for 4g, and its section 7 item 1.

**The earlier entry point for Fig. 4 panel work and the narrative rewrite was
`docs/2026-10-23_handoff.md`.** The 2026-11-03 handoff supersedes it for 4g
and the human text. It superseded `2026-10-22_handoff_evening.md`.
- **Its section 5 is the input to the rewrite.** It lists every Doc sentence
  that must change, updated by E47.
- **Its section 6 points to the text that is ready to use.**
- **Its section 7 lists the open decisions.**
- **THEN READ `docs/2026-11-02_human_section_reconciliation_proposal.md`**,
  which updates the handoff's sections 5 and 7. It is a proposal, and its
  recommendations are the author's to accept.
  - It retires "tracks OXPHOS but not MYC" and gives text to paste into the
    Doc.
  - It recommends E47 Display C as main 4g, which would reverse the E46-main
    entry in the figure map below.
  - It records that "trigger arm, not guardian" holds on `M_a` only: on the
    clean estimator `BCL2L1` also rises with MYC (its section 8).
- **The E47 result is `docs/2026-10-23_e47_result.md`.**
- **The figure map, as the brief and E48 left it.** Its 4g entry is
  SUPERSEDED: 4g is now E49's heatmap (see E49's entry above), and E48's three
  layouts are retired from the paper. The rest stands:
  - **Fig. 4g was the separation panel from `E48`.** The author was to choose
    one of three layouts at STOP 3 (brief decision 3);
  - **the pCR and METABRIC forests (`E45`) are MAIN.**
    - Once co-amplification moves to Extended Data they become 4h and 4i
      (brief 6.1).
    - If layout 3 is chosen, the author decides whether they stay;
  - **E46's quadrant map goes to Extended Data** (brief open item 1), as
    Fig. 6d, redrawn by `E56`;
  - **E47 Display C goes to Extended Data as E48's ED-b**, relabelled, with
    its values unchanged;
  - **for the paper, E48's ED-c and ED-d replace E47's Displays A and B**
    (declaration 12.3). E47 Display A had been placed in EXTENDED DATA by the
    author on 2026-10-23;
  - the Doc's 4h, co-amplification, becomes one Extended Data panel under
    decision 5. **Amended 2026-11-05: it is a main panel, Fig. 4x, drawn by
    E54.**

- **THE NARRATIVE OF RECORD IS THE AUTHOR'S GOOGLE DOC** "11.8.26_MK_Myc_Paper"
  (Nat Metab version), from 2026-10-22. It supersedes v5 in `docs/`, **and its
  panel letters are not v5's**:
  - Doc 4g is the configuration by quadrant;
  - Doc **4h is MYC / BCL2L1 / MCL1 co-amplification**, which v5 placed in ED
    Fig. 9;
  - Doc 4i is the pCR forest (v5's 4h);
  - Doc 4j is the METABRIC forest (v5's 4i).

  Read the Doc through the Drive connector. It is not in the repo. Name figure
  files by content, never by letter.
- **Fig. 4i and 4j are BUILT** by **`E45`**, a display-only script that fits
  nothing. Ratios are on a log axis. In place of the proposal's bracket, a
  contrast column runs on every row (7.2). **VERIFIED: 9 of 9 digests, 27 of 27
  assertions**, and the output is byte-identical to the dry run. **Caveat: those
  digests were committed AFTER the author's run**, so git cannot show they came
  first (7.6). The tracked figure is
  `docs/figures/2026-10-22_E45_outcome_forests.{png,pdf}`. **For the paper,
  `E57` redraws them as two separate panels** (the ruling of 2026-11-06,
  `cb3c7d3`).
- **Fig. 4g is BUILT by `E46`, run by the author and VERIFIED**: 10 of 10
  digests, committed before the run at `d91e156` (proposal 8.4, 8.5). It
  shows the quadrant configuration on BOTH instruments, side by side, in the
  main figure: the OXPHOS split by GSVA and by mitoPPS, with `M_a` the same
  in both. Rows are blocked by function. The tracked figure is
  `docs/figures/2026-10-22_E46_quadrant_configuration.{png,pdf}`. **On
  mitoPPS, Q2 sits well below Q4**, so "comparable OXPHOS-high groups" holds
  in one block of four. **E47's verdict is in (MIXED). 4g's place in the
  paper is the author's call**: see E47 result note 9.3.
- **The provenance of the six signature genes** is the author's note
  `docs/2026-10-22_six_gene_composite_provenance.md` (`ceb06c9`). Three of
  its facts are verified against their sources. Most of its exclusion
  reasons are reconstructions written on 2026-10-22; the `BCL2` ER rule and
  `BID`'s negative-control role were written beforehand. **Its open item 3,
  whether `BID` stays, is a decision for the author BEFORE E47 runs.**
- **`E47` is RUN and VERIFIED**: 22 of 22 digests, committed before the run
  at `9251355`. **Verdict MIXED.** The result note is
  `docs/2026-10-23_e47_result.md`; read its section 9 before writing any
  sentence about MYC and the signature. The declaration is
  `docs/2026-10-22_E47_mutual_adjustment_declaration.md` (`dcb58a7`, addendum
  4a.7 at `3b5e221`). It fits `readout ~ OX + MYC + spine` on 4g's own
  scores, sets and covariates.
  - **The verdict reads `beta_MYC|OX` on the signature (S6, `M_a`).** Its CI
    is above zero in 4 of 8 primary cells and below zero in TCGA/GSVA `m1`.
    S5 and S7 are also MIXED, and order S5 < S6 < S7 in all eight cells.
  - **Behind MIXED, one structure.**
    - Proliferation adjustment raises `beta_MYC|OX` in every block; with it,
      MYC adds in 3 of 4 blocks.
    - **On the MYC estimator least entangled with proliferation, MYC adds in
      all 8 cells (+0.15 to +0.42).**
    - mitoPPS exceeds GSVA in every comparison.
    - RETRIEVED, not declared: `PROLIF`'s own coefficient is negative,
      -0.22 to -0.38. That is the mechanism of N9's sign flips, and it is NOT
      licensed for the manuscript.
  - **The arms.** At fixed OXPHOS, MYC goes with the trigger arm (`BBC3` and
    `BID` positive at `m3` in all four blocks), not the guardian arm
    (`BCL2L1` not consistent). **RETRACTED for the paper (brief 1.4).** It
    holds on `M_a` only; see the trigger-arm paragraph at the top of this
    section.
  - **D1.** On GSVA, 4g's MYC step is mostly or entirely leftover OXPHOS; on
    mitoPPS, mostly not.
  - **Section 9 is LICENSED**: at fixed MYC activity, OXPHOS tracks `BCL2L1`
    more steeply than `BBC3` in all four cells. The differences are +0.08 to
    +0.19; TCGA's lower bounds are 0.014 and 0.015.
  - **NOT licensed:**
    - "MYC adds";
    - "OXPHOS only";
    - **the Doc's "MYC association lost after adjustment for OXPHOS"**, which
      must be withdrawn.
  - **Five falsifiers** for the emerging trigger-arm hypothesis are in
    section 10 of the result note. The hypothesis is retracted from the
    paper.
  - **4g's placement is the author's call** (result note 9.3). Display A
    answers the question directly.
- **TERMINOLOGY, from 2026-10-22 (the author's decision):**
  - "the configuration composite" is now **"the signature"**: the
    BCL-2-family signature, S6;
  - its +1/-1 are **"signature signs"**: the direction a gene moves with
    OXPHOS, NOT its function.

  Earlier notes and code (`config_comp`) keep the old name. In new text, call
  a MYC activity score an **estimator**, never a "signature". **The provenance
  of every gene choice, including that no numeric rule chose S6, is in the
  E47 declaration, section 4a.**
- **The Doc's 4h is now Fig. 4x, drawn by E54** (see E54's entry above). It is
  validation-study material.
  - **Decision 5 (2026-11-02) admits it**, for G2 co-amplification only; its
    amendment of 2026-11-05 makes it a main panel.
  - Brief section 9's verification is done: the snapshot and the committed
    note agree, 46 of 46.
- **Doc sentences that clash with the data or with this file's rules** are
  listed in section 5 of the evening handoff. The weightiest:
  - "comparable" OXPHOS-high groups;
  - "MYC expression" where MYC activity is meant (trap 4);
  - "BCL2L1 more steeply than BBC3". **Now LICENSED**:
    - by E47 section 9, at fixed MYC activity;
    - by E48's C3, with the other two held fixed;
  - **"or BCL2L1 expression"** as an adjustment. **E40 has no `BCL2L1`-only
    rung.** Its `m3` adds `BUFFER_c`, the mean z of `MCL1` and `BCL2L1`, on top
    of MYC.

  **The Doc is the author's. It was read, never edited.**

**Phase 2 is CONCLUDED as of 2026-09-11.** Start at
`docs/2026-09-11_phase2_conclusion.md` — it states the finding, what is not
settled, and four falsifiers, and it names the three scripts behind it (`E24`
superseded in part, `E25` the result, `E26` the composition check). **No phase
is currently open.**

Completed: **Phase 1, the correlation atlas** (`E00`–`E05`, aim doc
`docs/2026-08-31_aim.md`, plan approved 2026-08-31), then the exploratory
`E06`–`E23` line, then **Phase 2, the cross-species mitoPPS rank comparison**
(`E24`–`E26`), then **`E27`**, a companion presentation of Phase 2 on the native
MitoCarta level-1 partition with an internal null. `E27` re-reads `E25` and
changes nothing in the conclusion; see
`docs/2026-09-11_e27_category_readouts.md`. Then **`E28`**, a post-hoc gap
readout that **did not survive its own tests** and carries nothing on the human
side; it changes nothing in the conclusion either, and
`docs/2026-09-11_e28_gap_readout.md` section 8 says exactly what it may and may
not be used for. Then **`E29`**, a **map and explicitly not a readout**: it found
position matches between groups to be coincidental and established no opposing
programme, and it records that **mitoPPS is not a closed composition** (the
per-sample sum ranges, so the standard compositional objection does not apply).
See `docs/2026-09-12_e29_oxphos_neighbourhood.md`. Then **`E30`**, a second
comparator for `E28`'s gap which **failed the estimator ladder identically** -
two comparators, same numerator, 2 of 4 estimators each - demonstrating that the
failure belongs to the rank-gap statistic rather than the comparator. See
`docs/2026-09-12_e30_mitoribosome_and_nadh.md`.

Then **`E32`** (2026-09-28), outside any phase: which Menegollo bicluster this
arm's MYC-high, OXPHOS-high state corresponds to. The direction came from
Menegollo's published text, and the rule was committed before
`forkscale_models.rds` was opened. **Verdict MB2-ALIGNED**:
`rho(M_a, MB2 | MB1) = +0.688`, `rho(M_a, MB1 | MB2) = -0.195`, so conditioning
on MB2 *reverses* MB1's apparent MYC association. It survives proliferation
(+0.653 on `FELSHER__PROLIFSTRIP`, halved to +0.375 with `PROLIF_DISJOINT` as a
covariate). **It licenses one manuscript sentence and nothing about outcome.**
The sentence and the three qualifications that travel with it are in section 8
of `docs/2026-09-28_mb2_forkscale_result.md`. Note that the *marginals* already
existed in `forkscale_models.rds`; what `E32` adds is the partials.

**`E31` is not on `main`, and the gap is deliberate.** It is a BH3-mimetic
pharmacogenomics analysis (GDSC 8.5, PRISM 24Q2) on the unmerged branch
`bh3-mimetic-oxphos`, kept on origin so that nobody re-runs it without knowing
it has been done. Its declared primary is a precise null (an RB1 positive
control fires at p = 4e-13 in the same models) and is uninformative about the
manuscript, because `E18` found the BCL2L1-OXPHOS configuration absent from
cell lines. Read it with
`git show origin/bh3-mimetic-oxphos:docs/2026-09-28_bh3_mimetic_result.md`.

**`E33` and `E34` ARE on `main`**, merged 2026-10-01 at the author's request
once both were complete: declared alone before any data was opened, scripted,
run by the author, and verified object by object against digests recorded
*before* that run. Each kept its own branch on origin, undeleted, and each
arrived as a merge commit so its declaration / data note / unrun script /
result sequence survives as the audit trail.

- **`E33`**, merged from `e33-programme-specificity`. Is the BCL2-family death
  configuration specific to the MYC-coupled bicluster? **Verdict NEITHER
  SEPARATES.** Trigger arm `rho(arm, MB2 | MB1) = -0.140 [-0.200, -0.078]`
  against `rho(arm, MB1 | MB2) = +0.250`; guardian arm -0.350 against +0.517;
  unchanged on Block C, on the patient join and under the stricter reading. **It
  licenses no new sentence and WITHDRAWS one implication**: the manuscript may
  not state or imply that the standing pro-apoptotic configuration is confined
  to tumours whose respiration is MYC-coupled. A post-hoc observation that the
  configuration leans **MB1** is recorded with its own falsifier list and is not
  a finding. **Do not chain `E32` to `E33`**: the *state* is MB2-aligned, the
  *configuration* is not, and a property of the state is not a property of
  everything it correlates with. It is `docs/2026-09-30_e33_result.md`.
- **`E34`**, merged from `e34-quadrant-configuration`. The same configuration across
  the four MYC x OXPHOS quadrants, both cohorts, both instruments. **Verdict
  MIXED** — the rule demanded all 18 reading cells agree and 16 did; both
  dissents are mitoPPS within a subtype. What is new is the **Q2 cell**, MYC-low
  and OXPHOS-high, which `STATE` collapses into its level 1: it is large and
  carries the configuration with Q4 rather than Q1 (adjusted means +0.41 against
  +0.40 in TCGA, +0.24 against +0.37 in SCAN-B). **The MYC contrast is small but
  not zero**, and **proliferation adjustment enlarges it** by a uniform +0.12.
  Neither pre-written reading is licensed; the four-group description is.
  **The two OXPHOS contrasts are near-guaranteed by construction** — the
  configuration's genes and signs were read off OXPHOS in these same cohorts —
  so only the MYC contrasts carry information. It is
  `docs/2026-10-01_e34_result.md`.

**One `E33` input finding qualifies `E32`, which is on `main`.** The Menegollo
forkscale should be joined on the **aliquot** barcode, not the patient barcode:
six TCGA patients were profiled from a different vial, so `E32`'s patient-key
join paired those six across vials. On the aliquot join (n = 1,031; Block C
overlap 879) `E32`'s partials move by less than 0.01 and **its verdict and its
manuscript sentence stand** — but any future forkscale join should use the
aliquot.

**`E35`, `E36`, `E37` and `E38` are all on `main`**, each merged once complete,
run by the author and verified object by object against digests recorded
*before* that run — 25 of 25 for `E36`, 15 of 15 for `E37`, 30 of 30 for `E38`.
Branches kept undeleted on origin.

They continue one line of questioning and they finish it: `E34` asked where the
configuration sits, `E35` whether that depends on genomic burden, `E36` whether
it exists below the tumour range at all, and `E38` diagnosed why `E36` could not
answer. **The normal-tissue line is closed.**

**`E39`, `E40`, `E41`, `E42` and `E43` are on `main`**, committed straight to
the trunk rather than merged from branches, each verified object by object
against digests recorded *before* the author's run — **26 of 26 for `E41`, 16 of
16 for `E42`, 14 of 14 for `E43`, 9 of 9 for `E44`**. **The next free script
number is `E58`**.
- `E45` and `E46` are the Fig. 4 display scripts (4i/4j and 4g).
- `E47` is run and verified (MIXED).
- `E48` is run and verified, 33 of 33. It gives the separation reading and
  4g's three layouts.
- `E49` is run and verified, 24 of 24. It gives the signature's declared
  quantities, the effectors, and the new 4g in two versions.
- `E50` is display-only, run and verified (116 of 116 assertions; 11 of 12
  digests, the twelfth a label's wording). It redraws 4g, ED-e and ED-f in
  Fig. 1g's style and re-colours ED-g and ED-h. The STOP 3 rulings accept its
  displays; the rest is the author's Illustrator work, and no display script
  follows it for 4g.
- `E51` is a display script with one declared quantity, each cohort's *P*.
  It redraws ED-k as a main panel, with both cohorts overlaid. It is run and
  verified: 8 of 8 digests, 26 of 26 assertions.
- `E52` is a display script. It redraws E51's panel at a 26-mm axis, with
  the *P* labels in a key above the data area and no "adjusted" in the axis
  titles. It is run and verified: 8 of 8 digests, 38 of 38 assertions. Its
  panel replaces E51's for the paper.
- `E53` is a display script. It redraws E48's ED-a and ED-b (Extended Data
  Fig. 6a and 6b) at half width, with 6b's arrows removed and its labels
  edited. It is run and verified: 6 of 6 digests, 43 of 43 assertions. Its
  panels replace E48's ED-a and ED-b for the paper.
- `E54` is a display script. It draws Fig. 4x, the validation study's G2
  co-amplification odds ratios, from the snapshot in
  `data/from_validation_g2/`. It is run and verified twice; the second
  version, with smaller circles, is Fig. 4x.
- `E55` is a display script. It draws Extended Data Fig. 6c, E49's
  unadjusted correlations in Fig. 4g's grammar. It is run and verified:
  7 of 7 digests, 22 of 22 assertions.
- `E56` fits the effectors' quadrant means (declared, post hoc) and draws
  Extended Data Fig. 6d. It is run and verified: 11 of 11 digests, 25 of
  25 assertions.
- `E57` is a display script. It redraws E45's two outcome forests as two
  separate panels, Fig. 4's last two, at Fig. 4x's size and in its slates.
  Its first run is verified (7 of 7 digests, 75 of 75 assertions). The
  author then amended it: black text, and dashed connectors. The second
  version is run and verified (7 of 7 digests, 80 of 80 assertions), and its
  panels are Fig. 4's last two.

**`E31` is the only analysis still off `main`. No phase is open.**

They are **one analysis under one declaration**,
`docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md` (amended for
METABRIC at `ccadfa5` and `82659ff`), with a division of labour fixed in its
section 10: GSE25066 carries the ER-negative sign pair, METABRIC the
ER-positive half. **One result note covers all three**:
`docs/2026-10-03_e39_e40_e41_result.md`. It exists because the draft sentence
*"respiration did not itself predict response"* had **no fitted quantity behind
it** — only three numbers printed once in an `if (FALSE)` sandbox.

- **`E39`**, the identifiability gate. It **fits nothing**.
  `rho(OX, PROLIF_DISJOINT)` = 0.3575 / 0.3812 / 0.4934, all below 0.60, so the
  ladder is licensed in all three cohorts. **But it found something it was not
  looking for**: `rho(MYC, PROLIF_DISJOINT)` = 0.7797 / 0.7851 / **0.8236**, and
  the last is inside the declaration's own `>= 0.80 NOT IDENTIFIABLE` band. **In
  GSE164458 the MYC and proliferation terms are not separately identifiable**,
  so 5.3's prohibition on reading `m1 -> m2` as "what remains after
  proliferation is MYC" is **mandatory, not cautious**.
- **`E40`**, the ladder on pCR and DRFS, three neoadjuvant cohorts.
  **ONE QUANTITY MEETS THE DECLARED READING RULE**: the pooled pCR coefficient
  at `m1`, **-0.1616 [-0.281, -0.042]**, negative in all three cohorts, I2 = 0
  (and `m2`, `m3` likewise). **`m0` fails both conditions** — signs `+ - +`,
  pooled CI including zero, the only cohort-level `m0` excluding zero being the
  smallest cohort disagreeing in sign with the largest. **So the draft sentence
  refers to a quantity that is not readable in either direction and must be
  replaced or dropped.** Every VIF for OX is 1.00 to 1.43, so the
  MYC-proliferation collinearity **never reached the exposure**.
- **`E41`**, the METABRIC prognosis ladder — **1,505 ER-positive patients
  carrying 431 cause-specific events, against GSE25066's 40**. Nothing meets the
  rule. It licenses the **forkscale** statement: `rho(OX, MB1_forkscale)` = 0.416
  here against F3-pre's 0.529 (GSVA) and 0.418 (mitoPPS), so three readings on
  two instruments in two cohorts all sit inside the INTERMEDIATE band and
  **INTERMEDIATE is a property, not an instrument artefact**. `F3` was specified
  in `myc_human_validation`; **these rungs are NOT F3.**
- **`E42`**, the subtype scope check on `E40`'s primary finding. **It cannot
  promote or demote it** (13.2) and did neither. **The draft sentence stands
  UNQUALIFIED**: 13.6's third row is the only claim-changing one and it needs
  opposite signs *and* non-overlapping CIs — all three subtype summaries are
  negative and all three intervals overlap. **But 13.6 had no row for what
  happened.** The subtype whose pooled estimate excludes zero is
  **HR-positive/HER2-negative** (-0.2759 [-0.5299, -0.0219], I2 3.9%), not TNBC,
  which is **heterogeneous rather than underpowered** — 441 events, the best
  powered stratum in the study, at I2 62.1% with its three cohorts' signs
  disagreeing. **13.6 was NOT amended**: writing a row after seeing the result
  would be fitting a rule to an outcome, so the gap is recorded in section 5.2
  of `docs/2026-10-03_e42_result.md` and 13.6 carries a dated pointer with its
  content unchanged. `E42` also **refuted the `s508` explanation**:
  `rho(OX, subtype2) = -0.054`, and near zero was the one decisive direction,
  so **that divergence now has no explanation.** And **one of four subtypes was
  never fitted** (`HRneg_HER2pos`, minority class 33 against 7 parameters), so
  **the clause's generality over HER2-positive disease rests on nothing
  measured.**
- **`E43`**, a **RETRIEVAL and not a new analysis** (14.2). `E40` and `E41`
  carried `PROLIF_DISJOINT` in every rung from m1 and extracted **only the OX
  row** of each coefficient matrix, so the proliferation coefficient was
  computed and discarded. `E43` refits m1 from the saved frames — **the gate
  reproduced every stored OX coefficient and SE with a discrepancy of exactly
  0**, so these are the same models — and retrieves it. **Declaration 14.6 lands
  on row 1**: PROLIF is **positive on pCR in all three cohorts**, every interval
  excluding zero (pooled random **+0.5422 [0.307, 0.778]**, I2 66.9), and
  **above unity on survival** with intervals excluding the null in the pooled
  ladder (**HR 1.1336**) and both ER-positive strata. `docs/2026-10-08_e43_result.md`.

  **The gain is INFERENTIAL, and it is the reason to care.** PROLIF reproduces
  the canonical proliferation paradox — better pCR, worse survival — which is
  **exactly row 1 of section 9.1's sign-pair table**, inside the same models in
  which OX runs the other way. **So the models CAN detect a known sign pair in
  these data**, and OX's failure to show the canonical pattern is not a failure
  of the models, the spine, the cohorts or the endpoints. **It functions as a
  positive control; it was NOT designed as one** and no sentence may imply it
  was. **The asymmetry travels with it**: PROLIF's two halves are both solid and
  **OX's are not** — its survival limb is still `E42`'s boundary cell, and `E43`
  does not change that by one decimal.

  **`14.5`'s pre-stated expectation was WRONG IN THE CONSERVATIVE DIRECTION**
  and is recorded as such: it prepared for an attenuated or null coefficient,
  and the coefficient survives subtype adjustment in both arms.
- **`E44`**, the contrast between them, **TESTED rather than eyeballed** and
  computed with **NO REFIT** from `E43`'s saved fits. `delta = beta(PROLIF) -
  beta(OX)` from the same m1 fit, with `Var(delta) = Var(bP) + Var(bO) -
  2*Cov(bP, bO)`. **15.6 lands on row 1**: delta excludes zero on both
  endpoints with signs agreeing in all three pCR cohorts - pooled **+0.6577
  [0.4446, 0.8709]** at I2 0, a ratio of odds ratios of 1.930, and **+0.3160
  [0.1180, 0.5140]** in METABRIC's ER-positive IHC stratum. **So the two
  coefficients DIFFER, tested, on both endpoints.**
  `docs/2026-10-08_e44_result.md`.

  **IT REPLACES THE BASIS OF THE CLAIM, NOT THE CLAIM.** 15.7 measured before
  the computation that the eyeball non-overlap **fails** in the pooled ladder
  and in the declared `ER.Expr` sensitivity; the tested contrast **holds** in
  both. **So "the two intervals are again disjoint" MUST NOT be written on
  distant outcome** - it is true of the ER-positive IHC stratum and nothing
  else there - and the tested difference goes in its place.

  **AND 15.3 GOVERNS IT ABSOLUTELY. A significant contrast does not make either
  coefficient a finding.** The OX coefficient remains exactly as weak as the
  E39/E40/E41 result note's 7.2 records it, and **nothing in `E44` may be cited
  in support of a distant-outcome claim for respiration.** What it licenses is
  that the two DIFFER, not that either is established.

**Three things from this line that must travel with any sentence about it.**

1. **ON DISTANT OUTCOME NOTHING MEETS THE RULE — FOR THE OXPHOS EXPOSURE.**
   **Say which exposure**, because since `E43` the sentence is ambiguous without
   it: `PROLIF_DISJOINT` in those same models IS above unity with intervals
   excluding the null in the pooled METABRIC ladder and both ER-positive strata.
   What follows is about **OX**, in either cohort family, and
   **two near-misses are traps rather than results.** The only GSE25066 DRFS
   cells excluding zero are in **`s508`, the set with NO SUBTYPE TERM** — the
   pre-specified 465 and the 490 have none — so the coefficient crosses the line
   exactly when subtype adjustment is removed. **Those cells may not be cited.**
   And METABRIC's ER-positive `m1` excludes zero by **3.7e-04**, its own declared
   sensitivity **does not reproduce it** (p 0.186, and the 127 reclassified
   patients are 28 Basal and 32 Her2 leaving the stratum), and `m2` and `m3`
   include zero. **One boundary cell of thirty is not a finding**, and section
   9's ban on calling a p value a near miss cuts **both ways** around it.
2. **METABRIC DOES NOT REPLICATE `E40`, and no sentence may say it does** —
   declaration 10.5, whatever the sign. Different estimand: of 1,505 ER-positive
   patients **1,111 (73.8%) had endocrine therapy and 147 (9.8%) chemotherapy**,
   so that stratum is overwhelmingly a non-chemotherapy one and E40's endpoint
   does not exist in it. **The two spines also differ** (10.1a) — PAM50 there, a
   two-level ER/HER2 collapse in GSE25066 — so the two `m0` are not the same
   model.
3. **SINGLE-INSTRUMENT throughout.** mitoPPS exists in none of the four cohorts;
   three are arrays or panels and METABRIC is an array. **This line has not
   cleared the arm's two-instrument bar and must never be written as though it
   had.**

**Two traps this line added to the repo's stock, both measured.**

- **An alias is not automatically safe, and the guard belongs on the ROUTE, not
  the symbol.** Harmonising to METABRIC's pre-2018 build, `COX7A2L -> SCAF1` and
  `POLR1B -> RPA2` both send a gene to the **HGNC-approved symbol of a different
  gene** (58506 and 6118); taking the first would have put a splicing factor
  inside the OXPHOS exposure. Both are rejected, which is why OXPHOS subunits
  recover **69 of 89 and not 70**. **But a blanket ban on those symbols is also
  wrong**: `RPA2` is a legitimate *direct* member of `PROLIF_DISJOINT`, and
  banning it would have silently altered the proliferation score in every model.
  Declaration 10.3a carries the two-part rule.
- **NEVER DIFFERENCE TWO COEFFICIENTS FROM ONE MODEL WITHOUT THE COVARIANCE.**
  `Var(bP - bO)` is `Var(bP) + Var(bO) - 2*Cov(bP, bO)`, and the naive sum of
  squared SEs assumes the covariance is zero. In `E44`'s eight fits
  `Cov(bP, bO)` is **negative in all eight** (coefficient correlations -0.512
  to -0.312), so the correct SE is **14.3% to 22.9% LARGER** and **the naive
  version is ANTI-CONSERVATIVE, not merely wrong**: it would have given p
  0.000352 instead of 0.00176 in the primary stratum and 0.0196 instead of
  0.0434 in the pooled ladder. **Use a general linear contrast** - `L %*% beta`
  with `Var = L V L'` - so the covariance sign comes from the algebra. And
  **the difference of two separately POOLED estimates is not the pooled
  difference and has NO valid standard error**, because the pooled estimates
  come from overlapping patients; `E44` landed 0.6566 against the correct
  0.6577, which is luck rather than licence.
- **A MULTI-DF GVIF IS NOT ON THE SAME SCALE AS A 1-DF VIF, and section 5.3's
  threshold of 5 was set for OX, which is 1 df.** Comparing a raw multi-df GVIF
  against 5 overstates it. Fox and Monette's prescription is to compare
  **`GVIF^(1/(2*Df))`** against the **square root** of the 1-df threshold.
  `E43` flagged GSE194040's subtype GVIF 7.979 and treatment GVIF 7.269 as
  breaching the caveat; adjusted they are **1.414** (3 df) and **1.086** (12
  df), **no term in any model of the line exceeds `sqrt(5)` = 2.236**, and the
  maximum anywhere is 1.463. **So the line carries NO collinearity caveat**, and
  the raw figures concern the **spine** — I-SPY2's arm assignment by receptor
  status — rather than either exposure. The error was only visible because
  `14.7` required the full VIF vector.
- **`BBC3` is ABSENT from METABRIC.** It costs this ladder nothing, since the
  configuration is not fitted there — but `BBC3` is the trigger limb the
  guardian-switch finding is built around, so **any future METABRIC analysis
  touching the 6-gene signed configuration is blocked on it** and a five-gene
  stand-in is forbidden.

**`E40`'s 46 unexplained warnings are identified and the item is closed**: they
are `car::vif`'s "No intercept: vifs may not be sensible", spurious for a model
that has no intercept by design, at 40 fits times a double call. **No `E40`
number is affected.**

### `X01` is OFF the E-numbering, and that is the point

**`X01` is grant-preparatory and MANUSCRIPT-EXCLUDED.** It is deliberately not
an `E` number so the file listing itself shows the difference: the E-line is
manuscript-bound and reports whatever it shows under a declaration committed
before it ran, and `X01` is neither. Its governing note is
`docs/2026-10-03_X01_exploratory_leads_note.md` and its result is
`docs/2026-10-04_X01_result.md`.

**NO result from `X01` may enter any manuscript document, figure, legend or
supplementary file without its own declaration written and committed first**,
and **`E42` 13.2's prohibition on a subtype result becoming a finding is NOT
relaxed by it.** The vocabulary there is *suggests*, *is consistent with*,
*would be worth testing*.

**It surveyed METABRIC's ER-positive stratum by PAM50 and nothing cleared the
bar its own note set** — no subtype meets more than one of three criteria. Four
things it established, all negative or cautionary and all useful:

1. **`LumB` carries the ER-positive association and `LumA` does not** (LumB m1
   0.7532 [0.634, 0.895], the only interval excluding 1; LumA flat at both
   rungs). `Her2` scores three of three on the false-lead signature named in
   advance and is not a lead.
2. **The mechanism its own note predicted failed.** It named `LumA` as most
   differentiated and most oxidative; the association sits in the *more*
   proliferative luminal subtype. **This is the second mechanism this arm has
   proposed and had refused** — `E42`'s `rho(OX, subtype2) = -0.054` was the
   first — and **the pattern is worth more than either refutation.** No
   replacement mechanism was offered, deliberately.
3. **A criterion imported from the E-line inverted.** "m0 null and m1
   protective shows the adjustment doing work" is right for the E-line's
   mechanistic question and wrong for a biomarker one, where an unadjusted
   association is the deliverable and adjustment-survival is evidence of
   *incremental* value over grade and Ki67. On the reversed pair `LumB` ranks
   first where the criterion ranked it last. **Which criterion to apply
   changed; what the estimate is worth did not.**
4. **The predictive-marker hypothesis cannot be surveyed in METABRIC.** The one
   fittable chemotherapy group reverses sign (HR 1.348 [0.848, 2.142]) against
   0.828 to 0.880 elsewhere, but **confounding by indication is the leading
   explanation** — that group is the most proliferative (mean `PROLIF` -0.048
   against -0.147 to -0.383, *which is the indication*), has 5.08 y follow-up
   against up to 11.63 y where PH already fails for OX, and is depleted of
   `LumB`. **It cannot be tested from disk**: `grade`, `size`,
   `lymph_nodes_positive` and `age_at_diagnosis` were never carried into
   `E41`'s frame. That gap is one of the two blocks behind R rule 4.

- **`E35`**, merged from `e35-burden-coupling`, **complete and verified**. Is the
  OXPHOS-to-configuration coupling graded by genomic burden? **Reading
  SATURATED** — the declared expected outcome. On `Aneuploidy.Score`, m1, GSVA:
  `Q2 - Q1` +1.073 / +0.651 / +1.144 and `Q4 - Q3` +0.929 / +0.853 / +0.745
  across tertiles, all six intervals excluding zero and the lowest and highest
  overlapping. **11 of 12 readable blocks agree.** Adding `PROLIF_DISJOINT`
  moves it by at most 0.054, so it is not proliferation under a burden label.
  **It licenses one clause and reinstates no stratifier**, and it cannot
  demonstrate the mouse's condition — a cohort of tumours has no sample below
  the threshold. `docs/2026-10-01_e35_result.md`.
- **`E36`**, merged from `e36-normal-comparison`. Is the coupling present in
  normal breast — the human analogue of the mouse wild-type gland, and the one
  analysis that could come back against the model. **Verdict UNINTERPRETABLE.**
  Its declared positive control — the **Fe-S cluster assembly cytosolic half**,
  named before any fit, with zero genes in common with either the configuration
  or the OXPHOS arm — reproduces **69% to 96%** of the configuration's
  normal-to-tumour change, in 6 of 12 cells. **The verdict does not rest on the
  control alone**: independently, all six cells are SPLIT BY ADJUSTMENT.
  **Q1 LICENSES NOTHING** — not ABSENT, not WEAKER, not PRESENT IN NORMAL. The
  manuscript may not say the coupling is absent in normal breast, nor present,
  nor that the human data do or do not corroborate the mouse gland result.
  **Remember why that matters: every *unadjusted* cell reads ABSENT IN NORMAL**,
  which is the predicted result and what a less careful analysis would have
  reported. Only the declared control and the adjustment-agreement rule stopped
  it. **The declared `EPITHELIAL` adjustment failed its own validity check** —
  its luminal and basal halves separate the tissues in opposite directions (AUC
  0.875 against 0.134) and cancel to 0.445 — **and the gene list was not changed
  and no substitute introduced.**
  **Q2 does not depend on the control and STANDS.** Over 113 matched pairs the
  tumour exceeds its own normal in 70% on `M_a`, but the two
  proliferation-stripped CollecTRI regulons **reverse the sign** to 39.8%. By
  quadrant: **Q1 48.7% and Q2 47.1% against Q3 100% and Q4 90.6%** — MYC-low
  tumours have no more MYC activity than their own matched normal, so the
  threshold reconciliation **fails for them**. The comparison is **paired** (all
  113 normals come from patients already in the tumour set) and the normals come
  from **6 of 40 collection sites**, 87.6% from three, so every across-person
  comparison is reported on **three nested tumour sets** and the within-patient
  contrast is the one immune to it. `docs/2026-10-01_e36_result.md`.
  **Its section 3 carries a dated correction from `E38`**: the numbers stand,
  two interpretive sentences do not, and the verdict is unaffected.
- **`E37`**, the snapshot rebuild `E36` needed, **GATE PASS**. It rebuilds the
  expression layer over tumours **and** normals jointly, because DESeq2 size
  factors over 1,208 samples are not size factors over 1,095 and including
  normals moves every tumour value. Its gate re-derives E34 on the rebuilt
  snapshot, tumours only: **quadrant agreement 97.1% / 98.1%, all sixteen
  contrasts within 0.10 (largest 0.047), E34's four pooled labels unchanged.**
  **So E34 is robust to a renormalisation that moves every input value.** It
  licenses that one sentence and **nothing biological**.
  `docs/2026-10-01_e37_result.md`.

- **`E38`**, **a DIAGNOSTIC and explicitly not a hypothesis test**.
  It decomposes quantities `E36` already computed, **adds no claim and cannot
  overturn `E36`'s verdict**; two of its four legs are **second attempts after
  seeing a first result** and are labelled post-hoc throughout. Its arbitrary
  sign splits were named **gene by gene** before it ran, and the script asserts
  them. Four findings:
  **(1) The signed near-zero is NOT cancellation.** Unsigning the configuration
  moves its normal-side correlation from -0.011 to **-0.037**, and arbitrarily
  signed comparators do **not** collapse (+0.834, +0.457) even at balanced
  splits. The leading explanation for `E36`'s result is refuted.
  **(2) The compartment gap is present in normals and LARGER there**, and is not
  specific to the apoptotic machinery — Fe-S shows +1.022 against +0.754. No
  reading taken; E14's compartment-matched null was not rebuilt.
  **(3) A luminal-only epithelial score FAILS** `E36`'s criterion applied
  unchanged, by **+0.048 on one leg of three**. Not used; criterion not relaxed.
  **(4) `MYC` mRNA is DOWN in 85% of matched pairs** (median -1.394 log2) while
  MYC *activity* is up in 70% of the same pairs — **trap 4 in its sharpest
  form** — and the 82 targets `M_b__PROLIFSTRIP` removes are **4.9x enriched in
  `HALLMARK_MYC_TARGETS_V1`** and 19.6x in `YU_MYC_TARGETS_UP`, including
  `ODC1`, `NCL`, `NME1`, `TFRC`, `CCND1`, `CDK4`, `E2F1`. **So
  `M_b__PROLIFSTRIP` is not "MYC activity with growth removed" and must not be
  described as such.** It licenses **three methodological statements and no
  biology**: `docs/2026-10-02_e38_result.md` section 7.

**`results/joint_tcga_*.rds` is gitignored and must be rebuilt after a fresh
clone**, in this order: `E37`, then `E36`, then `E38` — about five minutes in
all. **`E36` and `E38` stop with a clear error if the snapshot is absent.** The
tracked figures in `docs/figures/` and the numbers written into the notes are
the durable record. See `docs/2026-10-02_handoff.md`.

**One thing the normal-tissue line established that outlives it.** Every
limitation in it reduces to **bulk tissue composition**, and **no composition
adjustment built from bulk markers has yet survived its own validity check** —
the declared nine-gene score failed, the luminal-only score failed by +0.048.
`E38` also names an instrument this repo does not have: **a signed control**,
matched arbitrary genes with matched signs, which would separate "signed" from
"apoptotic".

Named as decisions rather than drift, and still not done: MCbiclust / forkscale
beyond `E32`'s alignment check and `E41`'s two rungs (the Menegollo axis
proper), treatment as anything but a spine term, DepMap, causal or mediation
modelling, and anything that revisits the validation study's hypotheses.

**Two of these came off the list with the `E39`-`E41` line and the change is
narrow.** **Survival** is no longer undone, but only as `E41` fitted it —
METABRIC breast-cancer-specific and overall survival, on one instrument;
SCAN-B carries no endpoint and TCGA survival is still forbidden by plan
section 3. **METABRIC** is no longer undone, but only for that ladder and the
two forkscale rungs; its expression layer is untracked and its provenance rests
on the author's authority (Synapse `syn1757063`, nothing inside the file naming
it, md5 the only check).

Two items left open on purpose, recorded in
`docs/2026-09-11_handoff.md` (spent, kept for this): whether standalone
mouse-only figure reproductions are wanted — they would be a `myc_mouse` task
and cannot be built from here — and whether `data/from_myc_mouse/` should be
tracked rather than gitignored, since it is 511 KB and the convention it copied
exists for a 560 MB directory.

## Section 2 — the traps

Each is measured, not anticipated. Numbers are from the snapshot.

1. **Correlation is not the interaction.** The validation study found the
   `MYC x OXPHOS` *interaction* on apoptotic priming is null. A strong
   MYC–OXPHOS *correlation* is entirely compatible with that. They are different
   questions. Never present one as confirming or contradicting the other.
2. **Purity and immune infiltrate.** In TCGA `rho(OXPHOS subunits, purity) =
   0.214`, `rho(., leukocyte fraction) = -0.158`. Breast is the worst TCGA
   tissue for this: adipose is OXPHOS/FAO-high and infiltrate carries its own
   BCL2-family profile. **SCAN-B has no purity estimate.** Report raw and
   adjusted in TCGA, raw only in SCAN-B, and say so on the figure.
3. **Every MYC activity signature is entangled with proliferation, by wildly
   different amounts** — 1.5% (`MYC_UP.V1_UP`) to 47.6% (`YU_MYC_TARGETS_UP`),
   with the validation study's `M_a` at 14.8% and `HALLMARK_MYC_TARGETS_V1` at
   23.5%. Never report a MYC–OXPHOS correlation from one signature. Report the
   panel ordered by entanglement, and the proliferation-adjusted estimate.
4. **MYC mRNA is not MYC activity, and the difference is total.** In TCGA
   `rho(log2(MYC), OXPHOS subunits) = -0.032` against `+0.388` for the activity
   signature. Someone plotting MYC expression against OXPHOS sees nothing. Both
   are reported and the gap is itself a result.
5. **The four instruments disagree, by a lot.** GSVA-vs-mitoPPS agreement across
   arms runs 0.24 (lipid metabolism) to 0.94 (mtDNA). Instrument choice is not
   cosmetic: report all four, or justify the one.
6. **mitoPPS is blind to level by design.** It answers "is OXPHOS *prioritised*
   relative to other mitochondrial programmes", not "is OXPHOS high". Never
   compare mitoPPS values numerically across cohorts — only patterns.
7. **SCAN-B's symbols are a 2014 UCSC build.** 19 of the 89 `OXPHOS subunits`
   genes are pre-2018 ATP-synthase names; unharmonised the exposure covers 0.775
   instead of 0.989, and 70 of 89 genes is a Complex V with no F1 head and no
   c-ring. **Never score SCAN-B without `scanb_pheno.rds$symbol_map`.**
8. **mtDNA-encoded genes are held separately** and never pooled with
   nuclear-encoded subunits — expression-scale skew, and they behave differently:
   `rho(M_a, mtDNA-encoded OXPHOS) = +0.068` against `+0.388` for the nuclear
   subunits, and **negative on mitoPPS**.
9. **CICD is thin and must not be over-read.** 13 pro-death and 4 pro-survival
   human genes. It is the axis of most interest and the weakest measured. Score
   only what clears n >= 5; show the genes individually; never present a 4-gene
   GSVA score as a programme.
10. **Two Tang sets are near-transcriptome-wide**, and their sizes are easy to
    overstate. The CSVs carry one row per gene-per-evidence, so `Ferroptosis` is
    935 rows but **600 genes** and `Autophagy_dependent` 1,195 rows but **876** —
    4.8% and 3.3% of the matrix. Always count distinct genes. A correlation with
    either is still close to one with general expression: read them against a
    size-matched comparator before believing anything.
11. **The cell-death and MYC sets are human-native and must never be remapped.**
    Both carry first-class human columns and the upstream README says "Do NOT
    remap". See `data/genesets_celldeath_human/README.md`.

## Scale discipline — the most likely silent error

- **GSVA / ssGSEA** want log-scale input: VST, `kcdf = "Gaussian"`.
- **mitoPPS** wants linear DESeq2-normalised counts.

Opposite requirements. **They must not share an input object.** State the scale
in a comment at the top of every scoring block.

- **GSVA is cohort-relative.** Score all samples of a cohort in one run, and all
  sets of interest in the *same call* — the `.PIN_A`/`.PIN_B` half-matrix pins
  hold the gene universe, and two calls with different set collections are not
  comparable without them. **Never pool scores across cohorts**; compare
  correlations and patterns, not values.
- **mitoPPS baseline is composition-dependent.** It reports the *shape* of the
  mitochondrial programme, not its level.

### Two normalisations of TCGA-BRCA now exist. Say which one you read.

From `E37` there are **two** normalised TCGA expression layers in this repo, and
**a value from one is not comparable with a value from the other** — the DESeq2
size factors and the low-count gene filter are both computed over a different
sample set, so every value differs.

| snapshot | what | read by |
|---|---|---|
| `data/from_validation/` | 1,095 tumours, 18,115 genes | `E11`, `E33`, `E34`, `E35` and everything before them |
| `results/joint_tcga_*.rds` | 1,208 samples (1,095 tumour + 113 normal), 18,142 genes | **`E36` and `E38`, and nothing else unless separately declared** |

- **Every script states which snapshot it reads**, in a comment at the top of
  its scoring block, exactly as the scale rule above already requires.
- **No figure, table or sentence combines a quantity from one with a quantity
  from the other.** `E37`'s gate is the sole exception: it is a deliberate
  like-for-like comparison of the *same* quantities and is labelled as such.
- **Neither replaces the other.** Replacing `data/from_validation/` would
  invalidate `E11`, `E33`, `E34` and `E35` at a stroke.
- **One genomic annotation crosses, and it is not an expression value.**
  `E36` and `E37` carry `BUFFER_gistic`, and `E38` reads `MYC_amp`, from the
  frozen covariate table. These are per-patient GISTIC calls, not derived from
  any expression normalisation. **No expression value crosses between
  snapshots.**

## Gene sets — consume the snapshots, do not rebuild

Each directory under `data/` carries its own provenance README with a pinned
commit SHA. Consume as-is; do not rebuild here, do not edit in place. If a
source changes, re-snapshot and bump the SHA rather than patching.

| Input | Location |
|---|---|
| Validation-study matrices, scores, covariates | `data/from_validation/` (gitignored, ~563 MB) |
| Validation study's G2 result, for Fig. 4x | `data/from_validation_g2/` (tracked, 10 KB; md5 in its README) |
| Human MitoCarta 3.0 | `data/mitocarta_human/` |
| Cell death: pro/anti x apoptosis/CICD, + 15 Tang modalities | `data/genesets_celldeath_human/` |
| MYC signature compendium (16 sets) | `data/genesets_myc_human/` |
| CollecTRI regulons | `data/collectri_human/` |
| Felsher signature (library v1.0) | `data/genesets_from_library_human/` |
| Menegollo bicluster forkscale | `data/menegollo_biclusters/` |
| Curated metabolic genes | `data/genesets_metabolic_human/` |

**The library's `outputs/gmt/human/` tree is mouse-derived and must not be
loaded** — it is mouse-native sets pushed through `mouse_to_human()`, gitignored
upstream and unpinned by the tag. The rejection is of that tree only; tracked
raw inputs carrying their own native human data are a different object. **Check
the sheet, not the repository.**

### Standing conventions

- mtDNA-encoded protein-coding genes (13, `MT-` prefix) sit in their own
  synthetic pathway and are never pooled with nuclear-encoded OXPHOS subunits.
- MitoCarta's `OXPHOS` umbrella includes assembly factors; `OXPHOS subunits` is
  the narrower set. They are different — pick deliberately.

## R coding rules

Full file: `docs/R_CODING_INSTRUCTIONS.md`. These five cause the most damage.

1. **Never `print(n = X)` after `head()`.** `head()` may coerce a tibble to a
   data.frame, so `n` is read as `na.print`. Use `head(X) %>% print()`.
2. **Always `dplyr::count()`**, never bare `count()` — namespace conflicts.
3. **ASCII-only strings in scripts.**
4. **Save the analysis frame**, not only the summary tables — one row per
   analysis unit, with the join key, every exposure and covariate that entered
   a model, every outcome and time variable, every stratifier used or audited,
   **and the covariates a downstream question might plausibly need.** `E41`
   saved 26 tables and no frame; `X01` was then unable to stratify METABRIC at
   all until `51465b6` added it, and separately unable to weigh confounding by
   indication because `grade`, `size`, `lymph_nodes_positive` and
   `age_at_diagnosis` had not been carried across. **Not** the expression
   matrix — the analysis-unit table, which is trivial on disk. Appending a
   frame to an already-verified object changes no existing element, so recorded
   digests still hold; say so and say that a re-source is needed.
5. **Save the FULL coefficient matrix, the FULL VIF vector and the FITTED MODEL
   OBJECTS** — never a selected row. Declaration `14.7`. `E40` and `E41`
   extracted only `co["OX", ]`, which discarded every other term's estimate and
   forced `E43` to refit to recover one. **This is the third retrieval forced by
   selective extraction** and it is cheap to prevent. **`E43` is the first
   script in the arm to comply.** One caveat for verification: a fitted model is
   **not digest-stable** across runs — `coxph` and `glm` carry `terms`,
   `formula`, `model` and `family`, each holding an environment — so save a
   `fit_digests` table of `coef()` and `vcov()` per fit and exclude the fits
   themselves from any digest comparison.

No `renv`; packages are installed system-wide.

## Workflow — "Option A" (do not deviate)

- Claude Code **writes and edits** the numbered pipeline scripts. It does **not
  run them.** The author sources them in Positron interactively.
- Infrastructure (git, snapshots, provenance READMEs, editing this file,
  planning and result notes) Claude Code may execute directly.
- Every numbered script ends with an `if (FALSE) { ... }` sandbox block —
  skipped by `source()`, run line-by-line in Positron for inspection.
- Commit per verified phase. Git is the safety net.
- When in doubt, ask.

## Project structure

```
scripts/       numbered R pipeline, E00-E57. `E31` lives on an unmerged
               branch, so it is the one gap in this tree. `X01_*` is OFF the
               numbering on purpose: grant-preparatory, manuscript-excluded
docs/          the aim, the plan, dated notes
docs/figures/  tracked copies of the figures a note relies on
data/          snapshots, each with a provenance README
functions/     shared utilities
results/       intermediate .rds (gitignored, generated at runtime)
outputs/       figures and tables (gitignored, generated at runtime)
```

From `E31` on, a script writes its key figure twice: to `outputs/figures/` as
usual, and to `docs/figures/`, which is tracked. `outputs/` does not survive a
fresh clone, so the tracked copy and the numbers written into the note are the
durable record.

`results/` and `outputs/` are regenerable. `data/from_validation/` is
regenerable by re-copying from the validation repo at the pinned SHA.

## Git discipline

- `main` is the trunk. Feature branches off `main` as needed.
- Read-only git ops are always fine. Stop-and-check before anything destructive;
  never force-push a shared branch.
- **Never write to `myc_human_validation` or `myc_mouse` from this repo.**
