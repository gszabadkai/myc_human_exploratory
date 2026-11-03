---
date: 2026-11-03
status: >
  THE AUTHOR'S STOP 3 RULINGS, CONSOLIDATED. They replace the chat versions
  (items 1 to 6, amendments A to F), which never reached the repo. The display
  amendment (6f5bf4d) is already applied by E50 and is not repeated here.
  Record this file verbatim and commit it alone.
relates-to:
  - docs/2026-11-03_e49_result.md (4fd661b; Part 2 H2, Part 4 the text this replaces)
  - docs/2026-11-03_e50_result.md (E50's displays and their open choices)
  - docs/2026-11-03_E50_display_amendment.md (6f5bf4d)
  - docs/2026-11-03_signature_and_fig4g_brief.md (eb45a9b; its section 5 and 4.8 still govern)
---

# STOP 3 rulings, consolidated (the author, 2026-11-03)

## 1. The displays

1. **E50's displays are accepted, with the palette from Fig. 1g's source** (`#4A3525` / `#FAFAFA` / `#2A8A6D`, interpolated in CIELAB).
   - The measured values are not used.
   - They were extrapolated along a straight line in RGB, whereas Fig. 1g's scale is straight in CIELAB.
2. **E50's open choices are accepted:**
   - one ±0.6 scale shared by 4g, ED-e and ED-f;
   - the 1.6 mm gap before the Signature row;
   - white ink on 9 cells.
3. **The author does the following in Illustrator. No new script is needed.**
   - **Type sizes:** group names and the colour bar's title at 6 pt, as in Fig. 1g.
   - **Annotations on 4g:** the binding code, the signature signs, the § marks and any protein names.
   - **ED-g:** the "Proliferation minus OXPHOS" column is removed. It contains the coefficient that defines the signature and is read nowhere. It remains in E49's tables.
4. **Extended Data in the paper:** ED-a, ED-b, ED-e, ED-f, ED-g, ED-h, ED-i, ED-j and ED-k.
   - **ED-c and S-1 stay in the repo, out of the paper.** ED-c is superseded by ED-i and ED-j. S-1 is kept for referee responses on p53 and tumour composition.
   - **ED-d** stays superseded by ED-i.
   - **Every ED panel in the paper is cited:**
     - ED-a, ED-b, ED-h, ED-j and ED-k in the Results;
     - ED-g in the Methods;
     - ED-e, ED-f and ED-i from the 4g legend.

## 2. H2 (E49 result, Part 2)

1. **H2 is accepted.** The text states upfront that the signature was selected for tracking OXPHOS at fixed proliferation. That is the question the mouse posed, and the configuration the selection found is the finding.
2. **The signature's associations with OXPHOS and with proliferation are described as what the selection intended**, never as independent evidence.
3. **The independent evidence is:**
   - BBC3 and BCL2L1, named by the mouse before any selection;
   - MYC activity, which the selection never used;
   - the effectors, chosen by Fig. 1g's grouping.
4. **The H2 test is not declared.** It is recorded as the answer if a referee asks how much of the signature's proliferation coefficient the selection explains. If ever declared, it takes the next free number.

## 3. The text

This replaces E49's result note, Part 4. Measured against `docs/STYLE.md`.

**Abstract, human sentences:**

> Extending this relationship to 4,081 human breast tumours[8,9], a BCL-2 family signature that tracked OXPHOS at fixed proliferation showed BCL-xL rising more steeply than PUMA. MYC activity, which played no part in defining it, went with the same signature when measured with scores least correlated with proliferation.

**Mouse Fig. 4 paragraph, last sentence:**

> Consequently, established human tumours with high OXPHOS and MYC activity were predicted to carry raised pro-apoptotic signalling restrained by BCL-xL, a configuration distinct from that accompanying proliferation.

**Results, four paragraphs:**

> In tumours, however, OXPHOS rises with both MYC activity and proliferation (Extended Data Fig. Xa), so testing this prediction required OXPHOS to be separated from proliferation. We therefore searched the twelve BH3-only and anti-apoptotic genes of the BCL-2 family in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143) for those whose expression tracked OXPHOS at fixed proliferation. Six of the eleven mitochondrial members did, and together they formed a signature (Fig. 4g). BCL2L1, encoding BCL-xL, rose with the activators BBC3 and BID and the sensitisers BAD and BIK, whilst MCL1 fell.

> Its structure fitted the prediction, BCL2L1 rising more steeply with OXPHOS than BBC3, encoding PUMA, in both cohorts (Extended Data Fig. Xb). This is the reverse of the MYC-expressing gland, in which PUMA rose relative to BCL-xL with OXPHOS priority (Fig. 1i; Extended Data Fig. Xk), and resembles the mouse tumours that sustained high OXPHOS only alongside Bcl-xL (Fig. 3h,i). BAD rose more steeply than BCL2L1 (Fig. 4g), but BAD neutralises BCL-xL only when unphosphorylated, and survival signals phosphorylate it and hold it away from BCL-xL[Zha]. Its expression therefore does not measure its activity, unlike that of PUMA, which is controlled mainly through transcription[Nakano].

> MCL1 fell as BCL-xL rose, so BCL-xL was the only one of the two main guardians to rise alongside PUMA, which binds both[Chen]. MCL1's reported support of respiration in breast cancer stem cells depends on a mitochondrial function separable from its anti-apoptotic role[27]. This suggests that tumours with high OXPHOS rely on BCL-xL rather than MCL1, like persistent OXPHOS-high triple-negative breast cancer cells, which depend on BCL-xL and are insensitive to MCL1 inhibition[26].

> With OXPHOS, MYC activity and proliferation modelled together and intrinsic subtype adjusted for, the signature followed OXPHOS but not proliferation, as its selection intended (Fig. 4g). So did BBC3 and BCL2L1, the two genes named by the mouse before any selection (Extended Data Fig. Xh). MYC activity played no part in the selection, yet it went with the signature in the same direction as OXPHOS (Fig. 4g). This held on the two MYC scores least correlated with proliferation, and was weaker on scores that track proliferation closely (Extended Data Fig. Xa,j). The effectors BAX, BAK1 and BOK, outside the signature, also rose with OXPHOS and MYC activity in both cohorts (Fig. 4g).

**Clinical paragraph.** The opening sentence becomes:

> Having separated OXPHOS from proliferation on PUMA and BCL-xL, we next asked whether the two are also separable in their association with chemotherapy response.

Its other two sentences are unchanged: "Additional adjustment for MYC activity, and then for the mean expression of BCL2L1 and MCL1, did not materially alter the OXPHOS–response association (pooled odds ratios 0.87 and 0.88)." and "The transcriptional coupling of OXPHOS and proliferation therefore concealed a clinically relevant distinction."

**Closing paragraph, human sentences:**

> In human breast tumours, BCL-xL rose more steeply than PUMA as OXPHOS rose, and MYC activity went with the OXPHOS-associated signature. Proliferation, although correlated with OXPHOS, was not part of this configuration and went the other way on chemotherapy response.

**Methods notes.**
- Keep E49's three notes: the effectors outside the signature, the binding source, and the post hoc statement.
- Add the following, aligned to E47 declaration 4a.7, item 5, if its construction sentence differs:

> The signature comprised the BCL-2 family members listed in MitoCarta 3.0 whose partial Spearman correlation with OXPHOS at fixed proliferation, averaged over both cohorts, was at least 0.20 in magnitude. Each was entered with the sign of that correlation, and MYC activity played no part in the selection. Variants without BID, which met the threshold by the smallest margin, and with BMF, which met it but is not mitochondrial, gave similar associations (Extended Data Fig. Xg).

## 4. The 4g legend: edits to E50's draft

1. **After** "...every variable is rank-normalised within each cohort.", **add:** "Cells give the coefficient, rounded to one decimal; intervals are in Extended Data Fig. Xi."
2. **Replace** "MYC activity is scored on MYC_UP.V1_UP, a score with few proliferation genes; the other MYC activity scores are shown in Extended Data." **with:** "MYC activity is scored on MYC_UP.V1_UP, one of the two MYC activity scores least correlated with proliferation; the reference score is used in Extended Data Fig. Xf, and all four in Extended Data Fig. Xj."
3. **Replace** "The signature was selected for its association with OXPHOS at fixed proliferation (Methods)." **with:** "Bottom, the six-gene signature, selected for its association with OXPHOS at fixed proliferation (Methods)." This answers E50's point that the bottom row was described only by its label.
4. **Replace** "OXPHOS priority is shown in Extended Data." **with:** "OXPHOS priority is shown in Extended Data Fig. Xe."
5. **Add,** replacing E49's proposed additions (4.3): "Outside the signature, BCL2L11, BMF and BCL2 also fell with OXPHOS. At fixed OXPHOS and MYC activity, proliferation went with lower BAX and BOK and higher PMAIP1 and BCL2L11 in both cohorts; its association with BAK1 differed in sign between cohorts, and MCL1 was higher only in SCAN-B."
6. **Once the binding classes are entered (section 5),** delete "; blank, not yet verified".

## 5. The binding code

The author enters BID, BIK and BMF in `data/annotations/bh3_binding_code.csv`, from Chen et al. 2005's full binding data, with `verified_by` filled. They appear on 4g only through the author's Illustrator annotations.

## 6. Metrics

Measured on the four Results paragraphs:
- 16 sentences and 345 words;
- mean 21.6 words, median 20, range 13 to 30;
- 38% of sentences under 20 words.

The abstract's sentences are 25 and 22 words; the new Methods sentences are 30, 18 and 24. No sentence is under 12 words, and there are no em-dashes and no serial commas. Re-measure with `pack_text.py`.

## 7. What Claude Code does

1. **Record this file verbatim**, committed alone.
2. **Put section 3's sentences into** `docs/2026-11-03_e49_review_pack_src/pack_text.py`, the single source, and section 4's edits into the legend shown in the pack. Re-measure.
3. **Regenerate the E50 review pack.**
4. **Carry the text into** `docs/2026-11-03_e49_result.md`, Part 4, 4.1 and 4.3, stating that it supersedes the committed text.
5. **Update CLAUDE.md and the handoff,** including:
   - the Illustrator work (1.3);
   - the Extended Data list (1.4);
   - that the remaining display changes are the author's, outside the repo's scripts.

**No new display script is needed.**
