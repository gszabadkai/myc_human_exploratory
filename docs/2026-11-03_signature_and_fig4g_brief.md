---
date: 2026-11-03
status: >
  AUTHOR'S DECISION RECORD AND BRIEF. Written in the claude.ai Project chat
  with the author, from the E48 result note, the review pack, the provenance
  note and the Doc's Fig. 1. It restores the six-gene signature as the summary
  of Fig. 4g and specifies E49. E49 is NOT declared, written or run. Nothing
  here is a new result. Commit it as written, alone.
decisions: 2026-11-03, by the author (section 1)
supersedes: >
  - docs/2026-11-03_handoff.md sections 0 and 3 (the STOP-3 rulings on E48's
    layouts and on result note 4.1's text choices);
  - docs/2026-11-03_e48_result.md Part 4 as the paper's text (section 8 below
    replaces it);
  - in docs/2026-11-02_separation_line_brief.md: section 1.4's "the signature
    leaves the main text", section 6 (the 4g design) and section 8 (the text
    drafts). Its section 4 rules still govern, with the additions in section 5
    below.
relates-to:
  - docs/2026-11-03_e48_result.md (Part 2 holds all twelve genes' coefficients)
  - docs/2026-11-02_E48_separation_declaration.md (2941438; addendum 258a496)
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (addendum 4a.7: the signature's exact rule)
  - docs/2026-10-23_e47_result.md (MIXED on the reference score, which stands)
  - docs/2026-10-22_six_gene_composite_provenance.md (sections 4 and 5: why each gene is in or out)
  - docs/STYLE.md
next-action: >
  Claude Code checks structure, drafts the E49 declaration from section 6, and
  stops for the author's approval (see the session prompt).
---

# Fig. 4g around the signature, and the brief for E49

## 0. In brief

1. **The six-gene signature returns as Fig. 4g's summary.** Only its association with OXPHOS is circular. That number is shown and never read; everything else the signature carries is a finding.
2. **Fig. 4g becomes a heatmap in Fig. 1g's layout.**
   - Rows: the four groups of Fig. 1g (effectors, activators, sensitisers, anti-apoptotic).
   - Left of the heatmap: the BH3 binding code and the signature signs.
   - Columns: OXPHOS, MYC activity on two scores, and proliferation.
   - Beneath: a strip for the signature itself.
3. **PUMA against BCL-xL is the balance the text reads.** These are the two genes whose expression measures their activity. BAD rises most, but its activity is set by phosphorylation, so its transcript cannot be read as a load on BCL-xL. BID is treated the same way, because its activity is set by cleavage.
4. **MCL1 is answered without BAD.** PUMA binds both guardians. MCL1 falls as BCL-xL rises, so BCL-xL is the rising guardian able to hold PUMA. MCL1's link to OXPHOS in the literature is a non-apoptotic mitochondrial function, and OXPHOS-high persister cells depend on BCL-xL, not MCL1.
5. **E49 is light.** It recovers the signature's quantities from E47's saved fits and fits only the three effectors. It renders the displays and stops three times, as E48 did.

## 1. Decisions, the author, 2026-11-03

1. **The signature (S6) returns** as the main-figure summary, reversing the 2026-11-02 brief (1.4) and E48 (section 3, "not used"). The reasons:
   - the single-gene presentation was hard to read;
   - the circularity is confined to one number, which the new display shows grey and never uses as evidence.
2. **Fig. 4g is specified in section 7.** E48's three layouts retire from the paper and stay in the repo as E48's verified record.
3. **4g uses OXPHOS level** (GSVA); **priority (mitoPPS) goes to Extended Data.** This reverses "both measures every time" for this panel only, to match the clinical cohorts. The author may revisit it at STOP 3.
4. **The effectors BAX, BAK1 and BOK are fitted** and reported, so that 4g carries Fig. 1g's four groups. They stay **outside the signature**, whose definition is fixed (section 3).
5. **A human panel in Fig. 1i's format** goes to Extended Data: PUMA:BCL-xL against OXPHOS priority.
6. **E46 goes to Extended Data as built.** Its signature strip is now consistent with 4g, which closes open item 1 of the 2026-11-02 brief.
7. **Heatmap significance marks: render both versions** (section 6.8). The author chooses at STOP 3.
8. **Two claims from the chat are retracted** (section 4.8).
9. **"Nothing is discarded"** (E48 declaration, section 8) continues to apply.

## 2. Why the signature can return

**What is circular.** The six genes were chosen because they follow OXPHOS, so "the signature rises with OXPHOS" is true by selection. Its OXPHOS coefficient (+0.30 to +0.61, E47) is therefore drawn grey, labelled "defines the signature", and never read.

**What is not circular, and carries the messages:**

| quantity | value so far | status |
|---|---|---|
| composition: which genes follow OXPHOS, in which direction | E48 Part 2 | the finding the heatmap shows |
| proliferation's coefficient, at fixed OXPHOS and MYC | −0.22 to −0.38, 8 of 8 cells (E47 3.1) | RETRIEVED; **E49 declares it** |
| MYC activity's coefficient, at fixed OXPHOS | positive in 8 of 8 cells on the least-entangled score; MIXED on the reference score (E47) | declared; the MIXED verdict stands |
| BCL2L1 more steeply than BBC3 | E47 section 9 at `m1`; E48 C3 at `m3` | licensed |

**Raw against adjusted, for proliferation.**
- Unadjusted, the signature barely correlates with proliferation: −0.09 / −0.04 (reconciliation proposal 1.1, AD HOC).
- Its negative coefficient appears once OXPHOS and MYC activity, both correlated with proliferation, are held fixed.
- **Text wording:** "the signature did not follow it, and at fixed OXPHOS and MYC activity higher proliferation went with lower signature expression".
- E49 computes and displays the raw correlations, so that both halves of that sentence rest on declared output.

## 3. Why these genes

**The rule** (E47 declaration addendum 4a.7, item 5): the MitoCarta 3.0 members of the twelve whose mean partial correlation with OXPHOS, at fixed proliferation, over TCGA and SCAN-B, is at least 0.20 in size. Signs follow the direction of that association, and the rule never uses MYC.

**The six:** BAD, BIK, BBC3, BID and BCL2L1 enter positive; MCL1 enters negative.

**Why the other six are out** (provenance note, section 5):
- BMF is not mitochondrial;
- BIM's two cohorts disagree;
- BCL2 is oestrogen-driven, and was excluded before any fit;
- NOXA follows MYC activity and proliferation, not OXPHOS;
- BCL-w and A1 are flat.

**The threshold.** E47's variants S5 (without BID) and S7 (with BMF) give the same picture. E49 shows them in Extended Data.

**An observation, to be stated as one.** The selection used no binding information. Yet the genes that follow OXPHOS are BCL-xL with the activators PUMA and BID and the sensitisers BAD and BIK. MCL1 falls, and NOXA, its selective antagonist, does not rise.

## 4. BCL-xL rather than MCL1, and the BAD problem

**The sources were verified in PubMed on 2026-11-03.** Full references are in section 9.

### 4.1 BAD: what a reader sees first

A reader sees that BAD rises most and that BAD binds BCL-xL. The first reading is "BCL-xL is neutralised, so the tumour is closer to death". That is the textbook mechanism, so **the text must answer it in the sentence that reports BAD, or the next one.**

**The answer (Zha et al. 1996).** BAD neutralises BCL-xL only when it is unphosphorylated. Survival signals phosphorylate it, and phosphorylated BAD is held in the cytosol by 14-3-3 proteins, away from BCL-xL. So BAD's transcript does not measure active BAD.

### 4.2 PUMA: the gene whose expression measures its activity

PUMA was identified as a gene that p53 induces (Nakano & Vousden 2001), and it is controlled mainly through its transcription. **PUMA against BCL-xL is therefore the comparison that expression data can inform.** It is also the comparison the mouse made:
- in the MYC-expressing gland, PUMA:BCL-xL rose with OXPHOS priority (Fig. 1i, R² 0.73);
- the mouse tumours sustained high OXPHOS only alongside Bcl-xL (Fig. 3h,i).

### 4.3 The rule is applied evenly

BID's activity needs cleavage by caspase-8 (Li et al. 1998), so its transcript is treated as BAD's is. That does not help the argument: BID rises less steeply than BCL-xL. **BAD and BID carry a § in 4g.**

### 4.4 MCL1, without BAD

Chen et al. 2005 established the binding selectivities:
- PUMA and BIM bind every guardian;
- BAD binds BCL-2, BCL-xL and BCL-w, but not MCL1;
- NOXA binds only MCL1 and A1.

In the data, MCL1 falls and BCL-xL rises with OXPHOS, and PUMA, which binds both, rises less steeply than BCL-xL. **BCL-xL is therefore the only one of the two main guardians that rises alongside PUMA.**

### 4.5 MCL1's link to OXPHOS in the literature is a different function

- **Lee et al. 2017** (probably the Doc's ref. 27; check in Paperpile). MYC and MCL1 raise mitochondrial OXPHOS in chemoresistant TNBC stem cells. An MCL1 mutant that cannot enter mitochondria lost that effect but kept its anti-apoptotic function.
- **Perciavalle et al. 2012.** A matrix isoform of MCL1 supports respiration, whilst the outer-membrane isoform blocks apoptosis.

### 4.6 The functional dependence is on BCL-xL

**Andrzejewski et al. 2025** (the Doc's ref. 26) studied quiescent, OXPHOS-high TNBC cells that persist after doxorubicin.
- BCL-xL was the top survival dependency in a CRISPR screen.
- The cells were sensitive to BCL-xL inhibitors, and not to inhibitors of BCL-2 or MCL1.
- Their proliferating counterparts did not show this dependency.

### 4.7 The mouse, by direction only

- **Fig. 1i:** the human PUMA:BCL-xL slope against OXPHOS is the reverse of the gland's.
- **Fig. 3h,i:** Bcl-xL is what allowed high OXPHOS in the mouse tumours.
- **Proteomics:** the author's MYAZ proteomics note of 2026-10-23 (in the claude.ai Project, not in this repo) finds Mcl1 protein unchanged across all three OXPHOS constructs (51st to 80th percentile).

### 4.8 Retracted in the chat; do not use

1. **"Its [BAD's] rise increases reliance on BCL-xL rather than overriding it."** It is half true and does not answer the reader. If BAD outpaces BCL-xL, the spare BCL-xL shrinks.
2. **The BAD-based argument for BCL-xL over MCL1** ("the rising load includes BAD, which BCL-xL can absorb and MCL1 cannot"). It is replaced by 4.4.

### 4.9 Caveats and options

- **MCL1 protein** is short-lived and controlled after translation, so a fall in its mRNA is not a loss of MCL1 protection. The argument in 4.4 does not rely on MCL1 falling.
- **Optional, for the Discussion.** Phosphorylated BAD supports glucose-driven mitochondrial respiration (Danial et al. 2003, 2008). That was shown in liver and β-cells through glucokinase, so for breast it is a hypothesis only.
- **Not in E49: a protein check.** TCGA's protein arrays may include phospho-BAD, BCL-xL and MCL-1; check the panel. That would be a new data source and needs its own declaration.

## 5. Rules for interpreting and writing

These add to the 2026-11-02 brief's section 4, which still governs.

1. **The signature's OXPHOS association is never a finding.** The text says so once.
2. **Proliferation is worded at both levels**, raw and adjusted, as in section 2.
3. **The balance is PUMA against BCL-xL.** No sentence infers a net apoptotic load from the BAD, BID, BIM or MCL1 transcripts.
4. **BAD is reported with its explanation** in the same sentence or the next.
5. **Mouse and human are compared by direction only** (Fig. 1i; Fig. 3h,i). No coefficients sit side by side.
6. **Any pattern visible in 4g that the text's summary does not cover gets one descriptive clause**, in the text or the legend. Examples: MCL1's proliferation sign differs between cohorts, and BIM differs by cohort.
7. **Words that stay out:** the 2026-11-02 list (4.9), plus the two retractions in 4.8.

## 6. The E49 specification (to become the declaration)

### 6.1 Title and posture

**"E49: the OXPHOS-associated BCL-2 signature beside OXPHOS, MYC activity and proliferation; the three effectors; and Fig. 4g."**

EXPLORATORY and POST HOC.
- **Part A** recovers quantities from E47's verified saved fits (coefficients and `vcov`), on the E43, E44 and E48 precedents. It fits nothing.
- **Part B** fits the three effectors with E47's exact model.
- **The displays** are specified in 6.8.
- **There is no reading rule and no label.** Everything is reported and may be described. The declaration exists because E47 3.1 requires one before the signature's proliferation coefficient can carry a manuscript sentence.
- Confirm that E49 is still the next free script number.

### 6.2 Seen before writing: correct this list from your own record

**Seen:**
- the signature's proliferation coefficient (E47 3.1, as a range, and per cell if you printed it);
- its OXPHOS and MYC coefficients (E47 section 2 and the ladder);
- all twelve genes' coefficients (E48 Part 2);
- the signature's raw correlations with OXPHOS, MYC and proliferation (reconciliation proposal 1.1).

**Believed unseen:**
- the signature's within-model difference, and its interval;
- S5's and S7's proliferation coefficients;
- every effector quantity;
- the Fig. 1i-style display.

### 6.3 Readouts

- **S6, S5 and S7**, exactly as E47 built them. Their membership does not change.
- **The twelve genes**, re-read from E48's saved objects.
- **BAX, BAK1 and BOK**, fitted new (6.5).

### 6.4 Quantities

All are at `m3` unless stated. Each is computed on the reference and least-entangled MYC scores, on both OXPHOS measures, in both cohorts.

1. **S6, S5 and S7:** proliferation's coefficient, and the within-model difference (proliferation minus OXPHOS), from E47's `coef` and `vcov`, by E48's method.
2. **S6, S5 and S7:** MYC activity's coefficient at `m1` and `m3` on all four scores, re-read from E47. HALLMARK's `m3` cells are void.
3. **S6, S5 and S7:** OXPHOS's coefficient, displayed grey and never read.
4. **Raw Spearman correlations** of S6, S5, S7 and all fifteen genes with each OXPHOS measure, each MYC score and proliferation, per cohort, on E34's analysis sets. Descriptive.
5. **The effectors:** OXPHOS's, MYC's and proliferation's coefficients at `m1` and `m3`, and the within-model difference, as E48 computed them for the twelve.
6. **The Fig. 1i-style display.**
   - The added-variable relation of z(BBC3) − z(BCL2L1) to OXPHOS priority, at fixed reference MYC activity, proliferation and the spine, per cohort.
   - Draw decile means and the fitted slope.
   - **The slope must equal β_OX(BBC3) − β_OX(BCL2L1) from E47's `m3` priority fits, exactly.**
   - Cite C3's bootstrap interval for it.

### 6.5 The effector fits

- **Inputs:** the same expression matrices, normalisation and per-cohort Blom scoring as the twelve (E34's recipe), on E34's analysis sets.
- **The model:** `readout ~ OX + MYC + spine (+ PROLIF at m3)`, on all four MYC scores, both OXPHOS measures and both cohorts.
- **Gate:** the same code, refitting BBC3, must reproduce E47's stored coefficients and `vcov` exactly.
- **The effectors are never added to the signature,** whatever their association. The Methods says so.
- **Structure check first.** If an effector is missing from a matrix, record that and proceed without it.

### 6.6 The binding-code annotation

Create a small tracked table, `data/annotations/bh3_binding_code.csv`, with the columns `gene`, `binds_BCL2L1`, `binds_MCL1`, `source` and `verified_by`.

- **Verified from the abstract of Chen et al. 2005:**
  - BAD: binds BCL-xL, not MCL1;
  - PMAIP1 (NOXA): binds MCL1, not BCL-xL;
  - BBC3 (PUMA) and BCL2L11 (BIM): bind both.
- **BIK, BID and BMF: the author supplies their classes from the paper's full binding data.** Until then their cells are blank and nothing is drawn. **Never infer a class.**

### 6.7 Gates and assertions

As E48, and in addition:
- the S6, S5 and S7 re-read values equal E47's exactly;
- each difference recomputed from `coef` and `vcov` equals the stored one;
- the effector code reproduces BBC3 from E47 exactly;
- the twelve genes' values equal E48's saved objects exactly;
- the raw S6 correlations match the reconciliation proposal's values to 2 dp, as a sanity check;
- the Fig. 1i-style slope equals the E47 `m3` priority difference exactly;
- every drawn value equals its stored value, read back from `ggplot_build()`;
- colours come from the validated palette, and no text is below 5.5 pt;
- panels are 178 mm wide, and no statistical notation appears on any panel.

### 6.8 Displays

Their content is fixed now. Where versions are rendered, the author chooses after the run.

| display | content |
|---|---|
| **Fig. 4g** (main), two versions | the section 7 design: **(a)** unmarked; **(b)** a small dot on each cell whose 95% CI excludes zero |
| ED | the same heatmap on OXPHOS priority |
| ED | S5, S6 and S7 strips side by side |
| ED | raw correlations beside the adjusted coefficients, for the signature and the fifteen genes |
| ED | intervals for all fifteen genes: E48's ED-d extended with the effectors, the twelve asserted equal to E48's |
| ED | the Fig. 1i-style human panel, per cohort, on OXPHOS priority |

- **Kept from E48 unchanged:** ED-a (score correlations), ED-b (BCL2L1 against BBC3), ED-c (the MYC ladder) and S-1 (sensitivities).
- **Retired from the paper:** E48's layouts 1 to 3.
- **Kept, from E46:** the quadrant map, as built (decision 6).

### 6.9 Saved objects and digests

As E48: plotting frames, tables, assertions, and digests committed **before** the author's run. Fit objects are excluded from the digests.

## 7. The Fig. 4g design

The author's sketch, `2026-11-03_fig4g_fig1g_layout_SKETCH.png`, shows the intent. **It is not a source of values.**

**Size and type.** 178 mm wide, with no text below 5.5 pt. The legend is descriptive only.

**Rows**, in Fig. 1g's groups and order, with vertical group labels as in Fig. 1g:

| group | genes |
|---|---|
| Effectors | BAX, BAK1, BOK |
| Activators | BID, BBC3, BCL2L11 |
| Sensitisers | BAD, PMAIP1, BMF, BIK |
| Anti-apoptotic | MCL1, BCL2L2, BCL2A1, BCL2, BCL2L1 |

Protein names go in parentheses where they differ from the gene: PUMA, BIM, NOXA, BCL-w, A1 and BCL-xL.

**Left of the heatmap:**
- **"binds BCL-xL"** and **"binds MCL1":** filled or open circles, on BH3-only rows only, read from the annotation table and blank where unverified;
- **"signature sign":** + or −;
- **a § after BAD and BID.** The legend reads: "§ activity set mainly after translation (BAD: phosphorylation; BID: cleavage to tBID)".

**Heatmap columns:**
- OXPHOS | MYC activity, least-entangled | MYC activity, reference | Proliferation;
- each split into TCGA and SCAN-B, on OXPHOS level, at `m3`;
- OXPHOS and proliferation come from the reference-score fits, and each MYC column from its own score's fit. The legend says so.

**Colour.**
- Diverging blue to red, as in E46.
- Limits symmetric about zero, containing every value, rounded outward to 0.1.
- **Not Fig. 1g's brown to green.** The estimand differs: Fig. 1g shows z-scored group means, and 4g shows model coefficients.

**The signature strip** sits beneath, aligned to the eight columns. It shows S6's coefficient with its 95% CI in each column. The OXPHOS columns are grey and labelled "defines the signature".

## 8. Text drafts, to be finalised against E49

**Measured against `docs/STYLE.md`.** Brackets mark sentences that wait on E49. Figure citations and reference numbers are not counted.

**The last sentence of the mouse Fig. 4 paragraph:**

> Consequently, established human tumours with high OXPHOS and MYC activity were predicted to carry raised pro-apoptotic signalling restrained by BCL-xL, distinct from the consequences of proliferation.

**The signature paragraph:**

> In tumours, however, OXPHOS rises with both MYC activity and proliferation, so their associations can only be interpreted once the three are separated. In TCGA-BRCA (n = 938) and SCAN-B (n = 3,143) we therefore modelled each BCL-2 family gene on all three together, adjusting for intrinsic subtype (Fig. 4g). With the other two held fixed, six of the twelve BH3-only and anti-apoptotic genes followed OXPHOS and formed a signature (Fig. 4g). BCL2L1, encoding BCL-xL, rose with the activators BBC3 and BID and the sensitisers BAD and BIK, whilst MCL1 fell. Because the signature was defined by OXPHOS, that association is descriptive, whereas its associations with proliferation and MYC activity played no part in its definition. Although proliferation correlated with OXPHOS (ρ = 0.23–0.42), the signature did not follow it, and at fixed OXPHOS and MYC activity higher proliferation went with lower signature expression (Fig. 4g). MYC activity was associated with higher expression of the signature on the score least entangled with proliferation, and less so as the proliferation content of the score rose (Fig. 4g and Extended Data Fig. Xc). OXPHOS and MYC activity therefore went with the signature in one direction, and proliferation in the other. [One descriptive sentence on the effectors BAX, BAK1 and BOK (Fig. 4g).]

**The balance paragraph:**

> Within the signature, BCL2L1 rose more steeply with OXPHOS than BBC3, encoding PUMA, in both cohorts (Extended Data Fig. Xb). This is the reverse of the MYC-expressing gland, in which PUMA rose relative to BCL-xL with OXPHOS priority (Fig. 1i), and resembles the mouse tumours that sustained high OXPHOS only alongside Bcl-xL (Fig. 3h,i). BAD rose most steeply of all, but BAD neutralises BCL-xL only when unphosphorylated, and survival signals phosphorylate it and hold it away from BCL-xL[Zha]. Its expression therefore does not measure its activity, unlike that of PUMA, which is controlled mainly through transcription[Nakano].

**The MCL1 paragraph:**

> MCL1 fell as BCL-xL rose, so BCL-xL was the only one of the two main guardians to rise alongside PUMA, which binds both[Chen]. MCL1's reported support of respiration in breast cancer stem cells depends on a mitochondrial function separable from its anti-apoptotic role[27]. This suggests that tumours with high OXPHOS rely on BCL-xL rather than MCL1, like persistent OXPHOS-high triple-negative breast cancer cells, which depend on BCL-xL and are insensitive to MCL1 inhibition[26].

**Metrics for the three Results paragraphs**, without the bracketed sentence: 15 sentences, 342 words, mean 22.8, median 23, range 16 to 30, 27% under 20 words. No sentence is under 12 words, and there are no em-dashes and no serial commas.

**The abstract's human sentences:**

> Extending this relationship to 4,081 human breast tumours[8,9], OXPHOS was associated with a BCL-2 family signature in which BCL-xL rose more steeply than PUMA. MYC activity, once separated from proliferation, was associated with the same signature, whereas proliferation was associated with the opposite when OXPHOS was held fixed.

**The clinical paragraph.**
- Its opening becomes: "Having separated OXPHOS from proliferation on the BCL-2 family, we next asked whether the two are also separable in their association with chemotherapy response."
- Keep "Additional adjustment for MYC activity, and then for the mean expression of BCL2L1 and MCL1, did not materially alter the OXPHOS–response association (pooled odds ratios 0.87 and 0.88)."
- Keep "The transcriptional coupling of OXPHOS and proliferation therefore concealed a clinically relevant distinction."

**The closing paragraph, its human sentence:**

> In human breast tumours, OXPHOS and MYC activity were associated with a BCL-2 family signature in which BCL-xL outpaced PUMA, whereas proliferation ran the other way, on this signature as on chemotherapy response.

**Methods notes:**
- the signature's definition: E47 addendum 4a.7, item 5, with its construction sentence;
- the effectors were fitted for comparison with Fig. 1g and are outside the signature;
- the binding-code annotation and its source;
- the post hoc statement, as in E48 11.4.

## 9. References

**Verified in PubMed, 2026-11-03.** The author adds the new ones in Paperpile; the Doc's numbers then shift.

| reference | role |
|---|---|
| Chen L, et al. *Mol Cell* 2005;17:393–403. doi:10.1016/j.molcel.2004.12.030 | **new:** the binding code |
| Zha J, et al. *Cell* 1996;87:619–628. doi:10.1016/s0092-8674(00)81382-3 | **new:** BAD phosphorylation |
| Nakano K, Vousden KH. *Mol Cell* 2001;7:683–694. doi:10.1016/s1097-2765(01)00214-3 | **new:** PUMA's transcriptional control |
| Li H, et al. *Cell* 1998;94:491–501. doi:10.1016/s0092-8674(00)81590-1 | **new:** BID cleavage (the § note) |
| Lee KM, et al. *Cell Metab* 2017;26:633–647. doi:10.1016/j.cmet.2017.09.009 | probably ref. 27; check |
| Andrzejewski S, et al. *Cells* 2025;14:1557. doi:10.3390/cells14191557 | ref. 26 |
| Perciavalle RM, et al. *Nat Cell Biol* 2012;14:575–583. doi:10.1038/ncb2488 | optional, Discussion |
| Danial NN, et al. *Nature* 2003;424:952–956. doi:10.1038/nature01825 | optional, Discussion |
| Danial NN, et al. *Nat Med* 2008;14:144–153. doi:10.1038/nm1717 | optional, Discussion |

## 10. Superseded and carried

**Superseded** (front matter): the E48 handoff's STOP-3 rulings; E48 result Part 4 as the paper's text; the 2026-11-02 brief's 1.4 (on the signature), section 6 and section 8.

**Carried unchanged:**
- the 2026-11-02 brief's section 4;
- E47's MIXED verdict;
- E48's declared reading and every saved result;
- the co-amplification task (2026-11-02 brief, section 9);
- the REFEREE_QUESTIONS updates, now also a BAD entry and an MCL1 entry built from section 4;
- the Fig. 4 title;
- push only when the author asks.

## 11. Do not

- read any new E49 quantity before the declaration is committed: inspect structure only;
- refit E47 or E48, except the effector fits and their reproduction gate;
- change the signature's membership, or add the effectors to it;
- infer a binding class: annotate only verified ones;
- read the signature's OXPHOS coefficient as evidence;
- use the two retracted claims (4.8);
- set mouse and human coefficients side by side;
- choose between the two heatmap versions;
- download data, or edit the Google Doc;
- run numbered scripts except as redirected dry runs (Option A).
