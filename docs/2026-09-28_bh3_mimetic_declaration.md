---
date: 2026-09-28
status: DECLARATION. Committed BEFORE the script was written and before any
        model was fitted. Nothing below is a result.
posture: EXPLORATORY and POST-HOC. Nothing here is pre-registered. This is a
         declaration of order, not a registration: it fixes the estimand
         hierarchy, the directions, the population and the stopping rules in
         advance so that the ranking cannot be re-decided after the estimates
         are seen. The ordering is auditable in the git history of the
         bh3-mimetic-oxphos branch, as the STATE freeze was in the validation
         repo.
branch: bh3-mimetic-oxphos. NOT to be merged to main without being asked.
relates-to:
  - docs/2026-09-04_b4_addendum.md (E18 - the configuration is ABSENT from
    cell lines. This is why a null here is uninformative.)
  - docs/2026-09-04_b4_result.md (E17 / B4 - the CRISPR arm, reversed)
  - docs/2026-09-07_e22_bbc3_direct_myc_effect.md (beta2, the estimand that has
    held)
  - docs/2026-09-07_e23_guardian_balance.md (beta1, the estimand that failed)
  - data/raw/depmap_README.md (DepMap inputs, reached by symlink)
---

# Declaration - BH3 mimetic sensitivity on the OXPHOS axis in cell lines

**N3 throughout.** These are drug-response summary statistics and transcript
scores. No cell line is described as "primed".

---

## 1. The claim this comes from, stated so it can be held to

The manuscript argues that a high respiratory state is lethal on a MYC
background, and that tumours which retain a high respiratory state have escaped
either by **lowering the oxidative driver** or by **blocking the apoptotic
effector**. The escape-by-blocking route implies a dependency: a tumour holding
a high respiratory state needs the guardian that permits it, and should
therefore be selectively vulnerable to BCL-XL antagonism.

That is a statement about tumours. **This analysis asks only whether it shows
up in cancer cell line pharmacogenomics**, which is a different system, and
section 2 is the reason that difference is not a detail.

## 2. E18 already measured the thing this analysis depends on, and it is absent

`E18` (`docs/2026-09-04_b4_addendum.md`) asked whether the tumour transcript
configuration reproduces in CCLE. **It does not.** `BCL2L1` runs `+0.35` to
`+0.43` with OXPHOS in both tumour cohorts and `-0.11` to `+0.01` in CCLE; the
`BCL2L1 - BBC3` gap is **positive in 8 of 8 tumour cells and negative in 8 of 8
CCLE cells**. The verdict was 0 of 4 rulers in breast and 0 of 4 pan-cancer,
against a rule fixed before the numbers. The cheapest explanation - that the
tumour association is a purity or infiltrate artefact - was eliminated
independently by `E16` check 4 before `E18` was written.

**Two consequences, and they are asymmetric. Both go in the result note
whatever it says.**

- **A null here is uninformative.** The transcript configuration the prediction
  is read off does not exist in these cells. A cell-line system that does not
  carry the exposure cannot be asked to carry the dependency, so a null is a
  statement about cell lines and not about the model.
- **A positive here is interesting precisely because of that.** A dependency
  that appears without the transcript configuration behind it is not explained
  by the configuration, and would have to be explained some other way.

**This asymmetry is declared here, before the estimates, so that it cannot be
read as a rescue afterwards.** It is not a licence to accept a weak positive.
The guardrails in section 8 are what a positive must clear, and they are fixed
here too.

## 3. The estimand hierarchy - FIXED BEFORE ANY MODEL IS FITTED

**This order is not re-decided after the estimates are seen.** If a lower-ranked
model fires and the primary does not, that is reported as such and is not
promoted.

| Rank | Model | Role |
|---|---|---|
| **PRIMARY** | `LN_IC50 ~ OX + MY + lineage` | the **OXPHOS coefficient at fixed MYC**. Mirrors `beta2` in `E22` / `E23` - the one estimand in this project that has held |
| SECONDARY | `LN_IC50 ~ OX + lineage`, fitted **separately within MYC tertile 1 and tertile 3** | the conditionality claim, without assuming multiplicativity on the log-IC50 scale |
| TERTIARY | `LN_IC50 ~ OX * MY + lineage` | the product term. **Reported, never relied on** - see section 8 |
| contrast | `LN_IC50 ~ MY + lineage` | the `beta1` arm that failed in `E23`; also the most proliferation-contaminated |
| companion | PRIMARY + `PROLIF_DISJOINT` | is it growth rate? |
| companion | PRIMARY + RB1 loss | declared confounder, section 6 |
| mediation only | PRIMARY + `z(BCL2L1)` | **never in the primary**, section 7 |

Every model in this table is fitted and every one is saved, including the ones
that do not fire. The results object carries all of them or it is incomplete.

## 4. Drugs

**Navitoclax (ABT-263) is primary.** It antagonises BCL-XL, BCL-2 and BCL-W.

**Venetoclax (ABT-199) is the specificity control and is predicted NULL.** It is
BCL-2 selective. **If venetoclax moves with navitoclax, the effect is not
BCL-XL** and the reading collapses whatever the navitoclax estimate does.

Also fitted if present in the release, and which were found is recorded in the
data README before the script is run: **ABT-737** (BCL-XL / BCL-2 / BCL-W, the
tool compound; may be in GDSC1 rather than GDSC2, so both are loaded),
**WEHI-539** and **A-1331852** (BCL-XL selective), **S63845** and **AZD5991**
(MCL1 selective), **sabutoclax** (pan-BCL-2-family).

Drug names are not stable across releases. The script greps
`unique(DRUG_NAME)` and records what it found; nothing is assumed about
spelling.

## 5. Population

**Solid lineages only.** Haematological lines are BCL2-dependent and would
dominate any BCL2-family signal; `Lymphoid` and `Myeloid` are excluded on
`OncotreeLineage`.

**Pan-solid with lineage as a factor is the primary.** Breast is the stratum of
interest and **will be underpowered at roughly 50 lines** - that is stated here,
before the estimate, so that a wide breast interval is not read as a result in
either direction.

**Stop rule.** If navitoclax is absent from the release, or fewer than ~300
solid lines survive the join, the script stops and reports rather than
substituting something else.

**Join rate.** GDSC `SANGER_MODEL_ID` to DepMap `Model.csv$SangerModelID`.
Asserted, and **the run stops below 90%**.

## 6. Declared confounder - RB1 loss

Varkaris et al., *Nat Commun* 2025 (doi 10.1038/s41467-025-60238-x) showed in
GDSC that **RB1 alteration is the strongest sensitiser to navitoclax among all
drugs tested**. MYC-high lines may be enriched for RB1 loss, so an unadjusted
OXPHOS or MYC effect could be RB1. It is declared here as a confounder and
enters as a companion model, not as a post-hoc check.

## 7. BCL2L1 expression is a MEDIATOR, not a confounder

If the model is right the path is **OXPHOS -> more BCL-XL -> more dependency**.
Adjusting for `BCL2L1` is therefore **over-adjustment**: it would remove the
effect the model predicts. It is fitted **once**, as a mediation check, and
reported separately. **It never enters the primary.**

## 8. Guardrails - fixed here, not after

1. **The estimand hierarchy is not reordered after seeing results.**
2. **Every model fitted is reported**, including the failures. The results
   object must contain all of them.
3. **A null on the product term alone means nothing.** That estimand has already
   failed in Block C, Block B, Block G, H4 and N2. **A sixth null is not a
   finding and is not presented as one.**
4. **Claim only what at least two readouts support** - `LN_IC50` and `AUC`
   minimally, PRISM as a third platform if it joins.
5. **No per-gene FDR**, consistent with the rest of the project.
6. **Nothing is written to `myc_human_validation`.** Gene set definitions may be
   read from it; nothing is written.
7. **GSVA is cohort-relative.** Every score here is computed on this line panel
   in one run and **is never compared numerically to a TCGA, SCAN-B or mouse
   value.** Only signs and orderings cross cohorts. This is CLAUDE.md and it is
   restated because this analysis sits beside tumour numbers in the manuscript.

## 9. MYC estimator

**`M_a` (`FELSHER__MITOSTRIP`) is primary**, as in the manuscript.

Cell lines carry copy number, so **MYC expression and 8q24 copy number are
computable here as sensitivities** - which partly defuses the D2 estimator
entanglement the neoadjuvant cohorts could not address. **All three are
reported.**

## 10. Sign conventions - all three readouts point the same way

| readout | direction |
|---|---|
| GDSC `LN_IC50` | **lower = more sensitive** |
| GDSC `AUC` | **lower = more sensitive** |
| PRISM log fold change | **more negative = more killing** |

**So a NEGATIVE coefficient on `OX` supports the prediction on every readout.**
This is stated here and repeated in the result note so that nobody reads the
forest plot backwards.

## 11. The prediction, written before the estimates

- **Navitoclax: OXPHOS coefficient NEGATIVE** on `LN_IC50`.
- **Larger in the MYC-high tertile** than in the MYC-low tertile.
- **Venetoclax: NULL.**

## 12. What the result will mean - written before it is seen

**SUPPORTED.** The primary OXPHOS coefficient is negative for navitoclax on
**both** `LN_IC50` and `AUC`; venetoclax is null; it survives the proliferation
companion; and it is larger in the MYC-high tertile than the MYC-low one.

**NOT SHOWN.** Any one of: the primary is null; venetoclax moves the same way;
it dies when proliferation is added; or the two MYC strata are
indistinguishable - **which would be a real but non-MYC-conditional effect, and
a weaker claim than the one the manuscript makes.**

Either way the result note states **n**, **set coverage**, and **that this is
exploratory and post-hoc**.

## 13. What this analysis is not

- It is **not** a test of the `MYC x OXPHOS` interaction on apoptotic priming.
  That is pre-registered, null, and closed at `d3ac60e`. CLAUDE.md trap 1: a
  correlation and an interaction are different questions and neither confirms
  the other.
- It **does not reopen** `myc_human_validation`.
- It is **not** a tumour result and does not qualify one. It is a cell-line
  pharmacogenomic readout, and `E18` is the reason that distinction is load
  bearing here rather than ceremonial.
