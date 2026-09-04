---
date: 2026-09-04
status: DECLARED. Written and committed BEFORE the script was written or run.
posture: EXPLORATORY. Nothing in this repo is pre-registered, and this note is
         not a pre-registration. It is a declaration of direction and reading
         rules made in advance of the fit, so that the reading cannot be chosen
         after the numbers are seen.
relates-to:
  - docs/2026-09-03_human_arm_for_mouse_reconciliation.md (3.1, M2, V1, N2, N3)
  - docs/2026-09-04_e16_respiratory_rulers.md (the four rulers)
  - data/raw/depmap_README.md (the symlink, and what this repo cannot regenerate)
  - myc_human_validation @ d3ac60e, scripts/14_depmap_dependency.R (the template)
  - myc_human_validation @ d3ac60e, docs/2026-08-30_block_g_result_and_bcl2l11_declaration.md
decides:
  - The direction of the OX main effect on MCL1 and BCL2L1 dependency, per gene.
  - How each outcome is read, including the reverse direction, which is reported
    and is NOT scored as a pass.
  - That POLR2A is unavailable in 26Q1 and that POLR2B is a DECLARED substitute
    floor, named here before the fit.
  - B4-a to B4-e, and what voids what.
next-action: write scripts/E17_depmap_ox_dependency.R, then source it in Positron
---

# B4 - the OXPHOS main effect on guardian dependency. Declared before the fit.

**This is B4, an exploratory functional check.** It runs in
`myc_human_exploratory` and reads DepMap through the symlink created 2026-09-04.
It is not G3, does not depend on G3, and must not anticipate it.

**N3 applies to every sentence here.** These are transcript correlations and
gene-effect scores. Neither a tumour nor a cell line is described as "primed".

---

## 1. The question, and why it is a different estimand

Block G, in the closed validation study, fitted

```r
lm(Y ~ MYC * OX + PROLIF)          # and + lineage, pan-cancer
```

everywhere - breast, the 18-arm panel, pan-cancer, the BCL2L11 check, and the
drug arm - and extracted **only** `MYC:OX`. Verified 2026-09-04 by grep over the
forked script: every one of the seven `.tidy()` calls names `"MYC:OX"`, except
P3 which names `"MYC:OX:brst"`. **Nothing anywhere extracts `OX`.**

So the `OX` main effect was fitted in every one of those models and reported in
none of them. That coefficient is what B4 asks for.

**This is a different estimand, not a variant of a failed test.** The interaction
was null and stays null - Block G's own numbers are breast MCL1 +0.057 / +0.143
and BCL2L1 +0.097 / +0.042, and pan-cancer at n = 1,130 with CIs of about +/-
0.04. B4 does not revisit that and makes no interaction claim in either
direction. It asks whether OXPHOS has an effect **at all**.

**The prediction comes from the human tumour arm, independently of DepMap.**
Synthesis 3.1: in OXPHOS-high tumours `BCL2L1` transcript is up (+0.388) and
`MCL1` is down (-0.266), while MYC orders neither (+0.002, +0.022). **The
dependency prediction therefore carries no MYC term.** MYC stays in the model as
a covariate; it is not the exposure and nothing is stratified on it.

---

## 2. Direction, fixed now

**Chronos gene effect: 0 = no effect, -1 = median common essential. MORE
NEGATIVE = MORE ESSENTIAL.**

| gene | predicted `OX` coefficient | in words |
|---|---|---|
| **`BCL2L1`** | **negative** | OXPHOS-high lines are MORE BCL-XL-dependent |
| **`MCL1`** | **positive** | OXPHOS-high lines are LESS MCL1-dependent |

Nothing else has a declared direction. Every control gene is read as a floor or
as a specificity check, never as a hypothesis.

## 2.1 The two-sided problem, stated as script 14 states it

**A cell can need a protein more precisely because it has less of it.**
Transcript-up does not compel dependency-up, and the two quantities are not the
same thing.

The standard prior runs the other way - guardian expression tracks dependence on
it, which is the basis of the whole BH3-mimetic biomarker literature - and V1
assumes it. That assumption is what B4 tests. **The reverse is reported when it
occurs and is NOT scored as a pass.** If `BCL2L1` comes out positive and `MCL1`
negative, that is a clean, interesting, reportable result and it is a FAIL of
B4-a, written as such.

---

## 3. Primary and secondary, fixed now

| | scope | n | status |
|---|---|---|---|
| **primary** | Breast | **51** | the mouse is mammary and both human cohorts are breast |
| **secondary** | pan-cancer, lineage-adjusted | **1,140** | better powered, and the MYC/OXPHOS relationship need not be lineage-invariant |

**A breast null at n = 51 is uninformative, not evidence against.** That is fixed
here so it cannot be argued either way afterwards. The pan-cancer fit is
reported because saying "null at n = 51" without offering the powered comparison
would be a half-report.

Line counts confirmed against the files on disk 2026-09-04, and they reproduce
script 14's exactly: 96 breast models, **71** with expression, **53** with
CRISPR, **51** with both, **1,140** pan-cancer with both.

---

## 4. The models

```r
m_add <- lm(Y ~ MYC + OX + PROLIF)      # PRIMARY. V1's estimand
m_int <- lm(Y ~ MYC * OX + PROLIF)      # sensitivity only
```

`OX` is extracted from **both**, and both are reported. MYC, OX and PROLIF are
each z-scored within the fitted set, so in `m_int` the `OX` coefficient is the
effect **at mean MYC**, not at MYC = 0 in raw units.

**They should agree**, because Block G established `MYC:OX` is null on these
genes. A large disagreement is a diagnostic - it would say the interaction is
not negligible after all, in which case the additive `OX` is not interpretable
as a main effect - and it is reported as a diagnostic, **not as a result and not
as an interaction claim.** The line held here: `m_int` is fitted only as a
sensitivity, and only `OX` is extracted from it.

Pan-cancer adds `+ lineage`, as in script 14. Lineage adjustment is mandatory
there, not optional.

---

## 5. Rulers, genes, and one substitution declared in advance

**Four rulers**, all built or read in one CCLE run: `ox_gsva`, `ox_ppd`,
`ox_lvl`, `ox_rel`. E16 made `ox_rel` available and found it the most
infiltrate-robust ruler in TCGA; **cell lines have no infiltrate, so whether
that robustness was about infiltrate at all is itself worth knowing**, and it is
a reason to carry the ruler rather than a prediction about it.

`ox_rel` is built to E16's exact recipe, verified buildable in CCLE 2026-09-04:
numerator 87 of 89 (0.978, missing `ATP5MD` and `ATP5MPL`), denominator 1,036 of
1,047 (0.989), all 13 mtDNA-encoded genes present and in the **denominator**.

**Endpoints.** `DEP_PRIMARY` = `MCL1`, `BCL2L1`. `DEP_CONTROL` unchanged from
script 14 = `BCL2`, `BBC3`, `BCL2L11`, `BAX`, `BAK1`, `RPL3`, `POLR2A`. No new
endpoint genes.

### 5.1 POLR2A is not screened in 26Q1, and this is declared before the fit

Confirmed 2026-09-04 against `CRISPRGeneEffect.csv` through the symlink:
**`POLR2A` is absent.** Every other RNA Pol II subunit is present - `POLR2B`
through `POLR2M`, `POLR2J`, `POLR2J2` - but not `POLR2A`. The closed study hit
the same wall and recorded it: *"POLR2A is not screened in 26Q1, so the floor
rests on RPL3 alone."*

**B4-b is a hard requirement and it was written for two floor genes.** Leaving it
on one halves a criterion that voids B4-a on failure. So:

- `DEP_CONTROL` is **unchanged**. `POLR2A` is requested, found absent, and
  reported as not screened - exactly as script 14 does. It is not quietly
  dropped.
- **`POLR2B` is declared here, before the fit, as a substitute floor gene.** It
  is the second-largest RNA Pol II subunit, an unambiguous pan-essential, and
  fills the role `POLR2A` was named for. It is a FLOOR, not an endpoint, and no
  direction is predicted for it.
- **B4-b is reported both ways**: on `RPL3` alone, which is what the closed
  study could do, and on `RPL3` + `POLR2B`, which is the two-gene floor the
  criterion intends. **Both must pass.** If they disagree, B4-b fails.

---

## 6. Specificity is load-bearing here in a way it was not for the interaction

This is the substantive methodological difference between B4 and Block G, and it
is why the controls are not decoration.

**An interaction is robust to anything that shifts all dependencies together. A
main effect is not.** If OXPHOS-high lines simply grow differently - divide
faster, are more or less sensitive to any knockout - then `OX` picks that up on
every gene, and `MYC:OX` would not have. Three consequences, all fixed now:

1. **`RPL3` and `POLR2B` are a floor, not decoration.** Any signal there and the
   result is a growth-rate artefact. This voids B4-a; it does not qualify it.
2. **The full 18-arm panel runs for the primary genes.** If every mitochondrial
   arm predicts `BCL2L1` dependency, the effect is not respiratory.
3. **`PROLIF` stays in every model. Lineage adjustment is mandatory
   pan-cancer.** Neither is optional and neither is dropped to gain power.

## 6.1 The CCLE matched null does not exist, and that bounds the whole result

Script 14's `spec$no_null` still holds. **The TCGA expression-matched nulls are
not transferable to CCLE.** So the 18-arm panel is a **rank ordering, not a
calibrated comparison**, and:

> **A positive OXPHOS result is not reportable until the CCLE matched null is
> built.** The result note will report OXPHOS's rank position among the 18 arms
> and will attach no p-value to the arm comparison.

Building that null is out of scope for this session and becomes the next script
if B4-a passes.

---

## 7. PRISM, if present - and it is

Confirmed present through the symlink 2026-09-04. **PRISM Repurposing log2 fold
change: MORE NEGATIVE = MORE SENSITIVE**, the same direction as Chronos.

Compound coverage verified against the 24Q2 compound list, 2026-09-04:

| drug | target | present |
|---|---|---|
| `S63845`, `AMG-176`, `AZD5991` | MCL1 | **all three** |
| `venetoclax` | BCL2-selective | yes |
| `navitoclax` | BCL2 / BCL-XL / BCL-W | yes |
| `A-1331852`, `A-1155463` | **selective BCL-XL** | **ABSENT** |

**MCL1 is fully covered and is the counter-intuitive leg.** The field expects
MYC-driven breast cancer to be MCL1-dependent; the prediction here is `OX`
**positive** - OXPHOS-high lines LESS sensitive to MCL1 inhibition.

**BCL-XL is not covered.** Both selective inhibitors are absent. Navitoclax hits
BCL2, BCL-XL and BCL-W, and `BCL2` transcript also falls with OXPHOS (-0.117),
so navitoclax alone reads a mixture of two things moving the same way.

> **The interpretable quantity is navitoclax MINUS venetoclax on the same
> lines.** A navitoclax result is NEVER reported as a BCL-XL result.

Drug n is small (script 14 saw 26-33 breast lines with PRISM). A null there is
uninformative, and that is fixed here.

---

## 8. B4-a to B4-e, and how each outcome will be read

| | criterion | reading |
|---|---|---|
| **B4-a** | `BCL2L1` `OX` **negative** AND `MCL1` `OX` **positive**, breast, on **at least two of the four rulers** | **pass** both directions on >= 2 rulers; **partial** exactly one gene in its predicted direction on >= 2 rulers; **fail** otherwise, including the clean reverse |
| **B4-b** | neither direction present on `RPL3`, nor on `POLR2B` | **hard requirement. Failure voids B4-a**, which is then reported as a growth-rate artefact and not as a dependency result |
| **B4-c** | OXPHOS ranks above the other mitochondrial arms for the primary genes | reported as a **rank**, pending the null. No p-value is attached |
| **B4-d** | pan-cancer agrees in sign with breast | disagreement is **reported as lineage-dependence, not resolved**. Neither scope overrides the other |
| **B4-e** | the PRISM MCL1 leg agrees in sign with the CRISPR MCL1 result | a null at n ~ 30 is uninformative and is written as such |

**How "present" and "flat" are judged, fixed now so it is not chosen later.**
The criterion is the **sign of the estimate**, read across rulers and cohorts;
95% CIs are reported for every coefficient and it is stated whether each excludes
zero, but **no p-value threshold gates any verdict.** This follows the repo's
standing rule: consistency across rulers, genes and scopes is the evidence, and a
single cell of a grid is not. For B4-b specifically, "neither direction present"
means the floor genes' `OX` estimates are small relative to the primary genes'
and their CIs include zero; if a floor gene carries an estimate comparable to the
primary genes, B4-b fails whatever its p-value.

---

## 9. What would falsify B4

Written before the fit, as three named outcomes:

1. **Both primary genes flat** on both models and all four rulers, breast and
   pan-cancer. V1's dependency leg does not transfer to cell lines.
2. **The effect present on `RPL3` or `POLR2B`.** A growth-rate artefact. This is
   the most likely way a positive B4-a would be wrong, and it is why B4-b voids
   rather than qualifies.
3. **The effect as large on every other mitochondrial arm as on OXPHOS.** Not
   respiratory - a general mitochondrial-content effect.

## 10. What a positive still does not buy

- **Not a null overturned.** N2 is pre-registered and stands. B4 asks a
  different question of a different quantity in a different system.
- **Not an interaction.** No interaction claim is made in either direction.
  `m_int` is a sensitivity and only its `OX` term is read.
- **Not causation in tumours.** DepMap is a cell-line panel; CRISPR and PRISM are
  interventions on the gene, not on OXPHOS.
- **Not reportable as a specificity result** until the CCLE matched null exists.

## 11. Out of scope, named so it stays a decision

No MYC stratification. No Johnson-Neyman. No interaction claim. No new endpoint
genes - `POLR2B` is a declared floor substitute, not an endpoint. No CCLE matched
null this session. **Nothing written to `myc_human_validation`** - not results,
not tables, not a note; its `results/depmap_dependency.rds` is part of the frozen
record and is neither read nor rewritten.
