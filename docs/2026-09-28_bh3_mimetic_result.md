---
date: 2026-09-28
status: RESULT. The declared PRIMARY estimand is NULL for navitoclax on every
        readout, both GDSC screens and PRISM. The analysis is not merged and
        is not carried into the manuscript.
posture: EXPLORATORY AND POST-HOC. Nothing here is pre-registered. The estimand
         hierarchy was fixed in docs/2026-09-28_bh3_mimetic_declaration.md and
         COMMITTED AT 552ba4f, before the script was written and before any
         model was fitted.
branch: bh3-mimetic-oxphos. NOT merged. Left in place deliberately - see s.11.
script: scripts/E31_bh3_mimetic_oxphos.R, committed unrun at 1d9a75b, fixed at
        df79eba. THE NUMBERS BELOW ARE FROM A DRY RUN with output paths
        redirected to the scratchpad; the repo was left untouched and
        results/bh3_oxphos_gdsc.rds is NOT yet on disk. The dry-run caveat is
        discharged when the author sources the script and the run reproduces
        this note object by object.
relates-to:
  - docs/2026-09-28_bh3_mimetic_declaration.md (the declaration - read first)
  - docs/2026-09-04_b4_addendum.md (E18 - why a null here is uninformative)
  - docs/2026-09-04_b4_result.md (E17 / B4 - the CRISPR arm)
  - data/raw/gdsc_README.md (the inputs, with URLs and SHA-256)
---

# BH3-mimetic sensitivity on the OXPHOS axis - result

**EXPLORATORY AND POST-HOC.** Nothing here is pre-registered. Every direction,
every model and the order they rank in were fixed in the declaration and
committed before the script existed.

**N3 throughout.** These are drug-response summary statistics and transcript
scores. No cell line is described as "primed".

**THE RESULTS OBJECT AND THE RAW DATA ARE GITIGNORED AND ARE NOT ON ORIGIN.**
That is why every estimate that matters is written into this note with n,
estimate, 95% CI and p. `results/bh3_oxphos_gdsc.rds` will not survive a fresh
clone. This note will.

---

## 0. The answer in six lines

| | |
|---|---|
| **The PRIMARY estimand** | **NULL.** The OXPHOS coefficient at fixed MYC for navitoclax includes zero in **all 5** screen-by-readout cells: GDSC2 `LN_IC50` **-0.008 [-0.094, +0.077]**, GDSC2 `AUC` **+0.017 [-0.067, +0.102]**, GDSC1 `LN_IC50` **+0.023**, GDSC1 `AUC` **+0.016**, PRISM **-0.059 [-0.176, +0.057]** |
| **and the specificity control** | **Venetoclax is null too**, so the control cannot discriminate - there is nothing for it to discriminate between |
| **Is it a power failure?** | **No, not on the GDSC readouts.** The same cohort, the same models and the same endpoints detect **RB1 loss on navitoclax at -0.894 [-1.130, -0.658], p = 4e-13**, replicating Varkaris et al. in all four GDSC cells. The pipeline finds a real navitoclax sensitiser and finds nothing for OXPHOS |
| **Is it therefore evidence against the model?** | **NO, and this is the whole point.** `E18` established that the BCL2L1-OXPHOS configuration **is absent from cell lines**. A system that does not carry the exposure cannot be asked to carry the dependency |
| **The cleanest BCL-XL probe** | **`WEHI-539` runs the WRONG WAY** - **+0.066** and **+0.074**, null but positive on both readouts, i.e. OXPHOS-high lines nominally *less* sensitive |
| **Where structure does exist** | It is **MCL1-directed, not BCL-XL-directed**: `AZD5991` negative in 3 of 3 cells, `sabutoclax` 2 of 2, `S63845` 1 of 1. **Post-hoc, not declared, and recorded in s.7 as an observation only** |

**The verdict the script printed, on the declaration's own branch logic:**

> **NOT SHOWN** - the primary is null on at least one of the two GDSC readouts.
> **AND THIS IS THE UNINFORMATIVE BRANCH.**

---

## 1. How this relates to E18, stated before the estimates and repeated here

`E18` (`docs/2026-09-04_b4_addendum.md`) asked whether the tumour transcript
configuration reproduces in cell lines. **It does not.**

| | TCGA | SCAN-B | CCLE breast | CCLE pan |
|---|---|---|---|---|
| `BCL2L1` with OXPHOS | **+0.379** | **+0.327** | **-0.045** | **-0.062** |
| gap `BCL2L1 - BBC3` | positive **8 of 8** | | negative **8 of 8** | |

Verdict 0 of 4 rulers in breast and 0 of 4 pan-cancer, against a rule fixed
before the numbers. The cheapest escape - that the tumour association is a
purity or infiltrate artefact - was closed by `E16` check 4 *before* `E18` was
written: `BCL2L1` on `ox_gsva` moves **+0.419 -> +0.406** on adjustment.

**The consequence was declared in advance and is asymmetric.**

- **A null here is uninformative.** The dependency is read off a configuration
  that does not exist in these cells.
- **A positive here would have been interesting**, because it would not have
  been explained by that configuration.

**We are on the first branch.** This result therefore **does not bear on the
manuscript's claim** and must not be cited as evidence against it. What it does
close is a specific, cheap hope: that cell line pharmacogenomics would supply
easy support for the escape-by-blocking route. **It does not, and now that is
measured rather than assumed.**

---

## 2. What was measured

| | |
|---|---|
| **Analysis cohort** | **579 solid cell lines**, 18 lineages, all screened in GDSC and carrying DepMap 26Q1 expression |
| Join `SANGER_MODEL_ID` -> `SangerModelID` | **976 of 978, 99.8%** - declared floor was 90% |
| Solid lines after excluding `Lymphoid` and `Myeloid` | 805 |
| **Attrition, reported not absorbed** | of those 805, **601 carry DepMap expression** (lost 204); 7 lineages below 10 lines dropped (`Adrenal Gland`, `Biliary Tract`, `Fibroblast`, `Pleura`, `Prostate`, `Testis`, `Vulva/Vagina`), leaving **579** |
| **Breast** | **47 lines** - the stratum of interest, **declared underpowered before it was fitted**, and it is |
| MYC tertiles on `M_a` | 193 / 193 / 193 |
| RB1 loss | **77 of 579** (18 by `Deletion`, 60 by truncating mutation) |
| Models fitted | **798**, of which **0** failed to estimate |

### Set coverage, as script 13 reported it for the neoadjuvant cohorts

| set | n_defined | n_present | frac |
|---|---|---|---|
| `MYC` (`FELSHER__MITOSTRIP`, M_a) | 61 | **60** | 0.984 |
| `PROLIF_DISJOINT` | 318 | **317** | 0.997 |
| **`OXPHOS subunits`** | 89 | **87** | **0.978** |

The narrow subunit set, **not** the umbrella, which carries assembly factors.

### Which BH3 mimetics were found

Checked before the script was written. **Navitoclax and venetoclax are in both
GDSC screens**, so the specificity contrast does not rest on one. `WEHI-539` is
in GDSC2 and is a **second, cleanly BCL-XL-selective probe** - not a second look
at navitoclax, which is a triple inhibitor. **`A-1331852` and `A-1155463` are
absent from all three platforms**, recorded so their absence is not misread.

Two naming traps, both real: **GDSC1 spells venetoclax `Venotoclax`**, and
**GDSC2 spells ABT-737 `ABT737`** and omits it from the compound list entirely.
The script greps the loaded tables and records what it matched.

---

## 3. THE PRIMARY. `LN_IC50 ~ OX + MY + lineage`, the OXPHOS coefficient

`ox_gsva`, `M_a`, all 579 lines, lineage as a factor. Coefficients are in SD of
response per SD of OXPHOS - a linear rescaling that leaves **every p value and
every sign unchanged** and exists only so `AUC`, bounded in `[0, 1]`, can sit
beside `LN_IC50`, which spans about ten natural-log units.

**SIGN CONVENTION. `LN_IC50` lower = more sensitive. `AUC` lower = more
sensitive. PRISM LFC more negative = more killing. A NEGATIVE COEFFICIENT
SUPPORTS THE PREDICTION ON ALL THREE.** Do not read the forest backwards.

| drug | screen | readout | n | estimate | 95% CI | p |
|---|---|---|---|---|---|---|
| **navitoclax** | GDSC2 | `LN_IC50` | 576 | **-0.008** | **[-0.094, +0.077]** | 0.85 |
| **navitoclax** | GDSC2 | `AUC` | 576 | **+0.017** | **[-0.067, +0.102]** | 0.69 |
| **navitoclax** | GDSC1 | `LN_IC50` | 558 | **+0.023** | [-0.065, +0.112] | 0.60 |
| **navitoclax** | GDSC1 | `AUC` | 558 | **+0.016** | [-0.072, +0.103] | 0.72 |
| **navitoclax** | PRISM | `LFC` | 376 | **-0.059** | [-0.176, +0.057] | 0.32 |
| **venetoclax** | GDSC2 | `LN_IC50` | 571 | **-0.003** | [-0.080, +0.074] | 0.94 |
| **venetoclax** | GDSC2 | `AUC` | 571 | **+0.015** | [-0.069, +0.099] | 0.72 |
| **venetoclax** | GDSC1 | `LN_IC50` | 548 | -0.007 | [-0.097, +0.083] | 0.89 |
| **venetoclax** | GDSC1 | `AUC` | 548 | -0.014 | [-0.104, +0.077] | 0.76 |
| **venetoclax** | PRISM | `LFC` | 366 | +0.083 | [-0.035, +0.202] | 0.17 |
| **wehi-539** | GDSC2 | `LN_IC50` | 574 | **+0.066** | [-0.021, +0.153] | 0.13 |
| **wehi-539** | GDSC2 | `AUC` | 574 | **+0.074** | [-0.014, +0.162] | 0.097 |
| abt-737 | GDSC2 | `LN_IC50` | 570 | +0.002 | [-0.083, +0.087] | 0.97 |
| abt-737 | GDSC2 | `AUC` | 570 | -0.004 | [-0.088, +0.079] | 0.92 |
| abt-737 | PRISM | `LFC` | 377 | -0.091 | [-0.207, +0.025] | 0.12 |
| azd5991 | GDSC2 | `LN_IC50` | 449 | -0.093 | [-0.186, +0.000] | 0.050 |
| azd5991 | GDSC2 | `AUC` | 449 | -0.072 | [-0.169, +0.025] | 0.14 |
| azd5991 | PRISM | `LFC` | 440 | -0.095 | [-0.197, +0.007] | 0.068 |
| sabutoclax | GDSC2 | `LN_IC50` | 532 | -0.078 | [-0.165, +0.008] | 0.076 |
| sabutoclax | GDSC2 | `AUC` | 532 | -0.072 | [-0.159, +0.015] | 0.10 |
| s63845 | PRISM | `LFC` | 440 | -0.055 | [-0.158, +0.048] | 0.29 |
| obatoclax | GDSC1 | `LN_IC50` | 531 | **+0.100** | **[+0.013, +0.186]** | 0.025 |
| obatoclax | GDSC1 | `AUC` | 531 | **+0.101** | **[+0.015, +0.188]** | 0.022 |
| obatoclax | GDSC2 | `LN_IC50` | 490 | +0.035 | [-0.053, +0.124] | 0.43 |
| obatoclax | GDSC2 | `AUC` | 490 | +0.058 | [-0.030, +0.147] | 0.20 |
| obatoclax | PRISM | `LFC` | 371 | -0.010 | [-0.127, +0.106] | 0.86 |

**Not one navitoclax cell is negative with an interval excluding zero. Nor one
venetoclax cell. Nor one `WEHI-539` cell.** The only primary cells clearing zero
belong to obatoclax, in GDSC1 only, in the **wrong** direction, and not
reproduced in GDSC2 or PRISM - which is what a single-screen artefact looks
like.

### 3.1 How small an effect this rules out

The navitoclax standard errors are **0.043 to 0.045** on the GDSC cells, giving
interval half-widths of **0.085 to 0.089**. So the data exclude an OXPHOS effect
larger than **about 0.09 SD per SD** in either direction.

**Set that beside the positive control in the same models: RB1 loss moves
navitoclax `LN_IC50` by -0.894 SD.** The OXPHOS interval therefore excludes
anything above roughly **one tenth** of a known navitoclax sensitiser's effect.
**This is a precise null on the GDSC readouts, not an underpowered one.**

---

## 4. THE POSITIVE CONTROL - and why it changes what the null is worth

Read from the `comp_rb1` fit that already computed it: the same model, a
different coefficient. **Negative = RB1-null lines more sensitive.**

| drug | screen | readout | n | RB1 estimate | 95% CI | p |
|---|---|---|---|---|---|---|
| **navitoclax** | GDSC2 | `LN_IC50` | 576 | **-0.894** | **[-1.130, -0.658]** | **4.0e-13** |
| **navitoclax** | GDSC2 | `AUC` | 576 | **-0.774** | **[-1.011, -0.537]** | **2.9e-10** |
| **navitoclax** | GDSC1 | `LN_IC50` | 558 | **-0.679** | **[-0.939, -0.419]** | **4.0e-07** |
| **navitoclax** | GDSC1 | `AUC` | 558 | **-0.629** | **[-0.886, -0.372]** | **2.0e-06** |
| navitoclax | PRISM | `LFC` | 376 | -0.063 | [-0.431, +0.305] | 0.74 |
| wehi-539 | GDSC2 | `LN_IC50` | 574 | -0.677 | [-0.925, -0.429] | 1.2e-07 |
| abt-737 | GDSC2 | `LN_IC50` | 570 | -0.664 | [-0.906, -0.422] | 1.1e-07 |
| venetoclax | GDSC2 | `LN_IC50` | 571 | -0.493 | [-0.714, -0.272] | 1.4e-05 |

**This replicates Varkaris et al., Nat Commun 2025, in all four GDSC cells and
both screens**, and it does so on the BCL-XL-containing compounds
(navitoclax, `WEHI-539`, ABT-737) more strongly than on BCL-2-selective
venetoclax, which is the expected pattern.

**Three things follow.**

1. **The cohort, the models and the endpoints can detect a genuine navitoclax
   sensitiser**, at effect sizes an order of magnitude above the OXPHOS
   interval. The GDSC null is about OXPHOS.
2. **PRISM is the weaker platform.** RB1 is null there too, on a single-dose
   log fold change. **A PRISM null carries much less than a GDSC null**, and
   nothing in this note leans on one.
3. **The declared RB1 confounder confounds nothing here.** Adding it moves the
   navitoclax primary from **-0.008 to -0.006** (GDSC2 `LN_IC50`). It is a
   strong predictor of the outcome and is simply not entangled with OXPHOS in
   these lines.

---

## 5. Every other model in the hierarchy, navitoclax, `ox_gsva`, `M_a`

**Reported in the declared order, including the ones that do not fire.** The
hierarchy is not reordered after the estimates.

| rank | model | screen | readout | n | estimate | 95% CI |
|---|---|---|---|---|---|---|
| **PRIMARY** | primary | GDSC2 | `LN_IC50` | 576 | **-0.008** | [-0.094, +0.077] |
| SECONDARY | secondary, MYC T1 low | GDSC2 | `LN_IC50` | 191 | +0.024 | [-0.122, +0.171] |
| SECONDARY | secondary, MYC T3 high | GDSC2 | `LN_IC50` | 193 | +0.060 | [-0.101, +0.220] |
| TERTIARY | tertiary, `OX:MY` | GDSC2 | `LN_IC50` | 576 | +0.012 | [-0.070, +0.094] |
| contrast | contrast, `MY` term | GDSC2 | `LN_IC50` | 576 | +0.041 | [-0.041, +0.122] |
| companion | + `PROLIF_DISJOINT` | GDSC2 | `LN_IC50` | 576 | -0.002 | [-0.085, +0.081] |
| companion | + RB1 | GDSC2 | `LN_IC50` | 576 | -0.006 | [-0.087, +0.076] |
| mediation only | + `z(BCL2L1)` | GDSC2 | `LN_IC50` | 576 | -0.009 | [-0.095, +0.076] |
| focus stratum | Breast only | GDSC2 | `LN_IC50` | **47** | +0.058 | [-0.252, +0.368] |

**The conditionality claim fails in the direction as well as the magnitude.**
The declaration predicted the effect would be **larger in the MYC-high tertile**.
It is **+0.060 in T3 against +0.024 in T1** - both null, both positive, and the
high-MYC stratum is nominally the *less* sensitive one. On GDSC1 the same
pattern (+0.083 against +0.073); on PRISM it reverses (+0.008 against -0.083).
**No conditionality, and no stable direction to the non-conditionality.**

**The product term is null, and that means nothing.** Declaration guardrail 3:
that estimand has already failed in Block C, Block B, Block G, H4 and N2. **A
sixth null is not a finding and is not presented as one.**

**The mediation check moves nothing** (-0.008 to -0.009), which is what the
absence of a mediator looks like: `E18` already showed `BCL2L1` does not track
OXPHOS in these cells, so there is no path to mediate.

**The breast stratum is uninformative at n = 47**, exactly as declared. Its
interval, `[-0.252, +0.368]`, is four times the width of the pan-solid one. It
is reported, not read.

---

## 6. The ruler sensitivity - and it is a textbook trap 5

CLAUDE.md trap 5: the four instruments disagree by a lot, so a result that
lives on one of them is a **ruler result**. The primary model on all four:

| drug | cell | `ox_gsva` | `ox_ppd` | `ox_lvl` | `ox_rel` |
|---|---|---|---|---|---|
| navitoclax | GDSC2 `LN_IC50` | -0.008 | -0.063 | **-0.099** * | -0.006 |
| navitoclax | GDSC2 `AUC` | +0.017 | -0.063 | -0.079 | +0.011 |
| navitoclax | GDSC1 `LN_IC50` | +0.023 | **-0.096** * | -0.028 | -0.001 |
| navitoclax | GDSC1 `AUC` | +0.016 | **-0.105** * | -0.042 | -0.008 |
| navitoclax | PRISM `LFC` | -0.059 | +0.025 | -0.058 | -0.025 |
| venetoclax | GDSC1 `LN_IC50` | -0.007 | **-0.089** * | -0.076 | +0.006 |
| wehi-539 | GDSC2 `LN_IC50` | +0.066 | +0.013 | +0.043 | +0.019 |

`*` = 95% CI excludes zero.

**Read this as a warning, not as support.** Navitoclax clears zero on
**`ox_lvl` in GDSC2 and on `ox_ppd` in GDSC1** - **different rulers in different
screens**, with the primary ruler firing in neither, and `ox_rel` firing
nowhere. **That is the signature of a ruler artefact, not of an effect.** And
`ox_ppd` fires for **venetoclax** in the same screen, which on the declaration's
own logic would say the effect is not BCL-XL anyway.

Had the declaration named `ox_lvl`, this note would have opened with a
significant negative navitoclax coefficient. **It named `ox_gsva`, before the
data, which is the entire reason the ruler panel is a sensitivity here and not
a menu.**

---

## 7. What DOES show structure - MCL1, and it is post-hoc

**Not declared, not predicted, and recorded as an observation.** Reported
because the declaration says every model is reported, and because it has a shape
worth writing down rather than a p-value worth quoting.

On the primary model, the **MCL1-directed** compounds are negative in **every**
cell, while the BCL-XL-directed ones are not:

| drug | target | cells | all negative? | range |
|---|---|---|---|---|
| `AZD5991` | MCL1 | 3 | **yes, 3 of 3** | -0.095 to -0.072 |
| `sabutoclax` | pan, incl. MCL1 | 2 | **yes, 2 of 2** | -0.078 to -0.072 |
| `S63845` | MCL1 | 1 | **yes, 1 of 1** | -0.055 |
| **navitoclax** | **BCL-XL/2/W** | 5 | **no, 2 of 5** | -0.059 to +0.023 |
| **`WEHI-539`** | **BCL-XL** | 2 | **no, 0 of 2** | +0.066 to +0.074 |
| `venetoclax` | BCL-2 | 5 | no, 3 of 5 | -0.014 to +0.083 |

**On the primary ruler and primary estimator not one MCL1 cell clears zero**
(`AZD5991` GDSC2 `LN_IC50` p = 0.050, PRISM p = 0.068). They clear zero on
`myc_expr` and `myc_cn8q24` and in the MYC-high tertile. **So this is a
consistent sign across six cells and three platforms, not a demonstrated
effect**, and it is written that way.

**It is nonetheless coherent with `E18`, which is why it is recorded.** `E18`
found that **`MCL1` is one of the four genes whose OXPHOS configuration DOES
transfer to cell lines** - negative on 4 of 4 rulers in CCLE breast, as in both
tumour cohorts - while `BCL2L1` is the gene that does not. **The compound class
whose transcript configuration survives into culture is the compound class that
shows a gradient here.** That is a consistency worth one sentence and no more.

**It is not support for the manuscript.** The manuscript's claim is about
**BCL-XL**. On BCL-XL the cleanest probe runs the other way.

---

## 8. The MYC estimators - a standalone observation

Pairwise Spearman on the 579 lines:

| | `M_a` | MYC expression | 8q24 copy number |
|---|---|---|---|
| **`M_a`** | 1.000 | **+0.360** | **-0.037** |
| **MYC expression** | +0.360 | 1.000 | +0.231 |
| **8q24 copy number** | **-0.037** | +0.231 | 1.000 |

**MYC activity and 8q24 copy number are uncorrelated in these lines
(`rho = -0.037`).** Copy number tracks MYC *transcript* weakly (+0.231) and MYC
*activity* not at all.

This is CLAUDE.md trap 4 measured in a third system: MYC mRNA is not MYC
activity, and amplification is not either. The declaration hoped copy number
would "partly defuse the D2 estimator entanglement the neoadjuvant cohorts could
not address". **It does not - it is a third quantity, not a cleaner version of
the first.** Worth carrying forward; it is cheap evidence against ever using
8q24 amplification as a MYC-activity proxy.

All three estimators were fitted on the primary. They **do not change the
navitoclax conclusion** (GDSC2 `LN_IC50`: -0.008, -0.003, +0.002 across `M_a`,
expression and copy number).

---

## 9. The grid, stated plainly

**155 of 798 fitted cells have an interval excluding zero.** CLAUDE.md: multiple
comparisons are the default state here, and a p-value plucked from the grid is
not a result.

| rank | cells | clearing zero |
|---|---|---|
| PRIMARY | 156 | 31 |
| SECONDARY | 156 | 22 |
| TERTIARY | 78 | 4 |
| contrast | 78 | 27 |
| companion | 156 | 25 |
| mediation only | 78 | 7 |
| **positive control** | **78** | **39** |
| focus stratum (breast) | 18 | **0** |

These cells are heavily nested - the same lines, the same drugs, models that
contain one another - so the count is **not** interpretable as a false discovery
rate. It is here so that the 31 firing PRIMARY cells are read as what they are:
**mostly the non-declared drugs on non-declared rulers and estimators**, with
**navitoclax on the declared ruler and estimator clearing zero in none of them**.

**The positive control firing in 39 of 78 is the one count that should be high**,
and it is the reason section 4 can say the GDSC null is precise rather than
empty.

---

## 10. What this result is, and what it is not

**IT IS.** A well-powered null on the GDSC readouts. The declared primary
estimand is zero to within about a tenth of a known navitoclax sensitiser's
effect, in 579 solid lines, on two independent GDSC screens, with the
specificity control equally null, surviving proliferation and RB1 adjustment,
with no MYC conditionality in either direction, and with the cleanest
BCL-XL-selective probe trending the wrong way.

**IT IS NOT evidence against the manuscript.** `E18` measured the exposure
configuration in cell lines and found it absent. This was declared before the
estimates, in writing, at `552ba4f`. **The null is a statement about cell lines.**

**IT IS NOT a test of the `MYC x OXPHOS` interaction on apoptotic priming.**
That is pre-registered, null, and closed at `d3ac60e`. CLAUDE.md trap 1.

**IT DOES NOT reopen `myc_human_validation`.** Nothing was written there.

**IT IS NOT a tumour result** and it qualifies none.

### What would change the answer

1. **A cell line system that carries the exposure.** The blocking argument
   predicts a dependency *given* the configuration. Nothing here had it. A
   panel selected or engineered for high `BCL2L1` on a high respiratory state
   would be a test; an unselected panel is not.
2. **A combination readout.** `MCL1` routinely compensates for `BCL2L1` loss -
   the same limitation recorded in the `B4` addendum section 4. A navitoclax
   response conditioned on `MCL1` status, or a co-dependency, is the readout a
   single-agent screen cannot give, and section 7 is a reason to think it might
   be where the structure is.
3. **The MCL1 observation tested properly**, on a declared estimand, with its
   own control, rather than read off this grid.
4. **A tumour-side pharmacogenomic readout.** There is none, which is why this
   was attempted in cell lines at all.

---

## 11. The branch, and why it stays

**`bh3-mimetic-oxphos` is NOT merged and is NOT deleted.** The manuscript is
close to submission and this analysis carries nothing into it.

**It stays on origin so that it is discoverable.** An exploratory analysis that
was run and abandoned should remain visible, so that nobody runs it a second
time without knowing it has been run once. **This note is the record of what was
found, and the commit sequence is the record of the order it was found in:**

| commit | what |
|---|---|
| `552ba4f` | the declaration, **alone**, before any data or code |
| `894953f` | the data README, with URL, size and SHA-256 for every input |
| `1d9a75b` | the script, **before it had been run** |
| `df79eba` | two dry-run defects fixed, and the RB1 positive control added |
| this one | the result note and the figure |

---

## 12. Where the numbers live

| what | where |
|---|---|
| every model fitted, with n, estimate, SE, 95% CI and p | `results/bh3_oxphos_gdsc.rds`, `x$coefficients` - **gitignored** |
| the same, as a table | `outputs/tables/E31_bh3_coefficients.csv` - **gitignored** |
| set coverage | `x$coverage`, `outputs/tables/E31_bh3_coverage.csv` |
| the figure | `outputs/figures/E31_bh3_forest.png` (gitignored) **and `docs/figures/2026-09-28_E31_bh3_forest.png` (tracked)** |
| which drugs matched, and under which spelling | `x$drug_map` |
| the inputs, with checksums | `data/raw/gdsc_README.md`, `data/raw/depmap_README.md` |

**Because the object is gitignored, sections 3, 4 and 5 above carry the
estimates themselves.** Filter `x$coefficients` on `model`, `drug`, `screen`,
`readout`, `ruler == "ox_gsva"`, `myc_estimator == "myc_ma"` to reproduce them.

**Releases:** GDSC **8.5 (fitted 27Oct23)**, DepMap **26Q1**, PRISM Repurposing
**24Q2**, GSVA **2.6.2**.

**THE DRY-RUN CAVEAT.** Every number above comes from a dry run with output
paths redirected to the scratchpad. **`results/bh3_oxphos_gdsc.rds` is not yet
on disk and `outputs/` was not written.** The caveat is discharged when the
author sources `scripts/E31_bh3_mimetic_oxphos.R` in Positron and the run
reproduces this note object by object, as `E16` through `E30` each did.
