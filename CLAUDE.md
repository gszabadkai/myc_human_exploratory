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
found nothing supported. It is not reopened and it does not go in the paper.

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
number is `E45`, and `E31` is the only analysis still off `main`. No phase is
open.**

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
scripts/       numbered R pipeline, E00-E44. `E31` lives on an unmerged
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
