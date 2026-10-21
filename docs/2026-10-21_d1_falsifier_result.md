---
date: 2026-10-21
status: >
  RESULT NOTE, written late. It records a falsifier verdict computed on
  2026-09-01 by `scripts/E05b_d1_falsifier.R` and saved in
  `results/d1_falsifier.rds`, which no note, handoff or commit message had
  recorded. Nothing was re-run to write it: every number below is read from
  that saved object, or from the file or commit named beside it.
posture: EXPLORATORY. The verdict is a NEGATIVE for D1's proliferation-adjusted claim.
relates-to:
  - docs/2026-09-01_phase1_celldeath_findings.md (defines D1; now carries a pointer banner)
  - docs/2026-09-01_phase2_estimator_findings.md (S2, S3 and S6 dissect the raw D1 contrast)
  - scripts/E05b_d1_falsifier.R
  - results/d1_falsifier.rds (gitignored; `$built` = 2026-09-01 13:14:55 BST)
  - outputs/figures/E05b_d1_contrast_null.{png,pdf} (gitignored)
---

# D1 falsifier (E05b): FALSIFIED after proliferation adjustment

## 1. What D1 asserted

D1 is defined in `docs/2026-09-01_phase1_celldeath_findings.md`, section "D1.
The pro-death / pro-survival skew - and it is an OXPHOS finding, not a MYC one"
(L65). The finding, at L78-79:

> **MYC-high and OXPHOS-high tumours carry a pro-survival-skewed apoptotic
> transcriptome.**

The claim the note carried forward, at L99-102:

> **The defensible claim is: tumours with a high OXPHOS programme carry a
> pro-survival-skewed apoptotic transcriptome, and this is not explained by
> proliferation, purity or immune infiltrate. The apparent MYC version of the
> same statement is largely proliferation.**

Its evidence was a contrast of GSVA scores, pro-death minus pro-survival
(L86-91). Against OXPHOS this was -0.413 (TCGA) and -0.551 (SCAN-B) raw, and
-0.320 and -0.372 proliferation-adjusted.

## 2. What E05b tested

The note's falsifier list (L133-142, below the pointer banner added at L126;
added in `ec533a2`, 2026-09-01 00:43:54 +0100) named the test as item 2:
"**expression has not been matched between the two strata**, and that is the
single cheapest remaining test" (L139-140).

`scripts/E05b_d1_falsifier.R` ran that test as a null of the contrast itself
(script section 5):

- each draw takes an expression-matched gene set for the pro-death stratum and
  a non-overlapping one for the pro-survival stratum, then differences their
  mean per-gene rho;
- this repeats `$settings$draws` = 2000 times, over `$settings$bins` = 20
  expression bins, with `$settings$seed` = 20260831;
- it is done for three targets (MYC, OXPHOS subunits on GSVA, OXPHOS subunits
  on mitoPPS), raw and proliferation-adjusted (`PROLIF_DISJOINT`), in both
  cohorts;
- section 6 repeats it with every MitoCarta gene deleted from the strata and
  the background.

The genes entering the test (`$contrast_null`) number 493 pro-death and 571
pro-survival in TCGA, and 492 and 567 in SCAN-B. With MitoCarta deleted
(`$contrast_nomito`) they are 455 and 530, and 454 and 527.

**The statistic is not D1's.** `$rules$statistic` reads: "mean per-gene rho,
not the GSVA score correlation. The two diverge and the gene-level one is the
conservative choice". So E05b tests a gene-level version of D1's contrast.
D1's GSVA numbers (-0.320 / -0.372) were never themselves nulled.

## 3. The rule for calling it falsified

The rule is the same for every verdict row (script L265-268 and L277): **a
cell survives only if the observed contrast is negative AND outside the 95%
null in BOTH cohorts; otherwise it is FALSIFIED.**

How "pre-declared" it is, stated exactly:

- **The test** was declared in the note, in `ec533a2`, before the script
  existed. **The threshold** was not in the note. It was first written in the
  script.
- The rule is present in the script's first commit, `dafb301` (2026-09-01
  10:29:37 +0100). The two later commits to the script, `74cd1fa` and `482ab54`
  (the latter 13:02:54 +0100), did not touch the verdict logic.
- The object's `$built` field reads 2026-09-01 13:14:55 BST, after the last
  script commit. So the rule was committed before the only recorded run.
- `dafb301`'s message says "Machinery verified against real TCGA data before
  commit". Git cannot show whether any contrast was seen before the rule was
  committed.
- `$rules$declared_in` reads "docs/2026-09-01_phase1_celldeath_findings.md D1".

## 4. The verdict, as read from `results/d1_falsifier.rds$verdict`

Criterion on every row: "observed contrast negative and outside the 95% null,
BOTH cohorts".

| falsifier (verbatim) | result |
|---|---|
| B  contrast vs a matched-pair null, OXPHOS GSVA, unadjusted | SURVIVES |
| B  contrast vs a matched-pair null, OXPHOS GSVA, proliferation-adjusted | **FALSIFIED** |
| B  contrast vs a matched-pair null, OXPHOS mitoPPS, proliferation-adjusted | **FALSIFIED** |
| C  same, with every MitoCarta gene deleted | SURVIVES |
| D  the MYC axis, proliferation-adjusted (expected to FAIL - D1 says so) | FALSIFIED |

**The verdict is FALSIFIED for the OXPHOS contrast on both GSVA and mitoPPS
after adjustment for proliferation.** D1's claim that the skew "is not
explained by proliferation" does not survive its own declared falsifier.

### 4.1 The numbers behind it

Below are the proliferation-adjusted rows from `$contrast_null` (all genes)
and `$contrast_nomito` (MitoCarta deleted), rounded to 4 dp from the object.
"null 95%" is `null_lo` to `null_hi`.

| target | genes | cohort | observed | null 95% | outside_95 |
|---|---|---|---|---|---|
| OXPHOS GSVA | all | TCGA | -0.0345 | -0.0356 to 0.0338 | FALSE |
| OXPHOS GSVA | all | SCAN-B | -0.0402 | -0.0309 to 0.0217 | TRUE |
| OXPHOS mitoPPS | all | TCGA | -0.0255 | -0.0321 to 0.0294 | FALSE |
| OXPHOS mitoPPS | all | SCAN-B | -0.0246 | -0.0260 to 0.0171 | FALSE |
| MYC (control) | all | TCGA | -0.0024 | -0.0238 to 0.0195 | FALSE |
| MYC (control) | all | SCAN-B | -0.0093 | -0.0230 to 0.0124 | FALSE |
| OXPHOS GSVA | MitoCarta deleted | TCGA | -0.0373 | -0.0339 to 0.0342 | TRUE |
| OXPHOS GSVA | MitoCarta deleted | SCAN-B | -0.0442 | -0.0305 to 0.0222 | TRUE |
| OXPHOS mitoPPS | MitoCarta deleted | TCGA | -0.0268 | -0.0313 to 0.0299 | FALSE |
| OXPHOS mitoPPS | MitoCarta deleted | SCAN-B | -0.0273 | -0.0249 to 0.0180 | TRUE |

The unadjusted rows are all outside their nulls (`outside_95` TRUE), in both
cohorts, for all three targets, with and without MitoCarta genes. On OXPHOS
GSVA: -0.0442 TCGA and -0.0572 SCAN-B, with nulls -0.0364 to 0.0300 and
-0.0305 to 0.0243.

### 4.2 Three things that travel with the verdict

1. **The GSVA failure is TCGA alone, and narrow.** The observed -0.0345 sits
   inside a null whose lower bound is -0.0356. SCAN-B is outside. The mitoPPS
   failure is in both cohorts.
2. **Row C survives.** That row is proliferation-adjusted OXPHOS GSVA with
   every MitoCarta gene deleted, and it is outside the null in both cohorts.
   The verdict rows do not test mitoPPS with MitoCarta deleted. In
   `$contrast_nomito` that cell is outside in SCAN-B and inside in TCGA.
3. **The raw skew survives.** D1's unadjusted statement is not what fails. What
   fails is the qualifier "not explained by proliferation", on the gene-level
   statistic and under the rule above.

The MYC control behaved as D1 predicted: it was FALSIFIED, with observed
contrasts inside the null in both cohorts.

## 5. Why it went unrecorded

The script's last message is "Record the verdict in the D1 note BEFORE moving
on." (`scripts/E05b_d1_falsifier.R:350`). The D1 note was edited twice after
`$built`:

- `ed37713`, 2026-09-01 14:14:12 +0100;
- `9f2d45a`, 2026-09-01 15:52:34 +0100, which added the REINTERPRETED banner.

Neither edit recorded the verdict. Before this note, no file in `docs/`
mentioned `E05b` or `d1_falsifier`. `docs/handoff-2026-10-21.md`, written the
same day, flagged the gap first.

## 6. What depended on D1

E08's S2, S3 and S6 (`docs/2026-09-01_phase2_estimator_findings.md`)
dissect the **raw** contrast. They do not depend on its proliferation-adjusted
form.

- S2's values are -0.0442 / -0.0572 (all genes) and -0.0472 / -0.0620
  (MitoCarta deleted), at L327-329. They equal the unadjusted `observed`
  values in `$contrast_null` and `$contrast_nomito`, which survive.
- D1 itself was restated the same day ("**D1's number stands; D1's
  interpretation does not.**", L495).
- Its carried-forward part is the 44-gene machinery split, S6 → E10, E11, E14,
  E15 and the synthesis. That split rests on its own per-gene evidence.
- No later note, CLAUDE.md or the synthesis cites D1's proliferation-adjusted
  claim.

Nothing downstream depended on D1 holding, so no other result is affected.
