---
date: 2026-09-28
status: >
  RESULT, RUN AND VERIFIED. VERDICT = MB2-ALIGNED, on the rule fixed at
  a71e19d before any number was computed. It survives the proliferation
  companion. The analysis is not merged.
posture: EXPLORATORY. Nothing pre-registered. The direction came from
         Menegollo's published description, not from our data; the verdict rule
         was committed ALONE at a71e19d BEFORE forkscale_models.rds was opened.
branch: mb2-forkscale. NOT merged, NOT deleted.
script: scripts/E32_mb2_forkscale_alignment.R, committed unrun at 7312ed4,
        fixed at 6e9df6f. SOURCED BY THE AUTHOR 2026-09-28 at 16:47, and
        results/mb2_forkscale_alignment.rds, the CSV and the figure are on
        disk. The numbers below were written from a dry run beforehand; every
        one of them is now checked against the saved object and THE DRY-RUN
        CAVEAT IS DISCHARGED - see section 10 for exactly how, because the
        check differs from E31's.
relates-to:
  - docs/2026-09-28_mb2_forkscale_declaration.md (the declaration - read first)
  - docs/2026-09-28_mb2_forkscale_data.md (the inputs, and what F3-pre held)
  - myc_human_validation/docs/2026-08-29_G3_result_forkscale_availability.md (F3-pre)
  - Menegollo et al., Fig. 7 and Supplementary Fig. S7A
---

# Which bicluster programme is this arm's state? MB2. Result.

**EXPLORATORY.** Nothing here is pre-registered. **N3**: two exposure axes and a
transcript score; nothing is "primed". **No outcome or survival variable enters
this analysis at any point.** Nothing was written to `myc_human_validation`.

**The results object is gitignored, so the numbers live in this note.**

---

## 0. The answer

| | |
|---|---|
| **VERDICT** | **MB2-ALIGNED**, on the rule fixed at `a71e19d` before any number existed |
| **the discriminating pair** | `rho(M_a, MB2 \| MB1)` = **+0.688 [+0.648, +0.723]**; `rho(M_a, MB1 \| MB2)` = **-0.195 [-0.247, -0.141]**. n = 1,037. The two intervals **do not overlap** |
| **paired difference** | **+0.884 [+0.802, +0.959]**, bootstrapped on the same resamples - a genuine paired contrast, not two independent intervals compared by eye |
| **the sign flip is the striking part** | Conditioning does not merely shrink MB1, it **reverses** it. At fixed MB2, MYC activity runs **negative** on MB1 forkscale - and on `M_b` it is **-0.586** |
| **does it survive proliferation?** | **Yes, and this was the main threat.** On the proliferation-stripped estimator it is **+0.653 [+0.610, +0.691]**, barely moved. With `PROLIF_DISJOINT` as a covariate it roughly halves to **+0.375 [+0.317, +0.430]** but stays firmly positive |
| **all three estimators agree** | `M_a` +0.688, `M_b` +0.686, `M_c` (GISTIC) +0.324 on `MB2 \| MB1`; all three negative or null on `MB1 \| MB2` |
| **what it licenses** | **Exactly one sentence** (section 8). **Nothing about outcome.** |

---

## 1. What was already known, and what this adds

**This must be read before the numbers, because it bounds the claim.**

`forkscale_models.rds` was unread from 2026-08-30 until this session, and the
declaration was committed alone at `a71e19d` before it was opened. **But the
object already contained the MB2 marginals**, computed in the validation repo:
MB2 vs `M_a` **+0.739**, MB1 vs `M_a` **+0.429**. This was recorded at `fb98f33`,
in the commit that opened it, precisely so it could not look like a discovery
made afterwards.

**So the marginal direction was computable before this analysis existed.** What
is new here:

1. **The partials** - the declaration's discriminating quantity, not in that
   object, and the thing that turns "MB2 is bigger" into "MB1 reverses sign".
2. **This repo's own scoring on both sides**, so the contrast is internally
   consistent rather than assembled from two repos.
3. **Bootstrap intervals**, and the paired difference.
4. **`M_c`**, which F3-pre never computed for any fork.
5. **The proliferation companion**, which is the reason the verdict is worth
   anything at all (section 5).

---

## 2. What was measured

| | |
|---|---|
| forkscale patients | **1,037** (per-bicluster files, the primary source) |
| joining to this repo's 1,095 TCGA patients | **1,037 - 100%** |
| **overlap with Block C** | **885**, above the declaration's expected 880 and well above the 800 stop |
| `M_c` (GISTIC) missing | 52 of 1,037 |
| ER status missing | 45 of 1,037 |
| `rho(MB1, MB2 forkscale)` | **+0.712** - the Supplementary S7A diagonal, measured |

**Forkscale was recomputed here** as `pc1 / index`, plain form, from the pinned
snapshot (upstream `8fdbb34`, blob SHAs re-verified this session). `index` was
asserted to be a clean permutation of `1:1037` for both biclusters - without
that, forkscale is not the quantity the upstream scripts computed.

**Cross-check against the assembled 849-row frame: max |difference| = 0** for
both MB1 and MB2 across all 849 shared rows. The definition is confirmed, not
assumed.

**The log form** is `Inf` at `index == 1` because `log(1) = 0`. Exactly **one
MB1 and one MB2 sample** were set to `NA` explicitly rather than being left to
survive `complete.cases()`, which they would have done.

**Why the partial and not the marginal.** The two forkscales correlate at
**+0.712**. "Different programmes" does not mean orthogonal axes, so marginals
cannot separate them - which is the whole reason the declaration named the
partial in advance.

---

## 3. THE PRIMARY - all three estimators

Partial Spearman; 10,000 bootstrap resamples with the **ranking re-done inside
every resample**, not the coefficient alone.

| estimator | quantity | n | rho | 95% CI | CI excludes 0 |
|---|---|---|---|---|---|
| **`M_a`** | **MB2 \| MB1** | 1037 | **+0.688** | **[+0.648, +0.723]** | **yes** |
| **`M_a`** | **MB1 \| MB2** | 1037 | **-0.195** | **[-0.247, -0.141]** | **yes** |
| `M_a` | MB2 marginal | 1037 | +0.748 | [+0.714, +0.778] | yes |
| `M_a` | MB1 marginal | 1037 | +0.442 | [+0.391, +0.490] | yes |
| `M_b` | MB2 \| MB1 | 1037 | **+0.686** | [+0.649, +0.719] | yes |
| `M_b` | MB1 \| MB2 | 1037 | **-0.586** | [-0.622, -0.546] | yes |
| `M_b` | MB2 marginal | 1037 | +0.444 | [+0.391, +0.493] | yes |
| `M_b` | MB1 marginal | 1037 | -0.052 | [-0.109, +0.005] | **no** |
| `M_c` | MB2 \| MB1 | 985 | **+0.324** | [+0.267, +0.378] | yes |
| `M_c` | MB1 \| MB2 | 985 | -0.056 | [-0.117, +0.005] | **no** |
| `M_c` | MB2 marginal | 985 | +0.397 | [+0.341, +0.447] | yes |
| `M_c` | MB1 marginal | 985 | +0.249 | [+0.189, +0.306] | yes |

**All three point the same way**, on a signature score, a regulon score and a
copy-number call - three different kinds of measurement of MYC.

**`M_c` is GISTIC amplification**, a discrete call with 52 NA here. Its weaker
magnitude is a property of a four-level call against a continuous axis, not a
weaker finding, and **F3-pre stored no `M_c` for any fork, so it has no
cross-check**.

### 3.1 The sign flip, which is the real result

The declaration only asked whether MB2 was *larger*. What the data give is
stronger: **conditioning on MB2 reverses MB1**. Marginally MB1 looks
MYC-associated (+0.442); at fixed MB2 it is **-0.195**. On `M_b` the marginal is
already null (-0.052) and the partial is **-0.586**.

**Read plainly: the apparent MYC association of MB1 is borrowed from its +0.712
correlation with MB2.** Once that is removed, MB1 carries no MYC association and
on two of three estimators carries a negative one. **That is Menegollo's own
description of MB1 as Myc-independent, recovered from our data rather than
assumed** - and it is the cleanest agreement between the two studies anywhere in
this arm.

### 3.2 Robustness

| | MB2 \| MB1 | MB1 \| MB2 |
|---|---|---|
| primary, n = 1,037 | +0.688 [+0.648, +0.723] | -0.195 [-0.247, -0.141] |
| **log form** of forkscale, n = 1,035 | +0.670 [+0.629, +0.706] | -0.187 [-0.237, -0.134] |
| **Block C only**, n = 885 | +0.680 [+0.635, +0.720] | -0.182 [-0.239, -0.123] |

Neither the paper's preferred log variant nor the Block C restriction moves it.

---

## 4. Cross-check against F3-pre - reported, NOT reconciled

Two repos' scoring, two different n. The declaration says report the difference
and do not silently replace one with the other.

| fork | axis | F3-pre (n=885) | here (n=1037) | delta |
|---|---|---|---|---|
| MB1 | `M_a` | +0.429 | +0.442 | +0.012 |
| MB2 | `M_a` | +0.739 | +0.748 | +0.009 |
| MB1 | `M_b` | -0.031 | -0.052 | -0.021 |
| MB2 | `M_b` | +0.470 | +0.444 | -0.026 |
| MB1 | OXPHOS (GSVA) | +0.529 | +0.536 | +0.007 |
| MB2 | OXPHOS (GSVA) | +0.480 | +0.481 | +0.001 |
| MB1 | OXPHOS (mitoPPS) | +0.418 | +0.426 | +0.007 |
| MB2 | OXPHOS (mitoPPS) | +0.414 | +0.421 | +0.008 |
| MB1 | `PROLIF_DISJOINT` | +0.503 | +0.512 | +0.009 |
| MB2 | `PROLIF_DISJOINT` | +0.874 | +0.874 | 0.000 |

**Maximum |delta| = 0.026, across 10 quantities and two independent scoring
pipelines.** The two repos agree. That is a reassurance about the scoring, not a
result.

---

## 5. THE PROLIFERATION COMPANION - post-hoc, and the reason the verdict holds

**Not in the declaration.** Chosen after reading `forkscale_models.rds` and
declared at `fb98f33` **before the script was written**, for a specific reason:
**`MB2_forkscale` correlates with `PROLIF_DISJOINT` at +0.874, harder than its
+0.739 with `M_a`.** CLAUDE.md trap 3 requires a proliferation-adjusted
estimate; here the generic rule had a concrete target. If the partial had died,
"MYC tracks MB2" would have been "proliferation tracks MB2" - which is roughly
Menegollo's description of **MB1**.

**It did not die.**

| quantity | n | rho | 95% CI |
|---|---|---|---|
| **`M_a__PROLIFSTRIP` ~ MB2 \| MB1** | 1037 | **+0.653** | **[+0.610, +0.691]** |
| `M_a__PROLIFSTRIP` ~ MB1 \| MB2 | 1037 | -0.164 | [-0.217, -0.109] |
| **`M_a` ~ MB2 \| MB1 + `PROLIF`** | 1037 | **+0.375** | **[+0.317, +0.430]** |
| `M_a` ~ MB1 \| MB2 + `PROLIF` | 1037 | -0.147 | [-0.199, -0.091] |
| `PROLIF` ~ MB2 marginal | 1037 | +0.874 | [+0.858, +0.887] |
| `PROLIF` ~ MB1 marginal | 1037 | +0.512 | [+0.462, +0.557] |

**Two operations, two answers, both surviving, and they are not the same
operation.** Stripping the proliferation genes out of the MYC signature leaves
the partial **essentially intact** (+0.688 to +0.653). Holding proliferation
fixed as a covariate **roughly halves it** (+0.375) - a large attenuation that is
reported as such - but the interval is nowhere near zero, and the MB1 mirror
stays negative under both.

**The honest reading.** MB2 forkscale really is close to a proliferation axis,
and a substantial part of the raw `M_a`-MB2 association is shared with
proliferation. **What survives is a MYC association at fixed proliferation that
MB1 does not have, in either direction.** The contrast between the two
biclusters - which is the declaration's question - is not a proliferation
artefact.

---

## 6. Declared UNINFORMATIVE in advance - OXPHOS, and it is NOT used

Both biclusters carry increased mitochondrial biogenesis, so this comparison
discriminates nothing. Declaration P2 says it will not be presented as support
**whatever it returns**. It is computed and reported here and **it is not used
in section 8**.

| ruler | MB2 \| MB1 | MB1 \| MB2 | MB2 marginal | MB1 marginal |
|---|---|---|---|---|
| `ox_gsva` | +0.167 [+0.108, +0.225] | **+0.315** [+0.259, +0.370] | +0.481 | +0.536 |
| `ox_ppd` | +0.186 [+0.125, +0.245] | +0.197 [+0.138, +0.257] | +0.421 | +0.426 |

It behaves exactly as the declaration predicted: **no MB2 preference at all**.
On `ox_ppd` the two partials are indistinguishable; on `ox_gsva` the lean is
towards **MB1**.

**One observation, and the declaration bars it from being support.** The MB2
preference is specific to MYC and is not shared by the OXPHOS axis, so it is not
a generic "everything loads on MB2" artefact. **That is recorded as an
observation and is not counted as evidence**, because P2 was fixed in advance
and is not now reinterpreted because the numbers came out conveniently.

---

## 7. P3 - CONTEXT, never evidence

MB2_UF is **ER-negative only** in Menegollo's description. No test statistic is
computed here and none should be.

**ER by STATE level** (STATE read read-only from the frozen repo; this repo does
not build it):

| STATE | n | % ER-negative |
|---|---|---|
| `MYC_low` | 476 | **8.0%** |
| `MYC_high_OXPHOS_low` | 172 | 41.3% |
| `MYC_high_OXPHOS_high_unbuffered` | 126 | 27.8% |
| `MYC_high_OXPHOS_high_buffered` | 168 | **44.0%** |

The declaration quotes 41.2% against 27.9% from the STATE freeze; this table is
restricted to the 1,037 forkscale patients with known ER, which is why level 4
reads 44.0% rather than 41.2%. **The direction is the same and neither number is
a test.**

**ER by forkscale quadrant** (medians within this cohort):

| quadrant | n | % ER-negative |
|---|---|---|
| **MB1 low, MB2 high** | 113 | **74.3%** |
| MB1 high, MB2 high | 378 | 25.9% |
| MB1 low, MB2 low | 389 | 10.3% |
| **MB1 high, MB2 low** | 112 | **6.2%** |

**The quadrant that is high on MB2 and low on MB1 is 74.3% ER-negative; the
mirror quadrant is 6.2%.** That is a twelve-fold difference and it is exactly
the direction Menegollo's description implies. **It is consistency and it is
recorded as consistency. It is not evidence, it is not a test, and no p-value
is attached to it.**

---

## 8. What this licenses, and what it does not

**LICENSED - one sentence, in section 4 or the Discussion:**

> The MYC-high, OXPHOS-high state described here corresponds to the
> MYC-associated bicluster (MB2) rather than to the proliferation-driven one
> (MB1) whose prognostic behaviour has previously been reported.

This cites the companion paper as a **positive rather than as a rival**, and it
explains the ER-negative skew of the buffered cell without special pleading.

**Three qualifications that belong with that sentence if it is written:**

1. **The marginal direction was already computable** from `forkscale_models.rds`
   on 2026-08-30. The partials, the intervals, the sign flip and the
   proliferation companion are new; the bare direction is not.
2. **MB2 forkscale is close to a proliferation axis** (`rho` = +0.874 with
   `PROLIF_DISJOINT`). The association survives both proliferation treatments,
   but covariate adjustment halves it, and a reader who checks will find that.
   Say it rather than be caught by it.
3. **The OXPHOS half of the comparison discriminates nothing** and was declared
   uninformative before it was computed. The sentence rests on the MYC contrast
   alone.

**NOT LICENSED.**

- **Anything about outcome or prognosis.** No survival variable entered this
  analysis. Menegollo report no survival analysis for MB2 anywhere; every
  survival panel in that paper is MB1. **A correlation between two exposure axes
  licenses no prognostic claim**, and the fact that MB2 has no published
  survival analysis is a gap, not an invitation.
- **Any claim that the two programmes are biologically distinct.** That is
  Menegollo's result, not ours. This is an alignment check on one axis.
- **Any causal reading.** Cross-sectional and descriptive, 1,037 TCGA samples.

---

## 9. What would change this answer

1. **A scoring difference.** Two repos already agree to |delta| <= 0.026, so
   this is close to closed.
2. **A proliferation treatment that kills it.** Two were tried and both
   survived, but they are not the only two. An estimator that removed
   proliferation more aggressively than `__PROLIFSTRIP` would be a real test.
3. **METABRIC.** The snapshot carries METABRIC sort data but it is incomplete
   here (see `data/menegollo_biclusters/README.md`, "The METABRIC gap"), and
   MB3's sign convention differs between cohorts. A METABRIC replication of the
   partial would be the natural next step and is **not** done here.
4. **MB3 as a control axis.** Not used here at all. If the MB2 preference also
   appeared against MB3, the contrast would be less specific than it looks.

---

## 10. Where the numbers live

| what | where |
|---|---|
| every correlation, with n, rho, CI and what it was adjusted for | `results/mb2_forkscale_alignment.rds`, `x$correlations` - **gitignored** |
| the same, as a table | `outputs/tables/E32_forkscale_correlations.csv` - **gitignored** |
| the verdict and its three inputs | `x$verdict`, `x$verdict_inputs`, `x$cis_overlap` |
| the recomputed forkscales | `x$forkscale` |
| the F3-pre cross-check | `x$f3pre_crosscheck` |
| P3 context tables | `x$er_by_state`, `x$er_by_quadrant` |
| the figure | `outputs/figures/E32_mb1_vs_mb2_forkscale.png` (gitignored) **and `docs/figures/2026-09-28_E32_mb1_vs_mb2_forkscale.png` (tracked)** |

**Because the object is gitignored, sections 3 to 7 carry the estimates
themselves.**

**THE DRY-RUN CAVEAT IS DISCHARGED - by a different check from E31's, stated
so it is not mistaken for the stronger one.** The author sourced
`scripts/E32_mb2_forkscale_alignment.R` in Positron on 2026-09-28 at 16:47.

**The object-by-object comparison E16 to E31 each received was not possible
here**: the dry-run object lived in a session scratchpad that no longer exists,
so there is no second object to compare against. Two checks were run instead:

1. **Every number this note quotes was checked against the saved object.** All
   **26** quoted correlations match on n, rho and both interval bounds to three
   decimals; so do the four OXPHOS marginals, the verdict, the non-overlap flag,
   the paired difference, `rho(MB1, MB2)`, both cross-check maxima, all six
   counts, and every row of both ER tables.
2. **The tracked figure came back byte-identical** to the one committed at
   `b633dd4`. It embeds every plotted point for 1,037 patients and prints the
   verdict and both primary intervals in its subtitle, so this is an
   independent check on the data the dry run plotted, not only on the numbers
   the note happens to quote.

The bootstrap is seeded (`PROJECT_SEED`), so identical intervals are the
expected outcome, not a coincidence. **What this cannot rule out** is a
difference in an element of the saved object that this note does not quote and
the figure does not draw. That is a smaller gap than it sounds - the note
quotes every declared and companion quantity - but it is a gap, and it is why
this paragraph does not say "object by object".

---

## 11. The branch

**`mb2-forkscale` is NOT merged and NOT deleted**, whatever the verdict. The
commit sequence is the audit trail and it is the point of doing it this way:

| commit | what |
|---|---|
| `a71e19d` | the declaration, **alone**, before `forkscale_models.rds` was opened |
| `fb98f33` | the inputs, and what opening that object revealed |
| `7312ed4` | the script, **before it had been run** |
| `6e9df6f` | two dry-run defects fixed |
| `b633dd4` | the result note and the figure |
| this one | the dry-run caveat discharged against the author's run |
