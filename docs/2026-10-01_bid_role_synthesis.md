---
date: 2026-10-01
status: >
  SYNTHESIS. No new computation: every number is read from a committed result,
  a saved object on disk, or a read-only sweep of the two sibling repos. Written
  to decide what role BID can be given, after E34 made it conspicuous.
posture: >
  The human exploratory results are EXPLORATORY and hypothesis-generating. The
  validation repo's results are PRE-REGISTERED where its dated notes say so, and
  are flagged as such below. The mouse results are an INTERVENTION arm plus an
  n = 24 correlational arm, and the two are not interchangeable.
question: Does BID follow MYC rather than OXPHOS, and can that explain why
          MYC + OXPHOS kills in the mouse while raising OXPHOS alone does not?
answer: >
  Half of it holds. BID does follow MYC, and it is the best-evidenced MYC
  main-effect claim in this repo. It cannot carry the killing asymmetry,
  because that is an interaction and every interaction test on BID, in all
  three repos, is null or collapses. All three repos independently designated
  BID a negative control; only the human exploratory data refuse to let it be
  one.
relates-to:
  - docs/2026-10-01_e34_result.md, docs/2026-09-30_e33_result.md
  - docs/2026-09-03_human_arm_for_mouse_reconciliation.md (section 6; M3; N1; traps 9, 12)
  - myc_mouse docs/2026-09-02_myc_oxphos_priming_gate_model.md (read-only)
  - myc_human_validation docs/2026-08-29_block_c_result_H1_not_supported.md (read-only)
---

# What is the role of BID?

**Nothing here is a new analysis.** Every number is read from a committed
result, a saved object, or a read-only sweep. Nothing was written to
`myc_human_validation` or `myc_mouse`.

---

## 0. The answer in six lines

| | |
|---|---|
| **Does BID follow MYC?** | **Yes, and it is the single best-evidenced MYC main-effect claim in this repo.** Rank 1 of 44 genes vs MYC in TCGA, positive with every interval excluding zero in all 20 `E22` cells, and uniquely untouched by proliferation adjustment |
| **Instead of OXPHOS?** | **No.** BID follows both. In TCGA its OXPHOS contrasts are *larger* than its MYC ones; in SCAN-B its OXPHOS association collapses to +0.091. The MYC association replicates; the *dominance over OXPHOS* does not |
| **Is it the only one?** | **No.** `PMAIP1` also keeps a MYC association after conditioning on OXPHOS — in the opposite direction. The repo's standing phrasing is "**two genes, not a programme**" |
| **Can it explain the killing asymmetry?** | **No, not on present evidence.** That asymmetry is an **interaction**, and every interaction test on BID is null or collapses: mouse `bMX` p = 0.36, validation interaction 55% artefact, human `BID/BCL2` killed by its own spline falsifier |
| **The awkward fact** | **All three repos designated BID a negative control** — the mouse gate model by name, the validation study by pre-registration. In the mouse, a strong BID is a *falsification* condition |
| **What it can be** | A **MYC-responsive BH3-only transcript**, reported as one gene. Not a mediator, not a mechanism, not a programme |

---

## 1. Your statement, verified and refined

> *"BID is the only one which fits the hypothesis that a trigger follows MYC
> instead of OXPHOS."*

**Supported.** Among the four trigger-arm genes (`BBC3`, `BID`, `BIK`, `BAD`),
BID is the only one whose MYC association is positive with the interval
excluding zero in **all four** cohort x instrument cells of `E34`, on **both**
MYC contrasts — at OXPHOS-high (`Q4 - Q2`: +0.188, +0.343, +0.200, +0.253) and
at OXPHOS-low (`Q3 - Q1`: +0.313, +0.302, +0.172, +0.169). **BID rises with MYC
regardless of OXPHOS status**, which is the part your hypothesis needs.

`BAD` is its exact complement and sharpens the contrast: MYC-null in all four
cells, with the largest OXPHOS contrasts in the panel (+0.56 to +0.90).

**Refinement 1 — "instead of" is wrong.** BID follows *both* axes, and which
dominates is cohort-dependent. TCGA OXPHOS contrasts +0.35 to +0.60 exceed its
MYC contrasts; SCAN-B OXPHOS contrasts collapse to +0.001 to +0.137 and MYC
dominates. Mean |MYC| / |OXPHOS| per cell: 0.47, 0.88 (TCGA), 1.51, 4.90
(SCAN-B).

**Refinement 2 — "the only one" is wrong.** `PMAIP1` tracks MYC at least as
consistently (`Q4 - Q2` = -0.199, -0.330, -0.210, -0.265, every interval
excluding zero) and is OXPHOS-null. It simply runs **downward** and is not in
the trigger arm. `BCL2L2` behaves similarly in SCAN-B.

**The defensible form:** *among the four trigger-arm genes, BID is the only one
carrying a replicating MYC association present at both OXPHOS levels; `BAD` is
the purest OXPHOS gene; and `PMAIP1` is a second MYC-tracking BH3-only gene
running the opposite way.*

---

## 2. The case FOR BID, from this repo

This is stronger than I expected, and it is worth stating at full strength
before the counter-evidence.

### 2.1 It ranks first wherever MYC is the axis

| table | BID's rank |
|---|---|
| 15-gene BCL2 family vs MYC, pooled (`E08` D3) | **1 of 15** (+0.425 TCGA / +0.346 SCAN-B) |
| 44-gene machinery vs MYC, adjusted, TCGA (`E15`) | **1 of 44** (+0.357) |
| 44-gene machinery vs MYC, adjusted, SCAN-B (`E15`) | **4 of 44** (+0.243) |
| 12 priming genes vs MYC, adjusted (`E10` R5a) | **1 of 12** (+0.385 mean) |

### 2.2 The direct MYC effect at fixed OXPHOS is positive in all 20 cells

`E22` fitted `gene ~ MYC + OXPHOS` across 2 cohorts x 5 estimators x 2 rulers.
**BID's `beta1` is positive with every interval excluding zero in all 20**
(verified against `outputs/tables/E22_fit_a.csv`):

| cohort | dose (log2MYC) | signature | regulon | sig_clean | sig_entangled |
|---|---|---|---|---|---|
| TCGA, `ox_rel` | +0.297 | +0.335 | **+0.474** | +0.306 | +0.289 |
| TCGA, `ox_lvl` | +0.263 | +0.312 | **+0.481** | +0.280 | +0.277 |
| SCAN-B, `ox_rel` | +0.099 | +0.295 | **+0.522** | +0.227 | +0.232 |
| SCAN-B, `ox_lvl` | +0.088 | +0.312 | **+0.532** | +0.239 | +0.261 |

In `E22`'s specificity table BID is the extreme: **0 of 6 cells negative**, the
most positive range of the twelve (+0.093 to +0.530), against `BCL2L1`'s 6 of 6
negative.

### 2.3 It passes the test that killed the pre-specified genes

`E22`/trap 12 established that the sign of a MYC main effect on a BCL2-family
transcript is usually a function of the estimator's proliferation entanglement —
`BAD` and `BBC3` "fail it worst". **BID does not track entanglement on either
ruler** (spread 0.099 and 0.102, among the smallest of the twelve) and keeps its
sign across the whole 1.5%-23.5% range, +0.183 to +0.477.

### 2.4 It is the only gene proliferation adjustment leaves alone

From the `E34` object, mean `Q4 - Q2` over the four cells, m1 -> m3:

| gene | m1 | m3 | shift |
|---|---|---|---|
| **BID** | **+0.246** | **+0.257** | **+0.011** |
| BBC3 | +0.130 | +0.280 | +0.149 |
| BAD | -0.055 | +0.133 | +0.187 |
| BCL2L1 | -0.041 | +0.051 | +0.092 |
| PMAIP1 | -0.251 | -0.418 | -0.167 |

**Nothing else in the panel is stable.** `BBC3`'s and `BAD`'s apparent MYC
associations are largely *manufactured* by the adjustment — `BBC3` goes from
non-significant to significant in TCGA, `BAD` flips sign. This matters
especially because BID is a CollecTRI target of **E2F1**, so proliferation was
the obvious alternative explanation. It does not survive contact with the data.

### 2.5 It is the repo's one sanctioned MYC main-effect claim

The reconciliation document, section 6 — "**the one MYC main-effect claim in the
document**": *"`BID` (+0.12 to +0.28) and `PMAIP1`/NOXA (-0.05 to -0.21) keep a
MYC association after conditioning on OXPHOS, in both cohorts under all three
estimators."* Positive in all six cells, and neither gene is a member of any of
the three estimators, so it is not self-overlap.

---

## 3. The case AGAINST, from the same repo

Every item here is a committed, dated finding.

1. **The headline MYC number is pooling-inflated, and this is a formal dated
   qualification** (`docs/2026-09-01_phase1_celldeath_findings.md:193`,
   verified): *"QUALIFIED 2026-09-01 by `E08`. The MYC column above is inflated
   by pooling."* `BID` runs **0.42 pooled against 0.20-0.26 in the luminal
   strata**. The accompanying instruction is explicit: *"D3's BCL2-family claim
   should be made **against OXPHOS, not against MYC**."* Much of the headline
   MYC association is a between-subtype artefact.
2. **BID inverts in cell lines.** `E18`/CCLE: **0 of 4 rulers agree in sign with
   the tumours**, in breast (mean -0.225) and pan-cancer (-0.059). The B4
   addendum singles it out — `BCL2` and `BID` *"invert cleanly"*.
3. **Its OXPHOS association collapses between cohorts**: +0.329 TCGA to **+0.091**
   SCAN-B (adjusted). This is the named cause of the `BID/MCL1` ratio replication
   failure.
4. **It fails in SCAN-B Basal**: CIs include zero on 2 of 4 rulers. BID is not
   among the six genes clearing zero in all 8 Basal cells.
5. **It is purity-sensitive**: 2nd-largest shift of the twelve under purity +
   leukocyte adjustment, on two rulers. Not checkable in SCAN-B (trap 2).
6. **Its one interaction result was killed by its own falsifier.** `BID/BCL2`
   was one of five MYC x OXPHOS interactions replicating with the same sign in
   both cohorts — then F1 (splines) collapsed it from **0.081 to -0.003**, F2
   reversed cross-cohort replication to r = -0.33, F3 was non-monotone. Verdict:
   *"the apparent interaction is curvature in the main effects"* (N1).
7. **The two MYC-held genes move oppositely.** Proposition **M3**: *"MYC does not
   participate in the configuration... the two genes MYC holds on to move in
   opposite directions."* BID up, PMAIP1 down. **"It is two genes."**
8. **MYC *activity*, not MYC *dose*.** In SCAN-B the `dose` (log2MYC) estimator
   gives BID only **+0.099**, against +0.295 for the signature and +0.522 for
   the regulon. Trap 4 in its sharpest form: whatever BID responds to, it is not
   MYC transcript level.

---

## 4. The frozen validation study — BID was a PRE-REGISTERED NEGATIVE CONTROL

This is the one repo whose results carry pre-registered weight, and it reached
for BID deliberately.

- **Definition, pre-registered** in the plan dated 2026-08-27, before any model
  existed: *"Endpoint negatives: BID/BCL2L1, BAX/BCL2L1, BCL2L11/BCL2L1,
  BAK1/BCL2L1. The mouse says only PUMA/BCL-XL reverses. Human should show the
  same."* It survived the one revision that touched the panel: *"`BID`, `BAX`
  and `BAK1` remain valid negatives."*
- **Result**: `MYC:OXPHOS` on `BID/BCL2L1` = **+0.062, p = 0.012 / 0.015**.
- **But 55% of that is borrowed from the shared falling `BCL2L1` denominator.**
  The `log2(BID)` limb alone is **+0.0279**; the denominator contributes 0.0342.
  The decomposition is *"an algebraic identity... verified to machine precision,
  not an empirical finding"*.
- **It does not clear multiplicity**: Bonferroni over 114 coefficients is
  4.4e-4.
- **It was never carried into the pre-registered replication** — that named only
  `BCL2L11`, `BCL2L1` and `BBC3`.
- BID contributed to the failure of **decision clause 4** (*"endpoint negatives
  null"*) by moving when it was supposed to be silent.

**This is not evidence against the E34 finding.** The validation tested an
**interaction**; `E34` tests a **main effect**. Trap 1: they are different
questions. But it does mean BID's only pre-registered appearance is as a control
that misbehaved for an arithmetic reason.

---

## 5. The mouse — BID is the gene that must stay quiet

### 5.1 It is the designated negative control, by name

`docs/2026-09-02_myc_oxphos_priming_gate_model.md` §4.3, verified verbatim:
*"Negative-control endpoints unchanged, with the mouse's expected values now
attached: `BCL2L11:BCL2L1` (null), `PMAIP1:BCL2L1` (negative), `BAX:BCL2L1`
(weak, rides the denominator), **`BID:BCL2L1` (weak)**."*

And §5, Falsification: *"The model is dead if... `bMX` is as large for endpoints
the model does not name (`BAX`, **`BID`**, `NOXA`) as for the co-primaries, with
`BCL2L1` alone flat — the gate would then be a general apoptotic shift."*

**The role you are proposing for BID is, in the mouse's own framework, the
condition under which its model fails.** That is not fatal to the idea, but
promoting BID is an argument *against* the gate model, not a refinement of it.

### 5.2 And it behaved exactly as designated

| | `bMX` | p | pct expression-matched |
|---|---|---|---|
| `Mcl1` | +0.806 | 0.018 | 98.9 |
| `Bbc3` | +0.746 | 0.059 | 94.8 |
| **`Bid`** | **+0.342** | **0.362** | **76.6** |
| `Pmaip1` | -0.430 | 0.033 | 7.7 |
| `Bcl2l1` | -0.670 | 0.0026 | 2.1 |

Rank 5 of 10, upper-middle, nowhere near a tail. `Bid:Bcl2l1` = +0.588
(p 0.076), and **within-timepoint — the regime comparable to a human
cross-section — it collapses to +0.148, p 0.45**.

Elsewhere: gene-level MYC effect **+0.609 at 6W (padj 0.132, n.s.)**, +0.148 at
12W; genotype x time interaction -0.461 (p 0.23), in a panel the source doc
calls *non-corroborating*. **Low expression is flagged**: baseMean 148, with an
explicit note that a real half-log2 effect could be missed.

### 5.3 Under raised OXPHOS, mouse Bid goes DOWN

The construct-free orthotopic arm (PGC1a vs EV): **`Bid` -0.182 [-0.234,
-0.018]** — the opposite sign to the human. Both are post-hoc; the source doc
states *"none of which was pre-specified and all of which is reported without
per-gene FDR."*

### 5.4 The killing asymmetry already has an answer, and it is not BID

**Ndi1 and mtLbNOX do not kill MYAZ at baseline but do kill MYAZ + Tre-MYC.**
That is the MYC x OXPHOS interaction demonstrated **by intervention**, in the
form a human cross-section cannot show. The manipulation panel excludes **mass**
(Ndi1 and mtLbNOX expand nothing and kill), **ATP** (LbNOX makes none and
kills), and **membrane potential** (Ndi1 and mtLbNOX lower it, PGC1a raises it,
all three kill). The lethal variable is named as the **matrix NADH/NAD+ ratio**;
cytLbNOX killing "much less" also excludes the cytosolic relay.

**Bid appears nowhere in that chain.** Competing with an intervention result
using correlational transcript data is a weak position.

### 5.5 Two further constraints from the mouse

- The standing correction is *"respiration does nothing without MYC"*, **not**
  "respiration is protective" — the WT slope is -0.012 (p 0.95) once composition
  is adjusted. The guardian balance `Mcl1:Bcl2l1` *is* a genuine crossover; the
  trigger is not.
- `docs/2026-09-09_reprioritisation_narrative_v3.md` **forbids** setting a mouse
  `bMX` beside a human partial Spearman — *"different estimands, orthogonal by
  construction"*. That is precisely the juxtaposition this hypothesis wants.

---

## 6. The central tension

**Three independent efforts assigned BID the role of control, and in the human
exploratory data it refuses to play it.**

- The mouse gate model: negative control, and a strong BID is a falsification
  condition.
- The validation study: pre-registered endpoint negative.
- This repo: never the subject of its own analysis — only ever one row in a
  panel of 12, 15, 44 or 6.

And yet BID is rank 1 of 44 against MYC, positive in all 20 `E22` cells, and the
only gene in the panel that proliferation adjustment does not move.

**That is interesting, and it is also exactly the shape of a false positive**
reached by looking at a gene nobody pre-specified, in a repo whose own rule is
that no single cell of a 44 x 2 grid is a finding.

---

## 7. Decomposing the hypothesis

Your hypothesis has two parts, and they have very different standing.

**(a) BID follows MYC where the other trigger genes follow OXPHOS.**
**SUPPORTED**, with the two refinements in section 1. This is the repo's best
MYC main-effect claim and it survives entanglement, proliferation adjustment and
cross-cohort replication.

**(b) That explains why MYC + OXPHOS kills while OXPHOS alone does not.**
**NOT SUPPORTED.** The asymmetry in (b) is an **interaction**, and every
interaction estimand touching BID is null or collapses:

| test | result |
|---|---|
| mouse `bMX` on `Bid` | +0.342, **p 0.362**; within-timepoint p 0.78 |
| mouse `Bid:Bcl2l1` | +0.588, p 0.076; within-timepoint **+0.148, p 0.45** |
| validation `MYC:OXPHOS` on `BID/BCL2L1` (pre-registered) | +0.062, **55% denominator artefact**, limb +0.0279 |
| human `BID/BCL2` interaction | **collapsed by its own spline falsifier**, 0.081 -> -0.003 |

This is the same wall the whole project has hit: the MYC x OXPHOS interaction
estimand has failed in Block B, Block C, Block G, H4, N1 and N2.

### 7.1 The reading that would rescue it — and why it cannot be tested here

Your hypothesis does not actually *require* an interaction on BID. A
**conjunctive** model — death needs BID (supplied by MYC) **and** the
OXPHOS-driven remainder of the configuration — predicts main effects on
*different genes* and **no interaction on any single gene**. That is consistent
with everything above, including the mouse's own phenotypic observation that
*"death needs two inputs, proliferation needs roughly one"*.

**But it is untestable from transcripts**, for a reason this repo has already
committed to (N3): **BID must be cleaved to act**. Transcript abundance cannot
report tBID. The conjunction, if it exists, happens at a level no measurement in
any of the three repos reaches. The repo retired localisation language on
2026-09-02 for exactly this reason — MitoCarta membership for BID marks
*regulon membership, not where the protein sits*.

---

## 8. The mechanistic puzzle, which is genuinely open

**BID is not an annotated MYC target.** Verified in this repo's own snapshot:
CollecTRI carries **891 MYC->target edges and none to BID**. BID's annotated
regulators are TP53, HIF1A, ZBTB16, FOXO1, E2F1, TP73. BID is absent from the
Felsher signature — which *is* `M_a`, the estimator the association is measured
against. The twelve genes that *are* CollecTRI MYC targets are `BCL2L11`,
`PMAIP1`, `BBC3`, `BCL2`, `BCL2L1` — i.e. the ones with weaker or opposite MYC
associations.

Three facts that deepen rather than resolve this:

1. **The regulon estimator gives BID its largest effect** (+0.474 to +0.532),
   while BID is not a member of `M_b__MITOSTRIP`. So BID correlates best with a
   MYC-regulon-derived activity score without being in the regulon.
2. **The mouse chromatin data point the other way.** ChEA/Gray:
   **MYC -> Bid at FDR 1.8e-32** in luminal epithelium, and MYC reaches `Bbc3`
   *only* in the terminal end bud. On that reading BID is the default BH3-only
   gene MYC can touch in most chromatin contexts. The mouse repo flags this in
   capitals as *"a context statement, not a TF statement"*.
3. **BID is in `BILD_MYC_ONCOGENIC_SIGNATURE`** — one of the sixteen MYC sets,
   and a self-overlap hazard for that set alone. It is not one of the three
   claim estimators.

**An association this robust to a gene with no annotated MYC edge needs an
explanation.** It is currently unexplained, and that is a reason for caution,
not for excitement.

---

## 9. What BID would need to earn a role

In rough order of cost:

1. **A third cohort.** The repo already carries a declared, never-run falsifier
   aimed exactly here (reconciliation §10 item 4): *"The `BID`/`BBC3` numerator
   test named in 3.3b is the cheapest real falsifier left in the priming
   section, and it needs a third cohort."* A cohort behaving like TCGA should
   show `BID/MCL1` gaining; one behaving like SCAN-B should not.
2. **Resolve the pooling qualification.** The MYC claim must be shown to survive
   *within* subtype at full strength, not at the 0.20-0.26 the luminal strata
   give. `E34`'s `m2` already does part of this and should be read for BID
   specifically.
3. **Explain the missing MYC edge** (section 8) — or accept that the
   association is indirect and say so.
4. **Any protein-level readout.** tBID, or cleaved caspase-8. Without it, N3
   caps every sentence at "transcript".
5. **In the mouse, a direct test.** Bid was never a hypothesis gene there and at
   baseMean 148 with n = 24 the arm was underpowered to be one. "Never tested"
   is not "tested and refuted" — but it is also not support.

---

## 10. What must not be written

- **Never** that BID is "primed", or that any tumour is (N3).
- **Never** that BID mediates the MYC x OXPHOS interaction, or explains the
  mouse killing asymmetry. No interaction evidence survives.
- **Never** a mouse `bMX` beside a human partial Spearman — forbidden by the
  mouse repo, different estimands.
- **Never** "the pro-apoptotic genes rise with OXPHOS" — they do not move as a
  block, and BID and PMAIP1 move oppositely with MYC.
- **Never** the pooled `+0.425` MYC figure without its dated qualification.
- **Never** BID as a programme or a mechanism. The sanctioned phrasing is **two
  genes**.
- Nothing here licenses any treatment-selection claim.

---

## 11. Sources and epistemic status

| source | status |
|---|---|
| `E34` contrasts, m1/m3 shifts | EXPLORATORY; run and verified, 12/12 digests |
| `E22` `beta1` table | EXPLORATORY; verified against `outputs/tables/E22_fit_a.csv` |
| `E08` pooling qualification | EXPLORATORY; a **dated qualification**, verified in source |
| `E10`, `E15`, `E16`, `E18`, `E19`, `E33` | EXPLORATORY; read from committed notes and tables |
| reconciliation §6 MYC claim | EXPLORATORY; the document's **one** sanctioned MYC main-effect claim |
| CollecTRI MYC->BID absence | Verified by direct query of this repo's snapshot |
| validation `BID/BCL2L1` | **PRE-REGISTERED** as a negative control; result is 55% artefact |
| mouse gate model §4.3, §5, §3.3 | Verified verbatim read-only; BID is a **declared negative control** |
| mouse intervention panel | INTERVENTION — the strongest evidence in the project, and BID is absent from it |
| mouse orthotopic `Bid` -0.182 | **Explicitly post-hoc**, no per-gene FDR |
| mouse scripts 50/51 conclusions | **VOIDED** by the orthotopic identity correction — not used here |
| ChEA/Gray MYC->Bid | A **chromatin-context statement**, flagged as not a TF statement |

**Bottom line for the decision.** BID is worth one carefully qualified sentence
about a MYC-responsive BH3-only transcript, reported as one gene beside
`PMAIP1` running the other way. It is not worth a mechanism, and it cannot be
offered as the reason the mouse combination kills — that explanation is already
held, by intervention, by the matrix NADH/NAD+ ratio, and BID is not in it.
