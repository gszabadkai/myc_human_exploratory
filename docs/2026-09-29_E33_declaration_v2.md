---
date: 2026-09-29
status: >
  DECLARATION - pre-data. Nothing below has been computed. Exposure-side only;
  no outcome variable is read, and none is joined.
  SUPERSEDES an earlier draft of this note dated 2026-09-29 which framed limb 1
  as a test of whether the configuration is MB2-specific. That framing was
  dropped: E11 had already shown the configuration's ordering survives
  proliferation adjustment, and MB2 forkscale runs +0.874 with proliferation, so
  the expected verdict was UNRESOLVED and both branches were awkward. What
  remains is descriptive, plus one limb that is genuinely informative.
repo: myc_human_exploratory. `myc_human_validation` and `myc_mouse` read-only.
relates-to:
  - docs/2026-09-28_mb2_forkscale_result.md (E32, MB2-ALIGNED)
  - E10, E11, E14 (the configuration, the conditioning asymmetry, specificity)
  - myc_human_validation/docs/2026-08-30_STATE_frozen_and_H4_buffer_declaration.md
  - Menegollo et al., Fig. 7, Supplementary Figs. S1F and S7A
decides:
  - What is descriptive here and what is claim-relevant, so the two are not
    conflated when the result is written up.
  - Direction and verdict for limb 3, the only limb that could sharpen a
    manuscript sentence.
next-action: commit alone, then data note, then script unrun, then result.
---

# E33 — characterising the MB1 / MB2 distinction, and where the death signal sits

## 0. Purpose, stated honestly

**Primarily descriptive.** Section 4 and the Discussion currently carry five
overlapping variables — the MYC x OXPHOS quadrant, `STATE`, forkscale, UF/LF
assignment, and ER/PAM50 — and the manuscript describes the MB1/MB2 distinction
approximately because the quantities needed to describe it precisely were never
computed. This note computes them.

**One limb is claim-relevant**, limb 3. The rest sharpens naming, not status.

**What this cannot do, declared up front.** No OXPHOS-directed or BH3-directed
intervention exists in any cohort available here. **Nothing in this analysis can
support a treatment-selection claim**, and the result note must say so in those
words.

## 1. LIMB 1 — descriptive. Two quantities E32 never reported.

E32 partialled MB1 on **MB2**, not on proliferation, and reported ER by quadrant
but not MYC. Both gaps are one line each.

- **MYC activity by forkscale quadrant.** `M_a`, `M_b`, `M_c` and
  `FELSHER__PROLIFSTRIP`, median and IQR, in each of the four quadrants formed
  by median splits of MB1 and MB2 forkscale. Cell sizes stated (E32 gave
  n = 112 and n = 113 for the two off-diagonal cells).
- **`rho(M_a, MB1_forkscale | PROLIF_DISJOINT)`**, beside E32's
  `rho(M_a, MB1 | MB2)` = −0.195 and the marginal +0.442.

**No verdict, no prediction.** These exist so the manuscript can say what MB1_UF
and MB2_UF are rather than approximate it.

**State in the note, because it limits what these numbers can mean:** MB2
forkscale runs +0.874 with `PROLIF_DISJOINT`, so "conditioned on MB2" and
"conditioned on proliferation" are close to indistinguishable in these data. A
difference between the two partials is reportable; an interpretation of which
one is doing the work is not.

**Expect a gradient, not two states.** MB1's marginal with MYC is +0.442 and
MB2's is +0.739. MB1_UF is **not** MYC-low in absolute terms, and the write-up
must not say it is.

## 2. LIMB 2 — descriptive. The variable map.

The cross-tabulation that makes five variables legible at once, with the
**MYC x OXPHOS quadrant** as the shared row variable throughout.

Build the quadrant from the frozen constructor's own MYC and OXPHOS calls:
split `STATE` level 1 by the OXPHOS call it already computes, merge levels 3
and 4. **Assert before proceeding** that Q3 equals level 2 exactly and Q4 equals
levels 3 + 4 exactly.

```
Q1 MYC-low  OXPHOS-low     Q2 MYC-low  OXPHOS-high
Q3 MYC-high OXPHOS-low     Q4 MYC-high OXPHOS-high
```

Cross-tabulate the quadrant against: `STATE`; MB1 UF/LF/none; MB2 UF/LF/none;
PAM50; ER. Marginal counts only, every cell size stated.

**Do not conflate forkscale median splits with UF/LF assignments.** The first is
a cut on a continuous score; the second is an MCbiclust sample assignment with
an unassigned `none` group. Both appear here.

**Sanity check, one line:** recompute the MYC call on `FELSHER__PROLIFSTRIP` and
cross-tabulate against the frozen `M_a` call. E32 gives grounds to expect close
agreement (+0.688 against +0.653 on the same partial). **This is a sanity check
and nothing is gated on it.** If the calls diverge substantially that is a
finding worth reporting, not a reason to stop.

## 3. LIMB 3 — the claim-relevant one. Trigger arm against guardian arm.

**Why this one is different.** The therapeutic stratification rests on a premise
that has never been tested directly: that the standing pro-apoptotic signal is
confined to tumours whose respiration is MYC-coupled. The two arms of the
configuration separate that premise from its alternative.

Same statistic as E32, on two gene groups reported **separately and always**:

```
rho( arm , MB2_forkscale | MB1_forkscale )   against
rho( arm , MB1_forkscale | MB2_forkscale )
```

| arm | genes |
|---|---|
| **trigger** | `BBC3`, `BID`, `BIK`, `BAD` |
| **guardian** | `BCL2L1`, `MCL1`, reported individually as well as together |

Partial Spearman, 10,000 bootstrap resamples with ranking redone inside each.
Marginals alongside. TCGA n = 1,037; overlap with the Block C set reported
exactly (expect ~885). The composite configuration score is reported third, for
continuity with E11, and is **not** the readout the verdict turns on.

**Verdict, fixed now:**

| | condition | what it would mean |
|---|---|---|
| **TRIGGER MB2-SPECIFIC** | trigger arm `MB2\|MB1` positive, CI excludes zero, exceeds `MB1\|MB2`; guardian arm does not show the same separation | the standing death signal is confined to MYC-coupled tumours. Sharpens the Discussion premise |
| **BOTH ARMS SEPARATE** | both show it | the whole configuration is MB2-associated; descriptive, does not isolate the trigger |
| **NEITHER SEPARATES** | neither | the configuration is not programme-specific. Reportable, and consistent with E11's finding that it is an OXPHOS correlate |

**All three outcomes are acceptable and all three get written.** NEITHER is not
a failure: E11 already suggests it, and knowing it would stop the manuscript
implying a programme-specificity it does not have.

## 4. Declared uninformative in advance

- **OXPHOS against either forkscale.** Both biclusters carry mitochondrial
  biogenesis, so it discriminates nothing. Computed, reported, never presented
  as support, as in E32.
- **Any prognostic reading.** No outcome variable enters. Note that MB1_UF and
  MB2_UF would be expected to have poor outcome for **different** reasons,
  proliferation and MYC-coupling respectively, so survival could not separate
  them and is not attempted here.

## 5. Traps

- `forkscale = pc1 / index`, **plain form**; the log form is `Inf` where
  `index == 1` and `Inf` survives `complete.cases()`.
- **MB3 is not used.** If it ever enters, METABRIC negates `MB3.pc1`, TCGA does
  not.
- Use the **per-bicluster** files, not `TCGA.all.biclusters.RNAseq.Rdata`, which
  covers only 849 samples after an upstream join.
- `STATE` is not built in this repo; read the frozen constructor read-only.
- Reuse E32's forkscale construction verbatim so limbs 1 and 3 sit on the same
  quantities. Report any difference.

## 6. What a result licenses

**Limbs 1 and 2 license naming, not claims.** They let the manuscript describe
MB1_UF, MB2_UF and the quadrant precisely, and they supply one Extended Data
table.

**Limb 3, if TRIGGER MB2-SPECIFIC**, licenses one clause in the Discussion: that
the standing pro-apoptotic signal is a property of tumours whose respiration is
MYC-coupled. It does **not** license any statement about response to an
OXPHOS-directed or BH3-directed agent, because no such data exists in any cohort
here.

## 7. What this cannot do

Cross-sectional, descriptive, one cohort, bulk transcriptomes, and resting on a
bicluster assignment produced by another study. It cannot establish that the
barrier operated, only whether the configuration is distributed as the barrier
would predict.
