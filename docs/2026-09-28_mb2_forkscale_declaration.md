---
date: 2026-09-28
status: >
  DECLARATION - pre-data. Nothing below has been computed. Written and committed
  BEFORE `myc_human_validation/results/forkscale_models.rds` was opened, and
  before any MB2 quantity existed.
repo: myc_human_exploratory (nothing pre-registered here; the frozen validation
  repo is read-only for this work and nothing is written to it)
relates-to:
  - myc_human_validation/docs/2026-08-29_G3_result_forkscale_availability.md
    (F3-pre, section 4b - the declared diagnostic this extends)
  - myc_human_validation/results/forkscale_models.rds (computed 2026-08-30,
    UNREAD at the time of writing)
  - myc_human_validation/data/menegollo_biclusters/ (the pinned snapshot)
  - Menegollo et al., Fig. 7 and Supplementary Fig. S7A
  - section 4 draft, limb (i) and the Menegollo clause
decides:
  - The direction, the statistic and the verdict rule for the MB2 comparison.
  - Which comparison is declared UNINFORMATIVE in advance, so that it cannot be
    presented as support afterwards.
next-action: commit this note alone, then open forkscale_models.rds, then run.
---

# Which bicluster programme does this arm's state correspond to?

## 0. What this is, and what it is not

**It is an addendum to an exposure-side diagnostic.** F3-pre was declared on
2026-08-29 and computed on 2026-08-30. No outcome variable appears here, none is
joined, and nothing is written to `myc_human_validation`.

**It is not pre-registered.** F3-pre is declared; this extension to MB2 is not in
the original plan. It lives in the exploratory repo for that reason.

**It is written before the answer exists.** `forkscale_models.rds` has not been
opened, so the MB1 half of the question is as unknown to us as the MB2 half. The
ordering is the point and it is auditable from the commit sequence.

## 1. The question

Menegollo describe two biclusters that both carry increased mitochondrial
biogenesis but arise from different programmes:

| | cell state | biogenesis driven by | ER |
|---|---|---|---|
| **MB1_UF** | mature luminal (mL) | proliferation signalling, **without** Myc/miRNA signals | almost exclusively ER+ |
| **MB2_UF** | luminal progenitor (LP) | the same biogenesis, **with** Myc/miRNA signals | **ER-negative only** |

Every survival panel in that paper, Fig. 7E to 7L and Supplementary Fig. S7C to
S7E, is **MB1**. MB2 has no survival analysis anywhere in it.

This arm's state is MYC-high and OXPHOS-high. The question is whether it
corresponds to the MYC-associated programme (MB2) or to the proliferation-driven
one (MB1) whose prognostic behaviour has already been reported.

## 2. What F3-pre already covers, so this does not duplicate it

F3-pre section 4b declares, for **MB1 only**: correlation of `MB1.forkscale`
with OXPHOS on both instruments, with M-a, M-b and M-c, with PRIME and with
STATE, plus whether the Block C interaction estimate survives adding forkscale
as a covariate. Verdict thresholds fixed in advance at REDUNDANT for
abs(rho) >= 0.70 on both instruments and INDEPENDENT for <= 0.30 on both.

**New here is the MB2 side and the contrast between the two**, nothing else.

## 3. The direction is fixed by the published description, not by us

Menegollo state that MB1_UF's biogenesis runs **without** Myc signals and
MB2_UF's runs **with** them. The prediction therefore follows from their text
and could not be tuned to our data:

- **P1, the test.** MYC activity correlates with MB2 forkscale, and not or much
  less with MB1 forkscale.
- **P2, declared uninformative.** The OXPHOS axis correlates with both. Both
  biclusters carry mitochondrial biogenesis, so this comparison discriminates
  nothing. **It will not be presented as support whatever it returns.**
- **P3, context and not a test.** STATE level 4 is 41.2% ER-negative against
  27.9% in level 3 (STATE freeze section 5), which is consistent with MB2_UF
  being ER-negative only. Reported as consistency, never as evidence.

## 4. Statistic

**Partial Spearman.** Supplementary Fig. S7A plots the two forkscales against
each other on a clear positive diagonal, so different programmes does not mean
orthogonal axes and marginal correlations will not separate them. The
discriminating quantities are:

```
rho( M-a , MB2_forkscale | MB1_forkscale )
rho( M-a , MB1_forkscale | MB2_forkscale )
```

Marginals reported alongside. Bootstrap confidence intervals, percentile,
10,000 resamples. Spearman because forkscale is severely skewed by construction.

**Both MB1 and MB2 are computed fresh in this repo on this repo's own scores**,
so the two sides of the contrast are internally consistent. F3-pre's stored MB1
values are then a **cross-check, not an input**: if the MB1 marginal here differs
materially from the one in `forkscale_models.rds`, that is a scoring difference
between repos and is reported rather than reconciled silently.

**MYC estimator.** M-a primary, as everywhere in this arm. M-b and M-c reported,
since F3-pre declared all three on the MB1 side.

## 5. Verdict, fixed now

| Verdict | Condition |
|---|---|
| **MB2-ALIGNED** | `rho(M-a, MB2 \| MB1)` positive with CI excluding zero, **and** larger than `rho(M-a, MB1 \| MB2)` |
| **MB1-ALIGNED** | the reverse |
| **UNRESOLVED** | anything else, including both positive with overlapping intervals |

UNRESOLVED is an acceptable outcome. The manuscript sentence it would license
does not exist, and limb (i)'s Menegollo clause is then dropped rather than
softened.

## 6. Traps, carried from the G3 note

- `forkscale.log` is **`Inf`** where `index == 1`. The plain form
  `pc1 / index` carries the verdict, as F3-pre declared; the log form is reported
  alongside only.
- **MB3 is not used here.** If it ever enters, METABRIC negates `MB3.pc1` before
  forming forkscale and TCGA does not, and that must be handled explicitly.
- The pre-assembled `TCGA.all.biclusters.RNAseq.Rdata` covers only 849 samples
  because of an upstream join. **The per-bicluster files are the primary
  source**; the assembled frame is a cross-check.
- Overlap with the Block C analysis set is at least 880 of 1,037 and is to be
  reported exactly, not assumed.

## 7. What a verdict licenses

**MB2-ALIGNED licenses one sentence**, in section 4 or the Discussion: that the
state described here corresponds to the MYC-associated bicluster rather than to
the proliferation-driven one whose prognostic behaviour was previously reported.
That cites the companion paper as a positive rather than as a rival, and it
explains the ER-negative skew of the buffered cell without special pleading.

**MB1-ALIGNED is the informative failure.** It would mean this arm's axis is the
programme Menegollo already characterised, the published description of MB1 as
Myc-independent would be in tension with our data, and limb (i) would have to
cite them as prior work rather than as a contrast. That discrepancy would be
worth understanding rather than working around.

**Neither licenses anything about outcome.** No survival or prognostic claim
follows from a correlation between two exposure axes.

## 8. What this cannot do

It cannot test whether MB2 carries prognostic information, because Menegollo
report none and we compute none here. It cannot establish that the two
programmes are biologically distinct, which is their result and not ours. It is
an alignment check on one axis, run on 1,037 TCGA samples, and it is
cross-sectional and descriptive.
