# The orthotopic series — what it carries and how to score it

> # VOID — SAMPLE IDENTITY CORRECTION, 2026-09-10
>
> **Do not quote any part of this document.** The arm this plan treats as the
> respiratory manipulation **alone** was a **double-construct** arm — the
> respiratory manipulation **plus an anti-apoptotic transgene** — so its
> `Bcl2l1` readout was **construct transcript, not endogenous**.
>
> **What that voids, specifically:**
>
> - **§1(a), the licensing relationship.** "PGC1α overexpression alone raised
>   `Bcl2l1` 2.7× against EV … with no Bcl-xL construct present", and its claim
>   to be "the single cleanest argument" that the human `BCL2L1` coefficient is
>   causal. **Void.** In the **construct-free** arm of the same cohort,
>   endogenous `Bcl2l1` is **-0.113 log2 [-0.283, +0.174]** — null — with `Mcl1`
>   null beside it (-0.100 [-0.229, +0.047]).
> - **§1(b), the MYC ceiling as an in vivo demonstration.** It rests on the same
>   mislabelled arm. **Void as an intervention result.**
> - **§1(c), the abstract sentence and the gland-versus-human comparison.**
>   **Void**, and see synthesis 8.4 traps 3 and 13 before any version of it is
>   rewritten.
>
> **The reverse direction is void too, and it is the sharper half.** In the same
> cohort an anti-apoptotic construct moved **no respiratory ruler** — `ox_lvl`
> **+0.093 [-0.298, +0.581]**, `ox_rel` **-0.023**, `ox_mt` **-0.246**, all
> covering zero. So **"raising the guardian raises respiration" has no support
> either**, and **no claim that mouse tumours reproduce the human configuration
> survives**.
>
> **No human number in this repo changes.** What is void is what these mouse
> results were taken to *support*, never what the human data say.
>
> The body below is **preserved unaltered** as the record of what was planned and
> why. Read `docs/2026-09-10_mouse_corrections.md` first; the corrected mouse
> sources live in `myc_mouse` and are cited there by filename, never copied here.

**2026-09-07. Carry-over note. Written in chat, to be executed in Claude Code.**
Repo: `myc_mouse`, `experimental-cohorts` branch off `paper-final`. Numbered
scripts continue `paper-final`'s sequence.

**N3 throughout.** These are transcript associations. "Primed" is never written
of a transcript.

---

## 1. Why this dataset matters, stated precisely

It is the last piece of the causal spine, and it does three jobs no other
dataset in the paper can do.

**It cannot do a fourth, and that must be said first.** Every arm is MYAZ —
MMTV-Myc tumour derived, MYC-high throughout. There is no MYC-low arm, so this
dataset **cannot estimate the MYC × OXPHOS interaction**. The interaction lives
in the iMMEC rtTA-MYC ±dox × ±PGC1α design, on the death readout. Do not ask
this dataset for it.

What it does carry:

**(a) The licensing relationship by intervention.** In the preliminary look,
PGC1α overexpression alone raised `Bcl2l1` 2.7× against EV — 193 vs 71 median
CPM — with no Bcl-xL construct present. That is ρ(BCL2L1, OXPHOS) > 0 produced
by manipulating respiration, not observed by correlating it. It is the single
cleanest argument that the human β₂ on `BCL2L1` (+0.414 TCGA, +0.344 SCAN-B, at
fixed MYC, proliferation-robust) reflects a causal link rather than a
cross-sectional association.

**(b) The selection consequence — the MYC ceiling.** No PGC1α tumour exceeded
2000 `Myc` CPM; EV had two above 2600 and BclxL two above 2500. `Bcl2l1` in the
PGC1α arm converged into 181–227, a 1.25× spread against EV's 1.6× and BclxL's
1.8×. Read as a ceiling and a floor rather than shifted means, that is the in
vivo demonstration that M and O cannot both be high — the collider argument the
human arm depends on, shown by intervention.

**(c) A sentence already in the abstract that has never been tested.**

> *"Mouse tumours resemble human tumours rather than the gland they arose from,
> placing this reversal at the transition into malignancy rather than between
> species."*

The reversal is on the guardian pair. In the gland, `ox_rel → Mcl1:Bcl2l1` is
**+0.351** (mouse gate model, `docs/2026-09-02`, §3.4 within-timepoint, the
document's own cross-sectional transfer estimate). In human tumours the same
endpoint on the same ruler is **−0.31 to −0.46** (E23, at fixed MYC). These are
mouse *tumours*. Which side do they fall on? **This is the highest-value item in
the dataset and it can break a sentence that is currently written.**

---

## 2. What is known before any script runs

67 samples, 10 groups, `salmon_merged_gene_counts.tsv`. Median CPM from the
preliminary pass:

| group | n | `Bcl2l1` | `Ppargc1a` | `Myc` |
|---|---|---|---|---|
| EV | 7 | 71 | 0.0 | 1915 |
| BclxL | 7 | 326 | 0.3 | 2033 |
| Pgc1a | 7 | 193 | 302 | 1440 |
| NTEV | 7 | 82 | 0.2 | 2208 |
| KOEV | 7 | 62 | 0.1 | 2368 |
| NTPgc1a | 7 | 79 | 129 | 2286 |
| KOPgc1a | 7 | 79 | 0.2 | 2176 |
| MYAZ | 6 | 59 | 0.0 | 1370 |
| FaMY1 | 6 | 59 | 0.0 | 1280 |
| BoMY | 6 | 51 | 0.0 | 1293 |

**The p21-KO series (NT/KO × EV/Pgc1a) is set aside** as unrelated to the main
narrative. `Cdkn1a` runs 86–140 CPM across the KO groups, so the knockout is not
visible at transcript level — consistent with a CRISPR indel surviving NMD.
`KOPgc1a` has lost the transgene entirely (`Ppargc1a` 129 → 0.2), which is
interesting and is not this paper's business.

**MYAZ / FaMY1 / BoMY are unidentified.** Confirm what they are before scoring —
if they are primary vs fat-pad vs bone-derived lines they may be a second, quite
different comparison, and they must not be silently pooled into the vector
series.

---

## 3. GATE 1 — input audit, before any script is written

1. **Are salmon transcript-level quant files on disk, or only the merged gene
   matrix?** The merged matrix goes through `DESeqDataSetFromMatrix`, which
   silently discards the transcript-length offset that the `tximport` pathway
   carries. If the quant files exist, use `tximport`. If not, say so in the note
   and proceed — but say so.
2. **What are MYAZ, FaMY1, BoMY?** Names, provenance, whether they belong in the
   same cohort as the vector series.
3. **Which samples define the cohort?** This is not bookkeeping. `ox_rel` is a
   relative score and GSVA is cohort-relative, so including or excluding the p21
   series and the three cell lines changes every score. **Fix the scoring set in
   advance and state it.** Recommended: the vector series alone (EV / BclxL /
   Pgc1a, n = 21) as primary, with the wider set as a named sensitivity — but
   decide it before the numbers, not after.
4. **Composition.** All arms are established tumours in the same site, which is
   why this is the right system — but the fat-pad failure's transferable rule
   applies: *before testing any trend, check whether the composition axis and the
   biological axis are separable at all.* Run the loading table (epithelial,
   stromal, adipose, immune, proliferation markers against `ox_rel`, by group)
   and a positive control **first**. If the group contrast is confounded with
   cellularity, that changes what can be claimed.

Report and wait.

---

## 4. Standing rules that apply

- **Species = cohort.** Never pool these scores with the 6W/12W timeline or with
  human. mitoPPS values are not numerically comparable across cohorts.
- **`ox_rel`**: 88 nuclear OXPHOS subunits minus the rest of MitoCarta, 13
  mtDNA-encoded genes in the denominator.
- **Input object separation**: GSVA takes log-scale input (`kcdf="Gaussian"`);
  mitoPPS takes linear DESeq2-normalised counts. These two never share an input
  object.
- **No per-gene FDR** across the twelve. Report estimates with intervals.
- Option A workflow: numbered scripts written by Claude Code, sourced manually in
  Positron. `if (FALSE)` sandbox block in every script.

---

## 5. Pre-specify these three claims before fitting

**C1 — the licensing relationship.** `Bcl2l1` is higher in Pgc1a than in EV, in
the absence of a Bcl-xL construct. Directional, one-sided, n = 7 vs 7. Report
`Ppargc1a` and the OXPHOS score alongside to show the manipulation worked at the
level it claims to.

**C2 — the MYC ceiling.** This is a claim about the **upper tail**, not the
mean, and it must be tested as one. A shifted median and a ceiling are different
statements; the medians alone do not distinguish them. Pre-specify the statistic
(upper quantile contrast, or max, with a permutation null) and be explicit that
n = 7 per group makes this suggestive at best. **A ceiling claim at n = 7 is a
figure and an observation, not a test** — write it that way.

**C3 — which side of the reversal.** Fit `Mcl1:Bcl2l1` and its two members
against `ox_rel` across the vector series. Predicted **negative** if orthotopic
tumours behave like human tumours; **positive** if they behave like the gland.
Report both members separately, as E23 did, so it is visible which one carries
it. Pre-specify both readings as informative — this is the rare case where either
direction is a result.

**Failure readings, fixed now.** C1 or C3 reversing under the composition
adjustment from Gate 1 means the group contrast is confounded with cellularity
and the dataset contributes an exclusion, not a result. Write that as a negative
if it happens, in the fat-pad note's style.

---

## 6. What this will and will not license

**Will:** the causal version of the human β₂ finding; the in vivo selection
signature; and — depending on C3 — either support for or the withdrawal of the
abstract's "mouse tumours resemble human tumours" sentence.

**Will not:** anything about the MYC × OXPHOS interaction, anything about
apoptotic priming (N3), and nothing that pools with the human cohorts.

---

## 7. Deliverables

One numbered script continuing `paper-final`'s sequence, plus a dated note in
`docs/` in the house pattern, stating which manuscript sentence it licenses and
which it does not. The branch merges back into `paper-final` regardless of
outcome.
