---
date: 2026-10-01
status: >
  DECLARATION - pre-data. Nothing below has been computed. Exposure-side only;
  no outcome variable is read, and none is joined.
repo: myc_human_exploratory. `myc_human_validation` and `myc_mouse` read-only.
relates-to:
  - E11 (the conditioning asymmetry), E14 (compartment specificity)
  - E34 (OXPHOS-ordered in 16 of 18 reading cells; the MYC contrast four to
    seven times smaller and identical at both respiratory levels)
  - myc_mouse, the gate model (WT gland slope -0.012, p 0.95; the coupling
    appears only in Myc+ animals)
  - 2026-08-27_human_validation_plan.md (the Hannon pre-malignant series, the
    dataset this question actually wants, not in hand)
decides:
  - The exposure, the stratifier, the three possible readings and what each
    would mean, all fixed before any estimate exists.
  - That all three readings are reportable, and which one is expected.
next-action: commit alone, then data note, then script unrun, then result.
---

# E35 — is the OXPHOS-to-configuration coupling graded by genomic burden?

## 0. The question, and what it is NOT

**The tension this addresses.** In the mouse pre-tumour gland the respiratory
state and the apoptotic configuration are unrelated in wild-type animals (slope
−0.012, p 0.95) and coupled only where MYC is expressed. In human tumours E34
found the coupling present irrespective of MYC level, in both OXPHOS-high
quadrants equally (TCGA Q2 − Q1 +0.981 against Q4 − Q3 +0.845).

The leading reconciliation is that **the mouse contrasts a pre-tumour gland with
transformed tissue whilst the human contrasts tumours with tumours**, so the
human comparison may lie entirely above the condition the mouse manipulates. If
the condition is oncogenic or replicative drive rather than MYC specifically,
every tumour has some version of it and no tumour-versus-tumour contrast can
reveal it.

**This tests that from inside the tumour cohort**, using genomic burden as a
graded proxy for drive, without needing normal tissue.

**It is NOT a test of killing.** There is no death readout in TCGA. The endpoint
is the **coupling between respiratory state and apoptotic configuration**, which
is a transcript relationship. Killing is the mouse's variable and stays there.
Any sentence written from this must say "coupling", never "killing" and never
"lethality".

## 1. Exposure, stratifier, endpoint

| | |
|---|---|
| **endpoint** | the configuration composite, as E11 and E34 define it |
| **exposure** | OXPHOS status, the quadrant's own call. The coupling is the **Q2 − Q1 and Q4 − Q3 contrasts**, exactly as E34 computes them |
| **stratifier** | genomic burden, tertiles, from `tcga_brca_covariates.rds` |

**Burden proxies, all three reported, none promoted after the fact:**
`Aneuploidy.Score` (n = 1038), `FGA` (n = 1065), `genome_doublings` (n = 1020).
`Aneuploidy.Score` is primary because it has the clearest interpretation as
chromosomal instability; the other two are replications of the same question on
the same cohort, not independent evidence.

**TCGA only.** SCAN-B carries no CNV, so there is no replication cohort for
this. A result here is **single-cohort and must be written that way.**

## 2. The three readings, all reportable, fixed now

| reading | condition | meaning |
|---|---|---|
| **GRADED** | the OXPHOS contrasts increase monotonically across burden tertiles, with the low and high tertile CIs non-overlapping | the condition is dose-like and visible across human tumours. The strongest outcome, and it would reinstate a stratifier |
| **SATURATED** | the OXPHOS contrasts are present in all three tertiles and do not differ | every tumour lies above the condition. Supports the universal reading: the caution applies to the whole population the drug class targets |
| **ABSENT AT LOW BURDEN** | the OXPHOS contrasts do not clear zero in the lowest tertile but do in the highest | the condition is genuinely graded and some tumours lack it. **Strongest of the three** |

**SATURATED is the expected outcome** and it is recorded here so that it does not
read as a disappointment when it arrives. It is also the reading that most
directly supports the sentence the Discussion already wants.

## 3. Stratify, do not adjust

`Aneuploidy.Score` correlates with proliferation and with basal subtype, so the
circularity this analysis exists to escape partly follows it in. Therefore:

- **Tertile stratification, not a continuous interaction term.** No product term
  anywhere; that estimand has failed in Block B, Block C, Block G, H4, N1 and N2.
- **Within-subtype estimates are the result, not a robustness check**, as in E34.
  Report Luminal and Basal separately wherever cells allow, with a declared floor
  of n >= 20 per quadrant cell.
- **The proliferation companion is mandatory**, as in E34: the primary repeated
  with `PROLIF_DISJOINT` added, reported beside it, never chosen between.
- Purity and leukocyte fraction as in E11 and E34.

## 4. What this cannot do, stated before it runs

1. **Burden is not drive.** Aneuploidy, FGA and genome doublings index
   chromosomal instability, which is as much a consequence of drive as a measure
   of it. These are proxies and the note must call them that.
2. **Configuration is not lethality.** E34's ceiling applies unchanged. The
   sibling study's functional BH3 assays already came back null against a
   transcript-level expectation, so a transcript coupling and a functional
   liability can come apart.
3. **One cohort, no replication available.**
4. **It is not the right dataset.** The Hannon pre-malignant series is the
   graded-drive axis this question wants and it is not in hand. This is the
   best available substitute, not the intended test.

## 5. What each reading licenses

**SATURATED** licenses one clause: that the coupling is present across the range
of genomic burden represented in this cohort, consistent with every established
tumour lying above the condition. It supports the universal drug-class caution
and it reinstates no stratifier.

**GRADED or ABSENT AT LOW BURDEN** licenses a sentence that the coupling
strengthens with genomic burden, single-cohort, with the proxy caveat attached,
and it would make a burden-stratified therapeutic caution arguable rather than
universal.

**No reading licenses a treatment-selection claim**, and none licenses any
statement about killing. No cohort here carries an OXPHOS-directed intervention
or a death readout.

## 6. The companion analysis this belongs with

Genomic burden gives gradation **within** tumours. The TCGA solid-tissue-normal
comparison (~113 samples) gives the anchor **below** them, and is the human
analogue of the mouse wild-type gland. Together they are the human version of
the mouse's gland-to-tumour axis; separately each is half the argument. That
analysis needs its own declaration and carries its own confounder, since
normal-adjacent breast is largely stroma and adipose.
