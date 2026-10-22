---
date: 2026-10-22
status: Provenance reconstruction — no new analysis, no new claim
subject: Membership of the 6-gene signed configuration composite
relates-to:
  - docs/2026-09-03_human_arm_for_mouse_reconciliation.md (sections 2.5, 3.1, 5.1)
  - docs/2026-08-27_human_validation_plan.md (section 9, the ER confound)
  - 2026-10-02_E41_metabric_amendment_and_prompt.md (section A1)
  - 2026-10-02_human_arm_handoff.md (retired claims; limitations 4 and 8)
purpose: Section 4 narrative — record why the composite has these six genes
  and not the other six, before the section is finalised.
---

# The 6-gene signed configuration composite: where it came from, gene by gene

Reconstructed 2026-10-22 from project documents and prior sessions. **Nothing
here is a new result.** It records reasoning that is currently distributed
across several documents and nowhere stated in one place.

---

## 1. The object

**The 6-gene signed configuration composite:**

```
BBC3  BID  BIK  BAD  BCL2L1   entering POSITIVE
MCL1                          entering NEGATIVE
```

This is what scripts **E34** and **E36** compute.

- **E34** was the quadrant analysis: tumours split at the cohort median of MYC
  activity and of OXPHOS status into four groups, establishing that the
  configuration separates by respiratory state and not by MYC.
- **E36** was the tumour-versus-normal comparison against 113 TCGA solid-tissue
  normals.

**E10's twelve is a different, unsigned object** and is not used by either.

### 1.1 The only place the membership is written down — and it is a correction

Section **A1** of the 2 October METABRIC amendment records that the declaration
had been saying "the 12-gene configuration composite" throughout, including in
section 5, and that this was wrong: E34 and E36 use the signed six.

That entry **fixes the name. It does not argue the choice.** The rationale below
is reconstructed from the table the six were drawn off.

**This is the gap this note exists to close.**

---

## 2. The parent set — the twelve

The twelve are the author's own curated lists: **7 pro-apoptotic x 5
anti-apoptotic = 35 ratios**, no self-pair, no empty cell.

| side | genes |
|---|---|
| pro-apoptotic | `BCL2L11` `BMF` `PMAIP1` `BBC3` `BID` `BAD` `BIK` |
| anti-apoptotic | `BCL2` `BCL2L1` `MCL1` `BCL2L2` `BCL2A1` |

Two properties of the twelve bear on the selection:

1. **All twelve clear the 25th expression percentile in both cohorts.** Nothing
   in this set is a low-expression artefact. (`HRK` and `FASLG`, in the wider
   44-gene machinery set, do not clear it and are flagged everywhere.)
2. **Eleven of the twelve are MitoCarta outer-mitochondrial-membrane proteins.
   Only `BMF` is not.**

One curation correction on the record: **`BCL2L2` was supplied on both lists**
and was corrected on 2026-09-02 to anti-apoptotic only. Bcl-w is a multidomain
guardian, not a BH3-only sensitiser.

---

## 3. The table the six were cut from

Partial Spearman against the OXPHOS axis, adjusted for proliferation
(`PROLIF_DISJOINT`, the 318-gene E2F/G2M index built to share zero genes with
either axis), cross-cohort mean of TCGA and SCAN-B.

Source: `machinery_and_priming.rds$component_cor`; figure `E10_fig5`; tabulated
in reconciliation section 3.1.

| gene | side | OXPHOS TCGA | OXPHOS SCAN-B | OXPHOS mean | MYC mean | in the six? |
|---|---|---|---|---|---|---|
| **BAD** | pro | +0.525 | +0.482 | **+0.503** | +0.138 | **yes, +** |
| **BCL2L1** | anti | +0.426 | +0.350 | **+0.388** | +0.002 | **yes, +** |
| **BIK** | pro | +0.338 | +0.213 | **+0.276** | +0.029 | **yes, +** |
| **BBC3** | pro | +0.365 | +0.170 | **+0.267** | +0.136 | **yes, +** |
| **BID** | pro | +0.329 | +0.091 | **+0.210** | **+0.300** | **yes, +** |
| BCL2A1 | anti | -0.004 | -0.107 | -0.056 | +0.055 | no |
| BCL2L2 | anti | -0.145 | -0.008 | -0.077 | -0.076 | no |
| PMAIP1 | pro | -0.147 | -0.064 | -0.105 | **-0.209** | no |
| BCL2 | anti | -0.160 | -0.074 | -0.117 | -0.191 | no |
| BCL2L11 | pro | -0.314 | -0.025 | -0.170 | -0.101 | no |
| BMF | pro | -0.200 | -0.223 | -0.212 | -0.115 | no |
| **MCL1** | anti | -0.241 | -0.292 | **-0.266** | +0.022 | **yes, -** |

### 3.1 The implicit rule

**Take the genes at the two ends of the OXPHOS column, and give each the sign it
actually shows.** The composite is signed so that it rises when a tumour carries
the configuration, which is why MCL1 — the most negative gene — enters with a
minus.

A second criterion from reconciliation section 3.1 does most of the work.
**"Individually OXPHOS-led"** means `|rho(OXPHOS) - rho(MYC)| >= 0.20` in at
least one cohort. Exactly five genes meet it:

```
BCL2L1  +0.455 / +0.317
BAD     +0.363 / +0.368
BIK     +0.330 / +0.163
BBC3    +0.202 / +0.061
MCL1    -0.321 / -0.254   (leaning hard to OXPHOS while running the other way)
```

**The other seven are not.**

So: **the composite is the five OXPHOS-led genes, plus BID.** BID is the
exception and section 4.5 below is about it.

---

## 4. Gene by gene — the five positives and MCL1

### 4.1 BAD, +0.50

The largest positive coefficient in the family. Pro-apoptotic, MOM-resident,
individually OXPHOS-led at +0.363 / +0.368 — in **both** cohorts, which no other
member manages on that criterion. The cleanest member of the set.

### 4.2 BCL2L1 / BCL-XL, +0.39

**The load-bearing gene.** It is the anti-apoptotic member that *rises* with
respiration, which is the entire guardian argument.

- Individually OXPHOS-led, +0.455 / +0.317.
- **Purity is not the explanation.** The E16 check moved it +0.419 -> +0.406
  with purity and leukocyte fraction carried.
- **Subtype is not the explanation.** E19 stratified all twelve genes by
  Luminal / Basal using pooled scores with samples subset (never re-scored
  within stratum, which would produce a different axis). BCL2L1 held in both —
  branch 1 of the four pre-declared branches, the cleanest available reading.
- It carries the copy-number result. **G2: MYC-amplified tumours show excess
  BCL2L1 amplification, adjusted OR 1.77 (1.26-2.47).** Caveat that must travel
  with it: under `FRACTION_GENOME_ALTERED` instead of `ANEUPLOIDY_SCORE` it
  drops to 1.36 (0.984-1.88). Both covariates were pre-specified, so this is a
  disclosed sensitivity rather than a choice — but BCL2L1's pass is
  covariate-dependent and cannot be reported as a bare odds ratio.

**One known non-transfer.** E18 found BCL2L1 is the single gene in the
configuration that does not transfer to cell lines: +0.31 to +0.43 on 8 of 8
tumour-cohort cells, -0.12 to +0.01 on 8 of 8 CCLE cells.

### 4.3 BIK, +0.28

Individually OXPHOS-led at +0.330 / +0.163. In E18 — the CCLE configuration
check — BIK kept its sign on 4 of 4 rulers against both tumour cohorts,
alongside BBC3, MCL1 and BAD.

**Caution to carry:** in the orthotopic mouse tumours `Bik` sits at 12-20 counts
and was explicitly discounted there. It has no usable mouse counterpart.

### 4.4 BBC3 / PUMA, +0.27

The mouse's own trigger gene, so it is present for continuity with the mouse
arm. Individually OXPHOS-led at +0.202 / +0.061 — **the weakest of the five on
that criterion**, carried by TCGA.

It is one of the two genes appearing in the three Basal ratios that survived
every check (`BBC3/BCL2`, `BBC3/MCL1`, `BID/BCL2`), and `BBC3/MCL1` is the only
cell in the 420-cell grid carrying all three marks.

### 4.5 BID, +0.21 — the weak member

**This gene is weak in three independent ways and should be flagged in any
methods statement.**

1. **It is the only one of the six that is not OXPHOS-led.** It is MYC-led:
   +0.300 on MYC against +0.210 on OXPHOS. Reconciliation section 3.1 names BID
   and PMAIP1 as the only two genes whose MYC association survives conditioning
   on OXPHOS in both cohorts under all three estimators — and notes they move in
   opposite directions, so this is two genes, not a shift in balance.
2. **The largest cohort disagreement in the six**: +0.329 TCGA against +0.091
   SCAN-B, a 3.6-fold gap.
3. **It was pre-declared as a negative control.** The 27 August validation plan
   lists the endpoint negatives as BID/BCL2L1, BAX/BCL2L1, BCL2L11/BCL2L1,
   BAK1/BCL2L1, with the rationale: *the mouse says only PUMA/BCL-XL reverses;
   human should show the same.* BID was named in advance as the gene expected
   **not** to behave like PUMA. It then behaved distinctively. That is not
   cherry-picking, but it is the opposite of a confirmation.

**"BID as a mouse-human link" is on the explicit retired-claims list** (5
October), on the ground that no mouse Bid result exists.

### 4.6 MCL1, -0.27, entering negative

The most negative gene in the family, and individually OXPHOS-led while running
the other way (-0.321 / -0.254).

The reconciliation document calls the axis of the whole configuration
**"BCL2L1 up with MCL1 down"** — the two anti-apoptotic genes that move most, in
opposite directions. That is why MCL1 and not one of the other four down-genes
carries the minus: it is the one with a large, OXPHOS-specific coefficient.

G2 copy number: MYC x MCL1 adjusted OR 1.69 (1.17-2.44) at `+2`. Note it
**inverts at `>= +1`** (OR 0.503) and that the association is genuinely
basal-specific (Basal OR 3.12; Breslow-Day by TNBC p 0.022), unlike BCL2L1's,
which is subtype-broad.

---

## 5. Gene by gene — the six exclusions

### 5.1 BMF, -0.21 — compartment

**The strongest exclusion, and it is structural rather than statistical.** BMF
is the only one of the twelve that is **not** a MitoCarta outer-membrane
protein.

Dropping it keeps the composite compartment-clean. This matters because the
compartment argument — nuclear-encoded mitochondrial members rise along the
axis, cytosolic effectors run against it at the 0th percentile of thirty
MitoCarta leaf pathways subjected to the same split — would be undermined by a
signed score mixing the two sides.

Note BMF's coefficient (-0.21) is larger in magnitude than BCL2's or PMAIP1's,
so this exclusion is **not** explicable by effect size. Compartment is the
reason.

### 5.2 BCL2L11 / BIM, -0.17 — inconsistent, and historically tangled

Three grounds:

- Not OXPHOS-led.
- **The cohorts disagree badly**: -0.314 TCGA against -0.025 SCAN-B, essentially
  zero in the larger cohort.
- A tangled record. It was mistakenly relabelled a negative control when it was
  really a co-candidate pre-specified on independent protein evidence; the
  SCAN-B replication returned the effect with the **wrong sign** at n = 3,143;
  and the DepMap breast coefficient (+0.107 mitoPPS, p 0.00018; GSVA +0.040,
  p 0.055) was registered as **unresolved**, failing the two-instrument rule.

Including it would import all of that into the composite.

### 5.3 BCL2, -0.12 — the ER confound, ruled before any fit

**The cleanest documented exclusion.** Section 9 of the 27 August validation
plan, written before any model was fitted:

> BCL2 is a canonical estrogen-responsive gene and ER+ tumours are BCL2-high.
> Any BCL2-family ratio compared across MB2_UF and MB1_UF is partly an ER
> readout, not a MYC readout.

Mitigation 1 of 3: **restrict primary BCL-family measures to BCL2L1 and MCL1,
not BCL2** — with the added note that this also matches the mouse, which is a
BCL-XL story. BCL2 was additionally used as a deliberate negative control (G-d
in the DepMap block, which passed).

### 5.4 PMAIP1 / NOXA, -0.11 — small, and on the wrong axis

Small coefficient, and it is the second of the two MYC-held genes (-0.209 on MYC
against -0.105 on OXPHOS), so it indexes the wrong axis for a composite built on
respiratory state. In the mouse, NOXA was a clean negative at the 7.7th
percentile.

### 5.5 BCL2L2 / BCL-w, -0.08, and BCL2A1, -0.06 — flat

Both are labelled **"flat"** in the reconciliation document's own architecture
diagram. Neither is distinguishable from zero on the OXPHOS axis in either
cohort. Including flat genes only dilutes a mean-of-z composite.

BCL2L2 additionally carries the 2026-09-02 curation correction noted in section
2 above.

---

## 6. The scope problem — the six are NOT a mouse-human similarity selection

The composite is sometimes described as selecting for similarity between the
human and mouse patterns. **The selection does not do that**, and four things
work against reading it that way.

### 6.1 The membership is entirely human-derived

Every criterion in section 3 — the sign, the magnitude, the OXPHOS-versus-MYC
lead — comes from the TCGA / SCAN-B table. The mouse contributes the *motivation*
for looking at BBC3 and BCL2L1 and nothing else.

### 6.2 The mouse has no gene-level counterpart for most of the six

Verified directly against `interaction_sig_genes_p05.csv` on 2026-10-22:

```
rows: 357
BCL2-family / apoptosome / caspase hits: NONE
```

**Zero BCL2-family genes are in the 357-gene mouse interaction set**, and zero
apoptosome or effector-caspase genes. The mouse's gene-level apoptotic result is
`Bbc3` and `Bcl2l1`, with `Mcl1` cross-sectionally; the rest of the mouse cell
death work is pathway-level fGSEA on the Tang modalities.

So **four of the six have no mouse value to be similar to.**

### 6.3 Where both species have a value, they partly disagree

Against the orthotopic PGC1a intervention:

| gene | human, partial rho vs OXPHOS | orthotopic, PGC1a raises OXPHOS | agree? |
|---|---|---|---|
| `BCL2L1` | +0.414 / +0.344 | +1.54 | **yes** |
| `BAD` | +0.50 | +0.18 | yes |
| `BBC3` | +0.412 / +0.200 | **-0.73** | **no — opposite** |
| `MCL1` | -0.211 / -0.222 | -0.02 | weakly |
| `BID` | +0.453 / +0.189 | -0.01 | no |

**The PUMA fall is a mouse-tumour feature that human tumours do not show.**

### 6.4 The estimands are orthogonal — standing rule

Mouse `bMX` is a **MYC x OXPHOS interaction coefficient**. The human value is a
**partial Spearman at fixed MYC**. These are orthogonal by construction:

```
Cov(M.O, O) = Cov(M.O, M) = 0  exactly
```

One carries no information about the other. The standing rule forbids setting a
mouse interaction coefficient beside a human partial Spearman. A similarity
claim built across that boundary is not a claim about biology.

---

## 7. Sign-robustness — a reviewer will ask

From the E38 diagnostics, and already recorded as limitation 4 of the human-arm
handoff:

> **The normal-side zero is not robust to which gene carries the sign.**
> `config, -BIK` gives +0.529; six permutations span -0.276 to +0.529.

So **"MCL1 negative" is not a cosmetic choice** — it is doing a large share of
the work, and the question "why that gene?" has to have an answer ready. Section
4.6 is that answer.

Related, limitation 8: the genes that *fall* (MCL1, BMF, BCL2L11, BCL2, PMAIP1)
are conspicuously lymphocyte-expressed, leukocyte fraction is a coarse
adjustment, and BCL2 is also estrogen-responsive. **The immune alternative to the
compartment result is untested.** Four of those five are excluded from the
composite, which incidentally limits the exposure — but MCL1 is not, and it
carries the minus.

---

## 8. The defensible framing for section 4

**Human-derived, mouse-motivated.**

The six are the genes at the ends of the human OXPHOS column, with BBC3 and
BCL2L1 present because the mouse named them in advance, and the other four
present because the human data put them there.

That is what the analysis actually did, and it is defensible as written.

**If a genuine mouse-human similarity statement is wanted**, only BBC3 and
BCL2L1 can carry it — and BBC3 points the wrong way in the orthotopic tumours,
which leaves BCL2L1 alone.

---

## 9. Open items

1. **No document argues the six-gene membership.** A1 records it; this note
   reconstructs it. If the reconstruction above is accepted, it should be
   committed as the definition of record and A1 amended to point at it.
2. **No statement of why the five down-genes reduce to MCL1 alone.** Section 4.6
   gives the answer (largest magnitude, OXPHOS-led, forms the axis with BCL2L1);
   it has not previously been written down.
3. **BID's status inside the composite is unreconciled** with its pre-declared
   role as a negative control. Either the methods state that a pre-declared
   negative was included on human evidence, or BID comes out and the composite
   becomes five genes. **This is a decision, not an analysis.**
4. **The sign-permutation span (section 7) is not currently in the section 4
   text** and probably should be, given that the composite's headline property
   depends on it.
