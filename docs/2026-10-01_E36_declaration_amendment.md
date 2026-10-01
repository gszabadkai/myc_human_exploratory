---
date: 2026-10-01
status: >
  AMENDMENT to docs/2026-10-01_E36_declaration.md, PRE-DATA. Nothing has been
  computed. It amends the design for a constraint the original declaration did
  not anticipate, fixes the epithelial-content estimate by gene list, and
  records a disclosure about the positive control. Committed ALONE, before any
  rebuilt value or any fit exists.
amends: docs/2026-10-01_E36_declaration.md (committed d71fa45)
relates-to:
  - docs/2026-10-01_e36_data.md (the stop, and what established the pairing)
  - docs/2026-10-01_E37_declaration.md (the rebuild E36 waits on, and its gate)
decides:
  - Q2 becomes a WITHIN-PATIENT contrast; Q1 keeps its form but its
    normal-against-tumour comparison is bootstrapped by PATIENT.
  - The epithelial-content estimate, by gene list, before any fit.
  - That the positive control is pre-specified but NOT blind, and why.
---

# E36, amended: the comparison is paired

**Nothing here has been computed.** The original declaration stands except where
this amends it, and it is not rewritten — it is a dated record of what was
declared before the inputs were opened.

---

## 1. Why this amendment exists

The E36 data note established, from barcodes alone, that **all 113
solid-tissue-normal samples come from patients already in the 1,095-tumour
analysis set — 113 of 113.** The original declaration treats normal and tumour
as two groups. **They are not: they are 113 matched pairs plus 982 unmatched
tumours.** Treating them as independent would be wrong for every normal in the
set.

This is a better design than the one declared, not a worse one. A within-patient
contrast holds genotype, age, germline and collection constant, which no
across-person comparison can.

## 2. Q2 becomes a within-patient contrast

**Original (declaration section 4):** MYC activity distributions of normals and
of each tumour quadrant, as medians with IQR and the percentage of tumours in
each quadrant exceeding the normal median.

**Amended.** The distributional form is **kept and still reported**, because it
is what the premise is stated in. **Added, and primary for Q2:**

> **the within-patient difference, tumour minus that patient's own normal, on
> each of the three estimators (`M_a`, `M_b`, `FELSHER__PROLIFSTRIP`), over the
> 113 pairs.**

Reported as median and IQR of the paired difference, and as **the percentage of
the 113 pairs in which the tumour exceeds its own normal**. This is a stronger
test of the premise than comparing distributions across different people: if
tumour MYC does not exceed its own matched normal in a substantial fraction of
pairs, the threshold reconciliation fails, and the original declaration's
instruction to say so plainly applies unchanged.

**Still no verdict on Q2.** It remains descriptive.

## 3. Q1 keeps its form; its comparison is bootstrapped by patient

**Q1 is unchanged in form.** It compares two correlations computed *within*
tissues — `rho(OXPHOS, configuration | epithelial content)` in normals, and the
same in tumours — and a correlation computed within a tissue is not affected by
the pairing.

**What the pairing changes is the contrast between them.** The 113 patients
contribute to both correlations, so the two are not independent and no interval
on their difference may assume they are.

> **Any direct contrast of the normal and tumour correlations uses a
> PATIENT-LEVEL bootstrap: resample the 1,095 patients with replacement, carry
> each drawn patient's tumour and, where it exists, that patient's normal
> together into the resample, and recompute both correlations and their
> difference inside each one.**

**This is a change of resampling unit, not a redesign.** 10,000 resamples,
percentile intervals, as E32 and E33 used. The two correlations are also
reported with their own within-tissue intervals, as the original declaration
asks.

## 4. Epithelial content, fixed by gene list

The original declaration requires an estimate computable in both tissues and
forbids ABSOLUTE purity, which is tumour-only. **Declared now, before any fit:**

> **`EPITHELIAL` = mean per-gene z, over the joint snapshot, of:**
> **luminal epithelial** `EPCAM`, `CDH1`, `KRT8`, `KRT18`, `KRT19`;
> **basal / myoepithelial** `KRT5`, `KRT14`, `KRT17`, `TP63`.

**Nine genes, all present in the matrix** (checked; 5/5 and 4/4).

**Why both compartments.** Normal breast epithelium carries a luminal and a
basal/myoepithelial layer; most tumours are luminal. **A luminal-only marker
score would be part epithelial-content estimate and part subtype score**, and
would then adjust away some of the thing being measured. Including the
myoepithelial markers makes it a tissue-composition estimate rather than a
lineage one.

**No deconvolution reference is used.** Every such reference available here is
built on tumours and does not transfer to normal tissue, which is the same
objection that disqualifies ABSOLUTE.

**Reported with and without the adjustment**, as the original declaration
requires.

**A descriptive companion, not an adjustment**: the mean z of adipose
(`ADIPOQ`, `FABP4`, `PLIN1`, `LEP`, `CIDEC`) and fibroblast (`COL1A1`,
`COL1A2`, `FAP`, `PDGFRB`, `THY1`) markers, all present, reported **only** to
show that `EPITHELIAL` behaves as a composition estimate should — strongly
higher in tumours, and strongly anti-correlated with the stromal score. **It
enters no model.**

## 5. The positive control is pre-specified but NOT blind

The E36 data note named **Fe-S cluster assembly** as the declared primary
control, with **mitophagy** and **mitophagy PINK1/PRKN only** alongside, all
three before any fit. That protection stands and is the main thing.

**But these are not blind comparators, and that is recorded here rather than
left implicit.** All three are from E14, whose numbers are committed and have
been read. For Fe-S cluster assembly, the cytosolic half against OXPHOS is
**+0.118 in TCGA** (proliferation-adjusted; +0.152 raw) and **+0.159 in
SCAN-B**. The mitophagy halves are likewise on record.

**So the control is pre-specified, not naive.** What the advance naming buys is
that **one of the three cannot be promoted after the fact**; what it does not
buy is ignorance of how they behave in tumours. The result note must say this in
the same breath as the control's verdict.

## 6. Two checks before the figure is designed

Both descriptive, both to be reported whatever they show, because they decide
whether a per-quadrant comparison is a comparison or a description:

1. **How many of the 113 fall in each tumour quadrant** — strictly, how many of
   the 113 *patients'* tumours do. **If Q2 or Q4 holds a handful, the
   per-quadrant comparison is descriptive at best**, and the note says so rather
   than reporting four cells as if they were comparable.
2. **Whether the 113-pair subset is skewed by subtype or by collection site**
   (PAM50, and the TSS code in barcode positions 6-7), against the full
   1,095-tumour set.

## 7. One thing the stop established on its own

**The normals were excluded upstream for copy-number alignment, and E36 reads no
copy number.** The rule was `keep <- which(sty == "01")`, with the stated reason
*"so the expression and copy-number analyses stay aligned"*. **The rebuild is
therefore restoring samples dropped for an unrelated constraint, not reversing a
quality judgement about those samples.** That belongs beside any sentence about
where the normals came from.

## 8. Unchanged

Everything else in `docs/2026-10-01_E36_declaration.md` stands: the four
readings including **PRESENT IN NORMAL**, the positive control as a **stop
condition**, reporting the normal n at every step, and the limits in its
section 5 — field-adjacent tissue, n ≈ 113, unreplicated, transcript-level, and
the configuration not being lethality. **No reading licenses a
treatment-selection claim, and none concerns killing.**

**E36 does not start until E37's gate passes.**
