---
date: 2026-10-01
status: >
  DECLARATION - pre-data. Nothing below has been computed. This declares a
  SNAPSHOT REBUILD and the GATE it must pass before anything is read from it.
  No outcome variable is read, and none is joined.
repo: myc_human_exploratory. `myc_human_validation` and `myc_mouse` read-only.
relates-to:
  - docs/2026-10-01_E36_declaration.md and docs/2026-10-01_e36_data.md (the stop
    that produced this, and what opening the inputs established)
  - myc_human_validation/scripts/01_fetch_tcga_expression.R (read-only; the
    pipeline this reproduces, and the one line it changes)
  - docs/2026-10-01_e34_result.md (the analysis the gate re-derives)
decides:
  - That the rebuilt snapshot sits BESIDE data/from_validation/, never replacing
    it, and is read by E36 and by nothing else unless separately declared.
  - The gate: E34's headline quantities re-derived on the rebuilt snapshot,
    tumours only, with the thresholds for "material movement" fixed HERE,
    before any rebuilt value exists.
next-action: commit alone, then the E36 amendment, then the data note, then the
  script unrun, then the gate result.
---

# E37 — rebuild the expression layer with the normals, and gate it on E34

## 0. Why this exists

`E36` stopped because there are no solid-tissue-normal samples in any saved
object in this repo. They exist — 113 of them, in the cached Xena source in the
frozen repo — and they were dropped at the first step of the upstream pipeline
by `keep <- which(sty == "01")`, **for a documented reason that does not apply
to E36**: primary tumours only, so that expression and copy number stay aligned.
E36 reads no copy number. **This is restoring samples dropped for an unrelated
constraint, not reversing a quality judgement.**

`E36` cannot proceed without this, and this analysis is the only one in the
queue that can come back against the model.

## 1. What is rebuilt, and what is not

**The expression layer, over tumours AND normals jointly**, following
`myc_human_validation/scripts/01_fetch_tcga_expression.R` as closely as the
frozen repo permits, with **one deliberate change**: the sample filter admits
sample type `11` as well as `01`. Type `06` (metastatic, n = 7) stays excluded,
as upstream.

| step | treatment |
|---|---|
| source | the cached Xena matrix in the frozen repo, **read-only** |
| sample filter | `01` and `11`; `06` excluded. One column per patient per type |
| transform | invert the Xena log2, round to counts, exactly as upstream |
| genes | protein-coding, symbol-keyed, deterministic collapse, as upstream |
| normalisation | **DESeq2 size factors over the JOINT set**, which is the whole point |
| VST | over the joint set |
| scoring | GSVA and mitoPPS **once, over the joint set**, because both are cohort-relative |

**Nothing is written to `myc_human_validation`.** The source is read with
`gzip -dc`; the script is read with `git show`.

## 2. The snapshot sits BESIDE, and is scoped

**It does not replace `data/from_validation/`.** Replacing it would invalidate
E11, E33, E34 and E35 at a stroke and leave their figures unreproducible.

The rebuilt objects get **their own name, their own README, and a one-line
statement of what they are for**. `CLAUDE.md`'s scale-discipline section names
two normalisations in one repo as a trap; **the discipline that resolves it is
that every script states which snapshot it reads, and no result mixes the two.**

- **The rebuilt snapshot is read by `E36` and by nothing else unless separately
  declared.**
- Every script that reads it says so in a comment at the top of its scoring
  block, as `CLAUDE.md` already requires for scale.
- **No figure, table or sentence combines a quantity from one snapshot with a
  quantity from the other.**

## 3. THE GATE — E34 re-derived, tumours only

Before anything is read from the rebuilt snapshot for E36, **E34's headline
quantities are recomputed on it, restricted to the 1,095 tumours**, and compared
with the committed values in `docs/2026-10-01_e34_result.md`.

This is free, and it is valuable in both directions: if it lands in the same
place, the rebuild is sound **and E34 is robust to normalisation**; if it moves,
that is a finding about E34 and it is wanted now rather than from a reviewer.

### 3.1 What is compared

1. **Quadrant assignment**, GSVA and mitoPPS, against E34's.
2. **The declared assertions**: Q3 == `STATE` level 2, Q4 == levels 3 + 4, on
   both instruments. These must hold; failure is a stop regardless of anything
   else.
3. **The four contrasts on the configuration composite**, m1, pooled, both
   instruments: `Q2 - Q1`, `Q4 - Q3`, `Q4 - Q2`, `Q3 - Q1`.
4. **E34's per-cell reading labels** in the pooled cells (m1 and m3, both
   instruments).

### 3.2 The thresholds, fixed now

**PASS** requires all four:

| | threshold |
|---|---|
| quadrant agreement with E34, each instrument | **>= 95%** of tumours with a quadrant in both |
| the two OXPHOS contrasts | each within **0.10** of E34's value |
| the two MYC contrasts | each within **0.10** of E34's value |
| E34's pooled per-cell labels | **unchanged** in all four pooled cells |

0.10 is double `E16`'s `R1_MAX_DELTA` of 0.05, which this repo has used for
"materially different" since 2026-09-04, because a re-normalisation is a larger
perturbation than a ruler swap. Against OXPHOS contrasts of +0.85 to +0.98 it is
about a tenth; against MYC contrasts of 0.00 to 0.23 it is lenient in relative
terms, deliberately, because **what the gate protects is E34's conclusion — that
the MYC contrasts are several times smaller than the OXPHOS ones — and not its
third decimal.**

### 3.3 What happens on each outcome

| outcome | action |
|---|---|
| **PASS** | the rebuild is sound and E34 is robust to normalisation. Record both. **E36 proceeds.** |
| **FAIL on the assertions** | the quadrant is not what it is taken to be on the new scale. **Stop. E36 does not proceed.** |
| **FAIL on any threshold** | **STOP, and deal with that first.** It is a finding about E34, not a nuisance. The result note says what moved and by how much, and E36 waits |

**A FAIL is not a reason to adjust the thresholds.** They are fixed here,
before any rebuilt value exists, precisely so that a disappointing gate cannot
be renegotiated into a passing one.

## 4. What this cannot do

1. **It does not make the two snapshots comparable.** Values from the rebuilt
   snapshot are never compared numerically with values from
   `data/from_validation/` except in the gate itself, which is a deliberate
   like-for-like comparison of the *same* quantities and is reported as such.
2. **It does not re-derive the genomics.** Copy number, GISTIC, mutations and
   the frozen `STATE` inputs other than the scores are untouched. The frozen
   constructor is re-applied to new scores; it is not rebuilt.
3. **It is one cohort and one release.** The Xena matrix is as cached on
   2026-08-28 and is not re-downloaded.
4. **It establishes nothing biological.** It is infrastructure plus one
   robustness check.

## 5. What a PASS licenses

**One sentence about E34**: that its headline contrasts survive a change of
normalisation that moves every input value, which is a robustness check it has
not previously had.

**And it unblocks E36**, which is the point.

**It licenses nothing biological on its own**, nothing about killing, and
nothing about treatment selection.

## 6. What is already known, disclosed here

The E36 data note (`4802508`) was committed before this declaration and
established, from column headers and barcodes only:

- **113 normals**, sample type `11`, in the cached source.
- **All 113 from patients already in the 1,095-tumour set** — the design is
  inherently paired.
- The upstream exclusion rule, and its stated reason.

**No expression value, score or estimate has been read.** The gate's quantities
do not exist yet, and E34's committed values are the fixed comparator.
