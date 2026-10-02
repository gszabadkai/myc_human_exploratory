---
date: 2026-10-02
status: >
  CLOSES TWO ACTIONS against E11, by reading objects and committed text only.
  SECTION C, added 2026-10-02, records that the current draft and its
  limitation list are in no repository. One further item - the limitation-6
  sentence - is HELD, because the correction beside P3 that it cites was
  stopped before being written. See section C's last block.
  NOTHING WAS FITTED, RUN OR REFITTED, and no estimate here is new. (A) one of
  E11's two restatements was applied and one was not. (B) limitation 7 and
  E11's `$boot_ci` are DIFFERENT ESTIMANDS; the boot_ci result does not bear on
  limitation 7, and the E11 quantity that does bear on it agrees with it.
posture: EXPLORATORY. Nothing pre-registered. E11 itself carries no reading rule
         declared before it ran - see docs/2026-09-02_e11_prolif_adjusted.md.
relates-to:
  - docs/2026-09-02_e11_prolif_adjusted.md (P3 and P4, and the 2026-10-02 banner)
  - docs/2026-09-02_e10_machinery_and_priming.md (where P4 was applied)
  - docs/2026-09-02_paper_opening_human.md (the current draft)
  - results/prolif_adjusted_machinery.rds, results/curated_comparators.rds
---

# E11's two outstanding actions, closed

---

## A. "Correlates needs restating, in two places" - one done, one not

E11's verdict reads: *"The comparison is real and it survives every correction.
The word 'correlates' is where it needs restating, in two places."* The two
places are its own sections **P3** and **P4**.

### A.1 P4, the second restatement - APPLIED, in both places it needed to be

> **P4:** *"So E10 R1's 0.453 is a fact about mitochondrial localisation and
> expression, not about apoptosis."*

| where | state |
|---|---|
| `docs/2026-09-02_e10_machinery_and_priming.md`, lines 20-28 | **applied**, as a `> CORRECTION, 2026-09-02, from E11 section 3.1` block at the head of the note, before any of R1. P4's own claim that *"The E10 note now carries this correction at the top"* is accurate |
| `docs/2026-09-02_paper_opening_human.md`, claim ladder rung 4 | **applied**: *"...and OXPHOS therefore engages the death programme — **NOT SUPPORTED.**"* with the expression-matched and sub-compartment-matched nulls quoted |
| same draft, `## The two versions of "mitochondrial genes correlate more"` | **applied, and further than required**: it records that the plotted statistic was chosen after both were seen, that the plotted one has the larger z, and that **"the composition bound stays anchored on the rank split ... the honest word remains 'not separable'"** |

**Nothing is outstanding on P4.**

### A.2 P3, the first restatement - NOT APPLIED ANYWHERE

> **P3:** the finding is *not* "proliferation adjustment removes MYC's
> correlation with the machinery"; it is **"MYC never had one above background,
> and the correction is not what removes it."** Raw `z` for the 44 on MYC is
> **-0.95 (TCGA) / -0.81 (SCAN-B)** — at or below an expression-matched null
> **before any correction** — and **-0.86 / -0.74 (`__PROLIFSTRIP`)** and
> **-0.86 / -0.84 (`__BOTHSTRIP`)**, so it does not depend on the covariate
> being right. (`prolif_adjusted_machinery.rds$null_tests`, `set == "apoptotic
> machinery (44)"`, `adjustment == "raw"`, column `z_mean_abs`.)

**Searched: this argument appears in E11's own note and nowhere else.** No other
tracked document carries the phrasing or those z values.

**But it is a GAP, not a misstatement.** The draft does not say the wrong thing —
it does not make the claim at all:

- Rung 1 says *"correlate more strongly with OXPHOS than with MYC activity"* and
  **operationalises it in the same table cell** as the SD of the 44 per-gene rho
  (0.313 vs 0.191; 0.245 vs 0.157). The loose word is used and repaired on the
  spot.
- The draft's MYC-side argument is **rung 3**, the conditioning result —
  *"Conditioning MYC on OXPHOS takes it to zero (0.187 -> -0.043, 0.137 ->
  -0.058). MYC's ordering was inherited."*

**That is a different argument from P3's and it rests on more.** Rung 3 depends
on the conditioning being the right operation; P3 depends on nothing but the
raw numbers against a composition null, which is why E11 called it *"a cleaner
claim and a stronger one"*. **The draft currently uses the one that rests on
more and omits the one that rests on less.** Recording that is the whole of this
action; what to do about it is not decided here.

---

## B. Limitation 7 against `$boot_ci` - two different estimands

### B.0 One thing first: I could not find "limitation 7"

**No tracked file in any branch of any of the three repos contains a numbered
limitation 7 matching the description**, nor the phrase "rank argument".
Searched `*.md` across all seven exploratory branches and both single-branch
repos. The numbers it cites are traceable (below), so the reconciliation is
answerable regardless of where the list lives; **the pointer back into that list
is the one thing this note cannot supply.**

### B.1 What limitation 7's numbers are, and they come from two objects

| figure as cited | where it actually is | statistic |
|---|---|---|
| **-0.106 / -0.094** | `prolif_adjusted_machinery.rds$wide`, `OXPHOS_raw`, `mitocarta == FALSE` | **median** of the 24 cytosolic genes, TCGA -0.1062, SCAN-B -0.0944, **raw** |
| **"extreme of thirty comparators"** | `curated_comparators.rds$ladder_pct`, `pct_mean_rho` | **0.000 in both cohorts**, `n_reference_pathways = 30` — but computed on the **mean**, -0.097 / -0.068 |
| **"neither individually distinguishable from its own null"** | `curated_comparators.rds$half_tests`, `half == "cytosolic"`, `axis == "OXPHOS"`, `adjustment == "raw"` | **z = -1.07 (TCGA) / -0.60 (SCAN-B)** |

**All three statements are correct.** They are assembled from two objects and
two statistics — the level from E11 as a median, the rank from E14 as a mean.
Noted so the sentence can be sourced; it does not change the argument.

### B.2 The two estimands, side by side

| | **limitation 7's quantity** | **`$boot_ci`'s `split_gap` / `delta_rho_gap`** |
|---|---|---|
| **the contrast** | the cytosolic half's level against **thirty other pathways** | the mitochondrial-minus-cytosolic gap on the **OXPHOS axis** against **the same gap on the MYC axis** |
| **what varies** | which **pathway** you look at | which **axis** you look at |
| **the set** | the 24 cytosolic genes of the 44 | all 44, both halves, one set |
| **THE NULL** | an **expression-matched draw of genes** — i.e. composition | **zero** — i.e. the two axes split the compartments equally |
| **what the interval covers** | which **genes** were drawn (2,000 draws) | which **tumours** were drawn (1,000 resamples) |
| **adjustment** | raw | `PROLIF_DISJOINT`-adjusted (`cbind(rp)`, `E11:1019-1022`) |
| **result** | z **-1.07 / -0.60**; 0th percentile of 30 | **10 of 10 rows exclude their null**; `split_gap` TCGA **+0.266 [0.212, 0.331]**, SCAN-B **+0.352 [0.302, 0.395]** |

### B.3 THE ANSWER: different estimands, and `$boot_ci` does not bear on limitation 7

**`$boot_ci` contains no composition null.** Its null is zero — "OXPHOS orders
the compartments no differently from MYC". Clearing that null establishes that
**the OXPHOS axis separates the two compartments more than the MYC axis does,
and that this is not sampling noise over tumours.** It is silent on whether
either axis's separation is composition, because nothing in it is matched for
expression or sub-compartment.

**The E11 quantity that does ask limitation 7's question is `$split_null`, and
it agrees with limitation 7.** Same gene set, same axis, composition-matched
null, adjusted:

| cohort | null pool | observed | null mean +/- SD | z | pct of draws below |
|---|---|---|---|---|---|
| TCGA | all MitoCarta | 0.453 | 0.397 +/- 0.130 | **0.43** | 0.637 |
| TCGA | no OXPHOS/mitoribo (strict) | 0.453 | 0.289 +/- 0.142 | **1.16** | 0.877 |
| SCAN-B | all MitoCarta | 0.489 | 0.452 +/- 0.117 | **0.32** | 0.592 |
| SCAN-B | no OXPHOS/mitoribo (strict) | 0.489 | 0.362 +/- 0.135 | **0.94** | 0.827 |

**Maximum z = 1.16. Nothing clears.** This is the same conclusion limitation 7
states, reached inside E11 on the adjusted numbers.

> **So the two are not in tension and neither corrects the other.** A reader who
> sets `$boot_ci`'s ten-of-ten against limitation 7's "not individually
> distinguishable" is comparing a between-axes contrast under a zero null with a
> between-pathways rank under a composition null. **Limitation 7 stands as
> written. `$boot_ci` is a precision statement about a different question and
> must not be quoted in support of, or against, it.**

**The draft already behaves correctly on this**, which is why no text needs
changing: `docs/2026-09-02_paper_opening_human.md` anchors the composition bound
on the rank split, keeps "not separable", and records that the larger-z
statistic was chosen after both had been seen.

### B.4 What would change this answer

**A composition-matched null for `$boot_ci`'s contrast** — resample genes, not
tumours, matched on expression and sub-compartment, and ask whether the
OXPHOS-minus-MYC gap survives. **That object does not exist**, in E11 or
anywhere else in this repo. Until it does, the ten-of-ten result cannot be read
as evidence against composition.

---

## C. Two live documents sit outside version control

**Recorded 2026-10-02. Stated, not acted on.**

Two documents that the current work cites are **in no repository**:

| document | status |
|---|---|
| `2026-10-02_human_arm_handoff.md` | holds the limitation list, including **limitation 6** and **limitation 7** |
| `2026-10-02_section4_narrative_v2.md` | holds the current section-4 draft |

**Both exist only in the chat Project knowledge.** Verified 2026-10-02 across
`myc_human_exploratory`, `myc_human_validation` and `myc_mouse`:

- **not in any branch's tree** - all seven exploratory branches and both
  single-branch repos;
- **never added in any history** - `git log --all --diff-filter=A --name-only`
  returns nothing for either name, so they were not committed and later removed;
- **not on disk anywhere under `~/code`**, tracked or untracked.

### What follows from that, as fact

1. **The current draft and its limitation list are outside version control.**
   They have no commit history, no diff, and no provenance record. A number
   quoted from them cannot be traced to a version the way every number in
   `docs/` can.
2. **This is why the 2026-10-02 audit could not find "limitation 7".** That
   audit searched every tracked file in all three repos and reported the phrase
   absent. It was absent because the list is not in a repo - not because the
   limitation does not exist. Section B.0 above records the failed search; this
   section records the reason.
3. **It cuts against this arm's own convention.** `CLAUDE.md` places the aim,
   the plan, the dated notes and the tracked figures under `docs/`, on the
   stated ground that `outputs/` does not survive a fresh clone and that the
   durable record is the committed note. Every declaration, data note and result
   note in this arm - `E32` through `E38` - was committed before or alongside
   the work it governs, in the repo the work is about. **The draft and the
   limitation list are the exception.**

### What is NOT recorded here

- **No action is taken and none is proposed.** Where those two documents should
  live, and whether they should be committed, is not decided in this note.
- **The limitation-6 sentence is NOT filed here.** It was to read that `P3`, in
  its corrected magnitude form, answers limitation 6's objection to
  `M_b__PROLIFSTRIP`. **It is held** because the correction beside `P3` was
  stopped before being written: the step-1 read showed that `z_mean_signed`
  reverses the axis ordering - OXPHOS above MYC on `z_mean_abs` in 6 of 6 cells,
  MYC above OXPHOS on `z_mean_signed` in 6 of 6 - and that question is open.
  **Until P3's corrected form is settled, nothing may be filed that cites it.**
