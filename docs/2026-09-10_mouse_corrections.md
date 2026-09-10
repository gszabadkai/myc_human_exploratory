# Corrections from the mouse sample-identity error, and one input audit

**2026-09-10. EXPLORATORY. Nothing here is pre-registered.**

**No new analysis, no new script, no new numbers.** This note records six
corrections that follow from an error found and resolved in `myc_mouse`, and
answers an input audit for a cross-species mitoPPS **rank** comparison that has
not been scoped. **Every human number in this repo is unchanged.**

**N3 throughout.** Transcript associations. "Primed" appears nowhere.

---

## 0. The rule this note is written under

> **A document lives in the repo it is about.**
>
> The only mouse material that enters this repo is a **data artefact** under
> `data/from_myc_mouse/`, carrying a recorded **md5, source path and source
> commit**. Mouse documents are cited **by repo and filename** and are never
> copied in.

**This is not housekeeping, and 2026-09-10 supplied the demonstration.** Two
mouse documents were found sitting untracked in `docs/` here:

| file | lines | human md5 | `myc_mouse` md5 | |
|---|---|---|---|---|
| `2026-09-09_orthotopic_identity_correction.md` | 239 | `162ad760…` | `2c4bc4ed…` | **differ** |
| `2026-09-09_reprioritisation_narrative_v3.md` | 353 | `3e3c5fa2…` | `7478daab…` | **differ** |

Identical line counts, **different content**. The `myc_mouse` originals were
rewritten by commit **`97348a8` -- "The scoring set is 49, not 35, and I had
propagated the slip into three files."** The copies here predate it, so they
carry **a number that has since been corrected**. Had either been cited, this
repo would have imported a known-wrong figure and had no way to notice.

Neither was read beyond its opening lines and **no number in this note comes
from either.** Everything below was supplied directly or read from objects
already on disk here.

**Both were deleted from `docs/` on 2026-09-10.** Nothing was lost: both were
**untracked here with no git history in this repo**, and both are **committed in
`myc_mouse` at `97348a8`**, which is the corrected version. The mismatch was the
reason to delete rather than a reason to hesitate. **The authoritative copies are
cited by repo and filename in section 10 and are not reproduced here.**

---

## 1. V4 is ANSWERED, and it FAILED

**Synthesis 5.5 and the section 9 table are edited in place.** V4 was the
prediction that would move "the OXPHOS axis sets these levels" from association
toward cause: move OXPHOS, read the twelve.

**The mouse experiment intended to answer it was a double-construct arm** -- the
respiratory manipulation **plus** an anti-apoptotic transgene -- so its `BCL2L1`
readout was **construct transcript, not endogenous**. A later analysis found a
**construct-free** arm in the same cohort. There, with the respiratory programme
raised:

| endogenous | effect |
|---|---|
| `Bcl2l1` | **-0.113 log2 [-0.283, +0.174]** -- null |
| `Mcl1` | **-0.100 [-0.229, +0.047]** -- null beside it |

**Record V4 as answered and failed. The perturbation does not reproduce the
human `BCL2L1` association.**

**And read as exactly what it is: one mouse result at n = 7 v 7.** It is **not
evidence against the human finding**, which stands on its own cross-sectional
footing in two cohorts. What it removes is V4's availability as the outstanding
causal test -- the test was run and it did not deliver. Whether a differently
designed perturbation would is a separate question and nothing here answers it.

---

## 2. The orthotopic corroboration is VOID, everywhere it appears

**Searched, and it appears in exactly one tracked file:**
`docs/2026-09-07_orthotopic_analysis_plan.md`, in three places. No other
document or script in this repo mentions the orthotopic series.

| where | the void claim |
|---|---|
| §1(a) | PGC1a overexpression alone raised `Bcl2l1` 2.7x against EV, "with no Bcl-xL construct present"; called "the single cleanest argument" that the human `BCL2L1` coefficient is causal |
| §1(b) | the MYC ceiling read as "the in vivo demonstration that M and O cannot both be high, the collider argument shown by intervention" |
| §1(c) | the abstract sentence, and the gland/human comparison on the guardian pair |

**All void.** The arm labelled as the respiratory manipulation alone was not
alone. That file now carries a VOID banner; its body is preserved as the record
of what was planned and why, and is not to be quoted.

**The reverse direction is void too, and this is the sharper half.** In the same
cohort, an anti-apoptotic construct moved **no respiratory ruler**:

| ruler | effect |
|---|---|
| `ox_lvl` | +0.093 [-0.298, +0.581] |
| `ox_rel` | -0.023 |
| `ox_mt` | -0.246 |

**All three cover zero.** So "raising the guardian raises respiration" has no
support either, and no claim that mouse tumours reproduce the human
configuration survives.

---

## 3. The `BBC3` cross-species divergence is CLOSED, and the reason is the estimand

The mouse intervention shows `Bbc3` **falling** where the human coefficients
have it **rising**:

| | `BBC3` against the respiratory axis |
|---|---|
| human, TCGA, four rulers | **+0.293 to +0.366** |
| human, SCAN-B, four rulers | **+0.145 to +0.173** |
| mouse intervention | falls |

*(Human values from `docs/2026-09-04_e16_respiratory_rulers.md` §4.)*

**This is not an open conflict.** The mouse fall is a **survivor signature**.
Across three independent backgrounds it tracks how much transgene each tumour
**retained**, and goes flat where the transgene is lost:

| retained transgene | `Bbc3` |
|---|---|
| 302 CPM | **-0.731** |
| 129 CPM | **-0.345** |
| 0.17 CPM | **+0.044** |

**The two quantities are not the same estimand and cannot disagree.**

- The **mouse** quantity is a **between-arm contrast whose x-axis is retained
  transgene dose** -- an *outcome of selection*, with **no counterpart in a
  human tumour cohort**. Nothing in TCGA or SCAN-B has a transgene to retain.
- The **human** quantity is a **within-population partial correlation on
  unperturbed material**.

**Human values are unaffected and no human number changes.** What is recorded is
that the apparent divergence was a category error, not a discrepancy to
reconcile.

---

## 4. CDKN2A -- narrowly scoped, not excluded

An earlier instruction to exclude CDKN2A outright was **too broad**. `CDKN2A`
appears nowhere in this repo's notes or scripts, so there is no text to correct
-- this is the record of the correct position.

- **Not an anchor.** A p53/ARF-disabled cell line still shows the phenotype, so
  the route is **dispensable**. And Menegollo et al. 2024 already reports CDKN2A
  loss with the OXPHOS-high biclusters, so re-running it here **confirms prior
  work rather than advancing past it**.
- **But a valid instance.** It is one documented route past the gate, and **the
  only one with independent human evidence**. Cite the companion paper as
  support for that framing.
- **Do not run a CDKN2A-versus-OXPHOS analysis in this repo.**

---

## 5. The share ruler is already here -- `E16` tested it

**Recorded prominently so a future session does not propose porting a ruler this
repo has had since 2026-09-04.** `ox_rel` has been in **both cohorts since E16 /
G1**. Synthesis 2.3a is the full treatment.

| | |
|---|---|
| `BCL2L1 - BBC3` gap | **positive in 8 of 8** ruler-by-cohort cells |
| the gap on `ox_rel`, TCGA | **+0.049** |
| does the configuration track content or share? | **content, slightly** -- `BCL2L1` **+0.404** on `ox_lvl` against **+0.360** on `ox_rel` |

*(All four read back against `docs/2026-09-04_e16_respiratory_rulers.md` §4-§5.)*

`E20` then took the content-versus-share question further and returned
**ESTIMATOR-DEPENDENT**; `E21` showed the dissent there is a property of
composite **size**, not of regulons. Synthesis 3.11.

---

## 6. Two phenotypes, opposite valence, same variable -- scoping only

**Recorded; nothing is scoped.** The mouse work carries two phenotypes that move
in **opposite directions on the same variable**:

| axis | raised respiration is |
|---|---|
| transformation | **lethal** |
| growth and metastasis | **advantageous** |

**Anything in this repo touching outcome, recurrence or metastasis is on the
second axis, not the first.** The two must **not** be run together, and must not
be reported as one relationship. A single "OXPHOS-high is good/bad" sentence
spanning both is a sign the distinction has been lost.

Nothing in this repo currently sits on the second axis. Survival and treatment
are out of scope for phase 1 by the standing plan, which is why this is a
scoping note rather than a correction.

---

## 7. Input audit -- the mitoPPS rank comparison

A cross-species comparison of **pathway RANK** (never value) is planned: where
OXPHOS sits within the mitochondrial pathway panel in 6W mouse gland, 12W mouse
gland, and human tumours. **Nothing was built. These are the answers.**

### 7.1 A per-sample, full-panel human matrix EXISTS, in both cohorts

| cohort | object | element | dim |
|---|---|---|---|
| TCGA | `data/from_validation/tcga_brca_mito_scores.rds` | `$mitopps_universe` | **142 x 1,095** |
| SCAN-B | `results/scanb_scores.rds` | `$mitopps_universe` | **142 x 3,207** |

**Not only the OXPHOS column.** E16 read `ox_ppd` out of `$mitopps_arms`
(18 arms), but the full panel was retained alongside it.

**The 142 pathway rownames are identical, and in the same order, across both
cohorts.** TCGA additionally carries `$pathway_scores` (142 x 1,095), which is
the **mean-linear input** mitoPPS consumes rather than the score itself. SCAN-B
has **no** `$pathway_scores`; if the raw input is wanted there it must be
rebuilt.

### 7.2 Panel composition, per cohort

Identical in both, because both were scored from the same pinned pathway list.

| | |
|---|---|
| MitoCarta pathways before filtering | **150** |
| after filtering | **142** |
| the filter | **`min_genes = 3L`**, `.path_scores()`, `functions/mitopps.R:60-61` |
| the 8 dropped | `mtDNA modifications`, `OXA`, `Cytochrome C`, `Glycerol phosphate shuttle`, `Cholesterol-associated`, `Vitamin B1 metabolism`, `Vitamin B6 metabolism`, `Vitamin C metabolism` |

OXPHOS pathway sizes inside the panel:

| pathway | n |
|---|---|
| `OXPHOS` (umbrella, includes assembly factors) | **156** |
| `OXPHOS subunits` | **89** |
| `OXPHOS assembly factors` | 68 |
| `mtDNA-encoded OXPHOS subunits` | **13** |

**`MT-*` ARE stripped into a separate synthetic pathway, and this is the answer
that decides whether the panels are the same object across species.** Exactly
one pathway of the 142 contains `MT-` genes -- `mtDNA-encoded OXPHOS subunits`,
holding all 13 -- and **neither `OXPHOS subunits` nor the `OXPHOS` umbrella
contains a single one**. They are not left inside their complexes.

> **The mouse panel must be built the same way or the two are not comparable
> objects, and a rank comparison across them would be meaningless before it was
> even wrong.** This is the first thing to check on the mouse artefact.

**Two clarifications, so the check is not botched.**

**1. `89` and `86-87` are different quantities and their difference is not a
species difference.** The human **`OXPHOS subunits` = 89** is a **MitoCarta
pathway membership** -- a property of the catalogue. The **`ox_sub` = 86-87**
that appears in mouse scripts is the `ox_rel` **gene set after expression
filtering in a specific cohort** -- a property of that cohort's matrix. The
second is **always smaller** and is **cohort-dependent**; this repo's own
`E16` reports 87 mapping in one place against a nominal 88 in another for
exactly that reason. **Do not read 89 against 87 as mouse-versus-human.**

**2. Verify the `mt-` strip against the SAVED OBJECT, not the script.** The
`mt-*` strip was added to the mouse mitoPPS script **after its first version**,
so an old copy of the script shows no `mt-` handling while the saved object
already has it. Reading the script would give the wrong answer with no sign that
it was wrong. **The check, on `results/mitopps_scores.rds`:**

- does **`gene_to_pathway`** contain an **mtDNA pathway** at all, and
- does **any `mt-` gene sit inside `OXPHOS subunits`**?

Both must be answered from the object. This is the same failure mode as section
0 -- a stale copy that looks authoritative -- one level down, in code rather
than prose.

### 7.3 Input scale -- LINEAR, confirmed

Human mitoPPS was computed on **linear DESeq2-normalised counts**. Not log, not
VST. `E02` lines 72-74 assert `identical(scale, "linear_deseq2_normalised")` for
both cohorts before scoring; line 561 scores on it; the saved object records
`mitopps = "linear_deseq2_normalised"`. GSVA's log input never touches it --
the standing scale-discipline rule.

### 7.4 No mouse artefact is snapshotted here

`data/` holds `collectri_human`, `from_validation`, `genesets_celldeath_human`,
`genesets_from_library_human`, `genesets_metabolic_human`, `genesets_myc_human`,
`menegollo_biclusters`, `mitocarta_human` and `raw`. **There is no
`data/from_myc_mouse/` and nothing mouse-named anywhere under `data/`.**

**What would need copying in**, from `myc_mouse` `results/mitopps_scores.rds`:

| element | |
|---|---|
| `mitopps_scores` | samples x pathways, annotated -- the matrix the rank is read from |
| `mitopps_group_means` | |
| `gene_to_pathway` | needed to check 7.2's `MT-` question on the mouse side |
| `pathway_levels` | |
| `n_pathways` | the number to reconcile against 142 |

**Destination:** a new `data/from_myc_mouse/` with a provenance README recording
**md5, source path and source commit**, per section 0. **The mouse repo is not
added to this session's directories** and nothing is fetched unilaterally.

---

## 7a. Three traps carried in from the mouse branch notes

**Checked first: none of the three existed here under any wording.** Trap 5
notes that human confounds are human; trap 9 is about cross-species
separability. Neither is any of these. "Positive control" appears in `E16`,
`E19`, `E22` and `E23` only as a **reproduction** control, which is a different
object. **Added as traps 14, 15 and 16 in synthesis 8.4.**

Two of the three were addressed to this arm in the mouse notes -- *"tell that
session; it does not read this branch"* -- and had never arrived.

| trap | the rule | why it is live here |
|---|---|---|
| **14** | **A loading is a property of a SUBSET, not of a ruler.** Mouse fat-pad: `rho(ox_sub, Adipoq)` **+0.538 over all 30** samples, **+0.144 within the tested subset**, and a second ruler's loading **flipped sign** between them | this repo pre-specifies analyses **stratified by PAM50, TP53 status and purity**. An admissibility argument must be recomputed **inside each stratum**, never inherited from 1,095 tumours |
| **15** | **Carry a positive control through any composition or purity adjustment, and read it FIRST.** Mouse fat-pad: `Mki67` went **tau +0.313, p 0.045 -> -0.096, p 0.715** after adjustment, while the adjusters explained only **15%** of its variance | this repo adjusts for purity in TCGA and **cannot** in SCAN-B. The control is needed in both directions, and it is not the reproduction control already in use |
| **16** | **A rule keyed to ONE ruler can pass while the cohort says otherwise.** The pre-registered ruler is **itself a choice**. Mouse orthotopic: a rule passed on its named ruler while **two others in the same contrast excluded zero in the other direction** | **four rulers are in use here.** 2.3a and 3.11 both turn on which is read, and trap 12 is the estimator-side twin |

---

## 8. Flagged and NOT fixed -- the estimand audit

**A list only. Nothing here was touched.** These are places where a mouse
coefficient sits beside a human partial correlation. The audit predates this
correction and folding it in would get it half-done.

**THE SIX ARE NOT THE SAME KIND OF PROBLEM, and the audit should start with the
two that are.** A **juxtaposition** puts two numbers side by side; it can be
reworded and the surrounding claim survives. An **ARGUMENT** uses the comparison
as a *reason* -- and an argument built on comparing a mouse coefficient with a
human partial correlation **may not survive the audit at all.** Severity is
recorded so the audit session knows where to begin, not to pre-judge the outcome.

| # | severity | location | what sits together |
|---|---|---|---|
| **2** | **ARGUMENT -- start here** | synthesis ~1284, ~1294 | `bX_wt(Bcl2l1)` **+0.01**, a mouse slope at MYC-null, used as a **reason** the human **+0.39** cannot be explained away. Not a juxtaposition: the inference depends on the two being commensurable |
| **3** | **ARGUMENT -- start here** | synthesis 5.2 / ~789 | the same "`bX_wt` is flat" step, load-bearing in the reconciliation's own account of why the disagreement stands. **If the comparison is inadmissible, the explanation goes, not just the wording** |
| 1 | juxtaposition | synthesis 3.11 | gland **+0.351** and human **-0.31/-0.46** in one breath |
| 4 | juxtaposition | `docs/2026-09-04_e19_subtype.md` §3 | the per-gene `bX_wt` / `bX_myc` / `bMX` table beside human partial rhos |
| 6 | juxtaposition | `docs/2026-09-04_handoff.md` §5 item 0c, ~241 | `Bbc3`'s `bX_wt` listed as a thing to fetch from the mouse |
| 5 | **dormant** | `docs/2026-09-07_orthotopic_analysis_plan.md` §1(c) | **+0.351** against **-0.31 to -0.46**. Behind a VOID banner, so it cannot be quoted; fix it only if the file survives the decision in section 10 |

**Synthesis 8.4 trap 3 already forbids cross-species magnitude comparison, and
trap 13 now adds that sometimes the SIGN is not comparable either.** The audit is
whether they were obeyed, not whether the rules exist.

---

## 9. What this note changes, and what it does not

**Changes:** V4's status (answered, failed); the standing of the orthotopic
corroboration (void); the status of the `BBC3` divergence (closed, estimand
mismatch); the CDKN2A position (scoped, not excluded); the visibility of
`ox_rel` (already here); the phenotype-axis distinction (recorded).

**Does not change:** **any human number in this repo.** Not one coefficient,
interval, ruler or cohort value moves. The corrections are to what mouse results
were taken to *support*, never to what the human data say.

**And it does not touch N2.** `myc_human_validation` is frozen at `d3ac60e` and
was not read. The pre-registered `MYC x OXPHOS` null stands; nothing here reads
as overturning it, and V4's failure is about a *main-effect* perturbation
prediction, a different estimand entirely.

---

## 10. Where things live

| | |
|---|---|
| this note | `docs/2026-09-10_mouse_corrections.md` |
| edited in place | synthesis 5.5 and the section 9 table (V4); synthesis 8.4 (BBC3 estimand); synthesis 10 (CDKN2A); `docs/2026-09-08_handoff.md` §3 (ox_rel, phenotype axes) |
| VOID-bannered | `docs/2026-09-07_orthotopic_analysis_plan.md` |
| the mouse sources, cited never copied | `myc_mouse` `docs/2026-09-09_orthotopic_identity_correction.md`, `docs/2026-09-09_reprioritisation_narrative_v3.md`, `results/mitopps_scores.rds` |
| audit objects read | `data/from_validation/tcga_brca_mito_scores.rds`, `results/scanb_scores.rds`, `functions/mitopps.R`, `scripts/E02_score_cohorts.R` |
