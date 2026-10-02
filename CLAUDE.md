# CLAUDE.md — MYC / OXPHOS / cell death in human breast cancer (exploratory)

Rules and context for Claude Code sessions in this repo. Read this first, then
`docs/2026-08-31_aim.md`.

## What this repo is

An **openly exploratory** study of how MYC activity and OXPHOS relate in human
breast tumours, and how the cell-death programme sits on that plane. Two
cohorts, TCGA-BRCA (n=1,095) and SCAN-B / GSE202203 (n=3,207). Starting point:
Menegollo, Bentham et al., *Cancer Res* 2024 (CAN-23-3172).

## THIS IS EXPLORATORY. THAT IS THE POINT, AND IT IS ALSO THE RISK.

The most important rule, and the mirror image of the sibling repo's.

- **Nothing here is pre-registered.** Every finding is hypothesis-generating and
  must be labelled as such, in notes, in figures and in conversation.
- **Multiple comparisons are the default state, not an exception.** The atlas is
  21 MYC estimators x 4 instruments x ~30 arms x 8 strata. Any single cell of it
  is uninteresting. Report structure, gradients and reproducibility across
  cohorts — never a p-value plucked from the grid.
- **"Significant" is not a result here.** Consistency across cohorts,
  instruments and estimators is. If something holds in TCGA and SCAN-B, on GSVA
  and mitoPPS, across low- and high-entanglement MYC signatures — that is worth
  something. One cell of the grid is not.
- **When something looks real, say what would falsify it** — and write that down
  *before* the next analysis, so the exploratory phase can hand a real
  hypothesis to a confirmatory one.

### The sibling repo, and how to read it

`myc_human_validation` (`/Users/gs/code/myc_human_validation`, frozen at
`d3ac60e`) is a **completed, pre-registered** study of one narrow hypothesis. It
found nothing supported. It is not reopened and it does not go in the paper.

- Its `CLAUDE.md` forbids post-hoc hypotheses. That rule is correct *there* and
  **actively wrong here**. Do not import it.
- Its *results* are trustworthy in a way nothing here is, because they were
  declared first. **Do not import a conclusion from it without checking whether
  it was pre-registered.** Its dated notes say which.
- Read it read-only:
  `git -C /Users/gs/code/myc_human_validation show d3ac60e:<path>`.
  Never write to it. It is not attached to sessions here.

## THIS IS A HUMAN REPO

- Human MitoCarta 3.0, human gene symbols, human-native gene sets only.
- The mouse repo (`/Users/gs/code/myc_mouse`) is **not** attached.
  Read it read-only via `git -C ... show <ref>:<path>` if ever needed. It moved
  there from `/Users/gs/G/data/MK_myc_2022/myc_mouse`, which is now renamed
  `myc_mouse_OLD_do_not_use`; **never point anything at an `_OLD_` copy.**
  Dated notes before 2026-09-28 still give the old path, as written.
- **No ortholog function is called anywhere in this repo**, in either direction.
  The check is for *calls*, not for the word — comments asserting the rule are
  the reason it must be narrowed to a `(`, or the tripwire always fires:

  ```
  grep -rnE "(mouse_to_human|human_to_mouse|ortholog[s]?)[[:space:]]*\(" scripts/
  ```

  must return nothing. It catches `mouse_to_human(`, `ortholog (` and
  `convert_orthologs(`; it ignores prose. The cell-death and MYC sets are
  human-native — see their READMEs for why that is established rather than
  assumed.

## Current phase

**Phase 2 is CONCLUDED as of 2026-09-11.** Start at
`docs/2026-09-11_phase2_conclusion.md` — it states the finding, what is not
settled, and four falsifiers, and it names the three scripts behind it (`E24`
superseded in part, `E25` the result, `E26` the composition check). **No phase
is currently open.**

Completed: **Phase 1, the correlation atlas** (`E00`–`E05`, aim doc
`docs/2026-08-31_aim.md`, plan approved 2026-08-31), then the exploratory
`E06`–`E23` line, then **Phase 2, the cross-species mitoPPS rank comparison**
(`E24`–`E26`), then **`E27`**, a companion presentation of Phase 2 on the native
MitoCarta level-1 partition with an internal null. `E27` re-reads `E25` and
changes nothing in the conclusion; see
`docs/2026-09-11_e27_category_readouts.md`. Then **`E28`**, a post-hoc gap
readout that **did not survive its own tests** and carries nothing on the human
side; it changes nothing in the conclusion either, and
`docs/2026-09-11_e28_gap_readout.md` section 8 says exactly what it may and may
not be used for. Then **`E29`**, a **map and explicitly not a readout**: it found
position matches between groups to be coincidental and established no opposing
programme, and it records that **mitoPPS is not a closed composition** (the
per-sample sum ranges, so the standard compositional objection does not apply).
See `docs/2026-09-12_e29_oxphos_neighbourhood.md`. Then **`E30`**, a second
comparator for `E28`'s gap which **failed the estimator ladder identically** -
two comparators, same numerator, 2 of 4 estimators each - demonstrating that the
failure belongs to the rank-gap statistic rather than the comparator. See
`docs/2026-09-12_e30_mitoribosome_and_nadh.md`.

Then **`E32`** (2026-09-28), outside any phase: which Menegollo bicluster this
arm's MYC-high, OXPHOS-high state corresponds to. The direction came from
Menegollo's published text, and the rule was committed before
`forkscale_models.rds` was opened. **Verdict MB2-ALIGNED**:
`rho(M_a, MB2 | MB1) = +0.688`, `rho(M_a, MB1 | MB2) = -0.195`, so conditioning
on MB2 *reverses* MB1's apparent MYC association. It survives proliferation
(+0.653 on `FELSHER__PROLIFSTRIP`, halved to +0.375 with `PROLIF_DISJOINT` as a
covariate). **It licenses one manuscript sentence and nothing about outcome.**
The sentence and the three qualifications that travel with it are in section 8
of `docs/2026-09-28_mb2_forkscale_result.md`. Note that the *marginals* already
existed in `forkscale_models.rds`; what `E32` adds is the partials.

**`E31` is not on `main`, and the gap is deliberate.** It is a BH3-mimetic
pharmacogenomics analysis (GDSC 8.5, PRISM 24Q2) on the unmerged branch
`bh3-mimetic-oxphos`, kept on origin so that nobody re-runs it without knowing
it has been done. Its declared primary is a precise null (an RB1 positive
control fires at p = 4e-13 in the same models) and is uninformative about the
manuscript, because `E18` found the BCL2L1-OXPHOS configuration absent from
cell lines. Read it with
`git show origin/bh3-mimetic-oxphos:docs/2026-09-28_bh3_mimetic_result.md`.

**`E33` and `E34` ARE on `main`**, merged 2026-10-01 at the author's request
once both were complete: declared alone before any data was opened, scripted,
run by the author, and verified object by object against digests recorded
*before* that run. Each kept its own branch on origin, undeleted, and each
arrived as a merge commit so its declaration / data note / unrun script /
result sequence survives as the audit trail.

- **`E33`**, merged from `e33-programme-specificity`. Is the BCL2-family death
  configuration specific to the MYC-coupled bicluster? **Verdict NEITHER
  SEPARATES.** Trigger arm `rho(arm, MB2 | MB1) = -0.140 [-0.200, -0.078]`
  against `rho(arm, MB1 | MB2) = +0.250`; guardian arm -0.350 against +0.517;
  unchanged on Block C, on the patient join and under the stricter reading. **It
  licenses no new sentence and WITHDRAWS one implication**: the manuscript may
  not state or imply that the standing pro-apoptotic configuration is confined
  to tumours whose respiration is MYC-coupled. A post-hoc observation that the
  configuration leans **MB1** is recorded with its own falsifier list and is not
  a finding. **Do not chain `E32` to `E33`**: the *state* is MB2-aligned, the
  *configuration* is not, and a property of the state is not a property of
  everything it correlates with. It is `docs/2026-09-30_e33_result.md`.
- **`E34`**, merged from `e34-quadrant-configuration`. The same configuration across
  the four MYC x OXPHOS quadrants, both cohorts, both instruments. **Verdict
  MIXED** — the rule demanded all 18 reading cells agree and 16 did; both
  dissents are mitoPPS within a subtype. What is new is the **Q2 cell**, MYC-low
  and OXPHOS-high, which `STATE` collapses into its level 1: it is large and
  carries the configuration with Q4 rather than Q1 (adjusted means +0.41 against
  +0.40 in TCGA, +0.24 against +0.37 in SCAN-B). **The MYC contrast is small but
  not zero**, and **proliferation adjustment enlarges it** by a uniform +0.12.
  Neither pre-written reading is licensed; the four-group description is.
  **The two OXPHOS contrasts are near-guaranteed by construction** — the
  configuration's genes and signs were read off OXPHOS in these same cohorts —
  so only the MYC contrasts carry information. It is
  `docs/2026-10-01_e34_result.md`.

**One `E33` input finding qualifies `E32`, which is on `main`.** The Menegollo
forkscale should be joined on the **aliquot** barcode, not the patient barcode:
six TCGA patients were profiled from a different vial, so `E32`'s patient-key
join paired those six across vials. On the aliquot join (n = 1,031; Block C
overlap 879) `E32`'s partials move by less than 0.01 and **its verdict and its
manuscript sentence stand** — but any future forkscale join should use the
aliquot.

**`E35`, `E36`, `E37` and `E38` are all on `main`**, each merged once complete,
run by the author and verified object by object against digests recorded
*before* that run — 25 of 25 for `E36`, 15 of 15 for `E37`, 30 of 30 for `E38`.
Branches kept undeleted on origin. **The next free script number is `E39`, and
`E31` is the only analysis still off `main`.**

They continue one line of questioning and they finish it: `E34` asked where the
configuration sits, `E35` whether that depends on genomic burden, `E36` whether
it exists below the tumour range at all, and `E38` diagnosed why `E36` could not
answer. **The normal-tissue line is closed; no phase is open.**

- **`E35`**, merged from `e35-burden-coupling`, **complete and verified**. Is the
  OXPHOS-to-configuration coupling graded by genomic burden? **Reading
  SATURATED** — the declared expected outcome. On `Aneuploidy.Score`, m1, GSVA:
  `Q2 - Q1` +1.073 / +0.651 / +1.144 and `Q4 - Q3` +0.929 / +0.853 / +0.745
  across tertiles, all six intervals excluding zero and the lowest and highest
  overlapping. **11 of 12 readable blocks agree.** Adding `PROLIF_DISJOINT`
  moves it by at most 0.054, so it is not proliferation under a burden label.
  **It licenses one clause and reinstates no stratifier**, and it cannot
  demonstrate the mouse's condition — a cohort of tumours has no sample below
  the threshold. `docs/2026-10-01_e35_result.md`.
- **`E36`**, merged from `e36-normal-comparison`. Is the coupling present in
  normal breast — the human analogue of the mouse wild-type gland, and the one
  analysis that could come back against the model. **Verdict UNINTERPRETABLE.**
  Its declared positive control — the **Fe-S cluster assembly cytosolic half**,
  named before any fit, with zero genes in common with either the configuration
  or the OXPHOS arm — reproduces **69% to 96%** of the configuration's
  normal-to-tumour change, in 6 of 12 cells. **The verdict does not rest on the
  control alone**: independently, all six cells are SPLIT BY ADJUSTMENT.
  **Q1 LICENSES NOTHING** — not ABSENT, not WEAKER, not PRESENT IN NORMAL. The
  manuscript may not say the coupling is absent in normal breast, nor present,
  nor that the human data do or do not corroborate the mouse gland result.
  **Remember why that matters: every *unadjusted* cell reads ABSENT IN NORMAL**,
  which is the predicted result and what a less careful analysis would have
  reported. Only the declared control and the adjustment-agreement rule stopped
  it. **The declared `EPITHELIAL` adjustment failed its own validity check** —
  its luminal and basal halves separate the tissues in opposite directions (AUC
  0.875 against 0.134) and cancel to 0.445 — **and the gene list was not changed
  and no substitute introduced.**
  **Q2 does not depend on the control and STANDS.** Over 113 matched pairs the
  tumour exceeds its own normal in 70% on `M_a`, but the two
  proliferation-stripped CollecTRI regulons **reverse the sign** to 39.8%. By
  quadrant: **Q1 48.7% and Q2 47.1% against Q3 100% and Q4 90.6%** — MYC-low
  tumours have no more MYC activity than their own matched normal, so the
  threshold reconciliation **fails for them**. The comparison is **paired** (all
  113 normals come from patients already in the tumour set) and the normals come
  from **6 of 40 collection sites**, 87.6% from three, so every across-person
  comparison is reported on **three nested tumour sets** and the within-patient
  contrast is the one immune to it. `docs/2026-10-01_e36_result.md`.
  **Its section 3 carries a dated correction from `E38`**: the numbers stand,
  two interpretive sentences do not, and the verdict is unaffected.
- **`E37`**, the snapshot rebuild `E36` needed, **GATE PASS**. It rebuilds the
  expression layer over tumours **and** normals jointly, because DESeq2 size
  factors over 1,208 samples are not size factors over 1,095 and including
  normals moves every tumour value. Its gate re-derives E34 on the rebuilt
  snapshot, tumours only: **quadrant agreement 97.1% / 98.1%, all sixteen
  contrasts within 0.10 (largest 0.047), E34's four pooled labels unchanged.**
  **So E34 is robust to a renormalisation that moves every input value.** It
  licenses that one sentence and **nothing biological**.
  `docs/2026-10-01_e37_result.md`.

- **`E38`**, **a DIAGNOSTIC and explicitly not a hypothesis test**.
  It decomposes quantities `E36` already computed, **adds no claim and cannot
  overturn `E36`'s verdict**; two of its four legs are **second attempts after
  seeing a first result** and are labelled post-hoc throughout. Its arbitrary
  sign splits were named **gene by gene** before it ran, and the script asserts
  them. Four findings:
  **(1) The signed near-zero is NOT cancellation.** Unsigning the configuration
  moves its normal-side correlation from -0.011 to **-0.037**, and arbitrarily
  signed comparators do **not** collapse (+0.834, +0.457) even at balanced
  splits. The leading explanation for `E36`'s result is refuted.
  **(2) The compartment gap is present in normals and LARGER there**, and is not
  specific to the apoptotic machinery — Fe-S shows +1.022 against +0.754. No
  reading taken; E14's compartment-matched null was not rebuilt.
  **(3) A luminal-only epithelial score FAILS** `E36`'s criterion applied
  unchanged, by **+0.048 on one leg of three**. Not used; criterion not relaxed.
  **(4) `MYC` mRNA is DOWN in 85% of matched pairs** (median -1.394 log2) while
  MYC *activity* is up in 70% of the same pairs — **trap 4 in its sharpest
  form** — and the 82 targets `M_b__PROLIFSTRIP` removes are **4.9x enriched in
  `HALLMARK_MYC_TARGETS_V1`** and 19.6x in `YU_MYC_TARGETS_UP`, including
  `ODC1`, `NCL`, `NME1`, `TFRC`, `CCND1`, `CDK4`, `E2F1`. **So
  `M_b__PROLIFSTRIP` is not "MYC activity with growth removed" and must not be
  described as such.** It licenses **three methodological statements and no
  biology**: `docs/2026-10-02_e38_result.md` section 7.

**`results/joint_tcga_*.rds` is gitignored and must be rebuilt after a fresh
clone**, in this order: `E37`, then `E36`, then `E38` — about five minutes in
all. **`E36` and `E38` stop with a clear error if the snapshot is absent.** The
tracked figures in `docs/figures/` and the numbers written into the notes are
the durable record. See `docs/2026-10-02_handoff.md`.

**One thing the normal-tissue line established that outlives it.** Every
limitation in it reduces to **bulk tissue composition**, and **no composition
adjustment built from bulk markers has yet survived its own validity check** —
the declared nine-gene score failed, the luminal-only score failed by +0.048.
`E38` also names an instrument this repo does not have: **a signed control**,
matched arbitrary genes with matched signs, which would separate "signed" from
"apoptotic".

Named as decisions rather than drift, and still not done: MCbiclust / forkscale
beyond `E32`'s single alignment check (the Menegollo axis proper), survival,
treatment, METABRIC, DepMap, causal or mediation modelling, and anything that
revisits the validation study's hypotheses.

Two items left open on purpose, recorded in
`docs/2026-09-11_handoff.md` (spent, kept for this): whether standalone
mouse-only figure reproductions are wanted — they would be a `myc_mouse` task
and cannot be built from here — and whether `data/from_myc_mouse/` should be
tracked rather than gitignored, since it is 511 KB and the convention it copied
exists for a 560 MB directory.

## Section 2 — the traps

Each is measured, not anticipated. Numbers are from the snapshot.

1. **Correlation is not the interaction.** The validation study found the
   `MYC x OXPHOS` *interaction* on apoptotic priming is null. A strong
   MYC–OXPHOS *correlation* is entirely compatible with that. They are different
   questions. Never present one as confirming or contradicting the other.
2. **Purity and immune infiltrate.** In TCGA `rho(OXPHOS subunits, purity) =
   0.214`, `rho(., leukocyte fraction) = -0.158`. Breast is the worst TCGA
   tissue for this: adipose is OXPHOS/FAO-high and infiltrate carries its own
   BCL2-family profile. **SCAN-B has no purity estimate.** Report raw and
   adjusted in TCGA, raw only in SCAN-B, and say so on the figure.
3. **Every MYC activity signature is entangled with proliferation, by wildly
   different amounts** — 1.5% (`MYC_UP.V1_UP`) to 47.6% (`YU_MYC_TARGETS_UP`),
   with the validation study's `M_a` at 14.8% and `HALLMARK_MYC_TARGETS_V1` at
   23.5%. Never report a MYC–OXPHOS correlation from one signature. Report the
   panel ordered by entanglement, and the proliferation-adjusted estimate.
4. **MYC mRNA is not MYC activity, and the difference is total.** In TCGA
   `rho(log2(MYC), OXPHOS subunits) = -0.032` against `+0.388` for the activity
   signature. Someone plotting MYC expression against OXPHOS sees nothing. Both
   are reported and the gap is itself a result.
5. **The four instruments disagree, by a lot.** GSVA-vs-mitoPPS agreement across
   arms runs 0.24 (lipid metabolism) to 0.94 (mtDNA). Instrument choice is not
   cosmetic: report all four, or justify the one.
6. **mitoPPS is blind to level by design.** It answers "is OXPHOS *prioritised*
   relative to other mitochondrial programmes", not "is OXPHOS high". Never
   compare mitoPPS values numerically across cohorts — only patterns.
7. **SCAN-B's symbols are a 2014 UCSC build.** 19 of the 89 `OXPHOS subunits`
   genes are pre-2018 ATP-synthase names; unharmonised the exposure covers 0.775
   instead of 0.989, and 70 of 89 genes is a Complex V with no F1 head and no
   c-ring. **Never score SCAN-B without `scanb_pheno.rds$symbol_map`.**
8. **mtDNA-encoded genes are held separately** and never pooled with
   nuclear-encoded subunits — expression-scale skew, and they behave differently:
   `rho(M_a, mtDNA-encoded OXPHOS) = +0.068` against `+0.388` for the nuclear
   subunits, and **negative on mitoPPS**.
9. **CICD is thin and must not be over-read.** 13 pro-death and 4 pro-survival
   human genes. It is the axis of most interest and the weakest measured. Score
   only what clears n >= 5; show the genes individually; never present a 4-gene
   GSVA score as a programme.
10. **Two Tang sets are near-transcriptome-wide**, and their sizes are easy to
    overstate. The CSVs carry one row per gene-per-evidence, so `Ferroptosis` is
    935 rows but **600 genes** and `Autophagy_dependent` 1,195 rows but **876** —
    4.8% and 3.3% of the matrix. Always count distinct genes. A correlation with
    either is still close to one with general expression: read them against a
    size-matched comparator before believing anything.
11. **The cell-death and MYC sets are human-native and must never be remapped.**
    Both carry first-class human columns and the upstream README says "Do NOT
    remap". See `data/genesets_celldeath_human/README.md`.

## Scale discipline — the most likely silent error

- **GSVA / ssGSEA** want log-scale input: VST, `kcdf = "Gaussian"`.
- **mitoPPS** wants linear DESeq2-normalised counts.

Opposite requirements. **They must not share an input object.** State the scale
in a comment at the top of every scoring block.

- **GSVA is cohort-relative.** Score all samples of a cohort in one run, and all
  sets of interest in the *same call* — the `.PIN_A`/`.PIN_B` half-matrix pins
  hold the gene universe, and two calls with different set collections are not
  comparable without them. **Never pool scores across cohorts**; compare
  correlations and patterns, not values.
- **mitoPPS baseline is composition-dependent.** It reports the *shape* of the
  mitochondrial programme, not its level.

### Two normalisations of TCGA-BRCA now exist. Say which one you read.

From `E37` there are **two** normalised TCGA expression layers in this repo, and
**a value from one is not comparable with a value from the other** — the DESeq2
size factors and the low-count gene filter are both computed over a different
sample set, so every value differs.

| snapshot | what | read by |
|---|---|---|
| `data/from_validation/` | 1,095 tumours, 18,115 genes | `E11`, `E33`, `E34`, `E35` and everything before them |
| `results/joint_tcga_*.rds` | 1,208 samples (1,095 tumour + 113 normal), 18,142 genes | **`E36` and `E38`, and nothing else unless separately declared** |

- **Every script states which snapshot it reads**, in a comment at the top of
  its scoring block, exactly as the scale rule above already requires.
- **No figure, table or sentence combines a quantity from one with a quantity
  from the other.** `E37`'s gate is the sole exception: it is a deliberate
  like-for-like comparison of the *same* quantities and is labelled as such.
- **Neither replaces the other.** Replacing `data/from_validation/` would
  invalidate `E11`, `E33`, `E34` and `E35` at a stroke.
- **One genomic annotation crosses, and it is not an expression value.**
  `E36` and `E37` carry `BUFFER_gistic`, and `E38` reads `MYC_amp`, from the
  frozen covariate table. These are per-patient GISTIC calls, not derived from
  any expression normalisation. **No expression value crosses between
  snapshots.**

## Gene sets — consume the snapshots, do not rebuild

Each directory under `data/` carries its own provenance README with a pinned
commit SHA. Consume as-is; do not rebuild here, do not edit in place. If a
source changes, re-snapshot and bump the SHA rather than patching.

| Input | Location |
|---|---|
| Validation-study matrices, scores, covariates | `data/from_validation/` (gitignored, ~563 MB) |
| Human MitoCarta 3.0 | `data/mitocarta_human/` |
| Cell death: pro/anti x apoptosis/CICD, + 15 Tang modalities | `data/genesets_celldeath_human/` |
| MYC signature compendium (16 sets) | `data/genesets_myc_human/` |
| CollecTRI regulons | `data/collectri_human/` |
| Felsher signature (library v1.0) | `data/genesets_from_library_human/` |
| Menegollo bicluster forkscale | `data/menegollo_biclusters/` |
| Curated metabolic genes | `data/genesets_metabolic_human/` |

**The library's `outputs/gmt/human/` tree is mouse-derived and must not be
loaded** — it is mouse-native sets pushed through `mouse_to_human()`, gitignored
upstream and unpinned by the tag. The rejection is of that tree only; tracked
raw inputs carrying their own native human data are a different object. **Check
the sheet, not the repository.**

### Standing conventions

- mtDNA-encoded protein-coding genes (13, `MT-` prefix) sit in their own
  synthetic pathway and are never pooled with nuclear-encoded OXPHOS subunits.
- MitoCarta's `OXPHOS` umbrella includes assembly factors; `OXPHOS subunits` is
  the narrower set. They are different — pick deliberately.

## R coding rules

Full file: `docs/R_CODING_INSTRUCTIONS.md`. These three cause the most damage.

1. **Never `print(n = X)` after `head()`.** `head()` may coerce a tibble to a
   data.frame, so `n` is read as `na.print`. Use `head(X) %>% print()`.
2. **Always `dplyr::count()`**, never bare `count()` — namespace conflicts.
3. **ASCII-only strings in scripts.**

No `renv`; packages are installed system-wide.

## Workflow — "Option A" (do not deviate)

- Claude Code **writes and edits** the numbered pipeline scripts. It does **not
  run them.** The author sources them in Positron interactively.
- Infrastructure (git, snapshots, provenance READMEs, editing this file,
  planning and result notes) Claude Code may execute directly.
- Every numbered script ends with an `if (FALSE) { ... }` sandbox block —
  skipped by `source()`, run line-by-line in Positron for inspection.
- Commit per verified phase. Git is the safety net.
- When in doubt, ask.

## Project structure

```
scripts/       numbered R pipeline, E00-E38. `E31` lives on an unmerged
               branch, so it is the one gap in this tree
docs/          the aim, the plan, dated notes
docs/figures/  tracked copies of the figures a note relies on
data/          snapshots, each with a provenance README
functions/     shared utilities
results/       intermediate .rds (gitignored, generated at runtime)
outputs/       figures and tables (gitignored, generated at runtime)
```

From `E31` on, a script writes its key figure twice: to `outputs/figures/` as
usual, and to `docs/figures/`, which is tracked. `outputs/` does not survive a
fresh clone, so the tracked copy and the numbers written into the note are the
durable record.

`results/` and `outputs/` are regenerable. `data/from_validation/` is
regenerable by re-copying from the validation repo at the pinned SHA.

## Git discipline

- `main` is the trunk. Feature branches off `main` as needed.
- Read-only git ops are always fine. Stop-and-check before anything destructive;
  never force-push a shared branch.
- **Never write to `myc_human_validation` or `myc_mouse` from this repo.**
