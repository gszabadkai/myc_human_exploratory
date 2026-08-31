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
- The mouse repo (`/Users/gs/G/data/MK_myc_2022/myc_mouse`) is **not** attached.
  Read it read-only via `git -C ... show <ref>:<path>` if ever needed.
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

**Phase 1: the correlation atlas.** Scripts `E00`–`E05`. Aim doc:
`docs/2026-08-31_aim.md`. Plan as approved 2026-08-31.

Out of scope for phase 1, named so they stay decisions rather than drift:
MCbiclust / forkscale (the Menegollo axis proper — the obvious phase 2),
survival, treatment, METABRIC, DepMap, causal or mediation modelling, and
anything that revisits the validation study's hypotheses.

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
scripts/    numbered R pipeline, E00-E05
docs/       the aim, the plan, dated notes
data/       snapshots, each with a provenance README
functions/  shared utilities
results/    intermediate .rds (gitignored, generated at runtime)
outputs/    figures and tables (gitignored, generated at runtime)
```

`results/` and `outputs/` are regenerable. `data/from_validation/` is
regenerable by re-copying from the validation repo at the pinned SHA.

## Git discipline

- `main` is the trunk. Feature branches off `main` as needed.
- Read-only git ops are always fine. Stop-and-check before anything destructive;
  never force-push a shared branch.
- **Never write to `myc_human_validation` or `myc_mouse` from this repo.**
