# Snapshot from `myc_mouse`

The mouse mitoPPS artefact, copied once as a **read-only input** for a planned
cross-species comparison of **pathway RANK** (never value) in Phase 2.

**This is the only mouse material that enters this repo, and it enters as a data
artefact, never as a document.** Three mouse *documents* were found stale in
`docs/` during the week of 2026-09-10 and all three were deleted; one had drifted
in interpretation and not merely in a number. **A document lives in the repo it
is about.** See `docs/2026-09-10_mouse_corrections.md` section 0.

**Nothing is read across repos at runtime.** This is a copy, taken once, with the
source path, md5 and producing commit pinned below.

## Provenance

- **Source path:**
  `/Users/gs/G/data/MK_myc_2022/myc_mouse/results/mitopps_scores.rds`
- **md5:** `b8a125af4bc0d5909ad03f6e126e2890` (source and destination verified
  identical at copy)
- **Size:** 511,169 bytes
- **Copied:** 2026-09-10
- **`myc_mouse` HEAD at time of copy:** `129b0c0`
  (`129b0c080c4b19fd9e406967a77a5ad15611cda6`) on branch **`experimental-cohorts`**
- **Producing script:** `scripts/08_mitoPPS_analysis.R`
- **Last commit touching that script:** `0b9b28e`
  (`0b9b28eda687dc29e47ce9d27277a4133cd1d768`), **2026-07-24 12:29:01**,
  *"reconcile gene-symbol vintage (core set-based scripts)"*
- **Artefact mtime at source:** 2026-07-24 12:35:41 -- **later than the script
  commit**, so the object was built from the current version of its script and
  is not stale. Checked before copying.
- **The artefact's own `analysis_date`:** 2026-07-24

### The repo HEAD does NOT pin this file, and that matters

**`results/mitopps_scores.rds` is NOT tracked in `myc_mouse`.** That repo
gitignores `results/` wholesale (`.gitignore:10`), and unlike this repo it keeps
no tracked README beside it.

So **`129b0c0` records the state of the repo, not the state of this file.** What
determines the content is **the producing-script commit `0b9b28e`** together with
the mtime check above. If the artefact is ever refreshed, it is
`scripts/08_mitoPPS_analysis.R`'s history that says whether the content should
have changed -- the repo HEAD will have moved for unrelated reasons.

## Contents

One file. `mitopps_scores.rds` is a list of 21 elements; the five this repo needs:

| element | class | dim |
|---|---|---|
| `mitopps_scores` | `data.frame` | **24 x 148** -- 24 samples x (4 annotation + **144 pathway** columns) |
| `mitopps_group_means` | tibble | **576 x 4** -- `pathway`, `group`, `mean_mitopps`, `sd_mitopps` (144 x 4 groups) |
| `gene_to_pathway` | tibble | **3,881 x 2** |
| `pathway_levels` | tibble | **152 x 5** -- a lookup with tiering, wider than the scored panel |
| `n_pathways` | integer | **144** |

**Samples: 24. Group labels: `6W_neg`, `6W_pos`, `12W_neg`, `12W_pos`** -- the
four developmental groups, carried in `mitopps_scores$group`, with `timepoint`
(`6W` / `12W`) and `myc_status` (`neg` / `pos`) as separate columns. **The four
groups are separable and that is on record here.**

The remaining sixteen elements are carried unaltered and are named so a later
session knows what is available without opening the object: `raw_pathway_scores`,
`raw_group_means`, `raw_anova`, `mitopps_anova`, `raw_pairwise`,
`mitopps_pairwise`, `pathway_tier1_map`, `pca_raw`, `pca_mitopps`,
`n_mitocarta_genes_found`, `mtdna_genes_separated`, `mtdna_pathway_name`,
`apoptosis_pro_genes`, `apoptosis_anti_genes`, `analysis_date`, `description`.

## Rules

- **Do not edit in place.** If a refresh is wanted, re-copy from the source path
  and bump the md5, the HEAD and the producing-script commit above. Never patch
  the file here.
- **Read-only input.** No script in this repo writes to it, and nothing in this
  repo writes to `myc_mouse`.
- **Scores are cohort-relative and mitoPPS is composition-dependent.** Species is
  a cohort. **Never compare these values numerically with the human mitoPPS
  scores** -- the planned comparison is of **RANK within the panel**, which is
  why it is worth doing at all. CLAUDE.md traps 3 and 6.
- **The panel-composition check HAS NOW BEEN RUN, 2026-09-10, and it passed.**
  `scripts/E24_mitopps_rank_comparison.R` answered all three questions against
  the saved object, as required below: **144 against the human 142** with the
  **identical eight** pathways dropping under `min_genes = 3L` in both species,
  and the 13 `mt-` genes in exactly one synthetic pathway with **none inside
  `OXPHOS subunits` or the umbrella**. `OXPHOS subunits` is **89 in both**.
  The object's `description` claim was verified against `gene_to_pathway` rather
  than trusted. Full reconciliation table:
  `docs/2026-09-10_e24_mitopps_rank_comparison.md` section 1, and
  `results/mitopps_rank_comparison.rds$panel_reconciliation`. The paragraph
  below is kept because it fixes **how** the check must be run if it is ever
  re-run against a refreshed artefact.
  `docs/2026-09-10_mouse_corrections.md` section 7.2 fixes how to answer them:
  **against this saved object, never against the producing script**, because the
  `mt-` strip postdates the script's first version and an older copy of the
  script reads wrong with no sign that it is wrong. `gene_to_pathway` is where
  that check runs.
- **The object's own `description` field volunteers an answer** -- it states that
  `mt-*` genes were removed from the MitoCarta pathways and placed in a dedicated
  synthetic pathway. **That is a claim made by the script, recorded here as
  provenance, and it is not a substitute for the check.** Verify it against
  `gene_to_pathway` before relying on it.

## Reproduction check

One line. Run from the repo root; **`CURRENT` means the snapshot still matches
its source, `DRIFTED` means re-read this README before using the file.**

```sh
[ "$(md5 -q /Users/gs/G/data/MK_myc_2022/myc_mouse/results/mitopps_scores.rds)" \
  = "$(md5 -q data/from_myc_mouse/mitopps_scores.rds)" ] && echo CURRENT || echo DRIFTED
```

If it reports `DRIFTED`, do **not** silently re-copy: check
`git -C /Users/gs/G/data/MK_myc_2022/myc_mouse log --oneline -- scripts/08_mitoPPS_analysis.R`
first, because a changed artefact with an unchanged producing script means
something happened that this README cannot account for.

## Not copied, deliberately

Every other object under `myc_mouse` `results/`, and **every mouse document**.
Phase 2's rank comparison needs this one artefact. If a later phase needs
another, that is a decision with its own note and its own provenance entry here.

## Tracking

`mitopps_scores.rds` is **gitignored; this README is tracked** -- the convention
`data/from_validation/` already uses, and the reason it is regenerable is the
reproduction check above. **Note the size differs by three orders of magnitude
from the sibling's rationale**: `from_validation/` is ignored because it is
~560 MB, whereas this file is 511 KB and could be tracked without difficulty.
The convention was followed rather than the rationale re-derived; if this repo
would rather hold the artefact in git, that is a one-line `.gitignore` change and
a deliberate departure.
