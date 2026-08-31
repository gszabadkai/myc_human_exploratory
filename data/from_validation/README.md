# Snapshot from `myc_human_validation`

Derived objects copied from the validation study, which is **complete, closed
and not reopened**. This directory is the exploratory study's input layer.

**Nothing is read across repos at runtime.** These are copies, taken once, at a
pinned commit. That is the same contract `myc_human_validation/data/from_myc_mouse/`
uses for the mouse tables, and it exists for the same reason: a repo whose
inputs can change underneath it is not reproducible.

## Provenance

- Source repo: `/Users/gs/code/myc_human_validation`
  (`github.com/gszabadkai/myc_human_validation`)
- Pinned commit: **`d3ac60ec06c8e07df0626cbc063a52186ee167c5`** (`main`)
- Snapshot date: 2026-08-31

## Contents

| File | Scale / kind | Notes |
|---|---|---|
| `tcga_brca_vst.rds` | **LOG (VST)** | GSVA only. `$mat`, `$scale == "log_vst"` |
| `tcga_brca_linear.rds` | **LINEAR (DESeq2-normalised)** | mitoPPS and `log2(gene)` only |
| `tcga_brca_mito_scores.rds` | scores + **set definitions** | 18 arms x 4 instruments (`gsva_arms`, `mitopps_arms`, `content_arms`, `zmean_arms`), plus `arm_sets`, `covariate_sets`, `mito_paths`, `arm_universe_path`, `arm_summary` |
| `tcga_brca_myc_scores.rds` | scores | `M_a` (Felsher-61 GSVA), `M_b` (CollecTRI ULM), `M_c` (GISTIC amplification) |
| `tcga_brca_covariates.rds` | clinical + genomic | `PAM50`, `er_call`, `purity`, `leukocyte_fraction`, TP53/PIK3CA calls |
| `scanb_vst.rds` | **LOG (VST)** | GSVA only |
| `scanb_linear.rds` | **LINEAR** | mitoPPS and `log2(gene)` only |
| `scanb_pheno.rds` | phenotype + **`symbol_map`** | `PAM50`, `er_call`/`pgr_call`/`her2_call`, `age`, `NHG`; and the 2014-vintage symbol harmonisation map |
| `g1_overlap_audit.rds` | set definitions | `estimators_stripped$FELSHER` (61), `$COLLECTRI_MYC_ALL` (811) |

Cohorts: TCGA-BRCA 18,115 genes x 1,095 samples; SCAN-B (GSE202203) 18,153 x
3,207.

## Rules

- **Do not edit in place.** If the validation repo changes and a refresh is
  wanted, re-copy and bump the SHA above. Never patch a file here.
- **The two scales are opposite requirements and must not share an object.**
  GSVA wants log (VST, `kcdf = "Gaussian"`); mitoPPS wants linear
  DESeq2-normalised counts. Every script states which it is using.
- **The TCGA mitochondrial arms are reused, not recomputed.** Their values are
  already reported in the validation study; re-scoring risks a silently
  different answer for no gain. `E01` asserts three of them
  (`rho(M_a, OXPHOS subunits) = 0.388`, `Mitochondrial ribosome = 0.590`,
  `Fatty acid oxidation = -0.140`) as an end-to-end check that this snapshot is
  faithful.
- **SCAN-B must never be scored without `scanb_pheno.rds$symbol_map`.** It is
  annotated against a 2014 UCSC build, so 19 of the 89 `OXPHOS subunits` genes
  are pre-2018 ATP-synthase names; unharmonised the exposure covers 0.775
  instead of 0.989.

## Not copied, deliberately

`data/raw/` (2.8 GB of GDC/GEO/DepMap downloads), the neoadjuvant cohorts,
DepMap, GISTIC, MC3, the TCGA CDR, and every model object
(`block_c_models.rds`, `h4_outcome_models.rds`, `scanb_bim_replication.rds`, …).
Phase 1 needs none of them. If a later phase does, that is a decision with its
own note, and the raw data is re-downloadable from the URLs and checksums in the
validation repo's `data/*/README.md`.
