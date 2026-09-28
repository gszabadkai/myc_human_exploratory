---
date: 2026-09-28
status: DATA NOTE. Records which files were read, from where, and what was
        found in them. No analysis, no verdict.
posture: EXPLORATORY. Nothing pre-registered. The declaration was committed
         ALONE at a71e19d BEFORE forkscale_models.rds was opened; this note is
         the second commit and records the opening.
branch: mb2-forkscale. NOT to be merged.
relates-to:
  - docs/2026-09-28_mb2_forkscale_declaration.md (the declaration - read first)
  - data/menegollo_biclusters/README.md (the pinned snapshot's provenance)
  - myc_human_validation/results/forkscale_models.rds (F3-pre, read after a71e19d)
---

# MB2 forkscale alignment - inputs, and what opening F3-pre revealed

**Nothing was downloaded.** Every input was already on disk. Nothing was
written to `myc_human_validation`.

---

## 1. The bicluster snapshot - this repo's own tracked copy

**Used: `data/menegollo_biclusters/` in THIS repo**, not the validation repo's.
CLAUDE.md says each directory under `data/` carries its own provenance README
with a pinned SHA and is consumed as-is. This repo has such a copy, and it is
**tracked** (unlike `data/from_validation/`, which is gitignored), so it
survives a fresh clone.

**The two copies are byte-identical**, verified this session:

| file | SHA-256 | this repo vs `myc_human_validation` |
|---|---|---|
| `TCGA_MB1_RNASeq_data_nonorm.RData` | `d5a94bf50fb924c1b2eeb6338ee2f5c4293586fd99b2fe1d9c50430abcbfe7d7` | **identical** |
| `TCGA_MB2_RNASeq_data_nonorm.RData` | `fba6918017956ebaec780fcd5de85559a54e329d8f2ae2027f9fed7fe463f330` | **identical** |
| `TCGA.all.biclusters.RNAseq.Rdata` | `892c32996588956fa7d9128e1f631d0fd0c43703f704cb30e578a46124e6500d` | **identical** |

**The upstream pin was re-verified rather than taken on trust.**
`git hash-object` on the two per-bicluster files returns
`8422be1ae0f0da85c3ff03a74acc624906a5127b` and
`17e652b904ffca12f8cc9c33daa60e561da2d1e5`, which are **exactly the upstream
blob SHAs recorded in `data/menegollo_biclusters/README.md`** for the pinned
commit `8fdbb3437ae5537055d5d5429411bdb3b333c04a`
(`gszabadkai/Menegollo_Bentham`, `main`, 2024-06-18). So the files are
byte-for-byte the pinned upstream objects, confirmed now and not only when the
snapshot was taken on 2026-08-29.

**Primary source is the per-bicluster files**, each holding
`TCGA.MB{n}.RNAseq.df`, 1,037 x 6, with `X` a 28-character aliquot barcode and
the columns `MB{n}.pc1` and `MB{n}.index`.
**`TCGA.all.biclusters.RNAseq.Rdata` is a cross-check only** - it is 849 rows,
not 1,037, because an upstream `inner_join` to two files that are not in the
repo drops 188 samples. Forkscale is recomputed here as
`MB{k}_forkscale = MB{k}.pc1 / MB{k}.index`, the plain form, per the
declaration.

**MB3 is not read.** Its sign convention differs between TCGA and METABRIC and
nothing here needs it.

## 2. `forkscale_models.rds` - opened AFTER the declaration was committed

| | |
|---|---|
| path | `myc_human_validation/results/forkscale_models.rds` (read-only) |
| SHA-256 | `34803b54e9f2688d6b8564e5d1ec93e0c81c2875e87cf5e95b77f5a9dde32f89` |
| built | 2026-08-30 22:37:14 |
| declaration committed | **`a71e19d`**, before the first read |

### 2.1 What it contains, and the disclosure that matters

The declaration (section 0) says the object "has not been opened, so the MB1
half of the question is as unknown to us as the MB2 half".

**That understates what is in it. `x$f3pre$correlations` carries the MB2
marginals as well as the MB1 ones**, computed in the validation repo on
2026-08-30, n = 885:

| fork | vs `M_a` | vs `M_b` | vs `PROLIF_DISJOINT` | vs OXPHOS subunits (GSVA) |
|---|---|---|---|---|
| `MB1_forkscale` | **+0.429** | -0.031 | +0.503 | +0.529 |
| **`MB2_forkscale`** | **+0.739** | +0.470 | **+0.874** | +0.480 |

**This is recorded here, in the commit that opens the object, rather than in
the result note, so that it cannot look like a discovery made afterwards.**

Three consequences, all of which constrain what the result note may claim:

1. **The marginal direction the declaration predicts was already computable on
   2026-08-30.** MB2 `+0.739` against MB1 `+0.429`. Nobody had looked - the
   object was unread until this session, after `a71e19d` - so the declaration
   is still a genuine pre-data commitment by the people writing it. **But the
   result note must not be written as though the direction was unknowable.**
2. **What this analysis actually adds** is therefore: the **partials**, which
   are not in the object and are the declaration's discriminating quantity;
   **this repo's own scoring** on both sides of the contrast; **bootstrap
   intervals**; and **M_c**.
3. **`MB2_forkscale` tracks proliferation harder than it tracks MYC** -
   `+0.874` against `+0.739`. See section 5.

### 2.2 The cross-check the declaration asked for

F3-pre's stored MB1 marginals are a **cross-check, not an input**. They were
computed on the validation repo's scoring at n = 885; this analysis computes
its own on this repo's scoring at whatever n the join gives. **Differences are
reported, not reconciled**, per declaration section 4.

### 2.3 One discrepancy against the declaration's own account

The declaration (section 4) says "M-b and M-c reported, since F3-pre declared
all three on the MB1 side". **F3-pre's stored correlations contain `M_a` and
`M_b` only. There is no `M_c` row for any fork.** M_c is still computed here,
as the declaration asks, but the stored object offers **no M_c cross-check** and
the result note will say so.

## 3. This repo's TCGA quantities

All read from existing objects; nothing is re-scored. Key is the **12-character
patient barcode**; the forkscale files are joined on `substr(X, 1, 12)`.

| quantity | source |
|---|---|
| `ox_gsva` | `tcga_brca_mito_scores.rds$gsva_arms["OXPHOS subunits", ]` |
| `ox_ppd` | `tcga_brca_mito_scores.rds$mitopps_arms["OXPHOS subunits", ]` |
| **`M_a`** | `results/new_set_scores.rds$tcga_gsva_new["FELSHER__MITOSTRIP", ]` |
| `M_b` | `results/new_set_scores.rds$tcga_M_b_variants["M_b__MITOSTRIP", ]` |
| **`M_c`** | `tcga_brca_myc_scores.rds$estimators$M_c_call` |
| `PROLIF_DISJOINT` | `tcga_brca_mito_scores.rds$gsva_cov["PROLIF_DISJOINT", ]` |
| ER status | `tcga_brca_covariates.rds$covariates$er_call` |
| Block C flag | `tcga_brca_covariates.rds$covariates$complete_block_c` |

**`M_c` is GISTIC copy-number amplification, not MYC mRNA.** `M_c_call` takes
GISTIC values -1/0/1/2 (26 / 341 / 449 / 224) and **carries 55 NA**, the only
estimator in the panel that does. It is a discrete call, not a continuous
score; Spearman handles the ties, but its interval will be wider and its
effective n smaller than M_a's, and that is a property of the estimator and not
a finding. MYC mRNA (`log2MYC`) is a fourth, separate variable and is not one
of the three the declaration names.

**Counts.** TCGA scores cover 1,095 patients; `complete_block_c` is **938**;
the forkscale files carry **1,037**. The declaration expects the
forkscale-to-Block-C overlap to be at least 880 and the brief says to stop
below 800. The exact number is computed and asserted by the script, not
assumed here.

## 4. STATE - not in this repo, read from the frozen repo

**The brief lists STATE among "this repo's" quantities. It is not here.** No
script in `myc_human_exploratory` builds or loads it; it exists only in prose.
It lives in the frozen validation repo:

```
myc_human_validation/results/state_definition.rds  ->  $tcga$state$STATE_gsva
```

with levels `MYC_low` (520), `MYC_high_OXPHOS_low` (188),
`MYC_high_OXPHOS_high_unbuffered` (140),
`MYC_high_OXPHOS_high_buffered` (192), NA (55).

**Read read-only, and used only for the declaration's P3**, which is
**context and explicitly not a test**. Reading a pinned object from that repo is
not reopening it; nothing is written there. STATE is pre-registered there,
which is why it is read rather than rebuilt: a STATE rebuilt on this repo's
scoring would be a different variable wearing the same name.

## 5. The proliferation companion - a decision made AFTER reading section 2.1

**Declared in this note, before the script was written, and labelled post-hoc
because it is.**

The declaration fixes no proliferation adjustment. CLAUDE.md trap 3 requires a
proliferation-adjusted estimate for any MYC correlation in this repo, and
section 2.1 makes the reason concrete rather than generic: **`MB2_forkscale`
correlates with `PROLIF_DISJOINT` at `+0.874`, higher than its `+0.739` with
`M_a`**, and `M_a` is itself 14.8% proliferation-entangled. So the single most
likely alternative explanation of a positive result is that MB2 forkscale is
close to a proliferation axis - which is, near enough, Menegollo's own
description of **MB1**.

**The declared verdict is computed on the declared quantities and nothing
else.** The following are computed alongside, reported beside it, and **cannot
change it**:

- `rho(M_a, MB2 | MB1, PROLIF_DISJOINT)` and its MB1 mirror;
- the same partials on `FELSHER__PROLIFSTRIP`, the proliferation-stripped MYC
  estimator, which asks the question without the covariate.

**This ordering is auditable**: the companion is named here, in the commit that
records opening F3-pre, and the script is committed unrun in the next commit.

## 6. Script numbering - E32, not E31

`E31` is free on this branch but **is taken on the unmerged
`bh3-mimetic-oxphos` branch** by `E31_bh3_mimetic_oxphos.R`. Both branches are
deliberately kept and neither is merged, so reusing `E31` would put two
different scripts under one number. **This analysis is
`scripts/E32_mb2_forkscale_alignment.R`.**
