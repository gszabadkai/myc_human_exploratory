# Menegollo / Bentham bicluster snapshot (MCbiclust forks and forkscale)

The companion paper's MCbiclust output: per-sample bicluster PC1 values and sort
indices for the three mitochondrial biclusters MB1, MB2, MB3, in **TCGA-BRCA** and
**METABRIC**. This is the input to the F3 incremental-value stop gate (plan section 10)
and to Block D / script 15 (plan section 10, ED2).

**These are human data.** The companion paper is a human breast cancer study; nothing
here is mouse-derived and no ortholog projection is involved. This snapshot does not
touch the `myc_mouse` repo.

**Do not edit in place.** If the upstream repo changes, re-snapshot and bump the commit
SHA below rather than patching these files.

## Provenance

- Upstream: `github.com/gszabadkai/Menegollo_Bentham` (private; readable via `gh` as the
  repo owner). Code and data for Menegollo, Bentham et al., *Cancer Research* 2024
  (CAN-23-3172), the companion paper.
- Pinned commit: `8fdbb3437ae5537055d5d5429411bdb3b333c04a` (`main`, 2024-06-18,
  "Update README.md")
- Snapshot date: 2026-08-29
- Method: `gh api repos/.../git/blobs/<sha>` -> base64 decode. Every file was then
  re-hashed locally with `git hash-object` and matched to its upstream blob SHA, so
  these are byte-for-byte identical to the pinned commit.

| File | Upstream path (under `Input and output data/`) | Bytes | Blob SHA (verified) |
|---|---|---|---|
| `TCGA_MB1_RNASeq_data_nonorm.RData` | `TCGA analysis/` | 242,449 | `8422be1ae0f0da85c3ff03a74acc624906a5127b` |
| `TCGA_MB2_RNASeq_data_nonorm.RData` | `TCGA analysis/` | 241,801 | `17e652b904ffca12f8cc9c33daa60e561da2d1e5` |
| `TCGA_MB3_RNASeq_data_nonorm.RData` | `TCGA analysis/` | 241,236 | `fe82edaca89aa146ddbe10c9e88ad92cb8755bff` |
| `TCGA.all.biclusters.RNAseq.Rdata` | `TCGA analysis/` | 191,563 | `302f372120f3a522bf89b4b4fe9bb9b0f237dade` |
| `METABRIC_sort_data.RData` | `METABRIC analysis/` | 1,060,560 | `f33b29916f3cf0a1949381203b45b2c971f1c0a1` |
| `METABRIC_PC1_GSEA.RData` | `METABRIC analysis/` | 317,698 | `c2d6c7cc0061a96b8e978a9f4d60fd063b400f0f` |

To re-read any upstream file later without cloning:

```
gh api "repos/gszabadkai/Menegollo_Bentham/git/blobs/<blob-sha>" --jq '.content' \
  | tr -d '\n' | base64 -d > <dest>
```

## Contents

### TCGA, per bicluster - the primary source. Use these.

`TCGA_MB{1,2,3}_RNASeq_data_nonorm.RData`, each holding two data frames:

| Object | Dim | Columns |
|---|---|---|
| `TCGA.MB{n}.RNAseq.df` | 1037 x 6 | `X`, `MB{n}.pc1`, `MB{n}.index`, `PAM50`, `LumA`, `LumB` |
| `TCGA.MB{n}.CV.RNAseq.df` | 20531 x 2 | `Genes`, `CV` - the correlation vector |

- `X` is a **28-character TCGA aliquot barcode**, e.g. `TCGA-A2-A0YK-01A-22R-A109-07`.
  All 1,037 are sample type `01` (primary tumour); normals (`11`) and metastases (`06`)
  were removed upstream. Patient key is `substr(X, 1, 12)`; 1,037 rows, 1,037 unique
  patients, no duplicates.
- `MB{n}.index` is the MCbiclust `SampleSort` rank, a permutation of `1:1037`.
- `MB{n}.pc1` is `PC1VecFun` on the mean-centred bicluster gene block.
- The `CV` frames are on the old TCGA RNASeqV2 20,531-gene axis and their `Genes`
  column contains `"?"` for unannotated rows. **Not needed for F3 or Block D** - only
  `pc1` and `index` are.

### TCGA, assembled - cross-check only, NOT the primary source

`TCGA.all.biclusters.RNAseq.Rdata` -> `TCGA.all.biclusters.RNAseq.df2`, **849 x 213**.
It carries `MB{1,2,3}.forkscale`, the `.log` variants, the fork calls, plus histology,
purity, proliferation, and a large CNV/mutation panel.

**It is 849 rows, not 1,037.** The upstream `Fig_7ABC_...R` `inner_join`s it to
`consensus.csv` and `TCGA_molecular_profiles.csv`, neither of which is in the upstream
repo, so those 188 samples cannot be recovered from this snapshot. Recompute forkscale
from the per-bicluster files to keep all 1,037 and use this frame only to verify.

Verified 2026-08-29: recomputing `MB1.pc1 / MB1.index` from the per-bicluster file
reproduces this frame's `MB1.forkscale` with **max absolute difference 0** across all
849 shared rows. The definition below is confirmed, not assumed.

### METABRIC - incomplete here, see "The METABRIC gap"

| File | Objects |
|---|---|
| `METABRIC_sort_data.RData` | `mito.sort` (3), `ICT1.sort` (1), `random.sort` (2), and the matching `*.cv` / `*.top.seed` lists |
| `METABRIC_PC1_GSEA.RData` | `mito.pc1` (3), `ICT1.pc1` (1), `random.pc1` (2), each element length **1981**; `all.gsea` (6) |

## The forkscale definition

Transcribed from the upstream scripts. `MB1 = ICT1`, `MB2 = Mito2`, `MB3 = Mito1`
(upstream nomenclature note in the repo README).

```
forkscale      = pc1 / index
forkscale.log  = pc1 / log(index)      # the paper's preferred variant
```

TCGA, from `R scripts/TCGA further analysis - figures/Fig_7ABC_S7AB_...R`:

```r
MB1.forkscale = MB1.pc1 / MB1.index
MB2.forkscale = MB2.pc1 / MB2.index
MB3.forkscale = MB3.pc1 / MB3.index
```

METABRIC, from `R scripts/METABRIC MCbiclust analysis/compile_all_METABRIC_data_...R`:

```r
sort.df <- data.frame(sample = Complete_METABRIC_Clinical_Features_Data[, 1],
                      MB3.index = order(mito.sort[[1]]), MB3.pc1 = mito.pc1[[1]][order(mito.sort[[1]])],
                      MB2.index = order(mito.sort[[2]]), MB2.pc1 = mito.pc1[[2]][order(mito.sort[[2]])],
                      MB1.index = order(ICT1.sort[[1]]), MB1.pc1 = ICT1.pc1[[1]][order(ICT1.sort[[1]])])
MB1.forkscale = MB1.pc1 / MB1.index
MB3.pc1.rev   = MB3.pc1 * -1            # NOTE the sign flip
MB3.forkscale = MB3.pc1.rev / MB3.index
MB2.forkscale = MB2.pc1 / MB2.index
```

## Read-time traps (all verified 2026-08-29 against the snapshot)

- **`forkscale.log` is `Inf` for the sample at `index == 1`.** `log(1) == 0`. The
  published assembled frame carries exactly one `Inf` in `MB1.forkscale.log` and one in
  `MB3.forkscale.log`; `MB2` escapes only because its index-1 sample was dropped by the
  inner join. The full 1,037-row MB1 frame also contains an `index == 1` row. Any model
  on a `.log` variant must decide explicitly what to do with that sample - dropping it
  silently through `complete.cases()` will not happen, because `Inf` is not `NA`.
- **MB3's sign convention differs between cohorts.** METABRIC reverses `MB3.pc1` before
  forming forkscale (and negates `MB3.forkscale.log`); the TCGA script does neither.
  MB3 forkscale is therefore **not sign-comparable across the two cohorts** as stored.
  Fix the direction explicitly in any script that uses MB3, and say which convention.
  This matters: the plan uses MB3 as the ER-neutral control axis in Block D.
- **`forkscale` is heavily skewed by construction.** `pc1 / index` with `index` near 1
  gives huge values; observed MB1 range is `[-42.6, 8.8]` over 849 samples. This is why
  the paper prefers the log variant. Neither is symmetric; do not assume normality.
- **The assembled frame is 849, the per-bicluster frames are 1,037.** Say which you used.
- **Barcode length.** `X` is 28 characters (aliquot). Our covariate table is keyed on the
  12-character patient barcode. Join on `substr(X, 1, 12)`.
- **Cohort overlap with this repo.** All **1,037** Menegollo TCGA patients are present in
  our 1,095-patient covariate table (intersection 1,037; verified 2026-08-29). Our Block
  C/B analysis set is 938, so the overlap with it is **at least 880** and will be
  reported exactly by script 15.

## The METABRIC gap - one file is still missing

METABRIC forkscale is *almost* reconstructible from this snapshot:
`METABRIC_sort_data.RData` supplies the sort orders and `METABRIC_PC1_GSEA.RData`
supplies the PC1 vectors (n = 1,981). What is missing is the **sample identifier
vector**, which the upstream compile script takes from
`Complete_METABRIC_Clinical_Features_Data[, 1]` inside `METABRIC_DATA.RData`. The sort
and PC1 objects are positional indices into that frame, so without it the METABRIC
forkscale values cannot be attached to patients.

`METABRIC_DATA.RData` is not in the upstream repo; the upstream README links it on Google
Drive. Either that file, or the compile script's own output
(`METABRIC_starting_data_final_corr_groups.Rdata`, or `all.clinical.df` /
`METABRIC_data_w_forkscales`), closes the gap. See `docs/2026-08-29_G3_result_forkscale_availability.md`.

TCGA needs nothing further.
