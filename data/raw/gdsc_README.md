# GDSC release 8.5 - drug response, for the BH3-mimetic analysis

Provenance for the files consumed by `scripts/E31_bh3_mimetic_oxphos.R`.
Downloaded **2026-09-28**. Files live in `data/raw/gdsc/`, which is
**gitignored and not on origin**; this README is the record and the
re-download commands below are the reproduction route.

Follows the convention of `myc_human_validation/data/neoadjuvant/README.md`:
exact download URL, byte size and SHA-256 for every input.

---

## THE SOURCE MOVED, AND THE OLD URL IS GONE

**`cancerrxgene.org/downloads/bulk_download` now returns HTTP 410 Gone.** The
GDSC datasets have been folded into the Sanger DepMap, and the bulk download
page is a static notice pointing at Cell Model Passports. Anything written
before 2026 that names a `cancerrxgene.org` bulk path, or a
`cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/...` path, is stale.

**The data itself has NOT moved on.** The current release is still
**GDSC release 8.5, fitted 27 Oct 2023** - the same files, re-hosted. There is
no newer fitted dose-response release to be missing.

Current host: `https://cmp.cog.sanger.ac.uk/download/` (drug response) and
`https://cog.sanger.ac.uk/cmp/download/` (model and genomic annotation). Both
are plain HTTPS, **scriptable with `curl`, no browser and no challenge page** -
unlike the DepMap portal, which is Cloudflare-gated (see
`data/raw/depmap_README.md`).

The file list is discoverable from the Passports API rather than by scraping
the JavaScript page:

```sh
curl -s "https://api.cellmodelpassports.sanger.ac.uk/download_files?include=groups.files"
```

## Download

About 63 MB. Re-runs from a clean checkout.

```sh
mkdir -p data/raw/gdsc && cd data/raw/gdsc
B=https://cmp.cog.sanger.ac.uk/download
C=https://cog.sanger.ac.uk/cmp/download

# --- drug response, GDSC release 8.5 ---
curl -LO $B/GDSC2_fitted_dose_response_27Oct23.xlsx     # PRIMARY screen
curl -LO $B/GDSC1_fitted_dose_response_27Oct23.xlsx     # legacy, navitoclax replication
curl -LO $B/screened_compounds_rel_8.5.csv              # DRUG_ID -> name, target

# --- model and genomic annotation, for the declared RB1 confounder ---
curl -LO $C/model_list_20260921.csv
curl -LO $C/cnv_summary_20260804.zip      && unzip -o cnv_summary_20260804.zip
curl -LO $C/mutations_summary_20260724.zip && unzip -o mutations_summary_20260724.zip
```

## Checksums

`shasum -a 256 data/raw/gdsc/*`

| File | Bytes | SHA-256 |
|---|---|---|
| `GDSC2_fitted_dose_response_27Oct23.xlsx` | 22,112,537 | `0e99120a14a003dcbd3c27c020f1cfd35e38ee091067a2406f66c0d4dbc1c899` |
| `GDSC1_fitted_dose_response_27Oct23.xlsx` | 30,348,857 | `6dfad242263d93b6814cc51146b15346d3c69abe1aba834a0977ba9e3fd3b6a5` |
| `screened_compounds_rel_8.5.csv` | 46,414 | `d3cbb8b595980b36bf92d54af1da1dd995d2e2ce624c7cc70aacf2a754dce782` |
| `model_list_20260921.csv` | 957,555 | `419f65d1d6534454194f1a744550fab632cd1d7411f6b13411a885735aa37006` |
| `cnv_summary_20260804.zip` | 9,098,915 | `17c816f577e25503880b3d443870cd7a0c7300abee9214caffec3ebff7512277` |
| `cnv_summary_20260804.csv` (unzipped) | 70,156,911 | `9cefb776f8408a0db05a78d44b601a8ae9ba9b3906661703d92fa1053b188cdd` |
| `mutations_summary_20260724.zip` | 584,662 | `f8f97a7301ca604d69ffb1ed30d7d386901803f6a9bcda012084e0e6b1b44dce` |
| `mutations_summary_20260724.csv` (unzipped) | 2,241,637 | `998344a90cb7d9fc0e052c967a064e03b8e03eb80d5a73e6d9be3e09b76f06ef` |

**One byte of discrepancy, recorded rather than smoothed over.** The Passports
API reports `GDSC1_fitted_dose_response_27Oct23.xlsx` as 30,348,856 bytes; the
file that arrives is **30,348,857**. The other six sizes match the API exactly.
The archive opens cleanly (9 zip entries, as GDSC2 has) and parses to 333,161
rows. This is recorded as an upstream metadata quirk. **The SHA-256 above is
what the script should be held to**, not the API's number.

## What is NOT downloaded here, and why

**DepMap expression, lineage and PRISM are already on disk**, reached by the
existing symlink `data/raw/depmap`. See `data/raw/depmap_README.md` for their
provenance, checksums and the browser-only acquisition route. This analysis adds
no DepMap file and changes none.

| Input | Where | Release |
|---|---|---|
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | `data/raw/depmap/` (symlink) | 26Q1 |
| `Model.csv` - lineage **and `SangerModelID`** | `data/raw/depmap/` (symlink) | 26Q1 |
| `Repurposing_..._Primary_Data_Matrix.csv` | `data/raw/depmap/` (symlink) | 24Q2 |
| `Repurposing_..._Primary_Compound_List.csv` | `data/raw/depmap/` (symlink) | 24Q2 |

**`OmicsCNGene.csv` was NOT fetched, and RB1 is taken from Cell Model Passports
instead.** The declared RB1 confounder could come from either source. Passports
was chosen for two reasons, and the second is the stronger:

1. The DepMap portal needs a browser (Cloudflare); Passports does not.
2. **Passports keys on `model_id`, which IS `SANGER_MODEL_ID` - the GDSC join
   key.** Taking RB1 from DepMap would route the confounder through the same
   ACH-to-SIDM join as everything else and lose lines where that join fails.
   Passports joins to GDSC directly, with no intermediate key.

RB1 loss is therefore defined on two Passports tables: `cnv_summary` for
deletion (`cn_category`) and `mutations_summary` for truncating variants. The
exact rule is fixed in the script header, not here.

## Which BH3 mimetics are actually in the release

**Checked before the script was written, as the declaration requires**
(`docs/2026-09-28_bh3_mimetic_declaration.md` section 4). Counts are curve-fit
rows, i.e. lines screened, before any join or lineage filter.

| Compound | Declared role | GDSC2 | GDSC1 | PRISM 24Q2 | `PUTATIVE_TARGET` |
|---|---|---|---|---|---|
| **Navitoclax** (ABT-263) | **PRIMARY** | **967** | **928** | yes | BCL2, BCL-XL, BCL-W |
| **Venetoclax** (ABT-199) | **specificity control, predicted null** | **958** | **905** | yes | BCL2 |
| **WEHI-539** | BCL-XL selective | **962** | - | - | BCL-XL |
| ABT-737 | BCL-XL / BCL-2 / BCL-W tool compound | 957 | - | yes | BCL2, BCL-XL, BCL-W |
| AZD5991 | MCL1 selective | 713 | - | yes | MCL1 |
| Sabutoclax | pan-BCL-2-family | 892 | - | - | BCL2, BCL-XL, BFL1, MCL1 |
| Obatoclax Mesylate | pan, **not declared; added as found** | 785 | 862 | yes | BCL2, BCL-XL, BCL-W, MCL1 |
| S63845 | MCL1 selective | **absent** | **absent** | yes | MCL1 |
| A-1331852 | BCL-XL selective | **absent** | **absent** | **absent** | - |
| A-1155463 | BCL-XL selective | **absent** | **absent** | **absent** | - |

**Three things this table settles, and they were settled before any model was
fitted.**

1. **The primary and its control are both present, in both GDSC screens**, so
   the specificity contrast does not rest on one screen.
2. **`WEHI-539` is a second, cleanly BCL-XL-selective probe in GDSC2**, which
   the declaration asked for and did not assume. Navitoclax is a triple
   inhibitor; WEHI-539 is not. **They are not two looks at the same compound**,
   and the declaration's prediction applies to both.
3. **`A-1331852` and `A-1155463` are absent from all three platforms.** Recorded
   so that their absence from the results is not read as a failure.

### Two naming traps, both real and both in the data

- **GDSC1 spells venetoclax `Venotoclax`** (`DRUG_ID` 412, MGH site). GDSC2
  spells it `Venetoclax` (`DRUG_ID` 1909). A script matching one spelling
  silently loses the specificity control from an entire screen.
- **GDSC2 spells ABT-737 as `ABT737`**, with no hyphen, and it is **not** in
  `screened_compounds_rel_8.5.csv` under any grep for `ABT-737`. It is
  `DRUG_ID` 1910 and is found only in the dose-response table itself.

Both are why the script greps `unique(DRUG_NAME)` in the loaded tables and
records what it matched, rather than filtering on assumed spellings.

## The join, measured before the script was written

`GDSC$SANGER_MODEL_ID` to `Model.csv$SangerModelID`. **Not on cell line names.**

| | |
|---|---|
| `Model.csv` rows | 2,154, of which **1,220 carry a `SangerModelID`** |
| GDSC2 lines | 969 |
| **joining to DepMap** | **967, or 99.8%** - clears the declared 90% stop |
| of those, **solid** (not `Lymphoid`, not `Myeloid`) | **799** - clears the declared ~300 stop |
| of those, **Breast** | **52** - as the declaration predicted, underpowered |
| navitoclax lines that both join and are solid | **797** |

These are cohort counts, not estimates. The script re-derives every one of them
and asserts the two stop rules rather than trusting this table.

## Read-time traps

- **The GDSC tables are one row per curve fit**, keyed
  `NLME_CURVE_ID`, not one row per line. A line screened against a drug in both
  GDSC1 and GDSC2 appears in both files; `DATASET` distinguishes them.
- **`LN_IC50` is natural log of IC50 in micromolar**, and is **not** bounded by
  the screened range - fits extrapolate beyond `MAX_CONC`. `AUC` is bounded
  `[0, 1]`. **Both are lower = more sensitive.**
- **`Z_SCORE` is computed across the whole GDSC panel including haematological
  lines** and is not used here: this analysis fits its own models on a
  solid-only population, and a z-score standardised on a different population
  would smuggle that population back in.
- **The Passports `cnv_summary` is one row per model per gene**, 70 MB
  unzipped. Read only the `RB1` rows; do not load it whole into a data frame
  beside the expression matrix.
- **`cn_category` is a string**, with `Loss` and `Deletion` as distinct values.
  They are not interchangeable and the script's RB1 rule names which it uses.
