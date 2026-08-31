# Cell-death gene sets, human

The mouse arm's cell-death sets, **taken from their native human columns**. Two
sources, both curated, both carrying first-class human symbols.

## Why this is a legitimate human input, and not an ortholog round trip

The validation study (`myc_human_validation`) deliberately **refused** the
mouse Apoptosis-PRO/ANTI sets. Its script 07 section 5.1 gave two reasons, and
neither survives here:

1. *"those lists contain `Bbc3` and `Bcl2l1`, whose human orthologs are the
   numerator and denominator of PRIME"* - circularity, **only** because apoptotic
   priming was that study's endpoint. Here the cell-death programme **is the
   object of study**, so their presence is the point rather than a defect.
2. *"they are a MOUSE-curated gene set. Translating them by uppercasing symbols
   is exactly the ortholog projection CLAUDE.md forbids"* - **checked against the
   files, and false for both sources below.** Each carries a human symbol column
   and human Ensembl ids as first-class data. The validation repo's own rule
   applies: *"Check the sheet, not the repository."*

**No ortholog function is applied to anything in this directory, in either
direction.** Source A's upstream README says in as many words: *"Orthology done
via MyGene.info + HomoloGene (~93%), NOT biomaRt. Do NOT remap."* The mouse arm
took the `mouse_symbol` column from the same table; this study takes
`human_symbol`. That is not a round trip.

## Source A - `cell_death_genes_consolidated.csv`

- Upstream: `/Users/gs/G/data/MK_myc_2022/myc_mouse`
- Pinned commit: **`6a9c7dd513800a2a433934314a87d161ce98caa2`** (branch
  `paper-final`, 2026-08-18)
- Blob: `270d459115bcc2993344fc65d5691a1693173620`
- Path upstream: `data/cell_death_genes_consolidated.csv`
- Read read-only via `git -C <mouse repo> show`; the mouse repo is **not**
  attached to sessions and is never written to.

Curated from GO / KEGG / Reactome / MSigDB. **1,232 rows, 1,232 distinct
`human_symbol`, none missing.**

| Column | Values | Use here |
|---|---|---|
| `human_symbol`, `human_ensembl` | | **the symbols used; mouse columns ignored** |
| `effect` | pro-death 512, pro-survival 587, ambiguous 42, unclassified 91 | the pro/anti-death axis |
| `pathway` | apoptosis 1,215, CICD 13, both 4 | the CICD axis |
| `is_mitochondrial` | 90 TRUE | the `_MITO` strata |
| `is_core` | | reported, not used to filter |
| `family_pathway` | BCL2 family 7, BH3-only 7, effector caspase 3, apoptosome, death receptor, ... | gene-level overlay labels |
| `modality` | sparse (17 genes) | superseded by source B |
| `confidence`, `evidence_score`, `in_GO/KEGG/Reactome/Hallmark` | | reported, not tuned |

**The eight CDC strata** (`pathway == "both"` contributes to both, as the mouse
arm did), in human symbols:

```
                              all    MITO
pro-death     apoptosis       502      39
pro-survival  apoptosis       584      41
pro-death     CICD             13       5
pro-survival  CICD              4       2
```

**CICD is thin, and that is the axis of most interest.** At the emission floor
of n >= 5 (the upstream library's own rule) two strata cannot be scored:
`pro-survival CICD` (4) and `pro-survival CICD MITO` (2). They are reported as
named gene lists instead. This is why the mouse arm emitted six `CDC_*` sets and
not eight. `ambiguous` and `unclassified` effects are excluded from scores and
reported separately.

## Source B - `tang_modalities/` - the 15 regulated cell-death modalities

- Upstream: `/Users/gs/G/data/MK_myc_2022/mammary_geneset_library`
- Pinned tag: **`v1.0`** = `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b` (2026-06-20)
- Path upstream: `data/raw/from_myc_mouse/cell-death/*.csv`
- Primary reference: **Tang et al. 2024, doi:10.1016/j.csbj.2024.08.012**

Human symbols with `ENSG` ids, a PMID and a curator comment per gene.

> The library's own loader (`R/05_load_myc_mouse_assets.R`) records that its
> README describes these as "Originally human, mouse-converted via
> ortholog_table" but that *"inspection of the file content shows the gene
> column still holds HUMAN symbols (uppercase, with ENSG human Ensembl IDs in
> gene_id). The CSV content is the source of truth."* That is an independent
> confirmation of the same finding, made upstream and for the opposite purpose
> (they map human -> mouse; we keep human).

```
Autophagy_dependent_cell_death  1195     Necroptosis      83     Entotic_cell_death   17
Ferroptosis                      935     Pyroptosis       54     Disulfidptosis       16
Apoptosis                        610     Lysosome_dep     40     Alkaliptosis         15
                                         Immunogenic      34     Parthanatos          11
                                         MPT_driven_necr  33     Oxeiptosis           10
                                         Cuproptosis      27     NETotic_cell_death    9
```

**Two cautions, both structural.** `Autophagy_dependent_cell_death` (1,195) and
`Ferroptosis` (935) are 6.5% and 5% of the expression matrix - a correlation
with either is close to a correlation with general expression, and needs a
size-matched comparator before it is believed. At the other end `Parthanatos`
(11), `Oxeiptosis` (10) and `NETotic_cell_death` (9) fall below the n >= 15 GSVA
floor and are carried as gene lists, not scores.

## Rules

- **Do not edit in place.** Re-take from the pinned commits and bump the SHAs.
- **Do not remap.** No `human_to_mouse`, no `mouse_to_human`, no ortholog table.
- Sizes above are asserted on load. A changed count means the upstream file
  moved and the snapshot must be re-taken, not patched.
