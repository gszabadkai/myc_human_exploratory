# Human MitoCarta 3.0

The human mitochondrial gene inventory and pathway hierarchy for this repo. This is
the authoritative MitoCarta source for every mitochondrial gene set used here.

**Do not edit in place. If a newer MitoCarta release is adopted, re-snapshot and update
this file rather than patching the workbook.**

## Provenance

- File: `Human.MitoCarta3.0.xls`
- Upstream: Broad Institute, https://www.broadinstitute.org/mitocarta
  (Rath et al. 2021, Nucleic Acids Research, MitoCarta3.0)
- Snapshot source on this machine: `/Users/gs/Downloads/Human.MitoCarta3.0.xls`
- Source file mtime: 2025-11-25 00:02:43
- Snapshot date: 2026-08-28
- Method: byte-for-byte copy, verified by MD5 before and after
- MD5: `3c0bd24e362238216e142bc708e41286`
- Size: 10,179,584 bytes

### Why the Downloads copy and not the myc_mouse one

A second, content-identical copy exists at
`/Users/gs/G/data/MK_myc_2022/myc_mouse_main/data/Human_MitoCarta3_0.xls`. The two files
differ by 11 bytes of Excel metadata and have different MD5s, but their parsed contents
are identical: same 1,136 symbols, same 154 MitoPathway rows, identical pathway gene
strings (verified 2026-08-28).

The Downloads copy was taken because it is the original Broad download and sits outside
the mouse repo, so this human input has no cross-repo dependency at all.

## Contents

Four sheets:

| # | Sheet | Rows | Note |
|---|---|---|---|
| 1 | `Description` | - | Release notes |
| 2 | `A Human MitoCarta3.0` | 1,136 | The inventory. 1,136 unique symbols, 13 `MT-` prefixed |
| 3 | `B Human All Genes` | - | Full genome background with scoring columns |
| 4 | `C MitoPathways` | 154 rows / **149 pathways** | Pathway hierarchy; columns `MitoPathway`, `MitoPathways Hierarchy`, `Genes`. The last 5 rows are blank padding |

Sheet 2 columns of interest: `HumanGeneID`, `MouseOrthologGeneID`, `Symbol`, `Synonyms`,
`Description`, `MitoCarta3.0_List`.

`readxl::read_xls()` emits benign column-type warnings on sheets 2 and 3. Suppress them
at read time; they do not affect the symbol or pathway columns.

## Why this file exists rather than a gene set from the library

`mammary_geneset_library` v1.0 ships an `outputs/gmt/human/by_category/01_mitocarta_human.gmt`.
**It must not be used here.** That file is Mouse MitoCarta 3.0 projected through an
ortholog map: the library's only MitoCarta input is `data/raw/Mouse.MitoCarta3.0.xls`,
and `R/14_export_gmts.R` builds the "human" tree by running `mouse_to_human()` over the
mouse sets. It carries 1,083 genes against the 1,136 in the real human inventory.

It has human symbols and mouse provenance, which is exactly the failure mode the repo
split exists to prevent. See the "THIS IS A HUMAN REPO" section of `CLAUDE.md`.

## Standing convention

The 13 mtDNA-encoded protein-coding genes (`MT-` prefix) are held in a separate synthetic
"mtDNA-encoded OXPHOS subunits" pathway and are never pooled with nuclear-encoded OXPHOS
subunits - expression-scale skew. See `CLAUDE.md`.

## Read-time traps (verified 2026-08-28)

- **Sheet 4 has a leading unnamed index column.** `readxl` names it `"2"`, so
  `pw[[1]]` is a row number, not a pathway name. Address columns by name:
  `MitoPathway`, `MitoPathways Hierarchy`, `Genes`. Sheet 4 is 154 x 4.
- **Sheet 4 ends with 5 entirely blank padding rows.** 154 rows, 149 real
  pathways. Drop them at load with `filter(!is.na(MitoPathway))`. Leaving them
  in corrupts results two ways, and only one is loud: `MitoPathway == "OXPHOS"`
  returns 1 TRUE and 5 NA, so logical subsetting yields 6 elements (loud), and
  `strsplit()` then turns each NA into a phantom `"NA"` gene, inflating every
  pathway by exactly 5 (silent). Use `which()` for lookups - it is NA-safe.
- **`Genes` is a single comma-space-separated string per pathway**, not a list column.
  Split on `", "` before use.
- **`MitoPathways Hierarchy` encodes the tree with `>` separators**, e.g.
  `"Mitochondrial central dogma > mtDNA maintenance > mtDNA replication"`. Parent and
  child rows both exist, so pathway sets overlap by construction - do not treat the 149
  pathways as disjoint.
- Relevant `MitoPathway` values for this project include `OXPHOS` (169 genes),
  `OXPHOS subunits` (102), `OXPHOS assembly factors`, and `Complex I` through
  `Complex V`. Note the plan's primary OXPHOS measure is the *subunit* set, not the
  `OXPHOS` umbrella, which also contains assembly factors. Those two counts are
  post-padding-row-removal; with the blank rows left in they read as 174 and 107.
- `readxl::read_xls()` emits benign column-type warnings on sheets 2 and 3.
