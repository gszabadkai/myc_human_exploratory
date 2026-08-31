# Curated human metabolic gene list (`GS_metabolic`)

Snapshotted 2026-08-28 for the script 07 specificity work. Small enough to
commit, so this one is on `origin`.

| | |
|---|---|
| File | `gs_metabolic_human.csv` |
| Rows | 2,347 genes, one per row |
| Columns | `gene_name`, `gene_symbol`, `entrez_id`, `classification` |
| Classifications | 74 |
| SHA-256 | `891427debea8403a24e15a5d43eb309c4b32481ef18309bd21bd5d6606c50bb8` |

## Provenance

Extracted from the **`human` sheet** of
`data/raw/user_curated/GS_metabolic_genes_list.xlsx` in the geneset library:

```
repo   /Users/gs/G/data/MK_myc_2022/mammary_geneset_library
tag    v1.0
path   data/raw/user_curated/GS_metabolic_genes_list.xlsx
blob   21f79d9a04ae88dbeff8794fd5a4f34a51809047
```

Retrieved read-only with `git -C <repo> show v1.0:<path>`, per CLAUDE.md. Rows
with a missing symbol or classification were dropped (2,348 -> 2,347);
whitespace trimmed. Nothing else was changed.

## This is human-native, and that is the point

The xlsx has two sheets, `human` and `mouse orthologs`. The library reads the
**human** sheet and maps it *to* mouse (`R/04_load_user_curated_sets.R`, which
requires columns `Gene Symbol` and `Classification` from the `human` sheet).

So the human version is the **original**, not a round trip. This is not the
mouse-derived `outputs/gmt/human/` tree that CLAUDE.md rejects - that tree is
mouse-native sets pushed through `mouse_to_human()`, gitignored and unpinned by
the tag. This is a tracked raw input carrying its own human data. Different
objects; see the amended wording in CLAUDE.md.

## What it is for, and what it is NOT for

**NOT the specificity battery.** That is MitoCarta, matching the mouse arms
named in `myc_mouse/scripts/43_substrate_specificity_and_tradeoff.R`. See
`docs/2026-08-28_specificity_panel_proposal.md`.

**NOT the per-complex sets either. CORRECTED 2026-08-28.** This file was
snapshotted to supply complex-level resolution "which MitoCarta does not
provide". **That claim was wrong.** MitoCarta 3.0 Sheet 4 carries a full
hierarchy - `OXPHOS > Complex N > CN subunits / CN assembly factors` - giving
`CI subunits` 44, `CII` 4, `CIII` 10, `CIV` 21, `CV` 21, nesting exactly inside
`OXPHOS subunits`. Script 07 uses those.

**And the substitute was not merely redundant, it was wrong.** Read the members,
not the class name:

```
ATPase (n = 5)  =  ATP13A1  ATP13A2  ATP13A3  ATP13A4  ATP13A5
```

Those are P5-type cation-transporting ATPases - ER and lysosomal transporters,
`ATP13A2` being the Kufor-Rakeb Parkinson gene. **Not one is a subunit of the
mitochondrial ATP synthase.** Used as "Complex V" in script 07's first run, that
arm correlated 0.04 with OXPHOS subunits and -0.04 with Complex IV. The real ATP
synthase subunits are inside `Proton Transport` (49), mixed with the vacuolar
V-ATPase (`ATP6V*`) and the uncoupling proteins, so that class is not Complex V
either.

Three further problems in the same five sets: `Complex III` writes `UCRC` and
`UQCR` (legacy symbols for `UQCR10` and `UQCR11`) and so loses 20% of itself
against any current annotation; `Complex I` writes `NDUFA4L`; and `Complex IV`
writes its mtDNA subunits as `COX1/COX2/COX3`, which no `^MT-` strip catches.

**Do not repair those by resolving symbols in reverse.** `COX1` and `COX3`
resolve to `PTGS1` and `COX2` to `PTGS2` - the prostaglandin synthases,
pharmacology's COX-1 and COX-2 - which would inject two abundant inflammatory
genes into an OXPHOS set, invisibly. `PRODH2` resolves to `PRODH`, a different
gene. Script 07 leaves such symbols unresolved and names them in its
harmonisation report instead.

## What it IS for

**An independent second opinion, reported and never used.** Script 07 prints
this file's five complex classes against MitoCarta's on every run - the check
that would have caught the `ATPase` error before it was scored. The `ATPase` row
is expected to disagree completely, and that line is the permanent record of why
this file does not supply the sets.

It also carries the wider metabolic classifications (`Fatty Acid` 139,
`Krebs` 26, `Redox` 68, `Glutathione` 34, `Folate` 12, `Glycine` 7,
`Proton Transport` 49, `Ubiquinone` 13 ...), which are a different curation of
the same biology and are available if a specificity result needs a second view.
**Read the members before using any of them.**

## Read-time traps

- `classification` is a flat character column, not a hierarchy. `Complex I` and
  `Proton Transport` are siblings, and a gene appears **once**: the 74 classes
  partition the 2,347 genes rather than overlapping. Do not assume set overlap
  semantics from MSigDB or MitoCarta.
- Symbols are human HGNC as curated, and have not been re-harmonised against the
  MitoCarta Sheet-3 vocabulary the way G1 harmonised the Felsher and CollecTRI
  sets. **Run them through the same alias resolution before intersecting with
  anything**, or a deprecated symbol will silently drop (`KARS` -> `KARS1` cost
  G1 a gene until it was caught).
- Every gene appears exactly once, so class sizes sum to 2,347. A gene present
  in the expression matrix but absent here is not "not metabolic" - the list is a
  curation, not an exhaustive annotation.
