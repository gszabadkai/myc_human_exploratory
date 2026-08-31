# Human assets from mammary_geneset_library

Snapshot of the genuinely human-native assets from `mammary_geneset_library`.

**Consume as-is. Do not rebuild here, do not edit in place. If the library changes,
re-snapshot from a new tag rather than patching files.**

## Provenance

- Source repo: `mammary_geneset_library`, at
  `/Users/gs/G/data/MK_myc_2022/mammary_geneset_library`
- Tag: `v1.0` ("Gene set library v1.0 - used to finalise myc_mouse")
- Commit (dereferenced): `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b`
- Commit date: 2026-06-20
- Snapshot date: 2026-08-28
- Method: `git -C <library> show v1.0:<path>`, read-only. Nothing in the library was
  created, modified, or deleted, and no directory access was granted to it.

## Contents

```
data/genesets_from_library_human/
  README.md                            # this file
  felsher_integrative_signature.csv    # 67 genes, human ENSG + HGNC
```

### felsher_integrative_signature.csv

- Source path in the library: `data/raw/from_myc_mouse/felsher_integrative_signature.csv`
- Blob SHA at `v1.0`: `193a18f58cc9869acae6d2dcdcd70d2bed71e861`
- 67 rows, 67 unique symbols, no NA or blank
- Columns: unnamed ENSG rowname, `Symbol`, `P.Value`, `Adjusted.P.Value`,
  `Median.Pearson`, `Number.of.MYC.mouse.experiments.in.which.the.gene.is.differentially.expressed`

**This is a human file despite its directory name in the library.** All 67 rownames are
human Ensembl `ENSG` IDs and all 67 symbols are HGNC. The signature is derived from MYC
mouse experiments but is reported in human identifiers by its authors. The library's own
`R/05_load_myc_mouse_assets.R` records the same directory-name mislabelling for the Tang
cell-death CSVs and resolves it the same way: file content is the source of truth.

This is the M-a estimator input (plan section 7.1) and the subject of gate G1.

## What was deliberately NOT taken, and why

Recorded so the absence is not mistaken for an oversight.

**Nothing from `outputs/gmt/human/`.** The library ships a nine-category human GMT tree.
It was not snapshotted, for two independent reasons:

1. **It is mouse-derived.** The library README describes itself as "mouse-native ...
   parallel human-symbol versions", with "human export = HGNC via one-to-one ortholog
   mapping". `R/14_export_gmts.R` builds the tree by running `mouse_to_human()` over the
   mouse sets. Those files carry human symbols and mouse provenance. Loading them here
   would violate the "THIS IS A HUMAN REPO" rule while looking correct.
2. **The tag does not pin them.** `outputs/` is gitignored in the library, so no GMT is
   reachable by `git show` at `v1.0`. Only 64 files are tracked at that tag and none is a
   GMT. Anything taken from `outputs/` would be an unpinned copy of whatever is on disk
   that day - the same trap already documented in `data/from_myc_mouse/README.md` for the
   mouse CSVs.

Specifically:

- **`01_mitocarta_human.gmt` is not used.** It is Mouse MitoCarta 3.0 ortholog-projected,
  1,083 genes against 1,136 in the real human inventory. Human MitoCarta 3.0 is
  snapshotted separately at `data/mitocarta_human/`. See that README.
- **The Felsher set in `02_myc_signatures_human.gmt` is not used.** It is a
  human -> mouse -> human round trip and carries 62 genes against the 67 in the CSV above.
  The round trip drops `CD3EAP`, `CIRH1A`, `KARS`, `NPM1`, `RPLP0`, `VARS` and adds
  `EIF5AL1`, a mouse-paralog artefact. Losing NPM1 and RPLP0 from a MYC signature is not
  cosmetic, and the five-gene difference bears directly on the G1 decision threshold
  (~50 genes after MitoCarta stripping).

The remaining seven categories (mammary development, metabolism, proliferation, TF
targets, biogenesis discrimination, apoptosis, intersections) are mouse-derived on the
same basis. Script 07 needs a specificity panel that those categories would otherwise
supply; the source for it is an open question, not a settled one. Do not reach for the
library tree to fill it without deciding that explicitly.

## Reading a further library file later

Without granting directory access:

    git -C /Users/gs/G/data/MK_myc_2022/mammary_geneset_library show v1.0:<path>

This reaches tracked files only. `outputs/` and `results/` are gitignored in the library
and are not retrievable this way.
