# DepMap / CCLE - reached by SYMLINK, not held here

Provenance for the DepMap files consumed by `scripts/E17_depmap_ox_dependency.R`
(B4). **This repo does not hold these files.** `data/raw/depmap` is a symlink:

```
data/raw/depmap  ->  /Users/gs/code/myc_human_validation/data/raw/depmap
```

Created 2026-09-04. The target is ~1.2 GB, is **gitignored in that repo**, and
is not on any `origin`. The symlink itself is gitignored here - git would
otherwise track it as a link object that resolves only on this machine - so
**this file is the only record this repo has of the input.**

## THE CONSEQUENCE, STATED PLAINLY

**This repo now has an input it cannot regenerate from its own records.** Every
other directory under `data/` carries a pinned SHA or a checksum and can be
rebuilt from what is written here. This one cannot: the re-download URLs, the
figshare file ids, the checksums and the release-note archaeology all live in
`myc_human_validation/data/depmap/README.md`, and the files themselves live in
that repo's gitignored `data/raw/`.

Three things follow, and none of them is optional:

1. **If `myc_human_validation` is deleted or its `data/raw/` is cleared, B4's
   input is gone from here too**, and recovering it means re-downloading from
   the portal with a browser (below).
2. **Nothing here writes to that repo.** It is CLOSED, pre-registered, and
   frozen at `d3ac60e`. Reading a pinned raw data file through a symlink is not
   reopening a study; writing anything to it would be. Not results, not tables,
   not a note.
3. **The script asserts the link resolves before it reads anything.** See below.

## Why an assertion and not a graceful skip

`scripts/14_depmap_dependency.R` in the validation repo responds to an absent
optional file by skipping the section and saying so, which is right for PRISM.
**A broken symlink looks exactly like an absent directory.** Without a guard,
B4 would report "PRISM skipped, files absent" and carry on having silently lost
the CRISPR arm as well - a run that completes and answers a different question.
So `E17` stops, in this order, on: the link missing; the link present but its
target unresolvable; or any of the three REQUIRED files unreadable through it.
PRISM stays a skip, as in script 14, because it is genuinely optional.

## Files, and their two release pins

| File | Release | Required | Used for |
|---|---|---|---|
| `Model.csv` | **26Q1** | yes | lineage assignment (`OncotreeLineage == "Breast"`) |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | **26Q1** | yes | MYC / OXPHOS scoring. **log2(TPM+1)**, unstranded |
| `CRISPRGeneEffect.csv` | **26Q1** | yes | Chronos gene effect |
| `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` | **24Q2** | optional | drug sensitivity |
| `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` | **24Q2** | optional | compound name -> column id |
| *superseded* `CRISPRGeneEffect_24Q4_figshare_SUPERSEDED.csv` | 24Q4 | never | see the byte guard below |

**The three 26Q1 files must come from ONE release.** Model IDs are stable across
releases but the line set and the Chronos scaling are not, and mixing releases
produces a join that looks complete and is not.

**PRISM is 24Q2 and that is deliberate, not a stale copy.** Repurposing has its
own release cadence, `Repurposing Public 24Q2` (May 2024) is still the latest,
and 26Q1 contains no drug-sensitivity data at all. It joins on `ModelID`, which
is stable, so models added since 24Q2 simply have no PRISM row: **the drug arm
loses n, it does not gain bias.** Do not "fix" the gap by downgrading the rest.

**Two expression variants exist in 26Q1**, stranded and unstranded. The
**unstranded** file is used - it is continuous with every earlier DepMap release
and with the published CCLE analyses this is read beside. Only the unstranded
file is on disk. Do not mix them; the one used is recorded in the saved object.

## The byte guard, kept from script 14

`CRISPRGeneEffect_24Q4_figshare_SUPERSEDED.csv` is the 24Q4 matrix fetched from
figshare before the portal release was checked. figshare never mirrored past
24Q4, so it is not a 26Q1 file. It is renamed rather than deleted - but **a
re-download under the original name would resurrect it silently, and nothing
inside the file records its release.** Byte length is the only signal:

```r
CRISPR_24Q4_BYTES <- 428678699   # stop if CRISPRGeneEffect.csv is exactly this
```

Verified through the symlink 2026-09-04: `CRISPRGeneEffect.csv` is 440,646,050
bytes, so the guard passes.

## Acquisition - the omics files need a BROWSER

Established in the validation repo on 2026-08-30 by enumerating the figshare API
and reading the 26Q1 release notes, and carried forward here unchanged:

- **figshare's latest DepMap mirror is 24Q4 (December 2024).** It never mirrored
  25Qx or 26Q1. For the current release the portal is the only source.
- **The DepMap portal is behind a Cloudflare challenge.** `curl` of any portal
  download URL returns a ~5 KB HTML verification page, not data; saved under a
  `.csv` name it fails much later and very confusingly. The matching
  `storage.googleapis.com` paths return 403. **Use a browser** -
  https://depmap.org/portal/data_page/?tab=allData , release **26Q1**, three
  files: `Model.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
  `CRISPRGeneEffect.csv`.
- **PRISM Repurposing IS scriptable**, from figshare `25917643`
  (DOI `10.25452/figshare.plus.25917643`), file ids pinned:

  ```
  curl -L -o Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv \
    https://ndownloader.figshare.com/files/46630984
  curl -L -o Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv \
    https://ndownloader.figshare.com/files/46630981
  ```

## Checksums, copied so this repo carries them

Recorded in `myc_human_validation/data/depmap/README.md` and reproduced here
because that is the whole point of this file. `shasum -a 256 data/raw/depmap/*.csv`

| File | Bytes | SHA-256 |
|---|---|---|
| `Model.csv` | 697,455 | `ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8` |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 305,007,605 | `0377be80c525fde98cbd2c6e8b06bdf2a4014a9683eb70182c1f8649d711021a` |
| `CRISPRGeneEffect.csv` | 440,646,050 | `e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e` |
| `Repurposing_...Primary_Data_Matrix.csv` | 72,456,953 | `3b6554cfc6c765af53088a676edc7bce00ee7d84fe808b93bbfa892de607bc3d` |
| `Repurposing_...Primary_Compound_List.csv` | 719,567 | `7e78f5901c1a97d2baab0789ab89832e716388da4eacaa9f094e7d2f2f5a3463` |
| *superseded* `CRISPRGeneEffect_24Q4...csv` | 428,678,699 | `3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b` |

Byte lengths re-verified through the symlink 2026-09-04. The SHA-256 values are
carried over, not recomputed here - recomputing 745 MB to confirm a copy is a
copy is not a good use of the check, and the byte lengths all match.

## Read-time traps

- **Column names are `SYMBOL (ENTREZ)`** in the expression and CRISPR matrices.
  Strip the Entrez suffix; keep the first occurrence of a duplicated symbol.
- **The 26Q1 release notes misdescribe their own expression file.** They say
  `ProfileID, is_default_entry, ModelID, then genes`. The shipped file begins
  `<unnamed index>, SequencingID, ModelConditionID, ModelID, IsDefaultEntryForMC,
  IsDefaultEntryForModel, then "SYMBOL (ENTREZ)"` - no `ProfileID`, no
  `is_default_entry`, and the flags are `"Yes"`/`"No"` strings, not logicals.
  Trust the file. Locate columns by name with fallbacks; decide what is a gene
  by testing the column is numeric.
- **It is ONE ROW PER SEQUENCING RUN**, not per model. 1,775 rows collapse to
  1,719 models. Filter on **`IsDefaultEntryForModel`**, not
  `IsDefaultEntryForMC` - the latter is the model-*condition* default, a
  different and larger set.
- **`OmicsExpression...Logp1` is log2(TPM+1).** GSVA takes it as supplied;
  mitoPPS takes `2^x - 1`. Two objects that never meet.
- **Chronos sign:** 0 = no effect, -1 = median common essential. **More negative
  = more essential.**
- **PRISM log2 fold change:** more negative = more sensitive. Same direction.
- **The PRISM matrix is COMPOUNDS x CELL LINES**, the transpose of the other
  two. Read without transposing, the line intersection is empty and the drug arm
  reports "0 lines" rather than failing.

## Two declared deviations, carried forward verbatim and not re-argued

- **mitoPPS here runs on linear TPM, not linear DESeq2-normalised counts.**
  Declared in script 14's header. TPM is additionally length-normalised; mitoPPS
  is a COMPOSITION measure built from all-pairwise pathway ratios and is
  deliberately robust to total content, so the deviation is defensible - but it
  is a deviation. It makes the standing rule bite harder, not less: **CCLE
  mitoPPS values are NEVER comparable to TCGA or SCAN-B mitoPPS values. Only the
  pattern transfers.**
- **There is no expression-matched null in CCLE.** The TCGA nulls are not
  transferable. The 18-arm panel is a **rank ordering, not a calibrated
  p-value**, and a positive OXPHOS result is not reportable until that null is
  built in CCLE.

## What is NOT taken from the validation repo

Its `results/depmap_dependency.rds` - the Block G output - is part of the frozen
record and is neither read nor rewritten. B4 forks the *script* as a template
(read read-only with `git -C /Users/gs/code/myc_human_validation show
HEAD:scripts/14_depmap_dependency.R`) and asks a different question of the same
raw files: Block G extracted only the `MYC:OX` interaction and reported the `OX`
main effect nowhere.

Its `ARM_SETS` are also not reused. B4 rebuilds the 18 MitoCarta arms from
**this** repo's pinned `data/mitocarta_human/` so they match what `E16` used.
