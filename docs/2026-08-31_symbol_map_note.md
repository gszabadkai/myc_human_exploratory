# Note - MitoCarta's Synonyms column is not injective

2026-08-31. Found while running `E02` on SCAN-B; fixed in `E02` section 3.1.

## What stopped the run

```
4.1 symbol harmonisation
   matched   resolved unresolved
      3858         61        367
Error: ! the SCAN-B symbol map sends two inputs to the same row: QARS
```

## Why

SCAN-B is a 2014 UCSC build (CLAUDE.md trap 7), so post-2018 renames have to be
walked back through MitoCarta's `Synonyms` column. That column lists an alias
under **every** symbol it was ever associated with, not only under the gene that
was actually renamed:

```
QARS1 -> GLNRS|MSCCA|PRO2195|QARS
EPRS1 -> EARS|EPRS|GLUPRORS|HLD15|PARS|PIG32|QARS|QPRS
GARS  is listed under both GARS1 and GART
```

`EPRS1` carries the whole multi-tRNA-synthetase complex's abbreviations. So the
row `QARS` had two claimants and the guard fired, correctly.

The second, quieter defect was the cause. The candidate filter ended with
`setdiff(cand, genes)` - refuse any candidate that appears anywhere in the input
universe. Its purpose was to stop one gene stealing a row another input already
occupies. But it punishes a gene for **its own old name** being present
elsewhere in the universe:

- `HALLMARK_MYC_TARGETS_V1` spells the gene `EPRS1`; `MENSSEN_MYC_TARGETS`
  spells the same gene `EPRS`.
- `EPRS` is therefore in `genes`, so it was refused as a candidate for `EPRS1`.
- The only survivor was `QARS` - **the wrong gene**.

The exploratory universe is 4,286 symbols against the validation study's
narrower one, which is why this surfaced here and not there.

## The fix

A deprecated row belongs to exactly one gene. A candidate is admitted only when
**exactly one symbol of the gene universe being mapped claims it**:

```r
cand <- cand[vapply(cand, function(cc)
  length(intersect(.alias_owners[[cc]], genes)) == 1L, logical(1))]
```

`QARS` (claimed by `QARS1` and `EPRS1`) and `GARS` (claimed by `GARS1` and
`GART`) go to neither claimant and are reported with a new `contested` status.
That is the conservative direction and it costs `QARS1` and `GARS1` their
SCAN-B rows - 1 gene each out of `CDC_PROSURVIVAL_APOPTOSIS` (584),
`TANG_APOPTOSIS` (608) and the sets carrying `GARS1`.

The collision guard was rewritten to match. A row reached by both a synonym and
its own current name is one gene spelled two ways and is now **named, not
fatal**; only a row reached by two different synonym-resolvers stops the run.
The same guard now runs on TCGA as well, which previously had none.

## Measured effect

| | old rule | new rule |
|---|---|---|
| resolved | 61 | 66 |
| contested | - | 2 (`GARS1`, `QARS1`) |
| collisions | 1 fatal (`QARS`) | 8 same-gene, named |
| disagreements with `scanb_pheno$symbol_map` over 468 shared symbols | **3** | **0** |
| `OXPHOS subunits` coverage | 0.9888 | 0.9888 |

The middle row of that table is the point. The old rule had *already* silently
lost `H2AZ1 -> H2AFZ`, `POLR1G -> CD3EAP` and `VARS1 -> VARS`, for exactly the
same reason: the wider universe contained their old spellings. E02 would have
stopped at the validation-map agreement check even if `QARS` had not fired
first. **The agreement check is load-bearing and should stay.**

The 8 same-gene merges are harmless - every set is intersected against the
matrix, so a set naming a gene twice carries it once:

```
AES <- TLE5    DFNA5 <- GSDME   IARS <- IARS1   VARS <- VARS1
CD3EAP <- POLR1G   EPRS <- EPRS1   H2AFZ <- H2AZ1   KARS <- KARS1
```

TCGA is unaffected: modern symbols, 3,912 matched, 0 resolved, 0 contested,
0 collisions. The map is a genuine no-op there.

## Still true, and deliberate

The map is **forward only** (current symbol -> its alias). The reverse would let
`COX1`/`COX2`/`COX3` resolve to the prostaglandin synthases. So an input spelled
with an old name in a modern matrix - `EPRS` in `MENSSEN_MYC_TARGETS` against
TCGA - stays unresolved and drops from that set. E02 section 6's coverage table
is where that shows up; the `frac < 0.80` warning is the tripwire.
