# Snapshot: the validation study's G2 result (MYC co-amplification)

One read-only artefact copied from `myc_human_validation`: the saved result of
its gate **G2**, copy-number co-occurrence with MYC amplification in TCGA-BRCA.
It is here for one use, Fig. 4x, by the author's ruling of 2026-11-05
(`docs/2026-11-05_E54_coamplification_panel_ruling.md`, `8fde82a`), which
amends decision 5.

**Nothing is read across repos at runtime.** Scripts here read this copy, never
the source.

## Provenance

| item | value |
|---|---|
| file | `g2_cnv_cooccurrence.rds` (10,279 bytes) |
| md5 | `33da579ae1777c6bd246f5e8518d967d` |
| source | `/Users/gs/code/myc_human_validation/results/g2_cnv_cooccurrence.rds` |
| written | 2026-08-28 11:16:10 (its `$generated`; file time kept by `cp -p`) |
| written by | the validation repo's `scripts/05_gate_cnv_cooccurrence.R`, committed at `77968ac` (2026-08-28 11:11), before the result, and tracked at `d3ac60e` |
| its result note | `docs/2026-08-28_G2_result.md`, committed at `402b164` (2026-08-28 11:35), tracked at `d3ac60e` |
| copied | 2026-11-05, by Claude Code, read only from the source |

**Why the provenance is a checksum and a note, not a commit.**
- The validation repo ignores `results/`, so this file is in none of its
  commits, including the frozen `d3ac60e`. Its working tree holds the only
  other copy.
- So this file is **tracked here** (10 KB), unlike `data/from_validation/`.
  That makes it the only versioned copy anywhere.
- **What anchors it is the committed result note.** Every number that note
  prints was checked against this file on 2026-11-05 at the note's own
  precision: **46 of 46 match**.
  - section 2: the three primary tests (n, joint count, unadjusted and
    adjusted OR, interval, p);
  - section 4: MCL1 at any gain (0.503, 0.353 to 0.719);
  - section 5: the eight PAM50 odds ratios, both Basal intervals and the
    four Breslow-Day p values;
  - section 7: the covariate swap to fraction genome altered.

## What G2 pre-registered, and what it did not

From the G2 design note (`docs/2026-08-28_D4_scoping_and_G2_design.md`,
sections 2.3 to 2.9) and the validation audit
(`docs/handoff-2026-10-21_validation_audit.md`, in this repo):

- **Primary, pre-registered:** the three pooled tests at each gene's
  threshold, in ISAR calls, adjusted for `ANEUPLOIDY_SCORE`:
  - `MCL1` amplification (`+2`);
  - `BCL2L1` gain (`>= +1`);
  - `BBC3` loss (`<= -1`).

  The numeric pass rule was committed (`77968ac`) before the result
  (`402b164`). **Disclosed there:** crude odds ratios were seen before the
  thresholds were fixed.
- **Pre-specified alternative covariate:** `FRACTION_GENOME_ALTERED`. Both
  covariates are pre-specified, so **BCL2L1's result is reported under both,
  never as a bare odds ratio**: 1.77 (1.26-2.47) against 1.36 (0.984-1.88).
- **Secondary, declared but exploratory:** the PAM50 and TNBC strata, the
  Breslow-Day tests, the gene x threshold grid and the `BAX` regional
  control. The design says they are "secondary, BH-adjusted within family and
  labelled exploratory".

## Contents

`$freq_grid`, `$cooccurrence` (pooled tests, the grid, sources and covariates),
`$stratified` (PAM50 and TNBC strata, aneuploidy-adjusted), `$homogeneity`
(Breslow-Day), `$enrichment`, `$bbc3_deep`, `$verdict`, `$thresholds`,
`$samples`, `$generated`.

**n = 1,043 tumours** (ISAR, primary tumour samples). This repo's
`data/from_validation/tcga_brca_covariates.rds` holds the same calls for the
1,040 of them that also have expression. A recomputation there would not
reproduce these numbers exactly.

## The caveats G2 carries (its result note, section 9)

1. `MCL1` inverts at any gain (`>= +1`, OR 0.503). The association is with
   high-level amplification only.
2. `BCL2L1`'s pass depends on the covariate (above).
3. The basal enrichment of co-altered tumours is partly definitional, since
   MYC amplification is itself basal-enriched.
4. `BBC3` loss within Basal may be the 19q13 regional effect.

## Rules

- **Never patch this file.** A refresh is a new copy with a new md5 recorded
  here, and a dated note saying why.
- **Never write to the validation repo** (CLAUDE.md).
