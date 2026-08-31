# MYC signature compendium, human

The mouse arm scored MYC on a **compendium of signatures**, not one. This is
that compendium, in its original human form.

## Provenance

- Upstream: `/Users/gs/G/data/MK_myc_2022/mammary_geneset_library`
- Pinned tag: **`v1.0`** = `cbd8f16d2b0f95c5d4e86bed6aa112e42538a34b` (2026-06-20)
- Paths upstream: `data/raw/from_myc_mouse/myc_signature_genesets.gmx`,
  `data/raw/from_myc_mouse/felsher_integrative_signature.csv`
- Snapshot date: 2026-08-31

The upstream README describes the GMX as *"Felsher MYC compendium (**original
human GMX**)"*, and the file content agrees - uppercase human symbols
throughout. **Native human; no ortholog mapping is applied here.**

## Contents

`myc_signature_genesets.gmx` - **16 signatures**. GMX layout: row 1 is set
names, row 2 is the description line (all `NA` here), rows 3+ are genes, one
column per set, padded with `NA`.

Union 977 genes; **244 appear in >= 2 signatures and only 47 in >= 4.** The
signatures disagree substantially about what a MYC target is, and that
disagreement is measurable rather than a nuisance.

**Proliferation entanglement** - share of each signature that also sits in
`HALLMARK_E2F_TARGETS` + `HALLMARK_G2M_CHECKPOINT`, the covariate D7 was written
about. Measured 2026-08-31 against the TCGA matrix:

```
                              n   in matrix   % proliferation
MYC_UP.V1_UP                198       169        1.5    <- least entangled
BILD_MYC_ONCOGENIC_SIGNATURE 204      185        3.9
ALFANO_MYC_TARGETS          238       229        5.0
ELLWOOD_MYC_TARGETS_UP       13        13        7.7
PID_MYC_PATHWAY              24        24        8.3
MUHAR_MYC_SIGNATURE         100        90        9.0
MENSSEN_MYC_TARGETS          53        52       13.2
dang_myc_core_genes          51        50       13.7
SCHUHMACHER_MYC_TARGETS_UP   79        75       13.9
  [Felsher-61, the validation study's M_a]      14.8
COLLER_MYC_TARGETS_UP        25        24       16.0
DANG_REGULATED_BY_MYC_UP     71        70       16.9
PID_MYC_ACTIV_PATHWAY        78        77       19.2
HALLMARK_MYC_TARGETS_V2      57        57       19.3
DANG_MYC_TARGETS_UP         143       138       20.3
HALLMARK_MYC_TARGETS_V1     200       199       23.5
jung_myc_activity_signature  18        18       33.3
YU_MYC_TARGETS_UP            42        40       47.6    <- half proliferation
```

**A 30-fold range.** The validation study's `M_a` sits at 14.8% and the
compendium brackets it on both sides. Whether "MYC correlates with OXPHOS"
survives the low-entanglement signatures is the sharpest available test of
whether that correlation is MYC or is proliferation wearing MYC's name.

`felsher_integrative_signature.csv` - the integrated signature the validation
study's `M_a` derives from (67 genes raw, 61 after MitoCarta stripping; see
`../genesets_from_library_human/README.md`).

`HALLMARK_MYC_TARGETS_V1` is **not** in the GMX - only V2 is. V1 comes from
`msigdbr` at runtime, human-native either way.

## Rules

- **Do not edit in place**; re-take from the pinned tag and bump the SHA.
- **Never report a MYC-OXPHOS correlation from a single signature.** Report the
  panel, ordered by entanglement, alongside `log2(MYC)` and the
  proliferation-adjusted estimate.
- **MYC mRNA is not MYC activity.** In TCGA `rho(log2(MYC), OXPHOS subunits) =
  -0.032` against `+0.388` for the activity signature. Both are reported; the
  gap between them is itself a result.
