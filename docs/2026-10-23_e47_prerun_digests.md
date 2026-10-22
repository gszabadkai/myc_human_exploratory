---
date: 2026-10-23
status: >
  E47 WRITTEN and NOT YET RUN BY THE AUTHOR. These digests are recorded and
  committed BEFORE that run. No result is reported here; the result note
  follows the author's verified run.
relates-to:
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (dcb58a7, addendum 3b5e221)
  - scripts/E47_mutual_adjustment.R
---

# E47: digests recorded before the author's run

## 1. What the dry runs established

**There were two dry runs, every output redirected to the session scratchpad.**
Nothing was written in the repo.
- **They agree on 22 of 22 saved elements.** This excludes `built`, and
  `fits_S6`, which is not digest-stable (R rule 5).
- **33 of 33 assertions pass.**
- **There are no warnings.**
- **Run time is about 25 seconds.**

**The gates, all as the declaration requires:**

- **Selection rules (4a.6, 4a.7).** Each rule, applied to E10's
  `component_cor` at full precision, selects exactly its declared genes with
  exactly their signature signs.
  - `BID`'s margin in S6 is 0.0104.
  - `BBC3`'s margin in S5 is 0.0019.
  - `BMF` is the one gene of the twelve outside MitoCarta 3.0.
- **G1.** E34's quadrant labels are threshold splits of the continuous scores
  read here, in all four cohort-by-instrument blocks. The analysis sets are
  E34's, 938 and 3,143.
- **G2.** No fitted predictor contains any of the twelve. `M_b__FULL` and
  `M_b__PROLIFSTRIP` each contain five, which is why they are excluded.
- **G3.** `rho(MYC, PROLIF)` reaches 0.821 (TCGA) and 0.838 (SCAN-B) for the
  HALLMARK estimator.
  - **That estimator's `m3` cells are therefore VOID in both cohorts**:
    reported, never read.
  - No other estimator reaches 0.80; `M_a` is 0.696 / 0.615.
- **G4.** S6, both arms and all twelve genes, rebuilt from expression by E34's
  recipe, equal E34's saved values with a maximum difference of **0**. S5 and
  S7 are built by the same function.
- **VIF.** All 34 terms above the `sqrt(5)` caveat are MYC and `PROLIF` in
  TCGA's HALLMARK `m3` fits, which are exactly the cells G3 voids. **Over every
  cell the rule reads, the maximum `GVIF^(1/(2*Df))` is 1.776**, so no read
  cell carries a collinearity caveat.
- **Section 9.** The bootstrap's design matrix reproduces the reported OX
  coefficients exactly. The replicates are seeded, `PROJECT_SEED` plus the
  cell index, and identical across runs.
- **Displays A, B and C.** Every drawn point is one stored row, at its stored
  value, read back from the built plot.

## 2. Digests, recorded before the author's run

| element | digest |
|---|---|
| `verdict` | `f67c29fbb4e7ffe29c38e20eec5cd161` |
| `labels` | `b258c3c5f7562952a93bdef3f7f0ac82` |
| `primary_cells` | `faad5db6e76d8ae677891b11f6da3a76` |
| `ladder` | `0e5ab5716fb5d8ab1404645aef7956ea` |
| `results` | `baab2c53714f7d63a3d6d41ff9298ada` |
| `ratio` | `804613baeeb84cef7f5c591fac94bf8e` |
| `d1_leakage` | `f55878898645e7ce36e1ab2722a7083c` |
| `d2_rungs` | `4cdf69f68653134e3f21e92c74ae44e7` |
| `bcl2l1_bbc3` | `cc0a336eef81da2db8f36104be7e4279` |
| `steeper_licensed` | `bb73ad91bcb7e948250d465016f7bea8` |
| `boot_replicates` | `742e145c0c21e24b79b8093fac391a6b` |
| `coefficient_matrix` | `a7b196a288a474a693653efcd6204a08` |
| `vcov` | `152c3276674cfea4f4364a750c8ab664` |
| `vif` | `8742f159b01689f45df81752afaf8dd8` |
| `fit_digests` | `7fe71626f17096f2ce19bece096b9f27` |
| `frame` | `357afc907ea94997b8fdb718dfcfe440` |
| `gates` | `cf26e1de03f2b2d2e692c92b10a5bdc4` |
| `displays` | `604b6aeeebdfc2a716267066d43f8316` |
| `legend_draft` | `7c37084896aff3e255389ca42999a686` |
| `assertions` | `5c616aa95bbd0d168be91c66e9af47f1` |
| `figure_files` | `3df7421a787730480cb6ed8c4b9c4afd` |
| `spec` | `013d85dd424e67433539257fe4f10c8e` |

## 3. Verifying the author's run

```r
x  <- readRDS("results/e47_mutual_adjustment.rds")
el <- setdiff(names(x), c("built", "fits_S6"))
sapply(el, function(k) digest::digest(x[[k]]))   # compare with section 2
```

**The tracked figures** are
`docs/figures/2026-10-23_E47_display_{A_genes,B_estimator_ladder,C_bcl2l1_bbc3}.{png,pdf}`.
Commit them after the run is verified.
