---
date: 2026-11-06
status: >
  DECLARATION for E56's one new quantity, committed ALONE before any of it is
  computed: the quadrant means of the three effectors, BAX, BAK1 and BOK, on
  E34's model. EXPLORATORY and POST HOC: the effectors were never part of E34's
  configuration, and the author added them to the paper later (E49). No
  reading rule; nothing here is read as a finding. The values are drawn in
  Extended Data Fig. 6d beside E34's twelve genes.
relates-to:
  - docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md (dd5004a; question 2 asks for this fit, declared first)
  - docs/2026-10-01_e34_result.md and scripts/E34_quadrant_configuration.R (the model, repeated exactly)
  - scripts/E46_fig4_quadrant_heatmap.R (E46's panel, which Fig. 6d redraws)
  - docs/2026-11-03_e49_result.md (E49, which fitted the effectors on E47's model)
---

# E56 declaration: the effectors' quadrant means, for Fig. 6d

## 1. Why

**The author's ruling of 2026-11-05 (question 2)** asks Fig. 6d to carry Fig. 4g's fifteen genes in its four groups.

**E34, and E46 after it, computed quadrant means for twelve genes only.** The effectors BAX, BAK1 and BOK joined the paper later, with E49, and were never fitted by quadrant.

**So their quadrant means are a new quantity,** fitted here exactly as E34 fitted the other twelve, and declared before any is computed.

## 2. What is computed

**For each of BAX, BAK1 and BOK:**
- in TCGA-BRCA and in SCAN-B;
- with OXPHOS split by GSVA and by mitoPPS, as in E34's quadrants;
- **E34's m1 standardised group means** for the four quadrants, each with its standard error and 95% confidence interval.

That is 3 × 2 × 2 × 4 = **48 estimates.** Only m1 is fitted, because E46's heatmap, which Fig. 6d redraws, shows m1.

## 3. How: E34's path, verbatim

1. **Expression:** log2(linear DESeq2-normalised + 1), from the frozen snapshot (`data/from_validation/tcga_brca_linear.rds` and `scanb_linear.rds`, `d3ac60e`). SCAN-B's symbols resolve through `scanb_pheno.rds$symbol_map` with `functions/gene_matrix.R`'s resolver, as E34 resolves them (CLAUDE.md, trap 7).
2. **Tumours:** those of E34's saved analysis frames, in their order: TCGA-BRCA 1,095, SCAN-B 3,207.
3. **Scale:** Blom normal scores, `qnorm((rank - 3/8) / (n + 1/4))` scaled to unit variance, once per cohort over all its tumours (E34's `.blom()`).
4. **Model m1** (E34): `lm(gene ~ quadrant + PAM50 + purity + leukocyte fraction)` in TCGA-BRCA, and `lm(gene ~ quadrant + PAM50)` in SCAN-B, which has no purity or leukocyte estimate (trap 2). The covariates are E34's frames' own, already scored. Rows are complete cases, as E34 takes them.
5. **Group means** (E34's `.fit_groups()` and `.lincomb()`):
   - each quadrant's mean is the model's prediction averaged over its own rows, with the quadrant set to that quadrant;
   - its 95% interval comes from a′Va on the residual degrees of freedom.

   There is no product term.

## 4. Gates, which must pass before any effector is fitted

The script stops if any fails.
- **G1. The expression path reproduces E34's twelve genes.** Steps 1 to 3, applied to E34's twelve genes, reproduce the Blom scores in E34's saved frames, to 1e-12, in both cohorts.
- **G2. The model path reproduces E34's estimates.** Steps 4 and 5, applied to the twelve genes in E34's frames, reproduce E34's stored m1 group estimates (estimate, standard error and both limits) for both instruments and both cohorts, to 1e-12.
- **G3. The signature strip is E46's.** The same path reproduces the signature's m1 and m3 group estimates that E46 drew in its strip, to 1e-12.

## 5. What is not done

- **No reading rule, no label, no contrast, no test.** The effectors' means are drawn, not read.
- **No sentence may rest on them** without its own declaration.
- **Nothing else is refitted.** The twelve genes and the signature are E34's stored values, checked by G2 and G3.

## 6. What E49 leads one to expect, stated so that it can be wrong

E49 found BAX and BAK1 rising with OXPHOS level and with MYC activity in both cohorts, and BOK less so. On that, BAX's and BAK1's OXPHOS-high quadrants should sit above their OXPHOS-low ones, in both cohorts and on both instruments.

This expectation is written down so that a different pattern is seen as one. **It licenses nothing.**

## 7. What is saved (R rules 4 and 5)

- the 48 estimates, with the 192 gate values beside them;
- **the analysis frame:** E34's frames with the effectors' raw and Blom-scored columns added;
- **for each of the 12 effector fits** (3 genes × 2 cohorts × 2 instruments), the full coefficient matrix, the variance matrix, and a `fit_digests` table of `coef()` and `vcov()`. The fitted objects themselves are saved but not digested.

## 8. The display

**Fig. 6d is drawn by the same script, after the gates and the fit,** to the author's ruling (`dd5004a`):
- **Rows:** Fig. 4g's fifteen genes in its four groups.
- **Columns:** both OXPHOS splits at full width.
- **Colours:** 4g's colour code.
- **Keys:** the colour key reads "z-score"; the dot plot's key reads "Signature" and "proliferation adjusted".
