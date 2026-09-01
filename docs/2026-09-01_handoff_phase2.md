---
date: 2026-09-01
status: disposable - delete once used
purpose: resume mid-phase-2. Everything is run EXCEPT E09.
---

# Handoff - phase 2 in progress

`github.com/gszabadkai/myc_human_exploratory` (private), `main` level with
`origin`, tree clean.

**Every script has been written and run except `E09`.** `results/` holds twelve
artefacts, `outputs/figures/` twenty-five figures. The findings are in four
notes and several of them have been **corrected or reinterpreted** since they
were first written - see section 3, and read the corrections, not just the
headings.

---

## Resume prompt

Paste into a fresh session started in `/Users/gs/code/myc_human_exploratory`:

```
Read CLAUDE.md, then docs/2026-09-01_handoff_phase2.md.

This repo is EXPLORATORY. Nothing here is pre-registered. The sibling repo
myc_human_validation is a COMPLETED pre-registered study, frozen at d3ac60e,
not reopened, and its CLAUDE.md's no-post-hoc rule does NOT apply here.

E00-E08 are all written and run. E09 (Spearman against Pearson, bicor, Kendall
and a non-monotonicity probe) is written and has NEVER been run - that is where
I stopped and what I want to look at in detail.

I source scripts in Positron; you write and edit them, you do not run them.
```

---

## 1. THE NAMING CONTRACT - read this before touching anything

Phase 1's biggest error was reading an 18-member MYC panel as one kind of
object. It was not: `FELSHER_61` had been stripped of every MitoCarta gene by
the validation study and the other 17 had not, and nothing in any name said so.

**Every scored MYC estimator now carries an explicit suffix. There are no bare
names.** Defined in `E00`:

| suffix | definition |
|---|---|
| `__FULL` | the set exactly as distributed |
| `__MITOSTRIP` | minus `MITOCARTA_ALL` (1,136 genes) |
| `__PROLIFSTRIP` | minus HALLMARK E2F + G2M (327 genes) |
| `__BOTHSTRIP` | minus both |

18 signatures x 4 + 4 CollecTRI variants + `log2MYC` + `M_c_call` = **78
estimators**. Use the E00 constants rather than string literals:

```r
MYC_REF        <- "FELSHER__MITOSTRIP"        # the validation study's M_a
MYC_REF_FULL   <- "FELSHER__FULL"
MYC_LOW_ENTANG <- "MYC_UP.V1_UP__FULL"        # 1.5% proliferation
MYC_HALLMARK   <- "HALLMARK_MYC_TARGETS_V1__FULL"
MB_REF         <- "M_b__MITOSTRIP"
```

**And the second half of the same defect:** `PROLIF_DISJOINT` is `PROLIF_STD`
minus the 9 proliferation genes of `FELSHER_61` and nothing else. It is disjoint
from `M_a` **alone** and shares 97 genes with the FULL panel union. Never report
a proliferation-ADJUSTED rho without the `__PROLIFSTRIP` value beside it - E1
shows they disagree in sign.

---

## 2. What is on disk

```
scripts/   E00 setup   E01 load+anchor   E02 score   E03 atlas   E03b F3 falsifiers
           E04 figures E05 death axis    E05b D1 falsifier
           E06 estimator anatomy  E07 mtDNA anatomy  E08 strata+genes
           E09 correlation measures   <- WRITTEN, NEVER RUN
functions/ correlation_engine.R   mitopps.R   strata.R
```

| results/ | built by | holds |
|---|---|---|
| `frames.rds` | E01 | harmonised covariates |
| `scanb_scores.rds`, `new_set_scores.rds`, `set_definitions.rds` | E02 | 110 sets scored per cohort, the labelled panel, `strip_refs` |
| `correlation_atlas.rds` | E03 | **290,190 cells** + summaries |
| `mtdna_falsifiers.rds` | E03b | the F3 verdict over 77 estimators |
| `celldeath_axis.rds` | E05 | 125,604 cells, overlay, expression-matched null |
| `d1_falsifier.rds` | E05b | the D1 contrast null |
| `estimator_anatomy.rds` | E06 | variant table, regulon halves |
| `mtdna_anatomy.rds` | E07 | residual structure, bicistronic control |
| `strata_and_death_genes.rds` | E08 | strata, drivers, the 44 canonical genes |

`data/from_validation/` (563 MB) is gitignored and not on origin; re-copy from
`myc_human_validation @ d3ac60e` and `E01` stops with the instruction if absent.

---

## 3. Findings, and which ones were corrected

**Read the corrections.** Four findings changed materially after their first
write-up, and the notes carry the original claim with the correction beside it.

### The plane - `docs/2026-08-31_phase1_atlas_findings.md`

| | claim | status |
|---|---|---|
| **F1** | the MYC-OXPHOS correlation is not proliferation | **RESOLVED and strengthened.** With every MitoCarta and E2F/G2M gene removed, 16/18 (TCGA) and 15/18 (SCAN-B) signatures still give rho >= 0.2, median 0.427 / 0.416. Removing proliferation genes *raises* rho for the strongest signatures. |
| **F2** | MYC mRNA carries no OXPHOS signal | stands |
| **F3** | the two genomes diverge | **REFRAMED.** The mtDNA programme is FLAT, not suppressed - it falls behind a rising nuclear programme, furthest behind the mitoribosome. All three falsifiers survive over 77 estimators. **`MT-ND6` dropped as a technical artefact.** The mitoPPS framing is demoted. |
| **F4** | gsva and zmean agree at 0.994 | stands - there are three instruments, not four |
| **F5** | present in every subtype, weakest in ER-negative | extended by S4 |
| **F6** | the mitoribosome beats OXPHOS | stands, and F3 says it is the arm mtDNA falls furthest behind |
| **F7** | fatty-acid oxidation is the only negative arm | stands |

### The death axis - `docs/2026-09-01_phase1_celldeath_findings.md`

| | claim | status |
|---|---|---|
| **D0** | most of the death axis is mitochondrial content | stands - the frame for everything else |
| **D1** | pro-survival-skewed apoptotic transcriptome in OXPHOS-high tumours | **NUMBER STANDS, INTERPRETATION DOES NOT.** See S6. |
| **D2** | `TANG_FERROPTOSIS` fails its expression-matched null | stands - do not report its 0.398 |
| **D3** | BCL2 family splits by function | **QUALIFIED.** The MYC column is pooling: `BCL2` is -0.369 pooled and **-0.009 in LumA**. The claim belongs on the OXPHOS axis, which is stratum-stable. |
| **D4/D5** | CICD unscoreable; product terms empty | stand |

### Phase 2 - `docs/2026-09-01_phase2_estimator_findings.md`

| | |
|---|---|
| **E1** | adjusting for proliferation and removing proliferation genes disagree in sign. For `FELSHER` - the one signature the covariate is disjoint from - adjusting lowers rho 0.404 to 0.344 while stripping raises it to 0.416. |
| **E2** | internal coherence measures proliferation content (0.641/0.618), not quality. **`ELLWOOD` should be dropped** (agreement 0.202 against a next-lowest 0.485); `ALFANO` flagged unreliable. |
| **E3** | **`M_b` is the activated half of its own regulon, and that is the weak half.** Activated 736 genes give +0.246/+0.100 on OXPHOS; repressed 72 give **-0.518/-0.578**; M_b correlates 0.955 with the former and 0.068 with the latter, outnumbered 10.2:1. A balanced contrast gives 0.595/0.604. |
| **E4** | four FULL signatures contain a BCL2-family gene, so those E05 overlay cells are not independent. `M_b__MITOSTRIP` contains none. |
| **M1** | **the bicistronic control settles CO1 vs CO2.** Genes sharing one mRNA differ by <= 0.158; `CO1` vs `CO2` differ by 0.735/0.526 - 4.7x and 3.6x. Post-excision, not transcriptional. |
| **M2** | the per-gene mtDNA deviation replicates at 0.888; no positional or abundance gradient |
| **M3** | the reproducible residual axis is **PC2, not PC1** (loadings replicate 0.629 against -0.168). It is ~3.5% of gene-level variance - reproducible, not large. |
| **S1** | **BCL2 vs MYC is between-subtype.** Every stratum weaker than pooled; LumA -0.009. |
| **S2** | D1 is not self-overlap - only 3 of 1,086 genes are in the OXPHOS arm, and deleting all 79 MitoCarta genes makes the contrast stronger |
| **S3** | D1's magnitude lives in the tails: deleting the top 100 of ~500 per side removes 84%/72% |
| **S4** | the ER-negative gap is not range restriction and is arm-specific - **the mitoribosome shows none** (0.045/0.013) while FAO and ROS carry 0.19-0.32 |
| **S5** | `Luminal` exceeds both LumA and LumB because combining adds between-group variance |
| **S6** | **the canonical machinery reverses D1's direction.** Among the 44 `family_pathway` genes: pro-death median +0.048, pro-survival -0.225. The split that predicts sign is whether the gene acts at the mitochondrion (0.453), not its direction of effect (0.225). |

---

## 4. NEXT: E09, which is where I stopped

`E09_correlation_measures.R` is written, verified, and **has never been run.**
It answers "is Spearman hiding anything?":

- **spearman** (the study's default) against **pearson** on the log scale,
  **bicor** (Tukey biweight - robust like a rank method but keeps magnitude) and
  **kendall** on a 1,200-sample subsample, over 220 pairs
- a **non-monotonicity probe**: the R-squared gain of a natural spline over a
  straight line, on ranks, plus decile profiles that show the shape

**How to read the disagreements**, which is the point of the design:

- pearson departs from spearman, bicor sits WITH spearman -> a heavy-tail
  artefact, and Spearman is the safer reading
- bicor sits WITH pearson against spearman -> magnitude carries real information
  that ranks discard
- a large **spline gain** with a modest spearman -> a threshold or saturating
  curve, which is the thing a rank correlation structurally cannot report

Verified before commit: bicor holds at 0.615 where two leverage points pull
pearson to 0.79; the spline gain returns 0.79 on a U-shaped control and **0.012**
on the headline MYC-OXPHOS pair, which already suggests that relationship is
essentially monotone and Spearman is losing nothing there.

**Why Spearman was chosen, for the record:** it is invariant to monotone
transforms, which is what makes a log GSVA score and a linear mitoPPS score
comparable at all. CLAUDE.md names the log-vs-linear question as the most likely
silent error in this repo. That is a strong reason and E09 is a sensitivity
check on it, not a challenge to it.

---

## 5. Open, in rough priority order

1. **Run E09** and read it in detail. That was the stopping point.
2. **Score `COLLECTRI_MYC_STIM`** (739 genes, activation-only) as a fifth base in
   `E02`. It is already in the snapshot, has never been scored, and E3 predicts
   it lands near the activated half at about +0.25 / +0.10 on OXPHOS. Cheapest
   test of the sharpest phase-2 finding.
3. **Drop `ELLWOOD`** from the panel and recompute F1's entanglement slope; flag
   `ALFANO` rather than dropping it.
4. **Re-annotate the 44 canonical death genes from an independent source**
   (Reactome apoptosis sub-pathways) and check S6's module medians hold. `APAF1`
   at -0.452 is the named anomaly.
5. **The ER-negative fatty-acid-oxidation reversal** (-0.29 / -0.16 against ~0 in
   ER-positive) is the largest arm-level stratum effect in the study and nothing
   has been done with it.
6. A **stranded, total-RNA (non-polyA)** dataset. It resolves `MT-ND6` outright
   and tests whether `CO1`/`CO2` survives without polyA selection. Phase 1 and 2
   cannot do it.

**Still out of scope until explicitly reopened:** MCbiclust / forkscale,
survival, treatment, METABRIC, DepMap, causal or mediation modelling, and
anything that revisits the validation study's hypotheses.

---

## 6. Traps found the hard way in phase 2

Each of these cost a debugging cycle. They are in the scripts as comments; they
are here so the next session does not rediscover them.

1. **Stale saved artefacts.** `E07` reads `E03b`'s output and keys it by
   estimator name. After the relabelling that file still said `FELSHER_61` and
   the failure surfaced three sections later as a rename error. E07 now checks
   its dependency in section 1. **When an estimator name changes, sweep
   `results/` for it.**
2. **`slice(-seq_len(0))` returns NO rows.** `seq_len(0)` is `integer(0)`, so
   the negation is too. Use `filter(row_number() > k)`.
3. **A message that disagrees with the code beneath it is worse than no
   message.** Several printed `FELSHER_61` while filtering on `MYC_REF`.
   Interpolate the constant.
4. **Do not write prose that prejudges a test the script then runs.** E06's
   header asserted that mitochondrial contamination deflates F1; its own section
   3.1 showed it does not. The header now states both outcomes.
5. **`frac_prolif` of a `__PROLIFSTRIP` variant is 0 by construction.** Ordering
   or correlating against it is a correlation with a constant. Use
   `frac_prolif_full`, carried from the base signature's FULL variant.
6. **`thin` must be per cohort.** A set defined with 15 genes can land 14 in one
   matrix and 15 in the other. The atlas carries `thin_in_cohort`.

---

## 7. Rules that still bind

- **Option A.** Claude Code writes and edits the numbered scripts. It does not
  run them. Infrastructure - git, snapshots, provenance READMEs, `CLAUDE.md`,
  planning and result notes - it may execute directly.
- **Exploratory posture.** Nothing is pre-registered. Report structure,
  gradients and cross-cohort reproducibility; never a p-value plucked from a
  290,190-cell grid. When something looks real, write down what would falsify it
  *before* the next analysis - `E03b` and `E05b` are what that looks like.
- **Never write to `myc_human_validation` or `myc_mouse`.** Read-only via
  `git -C ... show <ref>:<path>`.
- **Human only.** The ortholog tripwire is narrowed to calls:
  `grep -rnE "(mouse_to_human|human_to_mouse|ortholog[s]?)[[:space:]]*\\(" scripts/`
  must return nothing.
- **Scale discipline.** GSVA wants log VST; mitoPPS wants linear DESeq2. They
  never share an input object. The correlation engine rank-transforms
  everything, so correlations are immune - the discipline binds where scores are
  *built*.
- **Never pool GSVA or mitoPPS values across cohorts.** Compare correlations and
  patterns.
- R: no `print(n = X)` after `head()`; always `dplyr::count()`; ASCII-only
  strings; every numbered script ends with an `if (FALSE) { ... }` sandbox.
- **Gene lists selected on a statistic are descriptions, not findings**, and
  must say so wherever they appear.
