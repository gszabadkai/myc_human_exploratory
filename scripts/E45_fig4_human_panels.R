# E45_fig4_human_panels.R
# =============================================================================
# FIG. 4 HUMAN PANELS, PART 1: THE TWO OUTCOME FORESTS. DISPLAY ONLY.
#
# EXPLORATORY. Nothing here is pre-registered, and NOTHING HERE IS ESTIMATED:
# every plotted number is read from a saved object, exponentiated and laid out.
# No declaration, on the E15 precedent for display scripts.
#
# Specification: docs/2026-10-22_fig4_human_panels_proposal.md section 3, as
# amended by the decisions in its section 7. Authority on what a panel may
# assert: docs/2026-10-04_handoff.md section 2. Nothing here relaxes it.
#
# =============================================================================
# THE NARRATIVE OF RECORD, AND ITS LETTERING
# =============================================================================
# The author's Google Doc "11.8.26_MK_Myc_Paper" is the final draft narrative
# (author, 2026-10-22). It supersedes v5 in docs/, and its lettering is NOT
# v5's:
#
#   Doc      v5        content
#   4g       4g        configuration by quadrant          - NOT here; part 2
#   4h       ED 9      MYC / BCL2L1 / MCL1 co-amplification - NOT here; it is
#                      validation-study material (proposal decision 5)
#   4i       4h        OXPHOS and proliferation on pCR    - BUILT HERE
#   4j       4i        the same on METABRIC BCSS          - BUILT HERE
#
# FILES ARE NAMED BY CONTENT, NOT BY LETTER, so the next relettering does not
# reach a filename.
#
# =============================================================================
# DECISIONS THIS SCRIPT IMPLEMENTS (proposal section 7)
# =============================================================================
#   - AXIS: ratios on a log axis with 1 marked, in both panels (author,
#     2026-10-22). The axes are NOT shared: one is an odds ratio, the other a
#     hazard ratio. A direction cue sits under each axis, because on both
#     panels proliferation lies right of 1 and OXPHOS left of it, and a higher
#     OR is benefit while a higher HR is harm.
#   - METABRIC ROWS: the ER-positive (IHC) primary, its two declared
#     comparators (ER.Expr sensitivity, all patients in one model) and
#     ER-negative (IHC). NO GSE25066 DRFS ROW: declaration 15.4 excluded it, no
#     PROLIF term for it exists, and the Doc examines long-term outcome in
#     METABRIC alone.
#   - NO BRACKET. The proposal put the pooled pCR contrast on a bracket between
#     the pooled OXPHOS and proliferation markers. That bracket would SPAN THE
#     WRONG QUANTITY: the pooled contrast is the meta-analysis of the
#     per-cohort contrasts (15.4), not the ratio of the two pooled odds ratios,
#     and the two differ (section 6 prints both). The contrast goes in a
#     right-hand column on EVERY row instead, which is also what puts "the
#     cohorts agree about the size of the gap" on the panel.
#
# =============================================================================
# WHAT THE PANELS MAY AND MAY NOT ASSERT (2026-10-04_handoff.md section 2)
# =============================================================================
#   - The pCR panel shows OXPHOS m0 beside m1. An m1-only forest BREAKS THE
#     LICENCE: m0 fails the reading rule in both directions, and the licensed
#     clause's second qualification is that it is reported beside m1.
#   - The annotated quantity is the TESTED CONTRAST (E44), not the gap between
#     whiskers. 15.3 binds: a significant contrast makes NEITHER coefficient a
#     finding. The METABRIC ER-positive OXPHOS marker is the boundary cell
#     (excludes 1 by 3.7e-04, unreproduced by its own sensitivity) and NO OX
#     VALUE IS PRINTED ON EITHER PANEL; its upper bound rounds to 1.000.
#   - No s508 cell, no GSE25066 DRFS row, nothing from X01.
#   - The two panels are NEVER on a shared axis and their contrasts are never
#     set side by side: a ratio of odds ratios is not a ratio of hazard ratios.
#   - METABRIC does not replicate the pCR analysis (10.5): different endpoint,
#     treatment mix and spine. Nothing here may be laid out to imply it does.
#   - SINGLE INSTRUMENT: GSVA only. mitoPPS exists in none of these cohorts.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits, refits and pools NOTHING. Pooled values are E40's and E43's
#     stored meta-analyses; the summary type is CHECKED against section 9's rule
#     (random effects at I2 >= 50), not chosen here.
#   - It reads E40's $drfs element NOWHERE.
#   - It does not build 4g (part 2, waiting on proposal decisions 1 and 2) or
#     the Doc's 4h co-amplification panel (decision 5).
#   - It does not interpret. Its legends are drafts and are DESCRIPTIVE ONLY.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

message("\nE45: Fig. 4 human panels, part 1 - the two outcome forests\n",
        strrep("=", 78))
message("DISPLAY ONLY. Nothing is fitted; every value is read from a saved object.")

PATH_E39 <- file.path(DIR_RESULTS, "e39_prepare_and_gate.rds")
PATH_E40 <- file.path(DIR_RESULTS, "e40_ladder.rds")
PATH_E41 <- file.path(DIR_RESULTS, "e41_metabric_ladder.rds")
PATH_E43 <- file.path(DIR_RESULTS, "e43_proliferation_coefficient.rds")
PATH_E44 <- file.path(DIR_RESULTS, "e44_contrast_test.rds")
PATH_E45 <- file.path(DIR_RESULTS, "e45_fig4_human_panels.rds")

# From E31 on, the key figure is written twice: outputs/ as usual, and the
# tracked docs/figures/, which survives a fresh clone.
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_STEM_DOC    <- "2026-10-22_E45_outcome_forests"
PATH_FIG_PCR    <- file.path(DIR_FIGURES, "E45_pcr_forest")
PATH_FIG_MB     <- file.path(DIR_FIGURES, "E45_metabric_forest")
PATH_FIG_BOTH   <- file.path(DIR_FIGURES, "E45_outcome_forests")
PATH_FIG_DOC    <- file.path(DIR_DOC_FIGURES, FIG_STEM_DOC)

# =============================================================================
# 0. CONSTANTS - fixed by the specification, none chosen on the data
# =============================================================================
PCR_ORDER <- c("GSE194040", "GSE164458", "GSE25066")     # top to bottom
PCR_N     <- c(GSE194040 = 988L, GSE164458 = 482L, GSE25066 = 465L)
# 465 is the PRE-SPECIFIED GSE25066 set (2026-10-04_handoff.md 2.2). 490 and
# 508 are sensitivity sets and never appear in a panel.

MB_ORDER <- c("METABRIC / ER-pos (IHC)", "METABRIC / ER-pos (ER.Expr)",
              "METABRIC / pooled", "METABRIC / ER-neg (IHC)")
MB_ROLE  <- c("primary", "declared sensitivity", "declared comparator",
              "other stratum")
MB_LABEL <- c("ER-positive", "ER-positive", "All patients", "ER-negative")
MB_ERDEF <- c("IHC, ", "expression, ", "", "IHC, ")    # on the second line
# "METABRIC / pooled" is ALL PATIENTS IN ONE MODEL, not a meta-analysis. It is
# labelled "All patients" so the word "pooled" means one thing per figure.
MB_FORBIDDEN <- c("METABRIC / ER-neg (ER.Expr)")         # goes to ED Fig. 8d

# Row geometry, SHARED by both panels so their rows align side by side.
Y_ROWS  <- c(4.6, 3.6, 2.6, 1.2)    # three rows, a rule, then the fourth
Y_SEP   <- 1.9
Y_CUE   <- 0.28                     # the direction cue, below every interval
Y_HEAD  <- 5.45
Y_LIM   <- c(0, 5.9)
OFF_PCR <- c("OXPHOS, m0" = 0.25, "OXPHOS, m1" = 0, "Proliferation, m1" = -0.25)
OFF_MB  <- c("OXPHOS, m1" = 0.14, "Proliferation, m1" = -0.14)

# The marker key. One factor for both panels, so the collected legend merges.
KEY_LEVELS <- c("OXPHOS, m0", "OXPHOS, m1", "Proliferation, m1")
KEY_TEXT   <- c("OXPHOS, m0" = "OXPHOS, no proliferation term (m0)",
                "OXPHOS, m1" = "OXPHOS, with proliferation term (m1)",
                "Proliferation, m1" = "Proliferation (m1)")
# Okabe-Ito bluish green and vermillion. Validated with the dataviz palette
# checker on a light surface: CVD separation dE 11.0 (deutan), normal-vision
# 25.8, both >= 3:1 against white - all checks pass. Shape is a second,
# colour-independent encoding of the term. Vermillion is NOT reused for
# "OXPHOS-high" here; E34's quadrant figure is not part of this panel set.
COL_OX <- "#009E73"
COL_PR <- "#D55E00"
KEY_COL   <- c("OXPHOS, m0" = COL_OX, "OXPHOS, m1" = COL_OX,
               "Proliferation, m1" = COL_PR)
KEY_FILL  <- c("OXPHOS, m0" = "white", "OXPHOS, m1" = COL_OX,
               "Proliferation, m1" = COL_PR)
KEY_SHAPE <- c("OXPHOS, m0" = 21, "OXPHOS, m1" = 21, "Proliferation, m1" = 24)

# Sized for the bottom row of a 180 mm (7.1 in) figure page at 1:1, so no
# text falls below 5.5 pt when the figure is placed unscaled.
BASE_PT <- 6.5
TXT_MM  <- BASE_PT / ggplot2::.pt          # geom_text size equal to 6.5 pt
SUB_MM  <- 5.5 / ggplot2::.pt

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.fmt_n   <- function(x) format(x, big.mark = ",", trim = TRUE)
.fmt_ci  <- function(e, lo, hi) sprintf("%.2f (%.2f, %.2f)", e, lo, hi)

# Section 9's rule, the only place a summary type is derived. It is CHECKED
# against each object's stored summary_is, never substituted for it.
.summary_for <- function(I2) ifelse(!is.na(I2) & I2 >= 50, "random effects",
                                    "fixed effect")

# Every assertion is recorded as well as enforced, so the saved object shows
# what was checked and by how much it passed.
ASSERT <- list()
.assert <- function(check, ok, n = NA_integer_, max_dev = NA_real_,
                    detail = "") {
  ASSERT[[length(ASSERT) + 1L]] <<- tibble::tibble(
    check = check, n_compared = as.integer(n), max_abs_dev = max_dev,
    pass = isTRUE(ok), detail = detail)
  .stop_if(ok, "ASSERTION FAILED: ", check,
           if (nzchar(detail)) paste0(" (", detail, ")") else "", ". STOP.")
  invisible(TRUE)
}
.dev <- function(a, b) if (length(a)) max(abs(a - b)) else NA_real_

# =============================================================================
# 1. Inputs - five saved objects, nothing recomputed from data
# =============================================================================
message("\n1. inputs")

for (p in c(PATH_E39, PATH_E40, PATH_E41, PATH_E43, PATH_E44)) {
  .stop_if(file.exists(p), "absent: ", p, ". Re-source the script that writes ",
           "it; the order is in docs/2026-10-04_handoff.md 5.2.")
}
e39 <- readRDS(PATH_E39)
e40 <- readRDS(PATH_E40)
e41 <- readRDS(PATH_E41)
e43 <- readRDS(PATH_E43)
e44 <- readRDS(PATH_E44)
message("   read: E39 $gate, E40 $pcr and $meta, E41 $exposure_corr, ",
        "E43 $prolif_focused and $meta, E44 $pcr, $metabric and $meta")
message("   NOT read: E40 $drfs. No DRFS cell can reach a panel.")

# =============================================================================
# 2. The pCR panel's values (the Doc's Fig. 4i)
# =============================================================================
message("\n2. pCR: OXPHOS m0 and m1 (E40), proliferation m1 (E43)")

OXP <- e40$pcr %>%
  dplyr::filter(cohort %in% PCR_ORDER, rung %in% c("m0", "m1"))
.assert("pCR rows come from the pre-specified set only",
        all(OXP$set == "primary") && nrow(OXP) == 6L, n = nrow(OXP),
        detail = paste(unique(OXP$set), collapse = ", "))
n_chk <- OXP %>% dplyr::distinct(cohort, n)
.assert("pCR cohort sizes are the declared analysis sets (GSE25066 is 465)",
        nrow(n_chk) == 3L &&
          all(n_chk$n == PCR_N[n_chk$cohort]), n = nrow(n_chk),
        detail = paste(n_chk$cohort, n_chk$n, collapse = "; "))

PRF <- e43$prolif_focused %>% dplyr::filter(set %in% PCR_ORDER)
.assert("E43 holds one m1 fit per pCR cohort", nrow(PRF) == 3L, n = nrow(PRF))

# OXPHOS m1 is in BOTH objects. They are the same fits (E43 14.4 reproduced
# them with a discrepancy of 0), so the panel may take OX from E40 and PROLIF
# from E43 only if that still holds exactly.
x40 <- OXP %>% dplyr::filter(rung == "m1") %>%
  dplyr::arrange(match(cohort, PCR_ORDER))
x43 <- PRF %>% dplyr::arrange(match(set, PCR_ORDER))
.assert("OXPHOS m1 in E40 equals OXPHOS in E43's m1 fits, exactly",
        identical(x40$cohort, x43$set) &&
          all(x40$estimate == x43$estimate_OX) &&
          all(x40$ci_lo == x43$ci_lo_OX) && all(x40$ci_hi == x43$ci_hi_OX),
        n = 3L, max_dev = max(.dev(x40$estimate, x43$estimate_OX),
                              .dev(x40$ci_lo, x43$ci_lo_OX),
                              .dev(x40$ci_hi, x43$ci_hi_OX)))
.assert("E43's pCR fits are logistic (odds-ratio scale)",
        all(x43$scale == "odds ratio"), n = 3L)

PCR_COH <- dplyr::bind_rows(
  OXP %>% dplyr::transmute(
    row_key = cohort,
    key = ifelse(rung == "m0", "OXPHOS, m0", "OXPHOS, m1"),
    term = "OX", rung, log_est = estimate, log_lo = ci_lo, log_hi = ci_hi,
    n, events, summary = NA_character_, I2 = NA_real_,
    source = paste0("e40_ladder.rds$pcr [set primary, rung ", rung, "]")),
  PRF %>% dplyr::transmute(
    row_key = set, key = "Proliferation, m1", term = "PROLIF", rung = "m1",
    log_est = estimate_PROLIF, log_lo = ci_lo_PROLIF, log_hi = ci_hi_PROLIF,
    n, events, summary = NA_character_, I2 = NA_real_,
    source = "e43_proliferation_coefficient.rds$prolif_focused"))

# --- the pooled row -------------------------------------------------------
.pooled <- function(m, key, term, rung, source) {
  .stop_if(nrow(m) == 1L, "expected one pooled row for ", key)
  rule <- .summary_for(m$I2)
  .assert(paste0("pooled ", key, ": stored summary type follows section 9 ",
                 "(I2 ", round(m$I2, 1), ")"),
          identical(m$summary_is, rule), n = 1L,
          detail = paste0("stored '", m$summary_is, "', rule '", rule, "'"))
  re <- rule == "random effects"
  tibble::tibble(
    row_key = "Pooled", key = key, term = term, rung = rung,
    log_est = if (re) m$re_estimate else m$fe_estimate,
    log_lo  = if (re) m$re_ci_lo    else m$fe_ci_lo,
    log_hi  = if (re) m$re_ci_hi    else m$fe_ci_hi,
    n = NA_integer_, events = NA_integer_, summary = rule, I2 = m$I2,
    source = paste0(source, " [", if (re) "re_" else "fe_", "*]"))
}

m40 <- e40$meta
m43 <- e43$meta
.assert("E40's pooled OXPHOS m1 equals E43's pooled OXPHOS, exactly",
        m40$fe_estimate[m40$rung == "m1"] == m43$fe_estimate[m43$term == "OX"] &&
          m40$re_estimate[m40$rung == "m1"] ==
            m43$re_estimate[m43$term == "OX"], n = 1L,
        max_dev = abs(m40$fe_estimate[m40$rung == "m1"] -
                        m43$fe_estimate[m43$term == "OX"]))

PCR_POOL <- dplyr::bind_rows(
  .pooled(m40 %>% dplyr::filter(rung == "m0"), "OXPHOS, m0", "OX", "m0",
          "e40_ladder.rds$meta [rung m0]"),
  .pooled(m40 %>% dplyr::filter(rung == "m1"), "OXPHOS, m1", "OX", "m1",
          "e40_ladder.rds$meta [rung m1]"),
  .pooled(m43 %>% dplyr::filter(term == "PROLIF"), "Proliferation, m1",
          "PROLIF", "m1", "e43_proliferation_coefficient.rds$meta [PROLIF]"))

PLOT_PCR <- dplyr::bind_rows(PCR_COH, PCR_POOL) %>%
  dplyr::mutate(
    panel = "pCR", scale = "odds ratio",
    ratio = exp(log_est), ratio_lo = exp(log_lo), ratio_hi = exp(log_hi),
    y_row = unname(stats::setNames(Y_ROWS, c(PCR_ORDER, "Pooled"))[row_key]),
    y_mark = y_row + unname(OFF_PCR[key]),
    key = factor(key, levels = KEY_LEVELS)) %>%
  dplyr::arrange(y_row, dplyr::desc(y_mark))
.stop_if(!anyNA(PLOT_PCR$y_mark), "a pCR row has no position")
.pr(PLOT_PCR %>% dplyr::transmute(row_key, key, summary,
    OR = round(ratio, 3), lo = round(ratio_lo, 3), hi = round(ratio_hi, 3)))

# =============================================================================
# 3. The METABRIC panel's values (the Doc's Fig. 4j)
# =============================================================================
message("\n3. METABRIC: OXPHOS and proliferation, from E43's m1 fits")

MBF <- e43$prolif_focused %>% dplyr::filter(set %in% MB_ORDER) %>%
  dplyr::arrange(match(set, MB_ORDER))
.assert("METABRIC rows are exactly the four specified, in order",
        identical(MBF$set, MB_ORDER), n = nrow(MBF),
        detail = paste(MBF$set, collapse = "; "))
.assert("no DRFS, s508, s490, GSE25066 or ER-neg (ER.Expr) row in METABRIC",
        !any(grepl("s508|s490|DRFS|GSE25066", MBF$set)) &&
          !any(MBF$set %in% MB_FORBIDDEN), n = nrow(MBF))
.assert("E43's METABRIC fits are Cox models (hazard-ratio scale)",
        all(MBF$scale == "hazard ratio"), n = nrow(MBF))

PLOT_MB <- dplyr::bind_rows(
  MBF %>% dplyr::transmute(row_key = set, key = "OXPHOS, m1", term = "OX",
    rung = "m1", log_est = estimate_OX, log_lo = ci_lo_OX, log_hi = ci_hi_OX,
    n, events),
  MBF %>% dplyr::transmute(row_key = set, key = "Proliferation, m1",
    term = "PROLIF", rung = "m1", log_est = estimate_PROLIF,
    log_lo = ci_lo_PROLIF, log_hi = ci_hi_PROLIF, n, events)) %>%
  dplyr::mutate(
    panel = "METABRIC", scale = "hazard ratio",
    summary = NA_character_, I2 = NA_real_,
    source = "e43_proliferation_coefficient.rds$prolif_focused",
    ratio = exp(log_est), ratio_lo = exp(log_lo), ratio_hi = exp(log_hi),
    y_row = unname(stats::setNames(Y_ROWS, MB_ORDER)[row_key]),
    y_mark = y_row + unname(OFF_MB[key]),
    key = factor(key, levels = KEY_LEVELS)) %>%
  dplyr::arrange(y_row, dplyr::desc(y_mark))
.stop_if(!anyNA(PLOT_MB$y_mark), "a METABRIC row has no position")
.pr(PLOT_MB %>% dplyr::transmute(row_key, key, n, events,
    HR = round(ratio, 3), lo = round(ratio_lo, 4), hi = round(ratio_hi, 4)))
message("   ", .sub(paste0(
  "The ER-positive (IHC) OXPHOS upper bound is ",
  sprintf("%.4f", PLOT_MB$ratio_hi[PLOT_MB$row_key == MB_ORDER[1] &
                                     PLOT_MB$term == "OX"]),
  ", which rounds to 1.000 at three decimals. It is the boundary cell and NO ",
  "OXPHOS VALUE IS PRINTED ON EITHER PANEL.")))

# =============================================================================
# 4. The contrast column - E44's tested PROLIF - OX, one entry per row
# =============================================================================
message("\n4. the contrast column (E44)")

# E44's components must still be E43's numbers, exactly, before its contrast is
# printed beside E43's markers.
c_all <- dplyr::bind_rows(e44$pcr, e44$metabric)
f_all <- e43$prolif_focused %>% dplyr::filter(set %in% c_all$set)
cj <- dplyr::inner_join(c_all, f_all, by = "set")
.assert("E44's contrast components equal E43's coefficients, exactly",
        nrow(cj) == nrow(c_all) &&
          all(cj$beta_A == cj$estimate_PROLIF) &&
          all(cj$beta_B == cj$estimate_OX), n = nrow(cj),
        max_dev = max(.dev(cj$beta_A, cj$estimate_PROLIF),
                      .dev(cj$beta_B, cj$estimate_OX)))

m44  <- e44$meta
rule <- .summary_for(m44$I2)
.assert(paste0("pooled pCR contrast: stored summary type follows section 9 ",
               "(I2 ", round(m44$I2, 1), ")"),
        nrow(m44) == 1L && identical(m44$summary_is, rule), n = 1L)
re44 <- rule == "random effects"

CONTRAST_COL <- dplyr::bind_rows(
  e44$pcr %>% dplyr::filter(set %in% PCR_ORDER) %>%
    dplyr::transmute(panel = "pCR", row_key = set, delta, ci_lo, ci_hi,
                     scale, summary = NA_character_, I2 = NA_real_,
                     source = "e44_contrast_test.rds$pcr"),
  tibble::tibble(panel = "pCR", row_key = "Pooled",
                 delta = if (re44) m44$re_estimate else m44$fe_estimate,
                 ci_lo = if (re44) m44$re_ci_lo else m44$fe_ci_lo,
                 ci_hi = if (re44) m44$re_ci_hi else m44$fe_ci_hi,
                 scale = "log odds ratio", summary = rule, I2 = m44$I2,
                 source = paste0("e44_contrast_test.rds$meta [",
                                 if (re44) "re_" else "fe_", "*]")),
  e44$metabric %>% dplyr::filter(set %in% MB_ORDER) %>%
    dplyr::transmute(panel = "METABRIC", row_key = set, delta, ci_lo, ci_hi,
                     scale, summary = NA_character_, I2 = NA_real_,
                     source = "e44_contrast_test.rds$metabric")) %>%
  dplyr::mutate(
    ratio = exp(delta), ratio_lo = exp(ci_lo), ratio_hi = exp(ci_hi),
    reads_as = ifelse(scale == "log odds ratio", "ratio of odds ratios",
                      "ratio of hazard ratios"),
    label = .fmt_ci(ratio, ratio_lo, ratio_hi),
    y_row = ifelse(panel == "pCR",
                   unname(stats::setNames(Y_ROWS, c(PCR_ORDER, "Pooled"))[row_key]),
                   unname(stats::setNames(Y_ROWS, MB_ORDER)[row_key])))
.assert("every panel row has exactly one contrast",
        nrow(CONTRAST_COL) == 8L && !anyNA(CONTRAST_COL$y_row) &&
          !any(duplicated(CONTRAST_COL[c("panel", "row_key")])),
        n = nrow(CONTRAST_COL))
.assert("pCR contrasts are on the odds scale and METABRIC's on the hazard scale",
        all(CONTRAST_COL$scale[CONTRAST_COL$panel == "pCR"] == "log odds ratio") &&
          all(CONTRAST_COL$scale[CONTRAST_COL$panel == "METABRIC"] ==
                "log hazard ratio"), n = nrow(CONTRAST_COL))

# Within ONE model the contrast IS the ratio of the two markers' ratios. This
# holds for every per-cohort and per-stratum row and is asserted; it does NOT
# hold for the pooled pCR row, which is why there is no bracket (section 6).
wm <- dplyr::bind_rows(PLOT_PCR, PLOT_MB) %>%
  dplyr::filter(row_key != "Pooled", rung == "m1") %>%
  dplyr::select(panel, row_key, term, ratio) %>%
  tidyr::pivot_wider(names_from = term, values_from = ratio) %>%
  dplyr::inner_join(CONTRAST_COL %>% dplyr::select(panel, row_key, cr = ratio),
                    by = c("panel", "row_key"))
.assert("within each model, the contrast equals PROLIF ratio / OX ratio",
        nrow(wm) == 7L && .dev(wm$cr, wm$PROLIF / wm$OX) < 1e-12, n = nrow(wm),
        max_dev = .dev(wm$cr, wm$PROLIF / wm$OX))
.pr(CONTRAST_COL %>% dplyr::transmute(panel, row_key, reads_as, label,
                                      summary, I2 = round(I2, 1)))

# =============================================================================
# 5. OXPHOS against proliferation, Spearman, in all four cohorts
# =============================================================================
# The Doc says "positively but moderately correlated". The four values are
# already recorded: three in E39's gate and the result note, METABRIC's in
# E41's $exposure_corr and 2026-10-03_e41_handoff.md:106. Read, not recomputed.
message("\n5. rho(OXPHOS, proliferation), four cohorts (read, not recomputed)")

RHO <- dplyr::bind_rows(
  e39$gate %>% dplyr::filter(pair == "OX vs PROLIF_DISJOINT") %>%
    dplyr::transmute(cohort, n, rho, source = "e39_prepare_and_gate.rds$gate"),
  e41$exposure_corr %>% dplyr::filter(pair == "OX vs PROLIF") %>%
    dplyr::transmute(cohort = "METABRIC", n = nrow(e41$frame), rho = spearman,
                     source = "e41_metabric_ladder.rds$exposure_corr"))
.assert("one OXPHOS-proliferation correlation per cohort, four cohorts",
        nrow(RHO) == 4L && !anyNA(RHO$rho), n = nrow(RHO))
.pr(RHO %>% dplyr::mutate(rho = round(rho, 4)))
message("   ", .sub(paste0(
  "Range ", sprintf("%.3f", min(RHO$rho)), " to ", sprintf("%.3f", max(RHO$rho)),
  ". GSE25066's gate value is on its 508-patient gate set, not the 465 ",
  "fitted; E39 computed the gate before the analysis sets were cut.")))

# =============================================================================
# 6. Why there is no bracket on the pooled pCR row
# =============================================================================
message("\n6. the bracket the proposal specified, and why it is not drawn")

pool_or <- PCR_POOL %>% dplyr::filter(rung == "m1")
gap_of_pooled <- pool_or$log_est[pool_or$term == "PROLIF"] -
                 pool_or$log_est[pool_or$term == "OX"]
pooled_delta  <- CONTRAST_COL$delta[CONTRAST_COL$row_key == "Pooled"]
BRACKET <- tibble::tibble(
  quantity = c("distance between the pooled markers (NOT an estimate)",
               "pooled contrast (the meta-analysis of per-cohort contrasts)"),
  log_scale = c(gap_of_pooled, pooled_delta),
  ratio = exp(c(gap_of_pooled, pooled_delta)))
.pr(BRACKET %>% dplyr::mutate(log_scale = round(log_scale, 4),
                              ratio = round(ratio, 3)))
message("   ", .sub(paste0(
  "A bracket between the pooled markers spans ", sprintf("%.2f", BRACKET$ratio[1]),
  " on the axis while its label would read ", sprintf("%.2f", BRACKET$ratio[2]),
  ". The first is the difference of two separately pooled estimates, ",
  "random-effects for proliferation and fixed-effect for OXPHOS, which 15.4 ",
  "forbids as an estimate and which has no valid standard error. The contrast ",
  "is printed in a column instead.")))

# =============================================================================
# 7. The figures
# =============================================================================
message("\n7. figures")

ROWS_PCR <- tibble::tibble(
  row_key = c(PCR_ORDER, "Pooled"), y = Y_ROWS,
  label = c(PCR_ORDER, "Pooled"),
  sub = c(sprintf("%s / %s", .fmt_n(PCR_N[PCR_ORDER]),
                  .fmt_n(x40$events[match(PCR_ORDER, x40$cohort)])), NA),
  face = c("plain", "plain", "plain", "bold"))
# The pooled row's second line is the three I2 values in marker-key order,
# drawn with plotmath so the source stays ASCII.
I2_POOL <- PCR_POOL$I2[match(KEY_LEVELS, as.character(PCR_POOL$key))]
I2_EXPR <- sprintf("italic(I)^2*' %s%%, %s%%, %s%%'",
                   round(I2_POOL[1]), round(I2_POOL[2]), round(I2_POOL[3]))

ROWS_MB <- tibble::tibble(
  row_key = MB_ORDER, y = Y_ROWS, label = MB_LABEL,
  sub = paste0(MB_ERDEF, .fmt_n(MBF$n), " / ", .fmt_n(MBF$events)),
  face = ifelse(MB_ROLE == "primary", "bold", "plain"))

.lab_panel <- function(rows, head, i2_expr = NULL) {
  p <- ggplot2::ggplot(rows) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = y + 0.17, label = label,
                                    fontface = face),
                       hjust = 0, size = TXT_MM) +
    ggplot2::geom_text(data = rows %>% dplyr::filter(!is.na(sub)),
                       ggplot2::aes(x = 0, y = y - 0.21, label = sub),
                       hjust = 0, size = SUB_MM, colour = "grey35") +
    ggplot2::annotate("text", x = 0, y = Y_HEAD, label = head, hjust = 0,
                      vjust = 0.5, size = SUB_MM, fontface = "bold",
                      lineheight = 0.95)
  if (!is.null(i2_expr)) {
    p <- p + ggplot2::annotate("text", x = 0,
                               y = rows$y[rows$row_key == "Pooled"] - 0.21,
                               label = i2_expr, parse = TRUE, hjust = 0,
                               size = SUB_MM, colour = "grey35")
  }
  p + ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::theme_void(base_size = BASE_PT)
}

.con_panel <- function(con, head, bold_key) {
  con <- con %>% dplyr::mutate(face = ifelse(row_key == bold_key, "bold",
                                             "plain"))
  ggplot2::ggplot(con) +
    ggplot2::geom_text(ggplot2::aes(x = 0.5, y = y_row, label = label,
                                    fontface = face),
                       hjust = 0.5, size = TXT_MM) +
    ggplot2::annotate("text", x = 0.5, y = Y_HEAD, label = head, hjust = 0.5,
                      size = SUB_MM, fontface = "bold", lineheight = 0.95) +
    ggplot2::geom_hline(yintercept = Y_SEP, colour = "grey75", linewidth = 0.25) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::theme_void(base_size = BASE_PT)
}

.forest <- function(pf, xlim, breaks, x_title, cue, title, key_levels,
                    connect_m0 = FALSE, legend = TRUE) {
  p <- ggplot2::ggplot(pf) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = Y_SEP, colour = "grey75", linewidth = 0.25)
  if (connect_m0) {
    # m0 to m1 joins TWO RUNGS OF ONE TERM, not two quadrants, so E34's rule
    # against connectors does not apply.
    cn <- pf %>% dplyr::filter(term == "OX") %>%
      dplyr::select(row_key, rung, ratio, y_mark) %>%
      tidyr::pivot_wider(names_from = rung, values_from = c(ratio, y_mark))
    p <- p + ggplot2::geom_segment(
      data = cn, ggplot2::aes(x = ratio_m0, xend = ratio_m1,
                              y = y_mark_m0, yend = y_mark_m1),
      colour = COL_OX, linewidth = 0.3, alpha = 0.6)
  }
  lo_a <- 1 / 1.04
  hi_a <- 1.04
  lo_b <- xlim[1] * 1.10
  hi_b <- xlim[2] / 1.10
  p +
    ggplot2::geom_segment(ggplot2::aes(x = ratio_lo, xend = ratio_hi,
                                       y = y_mark, yend = y_mark, colour = key),
                          linewidth = 0.45, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = ratio, y = y_mark, colour = key,
                                     fill = key, shape = key),
                        size = 1.7, stroke = 0.5) +
    # direction cues: description, not interpretation, so the legend rule
    # allows them
    ggplot2::annotate("segment", x = lo_a, xend = lo_b, y = Y_CUE, yend = Y_CUE,
                      linewidth = 0.3, colour = "grey30",
                      arrow = grid::arrow(length = grid::unit(1.1, "mm"),
                                          type = "closed")) +
    ggplot2::annotate("segment", x = hi_a, xend = hi_b, y = Y_CUE, yend = Y_CUE,
                      linewidth = 0.3, colour = "grey30",
                      arrow = grid::arrow(length = grid::unit(1.1, "mm"),
                                          type = "closed")) +
    # anchored at the arrow's inner end and running outward, so the text never
    # crosses the line at 1 however narrow the panel
    ggplot2::annotate("text", x = lo_a, y = Y_CUE + 0.21, hjust = 1,
                      label = cue[1], size = SUB_MM, colour = "grey30") +
    ggplot2::annotate("text", x = hi_a, y = Y_CUE + 0.21, hjust = 0,
                      label = cue[2], size = SUB_MM, colour = "grey30") +
    ggplot2::scale_x_continuous(transform = "log10", limits = xlim,
                                breaks = breaks, labels = format(breaks),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = KEY_COL, breaks = key_levels,
                                 labels = KEY_TEXT[key_levels], drop = FALSE,
                                 name = NULL) +
    ggplot2::scale_fill_manual(values = KEY_FILL, breaks = key_levels,
                               labels = KEY_TEXT[key_levels], drop = FALSE,
                               name = NULL) +
    ggplot2::scale_shape_manual(values = KEY_SHAPE, breaks = key_levels,
                                labels = KEY_TEXT[key_levels], drop = FALSE,
                                name = NULL) +
    ggplot2::labs(x = x_title, y = NULL, title = title) +
    ggplot2::theme_classic(base_size = BASE_PT) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.3),
      axis.text.x = ggplot2::element_text(size = BASE_PT - 1, colour = "black"),
      axis.title.x = ggplot2::element_text(size = BASE_PT - 0.5),
      plot.title = ggplot2::element_text(size = BASE_PT + 0.5, face = "bold",
                                         hjust = 0.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.key.width = grid::unit(3, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)) +
    # The key is removed through guides(), not legend.position: the combined
    # figure's closing `& theme(legend.position = ...)` overrides the latter,
    # which drew METABRIC's key a second time in the second dry run.
    if (legend) NULL else ggplot2::guides(colour = "none", fill = "none",
                                          shape = "none")
}

# Axis limits and breaks are fixed per panel. They are not shared: OR and HR
# are different quantities.
XLIM_PCR   <- c(0.5, 3.5)
BREAKS_PCR <- c(0.5, 0.75, 1, 1.5, 2, 3)
XLIM_MB    <- c(0.7, 1.6)
BREAKS_MB  <- c(0.75, 1, 1.25, 1.5)

# An interval outside the axis would be DROPPED by ggplot with only a warning,
# which would remove data from a panel silently. Stop instead.
for (a in list(list("pCR", PLOT_PCR, XLIM_PCR), list("METABRIC", PLOT_MB, XLIM_MB))) {
  .assert(paste0(a[[1]], " panel: every interval lies inside the axis limits"),
          all(a[[2]]$ratio_lo >= a[[3]][1]) && all(a[[2]]$ratio_hi <= a[[3]][2]),
          n = nrow(a[[2]]),
          detail = sprintf("limits %.2f to %.2f, data %.3f to %.3f", a[[3]][1],
                           a[[3]][2], min(a[[2]]$ratio_lo), max(a[[2]]$ratio_hi)))
}

# Widths in inches. The two panels together fill a 180 mm page with a gap.
W_PCR <- c(0.78, 1.85, 0.97)
W_MB  <- c(0.92, 1.42, 0.97)
W_GAP <- 0.12
H_PANEL <- 2.45

.build_pcr <- function(key_levels, legend = TRUE) {
  .lab_panel(ROWS_PCR, "Cohort\nn / pCR", I2_EXPR) +
    .forest(PLOT_PCR, XLIM_PCR, BREAKS_PCR,
            "Odds ratio for pCR\nper within-cohort SD",
            c("lower pCR", "higher pCR"), "Pathological complete response",
            key_levels, connect_m0 = TRUE, legend = legend) +
    .con_panel(CONTRAST_COL %>% dplyr::filter(panel == "pCR"),
               "Proliferation / OXPHOS\nratio of ORs (95% CI)", "Pooled") +
    patchwork::plot_layout(widths = W_PCR)
}
.build_mb <- function(key_levels, legend = TRUE) {
  .lab_panel(ROWS_MB, "METABRIC\nn / BC deaths") +
    .forest(PLOT_MB, XLIM_MB, BREAKS_MB,
            "Hazard ratio for breast-cancer death\nper within-cohort SD",
            c("longer survival", "shorter survival"),
            "Breast-cancer-specific survival", key_levels, legend = legend) +
    .con_panel(CONTRAST_COL %>% dplyr::filter(panel == "METABRIC"),
               "Proliferation / OXPHOS\nratio of HRs (95% CI)", MB_ORDER[1]) +
    patchwork::plot_layout(widths = W_MB)
}

FIG_PCR <- .build_pcr(KEY_LEVELS)
FIG_MB  <- .build_mb(KEY_LEVELS[2:3])
# In the combined figure only the pCR panel carries the key: it holds all three
# markers, and METABRIC's two are a subset with the same encoding. Collecting
# two non-identical keys drew two legends in the first dry run.
FIG_BOTH <- patchwork::wrap_plots(.build_pcr(KEY_LEVELS),
                                  patchwork::plot_spacer(),
                                  .build_mb(KEY_LEVELS[2:3], legend = FALSE),
                                  nrow = 1,
                                  widths = c(sum(W_PCR), W_GAP, sum(W_MB))) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

# --- every plotted value is its stored value, read back from the build ----
# This reads what ggplot will DRAW, after the log transform, and matches it to
# the plotting frame by marker position. It catches a mislabelled row, a
# dropped point and a wrong scale, none of which a check on the frame alone
# would see.
.check_drawn <- function(fig, pf, label) {
  fp  <- fig[[2]]                                   # the forest sub-plot
  b   <- ggplot2::ggplot_build(fp)
  geo <- vapply(fp$layers, function(l) class(l$geom)[1], character(1))
  xe  <- vapply(fp$layers, function(l) {
    m <- l$mapping$xend
    if (is.null(m)) "" else rlang::as_label(m)
  }, character(1))
  k_pt <- which(geo == "GeomPoint")
  k_ci <- which(geo == "GeomSegment" & xe == "ratio_hi")
  .stop_if(length(k_pt) == 1L && length(k_ci) == 1L,
           label, ": cannot identify the point and interval layers")
  dp <- b$data[[k_pt]]
  dc <- b$data[[k_ci]]
  mp <- match(round(dp$y, 9), round(pf$y_mark, 9))
  mc <- match(round(dc$y, 9), round(pf$y_mark, 9))
  .assert(paste0(label, ": every marker drawn is one frame row, none dropped"),
          nrow(dp) == nrow(pf) && !anyNA(mp) && !any(duplicated(mp)) &&
            nrow(dc) == nrow(pf) && !anyNA(mc) && !any(duplicated(mc)),
          n = nrow(dp))
  dev <- max(.dev(dp$x, log10(pf$ratio[mp])),
             .dev(dc$x, log10(pf$ratio_lo[mc])),
             .dev(dc$xend, log10(pf$ratio_hi[mc])))
  .assert(paste0(label, ": drawn positions equal the stored values"),
          dev < 1e-12, n = nrow(dp) + 2L * nrow(dc), max_dev = dev)
  dev_log <- max(.dev(log(pf$ratio), pf$log_est),
                 .dev(log(pf$ratio_lo), pf$log_lo),
                 .dev(log(pf$ratio_hi), pf$log_hi))
  .assert(paste0(label, ": plotted ratios are exp() of the stored log values"),
          dev_log < 1e-12, n = 3L * nrow(pf), max_dev = dev_log)
}
.check_drawn(FIG_PCR, PLOT_PCR, "pCR panel")
.check_drawn(FIG_MB,  PLOT_MB,  "METABRIC panel")

# The contrast labels are re-derived from the stored log values, independently
# of the column built above.
lab_chk <- CONTRAST_COL %>%
  dplyr::mutate(re = .fmt_ci(exp(delta), exp(ci_lo), exp(ci_hi)))
.assert("every contrast label is the stored contrast, exponentiated",
        identical(lab_chk$label, lab_chk$re), n = nrow(lab_chk))
# No OXPHOS or proliferation value is printed: the contrast column is the only
# text layer carrying numbers, and it is built from CONTRAST_COL alone.

.ensure_dir(DIR_FIGURES)
.ensure_dir(DIR_DOC_FIGURES)
.save_fig <- function(fig, stem, w, h) {
  ggplot2::ggsave(paste0(stem, ".png"), fig, width = w, height = h, dpi = 300,
                  bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), fig, width = w, height = h,
                  device = grDevices::cairo_pdf, bg = "white")
}
W_BOTH <- sum(W_PCR) + W_GAP + sum(W_MB)
.save_fig(FIG_PCR,  PATH_FIG_PCR,  sum(W_PCR), H_PANEL)
.save_fig(FIG_MB,   PATH_FIG_MB,   sum(W_MB),  H_PANEL)
.save_fig(FIG_BOTH, PATH_FIG_BOTH, W_BOTH, H_PANEL)
.save_fig(FIG_BOTH, PATH_FIG_DOC,  W_BOTH, H_PANEL)
message("   combined width ", sprintf("%.2f", W_BOTH), " in (",
        round(W_BOTH * 25.4), " mm), height ", H_PANEL, " in")
message("   ", PATH_FIG_PCR, ".{png,pdf}")
message("   ", PATH_FIG_MB, ".{png,pdf}")
message("   ", PATH_FIG_BOTH, ".{png,pdf}")
message("   ", PATH_FIG_DOC, ".{png,pdf}  (tracked)")

# =============================================================================
# 8. Legend drafts - DESCRIPTIVE ONLY (v5's legend rule)
# =============================================================================
# No interpretation, no verdict, and no "control": PROLIF was declared as a
# covariate. The model contents are read from E43's own $prolif_focused$terms
# rather than typed, so they cannot drift from the fits.
message("\n8. legend drafts")

terms_pcr <- stats::setNames(x43$terms, x43$set)
terms_mb  <- unique(MBF$terms)
.stop_if(length(terms_mb) == 1L, "the METABRIC fits do not share one formula")
pool_txt <- PCR_POOL %>% dplyr::mutate(t = sprintf("%s %.0f%%, %s", key, I2,
                                                   summary))

LEGEND <- list(
  pcr = paste0(
    "Odds ratios for pathological complete response (pCR) per within-cohort ",
    "standard deviation of the OXPHOS subunit score (GSVA) and of the ",
    "proliferation score, in three neoadjuvant chemotherapy cohorts (",
    paste(sprintf("%s, n = %s, %s pCR", PCR_ORDER, .fmt_n(PCR_N[PCR_ORDER]),
                  .fmt_n(x40$events[match(PCR_ORDER, x40$cohort)])),
          collapse = "; "),
    "). Open circle, OXPHOS in a logistic model with receptor subtype and ",
    "treatment arm where each varies within the cohort (m0); filled circle, ",
    "the same model with the proliferation score added (m1); triangle, ",
    "proliferation in that m1 model. Pooled, inverse-variance meta-analysis, ",
    "random effects (DerSimonian-Laird) where I2 >= 50% and fixed effect ",
    "otherwise; I2 is given for OXPHOS m0, OXPHOS m1 and proliferation (",
    paste(pool_txt$t, collapse = "; "), "). Right, the ratio of the ",
    "proliferation to the OXPHOS odds ratio within each m1 model, its variance ",
    "taken from that model's covariance matrix; the pooled value is the ",
    m44$summary_is, " meta-analysis of the three within-model ratios (I2 ",
    sprintf("%.0f", m44$I2), "%). Bars, 95% confidence intervals."),
  metabric = paste0(
    "Hazard ratios for breast-cancer-specific death per within-cohort ",
    "standard deviation of the OXPHOS subunit score (GSVA) and of the ",
    "proliferation score in METABRIC, from Cox models with PAM50 subtype, ",
    "treatment group, OXPHOS and proliferation (m1). Rows: ",
    paste(sprintf("%s (n = %s, %s events)", c("ER-positive by immunohistochemistry",
                                              "ER-positive by expression (ER.Expr)",
                                              "all patients in one model",
                                              "ER-negative by immunohistochemistry"),
                  .fmt_n(MBF$n), .fmt_n(MBF$events)), collapse = "; "),
    ". Circle, OXPHOS; triangle, proliferation. Right, the ratio of the ",
    "proliferation to the OXPHOS hazard ratio within each model, its variance ",
    "taken from that model's covariance matrix. The endpoint, treatment and ",
    "model differ from those in the pCR panel, and the axes are not shared. ",
    "Bars, 95% confidence intervals."),
  model_terms = list(pcr = terms_pcr, metabric = terms_mb))
message("   pCR:\n   ", .sub(LEGEND$pcr))
message("\n   METABRIC:\n   ", .sub(LEGEND$metabric))
message("\n   model formulas, from E43: ")
.pr(tibble::tibble(set = c(names(terms_pcr), "METABRIC (all four)"),
                   terms = c(unname(terms_pcr), terms_mb)))

# =============================================================================
# 9. Save - the plotting frames, per R rule 4
# =============================================================================
message("\n9. saving")

ASSERTIONS <- dplyr::bind_rows(ASSERT)
.pr(ASSERTIONS %>% dplyr::transmute(check, n_compared,
                                    max_abs_dev = signif(max_abs_dev, 3), pass))
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")

PLOT_FRAME <- dplyr::bind_rows(PLOT_PCR, PLOT_MB) %>%
  dplyr::mutate(key = as.character(key)) %>%
  dplyr::select(panel, row_key, key, term, rung, scale, log_est, log_lo,
                log_hi, ratio, ratio_lo, ratio_hi, n, events, summary, I2,
                y_row, y_mark, source)

readr::write_csv(PLOT_FRAME,   file.path(DIR_TABLES, "E45_plot_frame.csv"))
readr::write_csv(CONTRAST_COL, file.path(DIR_TABLES, "E45_contrast_column.csv"))
readr::write_csv(RHO,          file.path(DIR_TABLES, "E45_exposure_rho.csv"))
readr::write_csv(ASSERTIONS,   file.path(DIR_TABLES, "E45_assertions.csv"))

saveRDS(list(
  plot_frame      = PLOT_FRAME,
  contrast_column = CONTRAST_COL,
  rows            = list(pcr = ROWS_PCR, metabric = ROWS_MB),
  exposure_rho    = RHO,
  bracket_check   = BRACKET,
  assertions      = ASSERTIONS,
  legend_draft    = LEGEND,
  # FILE NAMES, NOT PATHS: a path differs between a redirected dry run and the
  # author's run, and would make this element fail a digest comparison for a
  # reason that has nothing to do with the figure. Figures themselves are not
  # digested: a PDF carries its creation time.
  figure_files    = list(
    outputs = basename(paste0(c(PATH_FIG_PCR, PATH_FIG_MB, PATH_FIG_BOTH),
                              rep(c(".png", ".pdf"), each = 3))),
    tracked_in_docs_figures = basename(paste0(PATH_FIG_DOC, c(".png", ".pdf")))),
  spec = list(
    posture = paste0(
      "EXPLORATORY and DISPLAY ONLY. Nothing is fitted, refitted or pooled; ",
      "every plotted value is read from E40, E43 or E44 and exponentiated."),
    narrative = paste0(
      "The author's Google Doc '11.8.26_MK_Myc_Paper' is the narrative of ",
      "record (2026-10-22). Its Fig. 4i is the pCR forest and its Fig. 4j the ",
      "METABRIC forest; v5 called them 4h and 4i. Files are named by content."),
    authority = paste0(
      "docs/2026-10-04_handoff.md section 2 governs what a panel may assert; ",
      "docs/2026-10-22_fig4_human_panels_proposal.md sections 3 and 7 specify ",
      "these panels."),
    axis = paste0(
      "Ratios on a log axis with 1 marked (author, 2026-10-22). Not shared ",
      "between panels: odds ratio and hazard ratio are different quantities."),
    no_bracket = paste0(
      "The pooled pCR contrast is the meta-analysis of the per-cohort ",
      "contrasts (15.4). A bracket between the pooled markers would span the ",
      "difference of two separately pooled estimates, which has no valid ",
      "standard error. See $bracket_check. The contrast is in a column on ",
      "every row instead."),
    rows = paste0(
      "METABRIC: the ER-positive IHC primary, ER.Expr sensitivity, all ",
      "patients, ER-negative IHC. No GSE25066 DRFS row (declaration 15.4; ",
      "the Doc examines long-term outcome in METABRIC alone). The ",
      "ER-negative ER.Expr row goes to ED Fig. 8d."),
    what_15_3_says = paste0(
      "A significant contrast does not make either coefficient a finding. ",
      "Nothing here may be cited in support of a distant-outcome claim for ",
      "respiration, and no OXPHOS value is printed on either panel."),
    instrument = "GSVA only. mitoPPS exists in none of these four cohorts.",
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  built = Sys.time()), PATH_E45)

message("   ", PATH_E45)
message("   4 tables in ", DIR_TABLES)
message("\nE45 part 1 done. Display only; nothing here is interpreted.")
message("A SIGNIFICANT CONTRAST DOES NOT MAKE EITHER COEFFICIENT A FINDING. ",
        "15.3.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E45)

  # every check, with how many values it compared and by how much it passed
  x$assertions %>% as.data.frame()

  # what is drawn, row by row, and where each number came from
  x$plot_frame %>% dplyr::transmute(panel, row_key, key, summary,
      ratio = round(ratio, 3), lo = round(ratio_lo, 3), hi = round(ratio_hi, 3),
      source) %>% as.data.frame()
  x$contrast_column %>% dplyr::select(panel, row_key, reads_as, label) %>%
    as.data.frame()

  # why there is no bracket on the pooled row
  x$bracket_check %>% as.data.frame()

  # OXPHOS against proliferation in all four cohorts
  x$exposure_rho %>% as.data.frame()

  # the legend drafts, descriptive only
  cat(x$legend_draft$pcr, "\n\n", x$legend_draft$metabric, "\n")

  # look at the panels
  FIG_BOTH
}
