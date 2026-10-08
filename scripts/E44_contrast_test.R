# E44_contrast_test.R
# =============================================================================
# THE RESPIRATION-PROLIFERATION CONTRAST, TESTED. IT MAKES NO DECISIONS.
#
# EXPLORATORY, POST-HOC and DESCRIPTIVE. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md,
# section 15, committed at af22110. Read 15.1 to 15.7 first.
#
# =============================================================================
# WHAT THIS COMPUTES (15.1)
# =============================================================================
# Section 4's draft states that the respiration and proliferation coefficients
# have non-overlapping intervals. THAT IS AN EYEBALL COMPARISON OF TWO
# INTERVALS, NOT A TEST OF THE DIFFERENCE BETWEEN THEM. This replaces it with
# the within-model linear contrast
#
#     delta = beta(PROLIF) - beta(OX)
#
# taken from the SAME m1 fit, so
#
#     Var(delta) = Var(bP) + Var(bO) - 2 * Cov(bP, bO)
#
# with the covariance read from that model's variance matrix. THE NAIVE SUM OF
# SQUARED STANDARD ERRORS IS WRONG: OX and PROLIF correlate at Spearman 0.36 to
# 0.49, so their coefficient estimates are correlated and the covariance term
# is not zero. Both are reported, so the size AND SIGN of the correction are
# visible rather than asserted (15.5).
#
# NO REFIT. This is computable only because declaration 14.7 required E43 to
# save the fitted objects. Every number comes from
# e43_proliferation_coefficient.rds$fits.
#
# =============================================================================
# WHAT A SIGNIFICANT CONTRAST DOES NOT DO - 15.3, VERBATIM
# =============================================================================
# "It does not make either coefficient a finding. The contrast can exclude zero
# while the OX coefficient alone remains exactly as weak as E39/E40/E41 section
# 7.2 records it - one boundary cell of thirty, unreproduced by its own
# sensitivity, with the adjacent rungs including the null.
#
# In particular, nothing in E44 may be cited in support of a distant-outcome
# claim for respiration. That prohibition is unchanged and is not relaxed by a
# contrast against a different term. What E44 can license is a statement that
# the two coefficients DIFFER, not that either is established."
#
# =============================================================================
# THE EYEBALL CLAIM IS NOT UNIVERSALLY TRUE (15.7)
# =============================================================================
# Measured before this script was written: the two intervals are disjoint in 4
# of the 8 saved fits - the three pCR cohorts and METABRIC's ER-positive IHC
# stratum, the last by a gap of 0.072. THEY OVERLAP in the pooled METABRIC
# ladder (by 0.031) and in the declared ER.Expr sensitivity (by 0.006), which
# are the two comparators 15.4 requires beside the primary. So section 4 may
# not write "again disjoint" without naming the stratum.
#
# AND OVERLAP IS THE STRICTER CRITERION. Two intervals can overlap while their
# difference still excludes zero, so this test may license a difference where
# the eyeball comparison cannot. The direction of the correction is NOT known
# in advance. This script recomputes the overlap alongside the contrast so the
# two criteria can be read against each other.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It does NOT refit anything. No model is estimated here.
#   - It does NOT compute the contrast as a difference of pooled estimates. The
#     pooled delta is the meta-analysis OF THE PER-COHORT DELTAS (15.4).
#   - It does NOT use the naive SE as a test statistic. Reported for comparison.
#   - It does NOT combine, pool or compare log-odds with log-hazard deltas, and
#     keeps the two endpoints in separate tables with their scales named (15.4).
#   - It does NOT touch GSE25066 DRFS (15.4), whose m1 is not among the saved
#     fits.
#   - It does NOT write or imply that a significant contrast strengthens the OX
#     coefficient on distant outcome (15.3).
#   - It does NOT classify the outcome against 15.6's table, which is printed
#     UNFILLED.
#   - It does NOT interpret. It prints numbers, never verdicts.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  # survival is attached because stats::vcov has NO METHOD for a coxph object
  # unless it is. Without this the extraction fails with "no applicable method
  # for 'vcov'" on five of the eight fits.
  library(survival)
})

message("\nE44: the respiration-proliferation contrast, tested\n", strrep("=", 78))
message("NO REFIT. Every number comes from E43's saved fits (14.7).")
message("A SIGNIFICANT CONTRAST DOES NOT MAKE EITHER COEFFICIENT A FINDING (15.3).")

PATH_E43 <- file.path(DIR_RESULTS, "e43_proliferation_coefficient.rds")
PATH_E44 <- file.path(DIR_RESULTS, "e44_contrast_test.rds")

# =============================================================================
# 0. CONSTANTS - from the declaration, none invented here
# =============================================================================
CI_LEVEL <- 0.95
CRIT     <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

TERM_A <- "PROLIF"      # the + limb of the contrast
TERM_B <- "OX"          # the - limb
N_FITS <- 8L

PCR_COHORTS <- c("GSE25066", "GSE194040", "GSE164458")
MB_PRIMARY  <- "METABRIC / ER-pos (IHC)"
MB_SETS     <- c("METABRIC / pooled", MB_PRIMARY,
                 "METABRIC / ER-pos (ER.Expr)",
                 "METABRIC / ER-neg (IHC)", "METABRIC / ER-neg (ER.Expr)")

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")

# =============================================================================
# 1. Inputs - E43's fitted objects, and NOTHING is refitted
# =============================================================================
message("\n1. inputs (E43's saved fits; NOTHING is refitted)")

.stop_if(file.exists(PATH_E43), "absent: ", PATH_E43,
         ". Source scripts/E43_retrieve_proliferation_coefficient.R first.")
e43 <- readRDS(PATH_E43)
.stop_if(!is.null(e43$fits),
         PATH_E43, " carries no $fits. Declaration 14.7 requires them and E44 ",
         "cannot be computed without them, because the CONTRAST NEEDS THE ",
         "COVARIANCE and only the fitted object carries the variance matrix. ",
         "Re-source E43.")
FITS <- e43$fits
.stop_if(length(FITS) == N_FITS, "$fits holds ", length(FITS), " models, not ",
         N_FITS)

# Every fit must carry BOTH terms and a usable variance matrix, asserted before
# anything is computed (prompt step A's precondition).
PRE <- dplyr::bind_rows(lapply(names(FITS), function(k) {
  f  <- FITS[[k]]
  cf <- tryCatch(stats::coef(f), error = function(e) NULL)
  vc <- tryCatch(stats::vcov(f), error = function(e) NULL)
  tibble::tibble(set = k, class = class(f)[1],
                 has_coef = !is.null(cf), has_vcov = !is.null(vc),
                 has_A = !is.null(cf) && TERM_A %in% names(cf),
                 has_B = !is.null(cf) && TERM_B %in% names(cf),
                 vcov_has_both = !is.null(vc) &&
                   all(c(TERM_A, TERM_B) %in% rownames(vc)))
}))
.pr(PRE)
.stop_if(all(PRE$has_coef & PRE$has_vcov & PRE$has_A & PRE$has_B &
             PRE$vcov_has_both),
         "a fit is missing a coefficient, a variance matrix or one of the two ",
         "terms: ", paste(PRE$set[!(PRE$has_coef & PRE$has_vcov & PRE$has_A &
                                    PRE$has_B & PRE$vcov_has_both)],
                          collapse = ", "), ". STOP.")
message("   ASSERTED: all ", N_FITS, " m1 fits carry ", TERM_A, ", ", TERM_B,
        " and a usable variance matrix")
message("   ", .sub(paste0(
  "GSE25066's DRFS m1 is NOT among these fits and is EXCLUDED by 15.4, ",
  "declared rather than left silent: nothing in that cohort met the reading ",
  "rule at 103 events and a contrast there would be uninformative. If it is ",
  "ever computed it is reported whatever it shows.")))

# =============================================================================
# 2. The contrast helper - a GENERAL linear contrast, not hand-arithmetic
# =============================================================================
# L %*% beta with Var = L %*% V %*% t(L). For L = c(1, -1) against the 2x2
# submatrix this gives Var(bP) + Var(bO) - 2*Cov(bP, bO) with the sign of the
# covariance handled by the algebra rather than by a typed minus sign.
.linear_contrast <- function(fit, terms, L, label) {
  b  <- stats::coef(fit)[terms]
  V  <- stats::vcov(fit)[terms, terms, drop = FALSE]
  est <- as.numeric(L %*% b)
  var <- as.numeric(L %*% V %*% t(L))
  .stop_if(var > 0, label, ": the contrast variance is ", var,
           ", which is not positive. The variance matrix is not usable. STOP.")
  se  <- sqrt(var)
  # the naive SE, for comparison ONLY (15.5). Never the test statistic.
  se_naive <- sqrt(sum(diag(V)))
  tibble::tibble(
    set = label,
    beta_A = unname(b[[terms[1]]]), beta_B = unname(b[[terms[2]]]),
    se_A = sqrt(V[terms[1], terms[1]]), se_B = sqrt(V[terms[2], terms[2]]),
    cov_AB = V[terms[1], terms[2]],
    cor_AB = V[terms[1], terms[2]] /
             sqrt(V[terms[1], terms[1]] * V[terms[2], terms[2]]),
    delta = est, se = se, var = var,
    se_naive = se_naive,
    se_ratio = se / se_naive,
    naive_understates = se > se_naive,
    ci_lo = est - CRIT * se, ci_hi = est + CRIT * se,
    excl_0 = (est - CRIT * se) > 0 | (est + CRIT * se) < 0,
    z = est / se, p = 2 * stats::pnorm(-abs(est / se)),
    # the eyeball criterion recomputed, so the two can be read against each
    # other (15.7). Overlap is the STRICTER criterion.
    a_lo = unname(b[[terms[1]]]) - CRIT * sqrt(V[terms[1], terms[1]]),
    a_hi = unname(b[[terms[1]]]) + CRIT * sqrt(V[terms[1], terms[1]]),
    b_lo = unname(b[[terms[2]]]) - CRIT * sqrt(V[terms[2], terms[2]]),
    b_hi = unname(b[[terms[2]]]) + CRIT * sqrt(V[terms[2], terms[2]])) %>%
    dplyr::mutate(intervals_disjoint = b_hi < a_lo | a_hi < b_lo,
                  interval_gap = pmax(a_lo - b_hi, b_lo - a_hi))
}

CONTRASTS <- dplyr::bind_rows(lapply(names(FITS), function(k)
  .linear_contrast(FITS[[k]], c(TERM_A, TERM_B), matrix(c(1, -1), nrow = 1), k)))
CONTRASTS <- CONTRASTS %>%
  dplyr::left_join(e43$model_info[c("set", "n", "events")], by = "set") %>%
  dplyr::mutate(scale = ifelse(set %in% PCR_COHORTS, "log odds ratio",
                               "log hazard ratio"))

# `scale` is assigned from the set name, which is a hard-coded assumption about
# which cohorts are logistic. ASSERTED against the fitted objects' CLASSES,
# which are authoritative, so the labelling cannot drift from the models.
# (An earlier version read a `family` column from e43$model_info, which has no
# such column; the assignment was a no-op that only emitted a warning.)
SCALE_CHK <- tibble::tibble(
  set = CONTRASTS$set, scale = CONTRASTS$scale,
  fit_class = vapply(CONTRASTS$set, function(k) class(FITS[[k]])[1],
                     character(1))) %>%
  dplyr::mutate(agrees = (fit_class == "glm"   & scale == "log odds ratio") |
                         (fit_class == "coxph" & scale == "log hazard ratio"))
.stop_if(all(SCALE_CHK$agrees),
         "the scale label disagrees with the fitted model class in: ",
         paste(SCALE_CHK$set[!SCALE_CHK$agrees], collapse = ", "),
         ". A log-odds delta would be reported as a log-hazard one or the ",
         "reverse. STOP.")

# =============================================================================
# 3. B. VERIFY against E43 before reporting any contrast - HARD STOP
# =============================================================================
# Every component coefficient extracted here must equal E43's stored value
# EXACTLY. These are the same objects, so anything other than exact equality
# means the extraction is reading the wrong rows.
message("\n3. B. verifying the components against E43's stored coefficients")

STORED <- e43$coefficient_matrix %>%
  dplyr::filter(term %in% c(TERM_A, TERM_B)) %>%
  dplyr::select(set, term, estimate, se) %>%
  tidyr::pivot_wider(names_from = term, values_from = c(estimate, se))
VER <- CONTRASTS %>%
  dplyr::select(set, beta_A, beta_B, se_A, se_B) %>%
  dplyr::left_join(STORED, by = "set") %>%
  dplyr::mutate(
    d_bA = abs(beta_A - .data[[paste0("estimate_", TERM_A)]]),
    d_bB = abs(beta_B - .data[[paste0("estimate_", TERM_B)]]),
    d_sA = abs(se_A  - .data[[paste0("se_", TERM_A)]]),
    d_sB = abs(se_B  - .data[[paste0("se_", TERM_B)]]),
    exact = d_bA == 0 & d_bB == 0 & d_sA == 0 & d_sB == 0)
.pr(VER %>% dplyr::transmute(set, d_beta_PROLIF = d_bA, d_beta_OX = d_bB,
                             d_se_PROLIF = d_sA, d_se_OX = d_sB, exact))
.stop_if(sum(is.na(VER$d_bA)) == 0L,
         "a fit has no stored counterpart in E43's coefficient matrix: ",
         paste(VER$set[is.na(VER$d_bA)], collapse = ", "))
.stop_if(all(VER$exact),
         "COMPONENT MISMATCH in ", sum(!VER$exact), " of ", nrow(VER),
         " fits: ", paste(VER$set[!VER$exact], collapse = ", "),
         ". The extraction is reading different numbers from the same objects, ",
         "so no contrast is reported. STOP.")
message("   ALL ", nrow(VER), " fits: every component coefficient and SE is ",
        "EXACTLY E43's stored value. The contrast is built on the reported ",
        "coefficients.")

# =============================================================================
# 4. A. THE CONTRAST, per fit - and what the covariance correction does
# =============================================================================
message("\n4. A. the contrast delta = beta(", TERM_A, ") - beta(", TERM_B, ")")

.show <- function(tb) .pr(tb %>% dplyr::transmute(set, n, events,
    OX = round(beta_B, 4),
    OX_ci = paste0("[", sprintf("%.4f", b_lo), ", ", sprintf("%.4f", b_hi), "]"),
    PROLIF = round(beta_A, 4),
    PR_ci = paste0("[", sprintf("%.4f", a_lo), ", ", sprintf("%.4f", a_hi), "]"),
    delta = round(delta, 4), se = round(se, 4),
    ci = paste0("[", sprintf("%.4f", ci_lo), ", ", sprintf("%.4f", ci_hi), "]"),
    excl_0, z = round(z, 3), p = signif(p, 3)))

message("\n   ENDPOINT 1: pCR, three cohorts. Scale: LOG ODDS RATIO.")
PCR <- CONTRASTS %>% dplyr::filter(set %in% PCR_COHORTS) %>%
  dplyr::arrange(match(set, PCR_COHORTS))
.show(PCR)

message("\n   ENDPOINT 2: METABRIC breast-cancer-specific survival. Scale: ",
        "LOG HAZARD RATIO.")
message("   ", .sub(paste0(
  "15.4 makes ", MB_PRIMARY, " the one computed directly, with the pooled ",
  "ladder and the ER.Expr sensitivity beside it. The ER-negative strata are ",
  "reported and are uninformative in BOTH terms.")))
MB <- CONTRASTS %>% dplyr::filter(set %in% MB_SETS) %>%
  dplyr::arrange(match(set, MB_SETS))
.show(MB)

message("\n   ", .sub(paste0(
  "THE TWO ENDPOINTS ARE ON DIFFERENT SCALES and are NOT combined, pooled or ",
  "compared as values (15.4, 4.2). Exponentiated, a pCR delta is a RATIO OF ",
  "ODDS RATIOS and a METABRIC delta a RATIO OF HAZARD RATIOS; they are never ",
  "set side by side.")))
.pr(CONTRASTS %>% dplyr::transmute(set, scale,
    exp_delta = round(exp(delta), 4),
    exp_ci = paste0("[", sprintf("%.4f", exp(ci_lo)), ", ",
                    sprintf("%.4f", exp(ci_hi)), "]"),
    reads_as = ifelse(scale == "log odds ratio", "ratio of odds ratios",
                      "ratio of hazard ratios")))

# --- 15.5: the covariance correction, made visible ------------------------
message("\n   4.1 WHY THE NAIVE CALCULATION IS WRONG, per fit (15.5)")
.pr(CONTRASTS %>% dplyr::transmute(set,
    cov_AB = signif(cov_AB, 3), cor_AB = round(cor_AB, 3),
    se_correct = round(se, 4), se_naive = round(se_naive, 4),
    ratio = round(se_ratio, 4),
    naive_too_small = naive_understates))
n_under <- sum(CONTRASTS$naive_understates)
message("   ", .sub(paste0(
  "The naive SE is sqrt(seP^2 + seO^2), which assumes Cov(bP, bO) = 0. The ",
  "correct SE is larger in ", n_under, " of ", nrow(CONTRASTS),
  " fits and smaller in ", nrow(CONTRASTS) - n_under,
  ". Where it is larger THE NAIVE CALCULATION IS ANTI-CONSERVATIVE, not ",
  "merely wrong: it would have produced too narrow an interval and too small ",
  "a p value. 15.7 recorded the expectation that this would be the common ",
  "case, because positively correlated predictors usually give negatively ",
  "correlated coefficient estimates. The numbers above are the measured ",
  "version; the naive SE is reported for comparison ONLY and is never the ",
  "test statistic.")))

# --- 15.7: the eyeball criterion, recomputed beside the test --------------
message("\n   4.2 THE EYEBALL CRITERION BESIDE THE TEST (15.7)")
.pr(CONTRASTS %>% dplyr::transmute(set,
    intervals_disjoint, interval_gap = round(interval_gap, 4),
    delta_excl_0 = excl_0,
    criteria_agree = intervals_disjoint == excl_0))
message("   ", .sub(paste0(
  "NON-OVERLAP IS THE STRICTER CRITERION: two intervals can overlap while ",
  "their difference still excludes zero, but not the reverse. So a row where ",
  "intervals_disjoint is FALSE and delta_excl_0 is TRUE is the expected kind ",
  "of disagreement, and a row the other way round would need explaining. ",
  "Section 4's draft rests on the eyeball column; 15.1 replaces it with the ",
  "delta column.")))

# =============================================================================
# 5. C. META-ANALYSE the three pCR deltas - of the DELTAS, not of the pooled
# =============================================================================
# 15.4: "The difference of two separately pooled estimates is NOT the pooled
# difference and must not be computed that way." The pooled delta below is the
# DerSimonian-Laird meta-analysis OF THE PER-COHORT DELTAS, each of which
# already carries its own covariance correction.
message("\n5. C. meta-analysis of the three pCR deltas (15.4)")

.meta <- function(est, se, label) {
  v <- se^2; w <- 1 / v; k <- length(est)
  if (k < 2L) {
    return(tibble::tibble(quantity = label, k = k,
      fe_estimate = est, fe_ci_lo = est - CRIT*se, fe_ci_hi = est + CRIT*se,
      fe_p = 2*stats::pnorm(-abs(est/se)),
      re_estimate = est, re_ci_lo = est - CRIT*se, re_ci_hi = est + CRIT*se,
      re_p = 2*stats::pnorm(-abs(est/se)),
      Q = NA_real_, Q_p = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  }
  fe <- sum(w*est)/sum(w); vfe <- 1/sum(w)
  Q  <- sum(w*(est-fe)^2); df <- k - 1
  Cc <- sum(w) - sum(w^2)/sum(w)
  tau2 <- max(0, (Q - df)/Cc)
  ws <- 1/(v + tau2); re <- sum(ws*est)/sum(ws); vre <- 1/sum(ws)
  tibble::tibble(quantity = label, k = k,
    fe_estimate = fe, fe_ci_lo = fe - CRIT*sqrt(vfe),
    fe_ci_hi = fe + CRIT*sqrt(vfe), fe_p = 2*stats::pnorm(-abs(fe/sqrt(vfe))),
    re_estimate = re, re_ci_lo = re - CRIT*sqrt(vre),
    re_ci_hi = re + CRIT*sqrt(vre), re_p = 2*stats::pnorm(-abs(re/sqrt(vre))),
    Q = Q, Q_p = stats::pchisq(Q, df, lower.tail = FALSE),
    tau2 = tau2, I2 = 100*max(0, (Q - df)/Q))
}

META <- .meta(PCR$delta, PCR$se, "delta = PROLIF - OX, pCR")
META$cohorts     <- paste(PCR$set, collapse = ", ")
META$signs       <- paste(sprintf("%s %s", PCR$set,
                                  ifelse(PCR$delta > 0, "+", "-")),
                          collapse = " ")
META$signs_agree <- length(unique(sign(PCR$delta))) == 1L
META$summary_is  <- ifelse(!is.na(META$I2) & META$I2 >= 50, "random effects",
                           "fixed effect")
.pr(META %>% dplyr::transmute(quantity, k,
    fe = round(fe_estimate, 4),
    fe_ci = paste0("[", sprintf("%.4f", fe_ci_lo), ", ",
                   sprintf("%.4f", fe_ci_hi), "]"),
    fe_p = signif(fe_p, 3),
    re = round(re_estimate, 4),
    re_ci = paste0("[", sprintf("%.4f", re_ci_lo), ", ",
                   sprintf("%.4f", re_ci_hi), "]"),
    re_p = signif(re_p, 3),
    Q = round(Q, 3), Q_p = signif(Q_p, 3), tau2 = round(tau2, 4),
    I2 = round(I2, 1), summary_is, signs_agree, signs))

# The forbidden calculation, shown as forbidden so nobody reaches for it.
WRONG <- e43$meta %>% dplyr::select(term, fe_estimate, re_estimate)
diff_of_pooled_fe <- WRONG$fe_estimate[WRONG$term == TERM_A] -
                     WRONG$fe_estimate[WRONG$term == TERM_B]
message("   ", .sub(paste0(
  "15.4 FORBIDS computing this as a difference of pooled estimates. For the ",
  "record, that calculation would give ", round(diff_of_pooled_fe, 4),
  " against the correct pooled delta of ", round(META$fe_estimate, 4),
  " (fixed) - and it has NO VALID STANDARD ERROR at all, because the two ",
  "pooled estimates come from overlapping patients and their covariance is ",
  "not available from the pooled summaries. IT IS NOT REPORTED AS AN ",
  "ESTIMATE and appears here only so that nobody reaches for it.")))
message("   ", .sub(paste0(
  "BOTH summaries are reported always; at I2 >= 50 the random-effects ",
  "estimate is the one read, with the fixed effect beside it (section 9).")))

# =============================================================================
# 6. 15.6's reading rule, printed UNFILLED. This script does NOT classify.
# =============================================================================
message("\n6. declaration 15.6, printed UNFILLED")

READING_RULE <- tibble::tribble(
  ~outcome, ~`what the Results may say`,
  "delta excludes zero on both endpoints, signs agreeing across the pCR cohorts",
    "the two coefficients differ, tested, on both endpoints",
  "delta excludes zero on pCR only",
    "the opposition is tested on treatment response and remains a described direction on distant outcome",
  "delta includes zero on an endpoint",
    "the non-overlap statement is withdrawn for that endpoint; the two intervals are reported and nothing is claimed about their difference")
.pr(READING_RULE)
message("   ", .sub(paste0(
  "This script does not classify the outcome against this table. No p value ",
  "is described as a near miss (section 9).")))

message("\n   AND 15.3, WHICH GOVERNS EVERY ROW OF IT:")
message("   ", .sub(paste0(
  "A SIGNIFICANT CONTRAST DOES NOT MAKE EITHER COEFFICIENT A FINDING. The ",
  "contrast can exclude zero while the OX coefficient alone remains exactly ",
  "as weak as E39/E40/E41 section 7.2 records it - one boundary cell of ",
  "thirty, unreproduced by its own sensitivity, with the adjacent rungs ",
  "including the null. NOTHING IN E44 MAY BE CITED IN SUPPORT OF A ",
  "DISTANT-OUTCOME CLAIM FOR RESPIRATION. What E44 can license is a statement ",
  "that the two coefficients DIFFER, not that either is established.")))

# =============================================================================
# 7. D. Save - the full objects, per 14.7
# =============================================================================
message("\n7. D. saving")

readr::write_csv(CONTRASTS,    file.path(DIR_TABLES, "E44_contrasts_all.csv"))
readr::write_csv(PCR,          file.path(DIR_TABLES, "E44_contrast_pcr.csv"))
readr::write_csv(MB,           file.path(DIR_TABLES, "E44_contrast_metabric.csv"))
readr::write_csv(META,         file.path(DIR_TABLES, "E44_contrast_pcr_meta.csv"))
readr::write_csv(VER,          file.path(DIR_TABLES, "E44_component_check.csv"))
readr::write_csv(PRE,          file.path(DIR_TABLES, "E44_fit_preconditions.csv"))
readr::write_csv(SCALE_CHK,    file.path(DIR_TABLES, "E44_scale_check.csv"))

saveRDS(list(
  contrasts        = CONTRASTS,
  pcr              = PCR,
  metabric         = MB,
  meta             = META,
  component_check  = VER,
  preconditions    = PRE,
  scale_check      = SCALE_CHK,
  reading_rule     = READING_RULE,
  diff_of_pooled_forbidden = list(
    value = diff_of_pooled_fe,
    why = paste0(
      "15.4 forbids it. Shown only so nobody reaches for it: it has NO VALID ",
      "STANDARD ERROR, because the two pooled estimates come from overlapping ",
      "patients and their covariance is not recoverable from pooled ",
      "summaries. The correct pooled delta is the meta-analysis of the ",
      "per-cohort deltas, each already covariance-corrected.")),
  spec = list(
    declaration = paste0("docs/2026-10-02_E39_respiratory_axis_decomposition",
                         "_declaration.md section 15, committed at af22110"),
    posture = paste0(
      "EXPLORATORY, POST-HOC and DESCRIPTIVE (15.2). A derived quantity from ",
      "models already declared and already reported, computed to replace a ",
      "descriptive statement the Results would otherwise make loosely. It is ",
      "not a hypothesis and carries no falsification criterion."),
    # 15.3, VERBATIM, as the prompt requires.
    what_a_significant_contrast_does_not_do = paste0(
      "It does not make either coefficient a finding. The contrast can ",
      "exclude zero while the OX coefficient alone remains exactly as weak as ",
      "E39/E40/E41 section 7.2 records it - one boundary cell of thirty, ",
      "unreproduced by its own sensitivity, with the adjacent rungs including ",
      "the null. In particular, nothing in E44 may be cited in support of a ",
      "distant-outcome claim for respiration. That prohibition is unchanged ",
      "and is not relaxed by a contrast against a different term. What E44 ",
      "can license is a statement that the two coefficients DIFFER, not that ",
      "either is established."),
    estimand = paste0(
      "delta = beta(PROLIF) - beta(OX) from the SAME m1 fit, with ",
      "Var(delta) = Var(bP) + Var(bO) - 2*Cov(bP, bO) and the covariance read ",
      "from that model's variance matrix. Computed as a general linear ",
      "contrast L %*% beta with Var = L V L', L = c(1, -1), so the sign of ",
      "the covariance is handled by the algebra rather than by hand."),
    naive_se = paste0(
      "15.5. The naive SE, sqrt(seP^2 + seO^2), assumes Cov(bP, bO) = 0 and ",
      "is reported BESIDE the correct one for comparison ONLY. It is never ",
      "the test statistic. Where the correct SE is LARGER the naive ",
      "calculation is ANTI-CONSERVATIVE, not merely wrong."),
    pooling = paste0(
      "15.4. The pooled pCR delta is the DerSimonian-Laird meta-analysis OF ",
      "THE PER-COHORT DELTAS. The difference of two separately pooled ",
      "estimates is NOT the pooled difference, has no valid standard error, ",
      "and must not be computed that way."),
    scales = paste0(
      "15.4 and 4.2. pCR deltas are log odds ratios and METABRIC deltas log ",
      "hazard ratios. They are NOT combined, pooled or compared as values, ",
      "and are kept in separate tables. Exponentiated, a pCR delta is a ratio ",
      "of odds ratios and a METABRIC delta a ratio of hazard ratios; they are ",
      "never set side by side."),
    eyeball = paste0(
      "15.7. The two intervals are disjoint in only 4 of the 8 saved fits - ",
      "the three pCR cohorts and METABRIC's ER-positive IHC stratum. They ",
      "OVERLAP in the pooled ladder and in the declared ER.Expr sensitivity, ",
      "so section 4 may not write 'again disjoint' without naming the ",
      "stratum. NON-OVERLAP IS THE STRICTER CRITERION, so this test may ",
      "license a difference where the eyeball comparison cannot; the ",
      "direction of the correction was not known in advance. Both criteria ",
      "are reported side by side in $contrasts."),
    excluded = paste0(
      "15.4. GSE25066 DRFS is excluded, declared rather than left silent: its ",
      "m1 is not among E43's saved fits, nothing in that cohort met the ",
      "reading rule at 103 events, and a contrast there would be ",
      "uninformative. If it is ever computed it is reported whatever it ",
      "shows."),
    no_refit = paste0(
      "NOTHING was refitted. Every number comes from ",
      "e43_proliferation_coefficient.rds$fits, and every component ",
      "coefficient and SE was verified EXACTLY equal to E43's stored value ",
      "before any contrast was reported."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  built = Sys.time()), PATH_E44)

message("   ", PATH_E44)
message("   7 tables in ", DIR_TABLES)
message("\nE44 done. Numbers only; nothing here is interpreted or classified.")
message("A SIGNIFICANT CONTRAST DOES NOT MAKE EITHER COEFFICIENT A FINDING. ",
        "15.3.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E44)

  # READ THIS FIRST, every time.
  cat(x$spec$what_a_significant_contrast_does_not_do, "\n")

  # B. the component check. Exact equality or the script would have stopped.
  x$component_check %>% as.data.frame()

  # A. the contrast, by endpoint. The scales are different; keep them apart.
  x$pcr %>% dplyr::transmute(set, n, events, OX = round(beta_B, 4),
      PROLIF = round(beta_A, 4), delta = round(delta, 4),
      ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4), excl_0,
      p = signif(p, 3)) %>% as.data.frame()
  x$metabric %>% dplyr::transmute(set, n, events, OX = round(beta_B, 4),
      PROLIF = round(beta_A, 4), delta = round(delta, 4),
      ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4), excl_0,
      p = signif(p, 3)) %>% as.data.frame()

  # 4.1 what the covariance correction did. Watch naive_too_small.
  x$contrasts %>% dplyr::transmute(set, cor_AB = round(cor_AB, 3),
      se = round(se, 4), se_naive = round(se_naive, 4),
      ratio = round(se_ratio, 4), naive_understates) %>% as.data.frame()

  # 4.2 the eyeball criterion beside the test. Disagreement is expected one way.
  x$contrasts %>% dplyr::transmute(set, intervals_disjoint,
      interval_gap = round(interval_gap, 4), delta_excl_0 = excl_0,
      criteria_agree = intervals_disjoint == excl_0) %>% as.data.frame()

  # C. pooled across the three pCR cohorts - of the DELTAS, not of the pooled.
  x$meta %>% as.data.frame()
  cat(x$diff_of_pooled_forbidden$why, "\n")

  # 15.6, UNFILLED. The script does not classify.
  x$reading_rule %>% as.data.frame()

  utils::str(x$spec)
}
