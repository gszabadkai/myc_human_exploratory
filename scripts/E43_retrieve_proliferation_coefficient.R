# E43_retrieve_proliferation_coefficient.R
# =============================================================================
# A RETRIEVAL, NOT A NEW ANALYSIS. IT MAKES NO DECISIONS.
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md,
# section 14, committed at 9fb7852. Read 14.1 to 14.7 first.
#
# =============================================================================
# WHAT WAS DISCARDED, AND WHY THIS IS A RETRIEVAL (14.1, 14.2)
# =============================================================================
# E40 and E41 carried PROLIF_DISJOINT in every rung from m1 onward and
# extracted only the OX row of each coefficient matrix, through a helper taking
# co["OX", ]. The proliferation coefficient, its SE, interval, p and VIF were
# COMPUTED AND DISCARDED AT EXTRACTION. The `terms` column records the formula
# as a string, so the term is documented as having been fitted while its
# estimate exists nowhere on disk. Neither object carries the fitted models.
#
# THE MODEL IS UNCHANGED, THE TERM WAS DECLARED, AND ONLY THE EXTRACTION IS
# NEW. 14.2 grounds that in the Block C precedent already established in this
# arm (section 2.1): "Those main effects are relevant to a constraint reading,
# and retrieving them is a one-line change ... when retrieved they inform E2's
# interpretation without redefining it."
#
# BUT THE READING IT SUPPORTS IS POST-HOC AND MUST BE WRITTEN AS SUCH. Section
# 9.1's sign-pair table anticipated better/worse, no-effect/worse and
# better/better. The observed worse-pCR-with-better-survival combination is NOT
# IN IT. That gap is recorded in 14.2 rather than repaired retrospectively.
#
# =============================================================================
# WHY m1 IS THE INTERPRETABLE RUNG (14.4)
# =============================================================================
# MYC is ABSENT from m1, so rho(MYC, PROLIF_DISJOINT) of 0.78 / 0.79 / 0.82 in
# the pCR cohorts and 0.71 in METABRIC DOES NOT BEAR on this coefficient. It is
# governed by rho(OX, PROLIF) of 0.36 to 0.49, inside the declaration's own
# separable band. m0, m2 and m3 are NOT retrieved.
#
# =============================================================================
# EXPECTED ATTENUATION, STATED BEFORE THE FIT (14.5)
# =============================================================================
# Subtype is in the spine, and subtype is substantially a proliferation
# classification - PAM50 LumA against LumB in METABRIC, hormone-receptor status
# in the neoadjuvant cohorts. Much of proliferation's prognostic and predictive
# signal is therefore absorbed before the PROLIF term is reached.
#
#   AN ATTENUATED OR NULL PROLIF COEFFICIENT IS EXPECTED AND DOES NOT REFUTE
#   THE PUBLISHED PATTERN. It would say that proliferation adds little within
#   subtype, not that the relationship is absent.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It retrieves m1 ONLY. NO m0, NO m2, NO m3 (14.4).
#   - It does NOT rescore anything and does NOT reload METABRIC_DATA.RData.
#     Every fit comes from a saved frame.
#   - It does NOT classify the outcome against 14.6's table. That table is
#     printed UNFILLED.
#   - It does NOT report one endpoint without the other (14.6).
#   - It does NOT interpret. It prints numbers, never verdicts.
#   - It writes NOTHING to myc_human_validation, FROZEN at d3ac60e, or to
#     myc_mouse. h4_outcome_models.rds is read READ-ONLY.
#
# PER 14.7 THIS SCRIPT SAVES THE FULL COEFFICIENT MATRIX, THE FULL VIF VECTOR
# AND THE FITTED MODEL OBJECTS. That is the standing fix, and this script is
# the first to apply it.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

message("\nE43: retrieving the proliferation coefficient\n", strrep("=", 78))
message("A RETRIEVAL, not a new analysis. The model is unchanged (14.2).")
message("m1 ONLY. The reproduction check is a HARD STOP (14.3).")

PATH_E39 <- file.path(DIR_RESULTS, "e39_prepare_and_gate.rds")
PATH_E40 <- file.path(DIR_RESULTS, "e40_ladder.rds")
PATH_E41 <- file.path(DIR_RESULTS, "e41_metabric_ladder.rds")
# READ-ONLY, from the FROZEN repo. Nothing is ever written there.
PATH_H4  <- paste0("/Users/gs/code/myc_human_validation/results/",
                   "h4_outcome_models.rds")
PATH_E43 <- file.path(DIR_RESULTS, "e43_proliferation_coefficient.rds")

# =============================================================================
# 0. CONSTANTS - from the declaration and from the stored objects
# =============================================================================
CI_LEVEL <- 0.95
CRIT     <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

COHORTS <- c("GSE25066", "GSE194040", "GSE164458")
PRIMARY <- "GSE25066"

# 14.3. ONE rung. The RHS is built exactly as E40 and E41 built it - OX, then
# the surviving spine, then PROLIF - because term order does not change a fit
# but using the original construction is what makes a 6-decimal reproduction
# check a real check rather than a coincidence.
RUNG <- "m1"
SPINE_CANDIDATES <- c("subtype_term", "treatment")

# 14.3's stop condition.
REPRO_DP  <- 6L
REPRO_TOL <- 10^(-REPRO_DP)

VIF_CAVEAT <- 5          # section 5.3

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.n_lev   <- function(v) dplyr::n_distinct(v[!is.na(v)])

# =============================================================================
# 1. Inputs - saved frames only. NOTHING is rescored (14.3).
# =============================================================================
message("\n1. inputs (saved frames only; NOTHING is rescored)")

for (p in c(PATH_E39, PATH_E40, PATH_E41, PATH_H4))
  .stop_if(file.exists(p), "absent: ", p)
e39 <- readRDS(PATH_E39)
e40 <- readRDS(PATH_E40)
e41 <- readRDS(PATH_E41)
h4  <- readRDS(PATH_H4)

.stop_if(!is.null(e41$frame),
         PATH_E41, " carries no $frame. It was added to E41's save block on ",
         "2026-10-03 (51465b6) and the object on disk predates that. ",
         "RE-SOURCE scripts/E41_metabric_ladder.R first. E43 will NOT reload ",
         "METABRIC_DATA.RData or rescore.")
.stop_if(!is.null(e39$frame), PATH_E39, " carries no $frame")

# --- the three pCR frames, rebuilt exactly as E40 rebuilt them -------------
PS  <- e39$prolif
D25 <- e39$frame %>%
  dplyr::left_join(PS[PS$cohort == PRIMARY, c("sample_id", "PROLIF")],
                   by = "sample_id")
.stop_if(sum(is.na(D25$PROLIF)) == 0L, "PROLIF did not join to every ",
         PRIMARY, " row")
FR <- stats::setNames(lapply(COHORTS, function(cn) {
  if (cn == PRIMARY) return(D25)
  h4$frames[[cn]] %>%
    dplyr::left_join(PS[PS$cohort == cn, c("sample_id", "PROLIF")],
                     by = "sample_id")
}), COHORTS)
for (cn in COHORTS) .stop_if(sum(is.na(FR[[cn]]$PROLIF)) == 0L,
                             cn, ": PROLIF has NA after the join")
# E40's subtype_term construction, reused UNCHANGED (14.3).
for (cn in COHORTS) {
  FR[[cn]]$subtype_term <- if (cn == PRIMARY) FR[[cn]]$subtype2 else
    factor(FR[[cn]]$subtype)
}
PCR_SETS <- stats::setNames(lapply(COHORTS, function(cn) {
  d <- FR[[cn]]
  if (cn == PRIMARY) d <- d[d$in_model, ]
  d
}), COHORTS)
message("   pCR frames: ",
        paste(sprintf("%s %d", COHORTS,
                      vapply(PCR_SETS, nrow, integer(1))), collapse = " | "))

# --- the METABRIC frame and its five analysis sets, as E41 built them ------
MB <- e41$frame
MB$subtype_term <- stats::relevel(factor(MB$subtype_pam50), ref = "LumA")
MB_SETS <- list(
  `pooled`             = MB,
  `ER-pos (IHC)`       = MB[!is.na(MB$er_primary)   & MB$er_primary   == "P", ],
  `ER-neg (IHC)`       = MB[!is.na(MB$er_primary)   & MB$er_primary   == "N", ],
  `ER-pos (ER.Expr)`   = MB[MB$er_sensitive == "P", ],
  `ER-neg (ER.Expr)`   = MB[MB$er_sensitive == "N", ])
message("   METABRIC sets: ",
        paste(sprintf("%s %d", names(MB_SETS),
                      vapply(MB_SETS, nrow, integer(1))), collapse = " | "))
message("   ", .sub(paste0(
  "Both subtype definitions are E40's and E41's, reused unchanged (14.3): ",
  "GSE25066's collapsed subtype2, GSE194040's four-level subtype, ",
  "GSE164458's constant TNBC, and METABRIC's PAM50 with LumA as reference ",
  "and NC kept as its own level.")))

# =============================================================================
# 2. Spine audit - constants omitted EXPLICITLY, as E40, E41 and E42 did
# =============================================================================
message("\n2. spine audit (5.2)")

.audit <- function(d, label) {
  tibble::tibble(set = label, n = nrow(d), term = SPINE_CANDIDATES,
                 n_levels = c(.n_lev(d$subtype_term), .n_lev(d$treatment))) %>%
    dplyr::mutate(constant = n_levels <= 1L,
                  action = ifelse(constant, "OMITTED (constant)", "fitted"))
}
AUDIT <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(COHORTS, function(cn)
    .audit(PCR_SETS[[cn]], cn))),
  dplyr::bind_rows(lapply(names(MB_SETS), function(k)
    .audit(MB_SETS[[k]], paste0("METABRIC / ", k)))))
.pr(AUDIT)
.spine_for <- function(label) {
  a <- AUDIT[AUDIT$set == label, ]
  a$term[!a$constant]
}
# The two constants E40 declared (5.2), asserted again rather than assumed.
.stop_if(.n_lev(PCR_SETS[[PRIMARY]]$treatment) == 1L,
         "treatment is NOT constant in ", PRIMARY, "; 5.2 says it is")
.stop_if(.n_lev(PCR_SETS[["GSE164458"]]$subtype_term) == 1L,
         "subtype is NOT constant in GSE164458; 5.2 says it is")
message("   ASSERTED: treatment constant in ", PRIMARY,
        ", subtype constant in GSE164458 - both OMITTED EXPLICITLY")

# Empty factor levels, dropped explicitly and reported, exactly as E41 did.
# An EMPTY level contributes nothing and changes no estimate; it is dropped so
# car::vif does not fail on aliased coefficients.
DROPPED <- list()
.drop_empty <- function(d, label, spine) {
  for (v in spine) {
    if (!is.factor(d[[v]])) next
    gone <- setdiff(levels(d[[v]]), unique(as.character(d[[v]])))
    if (length(gone)) {
      DROPPED[[length(DROPPED) + 1L]] <<- tibble::tibble(
        set = label, term = v, dropped = paste(gone, collapse = ", "),
        n_dropped = length(gone))
      d[[v]] <- droplevels(d[[v]])
    }
  }
  d
}

# =============================================================================
# 3. The estimator helpers - FULL matrix and FULL vector this time (14.3, 14.7)
# =============================================================================
# car::vif's "No intercept: vifs may not be sensible" arises on the COX fits
# only - a Cox model has no intercept by design and the VIF comes from the
# correlation matrix of the coefficient estimates, which needs none. It is
# counted and muffled deliberately, as E41 did. It does NOT arise on the
# logistic fits, which have an intercept.
.VIF_STATE <- new.env(parent = emptyenv())
.VIF_STATE$no_intercept <- 0L
OTHER_WARNINGS <- list()

.catch <- function(expr, label) {
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (grepl("No intercept", m, fixed = TRUE)) {
      .VIF_STATE$no_intercept <- .VIF_STATE$no_intercept + 1L
    } else {
      OTHER_WARNINGS[[length(OTHER_WARNINGS) + 1L]] <<- tibble::tibble(
        where = label, message = m)
    }
    invokeRestart("muffleWarning")
  })
}

# THE FULL coefficient matrix, every term, as a tidy table.
.coef_matrix <- function(fit, label, family) {
  co <- summary(fit)$coefficients
  se_col <- if (inherits(fit, "coxph")) "se(coef)" else 2L
  tibble::tibble(
    set = label, family = family, term = rownames(co),
    estimate = unname(co[, 1]), se = unname(co[, se_col]),
    statistic = unname(co[, ncol(co) - 1L]), p = unname(co[, ncol(co)])) %>%
    dplyr::mutate(ci_lo = estimate - CRIT * se,
                  ci_hi = estimate + CRIT * se,
                  excl_0 = ci_lo > 0 | ci_hi < 0)
}

# THE FULL VIF vector, every term, with the measure named.
.vif_vector <- function(fit, label) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2L) return(tibble::tibble(set = label, term = NA_character_,
                                     vif = NA_real_, vif_df = NA_real_,
                                     measure = "undefined (1 predictor)"))
  v <- .catch(tryCatch(car::vif(fit), error = function(e) NULL),
              paste0("car::vif / ", label))
  if (is.null(v)) return(tibble::tibble(set = label, term = NA_character_,
                                        vif = NA_real_, vif_df = NA_real_,
                                        measure = "car::vif failed"))
  if (is.matrix(v))
    tibble::tibble(set = label, term = rownames(v), vif = unname(v[, 1]),
                   vif_df = unname(v[, 2]),
                   measure = "GVIF (generalised; a 1-df term's GVIF == VIF)")
  else
    tibble::tibble(set = label, term = names(v), vif = unname(v),
                   vif_df = 1, measure = "VIF")
}

FITS <- list()          # SAVED this time, per 14.7

.fit_m1 <- function(d, label, family, time_var = NULL, event_var = NULL) {
  sp <- .spine_for(label)
  d  <- .drop_empty(d, label, sp)
  # E40's and E41's exact RHS order: OX, surviving spine, then PROLIF.
  rhs <- c("OX", sp, "PROLIF")
  resp <- if (family == "logistic") "pcr" else
    sprintf("survival::Surv(%s, %s)", time_var, event_var)
  f <- stats::reformulate(rhs, response = resp)
  fit <- .catch(
    if (family == "logistic")
      stats::glm(f, data = d, family = stats::binomial())
    else survival::coxph(f, data = d),
    paste0(tolower(family), " / ", label))
  FITS[[label]] <<- fit
  list(fit = fit,
       coefs = .coef_matrix(fit, label, family),
       vifs  = .vif_vector(fit, label),
       n = if (inherits(fit, "coxph")) fit$n else stats::nobs(fit),
       events = if (inherits(fit, "coxph")) fit$nevent else sum(fit$y == 1),
       omitted = paste(setdiff(SPINE_CANDIDATES, sp), collapse = ", "),
       terms = paste(attr(stats::terms(fit), "term.labels"), collapse = " + "))
}

# =============================================================================
# 4. A. REFIT m1 - three pCR cohorts and five METABRIC sets
# =============================================================================
message("\n4. A. refitting m1 (nothing rescored)")

RES <- list()
for (cn in COHORTS)
  RES[[cn]] <- .fit_m1(PCR_SETS[[cn]], cn, "logistic")
for (k in names(MB_SETS))
  RES[[paste0("METABRIC / ", k)]] <-
    .fit_m1(MB_SETS[[k]], paste0("METABRIC / ", k), "cox",
            "bcss_time_years", "bcss_event")
message("   fitted ", length(RES), " models: ", length(COHORTS),
        " logistic (pCR) + ", length(MB_SETS),
        " Cox (METABRIC breast-cancer-specific survival)")

COEFS <- dplyr::bind_rows(lapply(RES, `[[`, "coefs"))
VIFS  <- dplyr::bind_rows(lapply(RES, `[[`, "vifs"))
META_INFO <- dplyr::bind_rows(lapply(names(RES), function(k)
  tibble::tibble(set = k, n = RES[[k]]$n, events = RES[[k]]$events,
                 omitted = RES[[k]]$omitted, terms = RES[[k]]$terms)))

# =============================================================================
# 5. B. REPRODUCTION CHECK - A HARD STOP (14.3)
# =============================================================================
# Every refitted OX coefficient AND its SE must reproduce the stored value to
# at least 6 decimal places. A refit that does not reproduce the original is
# NOT the original model, and the PROLIF estimate cannot be read from it.
message("\n5. B. reproduction check against the stored OX coefficients")

STORED <- dplyr::bind_rows(
  e40$pcr %>% dplyr::filter(rung == RUNG) %>%
    dplyr::transmute(set = cohort, stored_estimate = estimate,
                     stored_se = se, stored_n = n, stored_events = events),
  e41$ladder %>% dplyr::filter(rung == RUNG, endpoint == "BCSS") %>%
    dplyr::transmute(set = paste0("METABRIC / ", set),
                     stored_estimate = estimate, stored_se = se,
                     stored_n = n, stored_events = events),
  e41$er_stratified %>% dplyr::filter(rung == RUNG) %>%
    dplyr::transmute(set = paste0("METABRIC / ", set),
                     stored_estimate = estimate, stored_se = se,
                     stored_n = n, stored_events = events))

REPRO <- COEFS %>% dplyr::filter(term == "OX") %>%
  dplyr::select(set, refit_estimate = estimate, refit_se = se) %>%
  dplyr::left_join(STORED, by = "set") %>%
  dplyr::left_join(META_INFO[c("set", "n", "events")], by = "set") %>%
  dplyr::mutate(d_estimate = abs(refit_estimate - stored_estimate),
                d_se       = abs(refit_se - stored_se),
                n_matches  = n == stored_n,
                ev_matches = events == stored_events,
                ok = d_estimate < REPRO_TOL & d_se < REPRO_TOL &
                     n_matches & ev_matches)
.pr(REPRO %>% dplyr::transmute(set, n, events,
    refit = signif(refit_estimate, 10), stored = signif(stored_estimate, 10),
    d_est = signif(d_estimate, 3), d_se = signif(d_se, 3),
    n_ev_match = n_matches & ev_matches, ok))

.stop_if(sum(is.na(REPRO$stored_estimate)) == 0L,
         "a refitted set has no stored counterpart: ",
         paste(REPRO$set[is.na(REPRO$stored_estimate)], collapse = ", "),
         ". The set labels do not line up with E40/E41 and the check is ",
         "not meaningful. STOP.")
.stop_if(nrow(REPRO) == length(RES),
         "the check covers ", nrow(REPRO), " of ", length(RES), " fits")
.stop_if(all(REPRO$ok),
         "REPRODUCTION FAILED in ", sum(!REPRO$ok), " of ", nrow(REPRO),
         " cells: ", paste(REPRO$set[!REPRO$ok], collapse = ", "),
         ". Largest discrepancy on the estimate ",
         signif(max(REPRO$d_estimate, na.rm = TRUE), 3), ", on the SE ",
         signif(max(REPRO$d_se, na.rm = TRUE), 3),
         ". A REFIT THAT DOES NOT REPRODUCE THE ORIGINAL IS NOT THE ORIGINAL ",
         "MODEL, so the PROLIF estimate CANNOT BE READ and is not reported ",
         "(declaration 14.3). STOP.")
message("   ALL ", nrow(REPRO), " cells reproduce to better than ",
        REPRO_DP, " decimal places, on the estimate AND the SE, with n and ",
        "events matching.")
message("   Largest discrepancy: estimate ",
        signif(max(REPRO$d_estimate), 3), ", SE ",
        signif(max(REPRO$d_se), 3), ".")
message("   ", .sub(paste0(
  "So these are the SAME MODELS E40 and E41 fitted, and the PROLIF ",
  "coefficient below is the one they computed and discarded (14.1).")))

# =============================================================================
# 6. C. EXTRACT - the full matrix, the full vector, then the focused table
# =============================================================================
message("\n6. C. the FULL coefficient matrix, every term (14.3, 14.7)")
.pr(COEFS %>% dplyr::transmute(set, family, term,
    estimate = round(estimate, 4), se = round(se, 4),
    ci = paste0("[", sprintf("%.4f", ci_lo), ", ",
                sprintf("%.4f", ci_hi), "]"),
    excl_0, p = signif(p, 3)))

message("\n   the FULL VIF vector, every term:")
.pr(VIFS %>% dplyr::transmute(set, term, vif = round(vif, 3), vif_df,
                              measure))

message("\n   D1. THE PROLIF TERM, with OX beside it in every row (14.3):")
PROLIF_TAB <- COEFS %>% dplyr::filter(term %in% c("OX", "PROLIF")) %>%
  dplyr::select(set, family, term, estimate, se, ci_lo, ci_hi, p, excl_0) %>%
  tidyr::pivot_wider(names_from = term,
                     values_from = c(estimate, se, ci_lo, ci_hi, p, excl_0)) %>%
  dplyr::left_join(META_INFO, by = "set") %>%
  dplyr::left_join(VIFS %>% dplyr::filter(term == "PROLIF") %>%
                     dplyr::select(set, vif_PROLIF = vif), by = "set") %>%
  dplyr::left_join(VIFS %>% dplyr::filter(term == "OX") %>%
                     dplyr::select(set, vif_OX = vif), by = "set") %>%
  dplyr::mutate(
    # On the logistic fits the coefficient is a log odds ratio; on the Cox
    # fits a log hazard ratio. The exponentiated column is labelled by family
    # rather than called "HR" everywhere, which would be wrong for three rows.
    effect_PROLIF = exp(estimate_PROLIF),
    effect_lo_PROLIF = exp(ci_lo_PROLIF), effect_hi_PROLIF = exp(ci_hi_PROLIF),
    scale = ifelse(family == "logistic", "odds ratio", "hazard ratio"))
.pr(PROLIF_TAB %>% dplyr::transmute(set, family, n, events,
    PROLIF = round(estimate_PROLIF, 4), se = round(se_PROLIF, 4),
    ci = paste0("[", sprintf("%.4f", ci_lo_PROLIF), ", ",
                sprintf("%.4f", ci_hi_PROLIF), "]"),
    excl_0 = excl_0_PROLIF, p = signif(p_PROLIF, 3),
    vif = round(vif_PROLIF, 3),
    OX = round(estimate_OX, 4), OX_excl_0 = excl_0_OX))
message("\n   the same rows on the exponentiated scale, labelled by family:")
.pr(PROLIF_TAB %>% dplyr::transmute(set, scale,
    PROLIF_effect = round(effect_PROLIF, 4),
    ci = paste0("[", sprintf("%.4f", effect_lo_PROLIF), ", ",
                sprintf("%.4f", effect_hi_PROLIF), "]"),
    OX_effect = round(exp(estimate_OX), 4)))

# A MULTI-DF GVIF IS NOT ON THE SAME SCALE AS A 1-DF VIF, and section 5.3's
# threshold of 5 was set for OX, which is 1 df. Fox and Monette's prescription
# is to compare GVIF^(1/(2*Df)) against the square root of the 1-df threshold.
# Comparing a raw multi-df GVIF against 5 overstates it: GSE194040's subtype
# term is 7.98 on 3 df, which is 1.414 adjusted, and its treatment term 7.27 on
# 12 df, which is 1.086.
# gvif_adj is NOT added to the saved $vif_vectors, so that table's digest is
# unchanged; it is a deterministic function of the two columns that are saved,
# vif and vif_df, and is recomputable as vif^(1/(2*vif_df)).
VIFS_ADJ <- VIFS %>%
  dplyr::mutate(gvif_adj = vif^(1 / (2 * vif_df)))
VIF_HI <- VIFS_ADJ %>%
  dplyr::filter(!is.na(gvif_adj), gvif_adj > sqrt(VIF_CAVEAT))
message("\n   collinearity, read on the df-ADJUSTED scale GVIF^(1/(2*Df)) ",
        "against sqrt(", VIF_CAVEAT, ") = ", round(sqrt(VIF_CAVEAT), 3), ":")
if (nrow(VIF_HI)) {
  .pr(VIF_HI %>% dplyr::transmute(set, term, gvif = round(vif, 3),
                                  df = vif_df, gvif_adj = round(gvif_adj, 3)))
} else {
  message("   NO term in any model exceeds it. The highest adjusted value ",
          "anywhere is ", signif(max(VIFS_ADJ$gvif_adj, na.rm = TRUE), 4),
          ", and the highest among the OX and PROLIF terms - both 1 df, so ",
          "GVIF == VIF - is ",
          signif(max(VIFS_ADJ$vif[VIFS_ADJ$term %in% c("OX", "PROLIF")]), 4),
          ".")
  message("   ", .sub(paste0(
    "The two largest RAW GVIFs are GSE194040's spine terms - subtype 7.98 on ",
    "3 df and treatment 7.27 on 12 df - which adjust to 1.414 and 1.086. ",
    "That is I-SPY2's arm assignment by receptor status, it concerns the ",
    "spine and not the exposures, and on the comparable scale it is modest.")))
}
message("   ", .sub(paste0(
  "14.4: MYC is absent from m1, so rho(MYC, PROLIF_DISJOINT) of 0.78 to 0.82 ",
  "in the pCR cohorts and 0.71 in METABRIC does NOT bear on this ",
  "coefficient. What governs it is rho(OX, PROLIF) of 0.36 to 0.49, inside ",
  "the declaration's separable band. The VIF vector above is the measured ",
  "version of that statement.")))

# =============================================================================
# 7. D. META-ANALYSE the PROLIF coefficient across the three pCR cohorts
# =============================================================================
# script 13's .meta(), reproduced verbatim as E40 and E42 did, so the pooling
# is the same DerSimonian-Laird estimator this arm has used throughout.
message("\n7. D. meta-analysis of PROLIF across the three pCR cohorts")

.meta <- function(est, se, label) {
  v <- se^2; w <- 1 / v; k <- length(est)
  if (k < 2L) {
    return(tibble::tibble(term = label, k = k,
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
  tibble::tibble(term = label, k = k,
    fe_estimate = fe, fe_ci_lo = fe - CRIT*sqrt(vfe),
    fe_ci_hi = fe + CRIT*sqrt(vfe), fe_p = 2*stats::pnorm(-abs(fe/sqrt(vfe))),
    re_estimate = re, re_ci_lo = re - CRIT*sqrt(vre),
    re_ci_hi = re + CRIT*sqrt(vre), re_p = 2*stats::pnorm(-abs(re/sqrt(vre))),
    Q = Q, Q_p = stats::pchisq(Q, df, lower.tail = FALSE),
    tau2 = tau2, I2 = 100*max(0, (Q - df)/Q))
}

PCR_COEFS <- COEFS %>% dplyr::filter(set %in% COHORTS)
META <- dplyr::bind_rows(lapply(c("PROLIF", "OX"), function(tm) {
  d <- PCR_COEFS[PCR_COEFS$term == tm, ]
  d <- d[match(COHORTS, d$set), ]
  m <- .meta(d$estimate, d$se, tm)
  m$cohorts     <- paste(d$set, collapse = ", ")
  m$signs       <- paste(sprintf("%s %s", d$set,
                                 ifelse(d$estimate > 0, "+", "-")),
                         collapse = " ")
  m$signs_agree <- length(unique(sign(d$estimate))) == 1L
  m$summary_is  <- ifelse(!is.na(m$I2) & m$I2 >= 50, "random effects",
                          "fixed effect")
  m
}))
.pr(META %>% dplyr::transmute(term, k,
    fe = round(fe_estimate, 4),
    fe_ci = paste0("[", sprintf("%.4f", fe_ci_lo), ", ",
                   sprintf("%.4f", fe_ci_hi), "]"),
    fe_p = signif(fe_p, 3),
    re = round(re_estimate, 4),
    re_ci = paste0("[", sprintf("%.4f", re_ci_lo), ", ",
                   sprintf("%.4f", re_ci_hi), "]"),
    Q = round(Q, 3), Q_p = signif(Q_p, 3), tau2 = round(tau2, 4),
    I2 = round(I2, 1), summary_is, signs_agree, signs))
message("   ", .sub(paste0(
  "BOTH summaries are reported always; at I2 >= 50% the random-effects ",
  "estimate is the summary and the fixed effect sits beside it, never ",
  "instead of it (section 9). OX is pooled here too, as the comparator the ",
  "PROLIF row is read against - and its m1 row must reproduce E40's pooled ",
  "-0.161555, which is asserted next.")))
ox_meta <- META[META$term == "OX", ]
e40_meta <- e40$meta[e40$meta$rung == RUNG, ]
.stop_if(abs(ox_meta$fe_estimate - e40_meta$fe_estimate) < REPRO_TOL,
         "the re-pooled OX fixed effect is ", signif(ox_meta$fe_estimate, 10),
         " against E40's stored ", signif(e40_meta$fe_estimate, 10),
         ". The pooling method has drifted. STOP.")
message("   ASSERTED: the re-pooled OX fixed effect reproduces E40's stored ",
        "value to better than ", REPRO_DP, " decimal places")

# =============================================================================
# 8. 14.6's reading rule, printed UNFILLED. This script does NOT classify.
# =============================================================================
message("\n8. declaration 14.6, printed UNFILLED")

READING_RULE <- tibble::tribble(
  ~outcome, ~`what may be written`,
  "PROLIF positive on pCR and above unity on survival, both intervals excluding the null",
    "the sign pair is present in these models; the Results may state the inversion as an internal observation",
  "the directions hold but one or both intervals include the null",
    "report the point estimates and directions; cite the published pattern for the relationship itself; do not claim the pair is demonstrated here",
  "a direction contradicts the published pattern",
    "report it, and withdraw the inversion framing from the Results")
.pr(READING_RULE)
message("   ", .sub(paste0(
  "THE THIRD ROW MAY NOT BE RESCUED BY REPORTING ONLY THE ENDPOINT THAT ",
  "COOPERATED. Both are retrieved and both are reported (14.6). This script ",
  "does not classify the outcome against this table.")))
message("   ", .sub(paste0(
  "AND 14.5, stated before the fit: subtype is in the spine and subtype is ",
  "substantially a proliferation classification, so much of proliferation's ",
  "signal is absorbed before the PROLIF term is reached. AN ATTENUATED OR ",
  "NULL PROLIF COEFFICIENT IS EXPECTED AND DOES NOT REFUTE THE PUBLISHED ",
  "PATTERN - it would say proliferation adds little WITHIN SUBTYPE, not that ",
  "the relationship is absent.")))

# =============================================================================
# 9. Warnings
# =============================================================================
message("\n9. warnings")
message("   ", .sub(paste0(
  "car::vif raised the spurious 'No intercept: vifs may not be sensible' ",
  .VIF_STATE$no_intercept, " time(s), counted and muffled deliberately. It ",
  "arises on the COX fits only - a Cox model has no intercept by design and ",
  "the VIF needs none - and NOT on the logistic fits, which have one. Same ",
  "artefact identified in E40 and E41.")))
OTHER_TAB <- if (length(OTHER_WARNINGS)) dplyr::bind_rows(OTHER_WARNINGS) else
  tibble::tibble(where = character(0), message = character(0))
if (nrow(OTHER_TAB)) {
  message("   OTHER warnings, reported and NOT muffled away:")
  .pr(OTHER_TAB %>% dplyr::count(where, message, name = "n"))
} else {
  message("   no other warning was raised.")
}
DROPPED_TAB <- if (length(DROPPED)) dplyr::bind_rows(DROPPED) else
  tibble::tibble(set = character(0), term = character(0),
                 dropped = character(0), n_dropped = integer(0))
if (nrow(DROPPED_TAB)) {
  message("   empty factor levels dropped EXPLICITLY (0 observations each, ",
          "so no estimate changes):")
  .pr(DROPPED_TAB)
}

# =============================================================================
# 10. E. Save - INCLUDING THE FITTED MODEL OBJECTS, per 14.7
# =============================================================================
message("\n10. E. saving")

readr::write_csv(COEFS,      file.path(DIR_TABLES, "E43_coefficient_matrix.csv"))
readr::write_csv(VIFS,       file.path(DIR_TABLES, "E43_vif_vectors.csv"))
readr::write_csv(PROLIF_TAB, file.path(DIR_TABLES, "E43_prolif_focused.csv"))
readr::write_csv(META,       file.path(DIR_TABLES, "E43_prolif_meta.csv"))
readr::write_csv(REPRO,      file.path(DIR_TABLES, "E43_reproduction_check.csv"))
readr::write_csv(AUDIT,      file.path(DIR_TABLES, "E43_spine_audit.csv"))

# 14.7 asks for the fitted objects, and they are saved - but a fitted model is
# NOT digest-stable across runs, so a verification pass must not read that as a
# failure. coxph and glm objects carry `terms`, `formula`, `model` and `family`,
# each of which holds an ENVIRONMENT whose address differs between sessions.
# Verified: across two runs the coefficients and the variance matrices are
# identical in every fit and only those four carriers differ.
# So $fits is EXCLUDED from a digest comparison and $fit_digests is the stable
# handle instead: a per-fit digest of coef() and vcov(), which are the numbers.
FIT_DIGESTS <- dplyr::bind_rows(lapply(names(FITS), function(k)
  tibble::tibble(set = k, class = class(FITS[[k]])[1],
                 n_coef = length(stats::coef(FITS[[k]])),
                 digest_coef = digest::digest(stats::coef(FITS[[k]])),
                 digest_vcov = digest::digest(stats::vcov(FITS[[k]])))))
message("\n   $fits is saved per 14.7 but is NOT digest-stable across runs ",
        "(model objects carry environments).")
message("   $fit_digests is the stable handle - coef() and vcov() per fit:")
.pr(FIT_DIGESTS %>% dplyr::transmute(set, class, n_coef,
    coef = substr(digest_coef, 1, 8), vcov = substr(digest_vcov, 1, 8)))

saveRDS(list(
  # 14.7's standing fix, applied here first: the FULL matrix, the FULL vector,
  # and the fitted objects themselves.
  coefficient_matrix = COEFS,
  vif_vectors        = VIFS,
  fits               = FITS,
  fit_digests        = FIT_DIGESTS,
  prolif_focused     = PROLIF_TAB,
  meta               = META,
  reproduction_check = REPRO,
  audit              = AUDIT,
  model_info         = META_INFO,
  dropped_levels     = DROPPED_TAB,
  other_warnings     = OTHER_TAB,
  vif_no_intercept_warnings = .VIF_STATE$no_intercept,
  reading_rule       = READING_RULE,
  frames_used        = list(
    pcr = "results/e39_prepare_and_gate.rds$frame + the FROZEN h4_outcome_models.rds$frames",
    metabric = "results/e41_metabric_ladder.rds$frame (added at 51465b6)"),
  spec = list(
    declaration = paste0("docs/2026-10-02_E39_respiratory_axis_decomposition",
                         "_declaration.md section 14, committed at 9fb7852"),
    # 14.2, verbatim.
    posture = paste0(
      "This is the Block C situation again. ",
      "2026-08-29_escape_reading_declaration.md fixed that the discarded main ",
      "effects were 'relevant to a constraint reading, and retrieving them is ",
      "a one-line change,' and that 'when retrieved they inform E2's ",
      "interpretation without redefining it.' The same holds here: the model ",
      "is unchanged, the term was declared, and only the extraction is new. ",
      "BUT THE READING IT SUPPORTS IS POST-HOC AND MUST BE WRITTEN AS SUCH. ",
      "Section 9.1's sign-pair table anticipated better/worse, ",
      "no-effect/worse and better/better. The observed ",
      "worse-pCR-with-better-survival combination is not in it. That gap is ",
      "recorded here rather than repaired retrospectively."),
    # 14.5, verbatim.
    attenuation = paste0(
      "Subtype is in the spine, and subtype is substantially a proliferation ",
      "classification - PAM50 LumA against LumB in METABRIC, ",
      "hormone-receptor status in the neoadjuvant cohorts. Much of ",
      "proliferation's prognostic and predictive signal is therefore ",
      "absorbed before the PROLIF term is reached. AN ATTENUATED OR NULL ",
      "PROLIF COEFFICIENT IS EXPECTED AND DOES NOT REFUTE THE PUBLISHED ",
      "PATTERN. It would say that proliferation adds little within subtype, ",
      "not that the relationship is absent."),
    # 14.6, verbatim.
    reading_rule = paste0(
      "PROLIF positive on pCR and above unity on survival, both intervals ",
      "excluding the null -> the sign pair is present in these models; the ",
      "Results may state the inversion as an internal observation. The ",
      "directions hold but one or both intervals include the null -> report ",
      "the point estimates and directions; cite the published pattern for ",
      "the relationship itself; do not claim the pair is demonstrated here. ",
      "A direction contradicts the published pattern -> report it, and ",
      "withdraw the inversion framing from the Results. THE THIRD ROW MAY ",
      "NOT BE RESCUED BY REPORTING ONLY THE ENDPOINT THAT COOPERATED. Both ",
      "are retrieved and both are reported."),
    rung = paste0(
      "m1 ONLY (14.4). MYC is absent from m1, so rho(MYC, PROLIF_DISJOINT) ",
      "of 0.78 / 0.79 / 0.82 in the pCR cohorts and 0.71 in METABRIC does ",
      "not bear on this coefficient; it is governed by rho(OX, PROLIF) of ",
      "0.36 to 0.49, inside the separable band. m0, m2 and m3 were NOT ",
      "retrieved."),
    reproduction = paste0(
      "14.3's stop condition. Every refitted OX coefficient AND its SE ",
      "reproduced the stored E40/E41 value to better than ", REPRO_DP,
      " decimal places, with n and events matching, and the re-pooled OX ",
      "fixed effect reproduced E40's stored pooled value. So these are the ",
      "same models E40 and E41 fitted and the PROLIF coefficient is the one ",
      "they computed and discarded. A refit that did not reproduce would ",
      "have stopped the script without reporting a PROLIF estimate."),
    no_rescoring = paste0(
      "Nothing was rescored and METABRIC_DATA.RData was not reloaded. The ",
      "pCR fits came from e39_prepare_and_gate.rds$frame plus the FROZEN ",
      "h4_outcome_models.rds$frames; the METABRIC fits from ",
      "e41_metabric_ladder.rds$frame. Both subtype definitions and both ",
      "spine audits are E40's and E41's, reused unchanged."),
    standing_fix = paste0(
      "14.7. This script saves the FULL coefficient matrix, the FULL VIF ",
      "vector and the FITTED MODEL OBJECTS. It is the first in the arm to do ",
      "so, and it exists because selective extraction has now forced three ",
      "retrievals. ONE CAVEAT FOR ANY VERIFICATION PASS: $fits is NOT ",
      "digest-stable across runs, because coxph and glm objects carry terms, ",
      "formula, model and family, each holding an environment whose address ",
      "differs between sessions. Across two dry runs the coefficients and ",
      "variance matrices were identical in every fit and only those carriers ",
      "differed. EXCLUDE $fits from a digest comparison and use $fit_digests, ",
      "which digests coef() and vcov() per fit - the numbers rather than the ",
      "wrapper."),
    instrument = paste0(
      "SINGLE-INSTRUMENT, as the whole line is. No coefficient is compared ",
      "across cohorts as a raw value (4.2). The pCR coefficients are log ",
      "odds ratios and the METABRIC ones log hazard ratios; the ",
      "exponentiated column is labelled by family, not called HR ",
      "throughout."),
    frozen = paste0("Nothing written to myc_human_validation (d3ac60e) or ",
                    "myc_mouse. h4_outcome_models.rds was read READ-ONLY.")),
  built = Sys.time()), PATH_E43)

message("   ", PATH_E43)
message("   6 tables in ", DIR_TABLES)
message("   FITTED MODEL OBJECTS saved: ", length(FITS), " (14.7)")

message("\nE43 done. Numbers only; nothing here is interpreted or classified.")
message("A RETRIEVAL. The model is unchanged (14.2).\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E43)

  # B. THE GATE. Read this before any coefficient below it.
  x$reproduction_check %>% dplyr::transmute(set, n, events,
      refit_estimate, stored_estimate, d_estimate, d_se, ok) %>%
    as.data.frame()

  # C. the focused PROLIF table, OX beside it.
  x$prolif_focused %>% dplyr::transmute(set, family, n, events,
      PROLIF = round(estimate_PROLIF, 4), ci_lo = round(ci_lo_PROLIF, 4),
      ci_hi = round(ci_hi_PROLIF, 4), excl_0 = excl_0_PROLIF,
      p = signif(p_PROLIF, 3), vif = round(vif_PROLIF, 3),
      OX = round(estimate_OX, 4)) %>% as.data.frame()

  # the FULL matrix and the FULL vector - what 14.7 exists to preserve.
  x$coefficient_matrix %>% as.data.frame()
  x$vif_vectors %>% as.data.frame()

  # D. pooled across the three pCR cohorts, PROLIF and OX side by side.
  x$meta %>% dplyr::select(term, k, fe_estimate, fe_ci_lo, fe_ci_hi, fe_p,
                           re_estimate, Q_p, I2, summary_is, signs_agree,
                           signs) %>% as.data.frame()

  # 14.6, UNFILLED. The script does not classify and neither does the object.
  x$reading_rule %>% as.data.frame()
  cat(x$spec$attenuation, "\n")
  cat(x$spec$reading_rule, "\n")

  # the fitted objects themselves, saved for the first time in this arm.
  names(x$fits)
  summary(x$fits[["GSE194040"]])
  summary(x$fits[["METABRIC / ER-pos (IHC)"]])

  # $fits is NOT digest-stable; $fit_digests is the handle that is.
  x$fit_digests %>% as.data.frame()

  x$dropped_levels %>% as.data.frame()
  x$other_warnings %>% as.data.frame()
  utils::str(x$spec)
}
