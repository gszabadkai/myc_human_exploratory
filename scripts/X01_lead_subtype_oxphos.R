# X01_lead_subtype_oxphos.R
# =============================================================================
# EXPLORATORY LEADS. NOT PART OF THE DECLARED E-LINE. MANUSCRIPT-EXCLUDED.
# =============================================================================
# Note: docs/2026-10-03_X01_exploratory_leads_note.md, committed at c42e3b0.
#
# POSTURE, from the note's front matter and printed as this script's first
# output so nobody reads a number before reading this:
#
#   No result in this note or in scripts/X01_lead_subtype_oxphos.R may enter
#   any manuscript document, figure, legend or supplementary file without its
#   own declaration written and committed first. It carries no reading rule, no
#   pass, no fail and no falsification criterion. It is a survey of where a
#   predictive-marker hypothesis might be worth building, and nothing more.
#
# WHY IT IS OFF THE E-NUMBERING. The E-numbered line is manuscript-bound:
# every E-script's result is reported whatever it shows, under a declaration
# committed before it ran. X01 is not. Keeping it off the numbering makes the
# difference visible in the file listing so nobody downstream mistakes it for
# part of the declared work.
#
# =============================================================================
# WHAT THIS CANNOT LICENSE, EVER (note section 5)
# =============================================================================
#   - No manuscript sentence, figure, legend or supplementary table.
#   - No claim that any subtype or treatment stratum "shows" or "demonstrates"
#     anything. The vocabulary here is SUGGESTS, IS CONSISTENT WITH, WOULD BE
#     WORTH TESTING.
#   - No comparison of any coefficient across cohorts as a raw value (E39
#     section 4.2 still holds).
#   - No reinstatement of anything on the retired list.
#   - No interaction test, no formal effect-modification test. Not because of
#     posture but because the event counts do not support one.
#
# AND the prohibition X01 does NOT relax: E42 section 13.2 forbids any subtype
# result becoming a finding in its own right. MOVING THE WORK INTO THIS NOTE
# DOES NOT RELAX IT.
#
# =============================================================================
# THE THREE CAUTIONS (note section 2)
# =============================================================================
# 1. WINNER'S CURSE. With a pooled METABRIC estimate of HR 0.896, individual
#    PAM50 subtypes scatter by sampling variation alone. The largest stratum
#    estimate will look impressive and will be partly the maximum of several
#    noisy draws. NO SUBTYPE IS THE LEAD MERELY FOR BEING THE LARGEST.
#    Section 7 reports the scatter beside the leads for this reason.
# 2. THE LUMINAL pCR ESTIMATE IS PROBABLY INFLATED. E42 declared the
#    detectable effect for HRpos_HER2neg at 0.312 log-odds; the observed
#    pooled estimate, 0.276, is SMALLER than the effect the stratum was
#    declared able to detect. The true luminal effect is likely nearer E40's
#    pooled -0.162 than -0.276.
# 3. THE s508 EXPLANATION ALREADY FAILED ONCE. E42 measured
#    rho(OX, subtype2) = -0.054 and refuted the proposed explanation. The
#    plausible story was wrong. Treat every mechanism the same way.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It does NOT rescore anything, does NOT reload METABRIC_DATA.RData and
#     does NOT refit anything already in e41_metabric_ladder.rds. It reads
#     that object's $frame and stratifies it.
#   - It fits m0 and m1 ONLY. NO m2, NO m3.
#   - NO interaction term of any kind.
#   - It does NOT describe the ER.Expr repeat as replication, validation or
#     confirmation. It is THE SAME COHORT REPARTITIONED.
#   - It does NOT compare any coefficient across cohorts as a raw value.
#   - It does NOT interpret. It prints numbers, never verdicts.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
#
# SINGLE-INSTRUMENT. METABRIC is an array; mitoPPS is unavailable.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

POSTURE <- paste0(
  "No result in this note or in scripts/X01_lead_subtype_oxphos.R may enter ",
  "any manuscript document, figure, legend or supplementary file without its ",
  "own declaration written and committed first. It carries no reading rule, ",
  "no pass, no fail and no falsification criterion. It is a survey of where a ",
  "predictive-marker hypothesis might be worth building, and nothing more.")

.sub <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n  ")

message("\n", strrep("=", 78))
message("X01 - EXPLORATORY LEADS. NOT PART OF THE DECLARED E-LINE.")
message("MANUSCRIPT-EXCLUDED. Read this before any number below.")
message(strrep("=", 78))
message("  ", .sub(POSTURE))
message(strrep("-", 78))
message("  ", .sub(paste0(
  "E42 section 13.2 forbids any subtype result becoming a finding in its own ",
  "right. MOVING THE WORK INTO THIS NOTE DOES NOT RELAX IT. Vocabulary here ",
  "is SUGGESTS, IS CONSISTENT WITH, WOULD BE WORTH TESTING - never shows, ",
  "demonstrates or establishes.")))
message(strrep("=", 78))

PATH_E41 <- file.path(DIR_RESULTS, "e41_metabric_ladder.rds")
PATH_X01 <- file.path(DIR_RESULTS, "x01_lead_subtype_oxphos.rds")

# =============================================================================
# 0. CONSTANTS - from E41 and from the note, none invented here
# =============================================================================
CI_LEVEL <- 0.95
CRIT     <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

# The rungs. m0 AND m1, because note section 3 item 3 asks whether the
# proliferation adjustment does work inside a subtype: a subtype where m0 is
# null and m1 is protective is a better lead than one where both are
# protective. m2 and m3 are NOT fitted.
RUNG_RHS <- list(m0 = "OX", m1 = c("OX", "PROLIF"))
RUNGS    <- names(RUNG_RHS)
SPINE_CANDIDATES <- c("subtype_term_x01", "treatment")

# E41's declared counts, re-asserted. The survival naming is counterintuitive
# (10.3a): the 623-event object is the CAUSE-SPECIFIC one and is primary.
N_ALL      <- 1981L
E_BCSS     <- 623L
N_ER_POS   <- 1505L
E_ER_POS   <- 431L
ER_PRIMARY <- "er_primary"
ER_SENS    <- "er_sensitive"

# The threshold, the SAME one E42 named, mapped to a Cox model. E42's rule was
# "minority class at least 5 AND events per parameter at least 5", with
# minority class used because a LOGISTIC fit is limited by the rarer outcome.
# A COX MODEL HAS NO MINORITY CLASS: the limiting count is the number of
# EVENTS. So the rule here is events >= 5 and events per parameter >= 5, and
# that mapping is stated rather than applied silently.
MIN_EVENTS <- 5L
EPV_FLOOR  <- 5
EPV_FLAG   <- 10

VIF_CAVEAT <- 5          # E39 section 5.3

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.n_lev   <- function(v) dplyr::n_distinct(v[!is.na(v)])

# =============================================================================
# 1. Inputs - E41's saved frame, NOTHING rescored
# =============================================================================
message("\n1. inputs (E41's saved frame; NOTHING is rescored or reloaded)")

.stop_if(file.exists(PATH_E41), "absent: ", PATH_E41,
         ". Source scripts/E41_metabric_ladder.R first.")
e41 <- readRDS(PATH_E41)
.stop_if(!is.null(e41$frame),
         PATH_E41, " carries no $frame. It was added to E41's save block on ",
         "2026-10-03 and the object on disk predates that. RE-SOURCE ",
         "scripts/E41_metabric_ladder.R, then re-source this script. X01 will ",
         "NOT reload METABRIC_DATA.RData or rescore, because E41 already did ",
         "that work.")

FR <- e41$frame
.stop_if(nrow(FR) == N_ALL, "$frame is ", nrow(FR), " rows, not ", N_ALL)
for (v in c("sample", "er_primary", "er_sensitive", "subtype_pam50",
            "treatment", "bcss_time_years", "bcss_event", "OX", "PROLIF")) {
  .stop_if(v %in% names(FR), "$frame is missing ", v)
}
# The endpoint, asserted again. E41's 10.3a: the object with the
# generic-sounding name is the ALL-CAUSE one; the 623-event one is
# cause-specific and is what this script uses.
.stop_if(sum(FR$bcss_event) == E_BCSS,
         "$frame carries ", sum(FR$bcss_event), " bcss events, not the ",
         E_BCSS, " E41 asserted. The WRONG ENDPOINT may have been saved. STOP.")
message("   $frame: ", nrow(FR), " patients | bcss events ", sum(FR$bcss_event),
        " (CAUSE-SPECIFIC, the primary endpoint; 10.3a's naming trap)")

# The ER-positive IHC stratum, which is what the note surveys.
DP <- FR[!is.na(FR[[ER_PRIMARY]]) & FR[[ER_PRIMARY]] == "P", ]
.stop_if(nrow(DP) == N_ER_POS && sum(DP$bcss_event) == E_ER_POS,
         "the ER-positive IHC stratum is ", nrow(DP), " / ",
         sum(DP$bcss_event), ", not E41's ", N_ER_POS, " / ", E_ER_POS)
DS <- FR[FR[[ER_SENS]] == "P", ]
message("   ER-positive (IHC, PRIMARY): ", nrow(DP), " patients, ",
        sum(DP$bcss_event), " events")
message("   ER-positive (ER.Expr, the INTERNAL CHECK): ", nrow(DS),
        " patients, ", sum(DS$bcss_event), " events")

# The pooled estimate every subtype is read against, taken FROM THE SAVED
# OBJECT rather than retyped.
POOLED <- e41$er_stratified %>%
  dplyr::filter(set == "ER-pos (IHC)", rung == "m1")
.stop_if(nrow(POOLED) == 1L, "could not find E41's ER-pos (IHC) m1 row")
message("   E41's pooled ER-positive m1, for reference: HR ",
        round(POOLED$HR, 4), " [", round(POOLED$hr_lo, 3), ", ",
        round(POOLED$hr_hi, 4), "], p ", signif(POOLED$p, 3))

# PAM50 is the stratifier here. It is NOT E41's spine subtype term renamed:
# inside a PAM50 stratum it is constant, which is the point.
FR$subtype_term_x01 <- factor(FR$subtype_pam50)
DP$subtype_term_x01 <- factor(DP$subtype_pam50)
DS$subtype_term_x01 <- factor(DS$subtype_pam50)

# =============================================================================
# 2. A. COUNTS FIRST - the stop-and-check, before any fit
# =============================================================================
message("\n2. A. COUNTS FIRST. Read this block before any coefficient.")

.counts <- function(d, label) {
  d %>% dplyr::group_by(subtype = as.character(subtype_term_x01)) %>%
    dplyr::summarise(set = label, n = dplyr::n(),
                     events = sum(bcss_event),
                     median_fu_y = stats::median(bcss_time_years),
                     treat_levels = .n_lev(treatment), .groups = "drop") %>%
    # params on the LARGER of the two rungs, m1, which is the binding one
    dplyr::mutate(
      params_m1 = 2L + ifelse(treat_levels <= 1L, 0L, treat_levels - 1L),
      epv       = events / params_m1,
      fit       = events >= MIN_EVENTS & epv >= EPV_FLOOR,
      epv_flag  = epv < EPV_FLAG) %>%
    dplyr::arrange(dplyr::desc(n))
}
COUNTS <- .counts(DP, "ER-pos (IHC)")
.pr(COUNTS %>% dplyr::mutate(median_fu_y = round(median_fu_y, 2),
                             epv = round(epv, 1)))

message("   ", .sub(paste0(
  "THE THRESHOLD IS E42's, MAPPED TO A COX MODEL AND THE MAPPING STATED. ",
  "E42 required a minority class of at least ", MIN_EVENTS, " and at least ",
  EPV_FLOOR, " events per parameter, using the MINORITY CLASS because a ",
  "logistic fit is limited by the rarer outcome. A COX MODEL HAS NO MINORITY ",
  "CLASS - its limiting count is EVENTS - so the rule here is events >= ",
  MIN_EVENTS, " and events per parameter >= ", EPV_FLOOR,
  ". Below ", EPV_FLAG, " events per parameter a stratum is FLAGGED but still ",
  "fitted. Parameters are counted on m1, the larger rung.")))

NOT_FIT <- COUNTS %>% dplyr::filter(!fit)
if (nrow(NOT_FIT)) {
  message("\n   NOT FITTED, below the threshold:")
  .pr(NOT_FIT %>% dplyr::transmute(subtype, n, events, treat_levels,
                                   params_m1, epv = round(epv, 1)))
  message("   ", .sub(paste0(
    "READ WHY EACH FAILS, because the reasons differ. A subtype can fail on ",
    "EVENT SCARCITY or on PARAMETER COUNT, and the remedy would differ. ",
    "Where events are comfortably above ", MIN_EVENTS,
    " but the treatment variable carries many levels, the stratum fails ",
    "because the SPINE is expensive, not because the data are thin. NO ",
    "REDUCED SPINE IS FITTED HERE - that would be a modelling decision this ",
    "survey does not make.")))
  .pr(NOT_FIT %>% dplyr::transmute(subtype, events, params_m1,
        epv = round(epv, 1),
        fails_on = dplyr::case_when(
          events < MIN_EVENTS ~ "event scarcity",
          events >= 2 * EPV_FLOOR & params_m1 >= 6L ~
            "PARAMETER COUNT - events alone would suffice for a cheaper spine",
          TRUE ~ "both event count and parameter count")))
}
FLAGGED <- COUNTS %>% dplyr::filter(fit, epv_flag)
if (nrow(FLAGGED)) {
  message("\n   FITTED BUT FLAGGED, under ", EPV_FLAG, " events per parameter:")
  .pr(FLAGGED %>% dplyr::transmute(subtype, n, events, params_m1,
                                   epv = round(epv, 1)))
}

# Approximate detectable HR per subtype. Schoenfeld's approximation for a
# STANDARDISED continuous covariate: se(log HR) ~ 1 / sqrt(events), so the
# detectable |log HR| is (z_{1-a/2} + z_{power}) / sqrt(events). IT IGNORES
# COVARIATE ADJUSTMENT, which inflates the variance, so the true detectable
# effect is LARGER. Labelled as an approximation throughout.
Z_A <- stats::qnorm(0.975); Z_B <- stats::qnorm(0.80)
DETECTABLE <- COUNTS %>%
  dplyr::transmute(subtype, n, events,
                   detectable_logHR = (Z_A + Z_B) / sqrt(events),
                   detectable_HR_protective = exp(-(Z_A + Z_B) / sqrt(events)),
                   with_2x_variance = (Z_A + Z_B) * sqrt(2) / sqrt(events))
message("\n   APPROXIMATE detectable effect per subtype (80% power, alpha ",
        "0.05, per 1 within-cohort SD of OX):")
.pr(DETECTABLE %>% dplyr::mutate(
  detectable_logHR = round(detectable_logHR, 3),
  detectable_HR_protective = round(detectable_HR_protective, 3),
  with_2x_variance = round(with_2x_variance, 3)))
message("   ", .sub(paste0(
  "APPROXIMATE, and an UNDERESTIMATE of what is needed: it ignores covariate ",
  "adjustment. For scale, E41's pooled ER-positive m1 is HR ",
  round(POOLED$HR, 3), " (log ", round(POOLED$estimate, 3), "). A NULL IN A ",
  "LOW-EVENT SUBTYPE IS UNINFORMATIVE AND IS NOT ABSENCE - it says the ",
  "survey could not look there, not that there is nothing to see.")))

# =============================================================================
# 3. B. SPINE AUDIT per subtype stratum
# =============================================================================
# Inside a PAM50 stratum the subtype term is constant BY CONSTRUCTION and is
# omitted EXPLICITLY. treatment is audited too and may be constant or
# near-constant in a small stratum. Empty factor levels are dropped explicitly
# and reported, as E41 does - an empty level contributes nothing, gives an NA
# coefficient and breaks car::vif on aliased coefficients.
message("\n3. B. spine audit, per subtype stratum")

DROPPED_LEVELS <- list()
.drop_empty <- function(d, label, spine) {
  for (v in spine) {
    if (!is.factor(d[[v]])) next
    gone <- setdiff(levels(d[[v]]), unique(as.character(d[[v]])))
    if (length(gone)) {
      DROPPED_LEVELS[[length(DROPPED_LEVELS) + 1L]] <<- tibble::tibble(
        set = label, term = v, dropped = paste(gone, collapse = ", "),
        n_dropped = length(gone), n_kept = nlevels(droplevels(d[[v]])))
      d[[v]] <- droplevels(d[[v]])
    }
  }
  d
}

.audit_one <- function(d, label) {
  tibble::tibble(set = label, n = nrow(d),
                 term = c("OX", SPINE_CANDIDATES),
                 n_levels = c(.n_lev(d$OX), .n_lev(d$subtype_term_x01),
                              .n_lev(d$treatment))) %>%
    dplyr::mutate(constant = n_levels <= 1L,
                  action = ifelse(constant, "OMITTED (constant)", "fitted"))
}
FIT_SUBTYPES <- COUNTS$subtype[COUNTS$fit]
AUDIT <- dplyr::bind_rows(lapply(FIT_SUBTYPES, function(st)
  .audit_one(DP[as.character(DP$subtype_term_x01) == st, ],
             paste0("ER-pos (IHC) / ", st))))
.pr(AUDIT)
st_rows <- AUDIT[AUDIT$term == "subtype_term_x01", ]
.stop_if(all(st_rows$constant),
         "the PAM50 term is NOT constant inside ", sum(!st_rows$constant),
         " stratum/strata. The split is wrong. STOP.")
message("   ASSERTED: the PAM50 term is constant in all ", nrow(st_rows),
        " fitted strata and is OMITTED EXPLICITLY in every one")
tr_const <- AUDIT[AUDIT$term == "treatment" & AUDIT$constant, ]
if (nrow(tr_const)) {
  message("   treatment is CONSTANT and therefore OMITTED EXPLICITLY in:")
  .pr(tr_const %>% dplyr::transmute(set, n, n_levels))
} else {
  message("   treatment varies in every fitted stratum and is fitted in each")
}

.spine_for <- function(label) {
  a <- AUDIT[AUDIT$set == label & AUDIT$term %in% SPINE_CANDIDATES, ]
  a$term[!a$constant]
}

# =============================================================================
# 4. The estimator helpers
# =============================================================================
# car::vif WILL raise "No intercept: vifs may not be sensible" on every Cox
# fit - a Cox model has no intercept BY DESIGN and the VIF comes from the
# correlation matrix of the coefficient estimates, which needs none. It is
# counted and muffled DELIBERATELY, as E41 did, and the count is reported.
.VIF_STATE <- new.env(parent = emptyenv())
.VIF_STATE$no_intercept_warnings <- 0L
OTHER_WARNINGS <- list()

# Per-fit convergence flag. A Cox spine level that has observations but ZERO
# EVENTS gives an infinite coefficient, and coxph says so in a warning. The
# OX coefficient remains estimable, but the fit's spine is partly unidentified
# and the row must not be read as clean.
#
# THE LEVEL IS NOT DROPPED. For a logistic fit an EMPTY level contributes
# nothing and dropping it is bookkeeping (E41 does exactly that). A zero-EVENT
# level in a Cox model is different: those patients are at risk and then
# censored, so they DO contribute to other levels' risk sets, and removing
# them would change the data rather than tidy it. So the fit is kept and
# FLAGGED.
.LAST_FIT <- new.env(parent = emptyenv())
.catch <- function(expr, label) {
  .LAST_FIT$warned <- FALSE
  .LAST_FIT$msg <- ""
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (!grepl("No intercept", m, fixed = TRUE)) {
      .LAST_FIT$warned <- TRUE
      .LAST_FIT$msg <- m
      OTHER_WARNINGS[[length(OTHER_WARNINGS) + 1L]] <<- tibble::tibble(
        where = label, message = m)
    }
    invokeRestart("muffleWarning")
  })
}

.vif_ox <- function(fit, label) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2L) return(list(vif = NA_real_, measure = "undefined (1 predictor)"))
  v <- withCallingHandlers(
    tryCatch(car::vif(fit), error = function(e) NULL),
    warning = function(w) {
      if (grepl("No intercept", conditionMessage(w), fixed = TRUE)) {
        .VIF_STATE$no_intercept_warnings <-
          .VIF_STATE$no_intercept_warnings + 1L
      } else {
        OTHER_WARNINGS[[length(OTHER_WARNINGS) + 1L]] <<- tibble::tibble(
          where = paste0("car::vif / ", label), message = conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    })
  if (is.null(v)) return(list(vif = NA_real_, measure = "car::vif failed"))
  if (is.matrix(v))
    list(vif = unname(v["OX", 1]),
         measure = "GVIF (generalised; OX has Df 1 so GVIF == VIF)")
  else list(vif = unname(v[["OX"]]), measure = "VIF")
}

.fit_rungs <- function(d, label, spine) {
  d <- .drop_empty(d, label, spine)
  om <- setdiff(SPINE_CANDIDATES, spine)
  dplyr::bind_rows(lapply(RUNGS, function(rg) {
    rhs <- c(RUNG_RHS[[rg]], spine)
    f <- stats::reformulate(rhs,
           response = "survival::Surv(bcss_time_years, bcss_event)")
    fit <- .catch(survival::coxph(f, data = d),
                  paste0("coxph / ", label, " / ", rg))
    cw <- .LAST_FIT$warned; cm <- .LAST_FIT$msg
    co <- summary(fit)$coefficients
    .stop_if("OX" %in% rownames(co), label, " ", rg, ": no OX coefficient")
    est <- co["OX", 1]; se <- co["OX", "se(coef)"]; p <- co["OX", ncol(co)]
    vf  <- .vif_ox(fit, paste(label, rg))
    tibble::tibble(set = label, rung = rg, n = fit$n, events = fit$nevent,
      estimate = est, se = se, HR = exp(est),
      ci_lo = est - CRIT * se, ci_hi = est + CRIT * se,
      hr_lo = exp(est - CRIT * se), hr_hi = exp(est + CRIT * se),
      excl_0 = (est - CRIT * se) > 0 | (est + CRIT * se) < 0,
      p = p, vif_ox = vf$vif, vif_measure = vf$measure,
      vif_above_caveat = isTRUE(vf$vif > VIF_CAVEAT),
      converge_warning = cw, converge_message = cm,
      terms = paste(attr(stats::terms(fit), "term.labels"), collapse = " + "),
      omitted = paste(om, collapse = ", "))
  }))
}

.show <- function(tb) .pr(tb %>% dplyr::transmute(set, rung, n, events,
    logHR = round(estimate, 4), HR = round(HR, 4),
    ci_HR = paste0("[", sprintf("%.4f", hr_lo), ", ",
                   sprintf("%.4f", hr_hi), "]"),
    excl_0, p = signif(p, 3), vif = round(vif_ox, 3),
    conv = ifelse(converge_warning, "WARNED", "")))

# =============================================================================
# 5. C. m0 AND m1 per PAM50 subtype, cause-specific survival
# =============================================================================
message("\n5. C. m0 and m1 per PAM50 subtype (ER-positive IHC, PRIMARY split)")
message("   m0: ~ OX + treatment      m1: + PROLIF")
message("   ", .sub(paste0(
  "BOTH rungs, because note section 3 item 3 asks whether the proliferation ",
  "adjustment DOES WORK inside a subtype: one where m0 is null and m1 is ",
  "protective would be a better lead than one where both are protective. ",
  "m2 and m3 are NOT fitted.")))

BY_SUBTYPE <- dplyr::bind_rows(lapply(FIT_SUBTYPES, function(st) {
  lab <- paste0("ER-pos (IHC) / ", st)
  .fit_rungs(DP[as.character(DP$subtype_term_x01) == st, ], lab,
             .spine_for(lab))
}))
.show(BY_SUBTYPE)

# --- the ER.Expr repartition: INTERNAL CHECK, NOT REPLICATION --------------
message("\n   C2. THE SAME BLOCK UNDER ER.Expr")
message("   ", .sub(paste0(
  "*** INTERNAL CHECK, NOT REPLICATION. *** This is THE SAME COHORT ",
  "REPARTITIONED, not a second cohort. The two ER calls disagree on 127 of ",
  "1,940 callable patients, so this is a cheap internal consistency check and ",
  "MUST NEVER be described as replication, validation or confirmation (note ",
  "section 3 item 2).")))
COUNTS_S <- .counts(DS, "ER-pos (ER.Expr)")
.pr(COUNTS_S %>% dplyr::mutate(median_fu_y = round(median_fu_y, 2),
                               epv = round(epv, 1)))
FIT_S <- COUNTS_S$subtype[COUNTS_S$fit]
AUDIT_S <- dplyr::bind_rows(lapply(FIT_S, function(st)
  .audit_one(DS[as.character(DS$subtype_term_x01) == st, ],
             paste0("ER-pos (ER.Expr) / ", st))))
.spine_for_s <- function(label) {
  a <- AUDIT_S[AUDIT_S$set == label & AUDIT_S$term %in% SPINE_CANDIDATES, ]
  a$term[!a$constant]
}
BY_SUBTYPE_S <- dplyr::bind_rows(lapply(FIT_S, function(st) {
  lab <- paste0("ER-pos (ER.Expr) / ", st)
  .fit_rungs(DS[as.character(DS$subtype_term_x01) == st, ], lab,
             .spine_for_s(lab))
}))
.show(BY_SUBTYPE_S)
if (!setequal(FIT_SUBTYPES, FIT_S)) {
  message("   NOTE: the two partitions do not make the same subtypes ",
          "fittable. IHC: ", paste(sort(FIT_SUBTYPES), collapse = ", "),
          " | ER.Expr: ", paste(sort(FIT_S), collapse = ", "))
}

# =============================================================================
# 6. D. THE TREATMENT SURVEY - cells first, and it may not be possible
# =============================================================================
# Note section 4. Cross-tabulate PAM50 by treatment and report cell sizes AND
# events BEFORE fitting anything. Reporting that the survey cannot be done is a
# legitimate output.
message("\n6. D. the treatment survey (note section 4). Cells FIRST.")

message("\n   PAM50 x treatment, cell sizes (ER-positive IHC):")
print(table(PAM50 = DP$subtype_term_x01, treatment = DP$treatment))
message("\n   PAM50 x treatment, EVENTS:")
print(stats::xtabs(bcss_event ~ subtype_term_x01 + treatment, data = DP))

TREAT_CELLS <- DP %>%
  dplyr::group_by(treatment = as.character(treatment)) %>%
  dplyr::summarise(n = dplyr::n(), events = sum(bcss_event),
                   pam50_levels = .n_lev(subtype_term_x01),
                   median_fu_y = stats::median(bcss_time_years),
                   .groups = "drop") %>%
  # m1 within a treatment group keeps PAM50 in the spine, so the parameter
  # count is OX + PROLIF + (PAM50 levels - 1).
  dplyr::mutate(params_m1 = 2L + ifelse(pam50_levels <= 1L, 0L,
                                        pam50_levels - 1L),
                epv = events / params_m1,
                fit = events >= MIN_EVENTS & epv >= EPV_FLOOR,
                epv_flag = epv < EPV_FLAG,
                has_chemo = grepl("CT", treatment),
                has_endocrine = grepl("HT", treatment)) %>%
  dplyr::arrange(dplyr::desc(events))
message("\n   treatment groups within the ER-positive stratum:")
.pr(TREAT_CELLS %>% dplyr::mutate(median_fu_y = round(median_fu_y, 2),
                                  epv = round(epv, 1)))

CHEMO <- TREAT_CELLS %>% dplyr::filter(has_chemo)
message("\n   ", .sub(paste0(
  "THE SURVEY'S OWN CONSTRAINT, reported before any fit. Note section 4's ",
  "question needs a chemotherapy arm to contrast with an endocrine one. ",
  "Of the ", nrow(CHEMO), " chemotherapy-containing groups here, carrying ",
  sum(CHEMO$events), " events between them, ", sum(CHEMO$fit),
  " clear the threshold. The endocrine-only groups carry ",
  sum(TREAT_CELLS$events[TREAT_CELLS$has_endocrine & !TREAT_CELLS$has_chemo]),
  " events. ANY CONTRAST BETWEEN THEM IS THEREFORE LOPSIDED, and that is a ",
  "property of METABRIC's treatment era, not of the exposure.")))

TREAT_FIT <- TREAT_CELLS$treatment[TREAT_CELLS$fit]
if (length(TREAT_FIT) == 0L) {
  BY_TREATMENT <- tibble::tibble()
  TREAT_NOTE <- paste0(
    "THE TREATMENT SURVEY CANNOT BE DONE. No treatment group in the ",
    "ER-positive stratum clears events >= ", MIN_EVENTS, " with at least ",
    EPV_FLOOR, " events per parameter once PAM50 is retained in the spine. ",
    "NOTHING WAS FITTED. Reporting that the survey cannot be done is the ",
    "output.")
  message("\n   ", .sub(TREAT_NOTE))
} else {
  message("\n   fitting m1 within ", length(TREAT_FIT),
          " treatment group(s), PAM50 RETAINED in the spine:")
  .treat_spine <- function(d) {
    sp <- character(0)
    if (.n_lev(d$subtype_term_x01) > 1L) sp <- c(sp, "subtype_term_x01")
    sp
  }
  BY_TREATMENT <- dplyr::bind_rows(lapply(TREAT_FIT, function(tg) {
    d   <- DP[as.character(DP$treatment) == tg, ]
    lab <- paste0("ER-pos (IHC) / treatment ", tg)
    sp  <- .treat_spine(d)
    d   <- .drop_empty(d, lab, sp)
    rhs <- c(RUNG_RHS$m1, sp)
    f <- stats::reformulate(rhs,
           response = "survival::Surv(bcss_time_years, bcss_event)")
    fit <- .catch(survival::coxph(f, data = d), paste0("coxph / ", lab))
    cw <- .LAST_FIT$warned; cm <- .LAST_FIT$msg
    co <- summary(fit)$coefficients
    est <- co["OX", 1]; se <- co["OX", "se(coef)"]; p <- co["OX", ncol(co)]
    vf  <- .vif_ox(fit, lab)
    tibble::tibble(treatment = tg, rung = "m1", n = fit$n,
      events = fit$nevent, has_chemo = grepl("CT", tg),
      estimate = est, se = se, HR = exp(est),
      ci_lo = est - CRIT*se, ci_hi = est + CRIT*se,
      hr_lo = exp(est - CRIT*se), hr_hi = exp(est + CRIT*se),
      excl_0 = (est - CRIT*se) > 0 | (est + CRIT*se) < 0,
      p = p, vif_ox = vf$vif,
      converge_warning = cw, converge_message = cm,
      terms = paste(attr(stats::terms(fit), "term.labels"), collapse = " + "))
  }))
  .pr(BY_TREATMENT %>% dplyr::transmute(treatment, n, events, has_chemo,
      logHR = round(estimate, 4), HR = round(HR, 4),
      ci_HR = paste0("[", sprintf("%.4f", hr_lo), ", ",
                     sprintf("%.4f", hr_hi), "]"),
      excl_0, p = signif(p, 3), vif = round(vif_ox, 3),
    conv = ifelse(converge_warning, "WARNED", "")))
  WARNED <- BY_TREATMENT[BY_TREATMENT$converge_warning, ]
  if (nrow(WARNED)) {
    message("   ", .sub(paste0(
      "A 'WARNED' ROW IS NOT A CLEAN FIT. coxph reported that the ",
      "log-likelihood converged before a spine coefficient did, which happens ",
      "when a PAM50 level inside that treatment group has observations but ",
      "ZERO EVENTS, giving it an infinite coefficient. The OX coefficient is ",
      "still estimable but the spine is partly unidentified. Affected: ",
      paste(WARNED$treatment, collapse = ", "),
      ". THE LEVEL WAS NOT DROPPED: unlike an EMPTY level, a zero-event level ",
      "still contributes its patients to other levels' risk sets, so removing ",
      "them would change the data rather than tidy it.")))
  }
  TREAT_NOTE <- paste0(
    "Fitted in ", length(TREAT_FIT), " of ", nrow(TREAT_CELLS),
    " treatment groups. The groups that could not be fitted are ",
    paste(setdiff(TREAT_CELLS$treatment, TREAT_FIT), collapse = ", "),
    ". These are SURVEY cells, not a test of effect modification - no ",
    "interaction was fitted and the event counts would not support one.")
  message("   ", .sub(TREAT_NOTE))
}

# =============================================================================
# 7. E. WINNER'S CURSE - the scatter, shown beside the leads
# =============================================================================
# Note section 2 caution 1. NOTHING IS TESTED HERE. The point is that the
# scatter is visible next to whichever subtype looks largest.
message("\n7. E. winner's curse: the scatter against E41's pooled estimate")

SCATTER <- BY_SUBTYPE %>%
  dplyr::filter(rung == "m1") %>%
  dplyr::transmute(set, n, events,
                   logHR = estimate, HR,
                   hr_lo, hr_hi,
                   pooled_logHR = POOLED$estimate, pooled_HR = POOLED$HR,
                   distance_from_pooled = estimate - POOLED$estimate,
                   ci_contains_pooled = hr_lo <= POOLED$HR &
                                        hr_hi >= POOLED$HR) %>%
  dplyr::arrange(logHR)
.pr(SCATTER %>% dplyr::transmute(set, n, events,
    HR = round(HR, 4),
    ci_HR = paste0("[", sprintf("%.4f", hr_lo), ", ",
                   sprintf("%.4f", hr_hi), "]"),
    pooled_HR = round(pooled_HR, 4),
    distance_logHR = round(distance_from_pooled, 4),
    ci_contains_pooled))

SPREAD <- tibble::tibble(
  n_subtypes = nrow(SCATTER),
  pooled_HR = POOLED$HR,
  lowest_HR = min(SCATTER$HR), highest_HR = max(SCATTER$HR),
  range_logHR = max(SCATTER$logHR) - min(SCATTER$logHR),
  furthest_below = min(SCATTER$distance_from_pooled),
  furthest_above = max(SCATTER$distance_from_pooled),
  all_cis_contain_pooled = all(SCATTER$ci_contains_pooled),
  n_cis_containing = sum(SCATTER$ci_contains_pooled))
.pr(SPREAD %>% dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                           ~ round(.x, 4))))
message("   ", .sub(paste0(
  "CAUTION 1, restated with these numbers in front of it. The ",
  nrow(SCATTER), " subtype estimates span ", round(SPREAD$range_logHR, 3),
  " on the log-hazard scale, and ", SPREAD$n_cis_containing, " of ",
  nrow(SCATTER), " intervals contain the pooled value. ",
  "SCATTER OF THIS KIND IS EXPECTED FROM SAMPLING VARIATION ALONE. The ",
  "largest estimate is partly the maximum of several noisy draws, so NO ",
  "SUBTYPE IS A LEAD MERELY FOR BEING THE LARGEST. Nothing here is tested.")))

# =============================================================================
# 8. Warnings
# =============================================================================
message("\n8. warnings")
message("   ", .sub(paste0(
  "car::vif raised the spurious 'No intercept: vifs may not be sensible' ",
  .VIF_STATE$no_intercept_warnings, " time(s), once per call, and it was ",
  "counted and muffled deliberately. A Cox model has NO INTERCEPT BY DESIGN ",
  "and the VIF comes from the correlation matrix of the coefficient ",
  "estimates, which needs none. This is the same artefact identified in E40 ",
  "and E41.")))
OTHER_TAB <- if (length(OTHER_WARNINGS)) dplyr::bind_rows(OTHER_WARNINGS) else
  tibble::tibble(where = character(0), message = character(0))
if (nrow(OTHER_TAB)) {
  message("   OTHER warnings, reported and NOT muffled away:")
  .pr(OTHER_TAB %>% dplyr::count(where, message, name = "n"))
} else {
  message("   no other warning was raised.")
}
DROPPED_TAB <- if (length(DROPPED_LEVELS)) dplyr::bind_rows(DROPPED_LEVELS) else
  tibble::tibble(set = character(0), term = character(0),
                 dropped = character(0))
if (nrow(DROPPED_TAB)) {
  message("   empty factor levels dropped EXPLICITLY (0 observations each, so ",
          "no estimate changes):")
  .pr(DROPPED_TAB)
}

# =============================================================================
# 9. F. Save - every table prefixed x01_ so it sorts away from the E-outputs
# =============================================================================
message("\n9. F. saving (all tables prefixed x01_)")

readr::write_csv(COUNTS,       file.path(DIR_TABLES, "x01_counts_pam50.csv"))
readr::write_csv(COUNTS_S,     file.path(DIR_TABLES, "x01_counts_pam50_erexpr.csv"))
readr::write_csv(DETECTABLE,   file.path(DIR_TABLES, "x01_detectable.csv"))
readr::write_csv(AUDIT,        file.path(DIR_TABLES, "x01_spine_audit.csv"))
readr::write_csv(BY_SUBTYPE,   file.path(DIR_TABLES, "x01_by_pam50.csv"))
readr::write_csv(BY_SUBTYPE_S, file.path(DIR_TABLES, "x01_by_pam50_erexpr.csv"))
readr::write_csv(TREAT_CELLS,  file.path(DIR_TABLES, "x01_treatment_cells.csv"))
if (nrow(BY_TREATMENT))
  readr::write_csv(BY_TREATMENT, file.path(DIR_TABLES, "x01_by_treatment.csv"))
readr::write_csv(SCATTER,      file.path(DIR_TABLES, "x01_scatter.csv"))

saveRDS(list(
  counts            = COUNTS,
  counts_erexpr     = COUNTS_S,
  detectable        = DETECTABLE,
  not_fitted        = NOT_FIT,
  flagged           = FLAGGED,
  audit             = AUDIT,
  audit_erexpr      = AUDIT_S,
  by_pam50          = BY_SUBTYPE,
  by_pam50_erexpr   = BY_SUBTYPE_S,
  treatment_cells   = TREAT_CELLS,
  by_treatment      = BY_TREATMENT,
  treatment_note    = TREAT_NOTE,
  scatter           = SCATTER,
  spread            = SPREAD,
  pooled_reference  = POOLED,
  dropped_levels    = DROPPED_TAB,
  other_warnings    = OTHER_TAB,
  vif_no_intercept_warnings = .VIF_STATE$no_intercept_warnings,
  spec = list(
    note = paste0("docs/2026-10-03_X01_exploratory_leads_note.md, committed ",
                  "at c42e3b0"),
    posture = POSTURE,
    not_an_e_number = paste0(
      "The E-numbered line is manuscript-bound: every E-script's result is ",
      "reported whatever it shows, under a declaration committed before it ",
      "ran. X01 is NOT. Keeping it off the numbering makes the difference ",
      "visible in the file listing so nobody downstream mistakes it for part ",
      "of the declared work."),
    # Note section 5, VERBATIM.
    cannot_license_1 = "No manuscript sentence, figure, legend or supplementary table.",
    cannot_license_2 = paste0(
      "No claim that any subtype or treatment stratum 'shows' or ",
      "'demonstrates' anything. The vocabulary here is suggests, is ",
      "consistent with, would be worth testing."),
    cannot_license_3 = paste0(
      "No comparison of any coefficient across cohorts as a raw value (E39 ",
      "section 4.2 still holds)."),
    cannot_license_4 = "No reinstatement of anything on the retired list.",
    cannot_license_5 = paste0(
      "No interaction test, no formal effect-modification test. Not because ",
      "of posture but because the event counts do not support one."),
    e42_prohibition_not_relaxed = paste0(
      "E42 section 13.2 forbids any subtype result becoming a finding in its ",
      "own right, and that prohibition is NOT relaxed by moving the work into ",
      "this note."),
    caution_1 = paste0(
      "WINNER'S CURSE. Individual PAM50 subtypes scatter by sampling ",
      "variation alone; the largest stratum estimate is partly the maximum of ",
      "several noisy draws. NO SUBTYPE IS THE LEAD MERELY FOR BEING THE ",
      "LARGEST. See $scatter and $spread."),
    caution_2 = paste0(
      "THE LUMINAL pCR ESTIMATE IS PROBABLY INFLATED. E42 declared the ",
      "detectable effect for HRpos_HER2neg at 0.312 log-odds; the observed ",
      "0.276 is SMALLER than that, and an underpowered analysis that clears ",
      "significance typically overstates magnitude."),
    caution_3 = paste0(
      "THE s508 EXPLANATION ALREADY FAILED ONCE. E42 measured ",
      "rho(OX, subtype2) = -0.054 and refuted it. Treat every mechanism here ",
      "the same way."),
    er_expr_is_not_replication = paste0(
      "The ER.Expr block is an INTERNAL CHECK on THE SAME COHORT ",
      "REPARTITIONED - the two calls disagree on 127 of 1,940 callable ",
      "patients. It must NEVER be described as replication, validation or ",
      "confirmation (note section 3 item 2)."),
    threshold = paste0(
      "E42's threshold, mapped to a Cox model with the mapping stated. E42 ",
      "used the MINORITY CLASS because a logistic fit is limited by the rarer ",
      "outcome; A COX MODEL HAS NO MINORITY CLASS and its limiting count is ",
      "EVENTS. So: events >= ", MIN_EVENTS, " and events per parameter >= ",
      EPV_FLOOR, ", flagged below ", EPV_FLAG,
      ", parameters counted on m1. The detectable-effect figures are ",
      "Schoenfeld approximations that IGNORE covariate adjustment, so the ",
      "true detectable effects are LARGER. A null in a low-event subtype is ",
      "UNINFORMATIVE and is NOT absence."),
    rungs = paste0(
      "m0 (~ OX + treatment) and m1 (+ PROLIF) only. BOTH, because note ",
      "section 3 item 3 asks whether the proliferation adjustment does work ",
      "inside a subtype. m2 and m3 are NOT fitted, and no interaction of any ",
      "kind is fitted."),
    no_rescoring = paste0(
      "E41's $frame is read and stratified. NOTHING was rescored, ",
      "METABRIC_DATA.RData was NOT reloaded, and nothing already in ",
      "e41_metabric_ladder.rds was refitted."),
    instrument = paste0(
      "SINGLE-INSTRUMENT. METABRIC is an array and mitoPPS is unavailable. ",
      "No coefficient is compared across cohorts as a raw value."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  built = Sys.time()), PATH_X01)

message("   ", PATH_X01)
message("   ", length(list.files(DIR_TABLES, pattern = "^x01_")),
        " x01_ tables in ", DIR_TABLES)

message("\n", strrep("=", 78))
message("X01 done. Numbers only; nothing here is interpreted.")
message("NOT PART OF THE DECLARED E-LINE. MANUSCRIPT-EXCLUDED.")
message("  ", .sub(POSTURE))
message(strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_X01)

  # READ THIS FIRST, every time.
  cat(x$spec$posture, "\n")
  cat(x$spec$e42_prohibition_not_relaxed, "\n")

  # A. the counts, which decide what could be asked at all.
  x$counts %>% dplyr::mutate(epv = round(epv, 1)) %>% as.data.frame()
  x$not_fitted %>% as.data.frame()
  x$flagged %>% as.data.frame()
  x$detectable %>% as.data.frame()

  # B. the spine audit.
  x$audit %>% as.data.frame()
  x$dropped_levels %>% as.data.frame()

  # C. m0 against m1 per subtype. Note section 3 item 3 is the m0-to-m1 step.
  x$by_pam50 %>% dplyr::transmute(set, rung, n, events, HR = round(HR, 4),
      lo = round(hr_lo, 4), hi = round(hr_hi, 4), excl_0, p = signif(p, 3),
      vif = round(vif_ox, 3)) %>% as.data.frame()

  # C2. the SAME COHORT REPARTITIONED. NOT replication.
  cat(x$spec$er_expr_is_not_replication, "\n")
  x$by_pam50_erexpr %>% dplyr::transmute(set, rung, n, events,
      HR = round(HR, 4), lo = round(hr_lo, 4), hi = round(hr_hi, 4),
      excl_0, p = signif(p, 3)) %>% as.data.frame()

  # D. the treatment survey, cells first.
  x$treatment_cells %>% dplyr::mutate(epv = round(epv, 1)) %>% as.data.frame()
  x$by_treatment %>% as.data.frame()
  cat(x$treatment_note, "\n")

  # E. the scatter, read BEFORE any subtype is called a lead.
  cat(x$spec$caution_1, "\n")
  x$scatter %>% as.data.frame()
  x$spread %>% as.data.frame()

  # warnings
  x$vif_no_intercept_warnings
  x$other_warnings %>% as.data.frame()

  utils::str(x$spec)
}
