# E42_pcr_subtype_scope.R
# =============================================================================
# A SCOPE CHECK ON E40's PRIMARY FINDING. IT MAKES NO DECISIONS.
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md,
# section 13, committed at dba8ee4. Read 13.1 to 13.7 first.
#
# =============================================================================
# WHAT THIS IS, AND THE FOUR THINGS IT CANNOT DO (declaration 13.2)
# =============================================================================
# E40's three-cohort m1 is the arm's primary finding: pooled log-odds
# -0.162 [-0.281, -0.042], I2 = 0.0%, signs agreeing in all three cohorts. That
# estimate is pooled WITHIN each cohort ACROSS subtypes, and the model's
# assumption of one common OX effect inside every subtype has never been tested.
# E42 BOUNDS THE SCOPE of a claim already being made. That is not subgroup
# searching, and the check cannot change the claim's status:
#
#   1. IT CANNOT PROMOTE OR DEMOTE THE PRIMARY FINDING. E40's three-cohort m1
#      stands as the primary result whatever the strata show.
#   2. NO SUBTYPE RESULT BECOMES A FINDING IN ITS OWN RIGHT. A strong stratum is
#      a stratum of a pooled estimate, never a separate discovery.
#   3. NO INTERACTION TERM AND NO FORMAL EFFECT-MODIFICATION TEST. Declared out,
#      as in section 7, at these event counts.
#   4. IT DOES NOT TOUCH METABRIC. The PAM50 split there is X01, a separate
#      activity with a separate posture, and is not part of this line.
#
# =============================================================================
# THE READING RULE, FIXED BEFORE THE FIT (declaration 13.6)
# =============================================================================
#   consistent across subtypes, CIs overlapping
#       -> the draft sentence is UNCHANGED; the scope statement is reported
#   carried by one subtype, others null but underpowered
#       -> UNCHANGED, with the stratum estimates beside it and the power
#          limitation stated
#   opposite signs with NON-OVERLAPPING CIs
#       -> the claim is QUALIFIED to the subtype where it holds
#
# THE THIRD ROW IS THE ONLY ONE THAT CHANGES THE CLAIM, and it requires
# non-overlapping intervals, not one stratum being significant and another not.
# DIFFERING p VALUES ACROSS STRATA OF DIFFERING SIZE ARE NOT EVIDENCE OF A
# DIFFERING EFFECT. This script does NOT classify the outcome.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits m1 ONLY. NO m0, NO m2, NO m3.
#   - NO interaction term of any kind.
#   - It does NOT rescore anything. E39's frames and scores are reused
#     unchanged, so the exposure is byte-identical to E40's.
#   - It does NOT refit any DRFS model (13.7).
#   - It does NOT touch METABRIC.
#   - NO p value is described as a near miss. NO per-gene FDR.
#   - It does NOT interpret. It prints numbers, never verdicts.
#   - It writes NOTHING to myc_human_validation, FROZEN at d3ac60e, or to
#     myc_mouse. h4_outcome_models.rds is read READ-ONLY.
#
# SINGLE-INSTRUMENT, as the whole line is. mitoPPS does not exist in these
# cohorts. No coefficient is compared across cohorts as a raw value (4.2).
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
})

message("\nE42: does the pCR finding hold across subtypes?\n", strrep("=", 78))
message("A SCOPE CHECK. It cannot promote or demote E40's primary finding.")
message("It fits m1 ONLY, within subtype, and classifies nothing.")

PATH_E39 <- file.path(DIR_RESULTS, "e39_prepare_and_gate.rds")
# READ-ONLY, from the FROZEN repo. Nothing is ever written there.
PATH_H4  <- paste0("/Users/gs/code/myc_human_validation/results/",
                   "h4_outcome_models.rds")
PATH_E42 <- file.path(DIR_RESULTS, "e42_pcr_subtype_scope.rds")

# =============================================================================
# 0. CONSTANTS - from the declaration, none invented here
# =============================================================================
CI_LEVEL <- 0.95
CRIT     <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

COHORTS <- c("GSE25066", "GSE194040", "GSE164458")
PRIMARY <- "GSE25066"

# Declaration 13.3. ONE rung, and the subtype term is gone by construction.
RUNG      <- "m1"
RUNG_RHS  <- c("OX", "PROLIF")          # plus treatment where it varies
SPINE_CANDIDATES <- c("subtype_term", "treatment")

# E40's primary finding, printed beside every stratum estimate (13.3, step E).
E40_POOLED       <- -0.1616
E40_POOLED_CI    <- c(-0.2811, -0.0422)
E40_POOLED_I2    <- 0.0
E40_TNBC_164458  <- -0.2053
E40_TNBC_CI      <- c(-0.4228, 0.0122)

# Declaration 13.4: the GSE25066 partition identity, on the 465.
N_PRIMARY   <- 465L
PART_HRPOS  <- 278L
PART_TNBC   <- 187L

# Declaration 13.5, the counts declared BEFORE the fit. Each is a STOP.
DECLARED <- tibble::tribble(
  ~cohort,      ~subtype,         ~n,    ~events,
  "GSE194040",  "HRneg_HER2pos",   89L,   56L,
  "GSE194040",  "HRpos_HER2neg",  379L,   64L,
  "GSE25066",   "HRpos_HER2neg",  278L,   30L,
  "GSE194040",  "HRpos_HER2pos",  156L,   57L,
  "GSE164458",  "TNBC",           482L,  236L,
  "GSE194040",  "TNBC",           364L,  142L,
  "GSE25066",   "TNBC",           187L,   63L)

# Declaration 13.5, the fitting threshold, NAMED THERE rather than chosen here.
# Events per parameter is computed on the MINORITY class, because a logistic
# fit is limited by whichever outcome is rarer - HRneg_HER2pos has 56 events
# and only 33 non-events.
MIN_CLASS_FLOOR <- 5L     # minority class must reach this
EPV_FLOOR       <- 5      # events per parameter must reach this to fit
EPV_FLAG        <- 10     # below this the stratum is FLAGGED in the output

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr  <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.n_lev <- function(v) dplyr::n_distinct(v[!is.na(v)])

# =============================================================================
# 1. Inputs - E39's frames and scores, REUSED UNCHANGED
# =============================================================================
message("\n1. inputs (E39's frames and scores, reused unchanged - NO rescoring)")

.stop_if(file.exists(PATH_E39), "absent: ", PATH_E39,
         ". Source scripts/E39_prepare_and_gate.R first.")
.stop_if(file.exists(PATH_H4), "absent: ", PATH_H4,
         ". myc_human_validation must be present and is read-only.")
e39 <- readRDS(PATH_E39)
h4  <- readRDS(PATH_H4)

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

# GSE25066 uses the collapsed subtype2 on the 465; the other two use subtype.
# Identical to E40's construction, so the strata are E40's strata.
for (cn in COHORTS) {
  FR[[cn]]$subtype_term <- if (cn == PRIMARY) FR[[cn]]$subtype2 else
    factor(FR[[cn]]$subtype)
}
MODEL <- stats::setNames(lapply(COHORTS, function(cn) {
  d <- FR[[cn]]
  if (cn == PRIMARY) d <- d[d$in_model, ]
  d[!is.na(d$subtype_term), ]
}), COHORTS)
.stop_if(nrow(MODEL[[PRIMARY]]) == N_PRIMARY,
         PRIMARY, "'s model set is ", nrow(MODEL[[PRIMARY]]), ", not ",
         N_PRIMARY)
message("   model sets: ",
        paste(sprintf("%s %d", COHORTS, vapply(MODEL, nrow, integer(1))),
              collapse = " | "))
message("   OX and MYC are E40's, from script 13's .z(); PROLIF is E39's. ",
        "NOTHING is rescored.")

# --- declaration 13.4, asserted again rather than assumed -------------------
dm <- MODEL[[PRIMARY]]
part_tab <- table(subtype2 = dm$subtype_term, er_primary = dm$er_primary)
.stop_if(all((dm$subtype_term == "HRpos_HER2neg") == (dm$er_primary == "P")),
         PRIMARY, ": subtype2 and er_primary do NOT partition identically. ",
         "Declaration 13.4 and 5.2 both say they do. STOP.")
.stop_if(sum(dm$subtype_term == "HRpos_HER2neg") == PART_HRPOS &&
         sum(dm$subtype_term == "TNBC") == PART_TNBC,
         PRIMARY, ": the partition is not the declared ", PART_HRPOS, " / ",
         PART_TNBC)
message("\n   DECLARATION 13.4, re-asserted:")
.pr(part_tab)
message("   ", .sub(paste0(
  "In ", PRIMARY, " subtype2 and er_primary partition the 465 IDENTICALLY, ",
  PART_HRPOS, " / ", PART_TNBC, " on both. So this cohort's subtype split and ",
  "its ER split ARE THE SAME SPLIT, and the TNBC stratum here is also its ",
  "ER-negative stratum. Asserted in E40, on the 465 - not on the 470, where ",
  "subtype still has four levels and the identity does not hold.")))

# =============================================================================
# 2. A. COUNTS FIRST - the stop-and-check, before any fit
# =============================================================================
# Declaration 13.5. Reported in full. The fits do not matter until this is read.
message("\n2. A. COUNTS FIRST (13.5). Read this block before any coefficient.")

COUNTS <- dplyr::bind_rows(lapply(COHORTS, function(cn) {
  MODEL[[cn]] %>%
    dplyr::group_by(subtype = as.character(subtype_term)) %>%
    dplyr::summarise(cohort = cn, n = dplyr::n(),
                     events = sum(pcr == 1), non_events = sum(pcr == 0),
                     rate = mean(pcr == 1),
                     treat_levels = .n_lev(treatment), .groups = "drop")
})) %>%
  dplyr::mutate(
    # treatment is omitted where constant (13.3, 5.2), so it costs no parameter
    params    = 2L + ifelse(treat_levels <= 1L, 0L, treat_levels - 1L),
    min_class = pmin(events, non_events),
    epv       = min_class / params,
    fit       = min_class >= MIN_CLASS_FLOOR & epv >= EPV_FLOOR,
    epv_flag  = epv < EPV_FLAG) %>%
  dplyr::select(cohort, subtype, n, events, non_events, rate, treat_levels,
                params, min_class, epv, fit, epv_flag) %>%
  dplyr::arrange(subtype, cohort)

.pr(COUNTS %>% dplyr::mutate(rate = round(rate, 3), epv = round(epv, 1)))

# Every count is a STOP against the declaration, not a warning.
chk <- DECLARED %>%
  dplyr::left_join(COUNTS[c("cohort", "subtype", "n", "events")],
                   by = c("cohort", "subtype"), suffix = c("_dec", "_got"))
.stop_if(nrow(chk) == nrow(DECLARED) && !anyNA(chk$n_got),
         "a declared stratum is missing from the data")
.stop_if(all(chk$n_dec == chk$n_got) && all(chk$events_dec == chk$events_got),
         "the strata do not match declaration 13.5's table. Got: ",
         paste(sprintf("%s/%s n=%d e=%d", chk$cohort, chk$subtype,
                       chk$n_got, chk$events_got), collapse = "; "))
.stop_if(nrow(COUNTS) == nrow(DECLARED),
         "there are ", nrow(COUNTS), " strata, not the declared ",
         nrow(DECLARED))
message("   ASSERTED: all ", nrow(DECLARED),
        " strata match declaration 13.5's declared n and events")
message("   ", .sub(paste0(
  "Events per parameter is computed on the MINORITY class, not on pCR events, ",
  "because a logistic fit is limited by whichever outcome is rarer: ",
  "HRneg_HER2pos has 56 events and only 33 non-events (13.5).")))

# --- which strata can be pooled at all --------------------------------------
BY_SUBTYPE <- COUNTS %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(n_cohorts = dplyr::n(), n = sum(n), events = sum(events),
                   min_class = sum(min_class),
                   poolable = dplyr::n() > 1L, .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_cohorts))
message("\n   pooled by subtype - and TWO subtypes cannot be pooled at all:")
.pr(BY_SUBTYPE)
message("   ", .sub(paste0(
  "HRneg_HER2pos and HRpos_HER2pos appear in GSE194040 ALONE (13.3, 13.5), so ",
  "neither can be pooled and neither has a second cohort to agree or disagree ",
  "with it. Their single estimate is reported as SINGLE, never as pooled.")))

# --- the fitting threshold, and the detectable effect -----------------------
NOT_FIT <- COUNTS %>% dplyr::filter(!fit)
message("\n   THE THRESHOLD, named in declaration 13.5 rather than chosen ",
        "here: a stratum is fitted where the minority class is at least ",
        MIN_CLASS_FLOOR, " AND events per parameter is at least ", EPV_FLOOR,
        ".")
message("   Strata below ", EPV_FLAG, " events per parameter are FLAGGED but ",
        "still fitted (Peduzzi's rule of thumb, not treated as a cliff; ",
        "Vittinghoff and McCulloch 2007 find 5 to 9 often acceptable).")
if (nrow(NOT_FIT)) {
  message("   NOT FITTED, below the threshold:")
  .pr(NOT_FIT %>% dplyr::transmute(cohort, subtype, min_class, params,
                                   epv = round(epv, 1)))
} else {
  message("   NOT FITTED: none. Every stratum clears the threshold.")
}
FLAGGED <- COUNTS %>% dplyr::filter(fit, epv_flag)
if (nrow(FLAGGED)) {
  message("   FITTED BUT FLAGGED, under ", EPV_FLAG, " events per parameter:")
  .pr(FLAGGED %>% dplyr::transmute(cohort, subtype, min_class, params,
                                   epv = round(epv, 1)))
}

# Approximate detectable effect per POOLED subtype stratum, on the log-odds
# scale, for a standardised continuous exposure at 80% power, alpha 0.05.
# se(beta) ~ 1 / sqrt(n * p * (1 - p)) for a 1-SD exposure in a logistic model;
# the detectable |beta| is (z_{1-a/2} + z_{power}) * se. This is an
# APPROXIMATION and is labelled as one - it ignores covariate adjustment, which
# inflates the variance, so the true detectable effect is LARGER.
Z_A <- stats::qnorm(0.975); Z_B <- stats::qnorm(0.80)
DETECTABLE <- COUNTS %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(n_cohorts = dplyr::n(), n = sum(n), events = sum(events),
                   .groups = "drop") %>%
  dplyr::mutate(p_event = events / n,
                se_approx = 1 / sqrt(n * p_event * (1 - p_event)),
                detectable_logodds = (Z_A + Z_B) * se_approx,
                detectable_OR = exp(detectable_logodds),
                detectable_logodds_x2var = (Z_A + Z_B) * se_approx * sqrt(2))
message("\n   APPROXIMATE detectable effect per POOLED subtype stratum ",
        "(80% power, alpha 0.05, per 1 within-cohort SD of OX):")
.pr(DETECTABLE %>% dplyr::transmute(subtype, n_cohorts, n, events,
      p_event = round(p_event, 3),
      detectable_logodds = round(detectable_logodds, 3),
      detectable_OR = round(detectable_OR, 3),
      with_2x_variance = round(detectable_logodds_x2var, 3)))
message("   ", .sub(paste0(
  "This is an APPROXIMATION and ignores covariate adjustment, which inflates ",
  "the variance, so the true detectable effect is LARGER - the last column ",
  "shows it at twofold variance inflation. For scale, E40's pooled estimate ",
  "is ", E40_POOLED, ". A NULL IN A LOW-EVENT STRATUM IS REPORTED AS ",
  "UNINFORMATIVE, NOT AS ABSENCE (13.5, section 8).")))

# =============================================================================
# 3. B. SPINE AUDIT per stratum - constants omitted EXPLICITLY
# =============================================================================
# Declaration 5.2 and 13.3. Inside a subtype stratum the subtype term is
# constant BY CONSTRUCTION and is omitted explicitly. treatment is audited too:
# it is constant across the whole of GSE25066 and may become constant inside a
# small stratum of the other two.
message("\n3. B. spine audit, per stratum (5.2, 13.3)")

STRATA <- COUNTS %>% dplyr::filter(fit) %>%
  dplyr::transmute(cohort, subtype,
                   label = paste(cohort, subtype, sep = " / "))

AUDIT <- dplyr::bind_rows(lapply(seq_len(nrow(STRATA)), function(i) {
  cn <- STRATA$cohort[i]; st <- STRATA$subtype[i]
  d  <- MODEL[[cn]]
  d  <- d[as.character(d$subtype_term) == st, ]
  tibble::tibble(cohort = cn, subtype = st, n = nrow(d),
                 term = c("OX", SPINE_CANDIDATES),
                 n_levels = c(.n_lev(d$OX), .n_lev(d$subtype_term),
                              .n_lev(d$treatment))) %>%
    dplyr::mutate(constant = n_levels <= 1L,
                  action = ifelse(constant, "OMITTED (constant)", "fitted"))
}))
.pr(AUDIT)

# The subtype term MUST be constant inside every stratum. That is the
# definition of a stratum, and asserting it guards against a mis-split.
st_rows <- AUDIT[AUDIT$term == "subtype_term", ]
.stop_if(all(st_rows$constant),
         "subtype_term is NOT constant inside ",
         sum(!st_rows$constant), " stratum/strata. The split is wrong. STOP.")
message("   ASSERTED: subtype_term is constant in all ", nrow(st_rows),
        " strata and is OMITTED EXPLICITLY in every one")
tr_const <- AUDIT[AUDIT$term == "treatment" & AUDIT$constant, ]
if (nrow(tr_const)) {
  message("   treatment is CONSTANT and therefore OMITTED EXPLICITLY in:")
  .pr(tr_const %>% dplyr::transmute(cohort, subtype, n, n_levels))
}
message("   ", .sub(paste0(
  "treatment is constant across the WHOLE of ", PRIMARY,
  " (taxane-anthracycline), so it is omitted in both of its strata for the ",
  "same reason E40 omitted it there. It varies in every GSE194040 and ",
  "GSE164458 stratum and is fitted in each.")))

.spine_for <- function(cn, st) {
  a <- AUDIT[AUDIT$cohort == cn & AUDIT$subtype == st &
             AUDIT$term %in% SPINE_CANDIDATES, ]
  a$term[!a$constant]
}

# =============================================================================
# 4. The estimator helpers
# =============================================================================
# VIF. OX is a single continuous term, so its Df is 1 and its GVIF IS the
# ordinary VIF. These are LOGISTIC fits WITH an intercept, so car::vif's
# "No intercept" warning does not arise here - any warning that does appear is
# captured and reported rather than muffled (declaration 13, conventions).
WARNINGS_SEEN <- list()

.vif_ox <- function(fit, label) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2L) return(list(vif = NA_real_, df = 1, err = "",
                           measure = "undefined (single predictor)"))
  err <- ""
  v <- withCallingHandlers(
    tryCatch(car::vif(fit),
             error = function(e) { err <<- conditionMessage(e); NULL }),
    warning = function(w) {
      WARNINGS_SEEN[[length(WARNINGS_SEEN) + 1L]] <<- tibble::tibble(
        where = paste0("car::vif / ", label), message = conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (is.null(v)) return(list(vif = NA_real_, df = NA_real_, err = err,
                              measure = "car::vif failed"))
  if (is.matrix(v))
    list(vif = unname(v["OX", 1]), df = unname(v["OX", 2]), err = "",
         measure = "GVIF (generalised; OX has Df 1 so GVIF == VIF)")
  else list(vif = unname(v[["OX"]]), df = 1, err = "", measure = "VIF")
}

# Separation check. A logistic fit that has separated reports a huge |beta| with
# a huge SE and a fitted probability at 0 or 1. It is DETECTED and REPORTED, not
# silently returned as a coefficient.
.separated <- function(fit) {
  p <- stats::fitted(fit)
  any(p < 1e-8 | p > 1 - 1e-8) ||
    any(abs(stats::coef(fit)[-1]) > 15, na.rm = TRUE)
}

.fit_stratum <- function(cn, st) {
  d  <- MODEL[[cn]]
  d  <- d[as.character(d$subtype_term) == st, ]
  sp <- .spine_for(cn, st)
  om <- setdiff(SPINE_CANDIDATES, sp)
  rhs <- c(RUNG_RHS, sp)
  f   <- stats::reformulate(rhs, response = "pcr")
  fit <- withCallingHandlers(
    stats::glm(f, data = d, family = stats::binomial()),
    warning = function(w) {
      WARNINGS_SEEN[[length(WARNINGS_SEEN) + 1L]] <<- tibble::tibble(
        where = paste0("glm / ", cn, " / ", st), message = conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  co <- summary(fit)$coefficients
  .stop_if("OX" %in% rownames(co), cn, " / ", st, ": no OX coefficient")
  est <- co["OX", 1]; se <- co["OX", 2]; p <- co["OX", 4]
  vf  <- .vif_ox(fit, paste(cn, st, sep = " / "))
  cnt <- COUNTS[COUNTS$cohort == cn & COUNTS$subtype == st, ]
  tibble::tibble(
    cohort = cn, subtype = st, rung = RUNG,
    n = stats::nobs(fit), events = sum(fit$y == 1),
    min_class = cnt$min_class, params = cnt$params, epv = cnt$epv,
    epv_flag = cnt$epv_flag,
    estimate = est, se = se,
    ci_lo = est - CRIT * se, ci_hi = est + CRIT * se,
    excl_0 = (est - CRIT * se) > 0 | (est + CRIT * se) < 0,
    OR = exp(est), or_lo = exp(est - CRIT * se), or_hi = exp(est + CRIT * se),
    p = p,
    vif_ox = vf$vif, vif_measure = vf$measure, vif_error = vf$err,
    separated = .separated(fit),
    terms = paste(attr(stats::terms(fit), "term.labels"), collapse = " + "),
    omitted = paste(om, collapse = ", "))
}

# script 13's .meta(), reproduced verbatim as E40 did, so the pooling is the
# same DerSimonian-Laird estimator this arm has used throughout.
.meta <- function(est, se, label) {
  v <- se^2; w <- 1 / v; k <- length(est)
  if (k < 2L) {
    return(tibble::tibble(subtype = label, k = k,
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
  tibble::tibble(subtype = label, k = k,
    fe_estimate = fe, fe_ci_lo = fe - CRIT*sqrt(vfe),
    fe_ci_hi = fe + CRIT*sqrt(vfe), fe_p = 2*stats::pnorm(-abs(fe/sqrt(vfe))),
    re_estimate = re, re_ci_lo = re - CRIT*sqrt(vre),
    re_ci_hi = re + CRIT*sqrt(vre), re_p = 2*stats::pnorm(-abs(re/sqrt(vre))),
    Q = Q, Q_p = stats::pchisq(Q, df, lower.tail = FALSE),
    tau2 = tau2, I2 = 100*max(0, (Q - df)/Q))
}

# =============================================================================
# 5. C. FIT m1 ONLY, within each subtype stratum of each cohort
# =============================================================================
message("\n5. C. m1 ONLY, logistic, within each subtype stratum")
message("   pcr ~ OX + PROLIF + treatment, with treatment omitted where ",
        "constant and subtype omitted everywhere (constant by construction)")

LADDER <- dplyr::bind_rows(lapply(seq_len(nrow(STRATA)), function(i)
  .fit_stratum(STRATA$cohort[i], STRATA$subtype[i])))

# The CI prints at 4 dp and excl_0 is COMPUTED, not eyeballed: at 3 dp a bound
# of -0.00052 prints as "-0.000", which reads as zero. E41 hit the same hazard.
# This function only PRINTS; every saved tibble carries the exact bounds.
.show <- function(tb) .pr(tb %>% dplyr::transmute(cohort, subtype, n, events,
    logOR = round(estimate, 4), se = round(se, 4),
    ci = paste0("[", sprintf("%.4f", ci_lo), ", ", sprintf("%.4f", ci_hi), "]"),
    excl_0, p = signif(p, 3), vif = round(vif_ox, 3),
    epv = round(epv, 1), flag = ifelse(epv_flag, "<10 EPV", "")))
.show(LADDER %>% dplyr::arrange(subtype, cohort))
message("   VIF measure: ", LADDER$vif_measure[!is.na(LADDER$vif_ox)][1])
.stop_if(!any(LADDER$separated),
         "a stratum fit SEPARATED: ",
         paste(sprintf("%s/%s", LADDER$cohort[LADDER$separated],
                       LADDER$subtype[LADDER$separated]), collapse = ", "),
         ". Its coefficient is not interpretable. STOP.")
message("   ASSERTED: no stratum fit separated")
message("   m1 ONLY. m0, m2 and m3 are NOT fitted here, and NO interaction ",
        "term of any kind is fitted (13.2, 13.3).")

# =============================================================================
# 6. D. META-ANALYSIS - WITHIN subtype, ACROSS cohorts
# =============================================================================
# Declaration 13.3. NOT within cohort across subtypes. Both fixed and random
# are reported always; random is the summary at I2 >= 50% (section 9).
message("\n6. D. meta-analysis WITHIN subtype ACROSS cohorts (13.3)")

SUBTYPES <- sort(unique(LADDER$subtype))
META <- dplyr::bind_rows(lapply(SUBTYPES, function(st) {
  d <- LADDER[LADDER$subtype == st, ]
  m <- .meta(d$estimate, d$se, st)
  m$cohorts      <- paste(d$cohort, collapse = ", ")
  m$signs        <- paste(sprintf("%s %s", d$cohort,
                                  ifelse(d$estimate > 0, "+", "-")),
                          collapse = " ")
  m$signs_agree  <- length(unique(sign(d$estimate))) == 1L
  m$summary_is   <- ifelse(!is.na(m$I2) & m$I2 >= 50, "random effects",
                           "fixed effect")
  m$pooled       <- m$k > 1L
  m
}))
.pr(META %>% dplyr::transmute(subtype, k, pooled,
    fe = round(fe_estimate, 4),
    fe_ci = paste0("[", sprintf("%.3f", fe_ci_lo), ", ",
                   sprintf("%.3f", fe_ci_hi), "]"),
    re = round(re_estimate, 4),
    re_ci = paste0("[", sprintf("%.3f", re_ci_lo), ", ",
                   sprintf("%.3f", re_ci_hi), "]"),
    Q = round(Q, 3), Q_p = signif(Q_p, 3), I2 = round(I2, 1),
    summary_is, signs_agree, signs))
SINGLE <- META[!META$pooled, ]
if (nrow(SINGLE)) {
  message("   ", .sub(paste0(
    "NOT POOLED, one cohort only: ",
    paste(sprintf("%s (%s)", SINGLE$subtype, SINGLE$cohorts), collapse = "; "),
    ". These rows are a SINGLE COHORT ESTIMATE reported in a pooled table and ",
    "must not be read as pooled. k = 1, so Q, I2 and tau2 are undefined and ",
    "the fixed and random columns are the same single estimate (13.3).")))
}

# =============================================================================
# 7. E. REPORT, NOT REFIT, the two strata already answered
# =============================================================================
# Declaration 13.4. And E40's pooled finding printed beside every stratum, so
# each is read against the finding it bounds.
message("\n7. E. what was already known (13.4), reported and NOT refitted")

ALREADY <- tibble::tibble(
  what = c(
    paste0("GSE164458's m1 IS a TNBC-only estimate: ", E40_TNBC_164458,
           " [", E40_TNBC_CI[1], ", ", E40_TNBC_CI[2],
           "]. Its subtype term dropped as CONSTANT in E40 (TNBC-only cohort),",
           " so E40 already reported a TNBC estimate for it."),
    paste0("In GSE25066 subtype2 and er_primary partition the 465 identically",
           " (", PART_HRPOS, " / ", PART_TNBC, ", re-asserted in section 1),",
           " so that cohort's subtype split and its ER split ARE THE SAME",
           " SPLIT.")))
for (i in seq_len(nrow(ALREADY))) message("   - ", .sub(ALREADY$what[i], 70))

AGAINST_E40 <- LADDER %>%
  dplyr::arrange(subtype, cohort) %>%
  dplyr::transmute(cohort, subtype, n, events,
                   stratum_logOR = round(estimate, 4),
                   stratum_ci = paste0("[", sprintf("%.3f", ci_lo), ", ",
                                       sprintf("%.3f", ci_hi), "]"),
                   E40_pooled = E40_POOLED,
                   E40_pooled_ci = paste0("[", E40_POOLED_CI[1], ", ",
                                          E40_POOLED_CI[2], "]"),
                   E40_pooled_I2 = E40_POOLED_I2)
message("\n   every stratum beside the finding it bounds:")
.pr(AGAINST_E40)
message("   ", .sub(paste0(
  "E40's three-cohort m1 is -0.162 [-0.281, -0.042], I2 = 0.0%, signs ",
  "agreeing in all three cohorts. IT STANDS AS THE PRIMARY RESULT WHATEVER ",
  "THESE STRATA SHOW. No row above promotes, demotes, confirms or weakens it, ",
  "and none is a finding in its own right (13.2).")))

# =============================================================================
# 8. F. ANCILLARY - rho(OX, subtype2) in GSE25066
# =============================================================================
# Declaration 13.7. Reported, NOT interpreted, and NO DRFS model is refitted.
message("\n8. F. ancillary: rho(OX, subtype2) in ", PRIMARY, " (13.7)")

st_bin <- as.integer(dm$subtype_term == "TNBC")   # 1 = TNBC, 0 = HRpos_HER2neg
RHO <- tibble::tibble(
  cohort = PRIMARY, n = nrow(dm),
  pair = "OX vs subtype2 (TNBC = 1, HRpos_HER2neg = 0)",
  spearman = stats::cor(dm$OX, st_bin, method = "spearman"),
  mean_OX_HRpos = mean(dm$OX[st_bin == 0]),
  mean_OX_TNBC  = mean(dm$OX[st_bin == 1]))
.pr(RHO %>% dplyr::mutate(spearman = round(spearman, 4),
                          mean_OX_HRpos = round(mean_OX_HRpos, 4),
                          mean_OX_TNBC = round(mean_OX_TNBC, 4)))

RHO_NOTE <- paste0(
  "DECLARATION 13.7, AND WHAT THIS NUMBER CANNOT SETTLE. The s508 divergence ",
  "- GSE25066's DRFS m1 moving from HR 0.856 [0.687, 1.067] with subtype in ",
  "the model to 0.776 [0.631, 0.955] without it - has TWO causes, not one. ",
  "s508 omits the subtype term AND adds 18 patients who have no callable ",
  "subtype, carrying 2 DRFS events. s490 keeps subtype on 490 patients and ",
  "gives 0.8557 [0.6893, 1.0622], all but identical to the primary's 0.8563, ",
  "so the 25 patients between 490 and 465 are immaterial: the live comparison ",
  "is s490 against s508, confounded between the dropped term and the 18 added ",
  "patients. THE 18 HAVE NO subtype2 VALUE AT ALL and cannot be in this ",
  "correlation. What would settle it is one fit NOT authorised here: m1 on ",
  "the 490 WITHOUT the subtype term. SO THIS CHECK IS NECESSARY AND NOT ",
  "SUFFICIENT: a substantial correlation is CONSISTENT WITH the proposed ",
  "explanation and does not demonstrate it, while a near-zero correlation ",
  "REFUTES it - the one direction in which a single number here is decisive. ",
  "No DRFS model was refitted and no sentence may claim the divergence is ",
  "explained.")
message("\n   ", .sub(RHO_NOTE))

# =============================================================================
# 9. Warnings, reported rather than muffled
# =============================================================================
message("\n9. warnings")
WARNINGS_TAB <- if (length(WARNINGS_SEEN)) dplyr::bind_rows(WARNINGS_SEEN) else
  tibble::tibble(where = character(0), message = character(0))
if (nrow(WARNINGS_TAB)) {
  message("   ", nrow(WARNINGS_TAB), " warning(s) were raised and are ",
          "REPORTED, not muffled. car::vif's 'No intercept' cannot arise ",
          "here - these are logistic fits WITH an intercept.")
  .pr(WARNINGS_TAB %>% dplyr::count(where, message, name = "n"))
} else {
  message("   none raised. car::vif's 'No intercept' warning cannot arise ",
          "here: these are logistic fits WITH an intercept, unlike E40's and ",
          "E41's Cox fits.")
}

# =============================================================================
# 10. G. Save
# =============================================================================
message("\n10. G. saving")

readr::write_csv(COUNTS,      file.path(DIR_TABLES, "E42_counts.csv"))
readr::write_csv(BY_SUBTYPE,  file.path(DIR_TABLES, "E42_by_subtype.csv"))
readr::write_csv(DETECTABLE,  file.path(DIR_TABLES, "E42_detectable.csv"))
readr::write_csv(AUDIT,       file.path(DIR_TABLES, "E42_spine_audit.csv"))
readr::write_csv(LADDER,      file.path(DIR_TABLES, "E42_m1_by_subtype.csv"))
readr::write_csv(META,        file.path(DIR_TABLES, "E42_meta_within_subtype.csv"))
readr::write_csv(AGAINST_E40, file.path(DIR_TABLES, "E42_against_e40.csv"))
readr::write_csv(RHO,         file.path(DIR_TABLES, "E42_rho_ox_subtype.csv"))

READING_RULE <- tibble::tribble(
  ~outcome,                                           ~`the draft sentence becomes`,
  "consistent across subtypes, CIs overlapping",      "unchanged; the scope statement is simply reported",
  "carried by one subtype, others null but underpowered", "unchanged, with the stratum estimates reported beside it and the power limitation stated",
  "opposite signs with NON-OVERLAPPING CIs",           "the claim is QUALIFIED to the subtype where it holds")
message("\n   Declaration 13.6, printed UNFILLED. This script does NOT ",
        "classify the outcome:")
.pr(READING_RULE)
message("   ", .sub(paste0(
  "The THIRD row is the only one that changes the claim, and it requires ",
  "NON-OVERLAPPING intervals rather than one stratum being significant and ",
  "another not. DIFFERING p VALUES ACROSS STRATA OF DIFFERING SIZE ARE NOT ",
  "EVIDENCE OF A DIFFERING EFFECT (13.6).")))

saveRDS(list(
  counts        = COUNTS,
  by_subtype    = BY_SUBTYPE,
  detectable    = DETECTABLE,
  not_fitted    = NOT_FIT,
  flagged       = FLAGGED,
  audit         = AUDIT,
  partition_tab = part_tab,
  m1_by_subtype = LADDER,
  meta          = META,
  already_known = ALREADY,
  against_e40   = AGAINST_E40,
  rho_ox_subtype = RHO,
  rho_note      = RHO_NOTE,
  reading_rule  = READING_RULE,
  warnings_seen = WARNINGS_TAB,
  spec = list(
    declaration = paste0("docs/2026-10-02_E39_respiratory_axis_decomposition",
                         "_declaration.md section 13, committed at dba8ee4"),
    posture = paste0("EXPLORATORY and POST-HOC. A SCOPE CHECK on a claim ",
                     "already being made, not a new analysis of it and not ",
                     "subgroup searching. No stratum was selected; every ",
                     "fittable stratum is reported whatever it shows."),
    cannot_1 = paste0("IT CANNOT PROMOTE OR DEMOTE THE PRIMARY FINDING. ",
                      "E40's three-cohort m1 (-0.162 [-0.281, -0.042], ",
                      "I2 = 0.0%, signs agreeing in all three cohorts) ",
                      "stands as the primary result whatever the strata ",
                      "show. E42 bounds its scope (13.2)."),
    cannot_2 = paste0("NO SUBTYPE RESULT BECOMES A FINDING IN ITS OWN RIGHT. ",
                      "A strong stratum is reported as a stratum of a pooled ",
                      "estimate, never as a separate discovery (13.2)."),
    cannot_3 = paste0("NO INTERACTION TERM AND NO FORMAL EFFECT-MODIFICATION ",
                      "TEST. Declared out, as in section 7, at these event ",
                      "counts. None was fitted (13.2)."),
    cannot_4 = paste0("IT DOES NOT TOUCH METABRIC. The PAM50 split there is ",
                      "X01, a separate activity with a separate posture, and ",
                      "is not part of this declared line (13.2)."),
    reading_rule = paste0(
      "13.6, verbatim. Consistent across subtypes with overlapping CIs -> the ",
      "draft sentence is UNCHANGED and the scope statement is simply ",
      "reported. Carried by one subtype with others null but underpowered -> ",
      "UNCHANGED, with the stratum estimates reported beside it and the power ",
      "limitation stated. Opposite signs with NON-OVERLAPPING CIs -> the ",
      "claim is QUALIFIED to the subtype where it holds. THE THIRD ROW IS THE ",
      "ONLY ONE THAT CHANGES THE CLAIM, and it requires non-overlapping ",
      "intervals rather than one stratum being significant and another not. ",
      "DIFFERING p VALUES ACROSS STRATA OF DIFFERING SIZE ARE NOT EVIDENCE ",
      "OF A DIFFERING EFFECT. This script does not classify the outcome."),
    rung = paste0("m1 ONLY: pcr ~ OX + PROLIF + treatment, fitted within each ",
                  "subtype stratum. subtype is omitted everywhere, constant ",
                  "by construction; treatment is omitted where constant, ",
                  "which is both strata of GSE25066. m0, m2 and m3 are NOT ",
                  "fitted - the ladder's decomposition work is complete and ",
                  "this asks a different question of one rung (13.3)."),
    meta = paste0("WITHIN subtype ACROSS cohorts, not within cohort across ",
                  "subtypes, so the TNBC estimate draws on all three ",
                  "cohorts. Fixed and random always, random the summary at ",
                  "I2 >= 50% (section 9). TWO OF THE FOUR SUBTYPES APPEAR IN ",
                  "GSE194040 ALONE and cannot be pooled at all; their single ",
                  "estimate is reported as single (13.3, 13.5)."),
    power = paste0("13.5. Events per parameter is computed on the MINORITY ",
                   "class, because a logistic fit is limited by whichever ",
                   "outcome is rarer - HRneg_HER2pos has 56 events and only ",
                   "33 non-events, giving 4.7 per parameter, and it is ",
                   "single-cohort besides. A stratum is fitted where the ",
                   "minority class reaches ", MIN_CLASS_FLOOR,
                   " and events per parameter reaches ", EPV_FLOOR,
                   "; below ", EPV_FLAG, " it is FLAGGED but still fitted. ",
                   "The detectable effect is an APPROXIMATION that ignores ",
                   "covariate adjustment, so the true value is LARGER. A ",
                   "NULL IN A LOW-EVENT STRATUM IS REPORTED AS ",
                   "UNINFORMATIVE, NOT AS ABSENCE (section 8)."),
    already_known = paste0(
      "13.4. GSE164458's m1 of -0.205 [-0.423, 0.012] IS already a TNBC-only ",
      "estimate - its subtype term dropped as constant in E40 because the ",
      "cohort is TNBC-only. And in GSE25066 subtype2 and er_primary ",
      "partition the 465 identically (278 / 187, asserted in E40 on the 465, ",
      "NOT on the 470 where subtype still has four levels), so that cohort's ",
      "subtype split and its ER split are the same split."),
    ancillary = "see $rho_note - the check is necessary and NOT sufficient",
    instrument = paste0("SINGLE-INSTRUMENT, as the whole line is. mitoPPS ",
                        "does not exist in these cohorts. No coefficient is ",
                        "compared across cohorts as a raw value (4.2)."),
    no_rescoring = paste0("E39's frames and scores are reused UNCHANGED, so ",
                          "the exposure is byte-identical to E40's. Nothing ",
                          "was rescored and no DRFS model was refitted."),
    frozen = paste0("Nothing written to myc_human_validation (d3ac60e) or ",
                    "myc_mouse. h4_outcome_models.rds was read READ-ONLY.")),
  built = Sys.time()), PATH_E42)

message("   ", PATH_E42)
message("   8 tables in ", DIR_TABLES)
message("\nE42 done. Numbers only; nothing here is interpreted or classified.")
message("E40's primary finding is unchanged by construction. 13.2.\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E42)

  # A. THE COUNTS. Read these before any coefficient.
  x$counts %>% dplyr::mutate(rate = round(rate, 3), epv = round(epv, 1)) %>%
    as.data.frame()
  x$by_subtype %>% as.data.frame()       # which subtypes can be pooled at all
  x$detectable %>% as.data.frame()       # approximate, and an underestimate
  x$not_fitted %>% as.data.frame()
  x$flagged %>% as.data.frame()          # fitted but under 10 EPV

  # B. the spine audit, and the partition 13.4 asserts.
  x$audit %>% as.data.frame()
  x$partition_tab

  # C. m1 within subtype. Read excl_0 and epv together, never p alone.
  x$m1_by_subtype %>% dplyr::transmute(cohort, subtype, n, events,
      logOR = round(estimate, 4), ci_lo = round(ci_lo, 3),
      ci_hi = round(ci_hi, 3), excl_0, p = signif(p, 3),
      epv = round(epv, 1), epv_flag, vif = round(vif_ox, 3)) %>%
    as.data.frame()

  # D. pooled WITHIN subtype ACROSS cohorts. Check `pooled` before reading I2.
  x$meta %>% dplyr::select(subtype, k, pooled, fe_estimate, re_estimate,
                           Q_p, I2, summary_is, signs_agree, signs) %>%
    as.data.frame()

  # E. every stratum beside the finding it bounds.
  x$against_e40 %>% as.data.frame()
  x$already_known %>% as.data.frame()

  # F. the ancillary correlation, and what it cannot settle.
  x$rho_ox_subtype %>% as.data.frame()
  cat(x$rho_note, "\n")

  # The reading rule, UNFILLED. The script does not classify.
  x$reading_rule %>% as.data.frame()
  x$warnings_seen %>% as.data.frame()

  utils::str(x$spec)
}
