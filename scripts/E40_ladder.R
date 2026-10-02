# E40_ladder.R
# =============================================================================
# THE DECLARED LADDER, BOTH ENDPOINTS. IT MAKES NO DECISIONS.
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md,
# amended with the gate result at 5a70707. Read it first. This script fits what
# section 5 specifies, reports every rung whatever it shows, and interprets
# nothing.
#
# Prepared frame and gate: results/e39_prepare_and_gate.rds (E39, cd54d9e).
#
# =============================================================================
# THE GATE PASSED, AND WHAT IT DID AND DID NOT LICENSE
# =============================================================================
# rho(OX, PROLIF_DISJOINT) = 0.358 / 0.381 / 0.493, all < 0.60, so declaration
# section 3 licenses the ladder as specified in all three cohorts.
#
# BUT rho(MYC, PROLIF_DISJOINT) = 0.780 / 0.785 / 0.824. Per declaration 5.3:
#   DECLARED OUT  "whether what remains after proliferation is MYC". That
#                 attribution is NOT available at this collinearity and must
#                 not be written from these numbers.
#   DECLARED IN   whether the OX coefficient is robust to adding a covariate
#                 strongly correlated with one already in the model.
# VIFs for OX are reported at every rung. Neither MYC nor PROLIF may be dropped
# to relieve the collinearity; it is recorded, not engineered away.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - NO model selection. Every rung is reported whatever it shows.
#   - NO ER x OX interaction, no three-way, no post-hoc term. Declared out.
#   - NO p value is described as a near miss, anywhere.
#   - NO per-gene FDR.
#   - It does NOT drop MYC or PROLIF to relieve the collinearity.
#   - It does NOT rescore or re-standardise OX (declaration 4.3).
#   - It does NOT touch METABRIC.
#   - It does NOT interpret. It prints numbers, never verdicts.
#   - It writes NOTHING to myc_human_validation, FROZEN at d3ac60e.
#
# SINGLE-INSTRUMENT THROUGHOUT. mitoPPS does not exist in these cohorts
# (amendment A1). This has not cleared the arm's two-instrument bar and must
# never be written as though it had.
#
# COVERAGE. GSE25066 carries 57 of 89 OXPHOS subunits (64.0%) on GPL96. NO
# coefficient from it may be compared numerically with a TCGA or SCAN-B one
# (declaration 4.1).
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

message("\nE40: the declared ladder, both endpoints\n", strrep("=", 78))
message("IT MAKES NO DECISIONS. Every rung is reported whatever it shows.")

PATH_E39 <- file.path(DIR_RESULTS, "e39_prepare_and_gate.rds")
PATH_H4  <- paste0("/Users/gs/code/myc_human_validation/results/",
                   "h4_outcome_models.rds")
PATH_E40 <- file.path(DIR_RESULTS, "e40_ladder.rds")

# =============================================================================
# 0. CONSTANTS - from the declaration, none invented here
# =============================================================================
CI_LEVEL <- 0.95
CRIT     <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

COHORTS <- c("GSE25066", "GSE194040", "GSE164458")
PRIMARY <- "GSE25066"

# Declaration section 5. Each rung adds one term to the same spine.
RUNG_ADDS <- list(m0 = character(0),
                  m1 = "PROLIF",
                  m2 = c("PROLIF", "MYC"),
                  m3 = c("PROLIF", "MYC", "BUFFER_c"))
RUNGS <- names(RUNG_ADDS)

SPINE_CANDIDATES <- c("subtype_term", "treatment")
VIF_CAVEAT       <- 5         # declaration 5.3

N_PRIMARY <- 465L
E_PRIMARY <- 103L
N_490     <- 490L
N_508     <- 508L

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")
.stop_if(file.exists(PATH_E39), "absent: ", PATH_E39,
         ". Source scripts/E39_prepare_and_gate.R first.")
e39 <- readRDS(PATH_E39)
h4  <- readRDS(PATH_H4)

# The GSE25066 frame already carries DRFS, PROLIF is joined on below.
PS <- e39$prolif
D25 <- e39$frame %>%
  dplyr::left_join(PS[PS$cohort == PRIMARY, c("sample_id", "PROLIF")],
                   by = "sample_id")
.stop_if(sum(is.na(D25$PROLIF)) == 0L, "PROLIF did not join to every GSE25066 row")

# The other two cohorts: H4 frame plus PROLIF.
FR <- stats::setNames(lapply(COHORTS, function(cn) {
  if (cn == PRIMARY) return(D25)
  h4$frames[[cn]] %>%
    dplyr::left_join(PS[PS$cohort == cn, c("sample_id", "PROLIF")],
                     by = "sample_id")
}), COHORTS)
for (cn in COHORTS) .stop_if(sum(is.na(FR[[cn]]$PROLIF)) == 0L,
                             cn, ": PROLIF has NA after the join")
message("   frames: ",
        paste(sprintf("%s %d", COHORTS, vapply(FR, nrow, integer(1))),
              collapse = " | "))

# =============================================================================
# 2. A. SPINE AUDIT - before any fit. Constant terms are omitted EXPLICITLY.
# =============================================================================
# Declaration 5.2: never left for R to drop silently.
message("\n2. A. spine audit")

# GSE25066 uses the collapsed subtype2 and its model frame is the 465.
for (cn in COHORTS) {
  FR[[cn]]$subtype_term <- if (cn == PRIMARY) FR[[cn]]$subtype2 else
    factor(FR[[cn]]$subtype)
}

.n_lev <- function(v) dplyr::n_distinct(v[!is.na(v)])

AUDIT <- dplyr::bind_rows(lapply(COHORTS, function(cn) {
  d <- FR[[cn]]
  if (cn == PRIMARY) d <- d[d$in_model, ]
  tibble::tibble(cohort = cn, n = nrow(d),
                 term = c("OX", SPINE_CANDIDATES),
                 n_levels = c(.n_lev(d$OX), .n_lev(d$subtype_term),
                              .n_lev(d$treatment))) %>%
    dplyr::mutate(constant = n_levels <= 1L,
                  action = ifelse(constant, "OMITTED (constant)", "fitted"))
}))
AUDIT %>% as.data.frame() %>% print(row.names = FALSE)

# The two expected cases are ASSERTED; anything else is reported, not assumed.
.is_const <- function(cn, tm) AUDIT$constant[AUDIT$cohort == cn & AUDIT$term == tm]
.stop_if(.is_const(PRIMARY, "treatment"),
         "treatment is NOT constant in ", PRIMARY, "; declaration 5.2 says it is")
.stop_if(.is_const("GSE164458", "subtype_term"),
         "subtype is NOT constant in GSE164458; declaration 5.2 says it is")
.stop_if(!.is_const(PRIMARY, "subtype_term"),
         PRIMARY, ": subtype2 is constant, which 5.2 does not expect")
.stop_if(!any(AUDIT$constant[AUDIT$term == "OX"]), "OX is constant somewhere")
unexpected <- AUDIT %>%
  dplyr::filter(constant, !(cohort == PRIMARY & term == "treatment"),
                !(cohort == "GSE164458" & term == "subtype_term"))
if (nrow(unexpected)) {
  message("   UNEXPECTED constant term(s), reported not assumed:")
  unexpected %>% as.data.frame() %>% print(row.names = FALSE)
} else {
  message("   the two expected constants are the only ones: treatment in ",
          PRIMARY, ", subtype in GSE164458")
}

# --- subtype2 vs er_primary, declaration 5.2 --------------------------------
dm <- D25[D25$in_model, ]
part_tab <- table(subtype2 = dm$subtype2, er_primary = dm$er_primary)
identical_partition <-
  all((dm$subtype2 == "HRpos_HER2neg") == (dm$er_primary == "P"))
.stop_if(identical_partition,
         "subtype2 and er_primary do NOT partition ", PRIMARY,
         " identically after the HER2 drop. Section 7's stratification assumes ",
         "they do, and section 5.2 says subtype2 is then redundant. STOP.")
.stop_if(nrow(dm) == N_PRIMARY, "the primary frame is ", nrow(dm), ", not ",
         N_PRIMARY)
.stop_if(sum(dm$drfs_event) == E_PRIMARY, "events: ", sum(dm$drfs_event),
         ", not ", E_PRIMARY)
message("   subtype2 and er_primary partition identically (",
        paste(as.vector(part_tab[part_tab > 0]), collapse = " / "),
        "); subtype2 is OMITTED inside each ER stratum")
message("   primary frame: ", nrow(dm), " patients, ", sum(dm$drfs_event),
        " events")

.spine_for <- function(cn) {
  a <- AUDIT[AUDIT$cohort == cn & AUDIT$term %in% SPINE_CANDIDATES, ]
  a$term[!a$constant]
}

# =============================================================================
# 3. The estimator helpers
# =============================================================================
# VIF. OX is always a single continuous term, so its Df is 1 and its GVIF IS
# the ordinary VIF - the two coincide at 1 df. The measure is named in the
# output either way. car::vif needs at least two terms; with one it is
# undefined and is reported as such rather than faked.
.vif_ox <- function(fit) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2L) return(c(vif_ox = NA_real_, vif_df = 1))
  v <- tryCatch(car::vif(fit), error = function(e) NULL)
  if (is.null(v)) return(c(vif_ox = NA_real_, vif_df = NA_real_))
  if (is.matrix(v)) c(vif_ox = unname(v["OX", 1]), vif_df = unname(v["OX", 2]))
  else c(vif_ox = unname(v[["OX"]]), vif_df = 1)
}
.vif_measure <- function(fit) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2L) return("undefined (single predictor)")
  v <- tryCatch(car::vif(fit), error = function(e) NULL)
  if (is.null(v)) "car::vif failed"
  else if (is.matrix(v)) "GVIF (generalised; OX has Df 1 so GVIF == VIF)"
  else "VIF"
}

.tidy_ox <- function(fit, cohort, rung, endpoint, set, n, events = NA_integer_) {
  co <- summary(fit)$coefficients
  .stop_if("OX" %in% rownames(co), cohort, " ", rung, ": no OX coefficient")
  est <- co["OX", 1]
  se  <- if (inherits(fit, "coxph")) co["OX", "se(coef)"] else co["OX", 2]
  p   <- co["OX", ncol(co)]
  vf  <- .vif_ox(fit)
  tibble::tibble(
    endpoint = endpoint, set = set, cohort = cohort, rung = rung,
    n = n, events = events,
    estimate = est, se = se,
    ci_lo = est - CRIT * se, ci_hi = est + CRIT * se, p = p,
    vif_ox = unname(vf["vif_ox"]), vif_df = unname(vf["vif_df"]),
    vif_measure = .vif_measure(fit),
    vif_above_caveat = isTRUE(unname(vf["vif_ox"]) > VIF_CAVEAT),
    terms = paste(attr(stats::terms(fit), "term.labels"), collapse = " + "))
}

# script 13's .meta(), reproduced verbatim so E40's pooling is the same
# DerSimonian-Laird estimator the arm has used throughout. metafor is not
# installed and is not needed.
.meta <- function(est, se, label) {
  v <- se^2; w <- 1 / v; k <- length(est)
  if (k < 2L) {
    return(tibble::tibble(rung = label, k = k,
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
  tibble::tibble(rung = label, k = k,
    fe_estimate = fe, fe_ci_lo = fe - CRIT*sqrt(vfe),
    fe_ci_hi = fe + CRIT*sqrt(vfe), fe_p = 2*stats::pnorm(-abs(fe/sqrt(vfe))),
    re_estimate = re, re_ci_lo = re - CRIT*sqrt(vre),
    re_ci_hi = re + CRIT*sqrt(vre), re_p = 2*stats::pnorm(-abs(re/sqrt(vre))),
    Q = Q, Q_p = stats::pchisq(Q, df, lower.tail = FALSE),
    tau2 = tau2, I2 = 100*max(0, (Q - df)/Q))
}

# =============================================================================
# 4. B. THE pCR LADDER - logistic, all three cohorts
# =============================================================================
message("\n4. B. the pCR ladder (logistic)")

PCR <- dplyr::bind_rows(lapply(COHORTS, function(cn) {
  d <- FR[[cn]]
  if (cn == PRIMARY) d <- d[d$in_model, ]
  sp <- .spine_for(cn)
  dplyr::bind_rows(lapply(RUNGS, function(rg) {
    rhs <- c("OX", sp, RUNG_ADDS[[rg]])
    f   <- stats::reformulate(rhs, response = "pcr")
    fit <- stats::glm(f, data = d, family = stats::binomial())
    .tidy_ox(fit, cn, rg, "pCR", "primary", n = stats::nobs(fit),
             events = sum(fit$y == 1))
  }))
}))
PCR %>% dplyr::transmute(cohort, rung, n, events,
    OX = round(estimate, 4), se = round(se, 4),
    ci = paste0("[", round(ci_lo, 3), ", ", round(ci_hi, 3), "]"),
    p = signif(p, 3), vif_ox = round(vif_ox, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   VIF measure: ", unique(PCR$vif_measure[!is.na(PCR$vif_ox)])[1])
message("   spine fitted per cohort: ",
        paste(sprintf("%s = OX + %s", COHORTS,
                      vapply(COHORTS, function(cn)
                        paste(.spine_for(cn), collapse = " + "), character(1))),
              collapse = " | "))

# =============================================================================
# 5. C. META-ANALYSIS of the pCR ladder, per rung
# =============================================================================
# Declaration section 9: BOTH fixed and random are reported, always. At
# I2 >= 50% the random-effects estimate is the summary and the fixed-effect one
# sits beside it, never instead of it. Sign agreement across the three cohorts
# is a condition for reading a signal at all.
message("\n5. C. meta-analysis, per rung, three cohorts")

META <- dplyr::bind_rows(lapply(RUNGS, function(rg) {
  d <- PCR[PCR$rung == rg, ]
  m <- .meta(d$estimate, d$se, rg)
  m$signs_agree <- length(unique(sign(d$estimate))) == 1L
  m$signs <- paste(sprintf("%s %s", d$cohort,
                           ifelse(d$estimate > 0, "+", "-")), collapse = " ")
  m$summary_is <- ifelse(!is.na(m$I2) & m$I2 >= 50, "random effects",
                         "fixed effect")
  m
}))
META %>% dplyr::transmute(rung, k,
    fe = round(fe_estimate, 4),
    fe_ci = paste0("[", round(fe_ci_lo, 3), ", ", round(fe_ci_hi, 3), "]"),
    re = round(re_estimate, 4),
    re_ci = paste0("[", round(re_ci_lo, 3), ", ", round(re_ci_hi, 3), "]"),
    Q = round(Q, 3), Q_p = signif(Q_p, 3), tau2 = round(tau2, 4),
    I2 = round(I2, 1), summary_is, signs_agree, signs) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. D. THE DRFS LADDER - Cox, GSE25066 only
# =============================================================================
message("\n6. D. the DRFS ladder (Cox), ", PRIMARY, " only")

.drfs_sets <- list(
  primary = D25[D25$in_model, ],                       # 465
  s490    = D25[!is.na(D25$subtype), ],                # 508 less the 18
  s508    = D25)                                       # all 508, no subtype
.stop_if(nrow(.drfs_sets$primary) == N_PRIMARY, "primary set is not ", N_PRIMARY)
.stop_if(nrow(.drfs_sets$s490) == N_490, "the 490 set is ",
         nrow(.drfs_sets$s490), ", not ", N_490)
.stop_if(nrow(.drfs_sets$s508) == N_508, "the 508 set is not ", N_508)

# The 490 keeps the UNCOLLAPSED 4-level subtype by construction - it is the 508
# less the 18 without a callable subtype, so the HER2-positive patients are in
# it. Its composition is reported so the author sees what was fitted.
message("   the 490's subtype composition, as fitted:")
.drfs_sets$s490 %>% dplyr::count(subtype, name = "n") %>%
  as.data.frame() %>% print(row.names = FALSE)

# The fits are kept, not just their coefficients: cox.zph in section E needs
# the model objects and refitting them there would be a second, separate fit.
FITS <- list()

.fit_drfs <- function(d, set, spine) {
  rows <- vector("list", length(RUNGS))
  FITS[[set]] <<- stats::setNames(vector("list", length(RUNGS)), RUNGS)
  for (i in seq_along(RUNGS)) {
    rg  <- RUNGS[[i]]
    rhs <- c("OX", spine, RUNG_ADDS[[rg]])
    f   <- stats::reformulate(rhs,
             response = "survival::Surv(drfs_time, drfs_event)")
    fit <- survival::coxph(f, data = d)
    FITS[[set]][[rg]] <<- fit
    rows[[i]] <- .tidy_ox(fit, PRIMARY, rg, "DRFS", set, n = fit$n,
                          events = fit$nevent)
  }
  dplyr::bind_rows(rows)
}

DRFS <- dplyr::bind_rows(
  .fit_drfs(.drfs_sets$primary, "primary", "subtype2"),
  .fit_drfs(.drfs_sets$s490,    "s490",    "subtype"),
  .fit_drfs(.drfs_sets$s508,    "s508",    character(0)))
DRFS %>% dplyr::transmute(set, rung, n, events,
    OX = round(estimate, 4), se = round(se, 4),
    HR = round(exp(estimate), 4),
    ci = paste0("[", round(exp(ci_lo), 3), ", ", round(exp(ci_hi), 3), "]"),
    p = signif(p, 3), vif_ox = round(vif_ox, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   treatment OMITTED throughout (constant in ", PRIMARY, ")")
message("   CI shown on the HR scale; estimate and se are on the log-hazard scale")

# =============================================================================
# 7. E. PROPORTIONAL HAZARDS - reported whether or not it holds
# =============================================================================
# Declaration 6.3. cox.zph on every DRFS model, global test and the OX term.
message("\n7. E. proportional hazards (cox.zph), every DRFS model")

PH <- dplyr::bind_rows(lapply(names(FITS), function(set) {
  dplyr::bind_rows(lapply(RUNGS, function(rg) {
    z <- survival::cox.zph(FITS[[set]][[rg]])
    tb <- z$table
    tibble::tibble(set = set, rung = rg,
                   global_chisq = unname(tb["GLOBAL", "chisq"]),
                   global_df = unname(tb["GLOBAL", "df"]),
                   global_p = unname(tb["GLOBAL", "p"]),
                   ox_chisq = unname(tb["OX", "chisq"]),
                   ox_df = unname(tb["OX", "df"]),
                   ox_p = unname(tb["OX", "p"])) %>%
      dplyr::mutate(ox_ph_holds = ox_p >= 0.05,
                    global_ph_holds = global_p >= 0.05)
  }))
}))
PH %>% dplyr::transmute(set, rung,
    global_p = signif(global_p, 3), global_ph_holds,
    ox_p = signif(ox_p, 3), ox_ph_holds) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- the RMST branch, and a dependency this repo does not have --------------
# Declaration 6.3: if PH fails for OX, report RMST at 3 y ALONGSIDE the hazard
# ratio. The script does NOT invent a way to do it:
#   - survRM2 is NOT installed (checked, not installed by this script), and
#   - OX is a standardised CONTINUOUS exposure. An RMST contrast needs two
#     arms, so reporting one would require dichotomising OX, and NO cut is
#     declared anywhere in the declaration.
# Making that cut here would be a modelling decision, and this script makes
# none. The branch therefore REPORTS what is required and what is missing.
PH_OX_FAILS <- PH %>% dplyr::filter(!ox_ph_holds)
RMST_NOTE <- if (nrow(PH_OX_FAILS) == 0L) {
  "PH holds for OX at every rung in every DRFS set; declaration 6.3's RMST branch is not triggered."
} else {
  paste0(
    "PH FAILS FOR OX in ", nrow(PH_OX_FAILS), " of ", nrow(PH),
    " DRFS models (", paste(sprintf("%s/%s", PH_OX_FAILS$set, PH_OX_FAILS$rung),
                            collapse = ", "),
    "). Declaration 6.3 requires RMST at 3 y ALONGSIDE the hazard ratio there. ",
    "THIS SCRIPT DID NOT COMPUTE IT, for two reasons, neither of which it may ",
    "resolve on its own: survRM2 is not installed in this environment, and OX ",
    "is a standardised CONTINUOUS exposure, so an RMST contrast would require ",
    "dichotomising it at a cut that the declaration does not specify. Choosing ",
    "that cut is a modelling decision and this script makes none.")
}
message("\n   ", paste(strwrap(RMST_NOTE, width = 74), collapse = "\n   "))

# =============================================================================
# 8. F. ER STRATIFICATION - descriptive, DRFS, primary frame
# =============================================================================
# Declaration section 7. subtype2 is constant inside each stratum (5.2) and is
# OMITTED explicitly. NO ER x OX interaction is fitted: it is declared out.
message("\n8. F. ER stratification (descriptive), DRFS")

# subtype2 is OMITTED in BOTH stratifications, for two DIFFERENT reasons
# (declaration 7.1, added after this script's first run stopped here):
#   er_primary    it is a relabelling of er_status_ihc, so it is CONSTANT
#                 inside the stratum and the omission is forced. Asserted.
#   er_sensitive  the two ER calls disagree on 54 patients, so subtype2
#                 genuinely VARIES inside the esr1 strata. It is omitted BY
#                 DECISION, so the sensitivity is the same model as the primary
#                 stratification and the two differ only in how patients are
#                 assigned. Asserting constancy here would assert something
#                 false, so the composition is REPORTED instead.
.fit_stratum <- function(d, er_var, level, set_label, assert_constant) {
  dd <- d[d[[er_var]] == level, ]
  n_lev <- dplyr::n_distinct(dd$subtype2[!is.na(dd$subtype2)])
  if (assert_constant) {
    .stop_if(n_lev <= 1L,
             "subtype2 is NOT constant inside ", er_var, " == ", level,
             "; declaration 5.2 says it is for er_primary")
  }
  dplyr::bind_rows(lapply(RUNGS, function(rg) {
    rhs <- c("OX", RUNG_ADDS[[rg]])   # subtype2 omitted - see the note above
    f <- stats::reformulate(rhs,
           response = "survival::Surv(drfs_time, drfs_event)")
    fit <- survival::coxph(f, data = dd)
    .tidy_ox(fit, PRIMARY, rg, "DRFS", set_label, n = fit$n,
             events = fit$nevent) %>%
      dplyr::mutate(er_variable = er_var, er_level = level,
                    subtype2_levels_in_stratum = n_lev,
                    subtype2_omitted_because =
                      if (assert_constant) "constant (forced)" else
                        "decision (declaration 7.1); it VARIES here")
  }))
}

dm <- D25[D25$in_model, ]
ER_STRAT <- dplyr::bind_rows(
  .fit_stratum(dm, "er_primary",   "P", "ER-positive (er_status_ihc)", TRUE),
  .fit_stratum(dm, "er_primary",   "N", "ER-negative (er_status_ihc)", TRUE),
  .fit_stratum(dm, "er_sensitive", "P", "ER-positive (esr1_status)",   FALSE),
  .fit_stratum(dm, "er_sensitive", "N", "ER-negative (esr1_status)",   FALSE))
ER_STRAT %>% dplyr::transmute(set, rung, n, events,
    OX = round(estimate, 4), HR = round(exp(estimate), 4),
    ci = paste0("[", round(exp(ci_lo), 3), ", ", round(exp(ci_hi), 3), "]"),
    p = signif(p, 3), vif_ox = round(vif_ox, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   subtype2 OMITTED inside every stratum, for two different reasons:")
message("     er_primary   : CONSTANT there (5.2). Asserted.")
message("     er_sensitive : it VARIES there. Omitted BY DECISION (7.1) so the")
message("                    sensitivity is the same model as the primary.")
message("   NO ER x OX interaction was fitted. It is declared out (section 7).")

# What the esr1 omission costs, reported where it is read (declaration 7.1).
ESR1_COST <- dm %>%
  dplyr::group_by(er_sensitive) %>%
  dplyr::summarise(n = dplyr::n(), events = sum(drfs_event),
                   n_HRpos_HER2neg = sum(subtype2 == "HRpos_HER2neg"),
                   n_TNBC = sum(subtype2 == "TNBC"), .groups = "drop")
DISCORDANT <- dm %>% dplyr::filter(er_primary != er_sensitive)
message("\n   WHAT THE esr1 OMISSION COSTS - subtype2 composition inside each",
        " esr1 stratum:")
ESR1_COST %>% as.data.frame() %>% print(row.names = FALSE)
message("   the two ER calls disagree on ", nrow(DISCORDANT), " of ", nrow(dm),
        " patients, carrying ", sum(DISCORDANT$drfs_event), " events")
message("   (declaration 7's 'roughly ten' is the NET MARGIN, corrected at ",
        "29aa4d6)")
DISCORDANT %>% dplyr::count(er_primary, er_sensitive, subtype2, name = "n") %>%
  as.data.frame() %>% print(row.names = FALSE)

# The declaration's own words, attached to the ER-positive output.
ER_POS_WORDS <- paste0(
  "At 40 events this stratum detects nothing below HR 1.56 to 1.87. A NULL IN ",
  "THE ER-POSITIVE STRATUM OF GSE25066 CARRIES NO INFORMATION and must be ",
  "reported in those words (declaration section 8).")
message("\n   ON THE ER-POSITIVE STRATUM:")
message("   ", paste(strwrap(ER_POS_WORDS, width = 74), collapse = "\n   "))

# =============================================================================
# 9. G. THE SIGN PAIR - one table, same 465 patients, no classification
# =============================================================================
message("\n9. G. the sign pair, ", PRIMARY, ", the same ", N_PRIMARY,
        " patients")

SIGN_PAIR <- PCR %>%
  dplyr::filter(cohort == PRIMARY) %>%
  dplyr::transmute(rung, n_pcr = n, pcr_events = events,
                   pcr_est = estimate, pcr_lo = ci_lo, pcr_hi = ci_hi,
                   pcr_p = p) %>%
  dplyr::left_join(
    DRFS %>% dplyr::filter(set == "primary") %>%
      dplyr::transmute(rung, n_drfs = n, drfs_events = events,
                       drfs_est = estimate, drfs_lo = ci_lo, drfs_hi = ci_hi,
                       drfs_p = p),
    by = "rung")
.stop_if(all(SIGN_PAIR$n_pcr == SIGN_PAIR$n_drfs),
         "the two endpoints are not on the same patients; the sign pair is ",
         "not interpretable (declaration 6.1)")
SIGN_PAIR %>% dplyr::transmute(rung, n = n_pcr,
    pCR_logOR = round(pcr_est, 4),
    pCR_ci = paste0("[", round(pcr_lo, 3), ", ", round(pcr_hi, 3), "]"),
    DRFS_logHR = round(drfs_est, 4),
    DRFS_ci = paste0("[", round(drfs_lo, 3), ", ", round(drfs_hi, 3), "]")) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Declaration 9.1, printed UNFILLED. The script does not classify.
READING_TABLE <- tibble::tribble(
  ~pCR,        ~DRFS,   ~`consistent with`,
  "better",    "worse", "the proliferation route",
  "no effect", "worse", "the metastasis route the mouse data predict",
  "better",    "better","an unbuffered apoptotic liability")
message("\n   Declaration 9.1, printed UNFILLED. NONE of these is a test of ",
        "anything,\n   and this script does not classify the result:")
READING_TABLE %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 10. H. Save
# =============================================================================
message("\n10. saving")

readr::write_csv(PCR,      file.path(DIR_TABLES, "E40_pcr_ladder.csv"))
readr::write_csv(META,     file.path(DIR_TABLES, "E40_pcr_meta.csv"))
readr::write_csv(DRFS,     file.path(DIR_TABLES, "E40_drfs_ladder.csv"))
readr::write_csv(PH,       file.path(DIR_TABLES, "E40_proportional_hazards.csv"))
readr::write_csv(ER_STRAT, file.path(DIR_TABLES, "E40_er_stratified.csv"))
readr::write_csv(SIGN_PAIR, file.path(DIR_TABLES, "E40_sign_pair.csv"))

saveRDS(list(
  audit         = AUDIT,
  partition_tab = part_tab,
  pcr           = PCR,
  meta          = META,
  drfs          = DRFS,
  ph            = PH,
  ph_ox_fails   = PH_OX_FAILS,
  rmst_note     = RMST_NOTE,
  er_stratified = ER_STRAT,
  esr1_cost     = ESR1_COST,
  discordant    = DISCORDANT,
  er_pos_words  = ER_POS_WORDS,
  sign_pair     = SIGN_PAIR,
  reading_table = READING_TABLE,
  rungs         = RUNG_ADDS,
  spec = list(
    declaration = paste0("docs/2026-10-02_E39_respiratory_axis_decomposition",
                         "_declaration.md, amended at 5a70707"),
    e39 = "results/e39_prepare_and_gate.rds (E39, cd54d9e)",
    posture = paste0("EXPLORATORY and POST-HOC. A descriptive decomposition, ",
                     "not a fifth hypothesis. Nothing here carries a pass, a ",
                     "fail or a falsification criterion. No rung was selected; ",
                     "every rung is reported whatever it shows."),
    instrument = paste0("SINGLE-INSTRUMENT throughout. mitoPPS does not exist ",
                        "in these cohorts (amendment A1); this has NOT cleared ",
                        "the arm's two-instrument bar and must never be ",
                        "written as though it had."),
    coverage_caveat = paste0("GSE25066 carries 57 of 89 OXPHOS subunits ",
                             "(64.0%) on GPL96. NO coefficient from it may be ",
                             "compared numerically with a TCGA or SCAN-B one ",
                             "(declaration 4.1)."),
    collinearity = paste0("rho(MYC, PROLIF_DISJOINT) = 0.780 / 0.785 / 0.824. ",
                          "Per declaration 5.3 the m1 -> m2 step does NOT ",
                          "decompose by variable: 'whether what remains after ",
                          "proliferation is MYC' is DECLARED OUT and must not ",
                          "be written. What is declared IN is whether the OX ",
                          "coefficient is robust to adding a covariate ",
                          "strongly correlated with one already in the model. ",
                          "VIF for OX is reported at every rung; above ", 
                          VIF_CAVEAT, " it is a caveat on that rung, not a ",
                          "reason to drop a term."),
    gate = paste0("rho(OX, PROLIF_DISJOINT) = 0.358 / 0.381 / 0.493, all ",
                  "< 0.60, so declaration section 3 licenses the ladder as ",
                  "specified in all three cohorts."),
    er_interaction = "NOT fitted. Declared out (section 7).",
    subtype2_in_strata = paste0(
      "subtype2 is omitted in BOTH stratifications for DIFFERENT reasons ",
      "(declaration 7.1). In the er_primary strata it is constant and the ",
      "omission is forced; asserted. In the esr1 strata it VARIES - the two ",
      "ER calls disagree on 54 of 465 patients carrying 11 events - and is ",
      "omitted BY DECISION so the sensitivity is the same model as the ",
      "primary. The esr1 estimates are therefore UNADJUSTED for subtype and ",
      "the reclassified patients carry their subtype imbalance into them."),
    rmst = RMST_NOTE,
    frozen = "Nothing written to myc_human_validation (d3ac60e)."),
  built = Sys.time()), PATH_E40)

message("   ", PATH_E40)
message("   6 tables in ", DIR_TABLES)
message("\nE40 done. Numbers only; nothing here is interpreted.\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E40)

  # The spine audit, and the partition the stratification assumes.
  x$audit %>% as.data.frame(); x$partition_tab

  # B. the pCR ladder. Watch the VIF for OX at m2 and m3.
  x$pcr %>% dplyr::transmute(cohort, rung, n, events, estimate = round(estimate, 4),
      ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3), p = signif(p, 3),
      vif_ox = round(vif_ox, 3), vif_above_caveat) %>% as.data.frame()

  # C. and whether the three cohorts agree in sign at each rung.
  x$meta %>% dplyr::select(rung, k, fe_estimate, re_estimate, Q_p, I2,
                           summary_is, signs_agree, signs) %>% as.data.frame()

  # D. the DRFS ladder, on the HR scale, with its two sensitivities.
  x$drfs %>% dplyr::transmute(set, rung, n, events, HR = round(exp(estimate), 4),
      lo = round(exp(ci_lo), 3), hi = round(exp(ci_hi), 3), p = signif(p, 3),
      vif_ox = round(vif_ox, 3)) %>% as.data.frame()

  # E. proportional hazards, and the RMST branch.
  x$ph %>% as.data.frame()
  x$ph_ox_fails %>% as.data.frame()
  cat(x$rmst_note, "\n")

  # F. ER strata. Read the ER-positive rows with x$er_pos_words in front of them.
  x$er_stratified %>% dplyr::transmute(set, rung, n, events,
      HR = round(exp(estimate), 4), lo = round(exp(ci_lo), 3),
      hi = round(exp(ci_hi), 3), p = signif(p, 3),
      subtype2_omitted_because) %>% as.data.frame()
  cat(x$er_pos_words, "\n")

  # What the esr1 omission costs: subtype2 varies there, and by how much.
  x$esr1_cost %>% as.data.frame()
  x$discordant %>% dplyr::count(er_primary, er_sensitive, subtype2,
                                name = "n") %>% as.data.frame()

  # G. the sign pair, unclassified.
  x$sign_pair %>% as.data.frame()
  x$reading_table %>% as.data.frame()

  utils::str(x$spec)
}
