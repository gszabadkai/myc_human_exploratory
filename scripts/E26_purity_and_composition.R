# scripts/E26_purity_and_composition.R
# =============================================================================
# E26 - does the MYC separation survive composition adjustment?
# =============================================================================
#
# THE QUESTION. E25 found that MYC-high and MYC-low human tumours order the
# mitochondrial panel differently, and that the two orderings match two
# different mouse groups. Every number there is UNADJUSTED. MYC-high tumours are
# more proliferative and plausibly more tumour-cell-pure, and adipose tissue is
# OXPHOS-high and FAO-high. So part of that separation could be COMPOSITION
# rather than programme. This script asks whether it is.
#
# TRAP 15 GOVERNS THIS SCRIPT AND IS THE REASON IT EXISTS
# -------------------------------------------------------
# From the mouse fat pad: Mki67 went tau +0.313, p 0.045 -> -0.096, p 0.715
# after composition adjustment, while the adjusters explained only 15% of its
# variance. An adjustment that destroys an association it barely explains is
# over-correcting, and the only way to know is to carry a quantity whose answer
# is not in doubt THROUGH the adjustment and READ IT FIRST.
#
#   POSITIVE CONTROL: proliferation. MYC drives proliferation. That is not in
#   doubt and it is not a composition artefact. If purity adjustment destroys
#   it, the adjustment is removing biology and NOTHING below it is readable.
#   It is read FIRST, before any OXPHOS number, and it is a GATE.
#
#   COMPOSITION TRACER: an adipocyte marker score. It is nearly pure
#   composition in bulk breast tissue. It shows how much traction each
#   adjustment actually has on composition. If purity and leukocyte fraction
#   leave it substantially intact, then THEY ARE NOT AN ADEQUATE COMPOSITION
#   CONTROL and the adipose-adjusted models become primary. That rule is fixed
#   below, before any number.
#
# WHY ATTENUATION IS READ AGAINST THE WHOLE PANEL, NOT ALONE
# ----------------------------------------------------------
# Any adjustment shrinks every coefficient somewhat. A number shrinking is not
# evidence. What would be evidence is the PRIMARY SET shrinking MORE THAN THE
# PANEL DOES - that is composition acting specifically on it. So the attenuation
# of all 142 pathways is computed and the primary's position in that
# distribution is the read. Same logic as E21's size-matched null.
#
# TCGA ONLY, AND THAT IS STRUCTURAL
# ----------------------------------
# TCGA carries purity, leukocyte fraction and ploidy. SCAN-B carries NONE and
# nothing is imputed (CLAUDE.md trap 2). So this arm cannot replicate, and
# SCAN-B is carried unadjusted as the cohort that cannot be checked. That
# asymmetry is a limitation of the data, not a choice, and the note says so.
#
# GENOME DOUBLINGS ARE NOT A COVARIATE HERE. The author ruled that out and it is
# not added. Ploidy is available and also not used; only the composition
# adjusters the question is about enter the models.
#
# WHAT IS UNCHANGED FROM E25
# --------------------------
# The estimand, the MYC groups, the panel, the mouse artefact. The residual
# machinery is proven against E25 in GATE 0: adjusting for an INTERCEPT ONLY
# must reproduce E25's percentiles exactly, because mitoPPS is centred at 1 per
# pathway so a group-mean ordering and a group-mean-residual ordering are the
# same object.
#
# RANK, NEVER VALUE across the species boundary. No pooling. TCGA and SCAN-B
# never averaged.
#
# EXPLORATORY, POST-HOC, DESCRIPTIVE. Nothing here is pre-registered and nothing
# here is a hypothesis test. Every reading rule below is fixed before any number.
#
# SPECIES: the human side is human. No ortholog function is called here or
# anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE26: does the MYC separation survive composition adjustment?\n",
        strrep("=", 78))

PATH_E26     <- file.path(DIR_RESULTS, "purity_and_composition.rds")
PATH_E26_TAB <- file.path(DIR_TABLES,  "E26_adjusted_contrasts.csv")
DIR_E26_FIG  <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_E25     <- file.path(DIR_RESULTS, "mitopps_group_rankings.rds")
PATH_MOUSE   <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")

set.seed(1)
NBOOT <- 5000L

# =============================================================================
# 0. CONSTANTS AND THE READING RULES, FIXED BEFORE ANY NUMBER
# =============================================================================

MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"
PRIMARY   <- "OXPHOS subunits"
SECONDARY <- "OXPHOS"
MTRIB     <- "Mitochondrial ribosome"
MTDNA     <- "mtDNA-encoded OXPHOS subunits"
READOUTS  <- c(PRIMARY, SECONDARY, MTRIB, MTDNA)

TERTILE   <- 1/3          # the E25 split, reused unchanged
PURITY_HI <- 0.65         # the high-purity subset threshold, fixed in advance

# --- the pre-declared thresholds ---------------------------------------------
CTRL_RETAIN <- 0.50   # the POSITIVE control must keep >= this share of its
                      # unadjusted magnitude, same sign, interval clear of zero
TRACER_LEFT <- 0.50   # if an adjustment leaves >= this share of the composition
                      # tracer's MYC difference, it is NOT an adequate
                      # composition control and the adipose models become primary
PANEL_TAIL  <- 0.05   # the primary is an attenuation OUTLIER if it falls outside
                      # the panel's central 90 per cent

# --- the adipocyte marker panel ----------------------------------------------
# DECLARED HERE AND NOT FROM A PINNED SNAPSHOT. It is a COVARIATE, never a
# scored set and never a claim-bearing object, which is why an in-script
# definition is acceptable where CLAUDE.md otherwise says consume the snapshots.
# Sixteen canonical adipocyte genes, all present in the TCGA matrix.
#
# WHAT IT CANNOT SEPARATE, stated up front: "less adipose tissue in the section"
# from "tumour cells expressing fewer adipocyte genes". In bulk breast tissue
# ADIPOQ, LEP, PLIN1 and FABP4 are overwhelmingly adipocyte-derived, so it is
# mostly composition - but "mostly" is the honest word and the note repeats it.
ADIPO_MARKERS <- c("ADIPOQ", "LEP", "PLIN1", "PLIN4", "FABP4", "CFD",
                   "CIDEC", "LIPE", "PPARG", "CD36", "AQP7", "GPD1",
                   "TRARG1", "THRSP", "LPL", "PNPLA2")

# --- the adjustment models ---------------------------------------------------
# M0 is the E25 reading. M1 is the adjustment the plan named. M2 adds the axis
# the exposure check says actually differs between the groups. M3 is both.
MODELS <- list(
  `M0 unadjusted`              = character(0),
  `M1 purity + leukocyte`      = c("purity", "leukocyte"),
  `M2 adipose`                 = c("adipose"),
  `M3 purity + leuk + adipose` = c("purity", "leukocyte", "adipose"))
MODEL_MAIN <- "M1 purity + leukocyte"

.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}

RULES <- list(
  question = paste(
    "does the MYC-high against MYC-low separation in the mitochondrial",
    "ordering survive adjustment for tumour composition?"),
  control_first = paste(
    "the POSITIVE control (proliferation) is read before any OXPHOS number and",
    "is a GATE: if it does not survive the adjustment, the adjustment is",
    "over-correcting and nothing below it is readable"),
  tracer = paste(
    "the composition tracer (adipocyte markers) measures how much traction an",
    "adjustment has on composition; an adjustment that leaves it intact is not",
    "a composition control"),
  panel_null = paste(
    "attenuation is read against the whole panel's attenuation distribution,",
    "never alone - every adjustment shrinks every coefficient"),
  tcga_only = paste(
    "SCAN-B has no purity estimate and nothing is imputed, so this arm cannot",
    "replicate; SCAN-B is carried unadjusted as the cohort that cannot be checked"),
  not_covariates = "genome doublings ruled out by the author; ploidy not used",
  rank_only = "no mitoPPS value crosses the species boundary",
  exploratory = "post-hoc and descriptive; not a test of any hypothesis",
  n3 = "transcript associations; 'primed' is never written of a transcript")

message("\n0. rules fixed before any number")
message("   G0   intercept-only residuals must reproduce E25's percentiles exactly")
message("   G1   POSITIVE CONTROL read FIRST: proliferation must keep >= ",
        round(100 * CTRL_RETAIN), " pct of its")
message("        unadjusted MYC difference, same sign, interval clear of zero")
message("   G2   TRACER: an adjustment leaving >= ", round(100 * TRACER_LEFT),
        " pct of the adipose MYC difference is")
message("        NOT an adequate composition control -> the adipose models become primary")
message("   R1   the primary's attenuation against the PANEL's distribution")
message("   R2   does the adjusted cross-species agreement survive?")
message("   TCGA ONLY. SCAN-B carried unadjusted; it structurally cannot replicate")

# =============================================================================
# 1. INPUTS
# =============================================================================

message("\n1. inputs")

stopifnot(file.exists(PATH_E25))
e25 <- readRDS(PATH_E25)
message("   E25 object loaded: ", length(e25), " elements, analysis date ",
        as.character(e25$analysis_date))

md5_now <- unname(tools::md5sum(PATH_MOUSE))
stopifnot(identical(md5_now, MOUSE_MD5))
mo <- readRDS(PATH_MOUSE)
message("   mouse artefact md5 matches: ", md5_now)

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
cvd  <- readRDS(here::here("data", "from_validation",
                          "tcga_brca_covariates.rds"))$covariates
tl   <- readRDS(PATH_TCGA_LINEAR)
stopifnot(identical(tl$scale, "linear_deseq2_normalised"))

UT   <- mito$mitopps_universe
ID_T <- colnames(UT)
shared <- e25$intersection$shared
stopifnot(all(shared %in% rownames(UT)), length(shared) == 142L)
UT <- UT[shared, , drop = FALSE]
message("   TCGA mitopps_universe: ", nrow(UT), " pathways x ", ncol(UT), " samples")

# --- covariates --------------------------------------------------------------
k        <- match(ID_T, cvd$patient)
stopifnot(!anyNA(k))
purity   <- cvd$purity[k]
leuk     <- cvd$leukocyte_fraction[k]
amp      <- cvd$MYC_amp[k]
names(purity) <- names(leuk) <- names(amp) <- ID_T

# --- the composition tracer --------------------------------------------------
stopifnot(all(ADIPO_MARKERS %in% rownames(tl$mat)))
LA      <- log2(tl$mat[ADIPO_MARKERS, ID_T, drop = FALSE] + 1)
adipose <- colMeans(t(scale(t(LA))))          # mean z across markers, per sample
stopifnot(identical(names(adipose), ID_T))
rm(tl); invisible(gc(verbose = FALSE))
message("   adipocyte tracer: ", length(ADIPO_MARKERS),
        " markers, all present; in-script definition, NOT a pinned snapshot")

# --- the positive control ----------------------------------------------------
PROLIF <- "PROLIF_DISJOINT"
stopifnot(PROLIF %in% rownames(mito$gsva_cov))
prolif <- as.numeric(mito$gsva_cov[PROLIF, ID_T])
names(prolif) <- ID_T
message("   positive control: ", PROLIF, " (GSVA on VST - never mixed with the ",
        "linear mitoPPS input)")

# --- the MYC groups, E25's split reused unchanged -----------------------------
est <- as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T]); names(est) <- ID_T
qq  <- stats::quantile(est, c(TERTILE, 1 - TERTILE), names = FALSE)
grp <- rep(NA_character_, length(est)); names(grp) <- ID_T
grp[est <= qq[1]] <- "MYC-low"
grp[est >= qq[2]] <- "MYC-high"
stopifnot(identical(
  as.integer(c(sum(grp == "MYC-low", na.rm = TRUE), sum(grp == "MYC-high", na.rm = TRUE))),
  as.integer(e25$myc_groups$n[e25$myc_groups$cohort == "TCGA"])))
message("   MYC groups reproduce E25 exactly: low ", sum(grp == "MYC-low", na.rm = TRUE),
        " | high ", sum(grp == "MYC-high", na.rm = TRUE))

# --- the analysis set: two MYC groups with complete adjusters ------------------
COV <- data.frame(id = ID_T, grp = grp, purity = purity, leukocyte = leuk,
                  adipose = adipose, prolif = prolif, amp = amp,
                  stringsAsFactors = FALSE)
COV$myc <- ifelse(COV$grp == "MYC-high", 1, ifelse(COV$grp == "MYC-low", 0, NA))
keep <- !is.na(COV$myc) & !is.na(COV$purity) & !is.na(COV$leukocyte)
FIT  <- COV[keep, , drop = FALSE]
UF   <- UT[, FIT$id, drop = FALSE]
message("   analysis set: ", nrow(FIT), " samples with a MYC group AND complete ",
        "purity and leukocyte fraction")
message("      of ", sum(!is.na(COV$myc)), " in the two MYC groups -- ",
        sum(!is.na(COV$myc)) - nrow(FIT), " dropped for missing covariates")
coverage <- tibble::tibble(
  step = c("TCGA samples", "in a MYC tertile group", "with complete adjusters"),
  n = c(length(ID_T), sum(!is.na(COV$myc)), nrow(FIT)))
print(as.data.frame(coverage), row.names = FALSE)

# =============================================================================
# 2. GATE 0 - the residual machinery must reproduce E25
# =============================================================================
# mitoPPS is centred at exactly 1 per pathway, so ranking group MEANS and
# ranking group-mean RESIDUALS from an intercept-only fit are the same object.
# If this does not reproduce E25's percentiles exactly, the machinery is wrong
# and every adjusted number below is a rebuild rather than an adjustment.

message("\n2. GATE 0: intercept-only residuals must reproduce E25")

.residualise <- function(M, D) {
  X <- if (length(D)) cbind(1, as.matrix(D)) else matrix(1, ncol(M), 1L)
  R <- t(qr.resid(qr(X), t(M)))
  dimnames(R) <- dimnames(M)
  R
}

.group_pct <- function(R, g, panel) {
  gm <- rowMeans(R[panel, g, drop = FALSE])
  stats::setNames(.pct(gm), panel)
}

# E25 used the WHOLE cohort's group means. Reproduce on exactly those samples.
id_lo_all <- names(grp)[which(grp == "MYC-low")]
id_hi_all <- names(grp)[which(grp == "MYC-high")]
R_int_all <- .residualise(UT, NULL)
g0 <- tibble::tibble(
  group = c("MYC-low", "MYC-high"),
  e26_pct = c(.group_pct(R_int_all, id_lo_all, shared)[[PRIMARY]],
              .group_pct(R_int_all, id_hi_all, shared)[[PRIMARY]]),
  e25_pct = c(
    e25$positions$pct[e25$positions$panel == "shared (142)" &
                        e25$positions$pathway == PRIMARY &
                        e25$positions$arm == "TCGA MYC-low"],
    e25$positions$pct[e25$positions$panel == "shared (142)" &
                        e25$positions$pathway == PRIMARY &
                        e25$positions$arm == "TCGA MYC-high"]))
g0$delta <- g0$e26_pct - g0$e25_pct
print(as.data.frame(g0), digits = 6, row.names = FALSE)
GATE_0 <- max(abs(g0$delta)) < 1e-9
if (!GATE_0) {
  stop("GATE 0 FAILED: the residual machinery does not reproduce E25. ",
       "Every adjusted number below would be a rebuild, not an adjustment.",
       call. = FALSE)
}
message("   GATE 0 PASSED: max |delta| = ", format(max(abs(g0$delta)), digits = 2))

# =============================================================================
# 3. PART 0 - the exposure check. How much does composition differ at all?
# =============================================================================
# This is not an outcome. It sizes the problem: an adjuster cannot change much
# if the groups do not differ on it.

message("\n3. PART 0: the exposure check - do the MYC groups differ in composition?")

.grp_diff <- function(v) {
  hi <- v[FIT$myc == 1]; lo <- v[FIT$myc == 0]
  tt <- stats::t.test(hi, lo)
  c(high = mean(hi), low = mean(lo), diff = mean(hi) - mean(lo),
    lo_ci = unname(tt$conf.int[1]), hi_ci = unname(tt$conf.int[2]))
}
exposure <- dplyr::bind_rows(lapply(
  c(purity = "purity", leukocyte = "leukocyte", adipose = "adipose"),
  function(v) { z <- .grp_diff(FIT[[v]]); tibble::tibble(
    covariate = v, mean_high = z[["high"]], mean_low = z[["low"]],
    diff = z[["diff"]], ci_lo = z[["lo_ci"]], ci_hi = z[["hi_ci"]]) }),
  .id = "which")
exposure$which <- NULL
print(as.data.frame(exposure), digits = 3, row.names = FALSE)

tracer_rho <- c(
  purity    = stats::cor(FIT$adipose, FIT$purity, method = "spearman"),
  leukocyte = stats::cor(FIT$adipose, FIT$leukocyte, method = "spearman"))
message("   tracer against the named adjusters: rho(adipose, purity) ",
        sprintf("%+.3f", tracer_rho[["purity"]]),
        " | rho(adipose, leukocyte) ", sprintf("%+.3f", tracer_rho[["leukocyte"]]))
message("   READ THIS BEFORE THE GATES: the adjusters the plan named capture ",
        "the tracer only partly.")

# =============================================================================
# 4. GATE 1 - THE POSITIVE CONTROL, READ FIRST
# =============================================================================

message("\n4. GATE 1: THE POSITIVE CONTROL - read before any OXPHOS number")

.scalar_fit <- function(y, adj) {
  D <- if (length(adj)) FIT[, adj, drop = FALSE] else NULL
  X <- if (is.null(D)) data.frame(myc = FIT$myc) else data.frame(myc = FIT$myc, D)
  fit <- stats::lm(y ~ ., data = X)
  ci  <- stats::confint(fit, "myc")
  r2_adj <- if (length(adj)) {
    summary(stats::lm(y ~ ., data = as.data.frame(FIT[, adj, drop = FALSE])))$r.squared
  } else NA_real_
  c(beta = unname(stats::coef(fit)[["myc"]]),
    lo = unname(ci[1]), hi = unname(ci[2]), r2_adjusters = r2_adj)
}

# `skip` names models whose adjuster set CONTAINS the outcome. Fitting the
# tracer on itself returns zero by construction, which is arithmetic rather than
# evidence, so those rows are labelled instead of fitted.
.control_table <- function(y, label, skip = character(0)) {
  base <- .scalar_fit(y, character(0))
  dplyr::bind_rows(lapply(names(MODELS), function(m) {
    if (m %in% skip) {
      return(tibble::tibble(quantity = label, model = m, beta = NA_real_,
                            ci_lo = NA_real_, ci_hi = NA_real_,
                            retained = 0, r2_adjusters = NA_real_,
                            clear_of_zero = NA,
                            note = "zero BY CONSTRUCTION - the outcome is an adjuster"))
    }
    z <- .scalar_fit(y, MODELS[[m]])
    tibble::tibble(quantity = label, model = m, beta = z[["beta"]],
                   ci_lo = z[["lo"]], ci_hi = z[["hi"]],
                   retained = z[["beta"]] / base[["beta"]],
                   r2_adjusters = z[["r2_adjusters"]],
                   clear_of_zero = (z[["lo"]] > 0 & z[["hi"]] > 0) |
                                   (z[["lo"]] < 0 & z[["hi"]] < 0),
                   note = "")
  }))
}

ADIPO_MODELS <- names(MODELS)[vapply(MODELS, function(a) "adipose" %in% a, TRUE)]
control <- .control_table(FIT$prolif,  "POSITIVE CONTROL: proliferation")
tracer  <- .control_table(FIT$adipose, "TRACER: adipocyte markers",
                          skip = ADIPO_MODELS)

message("   the positive control - MYC-high minus MYC-low in proliferation:")
print(as.data.frame(control[, -1]), digits = 3, row.names = FALSE)

ctrl_main <- control[control$model == MODEL_MAIN, ]
GATE_1 <- ctrl_main$retained >= CTRL_RETAIN &&
          ctrl_main$clear_of_zero &&
          sign(ctrl_main$beta) == sign(control$beta[control$model == "M0 unadjusted"])
message("\n   the adjusters explain ",
        sprintf("%.1f", 100 * ctrl_main$r2_adjusters),
        " pct of the control's variance ",
        "(the mouse fat-pad figure that triggered trap 15 was 15 pct)")
message("   GATE 1: ", if (GATE_1) {
  paste0("PASSED - the control keeps ", sprintf("%.0f", 100 * ctrl_main$retained),
         " pct of its unadjusted difference, same sign, interval clear of zero")
} else {
  paste0("*** FAILED - the control keeps only ",
         sprintf("%.0f", 100 * ctrl_main$retained),
         " pct. The adjustment is removing biology and NOTHING below is readable")
})

message("\n   the composition tracer - how much traction each adjustment has:")
print(as.data.frame(tracer[, -1]), digits = 3, row.names = FALSE)
trc_main <- tracer[tracer$model == MODEL_MAIN, ]
GATE_2_INADEQUATE <- abs(trc_main$retained) >= TRACER_LEFT
message("   GATE 2: ", if (GATE_2_INADEQUATE) {
  paste0("*** ", MODEL_MAIN, " leaves ", sprintf("%.0f", 100 * abs(trc_main$retained)),
         " pct of the composition difference standing. It is NOT an adequate",
         " composition control; the adipose models are primary")
} else {
  paste0(MODEL_MAIN, " removes ", sprintf("%.0f", 100 * (1 - abs(trc_main$retained))),
         " pct of the composition difference - adequate on its own terms")
})
MODEL_READ <- if (GATE_2_INADEQUATE) "M3 purity + leuk + adipose" else MODEL_MAIN
message("   the primary adjusted reading is therefore: ", MODEL_READ)

# =============================================================================
# 5. PART A - the adjusted contrast for every pathway
# =============================================================================

message("\n5. PART A: the adjusted MYC contrast, all ", length(shared), " pathways")

.panel_fit <- function(adj) {
  D <- if (length(adj)) as.matrix(FIT[, adj, drop = FALSE]) else NULL
  X <- cbind(1, myc = FIT$myc, D)
  Q <- qr(X)
  B <- t(qr.coef(Q, t(UF)))            # pathways x coefficients
  stats::setNames(B[, "myc"], rownames(UF))
}
beta <- lapply(MODELS, .panel_fit)
names(beta) <- names(MODELS)
stopifnot(all(vapply(beta, function(b) identical(names(b), shared), TRUE)))

# NOTE: the vector is hoisted OUT of tibble(). Inside tibble() the new `beta`
# column masks the outer `beta` list, so a second `beta[[m]]` would index a
# numeric vector by a string and fail. Same class as the E20 `says` bug.
contrasts_long <- dplyr::bind_rows(lapply(names(beta), function(m) {
  b <- unname(beta[[m]])
  tibble::tibble(model = m, pathway = shared, beta = b, pct = .pct(b))
}))

# attenuation of every pathway under every adjustment
attenuation <- dplyr::bind_rows(lapply(setdiff(names(beta), "M0 unadjusted"),
  function(m) {
    r <- unname(beta[[m]] / beta[["M0 unadjusted"]])
    tibble::tibble(model = m, pathway = shared, retained = r)
  }))

message("   attenuation across the panel, by model:")
attenuation %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    median_retained = stats::median(retained),
    q05 = unname(stats::quantile(retained, PANEL_TAIL)),
    q95 = unname(stats::quantile(retained, 1 - PANEL_TAIL)), .groups = "drop") %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

# =============================================================================
# 6. PART B - R1: the primary's attenuation against the panel's
# =============================================================================

message("\n6. PART B: R1 - is the primary attenuated MORE than the panel?")

r1 <- dplyr::bind_rows(lapply(setdiff(names(beta), "M0 unadjusted"), function(m) {
  a <- attenuation$retained[attenuation$model == m]
  names(a) <- attenuation$pathway[attenuation$model == m]
  dplyr::bind_rows(lapply(READOUTS, function(p) {
    tibble::tibble(model = m, pathway = p, retained = unname(a[[p]]),
                   panel_median = stats::median(a),
                   panel_pctile = 100 * mean(a <= a[[p]]),
                   outlier = a[[p]] < stats::quantile(a, PANEL_TAIL) |
                             a[[p]] > stats::quantile(a, 1 - PANEL_TAIL))
  }))
}))
print(as.data.frame(r1), digits = 3, row.names = FALSE)

# =============================================================================
# 7. PART C - the group orderings, adjusted
# =============================================================================

message("\n7. PART C: where the primary sits, adjusted and unadjusted")

positions <- dplyr::bind_rows(lapply(names(MODELS), function(m) {
  R <- .residualise(UF, if (length(MODELS[[m]])) FIT[, MODELS[[m]], drop = FALSE] else NULL)
  dplyr::bind_rows(lapply(c("MYC-low", "MYC-high"), function(g) {
    ids <- FIT$id[FIT$grp == g]
    v   <- .group_pct(R, ids, shared)
    bs <- vapply(seq_len(NBOOT), function(i) {
      kk <- sample(ids, replace = TRUE)
      .group_pct(R, kk, shared)[[PRIMARY]] }, 0)
    tibble::tibble(model = m, group = g, pathway = PRIMARY, pct = v[[PRIMARY]],
                   ci_lo = unname(stats::quantile(bs, 0.025)),
                   ci_hi = unname(stats::quantile(bs, 0.975)))
  }))
}))
print(as.data.frame(positions[, -3]), digits = 3, row.names = FALSE)

# THE COVERAGE EFFECT, surfaced rather than left to be mistaken for a
# discrepancy. E25's percentiles are on all 730 samples of the two MYC groups;
# every model here is on the 675 with complete adjusters. M0 is the SAME
# estimand on a SMALLER set, so any difference between these two rows is the 55
# dropped samples and nothing else.
coverage_effect <- tibble::tibble(
  set = c("all 730 (E25)", "all 730 (E25)",
          paste0("covariate-complete ", nrow(FIT), " (E26 M0)"),
          paste0("covariate-complete ", nrow(FIT), " (E26 M0)")),
  group = rep(c("MYC-low", "MYC-high"), 2),
  pct = c(g0$e25_pct[g0$group == "MYC-low"], g0$e25_pct[g0$group == "MYC-high"],
          positions$pct[positions$model == "M0 unadjusted" & positions$group == "MYC-low"],
          positions$pct[positions$model == "M0 unadjusted" & positions$group == "MYC-high"]))
message("\n   THE COVERAGE EFFECT - E25's full 730 against this script's ",
        nrow(FIT), " with complete adjusters:")
print(as.data.frame(coverage_effect), digits = 3, row.names = FALSE)
message("   the same estimand on a smaller set. Any difference here is the ",
        sum(!is.na(COV$myc)) - nrow(FIT), " dropped samples,")
message("   NOT an adjustment effect. Every adjusted row below is compared ",
        "against E26's M0, never against E25's.")

separation <- positions %>%
  dplyr::select(model, group, pct) %>%
  tidyr::pivot_wider(names_from = group, values_from = pct) %>%
  dplyr::mutate(separation = `MYC-high` - `MYC-low`)
separation$retained <-
  separation$separation / separation$separation[separation$model == "M0 unadjusted"]
message("\n   the separation, and how much of it each adjustment leaves:")
print(as.data.frame(separation), digits = 3, row.names = FALSE)

# =============================================================================
# 8. PART D - the high-purity subset
# =============================================================================

message("\n8. PART D: the high-purity subset (purity >= ", PURITY_HI, ")")

hp <- FIT$purity >= PURITY_HI
hp_n <- tibble::tibble(group = c("MYC-low", "MYC-high"),
                       n = c(sum(hp & FIT$myc == 0), sum(hp & FIT$myc == 1)))
print(as.data.frame(hp_n), row.names = FALSE)

R_hp <- .residualise(UF[, hp, drop = FALSE], NULL)
hp_pos <- dplyr::bind_rows(lapply(c("MYC-low", "MYC-high"), function(g) {
  ids <- FIT$id[hp & FIT$grp == g]
  v <- .group_pct(R_hp, ids, shared)
  bs <- vapply(seq_len(NBOOT), function(i)
    .group_pct(R_hp, sample(ids, replace = TRUE), shared)[[PRIMARY]], 0)
  tibble::tibble(subset = paste0("purity >= ", PURITY_HI), group = g,
                 pct = v[[PRIMARY]],
                 ci_lo = unname(stats::quantile(bs, 0.025)),
                 ci_hi = unname(stats::quantile(bs, 0.975)))
}))
print(as.data.frame(hp_pos), digits = 3, row.names = FALSE)
HP_SEP <- hp_pos$pct[hp_pos$group == "MYC-high"] - hp_pos$pct[hp_pos$group == "MYC-low"]
message("   separation in the high-purity subset: ", sprintf("%.1f", HP_SEP),
        " against ", sprintf("%.1f", separation$separation[separation$model == "M0 unadjusted"]),
        " unadjusted in the full analysis set")

# --- the DNA-level split, which composition cannot manufacture ----------------
message("\n   the DNA-level split, adjusted - a copy-number call, not expression:")
amp_fit <- COV[!is.na(COV$amp) & !is.na(COV$purity) & !is.na(COV$leukocyte), ]
UA <- UT[, amp_fit$id, drop = FALSE]
amp_pos <- dplyr::bind_rows(lapply(names(MODELS), function(m) {
  D <- if (length(MODELS[[m]])) amp_fit[, MODELS[[m]], drop = FALSE] else NULL
  R <- .residualise(UA, D)
  dplyr::bind_rows(lapply(c(FALSE, TRUE), function(a) {
    ids <- amp_fit$id[amp_fit$amp == a]
    tibble::tibble(model = m, group = if (a) "MYC-amplified" else "not amplified",
                   n = length(ids), pct = .group_pct(R, ids, shared)[[PRIMARY]])
  }))
}))
print(as.data.frame(amp_pos), digits = 3, row.names = FALSE)

# =============================================================================
# 9. PART E - R2: does the cross-species agreement survive?
# =============================================================================

message("\n9. PART E: R2 - the adjusted contrast against the mouse")

ms <- mo$mitopps_scores
MPATHS <- setdiff(names(ms), c("sample", "group", "timepoint", "myc_status"))
MM <- as.matrix(ms[, MPATHS, drop = FALSE]); rownames(MM) <- ms$sample
MOUSE_CONTRAST <- list(
  Myc_effect_6W  = c(a = "6W_neg",  b = "6W_pos"),
  Myc_effect_12W = c(a = "12W_neg", b = "12W_pos"),
  `Temporal_Myc-` = c(a = "6W_neg", b = "12W_neg"),
  `Temporal_Myc+` = c(a = "6W_pos", b = "12W_pos"))
mouse_diff <- lapply(MOUSE_CONTRAST, function(cn)
  (colMeans(MM[ms$group == cn[["b"]], shared, drop = FALSE]) -
   colMeans(MM[ms$group == cn[["a"]], shared, drop = FALSE])))
CONTRAST_KIND <- c(Myc_effect_6W = "genotype, within timepoint",
                   Myc_effect_12W = "genotype, within timepoint",
                   `Temporal_Myc-` = "TIMELINE - batch-confounded",
                   `Temporal_Myc+` = "TIMELINE - batch-confounded")

.rho_ci <- function(x, y) {
  r <- stats::cor(x, y, method = "spearman")
  bs <- vapply(seq_len(NBOOT), function(i) {
    kk <- sample.int(length(x), replace = TRUE)
    if (length(unique(y[kk])) < 3L) return(NA_real_)
    stats::cor(x[kk], y[kk], method = "spearman") }, 0)
  c(rho = r, lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
    hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)))
}

cross <- dplyr::bind_rows(lapply(names(MODELS), function(m) {
  dplyr::bind_rows(lapply(names(MOUSE_CONTRAST), function(mc) {
    z <- .rho_ci(unname(beta[[m]]), unname(mouse_diff[[mc]]))
    tibble::tibble(model = m, mouse = mc, kind = unname(CONTRAST_KIND[[mc]]),
                   rho = z[["rho"]], ci_lo = z[["lo"]], ci_hi = z[["hi"]],
                   agrees = z[["rho"]] > 0 & z[["lo"]] > 0)
  }))
}))
print(as.data.frame(cross), digits = 2, row.names = FALSE)

.survives <- function(m) {
  d <- cross[cross$model == m & cross$mouse %in%
               c("Myc_effect_6W", "Myc_effect_12W"), ]
  all(d$agrees)
}
SURVIVES <- vapply(names(MODELS), .survives, TRUE)
message("\n   agreement with the mouse Myc effect at BOTH timepoints:")
for (m in names(MODELS)) message("      ", sprintf("%-28s %s", m,
  if (SURVIVES[[m]]) "SURVIVES" else "does NOT survive"))

# =============================================================================
# 10. PART F - the figure. Insight only, no significance marks.
# =============================================================================

message("\n10. PART F: figure")
if (!dir.exists(DIR_E26_FIG)) dir.create(DIR_E26_FIG, recursive = TRUE)

f_dat <- dplyr::bind_rows(
  positions %>% dplyr::transmute(model, group, pct, ci_lo, ci_hi,
                                 panel = "primary percentile"),
  hp_pos %>% dplyr::transmute(model = paste0("purity >= ", PURITY_HI),
                              group, pct, ci_lo, ci_hi,
                              panel = "primary percentile")) %>%
  dplyr::mutate(model = factor(model, levels = c(names(MODELS),
                                                 paste0("purity >= ", PURITY_HI))))

p1 <- ggplot2::ggplot(f_dat, ggplot2::aes(model, pct, colour = group)) +
  ggplot2::geom_hline(yintercept = 50, linewidth = 0.3, colour = "grey70",
                      linetype = 2) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                          position = ggplot2::position_dodge(width = 0.4),
                          linewidth = 0.5) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.4), size = 2.4) +
  ggplot2::scale_colour_manual(values = c(`MYC-low` = "#0072B2",
                                          `MYC-high` = "#D55E00")) +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::labs(
    x = NULL, y = paste0(PRIMARY, " percentile in the group ordering"),
    colour = NULL,
    title = "Does the MYC separation survive composition adjustment?",
    subtitle = paste0("TCGA only. Positive control (proliferation) keeps ",
                      sprintf("%.0f", 100 * ctrl_main$retained),
                      " pct of its unadjusted difference under ", MODEL_MAIN, "."),
    caption = paste(
      "Exploratory, post-hoc, descriptive; for insight, not statistical evaluation.",
      "SCAN-B has no purity estimate and nothing is imputed, so this cannot replicate.",
      sep = "\n")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                 legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))

ggplot2::ggsave(file.path(DIR_E26_FIG, "E26_composition_adjustment.pdf"), p1,
                width = 7, height = 5)
message("   wrote E26_composition_adjustment.pdf")

# =============================================================================
# 11. PART G - the verdict
# =============================================================================

message("\n11. PART G: verdict")

sep_read <- separation$retained[separation$model == MODEL_READ]
VERDICT <- if (!GATE_1) {
  "UNREADABLE - the positive control did not survive; the adjustment is over-correcting"
} else if (SURVIVES[[MODEL_READ]] && sep_read >= CTRL_RETAIN) {
  "SURVIVES - the separation and the cross-species agreement both hold under adjustment"
} else if (!SURVIVES[[MODEL_READ]]) {
  "COMPOSITION ARTEFACT - the control survived but the cross-species agreement did not"
} else {
  "ATTENUATED - the agreement holds but the separation is substantially reduced"
}

verdict <- tibble::tibble(
  rule = c("GATE 0 machinery", "GATE 1 positive control", "GATE 2 tracer",
           "R1 primary vs panel attenuation", "R2 cross-species agreement",
           "high-purity subset", "VERDICT"),
  outcome = c(
    paste0("PASSED - intercept-only residuals reproduce E25 exactly (max |delta| ",
           format(max(abs(g0$delta)), digits = 2), ")"),
    paste0(if (GATE_1) "PASSED" else "FAILED", " - proliferation keeps ",
           sprintf("%.0f", 100 * ctrl_main$retained), " pct of its MYC difference under ",
           MODEL_MAIN, "; those adjusters explain ",
           sprintf("%.1f", 100 * ctrl_main$r2_adjusters), " pct of its variance"),
    paste0(MODEL_MAIN, " leaves ", sprintf("%.0f", 100 * abs(trc_main$retained)),
           " pct of the composition difference; primary adjusted model is ", MODEL_READ),
    paste0("the primary retains ",
           sprintf("%.0f", 100 * r1$retained[r1$model == MODEL_READ & r1$pathway == PRIMARY]),
           " pct against a panel median of ",
           sprintf("%.0f", 100 * r1$panel_median[r1$model == MODEL_READ & r1$pathway == PRIMARY]),
           " pct - panel percentile ",
           sprintf("%.0f", r1$panel_pctile[r1$model == MODEL_READ & r1$pathway == PRIMARY]),
           ", outlier: ", r1$outlier[r1$model == MODEL_READ & r1$pathway == PRIMARY]),
    paste0(MODEL_READ, ": ", if (SURVIVES[[MODEL_READ]]) "SURVIVES" else "does NOT survive",
           " agreement with the mouse Myc effect at both timepoints"),
    paste0("separation ", sprintf("%.1f", HP_SEP), " percentile points at n = ",
           sum(hp), ", against ",
           sprintf("%.1f", separation$separation[separation$model == "M0 unadjusted"]),
           " unadjusted"),
    VERDICT))
print(as.data.frame(verdict), row.names = FALSE, right = FALSE)
message("\n   AND WHATEVER IT SAYS, IT IS DESCRIPTIVE. NOT A TEST.")

NOTES <- c(
  paste("TCGA ONLY. SCAN-B has no purity estimate and nothing is imputed",
        "(trap 2), so this arm structurally cannot replicate. SCAN-B's E25",
        "numbers stay unadjusted and uncheckable, and that is a limitation of",
        "the data rather than a choice."),
  paste("Trap 15 governs the script: the positive control (proliferation) is",
        "read before any OXPHOS number and is a gate. The mouse fat-pad figure",
        "that created the trap was an adjustment destroying an association",
        "while explaining 15 pct of its variance."),
  paste("The adipocyte tracer is defined IN SCRIPT and is not a pinned",
        "snapshot. It is a covariate, never a scored set and never a",
        "claim-bearing object. It cannot separate 'less adipose tissue in the",
        "section' from 'tumour cells expressing fewer adipocyte genes', though",
        "in bulk breast tissue these markers are overwhelmingly adipocyte-derived."),
  paste("Attenuation is read against the panel's own attenuation distribution.",
        "Every adjustment shrinks every coefficient; only shrinking MORE than",
        "the panel is evidence of composition acting on the primary."),
  paste("Genome doublings are not a covariate here - the author ruled that out.",
        "Ploidy is available and deliberately unused."),
  paste("The adjusted models run on the", nrow(FIT), "samples with complete",
        "purity and leukocyte fraction, not E25's full 730. M0 is the same",
        "estimand on that smaller set and is the ONLY baseline the adjusted",
        "rows are compared against; the difference between E25's percentiles",
        "and E26's M0 is the dropped samples, not an adjustment effect."),
  paste("GATE 0 proves the residual machinery against E25: mitoPPS is centred",
        "at exactly 1 per pathway, so an intercept-only residual ordering and",
        "E25's group-mean ordering are the same object, and they match exactly."),
  "RANK ONLY across the species boundary. Mouse mitoPPS is not recomputed.",
  "Exploratory and post-hoc. The figure is for insight, not statistical evaluation.")

saveRDS(list(
  rules = RULES,
  settings = list(question = RULES$question, primary = PRIMARY,
                  readouts = READOUTS, tertile = TERTILE,
                  purity_hi = PURITY_HI, ctrl_retain = CTRL_RETAIN,
                  tracer_left = TRACER_LEFT, panel_tail = PANEL_TAIL,
                  models = MODELS, model_main = MODEL_MAIN,
                  adipo_markers = ADIPO_MARKERS, prolif = PROLIF,
                  nboot = NBOOT, seed = 1,
                  mouse_artefact = list(md5 = MOUSE_MD5, recomputed = FALSE)),
  coverage = coverage, coverage_effect = coverage_effect,
  gate_0 = g0, exposure = exposure, tracer_rho = tracer_rho,
  control = control, tracer = tracer, contrasts = contrasts_long,
  attenuation = attenuation, r1 = r1, positions = positions,
  separation = separation, high_purity = hp_pos, high_purity_n = hp_n,
  amp_positions = amp_pos, cross_species = cross, survives = SURVIVES,
  gate_1 = GATE_1, gate_2_inadequate = GATE_2_INADEQUATE,
  model_read = MODEL_READ, verdict = verdict, verdict_text = VERDICT,
  analysis_date = Sys.Date(), notes = NOTES), PATH_E26)

readr::write_csv(contrasts_long, PATH_E26_TAB)

message("\nE26: done.")
message("    ", PATH_E26)
message("    ", PATH_E26_TAB)
message("    ", file.path(DIR_E26_FIG, "E26_composition_adjustment.pdf"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E26)

  ## read in this order. the control is FIRST and it is a gate.
  x$gate_0 |> as.data.frame()
  x$control |> as.data.frame()
  x$tracer  |> as.data.frame()
  x$gate_1; x$gate_2_inadequate; x$model_read

  ## how much composition differs between the groups at all
  x$exposure |> as.data.frame()
  x$tracer_rho

  ## the verdict table, and the one-line verdict
  x$verdict |> as.data.frame()
  x$verdict_text

  ## R1: is the primary attenuated more than the panel is?
  subset(x$r1, pathway == "OXPHOS subunits") |> as.data.frame()

  ## the coverage effect - read this BEFORE comparing anything to E25
  x$coverage |> as.data.frame()
  x$coverage_effect |> as.data.frame()

  ## the primary's position under every model, and the separation it leaves
  x$positions  |> as.data.frame()
  x$separation |> as.data.frame()

  ## the high-purity subset and the DNA-level split
  x$high_purity_n |> as.data.frame()
  x$high_purity   |> as.data.frame()
  x$amp_positions |> as.data.frame()

  ## R2: does the cross-species agreement survive?
  x$cross_species |> as.data.frame()
  x$survives

  ## the panel-wide attenuation distribution the primary is read against
  tapply(x$attenuation$retained, x$attenuation$model, summary)

  ## every caveat the note must repeat
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")
}
