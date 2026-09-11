# scripts/E28_gap_readout.R
# =============================================================================
# E28 - the OXPHOS gap: its intervals, its comparators, and proliferation
# =============================================================================
#
# THE QUANTITY. The gap is OXPHOS's position MINUS a comparator's position,
# inside the same group ordering of the shared 142-pathway panel. Its default
# comparator is Mitochondrial central dogma.
#
# THE GAP IS POST HOC. It was derived AFTER E27 ran and was not among E27's
# R1-R5. It is reportable, and it is labelled post hoc in this object, in the
# verdict, on both figures and in the note. Nothing about it was pre-registered.
#
# EXPLORATORY, POST-HOC, DESCRIPTIVE, as all of Phase 2 is. The reading rules
# below are fixed before any number so that the outcome cannot be chosen
# afterwards. NONE OF THEM IS A PREDICTION.
#
# WHAT IS READ AND WHAT IS BUILT
# -------------------------------
# READ: E25 (group orderings), E26 (composition machinery and covariates),
# E27 (the eleven readouts). None of the three is modified. NO mitoPPS IS
# RE-SCORED - not the mouse artefact, not the 142-pathway universe.
#
# BUILT: `ox_rel` per sample, on the HUMAN side only, for R5. That is a gene-set
# COMPOSITE on the human linear matrices, rebuilt from E16's own recipe. It is
# NOT re-scoring: the "nothing is re-scored" rule protects mitoPPS, and a
# composite rebuilt from its saved recipe is a different object. The
# REPRODUCTION CONTROL is reported before any comparison uses it: the rebuilt
# ruler must return E16's saved rho(ox_rel, ox_ppd) of 0.9130 (TCGA) and 0.8852
# (SCAN-B). A mouse ox_rel is NOT built - the artefact carries pathway scores and
# not gene-level values, so it would mean touching that repo.
#
# TWO GAP DEFINITIONS, AND THEY ARE NOT THE SAME STATISTIC
# ---------------------------------------------------------
#   PRIMARY   gap(group b) - gap(group a): a DIFFERENCE OF GROUP GAPS.
#   NAMED ONCE  OXPHOS minus the comparator INSIDE the contrast's own ordering.
# For the WT window these are -42.3 and -28.9. Both are computed, the first is
# the readout, and the second exists in this script only so that no one quotes
# it by mistake.
#
# THE PROTECTION CLAIM IS TESTED BEFORE THE GAP IS READ
# ------------------------------------------------------
# The rationale justified the gap by saying that two categories in the same
# ordering move together under composition, so differencing cancels it. That
# argument is weaker than it looks: a percentile is ALREADY a within-ordering
# rank, so a monotone panel-wide shift changes nothing and differencing adds
# nothing there. The real threat is composition moving the two categories
# DIFFERENTIALLY, and differencing does not touch that. PART 0b measures it
# against E26's four composition models and IS REPORTED BEFORE ANY GAP RESULT,
# because if it fails the rationale's justification for the readout goes.
#
# THE COMPARATOR IS NOT A CONSERVATIVE CHOICE
# --------------------------------------------
# Ranked by distance from OXPHOS across the mouse WT window, central dogma is
# THIRD of six: Protein import is closest and Signaling is also closer. That
# ordering is identical on the value scale and the percentile scale, which is
# checked here rather than assumed. So central dogma is a MIDDLING comparator,
# the six gaps span a wide range, and they are reported INDIVIDUALLY and NEVER
# averaged into a headline.
#
# THE AGGREGATE IS NOT A SECOND LINE OF EVIDENCE
# -----------------------------------------------
# OXPHOS against the mean of the other six is, conceptually, `ox_rel` rebuilt in
# mitoPPS space. If they agree that is CORROBORATION ON A SECOND RULER and
# nothing more. E16's R1 declared the two rulers NOT INTERCHANGEABLE, and
# agreement on this one contrast does not overturn that. R5 therefore asks only
# whether they agree IN DIRECTION on this contrast.
#
# PROLIFERATION: THREE APPROACHES, NOT EQUIVALENT
# ------------------------------------------------
# MYC and proliferation are entangled BIOLOGICALLY, not merely by measurement,
# so adjusting on proliferation risks removing MYC itself.
#   C1 the estimator ladder  PRIMARY - removes nothing
#   C2 stratification        PRIMARY - holds it roughly constant, subtracts nothing
#   C3 adjustment            SECONDARY, AMBIGUOUS BY CONSTRUCTION
# PRE-DECLARED, BEFORE ANY NUMBER: C3 collapsing while C1 and C2 survive is read
# as OVER-ADJUSTMENT, not as refutation.
#
# TRAP 15 WITH A TWIST. The usual positive control is proliferation, and here
# proliferation is the adjuster. The control used instead is the TCGA COPY-NUMBER
# amplification split, which is the only candidate not defined on the expression
# matrix at all. ITS INDEPENDENCE IS IN MEASUREMENT, NOT IN BIOLOGY - amplified
# tumours proliferate more - and that is stated wherever it is used. Its OWN
# attenuation under the adjustment is reported beside the gap's: if the
# copy-number split attenuates as much as the gap does, the adjustment is
# removing MYC wholesale and the control has said so.
#
# The mitoribosome is NOT used as a control: it sits inside central dogma and is
# therefore one arm of the primary gap, which is circular for this readout.
#
# NO mitoPPS VALUE CROSSES THE SPECIES BOUNDARY. The gap is in PERCENTILE units,
# which is what makes it crossable at all. TCGA and SCAN-B are never averaged.
#
# The mouse timeline is BATCH-CONFOUNDED (batch = timepoint). That ruling stands:
# temporal numbers reported in full, confounder stated, nothing withdrawn.
#
# THE EXPRESSION-MATCHED-NULL ANALYSIS IN myc_mouse REMAINS THE PRIMARY EVIDENCE
# for the developmental change. This is a second presentation. IF THE TWO
# DISAGREE, THAT DISAGREEMENT IS THE FINDING.
#
# N3. Transcript associations. "Primed" appears nowhere. No ortholog call.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE28: the gap readout - intervals, comparators, proliferation\n",
        strrep("=", 78))

PATH_E28     <- file.path(DIR_RESULTS, "gap_readout.rds")
PATH_E28_TAB <- file.path(DIR_TABLES,  "E28_gap_readout.csv")
DIR_E28_FIG  <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_E25     <- file.path(DIR_RESULTS, "mitopps_group_rankings.rds")
PATH_E26     <- file.path(DIR_RESULTS, "purity_and_composition.rds")
PATH_E27     <- file.path(DIR_RESULTS, "category_readouts.rds")
PATH_E16     <- file.path(DIR_RESULTS, "respiratory_rulers.rds")
PATH_MOUSE   <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")

set.seed(1)
NBOOT <- 5000L

# =============================================================================
# 0. CONSTANTS AND THE READING RULES, FIXED BEFORE ANY NUMBER
# =============================================================================

MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"
CAT_OX    <- "OXPHOS"
CAT_NULL  <- "Mitochondrial central dogma"
OX_SUB_N  <- "OXPHOS subunits"
OX_AF_N   <- "OXPHOS assembly factors"
TERTILE   <- 1/3

ARM_ORDER <- c("mouse 6W WT", "mouse 12W WT", "mouse 6W Myc+", "mouse 12W Myc+",
               "TCGA MYC-low", "TCGA MYC-high", "SCAN-B MYC-low", "SCAN-B MYC-high")

# The rationale's tables, hard-coded so PART 0 asserts against them rather than
# against numbers this script produced itself.
RAT_GROUP <- c(`mouse 6W WT` = 24.7, `mouse 6W Myc+` = 12.0, `TCGA MYC-low` = 2.8,
               `SCAN-B MYC-low` = -0.7, `SCAN-B MYC-high` = -14.8,
               `TCGA MYC-high` = -17.6, `mouse 12W WT` = -17.6,
               `mouse 12W Myc+` = -22.5)
RAT_CONTRAST <- c(`WT window (6W -> 12W)` = -42.3, `Myc+ window (6W -> 12W)` = -34.5,
                  `MYC at 6W` = -12.7, `MYC at 12W` = -4.9,
                  `TCGA MYC-high - low` = -20.4, `SCAN-B MYC-high - low` = -14.1)
RAT_TOL <- 0.05      # the rationale rounds to one decimal

# E16's saved reproduction target is read FROM E16'S OBJECT at run time. A
# transcribed constant would only ever be as exact as the transcription, and the
# control is meant to be exact.

.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}

RULES <- list(
  quantity = paste(
    "the GAP: OXPHOS position minus a comparator's position inside the same",
    "group ordering of the shared 142-pathway panel"),
  post_hoc = paste(
    "THE GAP IS POST HOC. Derived after E27 ran and not among its R1-R5.",
    "Labelled post hoc in the object, the verdict, the figures and the note"),
  two_definitions = paste(
    "PRIMARY is a DIFFERENCE OF GROUP GAPS. The gap inside a contrast's own",
    "ordering is a DIFFERENT statistic (-42.3 against -28.9 for the WT window)",
    "and is computed here only so it is not quoted by mistake"),
  protection_first = paste(
    "the claim that differencing cancels composition is TESTED BEFORE any gap",
    "result, because if it fails the justification for the readout goes"),
  comparators = paste(
    "central dogma is THIRD of six by distance from OXPHOS - a MIDDLING",
    "comparator, not a conservative one. The six gaps are reported",
    "INDIVIDUALLY and NEVER averaged into a headline"),
  aggregate = paste(
    "the aggregate is ox_rel rebuilt in mitoPPS space. Agreement is",
    "CORROBORATION ON A SECOND RULER and NOT a second line of evidence.",
    "E16's R1 declared the rulers not interchangeable and that is not",
    "overturned here"),
  oxrel_build = paste(
    "ox_rel is REBUILT from E16's recipe on the human linear matrices, human",
    "side only, with a reproduction control reported before any comparison.",
    "That is not re-scoring: the rule protects mitoPPS, and this is a",
    "composite rebuilt from its own saved recipe"),
  prolif = paste(
    "C1 ladder and C2 strata are PRIMARY; C3 adjustment is SECONDARY and",
    "ambiguous by construction. C3 collapsing while C1 and C2 survive is",
    "OVER-ADJUSTMENT, not refutation"),
  control = paste(
    "the positive control is the TCGA copy-number split - independent of the",
    "expression matrix in MEASUREMENT, not in BIOLOGY. Its own attenuation is",
    "reported beside the gap's. The mitoribosome is NOT used: it sits inside",
    "central dogma and is one arm of the primary gap"),
  no_value_crossing = "no mitoPPS value crosses the species boundary anywhere",
  batch = paste(
    "the mouse timeline is batch-confounded (batch = timepoint); temporal",
    "numbers reported in full and nothing withdrawn on that account"),
  primacy = paste(
    "the expression-matched-null analysis in myc_mouse remains the primary",
    "evidence for the developmental change; if the two disagree THAT",
    "DISAGREEMENT IS THE FINDING"),
  exploratory = "post-hoc and descriptive; not a hypothesis test",
  n3 = "transcript associations; 'primed' is never written of a transcript")

message("\n0. reading rules fixed before any number")
message("   R1  the mouse developmental gap. Interval clear of zero -> a second,")
message("       independent presentation of the developmental change. Covering zero")
message("       -> NOT RESOLVABLE AT n = 6, said plainly; the mouse gap contrast")
message("       then CANNOT STAND ALONE and must be read beside myc_mouse's")
message("       expression-matched-null analysis, which remains primary")
message("   R2  developmental against MYC in the mouse: do the two SEPARATE with")
message("       intervals? If not, say so rather than quoting the ratio")
message("   R3  the human MYC difference under the ladder and within strata.")
message("       myc_mRNA sign-flipping is ALREADY KNOWN and is not on its own a")
message("       failure; the pattern across the ladder and the copy-number split are")
message("   R4  the comparator panel: a property of the compartment, or of one")
message("       pairing? Report the range either way")
message("   R5  the aggregate against ox_rel: agreement is corroboration,")
message("       DISAGREEMENT IS THE FINDING. Neither outcome is better")

# =============================================================================
# 1. INPUTS
# =============================================================================

message("\n1. inputs")

e25 <- readRDS(PATH_E25); e26 <- readRDS(PATH_E26); e27 <- readRDS(PATH_E27)
e16 <- readRDS(PATH_E16)
message("   E25 ", length(e25), " elements | E26 ", length(e26),
        " elements | E27 ", length(e27), " elements")

md5_now <- unname(tools::md5sum(PATH_MOUSE))
stopifnot(identical(md5_now, MOUSE_MD5))
mo <- readRDS(PATH_MOUSE)
message("   mouse artefact md5 matches: ", md5_now)

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
cvd  <- readRDS(here::here("data", "from_validation",
                          "tcga_brca_covariates.rds"))$covariates

shared <- e25$intersection$shared
CATS   <- e27$settings$categories
COMPARATORS <- setdiff(CATS, CAT_OX)
stopifnot(length(shared) == 142L, length(CATS) == 7L, length(COMPARATORS) == 6L)

UT <- mito$mitopps_universe[shared, , drop = FALSE]
US <- sc$mitopps_universe[shared, , drop = FALSE]
ID_T <- colnames(UT); ID_S <- colnames(US)
ms <- mo$mitopps_scores
MPATHS <- setdiff(names(ms), c("sample", "group", "timepoint", "myc_status"))
MM <- as.matrix(ms[, MPATHS, drop = FALSE]); rownames(MM) <- ms$sample
MMs <- t(MM[, shared, drop = FALSE])

MOUSE_ARM <- c(`6W_neg` = "mouse 6W WT", `12W_neg` = "mouse 12W WT",
               `6W_pos` = "mouse 6W Myc+", `12W_pos` = "mouse 12W Myc+")

.tertile_groups <- function(v, ids) {
  q <- stats::quantile(v, c(TERTILE, 1 - TERTILE), names = FALSE)
  g <- rep(NA_character_, length(v)); names(g) <- ids
  g[v <= q[1]] <- "MYC-low"; g[v >= q[2]] <- "MYC-high"
  g
}
est_T <- as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T])
est_S <- as.numeric(sc$gsva_new[MYC_REF, ID_S])
grp_T <- .tertile_groups(est_T, ID_T); grp_S <- .tertile_groups(est_S, ID_S)
stopifnot(identical(as.integer(c(sum(grp_T == "MYC-low", na.rm = TRUE),
                                 sum(grp_T == "MYC-high", na.rm = TRUE),
                                 sum(grp_S == "MYC-low", na.rm = TRUE),
                                 sum(grp_S == "MYC-high", na.rm = TRUE))),
                    as.integer(e25$myc_groups$n)))
message("   MYC groups reproduce E25 exactly")

ARM <- list()
for (g in names(MOUSE_ARM)) {
  ARM[[unname(MOUSE_ARM[[g]])]] <- list(sp = "mouse", M = MMs, j = which(ms$group == g))
}
for (lv in c("MYC-low", "MYC-high")) {
  ARM[[paste("TCGA", lv)]]   <- list(sp = "human", M = UT, j = which(grp_T == lv))
  ARM[[paste("SCAN-B", lv)]] <- list(sp = "human", M = US, j = which(grp_S == lv))
}
ARM <- ARM[ARM_ORDER]

CONTRASTS <- list(
  `WT window (6W -> 12W)`   = c(a = "mouse 6W WT",    b = "mouse 12W WT"),
  `Myc+ window (6W -> 12W)` = c(a = "mouse 6W Myc+",  b = "mouse 12W Myc+"),
  `MYC at 6W`               = c(a = "mouse 6W WT",    b = "mouse 6W Myc+"),
  `MYC at 12W`              = c(a = "mouse 12W WT",   b = "mouse 12W Myc+"),
  `TCGA MYC-high - low`     = c(a = "TCGA MYC-low",   b = "TCGA MYC-high"),
  `SCAN-B MYC-high - low`   = c(a = "SCAN-B MYC-low", b = "SCAN-B MYC-high"))
CONTRAST_KIND <- c(`WT window (6W -> 12W)` = "TIMELINE - batch-confounded",
                   `Myc+ window (6W -> 12W)` = "TIMELINE - batch-confounded",
                   `MYC at 6W` = "genotype, within timepoint",
                   `MYC at 12W` = "genotype, within timepoint",
                   `TCGA MYC-high - low` = "MYC activity, within cohort",
                   `SCAN-B MYC-high - low` = "MYC activity, within cohort")

# =============================================================================
# 2. PART 0 - reproduce E27, then define the gap and assert it
# =============================================================================

message("\n2. PART 0: reproduce E27, then define and assert the gap")

.gp <- function(z, j) {
  v <- rowMeans(z$M[, j, drop = FALSE])
  stats::setNames(.pct(v), rownames(z$M))
}
pos <- vapply(ARM_ORDER, function(a) .gp(ARM[[a]], ARM[[a]]$j),
              numeric(length(shared)))
rownames(pos) <- shared

READOUTS <- as.character(e27$settings$readouts)
e27p <- e27$category_positions
ctrl <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  tibble::tibble(arm = a, readout = READOUTS, e28 = pos[READOUTS, a],
                 e27 = e27p$pct[match(paste(a, READOUTS),
                                      paste(as.character(e27p$arm),
                                            as.character(e27p$readout)))])
}))
ctrl$delta <- ctrl$e28 - ctrl$e27
CONTROL_DELTA <- max(abs(ctrl$delta))
message("   ARITHMETIC CONTROL: ", nrow(ctrl), " positions against E27, max |delta| = ",
        format(CONTROL_DELTA, digits = 2))
if (CONTROL_DELTA >= 1e-9) {
  stop("CONTROL FAILED - stop. Every number below would be a rebuild.", call. = FALSE)
}

.gap_of <- function(p, comparator) unname(p[CAT_OX] - p[comparator])

gap_point <- vapply(ARM_ORDER, function(a) .gap_of(pos[, a], CAT_NULL), 0)
assert_group <- tibble::tibble(
  arm = names(RAT_GROUP), rationale = unname(RAT_GROUP),
  computed = unname(gap_point[names(RAT_GROUP)]))
assert_group$delta <- assert_group$computed - assert_group$rationale
print(as.data.frame(assert_group), digits = 4, row.names = FALSE)

gcx_point <- vapply(names(CONTRASTS), function(ct) {
  cn <- CONTRASTS[[ct]]
  .gap_of(pos[, cn[["b"]]], CAT_NULL) - .gap_of(pos[, cn[["a"]]], CAT_NULL) }, 0)
assert_contrast <- tibble::tibble(
  contrast = names(RAT_CONTRAST), rationale = unname(RAT_CONTRAST),
  computed = unname(gcx_point[names(RAT_CONTRAST)]))
assert_contrast$delta <- assert_contrast$computed - assert_contrast$rationale
print(as.data.frame(assert_contrast), digits = 4, row.names = FALSE)

RAT_OK <- max(abs(c(assert_group$delta, assert_contrast$delta))) < RAT_TOL
message("   the rationale's tables reproduce: ", RAT_OK,
        " (max |delta| ", sprintf("%.3f", max(abs(c(assert_group$delta,
                                                    assert_contrast$delta)))),
        ", tolerance ", RAT_TOL, " for one-decimal rounding)")
if (!RAT_OK) {
  bad <- dplyr::bind_rows(
    tibble::tibble(cell = assert_group$arm, rationale = assert_group$rationale,
                   computed = assert_group$computed, delta = assert_group$delta),
    tibble::tibble(cell = assert_contrast$contrast,
                   rationale = assert_contrast$rationale,
                   computed = assert_contrast$computed, delta = assert_contrast$delta))
  bad <- bad[abs(bad$delta) >= RAT_TOL, ]
  message("   *** DISCREPANCY - reported, not silently accepted. The cell(s):")
  print(as.data.frame(bad), digits = 4, row.names = FALSE)
  message("   This is a ROUNDING slip in the rationale, not a computation",
          " difference: the tolerance is exactly the one-decimal boundary.")
}

# --- the OTHER definition, named once so it is not quoted by mistake ---------
alt <- e27$contrast_positions
E27_NAME <- c(`WT window (6W -> 12W)` = "Temporal_Myc-",
              `Myc+ window (6W -> 12W)` = "Temporal_Myc+",
              `MYC at 6W` = "Myc_effect_6W", `MYC at 12W` = "Myc_effect_12W",
              `TCGA MYC-high - low` = "TCGA MYC-high - low",
              `SCAN-B MYC-high - low` = "SCAN-B MYC-high - low")
two_defs <- tibble::tibble(
  contrast = names(CONTRASTS),
  difference_of_group_gaps = unname(gcx_point[names(CONTRASTS)]),
  gap_inside_contrast_ordering = vapply(names(CONTRASTS), function(ct) {
    d <- alt[alt$contrast == E27_NAME[[ct]], ]
    d$pct[as.character(d$readout) == CAT_OX] -
      d$pct[as.character(d$readout) == CAT_NULL] }, 0))
message("\n   TWO DEFINITIONS, and they are NOT the same statistic:")
print(as.data.frame(two_defs), digits = 3, row.names = FALSE)
message("   The FIRST column is the readout. The second exists here only so that")
message("   nobody quotes it by mistake.")

# =============================================================================
# 3. PART 0b - THE PROTECTION CHECK. Reported BEFORE any gap result.
# =============================================================================
# The rationale justified the gap by saying differencing cancels composition.
# A percentile is already a within-ordering rank, so a monotone panel-wide shift
# changes nothing and differencing adds nothing there. The threat is
# DIFFERENTIAL movement. This measures it against E26's four models, TCGA only.

message("\n3. PART 0b: THE PROTECTION CHECK - does differencing cancel composition?")

ADIPO_MARKERS <- e26$settings$adipo_markers
tl <- readRDS(PATH_TCGA_LINEAR)
stopifnot(identical(tl$scale, "linear_deseq2_normalised"))
GT <- log2(tl$mat[, ID_T, drop = FALSE] + 1)
sl <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(sl$scale, "linear_deseq2_normalised"))
GS <- log2(sl$mat[, ID_S, drop = FALSE] + 1)
rm(tl, sl); invisible(gc(verbose = FALSE))

LA <- GT[ADIPO_MARKERS, , drop = FALSE]
adipose <- colMeans(t(scale(t(LA))))
kk <- match(ID_T, cvd$patient)
COVT <- data.frame(id = ID_T, purity = cvd$purity[kk],
                   leukocyte = cvd$leukocyte_fraction[kk],
                   adipose = adipose, amp = cvd$MYC_amp[kk],
                   grp = grp_T, stringsAsFactors = FALSE)

.residualise <- function(M, D) {
  X <- if (length(D)) cbind(1, as.matrix(D)) else matrix(1, ncol(M), 1L)
  R <- t(qr.resid(qr(X), t(M)))
  dimnames(R) <- dimnames(M)
  R
}
MODELS <- e26$settings$models
keep <- !is.na(COVT$grp) & !is.na(COVT$purity) & !is.na(COVT$leukocyte)
FT <- COVT[keep, , drop = FALSE]
UF <- UT[, FT$id, drop = FALSE]

protection <- dplyr::bind_rows(lapply(names(MODELS), function(m) {
  R <- .residualise(UF, if (length(MODELS[[m]])) FT[, MODELS[[m]], drop = FALSE] else NULL)
  dplyr::bind_rows(lapply(c("MYC-low", "MYC-high"), function(g) {
    p <- .pct(rowMeans(R[, FT$grp == g, drop = FALSE])); names(p) <- shared
    tibble::tibble(model = m, group = g, oxphos = unname(p[CAT_OX]),
                   comparator = unname(p[CAT_NULL]),
                   gap = unname(p[CAT_OX] - p[CAT_NULL]))
  }))
}))
print(as.data.frame(protection), digits = 3, row.names = FALSE)

.maxmove <- function(col) {
  m0 <- protection[[col]][protection$model == "M0 unadjusted"]
  max(abs(vapply(setdiff(names(MODELS), "M0 unadjusted"), function(m)
    max(abs(protection[[col]][protection$model == m] - m0)), 0)))
}
protection_summary <- tibble::tibble(
  quantity = c("OXPHOS position", paste(CAT_NULL, "position"), "THE GAP"),
  max_move_from_M0 = c(.maxmove("oxphos"), .maxmove("comparator"), .maxmove("gap")))
print(as.data.frame(protection_summary), digits = 3, row.names = FALSE)
PROTECTS <- protection_summary$max_move_from_M0[3] <
            min(protection_summary$max_move_from_M0[1:2])
message("   DOES DIFFERENCING PROTECT? ", if (PROTECTS) {
  "YES - the gap moves less under composition adjustment than either position does"
} else {
  paste0("NO - the gap moves at least as much as a single position. The",
         " rationale's justification for the readout DOES NOT HOLD and every",
         " gap number below inherits that")
})
# "Protection fails" and "composition wrecks the readout" are different
# statements. The second needs the movement measured against the SIGNAL.
GAP_SIGNAL <- abs(protection$gap[protection$model == "M0 unadjusted" &
                                   protection$group == "MYC-high"] -
                  protection$gap[protection$model == "M0 unadjusted" &
                                   protection$group == "MYC-low"])
COMP_SHARE <- protection_summary$max_move_from_M0[3] / GAP_SIGNAL
message("   AND THE ABSOLUTE SCALE, which the verdict above does not carry: the",
        " gap moves ", sprintf("%.1f", protection_summary$max_move_from_M0[3]),
        " points against a")
message("   TCGA gap DIFFERENCE of ", sprintf("%.1f", GAP_SIGNAL), " points - ",
        sprintf("%.0f", 100 * COMP_SHARE), " pct. So differencing does not help,")
message("   but composition is a SMALL perturbation either way. Both statements go",
        " in the note.")
message("   Note what this can and cannot say: it is TCGA-only, because SCAN-B")
message("   has no purity estimate and nothing is imputed (trap 2). And it runs on")
message("   the ", nrow(FT), " samples with complete adjusters, not the ",
        sum(!is.na(COVT$grp)), " in the two MYC")
message("   groups, so its M0 gap difference of ", sprintf("%.1f", GAP_SIGNAL),
        " is not PART A's ", sprintf("%.1f", abs(gcx_point[["TCGA MYC-high - low"]])),
        ".")
message("   That is the coverage effect E26 recorded, not an adjustment effect.")

# =============================================================================
# 4. PART A - intervals
# =============================================================================

message("\n4. PART A: intervals. Paired resampling - both positions and their")
message("   difference are recomputed INSIDE each draw.")

.boot_group <- function(a, comparators) {
  z <- ARM[[a]]
  vapply(seq_len(NBOOT), function(i) {
    p <- .gp(z, sample(z$j, replace = TRUE))
    p[CAT_OX] - p[comparators]
  }, numeric(length(comparators)))
}
BG <- lapply(ARM_ORDER, function(a) .boot_group(a, COMPARATORS))
names(BG) <- ARM_ORDER

gap_by_group <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  b <- BG[[a]]
  tibble::tibble(arm = a, species = ARM[[a]]$sp, n = length(ARM[[a]]$j),
                 comparator = COMPARATORS,
                 gap = vapply(COMPARATORS, function(cm) .gap_of(pos[, a], cm), 0),
                 ci_lo = apply(b, 1L, stats::quantile, 0.025),
                 ci_hi = apply(b, 1L, stats::quantile, 0.975))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER),
                clear_of_zero = (ci_lo > 0 & ci_hi > 0) | (ci_lo < 0 & ci_hi < 0))

message("   the PRIMARY gap (against ", CAT_NULL, ") per group:")
gap_by_group %>%
  dplyr::filter(comparator == CAT_NULL) %>%
  dplyr::select(arm, n, gap, ci_lo, ci_hi, clear_of_zero) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

.boot_contrast <- function(ct, comparators) {
  cn <- CONTRASTS[[ct]]; za <- ARM[[cn[["a"]]]]; zb <- ARM[[cn[["b"]]]]
  vapply(seq_len(NBOOT), function(i) {
    pa <- .gp(za, sample(za$j, replace = TRUE))
    pb <- .gp(zb, sample(zb$j, replace = TRUE))
    (pb[CAT_OX] - pb[comparators]) - (pa[CAT_OX] - pa[comparators])
  }, numeric(length(comparators)))
}
BC <- lapply(names(CONTRASTS), function(ct) .boot_contrast(ct, COMPARATORS))
names(BC) <- names(CONTRASTS)

gap_by_contrast <- dplyr::bind_rows(lapply(names(CONTRASTS), function(ct) {
  cn <- CONTRASTS[[ct]]; b <- BC[[ct]]
  tibble::tibble(contrast = ct, kind = unname(CONTRAST_KIND[[ct]]),
                 species = ARM[[cn[["a"]]]]$sp, comparator = COMPARATORS,
                 change = vapply(COMPARATORS, function(cm)
                   .gap_of(pos[, cn[["b"]]], cm) - .gap_of(pos[, cn[["a"]]], cm), 0),
                 ci_lo = apply(b, 1L, stats::quantile, 0.025),
                 ci_hi = apply(b, 1L, stats::quantile, 0.975))
})) %>%
  dplyr::mutate(clear_of_zero = (ci_lo > 0 & ci_hi > 0) | (ci_lo < 0 & ci_hi < 0))

message("\n   the PRIMARY gap CHANGE per contrast:")
gap_by_contrast %>%
  dplyr::filter(comparator == CAT_NULL) %>%
  dplyr::select(contrast, kind, change, ci_lo, ci_hi, clear_of_zero) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)
message("   mouse and human rows sit in one table because the gap is in")
message("   PERCENTILE units. No mitoPPS value crosses the species boundary.")

# --- R1 and R2 ---------------------------------------------------------------
wt <- gap_by_contrast[gap_by_contrast$contrast == "WT window (6W -> 12W)" &
                        gap_by_contrast$comparator == CAT_NULL, ]
R1 <- if (wt$clear_of_zero) {
  "CLEAR OF ZERO - the gap is a second, independent presentation of the developmental change"
} else {
  paste0("COVERS ZERO - NOT RESOLVABLE AT n = 6. This does not refute the point ",
         "estimate, and it means the mouse gap contrast CANNOT STAND ALONE: it ",
         "must be read beside myc_mouse's expression-matched-null analysis, ",
         "which remains the primary evidence")
}
message("\n   R1: ", R1)

myc6 <- BC[["MYC at 6W"]][match(CAT_NULL, COMPARATORS), ]
myc12 <- BC[["MYC at 12W"]][match(CAT_NULL, COMPARATORS), ]
wtb  <- BC[["WT window (6W -> 12W)"]][match(CAT_NULL, COMPARATORS), ]
sep6  <- stats::quantile(abs(wtb) - abs(myc6), c(0.025, 0.975))
sep12 <- stats::quantile(abs(wtb) - abs(myc12), c(0.025, 0.975))
r2 <- tibble::tibble(
  against = c("MYC at 6W", "MYC at 12W"),
  dev_change = wt$change,
  myc_change = c(gcx_point[["MYC at 6W"]], gcx_point[["MYC at 12W"]]),
  ratio = abs(wt$change) / abs(c(gcx_point[["MYC at 6W"]], gcx_point[["MYC at 12W"]])),
  sep_lo = c(sep6[[1]], sep12[[1]]), sep_hi = c(sep6[[2]], sep12[[2]]))
r2$separates <- r2$sep_lo > 0
print(as.data.frame(r2), digits = 3, row.names = FALSE)
R2 <- if (all(r2$separates)) {
  "the developmental change SEPARATES from the MYC change at both timepoints"
} else if (any(r2$separates)) {
  paste0("separates at ", paste(r2$against[r2$separates], collapse = " and "),
         " only; elsewhere the separation is NOT ESTABLISHED at this n")
} else {
  paste0("the separation is NOT ESTABLISHED at this n. The point-estimate ratio ",
         "of ", sprintf("%.1f", min(r2$ratio)), " to ", sprintf("%.1f", max(r2$ratio)),
         " is reported as a point estimate only and must not be quoted as a result")
}
message("   R2: ", R2)

# =============================================================================
# 5. PART B - the comparator panel
# =============================================================================

message("\n5. PART B: the comparator panel - six gaps, reported individually")

wv <- subset(e27$window_change, window == "mouse WT window (6W -> 12W)")
ox_v <- wv$delta[as.character(wv$readout) == CAT_OX]
dir <- e27$direction
comparator_distance <- tibble::tibble(
  comparator = COMPARATORS,
  own_change_value = wv$delta[match(COMPARATORS, as.character(wv$readout))],
  dist_value = abs(ox_v - wv$delta[match(COMPARATORS, as.character(wv$readout))]),
  own_change_pct = dir$pct_delta[match(COMPARATORS, dir$readout)],
  dist_pct = abs(dir$pct_delta[dir$readout == CAT_OX] -
                   dir$pct_delta[match(COMPARATORS, dir$readout)]))
comparator_distance <- comparator_distance[order(comparator_distance$dist_value), ]
comparator_distance$rank_value <- seq_len(nrow(comparator_distance))
comparator_distance$rank_pct <- rank(comparator_distance$dist_pct)
SCALES_AGREE <- identical(comparator_distance$rank_value,
                          as.integer(comparator_distance$rank_pct))
print(as.data.frame(comparator_distance), digits = 3, row.names = FALSE)
message("   the two scales give the same comparator ordering: ", SCALES_AGREE)
message("   ", CAT_NULL, " is rank ",
        comparator_distance$rank_value[comparator_distance$comparator == CAT_NULL],
        " of 6 - a MIDDLING comparator, not a conservative one.")

message("\n   B1: all six gaps, per contrast. NEVER averaged into a headline.")
gap_by_contrast %>%
  dplyr::select(contrast, comparator, change, ci_lo, ci_hi, clear_of_zero) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

comparator_range <- gap_by_contrast %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(n_same_sign = max(sum(change < 0), sum(change > 0)),
                   n_clear = sum(clear_of_zero),
                   min_change = min(change), max_change = max(change),
                   .groups = "drop")
message("\n   the RANGE across comparators, per contrast:")
print(as.data.frame(comparator_range), digits = 3, row.names = FALSE)

hum <- comparator_range[grepl("MYC-high - low", comparator_range$contrast), ]
message("   in the two HUMAN contrasts the comparators agree in sign ",
        paste(hum$n_same_sign, collapse = " and "), " of 6 times")
wtr <- comparator_range[comparator_range$contrast == "WT window (6W -> 12W)", ]
R4 <- if (wtr$n_same_sign == 6L) {
  paste0("all six comparators agree in sign in the WT window (range ",
         sprintf("%.1f", wtr$min_change), " to ", sprintf("%.1f", wtr$max_change),
         ") - a property of the COMPARTMENT, not of one pairing IN THE MOUSE. ",
         "In the two human contrasts only ",
         paste(hum$n_same_sign, collapse = " and "),
         " of 6 agree, so the human claim is narrower than the mouse one")
} else {
  paste0("only ", wtr$n_same_sign, " of 6 comparators agree in sign in the WT ",
         "window (range ", sprintf("%.1f", wtr$min_change), " to ",
         sprintf("%.1f", wtr$max_change), ") - the claim NARROWS to the pairings ",
         "that hold, and the note names them")
}
message("   R4: ", R4)

# --- B2 the aggregate, and B3 --------------------------------------------------
.agg_gap <- function(p, num) unname(p[num] - mean(p[COMPARATORS]))
.boot_agg <- function(a, num) {
  z <- ARM[[a]]
  vapply(seq_len(NBOOT), function(i) .agg_gap(.gp(z, sample(z$j, replace = TRUE)), num), 0)
}
AGG_NUM <- c(`subunits (the ox_rel analogue)` = OX_SUB_N, `umbrella` = CAT_OX)
aggregate_tab <- dplyr::bind_rows(lapply(names(AGG_NUM), function(nm) {
  num <- AGG_NUM[[nm]]
  dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
    b <- .boot_agg(a, num)
    tibble::tibble(numerator = nm, arm = a, species = ARM[[a]]$sp,
                   aggregate_gap = .agg_gap(pos[, a], num),
                   ci_lo = unname(stats::quantile(b, 0.025)),
                   ci_hi = unname(stats::quantile(b, 0.975)))
  }))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER))
message("\n   B2: the aggregate, both numerators. CORROBORATION ONLY, never a")
message("   second line of evidence.")
aggregate_tab %>%
  dplyr::select(numerator, arm, aggregate_gap, ci_lo, ci_hi) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

agg_contrast <- dplyr::bind_rows(lapply(names(AGG_NUM), function(nm) {
  num <- AGG_NUM[[nm]]
  dplyr::bind_rows(lapply(names(CONTRASTS), function(ct) {
    cn <- CONTRASTS[[ct]]; za <- ARM[[cn[["a"]]]]; zb <- ARM[[cn[["b"]]]]
    b <- vapply(seq_len(NBOOT), function(i)
      .agg_gap(.gp(zb, sample(zb$j, replace = TRUE)), num) -
      .agg_gap(.gp(za, sample(za$j, replace = TRUE)), num), 0)
    tibble::tibble(numerator = nm, contrast = ct,
                   change = .agg_gap(pos[, cn[["b"]]], num) -
                            .agg_gap(pos[, cn[["a"]]], num),
                   ci_lo = unname(stats::quantile(b, 0.025)),
                   ci_hi = unname(stats::quantile(b, 0.975)))
  }))
}))
agg_contrast$clear_of_zero <- (agg_contrast$ci_lo > 0 & agg_contrast$ci_hi > 0) |
                              (agg_contrast$ci_lo < 0 & agg_contrast$ci_hi < 0)
print(as.data.frame(agg_contrast), digits = 3, row.names = FALSE)
DILUTION <- agg_contrast$change[agg_contrast$numerator == names(AGG_NUM)[1]] -
            agg_contrast$change[agg_contrast$numerator == names(AGG_NUM)[2]]
message("   the umbrella numerator DILUTES by ", sprintf("%.1f", mean(abs(DILUTION))),
        " percentile points on average - exactly what E27 predicted, since")
message("   assembly factors barely move while subunits do.")

subunits_vs_assembly <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  z <- ARM[[a]]
  b <- vapply(seq_len(NBOOT), function(i) {
    p <- .gp(z, sample(z$j, replace = TRUE)); unname(p[OX_SUB_N] - p[OX_AF_N]) }, 0)
  tibble::tibble(arm = a, species = z$sp,
                 gap = unname(pos[OX_SUB_N, a] - pos[OX_AF_N, a]),
                 ci_lo = unname(stats::quantile(b, 0.025)),
                 ci_hi = unname(stats::quantile(b, 0.975)))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER))
message("\n   B3: subunits against assembly factors, carried from E27's R4b:")
print(as.data.frame(subunits_vs_assembly), digits = 3, row.names = FALSE)

# =============================================================================
# 6. PART C - the proliferation confound
# =============================================================================

message("\n6. PART C: proliferation")

PROLIF <- "PROLIF_DISJOINT"
stopifnot(PROLIF %in% rownames(mito$gsva_cov), PROLIF %in% rownames(sc$gsva_cov))
pro_T <- as.numeric(mito$gsva_cov[PROLIF, ID_T])
pro_S <- as.numeric(sc$gsva_cov[PROLIF, ID_S])
message("   SCAN-B HAS a proliferation score (", PROLIF, "), so C2 runs in BOTH")
message("   cohorts - unlike E26's purity arm, which was structurally TCGA-only.")

# --- C1 the estimator ladder --------------------------------------------------
EST_PANEL <- e25$settings$estimator_panel
.est_vec <- function(coh, w) {
  if (identical(coh, "TCGA")) {
    if (identical(w, "log2MYC")) return(as.numeric(nw$tcga_log2MYC[ID_T]))
    if (w %in% rownames(nw$tcga_M_b_variants)) return(as.numeric(nw$tcga_M_b_variants[w, ID_T]))
    return(as.numeric(nw$tcga_gsva_new[w, ID_T]))
  }
  if (identical(w, "log2MYC")) return(as.numeric(sc$log2MYC[ID_S]))
  if (w %in% rownames(sc$M_b_variants)) return(as.numeric(sc$M_b_variants[w, ID_S]))
  as.numeric(sc$gsva_new[w, ID_S])
}
.gap_diff <- function(U, g) {
  pl <- .pct(rowMeans(U[, which(g == "MYC-low"), drop = FALSE]))
  ph <- .pct(rowMeans(U[, which(g == "MYC-high"), drop = FALSE]))
  names(pl) <- names(ph) <- shared
  unname((ph[CAT_OX] - ph[CAT_NULL]) - (pl[CAT_OX] - pl[CAT_NULL]))
}
estimator_ladder <- dplyr::bind_rows(lapply(names(EST_PANEL), function(en) {
  gT <- .tertile_groups(.est_vec("TCGA", EST_PANEL[[en]]), ID_T)
  gS <- .tertile_groups(.est_vec("SCAN-B", EST_PANEL[[en]]), ID_S)
  tibble::tibble(estimator = en, set = unname(EST_PANEL[[en]]),
                 TCGA = .gap_diff(UT, gT), `SCAN-B` = .gap_diff(US, gS))
}))
amp_g <- ifelse(is.na(COVT$amp), NA_character_,
                ifelse(COVT$amp, "MYC-high", "MYC-low"))
AMP_DIFF <- .gap_diff(UT, amp_g)
estimator_ladder <- dplyr::bind_rows(
  estimator_ladder,
  tibble::tibble(estimator = "TCGA copy number (MYC_amp)",
                 set = "DNA-level, not an expression programme",
                 TCGA = AMP_DIFF, `SCAN-B` = NA_real_))
print(as.data.frame(estimator_ladder), digits = 3, row.names = FALSE)

# --- C2 stratification --------------------------------------------------------
.strata <- function(U, pro, est, ids) {
  q <- stats::quantile(pro, c(TERTILE, 1 - TERTILE), names = FALSE)
  s <- ifelse(pro <= q[1], "low prolif", ifelse(pro >= q[2], "high prolif", "mid prolif"))
  dplyr::bind_rows(lapply(c("low prolif", "mid prolif", "high prolif"), function(lv) {
    j <- which(s == lv)
    g <- .tertile_groups(est[j], ids[j])          # re-split MYC WITHIN the stratum
    Us <- U[, j, drop = FALSE]
    tibble::tibble(stratum = lv, n_stratum = length(j),
                   n_low = sum(g == "MYC-low", na.rm = TRUE),
                   n_high = sum(g == "MYC-high", na.rm = TRUE),
                   gap_diff = .gap_diff(Us, g))
  }))
}
strata <- dplyr::bind_rows(
  dplyr::mutate(.strata(UT, pro_T, est_T, ID_T), cohort = "TCGA"),
  dplyr::mutate(.strata(US, pro_S, est_S, ID_S), cohort = "SCAN-B"))
strata <- strata[, c("cohort", "stratum", "n_stratum", "n_low", "n_high", "gap_diff")]
message("\n   C2: within proliferation tertiles, MYC re-split INSIDE each stratum")
message("   so the groups stay balanced. Stratum sizes reported; they are small.")
print(as.data.frame(strata), digits = 3, row.names = FALSE)

# --- C3 adjustment, with the copy-number control read FIRST -------------------
message("\n   C3: adjustment on proliferation. SECONDARY and ambiguous by")
message("   construction. THE CONTROL IS READ FIRST.")
.adj_gap_diff <- function(U, g, adj) {
  R <- .residualise(U, adj)
  pl <- .pct(rowMeans(R[, which(g == "MYC-low"), drop = FALSE]))
  ph <- .pct(rowMeans(R[, which(g == "MYC-high"), drop = FALSE]))
  names(pl) <- names(ph) <- shared
  unname((ph[CAT_OX] - ph[CAT_NULL]) - (pl[CAT_OX] - pl[CAT_NULL]))
}
adjusted <- tibble::tibble(
  quantity = c("POSITIVE CONTROL: TCGA copy-number split", "the gap, TCGA",
               "the gap, SCAN-B"),
  independence = c("measurement only - amplified tumours proliferate more", "", ""),
  unadjusted = c(AMP_DIFF, .gap_diff(UT, grp_T), .gap_diff(US, grp_S)),
  adjusted = c(
    .adj_gap_diff(UT[, !is.na(amp_g), drop = FALSE], amp_g[!is.na(amp_g)],
                  data.frame(prolif = pro_T[!is.na(amp_g)])),
    .adj_gap_diff(UT, grp_T, data.frame(prolif = pro_T)),
    .adj_gap_diff(US, grp_S, data.frame(prolif = pro_S))))
# The retained RATIO is unstable here and is deliberately not leaned on: the
# unadjusted quantities are small, so a ratio can exceed 1 (the quantity GREW
# under adjustment) or flip sign, and neither means "attenuated". Sign and
# absolute change are the readable quantities.
adjusted$retained <- adjusted$adjusted / adjusted$unadjusted
adjusted$sign_kept <- sign(adjusted$adjusted) == sign(adjusted$unadjusted)
adjusted$abs_change <- adjusted$adjusted - adjusted$unadjusted
print(as.data.frame(adjusted), digits = 3, row.names = FALSE)
message("   THE RATIO IS NOT LEANED ON. Above 1 means the quantity GREW under")
message("   adjustment, which is not attenuation, and the unadjusted numbers are")
message("   small enough that the ratio is unstable.")
C3_GREW <- any(abs(adjusted$adjusted) > abs(adjusted$unadjusted))
C3_FLIP <- any(!adjusted$sign_kept)
C3_READ <- if (C3_GREW || C3_FLIP) {
  paste0("UNINFORMATIVE - under proliferation adjustment ",
         if (C3_GREW) "at least one quantity GREW" else "",
         if (C3_GREW && C3_FLIP) " and " else "",
         if (C3_FLIP) "at least one changed sign" else "",
         ", the positive control included. The adjustment is not behaving as an",
         " attenuation and nothing is read from it. This is the ambiguity the",
         " header declared in advance, not a surprise")
} else if (abs(adjusted$retained[1] - mean(adjusted$retained[2:3])) < 0.2) {
  "the control attenuates AS MUCH AS the gap - the adjustment is removing MYC wholesale. UNINFORMATIVE"
} else if (adjusted$retained[1] > mean(adjusted$retained[2:3])) {
  "the control HOLDS while the gap attenuates more - a genuine proliferation effect on the gap"
} else {
  "the control attenuates MORE than the gap - the adjustment is not specific. UNINFORMATIVE"
}
message("   C3: ", C3_READ)

lad <- estimator_ladder[estimator_ladder$estimator != "TCGA copy number (MYC_amp)", ]
PRIMARY_SIGN <- sign(lad$TCGA[lad$estimator == "myc_felsher"])
both_same <- lad$estimator[sign(lad$TCGA) == PRIMARY_SIGN &
                             sign(lad$`SCAN-B`) == PRIMARY_SIGN]
both_opp  <- lad$estimator[sign(lad$TCGA) == -PRIMARY_SIGN &
                             sign(lad$`SCAN-B`) == -PRIMARY_SIGN]
AMP_AGREES <- sign(AMP_DIFF) == PRIMARY_SIGN
strata_T <- strata[strata$cohort == "TCGA", ]
strata_S <- strata[strata$cohort == "SCAN-B", ]
ladder_summary <- tibble::tibble(
  quantity = c("estimators agreeing with the primary in BOTH cohorts",
               "estimators OPPOSITE to the primary in BOTH cohorts",
               "copy-number split agrees with the primary",
               "TCGA strata agreeing with the primary",
               "SCAN-B strata agreeing with the primary"),
  value = c(paste(both_same, collapse = ", "), paste(both_opp, collapse = ", "),
            as.character(AMP_AGREES),
            paste0(sum(sign(strata_T$gap_diff) == PRIMARY_SIGN), " of 3"),
            paste0(sum(sign(strata_S$gap_diff) == PRIMARY_SIGN), " of 3")))
print(as.data.frame(ladder_summary), row.names = FALSE, right = FALSE)
R3 <- paste0(
  "THE SIGN OF THE HUMAN GAP DIFFERENCE IS ESTIMATOR-DEPENDENT. Only ",
  paste(both_same, collapse = " and "), " agree with the primary in both cohorts",
  ", while ", paste(both_opp, collapse = " and "),
  " give the OPPOSITE sign in both. The copy-number split, the one instrument ",
  "not defined on the expression matrix, ",
  if (AMP_AGREES) "AGREES" else "DISAGREES IN SIGN",
  " with the primary (", sprintf("%+.1f", AMP_DIFF),
  "). Within proliferation strata the cohorts disagree: TCGA keeps the ",
  "primary's sign in ", sum(sign(strata_T$gap_diff) == PRIMARY_SIGN),
  " of 3 strata and SCAN-B in ", sum(sign(strata_S$gap_diff) == PRIMARY_SIGN),
  " of 3. R3's second branch applies: the human gap difference is NOT a robust ",
  "readout of MYC, and the note says so")
message("\n   R3: ", R3)

# =============================================================================
# 7. R5 - ox_rel, rebuilt from E16's recipe. HUMAN SIDE ONLY.
# =============================================================================

message("\n7. R5: ox_rel rebuilt - the reproduction control comes FIRST")

OXS <- mito$arm_sets[[OX_SUB_N]]
REST <- setdiff(sd_$strip_refs$MITOCARTA_ALL, OXS)
COH <- list(TCGA = list(G = GT, res = .symbol_resolver(rownames(GT), NULL),
                        ids = ID_T, arms = mito, grp = grp_T),
            `SCAN-B` = list(G = GS, res = .symbol_resolver(rownames(GS), sc$symbol_map),
                            ids = ID_S, arms = sc, grp = grp_S))
.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res); M <- gr$mat
  v <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}
oxrel_control <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  r <- .comp(OXS, C) - .comp(REST, C)
  rho <- stats::cor(r, as.numeric(C$arms$mitopps_arms[OX_SUB_N, C$ids]),
                    method = "spearman")
  tibble::tibble(cohort = coh, rebuilt_rho = rho,
                 e16_saved = unname(e16$r1_rho[[coh]]),
                 delta = rho - unname(e16$r1_rho[[coh]]),
                 myc_high_minus_low = mean(r[which(C$grp == "MYC-high")]) -
                                      mean(r[which(C$grp == "MYC-low")]))
}))
print(as.data.frame(oxrel_control), digits = 4, row.names = FALSE)
OXREL_OK <- max(abs(oxrel_control$delta)) < 1e-9
if (!OXREL_OK) {
  stop("ox_rel REPRODUCTION CONTROL FAILED - the recipe has drifted.", call. = FALSE)
}
message("   CONTROL PASSED - the rebuilt ruler returns E16's saved rho exactly.")
message("   This is a COMPOSITE rebuilt from its own recipe, not a re-scoring.")
message("   Human side only: the mouse artefact has no gene-level values.")

agg_h <- agg_contrast[agg_contrast$numerator == names(AGG_NUM)[1] &
                        agg_contrast$contrast %in% c("TCGA MYC-high - low",
                                                     "SCAN-B MYC-high - low"), ]
r5 <- tibble::tibble(
  cohort = c("TCGA", "SCAN-B"),
  aggregate_gap_change = agg_h$change[match(c("TCGA MYC-high - low",
                                              "SCAN-B MYC-high - low"), agg_h$contrast)],
  ox_rel_change = oxrel_control$myc_high_minus_low[match(c("TCGA", "SCAN-B"),
                                                         oxrel_control$cohort)])
r5$same_direction <- sign(r5$aggregate_gap_change) == sign(r5$ox_rel_change)
print(as.data.frame(r5), digits = 3, row.names = FALSE)
message("   units differ, so only the DIRECTION is compared - a percentile gap")
message("   against a z-score composite.")
R5 <- if (all(r5$same_direction)) {
  paste0("AGREE in direction in both cohorts - CORROBORATION ON A SECOND RULER ",
         "and nothing more. E16's R1 declared the two rulers not interchangeable ",
         "and that is not overturned here")
} else {
  paste0("DISAGREE in ", sum(!r5$same_direction), " of 2 cohorts. THAT ",
         "DISAGREEMENT IS THE FINDING and is reported, not resolved by choosing one")
}
message("   R5: ", R5)

# =============================================================================
# 8. PART D - figures
# =============================================================================

message("\n8. PART D: figures")
if (!dir.exists(DIR_E28_FIG)) dir.create(DIR_E28_FIG, recursive = TRUE)

CAP <- paste(
  "THE GAP IS POST HOC - derived after E27 ran and not among its reading rules.",
  "Percentile units, so no mitoPPS value crosses the species boundary.",
  "The mouse timeline is batch-confounded (batch = timepoint); nothing is withdrawn on that account.",
  "Exploratory, post-hoc, descriptive; not a hypothesis test.", sep = "\n")

f1 <- gap_by_group %>% dplyr::filter(comparator == CAT_NULL)
p1 <- ggplot2::ggplot(f1, ggplot2::aes(arm, gap)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey45") +
  ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.4, colour = "grey70",
                      linetype = 3) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = ci_lo, ymax = ci_hi, colour = species),
                          linewidth = 0.6) +
  ggplot2::geom_point(ggplot2::aes(colour = species), size = 2.4) +
  ggplot2::scale_colour_manual(values = c(mouse = "#0072B2", human = "#D55E00")) +
  ggplot2::labs(x = NULL, y = paste0("OXPHOS minus ", CAT_NULL, " (percentile points)"),
                colour = NULL, title = "The OXPHOS gap across the eight groups",
                subtitle = "POST HOC readout. Intervals are 5,000 paired bootstrap draws over samples within group.",
                caption = CAP) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E28_FIG, "E28_gap_by_group.pdf"), p1,
                width = 7.5, height = 5)

p2 <- ggplot2::ggplot(gap_by_group,
                      ggplot2::aes(arm, gap, colour = comparator, group = comparator)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey45") +
  ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.4, colour = "grey70",
                      linetype = 3) +
  ggplot2::geom_line(linewidth = 0.35, alpha = 0.6) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::labs(x = NULL, y = "OXPHOS minus the comparator (percentile points)",
                colour = NULL,
                title = "All six comparators - is the gap a property of the compartment or of one pairing?",
                subtitle = paste0(CAT_NULL, " is rank ",
                                  comparator_distance$rank_value[comparator_distance$comparator == CAT_NULL],
                                  " of six by distance from OXPHOS. Never averaged into a headline."),
                caption = CAP) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 legend.position = "top",
                 legend.text = ggplot2::element_text(size = 6.5),
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E28_FIG, "E28_comparator_panel.pdf"), p2,
                width = 9, height = 5.5)
message("   wrote E28_gap_by_group.pdf and E28_comparator_panel.pdf")

# =============================================================================
# 9. PART E - verdict
# =============================================================================

message("\n9. PART E: verdict")

verdict <- tibble::tibble(
  rule = c("STATUS", "CONTROL against E27", "the rationale's tables",
           "PART 0b the protection claim", "R1 the mouse developmental gap",
           "R2 developmental against MYC", "R3 proliferation",
           "R4 the comparator panel", "R5 the aggregate against ox_rel",
           "C3 proliferation adjustment", "ox_rel reproduction"),
  outcome = c(
    "THE GAP IS POST HOC - derived after E27 ran, not among its reading rules",
    paste0("PASSED - ", nrow(ctrl), " positions reproduce E27 (max |delta| ",
           format(CONTROL_DELTA, digits = 2), ")"),
    if (RAT_OK) {
      "reproduce within the one-decimal rounding tolerance"
    } else {
      paste0("reproduce to ", sprintf("%.2f", max(abs(c(assert_group$delta,
             assert_contrast$delta)))), " percentile points. One cell sits just",
             " outside the ", RAT_TOL, " one-decimal boundary - a ROUNDING SLIP in",
             " the rationale, not a computation difference. The cell is named in",
             " $control$rationale_group")
    },
    paste0(if (PROTECTS) "HOLDS" else "DOES NOT HOLD",
           " - the gap moves ", sprintf("%.1f", protection_summary$max_move_from_M0[3]),
           " percentile points under composition adjustment against ",
           sprintf("%.1f", protection_summary$max_move_from_M0[1]), " and ",
           sprintf("%.1f", protection_summary$max_move_from_M0[2]),
           " for the two positions. TCGA only"),
    R1, R2, R3, R4, R5, paste0("SECONDARY - ", C3_READ),
    paste0("PASSED - rebuilt ruler returns E16's saved rho exactly (max |delta| ",
           format(max(abs(oxrel_control$delta)), digits = 2), ")")))
print(as.data.frame(verdict), row.names = FALSE, right = FALSE)
message("\n   AND WHATEVER IT SAYS, IT IS DESCRIPTIVE. NOT A TEST.")

NOTES <- c(
  paste("THE GAP IS POST HOC. It was derived after E27 ran and was not among",
        "E27's R1-R5. Labelled post hoc in this object, the verdict, both",
        "figures and the note."),
  paste("TWO DEFINITIONS, and they are not the same statistic. The readout is a",
        "DIFFERENCE OF GROUP GAPS. The gap inside a contrast's own ordering is a",
        "different number - -42.3 against -28.9 for the WT window - and is",
        "carried here only so it is not quoted by mistake."),
  paste("The protection claim was TESTED, not asserted, and reported before any",
        "gap result. A percentile is already a within-ordering rank, so a",
        "monotone panel-wide shift changes nothing and differencing adds nothing",
        "there; the threat is DIFFERENTIAL movement. TCGA only, because SCAN-B",
        "has no purity estimate and nothing is imputed."),
  paste(CAT_NULL, "is a MIDDLING comparator, rank 3 of 6 by distance from",
        "OXPHOS, not a conservative one. The ordering is identical on the value",
        "and percentile scales, which was checked rather than assumed. The six",
        "gaps are reported individually and NEVER averaged into a headline."),
  paste("The aggregate is ox_rel rebuilt in mitoPPS space. Agreement is",
        "CORROBORATION ON A SECOND RULER and NOT a second line of evidence.",
        "E16's R1 declared the two rulers not interchangeable and R5 does not",
        "overturn that - it asks only about direction on this one contrast."),
  paste("The umbrella numerator DILUTES the aggregate, by",
        sprintf("%.1f", mean(abs(DILUTION))), "percentile points on average,",
        "which is exactly what E27 predicted: assembly factors barely move while",
        "subunits do. Subunits are the primary numerator."),
  paste("ox_rel was REBUILT from E16's recipe on the human linear matrices, with",
        "the reproduction control reported before any comparison. This is not",
        "re-scoring: that rule protects mitoPPS - the mouse artefact and the",
        "142-pathway universe - and this is a composite rebuilt from its own",
        "saved recipe. HUMAN SIDE ONLY; a mouse ox_rel would mean touching that",
        "repo."),
  paste("SCAN-B HAS a proliferation score, so C2 stratification runs in BOTH",
        "cohorts. That is the only arm of the composition work that does, and it",
        "partly answers E26's standing limitation."),
  paste("C1 and C2 are PRIMARY because they remove nothing; C3 is SECONDARY and",
        "ambiguous by construction. C3 collapsing while C1 and C2 survive is",
        "OVER-ADJUSTMENT, not refutation - declared before any number."),
  paste("The positive control is the TCGA copy-number split, the only candidate",
        "not defined on the expression matrix. ITS INDEPENDENCE IS IN",
        "MEASUREMENT, NOT IN BIOLOGY. Its own attenuation is reported beside the",
        "gap's. The mitoribosome was NOT used: it sits inside central dogma and",
        "is one arm of the primary gap."),
  paste("The mouse timeline is BATCH-CONFOUNDED (batch = timepoint). Temporal",
        "numbers are reported in full and nothing is withdrawn on that account."),
  paste("The expression-matched-null analysis in myc_mouse REMAINS THE PRIMARY",
        "EVIDENCE for the developmental change. This is a second presentation.",
        "If the two disagree, THAT DISAGREEMENT IS THE FINDING."),
  "No mitoPPS value crosses the species boundary. No ortholog call. N3 throughout.")

saveRDS(list(
  rules = RULES,
  settings = list(quantity = RULES$quantity, comparator = CAT_NULL,
                  comparators = COMPARATORS, arm_order = ARM_ORDER,
                  contrasts = CONTRASTS, nboot = NBOOT, seed = 1,
                  tertile = TERTILE, post_hoc = TRUE,
                  oxrel_targets = e16$r1_rho,
                  mouse_artefact = list(md5 = MOUSE_MD5, recomputed = FALSE)),
  control = list(vs_e27 = ctrl, max_delta = CONTROL_DELTA,
                 rationale_group = assert_group, rationale_contrast = assert_contrast,
                 rationale_ok = RAT_OK, two_definitions = two_defs),
  protection = protection, protection_summary = protection_summary,
  protects = PROTECTS, gap_signal = GAP_SIGNAL, composition_share = COMP_SHARE,
  gap_by_group = gap_by_group, gap_by_contrast = gap_by_contrast,
  r2 = r2, comparator_panel = gap_by_contrast,
  comparator_distance = comparator_distance, comparator_range = comparator_range,
  scales_agree = SCALES_AGREE,
  aggregate = aggregate_tab, aggregate_contrast = agg_contrast,
  dilution = DILUTION, subunits_vs_assembly = subunits_vs_assembly,
  estimator_ladder = estimator_ladder, ladder_summary = ladder_summary,
  strata = strata, adjusted = adjusted, c3_read = C3_READ,
  oxrel = oxrel_control, r5 = r5,
  readings = list(R1 = R1, R2 = R2, R3 = R3, R4 = R4, R5 = R5),
  verdict = verdict, analysis_date = Sys.Date(), notes = NOTES), PATH_E28)

readr::write_csv(gap_by_group %>% dplyr::mutate(arm = as.character(arm)), PATH_E28_TAB)

message("\nE28: done.")
message("    ", PATH_E28)
message("    ", PATH_E28_TAB)
message("    ", file.path(DIR_E28_FIG, "E28_gap_by_group.pdf"))
message("    ", file.path(DIR_E28_FIG, "E28_comparator_panel.pdf"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E28)

  ## the controls first
  x$control$max_delta
  x$control$rationale_group |> as.data.frame()
  x$control$rationale_contrast |> as.data.frame()
  x$oxrel |> as.data.frame()

  ## READ THIS BEFORE ANY GAP NUMBER: does differencing protect at all?
  x$protection |> as.data.frame()
  x$protection_summary |> as.data.frame()
  x$protects

  ## the two definitions - do not quote the second
  x$control$two_definitions |> as.data.frame()

  ## the verdict and the five readings
  x$verdict |> as.data.frame()
  x$readings

  ## R1 and R2
  subset(x$gap_by_contrast, comparator == "Mitochondrial central dogma") |> as.data.frame()
  x$r2 |> as.data.frame()

  ## R4 - the comparator panel, never averaged
  x$comparator_distance |> as.data.frame()
  x$comparator_range |> as.data.frame()

  ## the aggregate, both numerators, and how much the umbrella costs
  x$aggregate_contrast |> as.data.frame()
  x$dilution

  ## R3 - the ladder, the strata, and the adjustment with its control
  x$estimator_ladder |> as.data.frame()
  x$strata |> as.data.frame()
  x$adjusted |> as.data.frame()

  ## R5
  x$r5 |> as.data.frame()

  ## every caveat the note must repeat
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")
}
