# E48_separation_and_fig4g.R
# =============================================================================
# OXPHOS, PROLIFERATION AND MYC ACTIVITY ON BBC3, BCL2L1 AND MCL1, EACH WITH
# THE OTHERS HELD FIXED; AND FIG. 4G IN THREE LAYOUTS.
#
# EXPLORATORY and POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-11-02_E48_separation_declaration.md, committed ALONE
# at 2941438, before this script existed. The governing brief, the author's
# decision record, is docs/2026-11-02_separation_line_brief.md (cbf2e4c).
# Every choice below is the declaration's; section numbers cited are its
# sections. Read sections 2.9, 7, 8 and 11 first.
#
# =============================================================================
# WHAT THIS DOES
# =============================================================================
# PART A FITS NOTHING (sections 4-7). From E47's verified saved coefficients
# and variance matrices it recovers, for the gene fits only:
#   - proliferation's m3 coefficient: a post-hoc characterisation of a
#     declared covariate (brief 4.7), on the E43 precedent;
#   - the within-model difference, proliferation minus OXPHOS, at m3, as the
#     general linear contrast L beta with Var = L V L' (E44's orientation).
# It labels BBC3, BCL2L1 and MCL1 over their four M_a cells (7.1), applies the
# agreement check on the least-entangled score (7.2), and tallies the nine
# other genes descriptively (7.4). Section 8's closing paragraph: nothing
# computed is discarded; everything is saved, tabulated and displayed.
#
# PART B (MYC activity on pCR) IS UNAVAILABLE (section 9). It is recorded
# with its structure facts, and nothing is fitted.
#
# PART C IS REPORTED AND NEVER READ (section 10). It changes no label.
#   C1  TCGA, M_a, m3, with TP53 status added: 6 fits.
#   C3  "more steeply" at m3, by E47's bootstrap and E47's seeds: 4 cells.
#   C4  TCGA with SCAN-B's spine (PAM50 only): 12 fits.
#   C2 is NOT RUN (10.2): PROLIF_STD shares nine genes with M_a.
#
# DISPLAYS (section 12): Fig. 4g in three layouts (the author chooses), ED-a
# to ED-d, and S-1. EVERY E48 PANEL FOLLOWS E45's FILL CONVENTION: OPEN =
# without the proliferation term, FILLED = with it. E47's Displays A and B use
# the opposite convention; they are not redrawn here.
#
# =============================================================================
# TERMS
# =============================================================================
#   level, priority   OXPHOS subunits scored on GSVA (VST) and on mitoPPS
#                     (linear), named in plain words on every panel
#   estimator         a MYC activity score. Never a "signature"
#   M_a               FELSHER__MITOSTRIP, the reference estimator (primary)
#   clean             MYC_UP.V1_UP__MITOSTRIP, the least-entangled score
#   PROLIF            PROLIF_DISJOINT (318 genes), E47's term "PROLIF"
#   m1, m3            E47's rungs: OX + MYC + spine; m1 + PROLIF
#   difference        beta(PROLIF) - beta(OX) within one m3 fit
# No panel carries the notation beta, m1, m3, M_a or Delta (12.1, asserted).
#
# =============================================================================
# SNAPSHOT AND SCALE
# =============================================================================
# Reads E47's saved object: its coefficients, variance matrices and analysis
# frame, in which every score was Blom-scored ONCE PER COHORT over the whole
# cohort, as E34 scored. Also reads E45's and E39's saved objects, E40's and
# E43's objects for structure only (Part B), and TP53_status from the TCGA
# covariate snapshot data/from_validation/tcga_brca_covariates.rds (C1).
# IT RESCORES NOTHING. The joint snapshot is not read.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It does not refit E47 or amend its MIXED verdict. It reads no S6, S5,
#     S7 or arm fit and draws none (asserted).
#   - It fits no product term and tests no interaction.
#   - It reads no Part C result. C1, C3 and C4 change no label.
#   - It fits nothing for Part B.
#   - It chooses no layout for Fig. 4g. It renders all three.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

message("\nE48: OXPHOS, proliferation and MYC activity on BBC3, BCL2L1 and MCL1\n",
        strrep("=", 78))
message("EXPLORATORY, POST-HOC. Declaration 2941438; brief cbf2e4c.")
message("Part A fits nothing. Part C is reported and never read.")

PATH_E47 <- file.path(DIR_RESULTS, "e47_mutual_adjustment.rds")
PATH_E45 <- file.path(DIR_RESULTS, "e45_fig4_human_panels.rds")
PATH_E39 <- file.path(DIR_RESULTS, "e39_prepare_and_gate.rds")
PATH_E40 <- file.path(DIR_RESULTS, "e40_ladder.rds")
PATH_E43 <- file.path(DIR_RESULTS, "e43_proliferation_coefficient.rds")
PATH_E48  <- file.path(DIR_RESULTS, "e48_separation_and_fig4g.rds")

DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-03"
# Files are named by CONTENT, never by figure letter (CLAUDE.md).
FIG <- c(layout1 = "E48_separation_layout1_two_scores",
         layout2 = "E48_separation_layout2_three_scores",
         layout3 = "E48_separation_layout3_with_outcomes",
         ed_a    = "E48_score_correlations",
         ed_b    = "E48_bcl2l1_against_bbc3",
         ed_c    = "E48_myc_score_ladder",
         ed_d    = "E48_twelve_genes_three_coefficients",
         s1      = "E48_sensitivities_tp53_composition")

# =============================================================================
# 0. CONSTANTS - all from the declaration, none chosen here
# =============================================================================
COHORTS <- c("TCGA", "SCAN-B")
INSTR   <- c(gsva = "level", mitopps = "priority")
BLOCK_ORDER <- c("TCGA gsva", "TCGA mitopps", "SCAN-B gsva", "SCAN-B mitopps")

# MYC activity scores, in E47's order, which is the order of proliferation
# entanglement (section 5 of E47's declaration).
EST_ALL   <- c("clean", "M_a", "entangled", "regulon")
EST_READ  <- "M_a"
EST_AGREE <- "clean"
EST_NAME  <- c(clean = "MYC_UP.V1_UP (1.5%)",
               M_a = "Felsher, reference (14.8%)",
               entangled = "Hallmark MYC targets (23.5%)",
               regulon = "CollecTRI MYC regulon")

READ3  <- c("BBC3", "BCL2L1", "MCL1")                       # section 3
PRO    <- c("BAD", "BIK", "BBC3", "BID", "BCL2L11", "BMF", "PMAIP1")
ANTI   <- c("BCL2L1", "MCL1", "BCL2", "BCL2L2", "BCL2A1")
TWELVE <- c(PRO, ANTI)
NINE   <- setdiff(TWELVE, READ3)                             # 7.4
NOT_USED <- c("S6", "S5", "S7", "arm_trigger", "arm_guardian")

# E47's spines (section 4), and C4's: SCAN-B's spine applied to TCGA (10.4).
RUNG_COV <- list(
  TCGA     = list(m1 = c("PAM50", "purity", "leuko"),
                  m3 = c("PAM50", "purity", "leuko", "PROLIF")),
  `SCAN-B` = list(m1 = "PAM50", m3 = c("PAM50", "PROLIF")))
C4_COV <- c("PAM50", "PROLIF")
PAM50_LEVELS <- c("LumA", "LumB", "HER2", "Basal", "Normal")   # E34's order

EXPECT_N    <- c(TCGA = 938L, `SCAN-B` = 3143L)
TP53_EXPECT <- c(mutant = 325L, wildtype = 613L)               # 10.1
E39_DECLARED <- c(GSE25066 = 0.7797, GSE194040 = 0.7851, GSE164458 = 0.8236)

BAND_SEP <- 0.60    # E39's bands
BAND_ID  <- 0.80
TOL      <- 1e-10
CI_Q     <- 0.975
BOOT_N   <- 2000L   # E47 section 9, reused by C3

# Colours (12.1). OXPHOS and proliferation are E45's. MYC activity is the
# first of #0072B2, #CC79A7, #E69F00 that passes the dataviz palette checker
# against the other two, all pairs, light surface: the rule the declaration
# fixed. Run 2026-11-03, before this script was written:
#   "#009E73,#D55E00,#0072B2" --mode light --pairs all
#     lightness band PASS; chroma floor PASS; CVD separation worst all-pairs
#     dE 11.0 (deutan) PASS; normal-vision worst dE 18.7 PASS; contrast vs
#     surface all >= 3:1 PASS. ALL CHECKS PASS, no warning. -> #0072B2.
#   The same checker reproduces E45's recorded pair result (dE 11.0, 25.8).
COL_OX  <- "#009E73"
COL_PR  <- "#D55E00"
COL_MYC <- "#0072B2"
VALIDATED_PALETTE <- c(COL_OX, COL_PR, COL_MYC)

# Sizes: 178 mm wide (E45's width), no text below 5.5 pt (12.1).
BASE_PT <- 6.5
MIN_PT  <- 5.5
TXT_MM  <- BASE_PT / ggplot2::.pt
SUB_MM  <- MIN_PT / ggplot2::.pt
W_PAGE  <- 178 / 25.4

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.fmt2    <- function(e, lo, hi) sprintf("%.2f (%.2f, %.2f)", e, lo, hi)
.fmt3    <- function(e, lo, hi) sprintf("%+.3f [%+.3f, %+.3f]", e, lo, hi)
.dev     <- function(a, b) if (length(a)) max(abs(a - b)) else NA_real_
.band    <- function(r) ifelse(abs(r) >= BAND_ID, "0.80 or above",
                         ifelse(abs(r) >= BAND_SEP, "0.60 to 0.80",
                                "below 0.60"))
.fid     <- function(coh, inst, est, rung, rd)
  paste(coh, inst, est, rung, rd, sep = " | ")

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

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

for (p in c(PATH_E47, PATH_E45, PATH_E39, PATH_E40, PATH_E43, PATH_TCGA_COV)) {
  .stop_if(file.exists(p), "absent: ", p)
}
e47 <- readRDS(PATH_E47)
e45 <- readRDS(PATH_E45)
e39 <- readRDS(PATH_E39)
e40 <- readRDS(PATH_E40)
e43 <- readRDS(PATH_E43)
CV  <- readRDS(PATH_TCGA_COV)$covariates

NEED47 <- c("coefficient_matrix", "vcov", "results", "vif", "fit_digests",
            "frame", "gates", "displays", "bcl2l1_bbc3", "boot_replicates")
.assert("E47's saved object holds every element E48 reads",
        all(NEED47 %in% names(e47)), n = length(NEED47))
CM  <- e47$coefficient_matrix
VC  <- e47$vcov
R47 <- e47$results
G3  <- e47$gates$g3

# =============================================================================
# 2. Structure (section 6.1), asserted again before anything is read
# =============================================================================
message("\n2. structure")

ID13 <- unique(CM$fit_id[CM$rung %in% c("m1", "m3")])
st_ok <- vapply(ID13, function(f) {
  tm <- CM$term[CM$fit_id == f]
  v  <- VC[[f]]
  !is.null(v) && nrow(v) == ncol(v) && identical(rownames(v), tm) &&
    identical(colnames(v), tm)
}, logical(1))
.assert(paste0("every m1 and m3 fit has a square variance matrix named by its ",
               "coefficient terms, in order"), all(st_ok), n = length(ID13))
ID3 <- unique(CM$fit_id[CM$rung == "m3"])
.assert("PROLIF is in every m3 fit, in its coefficients and its variance matrix",
        all(vapply(ID3, function(f) "PROLIF" %in% CM$term[CM$fit_id == f] &&
                     "PROLIF" %in% rownames(VC[[f]]), logical(1))),
        n = length(ID3))
.assert("PROLIF is in none of the m1 fits",
        !any(CM$term[CM$rung == "m1"] == "PROLIF"),
        n = length(unique(CM$fit_id[CM$rung == "m1"])))

# =============================================================================
# 3. Part A: extraction from E47's saved fits (section 5), gates G-A1 and G-A2
# =============================================================================
message("\n3. Part A: extraction")

# Section 5 and section 8's closing paragraph: the twelve genes on M_a and on
# the least-entangled score; the three read genes on the regulon and on
# HALLMARK, whose m3 cells are void.
GRID_A <- dplyr::bind_rows(
  tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                     estimator = c(EST_READ, EST_AGREE), readout = TWELVE),
  tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                     estimator = c("regulon", "entangled"), readout = READ3))
.assert("Part A reads gene fits only: no signature, variant or arm",
        !any(GRID_A$readout %in% NOT_USED) && all(GRID_A$readout %in% TWELVE),
        n = nrow(GRID_A))

.get_fit <- function(fid) {
  r <- CM[CM$fit_id == fid, ]
  .stop_if(nrow(r) > 0L, "fit absent from E47: ", fid)
  list(fid = fid,
       b  = stats::setNames(unname(r$estimate), r$term),
       se = stats::setNames(unname(r$se), r$term),
       lo = stats::setNames(unname(r$ci_lo), r$term),
       hi = stats::setNames(unname(r$ci_hi), r$term),
       df = unique(r$df), n = unique(r$n), V = VC[[fid]])
}

# Quantity 2: the general linear contrast. L carries +1 on PROLIF and -1 on
# OX, so Var = L V L' takes the covariance from the algebra (CLAUDE.md: never
# difference two coefficients without it).
.contrast <- function(b, V, df) {
  L <- stats::setNames(numeric(length(b)), names(b))
  L[["PROLIF"]] <- 1
  L[["OX"]]     <- -1
  d  <- sum(L * b)
  se <- sqrt(as.numeric(t(L) %*% V[names(b), names(b)] %*% L))
  tq <- stats::qt(CI_Q, df)
  c(diff = d, diff_se = se, diff_lo = d - tq * se, diff_hi = d + tq * se,
    coef_cor = V["PROLIF", "OX"] / sqrt(V["PROLIF", "PROLIF"] * V["OX", "OX"]),
    naive_se = sqrt(V["PROLIF", "PROLIF"] + V["OX", "OX"]))
}

DIGS <- e47$fit_digests
.ga2 <- function(f) {
  dg <- DIGS[DIGS$fit_id == f$fid, ]
  tq <- stats::qt(CI_Q, f$df)
  tibble::tibble(
    fit_id = f$fid,
    coef_digest = identical(digest::digest(f$b), dg$digest_coef),
    vcov_digest = identical(digest::digest(f$V), dg$digest_vcov),
    se_dev = .dev(f$se, sqrt(diag(f$V))[names(f$se)]),
    ci_dev = max(.dev(f$lo, f$b - tq * f$se), .dev(f$hi, f$b + tq * f$se)))
}

EXT <- vector("list", nrow(GRID_A))
GA2 <- vector("list", nrow(GRID_A))
for (i in seq_len(nrow(GRID_A))) {
  g  <- GRID_A[i, ]
  f1 <- .get_fit(.fid(g$cohort, g$instrument, g$estimator, "m1", g$readout))
  f3 <- .get_fit(.fid(g$cohort, g$instrument, g$estimator, "m3", g$readout))
  ct <- .contrast(f3$b, f3$V, f3$df)
  # The same difference recomputed a second way, from the three variance
  # elements directly (assertion 2 of section 14.2).
  alt_d  <- f3$b[["PROLIF"]] - f3$b[["OX"]]
  alt_se <- sqrt(f3$V["PROLIF", "PROLIF"] + f3$V["OX", "OX"] -
                   2 * f3$V["PROLIF", "OX"])
  GA2[[i]] <- dplyr::bind_rows(.ga2(f1), .ga2(f3))
  EXT[[i]] <- tibble::tibble(
    cohort = g$cohort, instrument = g$instrument, estimator = g$estimator,
    readout = g$readout, n = f3$n, df_m3 = f3$df,
    void = isTRUE(R47$void[R47$fit_id == f3$fid]),
    b_OX_m1 = f1$b[["OX"]], ox_lo_m1 = f1$lo[["OX"]], ox_hi_m1 = f1$hi[["OX"]],
    b_MYC_m1 = f1$b[["MYC"]], myc_lo_m1 = f1$lo[["MYC"]],
    myc_hi_m1 = f1$hi[["MYC"]],
    b_OX_m3 = f3$b[["OX"]], ox_lo_m3 = f3$lo[["OX"]], ox_hi_m3 = f3$hi[["OX"]],
    b_MYC_m3 = f3$b[["MYC"]], myc_lo_m3 = f3$lo[["MYC"]],
    myc_hi_m3 = f3$hi[["MYC"]],
    b_PR = f3$b[["PROLIF"]], pr_se = f3$se[["PROLIF"]],
    pr_lo = f3$lo[["PROLIF"]], pr_hi = f3$hi[["PROLIF"]],
    diff = ct[["diff"]], diff_se = ct[["diff_se"]], diff_lo = ct[["diff_lo"]],
    diff_hi = ct[["diff_hi"]], coef_cor = ct[["coef_cor"]],
    naive_se = ct[["naive_se"]],
    alt_dev = max(abs(ct[["diff"]] - alt_d), abs(ct[["diff_se"]] - alt_se)))
}
EXT <- dplyr::bind_rows(EXT)
GA2 <- dplyr::bind_rows(GA2) %>% dplyr::distinct(fit_id, .keep_all = TRUE)

.assert(paste0("G-A2: every coefficient vector read equals its fit's coef(), by ",
               "E47's digest"), all(GA2$coef_digest), n = nrow(GA2))
.assert(paste0("G-A2: every variance matrix read equals its fit's vcov(), by ",
               "E47's digest"), all(GA2$vcov_digest), n = nrow(GA2))
.assert("G-A2: every SE equals sqrt(diag(vcov))", max(GA2$se_dev) < TOL,
        n = nrow(GA2), max_dev = max(GA2$se_dev))
.assert("G-A2: every CI equals one recomputed with qt(0.975, df)",
        max(GA2$ci_dev) < TOL, n = nrow(GA2), max_dev = max(GA2$ci_dev))
.assert(paste0("each difference recomputed independently from coef and vcov ",
               "equals the stored one (section 14.2, 2)"),
        max(EXT$alt_dev) < 1e-12, n = nrow(EXT), max_dev = max(EXT$alt_dev))

# G-A1: the re-read OXPHOS and MYC coefficients equal E47's $results exactly.
R47W <- R47 %>%
  dplyr::filter(rung %in% c("m1", "m3")) %>%
  dplyr::transmute(cohort, instrument, estimator, readout, rung,
                   b_OX = unname(b_OX), ox_lo = unname(ox_lo),
                   ox_hi = unname(ox_hi), b_MYC = unname(b_MYC),
                   myc_lo = unname(myc_lo), myc_hi = unname(myc_hi)) %>%
  tidyr::pivot_wider(names_from = rung,
                     values_from = c(b_OX, ox_lo, ox_hi, b_MYC, myc_lo, myc_hi))
CMP <- EXT %>% dplyr::left_join(R47W, suffix = c("", "_r47"),
                                by = c("cohort", "instrument", "estimator",
                                       "readout"))
GA1_COLS <- c("b_OX_m1", "ox_lo_m1", "ox_hi_m1", "b_MYC_m1", "myc_lo_m1",
              "myc_hi_m1", "b_OX_m3", "ox_lo_m3", "ox_hi_m3", "b_MYC_m3",
              "myc_lo_m3", "myc_hi_m3")
GA1_DEV <- max(vapply(GA1_COLS, function(cc)
  .dev(CMP[[cc]], CMP[[paste0(cc, "_r47")]]), numeric(1)))
.assert("G-A1: the re-read OXPHOS and MYC coefficients equal E47's $results",
        !anyNA(CMP[paste0(GA1_COLS, "_r47")]) && GA1_DEV < TOL,
        n = nrow(CMP) * length(GA1_COLS), max_dev = GA1_DEV)
.assert(paste0("G-A3: the only void cells read are HALLMARK's m3 cells, from ",
               "E47's $results"),
        all(EXT$void == (EXT$estimator == "entangled")), n = nrow(EXT))

# G-A4, nesting (section 6.2): MYC's and OXPHOS's m1-to-m3 shifts equal minus
# proliferation's m3 coefficient times that term's coefficient in the
# regression of PROLIF on m1's terms, in the same tumours.
FRM <- e47$frame %>%
  dplyr::mutate(PAM50 = factor(PAM50, levels = PAM50_LEVELS))
DSET <- lapply(stats::setNames(COHORTS, COHORTS), function(coh) {
  d <- FRM[FRM$cohort == coh & FRM$in_set, ]
  d$PAM50 <- droplevels(d$PAM50)
  d
})
for (coh in COHORTS) {
  .assert(paste0(coh, ": the analysis set is E34's, n ", EXPECT_N[[coh]]),
          nrow(DSET[[coh]]) == EXPECT_N[[coh]] &&
            !anyNA(DSET[[coh]][, c("PAM50", "PROLIF", READ3)]),
          n = nrow(DSET[[coh]]))
}
.design <- function(d, inst, est, cov) {
  dd <- data.frame(OX = d[[paste0("OX_", inst)]], MYC = d[[paste0("MYC_", est)]])
  for (cv in cov) dd[[cv]] <- d[[cv]]
  for (fc in intersect(cov, c("PAM50", "TP53_status"))) {
    dd[[fc]] <- droplevels(dd[[fc]])
  }
  .stop_if(!anyNA(dd), "NA inside a design: ", inst, " ", est)
  stats::model.matrix(stats::reformulate(c("OX", "MYC", cov)), data = dd)
}
NEST <- tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                           estimator = EST_ALL)
NEST$delta_OX <- NA_real_
NEST$delta_MYC <- NA_real_
NEST$identity_dev <- NA_real_
for (r in seq_len(nrow(NEST))) {
  coh <- NEST$cohort[r]
  inst <- NEST$instrument[r]
  est <- NEST$estimator[r]
  X1 <- .design(DSET[[coh]], inst, est, RUNG_COV[[coh]][["m1"]])
  dl <- stats::lm.fit(X1, DSET[[coh]]$PROLIF)$coefficients
  e  <- EXT[EXT$cohort == coh & EXT$instrument == inst & EXT$estimator == est, ]
  NEST$delta_OX[r]  <- dl[["OX"]]
  NEST$delta_MYC[r] <- dl[["MYC"]]
  NEST$identity_dev[r] <- max(
    abs((e$b_MYC_m3 - e$b_MYC_m1) + e$b_PR * dl[["MYC"]]),
    abs((e$b_OX_m3 - e$b_OX_m1) + e$b_PR * dl[["OX"]]))
}
NEST <- NEST %>% dplyr::mutate(delta_MYC_positive = delta_MYC > 0)
.assert(paste0("G-A4: every gene fit read is nested, m1 in m3, on one sample ",
               "(the omitted-variable identity holds)"),
        max(NEST$identity_dev) < TOL, n = nrow(EXT) * 2L,
        max_dev = max(NEST$identity_dev))
message("   delta (MYC's coefficient on proliferation, m1 terms) positive in ",
        sum(NEST$delta_MYC_positive), " of ", nrow(NEST),
        " cells; section 11 assumes it is positive.")

# =============================================================================
# 4. The primary cells, the labels and the agreement check (section 7)
# =============================================================================
message("\n4. the reading rule (section 7)")

PRIM <- EXT %>%
  dplyr::filter(estimator == EST_READ, readout %in% READ3) %>%
  dplyr::mutate(
    block = paste(cohort, instrument),
    ox_excl = ox_lo_m3 > 0 | ox_hi_m3 < 0,
    pr_opposite = ox_excl & ((b_OX_m3 > 0 & pr_hi < 0) |
                               (b_OX_m3 < 0 & pr_lo > 0)),
    diff_excl = diff_lo > 0 | diff_hi < 0,
    fit_id = .fid(cohort, instrument, EST_READ, "m3", readout)) %>%
  dplyr::arrange(match(readout, READ3), match(block, BLOCK_ORDER))
.assert("12 primary cells, M_a at m3, none void",
        nrow(PRIM) == 12L && !any(PRIM$void), n = nrow(PRIM))

# Section 7.1, transcribed. Read only for the three genes.
.label_gene <- function(rd) {
  p <- PRIM[PRIM$readout == rd, ]
  if (all(p$pr_opposite & p$diff_excl)) return("OPPOSITE")
  if (all(p$diff_excl) && length(unique(sign(p$diff))) == 1L) return("DISTINCT")
  "NOT SEPARATED"
}
# Section 7.2: proliferation's sign on the least-entangled score must equal
# its sign on M_a in all 12 cells. Sign only.
AGREE <- PRIM %>%
  dplyr::select(readout, cohort, instrument, block, b_PR_ref = b_PR) %>%
  dplyr::left_join(EXT %>% dplyr::filter(estimator == EST_AGREE,
                                         readout %in% READ3) %>%
                     dplyr::select(readout, cohort, instrument,
                                   b_PR_clean = b_PR),
                   by = c("readout", "cohort", "instrument")) %>%
  dplyr::mutate(agree = sign(b_PR_ref) == sign(b_PR_clean))
.assert("the agreement check covers all 12 cells", nrow(AGREE) == 12L &&
          !anyNA(AGREE$agree), n = nrow(AGREE))
LABELS <- tibble::tibble(
  readout = READ3,
  label = vapply(READ3, .label_gene, character(1)),
  cells_agree = vapply(READ3, function(r) sum(AGREE$agree[AGREE$readout == r]),
                       integer(1))) %>%
  dplyr::mutate(estimator_dependent = cells_agree < 4L,
                printed = ifelse(estimator_dependent,
                                 paste0(label, ", ESTIMATOR-DEPENDENT"), label))

# Section 8: rows 1-3 are exclusive on the labels; row 4 applies in addition.
TEXT_ROW <- list(
  row = if (all(LABELS$label == "OPPOSITE")) 1L else
    if (!any(LABELS$label == "NOT SEPARATED") &&
        any(LABELS$label == "DISTINCT")) 2L else 3L,
  row4_estimator_dependent = any(LABELS$estimator_dependent),
  genes_other_way = LABELS$readout[LABELS$label == "OPPOSITE"],
  genes_descriptive_sentence = LABELS$readout[LABELS$label != "OPPOSITE"])

# VIF (E47's, every read cell) and G3 bands, printed beside the cells.
VIF_READ <- e47$vif %>%
  dplyr::filter(fit_id %in% PRIM$fit_id, term %in% c("OX", "MYC", "PROLIF")) %>%
  dplyr::select(fit_id, term, gvif_adj) %>%
  tidyr::pivot_wider(names_from = term, values_from = gvif_adj,
                     names_prefix = "vif_")
G3_READ <- G3 %>%
  dplyr::filter(estimator == EST_READ) %>%
  dplyr::select(cohort, instrument, rho_MYC_PROLIF, rho_OX_PROLIF,
                band_MYC_PROLIF)
PRIM <- PRIM %>%
  dplyr::left_join(VIF_READ, by = "fit_id") %>%
  dplyr::left_join(G3_READ, by = c("cohort", "instrument"))
.assert("E47's VIF and G3 band are attached to every read cell",
        !anyNA(PRIM[, c("vif_OX", "vif_MYC", "vif_PROLIF", "band_MYC_PROLIF")]),
        n = nrow(PRIM))

message("   the twelve primary cells (M_a, m3); proliferation is a post-hoc ",
        "characterisation of a declared covariate:")
.pr(PRIM %>% dplyr::transmute(
  gene = readout, block,
  OXPHOS = .fmt3(b_OX_m3, ox_lo_m3, ox_hi_m3),
  proliferation = .fmt3(b_PR, pr_lo, pr_hi),
  difference = .fmt3(diff, diff_lo, diff_hi),
  vif_MYC = round(vif_MYC, 2), vif_PROLIF = round(vif_PROLIF, 2),
  MYC_PROLIF_band = band_MYC_PROLIF))
message("\n   labels (section 7.1) with the agreement check (7.2):")
.pr(LABELS)
message("   section 8: outcome row ", TEXT_ROW$row,
        if (TEXT_ROW$row4_estimator_dependent) ", plus row 4" else "",
        ". Genes needing a descriptive sentence: ",
        if (length(TEXT_ROW$genes_descriptive_sentence))
          paste(TEXT_ROW$genes_descriptive_sentence, collapse = ", ")
        else "none", ".")

# Section 7.4: the nine other genes, a descriptive tally fixed beforehand.
NINE_TALLY <- EXT %>%
  dplyr::filter(estimator == EST_READ, readout %in% NINE) %>%
  dplyr::mutate(both_opposite = (ox_lo_m3 > 0 & pr_hi < 0) |
                  (ox_hi_m3 < 0 & pr_lo > 0)) %>%
  dplyr::group_by(readout) %>%
  dplyr::summarise(cells = dplyr::n(), cells_opposite = sum(both_opposite),
                   .groups = "drop") %>%
  dplyr::mutate(counts_toward_n = cells_opposite == 4L) %>%
  dplyr::arrange(match(readout, NINE))
N_TALLY <- sum(NINE_TALLY$counts_toward_n)
.assert("the nine-gene tally covers nine genes in four cells each",
        nrow(NINE_TALLY) == 9L && all(NINE_TALLY$cells == 4L), n = 36L)
message("   section 7.4, descriptive: n = ", N_TALLY, " of 9")

# =============================================================================
# 5. Part B: unavailable (section 9)
# =============================================================================
message("\n5. Part B: unavailable")

.fit_like <- function(o) inherits(o, c("lm", "glm", "coxph"))
e40_fits <- any(vapply(e40, function(o) .fit_like(o) ||
                         (is.list(o) && !is.data.frame(o) && length(o) > 0L &&
                            any(vapply(o, .fit_like, logical(1)))),
                       logical(1)))
.assert("Part B: E40 saves no fitted model", !e40_fits, n = length(e40))
.assert("Part B: E40 saves no variance matrix",
        !any(grepl("vcov", names(e40))), n = length(e40))
E43_TERMS <- vapply(e43$fits, function(f)
  paste(attr(stats::terms(f), "term.labels"), collapse = " + "), character(1))
.assert("Part B: E43's fits are m1 only, with no MYC term",
        length(E43_TERMS) == 8L && !any(grepl("\\bMYC\\b", E43_TERMS)),
        n = length(E43_TERMS))
E39G <- e39$gate %>%
  dplyr::filter(pair == "MYC vs PROLIF_DISJOINT") %>%
  dplyr::select(cohort, n, rho) %>%
  dplyr::mutate(band = .band(rho))
.assert("ED-a: E39's three gate values equal the declaration's, to 4 dp",
        nrow(E39G) == 3L &&
          all(round(E39G$rho, 4) == unname(E39_DECLARED[E39G$cohort])),
        n = 3L,
        max_dev = .dev(round(E39G$rho, 4), unname(E39_DECLARED[E39G$cohort])))
PART_B <- list(
  status = "unavailable",
  reason = paste0(
    "No saved object holds the m2 fits' MYC coefficient or variance matrix: ",
    "E40 kept the OXPHOS row only, with no fitted model or vcov; E43 kept ",
    "the m1 fits only, with no MYC term. Under brief 5.9.3 nothing is ",
    "refitted; a refit needs its own declaration."),
  e40_elements = names(e40),
  e43_fit_terms = E43_TERMS,
  myc_score = paste0("M-a, the stripped Felsher estimator, from the ",
                     "validation study's H4 frames (E39 declaration, ",
                     "section 4)"),
  e39_gate = E39G,
  text = paste0("The optional pCR MYC sentence (brief section 8) leaves the ",
                "text. E40's MYC-adjusted OR, 0.87, stays."))

# =============================================================================
# 6. Part C: C1, C3 and C4 - reported, never read (section 10)
# =============================================================================
message("\n6. Part C: C1, C3, C4 (reported, never read)")

.fit_lm <- function(d, rd, inst, est, cov) {
  dd <- data.frame(y = d[[rd]], OX = d[[paste0("OX_", inst)]],
                   MYC = d[[paste0("MYC_", est)]])
  for (cv in cov) dd[[cv]] <- d[[cv]]
  for (fc in intersect(cov, c("PAM50", "TP53_status"))) {
    dd[[fc]] <- droplevels(dd[[fc]])
  }
  .stop_if(!anyNA(dd), "NA inside a fit: ", rd)
  stats::lm(stats::reformulate(c("OX", "MYC", cov), response = "y"), data = dd)
}
.tidy_fit <- function(f, id) {
  s  <- summary(f)$coefficients
  tq <- stats::qt(CI_Q, f$df.residual)
  tibble::tibble(fit_id = id, term = rownames(s), estimate = s[, 1],
                 se = s[, 2], t = s[, 3], p = s[, 4],
                 ci_lo = s[, 1] - tq * s[, 2], ci_hi = s[, 1] + tq * s[, 2],
                 n = stats::nobs(f), df = f$df.residual)
}
.vif_tidy <- function(f, id) {
  v <- car::vif(f)
  out <- if (is.matrix(v)) {
    tibble::tibble(term = rownames(v), gvif = v[, "GVIF"], df = v[, "Df"],
                   gvif_adj = v[, "GVIF^(1/(2*Df))"])
  } else {
    tibble::tibble(term = names(v), gvif = unname(v), df = 1,
                   gvif_adj = sqrt(unname(v)))
  }
  out %>% dplyr::mutate(fit_id = id, .before = 1)
}
# A refit must reproduce E47's stored fit before its variant is reported.
.vs_e47 <- function(f, fid) {
  s <- .get_fit(fid)
  b <- stats::coef(f)
  c(same_terms = setequal(names(b), names(s$b)),
    coef_dev = .dev(b[names(s$b)], s$b),
    vcov_dev = max(abs(stats::vcov(f)[names(s$b), names(s$b)] - s$V)))
}
.report <- function(f, spec, extra = list()) {
  b  <- stats::coef(f)
  V  <- stats::vcov(f)
  ct <- .contrast(b, V, f$df.residual)
  td <- .tidy_fit(f, "x")
  w  <- function(tm, col) td[[col]][td$term == tm]
  tibble::tibble(spec = spec, n = stats::nobs(f), df = f$df.residual,
                 b_OX = w("OX", "estimate"), ox_lo = w("OX", "ci_lo"),
                 ox_hi = w("OX", "ci_hi"),
                 b_PR = w("PROLIF", "estimate"), pr_lo = w("PROLIF", "ci_lo"),
                 pr_hi = w("PROLIF", "ci_hi"),
                 diff = ct[["diff"]], diff_lo = ct[["diff_lo"]],
                 diff_hi = ct[["diff_hi"]], !!!extra)
}

# --- C1, TP53 status (10.1) --------------------------------------------------
d_t <- DSET$TCGA
cvj <- CV[match(d_t$sample_id, CV$patient), ]
.assert("C1: every TCGA tumour joins the covariate table on the patient barcode",
        sum(!is.na(cvj$patient)) == EXPECT_N[["TCGA"]], n = nrow(d_t))
pam_cv <- sub("^BRCA_", "", cvj$PAM50)
pam_cv[pam_cv == "Her2"] <- "HER2"
.assert("C1: PAM50 agrees across the join, after the label prefix",
        identical(pam_cv, as.character(d_t$PAM50)), n = nrow(d_t))
tp <- cvj$TP53_status
.assert("C1: TP53 status is complete in the 938 tumours", !anyNA(tp),
        n = length(tp))
tp_tab <- table(factor(tp, levels = names(TP53_EXPECT)))
.assert("C1: 325 mutant and 613 wild-type, as the declaration records",
        identical(as.integer(tp_tab), unname(TP53_EXPECT)), n = length(tp),
        detail = paste(names(tp_tab), as.integer(tp_tab), collapse = ", "))
DSET$TCGA$TP53_status <- factor(tp, levels = c("wildtype", "mutant"))

C1 <- list(); C1_COEF <- list(); C1_VCOV <- list(); C1_VIF <- list()
C1_DIG <- list(); FITS_C1 <- list(); C1_GATE <- list()
for (inst in names(INSTR)) for (rd in READ3) {
  fid <- .fid("TCGA", inst, EST_READ, "m3", rd)
  g0 <- .vs_e47(.fit_lm(DSET$TCGA, rd, inst, EST_READ, RUNG_COV$TCGA$m3), fid)
  C1_GATE[[fid]] <- g0
  f1 <- .fit_lm(DSET$TCGA, rd, inst, EST_READ,
                c(RUNG_COV$TCGA$m3, "TP53_status"))
  id <- paste("C1", fid, "+ TP53_status", sep = " | ")
  td <- .tidy_fit(f1, id)
  C1[[id]] <- .report(f1, "with TP53 status", list(
    instrument = inst, readout = rd, estimator = EST_READ,
    b_TP53 = td$estimate[td$term == "TP53_statusmutant"],
    tp53_lo = td$ci_lo[td$term == "TP53_statusmutant"],
    tp53_hi = td$ci_hi[td$term == "TP53_statusmutant"]))
  C1_COEF[[id]] <- td
  C1_VCOV[[id]] <- stats::vcov(f1)
  C1_VIF[[id]]  <- .vif_tidy(f1, id)
  C1_DIG[[id]]  <- tibble::tibble(fit_id = id,
                                  digest_coef = digest::digest(stats::coef(f1)),
                                  digest_vcov = digest::digest(stats::vcov(f1)))
  FITS_C1[[id]] <- f1
}
C1_GATE <- dplyr::bind_rows(lapply(C1_GATE, function(z) tibble::as_tibble(as.list(z))))
.assert(paste0("C1 gate: without the TP53 term, the six TCGA m3 fits reproduce ",
               "E47's coefficients and variance matrices exactly"),
        all(C1_GATE$same_terms == 1) && max(C1_GATE$coef_dev) < TOL &&
          max(C1_GATE$vcov_dev) < TOL, n = nrow(C1_GATE),
        max_dev = max(C1_GATE$coef_dev, C1_GATE$vcov_dev))
C1 <- dplyr::bind_rows(C1)

# --- C4, TCGA with SCAN-B's spine (10.4) -------------------------------------
C4 <- list(); C4_COEF <- list(); C4_VCOV <- list(); C4_VIF <- list()
C4_DIG <- list(); FITS_C4 <- list(); C4_GATE <- list()
for (est in c(EST_READ, EST_AGREE)) for (inst in names(INSTR)) for (rd in READ3) {
  fid <- .fid("TCGA", inst, est, "m3", rd)
  g0 <- .vs_e47(.fit_lm(DSET$TCGA, rd, inst, est, RUNG_COV$TCGA$m3), fid)
  C4_GATE[[fid]] <- g0
  f4 <- .fit_lm(DSET$TCGA, rd, inst, est, C4_COV)
  id <- paste("C4", fid, "without purity and leukocyte fraction", sep = " | ")
  C4[[id]] <- .report(f4, "without purity and leukocyte fraction", list(
    instrument = inst, readout = rd, estimator = est))
  C4_COEF[[id]] <- .tidy_fit(f4, id)
  C4_VCOV[[id]] <- stats::vcov(f4)
  C4_VIF[[id]]  <- .vif_tidy(f4, id)
  C4_DIG[[id]]  <- tibble::tibble(fit_id = id,
                                  digest_coef = digest::digest(stats::coef(f4)),
                                  digest_vcov = digest::digest(stats::vcov(f4)))
  FITS_C4[[id]] <- f4
}
C4_GATE <- dplyr::bind_rows(lapply(C4_GATE, function(z) tibble::as_tibble(as.list(z))))
.assert(paste0("C4 gate: with purity and leukocyte fraction restored, the 12 ",
               "TCGA fits reproduce E47's m3 coefficients and variance ",
               "matrices exactly"),
        nrow(C4_GATE) == 12L && all(C4_GATE$same_terms == 1) &&
          max(C4_GATE$coef_dev) < TOL && max(C4_GATE$vcov_dev) < TOL,
        n = nrow(C4_GATE), max_dev = max(C4_GATE$coef_dev, C4_GATE$vcov_dev))
C4 <- dplyr::bind_rows(C4)

# The primary beside C1 and C4, TCGA only (S-1's rows).
SENS <- dplyr::bind_rows(
  EXT %>% dplyr::filter(cohort == "TCGA", readout %in% READ3,
                        estimator %in% c(EST_READ, EST_AGREE)) %>%
    dplyr::transmute(spec = "as fitted", instrument, readout, estimator,
                     n, df = df_m3, b_OX = b_OX_m3, ox_lo = ox_lo_m3,
                     ox_hi = ox_hi_m3, b_PR, pr_lo, pr_hi, diff, diff_lo,
                     diff_hi),
  C1 %>% dplyr::select(spec, instrument, readout, estimator, n, df, b_OX, ox_lo,
                       ox_hi, b_PR, pr_lo, pr_hi, diff, diff_lo, diff_hi),
  C4 %>% dplyr::select(spec, instrument, readout, estimator, n, df, b_OX, ox_lo,
                       ox_hi, b_PR, pr_lo, pr_hi, diff, diff_lo, diff_hi))
.assert("S-1 rows: 3 genes x 2 measures x (3 specs on M_a + 2 on the clean score)",
        nrow(SENS) == 30L, n = nrow(SENS))
message("   C1 and C4 (TCGA), REPORTED, NEVER READ:")
.pr(SENS %>% dplyr::arrange(readout, instrument, estimator, spec) %>%
      dplyr::transmute(readout, instrument, estimator, spec,
                       proliferation = .fmt3(b_PR, pr_lo, pr_hi),
                       difference = .fmt3(diff, diff_lo, diff_hi)))

# --- C3, "more steeply" at m3 (10.3): E47 section 11's code, copied ---------
# The only change from E47 is the design: RUNG_COV[[coh]][[rung]] instead of
# m1 alone. The seeds are E47's, PROJECT_SEED + k in E47's cell order, so
# each m3 replicate resamples the same tumours as E47's m1 replicate k.
SEC9 <- list()
BOOT_M3 <- list()
k <- 0L
for (coh in COHORTS) for (inst in names(INSTR)) {
  k <- k + 1L
  d <- DSET[[coh]]
  for (rung in c("m1", "m3")) {
    X  <- .design(d, inst, EST_READ, RUNG_COV[[coh]][[rung]])
    y1 <- d$BCL2L1
    y2 <- d$BBC3
    b1 <- stats::lm.fit(X, y1)$coefficients[["OX"]]
    b2 <- stats::lm.fit(X, y2)$coefficients[["OX"]]
    s1 <- .get_fit(.fid(coh, inst, EST_READ, rung, "BCL2L1"))
    s2 <- .get_fit(.fid(coh, inst, EST_READ, rung, "BBC3"))
    .assert(paste0("C3, ", coh, " ", inst, " ", rung, ": the bootstrap's design ",
                   "matrix reproduces E47's stored OX coefficients"),
            abs(b1 - s1$b[["OX"]]) < TOL && abs(b2 - s2$b[["OX"]]) < TOL,
            n = 2L, max_dev = max(abs(b1 - s1$b[["OX"]]),
                                  abs(b2 - s2$b[["OX"]])))
    set.seed(PROJECT_SEED + k)
    n <- nrow(X)
    reps <- numeric(BOOT_N)
    for (b in seq_len(BOOT_N)) {
      idx <- sample.int(n, n, replace = TRUE)
      Xb  <- X[idx, , drop = FALSE]
      reps[b] <- stats::lm.fit(Xb, y1[idx])$coefficients[["OX"]] -
                 stats::lm.fit(Xb, y2[idx])$coefficients[["OX"]]
    }
    .stop_if(!anyNA(reps), "a bootstrap replicate is NA: ", coh, " ", inst,
             " ", rung)
    ci <- stats::quantile(reps, c(0.025, 0.975), names = FALSE, type = 7)
    if (rung == "m1") {
      ref <- e47$boot_replicates[[paste(coh, inst)]]
      .assert(paste0("C3, ", coh, " ", inst, ": the copied code at m1 ",
                     "reproduces E47's 2,000 replicates exactly"),
              identical(reps, ref), n = BOOT_N, max_dev = .dev(reps, ref))
      r9 <- e47$bcl2l1_bbc3[e47$bcl2l1_bbc3$cohort == coh &
                              e47$bcl2l1_bbc3$instrument == inst, ]
      dev9 <- max(abs((b1 - b2) - r9$delta), abs(ci[1] - r9$ci_lo),
                  abs(ci[2] - r9$ci_hi))
      .assert(paste0("C3, ", coh, " ", inst, ": the copied code at m1 ",
                     "reproduces E47's section 9 row"), dev9 < TOL, n = 3L,
              max_dev = dev9)
    } else {
      BOOT_M3[[paste(coh, inst)]] <- reps
    }
    SEC9[[length(SEC9) + 1L]] <- tibble::tibble(
      cohort = coh, instrument = inst, rung = rung, n = n,
      seed = PROJECT_SEED + k,
      b_OX_BCL2L1 = s1$b[["OX"]], lo_BCL2L1 = s1$lo[["OX"]],
      hi_BCL2L1 = s1$hi[["OX"]],
      b_OX_BBC3 = s2$b[["OX"]], lo_BBC3 = s2$lo[["OX"]], hi_BBC3 = s2$hi[["OX"]],
      delta = b1 - b2, boot_se = stats::sd(reps), ci_lo = ci[1], ci_hi = ci[2],
      steeper = (b1 - b2) > 0 && ci[1] > 0)
  }
}
C3 <- dplyr::bind_rows(SEC9)
C3_HOLDS_M3 <- all(C3$steeper[C3$rung == "m3"])
message("   C3, REPORTED: 'more steeply' at m3 by E47 section 9's rule, applied ",
        "descriptively: ", if (C3_HOLDS_M3) "holds in all four cells" else
          "does not hold in all four cells", ".")
PART_C2 <- list(
  status = "NOT RUN, a declared omission (10.2)",
  reason = paste0(
    "PROLIF_STD is PROLIF_DISJOINT plus nine genes (CTPS1, DNMT1, HMGA1, NCL, ",
    "NOP56, PRMT5, RANBP1, TFDP1, UCK2), all nine in M_a ",
    "(FELSHER__MITOSTRIP). With M_a in the model, C2 would test the ",
    "disjoint construction, not proliferation."))

# =============================================================================
# 7. Descriptive correlations (quantity 4) for ED-a
# =============================================================================
message("\n7. correlations")

VARS <- c(OX_gsva = "OXPHOS level", OX_mitopps = "OXPHOS priority",
          MYC_clean = unname(EST_NAME[["clean"]]),
          MYC_M_a = unname(EST_NAME[["M_a"]]),
          MYC_entangled = unname(EST_NAME[["entangled"]]),
          MYC_regulon = unname(EST_NAME[["regulon"]]),
          PROLIF = "Proliferation")
CORR <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  m <- stats::cor(as.matrix(DSET[[coh]][, names(VARS)]), method = "spearman")
  ij <- which(lower.tri(m), arr.ind = TRUE)
  tibble::tibble(cohort = coh, var_row = rownames(m)[ij[, 1]],
                 var_col = colnames(m)[ij[, 2]], rho = m[ij])
})) %>% dplyr::mutate(band = .band(rho))
.rho <- function(coh, a, b) {
  r <- CORR$rho[CORR$cohort == coh & ((CORR$var_row == a & CORR$var_col == b) |
                                        (CORR$var_row == b & CORR$var_col == a))]
  .stop_if(length(r) == 1L, "pair not unique: ", a, " ", b)
  r
}
G3_CHK <- G3 %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    d_ox_myc = abs(rho_OX_MYC - .rho(cohort, paste0("OX_", instrument),
                                     paste0("MYC_", estimator))),
    d_myc_pr = abs(rho_MYC_PROLIF - .rho(cohort, paste0("MYC_", estimator),
                                         "PROLIF")),
    d_ox_pr = abs(rho_OX_PROLIF - .rho(cohort, paste0("OX_", instrument),
                                       "PROLIF"))) %>%
  dplyr::ungroup()
G3_DEV <- max(G3_CHK$d_ox_myc, G3_CHK$d_myc_pr, G3_CHK$d_ox_pr)
.assert("the g3 pairs, recomputed on E47's frame, reproduce E47's g3",
        G3_DEV < 1e-12, n = 3L * nrow(G3), max_dev = G3_DEV)
.assert("ED-a holds 21 pairs in each cohort", nrow(CORR) == 42L, n = 42L)

# =============================================================================
# 8. Displays (section 12), content fixed before any quantity was read
# =============================================================================
message("\n8. displays")

KEYS <- tibble::tribble(
  ~key,             ~text,                                           ~col,     ~fill,    ~shape,
  "ox_without",     "OXPHOS, without the proliferation term",        COL_OX,   "white",  21,
  "ox_with",        "OXPHOS, with the proliferation term",           COL_OX,   COL_OX,   21,
  "pr",             "Proliferation",                                 COL_PR,   COL_PR,   24,
  "myc_without",    "MYC activity, without the proliferation term",  COL_MYC,  "white",  22,
  "myc_with",       "MYC activity, with the proliferation term",     COL_MYC,  COL_MYC,  22,
  "bcl2l1_without", "BCL2L1, without the proliferation term",        COL_OX,   "white",  21,
  "bcl2l1_with",    "BCL2L1, with the proliferation term",           COL_OX,   COL_OX,   21,
  "bbc3_without",   "BBC3, without the proliferation term",          COL_OX,   "white",  23,
  "bbc3_with",      "BBC3, with the proliferation term",             COL_OX,   COL_OX,   23,
  "diff_without",   "Difference, without the proliferation term",    "black",  "white",  21,
  "diff_with",      "Difference, with the proliferation term",       "black",  "black",  21,
  "diff",           "Proliferation minus OXPHOS",                    "black",  "black",  21)
K_TEXT  <- stats::setNames(KEYS$text, KEYS$key)
K_COL   <- stats::setNames(KEYS$col, KEYS$key)
K_FILL  <- stats::setNames(KEYS$fill, KEYS$key)
K_SHAPE <- stats::setNames(KEYS$shape, KEYS$key)
CONN_PAIRS <- list(OX = c("ox_without", "ox_with"),
                   MYC = c("myc_without", "myc_with"),
                   BCL2L1 = c("bcl2l1_without", "bcl2l1_with"),
                   BBC3 = c("bbc3_without", "bbc3_with"),
                   DIFF = c("diff_without", "diff_with"))
CONN_COL <- c(OX = COL_OX, MYC = COL_MYC, BCL2L1 = COL_OX, BBC3 = COL_OX,
              DIFF = "black")

.theme_forest <- function() {
  ggplot2::theme_classic(base_size = BASE_PT) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.3),
      axis.text.x = ggplot2::element_text(size = MIN_PT, colour = "black"),
      axis.title.x = ggplot2::element_text(size = BASE_PT - 0.5),
      plot.title = ggplot2::element_text(size = BASE_PT, face = "bold",
                                         hjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = BASE_PT),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.title = ggplot2::element_text(size = BASE_PT - 0.5),
      plot.caption = ggplot2::element_text(size = MIN_PT),
      legend.key.width = grid::unit(3, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(2, 4, 2, 4))
}
.theme_text <- function() {
  ggplot2::theme_void(base_size = BASE_PT) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.title = ggplot2::element_text(size = BASE_PT - 0.5),
      strip.text = ggplot2::element_text(size = BASE_PT),
      plot.caption = ggplot2::element_text(size = MIN_PT),
      plot.title = ggplot2::element_text(size = BASE_PT, face = "bold"),
      plot.margin = ggplot2::margin(2, 1, 2, 1))
}

# The axis rule (12.1): symmetric about the null, containing every CI,
# rounded outward to 0.1 (on the log10 scale for ratios). A rule, not a value.
.sym_lim <- function(v) {
  L <- ceiling(max(abs(v)) * 10 - 1e-9) / 10
  c(-L, L)
}
.sym_breaks <- function(L) {
  st <- if (L <= 0.6) 0.2 else if (L <= 1.2) 0.4 else 0.5
  k <- floor(L / st + 1e-9)
  round(seq(-k, k) * st, 10)
}
.log_lim <- function(v) {
  L <- ceiling(max(abs(log10(v))) * 10 - 1e-9) / 10
  10^c(-L, L)
}
BRK_OR <- c(0.33, 0.5, 0.75, 1, 1.5, 2, 3)        # fixed sets, not data-chosen
BRK_HR <- c(0.67, 0.8, 1, 1.25, 1.5)

.connectors <- function(pf, id_cols) {
  out <- list()
  for (nm in names(CONN_PAIRS)) {
    a <- pf[as.character(pf$key) == CONN_PAIRS[[nm]][1], ]
    b <- pf[as.character(pf$key) == CONN_PAIRS[[nm]][2], ]
    if (!nrow(a) && !nrow(b)) next
    ida <- do.call(paste, c(a[id_cols], sep = " | "))
    idb <- do.call(paste, c(b[id_cols], sep = " | "))
    m <- match(ida, idb)
    .stop_if(nrow(a) == nrow(b) && !anyNA(m) && !any(duplicated(m)),
             "a connector cannot be paired: ", nm)
    out[[nm]] <- tibble::tibble(term = nm, id = ida, x = a$est, xend = b$est[m],
                                y = a$y_mark, yend = b$y_mark[m])
    if ("panel" %in% names(pf)) out[[nm]]$panel <- a$panel
  }
  dplyr::bind_rows(out)
}

# One forest builder for every E48 panel: E45's grammar.
.forest <- function(pf, cn, xlim, ylim, sep_y, y_cue, x_title, keys,
                    cue = c("lower expression", "higher expression"),
                    legend = TRUE, log_x = FALSE, breaks = NULL, size = 1.5,
                    facet = FALSE, title = NULL, legend_ncol = NULL) {
  null_x <- if (log_x) 1 else 0
  p <- ggplot2::ggplot(pf) +
    ggplot2::geom_vline(xintercept = null_x, colour = "grey45", linewidth = 0.3)
  if (length(sep_y)) {
    p <- p + ggplot2::geom_hline(yintercept = sep_y, colour = "grey85",
                                 linewidth = 0.25)
  }
  if (!is.null(cn) && nrow(cn)) {
    for (tm in unique(cn$term)) {
      p <- p + ggplot2::geom_segment(
        data = cn[cn$term == tm, ],
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        colour = unname(CONN_COL[tm]), linewidth = 0.3, alpha = 0.6,
        inherit.aes = FALSE)
    }
  }
  if (log_x) {
    a_in  <- c(1 / 1.04, 1.04)
    a_out <- c(xlim[1] * 1.10, xlim[2] / 1.10)
  } else {
    w <- diff(xlim)
    a_in  <- c(-0.025 * w, 0.025 * w)
    a_out <- c(xlim[1] + 0.06 * w, xlim[2] - 0.06 * w)
  }
  p <- p +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark,
                                       yend = y_mark, colour = key,
                                       alpha = void),
                          linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = est, y = y_mark, colour = key,
                                     fill = key, shape = key, alpha = void),
                        size = size, stroke = 0.45)
  # Direction cues describe the axis. A difference panel has no expression
  # direction, so it is drawn with cue = NULL and its axis title carries it.
  if (!is.null(cue)) {
    p <- p +
      ggplot2::annotate("segment", x = a_in[1], xend = a_out[1], y = y_cue,
                        yend = y_cue, linewidth = 0.3, colour = "grey30",
                        arrow = grid::arrow(length = grid::unit(1.1, "mm"),
                                            type = "closed")) +
      ggplot2::annotate("segment", x = a_in[2], xend = a_out[2], y = y_cue,
                        yend = y_cue, linewidth = 0.3, colour = "grey30",
                        arrow = grid::arrow(length = grid::unit(1.1, "mm"),
                                            type = "closed")) +
      ggplot2::annotate("text", x = a_in[1], y = y_cue + 0.5, hjust = 1,
                        label = cue[1], size = SUB_MM, colour = "grey30") +
      ggplot2::annotate("text", x = a_in[2], y = y_cue + 0.5, hjust = 0,
                        label = cue[2], size = SUB_MM, colour = "grey30")
  }
  xs <- if (log_x) {
    ggplot2::scale_x_continuous(transform = "log10", limits = xlim,
                                breaks = breaks, labels = format(breaks),
                                expand = ggplot2::expansion(mult = 0.02))
  } else {
    ggplot2::scale_x_continuous(limits = xlim, breaks = breaks,
                                expand = ggplot2::expansion(mult = 0.02))
  }
  p <- p + xs +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = K_COL, breaks = keys,
                                 labels = K_TEXT[keys], drop = FALSE,
                                 name = NULL) +
    ggplot2::scale_fill_manual(values = K_FILL, breaks = keys,
                               labels = K_TEXT[keys], drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = K_SHAPE, breaks = keys,
                                labels = K_TEXT[keys], drop = FALSE,
                                name = NULL) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.3),
                                guide = "none") +
    ggplot2::labs(x = x_title, title = title) +
    .theme_forest()
  if (facet) {
    p <- p + ggplot2::facet_grid(. ~ panel, labeller = ggplot2::label_parsed) +
      ggplot2::theme(panel.spacing.x = grid::unit(5, "mm"))
  }
  if (!legend) {
    p <- p + ggplot2::guides(colour = "none", fill = "none", shape = "none")
  } else if (!is.null(legend_ncol)) {
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE),
      fill = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE),
      shape = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE))
  }
  p
}
# Three ticks for a narrow panel, so the labels never touch.
.coarse_breaks <- function(L) {
  k <- floor(L / 0.4 + 1e-9) * 0.4
  if (k <= 0) round(c(-L, 0, L), 10) else round(c(-k, 0, k), 10)
}

# A text column sharing the forest's y scale.
.text_panel <- function(ylim, layers) {
  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    .theme_text()
  for (l in layers) {
    p <- p + ggplot2::geom_text(
      data = l$d, ggplot2::aes(x = x, y = y, label = label),
      hjust = l$hjust, size = l$size, parse = isTRUE(l$parse),
      fontface = if (is.null(l$face)) "plain" else l$face,
      colour = if (is.null(l$colour)) "black" else l$colour,
      lineheight = 0.9, inherit.aes = FALSE)
  }
  p
}
.lyr <- function(x, y, label, hjust = 0.5, size = TXT_MM, parse = FALSE,
                 face = "plain", colour = "black") {
  list(d = tibble::tibble(x = x, y = y, label = label), hjust = hjust,
       size = size, parse = parse, face = face, colour = colour)
}

# --- the checks every figure must pass (section 14.2) ------------------------
FIG_CHECKS <- list()
.check_points <- function(p, pf, label, log_x = FALSE) {
  b   <- ggplot2::ggplot_build(p)
  geo <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  xe  <- vapply(p$layers, function(l) {
    m <- l$mapping$xend
    if (is.null(m)) "" else rlang::as_label(m)
  }, character(1))
  pts <- dplyr::bind_rows(lapply(which(geo == "GeomPoint"), function(k) {
    dk <- b$data[[k]]
    tibble::tibble(x = dk$x, y = dk$y, PANEL = as.integer(dk$PANEL),
                   fill = as.character(dk$fill), colour = as.character(dk$colour),
                   alpha = if ("alpha" %in% names(dk)) dk$alpha else NA_real_)
  }))
  cis <- dplyr::bind_rows(lapply(which(geo == "GeomSegment" & xe == "hi"),
                                 function(k) {
    dk <- b$data[[k]]
    tibble::tibble(x = dk$x, xend = dk$xend, y = dk$y,
                   PANEL = as.integer(dk$PANEL),
                   colour = as.character(dk$colour))
  }))
  panel_f <- if ("panel" %in% names(pf)) as.integer(pf$panel) else
    rep(1L, nrow(pf))
  key_f <- paste(panel_f, round(pf$y_mark, 6))
  mp <- match(paste(pts$PANEL, round(pts$y, 6)), key_f)
  mc <- match(paste(cis$PANEL, round(cis$y, 6)), key_f)
  .assert(paste0(label, ": every point and interval drawn is one stored row, ",
                 "none dropped"),
          nrow(pts) == nrow(pf) && !anyNA(mp) && !any(duplicated(mp)) &&
            nrow(cis) == nrow(pf) && !anyNA(mc) && !any(duplicated(mc)),
          n = nrow(pts), detail = paste(nrow(pts), "points,", nrow(cis),
                                        "intervals,", nrow(pf), "stored"))
  tf  <- if (log_x) log10 else identity
  dev <- max(.dev(pts$x, tf(pf$est[mp])), .dev(cis$x, tf(pf$lo[mc])),
             .dev(cis$xend, tf(pf$hi[mc])))
  .assert(paste0(label, ": drawn positions equal the stored values"),
          dev < 1e-12, n = nrow(pts) + 2L * nrow(cis), max_dev = dev)
  open <- grepl("without", as.character(pf$key[mp]))
  .assert(paste0(label, ": open = without the proliferation term, filled = ",
                 "with it"),
          all(pts$fill[open] == "white") && all(pts$fill[!open] != "white"),
          n = nrow(pts))
  v <- pf$void[mp]
  .assert(paste0(label, ": no void cell is drawn as read"),
          all(pts$alpha[v] < 1) && all(is.na(pts$alpha[!v]) | pts$alpha[!v] == 1),
          n = sum(v))
  used <- unique(c(pts$colour, cis$colour))
  data_cols <- setdiff(used, c("black"))
  .assert(paste0(label, ": every data colour drawn is from the validated ",
                 "palette"), all(data_cols %in% VALIDATED_PALETTE),
          n = length(used), detail = paste(used, collapse = ", "))
  invisible(list(points = pts, colours = used))
}
.check_axis <- function(pf, xlim, label, log_x = FALSE) {
  v  <- c(pf$lo, pf$hi)
  ok <- all(v >= xlim[1] & v <= xlim[2])
  sym <- if (log_x) abs(log10(xlim[1]) + log10(xlim[2])) < 1e-9 else
    abs(xlim[1] + xlim[2]) < 1e-12
  L <- if (log_x) log10(xlim[2]) else xlim[2]
  rounded <- abs(L * 10 - round(L * 10)) < 1e-9
  .assert(paste0(label, ": the axis is symmetric about the null, contains ",
                 "every interval and is rounded outward to 0.1"),
          ok && sym && rounded, n = length(v),
          detail = sprintf("limits %.4f to %.4f", xlim[1], xlim[2]))
}
.text_sizes <- function(p) {
  th  <- p$theme
  els <- c("axis.text.x", "axis.title.x", "legend.text", "legend.title",
           "plot.title", "plot.subtitle", "plot.caption", "strip.text.x",
           "strip.text.y", "axis.text.y", "axis.title.y")
  es <- vapply(els, function(e) {
    z <- tryCatch(ggplot2::calc_element(e, th), error = function(err) NULL)
    if (is.null(z) || inherits(z, "element_blank") || is.null(z$size)) NA_real_
    else as.numeric(z$size)
  }, numeric(1))
  ls <- unlist(lapply(p$layers, function(l) {
    if (inherits(l$geom, c("GeomText", "GeomLabel"))) {
      s <- l$aes_params$size
      if (is.null(s)) s <- l$geom$default_aes$size
      as.numeric(s) * ggplot2::.pt
    } else NULL
  }))
  c(es, ls)
}
BAN_RX <- "(^|[^A-Za-z0-9_])(m1|m3)([^A-Za-z0-9_]|$)|M_a|beta|Delta"
.panel_strings <- function(p) {
  labs <- tryCatch(ggplot2::get_labs(p), error = function(e) list())
  out <- unlist(labs[intersect(names(labs), c("title", "subtitle", "caption",
                                              "x", "y"))])
  for (l in p$layers) {
    if (!is.null(l$aes_params$label)) {
      out <- c(out, as.character(l$aes_params$label))
    }
    m <- l$mapping$label
    if (!is.null(m)) {
      dd <- if (is.data.frame(l$data)) l$data else p$data
      v <- tryCatch(rlang::eval_tidy(m, dd), error = function(e) NULL)
      out <- c(out, as.character(v))
    }
  }
  for (s in p$scales$scales) {
    if (is.character(s$labels)) out <- c(out, s$labels)
  }
  if ("panel" %in% names(p$data)) out <- c(out, levels(factor(p$data$panel)))
  as.character(out)
}
.check_parts <- function(parts, label, width_in) {
  sz <- unlist(lapply(parts, .text_sizes))
  .assert(paste0(label, ": no text below 5.5 pt"),
          all(sz[!is.na(sz)] >= MIN_PT - 1e-9), n = sum(!is.na(sz)),
          detail = sprintf("smallest %.2f pt", min(sz, na.rm = TRUE)))
  st <- unlist(lapply(parts, .panel_strings))
  bad <- st[grepl(BAN_RX, st)]
  .assert(paste0(label, ": no beta, m1, m3, M_a or Delta in any label"),
          length(bad) == 0L, n = length(st),
          detail = paste(utils::head(bad, 3), collapse = " / "))
  .assert(paste0(label, ": 178 mm wide"), abs(width_in - W_PAGE) < 1e-9, n = 1L)
}
.compose <- function(parts, widths, heights = NULL, ncol = length(parts)) {
  patchwork::wrap_plots(parts, ncol = ncol,
                        widths = widths / sum(widths) * W_PAGE,
                        heights = heights) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

# Plain-word row labels, drawn with plotmath so the source stays ASCII.
GENE_LAB  <- c(BBC3 = 'italic("BBC3")~"(PUMA)"',
               BCL2L1 = 'italic("BCL2L1")~"(BCL-xL)"',
               MCL1 = 'italic("MCL1")')
BLOCK_LAB <- c("TCGA gsva" = '"TCGA" %.% "level"',
               "TCGA mitopps" = '"TCGA" %.% "priority"',
               "SCAN-B gsva" = '"SCAN-B" %.% "level"',
               "SCAN-B mitopps" = '"SCAN-B" %.% "priority"')
HEAD_RHO_OP <- 'bold(rho*", OXPHOS-")'
HEAD_RHO_MP <- 'bold(rho*", MYC-")'
HEAD_DIFF   <- 'bold("Proliferation" - "OXPHOS")'
X_TITLE_SD  <- "Association with expression (SD per SD of each score)"
YH1 <- 1.30
YH2 <- 0.72

# --- the long frame of one coefficient set ----------------------------------
.long <- function(E, keys, rows, id_cols, keep = character(0)) {
  kc <- unique(c(id_cols, keep))
  parts <- list(
    ox_without  = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         term = "OX", rung = "m1",
                                         est = b_OX_m1, lo = ox_lo_m1,
                                         hi = ox_hi_m1, void = FALSE),
    ox_with     = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         term = "OX", rung = "m3",
                                         est = b_OX_m3, lo = ox_lo_m3,
                                         hi = ox_hi_m3, void = void),
    pr          = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         term = "PROLIF", rung = "m3",
                                         est = b_PR, lo = pr_lo, hi = pr_hi,
                                         void = void),
    myc_without = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         term = "MYC", rung = "m1",
                                         est = b_MYC_m1, lo = myc_lo_m1,
                                         hi = myc_hi_m1, void = FALSE),
    myc_with    = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         term = "MYC", rung = "m3",
                                         est = b_MYC_m3, lo = myc_lo_m3,
                                         hi = myc_hi_m3, void = void))
  dplyr::bind_rows(lapply(names(keys), function(k)
    parts[[k]] %>% dplyr::mutate(key = k))) %>%
    dplyr::left_join(rows, by = id_cols) %>%
    dplyr::mutate(y_mark = y_row + unname(keys[key]),
                  key = factor(key, levels = names(keys)))
}

# =============================================================================
# 8a. Fig. 4g, layouts 1 and 2 (12.2)
# =============================================================================
GENE_PITCH <- 5.6
ROWS12 <- tidyr::expand_grid(readout = READ3, block = BLOCK_ORDER) %>%
  dplyr::mutate(i = match(readout, READ3), j = match(block, BLOCK_ORDER),
                top = -(i - 1) * GENE_PITCH, y_head = top - 0.45,
                y_row = top - 0.45 - j)
YL12   <- c(-17.70, 1.75)
SEP12  <- c(-5.05, -10.65)
YCUE12 <- -17.25
OFF_L1 <- c(ox_without = 0.22, ox_with = 0, pr = -0.22)
OFF_L2 <- c(ox_without = 0.30, ox_with = 0.15, pr = 0, myc_without = -0.15,
            myc_with = -0.30)
P12 <- PRIM %>% dplyr::select(-block) %>%
  dplyr::mutate(block = paste(cohort, instrument))
PF_L1 <- .long(P12, OFF_L1, ROWS12 %>% dplyr::select(readout, block, y_row),
               c("readout", "block"))
PF_L2 <- .long(P12, OFF_L2, ROWS12 %>% dplyr::select(readout, block, y_row),
               c("readout", "block"))
CN_L1 <- .connectors(PF_L1, c("readout", "block"))
CN_L2 <- .connectors(PF_L2, c("readout", "block"))
.assert(paste0("layouts 1-2: connectors join only the two rungs of one term in ",
               "one row"),
        nrow(CN_L1) == 12L && nrow(CN_L2) == 24L &&
          all(CN_L1$term == "OX") && all(CN_L2$term %in% c("OX", "MYC")),
        n = nrow(CN_L1) + nrow(CN_L2))
XL12 <- .sym_lim(c(PF_L2$lo, PF_L2$hi))
.check_axis(PF_L1, XL12, "layout 1")
.check_axis(PF_L2, XL12, "layout 2")

RHO_OP <- G3 %>% dplyr::filter(estimator == EST_READ) %>%
  dplyr::transmute(block = paste(cohort, instrument), rho = rho_OX_PROLIF)
RHO_MP <- G3 %>% dplyr::filter(estimator == EST_READ) %>%
  dplyr::distinct(cohort, rho = rho_MYC_PROLIF)
.assert("the MYC-proliferation rho does not depend on the OXPHOS measure",
        nrow(RHO_MP) == 2L, n = nrow(RHO_MP))
R12 <- ROWS12 %>%
  dplyr::left_join(RHO_OP, by = "block") %>%
  dplyr::mutate(cohort = sub(" .*$", "", block)) %>%
  dplyr::left_join(RHO_MP %>% dplyr::rename(rho_mp = rho), by = "cohort") %>%
  dplyr::left_join(P12 %>% dplyr::select(readout, block, diff, diff_lo, diff_hi),
                   by = c("readout", "block"))

.lab12 <- function() {
  heads <- ROWS12 %>% dplyr::distinct(readout, y_head)
  .text_panel(YL12, list(
    .lyr(0.02, heads$y_head, unname(GENE_LAB[heads$readout]), hjust = 0,
         size = TXT_MM, parse = TRUE),
    .lyr(0.12, ROWS12$y_row, unname(BLOCK_LAB[ROWS12$block]), hjust = 0,
         size = SUB_MM, parse = TRUE)))
}
.rho_col <- function(head, vals) {
  .text_panel(YL12, list(
    .lyr(0.5, YH1, head, size = SUB_MM, parse = TRUE, face = "bold"),
    .lyr(0.5, YH2, "proliferation", size = SUB_MM, face = "bold"),
    .lyr(0.5, R12$y_row, sprintf("%.2f", vals), size = TXT_MM)))
}
.diff_col12 <- function() {
  .text_panel(YL12, list(
    .lyr(0.5, YH1, HEAD_DIFF, size = SUB_MM, parse = TRUE, face = "bold"),
    .lyr(0.5, YH2, "(95% CI)", size = SUB_MM, face = "bold"),
    .lyr(0.5, R12$y_row, .fmt2(R12$diff, R12$diff_lo, R12$diff_hi),
         size = TXT_MM)))
}
FOR_L1 <- .forest(PF_L1, CN_L1, XL12, YL12, SEP12, YCUE12, X_TITLE_SD,
                  names(OFF_L1), breaks = .sym_breaks(XL12[2]), size = 1.6)
FOR_L2 <- .forest(PF_L2, CN_L2, XL12, YL12, SEP12, YCUE12, X_TITLE_SD,
                  names(OFF_L2), breaks = .sym_breaks(XL12[2]), size = 1.4,
                  legend_ncol = 3)
PARTS_L1 <- list(.lab12(), .rho_col(HEAD_RHO_OP, R12$rho), FOR_L1,
                 .diff_col12())
PARTS_L2 <- list(.lab12(), .rho_col(HEAD_RHO_OP, R12$rho),
                 .rho_col(HEAD_RHO_MP, R12$rho_mp), FOR_L2, .diff_col12())
W_L1 <- c(1.05, 0.52, 4.20, 1.25)
W_L2 <- c(1.05, 0.52, 0.52, 3.68, 1.25)
FIG_L1 <- .compose(PARTS_L1, W_L1)
FIG_L2 <- .compose(PARTS_L2, W_L2)
H_L12 <- 3.70
CK_L1 <- .check_points(FOR_L1, PF_L1, "layout 1")
CK_L2 <- .check_points(FOR_L2, PF_L2, "layout 2")
.assert("layout 2 draws all three validated colours",
        setequal(setdiff(CK_L2$colours, "black"), VALIDATED_PALETTE), n = 3L)
.check_parts(PARTS_L1, "layout 1", W_PAGE)
.check_parts(PARTS_L2, "layout 2", W_PAGE)

# =============================================================================
# 8b. Fig. 4g, layout 3: the whole argument in one panel (12.2)
# =============================================================================
ROWS3 <- tidyr::expand_grid(readout = READ3, cohort = COHORTS) %>%
  dplyr::mutate(r = dplyr::row_number(), i = match(readout, READ3),
                y_row = -r - (i - 1) * 0.4)
YL3A <- c(-8.65, 1.75)
SEP3A <- c(-2.7, -5.1)
YCUE3A <- -8.25
P3 <- PRIM %>% dplyr::filter(instrument == "gsva") %>% dplyr::select(-block)
PF_3A <- .long(P3, OFF_L1, ROWS3 %>% dplyr::select(readout, cohort, y_row),
               c("readout", "cohort"))
CN_3A <- .connectors(PF_3A, c("readout", "cohort"))
XL3A <- .sym_lim(c(PF_3A$lo, PF_3A$hi))
.check_axis(PF_3A, XL3A, "layout 3, BCL-2 family")

# The clinical facets read E45's plotting frame and change nothing.
PF45 <- e45$plot_frame %>% dplyr::mutate(key = as.character(key))
MAP45 <- c("OXPHOS, m0" = "ox_without", "OXPHOS, m1" = "ox_with",
           "Proliferation, m1" = "pr")
PCR_ROWS <- tibble::tibble(row_key = c("GSE194040", "GSE164458", "GSE25066",
                                       "Pooled"),
                           y_row = c(-1, -2, -3, -4.4))
MB_ROW <- "METABRIC / ER-pos (IHC)"
OFF_CL <- c(ox_without = 0.22, ox_with = 0, pr = -0.22)
OFF_MB <- c(ox_with = 0.11, pr = -0.11)
PF_3B <- PF45 %>% dplyr::filter(panel == "pCR") %>%
  dplyr::transmute(row_key, key = unname(MAP45[key]), term, rung,
                   est = ratio, lo = ratio_lo, hi = ratio_hi, void = FALSE) %>%
  dplyr::left_join(PCR_ROWS, by = "row_key") %>%
  dplyr::mutate(y_mark = y_row + unname(OFF_CL[key]),
                key = factor(key, levels = names(OFF_CL)))
PF_3C <- PF45 %>% dplyr::filter(panel == "METABRIC", row_key == MB_ROW) %>%
  dplyr::transmute(row_key, key = unname(MAP45[key]), term, rung,
                   est = ratio, lo = ratio_lo, hi = ratio_hi, void = FALSE) %>%
  dplyr::mutate(y_row = -1, y_mark = y_row + unname(OFF_MB[key]),
                key = factor(key, levels = names(OFF_MB)))
CHK45 <- dplyr::bind_rows(PF_3B, PF_3C) %>%
  dplyr::mutate(key45 = names(MAP45)[match(as.character(key), MAP45)]) %>%
  dplyr::left_join(PF45 %>% dplyr::select(row_key, key45 = key, ratio,
                                          ratio_lo, ratio_hi),
                   by = c("row_key", "key45"))
.assert("layout 3: the clinical values equal E45's plotting frame exactly",
        nrow(PF_3B) == 12L && nrow(PF_3C) == 2L && !anyNA(CHK45$ratio) &&
          identical(CHK45$est, CHK45$ratio) &&
          identical(CHK45$lo, CHK45$ratio_lo) &&
          identical(CHK45$hi, CHK45$ratio_hi),
        n = 3L * nrow(CHK45))
CN_3B <- .connectors(PF_3B, "row_key")
XL3B <- .log_lim(c(PF_3B$lo, PF_3B$hi))
XL3C <- .log_lim(c(PF_3C$lo, PF_3C$hi))
.check_axis(PF_3B, XL3B, "layout 3, pCR", log_x = TRUE)
.check_axis(PF_3C, XL3C, "layout 3, METABRIC", log_x = TRUE)
YL3B <- c(-6.20, 1.75)
YCUE3B <- -5.75
YL3C <- c(-2.80, 1.75)
YCUE3C <- -2.35

RHO45 <- e45$exposure_rho
CON45 <- e45$contrast_column
rho_cl <- function(k) {
  r <- RHO45$rho[RHO45$cohort == k]
  if (length(r) == 1L) sprintf("%.2f", r) else ""
}
LBL_3A <- ROWS3 %>%
  dplyr::left_join(RHO_OP %>% dplyr::filter(grepl("gsva$", block)) %>%
                     dplyr::mutate(cohort = sub(" .*$", "", block)) %>%
                     dplyr::select(cohort, rho), by = "cohort") %>%
  dplyr::left_join(P3 %>% dplyr::select(readout, cohort, diff, diff_lo,
                                        diff_hi), by = c("readout", "cohort"))
.facet_parts <- function(title, rows_y, row_lab, rho_txt, right_txt, right_h1,
                         right_h2, forest, ylim, parse_right_h1 = FALSE) {
  list(
    .text_panel(ylim, list(
      .lyr(0.02, (YH1 + YH2) / 2, title, hjust = 0, size = TXT_MM,
           face = "bold"),
      .lyr(0.04, rows_y, row_lab, hjust = 0, size = SUB_MM, parse = TRUE))),
    .text_panel(ylim, list(
      .lyr(0.5, YH1, HEAD_RHO_OP, size = SUB_MM, parse = TRUE, face = "bold"),
      .lyr(0.5, YH2, "proliferation", size = SUB_MM, face = "bold"),
      .lyr(0.5, rows_y, rho_txt, size = TXT_MM))),
    forest,
    .text_panel(ylim, list(
      .lyr(0.5, YH1, right_h1, size = SUB_MM, parse = parse_right_h1,
           face = "bold"),
      .lyr(0.5, YH2, right_h2, size = SUB_MM, face = "bold"),
      .lyr(0.5, rows_y, right_txt, size = TXT_MM))))
}
FOR_3A <- .forest(PF_3A, CN_3A, XL3A, YL3A, SEP3A, YCUE3A, X_TITLE_SD,
                  names(OFF_L1), breaks = .sym_breaks(XL3A[2]), size = 1.6)
FOR_3B <- .forest(PF_3B, CN_3B, XL3B, YL3B, -3.7, YCUE3B,
                  "Odds ratio for pCR per within-cohort SD", names(OFF_CL),
                  cue = c("lower pCR", "higher pCR"), legend = FALSE,
                  log_x = TRUE,
                  breaks = BRK_OR[BRK_OR >= XL3B[1] & BRK_OR <= XL3B[2]],
                  size = 1.6)
FOR_3C <- .forest(PF_3C, NULL, XL3C, YL3C, numeric(0), YCUE3C,
                  "Hazard ratio for breast-cancer death per within-cohort SD",
                  names(OFF_MB), cue = c("longer survival", "shorter survival"),
                  legend = FALSE, log_x = TRUE,
                  breaks = BRK_HR[BRK_HR >= XL3C[1] & BRK_HR <= XL3C[2]],
                  size = 1.6)
CONTR_PCR <- CON45$label[match(PCR_ROWS$row_key, CON45$row_key)]
CONTR_MB  <- CON45$label[CON45$row_key == MB_ROW]
.assert("layout 3: the outcome contrasts are E45's labels, unchanged",
        !anyNA(CONTR_PCR) && length(CONTR_MB) == 1L, n = 5L)
PARTS_3A <- .facet_parts(
  "BCL-2 family,\nOXPHOS level", LBL_3A$y_row,
  paste0('italic("', LBL_3A$readout, '") %.% "', LBL_3A$cohort, '"'),
  sprintf("%.2f", LBL_3A$rho),
  .fmt2(LBL_3A$diff, LBL_3A$diff_lo, LBL_3A$diff_hi),
  HEAD_DIFF, "(95% CI)", FOR_3A, YL3A, parse_right_h1 = TRUE)
PARTS_3B <- .facet_parts(
  "Chemotherapy response", PCR_ROWS$y_row,
  paste0('"', PCR_ROWS$row_key, '"'),
  vapply(PCR_ROWS$row_key, rho_cl, character(1)),
  CONTR_PCR, "Proliferation / OXPHOS", "ratio of ORs (95% CI)", FOR_3B, YL3B)
PARTS_3C <- .facet_parts(
  "Breast-cancer-specific\nsurvival", -1, '"METABRIC, ER+ (IHC)"',
  rho_cl("METABRIC"), CONTR_MB, "Proliferation / OXPHOS",
  "ratio of HRs (95% CI)", FOR_3C, YL3C)
PARTS_L3 <- c(PARTS_3A, PARTS_3B, PARTS_3C)
W_L3 <- c(1.25, 0.52, 3.98, 1.25)
FIG_L3 <- .compose(PARTS_L3, W_L3, ncol = 4,
                   heights = c(diff(YL3A), diff(YL3B), diff(YL3C)))
H_L3 <- 4.90
.check_points(FOR_3A, PF_3A, "layout 3, BCL-2 family")
.check_points(FOR_3B, PF_3B, "layout 3, pCR", log_x = TRUE)
.check_points(FOR_3C, PF_3C, "layout 3, METABRIC", log_x = TRUE)
.check_parts(PARTS_L3, "layout 3", W_PAGE)

# =============================================================================
# 8c. ED-a: correlations among the scores, with E39's three values (12.3)
# =============================================================================
BANDS <- c("below 0.60", "0.60 to 0.80", "0.80 or above")
BAND_FILL <- c("below 0.60" = "#F0F0F0", "0.60 to 0.80" = "#BDBDBD",
               "0.80 or above" = "#636363")
BAND_INK  <- c("below 0.60" = "black", "0.60 to 0.80" = "black",
               "0.80 or above" = "white")
TILE <- CORR %>%
  dplyr::mutate(x = factor(unname(VARS[var_col]), levels = unname(VARS)),
                y = factor(unname(VARS[var_row]), levels = rev(unname(VARS))),
                label = sprintf("%.2f", rho),
                band = factor(band, levels = BANDS))
TILE39 <- E39G %>%
  dplyr::mutate(x = factor("Reference MYC score\nvs proliferation"),
                y = factor(cohort, levels = rev(c("GSE25066", "GSE194040",
                                                  "GSE164458"))),
                label = sprintf("%.2f", rho), band = factor(band, levels = BANDS))
.tiles <- function(d, title, show_y = TRUE, legend = FALSE) {
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = band), colour = "white",
                       linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = band),
                       size = SUB_MM) +
    ggplot2::scale_fill_manual(values = BAND_FILL, breaks = BANDS, drop = FALSE,
                               name = "Spearman rho, absolute") +
    ggplot2::scale_colour_manual(values = BAND_INK, guide = "none") +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = BASE_PT) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = MIN_PT, angle = 40,
                                          hjust = 1, colour = "black"),
      axis.text.y = if (show_y) ggplot2::element_text(size = MIN_PT,
                                                      colour = "black") else
        ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = BASE_PT, face = "bold",
                                         hjust = 0.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.title = ggplot2::element_text(size = BASE_PT - 0.5),
      strip.text = ggplot2::element_text(size = BASE_PT),
      plot.caption = ggplot2::element_text(size = MIN_PT),
      plot.margin = ggplot2::margin(2, 3, 2, 3))
  if (!legend) p <- p + ggplot2::guides(fill = "none")
  p
}
PARTS_EDA <- list(
  .tiles(TILE %>% dplyr::filter(cohort == "TCGA"), "TCGA-BRCA (n = 938)",
         legend = TRUE),
  .tiles(TILE %>% dplyr::filter(cohort == "SCAN-B"), "SCAN-B (n = 3,143)",
         show_y = FALSE),
  .tiles(TILE39, "Neoadjuvant cohorts"))
W_EDA <- c(2.75, 1.85, 1.30)
FIG_EDA <- .compose(PARTS_EDA, W_EDA)
H_EDA <- 2.90
for (k in seq_along(PARTS_EDA)) {
  b <- ggplot2::ggplot_build(PARTS_EDA[[k]])
  dd <- list(TILE %>% dplyr::filter(cohort == "TCGA"),
             TILE %>% dplyr::filter(cohort == "SCAN-B"), TILE39)[[k]]
  .assert(paste0("ED-a panel ", k, ": every tile drawn carries its stored rho"),
          nrow(b$data[[1]]) == nrow(dd) &&
            identical(sort(as.character(b$data[[2]]$label)), sort(dd$label)),
          n = nrow(dd))
}
.check_parts(PARTS_EDA, "ED-a", W_PAGE)

# =============================================================================
# 8d. ED-b: E47 Display C relabelled, with C3's m3 result beside it (12.3)
# =============================================================================
DC <- e47$displays$C_coefficients
DD <- e47$displays$C_difference
S9 <- e47$bcl2l1_bbc3
chk_c <- max(
  .dev(DC$est[DC$gene == "BCL2L1"], S9$b_OX_BCL2L1),
  .dev(DC$est[DC$gene == "BBC3"], S9$b_OX_BBC3),
  .dev(DD$est, S9$delta), .dev(DD$lo, S9$ci_lo), .dev(DD$hi, S9$ci_hi))
.assert("ED-b: E47's display frames hold E47's section 9 values",
        chk_c == 0, n = 20L, max_dev = chk_c)
EDB_ROWS <- tibble::tibble(block = BLOCK_ORDER, y_row = -1.2 * (1:4))
OFF_EDB  <- c(bcl2l1_without = 0.30, bcl2l1_with = 0.12, bbc3_without = -0.12,
              bbc3_with = -0.30)
OFF_EDBD <- c(diff_without = 0.12, diff_with = -0.12)
C3M3 <- C3 %>% dplyr::filter(rung == "m3")
PF_EDB <- dplyr::bind_rows(
  DC %>% dplyr::transmute(block = paste(cohort, instrument),
                          key = ifelse(gene == "BCL2L1", "bcl2l1_without",
                                       "bbc3_without"),
                          est, lo, hi, source = "E47 Display C"),
  C3M3 %>% dplyr::transmute(block = paste(cohort, instrument),
                            key = "bcl2l1_with", est = b_OX_BCL2L1,
                            lo = lo_BCL2L1, hi = hi_BCL2L1,
                            source = "E47 m3 fit"),
  C3M3 %>% dplyr::transmute(block = paste(cohort, instrument),
                            key = "bbc3_with", est = b_OX_BBC3, lo = lo_BBC3,
                            hi = hi_BBC3, source = "E47 m3 fit")) %>%
  dplyr::left_join(EDB_ROWS, by = "block") %>%
  dplyr::mutate(y_mark = y_row + unname(OFF_EDB[key]), void = FALSE,
                key = factor(key, levels = names(OFF_EDB)))
PF_EDBD <- dplyr::bind_rows(
  DD %>% dplyr::transmute(block = paste(cohort, instrument),
                          key = "diff_without", est, lo, hi,
                          source = "E47 Display C"),
  C3M3 %>% dplyr::transmute(block = paste(cohort, instrument),
                            key = "diff_with", est = delta, lo = ci_lo,
                            hi = ci_hi, source = "C3")) %>%
  dplyr::left_join(EDB_ROWS, by = "block") %>%
  dplyr::mutate(y_mark = y_row + unname(OFF_EDBD[key]), void = FALSE,
                key = factor(key, levels = names(OFF_EDBD)))
chk_b <- PF_EDB %>% dplyr::filter(grepl("without", key)) %>%
  dplyr::mutate(gene = ifelse(key == "bcl2l1_without", "BCL2L1", "BBC3"),
                cohort = sub(" .*$", "", block), instrument = sub("^.* ", "", block)) %>%
  dplyr::left_join(DC %>% dplyr::select(cohort, instrument, gene, est0 = est,
                                        lo0 = lo, hi0 = hi),
                   by = c("cohort", "instrument", "gene"))
.assert("ED-b: the relabelled Display C values are identical to E47's display frame",
        identical(chk_b$est, chk_b$est0) && identical(chk_b$lo, chk_b$lo0) &&
          identical(chk_b$hi, chk_b$hi0), n = 3L * nrow(chk_b))
CN_EDB  <- .connectors(PF_EDB %>% dplyr::mutate(row = block), "row")
CN_EDBD <- .connectors(PF_EDBD %>% dplyr::mutate(row = block), "row")
YL_EDB <- c(-6.70, -0.35)
YCUE_EDB <- -6.30
XL_EDB  <- .sym_lim(c(PF_EDB$lo, PF_EDB$hi))
XL_EDBD <- .sym_lim(c(PF_EDBD$lo, PF_EDBD$hi))
.check_axis(PF_EDB, XL_EDB, "ED-b, coefficients")
.check_axis(PF_EDBD, XL_EDBD, "ED-b, difference")
FOR_EDB <- .forest(PF_EDB, CN_EDB, XL_EDB, YL_EDB, numeric(0), YCUE_EDB,
                   "OXPHOS association, MYC activity held fixed (SD per SD)",
                   names(OFF_EDB), breaks = .sym_breaks(XL_EDB[2]), size = 1.5,
                   legend_ncol = 2)
FOR_EDBD <- .forest(PF_EDBD, CN_EDBD, XL_EDBD, YL_EDB, numeric(0), YCUE_EDB,
                    "BCL2L1 minus BBC3 (bootstrap 95% CI)", names(OFF_EDBD),
                    cue = c("BBC3 steeper", "BCL2L1 steeper"),
                    breaks = .sym_breaks(XL_EDBD[2]), size = 1.5,
                    legend_ncol = 1)
PARTS_EDB <- list(
  .text_panel(YL_EDB, list(.lyr(0.04, EDB_ROWS$y_row,
                                unname(BLOCK_LAB[EDB_ROWS$block]), hjust = 0,
                                size = SUB_MM, parse = TRUE))),
  FOR_EDB, FOR_EDBD)
W_EDB <- c(0.95, 3.45, 2.60)
FIG_EDB <- .compose(PARTS_EDB, W_EDB)
H_EDB <- 2.45
.check_points(FOR_EDB, PF_EDB, "ED-b, coefficients")
.check_points(FOR_EDBD, PF_EDBD, "ED-b, difference")
.check_parts(PARTS_EDB, "ED-b", W_PAGE)

# =============================================================================
# 8e. ED-c: MYC activity on all four scores, with proliferation beside (12.3)
# =============================================================================
ROWS_C <- tidyr::expand_grid(estimator = EST_ALL, block = BLOCK_ORDER) %>%
  dplyr::mutate(s = match(estimator, EST_ALL), j = match(block, BLOCK_ORDER),
                top = -(s - 1) * 5.0, y_head = top - 0.30,
                y_row = top - 0.45 - j)
YL_C <- c(-21.50, 0.30)
SEP_C <- -(1:3) * 5.0 + 0.25
YCUE_C <- -21.05
OFF_C <- c(myc_without = 0.22, myc_with = 0, pr = -0.22)
PANEL_LAB <- unname(GENE_LAB[READ3])
E_C <- EXT %>% dplyr::filter(readout %in% READ3) %>%
  dplyr::mutate(block = paste(cohort, instrument))
PF_C <- .long(E_C, OFF_C, ROWS_C %>% dplyr::select(estimator, block, y_row),
              c("estimator", "block"), keep = "readout") %>%
  dplyr::mutate(panel = factor(unname(GENE_LAB[readout]), levels = PANEL_LAB))
.stop_if(nrow(PF_C) == 3L * 16L * 3L, "ED-c frame has the wrong size")
CN_C <- .connectors(PF_C %>% dplyr::mutate(row = paste(readout, estimator,
                                                       block)), "row")
XL_C <- .sym_lim(c(PF_C$lo, PF_C$hi))
.check_axis(PF_C, XL_C, "ED-c")
FOR_C <- .forest(PF_C, CN_C, XL_C, YL_C, SEP_C, YCUE_C, X_TITLE_SD,
                 names(OFF_C), breaks = .sym_breaks(XL_C[2]), size = 1.1,
                 facet = TRUE)
HEADS_C <- ROWS_C %>% dplyr::distinct(estimator, y_head)
LBL_C <- .text_panel(YL_C, list(
  .lyr(0.02, HEADS_C$y_head, unname(EST_NAME[HEADS_C$estimator]), hjust = 0,
       size = SUB_MM, face = "bold"),
  .lyr(0.02, HEADS_C$y_head[HEADS_C$estimator == "entangled"] - 0.62,
       "with proliferation: not identifiable, faded", hjust = 0, size = SUB_MM,
       colour = "grey35"),
  .lyr(0.10, ROWS_C$y_row, unname(BLOCK_LAB[ROWS_C$block]), hjust = 0,
       size = SUB_MM, parse = TRUE)))
PARTS_EDC <- list(LBL_C, FOR_C)
W_EDC <- c(1.55, 5.45)
FIG_EDC <- .compose(PARTS_EDC, W_EDC)
H_EDC <- 4.40
.check_points(FOR_C, PF_C, "ED-c")
.check_parts(PARTS_EDC, "ED-c", W_PAGE)

# =============================================================================
# 8f. ED-d: all twelve genes, three coefficients, both scores (12.3)
# =============================================================================
ROWS_D <- tidyr::expand_grid(readout = TWELVE, block = BLOCK_ORDER) %>%
  dplyr::mutate(g = match(readout, TWELVE), j = match(block, BLOCK_ORDER),
                top = -0.8 - (g - 1) * 4.6 - ifelse(readout %in% ANTI, 1.2, 0),
                y_row = top - j)
YL_D <- c(-58.70, 0.30)
YCUE_D <- -58.25
SEP_D <- ROWS_D %>% dplyr::distinct(g, top) %>%
  dplyr::filter(!g %in% c(1L, 8L)) %>% dplyr::pull(top) + 0.3
OFF_D <- c(ox_with = 0.25, pr = 0, myc_with = -0.25)
.edd <- function(est) {
  E <- EXT %>% dplyr::filter(estimator == est) %>%
    dplyr::mutate(block = paste(cohort, instrument))
  pf <- .long(E, OFF_D, ROWS_D %>% dplyr::select(readout, block, y_row),
              c("readout", "block"))
  pd <- E %>% dplyr::transmute(readout, block, key = "diff", term = "difference",
                               rung = "m3", est = diff, lo = diff_lo,
                               hi = diff_hi, void = FALSE) %>%
    dplyr::left_join(ROWS_D %>% dplyr::select(readout, block, y_row),
                     by = c("readout", "block")) %>%
    dplyr::mutate(y_mark = y_row, key = factor(key, levels = "diff"))
  list(coef = pf, diff = pd)
}
D_REF <- .edd(EST_READ)
D_CLN <- .edd(EST_AGREE)
.assert("ED-d: 12 genes x 4 blocks x 3 coefficients on each score",
        nrow(D_REF$coef) == 144L && nrow(D_CLN$coef) == 144L &&
          nrow(D_REF$diff) == 48L && nrow(D_CLN$diff) == 48L, n = 384L)
XL_D  <- .sym_lim(c(D_REF$coef$lo, D_REF$coef$hi, D_CLN$coef$lo,
                    D_CLN$coef$hi))
XL_DD <- .sym_lim(c(D_REF$diff$lo, D_REF$diff$hi, D_CLN$diff$lo,
                    D_CLN$diff$hi))
for (nm in c("ref", "clean")) {
  dd <- if (nm == "ref") D_REF else D_CLN
  .check_axis(dd$coef, XL_D, paste0("ED-d, ", nm, ", coefficients"))
  .check_axis(dd$diff, XL_DD, paste0("ED-d, ", nm, ", difference"))
}
.ford <- function(pf, title, xlim, keys, x_title, diff_panel = FALSE) {
  .forest(pf, NULL, xlim, YL_D, SEP_D, YCUE_D, x_title, keys,
          cue = if (diff_panel) NULL else c("lower expression",
                                            "higher expression"),
          breaks = if (diff_panel) .coarse_breaks(xlim[2]) else
            .sym_breaks(xlim[2]), size = 0.9, title = title)
}
FOR_D1 <- .ford(D_REF$coef, "Felsher, reference score", XL_D, names(OFF_D),
                "SD per SD of each score")
FOR_D2 <- .ford(D_REF$diff, "Difference", XL_DD, "diff",
                "Proliferation\nminus OXPHOS", diff_panel = TRUE)
FOR_D3 <- .ford(D_CLN$coef, "MYC_UP.V1_UP, least overlap", XL_D, names(OFF_D),
                "SD per SD of each score")
FOR_D4 <- .ford(D_CLN$diff, "Difference", XL_DD, "diff",
                "Proliferation\nminus OXPHOS", diff_panel = TRUE)
GENE_POS <- ROWS_D %>% dplyr::group_by(readout) %>%
  dplyr::summarise(y = mean(y_row), .groups = "drop")
LBL_D <- .text_panel(YL_D, list(
  .lyr(0.02, c(-0.35, ROWS_D$top[ROWS_D$readout == "BCL2L1"][1] + 0.55),
       c("Pro-apoptotic", "Anti-apoptotic"), hjust = 0, size = SUB_MM,
       face = "bold"),
  .lyr(0.02, GENE_POS$y,
       ifelse(GENE_POS$readout %in% READ3,
              paste0('bolditalic("', GENE_POS$readout, '")'),
              paste0('italic("', GENE_POS$readout, '")')),
       hjust = 0, size = TXT_MM, parse = TRUE),
  .lyr(0.40, ROWS_D$y_row, unname(BLOCK_LAB[ROWS_D$block]), hjust = 0,
       size = SUB_MM, parse = TRUE)))
PARTS_EDD <- list(LBL_D, FOR_D1, FOR_D2, FOR_D3, FOR_D4)
W_EDD <- c(1.35, 1.90, 0.95, 1.90, 0.95)
FIG_EDD <- .compose(PARTS_EDD, W_EDD)
H_EDD <- 7.80
.check_points(FOR_D1, D_REF$coef, "ED-d, reference, coefficients")
.check_points(FOR_D2, D_REF$diff, "ED-d, reference, difference")
.check_points(FOR_D3, D_CLN$coef, "ED-d, least overlap, coefficients")
.check_points(FOR_D4, D_CLN$diff, "ED-d, least overlap, difference")
.check_parts(PARTS_EDD, "ED-d", W_PAGE)

# =============================================================================
# 8g. S-1: C1 and C4 beside the primary, TCGA (12.3)
# =============================================================================
S1_SPECS <- tibble::tribble(
  ~estimator, ~spec,                                    ~r, ~label,
  "M_a",      "as fitted",                              1L, "Felsher, as fitted",
  "M_a",      "with TP53 status",                       2L, "Felsher, with TP53 status",
  "M_a",      "without purity and leukocyte fraction",  3L, "Felsher, without purity and leukocytes",
  "clean",    "as fitted",                              4L, "MYC_UP.V1_UP, as fitted",
  "clean",    "without purity and leukocyte fraction",  5L, "MYC_UP.V1_UP, without purity and leukocytes")
S1_GROUPS <- tidyr::expand_grid(readout = READ3, instrument = names(INSTR)) %>%
  dplyr::mutate(grp = dplyr::row_number(), top = -(grp - 1) * 6.3,
                y_head = top - 0.45)
S1_ROWS <- tidyr::expand_grid(S1_GROUPS, S1_SPECS) %>%
  dplyr::mutate(y_row = top - 0.45 - r)
YL_S1 <- c(-39.00, 0.30)
YCUE_S1 <- -38.55
SEP_S1 <- -(1:5) * 6.3 + 0.35
SENS_Y <- SENS %>% dplyr::left_join(
  S1_ROWS %>% dplyr::select(readout, instrument, estimator, spec, y_row),
  by = c("readout", "instrument", "estimator", "spec"))
.assert("S-1: every sensitivity row has its place", !anyNA(SENS_Y$y_row) &&
          nrow(SENS_Y) == 30L, n = nrow(SENS_Y))
.s1pf <- function(col_est, col_lo, col_hi, key) {
  SENS_Y %>% dplyr::transmute(readout, instrument, estimator, spec, key = key,
                              est = .data[[col_est]], lo = .data[[col_lo]],
                              hi = .data[[col_hi]], y_row, y_mark = y_row,
                              void = FALSE) %>%
    dplyr::mutate(key = factor(key, levels = key[1]))
}
PF_S1 <- list(ox = .s1pf("b_OX", "ox_lo", "ox_hi", "ox_with"),
              pr = .s1pf("b_PR", "pr_lo", "pr_hi", "pr"),
              diff = .s1pf("diff", "diff_lo", "diff_hi", "diff"))
XL_S1 <- lapply(PF_S1, function(pf) .sym_lim(c(pf$lo, pf$hi)))
for (nm in names(PF_S1)) .check_axis(PF_S1[[nm]], XL_S1[[nm]], paste0("S-1, ", nm))
FOR_S1 <- list(
  ox = .forest(PF_S1$ox, NULL, XL_S1$ox, YL_S1, SEP_S1, YCUE_S1,
               "OXPHOS (SD per SD)", "ox_with",
               breaks = .sym_breaks(XL_S1$ox[2]), size = 1.1,
               title = "OXPHOS"),
  pr = .forest(PF_S1$pr, NULL, XL_S1$pr, YL_S1, SEP_S1, YCUE_S1,
               "Proliferation (SD per SD)", "pr",
               breaks = .sym_breaks(XL_S1$pr[2]), size = 1.1,
               title = "Proliferation"),
  diff = .forest(PF_S1$diff, NULL, XL_S1$diff, YL_S1, SEP_S1, YCUE_S1,
                 "Proliferation minus OXPHOS", "diff", cue = NULL,
                 breaks = .sym_breaks(XL_S1$diff[2]), size = 1.1,
                 title = "Difference"))
LBL_S1 <- .text_panel(YL_S1, list(
  .lyr(0.02, S1_GROUPS$y_head,
       paste0('italic("', S1_GROUPS$readout, '") %.% "TCGA" %.% "',
              unname(INSTR[S1_GROUPS$instrument]), '"'),
       hjust = 0, size = TXT_MM, parse = TRUE),
  .lyr(0.06, S1_ROWS$y_row, S1_ROWS$label, hjust = 0, size = SUB_MM)))
PARTS_S1 <- list(LBL_S1, FOR_S1$ox, FOR_S1$pr, FOR_S1$diff)
W_S1 <- c(1.95, 1.70, 1.70, 1.65)
FIG_S1 <- .compose(PARTS_S1, W_S1)
H_S1 <- 5.50
for (nm in names(PF_S1)) .check_points(FOR_S1[[nm]], PF_S1[[nm]],
                                       paste0("S-1, ", nm))
.check_parts(PARTS_S1, "S-1", W_PAGE)

# --- save every figure: outputs/figures, and docs/figures (tracked) ----------
.ensure_dir(DIR_FIGURES)
.ensure_dir(DIR_DOC_FIGURES)
.save_fig <- function(fig, stem, h) {
  for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
    base <- file.path(dir, if (identical(dir, DIR_DOC_FIGURES))
      paste0(FIG_DATE, "_", stem) else stem)
    ggplot2::ggsave(paste0(base, ".png"), fig, width = W_PAGE, height = h,
                    units = "in", dpi = 300, bg = "white")
    ggplot2::ggsave(paste0(base, ".pdf"), fig, width = W_PAGE, height = h,
                    units = "in", device = grDevices::cairo_pdf, bg = "white")
  }
}
FIGS <- list(layout1 = list(FIG_L1, H_L12), layout2 = list(FIG_L2, H_L12),
             layout3 = list(FIG_L3, H_L3), ed_a = list(FIG_EDA, H_EDA),
             ed_b = list(FIG_EDB, H_EDB), ed_c = list(FIG_EDC, H_EDC),
             ed_d = list(FIG_EDD, H_EDD), s1 = list(FIG_S1, H_S1))
for (nm in names(FIGS)) .save_fig(FIGS[[nm]][[1]], FIG[[nm]], FIGS[[nm]][[2]])
message("   eight figures, 178 mm wide, to outputs/figures and docs/figures ",
        "(tracked)")

# =============================================================================
# 9. Legend drafts - DESCRIPTIVE ONLY (brief section 4)
# =============================================================================
N_TXT <- "TCGA-BRCA (n = 938) and SCAN-B (n = 3,143)"
MODEL_TXT <- paste0(
  "Coefficients (SD per SD of each score) are from linear models of each gene ",
  "on OXPHOS, MYC activity and intrinsic subtype, with tumour purity and ",
  "leukocyte fraction in TCGA-BRCA; variables are rank-normalised within each ",
  "cohort. OXPHOS was scored as subunit expression (level) or as its priority ",
  "within the mitochondrial programme (priority).")
LEGEND <- list(
  layout1 = paste0(
    "Associations of OXPHOS and proliferation with BBC3 (PUMA), BCL2L1 ",
    "(BCL-xL) and MCL1 expression in ", N_TXT, ". ", MODEL_TXT, " Open ",
    "circles, OXPHOS without the proliferation term; filled circles, OXPHOS ",
    "with it; triangles, proliferation. Left, Spearman correlation of OXPHOS ",
    "with proliferation. Right, the proliferation coefficient minus the OXPHOS ",
    "coefficient within the same model, with its 95% confidence interval from ",
    "the model covariance. Bars, 95% confidence intervals."),
  layout2 = paste0(
    "As layout 1, with MYC activity added: squares, MYC activity without ",
    "(open) and with (filled) the proliferation term. Second column, Spearman ",
    "correlation of MYC activity with proliferation."),
  layout3 = paste0(
    "Top, the BCL-2 family as in layout 1, OXPHOS level only. Middle, ",
    "pathological complete response in three neoadjuvant cohorts and their ",
    "pooled estimate: odds ratios per within-cohort SD from logistic models, ",
    "OXPHOS without (open) and with (filled) the proliferation term, and ",
    "proliferation (triangles); right, the ratio of the two odds ratios within ",
    "the same model. Bottom, breast-cancer-specific survival in ER-positive ",
    "(IHC) METABRIC: hazard ratios from Cox models with proliferation; right, ",
    "the ratio of hazard ratios. The BCL-2 family panel holds MYC activity ",
    "fixed; the outcome panels come from models without a MYC activity term. ",
    "The METABRIC correlation is over all 1,981 patients. Bars, 95% ",
    "confidence intervals."),
  ed_a = paste0(
    "Spearman correlations among the two OXPHOS scores, four MYC activity ",
    "scores and proliferation in ", N_TXT, ", shaded by band (below 0.60; ",
    "0.60 to 0.80; 0.80 or above). Each MYC activity score is labelled with ",
    "the percentage of its genes in proliferation gene sets (Hallmark E2F ",
    "targets and G2M checkpoint). Right, the reference MYC activity score ",
    "against proliferation in three neoadjuvant cohorts."),
  ed_b = paste0(
    "OXPHOS coefficients for BCL2L1 and BBC3 with MYC activity held fixed, ",
    "without (open) and with (filled) the proliferation term, in ", N_TXT,
    ", and the difference between them, with 95% percentile intervals from ",
    "2,000 paired bootstrap resamples of tumours. Bars on the coefficients, ",
    "model-based 95% confidence intervals."),
  ed_c = paste0(
    "MYC activity coefficients for BBC3, BCL2L1 and MCL1 on four MYC activity ",
    "scores, ordered by the percentage of their genes in proliferation gene ",
    "sets, without (open) and with (filled) the proliferation term, with the ",
    "proliferation coefficient (triangles), in ", N_TXT, ". Faded points are ",
    "cells in which the score and proliferation correlate at 0.80 or more; ",
    "they are shown and not interpreted. Bars, 95% confidence intervals."),
  ed_d = paste0(
    "All twelve BCL-2 family genes: coefficients of OXPHOS (circles), ",
    "proliferation (triangles) and MYC activity (squares), each with the ",
    "other two held fixed, and the proliferation minus OXPHOS difference, on ",
    "the reference MYC activity score (left) and the score with least overlap ",
    "with proliferation gene sets (right), in ", N_TXT, ". Bold, the three ",
    "genes read under the declared rule. No correction for multiple genes. ",
    "Bars, 95% confidence intervals."),
  s1 = paste0(
    "Sensitivity analyses in TCGA-BRCA (n = 938): OXPHOS and proliferation ",
    "coefficients and their difference for BBC3, BCL2L1 and MCL1, as fitted, ",
    "with TP53 mutation status added, and without tumour purity and leukocyte ",
    "fraction (the SCAN-B covariate set), on the reference MYC activity score ",
    "and the score with least overlap with proliferation gene sets. Bars, 95% ",
    "confidence intervals."))
BANNED <- c("independently of MYC", "MYC expression", "comparable", "invariant",
            "lost after adjustment", "signature", "primed", "priming")
.assert("no legend draft uses a word brief 4.9 keeps out",
        !any(vapply(BANNED, function(w) any(grepl(w, unlist(LEGEND),
                                                  ignore.case = TRUE)),
                    logical(1))), n = length(LEGEND))

# =============================================================================
# 10. Save, per R rules 4 and 5
# =============================================================================
message("\n10. saving")

ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("   ", nrow(ASSERTIONS), " assertions, all passed")

FRAME <- dplyr::bind_rows(DSET) %>%
  dplyr::mutate(PAM50 = as.character(PAM50),
                TP53_status = as.character(TP53_status)) %>%
  dplyr::select(sample_id, cohort, in_set, OX_gsva, OX_mitopps, MYC_clean,
                MYC_M_a, MYC_entangled, MYC_regulon, PROLIF, PAM50, purity,
                leuko, TP53_status, dplyr::all_of(TWELVE))

readr::write_csv(EXT,        file.path(DIR_TABLES, "E48_coefficients.csv"))
readr::write_csv(PRIM,       file.path(DIR_TABLES, "E48_primary_cells.csv"))
readr::write_csv(LABELS,     file.path(DIR_TABLES, "E48_labels.csv"))
readr::write_csv(AGREE,      file.path(DIR_TABLES, "E48_agreement.csv"))
readr::write_csv(NINE_TALLY, file.path(DIR_TABLES, "E48_nine_gene_tally.csv"))
readr::write_csv(NEST,       file.path(DIR_TABLES, "E48_nesting.csv"))
readr::write_csv(SENS,       file.path(DIR_TABLES, "E48_c1_c4_sensitivities.csv"))
readr::write_csv(C3,         file.path(DIR_TABLES, "E48_c3_more_steeply.csv"))
readr::write_csv(CORR,       file.path(DIR_TABLES, "E48_correlations.csv"))
readr::write_csv(ASSERTIONS, file.path(DIR_TABLES, "E48_assertions.csv"))

saveRDS(list(
  labels            = LABELS,
  primary_cells     = PRIM,
  agreement         = AGREE,
  text_row          = TEXT_ROW,
  nine_gene_tally   = NINE_TALLY,
  coefficients      = EXT,
  nesting           = NEST,
  gate_a2           = GA2,
  part_b            = PART_B,
  c1                = C1,
  c1_coefficient_matrix = dplyr::bind_rows(C1_COEF),
  c1_vcov           = C1_VCOV,
  c1_vif            = dplyr::bind_rows(C1_VIF),
  c1_gate           = C1_GATE,
  part_c2           = PART_C2,
  c3                = C3,
  c3_holds_m3       = C3_HOLDS_M3,
  c3_boot_replicates_m3 = BOOT_M3,
  c4                = C4,
  c4_coefficient_matrix = dplyr::bind_rows(C4_COEF),
  c4_vcov           = C4_VCOV,
  c4_vif            = dplyr::bind_rows(C4_VIF),
  c4_gate           = C4_GATE,
  sensitivities     = SENS,
  fit_digests       = dplyr::bind_rows(c(C1_DIG, C4_DIG)),
  correlations      = CORR,
  e39_gate          = E39G,
  frame             = FRAME,
  displays          = list(layout1 = PF_L1, layout1_connectors = CN_L1,
                           layout2 = PF_L2, layout2_connectors = CN_L2,
                           layout3_bcl2 = PF_3A, layout3_pcr = PF_3B,
                           layout3_metabric = PF_3C, rows_l12 = R12,
                           rows_l3 = LBL_3A, ed_a = TILE, ed_a_e39 = TILE39,
                           ed_b_coefficients = PF_EDB,
                           ed_b_difference = PF_EDBD, ed_c = PF_C,
                           ed_d_reference = D_REF, ed_d_least_overlap = D_CLN,
                           s1 = PF_S1,
                           axes = list(l12 = XL12, l3 = list(XL3A, XL3B, XL3C),
                                       ed_b = list(XL_EDB, XL_EDBD),
                                       ed_c = XL_C, ed_d = list(XL_D, XL_DD),
                                       s1 = XL_S1)),
  legend_draft      = LEGEND,
  assertions        = ASSERTIONS,
  # FILE NAMES, NOT PATHS, so this element digests the same in a redirected
  # dry run and in the author's run. Figures are not digested.
  figure_files      = list(
    outputs = paste0(rep(unname(FIG), each = 2), c(".png", ".pdf")),
    tracked_in_docs_figures = paste0(FIG_DATE, "_",
                                     rep(unname(FIG), each = 2),
                                     c(".png", ".pdf"))),
  spec = list(
    declaration = paste0(
      "docs/2026-11-02_E48_separation_declaration.md, committed alone at ",
      "2941438, before this script. Brief: ",
      "docs/2026-11-02_separation_line_brief.md, cbf2e4c."),
    posture = paste0("EXPLORATORY and POST HOC. Part A fits nothing. Part B ",
                     "unavailable. Part C reported, never read. Main effects ",
                     "only; no interaction."),
    reading_rule = paste0(
      "Section 7.1, transcribed in .label_gene(): over the 4 M_a m3 cells of ",
      "each of BBC3, BCL2L1, MCL1. OPPOSITE: in all 4 cells proliferation's ",
      "CI excludes zero with the sign opposite to OXPHOS's m3 coefficient in ",
      "that cell (undefined, so failed, if OXPHOS's CI includes zero), and ",
      "the difference's CI excludes zero. DISTINCT: not OPPOSITE, but the ",
      "difference's CI excludes zero in all 4 cells, with one sign. NOT ",
      "SEPARATED: anything else. Section 7.2: on MYC_UP.V1_UP at m3, ",
      "proliferation's sign must equal its sign on M_a in all 12 cells; any ",
      "disagreement prints ESTIMATOR-DEPENDENT beside that gene's label."),
    outcome_table = paste0(
      "Section 8. Rows 1-3 exclusive on the labels; row 4 in addition. ",
      "1, all three OPPOSITE, agreement holds: abstract 'opposite ",
      "associations'; brief section 8's sentence as drafted. 2, none NOT ",
      "SEPARATED and at least one DISTINCT: 'distinct associations'; name the ",
      "genes on which proliferation ran the other way. 3, any gene NOT ",
      "SEPARATED: that gene leaves the separation sentence (BCL2L1 takes the ",
      "reconciliation proposal's section 4 wording); the clinical paragraph ",
      "is unaffected; the abstract names the genes on which proliferation ran ",
      "the other way and does not say 'the BCL-2 family'. 4, any ",
      "ESTIMATOR-DEPENDENT: state the label with the dependence named, or not ",
      "at all. Always: any gene not OPPOSITE gets one descriptive sentence ",
      "giving the direction of its proliferation association in each cohort, ",
      "with no inference; the separation and closing sentences name the same ",
      "genes. No outcome licenses a causal, mediation, interaction or 'primed' ",
      "claim, or any link from the BCL-2 family to outcome."),
    nothing_discarded = paste0(
      "Section 8's closing paragraph: every coefficient and difference, for ",
      "all twelve genes, both MYC scores and C1, C3 and C4, is saved, ",
      "tabulated and displayed. The labels decide only the words used for ",
      "the three read genes."),
    palette = paste0("dataviz palette checker, 2026-11-03: #009E73, #D55E00, ",
                     "#0072B2, all pairs, light surface: ALL CHECKS PASS."),
    snapshot = paste0("E47's saved object (frame, coefficients, vcov); E45, ",
                      "E39, E40, E43 saved objects; data/from_validation/",
                      "tcga_brca_covariates.rds for TP53_status. Nothing ",
                      "rescored; joint snapshot not read."),
    not_digested = paste0("$fits_c1 and $fits_c4 hold fitted lm objects, ",
                          "which are not digest-stable; $fit_digests is the ",
                          "stable handle."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  fits_c1           = FITS_C1,
  fits_c4           = FITS_C4,
  built             = Sys.time()), PATH_E48)

message("   ", PATH_E48)
message("\nE48 done. Labels (section 7): ",
        paste(LABELS$readout, LABELS$printed, sep = " ", collapse = "; "),
        ". Outcome row ", TEXT_ROW$row,
        if (TEXT_ROW$row4_estimator_dependent) " plus row 4" else "", ".")
message("Part C is reported and never read. Part B is unavailable.\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E48)

  # the labels, the agreement check and the outcome row (sections 7 and 8)
  x$labels %>% as.data.frame()
  x$agreement %>% as.data.frame()
  str(x$text_row)

  # the twelve primary cells: M_a, m3
  x$primary_cells %>%
    dplyr::transmute(readout, block, b_OX_m3 = round(b_OX_m3, 3),
                     b_PR = round(b_PR, 3), pr_lo = round(pr_lo, 3),
                     pr_hi = round(pr_hi, 3), diff = round(diff, 3),
                     diff_lo = round(diff_lo, 3), diff_hi = round(diff_hi, 3)) %>%
    as.data.frame()

  # G-A4: the nesting identity, and delta's sign per cell
  x$nesting %>% as.data.frame()

  # the nine other genes (7.4), descriptive
  x$nine_gene_tally %>% as.data.frame()

  # Part C: reported, never read
  x$sensitivities %>% as.data.frame()
  x$c3 %>% as.data.frame()

  x$assertions %>% as.data.frame()
  cat(x$legend_draft$layout1, "\n")
}
