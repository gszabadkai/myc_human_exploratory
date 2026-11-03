# E49_signature_and_fig4g.R
# =============================================================================
# THE OXPHOS-ASSOCIATED BCL-2 SIGNATURE BESIDE OXPHOS, MYC ACTIVITY AND
# PROLIFERATION; THE THREE EFFECTORS; AND FIG. 4G AS A HEATMAP IN FIG. 1G'S
# LAYOUT.
#
# EXPLORATORY and POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-11-03_E49_signature_declaration.md, committed ALONE
# with data/annotations/ at 7bd7865, before this script existed. The governing
# brief, the author's decision record, is
# docs/2026-11-03_signature_and_fig4g_brief.md (eb45a9b), as amended by the
# author on 2026-11-03 (declaration section 3). Every choice below is the
# declaration's; section numbers cited are its sections.
#
# =============================================================================
# WHAT THIS DOES
# =============================================================================
# THERE IS NO READING RULE AND NO LABEL (section 8.1). Everything is reported.
#
# PART A FITS NOTHING (5.1). From E47's verified saved coefficients and
# variance matrices it recovers, for S6, S5 and S7:
#   - proliferation's m3 coefficient, and the difference (proliferation minus
#     OXPHOS, L V L') on the reference and least-entangled scores (6.1);
#   - MYC activity's coefficient at m1 and m3 on all four scores (6.2);
#   - OXPHOS's coefficient, which is drawn grey and NEVER READ (6.3, 8.2).
# It re-reads the twelve genes' values from E48's verified saved object.
# It evaluates the placement rule (8.3) and draws the branch it selects.
#
# PART B FITS ONLY THE THREE EFFECTORS (5.2), BAX, BAK1 AND BOK, with E47's
# exact model, behind two gates: the scoring code must reproduce the twelve
# genes in E47's frame (G-B1), and the fitting code must reproduce E47's 32
# BBC3 fits by digest (G-B2). The effectors are never added to any signature.
#
# Also: raw Spearman correlations (6.4), the effectors' expression percentiles
# (6.6), and the Fig. 1i-style display (6.8).
#
# DISPLAYS (section 9): Fig. 4g in two versions (unmarked; a dot where the CI
# excludes zero), and ED-e to ED-k. Every heatmap of model coefficients is
# rendered in both versions. The author chooses at STOP 3.
#
# =============================================================================
# TERMS
# =============================================================================
#   level, priority   OXPHOS subunits scored on GSVA (VST) and on mitoPPS
#                     (linear), named in plain words on every panel
#   estimator         a MYC activity score. Never a "signature"
#   M_a               FELSHER__MITOSTRIP, the reference score
#   clean             MYC_UP.V1_UP__MITOSTRIP, the least-entangled score
#   4g's model        the m3 fit on the least-entangled score (amendment 1)
#   PROLIF            PROLIF_DISJOINT (318 genes), E47's term "PROLIF"
#   m1, m3            E47's rungs: OX + MYC + spine; m1 + PROLIF
#   difference        beta(PROLIF) - beta(OX) within one m3 fit
#   the signature     S6 (E47 4a.7); variants S5 (no BID) and S7 (BMF -)
# No panel carries the notation beta, m1, m3, M_a or Delta (asserted).
#
# =============================================================================
# SNAPSHOT AND SCALE
# =============================================================================
# Reads E47's and E48's saved objects. Part B reads the snapshot matrices
# data/from_validation/{tcga_brca,scanb}_linear.rds, LINEAR DESeq2-normalised,
# and takes log2(x + 1), exactly as E47 and E34 did for the twelve; SCAN-B
# resolves through scanb_pheno.rds$symbol_map (trap 7). Each gene is
# Blom-scored ONCE over the whole cohort (1,095 / 3,207), then the fits subset
# to E34's analysis sets (938 / 3,143). Exposures, spine and sets are read
# from E47's frame and never rescored. The joint snapshot is not read.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It does not refit E47 or E48, except the effector fits and the BBC3
#     reproduction gate. E47's MIXED verdict stands.
#   - It reads no arm and draws none. It never reads the signature's OXPHOS
#     coefficient, or its difference, as evidence (8.2).
#   - It infers no binding class: only verified rows are drawn (section 10).
#   - It fits no product term and tests no interaction.
#   - It chooses no version of Fig. 4g. It renders both.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

message("\nE49: the signature beside OXPHOS, MYC activity and proliferation\n",
        strrep("=", 78))
message("EXPLORATORY, POST-HOC. Declaration 7bd7865; brief eb45a9b.")
message("No reading rule. Part A fits nothing. Part B fits the three effectors.")

PATH_E47  <- file.path(DIR_RESULTS, "e47_mutual_adjustment.rds")
PATH_E48  <- file.path(DIR_RESULTS, "e48_separation_and_fig4g.rds")
PATH_E49  <- file.path(DIR_RESULTS, "e49_signature_and_fig4g.rds")
PATH_SCS  <- file.path(DIR_RESULTS, "scanb_scores.rds")
PATH_SETS <- file.path(DIR_RESULTS, "set_definitions.rds")
PATH_BIND <- here::here("data", "annotations", "bh3_binding_code.csv")

DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-03"
# Files are named by CONTENT, never by figure letter (CLAUDE.md).
FIG <- c(fig4g_a  = "E49_fig4g_heatmap_unmarked",
         fig4g_b  = "E49_fig4g_heatmap_marked",
         ed_e_a   = "E49_heatmap_priority_unmarked",
         ed_e_b   = "E49_heatmap_priority_marked",
         ed_f_a   = "E49_heatmap_reference_score_unmarked",
         ed_f_b   = "E49_heatmap_reference_score_marked",
         ed_g     = "E49_signature_strips",
         ed_h_a   = "E49_raw_beside_adjusted_unmarked",
         ed_h_b   = "E49_raw_beside_adjusted_marked",
         ed_i     = "E49_fifteen_genes_intervals",
         ed_j     = "E49_myc_ladder_signatures_effectors",
         ed_k     = "E49_puma_bclxl_against_priority")

# =============================================================================
# 0. CONSTANTS - all from the declaration, none chosen here
# =============================================================================
COHORTS <- c("TCGA", "SCAN-B")
INSTR   <- c(gsva = "level", mitopps = "priority")
BLOCK_ORDER <- c("TCGA gsva", "TCGA mitopps", "SCAN-B gsva", "SCAN-B mitopps")

# MYC activity scores, in E47's order (E47 declaration section 5).
EST_ALL  <- c("clean", "M_a", "entangled", "regulon")
EST_4G   <- "clean"                      # 4g's model (amendment item 1)
EST_REF  <- "M_a"                        # the reference model
EST_NAME <- c(clean = "MYC_UP.V1_UP (1.5%)",
              M_a = "Felsher, reference (14.8%)",
              entangled = "Hallmark MYC targets (23.5%)",
              regulon = "CollecTRI MYC regulon")

# The signatures, E47's SIGS verbatim. MEMBERSHIP IS FIXED (section 4).
SIGS <- list(
  S6 = c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1),
  S5 = c(BBC3 = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1),
  S7 = c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, BMF = -1, MCL1 = -1))
SIG <- names(SIGS)
SIG_LAB <- c(S6 = "Signature (S6)", S5 = "S5, without BID",
             S7 = "S7, with BMF")

PRO    <- c("BAD", "BIK", "BBC3", "BID", "BCL2L11", "BMF", "PMAIP1")
ANTI   <- c("BCL2L1", "MCL1", "BCL2", "BCL2L2", "BCL2A1")
TWELVE <- c(PRO, ANTI)
EFF    <- c("BAX", "BAK1", "BOK")

# Fig. 1g's four groups and order (9.1).
GROUPS <- list(Effectors = c("BAX", "BAK1", "BOK"),
               Activators = c("BID", "BBC3", "BCL2L11"),
               Sensitisers = c("BAD", "PMAIP1", "BMF", "BIK"),
               `Anti-apoptotic` = c("MCL1", "BCL2L2", "BCL2A1", "BCL2",
                                    "BCL2L1"))
FIFTEEN  <- unlist(GROUPS, use.names = FALSE)
BH3_ONLY <- c(GROUPS$Activators, GROUPS$Sensitisers)
SECTION_GENES <- c("BAD", "BID")          # the section mark (9.1, 12)
PROTEIN <- c(BBC3 = "PUMA", BCL2L11 = "BIM", PMAIP1 = "NOXA",
             BCL2L2 = "BCL-w", BCL2A1 = "A1", BCL2L1 = "BCL-xL")

RUNG_COV <- list(
  TCGA     = list(m1 = c("PAM50", "purity", "leuko"),
                  m3 = c("PAM50", "purity", "leuko", "PROLIF")),
  `SCAN-B` = list(m1 = "PAM50", m3 = c("PAM50", "PROLIF")))
PAM50_LEVELS <- c("LumA", "LumB", "HER2", "Basal", "Normal")   # E34's order
EXPECT_N     <- c(TCGA = 938L, `SCAN-B` = 3143L)
EXPECT_COH   <- c(TCGA = 1095L, `SCAN-B` = 3207L)

# G-C1 (7.2): the reconciliation proposal's raw S6 correlations, 2 dp.
RECON_S6 <- tibble::tribble(
  ~cohort,  ~exposure,  ~rho_2dp,
  "TCGA",   "OX_gsva",   0.55,
  "SCAN-B", "OX_gsva",   0.40,
  "TCGA",   "MYC_M_a",   0.04,
  "SCAN-B", "MYC_M_a",   0.10,
  "TCGA",   "PROLIF",   -0.09,
  "SCAN-B", "PROLIF",   -0.04)

BAND_ID <- 0.80
TOL     <- 1e-10
CI_Q    <- 0.975
PCT_LOW <- 0.25            # E10's flag, reproduced for the twelve (6.6)

# Colours.
# - The heatmaps: E46's diverging scale, ColorBrewer RdBu's end points with
#   white at 0. The dataviz checker on the two poles, 2026-11-03:
#   "#2166AC,#B2182B" --mode light --pairs all: lightness band PASS; chroma
#   floor PASS; CVD worst all-pairs dE 21.1 (protan) PASS; normal-vision
#   dE 28.7 PASS; contrast PASS. ALL CHECKS PASS.
# - The signature's circular cells (8.2): a mid grey OFF the scale. Its least
#   CIE76 distance to any colour the scale can produce is 19.1 (checked
#   2026-11-03, before this script); asserted below to be at least 15.
# - The forests (ED-i, ED-j): E48's validated palette, unchanged.
FILL_LOW  <- "#2166AC"
FILL_HIGH <- "#B2182B"
GREY_CIRC <- "#A6A6A6"
COL_OX  <- "#009E73"
COL_PR  <- "#D55E00"
COL_MYC <- "#0072B2"
VALIDATED_PALETTE <- c(COL_OX, COL_PR, COL_MYC)

# Sizes: 178 mm wide, no text below 5.5 pt (6.7).
BASE_PT <- 6.5
MIN_PT  <- 5.5
TXT_MM  <- BASE_PT / ggplot2::.pt
SUB_MM  <- MIN_PT / ggplot2::.pt
W_MM    <- 178
W_PAGE  <- W_MM / 25.4
MINUS   <- "\u2212"
SECT    <- "\u00a7"

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.fmt3    <- function(e, lo, hi) sprintf("%+.3f [%+.3f, %+.3f]", e, lo, hi)
.dev     <- function(a, b) if (length(a)) max(abs(a - b)) else NA_real_
.fid     <- function(coh, inst, est, rung, rd)
  paste(coh, inst, est, rung, rd, sep = " | ")
.excl0   <- function(lo, hi) (lo > 0) | (hi < 0)

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

# E34's Blom, verbatim as E47 copies it, so a gene scored here is on exactly
# E34's scale.
.blom <- function(x) {
  n <- length(x)
  stopifnot(!anyNA(x))
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.375) / (n + 0.25))
  as.numeric(scale(z))
}

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

for (p in c(PATH_E47, PATH_E48, PATH_SCS, PATH_SETS, PATH_BIND, PATH_TCGA_MITO,
            PATH_SCANB_PHENO, PATH_TCGA_LINEAR, PATH_SCANB_LINEAR)) {
  .stop_if(file.exists(p), "absent: ", p)
}
e47 <- readRDS(PATH_E47)
x48 <- readRDS(PATH_E48)
NEED47 <- c("coefficient_matrix", "vcov", "results", "fit_digests", "frame",
            "gates")
.assert("E47's saved object holds every element E49 reads",
        all(NEED47 %in% names(e47)), n = length(NEED47))
NEED48 <- c("coefficients", "c3", "frame")
.assert("E48's saved object holds every element E49 reads",
        all(NEED48 %in% names(x48)), n = length(NEED48))
CM   <- e47$coefficient_matrix
VC   <- e47$vcov
R47  <- e47$results
DIGS <- e47$fit_digests
G3   <- e47$gates$g3

# The binding code (section 10): read, never inferred.
BIND <- utils::read.csv(PATH_BIND, stringsAsFactors = FALSE)
.assert("the binding-code table has its five declared columns",
        identical(names(BIND), c("gene", "binds_BCL2L1", "binds_MCL1", "source",
                                 "verified_by")), n = 5L)
.assert("the binding-code table covers exactly the seven BH3-only genes",
        setequal(BIND$gene, BH3_ONLY) && !any(duplicated(BIND$gene)),
        n = nrow(BIND))
.assert("binding classes are TRUE, FALSE or blank, and blank rows are blank",
        is.logical(BIND$binds_BCL2L1) && is.logical(BIND$binds_MCL1) &&
          all(is.na(BIND$binds_BCL2L1) == is.na(BIND$binds_MCL1)),
        n = nrow(BIND))
.assert("every verified row names its source and who verified it",
        all(nzchar(BIND$source[!is.na(BIND$binds_BCL2L1)])) &&
          all(nzchar(BIND$verified_by[!is.na(BIND$binds_BCL2L1)])),
        n = sum(!is.na(BIND$binds_BCL2L1)))
message("   binding code: ", sum(!is.na(BIND$binds_BCL2L1)), " verified rows (",
        paste(BIND$gene[!is.na(BIND$binds_BCL2L1)], collapse = ", "),
        "); blank, drawn blank: ",
        paste(BIND$gene[is.na(BIND$binds_BCL2L1)], collapse = ", "))

# =============================================================================
# 2. Structure (7.1), asserted again before anything is read
# =============================================================================
message("\n2. structure")

SIG_FITS <- unique(CM$fit_id[CM$readout %in% SIG])
.assert("E47 holds 144 signature fits: S6, S5, S7 x 3 rungs x 16",
        length(SIG_FITS) == 144L, n = length(SIG_FITS))
st_ok <- vapply(SIG_FITS, function(f) {
  tm <- CM$term[CM$fit_id == f]
  v  <- VC[[f]]
  !is.null(v) && nrow(v) == ncol(v) && identical(rownames(v), tm) &&
    identical(colnames(v), tm)
}, logical(1))
.assert(paste0("every signature fit has a square variance matrix named by its ",
               "terms, in order"), all(st_ok), n = length(SIG_FITS))
m3s <- SIG_FITS[grepl("\\| m3 \\|", SIG_FITS)]
.assert("PROLIF is in every signature m3 fit and in no other",
        all(vapply(m3s, function(f) "PROLIF" %in% CM$term[CM$fit_id == f],
                   logical(1))) &&
          !any(CM$term[CM$fit_id %in% setdiff(SIG_FITS, m3s)] == "PROLIF"),
        n = length(SIG_FITS))

# =============================================================================
# 3. Part A: the signatures, from E47's saved fits (5.1, 6.1-6.3)
# =============================================================================
message("\n3. Part A: the signatures")

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
# The general linear contrast, E48's: L carries +1 on PROLIF and -1 on OX, so
# Var = L V L' takes the covariance from the algebra.
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
.ga1 <- function(f) {
  dg <- DIGS[DIGS$fit_id == f$fid, ]
  tq <- stats::qt(CI_Q, f$df)
  tibble::tibble(
    fit_id = f$fid,
    coef_digest = identical(digest::digest(f$b), dg$digest_coef),
    vcov_digest = identical(digest::digest(f$V), dg$digest_vcov),
    se_dev = .dev(f$se, sqrt(diag(f$V))[names(f$se)]),
    ci_dev = max(.dev(f$lo, f$b - tq * f$se), .dev(f$hi, f$b + tq * f$se)))
}
# One row of quantities from an m1 and an m3 fit, with or without the
# difference (6.1 computes it on the reference and least-entangled scores).
.ext_row <- function(f1, f3, keys, with_diff, void) {
  ct <- .contrast(f3$b, f3$V, f3$df)
  alt_d  <- f3$b[["PROLIF"]] - f3$b[["OX"]]
  alt_se <- sqrt(f3$V["PROLIF", "PROLIF"] + f3$V["OX", "OX"] -
                   2 * f3$V["PROLIF", "OX"])
  na_if <- function(v) if (with_diff) v else NA_real_
  dplyr::bind_cols(keys, tibble::tibble(
    n = f3$n, df_m3 = f3$df, void = void,
    b_OX_m1 = f1$b[["OX"]], ox_lo_m1 = f1$lo[["OX"]], ox_hi_m1 = f1$hi[["OX"]],
    b_MYC_m1 = f1$b[["MYC"]], myc_lo_m1 = f1$lo[["MYC"]],
    myc_hi_m1 = f1$hi[["MYC"]],
    b_OX_m3 = f3$b[["OX"]], ox_lo_m3 = f3$lo[["OX"]], ox_hi_m3 = f3$hi[["OX"]],
    b_MYC_m3 = f3$b[["MYC"]], myc_lo_m3 = f3$lo[["MYC"]],
    myc_hi_m3 = f3$hi[["MYC"]],
    b_PR = f3$b[["PROLIF"]], pr_se = f3$se[["PROLIF"]],
    pr_lo = f3$lo[["PROLIF"]], pr_hi = f3$hi[["PROLIF"]],
    diff = na_if(ct[["diff"]]), diff_se = na_if(ct[["diff_se"]]),
    diff_lo = na_if(ct[["diff_lo"]]), diff_hi = na_if(ct[["diff_hi"]]),
    coef_cor = na_if(ct[["coef_cor"]]), naive_se = na_if(ct[["naive_se"]]),
    alt_dev = if (with_diff) max(abs(ct[["diff"]] - alt_d),
                                 abs(ct[["diff_se"]] - alt_se)) else 0))
}

GRID_S <- tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                             estimator = EST_ALL, readout = SIG)
EXT_SIG <- vector("list", nrow(GRID_S))
GA1 <- vector("list", nrow(GRID_S))
for (i in seq_len(nrow(GRID_S))) {
  g  <- GRID_S[i, ]
  f1 <- .get_fit(.fid(g$cohort, g$instrument, g$estimator, "m1", g$readout))
  f3 <- .get_fit(.fid(g$cohort, g$instrument, g$estimator, "m3", g$readout))
  GA1[[i]] <- dplyr::bind_rows(.ga1(f1), .ga1(f3))
  EXT_SIG[[i]] <- .ext_row(f1, f3, g,
                           with_diff = g$estimator %in% c(EST_4G, EST_REF),
                           void = isTRUE(R47$void[R47$fit_id == f3$fid]))
}
EXT_SIG <- dplyr::bind_rows(EXT_SIG)
GA1 <- dplyr::bind_rows(GA1) %>% dplyr::distinct(fit_id, .keep_all = TRUE)
.assert(paste0("G-A1: every signature coefficient vector read equals its fit's ",
               "coef(), by E47's digest"), all(GA1$coef_digest), n = nrow(GA1))
.assert(paste0("G-A1: every signature variance matrix read equals its fit's ",
               "vcov(), by E47's digest"), all(GA1$vcov_digest), n = nrow(GA1))
.assert("G-A1: every SE equals sqrt(diag(vcov)), every CI qt(0.975, df)",
        max(GA1$se_dev, GA1$ci_dev) < TOL, n = nrow(GA1),
        max_dev = max(GA1$se_dev, GA1$ci_dev))
.assert(paste0("G-A2: each signature difference recomputed from Var(P) + ",
               "Var(O) - 2 Cov(P, O) equals the L V L' contrast"),
        max(EXT_SIG$alt_dev) < 1e-12, n = sum(!is.na(EXT_SIG$diff)),
        max_dev = max(EXT_SIG$alt_dev))
.assert("the signature differences are on the reference and least-entangled scores only",
        all(!is.na(EXT_SIG$diff) == (EXT_SIG$estimator %in% c(EST_4G, EST_REF))),
        n = nrow(EXT_SIG))
# The re-read OXPHOS and MYC coefficients equal E47's $results exactly.
R47W <- R47 %>%
  dplyr::filter(rung %in% c("m1", "m3"), readout %in% SIG) %>%
  dplyr::transmute(cohort, instrument, estimator, readout, rung,
                   b_OX = unname(b_OX), b_MYC = unname(b_MYC)) %>%
  tidyr::pivot_wider(names_from = rung, values_from = c(b_OX, b_MYC))
CMPS <- EXT_SIG %>% dplyr::left_join(R47W, suffix = c("", "_r47"),
                                     by = c("cohort", "instrument", "estimator",
                                            "readout"))
DEV_S <- max(.dev(CMPS$b_OX_m1, CMPS$b_OX_m1_r47), .dev(CMPS$b_OX_m3, CMPS$b_OX_m3_r47),
             .dev(CMPS$b_MYC_m1, CMPS$b_MYC_m1_r47),
             .dev(CMPS$b_MYC_m3, CMPS$b_MYC_m3_r47))
.assert("the re-read signature OXPHOS and MYC coefficients equal E47's $results",
        DEV_S < TOL, n = 4L * nrow(CMPS), max_dev = DEV_S)
.assert("the only void signature cells are HALLMARK's m3 cells",
        all(EXT_SIG$void == (EXT_SIG$estimator == "entangled")), n = nrow(EXT_SIG))
.assert("membership is fixed: S6, S5 and S7 are E47's SIGS, and no effector is in any",
        identical(SIGS$S6, c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1,
                             MCL1 = -1)) &&
          !any(EFF %in% unlist(lapply(SIGS, names))), n = 3L)

# --- the placement rule (8.3), its inputs re-read from E47 -------------------
PLACEMENT <- EXT_SIG %>%
  dplyr::filter(readout == "S6", estimator %in% c("regulon", EST_4G)) %>%
  dplyr::select(cohort, instrument, estimator, b_MYC_m3) %>%
  tidyr::pivot_wider(names_from = estimator, values_from = b_MYC_m3) %>%
  dplyr::rename(myc_m3_regulon = regulon, myc_m3_least_entangled = clean) %>%
  dplyr::mutate(same_sign = sign(myc_m3_regulon) == sign(myc_m3_least_entangled))
KEEP_MYC <- all(PLACEMENT$same_sign)
.assert("the placement rule is evaluated in all four blocks", nrow(PLACEMENT) == 4L,
        n = 4L)
message("   placement rule (8.3), S6's m3 MYC coefficient, regulon against the ",
        "least-entangled score:")
.pr(PLACEMENT %>% dplyr::mutate(dplyr::across(dplyr::starts_with("myc"),
                                              ~ round(.x, 3))))
message("   -> ", if (KEEP_MYC) paste0("same sign in all four blocks: 4g keeps ",
                                       "the MYC column (4g's model).") else
          paste0("signs differ: 4g shows OXPHOS and proliferation from the ",
                 "reference model; MYC activity in Extended Data only."))

# =============================================================================
# 4. The twelve genes, re-read from E48 (G-A3)
# =============================================================================
message("\n4. the twelve genes, from E48")

G12 <- x48$coefficients %>%
  dplyr::filter(estimator %in% c(EST_4G, EST_REF), readout %in% TWELVE)
.assert("E48 holds the twelve on both scores in all four blocks",
        nrow(G12) == 96L && !anyNA(G12$b_PR), n = nrow(G12))
GA3 <- dplyr::bind_rows(lapply(seq_len(nrow(G12)), function(i) {
  g  <- G12[i, ]
  f1 <- .get_fit(.fid(g$cohort, g$instrument, g$estimator, "m1", g$readout))
  f3 <- .get_fit(.fid(g$cohort, g$instrument, g$estimator, "m3", g$readout))
  ct <- .contrast(f3$b, f3$V, f3$df)
  tibble::tibble(dev = max(
    abs(g$b_OX_m1 - f1$b[["OX"]]), abs(g$b_MYC_m1 - f1$b[["MYC"]]),
    abs(g$b_OX_m3 - f3$b[["OX"]]), abs(g$ox_lo_m3 - f3$lo[["OX"]]),
    abs(g$ox_hi_m3 - f3$hi[["OX"]]), abs(g$b_MYC_m3 - f3$b[["MYC"]]),
    abs(g$myc_lo_m3 - f3$lo[["MYC"]]), abs(g$myc_hi_m3 - f3$hi[["MYC"]]),
    abs(g$b_PR - f3$b[["PROLIF"]]), abs(g$pr_lo - f3$lo[["PROLIF"]]),
    abs(g$pr_hi - f3$hi[["PROLIF"]]), abs(g$diff - ct[["diff"]]),
    abs(g$diff_lo - ct[["diff_lo"]]), abs(g$diff_hi - ct[["diff_hi"]])))
}))
.assert("G-A3: E48's twelve-gene values equal E47's saved fits exactly",
        max(GA3$dev) < TOL, n = 14L * nrow(G12), max_dev = max(GA3$dev))

# =============================================================================
# 5. Part B: the effectors (5.2), behind G-B1 and G-B2
# =============================================================================
message("\n5. Part B: the effectors")

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(PATH_SCS)
sph  <- readRDS(PATH_SCANB_PHENO)
sd_  <- readRDS(PATH_SETS)
ID <- list(TCGA = colnames(mito$gsva_arms), `SCAN-B` = colnames(sc$gsva_arms))
.assert("the cohorts are E34's: 1,095 and 3,207 tumours",
        lengths(ID)[["TCGA"]] == EXPECT_COH[["TCGA"]] &&
          lengths(ID)[["SCAN-B"]] == EXPECT_COH[["SCAN-B"]], n = 2L)

# Percentiles (6.6), computed while each whole matrix is in memory:
#   pct_in_tumour_median: in each tumour, the gene's average-tie rank among all
#     genes divided by the number of genes; the median over the cohort;
#   pct_of_mean_e10: E10's definition, the rank of the gene's mean log2
#     expression among all genes' means, divided by the number of genes.
.percentiles <- function(G, rowmap) {
  n  <- nrow(G)
  e10 <- rank(rowMeans(G)) / n
  per <- matrix(NA_real_, length(rowmap), ncol(G))
  for (st in seq(1L, ncol(G), by = 256L)) {
    j <- st:min(ncol(G), st + 255L)
    C <- G[, j, drop = FALSE]
    for (i in seq_along(rowmap)) {
      v <- C[rowmap[[i]], ]
      less <- colSums(sweep(C, 2L, v, `<`))
      eq   <- colSums(sweep(C, 2L, v, `==`))
      per[i, j] <- (less + (eq + 1) / 2) / n
    }
  }
  tibble::tibble(gene = names(rowmap),
                 pct_in_tumour_median = apply(per, 1L, stats::median),
                 pct_of_mean_e10 = unname(e10[unname(rowmap)]))
}
# E47's .log2_twelve, generalised to any genes: log2(linear + 1), resolved
# through the cohort's symbol map, with every miss named.
.log2_genes <- function(path, ids, map, genes) {
  lin <- readRDS(path)
  .stop_if(identical(lin$scale, "linear_deseq2_normalised"),
           basename(path), " is not on the linear DESeq2-normalised scale")
  .stop_if(all(ids %in% colnames(lin$mat)), basename(path),
           ": a cohort sample is not a column")
  G  <- log2(lin$mat[, ids, drop = FALSE] + 1)
  rm(lin); invisible(gc(verbose = FALSE))
  rs <- .symbol_resolver(rownames(G), map)
  gr <- .gene_rows(genes, G, rs)
  .stop_if(length(gr$missing) == 0L, basename(path), ": did not resolve ",
           paste(gr$missing, collapse = ", "))
  rowmap <- vapply(genes, function(g) rs(g), character(1))
  pct <- .percentiles(G, rowmap)
  out <- list(mat = gr$mat[genes, , drop = FALSE], pct = pct)
  rm(G); invisible(gc(verbose = FALSE))
  out
}
GENES15 <- c(TWELVE, EFF)
GM <- list(TCGA     = .log2_genes(PATH_TCGA_LINEAR, ID$TCGA, NULL, GENES15),
           `SCAN-B` = .log2_genes(PATH_SCANB_LINEAR, ID$`SCAN-B`,
                                  sph$symbol_map, GENES15))
PCT <- dplyr::bind_rows(lapply(COHORTS, function(coh)
  GM[[coh]]$pct %>% dplyr::mutate(cohort = coh, .before = 1)))
.assert(paste0("the twelve clear E10's 25th expression percentile in both ",
               "cohorts, as E10 recorded"),
        all(PCT$pct_of_mean_e10[PCT$gene %in% TWELVE] >= PCT_LOW),
        n = sum(PCT$gene %in% TWELVE))

# G-B1: the scoring code, run on the twelve, reproduces E47's frame exactly.
FR47 <- e47$frame
GB1 <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  f <- FR47[FR47$cohort == coh, ]
  .stop_if(identical(f$sample_id, ID[[coh]]), coh,
           ": E47's frame is not in the score objects' order")
  tibble::tibble(cohort = coh, gene = TWELVE,
                 max_abs_diff = vapply(TWELVE, function(g)
                   max(abs(.blom(as.numeric(GM[[coh]]$mat[g, ])) - f[[g]])),
                   numeric(1)))
}))
.assert(paste0("G-B1: the gene-scoring code reproduces all twelve genes in ",
               "E47's frame, both cohorts"),
        all(GB1$max_abs_diff < TOL), n = nrow(GB1),
        max_dev = max(GB1$max_abs_diff))
message("   G-B1 passed: the effectors are scored by the same function.")

# The analysis frame (R rule 4): E47's frame plus the three effectors.
FRAME49 <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  f <- FR47[FR47$cohort == coh, ]
  for (g in EFF) f[[g]] <- .blom(as.numeric(GM[[coh]]$mat[g, ]))
  f
}))
.assert("the effectors are complete in both cohorts",
        !anyNA(FRAME49[, EFF]), n = 3L * nrow(FRAME49))

# G-B3, G2's rule: no effector in any fitted predictor set.
PRED_SETS <- list(
  MYC_UP.V1_UP__MITOSTRIP            = sd_$myc_sets[["MYC_UP.V1_UP__MITOSTRIP"]],
  FELSHER__MITOSTRIP                 = sd_$myc_sets[["FELSHER__MITOSTRIP"]],
  HALLMARK_MYC_TARGETS_V1__MITOSTRIP = sd_$myc_sets[["HALLMARK_MYC_TARGETS_V1__MITOSTRIP"]],
  M_b__MITOSTRIP                     = sd_$collectri_sets[["M_b__MITOSTRIP"]],
  `OXPHOS subunits`                  = mito$arm_sets[["OXPHOS subunits"]],
  PROLIF_DISJOINT                    = mito$covariate_sets[["PROLIF_DISJOINT"]])
.assert("G-B3: every predictor set is present, and none contains an effector",
        all(lengths(PRED_SETS) > 0L) &&
          !any(vapply(PRED_SETS, function(s) any(EFF %in% s), logical(1))),
        n = length(PRED_SETS))

# The fitting code: E47's .fit_one, on E47's frame, verbatim in its design.
DSET <- lapply(stats::setNames(COHORTS, COHORTS), function(coh) {
  d <- FRAME49[FRAME49$cohort == coh & FRAME49$in_set, ]
  d$PAM50 <- factor(d$PAM50, levels = PAM50_LEVELS)
  d
})
for (coh in COHORTS) {
  .assert(paste0(coh, ": the analysis set is E34's, n ", EXPECT_N[[coh]]),
          nrow(DSET[[coh]]) == EXPECT_N[[coh]] &&
            !anyNA(DSET[[coh]][, c("PAM50", "PROLIF", EFF)]),
          n = nrow(DSET[[coh]]))
}
.fit_e47 <- function(d, readout, inst, est, cov) {
  dd <- data.frame(y = d[[readout]], OX = d[[paste0("OX_", inst)]],
                   MYC = d[[paste0("MYC_", est)]])
  for (cv in cov) dd[[cv]] <- d[[cv]]
  if ("PAM50" %in% cov) dd$PAM50 <- droplevels(dd$PAM50)
  .stop_if(!anyNA(dd), "NA inside the analysis set: ", readout)
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
.as_fit <- function(f, fid) {
  td <- .tidy_fit(f, fid)
  list(fid = fid,
       b  = stats::setNames(td$estimate, td$term),
       se = stats::setNames(td$se, td$term),
       lo = stats::setNames(td$ci_lo, td$term),
       hi = stats::setNames(td$ci_hi, td$term),
       df = f$df.residual, n = stats::nobs(f), V = stats::vcov(f))
}

# G-B2: refitting BBC3 reproduces E47's 32 fits by digest.
GRID_B2 <- tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                              estimator = EST_ALL, rung = c("m1", "m3"))
GB2 <- dplyr::bind_rows(lapply(seq_len(nrow(GRID_B2)), function(i) {
  g   <- GRID_B2[i, ]
  fid <- .fid(g$cohort, g$instrument, g$estimator, g$rung, "BBC3")
  f   <- .fit_e47(DSET[[g$cohort]], "BBC3", g$instrument, g$estimator,
                  RUNG_COV[[g$cohort]][[g$rung]])
  dg  <- DIGS[DIGS$fit_id == fid, ]
  s   <- .get_fit(fid)
  tibble::tibble(fit_id = fid,
                 coef_digest = identical(digest::digest(stats::coef(f)),
                                         dg$digest_coef),
                 vcov_digest = identical(digest::digest(stats::vcov(f)),
                                         dg$digest_vcov),
                 coef_dev = .dev(stats::coef(f)[names(s$b)], s$b),
                 vcov_dev = max(abs(stats::vcov(f)[names(s$b), names(s$b)] - s$V)))
}))
.assert(paste0("G-B2: the fitting code reproduces E47's 32 BBC3 fits, ",
               "coefficients by digest"), all(GB2$coef_digest), n = nrow(GB2),
        max_dev = max(GB2$coef_dev))
.assert(paste0("G-B2: the fitting code reproduces E47's 32 BBC3 fits, variance ",
               "matrices by digest"), all(GB2$vcov_digest), n = nrow(GB2),
        max_dev = max(GB2$vcov_dev))
message("   G-B2 passed: the effectors are fitted by the same code.")

# The effector fits: 3 x 2 cohorts x 2 measures x 4 scores x {m1, m3} = 96.
VOID_M3 <- G3 %>% dplyr::distinct(cohort, estimator, void_m3)
.assert("G-B4: the m3 void flag depends on cohort and score only, and is HALLMARK's",
        nrow(VOID_M3) == 8L &&
          all(VOID_M3$void_m3 == (VOID_M3$estimator == "entangled")), n = 8L)
GRID_E <- tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                             estimator = EST_ALL, readout = EFF)
EFF_COEF <- list(); EFF_VCOV <- list(); EFF_VIF <- list(); EFF_DIG <- list()
FITS_EFF <- list(); EXT_EFF <- vector("list", nrow(GRID_E))
for (i in seq_len(nrow(GRID_E))) {
  g  <- GRID_E[i, ]
  fl <- list()
  for (rung in c("m1", "m3")) {
    fid <- .fid(g$cohort, g$instrument, g$estimator, rung, g$readout)
    f <- .fit_e47(DSET[[g$cohort]], g$readout, g$instrument, g$estimator,
                  RUNG_COV[[g$cohort]][[rung]])
    EFF_COEF[[fid]] <- .tidy_fit(f, fid)
    EFF_VCOV[[fid]] <- stats::vcov(f)
    EFF_VIF[[fid]]  <- .vif_tidy(f, fid)
    EFF_DIG[[fid]]  <- tibble::tibble(fit_id = fid,
                                      digest_coef = digest::digest(stats::coef(f)),
                                      digest_vcov = digest::digest(stats::vcov(f)))
    FITS_EFF[[fid]] <- f
    fl[[rung]] <- .as_fit(f, fid)
  }
  vd <- VOID_M3$void_m3[VOID_M3$cohort == g$cohort &
                          VOID_M3$estimator == g$estimator]
  EXT_EFF[[i]] <- .ext_row(fl$m1, fl$m3, g, with_diff = TRUE, void = vd)
}
EXT_EFF <- dplyr::bind_rows(EXT_EFF)
.assert("96 effector fits, each on its cohort's whole analysis set",
        length(FITS_EFF) == 96L &&
          all(vapply(FITS_EFF, stats::nobs, integer(1)) ==
                EXPECT_N[sub(" .*$", "", names(FITS_EFF))]),
        n = length(FITS_EFF))
.assert(paste0("G-A2 for the effectors: each difference recomputed from Var(P) ",
               "+ Var(O) - 2 Cov(P, O) equals L V L'"),
        max(EXT_EFF$alt_dev) < 1e-12, n = nrow(EXT_EFF),
        max_dev = max(EXT_EFF$alt_dev))
message("   96 effector fits. Their values are in the saved object; none is ",
        "printed here.")

# One table of every readout's m3 quantities, as the displays read them.
ALLX <- dplyr::bind_rows(
  EXT_SIG %>% dplyr::mutate(source = "E47, re-read"),
  G12 %>% dplyr::select(dplyr::any_of(names(EXT_SIG))) %>%
    dplyr::mutate(source = "E48, re-read"),
  EXT_EFF %>% dplyr::mutate(source = "E49, effector fit"))

# =============================================================================
# 6. Raw Spearman correlations (6.4), and G-C1
# =============================================================================
message("\n6. raw correlations")

EXPO <- c(OX_gsva = "OXPHOS level", OX_mitopps = "OXPHOS priority",
          MYC_clean = "MYC_UP.V1_UP", MYC_M_a = "Felsher, reference",
          MYC_entangled = "Hallmark MYC targets", MYC_regulon = "CollecTRI regulon",
          PROLIF = "Proliferation")
RAWC <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  d <- DSET[[coh]]
  tidyr::expand_grid(readout = c(SIG, FIFTEEN), exposure = names(EXPO)) %>%
    dplyr::mutate(cohort = coh, .before = 1) %>%
    dplyr::mutate(rho = mapply(function(r, e)
      stats::cor(d[[r]], d[[e]], method = "spearman"), readout, exposure))
}))
.assert("raw correlations: 18 readouts x 7 exposures x 2 cohorts",
        nrow(RAWC) == 252L && !anyNA(RAWC$rho), n = nrow(RAWC))
GC1 <- RECON_S6 %>%
  dplyr::left_join(RAWC %>% dplyr::filter(readout == "S6"),
                   by = c("cohort", "exposure")) %>%
  dplyr::mutate(dev = abs(rho - rho_2dp))
.assert(paste0("G-C1: the raw S6 correlations equal the reconciliation ",
               "proposal's, to 2 dp"), all(GC1$dev <= 0.005 + 1e-9),
        n = nrow(GC1), max_dev = max(GC1$dev))

# =============================================================================
# 7. The Fig. 1i-style display (6.8), and G-C2
# =============================================================================
message("\n7. the Fig. 1i-style relation")

FIG1I <- list(); FIG1I_RES <- list(); FIG1I_DEC <- list()
C3M3P <- x48$c3 %>% dplyr::filter(rung == "m3", instrument == "mitopps")
for (coh in COHORTS) {
  d   <- DSET[[coh]]
  cov <- RUNG_COV[[coh]][["m3"]]
  dd  <- data.frame(MYC = d$MYC_M_a)
  for (cv in cov) dd[[cv]] <- d[[cv]]
  dd$PAM50 <- droplevels(dd$PAM50)
  Xc <- stats::model.matrix(stats::reformulate(c("MYC", cov)), data = dd)
  y  <- d$BBC3 - d$BCL2L1
  yr <- stats::lm.fit(Xc, y)$residuals
  xr <- stats::lm.fit(Xc, d$OX_mitopps)$residuals
  slope <- sum(xr * yr) / sum(xr^2)
  b_bbc3   <- .get_fit(.fid(coh, "mitopps", EST_REF, "m3", "BBC3"))$b[["OX"]]
  b_bcl2l1 <- .get_fit(.fid(coh, "mitopps", EST_REF, "m3", "BCL2L1"))$b[["OX"]]
  .assert(paste0("G-C2, ", coh, ": the Fig. 1i-style slope equals E47's m3 ",
                 "priority OXPHOS coefficient for BBC3 minus BCL2L1"),
          abs(slope - (b_bbc3 - b_bcl2l1)) < TOL, n = 1L,
          max_dev = abs(slope - (b_bbc3 - b_bcl2l1)))
  c3 <- C3M3P[C3M3P$cohort == coh, ]
  .assert(paste0("G-C2, ", coh, ": C3's m3 priority quantity is minus this slope"),
          nrow(c3) == 1L && abs(slope + c3$delta) < TOL, n = 1L,
          max_dev = abs(slope + c3$delta))
  br  <- stats::quantile(xr, probs = seq(0, 1, 0.1), names = FALSE, type = 7)
  dec <- cut(xr, breaks = br, include.lowest = TRUE, labels = FALSE)
  FIG1I_RES[[coh]] <- tibble::tibble(sample_id = d$sample_id, cohort = coh,
                                     x_resid = xr, y_resid = yr, decile = dec)
  FIG1I_DEC[[coh]] <- FIG1I_RES[[coh]] %>%
    dplyr::group_by(cohort, decile) %>%
    dplyr::summarise(n = dplyr::n(), x = mean(x_resid), y = mean(y_resid),
                     .groups = "drop")
  FIG1I[[coh]] <- tibble::tibble(
    cohort = coh, n = nrow(d), slope = slope,
    slope_lo = -c3$ci_hi, slope_hi = -c3$ci_lo,
    partial_r = stats::cor(xr, yr),
    b_OX_BBC3 = b_bbc3, b_OX_BCL2L1 = b_bcl2l1,
    interval_source = "E48 C3: 2,000 paired bootstrap resamples, negated")
}
FIG1I     <- dplyr::bind_rows(FIG1I)
FIG1I_RES <- dplyr::bind_rows(FIG1I_RES)
FIG1I_DEC <- dplyr::bind_rows(FIG1I_DEC)
.assert("ten deciles per cohort", nrow(FIG1I_DEC) == 20L &&
          all(FIG1I_DEC$n > 0L), n = nrow(FIG1I_DEC))

# =============================================================================
# 8. Displays (section 9), content fixed before any quantity was read
# =============================================================================
message("\n8. displays")

FIG_CHECKS <- list()
.lin <- function(h) {
  m <- grDevices::col2rgb(h) / 255
  f <- function(u) ifelse(u <= 0.03928, u / 12.92, ((u + 0.055) / 1.055)^2.4)
  as.numeric(0.2126 * f(m[1, ]) + 0.7152 * f(m[2, ]) + 0.0722 * f(m[3, ]))
}
.ink <- function(h) ifelse(.lin(h) < 0.30, "white", "black")
.sym_lim <- function(v) {
  L <- ceiling(max(abs(v)) * 10 - 1e-9) / 10
  c(-L, L)
}
.sym_breaks <- function(L) {
  st <- if (L <= 0.6) 0.2 else if (L <= 1.2) 0.4 else 0.5
  k <- floor(L / st + 1e-9)
  round(seq(-k, k) * st, 10)
}
.scale_cols <- function(v, lim) {
  pal <- scales::div_gradient_pal(FILL_LOW, "white", FILL_HIGH, "Lab")
  toupper(pal(scales::rescale_mid(v, to = c(0, 1), from = lim, mid = 0)))
}

# The grey is off the scale (8.2, 7.2).
.lab <- function(h) grDevices::convertColor(t(grDevices::col2rgb(h) / 255),
                                            from = "sRGB", to = "Lab")
GRID_COLS <- .scale_cols(seq(-1, 1, length.out = 2001), c(-1, 1))
DE_GREY <- min(sqrt(rowSums((.lab(GRID_COLS) -
                               matrix(.lab(GREY_CIRC), length(GRID_COLS), 3,
                                      byrow = TRUE))^2)))
.assert("the circular cells' grey is off the diverging scale (least dE >= 15)",
        DE_GREY >= 15, n = length(GRID_COLS), max_dev = DE_GREY)

# --- plain-word labels, drawn with plotmath so the source stays ASCII -------
.gene_lab <- function(g, sect = TRUE) {
  base <- paste0('italic("', g, '")')
  if (g %in% names(PROTEIN)) base <- paste0(base, '~"(', PROTEIN[[g]], ')"')
  if (sect && g %in% SECTION_GENES) base <- paste0(base, '~"', SECT, '"')
  base
}
BLOCK_LAB <- c("TCGA gsva" = '"TCGA" %.% "level"',
               "TCGA mitopps" = '"TCGA" %.% "priority"',
               "SCAN-B gsva" = '"SCAN-B" %.% "level"',
               "SCAN-B mitopps" = '"SCAN-B" %.% "priority"')
TERM_TITLE <- c(OX = "OXPHOS", MYC = "MYC activity", PROLIF = "Proliferation",
                DIFF = "Proliferation minus OXPHOS")

# Every text element of the millimetre-coordinate panels is set explicitly:
# theme_void leaves some at rel(0.8) of the base size, which is below 5.5 pt.
.theme_hm <- function(legend_inside = c(0.985, 0.012)) {
  ggplot2::theme_void(base_size = BASE_PT) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = legend_inside,
      legend.justification = c(1, 0),
      legend.text = ggplot2::element_text(size = MIN_PT),
      legend.title = ggplot2::element_text(size = MIN_PT),
      plot.title = ggplot2::element_text(size = BASE_PT),
      plot.subtitle = ggplot2::element_text(size = MIN_PT),
      plot.caption = ggplot2::element_text(size = MIN_PT),
      strip.text = ggplot2::element_text(size = BASE_PT),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(0, 0, 0, 0))
}

# --- the heatmap, in Fig. 1g's layout (9.1). Coordinates are millimetres. ----
RH <- 4.4; GG <- 1.8; SH <- 7.6
X_GRP <- 3.2; X_GLINE <- 6.4; X_GENE <- 44.0
X_B1 <- 49.5; X_B2 <- 58.5; X_SG <- 69.0
X0 <- 76.5; BG <- 3.0; HDR_TOP <- -3.6; HDR_COH <- -8.4; Y_FIRST <- -13.0
ROWS15 <- tibble::tibble(readout = FIFTEEN,
                         group = rep(names(GROUPS), lengths(GROUPS))) %>%
  dplyr::mutate(k = dplyr::row_number(), g = match(group, names(GROUPS)),
                y = Y_FIRST - (k - 1) * RH - (g - 1) * GG)
Y_LAST  <- min(ROWS15$y)
Y_STRIP <- Y_LAST - RH / 2 - 3.2 - SH / 2

.cols <- function(terms, x0 = X0, width = W_MM - X0 - 3) {
  nb <- length(terms)
  cw <- (width - (nb - 1) * BG) / (2 * nb)
  tidyr::expand_grid(term = terms, cohort = COHORTS) %>%
    dplyr::mutate(b = match(term, terms), c = match(cohort, COHORTS),
                  x = x0 + (b - 1) * (2 * cw + BG) + (c - 0.5) * cw, cw = cw)
}
.cells <- function(model, inst, terms) {
  X <- ALLX %>% dplyr::filter(estimator == model, instrument == inst,
                              readout %in% FIFTEEN)
  .stop_if(nrow(X) == 30L, "a heatmap is missing cells: ", model, " ", inst)
  dplyr::bind_rows(lapply(terms, function(tm) {
    col <- switch(tm, OX = c("b_OX_m3", "ox_lo_m3", "ox_hi_m3"),
                  MYC = c("b_MYC_m3", "myc_lo_m3", "myc_hi_m3"),
                  PROLIF = c("b_PR", "pr_lo", "pr_hi"))
    X %>% dplyr::transmute(readout, cohort, term = tm, est = .data[[col[1]]],
                           lo = .data[[col[2]]], hi = .data[[col[3]]])
  }))
}
.strip_cells <- function(model, inst, terms, rd = "S6") {
  X <- EXT_SIG %>% dplyr::filter(estimator == model, instrument == inst,
                                 readout == rd)
  .stop_if(nrow(X) == 2L, "a strip is missing cells")
  dplyr::bind_rows(lapply(terms, function(tm) {
    col <- switch(tm, OX = c("b_OX_m3", "ox_lo_m3", "ox_hi_m3"),
                  MYC = c("b_MYC_m3", "myc_lo_m3", "myc_hi_m3"),
                  PROLIF = c("b_PR", "pr_lo", "pr_hi"),
                  DIFF = c("diff", "diff_lo", "diff_hi"))
    X %>% dplyr::transmute(readout, cohort, term = tm, est = .data[[col[1]]],
                           lo = .data[[col[2]]], hi = .data[[col[3]]],
                           grey = tm %in% c("OX", "DIFF"))
  }))
}
.cell_txt <- function(e, lo, hi) sprintf("%.2f\n(%.2f, %.2f)", e, lo, hi)

.heatmap <- function(model, inst, terms, marked, note) {
  cols  <- .cols(terms)
  cw    <- cols$cw[1]
  cells <- .cells(model, inst, terms) %>%
    dplyr::left_join(cols %>% dplyr::select(term, cohort, x),
                     by = c("term", "cohort")) %>%
    dplyr::left_join(ROWS15 %>% dplyr::select(readout, y), by = "readout")
  strip <- .strip_cells(model, inst, terms) %>%
    dplyr::left_join(cols %>% dplyr::select(term, cohort, x),
                     by = c("term", "cohort")) %>%
    dplyr::mutate(y = Y_STRIP, label = .cell_txt(est, lo, hi))
  lim <- .sym_lim(c(cells$est, strip$est[!strip$grey]))
  strip <- strip %>%
    dplyr::mutate(fill_hex = ifelse(grey, GREY_CIRC, .scale_cols(est, lim)),
                  ink = ifelse(grey, "grey15", .ink(fill_hex)))
  dots <- cells %>% dplyr::filter(.excl0(lo, hi))
  hdr_b <- cols %>% dplyr::group_by(term) %>%
    dplyr::summarise(x = mean(x), .groups = "drop") %>%
    dplyr::mutate(label = unname(TERM_TITLE[term]))
  bind <- BIND %>% dplyr::filter(!is.na(binds_BCL2L1)) %>%
    dplyr::inner_join(ROWS15 %>% dplyr::select(readout, y),
                      by = c("gene" = "readout"))
  bind_pts <- dplyr::bind_rows(
    bind %>% dplyr::transmute(gene, y, x = X_B1, binds = binds_BCL2L1),
    bind %>% dplyr::transmute(gene, y, x = X_B2, binds = binds_MCL1)) %>%
    dplyr::mutate(fill = ifelse(binds, "black", "white"))
  signs <- tibble::tibble(readout = names(SIGS$S6),
                          label = ifelse(SIGS$S6 > 0, "+", MINUS)) %>%
    dplyr::left_join(ROWS15 %>% dplyr::select(readout, y), by = "readout")
  grp <- ROWS15 %>% dplyr::group_by(group) %>%
    dplyr::summarise(top = max(y) + RH / 2 - 0.4, bot = min(y) - RH / 2 + 0.4,
                     mid = mean(y), .groups = "drop")
  H <- -(Y_STRIP - SH / 2 - 16)
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = .tag(cells, "cells"),
                       ggplot2::aes(x = x, y = y, fill = est),
                       width = cw - 0.6, height = RH - 0.6) +
    ggplot2::geom_tile(data = .tag(strip[!strip$grey, ], "strip_col"),
                       ggplot2::aes(x = x, y = y, fill = est),
                       width = cw - 0.6, height = SH) +
    ggplot2::geom_tile(data = .tag(strip[strip$grey, ], "strip_grey"),
                       ggplot2::aes(x = x, y = y),
                       fill = GREY_CIRC, width = cw - 0.6, height = SH) +
    ggplot2::geom_text(data = .tag(strip, "strip_text"),
                       ggplot2::aes(x = x, y = y, label = label, colour = ink),
                       size = SUB_MM, lineheight = 0.85) +
    ggplot2::scale_colour_identity() +
    ggplot2::annotate("text", x = mean(strip$x[strip$grey]),
                      y = Y_STRIP - SH / 2 - 2.2, label = note,
                      size = SUB_MM, colour = "grey25") +
    ggplot2::annotate("text", x = X_GENE, y = Y_STRIP, label = "Signature",
                      hjust = 1, size = TXT_MM, fontface = "bold") +
    ggplot2::geom_text(data = hdr_b, ggplot2::aes(x = x, y = HDR_TOP, label = label),
                       size = TXT_MM, fontface = "bold") +
    ggplot2::geom_text(data = cols, ggplot2::aes(x = x, y = HDR_COH, label = cohort),
                       size = SUB_MM, colour = "grey20") +
    ggplot2::annotate("text", x = c(X_B1, X_B2, X_SG), y = HDR_COH + 1.6,
                      label = c("binds\nBCL-xL", "binds\nMCL1", "signature\nsign"),
                      size = SUB_MM, fontface = "bold", lineheight = 0.85) +
    ggplot2::geom_text(data = ROWS15, ggplot2::aes(x = X_GENE, y = y,
                                                   label = vapply(readout,
                                                                  .gene_lab,
                                                                  character(1))),
                       hjust = 1, size = TXT_MM, parse = TRUE) +
    ggplot2::geom_point(data = .tag(bind_pts, "bind"), ggplot2::aes(x = x, y = y),
                        fill = bind_pts$fill, shape = 21, size = 1.7,
                        stroke = 0.4, colour = "black") +
    ggplot2::geom_text(data = signs, ggplot2::aes(x = X_SG, y = y, label = label),
                       size = TXT_MM + 0.6, fontface = "bold") +
    ggplot2::geom_segment(data = grp, ggplot2::aes(x = X_GLINE, xend = X_GLINE,
                                                   y = top, yend = bot),
                          linewidth = 0.4) +
    ggplot2::geom_text(data = grp, ggplot2::aes(x = X_GRP, y = mid, label = group),
                       angle = 90, size = SUB_MM + 0.15, fontface = "bold")
  if (marked) {
    p <- p + ggplot2::geom_point(data = .tag(dots, "dots"),
                                 ggplot2::aes(x = x, y = y),
                                 shape = 21, size = 0.9, fill = "black",
                                 colour = "white", stroke = 0.35) +
      ggplot2::annotate("point", x = X0 + 1.5, y = -H + 6.8, shape = 21,
                        size = 0.9, fill = "black", colour = "white",
                        stroke = 0.35) +
      ggplot2::annotate("text", x = X0 + 3.2, y = -H + 6.8, hjust = 0,
                        label = "95% CI excludes zero", size = SUB_MM)
  }
  p <- p +
    ggplot2::scale_fill_gradient2(
      low = FILL_LOW, mid = "white", high = FILL_HIGH, midpoint = 0,
      limits = lim, breaks = .sym_breaks(lim[2]),
      name = "SD per SD of each score, with the other two held fixed",
      guide = ggplot2::guide_colourbar(
        direction = "horizontal", theme = ggplot2::theme(
          legend.title.position = "top",
          legend.key.width = grid::unit(40, "mm"),
          legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_MM), ylim = c(-H, 0),
                         expand = FALSE, clip = "off") +
    .theme_hm()
  list(plot = p, cells = cells, strip = strip, dots = dots, bind = bind_pts,
       signs = signs, lim = lim, cols = cols, height_mm = H,
       model = model, instrument = inst, terms = terms, marked = marked)
}

# --- the checks every heatmap must pass (7.2) -------------------------------
# Each layer's data carries a .layer tag; a layer is found by its tag.
.tag <- function(d, tag) dplyr::mutate(d, .layer = tag)
.layer_of <- function(p, tag) {
  k <- which(vapply(p$layers, function(l) is.data.frame(l$data) &&
                      ".layer" %in% names(l$data) &&
                      identical(unique(as.character(l$data$.layer)), tag),
                    logical(1)))
  .stop_if(length(k) == 1L, "layer not unique: ", tag)
  k
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
                                              "x", "y", "fill"))])
  for (l in p$layers) {
    if (!is.null(l$aes_params$label)) out <- c(out, as.character(l$aes_params$label))
    m <- l$mapping$label
    if (!is.null(m)) {
      dd <- if (is.data.frame(l$data)) l$data else p$data
      v <- tryCatch(rlang::eval_tidy(m, dd), error = function(e) NULL)
      out <- c(out, as.character(v))
    }
  }
  for (s in p$scales$scales) {
    if (is.character(s$labels)) out <- c(out, s$labels)
    if (is.character(s$name)) out <- c(out, s$name)
  }
  as.character(out)
}
.check_text <- function(parts, label, width_in) {
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
.check_heatmap <- function(hm, label) {
  p <- hm$plot
  b <- ggplot2::ggplot_build(p)
  kc <- .layer_of(p, "cells")
  dc <- b$data[[kc]]
  key_s <- paste(round(hm$cells$x, 6), round(hm$cells$y, 6))
  m <- match(paste(round(dc$x, 6), round(dc$y, 6)), key_s)
  .assert(paste0(label, ": every cell drawn is one stored cell, none dropped"),
          nrow(dc) == nrow(hm$cells) && !anyNA(m) && !any(duplicated(m)),
          n = nrow(dc))
  .assert(paste0(label, ": every cell's colour is its stored value on the scale"),
          identical(toupper(dc$fill), .scale_cols(hm$cells$est[m], hm$lim)),
          n = nrow(dc))
  sv <- c(hm$cells$est, hm$strip$est[!hm$strip$grey])
  .assert(paste0(label, ": the limits are symmetric, contain every coloured value ",
                 "and are rounded outward to 0.1"),
          abs(hm$lim[1] + hm$lim[2]) < 1e-12 &&
            all(sv >= hm$lim[1] & sv <= hm$lim[2]) &&
            abs(hm$lim[2] * 10 - round(hm$lim[2] * 10)) < 1e-9,
          n = length(sv), detail = sprintf("limits +/- %.1f", hm$lim[2]))
  ks <- .layer_of(p, "strip_col")
  ds <- b$data[[ks]]
  .assert(paste0(label, ": the strip's coloured cells carry their stored values"),
          identical(toupper(ds$fill),
                    .scale_cols(hm$strip$est[!hm$strip$grey], hm$lim)),
          n = nrow(ds))
  kg <- .layer_of(p, "strip_grey")
  .assert(paste0(label, ": the circular cells are grey, off the scale"),
          all(toupper(b$data[[kg]]$fill) == GREY_CIRC) &&
            nrow(b$data[[kg]]) == sum(hm$strip$grey) &&
            all(hm$strip$term[hm$strip$grey] %in% c("OX", "DIFF")),
          n = nrow(b$data[[kg]]))
  kt <- .layer_of(p, "strip_text")
  .assert(paste0(label, ": the strip's text is its stored estimate and interval"),
          identical(as.character(b$data[[kt]]$label),
                    .cell_txt(hm$strip$est, hm$strip$lo, hm$strip$hi)),
          n = nrow(hm$strip))
  kb <- .layer_of(p, "bind")
  db <- b$data[[kb]]
  .assert(paste0(label, ": binding marks only on verified BH3-only rows, filled ",
                 "= binds"),
          nrow(db) == 2L * sum(!is.na(BIND$binds_BCL2L1)) &&
            all(hm$bind$gene %in% BH3_ONLY) &&
            !any(hm$bind$gene %in% BIND$gene[is.na(BIND$binds_BCL2L1)]) &&
            identical(as.character(db$fill), hm$bind$fill),
          n = nrow(db))
  .assert(paste0(label, ": the signature-sign column is S6's, from E47"),
          identical(hm$signs$readout, names(SIGS$S6)) &&
            identical(hm$signs$label, ifelse(SIGS$S6 > 0, "+", MINUS)),
          n = length(SIGS$S6))
  lab_all <- vapply(FIFTEEN, .gene_lab, character(1))
  .assert(paste0(label, ": the section mark is on BAD and BID only"),
          identical(names(lab_all)[grepl(SECT, lab_all, fixed = TRUE)],
                    SECTION_GENES[order(match(SECTION_GENES, FIFTEEN))]),
          n = length(FIFTEEN))
  if (hm$marked) {
    kd <- .layer_of(p, "dots")
    dd <- b$data[[kd]]
    exp_d <- hm$cells %>% dplyr::filter(.excl0(lo, hi))
    .assert(paste0(label, ": a dot exactly where the 95% CI excludes zero"),
            nrow(dd) == nrow(exp_d) &&
              setequal(paste(round(dd$x, 6), round(dd$y, 6)),
                       paste(round(exp_d$x, 6), round(exp_d$y, 6))),
            n = nrow(hm$cells))
  }
  .check_text(list(p), label, W_PAGE)
  invisible(TRUE)
}
.same_values <- function(a, b, label) {
  .assert(paste0(label, ": the two versions draw identical values"),
          identical(a$cells, b$cells) && identical(a$strip, b$strip) &&
            identical(a$lim, b$lim), n = nrow(a$cells) + nrow(a$strip))
}

# --- Fig. 4g, by the placement rule (8.3) ------------------------------------
NOTE_GREY <- "grey: defines the signature"
if (KEEP_MYC) {
  SPEC_4G <- list(model = EST_4G, terms = c("OX", "MYC", "PROLIF"))
} else {
  SPEC_4G <- list(model = EST_REF, terms = c("OX", "PROLIF"))
}
HM <- list()
HM$fig4g_a <- .heatmap(SPEC_4G$model, "gsva", SPEC_4G$terms, FALSE, NOTE_GREY)
HM$fig4g_b <- .heatmap(SPEC_4G$model, "gsva", SPEC_4G$terms, TRUE, NOTE_GREY)
HM$ed_e_a  <- .heatmap(SPEC_4G$model, "mitopps", SPEC_4G$terms, FALSE, NOTE_GREY)
HM$ed_e_b  <- .heatmap(SPEC_4G$model, "mitopps", SPEC_4G$terms, TRUE, NOTE_GREY)
HM$ed_f_a  <- .heatmap(EST_REF, "gsva", c("OX", "MYC", "PROLIF"), FALSE, NOTE_GREY)
HM$ed_f_b  <- .heatmap(EST_REF, "gsva", c("OX", "MYC", "PROLIF"), TRUE, NOTE_GREY)
for (nm in names(HM)) .check_heatmap(HM[[nm]], nm)
.same_values(HM$fig4g_a, HM$fig4g_b, "Fig. 4g")
.same_values(HM$ed_e_a, HM$ed_e_b, "ED-e")
.same_values(HM$ed_f_a, HM$ed_f_b, "ED-f")
.assert("Fig. 4g draws the branch the placement rule selected",
        identical(HM$fig4g_a$model, SPEC_4G$model) &&
          identical(HM$fig4g_a$terms, SPEC_4G$terms), n = 1L,
        detail = paste0("model ", SPEC_4G$model, "; ",
                        paste(SPEC_4G$terms, collapse = ", ")))

# --- ED-g: S5, S6 and S7 strips side by side (9.2) --------------------------
STRIP_TERMS <- c("OX", "MYC", "PROLIF", "DIFF")
EDG_ROWS <- tidyr::expand_grid(instrument = names(INSTR),
                               model = c(EST_4G, EST_REF), readout = SIG) %>%
  dplyr::mutate(p = match(instrument, names(INSTR)),
                r = (match(model, c(EST_4G, EST_REF)) - 1) * 3 +
                  match(readout, SIG),
                y = -16 - (p - 1) * 64 - (r - 1) * 8.4 -
                  ifelse(r > 3, 2.6, 0))
EDG <- dplyr::bind_rows(lapply(seq_len(nrow(EDG_ROWS)), function(i) {
  r <- EDG_ROWS[i, ]
  .strip_cells(r$model, r$instrument, STRIP_TERMS, rd = r$readout) %>%
    dplyr::mutate(instrument = r$instrument, model = r$model, y = r$y)
}))
EDG_COLS <- .cols(STRIP_TERMS, x0 = 50, width = W_MM - 50 - 3)
EDG <- EDG %>% dplyr::left_join(EDG_COLS %>% dplyr::select(term, cohort, x),
                                by = c("term", "cohort"))
LIM_G <- .sym_lim(EDG$est[!EDG$grey])
EDG <- EDG %>%
  dplyr::mutate(fill_hex = ifelse(grey, GREY_CIRC, .scale_cols(est, LIM_G)),
                ink = ifelse(grey, "grey15", .ink(fill_hex)),
                label = .cell_txt(est, lo, hi))
.assert("ED-g: 3 signatures x 2 models x 2 measures x 8 cells",
        nrow(EDG) == 96L && !anyNA(EDG$x), n = nrow(EDG))
CW_G <- EDG_COLS$cw[1]
H_G <- 140
EDG_HDR <- dplyr::bind_rows(lapply(names(INSTR), function(inst) {
  top <- -16 - (match(inst, names(INSTR)) - 1) * 64 + 8.4
  dplyr::bind_rows(
    EDG_COLS %>% dplyr::group_by(term) %>%
      dplyr::summarise(x = mean(x), .groups = "drop") %>%
      dplyr::mutate(y = top + 1.8, label = unname(TERM_TITLE[term]), bold = TRUE),
    EDG_COLS %>% dplyr::transmute(term, x, y = top - 2.2, label = cohort,
                                  bold = FALSE))
}))
EDG_RL <- EDG_ROWS %>% dplyr::mutate(
  label = paste0(SIG_LAB[readout], ", ",
                 ifelse(model == EST_4G, "least-entangled score",
                        "reference score")))
EDG_PANEL <- tibble::tibble(
  y = -16 - (0:1) * 64 + 13.2,
  label = paste0("OXPHOS ", unname(INSTR)))
P_EDG <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = .tag(EDG[!EDG$grey, ], "g_col"),
                     ggplot2::aes(x = x, y = y, fill = est),
                     width = CW_G - 0.6, height = 7.6) +
  ggplot2::geom_tile(data = .tag(EDG[EDG$grey, ], "g_grey"),
                     ggplot2::aes(x = x, y = y),
                     fill = GREY_CIRC, width = CW_G - 0.6, height = 7.6) +
  ggplot2::geom_text(data = .tag(EDG, "g_text"),
                     ggplot2::aes(x = x, y = y, label = label, colour = ink),
                     size = SUB_MM, lineheight = 0.85) +
  ggplot2::scale_colour_identity() +
  ggplot2::geom_text(data = EDG_HDR[EDG_HDR$bold, ],
                     ggplot2::aes(x = x, y = y, label = label),
                     size = TXT_MM, fontface = "bold") +
  ggplot2::geom_text(data = EDG_HDR[!EDG_HDR$bold, ],
                     ggplot2::aes(x = x, y = y, label = label),
                     size = SUB_MM, colour = "grey20") +
  ggplot2::geom_text(data = EDG_RL, ggplot2::aes(x = 48, y = y, label = label),
                     hjust = 1, size = SUB_MM) +
  ggplot2::geom_text(data = EDG_PANEL, ggplot2::aes(x = 2, y = y, label = label),
                     hjust = 0, size = TXT_MM, fontface = "bold") +
  ggplot2::annotate("text", x = 2, y = -H_G + 4.5, hjust = 0, size = SUB_MM,
                    colour = "grey25",
                    label = paste0("grey: contains the coefficient that defines ",
                                   "the signature")) +
  ggplot2::scale_fill_gradient2(
    low = FILL_LOW, mid = "white", high = FILL_HIGH, midpoint = 0,
    limits = LIM_G, breaks = .sym_breaks(LIM_G[2]),
    name = "SD per SD of each score, with the other two held fixed",
    guide = ggplot2::guide_colourbar(direction = "horizontal",
                                     theme = ggplot2::theme(
                                       legend.title.position = "top",
                                       legend.key.width = grid::unit(40, "mm"),
                                       legend.key.height = grid::unit(2.2, "mm")))) +
  ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_MM), ylim = c(-H_G, 0),
                       expand = FALSE, clip = "off") +
  .theme_hm(c(0.985, 0.005))
b_g <- ggplot2::ggplot_build(P_EDG)
kgc <- .layer_of(P_EDG, "g_col")
.assert("ED-g: every coloured strip cell carries its stored value on the scale",
        identical(toupper(b_g$data[[kgc]]$fill),
                  .scale_cols(EDG$est[!EDG$grey], LIM_G)),
        n = sum(!EDG$grey))
.assert("ED-g: OXPHOS and the difference are grey, the circular cells (8.2)",
        all(EDG$term[EDG$grey] %in% c("OX", "DIFF")) &&
          all(toupper(b_g$data[[.layer_of(P_EDG, "g_grey")]]$fill) ==
                GREY_CIRC), n = sum(EDG$grey))
.assert("ED-g: every cell's text is its stored estimate and interval",
        identical(as.character(b_g$data[[.layer_of(P_EDG, "g_text")]]$label),
                  .cell_txt(EDG$est, EDG$lo, EDG$hi)), n = nrow(EDG))
.check_text(list(P_EDG), "ED-g", W_PAGE)

# --- ED-h: raw correlations beside the adjusted coefficients (9.2) ----------
ROWS18 <- tibble::tibble(readout = c(SIG, FIFTEEN),
                         group = c(rep("Signature", 3),
                                   rep(names(GROUPS), lengths(GROUPS)))) %>%
  dplyr::mutate(k = dplyr::row_number(),
                g = match(group, c("Signature", names(GROUPS))),
                y = -32 - (k - 1) * 4.0 - (g - 1) * 1.6)
.row_lab18 <- function(rd) if (rd %in% SIG) paste0('"', SIG_LAB[[rd]], '"') else
  .gene_lab(rd)
H_H <- -(min(ROWS18$y) - 2 - 13)
# Left: raw rho, grouped by cohort (7 exposures each).
RAW_COLS <- tidyr::expand_grid(cohort = COHORTS, exposure = names(EXPO)) %>%
  dplyr::mutate(c = match(cohort, COHORTS), e = match(exposure, names(EXPO)),
                x = 36 + (c - 1) * (7 * 6.0 + 3) + (e - 0.5) * 6.0)
RAW_T <- RAWC %>% dplyr::inner_join(RAW_COLS, by = c("cohort", "exposure")) %>%
  dplyr::inner_join(ROWS18 %>% dplyr::select(readout, y), by = "readout")
LIM_R <- .sym_lim(RAW_T$rho)
RAW_T <- RAW_T %>% dplyr::mutate(fill_hex = .scale_cols(rho, LIM_R),
                                 ink = .ink(fill_hex),
                                 label = sprintf("%.2f", rho))
W_HL <- max(RAW_COLS$x) + 3 + 2
# Right: 4g's adjusted coefficients (4g's model, level), signatures included.
ADJ_COLS <- tidyr::expand_grid(cohort = COHORTS, term = c("OX", "MYC", "PROLIF")) %>%
  dplyr::mutate(c = match(cohort, COHORTS), t = match(term, c("OX", "MYC", "PROLIF")),
                x = 3 + (c - 1) * (3 * 6.4 + 3) + (t - 0.5) * 6.4)
ADJ_T <- dplyr::bind_rows(
  .cells(EST_4G, "gsva", c("OX", "MYC", "PROLIF")),
  dplyr::bind_rows(lapply(SIG, function(s)
    .strip_cells(EST_4G, "gsva", c("OX", "MYC", "PROLIF"), rd = s) %>%
      dplyr::select(-grey)))) %>%
  dplyr::mutate(grey = readout %in% SIG & term == "OX") %>%
  dplyr::inner_join(ADJ_COLS, by = c("cohort", "term")) %>%
  dplyr::inner_join(ROWS18 %>% dplyr::select(readout, y), by = "readout")
LIM_A <- .sym_lim(ADJ_T$est[!ADJ_T$grey])
ADJ_T <- ADJ_T %>% dplyr::mutate(
  fill_hex = ifelse(grey, GREY_CIRC, .scale_cols(est, LIM_A)),
  ink = ifelse(grey, "grey15", .ink(fill_hex)), label = sprintf("%.2f", est))
.assert("ED-h: 18 readouts x 14 raw cells and x 6 adjusted cells",
        nrow(RAW_T) == 252L && nrow(ADJ_T) == 108L, n = 360L)
# The right part takes the rest of the 178 mm, so the two parts sum exactly.
W_HR <- W_MM - W_HL
.stop_if(W_HR > max(ADJ_COLS$x) + 3.2, "ED-h: the right part is too narrow")
.edh_left <- function() {
  hdr <- RAW_COLS %>% dplyr::mutate(label = unname(EXPO[exposure]))
  coh <- RAW_COLS %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(x = mean(x), .groups = "drop")
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = .tag(RAW_T, "raw"),
                       ggplot2::aes(x = x, y = y, fill = rho),
                       width = 5.6, height = 3.6) +
    ggplot2::geom_text(data = RAW_T, ggplot2::aes(x = x, y = y, label = label,
                                                  colour = ink), size = SUB_MM) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x = x, y = -29.6, label = label),
                       angle = 90, hjust = 0, size = SUB_MM) +
    ggplot2::geom_text(data = coh, ggplot2::aes(x = x, y = -4.0, label = cohort),
                       size = TXT_MM, fontface = "bold") +
    ggplot2::geom_text(data = ROWS18, ggplot2::aes(
      x = 33, y = y, label = vapply(readout, .row_lab18, character(1))),
      hjust = 1, size = SUB_MM, parse = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = FILL_LOW, mid = "white", high = FILL_HIGH, midpoint = 0,
      limits = LIM_R, breaks = .sym_breaks(LIM_R[2]),
      name = "Spearman correlation, unadjusted",
      guide = ggplot2::guide_colourbar(direction = "horizontal",
                                       theme = ggplot2::theme(
                                         legend.title.position = "top",
                                         legend.key.width = grid::unit(30, "mm"),
                                         legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_HL), ylim = c(-H_H, 0),
                         expand = FALSE, clip = "off") +
    .theme_hm(c(0.98, 0.01))
}
.edh_right <- function(marked) {
  hdr <- ADJ_COLS %>% dplyr::mutate(label = unname(TERM_TITLE[term]))
  coh <- ADJ_COLS %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(x = mean(x), .groups = "drop")
  dots <- ADJ_T %>% dplyr::filter(!grey, .excl0(lo, hi))
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = .tag(ADJ_T[!ADJ_T$grey, ], "adj_col"),
                       ggplot2::aes(x = x, y = y, fill = est),
                       width = 6.0, height = 3.6) +
    ggplot2::geom_tile(data = ADJ_T[ADJ_T$grey, ], ggplot2::aes(x = x, y = y),
                       fill = GREY_CIRC, width = 6.0, height = 3.6) +
    ggplot2::geom_text(data = ADJ_T, ggplot2::aes(x = x, y = y, label = label,
                                                  colour = ink), size = SUB_MM) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x = x, y = -29.6, label = label),
                       angle = 90, hjust = 0, size = SUB_MM) +
    ggplot2::geom_text(data = coh, ggplot2::aes(x = x, y = -4.0, label = cohort),
                       size = TXT_MM, fontface = "bold")
  if (marked) {
    p <- p + ggplot2::geom_point(data = dots, ggplot2::aes(x = x + 2.55, y = y + 1.45),
                                 shape = 21, size = 0.6, fill = "black",
                                 colour = "white", stroke = 0.25)
  }
  p + ggplot2::scale_fill_gradient2(
    low = FILL_LOW, mid = "white", high = FILL_HIGH, midpoint = 0,
    limits = LIM_A, breaks = .sym_breaks(LIM_A[2]),
    name = "SD per SD, other two held fixed",
    guide = ggplot2::guide_colourbar(direction = "horizontal",
                                     theme = ggplot2::theme(
                                       legend.title.position = "top",
                                       legend.key.width = grid::unit(26, "mm"),
                                       legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_HR), ylim = c(-H_H, 0),
                         expand = FALSE, clip = "off") +
    .theme_hm(c(0.98, 0.01))
}
P_EDH_L <- .edh_left()
P_EDH_R <- list(a = .edh_right(FALSE), b = .edh_right(TRUE))
b_hl <- ggplot2::ggplot_build(P_EDH_L)
.assert("ED-h: every raw cell carries its stored correlation on its scale",
        identical(toupper(b_hl$data[[.layer_of(P_EDH_L, "raw")]]$fill),
                  .scale_cols(RAW_T$rho, LIM_R)), n = nrow(RAW_T))
for (v in names(P_EDH_R)) {
  b_hr <- ggplot2::ggplot_build(P_EDH_R[[v]])
  .assert(paste0("ED-h (", v, "): every adjusted cell carries its stored value"),
          identical(toupper(b_hr$data[[.layer_of(P_EDH_R[[v]], "adj_col")]]$fill),
                    .scale_cols(ADJ_T$est[!ADJ_T$grey], LIM_A)),
          n = sum(!ADJ_T$grey))
}
.assert("ED-h: the signatures' OXPHOS cells are grey, the circular cells",
        all(ADJ_T$readout[ADJ_T$grey] %in% SIG) &&
          all(ADJ_T$term[ADJ_T$grey] == "OX"), n = sum(ADJ_T$grey))
W_EDH <- c(W_HL, W_HR)
.check_text(list(P_EDH_L, P_EDH_R$a, P_EDH_R$b), "ED-h", W_PAGE)

# --- ED-i: intervals for all fifteen genes (9.2), E48's ED-d extended -------
KEYS <- tibble::tribble(
  ~key,          ~text,                                           ~col,    ~fill,   ~shape,
  "ox_with",     "OXPHOS, with the proliferation term",           COL_OX,  COL_OX,  21,
  "pr",          "Proliferation",                                 COL_PR,  COL_PR,  24,
  "myc_without", "MYC activity, without the proliferation term",  COL_MYC, "white", 22,
  "myc_with",    "MYC activity, with the proliferation term",     COL_MYC, COL_MYC, 22,
  "diff",        "Proliferation minus OXPHOS",                    "black", "black", 21)
K_TEXT  <- stats::setNames(KEYS$text, KEYS$key)
K_COL   <- stats::setNames(KEYS$col, KEYS$key)
K_FILL  <- stats::setNames(KEYS$fill, KEYS$key)
K_SHAPE <- stats::setNames(KEYS$shape, KEYS$key)
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
      plot.subtitle = ggplot2::element_text(size = MIN_PT),
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
      plot.subtitle = ggplot2::element_text(size = MIN_PT),
      plot.title = ggplot2::element_text(size = BASE_PT, face = "bold"),
      plot.margin = ggplot2::margin(2, 1, 2, 1))
}
.forest <- function(pf, cn, xlim, ylim, sep_y, y_cue, x_title, keys,
                    cue = c("lower expression", "higher expression"),
                    breaks = NULL, size = 1.0, facet = FALSE, title = NULL,
                    legend_ncol = NULL) {
  p <- ggplot2::ggplot(pf) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3)
  if (length(sep_y)) {
    p <- p + ggplot2::geom_hline(yintercept = sep_y, colour = "grey85",
                                 linewidth = 0.25)
  }
  if (!is.null(cn) && nrow(cn)) {
    p <- p + ggplot2::geom_segment(data = cn, ggplot2::aes(x = x, xend = xend,
                                                           y = y, yend = yend),
                                   colour = COL_MYC, linewidth = 0.3, alpha = 0.6,
                                   inherit.aes = FALSE)
  }
  w <- diff(xlim)
  a_in  <- c(-0.025 * w, 0.025 * w)
  a_out <- c(xlim[1] + 0.06 * w, xlim[2] - 0.06 * w)
  p <- p +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark,
                                       yend = y_mark, colour = key, alpha = void),
                          linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = est, y = y_mark, colour = key,
                                     fill = key, shape = key, alpha = void),
                        size = size, stroke = 0.45)
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
  p <- p +
    ggplot2::scale_x_continuous(limits = xlim, breaks = breaks,
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = K_COL, breaks = keys,
                                 labels = K_TEXT[keys], drop = FALSE, name = NULL) +
    ggplot2::scale_fill_manual(values = K_FILL, breaks = keys,
                               labels = K_TEXT[keys], drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = K_SHAPE, breaks = keys,
                                labels = K_TEXT[keys], drop = FALSE, name = NULL) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.3),
                                guide = "none") +
    ggplot2::labs(x = x_title, title = title) +
    .theme_forest()
  if (facet) {
    p <- p + ggplot2::facet_grid(. ~ panel, labeller = ggplot2::label_parsed) +
      ggplot2::theme(panel.spacing.x = grid::unit(5, "mm"))
  }
  if (!is.null(legend_ncol)) {
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE),
      fill = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE),
      shape = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE))
  }
  p
}
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
.coarse_breaks <- function(L) {
  k <- floor(L / 0.4 + 1e-9) * 0.4
  if (k <= 0) round(c(-L, 0, L), 10) else round(c(-k, 0, k), 10)
}
.check_points <- function(p, pf, label) {
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
                   PANEL = as.integer(dk$PANEL), colour = as.character(dk$colour))
  }))
  panel_f <- if ("panel" %in% names(pf)) as.integer(pf$panel) else rep(1L, nrow(pf))
  key_f <- paste(panel_f, round(pf$y_mark, 6))
  mp <- match(paste(pts$PANEL, round(pts$y, 6)), key_f)
  mc <- match(paste(cis$PANEL, round(cis$y, 6)), key_f)
  .assert(paste0(label, ": every point and interval drawn is one stored row"),
          nrow(pts) == nrow(pf) && !anyNA(mp) && !any(duplicated(mp)) &&
            nrow(cis) == nrow(pf) && !anyNA(mc) && !any(duplicated(mc)),
          n = nrow(pts), detail = paste(nrow(pts), "points,", nrow(cis),
                                        "intervals,", nrow(pf), "stored"))
  dev <- max(.dev(pts$x, pf$est[mp]), .dev(cis$x, pf$lo[mc]),
             .dev(cis$xend, pf$hi[mc]))
  .assert(paste0(label, ": drawn positions equal the stored values"),
          dev < 1e-12, n = nrow(pts) + 2L * nrow(cis), max_dev = dev)
  open <- grepl("without", as.character(pf$key[mp]))
  .assert(paste0(label, ": open = without the proliferation term, filled = with it"),
          all(pts$fill[open] == "white") && all(pts$fill[!open] != "white"),
          n = nrow(pts))
  v <- pf$void[mp]
  .assert(paste0(label, ": no void cell is drawn as read"),
          all(pts$alpha[v] < 1) && all(is.na(pts$alpha[!v]) | pts$alpha[!v] == 1),
          n = sum(v))
  data_cols <- setdiff(unique(c(pts$colour, cis$colour)), "black")
  .assert(paste0(label, ": every data colour drawn is from the validated palette"),
          all(data_cols %in% VALIDATED_PALETTE), n = length(data_cols))
}
.check_axis <- function(pf, xlim, label) {
  v <- c(pf$lo, pf$hi)
  .assert(paste0(label, ": the axis is symmetric, contains every interval and is ",
                 "rounded outward to 0.1"),
          all(v >= xlim[1] & v <= xlim[2]) && abs(xlim[1] + xlim[2]) < 1e-12 &&
            abs(xlim[2] * 10 - round(xlim[2] * 10)) < 1e-9, n = length(v),
          detail = sprintf("limits %.2f to %.2f", xlim[1], xlim[2]))
}
.long <- function(E, keys, rows, id_cols, keep = character(0)) {
  kc <- unique(c(id_cols, keep))
  parts <- list(
    ox_with     = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         est = b_OX_m3, lo = ox_lo_m3,
                                         hi = ox_hi_m3, void = void),
    pr          = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         est = b_PR, lo = pr_lo, hi = pr_hi,
                                         void = void),
    myc_without = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         est = b_MYC_m1, lo = myc_lo_m1,
                                         hi = myc_hi_m1, void = FALSE),
    myc_with    = E %>% dplyr::transmute(dplyr::across(dplyr::all_of(kc)),
                                         est = b_MYC_m3, lo = myc_lo_m3,
                                         hi = myc_hi_m3, void = void))
  dplyr::bind_rows(lapply(names(keys), function(k)
    parts[[k]] %>% dplyr::mutate(key = k))) %>%
    dplyr::left_join(rows, by = id_cols) %>%
    dplyr::mutate(y_mark = y_row + unname(keys[key]),
                  key = factor(key, levels = names(keys)))
}

ROWS_I <- tidyr::expand_grid(readout = FIFTEEN, block = BLOCK_ORDER) %>%
  dplyr::mutate(g = match(readout, FIFTEEN),
                grp = match(ROWS15$group[match(readout, ROWS15$readout)],
                            names(GROUPS)),
                j = match(block, BLOCK_ORDER),
                top = -0.8 - (g - 1) * 4.6 - (grp - 1) * 1.2, y_row = top - j)
YL_I <- c(min(ROWS_I$y_row) - 1.9, 0.3)
YCUE_I <- YL_I[1] + 0.45
SEP_I <- ROWS_I %>% dplyr::distinct(g, top) %>% dplyr::filter(g > 1L) %>%
  dplyr::pull(top) + 0.3
OFF_I <- c(ox_with = 0.25, pr = 0, myc_with = -0.25)
.edi <- function(est) {
  E <- ALLX %>% dplyr::filter(estimator == est, readout %in% FIFTEEN) %>%
    dplyr::mutate(block = paste(cohort, instrument))
  .stop_if(nrow(E) == 60L, "ED-i: missing rows on ", est)
  pf <- .long(E, OFF_I, ROWS_I %>% dplyr::select(readout, block, y_row),
              c("readout", "block"))
  pd <- E %>% dplyr::transmute(readout, block, key = "diff", est = diff,
                               lo = diff_lo, hi = diff_hi, void = FALSE) %>%
    dplyr::left_join(ROWS_I %>% dplyr::select(readout, block, y_row),
                     by = c("readout", "block")) %>%
    dplyr::mutate(y_mark = y_row, key = factor(key, levels = "diff"))
  list(coef = pf, diff = pd, E = E)
}
I_REF <- .edi(EST_REF)
I_CLN <- .edi(EST_4G)
for (nm in c("ref", "clean")) {
  E <- if (nm == "ref") I_REF$E else I_CLN$E
  e48 <- G12 %>% dplyr::filter(estimator == unique(E$estimator))
  chk <- E %>% dplyr::filter(readout %in% TWELVE) %>%
    dplyr::inner_join(e48, by = c("cohort", "instrument", "estimator", "readout"),
                      suffix = c("", "_48"))
  dv <- max(.dev(chk$b_OX_m3, chk$b_OX_m3_48), .dev(chk$b_PR, chk$b_PR_48),
            .dev(chk$b_MYC_m3, chk$b_MYC_m3_48), .dev(chk$diff, chk$diff_48))
  .assert(paste0("ED-i (", nm, "): the twelve equal E48's saved values"),
          nrow(chk) == 48L && dv == 0, n = 4L * nrow(chk), max_dev = dv)
}
XL_I  <- .sym_lim(c(I_REF$coef$lo, I_REF$coef$hi, I_CLN$coef$lo, I_CLN$coef$hi))
XL_ID <- .sym_lim(c(I_REF$diff$lo, I_REF$diff$hi, I_CLN$diff$lo, I_CLN$diff$hi))
for (nm in c("ref", "clean")) {
  dd <- if (nm == "ref") I_REF else I_CLN
  .check_axis(dd$coef, XL_I, paste0("ED-i, ", nm, ", coefficients"))
  .check_axis(dd$diff, XL_ID, paste0("ED-i, ", nm, ", difference"))
}
.fori <- function(pf, title, xlim, keys, x_title, diff_panel = FALSE) {
  .forest(pf, NULL, xlim, YL_I, SEP_I, YCUE_I, x_title, keys,
          cue = if (diff_panel) NULL else c("lower expression", "higher expression"),
          breaks = if (diff_panel) .coarse_breaks(xlim[2]) else .sym_breaks(xlim[2]),
          size = 0.9, title = title)
}
FOR_I <- list(
  .fori(I_REF$coef, "Felsher, reference score", XL_I, names(OFF_I),
        "SD per SD of each score"),
  .fori(I_REF$diff, "Difference", XL_ID, "diff", "Proliferation\nminus OXPHOS",
        diff_panel = TRUE),
  .fori(I_CLN$coef, "MYC_UP.V1_UP, least-entangled", XL_I, names(OFF_I),
        "SD per SD of each score"),
  .fori(I_CLN$diff, "Difference", XL_ID, "diff", "Proliferation\nminus OXPHOS",
        diff_panel = TRUE))
GENE_POS_I <- ROWS_I %>% dplyr::group_by(readout) %>%
  dplyr::summarise(y = mean(y_row), .groups = "drop")
GRP_POS_I <- ROWS_I %>% dplyr::mutate(group = ROWS15$group[match(readout, ROWS15$readout)]) %>%
  dplyr::group_by(group) %>% dplyr::summarise(y = max(top) + 0.25, .groups = "drop")
LBL_I <- .text_panel(YL_I, list(
  .lyr(0.02, GRP_POS_I$y, GRP_POS_I$group, hjust = 0, size = SUB_MM, face = "bold"),
  .lyr(0.02, GENE_POS_I$y, vapply(GENE_POS_I$readout, .gene_lab, character(1)),
       hjust = 0, size = SUB_MM, parse = TRUE),
  .lyr(0.48, ROWS_I$y_row, unname(BLOCK_LAB[ROWS_I$block]), hjust = 0,
       size = SUB_MM, parse = TRUE)))
PARTS_EDI <- c(list(LBL_I), FOR_I)
W_EDI <- c(1.55, 1.85, 0.90, 1.85, 0.90)
FIG_EDI <- patchwork::wrap_plots(PARTS_EDI, ncol = 5,
                                 widths = W_EDI / sum(W_EDI) * W_PAGE) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
H_EDI <- 9.40
.check_points(FOR_I[[1]], I_REF$coef, "ED-i, reference, coefficients")
.check_points(FOR_I[[2]], I_REF$diff, "ED-i, reference, difference")
.check_points(FOR_I[[3]], I_CLN$coef, "ED-i, least-entangled, coefficients")
.check_points(FOR_I[[4]], I_CLN$diff, "ED-i, least-entangled, difference")
.check_text(PARTS_EDI, "ED-i", W_PAGE)

# --- ED-j: the MYC ladder for the signatures and the effectors (9.2) --------
ROWS_J <- tidyr::expand_grid(estimator = EST_ALL, block = BLOCK_ORDER) %>%
  dplyr::mutate(s = match(estimator, EST_ALL), j = match(block, BLOCK_ORDER),
                top = -(s - 1) * 5.0, y_head = top - 0.30, y_row = top - 0.45 - j)
YL_J <- c(-21.50, 0.30)
SEP_J <- -(1:3) * 5.0 + 0.25
YCUE_J <- -21.05
OFF_J <- c(myc_without = 0.22, myc_with = 0, pr = -0.22)
.edj <- function(readouts, labs) {
  E <- ALLX %>% dplyr::filter(readout %in% readouts) %>%
    dplyr::mutate(block = paste(cohort, instrument))
  pf <- .long(E, OFF_J, ROWS_J %>% dplyr::select(estimator, block, y_row),
              c("estimator", "block"), keep = "readout") %>%
    dplyr::mutate(panel = factor(unname(labs[readout]), levels = unname(labs)))
  .stop_if(nrow(pf) == length(readouts) * 16L * 3L, "ED-j frame has the wrong size")
  pa <- pf[pf$key == "myc_without", ]
  pb <- pf[pf$key == "myc_with", ]
  m  <- match(paste(pa$readout, pa$estimator, pa$block),
              paste(pb$readout, pb$estimator, pb$block))
  cn <- tibble::tibble(x = pa$est, xend = pb$est[m], y = pa$y_mark,
                       yend = pb$y_mark[m], panel = pa$panel)
  list(pf = pf, cn = cn)
}
LAB_SIG <- c(S6 = '"Signature (S6)"', S5 = '"S5, without BID"',
             S7 = '"S7, with BMF"')
LAB_EFF <- stats::setNames(vapply(EFF, .gene_lab, character(1)), EFF)
J_SIG <- .edj(SIG, LAB_SIG)
J_EFF <- .edj(EFF, LAB_EFF)
XL_J <- .sym_lim(c(J_SIG$pf$lo, J_SIG$pf$hi, J_EFF$pf$lo, J_EFF$pf$hi))
.check_axis(J_SIG$pf, XL_J, "ED-j, signatures")
.check_axis(J_EFF$pf, XL_J, "ED-j, effectors")
FOR_JS <- .forest(J_SIG$pf, J_SIG$cn, XL_J, YL_J, SEP_J, YCUE_J,
                  "Association with expression (SD per SD of each score)",
                  names(OFF_J), breaks = .sym_breaks(XL_J[2]), size = 1.0,
                  facet = TRUE)
FOR_JE <- .forest(J_EFF$pf, J_EFF$cn, XL_J, YL_J, SEP_J, YCUE_J,
                  "Association with expression (SD per SD of each score)",
                  names(OFF_J), breaks = .sym_breaks(XL_J[2]), size = 1.0,
                  facet = TRUE)
HEADS_J <- ROWS_J %>% dplyr::distinct(estimator, y_head)
.lbl_j <- function() .text_panel(YL_J, list(
  .lyr(0.02, HEADS_J$y_head, unname(EST_NAME[HEADS_J$estimator]), hjust = 0,
       size = SUB_MM, face = "bold"),
  .lyr(0.02, HEADS_J$y_head[HEADS_J$estimator == "entangled"] - 0.62,
       "with proliferation: not identifiable, faded", hjust = 0, size = SUB_MM,
       colour = "grey35"),
  .lyr(0.10, ROWS_J$y_row, unname(BLOCK_LAB[ROWS_J$block]), hjust = 0,
       size = SUB_MM, parse = TRUE)))
PARTS_EDJ <- list(.lbl_j(), FOR_JS, .lbl_j(), FOR_JE)
W_EDJ <- c(1.55, 5.45)
FIG_EDJ <- patchwork::wrap_plots(PARTS_EDJ, ncol = 2,
                                 widths = W_EDJ / sum(W_EDJ) * W_PAGE) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
H_EDJ <- 8.60
.check_points(FOR_JS, J_SIG$pf, "ED-j, signatures")
.check_points(FOR_JE, J_EFF$pf, "ED-j, effectors")
.check_text(PARTS_EDJ, "ED-j", W_PAGE)

# --- ED-k: the Fig. 1i-style human panel (6.8) ------------------------------
XL_K <- range(FIG1I_DEC$x) * 1.15
YL_K <- range(FIG1I_DEC$y) * 1.30
.edk <- function(coh) {
  s <- FIG1I[FIG1I$cohort == coh, ]
  d <- FIG1I_DEC[FIG1I_DEC$cohort == coh, ]
  ln <- tibble::tibble(x = XL_K, y = s$slope * XL_K)
  txt <- sprintf("slope %.2f (95%% CI %.2f to %.2f)\npartial correlation %.2f",
                 s$slope, s$slope_lo, s$slope_hi, s$partial_r)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_line(data = ln, ggplot2::aes(x = x, y = y), colour = "black",
                       linewidth = 0.4) +
    ggplot2::geom_point(size = 1.4, colour = "black") +
    ggplot2::scale_x_continuous(limits = XL_K) +
    ggplot2::scale_y_continuous(limits = YL_K) +
    ggplot2::labs(title = paste0(coh, " (n = ", format(s$n, big.mark = ","), ")"),
                  subtitle = txt,
                  x = "OXPHOS priority, adjusted (SD)",
                  y = "PUMA relative to BCL-xL, adjusted (SD)") +
    ggplot2::theme_classic(base_size = BASE_PT) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = MIN_PT, colour = "black"),
                   axis.title = ggplot2::element_text(size = BASE_PT - 0.5),
                   plot.title = ggplot2::element_text(size = BASE_PT, face = "bold",
                                                      hjust = 0.5),
                   plot.caption = ggplot2::element_text(size = MIN_PT),
                   plot.subtitle = ggplot2::element_text(size = MIN_PT, hjust = 0.5,
                                                         lineheight = 0.95),
                   strip.text = ggplot2::element_text(size = BASE_PT),
                   legend.text = ggplot2::element_text(size = MIN_PT),
                   legend.title = ggplot2::element_text(size = MIN_PT),
                   plot.margin = ggplot2::margin(3, 5, 3, 3))
}
PARTS_EDK <- lapply(stats::setNames(COHORTS, COHORTS), .edk)
FIG_EDK <- patchwork::wrap_plots(PARTS_EDK, ncol = 2,
                                 widths = rep(W_PAGE / 2, 2))
H_EDK <- 2.85
for (coh in COHORTS) {
  b <- ggplot2::ggplot_build(PARTS_EDK[[coh]])
  geo <- vapply(PARTS_EDK[[coh]]$layers, function(l) class(l$geom)[1],
                character(1))
  pt <- b$data[[which(geo == "GeomPoint")]]
  d  <- FIG1I_DEC[FIG1I_DEC$cohort == coh, ]
  .assert(paste0("ED-k, ", coh, ": the ten points are the stored decile means"),
          nrow(pt) == 10L && .dev(pt$x, d$x) < 1e-12 && .dev(pt$y, d$y) < 1e-12,
          n = 10L)
  ln <- b$data[[which(geo == "GeomLine")]]
  sl <- diff(ln$y) / diff(ln$x)
  .assert(paste0("ED-k, ", coh, ": the drawn line has the stored slope"),
          abs(sl - FIG1I$slope[FIG1I$cohort == coh]) < 1e-10, n = 1L)
  .assert(paste0("ED-k, ", coh, ": every decile mean lies inside the axes"),
          all(d$x >= XL_K[1] & d$x <= XL_K[2] & d$y >= YL_K[1] & d$y <= YL_K[2]),
          n = 10L)
}
.check_text(PARTS_EDK, "ED-k", W_PAGE)

# --- save every figure: outputs/figures, and docs/figures (tracked) ----------
.ensure_dir(DIR_FIGURES)
.ensure_dir(DIR_DOC_FIGURES)
.save_fig <- function(fig, stem, h_in) {
  for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
    base <- file.path(dir, if (identical(dir, DIR_DOC_FIGURES))
      paste0(FIG_DATE, "_", stem) else stem)
    ggplot2::ggsave(paste0(base, ".png"), fig, width = W_PAGE, height = h_in,
                    units = "in", dpi = 300, bg = "white")
    ggplot2::ggsave(paste0(base, ".pdf"), fig, width = W_PAGE, height = h_in,
                    units = "in", device = grDevices::cairo_pdf, bg = "white")
  }
}
FIG_EDH <- list(
  a = patchwork::wrap_plots(list(P_EDH_L, P_EDH_R$a), ncol = 2,
                            widths = W_EDH / sum(W_EDH) * W_PAGE),
  b = patchwork::wrap_plots(list(P_EDH_L, P_EDH_R$b), ncol = 2,
                            widths = W_EDH / sum(W_EDH) * W_PAGE))
.assert("ED-h's two parts sum to 178 mm", abs(sum(W_EDH) - W_MM) < 1e-9, n = 1L,
        detail = sprintf("%.1f mm", sum(W_EDH)))
FIGS <- list(
  fig4g_a = list(HM$fig4g_a$plot, HM$fig4g_a$height_mm / 25.4),
  fig4g_b = list(HM$fig4g_b$plot, HM$fig4g_b$height_mm / 25.4),
  ed_e_a  = list(HM$ed_e_a$plot, HM$ed_e_a$height_mm / 25.4),
  ed_e_b  = list(HM$ed_e_b$plot, HM$ed_e_b$height_mm / 25.4),
  ed_f_a  = list(HM$ed_f_a$plot, HM$ed_f_a$height_mm / 25.4),
  ed_f_b  = list(HM$ed_f_b$plot, HM$ed_f_b$height_mm / 25.4),
  ed_g    = list(P_EDG, H_G / 25.4),
  ed_h_a  = list(FIG_EDH$a, H_H / 25.4),
  ed_h_b  = list(FIG_EDH$b, H_H / 25.4),
  ed_i    = list(FIG_EDI, H_EDI),
  ed_j    = list(FIG_EDJ, H_EDJ),
  ed_k    = list(FIG_EDK, H_EDK))
.assert("twelve figures, each with a declared file stem",
        setequal(names(FIGS), names(FIG)), n = length(FIGS))
for (nm in names(FIGS)) .save_fig(FIGS[[nm]][[1]], FIG[[nm]], FIGS[[nm]][[2]])
message("   twelve figures, 178 mm wide, to outputs/figures and docs/figures ",
        "(tracked)")

# =============================================================================
# 9. Legend drafts - DESCRIPTIVE ONLY (8.4, 9.1)
# =============================================================================
N_TXT <- "TCGA-BRCA (n = 938) and SCAN-B (n = 3,143)"
SECT_NOTE <- paste0(SECT, " BH3-only proteins whose activity is switched after ",
                    "translation (BAD: off by phosphorylation; BID: on by ",
                    "cleavage to tBID).")
BIND_NOTE <- paste0(
  "Left, binding of BCL-xL and MCL1 by BH3-only proteins, from Chen et al. ",
  "(2005): filled, binds; open, does not; blank, not yet verified. Signature ",
  "sign: the direction in which each of the six signature genes enters the ",
  "signature.")
MODEL_4G <- paste0(
  "Each row is one linear model of the gene's expression on OXPHOS, MYC ",
  "activity and proliferation, with intrinsic subtype, and with tumour purity ",
  "and leukocyte fraction in TCGA-BRCA; each coefficient is in SD per SD with ",
  "the other two held fixed, and every variable is rank-normalised within each ",
  "cohort. MYC activity is scored on MYC_UP.V1_UP, a score with few ",
  "proliferation genes; the other MYC activity scores are shown in Extended ",
  "Data.")
STRIP_NOTE <- paste0(
  "Bottom, the six-gene signature: its coefficient and 95% confidence ",
  "interval in each column; its association with OXPHOS, in grey, defines the ",
  "signature and is shown, not interpreted.")
LEGEND <- list(
  fig4g_a = paste0(
    "Associations of OXPHOS level, MYC activity and proliferation with the ",
    "expression of fifteen BCL-2 family genes, in the groups of Fig. 1g, in ",
    N_TXT, ". ", MODEL_4G, " ", BIND_NOTE, " ", SECT_NOTE, " ", STRIP_NOTE,
    " OXPHOS priority is shown in Extended Data."),
  fig4g_b = "As the unmarked version, with a dot on each cell whose 95% confidence interval excludes zero.",
  ed_e = paste0("As Fig. 4g, with OXPHOS scored as its priority within the ",
                "mitochondrial programme (mitoPPS)."),
  ed_f = paste0("As Fig. 4g, from models using the reference MYC activity score ",
                "(Felsher), which correlates with proliferation at 0.61 to 0.70."),
  ed_g = paste0(
    "The six-gene signature (S6) and its two variants, S5 (without BID) and S7 ",
    "(with BMF), on the MYC activity score with few proliferation genes and on ",
    "the reference score, with OXPHOS scored as level and as priority, in ",
    N_TXT, ". Each cell gives the coefficient (SD per SD) and its 95% ",
    "confidence interval. Grey cells contain the OXPHOS coefficient, which ",
    "defines the signature: OXPHOS itself, and proliferation minus OXPHOS."),
  ed_h = paste0(
    "Left, unadjusted Spearman correlations of the signature, its variants and ",
    "fifteen BCL-2 family genes with OXPHOS (level and priority), four MYC ",
    "activity scores and proliferation, in ", N_TXT, ". Right, the adjusted ",
    "coefficients of Fig. 4g. Grey, the signatures' association with OXPHOS, ",
    "which defines them. The marked version adds a dot where the 95% ",
    "confidence interval of an adjusted coefficient excludes zero."),
  ed_i = paste0(
    "Fifteen BCL-2 family genes, in the groups of Fig. 1g: coefficients of ",
    "OXPHOS (circles), proliferation (triangles) and MYC activity (squares), ",
    "each with the other two held fixed, and proliferation minus OXPHOS, on the ",
    "reference MYC activity score (left) and the score with few proliferation ",
    "genes (right), in ", N_TXT, ". ", SECT_NOTE, " No correction for multiple ",
    "genes. Bars, 95% confidence intervals."),
  ed_j = paste0(
    "MYC activity coefficients for the signature, its two variants and BAX, ",
    "BAK1 and BOK on four MYC activity scores, ordered by the percentage of ",
    "their genes in proliferation gene sets, without (open) and with (filled) ",
    "the proliferation term, with the proliferation coefficient (triangles), in ",
    N_TXT, ". Faded points are cells in which the score and proliferation ",
    "correlate at 0.80 or more; they are shown and not interpreted. Bars, 95% ",
    "confidence intervals."),
  ed_k = paste0(
    "BBC3 (PUMA) relative to BCL2L1 (BCL-xL) against OXPHOS priority, each ",
    "adjusted for the reference MYC activity score, proliferation and intrinsic ",
    "subtype (and purity and leukocyte fraction in TCGA-BRCA), in ", N_TXT,
    ". Points, means of the ten deciles of adjusted OXPHOS priority; line, the ",
    "fitted slope, with its 95% interval from 2,000 paired bootstrap resamples ",
    "of tumours."))
BANNED <- c("independently of MYC", "MYC expression", "comparable", "invariant",
            "lost after adjustment", "primed", "priming",
            "increases reliance on BCL-xL rather than overriding",
            "which BCL-xL can absorb and MCL1 cannot")
.assert("no legend draft uses a word the briefs keep out (4.9; 4.8)",
        !any(vapply(BANNED, function(w) any(grepl(w, unlist(LEGEND),
                                                  ignore.case = TRUE)),
                    logical(1))), n = length(LEGEND))
.assert("no legend draft calls a MYC activity score a signature",
        !any(grepl("MYC[^.;]{0,30}signature", unlist(LEGEND))), n = length(LEGEND))
.assert("the legend carries the author's section-mark note verbatim (9.1)",
        grepl(SECT_NOTE, LEGEND$fig4g_a, fixed = TRUE), n = 1L)

# =============================================================================
# 10. Save, per R rules 4 and 5
# =============================================================================
message("\n10. saving")

ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("   ", nrow(ASSERTIONS), " assertions, all passed")

FRAME_OUT <- FRAME49 %>%
  dplyr::select(sample_id, cohort, in_set, PAM50, purity, leuko, PROLIF,
                OX_gsva, OX_mitopps, MYC_clean, MYC_M_a, MYC_entangled,
                MYC_regulon, dplyr::all_of(SIG), dplyr::all_of(TWELVE),
                dplyr::all_of(EFF))

readr::write_csv(EXT_SIG,   file.path(DIR_TABLES, "E49_signature_coefficients.csv"))
readr::write_csv(EXT_EFF,   file.path(DIR_TABLES, "E49_effector_coefficients.csv"))
readr::write_csv(ALLX,      file.path(DIR_TABLES, "E49_all_readouts_m3.csv"))
readr::write_csv(RAWC,      file.path(DIR_TABLES, "E49_raw_correlations.csv"))
readr::write_csv(PCT,       file.path(DIR_TABLES, "E49_expression_percentiles.csv"))
readr::write_csv(PLACEMENT, file.path(DIR_TABLES, "E49_placement_rule.csv"))
readr::write_csv(FIG1I,     file.path(DIR_TABLES, "E49_fig1i_style.csv"))
readr::write_csv(FIG1I_DEC, file.path(DIR_TABLES, "E49_fig1i_style_deciles.csv"))
readr::write_csv(ASSERTIONS, file.path(DIR_TABLES, "E49_assertions.csv"))

saveRDS(list(
  signature_coefficients = EXT_SIG,
  genes12                = G12,
  effector_coefficients  = EXT_EFF,
  all_readouts_m3        = ALLX,
  effector_coefficient_matrix = dplyr::bind_rows(EFF_COEF),
  effector_vcov          = EFF_VCOV,
  effector_vif           = dplyr::bind_rows(EFF_VIF),
  fit_digests            = dplyr::bind_rows(EFF_DIG),
  gate_a1                = GA1,
  gate_a3                = GA3,
  gate_b1                = GB1,
  gate_b2                = GB2,
  placement_rule         = list(inputs = PLACEMENT, keep_myc = KEEP_MYC,
                                fig4g = SPEC_4G),
  raw_correlations       = RAWC,
  gate_c1                = GC1,
  expression_percentiles = PCT,
  fig1i                  = list(summary = FIG1I, deciles = FIG1I_DEC,
                                residuals = FIG1I_RES),
  binding_code           = BIND,
  frame                  = FRAME_OUT,
  displays               = list(
    heatmaps = lapply(HM, function(h) h[c("cells", "strip", "dots", "bind",
                                          "signs", "lim", "cols", "model",
                                          "instrument", "terms", "marked",
                                          "height_mm")]),
    ed_g = EDG, ed_h_raw = RAW_T, ed_h_adjusted = ADJ_T,
    ed_i_reference = I_REF[c("coef", "diff")],
    ed_i_least_entangled = I_CLN[c("coef", "diff")],
    ed_j_signatures = J_SIG, ed_j_effectors = J_EFF,
    ed_k = list(deciles = FIG1I_DEC, summary = FIG1I),
    axes = list(ed_g = LIM_G, ed_h = list(raw = LIM_R, adjusted = LIM_A),
                ed_i = list(XL_I, XL_ID), ed_j = XL_J, ed_k = list(XL_K, YL_K)),
    grey_off_scale_dE = DE_GREY),
  legend_draft           = LEGEND,
  assertions             = ASSERTIONS,
  # FILE NAMES, NOT PATHS, so this element digests the same in a redirected
  # dry run and in the author's run. Figures are not digested.
  figure_files = list(
    outputs = paste0(rep(unname(FIG), each = 2), c(".png", ".pdf")),
    tracked_in_docs_figures = paste0(FIG_DATE, "_", rep(unname(FIG), each = 2),
                                     c(".png", ".pdf"))),
  spec = list(
    declaration = paste0(
      "docs/2026-11-03_E49_signature_declaration.md, committed with ",
      "data/annotations/ at 7bd7865, before this script. Brief: ",
      "docs/2026-11-03_signature_and_fig4g_brief.md, eb45a9b."),
    posture = paste0("EXPLORATORY and POST HOC. No reading rule, no label. ",
                     "Part A fits nothing; Part B fits the three effectors ",
                     "only. Main effects only; no interaction."),
    placement_rule = paste0(
      "Amendment item 4: if S6's m3 MYC coefficient on the CollecTRI regulon ",
      "has the same sign as on the least-entangled score in all four blocks, ",
      "4g keeps the MYC column (4g's model, the least-entangled m3 fit); ",
      "otherwise 4g shows OXPHOS and proliferation only, from the reference ",
      "model, and MYC activity appears only in Extended Data. Its outcome was ",
      "already implied by E47 result 3.2 (declaration 2.4)."),
    circularity = paste0(
      "Section 8.2: the signatures' OXPHOS coefficient is drawn grey and never ",
      "read; their difference contains it and is flagged likewise."),
    amendment = paste0(
      "Declaration section 3, verbatim there: 4g has three exposure blocks ",
      "from one model per row on the least-entangled score; the strip uses ",
      "the same model; the reference-score heatmap goes to Extended Data; ",
      "the placement rule; the legend names the score; brief 8's two ",
      "sentences (declaration 14.1)."),
    palette = paste0("Heatmaps: E46's RdBu poles, #2166AC and #B2182B, white ",
                     "at 0 (checker: ALL CHECKS PASS). Grey #A6A6A6 for the ",
                     "circular cells, least dE to the scale >= 15. Forests: ",
                     "#009E73, #D55E00, #0072B2 (E48's validated palette)."),
    snapshot = paste0("E47's and E48's saved objects; for Part B, ",
                      "data/from_validation/{tcga_brca,scanb}_linear.rds, ",
                      "log2(x + 1), SCAN-B through scanb_pheno.rds$symbol_map. ",
                      "The joint snapshot is not read."),
    not_digested = paste0("$fits_effectors holds 96 fitted lm objects, which ",
                          "are not digest-stable; $fit_digests is the stable ",
                          "handle."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  fits_effectors = FITS_EFF,
  built          = Sys.time()), PATH_E49)

message("   ", PATH_E49)
message("\nE49 done. Fig. 4g: ",
        if (KEEP_MYC) "three blocks from 4g's model (placement rule kept MYC)." else
          "OXPHOS and proliferation from the reference model (placement rule).",
        "\nNo reading rule: everything above is reported, nothing labelled.\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E49)

  # the placement rule (8.3)
  x$placement_rule$inputs %>% as.data.frame()

  # the signatures at m3: proliferation and the difference (6.1)
  x$signature_coefficients %>%
    dplyr::filter(estimator %in% c("clean", "M_a")) %>%
    dplyr::transmute(readout, estimator, cohort, instrument,
                     prolif = .fmt3(b_PR, pr_lo, pr_hi),
                     difference = .fmt3(diff, diff_lo, diff_hi)) %>%
    as.data.frame()

  # the effectors at m3, both scores
  x$effector_coefficients %>%
    dplyr::filter(estimator %in% c("clean", "M_a")) %>%
    dplyr::transmute(readout, estimator, cohort, instrument,
                     ox = .fmt3(b_OX_m3, ox_lo_m3, ox_hi_m3),
                     myc = .fmt3(b_MYC_m3, myc_lo_m3, myc_hi_m3),
                     prolif = .fmt3(b_PR, pr_lo, pr_hi)) %>%
    as.data.frame()

  x$expression_percentiles %>% as.data.frame()
  x$fig1i$summary %>% as.data.frame()
  x$assertions %>% as.data.frame()
  cat(x$legend_draft$fig4g_a, "\n")
}
