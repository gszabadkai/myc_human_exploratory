# E47_mutual_adjustment.R
# =============================================================================
# MYC AND OXPHOS ON THE BCL-2-FAMILY SIGNATURE, EACH WITH THE OTHER HELD FIXED.
#
# EXPLORATORY and POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-22_E47_mutual_adjustment_declaration.md, committed
# ALONE at dcb58a7, with addendum 4a.7 at 3b5e221. Both were committed BEFORE
# this script existed. Every choice below is the declaration's, and section
# numbers cited are its sections. Read sections 0, 3, 4a, 7 and 10 first.
#
# =============================================================================
# THE QUESTION (section 0)
# =============================================================================
# Fig. 4g (E46) shows a MYC increment on the signature at both OXPHOS levels,
# most visibly on mitoPPS. Its quadrants are median splits, and a median split
# does not hold the other axis fixed. This script fits, on Fig. 4g's own scores,
# analysis sets and covariates,
#
#     readout ~ OX + MYC + spine            OLS, t-based 95% CI, as E34
#
# so that beta_MYC is MYC with OXPHOS held fixed and beta_OX is OXPHOS with MYC
# held fixed. MAIN EFFECTS ONLY: no product term (section 3).
#
# THE VERDICT READS beta_MYC|OX ON THE SIGNATURE (S6) ONLY (section 7).
# beta_OX|MYC is reported in every cell and NEVER read. The signature's genes
# and signs were selected on OXPHOS in these same cohorts, so beta_OX|MYC is
# near-guaranteed by construction (E34 1.2).
#
# =============================================================================
# TERMS (the declaration's "Terms used from here on")
# =============================================================================
#   S6              the six-gene BCL-2-family signature: "the configuration
#                   composite" of E33/E34, config_comp in E34's frames
#   S5, S7          the rule-based variants (4a.6): S5 drops BID; S7 adds BMF
#                   with sign -1
#   signature sign  the +1/-1 with which a gene enters a signature (4a.4)
#   estimator       a score of MYC activity. Never called a "signature" here
#
# S5 IS BIASED AGAINST MYC BY CONSTRUCTION (4a.7). Its rule keeps a gene only
# if the gene's OXPHOS correlation beats its MYC correlation by 0.20. An S5-S6
# disagreement therefore means "MYC's contribution is carried by BID", never
# that the contribution is spurious.
#
# =============================================================================
# SNAPSHOT AND SCALE
# =============================================================================
# Inputs are data/from_validation/ and E02's scores, as for E34. The joint
# snapshot is NOT read.
#   - OXPHOS: GSVA on VST, and mitoPPS on linear counts. READ, never rescored.
#   - MYC estimators: E02's scores. READ, never rescored.
#   - Genes: log2(linear DESeq2-normalised + 1). SCAN-B is resolved through
#     scanb_pheno.rds$symbol_map (trap 7).
# Every readout and every exposure is Blom-scored ONCE PER COHORT, over the
# whole cohort, as E34 scored its readouts and PROLIF. The fits then subset
# and never rescore.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits no product term and tests no interaction (section 3).
#   - It reads no verdict from beta_OX|MYC. It never prints beta_MYC/beta_OX
#     without saying that beta_OX is inflated by construction.
#   - It applies section 7's rule to S6 for the verdict, and to S5 and S7
#     DESCRIPTIVELY. It applies no rule to the arms or the twelve genes.
#   - It fits neither M_b__FULL nor M_b__PROLIFSTRIP, which each contain 5 of
#     the 12, and not log2 MYC mRNA either (trap 4).
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

message("\nE47: MYC and OXPHOS on the BCL-2-family signature, each held fixed\n",
        strrep("=", 78))
message("EXPLORATORY, POST-HOC. Declaration dcb58a7 + addendum 3b5e221.")
message("The verdict reads beta_MYC|OX on S6. beta_OX|MYC is never read.")

PATH_E34  <- file.path(DIR_RESULTS, "e34_quadrant_configuration.rds")
PATH_E10  <- file.path(DIR_RESULTS, "machinery_and_priming.rds")
PATH_NEW  <- file.path(DIR_RESULTS, "new_set_scores.rds")
PATH_SCS  <- file.path(DIR_RESULTS, "scanb_scores.rds")
PATH_SETS <- file.path(DIR_RESULTS, "set_definitions.rds")
PATH_E47  <- file.path(DIR_RESULTS, "e47_mutual_adjustment.rds")

DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_A <- "E47_display_A_genes"
FIG_B <- "E47_display_B_estimator_ladder"
FIG_C <- "E47_display_C_bcl2l1_bbc3"
FIG_DATE <- "2026-10-23"

# =============================================================================
# 0. CONSTANTS - all from the declaration, none chosen here
# =============================================================================
COHORTS <- c("TCGA", "SCAN-B")
INSTR   <- c(gsva = "GSVA", mitopps = "mitoPPS")
OX_ARM  <- "OXPHOS subunits"
PROLIF_SET <- "PROLIF_DISJOINT"

# Section 5. The primary estimator is Fig. 4g's axis; the ladder is mandatory
# (N9). Entanglement is the declaration's own figure for each.
EST <- c(clean     = "MYC_UP.V1_UP__MITOSTRIP",
         M_a       = "FELSHER__MITOSTRIP",
         entangled = "HALLMARK_MYC_TARGETS_V1__MITOSTRIP",
         regulon   = "M_b__MITOSTRIP")
EST_PRIMARY  <- "M_a"
EST_LABEL    <- c(clean = "MYC_UP.V1_UP (1.5%)", M_a = "M_a (14.8%)",
                  entangled = "HALLMARK (23.5%)", regulon = "CollecTRI regulon")
EST_EXCLUDED <- c("M_b__FULL", "M_b__PROLIFSTRIP")

# The twelve (4a.1), blocked by function as on Fig. 4g.
PRO    <- c("BAD", "BIK", "BBC3", "BID", "BCL2L11", "BMF", "PMAIP1")
ANTI   <- c("BCL2L1", "MCL1", "BCL2", "BCL2L2", "BCL2A1")
TWELVE <- c(PRO, ANTI)

# The signatures with their signature signs (4a.4, 4a.6), and the arms
# (section 4). S6 is E34's CONFIG_SIGN verbatim.
SIGS <- list(
  S6 = c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1),
  S5 = c(BBC3 = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1),
  S7 = c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, BMF = -1, MCL1 = -1))
ARMS <- list(arm_trigger  = c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1),
             arm_guardian = c(BCL2L1 = 1, MCL1 = -1))
READOUTS <- c(names(SIGS), names(ARMS), TWELVE)          # 17
READ_PRIMARY <- "S6"
READ_LABEL <- c(S6 = "Signature (S6)", S5 = "S5, without BID",
                S7 = "S7, with BMF", arm_trigger = "Trigger arm",
                arm_guardian = "Guardian arm",
                stats::setNames(TWELVE, TWELVE))

# 4a.7: the reconstructed rules, checked against E10's component_cor below.
RULE_THR <- 0.20

# Section 3. The spine of each rung, by cohort. SCAN-B has no purity or
# leukocyte estimate (trap 2). m0 is reported and never read.
RUNG_COV <- list(
  TCGA     = list(m0 = character(0), m1 = c("PAM50", "purity", "leuko"),
                  m3 = c("PAM50", "purity", "leuko", "PROLIF")),
  `SCAN-B` = list(m0 = character(0), m1 = "PAM50", m3 = c("PAM50", "PROLIF")))
RUNGS      <- c("m0", "m1", "m3")
RUNGS_READ <- c("m1", "m3")

# E34's m1 complete-case set, which is the analysis set here (section 3).
SET_COLS <- list(TCGA     = c("quad_gsva", "quad_mitopps", "PAM50", "purity",
                              "leuko"),
                 `SCAN-B` = c("quad_gsva", "quad_mitopps", "PAM50"))
EXPECT_N <- c(TCGA = 938L, `SCAN-B` = 3143L)

BAND_SEP <- 0.60   # E39's bands (G3)
BAND_ID  <- 0.80
G4_TOL   <- 1e-10
BOOT_N   <- 2000L  # section 9
CI_Q     <- 0.975

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.fmt_ci  <- function(e, lo, hi) sprintf("%+.3f [%+.3f, %+.3f]", e, lo, hi)

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

# E34's Blom, verbatim (scripts/E34_quadrant_configuration.R:410-415), so a
# score built here is on exactly E34's scale.
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

for (p in c(PATH_E34, PATH_E10, PATH_NEW, PATH_SCS, PATH_SETS, PATH_TCGA_MITO,
            PATH_SCANB_PHENO, PATH_TCGA_LINEAR, PATH_SCANB_LINEAR)) {
  .stop_if(file.exists(p), "absent: ", p)
}
e34  <- readRDS(PATH_E34)
FR   <- e34$analysis_frames
mito <- readRDS(PATH_TCGA_MITO)
nw   <- readRDS(PATH_NEW)
sc   <- readRDS(PATH_SCS)
sd_  <- readRDS(PATH_SETS)
sph  <- readRDS(PATH_SCANB_PHENO)
e10  <- readRDS(PATH_E10)

ID <- list(TCGA = colnames(mito$gsva_arms), `SCAN-B` = colnames(sc$gsva_arms))
.assert("E34's frames are the full cohorts, in the score objects' order",
        identical(FR$TCGA$sample_id, ID$TCGA) &&
          identical(FR$`SCAN-B`$sample_id, ID$`SCAN-B`),
        n = length(ID$TCGA) + length(ID$`SCAN-B`))

# =============================================================================
# 2. The selection rules (4a.6, 4a.7), checked against their source
# =============================================================================
# The table the six were read from is E10's component_cor: partial Spearman
# on PROLIF_DISJOINT, axes OXPHOS (GSVA) and MYC (M_a). Each declared rule is
# applied to it at FULL precision and must return exactly the declared genes
# with exactly their signature signs.
message("\n2. the selection rules, applied to E10's component_cor")

mc_symbols <- suppressWarnings(readxl::read_excel(
  PATH_MITOCARTA, sheet = "A Human MitoCarta3.0"))$Symbol
CC <- e10$component_cor %>%
  dplyr::filter(axis %in% c("OXPHOS", "MYC"), gene %in% TWELVE) %>%
  dplyr::select(axis, gene, cohort, rho) %>%
  tidyr::pivot_wider(names_from = c(axis, cohort), values_from = rho) %>%
  dplyr::mutate(
    ox_mean  = (OXPHOS_TCGA + `OXPHOS_SCAN-B`) / 2,
    d_tcga   = OXPHOS_TCGA - MYC_TCGA,
    d_scanb  = `OXPHOS_SCAN-B` - `MYC_SCAN-B`,
    mitocarta = gene %in% mc_symbols,
    rule_S6  = mitocarta & abs(ox_mean) >= RULE_THR,
    rule_S7  = abs(ox_mean) >= RULE_THR,
    rule_S5  = abs(d_tcga) >= RULE_THR | abs(d_scanb) >= RULE_THR)
.assert("component_cor holds all twelve genes on both axes in both cohorts",
        nrow(CC) == 12L && !anyNA(CC$ox_mean), n = nrow(CC))
.rule_ok <- function(sig, col) {
  sel <- CC$gene[CC[[col]]]
  setequal(sel, names(SIGS[[sig]])) &&
    all(sign(CC$ox_mean[match(names(SIGS[[sig]]), CC$gene)]) ==
          unname(SIGS[[sig]]))
}
for (s in names(SIGS)) {
  .assert(paste0(s, ": its declared rule selects exactly its genes, with ",
                 "exactly their signature signs"),
          .rule_ok(s, paste0("rule_", s)), n = 12L,
          detail = paste(CC$gene[CC[[paste0("rule_", s)]]], collapse = ", "))
}
.assert("BMF is the one gene of the twelve outside MitoCarta 3.0",
        identical(CC$gene[!CC$mitocarta], "BMF"), n = 12L)
MARGINS <- tibble::tibble(
  what = c("BID in S6 (|mean rho| - 0.20)", "BBC3 in S5 (TCGA |d| - 0.20)"),
  margin = c(abs(CC$ox_mean[CC$gene == "BID"]) - RULE_THR,
             abs(CC$d_tcga[CC$gene == "BBC3"]) - RULE_THR))
.pr(MARGINS %>% dplyr::mutate(margin = round(margin, 5)))

# =============================================================================
# 3. G2, self-overlap: none of the twelve inside any predictor
# =============================================================================
message("\n3. G2: self-overlap")

PRED_SETS <- list(
  MYC_UP.V1_UP__MITOSTRIP            = sd_$myc_sets[["MYC_UP.V1_UP__MITOSTRIP"]],
  FELSHER__MITOSTRIP                 = sd_$myc_sets[["FELSHER__MITOSTRIP"]],
  HALLMARK_MYC_TARGETS_V1__MITOSTRIP = sd_$myc_sets[["HALLMARK_MYC_TARGETS_V1__MITOSTRIP"]],
  M_b__MITOSTRIP                     = sd_$collectri_sets[["M_b__MITOSTRIP"]],
  `OXPHOS subunits`                  = mito$arm_sets[[OX_ARM]],
  PROLIF_DISJOINT                    = mito$covariate_sets[[PROLIF_SET]])
EXCL_SETS <- lapply(stats::setNames(EST_EXCLUDED, EST_EXCLUDED),
                    function(nm) sd_$collectri_sets[[nm]])
.assert("every predictor set is present and non-empty",
        all(lengths(PRED_SETS) > 0L) && all(lengths(EXCL_SETS) > 0L),
        n = length(PRED_SETS) + length(EXCL_SETS))
G2 <- dplyr::bind_rows(
  tibble::tibble(set = names(PRED_SETS), role = "fitted",
                 n = lengths(PRED_SETS)),
  tibble::tibble(set = names(EXCL_SETS), role = "excluded",
                 n = lengths(EXCL_SETS))) %>%
  dplyr::mutate(
    of_twelve = vapply(c(PRED_SETS, EXCL_SETS)[set],
                       function(s) length(intersect(TWELVE, s)), integer(1)),
    which = vapply(c(PRED_SETS, EXCL_SETS)[set],
                   function(s) paste(intersect(TWELVE, s), collapse = ", "),
                   character(1)))
.pr(G2)
.assert("G2: no fitted predictor contains any of the twelve",
        all(G2$of_twelve[G2$role == "fitted"] == 0L), n = sum(G2$role == "fitted"))
.assert("G2: both excluded estimators do contain some of the twelve",
        all(G2$of_twelve[G2$role == "excluded"] > 0L),
        n = sum(G2$role == "excluded"))

# =============================================================================
# 4. G4, the readouts rebuilt from expression, and S5 / S7 built the same way
# =============================================================================
# E34's recipe (scripts/E34_quadrant_configuration.R:248-262): z of each gene
# across the WHOLE cohort, times its signature sign, averaged, then Blom over
# the whole cohort. It must reproduce E34's saved S6, both arms and all twelve
# genes before S5 or S7 is built by it.
message("\n4. G4: the readouts rebuilt from expression")

.log2_twelve <- function(path, ids, map) {
  lin <- readRDS(path)
  .stop_if(identical(lin$scale, "linear_deseq2_normalised"),
           basename(path), " is not on the linear DESeq2-normalised scale")
  G  <- log2(lin$mat[, ids, drop = FALSE] + 1)
  gr <- .gene_rows(TWELVE, G, .symbol_resolver(rownames(G), map))
  .stop_if(length(gr$missing) == 0L, basename(path), ": did not resolve ",
           paste(gr$missing, collapse = ", "))
  gr$mat[TWELVE, , drop = FALSE]
}
GM <- list(TCGA     = .log2_twelve(PATH_TCGA_LINEAR, ID$TCGA, NULL),
           `SCAN-B` = .log2_twelve(PATH_SCANB_LINEAR, ID$`SCAN-B`,
                                   sph$symbol_map))
invisible(gc(verbose = FALSE))

.score <- function(M, signs) {
  X  <- M[names(signs), , drop = FALSE]
  gv <- apply(X, 1L, stats::var)
  .stop_if(all(gv > 0), "zero-variance gene: ",
           paste(names(gv)[!(gv > 0)], collapse = ", "))
  Z <- t(scale(t(X)))
  .blom(as.numeric(colMeans(sweep(Z, 1L, signs, "*"))))
}
BUILT <- lapply(COHORTS, function(coh) {
  M <- GM[[coh]]
  out <- c(lapply(SIGS, function(s) .score(M, s)),
           lapply(ARMS, function(s) .score(M, s)),
           lapply(stats::setNames(TWELVE, TWELVE),
                  function(g) .blom(as.numeric(M[g, ]))))
  out
})
names(BUILT) <- COHORTS

G4 <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  pairs <- c(S6 = "config_comp", arm_trigger = "arm_trigger",
             arm_guardian = "arm_guardian",
             stats::setNames(paste0("g_", TWELVE), TWELVE))
  tibble::tibble(cohort = coh, readout = names(pairs), e34_column = unname(pairs),
                 max_abs_diff = vapply(names(pairs), function(r)
                   max(abs(BUILT[[coh]][[r]] - FR[[coh]][[pairs[[r]]]])),
                   numeric(1)))
}))
.pr(G4 %>% dplyr::group_by(cohort) %>%
      dplyr::summarise(readouts = dplyr::n(), max_abs_diff = max(max_abs_diff),
                       .groups = "drop"))
.assert("G4: S6, both arms and all twelve genes rebuilt equal E34's saved values",
        all(G4$max_abs_diff < G4_TOL), n = nrow(G4),
        max_dev = max(G4$max_abs_diff))
message("   G4 passed, so S5 and S7 are built by the same function.")

# =============================================================================
# 5. Exposures, Blom-scored once per cohort
# =============================================================================
message("\n5. exposures")

.row <- function(M, r, ids) {
  .stop_if(r %in% rownames(M), "row absent: ", r)
  v <- as.numeric(M[r, ids])
  .stop_if(!anyNA(v), "NA in ", r)
  v
}
RAW <- list(
  TCGA = list(
    OX_gsva    = .row(mito$gsva_arms,    OX_ARM, ID$TCGA),
    OX_mitopps = .row(mito$mitopps_arms, OX_ARM, ID$TCGA),
    MYC_clean     = .row(nw$tcga_gsva_new, EST[["clean"]], ID$TCGA),
    MYC_M_a       = .row(nw$tcga_gsva_new, EST[["M_a"]], ID$TCGA),
    MYC_entangled = .row(nw$tcga_gsva_new, EST[["entangled"]], ID$TCGA),
    MYC_regulon   = .row(nw$tcga_M_b_variants, EST[["regulon"]], ID$TCGA),
    PROLIF     = .row(mito$gsva_cov, PROLIF_SET, ID$TCGA)),
  `SCAN-B` = list(
    OX_gsva    = .row(sc$gsva_arms,    OX_ARM, ID$`SCAN-B`),
    OX_mitopps = .row(sc$mitopps_arms, OX_ARM, ID$`SCAN-B`),
    MYC_clean     = .row(sc$gsva_new, EST[["clean"]], ID$`SCAN-B`),
    MYC_M_a       = .row(sc$gsva_new, EST[["M_a"]], ID$`SCAN-B`),
    MYC_entangled = .row(sc$gsva_new, EST[["entangled"]], ID$`SCAN-B`),
    MYC_regulon   = .row(sc$M_b_variants, EST[["regulon"]], ID$`SCAN-B`),
    PROLIF     = .row(sc$gsva_cov, PROLIF_SET, ID$`SCAN-B`)))
for (coh in COHORTS) {
  .assert(paste0(coh, ": PROLIF Blom-scored here equals E34's, exactly"),
          max(abs(.blom(RAW[[coh]]$PROLIF) - FR[[coh]]$PROLIF)) < G4_TOL,
          n = length(ID[[coh]]),
          max_dev = max(abs(.blom(RAW[[coh]]$PROLIF) - FR[[coh]]$PROLIF)))
}

# =============================================================================
# 6. The analysis frame (R rule 4), and G1
# =============================================================================
message("\n6. the analysis frame, and G1")

.frame <- function(coh) {
  f <- FR[[coh]]
  r <- RAW[[coh]]
  d <- tibble::tibble(
    sample_id = f$sample_id, cohort = coh,
    in_set = stats::complete.cases(f[, SET_COLS[[coh]]]),
    quad_gsva = f$quad_gsva, quad_mitopps = f$quad_mitopps,
    PAM50 = f$PAM50, purity = f$purity, leuko = f$leuko, PROLIF = f$PROLIF)
  for (nm in setdiff(names(r), "PROLIF")) {
    d[[paste0(nm, "_raw")]] <- r[[nm]]
    d[[nm]] <- .blom(r[[nm]])
  }
  d$S6 <- f$config_comp
  d$S5 <- BUILT[[coh]]$S5
  d$S7 <- BUILT[[coh]]$S7
  d$arm_trigger  <- f$arm_trigger
  d$arm_guardian <- f$arm_guardian
  for (g in TWELVE) d[[g]] <- f[[paste0("g_", g)]]
  d
}
DF <- lapply(stats::setNames(COHORTS, COHORTS), .frame)
for (coh in COHORTS) {
  .assert(paste0(coh, ": the analysis set is E34's m1 set, n ", EXPECT_N[[coh]]),
          sum(DF[[coh]]$in_set) == EXPECT_N[[coh]], n = nrow(DF[[coh]]),
          detail = paste("in_set", sum(DF[[coh]]$in_set)))
}

# G1: E34's quadrant labels must be a threshold split of the continuous scores
# read here, on every cohort and instrument.
G1 <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  dplyr::bind_rows(lapply(names(INSTR), function(inst) {
    d <- DF[[coh]]
    q <- as.character(d[[paste0("quad_", inst)]])
    ok <- !is.na(q)
    myc_hi <- grepl("MYC-high", q[ok])
    ox_hi  <- grepl("OXPHOS-high", q[ok])
    m <- d$MYC_M_a_raw[ok]
    o <- d[[paste0("OX_", inst, "_raw")]][ok]
    tibble::tibble(cohort = coh, instrument = inst, n_labelled = sum(ok),
                   myc_split_monotone = max(m[!myc_hi]) < min(m[myc_hi]),
                   ox_split_monotone  = max(o[!ox_hi]) < min(o[ox_hi]))
  }))
}))
.pr(G1)
.assert("G1: E34's quadrants are threshold splits of these continuous scores",
        all(G1$myc_split_monotone & G1$ox_split_monotone), n = nrow(G1))

# =============================================================================
# 7. G3, identifiability, with E39's bands
# =============================================================================
message("\n7. G3: identifiability")

.band <- function(r) ifelse(abs(r) >= BAND_ID, "NOT IDENTIFIABLE",
                     ifelse(abs(r) >= BAND_SEP, "not cleanly separable", "fine"))
G3 <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  d <- DF[[coh]][DF[[coh]]$in_set, ]
  dplyr::bind_rows(lapply(names(INSTR), function(inst) {
    dplyr::bind_rows(lapply(names(EST), function(e) {
      o <- d[[paste0("OX_", inst)]]
      m <- d[[paste0("MYC_", e)]]
      tibble::tibble(cohort = coh, instrument = inst, estimator = e,
                     rho_OX_MYC = stats::cor(o, m, method = "spearman"),
                     rho_MYC_PROLIF = stats::cor(m, d$PROLIF, method = "spearman"),
                     rho_OX_PROLIF = stats::cor(o, d$PROLIF, method = "spearman"))
    }))
  }))
})) %>%
  dplyr::mutate(band_OX_MYC = .band(rho_OX_MYC),
                band_MYC_PROLIF = .band(rho_MYC_PROLIF),
                band_OX_PROLIF = .band(rho_OX_PROLIF),
                void_m3 = abs(rho_MYC_PROLIF) >= BAND_ID)
.pr(G3 %>% dplyr::transmute(cohort, instrument, estimator,
                            OX_MYC = round(rho_OX_MYC, 3),
                            MYC_PROLIF = round(rho_MYC_PROLIF, 3),
                            OX_PROLIF = round(rho_OX_PROLIF, 3), void_m3))
VOID <- G3 %>% dplyr::distinct(cohort, estimator, void_m3)
.assert("G3: the m3 void flag depends on cohort and estimator only",
        nrow(VOID) == length(COHORTS) * length(EST), n = nrow(VOID))
message("   ", .sub(paste0(
  "m3 cells VOID by G3 (rho(MYC, PROLIF) >= 0.80; reported, never read): ",
  if (any(VOID$void_m3)) paste(VOID$cohort[VOID$void_m3],
                               VOID$estimator[VOID$void_m3], collapse = "; ")
  else "none")))

# =============================================================================
# 8. The fits: 2 cohorts x 2 instruments x 4 estimators x 3 rungs x 17 readouts
# =============================================================================
message("\n8. fits")

GRID <- tidyr::expand_grid(cohort = COHORTS, instrument = names(INSTR),
                           estimator = names(EST), rung = RUNGS,
                           readout = READOUTS)
GRID$fit_id <- paste(GRID$cohort, GRID$instrument, GRID$estimator, GRID$rung,
                     GRID$readout, sep = " | ")

.fit_one <- function(coh, inst, est, rung, readout) {
  d   <- DF[[coh]][DF[[coh]]$in_set, ]
  cov <- RUNG_COV[[coh]][[rung]]
  dd  <- data.frame(y = d[[readout]], OX = d[[paste0("OX_", inst)]],
                    MYC = d[[paste0("MYC_", est)]])
  for (cv in cov) dd[[cv]] <- d[[cv]]
  if ("PAM50" %in% cov) dd$PAM50 <- droplevels(dd$PAM50)
  .stop_if(!anyNA(dd), "NA inside the analysis set: ", coh, " ", readout)
  stats::lm(stats::reformulate(c("OX", "MYC", cov), response = "y"), data = dd)
}
.vif_tidy <- function(f) {
  v <- car::vif(f)
  if (is.matrix(v)) {
    tibble::tibble(term = rownames(v), gvif = v[, "GVIF"], df = v[, "Df"],
                   gvif_adj = v[, "GVIF^(1/(2*Df))"])
  } else {
    tibble::tibble(term = names(v), gvif = unname(v), df = 1,
                   gvif_adj = sqrt(unname(v)))
  }
}

COEF <- vector("list", nrow(GRID))
VCOV <- vector("list", nrow(GRID))
VIF  <- vector("list", nrow(GRID))
DIG  <- vector("list", nrow(GRID))
FITS_S6 <- list()
names(VCOV) <- GRID$fit_id
for (i in seq_len(nrow(GRID))) {
  g <- GRID[i, ]
  f <- .fit_one(g$cohort, g$instrument, g$estimator, g$rung, g$readout)
  s <- summary(f)$coefficients
  tq <- stats::qt(CI_Q, f$df.residual)
  COEF[[i]] <- tibble::tibble(
    fit_id = g$fit_id, cohort = g$cohort, instrument = g$instrument,
    estimator = g$estimator, rung = g$rung, readout = g$readout,
    term = rownames(s), estimate = s[, 1], se = s[, 2], t = s[, 3],
    p = s[, 4], ci_lo = s[, 1] - tq * s[, 2], ci_hi = s[, 1] + tq * s[, 2],
    n = stats::nobs(f), df = f$df.residual)
  VCOV[[i]] <- stats::vcov(f)
  VIF[[i]]  <- .vif_tidy(f) %>% dplyr::mutate(fit_id = g$fit_id, .before = 1)
  DIG[[i]]  <- tibble::tibble(fit_id = g$fit_id,
                              digest_coef = digest::digest(stats::coef(f)),
                              digest_vcov = digest::digest(stats::vcov(f)))
  if (g$readout == READ_PRIMARY) FITS_S6[[g$fit_id]] <- f
  if (i %% 204L == 0L) message("   ", i, " of ", nrow(GRID), " fits")
}
COEF <- dplyr::bind_rows(COEF)
VIF  <- dplyr::bind_rows(VIF)
FIT_DIGESTS <- dplyr::bind_rows(DIG)
.assert("816 fits, each with an OX and a MYC coefficient",
        nrow(GRID) == 816L &&
          all(table(COEF$fit_id[COEF$term %in% c("OX", "MYC")]) == 2L) &&
          !anyNA(COEF$estimate[COEF$term %in% c("OX", "MYC")]),
        n = nrow(GRID))
.assert("every fit uses its cohort's whole analysis set",
        all(COEF$n == EXPECT_N[COEF$cohort]), n = nrow(COEF))
.assert("48 S6 fits are kept as objects (section 13)",
        length(FITS_S6) == 48L, n = length(FITS_S6))

# The VIF rule (R rule 5): GVIF^(1/(2*Df)) against sqrt(5), every term.
VIF <- VIF %>% dplyr::mutate(above_caveat = gvif_adj > sqrt(5))
message("   largest GVIF^(1/(2*Df)) anywhere: ", round(max(VIF$gvif_adj), 3),
        " (caveat above ", round(sqrt(5), 3), "); terms above: ",
        sum(VIF$above_caveat))

# =============================================================================
# 9. The two coefficients of every fit, and the reading rule (section 7)
# =============================================================================
message("\n9. the reading rule")

RES <- COEF %>%
  dplyr::filter(term %in% c("OX", "MYC")) %>%
  dplyr::select(fit_id, cohort, instrument, estimator, rung, readout, term,
                estimate, ci_lo, ci_hi, p, n) %>%
  tidyr::pivot_wider(names_from = term,
                     values_from = c(estimate, ci_lo, ci_hi, p)) %>%
  dplyr::rename(b_OX = estimate_OX, ox_lo = ci_lo_OX, ox_hi = ci_hi_OX,
                b_MYC = estimate_MYC, myc_lo = ci_lo_MYC,
                myc_hi = ci_hi_MYC) %>%
  dplyr::left_join(VOID, by = c("cohort", "estimator")) %>%
  dplyr::mutate(void = rung == "m3" & void_m3,
                ratio_MYC_to_OX = b_MYC / b_OX) %>%
  dplyr::select(-void_m3)

# Section 7, transcribed. Applied in this order; the first match is the label.
# Counts are over NON-VOID cells.
.label <- function(rd) {
  r <- RES %>% dplyr::filter(readout == rd, rung %in% RUNGS_READ, !void)
  P <- r %>% dplyr::filter(estimator == EST_PRIMARY)
  if (all(P$myc_lo > 0)) {
    lad <- r %>% dplyr::filter(estimator != EST_PRIMARY) %>%
      dplyr::group_by(estimator) %>%
      dplyr::summarise(frac_pos = mean(b_MYC > 0), .groups = "drop")
    return(if (all(lad$frac_pos >= 0.75)) "MYC ADDS" else "ESTIMATOR-DEPENDENT")
  }
  by_i <- P %>% dplyr::group_by(instrument) %>%
    dplyr::summarise(all_above = all(myc_lo > 0), none_above = all(myc_lo <= 0),
                     .groups = "drop")
  split <- nrow(by_i) == 2L &&
    ((by_i$all_above[1] && by_i$none_above[2]) ||
       (by_i$all_above[2] && by_i$none_above[1]))
  if (split) return("INSTRUMENT-SPLIT")
  if (mean(P$myc_lo <= 0) >= 0.75) return("OXPHOS ONLY")
  "MIXED"
}
LABELS <- tibble::tibble(readout = names(SIGS),
                         label = vapply(names(SIGS), .label, character(1)),
                         role = c(S6 = "THE VERDICT", S5 = "descriptive",
                                  S7 = "descriptive")[names(SIGS)])
VERDICT <- LABELS$label[LABELS$readout == READ_PRIMARY]

PRIMARY_CELLS <- RES %>%
  dplyr::filter(readout == READ_PRIMARY, estimator == EST_PRIMARY,
                rung %in% RUNGS_READ) %>%
  dplyr::arrange(match(cohort, COHORTS), instrument, rung)
message("\n   S6, M_a: beta_MYC|OX (READ) and beta_OX|MYC (NOT read)")
.pr(PRIMARY_CELLS %>% dplyr::transmute(
  cohort, instrument, rung, n,
  beta_MYC = .fmt_ci(b_MYC, myc_lo, myc_hi), CI_above_0 = myc_lo > 0,
  beta_OX = .fmt_ci(b_OX, ox_lo, ox_hi), void))

LADDER <- RES %>%
  dplyr::filter(readout %in% names(SIGS), rung %in% RUNGS_READ, !void) %>%
  dplyr::group_by(readout, estimator) %>%
  dplyr::summarise(cells = dplyr::n(), positive = sum(b_MYC > 0),
                   ci_above_0 = sum(myc_lo > 0), ci_below_0 = sum(myc_hi < 0),
                   .groups = "drop") %>%
  dplyr::arrange(match(readout, names(SIGS)), match(estimator, names(EST)))
message("\n   the estimator ladder for beta_MYC|OX (N9):")
.pr(LADDER)

message("\n   ", strrep("-", 74))
message("   VERDICT, S6 (section 7): ", VERDICT)
message("   S5 (descriptive; BIASED AGAINST MYC BY CONSTRUCTION, 4a.7): ",
        LABELS$label[LABELS$readout == "S5"])
message("   S7 (descriptive): ", LABELS$label[LABELS$readout == "S7"])
if (any(LABELS$label[LABELS$readout != READ_PRIMARY] != VERDICT)) {
  message("   ", .sub(paste0(
    "A VARIANT'S LABEL DIFFERS FROM S6's. Section 7: the result note says so ",
    "in its first section, and no manuscript sentence states S6's verdict ",
    "without naming that dependence.")))
}
message("   ", strrep("-", 74))

RATIO <- RES %>%
  dplyr::filter(readout %in% names(SIGS), estimator == EST_PRIMARY,
                rung %in% RUNGS_READ) %>%
  dplyr::select(readout, cohort, instrument, rung, b_MYC, b_OX,
                ratio_MYC_to_OX)
message("   ", .sub(paste0(
  "beta_MYC / beta_OX is saved in $ratio as a DESCRIPTION ONLY. beta_OX is ",
  "inflated by construction, so the ratio UNDERSTATES MYC's share, and it is ",
  "never quoted without this sentence.")))

# =============================================================================
# 10. Diagnostics D1 (leakage) and D2 (rung movement), section 8
# =============================================================================
message("\n10. diagnostics")

QL <- c(Q1 = "Q1 MYC-low OXPHOS-low",  Q2 = "Q2 MYC-low OXPHOS-high",
        Q3 = "Q3 MYC-high OXPHOS-low", Q4 = "Q4 MYC-high OXPHOS-high")
D1 <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  d <- DF[[coh]][DF[[coh]]$in_set, ]
  dplyr::bind_rows(lapply(names(INSTR), function(inst) {
    q  <- as.character(d[[paste0("quad_", inst)]])
    ox <- d[[paste0("OX_", inst)]]
    mq <- vapply(QL, function(l) mean(ox[q == l]), numeric(1))
    b  <- RES$b_OX[RES$readout == READ_PRIMARY & RES$estimator == EST_PRIMARY &
                     RES$rung == "m1" & RES$cohort == coh &
                     RES$instrument == inst]
    obs <- e34$contrasts %>%
      dplyr::filter(subtype == "all", readout_col == "config_comp",
                    model == "m1", cohort == coh, instrument == inst)
    tibble::tibble(
      cohort = coh, instrument = inst,
      contrast = c("Q4 - Q2", "Q3 - Q1"),
      ox_gap = c(mq[["Q4"]] - mq[["Q2"]], mq[["Q3"]] - mq[["Q1"]]),
      beta_OX_m1 = b,
      leakage = c(mq[["Q4"]] - mq[["Q2"]], mq[["Q3"]] - mq[["Q1"]]) * b,
      e34_observed = obs$est[match(c("Q4 - Q2", "Q3 - Q1"), obs$contrast)])
  }))
}))
.pr(D1 %>% dplyr::mutate(dplyr::across(c(ox_gap, beta_OX_m1, leakage,
                                         e34_observed), ~ round(.x, 3))))
message("   ", .sub(paste0(
  "D1: 'leakage' is the MYC contrast the leftover OXPHOS inside each half ",
  "would produce on its own, beside E34's observed MYC contrast. It replaces ",
  "declaration 2(d)'s approximation with a measurement.")))

D2 <- RES %>%
  dplyr::filter(readout %in% c(names(SIGS), names(ARMS)),
                estimator == EST_PRIMARY) %>%
  dplyr::select(readout, cohort, instrument, rung, b_MYC, b_OX) %>%
  tidyr::pivot_wider(names_from = rung, values_from = c(b_MYC, b_OX))

# =============================================================================
# 11. Section 9: BCL2L1 against BBC3, paired tumour-level bootstrap
# =============================================================================
message("\n11. section 9: BCL2L1 against BBC3")

SEC9 <- list()
BOOT <- list()
k <- 0L
for (coh in COHORTS) for (inst in names(INSTR)) {
  k <- k + 1L
  d   <- DF[[coh]][DF[[coh]]$in_set, ]
  cov <- RUNG_COV[[coh]][["m1"]]
  dd  <- data.frame(OX = d[[paste0("OX_", inst)]], MYC = d$MYC_M_a)
  for (cv in cov) dd[[cv]] <- d[[cv]]
  if ("PAM50" %in% cov) dd$PAM50 <- droplevels(dd$PAM50)
  X  <- stats::model.matrix(stats::reformulate(c("OX", "MYC", cov)), data = dd)
  y1 <- d$BCL2L1
  y2 <- d$BBC3
  b1 <- stats::lm.fit(X, y1)$coefficients[["OX"]]
  b2 <- stats::lm.fit(X, y2)$coefficients[["OX"]]
  r1 <- RES %>% dplyr::filter(readout == "BCL2L1", estimator == EST_PRIMARY,
                              rung == "m1", cohort == coh, instrument == inst)
  r2 <- RES %>% dplyr::filter(readout == "BBC3", estimator == EST_PRIMARY,
                              rung == "m1", cohort == coh, instrument == inst)
  .assert(paste0("section 9, ", coh, " ", inst, ": the bootstrap's design ",
                 "matrix reproduces the reported OX coefficients"),
          abs(b1 - r1$b_OX) < G4_TOL && abs(b2 - r2$b_OX) < G4_TOL, n = 2L,
          max_dev = max(abs(b1 - r1$b_OX), abs(b2 - r2$b_OX)))
  set.seed(PROJECT_SEED + k)
  n <- nrow(X)
  reps <- numeric(BOOT_N)
  for (b in seq_len(BOOT_N)) {
    idx <- sample.int(n, n, replace = TRUE)
    Xb  <- X[idx, , drop = FALSE]
    reps[b] <- stats::lm.fit(Xb, y1[idx])$coefficients[["OX"]] -
               stats::lm.fit(Xb, y2[idx])$coefficients[["OX"]]
  }
  .stop_if(!anyNA(reps), "a bootstrap replicate is NA: ", coh, " ", inst)
  ci <- stats::quantile(reps, c(0.025, 0.975), names = FALSE, type = 7)
  SEC9[[k]] <- tibble::tibble(
    cohort = coh, instrument = inst, n = n, seed = PROJECT_SEED + k,
    b_OX_BCL2L1 = r1$b_OX, lo_BCL2L1 = r1$ox_lo, hi_BCL2L1 = r1$ox_hi,
    b_OX_BBC3 = r2$b_OX, lo_BBC3 = r2$ox_lo, hi_BBC3 = r2$ox_hi,
    delta = b1 - b2, boot_se = stats::sd(reps), ci_lo = ci[1], ci_hi = ci[2],
    steeper = (b1 - b2) > 0 && ci[1] > 0)
  BOOT[[paste(coh, inst)]] <- reps
}
SEC9 <- dplyr::bind_rows(SEC9)
STEEPER_LICENSED <- all(SEC9$steeper)
.pr(SEC9 %>% dplyr::transmute(
  cohort, instrument,
  BCL2L1 = .fmt_ci(b_OX_BCL2L1, lo_BCL2L1, hi_BCL2L1),
  BBC3 = .fmt_ci(b_OX_BBC3, lo_BBC3, hi_BBC3),
  delta = .fmt_ci(delta, ci_lo, ci_hi), steeper))
message("   ", .sub(paste0(
  "Section 9: 'more steeply' is licensed only if delta > 0 with its bootstrap ",
  "CI above zero in ALL FOUR cells. Licensed: ", STEEPER_LICENSED, ".")))

# =============================================================================
# 12. Displays A, B and C (section 8a), fixed before any fit
# =============================================================================
message("\n12. displays")

COL_COHORT <- c(TCGA = "#0072B2", `SCAN-B` = "#CC79A7")
SHAPE_INSTR <- c(gsva = 21, mitopps = 22)
COL_EST <- c(clean = "#009E73", M_a = "#0072B2", entangled = "#D55E00",
             regulon = "#CC79A7")
BASE_PT <- 6.5
TXT_MM  <- 6 / ggplot2::.pt
SUB_MM  <- 5.5 / ggplot2::.pt

# Row geometry shared by A and B: the scores, then the twelve by function.
ROW_Y <- c(S6 = 18.2, S5 = 17.2, S7 = 16.2, arm_trigger = 15.2,
           arm_guardian = 14.2,
           stats::setNames(seq(12.6, 6.6, by = -1), PRO),
           stats::setNames(5:1, ANTI))
Y_LIM <- c(0.4, 19.4)
SEP_Y <- c(13.4, 5.8)
ROWS <- tibble::tibble(readout = names(ROW_Y), y = unname(ROW_Y)) %>%
  dplyr::mutate(
    label = unname(READ_LABEL[readout]),
    face = ifelse(readout %in% TWELVE, "italic",
                  ifelse(readout == READ_PRIMARY, "bold", "plain")),
    sign = dplyr::case_when(
      readout %in% names(SIGS$S6) & SIGS$S6[readout] > 0 ~ "+",
      readout %in% names(SIGS$S6) & SIGS$S6[readout] < 0 ~ "-",
      readout == "BMF" ~ "-*",
      TRUE ~ ""))

.label_panel <- function() {
  blocks <- tibble::tibble(
    label = c("Scores", "Pro-apoptotic", "Anti-apoptotic"),
    y0 = c(min(ROW_Y[names(c(SIGS, ARMS))]), min(ROW_Y[PRO]),
           min(ROW_Y[ANTI])) - 0.4,
    y1 = c(max(ROW_Y[names(c(SIGS, ARMS))]), max(ROW_Y[PRO]),
           max(ROW_Y[ANTI])) + 0.4)
  ggplot2::ggplot(ROWS) +
    ggplot2::geom_segment(data = blocks, ggplot2::aes(x = 0.12, xend = 0.12,
                                                      y = y0, yend = y1),
                          linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_text(data = blocks, ggplot2::aes(x = 0.05, y = (y0 + y1) / 2,
                                                   label = label),
                       angle = 90, size = SUB_MM, colour = "grey25") +
    ggplot2::geom_text(ggplot2::aes(x = 0.18, y = y, label = label,
                                    fontface = face),
                       hjust = 0, size = TXT_MM) +
    ggplot2::geom_text(ggplot2::aes(x = 0.93, y = y, label = sign),
                       size = TXT_MM, fontface = "bold") +
    ggplot2::annotate("text", x = 0.995, y = max(ROW_Y) + 0.75, hjust = 1,
                      label = "signature\nsign", size = SUB_MM,
                      colour = "grey25", lineheight = 0.9) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::theme_void(base_size = BASE_PT)
}

.theme_forest <- function() {
  ggplot2::theme_bw(base_size = BASE_PT) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95",
                                               colour = "grey80"),
      strip.text = ggplot2::element_text(size = BASE_PT),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.title = ggplot2::element_text(size = BASE_PT - 0.5))
}

# Display A: both coefficients, M_a, m1 filled and m3 open.
OFF_A <- c("TCGA gsva m1" = 0.35, "TCGA gsva m3" = 0.25,
           "TCGA mitopps m1" = 0.15, "TCGA mitopps m3" = 0.05,
           "SCAN-B gsva m1" = -0.05, "SCAN-B gsva m3" = -0.15,
           "SCAN-B mitopps m1" = -0.25, "SCAN-B mitopps m3" = -0.35)
COEF_LEVELS <- c(OX = "OXPHOS, MYC held fixed", MYC = "MYC, OXPHOS held fixed")
PA <- dplyr::bind_rows(
  RES %>% dplyr::transmute(fit_id, cohort, instrument, estimator, rung, readout,
                           coef = "OX", est = b_OX, lo = ox_lo, hi = ox_hi),
  RES %>% dplyr::transmute(fit_id, cohort, instrument, estimator, rung, readout,
                           coef = "MYC", est = b_MYC, lo = myc_lo, hi = myc_hi)) %>%
  dplyr::filter(estimator == EST_PRIMARY, rung %in% RUNGS_READ) %>%
  dplyr::mutate(
    y_mark = unname(ROW_Y[readout]) +
      unname(OFF_A[paste(cohort, instrument, rung)]),
    panel = factor(unname(COEF_LEVELS[coef]), levels = unname(COEF_LEVELS)),
    cohort = factor(cohort, levels = COHORTS))
XL_A <- ceiling(max(abs(c(PA$lo, PA$hi))) * 10) / 10
.assert("Display A: no interval clipped by the shared axis",
        all(PA$lo >= -XL_A & PA$hi <= XL_A), n = nrow(PA),
        detail = sprintf("axis +/-%.1f", XL_A))

.forest_A <- function() {
  ggplot2::ggplot(PA) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = SEP_Y, colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark,
                                       yend = y_mark, colour = cohort),
                          linewidth = 0.3, show.legend = FALSE) +
    ggplot2::geom_point(data = PA %>% dplyr::filter(rung == "m1"),
                        ggplot2::aes(x = est, y = y_mark, colour = cohort,
                                     fill = cohort, shape = instrument),
                        size = 1.1, stroke = 0.35) +
    ggplot2::geom_point(data = PA %>% dplyr::filter(rung == "m3"),
                        ggplot2::aes(x = est, y = y_mark, colour = cohort,
                                     shape = instrument),
                        fill = "white", size = 1.1, stroke = 0.35,
                        show.legend = FALSE) +
    ggplot2::facet_grid(. ~ panel) +
    ggplot2::scale_colour_manual(values = COL_COHORT, name = "Cohort") +
    ggplot2::scale_fill_manual(values = COL_COHORT, name = "Cohort") +
    ggplot2::scale_shape_manual(values = SHAPE_INSTR, name = "OXPHOS by",
                                labels = INSTR) +
    ggplot2::scale_x_continuous(limits = c(-XL_A, XL_A)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::labs(x = "Coefficient, SD per SD (95% CI)") +
    .theme_forest()
}
FIG_A_PLOT <- (.label_panel() | .forest_A()) +
  patchwork::plot_layout(widths = c(1.25, 5.6)) +
  patchwork::plot_annotation(
    subtitle = paste0("MYC activity: M_a. Filled: m1 (spine). Open: m3 ",
                      "(spine + proliferation)."),
    theme = ggplot2::theme(plot.subtitle = ggplot2::element_text(
      size = BASE_PT)))

# Display B: beta_MYC|OX on all four estimators, by cohort and instrument.
OFF_B <- stats::setNames(seq(0.35, -0.35, length.out = 8),
                         paste(rep(names(EST), each = 2), RUNGS_READ))
BLOCKS <- c("TCGA gsva" = "TCGA, GSVA", "TCGA mitopps" = "TCGA, mitoPPS",
            "SCAN-B gsva" = "SCAN-B, GSVA", "SCAN-B mitopps" = "SCAN-B, mitoPPS")
PB <- RES %>%
  dplyr::filter(rung %in% RUNGS_READ) %>%
  dplyr::transmute(fit_id, cohort, instrument, estimator, rung, readout, void,
                   est = b_MYC, lo = myc_lo, hi = myc_hi) %>%
  dplyr::mutate(
    y_mark = unname(ROW_Y[readout]) + unname(OFF_B[paste(estimator, rung)]),
    panel = factor(unname(BLOCKS[paste(cohort, instrument)]),
                   levels = unname(BLOCKS)),
    estimator = factor(estimator, levels = names(EST)))
XL_B <- ceiling(max(abs(c(PB$lo, PB$hi))) * 10) / 10
.assert("Display B: no interval clipped by the shared axis",
        all(PB$lo >= -XL_B & PB$hi <= XL_B), n = nrow(PB),
        detail = sprintf("axis +/-%.1f", XL_B))
.forest_B <- function() {
  ggplot2::ggplot(PB) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = SEP_Y, colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark,
                                       yend = y_mark, colour = estimator,
                                       alpha = void),
                          linewidth = 0.3, show.legend = FALSE) +
    ggplot2::geom_point(data = PB %>% dplyr::filter(rung == "m1"),
                        ggplot2::aes(x = est, y = y_mark, colour = estimator,
                                     fill = estimator, alpha = void),
                        shape = 21, size = 1.0, stroke = 0.35) +
    ggplot2::geom_point(data = PB %>% dplyr::filter(rung == "m3"),
                        ggplot2::aes(x = est, y = y_mark, colour = estimator,
                                     alpha = void),
                        shape = 21, fill = "white", size = 1.0, stroke = 0.35,
                        show.legend = FALSE) +
    ggplot2::facet_grid(. ~ panel) +
    ggplot2::scale_colour_manual(values = COL_EST, labels = EST_LABEL,
                                 name = "MYC estimator") +
    ggplot2::scale_fill_manual(values = COL_EST, labels = EST_LABEL,
                               name = "MYC estimator") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.3),
                                guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-XL_B, XL_B)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::labs(x = "MYC with OXPHOS held fixed, SD per SD (95% CI)") +
    .theme_forest()
}
FIG_B_PLOT <- (.label_panel() | .forest_B()) +
  patchwork::plot_layout(widths = c(1.25, 5.6)) +
  patchwork::plot_annotation(
    subtitle = paste0("Filled: m1. Open: m3 (+ proliferation). Faded: m3 ",
                      "cells void under gate G3, reported and never read."),
    theme = ggplot2::theme(plot.subtitle = ggplot2::element_text(
      size = BASE_PT)))

# Display C: section 9.
C_Y <- c("TCGA gsva" = 4, "TCGA mitopps" = 3, "SCAN-B gsva" = 2,
         "SCAN-B mitopps" = 1)
PC1 <- dplyr::bind_rows(
  SEC9 %>% dplyr::transmute(cohort, instrument, gene = "BCL2L1",
                            est = b_OX_BCL2L1, lo = lo_BCL2L1, hi = hi_BCL2L1),
  SEC9 %>% dplyr::transmute(cohort, instrument, gene = "BBC3",
                            est = b_OX_BBC3, lo = lo_BBC3, hi = hi_BBC3)) %>%
  dplyr::mutate(y_mark = unname(C_Y[paste(cohort, instrument)]) +
                  ifelse(gene == "BCL2L1", 0.12, -0.12),
                gene = factor(gene, levels = c("BCL2L1", "BBC3")))
PC2 <- SEC9 %>% dplyr::transmute(cohort, instrument, est = delta, lo = ci_lo,
                                 hi = ci_hi,
                                 y_mark = unname(C_Y[paste(cohort, instrument)]))
.forest_C1 <- function() {
  ggplot2::ggplot(PC1) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark,
                                       yend = y_mark, colour = gene),
                          linewidth = 0.35, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = est, y = y_mark, colour = gene,
                                     shape = gene), size = 1.3, fill = "white",
                        stroke = 0.4) +
    ggplot2::scale_colour_manual(values = c(BCL2L1 = "grey10", BBC3 = "grey50"),
                                 name = NULL) +
    ggplot2::scale_shape_manual(values = c(BCL2L1 = 16, BBC3 = 24), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
    ggplot2::labs(x = "OXPHOS, MYC held fixed (m1, M_a)",
                  title = "beta_OX|MYC") +
    .theme_forest() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = BASE_PT,
                                                      face = "bold"))
}
.forest_C2 <- function() {
  ggplot2::ggplot(PC2) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark,
                                       yend = y_mark), linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(x = est, y = y_mark), size = 1.3) +
    ggplot2::scale_y_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
    ggplot2::labs(x = "BCL2L1 minus BBC3 (bootstrap 95% CI)",
                  title = "Difference") +
    .theme_forest() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = BASE_PT,
                                                      face = "bold"))
}
.label_C <- function() {
  ggplot2::ggplot(tibble::tibble(y = unname(C_Y), label = unname(BLOCKS[names(C_Y)]))) +
    ggplot2::geom_text(ggplot2::aes(x = 0.98, y = y, label = label), hjust = 1,
                       size = TXT_MM) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
    ggplot2::theme_void(base_size = BASE_PT)
}
FIG_C_PLOT <- (.label_C() | .forest_C1() | .forest_C2()) +
  patchwork::plot_layout(widths = c(1.0, 2.0, 1.6), guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

# --- every drawn point is its stored value, read back from the build -------
.check_points <- function(p, frame, label) {
  b   <- ggplot2::ggplot_build(p)
  geo <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  pts <- dplyr::bind_rows(lapply(which(geo == "GeomPoint"), function(k)
    tibble::as_tibble(b$data[[k]][, c("x", "y", "PANEL")])))
  panel_f <- if ("panel" %in% names(frame)) as.integer(frame$panel) else 1L
  key_f <- paste(panel_f, round(frame$y_mark, 6))
  key_p <- paste(as.integer(pts$PANEL), round(pts$y, 6))
  m <- match(key_p, key_f)
  .assert(paste0(label, ": every point drawn is one stored row, none dropped"),
          nrow(pts) == nrow(frame) && !anyNA(m) && !any(duplicated(m)),
          n = nrow(pts), detail = paste(nrow(pts), "drawn,", nrow(frame), "stored"))
  dev <- max(abs(pts$x - frame$est[m]))
  .assert(paste0(label, ": drawn positions equal the stored estimates"),
          dev < 1e-12, n = nrow(pts), max_dev = dev)
}
.check_points(.forest_A(), PA, "Display A")
.check_points(.forest_B(), PB, "Display B")
.check_points(.forest_C2(), PC2, "Display C, difference")
.check_points(.forest_C1(), PC1, "Display C, coefficients")

.ensure_dir(DIR_FIGURES)
.ensure_dir(DIR_DOC_FIGURES)
.save_fig <- function(fig, stem, w, h) {
  for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
    base <- file.path(dir, if (identical(dir, DIR_DOC_FIGURES))
      paste0(FIG_DATE, "_", stem) else stem)
    ggplot2::ggsave(paste0(base, ".png"), fig, width = w, height = h,
                    dpi = 300, bg = "white")
    ggplot2::ggsave(paste0(base, ".pdf"), fig, width = w, height = h,
                    device = grDevices::cairo_pdf, bg = "white")
  }
}
.save_fig(FIG_A_PLOT, FIG_A, 7.1, 5.4)
.save_fig(FIG_B_PLOT, FIG_B, 7.1, 5.6)
.save_fig(FIG_C_PLOT, FIG_C, 5.2, 2.1)
message("   ", FIG_A, ", ", FIG_B, ", ", FIG_C,
        " to outputs/figures and docs/figures (tracked)")

# =============================================================================
# 13. Legend drafts - DESCRIPTIVE ONLY
# =============================================================================
LEGEND <- list(
  A = paste0(
    "Coefficients from linear models of each readout on the OXPHOS subunit ",
    "score and MYC activity (M_a) together, in TCGA-BRCA (n = ",
    EXPECT_N[["TCGA"]], ") and SCAN-B (n = ", format(EXPECT_N[["SCAN-B"]],
                                                      big.mark = ","),
    "): left, OXPHOS with MYC held fixed; right, MYC with OXPHOS held fixed. ",
    "Rows: the six-gene BCL-2-family signature (S6) and two rule-based ",
    "variants (S5, without BID; S7, with BMF), the signature's trigger and ",
    "guardian arms, and the twelve genes, unsigned; signature signs at left ",
    "(* BMF enters S7 only). Colour, cohort; shape, OXPHOS scored by GSVA or ",
    "mitoPPS; filled, m1 (PAM50 subtype, tumour purity and leukocyte fraction ",
    "in TCGA; PAM50 subtype in SCAN-B, which has no purity estimate); open, m1 ",
    "with the proliferation score added (m3). Readouts and scores are ",
    "rank-normalised within each cohort, so coefficients are in SD per SD. ",
    "Bars, 95% confidence intervals. The OXPHOS coefficient is expected to be ",
    "positive by construction for the signatures, the arms and the six ",
    "signature genes, whose genes and signs were selected by association with ",
    "OXPHOS in these cohorts; S5 is additionally selected against association ",
    "with MYC. No correction for multiple genes is applied."),
  B = paste0(
    "MYC with OXPHOS held fixed, for each readout of Display A, on four ",
    "estimators of MYC activity ordered by their overlap with proliferation ",
    "(MYC_UP.V1_UP, 1.5%; M_a, 14.8%; HALLMARK MYC targets, 23.5%; the ",
    "CollecTRI MYC regulon), in each cohort and on each OXPHOS score. Filled, ",
    "m1; open, m3. Faded points are m3 cells in which the estimator and the ",
    "proliferation score correlate at 0.80 or more; they are shown but not ",
    "interpreted. Bars, 95% confidence intervals."),
  C = paste0(
    "OXPHOS coefficients with MYC activity held fixed (m1, M_a) for BCL2L1 ",
    "and BBC3, with model-based 95% confidence intervals, and their ",
    "difference with a 95% percentile interval from 2,000 paired bootstrap ",
    "resamples of tumours, in each cohort and on each OXPHOS score."))

# =============================================================================
# 14. Save, per R rules 4 and 5
# =============================================================================
message("\n14. saving")

ASSERTIONS <- dplyr::bind_rows(ASSERT)
.pr(ASSERTIONS %>% dplyr::transmute(check = substr(check, 1, 86), n_compared,
                                    max_abs_dev = signif(max_abs_dev, 3), pass))
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")

FRAME <- dplyr::bind_rows(DF) %>% dplyr::mutate(PAM50 = as.character(PAM50))
readr::write_csv(RES,     file.path(DIR_TABLES, "E47_coefficients_OX_MYC.csv"))
readr::write_csv(COEF,    file.path(DIR_TABLES, "E47_coefficient_matrix.csv"))
readr::write_csv(LABELS,  file.path(DIR_TABLES, "E47_labels.csv"))
readr::write_csv(SEC9,    file.path(DIR_TABLES, "E47_bcl2l1_bbc3.csv"))
readr::write_csv(D1,      file.path(DIR_TABLES, "E47_d1_leakage.csv"))
readr::write_csv(G3,      file.path(DIR_TABLES, "E47_g3_identifiability.csv"))
readr::write_csv(ASSERTIONS, file.path(DIR_TABLES, "E47_assertions.csv"))

saveRDS(list(
  verdict          = VERDICT,
  labels           = LABELS,
  primary_cells    = PRIMARY_CELLS,
  ladder           = LADDER,
  results          = RES,
  ratio            = RATIO,
  d1_leakage       = D1,
  d2_rungs         = D2,
  bcl2l1_bbc3      = SEC9,
  steeper_licensed = STEEPER_LICENSED,
  boot_replicates  = BOOT,
  coefficient_matrix = COEF,
  vcov             = VCOV,
  vif              = VIF,
  fit_digests      = FIT_DIGESTS,
  frame            = FRAME,
  gates            = list(rules = CC, margins = MARGINS, g1 = G1, g2 = G2,
                          g3 = G3, g4 = G4),
  displays         = list(A = PA, B = PB, C_coefficients = PC1,
                          C_difference = PC2),
  legend_draft     = LEGEND,
  assertions       = ASSERTIONS,
  # FILE NAMES, NOT PATHS, so this element digests the same in a redirected
  # dry run and in the author's run. Figures are not digested.
  figure_files     = list(
    outputs = paste0(rep(c(FIG_A, FIG_B, FIG_C), each = 2), c(".png", ".pdf")),
    tracked_in_docs_figures = paste0(FIG_DATE, "_",
                                     rep(c(FIG_A, FIG_B, FIG_C), each = 2),
                                     c(".png", ".pdf"))),
  spec = list(
    declaration = paste0(
      "docs/2026-10-22_E47_mutual_adjustment_declaration.md, committed alone ",
      "at dcb58a7; addendum 4a.7 at 3b5e221. Both before this script."),
    posture = "EXPLORATORY and POST-HOC (section 1). Main effects only.",
    model = paste0("readout ~ OX + MYC + spine, OLS, t-based 95% CI. m0 none; ",
                   "m1 PAM50 + purity + leukocyte (TCGA) or PAM50 (SCAN-B); ",
                   "m3 m1 + PROLIF_DISJOINT. m0 reported, never read."),
    verdict_rule = paste0(
      "Section 7, transcribed in .label(): on S6 and M_a over the non-void ",
      "m1/m3 cells - MYC ADDS / ESTIMATOR-DEPENDENT / INSTRUMENT-SPLIT / ",
      "OXPHOS ONLY / MIXED, first match wins. S5 and S7 labelled ",
      "descriptively. beta_OX|MYC never read."),
    s5_bias = paste0(
      "4a.7: S5 is selected against MYC association, so its beta_MYC|OX is ",
      "pushed towards zero by construction. An S5-S6 disagreement means ",
      "'MYC's contribution is carried by BID'."),
    ratio = paste0(
      "beta_MYC/beta_OX is a description only; beta_OX is inflated by ",
      "construction, so the ratio understates MYC's share."),
    sec9 = paste0("Paired tumour-level bootstrap, ", BOOT_N, " resamples, seed ",
                  "PROJECT_SEED + cell index; percentile 95% CI. Licensed only ",
                  "if delta > 0 with CI above zero in all four cells."),
    snapshot = "data/from_validation/ and E02's scores; joint snapshot not read.",
    not_digested = paste0(
      "$fits_S6 holds fitted lm objects, which are not digest-stable across ",
      "runs; $fit_digests (coef and vcov per fit) is the stable handle."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  fits_S6          = FITS_S6,
  built            = Sys.time()), PATH_E47)

message("   ", PATH_E47)
message("\nE47 done. VERDICT (S6): ", VERDICT,
        ". 'More steeply' (section 9) licensed: ", STEEPER_LICENSED, ".")
message("beta_OX|MYC is expected by construction and was not read.\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E47)

  # the verdict, and the variants' descriptive labels
  x$verdict
  x$labels %>% as.data.frame()

  # the eight primary cells: S6, M_a, m1 and m3
  x$primary_cells %>% dplyr::transmute(cohort, instrument, rung,
      b_MYC = round(b_MYC, 3), myc_lo = round(myc_lo, 3),
      myc_hi = round(myc_hi, 3), b_OX = round(b_OX, 3), void) %>%
    as.data.frame()

  # N9: does the MYC answer depend on how MYC activity is measured?
  x$ladder %>% as.data.frame()

  # D1: how much of Fig. 4g's MYC step is leftover OXPHOS
  x$d1_leakage %>% as.data.frame()

  # section 9: BCL2L1 against BBC3
  x$bcl2l1_bbc3 %>% as.data.frame()
  x$steeper_licensed

  # every gene, both coefficients, M_a, m1
  x$results %>% dplyr::filter(estimator == "M_a", rung == "m1") %>%
    dplyr::transmute(readout, cohort, instrument, b_MYC = round(b_MYC, 3),
                     b_OX = round(b_OX, 3)) %>% as.data.frame()

  x$assertions %>% as.data.frame()
  cat(x$legend_draft$A, "\n")
}
