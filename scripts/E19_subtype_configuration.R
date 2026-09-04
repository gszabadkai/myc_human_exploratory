# E19_subtype_configuration.R
# =============================================================================
# IS THE CONFIGURATION COMPOSITIONAL? The twelve genes inside Luminal and Basal.
#
# EXPLORATORY. Nothing here is pre-registered. The four branches in section 0b
# are fixed BEFORE any number is computed and are not re-decided after.
#
# =============================================================================
# WHY THIS EXISTS
# =============================================================================
# E18 established that BCL2L1 is the one gene in the 3.1 configuration that does
# not transfer to cell lines: +0.31 to +0.43 on 8 of 8 tumour cells, -0.12 to
# +0.01 on 8 of 8 CCLE cells, while BBC3, MCL1, BAD and BIK all keep their sign
# on 4 of 4 rulers against both cohorts.
#
# The B4 addendum eliminated purity (E16 check 4: BCL2L1 +0.419 -> +0.406 with
# purity and leukocyte fraction on 1,007 TCGA samples). SUBTYPE IS THE LAST
# CHEAP COMPOSITIONAL EXPLANATION and it has never been tested for this gene.
#
# Synthesis 3.7 carries Basal values for BBC3 and MCL1 because both appear in
# the three surviving ratios. BCL2L1 appears in none of them, so it has no
# within-subtype value anywhere in the repo. That is the gap this script fills.
#
# =============================================================================
# THE D3/S1 WARNING, WHICH IS WHY THE POOLED COLUMN STAYS IN
# =============================================================================
# E10 section 5.1: BCL2 against MYC is -0.369 pooled and -0.009 inside LumA,
# because the pooled value was reading the difference BETWEEN subtypes rather
# than anything within one. A pooled column that disagrees with BOTH of its
# strata is that same artefact, and it is only visible when all three are side
# by side. E10 ran that test for the 35 RATIOS and never for the 12 GENES.
#
# =============================================================================
# NO RE-SCORING WITHIN STRATUM. THE MOST LIKELY SILENT ERROR HERE.
# =============================================================================
# GSVA is cohort-relative and so is every z-score, so re-scoring inside Basal
# builds a DIFFERENT AXIS and answers a different question. The pooled scores
# are subsetted, exactly as E10's `.stratum_cors` does it. Section 5's
# reproduction check against E10's own `component_strata` is what proves this
# script did not quietly re-score.
#
# =============================================================================
# n AND POWER, FIXED NOW AND NOT ARGUED WITH LATER
# =============================================================================
# Basal is 171 TCGA and 317 SCAN-B. E10 records that a 171-sample stratum gives
# a 95% interval about +/- 0.15 wide on rho, so a Basal-versus-Luminal
# difference smaller than about 0.2 is not separable from sampling in TCGA.
# A BASAL NULL IS WEAK EVIDENCE AND IS WRITTEN AS SUCH. The informative
# comparison is BETWEEN strata, never against zero.
#
# ALL TWELVE ARE REPORTED. Selecting BCL2L1 after seeing them is the
# grid-of-cells trap. BBC3 and MCL1 are the internal reference precisely because
# synthesis 3.7 already carries their Basal values - a mismatch there is a bug,
# not a finding.
#
# =============================================================================
# THE MOUSE NUMBERS IN SECTION 9 ARE TRANSCRIBED, NOT COMPUTED
# =============================================================================
# myc_mouse is NOT attached and is never written to. Its per-gene slopes are
# transcribed here from the pinned ref, the way E16's header transcribes the
# mouse recipe, so this script has no runtime dependency on an unattached repo.
# Source and how to re-read it are given in section 9. NO CROSS-SPECIES VALUE
# COMPARISON IS MADE ANYWHERE: signs and orderings only. Twenty-four animals
# against 1,095 tumours is not a comparison of magnitudes, and a species is a
# cohort.
#
# SCALE: the genes are log2(linear DESeq2-normalised + 1); ox_gsva and ox_ppd
# are READ from their saved objects; ox_lvl and ox_rel are built by E16's
# recipe on the log matrix. Ranks make the mixture safe.
# N3: never "primed", of a transcript, a cell line or an animal.
# SPECIES: human throughout. No ortholog function is called.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))
source(here::here("functions", "strata.R"))

message("\nE19: the configuration inside Luminal and Basal\n", strrep("=", 78))

PATH_E19      <- file.path(DIR_RESULTS, "subtype_configuration.rds")
PATH_E19_CFG  <- file.path(DIR_TABLES,  "E19_subtype_configuration.csv")
PATH_E19_BTW  <- file.path(DIR_TABLES,  "E19_between_subtype.csv")

# =============================================================================
# 0. Constants
# =============================================================================
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
PRIMING_ALL  <- sort(unique(c(PRIMING_PRO, PRIMING_ANTI)))

STRATA_PRIMING <- c("all", "Luminal", "Basal")
MIN_STRATUM_N  <- 30L
RULERS         <- c("ox_gsva", "ox_ppd", "ox_lvl", "ox_rel")
ARM_PRIMARY    <- "OXPHOS subunits"
PROLIF_COV     <- "PROLIF_DISJOINT"

# E16's load-bearing quantity.
GAP_NUM <- "BCL2L1"; GAP_DEN <- "BBC3"
FOCUS   <- "BCL2L1"                       # the gene the branches are about
REFERENCE <- c("BBC3", "MCL1")            # the two with known Basal behaviour

# CCLE, section 8. Optional and SECONDARY - it only matters if section 7 gives
# branch 2 or 3, and it must never hold task 1 hostage to a symlink.
CCLE_MIN_LINES <- 25L

# =============================================================================
# 0b. THE FOUR BRANCHES. Fixed here, before any number.
# =============================================================================
# "Holds in a stratum" means the gene keeps the SIGN it has pooled, on at least
# 3 of the 4 rulers, in BOTH cohorts. It is a statement about sign, never about
# magnitude: at n = 171 a Basal magnitude is not separable from its pooled one.
HOLD_RULE <- paste0(
  "a gene HOLDS in a stratum if it keeps its pooled sign on >= 3 of 4 rulers ",
  "in BOTH cohorts. Sign only - a Basal magnitude is not separable at n = 171.")

BRANCHES <- tibble::tibble(
  branch = 1:4,
  holds_in = c("both strata", "Luminal only", "Basal only", "neither"),
  reading = c(
    paste("NOT COMPOSITIONAL. Tier 1 survives; the CCLE discrepancy stays",
          "unexplained and gets harder."),
    paste("THE DAMAGING BRANCH. MMTV-Myc is basal-like, so the mouse",
          "comparison belongs in Basal, where the association would not hold.",
          "It also explains CCLE, whose breast panel is basal-enriched.",
          "Tier 1 needs restating."),
    paste("The mouse comparison is valid, but CCLE's failure becomes MORE",
          "puzzling, since CCLE breast is basal-enriched."),
    "PURELY COMPOSITIONAL. Worst case for Tier 1."))

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
e10  <- readRDS(file.path(DIR_RESULTS, "machinery_and_priming.rds"))
fr   <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames

# DEPENDENCY GUARD. Section 5's reproduction check is the only thing tying this
# script's rulers to the ones already in the repo, so its input is not optional.
if (!"component_strata" %in% names(e10)) {
  stop("results/machinery_and_priming.rds carries no `component_strata`; ",
       "section 5's reproduction check cannot run and nothing below would be ",
       "known to be on E10's axes. Re-run E10.", call. = FALSE)
}

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)

tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
# log2(x + 1) is monotone, so every Spearman below is identical to one on the
# raw linear matrix. E10 asserts that; it is not re-asserted here.
GT <- log2(tcga_lin$mat[, ID_T, drop = FALSE] + 1)
GS <- log2(scanb_lin$mat[, ID_S, drop = FALSE] + 1)
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

RES_T <- .symbol_resolver(rownames(GT), NULL)
RES_S <- .symbol_resolver(rownames(GS), sc$symbol_map)

CV_T <- t(mito$gsva_cov[PROLIF_COV, ID_T, drop = FALSE])
CV_S <- t(sc$gsva_cov[PROLIF_COV, ID_S, drop = FALSE])
stopifnot(identical(rownames(CV_T), ID_T), identical(rownames(CV_S), ID_S))

COH <- list(
  TCGA     = list(G = GT, res = RES_T, ids = ID_T, cov = CV_T, arms = mito),
  `SCAN-B` = list(G = GS, res = RES_S, ids = ID_S, cov = CV_S, arms = sc))
message("   TCGA ", length(ID_T), " samples | SCAN-B ", length(ID_S),
        " samples")

OX_SUB   <- mito$arm_sets[[ARM_PRIMARY]]
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL
REST     <- setdiff(MITO_ALL, OX_SUB)
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL, length(OX_SUB) == 89L,
          !any(grepl("^MT-", OX_SUB)), all(OX_SUB %in% MITO_ALL))

# =============================================================================
# 2. The four rulers - E16 section 2, verbatim
# =============================================================================
# ox_ppd and ox_gsva are READ from their saved objects and are NOT rebuilt.
# ox_lvl and ox_rel are the mouse-recipe composites: mean per-gene z across
# samples, on the log matrix.
message("\n2. the four rulers")

.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res)
  M  <- gr$mat
  v  <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}

rulers <- lapply(COH, function(C) {
  M <- rbind(
    ox_gsva = as.numeric(C$arms$gsva_arms[ARM_PRIMARY, C$ids]),
    ox_ppd  = as.numeric(C$arms$mitopps_arms[ARM_PRIMARY, C$ids]),
    ox_lvl  = .comp(OX_SUB, C),
    ox_rel  = .comp(OX_SUB, C) - .comp(REST, C))
  colnames(M) <- C$ids
  stopifnot(!anyNA(M), identical(rownames(M), RULERS))
  M
})
message("   4 rulers x 2 cohorts built; ox_gsva and ox_ppd read, not rebuilt")

# =============================================================================
# 3. The strata
# =============================================================================
message("\n3. the strata")

STR <- list(TCGA = .build_strata(fr, "TCGA", ID_T),
            `SCAN-B` = .build_strata(fr, "SCAN-B", ID_S))
strata_n <- tibble::tibble(
  stratum = STRATA_PRIMING,
  TCGA = lengths(STR$TCGA[STRATA_PRIMING]),
  `SCAN-B` = lengths(STR$`SCAN-B`[STRATA_PRIMING]))
strata_n %>% as.data.frame() %>% print(row.names = FALSE)
# A stratum that silently fell below the floor would drop rows and turn a
# missing cell into a missing gene three sections downstream.
tiny <- strata_n$stratum[pmin(strata_n$TCGA, strata_n$`SCAN-B`) < MIN_STRATUM_N]
if (length(tiny)) {
  stop("stratum/strata below the ", MIN_STRATUM_N, " floor: ",
       paste(tiny, collapse = ", "), call. = FALSE)
}

GR <- lapply(COH, function(C) .gene_rows(PRIMING_ALL, C$G, C$res))
miss <- unique(unlist(lapply(GR, function(g) g$missing), use.names = FALSE))
if (length(miss)) {
  stop("these priming genes did not resolve to exactly one matrix row: ",
       paste(miss, collapse = ", "),
       "\nCheck the SCAN-B symbol map before going further.", call. = FALSE)
}
message("   all ", length(PRIMING_ALL), " priming genes present in both cohorts")

# =============================================================================
# 4. The 12 x 4 x 2 x 3 block
# =============================================================================
# `.cor_block` is E16's and E10's, copied verbatim: partial Spearman on ranks
# with Fisher-z intervals and the Bonett-Wright variance
# se = sqrt((1 + rho^2/2)/(n - 3 - k)). The plain 1/(n - 3) is the Pearson case
# and understates a rank correlation. Section 5 is what proves the copy has not
# drifted.
message("\n4. the twelve, on four rulers, in three strata")

.rank_rows <- function(M) {
  out <- t(apply(M, 1L, rank))
  if (nrow(M) == 1L) out <- matrix(out, nrow = 1L)
  dimnames(out) <- dimnames(M)
  out
}
.cor_block <- function(A, B, cov = NULL) {
  n <- ncol(A)
  stopifnot(identical(colnames(A), colnames(B)))
  RA <- .rank_rows(A); RB <- .rank_rows(B)
  k  <- 0L
  if (!is.null(cov)) {
    stopifnot(nrow(cov) == n)
    H  <- qr(cbind(`(Intercept)` = 1, apply(cov, 2L, rank)))
    k  <- ncol(cov)
    RA <- RA - t(qr.fitted(H, t(RA)))
    RB <- RB - t(qr.fitted(H, t(RB)))
  }
  R  <- suppressWarnings(stats::cor(t(RA), t(RB)))
  z  <- atanh(pmin(pmax(R, -0.999999999), 0.999999999))
  se <- sqrt((1 + R^2 / 2) / (n - 3 - k))
  tibble::tibble(
    ruler = rep(rownames(R), times = ncol(R)),
    gene  = rep(colnames(R), each  = nrow(R)),
    n = n, k_cov = k,
    rho   = as.vector(R),
    ci_lo = as.vector(tanh(z - 1.959964 * se)),
    ci_hi = as.vector(tanh(z + 1.959964 * se)))
}

# E10's `.stratum_cors`, generalised to this script's ruler matrix. THE POOLED
# SCORES ARE SUBSET, NEVER RE-SCORED - see the header.
.stratum_cors <- function(coh, st) {
  ids <- STR[[coh]][[st]]
  if (is.null(ids) || length(ids) < MIN_STRATUM_N) return(NULL)
  .cor_block(rulers[[coh]][, ids, drop = FALSE],
             GR[[coh]]$mat[PRIMING_ALL, ids, drop = FALSE],
             cov = COH[[coh]]$cov[ids, , drop = FALSE]) %>%
    dplyr::mutate(cohort = coh, stratum = st)
}

config <- dplyr::bind_rows(lapply(names(COH), function(coh)
  dplyr::bind_rows(lapply(STRATA_PRIMING, .stratum_cors, coh = coh)))) %>%
  dplyr::mutate(
    side = dplyr::if_else(gene %in% PRIMING_PRO, "pro-apoptotic",
                          "anti-apoptotic"),
    adjustment = PROLIF_COV,
    ci_excludes_0 = (ci_lo > 0) | (ci_hi < 0))
stopifnot(nrow(config) ==
            length(PRIMING_ALL) * length(RULERS) * 2L * length(STRATA_PRIMING))
message("   ", nrow(config), " rows = 12 genes x 4 rulers x 2 cohorts x 3 strata")

for (coh in names(COH)) {
  message("\n   ", coh, ", rho by stratum (ox_gsva | ox_ppd | ox_lvl | ox_rel):")
  config %>%
    dplyr::filter(cohort == coh) %>%
    dplyr::mutate(rho = round(rho, 3),
                  stratum = factor(stratum, levels = STRATA_PRIMING),
                  ruler = factor(ruler, levels = RULERS)) %>%
    dplyr::select(gene, side, stratum, ruler, rho) %>%
    tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
    dplyr::arrange(side, gene, stratum) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

# =============================================================================
# 5. THE REPRODUCTION CHECK. Not a result.
# =============================================================================
# E10's `component_strata` already carries this block for two of the four
# rulers - `OXPHOS` is ox_gsva and `OXPHOS_mitopps` is ox_ppd - on the same
# genes, the same strata and the same covariate. So those 144 cells MUST come
# back bit-equal. A failure here is a broken input or a re-scored axis, not a
# new result, and it is the only thing tying the two new rulers to the repo's
# existing numbers.
message("\n5. reproduction check against E10 (not a result)")

E10_MAP <- c(ox_gsva = "OXPHOS", ox_ppd = "OXPHOS_mitopps")
repro <- config %>%
  dplyr::filter(ruler %in% names(E10_MAP)) %>%
  dplyr::mutate(axis = unname(E10_MAP[ruler])) %>%
  dplyr::select(cohort, stratum, axis, gene, n, rho_new = rho,
                lo_new = ci_lo, hi_new = ci_hi) %>%
  dplyr::inner_join(
    e10$component_strata %>%
      dplyr::select(cohort, stratum, axis, gene, n_e10 = n, rho_e10 = rho,
                    lo_e10 = ci_lo, hi_e10 = ci_hi),
    by = c("cohort", "stratum", "axis", "gene")) %>%
  dplyr::mutate(d_rho = abs(rho_new - rho_e10),
                d_lo = abs(lo_new - lo_e10), d_hi = abs(hi_new - hi_e10))
stopifnot(nrow(repro) == length(PRIMING_ALL) * 2L * 2L * length(STRATA_PRIMING))
if (!all(repro$n == repro$n_e10) ||
    max(c(repro$d_rho, repro$d_lo, repro$d_hi)) > 1e-12) {
  bad <- repro %>% dplyr::arrange(dplyr::desc(d_rho)) %>% utils::head(5)
  stop("E19 does not reproduce E10's `component_strata`. Largest |d_rho| = ",
       signif(max(repro$d_rho), 4), " on ", bad$cohort[1], " / ",
       bad$stratum[1], " / ", bad$axis[1], " / ", bad$gene[1],
       ". Either an axis was re-scored inside the stratum, or the covariate ",
       "or the strata definition has moved. Do NOT proceed.", call. = FALSE)
}
message("   all ", nrow(repro), " shared cells reproduce E10 to 1e-12")

# =============================================================================
# 6. THE BETWEEN-SUBTYPE TEST, per gene
# =============================================================================
# E10's `$between_test` logic, applied per GENE rather than per ratio. A pooled
# value OUTSIDE the range of its own two strata is the D3/S1 signature: a
# difference between subtypes being read as an association within one.
message("\n6. the between-subtype test, per gene")

between_gene <- config %>%
  dplyr::select(cohort, ruler, gene, side, stratum, rho) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho) %>%
  dplyr::mutate(pooled_outside = all < pmin(Luminal, Basal) |
                                 all > pmax(Luminal, Basal),
                lum_minus_basal = Luminal - Basal,
                # HOW FAR outside. The flag is a SIGN test with no magnitude
                # threshold, so it fires whenever the pooled value misses a
                # narrow stratum range by any amount at all - including by
                # sampling noise. E10's motivating case, BCL2 against MYC, was
                # -0.369 pooled against -0.009 inside LumA: an excursion of
                # about 0.36. Without this column the flag cannot be told apart
                # from that, and it is what makes the flag readable.
                excursion = dplyr::case_when(
                  all < pmin(Luminal, Basal) ~ pmin(Luminal, Basal) - all,
                  all > pmax(Luminal, Basal) ~ all - pmax(Luminal, Basal),
                  TRUE                       ~ 0))
D3S1_EXCURSION <- 0.36    # BCL2 against MYC, E10 section 5.1. The comparator.

message("\n   how many of the 12 have the pooled value outside both strata?")
between_gene %>%
  dplyr::group_by(cohort, ruler) %>%
  dplyr::summarise(n_genes = dplyr::n(),
                   n_pooled_outside = sum(pooled_outside), .groups = "drop") %>%
  dplyr::mutate(ruler = factor(ruler, levels = RULERS)) %>%
  dplyr::arrange(cohort, ruler) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   HOW BIG are the excursions? The flag is a sign test; this is",
        " the magnitude.")
message("   E10's motivating case (BCL2 vs MYC, -0.369 pooled against -0.009",
        " in LumA) is about ", D3S1_EXCURSION, ".")
message("   largest excursion anywhere here: ",
        round(max(between_gene$excursion), 3), " (",
        with(between_gene, paste(gene[which.max(excursion)],
                                 cohort[which.max(excursion)],
                                 ruler[which.max(excursion)])), ")")
if (max(between_gene$excursion) < D3S1_EXCURSION / 3) {
  message("   -> EVERY excursion is under a third of the D3/S1 case. On these",
          " numbers the flag\n      separates almost nothing: it is firing on",
          " narrow stratum ranges, not on a\n      between-subtype effect.",
          " Read the magnitude, never the flag alone.")
}

message("\n   the genes flagged in BOTH cohorts on at least 3 of 4 rulers:")
flagged <- between_gene %>%
  dplyr::group_by(cohort, gene) %>%
  dplyr::summarise(n_flag = sum(pooled_outside), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cohort, values_from = n_flag,
                     names_prefix = "flags_") %>%
  dplyr::mutate(both = flags_TCGA >= 3L & `flags_SCAN-B` >= 3L)
flagged %>% dplyr::arrange(dplyr::desc(both), gene) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   ", FOCUS, " and the two reference genes, all three strata:")
between_gene %>%
  dplyr::filter(gene %in% c(FOCUS, REFERENCE)) %>%
  dplyr::mutate(ruler = factor(ruler, levels = RULERS),
                dplyr::across(c(all, Luminal, Basal, lum_minus_basal),
                              ~ round(.x, 3))) %>%
  dplyr::arrange(gene, cohort, ruler) %>%
  dplyr::mutate(excursion = round(excursion, 3)) %>%
  dplyr::select(gene, cohort, ruler, all, Luminal, Basal, lum_minus_basal,
                pooled_outside, excursion) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 6.1 the gap, within each stratum ---------------------------------------
# BCL2L1 minus BBC3, as a difference of partial rhos - exactly how E16's
# `the_gap` and E18's `ccle_gap` are built.
message("\n6.1 the gap, ", GAP_NUM, " minus ", GAP_DEN, ", by stratum:")
gap_strata <- config %>%
  dplyr::filter(gene %in% c(GAP_NUM, GAP_DEN)) %>%
  dplyr::select(cohort, stratum, ruler, gene, rho) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rho) %>%
  dplyr::mutate(gap = .data[[GAP_NUM]] - .data[[GAP_DEN]])
gap_strata %>%
  dplyr::mutate(stratum = factor(stratum, levels = STRATA_PRIMING),
                ruler = factor(ruler, levels = RULERS),
                gap = round(gap, 3)) %>%
  dplyr::select(cohort, stratum, ruler, gap) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = gap) %>%
  dplyr::arrange(cohort, stratum) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. THE BRANCH, on the rule fixed in section 0b
# =============================================================================
message("\n7. the branch\n", strrep("-", 78))
message("HOLD RULE (fixed before any number): ", HOLD_RULE)

.holds <- function(g, st) {
  d <- config %>% dplyr::filter(gene == g, stratum %in% c("all", st))
  pooled <- d %>% dplyr::filter(stratum == "all") %>%
    dplyr::select(cohort, ruler, rho_all = rho)
  d %>% dplyr::filter(stratum == st) %>%
    dplyr::select(cohort, ruler, rho) %>%
    dplyr::inner_join(pooled, by = c("cohort", "ruler")) %>%
    dplyr::mutate(same_sign = sign(rho) == sign(rho_all)) %>%
    dplyr::group_by(cohort) %>%
    dplyr::summarise(n_same = sum(same_sign), .groups = "drop") %>%
    dplyr::summarise(holds = all(n_same >= 3L),
                     detail = paste(sprintf("%s %d/4", cohort, n_same),
                                    collapse = ", ")) %>%
    dplyr::mutate(gene = g, stratum = st)
}

hold_tab <- dplyr::bind_rows(lapply(PRIMING_ALL, function(g)
  dplyr::bind_rows(lapply(c("Luminal", "Basal"), .holds, g = g)))) %>%
  dplyr::select(gene, stratum, holds, detail)
hold_wide <- hold_tab %>%
  dplyr::select(gene, stratum, holds) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = holds) %>%
  dplyr::mutate(branch = dplyr::case_when(
    Luminal &  Basal ~ 1L,
    Luminal & !Basal ~ 2L,
    !Luminal & Basal ~ 3L,
    TRUE             ~ 4L))
message("\n   every gene, which strata it holds in:")
hold_wide %>%
  dplyr::left_join(dplyr::distinct(dplyr::select(config, gene, side)),
                   by = "gene") %>%
  dplyr::arrange(branch, side, gene) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   ", FOCUS, " and the two reference genes, cohort by cohort:")
hold_tab %>% dplyr::filter(gene %in% c(FOCUS, REFERENCE)) %>%
  as.data.frame() %>% print(row.names = FALSE)

focus_branch <- hold_wide$branch[hold_wide$gene == FOCUS]
branch_row   <- BRANCHES[BRANCHES$branch == focus_branch, ]
message("\n   ", FOCUS, " -> BRANCH ", focus_branch, ": holds in ",
        branch_row$holds_in)
message("   ", branch_row$reading)

# The internal reference. A mismatch here is a bug, not a finding - synthesis
# 3.7 already carries BBC3 and MCL1 in Basal and section 5 has already asserted
# the ox_gsva cells bit-equal, so this is a readability check on top of that.
ref_ok <- all(hold_wide$Luminal[hold_wide$gene %in% REFERENCE] &
              hold_wide$Basal[hold_wide$gene %in% REFERENCE])
message("\n   internal reference (", paste(REFERENCE, collapse = " and "),
        ") hold in both strata: ", ref_ok,
        if (!ref_ok) " <- CHECK THIS, 3.7 says they should" else "")

# =============================================================================
# 8. CCLE breast composition. SECONDARY, and it is allowed to be absent.
# =============================================================================
# Only informative if section 7 gives branch 2 or 3. It is a COMPOSITION COUNT,
# not a classification: DepMap's `ModelSubtypeFeatures` is a curator free-text
# field with mixed conventions, and NO PAM50 CLASSIFIER IS BUILT HERE.
#
# Unlike E17, a missing symlink SKIPS rather than stops, because task 1 does not
# depend on it and a stop would cost the whole script. The skip is RECORDED.
message("\n8. CCLE breast composition (secondary)")

DIR_DEPMAP <- here::here("data", "raw", "depmap")
PATH_MODEL <- file.path(DIR_DEPMAP, "Model.csv")
PATH_EXPR  <- file.path(DIR_DEPMAP,
                        "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
ccle <- NULL; ccle_skip <- NULL

.link_ok <- (nzchar(Sys.readlink(DIR_DEPMAP)) || dir.exists(DIR_DEPMAP)) &&
  (!nzchar(Sys.readlink(DIR_DEPMAP)) || dir.exists(Sys.readlink(DIR_DEPMAP)))
if (!.link_ok || !file.exists(PATH_MODEL) || !file.exists(PATH_EXPR)) {
  ccle_skip <- paste0(
    "data/raw/depmap is not readable, so the CCLE composition count was NOT ",
    "made. This is a SKIP, not a null. See data/raw/depmap_README.md - the ",
    "symlink points into myc_human_validation's gitignored data/raw/ and this ",
    "repo cannot regenerate what it points at.")
  message("   SKIPPED - ", ccle_skip)
} else {
  MODEL <- data.table::fread(PATH_MODEL, data.table = FALSE)
  lin_col <- if ("OncotreeLineage" %in% names(MODEL)) "OncotreeLineage" else
    if ("lineage" %in% names(MODEL)) "lineage" else
      stop("Model.csv has no lineage column.", call. = FALSE)
  sub_col <- "ModelSubtypeFeatures"
  if (!sub_col %in% names(MODEL)) {
    ccle_skip <- paste0("Model.csv has no `", sub_col, "` column, so no ",
                        "subtype annotation is available without building a ",
                        "classifier - which this session does not do.")
    message("   SKIPPED - ", ccle_skip)
  } else {
    # Only the ModelID column of the expression file is read: this is a count,
    # and the 71 are defined the same way E17 and E18 define them.
    ex <- data.table::fread(PATH_EXPR, select = c("ModelID",
                                                  "IsDefaultEntryForModel"),
                            data.table = FALSE)
    keep <- tolower(trimws(as.character(ex$IsDefaultEntryForModel))) %in%
      c("yes", "true", "t", "1")
    expr_ids <- unique(as.character(ex$ModelID[keep]))

    b <- MODEL[!is.na(MODEL[[lin_col]]) & MODEL[[lin_col]] == "Breast", ]
    b$has_expr <- b$ModelID %in% expr_ids
    txt <- tolower(trimws(ifelse(is.na(b[[sub_col]]), "", b[[sub_col]])))
    b$call <- dplyr::case_when(
      grepl("basal", txt)   & !grepl("luminal", txt) ~ "basal",
      grepl("luminal", txt) & !grepl("basal", txt)   ~ "luminal",
      grepl("basal", txt)   &  grepl("luminal", txt) ~ "ambiguous",
      nzchar(txt)                                    ~ "unlabelled call",
      TRUE                                           ~ "no annotation")
    ccle_counts <- b %>%
      dplyr::filter(has_expr) %>%
      dplyr::count(call, name = "n_lines") %>%
      dplyr::arrange(dplyr::desc(n_lines))
    message("   ", sum(b$has_expr), " Breast lines with expression, by ",
            sub_col, " free text:")
    ccle_counts %>% as.data.frame() %>% print(row.names = FALSE)
    message("\n   the raw free-text values, so the classification is auditable:")
    b %>% dplyr::filter(has_expr) %>%
      dplyr::count(call, value = .data[[sub_col]], name = "n") %>%
      dplyr::arrange(call, dplyr::desc(n)) %>%
      as.data.frame() %>% print(row.names = FALSE)

    big <- ccle_counts$call[ccle_counts$call %in% c("basal", "luminal") &
                              ccle_counts$n_lines >= CCLE_MIN_LINES]
    ccle <- list(counts = ccle_counts,
                 n_expr = sum(b$has_expr),
                 strata_clearing_floor = big,
                 min_lines = CCLE_MIN_LINES,
                 annotation = paste(sub_col, "- a DepMap curator FREE-TEXT",
                                    "field with mixed conventions, not a",
                                    "molecular classifier. No PAM50 classifier",
                                    "was built."))

    # --- 8.1 the within-stratum fit, ONLY where the floor is cleared --------
    # The task's own rule: report the count, and fit BCL2L1 and MCL1 within a
    # stratum only if it carries >= CCLE_MIN_LINES. A stratum below the floor
    # is reported as below it and NOT fitted-and-caveated.
    #
    # SAME RULE AS SECTION 4: the 71-line breast scores are built once and the
    # stratum is a SUBSET of them. Re-scoring inside 27 lines would build a
    # different axis, and it would not be the axis E18 measured against.
    if (!length(big)) {
      message("\n   NO CCLE stratum clears ", CCLE_MIN_LINES, " lines, so no ",
              "within-stratum fit is made. Reported and STOPPED here rather ",
              "than fitted and caveated.")
    } else {
      message("\n8.1 CCLE within-stratum fit - ", paste(big, collapse = ", "),
              " clears ", CCLE_MIN_LINES, " lines")

      EXPR <- local({
        m <- data.table::fread(PATH_EXPR, data.table = FALSE)
        k <- tolower(trimws(as.character(m$IsDefaultEntryForModel))) %in%
          c("yes", "true", "t", "1")
        META <- c("V1", "", "ProfileID", "SequencingID", "ModelConditionID",
                  "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel")
        gc_ <- setdiff(names(m), META)
        M <- as.matrix(m[k, gc_, drop = FALSE])
        rownames(M) <- as.character(m$ModelID[k])
        colnames(M) <- trimws(sub("\\s*\\(\\d+\\)$", "", colnames(M)))
        M[, !duplicated(colnames(M)), drop = FALSE]
      })
      if (max(EXPR, na.rm = TRUE) > 40) {
        stop("CCLE expression is not log2(TPM+1).", call. = FALSE)
      }

      # E18 section 4, restricted to the breast cohort. ARM_SETS comes from the
      # snapshot, which E18 asserted its own MitoCarta rebuild set-identical to.
      ARM_SETS <- mito$arm_sets
      ids_b <- intersect(b$ModelID[b$has_expr], rownames(EXPR))
      E_LOG <- t(EXPR[ids_b, , drop = FALSE])
      E_LOG <- E_LOG[stats::complete.cases(E_LOG), , drop = FALSE]
      sets_def <- c(list(MYC = sd_$myc_sets[[MYC_REF]],
                         PROLIF = sd_$cov_sets[[PROLIF_COV]]), ARM_SETS)
      sets_g <- lapply(sets_def, function(g) intersect(g, rownames(E_LOG)))
      if (any(lengths(sets_g) / lengths(sets_def) < 0.80)) {
        stop("a gene set falls below 0.80 coverage in CCLE - a harmonisation ",
             "failure, not attrition.", call. = FALSE)
      }
      GSc <- GSVA::gsva(GSVA::gsvaParam(exprData = E_LOG, geneSets = sets_g,
                                        kcdf = "Gaussian", minSize = 3L,
                                        maxSize = Inf), verbose = FALSE)
      if (length(setdiff(names(sets_g), rownames(GSc)))) {
        stop("GSVA silently dropped a set in the CCLE block.", call. = FALSE)
      }
      E_LIN <- 2^E_LOG - 1; E_LIN[E_LIN < 0] <- 0
      Smat <- t(vapply(lapply(ARM_SETS, function(g) intersect(g, rownames(E_LOG))),
                       function(g) colMeans(E_LIN[g, , drop = FALSE]),
                       numeric(ncol(E_LIN))))
      if (any(!is.finite(Smat) | Smat <= 0)) {
        stop("a CCLE arm x line pathway mean is zero or non-finite; mitoPPS is ",
             "undefined. Do NOT add a floor silently.", call. = FALSE)
      }
      PPSc <- local({
        N <- ncol(Smat); P <- nrow(Smat); Bi <- 1 / Smat
        A <- (Smat %*% t(Bi)) / N
        o <- Smat * (((1 / A) %*% Bi) - Bi) / (P - 1)
        dimnames(o) <- dimnames(Smat); o
      })
      .compc <- function(g) {
        M <- E_LOG[intersect(g, rownames(E_LOG)), , drop = FALSE]
        v <- apply(M, 1L, stats::var)
        colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
      }
      .zc <- function(v) (v - mean(v)) / stats::sd(v)
      Rc <- rbind(ox_gsva = .zc(GSc[ARM_PRIMARY, ]),
                  ox_ppd  = .zc(PPSc[ARM_PRIMARY, ]),
                  ox_lvl  = .zc(.compc(ARM_SETS[[ARM_PRIMARY]])),
                  ox_rel  = .zc(.compc(ARM_SETS[[ARM_PRIMARY]]) -
                                  .compc(setdiff(MITO_ALL,
                                                 ARM_SETS[[ARM_PRIMARY]]))))
      colnames(Rc) <- ids_b
      PROLIFc <- .zc(GSc["PROLIF", ])

      call_of <- b$call[match(ids_b, b$ModelID)]
      ccle_fit <- dplyr::bind_rows(lapply(big, function(st) {
        ii <- ids_b[call_of == st]
        .cor_block(Rc[, ii, drop = FALSE],
                   t(EXPR[ii, c(FOCUS, "MCL1"), drop = FALSE]),
                   cov = matrix(PROLIFc[ii], ncol = 1L,
                                dimnames = list(ii, "PROLIF"))) %>%
          dplyr::mutate(stratum = st)
      })) %>%
        dplyr::mutate(ci_excludes_0 = (ci_lo > 0) | (ci_hi < 0))
      ccle$within_stratum <- ccle_fit
      ccle$not_fitted <- setdiff(c("basal", "luminal"), big)

      message("   scores built on all ", length(ids_b), " breast lines and ",
              "SUBSET to the stratum - not re-scored inside it")
      ccle_fit %>%
        dplyr::mutate(ruler = factor(ruler, levels = RULERS)) %>%
        dplyr::transmute(stratum, gene, ruler, n, rho = round(rho, 3),
                         ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                         ci_excludes_0) %>%
        dplyr::arrange(stratum, gene, ruler) %>%
        as.data.frame() %>% print(row.names = FALSE)
      message("\n   NOT FITTED, below the ", CCLE_MIN_LINES, "-line floor: ",
              paste(ccle$not_fitted, collapse = ", "),
              ". Reported as below the floor, NOT as null.")
      message("   n = ", min(ccle_fit$n), " gives a 95% interval about +/- 0.38",
              " wide, so this is the weakest\n   cell in the whole session and",
              " is written as such.")
    }
  }
}

# =============================================================================
# 9. THE MOUSE PER-GENE SLOPES. TRANSCRIBED, read-only, from the pinned ref.
# =============================================================================
# SOURCE, and how to re-read it without writing to that repo:
#   git -C /Users/gs/G/data/MK_myc_2022/myc_mouse show \
#       e348dd8:docs/2026-09-02_myc_oxphos_priming_gate_model.md
#
# WHAT IS AND IS NOT AVAILABLE AT THAT REF, checked 2026-09-04:
#   - `results/gate_model_verification.rds` is GITIGNORED upstream (.gitignore
#     lines 10-11 are `results/` and `outputs/`), so its `gene_gate` table -
#     which carries `bX_wt` and `bX_myc` for all 14 roster genes - CANNOT be
#     read at a pinned ref. Only what the doc prints is readable.
#   - the doc's section 3.3 prints per-gene `bMX` for all 14.
#   - the doc's section 7 prints per-gene WT and Myc+ OXPHOS slopes, epi + imm
#     adjusted, for Bcl2l1 and Mcl1 ONLY: "Both `Mcl1` slopes are negative in
#     absolute terms (`gene_gate`, epi + imm adjusted: WT -1.00, Myc+ -0.20);
#     the ratio rises because BCL-XL falls faster still (`Bcl2l1` WT +0.01,
#     Myc+ -0.66)."
#   - Bbc3's `bX_wt` is NOT printed anywhere in the doc and is therefore NOT
#     AVAILABLE. It is recorded as absent below, never inferred.
#
# CONSISTENCY CHECK, and it is the reason these transcribed numbers can be
# trusted: bX_myc should equal bX_wt + bMX. Bcl2l1: +0.01 + (-0.670) = -0.660
# against the printed -0.66. Mcl1: -1.00 + (+0.806) = -0.194 against -0.20.
# Both hold, so the two tables are on the same scale and the transcription is
# not mixing regimes. The check is asserted below, not just asserted here.
#
# NO CROSS-SPECIES VALUE COMPARISON. Signs and orderings only.
message("\n9. the mouse per-gene slopes (transcribed from the pinned ref)")

mouse_gene <- tibble::tribble(
  ~gene,     ~bX_wt, ~bX_myc, ~bMX,   ~p_bMX,  ~bMX_within_tp, ~p_within_tp,
  "Bbc3",    NA_real_, NA_real_, 0.746, 0.059,  NA_real_,       0.44,
  "Bcl2l1",  0.01,    -0.66,    -0.670, 0.0026, NA_real_,       0.0026,
  "Mcl1",   -1.00,    -0.20,     0.806, 0.018,  NA_real_,       0.027)
mouse_gene$bX_wt_available <- !is.na(mouse_gene$bX_wt)

chk <- mouse_gene %>% dplyr::filter(bX_wt_available) %>%
  dplyr::mutate(implied = bX_wt + bMX, d = abs(implied - bX_myc))
if (max(chk$d) > 0.011) {
  stop("the transcribed mouse numbers are not internally consistent: ",
       "bX_wt + bMX should equal bX_myc and the largest discrepancy is ",
       signif(max(chk$d), 3), ". Re-read the pinned doc before using them.",
       call. = FALSE)
}
message("   identity bX_myc = bX_wt + bMX holds on both available genes ",
        "(max |d| = ", signif(max(chk$d), 2), ", rounding only)")

mouse_composite <- tibble::tribble(
  ~endpoint,        ~definition,                  ~bX,     ~p_bX,  ~Mstar,
  "T (trigger)",    "log2(Bbc3) - log2(Bcl2l1)",  -0.006,  0.973,   0.007,
  "B (guardian)",   "log2(Mcl1) - log2(Bcl2l1)",  -0.541,  0.0040,  0.545)

mouse_gene %>%
  dplyr::select(gene, bX_wt, bX_myc, bMX, p_bMX, bX_wt_available) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   the composite endpoints, genotype estimator:")
mouse_composite %>% as.data.frame() %>% print(row.names = FALSE)

# THE QUESTION THIS ANSWERS, and it is answered on a SIGN, not a magnitude.
# The synthesis argues the human cohort sits below M*, so the human marginal
# should correspond to the mouse's slope at low MYC - which is bX_wt.
bx_bcl2l1  <- mouse_gene$bX_wt[mouse_gene$gene == "Bcl2l1"]
below_mstar_ok <- bx_bcl2l1 > 0.1     # "positive" in any usable sense
mouse_reading <- if (below_mstar_ok) {
  paste("bX_wt(Bcl2l1) is POSITIVE, so the human positive BCL2L1-OXPHOS",
        "association is what the mouse predicts below M*.")
} else {
  paste0(
    "bX_wt(Bcl2l1) = ", bx_bcl2l1, " - FLAT, not positive. The doc's own ",
    "words are 'a flat Bcl2l1'. AT LOW MYC THE MOUSE SAYS RESPIRATION DOES ",
    "NOTHING TO Bcl2l1, so the 'human sits below M*' argument does NOT ",
    "deliver the human association; it predicts approximately zero there. ",
    "What the two species DO agree on is the ORDERING of the two guardians ",
    "against respiration at low MYC - Bcl2l1 above Mcl1 in both - but in the ",
    "mouse that ordering is carried by Mcl1 FALLING (-1.00) and in the human ",
    "by BCL2L1 RISING. Same order, different limb.")
}
message("\n   ", mouse_reading)

# =============================================================================
# 10. Save
# =============================================================================
message("\n10. saving")

out <- list(
  configuration   = config,
  between_gene    = between_gene,
  flagged         = flagged,
  gap_strata      = gap_strata,
  holds           = hold_tab,
  holds_wide      = hold_wide,
  branches        = BRANCHES,
  focus_branch    = focus_branch,
  focus_reading   = branch_row$reading,
  reference_ok    = ref_ok,
  strata_n        = strata_n,
  repro           = repro,
  ccle            = ccle,
  ccle_skip       = ccle_skip,
  mouse_gene      = mouse_gene,
  mouse_composite = mouse_composite,
  mouse_reading   = mouse_reading,
  spec = list(
    what = paste("Is the 3.1 configuration compositional? The twelve genes on",
                 "four rulers inside Luminal and Basal, and the per-gene",
                 "between-subtype test E10 ran for ratios and never for genes."),
    not = paste("no re-scoring within stratum, no MYC stratification, no",
                "interaction, no Johnson-Neyman, no PAM50 classifier for CCLE,",
                "no re-fitting of B4 or E18."),
    hold_rule = HOLD_RULE,
    estimator = paste("partial Spearman on", PROLIF_COV, "- Fisher-z",
                      "intervals, Bonett-Wright variance. E16/E10 .cor_block",
                      "verbatim."),
    scoring = paste("POOLED scores, samples subset. GSVA and every z-score are",
                    "cohort-relative, so re-scoring inside a stratum would",
                    "build a different axis. Section 5 proves this did not",
                    "happen."),
    power = paste("Basal is 171 TCGA / 317 SCAN-B. A 171-sample stratum gives",
                  "a 95% interval about +/- 0.15 wide, so a Basal null is WEAK",
                  "evidence and the informative comparison is between strata,",
                  "not against zero."),
    mouse_source = paste("myc_mouse @ e348dd8,",
                         "docs/2026-09-02_myc_oxphos_priming_gate_model.md,",
                         "sections 3.3 and 7. TRANSCRIBED, read-only.",
                         "gene_gate's bX_wt is in a GITIGNORED results object",
                         "and is unavailable at a pinned ref for Bbc3."),
    no_cross_species_values = paste("signs and orderings only. 24 animals",
                                    "against 1,095 tumours is not a comparison",
                                    "of magnitudes; a species is a cohort."),
    genes = PRIMING_ALL, rulers = RULERS, strata = STRATA_PRIMING,
    n3 = "these are transcripts; the word 'primed' is used of nothing",
    seed = PROJECT_SEED),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E19)
utils::write.csv(as.data.frame(config), PATH_E19_CFG, row.names = FALSE)
utils::write.csv(as.data.frame(between_gene), PATH_E19_BTW, row.names = FALSE)
message("   ", PATH_E19)
message("   ", PATH_E19_CFG)
message("   ", PATH_E19_BTW)
message("\nE19 done.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E19)

  # The branch, and the rule behind it.
  x$focus_branch; cat(x$focus_reading, "\n")
  x$spec$hold_rule
  x$holds_wide %>% as.data.frame()

  # BCL2L1 and the two reference genes, all three strata, every ruler.
  x$configuration %>%
    dplyr::filter(gene %in% c("BCL2L1", "BBC3", "MCL1")) %>%
    dplyr::transmute(gene, cohort, stratum, ruler, n, rho = round(rho, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                     ci_excludes_0) %>%
    dplyr::arrange(gene, cohort, stratum, ruler) %>%
    as.data.frame()

  # The between-subtype test. A pooled value outside BOTH strata is the D3/S1
  # signature - a between-subtype difference read as a within-subtype one.
  x$flagged %>% as.data.frame()
  x$between_gene %>% dplyr::filter(pooled_outside) %>%
    dplyr::transmute(cohort, ruler, gene, all = round(all, 3),
                     Luminal = round(Luminal, 3), Basal = round(Basal, 3)) %>%
    as.data.frame()

  # The gap by stratum, against E16's pooled and E18's CCLE values.
  x$gap_strata %>% dplyr::mutate(gap = round(gap, 3)) %>% as.data.frame()
  readRDS(file.path(DIR_RESULTS, "respiratory_rulers.rds"))$the_gap
  readRDS(file.path(DIR_RESULTS, "ccle_configuration.rds"))$ccle_gap

  # The reproduction check. Not a result.
  summary(x$repro$d_rho)

  # CCLE composition - or the recorded reason there is none.
  x$ccle$counts %>% as.data.frame(); x$ccle_skip

  # The mouse. Signs and orderings only.
  x$mouse_gene %>% as.data.frame()
  cat(x$mouse_reading, "\n")
}
