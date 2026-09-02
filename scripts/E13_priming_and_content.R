# E13_priming_and_content.R
# =============================================================================
# RECONCILING THE PRIMING RATIOS WITH THE LOCALISATION RESULT.
#
# E11 concluded that the apoptotic machinery is ordered along OXPHOS by WHERE
# each protein sits, and that the ordering is no steeper than a gene set matched
# for expression and sub-mitochondrial compartment. A fair reading of that is
# "the OXPHOS score is a mitochondrial content index, and the death genes are
# along for the ride".
#
# THE PRIMING RATIOS CAN TEST THAT DIRECTLY, AND THE MACHINERY ANALYSIS COULD
# NOT. Eleven of the twelve priming genes are outer-mitochondrial-membrane
# proteins - only BMF is absent from MitoCarta. A ratio of two MOM transcripts
# CANCELS anything that scales with how many mitochondria a tumour has. So:
#
#   IF the OXPHOS association is mitochondrial abundance, every MOM/MOM ratio
#   is flat. IF it is not, they are not.
#
# That is a control on a different axis from the composition null: the null asks
# whether the machinery is special AMONG gene sets, this asks whether the SIGNAL
# is abundance at all.
#
# =============================================================================
# AND THE SAME DISCIPLINE APPLIES TO THE ANSWER
# =============================================================================
# If the BCL2-family MOM/MOM ratios are not flat, the next question is whether
# RANDOM MOM/MOM pairs are flat either. If they are not, the within-compartment
# resolution is a property of the OXPHOS axis and not of the BCL2 family, and
# the priming ratios are once again real, OXPHOS-carried, and unexceptional.
# Section 3 draws that null before the answer is looked at.
#
# SCALE: log2(linear + 1) for the ratios, which are differences on that scale.
# All correlations rank-based and adjusted for proliferation. SPECIES: human.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE13: are the priming ratios mitochondrial content?\n", strrep("=", 78))
PATH_E13   <- file.path(DIR_RESULTS, "priming_and_content.rds")
NULL_PAIRS <- 500L

# =============================================================================
# 1. Inputs, and the compartment of every priming gene
# =============================================================================
message("\n1. inputs")

sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)
mcs  <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
SUB  <- stats::setNames(mcs[["MitoCarta3.0_SubMitoLocalization"]], mcs$Symbol)
MOM_GENES <- names(SUB)[!is.na(SUB) & SUB == "MOM"]

PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
PRIMING_ALL  <- union(PRIMING_PRO, PRIMING_ANTI)

compartment <- tibble::tibble(
  gene = PRIMING_ALL,
  side = ifelse(PRIMING_ALL %in% PRIMING_PRO, "pro-apoptotic", "anti-apoptotic"),
  submito = ifelse(is.na(SUB[PRIMING_ALL]), "not in MitoCarta",
                   SUB[PRIMING_ALL]))
message("\n   where the 12 priming proteins sit:")
compartment %>% dplyr::count(submito, name = "genes") %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   THE LOCALISATION RULE CANNOT ORDER THESE - they share a",
        " compartment.\n   Whatever separates them on the OXPHOS axis is",
        " something else.")

RATIO_GRID <- tidyr::expand_grid(pro = PRIMING_PRO, anti = PRIMING_ANTI) %>%
  dplyr::mutate(ratio = paste0(pro, "/", anti),
                pro_mito = !is.na(SUB[pro]), anti_mito = !is.na(SUB[anti]),
                pair = dplyr::if_else(pro_mito & anti_mito,
                                      "both mitochondrial", "one of each"))
message("\n   ", sum(RATIO_GRID$pair == "both mitochondrial"), " of ",
        nrow(RATIO_GRID), " ratios cancel mitochondrial content by construction")

.pg <- function(RG, ry, Z) {
  H <- qr(cbind(1, Z))
  ry <- ry - qr.fitted(H, ry)
  RG <- RG - t(qr.fitted(H, t(RG)))
  as.numeric(stats::cor(t(RG), ry))
}

COH <- c("TCGA", "SCAN-B")
per_cohort <- list()
set.seed(PROJECT_SEED)
for (coh in COH) {
  lin <- readRDS(if (coh == "TCGA") PATH_TCGA_LINEAR else PATH_SCANB_LINEAR)
  ids <- if (coh == "TCGA") colnames(mito$gsva_arms) else colnames(sc$gsva_arms)
  L   <- lin$mat[, ids, drop = FALSE]; rm(lin); invisible(gc(verbose = FALSE))
  res <- if (coh == "TCGA") .symbol_resolver(rownames(L), NULL) else
         .symbol_resolver(rownames(L), sc$symbol_map)
  ox  <- if (coh == "TCGA") as.numeric(mito$gsva_arms["OXPHOS subunits", ids]) else
         as.numeric(sc$gsva_arms["OXPHOS subunits", ids])
  myc <- if (coh == "TCGA") as.numeric(nw$tcga_gsva_new[MYC_REF, ids]) else
         as.numeric(sc$gsva_new[MYC_REF, ids])
  pf  <- if (coh == "TCGA") as.numeric(mito$gsva_cov["PROLIF_DISJOINT", ids]) else
         as.numeric(sc$gsva_cov["PROLIF_DISJOINT", ids])

  G  <- log2(.gene_rows(PRIMING_ALL, L, res)$mat + 1)
  own <- match(rownames(G), rownames(L))
  RG <- t(apply(G, 1L, rank))
  rp <- rank(pf); rm_ <- rank(myc); ro <- rank(ox)

  # the 12 genes on their own
  gene_rho <- tibble::tibble(
    cohort = coh, gene = rownames(G),
    ox_rho  = .pg(RG, ro,  cbind(rp)),
    myc_rho = .pg(RG, rm_, cbind(rp)))

  # the 35 ratios, before and after conditioning on the other axis
  R  <- G[RATIO_GRID$pro, , drop = FALSE] - G[RATIO_GRID$anti, , drop = FALSE]
  rownames(R) <- RATIO_GRID$ratio
  RR <- t(apply(R, 1L, rank))
  ratio_rho <- RATIO_GRID %>% dplyr::mutate(
    cohort    = coh,
    ox_rho    = .pg(RR, ro,  cbind(rp)),
    ox_condM  = .pg(RR, ro,  cbind(rp, rm_)),
    myc_rho   = .pg(RR, rm_, cbind(rp)),
    myc_condO = .pg(RR, rm_, cbind(rp, ro)))

  # THE NULL: random MOM/MOM pairs, expression-matched, same 30-pair structure.
  bm  <- RATIO_GRID %>% dplyr::filter(pair == "both mitochondrial")
  gu  <- unique(c(bm$pro, bm$anti))
  own_bm <- match(gu, rownames(L))
  keep <- which(apply(L, 1L, stats::sd) > 0)
  pool <- union(intersect(which(rownames(L) %in% MOM_GENES), keep), own_bm)
  B <- .expression_bins(L, pool, n_bins = 6L)
  obs_bm <- mean(abs(ratio_rho$ox_rho[ratio_rho$pair == "both mitochondrial"]))
  nd <- vapply(seq_len(NULL_PAIRS), function(b) {
    d  <- .matched_draw(own_bm, B, exclude = own_bm)
    Gd <- log2(L[d, , drop = FALSE] + 1); rownames(Gd) <- gu
    Rd <- Gd[bm$pro, , drop = FALSE] - Gd[bm$anti, , drop = FALSE]
    mean(abs(.pg(t(apply(Rd, 1L, rank)), ro, cbind(rp))))
  }, numeric(1))

  per_cohort[[coh]] <- list(
    gene_rho = gene_rho, ratio_rho = ratio_rho,
    null = tibble::tibble(cohort = coh, observed = obs_bm,
                          null_mean = mean(nd), null_sd = stats::sd(nd),
                          z = (obs_bm - mean(nd)) / stats::sd(nd),
                          pct_below = mean(nd < obs_bm)),
    ox = ox, pf = pf, G = G[c("BAD", "MCL1"), , drop = FALSE])
  message("   ", coh, " done")
}
gene_rho  <- dplyr::bind_rows(lapply(per_cohort, `[[`, "gene_rho"))
ratio_rho <- dplyr::bind_rows(lapply(per_cohort, `[[`, "ratio_rho"))
mom_null  <- dplyr::bind_rows(lapply(per_cohort, `[[`, "null"))

# =============================================================================
# 2. Is the OXPHOS signal mitochondrial abundance?
# =============================================================================
message("\n2. would abundance explain it?")

message("\n   the 12 priming genes on the OXPHOS axis, by compartment:")
gene_rho %>% dplyr::left_join(compartment, by = "gene") %>%
  dplyr::select(cohort, gene, side, submito, ox_rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = ox_rho) %>%
  dplyr::arrange(dplyr::desc(TCGA)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   ELEVEN OF TWELVE ARE OUTER-MEMBRANE AND THEY SPAN THE WHOLE RANGE.",
        "\n   Mitochondrial abundance would push same-compartment transcripts",
        " the SAME way.")

message("\n   mean |rho| of the ratios, by whether content cancels:")
ratio_rho %>% dplyr::group_by(cohort, pair) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean_abs_OXPHOS = mean(abs(ox_rho)),
                   max_abs_OXPHOS  = max(abs(ox_rho)),
                   mean_abs_MYC    = mean(abs(myc_rho)), .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   `both mitochondrial` ratios would be FLAT if the signal were",
        " abundance.\n   They are not - they are the stronger group.")

message("\n   and do the ratios show the same asymmetry as the 44 genes?")
ratio_rho %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(OXPHOS = mean(abs(ox_rho)),
                   `OXPHOS, MYC removed` = mean(abs(ox_condM)),
                   MYC = mean(abs(myc_rho)),
                   `MYC, OXPHOS removed` = mean(abs(myc_condO)),
                   .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 3. But is it the BCL2 family, or the axis?
# =============================================================================
message("\n3. against random same-compartment pairs")
mom_null %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   If this does not replicate, the within-compartment resolution",
        " belongs to the\n   OXPHOS axis and not to the BCL2 family - which is",
        " the same verdict E11\n   reached for the localisation split.")

# =============================================================================
# 4. The figure
# =============================================================================
message("\n4. figure")

COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
theme_e13 <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey92",
                                                          colour = NA),
                 legend.position = "bottom",
                 plot.subtitle = ggplot2::element_text(size = 8.5),
                 plot.caption = ggplot2::element_text(size = 7.5,
                                                      colour = "grey35",
                                                      hjust = 0,
                                                      lineheight = 1.2))

gA <- gene_rho %>% dplyr::left_join(compartment, by = "gene") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                gene = stats::reorder(gene, ox_rho)) %>%
  ggplot2::ggplot(ggplot2::aes(ox_rho, gene, colour = cohort,
                               shape = submito)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(MOM = 16, `not in MitoCarta` = 1),
                              name = NULL) +
  ggplot2::labs(
    subtitle = paste("A.  Eleven of the twelve priming proteins sit in the",
                     "SAME membrane, yet they span the whole range."),
    x = "correlation with OXPHOS (proliferation removed)", y = NULL) +
  theme_e13

gB_dat <- dplyr::bind_rows(lapply(names(per_cohort), function(coh) {
  p <- per_cohort[[coh]]
  tibble::tibble(cohort = coh, ox = p$ox,
                 BAD = p$G["BAD", ], MCL1 = p$G["MCL1", ]) %>%
    tidyr::pivot_longer(c(BAD, MCL1), names_to = "gene", values_to = "expr")
})) %>% dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
gB_lab <- gene_rho %>% dplyr::filter(gene %in% c("BAD", "MCL1")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                lab = sprintf("rho = %+.2f", ox_rho))
gB <- ggplot2::ggplot(gB_dat, ggplot2::aes(ox, expr)) +
  ggplot2::geom_point(size = 0.3, alpha = 0.15, colour = "grey35") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                       colour = "#d7191c", linewidth = 0.6) +
  ggplot2::geom_text(data = gB_lab, ggplot2::aes(label = lab), x = -Inf,
                     y = Inf, hjust = -0.15, vjust = 1.5, size = 2.8,
                     inherit.aes = FALSE) +
  ggplot2::facet_grid(cohort ~ gene, scales = "free_y") +
  ggplot2::labs(
    subtitle = paste("B.  BAD and MCL1 are both outer-membrane BCL2-family",
                     "proteins - and they tilt opposite ways."),
    x = "OXPHOS score", y = "log2 expression") +
  theme_e13

gC_dat <- ratio_rho %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
gC_null <- mom_null %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                pair = "both mitochondrial",
                lo = null_mean - 1.96 * null_sd,
                hi = null_mean + 1.96 * null_sd)
gC <- ggplot2::ggplot(gC_dat, ggplot2::aes(abs(ox_rho), pair)) +
  ggplot2::geom_linerange(data = gC_null,
                          ggplot2::aes(xmin = lo, xmax = hi, y = pair),
                          inherit.aes = FALSE, linewidth = 6,
                          colour = "grey86") +
  ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.13,
                                                          width = 0, seed = 1),
                      size = 1.6, alpha = 0.8, colour = "#1b9e77") +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::labs(
    subtitle = paste("C.  But random same-compartment pairs do the same thing",
                     "(grey), and the excess does not replicate."),
    x = "|correlation of the ratio with OXPHOS|", y = NULL) +
  theme_e13

gE <- patchwork::wrap_plots(gA, gB, gC, ncol = 1, heights = c(1, 1.1, 0.7)) +
  patchwork::plot_annotation(
    title = "The priming ratios say the OXPHOS axis is not just counting mitochondria",
    caption = paste0(
      "A ratio of two transcripts from the same compartment cancels anything that scales with how many mitochondria a tumour\n",
      "has. Eleven of the twelve priming genes are outer-membrane proteins, so most of these ratios are that control - and they\n",
      "are not flat. B makes it concrete: BAD and MCL1 are both outer-membrane BCL2-family proteins, and across the OXPHOS axis\n",
      "one rises while the other falls. Abundance cannot produce that; it would move them together.\n",
      "C is the limit. Random outer-membrane pairs, expression-matched, behave the same way, and the BCL2-family excess over\n",
      "them is z = 1.2 in TCGA and z = -0.1 in SCAN-B - it does not replicate. So the OXPHOS axis resolves BETWEEN transcripts\n",
      "that share a compartment, which is a real property of the axis and not an artefact of mitochondrial mass; but that\n",
      "resolution is not a property of the apoptotic machinery in particular. EXPLORATORY - not pre-registered."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.caption = ggplot2::element_text(size = 7.5, colour = "grey35",
                                           hjust = 0, lineheight = 1.2)))
for (ext in c("png", "pdf")) {
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("E13_priming_and_content.", ext)),
                  gE, width = 9, height = 11, dpi = 300)
}
message("   E13_priming_and_content")

saveRDS(list(compartment = compartment, gene_rho = gene_rho,
             ratio_rho = ratio_rho, mom_null = mom_null,
             settings = list(null_pairs = NULL_PAIRS, seed = PROJECT_SEED),
             rules = list(
               control = paste("a MOM/MOM ratio cancels mitochondrial",
                               "abundance by construction, which is a control",
                               "the 44-gene analysis could not run because",
                               "those genes are spread across compartments"),
               verdict = paste("the ratios are NOT flat, so the OXPHOS axis is",
                               "not a mitochondrial mass index; but random",
                               "MOM/MOM pairs are not flat either and the",
                               "BCL2-family excess does not replicate, so the",
                               "within-compartment resolution belongs to the",
                               "axis and not to the death machinery")),
             built = Sys.time()), PATH_E13)
message("\nE13: done.  results/priming_and_content.rds + 1 figure")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E13)
  x$compartment %>% as.data.frame()
  x$mom_null %>% as.data.frame()

  # the ratios that cancel content and still correlate
  x$ratio_rho %>% dplyr::filter(pair == "both mitochondrial") %>%
    dplyr::select(cohort, ratio, ox_rho, ox_condM) %>%
    tidyr::pivot_wider(names_from = cohort, values_from = c(ox_rho, ox_condM)) %>%
    dplyr::arrange(dplyr::desc(abs(ox_rho_TCGA))) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

}
