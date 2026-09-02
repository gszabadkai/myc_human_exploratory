# E12_borrowing_explainer.R
# =============================================================================
# AN EXPLAINER, NOT AN ANALYSIS. It computes nothing E11 has not already
# established; it exists to make one sentence legible to someone who has not
# read the methods:
#
#   "MYC borrows its ordering of the apoptotic machinery from OXPHOS."
#
# It is a separate script because it is pedagogy, it runs in seconds, and it
# must not force a re-run of E11 every time a panel is reworded.
#
# =============================================================================
# THE SENTENCE IS TRUE ABOUT THE ORDERING AND FALSE ABOUT THE MAGNITUDE
# =============================================================================
# Said carelessly, "MYC borrows from OXPHOS" implies MYC has no relationship of
# its own with these genes. THAT IS NOT WHAT THE DATA SAY, and the explainer has
# to be built on the true version or it will be shot down by the first person
# who checks it.
#
# Decompose MYC's per-gene correlation into the part predictable from OXPHOS and
# the part left over:
#
#   BORROWED   rho(MYC, OXPHOS) x rho(OXPHOS, gene). Ordered by localisation by
#              construction, because it is a scaled copy of the OXPHOS column.
#              SD about 0.10 in TCGA.
#   MYC'S OWN  what remains after conditioning on OXPHOS. SD about 0.15 - LARGER
#              than the borrowed part, far above the ~0.03 sampling noise at
#              n = 1,095, and therefore real gene-to-gene structure. But its
#              correlation with mitochondrial localisation is -0.04. It is real
#              and it is unsorted.
#
# So MYC's observed split of 0.19 is a genuinely ordered component diluted by a
# larger unordered one, and conditioning on OXPHOS removes the ordered part and
# leaves the rest. The honest sentence is:
#
#   THE ONLY PART OF MYC'S ASSOCIATION WITH THESE GENES THAT IS ORGANISED BY
#   MITOCHONDRIAL LOCALISATION IS THE PART IT SHARES WITH OXPHOS. MYC HAS ITS
#   OWN, LARGER ASSOCIATION - IT IS SIMPLY NOT ABOUT THE MITOCHONDRION.
#
# That last clause is a loose end worth naming rather than hiding: something
# orders MYC's per-gene associations and this study has not asked what.
#
# SCALE: linear DESeq2-normalised, rank-based throughout. SPECIES: human.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE12: making the borrowing argument legible\n", strrep("=", 78))
PATH_E12 <- file.path(DIR_RESULTS, "borrowing_explainer.rds")

# =============================================================================
# 1. The three columns of numbers the explainer is built from
# =============================================================================
message("\n1. inputs")

e08  <- readRDS(file.path(DIR_RESULTS, "strata_and_death_genes.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)

CANON <- sort(e08$canonical$gene)
MC    <- sd_$strip_refs$MITOCARTA_ALL
ID_T  <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)

tcga_lin <- readRDS(PATH_TCGA_LINEAR); LT <- tcga_lin$mat[, ID_T, drop = FALSE]
scanb_lin <- readRDS(PATH_SCANB_LINEAR); LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

# Partial Spearman on 44 rows: rank, project out ranked covariates, correlate.
.pg <- function(RG, ry, Z) {
  H <- qr(cbind(1, Z))
  ry <- ry - qr.fitted(H, ry)
  RG <- RG - t(qr.fitted(H, t(RG)))
  as.numeric(stats::cor(t(RG), ry))
}

COH <- list(
  TCGA = list(L = LT, res = .symbol_resolver(rownames(LT), NULL),
              myc = as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T]),
              ox  = as.numeric(mito$gsva_arms["OXPHOS subunits", ID_T]),
              pf  = as.numeric(mito$gsva_cov["PROLIF_DISJOINT", ID_T])),
  `SCAN-B` = list(L = LS, res = .symbol_resolver(rownames(LS), sc$symbol_map),
              myc = as.numeric(sc$gsva_new[MYC_REF, ID_S]),
              ox  = as.numeric(sc$gsva_arms["OXPHOS subunits", ID_S]),
              pf  = as.numeric(sc$gsva_cov["PROLIF_DISJOINT", ID_S])))

axes <- list(); genes <- list()
for (coh in names(COH)) {
  C <- COH[[coh]]
  G <- .gene_rows(CANON, C$L, C$res)$mat
  RG <- t(apply(G, 1L, rank))
  rp <- rank(C$pf); rm_ <- rank(C$myc); ro <- rank(C$ox)
  H  <- qr(cbind(1, rp))
  em <- rm_ - qr.fitted(H, rm_)
  eo <- ro  - qr.fitted(H, ro)
  r_MO <- stats::cor(em, eo)
  r_OG  <- .pg(RG, ro,  cbind(rp))
  r_MG  <- .pg(RG, rm_, cbind(rp))
  r_MGo <- .pg(RG, rm_, cbind(rp, ro))
  axes[[coh]] <- tibble::tibble(cohort = coh, sample = colnames(G),
                                myc = as.numeric(em) / stats::sd(em),
                                ox  = as.numeric(eo) / stats::sd(eo))
  genes[[coh]] <- tibble::tibble(
    cohort = coh, gene = rownames(G),
    mitochondrial = rownames(G) %in% MC,
    r_OXPHOS = r_OG, r_MYC = r_MG, r_MYC_given_OXPHOS = r_MGo,
    borrowed = r_MO * r_OG, r_axes = r_MO)
}
COHORT_ORDER <- c("TCGA", "SCAN-B")
axes  <- dplyr::bind_rows(axes) %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_ORDER))
genes <- dplyr::bind_rows(genes) %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_ORDER))

message("\n   the decomposition, per cohort:")
decomp <- genes %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(
    `rho(MYC, OXPHOS)`    = dplyr::first(r_axes),
    `SD of OXPHOS column` = stats::sd(r_OXPHOS),
    `SD of MYC column`    = stats::sd(r_MYC),
    `SD of borrowed part` = stats::sd(borrowed),
    `SD of MYC's own`     = stats::sd(r_MYC_given_OXPHOS),
    .groups = "drop")
decomp %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   MYC's OWN part is LARGER than the borrowed part. The sentence is",
        " about the\n   ORDERING, not the size - which is why the next table is",
        " the one that matters.")

message("\n   median correlation of mitochondrial vs cytosolic genes:")
gaps <- genes %>%
  tidyr::pivot_longer(c(r_OXPHOS, r_MYC, r_MYC_given_OXPHOS, borrowed),
                      names_to = "column", values_to = "rho") %>%
  dplyr::group_by(cohort, column, mitochondrial) %>%
  dplyr::summarise(median = stats::median(rho), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = mitochondrial, values_from = median,
                     names_prefix = "mito_") %>%
  dplyr::mutate(gap = mito_TRUE - mito_FALSE,
                column = factor(column,
                  levels = c("r_OXPHOS", "borrowed", "r_MYC",
                             "r_MYC_given_OXPHOS"),
                  labels = c("OXPHOS", "what MYC borrows from OXPHOS",
                             "MYC, as observed", "MYC, after removing OXPHOS")))
gaps %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(cohort, column) %>% as.data.frame() %>% print(row.names = FALSE)
message("   The gap is the separation between mitochondrial and cytosolic",
        " genes. It is\n   OXPHOS's; MYC shows a diluted copy of it; and",
        " nothing is left once OXPHOS goes.")

# =============================================================================
# 2. The panel
# =============================================================================
# Three steps, in the order someone would ask them.
#   A  the two axes travel together, so anything tracking one appears to track
#      the other
#   B  the separation between mitochondrial and cytosolic genes belongs to
#      OXPHOS: big on OXPHOS, a diluted copy on MYC, gone once OXPHOS is removed
#   C  and yet MYC's leftover is not zero - it is real and it is unsorted
message("\n2. the panel")

MITO_COLS <- c(`TRUE` = "#d7191c", `FALSE` = "grey55")
theme_e12 <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey92",
                                                          colour = NA),
                 legend.position = "bottom",
                 plot.subtitle = ggplot2::element_text(size = 8.5),
                 plot.caption = ggplot2::element_text(size = 7.5,
                                                      colour = "grey35",
                                                      hjust = 0,
                                                      lineheight = 1.2))

ax_lab <- genes %>% dplyr::distinct(cohort, r_axes) %>%
  dplyr::mutate(lab = sprintf("rho = %.2f", r_axes))
pA <- ggplot2::ggplot(axes, ggplot2::aes(myc, ox)) +
  ggplot2::geom_point(size = 0.35, alpha = 0.18, colour = "#1f4e79") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                       colour = "black", linewidth = 0.5) +
  ggplot2::geom_text(data = ax_lab, ggplot2::aes(label = lab), x = -2.6, y = 2.6,
                     hjust = 0, size = 3, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::labs(
    subtitle = paste("A.  The two axes travel together. A tumour high in MYC",
                     "tends to be high in OXPHOS."),
    x = "MYC activity (ranked, proliferation removed)",
    y = "OXPHOS (ranked)") +
  theme_e12

pB <- gaps %>%
  dplyr::mutate(column = factor(column, levels = rev(levels(gaps$column)))) %>%
  ggplot2::ggplot(ggplot2::aes(y = column)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_segment(ggplot2::aes(x = mito_FALSE, xend = mito_TRUE,
                                     yend = column), linewidth = 0.6,
                        colour = "grey55") +
  ggplot2::geom_point(ggplot2::aes(x = mito_FALSE), colour = MITO_COLS[["FALSE"]],
                      size = 3) +
  ggplot2::geom_point(ggplot2::aes(x = mito_TRUE), colour = MITO_COLS[["TRUE"]],
                      size = 3) +
  ggplot2::geom_text(ggplot2::aes(x = (mito_FALSE + mito_TRUE) / 2,
                                  label = sprintf("%.2f", gap)),
                     vjust = -1.1, size = 2.7, colour = "grey25") +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::labs(
    subtitle = paste("B.  The separation between mitochondrial (red) and",
                     "cytosolic (grey) genes belongs to OXPHOS."),
    x = "median correlation of the 44 genes with that column", y = NULL) +
  theme_e12

pC <- genes %>%
  tidyr::pivot_longer(c(r_OXPHOS, r_MYC_given_OXPHOS), names_to = "column",
                      values_to = "rho") %>%
  dplyr::mutate(column = factor(column,
    levels = c("r_OXPHOS", "r_MYC_given_OXPHOS"),
    labels = c("OXPHOS", "MYC, after removing OXPHOS"))) %>%
  ggplot2::ggplot(ggplot2::aes(column, rho, colour = mitochondrial)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.13,
                                                          height = 0, seed = 1),
                      size = 1.5, alpha = 0.85) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_colour_manual(values = MITO_COLS,
                               labels = c(`FALSE` = "cytosolic",
                                          `TRUE` = "mitochondrial"),
                               name = NULL) +
  ggplot2::labs(
    subtitle = paste("C.  What is left of MYC is still spread out - it is real",
                     "- but the colours are now mixed."),
    x = NULL, y = "correlation of each gene with the column") +
  theme_e12

pE <- patchwork::wrap_plots(pA, pB, pC, ncol = 1, heights = c(1, 0.85, 1)) +
  patchwork::plot_annotation(
    title = "Why MYC looks like it orders the death machinery, and why it does not",
    caption = paste0(
      "Each of the 44 genes has a correlation with OXPHOS and a correlation with MYC. Because the two axes are themselves\n",
      "correlated (A), a gene that tracks OXPHOS is bound to track MYC a little as well - it inherits a scaled copy. B shows this\n",
      "happening: OXPHOS separates mitochondrial from cytosolic genes by about 0.41, the copy MYC inherits carries the same\n",
      "separation shrunk to about 0.13, MYC as observed shows about 0.11 - no more than it inherited - and after OXPHOS is\n",
      "taken out, nothing is left.\n",
      "C is the honest qualifier. MYC's leftover correlations are NOT zero - they spread as widely as before, far beyond the\n",
      "+/- 0.03 that sampling alone would give at these cohort sizes. What has gone is the SORTING: red and grey now interleave.\n",
      "So MYC does relate to these genes. That relationship is simply not about the mitochondrion, and this study has not asked\n",
      "what it is about. EXPLORATORY - not pre-registered."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.caption = ggplot2::element_text(size = 7.5, colour = "grey35",
                                           hjust = 0, lineheight = 1.2)))
for (ext in c("png", "pdf")) {
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("E12_borrowing_explainer.", ext)),
                  pE, width = 8.5, height = 11, dpi = 300)
}
message("   E12_borrowing_explainer")

saveRDS(list(genes = genes, gaps = gaps, decomp = decomp,
             rules = list(
               sentence = paste("'MYC borrows its ordering from OXPHOS' is true",
                                "about the ORDERING and false about the",
                                "MAGNITUDE. MYC's own component is larger than",
                                "the borrowed one; it is simply unsorted with",
                                "respect to mitochondrial localisation."),
               loose_end = paste("MYC's residual per-gene correlations have SD",
                                 "~0.15 against ~0.03 sampling noise, so",
                                 "something orders them. This study has not",
                                 "asked what.")),
             built = Sys.time()), PATH_E12)
message("\nE12: done.  results/borrowing_explainer.rds + 1 figure")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E12)
  x$decomp %>% as.data.frame()
  x$gaps %>% dplyr::arrange(cohort, column) %>% as.data.frame()

  # the genes MYC still relates to once OXPHOS is gone, which is the loose end
  x$genes %>%
    dplyr::select(cohort, gene, mitochondrial, r_MYC_given_OXPHOS) %>%
    tidyr::pivot_wider(names_from = cohort, values_from = r_MYC_given_OXPHOS) %>%
    dplyr::filter(sign(TCGA) == sign(`SCAN-B`)) %>%
    dplyr::arrange(dplyr::desc(abs(TCGA) + abs(`SCAN-B`))) %>%
    utils::head(15) %>% as.data.frame()

}
