# E46_fig4_quadrant_heatmap.R
# =============================================================================
# FIG. 4 HUMAN PANELS, PART 2: THE CONFIGURATION BY QUADRANT, ON BOTH
# INSTRUMENTS. DISPLAY ONLY.
#
# EXPLORATORY. Nothing here is pre-registered, and NOTHING HERE IS ESTIMATED:
# every plotted number is an adjusted group mean that E34 already fitted and
# saved. No declaration, on the E15 precedent for display scripts.
#
# Specification: docs/2026-10-22_fig4_human_panels_proposal.md section 2, as
# amended by the author's decisions in its section 8. Authority on what the
# panel may assert: E34's result note section 9 (docs/2026-10-01_e34_result.md)
# and docs/2026-10-04_handoff.md section 2. Nothing here relaxes either.
#
# SNAPSHOT: every value derives from data/from_validation/ through E34. The
# joint snapshot (results/joint_tcga_*.rds) is NOT read.
#
# =============================================================================
# THE AUTHOR'S DECISIONS (2026-10-22, proposal section 8)
# =============================================================================
#   1. ROWS ARE BLOCKED BY FUNCTION: pro-apoptotic, then anti-apoptotic. v5's
#      MitoCarta / cytosolic split cannot be drawn - 11 of the 12 genes are in
#      Human MitoCarta 3.0; only BMF is not. The order is FIXED, not
#      data-driven: composite members first in each block.
#   2. BOTH INSTRUMENTS GO IN THE MAIN FIGURE, side by side: the OXPHOS split
#      by GSVA and by mitoPPS. MYC activity (M_a) is the same in both; only the
#      OXPHOS call differs. The two assign the same quadrant to 85.0% of TCGA
#      and 82.7% of SCAN-B tumours (E34 8.1), so they are NOT the same
#      tumours in every column, and the legend says so.
#
# This replaces the proposal's "GSVA main, mitoPPS in ED Fig. 7". The reason
# the proposal gave for GSVA alone (one instrument through 4g-4j) is outweighed
# by the arm's two-instrument convention, and E34 2.2's observation that GSVA
# shows the SMALLER MYC contrast in every comparable cell is now visible on the
# page instead of only in the text.
#
# =============================================================================
# WHAT THE PANEL MAY AND MAY NOT ASSERT (E34 result note section 9)
# =============================================================================
#   - THE OXPHOS SEPARATION IS NEAR-GUARANTEED BY CONSTRUCTION: the composite's
#     genes and signs were read off OXPHOS correlations in these same two
#     cohorts (E34 1.2). The informative contrasts are the MYC ones, which is
#     why the MYC pairs sit side by side within each respiratory level.
#   - The MYC contrast is SMALL BUT NOT ZERO (6 of 8 pooled intervals exclude
#     zero). Nothing here may be read as "no MYC difference".
#   - NO PRODUCT TERM WAS FITTED. No line joins any two quadrants, in the strip
#     or anywhere else; a line would draw an interaction the analysis does not
#     have. Asserted on the built plot.
#   - Transcript abundance, not a functional readout; nothing is "primed".
#   - No treatment-selection, causal or prognostic claim.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing. Cells are E34's $groups, subtype "all", m1; the strip
#     adds m3. The only computation on per-sample data is section 3's CHECK of
#     the composite's signs, and nothing it computes is plotted.
#   - It does not draw E34's m2 (within subtype) or the arms; those are E34's
#     own figure and ED material.
#   - It does not interpret. The legend draft is DESCRIPTIVE ONLY.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})

message("\nE46: Fig. 4 human panels, part 2 - the configuration by quadrant\n",
        strrep("=", 78))
message("DISPLAY ONLY. Every value is an adjusted mean E34 already saved.")

PATH_E34 <- file.path(DIR_RESULTS, "e34_quadrant_configuration.rds")
PATH_E46 <- file.path(DIR_RESULTS, "e46_fig4_quadrant_panels.rds")

DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_STEM_DOC    <- "2026-10-22_E46_quadrant_configuration"
PATH_FIG_BOTH   <- file.path(DIR_FIGURES, "E46_quadrant_configuration")
PATH_FIG_GSVA   <- file.path(DIR_FIGURES, "E46_quadrant_configuration_gsva")
PATH_FIG_MPPS   <- file.path(DIR_FIGURES, "E46_quadrant_configuration_mitopps")
PATH_FIG_DOC    <- file.path(DIR_DOC_FIGURES, FIG_STEM_DOC)

# =============================================================================
# 0. CONSTANTS - fixed by the specification, none chosen on the data
# =============================================================================
COHORTS     <- c("TCGA", "SCAN-B")
INSTRUMENTS <- c(gsva = "GSVA", mitopps = "mitoPPS")
QLEV <- c(Q1 = "Q1 MYC-low OXPHOS-low",  Q2 = "Q2 MYC-low OXPHOS-high",
          Q3 = "Q3 MYC-high OXPHOS-low", Q4 = "Q4 MYC-high OXPHOS-high")

# Columns: within each cohort, OXPHOS-low (Q1, Q3) then OXPHOS-high (Q2, Q4),
# so each MYC pair sits side by side at fixed respiratory level.
X_COL <- list(TCGA     = c(Q1 = 1.00, Q3 = 2.00, Q2 = 3.25, Q4 = 4.25),
              `SCAN-B` = c(Q1 = 5.85, Q3 = 6.85, Q2 = 8.10, Q4 = 9.10))
X_LIM <- c(0.45, 9.65)

# Rows, blocked by function. Composite members first in each block, then the
# rest in E10's order. FIXED: nothing here is ordered on the data.
PRO  <- c("BAD", "BIK", "BBC3", "BID", "BCL2L11", "BMF", "PMAIP1")
ANTI <- c("BCL2L1", "MCL1", "BCL2", "BCL2L2", "BCL2A1")
Y_GENE <- c(stats::setNames(seq(12.6, 6.6, by = -1), PRO),
            stats::setNames(5:1, ANTI))
Y_OXLEV  <- 13.55
Y_COHORT <- 14.35
Y_LIM    <- c(0.45, 14.85)

# E33's composite, verbatim from E34 (scripts/E34_quadrant_configuration.R:141):
# mean z of BBC3, BID, BIK, BAD, BCL2L1 and minus MCL1. Section 3 CHECKS these
# signs against the saved per-sample composite rather than trusting the copy.
COMPOSITE_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)

# The diverging scale. White at 0, symmetric. ColorBrewer RdBu's end points,
# which ColorBrewer marks colour-blind safe; blue below the cohort mean, red
# above. Limits fixed at +/-0.6 (largest |gene mean| 0.50) and ASSERTED to
# contain every cell, so no value is squished silently.
FILL_LIM  <- c(-0.6, 0.6)
FILL_LOW  <- "#2166AC"
FILL_HIGH <- "#B2182B"
STRIP_LIM <- c(-0.8, 0.6)
DX_MODEL  <- c(m1 = -0.17, m3 = 0.17)

BASE_PT <- 6.5
TXT_MM  <- 6 / ggplot2::.pt
SUB_MM  <- 5.5 / ggplot2::.pt

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.pr      <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub     <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")
.fmt_n   <- function(x) format(x, big.mark = ",", trim = TRUE)

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
.dev <- function(a, b) if (length(a)) max(abs(a - b)) else NA_real_

# =============================================================================
# 1. Input - E34's saved object
# =============================================================================
message("\n1. input")

.stop_if(file.exists(PATH_E34), "absent: ", PATH_E34,
         ". Source scripts/E34_quadrant_configuration.R first.")
e34 <- readRDS(PATH_E34)
message("   read: E34 $groups, $contrasts, $quadrant, $analysis_frames")

G <- e34$groups %>%
  dplyr::filter(subtype == "all", model %in% c("m1", "m3")) %>%
  dplyr::mutate(q = names(QLEV)[match(quadrant, QLEV)])
.assert("every quadrant label is one of E34's four", !anyNA(G$q), n = nrow(G))

genes_in <- sub("^g_", "", unique(G$readout_col[grepl("^g_", G$readout_col)]))
.assert("E34 carries exactly the twelve genes the rows name",
        setequal(genes_in, c(PRO, ANTI)) && length(genes_in) == 12L,
        n = length(genes_in), detail = paste(sort(genes_in), collapse = ", "))
.assert("every composite member is one of the twelve",
        all(names(COMPOSITE_SIGN) %in% c(PRO, ANTI)), n = 6L)

# One row per cohort x instrument x model x quadrant x readout, no gaps.
cover <- G %>%
  dplyr::filter(readout_col %in% c("config_comp", paste0("g_", c(PRO, ANTI)))) %>%
  dplyr::count(cohort, instrument, model, q, readout_col)
.assert("every cell present exactly once: 2 cohorts x 2 instruments x 2 models x 4 quadrants x 13 readouts",
        nrow(cover) == 2L * 2L * 2L * 4L * 13L && all(cover$n == 1L),
        n = nrow(cover))

# Analysis-set sizes: one per cohort, identical across instruments and models.
nset <- G %>% dplyr::filter(readout_col == "config_comp") %>%
  dplyr::group_by(cohort, instrument, model) %>%
  dplyr::summarise(n = sum(n_cell), .groups = "drop")
N_COH <- nset %>% dplyr::distinct(cohort, n)
.assert("one analysis set per cohort, the same on both instruments and both models",
        nrow(N_COH) == 2L, n = nrow(nset),
        detail = paste(N_COH$cohort, N_COH$n, collapse = "; "))
N_COH <- stats::setNames(N_COH$n, N_COH$cohort)
.pr(nset)

# =============================================================================
# 2. The heatmap's cells (m1) and the strip's points (m1, m3)
# =============================================================================
message("\n2. the cells and the strip")

.xpos <- function(cohort, q) mapply(function(c, k) X_COL[[c]][[k]], cohort, q)

CELLS <- G %>%
  dplyr::filter(model == "m1", grepl("^g_", readout_col)) %>%
  dplyr::mutate(gene = sub("^g_", "", readout_col),
                block = ifelse(gene %in% PRO, "pro-apoptotic", "anti-apoptotic"),
                composite_sign = unname(COMPOSITE_SIGN[gene]),
                x = unname(.xpos(cohort, q)), y = unname(Y_GENE[gene])) %>%
  dplyr::select(instrument, cohort, q, quadrant, gene, block, composite_sign,
                n_cell, est, se, lo, hi, x, y) %>%
  dplyr::arrange(instrument, cohort, -y, x)
.assert("no heatmap cell lies outside the colour limits (nothing squished)",
        all(abs(CELLS$est) <= FILL_LIM[2]), n = nrow(CELLS),
        detail = sprintf("max |est| %.3f against %.2f", max(abs(CELLS$est)),
                         FILL_LIM[2]))

STRIP <- G %>%
  dplyr::filter(readout_col == "config_comp") %>%
  dplyr::mutate(x = unname(.xpos(cohort, q)) + unname(DX_MODEL[model])) %>%
  dplyr::select(instrument, cohort, model, q, quadrant, n_cell, est, se, lo,
                hi, x) %>%
  dplyr::arrange(instrument, cohort, x)
.assert("every strip interval lies inside the strip's axis",
        all(STRIP$lo >= STRIP_LIM[1]) && all(STRIP$hi <= STRIP_LIM[2]),
        n = nrow(STRIP), detail = sprintf("limits %.2f to %.2f, data %.3f to %.3f",
                                          STRIP_LIM[1], STRIP_LIM[2],
                                          min(STRIP$lo), max(STRIP$hi)))

# The proposal printed these GSVA values (section 2's table) from the same
# object. Reproducing them to 2 dp checks that the frame is reading the cells
# the note described, in the column order it described.
PROPOSAL_TAB <- tibble::tribble(
  ~readout,      ~model, ~cohort,  ~Q1,   ~Q3,   ~Q2,   ~Q4,
  "config_comp", "m1",   "TCGA",   -0.58, -0.44,  0.41,  0.40,
  "config_comp", "m1",   "SCAN-B", -0.36, -0.24,  0.23,  0.37,
  "config_comp", "m3",   "TCGA",   -0.65, -0.37,  0.35,  0.47,
  "config_comp", "m3",   "SCAN-B", -0.44, -0.19,  0.19,  0.45,
  "g_BAD",       "m1",   "TCGA",   -0.48, -0.50,  0.42,  0.29,
  "g_BAD",       "m1",   "SCAN-B", -0.33, -0.34,  0.38,  0.31,
  "g_BCL2L1",    "m1",   "TCGA",   -0.34, -0.33,  0.41,  0.26,
  "g_BCL2L1",    "m1",   "SCAN-B", -0.23, -0.18,  0.24,  0.21,
  "g_BID",       "m1",   "TCGA",   -0.44, -0.12,  0.16,  0.35,
  "g_BID",       "m1",   "SCAN-B", -0.16,  0.02, -0.05,  0.15,
  "g_PMAIP1",    "m1",   "TCGA",    0.18,  0.04,  0.08, -0.19,
  "g_PMAIP1",    "m1",   "SCAN-B",  0.09,  0.00,  0.09, -0.11,
  "g_BCL2L11",   "m1",   "TCGA",    0.37,  0.20, -0.39, -0.16,
  "g_BCL2L11",   "m1",   "SCAN-B", -0.02,  0.12, -0.02, -0.04) %>%
  tidyr::pivot_longer(c(Q1, Q3, Q2, Q4), names_to = "q", values_to = "printed")
pt <- PROPOSAL_TAB %>%
  dplyr::left_join(G %>% dplyr::filter(instrument == "gsva") %>%
                     dplyr::select(readout = readout_col, model, cohort, q, est),
                   by = c("readout", "model", "cohort", "q"))
.assert("the proposal's section-2 table is reproduced to 2 dp",
        !anyNA(pt$est) && all(abs(pt$est - pt$printed) <= 0.005 + 1e-9),
        n = nrow(pt), max_dev = .dev(pt$est, pt$printed))

# =============================================================================
# 3. CHECK the composite's signs against the saved per-sample composite
# =============================================================================
# The +/- column is printed on the panel, so it is checked, not copied. E34's
# composite is Blom(mean of sign x z(log expression)); the frames hold the
# Blom-scored genes, so the composite is NOT reconstructed exactly. The check is
# discriminating instead: the declared sign vector must correlate with the
# saved composite better than every vector with one sign flipped. Nothing
# computed here is plotted.
message("\n3. the composite's signs, checked against the saved composite")

SIGN_CHECK <- dplyr::bind_rows(lapply(COHORTS, function(coh) {
  f  <- e34$analysis_frames[[coh]]
  cc <- stats::complete.cases(f[c("config_comp", paste0("g_", names(COMPOSITE_SIGN)))])
  M  <- as.matrix(f[cc, paste0("g_", names(COMPOSITE_SIGN))])
  .rho <- function(s) stats::cor(f$config_comp[cc],
                                 as.numeric(M %*% s) / length(s),
                                 method = "spearman")
  flips <- vapply(seq_along(COMPOSITE_SIGN), function(i) {
    s <- COMPOSITE_SIGN; s[i] <- -s[i]; .rho(s)
  }, numeric(1))
  tibble::tibble(cohort = coh, n = sum(cc),
                 rho_declared = .rho(COMPOSITE_SIGN),
                 best_single_flip = max(flips),
                 flip_gene = names(COMPOSITE_SIGN)[which.max(flips)])
}))
.pr(SIGN_CHECK %>% dplyr::mutate(dplyr::across(dplyr::starts_with("rho"),
                                               ~ round(.x, 4)),
                                 best_single_flip = round(best_single_flip, 4)))
.assert("the declared signs fit the saved composite better than any single flip",
        all(SIGN_CHECK$rho_declared > SIGN_CHECK$best_single_flip) &&
          all(SIGN_CHECK$rho_declared > 0.95), n = nrow(SIGN_CHECK),
        detail = paste(sprintf("%s %.4f vs %.4f", SIGN_CHECK$cohort,
                               SIGN_CHECK$rho_declared,
                               SIGN_CHECK$best_single_flip), collapse = "; "))

# =============================================================================
# 4. The MYC contrasts, READ from E34 - for the text, not drawn
# =============================================================================
# The panel is read through its MYC pairs. These are E34's own contrasts, so
# the text can quote what the panel shows without anyone recomputing them.
message("\n4. E34's MYC contrasts (read, not recomputed; not drawn)")

MYC_CONTRASTS <- e34$contrasts %>%
  dplyr::filter(subtype == "all", readout_col == "config_comp",
                model %in% c("m1", "m3")) %>%
  dplyr::select(cohort, instrument, model, contrast, kind, est, lo, hi,
                ci_excludes_0) %>%
  dplyr::arrange(match(cohort, COHORTS), instrument, model, kind, contrast)
.pr(MYC_CONTRASTS %>% dplyr::filter(kind == "MYC") %>%
      dplyr::mutate(dplyr::across(c(est, lo, hi), ~ round(.x, 3))))
AGREE <- e34$quadrant$agreement
.pr(AGREE %>% dplyr::mutate(agreement = round(100 * agreement, 1)))

# =============================================================================
# 5. The figure
# =============================================================================
message("\n5. figure")

.heatmap <- function(inst, show_fill_legend = TRUE) {
  d <- CELLS %>% dplyr::filter(instrument == inst)
  hdr_coh <- tibble::tibble(
    x = vapply(COHORTS, function(c) mean(range(X_COL[[c]])), numeric(1)),
    label = paste0(COHORTS, " (n = ", .fmt_n(N_COH[COHORTS]), ")"))
  hdr_ox <- tibble::tibble(
    x = unlist(lapply(COHORTS, function(c) c(mean(X_COL[[c]][c("Q1", "Q3")]),
                                             mean(X_COL[[c]][c("Q2", "Q4")])))),
    label = rep(c("OXPHOS low", "OXPHOS high"), 2))
  # the thin rule between respiratory levels, inside each cohort
  rules <- tibble::tibble(
    x = vapply(COHORTS, function(c) mean(X_COL[[c]][c("Q3", "Q2")]), numeric(1)))
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = est),
                       width = 0.94, height = 0.92) +
    ggplot2::geom_segment(data = rules,
                          ggplot2::aes(x = x, xend = x, y = 0.5, yend = 13.1),
                          colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_text(data = hdr_coh, ggplot2::aes(x = x, y = Y_COHORT,
                                                    label = label),
                       size = TXT_MM, fontface = "bold") +
    ggplot2::geom_text(data = hdr_ox, ggplot2::aes(x = x, y = Y_OXLEV,
                                                   label = label),
                       size = SUB_MM, colour = "grey25") +
    ggplot2::scale_fill_gradient2(
      low = FILL_LOW, mid = "white", high = FILL_HIGH, midpoint = 0,
      limits = FILL_LIM, breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
      name = "Adjusted mean (within-cohort SD)",
      # title above the bar: beside it, it collided with the -0.6 label
      guide = ggplot2::guide_colourbar(theme = ggplot2::theme(
        legend.title.position = "top",
        legend.key.width = grid::unit(26, "mm"),
        legend.key.height = grid::unit(2, "mm")))) +
    ggplot2::scale_x_continuous(limits = X_LIM, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::labs(title = paste0("OXPHOS split by ", INSTRUMENTS[[inst]])) +
    ggplot2::theme_void(base_size = BASE_PT) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = BASE_PT + 0.5, face = "bold",
                                         hjust = 0.5,
                                         margin = ggplot2::margin(b = 2)),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.text = ggplot2::element_text(size = BASE_PT - 1))
  if (!show_fill_legend) p <- p + ggplot2::guides(fill = "none")
  p
}

.strip <- function(inst, show_y = TRUE, show_shape_legend = TRUE) {
  d <- STRIP %>% dplyr::filter(instrument == inst)
  brk <- unlist(lapply(COHORTS, function(c) X_COL[[c]][c("Q1", "Q3", "Q2", "Q4")]))
  # NO LINE JOINS ANY TWO QUADRANTS: no product term was fitted (E34).
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_linerange(ggplot2::aes(x = x, ymin = lo, ymax = hi),
                            colour = "grey15", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = est, shape = model),
                        colour = "grey15", fill = "white", size = 1.3,
                        stroke = 0.45) +
    ggplot2::scale_shape_manual(
      values = c(m1 = 16, m3 = 21),
      labels = c(m1 = "Composite, m1",
                 m3 = "Composite, m1 + proliferation (m3)"),
      name = NULL) +
    ggplot2::scale_x_continuous(limits = X_LIM, expand = c(0, 0),
                                breaks = unname(brk),
                                labels = rep(c("MYC\nlow", "MYC\nhigh"), 4)) +
    ggplot2::scale_y_continuous(limits = STRIP_LIM, expand = c(0, 0),
                                breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
    ggplot2::labs(x = NULL, y = if (show_y) "Composite\n(SD)" else NULL) +
    ggplot2::theme_classic(base_size = BASE_PT) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = BASE_PT - 1.5,
                                          colour = "grey25", lineheight = 0.9),
      axis.line.y = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.y = ggplot2::element_line(linewidth = 0.3),
      axis.text.y = if (show_y) ggplot2::element_text(size = BASE_PT - 1,
                                                      colour = "black")
                    else ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = BASE_PT - 0.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = BASE_PT - 0.5))
  if (!show_shape_legend) p <- p + ggplot2::guides(shape = "none")
  p
}

# The gene labels, the function blocks and the composite's sign column.
.labels <- function() {
  g <- tibble::tibble(gene = names(Y_GENE), y = unname(Y_GENE)) %>%
    dplyr::mutate(sign = dplyr::case_when(
      gene %in% names(COMPOSITE_SIGN) & COMPOSITE_SIGN[gene] > 0 ~ "+",
      gene %in% names(COMPOSITE_SIGN) & COMPOSITE_SIGN[gene] < 0 ~ "-",
      TRUE ~ ""))
  blocks <- tibble::tibble(
    label = c("Pro-apoptotic", "Anti-apoptotic"),
    y0 = c(min(Y_GENE[PRO]) - 0.4, min(Y_GENE[ANTI]) - 0.4),
    y1 = c(max(Y_GENE[PRO]) + 0.4, max(Y_GENE[ANTI]) + 0.4))
  ggplot2::ggplot(g) +
    ggplot2::geom_segment(data = blocks, ggplot2::aes(x = 0.16, xend = 0.16,
                                                      y = y0, yend = y1),
                          linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_text(data = blocks, ggplot2::aes(x = 0.07, y = (y0 + y1) / 2,
                                                   label = label),
                       angle = 90, size = SUB_MM, colour = "grey25") +
    ggplot2::geom_text(ggplot2::aes(x = 0.25, y = y, label = gene),
                       hjust = 0, size = TXT_MM, fontface = "italic") +
    ggplot2::geom_text(ggplot2::aes(x = 0.9, y = y, label = sign),
                       size = TXT_MM, fontface = "bold") +
    ggplot2::annotate("text", x = 0.9, y = Y_OXLEV, label = "sign",
                      size = SUB_MM, colour = "grey25") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    ggplot2::theme_void(base_size = BASE_PT)
}

W_LAB <- 0.85; W_HM <- 2.95; W_GAP <- 0.2
H_HM  <- 2.05; H_ST <- 0.85
.assemble <- function(insts) {
  plots <- list(.labels())
  widths <- W_LAB
  for (i in seq_along(insts)) {
    if (i > 1L) { plots <- c(plots, list(patchwork::plot_spacer())); widths <- c(widths, W_GAP) }
    plots <- c(plots, list(.heatmap(insts[i], show_fill_legend = i == 1L)))
    widths <- c(widths, W_HM)
  }
  plots <- c(plots, list(patchwork::plot_spacer()))
  for (i in seq_along(insts)) {
    if (i > 1L) plots <- c(plots, list(patchwork::plot_spacer()))
    plots <- c(plots, list(.strip(insts[i], show_y = i == 1L,
                                  show_shape_legend = i == 1L)))
  }
  patchwork::wrap_plots(plots, nrow = 2, byrow = TRUE, widths = widths,
                        heights = c(H_HM, H_ST)) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}
FIG_BOTH <- .assemble(c("gsva", "mitopps"))
FIG_GSVA <- .assemble("gsva")
FIG_MPPS <- .assemble("mitopps")

# --- every plotted value is its stored value, read back from the build ----
.layer_of <- function(p, geom) which(vapply(p$layers, function(l)
  class(l$geom)[1], character(1)) == geom)

.check_heatmap <- function(inst) {
  p  <- .heatmap(inst)
  b  <- ggplot2::ggplot_build(p)
  k  <- .layer_of(p, "GeomTile")
  .stop_if(length(k) == 1L, "cannot identify the tile layer")
  dt <- b$data[[k]]
  d  <- CELLS %>% dplyr::filter(instrument == inst)
  m  <- match(paste(round(dt$x, 6), round(dt$y, 6)),
              paste(round(d$x, 6), round(d$y, 6)))
  .assert(paste0(INSTRUMENTS[[inst]], " heatmap: 96 tiles drawn, one per cell"),
          nrow(dt) == 96L && nrow(d) == 96L && !anyNA(m) && !any(duplicated(m)),
          n = nrow(dt))
  # the colour drawn is the scale's colour for the stored value
  sc  <- b$plot$scales$get_scales("fill")
  exp_fill <- sc$map(d$est[m])
  .assert(paste0(INSTRUMENTS[[inst]], " heatmap: every tile's colour is the ",
                 "scale's colour for its stored value"),
          identical(as.character(dt$fill), as.character(exp_fill)),
          n = nrow(dt))
}
.check_strip <- function(inst) {
  p  <- .strip(inst)
  b  <- ggplot2::ggplot_build(p)
  geo <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  .assert(paste0(INSTRUMENTS[[inst]], " strip: no line joins two quadrants"),
          all(geo %in% c("GeomHline", "GeomLinerange", "GeomPoint")),
          n = length(geo), detail = paste(geo, collapse = ", "))
  dp <- b$data[[.layer_of(p, "GeomPoint")]]
  dl <- b$data[[.layer_of(p, "GeomLinerange")]]
  d  <- STRIP %>% dplyr::filter(instrument == inst)
  mp <- match(round(dp$x, 6), round(d$x, 6))
  ml <- match(round(dl$x, 6), round(d$x, 6))
  .assert(paste0(INSTRUMENTS[[inst]], " strip: 16 points and 16 intervals, ",
                 "one per quadrant x model"),
          nrow(dp) == 16L && nrow(dl) == 16L && !anyNA(c(mp, ml)) &&
            !any(duplicated(mp)) && !any(duplicated(ml)), n = nrow(dp))
  dev <- max(.dev(dp$y, d$est[mp]), .dev(dl$ymin, d$lo[ml]),
             .dev(dl$ymax, d$hi[ml]))
  .assert(paste0(INSTRUMENTS[[inst]], " strip: drawn values equal the stored ",
                 "adjusted means and intervals"),
          dev == 0, n = 3L * nrow(dp), max_dev = dev)
}
for (inst in names(INSTRUMENTS)) {
  .check_heatmap(inst)
  .check_strip(inst)
}

.ensure_dir(DIR_FIGURES)
.ensure_dir(DIR_DOC_FIGURES)
.save_fig <- function(fig, stem, w, h) {
  ggplot2::ggsave(paste0(stem, ".png"), fig, width = w, height = h, dpi = 300,
                  bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), fig, width = w, height = h,
                  device = grDevices::cairo_pdf, bg = "white")
}
H_FIG  <- H_HM + H_ST + 0.55          # plus titles and the collected legend
W_BOTH <- W_LAB + 2 * W_HM + W_GAP
W_ONE  <- W_LAB + W_HM
.save_fig(FIG_BOTH, PATH_FIG_BOTH, W_BOTH, H_FIG)
.save_fig(FIG_BOTH, PATH_FIG_DOC,  W_BOTH, H_FIG)
.save_fig(FIG_GSVA, PATH_FIG_GSVA, W_ONE,  H_FIG + 0.2)
.save_fig(FIG_MPPS, PATH_FIG_MPPS, W_ONE,  H_FIG + 0.2)
message("   combined ", sprintf("%.2f x %.2f", W_BOTH, H_FIG), " in (",
        round(W_BOTH * 25.4), " x ", round(H_FIG * 25.4), " mm)")
message("   ", PATH_FIG_BOTH, ".{png,pdf}")
message("   ", PATH_FIG_GSVA, ".{png,pdf}")
message("   ", PATH_FIG_MPPS, ".{png,pdf}")
message("   ", PATH_FIG_DOC, ".{png,pdf}  (tracked)")

# =============================================================================
# 6. Legend draft - DESCRIPTIVE ONLY
# =============================================================================
message("\n6. legend draft")

rng <- CELLS %>% dplyr::distinct(instrument, cohort, q, n_cell) %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(lo = min(n_cell), hi = max(n_cell), .groups = "drop")
rng <- stats::setNames(sprintf("%s to %s", .fmt_n(rng$lo), .fmt_n(rng$hi)),
                       rng$cohort)
agr <- stats::setNames(sprintf("%.1f%%", 100 * AGREE$agreement), AGREE$cohort)

LEGEND <- paste0(
  "Expression of twelve BCL-2-family genes (rows) and of a six-gene signed ",
  "composite (strip) across four groups of tumours defined by the cohort ",
  "medians of MYC activity (M_a signature, GSVA) and of the OXPHOS subunit ",
  "score, in TCGA-BRCA (n = ", .fmt_n(N_COH[["TCGA"]]), ") and SCAN-B (n = ",
  .fmt_n(N_COH[["SCAN-B"]]), "). Left, the OXPHOS score by GSVA; right, by ",
  "mitoPPS; MYC activity is the same in both. The two instruments assign the ",
  "same group to ", agr[["TCGA"]], " of TCGA and ", agr[["SCAN-B"]], " of ",
  "SCAN-B tumours. Columns: within each cohort, OXPHOS-low then OXPHOS-high ",
  "tumours, MYC-low then MYC-high within each; group sizes ", rng[["TCGA"]],
  " (TCGA) and ", rng[["SCAN-B"]], " (SCAN-B). Cells, adjusted group means ",
  "of Blom-transformed expression in within-cohort standard deviations, from ",
  "a linear model with group, PAM50 subtype, tumour purity and leukocyte ",
  "fraction in TCGA, and group and PAM50 subtype in SCAN-B, which has no ",
  "purity or leukocyte estimate (m1). Strip, the composite (mean z of BBC3, ",
  "BID, BIK, BAD and BCL2L1, and of MCL1 with its sign reversed; signs at ",
  "left): filled, m1; open, m1 with the proliferation score added (m3); ",
  "bars, 95% confidence intervals. The composite's genes and signs were ",
  "selected by their association with OXPHOS in these cohorts (Methods). ",
  "Rows are ordered by function, composite members first, and not by the ",
  "data.")
message("   ", .sub(LEGEND))

# =============================================================================
# 7. Save - the plotting frames, per R rule 4
# =============================================================================
message("\n7. saving")

ASSERTIONS <- dplyr::bind_rows(ASSERT)
.pr(ASSERTIONS %>% dplyr::transmute(check, n_compared,
                                    max_abs_dev = signif(max_abs_dev, 3), pass))
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")

readr::write_csv(CELLS,         file.path(DIR_TABLES, "E46_heatmap_cells.csv"))
readr::write_csv(STRIP,         file.path(DIR_TABLES, "E46_strip.csv"))
readr::write_csv(MYC_CONTRASTS, file.path(DIR_TABLES, "E46_contrasts_from_e34.csv"))
readr::write_csv(ASSERTIONS,    file.path(DIR_TABLES, "E46_assertions.csv"))

saveRDS(list(
  heatmap_cells  = CELLS,
  strip          = STRIP,
  contrasts_e34  = MYC_CONTRASTS,
  agreement_e34  = AGREE,
  sign_check     = SIGN_CHECK,
  analysis_n     = nset,
  assertions     = ASSERTIONS,
  legend_draft   = LEGEND,
  # FILE NAMES, NOT PATHS, so this element digests the same in a redirected
  # dry run and in the author's run. Figures are not digested.
  figure_files   = list(
    outputs = basename(paste0(c(PATH_FIG_BOTH, PATH_FIG_GSVA, PATH_FIG_MPPS),
                              rep(c(".png", ".pdf"), each = 3))),
    tracked_in_docs_figures = basename(paste0(PATH_FIG_DOC, c(".png", ".pdf")))),
  spec = list(
    posture = paste0(
      "EXPLORATORY and DISPLAY ONLY. Every plotted value is an adjusted group ",
      "mean from E34's $groups (subtype 'all', m1 cells, m1 and m3 strip). ",
      "Nothing is fitted."),
    decisions = paste0(
      "Author, 2026-10-22: rows blocked by function (pro-, then ",
      "anti-apoptotic; composite members first), and BOTH instruments in the ",
      "main figure, side by side. Proposal section 8."),
    snapshot = paste0(
      "data/from_validation/, through E34. The joint snapshot is not read."),
    circularity = paste0(
      "The composite's genes and signs were read off OXPHOS in these cohorts, ",
      "so the OXPHOS separation is near-guaranteed (E34 1.2). The informative ",
      "contrasts are the MYC ones, which sit side by side."),
    no_product = paste0(
      "No product term was fitted, and no line joins two quadrants; asserted ",
      "on the built strip."),
    instruments = paste0(
      "The two panels split OXPHOS differently and agree on the quadrant for ",
      "85.0% of TCGA and 82.7% of SCAN-B tumours, so a column is not the same ",
      "tumours in both panels."),
    sign_check = paste0(
      "The composite's printed signs are CHECKED against E34's saved ",
      "per-sample composite: the declared vector must correlate better than ",
      "any single flip. A check, not an estimate; nothing from it is plotted."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  built = Sys.time()), PATH_E46)

message("   ", PATH_E46)
message("   4 tables in ", DIR_TABLES)
message("\nE46 done. Display only; nothing here is interpreted.")
message("THE OXPHOS SEPARATION IS NEAR-GUARANTEED BY CONSTRUCTION. ",
        "Read the MYC pairs.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E46)

  # every check, with how much it compared
  x$assertions %>% as.data.frame()

  # the composite's signs, checked against the saved per-sample composite
  x$sign_check %>% as.data.frame()

  # what the MYC pairs show, as E34 estimated it - for the text
  x$contrasts_e34 %>% dplyr::filter(kind == "MYC") %>% as.data.frame()

  # one instrument's cells as a gene x column table
  x$heatmap_cells %>% dplyr::filter(instrument == "mitopps") %>%
    dplyr::transmute(gene, col = paste(cohort, q), est = round(est, 2)) %>%
    tidyr::pivot_wider(names_from = col, values_from = est) %>% as.data.frame()

  cat(x$legend_draft, "\n")

  FIG_BOTH
}
