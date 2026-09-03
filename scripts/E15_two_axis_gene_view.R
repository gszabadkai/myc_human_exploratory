# E15_two_axis_gene_view.R
# =============================================================================
# THE TWO AXES ON ONE ROW PER GENE.
#
# The 44 canonical apoptosis genes have been drawn against OXPHOS (E08 fig6,
# E10 fig1) and against MYC (E10 fig2) on SEPARATE panels with a shared row
# order, and as a cloud with one axis each way (E11 fig1). None of those shows
# WHICH GENE MOVES between the axes, which is the thing the author is asking
# about: the machinery's mitochondrial and cytosolic halves behave differently
# on OXPHOS, and the question is whether MYC sees the same split at all.
#
# This is a DISPLAY script. It computes no correlation. Every rho it draws was
# computed in E11 and is read from `results/prolif_adjusted_machinery.rds`;
# section 1 asserts that the split it recovers from those numbers is bit-equal
# to the value E14 published for the same 44 genes. The only new quantity is a
# per-gene DIFFERENCE between the two axes, and its one real hazard is below.
#
# =============================================================================
# THE HAZARD: A DIFFERENCE OF TWO CORRELATIONS AGAINST TWO DIFFERENT SCORES
# =============================================================================
# `rho(gene, OXPHOS) - rho(gene, MYC)` is not a difference in a common unit.
# The 44 per-gene correlations are SPREAD WIDER on the OXPHOS axis than on the
# MYC axis - E11 already reported this as the shape of its figure 1 cloud, and
# the SDs are printed in section 3 here. A gene therefore earns a positive gap
# partly by sitting on the more-spread axis, before any biology.
#
# So the gap is computed twice:
#   `gap`    rho_OXPHOS - rho_MYC. Interpretable, and the units are the ones
#            every other figure in the study uses.
#   `gap_z`  the same difference after dividing each axis by the SD of ITS OWN
#            44 values, within cohort and within adjustment. This asks only
#            whether a gene sits higher in the OXPHOS ORDERING than in the MYC
#            ordering, and is immune to the spread difference.
#
# If the regulon split of `gap` were pure scale artefact, `gap_z` would
# split at zero. HONESTY NOTE: this check was run while the figure was being
# designed, so its answer is known and is stated in section 3 and on figure 3.
# It is a caveat quantified, not a prediction tested. What it is NOT is a
# result: it is a check on whether figure 1 is allowed to be drawn at all.
#
# =============================================================================
# WHAT THE RED HALF IS, AND IT IS NOT A PLACE
# =============================================================================
# The study retired localisation language on 2026-09-02. MitoCarta is a
# PROTEOME catalogue and half the BCL2 family translocates - BAX is cytosolic
# until activated, BID must be cleaved - so at transcript level membership
# marks MEMBERSHIP OF THE NUCLEAR-ENCODED MITOCHONDRIAL REGULON, not where the
# protein sits. The sub-compartment labels (MOM, IMS, MIM) are carried in the
# saved table because E10 found a gradient by depth, but they are MitoCarta's
# annotation of the protein and they are not a claim about the transcript.
# "mitochondrial half" below always means "in MitoCarta 3.0".
#
# =============================================================================
# WHAT THE FIGURES MAY AND MAY NOT BE READ AS SAYING
# =============================================================================
#   MAY:  which genes carry the OXPHOS ordering, whether the two cohorts agree
#         on each one, and whether MitoCarta membership tracks the difference.
#   MAY NOT: anything about a single gene. These 44 were not selected on the
#         statistic shown, which is better than E08's driver lists, but 44
#         genes x 2 cohorts is still a grid and no cell of it is a finding.
#         Five of the 44 sit below the 25th expression percentile in at least
#         one cohort and are STARRED; their rho is largely quantisation noise
#         against a score.
#   MAY NOT: be read as MYC having no relation to these genes. E11 figure 2 is
#         the control that says the adjusted MYC axis still carries signal
#         (against the mitoribosome), and it must be read with this.
#
# EXPLORATORY. Nothing here is pre-registered.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

PATH_E08     <- file.path(DIR_RESULTS, "strata_and_death_genes.rds")
PATH_E11     <- file.path(DIR_RESULTS, "prolif_adjusted_machinery.rds")
PATH_E14     <- file.path(DIR_RESULTS, "curated_comparators.rds")
PATH_E15     <- file.path(DIR_RESULTS, "two_axis_gene_view.rds")
PATH_E15_CSV <- file.path(DIR_TABLES,  "E15_gene_rho_two_axes.csv")

AXES          <- c("MYC", "OXPHOS")
ADJ_LEVELS    <- c("raw", "adj. PROLIF_DISJOINT")
ADJ_MAIN      <- "adj. PROLIF_DISJOINT"
COHORT_LEVELS <- c("TCGA", "SCAN-B")
LOW_EXPR_PCT  <- 0.25
N_CANON       <- 44L

.ensure_dir(DIR_FIGURES); .ensure_dir(DIR_TABLES); .ensure_dir(DIR_RESULTS)

# =============================================================================
# 1. Read the three upstream objects, and prove they still say what they said
# =============================================================================
message("\n1. upstream objects")
for (p in c(PATH_E08, PATH_E11, PATH_E14))
  if (!file.exists(p)) stop("missing: ", p, "\nRun its script first.",
                            call. = FALSE)
e08 <- readRDS(PATH_E08)
e11 <- readRDS(PATH_E11)
e14 <- readRDS(PATH_E14)

MYC_REF     <- e11$settings$myc_axis
PROLIF_COV  <- e11$settings$prolif_covariate
stopifnot(identical(MYC_REF, "FELSHER__MITOSTRIP"),
          identical(PROLIF_COV, "PROLIF_DISJOINT"))
message("   MYC axis ", MYC_REF, " | proliferation covariate ", PROLIF_COV)

gene_tab <- e11$gene_tab %>%
  dplyr::filter(axis %in% AXES, adjustment %in% ADJ_LEVELS)
stopifnot(nrow(gene_tab) == N_CANON * 2L * 2L * 2L,
          dplyr::n_distinct(gene_tab$gene) == N_CANON,
          setequal(unique(gene_tab$cohort), COHORT_LEVELS),
          is.logical(gene_tab$mitocarta), !anyNA(gene_tab$rho))
CANON <- sort(unique(gene_tab$gene))
message("   E11 gene_tab: ", N_CANON, " genes x 2 cohorts x 2 axes x 2 ",
        "adjustments = ", nrow(gene_tab), " rows")

# THE CROSS-SCRIPT ANCHOR. E14 published the localisation split for these same
# 44 genes as the target of its comparator test. Recomputing it here from E11's
# per-gene numbers must reproduce it exactly, or the two scripts are looking at
# different objects and nothing below is comparable to the rest of the study.
split_here <- gene_tab %>%
  dplyr::group_by(cohort, axis, adjustment) %>%
  dplyr::summarise(split = stats::cor(rho, as.numeric(mitocarta),
                                      method = "spearman"), .groups = "drop")
anchor <- e14$comparator_splits %>%
  dplyr::filter(set == "apoptotic machinery (44)") %>%
  dplyr::select(cohort, axis, adjustment, split_E14 = split)
anchor_chk <- dplyr::inner_join(split_here, anchor,
                                by = c("cohort", "axis", "adjustment"))
stopifnot(nrow(anchor_chk) == 8L,
          max(abs(anchor_chk$split - anchor_chk$split_E14)) < 1e-10)
message("   anchor OK: the split recomputed here is bit-equal to E14's for all",
        " 8 cohort x axis x adjustment cells")

# --- MitoCarta sub-compartment, for the mitochondrial half only --------------
# Membership is the split every other script uses. The sub-compartment is
# carried only as a label, because E10 found a gradient by depth into the
# organelle and a per-gene figure is the natural place to see it.
mitocarta_sheet <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
stopifnot(all(c("Symbol", "MitoCarta3.0_SubMitoLocalization") %in%
                names(mitocarta_sheet)))
SUBMITO <- stats::setNames(
  mitocarta_sheet[["MitoCarta3.0_SubMitoLocalization"]], mitocarta_sheet$Symbol)

# --- expression rank, from E08, so the noisy genes can be starred ------------
expr <- e08$expr_rank %>% dplyr::filter(gene %in% CANON)
stopifnot(dplyr::n_distinct(expr$gene) == N_CANON,
          all(expr$low_expression == (expr$expr_pct < LOW_EXPR_PCT)))
low_expr <- expr %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(low_expression = any(low_expression),
                   min_expr_pct = min(expr_pct), .groups = "drop")
LOW_GENES <- sort(low_expr$gene[low_expr$low_expression])
message("   below the ", LOW_EXPR_PCT * 100, "th expression percentile in at ",
        "least one cohort, starred on every panel: ",
        paste(LOW_GENES, collapse = ", "))

# =============================================================================
# 2. The two-axis gene table
# =============================================================================
message("\n2. the gene table")

annot <- gene_tab %>%
  dplyr::distinct(gene, mitocarta, effect, cdc_module, in_prolif_cov,
                  in_oxphos_arm) %>%
  dplyr::mutate(submito = unname(SUBMITO[gene]))
stopifnot(nrow(annot) == N_CANON, !anyNA(annot$submito[annot$mitocarta]))
# `mito_class` is a FACTOR with the mitochondrial half first, because figure 2
# stacks the two halves and ggplot orders facet strips by the factor. Left as a
# character it sorts alphabetically, puts the cytosolic block on top, and runs
# the opposite way to figure 1's sort while claiming to share its row order.
MITO_LEVELS <- c("mitochondrial (MitoCarta 3.0)", "cytosolic (not in MitoCarta)")
annot <- annot %>%
  dplyr::mutate(
    compartment = dplyr::if_else(mitocarta, submito, "cytosolic"),
    mito_class  = factor(dplyr::if_else(mitocarta, MITO_LEVELS[1],
                                        MITO_LEVELS[2]), levels = MITO_LEVELS)) %>%
  dplyr::left_join(low_expr, by = "gene")
stopifnot(!anyNA(annot$mito_class))
# The sub-compartment counts are E14's, and a change here means the MitoCarta
# snapshot moved under the study rather than that a gene was reclassified.
CMP_EXPECT <- c(MOM = 13L, IMS = 5L, MIM = 2L)
got <- table(annot$compartment[annot$mitocarta])
stopifnot(setequal(names(got), names(CMP_EXPECT)),
          all(as.integer(got[names(CMP_EXPECT)]) == CMP_EXPECT))
message("   mitochondrial half: 13 MOM, 5 IMS, 2 MIM | cytosolic half: ",
        sum(!annot$mitocarta))

# Per-axis standardisation is WITHIN cohort and WITHIN adjustment, so the
# denominator is the spread of the very 44 numbers being differenced.
gene_axes <- gene_tab %>%
  dplyr::select(cohort, axis, adjustment, gene, rho) %>%
  dplyr::group_by(cohort, axis, adjustment) %>%
  dplyr::mutate(sd_axis = stats::sd(rho), rho_z = rho / sd_axis) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(annot, by = "gene")

pairs_tab <- gene_axes %>%
  dplyr::select(cohort, adjustment, gene, axis, rho, rho_z) %>%
  tidyr::pivot_wider(names_from = axis, values_from = c(rho, rho_z)) %>%
  dplyr::mutate(gap = rho_OXPHOS - rho_MYC, gap_z = rho_z_OXPHOS - rho_z_MYC) %>%
  dplyr::left_join(annot, by = "gene")
stopifnot(nrow(pairs_tab) == N_CANON * 2L * 2L, !anyNA(pairs_tab$gap))

# One row per gene: the cross-cohort mean that fixes the row order, and whether
# the two cohorts agree on the SIGN of the gap. Row order is taken from the
# ADJUSTED numbers and is reused on every figure, so the panels can be laid
# beside one another and read for movement.
gene_order <- pairs_tab %>%
  dplyr::filter(adjustment == ADJ_MAIN) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(mean_gap = mean(gap), mean_gap_z = mean(gap_z),
                   mean_OXPHOS = mean(rho_OXPHOS), mean_MYC = mean(rho_MYC),
                   gap_TCGA = gap[cohort == "TCGA"],
                   gap_SCANB = gap[cohort == "SCAN-B"],
                   .groups = "drop") %>%
  dplyr::mutate(agree = sign(gap_TCGA) == sign(gap_SCANB)) %>%
  dplyr::left_join(annot, by = "gene")
stopifnot(nrow(gene_order) == N_CANON)

# The label carries every flag a reader needs to distrust a row, because a row
# is where the distrust has to live: a caption is read once and a gene name is
# read forty-four times.
#   *     below the 25th expression percentile in at least one cohort
#   [ox]  in the 89-gene OXPHOS subunits arm, so partly correlated with itself
#   [p]   in the 318-gene proliferation covariate, so partly adjusted for itself
gene_order <- gene_order %>%
  dplyr::mutate(glab = paste0(gene,
                              ifelse(low_expression, " *", ""),
                              ifelse(in_oxphos_arm, " [ox]", ""),
                              ifelse(in_prolif_cov, " [p]", "")))
LAB_LEVELS <- gene_order %>% dplyr::arrange(mean_gap) %>% dplyr::pull(glab)
stopifnot(!anyDuplicated(LAB_LEVELS), length(LAB_LEVELS) == N_CANON)
GLAB <- stats::setNames(gene_order$glab, gene_order$gene)

.add_lab <- function(d) {
  d$glab <- factor(unname(GLAB[d$gene]), levels = LAB_LEVELS)
  stopifnot(!anyNA(d$glab))
  d
}
gene_axes <- .add_lab(gene_axes)
pairs_tab <- .add_lab(pairs_tab)
gene_order <- .add_lab(gene_order)

# =============================================================================
# 3. What the table says, printed
# =============================================================================
message("\n3. summaries")

message("\n   3.1 SD of the 44 per-gene correlations, per axis. This is the",
        "\n       denominator of gap_z and the reason gap_z exists:")
sd_tab <- gene_axes %>%
  dplyr::distinct(cohort, axis, adjustment, sd_axis) %>%
  dplyr::arrange(cohort, adjustment, axis)
sd_tab %>% dplyr::mutate(sd_axis = round(sd_axis, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   3.2 the localisation split, on each axis and on the difference:")
split_tab <- dplyr::bind_rows(
  split_here %>% dplyr::rename(value = split) %>%
    dplyr::mutate(quantity = axis) %>% dplyr::select(-axis),
  pairs_tab %>% dplyr::group_by(cohort, adjustment) %>%
    dplyr::summarise(
      `OXPHOS - MYC` = stats::cor(gap, as.numeric(mitocarta),
                                  method = "spearman"),
      `OXPHOS - MYC, standardised` = stats::cor(gap_z, as.numeric(mitocarta),
                                                method = "spearman"),
      .groups = "drop") %>%
    tidyr::pivot_longer(dplyr::starts_with("OXPHOS - MYC"),
                        names_to = "quantity", values_to = "value")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_LEVELS),
                adjustment = factor(adjustment, levels = ADJ_LEVELS),
                quantity = factor(quantity,
                  levels = c("MYC", "OXPHOS", "OXPHOS - MYC",
                             "OXPHOS - MYC, standardised"))) %>%
  dplyr::arrange(adjustment, quantity, cohort)
split_tab %>% dplyr::mutate(value = round(value, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   Read the two middle rows against each other. `OXPHOS - MYC` is",
        " essentially\n   the OXPHOS row unchanged: subtracting MYC takes",
        " almost nothing out of the\n   regulon ordering, because MYC",
        " barely carries any of it.")

message("\n   3.3 medians by half. All 44 genes, TP53 and BIRC5 INCLUDED",
        "\n       (E11's own s6_adj excludes them and so runs about 0.04",
        " higher):")
half_tab <- pairs_tab %>%
  dplyr::group_by(cohort, adjustment, mito_class) %>%
  dplyr::summarise(n = dplyr::n(),
                   med_MYC = stats::median(rho_MYC),
                   med_OXPHOS = stats::median(rho_OXPHOS),
                   med_gap = stats::median(gap),
                   med_gap_z = stats::median(gap_z), .groups = "drop")
half_tab %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   3.4 does the gap itself replicate across the two cohorts?")
rep_tab <- pairs_tab %>%
  dplyr::select(cohort, adjustment, gene, gap) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = gap) %>%
  dplyr::group_by(adjustment) %>%
  dplyr::summarise(spearman = stats::cor(TCGA, `SCAN-B`, method = "spearman"),
                   pearson = stats::cor(TCGA, `SCAN-B`),
                   n_sign_agree = sum(sign(TCGA) == sign(`SCAN-B`)),
                   .groups = "drop")
rep_tab %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
DISAGREE <- sort(gene_order$gene[!gene_order$agree])
message("   disagreeing in sign under ", ADJ_MAIN, ", drawn with a DOTTED",
        " bar on figure 1:\n   ", paste(DISAGREE, collapse = ", "))

message("\n   3.5 the ten largest and ten smallest gaps (adjusted, cross-cohort",
        " mean).\n       SELECTED ON THE STATISTIC SHOWN - a description of the",
        " ordering, not a\n       finding about any gene in it:")
ends <- dplyr::bind_rows(
  gene_order %>% dplyr::arrange(dplyr::desc(mean_gap)) %>% utils::head(10) %>%
    dplyr::mutate(end = "largest gap: ordered by OXPHOS, not by MYC"),
  gene_order %>% dplyr::arrange(mean_gap) %>% utils::head(10) %>%
    dplyr::mutate(end = "smallest gap: MYC as high or higher"))
ends %>%
  dplyr::transmute(end, gene, compartment, effect,
                   MYC = round(mean_MYC, 3), OXPHOS = round(mean_OXPHOS, 3),
                   gap = round(mean_gap, 3), agree) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4. Figures
# =============================================================================
message("\n4. figures")

COMP_COLS  <- c(`mitochondrial (MitoCarta 3.0)` = "#d7191c",
                `cytosolic (not in MitoCarta)`  = "grey55")
stopifnot(setequal(names(COMP_COLS), MITO_LEVELS))
AXIS_FILL  <- c(MYC = "white", OXPHOS = "#1f4e79")
theme_e15 <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey92",
                                                          colour = NA),
                 legend.position = "bottom",
                 legend.box = "vertical",
                 legend.spacing.y = ggplot2::unit(1, "pt"),
                 plot.caption = ggplot2::element_text(size = 7.4,
                                                      colour = "grey30",
                                                      hjust = 0))
.save <- function(p, nm, w, h) {
  for (ext in c("png", "pdf"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0(nm, ".", ext)), p,
                    width = w, height = h, dpi = 300)
  message("   ", nm); invisible(p)
}

# Numbers that appear in a caption are pasted from the objects above, never
# typed, so a caption cannot drift away from the figure it sits under.
.n <- function(d, ...) {
  v <- dplyr::filter(d, ...)
  stopifnot(nrow(v) == 1L)
  v
}
.sp <- function(coh, q, adj = ADJ_MAIN)
  .n(split_tab, cohort == coh, quantity == q, adjustment == adj)$value
.hf <- function(coh, cls, col, adj = ADJ_MAIN)
  .n(half_tab, cohort == coh, mito_class == cls, adjustment == adj)[[col]]
MITO_CLS <- "mitochondrial (MitoCarta 3.0)"
CYTO_CLS <- "cytosolic (not in MitoCarta)"

# --- FIG 1: the dumbbell, sorted by the gap ---------------------------------
# One row per gene, ordered by how far it moves between the axes. The row order
# is DATA-DRIVEN and the regulon colour is not used to sort, so the fact that
# the red rows collect at the top is a property of the numbers rather than of
# the layout. That is the whole reason not to facet by membership here.
f1_pts <- gene_axes %>% dplyr::filter(adjustment == ADJ_MAIN) %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_LEVELS),
                axis = factor(axis, levels = AXES))
f1_seg <- pairs_tab %>% dplyr::filter(adjustment == ADJ_MAIN) %>%
  dplyr::left_join(dplyr::select(gene_order, gene, agree), by = "gene") %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_LEVELS))
# Compartment is encoded TWICE and neither encoding is a theme hack. A
# vectorised `element_text(colour = ...)` would colour the 44 gene names, but
# ggplot2 warns that vectorised input is unsupported and its recycling order is
# not guaranteed - a figure that silently re-orders that vector would mislabel
# every row and still render. So membership is drawn as DATA: a square in the
# left margin of each panel, on the same colour scale as the bars.
X_RANGE <- range(c(f1_seg$rho_MYC, f1_seg$rho_OXPHOS))
X_MARK  <- X_RANGE[1] - diff(X_RANGE) * 0.075

CAP1 <- sprintf(paste0(
  "EACH BAR IS ONE GENE AND ITS LENGTH IS THE DIFFERENCE BETWEEN THE TWO AXES.",
  " The hollow point is\n",
  "the gene's correlation with MYC activity, the filled point its correlation",
  " with OXPHOS; rows are\n",
  "sorted by the cross-cohort mean of the difference, and the SAME ORDER is",
  " used in both panels and\n",
  "on figure 2. The square in the left margin and the bar itself both say the",
  " same thing: red = in\n",
  "MitoCarta 3.0. The sort does not know about the colour.\n",
  "THE ORDERING IS AN OXPHOS PROPERTY. Spearman of the value with MitoCarta",
  " membership across the\n",
  "44: OXPHOS %.2f / %.2f (TCGA / SCAN-B), MYC %.2f / %.2f, and the DIFFERENCE",
  " %.2f / %.2f - subtracting\n",
  "MYC leaves that ordering almost exactly where it was, because MYC hardly",
  " carries it.\n",
  "RED IS NOT A PLACE. MitoCarta is a proteome catalogue and half the BCL2",
  " family translocates, so at\n",
  "transcript level membership marks the NUCLEAR-ENCODED MITOCHONDRIAL REGULON",
  " rather than where the\n",
  "protein sits.\n",
  "Median difference: mitochondrial %+.2f / %+.2f, cytosolic %+.2f / %+.2f.\n",
  "READ THE LENGTHS WITH CARE. The 44 correlations spread wider on OXPHOS (SD",
  " %.2f / %.2f) than on\n",
  "MYC (%.2f / %.2f), so part of every bar is that difference in spread and",
  " not biology. Figure 3\n",
  "puts a number on how much: standardising each axis by its own SD takes the",
  " split from %.2f / %.2f\n",
  "to %.2f / %.2f - smaller, and still far above the MYC row.\n",
  "A DOTTED bar is a gene whose two cohorts disagree on the SIGN of the",
  " difference (%d of 44).\n",
  "* = below the 25th expression percentile in at least one cohort, so its rho",
  " is largely\n",
  "quantisation noise. [ox] = CYCS, which is one of the 89 OXPHOS-arm genes and",
  " is partly\n",
  "correlated with itself. [p] = in the proliferation covariate, so partly",
  " adjusted for itself.\n",
  "THIS FIGURE IS NOT EVIDENCE THAT MYC IS UNRELATED TO THESE GENES. E11",
  " figure 2 is the control\n",
  "that shows the adjusted MYC axis still tracks the mitoribosome; without it",
  " an empty MYC column\n",
  "and an emptied MYC score look identical."),
  .sp("TCGA", "OXPHOS"), .sp("SCAN-B", "OXPHOS"),
  .sp("TCGA", "MYC"), .sp("SCAN-B", "MYC"),
  .sp("TCGA", "OXPHOS - MYC"), .sp("SCAN-B", "OXPHOS - MYC"),
  .hf("TCGA", MITO_CLS, "med_gap"), .hf("SCAN-B", MITO_CLS, "med_gap"),
  .hf("TCGA", CYTO_CLS, "med_gap"), .hf("SCAN-B", CYTO_CLS, "med_gap"),
  .n(sd_tab, cohort == "TCGA", axis == "OXPHOS",
     adjustment == ADJ_MAIN)$sd_axis,
  .n(sd_tab, cohort == "SCAN-B", axis == "OXPHOS",
     adjustment == ADJ_MAIN)$sd_axis,
  .n(sd_tab, cohort == "TCGA", axis == "MYC", adjustment == ADJ_MAIN)$sd_axis,
  .n(sd_tab, cohort == "SCAN-B", axis == "MYC", adjustment == ADJ_MAIN)$sd_axis,
  .sp("TCGA", "OXPHOS - MYC"), .sp("SCAN-B", "OXPHOS - MYC"),
  .sp("TCGA", "OXPHOS - MYC, standardised"),
  .sp("SCAN-B", "OXPHOS - MYC, standardised"),
  length(DISAGREE))

g1 <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey20") +
  ggplot2::geom_segment(
    data = f1_seg,
    ggplot2::aes(x = rho_MYC, xend = rho_OXPHOS, y = glab, yend = glab,
                 colour = mito_class, linetype = agree),
    linewidth = 0.9, alpha = 0.85) +
  ggplot2::geom_point(
    data = f1_pts,
    ggplot2::aes(x = rho, y = glab, fill = axis),
    shape = 21, size = 2.1, colour = "grey15", stroke = 0.45) +
  ggplot2::geom_point(
    data = gene_order,
    ggplot2::aes(x = X_MARK, y = glab, colour = mito_class),
    shape = 15, size = 1.8, show.legend = FALSE) +
  ggplot2::facet_grid(. ~ cohort) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.035, 0.03))) +
  ggplot2::scale_colour_manual(values = COMP_COLS, name = NULL) +
  ggplot2::scale_linetype_manual(
    values = c(`TRUE` = "solid", `FALSE` = "dotted"),
    breaks = c("FALSE"), labels = c("cohorts disagree on the sign"),
    name = NULL) +
  ggplot2::scale_fill_manual(
    values = AXIS_FILL, name = NULL,
    labels = c(MYC = paste0("MYC activity (", MYC_REF, ")"),
               OXPHOS = "OXPHOS subunits")) +
  ggplot2::guides(fill = ggplot2::guide_legend(order = 1,
                    override.aes = list(size = 2.6)),
                  colour = ggplot2::guide_legend(order = 2),
                  linetype = ggplot2::guide_legend(order = 3)) +
  ggplot2::labs(
    title = "How far each apoptotic gene moves between the OXPHOS and MYC axes",
    subtitle = paste0("EXPLORATORY - not pre-registered | 44 canonical genes | ",
                      "partial Spearman on ", PROLIF_COV, " (318 genes)"),
    x = "per-gene Spearman rho with the axis", y = NULL,
    caption = CAP1) +
  theme_e15 +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.6),
                 legend.key.height = ggplot2::unit(10, "pt"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0),
                 legend.box.spacing = ggplot2::unit(4, "pt"))
.save(g1, "E15_fig1_gene_gap_dumbbell", 8.6, 10.4)

# --- FIG 2: the same 44 as a numbered heatmap -------------------------------
# The dumbbell shows shape; this shows values, before and after adjustment, in
# both cohorts, with the halves separated because by this point the reader has
# already been shown that the separation is not imposed.
COL_LEVELS <- c("MYC TCGA", "MYC SCAN-B", "OXPHOS TCGA", "OXPHOS SCAN-B",
                "OXPHOS - MYC TCGA", "OXPHOS - MYC SCAN-B")
COL_LABELS <- c("MYC\nTCGA", "MYC\nSCAN-B", "OXPHOS\nTCGA", "OXPHOS\nSCAN-B",
                "OX - MYC\nTCGA", "OX - MYC\nSCAN-B")
f2 <- dplyr::bind_rows(
  gene_axes %>% dplyr::transmute(cohort, adjustment, gene, glab, mito_class,
                                 quantity = axis, value = rho),
  pairs_tab %>% dplyr::transmute(cohort, adjustment, gene, glab, mito_class,
                                 quantity = "OXPHOS - MYC", value = gap)) %>%
  dplyr::mutate(
    col = factor(paste(quantity, cohort), levels = COL_LEVELS),
    adjustment = factor(adjustment, levels = ADJ_LEVELS))
stopifnot(!anyNA(f2$col), nrow(f2) == N_CANON * 2L * 2L * 3L)
LIM2 <- max(abs(f2$value))

CAP2 <- sprintf(paste0(
  "The same 44 rows in the same order as figure 1, now with the numbers on",
  " them and with the two\n",
  "halves separated. The right-hand pair of columns is the DIFFERENCE drawn in",
  " figure 1; it shares\n",
  "the colour scale with the two rho pairs, which is defensible because it is",
  " in the same units,\n",
  "but it is a difference of correlations and not a correlation.\n",
  "THE RAW AND ADJUSTED PANELS ARE NOT INTERCHANGEABLE. Trap 3: every MYC",
  " signature is entangled\n",
  "with proliferation, so the raw MYC columns are partly a proliferation",
  " readout and the adjusted\n",
  "ones are the estimate to quote. Proliferation adjustment barely moves the",
  " OXPHOS columns\n",
  "(split %.2f to %.2f in TCGA, %.2f to %.2f in SCAN-B) and rearranges the MYC",
  " ones more than it\n",
  "shrinks them - E11 figure 3 is that rearrangement drawn as crossings.\n",
  "* = below the 25th expression percentile in at least one cohort. [ox] =",
  " in the OXPHOS arm.\n",
  "[p] = in the proliferation covariate. Rows are not independent tests and no",
  " cell of this grid\n",
  "is a finding on its own."),
  .sp("TCGA", "OXPHOS", "raw"), .sp("TCGA", "OXPHOS"),
  .sp("SCAN-B", "OXPHOS", "raw"), .sp("SCAN-B", "OXPHOS"))

g2 <- ggplot2::ggplot(f2, ggplot2::aes(col, glab)) +
  ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white",
                     linewidth = 0.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", value),
                 colour = ifelse(abs(value) > 0.34, "white", "grey15")),
    size = 1.95) +
  ggplot2::facet_grid(mito_class ~ adjustment, scales = "free_y",
                      space = "free_y") +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white",
                                high = "#b2182b", midpoint = 0,
                                limits = c(-LIM2, LIM2),
                                name = "Spearman rho (right-hand pair: a difference of two)") +
  ggplot2::scale_x_discrete(labels = stats::setNames(COL_LABELS, COL_LEVELS)) +
  ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 12,
                                                  barheight = 0.5,
                                                  title.position = "top")) +
  ggplot2::labs(
    title = "The 44 canonical apoptosis genes against both axes, both cohorts",
    subtitle = paste0("EXPLORATORY - not pre-registered | row order fixed to ",
                      "figure 1 | partial Spearman on ", PROLIF_COV),
    x = NULL, y = NULL, caption = CAP2) +
  theme_e15 +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.4),
                 axis.text.x = ggplot2::element_text(size = 6.6),
                 panel.grid = ggplot2::element_blank(),
                 strip.text.y = ggplot2::element_text(size = 7))
.save(g2, "E15_fig2_gene_rho_heatmap", 9.6, 10.4)

# --- FIG 3: how much of the difference is the axes' spread ------------------
f3 <- split_tab %>%
  dplyr::mutate(quantity = factor(quantity, levels = rev(levels(quantity))))
CAP3 <- sprintf(paste0(
  "Each point is Spearman between the quantity and MitoCarta membership across",
  " the 44 genes - the\n",
  "study's anchor statistic, computed on four different quantities. Read it",
  " downwards.\n",
  "MYC separates the regulon half from the cytosolic half weakly (%.2f / %.2f",
  " adjusted). OXPHOS\n",
  "separates them strongly (%.2f / %.2f).\n",
  "The DIFFERENCE between the axes orders them just as strongly as OXPHOS",
  " alone, which is the same\n",
  "statement twice: MYC contributes almost nothing to the ordering, so removing",
  " it changes nothing.\n",
  "THE BOTTOM ROW IS THE CAVEAT MADE ARITHMETIC. Dividing each axis by the SD",
  " of its own 44 values\n",
  "before differencing removes the advantage OXPHOS gets from being the",
  " wider-spread axis. The split\n",
  "falls from %.2f / %.2f to %.2f / %.2f. So SOME of the gap in figures 1 and 2",
  " is the spread\n",
  "difference - and most of it is not, and what remains still sits well above",
  " the MYC row.\n",
  "This check was run while the figures were being designed, so its answer was",
  " known before it was\n",
  "written down. It is a caveat quantified, not a prediction tested, and it is",
  " not a result."),
  .sp("TCGA", "MYC"), .sp("SCAN-B", "MYC"),
  .sp("TCGA", "OXPHOS"), .sp("SCAN-B", "OXPHOS"),
  .sp("TCGA", "OXPHOS - MYC"), .sp("SCAN-B", "OXPHOS - MYC"),
  .sp("TCGA", "OXPHOS - MYC, standardised"),
  .sp("SCAN-B", "OXPHOS - MYC, standardised"))

g3 <- ggplot2::ggplot(f3, ggplot2::aes(value, quantity)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = value, yend = quantity),
                        colour = "grey75", linewidth = 0.5) +
  ggplot2::geom_point(ggplot2::aes(shape = cohort), size = 2.4,
                      colour = "#1f4e79") +
  ggplot2::facet_wrap(~ adjustment) +
  ggplot2::scale_shape_manual(values = c(TCGA = 16, `SCAN-B` = 1), name = NULL) +
  ggplot2::labs(
    title = "How much of the OXPHOS-minus-MYC difference is a difference in spread",
    subtitle = paste("EXPLORATORY - not pre-registered | Spearman with",
                     "MitoCarta membership across the 44 genes"),
    x = "split: Spearman of the quantity with MitoCarta membership", y = NULL,
    caption = CAP3) +
  theme_e15
.save(g3, "E15_fig3_gap_vs_spread", 8.6, 4.6)

# --- FIG 4: the two halves as bars, mean +/- SD -----------------------------
# The author asked for the difference as a bar chart. The two components are
# drawn beside it in the same units and on the same y axis, because a reader
# shown only the right-hand pair cannot tell whether the difference comes from
# OXPHOS moving or from MYC moving, and the answer is the whole point: the MYC
# bars are the same height as each other and the OXPHOS bars are opposite.
#
# THE ERROR BAR IS A SPREAD, NOT AN UNCERTAINTY. It is the SD of the 20 or 24
# genes in the group, so it is nearly as wide as the distance between the two
# means and the bars overlap heavily. That is the honest picture: these groups
# separate ON AVERAGE and not gene by gene. The individual genes are drawn on
# top so the SD can be read as what it is, and the SE of each mean - about a
# fifth of the SD - is in the caption and in `bar_tab`.
f4 <- dplyr::bind_rows(
  gene_axes %>% dplyr::filter(adjustment == ADJ_MAIN) %>%
    dplyr::transmute(cohort, gene, mito_class, quantity = axis, value = rho),
  pairs_tab %>% dplyr::filter(adjustment == ADJ_MAIN) %>%
    dplyr::transmute(cohort, gene, mito_class, quantity = "OXPHOS - MYC",
                     value = gap)) %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_LEVELS),
                quantity = factor(quantity,
                                  levels = c("MYC", "OXPHOS", "OXPHOS - MYC")))
stopifnot(nrow(f4) == N_CANON * 2L * 3L, !anyNA(f4$quantity))

bar_tab <- f4 %>%
  dplyr::group_by(cohort, quantity, mito_class) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(value), sd = stats::sd(value),
                   se = stats::sd(value) / sqrt(dplyr::n()), .groups = "drop") %>%
  dplyr::mutate(lo = mean - sd, hi = mean + sd)
message("\n   figure 4 as a table - mean +/- SD by half, under ", ADJ_MAIN, ":")
bar_tab %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# The difference of the two group means, and the same difference after each
# axis is standardised - the V2 caveat carried into this figure's units.
bar_diff <- dplyr::bind_rows(
  bar_tab %>% dplyr::select(cohort, quantity, mito_class, mean) %>%
    tidyr::pivot_wider(names_from = mito_class, values_from = mean) %>%
    dplyr::mutate(scale = "rho"),
  pairs_tab %>% dplyr::filter(adjustment == ADJ_MAIN) %>%
    dplyr::group_by(cohort, mito_class) %>%
    dplyr::summarise(m = mean(gap_z), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = mito_class, values_from = m) %>%
    dplyr::mutate(quantity = "OXPHOS - MYC", scale = "SD-standardised")) %>%
  dplyr::mutate(diff = .data[[MITO_LEVELS[1]]] - .data[[MITO_LEVELS[2]]]) %>%
  dplyr::select(cohort, quantity, scale, dplyr::all_of(MITO_LEVELS), diff)
message("\n   difference between the two group means:")
bar_diff %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Does the difference survive dropping every gene the figures flag? Three are
# partly self-correlated or self-adjusted, five are low-expression. Eight of 44.
FLAGGED <- sort(unique(c(annot$gene[annot$in_oxphos_arm | annot$in_prolif_cov],
                         LOW_GENES)))
bar_trim <- pairs_tab %>%
  dplyr::filter(adjustment == ADJ_MAIN, !gene %in% FLAGGED) %>%
  dplyr::group_by(cohort, mito_class) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(gap), sd = stats::sd(gap),
                   .groups = "drop")
message("\n   the same difference with the ", length(FLAGGED),
        " flagged genes deleted (", paste(FLAGGED, collapse = ", "), "):")
bar_trim %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

.bd <- function(coh, q, sc = "rho")
  .n(bar_diff, cohort == coh, quantity == q, scale == sc)$diff
.bt <- function(coh, cls, col)
  .n(bar_tab, cohort == coh, quantity == "OXPHOS - MYC",
     mito_class == cls)[[col]]
.tr <- function(coh, cls)
  .n(bar_trim, cohort == coh, mito_class == cls)$mean

CAP4 <- sprintf(paste0(
  "The right-hand pair is what was asked for: the per-gene difference between",
  " the axes, averaged\n",
  "within each half. The two left pairs are its components, on the same y axis",
  " so the heights are\n",
  "comparable.\n",
  "READ ACROSS. The GAP BETWEEN THE TWO BARS is %+.2f / %+.2f on OXPHOS (TCGA",
  " / SCAN-B) and only\n",
  "%+.2f / %+.2f on MYC - MYC separates the halves three to four times less,",
  " and its two bars sit on\n",
  "the same side of zero while the OXPHOS pair straddles it. Means subtract",
  " exactly, so the right-hand\n",
  "pair is the middle pair minus the left one and its separation is %+.2f /",
  " %+.2f.\n",
  "ON THE RANK STATISTIC THE SUBTRACTION COSTS NOTHING (figures 1 and 3: 0.45",
  " to 0.47 in TCGA, 0.49\n",
  "to 0.49 in SCAN-B) WHILE ON THIS MEAN SCALE IT COSTS %+.2f / %+.2f. The two",
  " are not in conflict -\n",
  "a rank statistic and a mean weight genes differently - but neither should be",
  " quoted as the other.\n",
  "Standardising each axis by the SD of its own 44 values first puts the",
  " separation at %+.2f / %+.2f,\n",
  "which is the figure-3 caveat in this figure's units.\n",
  "THE ERROR BAR IS +/- ONE SD OF THE GENES IN THE GROUP, NOT AN UNCERTAINTY",
  " ON THE MEAN. It is\n",
  "almost as wide as the distance between the means, and the bars overlap:",
  " these halves separate ON\n",
  "AVERAGE and not gene by gene. Every gene is drawn on top so that spread can",
  " be seen rather than\n",
  "inferred. The SE of each mean is about a fifth of the SD - for the",
  " difference, %.3f and %.3f\n",
  "(mitochondrial) against %.3f and %.3f (cytosolic).\n",
  "n = 20 mitochondrial and 24 cytosolic. Deleting all %d flagged genes - the",
  " three partly\n",
  "self-correlated or self-adjusted and the five below the 25th expression",
  " percentile - leaves the\n",
  "means at %+.3f / %+.3f (mitochondrial) and %+.3f / %+.3f (cytosolic), so",
  " the difference is not\n",
  "carried by them. NO TEST IS REPORTED. Two groups of a curated 44 in a study",
  " whose atlas is a grid\n",
  "of thousands of cells; the reading is the direction, its size relative to",
  " the spread, and that\n",
  "both cohorts show it."),
  .bd("TCGA", "OXPHOS"), .bd("SCAN-B", "OXPHOS"),
  .bd("TCGA", "MYC"), .bd("SCAN-B", "MYC"),
  .bd("TCGA", "OXPHOS - MYC"), .bd("SCAN-B", "OXPHOS - MYC"),
  .bd("TCGA", "MYC"), .bd("SCAN-B", "MYC"),
  .bd("TCGA", "OXPHOS - MYC", "SD-standardised"),
  .bd("SCAN-B", "OXPHOS - MYC", "SD-standardised"),
  .bt("TCGA", MITO_CLS, "se"), .bt("SCAN-B", MITO_CLS, "se"),
  .bt("TCGA", CYTO_CLS, "se"), .bt("SCAN-B", CYTO_CLS, "se"),
  length(FLAGGED),
  .tr("TCGA", MITO_CLS), .tr("SCAN-B", MITO_CLS),
  .tr("TCGA", CYTO_CLS), .tr("SCAN-B", CYTO_CLS))

g4 <- ggplot2::ggplot(bar_tab, ggplot2::aes(mito_class, mean)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_col(ggplot2::aes(fill = mito_class), width = 0.62,
                    alpha = 0.4, colour = NA) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi,
                                      colour = mito_class),
                         width = 0.16, linewidth = 0.55) +
  ggplot2::geom_point(data = f4, ggplot2::aes(mito_class, value,
                                              colour = mito_class),
                      position = ggplot2::position_jitter(width = 0.17,
                                                          height = 0,
                                                          seed = PROJECT_SEED),
                      size = 0.85, alpha = 0.65, show.legend = FALSE) +
  ggplot2::geom_point(ggplot2::aes(colour = mito_class), size = 2.3) +
  ggplot2::facet_grid(cohort ~ quantity) +
  ggplot2::scale_fill_manual(values = COMP_COLS, name = NULL) +
  ggplot2::scale_colour_manual(values = COMP_COLS, name = NULL) +
  ggplot2::scale_x_discrete(labels = c("mitochondrial\n(MitoCarta, n = 20)",
                                       "cytosolic\n(not in MitoCarta, n = 24)")) +
  ggplot2::labs(
    title = paste("OXPHOS separates the two halves in opposite directions;",
                  "MYC barely separates them"),
    subtitle = paste0("EXPLORATORY - not pre-registered | mean +/- 1 SD of the ",
                      "per-gene values | partial Spearman on ", PROLIF_COV),
    x = NULL, y = "per-gene Spearman rho (right column: a difference of two)",
    caption = CAP4) +
  theme_e15 +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6.8),
                 legend.position = "none")
.save(g4, "E15_fig4_gap_by_membership_bars", 8.4, 7.4)

# --- FIG 5: the difference alone, gene by gene, one panel per cohort --------
# Figure 1 with the two component points taken away. What is left is the bar
# length, which is the quantity the argument is about, and dropping the points
# buys enough room to put the gene names on the x axis and give each cohort its
# own panel.
#
# BOTH PANELS USE THE SAME GENE ORDER - the cross-cohort mean, as in figures 1
# and 2 - rather than each sorting on its own values. Per-panel sorting makes a
# tidier monotone descent and makes the two panels incomparable, which is the
# wrong trade in a study whose unit of evidence is cross-cohort agreement. The
# cost is that neither panel is monotone; that departure from a clean slope IS
# the replication being displayed.
f5 <- pairs_tab %>%
  dplyr::filter(adjustment == ADJ_MAIN) %>%
  dplyr::left_join(dplyr::select(gene_order, gene, agree), by = "gene") %>%
  dplyr::mutate(cohort = factor(cohort, levels = COHORT_LEVELS),
                glab = factor(as.character(glab), levels = rev(LAB_LEVELS)))
stopifnot(nrow(f5) == N_CANON * 2L, !anyNA(f5$glab))

CAP5 <- sprintf(paste0(
  "The bar length from figure 1, on its own, with the genes on the x axis and a",
  " panel per cohort.\n",
  "Above zero the gene tracks OXPHOS more than MYC; below zero, MYC more than",
  " OXPHOS.\n",
  "GENE ORDER IS THE CROSS-COHORT MEAN AND IS SHARED BY BOTH PANELS, so the",
  " two can be read against\n",
  "each other. Neither panel is therefore monotone, and where it breaks is",
  " where the cohorts differ.\n",
  "Mean +/- SD: mitochondrial %+.2f +/- %.2f and %+.2f +/- %.2f, cytosolic",
  " %+.2f +/- %.2f and\n",
  "%+.2f +/- %.2f (TCGA and SCAN-B). The SD is a spread across genes, not an",
  " uncertainty on the mean -\n",
  "figure 4 draws both together.\n",
  "A HOLLOW POINT is one of the %d genes whose two cohorts disagree on the",
  " sign: %s.\n",
  "* = below the 25th expression percentile in at least one cohort. [ox] =",
  " CYCS, in the OXPHOS arm.\n",
  "[p] = in the proliferation covariate. THE DIFFERENCE IS NOT IN A COMMON",
  " UNIT - the 44 rho values\n",
  "spread wider on OXPHOS than on MYC, and figure 3 puts a number on how much",
  " of the split that buys."),
  .bt("TCGA", MITO_CLS, "mean"), .bt("TCGA", MITO_CLS, "sd"),
  .bt("SCAN-B", MITO_CLS, "mean"), .bt("SCAN-B", MITO_CLS, "sd"),
  .bt("TCGA", CYTO_CLS, "mean"), .bt("TCGA", CYTO_CLS, "sd"),
  .bt("SCAN-B", CYTO_CLS, "mean"), .bt("SCAN-B", CYTO_CLS, "sd"),
  length(DISAGREE), paste(DISAGREE, collapse = ", "))

g5 <- ggplot2::ggplot(f5, ggplot2::aes(glab, gap)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_segment(ggplot2::aes(xend = glab, y = 0, yend = gap,
                                     colour = mito_class), linewidth = 0.75) +
  ggplot2::geom_point(ggplot2::aes(colour = mito_class, shape = agree),
                      size = 1.9) +
  ggplot2::facet_grid(cohort ~ .) +
  ggplot2::scale_colour_manual(values = COMP_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                              breaks = c("FALSE"),
                              labels = c("cohorts disagree on the sign"),
                              name = NULL) +
  ggplot2::labs(
    title = "The OXPHOS-minus-MYC difference, gene by gene, in each cohort",
    subtitle = paste0("EXPLORATORY - not pre-registered | 44 canonical genes | ",
                      "partial Spearman on ", PROLIF_COV,
                      " | gene order shared by both panels"),
    x = NULL,
    y = "per-gene rho with OXPHOS minus per-gene rho with MYC",
    caption = CAP5) +
  theme_e15 +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                     vjust = 0.5, size = 6.2),
                 panel.grid.major.x = ggplot2::element_blank(),
                 legend.key.height = ggplot2::unit(10, "pt"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))
.save(g5, "E15_fig5_difference_cleveland", 10.2, 6.6)

# =============================================================================
# 5. Save
# =============================================================================
message("\n5. save")
csv_out <- pairs_tab %>%
  dplyr::select(cohort, adjustment, gene, mitocarta, compartment, effect,
                cdc_module, low_expression, in_oxphos_arm, in_prolif_cov,
                rho_MYC, rho_OXPHOS, gap, rho_z_MYC, rho_z_OXPHOS, gap_z) %>%
  dplyr::arrange(adjustment, cohort, dplyr::desc(gap))
saveRDS(list(
  gene_axes = gene_axes, pairs_tab = pairs_tab, gene_order = gene_order,
  annot = annot, sd_tab = sd_tab, split_tab = split_tab, half_tab = half_tab,
  rep_tab = rep_tab, ends = ends, anchor_chk = anchor_chk,
  bar_tab = bar_tab, bar_diff = bar_diff, bar_trim = bar_trim,
  flagged = FLAGGED,
  low_expr = low_expr, disagree = DISAGREE, lab_levels = LAB_LEVELS,
  settings = list(myc_axis = MYC_REF, prolif_covariate = PROLIF_COV,
                  adjustments = ADJ_LEVELS, main_adjustment = ADJ_MAIN,
                  low_expr_pct = LOW_EXPR_PCT, n_canon = N_CANON,
                  reads = c(PATH_E08, PATH_E11, PATH_E14)),
  rules = list(
    computes_nothing = paste("this is a display script. Every rho is E11's,",
                             "read from gene_tab; section 1 asserts the",
                             "localisation split recomputed from them is",
                             "bit-equal to E14's published value for the same",
                             "44 genes in all 8 cells."),
    gap_units = paste("`gap` is rho_OXPHOS - rho_MYC and is NOT a difference",
                      "in a common unit: the 44 correlations spread wider on",
                      "OXPHOS than on MYC, so a gene earns part of its gap by",
                      "sitting on the wider axis. `gap_z` divides each axis by",
                      "the SD of its own 44 values first and is the version",
                      "immune to that. Both are saved; neither is reportable",
                      "without the other."),
    honesty = paste("the gap_z check was run while the figures were being",
                    "designed, so its answer was known before it was written",
                    "down. It is a caveat quantified, not a prediction tested."),
    row_order = paste("row order is the cross-cohort mean gap under",
                      "PROLIF_DISJOINT and is fixed across all three figures.",
                      "It is data-driven and does not know about MitoCarta,",
                      "which is why the mitochondrial rows collecting at the",
                      "top of figure 1 is a property of the numbers."),
    no_single_gene = paste("44 genes x 2 cohorts x 2 adjustments is a grid.",
                           "The 44 were not selected on the statistic shown,",
                           "which is better than E08's driver lists, but no",
                           "cell is a finding. Five genes sit below the 25th",
                           "expression percentile in at least one cohort and",
                           "are starred on every panel."),
    error_bar = paste("figure 4's error bar is +/- one SD of the 20 or 24 genes",
                      "in the group, NOT an uncertainty on the mean. It is",
                      "almost as wide as the distance between the two means:",
                      "the halves separate on average and not gene by gene.",
                      "Every gene is drawn on top for that reason, and the SE",
                      "is in `bar_tab`. No test is reported."),
    control = paste("none of these figures may be read as MYC being unrelated",
                    "to the machinery. E11 figure 2 is the control that shows",
                    "the adjusted MYC axis still tracks the mitoribosome; an",
                    "empty MYC column and an emptied MYC score look the same",
                    "here.")),
  built = Sys.time()), PATH_E15)
readr::write_csv(csv_out, PATH_E15_CSV)
message("\nE15: done.")
message("    results/two_axis_gene_view.rds")
message("    outputs/tables/E15_gene_rho_two_axes.csv")
message("    5 figures in outputs/figures/:")
message("      fig1 the dumbbell - one bar per gene, sorted by the difference")
message("      fig2 the same 44 as a numbered heatmap, raw and adjusted")
message("      fig3 how much of the difference is the axes' spread")
message("      fig4 the two halves as bars, mean +/- SD, with every gene on top")
message("      fig5 the difference alone, gene by gene, one panel per cohort")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E15)

  # The anchor: this must be bit-equal to E14's machinery split.
  x$anchor_chk %>% as.data.frame()

  # The four splits side by side - figure 3 as a table.
  x$split_tab %>%
    tidyr::pivot_wider(names_from = cohort, values_from = value) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # The halves.
  x$half_tab %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # Does the gap replicate?
  x$rep_tab %>% as.data.frame()
  x$disagree

  # The ordering, end to end.
  x$gene_order %>%
    dplyr::arrange(dplyr::desc(mean_gap)) %>%
    dplyr::transmute(gene, compartment, effect, agree,
                     MYC = round(mean_MYC, 3), OXPHOS = round(mean_OXPHOS, 3),
                     gap = round(mean_gap, 3), gap_z = round(mean_gap_z, 3)) %>%
    as.data.frame()

  # The mitochondrial half by sub-compartment: is there a depth gradient?
  x$pairs_tab %>%
    dplyr::filter(adjustment == "adj. PROLIF_DISJOINT", mitocarta) %>%
    dplyr::group_by(cohort, compartment) %>%
    dplyr::summarise(n = dplyr::n(),
                     med_OXPHOS = round(stats::median(rho_OXPHOS), 3),
                     med_MYC = round(stats::median(rho_MYC), 3),
                     med_gap = round(stats::median(gap), 3), .groups = "drop") %>%
    as.data.frame()

  # Figure 4 as a table: mean +/- SD, and the difference between the means.
  x$bar_tab %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$bar_diff %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$bar_trim %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$flagged

  # The starred genes, which must never be read as results.
  x$low_expr %>% dplyr::filter(low_expression) %>% as.data.frame()
}
