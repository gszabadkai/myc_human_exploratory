# E04_figures_correlation.R
# =============================================================================
# The figures for the MYC-OXPHOS plane. Reads E03's atlas and E02's scores;
# computes nothing new except what a plot needs to draw a line.
#
# Built to: docs/2026-08-31_phase1_plan.md section 2 (E04)
# Findings these figures carry: docs/2026-08-31_phase1_atlas_findings.md
#
# =============================================================================
# EVERY FIGURE IS LABELLED EXPLORATORY, ON THE FIGURE
# =============================================================================
# CLAUDE.md: findings are hypothesis-generating and must be labelled as such
# "in notes, in figures and in conversation". Not in a caption written later -
# on the panel, in the subtitle, every time. .lab_exploratory() does it and no
# figure here is built without it.
#
# TWO OTHER THINGS GO ON THE FIGURE RATHER THAN IN A METHODS PARAGRAPH:
#   - SCAN-B has no purity estimate (trap 2). Any panel showing an adjustment
#     says which cohort could be adjusted and which could not.
#   - mitoPPS is composition, not level (trap 6). Its panels say so, and no
#     figure puts a mitoPPS VALUE from one cohort beside one from the other.
#
# NO P-VALUE IS PLOTTED ANYWHERE. The atlas has 68,255 cells; a p-value on a
# figure invites the reader to pick one. Intervals are drawn instead.
#
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE04: figures for the MYC-OXPHOS plane\n", strrep("=", 78))

# =============================================================================
# 0. Constants and shared look
# =============================================================================
FIG_DPI <- 300
COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
INSTRUMENT_ORDER <- c("gsva", "zmean", "content", "mitopps")

.lab_exploratory <- function(extra = NULL) {
  paste(c("EXPLORATORY - not pre-registered; hypothesis-generating only",
          extra), collapse = " | ")
}

theme_atlas <- function(base = 10) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = base - 1),
      plot.title       = ggplot2::element_text(face = "bold", size = base + 2),
      plot.subtitle    = ggplot2::element_text(size = base - 1, colour = "grey30"),
      plot.caption     = ggplot2::element_text(size = base - 2, colour = "grey40",
                                               hjust = 0),
      legend.position  = "bottom")
}

.save_fig <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0(name, ".", ext)), p,
                    width = w, height = h, dpi = FIG_DPI, limitsize = FALSE)
  }
  message("   ", name, "  (", w, " x ", h, " in)")
  invisible(p)
}

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

a      <- readRDS(file.path(DIR_RESULTS, "correlation_atlas.rds"))
frames <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
sc     <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw     <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito   <- readRDS(PATH_TCGA_MITO)

A <- a$atlas %>%
  dplyr::mutate(instrument = factor(instrument,
                                    levels = c(INSTRUMENT_ORDER, "gene")),
                cohort = factor(cohort, levels = names(COHORT_COLS)))
ARM_ORDER <- A %>%
  dplyr::filter(cohort == "TCGA", instrument == "gsva",
                myc_estimator == "FELSHER_61", stratum == "all",
                adjusted == "raw", measure_class == "arm") %>%
  dplyr::arrange(dplyr::desc(rho)) %>% dplyr::pull(arm)
message("   atlas ", format(nrow(A), big.mark = ","), " cells | ",
        length(ARM_ORDER), " arms")

# =============================================================================
# 2. Per-sample scores, for the scatters only
# =============================================================================
# GSVA IS COHORT-RELATIVE. The two cohorts are drawn in separate panels with
# free scales and are never pooled into one cloud. What is comparable is the
# SLOPE and the spread, not the position.
message("\n2. per-sample scores for the scatters")

.plane <- function(coh, myc, arms_obj, ids) {
  tibble::tibble(cohort = coh, sample_id = ids,
                 myc = as.numeric(myc[ids]),
                 gsva    = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]),
                 zmean   = as.numeric(arms_obj$zmean_arms["OXPHOS subunits", ids]),
                 content = as.numeric(arms_obj$content_arms["OXPHOS subunits", ids]),
                 mitopps = as.numeric(arms_obj$mitopps_arms["OXPHOS subunits", ids]))
}
ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
plane <- dplyr::bind_rows(
  .plane("TCGA",   nw$tcga_gsva_new["FELSHER_61", ], mito, ID_T),
  .plane("SCAN-B", sc$gsva_new["FELSHER_61", ],      sc,   ID_S)) %>%
  dplyr::left_join(frames %>%
                     dplyr::select(cohort, sample_id, PAM50, ER) %>%
                     dplyr::mutate(cohort = as.character(cohort)),
                   by = c("cohort", "sample_id")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
plane_long <- plane %>%
  tidyr::pivot_longer(dplyr::all_of(INSTRUMENT_ORDER),
                      names_to = "instrument", values_to = "oxphos") %>%
  dplyr::mutate(instrument = factor(instrument, levels = INSTRUMENT_ORDER))
message("   ", nrow(plane), " samples on the plane")

# =============================================================================
# 3. Figure 1 - the plane itself
# =============================================================================
message("\n3. figures")

f1 <- ggplot2::ggplot(plane_long,
                      ggplot2::aes(myc, oxphos, colour = cohort)) +
  ggplot2::geom_point(size = 0.35, alpha = 0.25) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                       colour = "black", linewidth = 0.5) +
  ggplot2::facet_grid(instrument ~ cohort, scales = "free") +
  ggplot2::scale_colour_manual(values = COHORT_COLS, guide = "none") +
  ggplot2::labs(
    title = "MYC activity against OXPHOS subunits, four instruments",
    subtitle = .lab_exploratory("MYC = FELSHER_61 (GSVA); all samples, unadjusted"),
    x = "MYC activity (FELSHER_61, GSVA)", y = "OXPHOS subunits",
    caption = paste("Scores are COHORT-RELATIVE: panels have free scales and",
                    "the two cohorts are never pooled. mitoPPS is composition,",
                    "not level -\nits values are not comparable between cohorts,",
                    "only its pattern is. Line is OLS, shown to display",
                    "direction; the atlas reports Spearman.")) +
  theme_atlas()
.save_fig(f1, "E04_fig1_plane_by_instrument", 7.5, 9)

f1b <- ggplot2::ggplot(dplyr::filter(plane, !is.na(PAM50)),
                       ggplot2::aes(myc, gsva, colour = ER)) +
  ggplot2::geom_point(size = 0.4, alpha = 0.35) +
  ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm", formula = y ~ x,
                       se = TRUE, colour = "black", linewidth = 0.5) +
  ggplot2::facet_grid(cohort ~ PAM50, scales = "free") +
  ggplot2::scale_colour_manual(values = c(ERpos = "#7b3294", ERneg = "#008837"),
                               na.value = "grey70") +
  ggplot2::labs(
    title = "The same plane, split by PAM50 subtype and ER call",
    subtitle = .lab_exploratory("GSVA only; the correlation is present in every subtype"),
    x = "MYC activity (FELSHER_61, GSVA)", y = "OXPHOS subunits (GSVA)",
    caption = paste("TCGA HER2 (n = 78) and Normal (n = 36) are small; see the",
                    "intervals in figure 6 before reading them.")) +
  theme_atlas()
.save_fig(f1b, "E04_fig1b_plane_by_subtype", 11, 5.5)

# =============================================================================
# 4. Figure 2 - the arm x instrument heatmap
# =============================================================================
# Where "the mitoribosome beats OXPHOS" becomes visible, and where the mtDNA arm
# separates from every other arm on mitoPPS alone.
heat <- A %>%
  dplyr::filter(measure_class == "arm", stratum == "all", adjusted == "raw",
                myc_estimator == "FELSHER_61") %>%
  dplyr::mutate(arm = factor(arm, levels = rev(ARM_ORDER)))

f2 <- ggplot2::ggplot(heat, ggplot2::aes(instrument, arm, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", rho)), size = 2.6) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                midpoint = 0, limits = c(-0.7, 0.7),
                                name = "Spearman rho") +
  ggplot2::labs(
    title = "MYC activity against 18 mitochondrial arms, four instruments",
    subtitle = .lab_exploratory("FELSHER_61, all samples, unadjusted; arms ordered by TCGA GSVA"),
    x = NULL, y = NULL,
    caption = paste("The mitochondrial ribosome outranks OXPHOS subunits in",
                    "BOTH cohorts. mtDNA-encoded OXPHOS is the only arm whose",
                    "sign\ndepends on the instrument - negative on mitoPPS",
                    "(composition), ~0 on the other three (level).")) +
  theme_atlas()
.save_fig(f2, "E04_fig2_arm_instrument_heatmap", 9, 6.5)

# =============================================================================
# 5. Figure 3 - THE estimator panel, ordered by entanglement
# =============================================================================
# CLAUDE.md trap 3, and the sharpest test in phase 1. Read left to right:
# proliferation content rises. If rho rose with it, the correlation would be
# proliferation wearing MYC's name.
panel <- a$panel %>%
  dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                adjusted %in% c("raw", "prolif")) %>%
  dplyr::mutate(
    cohort = factor(cohort, levels = names(COHORT_COLS)),
    adjusted = factor(adjusted, levels = c("raw", "prolif"),
                      labels = c("unadjusted", "proliferation-adjusted")),
    myc_estimator = factor(myc_estimator,
      levels = a$est_meta %>% dplyr::filter(!is.na(frac_prolif)) %>%
        dplyr::arrange(frac_prolif) %>% dplyr::pull(myc_estimator)))

f3 <- ggplot2::ggplot(panel, ggplot2::aes(myc_estimator, rho, colour = cohort)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 0.2, linetype = 3, linewidth = 0.3,
                      colour = "grey50") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                           position = ggplot2::position_dodge(width = 0.55),
                           size = 0.28, linewidth = 0.45) +
  ggplot2::facet_wrap(~ adjusted, ncol = 1) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "Does the MYC-OXPHOS correlation track proliferation entanglement?",
    subtitle = .lab_exploratory(
      "18 signatures, ordered LEFT to RIGHT by % of genes in HALLMARK E2F + G2M (1.5% to 47.6%)"),
    x = NULL, y = "Spearman rho with OXPHOS subunits (GSVA)",
    caption = paste(
      "The least entangled signature (MYC_UP.V1_UP, 1.5%) gives one of the",
      "HIGHEST correlations and keeps 0.40 / 0.52 after adjustment.",
      "\n15 of 18 (SCAN-B) and 17 of 18 (TCGA) stay above 0.2 adjusted.",
      "Between-signature spread is far larger than any entanglement trend, so",
      "\nno single signature may be quoted - CLAUDE.md trap 3. Bars are 95%",
      "Fisher-z intervals; no p-value is shown or implied.")) +
  theme_atlas() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                     size = 7))
.save_fig(f3, "E04_fig3_estimator_panel_entanglement", 10, 7.5)

# =============================================================================
# 6. Figure 4 - the two genomes
# =============================================================================
# Sub-analysis (i). Top: the arms. Bottom: the 13 mtDNA genes individually,
# because they do not move together and an arm score averages that away.
nm <- a$nuclear_vs_mtdna %>%
  dplyr::filter(myc_estimator == "FELSHER_61", adjusted == "raw") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                instrument = factor(instrument, levels = INSTRUMENT_ORDER),
                arm = factor(arm, levels = rev(c(
                  "OXPHOS subunits", "CI subunits", "CII subunits",
                  "CIII subunits", "CIV subunits", "CV subunits",
                  "mtDNA-encoded OXPHOS"))))

f4a <- ggplot2::ggplot(nm, ggplot2::aes(rho, arm, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                           position = ggplot2::position_dodge(width = 0.5),
                           size = 0.25, linewidth = 0.45) +
  ggplot2::facet_wrap(~ instrument, nrow = 1) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "Nuclear-encoded and mtDNA-encoded OXPHOS are not regulated alike",
    subtitle = .lab_exploratory("FELSHER_61, all samples, unadjusted"),
    x = "Spearman rho", y = NULL,
    caption = paste("On mitoPPS - which reports the SHAPE of the mitochondrial",
                    "programme rather than its level - the mtDNA arm is the",
                    "only one that\nturns negative, in both cohorts",
                    "(-0.081 and -0.088). Whether that is biology or an",
                    "arithmetic consequence of the nuclear arm rising is\nNOT",
                    "yet excluded: see falsifier 3 in the findings note.")) +
  theme_atlas()
.save_fig(f4a, "E04_fig4a_nuclear_vs_mtdna_arms", 11, 4.5)

mtg <- a$mtdna_genes %>%
  dplyr::filter(myc_estimator == "FELSHER_61", adjusted == "raw") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
gene_order <- mtg %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::arrange(rho) %>% dplyr::pull(gene)
mtg$gene <- factor(mtg$gene, levels = gene_order)

f4b <- ggplot2::ggplot(mtg, ggplot2::aes(rho, gene, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                           position = ggplot2::position_dodge(width = 0.5),
                           size = 0.3, linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "The 13 mtDNA-encoded genes do not move together",
    subtitle = .lab_exploratory("FELSHER_61 against linear expression, all samples, unadjusted"),
    x = "Spearman rho", y = NULL,
    caption = paste("MT-CO2 (+0.240 / +0.243) and MT-ND5 (-0.160 / -0.143)",
                    "agree between cohorts to 0.003 and 0.017 and sit on the",
                    "SAME heavy-strand\npolycistron, so an opposite-signed",
                    "correlation cannot be mtDNA copy number or",
                    "transcription rate. In TCGA the split survives\npurity and",
                    "leukocyte adjustment. Ordered by TCGA rho.")) +
  theme_atlas()
.save_fig(f4b, "E04_fig4b_mtdna_genes", 7.5, 5.5)

# =============================================================================
# 7. Figure 5 - what the adjustments cost
# =============================================================================
# The plan asks for this explicitly, and it is the only honest way to show a
# purity adjustment that ONE cohort could have.
adj <- A %>%
  dplyr::filter(measure_class == "arm", stratum == "all",
                kind %in% c("signature (GSVA)", "CollecTRI regulon (ULM)",
                            "MYC mRNA")) %>%
  dplyr::select(cohort, instrument, myc_estimator, arm, adjusted, rho) %>%
  tidyr::pivot_wider(names_from = adjusted, values_from = rho) %>%
  tidyr::pivot_longer(dplyr::any_of(c("prolif", "purity_leuko")),
                      names_to = "adjustment", values_to = "rho_adj") %>%
  dplyr::filter(!is.na(rho_adj)) %>%
  dplyr::mutate(adjustment = factor(adjustment,
    levels = c("prolif", "purity_leuko"),
    labels = c("proliferation (PROLIF_DISJOINT)", "purity + leukocyte")))

f5 <- ggplot2::ggplot(adj, ggplot2::aes(raw, rho_adj, colour = instrument)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2) +
  ggplot2::geom_point(size = 0.5, alpha = 0.45) +
  ggplot2::facet_grid(cohort ~ adjustment) +
  ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
  ggplot2::labs(
    title = "What the adjustments cost",
    subtitle = .lab_exploratory("every arm x MYC estimator cell, all samples"),
    x = "unadjusted Spearman rho", y = "adjusted Spearman rho",
    caption = paste("SCAN-B HAS NO PURITY ESTIMATE and none is imputed, so the",
                    "right-hand column exists for TCGA only (CLAUDE.md",
                    "trap 2).\nPoints below the diagonal lost signal to the",
                    "adjustment; the mass staying near it is the point.")) +
  theme_atlas() +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2)))
.save_fig(f5, "E04_fig5_raw_vs_adjusted", 8, 6)

# =============================================================================
# 8. Figure 6 - cross-cohort reproducibility, and the strata
# =============================================================================
rep_arm <- a$reproducibility %>%
  dplyr::filter(measure_class == "arm", stratum == "all", adjusted == "raw") %>%
  dplyr::mutate(instrument = factor(instrument, levels = INSTRUMENT_ORDER))

f6a <- ggplot2::ggplot(rep_arm, ggplot2::aes(rho_tcga, rho_scanb,
                                             colour = instrument)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2) +
  ggplot2::geom_point(size = 0.7, alpha = 0.5) +
  ggplot2::facet_wrap(~ instrument, nrow = 1) +
  ggplot2::scale_colour_brewer(palette = "Dark2", guide = "none") +
  ggplot2::labs(
    title = "Does TCGA's pattern reappear in SCAN-B?",
    subtitle = .lab_exploratory(
      "18 arms x 20 shared MYC estimators, all samples, unadjusted"),
    x = "Spearman rho, TCGA", y = "Spearman rho, SCAN-B",
    caption = paste("This - not any single cell - is what the atlas is for.",
                    "Across the whole comparable grid the two cohorts agree in",
                    "sign in 85.7%\nof cells and their intervals overlap in",
                    "80.9%. Values are not pooled anywhere; only the",
                    "correlations are compared.")) +
  theme_atlas()
.save_fig(f6a, "E04_fig6a_cross_cohort", 11, 3.8)

STRATUM_ORDER <- c("all", "ERpos", "ERneg", "LumA", "LumB", "HER2", "Basal",
                   "Normal")
str_df <- A %>%
  dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                adjusted == "raw", myc_estimator == "FELSHER_61") %>%
  dplyr::mutate(stratum = factor(stratum, levels = rev(STRATUM_ORDER)))

f6b <- ggplot2::ggplot(str_df, ggplot2::aes(rho, stratum, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                           position = ggplot2::position_dodge(width = 0.5),
                           size = 0.3, linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n), x = -0.28),
                     position = ggplot2::position_dodge(width = 0.5),
                     size = 2.4, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "The correlation is present in every subtype",
    subtitle = .lab_exploratory("FELSHER_61 vs OXPHOS subunits, GSVA, unadjusted"),
    x = "Spearman rho", y = NULL,
    caption = paste("ER-negative is the weakest stratum in both cohorts.",
                    "TCGA HER2 (n = 78) and Normal (n = 36) carry intervals",
                    "wide enough that\nthey should not be read as differences.")) +
  theme_atlas()
.save_fig(f6b, "E04_fig6b_strata", 7, 5)

# =============================================================================
# 9. Figure 7 - instrument agreement (trap 5)
# =============================================================================
ia <- a$instrument_agreement %>%
  dplyr::mutate(pair = paste(instrument_a, "vs", instrument_b),
                cohort = factor(cohort, levels = names(COHORT_COLS)))
pair_order <- ia %>% dplyr::group_by(pair) %>%
  dplyr::summarise(m = stats::median(rho_across_arms), .groups = "drop") %>%
  dplyr::arrange(m) %>% dplyr::pull(pair)
ia$pair <- factor(ia$pair, levels = pair_order)

f7 <- ggplot2::ggplot(ia, ggplot2::aes(rho_across_arms, pair, colour = cohort)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.15,
                                                          seed = PROJECT_SEED),
                      size = 1.1, alpha = 0.7) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "The four instruments are not four independent measurements",
    subtitle = .lab_exploratory(
      "each point is one MYC estimator; agreement = rho of the two instruments' 18-arm answers"),
    x = "Spearman rho between the two instruments' arm rankings", y = NULL,
    caption = paste("gsva and zmean agree at a median of 0.994 - they are one",
                    "instrument, not two. Every other pair sits near 0.88 with",
                    "minima\naround 0.50, so disagreement is",
                    "estimator-specific. CLAUDE.md trap 5: report all four, or",
                    "justify the one.")) +
  theme_atlas()
.save_fig(f7, "E04_fig7_instrument_agreement", 8, 4)

# =============================================================================
# 10. Inventory
# =============================================================================
figures <- tibble::tribble(
  ~file,                                   ~shows,
  "E04_fig1_plane_by_instrument",          "the plane itself, 4 instruments x 2 cohorts",
  "E04_fig1b_plane_by_subtype",            "the plane by PAM50 and ER",
  "E04_fig2_arm_instrument_heatmap",       "18 arms x 4 instruments x 2 cohorts, rho as fill",
  "E04_fig3_estimator_panel_entanglement", "trap 3: rho against proliferation entanglement",
  "E04_fig4a_nuclear_vs_mtdna_arms",       "sub-analysis (i), arm level",
  "E04_fig4b_mtdna_genes",                 "sub-analysis (i), the 13 genes individually",
  "E04_fig5_raw_vs_adjusted",              "what proliferation and purity cost",
  "E04_fig6a_cross_cohort",                "TCGA vs SCAN-B over the comparable grid",
  "E04_fig6b_strata",                      "the correlation by ER and PAM50",
  "E04_fig7_instrument_agreement",         "trap 5, agreement of the answers")
saveRDS(list(figures = figures, arm_order = ARM_ORDER,
             cohort_cols = COHORT_COLS, built = Sys.time()),
        file.path(DIR_RESULTS, "figure_manifest.rds"))

message("\nE04: done. ", nrow(figures), " figures (png + pdf) in outputs/figures/")
message("    Every panel is labelled EXPLORATORY and no p-value is drawn.")
message("    NEXT: E05, the cell-death axis.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  print(f3)    # the sharpest test in phase 1
  print(f4b)   # the finding most likely to survive contact with phase 2
  print(f2)

  # the numbers behind figure 3, in case a panel looks wrong
  a$panel %>%
    dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva") %>%
    dplyr::select(cohort, adjusted, myc_estimator, frac_prolif, rho) %>%
    tidyr::pivot_wider(names_from = c(cohort, adjusted), values_from = rho) %>%
    dplyr::arrange(frac_prolif) %>% as.data.frame()

  # and behind figure 4b
  a$mtdna_genes %>%
    dplyr::filter(myc_estimator == "FELSHER_61") %>%
    tidyr::pivot_wider(id_cols = gene, names_from = c(cohort, adjusted),
                       values_from = rho) %>% as.data.frame()

  list.files(DIR_FIGURES, pattern = "^E04_")

}
