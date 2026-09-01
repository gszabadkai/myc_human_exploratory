# E08_strata_and_death_genes.R
# =============================================================================
# PHASE 2. Three questions that all turn out to be about STRATIFICATION and
# GENE-LEVEL resolution, so they share one script and one set of machinery.
#
#   Q-a (handoff candidate 4) THE ER-NEGATIVE GAP. ERneg is the weakest stratum
#       for MYC-OXPHOS in both cohorts - panel median 0.315 / 0.416 against
#       0.479 / 0.512 for ERpos - and nothing has been done with it. Is it
#       ER-negative biology, is it Basal, or is it a range restriction?
#
#   Q-b (author) WHAT DRIVES THE PRO-DEATH / PRO-SURVIVAL ANTICORRELATION.
#       D1's two strata are 502 and 584 genes and sit at roughly z = -2 and +2
#       against an expression-matched null. Is the contrast carried by a
#       minority of genes or spread across the whole curation, and does it track
#       any annotation the curation already carries?
#
#   Q-c (author) THE OVERLAY GENES BY SUBTYPE, AND WHETHER THEY COUNT. E05 put
#       the individual death genes on the plane using the WHOLE cohort. The
#       mouse arm shows luminal expansion, so the luminal compartment - LumA,
#       LumB, and LumA+LumB together - deserves its own reading, beside Basal
#       and the ER split. And a gene's rho means nothing if the gene is barely
#       expressed, so every gene is reported with its expression percentile.
#
# =============================================================================
# THE LUMINAL STRATUM IS NEW, AND IT IS NOT JUST LumA PLUS LumB
# =============================================================================
# functions/strata.R adds Luminal = LumA + LumB (696 TCGA, 2,436 SCAN-B). At
# gene level LumA alone (499 / 1,540) and LumB alone (197 / 896) are thin - a
# 197-sample stratum gives a 95% interval about +/- 0.14 wide on rho - and the
# combined compartment is the one the mouse expansion is about. All three are
# reported, so a difference BETWEEN LumA and LumB stays visible rather than
# being averaged away.
#
# =============================================================================
# THE RULE THAT GOVERNS Q-c, AND IT IS THE POINT OF THE EXPRESSION RANK
# =============================================================================
# A Spearman correlation of a gene that sits at the bottom of the expression
# range is a correlation of quantisation noise with a score. It will still have
# a small confidence interval at n = 3,207 and it will still look like a result.
# Every overlay gene therefore carries `expr_pct`, its percentile among all
# genes of that cohort's matrix, and anything below the 25th percentile is
# flagged. THE FLAG IS NOT A FILTER - low-expression genes are shown, marked,
# and not deleted, because a reader deciding what to believe needs to see them.
#
# SCALE: linear DESeq2-normalised for gene level; the scores as E02 built them.
# All correlations are rank-based. SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))
source(here::here("functions", "strata.R"))

message("\nE08: strata, and the death axis at gene resolution\n", strrep("=", 78))

PATH_E08     <- file.path(DIR_RESULTS, "strata_and_death_genes.rds")
PATH_E08_CSV <- file.path(DIR_TABLES,  "death_genes_by_stratum.csv")

MIN_STRATUM_N <- 30L
LOW_EXPR_PCT  <- 0.25    # flag, never filter

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

a      <- readRDS(file.path(DIR_RESULTS, "correlation_atlas.rds"))
sd_    <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc     <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw     <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
frames <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
mito   <- readRDS(PATH_TCGA_MITO)
myc_t  <- readRDS(PATH_TCGA_MYC)$estimators
cdc    <- readr::read_csv(PATH_CDC, show_col_types = FALSE, progress = FALSE)

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
STR_T <- .build_strata(frames, "TCGA", ID_T)
STR_S <- .build_strata(frames, "SCAN-B", ID_S)
tibble::tibble(stratum = names(STR_T), TCGA = lengths(STR_T),
               `SCAN-B` = lengths(STR_S)) %>%
  as.data.frame() %>% print(row.names = FALSE)

tcga_lin <- readRDS(PATH_TCGA_LINEAR); scanb_lin <- readRDS(PATH_SCANB_LINEAR)
LT <- tcga_lin$mat[, ID_T, drop = FALSE]; LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

.in_t <- function(g) intersect(unique(g), rownames(LT))
.in_s <- function(g) {
  h <- sc$symbol_map[unique(g)]
  intersect(unname(ifelse(is.na(h), unique(g), h)), rownames(LS))
}
.gene_rows <- function(genes, L, inf) {
  h <- vapply(genes, function(g) { r <- inf(g)
    if (length(r) != 1L) NA_character_ else r }, character(1))
  ok <- !is.na(h); M <- L[h[ok], , drop = FALSE]; rownames(M) <- genes[ok]
  list(mat = M, missing = genes[!ok])
}

# The two axes everything is read against, per cohort.
.axes <- function(gsva_new, arms_obj, mb, ids) {
  m <- rbind(MYC            = as.numeric(gsva_new["FELSHER_61", ids]),
             MYC_low_entang = as.numeric(gsva_new["MYC_UP.V1_UP", ids]),
             M_b            = as.numeric(mb[ids]),
             OXPHOS         = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]),
             OXPHOS_mitopps = as.numeric(arms_obj$mitopps_arms["OXPHOS subunits", ids]),
             MITORIBO       = as.numeric(arms_obj$gsva_arms["Mitochondrial ribosome", ids]))
  colnames(m) <- ids; m
}
AX_T <- .axes(nw$tcga_gsva_new, mito, stats::setNames(myc_t$M_b, myc_t$patient), ID_T)
AX_S <- .axes(sc$gsva_new,      sc,   sc$M_b,                                    ID_S)

COH <- list(
  TCGA     = list(L = LT, inf = .in_t, ax = AX_T, str = STR_T, ids = ID_T),
  `SCAN-B` = list(L = LS, inf = .in_s, ax = AX_S, str = STR_S, ids = ID_S))

# =============================================================================
# 2. Q-a: the ER-negative gap
# =============================================================================
# Read straight out of E03's atlas - no recomputation - across the whole
# estimator panel, so the gap is a property of the stratum and not of one
# signature. Luminal is absent from the atlas until E03 is re-run; LumA and
# LumB are there and are shown separately.
message("\n2. Q-a: the ER-negative gap")

gap_panel <- a$atlas %>%
  dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                adjusted == "raw", kind == "signature (GSVA)") %>%
  dplyr::group_by(cohort, stratum) %>%
  dplyr::summarise(n_samples = dplyr::first(n), n_sig = dplyr::n(),
                   median_rho = stats::median(rho),
                   min_rho = min(rho), max_rho = max(rho), .groups = "drop")
gap_panel %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Is it ERneg, or is it Basal wearing ERneg's name? Basal is ~72% of TCGA's
# ERneg and ~67% of SCAN-B's, so the two strata are not independent.
message("\n   ERneg is mostly Basal - the overlap, so the two are not read as",
        " independent:")
frames %>% dplyr::filter(!is.na(ER), !is.na(PAM50)) %>%
  dplyr::count(cohort, ER, basal = PAM50 == "Basal") %>%
  tidyr::pivot_wider(names_from = basal, values_from = n,
                     names_prefix = "basal_") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Is it range restriction? A stratum with less spread in either variable can
# only produce a smaller correlation, whatever the biology.
message("\n   spread of each axis within each stratum (IQR of the score):")
range_test <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  dplyr::bind_rows(lapply(names(C$str), function(st) {
    ids <- C$str[[st]]
    if (length(ids) < MIN_STRATUM_N) return(NULL)
    tibble::tibble(cohort = coh, stratum = st, n = length(ids),
                   iqr_MYC = stats::IQR(C$ax["MYC", ids]),
                   iqr_OXPHOS = stats::IQR(C$ax["OXPHOS", ids]))
  }))
}))
range_test %>%
  dplyr::left_join(dplyr::select(gap_panel, cohort, stratum, median_rho),
                   by = c("cohort", "stratum")) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# And is the gap OXPHOS-specific, or is the whole mitochondrial programme less
# MYC-coupled in ER-negative disease?
message("\n   every arm, ERpos minus ERneg (FELSHER_61, GSVA, raw):")
arm_gap <- a$atlas %>%
  dplyr::filter(instrument == "gsva", adjusted == "raw",
                myc_estimator == "FELSHER_61", measure_class == "arm",
                stratum %in% c("ERpos", "ERneg")) %>%
  dplyr::select(cohort, arm, stratum, rho) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho) %>%
  dplyr::mutate(gap = ERpos - ERneg)
arm_gap %>%
  tidyr::pivot_wider(id_cols = arm, names_from = cohort,
                     values_from = c(ERpos, ERneg, gap)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(gap_TCGA)) %>% as.data.frame() %>%
  print(row.names = FALSE)

# =============================================================================
# 3. Q-b: what drives the pro-death / pro-survival anticorrelation
# =============================================================================
message("\n3. Q-b: decomposing the D1 contrast")

PD <- sd_$cdc_sets$CDC_PRODEATH_APOPTOSIS
PS <- sd_$cdc_sets$CDC_PROSURVIVAL_APOPTOSIS
CDC_ANNOT <- cdc %>%
  dplyr::select(gene = human_symbol, effect, pathway, is_mitochondrial,
                family_pathway, dplyr::any_of(c("is_core", "confidence",
                                                "evidence_score")))

.gene_axis_rho <- function(coh, genes, strata = "all") {
  C <- COH[[coh]]
  gr <- .gene_rows(genes, C$L, C$inf)
  dplyr::bind_rows(lapply(strata, function(st) {
    ids <- C$str[[st]]
    if (is.null(ids) || length(ids) < MIN_STRATUM_N) return(NULL)
    .atlas_block(C$ax, gr$mat, ids, NULL, min_n = MIN_STRATUM_N) %>%
      dplyr::rename(axis = myc_estimator, gene = measure) %>%
      dplyr::mutate(cohort = coh, stratum = st)
  }))
}
contrast_genes <- dplyr::bind_rows(
  lapply(names(COH), .gene_axis_rho, genes = c(PD, PS), strata = "all")) %>%
  dplyr::left_join(CDC_ANNOT, by = "gene") %>%
  dplyr::filter(effect %in% c("pro-death", "pro-survival"))

message("\n   distribution of per-gene rho within each stratum of the curation:")
contrast_genes %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::group_by(cohort, axis, effect) %>%
  dplyr::summarise(n = dplyr::n(), median = stats::median(rho),
                   q25 = stats::quantile(rho, .25),
                   q75 = stats::quantile(rho, .75),
                   frac_positive = mean(rho > 0), .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   If `frac_positive` differs but the medians barely do, the contrast",
        " is a\n   SHIFT OF THE WHOLE DISTRIBUTION and not a few strong genes.")

# Is it a minority of genes? Sort each stratum by rho and ask what fraction of
# its mean is carried by the top decile.
message("\n   is the contrast carried by a minority? (share of the mean held by",
        " the top 10% of genes)")
concentration <- contrast_genes %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::group_by(cohort, effect) %>%
  dplyr::summarise(
    n = dplyr::n(), mean_rho = mean(rho),
    top10_share = { v <- sort(abs(rho), decreasing = TRUE)
                    sum(v[seq_len(ceiling(0.1 * length(v)))]) / sum(v) },
    .groups = "drop")
concentration %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   0.10 would mean perfectly even; 1.00 would mean one decile carries",
        " everything.")

# Does the contrast survive deleting the strongest contributors?
message("\n   the contrast after deleting the top N contributors from each side:")
.trimmed_contrast <- function(df, k) {
  df %>% dplyr::group_by(cohort, axis, effect) %>%
    dplyr::arrange(dplyr::desc(abs(rho)), .by_group = TRUE) %>%
    dplyr::slice(-seq_len(k)) %>%
    dplyr::summarise(m = mean(rho), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = effect, values_from = m) %>%
    dplyr::mutate(contrast = `pro-death` - `pro-survival`, dropped = k)
}
trimmed <- dplyr::bind_rows(lapply(c(0L, 10L, 25L, 50L, 100L),
                                   function(k) .trimmed_contrast(contrast_genes, k)))
trimmed %>% dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  dplyr::arrange(axis, cohort, dropped) %>% as.data.frame() %>%
  print(row.names = FALSE)

# Does it track anything the curation already knows?
message("\n   does the per-gene rho track the curation's own annotation?")
annot_test <- contrast_genes %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::group_by(cohort, effect) %>%
  dplyr::summarise(
    vs_is_mitochondrial = stats::cor(rho,
      as.numeric(is_mitochondrial %in% c(TRUE, "TRUE")), method = "spearman"),
    vs_evidence = if ("evidence_score" %in% names(contrast_genes) &&
                      sum(!is.na(evidence_score)) > 10)
      stats::cor(rho, evidence_score, method = "spearman",
                 use = "pairwise.complete.obs") else NA_real_,
    .groups = "drop")
annot_test %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the 12 strongest contributors on each side (SELECTED ON THE",
        " STATISTIC -\n   these are descriptions, not findings):")
contrast_genes %>% dplyr::filter(axis == "OXPHOS", cohort == "TCGA") %>%
  dplyr::group_by(effect) %>%
  dplyr::arrange(dplyr::desc(abs(rho)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  dplyr::select(effect, gene, rho, family_pathway, is_mitochondrial) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4. Q-c: the overlay genes, by stratum, with expression
# =============================================================================
message("\n4. Q-c: the overlay genes by subtype, ranked by expression")

CICD_GENES <- sort(unique(unlist(
  sd_$cdc_sets[grepl("CICD", names(sd_$cdc_sets))], use.names = FALSE)))
OVERLAY <- sort(unique(c(CICD_GENES, sd_$bcl2_family, sd_$family_labels$gene)))
STRATA_WANTED <- c("all", "ERpos", "ERneg", "Luminal", "LumA", "LumB", "Basal",
                   "HER2", "Normal")
message("   ", length(OVERLAY), " genes x ", length(STRATA_WANTED),
        " strata x 6 axes x 2 cohorts")

overlay_by_stratum <- dplyr::bind_rows(
  lapply(names(COH), .gene_axis_rho, genes = OVERLAY, strata = STRATA_WANTED))

# Expression, and the flag that says whether a rho is worth reading at all.
expr_rank <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  e <- rowMeans(C$L)
  pct <- rank(e) / length(e)
  gr <- .gene_rows(OVERLAY, C$L, C$inf)
  tibble::tibble(cohort = coh, gene = rownames(gr$mat),
                 expr_mean = e[rownames(gr$mat)],
                 expr_pct = pct[rownames(gr$mat)]) %>%
    dplyr::mutate(expr_decile = ceiling(10 * expr_pct),
                  low_expression = expr_pct < LOW_EXPR_PCT)
}))
overlay_by_stratum <- overlay_by_stratum %>%
  dplyr::left_join(expr_rank, by = c("cohort", "gene")) %>%
  dplyr::left_join(CDC_ANNOT, by = "gene") %>%
  dplyr::mutate(is_bcl2 = gene %in% sd_$bcl2_family,
                is_cicd = gene %in% CICD_GENES) %>%
  dplyr::select(cohort, gene, axis, stratum, n, rho, ci_lo, ci_hi, expr_mean,
                expr_pct, expr_decile, low_expression, effect, pathway,
                family_pathway, is_mitochondrial, is_bcl2, is_cicd)

message("\n   overlay genes below the ", LOW_EXPR_PCT * 100,
        "th expression percentile - flagged, not removed:")
expr_rank %>% dplyr::filter(low_expression) %>%
  dplyr::mutate(expr_pct = round(expr_pct, 3)) %>%
  dplyr::select(cohort, gene, expr_pct) %>%
  dplyr::arrange(cohort, expr_pct) %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   the BCL2 family against MYC, by stratum, TCGA:")
overlay_by_stratum %>%
  dplyr::filter(is_bcl2, axis == "MYC", cohort == "TCGA") %>%
  tidyr::pivot_wider(id_cols = c(gene, expr_decile), names_from = stratum,
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  dplyr::arrange(dplyr::desc(all)) %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   and against OXPHOS, TCGA:")
overlay_by_stratum %>%
  dplyr::filter(is_bcl2, axis == "OXPHOS", cohort == "TCGA") %>%
  tidyr::pivot_wider(id_cols = c(gene, expr_decile), names_from = stratum,
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  dplyr::arrange(dplyr::desc(all)) %>% as.data.frame() %>% print(row.names = FALSE)

# Which genes actually differ between the luminal and basal compartments?
message("\n   largest Luminal-minus-Basal difference, both cohorts agreeing on sign:")
lum_basal <- overlay_by_stratum %>%
  dplyr::filter(stratum %in% c("Luminal", "Basal"), axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(cohort, gene, axis, stratum, rho, low_expression) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho) %>%
  dplyr::mutate(diff = Luminal - Basal) %>%
  dplyr::select(cohort, gene, axis, Luminal, Basal, diff, low_expression) %>%
  tidyr::pivot_wider(id_cols = c(gene, axis), names_from = cohort,
                     values_from = c(Luminal, Basal, diff, low_expression)) %>%
  dplyr::filter(!is.na(diff_TCGA), !is.na(`diff_SCAN-B`),
                sign(diff_TCGA) == sign(`diff_SCAN-B`))
lum_basal %>%
  dplyr::arrange(dplyr::desc(abs(diff_TCGA) + abs(`diff_SCAN-B`))) %>%
  utils::head(15) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Figures
# =============================================================================
message("\n5. figures")

COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
theme_e08 <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey92",
                                                          colour = NA),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8, colour = "grey40",
                                                      hjust = 0))
.save <- function(p, nm, w, h) {
  for (ext in c("png", "pdf"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0(nm, ".", ext)), p,
                    width = w, height = h, dpi = 300)
  message("   ", nm); invisible(p)
}

g1 <- contrast_genes %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  ggplot2::ggplot(ggplot2::aes(rho, colour = effect)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::stat_ecdf(linewidth = 0.6) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_colour_manual(values = c(`pro-death` = "#d7191c",
                                          `pro-survival` = "#2c7bb6"), name = NULL) +
  ggplot2::labs(
    title = "What the D1 contrast looks like at gene resolution",
    subtitle = "EXPLORATORY - not pre-registered | per-gene Spearman rho, all samples",
    x = "per-gene Spearman rho", y = "cumulative fraction of the stratum",
    caption = paste("Two curves separated along their whole length is a shift",
                    "of the entire distribution - the annotation is doing work",
                    "across\nhundreds of genes. Two curves that meet except in",
                    "one tail would mean a handful of genes carry it.")) +
  theme_e08
.save(g1, "E08_fig1_d1_gene_distributions", 8, 5)

g2 <- overlay_by_stratum %>%
  dplyr::filter(is_bcl2, axis %in% c("MYC", "OXPHOS"),
                stratum %in% c("all", "Luminal", "Basal", "ERpos", "ERneg")) %>%
  dplyr::mutate(stratum = factor(stratum, levels = c("all", "ERpos", "ERneg",
                                                     "Luminal", "Basal")),
                gene = stats::reorder(gene, rho),
                cohort = factor(cohort, levels = names(COHORT_COLS))) %>%
  ggplot2::ggplot(ggplot2::aes(rho, gene, colour = cohort,
                               shape = low_expression)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6),
                      size = 1.5) +
  ggplot2::facet_grid(axis ~ stratum) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                              name = "bottom quartile of expression") +
  ggplot2::labs(
    title = "The BCL2 family by subtype",
    subtitle = paste("EXPLORATORY - not pre-registered | Luminal = LumA + LumB",
                     "(696 TCGA / 2,436 SCAN-B)"),
    x = "Spearman rho", y = NULL,
    caption = paste("Crosses mark genes in the bottom quartile of expression,",
                    "where a correlation is largely a correlation of",
                    "quantisation noise.\nThey are flagged rather than removed",
                    "so the reader can decide.")) +
  theme_e08
.save(g2, "E08_fig2_bcl2_by_stratum", 11, 6)

g3 <- gap_panel %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                stratum = stats::reorder(stratum, median_rho)) %>%
  ggplot2::ggplot(ggplot2::aes(median_rho, stratum, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = min_rho, xmax = max_rho),
                           position = ggplot2::position_dodge(width = 0.6),
                           size = 0.3, linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_samples), x = -0.02),
                     position = ggplot2::position_dodge(width = 0.6),
                     size = 2.2, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "The ER-negative gap, across the whole estimator panel",
    subtitle = paste("EXPLORATORY - not pre-registered | point is the median of",
                     "18 MYC signatures, bar is their full range"),
    x = "Spearman rho with OXPHOS subunits (GSVA, unadjusted)", y = NULL,
    caption = paste("The bar is the SPREAD ACROSS SIGNATURES, not a confidence",
                    "interval. ERneg is mostly Basal, so those two rows are not",
                    "independent\nreadings - see the overlap table in the",
                    "script output.")) +
  theme_e08
.save(g3, "E08_fig3_er_gap_panel", 8, 5)

# =============================================================================
# 6. Save
# =============================================================================
message("\n6. save")
saveRDS(list(gap_panel = gap_panel, range_test = range_test, arm_gap = arm_gap,
             contrast_genes = contrast_genes, concentration = concentration,
             trimmed = trimmed, annot_test = annot_test,
             overlay_by_stratum = overlay_by_stratum, expr_rank = expr_rank,
             lum_basal = lum_basal,
             settings = list(min_stratum_n = MIN_STRATUM_N,
                             low_expr_pct = LOW_EXPR_PCT,
                             strata = STRATA_WANTED),
             rules = list(
               luminal = paste("Luminal = LumA + LumB, added 2026-09-01;",
                               "reported beside LumA and LumB so a difference",
                               "between them stays visible"),
               erneg = paste("ERneg is ~70% Basal in both cohorts; the two",
                             "strata are not independent readings"),
               expression = paste("expr_pct is the gene's percentile in its",
                                  "own cohort's matrix. Below the 25th it is",
                                  "FLAGGED, never filtered"),
               selection = paste("the top-contributor lists are selected on the",
                                 "statistic they report and are descriptions,",
                                 "not findings")),
             built = Sys.time()), PATH_E08)
readr::write_csv(overlay_by_stratum, PATH_E08_CSV)
message("\nE08: done.")
message("    results/strata_and_death_genes.rds")
message("    outputs/tables/death_genes_by_stratum.csv")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E08)

  # Q-a
  x$gap_panel %>% as.data.frame()
  x$range_test %>% as.data.frame()
  x$arm_gap %>% dplyr::arrange(dplyr::desc(gap)) %>% as.data.frame()

  # Q-b: is it a shift or a tail?
  x$concentration %>% as.data.frame()
  x$trimmed %>% dplyr::filter(axis == "OXPHOS") %>% as.data.frame()

  # Q-c
  x$overlay_by_stratum %>%
    dplyr::filter(is_bcl2, axis == "OXPHOS", !low_expression) %>%
    tidyr::pivot_wider(id_cols = gene, names_from = c(cohort, stratum),
                       values_from = rho) %>% as.data.frame()

  x$lum_basal %>% dplyr::arrange(dplyr::desc(abs(diff_TCGA))) %>%
    utils::head(20) %>% as.data.frame()

  # which overlay genes should not be read at all
  x$expr_rank %>% dplyr::filter(low_expression) %>% as.data.frame()

}
