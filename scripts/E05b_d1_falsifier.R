# E05b_d1_falsifier.R
# =============================================================================
# THE FALSIFIER FOR D1, RUN BEFORE D1 IS BELIEVED.
#
# Declared in: docs/2026-09-01_phase1_celldeath_findings.md, section D1
# Same posture as E03b: this script exists to KILL D1 if D1 is wrong, and
# whatever it returns is recorded beside the finding rather than folded into it.
#
# D1 says: tumours with a high OXPHOS programme carry a pro-survival-skewed
# apoptotic transcriptome. The pro-death and pro-survival apoptosis strata (502
# and 584 genes, both only ~7% mitochondrial so D0 does not explain them) go
# opposite ways, contrast -0.413 / -0.551 against OXPHOS subunits, surviving
# proliferation adjustment at -0.320 / -0.372 and purity adjustment intact.
# Against MYC the same contrast collapses under proliferation adjustment, which
# is why D1 is stated as an OXPHOS finding.
#
# =============================================================================
# THE GAP THIS SCRIPT CLOSES, WHICH E05 LEFT OPEN
# =============================================================================
# E05's expression-matched null was computed against MYC ONLY - `.null_test`
# was called with FELSHER__MITOSTRIP as the target. **The OXPHOS axis, the one
# D1's surviving claim rests on, was never nulled at all.** That is the hole.
#
# And a null per stratum is not even the right test. D1 reports a CONTRAST, so
# the null has to be a null OF THE CONTRAST: draw an expression-matched set for
# each stratum, disjointly, as the real strata are, and take the difference.
# Section 5 does that.
#
# FOUR WAYS D1 COULD BE WRONG, ALL RUN HERE
# -----------------------------------------
#   A  The two strata differ in EXPRESSION, and the contrast follows from that
#      rather than from annotation. Section 2 measures it, section 5 nulls it.
#   B  The contrast is what any two same-size, same-expression gene sets would
#      give. Section 5.
#   C  It is carried by the mitochondrial minority after all - D0 leaking in
#      through 7% of the genes. Section 6 deletes every MitoCarta gene and
#      redoes it.
#   D  It is proliferation, as the MYC version of it turned out to be.
#      Every section runs raw AND proliferation-adjusted.
#
# NOTE ON WHAT IS BEING MEASURED. D1's headline number is a correlation of GSVA
# SCORES. Everything here is a mean of PER-GENE correlations, which is a
# different and more conservative statistic - E05 already showed the two diverge
# (the GSVA contrast is -0.41 where the gene-level null z was under 2). The
# gene-level contrast is the honest quantity and it is what the verdict uses.
#
# SCALE: linear DESeq2-normalised gene expression; every correlation is
# rank-based. SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE05b: the falsifier for D1\n", strrep("=", 78))

PATH_D1 <- file.path(DIR_RESULTS, "d1_falsifier.rds")
N_DRAWS <- 2000L
N_BINS  <- 20L

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

sd_   <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc    <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw    <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito  <- readRDS(PATH_TCGA_MITO)
dth   <- readRDS(file.path(DIR_RESULTS, "celldeath_axis.rds"))

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
tcga_lin <- readRDS(PATH_TCGA_LINEAR); scanb_lin <- readRDS(PATH_SCANB_LINEAR)
LT <- tcga_lin$mat[, ID_T, drop = FALSE]; LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

.in_t <- function(g) intersect(unique(g), rownames(LT))
.in_s <- function(g) {
  h <- sc$symbol_map[unique(g)]
  intersect(unname(ifelse(is.na(h), unique(g), h)), rownames(LS))
}

PRODEATH    <- sd_$cdc_sets$CDC_PRODEATH_APOPTOSIS
PROSURVIVAL <- sd_$cdc_sets$CDC_PROSURVIVAL_APOPTOSIS
stopifnot(length(intersect(PRODEATH, PROSURVIVAL)) == 0L)
message("   pro-death ", length(PRODEATH), " genes | pro-survival ",
        length(PROSURVIVAL), " | disjoint, as `effect` is exclusive")

MITO_GENES <- unique(unlist(mito$mito_paths, use.names = FALSE))

# The targets. MYC is carried so the contrast that DID collapse under
# proliferation adjustment stays visible beside the one that did not.
.targets <- function(gsva_new, arms_obj, ids) rbind(
  MYC            = as.numeric(gsva_new[MYC_REF, ids]),
  OXPHOS_gsva    = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]),
  OXPHOS_mitopps = as.numeric(arms_obj$mitopps_arms["OXPHOS subunits", ids]))
TGT_T <- .targets(nw$tcga_gsva_new, mito, ID_T)
TGT_S <- .targets(sc$gsva_new,      sc,   ID_S)
colnames(TGT_T) <- ID_T; colnames(TGT_S) <- ID_S

PROLIF_T <- matrix(mito$gsva_cov["PROLIF_DISJOINT", ID_T], ncol = 1,
                   dimnames = list(ID_T, "PROLIF_DISJOINT"))
PROLIF_S <- matrix(sc$gsva_cov["PROLIF_DISJOINT", ID_S], ncol = 1,
                   dimnames = list(ID_S, "PROLIF_DISJOINT"))

COHORTS <- list(
  TCGA   = list(L = LT, inf = .in_t, tgt = TGT_T, cov = PROLIF_T, ids = ID_T),
  `SCAN-B` = list(L = LS, inf = .in_s, tgt = TGT_S, cov = PROLIF_S, ids = ID_S))

# =============================================================================
# 2. FALSIFIER A - do the two strata differ in expression at all?
# =============================================================================
# The premise of the whole worry. If pro-survival genes are systematically more
# abundant than pro-death genes, a contrast between them is partly a contrast
# between expression levels. Measured, not assumed.
message("\n2. falsifier A: are the two strata matched on expression already?")

expression_profile <- dplyr::bind_rows(lapply(names(COHORTS), function(coh) {
  C <- COHORTS[[coh]]; expr <- log2(rowMeans(C$L) + 1)
  dplyr::bind_rows(
    tibble::tibble(cohort = coh, stratum = "pro-death",
                   value = expr[C$inf(PRODEATH)]),
    tibble::tibble(cohort = coh, stratum = "pro-survival",
                   value = expr[C$inf(PROSURVIVAL)]),
    tibble::tibble(cohort = coh, stratum = "background", value = expr))
}))
expression_profile %>% dplyr::group_by(cohort, stratum) %>%
  dplyr::summarise(n = dplyr::n(), median = round(stats::median(value), 3),
                   mean = round(mean(value), 3),
                   q25 = round(stats::quantile(value, .25), 3),
                   q75 = round(stats::quantile(value, .75), 3),
                   .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   Both strata sit ABOVE the background - curated genes are",
        " better expressed.\n   What matters for D1 is whether they differ",
        " FROM EACH OTHER; the null in\n   section 5 does not care either way,",
        " because it matches each separately.")

# =============================================================================
# 3. Per-gene correlations, every target, raw and proliferation-adjusted
# =============================================================================
message("\n3. per-gene rho over the whole matrix, 3 targets x 2 adjustments")

GRHO <- list()
for (coh in names(COHORTS)) {
  C <- COHORTS[[coh]]
  for (tg in rownames(C$tgt)) for (adj in c("raw", "prolif")) {
    GRHO[[paste(coh, tg, adj)]] <-
      .per_gene_rho(C$L, C$tgt[tg, ],
                    cov = if (adj == "prolif") C$cov else NULL)
  }
  message("   ", coh, " done")
}

# =============================================================================
# 4. Each stratum against its own expression-matched null - INCLUDING OXPHOS
# =============================================================================
# E05 ran this for MYC only. These are the cells that were missing.
message("\n4. each stratum against an expression-matched null, all targets")

set.seed(PROJECT_SEED)
.bins_for <- function(coh, key) {
  C <- COHORTS[[coh]]
  .expression_bins(C$L, which(is.finite(GRHO[[key]])), N_BINS)
}
stratum_null <- dplyr::bind_rows(lapply(names(GRHO), function(key) {
  parts <- strsplit(key, " ", fixed = TRUE)[[1]]
  coh <- parts[1]; tg <- parts[2]; adj <- parts[3]
  C <- COHORTS[[coh]]; g <- GRHO[[key]]; B <- .bins_for(coh, key)
  dplyr::bind_rows(lapply(list(`pro-death` = PRODEATH,
                               `pro-survival` = PROSURVIVAL), function(genes) {
    idx <- match(C$inf(genes), rownames(C$L)); idx <- idx[!is.na(idx)]
    idx <- idx[is.finite(g[idx])]
    obs <- mean(g[idx])
    d <- vapply(seq_len(N_DRAWS), function(i) mean(g[.matched_draw(idx, B)]),
                numeric(1))
    tibble::tibble(n = length(idx), observed = obs, null_mean = mean(d),
                   z = (obs - mean(d)) / stats::sd(d),
                   percentile = mean(d <= obs))
  }), .id = "stratum") %>%
    dplyr::mutate(cohort = coh, target = tg, adjusted = adj)
}))
stratum_null %>%
  dplyr::select(cohort, target, adjusted, stratum, n, observed, z, percentile) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. FALSIFIER B - the null OF THE CONTRAST, which is what D1 actually reports
# =============================================================================
# Each draw takes an expression-matched set for the pro-death stratum and a
# DISJOINT expression-matched set for the pro-survival stratum, and differences
# their mean per-gene rho. That is the same operation D1 performs, on genes that
# carry no annotation at all.
#
# If the observed contrast sits inside this distribution, D1 is what any two
# same-size, same-expression gene sets would have given, and it is dead.
message("\n5. falsifier B: the null of the contrast itself")

set.seed(PROJECT_SEED)
.contrast_null <- function(coh, tg, adj, drop_mito = FALSE) {
  key <- paste(coh, tg, adj); C <- COHORTS[[coh]]; g <- GRHO[[key]]
  bad <- if (drop_mito) match(C$inf(MITO_GENES), rownames(C$L)) else integer(0)
  bad <- bad[!is.na(bad)]
  keep <- setdiff(which(is.finite(g)), bad)
  B <- .expression_bins(C$L, keep, N_BINS)
  .idx <- function(genes) {
    i <- match(C$inf(genes), rownames(C$L)); i <- i[!is.na(i)]
    setdiff(i[is.finite(g[i])], bad)
  }
  i_pd <- .idx(PRODEATH); i_ps <- .idx(PROSURVIVAL)
  obs <- mean(g[i_pd]) - mean(g[i_ps])
  d <- vapply(seq_len(N_DRAWS), function(i) {
    a <- .matched_draw(i_pd, B)
    b <- .matched_draw(i_ps, B, exclude = a)
    mean(g[a]) - mean(g[b])
  }, numeric(1))
  tibble::tibble(cohort = coh, target = tg, adjusted = adj,
                 mito_dropped = drop_mito,
                 n_prodeath = length(i_pd), n_prosurvival = length(i_ps),
                 observed = obs, null_mean = mean(d), null_sd = stats::sd(d),
                 null_lo = unname(stats::quantile(d, 0.025)),
                 null_hi = unname(stats::quantile(d, 0.975)),
                 z = (obs - mean(d)) / stats::sd(d),
                 percentile = mean(d <= obs),
                 outside_95 = obs < stats::quantile(d, 0.025) |
                              obs > stats::quantile(d, 0.975))
}
grid <- expand.grid(cohort = names(COHORTS), target = rownames(TGT_T),
                    adjusted = c("raw", "prolif"),
                    stringsAsFactors = FALSE)
contrast_null <- dplyr::bind_rows(
  Map(function(a, b, c) .contrast_null(a, b, c),
      grid$cohort, grid$target, grid$adjusted))
contrast_null %>%
  dplyr::select(cohort, target, adjusted, observed, null_mean, null_lo, null_hi,
                z, outside_95) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. FALSIFIER C - delete every mitochondrial gene and do it again
# =============================================================================
# D0 says the death axis is mostly mitochondrial content. Both strata are only
# ~7% MitoCarta, but 7% of 500 genes is 35 genes that sit inside the arm being
# correlated against. This removes them from the strata AND from the background.
message("\n6. falsifier C: the same contrast with every MitoCarta gene deleted")

set.seed(PROJECT_SEED)
contrast_nomito <- dplyr::bind_rows(
  Map(function(a, b, c) .contrast_null(a, b, c, drop_mito = TRUE),
      grid$cohort, grid$target, grid$adjusted))
contrast_nomito %>%
  dplyr::select(cohort, target, adjusted, n_prodeath, n_prosurvival, observed,
                null_lo, null_hi, z, outside_95) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. Verdict
# =============================================================================
# Computed from the numbers, not argued. D1's surviving claim is the OXPHOS one
# under proliferation adjustment, so that is the cell the criterion names.
message("\n7. verdict")

.crit <- function(tab, tg, adj, mito) {
  x <- tab %>% dplyr::filter(target == tg, adjusted == adj,
                             mito_dropped == mito)
  all(x$outside_95) && all(x$observed < 0) && nrow(x) == 2L
}
verdict <- tibble::tibble(
  falsifier = c(
    "B  contrast vs a matched-pair null, OXPHOS GSVA, unadjusted",
    "B  contrast vs a matched-pair null, OXPHOS GSVA, proliferation-adjusted",
    "B  contrast vs a matched-pair null, OXPHOS mitoPPS, proliferation-adjusted",
    "C  same, with every MitoCarta gene deleted",
    "D  the MYC axis, proliferation-adjusted (expected to FAIL - D1 says so)"),
  criterion = "observed contrast negative and outside the 95% null, BOTH cohorts",
  result = c(
    ifelse(.crit(contrast_null,  "OXPHOS_gsva",    "raw",    FALSE), "SURVIVES", "FALSIFIED"),
    ifelse(.crit(contrast_null,  "OXPHOS_gsva",    "prolif", FALSE), "SURVIVES", "FALSIFIED"),
    ifelse(.crit(contrast_null,  "OXPHOS_mitopps", "prolif", FALSE), "SURVIVES", "FALSIFIED"),
    ifelse(.crit(contrast_nomito,"OXPHOS_gsva",    "prolif", TRUE),  "SURVIVES", "FALSIFIED"),
    ifelse(.crit(contrast_null,  "MYC",            "prolif", FALSE), "SURVIVES", "FALSIFIED")))
verdict %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   The last row is a CONTROL, not a test. D1 already says the MYC",
        " contrast is\n   proliferation; if it 'SURVIVES' here, D1's own",
        " framing is what needs revising.")

# =============================================================================
# 8. Figure and save
# =============================================================================
message("\n8. figure and save")

plotdat <- dplyr::bind_rows(contrast_null, contrast_nomito) %>%
  dplyr::mutate(
    panel = ifelse(mito_dropped, "MitoCarta genes deleted", "all genes"),
    cohort = factor(cohort, levels = c("TCGA", "SCAN-B")),
    lab = paste(target, adjusted, sep = " / "))
p <- ggplot2::ggplot(plotdat, ggplot2::aes(observed, lab, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = null_lo, xmax = null_hi),
                          height = 0.25, linewidth = 0.4, colour = "grey55",
                          position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::geom_point(size = 2.2,
                      position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::facet_wrap(~ panel) +
  ggplot2::scale_colour_manual(values = c(TCGA = "#1f4e79",
                                          `SCAN-B` = "#b8541a"), name = NULL) +
  ggplot2::labs(
    title = "D1 against a null of the contrast itself",
    subtitle = paste("EXPLORATORY - not pre-registered |", N_DRAWS,
                     "matched-pair draws; grey bars are the 95% null interval"),
    x = "pro-death minus pro-survival, mean per-gene Spearman rho", y = NULL,
    caption = paste("Each draw takes an expression-matched set for each",
                    "stratum, disjointly as the real strata are, and",
                    "differences them.\nA point inside its grey bar is a",
                    "contrast any two same-size, same-expression gene sets",
                    "would have given.")) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8,
                                                      colour = "grey40",
                                                      hjust = 0))
for (ext in c("png", "pdf"))
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("E05b_d1_contrast_null.", ext)),
                  p, width = 9, height = 4.5, dpi = 300)

saveRDS(list(expression_profile = expression_profile,
             stratum_null = stratum_null, contrast_null = contrast_null,
             contrast_nomito = contrast_nomito, verdict = verdict,
             settings = list(draws = N_DRAWS, bins = N_BINS,
                             seed = PROJECT_SEED),
             rules = list(
               purpose = paste("this script exists to KILL D1 if D1 is wrong;",
                               "the verdict is recorded beside the finding,",
                               "not folded into it"),
               declared_in = "docs/2026-09-01_phase1_celldeath_findings.md D1",
               gap_closed = paste("E05's expression-matched null was computed",
                                  "against MYC only; the OXPHOS axis, which is",
                                  "what D1 rests on, had never been nulled"),
               statistic = paste("mean per-gene rho, not the GSVA score",
                                 "correlation. The two diverge and the",
                                 "gene-level one is the conservative choice")),
             built = Sys.time()), PATH_D1)

message("\nE05b: done.")
message("    results/d1_falsifier.rds")
message("    outputs/figures/E05b_d1_contrast_null.{png,pdf}")
message("    Record the verdict in the D1 note BEFORE moving on.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_D1)
  x$verdict %>% as.data.frame()

  # the whole contrast grid, which is the answer
  x$contrast_null %>%
    dplyr::select(cohort, target, adjusted, observed, null_lo, null_hi, z,
                  outside_95) %>% as.data.frame()

  # and with the mitochondrial genes gone
  x$contrast_nomito %>%
    dplyr::select(cohort, target, adjusted, n_prodeath, observed, z,
                  outside_95) %>% as.data.frame()

  # the cells E05 never computed: each stratum against OXPHOS
  x$stratum_null %>% dplyr::filter(target != "MYC") %>%
    dplyr::select(cohort, target, adjusted, stratum, observed, z) %>%
    as.data.frame()

  # are the strata actually different in expression?
  x$expression_profile %>%
    ggplot2::ggplot(ggplot2::aes(value, colour = stratum)) +
    ggplot2::stat_ecdf() + ggplot2::facet_wrap(~ cohort) +
    ggplot2::labs(x = "log2 mean linear expression")

}
