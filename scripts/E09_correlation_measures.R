# E09_correlation_measures.R
# =============================================================================
# PHASE 2. Is Spearman hiding anything?
#
# From the 2026-09-01 handoff section 5, the author's third question: every
# correlation in this study is a Spearman. Is there a specific advantage to
# that, or would another measure reveal something else?
#
# =============================================================================
# WHY SPEARMAN WAS CHOSEN, AND IT IS A GOOD REASON
# =============================================================================
# It is INVARIANT TO MONOTONE TRANSFORMS. CLAUDE.md names the log-versus-linear
# scale question as the most likely silent error in this repo: GSVA needs log
# VST, mitoPPS needs linear DESeq2, and they must not share an input object. A
# Spearman between a log score and a linear one is well defined and carries no
# scale error. A Pearson between them is not the same quantity depending on
# which side was logged. That alone justifies the default.
#
# It is also robust to the heavy right tail of linear expression - the mtDNA
# genes are orders of magnitude above the median gene - and to the handful of
# extreme samples every tumour cohort has.
#
# WHAT IT CANNOT SEE, WHICH IS THE POINT OF THIS SCRIPT
# -----------------------------------------------------
#   1. NON-MONOTONE dependence. A threshold or saturating MYC-OXPHOS
#      relationship reads as a weakened rho; a U-shape reads as ~0. MYC's
#      biology is dose-dependent, so this is a real blind spot and not a
#      theoretical one. Section 4 tests it directly.
#   2. MAGNITUDE. Ranks discard it. If MYC-amplified tumours are genuinely
#      high-leverage rather than outliers, Pearson would weight them and
#      Spearman will not. Section 3 compares.
#
# FOUR MEASURES, AND WHAT EACH ADDS
#   spearman  rank, monotone, scale-free            - the study's default
#   pearson   linear on the LOG scale               - uses magnitude, fragile to tails
#   bicor     Tukey biweight midcorrelation         - uses magnitude, robust to tails
#   kendall   rank, concordance-based               - a consistency check on spearman
# plus a NON-MONOTONE probe: the gain in R-squared from a natural spline over a
# straight line, and the decile profile that shows the shape.
#
# A DISAGREEMENT BETWEEN spearman AND bicor IS INFORMATIVE (magnitude matters);
# a disagreement between spearman AND pearson WITH bicor agreeing with spearman
# is a tail artefact.
#
# SCALE: this script deliberately uses the LOG (VST) matrix for the gene-level
# work and the GSVA/mitoPPS scores as built, because Pearson and bicor are NOT
# scale-free and the comparison is only meaningful on a stated scale.
# SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE09: does the choice of correlation measure change the answer?\n",
        strrep("=", 78))

PATH_E09 <- file.path(DIR_RESULTS, "correlation_measures.rds")
KENDALL_N <- 1200L    # tau is O(n^2); a subsample is enough for a consistency check
SPLINE_DF <- 4L

# =============================================================================
# 1. Measures
# =============================================================================
message("\n1. measures")

# Tukey biweight midcorrelation, the WGCNA default, implemented here rather
# than adding a dependency E00 would have to check for. Robust like a rank
# method but it keeps magnitude: observations far from the median in MAD units
# are downweighted smoothly and those beyond 9 MAD get weight zero.
.bicor <- function(x, y, const = 9) {
  .u <- function(v) {
    m <- stats::median(v); s <- stats::median(abs(v - m))
    if (s == 0) return(NULL)
    u <- (v - m) / (const * s)
    w <- (1 - u^2)^2 * (abs(u) < 1)
    z <- (v - m) * w
    d <- sqrt(sum(z^2)); if (d == 0) return(NULL)
    z / d
  }
  a <- .u(x); b <- .u(y)
  if (is.null(a) || is.null(b)) return(NA_real_)
  sum(a * b)
}

# How much better does a smooth curve fit than a straight line? On RANKS, so
# the comparison is about SHAPE and not about the marginal distributions - a
# monotone-but-curved relationship is already fully captured by Spearman, and
# this is asking whether anything is left after that.
.nonmonotone_gain <- function(x, y, df = SPLINE_DF) {
  rx <- rank(x); ry <- rank(y)
  lin <- stats::lm(ry ~ rx)
  spl <- stats::lm(ry ~ splines::ns(rx, df = df))
  c(r2_linear = summary(lin)$r.squared,
    r2_spline = summary(spl)$r.squared,
    gain      = summary(spl)$r.squared - summary(lin)$r.squared)
}

.all_measures <- function(x, y, seed = PROJECT_SEED) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  set.seed(seed)
  i <- if (length(x) > KENDALL_N) sample.int(length(x), KENDALL_N) else seq_along(x)
  g <- .nonmonotone_gain(x, y)
  tibble::tibble(n = length(x),
                 spearman = stats::cor(x, y, method = "spearman"),
                 pearson  = stats::cor(x, y),
                 bicor    = .bicor(x, y),
                 kendall  = stats::cor(x[i], y[i], method = "kendall"),
                 r2_linear = g[["r2_linear"]], r2_spline = g[["r2_spline"]],
                 spline_gain = g[["gain"]])
}

# =============================================================================
# 2. Inputs and the pairs to test
# =============================================================================
message("\n2. inputs")

sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)
myc_t <- readRDS(PATH_TCGA_MYC)$estimators
dth  <- readRDS(file.path(DIR_RESULTS, "celldeath_axis.rds"))

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
tcga_lin <- readRDS(PATH_TCGA_LINEAR); scanb_lin <- readRDS(PATH_SCANB_LINEAR)
LT <- tcga_lin$mat[, ID_T, drop = FALSE]; LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

X_ARMS <- c("OXPHOS subunits", "Mitochondrial ribosome", "mtDNA-encoded OXPHOS",
            "TCA cycle", "Fatty acid oxidation")
X_MYC  <- c(MYC_REF, MYC_LOW_ENTANG, MYC_HALLMARK)

.side <- function(gsva_new, arms_obj, mb, L, ids) {
  myc <- rbind(gsva_new[X_MYC, ids, drop = FALSE],
               M_b = as.numeric(mb[ids]),
               log2MYC = as.numeric(log2(L["MYC", ids])))
  oth <- rbind(
    do.call(rbind, lapply(c("gsva", "mitopps", "content", "zmean"), function(i) {
      m <- arms_obj[[paste0(i, "_arms")]][X_ARMS, ids, drop = FALSE]
      rownames(m) <- paste0(i, "::", X_ARMS); m })),
    death = gsva_new["CDC_PROSURVIVAL_APOPTOSIS", ids],
    death2 = gsva_new["CDC_PRODEATH_APOPTOSIS", ids])
  rownames(oth)[rownames(oth) == "death"]  <- "CDC_PROSURVIVAL_APOPTOSIS"
  rownames(oth)[rownames(oth) == "death2"] <- "CDC_PRODEATH_APOPTOSIS"
  list(myc = myc, oth = oth)
}
S_T <- .side(nw$tcga_gsva_new, mito, nw$tcga_M_b_variants[MB_REF, ID_T],
             LT, ID_T)
S_S <- .side(sc$gsva_new, sc, sc$M_b_variants[MB_REF, ], LS, ID_S)
message("   ", nrow(S_T$myc), " MYC estimators x ", nrow(S_T$oth),
        " partners x 2 cohorts = ", 2 * nrow(S_T$myc) * nrow(S_T$oth), " pairs")

# =============================================================================
# 3. Every pair under every measure
# =============================================================================
message("\n3. computing")

.sweep <- function(S, coh) {
  dplyr::bind_rows(lapply(rownames(S$myc), function(m) {
    dplyr::bind_rows(lapply(rownames(S$oth), function(o) {
      .all_measures(S$myc[m, ], S$oth[o, ]) %>%
        dplyr::mutate(cohort = coh, myc_estimator = m, partner = o)
    }))
  }))
}
measures <- dplyr::bind_rows(.sweep(S_T, "TCGA"), .sweep(S_S, "SCAN-B")) %>%
  dplyr::select(cohort, myc_estimator, partner, n, spearman, pearson, bicor,
                kendall, r2_linear, r2_spline, spline_gain)

message("\n   the headline pair, all four measures:")
measures %>%
  dplyr::filter(myc_estimator == MYC_REF,
                partner %in% paste0(c("gsva", "mitopps", "content", "zmean"),
                                    "::OXPHOS subunits")) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(cohort, partner, spearman, pearson, bicor, kendall) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   do the measures ever disagree? (across all ",
        nrow(measures), " pairs)")
agreement <- tibble::tibble(
  pair = c("spearman vs pearson", "spearman vs bicor", "spearman vs kendall",
           "pearson vs bicor"),
  spearman_of_the_two = c(
    stats::cor(measures$spearman, measures$pearson, method = "spearman"),
    stats::cor(measures$spearman, measures$bicor,   method = "spearman"),
    stats::cor(measures$spearman, measures$kendall, method = "spearman"),
    stats::cor(measures$pearson,  measures$bicor,   method = "spearman")),
  max_abs_difference = c(
    max(abs(measures$spearman - measures$pearson)),
    max(abs(measures$spearman - measures$bicor)),
    max(abs(measures$spearman - measures$kendall)),
    max(abs(measures$pearson  - measures$bicor))))
agreement %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   where spearman and pearson disagree most:")
measures %>% dplyr::mutate(gap = pearson - spearman) %>%
  dplyr::arrange(dplyr::desc(abs(gap))) %>% utils::head(10) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(cohort, myc_estimator, partner, spearman, pearson, bicor, gap) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   If bicor sits with SPEARMAN, the gap is a tail artefact and",
        " Spearman is right.\n   If bicor sits with PEARSON, magnitude carries",
        " real information Spearman discards.")

# =============================================================================
# 4. THE QUESTION SPEARMAN CANNOT ANSWER: is the relationship monotone?
# =============================================================================
message("\n4. non-monotonicity")

message("\n   largest spline gain over a straight line (on ranks):")
measures %>% dplyr::arrange(dplyr::desc(spline_gain)) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  dplyr::select(cohort, myc_estimator, partner, spearman, r2_linear, r2_spline,
                spline_gain) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   A gain near zero means Spearman lost nothing. A large gain with a",
        " MODEST\n   spearman is the signature of a threshold or a saturating",
        " curve.")

# The decile profile, which shows the shape rather than scoring it.
.profile <- function(x, y, k = 10L) {
  b <- cut(rank(x, ties.method = "first"), breaks = k, labels = FALSE)
  tibble::tibble(decile = seq_len(k),
                 mean_y = as.numeric(tapply(y, b, mean)),
                 se_y = as.numeric(tapply(y, b, function(v)
                   stats::sd(v) / sqrt(length(v)))))
}
PROFILE_PAIRS <- list(
  c(MYC_REF, "gsva::OXPHOS subunits"),
  c(MYC_LOW_ENTANG, "gsva::OXPHOS subunits"),
  c(MYC_REF, "gsva::Mitochondrial ribosome"),
  c(MYC_REF, "gsva::mtDNA-encoded OXPHOS"),
  c(MYC_REF, "CDC_PROSURVIVAL_APOPTOSIS"),
  c("log2MYC", "gsva::OXPHOS subunits"))
profiles <- dplyr::bind_rows(lapply(list(TCGA = S_T, `SCAN-B` = S_S),
  function(S) dplyr::bind_rows(lapply(PROFILE_PAIRS, function(p)
    .profile(S$myc[p[1], ], S$oth[p[2], ]) %>%
      dplyr::mutate(myc_estimator = p[1], partner = p[2])))), .id = "cohort")

message("\n   decile profile of the headline pair (mean OXPHOS per MYC decile):")
profiles %>%
  dplyr::filter(myc_estimator == MYC_REF,
                partner == "gsva::OXPHOS subunits") %>%
  tidyr::pivot_wider(id_cols = decile, names_from = cohort, values_from = mean_y) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Figures
# =============================================================================
message("\n5. figures")

COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
mlong <- measures %>%
  tidyr::pivot_longer(c(pearson, bicor, kendall), names_to = "measure",
                      values_to = "value") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
p1 <- ggplot2::ggplot(mlong, ggplot2::aes(spearman, value, colour = cohort)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2) +
  ggplot2::geom_point(size = 1.1, alpha = 0.6) +
  ggplot2::facet_wrap(~ measure) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "Does the choice of correlation measure change the answer?",
    subtitle = paste("EXPLORATORY - not pre-registered |", nrow(measures),
                     "pairs: 5 MYC estimators x 22 partners x 2 cohorts"),
    x = "Spearman rho (the study's default)", y = "the other measure",
    caption = paste("Points on the diagonal mean the choice is immaterial.",
                    "bicor departing from spearman while pearson does not is",
                    "impossible;\npearson departing while bicor does not is a",
                    "heavy-tail artefact, and Spearman is the safer reading.")) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8, colour = "grey40",
                                                      hjust = 0))
for (ext in c("png", "pdf"))
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("E09_measure_agreement.", ext)),
                  p1, width = 9, height = 4, dpi = 300)

p2 <- ggplot2::ggplot(profiles, ggplot2::aes(decile, mean_y, colour = cohort)) +
  ggplot2::geom_line(linewidth = 0.5) +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = mean_y - 1.96 * se_y,
                                        ymax = mean_y + 1.96 * se_y),
                           size = 0.2) +
  ggplot2::facet_wrap(~ paste(myc_estimator, "vs", partner), scales = "free_y",
                      ncol = 3) +
  ggplot2::scale_x_continuous(breaks = c(1, 5, 10)) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "Is the relationship monotone? The shape Spearman compresses to one number",
    subtitle = "EXPLORATORY - not pre-registered | mean partner score per decile of the MYC estimator",
    x = "decile of the MYC estimator", y = "mean partner score",
    caption = paste("A straight rise is what a rank correlation assumes. A",
                    "plateau at either end is a threshold or a saturation, and",
                    "Spearman\nwould report it as a weaker linear-ish",
                    "association rather than as the shape it is.")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8, colour = "grey40",
                                                      hjust = 0))
for (ext in c("png", "pdf"))
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("E09_decile_profiles.", ext)),
                  p2, width = 9, height = 6, dpi = 300)

# =============================================================================
# 6. Save
# =============================================================================
saveRDS(list(measures = measures, agreement = agreement, profiles = profiles,
             settings = list(kendall_subsample = KENDALL_N,
                             spline_df = SPLINE_DF, seed = PROJECT_SEED),
             rules = list(
               why_spearman = paste("invariance to monotone transforms is what",
                                    "makes a log score and a linear score",
                                    "comparable at all; CLAUDE.md's scale",
                                    "discipline is the reason it is the default"),
               reading = paste("bicor with spearman against pearson = tail",
                               "artefact; bicor with pearson against spearman =",
                               "magnitude carries information"),
               blind_spot = paste("Spearman sees only monotone association; the",
                                  "spline gain and the decile profile are the",
                                  "parts it cannot report")),
             built = Sys.time()), PATH_E09)
message("\nE09: done.  results/correlation_measures.rds + 2 figures")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  m <- readRDS(PATH_E09)
  m$agreement %>% as.data.frame()

  # the whole sweep, wide
  m$measures %>%
    dplyr::select(cohort, myc_estimator, partner, spearman, pearson, bicor,
                  kendall) %>% as.data.frame()

  # is anything non-monotone?
  m$measures %>% dplyr::arrange(dplyr::desc(spline_gain)) %>%
    dplyr::select(cohort, myc_estimator, partner, spearman, spline_gain) %>%
    utils::head(20) %>% as.data.frame()

  # and the shapes
  m$profiles %>% dplyr::filter(partner == "gsva::OXPHOS subunits") %>%
    as.data.frame()

}
