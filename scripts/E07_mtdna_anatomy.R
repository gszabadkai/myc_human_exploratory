# E07_mtdna_anatomy.R
# =============================================================================
# PHASE 2. What organises the divergence WITHIN the mitochondrial genome?
#
# From the 2026-09-01 handoff section 5, candidate 1, and E03b's verdict:
# after conditioning on mtDNA content, MT-CO2 sits at +0.440 / +0.413 and
# MT-CO1 at -0.295 / -0.113, MT-ND1 at -0.198 / -0.131 and MT-ND5 at -0.266 /
# -0.337. MT-CO1 and MT-CO2 are ADJACENT on the same heavy-strand polycistron
# and go opposite ways. Their mutual correlation is unremarkable (0.77 / 0.78),
# so it is not obviously mismapping - but that was an absence of evidence, not
# a test.
#
# =============================================================================
# THE MOVE THAT MAKES THIS TRACTABLE
# =============================================================================
# Conditioning on mtDNA content removes what the 13 genes SHARE. What is left
# in each gene is its DEVIATION from the common mitochondrial-genome programme,
# and E03b showed those deviations are larger than the raw correlations (median
# spread 0.408 -> 0.717 in TCGA). So the object of study is the residual
# structure, and the questions are:
#
#   1. Is there ONE secondary axis, or several? (section 3, PCA of residuals)
#   2. Do MT-CO1, MT-ND1 and MT-ND5 sit on the SAME side of it, i.e. is CO1's
#      negative the same phenomenon as ND5's? (section 3)
#   3. Does that axis replicate between cohorts? (section 4)
#   4. Does it track POSITION on the genome, strand, or bicistronic grouping -
#      the things a processing or quantification explanation predicts?
#      (section 5)
#   5. What ELSE does it track, across the 18 arms and the MYC panel?
#      (section 6)
#
# MT-ND6 IS CARRIED BUT NEVER TRUSTED. E03b showed it correlates 0.12-0.36 with
# every other mtDNA gene in SCAN-B where every other pair sits at 0.6-0.94, and
# it is the only light-strand, non-polyadenylated protein gene. It is plotted so
# its anomaly stays visible and it is excluded from every summary statistic.
#
# SCALE: linear DESeq2-normalised; every correlation is rank-based.
# SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE07: the anatomy of the mtDNA divergence\n", strrep("=", 78))

PATH_E07 <- file.path(DIR_RESULTS, "mtdna_anatomy.rds")

# =============================================================================
# 0. The structural annotation the positional hypothesis needs
# =============================================================================
# ORDER along the heavy strand from the origin of replication, strand, and the
# two bicistronic pairs whose transcripts are not separated by a tRNA. Order and
# strand only - no coordinates are used, so a small coordinate error cannot
# change a result. Standard human chrM gene order; verify against the reference
# if it ever matters more than it does here.
#
# The two structural facts a processing explanation would lean on:
#   - ATP8/ATP6 and ND4L/ND4 are BICISTRONIC: one mRNA each, so the two members
#     of a pair cannot be transcribed apart. If the residual axis separates a
#     bicistronic pair, it is not transcription.
#   - ND6 is the only LIGHT-strand gene and the only one without a poly(A) tail.
MT_ANNOT <- tibble::tribble(
  ~gene,      ~order, ~strand, ~operon,      ~polyA,
  "MT-ND1",        1L, "H",    "ND1",        TRUE,
  "MT-ND2",        2L, "H",    "ND2",        TRUE,
  "MT-CO1",        3L, "H",    "CO1",        TRUE,
  "MT-CO2",        4L, "H",    "CO2",        TRUE,
  "MT-ATP8",       5L, "H",    "ATP8/ATP6",  TRUE,
  "MT-ATP6",       6L, "H",    "ATP8/ATP6",  TRUE,
  "MT-CO3",        7L, "H",    "CO3",        TRUE,
  "MT-ND3",        8L, "H",    "ND3",        TRUE,
  "MT-ND4L",       9L, "H",    "ND4L/ND4",   TRUE,
  "MT-ND4",       10L, "H",    "ND4L/ND4",   TRUE,
  "MT-ND5",       11L, "H",    "ND5",        TRUE,
  "MT-ND6",       12L, "L",    "ND6",        FALSE,
  "MT-CYB",       13L, "H",    "CYB",        TRUE)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

a    <- readRDS(file.path(DIR_RESULTS, "correlation_atlas.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)
f    <- readRDS(file.path(DIR_RESULTS, "mtdna_falsifiers.rds"))

# E03b's output is an INPUT here, and it is keyed by estimator name. If it was
# written before the 2026-09-01 relabelling it still says "FELSHER_61" where
# this script now looks for FELSHER__MITOSTRIP, and the failure downstream is a
# rename error three sections later that says nothing about the cause.
if (!MYC_REF %in% f$copy_number_test$myc_estimator) {
  stop("results/mtdna_falsifiers.rds does not contain the estimator '", MYC_REF,
       "'. It carries: ",
       paste(utils::head(sort(unique(f$copy_number_test$myc_estimator)), 4),
             collapse = ", "),
       ". That file predates the MYC estimator relabelling.\n",
       "  RE-RUN E03b, then this script.", call. = FALSE)
}

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
tcga_lin <- readRDS(PATH_TCGA_LINEAR); scanb_lin <- readRDS(PATH_SCANB_LINEAR)
LT <- tcga_lin$mat[, ID_T, drop = FALSE]; LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

MT_GENES <- sd_$arm_sets[["mtDNA-encoded OXPHOS"]]
stopifnot(setequal(MT_GENES, MT_ANNOT$gene))
MT_ANNOT <- MT_ANNOT[order(MT_ANNOT$order), ]
TRUSTED <- setdiff(MT_ANNOT$gene, "MT-ND6")

.in_t <- function(g) intersect(unique(g), rownames(LT))
.in_s <- function(g) {
  h <- sc$symbol_map[unique(g)]
  intersect(unname(ifelse(is.na(h), unique(g), h)), rownames(LS))
}
.mt_matrix <- function(L, inf) {
  h <- vapply(MT_ANNOT$gene, function(g) {
    r <- inf(g); if (length(r) != 1L) NA_character_ else r }, character(1))
  if (anyNA(h)) stop("mtDNA gene(s) absent: ",
                     paste(MT_ANNOT$gene[is.na(h)], collapse = ", "), call. = FALSE)
  M <- L[h, , drop = FALSE]; rownames(M) <- MT_ANNOT$gene; M
}
MT_T <- .mt_matrix(LT, .in_t); MT_S <- .mt_matrix(LS, .in_s)
message("   13 genes x ", ncol(MT_T), " (TCGA) and ", ncol(MT_S), " (SCAN-B)")

# =============================================================================
# 2. Content, and the residuals that are the object of study
# =============================================================================
# Content proxy: log2 of the mean linear expression of the 13, exactly as E03b
# used it. The residual matrix is ranks-with-content-projected-out, which is
# what makes a partial Spearman a plain Pearson on residuals downstream.
message("\n2. content and residuals")

.residuals <- function(M) {
  content <- log2(colMeans(M))
  Z <- cbind(`(Intercept)` = 1, content = rank(content))
  list(content = content, R = .resid_rows(.rank_rows(M), Z))
}
RES_T <- .residuals(MT_T); RES_S <- .residuals(MT_S)
message(sprintf("   content explains a median %.1f%% (TCGA) / %.1f%% (SCAN-B)",
        100 * stats::median(1 - apply(RES_T$R, 1, var) /
                              apply(.rank_rows(MT_T), 1, var)),
        100 * stats::median(1 - apply(RES_S$R, 1, var) /
                              apply(.rank_rows(MT_S), 1, var))))

# =============================================================================
# 3. Is there one secondary axis, and who is on which side of it?
# =============================================================================
message("\n3. the residual structure")

.resid_cor <- function(R) stats::cor(t(R))
CR_T <- .resid_cor(RES_T$R); CR_S <- .resid_cor(RES_S$R)
message("\n   residual correlation among the 13, TCGA (MT-ND6 shown, not trusted):")
print(round(CR_T, 2))

.pca <- function(R, trusted) {
  p <- stats::prcomp(t(R[trusted, , drop = FALSE]), center = TRUE, scale. = TRUE)
  list(pca = p,
       var = p$sdev^2 / sum(p$sdev^2),
       load = p$rotation[, 1:3, drop = FALSE],
       scores = p$x[, 1:3, drop = FALSE])
}
P_T <- .pca(RES_T$R, TRUSTED); P_S <- .pca(RES_S$R, TRUSTED)
message("\n   variance explained by the residual PCs (12 trusted genes):")
tibble::tibble(PC = paste0("PC", 1:3),
               TCGA = round(P_T$var[1:3], 3),
               `SCAN-B` = round(P_S$var[1:3], 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

loadings <- tibble::tibble(gene = TRUSTED) %>%
  dplyr::mutate(TCGA_PC1 = P_T$load[gene, 1], TCGA_PC2 = P_T$load[gene, 2],
                SCANB_PC1 = P_S$load[gene, 1], SCANB_PC2 = P_S$load[gene, 2]) %>%
  dplyr::left_join(MT_ANNOT, by = "gene")

# The MYC deviation each gene shows, from E03b, is the thing being explained.
dev <- f$copy_number_test %>%
  dplyr::filter(myc_estimator == MYC_REF, adjusted == "mtDNA_content") %>%
  dplyr::select(cohort, gene, dev = rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = dev) %>%
  dplyr::rename(dev_TCGA = TCGA, dev_SCANB = `SCAN-B`)
loadings <- loadings %>% dplyr::left_join(dev, by = "gene")

message("\n   residual PC loadings beside the MYC deviation they should explain:")
loadings %>% dplyr::arrange(dplyr::desc(dev_TCGA)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(gene, order, operon, TCGA_PC1, SCANB_PC1, TCGA_PC2, SCANB_PC2,
                dev_TCGA, dev_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)

axis_test <- tibble::tibble(
  comparison = c("PC1 loading vs MYC deviation, TCGA",
                 "PC1 loading vs MYC deviation, SCAN-B",
                 "PC2 loading vs MYC deviation, TCGA",
                 "PC2 loading vs MYC deviation, SCAN-B"),
  spearman = c(
    stats::cor(loadings$TCGA_PC1,  loadings$dev_TCGA,  method = "spearman"),
    stats::cor(loadings$SCANB_PC1, loadings$dev_SCANB, method = "spearman"),
    stats::cor(loadings$TCGA_PC2,  loadings$dev_TCGA,  method = "spearman"),
    stats::cor(loadings$SCANB_PC2, loadings$dev_SCANB, method = "spearman")))
message("\n   does a residual PC ORDER the genes the way MYC does?")
axis_test %>% dplyr::mutate(spearman = round(spearman, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4. Does the axis replicate between cohorts?
# =============================================================================
message("\n4. cross-cohort replication of the residual axis")
replication <- tibble::tibble(
  quantity = c("residual correlation matrix (66 trusted pairs)",
               "PC1 loadings", "PC2 loadings", "MYC deviation per gene"),
  spearman = c(
    stats::cor(CR_T[TRUSTED, TRUSTED][lower.tri(CR_T[TRUSTED, TRUSTED])],
               CR_S[TRUSTED, TRUSTED][lower.tri(CR_S[TRUSTED, TRUSTED])],
               method = "spearman"),
    stats::cor(loadings$TCGA_PC1, loadings$SCANB_PC1, method = "spearman"),
    stats::cor(loadings$TCGA_PC2, loadings$SCANB_PC2, method = "spearman"),
    stats::cor(loadings$dev_TCGA, loadings$dev_SCANB, method = "spearman")))
replication %>% dplyr::mutate(spearman = round(spearman, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   PC signs are arbitrary; a strong NEGATIVE here means the same axis",
        " with\n   the sign flipped, which is agreement, not disagreement.")

# =============================================================================
# 5. Position, strand and operon - what a processing explanation predicts
# =============================================================================
# If the divergence were transcript processing or 3' decay, the deviation should
# track ORDER along the polycistron. If it were quantification, it should
# concentrate on the overlapping and bicistronic genes.
#
# THE SHARPEST TEST IS THE BICISTRONIC PAIRS. ATP8/ATP6 and ND4L/ND4 are each a
# single mRNA. Two genes on one mRNA CANNOT be transcribed apart, so if the pair
# members sit far apart on the residual axis, no transcriptional explanation
# survives - and if they sit together while CO1 and CO2 sit apart, the CO1/CO2
# split is not explained by shared-transcript arithmetic either.
message("\n5. position, strand and the bicistronic pairs")

pos_test <- tibble::tibble(
  test = c("deviation vs order along the heavy strand, TCGA",
           "deviation vs order along the heavy strand, SCAN-B"),
  spearman = c(
    stats::cor(loadings$order, loadings$dev_TCGA,  method = "spearman"),
    stats::cor(loadings$order, loadings$dev_SCANB, method = "spearman")))
pos_test %>% dplyr::mutate(spearman = round(spearman, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the two bicistronic pairs, which share one mRNA each:")
loadings %>% dplyr::filter(operon %in% c("ATP8/ATP6", "ND4L/ND4")) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(gene, operon, dev_TCGA, dev_SCANB, TCGA_PC1, SCANB_PC1) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   and the adjacent pair that is NOT bicistronic, for contrast:")
loadings %>% dplyr::filter(gene %in% c("MT-CO1", "MT-CO2")) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(gene, operon, dev_TCGA, dev_SCANB, TCGA_PC1, SCANB_PC1) %>%
  as.data.frame() %>% print(row.names = FALSE)

within_operon <- tibble::tibble(
  pair = c("ATP8 vs ATP6", "ND4L vs ND4", "CO1 vs CO2"),
  TCGA_gap  = c(abs(diff(loadings$dev_TCGA[match(c("MT-ATP8","MT-ATP6"), loadings$gene)])),
                abs(diff(loadings$dev_TCGA[match(c("MT-ND4L","MT-ND4"), loadings$gene)])),
                abs(diff(loadings$dev_TCGA[match(c("MT-CO1","MT-CO2"), loadings$gene)]))),
  SCANB_gap = c(abs(diff(loadings$dev_SCANB[match(c("MT-ATP8","MT-ATP6"), loadings$gene)])),
                abs(diff(loadings$dev_SCANB[match(c("MT-ND4L","MT-ND4"), loadings$gene)])),
                abs(diff(loadings$dev_SCANB[match(c("MT-CO1","MT-CO2"), loadings$gene)]))),
  bicistronic = c(TRUE, TRUE, FALSE))
message("\n   gap in MYC deviation between adjacent genes:")
within_operon %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Expression level, the other thing that could organise it
expr_test <- tibble::tibble(gene = TRUSTED) %>%
  dplyr::mutate(expr_TCGA  = log2(rowMeans(MT_T)[gene] + 1),
                expr_SCANB = log2(rowMeans(MT_S)[gene] + 1)) %>%
  dplyr::left_join(dplyr::select(loadings, gene, dev_TCGA, dev_SCANB), by = "gene")
message("\n   deviation vs expression level: TCGA ",
        round(stats::cor(expr_test$expr_TCGA, expr_test$dev_TCGA,
                         method = "spearman"), 3),
        " | SCAN-B ",
        round(stats::cor(expr_test$expr_SCANB, expr_test$dev_SCANB,
                         method = "spearman"), 3))

# =============================================================================
# 6. What else does the secondary axis track?
# =============================================================================
# The residual PC1 is a per-SAMPLE score. Correlating it against the 18 arms and
# the MYC panel says what kind of tumour has a high one - which is the first
# real clue to what the axis is.
message("\n6. what the residual axis tracks across the rest of the atlas")

.axis_probe <- function(P, arms_obj, gsva_new, ids, coh) {
  AX <- rbind(residPC1 = P$scores[ids, 1], residPC2 = P$scores[ids, 2])
  colnames(AX) <- ids
  M <- rbind(arms_obj$gsva_arms[, ids, drop = FALSE],
             PROLIF = arms_obj$gsva_cov["PROLIF_DISJOINT", ids],
             MYC_ref = gsva_new[MYC_REF, ids],
             MYC_UP.V1_UP = gsva_new[MYC_LOW_ENTANG, ids])
  .atlas_block(AX, M, ids, NULL, min_n = 30L) %>%
    dplyr::rename(axis = myc_estimator, target = measure) %>%
    dplyr::mutate(cohort = coh) %>%
    dplyr::select(cohort, axis, target, rho, ci_lo, ci_hi)
}
axis_probe <- dplyr::bind_rows(
  .axis_probe(P_T, mito, nw$tcga_gsva_new, ID_T, "TCGA"),
  .axis_probe(P_S, sc,   sc$gsva_new,      ID_S, "SCAN-B"))
# BOTH axes, because which one matters is decided in section 4 and not here.
# In the 2026-09-01 run PC1 was the larger axis and did NOT replicate (loadings
# agree at -0.168); PC2 was smaller and did (0.629), and PC2 is the one that
# orders the genes the way MYC does. Printing only PC1 would show the axis that
# turned out to be cohort-specific.
for (ax in c("residPC1", "residPC2")) {
  message("\n   ", ax, ":")
  axis_probe %>% dplyr::filter(axis == ax) %>%
    tidyr::pivot_wider(id_cols = target, names_from = cohort,
                       values_from = rho) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::arrange(dplyr::desc(abs(TCGA))) %>% as.data.frame() %>%
    print(row.names = FALSE)
}

# =============================================================================
# 7. Figure and save
# =============================================================================
message("\n7. figure and save")

pl <- loadings %>%
  dplyr::select(gene, order, operon, dev_TCGA, dev_SCANB) %>%
  tidyr::pivot_longer(c(dev_TCGA, dev_SCANB), names_to = "cohort",
                      values_to = "dev") %>%
  dplyr::mutate(cohort = ifelse(cohort == "dev_TCGA", "TCGA", "SCAN-B"),
                bicistronic = operon %in% c("ATP8/ATP6", "ND4L/ND4"))
p <- ggplot2::ggplot(pl, ggplot2::aes(order, dev, colour = cohort)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_line(ggplot2::aes(group = cohort), linewidth = 0.4, alpha = 0.6) +
  ggplot2::geom_point(ggplot2::aes(shape = bicistronic), size = 2.4) +
  ggplot2::geom_text(data = dplyr::filter(pl, cohort == "TCGA"),
                     ggplot2::aes(label = sub("^MT-", "", gene)), size = 2.4,
                     vjust = -1.1, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(TCGA = "#1f4e79",
                                          `SCAN-B` = "#b8541a"), name = NULL) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                              name = "shares an mRNA") +
  ggplot2::scale_x_continuous(breaks = 1:13) +
  ggplot2::labs(
    title = "The mtDNA divergence is not positional",
    subtitle = paste("EXPLORATORY - not pre-registered | rho with MYC activity",
                     "after conditioning on mtDNA content"),
    x = "order along the heavy strand from the origin", y = "MYC deviation",
    caption = paste("If this were 3' decay or processing the line would trend.",
                    "The two bicistronic pairs (triangles) share one mRNA each",
                    "and\ncannot be transcribed apart, so the gap WITHIN a pair",
                    "is the scale of non-transcriptional variation to judge",
                    "CO1-vs-CO2 against.\nMT-ND6 (position 12) is the",
                    "light-strand, non-polyadenylated gene and is not trusted;",
                    "see E03b.")) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8, colour = "grey40",
                                                      hjust = 0))
for (ext in c("png", "pdf"))
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0("E07_mtdna_position.", ext)), p,
                  width = 8, height = 5, dpi = 300)

saveRDS(list(annot = MT_ANNOT, loadings = loadings, axis_test = axis_test,
             replication = replication, pos_test = pos_test,
             within_operon = within_operon, expr_test = expr_test,
             axis_probe = axis_probe,
             resid_cor = list(TCGA = CR_T, `SCAN-B` = CR_S),
             pca = list(TCGA = P_T$var, `SCAN-B` = P_S$var),
             rules = list(
               nd6 = paste("MT-ND6 is carried for visibility and excluded from",
                           "every summary statistic - see E03b falsifier 2"),
               annot = paste("gene ORDER and strand only; no coordinates are",
                             "used, so a coordinate error cannot change a",
                             "result"),
               residual = paste("conditioning on mtDNA content removes what the",
                                "13 share; the residual is each gene's",
                                "deviation from the common programme")),
             built = Sys.time()), PATH_E07)
message("\nE07: done.  results/mtdna_anatomy.rds")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E07)
  x$loadings %>% dplyr::arrange(dplyr::desc(dev_TCGA)) %>% as.data.frame()
  x$within_operon %>% as.data.frame()      # the decisive comparison
  x$replication %>% as.data.frame()
  round(x$resid_cor$TCGA, 2)
  round(x$resid_cor$`SCAN-B`, 2)
  x$axis_probe %>% dplyr::filter(axis == "residPC1") %>% as.data.frame()

}
