# E05_celldeath_axis.R
# =============================================================================
# THE DEATH AXIS ON THE MYC-OXPHOS PLANE.
#
# Built to: docs/2026-08-31_phase1_plan.md section 2 (E05), section 1.3
# The plane it sits on: docs/2026-08-31_phase1_atlas_findings.md
#
# Three layers, deliberately distinct because they are trustworthy to very
# different degrees:
#
#   LAYER 1  the CDC strata - pro-death / pro-survival x apoptosis / CICD,
#            with and without the mitochondrial restriction. 6 of 8 scored;
#            the other 2 are 4 and 2 genes and are shown as GENE LISTS.
#   LAYER 2  the 15 Tang modalities; 12 clear n >= 15 and are scored.
#   LAYER 3  the individual genes - the BCL2 family and the family_pathway
#            labels - as points on the plane rather than as scores.
#
# =============================================================================
# THE TWO SENTENCES THAT MUST NOT BE WRITTEN FROM THIS SCRIPT
# =============================================================================
# 1. "MYC and OXPHOS INTERACT to determine the death programme."
#    They may or may not. THIS SCRIPT DOES NOT TEST THAT. The validation study
#    (myc_human_validation @ d3ac60e) tested the MYC x OXPHOS interaction on
#    apoptotic priming properly, pre-registered, and found it NULL. A
#    correlation with a product term is not an interaction test - the product of
#    two z-scores is large in the both-HIGH corner and equally large in the
#    both-LOW corner. It is carried here because the plan asks for it, and it is
#    labelled `product` everywhere. CLAUDE.md trap 1.
#
# 2. "CICD tracks MYC."  CICD is 13 pro-death and 4 pro-survival human genes.
#    It is the axis of most interest and the weakest measured (trap 9). The
#    13-gene set is scored; the 4- and 2-gene sets are NOT, ever, and every
#    CICD gene is also shown individually so a reader can see how thin it is.
#
# 3. "The mitochondrial death machinery tracks OXPHOS."
#    CDC_PRODEATH_APOPTOSIS_MITO correlates +0.80 with OXPHOS subunits, and
#    CDC_PROSURVIVAL_APOPTOSIS_MITO correlates +0.80 as well - IN THE SAME
#    DIRECTION, which is the tell. The `_MITO` strata are MitoCarta-restricted
#    by construction, so they are made of the same transcripts as the arm they
#    are being correlated against. That correlation is close to tautological and
#    says nothing about death. What is informative is the CONTRAST between the
#    pro-death and pro-survival strata, not either one's level, and the
#    unrestricted strata are where the contrast can be read.
#
# AND THE FOURTH, WHICH IS TRAP 10: TANG_FERROPTOSIS (600 genes, 3.3% of the
# matrix), TANG_APOPTOSIS (608) and TANG_AUTOPHAGY_DEPENDENT_CELL_DEATH (876,
# 4.8%) are near-transcriptome-wide. A correlation with any of them is close to
# a correlation with general expression. Section 6 measures exactly that, with
# an expression-matched null, and no number from those three sets should be read
# before looking at it.
#
# SCALE: GSVA scores are log-derived, the gene-level work is linear. Everything
# is rank-transformed by the engine, so no comparison here carries a scale
# error - see functions/correlation_engine.R.
# SPECIES: human, natively. The death sets carry first-class human columns and
# no ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE05: the cell-death axis\n", strrep("=", 78))

# =============================================================================
# 0. Constants
# =============================================================================
PATH_DEATH     <- file.path(DIR_RESULTS, "celldeath_axis.rds")
PATH_DEATH_CSV <- file.path(DIR_TABLES,  "celldeath_axis.csv")

MIN_STRATUM_N <- 30L
NULL_DRAWS    <- 2000L    # expression-matched draws per set per cohort
NULL_BINS     <- 20L      # ventiles of mean linear expression
HUGE_SET      <- 500L     # the trap-10 threshold, from E02

# The three MYC estimators the product terms are built from: the least
# entangled, the validation study's, and the most widely used.
PRODUCT_MYC <- c("MYC_UP.V1_UP", "FELSHER_61", "HALLMARK_MYC_TARGETS_V1")
PRODUCT_ARM <- "OXPHOS subunits"

INSTRUMENTS <- c(gsva = "gsva_arms", mitopps = "mitopps_arms",
                 content = "content_arms", zmean = "zmean_arms")

.lab_exploratory <- function(extra = NULL)
  paste(c("EXPLORATORY - not pre-registered; hypothesis-generating only",
          extra), collapse = " | ")
COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")

theme_atlas <- function(base = 10) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey92",
                                                            colour = NA),
                   strip.text = ggplot2::element_text(face = "bold",
                                                      size = base - 1),
                   plot.title = ggplot2::element_text(face = "bold",
                                                      size = base + 2),
                   plot.subtitle = ggplot2::element_text(size = base - 1,
                                                         colour = "grey30"),
                   plot.caption = ggplot2::element_text(size = base - 2,
                                                        colour = "grey40",
                                                        hjust = 0),
                   legend.position = "bottom")
}
.save_fig <- function(p, name, w, h) {
  for (ext in c("png", "pdf"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0(name, ".", ext)), p,
                    width = w, height = h, dpi = 300, limitsize = FALSE)
  message("   ", name)
  invisible(p)
}

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

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
DEATH_SETS <- c(names(sd_$cdc_sets)[sd_$cdc_table$scored],
                names(sd_$tang_sets)[sd_$tang_table$scored])
stopifnot(all(DEATH_SETS %in% rownames(nw$tcga_gsva_new)),
          all(DEATH_SETS %in% rownames(sc$gsva_new)))
message("   ", length(DEATH_SETS), " death sets scored (",
        sum(sd_$cdc_table$scored), " CDC + ", sum(sd_$tang_table$scored),
        " Tang); ", sum(!sd_$cdc_table$scored), " CDC strata deliberately not")

# The sets that are NOT scored, carried as gene lists. Trap 9.
CICD_UNSCORED <- sd_$cdc_sets[!sd_$cdc_table$scored]
for (nmx in names(CICD_UNSCORED))
  message("   NOT SCORED, shown as genes: ", nmx, " = ",
          paste(CICD_UNSCORED[[nmx]], collapse = ", "))

# =============================================================================
# 2. The three axes a death score is read against
# =============================================================================
message("\n2. the partners")

.stack_arms <- function(obj) {
  ARMS <- rownames(obj$gsva_arms)
  do.call(rbind, lapply(names(INSTRUMENTS), function(i) {
    x <- obj[[INSTRUMENTS[[i]]]][ARMS, , drop = FALSE]
    rownames(x) <- paste0("arm.", i, "::", ARMS); x }))
}
MYC_SIGS <- sd_$myc_panel$signature

.myc_rows <- function(gsva_new, extra) rbind(
  gsva_new[MYC_SIGS, , drop = FALSE], do.call(rbind, extra))

MYC_T <- .myc_rows(nw$tcga_gsva_new,
                   list(M_b = myc_t$M_b[match(ID_T, myc_t$patient)],
                        log2MYC = nw$tcga_log2MYC[ID_T],
                        M_c_call = as.numeric(myc_t$M_c_call[match(ID_T, myc_t$patient)])))
MYC_S <- .myc_rows(sc$gsva_new, list(M_b = sc$M_b[ID_S], log2MYC = sc$log2MYC[ID_S]))
colnames(MYC_T) <- ID_T; colnames(MYC_S) <- ID_S

# The product term. z-scored WITHIN COHORT and over ALL samples, so the same
# coordinate is used in every stratum. It is NOT an interaction test - see the
# header - and it is large in the both-low corner as well as the both-high one.
.products <- function(mycm, armsm) {
  z <- function(v) as.numeric(scale(v))
  out <- list()
  for (m in PRODUCT_MYC) for (i in names(INSTRUMENTS)) {
    out[[paste0("product::", m, " x ", i)]] <-
      z(mycm[m, ]) * z(armsm[paste0("arm.", i, "::", PRODUCT_ARM), ])
  }
  do.call(rbind, out)
}
ARM_T <- .stack_arms(mito); ARM_S <- .stack_arms(sc)
PRD_T <- .products(MYC_T, ARM_T); PRD_S <- .products(MYC_S, ARM_S)
colnames(PRD_T) <- ID_T; colnames(PRD_S) <- ID_S

PART_T <- rbind(MYC_T, ARM_T, PRD_T)
PART_S <- rbind(MYC_S, ARM_S, PRD_S)
part_meta <- tibble::tibble(partner = rownames(PART_T)) %>%
  dplyr::mutate(
    partner_class = dplyr::case_when(
      grepl("^arm\\.", partner)     ~ "mito_arm",
      grepl("^product::", partner)  ~ "product",
      TRUE                          ~ "myc_estimator"),
    instrument = ifelse(partner_class == "mito_arm",
                        sub("^arm\\.([^:]+)::.*$", "\\1", partner), NA_character_),
    label = dplyr::case_when(
      partner_class == "mito_arm" ~ sub("^arm\\.[^:]+::", "", partner),
      partner_class == "product"  ~ sub("^product::", "", partner),
      TRUE                        ~ partner))
message("   ", nrow(PART_T), " partners: ",
        sum(part_meta$partner_class == "myc_estimator"), " MYC estimators + ",
        sum(part_meta$partner_class == "mito_arm"), " arm x instrument + ",
        sum(part_meta$partner_class == "product"), " product terms")

DEATH_T <- nw$tcga_gsva_new[DEATH_SETS, ID_T, drop = FALSE]
DEATH_S <- sc$gsva_new[DEATH_SETS, ID_S, drop = FALSE]

# =============================================================================
# 3. Strata and adjustments, identical to E03
# =============================================================================
message("\n3. strata and adjustments")

.strata <- function(coh, ids) {
  f <- frames[frames$cohort == coh, ]; f <- f[match(ids, f$sample_id), ]
  s <- list(all = ids,
            ERpos = ids[!is.na(f$ER) & f$ER == "ERpos"],
            ERneg = ids[!is.na(f$ER) & f$ER == "ERneg"])
  for (p in levels(f$PAM50)) s[[p]] <- ids[!is.na(f$PAM50) & f$PAM50 == p]
  s
}
STR_T <- .strata("TCGA", ID_T); STR_S <- .strata("SCAN-B", ID_S)

cov_t <- frames %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::select(sample_id, purity, leuko) %>%
  tibble::column_to_rownames("sample_id")
cov_t <- as.matrix(cov_t[ID_T, , drop = FALSE])
ADJ_T <- list(raw = NULL,
              prolif = matrix(mito$gsva_cov["PROLIF_DISJOINT", ID_T], ncol = 1,
                              dimnames = list(ID_T, "PROLIF_DISJOINT")),
              purity_leuko = cov_t)
ADJ_S <- list(raw = NULL,
              prolif = matrix(sc$gsva_cov["PROLIF_DISJOINT", ID_S], ncol = 1,
                              dimnames = list(ID_S, "PROLIF_DISJOINT")))

# =============================================================================
# 4. The grid
# =============================================================================
message("\n4. building the death grid")

.run <- function(coh, PART, DEATH, strata, adj) {
  out <- list()
  for (st in names(strata)) {
    ids <- strata[[st]]
    if (length(ids) < MIN_STRATUM_N) next
    for (aj in names(adj)) {
      b <- .atlas_block(PART, DEATH, ids, adj[[aj]], min_n = MIN_STRATUM_N)
      if (is.null(b)) next
      out[[paste(st, aj)]] <- b %>%
        dplyr::mutate(cohort = coh, stratum = st, adjusted = aj)
    }
  }
  dplyr::bind_rows(out)
}
death_grid <- dplyr::bind_rows(
  .run("TCGA",   PART_T, DEATH_T, STR_T, ADJ_T),
  .run("SCAN-B", PART_S, DEATH_S, STR_S, ADJ_S)) %>%
  dplyr::rename(partner = myc_estimator, death_set = measure) %>%
  dplyr::left_join(part_meta, by = "partner") %>%
  dplyr::left_join(sd_$coverage %>% dplyr::select(cohort, death_set = set,
                                                  set_frac = frac),
                   by = c("cohort", "death_set")) %>%
  dplyr::mutate(set_n = c(vapply(sd_$cdc_sets, length, integer(1)),
                          vapply(sd_$tang_sets, length, integer(1)))[death_set],
                huge = set_n >= HUGE_SET) %>%
  dplyr::select(cohort, death_set, set_n, huge, set_frac, partner,
                partner_class, label, instrument, stratum, adjusted, n, k_cov,
                rho, ci_lo, ci_hi, p)
message("   ", format(nrow(death_grid), big.mark = ","), " cells")

message("\n   death sets against MYC (FELSHER_61) and OXPHOS subunits (GSVA), raw:")
death_grid %>%
  dplyr::filter(stratum == "all", adjusted == "raw",
                partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(partner = ifelse(partner == "FELSHER_61", "MYC", "OXPHOS")) %>%
  tidyr::pivot_wider(id_cols = c(death_set, set_n, huge),
                     names_from = c(cohort, partner), values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(TCGA_MYC)) %>% as.data.frame() %>%
  print(row.names = FALSE)

# The contrast the CDC strata exist for. A level is confounded by how
# mitochondrial the set is; the DIFFERENCE between pro-death and pro-survival on
# the same axis is not, because both are drawn from the same curation.
message("\n   pro-death MINUS pro-survival, the contrast the strata exist for:")
cdc_contrast <- death_grid %>%
  dplyr::filter(stratum == "all", grepl("^CDC_", death_set),
                partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(axis = ifelse(partner == "FELSHER_61", "MYC", "OXPHOS"),
                effect = ifelse(grepl("PRODEATH", death_set), "prodeath",
                                "prosurvival"),
                stratum_kind = sub("^CDC_(PRODEATH|PROSURVIVAL)_", "", death_set)) %>%
  dplyr::select(cohort, adjusted, axis, stratum_kind, effect, rho) %>%
  tidyr::pivot_wider(names_from = effect, values_from = rho) %>%
  dplyr::filter(!is.na(prodeath), !is.na(prosurvival)) %>%
  dplyr::mutate(contrast = prodeath - prosurvival)
cdc_contrast %>% dplyr::filter(adjusted == "raw") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. CICD, gene by gene
# =============================================================================
# Trap 9. The 13-gene pro-death set is scored above; the 4- and 2-gene
# pro-survival sets never are. All of it is shown here at gene resolution, which
# is the only honest resolution for a 4-gene "programme".
message("\n5. CICD at gene resolution")

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
  ok <- !is.na(h)
  M <- L[h[ok], , drop = FALSE]; rownames(M) <- genes[ok]
  list(mat = M, missing = genes[!ok])
}

CICD_GENES <- sort(unique(unlist(sd_$cdc_sets[grepl("CICD", names(sd_$cdc_sets))],
                                 use.names = FALSE)))
OVERLAY_GENES <- sort(unique(c(CICD_GENES, sd_$bcl2_family,
                               sd_$family_labels$gene)))

.gene_partner_rho <- function(genes, L, inf, PART, ids, coh) {
  gr <- .gene_rows(genes, L, inf)
  if (length(gr$missing))
    message("   ", coh, ": ", length(gr$missing), " gene(s) absent - ",
            paste(utils::head(gr$missing, 8), collapse = ", "))
  keep <- c("FELSHER_61", "MYC_UP.V1_UP", "M_b", "log2MYC",
            "arm.gsva::OXPHOS subunits", "arm.mitopps::OXPHOS subunits",
            "arm.gsva::mtDNA-encoded OXPHOS")
  .atlas_block(PART[intersect(keep, rownames(PART)), , drop = FALSE],
               gr$mat, ids, NULL, min_n = MIN_STRATUM_N) %>%
    dplyr::rename(partner = myc_estimator, gene = measure) %>%
    dplyr::mutate(cohort = coh)
}
overlay <- dplyr::bind_rows(
  .gene_partner_rho(OVERLAY_GENES, LT, .in_t, PART_T, ID_T, "TCGA"),
  .gene_partner_rho(OVERLAY_GENES, LS, .in_s, PART_S, ID_S, "SCAN-B")) %>%
  dplyr::left_join(sd_$family_labels, by = "gene") %>%
  dplyr::mutate(is_cicd = gene %in% CICD_GENES,
                is_bcl2 = gene %in% sd_$bcl2_family)

message("\n   the CICD genes against MYC and OXPHOS (GSVA), both cohorts:")
overlay %>%
  dplyr::filter(is_cicd, partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(partner = ifelse(partner == "FELSHER_61", "MYC", "OXPHOS")) %>%
  tidyr::pivot_wider(id_cols = c(gene, effect, pathway), names_from = c(cohort, partner),
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(TCGA_MYC)) %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. Trap 10 - the near-transcriptome-wide sets against an expression-matched null
# =============================================================================
# The question is blunt: is a correlation with TANG_FERROPTOSIS (600 genes) any
# different from a correlation with 600 randomly chosen genes of the same
# expression? The null draws size- and expression-matched sets from the whole
# matrix, in ventiles of mean linear expression, without replacement.
#
# The observed statistic is the MEAN PER-GENE Spearman rho with MYC across the
# set - not the set's GSVA score, because a random set has no GSVA score here.
# That makes this a test of the set's GENES, which is the right question.
message("\n6. expression-matched null for every death set (trap 10)")

# Chunked so a 18,153 x 3,207 rank matrix never exists all at once. Spearman by
# construction: rank both sides, then a centred-and-scaled inner product.
.per_gene_rho <- function(L, myc, chunk = 2000L) {
  ry <- rank(myc); n <- length(ry)
  ry <- (ry - mean(ry)) / stats::sd(ry)
  out <- numeric(nrow(L)); names(out) <- rownames(L)
  for (s in seq(1, nrow(L), by = chunk)) {
    e  <- min(s + chunk - 1L, nrow(L))
    Rx <- .rank_rows(L[s:e, , drop = FALSE])
    sdx <- apply(Rx, 1L, stats::sd)
    Rx <- (Rx - rowMeans(Rx)) / sdx
    v <- as.numeric(Rx %*% ry) / (n - 1)
    v[sdx == 0] <- NA_real_          # a gene with no variance has no rho
    out[s:e] <- v
  }
  out
}
.null_test <- function(L, inf, myc, sets, coh) {
  message("   ", coh, ": per-gene rho over ", nrow(L), " genes")
  g_rho <- .per_gene_rho(L, myc)
  # A gene with zero variance cannot enter either the observed statistic or the
  # null, and leaving it in the background would bias the draws toward NA.
  ok <- which(is.finite(g_rho))
  if (length(ok) < nrow(L)) {
    message("   ", coh, ": ", nrow(L) - length(ok),
            " gene(s) with no variance excluded from set and null alike")
  }
  expr  <- rowMeans(L)
  bins  <- rep(NA_integer_, nrow(L))
  bins[ok] <- cut(rank(expr[ok], ties.method = "first"),
                  breaks = NULL_BINS, labels = FALSE)
  by_bin <- split(ok, bins[ok])
  set.seed(PROJECT_SEED)
  dplyr::bind_rows(lapply(names(sets), function(nmx) {
    idx <- match(inf(sets[[nmx]]), rownames(L)); idx <- idx[!is.na(idx)]
    idx <- idx[is.finite(g_rho[idx])]
    if (length(idx) < 5L) return(NULL)
    obs <- mean(g_rho[idx])
    want <- table(factor(bins[idx], levels = seq_len(NULL_BINS)))
    short <- vapply(seq_len(NULL_BINS), function(b)
      want[[b]] > length(by_bin[[as.character(b)]]), logical(1))
    if (any(short)) {
      stop("expression-matched null for '", nmx, "' in ", coh, ": ventile(s) ",
           paste(which(short), collapse = ", "), " hold fewer background genes ",
           "than the set needs. The null cannot be drawn without replacement.",
           call. = FALSE)
    }
    draws <- vapply(seq_len(NULL_DRAWS), function(d) {
      pick <- unlist(lapply(seq_len(NULL_BINS), function(b) {
        k <- want[[b]]
        if (k == 0L) return(integer(0))
        pool <- by_bin[[as.character(b)]]
        pool[sample.int(length(pool), k)]   # never sample(pool, k): length-1 trap
      }), use.names = FALSE)
      mean(g_rho[pick])
    }, numeric(1))
    tibble::tibble(cohort = coh, death_set = nmx, n_genes = length(idx),
                   observed = obs, null_mean = mean(draws),
                   null_sd = stats::sd(draws),
                   z = (obs - mean(draws)) / stats::sd(draws),
                   percentile = mean(draws <= obs))
  }))
}
ALL_DEATH_SETS <- c(sd_$cdc_sets, sd_$tang_sets)[DEATH_SETS]
null_test <- dplyr::bind_rows(
  .null_test(LT, .in_t, MYC_T["FELSHER_61", ], ALL_DEATH_SETS, "TCGA"),
  .null_test(LS, .in_s, MYC_S["FELSHER_61", ], ALL_DEATH_SETS, "SCAN-B")) %>%
  dplyr::mutate(huge = n_genes >= HUGE_SET)

message("\n   observed mean per-gene rho(MYC, gene) vs an expression-matched null:")
null_test %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(abs(z))) %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   READ THE THREE HUGE SETS ONLY THROUGH THIS TABLE (trap 10).")

# =============================================================================
# 7. Figures
# =============================================================================
message("\n7. figures")

plane_pts <- death_grid %>%
  dplyr::filter(stratum == "all", adjusted == "raw",
                partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(axis = ifelse(partner == "FELSHER_61", "MYC", "OXPHOS")) %>%
  tidyr::pivot_wider(id_cols = c(cohort, death_set, set_n, huge),
                     names_from = axis, values_from = rho) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                family = ifelse(grepl("^CDC_", death_set), "CDC stratum",
                                "Tang modality"))

g1 <- ggplot2::ggplot(plane_pts, ggplot2::aes(MYC, OXPHOS)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 3, linewidth = 0.3) +
  ggplot2::geom_point(ggplot2::aes(size = set_n, shape = huge, colour = family),
                      alpha = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = sub("^(CDC_|TANG_)", "", death_set)),
                     size = 2.2, hjust = -0.08, vjust = 0.4) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                              name = "n >= 500 (trap 10)") +
  ggplot2::scale_colour_manual(values = c(`CDC stratum` = "#5e3c99",
                                          `Tang modality` = "#e66101"), name = NULL) +
  ggplot2::scale_size_continuous(range = c(1, 5), name = "genes") +
  ggplot2::expand_limits(x = c(-0.2, 0.9)) +
  ggplot2::labs(
    title = "Where each death programme sits on the MYC-OXPHOS plane",
    subtitle = .lab_exploratory("rho with FELSHER_61 and with OXPHOS subunits, both GSVA, unadjusted"),
    x = "Spearman rho with MYC activity", y = "Spearman rho with OXPHOS subunits",
    caption = paste("Triangles are the three near-transcriptome-wide Tang sets",
                    "(600, 608 and 876 genes). Their position says little until",
                    "it is read\nagainst the expression-matched null - see the",
                    "null figure. NOTHING HERE IS AN INTERACTION TEST",
                    "(CLAUDE.md trap 1).")) +
  theme_atlas()
.save_fig(g1, "E05_fig1_death_on_the_plane", 10, 5.5)

forest <- death_grid %>%
  dplyr::filter(stratum == "all", adjusted %in% c("raw", "prolif"),
                partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(axis = ifelse(partner == "FELSHER_61", "vs MYC activity",
                              "vs OXPHOS subunits"),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                adjusted = factor(adjusted, levels = c("raw", "prolif"),
                                  labels = c("unadjusted", "proliferation-adjusted")),
                death_set = stats::reorder(death_set, rho))

g2 <- ggplot2::ggplot(forest, ggplot2::aes(rho, death_set, colour = cohort,
                                           shape = adjusted)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                           position = ggplot2::position_dodge(width = 0.7),
                           size = 0.22, linewidth = 0.4) +
  ggplot2::facet_wrap(~ axis) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
  ggplot2::labs(
    title = "Every scored death programme, against both axes",
    subtitle = .lab_exploratory("all samples; 6 CDC strata and 12 Tang modalities"),
    x = "Spearman rho", y = NULL,
    caption = paste("CDC_PROSURVIVAL_CICD (4 genes) and CDC_PROSURVIVAL_CICD_MITO",
                    "(2) are absent BY DESIGN - a 4-gene GSVA score is not a",
                    "programme\n(CLAUDE.md trap 9). They appear at gene",
                    "resolution in the CICD figure instead. The _MITO strata are",
                    "MitoCarta-restricted, so their high rho against\nOXPHOS is",
                    "close to tautological - read the pro-death vs pro-survival",
                    "CONTRAST, not either level.")) +
  theme_atlas()
.save_fig(g2, "E05_fig2_death_forest", 10, 6.5)

nullp <- null_test %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                death_set = stats::reorder(death_set, z))
g3 <- ggplot2::ggplot(nullp, ggplot2::aes(z, death_set, colour = cohort,
                                          shape = huge)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-2, 2), linetype = 3, linewidth = 0.3,
                      colour = "grey50") +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 2) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                              name = "n >= 500") +
  ggplot2::labs(
    title = "Is a death set's MYC correlation more than its genes' expression?",
    subtitle = .lab_exploratory(
      paste0(NULL_DRAWS, " size- and expression-matched draws, ", NULL_BINS,
             " ventiles, per set per cohort")),
    x = "z of observed mean per-gene rho(MYC, gene) against the matched null",
    y = NULL,
    caption = paste("A set near z = 0 correlates with MYC no more than any",
                    "same-size, same-expression set of genes would. This is the",
                    "only\nlegitimate way to read the three",
                    "near-transcriptome-wide Tang sets (triangles).")) +
  theme_atlas()
.save_fig(g3, "E05_fig3_expression_matched_null", 8, 6)

cicd <- overlay %>%
  dplyr::filter(is_cicd,
                partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(axis = ifelse(partner == "FELSHER_61", "vs MYC activity",
                              "vs OXPHOS subunits"),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                scored = ifelse(gene %in% sd_$cdc_sets$CDC_PRODEATH_CICD,
                                "in the scored 13-gene set", "NOT scored anywhere"),
                gene = stats::reorder(gene, rho))
g4 <- ggplot2::ggplot(cicd, ggplot2::aes(rho, gene, colour = cohort,
                                         shape = scored)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                           position = ggplot2::position_dodge(width = 0.6),
                           size = 0.28, linewidth = 0.45) +
  ggplot2::facet_wrap(~ axis) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
  ggplot2::labs(
    title = "CICD at gene resolution - the only honest resolution for it",
    subtitle = .lab_exploratory("13 pro-death and 4 pro-survival human genes, that is the whole axis"),
    x = "Spearman rho", y = NULL,
    caption = paste("CLAUDE.md trap 9: CICD is the axis of most interest and",
                    "the weakest measured. The pro-survival side is FOUR genes",
                    "and is never\nscored as a programme. Read the genes, not a",
                    "summary of them.")) +
  theme_atlas()
.save_fig(g4, "E05_fig4_cicd_genes", 9, 5.5)

bcl2 <- overlay %>%
  dplyr::filter(is_bcl2 | !is.na(family_pathway),
                partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits")) %>%
  dplyr::mutate(axis = ifelse(partner == "FELSHER_61", "MYC", "OXPHOS")) %>%
  tidyr::pivot_wider(id_cols = c(cohort, gene, family_pathway, effect, is_bcl2),
                     names_from = axis, values_from = rho) %>%
  dplyr::filter(!is.na(MYC), !is.na(OXPHOS)) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
g5 <- ggplot2::ggplot(bcl2, ggplot2::aes(MYC, OXPHOS)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(ggplot2::aes(colour = effect, size = is_bcl2), alpha = 0.75) +
  ggplot2::geom_text(data = dplyr::filter(bcl2, is_bcl2),
                     ggplot2::aes(label = gene), size = 2.2, hjust = -0.12) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.4),
                             name = "BCL2 family") +
  ggplot2::scale_colour_brewer(palette = "Set1", na.value = "grey70", name = NULL) +
  ggplot2::labs(
    title = "Individual death genes on the MYC-OXPHOS plane",
    subtitle = .lab_exploratory("every gene carrying a family_pathway label, plus the 15 BCL2-family genes"),
    x = "Spearman rho with MYC activity", y = "Spearman rho with OXPHOS subunits",
    caption = paste("Gene-level expression, linear scale, all samples,",
                    "unadjusted. Points, not scores - at these set sizes this is",
                    "what\n'overlaid on the correlation plots' can honestly",
                    "mean.")) +
  theme_atlas()
.save_fig(g5, "E05_fig5_genes_on_the_plane", 10, 5.5)

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. save")
saveRDS(list(
  death_grid = death_grid, cdc_contrast = cdc_contrast,
  overlay = overlay, null_test = null_test,
  part_meta = part_meta, death_sets = DEATH_SETS,
  cicd_unscored = CICD_UNSCORED, cicd_genes = CICD_GENES,
  settings = list(null_draws = NULL_DRAWS, null_bins = NULL_BINS,
                  huge_set = HUGE_SET, min_stratum_n = MIN_STRATUM_N,
                  product_myc = PRODUCT_MYC, product_arm = PRODUCT_ARM),
  rules = list(
    not_an_interaction = paste("the `product` partners are NOT an interaction",
                              "test. The validation study tested the MYC x",
                              "OXPHOS interaction on apoptotic priming,",
                              "pre-registered, and found it null. A product of",
                              "z-scores is large in the both-low corner too."),
    cicd = paste("pro-survival CICD is 4 genes and pro-survival CICD MITO is 2.",
                 "Neither is ever scored. Read the genes."),
    mito_strata = paste("the _MITO strata are MitoCarta-restricted, so their",
                        "correlation with a mitochondrial arm is close to",
                        "tautological - both _MITO strata sit near +0.80 with",
                        "OXPHOS subunits regardless of direction of effect.",
                        "Read the pro-death vs pro-survival CONTRAST, not the",
                        "level, and prefer the unrestricted strata for it."),
    huge = paste("TANG_FERROPTOSIS, TANG_APOPTOSIS and",
                 "TANG_AUTOPHAGY_DEPENDENT_CELL_DEATH are 3-5% of the matrix.",
                 "Read them only through null_test."),
    coverage = paste("TANG_FERROPTOSIS (0.70) and TANG_IMMUNOGENIC_CELL_DEATH",
                     "(0.62) are below the 0.80 coverage floor in BOTH cohorts",
                     "identically - the missing symbols are absent from the",
                     "gene sets' vocabulary, not from a platform"),
    exploratory = "nothing here is pre-registered; every cell is hypothesis-generating"),
  built = Sys.time()), PATH_DEATH)
readr::write_csv(death_grid, PATH_DEATH_CSV)

message("\nE05: done.")
message("    results/celldeath_axis.rds       grid, overlay, null")
message("    outputs/tables/celldeath_axis.csv ", format(nrow(death_grid),
        big.mark = ","), " rows")
message("    5 figures in outputs/figures/")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  d <- readRDS(PATH_DEATH)

  # the only legitimate first look at the three huge sets
  d$null_test %>% dplyr::arrange(dplyr::desc(abs(z))) %>% as.data.frame()

  # the CDC strata, which are the part of the death axis with a real hypothesis
  d$death_grid %>%
    dplyr::filter(grepl("^CDC_", death_set), stratum == "all",
                  partner %in% c("FELSHER_61", "arm.gsva::OXPHOS subunits",
                                 "arm.mitopps::OXPHOS subunits")) %>%
    tidyr::pivot_wider(id_cols = c(death_set, set_n),
                       names_from = c(cohort, partner, adjusted),
                       values_from = rho) %>% as.data.frame()
  # `partner`, not `label`: gsva and mitopps share the label "OXPHOS subunits"
  # and pivoting on it silently produces list-columns.

  # pro-death minus pro-survival, the contrast the strata exist for
  d$death_grid %>%
    dplyr::filter(stratum == "all", adjusted == "raw", partner == "FELSHER_61",
                  death_set %in% c("CDC_PRODEATH_APOPTOSIS",
                                   "CDC_PROSURVIVAL_APOPTOSIS",
                                   "CDC_PRODEATH_APOPTOSIS_MITO",
                                   "CDC_PROSURVIVAL_APOPTOSIS_MITO")) %>%
    dplyr::select(cohort, death_set, rho, ci_lo, ci_hi) %>% as.data.frame()

  # the product terms - and remember what they are not
  d$death_grid %>%
    dplyr::filter(partner_class == "product", stratum == "all",
                  adjusted == "raw") %>%
    tidyr::pivot_wider(id_cols = death_set, names_from = c(cohort, label),
                       values_from = rho) %>% as.data.frame()

  # BCL2 family, gene by gene
  d$overlay %>% dplyr::filter(is_bcl2) %>%
    tidyr::pivot_wider(id_cols = c(gene, effect), names_from = c(cohort, partner),
                       values_from = rho) %>% as.data.frame()

}
