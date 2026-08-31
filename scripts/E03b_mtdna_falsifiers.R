# E03b_mtdna_falsifiers.R
# =============================================================================
# THE FALSIFIERS FOR F3, RUN BEFORE F3 IS BELIEVED.
#
# Declared in: docs/2026-08-31_phase1_atlas_findings.md, section F3
# Reads:       E02's scores, E03's atlas, the linear matrices
#
# F3 says the two genomes are regulated differently: nuclear-encoded OXPHOS
# rises with MYC activity (+0.43 / +0.39) while the mtDNA-encoded arm does not
# (+0.06 / +0.07) and turns NEGATIVE on mitoPPS (-0.081 / -0.088); and within
# the 13 genes, MT-CO2 (+0.240 / +0.243) and MT-ND5 (-0.160 / -0.143) go
# opposite ways despite sitting on the same heavy-strand polycistron.
#
# The findings note named three ways that could be wrong. This script runs them.
#
#   FALSIFIER 3 (named as the most likely benign explanation, and the first to
#   run): the mitoPPS negative is an ARITHMETIC IDENTITY. mitoPPS is a ratio
#   against the rest of the mitochondrial programme. If the nuclear OXPHOS
#   programme simply rises, everything else must fall in composition terms with
#   no biology at all.
#     -> section 2 holds the denominator fixed: re-query the mtDNA arm against
#        universes with the OXPHOS pathways removed.
#     -> section 3 abandons mitoPPS entirely and asks the same question with a
#        transparent log ratio anyone can recompute in three lines.
#
#   FALSIFIER 1: mtDNA COPY NUMBER. If MYC-high tumours carry fewer genomes,
#   all 13 genes should move together. They do not - but the formal test is to
#   put mtDNA content in as a covariate and see whether the SPREAD among the 13
#   collapses.
#     -> section 4.
#
#   FALSIFIER 2: QUANTIFICATION. MT-ND5 and MT-ND6 overlap on opposite strands,
#   so multi-mapping could manufacture exactly this pattern.
#     -> section 5 asks whether the 13x13 correlation structure singles those
#        two out.
#
# THIS SCRIPT CAN KILL F3. That is what it is for. Whatever it returns is
# recorded beside the finding, not folded into it.
#
# SCALE: linear DESeq2-normalised throughout (mitoPPS and the ratios both need
# it). Correlations are Spearman and therefore rank-invariant.
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "mitopps.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE03b: the falsifiers for F3\n", strrep("=", 78))

PATH_FALSIFIERS <- file.path(DIR_RESULTS, "mtdna_falsifiers.rds")

# Sample floor for a block. Every block here uses a whole cohort, so this is a
# guard rather than a filter.
MIN_N <- 30L

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

a     <- readRDS(file.path(DIR_RESULTS, "correlation_atlas.rds"))
sd_   <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc    <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw    <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito  <- readRDS(PATH_TCGA_MITO)
myc_t <- readRDS(PATH_TCGA_MYC)$estimators

tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(tcga_lin$scale,  "linear_deseq2_normalised"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
LT <- tcga_lin$mat[, ID_T, drop = FALSE]
LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

arm_sets   <- sd_$arm_sets
mito_paths <- mito$mito_paths
arm_univ   <- mito$arm_universe_path
MT_GENES   <- arm_sets[["mtDNA-encoded OXPHOS"]]
NUC_OX     <- arm_sets[["OXPHOS subunits"]]

# SCAN-B vocabulary, exactly as E02 built it (2014 UCSC build; trap 7).
.in_t <- function(g) intersect(unique(g), rownames(LT))
.in_s <- function(g) {
  h <- sc$symbol_map[unique(g)]
  intersect(unname(ifelse(is.na(h), unique(g), h)), rownames(LS))
}

# The 13 rows, under their CANONICAL names whatever the matrix calls them. Built
# by explicit lookup rather than by relying on intersect() preserving order,
# and it refuses if any of the 13 is missing - a 12-gene mtDNA arm would change
# every number below without changing any label.
.mt_matrix <- function(L, inf) {
  h <- vapply(MT_GENES, function(g) {
    r <- inf(g)
    if (length(r) != 1L) NA_character_ else r
  }, character(1))
  if (anyNA(h)) {
    stop("mtDNA gene(s) absent from the matrix: ",
         paste(MT_GENES[is.na(h)], collapse = ", "), call. = FALSE)
  }
  M <- L[h, , drop = FALSE]
  rownames(M) <- MT_GENES
  M
}
message("   mtDNA genes: ", length(.in_t(MT_GENES)), " (TCGA) / ",
        length(.in_s(MT_GENES)), " (SCAN-B)")

# The MYC estimators, exactly as E03 used them.
MYC_SIGS <- sd_$myc_panel$signature
EST_T <- rbind(nw$tcga_gsva_new[MYC_SIGS, ID_T, drop = FALSE],
               M_b = myc_t$M_b[match(ID_T, myc_t$patient)],
               log2MYC = nw$tcga_log2MYC[ID_T])
EST_S <- rbind(sc$gsva_new[MYC_SIGS, ID_S, drop = FALSE],
               M_b = sc$M_b[ID_S], log2MYC = sc$log2MYC[ID_S])
colnames(EST_T) <- ID_T; colnames(EST_S) <- ID_S

# =============================================================================
# 2. FALSIFIER 3a - hold the mitoPPS denominator fixed
# =============================================================================
# The canonical mtDNA mitoPPS is queried against the 142-pathway universe with
# its own pathway held out. That universe is heavily OXPHOS-loaded: 22 of its
# pathways are OXPHOS, a complex, its subunits or its assembly factors.
#
# If the negative is an identity forced by nuclear OXPHOS rising, taking OXPHOS
# OUT OF THE DENOMINATOR must abolish it.
#
#   U_full     the canonical universe (the value E03 reports)
#   U_no_oxname  minus every pathway named for OXPHOS, a complex or a respirasome
#   U_no_oxgene  minus every pathway sharing ANY gene with `OXPHOS subunits`
#
# U_no_oxgene is the strict one: after it, no nuclear OXPHOS subunit appears
# anywhere in the denominator.
message("\n2. falsifier 3a: hold the mitoPPS denominator fixed")

OX_NAME_RX <- "OXPHOS|Complex [IV]+$|^C[IV]+ (subunits|assembly factors)$|Respirasome"

# nuc arrives ALREADY IN THE MATRIX'S VOCABULARY (arms_x is remapped upstream),
# so no symbol handling happens in here and neither cohort is a special case.
.universes <- function(paths_x, L, nuc) {
  S <- .path_scores(paths_x, L, MIN_SET_GENES)
  keep <- rownames(S)
  by_name <- setdiff(keep, grep(OX_NAME_RX, keep, value = TRUE))
  shares <- vapply(paths_x[keep], function(g) length(intersect(g, nuc)) > 0L,
                   logical(1))
  by_gene <- keep[!shares]
  list(full = S, no_oxname = S[by_name, , drop = FALSE],
       no_oxgene = S[by_gene, , drop = FALSE])
}

.requery <- function(paths_x, arms_x, L, label) {
  U <- .universes(paths_x, L, arms_x[["OXPHOS subunits"]])
  message("   ", label, " universe sizes: full ", nrow(U$full),
          " | no_oxname ", nrow(U$no_oxname), " | no_oxgene ", nrow(U$no_oxgene))
  S_arm <- .path_scores(arms_x[c("mtDNA-encoded OXPHOS", "OXPHOS subunits",
                                 "Mitochondrial ribosome")], L, MIN_SET_GENES)
  out <- list()
  for (u in names(U)) {
    Su <- U[[u]]
    for (arm in rownames(S_arm)) {
      hold <- arm_univ[[arm]]
      Sx <- if (!is.null(hold) && !is.na(hold) && hold %in% rownames(Su))
              Su[setdiff(rownames(Su), hold), , drop = FALSE] else Su
      out[[paste(u, arm)]] <- tibble::tibble(
        universe = u, n_paths = nrow(Sx), arm = arm,
        score = list(as.numeric(.mitopps_query(S_arm[arm, , drop = FALSE], Sx))))
    }
  }
  dplyr::bind_rows(out)
}

paths_t <- lapply(mito_paths, .in_t); paths_s <- lapply(mito_paths, .in_s)
arms_t  <- lapply(arm_sets,   .in_t); arms_s  <- lapply(arm_sets,   .in_s)
rq_t <- .requery(paths_t, arms_t, LT, "TCGA")
rq_s <- .requery(paths_s, arms_s, LS, "SCAN-B")

.rho_panel <- function(rq, EST, ids, cohort) {
  M <- do.call(rbind, lapply(seq_len(nrow(rq)), function(i) rq$score[[i]]))
  rownames(M) <- paste(rq$universe, rq$arm, sep = "|"); colnames(M) <- ids
  .atlas_block(EST, M, ids, NULL, min_n = MIN_N) %>%
    dplyr::mutate(cohort = cohort,
                  universe = sub("\\|.*$", "", measure),
                  arm = sub("^[^|]+\\|", "", measure)) %>%
    dplyr::select(cohort, universe, arm, myc_estimator, n, rho, ci_lo, ci_hi)
}
denominator_test <- dplyr::bind_rows(
  .rho_panel(rq_t, EST_T, ID_T, "TCGA"),
  .rho_panel(rq_s, EST_S, ID_S, "SCAN-B")) %>%
  dplyr::left_join(dplyr::select(a$est_meta, myc_estimator, frac_prolif),
                   by = "myc_estimator")

message("\n   FELSHER_61, mitoPPS rho, by what is left in the denominator:")
denominator_test %>%
  dplyr::filter(myc_estimator == "FELSHER_61") %>%
  tidyr::pivot_wider(id_cols = c(cohort, arm), names_from = universe,
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the mtDNA arm across all 20 estimators - how many stay negative:")
denominator_test %>%
  dplyr::filter(arm == "mtDNA-encoded OXPHOS") %>%
  dplyr::group_by(cohort, universe) %>%
  dplyr::summarise(n_est = dplyr::n(), n_negative = sum(rho < 0),
                   median_rho = round(stats::median(rho), 3),
                   min = round(min(rho), 3), max = round(max(rho), 3),
                   .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 3. FALSIFIER 3b - drop mitoPPS entirely and use a transparent ratio
# =============================================================================
# The strongest form of the test, because it needs no estimator at all: a log2
# ratio of mean linear expression that any reader can recompute.
#
#   vs_mito_rest   mtDNA 13  /  every other MitoCarta pathway gene
#   vs_mito_noox   mtDNA 13  /  MitoCarta pathway genes minus nuclear OXPHOS
#                               subunits and assembly factors  <- THE TEST
#   vs_mitoribo    mtDNA 13  /  mitochondrial ribosome  <- the F6 tension
#   vs_nuclear_ox  mtDNA 13  /  nuclear OXPHOS subunits
#   vs_all_genes   mtDNA 13  /  the whole matrix
#
# If `vs_mito_noox` is still negative, the composition-identity explanation is
# dead: there is no OXPHOS in that denominator to force it.
message("\n3. falsifier 3b: instrument-free log2 ratios")

.ratios <- function(L, inf) {
  invisible(.mt_matrix(L, inf))   # refuses here if any of the 13 is missing
  mt   <- inf(MT_GENES)           # the matrix's own names, for the setdiffs
  allm <- unique(unlist(lapply(mito_paths, inf), use.names = FALSE))
  nucox <- unique(c(inf(NUC_OX), inf(arm_sets[["OXPHOS assembly factors"]])))
  mrp  <- inf(arm_sets[["Mitochondrial ribosome"]])
  num  <- colMeans(L[mt, , drop = FALSE])
  .r <- function(den_genes) log2(num / colMeans(L[den_genes, , drop = FALSE]))
  rbind(vs_mito_rest  = .r(setdiff(allm, mt)),
        vs_mito_noox  = .r(setdiff(allm, c(mt, nucox))),
        vs_mitoribo   = .r(setdiff(mrp, mt)),
        vs_nuclear_ox = .r(setdiff(inf(NUC_OX), mt)),
        vs_all_genes  = .r(setdiff(rownames(L), mt)))
}
R_T <- .ratios(LT, .in_t); R_S <- .ratios(LS, .in_s)
stopifnot(all(is.finite(R_T)), all(is.finite(R_S)))
colnames(R_T) <- ID_T; colnames(R_S) <- ID_S

ratio_test <- dplyr::bind_rows(
  .atlas_block(EST_T, R_T, ID_T, NULL, min_n = MIN_N) %>%
    dplyr::mutate(cohort = "TCGA"),
  .atlas_block(EST_S, R_S, ID_S, NULL, min_n = MIN_N) %>%
    dplyr::mutate(cohort = "SCAN-B")) %>%
  dplyr::rename(ratio = measure) %>%
  dplyr::select(cohort, ratio, myc_estimator, n, rho, ci_lo, ci_hi)

message("\n   FELSHER_61 against each ratio:")
ratio_test %>% dplyr::filter(myc_estimator == "FELSHER_61") %>%
  tidyr::pivot_wider(id_cols = ratio, names_from = cohort, values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   all 20 estimators, how many negative per ratio:")
ratio_test %>% dplyr::group_by(cohort, ratio) %>%
  dplyr::summarise(n_est = dplyr::n(), n_negative = sum(rho < 0),
                   median_rho = round(stats::median(rho), 3), .groups = "drop") %>%
  tidyr::pivot_wider(id_cols = ratio, names_from = cohort,
                     values_from = c(n_negative, median_rho)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4. FALSIFIER 1 - mtDNA copy number
# =============================================================================
# If MYC-high tumours simply carry fewer mitochondrial genomes, the 13 genes
# should move TOGETHER. Content proxy = log2 mean linear expression of the 13.
# Putting it in as a covariate removes whatever the 13 share; what is left is
# WITHIN-mtDNA regulation. If the MT-CO2 / MT-ND5 split collapses, F3's per-gene
# claim was copy number plus quantification noise.
message("\n4. falsifier 1: mtDNA content as a covariate")

.gene_block <- function(L, inf, EST, ids, cohort) {
  Mg <- .mt_matrix(L, inf)[, ids, drop = FALSE]
  content <- matrix(log2(colMeans(Mg)), ncol = 1,
                    dimnames = list(ids, "mtDNA_content"))
  dplyr::bind_rows(
    .atlas_block(EST, Mg, ids, NULL,    min_n = MIN_N) %>%
      dplyr::mutate(adjusted = "raw"),
    .atlas_block(EST, Mg, ids, content, min_n = MIN_N) %>%
      dplyr::mutate(adjusted = "mtDNA_content")) %>%
    dplyr::mutate(cohort = cohort) %>%
    dplyr::rename(gene = measure) %>%
    dplyr::select(cohort, gene, myc_estimator, adjusted, n, rho, ci_lo, ci_hi)
}
copy_number_test <- dplyr::bind_rows(
  .gene_block(LT, .in_t, EST_T, ID_T, "TCGA"),
  .gene_block(LS, .in_s, EST_S, ID_S, "SCAN-B"))

message("\n   FELSHER_61, raw vs content-adjusted, both cohorts:")
copy_number_test %>% dplyr::filter(myc_estimator == "FELSHER_61") %>%
  tidyr::pivot_wider(id_cols = gene, names_from = c(cohort, adjusted),
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(TCGA_raw)) %>% as.data.frame() %>%
  print(row.names = FALSE)

spread <- copy_number_test %>%
  dplyr::group_by(cohort, myc_estimator, adjusted) %>%
  dplyr::summarise(spread = max(rho) - min(rho), .groups = "drop") %>%
  dplyr::group_by(cohort, adjusted) %>%
  dplyr::summarise(median_spread_over_13_genes = round(stats::median(spread), 3),
                   .groups = "drop")
message("\n   spread across the 13 genes (max rho - min rho), median over estimators:")
spread %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. FALSIFIER 2 - quantification
# =============================================================================
# MT-ND5 and MT-ND6 overlap on opposite strands. If multi-mapping manufactured
# the pattern, those two should be anomalous in the 13x13 correlation structure
# - unusually correlated with each other, or unusually decoupled from the rest.
message("\n5. falsifier 2: the 13 x 13 structure")

.gene_cor <- function(L, inf, ids) {
  stats::cor(t(.rank_rows(.mt_matrix(L, inf)[, ids, drop = FALSE])))
}
C_T <- .gene_cor(LT, .in_t, ID_T); C_S <- .gene_cor(LS, .in_s, ID_S)
.pair_rank <- function(C, g1, g2) {
  v <- C[lower.tri(C)]
  list(rho = C[g1, g2], percentile = mean(v <= C[g1, g2]))
}
qtab <- dplyr::bind_rows(lapply(list(TCGA = C_T, `SCAN-B` = C_S), function(C) {
  off <- C[lower.tri(C)]
  nd56 <- .pair_rank(C, "MT-ND5", "MT-ND6")
  tibble::tibble(
    median_pairwise = stats::median(off), min_pairwise = min(off),
    max_pairwise = max(off),
    ND5_ND6_rho = nd56$rho, ND5_ND6_percentile = nd56$percentile,
    ND5_mean_with_others = mean(C["MT-ND5", setdiff(colnames(C), "MT-ND5")]),
    ND6_mean_with_others = mean(C["MT-ND6", setdiff(colnames(C), "MT-ND6")]),
    CO2_mean_with_others = mean(C["MT-CO2", setdiff(colnames(C), "MT-CO2")]))
}), .id = "cohort")
qtab %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. Verdict
# =============================================================================
# Stated mechanically from the numbers, so the script cannot be read as arguing
# for the finding it was written to attack.
message("\n6. verdict")

.surv <- function(x) if (x) "SURVIVES" else "FALSIFIED"
v_noox <- ratio_test %>%
  dplyr::filter(ratio == "vs_mito_noox", myc_estimator == "FELSHER_61")
v_den <- denominator_test %>%
  dplyr::filter(arm == "mtDNA-encoded OXPHOS", universe == "no_oxgene",
                myc_estimator == "FELSHER_61")
v_cn <- copy_number_test %>%
  dplyr::filter(myc_estimator == "FELSHER_61", adjusted == "mtDNA_content",
                gene %in% c("MT-CO2", "MT-ND5")) %>%
  dplyr::select(cohort, gene, rho) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rho)

verdict <- tibble::tibble(
  falsifier = c("3a mitoPPS denominator (no OXPHOS gene in it)",
                "3b instrument-free ratio vs non-OXPHOS mitochondrion",
                "1 mtDNA copy number (content as covariate)"),
  criterion = c("mtDNA arm rho still < 0 in BOTH cohorts",
                "log2 ratio rho still < 0 in BOTH cohorts",
                "MT-CO2 still > 0 and MT-ND5 still < 0 in BOTH cohorts"),
  result = c(
    .surv(all(v_den$rho < 0) && dplyr::n_distinct(v_den$cohort) == 2L),
    .surv(all(v_noox$rho < 0) && dplyr::n_distinct(v_noox$cohort) == 2L),
    .surv(all(v_cn$`MT-CO2` > 0) && all(v_cn$`MT-ND5` < 0) &&
          nrow(v_cn) == 2L)))
verdict %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   Falsifier 2 (quantification) is not a pass/fail test - read the ",
        "13 x 13 table\n   in section 5 and judge whether ND5/ND6 look anomalous.")

# =============================================================================
# 7. Save
# =============================================================================
saveRDS(list(
  denominator_test = denominator_test,
  ratio_test       = ratio_test,
  copy_number_test = copy_number_test,
  gene_spread      = spread,
  gene_cor         = list(TCGA = C_T, `SCAN-B` = C_S),
  quantification   = qtab,
  verdict          = verdict,
  rules = list(
    purpose = paste("this script exists to KILL F3 if F3 is wrong; whatever it",
                    "returns is recorded beside the finding, not folded into it"),
    declared_in = "docs/2026-08-31_phase1_atlas_findings.md section F3",
    mitopps = paste("mitoPPS values are never compared across cohorts here -",
                    "only the SIGN and the pattern across denominators")),
  built = Sys.time()), PATH_FALSIFIERS)

message("\nE03b: done.")
message("    results/mtdna_falsifiers.rds")
message("    Record the verdict in the findings note BEFORE moving on.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  f <- readRDS(PATH_FALSIFIERS)
  f$verdict %>% as.data.frame()

  # the whole denominator series, every estimator
  f$denominator_test %>%
    dplyr::filter(arm == "mtDNA-encoded OXPHOS") %>%
    tidyr::pivot_wider(id_cols = c(myc_estimator, frac_prolif),
                       names_from = c(cohort, universe), values_from = rho) %>%
    dplyr::arrange(frac_prolif) %>% as.data.frame()

  # the ratios
  f$ratio_test %>%
    tidyr::pivot_wider(id_cols = myc_estimator, names_from = c(cohort, ratio),
                       values_from = rho) %>% as.data.frame()

  # per gene, what survives conditioning on mtDNA content
  f$copy_number_test %>%
    dplyr::filter(myc_estimator %in% c("FELSHER_61", "MYC_UP.V1_UP")) %>%
    tidyr::pivot_wider(id_cols = gene,
                       names_from = c(cohort, myc_estimator, adjusted),
                       values_from = rho) %>% as.data.frame()

  round(f$gene_cor$TCGA, 2)
  round(f$gene_cor$`SCAN-B`, 2)

}
