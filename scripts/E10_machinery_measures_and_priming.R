# E10_machinery_measures_and_priming.R
# =============================================================================
# PHASE 2. Three requests, all about the SAME 44 genes and the SAME two axes,
# so they share one script, one gene matrix and one set of measures.
#
#   Q-a (handoff open item 5.4) RE-ANNOTATE THE CANONICAL MACHINERY FROM AN
#       INDEPENDENT SOURCE. S6 says that among the 44 genes carrying a
#       `family_pathway` label the D1 direction reverses, and that the split
#       predicting the sign is WHERE THE GENE ACTS (rho 0.453) rather than its
#       annotated direction of effect (0.225). Both of those labels came from
#       the same CSV. If the finding is real it must survive being re-derived
#       from sources that never saw that curation: REACTOME for the module and
#       MITOCARTA 3.0 for the localisation. `APAF1` at -0.452 is the named
#       anomaly and this is the test that either explains or keeps it.
#
#   Q-b (author) THE SAME FIGURE AGAINST MYC. E08 fig6 is the 44 genes against
#       OXPHOS subunits. The MYC axis was computed in E08 but never plotted, and
#       the two axes are not interchangeable - F2 says MYC mRNA and MYC activity
#       disagree completely, and D3/S1 say the BCL2 column against MYC is
#       between-subtype pooling while the OXPHOS column is stratum-stable.
#
#   Q-c (author) THE BCL2-FAMILY PRIMING RATIOS - each of 7 pro-apoptotic genes
#       over each of 5 anti-apoptotic ones - against MYC and OXPHOS, as a
#       heatmap, and inside the luminal and basal compartments.
#
# =============================================================================
# EVERYTHING HERE IS SPEARMAN, AND E09 IS WHY
# =============================================================================
# An earlier version of this script carried a Pearson panel beside every
# Spearman one. E09 ran on 2026-09-02 and made that redundant: over 220 pairs
# the two measures correlate at 0.996, bicor sits WITH Spearman in all twelve of
# the largest disagreements - so every one of them is a heavy tail rather than
# magnitude carrying information ranks discard - and the mean departure is 0.009
# on GSVA. Nothing in the atlas is non-monotone either; the largest spline gain
# over a straight line is 0.051 and the headline pair is 0.012.
#
# So the Pearson panels were four extra figures that said what the Spearman
# panels said. THEY WERE REMOVED ON 2026-09-02, and the finding that they showed
# nothing is E09's, not this script's - see
# docs/2026-09-02_e09_correlation_measures.md, C1 to C4.
#
# THE ONE PLACE E09 SAYS THE MEASURE DOES MATTER IS mitoPPS, where the mean
# departure is 0.029 and reaches 0.093. A Pearson against a mitoPPS score is not
# to be reported anywhere in this study.
#
# SCALE. The gene matrix is log2(linear DESeq2-normalised + 1). Now that every
# correlation is rank-based this is invariant - Spearman on log2(x + 1) is the
# same number as Spearman on x, which is what lets section 4 assert against
# E08's values computed on the raw linear matrix. The transform survives because
# `expr_pct` is taken on it, and an expression percentile on a log scale is not
# dominated by a handful of enormous genes the way a linear one is.
#
# =============================================================================
# THE PRIMING RATIOS, AND THE ONE THING TO KNOW BEFORE READING THEM
# =============================================================================
# The author's lists:
#   pro   BCL2L11 (BIM), BMF, PMAIP1 (NOXA), BBC3 (PUMA), BID, BAD, BIK
#   anti  BCL2, BCL2L1 (Bcl-XL), MCL1, BCL2L2 (Bcl-w), BCL2A1 (A1)
#
# 7 x 5 = 35 ratios. BCL2L2 was originally supplied on both sides; the first run
# flagged it as canonically anti-apoptotic and the author confirmed on
# 2026-09-02 that the pro-side entry was an error. It is now anti only, so there
# is no degenerate self-pair and no empty cell in the heatmap.
#
# A RATIO IS ONLY INTERESTING IF IT BEATS ITS PARTS, and that is this section's
# built-in falsifier. log2(pro) - log2(anti) is a difference of two correlated
# variables; when the two genes co-express, the difference is mostly noise, and
# when they do not, the ratio's rho is close to the stronger component's. Every
# ratio is therefore reported beside `rho_pro`, `rho_anti`, the pro-anti
# co-expression, and `gain` = |rho_ratio| - max(|rho_pro|, |rho_anti|).
# GAIN AT OR BELOW ZERO MEANS THE RATIO ADDED NOTHING and the single gene is the
# better measurement.
#
# =============================================================================
# WHAT WOULD FALSIFY EACH CLAIM, WRITTEN BEFORE THE ANSWERS ARE SEEN
# =============================================================================
#   S6 SURVIVES if MitoCarta membership - a localisation call made with no
#   knowledge of apoptosis - predicts the sign about as well as the curation's
#   own `is_mitochondrial` did (0.453), AND better than direction of effect
#   (0.225). It FAILS if the independent localisation is uninformative, which
#   would mean the original 0.453 was an artefact of one curator's habit.
#
#   THE MODULE STORY SURVIVES if Reactome's own intrinsic/extrinsic partition
#   reproduces the sign ordering. It FAILS if Reactome scrambles it - and note
#   in advance that these two annotations answer DIFFERENT questions: a pathway
#   membership is about what a gene DOES, a localisation is about where it IS,
#   and APAF1 is exactly the gene where they come apart.
#
#   THE PRIMING RATIOS ARE WORTH REPORTING only where `gain` is positive in BOTH
#   cohorts. Anything else is a single gene wearing a ratio's name.
#
# EXPLORATORY. Nothing here is pre-registered. The gene lists in Q-c were chosen
# by the author on biological grounds BEFORE any of these numbers were seen,
# which is the one thing that makes 39 ratios a panel rather than a fishing
# expedition - but it is not a pre-registration and must not be read as one.
#
# SCALE: log2(linear + 1) at gene level; GSVA as built at axis level. Every
# correlation is Spearman. SPECIES: human, natively. No ortholog function is
# called.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))
source(here::here("functions", "strata.R"))

message("\nE10: the canonical machinery on both axes, and the priming ratios\n",
        strrep("=", 78))

PATH_E10       <- file.path(DIR_RESULTS, "machinery_and_priming.rds")
PATH_E10_CANON <- file.path(DIR_TABLES,  "E10_canonical_machinery.csv")
PATH_E10_PRIME <- file.path(DIR_TABLES,  "E10_priming_ratios.csv")

# The author's lists. BCL2L2 was originally given on BOTH sides and was
# confirmed on 2026-09-02 to be an error on the pro side - it is Bcl-w, a
# multidomain guardian, not a BH3-only sensitiser. It is anti-apoptotic here and
# nowhere else.
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")

# The priming ratios are read in the luminal and basal compartments as well as
# pooled. The mouse arm shows luminal expansion, so `Luminal` = LumA + LumB is
# the compartment the question is about; `Basal` is its contrast. `all` stays in
# so the pooled value is beside them and a between-subtype effect - which is
# what D3/S1 turned out to be for BCL2 against MYC - stays visible as the
# pooled column disagreeing with both.
STRATA_PRIMING <- c("all", "Luminal", "Basal")
MIN_STRATUM_N  <- 30L

# Percentile below which a gene-level correlation is flagged as unreadable.
# Same value E08 used, and for the same reason: below it a Spearman is largely
# a correlation of quantisation noise with a score.
LOW_EXPR_PCT <- 0.25

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

e08  <- readRDS(file.path(DIR_RESULTS, "strata_and_death_genes.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)

# DEPENDENCY GUARD. Handoff trap 1: a saved artefact that predates a rename
# fails three sections downstream as something unrecognisable. Check here.
if (!"canonical" %in% names(e08) || nrow(e08$canonical) != 44L) {
  stop("results/strata_and_death_genes.rds does not carry the 44-gene ",
       "`canonical` table this script re-annotates. Re-run E08.", call. = FALSE)
}
if (!all(c(MYC_REF, MYC_LOW_ENTANG) %in% rownames(nw$tcga_gsva_new))) {
  stop("the MYC estimators named by E00's naming contract are not in ",
       "new_set_scores.rds. Re-run E02.", call. = FALSE)
}

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)

tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
# log2(x + 1) is monotone, so every Spearman below is identical to the one E08
# computed on the raw linear matrix - section 4 asserts it. The transform is
# kept because `expr_pct` is taken on it, and an expression percentile on a log
# scale is not dominated by a handful of enormous genes the way a linear one is.
GT <- log2(tcga_lin$mat[, ID_T, drop = FALSE] + 1)
GS <- log2(scanb_lin$mat[, ID_S, drop = FALSE] + 1)
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

RES_T <- .symbol_resolver(rownames(GT), NULL)
RES_S <- .symbol_resolver(rownames(GS), sc$symbol_map)

# The axes. Identical construction to E08 section 1 so the two scripts are
# reading the same two columns of the same plane.
.axes <- function(gsva_new, arms_obj, mb, ids) {
  m <- rbind(
    MYC            = as.numeric(gsva_new[MYC_REF, ids]),
    MYC_low_entang = as.numeric(gsva_new[MYC_LOW_ENTANG, ids]),
    M_b_ref        = as.numeric(mb[ids]),
    OXPHOS         = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]),
    OXPHOS_mitopps = as.numeric(arms_obj$mitopps_arms["OXPHOS subunits", ids]),
    MITORIBO       = as.numeric(arms_obj$gsva_arms["Mitochondrial ribosome", ids]))
  colnames(m) <- ids
  m
}
AX_T <- .axes(nw$tcga_gsva_new, mito, nw$tcga_M_b_variants[MB_REF, ID_T], ID_T)
AX_S <- .axes(sc$gsva_new,      sc,   sc$M_b_variants[MB_REF, ],          ID_S)

# PROLIF_DISJOINT, the same covariate E11 and E13 use. It is disjoint from
# FELSHER__MITOSTRIP and from the OXPHOS-subunits arm, which is what lets the
# same adjustment be applied to both axes - see E11's header.
PROLIF_COV <- "PROLIF_DISJOINT"
CV_T <- t(mito$gsva_cov[PROLIF_COV, ID_T, drop = FALSE])
CV_S <- t(sc$gsva_cov[PROLIF_COV, ID_S, drop = FALSE])
stopifnot(identical(rownames(CV_T), ID_T), identical(rownames(CV_S), ID_S))

COH <- list(
  TCGA     = list(G = GT, res = RES_T, ax = AX_T, ids = ID_T, cov = CV_T),
  `SCAN-B` = list(G = GS, res = RES_S, ax = AX_S, ids = ID_S, cov = CV_S))
message("   TCGA ", length(ID_T), " samples | SCAN-B ", length(ID_S),
        " samples | ", nrow(AX_T), " axes")

frames <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
STR <- list(TCGA = .build_strata(frames, "TCGA", ID_T),
            `SCAN-B` = .build_strata(frames, "SCAN-B", ID_S))
tibble::tibble(stratum = STRATA_PRIMING,
               TCGA = lengths(STR$TCGA[STRATA_PRIMING]),
               `SCAN-B` = lengths(STR$`SCAN-B`[STRATA_PRIMING])) %>%
  as.data.frame() %>% print(row.names = FALSE)

CANON_GENES  <- sort(e08$canonical$gene)
PRIMING_ALL  <- sort(unique(c(PRIMING_PRO, PRIMING_ANTI)))
GENES_WANTED <- sort(unique(c(CANON_GENES, PRIMING_ALL)))
message("   ", length(CANON_GENES), " canonical machinery genes + ",
        length(PRIMING_ALL), " BCL2-family priming genes = ",
        length(GENES_WANTED), " distinct")

GR <- lapply(COH, function(C) .gene_rows(GENES_WANTED, C$G, C$res))
for (coh in names(GR)) {
  if (length(GR[[coh]]$missing)) {
    message("   MISSING in ", coh, ": ",
            paste(GR[[coh]]$missing, collapse = ", "))
  } else {
    message("   all ", length(GENES_WANTED), " genes present in ", coh)
  }
}
# Every one of these genes was present in both matrices when E08 ran. If that
# has changed, the sets below would silently shrink rather than fail, so stop
# here and name what went missing instead of reporting a smaller panel.
miss_all <- unique(unlist(lapply(GR, function(g) g$missing), use.names = FALSE))
if (length(miss_all)) {
  stop("these genes did not resolve to exactly one matrix row: ",
       paste(miss_all, collapse = ", "),
       "\nCheck the SCAN-B symbol map before going further.", call. = FALSE)
}

# =============================================================================
# 2. Measures, and the expression floor
# =============================================================================
message("\n2. measures")

# Spearman of every AXIS against every ITEM in one call, optionally partialled
# on a covariate. A partial Spearman is Pearson on the residuals of the ranks,
# which is what functions/correlation_engine.R does and what E11 and E13 do.
#
# EVERYTHING PLOTTED BY THIS SCRIPT IS ADJUSTED FOR PROLIFERATION, added
# 2026-09-02 so that its panels sit on the same footing as E11's and E13's. Two
# figures in one paper cannot be on different footings. The UNADJUSTED values
# are still computed for section 4's reproduction check against E08, which was
# built without a covariate, and are carried in the saved table.
#
# Fisher-z intervals with the Bonett-Wright variance, se =
# sqrt((1 + rho^2/2)/(n - 3 - k)), k covariates projected out; the plain
# 1/(n-3) is the Pearson case and understates a rank correlation.
#
# THE INTERVAL IS FOR DESCRIBING ONE CELL, NEVER FOR SELECTING ONE. This script
# emits 44 x 6 x 2 x 2 gene cells plus 39 x 6 x 2 x 2 ratio cells; a cell whose
# interval excludes zero is not a finding.
.rank_rows <- function(M) {
  out <- t(apply(M, 1L, rank))
  if (nrow(M) == 1L) out <- matrix(out, nrow = 1L)
  dimnames(out) <- dimnames(M)
  out
}
.cor_block <- function(A, B, cov = NULL) {
  n <- ncol(A)
  stopifnot(identical(colnames(A), colnames(B)))
  RA <- .rank_rows(A); RB <- .rank_rows(B)
  k  <- 0L
  if (!is.null(cov)) {
    stopifnot(nrow(cov) == n)
    H  <- qr(cbind(`(Intercept)` = 1, apply(cov, 2L, rank)))
    k  <- ncol(cov)
    RA <- RA - t(qr.fitted(H, t(RA)))
    RB <- RB - t(qr.fitted(H, t(RB)))
  }
  R  <- suppressWarnings(stats::cor(t(RA), t(RB)))
  z  <- atanh(pmin(pmax(R, -0.999999999), 0.999999999))
  se <- sqrt((1 + R^2 / 2) / (n - 3 - k))
  tibble::tibble(
    axis   = rep(rownames(R), times = ncol(R)),
    item   = rep(colnames(R), each  = nrow(R)),
    n      = n,
    k_cov  = k,
    rho    = as.vector(R),
    ci_lo  = as.vector(tanh(z - 1.959964 * se)),
    ci_hi  = as.vector(tanh(z + 1.959964 * se)))
}

# Expression percentile, per cohort, over the whole matrix. Carried so that a
# rho on a barely-expressed gene is visible as such. FLAG, NEVER FILTER.
expr_rank <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  e <- rowMeans(C$G)
  pct <- rank(e) / length(e)
  g <- rownames(GR[[coh]]$mat)
  tibble::tibble(cohort = coh, gene = g,
                 expr_log2mean = e[g], expr_pct = pct[g]) %>%
    dplyr::mutate(expr_decile = ceiling(10 * expr_pct),
                  low_expression = expr_pct < LOW_EXPR_PCT)
}))
message("   genes below the ", LOW_EXPR_PCT * 100, "th expression percentile:")
low <- expr_rank %>% dplyr::filter(low_expression)
if (nrow(low)) {
  low %>% dplyr::mutate(expr_pct = round(expr_pct, 3)) %>%
    dplyr::select(cohort, gene, expr_pct) %>% dplyr::arrange(cohort, expr_pct) %>%
    as.data.frame() %>% print(row.names = FALSE)
} else {
  message("   none - every gene in this script clears the floor in both cohorts")
}
message("\n   the 12 priming genes, by expression decile:")
expr_rank %>% dplyr::filter(gene %in% PRIMING_ALL) %>%
  tidyr::pivot_wider(id_cols = gene, names_from = cohort,
                     values_from = expr_decile) %>%
  dplyr::arrange(gene) %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 3. Q-a: re-annotating the 44 from sources that never saw the curation
# =============================================================================
# S6's two predictors both came out of cell_death_genes_consolidated.csv:
# `family_pathway` gave the module and `is_mitochondrial` gave the localisation.
# A single curator's habits could produce both. Two independent replacements:
#
#   MODULE        MSigDB C2:CP:REACTOME. Reactome partitions apoptosis into
#                 reaction-level sub-pathways - intrinsic, extrinsic, execution,
#                 NF-kB, p53-transcriptional - built from primary literature
#                 with no reference to this study's curation.
#   LOCALISATION  Human MitoCarta 3.0, already pinned in this repo. Its calls
#                 come from mass spectrometry, GFP imaging and targeting-signal
#                 prediction. It has never heard of apoptosis.
#
# THE TWO ANSWER DIFFERENT QUESTIONS and that is the point, not a nuisance:
# a Reactome pathway says what a gene DOES, MitoCarta says where it IS. APAF1
# is the gene where they come apart - apoptosome, therefore intrinsic by
# pathway, but cytosolic, therefore not mitochondrial by localisation.
message("\n3. Q-a: independent re-annotation of the 44")

reactome <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                             subcollection = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol, db_version) %>%
  dplyr::distinct()
MSIGDB_VERSION <- unique(reactome$db_version)
message("   MSigDB ", paste(MSIGDB_VERSION, collapse = "/"), " C2:CP:REACTOME - ",
        dplyr::n_distinct(reactome$gs_name), " sets")

# The anchor pathways, named one by one rather than pattern-matched, because a
# regex over set names is a rule nobody can check and it silently changes
# meaning when MSigDB adds a set.
REACTOME_MODULES <- list(
  `intrinsic / mitochondrial` = c(
    "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
    "REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS",
    "REACTOME_BH3_ONLY_PROTEINS_ASSOCIATE_WITH_AND_INACTIVATE_ANTI_APOPTOTIC_BCL_2_MEMBERS",
    "REACTOME_RELEASE_OF_APOPTOTIC_FACTORS_FROM_THE_MITOCHONDRIA",
    "REACTOME_CYTOCHROME_C_MEDIATED_APOPTOTIC_RESPONSE",
    "REACTOME_FORMATION_OF_APOPTOSOME",
    "REACTOME_ACTIVATION_OF_CASPASES_THROUGH_APOPTOSOME_MEDIATED_CLEAVAGE",
    "REACTOME_SMAC_XIAP_REGULATED_APOPTOTIC_RESPONSE",
    "REACTOME_APOPTOTIC_FACTOR_MEDIATED_RESPONSE"),
  `extrinsic / death receptor` = c(
    "REACTOME_DEATH_RECEPTOR_SIGNALING",
    "REACTOME_CASPASE_ACTIVATION_VIA_DEATH_RECEPTORS_IN_THE_PRESENCE_OF_LIGAND",
    "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
    "REACTOME_CASPASE_ACTIVATION_VIA_DEPENDENCE_RECEPTORS_IN_THE_ABSENCE_OF_LIGAND",
    "REACTOME_TNFR1_INDUCED_PROAPOPTOTIC_SIGNALING",
    "REACTOME_TNF_SIGNALING",
    "REACTOME_TNFS_BIND_THEIR_PHYSIOLOGICAL_RECEPTORS"),
  `execution phase` = c(
    "REACTOME_APOPTOTIC_EXECUTION_PHASE",
    "REACTOME_APOPTOTIC_CLEAVAGE_OF_CELLULAR_PROTEINS",
    "REACTOME_CASPASE_MEDIATED_CLEAVAGE_OF_CYTOSKELETAL_PROTEINS",
    "REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION",
    "REACTOME_APOPTOTIC_CLEAVAGE_OF_CELL_ADHESION_PROTEINS"),
  `NF-kB / survival` = c(
    "REACTOME_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY",
    "REACTOME_NF_KB_ACTIVATION_THROUGH_FADD_RIP_1_PATHWAY_MEDIATED_BY_CASPASE_8_AND_10",
    "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",
    "REACTOME_TNF_RECEPTOR_SUPERFAMILY_TNFSF_MEMBERS_MEDIATING_NON_CANONICAL_NF_KB_PATHWAY",
    "REACTOME_SUPPRESSION_OF_APOPTOSIS",
    "REACTOME_REGULATION_OF_APOPTOSIS"),
  `p53 transcriptional` = c(
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CASPASE_ACTIVATORS_AND_CASPASES",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DEATH_RECEPTORS_AND_LIGANDS",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_CYTOCHROME_C_RELEASE",
    paste0("REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_SEVERAL_ADDITIONAL_CELL_",
           "DEATH_GENES_WHOSE_SPECIFIC_ROLES_IN_P53_DEPENDENT_APOPTOSIS_",
           "REMAIN_UNCERTAIN")))

# FAIL LOUDLY IF MSIGDB MOVED. A renamed or retired set would otherwise just
# make a module smaller and shift a median with no message anywhere.
missing_sets <- setdiff(unlist(REACTOME_MODULES, use.names = FALSE),
                        unique(reactome$gs_name))
if (length(missing_sets)) {
  stop("these Reactome sets are not in the installed MSigDB (",
       paste(MSIGDB_VERSION, collapse = "/"), "): ",
       paste(missing_sets, collapse = ", "),
       "\nThe module definitions must be revisited, not patched around.",
       call. = FALSE)
}
message("   ", length(unlist(REACTOME_MODULES)), " anchor pathways across ",
        length(REACTOME_MODULES), " modules, all present")

# Hit counts per gene per module. A gene in several modules is normal - the
# Reactome hierarchy nests - so the assignment takes the module with most hits
# and breaks ties by a STATED priority, and the number of genes decided by the
# tie-break is printed rather than hidden.
REACTOME_PRIORITY <- names(REACTOME_MODULES)
hit_counts <- vapply(REACTOME_MODULES, function(sets)
  vapply(CANON_GENES, function(g)
    sum(reactome$gene_symbol == g & reactome$gs_name %in% sets), integer(1)),
  integer(length(CANON_GENES)))
rownames(hit_counts) <- CANON_GENES

reactome_annot <- tibble::tibble(
  gene = CANON_GENES,
  reactome_hits = rowSums(hit_counts),
  reactome_module = apply(hit_counts, 1L, function(v) {
    if (all(v == 0)) return("not in a Reactome apoptosis pathway")
    REACTOME_PRIORITY[which(v == max(v))][1]
  }),
  reactome_tied = apply(hit_counts, 1L, function(v)
    sum(v == max(v)) > 1L && max(v) > 0),
  reactome_intrinsic = unname(hit_counts[, "intrinsic / mitochondrial"] > 0),
  reactome_extrinsic = unname(hit_counts[, "extrinsic / death receptor"] > 0))
message("   ", sum(reactome_annot$reactome_hits == 0), " of 44 are in NO ",
        "Reactome apoptosis pathway; ", sum(reactome_annot$reactome_tied),
        " were assigned by the priority tie-break")
if (any(reactome_annot$reactome_hits == 0)) {
  message("   absent from Reactome apoptosis: ",
          paste(reactome_annot$gene[reactome_annot$reactome_hits == 0],
                collapse = ", "))
}

# MitoCarta 3.0: membership, and the sub-mitochondrial compartment where it has
# one. Read from the pinned xls, not from a derived list, so the compartment
# column travels with the membership call.
# NOTE the object name. tibble() evaluates its arguments in sequence and each
# new column shadows anything of the same name, so a column called `mitocarta`
# built beside a data frame called `mitocarta` makes the second reference read
# the logical vector rather than the sheet. The column name is the one used
# downstream, so the sheet is what gets renamed.
mitocarta_sheet <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
stopifnot(all(c("Symbol", "MitoCarta3.0_SubMitoLocalization") %in%
                colnames(mitocarta_sheet)))
mc_row <- match(CANON_GENES, mitocarta_sheet$Symbol)
mito_annot <- tibble::tibble(
  gene = CANON_GENES,
  mitocarta = !is.na(mc_row),
  submito = mitocarta_sheet[["MitoCarta3.0_SubMitoLocalization"]][mc_row]) %>%
  dplyr::mutate(submito = dplyr::if_else(is.na(submito),
                                         "(not in MitoCarta)", submito))

# The annotated table: the curation's own labels beside both replacements.
reannot <- e08$canonical %>%
  dplyr::select(gene, TCGA, `SCAN-B`, mean_rho, effect,
                cdc_module = module, cdc_acts_at_mito = acts_at_mito,
                family_pathway) %>%
  dplyr::left_join(reactome_annot, by = "gene") %>%
  dplyr::left_join(mito_annot, by = "gene")

message("\n   does the INDEPENDENT localisation agree with the curation's?")
print(table(CDC_is_mitochondrial = reannot$cdc_acts_at_mito,
            MitoCarta_3.0 = reannot$mitocarta))
disagree <- reannot$gene[reannot$cdc_acts_at_mito != reannot$mitocarta]
message("   disagreements: ",
        if (length(disagree)) paste(disagree, collapse = ", ") else "none")

message("\n   module medians of the OXPHOS rho - the curation's modules:")
reannot %>% dplyr::group_by(cdc_module) %>%
  dplyr::summarise(n = dplyr::n(), n_negative = sum(mean_rho < 0),
                   median = round(stats::median(mean_rho), 3), .groups = "drop") %>%
  dplyr::arrange(median) %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   the same, from REACTOME:")
reannot %>% dplyr::group_by(reactome_module) %>%
  dplyr::summarise(n = dplyr::n(), n_negative = sum(mean_rho < 0),
                   median = round(stats::median(mean_rho), 3), .groups = "drop") %>%
  dplyr::arrange(median) %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   and by MITOCARTA sub-compartment (MIM/IMS are small - read the n):")
reannot %>% dplyr::group_by(submito) %>%
  dplyr::summarise(n = dplyr::n(), n_negative = sum(mean_rho < 0),
                   median = round(stats::median(mean_rho), 3), .groups = "drop") %>%
  dplyr::arrange(median) %>% as.data.frame() %>% print(row.names = FALSE)

# THE TEST S6 EITHER PASSES OR FAILS.
predictors <- tibble::tibble(
  predictor = c("acts at the mitochondrion (CDC curation)",
                "in MitoCarta 3.0 (INDEPENDENT localisation)",
                "in a Reactome INTRINSIC pathway (INDEPENDENT module)",
                "in a Reactome EXTRINSIC pathway (INDEPENDENT module)",
                "annotated pro-death (CDC curation)"),
  source = c("the curation under test", "independent", "independent",
             "independent", "the curation under test"),
  spearman_with_rho = c(
    stats::cor(reannot$mean_rho, as.numeric(reannot$cdc_acts_at_mito),
               method = "spearman"),
    stats::cor(reannot$mean_rho, as.numeric(reannot$mitocarta),
               method = "spearman"),
    stats::cor(reannot$mean_rho, as.numeric(reannot$reactome_intrinsic),
               method = "spearman"),
    stats::cor(reannot$mean_rho, as.numeric(reannot$reactome_extrinsic),
               method = "spearman"),
    stats::cor(reannot$mean_rho, as.numeric(reannot$effect == "pro-death"),
               method = "spearman")))
message("\n   WHICH SPLIT PREDICTS THE SIGN OF THE OXPHOS CORRELATION?")
predictors %>% dplyr::mutate(spearman_with_rho = round(spearman_with_rho, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   S6 SURVIVES if the INDEPENDENT localisation row is close to the",
        " curation's\n   own (0.453) and both beat direction of effect (0.225).")

message("\n   APAF1, the named anomaly, under each annotation:")
reannot %>% dplyr::filter(gene %in% c("APAF1", "AIFM1", "CYCS", "BCL2A1")) %>%
  dplyr::select(gene, mean_rho, effect, cdc_module, cdc_acts_at_mito,
                reactome_module, mitocarta, submito) %>%
  dplyr::mutate(mean_rho = round(mean_rho, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4. Q-b: the 44 genes on both axes
# =============================================================================
message("\n4. the canonical machinery on both axes")

gene_cor_raw <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  .cor_block(C$ax, GR[[coh]]$mat[CANON_GENES, , drop = FALSE]) %>%
    dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(gene = item)

gene_cor <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  .cor_block(C$ax, GR[[coh]]$mat[CANON_GENES, , drop = FALSE], cov = C$cov) %>%
    dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(gene = item)

# --- THE REPRODUCTION CHECK --------------------------------------------------
# E08 computed these Spearman values on the RAW LINEAR matrix through the atlas
# engine, with no covariate. This script computes them on log2(linear + 1)
# through stats::cor. log2(x + 1) is monotone, so the two must agree to floating
# point. THE CHECK USES `gene_cor_raw` - the unadjusted pass - because E08 has
# no proliferation adjustment to reproduce. Everything plotted below uses the
# adjusted pass. If this check fails, one of the two scripts is not reading the
# plane it says it is and every number below is suspect.
check <- gene_cor_raw %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::select(cohort, gene, rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
  dplyr::inner_join(dplyr::select(e08$canonical, gene, e08_TCGA = TCGA,
                                  `e08_SCAN-B` = `SCAN-B`), by = "gene")
max_drift <- max(abs(c(check$TCGA - check$e08_TCGA,
                       check$`SCAN-B` - check$`e08_SCAN-B`)))
message("   reproduction of E08's OXPHOS Spearman over ", nrow(check),
        " genes: max |difference| = ", format(max_drift, scientific = TRUE,
                                              digits = 3))
if (!isTRUE(max_drift < 1e-6)) {
  stop("this script does not reproduce E08's per-gene Spearman (max drift ",
       max_drift, "). The gene matrix, the symbol map or the axis definition ",
       "has moved between the two scripts. Do not read anything below.",
       call. = FALSE)
}
message("   OK - the two scripts are reading the same plane")

# --- FLAGS THAT DECIDE WHETHER A CELL IS INDEPENDENT -------------------------
# Two ways a gene-versus-score correlation can be partly a correlation with
# itself, and both are present here:
#   CYCS is a member of the `OXPHOS subunits` arm (1 of 89), so its OXPHOS
#   column is 1/89 self-overlap. E08 fig6 did not mark it.
#   E4 established that some MYC signatures contain death genes. The reference
#   estimator FELSHER__MITOSTRIP contains NONE of these 44 and none of the 12
#   priming genes, which is why it is the MYC axis plotted; M_b__MITOSTRIP
#   contains seven of them and is tabled, not plotted.
OX_ARM   <- sd_$arm_sets[["OXPHOS subunits"]]
in_myc   <- function(est_name, genes) {
  s <- if (est_name %in% names(sd_$myc_sets)) sd_$myc_sets[[est_name]] else
       sd_$collectri_sets[[est_name]]
  if (is.null(s)) {
    stop("no gene list for estimator ", est_name, " in set_definitions.rds",
         call. = FALSE)
  }
  intersect(genes, s)
}
message("\n   self-overlap audit:")
message("     in the OXPHOS subunits arm: ",
        paste(intersect(GENES_WANTED, OX_ARM), collapse = ", "))
for (est in c(MYC_REF, MYC_LOW_ENTANG, MB_REF)) {
  ov <- in_myc(est, GENES_WANTED)
  message("     in ", est, ": ",
          if (length(ov)) paste(ov, collapse = ", ") else "none")
}

MB_OVERLAP <- in_myc(MB_REF, GENES_WANTED)

canon_wide <- gene_cor %>%
  dplyr::select(cohort, axis, gene, rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
  dplyr::mutate(
    mean_rho = (TCGA + `SCAN-B`) / 2,
    agree    = sign(TCGA) == sign(`SCAN-B`),
    in_oxphos_arm  = gene %in% OX_ARM,
    in_m_b         = gene %in% MB_OVERLAP) %>%
  dplyr::left_join(dplyr::select(reannot, gene, effect, cdc_module,
                                 cdc_acts_at_mito, reactome_module, mitocarta,
                                 submito, family_pathway), by = "gene") %>%
  # A rho on a gene near the bottom of the expression range is largely a
  # correlation of quantisation noise with a score. FLAGGED, NEVER FILTERED -
  # E08's rule, and the flag has to reach the figures or it does no work.
  dplyr::left_join(
    expr_rank %>% dplyr::group_by(gene) %>%
      dplyr::summarise(low_expr_any = any(low_expression), .groups = "drop"),
    by = "gene")

message("\n   the 44 against MYC (", MYC_REF, "), Spearman, ranked:")
canon_wide %>%
  dplyr::filter(axis == "MYC") %>%
  dplyr::arrange(mean_rho) %>%
  dplyr::mutate(dplyr::across(c(TCGA, `SCAN-B`, mean_rho), ~ round(.x, 3))) %>%
  dplyr::select(gene, TCGA, `SCAN-B`, mean_rho, agree, effect,
                cdc_acts_at_mito, cdc_module) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   do the two axes agree about these 44 genes?")
axis_agree <- canon_wide %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(axis, gene, mean_rho) %>%
  tidyr::pivot_wider(names_from = axis, values_from = mean_rho) %>%
  dplyr::summarise(spearman_of_the_two_axes =
                     stats::cor(MYC, OXPHOS, method = "spearman"),
                   median_MYC = stats::median(MYC),
                   median_OXPHOS = stats::median(OXPHOS))
axis_agree %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   does the S6 split hold on the MYC axis too?")
s6_by_axis <- canon_wide %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::group_by(axis) %>%
  dplyr::summarise(
    vs_mitocarta = stats::cor(mean_rho, as.numeric(mitocarta),
                              method = "spearman"),
    vs_pro_death = stats::cor(mean_rho, as.numeric(effect == "pro-death"),
                              method = "spearman"),
    median_mito  = stats::median(mean_rho[mitocarta]),
    median_nonmito = stats::median(mean_rho[!mitocarta]), .groups = "drop")
s6_by_axis %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Q-c: the BCL2-family priming ratios
# =============================================================================
# On log2 expression a ratio is a DIFFERENCE, so log2(pro + 1) - log2(anti + 1)
# is the log priming ratio and needs no separate arithmetic. Every pro x anti
# combination is formed and none is dropped.
#
# THE RATIO IS NOT AUTOMATICALLY A BETTER MEASUREMENT THAN ITS PARTS. Two genes
# that co-express strongly give a difference that is mostly the residual noise
# of both; two that do not give a difference whose correlation is close to the
# stronger component's. `gain` is the whole point of the table.
message("\n5. the BCL2-family priming ratios")

message("   pro  (", length(PRIMING_PRO), "): ",
        paste(PRIMING_PRO, collapse = ", "))
message("   anti (", length(PRIMING_ANTI), "): ",
        paste(PRIMING_ANTI, collapse = ", "))
# The lists are disjoint by construction now, but an overlap would silently
# produce a gene-over-itself column of zeros, so it is checked rather than
# assumed.
if (length(intersect(PRIMING_PRO, PRIMING_ANTI))) {
  stop("these genes are on BOTH priming lists: ",
       paste(intersect(PRIMING_PRO, PRIMING_ANTI), collapse = ", "),
       ". A gene over itself is the constant zero.", call. = FALSE)
}

RATIO_GRID <- tidyr::expand_grid(pro = PRIMING_PRO, anti = PRIMING_ANTI) %>%
  dplyr::mutate(ratio = paste0(pro, "/", anti))
message("   ", nrow(RATIO_GRID), " ratios (", length(PRIMING_PRO), " x ",
        length(PRIMING_ANTI), ")")

.ratio_matrix <- function(M) {
  R <- M[RATIO_GRID$pro, , drop = FALSE] - M[RATIO_GRID$anti, , drop = FALSE]
  rownames(R) <- RATIO_GRID$ratio
  R
}
ratio_cor <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  M <- GR[[coh]]$mat
  .cor_block(C$ax, .ratio_matrix(M), cov = C$cov) %>% dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(ratio = item)

# The components, and the co-expression that decides whether a ratio can add
# anything at all.
component_cor <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  .cor_block(C$ax, GR[[coh]]$mat[PRIMING_ALL, , drop = FALSE], cov = C$cov) %>%
    dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(gene = item) %>%
  dplyr::left_join(dplyr::select(expr_rank, cohort, gene, expr_decile,
                                 low_expression), by = c("cohort", "gene")) %>%
  dplyr::mutate(side = dplyr::if_else(gene %in% PRIMING_PRO,
                                      "pro-apoptotic", "anti-apoptotic"))

# Adjusted too, so `coexpr` describes the same residual space as the rho it is
# read beside.
coexpr <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M <- GR[[coh]]$mat
  v <- vapply(seq_len(nrow(RATIO_GRID)), function(i)
    .cor_block(M[RATIO_GRID$pro[i], , drop = FALSE],
               M[RATIO_GRID$anti[i], , drop = FALSE],
               cov = COH[[coh]]$cov)$rho, numeric(1))
  RATIO_GRID %>% dplyr::mutate(cohort = coh, coexpr = v)
}))

priming <- ratio_cor %>%
  dplyr::left_join(RATIO_GRID, by = "ratio") %>%
  dplyr::left_join(dplyr::select(component_cor, cohort, axis,
                                 pro = gene, rho_pro = rho),
                   by = c("cohort", "axis", "pro")) %>%
  dplyr::left_join(dplyr::select(component_cor, cohort, axis,
                                 anti = gene, rho_anti = rho),
                   by = c("cohort", "axis", "anti")) %>%
  dplyr::left_join(coexpr, by = c("cohort", "pro", "anti", "ratio")) %>%
  dplyr::mutate(best_component = pmax(abs(rho_pro), abs(rho_anti)),
                gain = abs(rho) - best_component) %>%
  dplyr::select(cohort, axis, ratio, pro, anti, n, rho, ci_lo, ci_hi,
                rho_pro, rho_anti, coexpr, best_component, gain)

message("\n   strongest ratios against OXPHOS subunits (Spearman, both cohorts",
        " agreeing in sign):")
priming %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::select(cohort, ratio, rho, gain) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = c(rho, gain)) %>%
  dplyr::filter(sign(rho_TCGA) == sign(`rho_SCAN-B`)) %>%
  dplyr::mutate(mean_rho = (rho_TCGA + `rho_SCAN-B`) / 2) %>%
  dplyr::arrange(dplyr::desc(abs(mean_rho))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the same against MYC (", MYC_REF, "):")
priming %>%
  dplyr::filter(axis == "MYC") %>%
  dplyr::select(cohort, ratio, rho, gain) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = c(rho, gain)) %>%
  dplyr::filter(sign(rho_TCGA) == sign(`rho_SCAN-B`)) %>%
  dplyr::mutate(mean_rho = (rho_TCGA + `rho_SCAN-B`) / 2) %>%
  dplyr::arrange(dplyr::desc(abs(mean_rho))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# THE FALSIFIER. A ratio earns its name only where it beats both of its parts
# in BOTH cohorts.
message("\n   DOES THE RATIO BEAT ITS PARTS? (gain > 0 in both cohorts)")
gain_summary <- priming %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(cohort, axis, ratio, gain) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = gain) %>%
  dplyr::mutate(both_positive = TCGA > 0 & `SCAN-B` > 0)
gain_summary %>% dplyr::group_by(axis) %>%
  dplyr::summarise(n_ratios = dplyr::n(),
                   n_gain_in_both = sum(both_positive),
                   median_gain_TCGA = round(stats::median(TCGA), 3),
                   `median_gain_SCAN-B` = round(stats::median(`SCAN-B`), 3),
                   .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   the ratios that DO beat both parts in both cohorts:")
winners <- gain_summary %>% dplyr::filter(both_positive) %>%
  dplyr::arrange(dplyr::desc(pmin(TCGA, `SCAN-B`)))
if (nrow(winners)) {
  winners %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    utils::head(20) %>% as.data.frame() %>% print(row.names = FALSE)
} else {
  message("   NONE. Every ratio is a single gene wearing a ratio's name.")
}

# --- 5.0b HOW MUCH OF A RATIO IS JUST ITS TWO GENES -------------------------
# THE NUMBER THE PAPER LEANS ON, AND UNTIL 2026-09-03 IT LIVED IN NO SCRIPT.
# It was computed once in an interactive session and typed into a note, which
# is exactly the way a figure and its text drift apart. It regenerates here.
#
# The model is deliberately crude: ratio rho ~ numerator identity + denominator
# identity, as two factors, with NO interaction term. Its R-squared is the
# fraction of the 35 values that is explained by WHICH TWO GENES were used,
# ignoring the pairing entirely. What the interaction term would have captured
# - pair-specific information, the thing a ratio is supposed to add - is
# whatever is left over.
#
# A high R-squared here is therefore a NEGATIVE result about ratios, and it is
# the reason the priming subsection is descriptive.
message("\n5.0b how much of each ratio is its two component genes?")
additive_fit <- dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
  dplyr::bind_rows(lapply(unique(priming$cohort), function(co) {
    d <- priming %>% dplyr::filter(axis == ax, cohort == co)
    # 35 rows, 7 numerator levels and 5 denominator levels: 11 parameters, so
    # this is not a saturated model dressed up as a finding.
    stopifnot(nrow(d) == 35L, dplyr::n_distinct(d$pro) == 7L,
              dplyr::n_distinct(d$anti) == 5L)
    .r2 <- function(f) summary(stats::lm(f, data = d))$r.squared
    tibble::tibble(axis = ax, cohort = co, n_ratios = nrow(d),
                   r2_additive = .r2(rho ~ pro + anti),
                   r2_numerator_only = .r2(rho ~ pro),
                   r2_denominator_only = .r2(rho ~ anti))
  }))))
additive_fit %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   READ r2_additive AS A NEGATIVE. It is the share of each ratio",
        " that is fixed by\n   WHICH genes were used, with the pairing",
        " ignored - so 1 minus it is the ceiling on\n   how much any",
        " pair-specific 'priming' signal could possibly be.")

message("\n   the components on their own (Spearman, mean over cohorts):")
component_cor %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::group_by(gene, side, axis) %>%
  dplyr::summarise(mean_rho = round(mean(rho), 3), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = axis, values_from = mean_rho) %>%
  dplyr::arrange(side, dplyr::desc(OXPHOS)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.1 THE SAME RATIOS INSIDE THE LUMINAL AND BASAL COMPARTMENTS -----------
#
# WHY THIS SPLIT AND NOT ANOTHER. The mouse arm shows luminal expansion, so the
# luminal compartment is the one the question is about. And D3/S1 are the
# warning: `BCL2` against MYC is -0.369 pooled and -0.009 inside LumA, because
# the pooled value was reading the difference BETWEEN subtypes rather than
# anything within one. A pooled column that disagrees with BOTH of its strata is
# that same artefact, and it is only visible when all three are side by side.
message("\n5.1 the priming ratios by compartment (Spearman)")

.stratum_cors <- function(coh, st, B) {
  ids <- STR[[coh]][[st]]
  if (is.null(ids) || length(ids) < MIN_STRATUM_N) return(NULL)
  .cor_block(COH[[coh]]$ax[, ids, drop = FALSE], B[, ids, drop = FALSE],
             cov = COH[[coh]]$cov[ids, , drop = FALSE]) %>%
    dplyr::mutate(cohort = coh, stratum = st)
}

ratio_strata <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  R <- .ratio_matrix(GR[[coh]]$mat)
  dplyr::bind_rows(lapply(STRATA_PRIMING, .stratum_cors, coh = coh, B = R))
})) %>% dplyr::rename(ratio = item)

component_strata <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M <- GR[[coh]]$mat[PRIMING_ALL, , drop = FALSE]
  dplyr::bind_rows(lapply(STRATA_PRIMING, .stratum_cors, coh = coh, B = M))
})) %>% dplyr::rename(gene = item) %>%
  dplyr::left_join(dplyr::distinct(dplyr::select(component_cor, gene, side)),
                   by = "gene")

priming_strata <- ratio_strata %>%
  dplyr::left_join(RATIO_GRID, by = "ratio") %>%
  dplyr::left_join(dplyr::select(component_strata, cohort, stratum, axis,
                                 pro = gene, rho_pro = rho),
                   by = c("cohort", "stratum", "axis", "pro")) %>%
  dplyr::left_join(dplyr::select(component_strata, cohort, stratum, axis,
                                 anti = gene, rho_anti = rho),
                   by = c("cohort", "stratum", "axis", "anti")) %>%
  dplyr::mutate(best_component = pmax(abs(rho_pro), abs(rho_anti)),
                gain = abs(rho) - best_component) %>%
  dplyr::select(cohort, stratum, axis, ratio, pro, anti, n, rho, ci_lo, ci_hi,
                rho_pro, rho_anti, best_component, gain)

message("\n   does the POOLED value sit between its two compartments, or",
        " outside both?")
message("   (outside both is the D3/S1 signature - a between-subtype effect",
        " read as a within one)")
between_test <- priming_strata %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(cohort, axis, ratio, stratum, rho) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho) %>%
  dplyr::mutate(pooled_outside = all < pmin(Luminal, Basal) |
                                 all > pmax(Luminal, Basal))
between_test %>% dplyr::group_by(cohort, axis) %>%
  dplyr::summarise(n_ratios = dplyr::n(),
                   n_pooled_outside_both = sum(pooled_outside),
                   .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the ratios where the two compartments most disagree",
        " (both cohorts agreeing on the sign of the difference):")
lum_basal_gap <- priming_strata %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"),
                stratum %in% c("Luminal", "Basal")) %>%
  dplyr::select(cohort, axis, ratio, stratum, rho) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho) %>%
  dplyr::mutate(diff = Luminal - Basal) %>%
  dplyr::select(cohort, axis, ratio, diff) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = diff) %>%
  dplyr::filter(sign(TCGA) == sign(`SCAN-B`))
lum_basal_gap %>%
  dplyr::arrange(dplyr::desc(abs(TCGA) + abs(`SCAN-B`))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::rename(`Luminal-Basal, TCGA` = TCGA,
                `Luminal-Basal, SCAN-B` = `SCAN-B`) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the 12 component genes by compartment (Spearman, OXPHOS):")
component_strata %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::select(cohort, gene, side, stratum, rho) %>%
  tidyr::pivot_wider(names_from = c(cohort, stratum), values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  dplyr::arrange(side, gene) %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   NOTE: Basal is 171 TCGA and 317 SCAN-B samples. A 171-sample",
        " stratum gives a\n   95% interval about +/- 0.15 wide on rho, so a",
        " Basal-versus-Luminal difference\n   smaller than about 0.2 is not",
        " separable from sampling in TCGA.")

# =============================================================================
# 6. Figures
# =============================================================================
message("\n6. figures")

COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
EFFECT_COLS <- c(`pro-death` = "#d7191c", `pro-survival` = "#2c7bb6")
theme_e10 <- ggplot2::theme_bw(base_size = 9) +
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

# ONE ROW ORDER FOR ALL FOUR MACHINERY PANELS, taken from the OXPHOS Spearman
# column - which is exactly E08 fig6. The figures are meant to be laid side by
# side and read for MOVEMENT; re-sorting each panel by its own values would
# make every panel look like a gradient and hide the only thing worth seeing.
ORDER_KEY <- canon_wide %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::select(gene, order_rho = mean_rho)
MODULE_ORDER <- ORDER_KEY %>%
  dplyr::left_join(dplyr::select(reannot, gene, cdc_module), by = "gene") %>%
  dplyr::group_by(cdc_module) %>%
  dplyr::summarise(m = stats::median(order_rho), .groups = "drop") %>%
  dplyr::arrange(m) %>% dplyr::pull(cdc_module)
GENE_ORDER <- ORDER_KEY %>% dplyr::arrange(order_rho) %>% dplyr::pull(gene)

.machinery_plot <- function(ax, xlab, title, caption, flag_genes, flag_label) {
  d <- canon_wide %>%
    dplyr::filter(axis == ax) %>%
    dplyr::left_join(ORDER_KEY, by = "gene") %>%
    dplyr::mutate(gene = factor(gene, levels = GENE_ORDER),
                  cdc_module = factor(cdc_module, levels = MODULE_ORDER))
  # The row labels carry the expression flag, so a reader cannot look at a
  # low-expression gene without being told it is one.
  lab <- stats::setNames(paste0(GENE_ORDER,
                                ifelse(GENE_ORDER %in% LOW_EXPR_GENES, " *", "")),
                         GENE_ORDER)
  med <- d %>% dplyr::group_by(cdc_module) %>%
    dplyr::summarise(m = stats::median(mean_rho), .groups = "drop")
  flagged <- dplyr::filter(d, gene %in% flag_genes)
  sub <- paste0("EXPLORATORY - not pre-registered | adjusted for proliferation | circle = in MitoCarta 3.0",
                if (nrow(flagged)) paste0(" | cross = ", flag_label) else "")
  p <- ggplot2::ggplot(d, ggplot2::aes(mean_rho, gene, colour = effect)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(data = med, ggplot2::aes(xintercept = m),
                        linetype = 2, linewidth = 0.4, colour = "grey40") +
    ggplot2::geom_segment(ggplot2::aes(x = TCGA, xend = `SCAN-B`,
                                       y = gene, yend = gene),
                          linewidth = 1.1, alpha = 0.5) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_point(data = dplyr::filter(d, mitocarta),
                        ggplot2::aes(x = mean_rho), shape = 1, size = 3.4,
                        colour = "black", stroke = 0.5) +
    ggplot2::facet_grid(cdc_module ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::scale_colour_manual(values = EFFECT_COLS, name = NULL) +
    ggplot2::scale_y_discrete(labels = lab, drop = TRUE) +
    ggplot2::labs(title = title, subtitle = sub, x = xlab, y = NULL,
                  caption = caption) +
    theme_e10 +
    ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, size = 7),
                   strip.placement = "outside")
  if (nrow(flagged)) {
    p <- p + ggplot2::geom_point(data = flagged, ggplot2::aes(x = mean_rho),
                                 shape = 4, size = 2.6, colour = "grey20",
                                 stroke = 0.7, show.legend = FALSE)
  }
  p
}

# THE CROSS MUST MEAN THE SAME THING AS THE CAPTION SAYS IT DOES, and what
# threatens independence is not the same on the two axes:
#   OXPHOS panels  a gene inside the `OXPHOS subunits` arm is correlated partly
#                  with itself. That is CYCS, 1 of 89.
#   MYC panels     FELSHER__MITOSTRIP contains none of these 44, so there is no
#                  self-overlap - but trap 3 says every MYC signature is
#                  entangled with proliferation, and two of the 44 are HALLMARK
#                  E2F/G2M genes. Those are the cells to distrust on this axis.
# M_b__MITOSTRIP contains seven of the 44 and is NOT plotted anywhere for that
# reason; the overlap is in the saved table as `in_m_b`.
FLAG_OXPHOS <- intersect(CANON_GENES, OX_ARM)
FLAG_MYC    <- intersect(CANON_GENES, sd_$strip_refs$PROLIF_REF)
LOW_EXPR_GENES <- sort(unique(expr_rank$gene[expr_rank$low_expression]))
message("   below the ", LOW_EXPR_PCT * 100, "th expression percentile in at",
        " least one cohort, starred on every panel: ",
        if (length(LOW_EXPR_GENES)) paste(LOW_EXPR_GENES, collapse = ", ") else
          "none")
message("\n   figure flags: OXPHOS panels cross ",
        paste(FLAG_OXPHOS, collapse = ", "), " | MYC panels cross ",
        paste(FLAG_MYC, collapse = ", "))

# Captions are HARD-WRAPPED at about 85 characters. ggplot does not wrap a
# caption, it clips it, and a clipped caption is a caption that lies by omission.
CAP_MITO <- paste0(
  "Circles mark the 20 genes MitoCarta 3.0 places in the mitochondrion - an INDEPENDENT\n",
  "localisation call, not the curation's own. The dashed line is the module median, and\n",
  "bars span the two cohorts, so a short bar is a gene that replicates. Rows and module\n",
  "order are FIXED to the OXPHOS panel in both figures, so the two can be laid side by\n",
  "side and read for movement. A STARRED gene sits below the 25th expression percentile\n",
  "in at least one cohort: its rho is largely a correlation of quantisation noise with a\n",
  "score and must not be read as a result.\n")

.save(.machinery_plot(
  "OXPHOS",
  "mean partial Spearman rho with OXPHOS subunits, across cohorts",
  "The canonical machinery against OXPHOS",
  paste0(CAP_MITO,
    "The cross is CYCS, 1 of the 89 genes in the OXPHOS subunits arm itself, so its value\n",
    "is partly a correlation with itself. This is the E08 fig6 panel, recomputed here and\n",
    "asserted identical to it."),
  FLAG_OXPHOS, "in the OXPHOS subunits arm itself"),
  "E10_fig1_machinery_oxphos", 9, 9)

.save(.machinery_plot(
  "MYC",
  paste0("mean partial Spearman rho with MYC activity (", MYC_REF,
         "), across cohorts"),
  "The canonical machinery against MYC",
  paste0(CAP_MITO,
    "FELSHER__MITOSTRIP contains none of these 44 genes, so there is no self-overlap on\n",
    "this axis. The crosses are the two that are HALLMARK E2F/G2M members - BIRC5 at the\n",
    "top of the range is survivin, and a MYC signature is entangled with proliferation by\n",
    "construction. This is NOT the OXPHOS panel re-drawn; the rows are in the same order,\n",
    "so compare the two gene by gene."),
  FLAG_MYC, "a HALLMARK E2F/G2M proliferation gene"),
  "E10_fig2_machinery_myc", 9, 9)

# --- WHICH RATIO CELLS ARE WORTH LOOKING AT ----------------------------------
# Author, 2026-09-03: put figure 4's test onto the heatmaps, so the two do not
# have to be read side by side. A cell is MARKED when it passes both of the
# things this script already computes separately:
#
#   `gain > 0`           |rho of the ratio| exceeds the larger |rho| of the two
#                        genes it is made of. That is figure 4's diagonal,
#                        evaluated per cell instead of drawn as a scatter.
#   `|rho| >= 0.30`      a hand-drawn effect-size floor. It is NOT a test, it is
#                        not derived from anything, and moving it moves the
#                        marks. It is |rho| and not rho so that a ratio running
#                        strongly the OTHER way is marked too - five of the
#                        marked cells are negative, and hiding them would make
#                        the mark mean "large and positive" while the caption
#                        said "large".
#
# TWO LEVELS OF MARK, because one cohort is not a result (CLAUDE.md section 2,
# and this script's own `ratios` rule):
#   *                    the cell passes IN THIS COHORT.
#   heavy black border   the same ratio passes in BOTH cohorts, at the same
#                        stratum and on the same axis. That is the only version
#                        of the mark this study's rules allow anyone to quote.
#
# The asterisk is deliberately the weaker-looking mark of the two.
MARK_MIN_ABS_RHO <- 0.30

.mark_cells <- function(d, keys) {
  d <- d %>% dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
    dplyr::mutate(mark = gain > 0 & abs(rho) >= MARK_MIN_ABS_RHO)
  both <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::summarise(n_marked = sum(mark),
                     n_cohorts = dplyr::n_distinct(cohort), .groups = "drop")
  # A pair of cohorts is required, not merely two marks: if a key ever appeared
  # in one cohort twice the count would pass while the agreement did not.
  stopifnot(all(both$n_cohorts == 2L))
  d %>% dplyr::left_join(
    both %>% dplyr::transmute(dplyr::across(dplyr::all_of(keys)),
                              mark_both = n_marked == 2L), by = keys)
}
priming_marked <- .mark_cells(priming, c("axis", "ratio"))
strata_marked  <- .mark_cells(priming_strata, c("stratum", "axis", "ratio"))
stopifnot(!anyNA(priming_marked$mark_both), !anyNA(strata_marked$mark_both))

# --- WHICH GENES ARE GOVERNED BY OXPHOS RATHER THAN BY MYC ------------------
# Author, 2026-09-03: among the marked cells, mark harder the ones whose TWO
# genes are BOTH governed more by OXPHOS than by MYC. The per-gene quantity is
# E15 figure 5's - rho(gene, OXPHOS) - rho(gene, MYC), both partial on
# PROLIF_DISJOINT - and it is RECOMPUTED HERE from `component_cor` rather than
# read from E15. Two reasons: E15 runs after this script, and E11's pipeline
# covers only 10 of these 12 genes, because BIK and BCL2L2 are not among the
# canonical 44. On the 10 it does cover the two pipelines agree to 0.000 in all
# 20 cohort-by-gene cells, checked 2026-09-03.
#
# THE THRESHOLD IS ON THE MAGNITUDE, |gap| >= 0.20, chosen by the author over
# the signed version. It admits MCL1, whose gap is -0.32 / -0.25: MCL1 leans
# hard to OXPHOS while running the OTHER WAY. That is not a loophole, it is the
# productive configuration - R5 found that a ratio only beats its own two genes
# when they move oppositely, so a denominator leaning negative is exactly what
# lets log2(pro/anti) add rather than cancel. Under the signed rule
# (gap >= +0.20) NOT ONE marked cell qualifies anywhere, which is the same
# structural tension seen from the other side; both counts are printed below.
MARK2_MIN_ABS_GAP <- 0.20

gene_gap <- component_cor %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(cohort, gene, side, axis, rho) %>%
  tidyr::pivot_wider(names_from = axis, values_from = rho) %>%
  dplyr::mutate(gap = OXPHOS - MYC)
stopifnot(nrow(gene_gap) == length(c(PRIMING_PRO, PRIMING_ANTI)) * 2L,
          !anyNA(gene_gap$gap))
gene_lean <- gene_gap %>%
  dplyr::group_by(gene, side) %>%
  dplyr::summarise(gap_TCGA = gap[cohort == "TCGA"],
                   gap_SCANB = gap[cohort == "SCAN-B"],
                   max_abs_gap = max(abs(gap)),
                   max_signed_gap = max(gap), .groups = "drop") %>%
  dplyr::mutate(oxphos_led = max_abs_gap >= MARK2_MIN_ABS_GAP,
                oxphos_led_signed = max_signed_gap >= MARK2_MIN_ABS_GAP) %>%
  dplyr::arrange(side, dplyr::desc(max_abs_gap))
LEAN_GENES <- sort(gene_lean$gene[gene_lean$oxphos_led])
LEAN_SIGNED <- sort(gene_lean$gene[gene_lean$oxphos_led_signed])

.mark2 <- function(d)
  dplyr::mutate(d, mark2 = mark & pro %in% LEAN_GENES & anti %in% LEAN_GENES,
                mark2_signed = mark & pro %in% LEAN_SIGNED &
                  anti %in% LEAN_SIGNED)
priming_marked <- .mark2(priming_marked)
strata_marked  <- .mark2(strata_marked)

message("\n   per-gene |rho(OXPHOS) - rho(MYC)|, the E15 figure 5 quantity:")
gene_lean %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   OXPHOS-led at |gap| >= ", MARK2_MIN_ABS_GAP, ": ",
        paste(LEAN_GENES, collapse = ", "))
message("   under the SIGNED rule it would be: ",
        paste(LEAN_SIGNED, collapse = ", "), " - and no marked cell qualifies")

message("\n   heatmap marks: `gain > 0` AND |rho| >= ", MARK_MIN_ABS_RHO)
mark_summary <- dplyr::bind_rows(
  priming_marked %>% dplyr::mutate(stratum = "all (figure 3)"),
  strata_marked %>% dplyr::mutate(stratum = paste0(stratum, " (figure 6)"))) %>%
  dplyr::group_by(figure = ifelse(grepl("figure 3", stratum), "fig3", "fig6"),
                  stratum, axis) %>%
  dplyr::summarise(n_cells = dplyr::n(), n_marked = sum(mark),
                   n_marked_both = sum(mark_both) / 2,
                   n_mark2 = sum(mark2), n_mark2_signed = sum(mark2_signed),
                   .groups = "drop")
mark_summary %>% as.data.frame() %>% print(row.names = FALSE)
message("   NOT ONE MYC CELL IS MARKED, in either heatmap, at any stratum.")
message("   ** cells - both genes OXPHOS-led - on figure 6:")
strata_marked %>% dplyr::filter(mark2) %>%
  dplyr::transmute(stratum, cohort, axis, ratio, rho = round(rho, 3),
                   gain = round(gain, 3), bordered = mark_both) %>%
  dplyr::arrange(stratum, cohort, ratio) %>%
  as.data.frame() %>% print(row.names = FALSE)
mark_both_list <- strata_marked %>%
  dplyr::filter(mark_both) %>%
  dplyr::distinct(stratum, axis, ratio) %>%
  dplyr::arrange(stratum, axis, ratio)
message("   marked in BOTH cohorts on figure 6 (heavy border):")
if (nrow(mark_both_list)) {
  mark_both_list %>% as.data.frame() %>% print(row.names = FALSE)
} else message("      none")
message("   marked in BOTH cohorts on figure 3 (pooled): ",
        sum(priming_marked$mark_both) / 2)
# WHY there are none, which is the pooled heatmap's actual result: the ratios
# that replicate as ratios are all small, and the large ones do not replicate.
pooled_two_way <- priming_marked %>%
  dplyr::group_by(axis, ratio) %>%
  dplyr::summarise(gain_both = sum(gain > 0) == 2L,
                   max_abs_rho = max(abs(rho)), .groups = "drop") %>%
  dplyr::group_by(axis) %>%
  dplyr::summarise(
    n_gain_both = sum(gain_both),
    # max() of an empty vector is -Inf with a warning; there are no OXPHOS
    # ratios in this group at all, and NA says that rather than pretending.
    max_abs_rho_of_those = if (any(gain_both)) max(max_abs_rho[gain_both]) else
      NA_real_,
    n_ratios_reaching_floor = sum(max_abs_rho >= MARK_MIN_ABS_RHO),
    .groups = "drop")
message("   pooled, and this is why: which ratios gain in BOTH cohorts, and how",
        " big they get")
pooled_two_way %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# A ^ on an axis label marks a gene that is OXPHOS-led on its own, so a reader
# can see WHERE a ** could ever appear - only at the crossing of a marked row
# and a marked column - and can see the crossings that stayed unstarred.
.axis_lab <- function(v)
  stats::setNames(paste0(v, ifelse(v %in% LEAN_GENES, " ^", "")), v)

# --- the priming ratios ------------------------------------------------------
g6dat <- priming_marked %>%
  dplyr::mutate(pro = factor(pro, levels = rev(PRIMING_PRO)),
                anti = factor(anti, levels = PRIMING_ANTI),
                cohort = factor(cohort, levels = names(COHORT_COLS)))
LIM <- max(abs(g6dat$rho))
g6 <- ggplot2::ggplot(g6dat, ggplot2::aes(anti, pro, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.1) +
  ggplot2::geom_tile(data = dplyr::filter(g6dat, mark_both),
                     ggplot2::aes(anti, pro), fill = NA, colour = "black",
                     linewidth = 0.85, inherit.aes = FALSE) +
  ggplot2::geom_text(data = dplyr::filter(g6dat, mark),
                     ggplot2::aes(anti, pro, label = ifelse(mark2, "**", "*")),
                     size = 4.2, fontface = "bold", nudge_x = 0.29,
                     nudge_y = 0.20, vjust = 0.75, inherit.aes = FALSE) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_x_discrete(labels = .axis_lab(PRIMING_ANTI)) +
  ggplot2::scale_y_discrete(labels = .axis_lab(PRIMING_PRO)) +
  ggplot2::scale_fill_gradient2(low = "#2c7bb6", mid = "grey96",
                                high = "#d7191c", midpoint = 0,
                                limits = c(-LIM, LIM), na.value = "grey85",
                                name = "partial Spearman rho of log2(pro/anti)") +
  ggplot2::labs(
    title = "Every BCL2-family priming ratio against MYC and against OXPHOS",
    # Hard-wrapped for the same reason the captions are: ggplot CLIPS a long
    # subtitle rather than wrapping it, and at 8in this one lost its last three
    # words for every run before 2026-09-03.
    subtitle = paste0("EXPLORATORY - not pre-registered | adjusted for ",
                      "proliferation | ", nrow(RATIO_GRID), " ratios\n",
                      "rows are the pro-apoptotic numerator, columns the ",
                      "anti-apoptotic denominator"),
    x = "anti-apoptotic (denominator)", y = "pro-apoptotic (numerator)",
    # paste0 does not insert spaces, so a rendered caption line may be split
    # across two source strings to keep the file inside 80 columns.
    caption = paste0(
      "RED means the ratio rises with the axis - a MORE PRIMED transcriptome\n",
      "in high-MYC or high-OXPHOS tumours. A COLUMN that is uniformly coloured\n",
      "is the denominator gene talking rather than priming, which is what the\n",
      "BCL2 column is: check `gain` before reading any cell as a ratio, because\n",
      "most of these ratios are beaten by one of the two genes they are made\n",
      "of.\n",
      "* MARKS A CELL THAT PASSES FIGURE 4'S TEST AND REACHES |rho| >= 0.30 -\n",
      "the ratio beats the stronger of its own two genes AND the effect is not\n",
      "small. 7 of 140 cells qualify and ALL SEVEN ARE OXPHOS; no MYC cell in\n",
      "this figure is marked. The 0.30 line is drawn by hand and is not a test.\n",
      "** UPGRADES A STARRED CELL whose TWO genes are each governed more by\n",
      "OXPHOS than by MYC on their own: |rho(OXPHOS) - rho(MYC)| >= 0.20 in at\n",
      "least one cohort. Those five genes carry a ^ on the axis - BAD, BBC3 and\n",
      "BIK above, BCL2L1 and MCL1 across - so a ** can only ever appear where a\n",
      "marked row crosses a marked column. Two do here, both TCGA and both with\n",
      "MCL1 underneath. MCL1 qualifies by running the OTHER WAY (-0.32 / -0.25)\n",
      "and that is the productive case, not a loophole: a ratio beats its own\n",
      "genes only when they move oppositely. Under a SIGNED rule no cell in\n",
      "either heatmap qualifies at all.\n",
      "A HEAVY BORDER would mark a ratio passing in BOTH cohorts. THERE ARE\n",
      "NONE HERE, and the two halves of that failure pull opposite ways: the 5\n",
      "ratios that beat their parts in both cohorts are ALL on MYC and their\n",
      "|rho| tops out at 0.27, while the 7 cells reaching 0.30 are ALL on OXPHOS\n",
      "and not one of them beats its parts in both cohorts. Figure 6 is where\n",
      "three cells clear both bars at once, and all three are Basal.")) +
  theme_e10 +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.key.width = ggplot2::unit(1.4, "cm"))
.save(g6, "E10_fig3_priming_ratio_heatmap", 8, 7)

# Does the ratio beat its parts? The falsifier, drawn.
g7dat <- priming %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                abs_rho = abs(rho),
                beats_parts = gain > 0)
g7 <- ggplot2::ggplot(g7dat, ggplot2::aes(best_component, abs_rho,
                                          colour = beats_parts)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::geom_point(size = 1.4, alpha = 0.75) +
  ggrepel::geom_text_repel(data = dplyr::filter(g7dat, gain > 0.045),
                           ggplot2::aes(label = ratio), size = 2.1,
                           max.overlaps = 12, min.segment.length = 0.2,
                           show.legend = FALSE, seed = PROJECT_SEED) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "#1b9e77",
                                          `FALSE` = "grey60"),
                               name = "the ratio beats both its genes") +
  ggplot2::labs(
    title = "Does a priming ratio measure more than its stronger half?",
    subtitle = paste("EXPLORATORY - not pre-registered | adjusted for",
                     "proliferation | each point is one of", nrow(RATIO_GRID),
                     "ratios"),
    x = "|partial rho| of the stronger of the two component genes",
    y = "|partial rho| of the log2 ratio",
    caption = paste0(
      "ABOVE the diagonal the ratio adds information; ON or BELOW it, the single\n",
      "gene is the better measurement and the ratio is that gene with noise added\n",
      "to it. Only ratios above the line IN BOTH COHORTS are worth reporting as\n",
      "ratios, and the median ratio sits below it on both axes.\n",
      "THIS TEST IS CARRIED ONTO THE HEATMAPS. A cell on figure 3 or figure 6 is\n",
      "starred when it is green here AND reaches |rho| >= 0.30, and is bordered\n",
      "when it does so in both cohorts. Being above the diagonal is necessary and\n",
      "not sufficient: many green points here are small effects.")) +
  theme_e10
.save(g7, "E10_fig4_ratio_vs_components", 8, 7)

# The 12 component genes on their own, which is what most of the heatmap is.
# Ordered by the mean OXPHOS Spearman, stated rather than inferred - reorder()
# on a column that mixes two axes would silently average them.
PRIMING_ORDER <- component_cor %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(m = mean(rho), .groups = "drop") %>%
  dplyr::arrange(m) %>% dplyr::pull(gene)
g8dat <- component_cor %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                gene = factor(gene, levels = PRIMING_ORDER))
g8 <- ggplot2::ggplot(g8dat, ggplot2::aes(rho, gene, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                          position = ggplot2::position_dodge(width = 0.5),
                          linewidth = 0.4, alpha = 0.6) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5),
                      size = 1.9) +
  ggplot2::facet_grid(side ~ axis, scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(
    title = "The 12 priming genes on their own, before any ratio is taken",
    subtitle = paste("EXPLORATORY - not pre-registered | adjusted for",
                     "proliferation | the author's pro and anti lists"),
    x = "partial Spearman rho with the axis", y = NULL,
    caption = paste0(
      "This is the panel the ratio heatmap is built out of. Where a pro and an\n",
      "anti gene move the SAME way their ratio cancels; where they move\n",
      "oppositely the ratio adds them. Neither list moves as a block - BCL2L1\n",
      "and BCL2 are both anti and sit at opposite ends of the OXPHOS axis - and\n",
      "that is why most of the ratios lose to their own genes. Every gene here\n",
      "clears the 25th expression percentile in both cohorts.")) +
  theme_e10
.save(g8, "E10_fig5_priming_components", 8, 6)

# --- the priming ratios by compartment ---------------------------------------
STRAT_COLS <- c(all = "grey30", Luminal = "#7b3294", Basal = "#008837")

g9dat <- strata_marked %>%
  dplyr::mutate(pro = factor(pro, levels = rev(PRIMING_PRO)),
                anti = factor(anti, levels = PRIMING_ANTI),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                stratum = factor(stratum, levels = STRATA_PRIMING))
LIM9 <- max(abs(g9dat$rho))
g9 <- ggplot2::ggplot(g9dat, ggplot2::aes(anti, pro, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2) +
  ggplot2::geom_tile(data = dplyr::filter(g9dat, mark_both),
                     ggplot2::aes(anti, pro), fill = NA, colour = "black",
                     linewidth = 0.85, inherit.aes = FALSE) +
  ggplot2::geom_text(data = dplyr::filter(g9dat, mark),
                     ggplot2::aes(anti, pro, label = ifelse(mark2, "**", "*")),
                     size = 4, fontface = "bold", nudge_x = 0.29,
                     nudge_y = 0.20, vjust = 0.75, inherit.aes = FALSE) +
  ggplot2::facet_grid(stratum ~ cohort + axis) +
  ggplot2::scale_x_discrete(labels = .axis_lab(PRIMING_ANTI)) +
  ggplot2::scale_y_discrete(labels = .axis_lab(PRIMING_PRO)) +
  ggplot2::scale_fill_gradient2(low = "#2c7bb6", mid = "grey96",
                                high = "#d7191c", midpoint = 0,
                                limits = c(-LIM9, LIM9), na.value = "grey85",
                                name = "Spearman rho of log2(pro/anti)") +
  ggplot2::labs(
    title = "The priming ratios inside the luminal and basal compartments",
    subtitle = paste0("EXPLORATORY - not pre-registered | adjusted for ",
                      "proliferation | Luminal = LumA + LumB (696 TCGA / 2,436 ",
                      "SCAN-B), Basal (171 / 317)"),
    x = "anti-apoptotic (denominator)", y = "pro-apoptotic (numerator)",
    caption = paste0(
      "* MARKS A CELL THAT PASSES FIGURE 4'S TEST AND REACHES |rho| >= 0.30: the ratio beats the stronger of the two genes it is\n",
      "made of AND the effect is not small. 22 of the 420 cells qualify and ALL 22 ARE OXPHOS - not one MYC cell is marked, at any\n",
      "stratum, in either cohort. The 0.30 line is drawn by hand, is not a test, and moving it moves the marks.\n",
      "** UPGRADES A STARRED CELL whose TWO genes are each governed more by OXPHOS than by MYC ON THEIR OWN - |rho(OXPHOS) - rho(MYC)|\n",
      ">= 0.20 in at least one cohort, which is the quantity E15 figure 5 draws. The five genes that qualify carry a ^ on the axis:\n",
      "BAD, BBC3, BIK as numerators, BCL2L1 and MCL1 as denominators. A ** can therefore only appear where a marked row crosses a\n",
      "marked column, and 7 of the 22 do - EVERY ONE OF THEM WITH MCL1 UNDERNEATH. MCL1 qualifies by running the OTHER WAY (-0.32 /\n",
      "-0.25), which is the productive case rather than a loophole: a ratio beats its own two genes only when they move oppositely,\n",
      "so a denominator leaning the other way is exactly what makes log2(pro/anti) add instead of cancel. Under a SIGNED rule\n",
      "(gap >= +0.20) not one cell in either heatmap qualifies, which is that same tension seen from the other side. The 0.20 line,\n",
      "like the 0.30 one, is drawn by hand and is not a test.\n",
      "A HEAVY BORDER is the only mark this study's rules allow anyone to quote: the same ratio passing in BOTH cohorts, at the same\n",
      "stratum and on the same axis. There are THREE, all Basal and all OXPHOS - BBC3/BCL2, BID/BCL2 and BBC3/MCL1 - and none at all\n",
      "in the pooled row or in Luminal. BBC3/MCL1 IN BASAL IS THE ONLY CELL IN THE FIGURE CARRYING ALL THREE MARKS - starred in both\n",
      "cohorts, bordered, and ** in both. Read that against where the marks are densest: Basal is 171 TCGA and 317 SCAN-B samples, so\n",
      "it is also where the intervals are widest, and a mark there is a weaker claim than the same mark would be in Luminal.\n",
      "The `all` row is the pooled value and is NOT an average of the two below it. A pooled cell that sits OUTSIDE the range of\n",
      "its own two compartments is reading a difference BETWEEN subtypes rather than anything within one - the D3/S1 artefact,\n",
      "where BCL2 against MYC is -0.369 pooled and -0.009 inside LumA. Basal is 171 TCGA samples, so a 95% interval there is\n",
      "about +/- 0.15 wide and a Luminal-Basal difference under roughly 0.2 is not separable from sampling in that cohort.")) +
  theme_e10 +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.key.width = ggplot2::unit(1.4, "cm"))
.save(g9, "E10_fig6_priming_by_compartment", 13, 10)

# --- FIG 8: the two axes as ONE number, and how that number must be formed ---
# Author, 2026-09-03: replace figure 3's two facet columns - MYC beside OXPHOS -
# with a single column holding their difference, so the cells where OXPHOS wins
# stand out instead of having to be compared across the page.
#
# THE QUESTION THE AUTHOR ASKED, AND IT IS THE RIGHT ONE: what happens to the
# negative correlations. There are two candidate differences and they are NOT
# the same quantity:
#
#   SIGNED     rho(OXPHOS) - rho(MYC). Keeps the DIRECTION of the relationship.
#              Its sign says which way the ratio moves, not which axis governs
#              it. `BCL2L11/BCL2L1` in TCGA is -0.106 on MYC and -0.515 on
#              OXPHOS: OXPHOS governs it five times over, and this quantity
#              scores it -0.409, the same deep blue a MYC-dominated cell would
#              get. On a diverging scale the OXPHOS-dominated cells therefore
#              land at BOTH ends and the figure cannot be read at a glance,
#              which is exactly what the new figure was supposed to fix.
#   WHICH AXIS |rho(OXPHOS)| - |rho(MYC)|. The difference OF THE MAGNITUDES,
#              not the magnitude of the difference. Its sign says WHICH AXIS
#              GOVERNS and nothing about direction, so one diverging scale is
#              unambiguous: red = OXPHOS, blue = MYC. This is the same
#              construction as `gain` and it is the panel to read.
#
# THEY DISAGREE IN SIGN ON 22 OF THE 70 CELLS, so this is not a quibble. Both
# are drawn, side by side, because the disagreement IS the answer to the
# question and asserting it in a caption would be weaker than showing it.
# Direction is not lost: it stays in figure 3, and the `*` `**` and border marks
# are carried here so the two figures line up cell for cell.
DIFF_LEVELS <- c(
  "SIGNED   rho(OXPHOS) - rho(MYC)\nred = the ratio RISES more with OXPHOS",
  "WHICH AXIS   |rho(OXPHOS)| - |rho(MYC)|\nred = OXPHOS GOVERNS it more")

diff_wide <- priming_marked %>%
  dplyr::select(cohort, axis, ratio, pro, anti, rho, mark, mark2, mark_both) %>%
  tidyr::pivot_wider(names_from = axis,
                     values_from = c(rho, mark, mark2, mark_both))
# Every mark in this study is an OXPHOS mark - section 6 asserts no MYC cell is
# ever marked - so carrying `_OXPHOS` here is a simplification, not a choice,
# and it stops being true silently if that ever changes.
stopifnot(nrow(diff_wide) == nrow(RATIO_GRID) * 2L,
          !any(diff_wide$mark_MYC), !any(diff_wide$mark2_MYC),
          !any(diff_wide$mark_both_MYC))
diff_tab <- diff_wide %>%
  dplyr::transmute(cohort, ratio, pro, anti,
                   rho_MYC, rho_OXPHOS,
                   signed = rho_OXPHOS - rho_MYC,
                   which_axis = abs(rho_OXPHOS) - abs(rho_MYC),
                   axes_discordant = sign(rho_OXPHOS) != sign(rho_MYC),
                   readings_disagree = sign(signed) != sign(which_axis),
                   mark = mark_OXPHOS, mark2 = mark2_OXPHOS,
                   mark_both = mark_both_OXPHOS)
# A caption number must come from the object, never from a typed literal, or
# the caption and the figure drift apart on the next run. This pulls exactly one
# row and stops if it does not.
.one <- function(d, ...) {
  v <- dplyr::filter(d, ...)
  stopifnot(nrow(v) == 1L)
  v
}
N_DISAGREE <- sum(diff_tab$readings_disagree)
N_DISCORD  <- sum(diff_tab$axes_discordant)
message("\n   figure 8: the two ways to difference the axes")
message("   the SIGNED and the WHICH-AXIS reading disagree in sign on ",
        N_DISAGREE, " of ", nrow(diff_tab), " cells")
message("   the two axes disagree about DIRECTION on ", N_DISCORD, " cells")
message("   the starred cells where the two readings flip:")
diff_tab %>% dplyr::filter(mark, readings_disagree) %>%
  dplyr::transmute(cohort, ratio, MYC = round(rho_MYC, 3),
                   OXPHOS = round(rho_OXPHOS, 3),
                   signed = round(signed, 3),
                   which_axis = round(which_axis, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

g11dat <- dplyr::bind_rows(
  diff_tab %>% dplyr::mutate(quantity = DIFF_LEVELS[1], value = signed),
  diff_tab %>% dplyr::mutate(quantity = DIFF_LEVELS[2], value = which_axis)) %>%
  dplyr::mutate(quantity = factor(quantity, levels = DIFF_LEVELS),
                pro = factor(pro, levels = rev(PRIMING_PRO)),
                anti = factor(anti, levels = PRIMING_ANTI),
                cohort = factor(cohort, levels = names(COHORT_COLS)))
LIM11 <- max(abs(g11dat$value))
FLIP <- diff_tab %>% dplyr::filter(mark, readings_disagree) %>%
  dplyr::mutate(lab = paste0(ratio, " in ", cohort)) %>% dplyr::pull(lab)

CAP11 <- sprintf(paste0(
  "ONE NUMBER PER RATIO INSTEAD OF TWO, and the point of the figure is that",
  " there are two ways to form\n",
  "it and only the right-hand one can be read at a glance.\n",
  "LEFT, the SIGNED difference. It keeps the DIRECTION of the relationship, so",
  " its sign says which way\n",
  "the ratio moves and NOT which axis governs it. RIGHT, the difference OF THE",
  " MAGNITUDES - not the\n",
  "magnitude of the difference. Its sign says only WHICH AXIS GOVERNS, so red",
  " means OXPHOS everywhere\n",
  "in that panel and direction is not encoded at all.\n",
  "THE TWO DISAGREE IN SIGN ON %d OF THE %d CELLS, so the choice is not",
  " cosmetic. %s is the clearest\n",
  "case: %.2f on MYC against %.2f on OXPHOS, which is OXPHOS governing it",
  " about five times over. The\n",
  "signed panel scores it %+.2f - the same deep blue a MYC-dominated cell",
  " would get - and the which-axis\n",
  "panel scores it %+.2f. Both starred cells that flip do so this way, and",
  " both are ratios whose\n",
  "correlation with OXPHOS is strongly NEGATIVE.\n",
  "SO: USE THE RIGHT-HAND PANEL TO FIND WHERE OXPHOS WINS, AND FIGURE 3 TO",
  " READ THE DIRECTION. The two\n",
  "cannot be one number. On %d of the 70 cells the axes point in opposite",
  " directions altogether, and\n",
  "there the signed difference ADDS the two correlations rather than",
  " differencing them.\n",
  "The * ** ^ and border marks are figure 3's, carried unchanged so the two",
  " figures line up cell for\n",
  "cell: * the ratio beats its own two genes and reaches |rho| 0.30 in this",
  " cohort; ** its two genes\n",
  "are each OXPHOS-led on their own; ^ on an axis marks those genes; a heavy",
  " border would mark a cell\n",
  "passing in both cohorts, and pooled there are none."),
  N_DISAGREE, nrow(diff_tab), FLIP[1],
  .one(diff_tab, ratio == "BCL2L11/BCL2L1", cohort == "TCGA")$rho_MYC,
  .one(diff_tab, ratio == "BCL2L11/BCL2L1", cohort == "TCGA")$rho_OXPHOS,
  .one(diff_tab, ratio == "BCL2L11/BCL2L1", cohort == "TCGA")$signed,
  .one(diff_tab, ratio == "BCL2L11/BCL2L1", cohort == "TCGA")$which_axis,
  N_DISCORD)

g11 <- ggplot2::ggplot(g11dat, ggplot2::aes(anti, pro, fill = value)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 2.1) +
  ggplot2::geom_tile(data = dplyr::filter(g11dat, mark_both),
                     ggplot2::aes(anti, pro), fill = NA, colour = "black",
                     linewidth = 0.85, inherit.aes = FALSE) +
  ggplot2::geom_text(data = dplyr::filter(g11dat, mark),
                     ggplot2::aes(anti, pro, label = ifelse(mark2, "**", "*")),
                     size = 4.2, fontface = "bold", nudge_x = 0.29,
                     nudge_y = 0.20, vjust = 0.75, inherit.aes = FALSE) +
  ggplot2::facet_grid(cohort ~ quantity) +
  ggplot2::scale_x_discrete(labels = .axis_lab(PRIMING_ANTI)) +
  ggplot2::scale_y_discrete(labels = .axis_lab(PRIMING_PRO)) +
  ggplot2::scale_fill_gradient2(low = "#2c7bb6", mid = "grey96",
                                high = "#d7191c", midpoint = 0,
                                limits = c(-LIM11, LIM11), na.value = "grey85",
                                name = paste("difference in partial Spearman",
                                             "rho - the two panels differ in",
                                             "WHAT the sign means")) +
  ggplot2::labs(
    title = "The two axes as one number per ratio, formed the two possible ways",
    subtitle = paste0("EXPLORATORY - not pre-registered | adjusted for ",
                      "proliferation | ", nrow(RATIO_GRID), " ratios\n",
                      "rows are the pro-apoptotic numerator, columns the ",
                      "anti-apoptotic denominator"),
    x = "anti-apoptotic (denominator)", y = "pro-apoptotic (numerator)",
    caption = CAP11) +
  theme_e10 +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.key.width = ggplot2::unit(1.4, "cm"),
                 strip.text.x = ggplot2::element_text(size = 7.5))
.save(g11, "E10_fig8_axis_difference_heatmap", 9.6, 8.2)

g10dat <- component_strata %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                stratum = factor(stratum, levels = STRATA_PRIMING),
                gene = factor(gene, levels = PRIMING_ORDER))
g10 <- ggplot2::ggplot(g10dat, ggplot2::aes(rho, gene, colour = stratum)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                          position = ggplot2::position_dodge(width = 0.7),
                          linewidth = 0.4, alpha = 0.6) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.7),
                      size = 1.7) +
  ggplot2::facet_grid(side ~ cohort + axis, scales = "free_y",
                      space = "free_y") +
  ggplot2::scale_colour_manual(values = STRAT_COLS, name = NULL) +
  ggplot2::labs(
    title = "The 12 priming genes by compartment",
    subtitle = paste("EXPLORATORY - not pre-registered | adjusted for",
                     "proliferation | 95% Fisher-z intervals; Basal is widest"),
    x = "partial Spearman rho with the axis", y = NULL,
    caption = paste0(
      "The intervals are here because the compartments have very different\n",
      "sizes and a Basal point is not the same kind of estimate as a Luminal\n",
      "one. Where the three colours stack, the correlation is a within-subtype\n",
      "property; where `all` sits outside both, it is between-subtype.")) +
  theme_e10
.save(g10, "E10_fig7_priming_components_by_compartment", 11, 6)

# =============================================================================
# 7. Save
# =============================================================================
message("\n7. save")
saveRDS(list(
  additive_fit = additive_fit,
  reannot = reannot, predictors = predictors,
  reactome_modules = REACTOME_MODULES, hit_counts = hit_counts,
  gene_cor = gene_cor, gene_cor_raw = gene_cor_raw,
  canon_wide = canon_wide, axis_agree = axis_agree,
  s6_by_axis = s6_by_axis,
  priming = priming, component_cor = component_cor, coexpr = coexpr,
  gain_summary = gain_summary, ratio_grid = RATIO_GRID,
  priming_strata = priming_strata, component_strata = component_strata,
  priming_marked = priming_marked, strata_marked = strata_marked,
  mark_summary = mark_summary, mark_both_list = mark_both_list,
  pooled_two_way = pooled_two_way, gene_gap = gene_gap, gene_lean = gene_lean,
  lean_genes = LEAN_GENES, diff_tab = diff_tab,
  between_test = between_test, lum_basal_gap = lum_basal_gap,
  expr_rank = expr_rank,
  settings = list(priming_pro = PRIMING_PRO, priming_anti = PRIMING_ANTI,
                  low_expr_pct = LOW_EXPR_PCT, msigdb_version = MSIGDB_VERSION,
                  strata = STRATA_PRIMING, min_stratum_n = MIN_STRATUM_N,
                  mark_min_abs_rho = MARK_MIN_ABS_RHO,
                  mark2_min_abs_gap = MARK2_MIN_ABS_GAP,
                  gene_scale = "log2(linear DESeq2-normalised + 1)",
                  axis_scale = "GSVA as built by E02",
                  measure = "spearman", covariate = PROLIF_COV,
                  myc_axis = MYC_REF, seed = PROJECT_SEED),
  rules = list(
    adjustment = paste("every value plotted by this script is a PARTIAL",
                       "Spearman on", PROLIF_COV, "- added 2026-09-02 so that",
                       "its panels sit on the same footing as E11's and E13's.",
                       "Two figures in one paper cannot be on different",
                       "footings. `gene_cor_raw` keeps the unadjusted pass,",
                       "which is what section 4 asserts against E08 - E08 has",
                       "no covariate to reproduce."),
    measure = paste("Spearman throughout. The Pearson panels this script",
                    "once carried were removed on 2026-09-02: E09 showed the",
                    "two measures correlate at 0.996 over 220 pairs, bicor",
                    "sides with Spearman in all twelve of the largest",
                    "disagreements, and the mean departure on GSVA is 0.009.",
                    "The one instrument where the measure does matter is",
                    "mitoPPS (0.029, up to 0.093) and a Pearson against a",
                    "mitoPPS score is not reported anywhere in this study."),
    independence = paste("Reactome supplies the module and MitoCarta 3.0 the",
                         "localisation. Neither saw the CDC curation, which is",
                         "the whole point of the re-annotation. They answer",
                         "different questions - what a gene does versus where",
                         "it is - and APAF1 is where they come apart."),
    self_overlap = paste("CYCS is 1 of the 89 genes in the OXPHOS subunits arm,",
                         "so its OXPHOS column is partly self-correlation.",
                         "FELSHER__MITOSTRIP contains none of the 44 or the 12,",
                         "so the MYC panels have no self-overlap;",
                         "M_b__MITOSTRIP contains seven and is tabled only."),
    bcl2l2 = paste("BCL2L2 was originally supplied on BOTH priming lists. The",
                   "first run flagged it as canonically anti-apoptotic - Bcl-w,",
                   "a multidomain guardian rather than a BH3-only sensitiser -",
                   "and the author confirmed on 2026-09-02 that the pro-side",
                   "entry was an error. It is anti-apoptotic only, so the grid",
                   "is 7 x 5 = 35 with no self-pair."),
    ratios = paste("A ratio is only worth reporting as a ratio where `gain` -",
                   "|rho of the ratio| minus the stronger component's |rho| -",
                   "is positive in BOTH cohorts. Otherwise it is a single gene",
                   "wearing a ratio's name."),
    marks = paste("the * on figures 3 and 6 is `gain > 0` AND |rho| >=",
                  MARK_MIN_ABS_RHO, "- figure 4's test plus a HAND-DRAWN",
                  "effect-size floor that is not a test and moves the marks if",
                  "it is moved. |rho| rather than rho, so a ratio running the",
                  "other way is marked too. The heavy border is the same ratio",
                  "passing in BOTH cohorts at the same stratum and axis, and",
                  "is the only version of the mark the `ratios` rule above",
                  "permits anyone to quote. All 22 marked cells on figure 6",
                  "are OXPHOS; the three bordered ones are all Basal; the",
                  "pooled heatmap has none."),
    marks2 = paste("** upgrades a starred cell whose two genes are BOTH",
                   "OXPHOS-led on their own: |rho(OXPHOS) - rho(MYC)| >=",
                   MARK2_MIN_ABS_GAP, "in at least one cohort, which is E15",
                   "figure 5's quantity recomputed from `component_cor` (E11's",
                   "pipeline misses BIK and BCL2L2; on the 10 genes both cover",
                   "they agree to 0.000). The threshold is on the MAGNITUDE by",
                   "the author's choice on 2026-09-03, which admits MCL1 at",
                   "-0.32 / -0.25. That is the productive configuration, not a",
                   "loophole - a ratio beats its parts only when they move",
                   "oppositely - and under the signed rule NOT ONE marked cell",
                   "qualifies. `mark2_signed` keeps that count. All 7 ** cells",
                   "on figure 6 have MCL1 as denominator, and BBC3/MCL1 in",
                   "Basal is the only cell carrying all three marks."),
    difference = paste("figure 8 collapses the two axes into one number per",
                       "ratio, and there are two ways to do it that disagree",
                       "in sign on 22 of 70 cells. SIGNED,",
                       "rho(OXPHOS) - rho(MYC), keeps direction and therefore",
                       "scatters the OXPHOS-dominated cells across both ends",
                       "of a diverging scale - BCL2L11/BCL2L1 in TCGA is",
                       "-0.11 on MYC and -0.52 on OXPHOS and scores -0.41.",
                       "WHICH-AXIS, |rho(OXPHOS)| - |rho(MYC)|, is the",
                       "difference of the magnitudes and its sign means only",
                       "which axis governs, so it is the panel to read for",
                       "that question. Direction stays in figure 3. On 13",
                       "cells the axes point opposite ways and the signed",
                       "difference ADDS them."),
    strata = paste("A POOLED value outside the range of its own two",
                   "compartments is a between-subtype effect, not a within-",
                   "subtype one; that is what D3/S1 found for BCL2 against",
                   "MYC. Basal is 171 TCGA samples and its intervals are",
                   "about +/- 0.15 wide."),
    selection = paste("no cell of this grid is a finding on its own. 44 x 6 x",
                      "2 gene cells plus 39 x 6 x 2 ratio cells; read",
                      "structure and cross-cohort agreement, never one",
                      "interval.")),
  built = Sys.time()), PATH_E10)
readr::write_csv(canon_wide, PATH_E10_CANON)
readr::write_csv(priming, PATH_E10_PRIME)

message("\nE10: done.")
message("    results/machinery_and_priming.rds")
message("    outputs/tables/E10_canonical_machinery.csv")
message("    outputs/tables/E10_priming_ratios.csv")
message("    8 figures in outputs/figures/:")
message("      fig1 the 44 vs OXPHOS - the E08 fig6 panel, re-derived")
message("      fig2 the 44 vs MYC - the new axis, same row order")
message("      fig3 the ", nrow(RATIO_GRID), "-cell priming-ratio heatmap")
message("      fig4 whether a ratio beats its stronger component")
message("      fig5 the 12 priming genes before any ratio is taken")
message("      fig6 the ratio heatmap again, split by luminal and basal")
message("           figs 3 and 6 carry figure 4's test as a * and a border,")
message("           and ** where both genes are OXPHOS-led (^ on the axes)")
message("      fig7 the 12 genes by compartment, with intervals")
message("      fig8 the two axes as ONE number, formed the two possible ways")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E10)

  # Q-a: did S6 survive an independent annotation?
  x$predictors %>% as.data.frame()
  x$reannot %>%
    dplyr::select(gene, mean_rho, effect, cdc_acts_at_mito, mitocarta, submito,
                  cdc_module, reactome_module) %>%
    dplyr::arrange(mean_rho) %>% as.data.frame()

  # the genes Reactome's apoptosis annotation does not carry at all
  x$reannot %>% dplyr::filter(reactome_hits == 0) %>%
    dplyr::select(gene, mean_rho, mitocarta, submito) %>% as.data.frame()

  # Q-b: MYC against OXPHOS for the same 44
  x$canon_wide %>%
    dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
    dplyr::select(axis, gene, mean_rho) %>%
    tidyr::pivot_wider(names_from = axis, values_from = mean_rho) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::arrange(OXPHOS) %>% as.data.frame()

  # Q-c: which ratios actually beat their parts
  x$gain_summary %>% dplyr::filter(both_positive) %>% as.data.frame()

  # and the heatmap's numbers, as a table
  x$priming %>%
    dplyr::filter(axis == "OXPHOS") %>%
    dplyr::select(cohort, ratio, rho, rho_pro, rho_anti, coexpr, gain) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

}
