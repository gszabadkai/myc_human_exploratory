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
#   Q-c (author) THE SAME TWO FIGURES UNDER PEARSON, and the BCL2-family
#       PRIMING RATIOS - every pro-apoptotic over every anti-apoptotic gene -
#       against MYC and OXPHOS under both measures, as a heatmap.
#
# =============================================================================
# WHAT MAKES THE PEARSON COMPARISON MEANINGFUL, AND IT IS A SCALE DECISION
# =============================================================================
# Spearman is invariant to monotone transforms; Pearson is not. So "the Pearson
# version of this figure" is undefined until the scale is stated. E09 made the
# same decision and it is repeated here:
#
#   GENE LEVEL   log2(linear DESeq2-normalised + 1). One matrix, both measures.
#                The +1 is a pseudocount, not a filter - every gene in this
#                script sits above the 25th expression percentile in both
#                cohorts (checked in section 2), so the offset is negligible
#                and no rho depends on it.
#   AXIS LEVEL   the GSVA scores as E02 built them. GSVA output is already a
#                roughly symmetric per-sample statistic, which is the scale
#                Pearson assumes.
#
# BECAUSE log2(x + 1) IS MONOTONE, THE SPEARMAN VALUES COMPUTED HERE ARE THE
# SAME NUMBERS E08 COMPUTED ON THE RAW LINEAR MATRIX. Section 4 asserts that
# against `canonical$mean_rho` rather than assuming it.
#
# mitoPPS is carried in the tables but NOT plotted under Pearson. It is a
# composition statistic on a linear scale with a long right tail (CLAUDE.md
# trap 6), and a Pearson on it would be measuring the tail.
#
# =============================================================================
# THE PRIMING RATIOS, AND THE ONE THING TO KNOW BEFORE READING THEM
# =============================================================================
# The author's lists, verbatim:
#   pro   BCL2L11 (BIM), BMF, PMAIP1 (NOXA), BBC3 (PUMA), BCL2L2, BID, BAD, BIK
#   anti  BCL2, BCL2L1 (Bcl-XL), MCL1, BCL2L2 (Bcl-w), BCL2A1 (A1)
#
# BCL2L2 APPEARS ON BOTH SIDES. It is canonically ANTI-apoptotic - Bcl-w, a
# multidomain guardian, not a BH3-only sensitiser - so its presence on the pro
# side is almost certainly a slip for another BH3-only gene. It is kept on both
# sides exactly as specified, because silently deleting an input is worse than
# reporting it; the degenerate BCL2L2/BCL2L2 cell is the only one dropped.
# 8 x 5 - 1 = 39 ratios. If the pro-side entry was meant to be HRK or BOK, both
# are already scored in section 4 and the ratio grid can be rebuilt in minutes.
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
# SCALE: log2(linear + 1) at gene level; GSVA as built at axis level.
# SPECIES: human, natively. No ortholog function is called.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE10: the canonical machinery, both axes, both measures, and the ",
        "priming ratios\n", strrep("=", 78))

PATH_E10       <- file.path(DIR_RESULTS, "machinery_and_priming.rds")
PATH_E10_CANON <- file.path(DIR_TABLES,  "E10_canonical_machinery.csv")
PATH_E10_PRIME <- file.path(DIR_TABLES,  "E10_priming_ratios.csv")

# The author's lists, verbatim. See the header for BCL2L2.
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BCL2L2", "BID", "BAD",
                  "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")

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
# THE ONE MATRIX BOTH MEASURES ARE COMPUTED FROM. log2(x + 1) is monotone, so
# every Spearman below is identical to the one E08 computed on the raw linear
# matrix - section 4 asserts it - while Pearson now has a stated scale.
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

COH <- list(
  TCGA     = list(G = GT, res = RES_T, ax = AX_T, ids = ID_T),
  `SCAN-B` = list(G = GS, res = RES_S, ax = AX_S, ids = ID_S))
message("   TCGA ", length(ID_T), " samples | SCAN-B ", length(ID_S),
        " samples | ", nrow(AX_T), " axes")

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

# Spearman and Pearson of every AXIS against every ITEM, in one call each.
# Fisher-z intervals: Bonett-Wright for Spearman, se = sqrt((1 + rho^2/2)/(n-3)),
# which is the same variance the atlas engine uses; the plain 1/(n-3) is the
# Pearson case and understates a rank correlation.
#
# THE INTERVAL IS FOR DESCRIBING ONE CELL, NEVER FOR SELECTING ONE. This script
# emits 44 x 6 x 2 x 2 gene cells plus 39 x 6 x 2 x 2 ratio cells; a cell whose
# interval excludes zero is not a finding.
.cor_block <- function(A, B, method) {
  n <- ncol(A)
  stopifnot(identical(colnames(A), colnames(B)))
  R  <- suppressWarnings(stats::cor(t(A), t(B), method = method))
  z  <- atanh(pmin(pmax(R, -0.999999999), 0.999999999))
  se <- if (method == "spearman") sqrt((1 + R^2 / 2) / (n - 3)) else
        matrix(1 / sqrt(n - 3), nrow(R), ncol(R))
  tibble::tibble(
    axis   = rep(rownames(R), times = ncol(R)),
    item   = rep(colnames(R), each  = nrow(R)),
    method = method,
    n      = n,
    rho    = as.vector(R),
    ci_lo  = as.vector(tanh(z - 1.959964 * se)),
    ci_hi  = as.vector(tanh(z + 1.959964 * se)))
}
.both_measures <- function(A, B) {
  dplyr::bind_rows(.cor_block(A, B, "spearman"), .cor_block(A, B, "pearson"))
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
# 4. Q-b and Q-c: the 44 genes on both axes, under both measures
# =============================================================================
message("\n4. the canonical machinery, both axes, both measures")

gene_cor <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  .both_measures(C$ax, GR[[coh]]$mat[CANON_GENES, , drop = FALSE]) %>%
    dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(gene = item)

# --- THE REPRODUCTION CHECK --------------------------------------------------
# E08 computed these Spearman values on the RAW LINEAR matrix through the atlas
# engine. This script computes them on log2(linear + 1) through stats::cor.
# log2(x + 1) is monotone, so the two must agree to floating point. If they do
# not, one of the two scripts is not reading the plane it says it is, and every
# number below is suspect.
check <- gene_cor %>%
  dplyr::filter(axis == "OXPHOS", method == "spearman") %>%
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
  dplyr::select(cohort, axis, method, gene, rho) %>%
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
  dplyr::filter(axis == "MYC", method == "spearman") %>%
  dplyr::arrange(mean_rho) %>%
  dplyr::mutate(dplyr::across(c(TCGA, `SCAN-B`, mean_rho), ~ round(.x, 3))) %>%
  dplyr::select(gene, TCGA, `SCAN-B`, mean_rho, agree, effect,
                cdc_acts_at_mito, cdc_module) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   do the two axes agree about these 44 genes?")
axis_agree <- canon_wide %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(method, axis, gene, mean_rho) %>%
  tidyr::pivot_wider(names_from = axis, values_from = mean_rho) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(spearman_of_the_two_axes =
                     stats::cor(MYC, OXPHOS, method = "spearman"),
                   median_MYC = stats::median(MYC),
                   median_OXPHOS = stats::median(OXPHOS), .groups = "drop")
axis_agree %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   does the S6 split hold on the MYC axis too?")
s6_by_axis <- canon_wide %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::group_by(axis, method) %>%
  dplyr::summarise(
    vs_mitocarta = stats::cor(mean_rho, as.numeric(mitocarta),
                              method = "spearman"),
    vs_pro_death = stats::cor(mean_rho, as.numeric(effect == "pro-death"),
                              method = "spearman"),
    median_mito  = stats::median(mean_rho[mitocarta]),
    median_nonmito = stats::median(mean_rho[!mitocarta]), .groups = "drop")
s6_by_axis %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   where do SPEARMAN and PEARSON disagree most over the 44?")
measure_gap <- canon_wide %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(axis, method, gene, TCGA, `SCAN-B`) %>%
  tidyr::pivot_longer(c(TCGA, `SCAN-B`), names_to = "cohort",
                      values_to = "rho") %>%
  tidyr::pivot_wider(names_from = method, values_from = rho) %>%
  dplyr::mutate(gap = pearson - spearman)
measure_gap %>% dplyr::arrange(dplyr::desc(abs(gap))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   Pearson ABOVE Spearman means a minority of samples with large",
        " expression\n   carry the association; BELOW means the rank",
        " relationship is the robust one.\n   E09 is the script that settles",
        " which reading to trust, with bicor as the arbiter.")

# =============================================================================
# 5. Q-c: the BCL2-family priming ratios
# =============================================================================
# On log2 expression a ratio is a DIFFERENCE, so log2(pro + 1) - log2(anti + 1)
# is the log priming ratio and needs no separate arithmetic. Every pro x anti
# combination is formed; BCL2L2/BCL2L2 is the only cell dropped, because a gene
# over itself is the constant zero.
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
if (length(intersect(PRIMING_PRO, PRIMING_ANTI))) {
  message("   NOTE: ", paste(intersect(PRIMING_PRO, PRIMING_ANTI),
                             collapse = ", "),
          " appears on BOTH sides of the author's lists.\n",
          "   BCL2L2 is Bcl-w, canonically ANTI-apoptotic. It is kept exactly",
          " as given;\n   only the degenerate gene-over-itself cell is dropped.")
}

RATIO_GRID <- tidyr::expand_grid(pro = PRIMING_PRO, anti = PRIMING_ANTI) %>%
  dplyr::filter(pro != anti) %>%
  dplyr::mutate(ratio = paste0(pro, "/", anti))
message("   ", nrow(RATIO_GRID), " ratios (", length(PRIMING_PRO), " x ",
        length(PRIMING_ANTI), " minus the self-pair)")

.ratio_matrix <- function(M) {
  R <- M[RATIO_GRID$pro, , drop = FALSE] - M[RATIO_GRID$anti, , drop = FALSE]
  rownames(R) <- RATIO_GRID$ratio
  R
}
ratio_cor <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  M <- GR[[coh]]$mat
  .both_measures(C$ax, .ratio_matrix(M)) %>% dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(ratio = item)

# The components, and the co-expression that decides whether a ratio can add
# anything at all.
component_cor <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  .both_measures(C$ax, GR[[coh]]$mat[PRIMING_ALL, , drop = FALSE]) %>%
    dplyr::mutate(cohort = coh)
})) %>% dplyr::rename(gene = item) %>%
  dplyr::left_join(dplyr::select(expr_rank, cohort, gene, expr_decile,
                                 low_expression), by = c("cohort", "gene")) %>%
  dplyr::mutate(side = dplyr::case_when(
    gene %in% PRIMING_PRO & gene %in% PRIMING_ANTI ~ "listed on both sides",
    gene %in% PRIMING_PRO ~ "pro-apoptotic (author's list)",
    TRUE ~ "anti-apoptotic (author's list)"))

coexpr <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M <- GR[[coh]]$mat
  dplyr::bind_rows(lapply(c("spearman", "pearson"), function(meth) {
    v <- vapply(seq_len(nrow(RATIO_GRID)), function(i)
      stats::cor(M[RATIO_GRID$pro[i], ], M[RATIO_GRID$anti[i], ],
                 method = meth), numeric(1))
    RATIO_GRID %>% dplyr::mutate(cohort = coh, method = meth, coexpr = v)
  }))
}))

priming <- ratio_cor %>%
  dplyr::left_join(RATIO_GRID, by = "ratio") %>%
  dplyr::left_join(dplyr::select(component_cor, cohort, method, axis,
                                 pro = gene, rho_pro = rho),
                   by = c("cohort", "method", "axis", "pro")) %>%
  dplyr::left_join(dplyr::select(component_cor, cohort, method, axis,
                                 anti = gene, rho_anti = rho),
                   by = c("cohort", "method", "axis", "anti")) %>%
  dplyr::left_join(coexpr, by = c("cohort", "method", "pro", "anti", "ratio")) %>%
  dplyr::mutate(best_component = pmax(abs(rho_pro), abs(rho_anti)),
                gain = abs(rho) - best_component) %>%
  dplyr::select(cohort, method, axis, ratio, pro, anti, n, rho, ci_lo, ci_hi,
                rho_pro, rho_anti, coexpr, best_component, gain)

message("\n   strongest ratios against OXPHOS subunits (Spearman, both cohorts",
        " agreeing in sign):")
priming %>%
  dplyr::filter(axis == "OXPHOS", method == "spearman") %>%
  dplyr::select(cohort, ratio, rho, gain) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = c(rho, gain)) %>%
  dplyr::filter(sign(rho_TCGA) == sign(`rho_SCAN-B`)) %>%
  dplyr::mutate(mean_rho = (rho_TCGA + `rho_SCAN-B`) / 2) %>%
  dplyr::arrange(dplyr::desc(abs(mean_rho))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the same against MYC (", MYC_REF, "):")
priming %>%
  dplyr::filter(axis == "MYC", method == "spearman") %>%
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
  dplyr::select(cohort, method, axis, ratio, gain) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = gain) %>%
  dplyr::mutate(both_positive = TCGA > 0 & `SCAN-B` > 0)
gain_summary %>% dplyr::group_by(axis, method) %>%
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

message("\n   the components on their own (Spearman, mean over cohorts):")
component_cor %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"), method == "spearman") %>%
  dplyr::group_by(gene, side, axis) %>%
  dplyr::summarise(mean_rho = round(mean(rho), 3), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = axis, values_from = mean_rho) %>%
  dplyr::arrange(side, dplyr::desc(OXPHOS)) %>%
  as.data.frame() %>% print(row.names = FALSE)

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
  dplyr::filter(axis == "OXPHOS", method == "spearman") %>%
  dplyr::select(gene, order_rho = mean_rho)
MODULE_ORDER <- ORDER_KEY %>%
  dplyr::left_join(dplyr::select(reannot, gene, cdc_module), by = "gene") %>%
  dplyr::group_by(cdc_module) %>%
  dplyr::summarise(m = stats::median(order_rho), .groups = "drop") %>%
  dplyr::arrange(m) %>% dplyr::pull(cdc_module)
GENE_ORDER <- ORDER_KEY %>% dplyr::arrange(order_rho) %>% dplyr::pull(gene)

.machinery_plot <- function(ax, meth, xlab, title, caption, flag_genes,
                            flag_label) {
  d <- canon_wide %>%
    dplyr::filter(axis == ax, method == meth) %>%
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
  sub <- paste0("EXPLORATORY - not pre-registered | circle = in MitoCarta 3.0",
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
  "order are FIXED to the OXPHOS Spearman panel in all four figures, so they can be laid\n",
  "side by side and read for movement. A STARRED gene sits below the 25th expression\n",
  "percentile in at least one cohort: its rho is largely a correlation of quantisation\n",
  "noise with a score and must not be read as a result.\n")

.save(.machinery_plot(
  "OXPHOS", "spearman",
  "mean Spearman rho with OXPHOS subunits (GSVA), across cohorts",
  "The canonical machinery against OXPHOS - Spearman",
  paste0(CAP_MITO,
    "The cross is CYCS, 1 of the 89 genes in the OXPHOS subunits arm itself, so its value\n",
    "is partly a correlation with itself. This is the E08 fig6 panel, recomputed here and\n",
    "asserted identical to it."),
  FLAG_OXPHOS, "in the OXPHOS subunits arm itself"),
  "E10_fig1_machinery_oxphos_spearman", 9, 9)

.save(.machinery_plot(
  "MYC", "spearman",
  paste0("mean Spearman rho with MYC activity (", MYC_REF,
         ", GSVA), across cohorts"),
  "The canonical machinery against MYC - Spearman",
  paste0(CAP_MITO,
    "FELSHER__MITOSTRIP contains none of these 44 genes, so there is no self-overlap on\n",
    "this axis. The crosses are the two that are HALLMARK E2F/G2M members - BIRC5 at the\n",
    "top of the range is survivin, and a MYC signature is entangled with proliferation by\n",
    "construction. This is NOT the OXPHOS panel re-drawn; the rows are in the same order,\n",
    "so compare the two gene by gene."),
  FLAG_MYC, "a HALLMARK E2F/G2M proliferation gene"),
  "E10_fig2_machinery_myc_spearman", 9, 9)

.save(.machinery_plot(
  "MYC", "pearson",
  paste0("mean Pearson r with MYC activity (", MYC_REF, ", GSVA), across cohorts"),
  "The canonical machinery against MYC - Pearson on log2 expression",
  paste0(CAP_MITO,
    "Pearson is NOT scale-free: this is log2(DESeq2-normalised + 1) against the GSVA\n",
    "score, and on any other gene-level scale it would be a different number. Where this\n",
    "differs from the Spearman panel, samples with extreme expression are pulling that\n",
    "gene - E09 is the script that adjudicates, with bicor as the arbiter."),
  FLAG_MYC, "a HALLMARK E2F/G2M proliferation gene"),
  "E10_fig3_machinery_myc_pearson", 9, 9)

.save(.machinery_plot(
  "OXPHOS", "pearson",
  "mean Pearson r with OXPHOS subunits (GSVA), across cohorts",
  "The canonical machinery against OXPHOS - Pearson on log2 expression",
  paste0(CAP_MITO,
    "Pearson is NOT scale-free: this is log2(DESeq2-normalised + 1) against the GSVA\n",
    "score. The cross is CYCS, inside the OXPHOS arm itself. If the localisation split\n",
    "holds under BOTH measures it is not a rank artefact; if it holds under only one,\n",
    "that is the finding."),
  FLAG_OXPHOS, "in the OXPHOS subunits arm itself"),
  "E10_fig4_machinery_oxphos_pearson", 9, 9)

# The measure disagreement, in one plot rather than four comparisons by eye.
g5dat <- canon_wide %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(axis, method, gene, TCGA, `SCAN-B`, mitocarta, effect) %>%
  tidyr::pivot_longer(c(TCGA, `SCAN-B`), names_to = "cohort",
                      values_to = "rho") %>%
  tidyr::pivot_wider(names_from = method, values_from = rho) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
g5 <- ggplot2::ggplot(g5dat, ggplot2::aes(spearman, pearson)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.2) +
  ggplot2::geom_point(ggplot2::aes(colour = effect, shape = mitocarta),
                      size = 1.7) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(g5dat, abs(pearson - spearman) > 0.045),
    ggplot2::aes(label = gene), size = 2.3, max.overlaps = 25,
    show.legend = FALSE, seed = PROJECT_SEED) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_colour_manual(values = EFFECT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                              name = "in MitoCarta 3.0") +
  ggplot2::labs(
    title = "Does the measure change the answer for the apoptotic machinery?",
    subtitle = paste("EXPLORATORY - not pre-registered | 44 genes x 2 axes x 2",
                     "cohorts; Pearson on log2(normalised + 1)"),
    x = "Spearman rho", y = "Pearson r",
    caption = paste0(
      "Points on the diagonal mean the choice of measure is immaterial for that\n",
      "gene. Labelled genes differ by more than 0.045. Pearson ABOVE the diagonal\n",
      "means a minority of high-expressing samples carry the association; BELOW\n",
      "means the rank relationship is the robust one.")) +
  theme_e10
.save(g5, "E10_fig5_spearman_vs_pearson", 8, 7)

# --- the priming ratios ------------------------------------------------------
g6dat <- priming %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(pro = factor(pro, levels = rev(PRIMING_PRO)),
                anti = factor(anti, levels = PRIMING_ANTI),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                method = factor(method, levels = c("spearman", "pearson")))
LIM <- max(abs(g6dat$rho))
g6 <- ggplot2::ggplot(g6dat, ggplot2::aes(anti, pro, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.1) +
  ggplot2::facet_grid(method ~ cohort + axis) +
  ggplot2::scale_fill_gradient2(low = "#2c7bb6", mid = "grey96",
                                high = "#d7191c", midpoint = 0,
                                limits = c(-LIM, LIM), na.value = "grey85",
                                name = "correlation of log2(pro/anti)") +
  ggplot2::labs(
    title = "Every BCL2-family priming ratio against MYC and against OXPHOS",
    subtitle = paste("EXPLORATORY - not pre-registered |", nrow(RATIO_GRID),
                     "ratios; rows are the pro-apoptotic numerator, columns the",
                     "anti-apoptotic denominator"),
    x = "anti-apoptotic (denominator)", y = "pro-apoptotic (numerator)",
    # paste0 does not insert spaces, so a rendered caption line may be split
    # across two source strings to keep the file inside 80 columns.
    caption = paste0(
      "BCL2L2 appears on both of the author's lists; it is Bcl-w and canonically ",
      "anti-apoptotic. The empty cell is BCL2L2/BCL2L2.\n",
      "RED means the ratio rises with the axis - a MORE PRIMED transcriptome in ",
      "high-MYC or high-OXPHOS tumours. A COLUMN that is\n",
      "uniformly coloured is the denominator gene talking rather than priming: ",
      "check `gain` in the table before reading any cell\n",
      "as a ratio. Only 4 to 6 of the 39 ratios beat both of their component ",
      "genes in both cohorts.")) +
  theme_e10 +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.key.width = ggplot2::unit(1.4, "cm"))
.save(g6, "E10_fig6_priming_ratio_heatmap", 13, 7)

# Does the ratio beat its parts? The falsifier, drawn.
g7dat <- priming %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                abs_rho = abs(rho))
g7 <- ggplot2::ggplot(g7dat, ggplot2::aes(best_component, abs_rho,
                                          colour = method)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  ggplot2::geom_point(size = 1.4, alpha = 0.75) +
  ggrepel::geom_text_repel(data = dplyr::filter(g7dat, gain > 0.045),
                           ggplot2::aes(label = ratio), size = 2.1,
                           max.overlaps = 12, min.segment.length = 0.2,
                           show.legend = FALSE, seed = PROJECT_SEED) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_colour_manual(values = c(spearman = "#1b9e77",
                                          pearson = "#7570b3"), name = NULL) +
  ggplot2::labs(
    title = "Does a priming ratio measure more than its stronger half?",
    subtitle = paste("EXPLORATORY - not pre-registered | each point is one of",
                     nrow(RATIO_GRID), "ratios"),
    x = "|rho| of the stronger of the two component genes",
    y = "|rho| of the log2 ratio",
    caption = paste0(
      "ABOVE the diagonal the ratio adds information; ON or BELOW it, the single\n",
      "gene is the better measurement and the ratio is that gene with noise added\n",
      "to it. Only ratios above the line IN BOTH COHORTS are worth reporting as\n",
      "ratios, and the median ratio sits below it on every axis and both measures.")) +
  theme_e10
.save(g7, "E10_fig7_ratio_vs_components", 8, 7)

# The 12 component genes on their own, which is what most of the heatmap is.
# Ordered by the mean OXPHOS Spearman, stated rather than inferred - reorder()
# on a column that mixes two axes would silently average them.
PRIMING_ORDER <- component_cor %>%
  dplyr::filter(axis == "OXPHOS", method == "spearman") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(m = mean(rho), .groups = "drop") %>%
  dplyr::arrange(m) %>% dplyr::pull(gene)
g8dat <- component_cor %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                gene = factor(gene, levels = PRIMING_ORDER))
g8 <- ggplot2::ggplot(g8dat, ggplot2::aes(rho, gene, colour = cohort,
                                          shape = method)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55),
                      size = 1.9) +
  ggplot2::facet_grid(side ~ axis, scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(spearman = 16, pearson = 1),
                              name = NULL) +
  ggplot2::labs(
    title = "The 12 priming genes on their own, before any ratio is taken",
    subtitle = paste("EXPLORATORY - not pre-registered | the author's pro and",
                     "anti lists; BCL2L2 appears on both"),
    x = "correlation with the axis", y = NULL,
    caption = paste0(
      "This is the panel the ratio heatmap is built out of. Where a pro and an\n",
      "anti gene move the SAME way their ratio cancels; where they move oppositely\n",
      "the ratio adds them. Every gene here clears the 25th expression percentile\n",
      "in both cohorts, so none of these cells carries the low-expression caveat.")) +
  theme_e10
.save(g8, "E10_fig8_priming_components", 8, 6)

# =============================================================================
# 7. Save
# =============================================================================
message("\n7. save")
saveRDS(list(
  reannot = reannot, predictors = predictors,
  reactome_modules = REACTOME_MODULES, hit_counts = hit_counts,
  gene_cor = gene_cor, canon_wide = canon_wide, axis_agree = axis_agree,
  s6_by_axis = s6_by_axis, measure_gap = measure_gap,
  priming = priming, component_cor = component_cor, coexpr = coexpr,
  gain_summary = gain_summary, ratio_grid = RATIO_GRID,
  expr_rank = expr_rank,
  settings = list(priming_pro = PRIMING_PRO, priming_anti = PRIMING_ANTI,
                  low_expr_pct = LOW_EXPR_PCT, msigdb_version = MSIGDB_VERSION,
                  gene_scale = "log2(linear DESeq2-normalised + 1)",
                  axis_scale = "GSVA as built by E02",
                  myc_axis = MYC_REF, seed = PROJECT_SEED),
  rules = list(
    scale = paste("Pearson is not scale-free. Every Pearson here is on",
                  "log2(normalised + 1) at gene level against the GSVA score;",
                  "on any other gene scale it is a different number. The",
                  "Spearman values are identical to E08's, asserted in",
                  "section 4."),
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
    bcl2l2 = paste("BCL2L2 is on both of the author's priming lists. It is",
                   "Bcl-w and canonically anti-apoptotic. Kept as given; only",
                   "BCL2L2/BCL2L2 is dropped."),
    ratios = paste("A ratio is only worth reporting as a ratio where `gain` -",
                   "|rho of the ratio| minus the stronger component's |rho| -",
                   "is positive in BOTH cohorts. Otherwise it is a single gene",
                   "wearing a ratio's name."),
    selection = paste("no cell of this grid is a finding on its own. 44 x 6 x 2",
                      "x 2 gene cells plus 39 x 6 x 2 x 2 ratio cells; read",
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
message("      fig1 the 44 vs OXPHOS, Spearman - the E08 fig6 panel, re-derived")
message("      fig2 the 44 vs MYC, Spearman - the new axis, same row order")
message("      fig3 the 44 vs MYC, Pearson")
message("      fig4 the 44 vs OXPHOS, Pearson")
message("      fig5 Spearman against Pearson for all 44, both axes")
message("      fig6 the ", nrow(RATIO_GRID),
        "-cell priming-ratio heatmap, both axes and both measures")
message("      fig7 whether a ratio beats its stronger component")
message("      fig8 the 12 priming genes before any ratio is taken")

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
    dplyr::filter(method == "spearman", axis %in% c("MYC", "OXPHOS")) %>%
    dplyr::select(axis, gene, mean_rho) %>%
    tidyr::pivot_wider(names_from = axis, values_from = mean_rho) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::arrange(OXPHOS) %>% as.data.frame()

  # Q-c: which ratios actually beat their parts
  x$gain_summary %>% dplyr::filter(both_positive) %>% as.data.frame()

  # and the heatmap's numbers, as a table
  x$priming %>%
    dplyr::filter(axis == "OXPHOS", method == "spearman") %>%
    dplyr::select(cohort, ratio, rho, rho_pro, rho_anti, coexpr, gain) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

}
