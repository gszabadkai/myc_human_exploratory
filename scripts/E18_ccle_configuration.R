# E18_ccle_configuration.R
# =============================================================================
# B4 ADDENDUM - DOES THE TRANSCRIPT CONFIGURATION REPRODUCE IN CELL LINES?
#
# This is an ADDENDUM to B4, not a new gate. It re-fits nothing, it re-scores
# no dependency, it revises no verdict. docs/2026-09-04_b4_declaration.md is
# dated, was correct, and is NOT amended by this script.
#
# =============================================================================
# THE QUESTION B4 NEVER ASKED
# =============================================================================
# B4 tested V1's DEPENDENCY prediction in DepMap and got the reverse on both
# genes, on all four rulers, in breast and pan-cancer.
#
# But B4 ASSUMED, and never checked, that the tumour transcript configuration
# holds in cell lines. Synthesis 3.1, TCGA, partial Spearman on
# PROLIF_DISJOINT:  BCL2L1 +0.388 with OXPHOS, MCL1 -0.266.  V1 reads a
# dependency off that configuration; B4 measured the dependency; nobody
# measured the configuration in the system where the dependency was measured.
#
# That is what this script measures, and it is the difference between two
# different results:
#
#   IF THE CONFIGURATION REPRODUCES (BCL2L1 positive, MCL1 negative, gap
#   positive) - then transcript and dependency genuinely run in OPPOSITE
#   directions in the same cells. Surplus BCL-XL is dispensable and withdrawn
#   MCL1 is load-bearing. That is script 14's own two-sided problem stated as a
#   finding rather than as a caveat: a cell can need a protein more precisely
#   because it has less of it. It is coherent - and it inverts the therapeutic
#   reading back toward MCL1 dependence, which is the field's default.
#
#   IF IT DOES NOT REPRODUCE - then B4 never tested V1. The dependency was
#   measured against an axis that orders these genes differently in culture
#   than in tumours, and the reversal says nothing about the tumour model.
#
# BOTH READINGS ARE LEGITIMATE AND NEITHER IS A RESCUE. The rule that decides
# between them is fixed in section 0b, BEFORE any number is computed, and it is
# not re-decided after.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - no CRISPR endpoint is fitted. The single lm() in section 5 is a
#     REPRODUCTION CONTROL whose only job is to prove these rulers are B4's
#     rulers and not a same-recipe rebuild. Its coefficient is asserted equal
#     to the saved value and is reported as a result NOWHERE.
#   - no MYC stratification, no interaction, no Johnson-Neyman. Still G3's.
#   - no CCLE matched null. B4-a failed; there is no positive arm result to
#     calibrate, and the declaration's 6.1 makes the null conditional on one.
#   - nothing is written to myc_human_validation. It is closed, pre-registered
#     and frozen at d3ac60e.
#
# =============================================================================
# THE ESTIMATOR, AND WHY IT IS A CORRELATION AND NOT A COEFFICIENT
# =============================================================================
# The tumour arm measures this configuration as a PARTIAL SPEARMAN adjusted for
# PROLIF_DISJOINT (E10, E16 check 3). This script uses the same estimator on
# the same 12 genes, adjusted for the CCLE PROLIF score, so the CCLE row and
# the tumour row are the same quantity computed the same way and can be read
# side by side. `.cor_block` below is E16's, copied verbatim.
#
# It is NOT the lm coefficient B4 fitted. B4's endpoint was a Chronos score and
# its exposure the same ruler; this endpoint is a transcript. Using B4's
# machinery here would put a rank measure and a least-squares measure in one
# table for no gain.
#
# =============================================================================
# SCALE, AND COHORT-RELATIVITY
# =============================================================================
# Gene-level endpoints are log2 of the LINEAR matrix. DepMap ships
# OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv as log2(TPM + 1), which is
# already that shape, so the 12 genes are read AS SUPPLIED and nothing is
# re-logged. The rulers are built by section 4 exactly as E17 built them:
# GSVA and the two composites on the log matrix, mitoPPS on 2^x - 1.
#
# GSVA IS COHORT-RELATIVE AND SO IS A Z-SCORE. The breast and pan-cancer
# versions of every table below are DIFFERENT QUANTITIES and are never compared
# numerically - only their signs and their orderings are. The same rule bars
# comparing any CCLE value to a TCGA or SCAN-B value. Sections 6 and 7 print
# the tumour rows beside the CCLE rows for exactly one purpose: to compare the
# SIGN and the ORDERING of a configuration, which is what "reproduces" means
# here and is stated as such in section 0b.
#
# N3: never "primed", of a transcript or of a cell line.
# SPECIES: human. Human MitoCarta, rebuilt from THIS repo's pinned workbook.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
suppressPackageStartupMessages(library(data.table))

message("\nE18: the CCLE configuration check - a B4 addendum\n", strrep("=", 78))

PATH_E17      <- file.path(DIR_RESULTS, "depmap_ox_dependency.rds")
PATH_E16      <- file.path(DIR_RESULTS, "respiratory_rulers.rds")
PATH_E18      <- file.path(DIR_RESULTS, "ccle_configuration.rds")
PATH_E18_CFG  <- file.path(DIR_TABLES,  "E18_ccle_configuration.csv")
PATH_E18_GAP  <- file.path(DIR_TABLES,  "E18_ccle_gap.csv")

if (!file.exists(PATH_E17)) {
  stop("results/depmap_ox_dependency.rds is absent. E18 is an addendum to B4 ",
       "and reads B4's saved object for its reproduction control and for the ",
       "coefficients it sets beside the configuration. Source E17 first.",
       call. = FALSE)
}
if (!file.exists(PATH_E16)) {
  stop("results/respiratory_rulers.rds is absent. E18 reads E16's tumour-side ",
       "configuration to put beside the CCLE one. Source E16 first.",
       call. = FALSE)
}
B4 <- readRDS(PATH_E17)
E16 <- readRDS(PATH_E16)

# =============================================================================
# 0. THE GUARD. Before any read. Verbatim from E17 - see its section 0 for why
#    a broken symlink must STOP rather than degrade to a skip.
# =============================================================================
DIR_DEPMAP <- here::here("data", "raw", "depmap")

DEPMAP_RELEASE <- "26Q1"
EXPR_STRANDED  <- FALSE
EXPR_FILE <- if (EXPR_STRANDED)
  "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv" else
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"

.link_help <- paste0(
  "\n\ndata/raw/depmap is a SYMLINK, not a copy:\n",
  "  data/raw/depmap -> /Users/gs/code/myc_human_validation/data/raw/depmap\n",
  "Recreate it with:\n",
  "  ln -s /Users/gs/code/myc_human_validation/data/raw/depmap ",
  "data/raw/depmap\n",
  "Provenance, checksums and the browser-only acquisition route: ",
  "data/raw/depmap_README.md")

.p <- DIR_DEPMAP
if (!nzchar(Sys.readlink(.p)) && !dir.exists(.p)) {
  stop("data/raw/depmap: symlink missing or unresolved.", .link_help,
       call. = FALSE)
}
if (nzchar(Sys.readlink(.p)) && !dir.exists(Sys.readlink(.p))) {
  stop("data/raw/depmap: symlink target does not resolve -> ",
       Sys.readlink(.p), .link_help, call. = FALSE)
}
for (f in c("Model.csv", EXPR_FILE, "CRISPRGeneEffect.csv")) {
  if (!file.exists(file.path(.p, f))) {
    stop("data/raw/depmap/", f, " unreadable through the symlink.", .link_help,
         call. = FALSE)
  }
}
message("\n0. symlink resolves -> ", Sys.readlink(.p))

PATH_MODEL  <- file.path(DIR_DEPMAP, "Model.csv")
PATH_EXPR   <- file.path(DIR_DEPMAP, EXPR_FILE)
PATH_CRISPR <- file.path(DIR_DEPMAP, "CRISPRGeneEffect.csv")

# KEPT FROM SCRIPT 14 AND FROM E17. The 24Q4 figshare CRISPR matrix would
# resurrect itself SILENTLY under the original name; byte length is the only
# signal that distinguishes it, and nothing inside the file records its release.
CRISPR_24Q4_BYTES <- 428678699
for (p in c(PATH_MODEL, PATH_EXPR, PATH_CRISPR)) {
  # An HTML Cloudflare challenge page saved under a .csv name is ~5 KB.
  if (file.size(p) < 1e5) {
    stop(basename(p), " is only ", file.size(p), " bytes - almost certainly the",
         " Cloudflare verification HTML page saved under a .csv name.",
         .link_help, call. = FALSE)
  }
}
if (identical(as.numeric(file.size(PATH_CRISPR)),
              as.numeric(CRISPR_24Q4_BYTES))) {
  stop("CRISPRGeneEffect.csv is exactly ", CRISPR_24Q4_BYTES, " bytes, which is",
       " the DepMap 24Q4 figshare file, not ", DEPMAP_RELEASE, ".", .link_help,
       call. = FALSE)
}
message("   byte guard passed")

# =============================================================================
# 0b. Constants, and THE READING RULE, fixed before any number
# =============================================================================
LINEAGE      <- "Breast"
ARM_PRIMARY  <- "OXPHOS subunits"
RULERS       <- c("ox_gsva", "ox_ppd", "ox_lvl", "ox_rel")
GSVA_MIN_SET <- 3L
MIN_SET_FRAC <- 0.80
MIN_LINES    <- 25L

# E16's twelve, unchanged. ALL TWELVE ARE REPORTED - taking BCL2L1 and MCL1
# alone after seeing them is the grid-of-cells trap, and the prompt for this
# addendum asks for the whole configuration for exactly that reason.
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
PRIMING_ALL  <- sort(unique(c(PRIMING_PRO, PRIMING_ANTI)))
PRIMARY_GENES <- c("BCL2L1", "MCL1")
# E16's load-bearing quantity: the anti-apoptotic guardian minus the BH3-only
# sensitiser that most closely tracks it on the axis. Positive in 8 of 8
# ruler-by-cohort cells in tumours.
GAP_NUM <- "BCL2L1"; GAP_DEN <- "BBC3"

# --- THE RULE. Fixed here, before section 6 runs. ---------------------------
# "The configuration reproduces" means, ON ONE RULER, all three of:
#     rho(BCL2L1, OX) > 0,  rho(MCL1, OX) < 0,  and gap = BCL2L1 - BBC3 > 0.
# It is a statement about SIGNS AND ORDERING, never about magnitudes, because
# a CCLE partial rho and a TCGA partial rho are not comparable quantities.
# The cohort-level verdict is the count of rulers, of 4, on which all three
# hold. 4 of 4 or 3 of 4 is "reproduces"; 1 or 0 is "does not"; 2 of 4 is
# recorded as SPLIT and is not resolved by a tie-break invented afterwards.
CFG_RULE <- paste0(
  "on one ruler: rho(BCL2L1) > 0 AND rho(MCL1) < 0 AND gap(BCL2L1 - BBC3) > 0",
  ". Cohort verdict: 4/4 or 3/4 = REPRODUCES, 2/4 = SPLIT, <= 1/4 = DOES NOT.")

# =============================================================================
# 1. Gene sets - REBUILT FROM THIS REPO'S PINNED MITOCARTA. E17 section 1,
#    unchanged, including the assertion that the rebuild is set-identical to
#    the snapshot's arm_sets.
# =============================================================================
message("\n1. gene sets, rebuilt from ", basename(PATH_MITOCARTA))

suppressWarnings({
  mc_inv <- readxl::read_xls(PATH_MITOCARTA, sheet = 2)
  mc_pw  <- readxl::read_xls(PATH_MITOCARTA, sheet = 4)
})
# Sheet 4 ends with 5 blank padding rows; leaving them in inflates every pathway
# by 5 phantom "NA" genes, silently. CLAUDE.md read-time trap.
mc_pw <- mc_pw %>% dplyr::filter(!is.na(MitoPathway))
stopifnot(length(unique(mc_inv$Symbol)) == EXPECT_MITOCARTA_ALL,
          nrow(mc_pw) == 149L)

.split_genes <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))
MTDNA_PATHWAY <- "mtDNA-encoded OXPHOS subunits"
mito_raw   <- stats::setNames(lapply(mc_pw$Genes, .split_genes), mc_pw$MitoPathway)
MTDNA_GENES <- grep("^MT-", unique(mc_inv$Symbol), value = TRUE)
stopifnot(length(MTDNA_GENES) == 13L, !MTDNA_PATHWAY %in% names(mito_raw))
mito_paths <- lapply(mito_raw, function(g) setdiff(g, MTDNA_GENES))
mito_paths[[MTDNA_PATHWAY]] <- MTDNA_GENES

ARM_PATHWAYS <- list(
  "OXPHOS subunits"          = "OXPHOS subunits",
  "OXPHOS umbrella"          = "OXPHOS",
  "OXPHOS assembly factors"  = "OXPHOS assembly factors",
  "Mitochondrial ribosome"   = "Mitochondrial ribosome",
  "Nucleotide metabolism"    = "Nucleotide metabolism",
  "ROS and glutathione"      = "ROS and glutathione metabolism",
  "TCA cycle"                = "TCA cycle",
  "Amino acid metabolism"    = "Amino acid metabolism",
  "Lipid metabolism"         = "Lipid metabolism",
  "Fatty acid oxidation"     = c("Fatty acid oxidation", "Carnitine shuttle"),
  "Folate and 1-C"           = "Folate and 1-C metabolism",
  "Glycine metabolism"       = "Glycine metabolism",
  "mtDNA-encoded OXPHOS"     = MTDNA_PATHWAY,
  "CI subunits"              = "CI subunits",
  "CII subunits"             = "CII subunits",
  "CIII subunits"            = "CIII subunits",
  "CIV subunits"             = "CIV subunits",
  "CV subunits"              = "CV subunits")
missing_paths <- setdiff(unlist(ARM_PATHWAYS), names(mito_paths))
if (length(missing_paths)) {
  stop("MitoPathway name(s) not in Human MitoCarta 3.0 sheet 4: ",
       paste(missing_paths, collapse = ", "), call. = FALSE)
}
ARM_SETS <- stats::setNames(
  lapply(ARM_PATHWAYS, function(p) unique(unlist(mito_paths[p], use.names = FALSE))),
  names(ARM_PATHWAYS))

snap_arms <- readRDS(PATH_TCGA_MITO)$arm_sets
same_arm  <- vapply(names(ARM_SETS),
                    function(a) setequal(ARM_SETS[[a]], snap_arms[[a]]), logical(1))
if (!all(same_arm)) {
  stop("arm(s) rebuilt from this repo's MitoCarta differ from the snapshot's ",
       "arm_sets: ", paste(names(same_arm)[!same_arm], collapse = ", "),
       call. = FALSE)
}
message("   18 arms rebuilt, all set-identical to the snapshot's arm_sets")

sd_       <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
MITO_ALL  <- sd_$strip_refs$MITOCARTA_ALL
MYC_SET   <- sd_$myc_sets[[MYC_REF]]           # FELSHER__MITOSTRIP, 61
PROLIF    <- sd_$cov_sets$PROLIF_DISJOINT      # 318
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL,
          length(MYC_SET) == EXPECT_FELSHER_STRIP, length(PROLIF) == 318L)
message("   MYC = ", MYC_REF, " (", length(MYC_SET), "), PROLIF_DISJOINT (",
        length(PROLIF), ")")

# =============================================================================
# 2. Cell lines
# =============================================================================
message("\n2. cell lines")

MODEL <- data.table::fread(PATH_MODEL, data.table = FALSE)
stopifnot("ModelID" %in% names(MODEL))
lin_col <- if ("OncotreeLineage" %in% names(MODEL)) "OncotreeLineage" else
  if ("lineage" %in% names(MODEL)) "lineage" else
    stop("Model.csv has no lineage column (looked for OncotreeLineage, ",
         "lineage). DepMap renamed it; fix `lin_col` rather than guessing.",
         call. = FALSE)
MODEL$lineage <- MODEL[[lin_col]]
breast <- MODEL$ModelID[!is.na(MODEL$lineage) & MODEL$lineage == LINEAGE]
message("   ", length(breast), " ", LINEAGE, " models of ", nrow(MODEL))

# =============================================================================
# 3. Expression. E17's reader, unchanged.
# =============================================================================
message("\n3. expression")

.is_true <- function(v, what) {
  if (is.logical(v)) return(v)
  if (is.numeric(v)) return(v == 1)
  if (is.character(v)) {
    x  <- tolower(trimws(v))
    ok <- c("yes", "true", "t", "1"); no <- c("no", "false", "f", "0", "")
    bad <- setdiff(unique(x), c(ok, no))
    if (length(bad)) {
      stop(what, " has unrecognised value(s): ",
           paste(utils::head(bad, 5), collapse = ", "),
           ". Refusing to guess.", call. = FALSE)
    }
    return(x %in% ok)
  }
  stop(what, " is ", class(v)[1], "; cannot be read as a flag.", call. = FALSE)
}

.read_expression <- function(path) {
  m <- data.table::fread(path, data.table = FALSE)
  if (!"ModelID" %in% names(m)) {
    stop("expression file has no ModelID column. Columns start: ",
         paste(utils::head(names(m), 8), collapse = ", "), call. = FALSE)
  }
  flag_col <- intersect(c("IsDefaultEntryForModel", "is_default_entry"),
                        names(m))[1]
  if (is.na(flag_col)) {
    stop("expression file has no MODEL-level default-entry flag. Without it, ",
         "models with several sequencing runs are counted more than once.",
         call. = FALSE)
  }
  # Deliberately NOT IsDefaultEntryForMC - see E17 section 3.
  keep <- .is_true(m[[flag_col]], flag_col)
  if (!any(keep)) stop("no rows have ", flag_col, " true.", call. = FALSE)
  message("   ", sum(keep), " default rows of ", nrow(m), " (flag: ", flag_col, ")")
  ids <- as.character(m$ModelID[keep])
  if (anyDuplicated(ids)) {
    stop(sum(duplicated(ids)), " ModelID(s) still duplicated after filtering. ",
         "Do not de-duplicate silently - find out why.", call. = FALSE)
  }
  META <- c("V1", "", "ProfileID", "SequencingID", "ModelConditionID", "ModelID",
            "IsDefaultEntryForMC", "IsDefaultEntryForModel", "is_default_entry")
  gene_cols <- setdiff(names(m), META)
  num_ok <- vapply(m[gene_cols], is.numeric, logical(1))
  if (!all(num_ok)) {
    stop("non-numeric column(s) survived the metadata filter: ",
         paste(utils::head(gene_cols[!num_ok], 5), collapse = ", "),
         ". Add them to META in .read_expression().", call. = FALSE)
  }
  M <- as.matrix(m[keep, gene_cols, drop = FALSE])
  rownames(M) <- ids
  colnames(M) <- trimws(sub("\\s*\\(\\d+\\)$", "", colnames(M)))
  dup <- duplicated(colnames(M))
  if (any(dup)) M <- M[, !dup, drop = FALSE]
  M
}

EXPR <- .read_expression(PATH_EXPR)
message("   expression: ", nrow(EXPR), " models x ", ncol(EXPR), " genes")
if (max(EXPR, na.rm = TRUE) > 40) {
  stop("expression matrix max is ", round(max(EXPR, na.rm = TRUE), 1),
       ", which is not log2(TPM+1). Check the file.", call. = FALSE)
}

# The 12 endpoints, read AS SUPPLIED. log2(TPM + 1) is already log of linear.
gene_missing <- setdiff(PRIMING_ALL, colnames(EXPR))
if (length(gene_missing)) {
  stop("BCL2-family gene(s) absent from the CCLE expression matrix: ",
       paste(gene_missing, collapse = ", "),
       ". The configuration cannot be compared to E16's twelve if it is not ",
       "twelve. Do not silently drop one.", call. = FALSE)
}
message("   all ", length(PRIMING_ALL), " BCL2-family genes present")

# The pan-cancer cohort is defined EXACTLY as B4 defined it - lines with both
# expression and CRISPR - so section 7's scores are the scores B4 fitted on and
# not a differently-composed pan-cancer cohort that happens to be about as big.
# Only the id column is read here; the reproduction control in section 5 reads
# the two endpoint columns and nothing else.
crispr_head <- data.table::fread(PATH_CRISPR, nrows = 0L, data.table = FALSE)
crispr_cols <- trimws(sub("\\s*\\(\\d+\\)$", "", names(crispr_head)))
crispr_ids  <- data.table::fread(PATH_CRISPR, select = 1L,
                                 data.table = FALSE)[[1]]
crispr_ids  <- as.character(crispr_ids)
message("   CRISPR: ", length(crispr_ids), " models (id column only)")

# =============================================================================
# 4. Scoring - E17 section 4, verbatim. One run per cohort, never pooled.
# =============================================================================
message("\n4. scoring")

.z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)

.pathway_means <- function(E, sets)
  t(vapply(sets, function(g) colMeans(E[g, , drop = FALSE]), numeric(ncol(E))))
.mitopps_universe <- function(S) {
  N <- ncol(S); P <- nrow(S)
  Bi <- 1 / S
  A  <- (S %*% t(Bi)) / N
  out <- S * (((1 / A) %*% Bi) - Bi) / (P - 1)
  dimnames(out) <- dimnames(S)
  out
}
.check_positive <- function(S, what) {
  bad <- which(!is.finite(S) | S <= 0, arr.ind = TRUE)
  if (nrow(bad)) {
    lab <- utils::head(sprintf("%s / %s", rownames(S)[bad[, 1]],
                               colnames(S)[bad[, 2]]), 6L)
    stop("mitoPPS (", what, "): ", nrow(bad), " arm x line pathway mean(s) are ",
         "zero or non-finite. First few: ", paste(lab, collapse = "; "),
         ". Do NOT add a floor silently.", call. = FALSE)
  }
}

# E16's composite: mean of per-gene z across samples, on a log matrix.
.comp <- function(genes, E) {
  M <- E[intersect(genes, rownames(E)), , drop = FALSE]
  v <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}

.score_cohort <- function(ids, what) {
  message("   scoring ", what, ": ", length(ids), " lines")
  E_LOG <- t(EXPR[ids, , drop = FALSE])
  E_LOG <- E_LOG[stats::complete.cases(E_LOG), , drop = FALSE]

  sets_defined <- c(list(MYC = MYC_SET, PROLIF = PROLIF), ARM_SETS)
  sets_gsva <- lapply(sets_defined, function(s) intersect(s, rownames(E_LOG)))
  cov_tab <- tibble::tibble(scope = what, set = names(sets_defined),
                            n_defined = lengths(sets_defined),
                            n_ccle = lengths(sets_gsva)) %>%
    dplyr::mutate(frac = round(n_ccle / n_defined, 3))
  lost <- cov_tab$set[cov_tab$frac < MIN_SET_FRAC]
  if (length(lost)) {
    stop("gene set(s) below ", MIN_SET_FRAC, " coverage in CCLE (", what, "): ",
         paste(lost, collapse = ", "),
         ". That is a symbol-harmonisation failure, not natural attrition.",
         call. = FALSE)
  }
  tiny <- cov_tab$set[cov_tab$n_ccle < GSVA_MIN_SET]
  if (length(tiny)) {
    stop("gene set(s) with fewer than ", GSVA_MIN_SET, " genes in CCLE: ",
         paste(tiny, collapse = ", "), call. = FALSE)
  }

  GS <- GSVA::gsva(GSVA::gsvaParam(exprData = E_LOG, geneSets = sets_gsva,
                                   kcdf = "Gaussian", minSize = GSVA_MIN_SET,
                                   maxSize = Inf), verbose = FALSE)
  miss <- setdiff(names(sets_gsva), rownames(GS))
  if (length(miss)) {
    stop("GSVA silently dropped ", length(miss), " set(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  E_LIN <- 2^E_LOG - 1
  E_LIN[E_LIN < 0] <- 0                 # floating point on an exact zero only
  S_LIN <- .pathway_means(E_LIN, lapply(ARM_SETS, function(s)
    intersect(s, rownames(E_LOG))))
  .check_positive(S_LIN, what)
  PPS <- .mitopps_universe(S_LIN)

  LVL <- t(vapply(ARM_SETS, function(s) .comp(s, E_LOG), numeric(ncol(E_LOG))))
  REL <- t(vapply(ARM_SETS, function(s)
    .comp(s, E_LOG) - .comp(setdiff(MITO_ALL, s), E_LOG), numeric(ncol(E_LOG))))
  rownames(LVL) <- rownames(REL) <- names(ARM_SETS)

  SC <- list(ox_gsva = t(apply(GS[names(ARM_SETS), , drop = FALSE], 1, .z)),
             ox_ppd  = t(apply(PPS, 1, .z)),
             ox_lvl  = t(apply(LVL, 1, .z)),
             ox_rel  = t(apply(REL, 1, .z)))
  for (r in RULERS) stopifnot(identical(colnames(SC[[r]]), ids))
  list(SC = SC, MYC = .z(GS["MYC", ]), PROLIF = .z(GS["PROLIF", ]),
       coverage = cov_tab, ids = ids)
}

lines_b <- intersect(breast, rownames(EXPR))
if (length(lines_b) < MIN_LINES) {
  stop("only ", length(lines_b), " ", LINEAGE, " lines have expression.",
       call. = FALSE)
}
message("   ", length(lines_b), " ", LINEAGE, " lines with expression - the ",
        "cohort B4 SCORED. B4 FITTED on the ", length(intersect(lines_b, crispr_ids)),
        " of them that also carry CRISPR; this check needs no CRISPR.")

SB <- .score_cohort(lines_b, LINEAGE)

pan_lines <- intersect(rownames(EXPR), crispr_ids)
message("   pan-cancer: ", length(pan_lines), " lines with expression AND ",
        "CRISPR - B4's pan-cancer scoring cohort, re-scored")
SP <- .score_cohort(pan_lines, "pan-cancer")

# =============================================================================
# 5. THE REPRODUCTION CONTROL. Not a result. Not a re-fit of B4.
# =============================================================================
# Two checks, and their only purpose is to establish that the rulers used below
# ARE B4's rulers rather than a same-recipe rebuild that might differ. Without
# this, section 6 would be measuring the configuration against an axis nobody
# had shown to be the axis B4 measured the dependency against - which is the
# very error this addendum exists to rule out.
#
# NO COEFFICIENT FROM THIS SECTION IS REPORTED AS A RESULT ANYWHERE.
message("\n5. reproduction control (not a result)")

# 5.1 The coverage tables must come back identical. This proves the expression
#     matrix, the default-entry filter and the 20 gene sets are the same.
cov_new <- dplyr::bind_rows(SB$coverage, SP$coverage) %>% as.data.frame()
cov_old <- as.data.frame(B4$coverage)
cov_new$scope[cov_new$scope == LINEAGE] <- "Breast"
if (!isTRUE(all.equal(cov_new, cov_old, check.attributes = FALSE))) {
  stop("the CCLE coverage table does not reproduce B4's. The expression ",
       "matrix, the release or a gene set has changed underneath, and the ",
       "rulers below would not be the rulers B4 fitted against.", call. = FALSE)
}
message("   coverage table reproduces B4's, all ", nrow(cov_old), " rows")

# 5.2 One saved coefficient must come back bit-equal. This proves the SCORES,
#     which the coverage table cannot: GSVA, mitoPPS and the two composites all
#     have to land in the same place for this to hold.
repro_gene  <- "MCL1"
repro_ruler <- "ox_gsva"
ix <- match(repro_gene, crispr_cols)
if (is.na(ix)) {
  stop(repro_gene, " is not a column of CRISPRGeneEffect.csv, so the ",
       "reproduction control cannot run. Do not proceed without it.",
       call. = FALSE)
}
Y_repro <- data.table::fread(PATH_CRISPR, select = c(1L, ix),
                             data.table = FALSE)
rownames(Y_repro) <- as.character(Y_repro[[1]])
dep_lines <- intersect(lines_b, rownames(Y_repro))
d_repro <- tibble::tibble(
  MYC    = SB$MYC[dep_lines],
  OX     = SB$SC[[repro_ruler]][ARM_PRIMARY, dep_lines],
  PROLIF = SB$PROLIF[dep_lines],
  Y      = Y_repro[dep_lines, 2L])
got  <- unname(stats::coef(stats::lm(Y ~ MYC + OX + PROLIF, data = d_repro))["OX"])
want <- B4$coefficients %>%
  dplyr::filter(label == "B4 breast", ruler == repro_ruler,
                model == "additive", gene == repro_gene,
                arm == ARM_PRIMARY) %>% dplyr::pull(estimate)
stopifnot(length(want) == 1L)
if (!isTRUE(all.equal(got, want, tolerance = 1e-10))) {
  stop("the reproduction control does NOT reproduce B4: ", repro_gene, " / ",
       repro_ruler, " OX is ", signif(got, 8), " here and ", signif(want, 8),
       " in results/depmap_ox_dependency.rds. The scores have drifted and ",
       "nothing below is comparable to B4.", call. = FALSE)
}
message("   B4's ", repro_gene, " / ", repro_ruler, " OX coefficient reproduces",
        " to 1e-10 (", signif(got, 6), ") on ", length(dep_lines), " lines")
rm(Y_repro)

# =============================================================================
# 6. THE CONFIGURATION. Partial Spearman, 12 genes x 4 rulers.
# =============================================================================
# E16's `.cor_block`, copied verbatim: partial Spearman on ranks with Fisher-z
# intervals and the Bonett-Wright variance se = sqrt((1 + rho^2/2)/(n - 3 - k)).
# The plain 1/(n - 3) is the Pearson case and understates a rank correlation.
# It is copied rather than sourced because E10 and E16 both define it inline and
# there is no shared file; section 6.1's comparison against E16 is what proves
# the copy has not drifted.
message("\n6. the configuration, ", LINEAGE, " (n = ", length(lines_b), ")")

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
    ruler = rep(rownames(R), times = ncol(R)),
    gene  = rep(colnames(R), each  = nrow(R)),
    n = n, k_cov = k,
    rho   = as.vector(R),
    ci_lo = as.vector(tanh(z - 1.959964 * se)),
    ci_hi = as.vector(tanh(z + 1.959964 * se)))
}

# One ruler matrix per cohort: the four rulers on the OXPHOS subunits arm.
.ruler_mat <- function(S) {
  M <- do.call(rbind, lapply(RULERS, function(r) S$SC[[r]][ARM_PRIMARY, S$ids]))
  rownames(M) <- RULERS; colnames(M) <- S$ids
  M
}
.gene_mat <- function(S) {
  M <- t(EXPR[S$ids, PRIMING_ALL, drop = FALSE])
  stopifnot(identical(colnames(M), S$ids))
  M
}

.configuration <- function(S, scope) {
  cfg <- .cor_block(.ruler_mat(S), .gene_mat(S),
                    cov = matrix(S$PROLIF[S$ids], ncol = 1L,
                                 dimnames = list(S$ids, "PROLIF"))) %>%
    dplyr::mutate(scope = scope, adjustment = "PROLIF (CCLE GSVA)",
                  side = dplyr::if_else(gene %in% PRIMING_PRO,
                                        "pro-apoptotic", "anti-apoptotic"))
  cfg
}

cfg_breast <- .configuration(SB, LINEAGE)
cfg_pan    <- .configuration(SP, "pan-cancer")
ccle_cfg   <- dplyr::bind_rows(cfg_breast, cfg_pan)

message("\n   ", LINEAGE, ", partial Spearman on PROLIF, rho:")
cfg_breast %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  dplyr::select(gene, side, ruler, rho) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
  dplyr::arrange(side, gene) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   pan-cancer (a DIFFERENT quantity - signs and ordering only):")
cfg_pan %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  dplyr::select(gene, side, ruler, rho) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
  dplyr::arrange(side, gene) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the two primary genes, with intervals:")
ccle_cfg %>%
  dplyr::filter(gene %in% PRIMARY_GENES) %>%
  dplyr::transmute(scope, gene, ruler, n,
                   rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3),
                   ci_excludes_0 = (ci_lo > 0) | (ci_hi < 0)) %>%
  dplyr::arrange(scope, gene, ruler) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 6.1 THE GAP, E16's load-bearing quantity -------------------------------
# BCL2L1 minus BBC3, per ruler, as a difference of partial rhos - exactly how
# E16's `the_gap` is built. Positive in 8 of 8 ruler-by-cohort cells in tumours.
ccle_gap <- ccle_cfg %>%
  dplyr::filter(gene %in% c(GAP_NUM, GAP_DEN)) %>%
  dplyr::select(scope, ruler, gene, rho) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rho) %>%
  dplyr::mutate(gap = .data[[GAP_NUM]] - .data[[GAP_DEN]])

message("\n6.1 the gap, ", GAP_NUM, " minus ", GAP_DEN, ":")
ccle_gap %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. THE TUMOUR SIDE, read from E16. Signs and ordering only.
# =============================================================================
# CLAUDE.md: never pool GSVA or mitoPPS values across cohorts. CCLE is a third
# cohort and these rows sit beside each other to be compared on SIGN and on
# ORDERING - which is what the section 0b rule asks - and on nothing else.
message("\n7. the tumour side (E16), for sign and ordering only")

tumour_cfg <- E16$twelve %>%
  dplyr::transmute(scope = cohort, gene, ruler, n, rho, ci_lo, ci_hi,
                   adjustment, side)
tumour_gap <- E16$the_gap %>%
  dplyr::transmute(scope = cohort, ruler, BCL2L1, BBC3, gap)

message("\n   the two primary genes, tumours and cell lines, rho:")
dplyr::bind_rows(
  tumour_cfg %>% dplyr::filter(gene %in% PRIMARY_GENES) %>%
    dplyr::transmute(scope, gene, ruler, rho),
  ccle_cfg %>% dplyr::filter(gene %in% PRIMARY_GENES) %>%
    dplyr::transmute(scope = paste0("CCLE ", scope), gene, ruler, rho)) %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
  dplyr::arrange(gene, scope) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the gap, tumours and cell lines:")
dplyr::bind_rows(
  tumour_gap %>% dplyr::transmute(scope, ruler, gap),
  ccle_gap   %>% dplyr::transmute(scope = paste0("CCLE ", scope), ruler, gap)) %>%
  dplyr::mutate(gap = round(gap, 3)) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = gap) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 7.1 the pan-cancer sensitivity ------------------------------------------
# B4's pan-cancer FIT dropped lines in lineages of fewer than 10, leaving 1,127
# of the 1,140 SCORED. This check runs on the scored cohort, which is the right
# one - no lineage term is fitted here - but the 13-line difference is reported
# so it cannot be mistaken for a cohort-definition artefact. The scores are NOT
# rebuilt: this is a subset of the same z-scored vectors.
lin_pan  <- MODEL$lineage[match(pan_lines, MODEL$ModelID)]
keep_lin <- lin_pan %in% names(which(table(lin_pan) >= 10L))
message("\n7.1 sensitivity: the ", sum(keep_lin), " of ", length(pan_lines),
        " lines B4 FITTED, same scores, no re-scoring")
SP_sub <- list(SC = SP$SC, MYC = SP$MYC, PROLIF = SP$PROLIF,
               ids = pan_lines[keep_lin])
cfg_pan_fit <- .configuration(SP_sub, "pan-cancer (B4 fit lines)")
cfg_pan_fit %>%
  dplyr::filter(gene %in% PRIMARY_GENES) %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  dplyr::select(gene, ruler, rho) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. THE READING, on the rule fixed in section 0b
# =============================================================================
# ONE THING THIS SCRIPT DELIBERATELY DOES NOT TEST, because it is already
# answered elsewhere: whether the TUMOUR association is a compartment artefact.
# Cell lines have no stroma or infiltrate, so that is the natural explanation
# for any non-reproduction - and E16 check 4 already refits the twelve on the
# 1,007 TCGA samples with purity and leukocyte fraction on top of PROLIF.
# BCL2L1 / ox_gsva moves +0.419 -> +0.406. The tumour association SURVIVES
# adjustment, so a non-reproduction here is not explained by purity and must
# not be written as if it were. CCLE carries no purity estimate at all.
# =============================================================================
message("\n8. the reading\n", strrep("-", 78))
message("RULE (fixed before any number): ", CFG_RULE)

.verdict <- function(cfg, gap, scope) {
  w <- cfg %>%
    dplyr::filter(gene %in% PRIMARY_GENES) %>%
    dplyr::select(ruler, gene, rho) %>%
    tidyr::pivot_wider(names_from = gene, values_from = rho) %>%
    dplyr::left_join(gap %>% dplyr::select(ruler, gap), by = "ruler") %>%
    dplyr::mutate(bcl2l1_pos = BCL2L1 > 0, mcl1_neg = MCL1 < 0, gap_pos = gap > 0,
                  all_three  = bcl2l1_pos & mcl1_neg & gap_pos)
  k <- sum(w$all_three)
  v <- if (k >= 3L) "REPRODUCES" else if (k == 2L) "SPLIT" else "DOES NOT REPRODUCE"
  list(table = dplyr::mutate(w, scope = scope), k = k, verdict = v)
}

v_breast <- .verdict(cfg_breast, dplyr::filter(ccle_gap, scope == LINEAGE), LINEAGE)
v_pan    <- .verdict(cfg_pan, dplyr::filter(ccle_gap, scope == "pan-cancer"),
                     "pan-cancer")

for (v in list(v_breast, v_pan)) {
  message("\n   ", unique(v$table$scope), ": ", v$k, " of 4 rulers carry all ",
          "three -> ", v$verdict)
  v$table %>%
    dplyr::transmute(ruler, BCL2L1 = round(BCL2L1, 3), MCL1 = round(MCL1, 3),
                     gap = round(gap, 3), bcl2l1_pos, mcl1_neg, gap_pos,
                     all_three) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

reading <- if (v_breast$verdict == "REPRODUCES") {
  paste0(
    "READING A - THE CONFIGURATION REPRODUCES IN BREAST CELL LINES.\n",
    "   Transcript and dependency genuinely run in OPPOSITE directions in the\n",
    "   same cells: the lines that carry more BCL2L1 and less MCL1 transcript\n",
    "   with OXPHOS are the lines that depend LESS on BCL-XL and MORE on MCL1.\n",
    "   That is script 14's two-sided problem as a finding - a cell can need a\n",
    "   protein more precisely because it has less of it - and it INVERTS the\n",
    "   therapeutic reading back toward MCL1 dependence, the field's default.\n",
    "   It is NOT a rescue of V1: V1's dependency prediction still fails.")
} else if (v_breast$verdict == "DOES NOT REPRODUCE") {
  paste0(
    "READING B - THE CONFIGURATION DOES NOT REPRODUCE. B4 NEVER TESTED V1.\n",
    "   The dependency was measured against an axis that orders these genes\n",
    "   differently in culture than in tumours, so B4's reversal says nothing\n",
    "   about the tumour model. It is NOT a rescue either: what it removes is\n",
    "   B4's standing as a test, not V1's exposure to one.")
} else {
  paste0(
    "SPLIT - 2 of 4 rulers carry the configuration and 2 do not.\n",
    "   The section 0b rule does NOT resolve this and no tie-break is invented\n",
    "   here. Report the split, name which rulers fall each way, and treat B4\n",
    "   as a test whose standing depends on a ruler choice - which is exactly\n",
    "   what synthesis 8.4 trap 9 says a ruler result must never be allowed to\n",
    "   do silently.")
}
message("\n", reading)

# =============================================================================
# 9. Save
# =============================================================================
message("\n9. saving")

out <- list(
  ccle_configuration = ccle_cfg,
  ccle_gap           = ccle_gap,
  ccle_pan_fit_lines = cfg_pan_fit,
  tumour_configuration = tumour_cfg,
  tumour_gap           = tumour_gap,
  verdict_breast     = v_breast$table,
  verdict_pan        = v_pan$table,
  verdicts = tibble::tibble(
    scope   = c(LINEAGE, "pan-cancer"),
    n_rulers_of_4 = c(v_breast$k, v_pan$k),
    verdict = c(v_breast$verdict, v_pan$verdict)),
  reading  = reading,
  coverage = dplyr::bind_rows(SB$coverage, SP$coverage),
  lines = list(breast_expression = length(lines_b),
               breast_fitted_by_b4 = length(dep_lines),
               pan_scored = length(pan_lines),
               pan_fitted_by_b4 = sum(keep_lin)),
  repro = list(gene = repro_gene, ruler = repro_ruler,
               got = got, want = want,
               note = paste("a reproduction control only - it proves these",
                            "rulers are B4's rulers. NOT a result and reported",
                            "as one nowhere.")),
  spec = list(
    what = paste("B4 ADDENDUM. Does the tumour transcript configuration",
                 "(synthesis 3.1) reproduce in CCLE? B4 assumed it and never",
                 "checked it."),
    not = paste("re-fits nothing, revises no verdict, amends no declaration.",
                "No MYC stratification, no interaction, no CCLE matched null."),
    rule = CFG_RULE,
    estimator = paste("partial Spearman on the CCLE PROLIF GSVA score,",
                      "Fisher-z intervals, Bonett-Wright variance.",
                      "E16 .cor_block verbatim."),
    endpoints = paste("the 12 BCL2-family genes as log2(TPM+1) AS SUPPLIED -",
                      "already log of linear, nothing is re-logged."),
    cohort_relative = paste("CCLE breast, CCLE pan-cancer, TCGA and SCAN-B are",
                            "four cohorts. Values are NEVER compared across",
                            "them; signs and orderings are."),
    genes = PRIMING_ALL,
    rulers = RULERS,
    releases = c(depmap = DEPMAP_RELEASE, expression = EXPR_FILE),
    depmap_source = paste("symlink data/raw/depmap ->",
                          "/Users/gs/code/myc_human_validation/data/raw/depmap",
                          "- see data/raw/depmap_README.md"),
    n3 = "these are transcripts; the word 'primed' is used of nothing",
    seed = PROJECT_SEED),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E18)
utils::write.csv(as.data.frame(ccle_cfg), PATH_E18_CFG, row.names = FALSE)
utils::write.csv(as.data.frame(ccle_gap), PATH_E18_GAP, row.names = FALSE)
message("   ", PATH_E18)
message("   ", PATH_E18_CFG)
message("   ", PATH_E18_GAP)
message("\nE18 done.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E18)

  # The reading, and the counts behind it.
  x$verdicts %>% as.data.frame()
  cat(x$reading, "\n")
  utils::str(x$lines)

  # The two primary genes, every ruler, both CCLE cohorts.
  x$ccle_configuration %>%
    dplyr::filter(gene %in% c("BCL2L1", "MCL1")) %>%
    dplyr::transmute(scope, gene, ruler, n, rho = round(rho, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3)) %>%
    as.data.frame()

  # The whole configuration - all twelve, so the two are read in context.
  x$ccle_configuration %>%
    dplyr::filter(scope == "Breast") %>%
    dplyr::mutate(rho = round(rho, 3)) %>%
    dplyr::select(gene, side, ruler, rho) %>%
    tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
    dplyr::arrange(side, gene) %>%
    as.data.frame()

  # The gap, cell lines beside tumours. SIGNS AND ORDERING ONLY - these are
  # four different cohorts and the values are not comparable.
  x$ccle_gap %>% as.data.frame()
  x$tumour_gap %>% as.data.frame()

  # The reproduction control. Not a result.
  utils::str(x$repro)

  # The 13 lines B4's pan-cancer fit dropped make no difference:
  dplyr::inner_join(
    x$ccle_configuration %>% dplyr::filter(scope == "pan-cancer") %>%
      dplyr::select(gene, ruler, rho_all = rho),
    x$ccle_pan_fit_lines %>% dplyr::select(gene, ruler, rho_fit = rho),
    by = c("gene", "ruler")) %>%
    dplyr::mutate(delta = round(rho_all - rho_fit, 4)) %>%
    as.data.frame()
}
