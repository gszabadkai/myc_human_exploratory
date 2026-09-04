# E17_depmap_ox_dependency.R
# =============================================================================
# B4 - DOES OXPHOS PREDICT GUARDIAN DEPENDENCY? The OX MAIN EFFECT in DepMap.
#
# Declared in advance: docs/2026-09-04_b4_declaration.md, committed BEFORE this
# script was written. Read it first. Every direction and every reading rule
# below was fixed there and none of them is re-decided here.
#
# =============================================================================
# THE ESTIMAND, AND WHY IT IS NOT BLOCK G REPEATED
# =============================================================================
# The closed validation study's Block G fitted
#
#     lm(Y ~ MYC * OX + PROLIF)        (+ lineage, pan-cancer)
#
# in breast, in the 18-arm panel, pan-cancer, the BCL2L11 check and the drug
# arm, and extracted ONLY `MYC:OX`. Verified 2026-09-04 by grep over the forked
# script: all seven .tidy() calls name "MYC:OX", except P3 which names
# "MYC:OX:brst". Nothing anywhere extracts `OX`.
#
# So the OX MAIN EFFECT was fitted in every one of those models and reported in
# none. That coefficient is what this script asks for.
#
# THE INTERACTION WAS NULL AND STAYS NULL. Block G's own numbers: breast MCL1
# +0.057 / +0.143, BCL2L1 +0.097 / +0.042; pan-cancer at n = 1,130 with CIs of
# about +/- 0.04. NOTHING HERE MAKES AN INTERACTION CLAIM IN EITHER DIRECTION.
# `m_int` below is fitted as a SENSITIVITY and only its `OX` term is read.
#
# THE PREDICTION COMES FROM THE HUMAN TUMOUR ARM, not from DepMap. Synthesis
# 3.1: in OXPHOS-high tumours BCL2L1 transcript is up (+0.388) and MCL1 is down
# (-0.266), while MYC orders neither (+0.002, +0.022). The dependency prediction
# therefore CARRIES NO MYC TERM. MYC is a covariate here, never the exposure,
# and nothing is stratified on it.
#
#   Chronos: 0 = no effect, -1 = median common essential. MORE NEGATIVE = MORE
#   ESSENTIAL. Declared directions:  BCL2L1 -> NEGATIVE OX.  MCL1 -> POSITIVE OX.
#
# THE TWO-SIDED PROBLEM, as script 14 states it: a cell can need a protein more
# precisely because it has less of it. The standard prior runs the other way and
# V1 assumes it. THE REVERSE IS REPORTED WHEN IT OCCURS AND IS NOT A PASS.
#
# =============================================================================
# WHY THE CONTROLS ARE LOAD-BEARING HERE AND WERE NOT FOR THE INTERACTION
# =============================================================================
# An interaction is robust to anything that shifts all dependencies together.
# A MAIN EFFECT IS NOT. If OXPHOS-high lines simply grow differently, `OX` picks
# that up on every gene and `MYC:OX` would not have. Hence:
#   - RPL3 and POLR2B are a FLOOR, not decoration. Signal there voids B4-a.
#   - the full 18-arm panel runs for the primary genes.
#   - PROLIF stays in every model; lineage adjustment is mandatory pan-cancer.
#
# POLR2A IS NOT SCREENED IN 26Q1. Every other RNA Pol II subunit is. The closed
# study hit the same wall and recorded it. DEP_CONTROL is UNCHANGED - POLR2A is
# requested and reported absent - and POLR2B is a SUBSTITUTE FLOOR declared in
# section 5.1 of the declaration, before the fit. It is a floor, not an endpoint.
#
# =============================================================================
# THE DATA IS REACHED BY SYMLINK, AND THE GUARD IS NOT OPTIONAL
# =============================================================================
# data/raw/depmap -> myc_human_validation/data/raw/depmap. See
# data/raw/depmap_README.md, including what this repo can no longer regenerate.
#
# A BROKEN SYMLINK READS AS AN ABSENT DIRECTORY, and script 14's structure
# responds to absence by SKIPPING. Unguarded, this script would report "PRISM
# skipped, files absent" and carry on having silently lost the CRISPR arm too -
# a run that completes and answers a different question. Section 0 therefore
# STOPS on a missing link, an unresolvable target, or any unreadable required
# file. PRISM absence stays a skip; it is genuinely optional.
#
# NOTHING IS WRITTEN TO myc_human_validation. It is closed, pre-registered and
# frozen at d3ac60e. Its results/depmap_dependency.rds is neither read nor
# rewritten.
#
# =============================================================================
# SCALE DISCIPLINE - READ BEFORE EDITING ANY SCORING BLOCK
# =============================================================================
# DepMap OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv is log2(TPM + 1).
#
#   GSVA          -> the file AS SUPPLIED. Log scale, kcdf = "Gaussian".
#   ox_lvl/ox_rel -> the file AS SUPPLIED. E16's composite is mean per-gene z of
#                    log2(linear + 1), and log2(TPM + 1) is that same shape.
#   mitoPPS       -> LINEAR, via 2^x - 1. It cannot see the log matrix.
#
# E_LOG and E_LIN are separate objects and are never rebuilt from one another.
#
# DECLARED DEVIATION, carried forward verbatim from script 14 and NOT re-argued:
# the plan's mitoPPS input is linear DESeq2-normalised counts; here it is linear
# TPM, because that is what DepMap ships. TPM is additionally length-normalised.
# mitoPPS is a COMPOSITION measure built from all-pairwise pathway ratios and is
# deliberately robust to total content, so the deviation is defensible - but it
# is a deviation and is recorded here and in the saved object. It also makes the
# standing CLAUDE.md rule bite harder, not less: CCLE mitoPPS values are NEVER
# comparable to TCGA mitoPPS values. Only the PATTERN transfers.
#
# GSVA IS COHORT-RELATIVE, and so is a z-score. Every line is scored in ONE run,
# and the pan-cancer fit RE-SCORES all four rulers rather than reusing the
# breast ones. CCLE scores are never pooled with TCGA or SCAN-B scores.
#
# N3: never "primed", of a transcript or of a cell line.
# SPECIES: human. Human MitoCarta, rebuilt from THIS repo's pinned workbook.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
suppressPackageStartupMessages(library(data.table))

message("\nE17: B4 - the OXPHOS main effect on guardian dependency\n",
        strrep("=", 78))

PATH_E17       <- file.path(DIR_RESULTS, "depmap_ox_dependency.rds")
PATH_E17_COEF  <- file.path(DIR_TABLES,  "E17_ox_coefficients.csv")
PATH_E17_ARMS  <- file.path(DIR_TABLES,  "E17_arm_panel.csv")

# =============================================================================
# 0. THE GUARD. Before any read.
# =============================================================================
DIR_DEPMAP <- here::here("data", "raw", "depmap")

DEPMAP_RELEASE <- "26Q1"
PRISM_RELEASE  <- "24Q2"
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
message("   3 required files readable through it")

PATH_MODEL  <- file.path(DIR_DEPMAP, "Model.csv")
PATH_EXPR   <- file.path(DIR_DEPMAP, EXPR_FILE)
PATH_CRISPR <- file.path(DIR_DEPMAP, "CRISPRGeneEffect.csv")
PATH_PRISM  <- file.path(DIR_DEPMAP,
  "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv")
PATH_PRISM_TREAT <- file.path(DIR_DEPMAP,
  "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")

# KEPT FROM SCRIPT 14. The 24Q4 CRISPR matrix was fetched from figshare before
# the portal release was checked; figshare never mirrored past 24Q4. It has been
# renamed with a _SUPERSEDED suffix, but a re-download under the original name
# would resurrect it SILENTLY, and nothing inside the file records its release.
# Byte length is the only signal that distinguishes it.
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
message("   byte guard passed (CRISPR ",
        format(file.size(PATH_CRISPR), big.mark = ","), " bytes, not the 24Q4 ",
        format(CRISPR_24Q4_BYTES, big.mark = ","), ")")

have_prism <- file.exists(PATH_PRISM) && file.exists(PATH_PRISM_TREAT)
message("   PRISM ", PRISM_RELEASE, ": ", if (have_prism) "present" else
        "ABSENT - the drug arm is SKIPPED, not faked")

# =============================================================================
# 0b. Constants
# =============================================================================
LINEAGE      <- "Breast"
DEP_PRIMARY  <- c("MCL1", "BCL2L1")
# UNCHANGED from script 14. POLR2A is requested even though it is known absent
# from 26Q1, so the absence is reported rather than quietly assumed.
DEP_CONTROL  <- c("BCL2", "BBC3", "BCL2L11", "BAX", "BAK1", "RPL3", "POLR2A")
# Declaration section 5.1. A FLOOR, not an endpoint; no direction is predicted.
DEP_FLOOR_SUB <- "POLR2B"
FLOOR_GENES   <- c("RPL3", "POLR2A", DEP_FLOOR_SUB)

DRUGS <- c("S63845", "AMG-176", "AZD5991",          # MCL1, all three present
           "A-1331852", "navitoclax", "A-1155463",  # BCL-XL; only navitoclax
           "venetoclax")                            # BCL2 specificity control
DRUG_DIFF <- c("navitoclax", "venetoclax")  # the ONLY interpretable BCL-XL read

ARM_PRIMARY  <- "OXPHOS subunits"
RULERS       <- c("ox_gsva", "ox_ppd", "ox_lvl", "ox_rel")
CI_LEVEL     <- 0.95
MIN_LINES    <- 25L
GSVA_MIN_SET <- 3L      # matches script 07/14; CII subunits is 4 genes and that
                        # is correct, not broken. Changing it would change the
                        # mitoPPS universe and so every other arm's value.
MIN_SET_FRAC <- 0.80    # below this is a harmonisation failure, not attrition

# The declared directions. Chronos: more negative = more essential.
DECLARED <- c(BCL2L1 = -1, MCL1 = +1)

# =============================================================================
# 1. Gene sets - REBUILT FROM THIS REPO'S PINNED MITOCARTA
# =============================================================================
# Script 14 read `arm_sets` from the validation repo's own script 07. Here they
# are rebuilt from `data/mitocarta_human/`, so the arms are the ones E16 used
# and B4 sits on this repo's inputs rather than on a neighbour's derived object.
#
# THE TWO BUILDS DO NOT DIFFER, and that was checked rather than assumed: the
# workbooks are byte-identical (MD5 3c0bd24e362238216e142bc708e41286 in both
# repos) and all 18 arms come out set-identical to the snapshot's `arm_sets`.
# The assertion below is that check, kept in the script.
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
# mtDNA-encoded transcripts are orders of magnitude more abundant, so a pathway
# containing even one has a severely inflated mean. They are removed from their
# canonical pathways and held in one synthetic pathway. Script 07 3a.
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

# The check that the rebuild is the same object E16 and the TCGA arms used.
snap_arms <- readRDS(PATH_TCGA_MITO)$arm_sets
same_arm  <- vapply(names(ARM_SETS),
                    function(a) setequal(ARM_SETS[[a]], snap_arms[[a]]), logical(1))
if (!all(same_arm)) {
  stop("arm(s) rebuilt from this repo's MitoCarta differ from the snapshot's ",
       "arm_sets: ", paste(names(same_arm)[!same_arm], collapse = ", "),
       ". The two builds have diverged and B4 would not be on E16's arms.",
       call. = FALSE)
}
message("   18 arms rebuilt, all ", sum(same_arm), " set-identical to the ",
        "snapshot's arm_sets")

sd_       <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
MITO_ALL  <- sd_$strip_refs$MITOCARTA_ALL
MYC_SET   <- sd_$myc_sets[[MYC_REF]]           # FELSHER__MITOSTRIP, 61
PROLIF    <- sd_$cov_sets$PROLIF_DISJOINT      # 318
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL,
          length(MYC_SET) == EXPECT_FELSHER_STRIP, length(PROLIF) == 318L)
message("   MYC = ", MYC_REF, " (", length(MYC_SET), "), PROLIF_DISJOINT (",
        length(PROLIF), "), MITOCARTA_ALL (", length(MITO_ALL), ")")

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
# 3. Expression, on two scales that never meet
# =============================================================================
message("\n3. expression")

.read_matrix <- function(path, what) {
  m <- data.table::fread(path, data.table = FALSE)
  rn <- as.character(m[[1]]); m <- as.matrix(m[, -1, drop = FALSE])
  rownames(m) <- rn
  colnames(m) <- trimws(sub("\\s*\\(\\d+\\)$", "", colnames(m)))
  dup <- duplicated(colnames(m))
  if (any(dup)) {
    message("   ", what, ": ", sum(dup), " duplicate symbols after stripping ",
            "Entrez ids - first occurrence kept")
    m <- m[, !dup, drop = FALSE]
  }
  m
}

# 26Q1's expression file is ONE ROW PER SEQUENCING RUN, not per model, and its
# own release notes misdescribe it - they name ProfileID and is_default_entry,
# neither of which exists, and the flags are "Yes"/"No" strings not logicals.
# So this reader trusts the FILE: it locates columns by name with fallbacks and
# decides what is a gene by testing the column is numeric.
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
           ". Refusing to guess - a wrong guess here either empties the matrix ",
           "or keeps every duplicate row.", call. = FALSE)
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
  # Deliberately NOT IsDefaultEntryForMC - the model-CONDITION default is a
  # different and larger set; in 26Q1 exactly one row differs between them.
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
    stop("non-numeric column(s) survived the metadata filter and would be ",
         "treated as genes: ",
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
# The scale check matters because the failure is silent: log2(TPM+1) is bounded
# well under 20. A linear matrix here would put GSVA on the wrong scale and
# mitoPPS on a doubly-exponentiated one.
if (max(EXPR, na.rm = TRUE) > 40) {
  stop("expression matrix max is ", round(max(EXPR, na.rm = TRUE), 1),
       ", which is not log2(TPM+1). Check the file.", call. = FALSE)
}

# =============================================================================
# 4. Scoring - the four rulers, in one run over one cohort
# =============================================================================
# `.score_cohort` builds every ruler for one set of lines. It is a function
# because the pan-cancer fit MUST re-score rather than reuse: GSVA is
# cohort-relative and so is a z-score, so a breast score and a pan-cancer score
# are different quantities and cannot be mixed.
message("\n4. scoring")

.z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)

# mitoPPS, same algorithm as E02/script 07: pathway means over a DECLARED
# universe, all pairwise ratios, each corrected by its across-sample mean.
# Changing the universe changes every value, which is why the 18 arms are fixed.
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
         "zero or non-finite, so the pairwise ratio is undefined. First few: ",
         paste(lab, collapse = "; "),
         ".\nDecide explicitly - drop the arm from the universe, or drop the ",
         "line - and record it. Do NOT add a floor silently: mitoPPS is a ",
         "composition measure and a floor moves every ratio in the universe.",
         call. = FALSE)
  }
}

# E16's composite, verbatim: mean of per-gene z across samples. Here the input is
# log2(TPM+1), which is the same shape as E16's log2(linear + 1).
.comp <- function(genes, E) {
  M <- E[intersect(genes, rownames(E)), , drop = FALSE]
  v <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}

.score_cohort <- function(ids, what) {
  message("   scoring ", what, ": ", length(ids), " lines")
  # LOG. GSVA and the two composites.
  E_LOG <- t(EXPR[ids, , drop = FALSE])
  E_LOG <- E_LOG[stats::complete.cases(E_LOG), , drop = FALSE]

  sets_defined <- c(list(MYC = MYC_SET, PROLIF = PROLIF), ARM_SETS)
  sets_gsva <- lapply(sets_defined, function(s) intersect(s, rownames(E_LOG)))
  cov_tab <- tibble::tibble(scope = what, set = names(sets_defined),
                            n_defined = lengths(sets_defined),
                            n_ccle = lengths(sets_gsva)) %>%
    dplyr::mutate(frac = round(n_ccle / n_defined, 3))
  # Two different failures, which the fraction separates and a count does not:
  # a set can be small because the pathway is small, or small because symbols
  # did not harmonise. CII subunits is the first kind - 4 of 4 is perfect.
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
         paste(tiny, collapse = ", "), ". Dropping an arm would change the ",
         "mitoPPS universe and so every other arm's value.", call. = FALSE)
  }

  GS <- GSVA::gsva(GSVA::gsvaParam(exprData = E_LOG, geneSets = sets_gsva,
                                   kcdf = "Gaussian", minSize = GSVA_MIN_SET,
                                   maxSize = Inf), verbose = FALSE)
  # GSVA drops sets below minSize SILENTLY, returning a smaller matrix; a
  # dropped arm would surface far downstream as a subscript error.
  miss <- setdiff(names(sets_gsva), rownames(GS))
  if (length(miss)) {
    stop("GSVA silently dropped ", length(miss), " set(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  # LINEAR. mitoPPS only. Built once from E_LOG and never mixed back.
  # NO FLOOR IS ADDED - .check_positive names the offending arm and line instead.
  E_LIN <- 2^E_LOG - 1
  E_LIN[E_LIN < 0] <- 0                 # floating point on an exact zero only
  S_LIN <- .pathway_means(E_LIN, lapply(ARM_SETS, function(s)
    intersect(s, rownames(E_LOG))))
  .check_positive(S_LIN, what)
  PPS <- .mitopps_universe(S_LIN)

  # ox_lvl and ox_rel, per ARM, on the LOG matrix. E16's recipe generalised the
  # way the mouse's relify() generalises it: an arm's relative score is its own
  # composite minus the composite of MitoCarta WITHOUT it. For the OXPHOS arm
  # that is exactly E16's ox_rel.
  LVL <- t(vapply(ARM_SETS, function(s) .comp(s, E_LOG), numeric(ncol(E_LOG))))
  REL <- t(vapply(ARM_SETS, function(s)
    .comp(s, E_LOG) - .comp(setdiff(MITO_ALL, s), E_LOG), numeric(ncol(E_LOG))))
  rownames(LVL) <- rownames(REL) <- names(ARM_SETS)

  # Every ruler z-scored WITHIN this cohort. Never across.
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
  stop("only ", length(lines_b), " ", LINEAGE, " lines have expression; below ",
       MIN_LINES, " this is not worth fitting.", call. = FALSE)
}
message("   ", length(lines_b), " ", LINEAGE, " lines with expression")

# =============================================================================
# 5. CRISPR
# =============================================================================
message("\n5. CRISPR gene effect")

CRISPR <- .read_matrix(PATH_CRISPR, "CRISPR")
message("   CRISPR: ", nrow(CRISPR), " models x ", ncol(CRISPR), " genes")

dep_lines <- intersect(lines_b, rownames(CRISPR))
message("   ", length(dep_lines), " ", LINEAGE,
        " lines have BOTH expression and CRISPR")
if (length(dep_lines) < MIN_LINES) {
  stop("only ", length(dep_lines), " lines carry both; below ", MIN_LINES,
       ". Report the gap; do not widen the lineage to manufacture n.",
       call. = FALSE)
}

GENES_WANTED <- c(DEP_PRIMARY, DEP_CONTROL, DEP_FLOOR_SUB)
missing_dep  <- setdiff(GENES_WANTED, colnames(CRISPR))
if (length(missing_dep)) {
  message("   NOT screened in ", DEPMAP_RELEASE, ": ",
          paste(missing_dep, collapse = ", "))
}
GENES_DEP <- intersect(GENES_WANTED, colnames(CRISPR))
# The primary endpoints and the substitute floor are not optional. If POLR2B
# were also absent the floor would be one gene and the declaration's B4-b could
# not be evaluated as written, so that stops rather than degrades.
if (!all(DEP_PRIMARY %in% GENES_DEP)) {
  stop("a PRIMARY endpoint is not screened: ",
       paste(setdiff(DEP_PRIMARY, GENES_DEP), collapse = ", "), call. = FALSE)
}
if (!DEP_FLOOR_SUB %in% GENES_DEP) {
  stop(DEP_FLOOR_SUB, " is not screened either, so the floor would rest on ",
       "RPL3 alone and B4-b cannot be evaluated as declared. Pick another ",
       "pan-essential and AMEND THE DECLARATION before re-running.",
       call. = FALSE)
}
message("   ", length(GENES_DEP), " endpoint/control genes screened; floor = ",
        paste(intersect(FLOOR_GENES, GENES_DEP), collapse = " + "))

# --- 5.1 the fitting machinery ----------------------------------------------
# BOTH models, every time. `OX` is extracted from each.
#
#   m_add  Y ~ MYC + OX + PROLIF          PRIMARY. V1's estimand.
#   m_int  Y ~ MYC * OX + PROLIF          SENSITIVITY. `OX` is the effect at
#                                         mean MYC, since MYC is z-scored.
#
# They should agree, because Block G established MYC:OX is null on these genes.
# A LARGE DISAGREEMENT IS A DIAGNOSTIC, NOT A RESULT: it would say the
# interaction is not negligible, in which case the additive `OX` is not
# interpretable as a main effect. It is reported as a diagnostic and NEVER as an
# interaction claim. Only `OX` is ever extracted from m_int.
.tidy <- function(m, term, label, ruler, model, extra = list()) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) {
    stop("term '", term, "' absent from fit '", label, "'", call. = FALSE)
  }
  crit <- stats::qt(1 - (1 - CI_LEVEL) / 2, m$df.residual)
  out <- tibble::tibble(
    label = label, term = term, ruler = ruler, model = model,
    n = stats::nobs(m),
    estimate = co[term, 1], se = co[term, 2],
    ci_lo = co[term, 1] - crit * co[term, 2],
    ci_hi = co[term, 1] + crit * co[term, 2],
    p = co[term, 4])
  if (length(extra)) out <- dplyr::bind_cols(out, tibble::as_tibble(extra))
  out
}

RES <- list(); .add <- function(x) if (!is.null(x)) RES[[length(RES) + 1L]] <<- x
# A SKIPPED FIT IS RECORDED, NEVER SILENT. This is the same principle as the
# symlink guard: an endpoint that quietly fails the n floor disappears from the
# results table looking exactly like an endpoint that was never asked for. The
# BCL-XL drug leg is the live case - navitoclax and venetoclax are in a
# different PRISM screen from the MCL1 compounds and cover fewer of these lines.
SKIPS <- list()
.skip <- function(label, ruler, extra, n_have, why)
  SKIPS[[length(SKIPS) + 1L]] <<- dplyr::bind_cols(
    tibble::tibble(label = label, ruler = ruler, n_complete = n_have,
                   reason = why), tibble::as_tibble(extra))

.fit_ox <- function(d, label, ruler, extra, lineage_term = FALSE) {
  n_have <- sum(stats::complete.cases(d))
  if (n_have < MIN_LINES) {
    .skip(label, ruler, extra, n_have,
          paste0("only ", n_have, " complete cases, below the MIN_LINES floor ",
                 "of ", MIN_LINES))
    return(invisible(NULL))
  }
  f_add <- if (lineage_term) Y ~ MYC + OX + PROLIF + lineage else
    Y ~ MYC + OX + PROLIF
  f_int <- if (lineage_term) Y ~ MYC * OX + PROLIF + lineage else
    Y ~ MYC * OX + PROLIF
  .add(.tidy(stats::lm(f_add, data = d), "OX", label, ruler, "additive", extra))
  .add(.tidy(stats::lm(f_int, data = d), "OX", label, ruler, "int-sensitivity",
             extra))
  invisible(NULL)
}

SB <- .score_cohort(lines_b, LINEAGE)

.frame_b <- function(rows, ruler, arm) tibble::tibble(
  MYC = SB$MYC[rows], OX = SB$SC[[ruler]][arm, rows], PROLIF = SB$PROLIF[rows])

# --- 5a. BREAST. PRIMARY. ----------------------------------------------------
# Deliberately minimal at n ~ 51: adding subtype here would spend the degrees of
# freedom the estimate needs.
message("\n5a. breast, primary (n = ", length(dep_lines), ")")
for (ruler in RULERS) {
  for (g in GENES_DEP) {
    d <- .frame_b(dep_lines, ruler, ARM_PRIMARY)
    d$Y <- CRISPR[dep_lines, g]
    .fit_ox(d, "B4 breast", ruler,
            list(gene = g, arm = ARM_PRIMARY, scope = LINEAGE))
  }
}

# --- 5b. the 18-arm panel, primary genes. SPECIFICITY. ----------------------
# THERE IS NO EXPRESSION-MATCHED NULL IN CCLE. The TCGA nulls are not
# transferable, so this is a RANK ORDERING, not a calibrated comparison, and no
# p-value is attached to it downstream. If OXPHOS is positive here, that null
# has to be BUILT in CCLE before the result is reportable. Named follow-up, not
# an optional extra. Declaration 6.1.
message("5b. the 18-arm panel, primary genes")
for (ruler in RULERS) {
  for (arm in names(ARM_SETS)) {
    for (g in DEP_PRIMARY) {
      d <- .frame_b(dep_lines, ruler, arm)
      d$Y <- CRISPR[dep_lines, g]
      .fit_ox(d, "B4 arm panel", ruler,
              list(gene = g, arm = arm, scope = LINEAGE))
    }
  }
}

# --- 5c. pan-cancer, lineage-adjusted. SECONDARY, better powered. ------------
# ~51 breast lines cannot support this on their own; the pan-cancer fit has 20x
# the n. It is SECONDARY because the mouse is mammary and both human cohorts are
# breast, and because the MYC/OXPHOS relationship need not be lineage-invariant.
# It is reported because a breast null at n = 51 is uninformative and saying so
# without offering the powered comparison would be a half-report.
#
# It RE-SCORES. GSVA is cohort-relative and so is a z-score, so a pan-cancer
# score is a different quantity from a breast-only one. This is the slow step.
pan_lines <- intersect(rownames(EXPR), rownames(CRISPR))
message("5c. pan-cancer (n = ", length(pan_lines), "), re-scoring all four ",
        "rulers - minutes, not seconds")
SP <- .score_cohort(pan_lines, "pan-cancer")
lin_pan  <- MODEL$lineage[match(pan_lines, MODEL$ModelID)]
keep_lin <- lin_pan %in% names(which(table(lin_pan) >= 10L))
message("   ", sum(keep_lin), " lines in lineages of >= 10")

for (ruler in RULERS) {
  for (g in GENES_DEP) {
    d <- tibble::tibble(MYC = SP$MYC[pan_lines],
                        OX = SP$SC[[ruler]][ARM_PRIMARY, pan_lines],
                        PROLIF = SP$PROLIF[pan_lines],
                        lineage = factor(lin_pan),
                        Y = CRISPR[pan_lines, g])[keep_lin, ]
    .fit_ox(droplevels(d), "B4 pan-cancer (secondary)", ruler,
            list(gene = g, arm = ARM_PRIMARY, scope = "pan-cancer"),
            lineage_term = TRUE)
  }
}

# =============================================================================
# 6. PRISM drug sensitivity
# =============================================================================
# PRISM Repurposing log2 fold change: MORE NEGATIVE = MORE SENSITIVE, the same
# direction as Chronos. Declared prediction, MCL1 leg: `OX` POSITIVE - OXPHOS-
# high lines LESS sensitive to MCL1 inhibition. That is the counter-intuitive
# leg; the field expects MYC-driven breast cancer to be MCL1-dependent.
#
# BCL-XL IS NOT COVERED. A-1331852 and A-1155463 are both absent from
# Repurposing. Navitoclax hits BCL2, BCL-XL and BCL-W, and BCL2 transcript ALSO
# falls with OXPHOS (-0.117), so navitoclax alone reads a mixture of two things
# moving the same way.
#
#   THE INTERPRETABLE QUANTITY IS NAVITOCLAX MINUS VENETOCLAX ON THE SAME LINES.
#   It is fitted below as its own endpoint. A NAVITOCLAX RESULT IS NEVER
#   REPORTED AS A BCL-XL RESULT.
message("\n6. drug sensitivity")

drug_res <- NULL; prism_note <- "not attempted"
if (!have_prism) {
  prism_note <- "SKIPPED - PRISM files absent through the symlink"
  message("   ", prism_note, ". Deferred, not dropped.")
} else {
  TREAT <- data.table::fread(PATH_PRISM_TREAT, data.table = FALSE)
  nm_col <- intersect(c("Drug.Name", "name", "Name"), names(TREAT))[1]
  id_col <- intersect(c("IDs", "column_name", "broad_id"), names(TREAT))[1]
  if (is.na(nm_col) || is.na(id_col)) {
    prism_note <- "SKIPPED - PRISM compound list has unexpected columns"
    message("   ", prism_note, " (", paste(utils::head(names(TREAT)),
            collapse = ", "), ") - skipped rather than guessed")
  } else {
    # ORIENTATION: the Repurposing matrix is COMPOUNDS x CELL LINES, the
    # transpose of the expression and CRISPR matrices. Read without transposing,
    # the line intersection is EMPTY and this section reports "0 lines" rather
    # than failing. Hence the assertion, not just the t().
    PR <- data.table::fread(PATH_PRISM, data.table = FALSE)
    rn <- PR[[1]]; PR <- as.matrix(PR[, -1, drop = FALSE]); rownames(PR) <- rn
    if (!all(grepl("^BRD", utils::head(rownames(PR), 20))) ||
        !all(grepl("^ACH-", utils::head(colnames(PR), 20)))) {
      stop("PRISM matrix orientation is not compounds x ModelIDs (rows start '",
           substr(rownames(PR)[1], 1, 12), "', columns start '",
           substr(colnames(PR)[1], 1, 12), "').", call. = FALSE)
    }
    PR <- t(PR)                                       # now ModelIDs x compounds
    # Drug.Name is UPPERCASE for some compounds and mixed for others.
    hit <- TREAT[tolower(TREAT[[nm_col]]) %in% tolower(DRUGS), , drop = FALSE]
    absent <- setdiff(tolower(DRUGS), tolower(hit[[nm_col]]))
    message("   ", nrow(hit), " of ", length(DRUGS), " named compounds present")
    if (length(absent)) {
      message("   NOT in ", PRISM_RELEASE, ": ", paste(absent, collapse = ", "),
              "\n   -> both SELECTIVE BCL-XL inhibitors. Read navitoclax only ",
              "as the difference below.")
    }
    pr_lines <- intersect(dep_lines, rownames(PR))
    message("   ", length(pr_lines), " ", LINEAGE, " lines with PRISM")
    prism_note <- paste0(nrow(hit), " of ", length(DRUGS), " compounds, ",
                         length(pr_lines), " ", LINEAGE, " lines")

    .drug_col <- function(nm) {
      k <- which(tolower(hit[[nm_col]]) == tolower(nm))
      if (!length(k)) return(NA_character_)
      cc <- as.character(hit[[id_col]][k])
      cc <- cc[cc %in% colnames(PR)]
      if (!length(cc)) NA_character_ else cc[1]
    }

    # Per-compound coverage, printed BEFORE any fit. The 24Q2 release is two
    # screens - REP.1M and REP.PRIMARY - and they do not cover the same lines,
    # so a compound can be present in the release and still be unfittable here.
    drug_cover <- dplyr::bind_rows(lapply(seq_len(nrow(hit)), function(k) {
      col <- as.character(hit[[id_col]][k])
      tibble::tibble(drug = as.character(hit[[nm_col]][k]),
                     screen = if ("screen" %in% names(hit)) as.character(hit$screen[k]) else NA_character_,
                     in_matrix = col %in% colnames(PR),
                     n_lines = if (col %in% colnames(PR))
                       sum(!is.na(PR[pr_lines, col])) else 0L)
    })) %>% dplyr::mutate(fittable = n_lines >= MIN_LINES)
    message("   per-compound coverage on those lines (floor ", MIN_LINES, "):")
    drug_cover %>% as.data.frame() %>% print(row.names = FALSE)
    if (any(!drug_cover$fittable)) {
      message("   NOT FITTABLE: ",
              paste(drug_cover$drug[!drug_cover$fittable], collapse = ", "),
              " - present in the release, too few lines here.")
    }

    if (nrow(hit) && length(pr_lines) >= MIN_LINES) {
      for (ruler in RULERS) {
        # each named compound
        for (k in seq_len(nrow(hit))) {
          col <- as.character(hit[[id_col]][k])
          if (!col %in% colnames(PR)) next
          d <- .frame_b(pr_lines, ruler, ARM_PRIMARY)
          d$Y <- PR[pr_lines, col]
          .fit_ox(d, "B4 drug", ruler,
                  list(gene = as.character(hit[[nm_col]][k]),
                       arm = ARM_PRIMARY, scope = LINEAGE))
        }
        # THE DERIVED ENDPOINT. Not a new gene - the only reading of the
        # BCL-XL/BCL-W component available in this release.
        c_nav <- .drug_col(DRUG_DIFF[1]); c_ven <- .drug_col(DRUG_DIFF[2])
        ex_diff <- list(gene = "navitoclax_minus_venetoclax",
                        arm = ARM_PRIMARY, scope = LINEAGE)
        if (is.na(c_nav) || is.na(c_ven)) {
          .skip("B4 drug", ruler, ex_diff, 0L,
                "navitoclax or venetoclax has no column in this PRISM release")
        } else {
          d <- .frame_b(pr_lines, ruler, ARM_PRIMARY)
          d$Y <- PR[pr_lines, c_nav] - PR[pr_lines, c_ven]
          .fit_ox(d, "B4 drug", ruler, ex_diff)
        }
      }
    } else {
      message("   too few lines or compounds to fit - reported, not faked")
    }
  }
}

coefs <- dplyr::bind_rows(RES)
skipped <- if (length(SKIPS)) dplyr::bind_rows(SKIPS) else
  tibble::tibble(label = character(), ruler = character(),
                 n_complete = integer(), reason = character())
message("\n   ", nrow(coefs), " OX coefficients extracted (",
        dplyr::n_distinct(coefs$label), " labels x ", length(RULERS),
        " rulers x 2 models)")
if (nrow(skipped)) {
  message("   ", nrow(skipped), " fit(s) SKIPPED and recorded - never silent:")
  skipped %>%
    dplyr::group_by(label, gene, reason) %>%
    dplyr::summarise(rulers = dplyr::n(), .groups = "drop") %>%
    as.data.frame() %>% print(row.names = FALSE)
} else {
  message("   no fits skipped")
}

# =============================================================================
# 7. B4-a to B4-e, on the rules fixed in the declaration
# =============================================================================
# THE CRITERION IS THE SIGN OF THE ESTIMATE, read across rulers and scopes.
# 95% CIs are reported for every coefficient and it is stated whether each
# excludes zero, but NO P-VALUE THRESHOLD GATES ANY VERDICT. Declaration 8.
message("\n7. the declared verdicts")

.est <- function(lab, g, mdl = "additive", arm = ARM_PRIMARY) {
  coefs %>% dplyr::filter(label == lab, gene == g, model == mdl, arm == !!arm)
}

# --- B4-a --------------------------------------------------------------------
b4a_tab <- .est("B4 breast", "BCL2L1") %>%
  dplyr::bind_rows(.est("B4 breast", "MCL1")) %>%
  dplyr::mutate(predicted = DECLARED[gene],
                as_declared = sign(estimate) == predicted,
                ci_excludes_0 = ci_lo > 0 | ci_hi < 0)
message("\n   B4-a: BCL2L1 OX negative and MCL1 OX positive, breast")
b4a_tab %>%
  dplyr::transmute(gene, ruler, n, est = sprintf("%+.3f", estimate),
                   ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
                   as_declared, ci_excludes_0) %>%
  as.data.frame() %>% print(row.names = FALSE)
n_rulers_ok <- b4a_tab %>% dplyr::group_by(gene) %>%
  dplyr::summarise(n_ok = sum(as_declared), .groups = "drop")
both_ok <- all(n_rulers_ok$n_ok >= 2L)
one_ok  <- any(n_rulers_ok$n_ok >= 2L)
b4a <- if (both_ok) "PASS" else if (one_ok) "PARTIAL" else "FAIL"
message("   rulers in the declared direction: ",
        paste(sprintf("%s %d/%d", n_rulers_ok$gene, n_rulers_ok$n_ok,
                      length(RULERS)), collapse = ", "),
        "  ->  B4-a ", b4a)

# --- B4-b. HARD REQUIREMENT. -------------------------------------------------
# Not a qualifier: an interaction is robust to a uniform shift in dependency and
# a main effect is not, so a floor signal means the result is growth rate.
floor_here <- intersect(FLOOR_GENES, GENES_DEP)
b4b_tab <- coefs %>%
  dplyr::filter(label == "B4 breast", model == "additive",
                gene %in% floor_here) %>%
  dplyr::mutate(ci_excludes_0 = ci_lo > 0 | ci_hi < 0)
prim_max <- max(abs(b4a_tab$estimate))
b4b_tab  <- b4b_tab %>%
  dplyr::mutate(vs_primary_max = round(abs(estimate) / prim_max, 2))
message("\n   B4-b: the floor (", paste(floor_here, collapse = " + "), ")")
b4b_tab %>%
  dplyr::transmute(gene, ruler, est = sprintf("%+.3f", estimate),
                   ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
                   ci_excludes_0, vs_primary_max) %>%
  as.data.frame() %>% print(row.names = FALSE)
b4b_rpl3  <- !any(b4b_tab$ci_excludes_0[b4b_tab$gene == "RPL3"])
b4b_both  <- !any(b4b_tab$ci_excludes_0)
b4b <- if (b4b_rpl3 && b4b_both) "PASS" else "FAIL"
message("   RPL3 alone: ", if (b4b_rpl3) "clean" else "SIGNAL",
        " | RPL3 + ", DEP_FLOOR_SUB, ": ", if (b4b_both) "clean" else "SIGNAL",
        "  ->  B4-b ", b4b)
if (b4b == "FAIL") {
  message("   B4-b VOIDS B4-a. The breast result is a growth-rate artefact and",
          "\n   is reported as one, not as a dependency result.")
}

# --- B4-c. A RANK, pending the null. ----------------------------------------
b4c_tab <- coefs %>%
  dplyr::filter(label == "B4 arm panel", model == "additive") %>%
  dplyr::group_by(ruler, gene) %>%
  dplyr::mutate(rank_by_signed = rank(estimate * DECLARED[gene] * -1)) %>%
  dplyr::ungroup()
b4c_rank <- b4c_tab %>% dplyr::filter(arm == ARM_PRIMARY) %>%
  dplyr::transmute(gene, ruler, oxphos_est = round(estimate, 3),
                   rank_of_18 = as.integer(rank_by_signed))
message("\n   B4-c: where OXPHOS subunits ranks among the 18 arms")
message("   (rank 1 = most extreme IN THE DECLARED DIRECTION for that gene)")
b4c_rank %>% as.data.frame() %>% print(row.names = FALSE)
message("   NO P-VALUE IS ATTACHED. There is no expression-matched null in ",
        "CCLE;\n   the TCGA nulls are not transferable. This is a rank ordering.")

# --- B4-d. Pan-cancer sign agreement. ---------------------------------------
b4d_tab <- coefs %>%
  dplyr::filter(label %in% c("B4 breast", "B4 pan-cancer (secondary)"),
                model == "additive", gene %in% DEP_PRIMARY,
                arm == ARM_PRIMARY) %>%
  dplyr::select(gene, ruler, scope, estimate) %>%
  tidyr::pivot_wider(names_from = scope, values_from = estimate) %>%
  dplyr::mutate(agree = sign(Breast) == sign(`pan-cancer`))
message("\n   B4-d: does pan-cancer agree in sign with breast?")
b4d_tab %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
  as.data.frame() %>% print(row.names = FALSE)
b4d <- if (all(b4d_tab$agree)) "AGREES" else
  paste0("DISAGREES on ", sum(!b4d_tab$agree), " of ", nrow(b4d_tab),
         " - reported as lineage-dependence, NOT resolved")
message("   B4-d: ", b4d)

# --- B4-e. PRISM MCL1 leg. ---------------------------------------------------
mcl1_drugs <- c("S63845", "AMG-176", "AZD5991")
b4e_tab <- coefs %>%
  dplyr::filter(label == "B4 drug", model == "additive") %>%
  dplyr::mutate(leg = dplyr::case_when(
    toupper(gene) %in% toupper(mcl1_drugs) ~ "MCL1",
    gene == "navitoclax_minus_venetoclax"  ~ "BCL-XL/W (difference)",
    TRUE                                   ~ "other"))
if (nrow(b4e_tab)) {
  message("\n   B4-e: PRISM. More negative = more sensitive, as Chronos.")
  b4e_tab %>%
    dplyr::transmute(gene, leg, ruler, n, est = sprintf("%+.3f", estimate),
                     ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi)) %>%
    dplyr::arrange(leg, gene, ruler) %>%
    as.data.frame() %>% print(row.names = FALSE)
  mcl1_crispr <- .est("B4 breast", "MCL1")$estimate
  mcl1_drug   <- b4e_tab$estimate[b4e_tab$leg == "MCL1"]
  b4e <- if (!length(mcl1_drug)) "NOT EVALUABLE" else
    sprintf("%d of %d MCL1 drug coefficients agree in sign with the CRISPR MCL1 result",
            sum(sign(mcl1_drug) == sign(stats::median(mcl1_crispr))),
            length(mcl1_drug))
  message("   B4-e: ", b4e)
  # The BCL-XL leg is separate, and its absence is a REPORTED absence.
  n_diff <- sum(b4e_tab$leg == "BCL-XL/W (difference)")
  if (!n_diff) {
    message("   THE BCL-XL/W LEG WAS NOT EVALUABLE. navitoclax minus ",
            "venetoclax could\n   not be fitted - see the skip table above. ",
            "The BCL-XL half of B4-e\n   is therefore ABSENT, not null.")
  }
  message("   navitoclax is NEVER read as a BCL-XL result. Only the difference.")
} else {
  b4e <- "NOT EVALUABLE - no PRISM fits"
  message("\n   B4-e: ", b4e)
}

# --- the additive/interaction agreement DIAGNOSTIC --------------------------
# Not a result, and not an interaction claim. If the two `OX` estimates diverge,
# the additive one is not interpretable as a main effect and that must be said.
model_agree <- coefs %>%
  dplyr::select(label, ruler, gene, arm, scope, model, estimate) %>%
  tidyr::pivot_wider(names_from = model, values_from = estimate) %>%
  dplyr::mutate(delta = `int-sensitivity` - additive)
message("\n   DIAGNOSTIC, not a result: |OX(additive) - OX(at mean MYC)|")
model_agree %>% dplyr::group_by(label) %>%
  dplyr::summarise(n = dplyr::n(), median_abs = round(stats::median(abs(delta)), 4),
                   max_abs = round(max(abs(delta)), 4), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   A large value would say MYC:OX is not negligible after all, and",
        "\n   that the additive OX is not readable as a main effect. It would",
        "\n   NOT be an interaction claim - Block G settled that and B4 does",
        "\n   not revisit it.")

verdicts <- tibble::tibble(
  rule = c("B4-a", "B4-b", "B4-c", "B4-d", "B4-e"),
  question = c(
    "BCL2L1 OX negative and MCL1 OX positive, breast, on >= 2 rulers",
    "neither direction on the floor (RPL3, POLR2B). HARD REQUIREMENT",
    "does OXPHOS rank above the other mitochondrial arms",
    "does pan-cancer agree in sign with breast",
    "does the PRISM MCL1 leg agree in sign with CRISPR MCL1"),
  verdict = c(b4a, b4b, "rank only, pending the CCLE null", b4d, b4e))
message("\n   THE FIVE VERDICTS")
verdicts %>% as.data.frame() %>% print(row.names = FALSE)
if (b4a %in% c("PASS", "PARTIAL") && b4b == "PASS") {
  message("\n   *** THE CCLE MATCHED NULL IS NOW REQUIRED. *** B4-a is not ",
          "FAIL and\n   the floor is clean, so B4-c's rank cannot be left as a",
          " rank: the\n   expression-matched null has to be BUILT in CCLE ",
          "before any OXPHOS\n   arm result is reportable. That is the next ",
          "script.")
} else {
  message("\n   The CCLE matched null is NOT required by this result: B4-a is ",
          b4a, "\n   and B4-b is ", b4b, ", so there is no positive OXPHOS ",
          "arm result to calibrate.")
}

# =============================================================================
# 8. Save
# =============================================================================
saveRDS(list(
  coefficients = coefs, skipped = skipped, verdicts = verdicts,
  b4a = b4a_tab, b4b = b4b_tab, b4c = b4c_tab, b4c_rank = b4c_rank,
  b4d = b4d_tab, b4e = b4e_tab, model_agree = model_agree,
  coverage = dplyr::bind_rows(SB$coverage, SP$coverage),
  arm_sets_n = lengths(ARM_SETS),
  lines = list(breast_models = length(breast),
               breast_expression = length(lines_b),
               breast_both = length(dep_lines),
               pan_both = length(pan_lines),
               pan_kept = sum(keep_lin)),
  genes = list(primary = DEP_PRIMARY, control = DEP_CONTROL,
               floor_substitute = DEP_FLOOR_SUB,
               not_screened = missing_dep),
  spec = list(
    estimand = paste("the OX MAIN EFFECT. Block G fitted MYC*OX everywhere and",
                     "extracted only MYC:OX; OX was reported nowhere."),
    no_interaction_claim = paste("m_int is a SENSITIVITY and only its OX term",
                                 "is read. No interaction claim in either",
                                 "direction; Block G settled that."),
    declared = DECLARED,
    declaration = "docs/2026-09-04_b4_declaration.md, committed before the fit",
    two_sided = paste("a cell can need a protein more because it has less of",
                      "it. The reverse is reported and is NOT a pass."),
    no_null = paste("NO expression-matched null exists in CCLE; the TCGA nulls",
                    "are not transferable. The 18-arm panel is a RANK ordering",
                    "and carries no p-value. A positive OXPHOS result is not",
                    "reportable until that null is built."),
    mitopps_deviation = paste("mitoPPS runs on linear TPM, not linear",
                              "DESeq2-normalised counts, because that is what",
                              "DepMap ships. A composition measure over",
                              "all-pairwise ratios is robust to total content,",
                              "so this is defensible - but CCLE mitoPPS values",
                              "are NEVER comparable to TCGA values. Only the",
                              "pattern transfers."),
    composite_scale = "mean per-gene z of log2(TPM+1); E16's recipe, CCLE input",
    cohort_relative = paste("GSVA and every z-score are cohort-relative, so the",
                            "pan-cancer fit RE-SCORES all four rulers rather",
                            "than reusing the breast ones."),
    polr2a = paste("POLR2A is not screened in DepMap", DEPMAP_RELEASE,
                   "- every other RNA Pol II subunit is. DEP_CONTROL is",
                   "unchanged and POLR2A is reported absent; POLR2B is a",
                   "SUBSTITUTE FLOOR declared before the fit."),
    prism = prism_note,
    prism_bclxl = paste("A-1331852 and A-1155463, the two SELECTIVE BCL-XL",
                        "inhibitors, are absent from Repurposing", PRISM_RELEASE,
                        ". Navitoclax hits BCL2/BCL-XL/BCL-W and BCL2",
                        "transcript also falls with OXPHOS, so navitoclax alone",
                        "reads a mixture. The interpretable quantity is",
                        "navitoclax MINUS venetoclax on the same lines. A",
                        "navitoclax result is NEVER a BCL-XL result."),
    n3 = "these are gene-effect scores; the word 'primed' is used of nothing",
    releases = c(depmap = DEPMAP_RELEASE, prism = PRISM_RELEASE,
                 expression = EXPR_FILE),
    depmap_source = paste("symlink data/raw/depmap ->", Sys.readlink(DIR_DEPMAP),
                          "- see data/raw/depmap_README.md. Nothing is written",
                          "to myc_human_validation."),
    myc_estimator = MYC_REF, covariate = "PROLIF_DISJOINT",
    rulers = RULERS, seed = PROJECT_SEED),
  built = Sys.time()), PATH_E17)

readr::write_csv(coefs, PATH_E17_COEF)
readr::write_csv(b4c_tab, PATH_E17_ARMS)

message("\nE17: done.")
message("    results/depmap_ox_dependency.rds")
message("    outputs/tables/E17_ox_coefficients.csv")
message("    outputs/tables/E17_arm_panel.csv")
message("    no figures - B4 is a coefficient, and every one of them is in a ",
        "table.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E17)

  # The five verdicts, and the line counts behind them.
  x$verdicts %>% as.data.frame()
  utils::str(x$lines)
  x$genes$not_screened          # POLR2A, expected

  # Every fit that was NOT made, and why. Read this before believing an
  # endpoint is absent because it was null.
  x$skipped %>% as.data.frame()

  # B4-a. The two declared directions, four rulers each.
  x$b4a %>%
    dplyr::transmute(gene, ruler, n, estimate = round(estimate, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                     as_declared, ci_excludes_0) %>%
    as.data.frame()

  # B4-b. The floor. A signal here voids B4-a.
  x$b4b %>%
    dplyr::transmute(gene, ruler, estimate = round(estimate, 3),
                     ci_excludes_0, vs_primary_max) %>%
    as.data.frame()

  # B4-c. The rank, and the whole 18-arm panel behind it.
  x$b4c_rank %>% as.data.frame()
  x$b4c %>%
    dplyr::filter(ruler == "ox_gsva", gene == "BCL2L1") %>%
    dplyr::arrange(estimate) %>%
    dplyr::transmute(arm, estimate = round(estimate, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3)) %>%
    as.data.frame()

  # B4-d and B4-e.
  x$b4d %>% as.data.frame()
  x$b4e %>%
    dplyr::transmute(gene, leg, ruler, n, estimate = round(estimate, 3)) %>%
    as.data.frame()

  # The diagnostic: do the additive and the at-mean-MYC estimates agree?
  # NOT an interaction claim - see $spec$no_interaction_claim.
  summary(abs(x$model_agree$delta))
  x$model_agree %>% dplyr::slice_max(abs(delta), n = 5) %>% as.data.frame()

  # Set coverage in CCLE, both scopes. Worst is CV subunits.
  x$coverage %>% dplyr::arrange(frac) %>% utils::head(8) %>% as.data.frame()

  # The standing limitations, in the object rather than only in a note.
  x$spec$no_null
  x$spec$prism_bclxl
}
