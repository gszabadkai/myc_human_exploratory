# E31_bh3_mimetic_oxphos.R
# =============================================================================
# BH3-MIMETIC SENSITIVITY ON THE OXPHOS AXIS IN CANCER CELL LINES
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered.
#
# The estimand hierarchy, the directions, the population, the confounders and
# the stopping rules were fixed in docs/2026-09-28_bh3_mimetic_declaration.md
# and COMMITTED BEFORE THIS FILE WAS WRITTEN. The ordering is auditable in the
# git history of the bh3-mimetic-oxphos branch. Read that note first; this
# header restates the rules it must obey, it does not invent them.
#
# =============================================================================
# THE QUESTION
# =============================================================================
# The manuscript argues that a high respiratory state is lethal on a MYC
# background, and that tumours retaining it escaped either by lowering the
# oxidative driver or by blocking the apoptotic effector. Escape-by-blocking
# implies a DEPENDENCY: a tumour holding a high respiratory state needs the
# guardian that permits it, and should be selectively vulnerable to BCL-XL
# antagonism.
#
# This script asks only whether that shows up in cell line pharmacogenomics.
#
# =============================================================================
# WHY A NULL HERE IS UNINFORMATIVE AND A POSITIVE IS INTERESTING
# =============================================================================
# E18 (docs/2026-09-04_b4_addendum.md) already measured the thing this analysis
# depends on. The tumour BCL2L1-OXPHOS configuration IS ABSENT FROM CELL LINES:
# BCL2L1 runs +0.35 to +0.43 with OXPHOS in both tumour cohorts and -0.11 to
# +0.01 in CCLE, and the BCL2L1 - BBC3 gap is positive in 8 of 8 tumour cells
# and negative in 8 of 8 CCLE cells. Verdict 0 of 4 rulers, breast and
# pan-cancer, on a rule fixed before the numbers. The purity escape was closed
# independently by E16 check 4 BEFORE E18 was written.
#
#   - A NULL HERE IS UNINFORMATIVE. A system that does not carry the exposure
#     cannot be asked to carry the dependency.
#   - A POSITIVE HERE IS INTERESTING, because it would not be explained by the
#     configuration and would need another explanation.
#
# This asymmetry is declared, not discovered. It is NOT a licence to accept a
# weak positive - section 0b's guardrails are what a positive must clear.
#
# =============================================================================
# SCALE. Two objects that never meet. CLAUDE.md, and E17/E18 section 4.
# =============================================================================
#   - GSVA takes the DepMap matrix AS SUPPLIED: it is log2(TPM + 1), which is
#     already log scale. kcdf = "Gaussian". NOTHING IS RE-LOGGED.
#   - mitoPPS takes 2^x - 1, i.e. LINEAR TPM.
#
# GSVA IS COHORT-RELATIVE, and so is every z-score here. This panel is a
# COHORT. Its scores are NEVER compared numerically to a TCGA, SCAN-B, mouse or
# CCLE-breast value - only signs and orderings cross cohorts. All sets are
# scored in ONE call on the FINAL analysis cohort, so the gene universe and the
# half-matrix pins are shared.
#
# mitoPPS here runs on linear TPM, not linear DESeq2 counts - the deviation
# declared in data/raw/depmap_README.md and carried forward unchanged. It makes
# the no-cross-cohort rule bite harder, not less.
#
# N3: these are drug-response summary statistics and transcripts. No cell line
# is described as "primed".
# SPECIES: human. Human MitoCarta, rebuilt from THIS repo's pinned workbook.
# NO ORTHOLOG FUNCTION IS CALLED. This script reads no mouse object.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
suppressPackageStartupMessages(library(data.table))

message("\nE31: BH3-mimetic sensitivity on the OXPHOS axis\n", strrep("=", 78))

PATH_E31     <- file.path(DIR_RESULTS, "bh3_oxphos_gdsc.rds")
PATH_E31_CO  <- file.path(DIR_TABLES,  "E31_bh3_coefficients.csv")
PATH_E31_COV <- file.path(DIR_TABLES,  "E31_bh3_coverage.csv")
PATH_E31_FIG <- file.path(DIR_FIGURES, "E31_bh3_forest.png")
# The figure is written TWICE, and the second path is deliberate. outputs/ is
# gitignored, so a figure written only there does not survive a fresh clone -
# and this analysis runs on a branch whose raw data and results object are also
# gitignored, which leaves the note and the figure as the only durable record.
# Writing the tracked copy from the SCRIPT rather than copying it by hand means
# the author's run regenerates it and git shows any difference.
PATH_E31_FIG_DOC <- here::here("docs", "figures", "2026-09-28_E31_bh3_forest.png")

# =============================================================================
# 0. THE GUARDS. Before any read.
# =============================================================================
# 0a. DepMap, reached by symlink. E17/E18 section 0 verbatim - see
#     data/raw/depmap_README.md for why a broken symlink must STOP rather than
#     degrade to a skip: an absent directory and a dead link look identical,
#     and the run would complete having silently answered a smaller question.
DIR_DEPMAP     <- here::here("data", "raw", "depmap")
DEPMAP_RELEASE <- "26Q1"
EXPR_FILE      <- "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
PRISM_RELEASE  <- "24Q2"

.link_help <- paste0(
  "\n\ndata/raw/depmap is a SYMLINK, not a copy:\n",
  "  data/raw/depmap -> /Users/gs/code/myc_human_validation/data/raw/depmap\n",
  "Recreate it with:\n",
  "  ln -s /Users/gs/code/myc_human_validation/data/raw/depmap data/raw/depmap\n",
  "Provenance, checksums and the browser-only acquisition route: ",
  "data/raw/depmap_README.md")

.p <- DIR_DEPMAP
if (!nzchar(Sys.readlink(.p)) && !dir.exists(.p)) {
  stop("data/raw/depmap: symlink missing or unresolved.", .link_help,
       call. = FALSE)
}
if (nzchar(Sys.readlink(.p)) && !dir.exists(Sys.readlink(.p))) {
  stop("data/raw/depmap: symlink target does not resolve -> ", Sys.readlink(.p),
       .link_help, call. = FALSE)
}
for (f in c("Model.csv", EXPR_FILE)) {
  if (!file.exists(file.path(.p, f))) {
    stop("data/raw/depmap/", f, " unreadable through the symlink.", .link_help,
         call. = FALSE)
  }
}
message("\n0a. symlink resolves -> ", Sys.readlink(.p))

PATH_MODEL <- file.path(DIR_DEPMAP, "Model.csv")
PATH_EXPR  <- file.path(DIR_DEPMAP, EXPR_FILE)
PATH_PRISM_M <- file.path(DIR_DEPMAP,
  "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv")
PATH_PRISM_C <- file.path(DIR_DEPMAP,
  "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")

for (p in c(PATH_MODEL, PATH_EXPR)) {
  # An HTML Cloudflare challenge page saved under a .csv name is ~5 KB.
  if (file.size(p) < 1e5) {
    stop(basename(p), " is only ", file.size(p), " bytes - almost certainly ",
         "the Cloudflare verification HTML page saved under a .csv name.",
         .link_help, call. = FALSE)
  }
}
message("    byte guard passed")

# 0b. GDSC. Byte-exact, because these files are STATIC: release 8.5 was fitted
#     27 Oct 2023 and there is no newer fitted dose-response release. Unlike a
#     GEO series matrix - which upstream regenerates, and which script 12 in
#     the validation repo therefore guards with a size FLOOR - a changed byte
#     count here means a different file, not a cosmetic rebuild.
DIR_GDSC <- here::here("data", "raw", "gdsc")
GDSC_RELEASE <- "8.5 (fitted 27Oct23)"

GDSC_FILES <- c(
  gdsc2      = "GDSC2_fitted_dose_response_27Oct23.xlsx",
  gdsc1      = "GDSC1_fitted_dose_response_27Oct23.xlsx",
  compounds  = "screened_compounds_rel_8.5.csv",
  model_list = "model_list_20260921.csv",
  cnv        = "cnv_summary_20260804.csv",
  mutations  = "mutations_summary_20260724.csv")
GDSC_BYTES <- c(gdsc2 = 22112537, gdsc1 = 30348857, compounds = 46414,
                model_list = 957555, cnv = 70156911, mutations = 2241637)

.gdsc_help <- paste0(
  "\n\ndata/raw/gdsc is gitignored and is NOT on origin. Re-download it with ",
  "the curl commands in data/raw/gdsc_README.md, which also carries the ",
  "byte size and SHA-256 of every file.\n",
  "NOTE: cancerrxgene.org/downloads/bulk_download now returns HTTP 410. The ",
  "current host is https://cmp.cog.sanger.ac.uk/download/ . The DATA has not ",
  "moved on - release 8.5 is still current - but every pre-2026 URL for it ",
  "is stale.")

for (k in names(GDSC_FILES)) {
  fp <- file.path(DIR_GDSC, GDSC_FILES[[k]])
  if (!file.exists(fp)) {
    stop("data/raw/gdsc/", GDSC_FILES[[k]], " is absent.", .gdsc_help,
         call. = FALSE)
  }
  if (!identical(as.numeric(file.size(fp)), as.numeric(GDSC_BYTES[[k]]))) {
    stop(GDSC_FILES[[k]], " is ", file.size(fp), " bytes, not the ",
         GDSC_BYTES[[k]], " recorded in data/raw/gdsc_README.md. These files ",
         "are STATIC; a changed size means a different file.", .gdsc_help,
         call. = FALSE)
  }
}
message("0b. GDSC ", GDSC_RELEASE, ": all ", length(GDSC_FILES),
        " files present and byte-exact")

PATH_GDSC2 <- file.path(DIR_GDSC, GDSC_FILES[["gdsc2"]])
PATH_GDSC1 <- file.path(DIR_GDSC, GDSC_FILES[["gdsc1"]])
PATH_CNV   <- file.path(DIR_GDSC, GDSC_FILES[["cnv"]])
PATH_MUT   <- file.path(DIR_GDSC, GDSC_FILES[["mutations"]])

# =============================================================================
# 0c. CONSTANTS AND THE RULES. Fixed here, before any number.
# =============================================================================
# --- population, and the two DECLARED STOP RULES ----------------------------
HAEM_LINEAGES <- c("Lymphoid", "Myeloid")   # BCL2-dependent; would dominate
MIN_JOIN_RATE <- 0.90                       # stop below this
MIN_SOLID     <- 300L                       # stop below this
MIN_LINEAGE   <- 10L                        # lineage kept as a factor level
LINEAGE_FOCUS <- "Breast"
MIN_FIT_N     <- 40L                        # a cell with fewer is not fitted

# --- the exposure -----------------------------------------------------------
ARM_PRIMARY <- "OXPHOS subunits"   # the narrow set, NOT the umbrella, which
                                   # carries assembly factors. CLAUDE.md.
RULER_PRIMARY <- "ox_gsva"
RULERS        <- c("ox_gsva", "ox_ppd", "ox_lvl", "ox_rel")
# THE RULER PANEL IS A SENSITIVITY ON THE PRIMARY MODEL ONLY, and it is an
# addition made at script-writing time, before any fit, recorded as such. The
# declaration fixed one OX. CLAUDE.md trap 5 says the four instruments disagree
# by a lot and that instrument choice must be reported or justified; running
# all four on the PRIMARY estimand answers that WITHOUT turning the hierarchy
# into a 4 x 7 x 3 x 7 grid, which trap 5's own sibling rule forbids. Every
# other model in the hierarchy is fitted on RULER_PRIMARY alone.

# --- the MYC estimators. M_a primary, two sensitivities. Declaration s.9. ----
# Cell lines carry copy number, so MYC dose is computable here in a way the
# neoadjuvant cohorts could not manage. All three are reported.
MYC_ESTIMATORS <- c("myc_ma", "myc_expr", "myc_cn8q24")
MYC_PRIMARY    <- "myc_ma"

# --- THE RB1 RULE. The declared confounder. Fixed before any fit. -----------
# Varkaris et al., Nat Commun 2025 (doi 10.1038/s41467-025-60238-x): RB1
# alteration is the strongest sensitiser to navitoclax in GDSC. MYC-high lines
# may be enriched for it, so an unadjusted OX or MY effect could be RB1.
#
# RB1 loss is TRUE if EITHER holds, on Cell Model Passports, keyed on
# model_id, which IS the GDSC join key:
#   (i)  cnv_summary cn_category == "Deletion" from any source, OR
#   (ii) mutations_summary effect in {nonsense, ess_splice, frameshift}.
# MISSENSE IS EXCLUDED: RB1 missense is not reliably loss of function, and
# every RB1 row in that table is flagged cancer_driver == "t", so the driver
# flag does not discriminate. "Loss" (one copy) is EXCLUDED from (i) and is
# NOT the same value as "Deletion"; both exist in the column.
RB1_CN_LOSS  <- "Deletion"
RB1_MUT_LOSS <- c("nonsense", "ess_splice", "frameshift")

# --- readouts, and the sign convention -------------------------------------
# ALL THREE POINT THE SAME WAY:
#   GDSC LN_IC50  lower = more sensitive
#   GDSC AUC      lower = more sensitive
#   PRISM LFC     more negative = more killing
# SO A NEGATIVE COEFFICIENT ON OX SUPPORTS THE PREDICTION ON EVERY READOUT.
# This sentence is repeated in the result note so nobody reads the forest
# backwards.
SIGN_NOTE <- paste(
  "LN_IC50 lower = more sensitive; AUC lower = more sensitive; PRISM LFC more",
  "negative = more killing. A NEGATIVE OX coefficient supports the prediction",
  "on all three.")

# --- the drugs. Declaration s.4. Availability recorded in the data README. --
# Spellings are NOT assumed: section 6 greps unique(DRUG_NAME) in the loaded
# tables and records what it matched. Two traps that are really in the data -
# GDSC1 spells venetoclax "Venotoclax", and GDSC2 spells ABT-737 "ABT737" and
# omits it from screened_compounds_rel_8.5.csv entirely.
DRUG_PATTERNS <- c(
  navitoclax  = "^navitoclax$|^abt-?263$",
  venetoclax  = "^venetoclax$|^venotoclax$|^abt-?199$",
  `wehi-539`  = "^wehi-?539$",
  `abt-737`   = "^abt-?737$",
  s63845      = "^s-?63845$",
  azd5991     = "^azd5991$",
  sabutoclax  = "^sabutoclax$",
  obatoclax   = "^obatoclax")
DRUG_PRIMARY <- "navitoclax"
DRUG_CONTROL <- "venetoclax"
# Declared and ABSENT FROM ALL THREE PLATFORMS - recorded so that their absence
# from the results is not read as a failure to find them:
DRUGS_ABSENT <- c("A-1331852", "A-1155463")

# --- THE ESTIMAND HIERARCHY. NOT REORDERED AFTER THE ESTIMATES. -------------
# Declaration section 3. If a lower-ranked model fires and the primary does
# not, that is REPORTED AS SUCH and is not promoted.
ESTIMANDS <- tibble::tribble(
  ~rank,           ~model,        ~formula_rhs,                  ~read,
  "PRIMARY",       "primary",     "OX + MY + lineage",           "OX",
  "SECONDARY",     "secondary",   "OX + lineage | MYC tertile",  "OX",
  "TERTIARY",      "tertiary",    "OX * MY + lineage",           "OX:MY",
  "contrast",      "contrast",    "MY + lineage",                "MY",
  "companion",     "comp_prolif", "OX + MY + PROLIF + lineage",  "OX",
  "companion",     "comp_rb1",    "OX + MY + RB1 + lineage",     "OX",
  "mediation only","mediation",   "OX + MY + z(BCL2L1) + lineage","OX",
  "positive control","rb1_control","OX + MY + RB1 + lineage",    "RB1")

# A NULL ON THE PRODUCT TERM ALONE MEANS NOTHING. That estimand has already
# failed in Block C, Block B, Block G, H4 and N2. A sixth null is not a finding
# and is not presented as one. Declaration s.8 guardrail 3.

# =============================================================================
# 1. Gene sets - REBUILT FROM THIS REPO'S PINNED MITOCARTA. E17/E18 section 1.
# =============================================================================
message("\n1. gene sets, rebuilt from ", basename(PATH_MITOCARTA))

suppressWarnings({
  mc_inv <- readxl::read_xls(PATH_MITOCARTA, sheet = 2)
  mc_pw  <- readxl::read_xls(PATH_MITOCARTA, sheet = 4)
})
# Sheet 4 ends with 5 blank padding rows; leaving them in inflates every
# pathway by 5 phantom "NA" genes, silently. CLAUDE.md read-time trap.
mc_pw <- mc_pw %>% dplyr::filter(!is.na(MitoPathway))
stopifnot(length(unique(mc_inv$Symbol)) == EXPECT_MITOCARTA_ALL,
          nrow(mc_pw) == 149L)

.split_genes <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))
MTDNA_PATHWAY <- "mtDNA-encoded OXPHOS subunits"
mito_raw    <- stats::setNames(lapply(mc_pw$Genes, .split_genes),
                               mc_pw$MitoPathway)
MTDNA_GENES <- grep("^MT-", unique(mc_inv$Symbol), value = TRUE)
stopifnot(length(MTDNA_GENES) == 13L, !MTDNA_PATHWAY %in% names(mito_raw))
# Standing convention: mtDNA-encoded genes are held separately and NEVER
# pooled with the nuclear-encoded subunits.
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
  lapply(ARM_PATHWAYS,
         function(p) unique(unlist(mito_paths[p], use.names = FALSE))),
  names(ARM_PATHWAYS))

snap_arms <- readRDS(PATH_TCGA_MITO)$arm_sets
same_arm  <- vapply(names(ARM_SETS),
                    function(a) setequal(ARM_SETS[[a]], snap_arms[[a]]),
                    logical(1))
if (!all(same_arm)) {
  stop("arm(s) rebuilt from this repo's MitoCarta differ from the snapshot's ",
       "arm_sets: ", paste(names(same_arm)[!same_arm], collapse = ", "),
       call. = FALSE)
}
message("   18 arms rebuilt, all set-identical to the snapshot's arm_sets")

sd_      <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL
MYC_SET  <- sd_$myc_sets[[MYC_REF]]        # FELSHER__MITOSTRIP - M_a, 61 genes
PROLIF   <- sd_$cov_sets$PROLIF_DISJOINT   # 318
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL,
          length(MYC_SET)  == EXPECT_FELSHER_STRIP,
          length(PROLIF)   == 318L)
message("   MYC = ", MYC_REF, " (", length(MYC_SET), " genes), ",
        "PROLIF_DISJOINT (", length(PROLIF), "), ",
        ARM_PRIMARY, " (", length(ARM_SETS[[ARM_PRIMARY]]), ")")

# =============================================================================
# 2. Cell lines, the join, and the two declared stop rules
# =============================================================================
message("\n2. cell lines and the join")

MODEL <- data.table::fread(PATH_MODEL, data.table = FALSE)
for (need in c("ModelID", "SangerModelID", "OncotreeLineage")) {
  if (!need %in% names(MODEL)) {
    stop("Model.csv has no ", need, " column. DepMap renamed it; fix the ",
         "name rather than guessing.", call. = FALSE)
  }
}
MODEL$SangerModelID[MODEL$SangerModelID == ""] <- NA_character_
message("   Model.csv: ", nrow(MODEL), " models, ",
        sum(!is.na(MODEL$SangerModelID)), " carrying a SangerModelID")

# --- GDSC dose response. BOTH screens: GDSC1 is where ABT-737 might have been,
#     and it carries an independent navitoclax screen. One row per curve fit,
#     NOT one row per line - a line screened in both files appears in both, and
#     DATASET is what distinguishes them.
.read_gdsc <- function(path, label) {
  d <- suppressMessages(readxl::read_excel(path, guess_max = 5000))
  need <- c("SANGER_MODEL_ID", "DRUG_NAME", "DRUG_ID", "LN_IC50", "AUC",
            "CELL_LINE_NAME")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop(basename(path), " lacks column(s): ", paste(miss, collapse = ", "),
         ". Columns are: ", paste(names(d), collapse = ", "), call. = FALSE)
  }
  d$DATASET_FILE <- label
  message("   ", label, ": ", nrow(d), " curve fits, ",
          length(unique(d$DRUG_NAME)), " drugs, ",
          length(unique(d$SANGER_MODEL_ID)), " lines")
  d
}
G2 <- .read_gdsc(PATH_GDSC2, "GDSC2")
G1 <- .read_gdsc(PATH_GDSC1, "GDSC1")
GD <- dplyr::bind_rows(G2, G1)

gdsc_lines <- unique(GD$SANGER_MODEL_ID)
join_rate  <- mean(gdsc_lines %in% MODEL$SangerModelID)
message("   GDSC lines: ", length(gdsc_lines), "; joining to DepMap: ",
        sum(gdsc_lines %in% MODEL$SangerModelID), " (",
        round(100 * join_rate, 1), "%)")
if (join_rate < MIN_JOIN_RATE) {
  stop("the GDSC-to-DepMap join rate is ", round(100 * join_rate, 1),
       "%, below the declared floor of ", round(100 * MIN_JOIN_RATE),
       "%. The declaration says STOP AND REPORT rather than proceed on a ",
       "join nobody has explained.", call. = FALSE)
}

# --- solid only. Haematological lines are BCL2-dependent and would dominate
#     any BCL2-family signal. Declaration s.5.
xw <- MODEL %>%
  dplyr::filter(!is.na(SangerModelID), !is.na(OncotreeLineage),
                OncotreeLineage != "",
                !OncotreeLineage %in% HAEM_LINEAGES) %>%
  dplyr::transmute(ModelID, SangerModelID, lineage = OncotreeLineage)
solid_ids <- intersect(xw$SangerModelID, gdsc_lines)
message("   solid GDSC lines joining to DepMap: ", length(solid_ids),
        " (excluded lineages: ", paste(HAEM_LINEAGES, collapse = ", "), ")")
if (length(solid_ids) < MIN_SOLID) {
  stop("only ", length(solid_ids), " solid lines survive the join, below the ",
       "declared floor of ", MIN_SOLID, ". STOP AND REPORT.", call. = FALSE)
}

# =============================================================================
# 3. Expression. E17/E18's reader, unchanged.
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
           paste(utils::head(bad, 5), collapse = ", "), ". Refusing to guess.",
           call. = FALSE)
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
  # Deliberately NOT IsDefaultEntryForMC - that is the model-CONDITION default,
  # a different and larger set. See E17 section 3.
  keep <- .is_true(m[[flag_col]], flag_col)
  if (!any(keep)) stop("no rows have ", flag_col, " true.", call. = FALSE)
  message("   ", sum(keep), " default rows of ", nrow(m), " (flag: ",
          flag_col, ")")
  ids <- as.character(m$ModelID[keep])
  if (anyDuplicated(ids)) {
    stop(sum(duplicated(ids)), " ModelID(s) still duplicated after filtering. ",
         "Do not de-duplicate silently - find out why.", call. = FALSE)
  }
  META <- c("V1", "", "ProfileID", "SequencingID", "ModelConditionID",
            "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel",
            "is_default_entry")
  gene_cols <- setdiff(names(m), META)
  num_ok <- vapply(m[gene_cols], is.numeric, logical(1))
  if (!all(num_ok)) {
    stop("non-numeric column(s) survived the metadata filter: ",
         paste(utils::head(gene_cols[!num_ok], 5), collapse = ", "),
         ". Add them to META in .read_expression().", call. = FALSE)
  }
  M <- as.matrix(m[keep, gene_cols, drop = FALSE])
  rownames(M) <- ids
  # Columns are "SYMBOL (ENTREZ)". Strip the suffix; keep the first occurrence
  # of a duplicated symbol; drop anything that strips to nothing.
  colnames(M) <- trimws(sub("\\s*\\(\\d+\\)$", "", colnames(M)))
  M <- M[, nzchar(colnames(M)) & !is.na(colnames(M)), drop = FALSE]
  dup <- duplicated(colnames(M))
  if (any(dup)) {
    message("   dropped ", sum(dup), " duplicated symbol column(s)")
    M <- M[, !dup, drop = FALSE]
  }
  M
}

EXPR <- .read_expression(PATH_EXPR)
message("   expression: ", nrow(EXPR), " models x ", ncol(EXPR), " genes")
if (max(EXPR, na.rm = TRUE) > 40) {
  stop("expression matrix max is ", round(max(EXPR, na.rm = TRUE), 1),
       ", which is not log2(TPM+1). Check the file.", call. = FALSE)
}

# --- THE ANALYSIS COHORT. Fixed here, and every score below is relative to it.
# ATTRITION IS REPORTED, NOT ABSORBED. The GDSC panel is largely Sanger lines
# and DepMap expression does not cover all of them, so requiring expression
# costs real lines. The count is printed and saved so a reader can see what the
# cohort is a cohort OF, rather than meeting only the surviving number.
n_solid_joined <- length(solid_ids)
xw <- xw %>% dplyr::filter(SangerModelID %in% solid_ids,
                           ModelID %in% rownames(EXPR))
# Two DepMap models can carry the same SangerModelID. Keep one, deterministically
# by sorted ModelID, rather than letting the join multiply rows silently.
dup_sanger <- xw$SangerModelID[duplicated(xw$SangerModelID)]
if (length(dup_sanger)) {
  message("   ", length(unique(dup_sanger)), " SangerModelID(s) map to more ",
          "than one DepMap ModelID; keeping the first by sorted ModelID")
  xw <- xw %>% dplyr::arrange(SangerModelID, ModelID) %>%
    dplyr::distinct(SangerModelID, .keep_all = TRUE)
}
# De-duplication happens AFTER the expression filter deliberately: where a
# Sanger line carries two DepMap models, the one with expression is the one
# worth keeping.
n_with_expr <- nrow(xw)
message("   of ", n_solid_joined, " solid joined lines, ", n_with_expr,
        " carry DepMap ", DEPMAP_RELEASE, " expression (lost ",
        n_solid_joined - n_with_expr, ")")
# Lineage as a factor: levels with fewer than MIN_LINEAGE lines are DROPPED
# from the cohort, as B4's pan-cancer fit did, rather than collapsed into an
# "Other" bucket that mixes unrelated tissues behind one coefficient.
lin_n  <- table(xw$lineage)
keep_l <- names(lin_n)[lin_n >= MIN_LINEAGE]
dropped_lineages <- sort(setdiff(names(lin_n), keep_l))
xw <- xw %>% dplyr::filter(lineage %in% keep_l)
# COHORT and sid are PARALLEL and both UNNAMED: COHORT indexes the DepMap
# matrices (ACH ids) and sid indexes GDSC and Cell Model Passports (SIDM ids).
# Keeping them as two aligned vectors rather than one named vector means the
# identity checks below compare values and not attributes.
COHORT <- xw$ModelID
sid    <- xw$SangerModelID
stopifnot(length(COHORT) == length(sid), !anyDuplicated(COHORT),
          !anyDuplicated(sid))
message("   ANALYSIS COHORT: ", length(COHORT), " solid lines with expression ",
        "and a GDSC screen, in ", length(keep_l), " lineages of >= ",
        MIN_LINEAGE)
if (length(dropped_lineages)) {
  message("   dropped ", length(dropped_lineages), " lineage(s) below ",
          MIN_LINEAGE, ": ", paste(dropped_lineages, collapse = ", "))
}
if (length(COHORT) < MIN_SOLID) {
  stop("only ", length(COHORT), " lines survive to the analysis cohort, below ",
       "the declared floor of ", MIN_SOLID, ". STOP AND REPORT.", call. = FALSE)
}
n_focus <- sum(xw$lineage == LINEAGE_FOCUS)
message("   ", LINEAGE_FOCUS, ": ", n_focus,
        " lines - the stratum of interest, and UNDERPOWERED as declared")

# =============================================================================
# 4. Scoring. ONE GSVA CALL ON THE ANALYSIS COHORT. E17/E18 section 4.
# =============================================================================
# GSVA IS COHORT-RELATIVE. All sets are scored in the SAME call so the gene
# universe and the half-matrix pins are shared, and the cohort is the FINAL
# analysis cohort, not a superset later subsetted - subsetting after scoring
# would leave scores relative to a population that is not the one fitted.
message("\n4. scoring the analysis cohort (one run, all sets)")

GSVA_MIN_SET <- 3L
MIN_SET_FRAC <- 0.80

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
.check_positive <- function(S) {
  bad <- which(!is.finite(S) | S <= 0, arr.ind = TRUE)
  if (nrow(bad)) {
    lab <- utils::head(sprintf("%s / %s", rownames(S)[bad[, 1]],
                               colnames(S)[bad[, 2]]), 6L)
    stop("mitoPPS: ", nrow(bad), " arm x line pathway mean(s) are zero or ",
         "non-finite. First few: ", paste(lab, collapse = "; "),
         ". Do NOT add a floor silently.", call. = FALSE)
  }
}
# E16's composite: mean of per-gene z across samples, on a LOG matrix.
.comp <- function(genes, E) {
  M <- E[intersect(genes, rownames(E)), , drop = FALSE]
  v <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}

E_LOG <- t(EXPR[COHORT, , drop = FALSE])          # log2(TPM+1), AS SUPPLIED
E_LOG <- E_LOG[stats::complete.cases(E_LOG), , drop = FALSE]
colnames(E_LOG) <- COHORT

sets_defined <- c(list(MYC = MYC_SET, PROLIF = PROLIF), ARM_SETS)
sets_gsva    <- lapply(sets_defined, function(s) intersect(s, rownames(E_LOG)))
coverage <- tibble::tibble(set = names(sets_defined),
                           n_defined = lengths(sets_defined),
                           n_present = lengths(sets_gsva)) %>%
  dplyr::mutate(frac = round(n_present / n_defined, 3))
lost <- coverage$set[coverage$frac < MIN_SET_FRAC]
if (length(lost)) {
  stop("gene set(s) below ", MIN_SET_FRAC, " coverage in the cohort: ",
       paste(lost, collapse = ", "),
       ". That is a symbol-harmonisation failure, not natural attrition.",
       call. = FALSE)
}
tiny <- coverage$set[coverage$n_present < GSVA_MIN_SET]
if (length(tiny)) {
  stop("gene set(s) with fewer than ", GSVA_MIN_SET, " genes: ",
       paste(tiny, collapse = ", "), call. = FALSE)
}
message("   set coverage, n_present / n_defined:")
coverage %>%
  dplyr::filter(set %in% c("MYC", "PROLIF", ARM_PRIMARY)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# GSVA's API changed at 1.50. Branch on the installed version rather than
# assuming - a wrong guess here is a stop, not a silently different score.
GSVA_VERSION <- as.character(utils::packageVersion("GSVA"))
if (utils::compareVersion(GSVA_VERSION, "1.50") >= 0) {
  GS <- GSVA::gsva(GSVA::gsvaParam(exprData = E_LOG, geneSets = sets_gsva,
                                   kcdf = "Gaussian", minSize = GSVA_MIN_SET,
                                   maxSize = Inf), verbose = FALSE)
} else {
  GS <- GSVA::gsva(E_LOG, sets_gsva, method = "gsva", kcdf = "Gaussian",
                   min.sz = GSVA_MIN_SET, verbose = FALSE)
}
miss <- setdiff(names(sets_gsva), rownames(GS))
if (length(miss)) {
  stop("GSVA silently dropped ", length(miss), " set(s): ",
       paste(miss, collapse = ", "), call. = FALSE)
}
message("   GSVA ", GSVA_VERSION, ": ", nrow(GS), " sets x ", ncol(GS), " lines")

E_LIN <- 2^E_LOG - 1                 # mitoPPS wants LINEAR
E_LIN[E_LIN < 0] <- 0                # floating point on an exact zero only
S_LIN <- .pathway_means(E_LIN, lapply(ARM_SETS,
                                      function(s) intersect(s, rownames(E_LOG))))
.check_positive(S_LIN)
PPS <- .mitopps_universe(S_LIN)

LVL <- t(vapply(ARM_SETS, function(s) .comp(s, E_LOG), numeric(ncol(E_LOG))))
REL <- t(vapply(ARM_SETS,
                function(s) .comp(s, E_LOG) - .comp(setdiff(MITO_ALL, s), E_LOG),
                numeric(ncol(E_LOG))))
rownames(LVL) <- rownames(REL) <- names(ARM_SETS)

SC <- list(ox_gsva = t(apply(GS[names(ARM_SETS), , drop = FALSE], 1, .z)),
           ox_ppd  = t(apply(PPS, 1, .z)),
           ox_lvl  = t(apply(LVL, 1, .z)),
           ox_rel  = t(apply(REL, 1, .z)))
for (r in RULERS) stopifnot(identical(colnames(SC[[r]]), COHORT))

# =============================================================================
# 5. The MYC estimators, RB1, and the mediator
# =============================================================================
message("\n5. MYC estimators, RB1 and the mediator")

# 5.1 M_a - the primary. GSVA on FELSHER__MITOSTRIP, z-scored on this cohort.
myc_ma <- .z(GS["MYC", ])

# 5.2 MYC transcript. log2(TPM+1) as supplied. CLAUDE.md trap 4 says MYC mRNA
#     is not MYC activity and that the difference is TOTAL - the gap is itself
#     a result, so both are reported and neither is presented as the other.
if (!"MYC" %in% rownames(E_LOG)) {
  stop("MYC is not a row of the expression matrix.", call. = FALSE)
}
myc_expr <- .z(E_LOG["MYC", ])

# 5.3 8q24 copy number, from Cell Model Passports. Keyed on model_id, which IS
#     SANGER_MODEL_ID. Heavily right-skewed (median 3, max 242), so log2.
#     A model may carry both a Sanger and a Broad row; Sanger is preferred
#     deterministically, falling back to Broad, rather than averaging two
#     pipelines into a number that is neither.
CNV <- data.table::fread(PATH_CNV, data.table = FALSE)
for (need in c("model_id", "symbol", "total_copy_number", "cn_category",
               "source")) {
  if (!need %in% names(CNV)) {
    stop("cnv_summary lacks column ", need, ". Columns: ",
         paste(names(CNV), collapse = ", "), call. = FALSE)
  }
}
.one_per_model <- function(d) {
  d %>% dplyr::mutate(.pref = dplyr::if_else(source == "Sanger", 1L, 2L)) %>%
    dplyr::arrange(model_id, .pref) %>%
    dplyr::distinct(model_id, .keep_all = TRUE)
}
cn_myc <- CNV %>% dplyr::filter(symbol == "MYC") %>% .one_per_model()
cn_rb1 <- CNV %>% dplyr::filter(symbol == "RB1") %>% .one_per_model()
rm(CNV)

myc_cn_raw <- cn_myc$total_copy_number[match(sid, cn_myc$model_id)]
myc_cn8q24 <- .z(log2(pmax(myc_cn_raw, 0.5)))   # 0.5 floor: log2(0) is -Inf
names(myc_cn8q24) <- COHORT
message("   MYC copy number present for ", sum(!is.na(myc_cn_raw)), " of ",
        length(sid), " lines")

MYC_MAT <- rbind(myc_ma = myc_ma, myc_expr = myc_expr, myc_cn8q24 = myc_cn8q24)
colnames(MYC_MAT) <- COHORT
stopifnot(identical(rownames(MYC_MAT), MYC_ESTIMATORS))
message("   MYC estimators, pairwise Spearman on this cohort:")
print(round(stats::cor(t(MYC_MAT), method = "spearman",
                       use = "pairwise.complete.obs"), 3))

# 5.4 RB1 loss. The rule fixed in 0c, applied - not decided here.
MUT <- data.table::fread(PATH_MUT, data.table = FALSE)
for (need in c("model_id", "gene_symbol", "effect")) {
  if (!need %in% names(MUT)) {
    stop("mutations_summary lacks column ", need, call. = FALSE)
  }
}
rb1_mut_ids <- unique(MUT$model_id[MUT$gene_symbol == "RB1" &
                                     MUT$effect %in% RB1_MUT_LOSS])
rb1_del_ids <- unique(cn_rb1$model_id[cn_rb1$cn_category %in% RB1_CN_LOSS])
rm(MUT)
RB1 <- as.integer(sid %in% union(rb1_mut_ids, rb1_del_ids))
names(RB1) <- COHORT
message("   RB1 loss: ", sum(RB1), " of ", length(RB1), " lines (",
        sum(sid %in% rb1_del_ids), " by ", RB1_CN_LOSS, ", ",
        sum(sid %in% rb1_mut_ids), " by truncating mutation)")

# 5.5 PROLIF, and the mediator. BCL2L1 IS A MEDIATOR, NOT A CONFOUNDER: if the
#     model is right the path is OXPHOS -> more BCL-XL -> more dependency, so
#     adjusting for it is OVER-adjustment. It is fitted ONCE and reported
#     separately. It NEVER enters the primary. Declaration s.7.
PROLIF_SC <- .z(GS["PROLIF", ])
if (!"BCL2L1" %in% rownames(E_LOG)) {
  stop("BCL2L1 is not a row of the expression matrix.", call. = FALSE)
}
BCL2L1 <- .z(E_LOG["BCL2L1", ])

# --- the MYC tertiles, on the PRIMARY estimator, over the analysis cohort ----
brk <- stats::quantile(MYC_MAT[MYC_PRIMARY, ], c(1/3, 2/3), na.rm = TRUE)
myc_tertile <- cut(MYC_MAT[MYC_PRIMARY, ], c(-Inf, brk, Inf),
                   labels = c("T1_low", "T2", "T3_high"))
names(myc_tertile) <- COHORT
message("   MYC tertiles on ", MYC_PRIMARY, ":")
print(table(myc_tertile))

COV <- tibble::tibble(
  ModelID = COHORT, SangerModelID = sid,
  lineage = xw$lineage[match(COHORT, xw$ModelID)],
  MY_ma = MYC_MAT["myc_ma", ], MY_expr = MYC_MAT["myc_expr", ],
  MY_cn8q24 = MYC_MAT["myc_cn8q24", ],
  PROLIF = PROLIF_SC, RB1 = RB1, BCL2L1 = BCL2L1,
  myc_tertile = as.character(myc_tertile))
for (r in RULERS) COV[[r]] <- SC[[r]][ARM_PRIMARY, COHORT]

# =============================================================================
# 6. The drug endpoints. Spellings are matched, never assumed.
# =============================================================================
message("\n6. drug endpoints")

.match_drugs <- function(names_seen, label) {
  nm <- unique(names_seen)
  out <- lapply(names(DRUG_PATTERNS), function(k) {
    hit <- nm[grepl(DRUG_PATTERNS[[k]], trimws(nm), ignore.case = TRUE)]
    if (!length(hit)) return(NULL)
    tibble::tibble(drug = k, matched = hit, source = label)
  })
  dplyr::bind_rows(out)
}
drug_map <- dplyr::bind_rows(
  .match_drugs(GD$DRUG_NAME[GD$DATASET_FILE == "GDSC2"], "GDSC2"),
  .match_drugs(GD$DRUG_NAME[GD$DATASET_FILE == "GDSC1"], "GDSC1"))
message("   BH3 mimetics matched in GDSC, by declared name:")
drug_map %>% as.data.frame() %>% print(row.names = FALSE)
if (!DRUG_PRIMARY %in% drug_map$drug) {
  stop(DRUG_PRIMARY, " is absent from the release. The declaration says STOP ",
       "AND REPORT rather than substitute something else.", call. = FALSE)
}
if (!DRUG_CONTROL %in% drug_map$drug) {
  stop(DRUG_CONTROL, " is absent from the release, so the specificity control ",
       "cannot be fitted. Without it a navitoclax effect cannot be attributed ",
       "to BCL-XL. STOP AND REPORT.", call. = FALSE)
}
message("   declared but ABSENT FROM ALL PLATFORMS: ",
        paste(DRUGS_ABSENT, collapse = ", "),
        " - recorded so their absence is not read as a failure")

# --- GDSC long table: one row per line per drug per screen -------------------
GD_LONG <- GD %>%
  dplyr::inner_join(drug_map %>% dplyr::select(drug, matched, source),
                    by = c("DRUG_NAME" = "matched",
                           "DATASET_FILE" = "source")) %>%
  dplyr::filter(SANGER_MODEL_ID %in% sid) %>%
  dplyr::transmute(screen = DATASET_FILE, drug, SangerModelID = SANGER_MODEL_ID,
                   LN_IC50 = as.numeric(LN_IC50), AUC = as.numeric(AUC)) %>%
  tidyr::pivot_longer(c(LN_IC50, AUC), names_to = "readout", values_to = "Y") %>%
  dplyr::filter(is.finite(Y))
# A line screened twice against the same drug in the same screen would silently
# double-weight. Collapse to the median and record how often it happens.
dup_curves <- GD_LONG %>%
  dplyr::count(screen, drug, readout, SangerModelID, name = "n_fits") %>%
  dplyr::filter(n_fits > 1L)
if (nrow(dup_curves)) {
  message("   ", nrow(dup_curves), " line x drug x readout cell(s) carry more ",
          "than one curve fit; collapsing to the median")
  GD_LONG <- GD_LONG %>%
    dplyr::group_by(screen, drug, readout, SangerModelID) %>%
    dplyr::summarise(Y = stats::median(Y), .groups = "drop")
}

# --- PRISM, the third PLATFORM. Not a third summary of the same curve. -------
# The matrix is COMPOUNDS x CELL LINES, the transpose of the DepMap omics
# files; read without transposing, the line intersection is empty and the arm
# reports "0 lines" rather than failing. Keyed on ACH ModelID, so it joins to
# the cohort directly and not through SANGER_MODEL_ID.
PRISM_LONG <- NULL
prism_map  <- NULL
if (file.exists(PATH_PRISM_M) && file.exists(PATH_PRISM_C)) {
  PC <- data.table::fread(PATH_PRISM_C, data.table = FALSE)
  if (!all(c("IDs", "Drug.Name") %in% names(PC))) {
    stop("the PRISM compound list lacks IDs or Drug.Name. Columns: ",
         paste(names(PC), collapse = ", "), call. = FALSE)
  }
  prism_map <- .match_drugs(PC$Drug.Name, "PRISM")
  if (nrow(prism_map)) {
    message("   BH3 mimetics matched in PRISM ", PRISM_RELEASE, ": ",
            paste(sort(unique(prism_map$drug)), collapse = ", "))
    brd <- PC %>%
      dplyr::inner_join(prism_map %>% dplyr::select(drug, matched),
                        by = c("Drug.Name" = "matched")) %>%
      dplyr::select(drug, IDs, screen, dose) %>% dplyr::distinct()
    PM <- data.table::fread(PATH_PRISM_M, data.table = FALSE)
    rn <- as.character(PM[[1]])
    keep_r <- rn %in% brd$IDs
    keep_c <- intersect(names(PM), COHORT)
    if (any(keep_r) && length(keep_c)) {
      sub <- PM[keep_r, c(names(PM)[1], keep_c), drop = FALSE]
      PRISM_LONG <- sub %>%
        tidyr::pivot_longer(-1, names_to = "ModelID", values_to = "Y") %>%
        dplyr::rename(IDs = 1) %>%
        dplyr::inner_join(brd %>% dplyr::select(drug, IDs), by = "IDs") %>%
        dplyr::filter(is.finite(Y)) %>%
        # One compound can appear under several BRD ids; median per line.
        dplyr::group_by(drug, ModelID) %>%
        dplyr::summarise(Y = stats::median(Y), .groups = "drop") %>%
        dplyr::mutate(screen = "PRISM", readout = "PRISM_LFC")
      message("   PRISM: ", length(unique(PRISM_LONG$ModelID)),
              " cohort lines x ", length(unique(PRISM_LONG$drug)), " drugs")
    }
    rm(PM)
  }
} else {
  message("   PRISM files absent - the third platform is SKIPPED, as script ",
          "14 skips it. The two GDSC readouts are summaries of the SAME ",
          "curve and are not two platforms.")
}

# --- one endpoint table -----------------------------------------------------
ENDPOINTS <- dplyr::bind_rows(
  GD_LONG %>% dplyr::mutate(ModelID = COHORT[match(SangerModelID, sid)]) %>%
    dplyr::select(screen, drug, readout, ModelID, Y),
  if (!is.null(PRISM_LONG)) PRISM_LONG %>%
    dplyr::select(screen, drug, readout, ModelID, Y))
message("   endpoint table: ", nrow(ENDPOINTS), " line x drug x readout rows, ",
        length(unique(ENDPOINTS$drug)), " drugs, ",
        length(unique(paste(ENDPOINTS$screen, ENDPOINTS$readout))),
        " screen x readout combinations")

# =============================================================================
# 7. THE MODELS. Every one in the hierarchy, every drug, every readout.
# =============================================================================
# EVERY MODEL FITTED IS SAVED, including the ones that do not fire. The
# declaration lists them and the results object must contain all of them.
#
# SCALE OF THE REPORTED COEFFICIENT. Y is z-scored within each
# screen x drug x readout cell, and OX and MY are z-scored over the cohort.
# This is a PRESENTATION choice and not an estimand change: standardising a
# response and a predictor is a linear rescaling, so it leaves every t
# statistic, every p value and every SIGN exactly as they were, and changes
# only the units. It is done because a forest plot cannot otherwise put AUC,
# which is bounded [0, 1], beside LN_IC50, which spans about ten natural-log
# units. sd_Y is saved per cell so raw units are recoverable by multiplication.
message("\n7. fitting")

.fit_one <- function(d, rhs, read_term, label) {
  d <- d[stats::complete.cases(d[, intersect(c("Y", "OX", "MY", "PROLIF",
                                               "RB1", "BCL2L1", "lineage"),
                                             names(d))]), , drop = FALSE]
  if (nrow(d) < MIN_FIT_N) return(NULL)
  d$lineage <- droplevels(factor(d$lineage))
  # A single surviving lineage makes the factor rank-deficient rather than
  # informative; drop the term instead of letting lm() silently alias it.
  if (nlevels(d$lineage) < 2L) rhs <- sub("\\s*\\+\\s*lineage", "", rhs)
  if (length(unique(d$Y)) < 3L) return(NULL)
  sd_Y <- stats::sd(d$Y)
  if (!is.finite(sd_Y) || sd_Y == 0) return(NULL)
  d$Y <- (d$Y - mean(d$Y)) / sd_Y
  # RB1 with no variation in a stratum is a constant column, not a covariate.
  if ("RB1" %in% all.vars(stats::as.formula(paste("~", rhs))) &&
      length(unique(d$RB1)) < 2L) return(NULL)
  fm <- stats::as.formula(paste("Y ~", rhs))
  fit <- try(stats::lm(fm, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  cf <- summary(fit)$coefficients
  if (!read_term %in% rownames(cf)) return(NULL)
  est <- cf[read_term, "Estimate"]; se <- cf[read_term, "Std. Error"]
  ci <- stats::confint(fit, read_term, level = 0.95)
  tibble::tibble(
    n = nrow(d), n_lineages = nlevels(d$lineage), term = read_term,
    estimate = est, se = se, ci_lo = ci[1, 1], ci_hi = ci[1, 2],
    p = cf[read_term, "Pr(>|t|)"], sd_Y_raw = sd_Y,
    ci_excludes_0 = (ci[1, 1] > 0) | (ci[1, 2] < 0))
}

# The seven models of the hierarchy, as fitted objects. `strata` marks the one
# that is fitted inside a MYC tertile rather than over the cohort.
MODELS <- list(
  primary     = list(rhs = "OX + MY + lineage",            read = "OX",    strata = NA),
  secondary   = list(rhs = "OX + lineage",                 read = "OX",    strata = c("T1_low", "T3_high")),
  tertiary    = list(rhs = "OX * MY + lineage",            read = "OX:MY", strata = NA),
  contrast    = list(rhs = "MY + lineage",                 read = "MY",    strata = NA),
  comp_prolif = list(rhs = "OX + MY + PROLIF + lineage",   read = "OX",    strata = NA),
  comp_rb1    = list(rhs = "OX + MY + RB1 + lineage",      read = "OX",    strata = NA),
  mediation   = list(rhs = "OX + MY + BCL2L1 + lineage",   read = "OX",    strata = NA),
  # THE POSITIVE CONTROL. Same fit as comp_rb1; the RB1 term is read instead of
  # the OX term. Added after the first dry run showed the primary to be null,
  # and it changes no estimand: it reads a coefficient that comp_rb1 already
  # computed, and it is about RB1, not about OXPHOS.
  #
  # It exists because a null is worth much more when the pipeline can be shown
  # to detect a real association of the expected kind. Varkaris et al. found
  # RB1 alteration to be the strongest sensitiser to navitoclax in GDSC. If
  # RB1 lands NEGATIVE here, this cohort, these models and these endpoints can
  # find a genuine navitoclax sensitiser - and the OXPHOS null is then a
  # statement about OXPHOS. If RB1 is ALSO null, the null is a statement about
  # the analysis and must be written as one.
  rb1_control = list(rhs = "OX + MY + RB1 + lineage",      read = "RB1",   strata = NA))
stopifnot(setequal(names(MODELS), ESTIMANDS$model))

cells <- ENDPOINTS %>% dplyr::distinct(screen, drug, readout)
message("   ", nrow(cells), " screen x drug x readout cells x ",
        length(MODELS), " models, on ", length(RULERS),
        " rulers for the primary and ", RULER_PRIMARY, " for the rest")

fits <- list(); k <- 0L
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  dd <- ENDPOINTS %>%
    dplyr::filter(screen == ce$screen, drug == ce$drug,
                  readout == ce$readout) %>%
    dplyr::inner_join(COV, by = "ModelID")
  for (myc_est in MYC_ESTIMATORS) {
    my_col <- c(myc_ma = "MY_ma", myc_expr = "MY_expr",
                myc_cn8q24 = "MY_cn8q24")[[myc_est]]
    # The ruler panel is a sensitivity on the PRIMARY model only, and only on
    # the primary MYC estimator - a 4-ruler x 3-estimator cross would be the
    # grid CLAUDE.md forbids.
    for (ruler in RULERS) {
      if (ruler != RULER_PRIMARY && myc_est != MYC_PRIMARY) next
      for (mn in names(MODELS)) {
        m <- MODELS[[mn]]
        if (ruler != RULER_PRIMARY && mn != "primary") next
        d0 <- dd %>%
          dplyr::transmute(Y, OX = .data[[ruler]], MY = .data[[my_col]],
                           PROLIF, RB1, BCL2L1, lineage,
                           myc_tertile)
        strata <- if (all(is.na(m$strata))) NA_character_ else m$strata
        for (st in strata) {
          d1 <- if (is.na(st)) d0 else d0[d0$myc_tertile == st, , drop = FALSE]
          r <- .fit_one(d1, m$rhs, m$read, mn)
          if (is.null(r)) {
            r <- tibble::tibble(n = nrow(d1), n_lineages = NA_integer_,
                                term = m$read, estimate = NA_real_,
                                se = NA_real_, ci_lo = NA_real_,
                                ci_hi = NA_real_, p = NA_real_,
                                sd_Y_raw = NA_real_, ci_excludes_0 = NA)
          }
          k <- k + 1L
          fits[[k]] <- dplyr::bind_cols(
            tibble::tibble(rank = ESTIMANDS$rank[match(mn, ESTIMANDS$model)],
                           model = mn, screen = ce$screen, drug = ce$drug,
                           readout = ce$readout, ruler = ruler,
                           myc_estimator = myc_est,
                           stratum = if (is.na(st)) "all" else st,
                           rhs = m$rhs), r)
        }
      }
    }
  }
}
COEF <- dplyr::bind_rows(fits)
message("   ", nrow(COEF), " model fits recorded, of which ",
        sum(is.na(COEF$estimate)), " did not estimate (n below ", MIN_FIT_N,
        ", or a constant column)")

# --- the primary, printed --------------------------------------------------
message("\n   THE PRIMARY: ", ARM_PRIMARY, " coefficient at fixed MYC, ",
        RULER_PRIMARY, ", ", MYC_PRIMARY)
message("   ", SIGN_NOTE)
COEF %>%
  dplyr::filter(model == "primary", ruler == RULER_PRIMARY,
                myc_estimator == MYC_PRIMARY) %>%
  dplyr::transmute(drug, screen, readout, n,
                   estimate = round(estimate, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3), p = signif(p, 2),
                   ci_excludes_0) %>%
  dplyr::arrange(drug, screen, readout) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   THE SPECIFICITY CONTRAST - if venetoclax moves with navitoclax, ",
        "the effect is NOT BCL-XL:")
COEF %>%
  dplyr::filter(model == "primary", ruler == RULER_PRIMARY,
                myc_estimator == MYC_PRIMARY,
                drug %in% c(DRUG_PRIMARY, DRUG_CONTROL)) %>%
  dplyr::transmute(drug, screen, readout, estimate = round(estimate, 3),
                   ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                   ci_excludes_0) %>%
  dplyr::arrange(screen, readout, drug) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   THE SECONDARY - the conditionality claim, by MYC tertile:")
COEF %>%
  dplyr::filter(model == "secondary", myc_estimator == MYC_PRIMARY,
                drug %in% c(DRUG_PRIMARY, DRUG_CONTROL)) %>%
  dplyr::transmute(drug, screen, readout, stratum, n,
                   estimate = round(estimate, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3)) %>%
  dplyr::arrange(drug, screen, readout, stratum) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- the breast stratum, declared underpowered before it was fitted ---------
message("\n   ", LINEAGE_FOCUS, " alone (n = ", n_focus,
        ") - DECLARED UNDERPOWERED. Reported, not read.")
focus_ids <- COV$ModelID[COV$lineage == LINEAGE_FOCUS]
focus_fits <- list(); j <- 0L
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  d0 <- ENDPOINTS %>%
    dplyr::filter(screen == ce$screen, drug == ce$drug,
                  readout == ce$readout, ModelID %in% focus_ids) %>%
    dplyr::inner_join(COV, by = "ModelID") %>%
    dplyr::transmute(Y, OX = .data[[RULER_PRIMARY]], MY = MY_ma, lineage)
  r <- .fit_one(d0, "OX + MY", "OX", "focus")
  if (!is.null(r)) {
    j <- j + 1L
    focus_fits[[j]] <- dplyr::bind_cols(
      tibble::tibble(rank = "focus stratum", model = "primary_breast",
                     screen = ce$screen, drug = ce$drug, readout = ce$readout,
                     ruler = RULER_PRIMARY, myc_estimator = MYC_PRIMARY,
                     stratum = LINEAGE_FOCUS, rhs = "OX + MY"), r)
  }
}
COEF_FOCUS <- dplyr::bind_rows(focus_fits)
if (nrow(COEF_FOCUS)) {
  COEF_FOCUS %>%
    dplyr::filter(drug %in% c(DRUG_PRIMARY, DRUG_CONTROL)) %>%
    dplyr::transmute(drug, screen, readout, n, estimate = round(estimate, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3)) %>%
    as.data.frame() %>% print(row.names = FALSE)
  COEF <- dplyr::bind_rows(COEF, COEF_FOCUS)
}

message("\n   THE POSITIVE CONTROL - RB1 loss on ", DRUG_PRIMARY,
        ". Varkaris et al. 2025 say this is the\n   strongest navitoclax ",
        "sensitiser in GDSC. NEGATIVE = RB1-null lines more sensitive.")
COEF %>%
  dplyr::filter(model == "rb1_control", ruler == RULER_PRIMARY,
                myc_estimator == MYC_PRIMARY) %>%
  dplyr::transmute(drug, screen, readout, n, estimate = round(estimate, 3),
                   ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                   p = signif(p, 2), ci_excludes_0) %>%
  dplyr::arrange(drug, screen, readout) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. THE READING, on the rules fixed in the declaration
# =============================================================================
message("\n8. the reading\n", strrep("-", 78))

.cell <- function(dr, rd, md = "primary", st = "all") {
  COEF %>% dplyr::filter(model == md, drug == dr, readout == rd,
                         ruler == RULER_PRIMARY, myc_estimator == MYC_PRIMARY,
                         stratum == st, screen %in% c("GDSC2", "PRISM"))
}
nav_ic50 <- .cell(DRUG_PRIMARY, "LN_IC50")
nav_auc  <- .cell(DRUG_PRIMARY, "AUC")
ven_ic50 <- .cell(DRUG_CONTROL, "LN_IC50")
ven_auc  <- .cell(DRUG_CONTROL, "AUC")

.neg <- function(x) isTRUE(nrow(x) == 1L && !is.na(x$estimate) &&
                            x$estimate < 0 && isTRUE(x$ci_excludes_0))
.null_ <- function(x) isTRUE(nrow(x) == 1L && !isTRUE(x$ci_excludes_0))

two_readouts   <- .neg(nav_ic50) && .neg(nav_auc)
control_is_null <- .null_(ven_ic50) && .null_(ven_auc)
prolif_survives <- {
  a <- .cell(DRUG_PRIMARY, "LN_IC50", "comp_prolif")
  b <- .cell(DRUG_PRIMARY, "AUC", "comp_prolif")
  .neg(a) && .neg(b)
}
t1 <- .cell(DRUG_PRIMARY, "LN_IC50", "secondary", "T1_low")
t3 <- .cell(DRUG_PRIMARY, "LN_IC50", "secondary", "T3_high")
myc_conditional <- isTRUE(nrow(t1) == 1L && nrow(t3) == 1L &&
                            !is.na(t1$estimate) && !is.na(t3$estimate) &&
                            t3$estimate < t1$estimate)

verdict <- if (two_readouts && control_is_null && prolif_survives &&
               myc_conditional) {
  paste0(
    "SUPPORTED, on the declaration's own definition. The primary OXPHOS\n",
    "   coefficient is negative for ", DRUG_PRIMARY, " on BOTH LN_IC50 and\n",
    "   AUC, ", DRUG_CONTROL, " is null, it survives the proliferation\n",
    "   companion, and it is larger in the MYC-high tertile.\n",
    "   E18 STILL APPLIES: the transcript configuration this would be read\n",
    "   off DOES NOT EXIST in these cells, so a positive here is not\n",
    "   explained by it and needs its own explanation.")
} else if (two_readouts && !control_is_null) {
  paste0(
    "NOT SHOWN - the specificity control moved with the primary.\n",
    "   ", DRUG_CONTROL, " is BCL-2 selective. If it tracks ",
    DRUG_PRIMARY, ", the\n   effect is NOT BCL-XL, whatever the primary ",
    "estimate does.")
} else if (two_readouts && !prolif_survives) {
  paste0("NOT SHOWN - the effect does not survive PROLIF_DISJOINT.\n",
         "   It is growth rate, or is not separable from it here.")
} else if (two_readouts && !myc_conditional) {
  paste0(
    "PARTIAL, AND THE WEAKER CLAIM. Both readouts carry a negative primary\n",
    "   coefficient, but the MYC strata are indistinguishable or run the\n",
    "   wrong way. That is a REAL BUT NON-MYC-CONDITIONAL effect, which is\n",
    "   NOT the claim the manuscript makes.")
} else {
  paste0(
    "NOT SHOWN - the primary is null on at least one of the two GDSC\n",
    "   readouts.\n",
    "   AND THIS IS THE UNINFORMATIVE BRANCH. E18 established that the\n",
    "   BCL2L1-OXPHOS configuration is ABSENT from cell lines, so a system\n",
    "   that does not carry the exposure cannot be asked to carry the\n",
    "   dependency. This null is a statement about cell lines, NOT about\n",
    "   the model. It is not evidence against the manuscript's claim and\n",
    "   must not be written as if it were.")
}
message("\n", verdict)
message("\n   NOTE ON 'TWO READOUTS': LN_IC50 and AUC are two summaries of ",
        "THE SAME\n   curve fit, not two platforms. They are not independent ",
        "evidence.\n   PRISM is the only genuinely independent platform here.")

# =============================================================================
# 9. The figure - one forest, per drug per readout, strata beside
# =============================================================================
message("\n9. figure")

# GDSC1 is kept as its OWN column rather than dropped. It is an independent
# screen of navitoclax and venetoclax on the same lines, and hiding it would
# make a two-column figure look like two platforms when LN_IC50 and AUC are
# two summaries of one curve.
FIG <- COEF %>%
  dplyr::filter(myc_estimator == MYC_PRIMARY, ruler == RULER_PRIMARY,
                model %in% c("primary", "secondary"), !is.na(estimate)) %>%
  dplyr::mutate(
    readout = paste(screen, sub("PRISM_", "", readout)),
    panel = dplyr::case_when(
      model == "primary"   ~ "all solid lines",
      stratum == "T1_low"  ~ "MYC tertile 1 (low)",
      stratum == "T3_high" ~ "MYC tertile 3 (high)"),
    panel = factor(panel, levels = c("all solid lines", "MYC tertile 1 (low)",
                                     "MYC tertile 3 (high)")),
    readout = factor(readout, levels = c("GDSC2 LN_IC50", "GDSC2 AUC",
                                         "GDSC1 LN_IC50", "GDSC1 AUC",
                                         "PRISM LFC")),
    drug = factor(drug, levels = rev(c(DRUG_PRIMARY, DRUG_CONTROL,
                                       setdiff(sort(unique(drug)),
                                               c(DRUG_PRIMARY, DRUG_CONTROL))))),
    is_key = drug %in% c(DRUG_PRIMARY, DRUG_CONTROL)) %>%
  dplyr::filter(!is.na(panel))

p <- ggplot2::ggplot(FIG, ggplot2::aes(x = estimate, y = drug,
                                       colour = is_key)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                         orientation = "y", width = 0, linewidth = 0.5) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::facet_grid(panel ~ readout, scales = "free_x") +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "#b2182b",
                                          `FALSE` = "grey45"),
                               guide = "none") +
  ggplot2::labs(
    title = paste0("OXPHOS coefficient on BH3-mimetic sensitivity, ",
                   length(COHORT), " solid cell lines"),
    subtitle = paste0("PRIMARY estimand: OX at fixed MYC (", MYC_PRIMARY,
                      "), ", RULER_PRIMARY, ", lineage as a factor.\n",
                      "NEGATIVE = more sensitive on all three readouts. ",
                      "Navitoclax and its venetoclax control in red.\n",
                      "EXPLORATORY AND POST-HOC. Not pre-registered."),
    x = "standardised coefficient on OXPHOS subunits (95% CI)", y = NULL,
    caption = paste0("GDSC ", GDSC_RELEASE, "; DepMap ", DEPMAP_RELEASE,
                     "; PRISM ", PRISM_RELEASE,
                     ". LN_IC50 and AUC summarise the SAME curve fit and are ",
                     "not independent platforms.")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(plot.title.position = "plot",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0),
                 strip.background = ggplot2::element_rect(fill = "grey93"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(PATH_E31_FIG, p, width = 13, height = 7.5, dpi = 200)
.ensure_dir(dirname(PATH_E31_FIG_DOC))
ggplot2::ggsave(PATH_E31_FIG_DOC, p, width = 13, height = 7.5, dpi = 200)
message("   ", PATH_E31_FIG)
message("   ", PATH_E31_FIG_DOC, "  (tracked - outputs/ is gitignored)")

# =============================================================================
# 10. Save
# =============================================================================
message("\n10. saving")

out <- list(
  coefficients = COEF,
  estimands    = ESTIMANDS,
  coverage     = coverage,
  drug_map     = dplyr::bind_rows(drug_map, prism_map),
  covariates   = COV,
  myc_estimator_cor = stats::cor(t(MYC_MAT), method = "spearman",
                                 use = "pairwise.complete.obs"),
  tertile_n    = as.data.frame(table(myc_tertile)),
  verdict      = verdict,
  lines = list(
    gdsc_lines       = length(gdsc_lines),
    join_rate        = join_rate,
    solid_joined     = n_solid_joined,
    with_expression  = n_with_expr,
    analysis_cohort  = length(COHORT),
    focus_stratum    = n_focus,
    n_lineages       = length(keep_l),
    dropped_lineages = dropped_lineages,
    rb1_loss         = sum(RB1)),
  spec = list(
    what = paste("Does a high respiratory state predict BH3-mimetic",
                 "sensitivity in cancer cell lines? PRIMARY estimand: the",
                 "OXPHOS coefficient at fixed MYC."),
    posture = paste("EXPLORATORY AND POST-HOC. Nothing pre-registered. The",
                    "estimand hierarchy was fixed in",
                    "docs/2026-09-28_bh3_mimetic_declaration.md and COMMITTED",
                    "BEFORE this script was written."),
    e18 = paste("E18 found the BCL2L1-OXPHOS configuration ABSENT from cell",
                "lines (0 of 4 rulers, breast and pan-cancer). A NULL HERE IS",
                "UNINFORMATIVE; a positive is interesting and would need its",
                "own explanation."),
    hierarchy_fixed = paste("The hierarchy is NOT reordered after the",
                            "estimates. A product term that fires while the",
                            "primary does not is reported as such, never",
                            "promoted."),
    signs = SIGN_NOTE,
    two_readouts = paste("LN_IC50 and AUC are two summaries of the SAME curve",
                         "fit, not two platforms. PRISM is the only",
                         "independent platform here."),
    rb1_rule = paste0("RB1 loss = cn_category '", RB1_CN_LOSS,
                      "' OR effect in {", paste(RB1_MUT_LOSS, collapse = ", "),
                      "}. Missense EXCLUDED; 'Loss' (one copy) EXCLUDED."),
    mediator = paste("BCL2L1 is a MEDIATOR, not a confounder. Fitted once,",
                     "reported separately, NEVER in the primary."),
    ruler_panel = paste("The 4-ruler panel is a sensitivity on the PRIMARY",
                        "model only, added at script-writing time before any",
                        "fit. Every other model uses", RULER_PRIMARY),
    scale = paste("Y z-scored within screen x drug x readout; OX and MY",
                  "z-scored over the cohort. A linear rescaling: every p and",
                  "every sign is unchanged. sd_Y_raw recovers raw units."),
    cohort_relative = paste("GSVA and every z-score are COHORT-RELATIVE. This",
                            "panel is a cohort. Its values are NEVER compared",
                            "numerically to TCGA, SCAN-B, mouse or CCLE-breast",
                            "values - only signs and orderings."),
    no_fdr = "No per-gene FDR, consistent with the rest of the project.",
    n3 = "drug-response statistics and transcripts; nothing is 'primed'",
    releases = c(gdsc = GDSC_RELEASE, depmap = DEPMAP_RELEASE,
                 prism = PRISM_RELEASE, expression = EXPR_FILE,
                 gsva = GSVA_VERSION),
    drugs_absent = DRUGS_ABSENT,
    seed = PROJECT_SEED),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E31)
utils::write.csv(as.data.frame(COEF), PATH_E31_CO, row.names = FALSE)
utils::write.csv(as.data.frame(coverage), PATH_E31_COV, row.names = FALSE)
message("   ", PATH_E31)
message("   ", PATH_E31_CO)
message("   ", PATH_E31_COV)
message("\nE31 done.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E31)

  cat(x$verdict, "\n")
  utils::str(x$lines)

  # The primary, every drug, both GDSC readouts and PRISM.
  x$coefficients %>%
    dplyr::filter(model == "primary", ruler == "ox_gsva",
                  myc_estimator == "myc_ma") %>%
    dplyr::transmute(drug, screen, readout, n, estimate = round(estimate, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                     p = signif(p, 2), ci_excludes_0) %>%
    dplyr::arrange(drug, screen, readout) %>%
    as.data.frame()

  # The four rulers on the primary estimand. CLAUDE.md trap 5: the instruments
  # disagree by a lot, and a result that lives on one of them is a ruler
  # result.
  x$coefficients %>%
    dplyr::filter(model == "primary", myc_estimator == "myc_ma",
                  drug %in% c("navitoclax", "venetoclax")) %>%
    dplyr::mutate(estimate = round(estimate, 3)) %>%
    dplyr::select(drug, screen, readout, ruler, estimate) %>%
    tidyr::pivot_wider(names_from = ruler, values_from = estimate) %>%
    as.data.frame()

  # The three MYC estimators on the primary. CLAUDE.md trap 4: MYC mRNA is not
  # MYC activity and the difference is total.
  x$coefficients %>%
    dplyr::filter(model == "primary", ruler == "ox_gsva",
                  drug %in% c("navitoclax", "venetoclax")) %>%
    dplyr::mutate(estimate = round(estimate, 3)) %>%
    dplyr::select(drug, screen, readout, myc_estimator, estimate) %>%
    tidyr::pivot_wider(names_from = myc_estimator, values_from = estimate) %>%
    as.data.frame()
  round(x$myc_estimator_cor, 3)

  # EVERY model of the hierarchy on the primary drug, including the failures.
  x$coefficients %>%
    dplyr::filter(drug == "navitoclax", ruler == "ox_gsva",
                  myc_estimator == "myc_ma", screen == "GDSC2") %>%
    dplyr::transmute(rank, model, stratum, readout, n,
                     estimate = round(estimate, 3), ci_lo = round(ci_lo, 3),
                     ci_hi = round(ci_hi, 3), ci_excludes_0) %>%
    dplyr::arrange(readout, model, stratum) %>%
    as.data.frame()

  # The product term. A NULL ON IT ALONE MEANS NOTHING - that estimand has
  # already failed in Block C, Block B, Block G, H4 and N2.
  x$coefficients %>%
    dplyr::filter(model == "tertiary", myc_estimator == "myc_ma") %>%
    dplyr::transmute(drug, screen, readout, estimate = round(estimate, 3),
                     ci_excludes_0) %>%
    as.data.frame()

  # Which drugs were matched, and under which spelling.
  x$drug_map %>% as.data.frame()

  # Set coverage, n_present / n_defined.
  x$coverage %>% as.data.frame()
}
