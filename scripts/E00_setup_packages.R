# E00_setup_packages.R
# =============================================================================
# Package loader, path constants and snapshot-shape expectations.
# Sourced by every numbered script. No side effects beyond attaching core
# packages.
#
# POLICY: this script CHECKS packages, it does not install them. Auto-installing
# behind a source() is how an analysis environment silently drifts between runs.
# A missing package produces one stop() naming everything that is absent.
#
# No renv; packages are installed system-wide. See CLAUDE.md.
# =============================================================================

# -----------------------------------------------------------------------------
# Package tiers
# -----------------------------------------------------------------------------
.pkg_core <- c("here", "dplyr", "tibble", "tidyr", "readr", "readxl",
               "stringr", "purrr", "ggplot2")

# Checked for presence, attached by the scripts that need them. GSVA is slow to
# attach and only E02 uses it.
.pkg_analysis <- c("GSVA", "decoupleR", "msigdbr")

.check_packages <- function(pkgs, tier) {
  have <- vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
  if (any(!have)) {
    stop("Missing ", tier, " packages: ", paste(pkgs[!have], collapse = ", "),
         "\nInstall them, then re-source this script. Nothing is auto-installed.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.check_packages(.pkg_core, "core")
.check_packages(.pkg_analysis, "analysis")

suppressPackageStartupMessages(
  invisible(lapply(.pkg_core, library, character.only = TRUE)))

# -----------------------------------------------------------------------------
# Paths - all relative to the project root via here::here(). No absolute paths.
# -----------------------------------------------------------------------------
DIR_DATA    <- here::here("data")
DIR_RESULTS <- here::here("results")
DIR_OUTPUTS <- here::here("outputs")
DIR_TABLES  <- here::here("outputs", "tables")
DIR_FIGURES <- here::here("outputs", "figures")

# The snapshot from myc_human_validation @ d3ac60e. Gitignored, ~563 MB.
# Provenance and the re-copy instruction: data/from_validation/README.md
DIR_SNAP <- here::here("data", "from_validation")

PATH_TCGA_VST    <- file.path(DIR_SNAP, "tcga_brca_vst.rds")
PATH_TCGA_LINEAR <- file.path(DIR_SNAP, "tcga_brca_linear.rds")
PATH_TCGA_MITO   <- file.path(DIR_SNAP, "tcga_brca_mito_scores.rds")
PATH_TCGA_MYC    <- file.path(DIR_SNAP, "tcga_brca_myc_scores.rds")
PATH_TCGA_COV    <- file.path(DIR_SNAP, "tcga_brca_covariates.rds")
PATH_SCANB_VST    <- file.path(DIR_SNAP, "scanb_vst.rds")
PATH_SCANB_LINEAR <- file.path(DIR_SNAP, "scanb_linear.rds")
PATH_SCANB_PHENO  <- file.path(DIR_SNAP, "scanb_pheno.rds")
PATH_G1           <- file.path(DIR_SNAP, "g1_overlap_audit.rds")

PATH_MITOCARTA <- here::here("data", "mitocarta_human", "Human.MitoCarta3.0.xls")
PATH_COLLECTRI <- here::here("data", "collectri_human",
                             "collectri_human_omnipath.tsv.gz")

# The death axis. See data/genesets_celldeath_human/README.md for why these are
# human-native and must never be remapped.
DIR_CELLDEATH  <- here::here("data", "genesets_celldeath_human")
PATH_CDC       <- file.path(DIR_CELLDEATH, "cell_death_genes_consolidated.csv")
DIR_TANG       <- file.path(DIR_CELLDEATH, "tang_modalities")

# The MYC estimator panel. data/genesets_myc_human/README.md.
DIR_MYCSETS    <- here::here("data", "genesets_myc_human")
PATH_MYC_GMX   <- file.path(DIR_MYCSETS, "myc_signature_genesets.gmx")
PATH_FELSHER   <- file.path(DIR_MYCSETS, "felsher_integrative_signature.csv")

.ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}
invisible(lapply(c(DIR_RESULTS, DIR_OUTPUTS, DIR_TABLES, DIR_FIGURES),
                 .ensure_dir))

PROJECT_SEED <- 20260831L

# =============================================================================
# Snapshot shape expectations
# =============================================================================
# Transcribed from the provenance READMEs. Scripts assert against these so a
# file swapped underneath the pipeline fails loudly instead of quietly changing
# a result. If one legitimately changes, re-snapshot and update BOTH the README
# and this block.
EXPECT_TCGA_GENES    <- 18115L
EXPECT_TCGA_SAMPLES  <- 1095L
EXPECT_SCANB_GENES   <- 18153L
EXPECT_SCANB_SAMPLES <- 3207L

EXPECT_CDC_ROWS      <- 1232L    # cell_death_genes_consolidated.csv
EXPECT_MYC_SETS      <- 16L      # signatures in the GMX
EXPECT_FELSHER_STRIP <- 61L      # g1$estimators_stripped$FELSHER

# The 15 Tang modalities, as DISTINCT GENES - not file rows. The CSVs carry one
# row per gene-per-evidence, so a gene with three PMIDs appears three times:
# Ferroptosis is 935 rows but 600 genes, Autophagy 1,195 rows but 876. Counting
# rows overstates every large set by up to 36%. A changed count means the
# upstream file moved and the snapshot must be re-taken, not patched.
EXPECT_TANG <- c(
  Alkaliptosis = 15L, Apoptosis = 608L,
  Autophagy_dependent_cell_death = 876L, Cuproptosis = 27L,
  Disulfidptosis = 16L, Entotic_cell_death = 17L, Ferroptosis = 600L,
  Immunogenic_cell_death = 34L, Lysosome_dependent_cell_death = 32L,
  MPT_driven_necrosis = 30L, NETotic_cell_death = 8L, Necroptosis = 83L,
  Oxeiptosis = 10L, Parthanatos = 10L, Pyroptosis = 51L)

# The three TCGA correlations that anchor the whole import. If E01 cannot
# reproduce these from the snapshot, nothing downstream is trustworthy.
# Source: myc_human_validation script 07 arm_summary, at d3ac60e.
EXPECT_ANCHOR <- c(
  `OXPHOS subunits`        =  0.388,
  `Mitochondrial ribosome` =  0.590,
  `Fatty acid oxidation`   = -0.140)
ANCHOR_TOL <- 0.001

# -----------------------------------------------------------------------------
# Shared helpers
# -----------------------------------------------------------------------------
.rho <- function(x, y) suppressWarnings(
  stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))

# Minimum genes present in the matrix for a set to be GSVA-scored at all.
MIN_SET_GENES <- 3L
# Floor below which a set is carried as a gene LIST rather than a score.
# CICD pro-survival (4 genes) sits under this deliberately - see CLAUDE.md
# trap 9. GSVA's own stability floor for a set is higher again; sets between
# these two are scored and flagged.
MIN_SCORE_N   <- 5L
MIN_GSVA_N    <- 15L

message("E00: setup complete (", length(.pkg_core), " core packages attached)")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  # Is the snapshot actually there? It is gitignored, so a fresh clone has an
  # empty data/from_validation/ and every script will stop in E01.
  file.exists(c(PATH_TCGA_VST, PATH_SCANB_VST, PATH_CDC, PATH_MYC_GMX))
  round(sum(file.size(list.files(DIR_SNAP, full.names = TRUE))) / 1e6, 0)

  # What the READMEs say the snapshot is pinned to.
  cat(readLines(file.path(DIR_SNAP, "README.md"), n = 20), sep = "\n")

}
