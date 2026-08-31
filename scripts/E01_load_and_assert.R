# E01_load_and_assert.R
# =============================================================================
# Load the snapshot, prove it is faithful, and build the one harmonised
# covariate frame both cohorts share. NOTHING IS SCORED HERE.
#
# Built to: docs/2026-08-31_phase1_plan.md sections 1.1, 2, 4
#
# =============================================================================
# WHY THIS SCRIPT EXISTS AT ALL
# =============================================================================
# data/from_validation/ is a COPY. Copies rot: a re-run of the validation
# pipeline, a half-finished cp, a file swapped by hand. So section 4 of the plan
# fixes an end-to-end anchor - three correlations whose values are already
# reported in the validation study - and this script asserts them.
#
#   rho(M_a, `OXPHOS subunits` GSVA)  =  0.388
#   rho(M_a, `Mitochondrial ribosome`)  =  0.590
#   rho(M_a, `Fatty acid oxidation`)  = -0.140
#
# If those do not come back, the import is not faithful and nothing downstream
# is trustworthy. That is a stop, not a warning.
#
# =============================================================================
# SCALE - stated here once, enforced everywhere
# =============================================================================
#   *_vst.rds     LOG (variance-stabilised)   -> GSVA only, kcdf = "Gaussian"
#   *_linear.rds  LINEAR (DESeq2-normalised)  -> mitoPPS and log2(gene) only
# These are OPPOSITE requirements and the two objects must never be confused.
# This script asserts both, in both cohorts.
#
# SPECIES: human. No ortholog function is called anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE01: load the snapshot and prove it is faithful\n", strrep("=", 78))

# =============================================================================
# 1. Is the snapshot present at all?
# =============================================================================
# data/from_validation/ is gitignored, so a fresh clone has an empty directory
# and every downstream script would fail somewhere less informative than here.
message("\n1. snapshot presence")

.needed <- c(PATH_TCGA_VST, PATH_TCGA_LINEAR, PATH_TCGA_MITO, PATH_TCGA_MYC,
             PATH_TCGA_COV, PATH_SCANB_VST, PATH_SCANB_LINEAR, PATH_SCANB_PHENO,
             PATH_G1)
missing <- .needed[!file.exists(.needed)]
if (length(missing)) {
  stop("the validation snapshot is missing ", length(missing), " file(s):\n  ",
       paste(basename(missing), collapse = "\n  "),
       "\n\ndata/from_validation/ is gitignored and is NOT on origin. Re-copy it",
       " from\n  /Users/gs/code/myc_human_validation/results/\n",
       "at the commit recorded in data/from_validation/README.md, then re-run.",
       call. = FALSE)
}
message("   all ", length(.needed), " snapshot files present (",
        round(sum(file.size(.needed)) / 1e6), " MB)")

tcga_vst   <- readRDS(PATH_TCGA_VST)
tcga_lin   <- readRDS(PATH_TCGA_LINEAR)
tcga_mito  <- readRDS(PATH_TCGA_MITO)
tcga_myc   <- readRDS(PATH_TCGA_MYC)$estimators
tcga_cov   <- readRDS(PATH_TCGA_COV)$covariates
scanb_vst  <- readRDS(PATH_SCANB_VST)
scanb_lin  <- readRDS(PATH_SCANB_LINEAR)
scanb_ph   <- readRDS(PATH_SCANB_PHENO)
g1         <- readRDS(PATH_G1)

# =============================================================================
# 2. Scale and shape
# =============================================================================
message("\n2. scale and shape")

.assert_scale <- function(obj, want, label) {
  if (!identical(obj$scale, want)) {
    stop(label, ": expected scale '", want, "', got '", obj$scale,
         "'. The wrong object has been snapshotted.", call. = FALSE)
  }
  invisible(TRUE)
}
.assert_scale(tcga_vst,  "log_vst",                   "TCGA VST")
.assert_scale(tcga_lin,  "linear_deseq2_normalised",  "TCGA linear")
.assert_scale(scanb_vst, "log_vst",                   "SCAN-B VST")
.assert_scale(scanb_lin, "linear_deseq2_normalised",  "SCAN-B linear")

.assert_dim <- function(m, g, s, label) {
  if (nrow(m) != g || ncol(m) != s) {
    stop(label, ": ", nrow(m), " x ", ncol(m), ", expected ", g, " x ", s,
         ". The snapshot does not match E00's expectations - re-snapshot and ",
         "update BOTH the README and E00.", call. = FALSE)
  }
}
.assert_dim(tcga_vst$mat,  EXPECT_TCGA_GENES,  EXPECT_TCGA_SAMPLES,  "TCGA VST")
.assert_dim(tcga_lin$mat,  EXPECT_TCGA_GENES,  EXPECT_TCGA_SAMPLES,  "TCGA linear")
.assert_dim(scanb_vst$mat, EXPECT_SCANB_GENES, EXPECT_SCANB_SAMPLES, "SCAN-B VST")
.assert_dim(scanb_lin$mat, EXPECT_SCANB_GENES, EXPECT_SCANB_SAMPLES, "SCAN-B linear")

# The two matrices of a cohort must describe the SAME genes and samples in the
# SAME order, or a per-gene comparison across instruments silently misaligns.
stopifnot(identical(dimnames(tcga_vst$mat),  dimnames(tcga_lin$mat)),
          identical(dimnames(scanb_vst$mat), dimnames(scanb_lin$mat)))

# The scale assertion proper. A VST matrix maxes in the tens; a linear one in
# the millions. If these ever cross, something has collapsed them.
for (x in list(list(tcga_vst, tcga_lin, "TCGA"), list(scanb_vst, scanb_lin, "SCAN-B"))) {
  mv <- max(x[[1]]$mat); ml <- max(x[[2]]$mat)
  if (mv > 100 || ml < 100) {
    stop(x[[3]], " scale assertion failed: VST max ", round(mv, 1),
         ", linear max ", round(ml, 1), ".", call. = FALSE)
  }
  message(sprintf("   %-7s %5d genes x %4d samples | VST max %.1f | linear max %s",
                  x[[3]], nrow(x[[1]]$mat), ncol(x[[1]]$mat), mv,
                  format(round(ml), big.mark = ",")))
}

# The instruments TCGA already carries, and the set definitions that travel
# with them. E02 reuses arm_sets rather than rebuilding it.
stopifnot(all(c("gsva_arms", "mitopps_arms", "content_arms", "zmean_arms",
                "arm_sets", "covariate_sets", "mito_paths",
                "arm_universe_path") %in% names(tcga_mito)))
message("   TCGA carries ", nrow(tcga_mito$gsva_arms), " arms on 4 instruments")
message("   arms: ", paste(utils::head(rownames(tcga_mito$gsva_arms), 4),
                           collapse = ", "), ", ...")

# =============================================================================
# 3. THE ANCHOR - is this snapshot the one the validation study reported on?
# =============================================================================
# Three correlations, recomputed here from the copied objects and compared with
# the values in myc_human_validation script 07's arm_summary at d3ac60e.
message("\n3. anchor: reproducing three published TCGA correlations")

pat <- colnames(tcga_mito$gsva_arms)
stopifnot(identical(sort(pat), sort(tcga_myc$patient)))
M_a <- tcga_myc$M_a[match(pat, tcga_myc$patient)]

anchor <- tibble::tibble(
  arm      = names(EXPECT_ANCHOR),
  expected = unname(EXPECT_ANCHOR),
  observed = vapply(names(EXPECT_ANCHOR),
                    function(a) .rho(tcga_mito$gsva_arms[a, ], M_a), numeric(1))) %>%
  dplyr::mutate(diff = abs(observed - expected))
anchor %>%
  dplyr::mutate(observed = round(observed, 4), diff = signif(diff, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (any(anchor$diff > ANCHOR_TOL)) {
  stop("the snapshot does NOT reproduce the validation study's reported ",
       "correlations (max |diff| = ", signif(max(anchor$diff), 3),
       ", tolerance ", ANCHOR_TOL, ").\n",
       "data/from_validation/ is not what its README says it is. Re-copy from ",
       "the pinned commit before going further - every number downstream would ",
       "otherwise be computed on unknown inputs.", call. = FALSE)
}
message("   all three reproduce within ", ANCHOR_TOL,
        ". The snapshot is faithful.")

# The two facts from the plan that motivate the whole study, printed so they
# are visible before any new analysis exists.
message(sprintf("\n   for orientation: rho(M_a, OXPHOS subunits) = %+.3f",
                anchor$observed[anchor$arm == "OXPHOS subunits"]))
message(sprintf("                    rho(M_a, mtDNA-encoded)   = %+.3f",
                .rho(tcga_mito$gsva_arms["mtDNA-encoded OXPHOS", ], M_a)))
message(sprintf("                    rho(log2(MYC), OXPHOS)    = %+.3f",
                .rho(tcga_mito$gsva_arms["OXPHOS subunits", ],
                     log2(tcga_lin$mat["MYC", pat]))))
message("   MYC ACTIVITY and MYC mRNA are different variables (CLAUDE.md trap 4),")
message("   and MYC does not track the mtDNA-encoded subunits (trap 8).")

# =============================================================================
# 4. The harmonised covariate frame
# =============================================================================
# One frame, both cohorts, only the columns that exist in both plus the
# TCGA-only confounders. Purity and leukocyte fraction are TCGA-only and stay
# NA in SCAN-B - NEVER imputed (CLAUDE.md trap 2).
message("\n4. harmonised covariates")

# PAM50 labels differ between cohorts ("BRCA_Her2" vs "HER2"). Harmonised to one
# vocabulary; the mapping is explicit rather than a regex, so a new level fails
# loudly instead of being silently dropped.
.PAM50_LEVELS <- c("LumA", "LumB", "HER2", "Basal", "Normal")
.harmonise_pam50 <- function(x, label) {
  y <- sub("^BRCA_", "", as.character(x))
  y <- ifelse(y == "Her2", "HER2", y)
  y[y %in% c("NA", "", "Unclassified")] <- NA_character_
  bad <- setdiff(unique(y[!is.na(y)]), .PAM50_LEVELS)
  if (length(bad)) {
    stop(label, ": unrecognised PAM50 level(s) -> ", paste(bad, collapse = ", "),
         ". Add them to .PAM50_LEVELS deliberately.", call. = FALSE)
  }
  factor(y, levels = .PAM50_LEVELS)
}

frame_tcga <- tibble::tibble(
  cohort    = "TCGA",
  sample_id = pat,
  PAM50     = .harmonise_pam50(tcga_cov$PAM50[match(pat, tcga_cov$patient)], "TCGA"),
  ER        = factor(dplyr::recode(tcga_cov$er_call[match(pat, tcga_cov$patient)],
                                   Positive = "ERpos", Negative = "ERneg",
                                   .default = NA_character_),
                     levels = c("ERpos", "ERneg")),
  purity    = tcga_cov$purity[match(pat, tcga_cov$patient)],
  leuko     = tcga_cov$leukocyte_fraction[match(pat, tcga_cov$patient)],
  age       = NA_real_,
  grade     = NA_character_)

sp <- scanb_ph$pheno
frame_scanb <- tibble::tibble(
  cohort    = "SCAN-B",
  sample_id = sp$sample_id,
  PAM50     = .harmonise_pam50(sp$PAM50, "SCAN-B"),
  ER        = factor(dplyr::recode(as.character(sp$er_call),
                                   `1` = "ERpos", `0` = "ERneg",
                                   .default = NA_character_),
                     levels = c("ERpos", "ERneg")),
  purity    = NA_real_,          # not available - see CLAUDE.md trap 2
  leuko     = NA_real_,
  age       = sp$age,
  grade     = as.character(sp$NHG))

frames <- dplyr::bind_rows(frame_tcga, frame_scanb) %>%
  dplyr::mutate(cohort = factor(cohort, levels = c("TCGA", "SCAN-B")))

stopifnot(identical(frame_tcga$sample_id,  colnames(tcga_vst$mat)),
          identical(frame_scanb$sample_id, colnames(scanb_vst$mat)))

# --- 4.1 stratum sizes, printed because they decide what is answerable -------
# SCAN-B is a lower FRACTION ER-negative than TCGA but has roughly DOUBLE the
# absolute number. For stratified work it is the stronger cohort, not the
# weaker one - the opposite of the intuition.
message("\n   PAM50:")
frames %>% dplyr::count(cohort, PAM50) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = n, values_fill = 0L) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   ER:")
frames %>% dplyr::count(cohort, ER) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = n, values_fill = 0L) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   NOTE: SCAN-B has the larger ABSOLUTE ER-negative and Basal strata,")
message("   despite the smaller fraction. Stratified reads are stronger there.")

# --- 4.2 the confounders, in the cohort that has them -----------------------
message(sprintf(
  "\n   TCGA purity: median %.2f (n = %d); leukocyte fraction: median %.3f",
  stats::median(frames$purity, na.rm = TRUE), sum(!is.na(frames$purity)),
  stats::median(frames$leuko, na.rm = TRUE)))
message("   SCAN-B has NEITHER. Every correlation is reported raw in both ",
        "cohorts and\n   additionally purity-adjusted in TCGA. Never imputed.")

# =============================================================================
# 5. Save
# =============================================================================
message("\n5. save")

out <- list(
  frames = frames,
  anchor = anchor,
  cohorts = list(
    TCGA   = list(n = EXPECT_TCGA_SAMPLES,  genes = EXPECT_TCGA_GENES),
    `SCAN-B` = list(n = EXPECT_SCANB_SAMPLES, genes = EXPECT_SCANB_GENES)),
  provenance = list(
    snapshot   = "data/from_validation/, myc_human_validation @ d3ac60e",
    anchor_src = "myc_human_validation script 07 arm_summary",
    pam50      = "harmonised to LumA/LumB/HER2/Basal/Normal; Unclassified -> NA",
    purity     = "TCGA only; NA in SCAN-B and never imputed",
    exploratory = paste("nothing in this repo is pre-registered; every finding",
                        "is hypothesis-generating")),
  built = Sys.time())

saveRDS(out, file.path(DIR_RESULTS, "frames.rds"))
readr::write_csv(frames, file.path(DIR_TABLES, "covariate_frame.csv"))

message("\nE01: done.")
message("    results/frames.rds")
message("    NEXT: E02, which scores SCAN-B and adds the death sets to both.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  f <- readRDS(file.path(DIR_RESULTS, "frames.rds"))

  # --- the anchor, which is the thing to read first ------------------------
  f$anchor %>% as.data.frame() %>% print(row.names = FALSE)

  # --- what is actually answerable, per stratum ----------------------------
  f$frames %>% dplyr::count(cohort, PAM50, ER) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- the full TCGA arm summary, already computed by the validation study --
  # This is the map of where MYC sits relative to every mitochondrial
  # programme, and it is the starting point of the whole study. Note
  # Mitochondrial ribosome ABOVE OXPHOS subunits, and mtDNA near zero.
  m <- readRDS(PATH_TCGA_MITO)
  m$arm_summary %>%
    dplyr::select(arm, n_genes, rho_M_a, rho_instruments, rho_purity, rho_leuko) %>%
    dplyr::arrange(dplyr::desc(rho_M_a)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # --- is purity going to be a problem? ------------------------------------
  # rho_purity ~0.2-0.3 on the mito arms. Not fatal, not ignorable, and
  # UNMEASURABLE in SCAN-B.
  plot(m$arm_summary$rho_purity, m$arm_summary$rho_M_a, pch = 16,
       xlab = "rho(arm, purity)", ylab = "rho(arm, MYC activity)")
  text(m$arm_summary$rho_purity, m$arm_summary$rho_M_a,
       m$arm_summary$arm, pos = 4, cex = 0.6)

}
