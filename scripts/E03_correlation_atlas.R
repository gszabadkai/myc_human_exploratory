# E03_correlation_atlas.R
# =============================================================================
# THE DATA LAYER. One tidy table, no plots, no verdicts.
#
# Built to: docs/2026-08-31_phase1_plan.md section 2 (E03)
#
#   cohort | instrument | myc_estimator | arm | stratum | n | rho | ci_lo |
#   ci_hi | p | adjusted
#
# 2 cohorts x 5 instruments x 21 MYC estimators x 31 measures x 8 strata x
# up to 3 adjustments. Roughly 70,000 cells.
#
# =============================================================================
# WHAT A CELL OF THIS TABLE IS WORTH: NOTHING
# =============================================================================
# CLAUDE.md, and it is the whole posture of the repo. Nothing here is
# pre-registered. Multiple comparisons are the default state, not an exception.
# A single cell of a 70,000-cell grid is uninteresting no matter what its p is.
#
# What the atlas is FOR is the three things a single cell cannot show:
#
#   1. GRADIENT   - does rho fall as the MYC signature's proliferation
#                   entanglement falls? (section 9; CLAUDE.md trap 3)
#   2. AGREEMENT  - do the four instruments say the same thing? (section 11;
#                   trap 5)
#   3. REPLICATION- does TCGA's pattern reappear in SCAN-B? (section 10)
#
# The p column exists so a cell can be DESCRIBED. It must never be used to
# SELECT one. E04 does not plot it.
#
# =============================================================================
# SCALE
# =============================================================================
# Every measure and estimator is rank-transformed by the engine before anything
# else happens, and ranks are invariant under monotone transforms. So mixing a
# log GSVA score with a linear mitoPPS score in one Spearman is well defined.
# The scale discipline in CLAUDE.md binds where scores are BUILT (E02).
#
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))
source(here::here("functions", "strata.R"))

message("\nE03: the correlation atlas\n", strrep("=", 78))

# =============================================================================
# 0. Constants
# =============================================================================
PATH_ATLAS     <- file.path(DIR_RESULTS, "correlation_atlas.rds")
PATH_ATLAS_CSV <- file.path(DIR_TABLES,  "correlation_atlas.csv")

# A stratum smaller than this is not emitted at all. TCGA Normal-like is 36 and
# survives; nothing else in either cohort is close to the floor.
MIN_STRATUM_N <- 30L

# The arms this study keeps coming back to. Used only to shorten printed
# summaries - the atlas itself carries all 18.
FOCUS_ARMS <- c("OXPHOS subunits", "mtDNA-encoded OXPHOS",
                "Mitochondrial ribosome", "TCA cycle", "Fatty acid oxidation")

# Descriptive, not inferential. Named here so it is a stated convention rather
# than a threshold discovered after looking.
REPLICATION_FLOOR <- 0.10

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

fr  <- readRDS(file.path(DIR_RESULTS, "frames.rds"))
mito <- readRDS(PATH_TCGA_MITO)                                    # TCGA arms
myc_t <- readRDS(PATH_TCGA_MYC)$estimators                         # M_a/M_b/M_c
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))        # E02
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))      # E02
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))     # E02

frames <- fr$frames
ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)
stopifnot(identical(ID_T, colnames(nw$tcga_gsva_new)),
          identical(ID_T, colnames(mito$mitopps_arms)),
          identical(ID_S, colnames(sc$mitopps_arms)),
          setequal(ID_T, frames$sample_id[frames$cohort == "TCGA"]),
          setequal(ID_S, frames$sample_id[frames$cohort == "SCAN-B"]))
message("   TCGA ", length(ID_T), " samples | SCAN-B ", length(ID_S))

# --- 1.1 the pin, verified rather than assumed -------------------------------
# E02 scored FELSHER__MITOSTRIP for TCGA in a SECOND GSVA call, months after
# snapshot's 18 arms were scored in the first. The .PIN_A/.PIN_B half-matrix
# sets are what is supposed to make those two calls walk the same gene universe.
# If they do, the new FELSHER__MITOSTRIP must BE the snapshot's M_a - not
# close to it. It is the same 61 genes, reconstructed by E02 rather than
# imported, so this also checks the variant labelling.
ma_snap <- myc_t$M_a[match(ID_T, myc_t$patient)]
f61_new <- nw$tcga_gsva_new[MYC_REF, ]
pin_diff <- max(abs(f61_new - ma_snap))
if (!is.finite(pin_diff) || pin_diff > 1e-10) {
  stop("the GSVA pin did not hold: ", MYC_REF, " rescored in E02 differs from ",
       "the snapshot's M_a by ", signif(pin_diff, 3), ". The two calls did not ",
       "walk the same gene universe, so TCGA's reused arms are NOT comparable ",
       "to the new sets. Do not build the atlas.", call. = FALSE)
}
message("   GSVA pin verified: ", MYC_REF,
        " rescored == snapshot M_a exactly (max abs diff ",
        signif(pin_diff, 3), ")")

# =============================================================================
# 2. The measures: 18 arms x 4 instruments, plus 13 mtDNA genes
# =============================================================================
# Rownames are "<instrument>::<measure>" and are decoded through meas_meta, not
# by splitting strings downstream.
message("\n2. measures")

ARMS <- rownames(mito$gsva_arms)
stopifnot(identical(ARMS, rownames(sc$gsva_arms)))
INSTRUMENTS <- c(gsva = "gsva_arms", mitopps = "mitopps_arms",
                 content = "content_arms", zmean = "zmean_arms")

.stack_arms <- function(obj) {
  m <- do.call(rbind, lapply(names(INSTRUMENTS), function(i) {
    x <- obj[[INSTRUMENTS[[i]]]][ARMS, , drop = FALSE]
    rownames(x) <- paste0(i, "::", ARMS)
    x
  }))
  m
}
MEAS_T <- .stack_arms(mito)
MEAS_S <- .stack_arms(sc)

# --- 2.1 the 13 mtDNA-encoded genes, individually ----------------------------
# CLAUDE.md trap 8 and plan sub-analysis (i). The arm is reported AND the genes
# are reported, because per gene they do not move together: in TCGA MT-ND6 sits
# at -0.154 and MT-CO2 at +0.243 against the same estimator. An arm score
# averages that away.
#
# Linear DESeq2-normalised expression. The engine ranks it, so linear vs log is
# immaterial here - see the header.
MT_GENES <- sd_$arm_sets[["mtDNA-encoded OXPHOS"]]
message("   ", length(MT_GENES), " mtDNA-encoded genes, carried individually")

tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(tcga_lin$scale,  "linear_deseq2_normalised"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))

mt_t <- MT_GENES
mt_s <- unname(ifelse(is.na(sc$symbol_map[MT_GENES]), MT_GENES,
                      sc$symbol_map[MT_GENES]))
if (!all(mt_t %in% rownames(tcga_lin$mat)) ||
    !all(mt_s %in% rownames(scanb_lin$mat))) {
  stop("mtDNA-encoded genes missing from a matrix: TCGA ",
       paste(setdiff(mt_t, rownames(tcga_lin$mat)), collapse = ","), " | SCAN-B ",
       paste(MT_GENES[!(mt_s %in% rownames(scanb_lin$mat))], collapse = ","),
       call. = FALSE)
}
g_t <- tcga_lin$mat[mt_t, ID_T, drop = FALSE]
g_s <- scanb_lin$mat[mt_s, ID_S, drop = FALSE]
rownames(g_t) <- rownames(g_s) <- paste0("gene::", MT_GENES)

MEAS_T <- rbind(MEAS_T, g_t)
MEAS_S <- rbind(MEAS_S, g_s)
rm(tcga_lin, scanb_lin, g_t, g_s); invisible(gc(verbose = FALSE))

meas_meta <- tibble::tibble(
  measure    = rownames(MEAS_T),
  instrument = sub("::.*$", "", rownames(MEAS_T)),
  arm        = sub("^[^:]+::", "", rownames(MEAS_T))) %>%
  dplyr::mutate(measure_class = ifelse(instrument == "gene", "mtdna_gene", "arm"))
stopifnot(identical(rownames(MEAS_T), rownames(MEAS_S)),
          !anyNA(MEAS_T), !anyNA(MEAS_S))
message("   ", nrow(MEAS_T), " measures (", length(ARMS), " arms x ",
        length(INSTRUMENTS), " instruments + ", length(MT_GENES), " genes)")

# =============================================================================
# 3. The estimators
# =============================================================================
# 18 MYC signatures x FOUR EXPLICITLY LABELLED VARIANTS (__FULL, __MITOSTRIP,
# __PROLIFSTRIP, __BOTHSTRIP - see the naming contract in E00), plus the four
# CollecTRI regulon variants scored by ULM, plus log2(MYC) (the mRNA, which
# CLAUDE.md trap 4 says is a different thing again). M_c_call - the copy-number
# call, -1/0/1/2 - is TCGA-only and is the only estimator carrying NA. The
# engine gives it its own complete-case set.
#
# EVERY ESTIMATOR CARRIES ITS strip_status INTO THE ATLAS. Phase 1's F1 was
# written on a panel that silently mixed one mito-stripped estimator with 17
# unstripped ones; a `strip_status` column is what stops that recurring.
message("\n3. estimators")

MYC_SIGS <- sd_$myc_panel$signature
MB_SIGS  <- rownames(nw$tcga_M_b_variants)
stopifnot(all(MYC_SIGS %in% rownames(nw$tcga_gsva_new)),
          all(MYC_SIGS %in% rownames(sc$gsva_new)),
          identical(MB_SIGS, rownames(sc$M_b_variants)),
          MYC_REF %in% MYC_SIGS, MB_REF %in% MB_SIGS)

EST_T <- rbind(nw$tcga_gsva_new[MYC_SIGS, ID_T, drop = FALSE],
               nw$tcga_M_b_variants[MB_SIGS, ID_T, drop = FALSE],
               log2MYC  = nw$tcga_log2MYC[ID_T],
               M_c_call = as.numeric(myc_t$M_c_call[match(ID_T, myc_t$patient)]))
EST_S <- rbind(sc$gsva_new[MYC_SIGS, ID_S, drop = FALSE],
               sc$M_b_variants[MB_SIGS, ID_S, drop = FALSE],
               log2MYC  = sc$log2MYC[ID_S])
colnames(EST_T) <- ID_T; colnames(EST_S) <- ID_S

est_meta <- dplyr::bind_rows(
  sd_$myc_panel %>%
    dplyr::select(myc_estimator = signature, base, strip_status, n_genes = n,
                  frac_prolif, frac_mito, thin),
  tibble::tibble(
    myc_estimator = MB_SIGS,
    base          = "M_b",
    strip_status  = sub("^.*__", "", MB_SIGS),
    n_genes       = vapply(sd_$collectri_sets[MB_SIGS], length, integer(1)),
    frac_prolif   = NA_real_, frac_mito = NA_real_, thin = FALSE),
  tibble::tibble(
    myc_estimator = c("log2MYC", "M_c_call"),
    base          = c("log2MYC", "M_c_call"),
    strip_status  = "NA",
    n_genes       = NA_integer_, frac_prolif = NA_real_, frac_mito = NA_real_,
    thin          = FALSE)) %>%
  dplyr::mutate(kind = dplyr::case_when(
    base == "M_b"      ~ "CollecTRI regulon (ULM)",
    base == "log2MYC"  ~ "MYC mRNA",
    base == "M_c_call" ~ "copy number (TCGA only)",
    TRUE               ~ "signature (GSVA)"))
message("   TCGA ", nrow(EST_T), " estimators | SCAN-B ", nrow(EST_S),
        " (M_c_call is TCGA-only)")
est_meta %>% dplyr::count(kind, strip_status) %>% as.data.frame() %>%
  print(row.names = FALSE)
message("   M_c_call missing in ", sum(is.na(EST_T["M_c_call", ])),
        " TCGA samples - it gets its own complete-case set, the other ",
        nrow(EST_T) - 1L, " keep their full n")

# =============================================================================
# 4. Covariates and adjustments
# =============================================================================
# raw          - nothing projected out
# prolif       - PROLIF_DISJOINT, the proliferation set built to be disjoint
#                from the MYC signatures (CLAUDE.md trap 3)
# purity_leuko - tumour purity and leukocyte fraction. TCGA ONLY. SCAN-B has no
#                purity estimate and it is never imputed (trap 2).
message("\n4. adjustments")

cov_t <- frames %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::select(sample_id, purity, leuko) %>%
  tibble::column_to_rownames("sample_id")
cov_t <- as.matrix(cov_t[ID_T, , drop = FALSE])

prolif_t <- matrix(mito$gsva_cov["PROLIF_DISJOINT", ID_T], ncol = 1,
                   dimnames = list(ID_T, "PROLIF_DISJOINT"))
prolif_s <- matrix(sc$gsva_cov["PROLIF_DISJOINT", ID_S], ncol = 1,
                   dimnames = list(ID_S, "PROLIF_DISJOINT"))

ADJ_T <- list(raw = NULL, prolif = prolif_t, purity_leuko = cov_t)
ADJ_S <- list(raw = NULL, prolif = prolif_s)
message("   TCGA: ", paste(names(ADJ_T), collapse = ", "))
message("   SCAN-B: ", paste(names(ADJ_S), collapse = ", "),
        " (no purity estimate exists - trap 2)")

# =============================================================================
# 5. Strata
# =============================================================================
message("\n5. strata")

# Strata come from functions/strata.R so every script cuts the cohorts
# identically. Luminal (LumA + LumB) was added 2026-09-01; it is additive and
# changes no stratum that existed before.
STR_T <- .build_strata(frames, "TCGA", ID_T)
STR_S <- .build_strata(frames, "SCAN-B", ID_S)
tibble::tibble(stratum = names(STR_T),
               TCGA    = vapply(STR_T, length, integer(1)),
               `SCAN-B` = vapply(STR_S, length, integer(1))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. The grid
# =============================================================================
message("\n6. building the atlas")

.run_cohort <- function(coh, EST, MEAS, strata, adjustments) {
  out <- list()
  for (st in names(strata)) {
    ids <- strata[[st]]
    if (length(ids) < MIN_STRATUM_N) {
      message("   ", coh, " / ", st, ": n = ", length(ids), " < ",
              MIN_STRATUM_N, ", not emitted")
      next
    }
    for (aj in names(adjustments)) {
      b <- .atlas_block(EST, MEAS, ids, adjustments[[aj]], min_n = MIN_STRATUM_N)
      if (is.null(b)) next
      out[[paste(st, aj)]] <- b %>%
        dplyr::mutate(cohort = coh, stratum = st, adjusted = aj)
    }
  }
  dplyr::bind_rows(out)
}

atlas <- dplyr::bind_rows(
  .run_cohort("TCGA",   EST_T, MEAS_T, STR_T, ADJ_T),
  .run_cohort("SCAN-B", EST_S, MEAS_S, STR_S, ADJ_S)) %>%
  dplyr::left_join(meas_meta, by = "measure") %>%
  dplyr::left_join(est_meta,  by = "myc_estimator") %>%
  dplyr::select(cohort, instrument, myc_estimator, base, strip_status, arm,
                measure_class, stratum, adjusted, n, k_cov, rho, ci_lo, ci_hi,
                p, frac_prolif, frac_mito, thin, kind, n_genes) %>%
  dplyr::arrange(cohort, instrument, myc_estimator, arm, stratum, adjusted)

message("   ", format(nrow(atlas), big.mark = ","), " cells")
atlas %>% dplyr::count(cohort, adjusted) %>% as.data.frame() %>%
  print(row.names = FALSE)

# =============================================================================
# 7. The anchor, recomputed FROM THE ATLAS
# =============================================================================
# E01 asserted these three against the snapshot. Asserting them again here, by
# reading them back out of the finished table, is what proves the atlas is
# INDEXED correctly - that the cell labelled (TCGA, gsva, FELSHER__MITOSTRIP,
# "OXPHOS subunits", all, raw) really is that correlation and not a neighbour's.
message("\n7. anchor, read back out of the atlas")

anchor <- atlas %>%
  dplyr::filter(cohort == "TCGA", instrument == "gsva",
                myc_estimator == MYC_REF, stratum == "all",
                adjusted == "raw", arm %in% names(EXPECT_ANCHOR)) %>%
  dplyr::transmute(arm, observed = rho,
                   expected = unname(EXPECT_ANCHOR[arm]),
                   diff = abs(observed - expected))
anchor %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)
if (nrow(anchor) != length(EXPECT_ANCHOR) || any(anchor$diff > ANCHOR_TOL)) {
  stop("the atlas does not reproduce the validation study's published ",
       "correlations. It is mis-indexed or the inputs moved. Nothing ",
       "downstream is trustworthy.", call. = FALSE)
}
message("   all ", nrow(anchor), " anchors reproduce within ", ANCHOR_TOL)

# =============================================================================
# 8. Sub-analysis (i): nuclear vs mtDNA-encoded OXPHOS
# =============================================================================
# Plan section 2, E03 (i). MYC does not regulate the two genomes alike, and on
# the composition instrument they point in OPPOSITE directions in TCGA. The arm
# summary and the per-gene table are both kept, because the 13 genes do not
# move together.
message("\n8. nuclear vs mtDNA-encoded")

NUC_MT <- c("OXPHOS subunits", "mtDNA-encoded OXPHOS", "CI subunits",
            "CII subunits", "CIII subunits", "CIV subunits", "CV subunits")
nuclear_vs_mtdna <- atlas %>%
  dplyr::filter(measure_class == "arm", arm %in% NUC_MT, stratum == "all") %>%
  dplyr::select(cohort, instrument, adjusted, myc_estimator, arm, n, rho,
                ci_lo, ci_hi)
nuclear_vs_mtdna %>%
  dplyr::filter(myc_estimator == MYC_REF, adjusted == "raw") %>%
  tidyr::pivot_wider(id_cols = c(cohort, arm), names_from = instrument,
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

mtdna_genes <- atlas %>%
  dplyr::filter(measure_class == "mtdna_gene", stratum == "all") %>%
  dplyr::select(cohort, adjusted, myc_estimator, gene = arm, n, rho, ci_lo,
                ci_hi, p)
message("\n   the 13 genes against ", MYC_REF, ", raw, both cohorts:")
mtdna_genes %>%
  dplyr::filter(myc_estimator == MYC_REF, adjusted == "raw") %>%
  tidyr::pivot_wider(id_cols = gene, names_from = cohort, values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(TCGA)) %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 9. Sub-analysis (ii): the estimator panel, ordered by entanglement
# =============================================================================
# Plan section 2, E03 (ii); CLAUDE.md trap 3. THE SHARPEST TEST IN PHASE 1.
#
# If rho(MYC, OXPHOS) decays monotonically as the signature's proliferation
# content falls, the correlation is proliferation wearing MYC's name. If it
# survives MYC_UP.V1_UP (1.5% entangled) and BILD (3.9%), it is not.
#
# `slope_vs_entanglement` is the Spearman correlation, ACROSS the 18 signatures,
# between each signature's frac_prolif and the rho it produces. It is a single
# number per (cohort, instrument, arm, adjustment) and it is the summary that
# actually answers the question - not any one signature's rho.
message("\n9. the estimator panel by proliferation entanglement")

panel <- atlas %>%
  dplyr::filter(kind == "signature (GSVA)", measure_class == "arm",
                stratum == "all") %>%
  dplyr::select(cohort, instrument, adjusted, arm, myc_estimator, base,
                strip_status, frac_prolif, frac_mito, thin, n, rho, ci_lo, ci_hi)

# GROUPED BY strip_status, because a slope computed across a mixture of FULL and
# stripped variants would be comparing sets of different composition as if they
# were the same estimator - which is exactly what F1 did by accident.
entanglement_slope <- panel %>%
  dplyr::group_by(cohort, instrument, adjusted, arm, strip_status) %>%
  dplyr::summarise(
    n_sig       = dplyr::n(),
    rho_min     = min(rho), rho_max = max(rho),
    rho_least   = rho[which.min(frac_prolif)],   # MYC_UP.V1_UP, 1.5%
    rho_most    = rho[which.max(frac_prolif)],   # YU_MYC_TARGETS_UP, 47.6%
    slope_vs_entanglement = stats::cor(rho, frac_prolif, method = "spearman"),
    .groups = "drop")

message("\n   OXPHOS subunits, raw - does rho track entanglement, per variant?")
entanglement_slope %>%
  dplyr::filter(arm == "OXPHOS subunits", adjusted == "raw",
                instrument == "gsva") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   OXPHOS subunits, GSVA, raw - every signature x every variant:")
panel %>%
  dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                adjusted == "raw") %>%
  tidyr::pivot_wider(id_cols = c(base, frac_prolif), names_from =
                       c(cohort, strip_status), values_from = rho) %>%
  dplyr::arrange(frac_prolif) %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   frac_prolif is the FULL variant's; a PROLIFSTRIP column is by",
        " construction 0%.")

# =============================================================================
# 10. Cross-cohort reproducibility
# =============================================================================
# CLAUDE.md: "'Significant' is not a result here. Consistency across cohorts,
# instruments and estimators is." This is that table.
#
# `replicated` is DESCRIPTIVE: same sign and |rho| >= REPLICATION_FLOOR in both
# cohorts. It is a stated convention, declared above before it was applied, not
# a threshold chosen after seeing the numbers. It is not a test and carries no
# error rate.
message("\n10. cross-cohort reproducibility")

reproducibility <- atlas %>%
  dplyr::filter(adjusted %in% c("raw", "prolif")) %>%
  dplyr::select(cohort, instrument, myc_estimator, arm, measure_class, stratum,
                adjusted, n, rho, ci_lo, ci_hi) %>%
  tidyr::pivot_wider(names_from = cohort,
                     values_from = c(n, rho, ci_lo, ci_hi)) %>%
  dplyr::rename(rho_tcga = `rho_TCGA`, rho_scanb = `rho_SCAN-B`,
                n_tcga = `n_TCGA`, n_scanb = `n_SCAN-B`) %>%
  dplyr::filter(!is.na(rho_tcga), !is.na(rho_scanb)) %>%
  dplyr::mutate(
    diff       = rho_scanb - rho_tcga,
    same_sign  = sign(rho_tcga) == sign(rho_scanb),
    ci_overlap = pmax(`ci_lo_TCGA`, `ci_lo_SCAN-B`) <=
                 pmin(`ci_hi_TCGA`, `ci_hi_SCAN-B`),
    replicated = same_sign & abs(rho_tcga) >= REPLICATION_FLOOR &
                 abs(rho_scanb) >= REPLICATION_FLOOR)

message("   ", format(nrow(reproducibility), big.mark = ","),
        " comparable cells; same sign in ",
        round(100 * mean(reproducibility$same_sign), 1), "%, CIs overlap in ",
        round(100 * mean(reproducibility$ci_overlap), 1), "%")

message("\n   agreement of the 18 arms between cohorts (", MYC_REF, ", all, raw):")
reproducibility %>%
  dplyr::filter(myc_estimator == MYC_REF, stratum == "all",
                adjusted == "raw", measure_class == "arm") %>%
  dplyr::group_by(instrument) %>%
  dplyr::summarise(n_arms = dplyr::n(),
                   rho_across_arms = round(stats::cor(rho_tcga, rho_scanb,
                                                      method = "spearman"), 3),
                   max_abs_diff = round(max(abs(diff)), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 11. Instrument agreement
# =============================================================================
# CLAUDE.md trap 5: the four instruments disagree, by a lot, and instrument
# choice is not cosmetic. This measures the disagreement AT THE LEVEL OF THE
# ANSWER - do two instruments rank the 18 arms the same way against the same
# MYC estimator? - rather than at the level of the raw scores.
message("\n11. instrument agreement (trap 5)")

.pairs <- utils::combn(names(INSTRUMENTS), 2)
instrument_agreement <- atlas %>%
  dplyr::filter(measure_class == "arm", stratum == "all", adjusted == "raw") %>%
  dplyr::select(cohort, myc_estimator, instrument, arm, rho) %>%
  tidyr::pivot_wider(names_from = instrument, values_from = rho) %>%
  dplyr::group_by(cohort, myc_estimator) %>%
  dplyr::group_modify(~ tibble::tibble(
    instrument_a = .pairs[1, ], instrument_b = .pairs[2, ],
    rho_across_arms = vapply(seq_len(ncol(.pairs)), function(j)
      stats::cor(.x[[.pairs[1, j]]], .x[[.pairs[2, j]]], method = "spearman"),
      numeric(1)))) %>%
  dplyr::ungroup()

instrument_agreement %>%
  dplyr::group_by(instrument_a, instrument_b) %>%
  dplyr::summarise(median_rho = round(stats::median(rho_across_arms), 3),
                   min_rho = round(min(rho_across_arms), 3),
                   max_rho = round(max(rho_across_arms), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 12. Coverage carried forward
# =============================================================================
# E02 warned that TANG_FERROPTOSIS (0.70) and TANG_IMMUNOGENIC_CELL_DEATH (0.62)
# sit below the 0.80 coverage floor IN BOTH COHORTS IDENTICALLY - which means
# the missing symbols are absent from the gene sets' own vocabulary, not from a
# platform. E05 consumes those sets; the flag travels with them.
low_coverage <- sd_$coverage %>% dplyr::filter(frac < 0.80) %>%
  dplyr::arrange(set, cohort)
message("\n12. sets below 0.80 coverage, carried forward to E05: ",
        paste(unique(low_coverage$set), collapse = ", "))

# =============================================================================
# 13. Save
# =============================================================================
message("\n13. save")

saveRDS(list(
  atlas               = atlas,
  reproducibility     = reproducibility,
  panel               = panel,
  entanglement_slope  = entanglement_slope,
  nuclear_vs_mtdna    = nuclear_vs_mtdna,
  mtdna_genes         = mtdna_genes,
  instrument_agreement = instrument_agreement,
  anchor              = anchor,
  low_coverage        = low_coverage,
  est_meta            = est_meta,
  meas_meta           = meas_meta,
  strata              = list(TCGA = lengths(STR_T), `SCAN-B` = lengths(STR_S)),
  settings = list(min_stratum_n = MIN_STRATUM_N,
                  replication_floor = REPLICATION_FLOOR,
                  ci = "Fisher-z, Bonett-Wright variance for Spearman",
                  adjustments = list(TCGA = names(ADJ_T), `SCAN-B` = names(ADJ_S))),
  rules = list(
    exploratory = paste("nothing here is pre-registered; every cell is",
                        "hypothesis-generating"),
    p = paste("the p column describes a cell, it never selects one. The grid",
              "is ~70,000 cells and no multiplicity correction would make a",
              "single one of them interesting."),
    scanb_purity = "SCAN-B has no purity estimate; purity_leuko is TCGA-only",
    mitopps = paste("mitoPPS is composition, not level. Never compare its",
                    "VALUES across cohorts - only patterns."),
    mtdna = paste("the mtDNA arm and the 13 genes are both reported; the genes",
                  "do not move together")),
  built = Sys.time()), PATH_ATLAS)

readr::write_csv(atlas, PATH_ATLAS_CSV)

message("\nE03: done.")
message("    results/correlation_atlas.rds        the atlas and its summaries")
message("    outputs/tables/correlation_atlas.csv ", format(nrow(atlas),
        big.mark = ","), " rows")
message("    NEXT: E04, the figures.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  a <- readRDS(PATH_ATLAS)

  # --- (ii) the question phase 1 exists to answer --------------------------
  # Does the MYC-OXPHOS correlation survive the least proliferation-entangled
  # signatures? Read DOWN the table: frac_prolif rises, rho should not track it.
  a$panel %>%
    dplyr::filter(arm == "OXPHOS subunits", adjusted == "raw") %>%
    tidyr::pivot_wider(id_cols = c(myc_estimator, frac_prolif),
                       names_from = c(cohort, instrument), values_from = rho) %>%
    dplyr::arrange(frac_prolif) %>% as.data.frame()

  # and the one-number summary of the same thing, over every arm
  a$entanglement_slope %>%
    dplyr::filter(adjusted == "raw", instrument == "gsva") %>%
    dplyr::arrange(slope_vs_entanglement) %>% as.data.frame()

  # --- what the proliferation adjustment costs -----------------------------
  a$atlas %>%
    dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                  stratum == "all", myc_estimator %in%
                    c(MYC_LOW_ENTANG, MYC_REF, MYC_HALLMARK,
                      "YU_MYC_TARGETS_UP__FULL", MB_REF, "log2MYC")) %>%
    tidyr::pivot_wider(id_cols = c(cohort, myc_estimator),
                       names_from = adjusted, values_from = rho) %>%
    as.data.frame()

  # --- (i) the two genomes -------------------------------------------------
  a$nuclear_vs_mtdna %>%
    dplyr::filter(myc_estimator == MYC_REF, adjusted == "raw") %>%
    tidyr::pivot_wider(id_cols = c(arm, cohort), names_from = instrument,
                       values_from = rho) %>% as.data.frame()

  a$mtdna_genes %>% dplyr::filter(adjusted == "raw",
                                  myc_estimator %in% c(MYC_REF, MB_REF)) %>%
    tidyr::pivot_wider(id_cols = gene, names_from = c(cohort, myc_estimator),
                       values_from = rho) %>% as.data.frame()

  # --- replication ---------------------------------------------------------
  a$reproducibility %>%
    dplyr::filter(stratum == "all", adjusted == "raw", measure_class == "arm") %>%
    dplyr::count(instrument, replicated) %>% as.data.frame()

  # the cells that do NOT replicate are the interesting ones
  a$reproducibility %>%
    dplyr::filter(stratum == "all", adjusted == "raw", !same_sign,
                  abs(rho_tcga) > 0.15 | abs(rho_scanb) > 0.15) %>%
    dplyr::arrange(dplyr::desc(abs(diff))) %>% head(20) %>% as.data.frame()

  # --- strata --------------------------------------------------------------
  a$atlas %>%
    dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                  myc_estimator == MYC_REF, adjusted == "raw") %>%
    dplyr::select(cohort, stratum, n, rho, ci_lo, ci_hi) %>% as.data.frame()

  # --- trap 5 --------------------------------------------------------------
  a$instrument_agreement %>%
    dplyr::filter(myc_estimator == MYC_REF) %>% as.data.frame()

}
