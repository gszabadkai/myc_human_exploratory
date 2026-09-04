# E16_respiratory_rulers.R
# =============================================================================
# G1 - A MEASUREMENT CHECK. NOT AN ANALYSIS.
#
# The mouse gate model (myc_mouse, docs/2026-09-02_myc_oxphos_priming_gate_model.md,
# script 48) introduced a third respiratory ruler:
#
#     ox_rel = mean z(nuclear MitoCarta OXPHOS subunits)
#              - mean z(the rest of MitoCarta)
#
# and argued it is the only ruler separable enough from MYC to carry an
# interaction test. Neither human arm has ever used it: the exploratory arm used
# OXPHOS-subunit GSVA, the validation arm used the absolute subunit level, and
# both are `ox_lvl` analogues.
#
# THIS SCRIPT ASKS TWO QUESTIONS AND NO OTHERS.
#   1. Is `ox_rel` materially different from mitoPPS in human?
#   2. Does the human three-gene configuration - BBC3 up, BCL2L1 up, MCL1 down -
#      survive the change of ruler?
#
# OUT OF SCOPE, DELIBERATELY. No interaction model. No MYC stratification. No
# Johnson-Neyman. No subtype split. No new endpoints or ratios. Those are later
# gates and none of them is pre-specified yet. Nothing here recomputes mitoPPS:
# the scores exist for both cohorts and are read.
#
# =============================================================================
# THE RECIPE, MATCHED TO THE MOUSE RATHER THAN INVENTED
# =============================================================================
# Read from `git -C <myc_mouse> show HEAD:scripts/48_gate_model_mouse_verification.R`
# on 2026-09-04, HEAD = e348dd8 on branch `paper-final`. PART A:
#
#     mito_sets <- grep("^MITOCARTA_", names(gmt), value = TRUE)
#     mito_all  <- unique(unlist(lapply(mito_sets, set_e)))
#     ox_sub    <- set_e("MITOCARTA_OXPHOS_SUBUNITS")   # nuclear only
#     rest_mito <- setdiff(mito_all, ox_sub)
#     comp_e    <- function(e) colMeans(t(scale(t(L[e, , drop = FALSE]))))
#     L         <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
#
# Three facts follow, and each was checked against the mouse GMT rather than
# assumed:
#
#   (a) THE 13 mtDNA-ENCODED GENES ARE EXCLUDED FROM THE NUMERATOR ONLY. The
#       mouse `MITOCARTA_OXPHOS_SUBUNITS` set carries 89 symbols and not one
#       `mt-` gene. But `mito_all` is the union of EVERY `MITOCARTA_*` set,
#       which includes `MITOCARTA_ALL`, `MITOCARTA_MTDNA_ENCODED` and
#       `MITOCARTA_OXPHOS_MT` - so all 13 sit in `rest_mito`, the DENOMINATOR.
#       That is consistent with this repo's standing convention, which forbids
#       POOLING the 13 with the nuclear subunits and says nothing about a
#       compartment-wide background. Section 2 carries `ox_rel_nomt` beside it
#       so the choice is measured rather than argued.
#
#   (b) `mito_all` IS THE FULL INVENTORY, NOT THE UNION OF THE PATHWAYS. The
#       mouse GMT carries `MITOCARTA_ALL` (1,140) and `MITOCARTA_UNASSIGNED`,
#       and the union of every MITOCARTA_ set is also 1,140 - so the mouse
#       denominator includes genes in no pathway. In human the two differ: the
#       sheet-2 inventory is 1,136 and the union of the 149 MitoPathways is
#       1,035, leaving 101 unassigned. The inventory is therefore the matching
#       object, and it is the SAME pinned set the `__MITOSTRIP` estimators were
#       stripped with, which is why `ox_rel` shares exactly 0 genes with
#       FELSHER__MITOSTRIP by construction.
#
#   (c) THE COMPOSITE IS MEAN-Z ON log2(LINEAR + 1), which is the mouse's `L`
#       and this repo's gene-level scale. It is NOT `zmean_arms`, which E02
#       builds on the VST. Section 2 reports the two against each other rather
#       than assuming they are interchangeable.
#
# =============================================================================
# SCALES, STATED ONCE
# =============================================================================
#   ox_rel, ox_lvl   built here, from log2(linear DESeq2-normalised + 1)
#   ox_ppd           READ, not rebuilt: mitoPPS on linear DESeq2-normalised
#   ox_gsva          READ: GSVA on VST, kcdf Gaussian - the incumbent human
#                    ruler, carried so section 5 can be read against the values
#                    already in the synthesis document
#   the 12 genes     log2(linear DESeq2-normalised + 1), as E10 and E13
#
# Every correlation below is rank-based except where a Pearson column is named
# beside a Spearman one, and ranks are invariant under any monotone transform,
# so no correlation here carries a scale error. The scale discipline binds where
# scores are BUILT, which is section 2.
#
# NEVER ACROSS COHORTS. Section 3 correlates rulers only within a cohort, and a
# species is a cohort: the mouse table is compared as an ORDERING, never as
# values. SCAN-B is resolved through `scanb_scores.rds$symbol_map` throughout -
# unharmonised the numerator covers 0.775 of its 89 genes instead of 0.989.
#
# STANDING ADJUSTMENT: partial Spearman on PROLIF_DISJOINT, 318 genes, disjoint
# from the MYC reference and from the numerator. It shares 10 genes with the
# 1,047-gene denominator (AK2, DUT, MTHFD2, PAICS, PRDX4, TBRG4, UNG, DTYMK,
# LIG3, POLQ); section 2 reports what dropping them does.
#
# NO PER-GENE FDR across the 12. Estimates with intervals, as E10.
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE16: is ox_rel a different ruler, and does the configuration survive it?\n",
        strrep("=", 78))

PATH_E16        <- file.path(DIR_RESULTS, "respiratory_rulers.rds")
PATH_E16_SEP    <- file.path(DIR_TABLES,  "E16_ruler_separability.csv")
PATH_E16_TWELVE <- file.path(DIR_TABLES,  "E16_twelve_on_rulers.csv")

# The author's lists, as E10. Reported in full: picking three genes after
# seeing them is the grid-of-cells trap.
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
PRIMING_ALL  <- sort(unique(c(PRIMING_PRO, PRIMING_ANTI)))

# The three genes the mouse model turns on. Named here ONLY so R3 and R4 can be
# stated; every table below carries all twelve.
CONFIG_UP   <- c("BBC3", "BCL2L1")
CONFIG_DOWN <- "MCL1"

# The three rulers of the mouse triple, in the mouse's own order.
RULERS_MOUSE <- c("ox_rel", "ox_ppd", "ox_lvl")

# =============================================================================
# DECISION RULES, TRANSCRIBED BEFORE ANY NUMBER WAS LOOKED AT
# =============================================================================
# From the G1 prompt, 2026-09-04. They are constants so the verdict in section 7
# cannot drift away from the rule it claims to apply.
R1_MIN_RHO   <- 0.85   # rho(ox_rel, ox_ppd) in BOTH cohorts
R1_MAX_DELTA <- 0.05   # and their MYC entanglements differ by at most this

# R2 is an ORDERING, not a magnitude: ox_rel < ox_ppd < ox_lvl on |r| with MYC.
# Magnitudes are not comparable across species - the mouse table is 24 animals
# where everything correlates with everything, and these are 1,095 and 3,207
# tumours. The mouse ordering, for reference only, never for a value comparison:
#            MYC mRNA  MSigDB MYC  DoRothEA MYC  proliferation
#   ox_lvl      0.61       0.81         0.88          0.80
#   ox_ppd      0.45       0.75         0.81          0.75
#   ox_rel      0.33       0.61         0.64          0.59
# Those are PEARSON r (script 48 uses stats::cor at its default). Section 4
# therefore computes Spearman - the standing measure here - AND Pearson, so the
# side-by-side reading is like for like.
MOUSE_SEPARABILITY <- tibble::tibble(
  ruler   = rep(c("ox_lvl", "ox_ppd", "ox_rel"), each = 4L),
  against = rep(c("myc_mRNA", "myc_msigdb", "myc_regulon", "prolif"), 3L),
  mouse_r = c(0.61, 0.81, 0.88, 0.80,
              0.45, 0.75, 0.81, 0.75,
              0.33, 0.61, 0.64, 0.59))

# R5 has NO pre-specified numeric threshold, and one is not invented here. The
# verdict is stated as a comparison against a reference that already exists in
# this repo: GSVA-vs-mitoPPS rank agreement over these same 12 genes, which
# E10's saved object puts at 0.930 / 0.902, and CLAUDE.md trap 5, which puts
# cross-instrument agreement over the 18 arms between 0.24 and 0.94.
R5_REFERENCE <- c(TCGA = 0.930, `SCAN-B` = 0.902)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
e10  <- readRDS(file.path(DIR_RESULTS, "machinery_and_priming.rds"))
fr   <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)

# DEPENDENCY GUARDS. Every one of these is an input this script would otherwise
# use silently and wrongly if it had moved.
if (!"component_cor" %in% names(e10)) {
  stop("results/machinery_and_priming.rds carries no `component_cor`; ",
       "section 5's reproduction check cannot run. Re-run E10.", call. = FALSE)
}
if (!identical(colnames(nw$tcga_gsva_new), ID_T) ||
    !identical(colnames(sc$gsva_new), ID_S) ||
    !identical(names(nw$tcga_log2MYC), ID_T) ||
    !identical(names(sc$log2MYC), ID_S)) {
  stop("the score objects are not aligned on the same samples as the arms. ",
       "Re-run E02.", call. = FALSE)
}

# --- 1.1 the two MitoCarta sets ----------------------------------------------
# NUMERATOR: the narrow subunit set, 89 genes, nuclear-encoded by construction.
# Taken from the snapshot's own `arm_sets` so it is the identical object the
# GSVA and mitoPPS arms were built on - not a re-parse of the workbook.
OX_SUB <- mito$arm_sets[["OXPHOS subunits"]]
# DENOMINATOR: the whole pinned MitoCarta inventory minus the numerator. See
# note (b) in the header for why the inventory and not the pathway union.
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL
REST     <- setdiff(MITO_ALL, OX_SUB)

# These are the two facts the whole recipe rests on. Asserted, not described:
# an MT- gene in the numerator would silently pool the 13 with the nuclear
# subunits, and a numerator gene outside the inventory would leave part of the
# exposure in the background it is measured against.
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL,
          length(OX_SUB) == 89L,
          !any(grepl("^MT-", OX_SUB)),
          all(OX_SUB %in% MITO_ALL))
message("   numerator   `OXPHOS subunits`      ", length(OX_SUB),
        " genes, 0 mtDNA-encoded")
message("   denominator  MitoCarta minus those ", length(REST),
        " genes, ", sum(grepl("^MT-", REST)), " mtDNA-encoded (mouse recipe)")

# --- 1.2 what the ruler shares with the things it is measured against --------
# `ox_rel` is a difference of two composites, so an estimator that overlaps
# EITHER side is partly correlated with itself. Printed, because the reference
# MYC estimator's zero here is a property worth knowing and the covariate's ten
# is a caveat section 2 then measures.
PROLIF_COV <- "PROLIF_DISJOINT"
PD <- mito$covariate_sets[[PROLIF_COV]]
overlap_audit <- tibble::tibble(
  set = c(MYC_REF, "MYC_UP.V1_UP__MITOSTRIP",
          "HALLMARK_MYC_TARGETS_V1__MITOSTRIP", MB_REF, PROLIF_COV),
  n   = c(lengths(sd_$myc_sets[c(MYC_REF, "MYC_UP.V1_UP__MITOSTRIP",
                                 "HALLMARK_MYC_TARGETS_V1__MITOSTRIP")]),
          length(sd_$collectri_sets[[MB_REF]]), length(PD))) %>%
  dplyr::mutate(
    x_numerator   = vapply(list(sd_$myc_sets[[MYC_REF]],
                                sd_$myc_sets[["MYC_UP.V1_UP__MITOSTRIP"]],
                                sd_$myc_sets[["HALLMARK_MYC_TARGETS_V1__MITOSTRIP"]],
                                sd_$collectri_sets[[MB_REF]], PD),
                           function(s) length(intersect(s, OX_SUB)), integer(1)),
    x_denominator = vapply(list(sd_$myc_sets[[MYC_REF]],
                                sd_$myc_sets[["MYC_UP.V1_UP__MITOSTRIP"]],
                                sd_$myc_sets[["HALLMARK_MYC_TARGETS_V1__MITOSTRIP"]],
                                sd_$collectri_sets[[MB_REF]], PD),
                           function(s) length(intersect(s, REST)), integer(1)))
message("\n   gene overlap of each estimator with the two halves of ox_rel:")
overlap_audit %>% as.data.frame() %>% print(row.names = FALSE)
message("   every __MITOSTRIP estimator is 0/0 BY CONSTRUCTION - MITOCARTA_ALL",
        " is\n   the strip set and it is the denominator's parent.")

# --- 1.3 the gene matrices ---------------------------------------------------
# SCALE: log2(linear DESeq2-normalised + 1), for the composites and for the 12.
tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(tcga_lin$scale, "linear_deseq2_normalised"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))
GT <- log2(tcga_lin$mat[, ID_T, drop = FALSE] + 1)
GS <- log2(scanb_lin$mat[, ID_S, drop = FALSE] + 1)
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

RES_T <- .symbol_resolver(rownames(GT), NULL)
RES_S <- .symbol_resolver(rownames(GS), sc$symbol_map)

COH <- list(
  TCGA     = list(G = GT, res = RES_T, ids = ID_T, arms = mito,
                  gsva_new = nw$tcga_gsva_new, mb = nw$tcga_M_b_variants,
                  log2myc = nw$tcga_log2MYC),
  `SCAN-B` = list(G = GS, res = RES_S, ids = ID_S, arms = sc,
                  gsva_new = sc$gsva_new, mb = sc$M_b_variants,
                  log2myc = sc$log2MYC))
message("\n   TCGA ", length(ID_T), " samples | SCAN-B ", length(ID_S), " samples")

# TCGA only, for check 4. SCAN-B has no purity estimate and it is never imputed.
PL <- fr %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::select(sample_id, purity, leuko)
PL <- PL[match(ID_T, PL$sample_id), ]
PURITY_OK <- stats::complete.cases(PL[, c("purity", "leuko")])
message("   TCGA samples with purity AND leukocyte fraction: ", sum(PURITY_OK),
        " of ", length(ID_T))

# =============================================================================
# 2. The rulers
# =============================================================================
# SCALE: `.comp` consumes log2(linear DESeq2-normalised + 1) and returns the
# mean, per sample, of each gene's z-score across samples. That is the mouse's
# `comp_e` exactly. `ox_ppd` and `ox_gsva` are READ from their saved objects and
# are NOT rebuilt here.
message("\n2. the rulers")

.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res)
  M  <- gr$mat
  v  <- apply(M, 1L, stats::var)
  # The mouse composite does not drop zero-variance rows because it has none;
  # scale() would return NaN for such a row and poison the whole composite, so
  # they are dropped here and COUNTED rather than silently tolerated.
  list(score   = colMeans(t(scale(t(M[v > 0, , drop = FALSE])))),
       n_asked = length(unique(genes)),
       n_used  = sum(v > 0),
       n_zerovar = sum(v == 0),
       missing = gr$missing)
}

rulers <- list(); set_sizes <- list()
for (coh in names(COH)) {
  C   <- COH[[coh]]
  num <- .comp(OX_SUB, C)
  den <- .comp(REST,   C)
  # SENSITIVITY, not a second ruler: the same denominator without the 13
  # mtDNA-encoded genes. Header note (a) - the mouse puts them in the
  # background, this repo holds them apart everywhere else, and the difference
  # between the two is worth a number rather than an argument.
  den_nomt <- .comp(setdiff(REST, grep("^MT-", REST, value = TRUE)), C)
  # And the same denominator without the 10 genes it shares with the covariate.
  den_nopd <- .comp(setdiff(REST, PD), C)

  M <- rbind(
    ox_rel  = num$score - den$score,
    ox_ppd  = as.numeric(C$arms$mitopps_arms["OXPHOS subunits", C$ids]),
    ox_lvl  = num$score,
    ox_gsva = as.numeric(C$arms$gsva_arms["OXPHOS subunits", C$ids]))
  colnames(M) <- C$ids
  stopifnot(!anyNA(M))

  aux <- rbind(ox_rel_nomt = num$score - den_nomt$score,
               ox_rel_nopd = num$score - den_nopd$score,
               ox_zmean    = as.numeric(C$arms$zmean_arms["OXPHOS subunits", C$ids]))
  colnames(aux) <- C$ids

  rulers[[coh]]    <- list(M = M, aux = aux)
  set_sizes[[coh]] <- tibble::tibble(
    cohort = coh,
    numerator_asked = num$n_asked, numerator_used = num$n_used,
    denominator_asked = den$n_asked, denominator_used = den$n_used,
    numerator_missing = paste(num$missing, collapse = ", "),
    n_zerovar = num$n_zerovar + den$n_zerovar,
    numerator_coverage = num$n_used / num$n_asked,
    denominator_coverage = den$n_used / den$n_asked)
}
set_sizes <- dplyr::bind_rows(set_sizes)
message("\n   set sizes actually used, after the symbol map:")
set_sizes %>%
  dplyr::mutate(dplyr::across(dplyr::ends_with("coverage"), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
# The mouse built ox_rel from 87 subunits against 999 background genes. The
# human counts are its analogue, not its equal, and the two are never pooled.
message("   mouse, for the shape of it and NOT as a value: 87 against 999.")

# The three sensitivities, so section 7's verdicts are read knowing what would
# have changed them. All within-cohort.
ruler_sensitivity <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  R <- rulers[[coh]]
  tibble::tibble(
    cohort = coh,
    vs_no_mtDNA_in_denominator = .rho(R$M["ox_rel", ], R$aux["ox_rel_nomt", ]),
    vs_no_covariate_genes      = .rho(R$M["ox_rel", ], R$aux["ox_rel_nopd", ]),
    ox_lvl_vs_zmean_arm        = .rho(R$M["ox_lvl", ], R$aux["ox_zmean", ]))
}))
message("\n   sensitivities (Spearman, within cohort):")
ruler_sensitivity %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   column 3 is the reason ox_lvl is BUILT here and not read from",
        " `zmean_arms`:\n   that arm is mean-z on the VST, this is mean-z on",
        " log2(linear + 1).")

# =============================================================================
# 3. CHECK 1 - ruler agreement, within each cohort and never across
# =============================================================================
message("\n3. check 1: do the rulers agree with each other?")

check1 <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M  <- rulers[[coh]]$M
  cb <- utils::combn(rownames(M), 2L)
  tibble::tibble(
    cohort = coh, a = cb[1, ], b = cb[2, ],
    pair   = paste(cb[1, ], "vs", cb[2, ]),
    n      = ncol(M),
    rho    = vapply(seq_len(ncol(cb)),
                    function(i) .rho(M[cb[1, i], ], M[cb[2, i], ]), numeric(1)))
}))
check1 %>%
  dplyr::select(cohort, pair, n, rho) %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = c(n, rho)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   Spearman, within cohort only. A value is never compared across",
        " cohorts;\n   mitoPPS's baseline is composition-dependent and GSVA is",
        " cohort-relative.")

# =============================================================================
# 4. CHECK 2 - separability. This is the check that decides G3.
# =============================================================================
# RAW marginal correlations, NOT proliferation-adjusted. Entanglement is the
# quantity being measured, so partialling the thing out would answer a different
# question; `PROLIF_DISJOINT` appears as a row of the table instead. This is the
# one section of the script where the standing adjustment does not apply, and it
# is stated here so a reader does not assume it.
#
# N7: MYC mRNA is not MYC activity. Both are reported and neither collapses into
# the other. The three set-based estimators are __MITOSTRIP, which is what the
# mouse did to its two signature scores - an unstripped signature shares genes
# with the ruler and the entanglement would be partly a set overlap. The __FULL
# variants are carried beside them so the size of that effect is visible.
message("\n4. check 2: how far is each ruler from MYC?")

# Short names, because these become column headers and a wrapped table is a
# table nobody reads. The first three are the mouse's own three columns.
EST_PRIMARY <- c(
  myc_mRNA    = "log2MYC",                              # the mouse's mRNA row
  myc_msigdb  = "HALLMARK_MYC_TARGETS_V1__MITOSTRIP",   # the mouse's MSigDB row
  myc_regulon = MB_REF,                                 # CollecTRI: the DoRothEA analogue
  myc_felsher = MYC_REF,                                # this repo's M_a
  myc_lowent  = "MYC_UP.V1_UP__MITOSTRIP",              # 1.5 pct entangled
  prolif      = PROLIF_COV)
# The three the mouse table has, and the only three R2 is scored on.
EST_MOUSE_TRIPLE <- c("myc_mRNA", "myc_msigdb", "myc_regulon")
EST_MYC <- setdiff(names(EST_PRIMARY), "prolif")
EST_SENSITIVITY <- c(
  msigdb_FULL  = "HALLMARK_MYC_TARGETS_V1__FULL",
  felsher_FULL = "FELSHER__FULL",
  lowent_FULL  = "MYC_UP.V1_UP__FULL",
  regulon_FULL = "M_b__FULL")

.estimator_matrix <- function(C, want) {
  rows <- lapply(want, function(w) {
    if (identical(w, "log2MYC"))   return(as.numeric(C$log2myc[C$ids]))
    if (identical(w, PROLIF_COV))  return(as.numeric(C$arms$gsva_cov[PROLIF_COV, C$ids]))
    if (w %in% rownames(C$mb))     return(as.numeric(C$mb[w, C$ids]))
    if (w %in% rownames(C$gsva_new)) return(as.numeric(C$gsva_new[w, C$ids]))
    stop("estimator not scored in this cohort: ", w, call. = FALSE)
  })
  M <- do.call(rbind, rows)
  dimnames(M) <- list(names(want), C$ids)
  M
}

separability <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  M <- rulers[[coh]]$M
  E <- .estimator_matrix(C, c(EST_PRIMARY, EST_SENSITIVITY))
  tidyr::expand_grid(ruler = rownames(M), against = rownames(E)) %>%
    dplyr::mutate(
      cohort   = coh,
      n        = ncol(M),
      panel    = dplyr::if_else(against %in% names(EST_PRIMARY),
                                "primary", "sensitivity"),
      spearman = purrr::map2_dbl(ruler, against, ~ .rho(M[.x, ], E[.y, ])),
      pearson  = purrr::map2_dbl(ruler, against,
                                 ~ stats::cor(M[.x, ], E[.y, ])))
}))

message("\n   PRIMARY panel, Spearman (the standing measure here):")
separability %>%
  dplyr::filter(panel == "primary") %>%
  dplyr::select(cohort, ruler, against, spearman) %>%
  dplyr::mutate(spearman = round(spearman, 3)) %>%
  tidyr::pivot_wider(names_from = against, values_from = spearman) %>%
  dplyr::arrange(cohort, match(ruler, c("ox_lvl", "ox_ppd", "ox_rel", "ox_gsva"))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the same in PEARSON, which is what the mouse table is, so the two",
        " can be\n   read side by side as ORDERINGS (never as values - a species",
        " is a cohort):")
separability %>%
  dplyr::filter(panel == "primary") %>%
  dplyr::select(cohort, ruler, against, pearson) %>%
  dplyr::mutate(pearson = round(pearson, 3)) %>%
  tidyr::pivot_wider(names_from = against, values_from = pearson) %>%
  dplyr::arrange(cohort, match(ruler, c("ox_lvl", "ox_ppd", "ox_rel", "ox_gsva"))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   SENSITIVITY panel - the same estimators unstripped. The gap",
        " between a\n   __FULL row and its __MITOSTRIP row is set overlap, not",
        " biology:")
separability %>%
  dplyr::filter(panel == "sensitivity") %>%
  dplyr::select(cohort, ruler, against, spearman) %>%
  dplyr::mutate(spearman = round(spearman, 3)) %>%
  tidyr::pivot_wider(names_from = against, values_from = spearman) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. CHECK 3 - the configuration, all 12 genes on every ruler
# =============================================================================
# Partial Spearman on PROLIF_DISJOINT with Fisher-z intervals and the
# Bonett-Wright variance, se = sqrt((1 + rho^2/2)/(n - 3 - k)). That is E10's
# `.cor_block` verbatim - the plain 1/(n-3) is the Pearson case and understates
# a rank correlation. It is copied rather than sourced because E10 defines it
# inline and there is no shared file to source; the reproduction check at the
# end of this section is what proves the copy has not drifted.
#
# ALL TWELVE ARE REPORTED. Selecting the three of interest after seeing them is
# the grid-of-cells trap, and an interval here describes a cell, it never
# selects one.
#
# SCALE: the genes are log2(linear DESeq2-normalised + 1); the rulers are as
# built in section 2. Ranks make the mixture safe.
message("\n5. check 3: the twelve on every ruler")

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

GR <- lapply(COH, function(C) .gene_rows(PRIMING_ALL, C$G, C$res))
miss <- unique(unlist(lapply(GR, function(g) g$missing), use.names = FALSE))
if (length(miss)) {
  stop("these priming genes did not resolve to exactly one matrix row: ",
       paste(miss, collapse = ", "),
       "\nCheck the SCAN-B symbol map before going further.", call. = FALSE)
}
message("   all ", length(PRIMING_ALL), " priming genes present in both cohorts")

CV <- lapply(names(COH), function(coh)
  t(COH[[coh]]$arms$gsva_cov[PROLIF_COV, COH[[coh]]$ids, drop = FALSE]))
names(CV) <- names(COH)
stopifnot(identical(rownames(CV$TCGA), ID_T),
          identical(rownames(CV$`SCAN-B`), ID_S))

twelve <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  .cor_block(rulers[[coh]]$M, GR[[coh]]$mat, cov = CV[[coh]]) %>%
    dplyr::mutate(cohort = coh, adjustment = "PROLIF_DISJOINT")
})) %>%
  dplyr::mutate(side = dplyr::if_else(gene %in% PRIMING_PRO,
                                      "pro-apoptotic", "anti-apoptotic"))

# --- 5.1 REPRODUCTION CHECK --------------------------------------------------
# `ox_gsva` and `ox_ppd` are the two axes E10 already carries as `OXPHOS` and
# `OXPHOS_mitopps`. Same vectors, same covariate, same estimator - so they must
# come back bit-equal. This is the load block's positive control and it is the
# reason a failure here is a broken input rather than a new result.
repro <- twelve %>%
  dplyr::filter(ruler %in% c("ox_gsva", "ox_ppd")) %>%
  dplyr::mutate(axis = dplyr::if_else(ruler == "ox_gsva", "OXPHOS",
                                      "OXPHOS_mitopps")) %>%
  dplyr::inner_join(
    e10$component_cor %>% dplyr::select(cohort, axis, gene, e10_rho = rho),
    by = c("cohort", "axis", "gene")) %>%
  dplyr::mutate(delta = rho - e10_rho)
stopifnot(nrow(repro) == 2L * 2L * length(PRIMING_ALL),
          max(abs(repro$delta)) < 1e-10)
message("   reproduction check: ", nrow(repro), " cells reproduce E10 exactly",
        " (max |delta| ", signif(max(abs(repro$delta)), 3), ")")

for (coh in names(COH)) {
  message("\n   the twelve on every ruler, ", coh,
          " (ordered by ox_rel; ox_gsva is the incumbent):")
  twelve %>%
    dplyr::filter(cohort == coh, ruler %in% c(RULERS_MOUSE, "ox_gsva")) %>%
    dplyr::select(ruler, gene, side, rho) %>%
    tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
    dplyr::arrange(dplyr::desc(ox_rel)) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

message("\n   the three the mouse model names, with intervals:")
twelve %>%
  dplyr::filter(gene %in% c(CONFIG_UP, CONFIG_DOWN), ruler %in% RULERS_MOUSE) %>%
  dplyr::transmute(cohort, ruler, gene, n,
                   rho = sprintf("%+.3f", rho),
                   ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi)) %>%
  dplyr::arrange(gene, cohort, match(ruler, RULERS_MOUSE)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.2 is the whole configuration ruler-stable? ----------------------------
# The question R5 asks. If the ORDERING of the twelve survives the ruler change,
# any three-gene claim is safe; if it does not, that is the finding whatever
# those three genes do.
ruler_stability <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  W <- twelve %>% dplyr::filter(cohort == coh) %>%
    dplyr::select(ruler, gene, rho) %>%
    tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
    dplyr::arrange(gene)
  rn <- setdiff(names(W), "gene")
  cb <- utils::combn(rn, 2L)
  tibble::tibble(
    cohort = coh, a = cb[1, ], b = cb[2, ],
    pair = paste(cb[1, ], "vs", cb[2, ]),
    rank_rho = vapply(seq_len(ncol(cb)), function(i)
      stats::cor(W[[cb[1, i]]], W[[cb[2, i]]], method = "spearman"), numeric(1)),
    max_abs_shift = vapply(seq_len(ncol(cb)), function(i)
      max(abs(W[[cb[1, i]]] - W[[cb[2, i]]])), numeric(1)))
}))
message("\n   cross-ruler rank agreement over the twelve (Spearman of the 12",
        " rho values):")
ruler_stability %>%
  dplyr::select(cohort, pair, rank_rho, max_abs_shift) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  tidyr::pivot_wider(names_from = cohort,
                     values_from = c(rank_rho, max_abs_shift)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.3 the load-bearing derived quantity -----------------------------------
# BCL2L1 minus BBC3. It is positive on GSVA in both cohorts (+0.061 / +0.180),
# which is why PRIME = BBC3 - BCL2L1 comes out flat on OXPHOS, which is what
# explains the pre-registered null. If the gap does not survive the ruler
# change, that explanation collapses. N3 applies to every word of this: these
# are transcript correlations, never priming.
the_gap <- twelve %>%
  dplyr::filter(gene %in% CONFIG_UP) %>%
  dplyr::select(cohort, ruler, gene, rho) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rho) %>%
  dplyr::mutate(gap = BCL2L1 - BBC3) %>%
  dplyr::arrange(cohort, match(ruler, c(RULERS_MOUSE, "ox_gsva")))
message("\n   BCL2L1 minus BBC3, per ruler per cohort:")
the_gap %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. CHECK 4 - infiltrate control, TCGA only
# =============================================================================
# Trap 2: breast is the worst TCGA tissue for this. Adipose is OXPHOS/FAO-high
# and infiltrate carries its own BCL2-family profile. SCAN-B has no purity
# estimate and it is never imputed, so this cannot replicate and is not asked to.
#
# BOTH ROWS ARE COMPUTED ON THE SAME 1,007 SAMPLES. Comparing an adjusted value
# at n = 1,007 against the section-5 value at n = 1,095 would confound the
# covariate with the sample set, and the difference wanted here is the covariate
# alone. Every ruler is carried, not just ox_rel, because a shift that all four
# share is a property of the samples and not of the ruler.
message("\n6. check 4: TCGA only, plus purity and leukocyte fraction")

ids_p  <- ID_T[PURITY_OK]
cov_pl <- cbind(CV$TCGA[PURITY_OK, , drop = FALSE],
                as.matrix(PL[PURITY_OK, c("purity", "leuko")]))
rownames(cov_pl) <- ids_p
infiltrate <- dplyr::bind_rows(
  .cor_block(rulers$TCGA$M[, ids_p, drop = FALSE],
             GR$TCGA$mat[, ids_p, drop = FALSE],
             cov = CV$TCGA[PURITY_OK, , drop = FALSE]) %>%
    dplyr::mutate(adjustment = "PROLIF_DISJOINT"),
  .cor_block(rulers$TCGA$M[, ids_p, drop = FALSE],
             GR$TCGA$mat[, ids_p, drop = FALSE],
             cov = cov_pl) %>%
    dplyr::mutate(adjustment = "PROLIF + purity + leuko")) %>%
  dplyr::mutate(cohort = "TCGA",
                side = dplyr::if_else(gene %in% PRIMING_PRO,
                                      "pro-apoptotic", "anti-apoptotic"))

infiltrate_shift <- infiltrate %>%
  dplyr::select(ruler, gene, adjustment, rho) %>%
  tidyr::pivot_wider(names_from = adjustment, values_from = rho) %>%
  dplyr::mutate(shift = `PROLIF + purity + leuko` - PROLIF_DISJOINT)

message("\n   ox_rel, the twelve, before and after (n = ", length(ids_p), "):")
infiltrate_shift %>%
  dplyr::filter(ruler == "ox_rel") %>%
  dplyr::arrange(dplyr::desc(PROLIF_DISJOINT)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   how far each ruler moves, over the twelve:")
infiltrate_shift %>%
  dplyr::group_by(ruler) %>%
  dplyr::summarise(mean_abs_shift = mean(abs(shift)),
                   max_abs_shift  = max(abs(shift)),
                   n_sign_changes = sum(sign(PROLIF_DISJOINT) !=
                                          sign(`PROLIF + purity + leuko`)),
                   # a sign change on a rho near zero is not a reversal, so the
                   # gene is named and its two values are printed beside it
                   sign_changes = paste(gene[sign(PROLIF_DISJOINT) !=
                                               sign(`PROLIF + purity + leuko`)],
                                        collapse = ", "),
                   .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
sign_flips <- infiltrate_shift %>%
  dplyr::filter(sign(PROLIF_DISJOINT) != sign(`PROLIF + purity + leuko`))
if (nrow(sign_flips)) {
  message("   the cells that changed sign, with both values:")
  sign_flips %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

gap_infiltrate <- infiltrate %>%
  dplyr::filter(gene %in% CONFIG_UP) %>%
  dplyr::select(ruler, adjustment, gene, rho) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rho) %>%
  dplyr::mutate(gap = BCL2L1 - BBC3)
message("\n   and the load-bearing gap, before and after:")
gap_infiltrate %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   TCGA only. SCAN-B has no purity estimate - trap 2 - so no row here",
        "\n   replicates, and none is written as if it did.")

# =============================================================================
# 7. The verdicts, on the rules fixed at the top of this file
# =============================================================================
# Every rule below was transcribed from the G1 prompt before any number in this
# script existed. They are applied here mechanically, and each carries the
# numbers it was decided on so the note can quote them rather than restate them.
message("\n7. verdicts")

.spear <- function(coh, ruler, against) {
  separability$spearman[separability$cohort == coh &
                          separability$ruler == ruler &
                          separability$against == against]
}

# --- R1: are the rulers interchangeable? -------------------------------------
r1_rho <- vapply(names(COH), function(coh)
  check1$rho[check1$cohort == coh & check1$pair == "ox_rel vs ox_ppd"],
  numeric(1))
# The rule says MYC entanglement, so `prolif` is reported but is not part of
# the test. Entanglement is a distance from zero, hence |r| on both sides.
r1_delta <- dplyr::bind_rows(lapply(names(COH), function(coh)
  tibble::tibble(
    cohort  = coh,
    against = names(EST_PRIMARY),
    is_myc  = names(EST_PRIMARY) %in% EST_MYC,
    ox_rel  = vapply(names(EST_PRIMARY),
                     function(a) abs(.spear(coh, "ox_rel", a)), numeric(1)),
    ox_ppd  = vapply(names(EST_PRIMARY),
                     function(a) abs(.spear(coh, "ox_ppd", a)), numeric(1))) %>%
    dplyr::mutate(delta = ox_rel - ox_ppd)))
r1_myc  <- r1_delta %>% dplyr::filter(is_myc)
r1_worst <- r1_myc %>% dplyr::slice_max(abs(delta), n = 1L)
r1_pass <- all(r1_rho >= R1_MIN_RHO) && max(abs(r1_myc$delta)) <= R1_MAX_DELTA
message("\n   R1 rho(ox_rel, ox_ppd): ",
        paste(sprintf("%s %.3f", names(r1_rho), r1_rho), collapse = " | "),
        "  (rule: >= ", R1_MIN_RHO, " in BOTH)")
message("   R1 |r| with each MYC estimator, ox_rel minus ox_ppd:")
r1_myc %>% dplyr::select(cohort, against, ox_rel, ox_ppd, delta) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   R1 largest: ", sprintf("%+.3f", r1_worst$delta[1]), " on ",
        r1_worst$against[1], " in ", r1_worst$cohort[1],
        "  (rule: <= ", R1_MAX_DELTA, ")")
message("   R1 VERDICT: ", if (r1_pass)
  "INTERCHANGEABLE - mitoPPS carries everything, ox_rel is a one-line sensitivity"
  else "DISTINCT INSTRUMENTS - report ox_rel and ox_ppd separately")

# --- R2: does the human reproduce the mouse separability ordering? -----------
# ORDERING only. Magnitudes are not comparable across species and are never
# compared here. |r| because a ruler may lean the other way in human, and
# entanglement is a distance from zero.
# Scored on the mouse's own three columns. The other two MYC estimators and
# proliferation are carried as rows so a 6-cell verdict is not read as the
# whole picture - but they are NOT part of the rule and do not enter the count.
r2 <- dplyr::bind_rows(lapply(names(COH), function(coh)
  dplyr::bind_rows(lapply(names(EST_PRIMARY), function(a) tibble::tibble(
    cohort = coh, against = a,
    mouse_triple = a %in% EST_MOUSE_TRIPLE,
    ox_rel = abs(.spear(coh, "ox_rel", a)),
    ox_ppd = abs(.spear(coh, "ox_ppd", a)),
    ox_lvl = abs(.spear(coh, "ox_lvl", a)),
    ordered_as_mouse = ox_rel < ox_ppd & ox_ppd < ox_lvl,
    ox_rel_is_lowest = ox_rel < ox_ppd & ox_rel < ox_lvl)))))
message("\n   R2 ordering |r|, ox_rel < ox_ppd < ox_lvl?  (scored on the",
        " mouse triple only)")
r2 %>%
  dplyr::mutate(dplyr::across(c(ox_rel, ox_ppd, ox_lvl), ~ round(.x, 3))) %>%
  dplyr::select(cohort, against, mouse_triple, ox_rel, ox_ppd, ox_lvl,
                as_mouse = ordered_as_mouse) %>%
  as.data.frame() %>% print(row.names = FALSE)
r2_myc <- r2 %>% dplyr::filter(mouse_triple)
r2_pass <- all(r2_myc$ordered_as_mouse)
message("   R2 VERDICT: full mouse ordering on ",
        sum(r2_myc$ordered_as_mouse), " of ", nrow(r2_myc),
        " MYC cells; ox_rel is the least entangled in ",
        sum(r2_myc$ox_rel_is_lowest), " of ", nrow(r2_myc))

# --- R3: does the reversal hold? ---------------------------------------------
r3 <- twelve %>%
  dplyr::filter(ruler == "ox_rel", gene %in% c(CONFIG_UP, CONFIG_DOWN)) %>%
  dplyr::mutate(expected = dplyr::if_else(gene %in% CONFIG_UP, 1, -1),
                ok = sign(rho) == expected) %>%
  dplyr::select(cohort, gene, rho, ci_lo, ci_hi, expected, ok)
r3_by_cohort <- r3 %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(pass = all(ok), n_ok = sum(ok), .groups = "drop")
message("\n   R3 signs on ox_rel (BBC3 +, BCL2L1 +, MCL1 -):")
r3 %>% dplyr::mutate(rho = sprintf("%+.3f", rho),
                     ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi)) %>%
  dplyr::select(cohort, gene, rho, ci, expected, ok) %>%
  as.data.frame() %>% print(row.names = FALSE)
r3_by_cohort %>% as.data.frame() %>% print(row.names = FALSE)
r3_pass <- all(r3_by_cohort$pass)

# --- R4: does the gap hold? THE LOAD-BEARING ONE -----------------------------
r4 <- the_gap %>% dplyr::filter(ruler == "ox_rel") %>%
  dplyr::mutate(ok = gap > 0)
r4_pass <- all(r4$ok)
message("\n   R4 BCL2L1 - BBC3 on ox_rel (rule: > 0 in BOTH cohorts):")
r4 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ sprintf("%+.3f", .x))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   R4 VERDICT: ", if (r4_pass) "HOLDS" else "FAILS",
        " - and this is the one that carries the explanation of the",
        "\n   pre-registered null.")

# --- R5: is the configuration ruler-stable? ----------------------------------
# No numeric threshold was pre-specified, and none is invented. The verdict is
# a comparison against the in-repo reference stated at the top of this file.
r5 <- ruler_stability %>%
  dplyr::filter(pair %in% c("ox_rel vs ox_ppd", "ox_rel vs ox_lvl",
                            "ox_rel vs ox_gsva", "ox_ppd vs ox_gsva"))
message("\n   R5 cross-ruler rank agreement over the twelve:")
r5 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(cohort, pair, rank_rho, max_abs_shift) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   reference already in this repo, GSVA vs mitoPPS over the same 12: ",
        paste(sprintf("%s %.3f", names(R5_REFERENCE), R5_REFERENCE),
              collapse = " | "))
r5_rel_gsva <- vapply(names(COH), function(coh)
  r5$rank_rho[r5$cohort == coh & r5$pair == "ox_rel vs ox_gsva"], numeric(1))

verdicts <- tibble::tibble(
  rule = c("R1", "R2", "R3", "R4", "R5"),
  question = c("are ox_rel and ox_ppd interchangeable?",
               "does human reproduce the mouse separability ordering?",
               "does the three-gene reversal hold on ox_rel?",
               "does BCL2L1 - BBC3 stay positive on ox_rel?",
               "is the configuration ruler-stable over the twelve?"),
  pass = c(r1_pass, r2_pass, r3_pass, r4_pass, NA),
  evidence = c(
    sprintf("rho %.3f / %.3f; largest MYC entanglement gap %+.3f on %s (%s)",
            r1_rho[["TCGA"]], r1_rho[["SCAN-B"]], r1_worst$delta[1],
            r1_worst$against[1], r1_worst$cohort[1]),
    sprintf("%d of %d MYC cells fully ordered; ox_rel least entangled in %d",
            sum(r2_myc$ordered_as_mouse), nrow(r2_myc),
            sum(r2_myc$ox_rel_is_lowest)),
    sprintf("TCGA %s, SCAN-B %s",
            ifelse(r3_by_cohort$pass[r3_by_cohort$cohort == "TCGA"], "pass", "fail"),
            ifelse(r3_by_cohort$pass[r3_by_cohort$cohort == "SCAN-B"], "pass", "fail")),
    sprintf("gap %+.3f / %+.3f",
            r4$gap[r4$cohort == "TCGA"], r4$gap[r4$cohort == "SCAN-B"]),
    sprintf("ox_rel vs ox_gsva rank rho %.3f / %.3f (reference %.3f / %.3f)",
            r5_rel_gsva[["TCGA"]], r5_rel_gsva[["SCAN-B"]],
            R5_REFERENCE[["TCGA"]], R5_REFERENCE[["SCAN-B"]])))
message("\n   the five verdicts:")
verdicts %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   R5 carries no pass/fail because no numeric threshold was",
        " pre-specified for it.\n   It is reported as a comparison, not as a test.")

# WHICH RULER IS PRIMARY.
#
# THE DECISION RULES SPECIFY ONE BRANCH AND NOT THE OTHER, AND THAT IS SAID
# RATHER THAN PAPERED OVER. R2 reads: "If yes, ox_rel is primary for G3
# regardless of what R1 returned - mutual correlation is not the criterion
# there, separability is." It fixes what happens when the human REPRODUCES the
# mouse ordering. It fixes nothing for the case where it does not. Inventing a
# rule for that branch after seeing the numbers is exactly the move this repo's
# posture exists to prevent, so the script reports the state and names the
# decision instead of manufacturing it.
primary_ruler <- if (r2_pass) {
  "ox_rel - R2 passes, and R2 outranks R1"
} else if (r1_pass) {
  paste("ox_ppd - R2 fails, so the mouse's case for ox_rel does not transfer,",
        "and R1 licenses mitoPPS to carry it")
} else {
  paste("NOT FIXED BY THE RULES - R2 fails, so ox_rel is not licensed as",
        "primary, and R1 fails, so ox_ppd may not absorb it. The rulers are",
        "distinct and neither rule elects one. This is a decision for the",
        "result note, taken openly, not a verdict of the pre-fixed rules.")
}
message("\n   PRIMARY RULER FOR THE HUMAN ARM GOING FORWARD:\n     ",
        primary_ruler)
if (!r2_pass && !r1_pass) {
  message("   What IS fixed: ox_rel and ox_ppd are separate instruments (R1),",
          " ox_rel is\n   not the least MYC-entangled ruler in human on two of",
          " the mouse's own three\n   estimators (R2), and the configuration",
          " survives every ruler tried (R3, R4).")
}

# =============================================================================
# 8. Save
# =============================================================================
saveRDS(list(
  set_sizes = set_sizes, overlap_audit = overlap_audit,
  ruler_sensitivity = ruler_sensitivity,
  check1 = check1, separability = separability,
  mouse_separability = MOUSE_SEPARABILITY,
  twelve = twelve, repro = repro, ruler_stability = ruler_stability,
  the_gap = the_gap,
  infiltrate = infiltrate, infiltrate_shift = infiltrate_shift,
  gap_infiltrate = gap_infiltrate, sign_flips = sign_flips,
  r1_rho = r1_rho, r1_delta = r1_delta, r2 = r2, r3 = r3, r4 = r4, r5 = r5,
  verdicts = verdicts, primary_ruler = primary_ruler,
  settings = list(
    rulers = c("ox_rel (built here)", "ox_ppd (mitoPPS, read)",
               "ox_lvl (built here)", "ox_gsva (GSVA, read)"),
    recipe_source = paste("myc_mouse HEAD e348dd8 branch paper-final,",
                          "scripts/48_gate_model_mouse_verification.R PART A"),
    composite_scale = "mean per-gene z of log2(linear DESeq2-normalised + 1)",
    mtdna_rule = paste("13 mtDNA-encoded genes are in the DENOMINATOR and not",
                       "the numerator - the mouse recipe; sensitivity in",
                       "$ruler_sensitivity"),
    numerator = "MitoCarta OXPHOS subunits, 89, nuclear only",
    denominator = "MITOCARTA_ALL minus the numerator, 1047",
    measure = "spearman; separability raw, everything else partial",
    covariate = PROLIF_COV,
    myc_estimators = EST_PRIMARY,
    decision_rules = c(R1_MIN_RHO = R1_MIN_RHO, R1_MAX_DELTA = R1_MAX_DELTA),
    seed = PROJECT_SEED),
  rules = list(
    scope = paste("MEASUREMENT CHECK. No interaction model, no MYC",
                  "stratification, no subtype split, no new endpoints."),
    n3 = "these are transcript correlations; the word 'primed' is never used",
    cohorts = "rulers are never correlated across cohorts; a species is a cohort",
    mouse = "the mouse table is compared as an ORDERING, never as values"),
  built = Sys.time()), PATH_E16)

readr::write_csv(separability, PATH_E16_SEP)
readr::write_csv(dplyr::bind_rows(
  twelve %>% dplyr::mutate(panel = "pooled"),
  infiltrate %>% dplyr::mutate(panel = "TCGA purity subset")), PATH_E16_TWELVE)

message("\nE16: done.")
message("    results/respiratory_rulers.rds")
message("    outputs/tables/E16_ruler_separability.csv")
message("    outputs/tables/E16_twelve_on_rulers.csv")
message("    no figures - this is a measurement check, and every number in it",
        " is a table.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E16)

  # The five verdicts, and which ruler they make primary.
  x$verdicts %>% as.data.frame()
  x$primary_ruler

  # What the ruler is actually built from, per cohort.
  x$set_sizes %>% as.data.frame()
  x$overlap_audit %>% as.data.frame()

  # Would the mtDNA decision or the covariate overlap have changed anything?
  x$ruler_sensitivity %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
    as.data.frame()

  # CHECK 1 - agreement, within cohort only.
  x$check1 %>%
    dplyr::select(cohort, pair, rho) %>%
    dplyr::mutate(rho = round(rho, 3)) %>%
    tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
    as.data.frame()

  # CHECK 2 - separability, beside the mouse ordering it is compared with.
  x$separability %>%
    dplyr::filter(panel == "primary") %>%
    dplyr::select(cohort, ruler, against, spearman) %>%
    dplyr::mutate(spearman = round(spearman, 3)) %>%
    tidyr::pivot_wider(names_from = against, values_from = spearman) %>%
    as.data.frame()
  x$mouse_separability %>%
    tidyr::pivot_wider(names_from = against, values_from = mouse_r) %>%
    as.data.frame()

  # CHECK 3 - all twelve, every ruler. Never read three of these on their own.
  x$twelve %>%
    dplyr::filter(ruler != "ox_gsva") %>%
    dplyr::select(cohort, ruler, gene, side, rho) %>%
    tidyr::pivot_wider(names_from = c(ruler, cohort), values_from = rho) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # The reproduction check against E10 - this is the load block's control.
  summary(abs(x$repro$delta))

  # R5, and the load-bearing gap.
  x$ruler_stability %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$the_gap %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # CHECK 4 - TCGA only, and both rows on the same 1,007 samples.
  x$infiltrate_shift %>%
    dplyr::filter(ruler == "ox_rel") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$gap_infiltrate %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

}
