# E20_share_vs_content.R
# =============================================================================
# A RE-READ OF E16, PLUS ONE NEW RULER. NOT A NEW ANALYSIS.
#
# E16's separability table (section 3 of docs/2026-09-04_e16_respiratory_rulers.md)
# correlates four respiratory rulers against five MYC estimators. It was read for
# ONE question - which ruler is least MYC-entangled - and it answers a second one
# directly, which nobody has asked it:
#
#     DOES MYC TRACK COMPARTMENT CONTENT, OR COMPARTMENT SHARE?
#
# The three rulers of the panel are THE SAME NUMERATOR WITH THREE DENOMINATORS,
# and that is the whole design:
#
#     ox_lvl        = comp(OXPHOS subunits)                        CONTENT
#     ox_rel        = comp(OXPHOS subunits) - comp(rest of MitoCarta)   SHARE
#     ox_nuc_mtrib  = comp(OXPHOS subunits) - comp(mitoribosome)    SHARE, but
#                     against the one denominator the mouse's developmental
#                     result is actually about.  <-- THE ONLY NEW COMPUTATION
#
# So the contrast rho(MYC, ox_lvl) - rho(MYC, ox_rel) isolates the DENOMINATOR
# and nothing else. Nothing is re-scored: mitoPPS and GSVA are read from their
# saved objects, and the two composites E16 already built are rebuilt to E16's
# own recipe and then CHECKED against E16's saved numbers (section 3). A
# reproduction failure there means this is not a re-read and the note may not be
# written.
#
# WHY IT MATTERS, stated so the reading cannot drift into it afterwards.
# Development lowers the OXPHOS SHARE of the mitochondrial compartment between
# 6W and 12W in the mouse while compartment CONTENT rises. If MYC in human
# tumours raises CONTENT rather than SHARE, then absolute respiration goes up
# and the human MYC-OXPHOS correlation appears WITHOUT the developmental
# reprioritisation ever being reversed. It was overwhelmed, not undone. That is
# a cleaner claim than "MYC reverses the trend" and it needs no reversal
# observed anywhere. Section 0 fixes what would license it and what would not,
# BEFORE any number exists.
#
# =============================================================================
# THE MAIN RISK, NAMED FIRST: E16 SECTION 3.1
# =============================================================================
# E16 established that separability from MYC is NOT a property of a ruler in
# human. It is a property of the RULER-ESTIMATOR PAIR: `ox_rel` buys separability
# from every signature estimator and LOSES it against the CollecTRI regulon, in
# both cohorts. That is CLAUDE.md trap 3 on a new axis.
#
# Therefore: EVERY claim in this script names its estimator, and a claim that
# does not hold across estimators is reported as ESTIMATOR-DEPENDENT, not as a
# finding. If the five estimators disagree about share-versus-content, THAT IS
# THE RESULT and section 6 stops there.
#
# =============================================================================
# SCALES, STATED ONCE
# =============================================================================
#   ox_lvl, ox_rel, ox_nuc_mtrib, mtrib_lvl
#                    BUILT here, from log2(linear DESeq2-normalised + 1), the
#                    mouse's `L` and E16's `.comp`, copied verbatim
#   ox_ppd           READ, not rebuilt: mitoPPS on linear DESeq2-normalised
#   ox_gsva          READ: GSVA on VST, kcdf Gaussian - the incumbent
#   estimators       READ: GSVA on VST (signatures), ULM (regulon), log2 of the
#                    linear matrix (MYC mRNA)
#
# Every correlation is rank-based, and ranks are invariant under any monotone
# transform, so no correlation here carries a scale error. The scale discipline
# binds where scores are BUILT, which is section 2.
#
# NEVER ACROSS COHORTS, and a species is a cohort. Rulers are correlated with
# estimators only WITHIN a cohort. The mouse is not touched by this script at
# all - not read, not compared, not cited as a value. SCAN-B is resolved through
# `scanb_scores.rds$symbol_map` throughout.
#
# ADJUSTMENT: the panel is reported BOTH raw and partial-Spearman on
# PROLIF_DISJOINT (318 genes). The VERDICT is scored on the ADJUSTED panel, as
# E16 section 4 adjusts; the raw panel is the literal re-read of E16 and is what
# section 3's reproduction check is run against.
#
# N3 THROUGHOUT. Every number here is a transcript correlation. The word
# "primed" appears nowhere as a description of a tumour, a cell or an animal.
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE20: does MYC track compartment CONTENT or compartment SHARE?\n",
        strrep("=", 78))

PATH_E20       <- file.path(DIR_RESULTS, "share_vs_content.rds")
PATH_E20_PANEL <- file.path(DIR_TABLES,  "E20_share_vs_content.csv")
PATH_E16       <- file.path(DIR_RESULTS, "respiratory_rulers.rds")

# =============================================================================
# 0. THE READING RULES, FIXED BEFORE ANY NUMBER IN THIS SCRIPT EXISTS
# =============================================================================
# These are constants so that section 6's verdict cannot drift away from the
# rule it claims to apply. Each one is either transcribed from a rule this repo
# already carries, or is stated here as new and says what it is for.

# The three rulers the verdict is scored on. `mtrib_lvl`, `ox_gsva` and `ox_ppd`
# are carried for context and are NOT part of the verdict.
PANEL_RULERS   <- c("ox_lvl", "ox_rel", "ox_nuc_mtrib")
CONTEXT_RULERS <- c("mtrib_lvl", "ox_gsva", "ox_ppd")

# ADJUDICATION FLOOR - new, and stated as new.
# An estimator can only speak to content-versus-share if it has a non-trivial
# association with the CONTENT ruler in the first place. The difference between
# two near-zero correlations has a direction and no meaning. E16 section 3.1
# made exactly this argument in prose about `myc_mRNA` (0.103 against 0.032 in
# TCGA; 0.048 against 0.011 in SCAN-B) and declined to lean on it. Here it is a
# number, fixed in advance rather than applied after the fact. Every estimator
# is REPORTED whatever it does; the floor governs only which ones ADJUDICATE.
ADJUDICATE_FLOOR <- 0.15   # |raw rho(estimator, ox_lvl)|, required in BOTH cohorts

# MATERIALITY - borrowed, not invented. E16's R1_MAX_DELTA is this repo's
# existing bar for "two entanglements differ by enough to matter", and the
# quantity here is the same kind of object: a difference between two of that
# table's cells.
MATERIAL_DELTA <- 0.05

# Bootstrap for the CONTRAST. The two rhos in `rho(est, ox_lvl) - rho(est,
# ox_rel)` are computed on the SAME samples against rulers that share their
# numerator and correlate at 0.94 / 0.92 (E16 check 1). Two independent Fisher-z
# intervals would therefore be badly wrong about their difference - far too
# wide - so the difference gets a PAIRED sample-level bootstrap and the marginal
# intervals are never used to judge it. Percentile interval, descriptive, not a
# test.
BOOT_B <- 2000L

# THE VERDICT RULE.
#   CONTENT  if, in EVERY adjudicating estimator-by-cohort cell,
#              rho(est, ox_lvl) > rho(est, ox_rel) AND
#              rho(est, ox_lvl) > rho(est, ox_nuc_mtrib)
#   SHARE    if the reverse holds in every such cell
#   ESTIMATOR-DEPENDENT otherwise - name which estimators disagree, and stop.
# MATERIAL is scored separately: it requires |delta| >= MATERIAL_DELTA AND a
# bootstrap interval excluding 0, again in every adjudicating cell.
VERDICT_RULE <- paste0(
  "CONTENT if rho(ox_lvl) exceeds BOTH share rulers in every adjudicating ",
  "estimator-by-cohort cell; SHARE if the reverse in every such cell; ",
  "ESTIMATOR-DEPENDENT otherwise. MATERIAL is separate: |delta| >= ",
  MATERIAL_DELTA, " and a bootstrap interval excluding 0, in every such cell.")

# THE FORMULATION, AND WHAT LICENSES IT. Fixed here so it cannot be reached for
# after seeing a number that only half supports it.
#   "Overwhelmed, not undone" is available ONLY if BOTH hold:
#     (a) the verdict is CONTENT and it is MATERIAL, and
#     (b) rho(est, ox_nuc_mtrib) is <= 0, or its interval covers 0, in every
#         adjudicating cell of BOTH cohorts.
#   If (a) holds and (b) does not - MYC raises the OXPHOS share against the
#   mitoribosome TOO - then the reprioritisation is being run backwards after
#   all, and this formulation is NOT available. That is the falsifier, and it is
#   written down before the panel is computed, which is what CLAUDE.md asks for.
FORMULATION_RULE <- paste0(
  "'overwhelmed, not undone' requires (a) CONTENT and MATERIAL, and (b) ",
  "rho(ox_nuc_mtrib) <= 0 or its CI covering 0, in every adjudicating cell of ",
  "both cohorts. (a) without (b) means the reprioritisation IS reversed and ",
  "the formulation is unavailable.")

message("\n0. rules fixed before any number")
message("   adjudication floor  |raw rho(est, ox_lvl)| >= ", ADJUDICATE_FLOOR,
        " in BOTH cohorts")
message("   materiality bar     |delta| >= ", MATERIAL_DELTA,
        "  (borrowed from E16's R1_MAX_DELTA)")
message("   contrast interval   paired sample-level bootstrap, B = ", BOOT_B,
        ", seed ", PROJECT_SEED)
message("   verdict rule        ", VERDICT_RULE)
message("   formulation rule    ", FORMULATION_RULE)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))

# THE DEPENDENCY THAT MAKES THIS A RE-READ. Without E16's saved separability
# table there is nothing to reproduce against, and this script would be a
# rebuild wearing a re-read's name.
if (!file.exists(PATH_E16)) {
  stop("results/respiratory_rulers.rds is missing. E20 is a RE-READ of E16 and ",
       "cannot run without it: section 3 reproduces E16's raw separability ",
       "cells exactly, and that check is what licenses the word 're-read'. ",
       "Run E16 first.", call. = FALSE)
}
e16 <- readRDS(PATH_E16)
if (!"separability" %in% names(e16)) {
  stop("results/respiratory_rulers.rds carries no `separability` table. ",
       "Re-run E16.", call. = FALSE)
}

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)

if (!identical(colnames(nw$tcga_gsva_new), ID_T) ||
    !identical(colnames(sc$gsva_new), ID_S) ||
    !identical(names(nw$tcga_log2MYC), ID_T) ||
    !identical(names(sc$log2MYC), ID_S)) {
  stop("the score objects are not aligned on the same samples as the arms. ",
       "Re-run E02.", call. = FALSE)
}

# --- 1.1 the three gene sets -------------------------------------------------
# All three come from the SNAPSHOT's own objects, so they are the identical
# sets the GSVA and mitoPPS arms were built on - not a re-parse of the workbook.
OX_SUB   <- mito$arm_sets[["OXPHOS subunits"]]        # numerator, all three rulers
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL              # the pinned inventory
REST     <- setdiff(MITO_ALL, OX_SUB)                 # ox_rel's denominator
MITORIB  <- mito$arm_sets[["Mitochondrial ribosome"]] # ox_nuc_mtrib's denominator

# The facts the new ruler rests on, asserted rather than described.
#   - the numerator carries no mtDNA-encoded gene, so nothing is pooled;
#   - the mitoribosome carries none either, so the new denominator does not
#     smuggle the 13 in through the back door;
#   - the two halves are DISJOINT, so neither is partly correlated with itself;
#   - the mitoribosome sits inside the same pinned inventory, so `ox_nuc_mtrib`
#     is a narrowing of `ox_rel`'s denominator and not a different universe.
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL,
          length(OX_SUB) == 89L, length(MITORIB) == 83L,
          !any(grepl("^MT-", OX_SUB)), !any(grepl("^MT-", MITORIB)),
          length(intersect(OX_SUB, MITORIB)) == 0L,
          all(MITORIB %in% MITO_ALL), all(MITORIB %in% REST))
message("   numerator    `OXPHOS subunits`       ", length(OX_SUB),
        " genes, 0 mtDNA-encoded")
message("   denominator1  MitoCarta minus those  ", length(REST),
        " genes  -> ox_rel        (E16's, unchanged)")
message("   denominator2 `Mitochondrial ribosome` ", length(MITORIB),
        " genes  -> ox_nuc_mtrib (NEW; a subset of denominator1)")

# --- 1.2 WHAT THE NEW RULER SHARES WITH WHAT IT IS MEASURED AGAINST ----------
# `ox_nuc_mtrib` is a difference of two composites, so an estimator overlapping
# EITHER half would be partly correlated with itself. Every `__MITOSTRIP`
# estimator is 0/0 by construction - MITOCARTA_ALL is the strip set and both
# halves live inside it - and this prints the number rather than asserting the
# argument. `PROLIF_DISJOINT` shares 10 genes with `ox_rel`'s 1,047-gene
# denominator (E16 section 1.2 measured that this moves nothing: Spearman
# 0.9999); against the mitoribosome the expected overlap is 0.
PROLIF_COV <- "PROLIF_DISJOINT"
PD <- mito$covariate_sets[[PROLIF_COV]]
EST_GENES <- list(
  FELSHER__MITOSTRIP                 = sd_$myc_sets[[MYC_REF]],
  MYC_UP.V1_UP__MITOSTRIP            = sd_$myc_sets[["MYC_UP.V1_UP__MITOSTRIP"]],
  HALLMARK_MYC_TARGETS_V1__MITOSTRIP = sd_$myc_sets[["HALLMARK_MYC_TARGETS_V1__MITOSTRIP"]],
  M_b__MITOSTRIP                     = sd_$collectri_sets[[MB_REF]],
  PROLIF_DISJOINT                    = PD)
overlap_audit <- tibble::tibble(
  set           = names(EST_GENES),
  n             = lengths(EST_GENES),
  x_numerator   = vapply(EST_GENES, function(s) length(intersect(s, OX_SUB)),
                         integer(1)),
  x_rest_mito   = vapply(EST_GENES, function(s) length(intersect(s, REST)),
                         integer(1)),
  x_mitoribosome = vapply(EST_GENES, function(s) length(intersect(s, MITORIB)),
                          integer(1)))
message("\n   gene overlap of each estimator with the halves of the rulers:")
overlap_audit %>% as.data.frame() %>% print(row.names = FALSE)
stopifnot(all(overlap_audit$x_mitoribosome == 0L))
message("   the mitoribosome column is 0 everywhere, INCLUDING the covariate:",
        "\n   the new ruler shares no gene with anything it is measured against.")

# --- 1.3 the gene matrices ---------------------------------------------------
# SCALE: log2(linear DESeq2-normalised + 1). E16's `L`, and the mouse's.
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

# =============================================================================
# 2. The rulers
# =============================================================================
# `.comp` is E16's, copied verbatim: the mean, per sample, of each gene's
# z-score across samples, on log2(linear + 1). That is the mouse's `comp_e`.
# The reproduction check in section 3 is what proves the copy has not drifted.
#
# NOTE ON THE RECIPE, so the departure is a decision and not a slip. The mouse's
# helper is `relify(s) = comp(s) - comp(mito_all \ s)` - a set against the REST
# OF THE COMPARTMENT. `ox_nuc_mtrib` is NOT `relify("MITOCARTA_MITOCHONDRIAL_
# RIBOSOME")`; it is a different object, and deliberately. What is matched to
# the mouse is the COMPOSITE recipe (`comp_e`), which is what "the same recipe"
# can mean here; the denominator is chosen to state the reprioritisation axis
# directly, OXPHOS against the mitoribosome, rather than each against the whole.
message("\n2. the rulers - one numerator, three denominators")

.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res)
  M  <- gr$mat
  v  <- apply(M, 1L, stats::var)
  # zero-variance rows would return NaN from scale() and poison the composite:
  # dropped, and COUNTED rather than silently tolerated.
  list(score   = colMeans(t(scale(t(M[v > 0, , drop = FALSE])))),
       n_asked = length(unique(genes)),
       n_used  = sum(v > 0),
       n_zerovar = sum(v == 0),
       missing = gr$missing)
}

rulers <- list(); set_sizes <- list()
for (coh in names(COH)) {
  C   <- COH[[coh]]
  num <- .comp(OX_SUB,  C)
  den <- .comp(REST,    C)
  mtr <- .comp(MITORIB, C)

  M <- rbind(
    ox_lvl       = num$score,                 # CONTENT
    ox_rel       = num$score - den$score,     # SHARE, vs the whole compartment
    ox_nuc_mtrib = num$score - mtr$score,     # SHARE, vs the mitoribosome (NEW)
    mtrib_lvl    = mtr$score,                 # context: the denominator alone
    ox_gsva      = as.numeric(C$arms$gsva_arms["OXPHOS subunits", C$ids]),
    ox_ppd       = as.numeric(C$arms$mitopps_arms["OXPHOS subunits", C$ids]))
  colnames(M) <- C$ids
  stopifnot(!anyNA(M), identical(rownames(M), c(PANEL_RULERS, CONTEXT_RULERS)))

  rulers[[coh]]    <- M
  set_sizes[[coh]] <- tibble::tibble(
    cohort = coh,
    numerator_used = num$n_used, numerator_coverage = num$n_used / num$n_asked,
    rest_used = den$n_used, rest_coverage = den$n_used / den$n_asked,
    mitorib_used = mtr$n_used, mitorib_coverage = mtr$n_used / mtr$n_asked,
    numerator_missing = paste(num$missing, collapse = ", "),
    mitorib_missing   = paste(mtr$missing, collapse = ", "),
    n_zerovar = num$n_zerovar + den$n_zerovar + mtr$n_zerovar)
}
set_sizes <- dplyr::bind_rows(set_sizes)
message("\n   set sizes actually used, after the symbol map:")
set_sizes %>%
  dplyr::select(-numerator_missing, -mitorib_missing) %>%
  dplyr::mutate(dplyr::across(dplyr::ends_with("coverage"), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
# Unharmonised, SCAN-B's numerator covers 0.775 instead of 0.989 (trap 7).
stopifnot(all(set_sizes$numerator_coverage >= 0.98),
          all(set_sizes$mitorib_coverage >= 0.95))

# --- 2.1 is the new ruler actually a new ruler? ------------------------------
# If `ox_nuc_mtrib` agrees with `ox_rel` at 0.99 it is `ox_rel` renamed and the
# panel has two columns, not three. Reported before anything is read off it.
ruler_agreement <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M  <- rulers[[coh]]
  cb <- utils::combn(rownames(M), 2L)
  tibble::tibble(
    cohort = coh, pair = paste(cb[1, ], "vs", cb[2, ]),
    rho = vapply(seq_len(ncol(cb)),
                 function(i) .rho(M[cb[1, i], ], M[cb[2, i], ]), numeric(1)))
}))
message("\n   ruler agreement, WITHIN cohort only (never across):")
ruler_agreement %>%
  dplyr::mutate(rho = round(rho, 3)) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 3. THE REPRODUCTION CHECK - what makes this a re-read
# =============================================================================
# `ox_lvl` and `ox_rel` are E16's own rulers, rebuilt here to E16's own recipe,
# and the raw Spearman panel in section 4 must therefore return E16's saved
# `$separability` cells BIT-EQUAL. `ox_gsva` and `ox_ppd` are read from the same
# saved objects E16 read and must reproduce too.
#
# This is the control for the whole script. If it fails, `.comp` has drifted or
# an input has moved, and every number below is a rebuild rather than a re-read.
message("\n3. reproduction check against E16")

EST_PRIMARY <- c(
  myc_mRNA    = "log2MYC",                              # the mouse's mRNA row
  myc_msigdb  = "HALLMARK_MYC_TARGETS_V1__MITOSTRIP",
  myc_regulon = MB_REF,                                 # CollecTRI
  myc_felsher = MYC_REF,                                # this repo's M_a
  myc_lowent  = "MYC_UP.V1_UP__MITOSTRIP",              # 1.5 pct entangled
  prolif      = PROLIF_COV)
EST_MYC <- setdiff(names(EST_PRIMARY), "prolif")

.estimator_matrix <- function(C, want) {
  rows <- lapply(want, function(w) {
    if (identical(w, "log2MYC"))     return(as.numeric(C$log2myc[C$ids]))
    if (identical(w, PROLIF_COV))    return(as.numeric(C$arms$gsva_cov[PROLIF_COV, C$ids]))
    if (w %in% rownames(C$mb))       return(as.numeric(C$mb[w, C$ids]))
    if (w %in% rownames(C$gsva_new)) return(as.numeric(C$gsva_new[w, C$ids]))
    stop("estimator not scored in this cohort: ", w, call. = FALSE)
  })
  M <- do.call(rbind, rows)
  dimnames(M) <- list(names(want), C$ids)
  M
}
ests <- lapply(COH, function(C) .estimator_matrix(C, EST_PRIMARY))

raw_rho <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M <- rulers[[coh]]; E <- ests[[coh]]
  tidyr::expand_grid(ruler = rownames(M), against = rownames(E)) %>%
    dplyr::mutate(cohort = coh,
                  rho = purrr::map2_dbl(ruler, against, ~ .rho(M[.x, ], E[.y, ])))
}))

repro <- raw_rho %>%
  dplyr::filter(ruler %in% c("ox_lvl", "ox_rel", "ox_gsva", "ox_ppd")) %>%
  dplyr::inner_join(
    e16$separability %>%
      dplyr::filter(panel == "primary") %>%
      dplyr::select(cohort, ruler, against, e16_rho = spearman),
    by = c("cohort", "ruler", "against")) %>%
  dplyr::mutate(delta = rho - e16_rho)
stopifnot(nrow(repro) == 2L * 4L * length(EST_PRIMARY))
if (max(abs(repro$delta)) >= 1e-10) {
  stop("E20 does not reproduce E16's separability cells (max |delta| = ",
       signif(max(abs(repro$delta)), 3), "). This is NOT a re-read: either ",
       "`.comp` has drifted from E16's copy or an input object has moved. ",
       "Stop and reconcile before writing anything.", call. = FALSE)
}
message("   ", nrow(repro), " of E16's raw separability cells reproduce exactly",
        " (max |delta| = ", signif(max(abs(repro$delta)), 3), ")")
message("   the two share rulers and the two read rulers all come back",
        " bit-equal,\n   so section 4's raw panel IS E16's table and the new",
        " ruler is the only\n   new number in this script.")

# =============================================================================
# 4. THE PANEL - rho with interval, raw and proliferation-adjusted
# =============================================================================
# `.cor_block` is E10's, copied verbatim as E16, E18 and E19 copy it: partial
# Spearman on ranks, Fisher-z intervals, Bonett-Wright variance
# se = sqrt((1 + rho^2/2)/(n - 3 - k)). The plain 1/(n-3) is the Pearson case
# and understates a rank correlation. Section 3 is what proves the copy is
# faithful.
#
# ALL FIVE MYC ESTIMATORS ARE REPORTED, plus proliferation, plus all six rulers.
# Selecting the pair that agrees after seeing them is the grid-of-cells trap.
#
# PROLIFERATION IS RAW ONLY. Partialling PROLIF_DISJOINT out of a correlation
# WITH PROLIF_DISJOINT is degenerate, so that cell is NA in the adjusted panel
# rather than printed as a near-zero that a reader might interpret.
message("\n4. the panel")

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

CV <- lapply(names(COH), function(coh)
  t(COH[[coh]]$arms$gsva_cov[PROLIF_COV, COH[[coh]]$ids, drop = FALSE]))
names(CV) <- names(COH)
stopifnot(identical(rownames(CV$TCGA), ID_T),
          identical(rownames(CV$`SCAN-B`), ID_S))

panel <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M <- rulers[[coh]]; E <- ests[[coh]]
  dplyr::bind_rows(
    .cor_block(M, E)              %>% dplyr::mutate(adjustment = "raw"),
    .cor_block(M, E, cov = CV[[coh]]) %>%
      dplyr::mutate(adjustment = "PROLIF_DISJOINT")) %>%
    dplyr::mutate(cohort = coh)
})) %>%
  dplyr::rename(against = gene) %>%
  # the degenerate cell, made NA rather than printed
  dplyr::mutate(dplyr::across(
    c(rho, ci_lo, ci_hi),
    ~ dplyr::if_else(against == "prolif" & adjustment == "PROLIF_DISJOINT",
                     NA_real_, .x))) %>%
  dplyr::mutate(
    role = dplyr::if_else(ruler %in% PANEL_RULERS, "panel", "context"),
    ruler = factor(ruler, levels = c(PANEL_RULERS, CONTEXT_RULERS)),
    against = factor(against, levels = names(EST_PRIMARY))) %>%
  dplyr::arrange(cohort, adjustment, ruler, against)

# The raw panel must still be E16's table, now with intervals attached.
chk <- panel %>%
  dplyr::filter(adjustment == "raw") %>%
  dplyr::mutate(ruler = as.character(ruler), against = as.character(against)) %>%
  dplyr::inner_join(raw_rho, by = c("cohort", "ruler", "against"),
                    suffix = c("", "_direct"))
stopifnot(nrow(chk) == 2L * 6L * length(EST_PRIMARY),
          max(abs(chk$rho - chk$rho_direct)) < 1e-12)

for (adj in c("raw", "PROLIF_DISJOINT")) {
  for (coh in names(COH)) {
    message("\n   ", adj, ", ", coh, ", Spearman. First three rows are the",
            " panel; the last three are context:")
    panel %>%
      dplyr::filter(adjustment == adj, cohort == coh) %>%
      dplyr::select(ruler, against, rho) %>%
      dplyr::mutate(rho = round(rho, 3)) %>%
      tidyr::pivot_wider(names_from = against, values_from = rho) %>%
      as.data.frame() %>% print(row.names = FALSE)
  }
}

message("\n   the three panel rulers with intervals, adjusted, MYC estimators",
        " only:")
panel %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT", role == "panel",
                against %in% EST_MYC) %>%
  dplyr::transmute(cohort, against, ruler, n,
                   rho = sprintf("%+.3f", rho),
                   ci  = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi)) %>%
  dplyr::arrange(against, cohort, ruler) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 4.1 which estimators adjudicate -----------------------------------------
# The floor from section 0, applied mechanically. Every estimator stays in the
# table; this decides only which ones the verdict is scored on.
adjudicate <- panel %>%
  dplyr::filter(adjustment == "raw", ruler == "ox_lvl", against %in% EST_MYC) %>%
  dplyr::select(cohort, against, rho) %>%
  dplyr::mutate(abs_rho = abs(rho)) %>%
  dplyr::group_by(against) %>%
  dplyr::summarise(min_abs_rho = min(abs_rho), .groups = "drop") %>%
  dplyr::mutate(adjudicates = min_abs_rho >= ADJUDICATE_FLOOR)
message("\n   which estimators adjudicate (|raw rho with ox_lvl| >= ",
        ADJUDICATE_FLOOR, " in BOTH cohorts):")
adjudicate %>%
  dplyr::mutate(min_abs_rho = round(min_abs_rho, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
EST_ADJ <- as.character(adjudicate$against[adjudicate$adjudicates])
message("   adjudicating: ", paste(EST_ADJ, collapse = ", "))
message("   excluded    : ",
        paste(setdiff(EST_MYC, EST_ADJ), collapse = ", "),
        "  (reported in every table, never scored)")

# =============================================================================
# 5. THE CONTRAST - share versus content, with a paired bootstrap
# =============================================================================
# Two differences per estimator per cohort:
#     content_minus_share  = rho(est, ox_lvl) - rho(est, ox_rel)
#     content_minus_mtrib  = rho(est, ox_lvl) - rho(est, ox_nuc_mtrib)
# The three rulers share their numerator, so these differences are differences
# of DEPENDENT correlations and the marginal Fisher-z intervals in section 4 say
# nothing useful about them. The interval below resamples SAMPLES and recomputes
# both rhos on the same resample, which is the paired quantity.
message("\n5. the contrast, with a paired bootstrap (B = ", BOOT_B, ")")

.boot_rho <- function(M, E, cov, B, seed) {
  n <- ncol(M)
  set.seed(seed)
  keep <- list(raw = array(NA_real_, c(B, nrow(M), nrow(E))),
               adj = array(NA_real_, c(B, nrow(M), nrow(E))))
  for (b in seq_len(B)) {
    i  <- sample.int(n, n, replace = TRUE)
    RA <- .rank_rows(M[, i, drop = FALSE])
    RB <- .rank_rows(E[, i, drop = FALSE])
    keep$raw[b, , ] <- suppressWarnings(stats::cor(t(RA), t(RB)))
    H  <- qr(cbind(`(Intercept)` = 1, apply(cov[i, , drop = FALSE], 2L, rank)))
    keep$adj[b, , ] <- suppressWarnings(
      stats::cor(t(RA - t(qr.fitted(H, t(RA)))),
                 t(RB - t(qr.fitted(H, t(RB))))))
  }
  for (nm in names(keep)) dimnames(keep[[nm]]) <- list(NULL, rownames(M), rownames(E))
  keep
}

message("   resampling ", BOOT_B, " times per cohort - the slow part of this",
        " script")
boots <- lapply(names(COH), function(coh)
  .boot_rho(rulers[[coh]], ests[[coh]], CV[[coh]], BOOT_B, PROJECT_SEED))
names(boots) <- names(COH)

.contrast <- function(coh, adj, est, a, b) {
  A <- boots[[coh]][[adj]][, a, est]
  B <- boots[[coh]][[adj]][, b, est]
  d <- A - B
  q <- stats::quantile(d, c(0.025, 0.975), names = FALSE)
  tibble::tibble(boot_lo = q[1], boot_hi = q[2], boot_sd = stats::sd(d))
}

contrast <- dplyr::bind_rows(lapply(names(COH), function(coh)
  dplyr::bind_rows(lapply(c("raw", "PROLIF_DISJOINT"), function(adj) {
    key <- if (adj == "raw") "raw" else "adj"
    P <- panel %>%
      dplyr::filter(cohort == coh, adjustment == adj) %>%
      dplyr::mutate(ruler = as.character(ruler), against = as.character(against))
    .get <- function(r, a) P$rho[P$ruler == r & P$against == a]
    dplyr::bind_rows(lapply(names(EST_PRIMARY), function(a) {
      if (a == "prolif" && adj != "raw") return(NULL)
      dplyr::bind_rows(
        tibble::tibble(contrast = "ox_lvl - ox_rel",
                       delta = .get("ox_lvl", a) - .get("ox_rel", a)) %>%
          dplyr::bind_cols(.contrast(coh, key, a, "ox_lvl", "ox_rel")),
        tibble::tibble(contrast = "ox_lvl - ox_nuc_mtrib",
                       delta = .get("ox_lvl", a) - .get("ox_nuc_mtrib", a)) %>%
          dplyr::bind_cols(.contrast(coh, key, a, "ox_lvl", "ox_nuc_mtrib")),
        # CONTEXT, NOT A VERDICT INPUT. `ox_nuc_mtrib` is a difference of two
        # composites, so a reader cannot tell from its sign alone whether MYC
        # is failing to raise OXPHOS or is raising the mitoribosome harder.
        # This row is that decomposition and nothing else; `role` keeps it out
        # of section 6.
        tibble::tibble(contrast = "mtrib_lvl - ox_lvl",
                       delta = .get("mtrib_lvl", a) - .get("ox_lvl", a)) %>%
          dplyr::bind_cols(.contrast(coh, key, a, "mtrib_lvl", "ox_lvl"))) %>%
        dplyr::mutate(cohort = coh, adjustment = adj, against = a, .before = 1L)
    }))
  })))) %>%
  dplyr::mutate(
    role          = dplyr::if_else(contrast == "mtrib_lvl - ox_lvl",
                                   "context", "verdict"),
    excludes_zero = (boot_lo > 0) | (boot_hi < 0),
    material      = excludes_zero & abs(delta) >= MATERIAL_DELTA,
    adjudicates   = against %in% EST_ADJ)

.print_contrast <- function(d) {
  d %>%
    dplyr::transmute(cohort, against, adj = adjudicates,
                     delta = sprintf("%+.3f", delta),
                     boot = sprintf("[%+.3f, %+.3f]", boot_lo, boot_hi),
                     material) %>%
    dplyr::arrange(against, cohort) %>%
    as.data.frame() %>% print(row.names = FALSE)
}
for (adj in c("raw", "PROLIF_DISJOINT")) {
  for (ct in c("ox_lvl - ox_rel", "ox_lvl - ox_nuc_mtrib",
               "mtrib_lvl - ox_lvl")) {
    message("\n   ", adj, "  |  ", ct,
            if (ct == "mtrib_lvl - ox_lvl") "   (CONTEXT, not a verdict input)"
            else "")
    .print_contrast(contrast %>%
                      dplyr::filter(adjustment == adj, contrast == ct))
  }
}

# =============================================================================
# 6. THE VERDICT, on the rules fixed in section 0
# =============================================================================
message("\n6. verdict")

VP <- panel %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT", role == "panel",
                against %in% EST_ADJ) %>%
  dplyr::mutate(ruler = as.character(ruler), against = as.character(against)) %>%
  dplyr::select(cohort, against, ruler, rho, ci_lo, ci_hi) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = c(rho, ci_lo, ci_hi))

VP <- VP %>%
  dplyr::mutate(
    content_wins = rho_ox_lvl > rho_ox_rel & rho_ox_lvl > rho_ox_nuc_mtrib,
    share_wins   = rho_ox_lvl < rho_ox_rel & rho_ox_lvl < rho_ox_nuc_mtrib)
message("\n   adjudicating cells, adjusted panel:")
VP %>%
  dplyr::transmute(cohort, against,
                   ox_lvl = sprintf("%+.3f", rho_ox_lvl),
                   ox_rel = sprintf("%+.3f", rho_ox_rel),
                   ox_nuc_mtrib = sprintf("%+.3f", rho_ox_nuc_mtrib),
                   content_wins, share_wins) %>%
  as.data.frame() %>% print(row.names = FALSE)

verdict <- if (all(VP$content_wins)) "CONTENT" else
  if (all(VP$share_wins)) "SHARE" else "ESTIMATOR-DEPENDENT"
dissent <- VP %>%
  dplyr::filter(!(if (verdict == "SHARE") share_wins else content_wins)) %>%
  dplyr::transmute(cell = paste0(against, " / ", cohort)) %>%
  dplyr::pull(cell)

material_all <- contrast %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT", adjudicates,
                role == "verdict") %>%
  dplyr::summarise(all_material = all(material),
                   min_abs_delta = min(abs(delta))) 

# (b) of the formulation rule: does the mitoribosome-relative ruler stay at or
# below zero? Its interval is the marginal Fisher-z one, which is the right
# object here - this is a statement about ONE correlation, not about a
# difference of two.
share_not_raised <- VP %>%
  dplyr::mutate(ok = rho_ox_nuc_mtrib <= 0 |
                  (ci_lo_ox_nuc_mtrib <= 0 & ci_hi_ox_nuc_mtrib >= 0))
formulation_ok <- verdict == "CONTENT" && isTRUE(material_all$all_material) &&
  all(share_not_raised$ok)

message("\n   VERDICT: ", verdict, "  (", nrow(VP), " adjudicating cells, ",
        length(EST_ADJ), " estimators x 2 cohorts)")
if (length(dissent))
  message("   cells that dissent: ", paste(dissent, collapse = "; "))
message("   MATERIAL in every adjudicating cell: ", material_all$all_material,
        "  (smallest |delta| ", sprintf("%.3f", material_all$min_abs_delta),
        " against a bar of ", MATERIAL_DELTA, ")")
message("   rho(ox_nuc_mtrib) <= 0 or interval covering 0, everywhere: ",
        all(share_not_raised$ok))
message("\n   'OVERWHELMED, NOT UNDONE' IS ",
        if (formulation_ok) "AVAILABLE" else "NOT AVAILABLE",
        " on these numbers.")
if (!formulation_ok) {
  message("   Report which estimators disagree and stop. Do not adjudicate a",
          " mixed\n   result with a tie-break invented afterwards - that is",
          " what section 0\n   exists to prevent.")
}

verdicts <- tibble::tibble(
  item = c("verdict", "adjudicating estimators", "dissenting cells",
           "material everywhere", "share not raised", "formulation available"),
  value = c(verdict, paste(EST_ADJ, collapse = ", "),
            if (length(dissent)) paste(dissent, collapse = "; ") else "none",
            as.character(material_all$all_material),
            as.character(all(share_not_raised$ok)),
            as.character(formulation_ok)))
message("")
verdicts %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. Save
# =============================================================================
saveRDS(list(
  set_sizes = set_sizes, overlap_audit = overlap_audit,
  ruler_agreement = ruler_agreement,
  repro = repro,
  panel = panel, adjudicate = adjudicate,
  contrast = contrast,
  verdict_panel = VP, verdicts = verdicts,
  verdict = verdict, dissent = dissent,
  formulation_available = formulation_ok,
  settings = list(
    what_this_is = paste("a RE-READ of E16's separability table for a different",
                         "question, plus ONE new ruler (ox_nuc_mtrib). Nothing",
                         "is re-scored; mitoPPS and GSVA are read."),
    rulers = c(ox_lvl = "comp(OXPHOS subunits) - CONTENT",
               ox_rel = "comp(OXPHOS subunits) - comp(MitoCarta rest) - SHARE",
               ox_nuc_mtrib = "comp(OXPHOS subunits) - comp(mitoribosome) - NEW",
               mtrib_lvl = "comp(mitoribosome) - context",
               ox_gsva = "GSVA on VST, read - the incumbent",
               ox_ppd = "mitoPPS, read"),
    panel_rulers = PANEL_RULERS, context_rulers = CONTEXT_RULERS,
    composite_scale = "mean per-gene z of log2(linear DESeq2-normalised + 1)",
    recipe_note = paste("comp_e matched to myc_mouse script 48 PART A. This is",
                        "NOT relify(): relify's denominator is the rest of the",
                        "compartment, ox_nuc_mtrib's is the mitoribosome alone,",
                        "chosen to state the reprioritisation axis directly."),
    myc_estimators = EST_PRIMARY,
    covariate = PROLIF_COV,
    measure = "spearman; panel reported raw AND partial; verdict on partial",
    decision_rules = c(ADJUDICATE_FLOOR = ADJUDICATE_FLOOR,
                       MATERIAL_DELTA = MATERIAL_DELTA, BOOT_B = BOOT_B),
    verdict_rule = VERDICT_RULE,
    formulation_rule = FORMULATION_RULE,
    seed = PROJECT_SEED),
  rules = list(
    scope = paste("RE-READ. No interaction model, no MYC stratification, no",
                  "subtype split, no new endpoint or ratio, no re-scoring."),
    n3 = "these are transcript correlations; the word 'primed' is never used",
    cohorts = "never pooled and never compared as values; a species is a cohort",
    trap9 = paste("separability is a property of the ruler-ESTIMATOR pair",
                  "(E16 3.1); every claim names its estimator"),
    mouse = "the mouse is not read, compared or cited as a value by this script"),
  built = Sys.time()), PATH_E20)

readr::write_csv(
  panel %>% dplyr::mutate(ruler = as.character(ruler),
                          against = as.character(against)),
  PATH_E20_PANEL)

message("\nE20: done.")
message("    results/share_vs_content.rds")
message("    outputs/tables/E20_share_vs_content.csv")
message("    no figures - this is a re-read, and every number in it is a table.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E20)

  # The answer, and the rules it was scored on.
  x$verdicts %>% as.data.frame()
  x$settings$verdict_rule
  x$settings$formulation_rule

  # THE CONTROL. If this is not 48 cells at max |delta| = 0, the script is a
  # rebuild rather than a re-read and nothing else should be read.
  nrow(x$repro); summary(abs(x$repro$delta))

  # Is the new ruler a new ruler, or ox_rel renamed?
  x$ruler_agreement %>%
    dplyr::mutate(rho = round(rho, 3)) %>%
    tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
    as.data.frame()

  # The panel, adjusted, all six rulers.
  x$panel %>%
    dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
    dplyr::select(cohort, ruler, against, rho) %>%
    dplyr::mutate(rho = round(rho, 3)) %>%
    tidyr::pivot_wider(names_from = against, values_from = rho) %>%
    as.data.frame()

  # ...and raw, which is E16's table with intervals attached.
  x$panel %>%
    dplyr::filter(adjustment == "raw") %>%
    dplyr::select(cohort, ruler, against, rho) %>%
    dplyr::mutate(rho = round(rho, 3)) %>%
    tidyr::pivot_wider(names_from = against, values_from = rho) %>%
    as.data.frame()

  # Which estimators were allowed to adjudicate, and why.
  x$adjudicate %>% as.data.frame()

  # The contrast. Read the bootstrap interval, NEVER the two marginal ones -
  # the rulers share a numerator and correlate above 0.9.
  x$contrast %>%
    dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
    dplyr::transmute(cohort, against, contrast, adjudicates,
                     delta = round(delta, 3),
                     boot = sprintf("[%+.3f, %+.3f]", boot_lo, boot_hi),
                     material) %>%
    as.data.frame()

  # The context rulers: is the mitoribosome really the MYC-nearer arm?
  x$panel %>%
    dplyr::filter(adjustment == "PROLIF_DISJOINT",
                  ruler %in% c("ox_lvl", "mtrib_lvl")) %>%
    dplyr::select(cohort, ruler, against, rho) %>%
    dplyr::mutate(rho = round(rho, 3)) %>%
    tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
    as.data.frame()

}
