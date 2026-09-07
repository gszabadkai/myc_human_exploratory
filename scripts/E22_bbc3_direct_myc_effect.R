# E22_bbc3_direct_myc_effect.R
# =============================================================================
# IS THERE A NEGATIVE DIRECT MYC EFFECT ON BBC3 AT FIXED OXPHOS?
#
# STAGE 1 ONLY - THE UNADJUSTED TWO-PREDICTOR FIT. THIS SCRIPT STOPS AT GATE 2.
# The adjusted fits (purity, PROLIF_DISJOINT), the within-stratum fits and the
# gated FOXO3 mediation are appended to THIS FILE after the gate clears. Do not
# add them before it does.
#
# =============================================================================
# WHAT THIS IS NOT - verbatim from the prompt, 2026-09-07
# =============================================================================
# This is a MAIN-EFFECTS model. It is a different estimand from the
# pre-registered MYC:OXPHOS interaction on PRIME, which is null and permanently
# closed at d3ac60e.
#
# - Do NOT include an M:O product term anywhere.
# - Do NOT compute PRIME.
# - Do NOT read, reference, or import anything from myc_human_validation.
# - Nothing produced here may be written as bearing on the registered
#   interaction result.
#
# =============================================================================
# PRE-SPECIFICATION - written before any coefficient in this script existed
# =============================================================================
# Primary quantities, per cohort, on Blom normal-scored variables:
#
#     beta1 = partial MYC    -> BBC3 at fixed OXPHOS    PREDICTED NEGATIVE
#     beta2 = partial OXPHOS -> BBC3 at fixed MYC       PREDICTED POSITIVE
#
# FALSIFICATION: beta1 >= 0, or an interval spanning zero with a point estimate
# below the magnitude needed to reconcile the observed marginal, means the
# FOXO3-repression arm does not transfer to human tumours. That is a reportable
# negative and is written as one.
#
# THE MOTIVATION. The human arm has reported that BBC3 tracks OXPHOS. Under an
# additive model, beta2 > 0 together with rho(MYC, OXPHOS) > 0 predicts
# rho(BBC3, MYC) > 0. The observed marginal is not positive, which requires a
# NEGATIVE direct MYC effect at fixed OXPHOS. The mouse gland predicts that sign
# (MYC represses Foxo3, which carries PUMA). This tests whether human delivers it.
#
# =============================================================================
# FOUR THINGS THE GATE 1 AUDIT SETTLED, AND THEY CHANGE HOW THIS IS READ
# =============================================================================
# (1) THE MARGINAL IS NOT ZERO, IT IS NEGATIVE - in 6 of 6 cohort-by-estimator
#     cells (TCGA -0.048 / -0.117 / -0.105; SCAN-B -0.057 / -0.004 / -0.158, for
#     dose / FELSHER / M_b). That STRENGTHENS the prediction: beta1 must be MORE
#     negative than reconciling to zero would need. The reconciliation target is
#     the observed negative value and never zero.
#
#     ANY MANUSCRIPT SENTENCE SAYING BBC3 IS "UNCORRELATED WITH MYC" IN HUMAN
#     TUMOURS IS WRONG AS WRITTEN and must read "weakly negatively correlated".
#     Correcting it is a deliverable of this script, not a footnote.
#
# (2) THE RECONCILIATION CHECK IS ALGEBRA, NOT EVIDENCE. For a two-predictor
#     OLS on variables of equal variance,
#         beta1 + beta2 * r(M,O) == r(Y,M)
#     identically, for any data whatsoever. Section 6 asserts it to 1e-10 and
#     reports it ONCE, labelled as ACCOUNTING - what it delivers is HOW NEGATIVE
#     beta1 has to be, not whether the model fits. The empirical weight sits on
#     beta1's SIGN, its replication across cohorts, estimators, strata and
#     rulers, and its specificity against the other eleven genes.
#
#     WHAT THE IDENTITY DOES NOT COVER, and it is the substantive test: beta1 is
#     NOT algebraically guaranteed to stay negative once purity,
#     PROLIF_DISJOINT and subtype enter. That is stage 2 and it is where the
#     result is actually decided.
#
# (3) THE DOSE ESTIMATOR'S RECONCILIATION ROW IS VACUOUS. The argument needs
#     rho(M, OXPHOS) > 0. For log2(MYC) it is -0.103 / -0.048 on ox_rel and
#     -0.016 / -0.003 on ox_lvl - at or below zero. With r(M,O) ~ 0 there is
#     nothing to reconcile: beta1 simply equals the marginal. That is CLAUDE.md
#     trap 4. Dose still carries the D2 sign rule; its identity row is flagged
#     VACUOUS and must not be read as agreement.
#
#     Because that fact is load-bearing well beyond this script - the collider
#     argument and the "MYC-high tumours are OXPHOS-high" premise both assume
#     rho(M, OXPHOS) > 0 - section 4 reports the FULL rho(M, OXPHOS) matrix in
#     its own right: 5 estimators x 2 rulers x 2 cohorts x 3 strata, with
#     intervals. It is reported and NOT chased here.
#
# (4) ELEVEN OF THE TWELVE OUTCOME GENES ARE INSIDE ox_rel's DENOMINATOR.
#     BAD, BBC3, BCL2, BCL2A1, BCL2L1, BCL2L11, BCL2L2, BID, BIK, MCL1 and
#     PMAIP1 are MitoCarta members; only BMF is not. So the outcome contributes
#     to its own predictor with weight -1/1047. Expected bias is tiny and it is
#     the OUTCOME rather than a covariate, so section 7 MEASURES it with a
#     leave-one-out rebuild of ox_rel per gene rather than arguing it away.
#     ox_lvl needs no leave-one-out: none of the twelve is an OXPHOS subunit,
#     and that is asserted.
#
# =============================================================================
# THE ESTIMATORS, AND WHY M_b SURVIVES THE CONTAMINATION CHECK
# =============================================================================
# BBC3 IS a CollecTRI MYC target. It is present in M_b__FULL and in
# M_b__PROLIFSTRIP, and it is REMOVED by the mitochondrial strip because BBC3 is
# a MitoCarta gene. The same holds for BCL2, BCL2L1, BCL2L11 and PMAIP1.
#
#   -> M_b__MITOSTRIP is structurally CLEAN for all twelve outcomes. Kept.
#   -> M_b__FULL and M_b__PROLIFSTRIP are DISQUALIFIED as the exposure for these
#      outcomes and are never fitted here. Asserted in section 1.
#
# D2 IS STILL OPEN, so three MYC estimators carry the claim and the rule is
# declared NOW, before any fit: the claim requires beta1 < 0 on ALL THREE, in
# BOTH cohorts. SIGN ONLY. M_b carries the widest standard error and is NOT
# dropped on effect size - the power table that would justify that was built for
# an interaction term, and a main effect needs roughly a quarter of the n.
#
# Two further signatures are carried as sensitivities and are NOT part of the
# rule: MYC_UP.V1_UP (1.5 pct proliferation-entangled, the clean one) and
# HALLMARK_MYC_TARGETS_V1 (23.5 pct, the entangled one).
#
# =============================================================================
# PURITY PROVENANCE - settled within this repo, as instructed
# =============================================================================
# No ingest or download script in this repository creates the purity column.
# `E01` line 202 reads `tcga_cov$purity` directly from the snapshot
# `data/from_validation/tcga_brca_covariates.rds`, and that directory's README
# lists the column without naming an algorithm.
#
#   -> IN THE NOTE AND IN METHODS: "purity, source undocumented in this repo."
#
# The column sits beside `ploidy` and `genome_doublings`, which is the triple
# ABSOLUTE emits. That is a provenance FINGERPRINT and is recorded as inference,
# never as documentation. `genome_doublings` is NOT a covariate here and is not
# added. SCAN-B has no purity estimate at all (0 of 3,207) and nothing is
# substituted for it - stage 2 bounds the SCAN-B fit using the TCGA
# with-minus-without delta instead.
#
# =============================================================================
# SCALES, STATED ONCE
# =============================================================================
#   the twelve genes   log2(linear DESeq2-normalised + 1), as E10, E16 and E20
#   ox_rel, ox_lvl     BUILT here from the same matrix, E16/E20's `.comp`
#   log2MYC            READ from the linear matrix
#   GSVA signatures    READ: GSVA on the VST, kcdf Gaussian
#   M_b__MITOSTRIP     READ: decoupleR ULM on the VST
#
# EVERY variable entering a model is then BLOM NORMAL-SCORED,
# qnorm((rank - 3/8)/(n + 1/4)), AND THEN SCALED TO UNIT VARIANCE.
#
# The second step is not cosmetic and the reason is worth recording. With no
# ties a Blom-scored vector is a permutation of one fixed vector, so every
# variable shares a variance and the section 6 identity is exact as it stands.
# TIES BREAK THAT: averaged ranks repeat values and shrink the variance a
# little, and SCAN-B has ties. Unstandardised, the variances there spread by
# 4.4e-4 and the identity would hold only to about 1e-4. Scaling to unit
# variance restores it to machine precision and changes nothing else - a
# two-predictor OLS on unit-variance variables has beta1 = (r_YM - r_YO r_MO) /
# (1 - r_MO^2) identically, which is what the pre-specification asks for.
# Section 3 counts the ties so their size is visible rather than absorbed.
#
# Blom-Pearson correlations are close to but not identical to the Spearman
# values in E10, E16 and E20. BOTH are reported side by side so this note cannot
# look like it contradicts those.
#
# BLOM SCORING IS POOLED, ONCE PER COHORT. Stage 2's within-stratum fits will
# SUBSET those pooled scores and will NOT re-score - E19's rule. One consequence
# is fixed here in advance: subsetting breaks the equal-variance property, so
# THE IDENTITY IS CLAIMED FOR THE POOLED FIT ONLY and is not computed within
# strata.
#
# NEVER POOLED ACROSS COHORTS - fitted separately, always. The mouse is not read,
# compared or cited by this script. N3 THROUGHOUT: these are transcript
# associations, and the word "primed" is never written of a transcript.
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))
source(here::here("functions", "strata.R"))

message("\nE22: a direct MYC effect on BBC3 at fixed OXPHOS? STAGE 1, to GATE 2\n",
        strrep("=", 78))

PATH_E22       <- file.path(DIR_RESULTS, "bbc3_direct_myc_effect.rds")
PATH_E22_FITA  <- file.path(DIR_TABLES,  "E22_fit_a.csv")
PATH_E22_RMO   <- file.path(DIR_TABLES,  "E22_rho_myc_oxphos.csv")

# =============================================================================
# 0. RULES FIXED BEFORE ANY COEFFICIENT EXISTS
# =============================================================================
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
PRIMING_ALL  <- sort(unique(c(PRIMING_PRO, PRIMING_ANTI)))
FOCUS_GENE   <- "BBC3"

# The three that carry the claim, and the two that do not.
EST_CLAIM <- c(dose      = "log2MYC",
               signature = "FELSHER__MITOSTRIP",
               regulon   = "M_b__MITOSTRIP")
EST_SENS  <- c(sig_clean     = "MYC_UP.V1_UP__MITOSTRIP",
               sig_entangled = "HALLMARK_MYC_TARGETS_V1__MITOSTRIP")
EST_ALL   <- c(EST_CLAIM, EST_SENS)

# Disqualified for THESE outcomes: they contain BBC3, BCL2, BCL2L1, BCL2L11 and
# PMAIP1. Named so the exclusion is a decision and section 1 can assert it.
EST_DISQUALIFIED <- c("M_b__FULL", "M_b__PROLIFSTRIP")

RULERS <- c("ox_rel", "ox_lvl")   # ox_rel primary; ox_lvl the ruler sensitivity
STRATA <- c("all", "Luminal", "Basal")

# THE D2 CONCORDANCE RULE, declared before fitting. SIGN ONLY.
D2_RULE <- paste0(
  "the claim requires beta1 < 0 on ALL THREE claim estimators (dose, signature, ",
  "regulon) in BOTH cohorts. SIGN ONLY - M_b carries the widest SE and is not ",
  "dropped on effect size. If they disagree, all three are reported and the ",
  "claim is UNSUPPORTED. No estimator is picked after the fact.")

# What "reconciled" would mean, so section 6 cannot drift into claiming more.
IDENTITY_NOTE <- paste0(
  "ACCOUNTING, NOT EVIDENCE. beta1 + beta2 * r(M,O) == r(Y,M) identically for a ",
  "two-predictor OLS on equal-variance variables. It cannot fail. What it ",
  "delivers is how negative beta1 must be to be consistent with the observed ",
  "marginal. A row with |r(M,O)| < ", 0.05, " is flagged VACUOUS.")
VACUOUS_RMO <- 0.05

# Basal is small and that is stated up front rather than discovered later.
BASAL_NOTE <- paste0(
  "TCGA Basal n = 171 (162 with purity). A 95% interval on rho there is about ",
  "+/- 0.15 wide, so a Basal null is WEAK EVIDENCE and is written as such, ",
  "never as failure to hold. The interval width is quoted beside every Basal ",
  "estimate.")

message("\n0. rules fixed before any coefficient")
message("   D2 rule      ", D2_RULE)
message("   identity     ", IDENTITY_NOTE)
message("   Basal        ", BASAL_NOTE)

# =============================================================================
# 1. Inputs, and the contamination audit as assertions
# =============================================================================
message("\n1. inputs")

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
fr   <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)

OX_SUB   <- mito$arm_sets[["OXPHOS subunits"]]
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL
REST     <- setdiff(MITO_ALL, OX_SUB)
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL, length(OX_SUB) == 89L,
          !any(grepl("^MT-", OX_SUB)))

# --- 1.1 what contaminates what ---------------------------------------------
EST_SETS <- c(
  lapply(EST_ALL[EST_ALL != "log2MYC"], function(nm)
    if (nm %in% names(sd_$myc_sets)) sd_$myc_sets[[nm]] else sd_$collectri_sets[[nm]]),
  lapply(EST_DISQUALIFIED, function(nm) sd_$collectri_sets[[nm]]),
  list(`ox_rel numerator` = OX_SUB, `ox_rel denominator` = REST))
names(EST_SETS) <- c(EST_ALL[EST_ALL != "log2MYC"], EST_DISQUALIFIED,
                     "ox_rel numerator", "ox_rel denominator")
contamination <- tibble::tibble(
  predictor = names(EST_SETS), n = lengths(EST_SETS),
  n_of_twelve = vapply(EST_SETS, function(s) length(intersect(PRIMING_ALL, s)),
                       integer(1)),
  which = vapply(EST_SETS, function(s)
    paste(intersect(PRIMING_ALL, s), collapse = ", "), character(1)))
message("\n   which of the twelve outcomes appear inside each predictor:")
contamination %>% as.data.frame() %>% print(row.names = FALSE)

# The three facts the design rests on, asserted rather than described.
.set_of <- function(nm) EST_SETS[[nm]]
stopifnot(
  # every FITTED estimator is clean of all twelve
  all(vapply(EST_ALL[EST_ALL != "log2MYC"],
             function(nm) length(intersect(PRIMING_ALL, .set_of(nm))) == 0L,
             logical(1))),
  # the two disqualified ones are NOT clean - which is why they are excluded
  all(vapply(EST_DISQUALIFIED,
             function(nm) FOCUS_GENE %in% .set_of(nm), logical(1))),
  # ox_lvl needs no leave-one-out
  length(intersect(PRIMING_ALL, OX_SUB)) == 0L)
LOO_NEEDED <- intersect(PRIMING_ALL, REST)
message("   FITTED estimators are clean of all twelve; ",
        paste(EST_DISQUALIFIED, collapse = " and "),
        " contain ", FOCUS_GENE, " and are EXCLUDED.")
message("   ox_lvl needs no leave-one-out (0 of the twelve are OXPHOS subunits);",
        "\n   ox_rel needs it for ", length(LOO_NEEDED), " of ",
        length(PRIMING_ALL), ": ", paste(LOO_NEEDED, collapse = ", "))

# --- 1.2 matrices and strata -------------------------------------------------
tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(tcga_lin$scale, "linear_deseq2_normalised"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))
GT <- log2(tcga_lin$mat[, ID_T, drop = FALSE] + 1)
GS <- log2(scanb_lin$mat[, ID_S, drop = FALSE] + 1)
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

COH <- list(
  TCGA     = list(G = GT, res = .symbol_resolver(rownames(GT), NULL), ids = ID_T,
                  gsva_new = nw$tcga_gsva_new, mb = nw$tcga_M_b_variants,
                  log2myc = nw$tcga_log2MYC),
  `SCAN-B` = list(G = GS, res = .symbol_resolver(rownames(GS), sc$symbol_map),
                  ids = ID_S, gsva_new = sc$gsva_new, mb = sc$M_b_variants,
                  log2myc = sc$log2MYC))

STR <- list(TCGA = .build_strata(fr, "TCGA", ID_T),
            `SCAN-B` = .build_strata(fr, "SCAN-B", ID_S))
strata_n <- tibble::tibble(stratum = STRATA,
                           TCGA = lengths(STR$TCGA[STRATA]),
                           `SCAN-B` = lengths(STR$`SCAN-B`[STRATA]))
message("\n   strata:")
strata_n %>% as.data.frame() %>% print(row.names = FALSE)

GR <- lapply(COH, function(C) .gene_rows(PRIMING_ALL, C$G, C$res))
miss <- unique(unlist(lapply(GR, function(g) g$missing), use.names = FALSE))
if (length(miss)) {
  stop("outcome gene(s) did not resolve to exactly one matrix row: ",
       paste(miss, collapse = ", "), call. = FALSE)
}
message("   all ", length(PRIMING_ALL), " outcome genes present in both cohorts")

# =============================================================================
# 2. The rulers, standard and leave-one-out
# =============================================================================
message("\n2. the rulers")

.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res)
  M  <- gr$mat
  v  <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}

rulers <- list(); rulers_loo <- list()
for (coh in names(COH)) {
  C <- COH[[coh]]
  num <- .comp(OX_SUB, C)
  M <- rbind(ox_rel = num - .comp(REST, C), ox_lvl = num)
  colnames(M) <- C$ids
  stopifnot(!anyNA(M))
  rulers[[coh]] <- M
  # ONE REBUILD PER OUTCOME GENE, for ox_rel only. The gene is dropped from the
  # DENOMINATOR; the numerator is untouched because no outcome is in it.
  L <- do.call(rbind, lapply(PRIMING_ALL, function(g)
    num - .comp(setdiff(REST, g), C)))
  dimnames(L) <- list(PRIMING_ALL, C$ids)
  rulers_loo[[coh]] <- L
}
message("   built ox_rel and ox_lvl, plus ", length(PRIMING_ALL),
        " leave-one-out ox_rel per cohort")

# =============================================================================
# 3. Blom normal scores
# =============================================================================
# qnorm((rank - 3/8)/(n + 1/4)). Every scored variable is then a PERMUTATION OF
# ONE FIXED VECTOR, so all of them share a variance exactly - which is what makes
# section 6's identity hold without further standardisation. Asserted below
# rather than assumed.
message("\n3. Blom normal scores")

# Blom normal score, THEN scaled to unit variance - see the header. The scaling
# matters only where there are ties, and it is applied unconditionally so that
# no cohort is treated differently from another.
.blom <- function(x) {
  n <- length(x)
  stopifnot(!anyNA(x))
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.375) / (n + 0.25))
  as.numeric(scale(z))
}
# How much tying there is, per variable. Reported, not assumed away: a variable
# with heavy ties is a variable whose normal score is not really continuous.
.tie_frac <- function(x) 1 - length(unique(x)) / length(x)

.estimator_vec <- function(C, w) {
  if (identical(w, "log2MYC"))     return(as.numeric(C$log2myc[C$ids]))
  if (w %in% rownames(C$mb))       return(as.numeric(C$mb[w, C$ids]))
  if (w %in% rownames(C$gsva_new)) return(as.numeric(C$gsva_new[w, C$ids]))
  stop("estimator not scored in this cohort: ", w, call. = FALSE)
}

B <- list(); ties <- list()
for (coh in names(COH)) {
  C <- COH[[coh]]
  E <- do.call(rbind, lapply(EST_ALL, function(w) .blom(.estimator_vec(C, w))))
  rownames(E) <- names(EST_ALL)
  Y <- do.call(rbind, lapply(PRIMING_ALL,
                             function(g) .blom(as.numeric(GR[[coh]]$mat[g, ]))))
  rownames(Y) <- PRIMING_ALL
  O <- do.call(rbind, lapply(RULERS, function(r) .blom(rulers[[coh]][r, ])))
  rownames(O) <- RULERS
  OL <- do.call(rbind, lapply(PRIMING_ALL,
                              function(g) .blom(rulers_loo[[coh]][g, ])))
  rownames(OL) <- PRIMING_ALL
  for (M in list(E, Y, O, OL)) colnames(M) <- C$ids
  colnames(E) <- colnames(Y) <- colnames(O) <- colnames(OL) <- C$ids
  B[[coh]] <- list(E = E, Y = Y, O = O, OL = OL)
  # THE PROPERTY THE IDENTITY NEEDS: unit variance on every scored vector.
  vs <- c(apply(E, 1L, stats::var), apply(Y, 1L, stats::var),
          apply(O, 1L, stats::var), apply(OL, 1L, stats::var))
  if (max(abs(vs - 1)) > 1e-10) {
    stop("Blom scores are not unit variance in ", coh, " (max |v - 1| = ",
         signif(max(abs(vs - 1)), 3), "). The identity in section 6 would not ",
         "hold; stop and reconcile.", call. = FALSE)
  }
  ties[[coh]] <- tibble::tibble(
    cohort = coh,
    variable = c(rownames(E), rownames(Y), "ox_rel", "ox_lvl"),
    kind = c(rep("estimator", nrow(E)), rep("outcome", nrow(Y)),
             "ruler", "ruler"),
    tie_frac = c(vapply(EST_ALL, function(w) .tie_frac(.estimator_vec(C, w)),
                        numeric(1)),
                 vapply(PRIMING_ALL, function(g)
                   .tie_frac(as.numeric(GR[[coh]]$mat[g, ])), numeric(1)),
                 .tie_frac(rulers[[coh]]["ox_rel", ]),
                 .tie_frac(rulers[[coh]]["ox_lvl", ])))
  message("   ", coh, ": ", nrow(E), " estimators, ", nrow(Y), " outcomes, ",
          nrow(O) + nrow(OL), " ruler vectors; all unit variance to 1e-10")
}

ties <- dplyr::bind_rows(ties)
tied <- ties %>% dplyr::filter(tie_frac > 0) %>% dplyr::arrange(dplyr::desc(tie_frac))
if (nrow(tied)) {
  message("\n   variables carrying ties (the reason for the unit-variance step):")
  tied %>% dplyr::mutate(tie_frac = round(tie_frac, 4)) %>%
    utils::head(10) %>% as.data.frame() %>% print(row.names = FALSE)
} else {
  message("\n   no ties in any scored variable")
}

# =============================================================================
# 4. rho(MYC, OXPHOS) IN ITS OWN RIGHT - reported, not chased
# =============================================================================
# Load-bearing well beyond this script. The collider argument and the
# "MYC-high tumours are OXPHOS-high" premise BOTH assume this is positive, and
# the dose estimator says it may not be. 5 estimators x 2 rulers x 2 cohorts x
# 3 strata, with intervals.
#
# Spearman with Fisher-z intervals and the Bonett-Wright variance - E10's
# `.cor_block`, copied verbatim as E16, E18, E19, E20 and E21 copy it. Computed
# on the RAW vectors, not the Blom scores, so these are directly comparable to
# every rho already in this repo.
message("\n4. rho(MYC, OXPHOS) - the matrix, reported in its own right")

.rank_rows <- function(M) {
  out <- t(apply(M, 1L, rank))
  if (nrow(M) == 1L) out <- matrix(out, nrow = 1L)
  dimnames(out) <- dimnames(M); out
}
.cor_block <- function(A, Bm, cov = NULL) {
  n <- ncol(A)
  stopifnot(identical(colnames(A), colnames(Bm)))
  RA <- .rank_rows(A); RB <- .rank_rows(Bm)
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
    a = rep(rownames(R), times = ncol(R)), b = rep(colnames(R), each = nrow(R)),
    n = n, k_cov = k, rho = as.vector(R),
    ci_lo = as.vector(tanh(z - 1.959964 * se)),
    ci_hi = as.vector(tanh(z + 1.959964 * se)))
}

rho_myc_ox <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  Eraw <- do.call(rbind, lapply(EST_ALL, function(w) .estimator_vec(C, w)))
  dimnames(Eraw) <- list(names(EST_ALL), C$ids)
  dplyr::bind_rows(lapply(STRATA, function(s) {
    ids <- STR[[coh]][[s]]
    .cor_block(Eraw[, ids, drop = FALSE],
               rulers[[coh]][RULERS, ids, drop = FALSE]) %>%
      dplyr::mutate(cohort = coh, stratum = s)
  }))
})) %>%
  dplyr::rename(estimator = a, ruler = b) %>%
  dplyr::mutate(ci_width = ci_hi - ci_lo,
                claims = estimator %in% names(EST_CLAIM),
                positive = ci_lo > 0, negative = ci_hi < 0,
                estimator = factor(estimator, levels = names(EST_ALL)),
                stratum = factor(stratum, levels = STRATA)) %>%
  dplyr::arrange(cohort, stratum, ruler, estimator)

for (r in RULERS) {
  message("\n   rho(MYC, ", r, ") with 95% intervals:")
  rho_myc_ox %>%
    dplyr::filter(ruler == r) %>%
    dplyr::transmute(cohort, stratum, estimator, n,
                     rho = sprintf("%+.3f", rho),
                     ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
                     sign = dplyr::case_when(positive ~ "> 0", negative ~ "< 0",
                                             TRUE ~ "spans 0")) %>%
    as.data.frame() %>% print(row.names = FALSE)
}
message("\n   READ THIS BEFORE USING 'MYC-high tumours are OXPHOS-high' ANYWHERE.")

# =============================================================================
# 5. FIT (a) - the unadjusted two-predictor model
# =============================================================================
# BBC3 ~ MYC + OXPHOS on Blom scores. NO COVARIATES, NO PRODUCT TERM, POOLED
# (never across cohorts). All twelve outcomes are fitted because the roster is
# the specificity control and picking BBC3 after seeing them would be the
# grid-of-cells trap; BBC3 is the gate deliverable and is reported first.
message("\n5. fit (a): Y ~ MYC + OXPHOS, Blom scores, no covariates")

.fit_a <- function(y, m, o) {
  d  <- data.frame(Y = y, M = m, O = o)
  f  <- stats::lm(Y ~ M + O, data = d)
  ci <- stats::confint(f)
  cf <- stats::coef(f)
  tibble::tibble(
    beta1 = unname(cf["M"]), b1_lo = ci["M", 1], b1_hi = ci["M", 2],
    beta2 = unname(cf["O"]), b2_lo = ci["O", 1], b2_hi = ci["O", 2],
    r_MO  = stats::cor(m, o),
    r_YM_observed  = stats::cor(y, m),
    r_YO_observed  = stats::cor(y, o),
    n = length(y), adj_r2 = summary(f)$adj.r.squared)
}

fit_a <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  Bc <- B[[coh]]
  tidyr::expand_grid(gene = PRIMING_ALL, estimator = names(EST_ALL),
                     ruler = RULERS) %>%
    purrr::pmap_dfr(function(gene, estimator, ruler)
      .fit_a(Bc$Y[gene, ], Bc$E[estimator, ], Bc$O[ruler, ]) %>%
        dplyr::mutate(cohort = coh, gene = gene, estimator = estimator,
                      ruler = ruler, .before = 1L))
})) %>%
  dplyr::mutate(
    side = dplyr::if_else(gene %in% PRIMING_PRO, "pro-apoptotic",
                          "anti-apoptotic"),
    claims = estimator %in% names(EST_CLAIM),
    b1_negative = b1_hi < 0, b1_positive = b1_lo > 0,
    estimator = factor(estimator, levels = names(EST_ALL)),
    ruler = factor(ruler, levels = RULERS))

message("\n   ", FOCUS_GENE, ", the gate deliverable - the three CLAIM",
        " estimators first:")
fit_a %>%
  dplyr::filter(gene == FOCUS_GENE, claims) %>%
  dplyr::transmute(cohort, ruler, estimator, n,
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi),
                   beta2 = sprintf("%+.3f", beta2),
                   b2_ci = sprintf("[%+.3f, %+.3f]", b2_lo, b2_hi)) %>%
  dplyr::arrange(ruler, cohort, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   ", FOCUS_GENE, ", the two SENSITIVITY signatures (not in the rule):")
fit_a %>%
  dplyr::filter(gene == FOCUS_GENE, !claims) %>%
  dplyr::transmute(cohort, ruler, estimator,
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi),
                   beta2 = sprintf("%+.3f", beta2)) %>%
  dplyr::arrange(ruler, cohort, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.1 the D2 concordance rule, applied mechanically ----------------------
d2 <- fit_a %>%
  dplyr::filter(gene == FOCUS_GENE, claims) %>%
  dplyr::group_by(ruler) %>%
  dplyr::summarise(n_cells = dplyr::n(),
                   n_beta1_lt0 = sum(beta1 < 0),
                   n_ci_excludes_0 = sum(b1_negative),
                   concordant = all(beta1 < 0), .groups = "drop")
message("\n   D2 concordance on the SIGN of beta1 (3 estimators x 2 cohorts):")
d2 %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. THE IDENTITY - reported ONCE, labelled ACCOUNTING
# =============================================================================
# beta1 + beta2 * r(M,O) == r(Y,M). It cannot fail; it is asserted to 1e-10 and
# then read for the only thing it says, which is HOW NEGATIVE beta1 has to be.
message("\n6. the identity - ACCOUNTING, NOT EVIDENCE")

identity_chk <- fit_a %>%
  dplyr::filter(gene == FOCUS_GENE) %>%
  dplyr::mutate(r_YM_predicted = beta1 + beta2 * r_MO,
                delta = r_YM_predicted - r_YM_observed,
                required_b1 = r_YM_observed - beta2 * r_MO,
                vacuous = abs(r_MO) < VACUOUS_RMO)
if (max(abs(identity_chk$delta)) > 1e-10) {
  stop("the algebraic identity did not hold to 1e-10 (max |delta| = ",
       signif(max(abs(identity_chk$delta)), 3), "). Something is wrong with the ",
       "Blom scoring or the fit, not with the biology.", call. = FALSE)
}
message("   identity holds to ", signif(max(abs(identity_chk$delta)), 3),
        " in all ", nrow(identity_chk), " cells - as it must.")
identity_chk %>%
  dplyr::filter(claims) %>%
  dplyr::transmute(cohort, ruler, estimator,
                   beta2_x_rMO = sprintf("%+.3f", beta2 * r_MO),
                   beta1 = sprintf("%+.3f", beta1),
                   r_MO = sprintf("%+.3f", r_MO),
                   observed = sprintf("%+.3f", r_YM_observed),
                   predicted = sprintf("%+.3f", r_YM_predicted),
                   vacuous) %>%
  dplyr::arrange(ruler, cohort, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   `vacuous` = |r(M,O)| < ", VACUOUS_RMO, ": there is nothing to",
        " reconcile and\n   beta1 is simply the marginal. Trap 4. Never read",
        " such a row as agreement.")

# --- 6.1 Blom-Pearson beside Spearman, so nothing looks contradicted --------
marginals <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]; Bc <- B[[coh]]
  Eraw <- do.call(rbind, lapply(EST_ALL, function(w) .estimator_vec(C, w)))
  dimnames(Eraw) <- list(names(EST_ALL), C$ids)
  tibble::tibble(cohort = coh, estimator = names(EST_ALL),
                 spearman = vapply(names(EST_ALL), function(e)
                   .rho(as.numeric(GR[[coh]]$mat[FOCUS_GENE, ]), Eraw[e, ]),
                   numeric(1)),
                 blom_pearson = vapply(names(EST_ALL), function(e)
                   stats::cor(Bc$Y[FOCUS_GENE, ], Bc$E[e, ]), numeric(1)))
})) %>% dplyr::mutate(gap = blom_pearson - spearman)
message("\n   marginal rho(", FOCUS_GENE, ", MYC): Spearman (the repo's standing",
        " measure)\n   beside Blom-Pearson (the one the identity uses):")
marginals %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   ALL NEGATIVE. Any manuscript sentence calling ", FOCUS_GENE,
        " 'uncorrelated with MYC'\n   in human tumours is wrong as written and",
        " must read 'weakly negatively correlated'.")

# =============================================================================
# 7. THE LEAVE-ONE-OUT MEASUREMENT
# =============================================================================
# Eleven of the twelve outcomes sit in ox_rel's denominator with weight -1/1047.
# This measures the consequence instead of arguing it away. ox_lvl is untouched
# and is not refitted here.
message("\n7. leave-one-out ox_rel: the outcome removed from its own predictor")

loo <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  Bc <- B[[coh]]
  tidyr::expand_grid(gene = PRIMING_ALL, estimator = names(EST_ALL)) %>%
    purrr::pmap_dfr(function(gene, estimator) {
      std <- .fit_a(Bc$Y[gene, ], Bc$E[estimator, ], Bc$O["ox_rel", ])
      out <- .fit_a(Bc$Y[gene, ], Bc$E[estimator, ], Bc$OL[gene, ])
      tibble::tibble(cohort = coh, gene = gene, estimator = estimator,
                     b1 = std$beta1, b1_loo = out$beta1,
                     b2 = std$beta2, b2_loo = out$beta2,
                     d_b1 = out$beta1 - std$beta1,
                     d_b2 = out$beta2 - std$beta2,
                     in_denominator = gene %in% LOO_NEEDED)
    })
}))
message("\n   worst shift over all ", nrow(loo), " cells: |d beta1| = ",
        signif(max(abs(loo$d_b1)), 3), ", |d beta2| = ",
        signif(max(abs(loo$d_b2)), 3))
message("   ", FOCUS_GENE, ", the three claim estimators:")
loo %>%
  dplyr::filter(gene == FOCUS_GENE, estimator %in% names(EST_CLAIM)) %>%
  dplyr::transmute(cohort, estimator,
                   beta1 = sprintf("%+.4f", b1), beta1_loo = sprintf("%+.4f", b1_loo),
                   beta2 = sprintf("%+.4f", b2), beta2_loo = sprintf("%+.4f", b2_loo),
                   d_b1 = sprintf("%+.4f", d_b1), d_b2 = sprintf("%+.4f", d_b2)) %>%
  as.data.frame() %>% print(row.names = FALSE)
loo_ok <- max(abs(loo$d_b1)) < 0.01 && max(abs(loo$d_b2)) < 0.01
message("   both shifts under 0.01 everywhere: ", loo_ok,
        "\n   PRIMARY REMAINS THE STANDARD ox_rel, for comparability with E10,",
        " E16, E20\n   and E21; the leave-one-out value is the measured",
        " sensitivity beside it.")

# =============================================================================
# 8. GATE 2
# =============================================================================
message("\n", strrep("=", 78), "\n8. GATE 2 - STOP. Report this and wait.\n",
        strrep("=", 78))

gate2 <- identity_chk %>%
  dplyr::filter(claims) %>%
  dplyr::transmute(cohort, ruler, estimator, n,
                   beta1 = round(beta1, 3), b1_lo = round(b1_lo, 3),
                   b1_hi = round(b1_hi, 3),
                   beta2 = round(beta2, 3), b2_lo = round(b2_lo, 3),
                   b2_hi = round(b2_hi, 3),
                   r_MO = round(r_MO, 3),
                   r_YM_observed = round(r_YM_observed, 3),
                   r_YM_predicted = round(r_YM_predicted, 3),
                   vacuous) %>%
  dplyr::arrange(ruler, cohort, estimator)
gate2 %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   D2 SIGN CONCORDANCE: ",
        paste(sprintf("%s %s", d2$ruler, ifelse(d2$concordant, "YES", "NO")),
              collapse = " | "))
message("   NOT RUN, AND NOT TO BE RUN UNTIL THE GATE CLEARS: the adjusted fits",
        "\n   (purity, PROLIF_DISJOINT), the within-stratum fits, the twelve-gene",
        "\n   specificity table as a REPORTED object, and the gated FOXO3",
        " mediation.")

# =============================================================================
# 9. Save
# =============================================================================
saveRDS(list(
  contamination = contamination, strata_n = strata_n, ties = ties,
  rho_myc_ox = rho_myc_ox,
  fit_a = fit_a, d2 = d2, identity_chk = identity_chk, marginals = marginals,
  loo = loo, loo_ok = loo_ok, gate2 = gate2,
  settings = list(
    stage = "STAGE 1 ONLY - stops at GATE 2. Adjusted, stratified and mediation
             fits are appended to this file after the gate clears.",
    not_this = paste("MAIN-EFFECTS model. A different estimand from the",
                     "pre-registered MYC:OXPHOS interaction on PRIME, which is",
                     "null and permanently closed. No product term, no PRIME,",
                     "nothing read from myc_human_validation, and nothing here",
                     "bears on the registered interaction result."),
    estimators_fitted = EST_ALL, estimators_claiming = EST_CLAIM,
    estimators_disqualified = EST_DISQUALIFIED,
    disqualification_reason = paste("M_b__FULL and M_b__PROLIFSTRIP contain",
                                    "BBC3, BCL2, BCL2L1, BCL2L11 and PMAIP1.",
                                    "BBC3 is a genuine CollecTRI MYC target;",
                                    "the MITOSTRIP variant removes it because",
                                    "BBC3 is a MitoCarta gene."),
    rulers = RULERS, strata = STRATA, focus = FOCUS_GENE,
    scoring = paste("Blom normal scores, qnorm((rank - 3/8)/(n + 1/4)), THEN",
                    "scaled to unit variance, pooled per cohort. The scaling",
                    "matters only where there are ties; without it SCAN-B's",
                    "variances spread by 4.4e-4 and the identity would hold",
                    "only to ~1e-4. See $ties."),
    d2_rule = D2_RULE, identity_note = IDENTITY_NOTE, basal_note = BASAL_NOTE,
    vacuous_rmo = VACUOUS_RMO,
    purity_provenance = paste("source UNDOCUMENTED IN THIS REPO. E01 reads the",
                              "column from data/from_validation/; that README",
                              "names the column and no algorithm. It sits beside",
                              "ploidy and genome_doublings, an ABSOLUTE",
                              "fingerprint - recorded as INFERENCE, not",
                              "documentation. genome_doublings is not a",
                              "covariate and is not added. SCAN-B has none and",
                              "nothing is substituted."),
    seed = PROJECT_SEED),
  rules = list(
    n3 = "transcript associations; 'primed' is never written of a transcript",
    cohorts = "fitted separately, never pooled",
    mouse = "the mouse is not read, compared or cited by this script",
    identity_scope = paste("the identity is claimed for the POOLED fit only.",
                           "Stage 2 subsets pooled Blom scores per E19's",
                           "no-re-scoring rule, which breaks equal variance,",
                           "so it is not computed within strata.")),
  built = Sys.time()), PATH_E22)

readr::write_csv(fit_a %>% dplyr::mutate(dplyr::across(where(is.factor), as.character)),
                 PATH_E22_FITA)
readr::write_csv(rho_myc_ox %>% dplyr::mutate(dplyr::across(where(is.factor), as.character)),
                 PATH_E22_RMO)

message("\nE22 stage 1: done.")
message("    results/bbc3_direct_myc_effect.rds")
message("    outputs/tables/E22_fit_a.csv")
message("    outputs/tables/E22_rho_myc_oxphos.csv")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E22)

  # THE GATE 2 TABLE. This is what gets reported and waited on.
  x$gate2 %>% as.data.frame()
  x$d2 %>% as.data.frame()

  # Is BBC3 really negatively correlated with MYC? Spearman and Blom-Pearson.
  x$marginals %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # The identity, and which rows are vacuous because r(M,O) is ~ 0.
  x$identity_chk %>%
    dplyr::transmute(cohort, ruler, estimator, r_MO = round(r_MO, 3),
                     observed = round(r_YM_observed, 3),
                     predicted = round(r_YM_predicted, 3),
                     beta1 = round(beta1, 3), vacuous) %>%
    as.data.frame()

  # rho(MYC, OXPHOS) IN ITS OWN RIGHT - read before using "MYC-high tumours are
  # OXPHOS-high" anywhere, in this repo or the manuscript.
  x$rho_myc_ox %>%
    dplyr::filter(ruler == "ox_rel") %>%
    dplyr::transmute(cohort, stratum, estimator, rho = round(rho, 3),
                     ci = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi), claims) %>%
    as.data.frame()

  # Did the outcome sitting in its own predictor matter? Expect no.
  summary(abs(x$loo$d_b1)); summary(abs(x$loo$d_b2)); x$loo_ok

  # All twelve on fit (a) - the specificity roster, NOT yet a reported result.
  x$fit_a %>%
    dplyr::filter(estimator == "signature", ruler == "ox_rel") %>%
    dplyr::transmute(cohort, gene, side, beta1 = round(beta1, 3),
                     beta2 = round(beta2, 3)) %>%
    dplyr::arrange(cohort, beta1) %>%
    as.data.frame()

  # What contaminates what, and why M_b__FULL is excluded.
  x$contamination %>% as.data.frame()

}
