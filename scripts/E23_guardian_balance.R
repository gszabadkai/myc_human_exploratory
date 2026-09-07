# E23_guardian_balance.R
# =============================================================================
# CONFIRMATORY: DOES MYC RAISE THE GUARDIAN BALANCE AT FIXED OXPHOS?
#
#     outcome  B = log2(MCL1) - log2(BCL2L1)
#     model    B ~ MYC + OXPHOS, Blom normal scores, unit variance
#     PREDICTED   beta1 > 0
#
# Same model, same ladder, same battery as E22. The endpoint only counts if it
# clears ALL of it - the battery is not a menu.
#
# =============================================================================
# WHAT THIS IS NOT - carried forward from E22, and it binds identically
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
# WHAT LICENSES THE PREDICTION - AND WHAT DOES NOT. READ THIS FIRST.
# =============================================================================
# LICENSED BY: `myc_mouse` docs/2026-09-02_myc_oxphos_priming_gate_model.md,
# read read-only and cited as CONCLUSIONS, never as values. There,
# `B = log2(Mcl1) - log2(Bcl2l1)` is the guardian balance, promoted to
# co-primary in section 4.3, and it is the model's most robust endpoint:
# section 3.3's endpoint menu puts `Mcl1:Bcl2l1` top on both the full model and
# the within-timepoint (cross-sectional) regime, and section 3.3's per-gene
# table has `Mcl1` at the 98.9th percentile of its expression-matched null and
# `Bcl2l1` at the 2.1st - the two halves moving oppositely is the endpoint's
# whole content.
#
# NOT LICENSED BY: E22's specificity table. That table is where the author
# NOTICED this endpoint - `BCL2L1` 6 of 6 negative, `MCL1` 0 of 6 - and noticing
# is not licensing. The distinction is stated here because it decides how the
# result may be written, and the next paragraph is why it matters so much.
#
# =============================================================================
# THIS IS NOT AN INDEPENDENT TEST IN THIS DATA, AND THAT IS SAID UP FRONT
# =============================================================================
# E22 already fitted both halves of this ratio, in the same samples, with the
# same estimators and the same rulers. At the top rung, across both rulers x
# three claim estimators x two cohorts (12 cells each):
#
#     BCL2L1   beta1 NEGATIVE in 12 of 12, every interval excluding zero
#     MCL1     beta1 POSITIVE in 12 of 12 point estimates
#
# A difference of those two is positive with near-certainty BEFORE anything is
# fitted. So `beta1 > 0` on the ratio is close to arithmetically assured here,
# and CONFIRMING IT IS WORTH ALMOST NOTHING. This is handoff trap 6i in a new
# guise: the same table asked a different question is not a second finding.
#
# TWO THINGS ARE STILL WORTH DOING, and they are the whole point of the script:
#
#   1. It puts the human arm on the MOUSE'S OWN ENDPOINT, so the cross-species
#      comparison is made on the mouse's terms rather than on a human-chosen
#      surrogate. That is a framing gain, not an evidential one.
#
#   2. THE BATTERY CAN STILL FAIL, and a failure IS informative. Section 0's
#      failure reading is pre-specified, and if it fires the endpoint is written
#      as a SECOND NEGATIVE - never softened into "attenuated".
#
# Section 6 measures point 1 directly: `beta1(ratio)` beside
# `beta1(MCL1) - beta1(BCL2L1)`. If the gap is negligible the ratio is a
# re-expression of two numbers already reported, and the note must say so.
#
# =============================================================================
# SCALES, STATED ONCE - identical to E22
# =============================================================================
#   MCL1, BCL2L1, B  log2(linear DESeq2-normalised + 1); B is their difference
#   ox_rel, ox_lvl   BUILT here from the same matrix, E16/E20/E22's `.comp`
#   log2MYC          READ from the linear matrix
#   GSVA signatures  READ: GSVA on the VST, kcdf Gaussian
#   M_b__MITOSTRIP   READ: decoupleR ULM on the VST
#
# Every model variable is Blom normal-scored, qnorm((rank - 3/8)/(n + 1/4)),
# THEN scaled to unit variance - E22's step, and for E22's reason: ties break
# the equal-variance property that scoring alone would give.
#
# BLOM SCORING IS POOLED, ONCE PER COHORT. Stratum fits SUBSET those scores and
# never re-score (E19's rule).
#
# NEVER POOLED ACROSS COHORTS - fitted separately, always. N3 THROUGHOUT: these
# are transcript associations, and "primed" is never written of a transcript.
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))
source(here::here("functions", "strata.R"))

message("\nE23: the guardian balance, MCL1:BCL2L1, against MYC at fixed OXPHOS\n",
        strrep("=", 78))

PATH_E23      <- file.path(DIR_RESULTS, "guardian_balance.rds")
PATH_E23_TAB  <- file.path(DIR_TABLES,  "E23_guardian_balance.csv")
PATH_E22      <- file.path(DIR_RESULTS, "bbc3_direct_myc_effect.rds")

# =============================================================================
# 0. PRE-SPECIFICATION - fixed before any coefficient in this script existed
# =============================================================================
NUM_GENE <- "MCL1"      # numerator of the balance
DEN_GENE <- "BCL2L1"    # denominator - the mouse's shared denominator
PAIR     <- c(NUM_GENE, DEN_GENE)
RATIO    <- "MCL1_minus_BCL2L1"

EST_CLAIM <- c(dose      = "log2MYC",
               signature = "FELSHER__MITOSTRIP",
               regulon   = "M_b__MITOSTRIP")
EST_SENS  <- c(sig_clean     = "MYC_UP.V1_UP__MITOSTRIP",              # 1.5 pct
               sig_entangled = "HALLMARK_MYC_TARGETS_V1__MITOSTRIP")   # 23.5 pct
EST_ALL   <- c(EST_CLAIM, EST_SENS)
# Proliferation entanglement, from CLAUDE.md trap 3. Used ONLY to order the
# three signature estimators for the failure reading - never as a weight.
ENTANGLE  <- c(sig_clean = 1.5, signature = 14.8, sig_entangled = 23.5)

RULERS <- c("ox_rel", "ox_lvl")
STRATA <- c("Luminal", "Basal")
MATERIAL_DELTA <- 0.05   # E16's R1_MAX_DELTA, as E20, E21 and E22 borrow it

# --- THE SUCCESS READING ------------------------------------------------------
# The endpoint counts ONLY if it clears the whole battery. Any one of these
# failing is a failure of the endpoint, not a reason to report the rest.
SUCCESS_RULE <- paste0(
  "beta1 > 0 with its interval excluding zero on ALL THREE claim estimators, ",
  "in BOTH cohorts, on BOTH rulers, at the TOP RUNG, POOLED and WITHIN ",
  "LUMINAL - and neither failure reading firing.")

# --- THE FAILURE READING, PRE-SPECIFIED --------------------------------------
# Two ways to fail, both taken from what E22 found for BBC3. If either fires,
# THE ENDPOINT IS WRITTEN AS A SECOND NEGATIVE - the same artefact as BBC3 -
# and NOT softened into "attenuated".
#
# F1 ENTANGLEMENT TRACKING. beta1 across the three SIGNATURE estimators, ordered
#    by proliferation entanglement (1.5, 14.8, 23.5), is MONOTONE in both
#    cohorts AND the spread exceeds MATERIAL_DELTA. Either direction counts:
#    a readout whose value is ordered by how much proliferation it carries is
#    measuring entanglement, not MYC.
#
#    HONESTLY LABELLED: with three points, monotone-by-chance is 1 in 3 per
#    cohort and 1 in 9 for both, so this is a HEURISTIC FLAG and not a test.
#    The spread condition is what stops it firing on three indistinguishable
#    values. Section 7 applies the identical rule to E22's BBC3 as a POSITIVE
#    CONTROL: it is known to be entanglement-tracking, so the rule must fire
#    there or the rule is not doing its job.
#
# F2 LUMINAL REVERSAL. beta1 <= 0 within Luminal in EITHER cohort at the top
#    rung. Luminal is the large, homogeneous stratum; a sign that survives
#    pooling but not Luminal is a composition artefact.
FAILURE_RULE <- paste0(
  "F1 entanglement tracking: beta1 monotone across sig_clean/signature/",
  "sig_entangled in BOTH cohorts with spread > ", MATERIAL_DELTA,
  ". F2 Luminal reversal: beta1 <= 0 within Luminal in either cohort at the ",
  "top rung. EITHER firing = a SECOND NEGATIVE, written as one, never as ",
  "'attenuated'.")

# beta2 IS NOT PREDICTED, and that is deliberate. The mouse's own statement for
# this endpoint is a CROSSOVER: `bX` is -0.541 at MYC-null and MYC reverses the
# sign, so the OXPHOS main effect pooled over a mixed-MYC tumour cohort has no
# predicted direction. It is reported, never scored.
BETA2_NOTE <- paste0(
  "beta2 is REPORTED AND NOT PREDICTED. The mouse's statement for this endpoint ",
  "is a crossover (bX = -0.541 at MYC-null, sign reversed by MYC), so pooled ",
  "over a mixed-MYC cohort its direction is not fixed in advance.")

# The estimand differs from the mouse's, and saying so is not a caveat but a
# definition. The mouse's `bMX` is an INTERACTION coefficient; `beta1` here is a
# MAIN EFFECT. What transfers is the PREDICTED SIGN on the guardian balance,
# not the estimand.
ESTIMAND_NOTE <- paste0(
  "the mouse's bMX is an INTERACTION coefficient; beta1 here is a MAIN EFFECT. ",
  "The predicted SIGN is inherited, the estimand is not. No interaction is ",
  "fitted anywhere in this script.")

message("\n0. pre-specification, fixed before any coefficient")
message("   outcome      ", RATIO, " = log2(", NUM_GENE, ") - log2(", DEN_GENE, ")")
message("   predicted    beta1 > 0")
message("   success      ", SUCCESS_RULE)
message("   failure      ", FAILURE_RULE)
message("   beta2        ", BETA2_NOTE)
message("   estimand     ", ESTIMAND_NOTE)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

if (!file.exists(PATH_E22)) {
  stop("results/bbc3_direct_myc_effect.rds is missing. E23 reproduces E22's ",
       "single-gene MCL1 and BCL2L1 fits as its control and cannot run without ",
       "it. Run E22 first.", call. = FALSE)
}
e22 <- readRDS(PATH_E22)
stopifnot("ladder_pooled" %in% names(e22))

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
fr   <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)

OX_SUB   <- mito$arm_sets[["OXPHOS subunits"]]
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL
REST     <- setdiff(MITO_ALL, OX_SUB)
PROLIF_COV <- "PROLIF_DISJOINT"
PD <- mito$covariate_sets[[PROLIF_COV]]

# --- 1.1 contamination, asserted as in E22 ----------------------------------
EST_SETS <- lapply(EST_ALL[EST_ALL != "log2MYC"], function(nm)
  if (nm %in% names(sd_$myc_sets)) sd_$myc_sets[[nm]] else sd_$collectri_sets[[nm]])
names(EST_SETS) <- EST_ALL[EST_ALL != "log2MYC"]
contamination <- tibble::tibble(
  predictor = c(names(EST_SETS), PROLIF_COV, "ox_rel numerator",
                "ox_rel denominator"),
  n = c(lengths(EST_SETS), length(PD), length(OX_SUB), length(REST)),
  which_of_pair = vapply(c(EST_SETS, list(PD, OX_SUB, REST)),
                         function(s) paste(intersect(PAIR, s), collapse = ", "),
                         character(1)))
message("\n   does either half of the ratio sit inside a predictor?")
contamination %>% as.data.frame() %>% print(row.names = FALSE)
# Both halves are MitoCarta members, so BOTH are in ox_rel's denominator and
# section 2 builds a pair-specific leave-one-out. Neither is in any fitted
# estimator (the __MITOSTRIP strip removed BCL2L1 from the regulon; MCL1 was
# never in it) nor in the covariate nor in the numerator - asserted, because a
# non-zero would mean the outcome is partly its own predictor.
stopifnot(
  all(vapply(EST_SETS, function(s) length(intersect(PAIR, s)) == 0L, logical(1))),
  length(intersect(PAIR, PD)) == 0L,
  length(intersect(PAIR, OX_SUB)) == 0L,
  setequal(intersect(PAIR, REST), PAIR))
message("   both halves are in ox_rel's DENOMINATOR (both are MitoCarta genes)",
        " and in\n   nothing else. Section 2 carries the pair-specific",
        " leave-one-out.")

tcga_lin  <- readRDS(PATH_TCGA_LINEAR); scanb_lin <- readRDS(PATH_SCANB_LINEAR)
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

GR <- lapply(COH, function(C) .gene_rows(PAIR, C$G, C$res))
miss <- unique(unlist(lapply(GR, function(g) g$missing), use.names = FALSE))
if (length(miss)) stop("did not resolve: ", paste(miss, collapse = ", "),
                       call. = FALSE)
message("\n   TCGA ", length(ID_T), " | SCAN-B ", length(ID_S),
        " samples; both halves present in both cohorts")

# =============================================================================
# 2. Rulers, including the PAIR-SPECIFIC leave-one-out
# =============================================================================
# Both halves of the outcome sit in ox_rel's denominator, so the leave-one-out
# here removes BOTH at once rather than one at a time. ox_lvl needs none -
# neither gene is an OXPHOS subunit, asserted in section 1.
message("\n2. the rulers")

.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res); M <- gr$mat
  v <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}
rulers <- lapply(COH, function(C) {
  num <- .comp(OX_SUB, C)
  M <- rbind(ox_rel     = num - .comp(REST, C),
             ox_lvl     = num,
             ox_rel_loo = num - .comp(setdiff(REST, PAIR), C))
  colnames(M) <- C$ids; stopifnot(!anyNA(M)); M
})
message("   ox_rel, ox_lvl, and ox_rel_loo (both halves dropped from the",
        " denominator)")

# =============================================================================
# 3. Blom scores - E22's function, verbatim
# =============================================================================
message("\n3. Blom normal scores, then unit variance")

.blom <- function(x) {
  n <- length(x); stopifnot(!anyNA(x))
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.375) / (n + 0.25))
  as.numeric(scale(z))
}
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
  raw_num <- as.numeric(GR[[coh]]$mat[NUM_GENE, ])
  raw_den <- as.numeric(GR[[coh]]$mat[DEN_GENE, ])
  # THE OUTCOME. The ratio is formed on the LOG scale, so it is a difference -
  # log2(MCL1) - log2(BCL2L1) - and only then normal-scored. Forming it after
  # scoring would be a different object.
  Y <- rbind(ratio = .blom(raw_num - raw_den),
             MCL1  = .blom(raw_num),
             BCL2L1 = .blom(raw_den))
  E <- do.call(rbind, lapply(EST_ALL, function(w) .blom(.estimator_vec(C, w))))
  rownames(E) <- names(EST_ALL)
  O <- do.call(rbind, lapply(rownames(rulers[[coh]]),
                             function(r) .blom(rulers[[coh]][r, ])))
  rownames(O) <- rownames(rulers[[coh]])
  colnames(Y) <- colnames(E) <- colnames(O) <- C$ids
  vs <- c(apply(Y, 1L, stats::var), apply(E, 1L, stats::var),
          apply(O, 1L, stats::var))
  if (max(abs(vs - 1)) > 1e-10) {
    stop("Blom scores are not unit variance in ", coh, call. = FALSE)
  }
  B[[coh]] <- list(Y = Y, E = E, O = O)
  ties[[coh]] <- tibble::tibble(
    cohort = coh, variable = c(RATIO, PAIR),
    tie_frac = c(.tie_frac(raw_num - raw_den), .tie_frac(raw_num),
                 .tie_frac(raw_den)))
  message("   ", coh, ": outcome + 2 halves, ", nrow(E), " estimators, ",
          nrow(O), " rulers; all unit variance")
}
ties <- dplyr::bind_rows(ties)
message("   tie fractions on the outcome and its halves:")
ties %>% dplyr::mutate(tie_frac = signif(tie_frac, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4. The ladder - identical construction to E22
# =============================================================================
# ONE FIXED SUBSET PER COHORT so every rung is comparable: TCGA the
# purity-complete set, SCAN-B everything. SCAN-B has no purity and nothing is
# substituted; its ladder stops one rung short.
message("\n4. the ladder")

PL <- fr %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::select(sample_id, purity, leuko)
PL <- PL[match(ID_T, PL$sample_id), ]
IDS_PUR <- ID_T[!is.na(PL$purity)]
COVV <- list(
  TCGA = list(prolif = stats::setNames(as.numeric(mito$gsva_cov[PROLIF_COV, ID_T]), ID_T),
              purity = stats::setNames(PL$purity, ID_T)),
  `SCAN-B` = list(prolif = stats::setNames(as.numeric(sc$gsva_cov[PROLIF_COV, ID_S]), ID_S)))
LADDER <- list(
  TCGA = list(subset = IDS_PUR,
              rungs = list(unadjusted = character(0), `+PROLIF` = "prolif",
                           `+PROLIF+purity` = c("prolif", "purity"))),
  `SCAN-B` = list(subset = ID_S,
                  rungs = list(unadjusted = character(0), `+PROLIF` = "prolif")))
TOP_RUNG <- c(TCGA = "+PROLIF+purity", `SCAN-B` = "+PROLIF")
RUNG_ORDER <- c("unadjusted", "+PROLIF", "+PROLIF+purity")
message("   TCGA n = ", length(IDS_PUR), " (purity-complete) | SCAN-B n = ",
        length(ID_S), "; SCAN-B has no purity rung and nothing is substituted")

.fit <- function(y, m, o, cov = NULL) {
  d <- data.frame(Y = y, M = m, O = o)
  if (!is.null(cov) && ncol(cov)) d <- cbind(d, cov)
  f <- stats::lm(Y ~ ., data = d); ci <- stats::confint(f); cf <- stats::coef(f)
  tibble::tibble(beta1 = unname(cf["M"]), b1_lo = ci["M", 1], b1_hi = ci["M", 2],
                 beta2 = unname(cf["O"]), b2_lo = ci["O", 1], b2_hi = ci["O", 2],
                 n = length(y))
}
.cov_frame <- function(coh, which, ids) {
  if (!length(which)) return(NULL)
  M <- do.call(cbind, lapply(which, function(w) .blom(COVV[[coh]][[w]][ids])))
  colnames(M) <- which; as.data.frame(M)
}
.ladder <- function(coh, ids, label, outcomes, rulers_use) {
  Bc <- B[[coh]]; L <- LADDER[[coh]]; k <- match(ids, colnames(Bc$Y))
  dplyr::bind_rows(lapply(names(L$rungs), function(rn)
    tidyr::expand_grid(outcome = outcomes, estimator = names(EST_ALL),
                       ruler = rulers_use) %>%
      purrr::pmap_dfr(function(outcome, estimator, ruler)
        .fit(Bc$Y[outcome, k], Bc$E[estimator, k], Bc$O[ruler, k],
             .cov_frame(coh, L$rungs[[rn]], ids)) %>%
          dplyr::mutate(cohort = coh, stratum = label, rung = rn,
                        outcome = outcome, estimator = estimator,
                        ruler = ruler, .before = 1L))))
}

OUTCOMES <- c("ratio", NUM_GENE, DEN_GENE)
pooled <- dplyr::bind_rows(lapply(names(COH), function(coh)
  .ladder(coh, LADDER[[coh]]$subset, "all", OUTCOMES, rownames(rulers$TCGA)))) %>%
  dplyr::mutate(claims = estimator %in% names(EST_CLAIM),
                b1_positive = b1_lo > 0, b1_negative = b1_hi < 0,
                rung = factor(rung, levels = RUNG_ORDER))

message("\n   THE RATIO, every rung, claim estimators, ox_rel:")
pooled %>%
  dplyr::filter(outcome == "ratio", ruler == "ox_rel", claims) %>%
  dplyr::transmute(cohort, estimator, rung, n,
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi),
                   beta2 = sprintf("%+.3f", beta2)) %>%
  dplyr::arrange(cohort, estimator, rung) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   THE TWO HALVES SEPARATELY, top rung, ox_rel - which member",
        " carries it:")
pooled %>%
  dplyr::filter(outcome != "ratio", ruler == "ox_rel", claims,
                rung == TOP_RUNG[cohort]) %>%
  dplyr::transmute(cohort, outcome, estimator,
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi)) %>%
  dplyr::arrange(cohort, outcome, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   THE RATIO on both rulers and the leave-one-out, top rung:")
pooled %>%
  dplyr::filter(outcome == "ratio", claims, rung == TOP_RUNG[cohort]) %>%
  dplyr::transmute(cohort, ruler, estimator,
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi)) %>%
  dplyr::arrange(cohort, ruler, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   sig_clean and sig_entangled on the ratio, top rung, ox_rel -",
        " REPORTED, NOT DROPPED:")
pooled %>%
  dplyr::filter(outcome == "ratio", ruler == "ox_rel", !claims,
                rung == TOP_RUNG[cohort]) %>%
  dplyr::transmute(cohort, estimator,
                   entanglement_pct = unname(ENTANGLE[estimator]),
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Strata
# =============================================================================
message("\n5. within stratum - pooled scores SUBSET, never re-scored")

STRAT_IDS <- list()
for (coh in names(COH)) for (s in STRATA)
  STRAT_IDS[[paste(coh, s)]] <- intersect(STR[[coh]][[s]], LADDER[[coh]]$subset)
strat_n <- tibble::tibble(key = names(STRAT_IDS), n = lengths(STRAT_IDS))
strat_n %>% as.data.frame() %>% print(row.names = FALSE)

strata_fits <- dplyr::bind_rows(lapply(names(COH), function(coh)
  dplyr::bind_rows(lapply(STRATA, function(s)
    .ladder(coh, STRAT_IDS[[paste(coh, s)]], s, OUTCOMES, RULERS))))) %>%
  dplyr::mutate(claims = estimator %in% names(EST_CLAIM),
                b1_positive = b1_lo > 0, ci_width = b1_hi - b1_lo,
                rung = factor(rung, levels = RUNG_ORDER))

message("\n   THE RATIO by stratum, top rung, ox_rel, with interval WIDTH:")
strata_fits %>%
  dplyr::filter(outcome == "ratio", ruler == "ox_rel",
                rung == TOP_RUNG[cohort]) %>%
  dplyr::transmute(cohort, stratum, estimator, n,
                   beta1 = sprintf("%+.3f", beta1),
                   b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi),
                   width = sprintf("%.3f", ci_width), b1_positive) %>%
  dplyr::arrange(cohort, stratum, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   A BASAL NULL IS WEAK EVIDENCE, NOT FAILURE TO HOLD. Quote the width.")

# =============================================================================
# 6. IS THE RATIO MORE THAN ITS TWO HALVES? - the honesty check
# =============================================================================
# The header says this endpoint is close to arithmetically assured given E22.
# This measures it. If `beta1(ratio)` is near `beta1(MCL1) - beta1(BCL2L1)`,
# the ratio is a RE-EXPRESSION of two numbers already reported and the note
# must say so rather than presenting it as a third result.
message("\n6. is the ratio more than its two halves?")

additivity <- pooled %>%
  dplyr::filter(rung == TOP_RUNG[cohort]) %>%
  dplyr::select(cohort, estimator, ruler, outcome, beta1) %>%
  tidyr::pivot_wider(names_from = outcome, values_from = beta1) %>%
  dplyr::mutate(halves_difference = .data[[NUM_GENE]] - .data[[DEN_GENE]],
                gap = ratio - halves_difference)
additivity %>%
  dplyr::transmute(cohort, ruler, estimator,
                   MCL1 = sprintf("%+.3f", .data[[NUM_GENE]]),
                   BCL2L1 = sprintf("%+.3f", .data[[DEN_GENE]]),
                   halves_diff = sprintf("%+.3f", halves_difference),
                   ratio = sprintf("%+.3f", ratio),
                   gap = sprintf("%+.3f", gap)) %>%
  dplyr::arrange(cohort, ruler, estimator) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   median |gap| = ", sprintf("%.3f", stats::median(abs(additivity$gap))),
        ". A small gap means the ratio ADDS NOTHING to the two",
        "\n   coefficients E22 already reported, and the note must say so.")

# --- 6.1 reproduction control against E22 -----------------------------------
# The two single-gene fits are cells E22 already carries, on the same samples,
# the same estimators and the same rulers. They must come back bit-equal.
repro <- pooled %>%
  dplyr::filter(outcome %in% PAIR, ruler %in% RULERS) %>%
  dplyr::mutate(rung = as.character(rung)) %>%
  dplyr::inner_join(
    e22$ladder_pooled %>%
      dplyr::mutate(ruler = as.character(ruler)) %>%
      dplyr::select(cohort, stratum, rung, gene, estimator, ruler,
                    e22_b1 = beta1, e22_b2 = beta2),
    by = c("cohort", "stratum", "rung", "outcome" = "gene", "estimator", "ruler")) %>%
  dplyr::mutate(d1 = beta1 - e22_b1, d2 = beta2 - e22_b2)
stopifnot(nrow(repro) > 0)
if (max(abs(c(repro$d1, repro$d2))) >= 1e-10) {
  stop("E23 does not reproduce E22's MCL1/BCL2L1 cells (max |delta| = ",
       signif(max(abs(c(repro$d1, repro$d2))), 3), "). An input has moved.",
       call. = FALSE)
}
message("   6.1 ", nrow(repro), " of E22's single-gene cells reproduce exactly",
        " (max |delta| = ", signif(max(abs(c(repro$d1, repro$d2))), 3), ")")

# =============================================================================
# 7. THE READINGS - success and the two failures, applied mechanically
# =============================================================================
message("\n7. the readings")

# --- F1, and its POSITIVE CONTROL -------------------------------------------
# The identical rule is applied to E22's BBC3, which is KNOWN to be
# entanglement-tracking. If it does not fire there, the rule is not working and
# its verdict here means nothing.
.f1 <- function(d, outcome_col, outcome_val, ruler_use) {
  x <- d %>%
    dplyr::filter(.data[[outcome_col]] == outcome_val, ruler == ruler_use,
                  estimator %in% names(ENTANGLE)) %>%
    dplyr::filter(rung == TOP_RUNG[cohort]) %>%
    dplyr::select(cohort, estimator, beta1)
  out <- x %>%
    dplyr::group_by(cohort) %>%
    dplyr::summarise(
      b = list(beta1[match(names(ENTANGLE), estimator)]),
      spread = diff(range(beta1)), .groups = "drop") %>%
    dplyr::mutate(monotone = vapply(b, function(v)
      all(diff(v) > 0) || all(diff(v) < 0), logical(1)))
  tibble::tibble(ruler = ruler_use,
                 n_monotone = sum(out$monotone),
                 min_spread = min(out$spread),
                 fires = all(out$monotone) && min(out$spread) > MATERIAL_DELTA)
}
f1_control <- dplyr::bind_rows(lapply(RULERS, function(r)
  .f1(e22$ladder_pooled %>% dplyr::mutate(ruler = as.character(ruler),
                                          rung = as.character(rung)),
      "gene", "BBC3", r))) %>% dplyr::mutate(target = "E22 BBC3 (control)")
f1_here <- dplyr::bind_rows(lapply(RULERS, function(r)
  .f1(pooled %>% dplyr::mutate(rung = as.character(rung)), "outcome", "ratio", r))) %>%
  dplyr::mutate(target = RATIO)
f1 <- dplyr::bind_rows(f1_control, f1_here) %>%
  dplyr::select(target, ruler, n_monotone, min_spread, fires)
message("\n   F1 entanglement tracking - the control row MUST fire:")
f1 %>% dplyr::mutate(min_spread = round(min_spread, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
if (!all(f1$fires[f1$target == "E22 BBC3 (control)"])) {
  message("   *** F1's POSITIVE CONTROL DID NOT FIRE ON BBC3. The rule is not",
          " detecting a\n   *** pattern known to be there, so its verdict on",
          " the ratio carries no weight.\n   *** Report the control failure",
          " and do not use F1 to adjudicate.")
}
F1_VALID <- all(f1$fires[f1$target == "E22 BBC3 (control)"])
F1_FIRES <- any(f1$fires[f1$target == RATIO])

# --- F2 Luminal reversal ----------------------------------------------------
lum <- strata_fits %>%
  dplyr::filter(outcome == "ratio", stratum == "Luminal", claims,
                rung == TOP_RUNG[cohort], ruler %in% RULERS)
f2 <- lum %>%
  dplyr::group_by(ruler) %>%
  dplyr::summarise(n_cells = dplyr::n(), n_b1_le0 = sum(beta1 <= 0),
                   min_beta1 = min(beta1),
                   fires = any(beta1 <= 0), .groups = "drop")
message("\n   F2 Luminal reversal (beta1 <= 0 in either cohort, top rung):")
f2 %>% dplyr::mutate(min_beta1 = round(min_beta1, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
F2_FIRES <- any(f2$fires)

# --- success ----------------------------------------------------------------
succ <- dplyr::bind_rows(lapply(RULERS, function(r) {
  p <- pooled %>% dplyr::filter(outcome == "ratio", ruler == r, claims,
                                rung == TOP_RUNG[cohort])
  l <- strata_fits %>% dplyr::filter(outcome == "ratio", ruler == r, claims,
                                     stratum == "Luminal",
                                     rung == TOP_RUNG[cohort])
  tibble::tibble(ruler = r,
                 pooled_ok = all(p$b1_positive), luminal_ok = all(l$b1_positive),
                 n_pooled = nrow(p), n_luminal = nrow(l))
}))
message("\n   success battery:")
succ %>% as.data.frame() %>% print(row.names = FALSE)
SUCCESS <- all(succ$pooled_ok) && all(succ$luminal_ok) && !F1_FIRES && !F2_FIRES

VERDICT <- if (SUCCESS) "CONFIRMED - beta1 > 0 on the guardian balance, whole battery clear" else
  if (F1_FIRES || F2_FIRES)
    "SECOND NEGATIVE - the same artefact as BBC3. NOT 'attenuated'." else
      "NOT CONFIRMED - the battery is not clear and no failure reading fired"
message("\n   VERDICT: ", VERDICT)
if (F1_FIRES) message("   F1 fired: beta1 tracks proliferation entanglement.")
if (F2_FIRES) message("   F2 fired: beta1 <= 0 within Luminal.")

verdicts <- tibble::tibble(
  item = c("verdict", "pooled battery clear", "Luminal battery clear",
           "F1 entanglement tracking fired", "F1 control valid",
           "F2 Luminal reversal fired", "median |ratio - halves difference|"),
  value = c(VERDICT, as.character(all(succ$pooled_ok)),
            as.character(all(succ$luminal_ok)), as.character(F1_FIRES),
            as.character(F1_VALID), as.character(F2_FIRES),
            sprintf("%.3f", stats::median(abs(additivity$gap)))))
message("")
verdicts %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. Save
# =============================================================================
saveRDS(list(
  contamination = contamination, ties = ties, strat_n = strat_n,
  pooled = pooled, strata_fits = strata_fits,
  additivity = additivity, repro = repro,
  f1 = f1, f2 = f2, succ = succ, verdicts = verdicts, verdict = VERDICT,
  f1_valid = F1_VALID, f1_fires = F1_FIRES, f2_fires = F2_FIRES,
  success = SUCCESS,
  settings = list(
    outcome = paste0(RATIO, " = log2(", NUM_GENE, ") - log2(", DEN_GENE, ")"),
    prediction = "beta1 > 0",
    licensed_by = paste("myc_mouse docs/2026-09-02_myc_oxphos_priming_gate_model.md,",
                        "read read-only and cited as conclusions only:",
                        "Mcl1:Bcl2l1 is the co-primary of section 4.3 and the",
                        "top endpoint of section 3.3's menu in both the full",
                        "and within-timepoint regimes."),
    not_licensed_by = paste("E22's specificity table. That is where the",
                            "endpoint was NOTICED. Noticing is not licensing."),
    not_independent = paste("E22 already fitted both halves in the same samples:",
                            "BCL2L1 negative 12 of 12 with intervals excluding",
                            "zero, MCL1 positive 12 of 12. beta1 > 0 on the",
                            "difference is close to arithmetically assured, so",
                            "confirming it is worth almost nothing. Section 6",
                            "measures how little the ratio adds. The value of",
                            "the script is that the BATTERY CAN FAIL."),
    estimand_note = ESTIMAND_NOTE, beta2_note = BETA2_NOTE,
    success_rule = SUCCESS_RULE, failure_rule = FAILURE_RULE,
    entanglement = ENTANGLE, material_delta = MATERIAL_DELTA,
    estimators = EST_ALL, rulers = RULERS, strata = STRATA,
    seed = PROJECT_SEED),
  rules = list(
    scope = paste("MAIN-EFFECTS model. No product term, no PRIME, nothing read",
                  "from myc_human_validation, and nothing here bears on the",
                  "registered interaction result."),
    n3 = "transcript associations; 'primed' is never written of a transcript",
    cohorts = "fitted separately, never pooled",
    failure_language = paste("if F1 or F2 fires the endpoint is a SECOND",
                             "NEGATIVE and is written as one - never softened",
                             "into 'attenuated'.")),
  built = Sys.time()), PATH_E23)

readr::write_csv(
  dplyr::bind_rows(pooled %>% dplyr::mutate(panel = "pooled"),
                   strata_fits %>% dplyr::mutate(panel = "stratum")) %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character)),
  PATH_E23_TAB)

message("\nE23: done.")
message("    results/guardian_balance.rds")
message("    outputs/tables/E23_guardian_balance.csv")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E23)

  # The verdict, and the rules it was scored on.
  x$verdicts %>% as.data.frame()
  x$settings$success_rule; x$settings$failure_rule
  x$settings$not_independent   # read this before quoting anything below

  # THE CONTROL: E22's MCL1 and BCL2L1 cells must come back bit-equal.
  nrow(x$repro); summary(abs(c(x$repro$d1, x$repro$d2)))

  # The ratio, every rung, claim estimators.
  x$pooled %>%
    dplyr::filter(outcome == "ratio", ruler == "ox_rel", claims) %>%
    dplyr::transmute(cohort, estimator, rung, beta1 = round(beta1, 3),
                     b1_ci = sprintf("[%+.3f, %+.3f]", b1_lo, b1_hi)) %>%
    as.data.frame()

  # Which member carries it - the whole reason both halves are reported.
  x$additivity %>%
    dplyr::filter(ruler == "ox_rel") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # F1 with its positive control. If the control row does not fire, F1's
  # verdict on the ratio carries no weight.
  x$f1 %>% as.data.frame()

  # F2, and the strata with widths.
  x$f2 %>% as.data.frame()
  x$strata_fits %>%
    dplyr::filter(outcome == "ratio", ruler == "ox_rel") %>%
    dplyr::transmute(cohort, stratum, estimator, rung, n,
                     beta1 = round(beta1, 3), width = round(ci_width, 3)) %>%
    as.data.frame()

  # The two sensitivity signatures, ordered by entanglement.
  x$pooled %>%
    dplyr::filter(outcome == "ratio", ruler == "ox_rel", !claims) %>%
    dplyr::transmute(cohort, estimator, rung, beta1 = round(beta1, 3)) %>%
    as.data.frame()

}
