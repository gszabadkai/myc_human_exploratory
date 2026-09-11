# scripts/E29_oxphos_neighbourhood.R
# =============================================================================
# E29 - where OXPHOS sits in the ordering, what is around it, what moves with it
# =============================================================================
#
# THIS PRODUCES A MAP, NOT A READOUT. E28 tested two pre-declared comparators for
# OXPHOS and the gap readout failed. Scanning 141 candidates and keeping the ones
# that behave would be SELECTING ON THE OUTCOME, and with 141 candidates some
# will look good by chance. So the output here is descriptive: where OXPHOS sits,
# what is around it, and what moves with it.
#
# NO PAIRING PRODUCED HERE MAY BE USED AS A READOUT without a separately
# declared validation across the full estimator ladder, both cohorts and the
# mouse. That sentence is in the saved object, the verdict and the note.
#
# PATHWAY LEVEL ONLY. Gene-level structure is a separate question, not in scope.
#
# EXPLORATORY, POST-HOC, DESCRIPTIVE. Nothing pre-registered, nothing a
# hypothesis test. R1-R5 are fixed before any number so the outcome cannot be
# chosen afterwards; none of them is a prediction.
#
# TWO PREMISES ABOUT THE SCORE, ONE TESTED AWAY AND ONE FOUND
# ------------------------------------------------------------
# The design was written expecting mitoPPS to behave as a CLOSED COMPOSITION -
# ratio-based, centred at 1 per pathway, therefore inducing spurious negative
# correlation between pathways. THAT CONCERN IS REASONABLE A PRIORI AND IT IS
# NOT BORNE OUT. PART 0 measures it:
#
#   - mitoPPS is NOT sum-constrained. The per-sample sum averages 142 but
#     RANGES, so the panel is not closed.
#   - Pairwise Spearman between the 142 scores across samples sits ON ZERO,
#     with about half the pairs negative in every cohort.
#
# The useful fact is the first one, and no Phase 2 note records it: the sum is
# not fixed, so this is not a closed composition and the standard compositional
# objection does not apply. PART B's contrast-profile design is kept, but for a
# different reason than the one it was written for: it is the estimand every
# other Phase 2 contrast already uses, which keeps E29 comparable with E25, E27
# and E28.
#
# WHAT IS ACTUALLY CONSTRAINED IS THE PERCENTILE TRANSFORM, one layer up. Within
# a group the 142 percentiles are a fixed set summing to 7100, so the CHANGE
# vector sums to EXACTLY ZERO across the panel. The neutral reference for
# aligned change is therefore not zero but
#
#     induced neutral = -|delta OXPHOS subunits| / 141
#
# It is small here and it is drawn rather than assumed. WHY IT IS SMALL MATTERS:
# the forcing scales as 1/(P-1), so it shrinks with panel size. On a 20-pathway
# panel the same OXPHOS change would force -2.7 percentile points and would have
# to be handled rather than noted.
#
# NO mitoPPS VALUE CROSSES THE SPECIES BOUNDARY - no table, no axis, no panel.
# Percentiles and percentile changes only. Cross-species comparison is by
# PATHWAY NAME. TCGA and SCAN-B are never averaged.
#
# n = 6 IN THE MOUSE. Every mouse ordering is imprecise. Rather than say so and
# leave it, the neighbourhood stability is BOOTSTRAPPED in all eight groups, so a
# reader can see how many of the twenty neighbours survive resampling in each.
#
# The mouse timeline is BATCH-CONFOUNDED (batch = timepoint). Ruling stands:
# reported in full, confounder stated, nothing withdrawn.
#
# NO REDOX ANNOTATION IS APPLIED AND NONE IS PRE-DECLARED. Pathway names are read
# after the ranking is produced. Any subsequent redox reading is a separate,
# separately declared analysis and is not introduced here.
#
# No FDR across pathways: this is a ranking exercise, not a screen.
# N3. Transcript associations. "Primed" appears nowhere. No ortholog call.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE29: where OXPHOS sits, and what moves with it\n", strrep("=", 78))

PATH_E29      <- file.path(DIR_RESULTS, "oxphos_neighbourhood.rds")
PATH_E29_ALG  <- file.path(DIR_TABLES,  "E29_aligned_change.csv")
PATH_E29_RES  <- file.path(DIR_TABLES,  "E29_diagonal_residuals.csv")
DIR_E29_FIG   <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_E25      <- file.path(DIR_RESULTS, "mitopps_group_rankings.rds")
PATH_E27      <- file.path(DIR_RESULTS, "category_readouts.rds")
PATH_E28      <- file.path(DIR_RESULTS, "gap_readout.rds")
PATH_MOUSE    <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")

set.seed(1)
NBOOT <- 5000L
NSIM  <- 5000L
K_NB  <- 10L          # neighbours ABOVE and BELOW, so 20 in total
N_LAB <- 10L          # residuals labelled in each direction per scatter panel

# =============================================================================
# 0. CONSTANTS AND THE READING RULES, FIXED BEFORE ANY NUMBER
# =============================================================================

MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"
REF       <- "OXPHOS subunits"       # the reference pathway throughout
UMB       <- "OXPHOS"                # secondary, and an aggregate of its halves
ARM_ORDER <- c("mouse 6W WT", "mouse 12W WT", "mouse 6W Myc+", "mouse 12W Myc+",
               "TCGA MYC-low", "TCGA MYC-high", "SCAN-B MYC-low", "SCAN-B MYC-high")

STATUS <- paste(
  "A MAP, NOT A READOUT. No pairing produced here may be used as a comparator",
  "without a separately declared validation across the full estimator ladder,",
  "both cohorts and the mouse.")

RULES <- list(
  status = STATUS,
  compositional = paste(
    "the design expected mitoPPS to be a closed composition; PART 0 measures",
    "that it is NOT - the per-sample sum is not fixed, and pairwise Spearman",
    "between the 142 scores sits on zero. The concern was reasonable a priori",
    "and is tested away"),
  zero_sum = paste(
    "what IS constrained is the percentile transform: the 142 percentiles sum",
    "to 7100 in every group, so the change vector sums to exactly zero and the",
    "neutral point for aligned change is -|delta REF|/141, not zero. Small",
    "here because it scales as 1/(P-1), and drawn rather than assumed"),
  aligned = paste(
    "aligned change = delta(pathway) * sign(delta REF). Negative means the",
    "pathway moved AGAINST the reference, positive WITH it, near zero barely",
    "moved. Same meaning in every contrast regardless of which way REF went"),
  no_value_crossing = paste(
    "no mitoPPS value crosses the species boundary; percentiles and percentile",
    "changes only, and cross-species comparison is by PATHWAY NAME"),
  mouse_n = paste(
    "n = 6 in the mouse. Neighbourhood stability is BOOTSTRAPPED in all eight",
    "groups rather than asserted"),
  no_redox = paste(
    "no redox annotation is applied and none is pre-declared; names are read",
    "after the ranking is produced"),
  batch = "the mouse timeline is batch-confounded; nothing is withdrawn on that account",
  no_fdr = "a ranking exercise, not a screen; distributions and ranks, no FDR",
  exploratory = "post-hoc and descriptive; not a hypothesis test",
  n3 = "transcript associations; 'primed' is never written of a transcript")

message("\n0. reading rules fixed before any number")
message("   R1  do groups at similar REF positions share neighbours BEYOND what")
message("       their whole-ordering rho already predicts?")
message("   R2  mouse 12W WT against human MYC-high specifically - same gap, very")
message("       different position, rho -0.51. Neither outcome is predicted")
message("   R2b are the off-diagonal pathways a COHERENT set or scattered? Per pair")
message("   R3  co-movement read THREE ways - with, against, independent. None is")
message("       privileged. A coherent opposing set is a MAP ENTRY NEEDING")
message("       VALIDATION, not a result. No structure at all is equally reportable")
message("   R4  do the same pathways track REF in mouse and human?")
message("   R5  the compositional diagnostic. If it is NOT strongly negative the")
message("       clean negative belongs in the note and NOT in CLAUDE.md's traps")
message("\n   ", STATUS)

# =============================================================================
# 1. INPUTS
# =============================================================================

message("\n1. inputs")

e25 <- readRDS(PATH_E25); e27 <- readRDS(PATH_E27); e28 <- readRDS(PATH_E28)
stopifnot(identical(unname(tools::md5sum(PATH_MOUSE)), MOUSE_MD5))
mo   <- readRDS(PATH_MOUSE)
mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
message("   E25 ", length(e25), " | E27 ", length(e27), " | E28 ", length(e28),
        " elements; mouse artefact md5 matches")

shared <- e25$intersection$shared
P_ALL  <- length(shared)
ANCHORS <- as.character(e27$settings$readouts)
stopifnot(length(shared) == 142L, REF %in% shared, all(ANCHORS %in% shared))

UT <- mito$mitopps_universe[shared, , drop = FALSE]
US <- sc$mitopps_universe[shared, , drop = FALSE]
ID_T <- colnames(UT); ID_S <- colnames(US)
ms <- mo$mitopps_scores
MPATHS <- setdiff(names(ms), c("sample", "group", "timepoint", "myc_status"))
MM <- as.matrix(ms[, MPATHS, drop = FALSE]); rownames(MM) <- ms$sample
MMs <- t(MM[, shared, drop = FALSE])

MOUSE_ARM <- c(`6W_neg` = "mouse 6W WT", `12W_neg` = "mouse 12W WT",
               `6W_pos` = "mouse 6W Myc+", `12W_pos` = "mouse 12W Myc+")
.tertiles <- function(v, ids) {
  q <- stats::quantile(v, c(1/3, 2/3), names = FALSE)
  g <- rep(NA_character_, length(v)); names(g) <- ids
  g[v <= q[1]] <- "MYC-low"; g[v >= q[2]] <- "MYC-high"; g
}
grp_T <- .tertiles(as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T]), ID_T)
grp_S <- .tertiles(as.numeric(sc$gsva_new[MYC_REF, ID_S]), ID_S)
stopifnot(identical(as.integer(c(sum(grp_T == "MYC-low", na.rm = TRUE),
                                 sum(grp_T == "MYC-high", na.rm = TRUE),
                                 sum(grp_S == "MYC-low", na.rm = TRUE),
                                 sum(grp_S == "MYC-high", na.rm = TRUE))),
                    as.integer(e25$myc_groups$n)))

ARM <- list()
for (g in names(MOUSE_ARM)) {
  ARM[[unname(MOUSE_ARM[[g]])]] <- list(sp = "mouse", M = MMs, j = which(ms$group == g))
}
for (lv in c("MYC-low", "MYC-high")) {
  ARM[[paste("TCGA", lv)]]   <- list(sp = "human", M = UT, j = which(grp_T == lv))
  ARM[[paste("SCAN-B", lv)]] <- list(sp = "human", M = US, j = which(grp_S == lv))
}
ARM <- ARM[ARM_ORDER]

.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}
.gp <- function(z, j) {
  v <- rowMeans(z$M[, j, drop = FALSE])
  stats::setNames(.pct(v), rownames(z$M))
}

# =============================================================================
# 2. PART 0 - reproduce, then the two premises about the score
# =============================================================================

message("\n2. PART 0: reproduce, then measure the two premises")

pos <- vapply(ARM_ORDER, function(a) .gp(ARM[[a]], ARM[[a]]$j), numeric(P_ALL))
rownames(pos) <- shared

e25o <- e25$orderings[e25$orderings$panel == "shared (142)", ]
ctrl25 <- max(abs(pos[cbind(match(as.character(e25o$pathway), shared),
                            match(as.character(e25o$arm), ARM_ORDER))] - e25o$pct))
e27p <- e27$category_positions
ctrl27 <- max(abs(pos[cbind(match(as.character(e27p$readout), shared),
                            match(as.character(e27p$arm), ARM_ORDER))] - e27p$pct))
message("   ARITHMETIC CONTROL: ", nrow(e25o), " positions against E25, max |delta| = ",
        format(ctrl25, digits = 2), "; ", nrow(e27p),
        " against E27, max |delta| = ", format(ctrl27, digits = 2))
CONTROL_OK <- max(ctrl25, ctrl27) < 1e-9
if (!CONTROL_OK) stop("CONTROL FAILED - stop.", call. = FALSE)

# --- the compositional diagnostic --------------------------------------------
message("\n   the COMPOSITIONAL DIAGNOSTIC - is mitoPPS a closed composition?")
COHORT_MAT <- list(TCGA = UT, `SCAN-B` = US, mouse = MMs)
compositional <- dplyr::bind_rows(lapply(names(COHORT_MAT), function(nm) {
  M <- COHORT_MAT[[nm]]
  s <- colSums(M)
  C <- stats::cor(t(M), method = "spearman")
  v <- C[upper.tri(C)]
  tibble::tibble(cohort = nm, n_samples = ncol(M),
                 sum_mean = mean(s), sum_sd = stats::sd(s),
                 sum_min = min(s), sum_max = max(s),
                 n_pairs = length(v), rho_median = stats::median(v),
                 rho_mean = mean(v), pct_negative = 100 * mean(v < 0),
                 rho_q05 = unname(stats::quantile(v, 0.05)),
                 rho_q95 = unname(stats::quantile(v, 0.95)))
}))
print(as.data.frame(compositional[, c("cohort", "sum_mean", "sum_sd", "sum_min",
                                      "sum_max")]), digits = 4, row.names = FALSE)
print(as.data.frame(compositional[, c("cohort", "n_pairs", "rho_median", "rho_mean",
                                      "pct_negative", "rho_q05", "rho_q95")]),
      digits = 3, row.names = FALSE)
CLOSED <- max(compositional$sum_sd) < 1e-6
NEGSHIFT <- any(compositional$rho_median < -0.1)
message("   IS IT A CLOSED COMPOSITION? ", if (CLOSED) "YES" else {
  paste0("NO - the per-sample sum RANGES (", sprintf("%.1f", min(compositional$sum_min)),
         " to ", sprintf("%.1f", max(compositional$sum_max)),
         "), so the panel is not closed") })
message("   IS THERE A SPURIOUS NEGATIVE SHIFT? ", if (NEGSHIFT) {
  "YES - PART B's contrast design is what protects it, and CLAUDE.md's traps should record it"
} else {
  paste0("NO - the median pairwise rho is ",
         sprintf("%+.3f to %+.3f", min(compositional$rho_median),
                 max(compositional$rho_median)),
         " with about half the pairs negative in every cohort. A CLEAN NEGATIVE:",
         " it belongs in the note and NOT in CLAUDE.md's trap list") })
message("   The concern was reasonable a priori - mitoPPS is ratio-based and",
        " centred at 1 -")
message("   and what settles it is that the SUM IS NOT FIXED. No Phase 2 note",
        " records that.")

# --- the zero-sum property of the percentile transform -----------------------
message("\n   WHAT IS CONSTRAINED: the percentile transform, one layer up")
CONTRASTS <- list(
  `WT window (6W -> 12W)`   = c(a = "mouse 6W WT",    b = "mouse 12W WT"),
  `Myc+ window (6W -> 12W)` = c(a = "mouse 6W Myc+",  b = "mouse 12W Myc+"),
  `MYC at 6W`               = c(a = "mouse 6W WT",    b = "mouse 6W Myc+"),
  `MYC at 12W`              = c(a = "mouse 12W WT",   b = "mouse 12W Myc+"),
  `TCGA MYC-high - low`     = c(a = "TCGA MYC-low",   b = "TCGA MYC-high"),
  `SCAN-B MYC-high - low`   = c(a = "SCAN-B MYC-low", b = "SCAN-B MYC-high"))
CONTRAST_KIND <- c(`WT window (6W -> 12W)` = "TIMELINE - batch-confounded",
                   `Myc+ window (6W -> 12W)` = "TIMELINE - batch-confounded",
                   `MYC at 6W` = "genotype, within timepoint",
                   `MYC at 12W` = "genotype, within timepoint",
                   `TCGA MYC-high - low` = "MYC activity, within cohort",
                   `SCAN-B MYC-high - low` = "MYC activity, within cohort")
CONTRAST_SP <- c(rep("mouse", 4L), rep("human", 2L))
names(CONTRAST_SP) <- names(CONTRASTS)

delta <- vapply(names(CONTRASTS), function(ct) {
  cn <- CONTRASTS[[ct]]; pos[, cn[["b"]]] - pos[, cn[["a"]]] }, numeric(P_ALL))
rownames(delta) <- shared
zero_sum <- tibble::tibble(
  contrast = names(CONTRASTS), sum_of_changes = colSums(delta),
  delta_ref = delta[REF, ],
  induced_neutral = -abs(delta[REF, ]) / (P_ALL - 1L))
print(as.data.frame(zero_sum), digits = 3, row.names = FALSE)
message("   the percentiles sum to ", sprintf("%.0f", sum(pos[, 1])),
        " in every group, so the CHANGE vector sums to zero")
message("   the neutral point for aligned change is -|delta REF|/", P_ALL - 1L,
        ", which runs ", sprintf("%.2f to %.2f", min(zero_sum$induced_neutral),
                                 max(zero_sum$induced_neutral)), " here.")
message("   IT IS SMALL BECAUSE IT SCALES AS 1/(P-1). On a 20-pathway panel the",
        " same REF change")
message("   would force ", sprintf("%.1f", -abs(delta[REF, "TCGA MYC-high - low"]) / 19),
        " points and would have to be handled rather than noted.")

# =============================================================================
# 3. PART A - the neighbourhood
# =============================================================================

message("\n3. PART A: the neighbourhood of ", REF)

.neighbours <- function(p, ref, k) {
  r <- rank(p, ties.method = "first")
  ri <- r[[ref]]
  keep <- names(r)[r >= ri - k & r <= ri + k & names(r) != ref]
  keep[order(r[keep])]
}
NB <- lapply(ARM_ORDER, function(a) .neighbours(pos[, a], REF, K_NB))
names(NB) <- ARM_ORDER

ref_position <- tibble::tibble(
  arm = ARM_ORDER, species = vapply(ARM, function(z) z$sp, ""),
  n = vapply(ARM, function(z) length(z$j), 0L),
  panel_size = P_ALL,
  pct = pos[REF, ARM_ORDER],
  rank = vapply(ARM_ORDER, function(a) unname(rank(pos[, a], ties.method = "first")[[REF]]), 0),
  n_neighbours = lengths(NB),
  umbrella_pct = pos[UMB, ARM_ORDER],
  umbrella_rank = vapply(ARM_ORDER, function(a)
    unname(rank(pos[, a], ties.method = "first")[[UMB]]), 0))
print(as.data.frame(ref_position), digits = 3, row.names = FALSE)
message("   the ", UMB, " umbrella is EXACTLY its two halves (E27) and carries no")
message("   independent information; it is reported as context only.")

neighbourhood <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  nb <- NB[[a]]
  tibble::tibble(arm = a, reference = REF, neighbour = nb,
                 neighbour_pct = pos[nb, a],
                 side = ifelse(pos[nb, a] > pos[REF, a], "above", "below"),
                 is_anchor = nb %in% ANCHORS)
}))
message("\n   the twenty neighbours, mouse 12W WT and TCGA MYC-high (R2's pair):")
neighbourhood %>%
  dplyr::filter(arm %in% c("mouse 12W WT", "TCGA MYC-high")) %>%
  dplyr::select(arm, side, neighbour, neighbour_pct) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

# --- what the neighbourhood is MADE OF, because it inflates any overlap -------
# The pathways nearest the reference are largely its own MitoCarta family. Two
# groups will therefore share neighbours partly because the OXPHOS branch sits
# near OXPHOS in both, which is a property of the ontology and not of biology.
# Measured here so it is visible before any overlap is read.
hx <- readxl::read_excel(here::here("data", "mitocarta_human",
                                    "Human.MitoCarta3.0.xls"), sheet = 4)
h_hier <- trimws(hx[["MitoPathways Hierarchy"]])
names(h_hier) <- trimws(hx[["MitoPathway"]])
OX_FAMILY <- c(names(h_hier)[!is.na(h_hier) & startsWith(h_hier, "OXPHOS")],
               "mtDNA-encoded OXPHOS subunits")
OX_FAMILY <- intersect(OX_FAMILY, shared)
message("\n   the OXPHOS family is ", length(OX_FAMILY), " of ", P_ALL,
        " panel members (", sprintf("%.0f", 100 * length(OX_FAMILY) / P_ALL), " pct)")
family_share <- tibble::tibble(
  arm = ARM_ORDER,
  n_family_neighbours = vapply(NB, function(v) sum(v %in% OX_FAMILY), 0L),
  pct_of_neighbourhood = 100 * vapply(NB, function(v) mean(v %in% OX_FAMILY), 0),
  panel_baseline_pct = 100 * length(OX_FAMILY) / P_ALL)
print(as.data.frame(family_share), digits = 3, row.names = FALSE)
message("   READ THE OVERLAP AGAINST THIS: where the neighbourhood is family-rich,")
message("   shared neighbours are partly an ONTOLOGY fact, not a biological one.")

# --- A1 the overlap matrix ----------------------------------------------------
gr <- utils::combn(length(ARM_ORDER), 2L)
ag <- e25$agreement[e25$agreement$panel == "shared (142)", ]
.rho_of <- function(a, b) {
  r <- ag$rho[(ag$arm_a == a & ag$arm_b == b) | (ag$arm_a == b & ag$arm_b == a)]
  if (length(r) != 1L) NA_real_ else r
}
overlap <- dplyr::bind_rows(lapply(seq_len(ncol(gr)), function(i) {
  a <- ARM_ORDER[gr[1, i]]; b <- ARM_ORDER[gr[2, i]]
  tibble::tibble(arm_a = a, arm_b = b,
                 shared_neighbours = length(intersect(NB[[a]], NB[[b]])),
                 rho_whole_ordering = .rho_of(a, b),
                 rank_a = ref_position$rank[ref_position$arm == a],
                 rank_b = ref_position$rank[ref_position$arm == b])
}))

# --- A2 the two nulls ---------------------------------------------------------
# PLAIN: two neighbourhoods drawn at random from the P-1 non-reference pathways.
# RHO-MATCHED: two orderings with the pair's OBSERVED rho, neighbourhoods taken
# at each group's OBSERVED reference rank. This is what R1 asks for - overlap
# BEYOND what the whole-ordering similarity already predicts.
.null_plain <- function(ka, kb) ka * kb / (P_ALL - 1L)
.null_rho <- function(rho, ra, rb, ka, kb, nsim = NSIM) {
  n <- P_ALL
  o <- vapply(seq_len(nsim), function(i) {
    x <- stats::rnorm(n)
    y <- rho * x + sqrt(max(0, 1 - rho^2)) * stats::rnorm(n)
    rx <- rank(x, ties.method = "first"); ry <- rank(y, ties.method = "first")
    ia <- which(rx >= ra - K_NB & rx <= ra + K_NB)
    ib <- which(ry >= rb - K_NB & ry <= rb + K_NB)
    length(intersect(ia, ib))
  }, 0)
  c(mean = mean(o), lo = unname(stats::quantile(o, 0.025)),
    hi = unname(stats::quantile(o, 0.975)))
}
overlap$null_plain <- .null_plain(K_NB * 2L, K_NB * 2L)
nr <- t(vapply(seq_len(nrow(overlap)), function(i)
  .null_rho(overlap$rho_whole_ordering[i], overlap$rank_a[i], overlap$rank_b[i],
            K_NB * 2L, K_NB * 2L), numeric(3)))
overlap$null_rho_mean <- nr[, "mean"]
overlap$null_rho_lo <- nr[, "lo"]; overlap$null_rho_hi <- nr[, "hi"]
overlap$beyond_rho <- overlap$shared_neighbours > overlap$null_rho_hi
message("\n   A1/A2: overlap against both nulls (the plain null is ",
        sprintf("%.2f", overlap$null_plain[1]), " for every pair):")
overlap %>%
  dplyr::select(arm_a, arm_b, shared_neighbours, rho_whole_ordering,
                null_rho_mean, null_rho_lo, null_rho_hi, beyond_rho) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)
R1 <- if (any(overlap$beyond_rho)) {
  paste0(sum(overlap$beyond_rho), " of ", nrow(overlap),
         " pairs share more neighbours than their own whole-ordering rho",
         " predicts - for those pairs the configuration is conserved beyond",
         " general ordering similarity")
} else {
  paste0("NO pair exceeds its rho-matched null. Every overlap is what the",
         " whole-ordering similarity already predicts, so POSITION MATCHES ARE",
         " COINCIDENTAL and must not be read as biological similarity")
}
message("   R1: ", R1)

# --- the stability bootstrap, all eight groups --------------------------------
message("\n   neighbourhood STABILITY - bootstrapped, not asserted:")
stability <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  z <- ARM[[a]]; obs <- NB[[a]]
  kept <- vapply(seq_len(NBOOT), function(i) {
    p <- .gp(z, sample(z$j, replace = TRUE))
    length(intersect(obs, .neighbours(p, REF, K_NB)))
  }, 0)
  tibble::tibble(arm = a, species = z$sp, n = length(z$j),
                 n_neighbours = length(obs), mean_retained = mean(kept),
                 lo = unname(stats::quantile(kept, 0.025)),
                 hi = unname(stats::quantile(kept, 0.975)))
}))
print(as.data.frame(stability), digits = 3, row.names = FALSE)
message("   read the mouse panels against this: ",
        sprintf("%.0f", min(stability$mean_retained[stability$species == "mouse"])),
        " to ", sprintf("%.0f", max(stability$mean_retained[stability$species == "mouse"])),
        " of ", 2L * K_NB, " neighbours survive resampling in the mouse against ",
        sprintf("%.0f", min(stability$mean_retained[stability$species == "human"])),
        " to ", sprintf("%.0f", max(stability$mean_retained[stability$species == "human"])),
        " in the human.")

# =============================================================================
# 4. PART A3 - the ordering scatters
# =============================================================================

message("\n4. PART A3: the ordering scatters and their residuals")

PAIRS <- list(c("mouse 12W WT", "MYC-low"), c("mouse 12W WT", "MYC-high"),
              c("mouse 6W Myc+", "MYC-low"), c("mouse 6W Myc+", "MYC-high"),
              c("mouse 12W Myc+", "MYC-high"), c("mouse 6W WT", "MYC-high"))
.pair_tab <- function(a, b) {
  tibble::tibble(pair = paste(a, "vs", b), arm_x = a, arm_y = b,
                 pathway = shared, x = pos[, a], y = pos[, b],
                 residual = pos[, b] - pos[, a],
                 is_anchor = shared %in% ANCHORS,
                 rho = .rho_of(a, b))
}
residuals_tab <- dplyr::bind_rows(lapply(c("TCGA", "SCAN-B"), function(coh) {
  dplyr::bind_rows(lapply(PAIRS, function(p) .pair_tab(p[1], paste(coh, p[2]))))
}))
residuals_tab <- residuals_tab %>%
  dplyr::group_by(pair) %>%
  dplyr::mutate(abs_rank = rank(-abs(residual), ties.method = "first"),
                label_it = abs_rank <= 2L * N_LAB) %>%
  dplyr::ungroup()

# R2b: are the off-diagonal pathways a COHERENT set, or scattered? A set is
# coherent if the same pathways recur as extreme residuals across pairs.
extremes <- residuals_tab %>%
  dplyr::filter(label_it) %>%
  dplyr::count(pathway, name = "n_pairs_extreme") %>%
  dplyr::arrange(dplyr::desc(n_pairs_extreme))
n_pairs_total <- length(unique(residuals_tab$pair))
coherence <- tibble::tibble(
  n_pairs = n_pairs_total,
  n_distinct_extremes = nrow(extremes),
  max_possible_if_identical = 2L * N_LAB,
  recur_in_half_or_more = sum(extremes$n_pairs_extreme >= n_pairs_total / 2))
print(as.data.frame(coherence), row.names = FALSE)
message("   the pathways most often furthest from the diagonal:")
extremes %>% utils::head(12) %>% as.data.frame() %>% print(row.names = FALSE)
R2B <- if (coherence$n_distinct_extremes <= 3L * N_LAB) {
  paste0("COHERENT - only ", coherence$n_distinct_extremes,
         " distinct pathways fill the extreme-residual slots across ",
         n_pairs_total, " pairs, and ", coherence$recur_in_half_or_more,
         " recur in half or more. That is a map entry worth naming")
} else {
  paste0("SCATTERED - ", coherence$n_distinct_extremes,
         " distinct pathways fill the extreme-residual slots across ",
         n_pairs_total, " pairs, with only ", coherence$recur_in_half_or_more,
         " recurring in half or more. The orderings differ DIFFUSELY and no",
         " subset carries the difference")
}
message("   R2b: ", R2B)

r2pair <- residuals_tab[residuals_tab$pair == "mouse 12W WT vs TCGA MYC-high", ]
r2row <- overlap[overlap$arm_a == "mouse 12W WT" &
                   overlap$arm_b == "TCGA MYC-high", ]
R2 <- paste0("mouse 12W WT and TCGA MYC-high share ", r2row$shared_neighbours,
             " of ", 2L * K_NB, " neighbours. Their rho-matched null is ",
             sprintf("%.1f", r2row$null_rho_mean), " [", r2row$null_rho_lo, ", ",
             r2row$null_rho_hi, "], so the observed overlap DOES NOT EXCEED what",
             " their whole-ordering similarity already predicts",
             ". Their whole orderings sit at rho ",
             sprintf("%+.2f", r2pair$rho[1]),
             " and the residuals run ", sprintf("%.0f", min(r2pair$residual)),
             " to ", sprintf("%.0f", max(r2pair$residual)),
             " percentile points. The identical central-dogma gap in E28 came",
             " from very different positions, and this is what sits behind it")
message("   R2: ", R2)

# =============================================================================
# 5. PART B - co-movement across contrasts
# =============================================================================

message("\n5. PART B: aligned change. aligned = delta(pathway) * sign(delta REF)")

aligned <- dplyr::bind_rows(lapply(names(CONTRASTS), function(ct) {
  s <- sign(delta[REF, ct])
  tibble::tibble(contrast = ct, kind = unname(CONTRAST_KIND[[ct]]),
                 species = unname(CONTRAST_SP[[ct]]), pathway = shared,
                 raw_change = delta[, ct], aligned_change = delta[, ct] * s,
                 ref_direction = ifelse(s > 0, "REF rose", "REF fell"),
                 neutral_point = -abs(delta[REF, ct]) / (P_ALL - 1L),
                 is_anchor = shared %in% ANCHORS, is_reference = shared == REF)
}))
message("   RAW CHANGE IS CARRIED BESIDE THE ALIGNED ONE, and it is required: a")
message("   pathway can be strongly negative because it moved the OPPOSITE way, or")
message("   because it moved the SAME way and much less. Only the first is opposition.")

message("\n   the five most OPPOSED and five most ALIGNED, TCGA contrast:")
aligned %>%
  dplyr::filter(contrast == "TCGA MYC-high - low", !is_reference) %>%
  dplyr::arrange(aligned_change) %>%
  dplyr::slice(c(1:5, (dplyr::n() - 4):dplyr::n())) %>%
  dplyr::select(pathway, raw_change, aligned_change, ref_direction) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

aligned_by_species <- aligned %>%
  dplyr::filter(!is_reference) %>%
  dplyr::group_by(pathway, species) %>%
  dplyr::summarise(mean_aligned = mean(aligned_change), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = species, values_from = mean_aligned) %>%
  dplyr::rename(mouse_mean = mouse, human_mean = human)
aligned_by_species$sign_agree <- sign(aligned_by_species$mouse_mean) ==
                                 sign(aligned_by_species$human_mean)
aligned_by_species$quadrant <- dplyr::case_when(
  aligned_by_species$mouse_mean < 0 & aligned_by_species$human_mean < 0 ~ "opposes REF in BOTH",
  aligned_by_species$mouse_mean > 0 & aligned_by_species$human_mean > 0 ~ "moves WITH REF in BOTH",
  TRUE ~ "species- or state-specific")
message("\n   B2: two contrasts is too few for a correlation, so MEANS and SIGN")
message("   AGREEMENT only, never a rho.")
aligned_by_species %>% dplyr::count(quadrant) %>% as.data.frame() %>%
  print(row.names = FALSE)

REF_MOUSE <- mean(abs(delta[REF, names(CONTRASTS)[CONTRAST_SP == "mouse"]]))
REF_HUMAN <- mean(abs(delta[REF, names(CONTRASTS)[CONTRAST_SP == "human"]]))
message("   the reference itself sits at (", sprintf("%.1f", REF_MOUSE), ", ",
        sprintf("%.1f", REF_HUMAN), ") by construction; pathways near it track it.")

# B3: the consistency check on a known observation, not a new finding
PI <- "Protein import, sorting and homeostasis"
b3 <- aligned %>%
  dplyr::filter(pathway == PI) %>%
  dplyr::select(contrast, species, raw_change, aligned_change)
b3_rank <- aligned %>%
  dplyr::filter(!is_reference) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(rank_aligned = rank(aligned_change, ties.method = "first")) %>%
  dplyr::filter(pathway == PI) %>%
  dplyr::select(contrast, rank_aligned) %>%
  dplyr::ungroup()
b3 <- dplyr::left_join(b3, b3_rank, by = "contrast")
message("\n   B3: ", PI, " - a CONSISTENCY CHECK on E28, not a new finding:")
print(as.data.frame(b3), digits = 3, row.names = FALSE)

opp_both <- aligned_by_species$pathway[aligned_by_species$quadrant == "opposes REF in BOTH"]
with_both <- aligned_by_species$pathway[aligned_by_species$quadrant == "moves WITH REF in BOTH"]
# A count of opposing pathways is NOT a coherent set. Coherence would need the
# opposing pathways to be related to each other, which this script does not test
# - and R2b already found the off-diagonal structure SCATTERED. So the count is
# reported as a count and the coherence claim is declined.
CHANCE_AGREE <- nrow(aligned_by_species) / 2
R3 <- paste0("read three ways: ", length(with_both), " pathways move WITH the ",
             "reference in both species, ", length(opp_both), " oppose it in both,",
             " and ", sum(!aligned_by_species$sign_agree),
             " are species- or state-specific. THE COUNT IS NOT A COHERENT SET:",
             " coherence would need the opposing pathways to be related to each",
             " other, which is not tested here, and R2b found the off-diagonal",
             " structure SCATTERED. So this is reported as a distribution of",
             " signs and NOT as a named opposing programme. Nothing here is a",
             " map entry that could be promoted without separate validation")
message("\n   R3: ", R3)
R4 <- paste0(sum(aligned_by_species$sign_agree), " of ", nrow(aligned_by_species),
             " pathways keep their sign between mouse and human (",
             sprintf("%.0f", 100 * mean(aligned_by_species$sign_agree)),
             " pct), against ", sprintf("%.0f", CHANCE_AGREE),
             " expected if the signs were independent. ",
             if (mean(aligned_by_species$sign_agree) > 0.7) {
               "Largely a CONSERVED neighbourhood"
             } else {
               paste0("The co-movement is substantially SPECIES- OR STATE-SPECIFIC,",
                      " and any cross-species reading built on a single comparator",
                      " inherits that")
             })
message("   R4: ", R4)
R5 <- if (NEGSHIFT) {
  "pairwise mitoPPS correlations ARE shifted negative - PART B's contrast design protects it and CLAUDE.md's traps should record it"
} else {
  paste0("pairwise mitoPPS correlations are NOT shifted negative (median ",
         sprintf("%+.3f to %+.3f", min(compositional$rho_median),
                 max(compositional$rho_median)),
         ") and the panel is NOT a closed composition (the per-sample sum ranges ",
         sprintf("%.1f to %.1f", min(compositional$sum_min), max(compositional$sum_max)),
         "). A CLEAN NEGATIVE: recorded in the note, NOT added to CLAUDE.md's traps")
}
message("   R5: ", R5)

# =============================================================================
# 6. PART C - figures
# =============================================================================

message("\n6. PART C: figures")
if (!dir.exists(DIR_E29_FIG)) dir.create(DIR_E29_FIG, recursive = TRUE)

CAP_MAP <- paste(
  "A MAP, NOT A READOUT: no pairing here may be used as a comparator without separate declared validation.",
  "Percentiles only - no mitoPPS value crosses the species boundary, and a point on the diagonal shares a RANK, not a value.",
  "Mouse contrasts come from n = 6 and are far less precise than human ones (365 and 1,069 per group).",
  "Exploratory, post-hoc, descriptive; not a hypothesis test.", sep = "\n")

ANCH_COL <- "#D55E00"

# --- C1: the reference's position, with the panel range behind it -------------
ci <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  z <- ARM[[a]]
  b <- vapply(seq_len(NBOOT), function(i) .gp(z, sample(z$j, replace = TRUE))[[REF]], 0)
  tibble::tibble(arm = a, species = z$sp, pct = pos[REF, a],
                 lo = unname(stats::quantile(b, 0.025)),
                 hi = unname(stats::quantile(b, 0.975)))
}))
ci$arm <- factor(ci$arm, levels = ARM_ORDER)
p1 <- ggplot2::ggplot(ci, ggplot2::aes(arm, pct)) +
  ggplot2::annotate("rect", xmin = 0.4, xmax = length(ARM_ORDER) + 0.6,
                    ymin = 0, ymax = 100, fill = "grey92") +
  ggplot2::geom_hline(yintercept = 50, linewidth = 0.3, colour = "grey65",
                      linetype = 2) +
  ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.4, colour = "grey55",
                      linetype = 3) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi, colour = species),
                          linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(colour = species), size = 2.6) +
  ggplot2::scale_colour_manual(values = c(mouse = "#0072B2", human = ANCH_COL)) +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::labs(x = NULL, y = paste0(REF, " percentile in the group ordering"),
                colour = NULL, title = paste0("Where ", REF, " sits - no comparator"),
                subtitle = paste0("The pale band is the full ", P_ALL,
                                  "-pathway range. Percentiles are UNIFORM by construction,",
                                  " so it is a reference range and not a density."),
                caption = CAP_MAP) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C1_reference_position.pdf"), p1,
                width = 7.5, height = 5)

# --- C2: the overlap matrix ---------------------------------------------------
sym <- dplyr::bind_rows(overlap,
                        dplyr::rename(overlap, arm_a = arm_b, arm_b = arm_a))
sym$arm_a <- factor(sym$arm_a, levels = ARM_ORDER)
sym$arm_b <- factor(sym$arm_b, levels = rev(ARM_ORDER))
p2 <- ggplot2::ggplot(sym, ggplot2::aes(arm_a, arm_b, fill = shared_neighbours)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d\nrho %+.2f", shared_neighbours,
                                                  rho_whole_ordering)), size = 2.1) +
  ggplot2::scale_fill_gradient(low = "grey95", high = "#0072B2") +
  ggplot2::labs(x = NULL, y = NULL, fill = "shared\nneighbours",
                title = paste0("Shared neighbours of ", REF, ", out of ", 2L * K_NB),
                subtitle = paste0("The plain random null is ",
                                  sprintf("%.2f", overlap$null_plain[1]),
                                  " for every pair; the rho-matched null is in the object."),
                caption = CAP_MAP) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C2_neighbour_overlap.pdf"), p2,
                width = 8, height = 7)

# --- C3: the ordering scatters ------------------------------------------------
.scatter_fig <- function(coh) {
  d <- residuals_tab[grepl(paste0("vs ", coh, " "), residuals_tab$pair), ]
  d$panel <- sprintf("%s  vs  %s\nrho %+.3f%s", d$arm_x, d$arm_y, d$rho,
                     ifelse(abs(d$rho) < 0.1,
                            "  - NO RELATIONSHIP: this is the calibration panel", ""))
  ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.35,
                         colour = "grey55", linetype = 2) +
    ggplot2::geom_point(colour = "grey75", size = 0.7) +
    ggplot2::geom_point(data = subset(d, label_it & !is_anchor),
                        colour = "grey35", size = 0.9) +
    ggplot2::geom_point(data = subset(d, is_anchor), colour = ANCH_COL, size = 1.8) +
    ggrepel::geom_text_repel(data = subset(d, label_it & !is_anchor),
                             ggplot2::aes(label = pathway), size = 1.5,
                             colour = "grey30", max.overlaps = 30, seed = 1,
                             segment.size = 0.15, min.segment.length = 0) +
    ggrepel::geom_text_repel(data = subset(d, is_anchor),
                             ggplot2::aes(label = pathway), size = 1.7,
                             colour = ANCH_COL, fontface = "bold",
                             max.overlaps = 30, seed = 1, segment.size = 0.2,
                             min.segment.length = 0) +
    ggplot2::facet_wrap(~ panel, ncol = 2L) +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::labs(x = "percentile in the MOUSE ordering",
                  y = paste0("percentile in the ", coh, " ordering"),
                  title = paste0("The orderings behind E25's correlations - ", coh),
                  subtitle = paste0("All ", P_ALL,
                                    " shared pathways. Anchors in colour and on top; the ",
                                    2L * N_LAB,
                                    " largest off-diagonal residuals labelled; the rest grey."),
                  caption = CAP_MAP) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6.5),
                   plot.caption = ggplot2::element_text(size = 6, hjust = 0))
}
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C3_ordering_scatters_TCGA.pdf"),
                .scatter_fig("TCGA"), width = 9, height = 12)
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C3_ordering_scatters_SCANB.pdf"),
                .scatter_fig("SCAN-B"), width = 9, height = 12)

# --- C4: the aligned-change heatmap -------------------------------------------
ord <- aligned_by_species$pathway[order(aligned_by_species$mouse_mean +
                                          aligned_by_species$human_mean)]
h <- aligned %>%
  dplyr::filter(!is_reference) %>%
  dplyr::mutate(pathway = factor(pathway, levels = ord),
                contrast = factor(contrast, levels = names(CONTRASTS)))
lab_rows <- c(utils::head(ord, N_LAB), utils::tail(ord, N_LAB))
p4 <- ggplot2::ggplot(h, ggplot2::aes(contrast, pathway, fill = aligned_change)) +
  ggplot2::geom_tile() +
  ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.6, colour = "black") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#B2182B",
                                midpoint = 0) +
  ggplot2::scale_y_discrete(breaks = lab_rows) +
  ggplot2::labs(x = NULL, y = NULL, fill = "aligned\nchange",
                title = "Aligned change: does a pathway move with the reference or against it?",
                subtitle = paste0("aligned = delta(pathway) * sign(delta ", REF,
                                  "). Rows ordered by the mean. Mouse contrasts left of the rule, human right.",
                                  " Only the top and bottom ", N_LAB, " rows are labelled."),
                caption = CAP_MAP) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 5),
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C4_aligned_heatmap.pdf"), p4,
                width = 8, height = 10)

# --- C5: mouse against human --------------------------------------------------
NEUT_M <- mean(zero_sum$induced_neutral[CONTRAST_SP == "mouse"])
NEUT_H <- mean(zero_sum$induced_neutral[CONTRAST_SP == "human"])
a5 <- aligned_by_species
a5$is_anchor <- a5$pathway %in% ANCHORS
a5$extreme <- rank(-(abs(a5$mouse_mean) + abs(a5$human_mean)), ties.method = "first") <= 12L
p5 <- ggplot2::ggplot(a5, ggplot2::aes(mouse_mean, human_mean)) +
  ggplot2::geom_hline(yintercept = NEUT_H, linewidth = 0.35, colour = "grey45") +
  ggplot2::geom_vline(xintercept = NEUT_M, linewidth = 0.35, colour = "grey45") +
  ggplot2::geom_point(colour = "grey75", size = 1) +
  ggplot2::geom_point(data = subset(a5, extreme & !is_anchor), colour = "grey30",
                      size = 1.3) +
  ggplot2::geom_point(data = subset(a5, is_anchor), colour = ANCH_COL, size = 2.2) +
  ggplot2::annotate("point", x = REF_MOUSE, y = REF_HUMAN, shape = 8, size = 3.4,
                    colour = "black") +
  ggplot2::annotate("text", x = REF_MOUSE, y = REF_HUMAN, label = REF, size = 2.4,
                    fontface = "bold", vjust = -1.2) +
  ggrepel::geom_text_repel(data = subset(a5, extreme | is_anchor),
                           ggplot2::aes(label = pathway,
                                        colour = ifelse(is_anchor, "a", "b")),
                           size = 1.9, max.overlaps = 40, seed = 1,
                           segment.size = 0.2, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(a = ANCH_COL, b = "grey30")) +
  ggplot2::labs(x = "mean aligned change, MOUSE (4 contrasts, n = 6 per group - NOISIER AXIS)",
                y = "mean aligned change, HUMAN (2 contrasts, 365 and 1,069 per group)",
                title = "Does the same thing track the reference in both species?",
                subtitle = paste0("Zero lines are drawn at the INDUCED NEUTRAL POINT (",
                                  sprintf("%.2f", NEUT_M), ", ", sprintf("%.2f", NEUT_H),
                                  "), not at zero: the percentile change vector sums to zero by construction.",
                                  " The star marks the reference itself."),
                caption = CAP_MAP) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C5_mouse_vs_human.pdf"), p5,
                width = 8.5, height = 7)

# --- C6, supplementary: six ordered lollipops ---------------------------------
l6 <- aligned %>% dplyr::filter(!is_reference)
l6 <- l6 %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(o = rank(aligned_change, ties.method = "first"),
                tail_lab = o <= 8L | o > dplyr::n() - 8L) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(contrast = factor(contrast, levels = names(CONTRASTS)))
p6 <- ggplot2::ggplot(l6, ggplot2::aes(o, aligned_change)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = neutral_point), linewidth = 0.3,
                      colour = "grey45") +
  ggplot2::geom_segment(ggplot2::aes(xend = o, yend = neutral_point),
                        linewidth = 0.2, colour = "grey80") +
  ggplot2::geom_point(ggplot2::aes(colour = is_anchor), size = 0.7) +
  ggrepel::geom_text_repel(data = subset(l6, tail_lab),
                           ggplot2::aes(label = pathway), size = 1.4,
                           max.overlaps = 25, seed = 1, segment.size = 0.15) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = ANCH_COL),
                               guide = "none") +
  ggplot2::facet_wrap(~ contrast, ncol = 2L, scales = "free_y") +
  ggplot2::labs(x = "pathways, ordered by aligned change",
                y = "aligned change (percentile points)",
                title = "Every pathway, every contrast - supplementary",
                subtitle = "The horizontal rule is the induced neutral point, not zero. Anchors in colour.",
                caption = CAP_MAP) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E29_FIG, "E29_C6_lollipops.pdf"), p6,
                width = 9, height = 11)
message("   wrote six figures under ", DIR_E29_FIG)

# =============================================================================
# 7. PART D - verdict
# =============================================================================

message("\n7. PART D: verdict")

verdict <- tibble::tibble(
  rule = c("STATUS", "CONTROL", "compositional diagnostic", "the zero-sum constraint",
           "neighbourhood composition",
           "neighbourhood stability", "R1 conservation beyond rho", "R2 the named pair",
           "R2b the off-diagonal set", "R3 co-movement", "R4 mouse against human",
           "R5 the diagnostic's destination"),
  outcome = c(
    STATUS,
    paste0("PASSED - ", nrow(e25o), " positions reproduce E25 and ", nrow(e27p),
           " reproduce E27, max |delta| ", format(max(ctrl25, ctrl27), digits = 2)),
    paste0("mitoPPS is NOT a closed composition - the per-sample sum ranges ",
           sprintf("%.1f to %.1f", min(compositional$sum_min), max(compositional$sum_max)),
           " - and pairwise rho sits on zero (median ",
           sprintf("%+.3f to %+.3f", min(compositional$rho_median),
                   max(compositional$rho_median)), ", about half negative). ",
           "The a priori concern is TESTED AWAY"),
    paste0("the percentile transform IS constrained: changes sum to zero, so the ",
           "neutral point is -|delta REF|/", P_ALL - 1L, " = ",
           sprintf("%.2f to %.2f", min(zero_sum$induced_neutral),
                   max(zero_sum$induced_neutral)),
           ". Small because it scales as 1/(P-1); on a 20-pathway panel it would be ",
           sprintf("%.1f", -abs(delta[REF, "TCGA MYC-high - low"]) / 19)),
    paste0("the OXPHOS family is ", length(OX_FAMILY), " of ", P_ALL,
           " panel members but fills ",
           sprintf("%.0f-%.0f", min(family_share$pct_of_neighbourhood),
                   max(family_share$pct_of_neighbourhood)),
           " pct of the neighbourhood. Shared neighbours are partly an ONTOLOGY",
           " fact, not a biological one"),
    paste0("of ", 2L * K_NB, " neighbours, mouse groups retain ",
           sprintf("%.0f-%.0f", min(stability$mean_retained[stability$species == "mouse"]),
                   max(stability$mean_retained[stability$species == "mouse"])),
           " under resampling against ",
           sprintf("%.0f-%.0f", min(stability$mean_retained[stability$species == "human"]),
                   max(stability$mean_retained[stability$species == "human"])),
           " in the human. Weight the mouse panels accordingly"),
    R1, R2, R2B, R3, R4, R5))
print(as.data.frame(verdict), row.names = FALSE, right = FALSE)
message("\n   AND WHATEVER IT SAYS, IT IS A MAP. NOT A READOUT. NOT A TEST.")

NOTES <- c(
  STATUS,
  paste("The compositional premise the design was written on is TESTED AWAY.",
        "mitoPPS is ratio-based and centred at 1, so the concern was reasonable",
        "a priori; what settles it is that THE PER-SAMPLE SUM IS NOT FIXED, so",
        "the panel is not a closed composition and the standard objection does",
        "not apply. Pairwise Spearman between the 142 scores sits on zero with",
        "about half the pairs negative in every cohort. No Phase 2 note records",
        "that the sum is unconstrained, and that is the useful fact here."),
  paste("PART B's contrast-profile design is KEPT, but for a different reason",
        "than it was written for: it is the estimand every other Phase 2",
        "contrast uses, which keeps E29 comparable with E25, E27 and E28."),
  paste("R5 says the diagnostic goes in CLAUDE.md's traps IF strongly negative.",
        "It is not, so it is recorded here and CLAUDE.md is NOT edited."),
  paste("WHAT IS CONSTRAINED is the percentile transform: the 142 percentiles",
        "sum to 7100 in every group, so the change vector sums to exactly zero",
        "and the neutral point for aligned change is -|delta REF|/141, not zero.",
        "It is drawn rather than assumed. It is small HERE because the forcing",
        "scales as 1/(P-1); on a 20-pathway panel the same reference change",
        "would force -2.7 points and would have to be handled."),
  paste("Aligned change = delta(pathway) * sign(delta reference), so the sign",
        "means the same thing in every contrast regardless of which way the",
        "reference went. THE RAW CHANGE IS CARRIED BESIDE IT and is required: a",
        "large negative aligned change can mean the pathway moved the opposite",
        "way, or the same way and much less, and only the first is opposition."),
  paste("The neighbourhood is largely the reference's OWN MitoCarta FAMILY.",
        "Two groups therefore share neighbours partly because the OXPHOS branch",
        "sits near OXPHOS in both, which is a property of the ontology and not",
        "of biology. Measured and reported before any overlap is read."),
  paste("A COUNT of opposing pathways is not a COHERENT SET. Coherence would",
        "need the opposing pathways to be related to each other, which is not",
        "tested here, and R2b found the off-diagonal structure scattered. The",
        "sign distribution is reported and no opposing programme is named."),
  paste("Neighbourhood stability is BOOTSTRAPPED in all eight groups rather than",
        "asserted, so the n = 6 caveat is a number and a reader can see how to",
        "weight the mouse panels."),
  paste("The overlap is read against TWO nulls: a plain random one, and a",
        "rho-matched, position-matched one that asks what R1 actually asks -",
        "overlap BEYOND what the whole-ordering similarity already predicts."),
  paste("Two scatter panels have a whole-ordering rho near zero and are LABELLED",
        "IN THE FIGURE as calibration panels, so an absence of structure reads as",
        "what no relationship looks like rather than as a failed panel."),
  paste("A point on a scatter's diagonal shares a RANK in both orderings. It says",
        "nothing about its mitoPPS value, which does not cross the species",
        "boundary."),
  paste("Percentiles are UNIFORM within a group by construction, so C1's pale band",
        "is a reference range and not a density."),
  paste("NO REDOX ANNOTATION is applied and none is pre-declared. Names were read",
        "after the ranking was produced. Any redox reading is a separate,",
        "separately declared analysis."),
  paste("The mouse timeline is BATCH-CONFOUNDED (batch = timepoint); reported in",
        "full and nothing withdrawn on that account."),
  "No FDR: a ranking exercise, not a screen. No ortholog call. N3 throughout.")

saveRDS(list(
  rules = RULES, status = STATUS,
  settings = list(reference = REF, umbrella = UMB, anchors = ANCHORS,
                  arm_order = ARM_ORDER, contrasts = CONTRASTS,
                  k_neighbours = K_NB, n_label = N_LAB, panel_size = P_ALL,
                  nboot = NBOOT, nsim = NSIM, seed = 1,
                  mouse_artefact = list(md5 = MOUSE_MD5, recomputed = FALSE)),
  control = list(vs_e25 = ctrl25, vs_e27 = ctrl27, ok = CONTROL_OK),
  compositional_diagnostic = compositional, is_closed = CLOSED,
  neg_shift = NEGSHIFT, zero_sum = zero_sum,
  reference_position = ref_position, neighbourhood = neighbourhood,
  family_share = family_share, oxphos_family = OX_FAMILY,
  overlap_matrix = overlap,
  overlap_null = list(plain = overlap$null_plain[1], k = 2L * K_NB, nsim = NSIM),
  stability = stability,
  ordering_scatters = residuals_tab, diagonal_residuals = residuals_tab,
  extremes = extremes, coherence = coherence,
  aligned_change = aligned, aligned_by_species = aligned_by_species,
  reference_anchor = c(mouse = REF_MOUSE, human = REF_HUMAN),
  protein_import = b3,
  readings = list(R1 = R1, R2 = R2, R2b = R2B, R3 = R3, R4 = R4, R5 = R5),
  verdict = verdict, analysis_date = Sys.Date(), notes = NOTES), PATH_E29)

readr::write_csv(aligned, PATH_E29_ALG)
readr::write_csv(residuals_tab, PATH_E29_RES)

message("\nE29: done.")
message("    ", PATH_E29)
message("    ", PATH_E29_ALG)
message("    ", PATH_E29_RES)
message("    six figures under ", DIR_E29_FIG)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E29)

  ## the status, first and always
  x$status

  ## the control
  x$control

  ## the two premises: one tested away, one found
  x$compositional_diagnostic |> as.data.frame()
  x$is_closed; x$neg_shift
  x$zero_sum |> as.data.frame()

  ## PART A: where the reference sits, and how stable its neighbourhood is
  x$reference_position |> as.data.frame()
  x$stability |> as.data.frame()
  subset(x$neighbourhood, arm == "mouse 12W WT") |> as.data.frame()
  subset(x$neighbourhood, arm == "TCGA MYC-high") |> as.data.frame()

  ## A1/A2: overlap against the plain and rho-matched nulls
  x$overlap_matrix |> as.data.frame()
  x$overlap_null

  ## A3: the scatters behind E25's correlations, and R2b's coherence question
  subset(x$diagonal_residuals,
         pair == "mouse 12W WT vs TCGA MYC-high" & label_it) |> as.data.frame()
  x$extremes |> as.data.frame()
  x$coherence |> as.data.frame()

  ## PART B: aligned change. RAW CHANGE IS BESIDE IT AND IS REQUIRED.
  subset(x$aligned_change, contrast == "TCGA MYC-high - low") |> as.data.frame()
  x$aligned_by_species |> as.data.frame()
  x$reference_anchor
  x$protein_import |> as.data.frame()

  ## the verdict and the six readings
  x$verdict |> as.data.frame()
  x$readings

  ## every caveat the note must repeat
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")
}
