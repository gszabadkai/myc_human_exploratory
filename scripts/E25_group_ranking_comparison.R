# scripts/E25_group_ranking_comparison.R
# =============================================================================
# E25 - the mitochondrial pathway ORDERING, group by group, mouse and human
# =============================================================================
#
# THE QUESTION, as the author set it: compare the RANKINGS in the four mouse
# groups (6W WT, 12W WT, 6W Myc+, 12W Myc+) and in MYC-high against MYC-low
# human tumours, in each cohort separately. Eight orderings of the same shared
# panel, compared as orderings.
#
# WHY THIS REPLACES E24's ESTIMAND, AND IT IS NOT A REFINEMENT
# ------------------------------------------------------------
# E24 ranked ONE pathway inside each SAMPLE and summarised the percentile by
# group. The mouse figures (fig02_reallocation_ranked.R,
# figS6_reallocation_temporal_ranked.R) do something else: they rank PATHWAYS by
# a GROUP CONTRAST in mitoPPS. The two statistics behave differently, and not
# symmetrically between the species. Measured, before this script was written:
#
#   rho(mean within-sample percentile, log set size)   mouse +0.165
#                                                      TCGA  +0.536
#                                                      SCAN-B +0.571
#   rho(group contrast, log set size)                  mouse -0.14 to +0.03
#
# So E24 compared two differently size-biased statistics across the species
# boundary. The mechanism is not mysterious: mitoPPS is centred at EXACTLY 1 for
# every pathway by construction (verified: colMeans of the mouse object are
# 1.000000 for all 144), so the first moment is size-free in both species and
# both implementations take the MEAN per set (mouse script 08 line 386, human
# .path_scores). What is NOT size-free is the cross-sample VARIANCE -
# rho(sd, log size) = -0.509 in the mouse - and because pairwise ratios are
# right-skewed, a rank transform converts that variance difference into a
# mean-rank difference. It scales with cohort heterogeneity, which is why human
# tumours show three times the mouse effect. Ranking pathways by a GROUP mean or
# a GROUP contrast does not do that.
#
# THE BATCH CONFOUNDER, AND THE AUTHOR'S DECISION ON IT
# -----------------------------------------------------
# The mouse time axis is batch-confounded: batch = timepoint. The mouse repo's
# own figS6 header says the individual 6W->12W changes are "DESCRIPTIVE, not
# claims" and that the clean quantity is their divergence; its evidence audit
# says the batch axis is genotype-independent and so cancels inside each
# timepoint's contrast.
#
# THE AUTHOR HAS RULED that the timeline comparisons STAY IN: reprioritisation
# is a timeline effect and DESeq2 normalisation is what handles the batch. This
# script therefore computes and reports every timeline comparison in full. The
# confounder is CARRIED as a stated caveat on the temporal quantities - in the
# saved object, in the verdict table and on the figure - and it is never used to
# withdraw a number. Genotype contrasts within a timepoint remain the cleaner
# axis and are labelled as such, not preferred silently.
#
# RANK, NEVER VALUE
# -----------------
# No mitoPPS value is compared across the species boundary anywhere. Every
# cross-species quantity in this script is an ordering, a rank correlation
# between orderings, or a percentile position within an ordering. Species is a
# cohort (CLAUDE.md trap 6). Mouse and human are never scored together, and TCGA
# and SCAN-B are never averaged.
#
# RE-RANKING, NEVER RE-SCORING
# ----------------------------
# The hierarchy-restricted and mtDNA-dropped panels are SUBSETS OF THE SAME
# SCORES. The scoring universe stays at 142 (human) and 144 (mouse) in every
# panel; only the set of pathways that enter the ordering changes. Nothing here
# recomputes mouse mitoPPS - the artefact is a read-only input with recorded
# provenance - and for symmetry nothing recomputes human mitoPPS either.
#
# INPUT SCALE: the saved mitoPPS universes, built from linear DESeq2-normalised
# expression on both sides. GSVA's log input is used ONLY to define the human
# MYC groups and never touches a mitoPPS number.
#
# EXPLORATORY, POST-HOC, DESCRIPTIVE. Nothing here is pre-registered and nothing
# here is a hypothesis test. The reading rules below are fixed before any number.
#
# SPECIES: the human side is human. No ortholog function is called here or
# anywhere in this repo. The MitoCarta hierarchy is read from the HUMAN
# MitoCarta 3.0 workbook for the human side and from the mouse artefact's own
# pathway_levels for the mouse side; that they agree is a GATE, not an
# assumption.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE25: the pathway ordering, group by group, mouse and human\n",
        strrep("=", 78))

PATH_E25     <- file.path(DIR_RESULTS, "mitopps_group_rankings.rds")
PATH_E25_TAB <- file.path(DIR_TABLES,  "E25_group_orderings.csv")
DIR_E25_FIG  <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_MOUSE   <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")
PATH_MOUSE_README <- here::here("data", "from_myc_mouse", "README.md")
PATH_MITOCARTA_H  <- here::here("data", "mitocarta_human", "Human.MitoCarta3.0.xls")
PATH_E24     <- file.path(DIR_RESULTS, "mitopps_rank_comparison.rds")

set.seed(1)
NBOOT <- 5000L

# =============================================================================
# 0. CONSTANTS AND THE READING RULES, FIXED BEFORE ANY NUMBER
# =============================================================================

MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"

PRIMARY   <- "OXPHOS subunits"
SECONDARY <- "OXPHOS"
MTDNA     <- "mtDNA-encoded OXPHOS subunits"
MTRIB     <- "Mitochondrial ribosome"
READOUTS  <- c(PRIMARY, SECONDARY, MTRIB, MTDNA)

CUT_SIZE  <- 0.5          # |rho(statistic, log size)| at or above this is structural
CUT_AGREE <- 0            # an ordering pair AGREES if rho > 0 with interval clear of 0

# The human MYC groups. Tertiles of the primary activity estimator, middle third
# dropped, so "high" and "low" are separated rather than adjacent. The split is
# made INSIDE each cohort and the cohorts are never pooled.
TERTILE   <- 1/3

# The estimator panel, ordered by proliferation entanglement (CLAUDE.md trap 3).
# The primary is the validation study's M_a. The other four are the sensitivity:
# no MYC-OXPHOS statement in this repo rests on one signature.
EST_PANEL <- c(
  myc_lowent  = "MYC_UP.V1_UP__MITOSTRIP",      # 1.5 pct entangled
  myc_felsher = MYC_REF,                        # 14.8 pct - M_a, THE PRIMARY
  myc_msigdb  = "HALLMARK_MYC_TARGETS_V1__MITOSTRIP",   # 23.5 pct
  myc_regulon = MB_REF,                         # CollecTRI - E21: size-dominated
  myc_mRNA    = "log2MYC")                      # trap 4 - not MYC activity
EST_PRIMARY_NAME <- "myc_felsher"

# Percentile within an ordering. 100 = most prioritised in that group.
.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}

RULES <- list(
  question = paste(
    "compare the pathway ORDERING across the four mouse groups and",
    "MYC-high against MYC-low human tumours, each cohort separately"),
  estimand = paste(
    "per group, the group-mean mitoPPS of each pathway, converted to a",
    "percentile within that group's ordering; and the group CONTRAST, which is",
    "the mouse figures' own statistic"),
  rank_only = paste(
    "no mitoPPS value crosses the species boundary; every cross-species",
    "quantity is an ordering, a rank correlation or a percentile"),
  rerank_only = paste(
    "the hierarchy and mtDNA panels are SUBSETS of the same scores; the",
    "scoring universe is never rebuilt and mouse mitoPPS is never recomputed"),
  batch = paste(
    "the mouse time axis is batch-confounded (batch = timepoint). The author",
    "has ruled the timeline comparisons stay in and that DESeq2 normalisation",
    "handles the batch. Temporal quantities are reported in full and carry the",
    "caveat; the genotype contrast within a timepoint is the cleaner axis"),
  no_pooling = "mouse and human never scored together; TCGA and SCAN-B never averaged",
  exploratory = "post-hoc and descriptive; not a test of any hypothesis",
  n3 = "transcript associations; 'primed' is never written of a transcript",
  visualisation = paste(
    "the figures are for general insight into the orderings and carry no",
    "significance marks; nothing is read off them statistically"))

message("\n0. rules fixed before any number")
message("   estimand  group-mean mitoPPS per pathway -> percentile within the group")
message("   G-A  the mouse contrasts must reproduce the artefact's own pairwise table")
message("   G-B  the shared panel and the MitoCarta hierarchy must agree across species")
message("   G-C  |rho(statistic, log size)| >= ", CUT_SIZE, " -> the statistic is structural")
message("   Q1   does the human MYC contrast order the panel like the mouse Myc effect?")
message("   Q2   which mouse ordering is each human ordering closest to? Report all of it")
message("   TIME timeline comparisons ARE included; the batch confounder is carried")

# =============================================================================
# 1. INPUTS
# =============================================================================

message("\n1. inputs")

if (!file.exists(PATH_MOUSE)) {
  stop("mouse artefact missing: ", PATH_MOUSE, call. = FALSE)
}
md5_now <- unname(tools::md5sum(PATH_MOUSE))
if (!identical(md5_now, MOUSE_MD5)) {
  stop("mouse artefact md5 has changed. README says ", MOUSE_MD5,
       ", file is ", md5_now, ". Re-read data/from_myc_mouse/README.md.",
       call. = FALSE)
}
stopifnot(any(grepl(MOUSE_MD5, readLines(PATH_MOUSE_README, warn = FALSE), fixed = TRUE)))
message("   mouse artefact md5 matches its README: ", md5_now)

mo   <- readRDS(PATH_MOUSE)
mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))

# --- the mouse side ----------------------------------------------------------
ms       <- mo$mitopps_scores
ANN_COLS <- c("sample", "group", "timepoint", "myc_status")
stopifnot(all(ANN_COLS %in% names(ms)))
MPATHS   <- setdiff(names(ms), ANN_COLS)
MM       <- as.matrix(ms[, MPATHS, drop = FALSE])      # samples x pathways
rownames(MM) <- ms$sample
stopifnot(!anyNA(MM), ncol(MM) == mo$n_pathways)

# mitoPPS is centred at exactly 1 per pathway. Assert it rather than trust it -
# the whole argument for this estimand rests on it.
CENTRE_DELTA <- max(abs(colMeans(MM) - 1))
stopifnot(CENTRE_DELTA < 1e-9)
message("   mouse: ", nrow(MM), " samples, ", ncol(MM), " pathways; ",
        "mitoPPS centred at 1 per pathway (max |delta| ",
        format(CENTRE_DELTA, digits = 2), ")")

MOUSE_ARM <- c(`6W_neg` = "mouse 6W WT", `12W_neg` = "mouse 12W WT",
               `6W_pos` = "mouse 6W Myc+", `12W_pos` = "mouse 12W Myc+")
stopifnot(all(ms$group %in% names(MOUSE_ARM)))
mouse_arm <- unname(MOUSE_ARM[ms$group])
message("   mouse groups: ",
        paste(sprintf("%s n=%d", MOUSE_ARM[names(table(ms$group))],
                      as.integer(table(ms$group))), collapse = " | "))

# --- the human side ----------------------------------------------------------
UT <- mito$mitopps_universe
US <- sc$mitopps_universe
stopifnot(is.matrix(UT), is.matrix(US), identical(rownames(UT), rownames(US)))
ID_T <- colnames(UT); ID_S <- colnames(US)
stopifnot(identical(ID_T, colnames(mito$gsva_arms)),
          identical(ID_S, colnames(sc$gsva_arms)))
message("   human mitopps_universe: TCGA ", length(ID_T), " samples | SCAN-B ",
        length(ID_S), " samples | ", nrow(UT), " pathways, identical rownames")

# The human universes are cohort-relative and must also be centred at 1.
stopifnot(max(abs(rowMeans(UT) - 1)) < 1e-9, max(abs(rowMeans(US) - 1)) < 1e-9)
message("   human mitoPPS centred at 1 per pathway in both cohorts")

# --- catalogue sizes ---------------------------------------------------------
HSZ <- vapply(mito$mito_paths, length, integer(1))     # human catalogue sizes
g2p <- mo$gene_to_pathway
MSZ <- c(table(g2p$Pathway))                           # mouse sizes, as detected

# =============================================================================
# 2. GATE A - the mouse contrasts must reproduce the artefact's own table
# =============================================================================
# If this fails, the snapshot is being read wrongly and nothing below it means
# anything. It is the control for the whole mouse side.

message("\n2. GATE A: reproducing the artefact's pairwise contrasts from its scores")

MOUSE_CONTRAST <- list(
  Myc_effect_6W  = c(a = "6W_neg",  b = "6W_pos"),
  Myc_effect_12W = c(a = "12W_neg", b = "12W_pos"),
  `Temporal_Myc-` = c(a = "6W_neg", b = "12W_neg"),
  `Temporal_Myc+` = c(a = "6W_pos", b = "12W_pos"))

.mouse_diff <- function(cn) {
  ia <- ms$group == cn[["a"]]; ib <- ms$group == cn[["b"]]
  colMeans(MM[ib, , drop = FALSE]) - colMeans(MM[ia, , drop = FALSE])
}
mouse_diff <- lapply(MOUSE_CONTRAST, .mouse_diff)

pwt <- mo$mitopps_pairwise
gate_a <- dplyr::bind_rows(lapply(names(MOUSE_CONTRAST), function(ct) {
  saved <- pwt[pwt$contrast == ct, ]
  saved <- saved$diff[match(MPATHS, saved$pathway)]
  tibble::tibble(contrast = ct, n_pathways = length(MPATHS),
                 max_abs_delta = max(abs(mouse_diff[[ct]] - saved)))
}))
print(as.data.frame(gate_a), digits = 3, row.names = FALSE)
GATE_A <- max(gate_a$max_abs_delta) < 1e-10
if (!GATE_A) {
  stop("GATE A FAILED: the recomputed mouse contrasts do not reproduce ",
       "$mitopps_pairwise. The artefact is being read wrongly.", call. = FALSE)
}
message("   GATE A PASSED: all four contrasts reproduce (max |delta| ",
        format(max(gate_a$max_abs_delta), digits = 2), ")")

# =============================================================================
# 3. GATE B - the shared panel, and the hierarchy, must agree across species
# =============================================================================

message("\n3. GATE B: the shared panel and the MitoCarta hierarchy")

shared <- sort(intersect(rownames(UT), MPATHS))
message("   shared: ", length(shared),
        " | human-only: ", length(setdiff(rownames(UT), MPATHS)),
        " | mouse-only: ", paste(sort(setdiff(MPATHS, rownames(UT))), collapse = ", "))

# E24 established this intersection. If its object is on disk, assert we match.
if (file.exists(PATH_E24)) {
  e24 <- readRDS(PATH_E24)
  stopifnot(identical(sort(e24$intersection$shared), shared))
  message("   matches E24's saved intersection exactly")
} else {
  message("   NOTE: E24's object is not on disk; the intersection is recomputed here")
}

# --- the hierarchy, each species from its OWN native source ------------------
# Each species reads its OWN native hierarchy. The full " > " path is compared,
# not just the depth, so this gate is about the shape of the tree and not only
# the number of levels in it.
hx <- readxl::read_excel(PATH_MITOCARTA_H, sheet = 4)
h_hier <- trimws(hx[["MitoPathways Hierarchy"]])
names(h_hier) <- trimws(hx[["MitoPathway"]])

pl <- mo$pathway_levels
LV_COLS <- c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3")
m_hier <- apply(pl[, LV_COLS], 1L, function(z) {
  z <- z[!is.na(z) & nzchar(z)]
  paste(trimws(z), collapse = " > ")
})
names(m_hier) <- trimws(pl$Pathway)

# The synthetic mtDNA pathway is in neither catalogue and has no place in the tree.
cat_shared <- setdiff(shared, MTDNA)
stopifnot(all(cat_shared %in% names(h_hier)), all(cat_shared %in% names(m_hier)))
path_agree <- sum(h_hier[cat_shared] == m_hier[cat_shared])
GATE_B <- path_agree == length(cat_shared)

h_lvl <- lengths(strsplit(h_hier, " > ", fixed = TRUE))
names(h_lvl) <- names(h_hier)
m_lvl <- lengths(strsplit(m_hier, " > ", fixed = TRUE))
names(m_lvl) <- names(m_hier)
hier <- tibble::tibble(
  level = 1:3,
  human = as.integer(table(factor(h_lvl[cat_shared], levels = 1:3))),
  mouse = as.integer(table(factor(m_lvl[cat_shared], levels = 1:3))))
print(as.data.frame(hier), row.names = FALSE)
if (!GATE_B) {
  stop("GATE B FAILED: the hierarchy PATH disagrees on ",
       length(cat_shared) - path_agree, " shared pathways. The level-restricted ",
       "panels are not comparable across species.", call. = FALSE)
}
message("   GATE B PASSED: the full hierarchy PATH agrees on all ", length(cat_shared),
        " shared catalogue pathways, each species read from its own source")
message("   (the synthetic ", MTDNA, " is in neither catalogue and has no place in the tree)")

LEVEL <- h_lvl[cat_shared]      # identical to the mouse's by the gate above

# --- the antichain: the author's "single pathways" ---------------------------
# A pathway is kept if NO OTHER member of the panel sits below it in the tree.
# This is the non-overlapping panel - no member contains another - and it is what
# a size-based measurement should be read on. It is NOT the level-3 slice: the
# primary set sits at level 2 and has no descendants, so it is a leaf of the
# INDUCED tree while being absent from any single level.
.antichain <- function(hh, members) {
  keep <- vapply(members, function(p) {
    !any(startsWith(hh[members], paste0(hh[[p]], " > ")))
  }, TRUE)
  sort(members[keep])
}
anti_h <- .antichain(h_hier, cat_shared)
anti_m <- .antichain(m_hier, cat_shared)
stopifnot(identical(anti_h, anti_m))
ANTICHAIN <- anti_h
message("   antichain (no member contains another): ", length(ANTICHAIN),
        " of ", length(cat_shared),
        " -- identical from both species' own hierarchies")
message("   ", PRIMARY, " is in the antichain: ", PRIMARY %in% ANTICHAIN,
        " | ", SECONDARY, ": ", SECONDARY %in% ANTICHAIN)

# =============================================================================
# 4. PANEL DEFINITIONS - subsets of the SAME scores, never a re-scoring
# =============================================================================

message("\n4. panels (re-ranking only; the scoring universe is unchanged)")

PANELS <- list(
  `shared (142)`         = shared,
  `disjoint (antichain)` = ANTICHAIN,
  `level-3 only`         = sort(names(LEVEL)[LEVEL == 3L]),
  `level-1 only`         = sort(names(LEVEL)[LEVEL == 1L]),
  `no mtDNA pathway`     = sort(setdiff(shared, MTDNA)))
PANEL_MAIN <- "shared (142)"
PANEL_LEAF <- "disjoint (antichain)"

for (nm in names(PANELS)) {
  message(sprintf("   %-18s %3d pathways", nm, length(PANELS[[nm]])))
}
message("   PRIMARY panel: ", PANEL_MAIN,
        " | the non-overlapping panel: ", PANEL_LEAF)
message("   a readout absent from a panel is reported as absent, never as an empty table")

# =============================================================================
# 5. THE HUMAN MYC GROUPS
# =============================================================================

message("\n5. the human MYC groups")

.estimator_vector <- function(coh, w) {
  if (identical(coh, "TCGA")) {
    if (identical(w, "log2MYC")) return(as.numeric(nw$tcga_log2MYC[ID_T]))
    if (w %in% rownames(nw$tcga_M_b_variants)) {
      return(as.numeric(nw$tcga_M_b_variants[w, ID_T]))
    }
    if (w %in% rownames(nw$tcga_gsva_new)) {
      return(as.numeric(nw$tcga_gsva_new[w, ID_T]))
    }
  } else {
    if (identical(w, "log2MYC")) return(as.numeric(sc$log2MYC[ID_S]))
    if (w %in% rownames(sc$M_b_variants)) return(as.numeric(sc$M_b_variants[w, ID_S]))
    if (w %in% rownames(sc$gsva_new))     return(as.numeric(sc$gsva_new[w, ID_S]))
  }
  stop("estimator not scored in this cohort: ", w, " (", coh, ")", call. = FALSE)
}

# Tertile split inside each cohort, middle third dropped.
.myc_groups <- function(v, ids) {
  stopifnot(length(v) == length(ids), !anyNA(v))
  q <- stats::quantile(v, c(TERTILE, 1 - TERTILE), names = FALSE)
  out <- rep(NA_character_, length(v))
  out[v <= q[1]] <- "MYC-low"
  out[v >= q[2]] <- "MYC-high"
  stats::setNames(out, ids)
}

est_T <- vapply(EST_PANEL, function(w) .estimator_vector("TCGA", w),
                numeric(length(ID_T)))
est_S <- vapply(EST_PANEL, function(w) .estimator_vector("SCAN-B", w),
                numeric(length(ID_S)))
rownames(est_T) <- ID_T; rownames(est_S) <- ID_S

grp_T <- .myc_groups(est_T[, EST_PRIMARY_NAME], ID_T)
grp_S <- .myc_groups(est_S[, EST_PRIMARY_NAME], ID_S)

myc_groups <- tibble::tibble(
  cohort = c("TCGA", "TCGA", "SCAN-B", "SCAN-B"),
  group  = c("MYC-low", "MYC-high", "MYC-low", "MYC-high"),
  n      = c(sum(grp_T == "MYC-low", na.rm = TRUE),
             sum(grp_T == "MYC-high", na.rm = TRUE),
             sum(grp_S == "MYC-low", na.rm = TRUE),
             sum(grp_S == "MYC-high", na.rm = TRUE)))
print(as.data.frame(myc_groups), row.names = FALSE)
message("   split on ", EST_PRIMARY_NAME, " (", EST_PANEL[[EST_PRIMARY_NAME]],
        "), tertiles, middle third dropped, inside each cohort")

# --- the TCGA amplification split: DNA-level, the closest analogue to a transgene
tmyc <- readRDS(here::here("data", "from_validation", "tcga_brca_myc_scores.rds"))
amp  <- tmyc$estimators$M_c_amp[match(ID_T, tmyc$estimators$patient)]
names(amp) <- ID_T
grp_T_amp <- ifelse(is.na(amp), NA_character_,
                    ifelse(amp, "MYC-amplified", "not amplified"))
message("   TCGA secondary, DNA-level: ", sum(grp_T_amp == "MYC-amplified", na.rm = TRUE),
        " amplified | ", sum(grp_T_amp == "not amplified", na.rm = TRUE),
        " not | ", sum(is.na(grp_T_amp)), " with no call")

# =============================================================================
# 6. PART A - the eight group orderings
# =============================================================================

message("\n6. PART A: the group orderings")

ARM_ORDER <- c("mouse 6W WT", "mouse 12W WT", "mouse 6W Myc+", "mouse 12W Myc+",
               "TCGA MYC-low", "TCGA MYC-high", "SCAN-B MYC-low", "SCAN-B MYC-high")

# group-mean mitoPPS per pathway, for every arm, on the shared panel
.group_means <- function() {
  out <- list()
  for (g in names(MOUSE_ARM)) {
    out[[unname(MOUSE_ARM[[g]])]] <- colMeans(MM[ms$group == g, shared, drop = FALSE])
  }
  for (lv in c("MYC-low", "MYC-high")) {
    out[[paste("TCGA", lv)]]   <- rowMeans(UT[shared, which(grp_T == lv), drop = FALSE])
    out[[paste("SCAN-B", lv)]] <- rowMeans(US[shared, which(grp_S == lv), drop = FALSE])
  }
  do.call(cbind, out)[shared, ARM_ORDER, drop = FALSE]
}
GMEAN <- .group_means()
stopifnot(!anyNA(GMEAN), identical(rownames(GMEAN), shared))

# the orderings: one percentile vector per arm per panel
orderings <- list()
for (pn in names(PANELS)) {
  p <- PANELS[[pn]]
  orderings[[pn]] <- apply(GMEAN[p, , drop = FALSE], 2L, .pct)
  rownames(orderings[[pn]]) <- p
}
message("   ", length(ARM_ORDER), " orderings x ", length(PANELS), " panels built")

ordering_long <- dplyr::bind_rows(lapply(names(orderings), function(pn) {
  M <- orderings[[pn]]
  tibble::tibble(panel = pn, pathway = rep(rownames(M), times = ncol(M)),
                 arm = rep(colnames(M), each = nrow(M)),
                 pct = as.numeric(M))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER))

# =============================================================================
# 7. PART B - GATE C, the size gate on BOTH estimands
# =============================================================================

message("\n7. PART B: GATE C, the size gate")

.size_of <- function(p, species) {
  s <- if (identical(species, "mouse")) MSZ[p] else HSZ[p]
  stopifnot(!anyNA(s))
  as.integer(s)
}

.rho_ci <- function(x, y) {
  r  <- stats::cor(x, y, method = "spearman")
  bs <- vapply(seq_len(NBOOT), function(i) {
    k <- sample.int(length(x), replace = TRUE)
    if (length(unique(y[k])) < 3L) return(NA_real_)
    stats::cor(x[k], y[k], method = "spearman")
  }, 0)
  c(rho = r,
    lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
    hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)))
}

# (a) the ORDERING estimand: does a group's ordering track set size?
size_ordering <- dplyr::bind_rows(lapply(names(PANELS), function(pn) {
  p <- PANELS[[pn]]
  dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
    sp <- if (grepl("^mouse", a)) "mouse" else "human"
    z  <- .rho_ci(orderings[[pn]][, a], log(.size_of(p, sp)))
    tibble::tibble(panel = pn, arm = a, statistic = "group ordering",
                   rho = z[["rho"]], ci_lo = z[["lo"]], ci_hi = z[["hi"]])
  }))
}))

# (b) the CONTRAST estimand: the mouse figures' own statistic
CONTRASTS <- c("Myc_effect_6W", "Myc_effect_12W", "Temporal_Myc-", "Temporal_Myc+",
               "TCGA MYC-high - low", "SCAN-B MYC-high - low")
CONTRAST_KIND <- c(Myc_effect_6W = "genotype, within timepoint",
                   Myc_effect_12W = "genotype, within timepoint",
                   `Temporal_Myc-` = "TIMELINE - batch-confounded",
                   `Temporal_Myc+` = "TIMELINE - batch-confounded",
                   `TCGA MYC-high - low` = "MYC activity, within cohort",
                   `SCAN-B MYC-high - low` = "MYC activity, within cohort")

.contrast_vec <- function(ct, p) {
  if (ct %in% names(MOUSE_CONTRAST)) return(mouse_diff[[ct]][p])
  if (identical(ct, "TCGA MYC-high - low")) {
    return(GMEAN[p, "TCGA MYC-high"] - GMEAN[p, "TCGA MYC-low"])
  }
  GMEAN[p, "SCAN-B MYC-high"] - GMEAN[p, "SCAN-B MYC-low"]
}

size_contrast <- dplyr::bind_rows(lapply(names(PANELS), function(pn) {
  p <- PANELS[[pn]]
  dplyr::bind_rows(lapply(CONTRASTS, function(ct) {
    sp <- if (ct %in% names(MOUSE_CONTRAST)) "mouse" else "human"
    z  <- .rho_ci(.contrast_vec(ct, p), log(.size_of(p, sp)))
    tibble::tibble(panel = pn, arm = ct, statistic = "group contrast",
                   rho = z[["rho"]], ci_lo = z[["lo"]], ci_hi = z[["hi"]])
  }))
}))

size_gate <- dplyr::bind_rows(size_ordering, size_contrast) %>%
  dplyr::mutate(structural = abs(rho) >= CUT_SIZE)

message("   on the PRIMARY panel (", PANEL_MAIN, "):")
size_gate %>%
  dplyr::filter(panel == PANEL_MAIN) %>%
  dplyr::select(statistic, arm, rho, ci_lo, ci_hi, structural) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

G_C_ORDER <- any(size_gate$structural[size_gate$statistic == "group ordering" &
                                        size_gate$panel == PANEL_MAIN])
G_C_CONTR <- any(size_gate$structural[size_gate$statistic == "group contrast" &
                                        size_gate$panel == PANEL_MAIN])
message("\n   GATE C, group ordering:  ",
        if (G_C_ORDER) "STRUCTURAL in at least one arm" else "clear in every arm")
message("   GATE C, group contrast:  ",
        if (G_C_CONTR) "STRUCTURAL in at least one arm" else "clear in every arm")
PANEL_FOR_READ <- if (G_C_CONTR) PANEL_LEAF else PANEL_MAIN
message("   the contrast comparison is read on: ", PANEL_FOR_READ)

# =============================================================================
# 8. PART C - the ordering agreement matrix, all eight against all eight
# =============================================================================

message("\n8. PART C: agreement between the orderings")

.pairwise_agree <- function(M, lab) {
  a <- colnames(M)
  gr <- utils::combn(length(a), 2L)
  dplyr::bind_rows(lapply(seq_len(ncol(gr)), function(i) {
    x <- M[, gr[1, i]]; y <- M[, gr[2, i]]
    z <- .rho_ci(x, y)
    tibble::tibble(panel = lab, arm_a = a[gr[1, i]], arm_b = a[gr[2, i]],
                   rho = z[["rho"]], ci_lo = z[["lo"]], ci_hi = z[["hi"]])
  }))
}

agree <- dplyr::bind_rows(lapply(names(PANELS), function(pn) {
  .pairwise_agree(orderings[[pn]], pn)
})) %>%
  dplyr::mutate(
    species_a = dplyr::if_else(grepl("^mouse", arm_a), "mouse", "human"),
    species_b = dplyr::if_else(grepl("^mouse", arm_b), "mouse", "human"),
    pair = dplyr::case_when(
      species_a != species_b ~ "cross-species",
      species_a == "mouse"   ~ "within mouse",
      TRUE                   ~ "within human"),
    agrees = rho > CUT_AGREE & ci_lo > CUT_AGREE)

message("   cross-species pairs on the primary panel:")
agree %>%
  dplyr::filter(panel == PANEL_MAIN, pair == "cross-species") %>%
  dplyr::select(arm_a, arm_b, rho, ci_lo, ci_hi, agrees) %>%
  as.data.frame() %>%
  print(digits = 2, row.names = FALSE)

message("\n   within-species pairs on the primary panel:")
agree %>%
  dplyr::filter(panel == PANEL_MAIN, pair != "cross-species") %>%
  dplyr::select(pair, arm_a, arm_b, rho, ci_lo, ci_hi) %>%
  as.data.frame() %>%
  print(digits = 2, row.names = FALSE)

# For each human arm, which mouse arm does its ordering sit closest to? The rule
# forbids assigning to the nearest when the intervals overlap.
closest <- dplyr::bind_rows(lapply(c("TCGA MYC-low", "TCGA MYC-high",
                                     "SCAN-B MYC-low", "SCAN-B MYC-high"),
  function(h) {
    d <- agree %>%
      dplyr::filter(panel == PANEL_MAIN, pair == "cross-species",
                    arm_a == h | arm_b == h) %>%
      dplyr::mutate(mouse_arm = dplyr::if_else(arm_a == h, arm_b, arm_a)) %>%
      dplyr::arrange(dplyr::desc(rho))
    best <- d[1, ]; second <- d[2, ]
    tibble::tibble(human_arm = h, nearest = best$mouse_arm, rho = best$rho,
                   ci_lo = best$ci_lo, ci_hi = best$ci_hi,
                   runner_up = second$mouse_arm, rho_2 = second$rho,
                   separated = best$ci_lo > second$ci_hi)
  }))
message("\n   nearest mouse ordering to each human ordering:")
print(as.data.frame(closest), digits = 2, row.names = FALSE)

# =============================================================================
# 9. PART D - where the readout pathways sit in each ordering
# =============================================================================

message("\n9. PART D: the readout pathways in each ordering")

# The ordering is a single number per group, so its uncertainty is the group's
# sampling uncertainty. Resample SAMPLES within the arm and re-rank.
.arm_samples <- function(a) {
  if (grepl("^mouse", a)) {
    g <- names(MOUSE_ARM)[MOUSE_ARM == a]
    return(list(kind = "mouse", idx = which(ms$group == g)))
  }
  if (grepl("^TCGA", a)) {
    return(list(kind = "TCGA", idx = which(grp_T == sub("^TCGA ", "", a))))
  }
  list(kind = "SCAN-B", idx = which(grp_S == sub("^SCAN-B ", "", a)))
}

.boot_pct <- function(a, p, path) {
  z <- .arm_samples(a)
  M <- switch(z$kind,
              mouse    = t(MM[z$idx, p, drop = FALSE]),
              TCGA     = UT[p, z$idx, drop = FALSE],
              `SCAN-B` = US[p, z$idx, drop = FALSE])
  n <- ncol(M)
  bs <- vapply(seq_len(NBOOT), function(i) {
    k <- sample.int(n, replace = TRUE)
    .pct(rowMeans(M[, k, drop = FALSE]))[match(path, p)]
  }, 0)
  c(lo = unname(stats::quantile(bs, 0.025)), hi = unname(stats::quantile(bs, 0.975)))
}

positions <- dplyr::bind_rows(lapply(names(PANELS), function(pn) {
  p <- PANELS[[pn]]
  rr <- intersect(READOUTS, p)
  dplyr::bind_rows(lapply(rr, function(path) {
    dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
      ci <- .boot_pct(a, p, path)
      tibble::tibble(panel = pn, pathway = path, arm = a,
                     pct = orderings[[pn]][path, a],
                     ci_lo = ci[["lo"]], ci_hi = ci[["hi"]])
    }))
  }))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER))

message("   ", PRIMARY, " on the primary panel:")
positions %>%
  dplyr::filter(panel == PANEL_MAIN, pathway == PRIMARY) %>%
  dplyr::select(arm, pct, ci_lo, ci_hi) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

message("\n   ", PRIMARY, " on ", PANEL_LEAF, " (no member contains another):")
if (PRIMARY %in% PANELS[[PANEL_LEAF]]) {
  positions %>%
    dplyr::filter(panel == PANEL_LEAF, pathway == PRIMARY) %>%
    dplyr::select(arm, pct, ci_lo, ci_hi) %>%
    as.data.frame() %>%
    print(digits = 3, row.names = FALSE)
} else {
  message("   ABSENT from this panel - it contains another member, so it is not",
          " an antichain element. Nothing is printed rather than an empty table.")
}

# which readouts exist in which panel, so no absence is ever silent
readout_membership <- tibble::tibble(
  pathway = rep(READOUTS, times = length(PANELS)),
  panel   = rep(names(PANELS), each = length(READOUTS)),
  present = as.logical(unlist(lapply(PANELS, function(p) READOUTS %in% p))))
message("\n   readout membership by panel:")
readout_membership %>%
  tidyr::pivot_wider(names_from = panel, values_from = present) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

message("\n   the other readouts on the primary panel:")
positions %>%
  dplyr::filter(panel == PANEL_MAIN, pathway != PRIMARY) %>%
  dplyr::select(pathway, arm, pct) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

# --- the estimator sensitivity: does the human position move with the ruler? --
est_sens <- dplyr::bind_rows(lapply(names(EST_PANEL), function(en) {
  gT <- .myc_groups(est_T[, en], ID_T)
  gS <- .myc_groups(est_S[, en], ID_S)
  p  <- PANELS[[PANEL_MAIN]]
  dplyr::bind_rows(
    tibble::tibble(estimator = en, cohort = "TCGA", group = c("MYC-low", "MYC-high"),
                   pct = c(.pct(rowMeans(UT[p, which(gT == "MYC-low"), drop = FALSE]))[match(PRIMARY, p)],
                           .pct(rowMeans(UT[p, which(gT == "MYC-high"), drop = FALSE]))[match(PRIMARY, p)])),
    tibble::tibble(estimator = en, cohort = "SCAN-B", group = c("MYC-low", "MYC-high"),
                   pct = c(.pct(rowMeans(US[p, which(gS == "MYC-low"), drop = FALSE]))[match(PRIMARY, p)],
                           .pct(rowMeans(US[p, which(gS == "MYC-high"), drop = FALSE]))[match(PRIMARY, p)])))
})) %>%
  tidyr::pivot_wider(names_from = group, values_from = pct) %>%
  dplyr::mutate(shift = `MYC-high` - `MYC-low`)
message("\n   estimator sensitivity, ", PRIMARY, " percentile (trap 3):")
print(as.data.frame(est_sens), digits = 3, row.names = FALSE)

# --- the TCGA amplification split, DNA-level ---------------------------------
p_main <- PANELS[[PANEL_MAIN]]
amp_pos <- tibble::tibble(
  group = c("not amplified", "MYC-amplified"),
  n     = c(sum(grp_T_amp == "not amplified", na.rm = TRUE),
            sum(grp_T_amp == "MYC-amplified", na.rm = TRUE)),
  pct   = c(.pct(rowMeans(UT[p_main, which(grp_T_amp == "not amplified"), drop = FALSE]))[match(PRIMARY, p_main)],
            .pct(rowMeans(UT[p_main, which(grp_T_amp == "MYC-amplified"), drop = FALSE]))[match(PRIMARY, p_main)]))
message("\n   TCGA DNA-level split, ", PRIMARY, " percentile:")
print(as.data.frame(amp_pos), digits = 3, row.names = FALSE)

# =============================================================================
# 10. PART E - the contrasts, and whether they order the panel alike
# =============================================================================

message("\n10. PART E: the contrasts (the mouse figures' own statistic)")

contrast_agree <- dplyr::bind_rows(lapply(names(PANELS), function(pn) {
  p <- PANELS[[pn]]
  hs <- c("TCGA MYC-high - low", "SCAN-B MYC-high - low")
  dplyr::bind_rows(lapply(hs, function(h) {
    dplyr::bind_rows(lapply(names(MOUSE_CONTRAST), function(mc) {
      z <- .rho_ci(.contrast_vec(mc, p), .contrast_vec(h, p))
      tibble::tibble(panel = pn, human = h, mouse = mc,
                     kind = unname(CONTRAST_KIND[[mc]]),
                     rho = z[["rho"]], ci_lo = z[["lo"]], ci_hi = z[["hi"]],
                     agrees = z[["rho"]] > CUT_AGREE & z[["lo"]] > CUT_AGREE)
    }))
  }))
}))

message("   cross-species agreement of the CONTRAST ordering, on ", PANEL_FOR_READ, ":")
contrast_agree %>%
  dplyr::filter(panel == PANEL_FOR_READ) %>%
  dplyr::select(human, mouse, kind, rho, ci_lo, ci_hi, agrees) %>%
  as.data.frame() %>%
  print(digits = 2, row.names = FALSE)

# where the readouts sit in each contrast's own ordering
contrast_pos <- dplyr::bind_rows(lapply(names(PANELS), function(pn) {
  p <- PANELS[[pn]]
  dplyr::bind_rows(lapply(CONTRASTS, function(ct) {
    v <- .pct(.contrast_vec(ct, p))
    rr <- intersect(READOUTS, p)
    tibble::tibble(panel = pn, contrast = ct, kind = unname(CONTRAST_KIND[[ct]]),
                   pathway = rr, pct = v[match(rr, p)])
  }))
}))
message("\n   ", PRIMARY, " inside each contrast's ordering, primary panel:")
contrast_pos %>%
  dplyr::filter(panel == PANEL_MAIN, pathway == PRIMARY) %>%
  dplyr::select(contrast, kind, pct) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

# =============================================================================
# 11. PART F - figures. INSIGHT ONLY. No significance marks, nothing read off.
# =============================================================================

message("\n11. PART F: figures (insight only)")

if (!dir.exists(DIR_E25_FIG)) dir.create(DIR_E25_FIG, recursive = TRUE)

CAP_BATCH <- paste(
  "Timeline comparisons are included by decision; the mouse time axis is",
  "batch-confounded (batch = timepoint) and DESeq2 normalisation is what",
  "handles it. Genotype contrasts within a timepoint are the cleaner axis.")
CAP_RANK <- paste(
  "RANK ONLY: no mitoPPS value is compared across species. Exploratory,",
  "post-hoc, descriptive; for insight, not for statistical evaluation.")

# --- figure 1: the readout pathways across all eight orderings ---------------
f1_dat <- positions %>%
  dplyr::filter(panel %in% c(PANEL_MAIN, PANEL_LEAF)) %>%
  dplyr::mutate(
    species = dplyr::if_else(grepl("^mouse", arm), "mouse", "human"),
    pathway = factor(pathway, levels = READOUTS))

p1 <- ggplot2::ggplot(f1_dat, ggplot2::aes(arm, pct, colour = species)) +
  ggplot2::geom_hline(yintercept = 50, linewidth = 0.3, colour = "grey70",
                      linetype = 2) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                          linewidth = 0.4, alpha = 0.7) +
  ggplot2::geom_point(size = 1.9) +
  ggplot2::facet_grid(pathway ~ panel) +
  ggplot2::scale_colour_manual(values = c(mouse = "#0072B2", human = "#D55E00")) +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::labs(
    x = NULL, y = "percentile within the group's ordering\n(100 = most prioritised)",
    title = "Where each pathway sits in each group's ordering",
    subtitle = paste("Eight orderings of the same shared panel. Interval =",
                     "bootstrap over samples within the group."),
    caption = paste(CAP_RANK, CAP_BATCH, sep = "\n")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))

ggplot2::ggsave(file.path(DIR_E25_FIG, "E25_group_orderings.pdf"), p1,
                width = 9, height = 9)

# --- figure 2: fig02 layout, the contrast ordering, mouse against human ------
f2_dat <- dplyr::bind_rows(lapply(names(MOUSE_CONTRAST), function(mc) {
  p <- PANELS[[PANEL_FOR_READ]]
  mv <- .pct(.contrast_vec(mc, p))
  dplyr::bind_rows(lapply(c("TCGA MYC-high - low", "SCAN-B MYC-high - low"),
    function(h) {
      tibble::tibble(mouse_contrast = mc, kind = unname(CONTRAST_KIND[[mc]]),
                     human = h, pathway = p, mouse_pct = mv,
                     human_pct = .pct(.contrast_vec(h, p)))
    }))
})) %>%
  dplyr::mutate(
    readout = dplyr::if_else(pathway %in% READOUTS, pathway, NA_character_),
    mouse_contrast = factor(mouse_contrast, levels = names(MOUSE_CONTRAST)))

p2 <- ggplot2::ggplot(f2_dat, ggplot2::aes(mouse_pct, human_pct)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3,
                       colour = "grey70", linetype = 2) +
  ggplot2::geom_point(colour = "grey55", size = 0.8, alpha = 0.7) +
  ggplot2::geom_point(data = subset(f2_dat, !is.na(readout)),
                      ggplot2::aes(fill = readout), shape = 21, size = 2.6,
                      colour = "black", stroke = 0.4) +
  ggplot2::facet_grid(human ~ mouse_contrast) +
  ggplot2::labs(
    x = "pathway rank in the MOUSE contrast (percentile)",
    y = "pathway rank in the HUMAN contrast (percentile)",
    fill = NULL,
    title = "Do MYC-high human tumours reorder the compartment like the mouse?",
    subtitle = paste0("Each point is a pathway on the ", PANEL_FOR_READ,
                      " panel. Ranks, not values, on both axes."),
    caption = paste(CAP_RANK, CAP_BATCH, sep = "\n")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))

ggplot2::ggsave(file.path(DIR_E25_FIG, "E25_contrast_rank_vs_rank.pdf"), p2,
                width = 10, height = 6)

message("   wrote E25_group_orderings.pdf and E25_contrast_rank_vs_rank.pdf")

# =============================================================================
# 12. PART G - the verdict
# =============================================================================

message("\n12. PART G: verdict")

.q1 <- function(h) {
  d <- contrast_agree %>%
    dplyr::filter(panel == PANEL_FOR_READ, human == h,
                  mouse %in% c("Myc_effect_6W", "Myc_effect_12W"))
  if (all(d$agrees)) return("AGREES with the mouse Myc effect at both timepoints")
  if (any(d$agrees)) {
    return(paste0("agrees at ", paste(sub("Myc_effect_", "", d$mouse[d$agrees]),
                                      collapse = " and "), " only"))
  }
  if (all(d$ci_hi < 0)) return("OPPOSITE to the mouse Myc effect")
  "NO agreement established - the interval covers zero"
}
Q1 <- vapply(c("TCGA MYC-high - low", "SCAN-B MYC-high - low"), .q1, "")

R2_AGREE <- identical(unname(Q1[[1]]), unname(Q1[[2]]))

verdict <- tibble::tibble(
  rule = c("GATE A mouse reproduction", "GATE B panel and hierarchy",
           "GATE C size, ordering", "GATE C size, contrast",
           "Q1 TCGA", "Q1 SCAN-B", "R2 cohort agreement", "TIME batch confounder"),
  outcome = c(
    paste0("PASSED - all four contrasts reproduce (max |delta| ",
           format(max(gate_a$max_abs_delta), digits = 2), ")"),
    paste0("PASSED - ", length(shared), " shared pathways; the full hierarchy PATH ",
           "agrees on all ", length(cat_shared), " catalogue members, and the ",
           length(ANTICHAIN), "-member antichain is identical from both species' ",
           "own sources"),
    if (G_C_ORDER) {
      "STRUCTURAL in at least one arm - read the level-restricted panels beside it"
    } else {
      "clear in every arm on the primary panel"
    },
    if (G_C_CONTR) {
      paste0("STRUCTURAL - the contrast comparison is read on ", PANEL_LEAF)
    } else {
      "clear in every arm - the contrast is size-free, unlike E24's estimand"
    },
    unname(Q1[[1]]), unname(Q1[[2]]),
    if (R2_AGREE) "the cohorts agree" else "THE COHORTS DISAGREE - that is the result",
    paste("carried, not used to withdraw anything. The two temporal contrasts",
          "are batch-confounded by design (batch = timepoint) and are reported",
          "in full by the author's decision")))
print(as.data.frame(verdict), row.names = FALSE, right = FALSE)

message("\n   AND WHATEVER IT SAYS, IT IS DESCRIPTIVE. NOT A TEST.")

NOTES <- c(
  paste("Estimand: group-mean mitoPPS per pathway, converted to a percentile",
        "within that group's ordering. This is the mouse figures' family of",
        "statistic, not E24's within-sample percentile."),
  paste("E24's estimand was size-biased asymmetrically - rho(mean within-sample",
        "percentile, log size) +0.165 mouse against +0.536 TCGA and +0.571",
        "SCAN-B - so it compared two differently biased statistics. The",
        "mechanism is the variance-to-rank conversion, not the score: mitoPPS is",
        "centred at exactly 1 per pathway in both species, asserted in section 1."),
  paste("Re-ranking, never re-scoring. The hierarchy and mtDNA panels are",
        "subsets of the same scores; the scoring universe stays 142 human and",
        "144 mouse. Mouse mitoPPS is never recomputed in this repo."),
  paste("The non-overlapping panel is the ANTICHAIN of", length(ANTICHAIN),
        "pathways, no member of which contains another. It is not a single",
        "level: the primary set sits at level 2 with no descendants in the",
        "panel, so it is an antichain member while being absent from the",
        "level-3 slice. The antichain was built independently from each",
        "species' own hierarchy and the two are identical."),
  paste("The MitoCarta hierarchy was read from the HUMAN workbook for the human",
        "side and the artefact's own pathway_levels for the mouse side. That",
        "they agree on all", length(cat_shared), "shared catalogue pathways is",
        "GATE B, not an assumption. No ortholog function is called."),
  paste("TIMELINE: included by the author's decision. The mouse time axis is",
        "batch-confounded (batch = timepoint); DESeq2 normalisation is what",
        "handles it, the confounder is stated wherever a temporal number",
        "appears, and no temporal number is withdrawn on its account."),
  paste("The human MYC split is tertiles of", EST_PANEL[[EST_PRIMARY_NAME]],
        "inside each cohort with the middle third dropped. Four other",
        "estimators are carried as the sensitivity because no MYC statement in",
        "this repo rests on one signature (trap 3). TCGA also carries the",
        "DNA-level amplification split, which SCAN-B cannot."),
  paste("PURITY IS NOT ADDRESSED HERE and the human numbers are unadjusted.",
        "That is E26's subject, and trap 15 governs it: a positive control is",
        "read FIRST, before any adjusted number."),
  "RANK ONLY. No mitoPPS value is saved in this object and none appears in the figures.",
  "Exploratory and post-hoc. The figures are for insight, not statistical evaluation.")

saveRDS(list(
  rules = RULES, settings = list(
    question = RULES$question, primary = PRIMARY, secondary = SECONDARY,
    readouts = READOUTS, cut_size = CUT_SIZE, tertile = TERTILE,
    nboot = NBOOT, seed = 1, estimator_panel = EST_PANEL,
    estimator_primary = EST_PRIMARY_NAME, arm_order = ARM_ORDER,
    mouse_artefact = list(path = "data/from_myc_mouse/mitopps_scores.rds",
                          md5 = MOUSE_MD5, recomputed = FALSE)),
  gate_a = gate_a, hierarchy = hier, panels = PANELS,
  antichain = ANTICHAIN, readout_membership = readout_membership,
  intersection = list(n_shared = length(shared), shared = shared,
                      human_only = setdiff(rownames(UT), MPATHS),
                      mouse_only = setdiff(MPATHS, rownames(UT))),
  myc_groups = myc_groups, amp_split = amp_pos,
  orderings = ordering_long, size_gate = size_gate,
  agreement = agree, closest = closest, positions = positions,
  estimator_sensitivity = est_sens,
  contrast_agreement = contrast_agree, contrast_positions = contrast_pos,
  panel_for_read = PANEL_FOR_READ, q1 = Q1, r2_agree = R2_AGREE,
  verdict = verdict, analysis_date = Sys.Date(), notes = NOTES), PATH_E25)

readr::write_csv(ordering_long %>% dplyr::mutate(arm = as.character(arm)),
                 PATH_E25_TAB)

message("\nE25: done.")
message("    ", PATH_E25)
message("    ", PATH_E25_TAB)
message("    ", file.path(DIR_E25_FIG, "E25_group_orderings.pdf"))
message("    ", file.path(DIR_E25_FIG, "E25_contrast_rank_vs_rank.pdf"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E25)

  ## the gates first, in the order they must be read
  x$gate_a |> as.data.frame()
  x$hierarchy |> as.data.frame()
  x$verdict |> as.data.frame()

  ## GATE C: is either estimand size-structural?
  subset(x$size_gate, panel == "shared (142)") |> as.data.frame()

  ## the size comparison that motivated the rebuild, E24 against E25
  ## (E24's within-sample percentile was +0.54 / +0.57 in human, +0.165 in mouse)
  subset(x$size_gate, statistic == "group contrast" & panel == "shared (142)") |>
    as.data.frame()

  ## the primary readout: where OXPHOS subunits sits in all eight orderings
  subset(x$positions, panel == "shared (142)" & pathway == "OXPHOS subunits") |>
    as.data.frame()
  subset(x$positions, panel == "level-3 leaves" & pathway == "OXPHOS subunits") |>
    as.data.frame()

  ## which mouse ordering is each human ordering nearest, and is it separated?
  x$closest |> as.data.frame()

  ## the cross-species question, on the contrast
  subset(x$contrast_agreement, panel == x$panel_for_read) |> as.data.frame()

  ## the timeline, carried with its caveat
  subset(x$contrast_positions,
         panel == "shared (142)" & pathway == "OXPHOS subunits") |> as.data.frame()

  ## trap 3: does the human answer move with the estimator?
  x$estimator_sensitivity |> as.data.frame()

  ## the DNA-level split TCGA can do and SCAN-B cannot
  x$amp_split |> as.data.frame()

  ## the notes, which carry every caveat the note must repeat
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")
}
