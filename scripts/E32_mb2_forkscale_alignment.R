# E32_mb2_forkscale_alignment.R
# =============================================================================
# WHICH BICLUSTER PROGRAMME DOES THIS ARM'S STATE CORRESPOND TO?
#
# EXPLORATORY. Nothing here is pre-registered.
#
# The direction, the statistic and the verdict rule were fixed in
# docs/2026-09-28_mb2_forkscale_declaration.md and COMMITTED ALONE AT a71e19d,
# BEFORE myc_human_validation/results/forkscale_models.rds was opened. The
# inputs and what opening that object revealed are recorded at fb98f33. Read
# both before this file; this header restates their rules, it does not invent
# them.
#
# =============================================================================
# THE QUESTION
# =============================================================================
# Menegollo describe two biclusters that both carry increased mitochondrial
# biogenesis from DIFFERENT programmes:
#
#   MB1_UF  mature luminal      biogenesis WITHOUT Myc signals   almost all ER+
#   MB2_UF  luminal progenitor  biogenesis WITH Myc signals      ER-NEGATIVE only
#
# Every survival panel in that paper is MB1. MB2 has none.
#
# This arm's state is MYC-high and OXPHOS-high. Which programme is it?
#
# THE DIRECTION COMES FROM THEIR PUBLISHED TEXT, NOT FROM US, so it could not
# have been tuned to our data: MYC activity should track MB2 and not MB1.
#
# =============================================================================
# THE STATISTIC, AND WHY IT IS A PARTIAL
# =============================================================================
# Supplementary Fig. S7A puts the two forkscales on a clear positive diagonal -
# measured here at rho = 0.712 - so "different programmes" does not mean
# "orthogonal axes", and MARGINAL CORRELATIONS CANNOT SEPARATE THEM. The
# discriminating quantities are the two partials:
#
#   rho( M_a , MB2_forkscale | MB1_forkscale )
#   rho( M_a , MB1_forkscale | MB2_forkscale )
#
# Partial SPEARMAN, because forkscale is severely skewed by construction
# (pc1 / index, with index near 1 giving huge values; MB1 observed near
# [-42.6, 8.8]). Pearson is reported as a companion only.
#
# Bootstrap percentile intervals, 10,000 resamples, and THE WHOLE PIPELINE IS
# BOOTSTRAPPED INCLUDING THE RANKING - not the coefficient alone.
#
# =============================================================================
# WHAT IS DECLARED UNINFORMATIVE IN ADVANCE
# =============================================================================
# OXPHOS against either forkscale. BOTH biclusters carry mitochondrial
# biogenesis, so the comparison discriminates NOTHING. It is computed, it is
# reported, and IT IS NOT PRESENTED AS SUPPORT IN EITHER DIRECTION whatever it
# returns. Declaration section 3, P2.
#
# =============================================================================
# THE PROLIFERATION COMPANION - POST-HOC, AND LABELLED AS SUCH
# =============================================================================
# NOT in the declaration. Chosen AFTER reading forkscale_models.rds, for a
# reason recorded at fb98f33 before this file was written: MB2_forkscale
# correlates with PROLIF_DISJOINT at +0.874, HIGHER than its +0.739 with M_a,
# and M_a is itself 14.8% proliferation-entangled (CLAUDE.md trap 3). So the
# likeliest alternative explanation of a positive result is that MB2 forkscale
# is close to a proliferation axis - which is roughly Menegollo's description
# of MB1.
#
# THE COMPANION CANNOT CHANGE THE VERDICT. The verdict is computed on the
# declared quantities and nothing else. The companion is reported beside it and
# bears on what a verdict LICENSES, which is a different thing.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - NO outcome or survival variable enters at any point. Not joined, not read.
#   - Nothing is written to myc_human_validation. It is read-only and frozen.
#   - No per-gene FDR, consistent with the rest of the project.
#   - MB3 is not used.
#   - Nothing is re-scored: every TCGA quantity is read from an existing object.
#
# N3: these are two exposure axes and a transcript score. Nothing is "primed".
# SPECIES: human throughout. NO ORTHOLOG FUNCTION IS CALLED.
# COHORT-RELATIVITY: GSVA and every z-score are cohort-relative. Everything here
# is ONE cohort (TCGA), so no pooling arises - and no value computed here is
# ever compared numerically to a SCAN-B, CCLE or mouse value.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE32: MB2 forkscale alignment\n", strrep("=", 78))

PATH_E32     <- file.path(DIR_RESULTS, "mb2_forkscale_alignment.rds")
PATH_E32_COR <- file.path(DIR_TABLES,  "E32_forkscale_correlations.csv")
PATH_E32_FIG <- file.path(DIR_FIGURES, "E32_mb1_vs_mb2_forkscale.png")
# Written twice: outputs/ is gitignored, and on this branch the note and the
# figure are the only durable record. Same reasoning as E31.
PATH_E32_FIG_DOC <- here::here("docs", "figures",
                               "2026-09-28_E32_mb1_vs_mb2_forkscale.png")

# =============================================================================
# 0. Constants and THE VERDICT RULE, transcribed from the declaration
# =============================================================================
DIR_MENEGOLLO <- here::here("data", "menegollo_biclusters")
PATH_MB1 <- file.path(DIR_MENEGOLLO, "TCGA_MB1_RNASeq_data_nonorm.RData")
PATH_MB2 <- file.path(DIR_MENEGOLLO, "TCGA_MB2_RNASeq_data_nonorm.RData")
PATH_ALL <- file.path(DIR_MENEGOLLO, "TCGA.all.biclusters.RNAseq.Rdata")
# Read-only, from the FROZEN repo. Nothing is written there, ever.
PATH_F3PRE <- "/Users/gs/code/myc_human_validation/results/forkscale_models.rds"
PATH_STATE <- "/Users/gs/code/myc_human_validation/results/state_definition.rds"

N_BOOT     <- 10000L
CI_LEVEL   <- 0.95
MIN_OVERLAP <- 800L     # the declared stop
EXPECT_OVERLAP_AT_LEAST <- 880L
EXPECT_FORK_N <- 1037L

MYC_ESTIMATORS <- c(M_a = "M_a", M_b = "M_b", M_c = "M_c")
MYC_PRIMARY    <- "M_a"
OX_RULERS      <- c("ox_gsva", "ox_ppd")
FORKS          <- c("MB1_forkscale", "MB2_forkscale")

# --- THE VERDICT RULE. Declaration section 5, transcribed, not reinterpreted.
# MB2-ALIGNED  rho(M_a, MB2|MB1) positive, CI excludes zero, AND larger than
#              rho(M_a, MB1|MB2)
# MB1-ALIGNED  the reverse
# UNRESOLVED   anything else, INCLUDING BOTH POSITIVE WITH OVERLAPPING CIs
#
# The overlap clause is listed under UNRESOLVED, so non-overlap of the two
# intervals is a REQUIREMENT of an aligned verdict and is implemented as one.
# UNRESOLVED IS AN ACCEPTABLE OUTCOME: the manuscript sentence simply does not
# get written, and limb (i)'s Menegollo clause is dropped rather than softened.
VERDICT_RULE <- paste0(
  "MB2-ALIGNED iff rho(M_a, MB2|MB1) > 0, its CI excludes zero, its point ",
  "estimate exceeds rho(M_a, MB1|MB2), and the two CIs DO NOT OVERLAP. ",
  "MB1-ALIGNED is the mirror. Everything else, including both positive with ",
  "overlapping intervals, is UNRESOLVED.")

message("\n0. the rule, fixed at a71e19d before any number:\n   ", VERDICT_RULE)

# =============================================================================
# 1. Forkscale, recomputed from the per-bicluster files
# =============================================================================
# THE PER-BICLUSTER FILES ARE PRIMARY. TCGA.all.biclusters.RNAseq.Rdata is 849
# rows, not 1,037, because an upstream inner_join to two files absent from the
# snapshot drops 188 samples. It is loaded ONLY as a cross-check.
message("\n1. forkscale, from the per-bicluster files")

.load_one <- function(path, want) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!want %in% ls(e)) {
    stop(basename(path), " does not contain ", want, ". It holds: ",
         paste(ls(e), collapse = ", "), call. = FALSE)
  }
  get(want, envir = e)
}

MB1 <- .load_one(PATH_MB1, "TCGA.MB1.RNAseq.df")
MB2 <- .load_one(PATH_MB2, "TCGA.MB2.RNAseq.df")
for (nm in c("MB1", "MB2")) {
  d <- get(nm)
  k <- paste0(nm, c(".pc1", ".index"))
  if (!all(c("X", k) %in% names(d))) {
    stop(nm, " lacks X / pc1 / index. Columns: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  }
  if (nrow(d) != EXPECT_FORK_N) {
    stop(nm, " has ", nrow(d), " rows, not ", EXPECT_FORK_N,
         ". The snapshot has changed underneath.", call. = FALSE)
  }
}
# index is the MCbiclust SampleSort rank, a permutation of 1:1037. If it is not,
# forkscale is not what the upstream scripts computed.
for (nm in c("MB1", "MB2")) {
  ix <- get(nm)[[paste0(nm, ".index")]]
  if (!setequal(ix, seq_len(EXPECT_FORK_N))) {
    stop(nm, ".index is not a permutation of 1:", EXPECT_FORK_N, call. = FALSE)
  }
}
message("   MB1 and MB2: ", EXPECT_FORK_N, " rows each, index a clean permutation")

# Patient key is the 12-character barcode; X is the 28-character aliquot.
FORK <- tibble::tibble(
  patient   = substr(MB1$X, 1, 12),
  aliquot   = MB1$X,
  MB1_pc1   = MB1$MB1.pc1,
  MB1_index = MB1$MB1.index)
m2 <- match(FORK$aliquot, MB2$X)
if (anyNA(m2)) {
  stop(sum(is.na(m2)), " MB1 aliquot(s) absent from the MB2 file. The two ",
       "per-bicluster frames are supposed to be the same 1,037 samples.",
       call. = FALSE)
}
FORK$MB2_pc1   <- MB2$MB2.pc1[m2]
FORK$MB2_index <- MB2$MB2.index[m2]
if (anyDuplicated(FORK$patient)) {
  stop(sum(duplicated(FORK$patient)), " duplicated patient barcode(s) in the ",
       "forkscale frame. Upstream says 1,037 rows, 1,037 unique patients.",
       call. = FALSE)
}

# THE DEFINITION, from the upstream Fig_7ABC script. Plain form carries the
# verdict; the log form is reported alongside only.
FORK$MB1_forkscale <- FORK$MB1_pc1 / FORK$MB1_index
FORK$MB2_forkscale <- FORK$MB2_pc1 / FORK$MB2_index
# forkscale.log is Inf where index == 1, because log(1) == 0. Inf is NOT NA and
# will survive complete.cases() and destroy a fit. Set it to NA explicitly and
# say how many, rather than letting it through.
.fs_log <- function(pc1, ix) {
  v <- pc1 / log(ix)
  v[!is.finite(v)] <- NA_real_
  v
}
FORK$MB1_forkscale_log <- .fs_log(FORK$MB1_pc1, FORK$MB1_index)
FORK$MB2_forkscale_log <- .fs_log(FORK$MB2_pc1, FORK$MB2_index)
message("   log form: ", sum(is.na(FORK$MB1_forkscale_log)), " MB1 and ",
        sum(is.na(FORK$MB2_forkscale_log)), " MB2 sample(s) set to NA at ",
        "index == 1 (Inf, not NA, and would otherwise survive complete.cases)")
message("   MB1 forkscale range [", paste(round(range(FORK$MB1_forkscale), 1),
        collapse = ", "), "], MB2 [",
        paste(round(range(FORK$MB2_forkscale), 1), collapse = ", "), "]")

# --- 1.1 cross-check against the assembled frame, on its 849 shared rows ------
ALLB <- .load_one(PATH_ALL, "TCGA.all.biclusters.RNAseq.df2")
xchk <- NULL
if (all(c("MB1.forkscale", "MB2.forkscale") %in% names(ALLB))) {
  key <- intersect(names(ALLB), c("X", "sample", "Sample"))[1]
  if (!is.na(key)) {
    ia <- match(ALLB[[key]], FORK$aliquot)
    ok <- !is.na(ia)
    xchk <- tibble::tibble(
      n_shared  = sum(ok),
      max_abs_d_MB1 = max(abs(FORK$MB1_forkscale[ia[ok]] -
                                ALLB$MB1.forkscale[ok]), na.rm = TRUE),
      max_abs_d_MB2 = max(abs(FORK$MB2_forkscale[ia[ok]] -
                                ALLB$MB2.forkscale[ok]), na.rm = TRUE))
    message("   cross-check vs the assembled 849-row frame: ", xchk$n_shared,
            " shared, max |diff| MB1 ", signif(xchk$max_abs_d_MB1, 3),
            ", MB2 ", signif(xchk$max_abs_d_MB2, 3))
  }
}
rm(ALLB)

# =============================================================================
# 2. This repo's TCGA quantities. Read, never re-scored.
# =============================================================================
message("\n2. TCGA quantities, read from existing objects")

mito <- readRDS(PATH_TCGA_MITO)
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
cv   <- readRDS(PATH_TCGA_COV)$covariates
my   <- readRDS(PATH_TCGA_MYC)$estimators
ID_T <- colnames(mito$gsva_arms)
stopifnot(length(ID_T) == EXPECT_TCGA_SAMPLES)

ARM <- "OXPHOS subunits"          # the narrow set, NOT the umbrella
PROLIF_COV <- "PROLIF_DISJOINT"
MYC_PROLIFSTRIP <- "FELSHER__PROLIFSTRIP"

i <- match(ID_T, cv$patient); j <- match(ID_T, my$patient)
TC <- tibble::tibble(
  patient = ID_T,
  ox_gsva = as.numeric(mito$gsva_arms[ARM, ID_T]),
  ox_ppd  = as.numeric(mito$mitopps_arms[ARM, ID_T]),
  M_a     = as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T]),
  M_b     = as.numeric(nw$tcga_M_b_variants[MB_REF, ID_T]),
  M_c     = as.numeric(my$M_c_call[j]),
  M_a_ps  = as.numeric(nw$tcga_gsva_new[MYC_PROLIFSTRIP, ID_T]),
  PROLIF  = as.numeric(mito$gsva_cov[PROLIF_COV, ID_T]),
  er_call = cv$er_call[i],
  block_c = cv$complete_block_c[i])
TC$ER <- factor(dplyr::recode(TC$er_call, Positive = "ERpos",
                              Negative = "ERneg", .default = NA_character_),
                levels = c("ERpos", "ERneg"))
message("   ", nrow(TC), " TCGA patients; M_c (GISTIC call) has ",
        sum(is.na(TC$M_c)), " NA - the only estimator that does")
message("   Block C complete: ", sum(TC$block_c, na.rm = TRUE))

# --- STATE, read-only from the FROZEN repo. Context only, never a test. -------
STATE <- NULL
if (file.exists(PATH_STATE)) {
  st <- readRDS(PATH_STATE)
  STATE <- tibble::tibble(patient = st$tcga$state$patient,
                          STATE = as.character(st$tcga$state$STATE_gsva))
  message("   STATE read from the frozen repo (", nrow(STATE),
          " patients) - P3 context, NOT a test")
} else {
  message("   STATE object absent - the P3 context table is skipped, and the ",
          "note must say so rather than quietly omitting it")
}

# =============================================================================
# 3. The join, and the declared stop
# =============================================================================
message("\n3. the join")

D <- FORK %>%
  dplyr::inner_join(TC, by = "patient")
message("   forkscale patients: ", nrow(FORK), "; joining to this repo's ",
        nrow(TC), " TCGA patients: ", nrow(D), " (",
        round(100 * nrow(D) / nrow(FORK), 1), "%)")
if (!is.null(STATE)) D <- D %>% dplyr::left_join(STATE, by = "patient")

n_block_c <- sum(D$block_c, na.rm = TRUE)
message("   of those, Block C complete: ", n_block_c,
        " (declaration expects at least ", EXPECT_OVERLAP_AT_LEAST, ")")
if (n_block_c < MIN_OVERLAP) {
  stop("the forkscale-to-Block-C overlap is ", n_block_c, ", below the ",
       "declared floor of ", MIN_OVERLAP, ". STOP AND REPORT rather than ",
       "proceed on a cohort nobody has explained.", call. = FALSE)
}
if (n_block_c < EXPECT_OVERLAP_AT_LEAST) {
  message("   NOTE: below the expected ", EXPECT_OVERLAP_AT_LEAST,
          " but above the ", MIN_OVERLAP, " stop. Reported, not silently used.")
}

# THE ANALYSIS SET. Everything below runs on the full forkscale-by-TCGA join,
# which is the widest set on which both axes exist. Block C is reported as a
# sensitivity, because it is the set Block C's estimates live on and the
# declaration asks for the overlap exactly.
IDS_ALL <- D$patient
IDS_BC  <- D$patient[which(D$block_c)]

# The S7A diagonal, measured rather than assumed.
rho_forks <- stats::cor(D$MB1_forkscale, D$MB2_forkscale, method = "spearman")
message("   MB1 vs MB2 forkscale, Spearman: ", round(rho_forks, 3),
        " - the Supp. S7A diagonal, and the reason the PARTIAL is the ",
        "discriminating quantity")

# =============================================================================
# 4. Partial Spearman, with a bootstrap that includes the ranking
# =============================================================================
message("\n4. partial Spearman, ", format(N_BOOT, big.mark = ","),
        " bootstrap resamples")

# Rank, then residualise on the ranked covariates, then correlate the residuals.
# This is E18's `.cor_block` logic and ppcor::pcor.test(method = "spearman")'s,
# written out so the bootstrap can re-rank inside every resample rather than
# ranking once outside it.
.partial_spearman <- function(x, y, Z = NULL) {
  ok <- stats::complete.cases(x, y, Z)
  if (sum(ok) < 10L) return(NA_real_)
  rx <- rank(x[ok]); ry <- rank(y[ok])
  if (is.null(Z)) return(stats::cor(rx, ry))
  RZ <- apply(as.matrix(Z)[ok, , drop = FALSE], 2L, rank)
  H  <- qr(cbind(1, RZ))
  stats::cor(qr.resid(H, rx), qr.resid(H, ry))
}
.n_complete <- function(x, y, Z = NULL) sum(stats::complete.cases(x, y, Z))

# One resample index set is shared by every statistic in a block, so the
# intervals and the differences between them are internally consistent.
set.seed(PROJECT_SEED)
BOOT_IX <- replicate(N_BOOT, sample.int(nrow(D), replace = TRUE), simplify = FALSE)

.ci <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, c((1 - CI_LEVEL) / 2, 1 - (1 - CI_LEVEL) / 2),
                         na.rm = TRUE))
}

# `spec` is a list of one-row descriptions; each carries the columns to use.
.run_block <- function(dat, boot_ix, x_col, y_col, z_cols, label, kind) {
  x <- dat[[x_col]]; y <- dat[[y_col]]
  Z <- if (length(z_cols)) as.matrix(dat[, z_cols, drop = FALSE]) else NULL
  est <- .partial_spearman(x, y, Z)
  bs  <- vapply(boot_ix, function(ix) {
    Zi <- if (is.null(Z)) NULL else Z[ix, , drop = FALSE]
    .partial_spearman(x[ix], y[ix], Zi)
  }, numeric(1))
  ci <- .ci(bs)
  # EVERYTHING IS COMPUTED BEFORE THE TIBBLE. Inside tibble() the new columns
  # `x` and `y` would shadow these vectors - data masking evaluates arguments
  # in order against the columns already created - and the helpers would then
  # silently receive the column NAMES instead of the data. The columns are also
  # named x_var / y_var so the shadowing cannot come back.
  n_cc <- .n_complete(x, y, Z)
  r_p  <- suppressWarnings(stats::cor(x, y, use = "complete.obs"))
  tibble::tibble(
    kind = kind, label = label, x_var = x_col, y_var = y_col,
    adjusted_for = if (length(z_cols)) paste(z_cols, collapse = " + ") else "-",
    n = n_cc, rho = est, ci_lo = ci[1], ci_hi = ci[2],
    ci_excludes_0 = isTRUE(ci[1] > 0) || isTRUE(ci[2] < 0),
    r_pearson = r_p)
}

SPEC <- list()
.add <- function(...) SPEC[[length(SPEC) + 1L]] <<- list(...)

# --- 4.1 THE DECLARED PRIMARY: the two partials, on all three estimators -----
for (m in MYC_ESTIMATORS) {
  .add(x = m, y = "MB2_forkscale", z = "MB1_forkscale",
       label = paste0(m, " ~ MB2 | MB1"), kind = "PRIMARY (declared)")
  .add(x = m, y = "MB1_forkscale", z = "MB2_forkscale",
       label = paste0(m, " ~ MB1 | MB2"), kind = "PRIMARY (declared)")
  # marginals, reported alongside as the declaration asks
  .add(x = m, y = "MB2_forkscale", z = character(0),
       label = paste0(m, " ~ MB2 marginal"), kind = "marginal (declared)")
  .add(x = m, y = "MB1_forkscale", z = character(0),
       label = paste0(m, " ~ MB1 marginal"), kind = "marginal (declared)")
}

# --- 4.2 DECLARED UNINFORMATIVE: OXPHOS against either fork ------------------
# Computed and reported. NOT presented as support in either direction, whatever
# it returns, because both biclusters carry biogenesis. Declaration P2.
for (r in OX_RULERS) {
  .add(x = r, y = "MB2_forkscale", z = "MB1_forkscale",
       label = paste0(r, " ~ MB2 | MB1"), kind = "UNINFORMATIVE (declared)")
  .add(x = r, y = "MB1_forkscale", z = "MB2_forkscale",
       label = paste0(r, " ~ MB1 | MB2"), kind = "UNINFORMATIVE (declared)")
  .add(x = r, y = "MB2_forkscale", z = character(0),
       label = paste0(r, " ~ MB2 marginal"), kind = "UNINFORMATIVE (declared)")
  .add(x = r, y = "MB1_forkscale", z = character(0),
       label = paste0(r, " ~ MB1 marginal"), kind = "UNINFORMATIVE (declared)")
}

# --- 4.3 THE PROLIFERATION COMPANION. POST-HOC. CANNOT CHANGE THE VERDICT. ---
.add(x = "M_a", y = "MB2_forkscale", z = c("MB1_forkscale", "PROLIF"),
     label = "M_a ~ MB2 | MB1 + PROLIF", kind = "companion (POST-HOC)")
.add(x = "M_a", y = "MB1_forkscale", z = c("MB2_forkscale", "PROLIF"),
     label = "M_a ~ MB1 | MB2 + PROLIF", kind = "companion (POST-HOC)")
# The same question asked without the covariate, on the proliferation-stripped
# estimator. CLAUDE.md: use __PROLIFSTRIP to ask the question cleanly.
.add(x = "M_a_ps", y = "MB2_forkscale", z = "MB1_forkscale",
     label = "M_a__PROLIFSTRIP ~ MB2 | MB1", kind = "companion (POST-HOC)")
.add(x = "M_a_ps", y = "MB1_forkscale", z = "MB2_forkscale",
     label = "M_a__PROLIFSTRIP ~ MB1 | MB2", kind = "companion (POST-HOC)")
# And what the forkscales are, against proliferation itself.
.add(x = "PROLIF", y = "MB2_forkscale", z = character(0),
     label = "PROLIF ~ MB2 marginal", kind = "companion (POST-HOC)")
.add(x = "PROLIF", y = "MB1_forkscale", z = character(0),
     label = "PROLIF ~ MB1 marginal", kind = "companion (POST-HOC)")

# --- 4.4 the log form, reported alongside only ------------------------------
.add(x = "M_a", y = "MB2_forkscale_log", z = "MB1_forkscale_log",
     label = "M_a ~ MB2 | MB1 (log form)", kind = "log form (alongside)")
.add(x = "M_a", y = "MB1_forkscale_log", z = "MB2_forkscale_log",
     label = "M_a ~ MB1 | MB2 (log form)", kind = "log form (alongside)")

COR <- dplyr::bind_rows(lapply(SPEC, function(s)
  .run_block(D, BOOT_IX, s$x, s$y, s$z, s$label, s$kind)))

message("   ", nrow(COR), " correlations computed")

# --- 4.5 the Block C sensitivity on the declared primary --------------------
D_BC <- D[which(D$block_c), , drop = FALSE]
set.seed(PROJECT_SEED + 1L)
BOOT_BC <- replicate(N_BOOT, sample.int(nrow(D_BC), replace = TRUE),
                     simplify = FALSE)
COR_BC <- dplyr::bind_rows(lapply(
  list(list(x = "M_a", y = "MB2_forkscale", z = "MB1_forkscale",
            label = "M_a ~ MB2 | MB1 (Block C)"),
       list(x = "M_a", y = "MB1_forkscale", z = "MB2_forkscale",
            label = "M_a ~ MB1 | MB2 (Block C)")),
  function(s) .run_block(D_BC, BOOT_BC, s$x, s$y, s$z, s$label,
                         "Block C sensitivity")))
COR <- dplyr::bind_rows(COR, COR_BC)

# =============================================================================
# 5. THE VERDICT, on the rule fixed at a71e19d
# =============================================================================
message("\n5. the verdict\n", strrep("-", 78))

.get <- function(lab) COR[match(lab, COR$label), ]
p2 <- .get("M_a ~ MB2 | MB1")
p1 <- .get("M_a ~ MB1 | MB2")

# The paired difference, bootstrapped on the SAME resamples, so it is a genuine
# paired contrast and not a comparison of two independent intervals. Reported
# as the stricter reading; the VERDICT uses the declaration's literal wording.
d_boot <- vapply(BOOT_IX, function(ix) {
  dd <- D[ix, , drop = FALSE]
  .partial_spearman(dd$M_a, dd$MB2_forkscale, as.matrix(dd[, "MB1_forkscale"])) -
    .partial_spearman(dd$M_a, dd$MB1_forkscale, as.matrix(dd[, "MB2_forkscale"]))
}, numeric(1))
d_est <- p2$rho - p1$rho
d_ci  <- .ci(d_boot)

.overlap <- function(a, b) !(a$ci_hi < b$ci_lo || b$ci_hi < a$ci_lo)
cis_overlap <- .overlap(p2, p1)

verdict <- if (p2$rho > 0 && isTRUE(p2$ci_lo > 0) && p2$rho > p1$rho &&
               !cis_overlap) {
  "MB2-ALIGNED"
} else if (p1$rho > 0 && isTRUE(p1$ci_lo > 0) && p1$rho > p2$rho &&
           !cis_overlap) {
  "MB1-ALIGNED"
} else {
  "UNRESOLVED"
}

message("RULE: ", VERDICT_RULE)
message("\n   M_a ~ MB2 | MB1 : rho = ", round(p2$rho, 3), "  [",
        round(p2$ci_lo, 3), ", ", round(p2$ci_hi, 3), "]  n = ", p2$n)
message("   M_a ~ MB1 | MB2 : rho = ", round(p1$rho, 3), "  [",
        round(p1$ci_lo, 3), ", ", round(p1$ci_hi, 3), "]  n = ", p1$n)
message("   difference (paired bootstrap): ", round(d_est, 3), "  [",
        round(d_ci[1], 3), ", ", round(d_ci[2], 3), "]")
message("   the two CIs overlap: ", cis_overlap)
message("\n   VERDICT: ", verdict)

reading <- switch(
  verdict,
  "MB2-ALIGNED" = paste0(
    "MB2-ALIGNED. MYC activity tracks the MB2 forkscale at fixed MB1 and not\n",
    "   the reverse. On the declaration this licenses ONE sentence: that the\n",
    "   state described here corresponds to the MYC-associated bicluster\n",
    "   rather than to the proliferation-driven one whose prognostic\n",
    "   behaviour was previously reported.\n",
    "   READ SECTION 6 BEFORE WRITING IT. MB2 forkscale is close to a\n",
    "   proliferation axis, and the companion says what survives adjustment."),
  "MB1-ALIGNED" = paste0(
    "MB1-ALIGNED - the INFORMATIVE FAILURE, and it is reported as a\n",
    "   discrepancy rather than worked around. This arm's axis would be the\n",
    "   programme Menegollo already characterised, the published description\n",
    "   of MB1 as Myc-independent would be in tension with our data, and\n",
    "   limb (i) would cite them as prior work rather than as a contrast."),
  "UNRESOLVED" = paste0(
    "UNRESOLVED, which the declaration names as an acceptable outcome.\n",
    "   The manuscript sentence does not get written and limb (i)'s Menegollo\n",
    "   clause is DROPPED rather than softened. The rule is not relaxed after\n",
    "   the fact to reach a conclusion."))
message("\n", reading)

# =============================================================================
# 6. The proliferation companion. POST-HOC. Reported, never decisive.
# =============================================================================
message("\n6. the proliferation companion (POST-HOC - cannot change section 5)")
message("   Declared at fb98f33 because MB2 forkscale tracks PROLIF_DISJOINT ",
        "harder\n   than it tracks MYC. If the partial dies here, 'MYC tracks ",
        "MB2' is\n   'proliferation tracks MB2' - which is close to ",
        "Menegollo's own MB1.")
COR %>%
  dplyr::filter(kind == "companion (POST-HOC)") %>%
  dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3), ci_excludes_0) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. Declared UNINFORMATIVE, and the other estimators
# =============================================================================
message("\n7. DECLARED UNINFORMATIVE IN ADVANCE - both biclusters carry ",
        "biogenesis.\n   Reported; NOT support in either direction whatever ",
        "it says.")
COR %>%
  dplyr::filter(kind == "UNINFORMATIVE (declared)") %>%
  dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the declared primary on all three estimators:")
COR %>%
  dplyr::filter(kind %in% c("PRIMARY (declared)", "marginal (declared)")) %>%
  dplyr::transmute(kind, label, n, rho = round(rho, 3),
                   ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                   ci_excludes_0) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. The F3-pre cross-check. Reported, NOT reconciled.
# =============================================================================
# Different repos' scoring. If the MB1 marginal here differs materially from
# the stored one, that IS the finding of this section - it is not repaired by
# replacing one with the other. Declaration section 4.
message("\n8. cross-check against F3-pre (a DIFFERENT repo's scoring)")

f3 <- NULL
if (file.exists(PATH_F3PRE)) {
  f3raw <- readRDS(PATH_F3PRE)$f3pre$correlations
  want <- tibble::tribble(
    ~fork,            ~axis,                       ~here,
    "MB1_forkscale",  "MYC (M_a)",                 "M_a ~ MB1 marginal",
    "MB2_forkscale",  "MYC (M_a)",                 "M_a ~ MB2 marginal",
    "MB1_forkscale",  "MYC (M_b)",                 "M_b ~ MB1 marginal",
    "MB2_forkscale",  "MYC (M_b)",                 "M_b ~ MB2 marginal",
    "MB1_forkscale",  "OXPHOS subunits (GSVA)",    "ox_gsva ~ MB1 marginal",
    "MB2_forkscale",  "OXPHOS subunits (GSVA)",    "ox_gsva ~ MB2 marginal",
    "MB1_forkscale",  "OXPHOS subunits (mitoPPS)", "ox_ppd ~ MB1 marginal",
    "MB2_forkscale",  "OXPHOS subunits (mitoPPS)", "ox_ppd ~ MB2 marginal",
    "MB1_forkscale",  "PROLIF_DISJOINT",           "PROLIF ~ MB1 marginal",
    "MB2_forkscale",  "PROLIF_DISJOINT",           "PROLIF ~ MB2 marginal")
  f3 <- want %>%
    dplyr::left_join(f3raw %>% dplyr::select(fork, axis, n_f3 = n,
                                             rho_f3 = rho_spearman),
                     by = c("fork", "axis")) %>%
    dplyr::left_join(COR %>% dplyr::select(here = label, n_here = n,
                                           rho_here = rho),
                     by = "here") %>%
    dplyr::mutate(delta = rho_here - rho_f3)
  f3 %>%
    dplyr::transmute(fork, axis, n_f3, rho_f3 = round(rho_f3, 3),
                     n_here, rho_here = round(rho_here, 3),
                     delta = round(delta, 3)) %>%
    as.data.frame() %>% print(row.names = FALSE)
  message("   max |delta| = ", round(max(abs(f3$delta), na.rm = TRUE), 3),
          ". REPORTED, NOT RECONCILED - these are two repos' scoring and ",
          "neither replaces the other.")
  message("   NOTE: F3-pre stored NO M_c row for any fork, so M_c has no ",
          "cross-check.")
} else {
  message("   forkscale_models.rds absent - no cross-check. Say so in the note.")
}

# =============================================================================
# 9. P3 - CONTEXT, NOT A TEST. ER status by STATE and by forkscale quadrant.
# =============================================================================
# MB2_UF is ER-negative only in Menegollo's description. Consistency is worth
# recording; it is NEVER evidence, and no test statistic is computed here.
message("\n9. P3 CONTEXT, NOT A TEST - ER status. No test statistic is ",
        "computed.")

er_by_state <- NULL
if (!is.null(STATE)) {
  er_by_state <- D %>%
    dplyr::filter(!is.na(STATE), !is.na(ER)) %>%
    dplyr::count(STATE, ER, name = "n") %>%
    tidyr::pivot_wider(names_from = ER, values_from = n, values_fill = 0L) %>%
    dplyr::mutate(n = ERpos + ERneg, pct_ERneg = round(100 * ERneg / n, 1))
  message("\n   ER by STATE level:")
  er_by_state %>% as.data.frame() %>% print(row.names = FALSE)
}

# Quadrants on the within-cohort medians of the two forkscales.
D$quadrant <- paste0(
  ifelse(D$MB1_forkscale > stats::median(D$MB1_forkscale), "MB1hi", "MB1lo"),
  "_",
  ifelse(D$MB2_forkscale > stats::median(D$MB2_forkscale), "MB2hi", "MB2lo"))
er_by_quad <- D %>%
  dplyr::filter(!is.na(ER)) %>%
  dplyr::count(quadrant, ER, name = "n") %>%
  tidyr::pivot_wider(names_from = ER, values_from = n, values_fill = 0L) %>%
  dplyr::mutate(n = ERpos + ERneg, pct_ERneg = round(100 * ERneg / n, 1))
message("\n   ER by forkscale quadrant (medians within this cohort):")
er_by_quad %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 10. The figure
# =============================================================================
# Every statistic above is rank-based, so the rank-rank panel is the faithful
# picture. The raw panel is shown beside it because forkscale's skew is a real
# property of the axis and hiding it behind ranks would be a different kind of
# dishonesty: MB1 spans roughly [-42.6, 8.8] and a raw scatter is dominated by
# a handful of low-index samples.
message("\n10. figure")

D$myc_tertile <- cut(D$M_a, stats::quantile(D$M_a, c(0, 1/3, 2/3, 1),
                                            na.rm = TRUE),
                     labels = c("MYC T1 (low)", "MYC T2", "MYC T3 (high)"),
                     include.lowest = TRUE)
PD <- D %>% dplyr::filter(!is.na(myc_tertile)) %>%
  dplyr::mutate(ER_shape = factor(ifelse(is.na(ER), "unknown",
                                         as.character(ER)),
                                  levels = c("ERpos", "ERneg", "unknown")),
                MB1_rank = rank(MB1_forkscale), MB2_rank = rank(MB2_forkscale))

COLS <- c("MYC T1 (low)" = "#4575b4", "MYC T2" = "#bdbdbd",
          "MYC T3 (high)" = "#b2182b")
SHAPES <- c(ERpos = 16, ERneg = 17, unknown = 4)

.base <- function(p) p +
  ggplot2::scale_colour_manual(values = COLS, name = "M_a tertile") +
  ggplot2::scale_shape_manual(values = SHAPES, name = "ER status") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom")

p_rank <- .base(ggplot2::ggplot(
  PD, ggplot2::aes(MB1_rank, MB2_rank, colour = myc_tertile,
                   shape = ER_shape)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.75)) +
  ggplot2::labs(subtitle = paste0("rank-rank - the scale every statistic ",
                                  "here is computed on\n(Spearman rho(MB1, ",
                                  "MB2) = ", round(rho_forks, 3), ")"),
                x = "MB1 forkscale, rank", y = "MB2 forkscale, rank")

# The raw values plotted raw are unreadable: pc1 / index with index near 1
# puts almost every sample on top of the origin and a handful of low-index
# samples at the edges. A signed log DISPLAY TRANSFORM keeps the sign and the
# skew visible without crushing the middle. It changes no statistic - every
# number in this script is rank-based - and the untransformed ranges are in
# the caption.
.slog <- function(v) sign(v) * log10(1 + abs(v))
PD$MB1_slog <- .slog(PD$MB1_forkscale)
PD$MB2_slog <- .slog(PD$MB2_forkscale)
p_raw <- .base(ggplot2::ggplot(
  PD, ggplot2::aes(MB1_slog, MB2_slog, colour = myc_tertile,
                   shape = ER_shape)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.75)) +
  ggplot2::labs(subtitle = paste0("raw forkscale = pc1 / index, on a signed ",
                                  "log display scale -\nseverely skewed BY ",
                                  "CONSTRUCTION, which is why Spearman is used"),
                x = "MB1 forkscale, sign(x) * log10(1 + |x|)",
                y = "MB2 forkscale, sign(x) * log10(1 + |x|)")

p <- patchwork::wrap_plots(p_rank, p_raw, nrow = 1, guides = "collect") +
  patchwork::plot_annotation(
    title = paste0("MB1 against MB2 forkscale, ", nrow(PD),
                   " TCGA-BRCA patients"),
    subtitle = paste0(
      "VERDICT on the rule fixed at a71e19d before any number: ", verdict,
      "\nM_a ~ MB2 | MB1 = ", round(p2$rho, 3), " [", round(p2$ci_lo, 3), ", ",
      round(p2$ci_hi, 3), "]   vs   M_a ~ MB1 | MB2 = ", round(p1$rho, 3),
      " [", round(p1$ci_lo, 3), ", ", round(p1$ci_hi, 3), "]",
      "\nEXPLORATORY. Not pre-registered. ER is CONTEXT, never evidence."),
    caption = paste0(
      "Forkscale recomputed as pc1 / index from the pinned Menegollo snapshot ",
      "(upstream 8fdbb34). Partial Spearman, ", format(N_BOOT, big.mark = ","),
      " bootstrap resamples, ranking inside every resample. Raw ranges: MB1 [",
      paste(round(range(PD$MB1_forkscale), 1), collapse = ", "), "], MB2 [",
      paste(round(range(PD$MB2_forkscale), 1), collapse = ", "),
      "]; the right panel's signed log is a DISPLAY transform only.\n",
      "OXPHOS against either forkscale was declared UNINFORMATIVE in advance ",
      "and is not shown here. No outcome variable enters this analysis."),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 6, hjust = 0),
      plot.title.position = "plot"))

ggplot2::ggsave(PATH_E32_FIG, p, width = 11, height = 6.2, dpi = 200)
.ensure_dir(dirname(PATH_E32_FIG_DOC))
ggplot2::ggsave(PATH_E32_FIG_DOC, p, width = 11, height = 6.2, dpi = 200)
message("   ", PATH_E32_FIG)
message("   ", PATH_E32_FIG_DOC, "  (tracked - outputs/ is gitignored)")

# =============================================================================
# 11. Save
# =============================================================================
message("\n11. saving")

out <- list(
  correlations = COR,
  verdict      = verdict,
  reading      = reading,
  verdict_inputs = tibble::tibble(
    quantity = c("M_a ~ MB2 | MB1", "M_a ~ MB1 | MB2", "difference (paired)"),
    rho = c(p2$rho, p1$rho, d_est),
    ci_lo = c(p2$ci_lo, p1$ci_lo, d_ci[1]),
    ci_hi = c(p2$ci_hi, p1$ci_hi, d_ci[2]),
    n = c(p2$n, p1$n, nrow(D))),
  cis_overlap  = cis_overlap,
  forkscale    = FORK,
  f3pre_crosscheck = f3,
  assembled_crosscheck = xchk,
  er_by_state  = er_by_state,
  er_by_quadrant = er_by_quad,
  rho_mb1_mb2  = rho_forks,
  counts = list(
    fork_patients = nrow(FORK),
    tcga_patients = nrow(TC),
    joined        = nrow(D),
    block_c_overlap = n_block_c,
    m_c_na        = sum(is.na(D$M_c)),
    er_na         = sum(is.na(D$ER))),
  spec = list(
    declaration = "docs/2026-09-28_mb2_forkscale_declaration.md, committed ALONE at a71e19d",
    data_note   = "docs/2026-09-28_mb2_forkscale_data.md, fb98f33",
    posture = paste("EXPLORATORY, not pre-registered. The verdict rule was",
                    "fixed BEFORE forkscale_models.rds was opened."),
    rule = VERDICT_RULE,
    already_known = paste("The MARGINALS were already in forkscale_models.rds,",
                          "computed 2026-08-30 (MB2 vs M_a +0.739, MB1 +0.429).",
                          "Unread until after a71e19d. What this adds is the",
                          "PARTIALS, this repo's scoring, bootstrap CIs and M_c."),
    uninformative = paste("OXPHOS vs either forkscale is DECLARED",
                          "UNINFORMATIVE - both biclusters carry biogenesis.",
                          "Never presented as support in either direction."),
    companion = paste("The proliferation companion is POST-HOC, chosen after",
                      "reading F3-pre, declared at fb98f33 before the script",
                      "was written. It CANNOT change the verdict."),
    forkscale = "pc1 / index, plain form; the log form is reported alongside only",
    log_inf = "forkscale.log is Inf at index == 1 and is set to NA explicitly",
    estimator = paste("partial Spearman: rank, residualise on ranked",
                      "covariates, correlate residuals. Bootstrap re-ranks",
                      "INSIDE every resample."),
    n_boot = N_BOOT,
    m_c = paste("M_c is GISTIC amplification (M_c_call), a discrete call with",
                "NAs - not MYC mRNA. F3-pre stored no M_c cross-check."),
    state = paste("STATE read read-only from the frozen validation repo;",
                  "this repo does not build it. P3 CONTEXT ONLY."),
    no_outcome = "NO outcome or survival variable enters at any point",
    no_fdr = "no per-gene FDR, consistent with the rest of the project",
    n3 = "two exposure axes and a transcript score; nothing is 'primed'",
    seed = PROJECT_SEED),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E32)
utils::write.csv(as.data.frame(COR), PATH_E32_COR, row.names = FALSE)
message("   ", PATH_E32)
message("   ", PATH_E32_COR)
message("\nE32 done.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E32)

  cat(x$verdict, "\n\n"); cat(x$reading, "\n")
  utils::str(x$counts)

  # The two quantities the verdict rests on, and their paired difference.
  x$verdict_inputs %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 3))) %>%
    as.data.frame()
  x$cis_overlap

  # The declared primary on all three estimators, partials and marginals.
  x$correlations %>%
    dplyr::filter(kind %in% c("PRIMARY (declared)", "marginal (declared)")) %>%
    dplyr::transmute(kind, label, n, rho = round(rho, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                     ci_excludes_0) %>%
    as.data.frame()

  # THE COMPANION. Post-hoc. If the partial dies here, "MYC tracks MB2" is
  # "proliferation tracks MB2".
  x$correlations %>%
    dplyr::filter(kind == "companion (POST-HOC)") %>%
    dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                     ci_hi = round(ci_hi, 3), ci_excludes_0) %>%
    as.data.frame()

  # Declared UNINFORMATIVE. Not support in either direction.
  x$correlations %>%
    dplyr::filter(kind == "UNINFORMATIVE (declared)") %>%
    dplyr::transmute(label, rho = round(rho, 3)) %>%
    as.data.frame()

  # The cross-check against F3-pre. Two repos' scoring; not reconciled.
  x$f3pre_crosscheck %>%
    dplyr::transmute(fork, axis, rho_f3 = round(rho_f3, 3),
                     rho_here = round(rho_here, 3), delta = round(delta, 3)) %>%
    as.data.frame()

  # P3 context. NEVER evidence.
  x$er_by_state %>% as.data.frame()
  x$er_by_quadrant %>% as.data.frame()
}
