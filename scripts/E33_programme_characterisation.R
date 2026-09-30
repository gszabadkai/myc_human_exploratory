# E33_programme_characterisation.R
# =============================================================================
# CHARACTERISING THE MB1 / MB2 DISTINCTION, AND WHERE THE DEATH SIGNAL SITS
#
# EXPLORATORY. Nothing here is pre-registered.
#
# The limbs, the readouts, the verdict and what is declared uninformative were
# fixed in docs/2026-09-29_E33_declaration_v2.md, COMMITTED ALONE AT 9559931
# BEFORE ANY DATA FILE WAS OPENED. The inputs, two findings that change what the
# declaration assumed, and the specification decisions it left open are in
# docs/2026-09-30_e33_data.md at 9b3ae75. Read both before this file; this
# header restates their rules, it does not invent them.
#
# =============================================================================
# THE THREE LIMBS
# =============================================================================
# LIMB 1 - DESCRIPTIVE, NO VERDICT. Two quantities E32 never reported:
#   (a) M_a, M_b, M_c and FELSHER__PROLIFSTRIP, median and IQR, in each of the
#       four quadrants formed by MEDIAN SPLITS of MB1 and MB2 forkscale.
#   (b) rho(M_a, MB1_forkscale | PROLIF_DISJOINT), beside E32's
#       rho(M_a, MB1 | MB2) = -0.195 and the marginal +0.442.
#   MB2 forkscale runs +0.874 with PROLIF_DISJOINT, so "conditioned on MB2" and
#   "conditioned on proliferation" are close to indistinguishable here. THE
#   DIFFERENCE BETWEEN THE TWO PARTIALS IS REPORTED; WHICH ONE IS DOING THE WORK
#   IS NOT INTERPRETED. Expect a gradient, not two states: MB1_UF is NOT
#   MYC-low in absolute terms and nothing here may say it is.
#
# LIMB 2 - DESCRIPTIVE, NO VERDICT. The variable map. The MYC x OXPHOS quadrant
#   from the FROZEN constructor's own calls, cross-tabulated against STATE,
#   MB1 UF/LF/none, MB2 UF/LF/none, PAM50 and ER. Counts only.
#
# LIMB 3 - THE CLAIM-RELEVANT ONE. Trigger arm against guardian arm:
#       rho( arm , MB2_forkscale | MB1_forkscale )   against
#       rho( arm , MB1_forkscale | MB2_forkscale )
#   trigger  BBC3, BID, BIK, BAD
#   guardian BCL2L1, MCL1 - together and individually
#   Partial Spearman, 10,000 bootstrap resamples, ranking redone inside each.
#   The composite configuration score is reported THIRD and is NOT what the
#   verdict turns on.
#
# =============================================================================
# FORKSCALE MEDIAN SPLITS ARE NOT UF/LF ASSIGNMENTS
# =============================================================================
# Both appear in this script and they are different variables:
#   fs_quad      a cut at the within-set median of a continuous score
#   MB{k}_fork   an MCbiclust sample assignment, UF / LF / none, where `none`
#                is every sample outside the top ~300 of that bicluster's sort
# Their names never overlap, and section 7 cross-tabulates one against the
# other so the difference is shown rather than asserted.
#
# =============================================================================
# DECISIONS TAKEN AFTER THE DECLARATION, BEFORE THIS SCRIPT (data note)
# =============================================================================
#   - ANALYSIS SET = the ALIQUOT join, n = 1,031, not 1,037. Six patients were
#     profiled from different vials; E32's patient join paired them across
#     vials. The patient join is kept to reproduce E32 exactly and as a
#     sensitivity. Data note section 3.
#   - UF/LF/none = the STORED calls (849 rows); the other 188 are an explicit
#     level `not in 849`. A bounded extension is a sensitivity. Section 2.
#   - ARMS ARE SIGNED BY THE CONFIGURATION, decided by the author 2026-09-30:
#     BBC3, BID, BIK, BAD, BCL2L1 carry +1, MCL1 carries -1. Each arm is E16's
#     `.comp` recipe on those signs. Section 6.
#   - "EXCEEDS" IS THE POINT ESTIMATE, as written. The paired difference is
#     reported as the stricter reading and is NOT the verdict. Section 7.
#   - A FOURTH OUTCOME the declaration does not name - guardian separates,
#     trigger does not - is reported literally and given no reading.
#
# =============================================================================
# DECLARED UNINFORMATIVE IN ADVANCE
# =============================================================================
# OXPHOS against either forkscale. Both biclusters carry mitochondrial
# biogenesis, so it discriminates nothing. Computed, reported, NEVER presented
# as support in either direction.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - NO outcome or survival variable enters at any point. Not joined, not read.
#     MB1_UF and MB2_UF would be expected to have poor outcome for DIFFERENT
#     reasons (proliferation, MYC-coupling), so survival could not separate them.
#   - NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM. No OXPHOS-directed or
#     BH3-directed intervention exists in any cohort available.
#   - Nothing is written to myc_human_validation or myc_mouse.
#   - No per-gene FDR. No causal language. MB3 and METABRIC are not read.
#   - STATE is not built here; the frozen constructor is read and called.
#
# N3: transcript scores and two exposure axes. Nothing is "primed".
# SPECIES: human throughout. NO ORTHOLOG FUNCTION IS CALLED.
# COHORT-RELATIVITY: one cohort (TCGA); nothing is pooled, and no value here is
# ever compared numerically with a SCAN-B, CCLE or mouse value.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE33: MB1 / MB2 characterisation, and trigger against guardian\n",
        strrep("=", 78))

PATH_E33      <- file.path(DIR_RESULTS, "e33_programme_characterisation.rds")
PATH_E33_COR  <- file.path(DIR_TABLES,  "E33_correlations.csv")
PATH_E33_FIG1 <- file.path(DIR_FIGURES, "E33_forkscale_quadrant_map.png")
PATH_E33_FIG2 <- file.path(DIR_FIGURES, "E33_trigger_vs_guardian.png")
# Written twice, as from E31 on: outputs/ is gitignored, and the note plus the
# tracked figure are the durable record.
DIR_DOC_FIGURES   <- here::here("docs", "figures")
PATH_E33_FIG1_DOC <- file.path(DIR_DOC_FIGURES,
                               "2026-09-30_E33_forkscale_quadrant_map.png")
PATH_E33_FIG2_DOC <- file.path(DIR_DOC_FIGURES,
                               "2026-09-30_E33_trigger_vs_guardian.png")

# =============================================================================
# 0. Constants and THE VERDICT RULE, transcribed from the declaration
# =============================================================================
DIR_MENEGOLLO <- here::here("data", "menegollo_biclusters")
PATH_MB1 <- file.path(DIR_MENEGOLLO, "TCGA_MB1_RNASeq_data_nonorm.RData")
PATH_MB2 <- file.path(DIR_MENEGOLLO, "TCGA_MB2_RNASeq_data_nonorm.RData")
PATH_ALL <- file.path(DIR_MENEGOLLO, "TCGA.all.biclusters.RNAseq.Rdata")
# Read-only, from the FROZEN repo. Nothing is written there, ever.
PATH_STATE <- "/Users/gs/code/myc_human_validation/results/state_definition.rds"
# E32's saved object: a CROSS-CHECK only, never an input. Named with here::here
# rather than DIR_RESULTS so it is always the repo's copy.
PATH_E32_IN <- here::here("results", "mb2_forkscale_alignment.rds")

N_BOOT   <- 10000L
CI_LEVEL <- 0.95

# Counts recorded in the data note. A mismatch means the inputs moved.
EXPECT_FORK_N      <- 1037L
EXPECT_ASSEMBLED_N <- 849L
EXPECT_ALIQUOT_N   <- 1031L
EXPECT_BLOCK_C_AL  <- 879L

MYC_PROLIFSTRIP <- "FELSHER__PROLIFSTRIP"
PROLIF_COV      <- "PROLIF_DISJOINT"
ARM_OX          <- "OXPHOS subunits"     # the narrow set, NOT the umbrella
OX_RULERS       <- c("ox_gsva", "ox_ppd")

# --- The arms, and the configuration's signs --------------------------------
# Signs from the reconciliation doc, section 3.1: "sensitisers and BCL2L1 up,
# MCL1 down" in OXPHOS-high tumours. The author chose on 2026-09-30, before
# this script existed, to sign every composite this way (data note section 6).
# An UNSIGNED guardian mean was rejected: its two genes run in opposite
# directions (+0.388 / -0.266 on OXPHOS), so it would partly cancel and be
# biased towards "does not separate" - i.e. towards the licensing verdict.
TRIGGER  <- c("BBC3", "BID", "BIK", "BAD")
GUARDIAN <- c("BCL2L1", "MCL1")
CONFIG_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)
stopifnot(setequal(names(CONFIG_SIGN), c(TRIGGER, GUARDIAN)))

# Readouts for limb 3. Column name -> display label. The two ARMS carry the
# verdict; the composite is third; the genes are reported, never decisive.
READ_ARMS  <- c(arm_trigger = "trigger arm", arm_guardian = "guardian arm")
READ_COMP  <- c(config_comp = "configuration composite")
READ_GENES <- stats::setNames(names(CONFIG_SIGN),
                              paste0("g_", names(CONFIG_SIGN)))
READOUTS   <- c(READ_ARMS, READ_COMP, READ_GENES)

# --- THE VERDICT RULE. Declaration section 3, transcribed, not reinterpreted.
# An arm SEPARATES iff rho(arm, MB2|MB1) > 0, its CI excludes zero, and its
# point estimate EXCEEDS rho(arm, MB1|MB2). "Exceeds" is the point estimate,
# as written; the declaration has no non-overlap clause (E32's had one). The
# paired difference is reported beside it as the STRICTER reading and is not
# the verdict.
#   TRIGGER MB2-SPECIFIC  trigger separates, guardian does not
#   BOTH ARMS SEPARATE    both separate - descriptive, does not isolate the trigger
#   NEITHER SEPARATES     neither - reportable, and what E11 predicts
# All three are acceptable. The fourth combination is not named by the
# declaration and is reported literally, with no reading.
VERDICT_RULE <- paste0(
  "An arm SEPARATES iff rho(arm, MB2|MB1) > 0, its 95% CI excludes zero, and ",
  "its point estimate exceeds rho(arm, MB1|MB2). TRIGGER MB2-SPECIFIC = ",
  "trigger separates, guardian does not. BOTH ARMS SEPARATE = both. NEITHER ",
  "SEPARATES = neither. The verdict turns on the two ARM composites only.")
V_TRIG    <- "TRIGGER MB2-SPECIFIC"
V_BOTH    <- "BOTH ARMS SEPARATE"
V_NEITHER <- "NEITHER SEPARATES"
V_OUTSIDE <- "OUTSIDE THE DECLARED OUTCOMES: GUARDIAN SEPARATES, TRIGGER DOES NOT"

message("\n0. the rule, fixed at 9559931 before any data was read:\n   ",
        VERDICT_RULE)

# =============================================================================
# 1. Forkscale - E32's construction, VERBATIM
# =============================================================================
# Copied from E32 section 1 without change, so limb 1 and E32 sit on the same
# quantities. Section 1.3 checks that against E32's saved object.
message("\n1. forkscale, from the per-bicluster files (E32's construction)")

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

# Menegollo's own PAM50 call, carried in the per-bicluster file. NOT the primary
# PAM50 here (data note section 5); saved for the agreement line only.
FORK$PAM50_menegollo <- MB1$PAM50

# --- 1.1 cross-check against the assembled frame, on its 849 shared rows ------
# E32 section 1.1, verbatim, except that the frame is kept: section 1.2 needs it.
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

# --- 1.2 THE UF / LF / none ASSIGNMENTS ---------------------------------------
# Located in the ASSEMBLED frame only ($MB1.fork, $MB2.fork; None/Lower/Upper),
# 849 rows. NOT the $MB{k}.forkscale.fork_0.4 / _0.8 columns, which are cuts on
# forkscale and are not used. Data note section 2.
#
# STORED CALLS ARE PRIMARY. The 188 absent from the 849 get an explicit level,
# `not in 849`, so that every cell size still sums to the analysis set.
#
# THE BOUNDED EXTENSION is a sensitivity. It does not guess the upstream rule.
# On the 849, membership is an index cut and UF/LF is a pc1 gap; a sample
# outside the 849 is assigned only where it lies beyond BOTH boundaries, where
# any rule consistent with the 849 gives the same answer. Anything else is NA.
message("\n1.2 UF / LF / none, from the assembled frame")

need <- c("X", "MB1.fork", "MB2.fork", "MB1.pc1", "MB1.index",
          "MB2.pc1", "MB2.index")
if (!all(need %in% names(ALLB))) {
  stop("the assembled frame lacks: ",
       paste(setdiff(need, names(ALLB)), collapse = ", "), call. = FALSE)
}
ia <- match(FORK$aliquot, ALLB$X)
FORK$in_849 <- !is.na(ia)
if (sum(FORK$in_849) != EXPECT_ASSEMBLED_N) {
  stop(sum(FORK$in_849), " forkscale samples in the assembled frame, not ",
       EXPECT_ASSEMBLED_N, call. = FALSE)
}

FORK_LEVELS     <- c("UF", "LF", "none", "not in 849")
FORK_EXT_LEVELS <- c("UF", "LF", "none")
FORK_BOUNDS <- list()
for (b in c("MB1", "MB2")) {
  raw <- as.character(ALLB[[paste0(b, ".fork")]])[ia]  # NA where not in 849
  if (!all(raw[!is.na(raw)] %in% c("None", "Lower", "Upper"))) {
    stop(b, ".fork has levels other than None / Lower / Upper", call. = FALSE)
  }
  lab <- unname(c(None = "none", Lower = "LF", Upper = "UF")[raw])
  ix  <- FORK[[paste0(b, "_index")]]
  pc  <- FORK[[paste0(b, "_pc1")]]
  s   <- !is.na(raw)
  # The assembled frame's pc1 and index must BE the per-bicluster values.
  if (!isTRUE(all.equal(ALLB[[paste0(b, ".pc1")]][ia[s]], pc[s])) ||
      !identical(as.integer(ALLB[[paste0(b, ".index")]][ia[s]]),
                 as.integer(ix[s]))) {
    stop(b, ": pc1 / index in the assembled frame differ from the ",
         "per-bicluster file", call. = FALSE)
  }
  a_max  <- max(ix[s & raw != "None"])     # highest index that is assigned
  n_min  <- min(ix[s & raw == "None"])     # lowest index that is none
  lo_max <- max(pc[s & raw == "Lower"])    # the LF side of the pc1 gap
  up_min <- min(pc[s & raw == "Upper"])    # the UF side
  if (a_max >= n_min) stop(b, ": membership is not an index cut", call. = FALSE)
  if (lo_max >= up_min) stop(b, ": UF/LF is not a pc1 threshold", call. = FALSE)

  stored <- ifelse(s, lab, "not in 849")
  ext <- ifelse(s, lab,
         ifelse(ix >= n_min, "none",
         ifelse(ix <= a_max,
                ifelse(pc <= lo_max, "LF",
                ifelse(pc >= up_min, "UF", NA_character_)),
                NA_character_)))
  FORK[[paste0(b, "_fork")]]     <- factor(stored, levels = FORK_LEVELS)
  FORK[[paste0(b, "_fork_ext")]] <- factor(ext, levels = FORK_EXT_LEVELS)

  n_undet <- sum(is.na(ext))
  FORK_BOUNDS[[b]] <- tibble::tibble(
    bicluster = b, assigned_index_max = a_max, none_index_min = n_min,
    LF_pc1_max = lo_max, UF_pc1_min = up_min,
    stored_UF = sum(stored == "UF"), stored_LF = sum(stored == "LF"),
    stored_none = sum(stored == "none"),
    not_in_849 = sum(stored == "not in 849"),
    ext_UF = sum(ext == "UF", na.rm = TRUE),
    ext_LF = sum(ext == "LF", na.rm = TRUE),
    ext_none = sum(ext == "none", na.rm = TRUE),
    ext_undetermined = n_undet)
  message("   ", b, ": stored UF ", sum(stored == "UF"), ", LF ",
          sum(stored == "LF"), ", none ", sum(stored == "none"),
          ", not in 849 ", sum(stored == "not in 849"),
          " | index cut ", a_max, " / ", n_min, ", pc1 gap [",
          round(lo_max, 3), ", ", round(up_min, 3), "] | extension leaves ",
          n_undet, " undetermined")
}
FORK_BOUNDS <- dplyr::bind_rows(FORK_BOUNDS)
rm(ALLB)

# --- 1.3 cross-check against E32's saved forkscale ---------------------------
e32 <- if (file.exists(PATH_E32_IN)) readRDS(PATH_E32_IN) else NULL
fs_vs_e32 <- NULL
if (!is.null(e32)) {
  f32 <- e32$forkscale
  fs_vs_e32 <- tibble::tibble(
    same_rows      = identical(f32$aliquot, FORK$aliquot),
    MB1_identical  = identical(f32$MB1_forkscale, FORK$MB1_forkscale),
    MB2_identical  = identical(f32$MB2_forkscale, FORK$MB2_forkscale),
    MB1_log_identical = identical(f32$MB1_forkscale_log, FORK$MB1_forkscale_log),
    MB2_log_identical = identical(f32$MB2_forkscale_log, FORK$MB2_forkscale_log))
  message("   against E32's saved forkscale: ",
          if (all(unlist(fs_vs_e32))) "IDENTICAL on every row and column" else
            "DIFFERS - see $forkscale_vs_e32, and the note must say so")
} else {
  message("   E32's saved object is absent; section 4.1 checks against the ",
          "values quoted in E32's note instead")
}

# =============================================================================
# 2. This repo's TCGA quantities, and the arm composites
# =============================================================================
message("\n2. TCGA quantities")

mito <- readRDS(PATH_TCGA_MITO)
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
cv   <- readRDS(PATH_TCGA_COV)$covariates
my   <- readRDS(PATH_TCGA_MYC)$estimators
ID_T <- colnames(mito$gsva_arms)
stopifnot(length(ID_T) == EXPECT_TCGA_SAMPLES)

# E32's TC block, plus the aliquot key, PAM50 and the arms.
i <- match(ID_T, cv$patient); j <- match(ID_T, my$patient)
TC <- tibble::tibble(
  patient = ID_T,
  aliquot_tcga = cv$aliquot_barcode[i],
  ox_gsva = as.numeric(mito$gsva_arms[ARM_OX, ID_T]),
  ox_ppd  = as.numeric(mito$mitopps_arms[ARM_OX, ID_T]),
  M_a     = as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T]),
  M_b     = as.numeric(nw$tcga_M_b_variants[MB_REF, ID_T]),
  M_c     = as.numeric(my$M_c_call[j]),
  M_a_ps  = as.numeric(nw$tcga_gsva_new[MYC_PROLIFSTRIP, ID_T]),
  PROLIF  = as.numeric(mito$gsva_cov[PROLIF_COV, ID_T]),
  er_call = cv$er_call[i],
  pam50_raw = cv$PAM50[i],
  block_c = cv$complete_block_c[i])
TC$ER <- factor(dplyr::recode(TC$er_call, Positive = "ERpos",
                              Negative = "ERneg", .default = NA_character_),
                levels = c("ERpos", "ERneg"))
TC$PAM50 <- factor(sub("^BRCA_", "", TC$pam50_raw),
                   levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))
if (anyNA(TC$aliquot_tcga)) stop("aliquot barcode missing for some TCGA patients",
                                 call. = FALSE)
message("   ", nrow(TC), " TCGA patients; M_c (GISTIC call) has ",
        sum(is.na(TC$M_c)), " NA; PAM50 has ", sum(is.na(TC$PAM50)), " NA")

# --- 2.1 the six genes, and the SIGNED arm composites ------------------------
# SCALE: log2(linear DESeq2-normalised + 1), this repo's gene-level scale (E10,
# E16, E23). The composites are E16's `.comp` recipe - per-gene z across the
# FULL cohort (1,095), averaged - with each gene multiplied by its sign in the
# configuration first. The z-scores are fixed here; the bootstrap re-ranks the
# composites inside every resample and never re-standardises them, which is
# how E32 treated its GSVA scores.
lin <- readRDS(PATH_TCGA_LINEAR)
if (!identical(lin$scale, "linear_deseq2_normalised")) {
  stop("the linear matrix is not on the linear DESeq2-normalised scale",
       call. = FALSE)
}
GX <- log2(lin$mat[, ID_T, drop = FALSE] + 1)
rm(lin); invisible(gc(verbose = FALSE))
gr <- .gene_rows(names(CONFIG_SIGN), GX, .symbol_resolver(rownames(GX), NULL))
if (length(gr$missing)) {
  stop("did not resolve: ", paste(gr$missing, collapse = ", "), call. = FALSE)
}
GM <- gr$mat[names(CONFIG_SIGN), , drop = FALSE]
rm(GX)
gv <- apply(GM, 1L, stats::var)
if (any(!(gv > 0))) {
  stop("zero-variance gene(s): ", paste(names(gv)[!(gv > 0)], collapse = ", "),
       call. = FALSE)
}
ZG <- t(scale(t(GM)))
.signed_mean <- function(genes) {
  as.numeric(colMeans(sweep(ZG[genes, , drop = FALSE], 1L,
                            CONFIG_SIGN[genes], "*")))
}
TC$arm_trigger  <- .signed_mean(TRIGGER)
TC$arm_guardian <- .signed_mean(GUARDIAN)
TC$config_comp  <- .signed_mean(names(CONFIG_SIGN))
# The genes themselves are UNSIGNED: MCL1's own row therefore has the opposite
# sign to its contribution to the guardian arm, by construction.
for (g in names(CONFIG_SIGN)) TC[[paste0("g_", g)]] <- as.numeric(GM[g, ])
message("   arms: trigger = mean z(", paste(TRIGGER, collapse = ", "),
        "); guardian = mean(z BCL2L1, -z MCL1); composite = signed mean of 6")

# =============================================================================
# 3. The join - on the ALIQUOT
# =============================================================================
message("\n3. the join")

# E32's join, verbatim, kept to reproduce E32 exactly (section 4.1).
D_PAT <- FORK %>%
  dplyr::inner_join(TC, by = "patient")
if (nrow(D_PAT) != EXPECT_FORK_N) {
  stop("the patient-key join gives ", nrow(D_PAT), ", not ", EXPECT_FORK_N,
       call. = FALSE)
}
# THE ANALYSIS SET: the same RNA aliquot on both sides of every correlation.
al_ok <- D_PAT$aliquot == D_PAT$aliquot_tcga
MISMATCH <- D_PAT[!al_ok, c("patient", "aliquot", "aliquot_tcga", "block_c",
                           "in_849")]
D <- D_PAT[al_ok, , drop = FALSE]
n_bc_al  <- sum(D$block_c, na.rm = TRUE)
n_bc_pat <- sum(D_PAT$block_c, na.rm = TRUE)
message("   forkscale samples: ", nrow(FORK), "; patient join: ", nrow(D_PAT),
        "; ALIQUOT join: ", nrow(D), " (", round(100 * nrow(D) / nrow(FORK), 1),
        "%)")
message("   different vial on the two sides: ", nrow(MISMATCH), " - ",
        paste(MISMATCH$patient, collapse = ", "))
message("   Block C overlap: ", n_bc_al, " on the aliquot join (", n_bc_pat,
        " on the patient join; the declaration expected ~885)")
if (nrow(D) != EXPECT_ALIQUOT_N || n_bc_al != EXPECT_BLOCK_C_AL) {
  stop("the aliquot join or its Block C overlap differs from the data note (",
       EXPECT_ALIQUOT_N, " / ", EXPECT_BLOCK_C_AL, "). The inputs have moved.",
       call. = FALSE)
}

rho_forks <- stats::cor(D$MB1_forkscale, D$MB2_forkscale, method = "spearman")
message("   MB1 vs MB2 forkscale, Spearman: ", round(rho_forks, 3))

# =============================================================================
# 4. The bootstrap machinery - E32's estimator, and one resample set per frame
# =============================================================================
message("\n4. partial Spearman, ", format(N_BOOT, big.mark = ","),
        " bootstrap resamples")

# E32's estimator verbatim: rank, residualise on the ranked covariates,
# correlate the residuals. Written out so the bootstrap re-ranks INSIDE every
# resample rather than ranking once outside it.
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

# Companion only. Forkscale is severely skewed by construction, so a Pearson
# is dominated by a few low-index samples; it is never read for the verdict.
.partial_pearson <- function(x, y, Z = NULL) {
  ok <- stats::complete.cases(x, y, Z)
  if (sum(ok) < 10L) return(NA_real_)
  if (is.null(Z)) return(stats::cor(x[ok], y[ok]))
  H <- qr(cbind(1, as.matrix(Z)[ok, , drop = FALSE]))
  stats::cor(qr.resid(H, x[ok]), qr.resid(H, y[ok]))
}

.ci <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, c((1 - CI_LEVEL) / 2, 1 - (1 - CI_LEVEL) / 2),
                         na.rm = TRUE))
}

# One statistic, point estimate and bootstrap. Returns the row AND the
# bootstrap vector, so paired differences can be taken on the SAME resamples
# without refitting.
.run_block <- function(dat, boot_ix, x_col, y_col, z_cols, label, kind, limb,
                       set) {
  x <- dat[[x_col]]; y <- dat[[y_col]]
  Z <- if (length(z_cols)) as.matrix(dat[, z_cols, drop = FALSE]) else NULL
  est <- .partial_spearman(x, y, Z)
  bs  <- vapply(boot_ix, function(ix) {
    Zi <- if (is.null(Z)) NULL else Z[ix, , drop = FALSE]
    .partial_spearman(x[ix], y[ix], Zi)
  }, numeric(1))
  # EVERYTHING IS COMPUTED BEFORE THE TIBBLE. Inside tibble() a new column
  # shadows an outer variable of the same name (the trap that stopped E32's
  # first run), so nothing below is named x, y or Z.
  ci   <- .ci(bs)
  n_cc <- .n_complete(x, y, Z)
  r_pp <- .partial_pearson(x, y, Z)
  adj  <- if (length(z_cols)) paste(z_cols, collapse = " + ") else "-"
  excl <- isTRUE(ci[1] > 0) || isTRUE(ci[2] < 0)
  row <- tibble::tibble(
    set = set, limb = limb, kind = kind, label = label,
    x_var = x_col, y_var = y_col, adjusted_for = adj,
    n = n_cc, rho = est, ci_lo = ci[1], ci_hi = ci[2], ci_excludes_0 = excl,
    r_pearson_partial = r_pp)
  list(row = row, boot = bs)
}

.run_spec <- function(spec, dat, boot_ix, set) {
  labs <- vapply(spec, function(s) s$label, character(1))
  if (anyDuplicated(labs)) stop("duplicated labels in a spec", call. = FALSE)
  res <- lapply(spec, function(s)
    .run_block(dat, boot_ix, s$x, s$y, s$z, s$label, s$kind, s$limb, set))
  list(cor  = dplyr::bind_rows(lapply(res, function(r) r$row)),
       boot = stats::setNames(lapply(res, function(r) r$boot), labs))
}

# a - b on the SAME resamples: a genuine paired contrast, not two intervals
# compared by eye.
.paired <- function(res, lab_a, lab_b, label, limb) {
  ra <- res$cor[match(lab_a, res$cor$label), ]
  rb <- res$cor[match(lab_b, res$cor$label), ]
  if (is.na(ra$label) || is.na(rb$label)) {
    stop("paired difference asks for a label that was not run: ",
         lab_a, " / ", lab_b, call. = FALSE)
  }
  ci  <- .ci(res$boot[[lab_a]] - res$boot[[lab_b]])
  dif <- ra$rho - rb$rho
  set_name <- ra$set
  tibble::tibble(set = set_name, limb = limb, label = label,
                 minuend = lab_a, subtrahend = lab_b, diff = dif,
                 ci_lo = ci[1], ci_hi = ci[2],
                 ci_excludes_0 = isTRUE(ci[1] > 0) || isTRUE(ci[2] < 0))
}

SPEC <- list()
.add <- function(...) SPEC[[length(SPEC) + 1L]] <<- list(...)

# One resample index set per frame, shared by every statistic on that frame.
set.seed(PROJECT_SEED)
BOOT_IX  <- replicate(N_BOOT, sample.int(nrow(D), replace = TRUE),
                      simplify = FALSE)
# E32's draw exactly: same seed, same n (1,037), same row order, nothing
# drawn in between. Its intervals should therefore reproduce to the digit.
set.seed(PROJECT_SEED)
BOOT_PAT <- replicate(N_BOOT, sample.int(nrow(D_PAT), replace = TRUE),
                      simplify = FALSE)
D_BC <- D[which(D$block_c), , drop = FALSE]
set.seed(PROJECT_SEED + 1L)
BOOT_BC  <- replicate(N_BOOT, sample.int(nrow(D_BC), replace = TRUE),
                      simplify = FALSE)

SET_MAIN <- "aliquot join (PRIMARY)"
SET_BC   <- "Block C (sensitivity)"
SET_PAT  <- "patient join (E32 reproduction and sensitivity)"

# =============================================================================
# 5. LIMB 1 - descriptive. No verdict, no prediction.
# =============================================================================
message("\n5. LIMB 1 - descriptive")

# --- 5.1 the forkscale MEDIAN-SPLIT quadrant ---------------------------------
# E32's rule: strictly greater than the within-set median is "high". Medians of
# THIS set (1,031), so the cells differ slightly from E32's 1,037.
MED_MB1 <- stats::median(D$MB1_forkscale)
MED_MB2 <- stats::median(D$MB2_forkscale)
FSQ <- c("MB1 low / MB2 low", "MB1 high / MB2 low",
         "MB1 low / MB2 high", "MB1 high / MB2 high")
D$MB1_split <- factor(ifelse(D$MB1_forkscale > MED_MB1, "MB1 high", "MB1 low"),
                      levels = c("MB1 high", "MB1 low"))
D$MB2_split <- factor(ifelse(D$MB2_forkscale > MED_MB2, "MB2 high", "MB2 low"),
                      levels = c("MB2 high", "MB2 low"))
D$fs_quad <- factor(paste0(ifelse(D$MB1_forkscale > MED_MB1, "MB1 high",
                                  "MB1 low"), " / ",
                           ifelse(D$MB2_forkscale > MED_MB2, "MB2 high",
                                  "MB2 low")),
                    levels = FSQ)

# --- 5.2 MYC activity by forkscale quadrant ----------------------------------
# Median and IQR, every cell size stated. `med_pct` is the median patient's
# percentile within the FULL TCGA cohort (1,095) on that estimator - GSVA is
# cohort-relative, so "MYC-low" can only mean low against this cohort, and this
# is the number that says whether a cell is.
MYC_VARS <- c(M_a = "M_a", M_b = "M_b", M_c = "M_c",
              FELSHER__PROLIFSTRIP = "M_a_ps")
for (v in MYC_VARS) {
  ref <- TC[[v]][!is.na(TC[[v]])]
  F_v <- stats::ecdf(ref)
  val <- D[[v]]
  D[[paste0(v, "_pct")]] <- ifelse(is.na(val), NA_real_, 100 * F_v(val))
}

.myc_by <- function(dat, grp) {
  dplyr::bind_rows(lapply(names(MYC_VARS), function(est) {
    v  <- MYC_VARS[[est]]
    vp <- paste0(v, "_pct")
    dat %>%
      dplyr::filter(!is.na(.data[[grp]])) %>%
      dplyr::group_by(group = .data[[grp]]) %>%
      dplyr::summarise(
        estimator = est,
        n_group = dplyr::n(),
        n       = sum(!is.na(.data[[v]])),
        med     = stats::median(.data[[v]], na.rm = TRUE),
        q25     = unname(stats::quantile(.data[[v]], 0.25, na.rm = TRUE)),
        q75     = unname(stats::quantile(.data[[v]], 0.75, na.rm = TRUE)),
        med_pct = stats::median(.data[[vp]], na.rm = TRUE),
        .groups = "drop") %>%
      dplyr::mutate(grouping = grp)
  })) %>%
    dplyr::select(grouping, group, estimator, n_group, n, med, q25, q75,
                  med_pct)
}

MYC_BY_FSQ <- .myc_by(D, "fs_quad")
message("\n   MYC activity by forkscale MEDIAN-SPLIT quadrant (n = ", nrow(D), "):")
MYC_BY_FSQ %>%
  dplyr::mutate(dplyr::across(c(med, q25, q75), ~ round(.x, 3)),
                med_pct = round(med_pct, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Added in the data note (section 8.2), not in the declaration: the same by the
# STORED fork assignment, which is what "MB1_UF" and "MB2_UF" literally name.
# Kept in a separate table so the two groupings are never read as one.
MYC_BY_FORK <- dplyr::bind_rows(.myc_by(D, "MB1_fork"), .myc_by(D, "MB2_fork"))
message("\n   MYC activity by STORED fork assignment (NOT the median split):")
MYC_BY_FORK %>%
  dplyr::filter(estimator == "M_a") %>%
  dplyr::mutate(dplyr::across(c(med, q25, q75), ~ round(.x, 3)),
                med_pct = round(med_pct, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.3 the two partials of MYC on MB1, and the marginals -------------------
L1 <- "LIMB 1"
.add(x = "M_a", y = "MB1_forkscale", z = "PROLIF", limb = L1,
     label = "M_a ~ MB1 | PROLIF", kind = "limb 1 NEW (declared)")
.add(x = "M_a", y = "MB1_forkscale", z = "MB2_forkscale", limb = L1,
     label = "M_a ~ MB1 | MB2", kind = "limb 1, E32's quantity on this set")
.add(x = "M_a", y = "MB2_forkscale", z = "MB1_forkscale", limb = L1,
     label = "M_a ~ MB2 | MB1", kind = "limb 1, E32's quantity on this set")
.add(x = "M_a", y = "MB1_forkscale", z = character(0), limb = L1,
     label = "M_a ~ MB1 marginal", kind = "limb 1 marginal")
.add(x = "M_a", y = "MB2_forkscale", z = character(0), limb = L1,
     label = "M_a ~ MB2 marginal", kind = "limb 1 marginal")
# What bounds the reading: MB2 forkscale against proliferation itself.
.add(x = "PROLIF", y = "MB2_forkscale", z = character(0), limb = L1,
     label = "PROLIF ~ MB2 marginal", kind = "limb 1 bound")
.add(x = "PROLIF", y = "MB1_forkscale", z = character(0), limb = L1,
     label = "PROLIF ~ MB1 marginal", kind = "limb 1 bound")

# =============================================================================
# 6. LIMB 3 - the claim-relevant one. Specified here, run with limb 1.
# =============================================================================
L3 <- "LIMB 3"
.kind3 <- function(col) {
  if (col %in% names(READ_ARMS)) "LIMB 3 ARM (declared) - THE VERDICT"
  else if (col %in% names(READ_COMP)) "LIMB 3 composite (declared, third, not the verdict)"
  else "LIMB 3 per gene (declared, reported, not the verdict)"
}
for (col in names(READOUTS)) {
  nm <- READOUTS[[col]]
  k  <- .kind3(col)
  .add(x = col, y = "MB2_forkscale", z = "MB1_forkscale", limb = L3,
       label = paste0(nm, " ~ MB2 | MB1"), kind = k)
  .add(x = col, y = "MB1_forkscale", z = "MB2_forkscale", limb = L3,
       label = paste0(nm, " ~ MB1 | MB2"), kind = k)
  .add(x = col, y = "MB2_forkscale", z = character(0), limb = L3,
       label = paste0(nm, " ~ MB2 marginal"), kind = "LIMB 3 marginal (declared)")
  .add(x = col, y = "MB1_forkscale", z = character(0), limb = L3,
       label = paste0(nm, " ~ MB1 marginal"), kind = "LIMB 3 marginal (declared)")
}

# --- the proliferation companion. POST-HOC (data note 8.1). CANNOT CHANGE THE
# VERDICT. As in limb 1, which of MB2 and proliferation does the work is NOT
# interpreted - MB2 runs +0.874 with PROLIF_DISJOINT.
for (col in c(names(READ_ARMS), names(READ_COMP))) {
  nm <- READOUTS[[col]]
  .add(x = col, y = "MB2_forkscale", z = c("MB1_forkscale", "PROLIF"),
       limb = L3, label = paste0(nm, " ~ MB2 | MB1 + PROLIF"),
       kind = "LIMB 3 companion (POST-HOC)")
  .add(x = col, y = "MB1_forkscale", z = c("MB2_forkscale", "PROLIF"),
       limb = L3, label = paste0(nm, " ~ MB1 | MB2 + PROLIF"),
       kind = "LIMB 3 companion (POST-HOC)")
}

# --- DECLARED UNINFORMATIVE: OXPHOS against either fork ------------------------
# Computed and reported. NOT presented as support in either direction, whatever
# it returns, because both biclusters carry biogenesis. Declaration section 4.
for (r in OX_RULERS) {
  .add(x = r, y = "MB2_forkscale", z = "MB1_forkscale", limb = "UNINFORMATIVE",
       label = paste0(r, " ~ MB2 | MB1"), kind = "UNINFORMATIVE (declared)")
  .add(x = r, y = "MB1_forkscale", z = "MB2_forkscale", limb = "UNINFORMATIVE",
       label = paste0(r, " ~ MB1 | MB2"), kind = "UNINFORMATIVE (declared)")
  .add(x = r, y = "MB2_forkscale", z = character(0), limb = "UNINFORMATIVE",
       label = paste0(r, " ~ MB2 marginal"), kind = "UNINFORMATIVE (declared)")
  .add(x = r, y = "MB1_forkscale", z = character(0), limb = "UNINFORMATIVE",
       label = paste0(r, " ~ MB1 marginal"), kind = "UNINFORMATIVE (declared)")
}

t0 <- Sys.time()
RES_MAIN <- .run_spec(SPEC, D, BOOT_IX, SET_MAIN)
message("   ", nrow(RES_MAIN$cor), " statistics on the primary set in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# --- sensitivities: Block C, and the patient join ---------------------------
SPEC_ARMS <- list()
for (col in names(READ_ARMS)) {
  nm <- READOUTS[[col]]
  SPEC_ARMS[[length(SPEC_ARMS) + 1L]] <- list(
    x = col, y = "MB2_forkscale", z = "MB1_forkscale", limb = L3,
    label = paste0(nm, " ~ MB2 | MB1"), kind = "LIMB 3 ARM - sensitivity")
  SPEC_ARMS[[length(SPEC_ARMS) + 1L]] <- list(
    x = col, y = "MB1_forkscale", z = "MB2_forkscale", limb = L3,
    label = paste0(nm, " ~ MB1 | MB2"), kind = "LIMB 3 ARM - sensitivity")
}
RES_BC <- .run_spec(SPEC_ARMS, D_BC, BOOT_BC, SET_BC)

# E32's two declared partials, on E32's own frame and E32's own resamples.
SPEC_E32 <- list(
  list(x = "M_a", y = "MB2_forkscale", z = "MB1_forkscale", limb = "E32 CHECK",
       label = "M_a ~ MB2 | MB1", kind = "E32 reproduction"),
  list(x = "M_a", y = "MB1_forkscale", z = "MB2_forkscale", limb = "E32 CHECK",
       label = "M_a ~ MB1 | MB2", kind = "E32 reproduction"))
RES_PAT <- .run_spec(c(SPEC_E32, SPEC_ARMS), D_PAT, BOOT_PAT, SET_PAT)

COR <- dplyr::bind_rows(RES_MAIN$cor, RES_BC$cor, RES_PAT$cor)

# =============================================================================
# 4.1 Does E32 reproduce? Patient join, E32's seed, E32's resamples
# =============================================================================
message("\n   E32 REPRODUCTION, patient join, E32's seed and resamples:")
E32_QUOTED <- tibble::tibble(
  label = c("M_a ~ MB2 | MB1", "M_a ~ MB1 | MB2"),
  rho = c(0.688, -0.195), ci_lo = c(0.648, -0.247), ci_hi = c(0.723, -0.141))
e32_here <- RES_PAT$cor %>% dplyr::filter(limb == "E32 CHECK") %>%
  dplyr::select(label, n, rho, ci_lo, ci_hi)
e32_repro <- e32_here %>%
  dplyr::left_join(E32_QUOTED %>% dplyr::rename(rho_note = rho,
                                                lo_note = ci_lo,
                                                hi_note = ci_hi),
                   by = "label") %>%
  dplyr::mutate(matches_note_3dp = round(rho, 3) == rho_note &
                  round(ci_lo, 3) == lo_note & round(ci_hi, 3) == hi_note)
if (!is.null(e32)) {
  saved <- e32$correlations %>%
    dplyr::filter(kind == "PRIMARY (declared)",
                  label %in% E32_QUOTED$label) %>%
    dplyr::select(label, rho_saved = rho, lo_saved = ci_lo, hi_saved = ci_hi)
  e32_repro <- e32_repro %>% dplyr::left_join(saved, by = "label") %>%
    dplyr::mutate(max_abs_diff_vs_saved = pmax(abs(rho - rho_saved),
                                               abs(ci_lo - lo_saved),
                                               abs(ci_hi - hi_saved)))
}
e32_repro %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5.4 LIMB 1 results
# =============================================================================
.get <- function(lab, set = SET_MAIN) {
  r <- COR[COR$label == lab & COR$set == set, ]
  if (nrow(r) != 1L) stop("expected one row for ", lab, " in ", set, call. = FALSE)
  r
}
L1_DIFF <- .paired(RES_MAIN, "M_a ~ MB1 | PROLIF", "M_a ~ MB1 | MB2",
                   "MB1 | PROLIF minus MB1 | MB2", L1)
message("\n   LIMB 1 partials (n = ", nrow(D), "). DESCRIPTIVE - no verdict:")
COR %>% dplyr::filter(set == SET_MAIN, limb == L1) %>%
  dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   paired difference, MB1|PROLIF - MB1|MB2: ", round(L1_DIFF$diff, 3),
        " [", round(L1_DIFF$ci_lo, 3), ", ", round(L1_DIFF$ci_hi, 3), "]",
        "\n   REPORTED. WHICH ONE IS DOING THE WORK IS NOT INTERPRETED: MB2 ",
        "runs ", round(.get("PROLIF ~ MB2 marginal")$rho, 3),
        " with PROLIF_DISJOINT.")

# =============================================================================
# 6.1 LIMB 3 - THE VERDICT, on the rule fixed at 9559931
# =============================================================================
message("\n6. LIMB 3 - the verdict\n", strrep("-", 78))

.pair_rows <- function(nm, res, set) {
  p2 <- res$cor[res$cor$label == paste0(nm, " ~ MB2 | MB1"), ]
  p1 <- res$cor[res$cor$label == paste0(nm, " ~ MB1 | MB2"), ]
  pd <- .paired(res, paste0(nm, " ~ MB2 | MB1"), paste0(nm, " ~ MB1 | MB2"),
                paste0(nm, ": MB2|MB1 - MB1|MB2"), L3)
  sep_as_written <- isTRUE(p2$rho > 0) && isTRUE(p2$ci_lo > 0) &&
    isTRUE(p2$rho > p1$rho)
  sep_strict <- sep_as_written && isTRUE(pd$ci_lo > 0)
  tibble::tibble(
    set = set, readout = nm,
    rho_MB2_given_MB1 = p2$rho, lo_MB2_given_MB1 = p2$ci_lo,
    hi_MB2_given_MB1 = p2$ci_hi,
    rho_MB1_given_MB2 = p1$rho, lo_MB1_given_MB2 = p1$ci_lo,
    hi_MB1_given_MB2 = p1$ci_hi,
    paired_diff = pd$diff, diff_lo = pd$ci_lo, diff_hi = pd$ci_hi,
    separates = sep_as_written, separates_strict = sep_strict,
    n = p2$n)
}
.verdict <- function(sT, sG) {
  if (sT && !sG) V_TRIG
  else if (sT && sG) V_BOTH
  else if (!sT && !sG) V_NEITHER
  else V_OUTSIDE
}

SEP_MAIN <- dplyr::bind_rows(lapply(READOUTS, .pair_rows, res = RES_MAIN,
                                    set = SET_MAIN))
SEP_BC   <- dplyr::bind_rows(lapply(READ_ARMS, .pair_rows, res = RES_BC,
                                    set = SET_BC))
SEP_PAT  <- dplyr::bind_rows(lapply(READ_ARMS, .pair_rows, res = RES_PAT,
                                    set = SET_PAT))
PAIRED <- dplyr::bind_rows(
  L1_DIFF,
  dplyr::bind_rows(lapply(READOUTS, function(nm)
    .paired(RES_MAIN, paste0(nm, " ~ MB2 | MB1"), paste0(nm, " ~ MB1 | MB2"),
            paste0(nm, ": MB2|MB1 - MB1|MB2"), L3))),
  dplyr::bind_rows(lapply(READ_ARMS, function(nm)
    .paired(RES_BC, paste0(nm, " ~ MB2 | MB1"), paste0(nm, " ~ MB1 | MB2"),
            paste0(nm, ": MB2|MB1 - MB1|MB2"), L3))),
  dplyr::bind_rows(lapply(READ_ARMS, function(nm)
    .paired(RES_PAT, paste0(nm, " ~ MB2 | MB1"), paste0(nm, " ~ MB1 | MB2"),
            paste0(nm, ": MB2|MB1 - MB1|MB2"), L3))),
  .paired(RES_PAT, "M_a ~ MB2 | MB1", "M_a ~ MB1 | MB2",
          "E32: M_a MB2|MB1 - MB1|MB2", "E32 CHECK"))

.sep_of <- function(tab, nm, col = "separates") tab[[col]][tab$readout == nm]
sT <- .sep_of(SEP_MAIN, "trigger arm"); sG <- .sep_of(SEP_MAIN, "guardian arm")
verdict <- .verdict(sT, sG)
verdict_strict <- .verdict(.sep_of(SEP_MAIN, "trigger arm", "separates_strict"),
                           .sep_of(SEP_MAIN, "guardian arm", "separates_strict"))
verdict_bc  <- .verdict(.sep_of(SEP_BC, "trigger arm"),
                        .sep_of(SEP_BC, "guardian arm"))
verdict_pat <- .verdict(.sep_of(SEP_PAT, "trigger arm"),
                        .sep_of(SEP_PAT, "guardian arm"))

.fmt_row <- function(r) paste0(
  sprintf("%-24s", r$readout),
  " MB2|MB1 ", sprintf("%+.3f [%+.3f, %+.3f]", r$rho_MB2_given_MB1,
                       r$lo_MB2_given_MB1, r$hi_MB2_given_MB1),
  "   MB1|MB2 ", sprintf("%+.3f [%+.3f, %+.3f]", r$rho_MB1_given_MB2,
                         r$lo_MB1_given_MB2, r$hi_MB1_given_MB2),
  "   diff ", sprintf("%+.3f [%+.3f, %+.3f]", r$paired_diff, r$diff_lo,
                      r$diff_hi),
  "   separates: ", r$separates)
message("RULE: ", VERDICT_RULE, "\n")
for (k in seq_len(nrow(SEP_MAIN))) {
  if (SEP_MAIN$readout[k] == READ_COMP[[1]]) message("   -- third, not the verdict:")
  if (SEP_MAIN$readout[k] == READ_GENES[[1]]) {
    message("   -- per gene, not the verdict (MCL1 UNSIGNED):")
  }
  message("   ", .fmt_row(SEP_MAIN[k, ]))
}
message("\n   VERDICT: ", verdict)
message("   the stricter reading (paired difference also excludes 0), NOT the ",
        "verdict: ", verdict_strict)
message("   sensitivities: Block C (n = ", nrow(D_BC), ") ", verdict_bc,
        "; patient join (n = ", nrow(D_PAT), ") ", verdict_pat)

reading <- switch(
  verdict,
  "TRIGGER MB2-SPECIFIC" = paste0(
    "TRIGGER MB2-SPECIFIC. The trigger arm tracks MB2 at fixed MB1 and not\n",
    "   the reverse, and the guardian arm does not. On the declaration this\n",
    "   licenses ONE Discussion clause: that the standing pro-apoptotic\n",
    "   transcript signal is a property of tumours whose respiration is\n",
    "   MYC-coupled. It licenses NOTHING about response to an OXPHOS-directed\n",
    "   or BH3-directed agent. Read the proliferation companion first."),
  "BOTH ARMS SEPARATE" = paste0(
    "BOTH ARMS SEPARATE. The whole configuration is MB2-associated. This is\n",
    "   descriptive and does NOT isolate the trigger, so the Discussion\n",
    "   premise is not sharpened by it."),
  "NEITHER SEPARATES" = paste0(
    "NEITHER SEPARATES. Neither arm tracks MB2 at fixed MB1 in the declared\n",
    "   sense: the configuration is not programme-specific. Reportable, NOT a\n",
    "   failure, and consistent with E11's finding that it is an OXPHOS\n",
    "   correlate. The manuscript should stop implying a programme-specificity\n",
    "   the configuration does not have."),
  paste0(
    "OUTSIDE THE DECLARED OUTCOMES. The guardian arm separates and the\n",
    "   trigger arm does not. The declaration did not name this outcome, so it\n",
    "   is reported literally and given NO licensing reading."))
message("\n", reading)

message("\n   the proliferation companion (POST-HOC - cannot change the verdict):")
COR %>% dplyr::filter(set == SET_MAIN,
                      kind == "LIMB 3 companion (POST-HOC)") %>%
  dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3), ci_excludes_0) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   DECLARED UNINFORMATIVE - both biclusters carry biogenesis. ",
        "Reported;\n   NOT support in either direction whatever it says.")
COR %>% dplyr::filter(set == SET_MAIN, kind == "UNINFORMATIVE (declared)") %>%
  dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                   ci_hi = round(ci_hi, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. LIMB 2 - the variable map. Counts only.
# =============================================================================
message("\n7. LIMB 2 - the variable map")

# --- 7.1 STATE, read-only, and its integrity contract ------------------------
st <- readRDS(PATH_STATE)
if (!identical(deparse(st$build_state), st$definition_source$build_state) ||
    !identical(environment(st$build_state), baseenv())) {
  stop("the frozen STATE constructor fails its own integrity contract",
       call. = FALSE)
}
.build_state <- st$build_state
S  <- st$tcga$state
LV <- levels(S$STATE_gsva)
if (!identical(LV, st$spec$levels)) stop("STATE levels moved", call. = FALSE)
# The frozen MYC IS this repo's M_a (data note section 4). Checked, not assumed.
rho_myc_frozen <- stats::cor(S$MYC, as.numeric(nw$tcga_gsva_new[MYC_REF,
                                                                 S$patient]),
                             method = "spearman")
message("   frozen MYC vs this repo's M_a, Spearman: ", round(rho_myc_frozen, 6))

# --- 7.2 the MYC x OXPHOS quadrant, from the constructor's OWN calls ---------
# The constructor computes its OXPHOS call internally and does not return it.
# Calling it a second time with myc := oxphos makes its MYC-high call BE its
# OXPHOS-high call: same median, same `>` tie rule, same complete-case set.
# No threshold is transcribed into this repo. Level 2 cannot occur in that
# probe (MYC-high but OXPHOS-low of the same vector), and it is asserted not to.
QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",
          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")
INSTR <- c(gsva = "OX_gsva", mitopps = "OX_mitopps")
myc_high <- S$STATE_gsva != LV[1]           # the constructor's MYC call
# The MYC call does not depend on the OXPHOS instrument; checked, not assumed.
if (!identical(myc_high, S$STATE_mitopps != LV[1])) {
  stop("STATE_gsva and STATE_mitopps disagree on the MYC call", call. = FALSE)
}
QUAD <- tibble::tibble(patient = S$patient,
                       STATE_gsva = S$STATE_gsva,
                       STATE_mitopps = S$STATE_mitopps,
                       myc_high_frozen = myc_high)
ASSERTS <- list()
for (inst in names(INSTR)) {
  ox <- S[[INSTR[[inst]]]]
  state_stored <- S[[paste0("STATE_", inst)]]
  rebuilt <- .build_state(S$MYC, ox, S$BUFFER_gistic)
  if (!identical(rebuilt, state_stored)) {
    stop("the constructor does not reproduce STATE_", inst, " from MYC, ",
         INSTR[[inst]], " and BUFFER_gistic", call. = FALSE)
  }
  probe <- .build_state(ifelse(is.na(S$MYC), NA_real_, ox), ox, S$BUFFER_gistic)
  if (any(probe == LV[2], na.rm = TRUE)) {
    stop("level 2 appeared in the myc := oxphos probe; the call is not what ",
         "it is taken to be", call. = FALSE)
  }
  if (!identical(is.na(probe), is.na(state_stored))) {
    stop("the probe's complete-case set differs from STATE's", call. = FALSE)
  }
  ox_high <- probe != LV[1]
  q <- dplyr::case_when(
    !myc_high & !ox_high ~ QLEV[1],
    !myc_high &  ox_high ~ QLEV[2],
     myc_high & !ox_high ~ QLEV[3],
     myc_high &  ox_high ~ QLEV[4])
  # THE DECLARED ASSERTIONS. Q3 and Q4 are built from the MYC call and the
  # EXTRACTED OXPHOS call, independently of STATE's levels 2-4, so these test
  # the extraction against the call the constructor made internally.
  a3 <- identical(which(q == QLEV[3]), which(state_stored == LV[2]))
  a4 <- identical(which(q == QLEV[4]), which(state_stored %in% LV[3:4]))
  if (!a3 || !a4) {
    stop("ASSERTION FAILED (", inst, "): Q3 == level 2 is ", a3,
         ", Q4 == levels 3 + 4 is ", a4, ". Something is wrong with the ",
         "constructor call. STOP.", call. = FALSE)
  }
  ASSERTS[[inst]] <- tibble::tibble(instrument = inst,
                                    Q3_equals_level2 = a3,
                                    Q4_equals_levels3_4 = a4)
  QUAD[[paste0("ox_high_", inst)]] <- ox_high
  QUAD[[paste0("quad_", inst)]]    <- factor(q, levels = QLEV)
  message("   ", inst, ": Q3 == level 2 TRUE; Q4 == levels 3 + 4 TRUE")
}
ASSERTS <- dplyr::bind_rows(ASSERTS)

QUAD_FULL <- QUAD %>% dplyr::count(quad_gsva, name = "n_gsva") %>%
  dplyr::full_join(QUAD %>% dplyr::count(quad_mitopps, name = "n_mitopps"),
                   by = c("quad_gsva" = "quad_mitopps")) %>%
  dplyr::rename(quadrant = quad_gsva)
message("\n   the quadrant, FULL cohort (", nrow(QUAD), "; cuts are the ",
        "constructor's, on its ", sum(!is.na(QUAD$quad_gsva)),
        " complete cases):")
QUAD_FULL %>% as.data.frame() %>% print(row.names = FALSE)
QUAD_INSTR <- table(gsva = QUAD$quad_gsva, mitopps = QUAD$quad_mitopps)
message("   GSVA and mitoPPS quadrant agree on ",
        round(100 * sum(diag(QUAD_INSTR)) / sum(QUAD_INSTR), 1), "%")

# --- 7.3 the sanity check: the MYC call on FELSHER__PROLIFSTRIP --------------
# One line, and NOTHING IS GATED ON IT. The same constructor, the same OXPHOS
# and BUFFER inputs, M_a replaced by the proliferation-stripped estimator.
ps <- as.numeric(nw$tcga_gsva_new[MYC_PROLIFSTRIP, S$patient])
myc_high_ps <- .build_state(ps, S$OX_gsva, S$BUFFER_gistic) != LV[1]
SANITY_TAB <- table(frozen_M_a = myc_high, PROLIFSTRIP = myc_high_ps)
.kappa <- function(tb) {
  n <- sum(tb); po <- sum(diag(tb)) / n
  pe <- sum(rowSums(tb) * colSums(tb)) / n^2
  (po - pe) / (1 - pe)
}
SANITY <- tibble::tibble(
  n = sum(SANITY_TAB),
  agreement = sum(diag(SANITY_TAB)) / sum(SANITY_TAB),
  kappa = .kappa(SANITY_TAB),
  frozen_high_ps_low = SANITY_TAB["TRUE", "FALSE"],
  frozen_low_ps_high = SANITY_TAB["FALSE", "TRUE"])
message("\n   SANITY CHECK (nothing gated): the MYC call on ", MYC_PROLIFSTRIP,
        " agrees with the frozen M_a call on ", round(100 * SANITY$agreement, 1),
        "% of ", SANITY$n, " (kappa ", round(SANITY$kappa, 3), "; ",
        SANITY$frozen_high_ps_low, " high -> low, ", SANITY$frozen_low_ps_high,
        " low -> high)")

# --- 7.4 the cross-tabulations, on the analysis set --------------------------
D <- D %>% dplyr::left_join(QUAD, by = "patient")
if (any(!D$patient %in% QUAD$patient)) stop("a joined patient has no STATE row",
                                             call. = FALSE)

# Every level shown, NA as an explicit level, row and column totals.
.lv <- function(v, na_lab) {
  f  <- if (is.factor(v)) v else factor(v)
  lv <- levels(f); ch <- as.character(f)
  if (anyNA(ch)) { ch[is.na(ch)] <- na_lab; lv <- c(lv, na_lab) }
  factor(ch, levels = lv)
}
.xtab <- function(dat, row, col, na_row = "no STATE (BUFFER NA)",
                  na_col = "NA") {
  rv <- .lv(dat[[row]], na_row); cvv <- .lv(dat[[col]], na_col)
  m  <- as.data.frame.matrix(table(rv, cvv))
  out <- tibble::as_tibble(m, rownames = row)
  out$total <- as.integer(rowSums(m))
  col_tot <- lapply(as.list(colSums(m)), as.integer)     # keeps the names
  tot <- tibble::as_tibble(c(stats::setNames(list("total"), row), col_tot,
                             list(total = as.integer(sum(m)))),
                           .name_repair = "minimal")
  if (!identical(names(tot), names(out))) stop("xtab totals misaligned",
                                               call. = FALSE)
  dplyr::bind_rows(out, tot)
}

XT <- list()
for (inst in names(INSTR)) {
  qv <- paste0("quad_", inst)
  XT[[paste0(inst, ": quadrant x STATE")]]    <- .xtab(D, qv, paste0("STATE_", inst))
  XT[[paste0(inst, ": quadrant x MB1 fork")]] <- .xtab(D, qv, "MB1_fork")
  XT[[paste0(inst, ": quadrant x MB2 fork")]] <- .xtab(D, qv, "MB2_fork")
  XT[[paste0(inst, ": quadrant x PAM50")]]    <- .xtab(D, qv, "PAM50")
  XT[[paste0(inst, ": quadrant x ER")]]       <- .xtab(D, qv, "ER")
}
# Forkscale median splits beside the fork assignments - two different
# variables, shown against each other so the difference is visible.
XT[["gsva: quadrant x forkscale median-split quadrant"]] <- .xtab(D, "quad_gsva", "fs_quad")
XT[["MB1 fork x MB1 forkscale median split"]] <- .xtab(D, "MB1_fork", "MB1_split", na_row = "NA")
XT[["MB2 fork x MB2 forkscale median split"]] <- .xtab(D, "MB2_fork", "MB2_split", na_row = "NA")
XT[["MB1 fork x MB2 fork"]] <- .xtab(D, "MB1_fork", "MB2_fork", na_row = "NA")
# Sensitivities and context.
XT[["gsva: quadrant x MB1 fork (bounded extension)"]] <-
  .xtab(D, "quad_gsva", "MB1_fork_ext", na_col = "undetermined")
XT[["gsva: quadrant x MB2 fork (bounded extension)"]] <-
  .xtab(D, "quad_gsva", "MB2_fork_ext", na_col = "undetermined")
XT[["gsva: quadrant x PAM50 (Menegollo's call)"]] <- .xtab(D, "quad_gsva", "PAM50_menegollo")
# The UF/LF labelling against Menegollo's descriptions (MB1_UF almost all ER+,
# MB2_UF ER-negative). CONTEXT, NOT A TEST.
XT[["MB1 fork x ER"]] <- .xtab(D, "MB1_fork", "ER", na_row = "NA")
XT[["MB2 fork x ER"]] <- .xtab(D, "MB2_fork", "ER", na_row = "NA")
XT[["MB1 fork x PAM50"]] <- .xtab(D, "MB1_fork", "PAM50", na_row = "NA")
XT[["MB2 fork x PAM50"]] <- .xtab(D, "MB2_fork", "PAM50", na_row = "NA")

for (nm in c("gsva: quadrant x STATE", "gsva: quadrant x MB1 fork",
             "gsva: quadrant x MB2 fork", "gsva: quadrant x PAM50",
             "gsva: quadrant x ER")) {
  message("\n   ", nm, " (n = ", nrow(D), ")")
  XT[[nm]] %>% as.data.frame() %>% print(row.names = FALSE)
}

pm_ok <- !is.na(D$PAM50) & !is.na(D$PAM50_menegollo)
PAM50_AGREE <- tibble::tibble(
  n_both = sum(pm_ok),
  agreement = mean(as.character(D$PAM50[pm_ok]) == D$PAM50_menegollo[pm_ok]),
  repo_na = sum(is.na(D$PAM50)))
message("\n   PAM50: this repo's call agrees with Menegollo's on ",
        round(100 * PAM50_AGREE$agreement, 1), "% of ", PAM50_AGREE$n_both)

# =============================================================================
# 8. Figures
# =============================================================================
message("\n8. figures")

# --- Figure 1: the forkscale plane, with MYC activity overlaid ---------------
# Rank-rank, because every statistic here is rank-based (E32's reasoning). The
# dashed lines are the MEDIAN SPLITS, not the UF/LF assignments.
PD1 <- D %>% dplyr::mutate(MB1_rank = rank(MB1_forkscale),
                           MB2_rank = rank(MB2_forkscale))
cut1 <- sum(D$MB1_forkscale <= MED_MB1) + 0.5
cut2 <- sum(D$MB2_forkscale <= MED_MB2) + 0.5
nR <- nrow(PD1)
ann <- MYC_BY_FSQ %>% dplyr::filter(estimator == "M_a") %>%
  dplyr::mutate(
    x = ifelse(grepl("^MB1 high", group), nR * 0.97, nR * 0.03),
    y = ifelse(grepl("MB2 high$", group), nR * 0.97, nR * 0.03),
    hj = ifelse(grepl("^MB1 high", group), 1, 0),
    vj = ifelse(grepl("MB2 high$", group), 1, 0),
    txt = paste0(group, "\nn = ", n_group, "\nmedian M_a at pct ",
                 round(med_pct)))

p1a <- ggplot2::ggplot(PD1, ggplot2::aes(MB1_rank, MB2_rank)) +
  ggplot2::geom_point(ggplot2::aes(colour = M_a_pct), size = 1.1, alpha = 0.85) +
  ggplot2::geom_vline(xintercept = cut1, linetype = "dashed", colour = "grey30") +
  ggplot2::geom_hline(yintercept = cut2, linetype = "dashed", colour = "grey30") +
  ggplot2::geom_label(data = ann, ggplot2::aes(x = x, y = y, label = txt,
                                                hjust = hj, vjust = vj),
                      size = 2.3, label.size = 0.2, alpha = 0.9,
                      inherit.aes = FALSE) +
  ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#E6E6E6",
                                  high = "#B2182B", midpoint = 50,
                                  limits = c(0, 100),
                                  name = "M_a, TCGA\npercentile") +
  ggplot2::coord_equal() +
  ggplot2::labs(subtitle = paste0("A. MYC activity (M_a) on the forkscale plane\n",
                                  "dashed = within-set MEDIAN SPLITS, not UF/LF"),
                x = "MB1 forkscale, rank", y = "MB2 forkscale, rank") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom")

QCOL <- c("Q1 MYC-low OXPHOS-low" = "#56B4E9", "Q2 MYC-low OXPHOS-high" = "#009E73",
          "Q3 MYC-high OXPHOS-low" = "#E69F00", "Q4 MYC-high OXPHOS-high" = "#D55E00",
          "no STATE (BUFFER NA)" = "#BDBDBD")
QSHP <- c("Q1 MYC-low OXPHOS-low" = 16, "Q2 MYC-low OXPHOS-high" = 17,
          "Q3 MYC-high OXPHOS-low" = 15, "Q4 MYC-high OXPHOS-high" = 18,
          "no STATE (BUFFER NA)" = 4)
PD1$quad_lab <- .lv(PD1$quad_gsva, "no STATE (BUFFER NA)")
p1b <- ggplot2::ggplot(PD1, ggplot2::aes(MB1_rank, MB2_rank, colour = quad_lab,
                                         shape = quad_lab)) +
  ggplot2::geom_point(size = 1.2, alpha = 0.85) +
  ggplot2::geom_vline(xintercept = cut1, linetype = "dashed", colour = "grey30") +
  ggplot2::geom_hline(yintercept = cut2, linetype = "dashed", colour = "grey30") +
  ggplot2::scale_colour_manual(values = QCOL, name = NULL, drop = TRUE) +
  ggplot2::scale_shape_manual(values = QSHP, name = NULL, drop = TRUE) +
  ggplot2::coord_equal() +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2),
                  shape = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::labs(subtitle = paste0("B. the MYC x OXPHOS quadrant (frozen ",
                                  "constructor, GSVA)\non the same plane"),
                x = "MB1 forkscale, rank", y = "MB2 forkscale, rank") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom")

f1 <- patchwork::wrap_plots(p1a, p1b, nrow = 1) +
  patchwork::plot_annotation(
    title = paste0("MB1 and MB2 forkscale, MYC activity and the MYC x OXPHOS ",
                   "quadrant - ", nR, " TCGA-BRCA tumours (aliquot join)"),
    subtitle = paste0(
      "LIMBS 1 AND 2 - DESCRIPTIVE, NO VERDICT. rho(M_a, MB1 | PROLIF) = ",
      sprintf("%+.3f", .get("M_a ~ MB1 | PROLIF")$rho), "; rho(M_a, MB1 | MB2) = ",
      sprintf("%+.3f", .get("M_a ~ MB1 | MB2")$rho), "; marginals MB1 ",
      sprintf("%+.3f", .get("M_a ~ MB1 marginal")$rho), ", MB2 ",
      sprintf("%+.3f", .get("M_a ~ MB2 marginal")$rho),
      ". A gradient, not two states.",
      "\nEXPLORATORY. Not pre-registered."),
    caption = paste0(
      "Forkscale = pc1 / index from the pinned Menegollo snapshot (upstream ",
      "8fdbb34), E32's construction. M_a percentile is within the full TCGA ",
      "cohort (1,095); GSVA is cohort-relative, so 'low' means low against this ",
      "cohort only.\nMB2 forkscale runs ",
      sprintf("%+.3f", .get("PROLIF ~ MB2 marginal")$rho),
      " with PROLIF_DISJOINT, so conditioning on MB2 and on proliferation are ",
      "close to indistinguishable here. STATE is read read-only from the frozen ",
      "validation repo. No outcome variable enters."),
    theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, hjust = 0),
                           plot.title.position = "plot"))

ggplot2::ggsave(PATH_E33_FIG1, f1, width = 11, height = 6.8, dpi = 200)
.ensure_dir(DIR_DOC_FIGURES)
ggplot2::ggsave(PATH_E33_FIG1_DOC, f1, width = 11, height = 6.8, dpi = 200)
message("   ", PATH_E33_FIG1, "\n   ", PATH_E33_FIG1_DOC, "  (tracked)")

# --- Figure 2: limb 3, paired partials, trigger beside guardian --------------
# Every readout carries BOTH partials, always - reporting only the one that
# fires is the failure mode this design exists to prevent. Marginals hollow.
GROUP_OF <- c(stats::setNames(rep("Trigger arm", 5),
                              c("trigger arm", TRIGGER)),
              stats::setNames(rep("Guardian arm", 3),
                              c("guardian arm", GUARDIAN)),
              stats::setNames("Configuration\n(third, not the verdict)",
                              "configuration composite"))
ROW_ORDER <- c("trigger arm", TRIGGER, "guardian arm", GUARDIAN,
               "configuration composite")
ROW_LAB <- c("trigger arm" = "TRIGGER ARM (signed composite)",
             BBC3 = "BBC3", BID = "BID", BIK = "BIK", BAD = "BAD",
             "guardian arm" = "GUARDIAN ARM (BCL2L1 - MCL1, in z)",
             BCL2L1 = "BCL2L1", MCL1 = "MCL1 (unsigned; enters the arm as -)",
             "configuration composite" = "configuration composite (signed, 6 genes)")
Q2LAB <- c("MB2 | MB1", "MB1 | MB2", "MB2 marginal", "MB1 marginal")
PD2 <- dplyr::bind_rows(lapply(ROW_ORDER, function(nm) {
  labs <- paste0(nm, " ~ ", Q2LAB)
  r <- COR[COR$set == SET_MAIN & COR$label %in% labs, ]
  r$quantity <- factor(sub(paste0(nm, " ~ "), "", r$label, fixed = TRUE),
                       levels = Q2LAB)
  r$readout <- nm
  r
})) %>%
  dplyr::mutate(
    row_i = match(readout, ROW_ORDER),
    offset = c(`MB2 | MB1` = 0.18, `MB1 | MB2` = -0.18,
               `MB2 marginal` = 0.18, `MB1 marginal` = -0.18)[as.character(quantity)],
    y = -(row_i) + offset,
    group = factor(GROUP_OF[readout], levels = unique(GROUP_OF)),
    partial = quantity %in% c("MB2 | MB1", "MB1 | MB2"))
PD2 <- PD2 %>% dplyr::mutate(
  pair = factor(ifelse(grepl("^MB2", as.character(quantity)), "MB2", "MB1"),
                levels = c("MB2", "MB1")))

y_breaks <- -seq_along(ROW_ORDER)
y_labels <- unname(ROW_LAB[ROW_ORDER])

PCOL <- c(MB2 = "#D55E00", MB1 = "#0072B2")
f2 <- ggplot2::ggplot(PD2, ggplot2::aes(x = rho, y = y, colour = pair)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  ggplot2::geom_errorbar(data = dplyr::filter(PD2, partial),
                         ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                         orientation = "y", width = 0.12, linewidth = 0.6) +
  ggplot2::geom_point(data = dplyr::filter(PD2, partial), size = 2.4,
                      shape = 16) +
  ggplot2::geom_point(data = dplyr::filter(PD2, !partial), size = 2.4,
                      shape = 1, stroke = 0.8) +
  ggplot2::scale_colour_manual(
    values = PCOL, name = NULL,
    labels = c(MB2 = "rho(readout, MB2 | MB1)  - filled; hollow = MB2 marginal",
               MB1 = "rho(readout, MB1 | MB2)  - filled; hollow = MB1 marginal")) +
  ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels,
                              expand = ggplot2::expansion(add = 0.5)) +
  ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  ggplot2::labs(x = "partial Spearman rho, 95% bootstrap CI (10,000 resamples)",
                y = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_blank(),
                 legend.position = "bottom", legend.direction = "vertical",
                 strip.text.y = ggplot2::element_text(angle = 0),
                 plot.title.position = "plot",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
.arm_line <- function(nm) {
  r <- SEP_MAIN[SEP_MAIN$readout == nm, ]
  sprintf("%s: MB2|MB1 %+.3f [%+.3f, %+.3f] vs MB1|MB2 %+.3f [%+.3f, %+.3f]",
          nm, r$rho_MB2_given_MB1, r$lo_MB2_given_MB1, r$hi_MB2_given_MB1,
          r$rho_MB1_given_MB2, r$lo_MB1_given_MB2, r$hi_MB1_given_MB2)
}
f2 <- f2 + ggplot2::labs(
  title = paste0("Trigger arm against guardian arm on MB2 and MB1 forkscale - ",
                 nR, " TCGA-BRCA tumours"),
  subtitle = paste0("VERDICT on the rule fixed at 9559931 before any data was ",
                    "read: ", verdict, "\n", .arm_line("trigger arm"), "\n",
                    .arm_line("guardian arm"),
                    "\nEXPLORATORY. Not pre-registered. Both partials shown for ",
                    "every readout."),
  caption = paste0(
    "Arms: E16's mean-z recipe on log2(linear + 1), each gene signed by its ",
    "direction in the configuration (MCL1 -1). Genes shown UNSIGNED. The verdict ",
    "turns on the two ARM rows only.\nNo outcome variable enters. NOTHING HERE ",
    "LICENSES A TREATMENT-SELECTION CLAIM: no OXPHOS-directed or BH3-directed ",
    "intervention exists in any cohort available."))

ggplot2::ggsave(PATH_E33_FIG2, f2, width = 8.5, height = 7.2, dpi = 200)
ggplot2::ggsave(PATH_E33_FIG2_DOC, f2, width = 8.5, height = 7.2, dpi = 200)
message("   ", PATH_E33_FIG2, "\n   ", PATH_E33_FIG2_DOC, "  (tracked)")

# =============================================================================
# 9. Save
# =============================================================================
message("\n9. saving")

out <- list(
  verdict        = verdict,
  reading        = reading,
  verdict_rule   = VERDICT_RULE,
  verdict_inputs = SEP_MAIN %>% dplyr::filter(readout %in% READ_ARMS),
  verdict_strict = verdict_strict,
  verdict_sensitivity = tibble::tibble(
    set = c(SET_BC, SET_PAT), n = c(nrow(D_BC), nrow(D_PAT)),
    verdict = c(verdict_bc, verdict_pat)),
  correlations   = COR,
  paired_differences = PAIRED,
  limb1 = list(
    myc_by_forkscale_quadrant = MYC_BY_FSQ,
    myc_by_fork_assignment    = MYC_BY_FORK,
    partials   = COR %>% dplyr::filter(set == SET_MAIN, limb == L1),
    difference = L1_DIFF,
    medians    = c(MB1 = MED_MB1, MB2 = MED_MB2)),
  limb2 = list(
    crosstabs = XT,
    quadrant_full_cohort = QUAD_FULL,
    quadrant_instrument_table = QUAD_INSTR,
    asserts   = ASSERTS,
    sanity_prolifstrip = list(table = SANITY_TAB, summary = SANITY),
    pam50_agreement = PAM50_AGREE,
    rho_frozen_myc_vs_M_a = rho_myc_frozen),
  limb3 = list(
    separation = dplyr::bind_rows(SEP_MAIN, SEP_BC, SEP_PAT),
    companion  = COR %>% dplyr::filter(kind == "LIMB 3 companion (POST-HOC)"),
    uninformative = COR %>% dplyr::filter(kind == "UNINFORMATIVE (declared)")),
  e32_reproduction = e32_repro,
  forkscale_vs_e32 = fs_vs_e32,
  assembled_crosscheck = xchk,
  fork_assignments = FORK_BOUNDS,
  join = list(
    fork_samples = nrow(FORK), patient_join = nrow(D_PAT),
    aliquot_join = nrow(D), match_rate = nrow(D) / nrow(FORK),
    different_vial = MISMATCH,
    block_c_aliquot = n_bc_al, block_c_patient = n_bc_pat),
  rho_mb1_mb2 = rho_forks,
  analysis_frame = D %>% dplyr::select(
    patient, aliquot, in_849, MB1_forkscale, MB2_forkscale, MB1_split,
    MB2_split, fs_quad, MB1_fork, MB2_fork, MB1_fork_ext, MB2_fork_ext,
    M_a, M_b, M_c, M_a_ps, PROLIF, ox_gsva, ox_ppd, arm_trigger, arm_guardian,
    config_comp, dplyr::starts_with("g_"), STATE_gsva, STATE_mitopps,
    quad_gsva, quad_mitopps, ER, PAM50, PAM50_menegollo, block_c),
  counts = list(
    analysis_set = nrow(D), block_c = nrow(D_BC),
    m_c_na = sum(is.na(D$M_c)), er_na = sum(is.na(D$ER)),
    pam50_na = sum(is.na(D$PAM50)), state_na = sum(is.na(D$quad_gsva))),
  spec = list(
    declaration = "docs/2026-09-29_E33_declaration_v2.md, committed ALONE at 9559931",
    data_note   = "docs/2026-09-30_e33_data.md, 9b3ae75",
    posture = paste("EXPLORATORY, not pre-registered. The verdict rule was",
                    "fixed BEFORE any data file was opened."),
    arms = paste("E16 .comp recipe on log2(linear + 1), each gene signed by",
                 "its direction in the configuration: BBC3, BID, BIK, BAD,",
                 "BCL2L1 +1, MCL1 -1. Chosen by the author 2026-09-30, before",
                 "the script. Genes are reported UNSIGNED."),
    exceeds = paste("'Exceeds' is the point estimate, as written. The paired",
                    "difference is the stricter reading and is NOT the verdict."),
    analysis_set = paste("ALIQUOT join, n = 1,031; six patients were profiled",
                         "from different vials and are excluded. The patient",
                         "join (1,037) reproduces E32 and is a sensitivity."),
    fork = paste("UF/LF/none = the STORED calls from the assembled frame (849);",
                 "the other 188 are 'not in 849'. The bounded extension is a",
                 "sensitivity. NOT the forkscale.fork_0.4/0.8 columns."),
    split_vs_fork = paste("forkscale MEDIAN SPLITS (fs_quad, MB{k}_split) and",
                          "UF/LF ASSIGNMENTS (MB{k}_fork) are different",
                          "variables and are never merged."),
    quadrant = paste("the frozen constructor called twice: once to reproduce",
                     "STATE, once with myc := oxphos to extract its OXPHOS",
                     "call. Q3 == level 2 and Q4 == levels 3 + 4 asserted."),
    uninformative = paste("OXPHOS vs either forkscale is DECLARED",
                          "UNINFORMATIVE - both biclusters carry biogenesis."),
    companion = paste("The limb 3 proliferation companion is POST-HOC",
                      "(data note 8.1) and CANNOT change the verdict."),
    limb1_bound = paste("MB2 forkscale runs ~+0.874 with PROLIF_DISJOINT; the",
                        "difference between the MB1 partials is reported, which",
                        "one is doing the work is NOT interpreted."),
    no_outcome  = "NO outcome or survival variable enters at any point",
    no_treatment = paste("NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM: no",
                         "OXPHOS- or BH3-directed intervention exists in any",
                         "cohort available."),
    no_fdr = "no per-gene FDR, consistent with the rest of the project",
    n3 = "transcript scores and exposure axes; nothing is 'primed'",
    n_boot = N_BOOT, seed = PROJECT_SEED),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E33)
utils::write.csv(as.data.frame(COR), PATH_E33_COR, row.names = FALSE)
message("   ", PATH_E33)
message("   ", PATH_E33_COR)
message("\nE33 done. VERDICT: ", verdict, "\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E33)

  cat(x$verdict, "\n\n"); cat(x$reading, "\n")
  x$verdict_strict
  x$verdict_sensitivity %>% as.data.frame()
  utils::str(x$counts)

  # LIMB 3: both partials, both arms, always. Then the composite and genes.
  x$limb3$separation %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # The proliferation companion. POST-HOC; cannot change the verdict.
  x$limb3$companion %>%
    dplyr::transmute(set, label, n, rho = round(rho, 3),
                     ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3)) %>%
    as.data.frame()

  # Declared UNINFORMATIVE. Not support in either direction.
  x$limb3$uninformative %>%
    dplyr::transmute(label, rho = round(rho, 3)) %>% as.data.frame()

  # LIMB 1: MYC by forkscale MEDIAN-SPLIT quadrant, and by STORED fork.
  x$limb1$myc_by_forkscale_quadrant %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$limb1$myc_by_fork_assignment %>%
    dplyr::filter(estimator == "M_a") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$limb1$partials %>%
    dplyr::transmute(label, n, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
                     ci_hi = round(ci_hi, 3)) %>%
    as.data.frame()
  x$limb1$difference %>% as.data.frame()

  # LIMB 2: the variable map.
  names(x$limb2$crosstabs)
  x$limb2$crosstabs[["gsva: quadrant x MB2 fork"]] %>% as.data.frame()
  x$limb2$asserts %>% as.data.frame()
  x$limb2$sanity_prolifstrip$summary %>% as.data.frame()

  # The join, and the E32 reproduction.
  x$join$different_vial %>% as.data.frame()
  x$e32_reproduction %>% as.data.frame()
  x$fork_assignments %>% as.data.frame()
}
