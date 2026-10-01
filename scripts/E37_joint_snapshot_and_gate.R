# E37_joint_snapshot_and_gate.R
# =============================================================================
# REBUILD THE EXPRESSION LAYER WITH THE NORMALS, AND GATE IT ON E34
#
# EXPLORATORY. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-01_E37_declaration.md, COMMITTED ALONE AT 205ef8e
# before any rebuilt value existed. Data note, with the source audit, the counts
# and the column-key decision: docs/2026-10-01_e37_data.md at cc5dfe9. Read both
# first; this header restates their rules, it does not invent them.
#
# =============================================================================
# WHY THIS EXISTS
# =============================================================================
# E36 stopped: there are no solid-tissue-normal samples in any saved object in
# this repo. They exist - 113 of them - and were dropped at the first step of
# the upstream pipeline by `keep <- which(sty == "01")`, for a documented reason
# that does NOT apply to E36: primary tumours only, so expression and copy
# number stay aligned. E36 reads no copy number. This restores samples dropped
# for an unrelated constraint; it does not reverse a quality judgement.
#
# =============================================================================
# THE SNAPSHOT SITS BESIDE, AND IS SCOPED
# =============================================================================
# It does NOT replace data/from_validation/. Replacing it would invalidate E11,
# E33, E34 and E35 at a stroke. CLAUDE.md names two normalisations in one repo
# as a trap; the discipline that resolves it is that EVERY SCRIPT STATES WHICH
# SNAPSHOT IT READS and NO RESULT MIXES THE TWO.
#
#   THIS SCRIPT WRITES THE JOINT SNAPSHOT. IT IS READ BY E36 AND BY NOTHING
#   ELSE UNLESS SEPARATELY DECLARED.
#
# =============================================================================
# THE ONE DELIBERATE CHANGE
# =============================================================================
# The sample filter admits sample type 11 as well as 01. Type 06 (metastatic,
# n = 7) stays excluded, as upstream. Everything else follows
# myc_human_validation/scripts/01_fetch_tcga_expression.R, which is READ-ONLY.
#
# The log2-to-counts inversion is NOT repeated. Upstream cached its own
# intermediate BEFORE the sample filter, holding integer counts whose inversion
# was verified against four GDC per-file downloads. This script inherits that
# verified result (data note section 1).
#
# =============================================================================
# THE GATE - fixed in the declaration, transcribed here, not reinterpreted
# =============================================================================
# E34's headline quantities are recomputed on the rebuilt snapshot, RESTRICTED
# TO THE 1,095 TUMOURS, and compared with the committed values in
# docs/2026-10-01_e34_result.md. PASS requires ALL FOUR:
#
#   quadrant agreement with E34, each instrument   >= 95%
#   the two OXPHOS contrasts                       each within 0.10
#   the two MYC contrasts                          each within 0.10
#   E34's pooled per-cell labels                   unchanged (all four cells)
#
# plus the declared assertions (Q3 == STATE level 2, Q4 == levels 3 + 4, both
# instruments), which are a STOP regardless of the numeric thresholds.
#
# A FAIL IS NOT A REASON TO ADJUST THE THRESHOLDS. They were fixed before any
# rebuilt value existed. On a FAIL the script reports what moved and by how
# much, and E36 does not proceed.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - Nothing is written to myc_human_validation or myc_mouse. The cached counts
#     object is read with readRDS from the frozen repo and nothing else.
#   - NO outcome or survival variable enters. No causal language.
#   - It establishes nothing biological. It is infrastructure plus one
#     robustness check.
#   - It does NOT make the two snapshots comparable. The gate is a deliberate
#     like-for-like comparison of the SAME quantities and is the only place the
#     two are set side by side.
#
# N3: transcript scores. Nothing is "primed".
# SPECIES: human throughout. NO ORTHOLOG FUNCTION IS CALLED.
# SCALE: two matrices, opposite requirements, saved separately and never merged.
#   joint_tcga_linear.rds  LINEAR DESeq2-normalised - mitoPPS only, never logged
#   joint_tcga_vst.rds     VST, log scale - GSVA only, kcdf Gaussian
# COHORT-RELATIVITY: GSVA and mitoPPS are scored ONCE over the joint 1,208, which
# is the whole point. No score here is ever compared numerically with a score
# from data/from_validation/, SCAN-B, CCLE or mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "mitopps.R"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

message("\nE37: joint snapshot with normals, and the E34 gate\n", strrep("=", 78))

PATH_E37_LINEAR <- file.path(DIR_RESULTS, "joint_tcga_linear.rds")
PATH_E37_VST    <- file.path(DIR_RESULTS, "joint_tcga_vst.rds")
PATH_E37_SCORES <- file.path(DIR_RESULTS, "joint_tcga_scores.rds")
PATH_E37_GATE   <- file.path(DIR_RESULTS, "e37_gate.rds")
PATH_E37_README <- file.path(DIR_RESULTS, "joint_snapshot_README.md")
PATH_E37_CSV    <- file.path(DIR_TABLES,  "E37_gate_comparison.csv")

# =============================================================================
# 0. Constants, and THE GATE, transcribed from the declaration
# =============================================================================
# READ-ONLY inputs from the frozen repo. Nothing is written there, ever.
PATH_CACHED_COUNTS <- paste0("/Users/gs/code/myc_human_validation/data/raw/",
                             "tcga_brca_expression/tcga_brca_star_counts.rds")
PATH_STATE <- "/Users/gs/code/myc_human_validation/results/state_definition.rds"

KEEP_TYPES   <- c("01", "11")    # THE ONE CHANGE. 06 stays excluded, as upstream.
EXPECT_TUMOUR <- 1095L
EXPECT_NORMAL <- 113L
EXPECT_GENES_IN <- 60660L
EXPECT_COLS_IN  <- 1226L

ARM_OX     <- "OXPHOS subunits"
PROLIF_COV <- "PROLIF_DISJOINT"
MYC_WANT   <- c(MYC_REF, "FELSHER__PROLIFSTRIP")

TRIGGER  <- c("BBC3", "BID", "BIK", "BAD")
GUARDIAN <- c("BCL2L1", "MCL1")
CONFIG_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)

QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",
          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")
INSTR        <- c(gsva = "OX_gsva", mitopps = "OX_mitopps")
PAM50_LEVELS <- c("LumA", "LumB", "HER2", "Basal", "Normal")

CONTR <- tibble::tribble(
  ~contrast,  ~hi,     ~lo,     ~kind,
  "Q2 - Q1",  QLEV[2], QLEV[1], "OXPHOS",
  "Q4 - Q3",  QLEV[4], QLEV[3], "OXPHOS",
  "Q4 - Q2",  QLEV[4], QLEV[2], "MYC",
  "Q3 - Q1",  QLEV[3], QLEV[1], "MYC")

# --- E34's committed values, transcribed from docs/2026-10-01_e34_result.md --
# Section 4, m1 and m3, pooled, TCGA. These are the FIXED comparator; they are
# not recomputed from E34's object, which may or may not be on disk, though the
# object is cross-checked against them in section 8 when it is present.
E34_REF <- tibble::tribble(
  ~instrument, ~model, ~contrast,  ~e34,
  "gsva",      "m1",   "Q4 - Q2",  -0.002,
  "gsva",      "m1",   "Q3 - Q1",  +0.134,
  "gsva",      "m1",   "Q2 - Q1",  +0.981,
  "gsva",      "m1",   "Q4 - Q3",  +0.845,
  "gsva",      "m3",   "Q4 - Q2",  +0.115,
  "gsva",      "m3",   "Q3 - Q1",  +0.278,
  "gsva",      "m3",   "Q2 - Q1",  +1.005,
  "gsva",      "m3",   "Q4 - Q3",  +0.841,
  "mitopps",   "m1",   "Q4 - Q2",  +0.185,
  "mitopps",   "m1",   "Q3 - Q1",  +0.163,
  "mitopps",   "m1",   "Q2 - Q1",  +0.723,
  "mitopps",   "m1",   "Q4 - Q3",  +0.744,
  "mitopps",   "m3",   "Q4 - Q2",  +0.314,
  "mitopps",   "m3",   "Q3 - Q1",  +0.308,
  "mitopps",   "m3",   "Q2 - Q1",  +0.751,
  "mitopps",   "m3",   "Q4 - Q3",  +0.757)
# E34 read all four pooled cells as OXPHOS-ORDERED.
E34_LABELS <- c("gsva.m1" = "OXPHOS-ORDERED", "gsva.m3" = "OXPHOS-ORDERED",
                "mitopps.m1" = "OXPHOS-ORDERED", "mitopps.m3" = "OXPHOS-ORDERED")

GATE_MIN_AGREE <- 0.95
GATE_MAX_DELTA <- 0.10
GATE_RULE <- paste0(
  "PASS iff quadrant agreement with E34 >= ", GATE_MIN_AGREE * 100,
  "% on EACH instrument, AND all four contrasts within ", GATE_MAX_DELTA,
  " of E34's committed values in every pooled cell, AND E34's pooled per-cell ",
  "labels unchanged in all four. The declared assertions are a STOP regardless.")

message("\n0. the gate, fixed at 205ef8e before any rebuilt value existed:\n   ",
        GATE_RULE)

# =============================================================================
# 1. The source - read-only, and the inversion inherited
# =============================================================================
message("\n1. source (READ-ONLY, frozen repo)")
if (!file.exists(PATH_CACHED_COUNTS)) {
  stop("the cached upstream counts object is not at ", PATH_CACHED_COUNTS,
       ". Without it this rebuild cannot start, and re-fetching is not ",
       "authorised.", call. = FALSE)
}
src <- readRDS(PATH_CACHED_COUNTS)
for (el in c("counts", "annotation", "averaged_col")) {
  if (is.null(src[[el]])) stop("cached object lacks $", el, call. = FALSE)
}
counts_raw   <- src$counts
gene_ann     <- src$annotation
averaged_col <- src$averaged_col
if (!identical(dim(counts_raw), c(EXPECT_GENES_IN, EXPECT_COLS_IN))) {
  stop("cached counts are ", paste(dim(counts_raw), collapse = " x "), ", not ",
       EXPECT_GENES_IN, " x ", EXPECT_COLS_IN, ". The source has moved.",
       call. = FALSE)
}
if (!identical(storage.mode(counts_raw), "integer")) {
  stop("cached counts are not integer; the upstream inversion is not what it ",
       "is taken to be", call. = FALSE)
}
message("   ", nrow(counts_raw), " genes x ", ncol(counts_raw),
        " samples, integer; inversion inherited from upstream and NOT repeated")
message("   Xena SHA-256 recorded upstream: ", substr(src$sha256, 1, 16), "...")

# =============================================================================
# 2. Samples - upstream's rule, with sample type 11 admitted
# =============================================================================
# SAMPLE TYPE IS POSITIONS 14-15 OF THE BARCODE, read directly and never
# inferred. Upstream's tie-break is reproduced exactly: drop an averaged column
# only where that patient has a clean alternative, then keep the lowest vial.
message("\n2. samples")

bc   <- colnames(counts_raw)
pat  <- substr(bc, 1, 12)
sty  <- substr(bc, 14, 15)
vial <- substr(bc, 16, 16)
avg  <- bc %in% averaged_col
message("   types present: ",
        paste(sprintf("%s=%d", names(table(sty)), table(sty)), collapse = " "))
if (length(averaged_col) > 0L && !any(avg)) {
  stop("averaged columns were flagged upstream but none match these column ",
       "names. The flag list has been lost; the averaged aliquots would be ",
       "selected silently.", call. = FALSE)
}

keep_idx <- integer(0)
per_type <- list()
for (ty in KEEP_TYPES) {
  k  <- which(sty == ty)
  bk <- bc[k]; pk <- pat[k]; vk <- vial[k]; ak <- avg[k]
  n_in <- length(k)
  # drop an averaged column only where the patient has a clean alternative
  drop_avg <- ak & (pk %in% pk[!ak])
  k <- k[!drop_avg]; pk <- pk[!drop_avg]; vk <- vk[!drop_avg]
  n_avg <- sum(drop_avg)
  # then keep the lowest vial per patient
  o <- order(pk, vk); k <- k[o]; pk <- pk[o]
  first <- !duplicated(pk)
  n_dup <- sum(!first)
  k <- k[first]
  keep_idx <- c(keep_idx, k)
  per_type[[ty]] <- tibble::tibble(sample_type = ty, columns_in = n_in,
                                   averaged_dropped = n_avg,
                                   duplicates_dropped = n_dup,
                                   surviving = length(k))
  message("   type ", ty, ": ", n_in, " in, -", n_avg, " averaged, -", n_dup,
          " duplicate vial, ", length(k), " surviving")
}
SAMPLE_STEPS <- dplyr::bind_rows(per_type)
counts_raw <- counts_raw[, keep_idx, drop = FALSE]

n_t <- sum(substr(colnames(counts_raw), 14, 15) == "01")
n_n <- sum(substr(colnames(counts_raw), 14, 15) == "11")
if (n_t != EXPECT_TUMOUR || n_n != EXPECT_NORMAL) {
  stop("after filtering: ", n_t, " tumours and ", n_n, " normals, not ",
       EXPECT_TUMOUR, " and ", EXPECT_NORMAL, ". The data note's counts do not ",
       "reproduce; stop rather than build on a cohort nobody has explained.",
       call. = FALSE)
}

# THE COLUMN KEY. Patient barcodes are NOT unique once normals are in - 113
# patients carry both - so columns are keyed by the 16-character sample id and
# the pairing travels in sample_map rather than being re-derived downstream.
sample_map <- tibble::tibble(
  sample_id   = substr(colnames(counts_raw), 1, 16),
  patient     = substr(colnames(counts_raw), 1, 12),
  sample_type = substr(colnames(counts_raw), 14, 15),
  aliquot     = colnames(counts_raw))
sample_map$tissue <- ifelse(sample_map$sample_type == "01", "tumour", "normal")
if (anyDuplicated(sample_map$sample_id)) {
  stop("16-character sample ids are not unique", call. = FALSE)
}
colnames(counts_raw) <- sample_map$sample_id
message("   joint set: ", ncol(counts_raw), " columns (", n_t, " tumour, ",
        n_n, " normal); ", sum(duplicated(sample_map$patient)),
        " patients contribute both")

# =============================================================================
# 3. Genes - upstream's rules, applied to the joint set
# =============================================================================
message("\n3. genes")
common <- intersect(rownames(counts_raw), gene_ann$gene_id)
message("   gene ids matched to annotation: ", length(common), " of ",
        nrow(counts_raw))
if (length(common) < 0.95 * nrow(counts_raw)) {
  stop("fewer than 95% of gene ids matched the annotation", call. = FALSE)
}
counts_raw <- counts_raw[common, , drop = FALSE]
gene_ann   <- gene_ann[match(common, gene_ann$gene_id), ]
stopifnot(identical(rownames(counts_raw), gene_ann$gene_id))

is_pc <- gene_ann$gene_type == "protein_coding" &
         !is.na(gene_ann$gene_name) & gene_ann$gene_name != ""
message("   protein-coding: ", sum(is_pc), " of ", nrow(gene_ann))
counts <- counts_raw[is_pc, , drop = FALSE]
sym    <- gene_ann$gene_name[is_pc]
rm(counts_raw); invisible(gc(verbose = FALSE))

# Collapse duplicate symbols by HIGHEST MEAN COUNT, not by summing - upstream's
# rule, and summing would inflate genes carrying several ENSG ids.
if (any(duplicated(sym))) {
  mu  <- rowMeans(counts)
  ord <- order(sym, -mu)
  counts <- counts[ord, , drop = FALSE]; sym <- sym[ord]
  first <- !duplicated(sym)
  message("   ", sum(!first), " duplicate symbol row(s) collapsed (highest mean)")
  counts <- counts[first, , drop = FALSE]; sym <- sym[first]
}
rownames(counts) <- sym
stopifnot(!any(duplicated(rownames(counts))))

# The low-count filter is computed over the JOINT set, so it differs from
# upstream by construction. That is the point, and the gate is what catches it.
keep_g <- rowSums(counts >= 10) >= 10
message("   low-count filter (>=10 counts in >=10 samples): ", sum(keep_g),
        " kept, ", sum(!keep_g), " dropped")
counts <- counts[keep_g, , drop = FALSE]
message("   final matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# =============================================================================
# 4. The two matrices - DESeq2 size factors over the JOINT set
# =============================================================================
# The whole point of the rebuild. Size factors over 1,208 are not size factors
# over 1,095, so every tumour value moves - which the gate measures.
message("\n4. normalisation, over the joint set")
cd  <- data.frame(sample_id = colnames(counts), row.names = colnames(counts))
dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)

mat_linear <- DESeq2::counts(dds, normalized = TRUE)
stopifnot(min(mat_linear) >= 0)
message("   LINEAR (mitoPPS only, never logged): range ",
        sprintf("%.1f", min(mat_linear)), " to ", sprintf("%.0f", max(mat_linear)))

vsd <- DESeq2::vst(dds, blind = TRUE)
mat_vst <- SummarizedExperiment::assay(vsd)
message("   VST (GSVA only, kcdf Gaussian): range ",
        sprintf("%.2f", min(mat_vst)), " to ", sprintf("%.2f", max(mat_vst)))
stopifnot(identical(colnames(mat_linear), colnames(mat_vst)),
          identical(rownames(mat_linear), rownames(mat_vst)))

saveRDS(list(mat = mat_linear, scale = "linear_deseq2_normalised",
             consumer = "mitoPPS only", sample_map = sample_map,
             snapshot = "JOINT tumour+normal, E37", built = Sys.time()),
        PATH_E37_LINEAR)
saveRDS(list(mat = mat_vst, scale = "vst_log", consumer = "GSVA only",
             sample_map = sample_map,
             snapshot = "JOINT tumour+normal, E37", built = Sys.time()),
        PATH_E37_VST)
message("   saved ", basename(PATH_E37_LINEAR), " and ", basename(PATH_E37_VST))

# =============================================================================
# 5. Scoring - every definition read from objects already in this repo
# =============================================================================
message("\n5. scoring over the joint set (both instruments are cohort-relative)")

sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
mref <- readRDS(PATH_TCGA_MITO)          # for mito_paths and arm_universe_path only
arm_sets  <- sd_$arm_sets
cov_sets  <- sd_$cov_sets
myc_sets  <- sd_$myc_sets[MYC_WANT]
if (any(vapply(myc_sets, is.null, logical(1)))) {
  stop("MYC set(s) absent from set_definitions: ",
       paste(MYC_WANT[vapply(myc_sets, is.null, logical(1))], collapse = ", "),
       call. = FALSE)
}
paths     <- mref$mito_paths
arm_univ  <- mref$arm_universe_path

# --- 5.1 GSVA, universe pinned, all sets in ONE call -------------------------
# CLAUDE.md: score all sets of interest in the SAME call; the .PIN half-matrix
# sets hold the gene universe so two calls cannot drift apart.
.gsva_batch <- function(sets, E, label) {
  n_ok <- vapply(sets, length, integer(1))
  if (any(n_ok < MIN_SET_GENES)) {
    stop("GSVA batch '", label, "': set(s) below the size floor -> ",
         paste(names(sets)[n_ok < MIN_SET_GENES], collapse = ", "), call. = FALSE)
  }
  wanted <- names(sets)
  syms <- rownames(E)
  sets[[".PIN_A"]] <- syms[c(TRUE, FALSE)]
  sets[[".PIN_B"]] <- syms[c(FALSE, TRUE)]
  par <- GSVA::gsvaParam(exprData = E, geneSets = sets, kcdf = "Gaussian",
                         minSize = MIN_SET_GENES, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)
  dropped <- setdiff(wanted, rownames(s))
  if (length(dropped)) {
    stop("GSVA silently dropped set(s) in '", label, "': ",
         paste(utils::head(dropped, 10), collapse = ", "), call. = FALSE)
  }
  s[wanted, , drop = FALSE]
}
.in_matrix <- function(g) intersect(unique(g), rownames(mat_vst))
all_sets <- c(lapply(arm_sets, .in_matrix), lapply(cov_sets, .in_matrix),
              lapply(myc_sets, .in_matrix))
all_sets <- all_sets[vapply(all_sets, length, integer(1)) >= MIN_SET_GENES]
message("   GSVA: ", length(all_sets), " sets x ", ncol(mat_vst),
        " samples; this takes a few minutes")
gsva_all <- .gsva_batch(all_sets, mat_vst, "joint")
message("   scored ", nrow(gsva_all), " sets")

# --- 5.2 mitoPPS, E02's construction, from functions/mitopps.R ---------------
message("   mitoPPS (linear DESeq2-normalised)")
paths_m <- lapply(paths, .in_matrix)
S_univ  <- .path_scores(paths_m, mat_linear, MIN_SET_GENES)
if (min(S_univ) <= 0) {
  stop("a MitoPathway score is zero or negative in some sample; the pairwise ",
       "ratio is undefined", call. = FALSE)
}
mpps_univ <- .mitopps_universe(S_univ)
message(sprintf("   universe: %d pathways; global mean %.4f (should be ~1)",
                nrow(S_univ), mean(mpps_univ)))
arm_m   <- lapply(arm_sets, .in_matrix)
S_arms  <- .path_scores(arm_m, mat_linear, MIN_SET_GENES)
# Each arm is queried against the universe with its OWN pathway held out, so an
# arm that IS a universe pathway reproduces its canonical value exactly.
mitopps_arms <- t(vapply(rownames(S_arms), function(a) {
  hold <- arm_univ[[a]]
  U <- if (is.null(hold) || is.na(hold)) S_univ
       else S_univ[setdiff(rownames(S_univ), hold), , drop = FALSE]
  as.numeric(.mitopps_query(S_arms[a, , drop = FALSE], U))
}, numeric(ncol(mat_linear))))
colnames(mitopps_arms) <- colnames(mat_linear)

saveRDS(list(gsva = gsva_all, mitopps_arms = mitopps_arms,
             mitopps_universe = mpps_univ, sample_map = sample_map,
             sets_used = names(all_sets),
             snapshot = "JOINT tumour+normal, E37",
             note = "cohort-relative over the JOINT 1,208; never compared with data/from_validation/",
             built = Sys.time()), PATH_E37_SCORES)
message("   saved ", basename(PATH_E37_SCORES))

# =============================================================================
# 6. THE QUADRANT, tumours only, from the frozen constructor
# =============================================================================
message("\n6. the quadrant (tumours only; the frozen constructor, re-applied)")

st <- readRDS(PATH_STATE)
if (!identical(deparse(st$build_state), st$definition_source$build_state) ||
    !identical(environment(st$build_state), baseenv())) {
  stop("the frozen STATE constructor fails its own integrity contract",
       call. = FALSE)
}
.build_state <- st$build_state
LV <- st$spec$levels
S0 <- st$tcga$state

TUM <- sample_map$sample_id[sample_map$tissue == "tumour"]
tum_pat <- sample_map$patient[match(TUM, sample_map$sample_id)]
si <- match(tum_pat, S0$patient)
if (anyNA(si)) stop("a rebuilt tumour has no frozen STATE row", call. = FALSE)

# The constructor's inputs, rebuilt: M_a, the OXPHOS arm on both instruments,
# and BUFFER_gistic carried unchanged from the frozen object (it is copy number
# and this rebuild does not touch genomics).
QD <- tibble::tibble(
  sample_id  = TUM,
  patient    = tum_pat,
  MYC        = as.numeric(gsva_all[MYC_REF, TUM]),
  OX_gsva    = as.numeric(gsva_all[ARM_OX, TUM]),
  OX_mitopps = as.numeric(mitopps_arms[ARM_OX, TUM]),
  BUFFER     = S0$BUFFER_gistic[si])

myc_high <- NULL
ASSERTS <- list()
for (inst in names(INSTR)) {
  ox <- QD[[INSTR[[inst]]]]
  state_new <- .build_state(QD$MYC, ox, QD$BUFFER)
  if (is.null(myc_high)) {
    myc_high <- state_new != LV[1]
  } else if (!identical(myc_high, state_new != LV[1])) {
    stop("the two instruments disagree on the MYC call", call. = FALSE)
  }
  probe <- .build_state(ifelse(is.na(QD$MYC), NA_real_, ox), ox, QD$BUFFER)
  if (any(probe == LV[2], na.rm = TRUE)) {
    stop("level 2 appeared in the myc := oxphos probe", call. = FALSE)
  }
  if (!identical(is.na(probe), is.na(state_new))) {
    stop("the probe's complete-case set differs from STATE's", call. = FALSE)
  }
  ox_high <- probe != LV[1]
  q <- dplyr::case_when(
    !myc_high & !ox_high ~ QLEV[1],
    !myc_high &  ox_high ~ QLEV[2],
     myc_high & !ox_high ~ QLEV[3],
     myc_high &  ox_high ~ QLEV[4])
  # THE DECLARED ASSERTIONS - a STOP regardless of the numeric thresholds.
  a3 <- identical(which(q == QLEV[3]), which(state_new == LV[2]))
  a4 <- identical(which(q == QLEV[4]), which(state_new %in% LV[3:4]))
  if (!a3 || !a4) {
    stop("ASSERTION FAILED (", inst, "): Q3 == level 2 is ", a3,
         ", Q4 == levels 3 + 4 is ", a4, ". The quadrant is not what it is ",
         "taken to be on the rebuilt scale. STOP; E36 does not proceed.",
         call. = FALSE)
  }
  ASSERTS[[inst]] <- tibble::tibble(instrument = inst, Q3_equals_level2 = a3,
                                    Q4_equals_levels3_4 = a4)
  QD[[paste0("quad_", inst)]] <- factor(q, levels = QLEV)
  message("   ", inst, ": Q3 == level 2 TRUE; Q4 == levels 3 + 4 TRUE")
}
ASSERTS <- dplyr::bind_rows(ASSERTS)

# =============================================================================
# 7. The configuration, and the model frame - E34's construction
# =============================================================================
message("\n7. the configuration composite and the model frame")

gm <- mat_linear[, TUM, drop = FALSE]
miss <- setdiff(names(CONFIG_SIGN), rownames(gm))
if (length(miss)) stop("configuration gene(s) absent: ",
                       paste(miss, collapse = ", "), call. = FALSE)
GL <- log2(gm[names(CONFIG_SIGN), , drop = FALSE] + 1)
gv <- apply(GL, 1L, stats::var)
if (any(!(gv > 0))) stop("zero-variance configuration gene", call. = FALSE)
ZG <- t(scale(t(GL)))
config_raw <- as.numeric(
  colMeans(sweep(ZG, 1L, CONFIG_SIGN[rownames(ZG)], "*")))

.blom <- function(x) {
  n <- length(x); stopifnot(!anyNA(x))
  as.numeric(scale(stats::qnorm((rank(x, ties.method = "average") - 0.375) /
                                  (n + 0.25))))
}
.blom_na <- function(x) {
  out <- rep(NA_real_, length(x)); ok <- !is.na(x)
  if (any(ok)) out[ok] <- .blom(x[ok]); out
}
fr <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
fr <- fr[fr$cohort == "TCGA", ]
fi <- match(QD$patient, fr$sample_id)
# Computed before the tibble: inside tibble() a new column shadows an outer
# object of the same name - the trap that stopped E32's first run.
D <- tibble::tibble(
  sample_id = QD$sample_id,
  quad_gsva = QD$quad_gsva, quad_mitopps = QD$quad_mitopps,
  PAM50  = factor(as.character(fr$PAM50[fi]), levels = PAM50_LEVELS),
  purity = .blom_na(fr$purity[fi]),
  leuko  = .blom_na(fr$leuko[fi]),
  PROLIF = .blom(as.numeric(gsva_all[PROLIF_COV, TUM])),
  config = .blom(config_raw))
message("   ", nrow(D), " tumours; PAM50 NA ", sum(is.na(D$PAM50)),
        "; purity and leuko present ",
        sum(!is.na(D$purity) & !is.na(D$leuko)))

# =============================================================================
# 8. THE GATE
# =============================================================================
message("\n8. the gate\n", strrep("-", 78))

.fit_groups <- function(dat, covars) {
  cols <- c("config", "quad", covars)
  dd <- as.data.frame(dat[stats::complete.cases(dat[, cols]), cols, drop = FALSE])
  dd$quad <- droplevels(factor(dd$quad, levels = QLEV))
  for (k in covars) if (is.factor(dd[[k]])) dd[[k]] <- droplevels(dd[[k]])
  f <- stats::lm(stats::reformulate(c("quad", covars), response = "config"),
                 data = dd)
  tt <- stats::delete.response(stats::terms(f))
  A <- vapply(QLEV, function(q) {
    dq <- dd; dq$quad <- factor(rep(q, nrow(dd)), levels = levels(dd$quad))
    colMeans(stats::model.matrix(tt, data = dq, xlev = f$xlevels,
                                 contrasts.arg = f$contrasts))
  }, numeric(length(stats::coef(f))))
  list(n = nrow(dd), df = f$df.residual, A = A, b = stats::coef(f),
       V = stats::vcov(f))
}
.lincomb <- function(fit, a) {
  est <- sum(a * fit$b); se <- sqrt(drop(t(a) %*% fit$V %*% a))
  tq <- stats::qt(0.975, fit$df)
  c(est = est, lo = est - tq * se, hi = est + tq * se)
}
.label_cell <- function(ct) {
  g <- function(nm, col) ct[[col]][ct$contrast == nm]
  ox1 <- g("Q2 - Q1", "est"); ox2 <- g("Q4 - Q3", "est")
  mp  <- g("Q4 - Q2", "est"); mc  <- g("Q3 - Q1", "est")
  ox_pos <- ox1 > 0 && g("Q2 - Q1", "lo") > 0 && ox2 > 0 && g("Q4 - Q3", "lo") > 0
  ox_min <- min(ox1, ox2)
  if (ox_pos && abs(mp) < ox_min && abs(mc) < ox_min) return("OXPHOS-ORDERED")
  if (mp > 0 && g("Q4 - Q2", "lo") > 0 && mp >= ox_min) return("MYC-DEPENDENT")
  "neither"
}

COVARS <- list(m1 = c("PAM50", "purity", "leuko"),
               m3 = c("PAM50", "purity", "leuko", "PROLIF"))
NEW <- list(); LAB <- character(0)
for (inst in names(INSTR)) {
  d0 <- D; d0$quad <- d0[[paste0("quad_", inst)]]
  for (mod in names(COVARS)) {
    fit <- .fit_groups(d0, COVARS[[mod]])
    C <- t(vapply(seq_len(nrow(CONTR)), function(i)
      .lincomb(fit, fit$A[, CONTR$hi[i]] - fit$A[, CONTR$lo[i]]), numeric(3)))
    ct <- tibble::tibble(instrument = inst, model = mod,
                         contrast = CONTR$contrast, kind = CONTR$kind,
                         n = fit$n, est = C[, "est"], lo = C[, "lo"],
                         hi = C[, "hi"])
    NEW[[length(NEW) + 1L]] <- ct
    LAB[paste(inst, mod, sep = ".")] <- .label_cell(ct)
  }
}
NEW <- dplyr::bind_rows(NEW)

GATE_CONTRASTS <- NEW %>%
  dplyr::left_join(E34_REF, by = c("instrument", "model", "contrast")) %>%
  dplyr::mutate(delta = est - e34, within = abs(delta) <= GATE_MAX_DELTA)

# --- quadrant agreement against E34 -----------------------------------------
PATH_E34_OBJ <- file.path(DIR_RESULTS, "e34_quadrant_configuration.rds")
AGREE <- NULL; e34_ref_check <- NULL
if (file.exists(PATH_E34_OBJ)) {
  e34 <- readRDS(PATH_E34_OBJ)
  af  <- e34$analysis_frames$TCGA
  j   <- match(QD$patient, af$sample_id)
  AGREE <- dplyr::bind_rows(lapply(names(INSTR), function(inst) {
    old <- af[[paste0("quad_", inst)]][j]
    new <- QD[[paste0("quad_", inst)]]
    ok  <- !is.na(old) & !is.na(new)
    tibble::tibble(instrument = inst, n_comparable = sum(ok),
                   agreement = mean(as.character(old[ok]) ==
                                      as.character(new[ok])))
  }))
  # Cross-check the transcribed comparator against E34's own object.
  e34_ref_check <- e34$contrasts %>%
    dplyr::filter(readout_col == "config_comp", subset == "all",
                  model %in% c("m1", "m3"), cohort == "TCGA") %>%
    dplyr::transmute(instrument, model, contrast, obj = round(est, 3)) %>%
    dplyr::left_join(E34_REF, by = c("instrument", "model", "contrast")) %>%
    dplyr::mutate(matches_note = abs(obj - e34) < 1e-9)
} else {
  message("   NOTE: E34's object is absent, so quadrant agreement cannot be ",
          "computed.\n   The contrast thresholds still apply and the gate ",
          "FAILS for want of the check.")
}

labels_unchanged <- identical(LAB[names(E34_LABELS)], E34_LABELS)
agree_ok  <- !is.null(AGREE) && all(AGREE$agreement >= GATE_MIN_AGREE)
within_ok <- all(GATE_CONTRASTS$within)
gate <- if (agree_ok && within_ok && labels_unchanged) "PASS" else "FAIL"

message("RULE: ", GATE_RULE, "\n")
if (!is.null(AGREE)) {
  message("   quadrant agreement with E34:")
  AGREE %>% dplyr::mutate(agreement = round(agreement, 4)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}
message("\n   contrasts, rebuilt against E34's committed values:")
GATE_CONTRASTS %>%
  dplyr::transmute(instrument, model, contrast, kind,
                   rebuilt = sprintf("%+.3f", est), e34 = sprintf("%+.3f", e34),
                   delta = sprintf("%+.3f", delta), within) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   pooled per-cell labels: ",
        paste(sprintf("%s=%s", names(LAB), LAB), collapse = "; "))
message("   labels unchanged from E34: ", labels_unchanged)
message("\n   GATE: ", gate)

gate_text <- if (identical(gate, "PASS")) paste0(
  "PASS. The rebuild is sound AND E34 is robust to normalisation - its\n",
  "   headline contrasts survive a change of size factors and gene filter that\n",
  "   moves every input value. E36 PROCEEDS.") else paste0(
  "FAIL. STOP. E36 DOES NOT PROCEED.\n",
  "   This is a finding about E34, not a nuisance, and the thresholds are NOT\n",
  "   adjusted - they were fixed at 205ef8e before any rebuilt value existed.\n",
  "   The table above says what moved and by how much; that is what gets\n",
  "   written up, and E36 waits.")
message("\n", gate_text)

# =============================================================================
# 9. Save, and the snapshot README
# =============================================================================
message("\n9. saving")

out <- list(
  gate = gate, gate_text = gate_text, gate_rule = GATE_RULE,
  thresholds = list(min_agreement = GATE_MIN_AGREE, max_delta = GATE_MAX_DELTA),
  contrasts = GATE_CONTRASTS, agreement = AGREE,
  labels_new = LAB, labels_e34 = E34_LABELS,
  labels_unchanged = labels_unchanged,
  e34_reference_check = e34_ref_check,
  asserts = ASSERTS, sample_steps = SAMPLE_STEPS, sample_map = sample_map,
  counts = list(genes = nrow(mat_linear), samples = ncol(mat_linear),
                tumours = n_t, normals = n_n,
                patients_with_both = sum(duplicated(sample_map$patient))),
  spec = list(
    declaration = "docs/2026-10-01_E37_declaration.md, committed ALONE at 205ef8e",
    data_note   = "docs/2026-10-01_e37_data.md, cc5dfe9",
    change = "the sample filter admits type 11 as well as 01; type 06 stays excluded",
    inversion = paste("the log2-to-counts inversion is INHERITED from the",
                      "upstream cache, not repeated; it was verified upstream",
                      "against four GDC per-file downloads"),
    scope = paste("this snapshot sits BESIDE data/from_validation/ and is read",
                  "by E36 and by nothing else unless separately declared. NO",
                  "result mixes the two."),
    key = paste("columns keyed by the 16-character sample id; patient barcodes",
                "are not unique once normals are in"),
    joint = paste("DESeq2 size factors, the low-count filter, GSVA and mitoPPS",
                  "are ALL computed over the joint 1,208, which is the point"),
    no_outcome = "NO outcome or survival variable enters at any point",
    n3 = "transcript scores; nothing is 'primed'"),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
saveRDS(out, PATH_E37_GATE)
utils::write.csv(as.data.frame(GATE_CONTRASTS), PATH_E37_CSV, row.names = FALSE)

writeLines(c(
  "# Joint TCGA snapshot (tumour + solid-tissue-normal) - E37",
  "",
  "**What it is for: read by E36 and by nothing else unless separately declared.**",
  "",
  paste0("Built ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
         " by scripts/E37_joint_snapshot_and_gate.R."),
  "",
  "## Provenance",
  "",
  paste0("- Source: `", PATH_CACHED_COUNTS, "`, READ-ONLY, in the frozen repo."),
  paste0("- Upstream Xena SHA-256: `", src$sha256, "`"),
  "- The log2-to-counts inversion is INHERITED from that cache, not repeated.",
  "- Pipeline follows `myc_human_validation/scripts/01_fetch_tcga_expression.R`",
  "  with ONE change: the sample filter admits type 11 as well as 01.",
  "  Type 06 (metastatic) stays excluded.",
  "",
  "## Contents",
  "",
  paste0("- ", nrow(mat_linear), " genes x ", ncol(mat_linear), " samples (",
         n_t, " tumour, ", n_n, " normal; ",
         sum(duplicated(sample_map$patient)), " patients contribute both)."),
  "- Columns keyed by the 16-character sample id. Patient barcodes are NOT",
  "  unique here; `sample_map` carries patient, sample_id, sample_type, tissue.",
  "- `joint_tcga_linear.rds` LINEAR DESeq2-normalised. mitoPPS only. Never log it.",
  "- `joint_tcga_vst.rds` VST, log scale. GSVA only, kcdf Gaussian.",
  "- `joint_tcga_scores.rds` GSVA (arms, covariates, MYC sets) and mitoPPS arms.",
  "",
  "## The discipline that makes two normalisations safe",
  "",
  "This snapshot does NOT replace `data/from_validation/`, which E11, E33, E34",
  "and E35 rest on. Every script states which snapshot it reads, and **no",
  "result mixes the two**. Size factors, the low-count filter and both scoring",
  "instruments are computed over the joint 1,208, so no value here is",
  "comparable with a value from the other snapshot.",
  "",
  paste0("## Gate against E34: ", gate),
  "",
  "See `docs/2026-10-01_e37_data.md` and the E37 result note."), PATH_E37_README)

message("   ", PATH_E37_GATE, "\n   ", PATH_E37_CSV, "\n   ", PATH_E37_README)
message("\nE37 done. GATE: ", gate, "\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E37_GATE)

  cat(x$gate, "\n\n"); cat(x$gate_text, "\n")
  utils::str(x$counts)

  # The gate's two numeric legs.
  x$agreement %>% as.data.frame()
  x$contrasts %>%
    dplyr::transmute(instrument, model, contrast, kind, est = round(est, 3),
                     e34, delta = round(delta, 3), within) %>%
    as.data.frame()

  # The labels, and the declared assertions.
  x$labels_new; x$labels_e34; x$labels_unchanged
  x$asserts %>% as.data.frame()

  # Was the transcribed E34 comparator right? (needs E34's object on disk)
  x$e34_reference_check %>% as.data.frame()

  # How the sample filter went, step by step.
  x$sample_steps %>% as.data.frame()
  x$sample_map %>% dplyr::count(tissue) %>% as.data.frame()
}
