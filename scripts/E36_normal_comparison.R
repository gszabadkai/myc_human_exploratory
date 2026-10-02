# E36_normal_comparison.R
# =============================================================================
# IS THE OXPHOS-TO-CONFIGURATION COUPLING PRESENT IN NORMAL BREAST?
#
# EXPLORATORY. Nothing here is pre-registered.
#
# Declaration:  docs/2026-10-01_E36_declaration.md, COMMITTED ALONE at d71fa45
#               before any sample was opened.
# Amendment:    docs/2026-10-01_E36_declaration_amendment.md, ALONE at be64f9a.
#               The comparison is PAIRED; Q2 becomes within-patient; Q1's
#               contrast is bootstrapped by PATIENT; EPITHELIAL fixed by gene.
# Data note 1:  docs/2026-10-01_e36_data.md at 4802508 - the STOP, and the
#               positive control named before any fit.
# Data note 2:  docs/2026-10-01_e36_data_part2.md - M_b, the two pre-figure
#               checks, and THE THREE COMPARATORS forced by the site skew.
# Read all four before this file. This header restates their rules; it invents
# only the numeric thresholds in section 0, which are fixed HERE, before any
# Q1 or Q2 value has been computed.
#
# THIS IS THE LOAD-BEARING ANALYSIS. It can come back against the model, and
# PRESENT IN NORMAL is a declared reading, not a failure mode.
#
# =============================================================================
# THE SNAPSHOT - SAY IT ONCE, LOUDLY
# =============================================================================
# This script reads the JOINT snapshot built by E37 and NOTHING ELSE:
#   results/joint_tcga_linear.rds   LINEAR DESeq2-normalised, mitoPPS only
#   results/joint_tcga_vst.rds      VST, log scale, GSVA and ULM only
#   results/joint_tcga_scores.rds   GSVA arms + covariates, mitoPPS arms
# Size factors, the low-count filter and both instruments were computed over
# the JOINT 1,208. NO VALUE HERE IS COMPARABLE WITH A VALUE FROM
# data/from_validation/, and no result mixes the two. E37 section 2 declared
# this; E34's and E35's numbers are NOT a comparator for anything below.
#
# =============================================================================
# WHAT IS FORBIDDEN
# =============================================================================
#   - No outcome, survival or treatment variable. None is read or joined.
#   - No product term. The MYC x OXPHOS interaction estimand is not estimated.
#   - No causal language, no killing language. N3: these are transcript scores;
#     nothing is "primed", in any sample of any tissue.
#   - No per-gene FDR, and no p-value plucked from the grid.
#   - Nothing is written to myc_human_validation or myc_mouse.
#   - No PAM50 and no purity on a normal sample - neither exists for one.
#     ABSOLUTE purity is tumour-only; that is why EPITHELIAL exists.
#
# =============================================================================
# THE TWO QUESTIONS
# =============================================================================
# Q1 (PRIMARY, carries the verdict)
#   Partial Spearman of OXPHOS against the configuration composite, computed
#   WITHIN each tissue, then contrasted. A correlation within a tissue is
#   unaffected by the pairing; the CONTRAST between them is not, because the
#   same 113 patients contribute to both - hence the patient-level bootstrap.
#
# Q2 (DESCRIPTIVE, NO VERDICT)
#   Where tumour MYC sits relative to normal breast. Primary is the
#   WITHIN-PATIENT difference over the pairs; the distributional form the
#   premise is stated in is kept and also reported.
#
# =============================================================================
# THE THREE COMPARATORS - data note 2 section 4, forced by the site skew
# =============================================================================
# The 113 normals come from 6 of 40 collection sites, 87.6% from the top 3,
# against 20.8% of the other 982 tumours (Fisher p = 1e-05). Site BH alone
# supplies 73 of 113. So every ACROSS-PERSON comparison is computed on three
# nested tumour sets and ALL THREE ARE ALWAYS REPORTED:
#   all      1,095 tumours. The basis E34 and E35 used. Confounded with site.
#   site     the 394 at the 6 normal-contributing sites. Same places, not the
#            same proportions - BH is 38% of these and 65% of the normals.
#   paired   the 113 patient-matched. Tightest; its interval is as wide as the
#            normal side's.
# A READING IS TAKEN ONLY IF IT HOLDS IN ALL THREE. If they disagree, THE
# DISAGREEMENT IS THE RESULT and is reported as such, never resolved by
# choosing one. Site-matching cannot remove site, because the normals are
# concentrated inside the matched set too.
#
# =============================================================================
# THE POSITIVE CONTROL - A STOP CONDITION, NOT A ROBUSTNESS CHECK
# =============================================================================
# DECLARED PRIMARY, at 4802508 before any fit: Fe-S cluster assembly, CYTOSOLIC
# HALF. Reported alongside and NOT promotable into its place: the mitochondrial
# half (not a control - it is mitochondrial, so coupling to OXPHOS is close to
# guaranteed; it shows the instrument responds), mitophagy, and mitophagy
# PINK1/PRKN only.
#
# The control substitutes its half for the configuration and asks Q1 again,
# unchanged - same exposure, same covariates, same three comparators.
#
# PRE-SPECIFIED BUT NOT BLIND, disclosed at be64f9a and restated here: in
# TUMOURS the Fe-S cytosolic half against OXPHOS is +0.152 raw / +0.118
# proliferation-adjusted in TCGA and +0.214 / +0.159 in SCAN-B (E14, checked
# against results/curated_comparators.rds). The control asks whether that weak
# coupling CHANGES BETWEEN TISSUES AS MUCH AS THE CONFIGURATION'S DOES. It does
# NOT ask whether it is zero.
#
# =============================================================================
# SCALE DISCIPLINE - CLAUDE.md, and the two instruments never share an object
# =============================================================================
#   GSVA and ULM    read joint_tcga_vst.rds        (VST, log scale)
#   mitoPPS         read joint_tcga_linear.rds     (linear DESeq2-normalised)
#   all composites  log2(linear + 1), then z ACROSS THE JOINT 1,208 so both
#                   tissues sit on one scale. Spearman is invariant to that
#                   choice, so it changes no correlation; it matters only for
#                   the descriptive means, and it is stated where they appear.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

message("\nE36: is the coupling present in normal breast?\n", strrep("=", 78))

PATH_JOINT_LINEAR <- file.path(DIR_RESULTS, "joint_tcga_linear.rds")
PATH_JOINT_VST    <- file.path(DIR_RESULTS, "joint_tcga_vst.rds")
PATH_JOINT_SCORES <- file.path(DIR_RESULTS, "joint_tcga_scores.rds")
PATH_STATE   <- "/Users/gs/code/myc_human_validation/results/state_definition.rds"
PATH_E14     <- file.path(DIR_RESULTS, "curated_comparators.rds")

PATH_E36         <- file.path(DIR_RESULTS, "e36_normal_comparison.rds")
PATH_E36_CSV_Q1  <- file.path(DIR_TABLES,  "E36_q1_coupling.csv")
PATH_E36_CSV_Q2  <- file.path(DIR_TABLES,  "E36_q2_myc_level.csv")
PATH_E36_FIG     <- file.path(DIR_FIGURES, "E36_normal_comparison.png")
PATH_E36_FIG_PDF <- file.path(DIR_FIGURES, "E36_normal_comparison.pdf")
DIR_DOC_FIGURES  <- here::here("docs", "figures")
PATH_E36_FIG_DOC     <- file.path(DIR_DOC_FIGURES, "E36_normal_comparison.png")
PATH_E36_FIG_DOC_PDF <- file.path(DIR_DOC_FIGURES, "E36_normal_comparison.pdf")

# =============================================================================
# 0. CONSTANTS, AND THE READING RULE - fixed HERE, before any value exists
# =============================================================================
N_BOOT   <- 10000L
CI_LEVEL <- 0.95
SEED     <- 20261001L

ARM_OX     <- "OXPHOS subunits"
PROLIF_COV <- "PROLIF_DISJOINT"
INSTR      <- c(gsva = "OX_gsva", mitopps = "OX_mitopps")

# The configuration: E33's composite, verbatim, signs included.
CONFIG_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)
TRIGGER  <- c("BBC3", "BID", "BIK", "BAD")
GUARDIAN <- c("BCL2L1", "MCL1")

# EPITHELIAL, fixed by gene list in the amendment, section 4. Both compartments
# deliberately: a luminal-only score would be part composition estimate and
# part subtype score, and would adjust away some of the thing being measured.
EPI_GENES <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19",   # luminal
               "KRT5", "KRT14", "KRT17", "TP63")            # basal / myoepi
# Descriptive companions ONLY. They enter NO model.
ADIPOSE_GENES    <- c("ADIPOQ", "FABP4", "PLIN1", "LEP", "CIDEC")
FIBROBLAST_GENES <- c("COL1A1", "COL1A2", "FAP", "PDGFRB", "THY1")

QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",
          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")

COMPARATOR_LEV <- c("all", "site", "paired")

# Q2's estimators. All four CollecTRI variants, never one under the bare name
# M_b (E05's rule). M_b__MITOSTRIP is the CONSTRUCTION counterpart of M_a -
# both are mito-stripped - which is a statement about construction and NOT a
# promotion; none of the six is designated primary.
MB_VARIANTS <- c("M_b__FULL", "M_b__MITOSTRIP", "M_b__PROLIFSTRIP",
                 "M_b__BOTHSTRIP")

# --- the readouts ------------------------------------------------------------
# PRIMARY carries the verdict. CONTROL is the declared stop condition.
# COMPANION is reported always and can never be promoted into either place.
READOUTS <- tibble::tribble(
  ~readout,              ~role,       ~signed,
  "config",              "PRIMARY",   TRUE,
  "fes_cyto",            "CONTROL",   FALSE,
  "fes_mito",            "COMPANION", FALSE,
  "mitophagy_cyto",      "COMPANION", FALSE,
  "mitophagy_mito",      "COMPANION", FALSE,
  "mitophagy_pp_cyto",   "COMPANION", FALSE,
  "mitophagy_pp_mito",   "COMPANION", FALSE)

# Adjustments. The DECLARED one is EPITHELIAL, reported with and without
# (declaration section 2a). The proliferation-adjusted column is a COMPANION:
# it is standard everywhere else in this repo and is reported always, but the
# VERDICT READS THE DECLARED ADJUSTMENT ONLY. If the companion disagrees with
# the declared one, that disagreement is reported as a qualification in the
# same breath as the verdict - it is not allowed to replace it, and it is not
# allowed to be omitted.
ADJUST <- list(
  raw        = character(0),
  epithelial = "EPITHELIAL",
  epi_prolif = c("EPITHELIAL", "PROLIF"))
ADJUST_DECLARED <- c("raw", "epithelial")

# Bootstrapped: the PRIMARY and the declared CONTROL, on all three adjustments.
# The COMPANION readouts get point estimates only - they never carry a verdict
# and 10,000 patient-level resamples across the whole grid would cost minutes
# for numbers that are read descriptively. Declared here so the asymmetry is a
# decision rather than something noticed in the output.
BOOT_ROLES <- c("PRIMARY", "CONTROL")

# --- THE READING RULE, operationalised ---------------------------------------
# The declaration fixes four readings in words. The thresholds below turn them
# into arithmetic, and they are fixed before any Q1 value exists.
#
#   "clears zero"        the 95% patient-level bootstrap CI excludes 0
#   "materially smaller" the CI on (rho_tumour - rho_normal) excludes 0
#
# ABSENT IN NORMAL   normal CI includes 0, tumour CI excludes 0, AND the
#                    difference CI excludes 0.
# WEAKER IN NORMAL   both CIs exclude 0, AND the difference CI excludes 0.
# PRESENT IN NORMAL  normal CI excludes 0, AND the difference CI INCLUDES 0.
# INDETERMINATE      normal CI includes 0 AND the difference CI includes 0.
#
# INDETERMINATE IS A FIFTH READING AND IT IS NOT IN THE DECLARATION. It is
# added here, before any value, because the declaration's own mitigation (c)
# says "a wide interval around zero is not evidence of absence" and at n = 113
# that case is likely. Calling it ABSENT would be exactly the error mitigation
# (c) warns against, so it gets its own name instead of being folded into the
# predicted outcome. THIS IS A GAP IN THE DECLARATION, CLOSED IN THE
# CONSERVATIVE DIRECTION, and the result note must say so.
#
# UNINTERPRETABLE    OVERRIDES EVERYTHING. The declared control's difference CI
#                    excludes 0 AND |diff_control| >= CONTROL_FRAC *
#                    |diff_config|. "Differs as much as apoptosis does", with a
#                    factor-of-two allowance.
CONTROL_FRAC <- 0.5
#
# ONE DISCLOSURE ABOUT THIS RULE, AND IT MATTERS. As first written, the
# override was unconditional: control CI excludes 0 AND the ratio >= 0.5. That
# is defective, and a dry run is what exposed it. If the configuration's
# difference is near zero - PRESENT IN NORMAL, the coupling unchanged between
# tissues - then the denominator is near zero and ANY control movement makes
# the ratio large. But a control that moves where the configuration does not is
# not evidence of composition; it is evidence that the configuration is
# SPECIFICALLY unchanged. The declaration's own words presuppose otherwise:
# "the apoptotic result is uninterpretable" requires an apoptotic result to
# invalidate.
#
# So the override is CONDITIONAL on the cell having claimed a change - ABSENT
# or WEAKER IN NORMAL. The fix is argued from the declaration's text, not from
# any value.
#
# FULL DISCLOSURE OF THE SEQUENCE, because the appearance matters more than the
# arithmetic here: the defect was noticed in a redirected dry run whose output
# I had already seen, and that dry run returned UNINTERPRETABLE. The fix was
# then checked against it and IS CELL-FOR-CELL NEUTRAL - six cells fired before
# and the same six fire after, because every cell with an inflated ratio
# already had a control interval crossing zero and so never fired. The verdict
# does not move. This is recorded so that nobody has to reconstruct whether a
# threshold was renegotiated after a disappointing result: the threshold
# CONTROL_FRAC was not touched, and the change that was made is verifiable as
# inert on this data by re-running with the condition removed.

# The verdict is taken over 3 comparators x 2 instruments = 6 cells, and all
# six must agree. Anything else is MIXED and the mixture is the result.
#
# EACH CELL IS READ ON BOTH DECLARED ADJUSTMENTS AND THEY MUST AGREE. The
# declaration says to report the comparison "with and without" the epithelial
# estimate and does not say which carries the reading, so requiring both
# rather than choosing one is the only option that cannot be settled after the
# fact. A cell whose two declared adjustments disagree is SPLIT BY ADJUSTMENT
# and no reading is taken from it.
N_CELLS <- 6L

READING_RULE <- paste0(
  "Per cell (comparator x instrument), on EACH declared adjustment (raw and ",
  "EPITHELIAL-adjusted), which must agree: ABSENT IN NORMAL iff the normal CI ",
  "includes 0, the tumour CI excludes 0 and the difference CI excludes 0; ",
  "WEAKER IN NORMAL iff both exclude 0 and the difference CI excludes 0; ",
  "PRESENT IN NORMAL iff the normal CI excludes 0 and the difference CI ",
  "includes 0; INDETERMINATE iff the normal CI and the difference CI both ",
  "include 0. UNINTERPRETABLE overrides all of it iff the Fe-S cytosolic ",
  "difference CI excludes 0 and is at least ", CONTROL_FRAC, " times the ",
  "configuration's difference in absolute size, in the same cell. A reading ",
  "is taken only if all ", N_CELLS, " cells agree; otherwise MIXED, and the ",
  "mixture is the result.")

message("\n0. the reading rule, fixed before any value exists:\n   ",
        paste(strwrap(READING_RULE, width = 74), collapse = "\n   "))
set.seed(SEED)

# =============================================================================
# 1. The joint snapshot, and nothing else
# =============================================================================
message("\n1. reading the joint snapshot (E37)")

for (p in c(PATH_JOINT_LINEAR, PATH_JOINT_VST, PATH_JOINT_SCORES)) {
  if (!file.exists(p)) {
    stop("the joint snapshot is absent: ", basename(p),
         ". Source scripts/E37_joint_snapshot_and_gate.R first - E36 cannot ",
         "start without it, and results/ does not survive a fresh clone.",
         call. = FALSE)
  }
}

jl <- readRDS(PATH_JOINT_LINEAR)   # LINEAR DESeq2-normalised. Never logged here
jv <- readRDS(PATH_JOINT_VST)      # VST, log scale
js <- readRDS(PATH_JOINT_SCORES)

# The scale fields are the guard against a swapped object, not decoration.
stopifnot(identical(jl$scale, "linear_deseq2_normalised"),
          identical(jv$scale, "vst_log"),
          identical(rownames(jl$mat), rownames(jv$mat)),
          identical(colnames(jl$mat), colnames(jv$mat)))

smap <- jv$sample_map
stopifnot(nrow(smap) == ncol(jv$mat), !anyDuplicated(smap$sample_id),
          identical(smap$sample_id, colnames(jv$mat)))
N_TUM <- sum(smap$tissue == "tumour"); N_NOR <- sum(smap$tissue == "normal")
message("   ", nrow(smap), " samples: ", N_TUM, " tumour, ", N_NOR, " normal")
message("   genes ", nrow(jv$mat), "; linear max ",
        signif(max(jl$mat), 3), ", VST max ", signif(max(jv$mat), 3),
        "  (the cheap check that the two are not swapped)")

sd_ <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
e14 <- readRDS(PATH_E14)

# =============================================================================
# 2. The quadrant - rebuilt, because E37 did not persist it, and RE-ASSERTED
# =============================================================================
# E37's object holds the gate and its audit, nothing sample-level beyond
# sample_map. So the frozen constructor is called again, the same way, and the
# two declared assertions become a precondition of THIS run rather than
# something inherited from a previous script's say-so.
#
# NORMALS ARE NEVER QUADRANTED AND CANNOT BE: STATE takes BUFFER_gistic, which
# is copy number, and a normal sample has none.
message("\n2. the quadrant (tumours only; the frozen constructor, re-applied)")

st <- readRDS(PATH_STATE)
if (!identical(deparse(st$build_state), st$definition_source$build_state) ||
    !identical(environment(st$build_state), baseenv())) {
  stop("the frozen STATE constructor fails its own integrity contract",
       call. = FALSE)
}
.build_state <- st$build_state
LV <- st$spec$levels
S0 <- st$tcga$state

TUM <- smap$sample_id[smap$tissue == "tumour"]
tum_pat <- smap$patient[match(TUM, smap$sample_id)]
si <- match(tum_pat, S0$patient)
if (anyNA(si)) stop("a tumour has no frozen STATE row", call. = FALSE)

QD <- tibble::tibble(
  sample_id  = TUM,
  patient    = tum_pat,
  MYC        = as.numeric(js$gsva[MYC_REF, TUM]),
  OX_gsva    = as.numeric(js$gsva[ARM_OX, TUM]),
  OX_mitopps = as.numeric(js$mitopps_arms[ARM_OX, TUM]),
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
  a3 <- identical(which(q == QLEV[3]), which(state_new == LV[2]))
  a4 <- identical(which(q == QLEV[4]), which(state_new %in% LV[3:4]))
  if (!a3 || !a4) {
    stop("ASSERTION FAILED (", inst, "): Q3 == level 2 is ", a3,
         ", Q4 == levels 3 + 4 is ", a4, ". STOP.", call. = FALSE)
  }
  ASSERTS[[inst]] <- tibble::tibble(instrument = inst, Q3_equals_level2 = a3,
                                    Q4_equals_levels3_4 = a4)
  QD[[paste0("quad_", inst)]] <- factor(q, levels = QLEV)
  message("   ", inst, ": Q3 == level 2 TRUE; Q4 == levels 3 + 4 TRUE")
}
ASSERTS <- dplyr::bind_rows(ASSERTS)
message("   quadranted tumours ", sum(!is.na(QD$quad_gsva)), " of ", nrow(QD),
        " (", sum(is.na(QD$quad_gsva)), " have no BUFFER_gistic)")

# =============================================================================
# 3. The composites - log2(LINEAR + 1), z across the JOINT 1,208
# =============================================================================
# SCALE: every composite below is built from the LINEAR DESeq2-normalised
# matrix, logged here and nowhere else. The z is taken over all 1,208 samples
# so the two tissues sit on one scale. Spearman is invariant to that choice, so
# it changes NO correlation - it matters only for the descriptive means, where
# it is stated.
message("\n3. the composites (log2(linear + 1), z over the joint 1,208)")

GL_ALL <- log2(jl$mat + 1)

# Mean per-gene z over a set, with every gene entering +1 unless signs given.
.composite <- function(genes, signs = NULL, label) {
  miss <- setdiff(genes, rownames(GL_ALL))
  if (length(miss)) {
    stop("composite '", label, "' is missing gene(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  G <- GL_ALL[genes, , drop = FALSE]
  gv <- apply(G, 1L, stats::var)
  if (any(!(gv > 0))) {
    stop("composite '", label, "' has a zero-variance gene: ",
         paste(genes[!(gv > 0)], collapse = ", "), call. = FALSE)
  }
  Z <- t(scale(t(G)))
  w <- if (is.null(signs)) stats::setNames(rep(1, length(genes)), genes) else
    signs[rownames(Z)]
  as.numeric(colMeans(sweep(Z, 1L, w[rownames(Z)], "*")))
}

# --- the control's two halves, split by MitoCarta membership as in E14 -------
mito_inventory <- unique(readxl::read_excel(
  PATH_MITOCARTA, sheet = "A Human MitoCarta3.0")$Symbol)
.halves <- function(set_name) {
  g <- intersect(e14$comparators[[set_name]], rownames(GL_ALL))
  list(mito = intersect(g, mito_inventory),
       cyto = setdiff(g, mito_inventory))
}
FES <- .halves("Fe-S cluster assembly")
MPG <- .halves("mitophagy")
MPP <- .halves("mitophagy, PINK1/PRKN only")

# The control must not be partly made of the exposure, nor of the readout.
# Checked, not assumed - data note 2 section 5 records both as zero.
for (nm in c("Fe-S cluster assembly", "mitophagy",
             "mitophagy, PINK1/PRKN only")) {
  g <- intersect(e14$comparators[[nm]], rownames(GL_ALL))
  ov_cfg <- intersect(g, names(CONFIG_SIGN))
  ov_ox  <- intersect(g, sd_$arm_sets[[ARM_OX]])
  if (length(ov_cfg) || length(ov_ox)) {
    stop("control '", nm, "' overlaps the readout (",
         paste(ov_cfg, collapse = ","), ") or the exposure (",
         paste(ov_ox, collapse = ","), ")", call. = FALSE)
  }
}
message("   Fe-S cluster assembly: ", length(FES$mito), " mitochondrial / ",
        length(FES$cyto), " cytosolic; zero overlap with the configuration ",
        "and with the OXPHOS arm")
message("   mitophagy ", length(MPG$mito), "/", length(MPG$cyto),
        "; PINK1/PRKN only ", length(MPP$mito), "/", length(MPP$cyto))

COMP <- tibble::tibble(
  sample_id  = smap$sample_id,
  patient    = smap$patient,
  tissue     = factor(smap$tissue, levels = c("normal", "tumour")),
  TSS        = substr(smap$patient, 6L, 7L),
  config     = .composite(names(CONFIG_SIGN), CONFIG_SIGN, "configuration"),
  trigger    = .composite(TRIGGER, NULL, "trigger arm"),
  guardian   = .composite(GUARDIAN, c(BCL2L1 = 1, MCL1 = -1), "guardian arm"),
  fes_cyto   = .composite(FES$cyto, NULL, "Fe-S cytosolic"),
  fes_mito   = .composite(FES$mito, NULL, "Fe-S mitochondrial"),
  mitophagy_cyto    = .composite(MPG$cyto, NULL, "mitophagy cytosolic"),
  mitophagy_mito    = .composite(MPG$mito, NULL, "mitophagy mitochondrial"),
  mitophagy_pp_cyto = .composite(MPP$cyto, NULL, "PINK1/PRKN cytosolic"),
  mitophagy_pp_mito = .composite(MPP$mito, NULL, "PINK1/PRKN mitochondrial"),
  EPITHELIAL = .composite(EPI_GENES, NULL, "EPITHELIAL"),
  ADIPOSE    = .composite(ADIPOSE_GENES, NULL, "adipose companion"),
  FIBROBLAST = .composite(FIBROBLAST_GENES, NULL, "fibroblast companion"),
  OX_gsva    = as.numeric(js$gsva[ARM_OX, smap$sample_id]),
  OX_mitopps = as.numeric(js$mitopps_arms[ARM_OX, smap$sample_id]),
  PROLIF     = as.numeric(js$gsva[PROLIF_COV, smap$sample_id]),
  M_a        = as.numeric(js$gsva[MYC_REF, smap$sample_id]),
  PROLIFSTRIP = as.numeric(js$gsva["FELSHER__PROLIFSTRIP", smap$sample_id]))

# --- M_b, all four CollecTRI variants ----------------------------------------
# SCALE: ULM reads the VST matrix. And ULM IS NOT COHORT-RELATIVE - a score
# depends only on its own sample's expression across genes, demonstrated in
# data note 2 section 1 (identical over 1,208 and over the 1,095 tumours, max
# deviation 0, r = 1). So M_b carries absolute information that M_a, being
# GSVA, cannot. It is still NOT comparable with M_b in data/from_validation/,
# because the VST and the size factors were computed over the joint set.
message("   M_b: four CollecTRI variants by ULM on the VST matrix")
ct <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE, progress = FALSE)
# SIGN RULE, identical to E02 section 3.1b and to the validation study's
# script 06 at d3ac60e: mor +1 for a stimulatory edge, -1 for an inhibitory
# one, an edge flagged BOTH takes +1. Reproduced, not improved on.
myc_net <- ct %>%
  dplyr::filter(source_genesymbol == "MYC", !is.na(target_genesymbol),
                target_genesymbol != "") %>%
  dplyr::transmute(source = "MYC", target = target_genesymbol,
                   mor = dplyr::if_else(as.logical(is_stimulation), 1, -1),
                   likelihood = 1) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)

.ulm <- function(targets, E, label) {
  net <- myc_net %>%
    dplyr::filter(target %in% targets, target %in% rownames(E))
  message("      ", sprintf("%-20s %3d targets in the matrix", label,
                            nrow(net)))
  res <- decoupleR::run_ulm(mat = E, network = net, .source = source,
                            .target = target, .mor = mor, minsize = 5L)
  v <- res %>% dplyr::filter(statistic == "ulm", source == "MYC") %>%
    dplyr::select(condition, score)
  stats::setNames(v$score, v$condition)[colnames(E)]
}
for (nm in MB_VARIANTS) {
  COMP[[nm]] <- as.numeric(.ulm(sd_$collectri_sets[[nm]], jv$mat, nm))
}

# The quadrant joins on sample_id; normals stay NA, by construction.
COMP <- COMP %>%
  dplyr::left_join(QD %>% dplyr::select(sample_id, quad_gsva, quad_mitopps),
                   by = "sample_id")
stopifnot(all(is.na(COMP$quad_gsva[COMP$tissue == "normal"])))
message("   ", nrow(COMP), " rows; normals carry no quadrant, as required")

# =============================================================================
# 4. The three comparators, and the PATIENT-LEVEL bootstrap
# =============================================================================
# Amendment section 3: "resample the patients with replacement, carry each
# drawn patient's tumour and, where it exists, that patient's normal TOGETHER
# into the resample, and recompute both correlations and their difference
# inside each one." The resampling unit is the PATIENT, never the sample,
# because the same 113 patients contribute to both correlations.
message("\n4. the three comparators and the patient-level resamples")

NOR_PAT   <- COMP$patient[COMP$tissue == "normal"]
NOR_SITES <- sort(unique(substr(NOR_PAT, 6L, 7L)))
tum_rows  <- which(COMP$tissue == "tumour")
nor_rows  <- which(COMP$tissue == "normal")

# patient -> row index, one per tissue. Unique by construction (one sample per
# patient per type), asserted rather than trusted.
stopifnot(!anyDuplicated(COMP$patient[tum_rows]),
          !anyDuplicated(COMP$patient[nor_rows]))
T_IX <- stats::setNames(tum_rows, COMP$patient[tum_rows])
N_IX <- stats::setNames(nor_rows, COMP$patient[nor_rows])

POOL <- list(
  all    = names(T_IX),
  site   = names(T_IX)[substr(names(T_IX), 6L, 7L) %in% NOR_SITES],
  paired = intersect(names(T_IX), NOR_PAT))
stopifnot(identical(sort(POOL$paired), sort(NOR_PAT)),
          all(POOL$paired %in% POOL$site), all(POOL$site %in% POOL$all))

COMPARATOR_NOTE <- c(
  all    = "all tumours; the basis E34 and E35 used; confounded with site",
  site   = paste0("tumours at the ", length(NOR_SITES),
                  " sites that contributed a normal (",
                  paste(NOR_SITES, collapse = ", "), ")"),
  paired = "the patient-matched tumours; identical patients, hence sites")

for (cm in COMPARATOR_LEV) {
  message("   ", sprintf("%-7s n_tumour %5d  n_normal %4d   %s", cm,
          length(POOL[[cm]]),
          sum(POOL[[cm]] %in% names(N_IX)), COMPARATOR_NOTE[[cm]]))
}

# One resample set per comparator, drawn once and reused across every readout,
# adjustment and instrument, so paired differences are taken on the SAME
# resamples without refitting.
message("   drawing ", format(N_BOOT, big.mark = ","),
        " patient-level resamples per comparator")
BOOT <- lapply(COMPARATOR_LEV, function(cm) {
  P  <- POOL[[cm]]
  np <- length(P)
  lapply(seq_len(N_BOOT), function(b) {
    d  <- P[sample.int(np, np, replace = TRUE)]
    nn <- N_IX[d]
    list(t = unname(T_IX[d]), n = unname(nn[!is.na(nn)]))
  })
})
names(BOOT) <- COMPARATOR_LEV
message("   mean normal rows per resample: ",
        paste(sprintf("%s %.0f", COMPARATOR_LEV,
                      vapply(COMPARATOR_LEV, function(cm)
                        mean(vapply(BOOT[[cm]], function(z) length(z$n),
                                    integer(1))), numeric(1))),
              collapse = " | "))

# =============================================================================
# 5. Q1 - the coupling, within each tissue, then contrasted
# =============================================================================
# E32's and E33's estimator, verbatim: rank, residualise on the ranked
# covariates, correlate the residuals. Written out so the bootstrap re-ranks
# INSIDE every resample rather than ranking once outside it.
message("\n5. Q1: partial Spearman of OXPHOS against each readout")

.partial_spearman <- function(x, y, Z = NULL) {
  ok <- stats::complete.cases(x, y, Z)
  if (sum(ok) < 10L) return(NA_real_)
  rx <- rank(x[ok]); ry <- rank(y[ok])
  if (is.null(Z)) return(stats::cor(rx, ry))
  RZ <- apply(as.matrix(Z)[ok, , drop = FALSE], 2L, rank)
  H  <- qr(cbind(1, RZ))
  stats::cor(qr.resid(H, rx), qr.resid(H, ry))
}
.ci <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, c((1 - CI_LEVEL) / 2, 1 - (1 - CI_LEVEL) / 2)))
}
.excl0 <- function(lo, hi) isTRUE(lo > 0) || isTRUE(hi < 0)

.block <- function(cm, inst, readout, adj_name, do_boot) {
  ox_col <- INSTR[[inst]]
  zc     <- ADJUST[[adj_name]]
  pt     <- POOL[[cm]]
  rt     <- unname(T_IX[pt])
  rn     <- unname(N_IX[intersect(pt, names(N_IX))])
  .ps <- function(rows) .partial_spearman(
    COMP[[ox_col]][rows], COMP[[readout]][rows],
    if (length(zc)) as.matrix(COMP[rows, zc, drop = FALSE]) else NULL)
  est_t <- .ps(rt); est_n <- .ps(rn)
  if (do_boot) {
    bs <- BOOT[[cm]]
    bt <- vapply(bs, function(z) .ps(z$t), numeric(1))
    bn <- vapply(bs, function(z) .ps(z$n), numeric(1))
    bd <- bt - bn
    ci_t <- .ci(bt); ci_n <- .ci(bn); ci_d <- .ci(bd)
  } else {
    ci_t <- ci_n <- ci_d <- c(NA_real_, NA_real_)
  }
  tibble::tibble(
    comparator = cm, instrument = inst, readout = readout,
    adjustment = adj_name,
    adjusted_for = if (length(zc)) paste(zc, collapse = " + ") else "-",
    n_tumour = sum(stats::complete.cases(
      COMP[rt, c(ox_col, readout, zc), drop = FALSE])),
    n_normal = sum(stats::complete.cases(
      COMP[rn, c(ox_col, readout, zc), drop = FALSE])),
    rho_normal = est_n, n_lo = ci_n[1], n_hi = ci_n[2],
    rho_tumour = est_t, t_lo = ci_t[1], t_hi = ci_t[2],
    diff = est_t - est_n, d_lo = ci_d[1], d_hi = ci_d[2],
    normal_excludes_0 = .excl0(ci_n[1], ci_n[2]),
    tumour_excludes_0 = .excl0(ci_t[1], ci_t[2]),
    diff_excludes_0   = .excl0(ci_d[1], ci_d[2]),
    bootstrapped = do_boot)
}

GRID <- tidyr::expand_grid(
  comparator = COMPARATOR_LEV, instrument = names(INSTR),
  readout = READOUTS$readout, adjustment = names(ADJUST)) %>%
  dplyr::left_join(READOUTS, by = "readout") %>%
  dplyr::mutate(do_boot = role %in% BOOT_ROLES)
message("   ", nrow(GRID), " cells, of which ", sum(GRID$do_boot),
        " bootstrapped (PRIMARY and the declared CONTROL only)")

Q1 <- vector("list", nrow(GRID))
for (i in seq_len(nrow(GRID))) {
  g <- GRID[i, ]
  if (g$do_boot) {
    message("      ", sprintf("%-7s %-8s %-18s %-11s", g$comparator,
                              g$instrument, g$readout, g$adjustment))
  }
  Q1[[i]] <- .block(g$comparator, g$instrument, g$readout, g$adjustment,
                    g$do_boot)
}
Q1 <- dplyr::bind_rows(Q1) %>%
  dplyr::left_join(READOUTS, by = "readout") %>%
  dplyr::mutate(
    comparator = factor(comparator, levels = COMPARATOR_LEV),
    declared_adjustment = adjustment %in% ADJUST_DECLARED)

# =============================================================================
# 6. THE READING - applied mechanically, from the rule fixed in section 0
# =============================================================================
message("\n6. the reading")

# The COMPANION readouts, printed because they are declared-reported and
# because the LEVEL of correlation they show in normals bounds what any
# normal-side number means. Unsigned composites in a tissue whose variance is
# dominated by one composition axis all ride that axis together.
.companion_range <- function(tis) {
  v <- Q1 %>%
    dplyr::filter(role == "COMPANION", adjustment == "raw") %>%
    dplyr::pull(!!rlang::sym(paste0("rho_", tis)))
  c(min(v, na.rm = TRUE), max(v, na.rm = TRUE))
}
CMP_RANGE <- list(normal = .companion_range("normal"),
                  tumour = .companion_range("tumour"))
message("   the five COMPANION composites against OXPHOS, unadjusted:")
message("     in normals ", sprintf("%+.3f to %+.3f", CMP_RANGE$normal[1],
                                    CMP_RANGE$normal[2]),
        " | in tumours ", sprintf("%+.3f to %+.3f", CMP_RANGE$tumour[1],
                                  CMP_RANGE$tumour[2]))

.read_cell <- function(r) {
  if (is.na(r$rho_normal) || is.na(r$diff)) return(NA_character_)
  if (!r$normal_excludes_0 && r$tumour_excludes_0 && r$diff_excludes_0)
    return("ABSENT IN NORMAL")
  if (r$normal_excludes_0 && r$tumour_excludes_0 && r$diff_excludes_0)
    return("WEAKER IN NORMAL")
  if (r$normal_excludes_0 && !r$diff_excludes_0)
    return("PRESENT IN NORMAL")
  if (!r$normal_excludes_0 && !r$diff_excludes_0)
    return("INDETERMINATE")
  # The remaining combination: normal clears zero, tumour does not, and the
  # difference does. Not anticipated by the declaration; named, never folded.
  "COUPLING IN NORMAL ONLY"
}

CFG  <- Q1 %>% dplyr::filter(readout == "config", declared_adjustment)
CTRL <- Q1 %>% dplyr::filter(readout == "fes_cyto", declared_adjustment)

READ <- CFG %>%
  dplyr::rowwise() %>%
  dplyr::mutate(reading = .read_cell(dplyr::pick(dplyr::everything()))) %>%
  dplyr::ungroup() %>%
  dplyr::select(comparator, instrument, adjustment, rho_normal, n_lo, n_hi,
                rho_tumour, t_lo, t_hi, diff, d_lo, d_hi, reading)

# --- the control's override, matched cell by cell ---------------------------
# CONDITIONAL ON A CLAIMED CHANGE - see the disclosure in section 0.
CTRL_CHK <- CTRL %>%
  dplyr::select(comparator, instrument, adjustment,
                c_diff = diff, c_lo = d_lo, c_hi = d_hi,
                c_excl = diff_excludes_0) %>%
  dplyr::left_join(CFG %>% dplyr::select(comparator, instrument, adjustment,
                                         cfg_diff = diff),
                   by = c("comparator", "instrument", "adjustment")) %>%
  dplyr::left_join(READ %>% dplyr::select(comparator, instrument, adjustment,
                                          reading),
                   by = c("comparator", "instrument", "adjustment")) %>%
  dplyr::mutate(
    ratio = abs(c_diff) / abs(cfg_diff),
    claimed_change = reading %in% c("ABSENT IN NORMAL", "WEAKER IN NORMAL"),
    control_fails  = c_excl & (ratio >= CONTROL_FRAC) & claimed_change)

UNINTERPRETABLE <- any(CTRL_CHK$control_fails, na.rm = TRUE)

# Reported whatever the override does: a control that moves where the
# configuration does NOT is not a composition problem, it is the opposite, and
# it is named rather than discarded.
CTRL_ASYMMETRY <- CTRL_CHK %>%
  dplyr::filter(c_excl, !claimed_change) %>%
  dplyr::select(comparator, instrument, adjustment, c_diff, cfg_diff, reading)

# --- per cell, both declared adjustments must agree -------------------------
CELL <- READ %>%
  dplyr::group_by(comparator, instrument) %>%
  dplyr::summarise(
    readings = paste(sort(unique(reading)), collapse = " / "),
    agreed   = dplyr::n_distinct(reading) == 1L,
    reading  = dplyr::if_else(agreed, dplyr::first(reading),
                              "SPLIT BY ADJUSTMENT"),
    .groups = "drop")

if (UNINTERPRETABLE) {
  VERDICT <- "UNINTERPRETABLE"
  VERDICT_TEXT <- paste0(
    "UNINTERPRETABLE. The declared positive control - the Fe-S cluster ",
    "assembly CYTOSOLIC half - changes between normal and tumour by at least ",
    CONTROL_FRAC, " of the configuration's change, with a difference interval ",
    "excluding zero, in ", sum(CTRL_CHK$control_fails, na.rm = TRUE), " of ",
    nrow(CTRL_CHK), " cells. The analysis is measuring composition and the ",
    "apoptotic result is uninterpretable. This is the stop condition the ",
    "declaration set, and NO READING IS TAKEN.")
} else if (all(CELL$agreed) && dplyr::n_distinct(CELL$reading) == 1L &&
           nrow(CELL) == N_CELLS) {
  VERDICT <- unique(CELL$reading)
  VERDICT_TEXT <- paste0(
    VERDICT, ". All ", N_CELLS, " cells agree, on both declared adjustments, ",
    "and the declared control is clean.")
} else {
  VERDICT <- "MIXED"
  VERDICT_TEXT <- paste0(
    "MIXED. The ", N_CELLS, " cells do not all agree, so no single reading is ",
    "taken and THE MIXTURE IS THE RESULT: ",
    paste(sprintf("%s/%s %s", CELL$comparator, CELL$instrument, CELL$reading),
          collapse = "; "), ". The declared control is clean.")
}

message("\n   ", paste(strwrap(VERDICT_TEXT, width = 74),
                       collapse = "\n   "))

# =============================================================================
# 7. Q2 - where tumour MYC sits relative to normal breast. NO VERDICT.
# =============================================================================
# Amendment section 2: the WITHIN-PATIENT difference is primary; the
# distributional form the premise is stated in is kept and also reported.
message("\n7. Q2: MYC level, normal against tumour. DESCRIPTIVE - no verdict")

MYC_EST <- c("M_a", "PROLIFSTRIP", MB_VARIANTS)
EST_NOTE <- c(
  M_a              = "GSVA, cohort-relative: a RANK within a set that is 91% tumour",
  PROLIFSTRIP      = "GSVA, cohort-relative: the same caveat",
  M_b__FULL        = "ULM, NOT cohort-relative: independent of the mixture",
  M_b__MITOSTRIP   = "ULM, NOT cohort-relative; the construction counterpart of M_a",
  M_b__PROLIFSTRIP = "ULM, NOT cohort-relative",
  M_b__BOTHSTRIP   = "ULM, NOT cohort-relative")

.wilson <- function(k, n) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2); p <- k / n
  d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}

# --- 7.1 the paired difference, tumour minus that patient's OWN normal -------
PAIR <- tibble::tibble(patient = POOL$paired) %>%
  dplyr::mutate(row_t = unname(T_IX[patient]), row_n = unname(N_IX[patient]))
stopifnot(nrow(PAIR) == N_NOR, !anyNA(PAIR$row_t), !anyNA(PAIR$row_n))

Q2_PAIRED <- dplyr::bind_rows(lapply(MYC_EST, function(est) {
  d <- COMP[[est]][PAIR$row_t] - COMP[[est]][PAIR$row_n]
  k <- sum(d > 0); n <- length(d); w <- .wilson(k, n)
  tibble::tibble(
    estimator = est, instrument_free = TRUE, n_pairs = n,
    median_diff = stats::median(d),
    iqr_lo = unname(stats::quantile(d, 0.25)),
    iqr_hi = unname(stats::quantile(d, 0.75)),
    n_tumour_higher = k, pct_tumour_higher = 100 * k / n,
    pct_lo = 100 * w[1], pct_hi = 100 * w[2],
    note = EST_NOTE[[est]])
}))

message("   7.1 within-patient difference over ", N_NOR, " pairs:")
Q2_PAIRED %>%
  dplyr::transmute(estimator, n_pairs, median_diff = round(median_diff, 3),
                   IQR = paste0("[", round(iqr_lo, 2), ", ",
                                round(iqr_hi, 2), "]"),
                   pct_higher = round(pct_tumour_higher, 1),
                   CI = paste0("[", round(pct_lo, 1), ", ",
                               round(pct_hi, 1), "]")) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 7.2 the distributional form, by quadrant, on both unpaired comparators --
# The normal median is the threshold. For the GSVA estimators it is a position
# in a tumour-dominated pooled distribution, NOT an independent reference
# point - data note 2 section 1.1. For the ULM estimators it is not.
Q2_DIST <- dplyr::bind_rows(lapply(c("all", "site"), function(cm) {
  pt <- POOL[[cm]]
  dplyr::bind_rows(lapply(MYC_EST, function(est) {
    nmed <- stats::median(COMP[[est]][nor_rows])
    dplyr::bind_rows(lapply(names(INSTR), function(inst) {
      qc <- paste0("quad_", inst)
      dd <- COMP[COMP$patient %in% pt & COMP$tissue == "tumour", ]
      dplyr::bind_rows(
        tibble::tibble(comparator = cm, estimator = est, instrument = inst,
                       group = "normal", n = length(nor_rows),
                       median = nmed,
                       iqr_lo = unname(stats::quantile(COMP[[est]][nor_rows], .25)),
                       iqr_hi = unname(stats::quantile(COMP[[est]][nor_rows], .75)),
                       pct_above_normal_median = NA_real_),
        dd %>% dplyr::filter(!is.na(.data[[qc]])) %>%
          dplyr::group_by(group = as.character(.data[[qc]])) %>%
          dplyr::summarise(
            n = dplyr::n(),
            median = stats::median(.data[[est]]),
            iqr_lo = unname(stats::quantile(.data[[est]], .25)),
            iqr_hi = unname(stats::quantile(.data[[est]], .75)),
            pct_above_normal_median = 100 * mean(.data[[est]] > nmed),
            .groups = "drop") %>%
          dplyr::mutate(comparator = cm, estimator = est, instrument = inst))
    }))
  }))
})) %>% dplyr::select(comparator, estimator, instrument, group, n, median,
                      iqr_lo, iqr_hi, pct_above_normal_median)

message("\n   7.2 % of tumours above the normal median, by quadrant",
        " (comparator 'all', GSVA):")
Q2_DIST %>%
  dplyr::filter(comparator == "all", instrument == "gsva", group != "normal") %>%
  dplyr::transmute(estimator, group = substr(group, 1, 2), n,
                   pct = round(pct_above_normal_median, 1)) %>%
  tidyr::pivot_wider(names_from = group, values_from = c(n, pct)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 7.3 the paired difference BY QUADRANT - descriptive, n printed ----------
# Data note 2 section 2: this runs n = 17 to 39 and is DESCRIPTIVE. The counts
# are carried in the table and printed in the figure panel rather than implied.
Q2_PAIRED_QUAD <- dplyr::bind_rows(lapply(names(INSTR), function(inst) {
  qc <- paste0("quad_", inst)
  dplyr::bind_rows(lapply(MYC_EST, function(est) {
    tibble::tibble(
      instrument = inst, estimator = est,
      quad = as.character(COMP[[qc]][PAIR$row_t]),
      d = COMP[[est]][PAIR$row_t] - COMP[[est]][PAIR$row_n]) %>%
      dplyr::filter(!is.na(quad)) %>%
      dplyr::group_by(instrument, estimator, quad) %>%
      dplyr::summarise(n_pairs = dplyr::n(),
                       median_diff = stats::median(d),
                       pct_tumour_higher = 100 * mean(d > 0),
                       .groups = "drop")
  }))
}))
message("\n   7.3 by quadrant (DESCRIPTIVE, n = ",
        min(Q2_PAIRED_QUAD$n_pairs), " to ", max(Q2_PAIRED_QUAD$n_pairs), ")")

# =============================================================================
# 8. EPITHELIAL, as a composition estimate. Companion - it enters no model.
# =============================================================================
# Amendment section 4: reported ONLY to show that EPITHELIAL behaves as a
# composition estimate should - higher in tumours, anti-correlated with stroma.
message("\n8. does EPITHELIAL behave as a composition estimate?")

EPI_CHK <- tibble::tibble(
  score = c("EPITHELIAL", "ADIPOSE", "FIBROBLAST"),
  median_normal = vapply(c("EPITHELIAL", "ADIPOSE", "FIBROBLAST"),
                         function(s) stats::median(COMP[[s]][nor_rows]),
                         numeric(1)),
  median_tumour = vapply(c("EPITHELIAL", "ADIPOSE", "FIBROBLAST"),
                         function(s) stats::median(COMP[[s]][tum_rows]),
                         numeric(1)),
  auc_tumour_higher = vapply(c("EPITHELIAL", "ADIPOSE", "FIBROBLAST"),
    function(s) {
      r <- rank(COMP[[s]]); nt <- length(tum_rows); nn <- length(nor_rows)
      (sum(r[tum_rows]) - nt * (nt + 1) / 2) / (nt * nn)
    }, numeric(1)),
  rho_with_epithelial = vapply(c("EPITHELIAL", "ADIPOSE", "FIBROBLAST"),
    function(s) stats::cor(COMP$EPITHELIAL, COMP[[s]], method = "spearman"),
    numeric(1)))
EPI_CHK %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   AUC is P(a tumour scores higher than a normal); 0.5 is no",
        " separation.")

# The amendment's validity criterion, applied mechanically: EPITHELIAL should
# be "strongly higher in tumours, and strongly anti-correlated with the
# stromal score". Both are checked, and the answer is recorded whichever way
# it comes out, because that is what a declared validity check is for.
EPI_VALID <- list(
  auc = EPI_CHK$auc_tumour_higher[EPI_CHK$score == "EPITHELIAL"],
  rho_adipose = EPI_CHK$rho_with_epithelial[EPI_CHK$score == "ADIPOSE"],
  rho_fibroblast = EPI_CHK$rho_with_epithelial[EPI_CHK$score == "FIBROBLAST"])
EPI_VALID$higher_in_tumour <- EPI_VALID$auc > 0.5
EPI_VALID$anticorrelated_with_stroma <-
  EPI_VALID$rho_adipose < 0 && EPI_VALID$rho_fibroblast < 0
EPI_VALID$passes <- EPI_VALID$higher_in_tumour &&
  EPI_VALID$anticorrelated_with_stroma

# --- the diagnostic that explains the answer, whichever way it went ---------
# DIAGNOSTIC ONLY. These two halves ENTER NO MODEL and are not an alternative
# adjustment - substituting one after seeing EPITHELIAL behave is exactly the
# move the amendment's gene list exists to prevent. They are here to say WHY
# the composite does what it does, which a single AUC cannot.
LUMINAL <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19")
BASAL   <- c("KRT5", "KRT14", "KRT17", "TP63")
.auc <- function(v) {
  r <- rank(v); nt <- length(tum_rows); nn <- length(nor_rows)
  (sum(r[tum_rows]) - nt * (nt + 1) / 2) / (nt * nn)
}
EPI_HALVES <- tibble::tibble(
  half = c("luminal (5)", "basal / myoepithelial (4)", "both (9), DECLARED"),
  auc_tumour_higher = c(.auc(.composite(LUMINAL, NULL, "luminal")),
                        .auc(.composite(BASAL, NULL, "basal")),
                        .auc(COMP$EPITHELIAL)))
EPI_GENE_AUC <- tibble::tibble(
  gene = EPI_GENES,
  compartment = dplyr::if_else(EPI_GENES %in% LUMINAL, "luminal", "basal"),
  auc_tumour_higher = vapply(EPI_GENES,
    function(g) .auc(as.numeric(GL_ALL[g, ])), numeric(1)))

message("\n   the amendment's validity criterion for EPITHELIAL: ",
        if (EPI_VALID$passes) "PASSES" else "FAILS")
message("     higher in tumours: ", EPI_VALID$higher_in_tumour,
        " (AUC ", round(EPI_VALID$auc, 3), ")")
message("     anti-correlated with both stromal scores: ",
        EPI_VALID$anticorrelated_with_stroma,
        " (adipose ", round(EPI_VALID$rho_adipose, 3),
        ", fibroblast ", round(EPI_VALID$rho_fibroblast, 3), ")")
message("   DIAGNOSTIC, entering no model - the two halves separately:")
EPI_HALVES %>%
  dplyr::mutate(auc_tumour_higher = round(auc_tumour_higher, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
EPI_GENE_AUC %>%
  dplyr::mutate(auc_tumour_higher = round(auc_tumour_higher, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
if (!EPI_VALID$passes) {
  message("   -> THE DECLARED ADJUSTMENT FAILS ITS OWN VALIDITY CHECK. The ",
          "EPITHELIAL-adjusted\n      columns may NOT be read as ",
          "'composition-adjusted'. The gene list is NOT changed - it was ",
          "fixed in the amendment at be64f9a -\n      and no substitute is ",
          "introduced. The failure is reported as a limit on what the ",
          "adjusted columns mean.")
}

# =============================================================================
# 8.1 The qualifications that travel with the verdict, whatever it is
# =============================================================================
# Assembled here because the EPITHELIAL check runs after the reading. Each one
# is attached to VERDICT_TEXT so that it CANNOT be separated from the verdict
# when either is quoted.
QUALIFICATIONS <- character(0)

if (!EPI_VALID$passes) {
  QUALIFICATIONS <- c(QUALIFICATIONS, paste0(
    "THE DECLARED ADJUSTMENT FAILS ITS OWN VALIDITY CHECK. The amendment ",
    "required EPITHELIAL to be strongly higher in tumours and strongly ",
    "anti-correlated with the stromal score. It is neither: AUC ",
    round(EPI_VALID$auc, 3), " (0.5 is no separation), rho with adipose ",
    round(EPI_VALID$rho_adipose, 3), " and with fibroblast ",
    round(EPI_VALID$rho_fibroblast, 3), ". The luminal and basal halves ",
    "separate the tissues in OPPOSITE directions and largely cancel - the ",
    "diagnostic in section 8 gives both. So the EPITHELIAL-adjusted columns ",
    "are NOT 'composition-adjusted'; they are adjusted for a score that does ",
    "not track composition. The gene list was NOT changed and no substitute ",
    "was introduced. Note also that the adjustment MOVES THE NORMAL ",
    "CORRELATION A LONG WAY, which is itself a reason not to read it as a ",
    "composition control."))
}
if (nrow(CTRL_ASYMMETRY)) {
  QUALIFICATIONS <- c(QUALIFICATIONS, paste0(
    "In ", nrow(CTRL_ASYMMETRY), " cell(s) the control's difference interval ",
    "excludes zero while the configuration's does not. That is NOT a ",
    "composition failure - it is the opposite - and it is reported rather ",
    "than discarded."))
}
QUALIFICATIONS <- c(QUALIFICATIONS, paste0(
  "The two stromal companions, which enter no model, separate the tissues ",
  "strongly (adipose AUC ",
  round(EPI_CHK$auc_tumour_higher[EPI_CHK$score == "ADIPOSE"], 3),
  ", fibroblast ",
  round(EPI_CHK$auc_tumour_higher[EPI_CHK$score == "FIBROBLAST"], 3),
  "). So composition DOES differ enormously between the two tissues; it is ",
  "the declared estimate of it that fails to capture the difference."))
if (all(CELL$reading == "SPLIT BY ADJUSTMENT")) {
  QUALIFICATIONS <- c(QUALIFICATIONS, paste0(
    "NO READING WOULD HAVE BEEN TAKEN EVEN IF THE CONTROL HAD BEEN CLEAN. All ",
    N_CELLS, " cells are SPLIT BY ADJUSTMENT - the unadjusted and the ",
    "EPITHELIAL-adjusted readings disagree in every one of them - so the rule ",
    "in section 0 yields no reading by a route entirely independent of the ",
    "stop condition. The verdict therefore does not rest on the control ",
    "override or on its threshold."))
}
QUALIFICATIONS <- c(QUALIFICATIONS, paste0(
  "THE LEVEL OF CORRELATION IN NORMAL TISSUE BOUNDS EVERY NORMAL-SIDE NUMBER. ",
  "Unadjusted, the five COMPANION composites correlate with OXPHOS at ",
  sprintf("%+.3f to %+.3f", CMP_RANGE$normal[1], CMP_RANGE$normal[2]),
  " in normals against ",
  sprintf("%+.3f to %+.3f", CMP_RANGE$tumour[1], CMP_RANGE$tumour[2]),
  " in tumours. In normal breast almost everything correlates with almost ",
  "everything, which is what a bulk tissue looks like when its variance is ",
  "dominated by a single composition axis. Against that background the ",
  "configuration's own normal-side value is NOT a low correlation in a quiet ",
  "tissue; it is a near-zero value in a tissue where every unsigned composite ",
  "sits far above it. Whether the signed construction of the configuration (it ",
  "carries -MCL1) is what cancels that common component is NOT established ",
  "here and is a question for a separate analysis."))
if (UNINTERPRETABLE) {
  QUALIFICATIONS <- c(QUALIFICATIONS, paste0(
    "THE STOP CONDITION BINDS Q1 ONLY. Q2 does not depend on the control, on ",
    "the configuration or on any adjustment - it is a within-patient ",
    "difference in MYC activity - so Q2 STANDS and is reported in full."))
}

VERDICT_TEXT <- paste0(VERDICT_TEXT, "\n\nQUALIFICATIONS THAT TRAVEL WITH ",
                       "THIS VERDICT:\n",
                       paste0("  (", seq_along(QUALIFICATIONS), ") ",
                              QUALIFICATIONS, collapse = "\n"))

message("\n8.1 qualifications attached to the verdict: ",
        length(QUALIFICATIONS))
for (q in QUALIFICATIONS) {
  message("   - ", paste(strwrap(q, width = 72), collapse = "\n     "))
}

# =============================================================================
# 9. The figure - three panels, and every n printed
# =============================================================================
message("\n9. the figure")

.ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
.ensure_dir(DIR_FIGURES); .ensure_dir(DIR_TABLES); .ensure_dir(DIR_DOC_FIGURES)

CMP_LAB <- c(all = paste0("all tumours (n = ", length(POOL$all), ")"),
             site = paste0("site-matched (n = ", length(POOL$site), ")"),
             paired = paste0("patient-matched (n = ", length(POOL$paired), ")"))
ADJ_LAB <- c(raw = "unadjusted", epithelial = "+ EPITHELIAL",
             epi_prolif = "+ EPITHELIAL + PROLIF (companion)")

# Subtitles are wrapped explicitly. An unwrapped ggplot subtitle is CLIPPED at
# the device edge rather than reflowed, which is how E34 lost the right-hand
# end of its caption; the fix there was strwrap and it is the fix here.
.sub <- function(x, width = 118) paste(strwrap(x, width = width),
                                       collapse = "\n")

base_theme <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
    strip.text = ggplot2::element_text(size = 8),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold", size = 10),
    plot.subtitle = ggplot2::element_text(size = 8),
    plot.title.position = "plot")

# --- Panel A: Q1, the two within-tissue correlations, side by side -----------
pA_dat <- CFG %>%
  tidyr::pivot_longer(
    cols = c(rho_normal, rho_tumour),
    names_to = "tissue", values_to = "rho") %>%
  dplyr::mutate(
    lo = dplyr::if_else(tissue == "rho_normal", n_lo, t_lo),
    hi = dplyr::if_else(tissue == "rho_normal", n_hi, t_hi),
    tissue = factor(dplyr::if_else(tissue == "rho_normal",
                                   "normal (n = 113)", "tumour"),
                    levels = c("normal (n = 113)", "tumour")),
    cmp = factor(CMP_LAB[as.character(comparator)], levels = CMP_LAB),
    adj = factor(ADJ_LAB[adjustment], levels = ADJ_LAB))

pA <- ggplot2::ggplot(pA_dat,
    ggplot2::aes(x = rho, y = cmp, colour = tissue)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = lo, xmax = hi),
    position = ggplot2::position_dodge(width = 0.6), linewidth = 0.5) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6),
                      size = 1.9) +
  ggplot2::facet_grid(adj ~ instrument) +
  ggplot2::scale_colour_manual(values = c("normal (n = 113)" = "#1f78b4",
                                          "tumour" = "#c2272d"), name = NULL) +
  ggplot2::labs(
    title = "A. Q1: the OXPHOS-to-configuration coupling, within each tissue",
    subtitle = .sub(paste0("Partial Spearman, ", format(N_BOOT, big.mark = ","),
                      " PATIENT-level bootstrap resamples. All three tumour ",
                      "comparators, always. Normals carry no quadrant and ",
                      "none is used here.")),
    x = "partial Spearman (OXPHOS vs the configuration composite)", y = NULL) +
  base_theme

# --- Panel B: the difference, and the declared control beside it -------------
pB_dat <- Q1 %>%
  dplyr::filter(declared_adjustment,
                readout %in% c("config", "fes_cyto", "fes_mito")) %>%
  dplyr::mutate(
    ro = factor(c(config = "configuration (PRIMARY)",
                  fes_cyto = "Fe-S cytosolic (CONTROL)",
                  fes_mito = "Fe-S mitochondrial (not a control)")[readout],
                levels = c("configuration (PRIMARY)",
                           "Fe-S cytosolic (CONTROL)",
                           "Fe-S mitochondrial (not a control)")),
    cmp = factor(CMP_LAB[as.character(comparator)], levels = CMP_LAB),
    adj = factor(ADJ_LAB[adjustment], levels = ADJ_LAB))

pB <- ggplot2::ggplot(pB_dat, ggplot2::aes(x = diff, y = cmp, colour = ro)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = d_lo, xmax = d_hi),
    position = ggplot2::position_dodge(width = 0.65), linewidth = 0.5,
    na.rm = TRUE) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.65),
                      size = 1.9) +
  ggplot2::facet_grid(adj ~ instrument) +
  ggplot2::scale_colour_manual(
    values = c("configuration (PRIMARY)" = "#1b1b1b",
               "Fe-S cytosolic (CONTROL)" = "#e08214",
               "Fe-S mitochondrial (not a control)" = "grey60"),
    name = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3)) +
  ggplot2::labs(
    title = "B. tumour minus normal, with the declared positive control",
    subtitle = .sub(paste0("The control is PRE-SPECIFIED BUT NOT BLIND. ",
                      "Fe-S mitochondrial is NOT a control - it is ",
                      "mitochondrial, so coupling to OXPHOS is near-guaranteed ",
                      "- and carries no interval, being a companion.")),
    x = "difference in partial Spearman (tumour - normal)", y = NULL) +
  base_theme

# --- Panel C: Q2, the within-patient difference ------------------------------
pC_dat <- Q2_PAIRED %>%
  dplyr::mutate(
    est = factor(estimator, levels = rev(MYC_EST)),
    lab = sprintf("%.0f%% [%.0f, %.0f]", pct_tumour_higher, pct_lo, pct_hi),
    kind = dplyr::if_else(grepl("^M_b", estimator),
                          "ULM (not cohort-relative)",
                          "GSVA (cohort-relative)"))

pC <- ggplot2::ggplot(pC_dat,
    ggplot2::aes(x = median_diff, y = est, colour = kind)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = iqr_lo, xmax = iqr_hi),
                          linewidth = 0.5) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(ggplot2::aes(label = lab), hjust = 0, nudge_y = 0.26,
                     size = 2.5, show.legend = FALSE) +
  ggplot2::scale_colour_manual(
    values = c("GSVA (cohort-relative)" = "#6a51a3",
               "ULM (not cohort-relative)" = "#238b45"), name = NULL) +
  ggplot2::labs(
    title = paste0("C. Q2: tumour MYC minus that patient's OWN normal, ",
                   N_NOR, " pairs. DESCRIPTIVE - no verdict"),
    subtitle = .sub(paste0("Point is the median paired difference, bar the ",
                      "IQR. The label is the % of pairs in which the tumour ",
                      "exceeds its own normal, with a Wilson interval. ",
                      "Immune to the collection-site skew by construction, ",
                      "and unaffected by the Q1 stop condition.")),
    x = "median within-patient difference (tumour - own normal)", y = NULL) +
  base_theme

CAPTION <- paste(strwrap(paste0(
  "E36. EXPLORATORY, nothing pre-registered. TCGA-BRCA only; SCAN-B carries ",
  "no normals, so this is UNREPLICATED. Joint snapshot built by E37 over ",
  nrow(smap), " samples (", N_TUM, " tumour, ", N_NOR, " normal); no value ",
  "here is comparable with one from data/from_validation/. The normals are ",
  "FIELD-ADJACENT tissue from breasts that grew a tumour, and all ", N_NOR,
  " come from patients whose tumour is in the same analysis - the design is ",
  "PAIRED. They come from ", length(NOR_SITES), " of 40 collection sites, ",
  "which is why every across-person comparison is shown on three nested ",
  "tumour sets and why the within-patient contrast in panel C is the one ",
  "immune to it. Cross-sectional, bulk, transcript-level: N3, nothing here ",
  "is 'primed', and the configuration is not lethality. No product term; no ",
  "outcome variable. VERDICT: ", VERDICT, "."), width = 164), collapse = "\n")

fig <- patchwork::wrap_plots(pA, pB, pC, ncol = 1,
                             heights = c(1.15, 1.15, 0.95)) +
  patchwork::plot_annotation(caption = CAPTION,
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(size = 5.6, hjust = 0,
                                           lineheight = 1.15),
      plot.caption.position = "plot"))

for (p in c(PATH_E36_FIG, PATH_E36_FIG_DOC)) {
  ggplot2::ggsave(p, fig, width = 9.2, height = 12.6, dpi = 300)
}
for (p in c(PATH_E36_FIG_PDF, PATH_E36_FIG_DOC_PDF)) {
  ggplot2::ggsave(p, fig, width = 9.2, height = 12.6,
                  device = grDevices::cairo_pdf)
}
message("   ", PATH_E36_FIG, "\n   ", PATH_E36_FIG_DOC, "  (tracked)")
message("   ", PATH_E36_FIG_PDF, "\n   ", PATH_E36_FIG_DOC_PDF, "  (tracked)")

# =============================================================================
# 10. Save
# =============================================================================
message("\n10. saving")

readr::write_csv(Q1, PATH_E36_CSV_Q1)
readr::write_csv(Q2_DIST, PATH_E36_CSV_Q2)

saveRDS(list(
  verdict        = VERDICT,
  verdict_text   = VERDICT_TEXT,
  reading_rule   = READING_RULE,
  control_frac   = CONTROL_FRAC,
  q1             = Q1,
  reading        = READ,
  cell           = CELL,
  control_check  = CTRL_CHK,
  control_asymmetry = CTRL_ASYMMETRY,
  uninterpretable = UNINTERPRETABLE,
  qualifications = QUALIFICATIONS,
  epi_valid      = EPI_VALID,
  epi_halves     = EPI_HALVES,
  epi_gene_auc   = EPI_GENE_AUC,
  companion_range = CMP_RANGE,
  q2_paired      = Q2_PAIRED,
  q2_paired_quad = Q2_PAIRED_QUAD,
  q2_dist        = Q2_DIST,
  epi_check      = EPI_CHK,
  asserts        = ASSERTS,
  pools          = lapply(POOL, length),
  pool_patients  = POOL,
  normal_sites   = NOR_SITES,
  sets = list(
    config = CONFIG_SIGN, trigger = TRIGGER, guardian = GUARDIAN,
    epithelial = EPI_GENES, adipose = ADIPOSE_GENES,
    fibroblast = FIBROBLAST_GENES,
    fes_mito = FES$mito, fes_cyto = FES$cyto,
    mitophagy_mito = MPG$mito, mitophagy_cyto = MPG$cyto,
    mitophagy_pp_mito = MPP$mito, mitophagy_pp_cyto = MPP$cyto),
  spec = list(
    snapshot = "JOINT tumour+normal, built by E37. NEVER data/from_validation/",
    n_boot = N_BOOT, ci_level = CI_LEVEL, seed = SEED,
    instruments = INSTR, adjustments = ADJUST,
    declared_adjustments = ADJUST_DECLARED,
    boot_roles = BOOT_ROLES, readouts = READOUTS,
    comparator_note = COMPARATOR_NOTE,
    myc_estimators = MYC_EST, estimator_note = EST_NOTE,
    declaration = "docs/2026-10-01_E36_declaration.md at d71fa45",
    amendment   = "docs/2026-10-01_E36_declaration_amendment.md at be64f9a",
    data_note_1 = "docs/2026-10-01_e36_data.md at 4802508",
    data_note_2 = "docs/2026-10-01_e36_data_part2.md",
    control_primary = "Fe-S cluster assembly, CYTOSOLIC half (declared 4802508)",
    control_not_blind = paste0(
      "In tumours E14 recorded the Fe-S cytosolic half against OXPHOS at ",
      "+0.152 raw / +0.118 prolif-adjusted (TCGA) and +0.214 / +0.159 ",
      "(SCAN-B). The control asks whether that CHANGES between tissues as ",
      "much as the configuration's does, not whether it is zero."),
    limits = paste0(
      "Field-adjacent normal tissue, n = ", N_NOR, ", unreplicated (SCAN-B ",
      "has no normals), cross-sectional, bulk, transcript-level. The ",
      "configuration is not lethality. Site-matching cannot remove site, ",
      "because the normals are concentrated inside the matched set too."),
    indeterminate_is_new = paste0(
      "INDETERMINATE is a FIFTH reading, not in the declaration. Added in the ",
      "script before any value, because the declaration's mitigation (c) says ",
      "a wide interval around zero is not evidence of absence. A gap in the ",
      "declaration, closed in the conservative direction.")),
  built = Sys.time()), PATH_E36)

message("   ", PATH_E36, "\n   ", PATH_E36_CSV_Q1, "\n   ", PATH_E36_CSV_Q2)
message("\nE36 done. VERDICT: ", VERDICT, "\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E36)

  cat(x$verdict, "\n\n"); cat(x$verdict_text, "\n")

  # The reading, cell by cell, on both declared adjustments.
  x$reading %>%
    dplyr::transmute(comparator, instrument, adjustment,
                     normal = round(rho_normal, 3),
                     n_ci = paste0("[", round(n_lo, 2), ",", round(n_hi, 2), "]"),
                     tumour = round(rho_tumour, 3),
                     t_ci = paste0("[", round(t_lo, 2), ",", round(t_hi, 2), "]"),
                     diff = round(diff, 3),
                     d_ci = paste0("[", round(d_lo, 2), ",", round(d_hi, 2), "]"),
                     reading) %>%
    as.data.frame()

  # Did the six cells agree?
  x$cell %>% as.data.frame()

  # THE STOP CONDITION. control_fails TRUE anywhere means UNINTERPRETABLE.
  x$control_check %>%
    dplyr::transmute(comparator, instrument, adjustment,
                     control_diff = round(c_diff, 3),
                     c_ci = paste0("[", round(c_lo, 2), ",", round(c_hi, 2), "]"),
                     config_diff = round(cfg_diff, 3),
                     ratio = round(ratio, 2), control_fails) %>%
    as.data.frame()

  # The companions, point estimates only - never a verdict.
  x$q1 %>%
    dplyr::filter(role == "COMPANION", adjustment == "epithelial") %>%
    dplyr::transmute(comparator, instrument, readout,
                     normal = round(rho_normal, 3),
                     tumour = round(rho_tumour, 3),
                     diff = round(diff, 3)) %>%
    as.data.frame()

  # Q2, the primary form and the distributional one.
  x$q2_paired %>% as.data.frame()
  x$q2_dist %>% dplyr::filter(comparator == "all", estimator == "M_a") %>%
    as.data.frame()
  x$q2_paired_quad %>% dplyr::filter(estimator == "M_a") %>% as.data.frame()

  # Does EPITHELIAL behave as a composition estimate? (It is a DECLARED check,
  # and its answer is a limit on what the adjusted columns mean.)
  x$epi_check %>% as.data.frame()
  utils::str(x$epi_valid)
  x$epi_halves %>% as.data.frame()
  x$epi_gene_auc %>% as.data.frame()

  # The qualifications are part of the verdict, not a footnote to it.
  for (q in x$qualifications) cat("-", q, "\n\n")
  x$control_asymmetry %>% as.data.frame()

  # The declared assertions, and the three pools.
  x$asserts %>% as.data.frame()
  utils::str(x$pools)
  x$normal_sites
}
