# E34_quadrant_configuration.R
# =============================================================================
# THE APOPTOTIC CONFIGURATION ACROSS THE MYC x OXPHOS QUADRANTS
#
# EXPLORATORY. Nothing here is pre-registered.
#
# Declaration: docs/2026-09-30_E34_declaration.md, COMMITTED ALONE AT 688539d
# before any data file was opened. Data note, with every decision the
# declaration left open: docs/2026-09-30_e34_data.md at 1345ec5. Read both
# before this file; this header restates their rules, it does not invent them.
#
# LARGELY CONFIRMATORY AND PRESENTATIONAL. E11 showed, on two continuous axes,
# that the apoptotic machinery is ordered by OXPHOS and not by MYC. This puts
# that finding into four groups a reader can see, and isolates the one cell the
# continuous treatment cannot: Q2, MYC-low and OXPHOS-high, which STATE
# collapses into its level 1. If the verdict is OXPHOS-ORDERED it CONFIRMS E11
# IN A GROUP-LEVEL FORM and is written as confirmation, NEVER as discovery.
#
# THE CIRCULARITY, from the data note 1.2 and repeated here because it bounds
# every reading: the configuration's genes AND their signs were read off OXPHOS
# correlations IN THESE SAME TWO COHORTS. A composite built that way separates
# OXPHOS-high from OXPHOS-low largely BY CONSTRUCTION. So the two OXPHOS
# contrasts are close to guaranteed and confirm almost nothing.
# THE INFORMATIVE CONTRASTS ARE THE TWO MYC CONTRASTS.
#
# =============================================================================
# THE VARIABLE
# =============================================================================
#   Q1 MYC-low  OXPHOS-low      Q2 MYC-low  OXPHOS-high
#   Q3 MYC-high OXPHOS-low      Q4 MYC-high OXPHOS-high
#
# From the FROZEN constructor's own calls: STATE level 1 split by the OXPHOS
# call it already computes, levels 3 and 4 merged. E33's construction, COPIED
# VERBATIM - E33 is on an unmerged branch and is never sourced, so E34 stands
# whether or not E33 is merged.
#   TCGA    the stored frozen STATE, reproduced from its own stored inputs. 55
#           patients have no quadrant because that object takes its complete
#           cases from GISTIC BUFFER, which the quadrant never uses. Kept, not
#           repaired, so TCGA is identical to E33 and to the frozen object.
#   SCAN-B  no stored STATE exists. The same constructor is called on SCAN-B's
#           own scores with a CONSTANT FALSE buffer - the quadrant does not use
#           BUFFER, so a constant contributes no NA and level 4 is empty. The
#           constructor takes its own within-cohort median; there is no argument
#           through which a TCGA threshold could enter.
# ASSERTED IN BOTH COHORTS ON BOTH INSTRUMENTS, AND THE SCRIPT STOPS ON FAILURE:
# Q3 == STATE level 2 exactly, Q4 == levels 3 + 4 exactly.
#
# =============================================================================
# READOUTS - data note 1.1. E11 DEFINES NO COMPOSITE; this is E33's, verbatim
# =============================================================================
#   PRIMARY   configuration composite = mean z of BBC3, BID, BIK, BAD, BCL2L1
#             and MINUS MCL1, on log2(linear + 1), z across the whole cohort
#   ALWAYS    trigger arm (BBC3, BID, BIK, BAD); guardian arm (BCL2L1, -MCL1)
#   REPORTED  the twelve transcripts, UNSIGNED. Never the verdict. No FDR.
#
# =============================================================================
# THE FOUR CONTRASTS - high minus low. NO PRODUCT TERM ANYWHERE.
# =============================================================================
#   Q4 - Q2   MYC at OXPHOS-high   THE PRIMARY   predicted: none, or much smaller
#   Q3 - Q1   MYC at OXPHOS-low    control       predicted: none
#   Q2 - Q1   OXPHOS at MYC-low                  predicted: positive
#   Q4 - Q3   OXPHOS at MYC-high                 predicted: positive
# The MYC x OXPHOS product-term estimand failed in Block C, Block B, Block G, H4
# and N2. It is not fitted here. Four group means stay four group means.
#
# =============================================================================
# MODELS - all three reported, never chosen between
# =============================================================================
#   m1  readout ~ quadrant + PAM50 + purity + leukocyte fraction
#   m2  m1 without PAM50, fitted WITHIN each PAM50 subtype - the result
#   m3  m1 + PROLIF_DISJOINT - mandatory companion
# SCAN-B HAS NO PURITY OR LEUKOCYTE ESTIMATE (CLAUDE.md trap 2), so there m1 is
# quadrant + PAM50 and m2 is quadrant alone. Forced by the data, and said.
# Blom normal scores, pooled ONCE PER COHORT; every fit subsets them and NEVER
# re-scores (E23's rule). Group estimates are standardised means from lm.
#
# =============================================================================
# THE READING RULE - declaration section 4, made operational in data note 5
# =============================================================================
# Applied in every READING CELL: 2 cohorts x 2 instruments x {m1, m3, m2 in
# each READABLE subtype, readable iff all four quadrant cells hold n >= 20}.
# ALL READING CELLS MUST AGREE (the author's decision, 2026-09-30).
#   OXPHOS-ORDERED  Q2-Q1 > 0 and Q4-Q3 > 0, both CIs excluding zero, AND
#                   |Q4-Q2| and |Q3-Q1| each smaller than both OXPHOS contrasts
#   MYC-DEPENDENT   Q4-Q2 > 0, CI excluding zero, AND Q4-Q2 at least as large as
#                   the smaller of the two OXPHOS contrasts
#   MIXED           anything else, INCLUDING ANY DISAGREEMENT between cohorts,
#                   instruments or models
# THE VERDICT READS THE CONFIGURATION COMPOSITE ONLY.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - NO outcome or survival variable enters. NO causal language.
#   - NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM: no cohort available
#     carries an OXPHOS-directed or BH3-directed intervention.
#   - Nothing is written to myc_human_validation or myc_mouse.
#   - No per-gene FDR. No product term. No value compared across cohorts.
#
# N3: transcript scores. Nothing is "primed".
# SPECIES: human throughout. NO ORTHOLOG FUNCTION IS CALLED.
# SCALE: genes on log2(linear DESeq2-normalised + 1). The quadrant's inputs are
# READ, never re-scored: GSVA on VST for M_a and ox_gsva, mitoPPS on linear for
# ox_ppd. GSVA is cohort-relative and nothing is pooled across cohorts.
# SCAN-B symbols resolve through scanb_pheno.rds$symbol_map, as CLAUDE.md
# requires, and the E02 map is asserted to agree.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE34: the configuration across the MYC x OXPHOS quadrants\n",
        strrep("=", 78))

PATH_E34         <- file.path(DIR_RESULTS, "e34_quadrant_configuration.rds")
PATH_E34_CSV     <- file.path(DIR_TABLES,  "E34_contrasts.csv")
PATH_E34_FIG     <- file.path(DIR_FIGURES, "E34_configuration_by_quadrant.png")
PATH_E34_FIG_PDF <- file.path(DIR_FIGURES, "E34_configuration_by_quadrant.pdf")
# Written twice, as from E31 on: outputs/ is gitignored, so the note plus the
# tracked figure are the durable record.
DIR_DOC_FIGURES      <- here::here("docs", "figures")
PATH_E34_FIG_DOC     <- file.path(DIR_DOC_FIGURES,
                                  "2026-10-01_E34_configuration_by_quadrant.png")
PATH_E34_FIG_DOC_PDF <- file.path(DIR_DOC_FIGURES,
                                  "2026-10-01_E34_configuration_by_quadrant.pdf")

# =============================================================================
# 0. Constants, and the reading rule transcribed
# =============================================================================
# Read-only, from the FROZEN repo. Nothing is written there, ever.
PATH_STATE <- "/Users/gs/code/myc_human_validation/results/state_definition.rds"

ARM_OX      <- "OXPHOS subunits"    # the narrow set, NOT the umbrella
PROLIF_COV  <- "PROLIF_DISJOINT"
MIN_CELL_M2 <- 20L                  # the author's decision, data note 5

# E33's arms and signs, verbatim. The signs come from the reconciliation doc
# section 3.1 - "sensitisers and BCL2L1 up, MCL1 down" in OXPHOS-high tumours.
TRIGGER  <- c("BBC3", "BID", "BIK", "BAD")
GUARDIAN <- c("BCL2L1", "MCL1")
CONFIG_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)
stopifnot(setequal(names(CONFIG_SIGN), c(TRIGGER, GUARDIAN)))
# E10's twelve, verbatim (E10 line 117). Reported UNSIGNED; never the verdict.
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
TWELVE <- c(PRIMING_PRO, PRIMING_ANTI)
stopifnot(all(names(CONFIG_SIGN) %in% TWELVE), !anyDuplicated(TWELVE))

READ_PRIMARY <- c(config_comp = "configuration composite")
READ_ARMS    <- c(arm_trigger = "trigger arm", arm_guardian = "guardian arm")
READ_GENES   <- stats::setNames(TWELVE, paste0("g_", TWELVE))
READOUTS     <- c(READ_PRIMARY, READ_ARMS, READ_GENES)

QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",
          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")
INSTR        <- c(gsva = "OX_gsva", mitopps = "OX_mitopps")
PAM50_LEVELS <- c("LumA", "LumB", "HER2", "Basal", "Normal")

# The four contrasts, high minus low. `kind` is what the reading rule reads.
CONTR <- tibble::tribble(
  ~contrast,  ~hi,     ~lo,     ~kind,    ~role,
  "Q4 - Q2",  QLEV[4], QLEV[2], "MYC",    "PRIMARY: MYC at OXPHOS-high",
  "Q3 - Q1",  QLEV[3], QLEV[1], "MYC",    "control: MYC at OXPHOS-low",
  "Q2 - Q1",  QLEV[2], QLEV[1], "OXPHOS", "OXPHOS at MYC-low",
  "Q4 - Q3",  QLEV[4], QLEV[3], "OXPHOS", "OXPHOS at MYC-high")

V_OX    <- "OXPHOS-ORDERED"
V_MYC   <- "MYC-DEPENDENT"
V_MIXED <- "MIXED"
READING_RULE <- paste0(
  "In every reading cell (2 cohorts x 2 instruments x {m1, m3, m2 in each ",
  "readable subtype, all four cells n >= ", MIN_CELL_M2, "}): OXPHOS-ORDERED ",
  "iff Q2-Q1 > 0 and Q4-Q3 > 0 with CIs excluding zero AND |Q4-Q2| and |Q3-Q1| ",
  "each smaller than both OXPHOS contrasts; MYC-DEPENDENT iff Q4-Q2 > 0 with ",
  "CI excluding zero AND Q4-Q2 >= the smaller OXPHOS contrast. OXPHOS-ORDERED ",
  "or MYC-DEPENDENT only if EVERY reading cell agrees; otherwise MIXED. The ",
  "verdict reads the configuration composite only.")

message("\n0. the rule, fixed at 688539d / 1345ec5 before any readout existed:\n   ",
        READING_RULE)

# =============================================================================
# 1. The frozen constructor, read-only, and its integrity contract
# =============================================================================
message("\n1. the frozen STATE constructor")
st <- readRDS(PATH_STATE)
if (!identical(deparse(st$build_state), st$definition_source$build_state) ||
    !identical(environment(st$build_state), baseenv())) {
  stop("the frozen STATE constructor fails its own integrity contract",
       call. = FALSE)
}
.build_state <- st$build_state
LV <- st$spec$levels
if (!identical(LV, levels(st$tcga$state$STATE_gsva))) {
  stop("STATE levels moved", call. = FALSE)
}
message("   integrity contract holds; the constructor's environment is baseenv()")

# =============================================================================
# 2. Inputs
# =============================================================================
message("\n2. inputs")

frames <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
mito   <- readRDS(PATH_TCGA_MITO)
nw     <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sc     <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
sph    <- readRDS(PATH_SCANB_PHENO)
ID_T   <- colnames(mito$gsva_arms)
ID_S   <- colnames(sc$gsva_arms)
stopifnot(length(ID_T) == EXPECT_TCGA_SAMPLES,
          length(ID_S) == EXPECT_SCANB_SAMPLES)

# --- 2.1 the twelve genes, on log2(linear DESeq2-normalised + 1) --------------
# SCAN-B resolves through scanb_pheno.rds$symbol_map (CLAUDE.md). The E02 map
# inside scanb_scores.rds is asserted to resolve the same twelve to the same
# rows, so the mandated map is used and the one the rest of the repo uses is
# checked against it rather than assumed equivalent (data note section 3).
.log2_genes <- function(path, ids, map) {
  lin <- readRDS(path)
  if (!identical(lin$scale, "linear_deseq2_normalised")) {
    stop(basename(path), " is not on the linear DESeq2-normalised scale",
         call. = FALSE)
  }
  G  <- log2(lin$mat[, ids, drop = FALSE] + 1)
  gr <- .gene_rows(TWELVE, G, .symbol_resolver(rownames(G), map))
  if (length(gr$missing)) {
    stop(basename(path), ": did not resolve ",
         paste(gr$missing, collapse = ", "), call. = FALSE)
  }
  list(mat = gr$mat[TWELVE, , drop = FALSE], rownames = rownames(lin$mat))
}
gt <- .log2_genes(PATH_TCGA_LINEAR, ID_T, NULL)
gs <- .log2_genes(PATH_SCANB_LINEAR, ID_S, sph$symbol_map)
.resolved <- function(rn, map) vapply(TWELVE, function(g) {
  r <- .symbol_resolver(rn, map)(g)
  if (length(r) == 1L) r else NA_character_
}, character(1))
if (!identical(.resolved(gs$rownames, sph$symbol_map),
               .resolved(gs$rownames, sc$symbol_map))) {
  stop("the pheno and E02 symbol maps resolve the twelve genes differently",
       call. = FALSE)
}
message("   twelve genes resolved in both cohorts; SCAN-B through the pheno ",
        "map,\n   identical under the E02 map")

# --- 2.2 E33's composites, verbatim. z across the WHOLE cohort ----------------
.composites <- function(GM) {
  M  <- GM[names(CONFIG_SIGN), , drop = FALSE]
  gv <- apply(M, 1L, stats::var)
  if (any(!(gv > 0))) {
    stop("zero-variance gene(s): ",
         paste(names(gv)[!(gv > 0)], collapse = ", "), call. = FALSE)
  }
  ZG <- t(scale(t(M)))
  .signed_mean <- function(genes) as.numeric(
    colMeans(sweep(ZG[genes, , drop = FALSE], 1L, CONFIG_SIGN[genes], "*")))
  out <- list(config_comp  = .signed_mean(names(CONFIG_SIGN)),
              arm_trigger  = .signed_mean(TRIGGER),
              arm_guardian = .signed_mean(GUARDIAN))
  for (g in TWELVE) out[[paste0("g_", g)]] <- as.numeric(GM[g, ])
  out
}
RAW <- list(TCGA = .composites(gt$mat), `SCAN-B` = .composites(gs$mat))
rm(gt, gs); invisible(gc(verbose = FALSE))
message("   composite = mean z(BBC3, BID, BIK, BAD, BCL2L1, -MCL1); ",
        "arms as E33;\n   the twelve are UNSIGNED")

# --- 2.3 covariates, from frames.rds - E11's and E19's source -----------------
.frame_for <- function(coh, ids) {
  f <- frames[frames$cohort == coh, ]
  f <- f[match(ids, f$sample_id), ]
  if (anyNA(f$sample_id)) stop(coh, ": ids missing from frames.rds", call. = FALSE)
  f
}
FT <- .frame_for("TCGA", ID_T)
FS <- .frame_for("SCAN-B", ID_S)

# =============================================================================
# 3. THE QUADRANT - E33's construction, verbatim, wrapped per cohort
# =============================================================================
message("\n3. the quadrant")

.build_quadrant <- function(cohort, S, buffer, stored = NULL) {
  states <- list()
  for (inst in names(INSTR)) {
    rebuilt <- .build_state(S$MYC, S[[INSTR[[inst]]]], buffer)
    if (!is.null(stored) && !identical(rebuilt, stored[[inst]])) {
      stop(cohort, ": the constructor does not reproduce STATE_", inst,
           " from its own stored inputs", call. = FALSE)
    }
    states[[inst]] <- rebuilt
  }
  myc_high <- states$gsva != LV[1]        # the constructor's own MYC call
  # The MYC call cannot depend on the OXPHOS instrument. Checked, not assumed.
  if (!identical(myc_high, states$mitopps != LV[1])) {
    stop(cohort, ": the two instruments disagree on the MYC call", call. = FALSE)
  }
  Q <- tibble::tibble(sample_id = S$sample_id, myc_high = myc_high)
  asserts <- list()
  for (inst in names(INSTR)) {
    ox           <- S[[INSTR[[inst]]]]
    state_stored <- states[[inst]]
    # The constructor computes its OXPHOS call internally and does not return
    # it. Calling it again with myc := oxphos makes its MYC-high call BE its
    # OXPHOS-high call: same median, same `>` tie rule, same complete cases.
    # No threshold is transcribed into this repo, so none can drift.
    probe <- .build_state(ifelse(is.na(S$MYC), NA_real_, ox), ox, buffer)
    if (any(probe == LV[2], na.rm = TRUE)) {
      stop(cohort, ": level 2 appeared in the myc := oxphos probe; the call is ",
           "not what it is taken to be", call. = FALSE)
    }
    if (!identical(is.na(probe), is.na(state_stored))) {
      stop(cohort, ": the probe's complete-case set differs from STATE's",
           call. = FALSE)
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
      stop("ASSERTION FAILED (", cohort, ", ", inst, "): Q3 == level 2 is ", a3,
           ", Q4 == levels 3 + 4 is ", a4, ". STOP.", call. = FALSE)
    }
    asserts[[inst]] <- tibble::tibble(
      cohort = cohort, instrument = inst,
      Q3_equals_level2 = a3, Q4_equals_levels3_4 = a4,
      level3 = sum(state_stored == LV[3], na.rm = TRUE),
      level4 = sum(state_stored == LV[4], na.rm = TRUE))
    Q[[paste0("quad_", inst)]] <- factor(q, levels = QLEV)
    message("   ", cohort, " ", inst,
            ": Q3 == level 2 TRUE; Q4 == levels 3 + 4 TRUE")
  }
  list(Q = Q, asserts = dplyr::bind_rows(asserts))
}

# --- TCGA: the stored frozen STATE, reproduced from its own inputs ------------
ST <- st$tcga$state
if (!setequal(ST$patient, ID_T)) {
  stop("the frozen STATE's patients differ from this repo's TCGA ids",
       call. = FALSE)
}
ST$sample_id <- ST$patient
rho_myc_frozen <- stats::cor(ST$MYC,
                             as.numeric(nw$tcga_gsva_new[MYC_REF, ST$patient]),
                             method = "spearman")
rho_ox_frozen  <- stats::cor(ST$OX_gsva,
                             as.numeric(mito$gsva_arms[ARM_OX, ST$patient]),
                             method = "spearman")
message("   TCGA: frozen MYC vs this repo's M_a, Spearman ",
        round(rho_myc_frozen, 6), "; frozen OX_gsva vs this repo's, ",
        round(rho_ox_frozen, 6))
QT <- .build_quadrant("TCGA", ST, ST$BUFFER_gistic,
                      stored = list(gsva = ST$STATE_gsva,
                                    mitopps = ST$STATE_mitopps))

# --- SCAN-B: the same constructor, SCAN-B's own scores, a CONSTANT buffer -----
SS <- tibble::tibble(
  sample_id  = ID_S,
  MYC        = as.numeric(sc$gsva_new[MYC_REF, ID_S]),
  OX_gsva    = as.numeric(sc$gsva_arms[ARM_OX, ID_S]),
  OX_mitopps = as.numeric(sc$mitopps_arms[ARM_OX, ID_S]))
QS <- .build_quadrant("SCAN-B", SS, rep(FALSE, nrow(SS)))
if (any(QS$asserts$level4 != 0L)) {
  stop("SCAN-B level 4 is not empty under a constant FALSE buffer", call. = FALSE)
}
ASSERTS <- dplyr::bind_rows(QT$asserts, QS$asserts)

# --- 3.1 cell sizes, NA, and instrument agreement ----------------------------
.cells <- function(cohort, Q) dplyr::bind_rows(lapply(names(INSTR), function(inst) {
  v <- Q[[paste0("quad_", inst)]]
  tibble::tibble(cohort = cohort, instrument = inst, quadrant = QLEV,
                 n = as.integer(table(factor(v, levels = QLEV))),
                 n_no_quadrant = sum(is.na(v)))
}))
CELLS <- dplyr::bind_rows(.cells("TCGA", QT$Q), .cells("SCAN-B", QS$Q))
.agree <- function(cohort, Q) {
  ok <- !is.na(Q$quad_gsva) & !is.na(Q$quad_mitopps)
  tb <- table(gsva = Q$quad_gsva[ok], mitopps = Q$quad_mitopps[ok])
  list(summary = tibble::tibble(cohort = cohort, n = sum(ok),
                                agreement = sum(diag(tb)) / sum(tb)),
       table = tb)
}
AGREE <- list(TCGA = .agree("TCGA", QT$Q), `SCAN-B` = .agree("SCAN-B", QS$Q))
message("\n   quadrant cell sizes, full cohorts:")
CELLS %>% tidyr::pivot_wider(names_from = quadrant, values_from = n) %>%
  as.data.frame() %>% print(row.names = FALSE)
for (coh in names(AGREE)) {
  message("   ", coh, ": GSVA and mitoPPS agree on ",
          round(100 * AGREE[[coh]]$summary$agreement, 1), "% of ",
          AGREE[[coh]]$summary$n, " with both")
}

# =============================================================================
# 4. Blom normal scores, pooled once per cohort
# =============================================================================
# E22's scorer, verbatim, plus an NA-aware wrapper: TCGA purity and leukocyte
# fraction have gaps, so their present values are scored among themselves and
# NA stays NA. Scored ONCE PER COHORT over the whole cohort; every fit below,
# including the within-subtype fits, SUBSETS these and never re-scores.
message("\n4. Blom normal scores, once per cohort")

.blom <- function(x) {
  n <- length(x)
  stopifnot(!anyNA(x))
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.375) / (n + 0.25))
  as.numeric(scale(z))
}
.blom_na <- function(x) {
  out <- rep(NA_real_, length(x))
  ok  <- !is.na(x)
  if (any(ok)) out[ok] <- .blom(x[ok])
  out
}

.frame <- function(cohort, ids, Q, FR, raw, prolif) {
  qi <- match(ids, Q$sample_id)
  if (anyNA(qi)) stop(cohort, ": quadrant rows missing", call. = FALSE)
  # Everything is computed BEFORE the tibble: inside tibble() a new column
  # shadows an outer object of the same name, which is the trap that stopped
  # E32's first run. No column below reuses the name of anything it needs.
  dat <- tibble::tibble(
    sample_id    = ids,
    quad_gsva    = Q$quad_gsva[qi],
    quad_mitopps = Q$quad_mitopps[qi],
    PAM50        = factor(as.character(FR$PAM50), levels = PAM50_LEVELS),
    purity       = .blom_na(FR$purity),
    leuko        = .blom_na(FR$leuko),
    PROLIF       = .blom(prolif))
  for (r in names(READOUTS)) {
    dat[[paste0("raw_", r)]] <- raw[[r]]
    dat[[r]] <- .blom(raw[[r]])
  }
  dat
}
DAT <- list(
  TCGA = .frame("TCGA", ID_T, QT$Q, FT, RAW$TCGA,
                as.numeric(mito$gsva_cov[PROLIF_COV, ID_T])),
  `SCAN-B` = .frame("SCAN-B", ID_S, QS$Q, FS, RAW$`SCAN-B`,
                    as.numeric(sc$gsva_cov[PROLIF_COV, ID_S])))
if (!all(is.na(DAT$`SCAN-B`$purity)) || !all(is.na(DAT$`SCAN-B`$leuko))) {
  stop("SCAN-B carries a purity or leukocyte estimate; trap 2 says it does not",
       call. = FALSE)
}
message("   TCGA: PAM50 NA ", sum(is.na(DAT$TCGA$PAM50)),
        ", purity and leuko both present ",
        sum(!is.na(DAT$TCGA$purity) & !is.na(DAT$TCGA$leuko)),
        "\n   SCAN-B: PAM50 NA ", sum(is.na(DAT$`SCAN-B`$PAM50)),
        ", NO purity or leukocyte estimate (trap 2) - forced, not a choice")

# =============================================================================
# 5. The models - standardised group means. NO PRODUCT TERM.
# =============================================================================
message("\n5. models")

COV <- list(
  TCGA     = list(m1 = c("PAM50", "purity", "leuko"),
                  m2 = c("purity", "leuko"),
                  m3 = c("PAM50", "purity", "leuko", "PROLIF")),
  `SCAN-B` = list(m1 = "PAM50",
                  m2 = character(0),
                  m3 = c("PAM50", "PROLIF")))
MODEL_LABEL <- c(m1 = "m1 (PAM50-adjusted)",
                 m2 = "m2 (within subtype)",
                 m3 = "m3 (m1 + PROLIF_DISJOINT)")

# One fit. The group estimate for quadrant q is the model's prediction averaged
# over its OWN analysis rows with quadrant set to q - a standardised mean, and a
# linear combination a'b whose 95% CI is a'Va on the residual df. The four
# contrasts are differences of those combinations; with no product term they
# equal coefficient differences, but they are computed generally so the code
# says what it means.
.fit_groups <- function(dat, y, covars) {
  cols <- c(y, "quad", covars)
  dd <- as.data.frame(dat[stats::complete.cases(dat[, cols]), cols, drop = FALSE])
  names(dd)[1] <- "yv"
  cell_n <- as.integer(table(factor(dd$quad, levels = QLEV)))
  names(cell_n) <- QLEV
  empty <- list(ok = FALSE, n = nrow(dd), cell_n = cell_n, df = NA_real_,
                A = NULL, b = NULL, V = NULL)
  dd$quad <- droplevels(factor(dd$quad, levels = QLEV))
  for (cv in covars) if (is.factor(dd[[cv]])) dd[[cv]] <- droplevels(dd[[cv]])
  if (nlevels(dd$quad) < 2L) return(empty)
  f <- tryCatch(stats::lm(stats::reformulate(c("quad", covars), response = "yv"),
                          data = dd), error = function(e) NULL)
  if (is.null(f) || f$df.residual < 1L || anyNA(stats::coef(f))) return(empty)
  tt <- stats::delete.response(stats::terms(f))
  A <- vapply(QLEV, function(q) {
    if (!q %in% levels(dd$quad)) return(rep(NA_real_, length(stats::coef(f))))
    dq <- dd
    dq$quad <- factor(rep(q, nrow(dd)), levels = levels(dd$quad))
    colMeans(stats::model.matrix(tt, data = dq, xlev = f$xlevels,
                                 contrasts.arg = f$contrasts))
  }, numeric(length(stats::coef(f))))
  list(ok = TRUE, n = nrow(dd), cell_n = cell_n, df = f$df.residual,
       A = A, b = stats::coef(f), V = stats::vcov(f))
}
.lincomb <- function(fit, a) {
  if (!fit$ok || anyNA(a)) {
    return(c(est = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  est <- sum(a * fit$b)
  se  <- sqrt(drop(t(a) %*% fit$V %*% a))
  tq  <- stats::qt(0.975, fit$df)
  c(est = est, se = se, lo = est - tq * se, hi = est + tq * se)
}

# The label of ONE reading cell, from its four contrasts. Transcribed from the
# data note section 5 and not reinterpreted.
.read_cell <- function(ct) {
  g <- function(nm, col) ct[[col]][ct$contrast == nm]
  ox1 <- g("Q2 - Q1", "est"); ox2 <- g("Q4 - Q3", "est")
  mp  <- g("Q4 - Q2", "est"); mc  <- g("Q3 - Q1", "est")
  if (anyNA(c(ox1, ox2, mp, mc))) return("not estimable")
  ox_pos <- ox1 > 0 && g("Q2 - Q1", "lo") > 0 &&
            ox2 > 0 && g("Q4 - Q3", "lo") > 0
  ox_min <- min(ox1, ox2)
  if (ox_pos && abs(mp) < ox_min && abs(mc) < ox_min) return(V_OX)
  if (mp > 0 && g("Q4 - Q2", "lo") > 0 && mp >= ox_min) return(V_MYC)
  "neither"
}

GROUPS <- list(); CONTRASTS <- list(); FITS <- list(); k <- 0L
for (coh in names(DAT)) {
  for (inst in names(INSTR)) {
    d0 <- DAT[[coh]]
    d0$quad <- d0[[paste0("quad_", inst)]]
    SETS <- list(list(model = "m1", subtype = "all", rows = rep(TRUE, nrow(d0))),
                 list(model = "m3", subtype = "all", rows = rep(TRUE, nrow(d0))))
    for (s in PAM50_LEVELS) {
      SETS[[length(SETS) + 1L]] <- list(
        model = "m2", subtype = s,
        rows = !is.na(d0$PAM50) & d0$PAM50 == s)
    }
    for (S_ in SETS) {
      dsub <- d0[S_$rows, , drop = FALSE]
      for (r in names(READOUTS)) {
        fit <- .fit_groups(dsub, r, COV[[coh]][[S_$model]])
        G <- t(vapply(QLEV, function(q) {
          if (!fit$ok) return(c(est = NA_real_, se = NA_real_,
                                lo = NA_real_, hi = NA_real_))
          .lincomb(fit, fit$A[, q])
        }, numeric(4)))
        C <- t(vapply(seq_len(nrow(CONTR)), function(i) {
          if (!fit$ok) return(c(est = NA_real_, se = NA_real_,
                                lo = NA_real_, hi = NA_real_))
          .lincomb(fit, fit$A[, CONTR$hi[i]] - fit$A[, CONTR$lo[i]])
        }, numeric(4)))
        # Computed before the tibbles, and no column reuses an outer name.
        ct_lab <- .read_cell(tibble::tibble(contrast = CONTR$contrast,
                                            est = C[, "est"], lo = C[, "lo"],
                                            hi = C[, "hi"]))
        is_readable <- fit$ok &&
          (S_$model != "m2" || min(fit$cell_n) >= MIN_CELL_M2)
        key <- tibble::tibble(cohort = coh, instrument = inst,
                              readout = unname(READOUTS[[r]]), readout_col = r,
                              model = S_$model, subtype = S_$subtype)
        k <- k + 1L
        GROUPS[[k]] <- dplyr::bind_cols(
          key[rep(1L, length(QLEV)), ],
          tibble::tibble(quadrant = QLEV, n_cell = unname(fit$cell_n[QLEV]),
                         est = G[, "est"], se = G[, "se"],
                         lo = G[, "lo"], hi = G[, "hi"]))
        CONTRASTS[[k]] <- dplyr::bind_cols(
          key[rep(1L, nrow(CONTR)), ],
          tibble::tibble(contrast = CONTR$contrast, kind = CONTR$kind,
                         role = CONTR$role, est = C[, "est"], se = C[, "se"],
                         lo = C[, "lo"], hi = C[, "hi"],
                         ci_excludes_0 = (C[, "lo"] > 0) | (C[, "hi"] < 0)))
        FITS[[k]] <- dplyr::bind_cols(
          key, tibble::tibble(n = fit$n, df = fit$df, fitted = fit$ok,
                              min_cell = min(fit$cell_n), readable = is_readable,
                              label = ct_lab))
      }
    }
  }
}
GROUPS    <- dplyr::bind_rows(GROUPS)
CONTRASTS <- dplyr::bind_rows(CONTRASTS)
FITS      <- dplyr::bind_rows(FITS)
message("   ", nrow(FITS), " fits: 2 cohorts x 2 instruments x ",
        length(READOUTS), " readouts x (m1, m3, m2 in ", length(PAM50_LEVELS),
        " subtypes)")

# =============================================================================
# 6. THE VERDICT - on the configuration composite only
# =============================================================================
message("\n6. the verdict\n", strrep("-", 78))

.overall <- function(fits) {
  # If any cohort x instrument has no readable subtype, m2 cannot contribute
  # there and the verdict is MIXED for want of the declared result.
  no_m2 <- fits %>% dplyr::filter(model == "m2") %>%
    dplyr::group_by(cohort, instrument) %>%
    dplyr::summarise(none = sum(readable) == 0L, .groups = "drop")
  if (nrow(no_m2) == 0L || any(no_m2$none)) return(V_MIXED)
  cells <- fits %>% dplyr::filter(readable)
  if (!nrow(cells)) return(V_MIXED)
  if (all(cells$label == V_OX))  return(V_OX)
  if (all(cells$label == V_MYC)) return(V_MYC)
  V_MIXED
}
VERDICT_CELLS <- FITS %>% dplyr::filter(readout_col == "config_comp")
verdict <- .overall(VERDICT_CELLS)

# The same rule on every other readout. REPORTED, NEVER THE VERDICT.
READOUT_SUMMARY <- FITS %>%
  dplyr::group_by(readout, readout_col) %>%
  dplyr::group_modify(~ tibble::tibble(
    overall              = .overall(.x),
    cells_read           = sum(.x$readable),
    cells_oxphos_ordered = sum(.x$readable & .x$label == V_OX),
    cells_myc_dependent  = sum(.x$readable & .x$label == V_MYC))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(role = dplyr::case_when(
    readout_col == "config_comp"      ~ "PRIMARY - THE VERDICT",
    readout_col %in% names(READ_ARMS) ~ "arm - reported always, not the verdict",
    TRUE                              ~ "gene - reported, not the verdict")) %>%
  dplyr::arrange(match(readout_col, names(READOUTS)))

message("RULE: ", READING_RULE, "\n")
message("   the configuration composite, m1 and m3, pooled:")
CONTRASTS %>%
  dplyr::filter(readout_col == "config_comp", model %in% c("m1", "m3"),
                subtype == "all") %>%
  dplyr::transmute(cohort, instrument, model, contrast, kind,
                   est = sprintf("%+.3f", est),
                   ci = sprintf("[%+.3f, %+.3f]", lo, hi)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   every reading cell for the verdict (configuration composite):")
VERDICT_CELLS %>%
  dplyr::transmute(cohort, instrument, model, subtype, n, min_cell, readable,
                   label) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   VERDICT: ", verdict)

reading <- switch(
  verdict,
  "OXPHOS-ORDERED" = paste0(
    "OXPHOS-ORDERED - the PREDICTED outcome, and it CONFIRMS E11 IN A\n",
    "   GROUP-LEVEL FORM. It is NOT a new finding and must not be written as\n",
    "   one. The configuration is present in MYC-low, OXPHOS-high tumours to\n",
    "   the same degree as in MYC-high ones, so the human transcriptome gives\n",
    "   NO SUPPORT to the claim that the standing pro-apoptotic state is\n",
    "   MYC-restricted. That claim rests on the mouse, where MYC-dependence is\n",
    "   shown by intervention, and the Discussion must say so: a LIMITATION TO\n",
    "   STATE, NOT A CONTRADICTION - transcript abundance is not a functional\n",
    "   readout. And the two OXPHOS contrasts are near-guaranteed by\n",
    "   construction (the circularity, data note 1.2); the MYC contrasts are\n",
    "   what carried information."),
  "MYC-DEPENDENT" = paste0(
    "MYC-DEPENDENT - surprising given E11 and E33, and it holds only because\n",
    "   EVERY reading cell agrees, within-subtype and proliferation-adjusted\n",
    "   included. The human data would corroborate the MYC restriction."),
  paste0(
    "MIXED - the rule is not met in every reading cell. The cell table above\n",
    "   says where it fails. The rule is NOT relaxed after the fact, and no\n",
    "   single agreeing cell is promoted."))
message("\n", reading)

message("\n   the same rule on every readout (REPORTED, NOT THE VERDICT):")
READOUT_SUMMARY %>%
  dplyr::select(readout, role, overall, cells_read, cells_oxphos_ordered,
                cells_myc_dependent) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. The figure - a candidate main-figure panel
# =============================================================================
# m1 standardised means with 95% CIs on the GSVA quadrant. MYC status on x,
# OXPHOS status as colour, so each OXPHOS contrast reads as a within-pair
# difference and each MYC contrast reads across the pairs. NO CONNECTING LINES:
# a line between the pairs would draw a product term this analysis does not fit.
# Okabe-Ito blue and vermillion; worst adjacent CVD delta E 21.9, all checks pass.
message("\n7. figure")

FIG_READ <- c(config_comp  = "Configuration composite (primary)",
              arm_trigger  = "Trigger arm (BBC3, BID, BIK, BAD)",
              arm_guardian = "Guardian arm (BCL2L1 minus MCL1)")
n_m1 <- FITS %>%
  dplyr::filter(model == "m1", instrument == "gsva",
                readout_col == "config_comp")
COH_LAB <- stats::setNames(
  paste0(n_m1$cohort, " (n = ", format(n_m1$n, big.mark = ",", trim = TRUE), ")"),
  n_m1$cohort)

PF <- GROUPS %>%
  dplyr::filter(instrument == "gsva", model == "m1", subtype == "all",
                readout_col %in% names(FIG_READ)) %>%
  dplyr::mutate(
    MYC = factor(ifelse(quadrant %in% QLEV[1:2], "MYC low", "MYC high"),
                 levels = c("MYC low", "MYC high")),
    OXPHOS = factor(ifelse(quadrant %in% QLEV[c(2, 4)], "OXPHOS high",
                           "OXPHOS low"),
                    levels = c("OXPHOS low", "OXPHOS high")),
    panel = factor(unname(FIG_READ[readout_col]), levels = unname(FIG_READ)),
    coh   = factor(unname(COH_LAB[cohort]), levels = unname(COH_LAB)))
y_lo <- min(PF$lo, na.rm = TRUE) - 0.20
y_hi <- max(PF$hi, na.rm = TRUE) + 0.34
PF$y_n <- y_lo + 0.05

# The three contrasts the figure is about, printed into each panel. The full
# table is in the note; these are the ones the reading rule turns on.
ANN <- CONTRASTS %>%
  dplyr::filter(instrument == "gsva", model == "m1", subtype == "all",
                readout_col %in% names(FIG_READ)) %>%
  dplyr::group_by(cohort, readout_col) %>%
  dplyr::summarise(
    txt = paste0(
      "MYC | OXPHOS-high  ", sprintf("%+.2f [%+.2f, %+.2f]",
                                     est[contrast == "Q4 - Q2"],
                                     lo[contrast == "Q4 - Q2"],
                                     hi[contrast == "Q4 - Q2"]),
      "\nOXPHOS | MYC-low   ", sprintf("%+.2f [%+.2f, %+.2f]",
                                       est[contrast == "Q2 - Q1"],
                                       lo[contrast == "Q2 - Q1"],
                                       hi[contrast == "Q2 - Q1"]),
      "\nOXPHOS | MYC-high  ", sprintf("%+.2f [%+.2f, %+.2f]",
                                       est[contrast == "Q4 - Q3"],
                                       lo[contrast == "Q4 - Q3"],
                                       hi[contrast == "Q4 - Q3"])),
    .groups = "drop") %>%
  dplyr::mutate(panel = factor(unname(FIG_READ[readout_col]),
                               levels = unname(FIG_READ)),
                coh   = factor(unname(COH_LAB[cohort]),
                               levels = unname(COH_LAB)))

DODGE <- ggplot2::position_dodge(width = 0.55)
OXCOL <- c("OXPHOS low" = "#0072B2", "OXPHOS high" = "#D55E00")
fig <- ggplot2::ggplot(PF, ggplot2::aes(x = MYC, y = est, colour = OXPHOS)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.16,
                         linewidth = 0.55, position = DODGE) +
  ggplot2::geom_point(size = 2.4, position = DODGE) +
  ggplot2::geom_text(ggplot2::aes(y = y_n, label = paste0("n = ", n_cell)),
                     position = DODGE, size = 2.1, colour = "grey35",
                     show.legend = FALSE) +
  ggplot2::geom_text(data = ANN,
                     ggplot2::aes(x = 0.42, y = y_hi, label = txt),
                     inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.0,
                     family = "mono", colour = "grey25", lineheight = 1.0) +
  ggplot2::facet_grid(panel ~ coh, switch = "y") +
  ggplot2::scale_colour_manual(values = OXCOL, name = NULL) +
  ggplot2::scale_y_continuous(limits = c(y_lo, y_hi),
                              expand = ggplot2::expansion(mult = 0)) +
  ggplot2::labs(
    x = NULL,
    y = "Adjusted mean, within-cohort SD units (Blom-scored), 95% CI",
    title = "The apoptotic configuration across the MYC x OXPHOS quadrants",
    subtitle = paste0(
      "Verdict on the rule fixed at 688539d before any readout: ", verdict,
      if (identical(verdict, V_OX)) " - confirms E11 in a group-level form" else "",
      "\nEXPLORATORY; not pre-registered. No outcome variable enters."),
    caption = paste0(
      "Standardised means from m1: readout ~ quadrant + PAM50 + purity + ",
      "leukocyte fraction. SCAN-B HAS NO PURITY OR LEUKOCYTE ESTIMATE, so ",
      "there m1 is quadrant + PAM50.\nQuadrant from the frozen STATE ",
      "constructor's own MYC and OXPHOS (GSVA) calls, each cohort on its own ",
      "medians. NO PRODUCT TERM is fitted, which is why no line joins the ",
      "pairs.\nComposite and arms are E33's signed mean-z on log2(linear + 1). ",
      "The two OXPHOS contrasts are near-guaranteed: the configuration's genes ",
      "and signs were chosen on OXPHOS in these cohorts.\nValues are ",
      "cohort-relative and are NEVER compared numerically between cohorts. ",
      "NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM: no cohort available ",
      "carries an OXPHOS- or BH3-directed intervention."),
    colour = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(size = 8),
    strip.placement = "outside",
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.caption = ggplot2::element_text(size = 5.6, hjust = 0, lineheight = 1.1),
    plot.caption.position = "plot")

.ensure_dir(DIR_DOC_FIGURES)
for (p in c(PATH_E34_FIG, PATH_E34_FIG_DOC)) {
  ggplot2::ggsave(p, fig, width = 8.4, height = 8.0, dpi = 300)
}
for (p in c(PATH_E34_FIG_PDF, PATH_E34_FIG_DOC_PDF)) {
  ggplot2::ggsave(p, fig, width = 8.4, height = 8.0, device = grDevices::cairo_pdf)
}
message("   ", PATH_E34_FIG, "\n   ", PATH_E34_FIG_DOC, "  (tracked)")
message("   ", PATH_E34_FIG_PDF, "\n   ", PATH_E34_FIG_DOC_PDF, "  (tracked)")

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. saving")

out <- list(
  verdict         = verdict,
  reading         = reading,
  reading_rule    = READING_RULE,
  verdict_cells   = VERDICT_CELLS,
  readout_summary = READOUT_SUMMARY,
  groups          = GROUPS,
  contrasts       = CONTRASTS,
  fits            = FITS,
  quadrant = list(cells = CELLS, asserts = ASSERTS,
                  agreement = dplyr::bind_rows(lapply(AGREE, `[[`, "summary")),
                  agreement_tables = lapply(AGREE, `[[`, "table"),
                  rho_frozen_myc_vs_M_a = rho_myc_frozen,
                  rho_frozen_ox_vs_repo = rho_ox_frozen),
  analysis_frames = lapply(DAT, function(d)
    d %>% dplyr::select(sample_id, quad_gsva, quad_mitopps, PAM50, purity,
                        leuko, PROLIF, dplyr::all_of(names(READOUTS)))),
  counts = list(
    tcga_n = length(ID_T), scanb_n = length(ID_S),
    tcga_no_quadrant = sum(is.na(DAT$TCGA$quad_gsva)),
    scanb_no_quadrant = sum(is.na(DAT$`SCAN-B`$quad_gsva)),
    tcga_pam50_na = sum(is.na(DAT$TCGA$PAM50)),
    scanb_pam50_na = sum(is.na(DAT$`SCAN-B`$PAM50)),
    tcga_purity_ok = sum(!is.na(DAT$TCGA$purity) & !is.na(DAT$TCGA$leuko))),
  spec = list(
    declaration = "docs/2026-09-30_E34_declaration.md, committed ALONE at 688539d",
    data_note   = "docs/2026-09-30_e34_data.md, 1345ec5",
    posture = paste("EXPLORATORY, not pre-registered. LARGELY CONFIRMATORY AND",
                    "PRESENTATIONAL: E11 in a group-level form. What is new is",
                    "the Q2 cell and the group-level form."),
    rule = READING_RULE,
    divergence = paste("The declaration names the primary as 'the configuration",
                       "composite, as E11 defines it'. E11 DEFINES NO COMPOSITE",
                       "- it works gene by gene over 44 genes. The primary here",
                       "is E33's signed composite, chosen by the author before",
                       "any data was read. Data note 1.1."),
    circularity = paste("The configuration's genes AND signs were read off",
                        "OXPHOS correlations in these same two cohorts, so the",
                        "two OXPHOS contrasts are close to guaranteed and",
                        "confirm almost nothing. The MYC contrasts are the",
                        "informative ones. Data note 1.2."),
    quadrant = paste("the frozen constructor called twice per cohort: once to",
                     "reproduce STATE, once with myc := oxphos to extract its",
                     "own OXPHOS call. Q3 == level 2 and Q4 == levels 3 + 4",
                     "asserted in both cohorts on both instruments."),
    scanb_quadrant = paste("no stored STATE exists for SCAN-B; the constructor",
                           "is called on SCAN-B's own scores with a CONSTANT",
                           "FALSE buffer, which the quadrant never uses, so its",
                           "cuts are SCAN-B's own medians and level 4 is empty."),
    tcga_55 = paste("55 TCGA patients have no quadrant because the frozen STATE",
                    "takes its complete cases from GISTIC BUFFER. Kept, not",
                    "repaired, so TCGA is identical to E33 and to the frozen",
                    "object."),
    scale = paste("Blom normal scores, qnorm((rank - 3/8)/(n + 1/4)), then unit",
                  "variance, pooled ONCE PER COHORT; fits subset and never",
                  "re-score. Estimates are in within-cohort SD units and are",
                  "NEVER compared numerically across cohorts."),
    models = paste("m1 quadrant + PAM50 + purity + leuko; m2 the same without",
                   "PAM50 within each subtype; m3 m1 + PROLIF_DISJOINT.",
                   "SCAN-B has no purity or leukocyte estimate (trap 2), so",
                   "there m1 is quadrant + PAM50 and m2 is quadrant alone."),
    no_product = "NO PRODUCT TERM anywhere. Four group means stay four group means.",
    m2_gate = paste("a subtype is READ in m2 only if all four quadrant cells",
                    "hold n >=", MIN_CELL_M2, "- the author's decision. Others",
                    "are reported and marked not read."),
    symbols = paste("SCAN-B resolves through scanb_pheno.rds$symbol_map, as",
                    "CLAUDE.md requires; the E02 map is asserted to agree."),
    no_outcome = "NO outcome or survival variable enters at any point",
    no_treatment = paste("NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM: no",
                         "cohort available carries an OXPHOS-directed or",
                         "BH3-directed intervention."),
    no_fdr = "no per-gene FDR, consistent with the rest of the project",
    n3 = "transcript scores; nothing is 'primed'"),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E34)
utils::write.csv(as.data.frame(CONTRASTS), PATH_E34_CSV, row.names = FALSE)
message("   ", PATH_E34)
message("   ", PATH_E34_CSV)
message("\nE34 done. VERDICT: ", verdict, "\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E34)

  cat(x$verdict, "\n\n"); cat(x$reading, "\n")
  utils::str(x$counts)

  # The quadrant: cell sizes, the declared assertions, instrument agreement.
  x$quadrant$cells %>%
    tidyr::pivot_wider(names_from = quadrant, values_from = n) %>%
    as.data.frame()
  x$quadrant$asserts %>% as.data.frame()
  x$quadrant$agreement %>% as.data.frame()
  x$quadrant$agreement_tables$TCGA

  # THE PRIMARY. All four contrasts, both cohorts, both instruments, m1 and m3.
  x$contrasts %>%
    dplyr::filter(readout_col == "config_comp", subtype == "all") %>%
    dplyr::transmute(cohort, instrument, model, contrast, kind,
                     est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3),
                     ci_excludes_0) %>%
    as.data.frame()

  # m2 - the result, not a check. Which subtypes were READ, and what they say.
  x$fits %>%
    dplyr::filter(readout_col == "config_comp", model == "m2") %>%
    dplyr::transmute(cohort, instrument, subtype, n, min_cell, readable, label) %>%
    as.data.frame()

  # Every reading cell the verdict is made of.
  x$verdict_cells %>%
    dplyr::transmute(cohort, instrument, model, subtype, n, readable, label) %>%
    as.data.frame()

  # The four group means the figure draws.
  x$groups %>%
    dplyr::filter(readout_col == "config_comp", model == "m1",
                  instrument == "gsva", subtype == "all") %>%
    dplyr::transmute(cohort, quadrant, n_cell, est = round(est, 3),
                     lo = round(lo, 3), hi = round(hi, 3)) %>%
    as.data.frame()

  # The arms, then the twelve. REPORTED, NEVER THE VERDICT.
  x$readout_summary %>% as.data.frame()
  x$contrasts %>%
    dplyr::filter(model == "m1", subtype == "all", instrument == "gsva",
                  contrast == "Q4 - Q2", grepl("^g_", readout_col)) %>%
    dplyr::transmute(cohort, readout, est = round(est, 3), lo = round(lo, 3),
                     hi = round(hi, 3), ci_excludes_0) %>%
    as.data.frame()
}
