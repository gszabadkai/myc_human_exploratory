# E41_metabric_ladder.R
# =============================================================================
# THE DECLARED LADDER IN METABRIC. IT MAKES NO DECISIONS.
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md,
# amended for METABRIC at ccadfa5. Read sections 4.5, 10.1, 10.1a, 10.2, 10.3,
# 10.3a, 10.4, 10.5 and 10.6 first. This script fits what section 5 specifies
# with the METABRIC spine of 10.1a, reports every rung whatever it shows, and
# interprets nothing.
#
# =============================================================================
# THE ESTIMAND IS NOT E40's, AND THIS IS THE FIRST THING TO READ
# =============================================================================
# E40 measured RESPONSE TO NEOADJUVANT CHEMOTHERAPY. This script measures
# PROGNOSIS UNDER MIXED, LARGELY HISTORICAL CARE. Of the 1,505 ER-positive
# patients, 1,111 (73.8%) had endocrine therapy and 147 (9.8%) chemotherapy, so
# the stratum METABRIC is here for is overwhelmingly a NON-CHEMOTHERAPY
# stratum. Per declaration 10.5:
#   NO SENTENCE MAY SAY THE pCR FINDING WAS REPLICATED, CONFIRMED OR EXTENDED
#   IN METABRIC, whatever the sign.
# The two arms sit beside each other as separate statements. They do not even
# share a spine (10.1a), so their m0 are not the same model.
#
# =============================================================================
# THREE TRAPS THAT FAIL SILENTLY (declaration 10.3a)
# =============================================================================
# 1. THE SURVIVAL OBJECTS ARE NAMED COUNTERINTUITIVELY.
#    Clinical_Overall_Survival_Data_from_METABRIC          = 888 events, ALL-CAUSE
#    Complete_METABRIC_Clinical_Survival_Data_from_METABRIC = 623 events, CAUSE-SPECIFIC
#    10.1 makes breast-cancer-specific survival PRIMARY, so the object with the
#    generic-sounding name is the primary one. Both counts are asserted below
#    and either mismatch is a hard stop.
# 2. SYMBOL HARMONISATION IS DIRECTIONAL, AND AN ALIAS IS NOT AUTOMATICALLY
#    SAFE. The map is applied CURRENT -> LEGACY and never reversed. An alias is
#    REJECTED where its target is the HGNC-approved symbol of a DIFFERENT gene:
#      COX7A2L -> SCAF1  REJECTED  SCAF1 is approved for SR-related CTD
#                                  associated factor 1 (58506); COX7A2L is 9167
#      POLR1B  -> RPA2   REJECTED  RPA2 is approved for replication protein A2
#                                  (6118); POLR1B is 84172
#      PRP4K   -> PRPF4B accepted  one gene, Entrez 8899
#    Taking SCAF1 would have put a splicing factor inside the OXPHOS exposure.
#    All three decisions are named gene by gene below and ASSERTED.
# 3. COX8C IS PRESENT IN METABRIC AND ABSENT FROM TCGA AND SCAN-B, so this gene
#    set is NOT a strict subset of the others. Actual counts only.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - NO model selection. Every rung is reported whatever it shows.
#   - NO ER x OX interaction, no three-way, no post-hoc term. Declared out.
#   - NO p value is described as a near miss, anywhere.
#   - NO per-gene FDR.
#   - It does NOT use MB3 in any form (declaration 10.2, 10.3).
#   - It does NOT substitute ATP5EP2 for ATP5F1E, or SCAF1 for COX7A2L, or
#     RPA2 for POLR1B.
#   - It does NOT reverse the symbol map.
#   - It does NOT reuse any score from another cohort, and compares no METABRIC
#     score numerically with a TCGA, SCAN-B or GSE25066 value.
#   - It does NOT fit the 6-gene signed configuration composite. BBC3 is absent
#     from METABRIC (4.5) and a five-gene stand-in is forbidden.
#   - It does NOT interpret. It prints numbers, never verdicts.
#   - It writes NOTHING to myc_human_validation, FROZEN at d3ac60e, or to
#     myc_mouse.
#
# SINGLE-INSTRUMENT THROUGHOUT. METABRIC is an array, so mitoPPS is unavailable
# (10.5). This has not cleared the arm's two-instrument bar and must never be
# written as though it had.
#
# COVERAGE. 69 of 89 OXPHOS subunits (77.5%). 10.4 measured the restriction
# rather than assuming it: Spearman 0.9960 (TCGA) and 0.9954 (SCAN-B) against
# the full set, largest movement in any manuscript reading -0.0174. It is
# treated as cosmetic. NO coefficient is compared across cohorts as a raw
# value (4.2).
#
# SCALE. METABRIC is ALREADY log2 and per-gene median-centred. GSVA wants
# log-scale input with kcdf = "Gaussian", so it is fed as-is; the per-gene
# centring is NOT undone because GSVA's per-gene ECDF is invariant to a
# per-gene location shift (4.5). mitoPPS, which would want linear counts, is
# not computed anywhere here.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

message("\nE41: the declared ladder in METABRIC\n", strrep("=", 78))
message("IT MAKES NO DECISIONS. Every rung is reported whatever it shows.")
message("THE ESTIMAND IS NOT E40's. See the header and declaration 10.5.")

PATH_MB_EXPR <- file.path(DIR_DATA, "menegollo_biclusters",
                          "METABRIC_DATA.RData")
PATH_MB_FORK <- file.path(DIR_DATA, "menegollo_biclusters",
                          "METABRIC_starting_data_final_corr_groups.Rdata")
PATH_SETDEFS <- file.path(DIR_RESULTS, "set_definitions.rds")
PATH_SCANB   <- file.path(DIR_DATA, "from_validation", "scanb_pheno.rds")
# READ-ONLY, from the FROZEN repo. Nothing is ever written there. This is read
# only to assert that M-a is the same 61 genes script 13 used.
PATH_G1      <- paste0("/Users/gs/code/myc_human_validation/results/",
                       "g1_overlap_audit.rds")
PATH_E41     <- file.path(DIR_RESULTS, "e41_metabric_ladder.rds")

# =============================================================================
# 0. CONSTANTS - every one from the declaration or from a recorded read, none
#    invented here. Each is a STOP, not a warning.
# =============================================================================
CI_LEVEL <- 0.95
CRIT     <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

COHORT <- "METABRIC"

# Declaration section 5. Each rung adds one term to the same spine.
RUNG_ADDS <- list(m0 = character(0),
                  m1 = "PROLIF",
                  m2 = c("PROLIF", "MYC"),
                  m3 = c("PROLIF", "MYC", "BUFFER_c"))
RUNGS <- names(RUNG_ADDS)

# Declaration 10.2. SECONDARY and DESCRIPTIVE, and NEITHER IS F3.
FORK_ADDS <- list(f_mb1 = "MB1_forkscale", f_mb2 = "MB2_forkscale")

SPINE_CANDIDATES <- c("subtype_term", "treatment")
VIF_CAVEAT       <- 5              # declaration 5.3

# script 13's own constants, reproduced so the GSVA call cannot drift (4.4).
GSVA_MIN_SET     <- 3L
NA_GENE_MAX_FRAC <- 0.05           # A5

# --- the file, as recorded in 4.5 and 10.3 ----------------------------------
N_ALL       <- 1981L
N_SYMBOLS   <- 23043L
MAX_DAYS    <- 9218L
DAYS_IN_YR  <- 365.25
E_BCSS      <- 623L                # cause-specific, the PRIMARY endpoint
E_OS        <- 888L                # all-cause, the declared sensitivity
OBJ_BCSS    <- "Complete_METABRIC_Clinical_Survival_Data_from_METABRIC"
OBJ_OS      <- "Clinical_Overall_Survival_Data_from_METABRIC"
OBJ_FEAT    <- "Complete_METABRIC_Clinical_Features_Data"
OBJ_EXPR    <- "METABRIC.data"
COL_SYMBOL  <- "external_gene_name"

# --- gene sets and their declared recovery (4.5) ----------------------------
SET_OX <- "OXPHOS subunits"
SET_PD <- "PROLIF_DISJOINT"
SET_MA <- "FELSHER__MITOSTRIP"     # "Felsher M-a", script 13's MYC estimator
N_SET     <- c(89L, 318L, 61L)
N_RECOVER <- c(69L, 293L, 59L)     # 69 NOT 70 - see ALIAS_RULE below
names(N_SET) <- names(N_RECOVER) <- c(SET_OX, SET_PD, SET_MA)
MA_ABSENT <- c("CTPS1", "POLR1B")
BUFFER_GENES <- c("MCL1", "BCL2L1")   # script 13's BUFFER_c, 2026-08-30 s6.1

# --- the alias rule, named GENE BY GENE and asserted (10.3a) ----------------
# An alias is REJECTED where its target is the HGNC-approved symbol of a
# DIFFERENT gene. These three decisions are hard-coded rather than looked up at
# runtime, so that the script carries its own audit trail and does not depend
# on whichever org.Hs.eg.db build happens to be installed.
ALIAS_ACCEPT <- c(PRP4K = "PRPF4B")
ALIAS_REJECT <- tibble::tribble(
  ~gene,     ~alias,   ~entrez_gene, ~entrez_alias, ~alias_is_approved_for,
  "COX7A2L", "SCAF1",  "9167",       "58506",   "SR-related CTD associated factor 1",
  "POLR1B",  "RPA2",   "84172",      "6118",    "replication protein A2")
# ATP5F1E: absent as ATP5F1E, as ATP5E and under every HGNC alias. The only
# ^ATP5E string in the matrix is the PSEUDOGENE ATP5EP2, which must NOT be
# substituted - array pseudogene probes carry cross-hybridisation noise.
EXCLUDED_GENES <- c("ATP5F1E", "COX7A2L", "POLR1B")
FORBIDDEN_SUBS <- c("ATP5EP2", "SCAF1", "RPA2")

# --- covariates (10.1a) -----------------------------------------------------
ER_PRIMARY   <- "ER_IHC_status"    # IHC. PRIMARY.
ER_SENSITIVE <- "ER.Expr"          # expression-derived. The sensitivity.
SUBTYPE_COL  <- "NOT_IN_OSLOVAL_Pam50Subtype"
TREAT_COL    <- "Treatment"
PAM50_REF    <- "LumA"             # reference level, fixed so output is stable
TREAT_REF    <- "NONE"
ER_IHC_N  <- c(pos = 1505L, neg = 435L)
ER_IHC_NA <- 41L
ER_EXP_N  <- c(`+` = 1512L, `-` = 469L)
ER_DISCORD        <- 127L
ER_DISCORD_NEGPOS <- 54L           # IHC negative, array positive
ER_DISCORD_POSNEG <- 73L           # IHC positive, array negative
PAM50_N <- c(LumA = 719L, LumB = 490L, Basal = 328L, Her2 = 238L,
             Normal = 200L, NC = 6L)
TREAT_N <- c(`HT/RT` = 609L, HT = 419L, NONE = 304L, RT = 231L,
             `CT/RT` = 170L, `CT/HT/RT` = 165L, CT = 52L, `CT/HT` = 31L)

# --- forkscale (10.2, 10.3) -------------------------------------------------
FORK_COLS <- c(MB1 = "MB1.forkscale", MB2 = "MB2.forkscale")
FORK_OBJ  <- "all.clinical.df"
# 10.3 trap 1: one Inf per .log column, at MB1 row 290, MB2 row 85, MB3 row
# 848. This script uses the UN-LOGGED forkscale, so no row is dropped. The
# rows are carried anyway and finiteness is asserted, because a silent Inf in
# the un-logged column would be a different and unrecorded defect.
FORK_INF_LOG_ROWS <- c(MB1 = 290L, MB2 = 85L, MB3 = 848L)
# F3-pre's figures, for a comparator only. F3 was specified in
# myc_human_validation, FROZEN at d3ac60e. THESE RUNGS ARE NOT F3.
F3PRE_RHO <- c(GSVA = 0.529, mitoPPS = 0.418)

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
.z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)
.pr <- function(d) print(as.data.frame(d), row.names = FALSE)
.sub <- function(s, w = 74) paste(strwrap(s, width = w), collapse = "\n   ")

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

for (p in c(PATH_MB_EXPR, PATH_MB_FORK, PATH_SETDEFS, PATH_SCANB, PATH_G1)) {
  .stop_if(file.exists(p), "absent: ", p)
}
sd_   <- readRDS(PATH_SETDEFS)
smap  <- readRDS(PATH_SCANB)$symbol_map
g1    <- readRDS(PATH_G1)

GS <- list(sd_$arm_sets[[SET_OX]], sd_$cov_sets[[SET_PD]], sd_$myc_sets[[SET_MA]])
names(GS) <- c(SET_OX, SET_PD, SET_MA)
for (nm in names(GS)) {
  .stop_if(length(GS[[nm]]) == N_SET[[nm]],
           nm, " is ", length(GS[[nm]]), " genes, not the declared ",
           N_SET[[nm]], ". The on-disk set has moved; stop.")
}
# M-a must be the SAME 61 genes script 13 scored, read from the frozen repo.
.stop_if(identical(sort(GS[[SET_MA]]), sort(g1$estimators_stripped$FELSHER)),
         SET_MA, " does not reproduce the frozen g1 estimators_stripped$FELSHER. ",
         "The MYC estimator would differ from the one E40 used. STOP.")
# D7: PROLIF_DISJOINT is defined against M-a and must share no gene with it.
.stop_if(length(intersect(GS[[SET_MA]], GS[[SET_PD]])) == 0L,
         "M-a and ", SET_PD, " intersect; D7 says they are disjoint")
message("   gene sets: ", paste(sprintf("%s %d", names(GS),
        vapply(GS, length, integer(1))), collapse = " | "))
message("   M-a verified identical to the frozen g1 set (61 genes), and ",
        "disjoint from ", SET_PD)
message("   symbol_map: ", length(smap), " entries, applied CURRENT -> LEGACY")

# =============================================================================
# 2. A. PREPARE - assertions throughout. Any failure stops.
# =============================================================================
message("\n2. A. prepare")

# --- 2.1 the expression matrix ----------------------------------------------
mb <- new.env(parent = emptyenv())
load(PATH_MB_EXPR, envir = mb)
for (o in c(OBJ_EXPR, OBJ_FEAT, OBJ_BCSS, OBJ_OS)) {
  .stop_if(exists(o, envir = mb, inherits = FALSE),
           PATH_MB_EXPR, " does not carry ", o)
}
raw  <- get(OBJ_EXPR, envir = mb)
feat <- get(OBJ_FEAT, envir = mb)
bcss <- get(OBJ_BCSS, envir = mb)
osv  <- get(OBJ_OS,   envir = mb)

sym <- as.character(raw[[COL_SYMBOL]])
.stop_if(length(sym) == N_SYMBOLS, "symbol column is ", length(sym),
         " long, not the declared ", N_SYMBOLS)
.stop_if(anyDuplicated(sym) == 0L,
         "duplicated gene symbols in METABRIC; 4.5 records 0 and says NO ",
         "probe collapse is needed. One is needed after all. STOP.")
.stop_if(sum(is.na(sym) | sym == "") == 0L, "empty or NA gene symbols")

# Sample columns: MB.0000 -> MB-0000. make.names mangling, undone.
samp_cols <- setdiff(colnames(raw), COL_SYMBOL)
.stop_if(length(samp_cols) == N_ALL, "expression carries ", length(samp_cols),
         " sample columns, not ", N_ALL)
ids <- gsub("\\.", "-", samp_cols)
.stop_if(anyDuplicated(ids) == 0L, "sample ids are not unique after un-mangling")
.stop_if(sum(ids %in% feat$sample) == N_ALL,
         "only ", sum(ids %in% feat$sample), " of ", N_ALL,
         " expression columns join to ", OBJ_FEAT, "$sample. 4.5 declares ",
         N_ALL, " of ", N_ALL, ".")

M <- as.matrix(raw[, samp_cols, drop = FALSE])
rownames(M) <- sym
colnames(M) <- ids
storage.mode(M) <- "double"
rm(raw); invisible(gc())
message("   expression: ", nrow(M), " genes x ", ncol(M), " samples; ",
        N_ALL, " of ", N_ALL, " ids join to ", OBJ_FEAT)
message("   ALREADY log2 and per-gene median-centred; NOT transformed, ",
        "NOT un-centred (4.5)")

# --- 2.2 symbol harmonisation, CURRENT -> LEGACY, never reversed ------------
# Returns the MATRIX ROWNAMES a set resolves to, with every decision reported.
UNIVERSE <- rownames(M)

.resolve <- function(genes, label) {
  direct <- intersect(genes, UNIVERSE)
  left   <- setdiff(genes, direct)
  # the map is applied one way only and the RESULT is checked against the matrix
  cand   <- ifelse(left %in% names(smap), unname(smap[left]), NA_character_)
  ok_map <- !is.na(cand) & cand %in% UNIVERSE
  left2  <- left[!ok_map]
  # the alias step: ONLY the accepted mapping, and only where it lands
  ok_al  <- left2 %in% names(ALIAS_ACCEPT) &
            unname(ALIAS_ACCEPT[left2]) %in% UNIVERSE
  absent <- sort(left2[!ok_al])

  # Per-row PROVENANCE, so the guards below can distinguish "this symbol is a
  # legitimate member of the set" from "this symbol was substituted in". A
  # global ban on a symbol would be wrong: see the RPA2 case in 10.3a.
  prov <- dplyr::bind_rows(
    tibble::tibble(set_gene = direct,        matrix_row = direct,
                   route = "direct"),
    tibble::tibble(set_gene = left[ok_map],  matrix_row = cand[ok_map],
                   route = "map"),
    tibble::tibble(set_gene = left2[ok_al],
                   matrix_row = unname(ALIAS_ACCEPT[left2[ok_al]]),
                   route = "alias"))
  rows <- prov$matrix_row
  .stop_if(anyDuplicated(rows) == 0L,
           label, ": a gene resolves to a matrix row already taken by another ",
           "(", paste(rows[duplicated(rows)], collapse = ", "),
           "). That would double-count it. STOP.")

  # A REJECTED gene must not be resolved by any route.
  for (i in seq_len(nrow(ALIAS_REJECT))) {
    g <- ALIAS_REJECT$gene[i]; a <- ALIAS_REJECT$alias[i]
    .stop_if(!(g %in% prov$set_gene[prov$route != "direct"]),
             label, ": ", g, " was resolved via ",
             prov$matrix_row[prov$set_gene == g][1], ". Its only METABRIC ",
             "route is ", a, ", the approved symbol of ",
             ALIAS_REJECT$alias_is_approved_for[i], " (Entrez ",
             ALIAS_REJECT$entrez_alias[i], "), NOT of ", g, " (",
             ALIAS_REJECT$entrez_gene[i], "). 10.3a rejects it. STOP.")
  }
  # A FORBIDDEN SUBSTITUTE must not appear as a substitute. It MAY appear by
  # the direct route, because then it is the gene the set actually names.
  for (f in FORBIDDEN_SUBS) {
    sub_in <- prov$route != "direct" & prov$matrix_row == f
    .stop_if(!any(sub_in),
             label, ": ", f, " entered as a SUBSTITUTE for ",
             paste(prov$set_gene[sub_in], collapse = ", "),
             " (10.3a forbids it). STOP.")
  }
  list(rows = rows, prov = prov, n_direct = length(direct),
       n_map = sum(ok_map), n_alias = sum(ok_al), absent = absent,
       alias_used = if (any(ok_al))
         paste(sprintf("%s -> %s", left2[ok_al],
                       unname(ALIAS_ACCEPT[left2[ok_al]])),
               collapse = ", ") else "",
       direct_but_also_forbidden =
         intersect(direct, c(FORBIDDEN_SUBS, ALIAS_REJECT$alias)))
}

RES <- lapply(names(GS), function(nm) .resolve(GS[[nm]], nm))
names(RES) <- names(GS)

RECOVERY <- dplyr::bind_rows(lapply(names(RES), function(nm) {
  r <- RES[[nm]]
  tibble::tibble(set = nm, n_set = N_SET[[nm]], direct = r$n_direct,
                 via_map = r$n_map, via_alias = r$n_alias,
                 total = length(r$rows),
                 pct = round(100 * length(r$rows) / N_SET[[nm]], 1),
                 absent = length(r$absent), alias_used = r$alias_used)
}))
.pr(RECOVERY)
for (nm in names(RES)) {
  .stop_if(length(RES[[nm]]$rows) == N_RECOVER[[nm]],
           nm, " recovers ", length(RES[[nm]]$rows), " of ", N_SET[[nm]],
           ", not the declared ", N_RECOVER[[nm]],
           ". 4.5's figure and this run disagree; find the cause before ",
           "anything is scored.")
  .stop_if(length(RES[[nm]]$rows) >= GSVA_MIN_SET, nm, ": below the GSVA floor")
}
.stop_if(identical(sort(RES[[SET_MA]]$absent), sort(MA_ABSENT)),
         "M-a's absent genes are ", paste(RES[[SET_MA]]$absent, collapse = ", "),
         ", not the declared ", paste(MA_ABSENT, collapse = ", "))
message("   ASSERTED: 69 of 89, 293 of 318, 59 of 61 (4.5)")
EXCLUDED_TAB <- tibble::tibble(
  gene = EXCLUDED_GENES,
  reason = c(paste0("absent as ATP5F1E, as ATP5E and under every alias; the ",
                    "only ^ATP5E string is the PSEUDOGENE ATP5EP2, NOT substituted"),
             paste0("its only METABRIC route is SCAF1, the approved symbol of ",
                    "SR-related CTD associated factor 1 (58506) - REJECTED"),
             paste0("its only METABRIC route is RPA2, the approved symbol of ",
                    "replication protein A2 (6118) - REJECTED")))
message("   the three EXCLUDED genes, named gene by gene (10.3a):")
for (i in seq_len(nrow(EXCLUDED_TAB)))
  message("     ", EXCLUDED_TAB$gene[i], ": ",
          .sub(EXCLUDED_TAB$reason[i], 66))
message("   ASSERTED: no rejected gene was resolved by any route, and none of ",
        paste(FORBIDDEN_SUBS, collapse = ", "),
        " entered as a SUBSTITUTE for anything")
# The guard is on the ROUTE, not on the symbol. A symbol that is a forbidden
# substitute for one gene can still be a legitimate direct member of a set in
# its own right, and banning it outright would silently alter that set.
DIRECT_COLLISIONS <- dplyr::bind_rows(lapply(names(RES), function(nm)
  if (length(RES[[nm]]$direct_but_also_forbidden))
    tibble::tibble(set = nm, symbol = RES[[nm]]$direct_but_also_forbidden,
                   route = "direct") else NULL))
if (nrow(DIRECT_COLLISIONS)) {
  message("   ", .sub(paste0(
    "AND a symbol that is forbidden AS A SUBSTITUTE is still allowed as a ",
    "DIRECT member, because then it is the gene the set actually names. ",
    "Here that applies to: ",
    paste(sprintf("%s in %s", DIRECT_COLLISIONS$symbol, DIRECT_COLLISIONS$set),
          collapse = "; "),
    ". RPA2 is replication protein A2, a proliferation gene, so its presence ",
    "in PROLIF_DISJOINT is correct and is NOT the POLR1B substitution. A ",
    "global ban on the symbol would have removed a legitimate member and ",
    "changed the proliferation score (10.3a).")))
} else {
  DIRECT_COLLISIONS <- tibble::tibble(set = character(0), symbol = character(0),
                                      route = character(0))
}
message("   the 20 OXPHOS subunits absent: ",
        paste(RES[[SET_OX]]$absent, collapse = ", "))

# BUFFER_c's two genes, read from the matrix directly as script 13 does
.stop_if(all(BUFFER_GENES %in% UNIVERSE),
         "BUFFER_c needs ", paste(BUFFER_GENES, collapse = " and "),
         "; missing: ", paste(setdiff(BUFFER_GENES, UNIVERSE), collapse = ", "),
         ". Fit m0 to m2 only and say so, rather than substituting anything.")
message("   BUFFER_c genes present: ", paste(BUFFER_GENES, collapse = ", "),
        " - so every rung of the ladder is buildable")
# 4.5: BBC3 is absent, and a five-gene stand-in is forbidden. Asserted so that
# nobody later adds the configuration to this script without seeing it.
BBC3_PRESENT <- "BBC3" %in% UNIVERSE
message("   BBC3 present: ", BBC3_PRESENT, " - the 6-gene signed ",
        "configuration is NOT fitted here and a 5-gene stand-in is ",
        "FORBIDDEN (4.5)")

# --- 2.3 survival, asserted BEFORE use --------------------------------------
# The naming is counterintuitive. 10.3a. Either mismatch is a hard stop.
.stop_if(sum(osv$status == 1L) == E_OS,
         OBJ_OS, " has ", sum(osv$status == 1L), " events, not the declared ",
         E_OS, ". The two survival objects may have been swapped. STOP.")
.stop_if(sum(bcss$status == 1L) == E_BCSS,
         OBJ_BCSS, " has ", sum(bcss$status == 1L),
         " events, not the declared ", E_BCSS, ". STOP.")
.stop_if(identical(osv$time, bcss$time),
         "the two survival objects' time vectors differ; 10.3a records them ",
         "as identical and differing only in status")
.stop_if(max(bcss$time) == MAX_DAYS, "max time is ", max(bcss$time),
         " not ", MAX_DAYS, "; the unit may not be days")
.stop_if(all(bcss$time > 0), "non-positive survival times")
.stop_if(setequal(bcss$sample, feat$sample),
         "the survival and clinical sample sets differ")
message("   survival asserted BEFORE use (10.3a):")
message("     PRIMARY   ", OBJ_BCSS, " = ", E_BCSS,
        " events, CAUSE-SPECIFIC")
message("     sensitivity ", OBJ_OS, " = ", E_OS, " events, ALL-CAUSE")
message("     time in DAYS, max ", MAX_DAYS, " (", round(MAX_DAYS/DAYS_IN_YR, 1),
        " y); converted to YEARS")

# --- 2.4 covariates and the two ER calls (10.1a) ----------------------------
CL <- feat %>%
  dplyr::transmute(
    sample,
    er_primary   = dplyr::case_when(.data[[ER_PRIMARY]] == "pos" ~ "P",
                                    .data[[ER_PRIMARY]] == "neg" ~ "N",
                                    TRUE ~ NA_character_),
    er_sensitive = dplyr::case_when(.data[[ER_SENSITIVE]] == "+" ~ "P",
                                    .data[[ER_SENSITIVE]] == "-" ~ "N",
                                    TRUE ~ NA_character_),
    subtype_pam50 = .data[[SUBTYPE_COL]],
    her2_expr     = Her2.Expr,
    treatment_raw = .data[[TREAT_COL]])

.stop_if(sum(CL$er_primary == "P", na.rm = TRUE) == ER_IHC_N[["pos"]] &&
         sum(CL$er_primary == "N", na.rm = TRUE) == ER_IHC_N[["neg"]] &&
         sum(is.na(CL$er_primary)) == ER_IHC_NA,
         ER_PRIMARY, " gives ", sum(CL$er_primary == "P", na.rm = TRUE), "/",
         sum(CL$er_primary == "N", na.rm = TRUE), " with ",
         sum(is.na(CL$er_primary)), " missing, not the declared ",
         ER_IHC_N[["pos"]], "/", ER_IHC_N[["neg"]], " with ", ER_IHC_NA)
.stop_if(sum(CL$er_sensitive == "P") == ER_EXP_N[["+"]] &&
         sum(CL$er_sensitive == "N") == ER_EXP_N[["-"]] &&
         sum(is.na(CL$er_sensitive)) == 0L,
         ER_SENSITIVE, " does not give the declared ", ER_EXP_N[["+"]], "/",
         ER_EXP_N[["-"]], " with none missing")
# table() carries a class and dimnames, so identical() against a named
# integer vector is ALWAYS FALSE and the assertion would never fire. Compared
# name by name instead.
.counts_match <- function(v, declared) {
  tb  <- table(v, useNA = "no")
  got <- stats::setNames(as.integer(tb), names(tb))
  identical(sort(names(got)), sort(names(declared))) &&
    all(got[names(declared)] == declared)
}
.stop_if(.counts_match(CL$subtype_pam50, PAM50_N),
         "the PAM50 composition differs from 10.1a's: got ",
         paste(sprintf("%s %d", names(table(CL$subtype_pam50)),
                       as.integer(table(CL$subtype_pam50))), collapse = " "))
.stop_if(.counts_match(CL$treatment_raw, TREAT_N),
         "the Treatment composition differs from 10.1a's: got ",
         paste(sprintf("%s %d", names(table(CL$treatment_raw)),
                       as.integer(table(CL$treatment_raw))), collapse = " "))

# The discordance, recorded BEFORE the sensitivity is read (10.1a).
DISCORD <- CL %>% dplyr::filter(!is.na(er_primary), er_primary != er_sensitive)
D_NEGPOS <- sum(DISCORD$er_primary == "N" & DISCORD$er_sensitive == "P")
D_POSNEG <- sum(DISCORD$er_primary == "P" & DISCORD$er_sensitive == "N")
.stop_if(nrow(DISCORD) == ER_DISCORD && D_NEGPOS == ER_DISCORD_NEGPOS &&
         D_POSNEG == ER_DISCORD_POSNEG,
         "the ER discordance is ", nrow(DISCORD), " (", D_NEGPOS, " / ",
         D_POSNEG, "), not the declared ", ER_DISCORD, " (",
         ER_DISCORD_NEGPOS, " / ", ER_DISCORD_POSNEG, ")")
message("   the two ER calls (10.1a):")
.pr(tibble::tibble(
  call = c(paste0(ER_PRIMARY, "  PRIMARY"), paste0(ER_SENSITIVE, "  sensitivity")),
  positive = c(ER_IHC_N[["pos"]], ER_EXP_N[["+"]]),
  negative = c(ER_IHC_N[["neg"]], ER_EXP_N[["-"]]),
  missing  = c(ER_IHC_NA, 0L)))
message("   they disagree on ", nrow(DISCORD), " of ",
        sum(!is.na(CL$er_primary)), " callable patients: ", D_NEGPOS,
        " IHC-negative/array-positive, ", D_POSNEG, " the other way")
message("   ", .sub(paste0(
  "The ", ER_IHC_NA, " IHC-missing patients fall OUT of the primary ER ",
  "stratification entirely, so its two strata sum to ",
  sum(!is.na(CL$er_primary)), " and not ", N_ALL, ". They do NOT fall out of ",
  "the pooled ladder, which does not condition on ER. Under ", ER_SENSITIVE,
  " nothing falls out, so the sensitivity's strata sum to ", N_ALL,
  ": a ", ER_IHC_NA, "-patient difference in the analysis set travels with ",
  "the change of definition and must not be read as an effect of it.")))

# --- 2.5 the spine's subtype term: PAM50, NC kept as its own level ----------
# 10.1a. No patient is dropped for subtype. The two-level analogue is
# CROSS-TABBED and NOT FITTED.
CL <- CL %>% dplyr::mutate(
  subtype_term = stats::relevel(factor(subtype_pam50), ref = PAM50_REF),
  treatment    = stats::relevel(factor(treatment_raw), ref = TREAT_REF),
  # the analogue, built and reported only
  subtype2_analogue = dplyr::case_when(
    is.na(er_primary) | is.na(her2_expr) ~ NA_character_,
    er_primary == "P" & her2_expr == "-" ~ "HRpos_HER2neg",
    er_primary == "N" & her2_expr == "-" ~ "TNBC",
    TRUE ~ "HER2pos"))
.stop_if(PAM50_REF %in% levels(CL$subtype_term) &&
         levels(CL$subtype_term)[1] == PAM50_REF,
         "the PAM50 reference level is not ", PAM50_REF)
.stop_if(sum(is.na(CL$subtype_term)) == 0L,
         "PAM50 has missing values; 10.1a records it complete at ", N_ALL)
message("   spine subtype term = PAM50, ", nlevels(CL$subtype_term),
        " levels, reference ", PAM50_REF, ", NC KEPT as its own level")
.pr(CL %>% dplyr::count(subtype_term, name = "n"))
message("   the two-level analogue, CROSS-TABBED and NOT FITTED (10.1a):")
print(table(analogue = CL$subtype2_analogue, PAM50 = CL$subtype_pam50,
            useNA = "ifany"))
message("   ", .sub(paste0(
  "It is not fitted because it would be near-constant inside each ER stratum, ",
  "the collapse 7.1 documents for GSE25066, and because HER2_IHC_status ",
  "carries only 821 non-missing of ", N_ALL, " so it would have to come from ",
  "Her2.Expr. THE METABRIC SPINE THEREFORE DIFFERS FROM GSE25066's AND THE ",
  "TWO COHORTS' m0 ARE NOT THE SAME MODEL (10.1a).")))

# --- 2.6 the analysis frame -------------------------------------------------
DM <- CL %>%
  dplyr::left_join(bcss %>% dplyr::transmute(sample, bcss_time = time,
                                             bcss_event = status),
                   by = "sample") %>%
  dplyr::left_join(osv %>% dplyr::transmute(sample, os_time = time,
                                            os_event = status),
                   by = "sample") %>%
  dplyr::mutate(bcss_time_years = bcss_time / DAYS_IN_YR,
                os_time_years   = os_time   / DAYS_IN_YR)
.stop_if(nrow(DM) == N_ALL, "the frame is ", nrow(DM), " rows, not ", N_ALL)
.stop_if(sum(is.na(DM$bcss_event)) == 0L && sum(is.na(DM$os_event)) == 0L,
         "survival did not join to every row")
.stop_if(sum(DM$bcss_event) == E_BCSS && sum(DM$os_event) == E_OS,
         "the joined event counts are not ", E_BCSS, " and ", E_OS)

rm(list = ls(envir = mb), envir = mb); rm(mb); invisible(gc())

# =============================================================================
# 3. B. SCORE - ONE GSVA call, script 13's parameters exactly (4.4)
# =============================================================================
# Declaration 4.4 and 4.5. log2 input, kcdf Gaussian, universe pinned with
# .PIN_A/.PIN_B, minSize 3L, maxSize Inf, and a silently dropped set is a hard
# stop. METABRIC's matrix is already log2 and per-gene median-centred, which is
# what GSVA wants and is NOT undone. mitoPPS is not computed anywhere: it would
# want linear counts, and METABRIC is an array (10.5).
message("\n3. B. scoring, one GSVA call")

.impute <- function(Mx, label) {
  na_g <- rowSums(is.na(Mx))
  drop <- na_g > NA_GENE_MAX_FRAC * ncol(Mx)
  if (any(drop)) {
    message("      ", label, ": dropped ", sum(drop), " gene(s) above ",
            100 * NA_GENE_MAX_FRAC, "% NA")
    Mx <- Mx[!drop, , drop = FALSE]
  }
  n_imp <- sum(is.na(Mx))
  if (n_imp) {
    med <- apply(Mx, 1L, stats::median, na.rm = TRUE)
    idx <- which(is.na(Mx), arr.ind = TRUE)
    Mx[idx] <- med[idx[, "row"]]
    message("      ", label, ": imputed ", n_imp, " value(s) at the gene median")
  }
  .stop_if(!anyNA(Mx), label, ": NA survived .impute()")
  Mx
}
MI <- .impute(M, COHORT)
# .impute() may drop rows, which would invalidate a resolved set. Re-checked
# rather than assumed.
for (nm in names(RES)) {
  .stop_if(all(RES[[nm]]$rows %in% rownames(MI)),
           nm, ": .impute() removed ", sum(!RES[[nm]]$rows %in% rownames(MI)),
           " of its genes. Recovery is no longer ", N_RECOVER[[nm]], ". STOP.")
}
.stop_if(all(BUFFER_GENES %in% rownames(MI)),
         ".impute() removed a BUFFER_c gene")
# M is ~365 MB and is dead once MI exists; UNIVERSE was taken from it before
# .impute() and the resolved sets have just been re-checked against MI.
rm(M); invisible(gc())

sets <- list(RES[[SET_OX]]$rows, RES[[SET_PD]]$rows, RES[[SET_MA]]$rows)
names(sets) <- c("OX", "PROLIF", "MYC")
u <- rownames(MI)
sets[[".PIN_A"]] <- u[c(TRUE, FALSE)]
sets[[".PIN_B"]] <- u[c(FALSE, TRUE)]

par <- GSVA::gsvaParam(exprData = MI, geneSets = sets, kcdf = "Gaussian",
                       minSize = GSVA_MIN_SET, maxSize = Inf)
S <- GSVA::gsva(par, verbose = FALSE)
for (nm in c("OX", "PROLIF", "MYC")) {
  .stop_if(nm %in% rownames(S), "GSVA silently dropped ", nm)
}
.stop_if(identical(colnames(S), colnames(MI)), "GSVA reordered the samples")
message("   GSVA: ", nrow(S), " rows (3 sets + 2 pins) x ", ncol(S),
        " samples, kcdf Gaussian, minSize ", GSVA_MIN_SET)

SCORES <- tibble::tibble(
  sample   = colnames(S),
  OX_raw     = as.numeric(S["OX", ]),
  PROLIF_raw = as.numeric(S["PROLIF", ]),
  MYC_raw    = as.numeric(S["MYC", ]),
  z_MCL1     = .z(as.numeric(MI["MCL1", ])),
  z_BCL2L1   = .z(as.numeric(MI["BCL2L1", ]))) %>%
  # Declaration 4.2: every exposure standardised to WITHIN-COHORT SD, so
  # coefficients read per within-cohort SD and never as raw values across
  # cohorts. Unlike E39/E40, OX is scored here rather than inherited from
  # script 13's .z(), so all three are standardised here.
  dplyr::mutate(OX = .z(OX_raw), PROLIF = .z(PROLIF_raw), MYC = .z(MYC_raw),
                BUFFER_c = (z_MCL1 + z_BCL2L1) / 2)
.stop_if(all(abs(vapply(SCORES[c("OX","PROLIF","MYC")], stats::sd, numeric(1))
                 - 1) < 1e-9), "an exposure is not at sd = 1")

DM <- DM %>% dplyr::left_join(SCORES, by = "sample")
.stop_if(sum(is.na(DM[c("OX", "PROLIF", "MYC", "BUFFER_c")])) == 0L,
         "a score did not join to every patient")
message("   standardised to within-cohort SD (4.2): ",
        paste(sprintf("%s sd %.6f", c("OX","PROLIF","MYC"),
              vapply(DM[c("OX","PROLIF","MYC")], stats::sd, numeric(1))),
              collapse = " | "))
message("   BUFFER_c = mean(z(MCL1), z(BCL2L1)), script 13's definition. ",
        "NO result from m3 may be described as about 'buffering capacity' (5.1)")
CORR_EXPOSURES <- tibble::tibble(
  pair = c("OX vs PROLIF", "OX vs MYC", "MYC vs PROLIF"),
  spearman = c(stats::cor(DM$OX, DM$PROLIF, method = "spearman"),
               stats::cor(DM$OX, DM$MYC, method = "spearman"),
               stats::cor(DM$MYC, DM$PROLIF, method = "spearman")))
message("   the exposures against each other (5.3's collinearity question):")
.pr(CORR_EXPOSURES %>% dplyr::mutate(spearman = round(spearman, 3)))

# =============================================================================
# 4. Stratum descriptions - n, events, median follow-up, BEFORE any fit
# =============================================================================
message("\n4. n, events and follow-up for every stratum, BEFORE any fit")

.describe <- function(d, label) {
  tibble::tibble(stratum = label, n = nrow(d),
    bcss_events = sum(d$bcss_event), os_events = sum(d$os_event),
    median_fu_y = stats::median(d$bcss_time_years),
    median_fu_censored_y =
      if (any(d$bcss_event == 0L))
        stats::median(d$bcss_time_years[d$bcss_event == 0L]) else NA_real_,
    max_fu_y = max(d$bcss_time_years))
}
STRATA <- dplyr::bind_rows(
  .describe(DM, "pooled, all"),
  .describe(DM[!is.na(DM$er_primary) & DM$er_primary == "P", ],
            "ER-positive (ER_IHC_status) PRIMARY"),
  .describe(DM[!is.na(DM$er_primary) & DM$er_primary == "N", ],
            "ER-negative (ER_IHC_status)"),
  .describe(DM[DM$er_sensitive == "P", ], "ER-positive (ER.Expr) sensitivity"),
  .describe(DM[DM$er_sensitive == "N", ], "ER-negative (ER.Expr) sensitivity"))
.pr(STRATA %>% dplyr::mutate(dplyr::across(dplyr::ends_with("_y"),
                                           ~ round(.x, 2))))
ER_POS_ROW <- which(STRATA$stratum == "ER-positive (ER_IHC_status) PRIMARY")
.stop_if(length(ER_POS_ROW) == 1L, "the ER-positive primary stratum row is ",
         "not identifiable in STRATA by name")
message("   ", .sub(paste0(
  "THE ER-POSITIVE STRATUM IS WHAT METABRIC IS HERE FOR (10.6): ",
  STRATA$n[ER_POS_ROW], " patients carrying ", STRATA$bcss_events[ER_POS_ROW],
  " cause-specific events, against GSE25066's 40. ",
  "Median follow-up ", round(STRATA$median_fu_y[ER_POS_ROW], 2),
  " y, maximum ", round(STRATA$max_fu_y[ER_POS_ROW], 2), " y.")))
message("   median_fu_y is the median of ALL times; median_fu_censored_y is ",
        "the median among the censored. Neither is a reverse-KM estimate.")

# =============================================================================
# 5. Spine audit - before any fit. Constant terms are omitted EXPLICITLY.
# =============================================================================
# Declaration 5.2: never left for R to drop silently.
message("\n5. spine audit (5.2)")

.n_lev <- function(v) dplyr::n_distinct(v[!is.na(v)])
.audit <- function(d, label) {
  tibble::tibble(set = label,
                 term = c("OX", SPINE_CANDIDATES),
                 n_levels = c(.n_lev(d$OX), .n_lev(d$subtype_term),
                              .n_lev(d$treatment))) %>%
    dplyr::mutate(constant = n_levels <= 1L)
}
AUDIT <- dplyr::bind_rows(
  .audit(DM, "pooled"),
  .audit(DM[!is.na(DM$er_primary) & DM$er_primary == "P", ], "ER-pos (IHC)"),
  .audit(DM[!is.na(DM$er_primary) & DM$er_primary == "N", ], "ER-neg (IHC)"),
  .audit(DM[DM$er_sensitive == "P", ], "ER-pos (ER.Expr)"),
  .audit(DM[DM$er_sensitive == "N", ], "ER-neg (ER.Expr)"))
.pr(AUDIT)
CONSTANTS <- AUDIT %>% dplyr::filter(constant)
if (nrow(CONSTANTS)) {
  message("   CONSTANT and therefore OMITTED EXPLICITLY:")
  .pr(CONSTANTS)
} else {
  message("   nothing is constant anywhere: the full spine is fitted in every ",
          "set, including inside each ER stratum. This is what PAM50 buys ",
          "over the two-level analogue (10.1a).")
}
.spine_for <- function(label) {
  a <- AUDIT[AUDIT$set == label & AUDIT$term %in% SPINE_CANDIDATES, ]
  a$term[!a$constant]
}

# =============================================================================
# 6. The estimator helpers
# =============================================================================
# VIF. OX is always a single continuous term, so its Df is 1 and its GVIF IS
# the ordinary VIF - the two coincide at 1 df. car::vif is called ONCE per fit
# and both the value and the measure name come from that one call. (E40 called
# it twice per fit; that defect is not reproduced here.)
#
# TWO THINGS THE DRY RUN EXPOSED, both handled here rather than left to
# accumulate:
#
# 1. car::vif WARNS "No intercept: vifs may not be sensible" on every Cox fit,
#    because a Cox model has no intercept BY DESIGN. The warning is spurious
#    here: the VIF is computed from the correlation matrix of the coefficient
#    estimates, which needs no intercept. IT IS COUNTED AND REPORTED ONCE,
#    deliberately, rather than suppressed silently or allowed to pile up.
#    THIS IS THE SOURCE OF E40's 46 UNIDENTIFIED WARNINGS: 40 fits there, one
#    warning per car::vif call, and E40 called it twice per fit.
# 2. car::vif ERRORS with "there are aliased coefficients in the model" when a
#    factor carries a level with zero observations in the subset being fitted -
#    which happens in every ER-negative stratum, where PAM50 has 5 or 4 of its
#    6 levels. The fix is droplevels() applied per fitting set and REPORTED
#    (see .fit_ladder); an empty level contributes nothing and would otherwise
#    give an NA coefficient silently.
.VIF_STATE <- new.env(parent = emptyenv())
.VIF_STATE$no_intercept_warnings <- 0L

.vif_ox <- function(fit) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2L) return(list(vif = NA_real_, df = 1, err = "",
                           measure = "undefined (single predictor)"))
  err <- ""
  v <- withCallingHandlers(
    tryCatch(car::vif(fit),
             error = function(e) { err <<- conditionMessage(e); NULL }),
    warning = function(w) {
      if (grepl("No intercept", conditionMessage(w), fixed = TRUE)) {
        .VIF_STATE$no_intercept_warnings <-
          .VIF_STATE$no_intercept_warnings + 1L
        invokeRestart("muffleWarning")
      }
    })
  if (is.null(v)) return(list(vif = NA_real_, df = NA_real_, err = err,
                              measure = "car::vif failed"))
  if (is.matrix(v))
    list(vif = unname(v["OX", 1]), df = unname(v["OX", 2]), err = "",
         measure = "GVIF (generalised; OX has Df 1 so GVIF == VIF)")
  else list(vif = unname(v[["OX"]]), df = 1, err = "", measure = "VIF")
}

.tidy_ox <- function(fit, rung, endpoint, set, omitted = "") {
  co <- summary(fit)$coefficients
  .stop_if("OX" %in% rownames(co), set, " ", rung, ": no OX coefficient")
  est <- co["OX", 1]
  se  <- co["OX", "se(coef)"]
  p   <- co["OX", ncol(co)]
  vf  <- .vif_ox(fit)
  tibble::tibble(
    endpoint = endpoint, set = set, rung = rung,
    n = fit$n, events = fit$nevent,
    estimate = est, se = se, HR = exp(est),
    ci_lo = est - CRIT * se, ci_hi = est + CRIT * se,
    hr_lo = exp(est - CRIT * se), hr_hi = exp(est + CRIT * se), p = p,
    vif_ox = vf$vif, vif_df = vf$df, vif_measure = vf$measure,
    vif_error = vf$err,
    vif_above_caveat = isTRUE(vf$vif > VIF_CAVEAT),
    terms = paste(attr(stats::terms(fit), "term.labels"), collapse = " + "),
    omitted = omitted)
}

# The fits are kept, not just their coefficients: cox.zph in section F needs
# the model objects and refitting them there would be a second, separate fit.
FITS <- list()

# EMPTY FACTOR LEVELS ARE DROPPED EXPLICITLY AND REPORTED, never left for R
# to alias away silently. A level with zero observations in this fitting set
# contributes no information, gives an NA coefficient, and breaks car::vif with
# "there are aliased coefficients in the model" - which is why every
# ER-negative VIF came back NA before this was added. Dropping an EMPTY level
# changes no estimate; it is bookkeeping, not a modelling decision.
DROPPED_LEVELS <- list()

.drop_empty <- function(d, set, spine) {
  for (v in spine) {
    if (!is.factor(d[[v]])) next
    gone <- setdiff(levels(d[[v]]), unique(as.character(d[[v]])))
    if (length(gone)) {
      DROPPED_LEVELS[[length(DROPPED_LEVELS) + 1L]] <<- tibble::tibble(
        set = set, term = v, dropped = paste(gone, collapse = ", "),
        n_dropped = length(gone), n_kept = nlevels(droplevels(d[[v]])))
      message("   ", set, ": ", v, " has ", length(gone),
              " EMPTY level(s) here - ", paste(gone, collapse = ", "),
              " - dropped explicitly (0 observations each, so no estimate ",
              "changes)")
      d[[v]] <- droplevels(d[[v]])
    }
  }
  d
}

.fit_ladder <- function(d, set, endpoint, spine, adds = RUNG_ADDS,
                        omitted = "") {
  d  <- .drop_empty(d, set, spine)
  tv <- if (endpoint == "BCSS") "bcss_time_years" else "os_time_years"
  ev <- if (endpoint == "BCSS") "bcss_event" else "os_event"
  key <- paste(endpoint, set, sep = "|")
  FITS[[key]] <<- stats::setNames(vector("list", length(adds)), names(adds))
  rows <- lapply(names(adds), function(rg) {
    rhs <- c("OX", spine, adds[[rg]])
    f <- stats::reformulate(rhs,
           response = sprintf("survival::Surv(%s, %s)", tv, ev))
    fit <- survival::coxph(f, data = d)
    FITS[[key]][[rg]] <<- fit
    .tidy_ox(fit, rg, endpoint, set, omitted)
  })
  dplyr::bind_rows(rows)
}

# se is in the saved object and the csv; it is left out of the printed view so
# the table fits a terminal without wrapping onto a second block.
.show <- function(tb) .pr(tb %>% dplyr::transmute(set, rung, n, events,
    logHR = round(estimate, 4), HR = round(HR, 4),
    ci_HR = paste0("[", round(hr_lo, 3), ", ", round(hr_hi, 3), "]"),
    p = signif(p, 3), vif = round(vif_ox, 3)))

# =============================================================================
# 7. C. THE LADDER - Cox, breast-cancer-specific survival PRIMARY
# =============================================================================
message("\n7. C. the ladder, Cox, ", OBJ_BCSS, " (", E_BCSS,
        " events) PRIMARY")

SP_POOLED <- .spine_for("pooled")
LADDER <- .fit_ladder(DM, "pooled", "BCSS", SP_POOLED,
                      omitted = paste(setdiff(SPINE_CANDIDATES, SP_POOLED),
                                      collapse = ", "))
.show(LADDER)
message("   spine = OX + ", paste(SP_POOLED, collapse = " + "),
        "; CI shown on the HR scale, estimate and se on the log-hazard scale")
message("   VIF measure: ", LADDER$vif_measure[1])
message("   ", .sub(paste0(
  "What the steps answer (5, 5.3). m0 is the unconditional respiratory ",
  "coefficient. m0 -> m1 is how much of it is proliferation. m1 -> m2 is ",
  "DECLARED NOT a decomposition by variable: at this collinearity ",
  "'whether what remains after proliferation is MYC' is DECLARED OUT and ",
  "must not be written. What is declared IN is whether the OX coefficient is ",
  "robust to adding a covariate strongly correlated with one already in the ",
  "model. The same applies to m2 -> m3.")))

# --- the declared sensitivity: overall survival, same ladder ----------------
message("\n   C2. declared sensitivity: ", OBJ_OS, " (", E_OS,
        " events, ALL-CAUSE)")
message("   ", .sub(paste0(
  "10.1 makes cause-specific survival primary because twenty years of ",
  "follow-up in an older cohort loads OS with non-cancer death. OS is ",
  "reported BESIDE it, never instead of it.")))
LADDER_OS <- .fit_ladder(DM, "pooled", "OS", SP_POOLED,
                         omitted = paste(setdiff(SPINE_CANDIDATES, SP_POOLED),
                                         collapse = ", "))
.show(LADDER_OS)

# =============================================================================
# 8. D. FORKSCALE RUNGS - MB1 SECONDARY, MB2 DESCRIPTIVE, NEITHER IS F3
# =============================================================================
message("\n8. D. forkscale rungs (10.2). NEITHER IS F3.")

fk <- new.env(parent = emptyenv())
load(PATH_MB_FORK, envir = fk)
.stop_if(exists(FORK_OBJ, envir = fk, inherits = FALSE),
         PATH_MB_FORK, " does not carry ", FORK_OBJ)
acd <- get(FORK_OBJ, envir = fk)
.stop_if(all(c("sample", FORK_COLS) %in% names(acd)),
         FORK_OBJ, " is missing: ",
         paste(setdiff(c("sample", FORK_COLS), names(acd)), collapse = ", "))
.stop_if(nrow(acd) == N_ALL, FORK_OBJ, " is ", nrow(acd), " rows, not ", N_ALL)
FORK <- tibble::tibble(sample = acd$sample,
                       MB1_forkscale = acd[[FORK_COLS[["MB1"]]]],
                       MB2_forkscale = acd[[FORK_COLS[["MB2"]]]])
rm(list = ls(envir = fk), envir = fk); rm(fk, acd); invisible(gc())

# 10.3 trap 1. The UN-LOGGED forkscale is used, so no row is dropped. The .log
# Inf rows are reported anyway, and finiteness here is ASSERTED rather than
# assumed: a silent Inf in the un-logged column would be a different defect.
message("   using the UN-LOGGED ", paste(FORK_COLS, collapse = " and "),
        " as stored, so NO row is dropped")
message("   10.3's one Inf per .log column is at MB1 row ",
        FORK_INF_LOG_ROWS[["MB1"]], ", MB2 row ", FORK_INF_LOG_ROWS[["MB2"]],
        ", MB3 row ", FORK_INF_LOG_ROWS[["MB3"]],
        " - no .log variant is used here")
for (v in c("MB1_forkscale", "MB2_forkscale")) {
  .stop_if(all(is.finite(FORK[[v]])),
           v, " carries ", sum(!is.finite(FORK[[v]])),
           " non-finite value(s) in the UN-LOGGED column. 10.3 records the Inf ",
           "as confined to the .log variants. STOP and find it.")
}
.stop_if(!any(grepl("MB3", names(FORK))), "an MB3 column entered FORK")
message("   ASSERTED: finite throughout, and NO MB3 column is carried (10.2)")

DM <- DM %>% dplyr::left_join(FORK, by = "sample")
.stop_if(sum(is.na(DM$MB1_forkscale)) == 0L &&
         sum(is.na(DM$MB2_forkscale)) == 0L,
         "forkscale did not join to every patient")

# 10.3 trap 3: forkscale is severely skewed by construction, so SPEARMAN and
# never Pearson for any correlation involving it.
FORK_RHO <- tibble::tibble(
  pair = c("OX vs MB1_forkscale", "OX vs MB2_forkscale",
           "MB1_forkscale vs MB2_forkscale", "PROLIF vs MB1_forkscale",
           "MYC vs MB1_forkscale"),
  spearman = c(stats::cor(DM$OX, DM$MB1_forkscale, method = "spearman"),
               stats::cor(DM$OX, DM$MB2_forkscale, method = "spearman"),
               stats::cor(DM$MB1_forkscale, DM$MB2_forkscale, method = "spearman"),
               stats::cor(DM$PROLIF, DM$MB1_forkscale, method = "spearman"),
               stats::cor(DM$MYC, DM$MB1_forkscale, method = "spearman")))
.pr(FORK_RHO %>% dplyr::mutate(spearman = round(spearman, 3)))
message("   SPEARMAN only; forkscale is severely skewed (10.3 trap 3)")
message("   ", .sub(paste0(
  "F3-pre's comparator figures are Spearman ", F3PRE_RHO[["GSVA"]],
  " (GSVA) and ", F3PRE_RHO[["mitoPPS"]], " (mitoPPS), which left the ",
  "question INTERMEDIATE on both instruments. The METABRIC figure sits beside ",
  "them as a third reading in the cohort the forks were derived from. ",
  "F3 WAS SPECIFIED IN myc_human_validation, FROZEN AT d3ac60e. THESE RUNGS ",
  "ARE NOT F3 and must be named as separate exploratory analyses in every ",
  "sentence that mentions them.")))

# The rungs themselves, added to the full ladder's top rung.
FORK_RUNGS <- list(
  f_mb1 = c(RUNG_ADDS$m3, FORK_ADDS$f_mb1),
  f_mb2 = c(RUNG_ADDS$m3, FORK_ADDS$f_mb2))
LADDER_FORK <- .fit_ladder(DM, "fork", "BCSS", SP_POOLED, adds = FORK_RUNGS)
.show(LADDER_FORK)
message("   f_mb1 is m3 + MB1_forkscale, DECLARED SECONDARY")
message("   f_mb2 is m3 + MB2_forkscale, DECLARED DESCRIPTIVE - a FIRST ",
        "REPORT carried as a ride-along, not a result the arm rests on (10.2)")
message("   ", .sub(paste0(
  "Two limits on the MB2 rung, recorded before it ran (10.2). E32 placed MYC ",
  "activity on MB2 at +0.688 and E33 placed the configuration on MB1, so an ",
  "MB2 result contributes to the companion paper's framework and not to this ",
  "one's argument. And the MB2 FORK CONTRAST is confounded with ER by ",
  "construction, so ANY MB2 fork contrast is reported WITHIN ER STRATA ONLY, ",
  "never pooled - see section 9.")))

# =============================================================================
# 9. E. ER STRATIFICATION - the point of this cohort
# =============================================================================
# Declaration 7 and 10.1a. Descriptive, not a test. NO ER x OX interaction:
# declared out, and not fitted here either.
message("\n9. E. ER stratification - what METABRIC is here for")

# fit_label is separate from label because label is the key into AUDIT, while
# fit_label is the key into FITS. Reusing one label for both would make the
# fork fits overwrite this stratum's m0-m3 fits and silently drop them from the
# proportional-hazards table.
.fit_stratum <- function(er_var, level, label, adds = RUNG_ADDS,
                         fit_label = label) {
  d <- DM[!is.na(DM[[er_var]]) & DM[[er_var]] == level, ]
  sp <- .spine_for(label)
  om <- setdiff(SPINE_CANDIDATES, sp)
  if (length(om))
    message("   ", label, ": OMITTING ", paste(om, collapse = ", "),
            " - constant in this stratum (5.2)")
  .stop_if(!(paste("BCSS", fit_label, sep = "|") %in% names(FITS)),
           "fit label '", fit_label, "' is already taken; it would overwrite ",
           "an earlier set of fits")
  .fit_ladder(d, fit_label, "BCSS", sp,
              adds = adds, omitted = paste(om, collapse = ", "))
}

ER_STRAT <- dplyr::bind_rows(
  .fit_stratum("er_primary",   "P", "ER-pos (IHC)"),
  .fit_stratum("er_primary",   "N", "ER-neg (IHC)"),
  .fit_stratum("er_sensitive", "P", "ER-pos (ER.Expr)"),
  .fit_stratum("er_sensitive", "N", "ER-neg (ER.Expr)"))
.show(ER_STRAT)
message("   ", .sub(paste0(
  "The ER-positive IHC stratum carries ", STRATA$bcss_events[ER_POS_ROW],
  " cause-specific events in ", STRATA$n[ER_POS_ROW], " patients. The ",
  ER_IHC_NA, " IHC-missing patients are in NEITHER primary stratum, so the ",
  "two sum to ", sum(!is.na(DM$er_primary)), " and not ", N_ALL,
  "; under ER.Expr nothing falls out and the two sum to ", N_ALL,
  ". The two ER calls disagree on ", ER_DISCORD, " of ",
  sum(!is.na(DM$er_primary)), " callable patients (10.1a), so the ",
  "sensitivity reclassifies ", ER_DISCORD, " patients and is not a ",
  "re-reading of the same strata.")))
message("   NO ER x OX interaction is fitted. Declared out (section 7).")

# --- the MB2 fork rung, WITHIN ER STRATA ONLY (10.2) -----------------------
message("\n   E2. the MB2 forkscale rung, WITHIN ER STRATA ONLY (10.2)")
message("   ", .sub(paste0(
  "Menegollo's Figure 3C labels MB2_UF the Myc/miRNA ER-negative state and ",
  "the lower fork ER-positive, so upper-versus-lower is largely ",
  "ER-negative-versus-positive. A POOLED MB2 CONTRAST IS THEREFORE FORBIDDEN. ",
  "That reading of Figure 3C is AUTHOR-SOURCED and unverified in this repo ",
  "(10.2): the Menegollo PDF is in no repository, so no script can check it.")))
ER_FORK <- dplyr::bind_rows(
  .fit_stratum("er_primary", "P", "ER-pos (IHC)", adds = FORK_RUNGS,
               fit_label = "ER-pos (IHC) fork"),
  .fit_stratum("er_primary", "N", "ER-neg (IHC)", adds = FORK_RUNGS,
               fit_label = "ER-neg (IHC) fork"))
.show(ER_FORK)
ER_FORK_RHO <- dplyr::bind_rows(lapply(c("P", "N"), function(lv) {
  d <- DM[!is.na(DM$er_primary) & DM$er_primary == lv, ]
  tibble::tibble(stratum = paste0("ER-", ifelse(lv == "P", "pos", "neg"),
                                  " (IHC)"), n = nrow(d),
    rho_OX_MB1 = stats::cor(d$OX, d$MB1_forkscale, method = "spearman"),
    rho_OX_MB2 = stats::cor(d$OX, d$MB2_forkscale, method = "spearman"))
}))
.pr(ER_FORK_RHO %>% dplyr::mutate(dplyr::across(dplyr::starts_with("rho"),
                                                ~ round(.x, 3))))

# =============================================================================
# 10. F. PROPORTIONAL HAZARDS - reported whether or not it holds
# =============================================================================
# Declaration 6.3. cox.zph on EVERY model, global test and the OX term.
# Global failures are EXPECTED over 25 years of follow-up and are NOT a reason
# to change the model.
message("\n10. F. proportional hazards (cox.zph), every model")

PH <- dplyr::bind_rows(lapply(names(FITS), function(key) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  dplyr::bind_rows(lapply(names(FITS[[key]]), function(rg) {
    z  <- survival::cox.zph(FITS[[key]][[rg]])
    tb <- z$table
    tibble::tibble(endpoint = parts[1], set = parts[2], rung = rg,
                   global_chisq = unname(tb["GLOBAL", "chisq"]),
                   global_df    = unname(tb["GLOBAL", "df"]),
                   global_p     = unname(tb["GLOBAL", "p"]),
                   ox_chisq     = unname(tb["OX", "chisq"]),
                   ox_df        = unname(tb["OX", "df"]),
                   ox_p         = unname(tb["OX", "p"])) %>%
      dplyr::mutate(ox_ph_holds = ox_p >= 0.05,
                    global_ph_holds = global_p >= 0.05)
  }))
}))
.pr(PH %>% dplyr::transmute(endpoint, set, rung,
    global_p = signif(global_p, 3), global_ph_holds,
    ox_p = signif(ox_p, 3), ox_ph_holds))
message("   ", .sub(paste0(
  "Declaration 6.3 requires PH be reported whether or not it holds. GLOBAL ",
  "failures are EXPECTED over up to ", round(MAX_DAYS/DAYS_IN_YR, 1),
  " years of follow-up and are NOT a reason to change the model. ",
  sum(!PH$global_ph_holds), " of ", nrow(PH), " models fail globally and ",
  sum(!PH$ox_ph_holds), " of ", nrow(PH), " fail FOR OX.")))

# --- the RMST branch, mirroring E40 exactly --------------------------------
# Declaration 6.3 asks for RMST at 5 and 10 y ALONGSIDE the hazard ratio where
# PH fails FOR OX. This script does NOT invent a way to do it. The two blockers
# E40 recorded apply here unchanged:
#   - survRM2 is NOT installed (checked, not installed by this script), and
#   - OX is a standardised CONTINUOUS exposure. An RMST contrast needs two
#     arms, so reporting one would require dichotomising OX, and NO cut is
#     declared anywhere in the declaration.
# Making that cut here would be a modelling decision, and this script makes
# none. The branch REPORTS what is required and what is missing.
PH_OX_FAILS <- PH %>% dplyr::filter(!ox_ph_holds)
HAVE_SURVRM2 <- requireNamespace("survRM2", quietly = TRUE)
RMST_NOTE <- if (nrow(PH_OX_FAILS) == 0L) {
  paste0("PH holds for OX at every rung in every model; declaration 6.3's ",
         "RMST branch is not triggered.")
} else {
  paste0(
    "PH FAILS FOR OX in ", nrow(PH_OX_FAILS), " of ", nrow(PH),
    " models (", paste(sprintf("%s/%s/%s", PH_OX_FAILS$endpoint,
                               PH_OX_FAILS$set, PH_OX_FAILS$rung),
                       collapse = ", "),
    "). Declaration 6.3 requires RMST at 5 and 10 y ALONGSIDE the hazard ",
    "ratio there. THIS SCRIPT DID NOT COMPUTE IT, for the two reasons E40 ",
    "recorded, neither of which it may resolve on its own: survRM2 is ",
    if (HAVE_SURVRM2) "installed, BUT " else "not installed in this ",
    if (HAVE_SURVRM2) "" else "environment, and ",
    "OX is a standardised CONTINUOUS exposure, so an RMST contrast would ",
    "require dichotomising it at a cut that the declaration does not ",
    "specify. Choosing that cut is a modelling decision and this script ",
    "makes none. The hazard ratio is reported with the PH result beside it, ",
    "and the limitation is recorded rather than engineered away.")
}
message("\n   ", .sub(RMST_NOTE))

# =============================================================================
# 10b. The two car::vif artefacts, reported once rather than left to pile up
# =============================================================================
ALL_ROWS <- dplyr::bind_rows(LADDER, LADDER_OS, LADDER_FORK, ER_STRAT, ER_FORK)
VIF_NOTE <- paste0(
  "car::vif raised the spurious warning 'No intercept: vifs may not be ",
  "sensible' ", .VIF_STATE$no_intercept_warnings, " time(s), once per call. ",
  "A Cox model has NO INTERCEPT BY DESIGN and the VIF is computed from the ",
  "correlation matrix of the coefficient estimates, which needs none, so the ",
  "warning carries no information here. It was counted and muffled ",
  "deliberately. THIS IS THE SOURCE OF E40's 46 UNIDENTIFIED WARNINGS: 40 ",
  "fits, one warning per car::vif call, and E40 called it twice per fit.")
message("\n10b. ", .sub(VIF_NOTE))
VIF_FAILED <- ALL_ROWS %>% dplyr::filter(is.na(vif_ox))
if (nrow(VIF_FAILED)) {
  message("   car::vif still returned NO value for ", nrow(VIF_FAILED),
          " of ", nrow(ALL_ROWS), " fits:")
  .pr(VIF_FAILED %>% dplyr::count(set, rung, vif_measure, vif_error,
                                  name = "n"))
} else {
  message("   car::vif returned a value for all ", nrow(ALL_ROWS), " fits.")
}
if (length(DROPPED_LEVELS)) {
  message("   empty factor levels dropped explicitly, by set:")
  .pr(dplyr::bind_rows(DROPPED_LEVELS))
} else {
  message("   no spine factor had an empty level in any fitting set.")
}

# =============================================================================
# 11. G. Save
# =============================================================================
message("\n11. G. saving")

readr::write_csv(RECOVERY,    file.path(DIR_TABLES, "E41_gene_recovery.csv"))
readr::write_csv(STRATA,      file.path(DIR_TABLES, "E41_strata.csv"))
readr::write_csv(AUDIT,       file.path(DIR_TABLES, "E41_spine_audit.csv"))
readr::write_csv(LADDER,      file.path(DIR_TABLES, "E41_ladder_bcss.csv"))
readr::write_csv(LADDER_OS,   file.path(DIR_TABLES, "E41_ladder_os.csv"))
readr::write_csv(LADDER_FORK, file.path(DIR_TABLES, "E41_fork_rungs.csv"))
readr::write_csv(ER_STRAT,    file.path(DIR_TABLES, "E41_er_stratified.csv"))
readr::write_csv(ER_FORK,     file.path(DIR_TABLES, "E41_er_fork_rungs.csv"))
readr::write_csv(PH,          file.path(DIR_TABLES, "E41_proportional_hazards.csv"))
readr::write_csv(FORK_RHO,    file.path(DIR_TABLES, "E41_forkscale_rho.csv"))

saveRDS(list(
  recovery      = RECOVERY,
  excluded      = EXCLUDED_TAB,
  dropped_levels = if (length(DROPPED_LEVELS))
                     dplyr::bind_rows(DROPPED_LEVELS) else
                     tibble::tibble(set = character(0), term = character(0),
                                    dropped = character(0)),
  alias_reject  = ALIAS_REJECT,
  alias_accept  = ALIAS_ACCEPT,
  direct_collisions = DIRECT_COLLISIONS,
  er_calls      = tibble::tibble(
                    call = c(ER_PRIMARY, ER_SENSITIVE),
                    role = c("PRIMARY", "sensitivity"),
                    positive = c(ER_IHC_N[["pos"]], ER_EXP_N[["+"]]),
                    negative = c(ER_IHC_N[["neg"]], ER_EXP_N[["-"]]),
                    missing = c(ER_IHC_NA, 0L)),
  er_discordant = DISCORD,
  strata        = STRATA,
  audit         = AUDIT,
  exposure_corr = CORR_EXPOSURES,
  ladder        = LADDER,
  ladder_os     = LADDER_OS,
  ladder_fork   = LADDER_FORK,
  fork_rho      = FORK_RHO,
  er_stratified = ER_STRAT,
  er_fork       = ER_FORK,
  er_fork_rho   = ER_FORK_RHO,
  ph            = PH,
  ph_ox_fails   = PH_OX_FAILS,
  rmst_note     = RMST_NOTE,
  vif_note      = VIF_NOTE,
  vif_failed    = VIF_FAILED,
  rungs         = RUNG_ADDS,
  fork_rungs    = FORK_RUNGS,
  spec = list(
    declaration = paste0("docs/2026-10-02_E39_respiratory_axis_decomposition",
                         "_declaration.md, amended for METABRIC at ccadfa5"),
    posture = paste0("EXPLORATORY and POST-HOC. A descriptive decomposition, ",
                     "not a hypothesis test. Nothing here carries a pass, a ",
                     "fail or a falsification criterion. No rung was ",
                     "selected; every rung is reported whatever it shows."),
    provenance = paste0(
      "Expression and clinical from data/menegollo_biclusters/",
      "METABRIC_DATA.RData, untracked. Source recorded on the AUTHOR'S ",
      "AUTHORITY as Synapse syn1757063, collection syn1688369, per the ",
      "Methods of Menegollo et al. 2024. NOTHING INSIDE THE FILE NAMES ",
      "SYNAPSE: no README object, no date stamp, no identifier, no checksum. ",
      "The three clinical tables carry a readr col_spec, which is CONSISTENT ",
      "WITH a delimited-text distribution and IS NOT EVIDENCE of one. The ",
      "md5 in data/menegollo_biclusters/README.md remains the only check. ",
      "Forkscale from METABRIC_starting_data_final_corr_groups.Rdata, also ",
      "untracked. The cBioPortal datahub fetch that 10.1 once specified was ",
      "NEVER PERFORMED (4.5)."),
    estimand = paste0(
      "NOT E40's. E40 measured RESPONSE TO NEOADJUVANT CHEMOTHERAPY; this ",
      "measures PROGNOSIS UNDER MIXED, LARGELY HISTORICAL CARE. Of 1,505 ",
      "ER-positive patients, 1,111 (73.8%) had endocrine therapy and 147 ",
      "(9.8%) chemotherapy; of 435 ER-negative, 270 (62.1%) had ",
      "chemotherapy. NO SENTENCE MAY SAY THE pCR FINDING WAS REPLICATED, ",
      "CONFIRMED OR EXTENDED IN METABRIC, whatever the sign (10.5). The two ",
      "arms do not share a spine, so their m0 are not the same model."),
    traps = paste0(
      "10.3a, all three. (1) The survival objects are named ",
      "counterintuitively: ", OBJ_OS, " is the 888-event ALL-CAUSE endpoint ",
      "and ", OBJ_BCSS, " the 623-event CAUSE-SPECIFIC one, which is the ",
      "PRIMARY. Both counts asserted before any fit. (2) Symbol ",
      "harmonisation is CURRENT -> LEGACY and never reversed, and an alias ",
      "is REJECTED where its target is the approved symbol of a different ",
      "gene: COX7A2L -> SCAF1 (58506, SR-related CTD associated factor 1) ",
      "and POLR1B -> RPA2 (6118, replication protein A2) are both rejected; ",
      "PRP4K -> PRPF4B (8899, one gene) is accepted. THE GUARD IS ON THE ",
      "ROUTE, NOT ON THE SYMBOL: RPA2 is a legitimate DIRECT member of ",
      "PROLIF_DISJOINT (replication protein A2 is a proliferation gene), so ",
      "banning the symbol outright would have removed a real member and ",
      "changed the proliferation score. ATP5F1E is excluded and ",
      "ATP5EP2, a pseudogene, is NOT substituted. (3) COX8C is present here ",
      "and absent from TCGA and SCAN-B, so this gene set is not a strict ",
      "subset of the others."),
    instrument_check = paste0(
      "10.4. The 69-gene restriction was measured, not assumed: against the ",
      "full 89-gene set it correlates at Spearman 0.9960 in TCGA (n=1,095) ",
      "and 0.9954 in SCAN-B (n=3,207), with the largest movement in any ",
      "manuscript reading -0.0174. 8 of the 10 movements are toward zero; ",
      "the two exceptions are the MitoCarta half of the 44-gene comparator, ",
      "by +0.0002 and +0.0021. The restriction is treated as COSMETIC, ",
      "subject to 4.2's standing rule that no coefficient is compared across ",
      "cohorts as a raw value. An earlier 70-gene version of this check ",
      "accepted COX7A2L -> SCAF1 and is superseded."),
    instrument = paste0(
      "SINGLE-INSTRUMENT THROUGHOUT. METABRIC is an array, so mitoPPS is ",
      "unavailable (10.5). This has NOT cleared the arm's two-instrument bar ",
      "and must never be written as though it had."),
    spine = paste0(
      "10.1a. The subtype term is PAM50 (", SUBTYPE_COL, "), complete at ",
      N_ALL, ", reference ", PAM50_REF, ", with the 6 NC patients KEPT as ",
      "their own level so no patient is dropped for subtype. The arm's ",
      "two-level HRpos_HER2neg/TNBC analogue is CROSS-TABBED AND NOT FITTED: ",
      "it would be near-constant inside each ER stratum, and HER2_IHC_status ",
      "carries only 821 non-missing. THE METABRIC SPINE THEREFORE DIFFERS ",
      "FROM GSE25066's AND THE TWO COHORTS' m0 ARE NOT THE SAME MODEL."),
    er = paste0(
      "10.1a. ", ER_PRIMARY, " is PRIMARY (", ER_IHC_N[["pos"]], "/",
      ER_IHC_N[["neg"]], ", ", ER_IHC_NA, " missing); ", ER_SENSITIVE,
      " is the declared sensitivity (", ER_EXP_N[["+"]], "/", ER_EXP_N[["-"]],
      ", none missing). They disagree on ", ER_DISCORD, " of 1,940 callable ",
      "patients, ", ER_DISCORD_NEGPOS, " IHC-negative/array-positive and ",
      ER_DISCORD_POSNEG, " the other way. The ", ER_IHC_NA,
      " IHC-missing fall OUT of the primary stratification entirely, so the ",
      "pooled set and the sum of the primary strata differ by ", ER_IHC_NA,
      " patients; under the sensitivity nothing falls out. That difference ",
      "travels with the change of definition and is not an effect of it. NO ",
      "ER x OX interaction was fitted: declared out (section 7)."),
    collinearity = paste0(
      "5.3. The m1 -> m2 step does NOT decompose by variable: 'whether what ",
      "remains after proliferation is MYC' is DECLARED OUT and must not be ",
      "written. What is declared IN is whether the OX coefficient is robust ",
      "to adding a covariate strongly correlated with one already in the ",
      "model. The same applies to m2 -> m3. VIF for OX is reported at every ",
      "rung; above ", VIF_CAVEAT, " it is a caveat on that rung, not a ",
      "reason to drop a term. Neither MYC nor PROLIF may be dropped."),
    buffer = paste0(
      "5.1. BUFFER_c is the mean of z(MCL1) and z(BCL2L1). NO result from m3 ",
      "may be described as being about 'buffering capacity'; it is about ",
      "adjustment for two anti-apoptotic transcripts."),
    forkscale = paste0(
      "10.2. MB1_forkscale is DECLARED SECONDARY and MB2_forkscale ",
      "DESCRIPTIVE, a first report carried as a ride-along. NEITHER IS F3 - ",
      "F3 was specified in myc_human_validation, FROZEN at d3ac60e, and ",
      "these must be named separate exploratory analyses in every sentence ",
      "that mentions them. F3-pre's comparators are Spearman 0.529 (GSVA) ",
      "and 0.418 (mitoPPS). MB3 IS NOT USED IN ANY FORM, so 10.3's ",
      "unverified MB3.forkscale == MB3.pc1.rev / MB3.index stays open. The ",
      "UN-LOGGED forkscale is used, so 10.3's one Inf per .log column (MB1 ",
      "row 290, MB2 row 85, MB3 row 848) drops no row; finiteness of the ",
      "un-logged columns is asserted anyway. Spearman only: forkscale is ",
      "severely skewed by construction. ANY MB2 FORK CONTRAST IS REPORTED ",
      "WITHIN ER STRATA ONLY, never pooled, because Figure 3C makes the ",
      "contrast largely ER-negative-versus-positive - a reading that is ",
      "AUTHOR-SOURCED and unverified in this repo."),
    configuration = paste0(
      "4.5. BBC3 is ABSENT from METABRIC. The 6-gene signed configuration ",
      "composite is NOT fitted anywhere here, and a five-gene stand-in is ",
      "FORBIDDEN. Any future METABRIC analysis touching the configuration is ",
      "blocked on BBC3."),
    rmst = RMST_NOTE,
    vif = VIF_NOTE,
    scale = paste0(
      "METABRIC is ALREADY log2 and per-gene median-centred. GSVA is fed it ",
      "as-is with kcdf = 'Gaussian'; the per-gene centring is NOT undone ",
      "because GSVA's per-gene ECDF is invariant to a per-gene location ",
      "shift (4.5). mitoPPS, which would want linear counts, is not computed ",
      "anywhere here. All exposures standardised to WITHIN-COHORT SD (4.2)."),
    frozen = paste0("Nothing written to myc_human_validation (d3ac60e) or ",
                    "myc_mouse. g1_overlap_audit.rds was read ONLY to assert ",
                    "that M-a is the same 61 genes script 13 scored.")),
  built = Sys.time()), PATH_E41)

message("   ", PATH_E41)
message("   10 tables in ", DIR_TABLES)
message("\nE41 done. Numbers only; nothing here is interpreted.")
message("THE ESTIMAND IS NOT E40's. Declaration 10.5.\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E41)

  # A. what was recovered, and the three genes that were refused.
  x$recovery %>% as.data.frame()
  x$alias_reject %>% as.data.frame()
  x$excluded %>% as.data.frame()

  # The two ER calls and the 127 they disagree on, before reading any stratum.
  x$er_calls %>% as.data.frame()
  x$er_discordant %>% dplyr::count(er_primary, er_sensitive, subtype_pam50,
                                   name = "n") %>% as.data.frame()

  # n, events and follow-up per stratum. Row 2 is what METABRIC is here for.
  x$strata %>% as.data.frame()
  x$audit %>% as.data.frame()

  # C. the ladder on cause-specific survival. Watch the m0 -> m1 step, and the
  # VIF for OX at m2 and m3.
  x$ladder %>% dplyr::transmute(rung, n, events, HR = round(HR, 4),
      lo = round(hr_lo, 3), hi = round(hr_hi, 3), p = signif(p, 3),
      vif_ox = round(vif_ox, 3), vif_above_caveat) %>% as.data.frame()

  # C2. the same ladder on all-cause survival, the declared sensitivity.
  x$ladder_os %>% dplyr::transmute(rung, n, events, HR = round(HR, 4),
      lo = round(hr_lo, 3), hi = round(hr_hi, 3), p = signif(p, 3)) %>%
    as.data.frame()

  # D. the forkscale rungs and the correlations beside them. NEITHER IS F3.
  x$fork_rho %>% as.data.frame()
  x$ladder_fork %>% dplyr::transmute(rung, n, events, HR = round(HR, 4),
      lo = round(hr_lo, 3), hi = round(hr_hi, 3), p = signif(p, 3),
      vif_ox = round(vif_ox, 3)) %>% as.data.frame()

  # E. the ER strata. The ER-positive IHC rows are the reason for the cohort.
  x$er_stratified %>% dplyr::transmute(set, rung, n, events,
      HR = round(HR, 4), lo = round(hr_lo, 3), hi = round(hr_hi, 3),
      p = signif(p, 3), omitted) %>% as.data.frame()

  # E2. the MB2 rung WITHIN strata only, never pooled.
  x$er_fork %>% dplyr::transmute(set, rung, n, events, HR = round(HR, 4),
      lo = round(hr_lo, 3), hi = round(hr_hi, 3), p = signif(p, 3)) %>%
    as.data.frame()
  x$er_fork_rho %>% as.data.frame()

  # F. proportional hazards, and the RMST branch that was NOT computed.
  x$ph %>% as.data.frame()
  x$ph_ox_fails %>% as.data.frame()
  cat(x$rmst_note, "\n")

  # The exposures against each other, 5.3's question.
  x$exposure_corr %>% as.data.frame()

  utils::str(x$spec)
}
