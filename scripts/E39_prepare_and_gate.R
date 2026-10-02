# E39_prepare_and_gate.R
# =============================================================================
# PREPARE THE GSE25066 ANALYSIS FRAME, SCORE PROLIFERATION, AND RUN THE
# IDENTIFIABILITY GATE. IT FITS NOTHING.
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md,
# COMMITTED AT 536c29f together with the amendment that resolved its section 3 /
# section 4 ambiguity toward three cohorts. Read it first; this header restates
# its rules and invents none.
#
# =============================================================================
# THE SCRIPT STOPS AT THE GATE. THAT IS THE POINT.
# =============================================================================
# Declaration section 3 makes the identifiability gate a STOP-AND-CHECK: its
# output is read by the author BEFORE E40 is written, because rho(OXPHOS,
# PROLIF_DISJOINT) decides whether the ladder is separable at all.
#
#   >= 0.80   the ladder is not identifiable in that cohort. m0 only there.
#   0.60-0.80 report the ladder with VIFs beside every coefficient, and state
#             that the rungs are not cleanly separable.
#   < 0.60    report the ladder as specified.
#
# The band is read PER COHORT. A cohort may be identifiable where another is not
# and the consequence then applies to that cohort alone.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits NO outcome model. No glm, no coxph, no survfit. E40 does that.
#   - It does NOT rescore OXPHOS. OX is reused from the H4 frames so the
#     exposure is byte-identical to the one H4 used (declaration section 4).
#   - It does NOT re-standardise OX, which already arrives at sd = 1 from
#     script 13's .z() (declaration section 4.3). Only PROLIF is standardised.
#   - It does NOT use mitoPPS, which does not exist in these cohorts
#     (amendment A1). THIS ARM IS SINGLE-INSTRUMENT THROUGHOUT and must never
#     be written as though it had cleared the two-instrument bar.
#   - It does NOT touch METABRIC. That is declared in, and is a later phase.
#   - It writes NOTHING to myc_human_validation, which is FROZEN at d3ac60e.
#   - No per-gene FDR anywhere.
#
# =============================================================================
# SCALE - CLAUDE.md
# =============================================================================
# Every neoadjuvant expression matrix is log2 (script 12). GSVA therefore takes
# kcdf = "Gaussian". The call below matches scripts/13_outcome_models.R EXACTLY
# and NOT E02's - the two differ, and script 13 is the one that produced the OX
# and MYC already sitting in the H4 frames.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(GSVA)
})

message("\nE39: prepare the frame, score proliferation, run the gate\n",
        strrep("=", 78))
message("IT FITS NOTHING. The gate is a stop-and-check read before E40.")

# READ-ONLY inputs from the FROZEN repo. Nothing is ever written there.
PATH_NEO <- paste0("/Users/gs/code/myc_human_validation/results/",
                   "neoadjuvant_cohorts.rds")
PATH_H4  <- paste0("/Users/gs/code/myc_human_validation/results/",
                   "h4_outcome_models.rds")

PATH_E39      <- file.path(DIR_RESULTS, "e39_prepare_and_gate.rds")
PATH_E39_GATE <- file.path(DIR_TABLES,  "E39_identifiability_gate.csv")

# =============================================================================
# 0. CONSTANTS - every one from the declaration, none invented here
# =============================================================================
COHORTS    <- c("GSE25066", "GSE194040", "GSE164458")
PRIMARY    <- "GSE25066"
PROLIF_SET <- "PROLIF_DISJOINT"

# script 13's own constants, reproduced so the call cannot drift.
GSVA_MIN_SET     <- 3L
NA_GENE_MAX_FRAC <- 0.05        # A5

# The DRFS columns. "even" is GEO's own misspelling of "event" and is kept
# verbatim: renaming it here would hide which upstream column was read.
COL_DRFS_EVENT <- "drfs_1_event_0_censored"
COL_DRFS_TIME  <- "drfs_even_time_years"
# The column NAME says 1 = event and 0 = censored. That is asserted against the
# data, not trusted. A silent reversal here would invert every hazard ratio in
# E40 with no error thrown.
DRFS_EVENT_STR    <- "1"
DRFS_CENSORED_STR <- "0"

# Declared counts. Every one is a STOP, not a warning.
N_ALL        <- 508L
N_EVENTS_ALL <- 111L
MED_T_ALL    <- 2.7377
MAX_T_ALL    <- 7.439
TOL_T        <- 1e-3
N_470        <- 470L
N_EVENTS_470 <- 104L
N_HER2_DROP  <- 5L              # declaration 5.2

SUBTYPE_KEEP <- c("HRpos_HER2neg", "TNBC")
SUBTYPE_DROP <- c("HRneg_HER2pos", "HRpos_HER2pos")

# Declaration section 7, the ER counts in the 470.
ER_PRIMARY   <- "er_status_ihc"
ER_SENSITIVE <- "esr1_status"
ER_470 <- tibble::tribble(
  ~which,         ~level, ~n,   ~events,
  "er_status_ihc", "P",   279L, 40L,
  "er_status_ihc", "N",   191L, 64L,
  "esr1_status",   "P",   269L, 37L,
  "esr1_status",   "N",   201L, 67L)

# Declaration section 4, the recovery figures.
RECOVERY <- c(GSE25066 = 293L, GSE194040 = 305L, GSE164458 = 311L)

# Declaration section 3, the bands.
GATE_BANDS <- c(not_identifiable = 0.80, not_separable = 0.60)
.band <- function(r) {
  a <- abs(r)
  if (is.na(a)) NA_character_
  else if (a >= GATE_BANDS[["not_identifiable"]]) ">= 0.80 NOT IDENTIFIABLE"
  else if (a >= GATE_BANDS[["not_separable"]])    "0.60-0.80 not cleanly separable"
  else                                            "< 0.60 ladder as specified"
}

.z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs (both READ-ONLY, from the frozen repo)")

for (p in c(PATH_NEO, PATH_H4)) {
  .stop_if(file.exists(p), "absent: ", p,
           ". myc_human_validation must be present and is read-only.")
}
neo <- readRDS(PATH_NEO)
h4  <- readRDS(PATH_H4)
sd_ <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))

.stop_if(all(COHORTS %in% names(neo$cohorts)), "a cohort is missing from ",
         basename(PATH_NEO))
.stop_if(all(COHORTS %in% names(h4$frames)), "a cohort is missing from ",
         basename(PATH_H4), "$frames")

PD <- sd_$cov_sets[[PROLIF_SET]]
.stop_if(length(PD) == 318L, PROLIF_SET, " is ", length(PD),
         " genes, not 318. The on-disk set differs from the declaration.")
message("   ", PROLIF_SET, ": ", length(PD), " genes")
message("   OX is reused from the H4 frames and is NOT rescored.")

# =============================================================================
# 2. The GSE25066 analysis frame
# =============================================================================
message("\n2. the GSE25066 analysis frame")

ph <- neo$cohorts[[PRIMARY]]$pheno
fr <- h4$frames[[PRIMARY]]
.stop_if(nrow(ph) == N_ALL && nrow(fr) == N_ALL,
         "expected ", N_ALL, " rows in both pheno and frame; got ",
         nrow(ph), " and ", nrow(fr))

# --- 2.1 the join -----------------------------------------------------------
matched <- sum(fr$sample_id %in% ph$sample_id)
.stop_if(matched == N_ALL, "join on sample_id matched ", matched, " of ",
         N_ALL, "; the key is not what it is taken to be")
message("   join on sample_id: ", matched, " of ", N_ALL)

# --- 2.2 DRFS coercion, with the direction ASSERTED -------------------------
# The column name says 1 = event, 0 = censored. Checked against the data before
# anything is coerced: as.integer() on an unexpected encoding would silently
# invert every hazard ratio in E40 and throw no error.
.stop_if(all(c(COL_DRFS_EVENT, COL_DRFS_TIME) %in% names(ph)),
         "a DRFS column is absent from pheno")
ev_raw <- ph[[COL_DRFS_EVENT]]
.stop_if(is.character(ev_raw), COL_DRFS_EVENT, " is ", class(ev_raw)[1],
         ", not character. The encoding has changed; stop and look.")
ev_levels <- sort(unique(ev_raw))
.stop_if(identical(ev_levels, sort(c(DRFS_EVENT_STR, DRFS_CENSORED_STR))),
         COL_DRFS_EVENT, " has values {", paste(ev_levels, collapse = ", "),
         "}, not {", DRFS_CENSORED_STR, ", ", DRFS_EVENT_STR,
         "}. The event encoding is NOT what the column name says. STOP.")
n_ev_str  <- sum(ev_raw == DRFS_EVENT_STR)
n_cen_str <- sum(ev_raw == DRFS_CENSORED_STR)
.stop_if(n_ev_str == N_EVENTS_ALL,
         "'", DRFS_EVENT_STR, "' occurs ", n_ev_str, " times, not ",
         N_EVENTS_ALL, ". If the encoding were reversed this would be ",
         n_cen_str, ". STOP - do not guess which way round it is.")
message("   DRFS encoding asserted: '", DRFS_EVENT_STR, "' = event (",
        n_ev_str, "), '", DRFS_CENSORED_STR, "' = censored (", n_cen_str, ")")

tm_raw <- ph[[COL_DRFS_TIME]]
.stop_if(is.character(tm_raw), COL_DRFS_TIME, " is not character")
tm <- suppressWarnings(as.numeric(tm_raw))
.stop_if(!anyNA(tm), "coercing ", COL_DRFS_TIME,
         " to numeric produced ", sum(is.na(tm)), " NA")

D <- ph %>%
  dplyr::transmute(
    sample_id,
    drfs_event = as.integer(.data[[COL_DRFS_EVENT]] == DRFS_EVENT_STR),
    drfs_time  = suppressWarnings(as.numeric(.data[[COL_DRFS_TIME]])),
    er_primary   = as.character(.data[[ER_PRIMARY]]),
    er_sensitive = as.character(.data[[ER_SENSITIVE]])) %>%
  dplyr::inner_join(fr, by = "sample_id")
.stop_if(nrow(D) == N_ALL, "the joined frame is ", nrow(D), " rows, not ", N_ALL)

.stop_if(sum(stats::complete.cases(D$drfs_time, D$drfs_event)) == N_ALL,
         "DRFS is not complete for all ", N_ALL)
.stop_if(sum(D$drfs_event) == N_EVENTS_ALL,
         "events after coercion: ", sum(D$drfs_event), ", not ", N_EVENTS_ALL)
med_t <- stats::median(D$drfs_time); max_t <- max(D$drfs_time)
.stop_if(abs(med_t - MED_T_ALL) < TOL_T, "median follow-up ", signif(med_t, 6),
         ", not ", MED_T_ALL)
.stop_if(abs(max_t - MAX_T_ALL) < TOL_T, "max follow-up ", signif(max_t, 6),
         ", not ", MAX_T_ALL)
message("   DRFS: ", N_ALL, " complete, ", sum(D$drfs_event),
        " events, median ", round(med_t, 4), " y, max ", round(max_t, 4), " y")

# --- 2.3 the 470 ------------------------------------------------------------
MODEL_VARS <- c("pcr", "MYC", "OX", "BUFFER_c", "subtype")
D$in_470 <- stats::complete.cases(D[, MODEL_VARS])
.stop_if(sum(D$in_470) == N_470, "the pCR model set is ", sum(D$in_470),
         ", not ", N_470)
.stop_if(sum(D$drfs_event[D$in_470]) == N_EVENTS_470,
         "events in the 470: ", sum(D$drfs_event[D$in_470]), ", not ",
         N_EVENTS_470)
message("   the 470: ", sum(D$in_470), " complete on ",
        paste(MODEL_VARS, collapse = " + "), ", ",
        sum(D$drfs_event[D$in_470]), " events")

# --- 2.4 ER counts in the 470, asserted against declaration section 7 --------
er_obs <- dplyr::bind_rows(lapply(c(ER_PRIMARY, ER_SENSITIVE), function(w) {
  col <- if (w == ER_PRIMARY) "er_primary" else "er_sensitive"
  D[D$in_470, ] %>% dplyr::group_by(which = w, level = .data[[col]]) %>%
    dplyr::summarise(n = dplyr::n(), events = sum(drfs_event), .groups = "drop")
}))
chk <- ER_470 %>% dplyr::left_join(er_obs, by = c("which", "level"),
                                   suffix = c("_dec", "_obs"))
.stop_if(all(chk$n_dec == chk$n_obs) && all(chk$events_dec == chk$events_obs),
         "ER counts in the 470 do not match declaration section 7")
message("   ER counts in the 470 match declaration section 7 (8 of 8 cells)")

# --- 2.5 collapse subtype, dropping the HER2-positive patients --------------
# Declaration 5.2: at 104 events a factor level of n = 1 is not fittable. The 5
# are DROPPED, not pooled, and the drop is reported.
sub_470 <- D[D$in_470, ] %>% dplyr::count(subtype, name = "n")
n_drop <- sum(sub_470$n[sub_470$subtype %in% SUBTYPE_DROP])
.stop_if(n_drop == N_HER2_DROP, "HER2-positive patients in the 470: ", n_drop,
         ", not ", N_HER2_DROP)
D$in_model <- D$in_470 & D$subtype %in% SUBTYPE_KEEP
D$subtype2 <- factor(ifelse(D$in_model, as.character(D$subtype), NA_character_),
                     levels = SUBTYPE_KEEP)
n_model <- sum(D$in_model)
message("   subtype collapsed to ", paste(SUBTYPE_KEEP, collapse = " vs "),
        "; DROPPED ", n_drop, " HER2-positive (",
        paste(sprintf("%s n=%d", sub_470$subtype[sub_470$subtype %in% SUBTYPE_DROP],
                      sub_470$n[sub_470$subtype %in% SUBTYPE_DROP]),
              collapse = ", "), ")")
message("   analysis frame after the drop: ", n_model, " patients, ",
        sum(D$drfs_event[D$in_model]), " events")

SET_SIZES <- tibble::tibble(
  set = c("all", "pCR model set (the 470)", "after HER2 drop"),
  n = c(N_ALL, sum(D$in_470), n_model),
  events = c(sum(D$drfs_event), sum(D$drfs_event[D$in_470]),
             sum(D$drfs_event[D$in_model])))

# =============================================================================
# 3. Score PROLIF_DISJOINT into all three cohorts
# =============================================================================
# Declaration section 4 and 4.4. The call matches scripts/13_outcome_models.R
# EXACTLY: .impute() first, universe pinned, kcdf Gaussian, minSize 3L, and a
# silently dropped set is a hard stop.
message("\n3. scoring ", PROLIF_SET, " into all three cohorts")

# script 13's A5 rule, reproduced. Called even where it is a no-op so that the
# provenance is identical and the script does not rely on a property of one
# cohort - it is a no-op in GSE25066 (0 NA) and may not be in the other two.
.impute <- function(M, label) {
  na_g <- rowSums(is.na(M))
  drop <- na_g > NA_GENE_MAX_FRAC * ncol(M)
  if (any(drop)) {
    message("      ", label, ": dropped ", sum(drop), " gene(s) above ",
            100 * NA_GENE_MAX_FRAC, "% NA")
    M <- M[!drop, , drop = FALSE]
  }
  n_imp <- sum(is.na(M))
  if (n_imp) {
    med <- apply(M, 1L, stats::median, na.rm = TRUE)
    idx <- which(is.na(M), arr.ind = TRUE)
    M[idx] <- med[idx[, "row"]]
    message("      ", label, ": imputed ", n_imp, " value(s) at the gene median")
  }
  .stop_if(!anyNA(M), label, ": NA survived .impute()")
  M
}

.score_prolif <- function(cn) {
  M  <- .impute(neo$cohorts[[cn]]$expr, cn)
  ph <- neo$cohorts[[cn]]$pheno
  .stop_if(identical(colnames(M), ph$sample_id),
           cn, ": expr colnames are not pheno$sample_id. The join key differs ",
           "from script 13's and the scores would be misaligned.")
  present <- intersect(PD, rownames(M))
  .stop_if(length(present) == RECOVERY[[cn]],
           cn, ": ", PROLIF_SET, " recovery is ", length(present), " of ",
           length(PD), ", not the declared ", RECOVERY[[cn]],
           ". The on-disk set differs from the audit; find the cause before ",
           "scoring proceeds.")
  .stop_if(length(present) >= GSVA_MIN_SET, cn, ": below the GSVA floor")

  sym  <- rownames(M)
  sets <- list(present)
  names(sets) <- PROLIF_SET
  sets[[".PIN_A"]] <- sym[c(TRUE, FALSE)]
  sets[[".PIN_B"]] <- sym[c(FALSE, TRUE)]

  par <- GSVA::gsvaParam(exprData = M, geneSets = sets, kcdf = "Gaussian",
                         minSize = GSVA_MIN_SET, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)
  .stop_if(PROLIF_SET %in% rownames(s),
           cn, ": GSVA silently dropped ", PROLIF_SET)
  message("      ", sprintf("%-10s %3d of %d genes (%.1f%%)", cn,
          length(present), length(PD), 100 * length(present) / length(PD)))
  tibble::tibble(cohort = cn, sample_id = colnames(s),
                 PROLIF_raw = as.numeric(s[PROLIF_SET, ]),
                 n_genes = length(present))
}

PS <- dplyr::bind_rows(lapply(COHORTS, .score_prolif))

# --- standardise PROLIF within cohort; OX is NOT touched --------------------
# Declaration 4.3: OX arrives at sd = 1 from script 13's .z(). Re-standardising
# would be arithmetically a no-op and would break the byte-identity section 4
# requires, so it is left exactly as it arrives.
PS <- PS %>% dplyr::group_by(cohort) %>%
  dplyr::mutate(PROLIF = .z(PROLIF_raw)) %>% dplyr::ungroup()

OX_AUDIT <- dplyr::bind_rows(lapply(COHORTS, function(cn) {
  v <- h4$frames[[cn]]$OX
  tibble::tibble(cohort = cn, n = length(v), ox_mean = mean(v), ox_sd = sd(v))
}))
.stop_if(all(abs(OX_AUDIT$ox_sd - 1) < 1e-8),
         "OX does not arrive standardised in every cohort; declaration 4.3 ",
         "assumes it does")
message("   OX ARRIVED STANDARDISED (sd = 1 in all three) and was NOT ",
        "re-standardised:")
OX_AUDIT %>%
  dplyr::mutate(ox_mean = signif(ox_mean, 3), ox_sd = round(ox_sd, 6)) %>%
  as.data.frame() %>% print(row.names = FALSE)

RECOVERY_TAB <- PS %>% dplyr::distinct(cohort, n_genes) %>%
  dplyr::mutate(n_declared = length(PD),
                pct = round(100 * n_genes / n_declared, 1))

# =============================================================================
# 4. THE IDENTIFIABILITY GATE - declaration section 3
# =============================================================================
message("\n4. the identifiability gate, per cohort, all three")

GATE <- dplyr::bind_rows(lapply(COHORTS, function(cn) {
  d <- h4$frames[[cn]] %>%
    dplyr::select(sample_id, OX, MYC) %>%
    dplyr::inner_join(PS[PS$cohort == cn, c("sample_id", "PROLIF")],
                      by = "sample_id")
  .stop_if(nrow(d) == nrow(h4$frames[[cn]]),
           cn, ": the PROLIF join lost rows (", nrow(d), " of ",
           nrow(h4$frames[[cn]]), ")")
  .rho <- function(a, b) stats::cor(d[[a]], d[[b]], method = "spearman",
                                    use = "complete.obs")
  tibble::tibble(
    cohort = cn, n = nrow(d),
    pair = c("OX vs PROLIF_DISJOINT", "OX vs MYC", "MYC vs PROLIF_DISJOINT"),
    rho  = c(.rho("OX", "PROLIF"), .rho("OX", "MYC"), .rho("MYC", "PROLIF")))
})) %>% dplyr::rowwise() %>% dplyr::mutate(band = .band(rho)) %>%
  dplyr::ungroup()

GATE %>% dplyr::mutate(rho = round(rho, 4)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   THE BAND THAT DECIDES E40 is rho(OX, PROLIF_DISJOINT), per cohort:")
GATE %>% dplyr::filter(pair == "OX vs PROLIF_DISJOINT") %>%
  dplyr::transmute(cohort, n, rho = round(rho, 4), band) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 5. Save - every number reported above comes from one of these objects
# =============================================================================
message("\n5. saving")

readr::write_csv(GATE, PATH_E39_GATE)
saveRDS(list(
  frame        = D,
  set_sizes    = SET_SIZES,
  er_counts    = er_obs,
  er_declared  = ER_470,
  subtype_470  = sub_470,
  her2_dropped = n_drop,
  prolif       = PS,
  recovery     = RECOVERY_TAB,
  ox_audit     = OX_AUDIT,
  gate         = GATE,
  gate_bands   = GATE_BANDS,
  spec = list(
    declaration = paste0("docs/2026-10-02_E39_respiratory_axis_decomposition",
                         "_declaration.md at 536c29f"),
    posture = paste0("EXPLORATORY and POST-HOC. A descriptive decomposition, ",
                     "not a fifth hypothesis. Nothing here carries a pass, a ",
                     "fail or a falsification criterion."),
    fits_nothing = paste0("E39 fits no outcome model. The gate is a ",
                          "stop-and-check read before E40 is written."),
    instrument = paste0("SINGLE-INSTRUMENT throughout. mitoPPS does not exist ",
                        "in these cohorts (amendment A1); this has not ",
                        "cleared the arm's two-instrument bar and must never ",
                        "be written as though it had."),
    oxphos = paste0("OX reused from the H4 frames, NOT rescored and NOT ",
                    "re-standardised - it arrives at sd = 1 from script 13's ",
                    ".z(). Only PROLIF_DISJOINT is standardised."),
    gsva = paste0("scripts/13_outcome_models.R's call, not E02's: ",
                  ".impute() then gsvaParam(kcdf = 'Gaussian', minSize = 3L, ",
                  "maxSize = Inf) on log2 input, universe pinned with ",
                  ".PIN_A/.PIN_B, a silently dropped set a hard stop."),
    coverage_caveat = paste0("GSE25066 carries 57 of 89 OXPHOS subunits ",
                             "(64.0%) on GPL96. NO coefficient from it may be ",
                             "compared numerically with a TCGA or SCAN-B one ",
                             "(declaration 4.1)."),
    drfs_encoding = paste0("'", DRFS_EVENT_STR, "' = event, '",
                           DRFS_CENSORED_STR, "' = censored, ASSERTED against ",
                           "the data and not inferred from the column name."),
    frozen = "Nothing written to myc_human_validation (d3ac60e)."),
  built = Sys.time()), PATH_E39)

message("   ", PATH_E39, "\n   ", PATH_E39_GATE)
message("\nE39 done. NO MODEL WAS FITTED. Read the gate before E40.\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E39)

  # THE GATE. This is what decides whether E40's ladder is reportable past m0.
  x$gate %>% dplyr::mutate(rho = round(rho, 4)) %>% as.data.frame()
  x$gate %>% dplyr::filter(pair == "OX vs PROLIF_DISJOINT") %>% as.data.frame()

  # The analysis sets, and the HER2 drop.
  x$set_sizes %>% as.data.frame()
  x$subtype_470 %>% as.data.frame(); x$her2_dropped

  # ER, observed against declaration section 7.
  x$er_counts %>% as.data.frame()
  x$er_declared %>% as.data.frame()

  # Proliferation: recovery, and that OX was left alone.
  x$recovery %>% as.data.frame()
  x$ox_audit %>% as.data.frame()
  x$prolif %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(n = dplyr::n(), mean = mean(PROLIF),
                     sd = sd(PROLIF), .groups = "drop") %>% as.data.frame()

  # The frame itself.
  utils::head(x$frame) %>% print()
  x$frame %>% dplyr::count(in_470, in_model, name = "n") %>% as.data.frame()

  utils::str(x$spec)
}
