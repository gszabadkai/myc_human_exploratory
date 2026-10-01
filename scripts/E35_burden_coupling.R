# E35_burden_coupling.R
# =============================================================================
# IS THE OXPHOS-TO-CONFIGURATION COUPLING GRADED BY GENOMIC BURDEN?
#
# EXPLORATORY. Nothing here is pre-registered.
#
# Declaration: docs/2026-10-01_E35_declaration.md, COMMITTED ALONE AT 74507b1
# before any burden variable was opened. Data note, with the two feasibility
# constraints the declaration did not anticipate and every decision it left
# open: docs/2026-10-01_e35_data.md at f7f19e3. Read both first; this header
# restates their rules, it does not invent them.
#
# =============================================================================
# THE QUESTION, AND WHAT IT IS NOT
# =============================================================================
# In the mouse pre-tumour gland the respiratory state and the apoptotic
# configuration are unrelated in wild-type animals (slope -0.012, p 0.95) and
# coupled only where MYC is expressed. E34 found the coupling in human tumours
# present irrespective of MYC level. The leading reconciliation is that the
# mouse contrasts a pre-tumour gland with transformed tissue whilst the human
# contrasts tumours with tumours, so the human comparison may lie entirely above
# the condition the mouse manipulates. This tests that from inside the tumour
# cohort, using genomic burden as a graded proxy for drive.
#
# THE ENDPOINT IS THE COUPLING, NOT KILLING. There is no death readout in TCGA.
# The coupling is a transcript relationship: the two OXPHOS contrasts on the
# configuration composite. EVERY SENTENCE WRITTEN FROM THIS SAYS "COUPLING",
# NEVER "KILLING" AND NEVER "LETHALITY". Killing is the mouse's variable and
# stays there.
#
#   Q2 - Q1   the coupling where MYC is low
#   Q4 - Q3   the coupling where MYC is high
#
# =============================================================================
# STRATIFY, DO NOT ADJUST. NO PRODUCT TERM ANYWHERE.
# =============================================================================
# Burden correlates with proliferation (+0.375 to +0.568) and with subtype, so
# the circularity this analysis exists to escape partly follows it in. The two
# contrasts are therefore computed WITHIN each burden stratum. There is no
# continuous interaction term and no burden x quadrant product term: that
# estimand has failed in Block B, Block C, Block G, H4, N1 and N2.
#
# =============================================================================
# THE THREE PROXIES - ALL REPORTED, NONE PROMOTED AFTER THE FACT
# =============================================================================
#   Aneuploidy.Score  PRIMARY, tertiles within the analysis set
#   FGA               tertiles - the same question on the same cohort
#   genome_doublings  ITS OWN THREE NATIVE LEVELS, NOT TERTILES (data note 1.1);
#                     level 2 is NOT READABLE (min quadrant cell 4) and is
#                     reported as such, never estimated
# They correlate 0.50 to 0.80 with each other. They are replications, NOT
# independent evidence, and the result note says so.
#
# =============================================================================
# MODELS - all reported, never chosen between
# =============================================================================
#   m1   configuration ~ quadrant + PAM50 + purity + leukocyte fraction
#   m3   m1 + PROLIF_DISJOINT                  mandatory companion
#   m2   m1 without PAM50, within each PAM50 subtype - THE DECLARED RESULT, and
#        readable in only 1 of 15 cells. Reported honestly as that.
#   m2L  the same within COMBINED LUMINAL (LumA + LumB) - the companion the
#        author added on 2026-10-01, before this script, because m2 as declared
#        cannot carry the result. Readable in 3 of 3 strata on both continuous
#        proxies. functions/strata.R already defines this stratum for exactly
#        this underpowering. IT DOES NOT REPLACE m2; both are reported.
#
# READABILITY FLOOR: n >= 20 in all four quadrant cells. Below it a cell is
# REPORTED AS NOT READABLE, never estimated.
#
# =============================================================================
# THE READING RULE (declaration section 2; data note section 6)
# =============================================================================
# Applied to the CONFIGURATION COMPOSITE on the PRIMARY proxy:
#   GRADED               both OXPHOS contrasts increase monotonically across
#                        strata in the same direction, AND both low/high CIs
#                        do not overlap
#   SATURATED            both present (CI excludes zero) in every stratum, and
#                        low/high CIs overlap
#   ABSENT AT LOW BURDEN both fail to clear zero in the lowest stratum and
#                        clear it in the highest
#   NONE OF THE THREE    anything else, reported literally, no reading invented
# SATURATED IS THE DECLARED EXPECTED OUTCOME and is written as the expected
# result - not as a disappointment and not as a failure.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - TCGA ONLY. SCAN-B has no CNV, so NO REPLICATION COHORT EXISTS.
#   - NO outcome or survival variable enters. NO causal language.
#   - Burden is a PROXY for drive, not a measure of it.
#   - Configuration is NOT lethality.
#   - NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM.
#   - Nothing is written to myc_human_validation or myc_mouse.
#   - No per-gene FDR. No product term.
#
# N3: transcript scores. Nothing is "primed".
# SPECIES: human throughout. NO ORTHOLOG FUNCTION IS CALLED.
# SCALE: genes on log2(linear DESeq2-normalised + 1); the quadrant's inputs are
# READ, never re-scored. GSVA is cohort-relative; one cohort here, nothing
# pooled, and no value is compared numerically with any other cohort.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE35: is the OXPHOS-to-configuration coupling graded by genomic burden?\n",
        strrep("=", 78))

PATH_E35         <- file.path(DIR_RESULTS, "e35_burden_coupling.rds")
PATH_E35_CSV     <- file.path(DIR_TABLES,  "E35_coupling_by_burden.csv")
PATH_E35_FIG     <- file.path(DIR_FIGURES, "E35_coupling_by_burden.png")
PATH_E35_FIG_PDF <- file.path(DIR_FIGURES, "E35_coupling_by_burden.pdf")
DIR_DOC_FIGURES      <- here::here("docs", "figures")
PATH_E35_FIG_DOC     <- file.path(DIR_DOC_FIGURES,
                                  "2026-10-01_E35_coupling_by_burden.png")
PATH_E35_FIG_DOC_PDF <- file.path(DIR_DOC_FIGURES,
                                  "2026-10-01_E35_coupling_by_burden.pdf")

# =============================================================================
# 0. Constants and the reading rule
# =============================================================================
# Read-only, from the FROZEN repo. Nothing is written there, ever.
PATH_STATE <- "/Users/gs/code/myc_human_validation/results/state_definition.rds"

ARM_OX      <- "OXPHOS subunits"     # the narrow set, NOT the umbrella
PROLIF_COV  <- "PROLIF_DISJOINT"
MIN_CELL    <- 20L                   # the declared floor

# E33's arms and signs, as E34 adopted them. E11 defines NO composite - the
# declaration's wording is inherited and the divergence is recorded in both
# data notes.
TRIGGER  <- c("BBC3", "BID", "BIK", "BAD")
GUARDIAN <- c("BCL2L1", "MCL1")
CONFIG_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
TWELVE <- c(PRIMING_PRO, PRIMING_ANTI)
stopifnot(setequal(names(CONFIG_SIGN), c(TRIGGER, GUARDIAN)),
          all(names(CONFIG_SIGN) %in% TWELVE), !anyDuplicated(TWELVE))

READ_PRIMARY <- c(config_comp = "configuration composite")
READ_ARMS    <- c(arm_trigger = "trigger arm", arm_guardian = "guardian arm")
READ_GENES   <- stats::setNames(TWELVE, paste0("g_", TWELVE))
READOUTS     <- c(READ_PRIMARY, READ_ARMS, READ_GENES)

QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",
          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")
INSTR        <- c(gsva = "OX_gsva", mitopps = "OX_mitopps")
PAM50_LEVELS <- c("LumA", "LumB", "HER2", "Basal", "Normal")

# The burden proxies. `tertile` FALSE means stratify on the native levels.
PROXIES <- tibble::tribble(
  ~proxy,              ~label,              ~tertile, ~role,
  "Aneuploidy.Score",  "Aneuploidy score",  TRUE,     "PRIMARY",
  "FGA",               "Fraction genome altered", TRUE, "replication, same cohort",
  "genome_doublings",  "Genome doublings",  FALSE,    "replication, same cohort")

# The four contrasts, high minus low. THE READING USES THE OXPHOS PAIR ONLY.
CONTR <- tibble::tribble(
  ~contrast,  ~hi,     ~lo,     ~kind,    ~role,
  "Q2 - Q1",  QLEV[2], QLEV[1], "OXPHOS", "THE COUPLING where MYC is low",
  "Q4 - Q3",  QLEV[4], QLEV[3], "OXPHOS", "THE COUPLING where MYC is high",
  "Q4 - Q2",  QLEV[4], QLEV[2], "MYC",    "context only, not the reading",
  "Q3 - Q1",  QLEV[3], QLEV[1], "MYC",    "context only, not the reading")
COUPLING <- c("Q2 - Q1", "Q4 - Q3")

R_GRADED <- "GRADED"
R_SAT    <- "SATURATED"
R_ABSENT <- "ABSENT AT LOW BURDEN"
R_NONE   <- "NONE OF THE THREE"
READING_RULE <- paste0(
  "On the configuration composite, both OXPHOS contrasts (Q2-Q1, Q4-Q3): ",
  "GRADED iff both increase monotonically across strata in the same direction ",
  "AND both lowest/highest CIs do not overlap; SATURATED iff both have CIs ",
  "excluding zero in every stratum AND both lowest/highest CIs overlap; ",
  "ABSENT AT LOW BURDEN iff both fail to clear zero in the lowest stratum and ",
  "clear it in the highest; otherwise NONE OF THE THREE, reported literally. ",
  "SATURATED is the declared expected outcome.")

message("\n0. the rule, fixed at 74507b1 / f7f19e3 before any estimate:\n   ",
        READING_RULE)

# =============================================================================
# 1. The frozen constructor, read-only
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
message("   integrity contract holds; environment is baseenv()")

# =============================================================================
# 2. Inputs - TCGA only
# =============================================================================
message("\n2. inputs (TCGA only; SCAN-B has no CNV, so no replication cohort)")

frames <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
mito   <- readRDS(PATH_TCGA_MITO)
nw     <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
cv     <- readRDS(PATH_TCGA_COV)$covariates
ID_T   <- colnames(mito$gsva_arms)
stopifnot(length(ID_T) == EXPECT_TCGA_SAMPLES)

# --- 2.1 the twelve genes, log2(linear DESeq2-normalised + 1) ----------------
lin <- readRDS(PATH_TCGA_LINEAR)
if (!identical(lin$scale, "linear_deseq2_normalised")) {
  stop("the linear matrix is not on the linear DESeq2-normalised scale",
       call. = FALSE)
}
G  <- log2(lin$mat[, ID_T, drop = FALSE] + 1)
gr <- .gene_rows(TWELVE, G, .symbol_resolver(rownames(G), NULL))
if (length(gr$missing)) {
  stop("did not resolve: ", paste(gr$missing, collapse = ", "), call. = FALSE)
}
GM <- gr$mat[TWELVE, , drop = FALSE]
rm(lin, G); invisible(gc(verbose = FALSE))

# --- 2.2 E33's composites, as E34 built them. z across the WHOLE cohort ------
gv <- apply(GM[names(CONFIG_SIGN), , drop = FALSE], 1L, stats::var)
if (any(!(gv > 0))) stop("zero-variance gene(s) in the composite", call. = FALSE)
ZG <- t(scale(t(GM[names(CONFIG_SIGN), , drop = FALSE])))
.signed_mean <- function(genes) as.numeric(
  colMeans(sweep(ZG[genes, , drop = FALSE], 1L, CONFIG_SIGN[genes], "*")))
RAW <- list(config_comp  = .signed_mean(names(CONFIG_SIGN)),
            arm_trigger  = .signed_mean(TRIGGER),
            arm_guardian = .signed_mean(GUARDIAN))
for (g in TWELVE) RAW[[paste0("g_", g)]] <- as.numeric(GM[g, ])
message("   composite = mean z(BBC3, BID, BIK, BAD, BCL2L1, -MCL1); arms as E33")

# --- 2.3 covariates and the three burden proxies ----------------------------
fr <- frames[frames$cohort == "TCGA", ]
fr <- fr[match(ID_T, fr$sample_id), ]
if (anyNA(fr$sample_id)) stop("TCGA ids missing from frames.rds", call. = FALSE)
ci <- match(ID_T, cv$patient)
if (anyNA(ci)) stop("TCGA ids missing from the covariate table", call. = FALSE)
for (p in PROXIES$proxy) {
  if (!p %in% names(cv)) stop("burden proxy absent: ", p, call. = FALSE)
}
message("   burden proxies read: ",
        paste(sprintf("%s (NA %d)", PROXIES$proxy,
                      vapply(PROXIES$proxy, function(p) sum(is.na(cv[[p]][ci])),
                             integer(1))), collapse = ", "))

# =============================================================================
# 3. THE QUADRANT - E34's construction, copied verbatim
# =============================================================================
message("\n3. the quadrant")

S <- st$tcga$state
if (!setequal(S$patient, ID_T)) {
  stop("the frozen STATE's patients differ from this repo's TCGA ids",
       call. = FALSE)
}
myc_high <- S$STATE_gsva != LV[1]
if (!identical(myc_high, S$STATE_mitopps != LV[1])) {
  stop("the two instruments disagree on the MYC call", call. = FALSE)
}
QUAD <- tibble::tibble(sample_id = S$patient)
ASSERTS <- list()
for (inst in names(INSTR)) {
  ox <- S[[INSTR[[inst]]]]
  state_stored <- S[[paste0("STATE_", inst)]]
  rebuilt <- .build_state(S$MYC, ox, S$BUFFER_gistic)
  if (!identical(rebuilt, state_stored)) {
    stop("the constructor does not reproduce STATE_", inst, call. = FALSE)
  }
  # The constructor computes its OXPHOS call internally and does not return it.
  # Calling it again with myc := oxphos makes its MYC-high call BE its
  # OXPHOS-high call: same median, same `>` tie rule, same complete cases.
  probe <- .build_state(ifelse(is.na(S$MYC), NA_real_, ox), ox, S$BUFFER_gistic)
  if (any(probe == LV[2], na.rm = TRUE)) {
    stop("level 2 appeared in the myc := oxphos probe", call. = FALSE)
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
  # THE DECLARED ASSERTIONS - re-asserted here, and the script stops on failure.
  a3 <- identical(which(q == QLEV[3]), which(state_stored == LV[2]))
  a4 <- identical(which(q == QLEV[4]), which(state_stored %in% LV[3:4]))
  if (!a3 || !a4) {
    stop("ASSERTION FAILED (", inst, "): Q3 == level 2 is ", a3,
         ", Q4 == levels 3 + 4 is ", a4, ". STOP.", call. = FALSE)
  }
  ASSERTS[[inst]] <- tibble::tibble(instrument = inst, Q3_equals_level2 = a3,
                                    Q4_equals_levels3_4 = a4)
  QUAD[[paste0("quad_", inst)]] <- factor(q, levels = QLEV)
  message("   ", inst, ": Q3 == level 2 TRUE; Q4 == levels 3 + 4 TRUE")
}
ASSERTS <- dplyr::bind_rows(ASSERTS)

# =============================================================================
# 4. Blom normal scores, pooled once over the whole cohort
# =============================================================================
# E22's scorer. Scored ONCE over all 1,095; every stratum below SUBSETS these
# and never re-scores (E23's rule). Estimates are therefore in within-cohort SD
# units on a rank-based scale and are never compared with any other cohort.
message("\n4. Blom normal scores, once over the whole cohort")

.blom <- function(x) {
  n <- length(x)
  stopifnot(!anyNA(x))
  as.numeric(scale(stats::qnorm((rank(x, ties.method = "average") - 0.375) /
                                  (n + 0.25))))
}
.blom_na <- function(x) {
  out <- rep(NA_real_, length(x)); ok <- !is.na(x)
  if (any(ok)) out[ok] <- .blom(x[ok])
  out
}
# Everything is computed BEFORE the tibble: inside tibble() a new column shadows
# an outer object of the same name, the trap that stopped E32's first run.
qi <- match(ID_T, QUAD$sample_id)
DAT <- tibble::tibble(
  sample_id    = ID_T,
  quad_gsva    = QUAD$quad_gsva[qi],
  quad_mitopps = QUAD$quad_mitopps[qi],
  PAM50  = factor(as.character(fr$PAM50), levels = PAM50_LEVELS),
  purity = .blom_na(fr$purity),
  leuko  = .blom_na(fr$leuko),
  PROLIF = .blom(as.numeric(mito$gsva_cov[PROLIF_COV, ID_T])))
DAT$Luminal <- !is.na(DAT$PAM50) & DAT$PAM50 %in% c("LumA", "LumB")
for (p in PROXIES$proxy) DAT[[p]] <- as.numeric(cv[[p]][ci])
for (r in names(READOUTS)) DAT[[r]] <- .blom(RAW[[r]])
message("   PAM50 NA ", sum(is.na(DAT$PAM50)), "; purity and leuko both present ",
        sum(!is.na(DAT$purity) & !is.na(DAT$leuko)))

# =============================================================================
# 5. The strata - tertiles within the analysis set, or native levels
# =============================================================================
# Tertile cuts are computed on the ANALYSIS SET's complete cases for that
# proxy, within cohort. Aneuploidy.Score is an integer with 36 levels, so its
# tertiles are UNEQUAL because the cuts fall on ties. Reported, not corrected.
message("\n5. burden strata")

.base_rows <- function(dat, inst) {
  !is.na(dat[[paste0("quad_", inst)]]) & !is.na(dat$PAM50) &
    !is.na(dat$purity) & !is.na(dat$leuko)
}
.strata_for <- function(dat, proxy, use_tertile, rows) {
  v <- dat[[proxy]]
  ok <- rows & !is.na(v)
  out <- rep(NA_character_, nrow(dat))
  if (use_tertile) {
    cut <- unname(stats::quantile(v[ok], c(1/3, 2/3)))
    out[ok] <- as.character(cut(v[ok], c(-Inf, cut, Inf),
                                labels = c("T1", "T2", "T3")))
    lv <- c("T1", "T2", "T3")
  } else {
    out[ok] <- as.character(v[ok])
    lv <- sort(unique(out[ok]))
  }
  list(stratum = factor(out, levels = lv), levels = lv,
       cuts = if (use_tertile) cut else NA_real_)
}

# =============================================================================
# 6. The models - standardised means within each stratum. NO PRODUCT TERM.
# =============================================================================
message("\n6. models")

COVARS <- list(m1  = c("PAM50", "purity", "leuko"),
               m3  = c("PAM50", "purity", "leuko", "PROLIF"),
               m2  = c("purity", "leuko"),
               m2L = c("purity", "leuko"))
MODEL_LABEL <- c(m1  = "m1 (PAM50-adjusted)",
                 m3  = "m3 (m1 + PROLIF_DISJOINT)",
                 m2  = "m2 (within PAM50 subtype) - the declared result",
                 m2L = "m2L (within combined Luminal) - the companion")

# One fit. The group estimate for quadrant q is the model's prediction averaged
# over its OWN analysis rows with quadrant set to q - a standardised mean, a
# linear combination a'b whose 95% CI is a'Va on the residual df. The contrasts
# are differences of those combinations.
.fit_groups <- function(dat, y, covars) {
  cols <- c(y, "quad", covars)
  dd <- as.data.frame(dat[stats::complete.cases(dat[, cols]), cols, drop = FALSE])
  names(dd)[1] <- "yv"
  cell_n <- as.integer(table(factor(dd$quad, levels = QLEV)))
  names(cell_n) <- QLEV
  empty <- list(ok = FALSE, n = nrow(dd), cell_n = cell_n, df = NA_real_,
                A = NULL, b = NULL, V = NULL)
  dd$quad <- droplevels(factor(dd$quad, levels = QLEV))
  for (cvn in covars) if (is.factor(dd[[cvn]])) dd[[cvn]] <- droplevels(dd[[cvn]])
  if (nlevels(dd$quad) < 2L) return(empty)
  keep <- covars[vapply(covars, function(k)
    if (is.factor(dd[[k]])) nlevels(dd[[k]]) > 1L else
      stats::var(dd[[k]], na.rm = TRUE) > 0, logical(1))]
  f <- tryCatch(stats::lm(stats::reformulate(c("quad", keep), response = "yv"),
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

CONTRASTS <- list(); FITS <- list(); STRATA_INFO <- list(); k <- 0L
for (inst in names(INSTR)) {
  rows0 <- .base_rows(DAT, inst)
  for (pi in seq_len(nrow(PROXIES))) {
    proxy <- PROXIES$proxy[pi]
    sinfo <- .strata_for(DAT, proxy, PROXIES$tertile[pi], rows0)
    STRATA_INFO[[length(STRATA_INFO) + 1L]] <- tibble::tibble(
      instrument = inst, proxy = proxy, tertile = PROXIES$tertile[pi],
      levels = paste(sinfo$levels, collapse = "/"),
      cuts = paste(signif(sinfo$cuts, 4), collapse = ", "),
      n = sum(!is.na(sinfo$stratum)))
    d0 <- DAT
    d0$quad    <- d0[[paste0("quad_", inst)]]
    d0$stratum <- sinfo$stratum
    for (lev in sinfo$levels) {
      in_lev <- !is.na(d0$stratum) & d0$stratum == lev
      SETS <- list(list(model = "m1",  sub = "all",     rows = in_lev),
                   list(model = "m3",  sub = "all",     rows = in_lev),
                   list(model = "m2L", sub = "Luminal", rows = in_lev & d0$Luminal))
      for (s in PAM50_LEVELS) {
        SETS[[length(SETS) + 1L]] <- list(
          model = "m2", sub = s,
          rows = in_lev & !is.na(d0$PAM50) & d0$PAM50 == s)
      }
      for (S_ in SETS) {
        dsub <- d0[S_$rows, , drop = FALSE]
        for (r in names(READOUTS)) {
          fit <- .fit_groups(dsub, r, COVARS[[S_$model]])
          C <- t(vapply(seq_len(nrow(CONTR)), function(i) {
            if (!fit$ok) return(c(est = NA_real_, se = NA_real_,
                                  lo = NA_real_, hi = NA_real_))
            .lincomb(fit, fit$A[, CONTR$hi[i]] - fit$A[, CONTR$lo[i]])
          }, numeric(4)))
          is_readable <- fit$ok && min(fit$cell_n) >= MIN_CELL
          key <- tibble::tibble(
            instrument = inst, proxy = proxy, stratum = lev,
            model = S_$model, subset = S_$sub,
            readout = unname(READOUTS[[r]]), readout_col = r)
          k <- k + 1L
          CONTRASTS[[k]] <- dplyr::bind_cols(
            key[rep(1L, nrow(CONTR)), ],
            tibble::tibble(contrast = CONTR$contrast, kind = CONTR$kind,
                           est = C[, "est"], se = C[, "se"],
                           lo = C[, "lo"], hi = C[, "hi"],
                           ci_excludes_0 = (C[, "lo"] > 0) | (C[, "hi"] < 0),
                           readable = is_readable))
          FITS[[k]] <- dplyr::bind_cols(
            key, tibble::tibble(n = fit$n, df = fit$df, fitted = fit$ok,
                                min_cell = min(fit$cell_n),
                                readable = is_readable,
                                Q1 = fit$cell_n[QLEV[1]], Q2 = fit$cell_n[QLEV[2]],
                                Q3 = fit$cell_n[QLEV[3]], Q4 = fit$cell_n[QLEV[4]]))
        }
      }
    }
  }
}
CONTRASTS   <- dplyr::bind_rows(CONTRASTS)
FITS        <- dplyr::bind_rows(FITS)
STRATA_INFO <- dplyr::bind_rows(STRATA_INFO)
message("   ", nrow(FITS), " fits across ", nrow(PROXIES), " proxies x 2 ",
        "instruments x strata x {m1, m3, m2L, m2 x 5} x ", length(READOUTS),
        " readouts")

# =============================================================================
# 7. THE READING - configuration composite, primary proxy
# =============================================================================
message("\n7. the reading\n", strrep("-", 78))

# One (instrument, proxy, model, subset) block -> one reading. Needs every
# stratum of that block to be readable; otherwise it is not read at all.
.read_block <- function(tab, levs) {
  g <- function(lev, nm, col) {
    v <- tab[[col]][tab$stratum == lev & tab$contrast == nm]
    if (length(v) == 1L) v else NA_real_
  }
  if (!all(levs %in% tab$stratum)) return("not readable")
  if (!all(tab$readable)) return("not readable")
  lo_lev <- levs[1]; hi_lev <- levs[length(levs)]
  res <- vapply(COUPLING, function(nm) {
    e <- vapply(levs, function(l) g(l, nm, "est"), numeric(1))
    if (anyNA(e)) return(NA_character_)
    excl <- vapply(levs, function(l) {
      lo <- g(l, nm, "lo"); hi <- g(l, nm, "hi")
      isTRUE(lo > 0) || isTRUE(hi < 0)
    }, logical(1))
    overlap <- !(g(hi_lev, nm, "lo") > g(lo_lev, nm, "hi") ||
                 g(lo_lev, nm, "lo") > g(hi_lev, nm, "hi"))
    mono_up <- all(diff(e) > 0)
    if (mono_up && !overlap) return(R_GRADED)
    if (all(excl) && overlap) return(R_SAT)
    if (!excl[1] && excl[length(excl)]) return(R_ABSENT)
    R_NONE
  }, character(1))
  if (anyNA(res)) return("not readable")
  if (length(unique(res)) == 1L) unique(res) else R_NONE
}

READINGS <- CONTRASTS %>%
  dplyr::filter(readout_col == "config_comp", kind == "OXPHOS") %>%
  dplyr::group_by(instrument, proxy, model, subset) %>%
  dplyr::group_modify(~ {
    levs <- STRATA_INFO$levels[STRATA_INFO$instrument == .y$instrument &
                                 STRATA_INFO$proxy == .y$proxy][1]
    levs <- strsplit(levs, "/")[[1]]
    tibble::tibble(reading = .read_block(.x, levs),
                   strata_readable = sum(!duplicated(.x$stratum) & .x$readable),
                   strata_total = length(levs))
  }) %>%
  dplyr::ungroup()

PRIMARY_PROXY <- PROXIES$proxy[PROXIES$role == "PRIMARY"]
headline <- READINGS %>%
  dplyr::filter(proxy == PRIMARY_PROXY, model == "m1", instrument == "gsva") %>%
  dplyr::pull(reading)
reading <- if (length(headline) == 1L) headline else "not readable"

message("RULE: ", READING_RULE, "\n")
message("   the coupling by burden stratum, configuration composite, m1:")
CONTRASTS %>%
  dplyr::filter(readout_col == "config_comp", kind == "OXPHOS", model == "m1",
                subset == "all") %>%
  dplyr::transmute(instrument, proxy, stratum, contrast,
                   est = sprintf("%+.3f", est),
                   ci = sprintf("[%+.3f, %+.3f]", lo, hi), readable) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   readings, every block (the headline is ", PRIMARY_PROXY,
        ", m1, gsva):")
READINGS %>% as.data.frame() %>% print(row.names = FALSE)
message("\n   READING: ", reading)

reading_text <- switch(
  reading,
  "SATURATED" = paste0(
    "SATURATED - THE DECLARED EXPECTED OUTCOME, and it is written as the\n",
    "   expected result, not as a disappointment and not as a failure. The\n",
    "   coupling is present across the range of genomic burden represented in\n",
    "   this cohort, which is consistent with every established tumour lying\n",
    "   above the condition the mouse manipulates. It licenses ONE clause, and\n",
    "   it reinstates NO stratifier. Burden is a PROXY for drive, not a\n",
    "   measure of it; the configuration is NOT lethality; one cohort, and no\n",
    "   replication is possible because SCAN-B has no CNV."),
  "GRADED" = paste0(
    "GRADED - the coupling strengthens with genomic burden. Single cohort,\n",
    "   with the proxy caveat attached. It would make a burden-stratified\n",
    "   caution arguable rather than universal. Read the within-subtype rows\n",
    "   before believing it: burden is close to a subtype stratum here."),
  "ABSENT AT LOW BURDEN" = paste0(
    "ABSENT AT LOW BURDEN - the strongest of the three readings. The\n",
    "   coupling does not clear zero in the lowest burden stratum and does in\n",
    "   the highest, so some tumours lack the condition. Single cohort, proxy\n",
    "   exposure, and the subtype confound must be read before this is used."),
  "NONE OF THE THREE" = paste0(
    "NONE OF THE THREE declared readings fits. Reported literally, with no\n",
    "   reading invented for it. The rule is not relaxed after the fact."),
  paste0(
    "NOT READABLE on the primary block. The stratum cells did not clear the\n",
    "   declared floor of n >= ", MIN_CELL, ". Reported as such."))
message("\n", reading_text)

message("\n   the arms, m1, primary proxy (reported always, not the reading):")
CONTRASTS %>%
  dplyr::filter(proxy == PRIMARY_PROXY, model == "m1", subset == "all",
                instrument == "gsva", kind == "OXPHOS",
                readout_col %in% names(READ_ARMS)) %>%
  dplyr::transmute(readout, stratum, contrast, est = sprintf("%+.3f", est),
                   ci = sprintf("[%+.3f, %+.3f]", lo, hi)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   m2 as declared - readable blocks out of total:")
READINGS %>% dplyr::filter(model == "m2") %>%
  dplyr::summarise(readable_blocks = sum(reading != "not readable"),
                   total_blocks = dplyr::n()) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   m2L, the combined-Luminal companion:")
READINGS %>% dplyr::filter(model == "m2L") %>%
  dplyr::transmute(instrument, proxy, reading, strata_readable, strata_total) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. The figure
# =============================================================================
# The two OXPHOS contrasts with CIs against burden stratum, faceted by proxy,
# with the pooled m1 row and the combined-Luminal companion row. GSVA quadrant.
# NO line joins the strata: they are not exchangeable groups (data note 3).
message("\n8. figure")

FIG <- CONTRASTS %>%
  dplyr::filter(readout_col == "config_comp", kind == "OXPHOS",
                instrument == "gsva", model %in% c("m1", "m2L"),
                readable) %>%
  dplyr::left_join(PROXIES %>% dplyr::select(proxy, label), by = "proxy") %>%
  dplyr::mutate(
    panel = factor(label, levels = PROXIES$label),
    row = factor(ifelse(model == "m1", "All tumours (m1)",
                        "Luminal only (m2L)"),
                 levels = c("All tumours (m1)", "Luminal only (m2L)")),
    coupling = factor(ifelse(contrast == "Q2 - Q1",
                             "Q2 - Q1   coupling at MYC-low",
                             "Q4 - Q3   coupling at MYC-high"),
                      levels = c("Q2 - Q1   coupling at MYC-low",
                                 "Q4 - Q3   coupling at MYC-high")))
# The caption is wrapped by strwrap so its width is a property of the code
# rather than of my counting - the defect that took three passes in E34.
CAPTION <- paste(strwrap(paste0(
  "Standardised means from lm within each burden stratum: configuration ~ ",
  "quadrant + PAM50 + purity + leukocyte fraction (m1), and the same without ",
  "PAM50 within combined Luminal (m2L). The configuration is E33's signed ",
  "mean-z on log2(linear + 1). STRATIFIED, NOT ADJUSTED: there is no product ",
  "term and no continuous interaction anywhere, and no line joins the strata ",
  "because they are not exchangeable groups - burden tracks proliferation ",
  "(rho +0.375 to +0.568) and subtype, so the lowest tertile is three quarters ",
  "LumA. THE ENDPOINT IS THE COUPLING, NOT KILLING: there is no death readout ",
  "in TCGA. Burden is a PROXY for drive, not a measure of it, and the ",
  "configuration is not lethality. TCGA only; SCAN-B carries no CNV, so no ",
  "replication cohort exists. Genome doublings is shown on its three native ",
  "levels, not tertiles, and its level 2 is omitted as not readable (minimum ",
  "quadrant cell 4 against a declared floor of 20). Only strata clearing that ",
  "floor are drawn. NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM."),
  width = 190), collapse = "\n")

CCOL <- c("#0072B2", "#D55E00")
names(CCOL) <- levels(FIG$coupling)
DW <- 0.5
fig <- ggplot2::ggplot(FIG, ggplot2::aes(x = stratum, y = est, colour = coupling)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.16,
                         linewidth = 0.55,
                         position = ggplot2::position_dodge(width = DW)) +
  ggplot2::geom_point(size = 2.3,
                      position = ggplot2::position_dodge(width = DW)) +
  ggplot2::facet_grid(row ~ panel, scales = "free_x", space = "free_x") +
  ggplot2::scale_colour_manual(values = CCOL, name = NULL) +
  ggplot2::labs(
    x = "Genomic burden stratum (tertiles; genome doublings on its native levels)",
    y = "OXPHOS contrast on the configuration, within-cohort SD units, 95% CI",
    title = "Is the OXPHOS-to-configuration coupling graded by genomic burden?",
    subtitle = paste0(
      "Reading on the rule fixed at 74507b1 before any estimate: ", reading,
      if (identical(reading, R_SAT)) " - the declared expected outcome" else "",
      "\nTCGA only; SCAN-B has no CNV, so no replication cohort exists. ",
      "EXPLORATORY; not pre-registered."),
    caption = CAPTION, colour = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
    strip.text = ggplot2::element_text(size = 8),
    legend.position = "bottom", plot.title.position = "plot",
    plot.caption = ggplot2::element_text(size = 5.8, hjust = 0, lineheight = 1.1),
    plot.caption.position = "plot")

.ensure_dir(DIR_DOC_FIGURES)
for (p in c(PATH_E35_FIG, PATH_E35_FIG_DOC)) {
  ggplot2::ggsave(p, fig, width = 9.2, height = 6.4, dpi = 300)
}
for (p in c(PATH_E35_FIG_PDF, PATH_E35_FIG_DOC_PDF)) {
  ggplot2::ggsave(p, fig, width = 9.2, height = 6.4,
                  device = grDevices::cairo_pdf)
}
message("   ", PATH_E35_FIG, "\n   ", PATH_E35_FIG_DOC, "  (tracked)")

# =============================================================================
# 9. Save
# =============================================================================
message("\n9. saving")

out <- list(
  reading       = reading,
  reading_text  = reading_text,
  reading_rule  = READING_RULE,
  readings      = READINGS,
  contrasts     = CONTRASTS,
  fits          = FITS,
  strata        = STRATA_INFO,
  asserts       = ASSERTS,
  burden_summary = tibble::tibble(
    proxy = PROXIES$proxy, role = PROXIES$role, tertile = PROXIES$tertile,
    n_nonmissing = vapply(PROXIES$proxy,
                          function(p) sum(!is.na(DAT[[p]])), integer(1))),
  burden_cor = stats::cor(as.matrix(DAT[, PROXIES$proxy]),
                          use = "pairwise.complete.obs", method = "spearman"),
  burden_vs_prolif = vapply(PROXIES$proxy, function(p)
    stats::cor(DAT[[p]], DAT$PROLIF, method = "spearman",
               use = "complete.obs"), numeric(1)),
  counts = list(tcga_n = length(ID_T), pam50_na = sum(is.na(DAT$PAM50)),
                purity_ok = sum(!is.na(DAT$purity) & !is.na(DAT$leuko))),
  spec = list(
    declaration = "docs/2026-10-01_E35_declaration.md, committed ALONE at 74507b1",
    data_note   = "docs/2026-10-01_e35_data.md, f7f19e3",
    endpoint = paste("THE COUPLING, NOT KILLING. There is no death readout in",
                     "TCGA. Every sentence says coupling; never killing, never",
                     "lethality."),
    rule = READING_RULE,
    expected = "SATURATED is the declared expected outcome",
    no_product = paste("NO product term and no continuous interaction. The two",
                       "contrasts are computed WITHIN each burden stratum."),
    proxies = paste("Aneuploidy.Score primary; FGA and genome_doublings are",
                    "replications of the same question on the same cohort, NOT",
                    "independent evidence. They correlate 0.50 to 0.80."),
    gd = paste("genome_doublings has THREE values and cannot be tertiled; it is",
               "stratified on its native levels and level 2 is not readable",
               "(min quadrant cell 4). Data note 1.1."),
    m2 = paste("m2 within each PAM50 subtype is the DECLARED result and is",
               "readable in 1 of 15 cells. m2L, within combined Luminal, was",
               "added by the author before the script because m2 cannot carry",
               "it; it does NOT replace m2 and both are reported. Data note 1.2."),
    circularity = paste("burden correlates with PROLIF_DISJOINT at +0.375 to",
                        "+0.568 and tracks subtype so hard that a tertile is",
                        "close to a subtype stratum; the quadrant distribution",
                        "also shifts across burden, so the strata are NOT",
                        "exchangeable groups."),
    burden_is_proxy = "burden is a PROXY for drive, not a measure of it",
    not_lethality = "the configuration is NOT lethality",
    single_cohort = "TCGA only; SCAN-B has no CNV, so NO replication is possible",
    no_outcome = "NO outcome or survival variable enters at any point",
    no_treatment = paste("NOTHING HERE LICENSES A TREATMENT-SELECTION CLAIM and",
                         "nothing licenses any statement about killing."),
    no_fdr = "no per-gene FDR", n3 = "transcript scores; nothing is 'primed'"),
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

saveRDS(out, PATH_E35)
utils::write.csv(as.data.frame(CONTRASTS), PATH_E35_CSV, row.names = FALSE)
message("   ", PATH_E35, "\n   ", PATH_E35_CSV)
message("\nE35 done. READING: ", reading, "\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E35)

  cat(x$reading, "\n\n"); cat(x$reading_text, "\n")
  utils::str(x$counts)

  # The strata actually used, and the tertile cuts.
  x$strata %>% as.data.frame()

  # THE READING, every block. The headline is Aneuploidy.Score / m1 / gsva.
  x$readings %>% as.data.frame()

  # The coupling itself: both OXPHOS contrasts, every stratum, every proxy, m1.
  x$contrasts %>%
    dplyr::filter(readout_col == "config_comp", kind == "OXPHOS",
                  model == "m1", subset == "all") %>%
    dplyr::transmute(instrument, proxy, stratum, contrast,
                     est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3),
                     ci_excludes_0, readable) %>%
    as.data.frame()

  # m3, the mandatory proliferation companion, beside it.
  x$contrasts %>%
    dplyr::filter(readout_col == "config_comp", kind == "OXPHOS",
                  model == "m3", subset == "all") %>%
    dplyr::transmute(instrument, proxy, stratum, contrast,
                     est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3)) %>%
    as.data.frame()

  # m2 as declared (1 of 15 readable) and m2L, the Luminal companion.
  x$fits %>%
    dplyr::filter(readout_col == "config_comp", model %in% c("m2", "m2L"),
                  instrument == "gsva") %>%
    dplyr::transmute(proxy, stratum, model, subset, n, min_cell, readable) %>%
    as.data.frame()

  # The arms, then the twelve. Reported, never the reading.
  x$contrasts %>%
    dplyr::filter(kind == "OXPHOS", model == "m1", subset == "all",
                  instrument == "gsva", proxy == "Aneuploidy.Score",
                  readout_col %in% c("arm_trigger", "arm_guardian")) %>%
    dplyr::transmute(readout, stratum, contrast, est = round(est, 3),
                     lo = round(lo, 3), hi = round(hi, 3)) %>%
    as.data.frame()

  # The circularity, as measured.
  x$burden_cor
  round(x$burden_vs_prolif, 3)
  x$asserts %>% as.data.frame()
}
