# E56_quadrant_configuration_fifteen_genes.R
# =============================================================================
# EXTENDED DATA FIG. 6d: E46'S QUADRANT PANEL WITH FIG. 4g'S FIFTEEN GENES,
# GROUPS AND COLOURS. FIRST IT FITS THE ONE NEW QUANTITY, THE THREE EFFECTORS'
# QUADRANT MEANS ON E34'S MODEL, BEHIND THREE GATES.
#
# Specification: the author's ruling of 2026-11-05, verbatim in
# docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md (dd5004a), and the declaration
# of the new quantity, docs/2026-11-06_E56_effector_quadrant_means_declaration.md,
# committed alone at 922f2ad (2026-11-05 23:55 GMT) before this script.
#
# EXPLORATORY and POST HOC. The effectors' means are drawn, not read: no
# reading rule, no label, no test (declaration, section 5).
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNG is compared byte for
# byte with the dry run. It reads the two linear expression snapshots, so it
# takes a minute.
#
# =============================================================================
# WHAT IT DOES
# =============================================================================
#   1. Checks E34's saved object (by the serialize-version-3 md5 its result
#      note records) and E46's (by digest::digest, as its record has it).
#   2. GATE G1: E34's expression path (log2 of the linear snapshot + 1, SCAN-B
#      through symbol_map, Blom scores once per cohort) reproduces the twelve
#      genes in E34's saved frames, to 1e-12.
#   3. GATE G2: E34's model path (.fit_groups, .lincomb, verbatim) reproduces
#      E34's stored m1 group means for the twelve, to 1e-12. GATE G3: it
#      reproduces E46's signature strip, m1 and m3, to 1e-12.
#   4. Fits BAX, BAK1 and BOK the same way: 48 m1 group means.
#   5. Draws Fig. 6d in Fig. 4g's grammar (E50's constants and helpers):
#      - rows, Fig. 4g's fifteen genes in its four groups (the ruling);
#      - columns, both OXPHOS splits at full width, each cohort's four
#        quadrants in E46's order (question 3);
#      - cells, the m1 group means, in Fig. 1g's palette at +/-0.6, colour only,
#        as E46 drew them; the colour key reads "z-score" (the ruling);
#      - below, E46's signature strip, m1 filled and m3 open, keyed "Signature"
#        and "proliferation adjusted" (question 4), its circles 1.05 mm across
#        as in Fig. 4x.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It refits nothing of E34's: the twelve genes and the signature are E34's
#     and E46's stored values, which the gates reproduce.
#   - It reads nothing across repos at runtime; it writes NOTHING to
#     myc_human_validation (d3ac60e) or myc_mouse.
#   - It changes no text of the paper. Its legend draft is Claude Code's.
#   - It leaves no raw file bytes in the workspace: the PDF is read back inside
#     a function that returns only ASCII and numbers (E51's lesson), and the
#     expression matrices are read inside a function that returns only the
#     fifteen genes' rows.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
for (pkg in c("digest", "png", "scales", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E56 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE56: Extended Data Fig. 6d, the quadrants with fifteen genes\n", strrep("=", 78))
message("Fits the effectors on E34's model behind three gates, then draws. ",
        "Specification dd5004a; declaration 922f2ad.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_E34        <- file.path(DIR_RESULTS, "e34_quadrant_configuration.rds")
PATH_E46        <- file.path(DIR_RESULTS, "e46_fig4_quadrant_panels.rds")
PATH_E34_SCRIPT <- here::here("scripts", "E34_quadrant_configuration.R")
PATH_E50_SCRIPT <- here::here("scripts", "E50_fig4g_in_fig1g_style.R")
PATH_RULING     <- here::here("docs", "2026-11-05_E55_E56_ed_fig6cd_ruling.md")
PATH_DECL       <- here::here("docs", "2026-11-06_E56_effector_quadrant_means_declaration.md")
PATH_E56        <- file.path(DIR_RESULTS, "e56_quadrant_configuration_fifteen_genes.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-06"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- "E56_quadrant_configuration_fifteen_genes"

# =============================================================================
# 0. Constants
# =============================================================================
# E34's elements, with the md5 its result note records (docs/2026-10-01_e34_result.md:
# serialize(x, version = 3), then tools::md5sum), which the author's run reproduced.
E34_MD5 <- c(groups          = "86763123ea2d0a725d14ccd4e0e3fc51",
             analysis_frames = "530fd16facc582ef45038cf9f7f2b0b0")
# E46's elements, with the digests its record gives for its final version
# (docs/2026-10-22_fig4_human_panels_proposal.md, the second table).
E46_DIGESTS <- c(heatmap_cells = "326f8533bb96fa2214f141f2f471dad9",
                 strip         = "a65aacc6636ededef218ee137062dde4")

COHORTS <- c("TCGA", "SCAN-B")
INSTRS  <- c("gsva", "mitopps")
HALF_TITLE <- c(gsva = "OXPHOS split by GSVA", mitopps = "OXPHOS split by mitoPPS")  # E46's
EFFECTORS  <- c("BAX", "BAK1", "BOK")

# --- E34's definitions, verbatim (checked against its script below) -------------
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
TWELVE <- c(PRIMING_PRO, PRIMING_ANTI)
QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",
          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")
COV <- list(
  TCGA     = list(m1 = c("PAM50", "purity", "leuko"),
                  m2 = c("purity", "leuko"),
                  m3 = c("PAM50", "purity", "leuko", "PROLIF")),
  `SCAN-B` = list(m1 = "PAM50",
                  m2 = character(0),
                  m3 = c("PAM50", "PROLIF")))
# E46's columns: within each cohort, OXPHOS-low then OXPHOS-high, MYC-low then
# MYC-high within each.
COL_Q <- c(QLEV[1], QLEV[3], QLEV[2], QLEV[4])

# --- Fig. 4g's rows, palette and geometry: E50's lines, verbatim -------------------
GROUPS <- list(Effectors = c("BAX", "BAK1", "BOK"),
               Activators = c("BID", "BBC3", "BCL2L11"),
               Sensitisers = c("BAD", "PMAIP1", "BMF", "BIK"),
               `Anti-apoptotic` = c("MCL1", "BCL2L2", "BCL2A1", "BCL2",
                                    "BCL2L1"))
FIFTEEN   <- unlist(GROUPS, use.names = FALSE)
DIV_COLS <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")
MM       <- 72 / 25.4
CELL_H   <- 8.46
GAP      <- 0.8 * MM     # Fig. 1g: panel.spacing.y = 0.8 mm (2.27 pt)
FONT     <- "Helvetica"
PT_TXT   <- 5            # values and every label (item 3)
PT_HEAD  <- 6            # the score headers (item 3)
CAP      <- 0.718        # Helvetica's cap height, em (AFM)
DESC     <- 0.207        # Helvetica's descender, em (AFM)
LAB_GAP  <- 1.2          # gene symbols end 1.2 pt left of the cells (Fig. 1g)
COH_GAP  <- 1.2          # cohort labels' baseline 1.2 pt above the cells
GRP_GAP  <- 2.4          # group names start 2.4 pt right of the cells
MARGIN   <- c(top = 1, right = 2.5, bottom = 1, left = 1) * MM  # plot.margin
BAR_W    <- 2 * MM       # Fig. 1g's colour bar
BAR_H    <- 18 * MM
KEY_GAP  <- 1 * MM       # 1 mm past the wider of the first two group names
TITLE_GAP <- 2.4         # bar title's descender to the bar (Fig. 1g)
TICK_LEN <- 0.6 * MM     # 1.70 pt, Fig. 1g's
TICK_LWD <- 0.375        # pt
TICK_COL <- "#4D4D4D"    # grey30
TICKLAB_GAP <- 3.0       # bar to tick labels (Fig. 1g: 2.99 pt)
N_STOPS  <- 41L          # gradient stops on the bar, each the palette's colour
PNG_DPI  <- 600
MINUS_U  <- intToUtf8(8722L) # the minus sign, in the PNG (U+2212)
SOFT_HY  <- intToUtf8(173L)  # the PDF's hyphen; its encoding prints "-" as minus (U+00AD)
HJ <- c(start = 0, middle = 0.5, end = 1)

# --- Fig. 6d's own layout, in points ----------------------------------------------------
L_SCALE   <- 0.6                 # Fig. 4g's scale (the ruling: the same colour code)
BAR_TITLE <- "z-score"           # the ruling
KEY_TEXT  <- c(m1 = "Signature", m3 = "proliferation adjusted")   # the ruling, question 4
STRIP_TITLE <- "Signature (SD)"  # E46's y axis title for the strip
W_PAGE_PT <- 510                 # 179.9 mm: Extended Data's full width, in whole points
CGAP      <- 3 * GAP             # between the two cohorts of a split
HGAP      <- 8 * MM              # between the two splits
HEAD_STEP <- 2.4 * MM            # header lines, baseline to baseline
TITLE_STEP <- 3.0 * MM           # the split titles to the cohort headers, a little more
S_GAP     <- 3 * MM              # the heatmap to the strip
S_H       <- 24 * MM             # the strip's height
S_TICKS   <- c(-0.6, -0.3, 0, 0.3, 0.6)   # E46's
AX_GAP    <- 1.5                 # pt: a strip's y axis left of its first column
AX_TICK   <- 1 * MM
AX_LWD_PT <- 0.53                # Fig. 1's axis line (NM_LINE)
CI_PT     <- 0.6                 # pt: the strip's intervals
NULL_COL  <- "#737373"; NULL_PT <- 0.5    # grey45, the line at zero
CIRCLE_MM <- 1.05                # as Fig. 4x (the author's point-size ruling)
CIRC_LW_PT <- 0.425              # its outline: E54's 0.3 stroke
DOT_DX    <- 0.2                 # m1 left and m3 right of a column's centre, in columns
LINE5     <- 5 * 1.15            # pt: two lines of 5-pt text, baseline to baseline

.stop_if <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
ASSERT <- list()
.assert <- function(check, ok, n = NA_integer_, max_dev = NA_real_,
                    detail = "") {
  ASSERT[[length(ASSERT) + 1L]] <<- tibble::tibble(
    check = check, n_compared = as.integer(n), max_abs_dev = max_dev,
    pass = isTRUE(ok), detail = detail)
  .stop_if(ok, "ASSERTION FAILED: ", check,
           if (nzchar(detail)) paste0(" (", detail, ")") else "", ". STOP.")
  invisible(TRUE)
}

.assert("the panel is the ruling's, and the new fit is declared in its own committed file",
        file.exists(PATH_RULING) && file.exists(PATH_DECL) &&
          any(grepl("Fit the effectors", readLines(PATH_RULING, warn = FALSE), fixed = TRUE)) &&
          any(grepl("Both splits, full width", readLines(PATH_RULING, warn = FALSE), fixed = TRUE)) &&
          any(grepl("48 estimates", readLines(PATH_DECL, warn = FALSE), fixed = TRUE)), n = 2L)
E34_SRC <- readLines(PATH_E34_SCRIPT, warn = FALSE)
E50_SRC <- readLines(PATH_E50_SCRIPT, warn = FALSE)
E34_LINES <- c(
  'PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")',
  'PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")',
  'QLEV <- c("Q1 MYC-low OXPHOS-low", "Q2 MYC-low OXPHOS-high",',
  '          "Q3 MYC-high OXPHOS-low", "Q4 MYC-high OXPHOS-high")',
  '  TCGA     = list(m1 = c("PAM50", "purity", "leuko"),',
  '                  m3 = c("PAM50", "purity", "leuko", "PROLIF")),',
  '  `SCAN-B` = list(m1 = "PAM50",', '                  m3 = c("PAM50", "PROLIF")))')
E50_LINES <- c(
  'GROUPS <- list(Effectors = c("BAX", "BAK1", "BOK"),',
  '               Activators = c("BID", "BBC3", "BCL2L11"),',
  '               Sensitisers = c("BAD", "PMAIP1", "BMF", "BIK"),',
  '               `Anti-apoptotic` = c("MCL1", "BCL2L2", "BCL2A1", "BCL2",',
  'DIV_COLS <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")',
  'CELL_H   <- 8.46',
  'GAP      <- 0.8 * MM     # Fig. 1g: panel.spacing.y = 0.8 mm (2.27 pt)',
  'LAB_GAP  <- 1.2          # gene symbols end 1.2 pt left of the cells (Fig. 1g)',
  'GRP_GAP  <- 2.4          # group names start 2.4 pt right of the cells',
  'BAR_W    <- 2 * MM       # Fig. 1g\'s colour bar', 'BAR_H    <- 18 * MM',
  'KEY_GAP  <- 1 * MM       # 1 mm past the wider of the first two group names',
  'TICK_LEN <- 0.6 * MM     # 1.70 pt, Fig. 1g\'s',
  'TICKLAB_GAP <- 3.0       # bar to tick labels (Fig. 1g: 2.99 pt)')
.assert("E34's definitions and Fig. 4g's rows, palette and geometry are their scripts', line for line",
        all(E34_LINES %in% E34_SRC) && all(E50_LINES %in% E50_SRC),
        n = length(E34_LINES) + length(E50_LINES),
        detail = paste(utils::head(c(setdiff(E34_LINES, E34_SRC), setdiff(E50_LINES, E50_SRC)), 2),
                       collapse = " / "))

# =============================================================================
# 1. E34's and E46's objects, checked against their records
# =============================================================================
message("\n1. E34's and E46's saved objects")
.stop_if(file.exists(PATH_E34) && file.exists(PATH_E46), "E34's or E46's saved object is missing. STOP.")
e34 <- readRDS(PATH_E34)
e46 <- readRDS(PATH_E46)
.md5_v3 <- function(o) {
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb")
  serialize(o, con, version = 3)
  close(con)
  unname(tools::md5sum(f))
}
GOT34 <- vapply(names(E34_MD5), function(k) .md5_v3(e34[[k]]), character(1))
GOT46 <- vapply(names(E46_DIGESTS), function(k) digest::digest(e46[[k]]), character(1))
.assert(paste0("E34's and E46's objects are the verified ones, each by the checksum ",
               "its own record uses"),
        identical(unname(GOT34), unname(E34_MD5)) && identical(unname(GOT46), unname(E46_DIGESTS)),
        n = 4L)
FR34 <- e34$analysis_frames
G34  <- as.data.frame(e34$groups)

# =============================================================================
# 2. Gate G1: the expression path reproduces E34's twelve genes
# =============================================================================
message("\n2. Gate G1: the expression path")
sph <- readRDS(PATH_SCANB_PHENO)
# The linear matrix is read inside a function that returns only the requested
# genes' rows, on log2(linear + 1), for the tumours given, in their order.
.log2_rows <- function(path, ids, map, genes) {
  lin <- readRDS(path)
  .stop_if(identical(lin$scale, "linear_deseq2_normalised"), basename(path),
           " is not on the linear DESeq2-normalised scale. STOP.")
  .stop_if(all(ids %in% colnames(lin$mat)), basename(path), ": a tumour is missing. STOP.")
  G  <- log2(lin$mat[, ids, drop = FALSE] + 1)
  gr <- .gene_rows(genes, G, .symbol_resolver(rownames(G), map))
  .stop_if(!length(gr$missing), basename(path), ": did not resolve ",
           paste(gr$missing, collapse = ", "), ". STOP.")
  gr$mat[genes, , drop = FALSE]
}
GENES_ALL <- c(TWELVE, EFFECTORS)
LOG2 <- list(
  TCGA     = .log2_rows(PATH_TCGA_LINEAR, FR34$TCGA$sample_id, NULL, GENES_ALL),
  `SCAN-B` = .log2_rows(PATH_SCANB_LINEAR, FR34$`SCAN-B`$sample_id, sph$symbol_map, GENES_ALL))
invisible(gc(verbose = FALSE))
# E34's scorer, verbatim.
.blom <- function(x) {
  n <- length(x)
  stopifnot(!anyNA(x))
  z <- stats::qnorm((rank(x, ties.method = "average") - 0.375) / (n + 0.25))
  as.numeric(scale(z))
}
dev_g1 <- max(vapply(COHORTS, function(coh) max(vapply(TWELVE, function(g)
  max(abs(.blom(as.numeric(LOG2[[coh]][g, ])) - FR34[[coh]][[paste0("g_", g)]])),
  numeric(1))), numeric(1)))
.assert(paste0("G1: log2(linear + 1), Blom-scored once per cohort, reproduces the ",
               "twelve genes in E34's saved frames, in both cohorts"),
        dev_g1 < 1e-12 && nrow(FR34$TCGA) == 1095L && nrow(FR34$`SCAN-B`) == 3207L,
        n = 2L * length(TWELVE), max_dev = dev_g1)
FRX <- stats::setNames(lapply(COHORTS, function(coh) {
  f <- FR34[[coh]]
  for (e in EFFECTORS) {
    f[[paste0("raw_g_", e)]] <- as.numeric(LOG2[[coh]][e, ])
    f[[paste0("g_", e)]] <- .blom(as.numeric(LOG2[[coh]][e, ]))
  }
  f
}), COHORTS)
rm(LOG2)

# =============================================================================
# 3. Gates G2 and G3: the model path reproduces E34's and E46's estimates
# =============================================================================
message("\n3. Gates G2 and G3: the model path")
# E34's group means, verbatim.
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
.groups_of <- function(coh, inst, y, model) {
  d0 <- FRX[[coh]]
  d0$quad <- d0[[paste0("quad_", inst)]]
  fit <- .fit_groups(d0, y, COV[[coh]][[model]])
  .stop_if(fit$ok, "a fit failed: ", coh, " ", inst, " ", y, " ", model, ". STOP.")
  G <- t(vapply(QLEV, function(q) .lincomb(fit, fit$A[, q]), numeric(4)))
  list(fit = fit, tab = tibble::tibble(cohort = coh, instrument = inst, readout_col = y,
                                       model = model, quadrant = QLEV,
                                       n_cell = unname(fit$cell_n[QLEV]),
                                       est = G[, "est"], se = G[, "se"], lo = G[, "lo"],
                                       hi = G[, "hi"]))
}
KEYS <- tidyr::expand_grid(cohort = COHORTS, instrument = INSTRS)
G2_NEW <- dplyr::bind_rows(lapply(seq_len(nrow(KEYS)), function(i)
  dplyr::bind_rows(lapply(paste0("g_", TWELVE), function(y)
    .groups_of(KEYS$cohort[i], KEYS$instrument[i], y, "m1")$tab))))
G2_OLD <- G34 %>% dplyr::filter(model == "m1", subtype == "all", readout_col %in% paste0("g_", TWELVE)) %>%
  dplyr::select(cohort, instrument, readout_col, quadrant, n_cell, est, se, lo, hi)
G2J <- dplyr::inner_join(G2_NEW, G2_OLD, by = c("cohort", "instrument", "readout_col", "quadrant"),
                         suffix = c("", "_e34"))
dev_g2 <- max(abs(c(G2J$est - G2J$est_e34, G2J$se - G2J$se_e34, G2J$lo - G2J$lo_e34,
                    G2J$hi - G2J$hi_e34)))
.assert(paste0("G2: E34's model path reproduces its stored m1 group means for the ",
               "twelve genes (estimate, SE, both limits), both splits, both cohorts"),
        nrow(G2J) == 192L && nrow(G2_OLD) == 192L && dev_g2 < 1e-12 &&
          identical(G2J$n_cell, G2J$n_cell_e34), n = 4L * nrow(G2J), max_dev = dev_g2)
G3_NEW <- dplyr::bind_rows(lapply(seq_len(nrow(KEYS)), function(i)
  dplyr::bind_rows(lapply(c("m1", "m3"), function(m)
    .groups_of(KEYS$cohort[i], KEYS$instrument[i], "config_comp", m)$tab))))
STRIP <- as.data.frame(e46$strip)
G3J <- dplyr::inner_join(G3_NEW, STRIP %>% dplyr::select(cohort, instrument, model, quadrant,
                                                          est46 = est, se46 = se, lo46 = lo, hi46 = hi),
                         by = c("cohort", "instrument", "model", "quadrant"))
dev_g3 <- max(abs(c(G3J$est - G3J$est46, G3J$se - G3J$se46, G3J$lo - G3J$lo46, G3J$hi - G3J$hi46)))
.assert(paste0("G3: the same path reproduces E46's signature strip, m1 and m3, both ",
               "splits, both cohorts"),
        nrow(G3J) == 32L && nrow(STRIP) == 32L && dev_g3 < 1e-12, n = 4L * nrow(G3J),
        max_dev = dev_g3)
HM46 <- as.data.frame(e46$heatmap_cells)
H46J <- dplyr::inner_join(HM46 %>% dplyr::mutate(readout_col = paste0("g_", gene)),
                          G2_OLD %>% dplyr::select(cohort, instrument, readout_col, quadrant, est34 = est),
                          by = c("cohort", "instrument", "readout_col", "quadrant"))
.assert("E46's heatmap cells are E34's stored m1 means for the twelve genes",
        nrow(H46J) == 192L && max(abs(H46J$est - H46J$est34)) < 1e-12, n = 192L,
        max_dev = max(abs(H46J$est - H46J$est34)))

# =============================================================================
# 4. The new quantity: the effectors' m1 group means (declaration, section 2)
# =============================================================================
message("\n4. The effectors, fitted on E34's m1")
EFF_FITS <- list()
EFF <- dplyr::bind_rows(lapply(seq_len(nrow(KEYS)), function(i)
  dplyr::bind_rows(lapply(EFFECTORS, function(e) {
    r <- .groups_of(KEYS$cohort[i], KEYS$instrument[i], paste0("g_", e), "m1")
    EFF_FITS[[paste(KEYS$cohort[i], KEYS$instrument[i], e)]] <<- r$fit
    r$tab %>% dplyr::mutate(gene = e)
  }))))
.assert("48 effector means, each finite, each interval around its estimate",
        nrow(EFF) == 48L && all(is.finite(c(EFF$est, EFF$lo, EFF$hi))) &&
          all(EFF$lo < EFF$est & EFF$est < EFF$hi), n = 48L)
N_COH <- vapply(COHORTS, function(coh) sum(G2_OLD$n_cell[G2_OLD$cohort == coh &
                                                           G2_OLD$instrument == "gsva" &
                                                           G2_OLD$readout_col == "g_BAD"]),
                numeric(1))
.assert(paste0("each effector fit uses the same tumours as E34's fits: n = 938 in ",
               "TCGA-BRCA and 3,143 in SCAN-B, with the same quadrant counts"),
        all(vapply(EFF_FITS, function(f) f$n, numeric(1)) ==
              N_COH[sub(" .*$", "", names(EFF_FITS))]) &&
          identical(unname(N_COH), c(938, 3143)), n = length(EFF_FITS),
        detail = paste(N_COH, collapse = ", "))
# R rules 4 and 5: the full coefficient vectors, their variance matrices and a
# digest table; the lm objects are refitted alike and saved, not digested.
.lm_of <- function(coh, inst, y) {
  d0 <- FRX[[coh]]
  d0$quad <- d0[[paste0("quad_", inst)]]
  covars <- COV[[coh]][["m1"]]
  cols <- c(y, "quad", covars)
  dd <- as.data.frame(d0[stats::complete.cases(d0[, cols]), cols, drop = FALSE])
  names(dd)[1] <- "yv"
  dd$quad <- droplevels(factor(dd$quad, levels = QLEV))
  for (cv in covars) if (is.factor(dd[[cv]])) dd[[cv]] <- droplevels(dd[[cv]])
  stats::lm(stats::reformulate(c("quad", covars), response = "yv"), data = dd)
}
EFF_LM <- list()
for (k in names(EFF_FITS)) {
  p <- strsplit(k, " ")[[1]]
  EFF_LM[[k]] <- .lm_of(paste(p[seq_len(length(p) - 2L)], collapse = " "),
                        p[length(p) - 1L], paste0("g_", p[length(p)]))
}
.assert("the saved lm objects are the fits the means come from: identical coefficients",
        all(vapply(names(EFF_FITS), function(k)
          identical(stats::coef(EFF_LM[[k]]), EFF_FITS[[k]]$b), logical(1))),
        n = length(EFF_FITS))
FIT_DIGESTS <- tibble::tibble(
  fit = names(EFF_FITS),
  coef_digest = vapply(EFF_FITS, function(f) digest::digest(f$b), character(1)),
  vcov_digest = vapply(EFF_FITS, function(f) digest::digest(f$V), character(1)),
  n = vapply(EFF_FITS, function(f) f$n, numeric(1)),
  df = vapply(EFF_FITS, function(f) f$df, numeric(1)))
for (e in EFFECTORS) {
  r <- EFF %>% dplyr::filter(gene == e)
  message(sprintf("   %-4s means from %+.2f to %+.2f across the 16 cells", e, min(r$est), max(r$est)))
}

# =============================================================================
# 5. The values drawn, and Fig. 1g's palette (E50's helpers, verbatim)
# =============================================================================
message("\n5. values and palette")
VAL <- dplyr::bind_rows(
  HM46 %>% dplyr::transmute(gene, cohort, instrument, quadrant, est, source = "E34, via E46"),
  EFF %>% dplyr::transmute(gene, cohort, instrument, quadrant, est, source = "E56, declared"))
.assert("240 cells: fifteen genes, two splits, two cohorts, four quadrants, each once",
        nrow(VAL) == 240L && setequal(VAL$gene, FIFTEEN) &&
          nrow(dplyr::distinct(VAL, gene, cohort, instrument, quadrant)) == 240L, n = 240L)
.fig1g_scale <- function(L, ...) {
  ggplot2::scale_fill_gradientn(
    colours = unname(DIV_COLS[c("neg", "zero", "pos")]),
    values = scales::rescale(c(-L, 0, L)), limits = c(-L, L),
    oob = scales::squish, ...)
}
.fill_hex <- function(v, L) {
  pal <- scales::pal_gradient_n(unname(DIV_COLS[c("neg", "zero", "pos")]),
                                values = scales::rescale(c(-L, 0, L)))
  toupper(pal(scales::rescale(scales::oob_squish(v, c(-L, L)), from = c(-L, L))))
}
.ggplot_fill <- function(v, L) {
  d <- data.frame(x = seq_along(v), y = 1, v = v)
  b <- ggplot2::ggplot_build(ggplot2::ggplot(d, ggplot2::aes(x, y, fill = v)) +
                               ggplot2::geom_tile() + .fig1g_scale(L))
  dd <- b$data[[1]]
  toupper(dd$fill[order(dd$x)])
}
lab_num <- function(x) {
  out <- format(x, big.mark = ",", scientific = FALSE, trim = TRUE,
                drop0trailing = TRUE)
  out[is.na(x)] <- NA_character_
  out
}
lab_signed <- function(x) {
  out <- lab_num(x)
  pos <- !is.na(x) & x > 0
  out[pos] <- paste0("+", out[pos])
  out
}
.to_pdf  <- function(u) gsub(MINUS_U, "-", gsub("-", SOFT_HY, u, fixed = TRUE),
                             fixed = TRUE)
.text_w <- function(s, size, face) {
  n <- length(s)
  size <- rep_len(size, n)
  face <- rep_len(face, n)
  prev <- grDevices::dev.cur()
  grDevices::pdf(NULL, family = FONT)
  cur <- grDevices::dev.cur()
  on.exit({
    grDevices::dev.off(cur)
    if (prev > 1L) grDevices::dev.set(prev)
  }, add = TRUE)
  vapply(seq_len(n), function(i) grid::convertWidth(grid::grobWidth(
    grid::textGrob(s[i], gp = grid::gpar(fontfamily = FONT, fontsize = size[i],
                                         fontface = face[i]))),
    "bigpts", valueOnly = TRUE), numeric(1))
}
.txt <- function(role, x, y, anchor, size, face, label_uni,
                 ink = "#000000", key = NA_character_) {
  tibble::tibble(role = role, x = x, y = y, anchor = anchor, size = size,
                 face = as.integer(face), label_uni = label_uni, ink = ink,
                 key = key)
}
BEYOND <- VAL %>% dplyr::filter(abs(est) > L_SCALE)
.assert(paste0("every cell's colour is the one Fig. 1g's own scale call gives its ",
               "value at +/-0.6"),
        identical(.fill_hex(c(-L_SCALE, 0, L_SCALE), L_SCALE), unname(DIV_COLS)) &&
          identical(.fill_hex(VAL$est, L_SCALE), .ggplot_fill(VAL$est, L_SCALE)),
        n = nrow(VAL), detail = sprintf("%d beyond +/-0.6", nrow(BEYOND)))

# =============================================================================
# 6. The layout, in points from the top left
# =============================================================================
message("\n6. layout")
.layout_6d <- function(VAL, STRIP) {
  L <- L_SCALE
  rows <- tibble::tibble(gene = FIFTEEN, group = rep(names(GROUPS), lengths(GROUPS))) %>%
    dplyr::mutate(k = dplyr::row_number(), g = match(group, names(GROUPS)))
  lab_w <- max(.text_w(rows$gene, PT_TXT, 3L))
  ytick_lab <- sub("^-", MINUS_U, lab_signed(S_TICKS))
  ytick_w <- max(.text_w(.to_pdf(ytick_lab), PT_TXT, 1L))
  # The strip's y axis, its ticks, their labels and the rotated title share the
  # left margin with the gene symbols.
  y_zone <- AX_GAP + AX_TICK + LAB_GAP + ytick_w + 2 + (CAP + DESC) * PT_HEAD
  LEFT <- max(lab_w + LAB_GAP, y_zone)
  X0 <- MARGIN[["left"]] + LEFT
  grp_w12 <- max(.text_w(c("Effectors", "Activators"), PT_TXT, 1L))
  bar_tick_lab <- sub("^-", MINUS_U, lab_signed(c(-L, -L / 2, 0, L / 2, L)))
  bar_tick_w <- max(.text_w(.to_pdf(bar_tick_lab), PT_TXT, 1L))
  RIGHT <- GRP_GAP + grp_w12 + KEY_GAP + BAR_W + TICKLAB_GAP + bar_tick_w
  CW <- (W_PAGE_PT - MARGIN[["left"]] - LEFT - RIGHT - MARGIN[["right"]] -
           4 * GAP - 2 * CGAP - HGAP) / 16
  HALF_W <- 8 * CW + 2 * GAP + CGAP
  cols <- tidyr::expand_grid(instrument = INSTRS, cohort = COHORTS, quadrant = COL_Q) %>%
    dplyr::mutate(h = match(instrument, INSTRS), c = match(cohort, COHORTS),
                  j = match(quadrant, COL_Q),
                  x = X0 + (h - 1) * (HALF_W + HGAP) + (c - 1) * (4 * CW + GAP + CGAP) +
                    (j - 1) * CW + (j >= 3) * GAP)
  X1 <- max(cols$x) + CW
  # Down: the split titles, the cohort headers, the OXPHOS headers, the cells.
  Y_T1 <- MARGIN[["top"]] + CAP * PT_HEAD
  Y_T2 <- Y_T1 + TITLE_STEP
  Y_T3 <- Y_T2 + HEAD_STEP
  Y0 <- Y_T3 + DESC * PT_TXT + COH_GAP
  rows <- rows %>% dplyr::mutate(y = Y0 + (k - 1) * CELL_H + (g - 1) * GAP)
  YH1 <- max(rows$y) + CELL_H
  rects <- VAL %>%
    dplyr::inner_join(rows %>% dplyr::select(gene, y), by = "gene") %>%
    dplyr::inner_join(cols %>% dplyr::select(instrument, cohort, quadrant, x),
                      by = c("instrument", "cohort", "quadrant")) %>%
    dplyr::arrange(match(gene, FIFTEEN), match(instrument, INSTRS), match(cohort, COHORTS),
                   match(quadrant, COL_Q)) %>%
    dplyr::mutate(w = CW, h = CELL_H, fill = .fill_hex(est, L))
  # The strip: one band under the cells, its own y axis left of each split.
  YLIM <- ceiling(max(abs(c(STRIP$lo, STRIP$hi))) * 10 - 1e-9) / 10
  Y_S0 <- YH1 + S_GAP
  Y_S1 <- Y_S0 + S_H
  .ys <- function(v) Y_S0 + (YLIM - v) / (2 * YLIM) * S_H
  dots <- STRIP %>% dplyr::select(-x) %>%           # E46's own x positions are not used
    dplyr::inner_join(cols %>% dplyr::select(instrument, cohort, quadrant, x),
                      by = c("instrument", "cohort", "quadrant")) %>%
    dplyr::mutate(xd = x + CW / 2 + ifelse(model == "m1", -1, 1) * DOT_DX * CW,
                  yd = .ys(est), y_lo = .ys(lo), y_hi = .ys(hi),
                  filled = model == "m1")
  halves <- cols %>% dplyr::group_by(instrument) %>%
    dplyr::summarise(x0 = min(x), x1 = max(x) + CW, .groups = "drop") %>%
    dplyr::arrange(match(instrument, INSTRS))
  axes <- halves %>% dplyr::mutate(xa = x0 - AX_GAP, ya0 = .ys(YLIM), ya1 = .ys(-YLIM))
  yticks <- tidyr::expand_grid(instrument = INSTRS, v = S_TICKS) %>%
    dplyr::left_join(axes %>% dplyr::select(instrument, xa), by = "instrument") %>%
    dplyr::mutate(x0 = xa - AX_TICK, x1 = xa, y = .ys(v))
  zeros <- halves %>% dplyr::mutate(y = .ys(0))
  # The MYC labels under each column, on two lines.
  Y_M1 <- Y_S1 + 1.5 * MM + CAP * PT_TXT
  Y_M2 <- Y_M1 + LINE5
  # The strip's key, under the first split, from the cells' left edge.
  Y_K <- Y_M2 + DESC * PT_TXT + 2.5 * MM + CAP * PT_TXT
  CIRC_PT <- CIRCLE_MM * MM
  key_x1 <- X0 + CIRC_PT / 2
  key_t1 <- X0 + CIRC_PT + 0.8 * MM
  key_x2 <- key_t1 + .text_w(KEY_TEXT[["m1"]], PT_TXT, 1L) + 3 * MM + CIRC_PT / 2
  key_t2 <- key_x2 + CIRC_PT / 2 + 0.8 * MM
  keys <- tibble::tibble(model = c("m1", "m3"), xd = c(key_x1, key_x2),
                         yd = Y_K - CAP * PT_TXT / 2, filled = c(TRUE, FALSE))
  # The colour bar, placed as E50 places 4g's.
  XG <- X1 + GRP_GAP
  XB <- XG + grp_w12 + KEY_GAP
  Y_TITLE <- Y0 + CAP * PT_TXT
  YB_TOP  <- Y_TITLE + DESC * PT_TXT + TITLE_GAP
  YB_BOT  <- YB_TOP + BAR_H
  TICKS   <- c(-L, -L / 2, 0, L / 2, L)
  y_tick  <- YB_BOT - (TICKS + L) / (2 * L) * BAR_H
  grp <- rows %>% dplyr::group_by(group) %>%
    dplyr::summarise(mid = (min(y) + max(y) + CELL_H) / 2, .groups = "drop") %>%
    dplyr::arrange(match(group, names(GROUPS)))
  coh_mid <- cols %>% dplyr::group_by(instrument, cohort) %>%
    dplyr::summarise(x = (min(x) + max(x) + CW) / 2, .groups = "drop")
  ox_mid <- cols %>% dplyr::mutate(ox = ifelse(j <= 2, "OXPHOS low", "OXPHOS high")) %>%
    dplyr::group_by(instrument, cohort, ox) %>%
    dplyr::summarise(x = (min(x) + max(x) + CW) / 2, .groups = "drop")
  n_lab <- c(TCGA = paste0("TCGA (n = ", format(N_COH[["TCGA"]], big.mark = ","), ")"),
             `SCAN-B` = paste0("SCAN-B (n = ", format(N_COH[["SCAN-B"]], big.mark = ","), ")"))
  myc <- cols %>% dplyr::mutate(xc = x + CW / 2,
                                lvl = ifelse(grepl("MYC-low", quadrant), "low", "high"))
  y_mid <- (Y_S0 + Y_S1) / 2
  x_ytitle <- axes$xa[1] - AX_TICK - LAB_GAP - ytick_w - 2 - DESC * PT_HEAD
  texts <- dplyr::bind_rows(
    .txt("gene", X0 - LAB_GAP, rows$y + CELL_H / 2 + CAP * PT_TXT / 2, "end", PT_TXT, 3L,
         rows$gene, key = rows$gene),
    .txt("group", XG, grp$mid + CAP * PT_TXT / 2, "start", PT_TXT, 1L, grp$group,
         key = grp$group),
    .txt("split", (halves$x0 + halves$x1) / 2, Y_T1, "middle", PT_HEAD, 1L,
         unname(HALF_TITLE[halves$instrument]), key = halves$instrument),
    .txt("cohort", coh_mid$x, Y_T2, "middle", PT_HEAD, 1L, unname(n_lab[coh_mid$cohort]),
         key = paste(coh_mid$instrument, coh_mid$cohort)),
    .txt("oxphos", ox_mid$x, Y_T3, "middle", PT_TXT, 1L, ox_mid$ox,
         key = paste(ox_mid$instrument, ox_mid$cohort, ox_mid$ox)),
    .txt("myc", myc$xc, Y_M1, "middle", PT_TXT, 1L, "MYC",
         key = paste(myc$instrument, myc$cohort, myc$quadrant, "1")),
    .txt("myc", myc$xc, Y_M2, "middle", PT_TXT, 1L, myc$lvl,
         key = paste(myc$instrument, myc$cohort, myc$quadrant, "2")),
    .txt("ytick", axes$xa[1] - AX_TICK - LAB_GAP, .ys(S_TICKS) + CAP * PT_TXT / 2, "end",
         PT_TXT, 1L, ytick_lab, key = sprintf("%+.1f", S_TICKS)),
    .txt("ytitle", x_ytitle, y_mid, "middle", PT_HEAD, 1L, STRIP_TITLE),
    .txt("key", c(key_t1, key_t2), Y_K, "start", PT_TXT, 1L, unname(KEY_TEXT),
         key = names(KEY_TEXT)),
    .txt("bar_title", XB, Y_TITLE, "start", PT_TXT, 1L, BAR_TITLE),
    .txt("tick", XB + BAR_W + TICKLAB_GAP, y_tick + CAP * PT_TXT / 2, "start", PT_TXT, 1L,
         bar_tick_lab, key = sprintf("%+.2f", TICKS))) %>%
    dplyr::mutate(rot = ifelse(role == "ytitle", 90, 0),
                  label_pdf = .to_pdf(label_uni),
                  width = .text_w(label_pdf, size, face),
                  left = ifelse(rot == 90, x - CAP * size, x - unname(HJ[anchor]) * width),
                  right = ifelse(rot == 90, x + DESC * size, x - unname(HJ[anchor]) * width + width),
                  top = ifelse(rot == 90, y - width / 2, y - CAP * size),
                  bottom = ifelse(rot == 90, y + width / 2, y + DESC * size))
  bar <- tibble::tibble(x = XB, y = YB_TOP, w = BAR_W, h = BAR_H)
  stops <- tibble::tibble(offset = seq(0, 1, length.out = N_STOPS)) %>%
    dplyr::mutate(value = -L + offset * 2 * L, colour = .fill_hex(value, L))
  bticks <- dplyr::bind_rows(
    tibble::tibble(edge = "left", value = TICKS, x0 = XB, x1 = XB + TICK_LEN,
                   y0 = y_tick, y1 = y_tick),
    tibble::tibble(edge = "right", value = TICKS, x0 = XB + BAR_W - TICK_LEN,
                   x1 = XB + BAR_W, y0 = y_tick, y1 = y_tick))
  W <- max(texts$right, bar$x + bar$w) + MARGIN[["right"]]
  H <- max(rects$y + rects$h, YB_BOT, texts$bottom, Y_K + DESC * PT_TXT) + MARGIN[["bottom"]]
  list(L = L, CW = CW, YLIM = YLIM, rows = rows, cols = cols, rects = rects, dots = dots,
       keys = keys, axes = axes, yticks = yticks, zeros = zeros, texts = texts, bar = bar,
       stops = stops, bticks = bticks, size = c(W = W, H = H), x0 = X0, x1 = X1, y0 = Y0,
       strip = c(top = Y_S0, bottom = Y_S1))
}
LAY <- .layout_6d(VAL, STRIP)
r <- LAY$rects
tx <- LAY$texts
.assert(paste0("240 cells, Fig. 4g's rows (8.46 pt, 0.8 mm between groups) in its ",
               "four groups and order; columns in E46's order"),
        nrow(r) == 240L && all(r$h == CELL_H) && all(abs(r$w - LAY$CW) < 1e-9) &&
          identical(LAY$rows$gene[order(LAY$rows$y)], FIFTEEN) &&
          identical(LAY$cols$quadrant[LAY$cols$instrument == "gsva" & LAY$cols$cohort == "TCGA"],
                    COL_Q), n = 240L,
        detail = sprintf("columns %.2f mm", LAY$CW / MM))
.assert("the page fits Extended Data's full width: at most 180 mm",
        LAY$size[["W"]] / MM <= 180, n = 1L,
        detail = sprintf("%.1f x %.1f mm", LAY$size[["W"]] / MM, LAY$size[["H"]] / MM))
ov <- outer(seq_len(nrow(tx)), seq_len(nrow(tx)), function(i, j)
  tx$left[i] < tx$right[j] & tx$left[j] < tx$right[i] &
    tx$top[i] < tx$bottom[j] & tx$top[j] < tx$bottom[i])
n_ov <- sum(ov[upper.tri(ov)])
bb <- LAY$bar
hit_bar <- tx$left < bb$x + bb$w & bb$x < tx$right & tx$top < bb$y + bb$h & bb$y < tx$bottom
cc <- outer(seq_len(nrow(tx)), seq_len(nrow(r)), function(i, k)
  tx$left[i] < r$x[k] + r$w[k] & r$x[k] < tx$right[i] &
    tx$top[i] < r$y[k] + r$h[k] & r$y[k] < tx$bottom[i])
.assert(paste0("no two labels overlap; no label touches a cell or the colour bar; ",
               "everything lies inside the page"),
        n_ov == 0L && !any(hit_bar) && sum(cc) == 0L && all(tx$left >= 0) &&
          all(tx$right <= LAY$size[["W"]]) && all(tx$top >= 0) &&
          all(tx$bottom <= LAY$size[["H"]]), n = nrow(tx),
        detail = sprintf("%d overlaps, %d label-cell contacts", n_ov, sum(cc)))
dd <- LAY$dots
.assert(paste0("the strip: E46's 32 signature means, each inside its axis, m1 left ",
               "of m3 in its column"),
        nrow(dd) == 32L && all(dd$y_hi >= LAY$strip[["top"]] & dd$y_lo <= LAY$strip[["bottom"]]) &&
          all(tapply(dd$xd[dd$model == "m1"], paste(dd$instrument, dd$cohort, dd$quadrant)[dd$model == "m1"], min) <
                tapply(dd$xd[dd$model == "m3"], paste(dd$instrument, dd$cohort, dd$quadrant)[dd$model == "m3"], min)),
        n = 32L, detail = sprintf("axis +/-%.1f", LAY$YLIM))
message(sprintf("   page %.1f x %.1f mm; columns %.2f mm; strip axis +/-%.1f",
                LAY$size[["W"]] / MM, LAY$size[["H"]] / MM, LAY$CW / MM, LAY$YLIM))

# =============================================================================
# 7. The files: the PNG, then the PDF, each read back
# =============================================================================
message("\n7. The files")
.draw <- function(lay, variant) {
  H <- lay$size[["H"]]
  bp <- function(v) grid::unit(v, "bigpts")
  lwd <- function(pt) pt * 96 / 72
  grid::grid.newpage()
  r <- lay$rects
  grid::grid.rect(x = bp(r$x), y = bp(H - r$y), width = bp(r$w), height = bp(r$h),
                  just = c("left", "top"), gp = grid::gpar(fill = r$fill, col = NA))
  b <- lay$bar
  grid::grid.rect(x = bp(b$x), y = bp(H - b$y), width = bp(b$w), height = bp(b$h),
                  just = c("left", "top"),
                  gp = grid::gpar(col = NA, fill = grid::linearGradient(
                    colours = lay$stops$colour, stops = lay$stops$offset,
                    x1 = 0.5, y1 = 0, x2 = 0.5, y2 = 1)))
  tk <- lay$bticks
  grid::grid.segments(x0 = bp(tk$x0), y0 = bp(H - tk$y0), x1 = bp(tk$x1), y1 = bp(H - tk$y1),
                      gp = grid::gpar(col = TICK_COL, lwd = lwd(TICK_LWD), lineend = "butt"))
  z <- lay$zeros
  grid::grid.segments(x0 = bp(z$x0), y0 = bp(H - z$y), x1 = bp(z$x1), y1 = bp(H - z$y),
                      gp = grid::gpar(col = NULL_COL, lwd = lwd(NULL_PT), lineend = "butt"))
  a <- lay$axes
  grid::grid.segments(x0 = bp(a$xa), y0 = bp(H - a$ya0), x1 = bp(a$xa), y1 = bp(H - a$ya1),
                      gp = grid::gpar(col = "black", lwd = lwd(AX_LWD_PT), lineend = "butt"))
  yt <- lay$yticks
  grid::grid.segments(x0 = bp(yt$x0), y0 = bp(H - yt$y), x1 = bp(yt$x1), y1 = bp(H - yt$y),
                      gp = grid::gpar(col = "black", lwd = lwd(AX_LWD_PT), lineend = "butt"))
  d <- lay$dots
  grid::grid.segments(x0 = bp(d$xd), y0 = bp(H - d$y_lo), x1 = bp(d$xd), y1 = bp(H - d$y_hi),
                      gp = grid::gpar(col = "black", lwd = lwd(CI_PT), lineend = "butt"))
  circ <- dplyr::bind_rows(d %>% dplyr::select(xd, yd, filled), lay$keys %>% dplyr::select(xd, yd, filled))
  r_path <- (CIRCLE_MM * MM - CIRC_LW_PT) / 2
  grid::grid.circle(x = bp(circ$xd), y = bp(H - circ$yd), r = bp(r_path),
                    gp = grid::gpar(fill = ifelse(circ$filled, "black", "white"), col = "black",
                                    lwd = lwd(CIRC_LW_PT)))
  tx <- lay$texts
  lab <- if (identical(variant, "pdf")) tx$label_pdf else tx$label_uni
  for (i in seq_len(nrow(tx))) {
    vj <- if (tx$rot[i] == 90) 0 else 0
    grid::grid.text(lab[i], x = bp(tx$x[i]), y = bp(H - tx$y[i]),
                    hjust = if (tx$rot[i] == 90) 0.5 else unname(HJ[tx$anchor[i]]),
                    vjust = vj, rot = tx$rot[i],
                    gp = grid::gpar(fontfamily = FONT, fontsize = tx$size[i],
                                    fontface = tx$face[i], col = tx$ink[i]))
  }
  invisible(TRUE)
}
.write_pdf <- function(lay, file, title) {
  grDevices::pdf(file, width = lay$size[["W"]] / 72, height = lay$size[["H"]] / 72,
                 family = FONT, useDingbats = FALSE, compress = FALSE, title = title)
  on.exit(grDevices::dev.off(), add = TRUE)
  .draw(lay, "pdf")
}
.write_png <- function(lay, file) {
  grDevices::png(file, width = lay$size[["W"]] / 72, height = lay$size[["H"]] / 72,
                 units = "in", res = PNG_DPI, bg = "white", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  .draw(lay, "uni")
}
FILES <- c(
  png_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM, ".png")),
  png_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM, ".png")),
  pdf_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM, ".pdf")),
  pdf_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM, ".pdf")))
for (f in FILES[c("png_outputs", "png_docs")]) .write_png(LAY, f)
for (f in FILES[c("pdf_outputs", "pdf_docs")]) .write_pdf(LAY, f, FIG_STEM)

.md5 <- function(f) unname(tools::md5sum(f))
.pdf_md5 <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  s <- rawToChar(b[b != as.raw(0)])
  s <- gsub("(/(CreationDate|ModDate) \\(D:)[0-9]+", "\\10", s, useBytes = TRUE)
  digest::digest(s, algo = "md5", serialize = FALSE)
}
.rgb3 <- function(h) {
  m <- grDevices::col2rgb(h) / 255
  sprintf("%.3f %.3f %.3f", m[1, ], m[2, ], m[3, ])
}
# The PDF, read back INSIDE a function (E51's lesson): rectangles with their
# fills (E50's reader), every string with its font, size, angle and place, the
# straight two-point lines, and the circles. Only ASCII and numbers come back.
.read_pdf <- function(file) {
  b <- readBin(file, "raw", file.info(file)$size)
  b <- b[b != as.raw(0)]
  b[b == as.raw(0xAD)] <- as.raw(0x7E)
  b[as.integer(b) > 127L] <- as.raw(0x3F)
  txt <- rawToChar(b)
  ln <- trimws(strsplit(txt, "\n", fixed = TRUE)[[1]])
  fm <- regmatches(ln, regexec("/Name /F([0-9]+) /BaseFont /([A-Za-z-]+)", ln))
  fm <- do.call(rbind, fm[lengths(fm) == 3L])
  fonts <- stats::setNames(fm[, 3], paste0("F", fm[, 2]))
  mb <- regmatches(txt, regexpr("/MediaBox \\[[^]]*\\]", txt))
  mb_num <- as.numeric(strsplit(gsub("[^0-9. ]", "", mb), " +")[[1]])
  mb_num <- mb_num[!is.na(mb_num)]
  str_re <- "\\((?:[^()\\\\]|\\\\.)*\\)"
  fill <- NA_character_
  rects <- list()
  texts <- list()
  for (i in seq_along(ln)) {
    l <- ln[i]
    m <- regmatches(l, regexec("^([0-9.]+) ([0-9.]+) ([0-9.]+) scn$", l))[[1]]
    if (length(m) == 4L) fill <- paste(m[2:4], collapse = " ")
    if (grepl("/Pattern cs", l, fixed = TRUE)) fill <- "pattern"
    m <- regmatches(l, regexec("^([0-9.-]+) ([0-9.-]+) ([0-9.-]+) ([0-9.-]+) re$", l))[[1]]
    if (length(m) == 5L && i < length(ln) && ln[i + 1L] == "f") {
      rects[[length(rects) + 1L]] <- tibble::tibble(
        x = as.numeric(m[2]), y = as.numeric(m[3]), w = as.numeric(m[4]),
        h = as.numeric(m[5]), fill = fill)
    }
    m <- regmatches(l, regexec(paste0(
      "^/F([0-9]+) 1 Tf ([0-9.-]+) ([0-9.-]+) [0-9.-]+ [0-9.-]+ ([0-9.-]+) ",
      "([0-9.-]+) Tm (.*) T[jJ]$"), l))[[1]]
    if (length(m) == 7L) {
      parts <- regmatches(m[7], gregexpr(str_re, m[7], perl = TRUE))[[1]]
      body <- gsub("\\\\([()\\\\])", "\\1", sub("^\\((.*)\\)$", "\\1", parts))
      a <- as.numeric(m[3]); bb <- as.numeric(m[4])
      texts[[length(texts) + 1L]] <- tibble::tibble(
        font = unname(fonts[paste0("F", m[2])]), size = sqrt(a^2 + bb^2),
        angle = atan2(bb, a) * 180 / pi, x = as.numeric(m[5]), y = as.numeric(m[6]),
        label = paste(body, collapse = ""))
    }
  }
  num <- "(?<![0-9.-])(-?[0-9.]+)"
  sg <- regmatches(txt, gregexpr(paste0(num, " (-?[0-9.]+) m\\s+(-?[0-9.]+) (-?[0-9.]+) l\\s+S"),
                                 txt, perl = TRUE))[[1]]
  segments <- matrix(vapply(strsplit(sg, "\\s+"), function(v) as.numeric(v[c(1, 2, 4, 5)]),
                            numeric(4)), ncol = 4L, byrow = TRUE,
                     dimnames = list(NULL, c("x0", "y0", "x1", "y1")))
  circ_re <- paste0("(-?[0-9.]+) (-?[0-9.]+) m\\s+((?:-?[0-9.]+ -?[0-9.]+ -?[0-9.]+ ",
                    "-?[0-9.]+ -?[0-9.]+ -?[0-9.]+ c\\s+){4})B")
  cs <- regmatches(txt, gregexpr(circ_re, txt, perl = TRUE))[[1]]
  circles <- do.call(rbind, lapply(cs, function(s) {
    v <- as.numeric(regmatches(s, gregexpr("-?[0-9.]+", s))[[1]])
    c(d = max(v[c(TRUE, FALSE)]) - min(v[c(TRUE, FALSE)]))
  }))
  list(mediabox = mb_num, fonts = fonts, rects = dplyr::bind_rows(rects),
       texts = dplyr::bind_rows(texts), segments = segments, circles = circles,
       n_images = sum(grepl("/Subtype /Image", ln, fixed = TRUE)),
       n_type3 = sum(grepl("/Subtype /Type3", ln, fixed = TRUE)),
       n_shadings = sum(grepl("/ShadingType 2", ln, fixed = TRUE)),
       font_file = any(grepl("/FontFile", ln, fixed = TRUE)))
}

.assert("the two PNGs are identical", .md5(FILES[["png_outputs"]]) == .md5(FILES[["png_docs"]]),
        n = 2L)
.assert("the two PDFs are identical, timestamps aside",
        .pdf_md5(FILES[["pdf_outputs"]]) == .pdf_md5(FILES[["pdf_docs"]]), n = 2L)
W <- LAY$size[["W"]]
H <- LAY$size[["H"]]
P <- .read_pdf(FILES[["pdf_docs"]])
pc <- P$rects %>% dplyr::filter(fill != "pattern")
.assert(paste0("in the PDF, 240 cells, each in its size, colour and place, and the ",
               "colour bar's gradient"),
        nrow(pc) == 240L && all(abs(pc$w - LAY$CW) < 0.006) && all(abs(pc$h - CELL_H) < 0.006) &&
          identical(pc$fill, .rgb3(r$fill)) && all(abs(pc$x - r$x) < 0.006) &&
          all(abs(pc$y - (H - r$y - r$h)) < 0.006) &&
          sum(P$rects$fill == "pattern") == 1L && P$n_shadings == 1L, n = nrow(pc))
pt <- P$texts
exp_font <- ifelse(tx$face == 3L, "Helvetica-Oblique", "Helvetica")
exp_lab <- gsub(SOFT_HY, "~", tx$label_pdf, fixed = TRUE)
.assert(paste0("in the PDF every label is live Helvetica text (no image, no embedded ",
               "or Type 3 font), at its size and angle, oblique for the gene symbols only"),
        nrow(pt) == nrow(tx) && identical(pt$label, exp_lab) && identical(pt$font, exp_font) &&
          all(abs(pt$size - tx$size) < 0.01) && all(abs(pt$angle - tx$rot) < 0.05) &&
          P$n_images == 0L && P$n_type3 == 0L && !P$font_file, n = nrow(pt))
.assert("the PDF's page is the panel's size, in whole points, within 180 mm",
        identical(as.integer(P$mediabox[3:4]), as.integer(floor(c(W, H)))) &&
          P$mediabox[3] * 25.4 / 72 <= 180, n = 1L, detail = paste(P$mediabox, collapse = " "))
CI_EXP <- sort(abs(LAY$dots$y_hi - LAY$dots$y_lo))
SG <- P$segments
v_len <- abs(SG[abs(SG[, "x1"] - SG[, "x0"]) < 1e-9, "y1"] - SG[abs(SG[, "x1"] - SG[, "x0"]) < 1e-9, "y0"])
ci_pdf <- sort(v_len[v_len < 0.95 * S_H & v_len > 0])
dev_ci <- if (length(ci_pdf) == 32L) max(abs(ci_pdf - CI_EXP)) else Inf
.assert("in the PDF the strip's 32 intervals have their lengths",
        is.finite(dev_ci) && dev_ci <= 0.02, n = 32L, max_dev = dev_ci)
circ_outer_mm <- (P$circles[, "d"] + CIRC_LW_PT) * 25.4 / 72
.assert("in the PDF the 34 circles, 32 means and the key's 2, are 1.05 mm across, outline included",
        !is.null(P$circles) && nrow(P$circles) == 34L && all(abs(circ_outer_mm - CIRCLE_MM) <= 0.01),
        n = 34L, max_dev = max(abs(circ_outer_mm - CIRCLE_MM)))
img <- png::readPNG(FILES[["png_docs"]])
s_px <- PNG_DPI / 72
px_x <- floor((r$x + 1.5) * s_px) + 1L
px_y <- floor((r$y + 1.5) * s_px) + 1L
px_col <- grDevices::rgb(round(img[cbind(px_y, px_x, 1L)] * 255),
                         round(img[cbind(px_y, px_x, 2L)] * 255),
                         round(img[cbind(px_y, px_x, 3L)] * 255), maxColorValue = 255)
.assert(paste0("the PNG is the panel at ", PNG_DPI, " dpi, and every cell shows its ",
               "colour exactly"),
        abs(dim(img)[2] - W / 72 * PNG_DPI) <= 1 && abs(dim(img)[1] - H / 72 * PNG_DPI) <= 1 &&
          identical(px_col, r$fill), n = nrow(r),
        detail = sprintf("%d x %d px", dim(img)[2], dim(img)[1]))
rm(img)
FIGURE_FILES <- tibble::tibble(
  copy = names(FILES), file = basename(FILES),
  checksum = c(.md5(FILES[["png_outputs"]]), .md5(FILES[["png_docs"]]),
               .pdf_md5(FILES[["pdf_outputs"]]), .pdf_md5(FILES[["pdf_docs"]])),
  checksum_is = c("md5", "md5", "md5, timestamps zeroed", "md5, timestamps zeroed"))

# =============================================================================
# 8. The legend draft, Claude Code's, for the author
# =============================================================================
LEGEND_DRAFT <- paste0(
  "Fifteen BCL-2 family genes across the four MYC x OXPHOS groups of tumours, in ",
  "TCGA-BRCA (n = 938) and SCAN-B (n = 3,143). Left, the OXPHOS score by GSVA; ",
  "right, by mitoPPS; MYC activity is the same in both. Columns: within each ",
  "cohort, OXPHOS-low then OXPHOS-high tumours, MYC-low then MYC-high within ",
  "each. Cells, z-scores: adjusted group means of Blom-transformed expression in ",
  "within-cohort standard deviations, from a linear model with group, PAM50 ",
  "subtype, tumour purity and leukocyte fraction in TCGA-BRCA, and group and ",
  "PAM50 subtype in SCAN-B. The three effectors (BAX, BAK1, BOK) were fitted ",
  "for this figure on the same model, post hoc. Rows and colours as in Fig. 4g. ",
  "Strip, the signature: filled, the same model; open, proliferation adjusted; ",
  "bars, 95% confidence intervals.")
.assert("the legend draft is ASCII, as the script is", !grepl("[^ -~]", LEGEND_DRAFT), n = 1L)

# =============================================================================
# 9. Saved: the object, three tables and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

SPEC <- list(
  ruling = "docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md (dd5004a)",
  declaration = "docs/2026-11-06_E56_effector_quadrant_means_declaration.md (922f2ad)",
  sources = list(e34 = E34_MD5, e46 = E46_DIGESTS),
  model = "E34's m1: quad + PAM50 (+ purity + leuko in TCGA); Blom-scored log2(linear + 1)",
  genes = FIFTEEN, groups = GROUPS, effectors = EFFECTORS, columns = COL_Q,
  scale = L_SCALE, palette = DIV_COLS, bar_title = BAR_TITLE, key = KEY_TEXT,
  page_pt = LAY$size, column_pt = LAY$CW, strip_axis = LAY$YLIM, circle_mm = CIRCLE_MM)
GATES <- list(g1_max_dev = dev_g1, g2 = G2J, g2_max_dev = dev_g2, g3 = G3J, g3_max_dev = dev_g3)
OUT <- list(
  effector_means = EFF,
  gates = GATES,
  values = VAL,
  frame = lapply(FRX, function(f) f %>% dplyr::select(sample_id, quad_gsva, quad_mitopps,
                                                     PAM50, purity, leuko, PROLIF,
                                                     dplyr::starts_with("raw_g_"),
                                                     dplyr::all_of(paste0("g_", EFFECTORS)))),
  coefficients = lapply(EFF_FITS, function(f) list(b = f$b, V = f$V, df = f$df, n = f$n,
                                                   cell_n = f$cell_n)),
  fit_digests = FIT_DIGESTS,
  layout = LAY,
  legend_draft = LEGEND_DRAFT,
  assertions = ASSERTIONS,
  figure_files = FIGURE_FILES,
  spec = SPEC,
  fits = EFF_LM,
  run_info = list(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  R = R.version.string,
                  bitmap_type_default = getOption("bitmapType")))
saveRDS(OUT, PATH_E56)
utils::write.csv(EFF, file.path(OUT_TABLES, "E56_effector_quadrant_means.csv"), row.names = FALSE)
utils::write.csv(VAL, file.path(OUT_TABLES, "E56_heatmap_values.csv"), row.names = FALSE)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E56_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), c("run_info", "fits")),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info and the lm objects excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-15s %s", k, DIGESTS[[k]]))
message("\n   the legend draft (Claude Code's):\n   ", LEGEND_DRAFT)
message("\nE56 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figure.
  library(dplyr)
  x <- readRDS(here::here("results", "e56_quadrant_configuration_fifteen_genes.rds"))
  x$effector_means %>% dplyr::select(gene, cohort, instrument, quadrant, est, lo, hi) %>% print(n = 48)
  x$gates[c("g1_max_dev", "g2_max_dev", "g3_max_dev")] %>% str()
  x$assertions %>% print(n = Inf)
  cat(x$legend_draft, "\n")
  img <- png::readPNG(here::here("docs", "figures",
                                 "2026-11-06_E56_quadrant_configuration_fifteen_genes.png"))
  grid::grid.newpage()
  grid::grid.raster(img)
}
