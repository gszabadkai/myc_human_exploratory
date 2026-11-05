# E54_myc_coamplification_panel.R
# =============================================================================
# FIG. 4x: BCL2L1 GAIN AND MCL1 AMPLIFICATION WITH MYC AMPLIFICATION IN
# TCGA-BRCA, AS ADJUSTED ODDS RATIOS: ALL TUMOURS, THEN BY PAM50 SUBTYPE.
#
# A DISPLAY SCRIPT. It fits nothing. It draws the validation study's saved G2
# result (myc_human_validation, gate G2, 2026-08-28), snapshotted into
# data/from_validation_g2/ with its provenance (c1dceed). It checks the
# snapshot's md5, and checks every number it draws that the committed G2
# result note prints (myc_human_validation d3ac60e,
# docs/2026-08-28_G2_result.md) at the note's precision.
#
# Specification: the author's ruling of 2026-11-05, verbatim in
# docs/2026-11-05_E54_coamplification_panel_ruling.md, committed alone at
# 8fde82a before this script was written. It amends decision 5 of 2026-11-02:
# the co-amplification result becomes a main panel, Fig. 4x.
#
# AMENDED after the first run (verified, ec9b590) by the author's ruling of
# 2026-11-05, docs/2026-11-05_E54_point_size_ruling.md (6550011): the circles
# are 1.05 mm across, outline included, in place of 1.6 mm. Nothing else
# changes.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNG is compared byte for
# byte with the dry run.
#
# =============================================================================
# WHAT IT DRAWS, AND WHY
# =============================================================================
#   - Rows, for BCL2L1 gain (GISTIC >= +1) and then MCL1 amplification (+2),
#     each against MYC amplification (+2), as G2 defined them (the ruling,
#     question 2):
#     - all tumours, the pre-registered primary test, under both pre-specified
#       covariates: aneuploidy score (filled) and fraction genome altered
#       (open). G2 says BCL2L1's result must never be shown as a bare odds
#       ratio, because its interval crosses 1 under the second;
#     - Basal, Her2, LumA and LumB, adjusted for aneuploidy score, which G2
#       declared secondary and exploratory. The Normal stratum (n = 35) is
#       left out, as the ruling's preview said.
#   - The x axis is 26 mm long (the ruling, question 1), on a log2 scale from
#     0.125 to 8, so that 1 sits in the middle. Ticks are labelled at 0.25, 1
#     and 4, with unlabelled ticks at each doubling between.
#   - A key above the data area names the two covariates, as E52's key names
#     the cohorts.
#   - The marks are TCGA-BRCA's dark slate (#243642, E51's ruling); the line at
#     1 is E48's grey. Type is 5 pt, with 6 pt for the gene headers and the
#     axis title.
#   - The circles are 1.05 mm across, outline included (the point-size ruling):
#     size 1.0 with a 0.3 outline, a little heavier than E52's, so that the
#     open circles still read as open. The PDF read-back measures them.
#   - The PDF carries live Helvetica, for Illustrator; the PNG is a preview.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing, recomputes nothing and reads no expression data.
#   - It reads nothing across repos: only this repo's snapshot.
#   - It changes no text of the paper. The author's sentence waits on a text
#     ruling (the ruling, question 4). Its legend draft is Claude Code's, for
#     that ruling.
#   - It leaves no raw file bytes in the workspace: the PDF is read back inside
#     a function that returns only ASCII and numbers (E51's lesson).
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
for (pkg in c("digest", "png", "patchwork", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E54 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE54: Fig. 4x, co-amplification with MYC in TCGA-BRCA\n", strrep("=", 78))
message("Draws G2's saved result and fits nothing. Specification 8fde82a.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_G2         <- here::here("data", "from_validation_g2", "g2_cnv_cooccurrence.rds")
PATH_RULING     <- here::here("docs", "2026-11-05_E54_coamplification_panel_ruling.md")
PATH_RULING_SIZE <- here::here("docs", "2026-11-05_E54_point_size_ruling.md")
PATH_E54        <- file.path(DIR_RESULTS, "e54_myc_coamplification_panel.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-05"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- "E54_myc_coamplification_odds_ratios"

# =============================================================================
# 0. Constants
# =============================================================================
G2_MD5 <- "33da579ae1777c6bd246f5e8518d967d"   # data/from_validation_g2/README.md

# What the committed G2 result note prints for every number drawn here (its
# sections 2, 5 and 7). Odds ratios are printed to two decimals and interval
# bounds to three significant figures. The note gives intervals for the pooled
# rows and for Basal only.
NOTE <- tibble::tribble(
  ~gene,    ~row,          ~covariate,   ~or,  ~lo,   ~hi,
  "BCL2L1", "All tumours", "aneuploidy", 1.77, 1.26,  2.47,
  "BCL2L1", "All tumours", "fga",        1.36, 0.984, 1.88,
  "BCL2L1", "Basal",       "aneuploidy", 2.69, 1.33,  5.62,
  "BCL2L1", "Her2",        "aneuploidy", 0.99, NA,    NA,
  "BCL2L1", "LumA",        "aneuploidy", 1.70, NA,    NA,
  "BCL2L1", "LumB",        "aneuploidy", 0.90, NA,    NA,
  "MCL1",   "All tumours", "aneuploidy", 1.69, 1.17,  2.44,
  "MCL1",   "All tumours", "fga",        1.60, 1.10,  2.31,
  "MCL1",   "Basal",       "aneuploidy", 3.12, 1.53,  6.56,
  "MCL1",   "Her2",        "aneuploidy", 1.22, NA,    NA,
  "MCL1",   "LumA",        "aneuploidy", 0.68, NA,    NA,
  "MCL1",   "LumB",        "aneuploidy", 1.12, NA,    NA)

GENES <- tibble::tribble(
  ~gene,    ~rule, ~header,
  "BCL2L1", "ge1", 'italic("BCL2L1")~"gain"',
  "MCL1",   "eq2", 'italic("MCL1")~"amplification"')
SUBTYPES <- c("Basal", "Her2", "LumA", "LumB")

# --- size and type -----------------------------------------------------------------
AXIS_MM  <- 26                  # the x axis (the ruling, question 1)
X_LIM    <- c(0.125, 8)         # log2: 1 in the middle
X_BREAKS <- c(0.25, 1, 4)
X_MINOR  <- c(0.125, 0.5, 2, 8)
PITCH    <- 2.2                 # mm per row slot
PT_EPS   <- 0.001               # pdf() truncates the page box to whole points
PT_TXT   <- 5
PT_TITLE <- 6
CAP_EM   <- 0.718               # Helvetica's capitals above the baseline, in em (_style.R)
COL_MARK <- "#243642"           # TCGA-BRCA's dark slate (E51's ruling, question 3)
NULL_COL <- "grey45"; NULL_LW <- 0.3      # E48's reference line
CI_LW    <- 0.4                           # E48's intervals
PT_SIZE  <- 1.0; PT_STROKE <- 0.3         # the point-size ruling (6550011)
# A circle's width, outline included, in mm: R draws pch 21 with a diameter of
# 0.75 x its font size, which ggplot2 sets at size x .pt + stroke x .stroke / 2,
# and an outline of stroke x .stroke / 2 line widths, each 0.75 pt.
CIRCLE_MM <- (0.75 * (PT_SIZE * ggplot2::.pt + PT_STROKE * ggplot2::.stroke / 2) +
                0.75 * PT_STROKE * ggplot2::.stroke / 2) * 25.4 / 72
M_EDGE   <- 1                   # mm: the page's margin
INDENT   <- 1.5                 # mm: a subtype's label under its gene's header
G_LAB    <- 1.5                 # mm: the longest row label to the data area
KEY_ABOVE <- 1.2                # mm: the key's lower baseline above the data area (E52)
KEY_STEP <- 2.8                 # mm: the key's rows, baseline to baseline (E52)
KEY_GAP  <- 0.8                 # mm: a key symbol to its text
KEY_TEXT <- c(aneuploidy = "adjusted for aneuploidy score",
              fga = "adjusted for fraction genome altered")
X_TITLE  <- "Odds ratio for co-occurrence\nwith MYC amplification"
TITLE_GAP <- 2.4                # mm, baseline to baseline (the Prism format; E52)
TITLE_LH <- TITLE_GAP / (1.2 * PT_TITLE * 25.4 / 72)   # E52's finding: R spaces lines 1.2 x
PNG_DPI  <- 600
TOL_PT   <- 0.02                # the PDF writes coordinates to 0.01 pt
SOFT_HY  <- intToUtf8(173L)     # the PDF's hyphen: its encoding prints "-" as minus

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

.assert("the panel is the ruling's: a main panel with a 26-mm axis, as its own file records",
        file.exists(PATH_RULING) && AXIS_MM == 26 &&
          any(grepl("Main panel, Fig. 4x", readLines(PATH_RULING, warn = FALSE),
                    fixed = TRUE)), n = 1L)
.assert("the circles are the point-size ruling's, 1.05 mm across, as its own file records",
        file.exists(PATH_RULING_SIZE) && round(CIRCLE_MM, 2) == 1.05 &&
          any(grepl("1.05 mm across", readLines(PATH_RULING_SIZE, warn = FALSE),
                    fixed = TRUE)), n = 1L, detail = sprintf("%.3f mm", CIRCLE_MM))

# =============================================================================
# 1. The snapshot, checked by its md5
# =============================================================================
message("\n1. G2's saved result")
.stop_if(file.exists(PATH_G2), "the G2 snapshot is missing: see ",
         "data/from_validation_g2/README.md. STOP.")
.assert("the G2 snapshot is the one recorded in its README, by md5",
        identical(unname(tools::md5sum(PATH_G2)), G2_MD5), n = 1L)
g2 <- readRDS(PATH_G2)

# =============================================================================
# 2. The rows drawn, from the snapshot, checked against the committed note
# =============================================================================
message("\n2. The rows")
CO <- as.data.frame(g2$cooccurrence)
ST <- as.data.frame(g2$stratified)
.one <- function(d, what) {
  .stop_if(nrow(d) == 1L, "not exactly one row: ", what)
  d
}
ROWS <- dplyr::bind_rows(lapply(seq_len(nrow(GENES)), function(i) {
  g <- GENES$gene[i]; r <- GENES$rule[i]
  pooled <- function(tier, cov) {
    d <- .one(CO[CO$source == "A_ISAR" & CO$tier == tier & CO$partner == g &
                   CO$partner_rule == r & CO$exposure_rule == "eq2", ],
              paste(g, tier))
    tibble::tibble(gene = g, rule = r, row = "All tumours", covariate = cov,
                   n = d$n, n_both = d$n_both, or = d$adj_or, lo = d$adj_lo,
                   hi = d$adj_hi, p = d$adj_p, aneu_covariate = d$aneu_covariate,
                   tier = tier)
  }
  strata <- dplyr::bind_rows(lapply(SUBTYPES, function(z) {
    d <- .one(ST[ST$strat_var == "pam50" & ST$stratum == paste0("BRCA_", z) &
                   ST$partner == g & ST$partner_rule == r, ], paste(g, z))
    tibble::tibble(gene = g, rule = r, row = z, covariate = "aneuploidy",
                   n = d$n, n_both = d$n_both, or = d$adj_or, lo = d$adj_lo,
                   hi = d$adj_hi, p = d$adj_p, aneu_covariate = d$aneu_covariate,
                   tier = "pam50 stratum")
  }))
  dplyr::bind_rows(pooled("primary", "aneuploidy"),
                   pooled("sensitivity_fga", "fga"), strata)
}))
.assert(paste0("the rows are the ruling's: for each gene at G2's threshold, all ",
               "tumours under both covariates, then the four PAM50 subtypes"),
        nrow(ROWS) == 12L && all(ROWS$n[ROWS$row == "All tumours"] == 1043L) &&
          all(ROWS$aneu_covariate[ROWS$covariate == "aneuploidy"] == "ANEUPLOIDY_SCORE") &&
          all(ROWS$aneu_covariate[ROWS$covariate == "fga"] == "FRACTION_GENOME_ALTERED"),
        n = 12L, detail = paste(unique(ROWS$aneu_covariate), collapse = ", "))
.assert("every interval contains its estimate, and all are finite",
        all(is.finite(c(ROWS$or, ROWS$lo, ROWS$hi))) &&
          all(ROWS$lo < ROWS$or & ROWS$or < ROWS$hi), n = 12L)

CHK <- dplyr::left_join(NOTE, ROWS, by = c("gene", "row", "covariate"),
                        suffix = c("_note", "")) %>%
  dplyr::mutate(or_ok = round(or, 2) == or_note,
                lo_ok = is.na(lo_note) | signif(lo, 3) == lo_note,
                hi_ok = is.na(hi_note) | signif(hi, 3) == hi_note)
n_checked <- nrow(CHK) + sum(!is.na(CHK$lo_note)) + sum(!is.na(CHK$hi_note))
.assert(paste0("every number drawn that the committed G2 note prints matches it, ",
               "at the note's precision"),
        nrow(CHK) == 12L && all(CHK$or_ok & CHK$lo_ok & CHK$hi_ok), n = n_checked,
        detail = paste(CHK$gene, CHK$row, CHK$covariate)[!(CHK$or_ok & CHK$lo_ok & CHK$hi_ok)] %>%
          paste(collapse = "; "))
.assert("the axis, 0.125 to 8, holds every interval",
        min(ROWS$lo) > X_LIM[1] && max(ROWS$hi) < X_LIM[2], n = 12L,
        detail = sprintf("%.3f to %.2f", min(ROWS$lo), max(ROWS$hi)))
for (i in seq_len(nrow(ROWS))) {
  r <- ROWS[i, ]
  message(sprintf("   %-6s %-11s %-10s n = %4d: OR %.2f [%.2f, %.2f]", r$gene, r$row,
                  r$covariate, r$n, r$or, r$lo, r$hi))
}

# =============================================================================
# 3. Fig. 1's style helpers, copied unchanged from myc_mouse 2ebf646,
#    figures/natmetab/_style.R
# =============================================================================
NM_FONT      <- "Helvetica"  # throughout: R panels, Prism graphs, Illustrator labels
NM_TXT_MIN   <- 5            # pt; the floor -- nothing on a panel may be smaller
NM_TXT       <- getOption("natmetab.txt", 6)        # pt; tick labels, keys, on-panel labels, P values
NM_TXT_TITLE <- getOption("natmetab.txt_title", 7)  # pt; axis titles and strips (Nature's maximum)
NM_LINE      <- 0.25         # ggplot linewidth: axes, ticks, brackets, box outlines (0.53 pt printed)

pt <- function(x) x / ggplot2::.pt   # points -> ggplot text size units

lab_num <- function(x) {
  out <- format(x, big.mark = ",", scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  out[is.na(x)] <- NA_character_
  out
}

theme_nm <- function(legend = "none", txt = NM_TXT, txt_title = NM_TXT_TITLE, line = NM_LINE) {
  stopifnot(txt >= NM_TXT_MIN, txt_title >= txt, txt_title <= 7)
  ggplot2::theme_classic(base_size = txt_title, base_family = NM_FONT) +
    ggplot2::theme(
      text              = ggplot2::element_text(colour = "black"),
      axis.text         = ggplot2::element_text(size = txt, colour = "black"),
      axis.title        = ggplot2::element_text(size = txt_title, colour = "black"),
      axis.line         = ggplot2::element_line(linewidth = line, colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = line, colour = "black"),
      axis.ticks.length = ggplot2::unit(1, "mm"),
      axis.minor.ticks.length = ggplot2::rel(0.5),
      legend.text       = ggplot2::element_text(size = txt),
      legend.title      = ggplot2::element_text(size = txt),
      legend.key.size   = ggplot2::unit(3, "mm"),
      legend.key.spacing.y = ggplot2::unit(0.4, "mm"),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(2, "mm"),     # key to plot; the 7-pt default is 2.5 mm
      legend.position   = legend,
      strip.text        = ggplot2::element_text(size = txt_title),
      strip.background  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_blank(),
      plot.subtitle     = ggplot2::element_blank(),
      plot.caption      = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(1, 2.5, 1, 1, "mm"))  # right: end labels overhang the axis
}

# =============================================================================
# 4. The panel
# =============================================================================
message("\n4. The panel")
# The printed size of a label at size_pt, in mm, with the PDF device's
# Helvetica metrics: width, height above its baseline, and descent.
.txt_mm <- function(labels, size_pt, parse = FALSE) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  t(vapply(labels, function(s) {
    lab <- if (parse) parse(text = s) else s
    g <- grid::textGrob(lab, gp = grid::gpar(fontsize = size_pt, fontfamily = NM_FONT))
    c(w = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
      h = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE),
      d = grid::convertHeight(grid::grobDescent(g), "mm", valueOnly = TRUE))
  }, numeric(3)))
}
CAP_MM <- function(size_pt) CAP_EM * size_pt * 25.4 / 72
M_ROW  <- .txt_mm(c("All tumours", SUBTYPES), PT_TXT)
M_HEAD <- .txt_mm(GENES$header, PT_TITLE, parse = TRUE)
M_KEY  <- .txt_mm(KEY_TEXT, PT_TXT)

# Row slots, from the top: a gene's header, its two pooled rows, its four
# subtypes; then the next gene. y is minus the slot's number.
SLOTS <- dplyr::bind_rows(lapply(seq_len(nrow(GENES)), function(i) {
  k0 <- (i - 1) * 7
  tibble::tibble(gene = GENES$gene[i],
                 row = c("header", "All tumours", "All tumours", SUBTYPES),
                 covariate = c(NA, "aneuploidy", "fga", rep("aneuploidy", 4)),
                 slot = k0 + 1:7)
})) %>% dplyr::mutate(y = -slot)
N_SLOT <- nrow(SLOTS)
PANEL_H <- N_SLOT * PITCH                               # mm
Y_LIM <- c(-N_SLOT - 0.5, -0.5)
FR <- dplyr::inner_join(ROWS, SLOTS, by = c("gene", "row", "covariate")) %>%
  dplyr::mutate(filled = covariate == "aneuploidy")
.assert("every row drawn has its slot, and every slot but the headers has one row",
        nrow(FR) == 12L && setequal(FR$slot, SLOTS$slot[SLOTS$row != "header"]), n = 12L)

# Across, in mm from the page's left: headers at the margin, subtype labels
# indented, then the data area.
X_HEAD <- M_EDGE
X_ITEM <- M_EDGE + INDENT
X_DATA <- X_ITEM + max(M_ROW[, "w"]) + G_LAB            # the data area's left edge
MM_PER_DOUBLING <- AXIS_MM / log2(X_LIM[2] / X_LIM[1])
.x_at <- function(mm) X_LIM[1] * 2^((mm - X_DATA) / MM_PER_DOUBLING)   # page mm -> odds ratio
X_NULL_MM <- X_DATA + log2(1 / X_LIM[1]) * MM_PER_DOUBLING
.assert(paste0("each gene's header, which runs into its empty slot, ends at least ",
               "1 mm before the line at 1"),
        all(X_HEAD + M_HEAD[, "w"] <= X_NULL_MM - 1), n = nrow(GENES),
        detail = sprintf("%.1f mm against %.1f", max(X_HEAD + M_HEAD[, "w"]), X_NULL_MM))

# Down: label baselines with their capitals centred on their slots. A plain
# string's baseline sits at y with vjust = 0 (E53's finding); a plotmath
# expression's lowest point does, so a header is drawn its descent lower (E52).
.mm_y <- function(mm) mm / PITCH                        # mm -> slot units
LABS <- dplyr::bind_rows(
  tibble::tibble(kind = "header", label = GENES$header, parse = TRUE, size_pt = PT_TITLE,
                 x_mm = X_HEAD,
                 y = SLOTS$y[SLOTS$row == "header"] - .mm_y(CAP_MM(PT_TITLE) / 2) -
                   .mm_y(unname(M_HEAD[, "d"]))),
  tibble::tibble(kind = "pooled", label = "All tumours", parse = FALSE, size_pt = PT_TXT,
                 x_mm = X_ITEM,
                 y = vapply(GENES$gene, function(g)
                   mean(SLOTS$y[SLOTS$gene == g & SLOTS$row == "All tumours"]), numeric(1)) -
                   .mm_y(CAP_MM(PT_TXT) / 2)),
  tibble::tibble(kind = "subtype", label = SLOTS$row[SLOTS$row %in% SUBTYPES], parse = FALSE,
                 size_pt = PT_TXT, x_mm = X_ITEM,
                 y = SLOTS$y[SLOTS$row %in% SUBTYPES] - .mm_y(CAP_MM(PT_TXT) / 2))) %>%
  dplyr::mutate(x = .x_at(x_mm), size_gg = pt(size_pt))

# The key, above the data area, as E52's: the aneuploidy row on top.
KEY <- tibble::tibble(covariate = c("aneuploidy", "fga"), text = unname(KEY_TEXT),
                      base_mm = KEY_ABOVE + c(1, 0) * KEY_STEP) %>%
  dplyr::mutate(y_text = Y_LIM[2] + .mm_y(base_mm),
                y_mark = Y_LIM[2] + .mm_y(base_mm + CAP_MM(PT_TXT) / 2),
                x_mark = .x_at(M_EDGE + CIRCLE_MM / 2),
                x_text = .x_at(M_EDGE + CIRCLE_MM + KEY_GAP),
                fill = ifelse(covariate == "aneuploidy", COL_MARK, "white"))
TOP_MM <- max(KEY$base_mm) + max(M_KEY[, "h"]) + M_EDGE   # the margin kept above the data area
.assert("the key is on the page and to the left of nothing it could cover",
        M_EDGE + CIRCLE_MM + KEY_GAP + max(M_KEY[, "w"]) <= X_DATA + AXIS_MM, n = 2L,
        detail = sprintf("%.1f mm wide", CIRCLE_MM + KEY_GAP + max(M_KEY[, "w"])))

# The PDF writes the hyphen in "co-occurrence" as a soft hyphen; the PNG as a
# hyphen (E51 to E53).
.x_title <- function(variant) if (variant == "pdf") gsub("-", SOFT_HY, X_TITLE, fixed = TRUE) else X_TITLE
.panel <- function(variant) {
  ggplot2::ggplot(FR) +
    ggplot2::geom_vline(xintercept = 1, colour = NULL_COL, linewidth = NULL_LW) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y, yend = y),
                          colour = COL_MARK, linewidth = CI_LW) +
    ggplot2::geom_point(ggplot2::aes(x = or, y = y, fill = filled), shape = 21,
                        size = PT_SIZE, stroke = PT_STROKE, colour = COL_MARK) +
    ggplot2::geom_point(data = KEY, ggplot2::aes(x = x_mark, y = y_mark, fill = covariate == "aneuploidy"),
                        shape = 21, size = PT_SIZE, stroke = PT_STROKE, colour = COL_MARK,
                        inherit.aes = FALSE) +
    ggplot2::geom_text(data = KEY, ggplot2::aes(x = x_text, y = y_text, label = text),
                       hjust = 0, vjust = 0, size = pt(PT_TXT), family = NM_FONT,
                       colour = "black", inherit.aes = FALSE) +
    ggplot2::geom_text(data = LABS[!LABS$parse, ], ggplot2::aes(x = x, y = y, label = label,
                                                                 size = size_gg),
                       hjust = 0, vjust = 0, family = NM_FONT, colour = "black",
                       inherit.aes = FALSE) +
    ggplot2::geom_text(data = LABS[LABS$parse, ], ggplot2::aes(x = x, y = y, label = label,
                                                                size = size_gg),
                       hjust = 0, vjust = 0, family = NM_FONT, colour = "black", parse = TRUE,
                       inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = COL_MARK, `FALSE` = "white"), guide = "none") +
    ggplot2::scale_size_identity() +
    # oob_keep: the labels and the key sit outside the axes' limits.
    ggplot2::scale_x_continuous(transform = "log2", limits = X_LIM, breaks = X_BREAKS,
                                minor_breaks = X_MINOR, labels = lab_num(X_BREAKS),
                                guide = ggplot2::guide_axis(minor.ticks = TRUE),
                                expand = c(0, 0), oob = scales::oob_keep) +
    ggplot2::scale_y_continuous(limits = Y_LIM, breaks = NULL, expand = c(0, 0),
                                oob = scales::oob_keep) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = .x_title(variant), y = NULL) +
    theme_nm(legend = "none", txt = PT_TXT, txt_title = PT_TITLE) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(lineheight = TITLE_LH),
                   plot.margin = ggplot2::margin(TOP_MM, M_EDGE, M_EDGE, X_DATA, "mm"))
}
# patchwork fixes the data area at 26 x PANEL_H mm; its own outer margin is set
# to 0, so that the plot's margins frame the page (E53's finding).
.wrap <- function(p) {
  patchwork::wrap_plots(p) +
    patchwork::plot_layout(widths = grid::unit(AXIS_MM, "mm"),
                           heights = grid::unit(PANEL_H, "mm")) +
    patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0)))
}
PANEL <- list(pdf = .wrap(.panel("pdf")), png = .wrap(.panel("png")))

B <- ggplot2::ggplot_build(PANEL$pdf[[1]])
geo <- vapply(PANEL$pdf[[1]]$layers, function(l) class(l$geom)[1], character(1))
pts <- B$data[[which(geo == "GeomPoint")[1]]]
ci  <- B$data[[which(geo == "GeomSegment")[1]]]
k1 <- function(x, y) paste(sprintf("%.10f", x), sprintf("%.10f", y))
.assert(paste0("every point and interval drawn is its stored odds ratio and ",
               "interval, in its slot (on the log2 scale)"),
        nrow(pts) == 12L && nrow(ci) == 12L &&
          setequal(k1(pts$x, pts$y), k1(log2(FR$or), FR$y)) &&
          setequal(paste(k1(ci$x, ci$y), sprintf("%.10f", ci$xend)),
                   paste(k1(log2(FR$lo), FR$y), sprintf("%.10f", log2(FR$hi)))) &&
          sum(pts$fill == COL_MARK) == 10L, n = 24L)

.natural <- function(pw) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- patchwork::patchworkGrob(pw)
  c(width  = sum(grid::convertWidth(g$widths, "mm", valueOnly = TRUE)),
    height = sum(grid::convertHeight(g$heights, "mm", valueOnly = TRUE)))
}
NEED <- .natural(PANEL$pdf)
PAGE_PT <- as.integer(ceiling(NEED * 72 / 25.4 - 1e-9))
names(PAGE_PT) <- c("w", "h")
PAGE <- (PAGE_PT + PT_EPS) * 25.4 / 72
over <- PAGE_PT * 25.4 / 72 - NEED
.assert("the page is the panel's natural size, rounded up to whole points",
        all(over >= -1e-9) && all(over < 25.4 / 72), n = 2L,
        detail = sprintf("%.3f and %.3f mm over it", over[1], over[2]))
message(sprintf("   data area %.0f x %.1f mm; page %d x %d pt (%.1f x %.1f mm)", AXIS_MM,
                PANEL_H, PAGE_PT[["w"]], PAGE_PT[["h"]], PAGE_PT[["w"]] * 25.4 / 72,
                PAGE_PT[["h"]] * 25.4 / 72))

# =============================================================================
# 5. The files: the PNG first, then the PDF (base pdf(), live Helvetica)
# =============================================================================
message("\n5. The files")
FILES <- c(
  png_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM, ".png")),
  png_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM, ".png")),
  pdf_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM, ".pdf")),
  pdf_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM, ".pdf")))
for (f in FILES[c("png_outputs", "png_docs")]) {
  ggplot2::ggsave(f, PANEL$png, width = PAGE[["w"]], height = PAGE[["h"]], units = "mm",
                  dpi = PNG_DPI, device = "png", type = "cairo", bg = "white")
}
.pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(filename, width = width, height = height, family = NM_FONT,
                 useDingbats = FALSE, compress = FALSE)
}
for (f in FILES[c("pdf_outputs", "pdf_docs")]) {
  ggplot2::ggsave(f, PANEL$pdf, width = PAGE[["w"]], height = PAGE[["h"]], units = "mm",
                  device = .pdf_device)
}

.md5 <- function(f) unname(tools::md5sum(f))
.pdf_md5 <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  s <- rawToChar(b[b != as.raw(0)])
  s <- gsub("(/(CreationDate|ModDate) \\(D:)[0-9]+", "\\10", s, useBytes = TRUE)
  digest::digest(s, algo = "md5", serialize = FALSE)
}
.pdf_text <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  rawToChar(b[b != as.raw(0)])
}
# The PDF is read back INSIDE a function, so that its raw bytes, which are not
# valid UTF-8, never reach the workspace (E51's lesson). Only ASCII and numbers
# come back: the page box, the fonts, every string the PDF shows with its
# place, size and angle, and every straight two-point line.
.pdf_read_back <- function(f) {
  txt <- .pdf_text(f)
  mb <- regmatches(txt, regexpr("/MediaBox \\[[^]]*\\]", txt, useBytes = TRUE))
  mb_num <- as.numeric(strsplit(gsub("[^0-9. ]", "", mb), " +")[[1]])
  mb_num <- mb_num[!is.na(mb_num)]
  fonts <- unique(regmatches(txt, gregexpr("/BaseFont /[A-Za-z-]+", txt, useBytes = TRUE))[[1]])
  str_re <- "\\((?:[^()\\\\]|\\\\.)*\\)"
  body   <- function(v) sub("^\\((.*)\\)$", "\\1", v, useBytes = TRUE)
  decode <- function(v) {
    v <- gsub("\\\\([()\\\\])", "\\1", v, useBytes = TRUE)
    v <- gsub(rawToChar(as.raw(0xAD)), "<SHY>", v, fixed = TRUE, useBytes = TRUE)
    iconv(v, from = "latin1", to = "ASCII", sub = "byte")
  }
  bt <- regmatches(txt, gregexpr("(?s)(?<![A-Za-z])BT\\s.*?\\sET(?![A-Za-z])", txt,
                                 perl = TRUE, useBytes = TRUE))[[1]]
  # Each font object names its resource and its font (/Name /F4 /BaseFont
  # /Helvetica-Oblique); each text block selects one resource (/F4 1 Tf).
  fm <- regmatches(txt, gregexpr("/Name /F[0-9]+ /BaseFont /[A-Za-z-]+", txt,
                                 useBytes = TRUE))[[1]]
  font_of <- stats::setNames(sub(".*/BaseFont /", "", fm), sub("^/Name /(F[0-9]+) .*$", "\\1", fm))
  num6 <- paste(rep("(-?[0-9.]+)", 6), collapse = " ")
  placed <- do.call(rbind, lapply(bt, function(b) {
    m <- regmatches(b, regexec(paste0(num6, " Tm"), b, useBytes = TRUE))[[1]]
    s <- regmatches(b, gregexpr(str_re, b, perl = TRUE, useBytes = TRUE))[[1]]
    f <- regmatches(b, regexec("/(F[0-9]+) [0-9.]+ Tf", b, useBytes = TRUE))[[1]]
    if (length(m) != 7L || length(s) == 0L) return(NULL)
    a <- as.numeric(m[2]); bb <- as.numeric(m[3])
    data.frame(s = decode(paste(body(s), collapse = "")),
               x = as.numeric(m[6]), y = as.numeric(m[7]),
               size = sqrt(a^2 + bb^2), angle = atan2(bb, a) * 180 / pi,
               font = if (length(f) == 2L && f[2] %in% names(font_of)) font_of[[f[2]]] else NA_character_,
               stringsAsFactors = FALSE)
  }))
  num <- "(?<![0-9.-])(-?[0-9.]+)"
  sg <- regmatches(txt, gregexpr(paste0(num, " (-?[0-9.]+) m\\s+(-?[0-9.]+) (-?[0-9.]+) l\\s+S"),
                                 txt, perl = TRUE, useBytes = TRUE))[[1]]
  segments <- matrix(vapply(strsplit(sg, "\\s+"), function(v) as.numeric(v[c(1, 2, 4, 5)]),
                            numeric(4)), ncol = 4L, byrow = TRUE,
                     dimnames = list(NULL, c("x0", "y0", "x1", "y1")))
  # Every circle (pch 21): a moveto, four Bezier curves and a fill-and-stroke
  # (B), drawn at the line width (w) last set before it. Its width and height
  # are its path's; its outline adds one line width.
  circ_re <- paste0("(-?[0-9.]+) (-?[0-9.]+) m\\s+((?:-?[0-9.]+ -?[0-9.]+ -?[0-9.]+ ",
                    "-?[0-9.]+ -?[0-9.]+ -?[0-9.]+ c\\s+){4})B")
  cm <- gregexpr(circ_re, txt, perl = TRUE, useBytes = TRUE)
  cs <- regmatches(txt, cm)[[1]]
  c_pos <- as.integer(cm[[1]])
  wm <- gregexpr("(?<![0-9.])[0-9.]+ w(?=\\s)", txt, perl = TRUE, useBytes = TRUE)
  w_pos <- as.integer(wm[[1]])
  w_val <- as.numeric(sub(" w$", "", regmatches(txt, wm)[[1]]))
  circles <- do.call(rbind, lapply(seq_along(cs), function(i) {
    v <- as.numeric(regmatches(cs[i], gregexpr("-?[0-9.]+", cs[i]))[[1]])
    xs <- v[c(TRUE, FALSE)]; ys <- v[c(FALSE, TRUE)]
    k <- sum(w_pos < c_pos[i])
    c(d = max(xs) - min(xs), h = max(ys) - min(ys), lw = if (k > 0) w_val[k] else NA_real_)
  }))
  list(mediabox = mb_num, fonts = fonts,
       font_file = grepl("/FontFile", txt, fixed = TRUE, useBytes = TRUE),
       placed = placed, n_text_blocks = length(bt), segments = segments,
       circles = circles)
}
.mm2pt <- function(v) v * 72 / 25.4

.assert("the two PNGs are identical", .md5(FILES[["png_outputs"]]) == .md5(FILES[["png_docs"]]),
        n = 2L)
.assert("the two PDFs are identical, timestamps aside",
        .pdf_md5(FILES[["pdf_outputs"]]) == .pdf_md5(FILES[["pdf_docs"]]), n = 2L)
PNG <- png::readPNG(FILES[["png_docs"]])
.assert("the PNG is the page's size at 600 dpi",
        abs(dim(PNG)[1] - PAGE[["h"]] / 25.4 * PNG_DPI) < 1 &&
          abs(dim(PNG)[2] - PAGE[["w"]] / 25.4 * PNG_DPI) < 1, n = 1L,
        detail = paste(dim(PNG)[1:2], collapse = " x "))
rm(PNG)

RB <- .pdf_read_back(FILES[["pdf_docs"]])
.assert("the PDF's page is the panel's size, in whole points",
        identical(as.integer(RB$mediabox[3:4]), unname(PAGE_PT)), n = 1L,
        detail = paste(RB$mediabox, collapse = " "))
.assert("the PDF's text is live Helvetica, referenced and not embedded, and every string was found",
        "/BaseFont /Helvetica" %in% RB$fonts &&
          all(RB$fonts %in% c("/BaseFont /Helvetica", "/BaseFont /Helvetica-Oblique",
                              "/BaseFont /Helvetica-Bold", "/BaseFont /Helvetica-BoldOblique",
                              "/BaseFont /Symbol")) &&
          !RB$font_file && nrow(RB$placed) == RB$n_text_blocks,
        n = nrow(RB$placed), detail = paste(sub("/BaseFont /", "", RB$fonts), collapse = ", "))
TX <- RB$placed
TX_WORDS <- TX[nzchar(trimws(TX$s)), ]                  # plotmath's spaces are their own strings
TITLE_LINES <- gsub(SOFT_HY, "<SHY>", strsplit(.x_title("pdf"), "\n", fixed = TRUE)[[1]],
                    fixed = TRUE)
WANT <- unname(c("BCL2L1", "gain", "MCL1", "amplification", rep(c("All tumours", SUBTYPES), 2),
                 lab_num(X_BREAKS), TITLE_LINES, KEY_TEXT))
.assert("the PDF shows exactly the strings it should, and no others",
        identical(sort(TX_WORDS$s), sort(WANT)), n = length(WANT),
        detail = paste(utils::head(c(setdiff(TX_WORDS$s, WANT), setdiff(WANT, TX_WORDS$s)), 3),
                       collapse = " / "))
BIG <- c("BCL2L1", "gain", "MCL1", "amplification", TITLE_LINES)
.assert("the gene headers and the axis title are 6 pt and everything else 5 pt; nothing is rotated",
        all(abs(TX_WORDS$size[TX_WORDS$s %in% BIG] - PT_TITLE) < 0.01) &&
          all(abs(TX_WORDS$size[!TX_WORDS$s %in% BIG] - PT_TXT) < 0.01) &&
          all(abs(TX_WORDS$angle) < 1e-6), n = nrow(TX_WORDS))
.assert("the two gene names are in Helvetica's oblique, and every other word in Helvetica",
        setequal(TX_WORDS$s[TX_WORDS$font %in% "Helvetica-Oblique"], c("BCL2L1", "MCL1")) &&
          sum(TX_WORDS$font %in% "Helvetica-Oblique") == 2L &&
          all(TX_WORDS$font[!TX_WORDS$s %in% c("BCL2L1", "MCL1")] %in% "Helvetica"),
        n = nrow(TX_WORDS), detail = paste(unique(TX_WORDS$font), collapse = ", "))

SG <- RB$segments
H_SG <- SG[abs(SG[, "y1"] - SG[, "y0"]) < 1e-9, , drop = FALSE]
V_SG <- SG[abs(SG[, "x1"] - SG[, "x0"]) < 1e-9, , drop = FALSE]
len_h <- abs(H_SG[, "x1"] - H_SG[, "x0"])
len_v <- abs(V_SG[, "y1"] - V_SG[, "y0"])
AX <- H_SG[abs(len_h - .mm2pt(AXIS_MM)) <= TOL_PT, , drop = FALSE]
NL <- V_SG[abs(len_v - .mm2pt(PANEL_H)) <= TOL_PT, , drop = FALSE]
.assert(paste0("in the PDF the x axis measures 26 mm and the line at 1 spans the ",
               "data area at its middle"),
        nrow(AX) == 1L && nrow(NL) == 1L &&
          abs(NL[1, "x0"] - (min(AX[1, c("x0", "x1")]) + .mm2pt(AXIS_MM / 2))) <= TOL_PT,
        n = 2L)
AX_X0 <- min(AX[1, c("x0", "x1")])
# Each interval, read from the PDF: its ends must sit where its odds ratios
# fall on the 26-mm log2 axis.
CI_PDF <- H_SG[len_h < .mm2pt(AXIS_MM) - 1, , drop = FALSE]
want_x0 <- sort(AX_X0 + .mm2pt(log2(FR$lo / X_LIM[1]) * MM_PER_DOUBLING))
want_x1 <- sort(AX_X0 + .mm2pt(log2(FR$hi / X_LIM[1]) * MM_PER_DOUBLING))
got_x0 <- sort(pmin(CI_PDF[, "x0"], CI_PDF[, "x1"]))
got_x1 <- sort(pmax(CI_PDF[, "x0"], CI_PDF[, "x1"]))
dev_ci <- if (nrow(CI_PDF) == 12L) max(abs(c(got_x0 - want_x0, got_x1 - want_x1))) else Inf
.assert("in the PDF every interval's ends sit where its odds ratios fall on the axis",
        is.finite(dev_ci) && dev_ci <= TOL_PT, n = 24L, max_dev = dev_ci)
y_t <- vapply(TITLE_LINES, function(s) TX_WORDS$y[TX_WORDS$s == s], numeric(1))
dev_t <- abs(unname(y_t[1] - y_t[2]) - .mm2pt(TITLE_GAP))
.assert("in the PDF the axis title's two lines are 2.4 mm apart",
        dev_t <= TOL_PT, n = 1L, max_dev = dev_t)
y_key <- vapply(KEY_TEXT, function(s) TX_WORDS$y[TX_WORDS$s == s], numeric(1))
y_top <- max(V_SG[abs(len_v - .mm2pt(PANEL_H)) <= TOL_PT, c("y0", "y1")])
dev_k <- abs(c(unname(y_key[["fga"]]) - y_top - .mm2pt(KEY_ABOVE),
               unname(y_key[["aneuploidy"]] - y_key[["fga"]]) - .mm2pt(KEY_STEP)))
.assert("in the PDF the key lies above the data area: 1.2 mm, then 2.8 mm, baseline to baseline",
        max(dev_k) <= TOL_PT, n = 2L, max_dev = max(dev_k))
y_sub <- TX_WORDS$y[TX_WORDS$s %in% SUBTYPES]
CIRC <- RB$circles
outer_mm <- if (!is.null(CIRC)) (CIRC[, "d"] + CIRC[, "lw"]) * 25.4 / 72 else Inf
.assert(paste0("in the PDF the 14 circles, 12 marks and the key's 2, are 1.05 mm ",
               "across, outline included (the point-size ruling)"),
        !is.null(CIRC) && nrow(CIRC) == 14L && all(abs(outer_mm - CIRCLE_MM) <= 0.01) &&
          all(round(outer_mm, 2) == 1.05) && all(abs(CIRC[, "d"] - CIRC[, "h"]) <= TOL_PT),
        n = if (is.null(CIRC)) 0L else nrow(CIRC), max_dev = max(abs(outer_mm - CIRCLE_MM)),
        detail = sprintf("%.3f to %.3f mm", min(outer_mm), max(outer_mm)))
.assert("in the PDF the rows run top to bottom: BCL2L1's, then MCL1's",
        all(diff(TX_WORDS$y[TX_WORDS$s %in% c("BCL2L1", "MCL1")]) < 0) &&
          length(y_sub) == 8L, n = 10L)

FIGURE_FILES <- tibble::tibble(
  copy = names(FILES), file = basename(FILES),
  checksum = c(.md5(FILES[["png_outputs"]]), .md5(FILES[["png_docs"]]),
               .pdf_md5(FILES[["pdf_outputs"]]), .pdf_md5(FILES[["pdf_docs"]])),
  checksum_is = c("md5", "md5", "md5, timestamps zeroed", "md5, timestamps zeroed"))

# =============================================================================
# 6. The legend draft, Claude Code's, for the author's text ruling
# =============================================================================
.n <- function(v) format(v, big.mark = ",", trim = TRUE)
NZ <- stats::setNames(FR$n[FR$row %in% SUBTYPES & FR$gene == "BCL2L1"],
                      FR$row[FR$row %in% SUBTYPES & FR$gene == "BCL2L1"])
.assert("each subtype's n is the same for both genes",
        identical(unname(NZ[SUBTYPES]),
                  unname(stats::setNames(FR$n[FR$row %in% SUBTYPES & FR$gene == "MCL1"],
                                         FR$row[FR$row %in% SUBTYPES & FR$gene == "MCL1"])[SUBTYPES])),
        n = 4L)
LEGEND_DRAFT <- paste0(
  "BCL2L1 gain and MCL1 amplification in MYC-amplified tumours of TCGA-BRCA ",
  "(n = ", .n(1043L), "): odds ratios against tumours without MYC amplification, ",
  "from logistic models adjusted for aneuploidy score (filled) or, in all tumours, ",
  "for fraction genome altered (open), with 95% confidence intervals, on a log ",
  "scale. Below each gene, PAM50 subtypes, adjusted for aneuploidy score (Basal, ",
  "n = ", .n(NZ[["Basal"]]), "; Her2, n = ", .n(NZ[["Her2"]]), "; LumA, n = ",
  .n(NZ[["LumA"]]), "; LumB, n = ", .n(NZ[["LumB"]]), "). Copy number from GISTIC ",
  "calls: MYC and MCL1 amplification, +2; BCL2L1 gain, +1 or more. The tests in ",
  "all tumours were pre-registered in the validation study; the subtype estimates ",
  "were declared there as secondary and exploratory.")
.assert("the legend draft is ASCII, as the script is",
        !grepl("[^ -~]", LEGEND_DRAFT), n = 1L)

# =============================================================================
# 7. Saved: the object, two tables and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

SPEC <- list(
  ruling = "docs/2026-11-05_E54_coamplification_panel_ruling.md (8fde82a)",
  source = list(file = "data/from_validation_g2/g2_cnv_cooccurrence.rds", md5 = G2_MD5,
                note = "myc_human_validation d3ac60e, docs/2026-08-28_G2_result.md"),
  genes = GENES, subtypes = SUBTYPES, axis_mm = AXIS_MM,
  x = list(limits = X_LIM, breaks = X_BREAKS, minor = X_MINOR, transform = "log2"),
  pitch_mm = PITCH, panel_height_mm = PANEL_H, page_pt = PAGE_PT,
  type_pt = c(text = PT_TXT, titles = PT_TITLE), colour = COL_MARK,
  marks = list(null = c(NULL_COL, NULL_LW), interval = CI_LW, point = c(PT_SIZE, PT_STROKE),
               circle_mm = CIRCLE_MM,
               axes = NM_LINE),
  key = KEY_TEXT, x_title = X_TITLE, title_gap_mm = TITLE_GAP)
LAYOUT <- list(slots = SLOTS, labels = LABS, key = KEY,
               mm = c(data_left = X_DATA, null_line = X_NULL_MM, top_margin = TOP_MM,
                      per_doubling = MM_PER_DOUBLING),
               natural_mm = NEED, page_mm = PAGE)

OUT <- list(
  frame = FR,
  note_checks = CHK,
  layout = LAYOUT,
  legend_draft = LEGEND_DRAFT,
  assertions = ASSERTIONS,
  figure_files = FIGURE_FILES,
  spec = SPEC,
  run_info = list(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  R = R.version.string,
                  bitmap_type_default = getOption("bitmapType")))
saveRDS(OUT, PATH_E54)
utils::write.csv(FR %>% dplyr::select(gene, rule, row, covariate, n, n_both, or, lo, hi, p, tier),
                 file.path(OUT_TABLES, "E54_odds_ratios.csv"), row.names = FALSE)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E54_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), "run_info"),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-13s %s", k, DIGESTS[[k]]))
message("\n   the legend draft (Claude Code's, for the text ruling):\n   ", LEGEND_DRAFT)
message("\nE54 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figure.
  library(dplyr)
  x <- readRDS(here::here("results", "e54_myc_coamplification_panel.rds"))
  x$frame %>% dplyr::select(gene, row, covariate, n, n_both, or, lo, hi, p) %>% print()
  x$note_checks %>% print()
  x$assertions %>% print(n = Inf)
  cat(x$legend_draft, "\n")
  img <- png::readPNG(here::here("docs", "figures",
                                 "2026-11-05_E54_myc_coamplification_odds_ratios.png"))
  grid::grid.newpage()
  grid::grid.raster(img)
}
