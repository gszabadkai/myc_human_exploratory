# E53_half_width_ed_panels.R
# =============================================================================
# EXTENDED DATA FIG. 6a AND 6b (E48'S ED-a AND ED-b) REDRAWN AT HALF WIDTH,
# AT MOST 88 MM, WITH 6b'S ARROWS REMOVED AND ITS LABELS EDITED.
#
# A DISPLAY SCRIPT. It fits nothing and computes no new quantity. It reads
# E48's verified saved object (digests committed at 1ccf805 before E48's run,
# which reproduced all 33) and draws its stored display frames again.
#
# Specification: the author's ruling of 2026-11-05, verbatim in
# docs/2026-11-05_E53_half_width_ed_panels_ruling.md, committed alone at
# 7325ff1 before this script was written.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNGs are compared byte for
# byte with the dry run.
#
# =============================================================================
# WHAT CHANGES FROM E48, AND WHY
# =============================================================================
# Both panels:
#   - Each is at most 88 mm wide (the ruling), so that two fit the 180-mm
#     Extended Data page (myc_mouse 2ebf646, FIGURE_RULES.md, "Page").
#   - The type follows Fig. 1's figure rules at this size: 5 pt, the floor,
#     for labels, values, ticks and key text; 6 pt for titles. Titles are
#     plain. E48 drew them bold at 6.5 pt, and its other text at 5.5 pt.
#   - Each page is the panel's natural size, rounded up to whole points.
#   - The PDFs come from base pdf(), with live Helvetica for Illustrator.
#     E48 wrote cairo_pdf.
# 6a (E48's ED-a):
#   - The TCGA-BRCA and SCAN-B triangles stay side by side, with 4.6-mm cells
#     (the ruling, question 1).
#   - The three neoadjuvant values move to a row of tiles under the colour
#     key, titled "Neoadjuvant cohorts, reference MYC score vs proliferation".
#   - The bands, their greys and their inks are E48's.
#   - It is drawn on a millimetre canvas, so that every element sits where
#     the script computes it, and the checks read the cells back from the PDF.
#   - Every label is set on its baseline, with its capitals centred on its
#     row, tile or key square.
# 6b (E48's ED-b):
#   - The direction arrows and their words are removed (the ruling).
#   - The row labels read "TCGA: OXPHOS level" and so on (question 3).
#   - The key reads "BCL2L1" and "proliferation adjusted", and the same for
#     "BBC3" and "Difference" (question 4), in three rows. The preview's two
#     closing lines ("open = ...", "filled = ...") glossed the symbols for
#     the author and are not drawn: the legend draft says it.
#   - The x axes start at -0.1 and run to 0.6 and 0.3 (question 2). Ticks are
#     labelled at 0, 0.2 and so on, and both panels share one scale, so the
#     difference reads at the coefficients' scale.
#   - The axis titles keep E48's words, each on two lines 2.4 mm apart.
#   - The marks are E48's: colours, shapes, sizes, offsets, intervals, the
#     zero line and the connectors. Its axes take Fig. 1's 0.53-pt line.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing and reads no expression data.
#   - It changes no text of the paper. Its legend drafts are E48's, with one
#     word changed in 6a's: the neoadjuvant values are now "Below", not
#     "Right".
#   - It leaves no raw file bytes in the workspace: each PDF is read back
#     inside a function that returns only ASCII and numbers (E51's lesson).
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
for (pkg in c("digest", "png", "patchwork", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E53 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE53: Extended Data Fig. 6a and 6b at half width\n", strrep("=", 78))
message("Reads E48's verified object and fits nothing. Specification 7325ff1.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_E48        <- file.path(DIR_RESULTS, "e48_separation_and_fig4g.rds")
PATH_E48_SCRIPT <- here::here("scripts", "E48_separation_and_fig4g.R")
PATH_RULING     <- here::here("docs", "2026-11-05_E53_half_width_ed_panels_ruling.md")
PATH_E53        <- file.path(DIR_RESULTS, "e53_half_width_ed_panels.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-05"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- c(a = "E53_score_correlations_half_width",
              b = "E53_bcl2l1_against_bbc3_half_width")

# =============================================================================
# 0. Constants
# =============================================================================
# The elements of E48's object read here, with the digests recorded in
# docs/2026-11-03_e48_prerun_digests.md, committed at 1ccf805 before the
# author's run. That run reproduced all 33.
E48_DIGESTS <- c(
  displays     = "c8397b9106fbe770476a7b10bb513f14",
  legend_draft = "b9dca928c086b687785da7436aca297e",
  assertions   = "f6a6ee8b7d1edcbfad1a277df3b6ef56")

W_MAX_MM <- 88                  # half width (the ruling): two fit the 180-mm page
PT_EPS   <- 0.001               # pdf() truncates the page box to whole points
PT_TXT   <- 5                   # labels, values, ticks and key text
PT_TITLE <- 6                   # titles
CAP_EM   <- 0.718               # Helvetica's capitals above the baseline, in em (_style.R)
PNG_DPI  <- 600
SOFT_HY  <- intToUtf8(173L)     # the PDF's hyphen: its encoding prints "-" as minus
TOL_PT   <- 0.02                # the PDF writes coordinates to 0.01 pt
TOL_SIZE <- 0.01                # and text matrices to 0.01, so a rotated size reads 4.997

# --- 6a: E48's bands, greys and inks (E48, lines 1442 to 1446) --------------------
BANDS     <- c("below 0.60", "0.60 to 0.80", "0.80 or above")
BAND_FILL <- c("below 0.60" = "#F0F0F0", "0.60 to 0.80" = "#BDBDBD",
               "0.80 or above" = "#636363")
BAND_INK  <- c("below 0.60" = "black", "0.60 to 0.80" = "black",
               "0.80 or above" = "white")
BAND_SEP  <- 0.60               # E39's bands, as E48 applies them
BAND_ID   <- 0.80
KEY_TITLE   <- "Spearman rho, absolute"                                    # E48's
STRIP_TITLE <- "Neoadjuvant cohorts, reference MYC score vs proliferation"  # the ruling
TITLE_A <- c(TCGA = "TCGA-BRCA (n = 938)", `SCAN-B` = "SCAN-B (n = 3,143)")  # E48's
# --- 6a: the layout, in mm ----------------------------------------------------------
CELL       <- 4.6               # the cells (the ruling, question 1)
TILE_LW    <- 0.35              # the white line between cells (E48: 0.5, at full width)
M_EDGE     <- 1                 # the page's margin
G_ROWLAB   <- 1                 # a row label to its triangle
G_TRI      <- 3                 # between the two triangles
G_TITLE    <- 1.2               # a triangle's title baseline to its top
G_COLLAB   <- 1                 # a triangle's foot to its column labels
COL_ROT    <- 40                # the column labels' angle, as E48's
G_KEY      <- 2.5               # the lowest column label to the colour key
KEY_SQ     <- 3                 # the colour key's squares
G_KEY_T    <- 2                 # the key's title to its first square
G_SQ_TXT   <- 1                 # a square to its label
G_KEY_IT   <- 3                 # a key label to the next square
G_STRIP    <- 2.5               # the colour key to the neoadjuvant title
G_STRIP_T  <- 1.2               # the neoadjuvant title's baseline to its tiles
G_GSE      <- 1                 # a cohort's name to its tile
G_STRIP_IT <- 4                 # a tile to the next cohort's name
LAB_CLEAR  <- 0.3               # mm kept between two rotated column labels

# --- 6b: E48's marks, with the ruling's texts -----------------------------------------
COL_OX <- "#009E73"             # E48, line 159
KEY_B <- tibble::tribble(
  ~key,             ~text,                    ~col,    ~fill,   ~shape,
  "bcl2l1_without", "BCL2L1",                 COL_OX,  "white", 21,
  "bcl2l1_with",    "proliferation adjusted", COL_OX,  COL_OX,  21,
  "bbc3_without",   "BBC3",                   COL_OX,  "white", 23,
  "bbc3_with",      "proliferation adjusted", COL_OX,  COL_OX,  23,
  "diff_without",   "Difference",             "black", "white", 21,
  "diff_with",      "proliferation adjusted", "black", "black", 21)   # the ruling, question 4
BLOCK_ORDER <- c("TCGA gsva", "TCGA mitopps", "SCAN-B gsva", "SCAN-B mitopps")  # E48's
ROW_LAB_B <- c(`TCGA gsva` = "TCGA: OXPHOS level", `TCGA mitopps` = "TCGA: OXPHOS priority",
               `SCAN-B gsva` = "SCAN-B: OXPHOS level",
               `SCAN-B mitopps` = "SCAN-B: OXPHOS priority")         # the ruling, question 3
X_TITLE_B <- c(coef = "OXPHOS association, MYC activity\nheld fixed (SD per SD)",
               diff = "BCL2L1 minus BBC3\n(bootstrap 95% CI)")          # E48's words
X_LIM_B  <- list(coef = c(-0.1, 0.6), diff = c(-0.1, 0.3))            # the ruling, question 2
X_STEP_B <- 0.2
MM_PER_X <- 54                  # mm per unit, on both x axes: one scale
PANEL_H_B <- 35                 # mm: the panels' height, 7 mm per unit of E48's rows
Y_PAD_B  <- 0.4                 # E48's row units kept above the top mark and below the lowest
ZERO_COL <- "grey45"; ZERO_LW <- 0.3      # E48's zero line (line 938)
CI_LW    <- 0.4                           # E48's intervals (line 964)
CONN_LW  <- 0.3; CONN_ALPHA <- 0.6        # E48's connectors (line 948)
PT_SIZE  <- 1.5; PT_STROKE <- 0.45        # E48's points (lines 967 and 1568)
TITLE_GAP <- 2.4                # mm, baseline to baseline in a two-line axis title (as E52)
# R spaces lines of text 1.2 x the type size apart at lineheight 1, so the
# lineheight that gives the gap is 0.945 (E52's finding).
TITLE_LH <- TITLE_GAP / (1.2 * PT_TITLE * 25.4 / 72)

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

.assert("the size is the ruling's, half width at 88 mm, as recorded in the ruling's own file",
        file.exists(PATH_RULING) && W_MAX_MM == 88 &&
          any(grepl("half width is 88 mm", readLines(PATH_RULING, warn = FALSE),
                    fixed = TRUE)), n = 1L)
# E48's script is read as text, to show that the colours, shapes and marks
# copied above are the ones it drew with.
E48_SRC <- readLines(PATH_E48_SCRIPT, warn = FALSE)
E48_LINES <- c(
  'COL_OX  <- "#009E73"',
  'BAND_FILL <- c("below 0.60" = "#F0F0F0", "0.60 to 0.80" = "#BDBDBD",',
  '               "0.80 or above" = "#636363")',
  'BAND_INK  <- c("below 0.60" = "black", "0.60 to 0.80" = "black",',
  '               "0.80 or above" = "white")',
  '  "bcl2l1_without", "BCL2L1, without the proliferation term",        COL_OX,   "white",  21,',
  '  "bcl2l1_with",    "BCL2L1, with the proliferation term",           COL_OX,   COL_OX,   21,',
  '  "bbc3_without",   "BBC3, without the proliferation term",          COL_OX,   "white",  23,',
  '  "bbc3_with",      "BBC3, with the proliferation term",             COL_OX,   COL_OX,   23,',
  '  "diff_without",   "Difference, without the proliferation term",    "black",  "white",  21,',
  '  "diff_with",      "Difference, with the proliferation term",       "black",  "black",  21,',
  'BLOCK_ORDER <- c("TCGA gsva", "TCGA mitopps", "SCAN-B gsva", "SCAN-B mitopps")',
  '    ggplot2::geom_vline(xintercept = null_x, colour = "grey45", linewidth = 0.3)',
  '        colour = unname(CONN_COL[tm]), linewidth = 0.3, alpha = 0.6,',
  '                          linewidth = 0.4, show.legend = FALSE) +',
  '                        size = size, stroke = 0.45)',
  '                   names(OFF_EDB), breaks = .sym_breaks(XL_EDB[2]), size = 1.5,',
  '    ggplot2::scale_fill_manual(values = BAND_FILL, breaks = BANDS, drop = FALSE,',
  '                               name = "Spearman rho, absolute") +',
  '      axis.text.x = ggplot2::element_text(size = MIN_PT, angle = 40,')
.assert(paste0("the colours, shapes, marks and key title copied here are E48's, as ",
               "its script writes them"),
        all(E48_LINES %in% E48_SRC), n = length(E48_LINES),
        detail = paste(utils::head(setdiff(E48_LINES, E48_SRC), 2), collapse = " / "))

# =============================================================================
# 1. E48's object, checked against its committed digests
# =============================================================================
message("\n1. E48's saved object")
.stop_if(file.exists(PATH_E48), "E48's saved object is missing: source ",
         "scripts/E48_separation_and_fig4g.R first. STOP.")
e48 <- readRDS(PATH_E48)
GOT48 <- vapply(names(E48_DIGESTS), function(k) digest::digest(e48[[k]]),
                character(1))
.assert(paste0("E48's object is the verified one: every element read here ",
               "matches its digest committed before E48's run (1ccf805)"),
        identical(unname(GOT48), unname(E48_DIGESTS)), n = length(GOT48),
        detail = paste(names(GOT48)[GOT48 != E48_DIGESTS], collapse = ", "))
.assert("E48's 175 assertions all passed",
        nrow(e48$assertions) == 175L && all(e48$assertions$pass),
        n = nrow(e48$assertions))

# =============================================================================
# 2. The stored display frames
# =============================================================================
message("\n2. The stored display frames")
TILE   <- as.data.frame(e48$displays$ed_a)
TILE39 <- as.data.frame(e48$displays$ed_a_e39)
PF_B   <- as.data.frame(e48$displays$ed_b_coefficients)
PF_BD  <- as.data.frame(e48$displays$ed_b_difference)

.band <- function(r) ifelse(abs(r) >= BAND_ID, "0.80 or above",
                            ifelse(abs(r) >= BAND_SEP, "0.60 to 0.80", "below 0.60"))
VAR_ORDER <- levels(TILE$x)
COLS_A <- VAR_ORDER[1:6]
ROWS_A <- VAR_ORDER[2:7]
GSE    <- rev(levels(TILE39$y))
.assert(paste0("6a: 21 pairs in each cohort, each below the diagonal of E48's ",
               "seven scores"),
        nrow(TILE) == 42L && all(table(TILE$cohort)[c("TCGA", "SCAN-B")] == 21L) &&
          length(VAR_ORDER) == 7L && identical(levels(TILE$y), rev(VAR_ORDER)) &&
          all(match(as.character(TILE$x), VAR_ORDER) <
                match(as.character(TILE$y), VAR_ORDER)), n = 42L)
.assert(paste0("6a: every value is its rho to two decimals, in its band, and the ",
               "three neoadjuvant cohorts are E48's, in its order"),
        identical(TILE$label, sprintf("%.2f", TILE$rho)) &&
          identical(as.character(TILE$band), .band(TILE$rho)) &&
          identical(TILE39$label, sprintf("%.2f", TILE39$rho)) &&
          identical(as.character(TILE39$band), .band(TILE39$rho)) &&
          identical(GSE, c("GSE25066", "GSE194040", "GSE164458")) &&
          setequal(TILE39$cohort, GSE), n = 45L)

ROW_Y_B <- stats::setNames(-1.2 * seq_along(BLOCK_ORDER), BLOCK_ORDER)  # E48's rows
.assert(paste0("6b: 16 coefficients and 8 differences in E48's four rows, keyed ",
               "as E48 keyed them"),
        nrow(PF_B) == 16L && nrow(PF_BD) == 8L &&
          setequal(PF_B$block, BLOCK_ORDER) && setequal(PF_BD$block, BLOCK_ORDER) &&
          all(abs(PF_B$y_row - ROW_Y_B[PF_B$block]) < 1e-12) &&
          all(abs(PF_BD$y_row - ROW_Y_B[PF_BD$block]) < 1e-12) &&
          setequal(as.character(PF_B$key), KEY_B$key[1:4]) &&
          setequal(as.character(PF_BD$key), KEY_B$key[5:6]) &&
          !any(PF_B$void) && !any(PF_BD$void), n = 24L)
.assert("6b: every interval contains its estimate",
        all(PF_B$lo < PF_B$est & PF_B$est < PF_B$hi) &&
          all(PF_BD$lo < PF_BD$est & PF_BD$est < PF_BD$hi), n = 24L)
.assert("6b: the ruling's axes hold every interval: -0.1 to 0.6, and -0.1 to 0.3",
        min(PF_B$lo) > X_LIM_B$coef[1] && max(PF_B$hi) < X_LIM_B$coef[2] &&
          min(PF_BD$lo) > X_LIM_B$diff[1] && max(PF_BD$hi) < X_LIM_B$diff[2], n = 24L)

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

.nm_axis <- function(step, top, bottom, minor) {
  stopifnot(is.numeric(step), length(step) == 1L, step > 0, top > bottom)
  unit <- if (minor) step / 2 else step
  # Everything is rounded to 10 decimals: in floating point -0.2 + 3 * 0.2 is
  # 0.4000000000000001, a hair past an axis that ends at 0.4, and ggplot then drops
  # that end tick as out of range (caught by check_panel on Fig. 1f, 2026-10-11).
  r  <- function(x) round(x, 10)
  lo <- r(floor(bottom / unit + 1e-9) * unit)
  hi <- r(ceiling(top / unit - 1e-9) * unit)
  list(limits = c(lo, hi),
       breaks = r(seq(ceiling(lo / step - 1e-9) * step, hi + 1e-9 * step, by = step)),
       minor_breaks = if (minor) r(seq(lo, hi + 1e-9 * unit, by = unit)) else NULL,
       guide = ggplot2::guide_axis(minor.ticks = minor))
}

scale_x_nm <- function(step, top, bottom = 0, minor = FALSE, labels = lab_num, ...) {
  a <- .nm_axis(step, top, bottom, minor)
  ggplot2::scale_x_continuous(limits = a$limits, breaks = a$breaks, minor_breaks = a$minor_breaks,
                              guide = a$guide, expand = c(0, 0), labels = labels, ...)
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
# 4. Fig. 6a: the score correlations, on a millimetre canvas
# =============================================================================
message("\n4. Fig. 6a")
# The printed size of plain text at size_pt, in mm, with the PDF device's
# Helvetica metrics: width, height above the baseline, and how far its
# descenders reach below it.
.txt_mm <- function(labels, size_pt) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  t(vapply(labels, function(s) {
    g <- grid::textGrob(s, gp = grid::gpar(fontsize = size_pt, fontfamily = NM_FONT))
    c(w = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
      h = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE),
      d = grid::convertHeight(grid::grobDescent(g), "mm", valueOnly = TRUE))
  }, numeric(3)))
}
CAP_MM <- function(size_pt) CAP_EM * size_pt * 25.4 / 72
.pdf_lab <- function(s) gsub("-", SOFT_HY, s, fixed = TRUE)

TTL_PDF <- .pdf_lab(TITLE_A)
M_ROW  <- .txt_mm(ROWS_A, PT_TXT)
M_COL  <- .txt_mm(COLS_A, PT_TXT)
M_TTL  <- .txt_mm(TTL_PDF, PT_TITLE)
M_BAND <- .txt_mm(BANDS, PT_TXT)
M_KT   <- .txt_mm(KEY_TITLE, PT_TITLE)
M_ST   <- .txt_mm(STRIP_TITLE, PT_TITLE)
M_GSE  <- .txt_mm(GSE, PT_TXT)
M_VAL  <- .txt_mm(unique(c(TILE$label, TILE39$label)), PT_TXT)
.assert("6a: every value fits inside its cell, within the white lines",
        max(M_VAL[, "w"]) < CELL - TILE_LW * ggplot2::.pt * 25.4 / 72, n = nrow(M_VAL),
        detail = sprintf("widest %.2f mm in a %.1f-mm cell", max(M_VAL[, "w"]), CELL))

TH <- COL_ROT * pi / 180
# Across: the row labels, right-aligned, then the two triangles.
X_RL  <- M_EDGE + max(M_ROW[, "w"])
X_TRI <- c(TCGA = X_RL + G_ROWLAB)
X_TRI[["SCAN-B"]] <- X_TRI[["TCGA"]] + length(COLS_A) * CELL + G_TRI
COL_X <- lapply(X_TRI, function(x0) x0 + (seq_along(COLS_A) - 0.5) * CELL)
# Down, as distances from the page's top.
O_TITLE  <- M_EDGE + CAP_MM(PT_TITLE)                  # the titles' baseline
O_MTOP   <- O_TITLE + G_TITLE                          # the triangles' top
O_MBOT   <- O_MTOP + length(ROWS_A) * CELL             # their foot
O_ANCHOR <- O_MBOT + G_COLLAB                          # where each column label ends
# A label rotated by COL_ROT ends at its anchor, with its capitals centred on
# it (vjust = 0.5). Its descenders hang d further below its baseline. So it
# reaches w sin + (h / 2 + d) cos below the anchor, w cos + (h / 2) sin to
# its left, and (h / 2 + d) sin to its right.
DROP_COL <- M_COL[, "w"] * sin(TH) + (M_COL[, "h"] / 2 + M_COL[, "d"]) * cos(TH)
O_LOW    <- O_ANCHOR + max(DROP_COL)
O_KEY    <- O_LOW + G_KEY + KEY_SQ / 2                 # the colour key's centre line
O_STITLE <- O_KEY + KEY_SQ / 2 + G_STRIP + CAP_MM(PT_TITLE)
O_STILE  <- O_STITLE + G_STRIP_T                       # the neoadjuvant tiles' top
H_NEED_A <- O_STILE + CELL + M_EDGE

# The colour key and the neoadjuvant row, each from the left margin.
KX_SQ <- numeric(length(BANDS)); KX_LAB <- numeric(length(BANDS))
x <- M_EDGE + M_KT[1, "w"] + G_KEY_T
for (i in seq_along(BANDS)) {
  KX_SQ[i]  <- x
  KX_LAB[i] <- x + KEY_SQ + G_SQ_TXT
  x <- KX_LAB[i] + M_BAND[i, "w"] + G_KEY_IT
}
SX_LAB <- numeric(length(GSE)); SX_TILE <- numeric(length(GSE))
x <- M_EDGE
for (i in seq_along(GSE)) {
  SX_LAB[i]  <- x + M_GSE[i, "w"]                      # the name's right edge
  SX_TILE[i] <- SX_LAB[i] + G_GSE
  x <- SX_TILE[i] + CELL + G_STRIP_IT
}

ANCHOR_X <- unlist(COL_X, use.names = FALSE)
LAB_W <- rep(M_COL[, "w"], length(COL_X))
LAB_H <- rep(M_COL[, "h"], length(COL_X))
LAB_D <- rep(M_COL[, "d"], length(COL_X))
RIGHT_A <- max(X_TRI[["SCAN-B"]] + length(COLS_A) * CELL,
               ANCHOR_X + (LAB_H / 2 + LAB_D) * sin(TH),
               X_TRI + length(COLS_A) * CELL / 2 + M_TTL[, "w"] / 2,
               KX_LAB[length(BANDS)] + M_BAND[length(BANDS), "w"],
               SX_TILE[length(GSE)] + CELL, M_EDGE + M_ST[1, "w"])
LEFT_A <- min(ANCHOR_X - LAB_W * cos(TH) - LAB_H / 2 * sin(TH),
              X_TRI + length(COLS_A) * CELL / 2 - M_TTL[, "w"] / 2)
W_NEED_A <- RIGHT_A + M_EDGE
.assert("6a: nothing reaches into the left margin", LEFT_A >= M_EDGE - 1e-9, n = 1L,
        detail = sprintf("leftmost %.2f mm", LEFT_A))
# Two neighbouring rotated labels are parallel, and the right one lies on the
# left one's descender side. They cannot touch if their centre lines are
# further apart than the left one's lower half with its descenders and the
# right one's upper half.
SEP <- diff(ANCHOR_X) * sin(TH)
NEED_SEP <- utils::head(LAB_H / 2 + LAB_D, -1) + utils::tail(LAB_H / 2, -1) + LAB_CLEAR
.assert("6a: no column label touches its neighbour",
        all(SEP >= NEED_SEP), n = length(SEP),
        detail = sprintf("closest %.2f mm apart, %.2f needed", min(SEP), max(NEED_SEP)))

PAGE_PT_A <- as.integer(ceiling(c(w = W_NEED_A, h = H_NEED_A) * 72 / 25.4 - 1e-9))
names(PAGE_PT_A) <- c("w", "h")
PAGE_A <- (PAGE_PT_A + PT_EPS) * 25.4 / 72
.assert("6a: half width: the page is at most 88 mm wide",
        PAGE_PT_A[["w"]] * 25.4 / 72 <= W_MAX_MM, n = 1L,
        detail = sprintf("%d pt, %.2f mm", PAGE_PT_A[["w"]], PAGE_PT_A[["w"]] * 25.4 / 72))
.Y <- function(o) PAGE_A[["h"]] - o                     # a distance from the top, as y

# The cells.
ROW_YC <- stats::setNames(.Y(O_MTOP + (seq_along(ROWS_A) - 0.5) * CELL), ROWS_A)
CELLS_A <- dplyr::bind_rows(
  tibble::tibble(part = "triangle", cohort = TILE$cohort,
                 row = as.character(TILE$y), col = as.character(TILE$x),
                 rho = TILE$rho, label = TILE$label, band = as.character(TILE$band),
                 xc = vapply(seq_len(nrow(TILE)), function(i)
                   COL_X[[TILE$cohort[i]]][match(as.character(TILE$x[i]), COLS_A)],
                   numeric(1)),
                 yc = unname(ROW_YC[as.character(TILE$y)])),
  tibble::tibble(part = "neoadjuvant", cohort = TILE39$cohort, row = NA_character_,
                 col = NA_character_, rho = TILE39$rho, label = TILE39$label,
                 band = as.character(TILE39$band),
                 xc = SX_TILE[match(TILE39$cohort, GSE)] + CELL / 2,
                 yc = .Y(O_STILE + CELL / 2)))
KEYSQ <- tibble::tibble(band = BANDS, xc = KX_SQ + KEY_SQ / 2, yc = .Y(O_KEY))

# The text. For a plain string, vjust = 0 sets its baseline at y. (For a
# plotmath expression it sets its lowest point there, as E52 found; there is no
# expression here.)
.on_base <- function(id, lab_pdf, lab_png, size_pt, x, y_base, hjust) {
  tibble::tibble(id = id, label_pdf = lab_pdf, label_png = lab_png, size_pt = size_pt,
                 x = x, y = y_base, y_base = y_base, hjust = hjust, vjust = 0, angle = 0)
}
TEXT_A <- dplyr::bind_rows(
  .on_base("title", TTL_PDF, unname(TITLE_A), PT_TITLE,
           X_TRI + length(COLS_A) * CELL / 2, .Y(O_TITLE), 0.5),
  .on_base("row", ROWS_A, ROWS_A, PT_TXT, X_RL, unname(ROW_YC) - CAP_MM(PT_TXT) / 2, 1),
  .on_base("key_title", KEY_TITLE, KEY_TITLE, PT_TITLE, M_EDGE,
           .Y(O_KEY) - CAP_MM(PT_TITLE) / 2, 0),
  .on_base("key", BANDS, BANDS, PT_TXT, KX_LAB, .Y(O_KEY) - CAP_MM(PT_TXT) / 2, 0),
  .on_base("strip_title", STRIP_TITLE, STRIP_TITLE, PT_TITLE, M_EDGE, .Y(O_STITLE), 0),
  .on_base("cohort", GSE, GSE, PT_TXT, SX_LAB, .Y(O_STILE + CELL / 2) - CAP_MM(PT_TXT) / 2, 1),
  tibble::tibble(id = "column", label_pdf = rep(COLS_A, length(COL_X)),
                 label_png = rep(COLS_A, length(COL_X)), size_pt = PT_TXT,
                 x = ANCHOR_X, y = .Y(O_ANCHOR), y_base = NA_real_,
                 hjust = 1, vjust = 0.5, angle = COL_ROT)) %>%
  dplyr::mutate(size_gg = pt(size_pt))

.panel_a <- function(variant) {
  tx <- TEXT_A
  tx$label <- tx[[paste0("label_", variant)]]
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = CELLS_A, ggplot2::aes(x = xc, y = yc, fill = band),
                       width = CELL, height = CELL, colour = "white", linewidth = TILE_LW) +
    ggplot2::geom_tile(data = KEYSQ, ggplot2::aes(x = xc, y = yc, fill = band),
                       width = KEY_SQ, height = KEY_SQ) +
    ggplot2::geom_text(data = CELLS_A, ggplot2::aes(x = xc, y = yc, label = label,
                                                    colour = band),
                       size = pt(PT_TXT), family = NM_FONT) +
    ggplot2::geom_text(data = tx, ggplot2::aes(x = x, y = y, label = label, hjust = hjust,
                                               vjust = vjust, angle = angle,
                                               size = size_gg),
                       family = NM_FONT, colour = "black") +
    ggplot2::scale_fill_manual(values = BAND_FILL, guide = "none") +
    ggplot2::scale_colour_manual(values = BAND_INK, guide = "none") +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(limits = c(0, PAGE_A[["w"]]), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, PAGE_A[["h"]]), expand = c(0, 0)) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0), legend.position = "none")
}
PANEL_A <- list(pdf = .panel_a("pdf"), png = .panel_a("png"))

BUILD_A <- ggplot2::ggplot_build(PANEL_A$pdf)
vals <- BUILD_A$data[[3]]
key_a <- function(x, y, l) paste(sprintf("%.4f", x), sprintf("%.4f", y), l)
.assert(paste0("6a: every cell drawn carries its stored value, in its place, with ",
               "its band's fill"),
        nrow(BUILD_A$data[[1]]) == 45L && nrow(vals) == 45L &&
          setequal(key_a(vals$x, vals$y, vals$label),
                   key_a(CELLS_A$xc, CELLS_A$yc, CELLS_A$label)) &&
          identical(toupper(BUILD_A$data[[1]]$fill), toupper(unname(BAND_FILL[CELLS_A$band]))),
        n = 45L)
message(sprintf("   page %d x %d pt (%.1f x %.1f mm)", PAGE_PT_A[["w"]], PAGE_PT_A[["h"]],
                PAGE_PT_A[["w"]] * 25.4 / 72, PAGE_PT_A[["h"]] * 25.4 / 72))

# =============================================================================
# 5. Fig. 6b: BCL2L1 against BBC3, two panels on one scale
# =============================================================================
message("\n5. Fig. 6b")
K_COL   <- stats::setNames(KEY_B$col, KEY_B$key)
K_FILL  <- stats::setNames(KEY_B$fill, KEY_B$key)
K_SHAPE <- stats::setNames(KEY_B$shape, KEY_B$key)
K_TEXT  <- stats::setNames(KEY_B$text, KEY_B$key)
FR_B  <- PF_B %>% dplyr::mutate(key = factor(as.character(key), levels = KEY_B$key))
FR_BD <- PF_BD %>% dplyr::mutate(key = factor(as.character(key), levels = KEY_B$key))
YL_B <- c(min(FR_B$y_mark, FR_BD$y_mark) - Y_PAD_B, max(FR_B$y_mark, FR_BD$y_mark) + Y_PAD_B)

# E48's connectors: in each row, a gene's estimate without the proliferation
# term joined to its estimate with it.
.conn <- function(fr, pairs) {
  dplyr::bind_rows(lapply(pairs, function(pr) {
    a <- fr[as.character(fr$key) == pr[1], ]
    b <- fr[as.character(fr$key) == pr[2], ]
    m <- match(a$block, b$block)
    .stop_if(nrow(a) == 4L && nrow(b) == 4L && !anyNA(m), "a connector cannot be paired")
    tibble::tibble(pair = pr[1], block = a$block, x = a$est, xend = b$est[m],
                   y = a$y_mark, yend = b$y_mark[m])
  }))
}
CN_B  <- .conn(FR_B, list(c("bcl2l1_without", "bcl2l1_with"), c("bbc3_without", "bbc3_with")))
CN_BD <- .conn(FR_BD, list(c("diff_without", "diff_with")))

.forest_b <- function(fr, cn, conn_col, xlim, x_title, y_labels, legend) {
  p <- ggplot2::ggplot(fr) +
    ggplot2::geom_vline(xintercept = 0, colour = ZERO_COL, linewidth = ZERO_LW) +
    ggplot2::geom_segment(data = cn, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                          colour = conn_col, linewidth = CONN_LW, alpha = CONN_ALPHA,
                          inherit.aes = FALSE) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = y_mark, yend = y_mark,
                                       colour = key),
                          linewidth = CI_LW, show.legend = FALSE) +
    # show.legend = TRUE: since ggplot2 3.5 a layer draws key symbols only for
    # the levels in its own data unless told to, and the key, drawn from the
    # left panel, also shows the difference's two symbols.
    ggplot2::geom_point(ggplot2::aes(x = est, y = y_mark, colour = key, fill = key,
                                     shape = key),
                        size = PT_SIZE, stroke = PT_STROKE, show.legend = legend) +
    scale_x_nm(step = X_STEP_B, top = xlim[2], bottom = xlim[1], minor = TRUE) +
    ggplot2::scale_y_continuous(limits = YL_B, breaks = unname(ROW_Y_B),
                                labels = if (is.null(y_labels)) NULL else y_labels,
                                expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = K_COL, breaks = KEY_B$key, labels = K_TEXT,
                                 drop = FALSE, name = NULL) +
    ggplot2::scale_fill_manual(values = K_FILL, breaks = KEY_B$key, labels = K_TEXT,
                               drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = K_SHAPE, breaks = KEY_B$key, labels = K_TEXT,
                                drop = FALSE, name = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = x_title, y = NULL) +
    theme_nm(legend = "bottom", txt = PT_TXT, txt_title = PT_TITLE) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(lineheight = TITLE_LH),
                   legend.key.spacing.x = grid::unit(3, "mm"),
                   legend.justification = "left")
  if (is.null(y_labels)) {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                            plot.margin = ggplot2::margin(1, 1, 1, 1.5, "mm"))
  } else {
    p <- p + ggplot2::theme(plot.margin = ggplot2::margin(1, 2, 1, 1, "mm"))
  }
  if (legend) {
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE),
      fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
  } else {
    p <- p + ggplot2::guides(colour = "none", fill = "none", shape = "none")
  }
  p
}
W_PANELS_B <- c(coef = diff(X_LIM_B$coef), diff = diff(X_LIM_B$diff)) * MM_PER_X
.panel_b <- function(variant) {
  rows <- unname(ROW_LAB_B[BLOCK_ORDER])
  if (variant == "pdf") rows <- .pdf_lab(rows)
  left  <- .forest_b(FR_B, CN_B, COL_OX, X_LIM_B$coef, X_TITLE_B[["coef"]], rows, TRUE)
  right <- .forest_b(FR_BD, CN_BD, "black", X_LIM_B$diff, X_TITLE_B[["diff"]], NULL, FALSE)
  # patchwork's own outer margin (5.5 pt each side, from its default theme) is
  # set to 0: the two plots' margins already frame the page.
  patchwork::wrap_plots(left, right, ncol = 2) +
    patchwork::plot_layout(widths = grid::unit(W_PANELS_B, "mm"),
                           heights = grid::unit(PANEL_H_B, "mm"), guides = "collect") +
    patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))) &
    ggplot2::theme(legend.position = "bottom", legend.justification = "left")
}
PANEL_B <- list(pdf = .panel_b("pdf"), png = .panel_b("png"))

# What each panel draws: its points and intervals, against the stored frames.
.check_drawn <- function(p, fr, label) {
  b <- ggplot2::ggplot_build(p)
  geo <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  pts <- b$data[[which(geo == "GeomPoint")]]
  ci  <- b$data[[which(geo == "GeomSegment")[2]]]
  k1 <- function(x, y) paste(sprintf("%.10f", x), sprintf("%.10f", y))
  ok <- nrow(pts) == nrow(fr) && nrow(ci) == nrow(fr) &&
    setequal(k1(pts$x, pts$y), k1(fr$est, fr$y_mark)) &&
    setequal(paste(k1(ci$x, ci$y), sprintf("%.10f", ci$xend)),
             paste(k1(fr$lo, fr$y_mark), sprintf("%.10f", fr$hi)))
  .assert(paste0("6b, ", label, ": every point and interval drawn is its stored ",
                 "estimate and interval, in its row"), ok, n = 2L * nrow(fr))
}
.check_drawn(PANEL_B$pdf[[1]], FR_B, "coefficients")
.check_drawn(PANEL_B$pdf[[2]], FR_BD, "difference")
KEY_DRAWN <- as.data.frame(ggplot2::get_guide_data(PANEL_B$pdf[[1]], "colour"))
pt_layer <- PANEL_B$pdf[[1]]$layers[[which(vapply(PANEL_B$pdf[[1]]$layers, function(l)
  class(l$geom)[1], character(1)) == "GeomPoint")]]
.assert(paste0("6b: the key holds all six entries, in the ruling's words, colours ",
               "and shapes, and draws a symbol for each"),
        nrow(KEY_DRAWN) == 6L && identical(as.character(KEY_DRAWN$.label), KEY_B$text) &&
          identical(toupper(KEY_DRAWN$colour), toupper(KEY_B$col)) &&
          isTRUE(pt_layer$show.legend), n = 6L)

.natural <- function(pw) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- patchwork::patchworkGrob(pw)
  c(width  = sum(grid::convertWidth(g$widths, "mm", valueOnly = TRUE)),
    height = sum(grid::convertHeight(g$heights, "mm", valueOnly = TRUE)))
}
NEED_B <- .natural(PANEL_B$pdf)
PAGE_PT_B <- as.integer(ceiling(NEED_B * 72 / 25.4 - 1e-9))
names(PAGE_PT_B) <- c("w", "h")
PAGE_B <- (PAGE_PT_B + PT_EPS) * 25.4 / 72
over_b <- PAGE_PT_B * 25.4 / 72 - NEED_B
.assert("6b: the page is the panel's natural size, rounded up to whole points",
        all(over_b >= -1e-9) && all(over_b < 25.4 / 72), n = 2L,
        detail = sprintf("%.3f and %.3f mm over it", over_b[1], over_b[2]))
.assert("6b: half width: the page is at most 88 mm wide",
        PAGE_PT_B[["w"]] * 25.4 / 72 <= W_MAX_MM, n = 1L,
        detail = sprintf("%d pt, %.2f mm", PAGE_PT_B[["w"]], PAGE_PT_B[["w"]] * 25.4 / 72))
message(sprintf("   page %d x %d pt (%.1f x %.1f mm)", PAGE_PT_B[["w"]], PAGE_PT_B[["h"]],
                PAGE_PT_B[["w"]] * 25.4 / 72, PAGE_PT_B[["h"]] * 25.4 / 72))

# =============================================================================
# 6. The files: PNGs first, then the PDFs (base pdf(), live Helvetica)
# =============================================================================
message("\n6. The files")
FILES <- unlist(lapply(names(FIG_STEM), function(k) c(
  stats::setNames(file.path(OUT_FIGURES, paste0(FIG_STEM[[k]], ".png")), paste0(k, "_png_outputs")),
  stats::setNames(file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM[[k]], ".png")),
                  paste0(k, "_png_docs")),
  stats::setNames(file.path(OUT_FIGURES, paste0(FIG_STEM[[k]], ".pdf")), paste0(k, "_pdf_outputs")),
  stats::setNames(file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM[[k]], ".pdf")),
                  paste0(k, "_pdf_docs")))))
PAGES <- list(a = PAGE_A, b = PAGE_B)
PANELS <- list(a = PANEL_A, b = PANEL_B)
for (k in names(FIG_STEM)) {
  for (f in FILES[paste0(k, c("_png_outputs", "_png_docs"))]) {
    ggplot2::ggsave(f, PANELS[[k]]$png, width = PAGES[[k]][[1]], height = PAGES[[k]][[2]],
                    units = "mm", dpi = PNG_DPI, device = "png", type = "cairo", bg = "white")
  }
}
.pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(filename, width = width, height = height, family = NM_FONT,
                 useDingbats = FALSE, compress = FALSE)
}
for (k in names(FIG_STEM)) {
  for (f in FILES[paste0(k, c("_pdf_outputs", "_pdf_docs"))]) {
    ggplot2::ggsave(f, PANELS[[k]]$pdf, width = PAGES[[k]][[1]], height = PAGES[[k]][[2]],
                    units = "mm", device = .pdf_device)
  }
}

# --- read back what was written -----------------------------------------------------
.md5 <- function(f) unname(tools::md5sum(f))
# A PDF's checksum with its two timestamps zeroed, which is how the dry run's
# and the author's run's PDFs are compared.
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
# place, size and angle, every straight two-point line, and every rectangle.
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
  # pdf() writes every string in its own BT ... ET block, after its text
  # matrix (a b c d x y Tm): (x, y) starts its baseline, sqrt(a^2 + b^2) is its
  # size in points, and atan2(b, a) its angle.
  bt <- regmatches(txt, gregexpr("(?s)(?<![A-Za-z])BT\\s.*?\\sET(?![A-Za-z])", txt,
                                 perl = TRUE, useBytes = TRUE))[[1]]
  num6 <- paste(rep("(-?[0-9.]+)", 6), collapse = " ")
  placed <- do.call(rbind, lapply(bt, function(b) {
    m <- regmatches(b, regexec(paste0(num6, " Tm"), b, useBytes = TRUE))[[1]]
    s <- regmatches(b, gregexpr(str_re, b, perl = TRUE, useBytes = TRUE))[[1]]
    if (length(m) != 7L || length(s) == 0L) return(NULL)
    a <- as.numeric(m[2]); bb <- as.numeric(m[3])
    data.frame(s = decode(paste(body(s), collapse = "")),
               x = as.numeric(m[6]), y = as.numeric(m[7]),
               size = sqrt(a^2 + bb^2), angle = atan2(bb, a) * 180 / pi,
               stringsAsFactors = FALSE)
  }))
  num <- "(?<![0-9.-])(-?[0-9.]+)"
  sg <- regmatches(txt, gregexpr(paste0(num, " (-?[0-9.]+) m\\s+(-?[0-9.]+) (-?[0-9.]+) l\\s+S"),
                                 txt, perl = TRUE, useBytes = TRUE))[[1]]
  segments <- matrix(vapply(strsplit(sg, "\\s+"), function(v) as.numeric(v[c(1, 2, 4, 5)]),
                            numeric(4)), ncol = 4L, byrow = TRUE,
                     dimnames = list(NULL, c("x0", "y0", "x1", "y1")))
  rc <- regmatches(txt, gregexpr(paste0(num, " (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) re"),
                                 txt, perl = TRUE, useBytes = TRUE))[[1]]
  rects <- matrix(vapply(strsplit(rc, "\\s+"), function(v) as.numeric(v[1:4]), numeric(4)),
                  ncol = 4L, byrow = TRUE, dimnames = list(NULL, c("x", "y", "w", "h")))
  list(mediabox = mb_num, fonts = fonts,
       font_file = grepl("/FontFile", txt, fixed = TRUE, useBytes = TRUE),
       placed = placed, n_text_blocks = length(bt), segments = segments, rects = rects)
}
.mm2pt <- function(v) v * 72 / 25.4
LIVE_FONTS <- c("/BaseFont /Helvetica", "/BaseFont /Helvetica-Oblique",
                "/BaseFont /Helvetica-Bold", "/BaseFont /Helvetica-BoldOblique",
                "/BaseFont /Symbol")
RB <- list()
for (k in names(FIG_STEM)) {
  .assert(paste0(k, ": the two PNGs are identical"),
          .md5(FILES[[paste0(k, "_png_outputs")]]) == .md5(FILES[[paste0(k, "_png_docs")]]),
          n = 2L)
  .assert(paste0(k, ": the two PDFs are identical, timestamps aside"),
          .pdf_md5(FILES[[paste0(k, "_pdf_outputs")]]) == .pdf_md5(FILES[[paste0(k, "_pdf_docs")]]),
          n = 2L)
  png_k <- png::readPNG(FILES[[paste0(k, "_png_docs")]])
  .assert(paste0(k, ": the PNG is the page's size at 600 dpi"),
          abs(dim(png_k)[1] - PAGES[[k]][[2]] / 25.4 * PNG_DPI) < 1 &&
            abs(dim(png_k)[2] - PAGES[[k]][[1]] / 25.4 * PNG_DPI) < 1, n = 1L,
          detail = paste(dim(png_k)[1:2], collapse = " x "))
  RB[[k]] <- .pdf_read_back(FILES[[paste0(k, "_pdf_docs")]])
  .assert(paste0(k, ": the PDF's page is the panel's size, in whole points"),
          identical(as.integer(RB[[k]]$mediabox[3:4]),
                    unname(as.integer(list(a = PAGE_PT_A, b = PAGE_PT_B)[[k]]))), n = 1L,
          detail = paste(RB[[k]]$mediabox, collapse = " "))
  .assert(paste0(k, ": the PDF's text is live Helvetica, referenced and not embedded, ",
                 "and every string was found with its place"),
          "/BaseFont /Helvetica" %in% RB[[k]]$fonts && all(RB[[k]]$fonts %in% LIVE_FONTS) &&
            !RB[[k]]$font_file && nrow(RB[[k]]$placed) == RB[[k]]$n_text_blocks,
          n = nrow(RB[[k]]$placed),
          detail = paste(sub("/BaseFont /", "", RB[[k]]$fonts), collapse = ", "))
}
rm(png_k)

# --- 6a, read from its PDF ------------------------------------------------------------
TX_A <- RB$a$placed
.shy <- function(s) gsub(SOFT_HY, "<SHY>", s, fixed = TRUE)
WANT_A <- unname(c(.shy(TTL_PDF), ROWS_A, rep(COLS_A, 2), CELLS_A$label, KEY_TITLE,
                   BANDS, STRIP_TITLE, GSE))
.assert(paste0("6a: the PDF shows exactly the strings it should: the titles, the ",
               "labels, the 45 values, the key and the neoadjuvant row"),
        identical(sort(TX_A$s), sort(WANT_A)), n = length(WANT_A),
        detail = paste(utils::head(c(setdiff(TX_A$s, WANT_A), setdiff(WANT_A, TX_A$s)), 3),
                       collapse = " / "))
BIG_A <- c(.shy(TTL_PDF), KEY_TITLE, STRIP_TITLE)
.assert("6a: the four titles are 6 pt and everything else 5 pt",
        all(abs(TX_A$size[TX_A$s %in% BIG_A] - PT_TITLE) < TOL_SIZE) &&
          sum(TX_A$s %in% BIG_A) == 4L &&
          all(abs(TX_A$size[!TX_A$s %in% BIG_A] - PT_TXT) < TOL_SIZE), n = nrow(TX_A),
        max_dev = max(abs(TX_A$size - ifelse(TX_A$s %in% BIG_A, PT_TITLE, PT_TXT))))
ROT_A <- abs(TX_A$angle) > 1e-6
.assert("6a: the twelve column labels are at 40 degrees, and no other text is rotated",
        sum(ROT_A) == 12L && all(abs(TX_A$angle[ROT_A] - COL_ROT) < 0.05) &&
          all(TX_A$s[ROT_A] %in% COLS_A), n = nrow(TX_A))
RC_A <- RB$a$rects
CELL_PT <- .mm2pt(CELL)
n_cell <- sum(abs(RC_A[, "w"] - CELL_PT) <= TOL_PT & abs(RC_A[, "h"] - CELL_PT) <= TOL_PT)
n_sq   <- sum(abs(RC_A[, "w"] - .mm2pt(KEY_SQ)) <= TOL_PT &
                abs(RC_A[, "h"] - .mm2pt(KEY_SQ)) <= TOL_PT)
.assert(paste0("6a: in the PDF, 45 cells measure 4.6 mm (", sprintf("%.2f", CELL_PT),
               " pt) and the key's three squares 3 mm"),
        n_cell == 45L && n_sq == 3L, n = 48L,
        detail = sprintf("%d cells, %d squares", n_cell, n_sq))
rl <- TX_A[TX_A$s %in% ROWS_A & !ROT_A, ]                 # the column labels are rotated
dev_rl <- max(abs(rl$y - .mm2pt(ROW_YC[rl$s] - CAP_MM(PT_TXT) / 2)))
.assert("6a: in the PDF, each row label's capitals are centred on its row",
        nrow(rl) == 6L && dev_rl <= TOL_PT, n = 6L, max_dev = dev_rl)

# --- 6b, read from its PDF ------------------------------------------------------------
TX_B <- RB$b$placed
ROWS_B_PDF <- .shy(.pdf_lab(unname(ROW_LAB_B[BLOCK_ORDER])))
TICKS_B <- c(lab_num(c(0, 0.2, 0.4, 0.6)), lab_num(c(0, 0.2)))
TITLE_LINES_B <- unlist(strsplit(X_TITLE_B, "\n", fixed = TRUE), use.names = FALSE)
WANT_B <- unname(c(ROWS_B_PDF, TICKS_B, TITLE_LINES_B, KEY_B$text))
.assert(paste0("6b: the PDF shows exactly the strings it should: the four row ",
               "labels, the ticks, the two axis titles and the six key entries, so no ",
               "arrow's words and no old key text"),
        identical(sort(TX_B$s), sort(WANT_B)), n = length(WANT_B),
        detail = paste(utils::head(c(setdiff(TX_B$s, WANT_B), setdiff(WANT_B, TX_B$s)), 3),
                       collapse = " / "))
.assert("6b: the axis titles are 6 pt and everything else 5 pt; nothing is rotated",
        all(abs(TX_B$size[TX_B$s %in% TITLE_LINES_B] - PT_TITLE) < TOL_SIZE) &&
          all(abs(TX_B$size[!TX_B$s %in% TITLE_LINES_B] - PT_TXT) < TOL_SIZE) &&
          all(abs(TX_B$angle) < 1e-6), n = nrow(TX_B),
        max_dev = max(abs(TX_B$size - ifelse(TX_B$s %in% TITLE_LINES_B, PT_TITLE, PT_TXT))))
y_rows <- vapply(ROWS_B_PDF, function(s) TX_B$y[TX_B$s == s], numeric(1))
.assert("6b: the row labels run from TCGA's OXPHOS level at the top to SCAN-B's priority",
        all(diff(y_rows) < 0), n = 4L)
SG_B <- RB$b$segments
H_SG <- SG_B[abs(SG_B[, "y1"] - SG_B[, "y0"]) < 1e-9, , drop = FALSE]
V_SG <- SG_B[abs(SG_B[, "x1"] - SG_B[, "x0"]) < 1e-9, , drop = FALSE]
len_h <- abs(H_SG[, "x1"] - H_SG[, "x0"])
len_v <- abs(V_SG[, "y1"] - V_SG[, "y0"])
AX_PT <- .mm2pt(W_PANELS_B)
ax_l <- H_SG[abs(len_h - AX_PT[["coef"]]) <= TOL_PT, , drop = FALSE]
ax_r <- H_SG[abs(len_h - AX_PT[["diff"]]) <= TOL_PT, , drop = FALSE]
zero <- V_SG[abs(len_v - .mm2pt(PANEL_H_B)) <= TOL_PT, , drop = FALSE]
dev_ax <- if (nrow(ax_l) == 1L && nrow(ax_r) == 1L) {
  max(abs(c(len_h[abs(len_h - AX_PT[["coef"]]) <= TOL_PT] - AX_PT[["coef"]],
            len_h[abs(len_h - AX_PT[["diff"]]) <= TOL_PT] - AX_PT[["diff"]])))
} else Inf
.assert(paste0("6b: in the PDF the two x axes measure ", sprintf("%.1f", W_PANELS_B[["coef"]]),
               " and ", sprintf("%.1f", W_PANELS_B[["diff"]]), " mm, 54 mm per unit, ",
               "so both panels share one scale"),
        is.finite(dev_ax) && dev_ax <= TOL_PT, n = 2L, max_dev = dev_ax)
x0_l <- min(ax_l[, c("x0", "x1")]); x0_r <- min(ax_r[, c("x0", "x1")])
dev_zero <- if (nrow(zero) == 2L) {
  max(abs(sort(zero[, "x0"]) - c(x0_l, x0_r) - .mm2pt(0.1 * MM_PER_X)))
} else Inf
.assert("6b: in the PDF each zero line spans its panel, 0.1 of a unit from the axis's start",
        is.finite(dev_zero) && dev_zero <= TOL_PT, n = 2L, max_dev = dev_zero)
y_title2 <- vapply(TITLE_LINES_B, function(s) TX_B$y[TX_B$s == s], numeric(1))
dev_tl <- max(abs(c(y_title2[1] - y_title2[2], y_title2[3] - y_title2[4]) -
                    .mm2pt(TITLE_GAP)))
.assert("6b: in the PDF each axis title's two lines are 2.4 mm apart",
        dev_tl <= TOL_PT, n = 2L, max_dev = dev_tl)
y_key <- TX_B$y[TX_B$s %in% KEY_B$text]
.assert("6b: in the PDF the key lies below the axis titles, in three rows",
        max(y_key) < min(y_title2) && length(unique(round(y_key, 2))) == 3L,
        n = length(y_key))

FIGURE_FILES <- tibble::tibble(
  copy = names(FILES), file = basename(FILES),
  checksum = vapply(names(FILES), function(nm) {
    if (grepl("_pdf_", nm)) .pdf_md5(FILES[[nm]]) else .md5(FILES[[nm]])
  }, character(1), USE.NAMES = FALSE),
  checksum_is = ifelse(grepl("_pdf_", names(FILES)), "md5, timestamps zeroed", "md5"))

# =============================================================================
# 7. The legend drafts, for the author: E48's, with one word changed in 6a's
# =============================================================================
LEG_A_E48 <- e48$legend_draft$ed_a
LEGEND_DRAFTS <- list(
  fig_6a = sub("Right, the reference", "Below, the reference", LEG_A_E48, fixed = TRUE),
  fig_6b = e48$legend_draft$ed_b)
.assert(paste0("the legend drafts are E48's, except that 6a's neoadjuvant values ",
               "are 'Below', not 'Right'"),
        lengths(gregexpr("Right, the reference", LEG_A_E48, fixed = TRUE)) == 1L &&
          !identical(LEGEND_DRAFTS$fig_6a, LEG_A_E48) &&
          identical(sub("Below, the reference", "Right, the reference",
                        LEGEND_DRAFTS$fig_6a, fixed = TRUE), LEG_A_E48) &&
          identical(LEGEND_DRAFTS$fig_6b, e48$legend_draft$ed_b), n = 2L)
.assert("the legend drafts are ASCII, as the script is",
        !any(grepl("[^ -~]", unlist(LEGEND_DRAFTS))), n = 2L)

# =============================================================================
# 8. Saved: the object, a table and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

SPEC <- list(
  ruling = "docs/2026-11-05_E53_half_width_ed_panels_ruling.md (7325ff1)",
  source = "E48's saved object, its display frames (displays$ed_a, ed_a_e39, ed_b_coefficients, ed_b_difference)",
  width_max_mm = W_MAX_MM, type_pt = c(text = PT_TXT, titles = PT_TITLE),
  fig_6a = list(page_pt = PAGE_PT_A, cell_mm = CELL, tile_line = TILE_LW,
                bands = BAND_FILL, inks = BAND_INK, column_angle = COL_ROT,
                key_title = KEY_TITLE, strip_title = STRIP_TITLE, titles = TITLE_A),
  fig_6b = list(page_pt = PAGE_PT_B, keys = KEY_B, rows = ROW_LAB_B,
                x_titles = X_TITLE_B, x_limits = X_LIM_B, x_step = X_STEP_B,
                mm_per_unit = MM_PER_X, panel_height_mm = PANEL_H_B,
                title_gap_mm = TITLE_GAP, title_lineheight = TITLE_LH,
                marks = list(zero = c(ZERO_COL, ZERO_LW), interval = CI_LW,
                             connector = c(CONN_LW, CONN_ALPHA),
                             point = c(PT_SIZE, PT_STROKE), axes = NM_LINE)),
  source_style = "myc_mouse 2ebf646: figures/natmetab/_style.R; E48's marks and colours")
LAYOUT <- list(
  fig_6a = list(cells = CELLS_A, key_squares = KEYSQ, text = TEXT_A,
                need_mm = c(w = W_NEED_A, h = H_NEED_A), page_mm = PAGE_A,
                label_mm = list(rows = M_ROW, columns = M_COL, values = M_VAL),
                column_label_separation_mm = SEP),
  fig_6b = list(connectors = list(coefficients = CN_B, difference = CN_BD),
                y_limits = YL_B, panel_widths_mm = W_PANELS_B, natural_mm = NEED_B,
                page_mm = PAGE_B))
FRAMES <- list(ed_a = TILE, ed_a_e39 = TILE39, ed_b_coefficients = PF_B,
               ed_b_difference = PF_BD)

OUT <- list(
  frames = FRAMES,
  layout = LAYOUT,
  legend_drafts = LEGEND_DRAFTS,
  assertions = ASSERTIONS,
  figure_files = FIGURE_FILES,
  spec = SPEC,
  run_info = list(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  R = R.version.string,
                  bitmap_type_default = getOption("bitmapType")))
saveRDS(OUT, PATH_E53)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E53_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), "run_info"),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-13s %s", k, DIGESTS[[k]]))
message("\n   the legend drafts:\n   6a: ", LEGEND_DRAFTS$fig_6a, "\n   6b: ",
        LEGEND_DRAFTS$fig_6b)
message("\nE53 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figures, so
  # it also works in a fresh session, not only straight after source().
  library(dplyr)
  x <- readRDS(here::here("results", "e53_half_width_ed_panels.rds"))
  x$assertions %>% print(n = Inf)
  x$figure_files %>% print()
  x$spec$fig_6a$page_pt
  x$spec$fig_6b$page_pt
  cat(x$legend_drafts$fig_6a, "\n\n", x$legend_drafts$fig_6b, "\n")
  # the panels as written: the PNG previews (the PDFs beside them are the files
  # for Illustrator)
  for (stem in c("E53_score_correlations_half_width", "E53_bcl2l1_against_bbc3_half_width")) {
    img <- png::readPNG(here::here("docs", "figures", paste0("2026-11-05_", stem, ".png")))
    grid::grid.newpage()
    grid::grid.raster(img)
  }
}
