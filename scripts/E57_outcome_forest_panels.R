# E57_outcome_forest_panels.R
# =============================================================================
# FIG. 4'S LAST TWO MAIN PANELS: THE OUTCOME FORESTS, AS TWO SEPARATE PANELS.
# PATHOLOGICAL COMPLETE RESPONSE IN THREE NEOADJUVANT COHORTS, AND
# BREAST-CANCER-SPECIFIC SURVIVAL IN METABRIC.
#
# A DISPLAY SCRIPT. It fits nothing. Every value it draws is E45's, read from
# E45's saved object and checked against the digests recorded before E45's run
# (docs/2026-10-22_fig4_human_panels_proposal.md, section 7.6). E45 read them
# from E40, E43 and E44 and exponentiated them; nothing is recomputed here.
#
# Specification: the author's ruling of 2026-11-06, verbatim in
# docs/2026-11-06_E57_outcome_forest_panels_ruling.md, committed alone at
# cb3c7d3 before this script was written.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNGs are compared byte for
# byte with the dry run.
#
# =============================================================================
# WHAT IT DRAWS, AND WHY
# =============================================================================
#   - E45's two forests, each as its own panel and file (the ruling's message).
#   - The colours (the ruling, question 1): OXPHOS in Fig. 4x's dark slate
#     (#243642), open without the proliferation term and filled with it;
#     proliferation in mid slate (#8FA3B0). They are the ends of Fig. 1's
#     unsigned palette and the PUMA/BCL-xL panel's two colours. E45's Okabe-Ito
#     green and orange are the paper's PGC-1a construct and MYC-positive 6-week
#     mice in Figs 1 to 3.
#   - Every mark is a circle; proliferation's differs by colour alone (the
#     message).
#   - No bold anywhere (the message). The PDF read-back checks the font of
#     every string.
#   - The keys (the message and question 4): "OXPHOS", "OXPHOS, adjusted for
#     proliferation" and "Proliferation, adjusted for OXPHOS". METABRIC's panel
#     has no model without proliferation, so its key has the last two.
#   - Fig. 4x's grammar (question 2): a 26-mm axis; 5-pt text with 6-pt
#     titles; its line weights; 1.05-mm circles; the key above the data area,
#     as in Fig. 4x and E52.
#   - The direction cues stay (question 3), at 5 pt.
#   - Otherwise E45's content: the rows and their labels, the column headers,
#     the axes (log scale, E45's limits and breaks), the line from OXPHOS
#     without to with proliferation, the rule before the fourth row, and the
#     contrast column. E45's row geometry is kept, scaled to 5.6 mm per unit.
#
# =============================================================================
# CHOICES THE RULING LEAVES OPEN (each the author's to reverse)
# =============================================================================
#   - Each page is the panel's natural size. Both data areas are 33.0 mm tall
#     and sit the same height above their pages' bottom edges, so the two
#     panels, aligned by their bottom edges, have their axes and rows in line.
#   - The title is kept, plain, at 6 pt, at the page's left margin above the
#     key. E45 centred it over the forest; the key is now where it was.
#   - The tick labels drop their trailing zeros, as Fig. 4x's do ("0.5", not
#     "0.50"): at 26 mm, "1.50" and "2.00" would touch.
#   - The line from OXPHOS without to with proliferation is dark slate at
#     0.43 pt, solid. E45 drew it at 60% opacity, which here would make a grey
#     close to proliferation's mid slate.
#   - The "2" of I-squared is set at the full 5 pt, raised, as E52 sets its
#     exponents: nothing on a panel is smaller than 5 pt. Its top reaches the
#     line above, so the pooled row's two lines are 0.4 mm further apart than
#     the other rows', still centred on the row.
#   - E45's greys stay: the second lines and I-squared in grey35, the cues in
#     grey30, the rule in grey75.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits, refits and pools nothing, and recomputes no value.
#   - It prints no OXPHOS or proliferation ratio: the contrast column is the
#     only text that carries numbers, as in E45 (2026-10-04_handoff.md 2).
#   - It changes no text of the paper. Its legend drafts are E45's, edited for
#     the new keys, for the author's text ruling.
#   - It leaves no raw file bytes in the workspace: each PDF is read back inside
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
    stop("E57 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE57: Fig. 4's outcome forests, as two panels\n", strrep("=", 78))
message("Draws E45's saved values and fits nothing. Specification cb3c7d3.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_E45        <- file.path(DIR_RESULTS, "e45_fig4_human_panels.rds")
PATH_RULING     <- here::here("docs", "2026-11-06_E57_outcome_forest_panels_ruling.md")
PATH_E45_SCRIPT <- here::here("scripts", "E45_fig4_human_panels.R")
PATH_E52_SCRIPT <- here::here("scripts", "E52_puma_bclxl_panel_26mm.R")
PATH_E54_SCRIPT <- here::here("scripts", "E54_myc_coamplification_panel.R")
PATH_E57        <- file.path(DIR_RESULTS, "e57_outcome_forest_panels.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-06"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- c(pcr = "E57_pcr_forest", metabric = "E57_metabric_survival_forest")
HERE_ENV <- environment()

# =============================================================================
# 0. Constants
# =============================================================================
# E45's saved object: the digests recorded before its run (the proposal,
# section 7.6), for the five elements read here.
E45_DIGESTS <- c(plot_frame      = "c571093f601b47b2ac4aae3f27a8c468",
                 contrast_column = "1eecec69b3d3fc480296ee010f515bdc",
                 rows            = "bd0aaba03ffb1ceebd3a7fe6ccbdf6ea",
                 assertions      = "047181c1dc163720a244ee04dae65b01",
                 legend_draft    = "615da6d3195c62b9b1120204a41c8d47")

# --- Fig. 4x's size, type and marks: E54's values, checked in section 1 -------
AXIS_MM  <- 26                  # the x axis (question 2)
PT_EPS   <- 0.001               # pdf() truncates the page box to whole points
PT_TXT   <- 5                   # labels, key, ticks, cues, contrast column
PT_TITLE <- 6                   # the panel's title and the axis title
CAP_EM   <- 0.718               # Helvetica's capitals above the baseline, in em
COL_MARK <- "#243642"           # Fig. 4x's dark slate
NULL_COL <- "grey45"; NULL_LW <- 0.3      # the line at 1
CI_LW    <- 0.4                           # the intervals
PT_SIZE  <- 1.0; PT_STROKE <- 0.3         # 1.05-mm circles, outline included
CIRCLE_MM <- (0.75 * (PT_SIZE * ggplot2::.pt + PT_STROKE * ggplot2::.stroke / 2) +
                0.75 * PT_STROKE * ggplot2::.stroke / 2) * 25.4 / 72
M_EDGE   <- 1                   # mm: the page's margin
G_LAB    <- 1.5                 # mm: the longest row label to the data area
KEY_ABOVE <- 1.2                # mm: the key's lowest baseline above the data area
KEY_STEP <- 2.8                 # mm: the key's rows, baseline to baseline
KEY_GAP  <- 0.8                 # mm: a key symbol to its text
TITLE_GAP <- 2.4                # mm: the axis title's lines, baseline to baseline
TITLE_LH <- TITLE_GAP / (1.2 * PT_TITLE * 25.4 / 72)   # R spaces lines 1.2 x the size
PNG_DPI  <- 600
TOL_PT   <- 0.02                # the PDF writes coordinates to 0.01 pt
SOFT_HY  <- intToUtf8(173L)     # the PDF's hyphen: its encoding prints "-" as minus

# --- the colours (question 1) ---------------------------------------------------
COL_OX <- COL_MARK              # OXPHOS: dark slate
COL_PR <- "#8FA3B0"             # proliferation: mid slate (checked against E52)

# --- E45's content: its values, checked in section 1 ----------------------------
PCR_ORDER <- c("GSE194040", "GSE164458", "GSE25066")     # top to bottom
MB_ORDER <- c("METABRIC / ER-pos (IHC)", "METABRIC / ER-pos (ER.Expr)",
              "METABRIC / pooled", "METABRIC / ER-neg (IHC)")
Y_ROWS  <- c(4.6, 3.6, 2.6, 1.2)    # three rows, a rule, then the fourth
Y_SEP   <- 1.9
Y_CUE   <- 0.28
Y_HEAD  <- 5.45
Y_LIM   <- c(0, 5.9)
KEY_LEVELS <- c("OXPHOS, m0", "OXPHOS, m1", "Proliferation, m1")
XLIM_PCR   <- c(0.5, 3.5)
BREAKS_PCR <- c(0.5, 0.75, 1, 1.5, 2, 3)
XLIM_MB    <- c(0.7, 1.6)
BREAKS_MB  <- c(0.75, 1, 1.25, 1.5)
CUE_IN  <- 1.04                 # a cue's arrow runs from 1.04 (or 1/1.04) ...
CUE_OUT <- 1.10                 # ... to 1.10 inside the axis's end (E45)
TEXT <- list(
  pcr = list(title    = "Pathological complete response",
             head     = "Cohort\nn / pCR",
             con_head = "Proliferation / OXPHOS\nratio of ORs (95% CI)",
             x_title  = "Odds ratio for pCR\nper within-cohort SD",
             cue      = c("lower pCR", "higher pCR")),
  metabric = list(title    = "Breast-cancer-specific survival",
                  head     = "METABRIC\nn / BC deaths",
                  con_head = "Proliferation / OXPHOS\nratio of HRs (95% CI)",
                  x_title  = "Hazard ratio for breast-cancer death\nper within-cohort SD",
                  cue      = c("longer survival", "shorter survival")))
SUB_COL  <- "grey35"            # E45's second lines
CUE_COL  <- "grey30"            # E45's cues
RULE_COL <- "grey75"; RULE_LW <- 0.25     # E45's rule

# --- the keys (the message; question 4) ------------------------------------------
KEY_TEXT <- c("OXPHOS, m0"        = "OXPHOS",
              "OXPHOS, m1"        = "OXPHOS, adjusted for proliferation",
              "Proliferation, m1" = "Proliferation, adjusted for OXPHOS")
KEY_FILL <- c("OXPHOS, m0" = "white", "OXPHOS, m1" = COL_OX, "Proliferation, m1" = COL_PR)
KEY_LINE <- c("OXPHOS, m0" = COL_OX, "OXPHOS, m1" = COL_OX, "Proliferation, m1" = COL_PR)

# --- the geometry at this size ------------------------------------------------------
S_MM      <- 5.6                # mm per unit of E45's row geometry: rows 5.6 mm apart
PANEL_H   <- Y_LIM[2] * S_MM    # the data area's height: 33.04 mm
CON_GAP   <- 2.0                # mm: the data area to the contrast column
LINE_STEP <- 1.2 * PT_TXT * 25.4 / 72     # a two-line label, baseline to baseline (2.12 mm)
I2_STEP   <- LINE_STEP + 0.4    # the pooled row's lines: its raised 2 reaches 2.12 mm
TITLE_ABOVE <- 3.2              # mm: the title's baseline above the key's top baseline
CUE_UP    <- 1.0                # mm: a cue's baseline above its arrow
ARROW_MM  <- 0.8                # mm: the cues' arrowheads (E45: 1.1 mm at its larger type)
CUE_LW    <- 0.25               # the cues' arrows, at the axes' weight (E45: 0.3)
CONN_LW   <- 0.2                # OXPHOS without to with proliferation: 0.43 pt

# --- helpers: E54's, unchanged (checked in section 1) ------------------------------
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

# Fig. 1's style helpers, as E54 copied them from myc_mouse 2ebf646,
# figures/natmetab/_style.R
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

# patchwork fixes the data area at 26 x PANEL_H mm; its own outer margin is set
# to 0, so that the plot's margins frame the page (E53's finding).
.wrap <- function(p) {
  patchwork::wrap_plots(p) +
    patchwork::plot_layout(widths = grid::unit(AXIS_MM, "mm"),
                           heights = grid::unit(PANEL_H, "mm")) +
    patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0)))
}
.natural <- function(pw) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- patchwork::patchworkGrob(pw)
  c(width  = sum(grid::convertWidth(g$widths, "mm", valueOnly = TRUE)),
    height = sum(grid::convertHeight(g$heights, "mm", valueOnly = TRUE)))
}
.pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(filename, width = width, height = height, family = NM_FONT,
                 useDingbats = FALSE, compress = FALSE)
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

# --- helpers new here -----------------------------------------------------------
# The paint at every circle and every two-point line, read from the PDF's
# content in order, as a PDF reader applies it: fill (scn), stroke (SCN) and
# line width (w) are set, saved with q and restored with Q. Inside a function,
# returning only numbers and hex colours (E51's lesson).
.pdf_paint <- function(f) {
  txt <- .pdf_text(f)
  find <- function(re, type) {
    m <- gregexpr(re, txt, perl = TRUE, useBytes = TRUE)[[1]]
    if (m[1] < 0) return(NULL)
    data.frame(pos = as.integer(m), type = type,
               s = regmatches(txt, list(m))[[1]], stringsAsFactors = FALSE)
  }
  n1 <- "-?[0-9.]+"
  circ_re <- paste0("(?<![0-9.-])", n1, " ", n1, " m\\s+(?:", paste(rep(n1, 6), collapse = " "),
                    " c\\s+){4}B")
  seg_re <- paste0("(?<![0-9.-])", n1, " ", n1, " m\\s+", n1, " ", n1, " l\\s+S")
  ev <- rbind(find("(?<=\\s)q(?=\\s)", "q"), find("(?<=\\s)Q(?=\\s)", "Q"),
              find("(?<![0-9.])[0-9.]+ [0-9.]+ [0-9.]+ scn", "fill"),
              find("(?<![0-9.])[0-9.]+ [0-9.]+ [0-9.]+ SCN", "stroke"),
              find("(?<![0-9.])[0-9.]+ w(?=\\s)", "lw"),
              find(circ_re, "circle"), find(seg_re, "segment"))
  ev <- ev[order(ev$pos), ]
  hex <- function(s) {
    v <- as.numeric(strsplit(s, " ", fixed = TRUE)[[1]][1:3])
    grDevices::rgb(v[1], v[2], v[3])
  }
  nums <- function(s) as.numeric(regmatches(s, gregexpr(n1, s))[[1]])
  st <- list(fill = NA_character_, stroke = NA_character_, lw = NA_real_)
  stack <- list(); circ <- list(); segs <- list()
  for (i in seq_len(nrow(ev))) {
    t <- ev$type[i]; s <- ev$s[i]
    if (t == "q") {
      stack[[length(stack) + 1L]] <- st
    } else if (t == "Q") {
      if (length(stack)) { st <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL }
    } else if (t == "fill") {
      st$fill <- hex(s)
    } else if (t == "stroke") {
      st$stroke <- hex(s)
    } else if (t == "lw") {
      st$lw <- as.numeric(sub(" w$", "", s))
    } else if (t == "circle") {
      v <- nums(s); xs <- v[c(TRUE, FALSE)]; ys <- v[c(FALSE, TRUE)]
      circ[[length(circ) + 1L]] <- data.frame(
        cx = (max(xs) + min(xs)) / 2, cy = (max(ys) + min(ys)) / 2,
        d = max(xs) - min(xs), h = max(ys) - min(ys), lw = st$lw,
        fill = st$fill, stroke = st$stroke, stringsAsFactors = FALSE)
    } else if (t == "segment") {
      v <- nums(s)
      segs[[length(segs) + 1L]] <- data.frame(x0 = v[1], y0 = v[2], x1 = v[3], y1 = v[4],
                                              lw = st$lw, stroke = st$stroke,
                                              stringsAsFactors = FALSE)
    }
  }
  list(circles = do.call(rbind, circ), segments = do.call(rbind, segs))
}
# A ggplot2 linewidth, as the PDF writes it: lwd = linewidth x .pt, at 0.75 pt each,
# to the PDF's two decimals.
.w_pt <- function(lw) lw * ggplot2::.pt * 0.75
.lw_is <- function(w, lw) !is.na(w) & abs(w - .w_pt(lw)) < 0.006
.hex <- function(col) toupper(grDevices::rgb(t(grDevices::col2rgb(col)), maxColorValue = 255))
.hy <- function(s) gsub("-", SOFT_HY, s, fixed = TRUE)   # the PDF's hyphen
.P <- function(p) c(pcr = "pCR", metabric = "METABRIC")[[p]]   # a panel's name in a check
.shy <- function(s) gsub("-", "<SHY>", s, fixed = TRUE)  # ... as the read-back decodes it

# =============================================================================
# 1. The ruling, and what is copied from E45, E52 and E54
# =============================================================================
message("\n1. The ruling and the copies")
RULING <- readLines(PATH_RULING, warn = FALSE)
.assert("the panels are the ruling's: its four answers, as its own file records",
        all(vapply(c("> Mid slate (Recommended):", "> 26-mm axis (Recommended):",
                     "> Keep them (Recommended):", "> Adjusted for OXPHOS (Recommended):"),
                   function(s) any(startsWith(RULING, s)), logical(1))), n = 4L)
.assert("the keys are the ruling's words",
        any(grepl("keys again OXPHOS; OXPHOS, adjusted for proliferation", RULING, fixed = TRUE)) &&
          any(grepl("'Proliferation, adjusted for OXPHOS'", RULING, fixed = TRUE)) &&
          identical(unname(KEY_TEXT), c("OXPHOS", "OXPHOS, adjusted for proliferation",
                                        "Proliferation, adjusted for OXPHOS")), n = 3L)

# A script's top-level assignments of the named objects, evaluated in the order
# it makes them, in a fresh environment: what that script's objects are,
# without running it.
.from_script <- function(path, names) {
  ex <- parse(path, keep.source = FALSE)
  env <- new.env(parent = globalenv())
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% names) eval(e, env)
  }
  env
}
# Two objects are the same if their values are, or, for functions, if their
# arguments and bodies are, token for token (comments and layout aside).
.same <- function(a, b) {
  if (is.function(a) && is.function(b)) {
    identical(deparse(formals(a), control = NULL), deparse(formals(b), control = NULL)) &&
      identical(deparse(body(a), control = NULL), deparse(body(b), control = NULL))
  } else {
    identical(a, b)
  }
}
.copied <- function(path, names) {
  env <- .from_script(path, names)
  vapply(names, function(k) exists(k, envir = env, inherits = FALSE) &&
           .same(get(k, envir = env), get(k, envir = HERE_ENV)), logical(1))
}

E54_COPIED <- c("AXIS_MM", "PT_EPS", "PT_TXT", "PT_TITLE", "CAP_EM", "COL_MARK", "NULL_COL",
                "NULL_LW", "CI_LW", "PT_SIZE", "PT_STROKE", "CIRCLE_MM", "M_EDGE", "G_LAB",
                "KEY_ABOVE", "KEY_STEP", "KEY_GAP", "TITLE_GAP", "TITLE_LH", "PNG_DPI", "TOL_PT",
                "SOFT_HY", ".stop_if", ".assert", "NM_FONT", "NM_TXT_MIN", "NM_TXT", "NM_TXT_TITLE",
                "NM_LINE", "pt", "lab_num", "theme_nm", ".txt_mm", "CAP_MM", ".wrap", ".natural",
                ".pdf_device", ".md5", ".pdf_md5", ".pdf_text", ".pdf_read_back", ".mm2pt")
same54 <- .copied(PATH_E54_SCRIPT, E54_COPIED)
.assert("Fig. 4x's sizes, type, marks and helpers are E54's, object for object",
        all(same54), n = length(same54), detail = paste(names(same54)[!same54], collapse = ", "))

E45_COPIED <- c("PCR_ORDER", "MB_ORDER", "Y_ROWS", "Y_SEP", "Y_CUE", "Y_HEAD", "Y_LIM",
                "KEY_LEVELS", "XLIM_PCR", "BREAKS_PCR", "XLIM_MB", "BREAKS_MB")
same45 <- .copied(PATH_E45_SCRIPT, E45_COPIED)
.assert("the rows, their geometry, the keys' order and the axes are E45's, value for value",
        all(same45), n = length(same45), detail = paste(names(same45)[!same45], collapse = ", "))
E45_SRC <- readLines(PATH_E45_SCRIPT, warn = FALSE)
in45 <- c(vapply(unlist(TEXT, use.names = FALSE),
                 function(s) any(grepl(deparse(s), E45_SRC, fixed = TRUE)), logical(1)),
          vapply(c("lo_a <- 1 / 1.04", "hi_a <- 1.04", "lo_b <- xlim[1] * 1.10",
                   "hi_b <- xlim[2] / 1.10", 'colour = "grey35"', 'colour = "grey30"',
                   'colour = "grey75", linewidth = 0.25'),
                 function(s) any(grepl(s, E45_SRC, fixed = TRUE)), logical(1)))
.assert("every title, header, axis title and cue is E45's, character for character, as are its greys and the cues' ends",
        all(in45) && CUE_IN == 1.04 && CUE_OUT == 1.10 && SUB_COL == "grey35" &&
          CUE_COL == "grey30" && RULE_COL == "grey75" && RULE_LW == 0.25, n = length(in45),
        detail = paste(names(in45)[!in45], collapse = " / "))

COH52 <- get("COH_COL", envir = .from_script(PATH_E52_SCRIPT, "COH_COL"))
.assert("dark slate and mid slate are the PUMA/BCL-xL panel's two colours (E52)",
        identical(unname(COH52), c(COL_OX, COL_PR)), n = 2L)

# =============================================================================
# 2. E45's saved values, checked against their recorded digests
# =============================================================================
message("\n2. E45's saved values")
.stop_if(file.exists(PATH_E45), "E45's saved object is missing: re-source ",
         "scripts/E45_fig4_human_panels.R. STOP.")
e45 <- readRDS(PATH_E45)
GOT45 <- vapply(names(E45_DIGESTS), function(k) digest::digest(e45[[k]]), character(1))
.assert("E45's saved object is the one recorded before its run, in the five elements read here",
        identical(GOT45, E45_DIGESTS), n = 5L,
        detail = paste(names(GOT45)[GOT45 != E45_DIGESTS], collapse = ", "))
.assert("E45's own 27 assertions all passed",
        nrow(e45$assertions) == 27L && all(e45$assertions$pass), n = 27L)
PF <- as.data.frame(e45$plot_frame)
CC <- as.data.frame(e45$contrast_column)
.assert(paste0("E45's frame: 12 pCR marks (four rows, three keys) and 8 METABRIC marks ",
               "(four rows, two keys), in E45's rows"),
        sum(PF$panel == "pCR") == 12L && sum(PF$panel == "METABRIC") == 8L &&
          setequal(PF$key[PF$panel == "pCR"], KEY_LEVELS) &&
          setequal(PF$key[PF$panel == "METABRIC"], KEY_LEVELS[2:3]) &&
          setequal(PF$row_key[PF$panel == "pCR"], c(PCR_ORDER, "Pooled")) &&
          setequal(PF$row_key[PF$panel == "METABRIC"], MB_ORDER), n = 20L)
.assert("every ratio and bound drawn is exp() of E45's stored log value, exactly",
        identical(PF$ratio, exp(PF$log_est)) && identical(PF$ratio_lo, exp(PF$log_lo)) &&
          identical(PF$ratio_hi, exp(PF$log_hi)), n = 60L)
.assert("every contrast printed is E45's stored contrast, exponentiated, to two decimals",
        nrow(CC) == 8L &&
          identical(CC$label, sprintf("%.2f (%.2f, %.2f)", exp(CC$delta), exp(CC$ci_lo),
                                      exp(CC$ci_hi))), n = 8L)
for (p in c("pCR", "METABRIC")) {
  d <- PF[PF$panel == p, ]
  lim <- if (p == "pCR") XLIM_PCR else XLIM_MB
  .assert(paste0(p, ": every interval lies inside E45's axis limits"),
          all(d$ratio_lo > lim[1]) && all(d$ratio_hi < lim[2]), n = nrow(d),
          detail = sprintf("%.3f to %.3f", min(d$ratio_lo), max(d$ratio_hi)))
}

# =============================================================================
# 3. The layout of each panel, in mm
# =============================================================================
# Across, from the page's left edge: the label column, the data area (26 mm),
# the contrast column. Up, from the axis line: E45's geometry times S_MM, so a
# row's centre is E45's y times 5.6 mm. Every string is measured with the PDF
# device's metrics, in the form the PDF prints it.
message("\n3. The layout")
CAP <- CAP_MM(PT_TXT)
.measure <- function(d) {
  m <- t(vapply(seq_len(nrow(d)), function(i) {
    .txt_mm(if (d$parse[i]) d$label[i] else .hy(d$label[i]), d$size_pt[i],
            parse = d$parse[i])[1, ]
  }, numeric(3)))
  d$w <- unname(m[, 1]); d$h <- unname(m[, 2]); d$d <- unname(m[, 3])
  d
}
.text <- function(role, label, base, colour = "black", size_pt = PT_TXT, hjust = 0,
                  parse = FALSE) {
  tibble::tibble(role = role, label = label, base = base, colour = colour,
                 size_pt = size_pt, hjust = hjust, parse = parse)
}

.layout <- function(p) {
  tx   <- TEXT[[p]]
  fp   <- c(pcr = "pCR", metabric = "METABRIC")[[p]]
  xlim <- if (p == "pcr") XLIM_PCR else XLIM_MB
  brk  <- if (p == "pcr") BREAKS_PCR else BREAKS_MB
  keys <- if (p == "pcr") KEY_LEVELS else KEY_LEVELS[2:3]
  rows <- as.data.frame(e45$rows[[p]])
  mk   <- PF[PF$panel == fp, ]
  con  <- CC[CC$panel == fp, ]

  # The label column: each row's two lines, centred on the row; the header's
  # two lines centred on E45's header height.
  yc  <- rows$y * S_MM
  step <- ifelse(rows$row_key == "Pooled", I2_STEP, LINE_STEP)
  b1  <- yc + (step - CAP) / 2
  hb1 <- Y_HEAD * S_MM + (LINE_STEP - CAP) / 2
  head <- strsplit(tx$head, "\n", fixed = TRUE)[[1]]
  lab <- dplyr::bind_rows(
    .text("row", rows$label, b1),
    .text("row_n", rows$sub[!is.na(rows$sub)], b1[!is.na(rows$sub)] - LINE_STEP, SUB_COL),
    .text("head", head, hb1 - c(0, LINE_STEP)))
  if (p == "pcr") {
    pool <- mk[mk$row_key == "Pooled", ]
    i2 <- round(pool$I2[match(KEY_LEVELS, pool$key)])
    lab <- dplyr::bind_rows(lab, .text(
      "i2", sprintf("italic(I)^displaystyle('2')*' %s%%, %s%%, %s%%'", i2[1], i2[2], i2[3]),
      b1[rows$row_key == "Pooled"] - I2_STEP, SUB_COL, parse = TRUE))
  }
  lab <- .measure(lab)
  # The contrast column: its header's two lines and each row's contrast,
  # centred; a contrast's capitals are centred on its row.
  chead <- strsplit(tx$con_head, "\n", fixed = TRUE)[[1]]
  cn <- .measure(dplyr::bind_rows(
    .text("con_head", chead, hb1 - c(0, LINE_STEP), hjust = 0.5),
    .text("contrast", con$label, con$y_row * S_MM - CAP / 2, hjust = 0.5)))
  x_data <- M_EDGE + max(lab$w) + G_LAB
  w_con  <- max(cn$w)
  x_con  <- x_data + AXIS_MM + CON_GAP + w_con / 2
  x_conr <- x_data + AXIS_MM + CON_GAP + w_con
  lab$x_mm <- M_EDGE
  cn$x_mm  <- x_con
  x_at  <- function(mm) xlim[1] * (xlim[2] / xlim[1])^((mm - x_data) / AXIS_MM)
  mm_of <- function(r) x_data + log(r / xlim[1]) / log(xlim[2] / xlim[1]) * AXIS_MM

  # The key above the data area, top to bottom in the marks' order within a row;
  # the title above it.
  nk <- length(keys)
  kt <- .measure(.text("key", unname(KEY_TEXT[keys]),
                       PANEL_H + KEY_ABOVE + (rev(seq_len(nk)) - 1) * KEY_STEP))
  kt$x_mm <- M_EDGE + CIRCLE_MM + KEY_GAP
  key_mk <- tibble::tibble(key = keys, x = x_at(M_EDGE + CIRCLE_MM / 2),
                           y = kt$base + CAP / 2, x_mm = M_EDGE + CIRCLE_MM / 2,
                           fill = unname(KEY_FILL[keys]), line = unname(KEY_LINE[keys]))
  ttl <- .measure(.text("title", tx$title, max(kt$base) + TITLE_ABOVE, size_pt = PT_TITLE))
  ttl$x_mm <- M_EDGE

  # The direction cues: arrows from 1.04 (or 1/1.04) to 1.10 inside the axis's
  # ends, each text anchored at its arrow's inner end and running outward.
  y_cue <- Y_CUE * S_MM
  cue <- tibble::tibble(x = c(1 / CUE_IN, CUE_IN), xend = c(xlim[1] * CUE_OUT, xlim[2] / CUE_OUT),
                        y = y_cue)
  cue_t <- .measure(.text("cue", tx$cue, y_cue + CUE_UP, CUE_COL, hjust = c(1, 0)))
  cue_t$x_mm <- mm_of(c(1 / CUE_IN, CUE_IN))

  # The marks, the rule before the fourth row, and (pCR) the line from each
  # row's OXPHOS without proliferation to OXPHOS with it.
  marks <- tibble::tibble(row_key = mk$row_key, key = mk$key, term = mk$term, rung = mk$rung,
                          ratio = mk$ratio, lo = mk$ratio_lo, hi = mk$ratio_hi,
                          y = mk$y_mark * S_MM, fill = unname(KEY_FILL[mk$key]),
                          line = unname(KEY_LINE[mk$key])) %>%
    dplyr::mutate(x_mm = mm_of(ratio), lo_mm = mm_of(lo), hi_mm = mm_of(hi))
  rule <- tibble::tibble(x = xlim[1], xend = x_at(x_conr), y = Y_SEP * S_MM)
  conn <- tibble::tibble(row_key = character(), x = numeric(), y = numeric(),
                         xend = numeric(), yend = numeric())
  if (p == "pcr") {
    m0 <- marks[marks$key == "OXPHOS, m0", ]
    m1 <- marks[marks$key == "OXPHOS, m1", ]
    m1 <- m1[match(m0$row_key, m1$row_key), ]
    conn <- tibble::tibble(row_key = m0$row_key, x = m0$ratio, y = m0$y,
                           xend = m1$ratio, yend = m1$y)
  }

  text <- dplyr::bind_rows(lab, cn, kt, ttl, cue_t) %>%
    dplyr::mutate(x = x_at(x_mm), y = base - ifelse(parse, d, 0),
                  x0_mm = x_mm - hjust * w, x1_mm = x0_mm + w,
                  y0_mm = base - d, y1_mm = ifelse(parse, base - d + h, base + h))
  top <- max(ttl$base) - PANEL_H + ttl$h + M_EDGE
  list(p = p, fp = fp, xlim = xlim, breaks = brk, keys = keys, x_title = tx$x_title,
       marks = marks, conn = conn, key_marks = key_mk, cue = cue, rule = rule, text = text,
       x_at = x_at, mm_of = mm_of,
       mm = c(data_left = x_data, con_centre = x_con, con_right = x_conr, top = top,
              right = CON_GAP + w_con + M_EDGE, label_w = max(lab$w), con_w = w_con))
}
LAY <- list(pcr = .layout("pcr"), metabric = .layout("metabric"))

# --- the layout's own checks, in mm ------------------------------------------------
.boxes_overlap <- function(a, b, gap = 0) {
  a$x0 < b$x1 + gap & b$x0 < a$x1 + gap & a$y0 < b$y1 + gap & b$y0 < a$y1 + gap
}
for (L in LAY) {
  P <- .P(L$p)
  tb <- data.frame(x0 = L$text$x0_mm, x1 = L$text$x1_mm, y0 = L$text$y0_mm, y1 = L$text$y1_mm)
  pairs <- utils::combn(nrow(tb), 2)
  hit <- vapply(seq_len(ncol(pairs)), function(k)
    .boxes_overlap(tb[pairs[1, k], ], tb[pairs[2, k], ], gap = 0.2), logical(1))
  .assert(paste0(P, ": no two strings come within 0.2 mm of each other"),
          !any(hit), n = ncol(pairs),
          detail = paste(utils::head(apply(pairs[, hit, drop = FALSE], 2, function(k)
            paste(L$text$label[k], collapse = " & ")), 2), collapse = " / "))
  # Marks and intervals as boxes: circles at their diameter, intervals at their
  # printed width; the cue arrows with their heads (0.4 mm either side).
  half_ci <- .w_pt(CI_LW) * 25.4 / 72 / 2
  mb <- rbind(
    data.frame(x0 = L$marks$x_mm - CIRCLE_MM / 2, x1 = L$marks$x_mm + CIRCLE_MM / 2,
               y0 = L$marks$y - CIRCLE_MM / 2, y1 = L$marks$y + CIRCLE_MM / 2),
    data.frame(x0 = L$marks$lo_mm, x1 = L$marks$hi_mm, y0 = L$marks$y - half_ci,
               y1 = L$marks$y + half_ci),
    data.frame(x0 = pmin(L$mm_of(L$cue$x), L$mm_of(L$cue$xend)),
               x1 = pmax(L$mm_of(L$cue$x), L$mm_of(L$cue$xend)),
               y0 = L$cue$y - ARROW_MM * sin(pi / 6), y1 = L$cue$y + ARROW_MM * sin(pi / 6)))
  hit2 <- outer(seq_len(nrow(tb)), seq_len(nrow(mb)), Vectorize(function(i, j)
    .boxes_overlap(tb[i, ], mb[j, ], gap = 0.2)))
  .assert(paste0(P, ": no string comes within 0.2 mm of a mark, an interval or a cue's arrow"),
          !any(hit2), n = length(hit2),
          detail = paste(utils::head(L$text$label[row(hit2)[hit2]], 3), collapse = " / "))
  .assert(paste0(P, ": the label column ends 1.5 mm before the data area, and the contrast ",
                 "column starts 2 mm after it"),
          max(L$text$x1_mm[L$text$role %in% c("row", "row_n", "head", "i2")]) <=
            L$mm[["data_left"]] - G_LAB + 1e-9 &&
            min(L$text$x0_mm[L$text$role %in% c("con_head", "contrast")]) >=
            L$mm[["data_left"]] + AXIS_MM + CON_GAP - 1e-9, n = nrow(L$text))
  page_w <- L$mm[["data_left"]] + AXIS_MM + L$mm[["right"]]
  w_tick <- .txt_mm(lab_num(L$breaks), PT_TXT)[, "w"]
  x_tick <- L$mm_of(L$breaks)
  gaps <- (x_tick[-1] - w_tick[-1] / 2) - (x_tick[-length(x_tick)] + w_tick[-length(w_tick)] / 2)
  .assert(paste0(P, ": the tick labels are at least 0.5 mm apart"),
          all(gaps >= 0.5), n = length(gaps), detail = sprintf("smallest %.2f mm", min(gaps)))
  w_xt <- .txt_mm(.hy(strsplit(L$x_title, "\n", fixed = TRUE)[[1]]), PT_TITLE)[, "w"]
  .assert(paste0(P, ": the axis title, centred on the axis, the key and the title all lie ",
                 "inside the page's margins"),
          L$mm[["data_left"]] + AXIS_MM / 2 - max(w_xt) / 2 >= M_EDGE &&
            L$mm[["data_left"]] + AXIS_MM / 2 + max(w_xt) / 2 <= page_w - M_EDGE &&
            all(L$text$x1_mm[L$text$role %in% c("key", "title")] <= page_w - M_EDGE) &&
            all(L$text$x0_mm >= M_EDGE - 1e-9), n = 3L,
          detail = sprintf("page %.1f mm wide", page_w))
  message(sprintf("   %-8s label column %.1f mm, data area %d x %.2f mm, contrast column %.1f mm",
                  L$p, L$mm[["label_w"]], AXIS_MM, PANEL_H, L$mm[["con_w"]]))
}
.assert("the two panels share E45's row heights, cue height and rule height",
        identical(e45$rows$pcr$y, e45$rows$metabric$y) &&
          identical(LAY$pcr$cue$y, LAY$metabric$cue$y) &&
          identical(LAY$pcr$rule$y, LAY$metabric$rule$y), n = 6L)

# =============================================================================
# 4. The panels
# =============================================================================
message("\n4. The panels")
.panel <- function(L, variant) {
  hy <- function(s) if (variant == "pdf") .hy(s) else s
  tp <- L$text[!L$text$parse, ]
  tp$label <- hy(tp$label)
  te <- L$text[L$text$parse, ]
  g <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 1, colour = NULL_COL, linewidth = NULL_LW) +
    ggplot2::geom_segment(data = L$rule, ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          colour = RULE_COL, linewidth = RULE_LW)
  if (nrow(L$conn)) {
    g <- g + ggplot2::geom_segment(data = L$conn,
                                   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                   colour = COL_OX, linewidth = CONN_LW)
  }
  g <- g +
    ggplot2::geom_segment(data = L$marks, ggplot2::aes(x = lo, xend = hi, y = y, yend = y,
                                                       colour = line), linewidth = CI_LW) +
    ggplot2::geom_point(data = L$marks, ggplot2::aes(x = ratio, y = y, colour = line, fill = fill),
                        shape = 21, size = PT_SIZE, stroke = PT_STROKE) +
    ggplot2::geom_point(data = L$key_marks, ggplot2::aes(x = x, y = y, colour = line, fill = fill),
                        shape = 21, size = PT_SIZE, stroke = PT_STROKE) +
    ggplot2::geom_segment(data = L$cue, ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          colour = CUE_COL, linewidth = CUE_LW, arrow.fill = CUE_COL,
                          arrow = grid::arrow(length = grid::unit(ARROW_MM, "mm"),
                                              type = "closed")) +
    ggplot2::geom_text(data = tp, ggplot2::aes(x = x, y = y, label = label, hjust = hjust,
                                               size = pt(size_pt), colour = colour),
                       vjust = 0, family = NM_FONT)
  if (nrow(te)) {
    g <- g + ggplot2::geom_text(data = te, ggplot2::aes(x = x, y = y, label = label, hjust = hjust,
                                                        size = pt(size_pt), colour = colour),
                                vjust = 0, family = NM_FONT, parse = TRUE)
  }
  g +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_identity() +
    # oob_keep: the labels, the key and the columns sit outside the axes' limits.
    ggplot2::scale_x_continuous(transform = "log10", limits = L$xlim, breaks = L$breaks,
                                labels = lab_num(L$breaks), expand = c(0, 0),
                                oob = scales::oob_keep) +
    ggplot2::scale_y_continuous(limits = c(0, PANEL_H), breaks = NULL, expand = c(0, 0),
                                oob = scales::oob_keep) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = hy(L$x_title), y = NULL) +
    theme_nm(legend = "none", txt = PT_TXT, txt_title = PT_TITLE) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(lineheight = TITLE_LH),
                   plot.margin = ggplot2::margin(L$mm[["top"]], L$mm[["right"]], M_EDGE,
                                                 L$mm[["data_left"]], "mm"))
}
PANEL <- lapply(LAY, function(L) list(pdf = .wrap(.panel(L, "pdf")), png = .wrap(.panel(L, "png"))))

# Each page is the panel's natural size in whole points: the fraction of a point
# that rounding up would leave is added to the top and right margins, so that
# nothing is spread above and below the panel and the two axes sit at the same
# height above their pages' bottom edges.
PAGE_PT <- list(); PAGE <- list()
for (p in names(LAY)) {
  need <- .natural(PANEL[[p]]$pdf) * 72 / 25.4
  pad  <- (ceiling(need - 1e-9) - need) * 25.4 / 72
  LAY[[p]]$mm[["top"]]   <- LAY[[p]]$mm[["top"]] + pad[[2]]
  LAY[[p]]$mm[["right"]] <- LAY[[p]]$mm[["right"]] + pad[[1]]
  PANEL[[p]] <- list(pdf = .wrap(.panel(LAY[[p]], "pdf")), png = .wrap(.panel(LAY[[p]], "png")))
  got <- .natural(PANEL[[p]]$pdf) * 72 / 25.4
  pp <- as.integer(round(got))
  names(pp) <- c("w", "h")
  .assert(paste0(.P(p), ": the page is the panel's natural size in whole points, its top ",
                 "and right margins taking the fraction of a point"),
          all(abs(got - pp) < 1e-6) && all(pad >= 0) && all(pad < 25.4 / 72), n = 2L,
          max_dev = max(abs(got - pp)), detail = sprintf("%.3f and %.3f mm added", pad[1], pad[2]))
  PAGE_PT[[p]] <- pp
  PAGE[[p]] <- (pp + PT_EPS) * 25.4 / 72
  message(sprintf("   %-8s page %d x %d pt (%.1f x %.1f mm)", p, pp[["w"]], pp[["h"]],
                  pp[["w"]] * 25.4 / 72, pp[["h"]] * 25.4 / 72))
}

# Every point and interval ggplot will draw is its stored value, in its place.
for (L in LAY) {
  gp <- PANEL[[L$p]]$pdf[[1]]
  B  <- ggplot2::ggplot_build(gp)
  is_marks <- vapply(gp$layers, function(l) identical(l$data, L$marks), logical(1))
  geo <- vapply(gp$layers, function(l) class(l$geom)[1], character(1))
  pts <- B$data[[which(is_marks & geo == "GeomPoint")]]
  ci  <- B$data[[which(is_marks & geo == "GeomSegment")]]
  dev <- max(abs(c(pts$x - log10(L$marks$ratio), pts$y - L$marks$y,
                   ci$x - log10(L$marks$lo), ci$xend - log10(L$marks$hi), ci$y - L$marks$y)))
  .assert(paste0(.P(L$p), ": every point and interval drawn is its stored ratio and ",
                 "interval, at its row's height (on the log scale)"),
          nrow(pts) == nrow(L$marks) && nrow(ci) == nrow(L$marks) && dev < 1e-12,
          n = 3L * nrow(L$marks), max_dev = dev)
}

# =============================================================================
# 5. The files: the PNGs first, then the PDFs (base pdf(), live Helvetica)
# =============================================================================
message("\n5. The files")
for (d in c(OUT_TABLES, OUT_FIGURES, DIR_DOC_FIGURES, dirname(PATH_E57))) .ensure_dir(d)
FILES <- lapply(stats::setNames(names(LAY), names(LAY)), function(p) c(
  png_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM[[p]], ".png")),
  png_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM[[p]], ".png")),
  pdf_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM[[p]], ".pdf")),
  pdf_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM[[p]], ".pdf"))))
for (p in names(LAY)) {
  for (f in FILES[[p]][c("png_outputs", "png_docs")]) {
    ggplot2::ggsave(f, PANEL[[p]]$png, width = PAGE[[p]][["w"]], height = PAGE[[p]][["h"]],
                    units = "mm", dpi = PNG_DPI, device = "png", type = "cairo", bg = "white")
  }
  for (f in FILES[[p]][c("pdf_outputs", "pdf_docs")]) {
    ggplot2::ggsave(f, PANEL[[p]]$pdf, width = PAGE[[p]][["w"]], height = PAGE[[p]][["h"]],
                    units = "mm", device = .pdf_device)
  }
}

# =============================================================================
# 6. What the files show, read back
# =============================================================================
message("\n6. The files, read back")
AXIS_AT <- list()
for (L in LAY) {
  p <- L$p; P <- .P(p); FL <- FILES[[p]]
  .assert(paste0(P, ": the two PNGs are identical"),
          .md5(FL[["png_outputs"]]) == .md5(FL[["png_docs"]]), n = 2L)
  .assert(paste0(P, ": the two PDFs are identical, timestamps aside"),
          .pdf_md5(FL[["pdf_outputs"]]) == .pdf_md5(FL[["pdf_docs"]]), n = 2L)

  RB <- .pdf_read_back(FL[["pdf_docs"]])
  PA <- .pdf_paint(FL[["pdf_docs"]])
  .assert(paste0(P, ": the PDF's page is the panel's size, in whole points"),
          identical(as.integer(RB$mediabox[3:4]), unname(PAGE_PT[[p]])), n = 1L,
          detail = paste(RB$mediabox, collapse = " "))

  # --- the text ---------------------------------------------------------------
  TX <- RB$placed
  TX$s <- trimws(TX$s)
  TX <- TX[nzchar(TX$s), ]                    # plotmath's spaces are their own strings
  xt_lines <- strsplit(L$x_title, "\n", fixed = TRUE)[[1]]
  plain <- L$text[!L$text$parse, ]
  want <- c(.shy(plain$label), lab_num(L$breaks), .shy(xt_lines))
  if (p == "pcr") {
    i2 <- L$text$label[L$text$parse]
    want <- c(want, "I", "2", sub("^.*\\*' (.*)'$", "\\1", i2))
  }
  .assert(paste0(P, ": the PDF shows exactly the strings it should, and no others"),
          identical(sort(TX$s), sort(want)), n = length(want),
          detail = paste(utils::head(c(setdiff(TX$s, want), setdiff(want, TX$s)), 3),
                         collapse = " / "))
  .assert(paste0(P, ": the PDF's text is live Helvetica, referenced and not embedded, every ",
                 "string was found, and NO STRING IS BOLD"),
          !RB$font_file && nrow(RB$placed) == RB$n_text_blocks &&
            all(TX$font %in% c("Helvetica", "Helvetica-Oblique")) &&
            !any(grepl("Bold", TX$font)), n = nrow(TX),
          detail = paste(unique(TX$font), collapse = ", "))
  .assert(paste0(P, if (p == "pcr") ": only the I of I-squared is oblique" else
                   ": no string is oblique"),
          if (p == "pcr") identical(TX$s[TX$font == "Helvetica-Oblique"], "I") else
            !any(TX$font == "Helvetica-Oblique"), n = nrow(TX))
  big <- .shy(c(L$text$label[L$text$role == "title"], xt_lines))
  .assert(paste0(P, ": the title and the axis title are 6 pt, everything else 5 pt, and ",
                 "nothing is rotated"),
          all(abs(TX$size[TX$s %in% big] - PT_TITLE) < 0.01) &&
            all(abs(TX$size[!TX$s %in% big] - PT_TXT) < 0.01) &&
            sum(TX$s %in% big) == 3L && all(abs(TX$angle) < 1e-6), n = nrow(TX))

  # --- the axis, the line at 1, and the data area's origin in the PDF -------------
  SG <- PA$segments
  hz <- abs(SG$y1 - SG$y0) < 1e-9
  vt <- abs(SG$x1 - SG$x0) < 1e-9
  len <- sqrt((SG$x1 - SG$x0)^2 + (SG$y1 - SG$y0)^2)
  ax <- which(hz & abs(len - .mm2pt(AXIS_MM)) <= TOL_PT & SG$stroke == "#000000")
  .assert(paste0(P, ": in the PDF the x axis measures 26 mm"), length(ax) == 1L, n = 1L)
  AX_X0 <- min(SG$x0[ax], SG$x1[ax]); AX_Y <- SG$y0[ax]
  AXIS_AT[[p]] <- AX_Y
  X_PT <- function(mm) AX_X0 + .mm2pt(mm - L$mm[["data_left"]])   # page mm -> PDF x
  Y_PT <- function(mm) AX_Y + .mm2pt(mm)                          # mm above the axis -> PDF y
  nl <- which(vt & abs(len - .mm2pt(PANEL_H)) <= TOL_PT)
  .assert(paste0(P, ": the line at 1 spans the data area at 1's place, grey45 at 0.64 pt"),
          length(nl) == 1L && abs(SG$x0[nl] - X_PT(L$mm_of(1))) <= TOL_PT &&
            abs(min(SG$y0[nl], SG$y1[nl]) - AX_Y) <= TOL_PT &&
            SG$stroke[nl] == .hex(NULL_COL) && .lw_is(SG$lw[nl], NULL_LW), n = 1L)

  # --- every string in its place ----------------------------------------------------
  ex <- plain %>% dplyr::mutate(s = .shy(label), x_pt = X_PT(x0_mm), y_pt = Y_PT(base))
  dev_t <- unlist(lapply(unique(ex$s), function(s) {
    e <- ex[ex$s == s, ]; g <- TX[TX$s == s, ]
    if (nrow(e) != nrow(g)) return(Inf)
    e <- e[order(-e$y_pt), ]; g <- g[order(-g$y), ]
    c(abs(e$x_pt - g$x), abs(e$y_pt - g$y))
  }))
  .assert(paste0(P, ": in the PDF every label, header, contrast, key line, title and cue ",
                 "sits at its place, to 0.02 pt"),
          all(dev_t <= TOL_PT), n = nrow(ex), max_dev = max(dev_t))
  w_tick <- .txt_mm(lab_num(L$breaks), PT_TXT)[, "w"]
  dev_k <- vapply(seq_along(L$breaks), function(i) {
    g <- TX$x[TX$s == lab_num(L$breaks)[i] & TX$y < AX_Y]   # below the axis: I-squared has a 2 too
    if (length(g) != 1L) return(Inf)
    abs(g - (X_PT(L$mm_of(L$breaks[i])) - .mm2pt(w_tick[i] / 2)))
  }, numeric(1))
  .assert(paste0(P, ": every tick label is centred under its tick"),
          all(dev_k <= TOL_PT), n = length(dev_k), max_dev = max(dev_k))
  y_t <- vapply(.shy(xt_lines), function(s) TX$y[TX$s == s][1], numeric(1))
  dev_xt <- abs(unname(y_t[1] - y_t[2]) - .mm2pt(TITLE_GAP))
  .assert(paste0(P, ": the axis title's two lines are 2.4 mm apart"),
          dev_xt <= TOL_PT, n = 1L, max_dev = dev_xt)
  if (p == "pcr") {
    i2_base <- L$text$base[L$text$parse]
    g_i <- TX[TX$s == "I", ]; g_r <- TX[TX$s == sub("^.*\\*' (.*)'$", "\\1", i2), ]
    g_2 <- TX[TX$s == "2" & TX$y > AX_Y, ]   # the tick label 2 is below the axis
    dev_i2 <- c(abs(g_i$x - X_PT(M_EDGE)), abs(g_i$y - Y_PT(i2_base)), abs(g_r$y - Y_PT(i2_base)))
    .assert(paste0("pCR: I-squared starts at the margin on its row's second line, with its ",
                   "2 raised at the full 5 pt"),
            all(dev_i2 <= TOL_PT) && g_2$y > g_i$y + 1 && abs(g_2$size - PT_TXT) < 0.01, n = 3L,
            max_dev = max(dev_i2))
  }

  # --- the circles: count, size, place and paint ----------------------------------------
  CI_ <- PA$circles
  want_c <- dplyr::bind_rows(
    L$marks %>% dplyr::transmute(what = "mark", key, cx = X_PT(x_mm), cy = Y_PT(y), fill, line),
    L$key_marks %>% dplyr::transmute(what = "key", key, cx = X_PT(x_mm), cy = Y_PT(y), fill, line))
  outer_mm <- (CI_$d + CI_$lw) * 25.4 / 72
  .assert(paste0(P, ": in the PDF the ", nrow(want_c), " circles (", nrow(L$marks), " marks and ",
                 "the key's ", nrow(L$key_marks), ") are 1.05 mm across, outline included"),
          !is.null(CI_) && nrow(CI_) == nrow(want_c) && all(round(outer_mm, 2) == 1.05) &&
            all(abs(CI_$d - CI_$h) <= TOL_PT), n = nrow(want_c),
          max_dev = max(abs(outer_mm - CIRCLE_MM)))
  mt <- vapply(seq_len(nrow(want_c)), function(i) {
    k <- which(abs(CI_$cx - want_c$cx[i]) <= TOL_PT & abs(CI_$cy - want_c$cy[i]) <= TOL_PT)
    if (length(k) != 1L) return(NA_integer_)
    k
  }, integer(1))
  .assert(paste0(P, ": every circle sits at its ratio and row, or its key line, to 0.02 pt"),
          !anyNA(mt) && !any(duplicated(mt)), n = nrow(want_c))
  .assert(paste0(P, ": every circle has its key's paint: OXPHOS dark slate, open (white) ",
                 "without proliferation and filled with it; proliferation mid slate"),
          identical(CI_$fill[mt], .hex(want_c$fill)) && identical(CI_$stroke[mt], .hex(want_c$line)),
          n = 2L * nrow(want_c),
          detail = paste(unique(paste(CI_$fill, CI_$stroke)), collapse = "; "))

  # --- the intervals, the rule, the connector and the cues -----------------------------
  ci_rows <- which(hz & SG$stroke %in% .hex(c(COL_OX, COL_PR)) & .lw_is(SG$lw, CI_LW))
  mi <- vapply(seq_len(nrow(L$marks)), function(i) {
    k <- ci_rows[abs(SG$y0[ci_rows] - Y_PT(L$marks$y[i])) <= TOL_PT]
    if (length(k) != 1L) return(NA_integer_)
    k
  }, integer(1))
  dev_ci <- if (anyNA(mi)) Inf else max(abs(c(
    pmin(SG$x0[mi], SG$x1[mi]) - X_PT(L$marks$lo_mm),
    pmax(SG$x0[mi], SG$x1[mi]) - X_PT(L$marks$hi_mm))))
  .assert(paste0(P, ": in the PDF every interval's ends sit where its bounds fall on the axis, ",
                 "in its key's colour at 0.85 pt"),
          length(ci_rows) == nrow(L$marks) && dev_ci <= TOL_PT &&
            identical(SG$stroke[mi], .hex(L$marks$line)), n = 2L * nrow(L$marks), max_dev = dev_ci)
  ru <- which(hz & SG$stroke == .hex(RULE_COL))
  .assert(paste0(P, ": the rule before the fourth row runs from the data area to the contrast ",
                 "column's end, grey75 at 0.53 pt"),
          length(ru) == 1L && abs(SG$y0[ru] - Y_PT(Y_SEP * S_MM)) <= TOL_PT &&
            abs(min(SG$x0[ru], SG$x1[ru]) - AX_X0) <= TOL_PT &&
            abs(max(SG$x0[ru], SG$x1[ru]) - X_PT(L$mm[["con_right"]])) <= TOL_PT &&
            .lw_is(SG$lw[ru], RULE_LW), n = 1L)
  dg <- which(!hz & !vt)
  if (p == "pcr") {
    dev_cn <- vapply(seq_len(nrow(L$conn)), function(i) {
      a <- c(X_PT(L$mm_of(L$conn$x[i])), Y_PT(L$conn$y[i]))
      b <- c(X_PT(L$mm_of(L$conn$xend[i])), Y_PT(L$conn$yend[i]))
      k <- dg[abs(SG$y0[dg] - a[2]) <= TOL_PT | abs(SG$y1[dg] - a[2]) <= TOL_PT]
      if (length(k) != 1L) return(Inf)
      ends <- rbind(c(SG$x0[k], SG$y0[k]), c(SG$x1[k], SG$y1[k]))
      min(max(abs(c(ends[1, ] - a, ends[2, ] - b))), max(abs(c(ends[1, ] - b, ends[2, ] - a))))
    }, numeric(1))
    .assert(paste0("pCR: each row's line from OXPHOS without to with proliferation joins its ",
                   "two circles' centres, dark slate at 0.43 pt"),
            length(dg) == nrow(L$conn) && all(dev_cn <= TOL_PT) &&
              all(SG$stroke[dg] == .hex(COL_OX)) && all(.lw_is(SG$lw[dg], CONN_LW)),
            n = nrow(L$conn), max_dev = max(dev_cn))
  } else {
    .assert("METABRIC: no diagonal line", length(dg) == 0L, n = 1L)
  }
  cu <- which(hz & SG$stroke == .hex(CUE_COL))
  dev_cu <- if (length(cu) == 2L) {
    o <- order(pmin(SG$x0[cu], SG$x1[cu]))
    c(abs(SG$y0[cu] - Y_PT(Y_CUE * S_MM)),
      abs(sort(c(SG$x0[cu], SG$x1[cu])) -
            sort(X_PT(L$mm_of(c(L$cue$x, L$cue$xend))))))
  } else Inf
  .assert(paste0(P, ": the two cues' arrows run from 1/1.04 and 1.04 to 1.10 inside the axis's ",
                 "ends, grey30 at 0.53 pt"),
          length(cu) == 2L && all(dev_cu <= TOL_PT) && all(.lw_is(SG$lw[cu], CUE_LW)), n = 2L,
          max_dev = max(dev_cu))
  ys_key <- sort(vapply(.shy(L$text$label[L$text$role == "key"]),
                        function(s) TX$y[TX$s == s][1], numeric(1)))
  top_pt <- Y_PT(PANEL_H)
  dev_kk <- abs(c(ys_key[1] - top_pt - .mm2pt(KEY_ABOVE), diff(ys_key) - .mm2pt(KEY_STEP)))
  .assert(paste0(P, ": the key lies above the data area, 1.2 mm then 2.8 mm baseline to baseline"),
          all(dev_kk <= TOL_PT), n = length(dev_kk), max_dev = max(dev_kk))

  # --- the PNG: its size, and every circle's centre in its fill -----------------------
  PNG <- png::readPNG(FL[["png_docs"]])
  h_pt <- PAGE_PT[[p]][["h"]] + PT_EPS
  .assert(paste0(P, ": the PNG is the page's size at 600 dpi"),
          abs(dim(PNG)[1] - PAGE[[p]][["h"]] / 25.4 * PNG_DPI) < 1 &&
            abs(dim(PNG)[2] - PAGE[[p]][["w"]] / 25.4 * PNG_DPI) < 1, n = 1L,
          detail = paste(dim(PNG)[1:2], collapse = " x "))
  px <- vapply(seq_len(nrow(CI_)), function(i) {
    r <- floor((h_pt - CI_$cy[i]) / 72 * PNG_DPI) + 1
    k <- floor(CI_$cx[i] / 72 * PNG_DPI) + 1
    v <- PNG[r, k, 1:3]
    toupper(grDevices::rgb(v[1], v[2], v[3]))
  }, character(1))
  .assert(paste0(P, ": in the PNG every circle's centre shows its fill exactly"),
          identical(px, CI_$fill), n = nrow(CI_),
          detail = paste(utils::head(unique(px[px != CI_$fill]), 3), collapse = ", "))
  rm(PNG)
}

.assert(paste0("the two axes sit at the same height above their pages' bottom edges, so ",
               "the panels aligned by those edges have their axes and rows in line"),
        abs(AXIS_AT$pcr - AXIS_AT$metabric) <= TOL_PT, n = 2L,
        max_dev = abs(AXIS_AT$pcr - AXIS_AT$metabric))

FIGURE_FILES <- dplyr::bind_rows(lapply(names(FILES), function(p) tibble::tibble(
  panel = p, copy = names(FILES[[p]]), file = basename(FILES[[p]]),
  checksum = c(.md5(FILES[[p]][["png_outputs"]]), .md5(FILES[[p]][["png_docs"]]),
               .pdf_md5(FILES[[p]][["pdf_outputs"]]), .pdf_md5(FILES[[p]][["pdf_docs"]])),
  checksum_is = c("md5", "md5", "md5, timestamps zeroed", "md5, timestamps zeroed"))))

# =============================================================================
# 7. The legend drafts: E45's, edited for the new keys, for the text ruling
# =============================================================================
message("\n7. The legend drafts")
.edit <- function(s, old, new) {
  .stop_if(lengths(regmatches(s, gregexpr(old, s, fixed = TRUE))) == 1L,
           "legend edit: '", substr(old, 1, 50), "' is not in E45's draft exactly once. STOP.")
  sub(old, new, s, fixed = TRUE)
}
pool <- PF[PF$panel == "pCR" & PF$row_key == "Pooled", ]
pool <- pool[match(KEY_LEVELS, pool$key), ]
I2_OLD <- paste0("I2 is given for OXPHOS m0, OXPHOS m1 and proliferation (",
                 paste(sprintf("%s %.0f%%, %s", pool$key, pool$I2, pool$summary), collapse = "; "),
                 ").")
I2_NEW <- paste0("I2 is given for the three, in that order (",
                 paste(sprintf("%.0f%%, %s", pool$I2, pool$summary), collapse = "; "), ").")
LEGEND_DRAFT <- list(
  pcr = e45$legend_draft$pcr %>%
    .edit(paste0("Open circle, OXPHOS in a logistic model with receptor subtype and treatment ",
                 "arm where each varies within the cohort (m0); filled circle, the same model ",
                 "with the proliferation score added (m1); triangle, proliferation in that m1 ",
                 "model."),
          paste0("Open circles, OXPHOS, from a logistic model with receptor subtype and ",
                 "treatment arm where each varies within the cohort; dark circles, OXPHOS, ",
                 "adjusted for proliferation, from the same model with the proliferation score ",
                 "added; light circles, proliferation, adjusted for OXPHOS, from that model.")) %>%
    .edit(I2_OLD, I2_NEW) %>%
    .edit("odds ratio within each m1 model", "odds ratio within each adjusted model"),
  metabric = e45$legend_draft$metabric %>%
    .edit("OXPHOS and proliferation (m1). Rows:", "OXPHOS and proliferation. Rows:") %>%
    .edit("Circle, OXPHOS; triangle, proliferation.",
          paste0("Dark circles, OXPHOS, adjusted for proliferation; light circles, ",
                 "proliferation, adjusted for OXPHOS.")),
  model_terms = e45$legend_draft$model_terms)
.assert("the legend drafts are ASCII and name no triangle, m0 or m1",
        !any(grepl("[^ -~]", unlist(LEGEND_DRAFT[c("pcr", "metabric")]))) &&
          !any(grepl("triangle|\\bm0\\b|\\bm1\\b", unlist(LEGEND_DRAFT[c("pcr", "metabric")]))),
        n = 2L)

# =============================================================================
# 8. Saved: the object, two tables and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

DRAWN <- dplyr::bind_rows(lapply(LAY, function(L) L$marks %>%
  dplyr::mutate(panel = L$p, .before = 1) %>% dplyr::select(-fill, -line) %>%
  dplyr::mutate(paint = ifelse(key == "OXPHOS, m0", "open dark slate",
                               ifelse(key == "OXPHOS, m1", "dark slate", "mid slate")),
                key_text = unname(KEY_TEXT[key]))))
CONTRAST <- dplyr::bind_rows(lapply(LAY, function(L) L$text %>%
  dplyr::filter(role == "contrast") %>%
  dplyr::transmute(panel = L$p, label, base_mm = base)))
LAYOUT <- lapply(LAY, function(L) list(
  mm = L$mm, page_pt = PAGE_PT[[L$p]], page_mm = PAGE[[L$p]],
  text = L$text %>% dplyr::select(role, label, size_pt, colour, hjust, parse, x_mm, base,
                                  w, h, d, x0_mm, x1_mm, y0_mm, y1_mm),
  key_marks = L$key_marks, cue = L$cue, rule = L$rule, conn = L$conn))
SPEC <- list(
  ruling = "docs/2026-11-06_E57_outcome_forest_panels_ruling.md (cb3c7d3)",
  source = list(object = "results/e45_fig4_human_panels.rds", digests = E45_DIGESTS,
                recorded_in = "docs/2026-10-22_fig4_human_panels_proposal.md, section 7.6"),
  colours = c(oxphos = COL_OX, proliferation = COL_PR, open = "white"),
  keys = KEY_TEXT, axis_mm = AXIS_MM, panel_height_mm = PANEL_H, row_unit_mm = S_MM,
  x = list(pcr = list(limits = XLIM_PCR, breaks = BREAKS_PCR),
           metabric = list(limits = XLIM_MB, breaks = BREAKS_MB), transform = "log10"),
  type_pt = c(text = PT_TXT, titles = PT_TITLE),
  marks = list(null = c(NULL_COL, NULL_LW), interval = CI_LW, point = c(PT_SIZE, PT_STROKE),
               circle_mm = CIRCLE_MM, connector = CONN_LW, rule = c(RULE_COL, RULE_LW),
               cue = c(CUE_COL, CUE_LW, ARROW_MM), axes = NM_LINE),
  text = TEXT, greys = c(second_lines = SUB_COL, cues = CUE_COL))

OUT <- list(
  frame = DRAWN,
  contrast = CONTRAST,
  layout = LAYOUT,
  legend_draft = LEGEND_DRAFT,
  assertions = ASSERTIONS,
  figure_files = FIGURE_FILES,
  spec = SPEC,
  run_info = list(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  R = R.version.string,
                  bitmap_type_default = getOption("bitmapType")))
saveRDS(OUT, PATH_E57)
utils::write.csv(DRAWN, file.path(OUT_TABLES, "E57_drawn_values.csv"), row.names = FALSE)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E57_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), "run_info"),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-13s %s", k, DIGESTS[[k]]))
message("\nE57 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figures.
  library(dplyr)
  x <- readRDS(here::here("results", "e57_outcome_forest_panels.rds"))
  x$frame %>% dplyr::select(panel, row_key, key_text, ratio, lo, hi, paint) %>% print(n = Inf)
  x$contrast %>% print()
  x$assertions %>% print(n = Inf)
  cat(x$legend_draft$pcr, "\n\n", x$legend_draft$metabric, "\n")
  for (stem in c("E57_pcr_forest", "E57_metabric_survival_forest")) {
    img <- png::readPNG(here::here("docs", "figures", paste0("2026-11-06_", stem, ".png")))
    grid::grid.newpage()
    grid::grid.raster(img)
  }
}
