# E55_raw_correlations_heatmap.R
# =============================================================================
# EXTENDED DATA FIG. 6c: THE UNADJUSTED SPEARMAN CORRELATIONS OF FIG. 4g'S
# FIFTEEN GENES AND THE SIGNATURE WITH ALL SEVEN MEASURES, IN BOTH COHORTS,
# DRAWN AS FIG. 4g IS: ITS ROWS, GROUPS, HEATMAP SIZE AND COLOURS.
#
# A DISPLAY SCRIPT. It fits nothing and computes no new quantity. It reads
# E49's verified saved correlations (the left side of E49's raw-beside-adjusted
# heatmap, ED-h; digests committed at a63ef6c) and draws them with Fig. 4g's
# grammar, taken from E50 (run and verified; its constants are checked against
# its script here).
#
# Specification: the author's ruling of 2026-11-05, verbatim in
# docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md, committed alone at dd5004a
# before this script was written.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNG is compared byte for
# byte with the dry run.
#
# =============================================================================
# WHAT IT DRAWS, AND WHY
# =============================================================================
#   - Rows: Fig. 4g's fifteen genes in its four groups and order, then one
#     signature row (S6, "Signature") at the bottom (the ruling). E49 drew three
#     signatures at the top; S5 and S7 are left out.
#   - Columns: all seven measures of E49's left side (the ruling, question 1),
#     in two cohort blocks, TCGA then SCAN-B, as E49 drew them: OXPHOS level,
#     OXPHOS priority, the four MYC activity scores and proliferation. The MYC
#     scores carry Fig. 6a's names, with the share of their genes in
#     proliferation gene sets.
#   - The heatmap takes Fig. 4g's size: its fourteen columns share 4g's width,
#     and its rows keep 4g's height and gaps. The cells are colour only; at
#     this width they cannot hold printed values (the ruling, question 1).
#   - The colours are 4g's: Fig. 1g's palette at +/-0.6. Correlations beyond
#     +/-0.6 take the end colour; the script lists them.
#   - Labels and the colour bar are placed as E50 places 4g's. The column
#     labels stand vertically above their columns, as in E49.
#   - The PDF carries live Helvetica, for Illustrator; the PNG is a preview.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing and reads no expression data.
#   - It changes no text of the paper. Its legend draft is Claude Code's.
#   - It leaves no raw file bytes in the workspace: the PDF is read back inside
#     a function that returns only ASCII and numbers (E51's lesson).
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
for (pkg in c("digest", "png", "scales", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E55 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE55: Extended Data Fig. 6c, the unadjusted correlations\n", strrep("=", 78))
message("Reads E49's verified object and fits nothing. Specification dd5004a.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_E49        <- file.path(DIR_RESULTS, "e49_signature_and_fig4g.rds")
PATH_E48        <- file.path(DIR_RESULTS, "e48_separation_and_fig4g.rds")
PATH_E50_SCRIPT <- here::here("scripts", "E50_fig4g_in_fig1g_style.R")
PATH_RULING     <- here::here("docs", "2026-11-05_E55_E56_ed_fig6cd_ruling.md")
PATH_E55        <- file.path(DIR_RESULTS, "e55_raw_correlations_heatmap.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-05"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- "E55_raw_correlations_heatmap"

# =============================================================================
# 0. Constants
# =============================================================================
# E49's elements read here, with the digests in docs/2026-11-03_e49_prerun_digests.md
# (a63ef6c), which E49's run reproduced; and E48's display frames, for Fig. 6a's
# names of the measures (docs/2026-11-03_e48_prerun_digests.md, 1ccf805).
E49_DIGESTS <- c(raw_correlations = "f9f32a13550b07075b1690cb7d0a6a36",
                 legend_draft     = "6b7639c40dcde9ae894183f76b367829",
                 assertions       = "2710612444bdbe6592173817b73c2d86")
E48_DIGESTS <- c(displays = "c8397b9106fbe770476a7b10bb513f14")

COHORTS   <- c("TCGA", "SCAN-B")
EXPOSURES <- c("OX_gsva", "OX_mitopps", "MYC_clean", "MYC_M_a", "MYC_entangled",
               "MYC_regulon", "PROLIF")                 # E49's order, and E48's
# The colour bar's title, in Fig. 6a's words ("Spearman rho, absolute"). On one
# line it is 11.2 mm wide and would widen the page past Fig. 4g's 87.6 mm, since
# 4g's "SD per SD" is 8.4 mm; on two lines it is 8.2 mm.
BAR_TITLE <- c("Spearman", "rho")
TITLE_LEAD <- 1.15 * 5                                    # pt: its two lines, 5 pt at 1.15
W_MAX_MM  <- 88                                           # half width, as Fig. 4g's 87.6 mm
L_SCALE   <- 0.6                                         # Fig. 4g's scale (E50)
COL_GAP   <- 1.2                                         # pt: a column label's start above the cells

# --- Fig. 4g's rows, palette and geometry: E50's lines, verbatim -------------------
GROUPS <- list(Effectors = c("BAX", "BAK1", "BOK"),
               Activators = c("BID", "BBC3", "BCL2L11"),
               Sensitisers = c("BAD", "PMAIP1", "BMF", "BIK"),
               `Anti-apoptotic` = c("MCL1", "BCL2L2", "BCL2A1", "BCL2",
                                    "BCL2L1"))
FIFTEEN   <- unlist(GROUPS, use.names = FALSE)
SIG_ROW   <- "S6"
SIG_LABEL <- "Signature"
DIV_COLS <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")
MM       <- 72 / 25.4
CELL_W   <- 26.71        # item 3; Fig. 1g's tiles, measured from its code
CELL_H   <- 8.46
GAP      <- 0.8 * MM     # Fig. 1g: panel.spacing.y = 0.8 mm (2.27 pt)
GAP_SIG  <- 2 * GAP      # "a larger gap" before the signature: twice (a choice)
FONT     <- "Helvetica"
PT_TXT   <- 5            # values and every label (item 3)
PT_HEAD  <- 6            # the score headers (item 3)
CAP      <- 0.718        # Helvetica's cap height, em (AFM)
DESC     <- 0.207        # Helvetica's descender, em (AFM)
LAB_GAP  <- 1.2          # gene symbols end 1.2 pt left of the cells (Fig. 1g)
COH_GAP  <- 1.2          # cohort labels' baseline 1.2 pt above the cells
GRP_GAP  <- 2.4          # group names start 2.4 pt right of the cells
HEAD_GAP <- 2.4          # score header above the cohort labels (strip margin)
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
MINUS_U  <- intToUtf8(8722L) # the minus sign, in the PNG (E50: "\u2212")
SOFT_HY  <- intToUtf8(173L)  # the PDF's hyphen; its encoding prints "-" as minus (E50: "\u00ad")
HJ <- c(start = 0, middle = 0.5, end = 1)

# 6c's cells: Fig. 4g's heatmap is six columns of CELL_W in three blocks; 6c
# puts fourteen columns, in two blocks, into the same width.
HEAT_W_4G <- 6 * CELL_W + 2 * GAP
HEAT_H_4G <- 16 * CELL_H + 3 * GAP + GAP_SIG
CELL_W6   <- (HEAT_W_4G - GAP) / (length(COHORTS) * length(EXPOSURES))

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

.assert("the panel is the ruling's: all seven columns, as its own file records",
        file.exists(PATH_RULING) &&
          any(grepl("All seven columns", readLines(PATH_RULING, warn = FALSE),
                    fixed = TRUE)), n = 1L)
E50_SRC <- readLines(PATH_E50_SCRIPT, warn = FALSE)
E50_LINES <- c(
  'GROUPS <- list(Effectors = c("BAX", "BAK1", "BOK"),',
  '               Activators = c("BID", "BBC3", "BCL2L11"),',
  '               Sensitisers = c("BAD", "PMAIP1", "BMF", "BIK"),',
  '               `Anti-apoptotic` = c("MCL1", "BCL2L2", "BCL2A1", "BCL2",',
  'SIG_ROW   <- "S6"', 'SIG_LABEL <- "Signature"',
  'DIV_COLS <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")',
  'CELL_W   <- 26.71        # item 3; Fig. 1g\'s tiles, measured from its code',
  'CELL_H   <- 8.46',
  'GAP      <- 0.8 * MM     # Fig. 1g: panel.spacing.y = 0.8 mm (2.27 pt)',
  'GAP_SIG  <- 2 * GAP      # "a larger gap" before the signature: twice (a choice)',
  'PT_TXT   <- 5            # values and every label (item 3)',
  'PT_HEAD  <- 6            # the score headers (item 3)',
  'LAB_GAP  <- 1.2          # gene symbols end 1.2 pt left of the cells (Fig. 1g)',
  'GRP_GAP  <- 2.4          # group names start 2.4 pt right of the cells',
  'HEAD_GAP <- 2.4          # score header above the cohort labels (strip margin)',
  'MARGIN   <- c(top = 1, right = 2.5, bottom = 1, left = 1) * MM  # plot.margin',
  'BAR_W    <- 2 * MM       # Fig. 1g\'s colour bar', 'BAR_H    <- 18 * MM',
  'KEY_GAP  <- 1 * MM       # 1 mm past the wider of the first two group names',
  'TITLE_GAP <- 2.4         # bar title\'s descender to the bar (Fig. 1g)',
  'TICK_LEN <- 0.6 * MM     # 1.70 pt, Fig. 1g\'s', 'TICK_LWD <- 0.375        # pt',
  'TICK_COL <- "#4D4D4D"    # grey30',
  'TICKLAB_GAP <- 3.0       # bar to tick labels (Fig. 1g: 2.99 pt)',
  'N_STOPS  <- 41L          # gradient stops on the bar, each the palette\'s colour')
.assert("Fig. 4g's rows, palette and geometry are E50's, line for line",
        all(E50_LINES %in% E50_SRC), n = length(E50_LINES),
        detail = paste(utils::head(setdiff(E50_LINES, E50_SRC), 2), collapse = " / "))

# =============================================================================
# 1. E49's and E48's objects, checked against their committed digests
# =============================================================================
message("\n1. E49's and E48's saved objects")
.stop_if(file.exists(PATH_E49) && file.exists(PATH_E48),
         "E49's or E48's saved object is missing. STOP.")
e49 <- readRDS(PATH_E49)
e48 <- readRDS(PATH_E48)
GOT49 <- vapply(names(E49_DIGESTS), function(k) digest::digest(e49[[k]]), character(1))
GOT48 <- vapply(names(E48_DIGESTS), function(k) digest::digest(e48[[k]]), character(1))
.assert(paste0("E49's and E48's objects are the verified ones: every element read ",
               "here matches its digest committed before its run"),
        identical(unname(GOT49), unname(E49_DIGESTS)) &&
          identical(unname(GOT48), unname(E48_DIGESTS)), n = 4L)
.assert("E49's 192 assertions all passed",
        nrow(e49$assertions) == 192L && all(e49$assertions$pass), n = 192L)

# =============================================================================
# 2. The values: E49's correlations, as E49 drew them
# =============================================================================
message("\n2. values")
RAW <- as.data.frame(e49$raw_correlations)
V <- RAW %>% dplyr::filter(readout %in% c(FIFTEEN, SIG_ROW), exposure %in% EXPOSURES,
                           cohort %in% COHORTS)
DRAWN49 <- as.data.frame(e49$displays$ed_h_raw) %>%
  dplyr::select(cohort, readout, exposure, rho_drawn = rho)
VJ <- dplyr::inner_join(V, DRAWN49, by = c("cohort", "readout", "exposure"))
.assert(paste0("224 correlations, sixteen rows by seven measures by two cohorts, ",
               "each the one E49 drew in ED-h"),
        nrow(V) == 224L && nrow(VJ) == 224L && identical(VJ$rho, VJ$rho_drawn) &&
          !anyNA(V$rho) && all(abs(V$rho) <= 1), n = 224L)
EXPO_LABEL <- stats::setNames(levels(e48$displays$ed_a$x), EXPOSURES)
.assert(paste0("the measures carry Fig. 6a's names, from E48's frame, in E49's ",
               "order"),
        length(levels(e48$displays$ed_a$x)) == 7L &&
          identical(unname(EXPO_LABEL), c("OXPHOS level", "OXPHOS priority",
                                          "MYC_UP.V1_UP (1.5%)", "Felsher, reference (14.8%)",
                                          "Hallmark MYC targets (23.5%)",
                                          "CollecTRI MYC regulon", "Proliferation")),
        n = 7L)
BEYOND <- V %>% dplyr::filter(abs(rho) > L_SCALE) %>% dplyr::arrange(dplyr::desc(abs(rho)))
message(sprintf("   %d correlations lie beyond +/-%.1f and take the end colour: %s",
                nrow(BEYOND), L_SCALE,
                paste(sprintf("%s %s %s %+.2f", BEYOND$readout, BEYOND$exposure,
                              BEYOND$cohort, BEYOND$rho), collapse = "; ")))

# =============================================================================
# 3. Fig. 1g's palette and E50's helpers, verbatim
# =============================================================================
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
.assert(paste0("the palette's anchors are Fig. 1g's, and every cell's colour is ",
               "the one Fig. 1g's own scale call gives its value at +/-0.6 ",
               "(squished beyond it)"),
        identical(.fill_hex(c(-L_SCALE, 0, L_SCALE), L_SCALE), unname(DIV_COLS)) &&
          identical(.fill_hex(V$rho, L_SCALE), .ggplot_fill(V$rho, L_SCALE)) &&
          all(.fill_hex(BEYOND$rho, L_SCALE) %in% unname(DIV_COLS[c("neg", "pos")])),
        n = nrow(V), detail = sprintf("%d squished", nrow(BEYOND)))

# =============================================================================
# 4. The layout, in points from the top left, as E50 lays out 4g
# =============================================================================
message("\n4. layout")
.layout_6c <- function(V) {
  rows <- tibble::tibble(readout = c(FIFTEEN, SIG_ROW),
                         group = c(rep(names(GROUPS), lengths(GROUPS)), SIG_LABEL)) %>%
    dplyr::mutate(k = dplyr::row_number(),
                  g = match(group, c(names(GROUPS), SIG_LABEL)),
                  label = ifelse(readout == SIG_ROW, SIG_LABEL, readout),
                  face = ifelse(readout == SIG_ROW, 1L, 3L),
                  role = ifelse(readout == SIG_ROW, "signature", "gene"))
  lab_w <- max(.text_w(.to_pdf(rows$label), PT_TXT, rows$face))
  col_w <- .text_w(.to_pdf(unname(EXPO_LABEL)), PT_TXT, 1L)
  X0 <- MARGIN[["left"]] + lab_w + LAB_GAP
  Y_COH <- MARGIN[["top"]] + CAP * PT_HEAD              # the cohort headers' baseline
  Y_LAB_TOP <- Y_COH + DESC * PT_HEAD + HEAD_GAP         # the longest column label's top
  Y0 <- Y_LAB_TOP + max(col_w) + COL_GAP                 # the cells' top
  rows <- rows %>%
    dplyr::mutate(y = Y0 + (k - 1) * CELL_H + pmin(g - 1L, 3L) * GAP +
                    (g == 5L) * GAP_SIG)
  cols <- tidyr::expand_grid(cohort = COHORTS, exposure = EXPOSURES) %>%
    dplyr::mutate(c = match(cohort, COHORTS), e = match(exposure, EXPOSURES),
                  x = X0 + (c - 1) * (length(EXPOSURES) * CELL_W6 + GAP) +
                    (e - 1) * CELL_W6)
  X1 <- max(cols$x) + CELL_W6
  rects <- V %>%
    dplyr::inner_join(rows %>% dplyr::select(readout, y), by = "readout") %>%
    dplyr::inner_join(cols %>% dplyr::select(cohort, exposure, x),
                      by = c("cohort", "exposure")) %>%
    dplyr::arrange(match(readout, rows$readout), match(cohort, COHORTS),
                   match(exposure, EXPOSURES)) %>%
    dplyr::mutate(w = CELL_W6, h = CELL_H, fill = .fill_hex(rho, L_SCALE))
  grp <- rows %>% dplyr::filter(g <= 4L) %>% dplyr::group_by(group) %>%
    dplyr::summarise(mid = (min(y) + max(y) + CELL_H) / 2, .groups = "drop") %>%
    dplyr::arrange(match(group, names(GROUPS)))
  blocks <- cols %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(x = (min(x) + max(x) + CELL_W6) / 2, .groups = "drop") %>%
    dplyr::arrange(match(cohort, COHORTS))
  # The colour bar, placed as E50 places 4g's.
  L <- L_SCALE
  XG <- X1 + GRP_GAP
  XB <- XG + max(.text_w(c("Effectors", "Activators"), PT_TXT, 1L)) + KEY_GAP
  Y_TITLE <- Y0 + CAP * PT_TXT
  YB_TOP  <- Y_TITLE + TITLE_LEAD + DESC * PT_TXT + TITLE_GAP
  YB_BOT  <- YB_TOP + BAR_H
  TICKS   <- c(-L, -L / 2, 0, L / 2, L)
  y_tick  <- YB_BOT - (TICKS + L) / (2 * L) * BAR_H
  # A vertical column label reads upwards from just above its column. Its
  # baseline runs up the column, CAP / 2 right of the column's centre, so that
  # its capitals are centred on it.
  texts <- dplyr::bind_rows(
    .txt(rows$role, X0 - LAB_GAP, rows$y + CELL_H / 2 + CAP * PT_TXT / 2, "end",
         PT_TXT, rows$face, rows$label, key = rows$readout),
    .txt("group", XG, grp$mid + CAP * PT_TXT / 2, "start", PT_TXT, 1L,
         grp$group, key = grp$group),
    .txt("cohort", blocks$x, Y_COH, "middle", PT_HEAD, 1L, blocks$cohort,
         key = blocks$cohort),
    .txt("column", cols$x + CELL_W6 / 2 + CAP * PT_TXT / 2, Y0 - COL_GAP, "start",
         PT_TXT, 1L, unname(EXPO_LABEL[cols$exposure]),
         key = paste(cols$cohort, cols$exposure)),
    .txt("bar_title", XB, Y_TITLE + c(0, TITLE_LEAD), "start", PT_TXT, 1L, BAR_TITLE),
    .txt("tick", XB + BAR_W + TICKLAB_GAP, y_tick + CAP * PT_TXT / 2, "start",
         PT_TXT, 1L, sub("^-", MINUS_U, lab_signed(TICKS)),
         key = sprintf("%+.2f", TICKS))) %>%
    dplyr::mutate(rot = ifelse(role == "column", 90, 0),
                  label_pdf = .to_pdf(label_uni),
                  width = .text_w(label_pdf, size, face),
                  left = ifelse(rot == 90, x - CAP * size, x - unname(HJ[anchor]) * width),
                  right = ifelse(rot == 90, x + DESC * size,
                                 x - unname(HJ[anchor]) * width + width),
                  top = ifelse(rot == 90, y - width, y - CAP * size),
                  bottom = ifelse(rot == 90, y, y + DESC * size))
  bar <- tibble::tibble(x = XB, y = YB_TOP, w = BAR_W, h = BAR_H)
  stops <- tibble::tibble(offset = seq(0, 1, length.out = N_STOPS)) %>%
    dplyr::mutate(value = -L + offset * 2 * L, colour = .fill_hex(value, L))
  ticks <- dplyr::bind_rows(
    tibble::tibble(edge = "left", value = TICKS, x0 = XB, x1 = XB + TICK_LEN,
                   y0 = y_tick, y1 = y_tick),
    tibble::tibble(edge = "right", value = TICKS, x0 = XB + BAR_W - TICK_LEN,
                   x1 = XB + BAR_W, y0 = y_tick, y1 = y_tick))
  W <- max(texts$right, bar$x + bar$w) + MARGIN[["right"]]
  H <- max(rects$y + rects$h, YB_BOT, texts$bottom) + MARGIN[["bottom"]]
  list(L = L, rows = rows, cols = cols, rects = rects, texts = texts, bar = bar,
       stops = stops, ticks = ticks, size = c(W = W, H = H), x0 = X0, x1 = X1,
       y0 = Y0)
}
LAY <- .layout_6c(V)
r <- LAY$rects
tx <- LAY$texts
.assert(paste0("224 cells; the heatmap has Fig. 4g's width and height, with its ",
               "rows 8.46 pt tall and its fourteen columns sharing 4g's width"),
        nrow(r) == 224L && all(r$w == CELL_W6) && all(r$h == CELL_H) &&
          abs((LAY$x1 - LAY$x0) - HEAT_W_4G) < 1e-9 &&
          abs((max(r$y + r$h) - LAY$y0) - HEAT_H_4G) < 1e-9, n = 224L,
        detail = sprintf("%.2f x %.2f mm; columns %.2f mm", (LAY$x1 - LAY$x0) / MM,
                         HEAT_H_4G / MM, CELL_W6 / MM))
ry <- LAY$rows %>% dplyr::arrange(k) %>%
  dplyr::mutate(gap = dplyr::lead(y) - (y + CELL_H), g_next = dplyr::lead(g))
within  <- ry$gap[!is.na(ry$gap) & ry$g == ry$g_next]
between <- ry$gap[!is.na(ry$gap) & ry$g != ry$g_next & ry$g_next <= 4L]
before  <- ry$gap[!is.na(ry$gap) & ry$g_next == 5L]
.assert(paste0("Fig. 4g's rows: no gap within a group, 0.8 mm between groups, ",
               "twice that before the signature, which comes last"),
        length(within) == 11L && all(abs(within) < 1e-9) && length(between) == 3L &&
          all(abs(between - GAP) < 1e-9) && length(before) == 1L &&
          abs(before - GAP_SIG) < 1e-9 &&
          identical(LAY$rows$readout[order(LAY$rows$y)], c(FIFTEEN, SIG_ROW)), n = 16L)
cx <- LAY$cols %>% dplyr::arrange(c, e) %>%
  dplyr::mutate(gap = dplyr::lead(x) - (x + CELL_W6), c_next = dplyr::lead(c))
.assert(paste0("the seven measures abut within each cohort's block, TCGA then ",
               "SCAN-B, with 0.8 mm between the blocks"),
        sum(cx$c == cx$c_next, na.rm = TRUE) == 12L &&
          all(abs(cx$gap[!is.na(cx$gap) & cx$c == cx$c_next]) < 1e-9) &&
          abs(cx$gap[!is.na(cx$gap) & cx$c != cx$c_next] - GAP) < 1e-9 &&
          identical(cx$exposure, rep(EXPOSURES, 2)) &&
          identical(cx$cohort, rep(COHORTS, each = 7)), n = 14L)
exp_size <- ifelse(tx$role == "cohort", PT_HEAD, PT_TXT)
.assert(paste0("Helvetica 5 pt, with 6 pt for the cohort headers; oblique for ",
               "the gene symbols and nothing else; the fourteen column labels vertical"),
        all(tx$size == exp_size) && all((tx$face == 3L) == (tx$role == "gene")) &&
          sum(tx$rot == 90) == 14L && all(tx$role[tx$rot == 90] == "column"),
        n = nrow(tx))
ov <- outer(seq_len(nrow(tx)), seq_len(nrow(tx)), function(i, j)
  tx$left[i] < tx$right[j] & tx$left[j] < tx$right[i] &
    tx$top[i] < tx$bottom[j] & tx$top[j] < tx$bottom[i])
n_ov <- sum(ov[upper.tri(ov)])
bb <- LAY$bar
hit_bar <- tx$left < bb$x + bb$w & bb$x < tx$right & tx$top < bb$y + bb$h & bb$y < tx$bottom
cc <- outer(seq_len(nrow(tx)), seq_len(nrow(r)), function(i, k)
  tx$left[i] < r$x[k] + r$w[k] & r$x[k] < tx$right[i] &
    tx$top[i] < r$y[k] + r$h[k] & r$y[k] < tx$bottom[i])
.assert(paste0("no two labels overlap; no label touches a cell or the colour ",
               "bar; everything lies inside the page"),
        n_ov == 0L && !any(hit_bar) && sum(cc) == 0L && all(tx$left >= 0) &&
          all(tx$right <= LAY$size[["W"]]) && all(tx$top >= 0) &&
          all(tx$bottom <= LAY$size[["H"]]), n = nrow(tx),
        detail = sprintf("%d overlaps, %d label-cell contacts", n_ov, sum(cc)))
tk <- tx %>% dplyr::filter(role == "tick")
.assert(paste0("the colour bar is 4g's: 2 x 18 mm, ticked at -0.6, -0.3, 0, ",
               "+0.3 and +0.6, through Fig. 1g's palette; titled \"Spearman\" over ",
               "\"rho\""),
        abs(bb$w - 2 * MM) < 1e-9 && abs(bb$h - 18 * MM) < 1e-9 &&
          identical(tk$label_uni, sub("^-", MINUS_U, lab_signed(c(-0.6, -0.3, 0, 0.3, 0.6)))) &&
          identical(tx$label_uni[tx$role == "bar_title"], BAR_TITLE) &&
          identical(LAY$stops$colour[c(1L, (N_STOPS + 1L) / 2L, N_STOPS)], unname(DIV_COLS)),
        n = nrow(tk))
.assert("the page fits half width, as Fig. 4g's does: at most 88 mm wide",
        LAY$size[["W"]] / MM <= W_MAX_MM, n = 1L,
        detail = sprintf("%.2f mm", LAY$size[["W"]] / MM))
message(sprintf("   page %.1f x %.1f mm; heatmap %.1f x %.1f mm (4g's); columns %.2f mm",
                LAY$size[["W"]] / MM, LAY$size[["H"]] / MM, HEAT_W_4G / MM,
                HEAT_H_4G / MM, CELL_W6 / MM))

# =============================================================================
# 5. The files: the PNG, then the PDF, each read back
# =============================================================================
message("\n5. The files")
.draw <- function(lay, variant) {
  H <- lay$size[["H"]]
  bp <- function(v) grid::unit(v, "bigpts")
  grid::grid.newpage()
  r <- lay$rects
  grid::grid.rect(x = bp(r$x), y = bp(H - r$y), width = bp(r$w),
                  height = bp(r$h), just = c("left", "top"),
                  gp = grid::gpar(fill = r$fill, col = NA))
  b <- lay$bar
  grid::grid.rect(x = bp(b$x), y = bp(H - b$y), width = bp(b$w),
                  height = bp(b$h), just = c("left", "top"),
                  gp = grid::gpar(col = NA, fill = grid::linearGradient(
                    colours = lay$stops$colour, stops = lay$stops$offset,
                    x1 = 0.5, y1 = 0, x2 = 0.5, y2 = 1)))
  tk <- lay$ticks
  grid::grid.segments(x0 = bp(tk$x0), y0 = bp(H - tk$y0), x1 = bp(tk$x1),
                      y1 = bp(H - tk$y1),
                      gp = grid::gpar(col = TICK_COL, lwd = TICK_LWD * 96 / 72,
                                      lineend = "butt"))
  tx <- lay$texts
  lab <- if (identical(variant, "pdf")) tx$label_pdf else tx$label_uni
  for (i in seq_len(nrow(tx))) {
    grid::grid.text(lab[i], x = bp(tx$x[i]), y = bp(H - tx$y[i]),
                    hjust = unname(HJ[tx$anchor[i]]), vjust = 0, rot = tx$rot[i],
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
# The PNG through cairo, which writes each colour exactly (E50's finding).
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
# The PDF is read back INSIDE a function: E50's reader, which follows the fill
# colour into each rectangle, with E53's text reader, which handles rotated
# text, escaped parentheses and kerned arrays. Bytes above 127 are mapped to
# ASCII before anything else (the soft hyphen to "~"), so only ASCII and
# numbers come back (E51's lesson).
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
    m <- regmatches(l, regexec(
      "^([0-9.-]+) ([0-9.-]+) ([0-9.-]+) ([0-9.-]+) re$", l))[[1]]
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
        label = paste(body, collapse = ""), fill = fill)
    }
  }
  list(mediabox = mb_num, fonts = fonts, rects = dplyr::bind_rows(rects),
       texts = dplyr::bind_rows(texts),
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
pb <- P$rects %>% dplyr::filter(fill == "pattern")
.assert(paste0("in the PDF, 224 cells, each in its size, colour and place, and ",
               "the colour bar's gradient"),
        nrow(pc) == 224L && all(abs(pc$w - CELL_W6) < 0.006) &&
          all(abs(pc$h - CELL_H) < 0.006) && identical(pc$fill, .rgb3(r$fill)) &&
          all(abs(pc$x - r$x) < 0.006) && all(abs(pc$y - (H - r$y - r$h)) < 0.006) &&
          nrow(pb) == 1L && P$n_shadings == 1L, n = nrow(pc))
pt <- P$texts
exp_font <- ifelse(tx$face == 3L, "Helvetica-Oblique", "Helvetica")
exp_lab <- gsub(SOFT_HY, "~", tx$label_pdf, fixed = TRUE)
# A label's start: the anchor less its justification, or, standing vertically,
# the anchor itself.
exp_x <- ifelse(tx$rot == 90, tx$x, tx$left)
.assert(paste0("in the PDF every label is live Helvetica text (no image, no ",
               "embedded or Type 3 font), at its size, angle and place, oblique for ",
               "the gene symbols only"),
        nrow(pt) == nrow(tx) && identical(pt$label, exp_lab) &&
          identical(pt$font, exp_font) && all(abs(pt$size - tx$size) < 0.01) &&
          all(abs(pt$angle - tx$rot) < 0.05) && all(abs(pt$y - (H - tx$y)) < 0.006) &&
          all(abs(pt$x - exp_x) < 0.006) && P$n_images == 0L && P$n_type3 == 0L &&
          !P$font_file, n = nrow(pt), max_dev = max(abs(pt$x - exp_x)))
.assert("the PDF's page is the panel's size, in whole points",
        identical(as.integer(P$mediabox[3:4]), as.integer(floor(c(W, H)))), n = 1L,
        detail = paste(P$mediabox, collapse = " "))
img <- png::readPNG(FILES[["png_docs"]])
s_px <- PNG_DPI / 72
px_x <- floor((r$x + 2) * s_px) + 1L
px_y <- floor((r$y + 1.5) * s_px) + 1L
px_col <- grDevices::rgb(round(img[cbind(px_y, px_x, 1L)] * 255),
                         round(img[cbind(px_y, px_x, 2L)] * 255),
                         round(img[cbind(px_y, px_x, 3L)] * 255), maxColorValue = 255)
.assert(paste0("the PNG is the panel at ", PNG_DPI, " dpi, and every cell shows ",
               "its colour exactly"),
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
# 6. The legend draft, Claude Code's, for the author
# =============================================================================
.assert("E49's legend for ED-h gives the cohorts' n, which 6c's draft repeats",
        grepl("TCGA-BRCA (n = 938) and SCAN-B (n = 3,143)", e49$legend_draft$ed_h,
              fixed = TRUE), n = 1L)
LEGEND_DRAFT <- paste0(
  "Unadjusted Spearman correlations of fifteen BCL-2 family genes and the ",
  "signature with OXPHOS (level and priority), four MYC activity scores and ",
  "proliferation, in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143). Each MYC ",
  "activity score is labelled with the percentage of its genes in proliferation ",
  "gene sets. Rows and colours as in Fig. 4g; colours saturate at +/-0.6",
  if (nrow(BEYOND)) paste0(", which ", nrow(BEYOND), " correlations exceed (",
                           paste(sprintf("%s with the %s, %s, %+.2f", BEYOND$readout,
                                         EXPO_LABEL[BEYOND$exposure], BEYOND$cohort,
                                         BEYOND$rho), collapse = "; "), ")") else "",
  ".")
.assert("the legend draft is ASCII, as the script is", !grepl("[^ -~]", LEGEND_DRAFT), n = 1L)

# =============================================================================
# 7. Saved: the object, two tables and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

SPEC <- list(
  ruling = "docs/2026-11-05_E55_E56_ed_fig6cd_ruling.md (dd5004a)",
  source = "E49's raw_correlations (a63ef6c), as drawn in ED-h; E48's names for the measures",
  rows = c(FIFTEEN, SIG_ROW), groups = GROUPS, exposures = EXPO_LABEL,
  cohorts = COHORTS, scale = L_SCALE, palette = DIV_COLS, bar_title = BAR_TITLE,
  cell_pt = c(w = CELL_W6, h = CELL_H), heatmap_pt = c(w = HEAT_W_4G, h = HEAT_H_4G),
  page_pt = LAY$size, type_pt = c(text = PT_TXT, cohort_headers = PT_HEAD),
  geometry_source = "E50 (Fig. 4g), its constants checked line for line")
OUT <- list(
  values = V %>% dplyr::arrange(match(readout, c(FIFTEEN, SIG_ROW)),
                                match(cohort, COHORTS), match(exposure, EXPOSURES)),
  beyond_scale = BEYOND,
  layout = LAY,
  legend_draft = LEGEND_DRAFT,
  assertions = ASSERTIONS,
  figure_files = FIGURE_FILES,
  spec = SPEC,
  run_info = list(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  R = R.version.string,
                  bitmap_type_default = getOption("bitmapType")))
saveRDS(OUT, PATH_E55)
utils::write.csv(OUT$values, file.path(OUT_TABLES, "E55_raw_correlations.csv"), row.names = FALSE)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E55_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), "run_info"),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-13s %s", k, DIGESTS[[k]]))
message("\n   the legend draft (Claude Code's):\n   ", LEGEND_DRAFT)
message("\nE55 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figure.
  library(dplyr)
  x <- readRDS(here::here("results", "e55_raw_correlations_heatmap.rds"))
  x$beyond_scale %>% print()
  x$assertions %>% print(n = Inf)
  cat(x$legend_draft, "\n")
  img <- png::readPNG(here::here("docs", "figures",
                                 "2026-11-05_E55_raw_correlations_heatmap.png"))
  grid::grid.newpage()
  grid::grid.raster(img)
}
