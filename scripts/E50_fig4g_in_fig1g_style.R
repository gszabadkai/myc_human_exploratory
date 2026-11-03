# E50_fig4g_in_fig1g_style.R
# =============================================================================
# FIG. 4G, ED-E AND ED-F REDRAWN IN FIG. 1G'S STYLE, AS A PDF AND AN SVG THAT
# STAY EDITABLE IN ILLUSTRATOR; ED-G AND ED-H RE-COLOURED. NO GREY ANYWHERE.
#
# DISPLAY-ONLY. It fits nothing and computes no new quantity. It reads E49's
# verified saved object, checked against the digests committed before E49's
# run (a63ef6c), and draws what is stored there.
#
# Specification: the author's further amendment to the STOP 3 rulings,
# verbatim in docs/2026-11-03_E50_display_amendment.md, committed alone at
# 6f5bf4d before this script was written. "Item" below means its items.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. An earlier cairo_pdf in the same
# session can move the text in a later cairo PNG (section 6), and the PNGs are
# compared byte for byte with the dry run.
#
# =============================================================================
# WHAT THIS DOES
# =============================================================================
# Item 1, NO GREY. Every cell is drawn in its actual colour, including the
#   signature's OXPHOS coefficient and the difference that contains it. No
#   text reads either. In each legend the grey sentence is replaced by the
#   author's: "The signature was selected for its association with OXPHOS at
#   fixed proliferation (Methods)."
# Item 2, FIG. 1G'S PALETTE. Its source IS in the mouse repo, read with git
#   show at myc_mouse 2ebf646:
#     figures/natmetab/_style.R              div_cols, scale_fill_div(),
#                                            lab_num(), lab_signed()
#     figures/natmetab/fig1g_bcl2_heatmap.R  the panel
#   So the source's anchors are used, not the values measured from the PDF:
#     source     #4A3525  #FAFAFA  #2A8A6D   used
#     measured   #462915  #F9F9F9  #318866   recorded, not used
#   Fig. 1g's scale is ggplot2's scale_fill_gradientn(), which interpolates in
#   CIELAB, not in RGB. The source is followed, and every cell's colour is
#   asserted equal to the colour Fig. 1g's own scale call gives its value.
#   The limits are symmetric and rounded outward to 0.1, as in E49.
# Item 3, 4G, ED-E AND ED-F IN FIG. 1G'S STYLE, with no annotation columns.
#   The author adds the binding code, the signature signs, the section marks
#   and any protein names in Illustrator. data/annotations/bh3_binding_code.csv
#   stays the record behind them; this script does not read it.
# Item 4, OUTPUT. For 4g, ED-e and ED-f:
#   - a PDF from base pdf() with Helvetica. Helvetica is a standard PDF font,
#     referenced rather than embedded, so the text stays live and editable.
#     This is Fig. 1g's save_panel(), which gives the reason cairo_pdf is not
#     used;
#   - an SVG whose text is <text> elements;
#   - a PNG preview at 600 dpi, through the cairo device, which writes every
#     colour exactly. The default quartz device colour-manages the PNG and
#     darkens the palette (#624E3F is written as #4F3D30).
#   ED-g and ED-h are re-coloured only, and written as E49 wrote them, except
#   that their PNGs also go through cairo.
#
# =============================================================================
# HOW FIG. 1G'S DRAWING CODE IS REUSED
# =============================================================================
# Fig. 1g is a ggplot. Here the panel is laid out ONCE, in PostScript points,
# and that one layout is written three ways: grid draws the PDF and the PNG,
# and a short writer emits the SVG. There are two reasons:
#   - the cells must be exactly 26.71 x 8.46 pt, the size Fig. 1g's code
#     draws at 66 x 56 mm (item 3);
#   - the SVG must carry live text. svglite, the device that would write one
#     from a ggplot, is not installed, and base svg() outlines every letter.
# What is taken from Fig. 1g's code, unchanged:
#   - the palette and the scale call: anchors, values, limits, squish;
#   - the value labels: sprintf("%+.1f"), with "0" for +0.0 and -0.0;
#   - the tick labels: lab_signed(), at -L, -L/2, 0, +L/2 and +L, as Fig. 1g's
#     five ticks;
#   - the colour bar:
#     - 2 x 18 mm, with no frame;
#     - ticked in grey30 on both edges;
#     - its title top-aligned with the cells;
#     - placed 1 mm past the wider of the first two group names, as Fig. 1g's
#       key_spacing places its key;
#   - the gap between groups, 0.8 mm (2.27 pt);
#   - the offsets of the labels from the cells:
#     - gene symbols end 1.2 pt to the left;
#     - column labels sit 1.2 pt above;
#     - group names start 2.4 pt to the right, centred on their group;
#   - Helvetica throughout, from base pdf(family = "Helvetica", useDingbats =
#     FALSE).
#     - In the PDF, a hyphen inside a word is written as a soft hyphen,
#       because the PDF encoding prints "-" as a minus sign.
#     - The SVG and the PNG carry a true minus sign (U+2212) and a hyphen.
# The offsets and the bar's geometry were measured on 2026-11-03. Fig. 1g's
# code was rebuilt with dummy values and written to an uncompressed PDF, and
# the positions were read from it.
#
# =============================================================================
# CHOICES THE AMENDMENT LEAVES OPEN, made here and reported
# =============================================================================
#   - No cell borders (item 3). Fig. 1g's tiles carry a 0.64 pt white border.
#   - ONE colour scale for 4g, ED-e and ED-f, by item 2's rule over all three
#     panels' values. It is +/-0.6, so a colour means the same value in each.
#     On its own values ED-e would be +/-0.5.
#   - The gap before the signature row is twice the group gap, 1.6 mm.
#   - The ink of a value is black, or white where white has the higher
#     contrast with the cell.
#     - Fig. 1g needs only black, because its limits are +/-2.
#     - At +/-0.6 the brown end is too dark for black.
#     - Every value's contrast is asserted at least 4.5.
#   - Type follows item 3's two sizes: 6 pt for the score headers, 5 pt for
#     everything else. That includes the group names and the bar's title,
#     which Fig. 1g sets at 6 pt.
#   - The score header sits 2.4 pt above the cohort labels, Fig. 1g's strip
#     margin. Fig. 1g has no second header tier.
#   - The PDF and the SVG have no background rectangle. The PNG is on white.
#   - ED-h is drawn once, unmarked, because there are no dots (item 3).
#   - ED-g and ED-h keep E49's layout, labels and text. Each takes the same
#     ink rule, and limits by the same rule on its own values.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing, refits nothing and reads no expression data.
#   - It does not draw the binding code, the signature signs or the section
#     marks, and it does not read the binding-code table (item 5).
#   - It reads no E49 quantity as evidence.
#   - It writes no text except the legend drafts, in which only the grey
#     sentences and the marked version's sentence change.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages(library(dplyr))
for (pkg in c("xml2", "png", "digest", "scales", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E50 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE50: Fig. 4g, ED-e and ED-f in Fig. 1g's style; ED-g and ED-h ",
        "re-coloured\n", strrep("=", 78))
message("DISPLAY-ONLY. Reads E49's verified object and fits nothing. ",
        "Specification 6f5bf4d.")

PATH_E49 <- file.path(DIR_RESULTS, "e49_signature_and_fig4g.rds")
PATH_E50 <- file.path(DIR_RESULTS, "e50_fig4g_in_fig1g_style.rds")
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-03"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG <- c(fig4g = "E50_bcl2_family_heatmap_level",
         ed_e  = "E50_bcl2_family_heatmap_priority",
         ed_f  = "E50_bcl2_family_heatmap_reference_score",
         ed_g  = "E50_signature_strips",
         ed_h  = "E50_raw_beside_adjusted")

# =============================================================================
# 0. Constants
# =============================================================================
# The elements of E49's object read here, with the digests recorded in
# docs/2026-11-03_e49_prerun_digests.md section 3, committed at a63ef6c before
# the author's run. That run reproduced all 24.
E49_DIGESTS <- c(
  all_readouts_m3  = "7023883b873b17020e800838db591684",
  placement_rule   = "82d62a1e160e2ff889e0aa854872690f",
  raw_correlations = "f9f32a13550b07075b1690cb7d0a6a36",
  displays         = "419f5c39819c1f951129a794948b6294",
  legend_draft     = "6b7639c40dcde9ae894183f76b367829",
  assertions       = "2710612444bdbe6592173817b73c2d86")

COHORTS <- c("TCGA", "SCAN-B")
TERMS   <- c("OX", "MYC", "PROLIF")
TERM_TITLE <- c(OX = "OXPHOS", MYC = "MYC activity", PROLIF = "Proliferation")
EST_4G  <- "clean"      # MYC_UP.V1_UP__MITOSTRIP, the least-entangled score
EST_REF <- "M_a"        # FELSHER__MITOSTRIP, the reference score
INSTR   <- c(gsva = "level", mitopps = "priority")
COEF_COL <- c(OX = "b_OX_m3", MYC = "b_MYC_m3", PROLIF = "b_PR")

# Fig. 1g's four groups and order (E49 declaration 9.1), asserted below to be
# the rows E49 drew.
GROUPS <- list(Effectors = c("BAX", "BAK1", "BOK"),
               Activators = c("BID", "BBC3", "BCL2L11"),
               Sensitisers = c("BAD", "PMAIP1", "BMF", "BIK"),
               `Anti-apoptotic` = c("MCL1", "BCL2L2", "BCL2A1", "BCL2",
                                    "BCL2L1"))
FIFTEEN   <- unlist(GROUPS, use.names = FALSE)
SIG_ROW   <- "S6"
SIG_LABEL <- "Signature"

# The three panels drawn in Fig. 1g's style, and the E49 display each redraws.
PANELS <- tibble::tribble(
  ~panel,  ~estimator, ~instrument, ~e49_display,
  "fig4g", EST_4G,     "gsva",      "fig4g_a",
  "ed_e",  EST_4G,     "mitopps",   "ed_e_a",
  "ed_f",  EST_REF,    "gsva",      "ed_f_a")

# --- Fig. 1g's palette (item 2) ----------------------------------------------
# myc_mouse 2ebf646, figures/natmetab/_style.R, verbatim:
#   div_cols <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")
DIV_COLS <- c(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")
# The author's measurement from the Fig. 1 PDF (item 2). Recorded, not used,
# because the source exists.
MEASURED_COLS <- c(neg = "#462915", zero = "#F9F9F9", pos = "#318866")
E49_GREY <- "#A6A6A6"   # E49's grey for the circular cells; asserted absent

# --- geometry, in PostScript points (1/72 in; grid's "bigpts") ---------------
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
MINUS_U  <- "\u2212"     # the minus sign, in the SVG and the PNG
SOFT_HY  <- "\u00ad"     # the PDF's hyphen (its encoding prints "-" as minus)
HJ <- c(start = 0, middle = 0.5, end = 1)

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

# =============================================================================
# 1. E49's object, checked against its committed digests
# =============================================================================
message("\n1. E49's saved object")
.stop_if(file.exists(PATH_E49), "E49's saved object is missing: source ",
         "scripts/E49_signature_and_fig4g.R first. STOP.")
e49 <- readRDS(PATH_E49)
GOT49 <- vapply(names(E49_DIGESTS), function(k) digest::digest(e49[[k]]),
                character(1))
.assert(paste0("E49's object is the verified one: every element read here ",
               "matches its digest committed before E49's run (a63ef6c)"),
        identical(unname(GOT49), unname(E49_DIGESTS)), n = length(GOT49),
        detail = paste(names(GOT49)[GOT49 != E49_DIGESTS], collapse = ", "))
.assert("E49's 192 assertions all passed",
        nrow(e49$assertions) == 192L && all(e49$assertions$pass),
        n = nrow(e49$assertions))
.assert(paste0("4g is E49's branch: three blocks, all from the least-entangled ",
               "model (placement rule)"),
        identical(e49$placement_rule$fig4g$model, EST_4G) &&
          identical(e49$placement_rule$fig4g$terms, TERMS), n = 1L)
E49_ROWS <- e49$displays$heatmaps$fig4g_a$cells %>%
  dplyr::distinct(readout, y) %>% dplyr::arrange(dplyr::desc(y))
.assert(paste0("the rows are Fig. 1g's groups and order, as E49 drew them ",
               "(declaration 9.1)"),
        identical(E49_ROWS$readout, FIFTEEN), n = length(FIFTEEN))

# =============================================================================
# 2. The values of 4g, ED-e and ED-f, and the shared scale
# =============================================================================
message("\n2. values")
.panel_values <- function(est, inst) {
  X <- e49$all_readouts_m3 %>%
    dplyr::filter(estimator == est, instrument == inst,
                  readout %in% c(FIFTEEN, SIG_ROW))
  .stop_if(nrow(X) == 32L, "E49 has not 32 rows for ", est, " / ", inst)
  dplyr::bind_rows(lapply(TERMS, function(tm) X %>%
    dplyr::transmute(readout, cohort, term = tm,
                     value = .data[[COEF_COL[[tm]]]], void)))
}
VALS <- lapply(stats::setNames(PANELS$panel, PANELS$panel), function(p) {
  r <- PANELS[PANELS$panel == p, ]
  .panel_values(r$estimator, r$instrument)
})
for (p in PANELS$panel) {
  h <- e49$displays$heatmaps[[PANELS$e49_display[PANELS$panel == p]]]
  shown <- dplyr::bind_rows(
    h$cells %>% dplyr::select(readout, cohort, term, est),
    h$strip %>% dplyr::select(readout, cohort, term, est))
  j <- VALS[[p]] %>% dplyr::inner_join(shown, by = c("readout", "cohort", "term"))
  .assert(paste0(p, ": the 96 values are those E49 drew, the signature's ",
                 "grey cells included, and none is void"),
          nrow(VALS[[p]]) == 96L && nrow(j) == 96L &&
            identical(j$value, j$est) && !any(VALS[[p]]$void) &&
            !anyNA(VALS[[p]]$value), n = nrow(j))
}

.lim_rule <- function(v) ceiling(max(abs(v)) * 10 - 1e-9) / 10   # E49's
ALL_V <- unlist(lapply(VALS, `[[`, "value"), use.names = FALSE)
L_HM  <- .lim_rule(ALL_V)
L_OWN <- vapply(VALS, function(v) .lim_rule(v$value), numeric(1))
.assert(paste0("4g, ED-e and ED-f share one scale, symmetric, rounded outward ",
               "to 0.1 and containing every value: +/-0.6, as item 2 states ",
               "for the present values"),
        isTRUE(all.equal(L_HM, 0.6)) && all(abs(ALL_V) <= L_HM),
        n = length(ALL_V), max_dev = max(abs(ALL_V)),
        detail = paste0("own limits: ", paste(names(L_OWN), sprintf("%.1f", L_OWN),
                                              collapse = ", ")))
message(sprintf("   one scale for the three panels: +/-%.1f", L_HM))

# =============================================================================
# 3. Fig. 1g's palette and labels (item 2)
# =============================================================================
message("\n3. palette")
# Fig. 1g's scale_fill_div(), whose colours are those of this call.
.fig1g_scale <- function(L, ...) {
  ggplot2::scale_fill_gradientn(
    colours = unname(DIV_COLS[c("neg", "zero", "pos")]),
    values = scales::rescale(c(-L, 0, L)), limits = c(-L, L),
    oob = scales::squish, ...)
}
# The same mapping as a function of the value, for the layouts below.
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
# Fig. 1g's number formats, verbatim from _style.R.
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
# Fig. 1g's tile label, verbatim: sprintf("%+.1f"), with "0" for a zero.
.val_ascii <- function(v) {
  s <- sprintf("%+.1f", v)
  s[s %in% c("+0.0", "-0.0")] <- "0"
  s
}
.val_uni <- function(v) sub("^-", MINUS_U, .val_ascii(v))
.to_pdf  <- function(u) gsub(MINUS_U, "-", gsub("-", SOFT_HY, u, fixed = TRUE),
                             fixed = TRUE)
# Ink: black, or white where white has the higher WCAG contrast.
.lin1 <- function(u) ifelse(u <= 0.04045, u / 12.92, ((u + 0.055) / 1.055)^2.4)
.lum <- function(h) {
  m <- grDevices::col2rgb(h) / 255
  as.numeric(0.2126 * .lin1(m[1, ]) + 0.7152 * .lin1(m[2, ]) +
               0.0722 * .lin1(m[3, ]))
}
.contrast <- function(a, b) {
  la <- .lum(a)
  lb <- .lum(b)
  (pmax(la, lb) + 0.05) / (pmin(la, lb) + 0.05)
}
.ink <- function(h) ifelse(.contrast(h, "#FFFFFF") > .contrast(h, "#000000"),
                           "#FFFFFF", "#000000")

.assert(paste0("the palette's three anchors are Fig. 1g's own (myc_mouse ",
               "2ebf646, _style.R div_cols), which item 2 puts before the ",
               "measured values, and the scale returns them at -L, 0 and +L"),
        identical(unname(DIV_COLS), c("#4A3525", "#FAFAFA", "#2A8A6D")) &&
          identical(.fill_hex(c(-L_HM, 0, L_HM), L_HM), unname(DIV_COLS)) &&
          identical(.ggplot_fill(c(-L_HM, 0, L_HM), L_HM), unname(DIV_COLS)),
        n = 3L)
.assert(paste0("every cell's colour is the one Fig. 1g's own scale call gives ",
               "its value (ggplot2::scale_fill_gradientn, CIELAB)"),
        identical(.fill_hex(ALL_V, L_HM), .ggplot_fill(ALL_V, L_HM)),
        n = length(ALL_V))
PAL_TABLE <- tibble::tibble(
  anchor = c("negative limit", "zero", "positive limit"),
  source_used = unname(DIV_COLS), measured_not_used = unname(MEASURED_COLS),
  max_channel_difference = vapply(1:3, function(i)
    max(abs(grDevices::col2rgb(DIV_COLS[i]) -
              grDevices::col2rgb(MEASURED_COLS[i]))), numeric(1)))

# =============================================================================
# 4. The layout of 4g, ED-e and ED-f (item 3), in points from the top left
# =============================================================================
message("\n4. layouts")
# Text widths from Helvetica's AFM metrics, through the PDF device, as Fig.
# 1g's style file measures them. Strings are the PDF's.
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

.layout <- function(V, L, panel) {
  rows <- tibble::tibble(
    readout = c(FIFTEEN, SIG_ROW),
    group = c(rep(names(GROUPS), lengths(GROUPS)), SIG_LABEL)) %>%
    dplyr::mutate(k = dplyr::row_number(),
                  g = match(group, c(names(GROUPS), SIG_LABEL)),
                  label = ifelse(readout == SIG_ROW, SIG_LABEL, readout),
                  face = ifelse(readout == SIG_ROW, 1L, 3L),
                  role = ifelse(readout == SIG_ROW, "signature", "gene"))
  lab_w <- max(.text_w(.to_pdf(rows$label), PT_TXT, rows$face))
  X0 <- MARGIN[["left"]] + lab_w + LAB_GAP
  Y_SCORE <- MARGIN[["top"]] + CAP * PT_HEAD
  Y_COH <- Y_SCORE + DESC * PT_HEAD + HEAD_GAP + CAP * PT_TXT
  Y0 <- Y_COH + COH_GAP
  rows <- rows %>%
    dplyr::mutate(y = Y0 + (k - 1) * CELL_H + pmin(g - 1L, 3L) * GAP +
                    (g == 5L) * GAP_SIG)
  cols <- tidyr::expand_grid(term = TERMS, cohort = COHORTS) %>%
    dplyr::mutate(b = match(term, TERMS), c = match(cohort, COHORTS),
                  x = X0 + (b - 1) * (2 * CELL_W + GAP) + (c - 1) * CELL_W)
  X1 <- max(cols$x) + CELL_W
  rects <- V %>%
    dplyr::inner_join(rows %>% dplyr::select(readout, y), by = "readout") %>%
    dplyr::inner_join(cols %>% dplyr::select(term, cohort, x),
                      by = c("term", "cohort")) %>%
    dplyr::arrange(match(readout, rows$readout), match(term, TERMS),
                   match(cohort, COHORTS)) %>%
    dplyr::mutate(w = CELL_W, h = CELL_H, fill = .fill_hex(value, L),
                  ink = .ink(fill), contrast = .contrast(fill, ink))
  .stop_if(nrow(rects) == 96L, panel, ": a cell is missing")
  grp <- rows %>% dplyr::filter(g <= 4L) %>% dplyr::group_by(group) %>%
    dplyr::summarise(mid = (min(y) + max(y) + CELL_H) / 2, .groups = "drop") %>%
    dplyr::arrange(match(group, names(GROUPS)))
  blocks <- cols %>% dplyr::group_by(term) %>%
    dplyr::summarise(x = min(x) + CELL_W, .groups = "drop") %>%
    dplyr::arrange(match(term, TERMS))
  # The colour bar, placed as Fig. 1g places its key.
  XG <- X1 + GRP_GAP
  XB <- XG + max(.text_w(c("Effectors", "Activators"), PT_TXT, 1L)) + KEY_GAP
  Y_TITLE <- Y0 + CAP * PT_TXT
  YB_TOP  <- Y_TITLE + DESC * PT_TXT + TITLE_GAP
  YB_BOT  <- YB_TOP + BAR_H
  TICKS   <- c(-L, -L / 2, 0, L / 2, L)
  y_tick  <- YB_BOT - (TICKS + L) / (2 * L) * BAR_H
  texts <- dplyr::bind_rows(
    .txt("value", rects$x + CELL_W / 2, rects$y + CELL_H / 2 + CAP * PT_TXT / 2,
         "middle", PT_TXT, 1L, .val_uni(rects$value), rects$ink,
         paste(rects$readout, rects$term, rects$cohort)),
    .txt(rows$role, X0 - LAB_GAP, rows$y + CELL_H / 2 + CAP * PT_TXT / 2, "end",
         PT_TXT, rows$face, rows$label, key = rows$readout),
    .txt("group", XG, grp$mid + CAP * PT_TXT / 2, "start", PT_TXT, 1L,
         grp$group, key = grp$group),
    .txt("cohort", cols$x + CELL_W / 2, Y_COH, "middle", PT_TXT, 1L,
         cols$cohort, key = paste(cols$term, cols$cohort)),
    .txt("score", blocks$x, Y_SCORE, "middle", PT_HEAD, 1L,
         unname(TERM_TITLE[blocks$term]), key = blocks$term),
    .txt("bar_title", XB, Y_TITLE, "start", PT_TXT, 1L, "SD per SD"),
    .txt("tick", XB + BAR_W + TICKLAB_GAP, y_tick + CAP * PT_TXT / 2, "start",
         PT_TXT, 1L, sub("^-", MINUS_U, lab_signed(TICKS)),
         key = sprintf("%+.2f", TICKS))) %>%
    dplyr::mutate(label_pdf = .to_pdf(label_uni),
                  width = .text_w(label_pdf, size, face),
                  left = x - unname(HJ[anchor]) * width,
                  right = left + width,
                  top = y - CAP * size, bottom = y + DESC * size)
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
  list(panel = panel, L = L, rows = rows, cols = cols, rects = rects,
       texts = texts, bar = bar, stops = stops, ticks = ticks,
       size = c(W = W, H = H), x0 = X0, x1 = X1, y0 = Y0)
}

.check_layout <- function(lay, label) {
  r <- lay$rects
  tx <- lay$texts
  rows <- lay$rows
  .assert(paste0(label, ": 96 cells, each 26.71 x 8.46 pt, with no border"),
          nrow(r) == 96L && all(r$w == CELL_W) && all(r$h == CELL_H),
          n = nrow(r))
  ry <- rows %>% dplyr::arrange(k) %>%
    dplyr::mutate(gap = dplyr::lead(y) - (y + CELL_H), g_next = dplyr::lead(g))
  within  <- ry$gap[!is.na(ry$gap) & ry$g == ry$g_next]
  between <- ry$gap[!is.na(ry$gap) & ry$g != ry$g_next & ry$g_next <= 4L]
  before  <- ry$gap[!is.na(ry$gap) & ry$g_next == 5L]
  .assert(paste0(label, ": no gap within a group; about 2.3 pt (0.8 mm) ",
                 "between groups; a larger gap before the signature"),
          length(within) == 11L && all(abs(within) < 1e-9) &&
            length(between) == 3L && all(abs(between - GAP) < 1e-9) &&
            all(abs(between - 2.3) < 0.05) && length(before) == 1L &&
            before > GAP + 1e-9, n = 15L,
          detail = sprintf("between groups %.3f pt; before the signature %.3f pt",
                           GAP, before))
  cx <- lay$cols %>% dplyr::arrange(b, c) %>%
    dplyr::mutate(gap = dplyr::lead(x) - (x + CELL_W), b_next = dplyr::lead(b))
  in_block <- cx$gap[!is.na(cx$gap) & cx$b == cx$b_next]
  across   <- cx$gap[!is.na(cx$gap) & cx$b != cx$b_next]
  .assert(paste0(label, ": the two cohorts abut within a block; about 2.3 pt ",
                 "(0.8 mm) between the three blocks"),
          length(in_block) == 3L && all(abs(in_block) < 1e-9) &&
            length(across) == 2L && all(abs(across - GAP) < 1e-9), n = 5L)
  .assert(paste0(label, ": rows in Fig. 1g's groups and order with the ",
                 "signature last; columns OXPHOS, MYC activity, ",
                 "proliferation, each TCGA then SCAN-B"),
          identical(rows$readout[order(rows$y)], c(FIFTEEN, SIG_ROW)) &&
            identical(cx$term, rep(TERMS, each = 2)) &&
            identical(cx$cohort, rep(COHORTS, 3)) &&
            identical(tx$label_uni[tx$role == "score"], unname(TERM_TITLE)) &&
            identical(tx$label_uni[tx$role == "cohort"], rep(COHORTS, 3)) &&
            identical(tx$label_uni[tx$role == "group"], names(GROUPS)) &&
            identical(tx$label_uni[tx$role == "signature"], SIG_LABEL),
          n = 16L)
  exp_size <- ifelse(tx$role == "score", PT_HEAD, PT_TXT)
  .assert(paste0(label, ": Helvetica 5 pt for values and labels, 6 pt for the ",
                 "score headers; oblique for the gene symbols and nothing else"),
          all(tx$size == exp_size) && all((tx$face == 3L) == (tx$role == "gene")),
          n = nrow(tx))
  vt <- tx %>% dplyr::filter(role == "value")
  m <- match(vt$key, paste(r$readout, r$term, r$cohort))
  printed <- as.numeric(sub(MINUS_U, "-", vt$label_uni, fixed = TRUE))
  dev <- abs(printed - r$value[m])
  .assert(paste0(label, ": every cell prints its stored value, signed, to one ",
                 "decimal, and \"0\" when it rounds to zero"),
          nrow(vt) == 96L && !anyNA(m) && !any(duplicated(m)) &&
            identical(vt$label_uni, .val_uni(r$value[m])) &&
            all(dev <= 0.05 + 1e-12) &&
            all(grepl(paste0("^([+", MINUS_U, "][0-9]+\\.[0-9]|0)$"),
                      vt$label_uni)) &&
            all((vt$label_uni == "0") ==
                  (sprintf("%.1f", abs(r$value[m])) == "0.0")),
          n = nrow(vt), max_dev = max(dev))
  .assert(paste0(label, ": each value is centred in its cell and reads at a ",
                 "contrast of 4.5 or more"),
          all(abs(vt$x - (r$x[m] + CELL_W / 2)) < 1e-9) &&
            all(abs(vt$y - (r$y[m] + CELL_H / 2 + CAP * PT_TXT / 2)) < 1e-9) &&
            all(r$contrast >= 4.5) && all(vt$ink == r$ink[m]),
          n = nrow(vt), detail = sprintf("lowest contrast %.2f", min(r$contrast)))
  inside <- vt$left >= r$x[m] + 1 & vt$right <= r$x[m] + CELL_W - 1 &
    vt$top >= r$y[m] & vt$bottom <= r$y[m] + CELL_H
  ct <- tx %>% dplyr::filter(role %in% c("cohort", "score"))
  .assert(paste0(label, ": every value fits inside its cell; each cohort label ",
                 "fits its column and each score header its block"),
          all(inside) &&
            all(ct$width[ct$role == "cohort"] <= CELL_W - 2) &&
            all(ct$width[ct$role == "score"] <= 2 * CELL_W - 2),
          n = nrow(vt) + nrow(ct))
  # Text boxes run from the cap height to the descender.
  ov <- outer(seq_len(nrow(tx)), seq_len(nrow(tx)), function(i, j)
    tx$left[i] < tx$right[j] & tx$left[j] < tx$right[i] &
      tx$top[i] < tx$bottom[j] & tx$top[j] < tx$bottom[i])
  n_ov <- sum(ov[upper.tri(ov)])
  bb <- lay$bar
  hit_bar <- tx$left < bb$x + bb$w & bb$x < tx$right &
    tx$top < bb$y + bb$h & bb$y < tx$bottom
  lab <- tx[tx$role != "value", ]
  cc <- outer(seq_len(nrow(lab)), seq_len(nrow(r)), function(i, k)
    lab$left[i] < r$x[k] + CELL_W & r$x[k] < lab$right[i] &
      lab$top[i] < r$y[k] + CELL_H & r$y[k] < lab$bottom[i])
  cell_hits <- sum(cc)
  W <- lay$size[["W"]]
  H <- lay$size[["H"]]
  .assert(paste0(label, ": no two labels overlap; no label touches a cell or ",
                 "the colour bar; everything lies inside the panel"),
          n_ov == 0L && !any(hit_bar) && cell_hits == 0L &&
            all(tx$left >= 0) && all(tx$right <= W) && all(tx$top >= 0) &&
            all(tx$bottom <= H) && lay$bar$x + lay$bar$w <= W &&
            lay$bar$y + lay$bar$h <= H,
          n = nrow(tx), detail = sprintf("%d overlaps, %d label-cell contacts",
                                         n_ov, cell_hits))
  tk <- tx %>% dplyr::filter(role == "tick")
  .assert(paste0(label, ": a vertical colour bar on the right, 2 x 18 mm, ",
                 "titled \"SD per SD\", ticked at -L, -L/2, 0, +L/2 and +L ",
                 "with signed labels"),
          lay$bar$x > lay$x1 && abs(lay$bar$w - 2 * MM) < 1e-9 &&
            abs(lay$bar$h - 18 * MM) < 1e-9 &&
            identical(tx$label_uni[tx$role == "bar_title"], "SD per SD") &&
            identical(tk$label_uni, sub("^-", MINUS_U, lab_signed(
              c(-lay$L, -lay$L / 2, 0, lay$L / 2, lay$L)))) &&
            identical(tk$label_uni[c(1, 3, 5)],
                      c(paste0(MINUS_U, lab_num(lay$L)), "0",
                        paste0("+", lab_num(lay$L)))),
          n = nrow(tk))
  st <- lay$stops
  .assert(paste0(label, ": the bar runs from the negative limit at the bottom ",
                 "to the positive limit at the top, through the palette, with ",
                 "its three anchors at -L, 0 and +L"),
          identical(st$colour[c(1L, (N_STOPS + 1L) / 2L, N_STOPS)],
                    unname(DIV_COLS)) &&
            identical(st$colour, .fill_hex(st$value, lay$L)),
          n = nrow(st))
  invisible(TRUE)
}

LAY <- list()
for (p in PANELS$panel) {
  LAY[[p]] <- .layout(VALS[[p]], L_HM, p)
  .check_layout(LAY[[p]], p)
}
.assert("4g, ED-e and ED-f share one layout: the same size and positions",
        identical(LAY$fig4g$size, LAY$ed_e$size) &&
          identical(LAY$fig4g$size, LAY$ed_f$size) &&
          identical(LAY$fig4g$rects[c("x", "y")], LAY$ed_e$rects[c("x", "y")]) &&
          identical(LAY$fig4g$rects[c("x", "y")], LAY$ed_f$rects[c("x", "y")]),
        n = 3L, detail = sprintf("%.1f x %.1f mm", LAY$fig4g$size[["W"]] / MM,
                                 LAY$fig4g$size[["H"]] / MM))

# =============================================================================
# 5. Writing 4g, ED-e and ED-f (item 4), and reading each file back
# =============================================================================
message("\n5. writing PDF, SVG and PNG, then reading them back")
.ensure_dir(DIR_FIGURES)
.ensure_dir(DIR_DOC_FIGURES)
.ensure_dir(DIR_TABLES)

# grid, for the PDF and the PNG. Coordinates in big points from the bottom.
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
                    hjust = unname(HJ[tx$anchor[i]]), vjust = 0,
                    gp = grid::gpar(fontfamily = FONT, fontsize = tx$size[i],
                                    fontface = tx$face[i], col = tx$ink[i]))
  }
  invisible(TRUE)
}
.write_pdf <- function(lay, file, title) {
  grDevices::pdf(file, width = lay$size[["W"]] / 72,
                 height = lay$size[["H"]] / 72, family = FONT,
                 useDingbats = FALSE, compress = FALSE, title = title)
  on.exit(grDevices::dev.off(), add = TRUE)
  .draw(lay, "pdf")
}
# The PNG through cairo, which writes each colour exactly. The default quartz
# device colour-manages it: #624E3F is written as #4F3D30 (checked 2026-11-03).
.write_png <- function(lay, file) {
  grDevices::png(file, width = lay$size[["W"]] / 72,
                 height = lay$size[["H"]] / 72, units = "in", res = PNG_DPI,
                 bg = "white", type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  .draw(lay, "uni")
}
# The SVG: user units are points, text is <text>, groups are named.
.f3 <- function(v) formatC(v, format = "f", digits = 3)
.esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}
.svg_text <- function(d) {
  if (!nrow(d)) return(character())
  sprintf(paste0('<text x="%s" y="%s" text-anchor="%s" font-family="%s"%s ',
                 'font-size="%s" fill="%s">%s</text>'),
          .f3(d$x), .f3(d$y), d$anchor,
          ifelse(d$face == 3L, "Helvetica-Oblique, Helvetica", "Helvetica"),
          ifelse(d$face == 3L, ' font-style="italic"', ""),
          format(d$size), d$ink, .esc(d$label_uni))
}
SVG_GROUPS <- list(values = "value", `row-labels` = c("gene", "signature"),
                   `group-names` = "group", `column-headers` = c("score", "cohort"))
.write_svg <- function(lay, file) {
  W <- lay$size[["W"]]
  H <- lay$size[["H"]]
  r <- lay$rects
  b <- lay$bar
  tk <- lay$ticks
  tx <- lay$texts
  st <- lay$stops
  out <- c(
    '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
    sprintf(paste0('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ',
                   'width="%spt" height="%spt" viewBox="0 0 %s %s">'),
            .f3(W), .f3(H), .f3(W), .f3(H)),
    '<defs>',
    '<linearGradient id="colour-bar-gradient" x1="0" y1="1" x2="0" y2="0">',
    sprintf('<stop offset="%s" stop-color="%s"/>',
            formatC(st$offset, format = "f", digits = 4), st$colour),
    '</linearGradient>',
    '</defs>',
    '<g id="cells">',
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s"/>',
            .f3(r$x), .f3(r$y), .f3(r$w), .f3(r$h), r$fill),
    '</g>')
  for (g in names(SVG_GROUPS)) {
    out <- c(out, sprintf('<g id="%s">', g),
             .svg_text(tx[tx$role %in% SVG_GROUPS[[g]], ]), '</g>')
  }
  out <- c(out, '<g id="colour-bar">',
           sprintf(paste0('<rect x="%s" y="%s" width="%s" height="%s" ',
                          'fill="url(#colour-bar-gradient)"/>'),
                   .f3(b$x), .f3(b$y), .f3(b$w), .f3(b$h)),
           sprintf(paste0('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" ',
                          'stroke-width="%s"/>'),
                   .f3(tk$x0), .f3(tk$y0), .f3(tk$x1), .f3(tk$y1), TICK_COL,
                   format(TICK_LWD)),
           .svg_text(tx[tx$role %in% c("bar_title", "tick"), ]), '</g>',
           '</svg>')
  con <- file(file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(paste0(paste(out, collapse = "\n"), "\n"))), con)
}

# --- reading the files back --------------------------------------------------
.read_svg <- function(file) {
  doc <- xml2::read_xml(file)
  xml2::xml_ns_strip(doc)
  rc <- xml2::xml_find_all(doc, "//g[@id='cells']/rect")
  tx <- xml2::xml_find_all(doc, "//text")
  sp <- xml2::xml_find_all(doc, "//linearGradient/stop")
  num <- function(n, a) as.numeric(xml2::xml_attr(n, a))
  list(
    rects = tibble::tibble(x = num(rc, "x"), y = num(rc, "y"),
                           w = num(rc, "width"), h = num(rc, "height"),
                           fill = xml2::xml_attr(rc, "fill")),
    texts = tibble::tibble(
      group = xml2::xml_attr(xml2::xml_find_first(tx, "parent::*"), "id"),
      x = num(tx, "x"), y = num(tx, "y"),
      anchor = xml2::xml_attr(tx, "text-anchor"),
      family = xml2::xml_attr(tx, "font-family"),
      style = xml2::xml_attr(tx, "font-style"),
      size = num(tx, "font-size"), fill = xml2::xml_attr(tx, "fill"),
      label = xml2::xml_text(tx)),
    stops = tibble::tibble(offset = num(sp, "offset"),
                           colour = xml2::xml_attr(sp, "stop-color")),
    n_paths = length(xml2::xml_find_all(doc, "//path")),
    width = xml2::xml_attr(doc, "width"), height = xml2::xml_attr(doc, "height"))
}
# The PDF is written uncompressed, so its drawing operators are plain text.
# The soft hyphen (byte 0xAD) is read as "~"; other bytes above 127 as "?".
.read_pdf <- function(file) {
  b <- readBin(file, "raw", file.info(file)$size)
  b[b == as.raw(0xAD)] <- as.raw(0x7E)
  b[as.integer(b) > 127L] <- as.raw(0x3F)
  ln <- trimws(strsplit(rawToChar(b), "\n", fixed = TRUE)[[1]])
  fm <- regmatches(ln, regexec("/Name /F([0-9]+) /BaseFont /([A-Za-z-]+)", ln))
  fm <- do.call(rbind, fm[lengths(fm) == 3L])
  fonts <- stats::setNames(fm[, 3], paste0("F", fm[, 2]))
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
      "^/F([0-9]+) 1 Tf ([0-9.]+) [0-9.-]+ [0-9.-]+ [0-9.]+ ([0-9.-]+) ",
      "([0-9.-]+) Tm (.*) T[jJ]$"), l))[[1]]
    if (length(m) == 6L) {
      parts <- regmatches(m[6], gregexpr("\\(([^)]*)\\)", m[6]))[[1]]
      texts[[length(texts) + 1L]] <- tibble::tibble(
        font = unname(fonts[paste0("F", m[2])]), size = as.numeric(m[3]),
        x = as.numeric(m[4]), y = as.numeric(m[5]),
        label = paste(substr(parts, 2L, nchar(parts) - 1L), collapse = ""),
        fill = fill)
    }
  }
  list(fonts = fonts, rects = dplyr::bind_rows(rects),
       texts = dplyr::bind_rows(texts),
       n_images = sum(grepl("/Subtype /Image", ln, fixed = TRUE)),
       n_type3 = sum(grepl("/Subtype /Type3", ln, fixed = TRUE)),
       n_shadings = sum(grepl("/ShadingType 2", ln, fixed = TRUE)))
}
.png_dims <- function(file) {
  b <- as.integer(readBin(file, "raw", 24L))
  c(w = sum(b[17:20] * 256^(3:0)), h = sum(b[21:24] * 256^(3:0)))
}
.png_hex <- function(img) {
  grDevices::rgb(round(img[, , 1] * 255), round(img[, , 2] * 255),
                 round(img[, , 3] * 255), maxColorValue = 255)
}
# Every colour drawn must appear in the PNG, exactly, on at least 100 pixels.
.png_has_colours <- function(file, cols) {
  px <- table(.png_hex(png::readPNG(file)))
  got <- as.integer(px[toupper(cols)])
  all(!is.na(got) & got >= 100L)
}
.rgb3 <- function(h) {
  m <- grDevices::col2rgb(h) / 255
  sprintf("%.3f %.3f %.3f", m[1, ], m[2, ], m[3, ])
}

.check_files <- function(lay, f_pdf, f_svg, f_png, label) {
  W <- lay$size[["W"]]
  H <- lay$size[["H"]]
  r <- lay$rects
  tx <- lay$texts
  # SVG
  s <- .read_svg(f_svg)
  sx <- s$texts
  ord <- unlist(lapply(c(SVG_GROUPS, list(`colour-bar` = c("bar_title", "tick"))),
                       function(rl) which(tx$role %in% rl)), use.names = FALSE)
  tl <- tx[ord, ]
  .assert(paste0(label, " (SVG): 96 cells, each 26.71 x 8.46 pt, in its colour ",
                 "and place; no path, so no outlined text"),
          nrow(s$rects) == 96L && all(abs(s$rects$w - CELL_W) < 5e-4) &&
            all(abs(s$rects$h - CELL_H) < 5e-4) &&
            identical(s$rects$fill, r$fill) &&
            all(abs(s$rects$x - r$x) < 5e-4) && all(abs(s$rects$y - r$y) < 5e-4) &&
            s$n_paths == 0L &&
            identical(s$width, paste0(.f3(W), "pt")) &&
            identical(s$height, paste0(.f3(H), "pt")),
          n = nrow(s$rects))
  .assert(paste0(label, " (SVG): every label is live text, in Helvetica at its ",
                 "size, oblique for the gene symbols only, in its place"),
          nrow(sx) == nrow(tl) && identical(sx$label, tl$label_uni) &&
            all(sx$size == tl$size) && identical(sx$anchor, tl$anchor) &&
            identical(sx$family, ifelse(tl$face == 3L,
                                        "Helvetica-Oblique, Helvetica",
                                        "Helvetica")) &&
            identical(!is.na(sx$style) & sx$style == "italic", tl$face == 3L) &&
            identical(sx$fill, tl$ink) &&
            all(abs(sx$x - tl$x) < 5e-4) && all(abs(sx$y - tl$y) < 5e-4),
          n = nrow(sx))
  sv <- sx %>% dplyr::filter(group == "values")
  pr <- as.numeric(sub(MINUS_U, "-", sv$label, fixed = TRUE))
  .assert(paste0(label, " (SVG): every printed value equals its stored ",
                 "value, rounded to one decimal"),
          nrow(sv) == 96L && identical(sv$label, .val_uni(r$value)) &&
            all(abs(pr - r$value) <= 0.05 + 1e-12), n = nrow(sv),
          max_dev = max(abs(pr - r$value)))
  .assert(paste0(label, " (SVG): the bar's gradient is the palette, its ",
                 "anchors Fig. 1g's"),
          nrow(s$stops) == N_STOPS && identical(s$stops$colour, lay$stops$colour) &&
            identical(s$stops$colour[c(1L, (N_STOPS + 1L) / 2L, N_STOPS)],
                      unname(DIV_COLS)), n = nrow(s$stops))
  # PDF
  p <- .read_pdf(f_pdf)
  pc <- p$rects %>% dplyr::filter(fill != "pattern")
  pb <- p$rects %>% dplyr::filter(fill == "pattern")
  .assert(paste0(label, " (PDF): 96 cells, each 26.71 x 8.46 pt, in its colour ",
                 "and place, and the colour bar"),
          nrow(pc) == 96L && all(abs(pc$w - CELL_W) < 0.006) &&
            all(abs(pc$h - CELL_H) < 0.006) && identical(pc$fill, .rgb3(r$fill)) &&
            all(abs(pc$x - r$x) < 0.006) &&
            all(abs(pc$y - (H - r$y - r$h)) < 0.006) &&
            nrow(pb) == 1L && p$n_shadings == 1L, n = nrow(pc))
  pt <- p$texts
  tl <- tx
  exp_font <- ifelse(tl$face == 3L, "Helvetica-Oblique", "Helvetica")
  exp_lab <- gsub(SOFT_HY, "~", tl$label_pdf, fixed = TRUE)
  .assert(paste0(label, " (PDF): every label is live Helvetica text (no image, ",
                 "no Type 3 font), at its size, oblique for the gene symbols ",
                 "only, in its place and ink"),
          nrow(pt) == nrow(tl) && identical(pt$label, exp_lab) &&
            identical(pt$font, exp_font) && all(pt$size == tl$size) &&
            all(abs(pt$y - (H - tl$y)) < 0.006) &&
            all(abs(pt$x - tl$left) < 0.006) &&
            identical(pt$fill, .rgb3(tl$ink)) &&
            all(pt$font %in% c("Helvetica", "Helvetica-Oblique")) &&
            p$n_images == 0L && p$n_type3 == 0L,
          n = nrow(pt), max_dev = max(abs(pt$x - tl$left)))
  pv <- pt[tl$role == "value", ]
  prp <- as.numeric(pv$label)
  .assert(paste0(label, " (PDF): every printed value equals its stored value, ",
                 "rounded to one decimal"),
          nrow(pv) == 96L && identical(pv$label, .val_ascii(r$value)) &&
            all(abs(prp - r$value) <= 0.05 + 1e-12), n = nrow(pv),
          max_dev = max(abs(prp - r$value)))
  d <- .png_dims(f_png)
  img <- png::readPNG(f_png)
  s_px <- PNG_DPI / 72
  # One pixel per cell, 2 pt in from its left edge and 1.5 pt down from its
  # top, clear of the printed value and of the cell's anti-aliased edge.
  px_x <- floor((r$x + 2) * s_px) + 1L
  px_y <- floor((r$y + 1.5) * s_px) + 1L
  px_col <- grDevices::rgb(round(img[cbind(px_y, px_x, 1L)] * 255),
                           round(img[cbind(px_y, px_x, 2L)] * 255),
                           round(img[cbind(px_y, px_x, 3L)] * 255),
                           maxColorValue = 255)
  .assert(paste0(label, " (PNG): the preview is the panel at ", PNG_DPI,
                 " dpi, and every cell shows its colour exactly"),
          all(abs(d - c(W, H) / 72 * PNG_DPI) <= 1) &&
            identical(px_col, r$fill), n = nrow(r),
          detail = sprintf("%d x %d px", d[["w"]], d[["h"]]))
  list(svg = list(n_rects = nrow(s$rects), n_texts = nrow(sx),
                  stops = nrow(s$stops)),
       pdf = list(n_cells = nrow(pc), n_texts = nrow(pt), fonts = sort(unique(pt$font)),
                  max_x_dev = max(abs(pt$x - tl$left))),
       png = d)
}

.fig_base <- function(dir, stem) {
  file.path(dir, if (identical(dir, DIR_DOC_FIGURES))
    paste0(FIG_DATE, "_", stem) else stem)
}
.same_in_both <- function(stem, ext) {
  identical(unname(tools::md5sum(paste0(.fig_base(DIR_FIGURES, stem), ext))),
            unname(tools::md5sum(paste0(.fig_base(DIR_DOC_FIGURES, stem), ext))))
}
OUT_CHECKS <- list()
for (p in PANELS$panel) {
  stem <- FIG[[p]]
  for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
    base <- .fig_base(dir, stem)
    .write_pdf(LAY[[p]], paste0(base, ".pdf"), stem)
    .write_svg(LAY[[p]], paste0(base, ".svg"))
    .write_png(LAY[[p]], paste0(base, ".png"))
    OUT_CHECKS[[paste0(p, if (identical(dir, DIR_DOC_FIGURES)) "_docs" else "")]] <-
      .check_files(LAY[[p]], paste0(base, ".pdf"), paste0(base, ".svg"),
                   paste0(base, ".png"), paste0(p, if (identical(dir, DIR_DOC_FIGURES))
                     " [docs/figures]" else " [outputs/figures]"))
  }
}
.assert(paste0("4g, ED-e and ED-f: each SVG and PNG is identical in ",
               "outputs/figures and docs/figures"),
        all(vapply(PANELS$panel, function(p) .same_in_both(FIG[[p]], ".svg") &&
                     .same_in_both(FIG[[p]], ".png"), logical(1))), n = 6L)
message("   4g, ED-e and ED-f: PDF, SVG and PNG, to outputs/figures and ",
        "docs/figures (tracked)")

# =============================================================================
# 6. ED-g and ED-h, re-coloured only (item 4)
# =============================================================================
# E49's display helpers and layout constants, copied unchanged from
# scripts/E49_signature_and_fig4g.R (section 8), so that these two panels are
# drawn as E49 drew them except for the colour, the grey and the dots.
message("\n6. ED-g and ED-h, re-coloured")
BASE_PT <- 6.5
MIN_PT  <- 5.5
TXT_MM  <- BASE_PT / ggplot2::.pt
SUB_MM  <- MIN_PT / ggplot2::.pt
W_MM    <- 178
W_PAGE  <- W_MM / 25.4
SECT    <- "\u00a7"
BG      <- 3.0
SIG     <- c("S6", "S5", "S7")
SIG_LAB <- c(S6 = "Signature (S6)", S5 = "S5, without BID",
             S7 = "S7, with BMF")
SECTION_GENES <- c("BAD", "BID")
PROTEIN <- c(BBC3 = "PUMA", BCL2L11 = "BIM", PMAIP1 = "NOXA",
             BCL2L2 = "BCL-w", BCL2A1 = "A1", BCL2L1 = "BCL-xL")
TERM_TITLE49 <- c(TERM_TITLE, DIFF = "Proliferation minus OXPHOS")
EXPO <- c(OX_gsva = "OXPHOS level", OX_mitopps = "OXPHOS priority",
          MYC_clean = "MYC_UP.V1_UP", MYC_M_a = "Felsher, reference",
          MYC_entangled = "Hallmark MYC targets", MYC_regulon = "CollecTRI regulon",
          PROLIF = "Proliferation")

.sym_breaks <- function(L) {
  st <- if (L <= 0.6) 0.2 else if (L <= 1.2) 0.4 else 0.5
  k <- floor(L / st + 1e-9)
  round(seq(-k, k) * st, 10)
}
.cols49 <- function(terms, x0, width) {
  nb <- length(terms)
  cw <- (width - (nb - 1) * BG) / (2 * nb)
  tidyr::expand_grid(term = terms, cohort = COHORTS) %>%
    dplyr::mutate(b = match(term, terms), c = match(cohort, COHORTS),
                  x = x0 + (b - 1) * (2 * cw + BG) + (c - 0.5) * cw, cw = cw)
}
.cell_txt <- function(e, lo, hi) sprintf("%.2f\n(%.2f, %.2f)", e, lo, hi)
.gene_lab <- function(g, sect = TRUE) {
  base <- paste0('italic("', g, '")')
  if (g %in% names(PROTEIN)) base <- paste0(base, '~"(', PROTEIN[[g]], ')"')
  if (sect && g %in% SECTION_GENES) base <- paste0(base, '~"', SECT, '"')
  base
}
.row_lab18 <- function(rd) if (rd %in% SIG) paste0('"', SIG_LAB[[rd]], '"') else
  .gene_lab(rd)
.theme_hm <- function(legend_inside = c(0.985, 0.012)) {
  ggplot2::theme_void(base_size = BASE_PT) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = legend_inside,
      legend.justification = c(1, 0),
      legend.text = ggplot2::element_text(size = MIN_PT),
      legend.title = ggplot2::element_text(size = MIN_PT),
      plot.title = ggplot2::element_text(size = BASE_PT),
      plot.subtitle = ggplot2::element_text(size = MIN_PT),
      plot.caption = ggplot2::element_text(size = MIN_PT),
      strip.text = ggplot2::element_text(size = BASE_PT),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(0, 0, 0, 0))
}
.tag <- function(d, tag) dplyr::mutate(d, .layer = tag)
.layer_of <- function(p, tag) {
  k <- which(vapply(p$layers, function(l) is.data.frame(l$data) &&
                      ".layer" %in% names(l$data) &&
                      identical(unique(as.character(l$data$.layer)), tag),
                    logical(1)))
  .stop_if(length(k) == 1L, "layer not unique: ", tag)
  k
}
.text_sizes <- function(p) {
  th  <- p$theme
  els <- c("axis.text.x", "axis.title.x", "legend.text", "legend.title",
           "plot.title", "plot.subtitle", "plot.caption", "strip.text.x",
           "strip.text.y", "axis.text.y", "axis.title.y")
  es <- vapply(els, function(e) {
    z <- tryCatch(ggplot2::calc_element(e, th), error = function(err) NULL)
    if (is.null(z) || inherits(z, "element_blank") || is.null(z$size)) NA_real_
    else as.numeric(z$size)
  }, numeric(1))
  ls <- unlist(lapply(p$layers, function(l) {
    if (inherits(l$geom, c("GeomText", "GeomLabel"))) {
      s <- l$aes_params$size
      if (is.null(s)) s <- l$geom$default_aes$size
      as.numeric(s) * ggplot2::.pt
    } else NULL
  }))
  c(es, ls)
}
BAN_RX <- "(^|[^A-Za-z0-9_])(m1|m3)([^A-Za-z0-9_]|$)|M_a|beta|Delta"
.panel_strings <- function(p) {
  labs <- tryCatch(ggplot2::get_labs(p), error = function(e) list())
  out <- unlist(labs[intersect(names(labs), c("title", "subtitle", "caption",
                                              "x", "y", "fill"))])
  for (l in p$layers) {
    if (!is.null(l$aes_params$label)) out <- c(out, as.character(l$aes_params$label))
    m <- l$mapping$label
    if (!is.null(m)) {
      dd <- if (is.data.frame(l$data)) l$data else p$data
      v <- tryCatch(rlang::eval_tidy(m, dd), error = function(e) NULL)
      out <- c(out, as.character(v))
    }
  }
  for (s in p$scales$scales) {
    if (is.character(s$labels)) out <- c(out, s$labels)
    if (is.character(s$name)) out <- c(out, s$name)
  }
  as.character(out)
}
.check_text <- function(parts, label, width_in) {
  sz <- unlist(lapply(parts, .text_sizes))
  .assert(paste0(label, ": no text below 5.5 pt"),
          all(sz[!is.na(sz)] >= MIN_PT - 1e-9), n = sum(!is.na(sz)),
          detail = sprintf("smallest %.2f pt", min(sz, na.rm = TRUE)))
  st <- unlist(lapply(parts, .panel_strings))
  bad <- st[grepl(BAN_RX, st)]
  .assert(paste0(label, ": no beta, m1, m3, M_a or Delta in any label"),
          length(bad) == 0L, n = length(st),
          detail = paste(utils::head(bad, 3), collapse = " / "))
  .assert(paste0(label, ": 178 mm wide"), abs(width_in - W_PAGE) < 1e-9, n = 1L)
}
.no_grey <- function(fills) !any(toupper(fills) == E49_GREY)

# --- ED-g: the three signatures, both models, both measures ------------------
STRIP_TERMS <- c("OX", "MYC", "PROLIF", "DIFF")
STRIP_COL <- c(OX = "b_OX_m3", MYC = "b_MYC_m3", PROLIF = "b_PR", DIFF = "diff")
EDG_ROWS <- tidyr::expand_grid(instrument = names(INSTR),
                               model = c(EST_4G, EST_REF), readout = SIG) %>%
  dplyr::mutate(p = match(instrument, names(INSTR)),
                r = (match(model, c(EST_4G, EST_REF)) - 1) * 3 +
                  match(readout, SIG),
                y = -16 - (p - 1) * 64 - (r - 1) * 8.4 -
                  ifelse(r > 3, 2.6, 0))
EDG_COLS <- .cols49(STRIP_TERMS, x0 = 50, width = W_MM - 50 - 3)
EDG <- e49$displays$ed_g %>%
  dplyr::select(readout, cohort, term, est, lo, hi, instrument, model, x, y,
                label)
EDG_EXP <- EDG %>%
  dplyr::select(readout, cohort, term, instrument, model) %>%
  dplyr::left_join(EDG_ROWS %>% dplyr::select(instrument, model, readout, y),
                   by = c("instrument", "model", "readout")) %>%
  dplyr::left_join(EDG_COLS %>% dplyr::select(term, cohort, x),
                   by = c("term", "cohort"))
EDG_SRC <- dplyr::bind_rows(lapply(STRIP_TERMS, function(tm)
  e49$all_readouts_m3 %>%
    dplyr::filter(readout %in% SIG, estimator %in% c(EST_4G, EST_REF)) %>%
    dplyr::transmute(readout, cohort, term = tm, instrument, model = estimator,
                     src = .data[[STRIP_COL[[tm]]]])))
EDG_CHK <- EDG %>% dplyr::left_join(EDG_SRC, by = c("readout", "cohort", "term",
                                                    "instrument", "model"))
.assert(paste0("ED-g: 96 cells, each where E49 drew it, each value E49's ",
               "stored coefficient, each text E49's"),
        nrow(EDG) == 96L && isTRUE(all.equal(EDG$x, EDG_EXP$x, tolerance = 0)) &&
          isTRUE(all.equal(EDG$y, EDG_EXP$y, tolerance = 0)) &&
          identical(EDG_CHK$est, EDG_CHK$src) &&
          identical(EDG$label, .cell_txt(EDG$est, EDG$lo, EDG$hi)),
        n = nrow(EDG))
LIM_G <- .lim_rule(EDG$est)
EDG <- EDG %>% dplyr::mutate(fill_hex = .fill_hex(est, LIM_G),
                             ink = .ink(fill_hex))
CW_G <- EDG_COLS$cw[1]
H_G  <- 140
EDG_HDR <- dplyr::bind_rows(lapply(names(INSTR), function(inst) {
  top <- -16 - (match(inst, names(INSTR)) - 1) * 64 + 8.4
  dplyr::bind_rows(
    EDG_COLS %>% dplyr::group_by(term) %>%
      dplyr::summarise(x = mean(x), .groups = "drop") %>%
      dplyr::mutate(y = top + 1.8, label = unname(TERM_TITLE49[term]), bold = TRUE),
    EDG_COLS %>% dplyr::transmute(term, x, y = top - 2.2, label = cohort,
                                  bold = FALSE))
}))
EDG_RL <- EDG_ROWS %>% dplyr::mutate(
  label = paste0(SIG_LAB[readout], ", ",
                 ifelse(model == EST_4G, "least-entangled score",
                        "reference score")))
EDG_PANEL <- tibble::tibble(
  y = -16 - (0:1) * 64 + 13.2,
  label = paste0("OXPHOS ", unname(INSTR)))
P_EDG <- ggplot2::ggplot() +
  ggplot2::geom_tile(data = .tag(EDG, "g_col"),
                     ggplot2::aes(x = x, y = y, fill = est),
                     width = CW_G - 0.6, height = 7.6) +
  ggplot2::geom_text(data = .tag(EDG, "g_text"),
                     ggplot2::aes(x = x, y = y, label = label, colour = ink),
                     size = SUB_MM, lineheight = 0.85) +
  ggplot2::scale_colour_identity() +
  ggplot2::geom_text(data = EDG_HDR[EDG_HDR$bold, ],
                     ggplot2::aes(x = x, y = y, label = label),
                     size = TXT_MM, fontface = "bold") +
  ggplot2::geom_text(data = EDG_HDR[!EDG_HDR$bold, ],
                     ggplot2::aes(x = x, y = y, label = label),
                     size = SUB_MM, colour = "grey20") +
  ggplot2::geom_text(data = EDG_RL, ggplot2::aes(x = 48, y = y, label = label),
                     hjust = 1, size = SUB_MM) +
  ggplot2::geom_text(data = EDG_PANEL, ggplot2::aes(x = 2, y = y, label = label),
                     hjust = 0, size = TXT_MM, fontface = "bold") +
  .fig1g_scale(
    LIM_G, breaks = .sym_breaks(LIM_G),
    name = "SD per SD of each score, with the other two held fixed",
    guide = ggplot2::guide_colourbar(direction = "horizontal",
                                     theme = ggplot2::theme(
                                       legend.title.position = "top",
                                       legend.key.width = grid::unit(40, "mm"),
                                       legend.key.height = grid::unit(2.2, "mm")))) +
  ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_MM), ylim = c(-H_G, 0),
                       expand = FALSE, clip = "off") +
  .theme_hm(c(0.985, 0.005))
b_g <- ggplot2::ggplot_build(P_EDG)
f_g <- b_g$data[[.layer_of(P_EDG, "g_col")]]$fill
.assert(paste0("ED-g: every cell is drawn in its actual colour on Fig. 1g's ",
               "palette, OXPHOS and the difference included; nothing is grey"),
        identical(toupper(f_g), .fill_hex(EDG$est, LIM_G)) &&
          identical(toupper(f_g), .ggplot_fill(EDG$est, LIM_G)) && .no_grey(f_g),
        n = length(f_g))
.assert("ED-g: the limits are symmetric, rounded outward to 0.1 and contain every value",
        all(abs(EDG$est) <= LIM_G) &&
          abs(LIM_G * 10 - round(LIM_G * 10)) < 1e-9, n = nrow(EDG),
        detail = sprintf("+/-%.1f", LIM_G))
.assert("ED-g: every cell's text is its stored estimate and interval",
        identical(as.character(b_g$data[[.layer_of(P_EDG, "g_text")]]$label),
                  .cell_txt(EDG$est, EDG$lo, EDG$hi)), n = nrow(EDG))
.check_text(list(P_EDG), "ED-g", W_PAGE)

# --- ED-h: raw correlations beside the adjusted coefficients -----------------
ROWS18 <- tibble::tibble(readout = c(SIG, FIFTEEN),
                         group = c(rep("Signature", 3),
                                   rep(names(GROUPS), lengths(GROUPS)))) %>%
  dplyr::mutate(k = dplyr::row_number(),
                g = match(group, c("Signature", names(GROUPS))),
                y = -32 - (k - 1) * 4.0 - (g - 1) * 1.6)
H_H <- -(min(ROWS18$y) - 2 - 13)
RAW_COLS <- tidyr::expand_grid(cohort = COHORTS, exposure = names(EXPO)) %>%
  dplyr::mutate(c = match(cohort, COHORTS), e = match(exposure, names(EXPO)),
                x = 36 + (c - 1) * (7 * 6.0 + 3) + (e - 0.5) * 6.0)
RAW_T <- e49$displays$ed_h_raw %>%
  dplyr::select(cohort, readout, exposure, rho, x, y, label)
RAW_EXP <- RAW_T %>% dplyr::select(cohort, readout, exposure) %>%
  dplyr::left_join(RAW_COLS %>% dplyr::select(cohort, exposure, x),
                   by = c("cohort", "exposure")) %>%
  dplyr::left_join(ROWS18 %>% dplyr::select(readout, y), by = "readout")
RAW_SRC <- RAW_T %>% dplyr::select(cohort, readout, exposure) %>%
  dplyr::left_join(e49$raw_correlations, by = c("cohort", "readout", "exposure"))
.assert(paste0("ED-h (raw): 252 cells, each where E49 drew it, each value ",
               "E49's stored correlation, each text E49's"),
        nrow(RAW_T) == 252L && isTRUE(all.equal(RAW_T$x, RAW_EXP$x, tolerance = 0)) &&
          isTRUE(all.equal(RAW_T$y, RAW_EXP$y, tolerance = 0)) &&
          identical(RAW_T$rho, RAW_SRC$rho) &&
          identical(RAW_T$label, sprintf("%.2f", RAW_T$rho)), n = nrow(RAW_T))
LIM_R <- .lim_rule(RAW_T$rho)
RAW_T <- RAW_T %>% dplyr::mutate(fill_hex = .fill_hex(rho, LIM_R),
                                 ink = .ink(fill_hex))
W_HL <- max(RAW_COLS$x) + 3 + 2
ADJ_COLS <- tidyr::expand_grid(cohort = COHORTS, term = TERMS) %>%
  dplyr::mutate(c = match(cohort, COHORTS), t = match(term, TERMS),
                x = 3 + (c - 1) * (3 * 6.4 + 3) + (t - 0.5) * 6.4)
ADJ_T <- e49$displays$ed_h_adjusted %>%
  dplyr::select(readout, cohort, term, est, lo, hi, x, y, label)
ADJ_EXP <- ADJ_T %>% dplyr::select(readout, cohort, term) %>%
  dplyr::left_join(ADJ_COLS %>% dplyr::select(cohort, term, x),
                   by = c("cohort", "term")) %>%
  dplyr::left_join(ROWS18 %>% dplyr::select(readout, y), by = "readout")
ADJ_SRC <- dplyr::bind_rows(lapply(TERMS, function(tm)
  e49$all_readouts_m3 %>%
    dplyr::filter(estimator == EST_4G, instrument == "gsva") %>%
    dplyr::transmute(readout, cohort, term = tm, src = .data[[COEF_COL[[tm]]]])))
ADJ_CHK <- ADJ_T %>% dplyr::left_join(ADJ_SRC, by = c("readout", "cohort", "term"))
.assert(paste0("ED-h (adjusted): 108 cells, each where E49 drew it, each ",
               "value 4g's stored coefficient, each text E49's"),
        nrow(ADJ_T) == 108L && isTRUE(all.equal(ADJ_T$x, ADJ_EXP$x, tolerance = 0)) &&
          isTRUE(all.equal(ADJ_T$y, ADJ_EXP$y, tolerance = 0)) &&
          identical(ADJ_CHK$est, ADJ_CHK$src) &&
          identical(ADJ_T$label, sprintf("%.2f", ADJ_T$est)), n = nrow(ADJ_T))
LIM_A <- .lim_rule(ADJ_T$est)
ADJ_T <- ADJ_T %>% dplyr::mutate(fill_hex = .fill_hex(est, LIM_A),
                                 ink = .ink(fill_hex))
W_HR <- W_MM - W_HL
.stop_if(W_HR > max(ADJ_COLS$x) + 3.2, "ED-h: the right part is too narrow")
.edh_left <- function() {
  hdr <- RAW_COLS %>% dplyr::mutate(label = unname(EXPO[exposure]))
  coh <- RAW_COLS %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(x = mean(x), .groups = "drop")
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = .tag(RAW_T, "raw"),
                       ggplot2::aes(x = x, y = y, fill = rho),
                       width = 5.6, height = 3.6) +
    ggplot2::geom_text(data = .tag(RAW_T, "raw_text"),
                       ggplot2::aes(x = x, y = y, label = label, colour = ink),
                       size = SUB_MM) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x = x, y = -29.6, label = label),
                       angle = 90, hjust = 0, size = SUB_MM) +
    ggplot2::geom_text(data = coh, ggplot2::aes(x = x, y = -4.0, label = cohort),
                       size = TXT_MM, fontface = "bold") +
    ggplot2::geom_text(data = ROWS18, ggplot2::aes(
      x = 33, y = y, label = vapply(readout, .row_lab18, character(1))),
      hjust = 1, size = SUB_MM, parse = TRUE) +
    .fig1g_scale(
      LIM_R, breaks = .sym_breaks(LIM_R),
      name = "Spearman correlation, unadjusted",
      guide = ggplot2::guide_colourbar(direction = "horizontal",
                                       theme = ggplot2::theme(
                                         legend.title.position = "top",
                                         legend.key.width = grid::unit(30, "mm"),
                                         legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_HL), ylim = c(-H_H, 0),
                         expand = FALSE, clip = "off") +
    .theme_hm(c(0.98, 0.01))
}
.edh_right <- function() {
  hdr <- ADJ_COLS %>% dplyr::mutate(label = unname(TERM_TITLE[term]))
  coh <- ADJ_COLS %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(x = mean(x), .groups = "drop")
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = .tag(ADJ_T, "adj_col"),
                       ggplot2::aes(x = x, y = y, fill = est),
                       width = 6.0, height = 3.6) +
    ggplot2::geom_text(data = .tag(ADJ_T, "adj_text"),
                       ggplot2::aes(x = x, y = y, label = label, colour = ink),
                       size = SUB_MM) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x = x, y = -29.6, label = label),
                       angle = 90, hjust = 0, size = SUB_MM) +
    ggplot2::geom_text(data = coh, ggplot2::aes(x = x, y = -4.0, label = cohort),
                       size = TXT_MM, fontface = "bold") +
    .fig1g_scale(
      LIM_A, breaks = .sym_breaks(LIM_A),
      name = "SD per SD, other two held fixed",
      guide = ggplot2::guide_colourbar(direction = "horizontal",
                                       theme = ggplot2::theme(
                                         legend.title.position = "top",
                                         legend.key.width = grid::unit(26, "mm"),
                                         legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, W_HR), ylim = c(-H_H, 0),
                         expand = FALSE, clip = "off") +
    .theme_hm(c(0.98, 0.01))
}
P_EDH_L <- .edh_left()
P_EDH_R <- .edh_right()
b_hl <- ggplot2::ggplot_build(P_EDH_L)
b_hr <- ggplot2::ggplot_build(P_EDH_R)
f_hl <- b_hl$data[[.layer_of(P_EDH_L, "raw")]]$fill
f_hr <- b_hr$data[[.layer_of(P_EDH_R, "adj_col")]]$fill
.assert("ED-h (raw): every cell is drawn in its actual colour on Fig. 1g's palette",
        identical(toupper(f_hl), .fill_hex(RAW_T$rho, LIM_R)) &&
          identical(toupper(f_hl), .ggplot_fill(RAW_T$rho, LIM_R)) && .no_grey(f_hl),
        n = length(f_hl))
.assert(paste0("ED-h (adjusted): every cell is drawn in its actual colour on ",
               "Fig. 1g's palette, the signatures' OXPHOS cells included; ",
               "nothing is grey"),
        identical(toupper(f_hr), .fill_hex(ADJ_T$est, LIM_A)) &&
          identical(toupper(f_hr), .ggplot_fill(ADJ_T$est, LIM_A)) &&
          .no_grey(f_hr) && sum(ADJ_T$readout %in% SIG & ADJ_T$term == "OX") == 6L,
        n = length(f_hr))
.assert("ED-h: both limits symmetric, rounded outward to 0.1, containing every value",
        all(abs(RAW_T$rho) <= LIM_R) && all(abs(ADJ_T$est) <= LIM_A) &&
          abs(LIM_R * 10 - round(LIM_R * 10)) < 1e-9 &&
          abs(LIM_A * 10 - round(LIM_A * 10)) < 1e-9, n = 360L,
        detail = sprintf("raw +/-%.1f; adjusted +/-%.1f", LIM_R, LIM_A))
.assert("ED-h: every cell's text is its stored value, as E49 printed it",
        identical(as.character(b_hl$data[[.layer_of(P_EDH_L, "raw_text")]]$label),
                  RAW_T$label) &&
          identical(as.character(b_hr$data[[.layer_of(P_EDH_R, "adj_text")]]$label),
                    ADJ_T$label), n = 360L)
W_EDH <- c(W_HL, W_HR)
.assert("ED-h's two parts sum to 178 mm", abs(sum(W_EDH) - W_MM) < 1e-9, n = 1L)
.check_text(list(P_EDH_L, P_EDH_R), "ED-h", W_PAGE)

# --- save ED-g and ED-h with E49's devices -----------------------------------
# Every PNG is written BEFORE any cairo_pdf. In one R session, a cairo PNG
# drawn after a cairo PNG and a cairo_pdf sets ggplot's two-line cell text
# higher than the first PNG did (found in the dry runs, 2026-11-03), so the two
# copies of ED-g differed. A cairo_pdf drawn first, or none, leaves the PNGs
# alike. Source this script in a freshly restarted R session.
FIG_EDH <- patchwork::wrap_plots(list(P_EDH_L, P_EDH_R), ncol = 2,
                                 widths = W_EDH / sum(W_EDH) * W_PAGE)
RECOL <- list(ed_g = list(P_EDG, H_G / 25.4), ed_h = list(FIG_EDH, H_H / 25.4))
for (nm in names(RECOL)) for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
  ggplot2::ggsave(paste0(.fig_base(dir, FIG[[nm]]), ".png"), RECOL[[nm]][[1]],
                  width = W_PAGE, height = RECOL[[nm]][[2]], units = "in",
                  dpi = 300, bg = "white", type = "cairo")
}
for (nm in names(RECOL)) for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
  ggplot2::ggsave(paste0(.fig_base(dir, FIG[[nm]]), ".pdf"), RECOL[[nm]][[1]],
                  width = W_PAGE, height = RECOL[[nm]][[2]], units = "in",
                  device = grDevices::cairo_pdf, bg = "white")
}
.assert("ED-g and ED-h: each PNG is identical in outputs/figures and docs/figures",
        all(vapply(names(RECOL), function(nm) .same_in_both(FIG[[nm]], ".png"),
                   logical(1))), n = 2L)
for (dir in c(DIR_FIGURES, DIR_DOC_FIGURES)) {
  pre <- if (identical(dir, DIR_DOC_FIGURES)) paste0(FIG_DATE, "_") else ""
  .assert(paste0("ED-g and ED-h (PNG, ", basename(dir), "): every cell's ",
                 "colour appears exactly in the PNG"),
          .png_has_colours(file.path(dir, paste0(pre, FIG[["ed_g"]], ".png")),
                           EDG$fill_hex) &&
            .png_has_colours(file.path(dir, paste0(pre, FIG[["ed_h"]], ".png")),
                             c(RAW_T$fill_hex, ADJ_T$fill_hex)),
          n = nrow(EDG) + nrow(RAW_T) + nrow(ADJ_T))
}
message("   ED-g and ED-h: PNG and PDF, 178 mm wide")

# =============================================================================
# 7. Legend drafts (item 1). Only the grey sentences, and the marked
#    version's sentence, change.
# =============================================================================
message("\n7. legend drafts")
L49 <- e49$legend_draft
NEW_SENTENCE <- paste0("The signature was selected for its association with ",
                       "OXPHOS at fixed proliferation (Methods).")
GREY_4G <- paste0("Bottom, the six-gene signature: its coefficient and 95% ",
                  "confidence interval in each column; its association with ",
                  "OXPHOS, in grey, defines the signature and is shown, not ",
                  "interpreted.")
GREY_G <- paste0("Grey cells contain the OXPHOS coefficient, which defines ",
                 "the signature: OXPHOS itself, and proliferation minus OXPHOS.")
GREY_H <- "Grey, the signatures' association with OXPHOS, which defines them."
DOTS_H <- paste0(" The marked version adds a dot where the 95% confidence ",
                 "interval of an adjusted coefficient excludes zero.")
.n_fixed <- function(s, pat) lengths(regmatches(s, gregexpr(pat, s, fixed = TRUE)))
.assert(paste0("each sentence replaced occurs exactly once in E49's draft: ",
               "4g's, ED-g's and ED-h's grey sentences, and ED-h's marked ",
               "version"),
        .n_fixed(L49$fig4g_a, GREY_4G) == 1L && .n_fixed(L49$ed_g, GREY_G) == 1L &&
          .n_fixed(L49$ed_h, GREY_H) == 1L && .n_fixed(L49$ed_h, DOTS_H) == 1L,
        n = 4L)
LEGEND <- list(
  fig4g = sub(GREY_4G, NEW_SENTENCE, L49$fig4g_a, fixed = TRUE),
  ed_e  = L49$ed_e,
  ed_f  = L49$ed_f,
  ed_g  = sub(GREY_G, NEW_SENTENCE, L49$ed_g, fixed = TRUE),
  ed_h  = sub(DOTS_H, "", sub(GREY_H, NEW_SENTENCE, L49$ed_h, fixed = TRUE),
              fixed = TRUE))
BANNED <- c("independently of MYC", "MYC expression", "comparable", "invariant",
            "lost after adjustment", "primed", "priming",
            "increases reliance on BCL-xL rather than overriding",
            "which BCL-xL can absorb and MCL1 cannot")
.assert(paste0("the author's sentence stands in 4g's, ED-g's and ED-h's ",
               "legends; no legend mentions grey, a dot or a marked version"),
        all(vapply(LEGEND[c("fig4g", "ed_g", "ed_h")], function(s)
          .n_fixed(s, NEW_SENTENCE) == 1L, logical(1))) &&
          !any(grepl("grey|dot|marked", unlist(LEGEND), ignore.case = TRUE)),
        n = length(LEGEND))
.assert(paste0("nothing else in the legends changed, and none uses a word ",
               "the briefs keep out"),
        identical(LEGEND$ed_e, L49$ed_e) && identical(LEGEND$ed_f, L49$ed_f) &&
          identical(sub(NEW_SENTENCE, GREY_4G, LEGEND$fig4g, fixed = TRUE),
                    L49$fig4g_a) &&
          identical(sub(NEW_SENTENCE, GREY_G, LEGEND$ed_g, fixed = TRUE), L49$ed_g) &&
          identical(paste0(sub(NEW_SENTENCE, GREY_H, LEGEND$ed_h, fixed = TRUE),
                           DOTS_H), L49$ed_h) &&
          !any(vapply(BANNED, function(w) any(grepl(w, unlist(LEGEND),
                                                    ignore.case = TRUE)),
                      logical(1))) &&
          !any(grepl("MYC[^.;]{0,30}signature", unlist(LEGEND))),
        n = length(LEGEND))

# =============================================================================
# 8. Save, per R rules 4 and 5
# =============================================================================
message("\n8. saving")
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("   ", nrow(ASSERTIONS), " assertions, all passed")

CELLS_OUT <- dplyr::bind_rows(lapply(names(LAY), function(p) LAY[[p]]$rects %>%
  dplyr::transmute(panel = p, readout, cohort, term, value,
                   printed = .val_ascii(value), fill, ink,
                   contrast = round(contrast, 3))))
readr::write_csv(CELLS_OUT, file.path(DIR_TABLES, "E50_heatmap_cells.csv"))
readr::write_csv(ASSERTIONS, file.path(DIR_TABLES, "E50_assertions.csv"))

saveRDS(list(
  values        = CELLS_OUT,
  layouts       = LAY,
  palette       = list(anchors = DIV_COLS, measured_not_used = MEASURED_COLS,
                       measured_vs_source = PAL_TABLE,
                       interpolation = paste0("CIELAB, as ggplot2's ",
                                              "scale_fill_gradientn in Fig. 1g"),
                       ink = paste0("black, or white where white has the ",
                                    "higher WCAG contrast")),
  axes          = list(heatmaps = L_HM, heatmaps_own = L_OWN, ed_g = LIM_G,
                       ed_h = list(raw = LIM_R, adjusted = LIM_A)),
  ed_g          = EDG,
  ed_h_raw      = RAW_T,
  ed_h_adjusted = ADJ_T,
  output_checks = OUT_CHECKS,
  legend_draft  = LEGEND,
  assertions    = ASSERTIONS,
  # FILE NAMES, NOT PATHS, so this element digests the same in a redirected
  # dry run and in the author's run.
  figure_files  = list(
    heatmaps = paste0(rep(unname(FIG[PANELS$panel]), each = 3),
                      c(".pdf", ".svg", ".png")),
    re_coloured = paste0(rep(unname(FIG[c("ed_g", "ed_h")]), each = 2),
                         c(".png", ".pdf")),
    tracked_prefix = paste0(FIG_DATE, "_")),
  spec          = list(
    specification = paste0("docs/2026-11-03_E50_display_amendment.md, the ",
                           "author's further amendment to the STOP 3 rulings, ",
                           "verbatim, committed alone at 6f5bf4d."),
    posture = "DISPLAY-ONLY. Reads E49's verified object; fits nothing.",
    e49_gate = paste0("Six elements of E49's object, by the digests committed ",
                      "at a63ef6c."),
    palette_source = paste0("myc_mouse 2ebf646, figures/natmetab/_style.R ",
                            "(div_cols, scale_fill_div, lab_num, lab_signed) ",
                            "and figures/natmetab/fig1g_bcl2_heatmap.R."),
    frozen = "Nothing written to myc_human_validation (d3ac60e) or myc_mouse."),
  run_info      = list(built = Sys.time(), bitmap_type = getOption("bitmapType"),
                       r_version = R.version.string)), PATH_E50)

message("   ", PATH_E50)
message("\nE50 done. One scale for 4g, ED-e and ED-f: +/-", sprintf("%.1f", L_HM),
        "; ED-g +/-", sprintf("%.1f", LIM_G), "; ED-h raw +/-",
        sprintf("%.1f", LIM_R), ", adjusted +/-", sprintf("%.1f", LIM_A), ".\n",
        strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E50)

  # the palette: the source used, the measurement recorded
  x$palette$measured_vs_source %>% as.data.frame()
  x$axes

  # 4g's cells as printed, with their colours and ink
  x$values %>% dplyr::filter(panel == "fig4g") %>% as.data.frame()

  # the size of the three panels, in mm
  x$layouts$fig4g$size / (72 / 25.4)

  x$assertions %>% as.data.frame()
  cat(x$legend_draft$fig4g, "\n")
}
