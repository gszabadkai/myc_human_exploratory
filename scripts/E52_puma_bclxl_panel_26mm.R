# E52_puma_bclxl_panel_26mm.R
# =============================================================================
# E51'S PANEL AT THE SIZE OF FIG. 4'S PRISM GRAPHS: A 26-MM AXIS, AXIS TITLES
# WITHOUT "ADJUSTED", AND THE TWO P LABELS IN A KEY ABOVE THE DATA AREA.
#
# A DISPLAY SCRIPT. It fits nothing and computes nothing that E51 did not. It
# reads E49's verified saved object (digests committed at a63ef6c) and repeats
# E51's computation of each cohort's P exactly. It asserts that the slopes,
# deciles, fitted lines and legend draft it draws are E51's verified ones, by
# their digests (E51's pre-run note, 1b0be30; verified at beead93).
#
# Specification: the author's ruling of 2026-11-05, verbatim in
# docs/2026-11-05_E52_prism_size_ruling.md, committed alone at 9c2ca44 before
# this script was written. E51's ruling (86f38aa) stands except where it
# changes it.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNGs are compared byte for
# byte with the dry run.
#
# =============================================================================
# WHAT CHANGES FROM E51, AND WHY
# =============================================================================
#   - The data area is 26 x 26 mm, the axis length of Fig. 4's Prism graphs
#     (the ruling). E51's was 46.75 mm, at Fig. 4g's height. The axis lines'
#     length is read back from the PDF and asserted.
#   - The page is the panel's natural size, rounded up to whole points. This is
#     the Prism format's rule that the saved page is the panel, with no margin
#     to trim (myc_mouse 2ebf646, figures/natmetab/FIGURE_RULES.md).
#   - The axis titles lose "adjusted" (the ruling), and the legend carries the
#     adjustment.
#     - The x title, 30.2 mm at 6 pt, is longer than the 26-mm axis. It is set
#       on two lines, "OXPHOS priority" and "(mitoPPS, SD)", 2.4 mm baseline to
#       baseline, as the Prism format asks ("6 pt at 1.15") and as Fig. 1h sets
#       it (fig1h_bcl2_fade.R, LEAD <- 2.4). R spaces lines 1.2 x the type size
#       apart at lineheight 1, so the lineheight is 0.945, not 1.15.
#     - The y title, "BBC3:BCL2L1 (SD)" with the gene names in italics, is
#       18.6 mm and stays on one line.
#   - The two P labels, with the cohort names and their key lines, form a key
#     above the data area, left-aligned with it (the ruling, question 1). At
#     26 mm, inside the data area, they would sit on the top SCAN-B points.
#     TCGA-BRCA's row is on top.
#   - The axes and ticks stay at 0.53 pt, NM_LINE (the ruling, question 2).
# Unchanged from E51:
#   - the data, the P values, the colours, the points and the lines;
#   - the type sizes: 5 pt for ticks and P, 6 pt for titles;
#   - Fig. 1i's helpers, copied from myc_mouse 2ebf646;
#   - the PDF with live Helvetica, beside a cairo PNG preview.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing and reads no expression data.
#   - It changes no text of the paper. Its legend draft is E51's, unchanged.
#   - It leaves no raw file bytes in the workspace. The PDF is read back inside
#     a function that returns only ASCII (E51's lesson).
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
for (pkg in c("digest", "png", "patchwork", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E52 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE52: E51's panel at the Prism graphs' size (a 26-mm axis)\n",
        strrep("=", 78))
message("Reads E49's verified object and fits nothing. Specifications 86f38aa, 9c2ca44.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_E49        <- file.path(DIR_RESULTS, "e49_signature_and_fig4g.rds")
PATH_RULING     <- here::here("docs", "2026-11-05_E52_prism_size_ruling.md")
PATH_RULING_E51 <- here::here("docs", "2026-11-05_E51_main_panel_ruling.md")
PATH_E52        <- file.path(DIR_RESULTS, "e52_puma_bclxl_panel_26mm.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-05"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- "E52_puma_bclxl_against_priority_26mm"

# =============================================================================
# 0. Constants
# =============================================================================
# The elements of E49's object read here, with the digests recorded in
# docs/2026-11-03_e49_prerun_digests.md, committed at a63ef6c before the
# author's run. That run reproduced all 24.
E49_DIGESTS <- c(
  fig1i           = "f9a71873b59f3d12a9dad7641c4c0eda",
  all_readouts_m3 = "7023883b873b17020e800838db591684",
  assertions      = "2710612444bdbe6592173817b73c2d86")
# E51's saved elements that E52 draws again, with the digests in E51's pre-run
# note (1b0be30), which E51's run reproduced (beead93).
E51_DIGESTS <- c(
  slopes       = "dec549b6e83e26be49aaaf9d49ce34a5",
  deciles      = "85136d7800b63e5dd8d69b064c90acd8",
  lines        = "ae9cdfec16f9f8682d3eea90a37f70fe",
  legend_draft = "8ccf7b2f800aad901bc5d4cc199f9288")

COHORTS   <- c("TCGA", "SCAN-B")
COH_LABEL <- c(TCGA = "TCGA-BRCA", `SCAN-B` = "SCAN-B")
COH_COL   <- c(TCGA = "#243642", `SCAN-B` = "#8FA3B0")   # E51's ruling, question 3
DRAW_ORDER <- c("SCAN-B", "TCGA")                        # TCGA-BRCA on top
# ED-k's model (E49 declaration 6.8): OXPHOS priority, the reference MYC score.
MODEL <- list(instrument = "mitopps", estimator = "M_a",
              genes = c("BBC3", "BCL2L1"))

# --- size and type ---------------------------------------------------------------
SIDE     <- 26                  # mm: the data area, the Prism graphs' axis length (the ruling)
PT_EPS   <- 0.001               # pdf() truncates the page box to whole points
PT_TXT   <- 5                   # tick labels and P values (Fig. 1d to 1i)
PT_TITLE <- 6                   # axis titles
LINE_W   <- 0.55                # fig1i_coupling.R: the solid line
PT_SIZE  <- 1.4                 # fig1i_coupling.R: the points
PT_STROKE <- 0.2
PT_EDGE  <- "grey20"
ZERO_COL <- "grey85"
X_STEP   <- 1
Y_STEP   <- 0.2
KEY_LEN  <- 2.5                 # mm: the key line before each P label
KEY_GAP  <- 0.8                 # mm: key line to text
KEY_ABOVE <- 1.2                # mm: the key's lower baseline above the data area's top
LAB_STEP <- 2.8                 # mm: the key's rows, baseline to baseline
CAP_EM   <- 0.718               # Helvetica's capitals above the baseline, in em (_style.R)
PNG_DPI  <- 600
SOFT_HY  <- intToUtf8(173L)     # the PDF's hyphen: its encoding prints "-" as minus
MINUS_U  <- intToUtf8(8722L)    # the minus sign, in the PNG
X_TITLE  <- "OXPHOS priority\n(mitoPPS, SD)"
X_TITLE_GAP <- 2.4              # mm, baseline to baseline (the Prism format; Fig. 1h's LEAD)
# R spaces lines of text 1.2 x the type size apart at lineheight 1 (the
# devices' character height), so the lineheight that gives the gap is 0.945.
X_TITLE_LH <- X_TITLE_GAP / (1.2 * PT_TITLE * 25.4 / 72)
Y_TITLE  <- quote(italic("BBC3:BCL2L1") ~ "(SD)")
AXIS_TOL_PT <- 0.02             # the PDF writes coordinates to 0.01 pt

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

.assert("the colours are E51's ruling's, as recorded in that ruling's own file",
        file.exists(PATH_RULING_E51) &&
          all(vapply(COH_COL, function(h) any(grepl(h, readLines(PATH_RULING_E51,
                                                    warn = FALSE), fixed = TRUE)),
                     logical(1))), n = length(COH_COL))
.assert("the size is the ruling's, a 26-mm axis, as recorded in the ruling's own file",
        file.exists(PATH_RULING) && SIDE == 26 &&
          any(grepl("axis length should be 26 mm", readLines(PATH_RULING, warn = FALSE),
                    fixed = TRUE)), n = 1L)

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

# =============================================================================
# 2. The relation as stored, and each cohort's P
# =============================================================================
message("\n2. The stored relation, and the P values")
SUM <- as.data.frame(e49$fig1i$summary)
DEC <- as.data.frame(e49$fig1i$deciles)
RES <- as.data.frame(e49$fig1i$residuals)
.assert("E49 stored both cohorts, with ten deciles each",
        setequal(SUM$cohort, COHORTS) && nrow(SUM) == 2L &&
          all(as.integer(table(DEC$cohort)[COHORTS]) == 10L),
        n = nrow(DEC))

FITS <- e49$all_readouts_m3 %>%
  dplyr::filter(instrument == MODEL$instrument, estimator == MODEL$estimator,
                readout %in% MODEL$genes) %>%
  dplyr::select(cohort, readout, n, df_m3, b_OX_m3)
.assert(paste0("ED-k's model: in each cohort the BBC3 and BCL2L1 fits share ",
               "one design, with the same n and residual degrees of freedom"),
        nrow(FITS) == 4L &&
          all(vapply(COHORTS, function(k) {
            f <- FITS[FITS$cohort == k, ]
            nrow(f) == 2L && length(unique(f$n)) == 1L && length(unique(f$df_m3)) == 1L
          }, logical(1))), n = nrow(FITS))
dev_gc2 <- max(vapply(COHORTS, function(k) {
  f <- FITS[FITS$cohort == k, ]
  abs((f$b_OX_m3[f$readout == "BBC3"] - f$b_OX_m3[f$readout == "BCL2L1"]) -
        SUM$slope[SUM$cohort == k])
}, numeric(1)))
.assert(paste0("the stored slope is BBC3's OXPHOS priority coefficient minus ",
               "BCL2L1's, in each cohort (E49's G-C2)"),
        dev_gc2 < 1e-12, n = 2L, max_dev = dev_gc2)

SLOPES <- dplyr::bind_rows(lapply(COHORTS, function(k) {
  r  <- RES[RES$cohort == k, ]
  s  <- SUM[SUM$cohort == k, ]
  f  <- FITS[FITS$cohort == k & FITS$readout == "BBC3", ]
  x  <- r$x_resid
  y  <- r$y_resid
  b  <- sum(x * y) / sum(x^2)
  e  <- y - b * x
  df <- f$df_m3
  se <- sqrt(sum(e^2) / df / sum(x^2))
  rr <- stats::cor(x, y)
  tibble::tibble(
    cohort = k, n = nrow(r), n_model = f$n, df = df,
    slope = s$slope, slope_lo = s$slope_lo, slope_hi = s$slope_hi,
    interval_source = s$interval_source, partial_r = s$partial_r,
    slope_from_residuals = b, r_from_residuals = rr,
    mean_x_resid = mean(x), mean_y_resid = mean(y),
    se = se, t = b / se, t_from_r = rr * sqrt(df / (1 - rr^2)),
    p = 2 * stats::pt(-abs(b / se), df))
}))

.assert("each cohort's residuals number the model's n",
        all(SLOPES$n == SLOPES$n_model) && all(SLOPES$n == SUM$n[match(SLOPES$cohort, SUM$cohort)]),
        n = 2L)
dev_mean <- max(abs(c(SLOPES$mean_x_resid, SLOPES$mean_y_resid)))
.assert("both residuals have mean 0 in each cohort, so the lines pass through the origin",
        dev_mean < 1e-10, n = 4L, max_dev = dev_mean)
dev_slope <- max(abs(SLOPES$slope_from_residuals - SLOPES$slope))
.assert(paste0("the slope of the residual on the residual equals the stored ",
               "model coefficient (Frisch-Waugh-Lovell)"),
        dev_slope < 1e-10, n = 2L, max_dev = dev_slope)
dev_r <- max(abs(SLOPES$r_from_residuals - SLOPES$partial_r))
.assert("the residuals' correlation equals the stored partial correlation",
        dev_r < 1e-10, n = 2L, max_dev = dev_r)
dev_t <- max(abs(SLOPES$t - SLOPES$t_from_r))
.assert(paste0("the slope's t, from its standard error, equals the partial ",
               "correlation's t, on the model's residual degrees of freedom"),
        dev_t < 1e-8, n = 2L, max_dev = dev_t)
.assert("each stored bootstrap interval contains its slope",
        all(SLOPES$slope_lo < SLOPES$slope & SLOPES$slope < SLOPES$slope_hi), n = 2L)
.assert("each P lies strictly between 0 and 1",
        all(SLOPES$p > 0 & SLOPES$p < 1), n = 2L)

RECOMP <- RES %>% dplyr::group_by(cohort, decile) %>%
  dplyr::summarise(n2 = dplyr::n(), x2 = mean(x_resid), y2 = mean(y_resid),
                   .groups = "drop") %>%
  dplyr::inner_join(DEC, by = c("cohort", "decile"))
dev_dec <- max(abs(c(RECOMP$x2 - RECOMP$x, RECOMP$y2 - RECOMP$y)))
.assert("every stored decile mean is the mean of its tumours' residuals",
        nrow(RECOMP) == 20L && all(RECOMP$n2 == RECOMP$n) && dev_dec < 1e-12,
        n = nrow(RECOMP), max_dev = dev_dec)

.assert("the slopes and P values are E51's verified ones, by digest",
        identical(digest::digest(SLOPES), E51_DIGESTS[["slopes"]]), n = 1L)

for (k in COHORTS) {
  s <- SLOPES[SLOPES$cohort == k, ]
  message(sprintf("   %-9s n = %5d, df = %5d: slope %+.3f [%+.3f, %+.3f]; partial r %+.3f; t = %.3f; P = %.2g",
                  COH_LABEL[[k]], s$n, s$df, s$slope, s$slope_lo, s$slope_hi,
                  s$partial_r, s$t, s$p))
}

# =============================================================================
# 3. Fig. 1i's style helpers, copied unchanged from myc_mouse 2ebf646,
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
lab_signed <- function(x) {
  out <- lab_num(x)
  pos <- !is.na(x) & x > 0
  out[pos] <- paste0("+", out[pos])
  out
}

fmt_p <- function(p, exponent = c("full", "reduced")) {
  style <- c(full = "displaystyle", reduced = "scriptstyle")[[match.arg(exponent)]]
  stopifnot(is.numeric(p), !anyNA(p), all(p >= 0 & p <= 1))
  vapply(p, function(v) {
    if (v >= 0.001) {
      return(sprintf("italic(P)=='%s'", formatC(signif(v, 2), format = "fg", digits = 2, flag = "#")))
    }
    if (v < 1e-300) return(sprintf("italic(P)<10^%s('-300')", style))
    e <- floor(log10(v))
    m <- signif(v / 10^e, 2)
    if (m >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("italic(P)=='%s'%%*%%10^%s('%d')", format(m), style, e)
  }, character(1))
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

scale_y_nm <- function(step, top, bottom = 0, minor = FALSE, labels = lab_num, ...) {
  a <- .nm_axis(step, top, bottom, minor)
  ggplot2::scale_y_continuous(limits = a$limits, breaks = a$breaks, minor_breaks = a$minor_breaks,
                              guide = a$guide, expand = c(0, 0), labels = labels, ...)
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
# 4. The panel
# =============================================================================
message("\n4. The panel")
DEC_DRAW <- DEC %>%
  dplyr::mutate(cohort = factor(cohort, levels = DRAW_ORDER)) %>%
  dplyr::arrange(cohort, decile)
SEG <- dplyr::bind_rows(lapply(COHORTS, function(k) {
  d <- DEC[DEC$cohort == k, ]
  b <- SLOPES$slope[SLOPES$cohort == k]
  tibble::tibble(cohort = k, x = min(d$x), xend = max(d$x),
                 y = b * min(d$x), yend = b * max(d$x), slope = b)
})) %>%
  dplyr::mutate(cohort = factor(cohort, levels = DRAW_ORDER)) %>%
  dplyr::arrange(cohort)
.assert("the points and lines drawn are E51's verified ones, by digest",
        identical(digest::digest(DEC_DRAW), E51_DIGESTS[["deciles"]]) &&
          identical(digest::digest(SEG), E51_DIGESTS[["lines"]]), n = 2L)

# Both axes are E51's: symmetric about 0, ending on a minor tick past every
# point and every line end.
X_TOP <- ceiling(max(abs(c(DEC$x, SEG$x, SEG$xend))) / (X_STEP / 2) - 1e-9) * (X_STEP / 2)
Y_TOP <- ceiling(max(abs(c(DEC$y, SEG$y, SEG$yend))) / (Y_STEP / 2) - 1e-9) * (Y_STEP / 2)
XA <- .nm_axis(X_STEP, X_TOP, -X_TOP, TRUE)
YA <- .nm_axis(Y_STEP, Y_TOP, -Y_TOP, TRUE)
XL <- XA$limits
YL <- YA$limits
.assert("every point and line end lies inside both axes",
        all(DEC$x > XL[1] & DEC$x < XL[2] & DEC$y > YL[1] & DEC$y < YL[2]) &&
          all(c(SEG$x, SEG$xend) >= XL[1] & c(SEG$x, SEG$xend) <= XL[2] &
                c(SEG$y, SEG$yend) >= YL[1] & c(SEG$y, SEG$yend) <= YL[2]),
        n = nrow(DEC) + 2L * nrow(SEG))
.assert("the axes are E51's: x from -2 to +2, and y from -0.3 to +0.3, as the ruling keeps it",
        isTRUE(all.equal(XL, c(-2, 2))) && isTRUE(all.equal(YL, c(-0.3, 0.3))), n = 2L)

# The P labels, one per cohort, for each output, exactly as E51 wrote them. The
# PDF writes a soft hyphen in "TCGA-BRCA"; the PNG writes a hyphen and a true
# minus in an exponent.
P_LAB <- fmt_p(SLOPES$p[match(COHORTS, SLOPES$cohort)])
LABELS <- list(
  pdf = sprintf("'%s'~~%s", gsub("-", SOFT_HY, COH_LABEL[COHORTS], fixed = TRUE), P_LAB),
  png = sprintf("'%s'~~%s", COH_LABEL[COHORTS],
                gsub("('-", paste0("('", MINUS_U), P_LAB, fixed = TRUE)))

# The printed size of a plotmath label at PT_TXT, in mm, with the PDF device's
# Helvetica metrics (as _style.R's label_mm() measures a P label): its width, its
# height, and how far it reaches below its baseline (the round letters'
# overshoot).
.label_mm <- function(label) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- grid::textGrob(parse(text = label),
                      gp = grid::gpar(fontsize = PT_TXT, fontfamily = NM_FONT))
  c(w = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
    h = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE),
    d = grid::convertHeight(grid::grobDescent(g), "mm", valueOnly = TRUE))
}
LAB_MM <- t(vapply(LABELS$pdf, .label_mm, numeric(3)))

# The key: one row per cohort above the data area, TCGA-BRCA's on top, each a key
# line in its cohort's colour and the label in black, left-aligned with the data
# area. The rows are set on their baselines, LAB_STEP apart, and each key line
# is level with the middle of its row's capitals, so that the exponent in
# SCAN-B's row does not lift its line. The margin above the data area holds the
# key and keeps theme_nm's 1 mm above it; theme_nm's other three margins are kept.
M_NM <- theme_nm()$plot.margin
.assert(paste0("theme_nm's margins are 1, 2.5, 1 and 1 mm (top, right, bottom, ",
               "left); E52 replaces only the top, keeping its 1 mm above the key"),
        identical(as.numeric(M_NM), c(1, 2.5, 1, 1)) && all(grid::unitType(M_NM) == "mm"),
        n = 4L)
ROW     <- rev(seq_along(COHORTS)) - 1                   # from the bottom: SCAN-B 0, TCGA-BRCA 1
BASE_MM <- KEY_ABOVE + ROW * LAB_STEP                    # each baseline, above the data area
ASC_MM  <- LAB_MM[, "h"] - LAB_MM[, "d"]                 # each label above its baseline
CAP_MM  <- CAP_EM * PT_TXT * 25.4 / 72                   # the capitals' height, 1.27 mm
W_KEY   <- KEY_LEN + KEY_GAP + max(LAB_MM[, "w"])        # mm
TOP_MM  <- max(BASE_MM + ASC_MM) + as.numeric(M_NM)[1]   # mm
.assert("the key fits within the data area's width, left-aligned with it",
        W_KEY <= SIDE, n = 1L,
        detail = sprintf("%.2f mm of %.0f", W_KEY, SIDE))
lower <- order(ROW)
.assert("each key row clears the row above it: no label reaches the next baseline",
        all(ASC_MM[lower][-length(lower)] + LAB_MM[lower, "d"][-1] < LAB_STEP),
        n = length(COHORTS) - 1L,
        detail = sprintf("%.2f mm of %.1f",
                         max(ASC_MM[lower][-length(lower)] + LAB_MM[lower, "d"][-1]),
                         LAB_STEP))

.panel <- function(variant, key = NULL) {
  ax <- if (variant == "png") function(x) sub("^-", MINUS_U, lab_signed(x)) else lab_signed
  p <- ggplot2::ggplot(DEC_DRAW, ggplot2::aes(x, y)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = NM_LINE, colour = ZERO_COL) +
    ggplot2::geom_vline(xintercept = 0, linewidth = NM_LINE, colour = ZERO_COL) +
    ggplot2::geom_segment(data = SEG, ggplot2::aes(x = x, xend = xend, y = y,
                                                   yend = yend, colour = cohort),
                          linewidth = LINE_W, inherit.aes = FALSE) +
    ggplot2::geom_point(ggplot2::aes(fill = cohort), shape = 21, size = PT_SIZE,
                        stroke = PT_STROKE, colour = PT_EDGE) +
    ggplot2::scale_colour_manual(values = COH_COL, guide = "none") +
    ggplot2::scale_fill_manual(values = COH_COL, guide = "none") +
    scale_x_nm(step = X_STEP, top = XL[2], bottom = XL[1], minor = TRUE, labels = ax) +
    # oob_keep: the key sits above the y axis's top, and the default would drop it.
    scale_y_nm(step = Y_STEP, top = YL[2], bottom = YL[1], minor = TRUE, labels = ax,
               oob = scales::oob_keep) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = X_TITLE, y = Y_TITLE) +
    theme_nm(legend = "none", txt = PT_TXT, txt_title = PT_TITLE) +
    ggplot2::theme(plot.margin = ggplot2::margin(TOP_MM, 2.5, 1, 1, "mm"),
                   axis.title.x = ggplot2::element_text(lineheight = X_TITLE_LH))
  if (!is.null(key)) {
    for (i in seq_along(COHORTS)) {
      k <- COHORTS[i]
      p <- p +
        ggplot2::annotate("segment", x = key$x_key[1], xend = key$x_key[2],
                          y = key$y_line[i], yend = key$y_line[i], colour = COH_COL[[k]],
                          linewidth = LINE_W) +
        # vjust = 0 sets the label's lowest point at y_text, its overshoot below
        # the baseline, so that the baseline itself falls where it is meant to.
        ggplot2::annotate("text", x = key$x_text, y = key$y_text[i],
                          label = LABELS[[variant]][i], parse = TRUE, hjust = 0,
                          vjust = 0, size = pt(PT_TXT), family = NM_FONT,
                          colour = "black")
    }
  }
  p
}

.fixed <- function(p, w_mm, h_mm) {
  patchwork::wrap_plots(p) +
    patchwork::plot_layout(widths = grid::unit(w_mm, "mm"),
                           heights = grid::unit(h_mm, "mm"))
}
.natural <- function(pw) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- patchwork::patchworkGrob(pw)
  c(width  = sum(grid::convertWidth(g$widths, "mm", valueOnly = TRUE)),
    height = sum(grid::convertHeight(g$heights, "mm", valueOnly = TRUE)))
}

# The key's places, in data units: the data area is SIDE mm on each side.
MM_X <- diff(XL) / SIDE
MM_Y <- diff(YL) / SIDE
KEY <- list(
  x_key  = XL[1] + c(0, KEY_LEN) * MM_X,
  x_text = XL[1] + (KEY_LEN + KEY_GAP) * MM_X,
  y_text = unname(YL[2] + (BASE_MM - LAB_MM[, "d"]) * MM_Y),
  y_line = YL[2] + (BASE_MM + CAP_MM / 2) * MM_Y)

PANEL <- list(pdf = .fixed(.panel("pdf", KEY), SIDE, SIDE),
              png = .fixed(.panel("png", KEY), SIDE, SIDE))

# The page is the panel's natural size, rounded up to whole points: pdf()
# truncates its page box to whole points, so the page is set PT_EPS over.
NEED <- .natural(PANEL$pdf)
PAGE_W_PT <- as.integer(ceiling(NEED[["width"]]  * 72 / 25.4 - 1e-9))
PAGE_H_PT <- as.integer(ceiling(NEED[["height"]] * 72 / 25.4 - 1e-9))
PAGE_W <- (PAGE_W_PT + PT_EPS) * 25.4 / 72
PAGE_H <- (PAGE_H_PT + PT_EPS) * 25.4 / 72
over_mm <- c(PAGE_W_PT, PAGE_H_PT) * 25.4 / 72 - NEED
.assert("the page is the panel's natural size, rounded up to whole points",
        all(over_mm >= -1e-9) && all(over_mm < 25.4 / 72), n = 2L,
        detail = sprintf("%.3f and %.3f mm over it", over_mm[1], over_mm[2]))
message(sprintf("   data area %.0f mm square; page %d x %d pt (%.2f x %.2f mm)",
                SIDE, PAGE_W_PT, PAGE_H_PT, PAGE_W_PT * 25.4 / 72, PAGE_H_PT * 25.4 / 72))

# =============================================================================
# 5. The files: PNGs first, then the PDFs (base pdf(), live Helvetica)
# =============================================================================
message("\n5. The files")
FILES <- c(
  png_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM, ".png")),
  png_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM, ".png")),
  pdf_outputs = file.path(OUT_FIGURES, paste0(FIG_STEM, ".pdf")),
  pdf_docs    = file.path(DIR_DOC_FIGURES, paste0(FIG_DATE, "_", FIG_STEM, ".pdf")))
for (f in FILES[c("png_outputs", "png_docs")]) {
  ggplot2::ggsave(f, PANEL$png, width = PAGE_W, height = PAGE_H, units = "mm",
                  dpi = PNG_DPI, device = "png", type = "cairo", bg = "white")
}
.pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(filename, width = width, height = height, family = NM_FONT,
                 useDingbats = FALSE, compress = FALSE)
}
for (f in FILES[c("pdf_outputs", "pdf_docs")]) {
  ggplot2::ggsave(f, PANEL$pdf, width = PAGE_W, height = PAGE_H, units = "mm",
                  device = .pdf_device)
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
.assert("the two PNGs are identical", .md5(FILES[["png_outputs"]]) == .md5(FILES[["png_docs"]]),
        n = 2L)
.assert("the two PDFs are identical, timestamps aside",
        .pdf_md5(FILES[["pdf_outputs"]]) == .pdf_md5(FILES[["pdf_docs"]]), n = 2L)

PNG <- png::readPNG(FILES[["png_docs"]])
dims_ok <- abs(dim(PNG)[1] - PAGE_H / 25.4 * PNG_DPI) < 1 &&
  abs(dim(PNG)[2] - PAGE_W / 25.4 * PNG_DPI) < 1
.assert("the PNG is the page's size at 600 dpi", dims_ok, n = 1L,
        detail = paste(dim(PNG)[1:2], collapse = " x "))
.hexes <- function(a) {
  rgb <- a[, , 1:3, drop = FALSE]
  unique(grDevices::rgb(rgb[, , 1], rgb[, , 2], rgb[, , 3]))
}
PNG_HEX <- toupper(.hexes(PNG))
.assert("both cohorts' colours appear exactly in the PNG",
        all(toupper(COH_COL) %in% PNG_HEX), n = length(COH_COL))

# The PDF is read back INSIDE a function, so that its raw bytes, which are not
# valid UTF-8 (the binary header, soft hyphens, Symbol's multiplication sign),
# never reach the workspace. Left there, they crashed Positron's R session each
# time it listed the workspace after a run (E51 result note). Only ASCII and
# numbers come back: the page box, the fonts, every string the PDF shows with
# its place, and every straight two-point line.
.pdf_read_back <- function(f) {
  txt <- .pdf_text(f)
  mb <- regmatches(txt, regexpr("/MediaBox \\[[^]]*\\]", txt, useBytes = TRUE))
  mb_num <- as.numeric(strsplit(gsub("[^0-9. ]", "", mb), " +")[[1]])
  mb_num <- mb_num[!is.na(mb_num)]
  fonts <- unique(regmatches(txt, gregexpr("/BaseFont /[A-Za-z-]+", txt, useBytes = TRUE))[[1]])
  # Every string the PDF shows: plain strings (Tj) and kerned arrays (TJ), whose
  # pieces are joined. Escapes are removed, the soft hyphen's byte is named, and
  # any other byte above 127 is written as <xx>.
  str_re <- "\\((?:[^()\\\\]|\\\\.)*\\)"
  body   <- function(v) sub("^\\((.*)\\)$", "\\1", v, useBytes = TRUE)
  decode <- function(v) {
    v <- gsub("\\\\([()\\\\])", "\\1", v, useBytes = TRUE)
    v <- gsub(rawToChar(as.raw(0xAD)), "<SHY>", v, fixed = TRUE, useBytes = TRUE)
    iconv(v, from = "latin1", to = "ASCII", sub = "byte")
  }
  tj  <- regmatches(txt, gregexpr(paste0(str_re, " Tj"), txt, perl = TRUE,
                                  useBytes = TRUE))[[1]]
  tj  <- body(sub(" Tj$", "", tj, useBytes = TRUE))
  arr <- regmatches(txt, gregexpr(paste0("\\[(?:", str_re, "|[ 0-9.-])*\\] TJ"),
                                  txt, perl = TRUE, useBytes = TRUE))[[1]]
  arr <- vapply(arr, function(a) {
    parts <- regmatches(a, gregexpr(str_re, a, perl = TRUE, useBytes = TRUE))[[1]]
    paste(body(parts), collapse = "")
  }, character(1), USE.NAMES = FALSE)
  shown <- decode(c(tj, arr))
  # Each string's place. pdf() writes every string in its own BT ... ET block,
  # after its text matrix (a b c d x y Tm), so (x, y) is the start of its
  # baseline, and b is 0 unless the text is rotated.
  bt <- regmatches(txt, gregexpr("(?s)(?<![A-Za-z])BT\\s.*?\\sET(?![A-Za-z])", txt,
                                 perl = TRUE, useBytes = TRUE))[[1]]
  num6 <- paste(rep("(-?[0-9.]+)", 6), collapse = " ")
  placed <- do.call(rbind, lapply(bt, function(b) {
    m <- regmatches(b, regexec(paste0(num6, " Tm"), b, useBytes = TRUE))[[1]]
    s <- regmatches(b, gregexpr(str_re, b, perl = TRUE, useBytes = TRUE))[[1]]
    if (length(m) != 7L || length(s) == 0L) return(NULL)
    data.frame(s = decode(paste(body(s), collapse = "")),
               x = as.numeric(m[6]), y = as.numeric(m[7]),
               rotated = as.numeric(m[3]) != 0, stringsAsFactors = FALSE)
  }))
  # Every straight line drawn from one point to one other (x0 y0 m x1 y1 l S).
  sg <- regmatches(txt, gregexpr(paste0("(?<![0-9.-])(-?[0-9.]+) (-?[0-9.]+) m\\s+",
                                        "(-?[0-9.]+) (-?[0-9.]+) l\\s+S"),
                                 txt, perl = TRUE, useBytes = TRUE))[[1]]
  segments <- matrix(vapply(strsplit(sg, "\\s+"), function(v) as.numeric(v[c(1, 2, 4, 5)]),
                            numeric(4)), ncol = 4L, byrow = TRUE,
                     dimnames = list(NULL, c("x0", "y0", "x1", "y1")))
  list(mediabox = mb_num, fonts = fonts,
       font_file = grepl("/FontFile", txt, fixed = TRUE, useBytes = TRUE),
       shown = shown, placed = placed, n_text_blocks = length(bt),
       segments = segments)
}
RB <- .pdf_read_back(FILES[["pdf_docs"]])
.assert("the PDF's page is the panel's size, in whole points",
        identical(as.integer(RB$mediabox[3:4]), c(PAGE_W_PT, PAGE_H_PT)), n = 1L,
        detail = paste(RB$mediabox, collapse = " "))
.assert(paste0("the PDF's text is live: Helvetica, its oblique and Symbol, ",
               "referenced and not embedded"),
        all(c("/BaseFont /Helvetica", "/BaseFont /Helvetica-Oblique") %in% RB$fonts) &&
          all(RB$fonts %in% c("/BaseFont /Helvetica", "/BaseFont /Helvetica-Oblique",
                              "/BaseFont /Helvetica-Bold", "/BaseFont /Helvetica-BoldOblique",
                              "/BaseFont /Symbol")) &&
          !RB$font_file,
        n = length(RB$fonts), detail = paste(sub("/BaseFont /", "", RB$fonts), collapse = ", "))
mant <- vapply(P_LAB, function(l) sub(".*=='([^']+)'.*", "\\1", l), character(1),
               USE.NAMES = FALSE)
X_LINES <- strsplit(X_TITLE, "\n", fixed = TRUE)[[1]]
want <- c("TCGA<SHY>BRCA", "SCAN<SHY>B", "P", mant, X_LINES, "BBC3:BCL2L1", "(SD)")
.assert(paste0("the PDF shows each cohort's name, with a soft hyphen, each P ",
               "value as computed, and both axis titles"),
        all(want %in% RB$shown), n = length(want),
        detail = paste(setdiff(want, RB$shown), collapse = ", "))
.assert("no text in the PDF contains \"adjusted\" (the ruling)",
        !any(grepl("adjusted", RB$shown, fixed = TRUE)), n = length(RB$shown))
.assert("every string the PDF shows was found with its place",
        nrow(RB$placed) == length(RB$shown) && RB$n_text_blocks == length(RB$shown) &&
          setequal(RB$placed$s, RB$shown), n = length(RB$shown))

# The data area, measured in the PDF. Its axis lines and its two zero lines span
# it, so they are the only straight lines 26 mm long, and nothing is longer.
SIDE_PT <- SIDE * 72 / 25.4
KEY_PT  <- KEY_LEN * 72 / 25.4
SG <- RB$segments
H_SG <- SG[abs(SG[, "y1"] - SG[, "y0"]) < 1e-9, , drop = FALSE]
V_SG <- SG[abs(SG[, "x1"] - SG[, "x0"]) < 1e-9, , drop = FALSE]
len_h <- abs(H_SG[, "x1"] - H_SG[, "x0"])
len_v <- abs(V_SG[, "y1"] - V_SG[, "y0"])
SPAN_H <- H_SG[abs(len_h - SIDE_PT) <= AXIS_TOL_PT, , drop = FALSE]
SPAN_V <- V_SG[abs(len_v - SIDE_PT) <= AXIS_TOL_PT, , drop = FALSE]
AREA <- c(x0 = min(SPAN_V[, "x0"]), y0 = min(SPAN_H[, "y0"]))
AREA <- c(AREA, x1 = AREA[["x0"]] + SIDE_PT, y1 = AREA[["y0"]] + SIDE_PT)
dev_side <- max(abs(c(len_h[abs(len_h - SIDE_PT) <= AXIS_TOL_PT],
                      len_v[abs(len_v - SIDE_PT) <= AXIS_TOL_PT]) - SIDE_PT))
.assert(paste0("in the PDF the data area is 26 mm square: two horizontal and two ",
               "vertical lines (the axis lines and the zero lines) measure ",
               sprintf("%.2f", SIDE_PT), " pt, start on its edges, and no ",
               "straight line is longer"),
        nrow(SPAN_H) == 2L && nrow(SPAN_V) == 2L &&
          max(len_h) <= SIDE_PT + AXIS_TOL_PT && max(len_v) <= SIDE_PT + AXIS_TOL_PT &&
          all(abs(pmin(SPAN_H[, "x0"], SPAN_H[, "x1"]) - AREA[["x0"]]) <= AXIS_TOL_PT) &&
          all(abs(pmin(SPAN_V[, "y0"], SPAN_V[, "y1"]) - AREA[["y0"]]) <= AXIS_TOL_PT),
        n = 4L, max_dev = dev_side)
TX <- RB$placed
.y_of <- function(s) TX$y[TX$s == s & !TX$rotated]
.mm2pt <- function(v) v * 72 / 25.4
Y_T <- .y_of("TCGA<SHY>BRCA")
Y_S <- .y_of("SCAN<SHY>B")
dev_base <- if (length(Y_T) == 1L && length(Y_S) == 1L) {
  max(abs(c(Y_S - AREA[["y1"]] - .mm2pt(KEY_ABOVE), Y_T - Y_S - .mm2pt(LAB_STEP))))
} else Inf
.assert(paste0("in the PDF the key's rows lie above the data area, TCGA-BRCA's on ",
               "top: the lower baseline 1.2 mm above the data area, and the two ",
               "baselines 2.8 mm apart"),
        is.finite(dev_base) && dev_base <= AXIS_TOL_PT, n = 2L, max_dev = dev_base)
KEY_SG <- H_SG[abs(len_h - KEY_PT) <= AXIS_TOL_PT, , drop = FALSE]
KEY_SG <- KEY_SG[order(-KEY_SG[, "y0"]), , drop = FALSE]      # the top row's first
dev_line <- if (nrow(KEY_SG) == 2L && is.finite(dev_base)) {
  max(abs(KEY_SG[, "y0"] - c(Y_T, Y_S) - .mm2pt(CAP_MM / 2)))
} else Inf
.assert(paste0("in the PDF each key line starts at the data area's left edge, ",
               "level with the middle of its row's capitals"),
        nrow(KEY_SG) == 2L && dev_line <= AXIS_TOL_PT &&
          all(abs(pmin(KEY_SG[, "x0"], KEY_SG[, "x1"]) - AREA[["x0"]]) <= AXIS_TOL_PT),
        n = 2L, max_dev = dev_line)
Y_X1 <- .y_of(X_LINES[1])
Y_X2 <- .y_of(X_LINES[2])
GAP_X <- if (length(Y_X1) == 1L && length(Y_X2) == 1L) Y_X1 - Y_X2 else NA_real_
GAP_PT <- X_TITLE_GAP * 72 / 25.4
.assert(paste0("in the PDF the x title's two lines lie below the data area, ",
               "2.4 mm (", sprintf("%.2f", GAP_PT), " pt) apart, baseline to baseline"),
        is.finite(GAP_X) && abs(GAP_X - GAP_PT) <= AXIS_TOL_PT &&
          Y_X1 < AREA[["y0"]], n = 1L, max_dev = abs(GAP_X - GAP_PT))

FIGURE_FILES <- tibble::tibble(
  copy = names(FILES), file = basename(FILES),
  checksum = c(.md5(FILES[["png_outputs"]]), .md5(FILES[["png_docs"]]),
               .pdf_md5(FILES[["pdf_outputs"]]), .pdf_md5(FILES[["pdf_docs"]])),
  checksum_is = c("md5", "md5", "md5, timestamps zeroed", "md5, timestamps zeroed"))

# =============================================================================
# 6. The legend draft, for the author: E51's, unchanged
# =============================================================================
# It already states the adjustment that leaves the axis titles: "each adjusted
# for MYC activity on the reference score, proliferation and intrinsic subtype,
# and for tumour purity and leukocyte fraction in TCGA-BRCA".
.num <- function(v) sprintf("%.2f", v)
.n   <- function(v) format(v, big.mark = ",", trim = TRUE)
S_T <- SLOPES[SLOPES$cohort == "TCGA", ]
S_S <- SLOPES[SLOPES$cohort == "SCAN-B", ]
LEGEND_DRAFT <- paste0(
  "BBC3 (PUMA) relative to BCL2L1 (BCL-xL) against OXPHOS priority (mitoPPS), ",
  "each adjusted for MYC activity on the reference score, proliferation and ",
  "intrinsic subtype, and for tumour purity and leukocyte fraction in ",
  "TCGA-BRCA; TCGA-BRCA (n = ", .n(S_T$n), "), dark, and SCAN-B (n = ",
  .n(S_S$n), "), light. Each cohort is standardised within itself, and the ",
  "two are not pooled. Points, means of the ten deciles of adjusted OXPHOS ",
  "priority; lines, the fitted slopes, ", .num(S_T$slope), " (95% CI ",
  .num(S_T$slope_lo), " to ", .num(S_T$slope_hi), ") in TCGA-BRCA and ",
  .num(S_S$slope), " (", .num(S_S$slope_lo), " to ", .num(S_S$slope_hi),
  ") in SCAN-B, with intervals from 2,000 paired bootstrap resamples of ",
  "tumours. P, two-sided t-test of the slope.")
.assert("the legend draft is E51's verified one, by digest",
        identical(digest::digest(LEGEND_DRAFT), E51_DIGESTS[["legend_draft"]]), n = 1L)
.assert("the legend draft states the adjustment that leaves the axis titles",
        grepl("each adjusted for MYC activity", LEGEND_DRAFT, fixed = TRUE), n = 1L)

# =============================================================================
# 7. Saved: the object, two tables and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

SPEC <- list(
  ruling = c(e52 = "docs/2026-11-05_E52_prism_size_ruling.md (9c2ca44)",
             e51 = "docs/2026-11-05_E51_main_panel_ruling.md (86f38aa), where E52's does not change it"),
  model = MODEL, colours = COH_COL, draw_order = DRAW_ORDER,
  test = "two-sided t-test of the slope, on the model's residual degrees of freedom",
  side_mm = SIDE, page_pt = c(width = PAGE_W_PT, height = PAGE_H_PT),
  type_pt = c(text = PT_TXT, titles = PT_TITLE), line = LINE_W, axis_line = NM_LINE,
  points = list(shape = 21L, size = PT_SIZE, stroke = PT_STROKE, edge = PT_EDGE),
  axes = list(x = XL, y = YL, x_step = X_STEP, y_step = Y_STEP),
  x_title = X_TITLE, x_title_gap_mm = X_TITLE_GAP, x_title_lineheight = X_TITLE_LH,
  y_title = deparse(Y_TITLE),
  key = list(place = "above the data area, left-aligned with it",
             row_order = COHORTS, rows_set_on = "baselines", line_mm = KEY_LEN,
             gap_mm = KEY_GAP, baseline_above_mm = KEY_ABOVE,
             baseline_step_mm = LAB_STEP, line_at = "half the capitals' height",
             cap_em = CAP_EM, margin_above_mm = as.numeric(M_NM)[1]),
  source_style = "myc_mouse 2ebf646: figures/natmetab/_style.R, fig1i_coupling.R")
LAYOUT <- list(labels = LABELS, label_mm = LAB_MM, key_pos = KEY, key_width_mm = W_KEY,
               top_margin_mm = TOP_MM, natural_mm = NEED, rest_mm = NEED - SIDE,
               pdf_area_pt = AREA, pdf_key_y_pt = c(TCGA = Y_T, `SCAN-B` = Y_S),
               pdf_x_title_gap_pt = GAP_X)

OUT <- list(
  slopes = SLOPES,
  deciles = DEC_DRAW,
  lines = SEG,
  layout = LAYOUT,
  legend_draft = LEGEND_DRAFT,
  assertions = ASSERTIONS,
  figure_files = FIGURE_FILES,
  spec = SPEC,
  run_info = list(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  R = R.version.string,
                  bitmap_type_default = getOption("bitmapType")))
saveRDS(OUT, PATH_E52)
utils::write.csv(SLOPES, file.path(OUT_TABLES, "E52_slopes_and_p.csv"), row.names = FALSE)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E52_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), "run_info"),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-13s %s", k, DIGESTS[[k]]))
message("\n   the legend draft (E51's, unchanged):\n   ", LEGEND_DRAFT)
message("\nE52 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figure, so it
  # also works in a fresh session, not only straight after source(). PANEL
  # exists only in the session that sourced the script.
  library(dplyr)
  x <- readRDS(here::here("results", "e52_puma_bclxl_panel_26mm.rds"))
  x$slopes %>% dplyr::select(cohort, n, df, slope, slope_lo, slope_hi, partial_r, t, p) %>% print()
  x$assertions %>% print(n = Inf)
  x$figure_files %>% print()
  x$spec$page_pt
  x$layout$pdf_area_pt
  cat(x$legend_draft, "\n")
  # the panel as written: the PNG preview (the PDF beside it is the file for
  # Illustrator). Straight after source(), print(PANEL$png) draws it afresh.
  img <- png::readPNG(here::here("docs", "figures",
                                 "2026-11-05_E52_puma_bclxl_against_priority_26mm.png"))
  grid::grid.newpage()
  grid::grid.raster(img)
}
