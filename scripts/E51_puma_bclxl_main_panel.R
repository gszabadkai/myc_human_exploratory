# E51_puma_bclxl_main_panel.R
# =============================================================================
# FIG. 1I'S HUMAN ANALOGUE AS A MAIN PANEL: PUMA RELATIVE TO BCL-XL AGAINST
# OXPHOS PRIORITY, WITH TCGA-BRCA AND SCAN-B OVERLAID IN ONE PANEL AND EACH
# COHORT'S P VALUE INSIDE IT.
#
# A DISPLAY SCRIPT WITH ONE DECLARED QUANTITY. It fits nothing and refits
# nothing. It reads E49's verified saved object, checked against the digests
# committed before E49's run (a63ef6c), and draws ED-k's relation again. The
# one quantity not computed before is each cohort's P value, the two-sided
# t-test of the drawn slope in its own model. It is computed exactly from the
# per-tumour residuals E49 saved and that model's residual degrees of freedom.
#
# Specification: the author's ruling of 2026-11-05, verbatim in
# docs/2026-11-05_E51_main_panel_ruling.md, committed alone at 86f38aa before
# this script was written.
#
# SOURCE IT IN A FRESHLY RESTARTED R SESSION. Its PNGs are compared byte for
# byte with the dry run.
#
# =============================================================================
# WHAT THIS DRAWS
# =============================================================================
# The relation is E49's (result note 1.7; declaration 6.8), unchanged:
#   - y: Blom(BBC3) - Blom(BCL2L1), PUMA relative to BCL-xL;
#   - x: OXPHOS priority (mitoPPS);
#   - each adjusted for the reference MYC activity score (Felsher, M_a),
#     proliferation (PROLIF_DISJOINT) and intrinsic subtype, and for tumour
#     purity and leukocyte fraction in TCGA-BRCA;
#   - the slope is E47's m3 priority coefficient difference, BBC3 minus
#     BCL2L1, and its 95% interval is E48 C3's 2,000 paired bootstrap
#     resamples, negated. Both are read as stored.
# The ruling's four choices:
#   1. P is the two-sided t-test of the drawn slope in its own model. By the
#      Frisch-Waugh-Lovell theorem the slope of the residual on the residual
#      is that model's coefficient, and its t statistic, on the model's
#      residual degrees of freedom, is the model's own. It is also the test of
#      the partial correlation, t = r sqrt(df / (1 - r^2)). Both routes are
#      computed here and asserted equal.
#   2. The model is ED-k's, on the reference MYC score, as declared.
#   3. The colours are dark slate #243642 for TCGA-BRCA and mid slate #8FA3B0
#      for SCAN-B, the two ends of Fig. 1's unsigned palette (myc_mouse
#      2ebf646, figures/natmetab/_style.R, seq_cols).
#   4. ED-b stays in Extended Data. Nothing here touches it.
#
# =============================================================================
# FIG. 1I'S STYLE, AND WHAT IS TAKEN FROM ITS CODE
# =============================================================================
# Fig. 1i is myc_mouse 2ebf646, figures/natmetab/fig1i_coupling.R, under
# figures/natmetab/FIGURE_RULES.md. Both were read with git show; nothing is
# written to that repo.
# Taken unchanged from _style.R: NM_FONT, NM_TXT_MIN, NM_TXT, NM_TXT_TITLE,
#   NM_LINE, pt(), lab_num(), lab_signed(), fmt_p(), .nm_axis(), scale_y_nm(),
#   scale_x_nm() and theme_nm(), with the sizes Fig. 1d to 1i use: 5 pt for
#   tick labels and P values, 6 pt for axis titles.
# Taken from fig1i_coupling.R:
#   - zero lines at x = 0 and y = 0, grey85, at NM_LINE;
#   - points of shape 21, size 1.4, stroke 0.2, outlined grey20, filled by
#     group;
#   - the fitted lines as segments over each group's x range, solid, 0.55;
#   - signed axes with minor ticks, ending on ticks;
#   - a square data area, as large as the page height allows, measured on a
#     30-mm square as Fig. 1i measures it.
# The figure rules followed: Helvetica throughout; P exact to two significant
# figures, italic capital P, x 10^k below 0.001 with the exponent at full
# size; a soft hyphen for a true hyphen in the PDF; all text in black.
#
# =============================================================================
# CHOICES THE RULING LEAVES OPEN, made here and reported
# =============================================================================
#   - The page is 165 pt (58.21 mm) tall, the height of Fig. 4g's PDF page
#     (E50 wrote 248 x 165 pt), so that the two can sit in one row. Fig. 1i is
#     56 mm. Its width is the panel's natural width, rounded up to a whole
#     point. Base pdf() writes the page box in whole points, truncating, so each
#     size is given 0.001 pt more than its whole number.
#   - The points are the ten decile means per cohort that E49 stored, as ED-k
#     drew them. The axes do not share a scale: y spans about +/-0.2 and x
#     about +/-1.7, so Fig. 1i's coord_equal() is not used.
#   - The fitted lines pass through the origin, because both residuals have
#     mean 0, and run over each cohort's decile range, as Fig. 1i draws its
#     lines over each group's range.
#   - Both lines are solid at 0.55. The cohorts are equals, so neither is
#     dashed, whereas Fig. 1i dashes WT.
#   - SCAN-B is drawn first and TCGA-BRCA over it.
#   - The P labels sit at the top right of the data area, where no data fall,
#     in black. Each begins with a short line in its cohort's colour, so the
#     labels are the key, and no other key is drawn. There is no title and no
#     R^2.
#   - The axis titles follow Fig. 1i's: "OXPHOS priority (mitoPPS, SD),
#     adjusted" and, on two lines, "BBC3:BCL2L1 (SD)," and "adjusted", with
#     the gene names in italics. They are live text, for the author to edit.
#   - The PDF writes R's minus sign and soft hyphens, which print as a minus
#     and a hyphen. The PNG preview writes a true minus sign (U+2212) and
#     hyphens.
#   - Like Fig. 1i's, the PDF has a white background.
#
# =============================================================================
# WHAT THIS SCRIPT DOES NOT DO
# =============================================================================
#   - It fits nothing and reads no expression data.
#   - It changes no text of the paper. Its legend draft is for the author.
#   - It writes NOTHING to myc_human_validation (d3ac60e) or myc_mouse.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
for (pkg in c("digest", "png", "patchwork", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("E51 needs the package ", pkg, ". Nothing is auto-installed.",
         call. = FALSE)
  }
}

message("\nE51: PUMA relative to BCL-xL against OXPHOS priority, both cohorts ",
        "in one panel\n", strrep("=", 78))
message("Reads E49's verified object and fits nothing. Specification 86f38aa.")

# --- where it reads and writes: one block, so that a dry run can redirect it --
PATH_E49        <- file.path(DIR_RESULTS, "e49_signature_and_fig4g.rds")
PATH_RULING     <- here::here("docs", "2026-11-05_E51_main_panel_ruling.md")
PATH_E51        <- file.path(DIR_RESULTS, "e51_puma_bclxl_main_panel.rds")
OUT_TABLES      <- DIR_TABLES
OUT_FIGURES     <- DIR_FIGURES
DIR_DOC_FIGURES <- here::here("docs", "figures")
FIG_DATE <- "2026-11-05"
# Files are named by content, never by figure letter (CLAUDE.md).
FIG_STEM <- "E51_puma_bclxl_against_priority_overlaid"

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

COHORTS   <- c("TCGA", "SCAN-B")
COH_LABEL <- c(TCGA = "TCGA-BRCA", `SCAN-B` = "SCAN-B")
COH_COL   <- c(TCGA = "#243642", `SCAN-B` = "#8FA3B0")   # the ruling, question 3
DRAW_ORDER <- c("SCAN-B", "TCGA")                        # TCGA-BRCA on top
# ED-k's model (E49 declaration 6.8): OXPHOS priority, the reference MYC score.
MODEL <- list(instrument = "mitopps", estimator = "M_a",
              genes = c("BBC3", "BCL2L1"))

# --- size and type ---------------------------------------------------------------
PAGE_H_PT <- 165L               # Fig. 4g's PDF page height (E50: /MediaBox [0 0 248 165])
PT_EPS   <- 0.001               # pdf() truncates the page box to whole points
PAGE_H   <- (PAGE_H_PT + PT_EPS) * 25.4 / 72   # 58.21 mm
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
LAB_INSET <- 1.0                # mm: the labels' block from the data area's top right
LAB_STEP <- 2.8                 # mm: the two label rows, centre to centre
PNG_DPI  <- 600
SOFT_HY  <- intToUtf8(173L)     # the PDF's hyphen: its encoding prints "-" as minus
MINUS_U  <- intToUtf8(8722L)    # the minus sign, in the PNG
X_TITLE  <- "OXPHOS priority (mitoPPS, SD), adjusted"
Y_TITLE  <- quote(atop(italic("BBC3:BCL2L1") ~ "(SD),", "adjusted"))

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

.assert("the colours are the ruling's, as recorded in the ruling's own file",
        file.exists(PATH_RULING) &&
          all(vapply(COH_COL, function(h) any(grepl(h, readLines(PATH_RULING,
                                                    warn = FALSE), fixed = TRUE)),
                     logical(1))), n = length(COH_COL))

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

# Both axes are symmetric about 0 and end on a minor tick, past every point and
# every line end.
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

# The P labels, one per cohort, for each output. The PDF writes a soft hyphen in
# "TCGA-BRCA"; the PNG writes a hyphen and a true minus in an exponent.
P_LAB <- fmt_p(SLOPES$p[match(COHORTS, SLOPES$cohort)])
LABELS <- list(
  pdf = sprintf("'%s'~~%s", gsub("-", SOFT_HY, COH_LABEL[COHORTS], fixed = TRUE), P_LAB),
  png = sprintf("'%s'~~%s", COH_LABEL[COHORTS],
                gsub("('-", paste0("('", MINUS_U), P_LAB, fixed = TRUE)))

# The printed size of a plotmath label at PT_TXT, in mm, with the PDF device's
# Helvetica metrics (as _style.R's label_mm() measures a P label).
.label_mm <- function(label) {
  grDevices::pdf(NULL, family = NM_FONT)
  on.exit(grDevices::dev.off())
  g <- grid::textGrob(parse(text = label),
                      gp = grid::gpar(fontsize = PT_TXT, fontfamily = NM_FONT))
  c(w = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
    h = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE))
}
LAB_MM <- t(vapply(LABELS$pdf, .label_mm, numeric(2)))

.panel <- function(variant, lab = NULL) {
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
    scale_y_nm(step = Y_STEP, top = YL[2], bottom = YL[1], minor = TRUE, labels = ax) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = X_TITLE, y = Y_TITLE) +
    theme_nm(legend = "none", txt = PT_TXT, txt_title = PT_TITLE)
  if (!is.null(lab)) {
    for (i in seq_along(COHORTS)) {
      k <- COHORTS[i]
      p <- p +
        ggplot2::annotate("segment", x = lab$x_key[1], xend = lab$x_key[2],
                          y = lab$y[i], yend = lab$y[i], colour = COH_COL[[k]],
                          linewidth = LINE_W) +
        ggplot2::annotate("text", x = lab$x_text, y = lab$y[i],
                          label = LABELS[[variant]][i], parse = TRUE, hjust = 0,
                          vjust = 0.5, size = pt(PT_TXT), family = NM_FONT,
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

# The square data area as large as PAGE_H allows, as fig1i_coupling.R finds it:
# measured with a 30-mm square, the rest of the panel is fixed.
REST   <- .natural(.fixed(.panel("pdf"), 30, 30)) - 30
SIDE   <- PAGE_H - PT_EPS * 25.4 / 72 - REST[["height"]]
PAGE_W_PT <- as.integer(ceiling((SIDE + REST[["width"]]) * 72 / 25.4 - 1e-9))
PAGE_W <- (PAGE_W_PT + PT_EPS) * 25.4 / 72

# The labels' places, in data units, now that the data area's size is known.
MM_X <- diff(XL) / SIDE
MM_Y <- diff(YL) / SIDE
LAB <- list(
  x_text = XL[2] - (LAB_INSET + max(LAB_MM[, "w"])) * MM_X,
  y      = YL[2] - (LAB_INSET + max(LAB_MM[, "h"]) / 2 + (seq_along(COHORTS) - 1) * LAB_STEP) * MM_Y)
LAB$x_key <- LAB$x_text - c(KEY_GAP + KEY_LEN, KEY_GAP) * MM_X
BOX <- c(x0 = LAB$x_key[1], x1 = XL[2] - LAB_INSET * MM_X,
         y0 = min(LAB$y) - max(LAB_MM[, "h"]) / 2 * MM_Y, y1 = YL[2] - LAB_INSET * MM_Y)
# Clearance around a point: its radius, about 0.7 mm at size 1.4, and 0.3 mm more.
PAD_X <- 1.0 * MM_X
PAD_Y <- 1.0 * MM_Y
hit_pt <- with(DEC, x > BOX[["x0"]] - PAD_X & x < BOX[["x1"]] + PAD_X &
                    y > BOX[["y0"]] - PAD_Y & y < BOX[["y1"]] + PAD_Y)
xs <- seq(BOX[["x0"]], BOX[["x1"]], length.out = 200)
hit_ln <- unlist(lapply(seq_len(nrow(SEG)), function(i) {
  s <- SEG[i, ]
  inr <- xs >= s$x & xs <= s$xend
  ys <- s$slope * xs[inr]
  ys > BOX[["y0"]] - PAD_Y & ys < BOX[["y1"]] + PAD_Y
}))
.assert("the P labels and their key lines overlap no point and no fitted line",
        !any(hit_pt) && !any(hit_ln), n = nrow(DEC) + nrow(SEG))
.assert("the labels' block lies inside the data area",
        BOX[["x0"]] > XL[1] && BOX[["x1"]] <= XL[2] && BOX[["y0"]] > YL[1] &&
          BOX[["y1"]] <= YL[2], n = 1L)

PANEL <- list(pdf = .fixed(.panel("pdf", LAB), SIDE, SIDE),
              png = .fixed(.panel("png", LAB), SIDE, SIDE))
NEED <- .natural(PANEL$pdf)
.assert("the page's whole points hold the whole panel, and it is Fig. 4g's height",
        NEED[["width"]] <= PAGE_W_PT * 25.4 / 72 + 1e-9 &&
          NEED[["height"]] <= PAGE_H_PT * 25.4 / 72 + 1e-9 && PAGE_H_PT == 165L,
        n = 1L, max_dev = max(NEED - c(PAGE_W_PT, PAGE_H_PT) * 25.4 / 72))
message(sprintf("   data area %.2f mm square; page %d x %d pt (%.2f x %.2f mm)",
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

PDF_TXT <- .pdf_text(FILES[["pdf_docs"]])
mb <- regmatches(PDF_TXT, regexpr("/MediaBox \\[[^]]*\\]", PDF_TXT, useBytes = TRUE))
mb_num <- as.numeric(strsplit(gsub("[^0-9. ]", "", mb), " +")[[1]])
mb_num <- mb_num[!is.na(mb_num)]
.assert("the PDF's page is the panel's size, in whole points",
        identical(as.integer(mb_num[3:4]), c(PAGE_W_PT, PAGE_H_PT)), n = 1L,
        detail = paste(mb_num, collapse = " "))
fonts <- unique(regmatches(PDF_TXT, gregexpr("/BaseFont /[A-Za-z-]+", PDF_TXT, useBytes = TRUE))[[1]])
.assert(paste0("the PDF's text is live: Helvetica, its oblique and Symbol, ",
               "referenced and not embedded"),
        all(c("/BaseFont /Helvetica", "/BaseFont /Helvetica-Oblique") %in% fonts) &&
          all(fonts %in% c("/BaseFont /Helvetica", "/BaseFont /Helvetica-Oblique",
                           "/BaseFont /Helvetica-Bold", "/BaseFont /Helvetica-BoldOblique",
                           "/BaseFont /Symbol")) &&
          !grepl("/FontFile", PDF_TXT, fixed = TRUE, useBytes = TRUE),
        n = length(fonts), detail = paste(sub("/BaseFont /", "", fonts), collapse = ", "))
# Every string the PDF shows: plain strings (Tj) and kerned arrays (TJ), whose
# pieces are joined. Escapes are removed, and the soft hyphen's byte is named.
STR_RE <- "\\((?:[^()\\\\]|\\\\.)*\\)"
.body  <- function(s) sub("^\\((.*)\\)$", "\\1", s, useBytes = TRUE)
tj  <- regmatches(PDF_TXT, gregexpr(paste0(STR_RE, " Tj"), PDF_TXT, perl = TRUE,
                                    useBytes = TRUE))[[1]]
tj  <- .body(sub(" Tj$", "", tj, useBytes = TRUE))
arr <- regmatches(PDF_TXT, gregexpr(paste0("\\[(?:", STR_RE, "|[ 0-9.-])*\\] TJ"),
                                    PDF_TXT, perl = TRUE, useBytes = TRUE))[[1]]
arr <- vapply(arr, function(a) {
  parts <- regmatches(a, gregexpr(STR_RE, a, perl = TRUE, useBytes = TRUE))[[1]]
  paste(.body(parts), collapse = "")
}, character(1), USE.NAMES = FALSE)
shown <- gsub("\\\\([()\\\\])", "\\1", c(tj, arr), useBytes = TRUE)
shown <- gsub(rawToChar(as.raw(0xAD)), "<SHY>", shown, fixed = TRUE, useBytes = TRUE)
mant <- vapply(P_LAB, function(l) sub(".*=='([^']+)'.*", "\\1", l), character(1),
               USE.NAMES = FALSE)
want <- c("TCGA<SHY>BRCA", "SCAN<SHY>B", "P", mant, X_TITLE, "BBC3:BCL2L1")
.assert(paste0("the PDF shows each cohort's name, with a soft hyphen, each P ",
               "value as computed, and both axis titles"),
        all(want %in% shown), n = length(want),
        detail = paste(setdiff(want, shown), collapse = ", "))

FIGURE_FILES <- tibble::tibble(
  copy = names(FILES), file = basename(FILES),
  checksum = c(.md5(FILES[["png_outputs"]]), .md5(FILES[["png_docs"]]),
               .pdf_md5(FILES[["pdf_outputs"]]), .pdf_md5(FILES[["pdf_docs"]])),
  checksum_is = c("md5", "md5", "md5, timestamps zeroed", "md5, timestamps zeroed"))

# =============================================================================
# 6. The legend draft, for the author
# =============================================================================
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
.assert("the legend draft is ASCII, as the script is",
        !grepl("[^ -~]", LEGEND_DRAFT), n = 1L)

# =============================================================================
# 7. Saved: the object, two tables and the digests
# =============================================================================
ASSERTIONS <- dplyr::bind_rows(ASSERT)
.stop_if(all(ASSERTIONS$pass), "an assertion is recorded as failed")
message("\n   ", nrow(ASSERTIONS), " assertions, all passed")

SPEC <- list(
  ruling = "docs/2026-11-05_E51_main_panel_ruling.md (86f38aa)",
  model = MODEL, colours = COH_COL, draw_order = DRAW_ORDER,
  test = "two-sided t-test of the slope, on the model's residual degrees of freedom",
  page_pt = c(width = PAGE_W_PT, height = PAGE_H_PT), side_mm = SIDE,
  type_pt = c(text = PT_TXT, titles = PT_TITLE), line = LINE_W,
  points = list(shape = 21L, size = PT_SIZE, stroke = PT_STROKE, edge = PT_EDGE),
  axes = list(x = XL, y = YL, x_step = X_STEP, y_step = Y_STEP),
  x_title = X_TITLE, y_title = deparse(Y_TITLE),
  source_style = "myc_mouse 2ebf646: figures/natmetab/_style.R, fig1i_coupling.R")
LAYOUT <- list(labels = LABELS, label_mm = LAB_MM, label_pos = LAB, label_box = BOX,
               rest_mm = REST, natural_mm = NEED)

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
saveRDS(OUT, PATH_E51)
utils::write.csv(SLOPES, file.path(OUT_TABLES, "E51_slopes_and_p.csv"), row.names = FALSE)
utils::write.csv(ASSERTIONS, file.path(OUT_TABLES, "E51_assertions.csv"), row.names = FALSE)

DIGESTS <- vapply(setdiff(names(OUT), "run_info"),
                  function(k) digest::digest(OUT[[k]]), character(1))
message("\n   digests of the saved elements (run_info excluded):")
for (k in names(DIGESTS)) message(sprintf("     %-13s %s", k, DIGESTS[[k]]))
message("\n   the legend draft:\n   ", LEGEND_DRAFT)
message("\nE51 done.")

# =============================================================================
# SANDBOX -- run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {
  # Self-contained: it reads only the saved object and the tracked figure, so it
  # also works in a fresh session (after a restart, say), not only straight after
  # source(). PANEL exists only in the session that sourced the script.
  library(dplyr)
  x <- readRDS(here::here("results", "e51_puma_bclxl_main_panel.rds"))
  x$slopes %>% dplyr::select(cohort, n, df, slope, slope_lo, slope_hi, partial_r, t, p) %>% print()
  x$assertions %>% print(n = Inf)
  x$figure_files %>% print()
  cat(x$legend_draft, "\n")
  # the panel as written: the PNG preview (the PDF beside it is the file for
  # Illustrator). Straight after source(), print(PANEL$png) draws it afresh.
  img <- png::readPNG(here::here("docs", "figures",
                                 "2026-11-05_E51_puma_bclxl_against_priority_overlaid.png"))
  grid::grid.newpage()
  grid::grid.raster(img)
}
