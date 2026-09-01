# functions/correlation_engine.R
# =============================================================================
# The ranked partial-correlation engine behind E03 (the atlas) and E05 (the
# death axis). Sourced, never run on its own.
#
# WHY RANKS ONCE, THEN LINEAR ALGEBRA
# -----------------------------------
# Spearman rho is Pearson r on ranks. A partial Spearman is Pearson r on the
# residuals of those ranks. So the whole grid - tens of thousands of cells -
# reduces to: rank each row once per (cohort, stratum, complete-case set),
# project out the covariates once, then ONE cor() call per block. Calling
# cor(method = "spearman") cell by cell would re-rank the same vectors
# thousands of times for an identical answer.
#
# SCALE IS NOT AT ISSUE HERE, AND THAT IS WORTH SAYING OUT LOUD. Every measure
# entering this engine is rank-transformed, and ranks are invariant under any
# monotone transform. So a Spearman correlation of a LOG score against a LINEAR
# one is well defined and carries no scale error. The scale discipline in
# CLAUDE.md binds where scores are BUILT (E02), not where they are correlated.
#
# CONFIDENCE INTERVALS. Fisher-z with the Bonett-Wright variance for Spearman,
# se = sqrt((1 + rho^2 / 2) / (n - 3 - k)), k = number of covariates projected
# out. The plain 1/(n-3) understates the variance of a rank correlation.
#
# THE p COLUMN IS NOT A RESULT. It is carried so a cell can be described, never
# so a cell can be selected. CLAUDE.md: the atlas is a grid of tens of thousands
# of cells and any single one of them is uninteresting. Report structure,
# gradients and cross-cohort reproducibility.
# =============================================================================

# Rank every ROW of a matrix (measures x samples), ties averaged.
.rank_rows <- function(M) {
  out <- t(apply(M, 1L, rank))
  if (nrow(M) == 1L) out <- matrix(out, nrow = 1L, dimnames = dimnames(M))
  dimnames(out) <- dimnames(M)
  out
}

# Residualise every ROW of R (rows x n) on the design Z (n x p, intercept
# included by the caller). Returns the residuals, same shape as R.
.resid_rows <- function(R, Z) {
  R - t(qr.fitted(qr(Z), t(R)))
}

# One block of the atlas: every estimator against every measure, on one sample
# set, under one adjustment.
#
#   EST   estimators x samples. MAY carry NA (TCGA's M_c_call does).
#   MEAS  measures   x samples. Must not carry NA.
#   ids   the samples of this stratum, as column names.
#   cov   samples x k covariate matrix (rownames = sample ids), or NULL.
#   min_n cells below this are not emitted at all.
#
# Estimators are grouped by their MISSINGNESS PATTERN and each group gets its
# own complete-case sample set and its own ranking. Dropping every sample that
# any estimator is missing would silently shrink n for the 20 estimators that
# are complete, so that is not done.
.atlas_block <- function(EST, MEAS, ids, cov = NULL, min_n = 30L) {
  stopifnot(all(ids %in% colnames(EST)), all(ids %in% colnames(MEAS)))
  E <- EST[, ids, drop = FALSE]
  M <- MEAS[, ids, drop = FALSE]
  if (anyNA(M)) {
    stop("the measure matrix carries NA - correlations would be silently ",
         "pairwise and n would not mean what the column says", call. = FALSE)
  }
  Z0 <- if (is.null(cov)) NULL else cov[ids, , drop = FALSE]
  ok_cov <- if (is.null(Z0)) rep(TRUE, length(ids)) else stats::complete.cases(Z0)

  na_key <- apply(is.na(E), 1L, function(v) paste(which(v), collapse = ","))
  out <- list()
  for (key in unique(na_key)) {
    rows <- which(na_key == key)
    ok   <- ok_cov & !apply(is.na(E[rows, , drop = FALSE]), 2L, any)
    n    <- sum(ok)
    if (n < min_n) next

    Er <- .rank_rows(E[rows, ok, drop = FALSE])
    Mr <- .rank_rows(M[, ok, drop = FALSE])
    k  <- 0L
    if (!is.null(Z0)) {
      Zr <- apply(Z0[ok, , drop = FALSE], 2L, rank)
      Z  <- cbind(`(Intercept)` = 1, Zr)
      k  <- ncol(Zr)
      Er <- .resid_rows(Er, Z)
      Mr <- .resid_rows(Mr, Z)
    }
    # A measure or estimator that is constant within this stratum has no
    # correlation to report. Named, not silently turned into NA.
    sdE <- apply(Er, 1L, stats::sd); sdM <- apply(Mr, 1L, stats::sd)
    if (any(sdE == 0) || any(sdM == 0)) {
      warning("constant within a stratum, dropped: ",
              paste(c(rownames(Er)[sdE == 0], rownames(Mr)[sdM == 0]),
                    collapse = ", "), call. = FALSE)
      Er <- Er[sdE > 0, , drop = FALSE]; Mr <- Mr[sdM > 0, , drop = FALSE]
      if (!nrow(Er) || !nrow(Mr)) next
    }

    rho <- stats::cor(t(Er), t(Mr))
    z   <- atanh(pmin(pmax(rho, -0.999999999), 0.999999999))
    se  <- sqrt((1 + rho^2 / 2) / (n - 3 - k))
    out[[key]] <- tibble::tibble(
      myc_estimator = rep(rownames(rho), times = ncol(rho)),
      measure       = rep(colnames(rho), each  = nrow(rho)),
      n             = n,
      k_cov         = k,
      rho           = as.vector(rho),
      ci_lo         = as.vector(tanh(z - 1.959964 * se)),
      ci_hi         = as.vector(tanh(z + 1.959964 * se)),
      p             = as.vector(2 * stats::pnorm(-abs(z / se))))
  }
  if (!length(out)) return(NULL)
  dplyr::bind_rows(out)
}

# Spearman of EVERY row of L against one vector, optionally partialled on
# covariates. Chunked so an 18,000 x 3,200 rank matrix never exists all at once.
#
# This is the background an expression-matched null is drawn against, so it must
# cover the WHOLE matrix, not just a gene set. A gene with no variance has no
# correlation and is returned NA rather than NaN - the caller must drop those
# from the set and the background alike, or the draws are biased.
#
# cov = NULL reproduces a plain Spearman exactly; with cov it is the same
# partial Spearman .atlas_block computes (Pearson on residuals of ranks).
.per_gene_rho <- function(L, y, cov = NULL, chunk = 2000L) {
  stopifnot(length(y) == ncol(L))
  n <- length(y)
  Z <- NULL
  if (!is.null(cov)) {
    stopifnot(nrow(cov) == n, !anyNA(cov))
    Z <- cbind(`(Intercept)` = 1, apply(cov, 2L, rank))
  }
  ry <- rank(y)
  if (!is.null(Z)) ry <- as.numeric(.resid_rows(matrix(ry, nrow = 1L), Z))
  ry <- (ry - mean(ry)) / stats::sd(ry)
  out <- numeric(nrow(L)); names(out) <- rownames(L)
  for (s in seq(1, nrow(L), by = chunk)) {
    e  <- min(s + chunk - 1L, nrow(L))
    Rx <- .rank_rows(L[s:e, , drop = FALSE])
    if (!is.null(Z)) Rx <- .resid_rows(Rx, Z)
    sdx <- apply(Rx, 1L, stats::sd)
    Rx <- (Rx - rowMeans(Rx)) / sdx
    v <- as.numeric(Rx %*% ry) / (n - 1)
    v[sdx == 0] <- NA_real_
    out[s:e] <- v
  }
  out
}

# Split a background into equal-count expression ventiles, ready for matched
# draws. Returns bin index per row (NA where excluded) and the pools.
.expression_bins <- function(L, keep, n_bins = 20L) {
  expr <- rowMeans(L)
  bins <- rep(NA_integer_, nrow(L))
  bins[keep] <- cut(rank(expr[keep], ties.method = "first"),
                    breaks = n_bins, labels = FALSE)
  list(bins = bins, pools = split(keep, bins[keep]), n_bins = n_bins)
}

# One expression-matched draw for a set of row indices. `exclude` lets a paired
# draw keep its two halves disjoint, as the real strata are.
.matched_draw <- function(idx, B, exclude = integer(0)) {
  want <- table(factor(B$bins[idx], levels = seq_len(B$n_bins)))
  unlist(lapply(seq_len(B$n_bins), function(b) {
    k <- want[[b]]
    if (k == 0L) return(integer(0))
    pool <- setdiff(B$pools[[as.character(b)]], exclude)
    if (length(pool) < k) {
      stop("expression ventile ", b, " holds ", length(pool),
           " usable background genes but the set needs ", k,
           "; the matched draw cannot be made without replacement.",
           call. = FALSE)
    }
    pool[sample.int(length(pool), k)]
  }), use.names = FALSE)
}

