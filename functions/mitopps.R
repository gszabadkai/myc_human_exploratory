# functions/mitopps.R
# =============================================================================
# mitoPPS - the composition instrument. Sourced by E02 (which builds the scores)
# and E03b (which tests whether one of them is an arithmetic artefact).
# Defined once, here, so the two cannot drift apart.
#
# WHAT IT MEASURES, AND WHAT IT CANNOT
# ------------------------------------
# A pathway's mitoPPS is the mean, over every OTHER pathway in the universe, of
# its pairwise expression ratio, each ratio corrected by that ratio's own
# cross-sample mean. It answers "is this programme PRIORITISED relative to the
# rest of the mitochondrion", not "is this programme high".
#
# It is therefore BLIND TO LEVEL BY DESIGN (CLAUDE.md trap 6), its baseline is
# composition-dependent, and its values are never comparable between cohorts -
# only its patterns are.
#
# THE COROLLARY THAT E03b EXISTS TO TEST: because the scores are ratios against
# a shared denominator, a programme can be pushed negative purely because
# ANOTHER programme in the universe rose. That is an arithmetic consequence, not
# biology, and telling the two apart requires holding the denominator fixed.
#
# INPUT SCALE: linear DESeq2-normalised. A pathway score is the mean linear
# expression of its genes and must be strictly positive - the ratio is otherwise
# undefined, and both functions refuse rather than return Inf.
# =============================================================================

# Every pathway of the universe against every other. S: pathways x samples.
.mitopps_universe <- function(S) {
  stopifnot(is.matrix(S), nrow(S) > 1L)
  if (min(S) <= 0) {
    stop("non-positive pathway score reached .mitopps_universe", call. = FALSE)
  }
  N <- ncol(S); P <- nrow(S)
  Bi <- 1 / S
  A  <- (S %*% t(Bi)) / N
  out <- S * (((1 / A) %*% Bi) - Bi) / (P - 1)
  dimnames(out) <- dimnames(S)
  out
}

# Query pathways against a HELD-OUT universe. Sq and Su share samples but not
# rows; when Su is the universe with Sq's own pathway removed, the result is
# identical to that pathway's .mitopps_universe value. E02 asserts that identity
# and it is what proves the query form is right.
.mitopps_query <- function(Sq, Su) {
  stopifnot(is.matrix(Sq), is.matrix(Su), ncol(Sq) == ncol(Su))
  if (min(Sq) <= 0 || min(Su) <= 0) {
    stop("non-positive pathway score reached .mitopps_query", call. = FALSE)
  }
  N <- ncol(Su); P <- nrow(Su)
  Bi  <- 1 / Su
  A   <- (Sq %*% t(Bi)) / N
  out <- Sq * ((1 / A) %*% Bi) / P
  dimnames(out) <- list(rownames(Sq), colnames(Sq))
  out
}

# Mean LINEAR expression per set. The pathway score mitoPPS consumes.
.path_scores <- function(sets, L, min_genes = 3L) {
  keep <- sets[vapply(sets, length, integer(1)) >= min_genes]
  out  <- t(vapply(keep, function(g) colMeans(L[g, , drop = FALSE]),
                   numeric(ncol(L))))
  rownames(out) <- names(keep); colnames(out) <- colnames(L)
  out
}
