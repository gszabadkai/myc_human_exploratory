# functions/gene_matrix.R
# =============================================================================
# Pulling a named set of HUMAN GENE SYMBOLS out of a cohort's expression matrix,
# through that cohort's symbol map, with the misses NAMED rather than dropped.
#
# WHY THIS IS NOT A ONE-LINER. SCAN-B's symbols are a 2014 UCSC build
# (CLAUDE.md trap 7). A gene asked for by its current symbol is simply ABSENT
# from that matrix unless it is resolved through scanb_pheno.rds$symbol_map -
# and absence is silent. `M[genes, ]` on a symbol that is not a rowname errors,
# but `intersect()` on the same list quietly returns a shorter set and the
# analysis continues with fewer genes than it claims. That is the failure mode
# this file exists to make impossible: every miss comes back in `$missing` and
# the caller is expected to print it.
#
# A symbol that resolves to MORE than one row is also a miss, not a guess.
# E02's symbol map is built to be injective, but a map is an input and inputs
# change; picking the first row of two would be a silent wrong answer.
#
# E08 carries its own copy of this logic, written before this file existed. The
# two are identical in behaviour and must stay so - the next edit to E08 should
# delete its local .in_t / .in_s / .gene_rows and source this file instead.
# =============================================================================

# Build a resolver for one cohort: a function from one requested symbol to the
# matrix rownames it maps to (zero, one, or - a problem - several).
#
#   rn   the matrix rownames
#   map  named character vector, current symbol -> matrix symbol, or NULL for a
#        cohort that needs no mapping (TCGA)
.symbol_resolver <- function(rn, map = NULL) {
  force(rn); force(map)
  function(g) {
    g <- unique(g)
    h <- if (is.null(map)) g else {
      m <- map[g]
      unname(ifelse(is.na(m), g, m))
    }
    intersect(h, rn)
  }
}

# genes x samples submatrix, ROWS NAMED BY THE REQUESTED SYMBOL rather than by
# whatever the matrix calls it, so downstream joins are on one vocabulary.
.gene_rows <- function(genes, L, resolve) {
  genes <- unique(genes)
  h <- vapply(genes, function(g) {
    r <- resolve(g)
    if (length(r) != 1L) NA_character_ else r
  }, character(1))
  ok <- !is.na(h)
  M <- L[h[ok], , drop = FALSE]
  rownames(M) <- genes[ok]
  list(mat = M, missing = genes[!ok])
}
