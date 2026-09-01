# functions/strata.R
# =============================================================================
# The sample strata, defined once. E03, E05, E07 and E08 all cut the cohorts the
# same way, and three copies of that definition is how they stop agreeing.
#
# LUMINAL (LumA + LumB) is a deliberate addition, made 2026-09-01. The mouse arm
# shows luminal expansion, so the combined luminal compartment is a question in
# its own right and not merely the sum of two answers - LumA is 499/1540 samples
# and LumB 197/896, and each alone is underpowered for a gene-level breakdown
# where the combined 696/2436 is not.
#
# ADDING IT IS PURELY ADDITIVE. Every stratum that existed before is defined
# exactly as before, so re-running E03 or E05 adds rows and changes no cell that
# either findings note quotes.
#
# PAM50 is NA for 114 TCGA and 64 SCAN-B samples and ER is NA for 51 and 17.
# Those samples are in `all` and in nothing else. They are never imputed.
# =============================================================================

.build_strata <- function(frames, cohort, ids) {
  f <- frames[frames$cohort == cohort, ]
  f <- f[match(ids, f$sample_id), ]
  s <- list(
    all     = ids,
    ERpos   = ids[!is.na(f$ER) & f$ER == "ERpos"],
    ERneg   = ids[!is.na(f$ER) & f$ER == "ERneg"],
    Luminal = ids[!is.na(f$PAM50) & f$PAM50 %in% c("LumA", "LumB")])
  for (p in levels(f$PAM50)) s[[p]] <- ids[!is.na(f$PAM50) & f$PAM50 == p]
  s
}
