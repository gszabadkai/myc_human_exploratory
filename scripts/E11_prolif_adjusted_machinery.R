# E11_prolif_adjusted_machinery.R
# =============================================================================
# PHASE 2. ONE CLAIM, TESTED PROPERLY, AND ITS CONTROL.
#
#   THE CLAIM (author, 2026-09-02): after correcting for proliferation, OXPHOS
#   correlates with the apoptotic machinery and MYC does not.
#
# E10 found the raw version of half of this - the localisation split predicts
# the sign of a gene's OXPHOS correlation at 0.453 and its MYC correlation at
# 0.140 - but E10 adjusted for nothing. This script is the adjusted test, and
# more importantly it is the CONTROL that decides whether the adjusted result
# means anything at all.
#
# =============================================================================
# THE CONTROL IS THE POINT. WITHOUT IT THE CLAIM IS UNFALSIFIABLE
# =============================================================================
# There is a boring explanation for "MYC stops correlating with the machinery
# after proliferation adjustment", and it has nothing to do with apoptosis:
#
#   A MYC ACTIVITY SCORE IS PARTLY A PROLIFERATION SCORE (trap 3, 1.5% to 47.6%
#   entanglement across the panel). Partial out proliferation and it is
#   possible that what remains correlates with NOTHING. "MYC does not correlate
#   with the machinery" would then be a special case of "adjusted MYC does not
#   correlate with anything", which is a statement about the estimator, not
#   about apoptosis.
#
# So every set is measured against three others on the same axis, under the
# same adjustment, in the same samples:
#
#   MITOCHONDRIAL RIBOSOME (83 genes)   THE POSITIVE CONTROL. F6 says it is the
#       arm MYC tracks most strongly (0.590 raw). It shares no gene with the 44,
#       none with the OXPHOS arm, and none with the proliferation covariate. If
#       adjusted MYC still tracks these, adjusted MYC still carries signal and
#       the machinery's collapse is about the machinery.
#   OXPHOS SUBUNITS (89 genes)          A second positive control FOR THE MYC
#       AXIS ONLY. On the OXPHOS axis it is the score's own genes and its value
#       there is self-correlation - reported, marked, never read as a result.
#   EXPRESSION-MATCHED BACKGROUND       what a set of this size and expression
#       profile scores by construction. 2,000 draws over 20 ventiles.
#
# A claim of the form "X correlates with the axis" is a claim about X ABOVE this
# background, and that is what the z column reports.
#
# =============================================================================
# WHAT "CORRECTING FOR PROLIFERATION" MEANS, AND WHY BOTH WAYS ARE DONE
# =============================================================================
# E1 established that the two available corrections DISAGREE IN SIGN, so a
# result under one of them is not a result.
#
#   PARTIAL CORRELATION on a proliferation score. Applies IDENTICALLY to both
#       axes, which is what makes MYC and OXPHOS comparable at all - this is the
#       primary correction here for that reason alone. The covariate is
#       PROLIF_DISJOINT (318 genes), and it is genuinely clean for this
#       comparison: 0 genes shared with FELSHER__MITOSTRIP, 0 with the 89
#       OXPHOS-arm genes, 0 with the 83 mitoribosome genes. PROLIF_STD (327) is
#       carried beside it and is NOT clean - it shares 9 genes with the MYC
#       reference - so where the two disagree, PROLIF_DISJOINT is the reading.
#   REMOVING THE GENES from the estimator. FELSHER__PROLIFSTRIP (57 genes) and
#       FELSHER__BOTHSTRIP (52) are scored and read raw. This has NO OXPHOS
#       counterpart - the OXPHOS arm contains no E2F/G2M gene to remove - so it
#       is a one-sided check on the MYC half, not a symmetric comparison.
#
# TWO OF THE 44 ARE IN THE PROLIFERATION COVARIATE - `TP53` and `BIRC5`. Their
# adjusted values are partly adjustments for themselves. They are marked on
# every figure and every summary is reported with and without them.
#
# =============================================================================
# WHAT WOULD FALSIFY THE CLAIM, WRITTEN BEFORE THE ANSWER IS SEEN
# =============================================================================
#   IT SURVIVES if, after adjustment: the machinery's z against its matched null
#   stays large on OXPHOS and falls to about zero on MYC, WHILE the mitoribosome
#   control stays large on BOTH. That combination cannot be explained by the MYC
#   estimator losing its signal, because the control kept it.
#
#   IT FAILS, and the finding is about the estimator rather than about
#   apoptosis, if the mitoribosome control collapses on MYC too.
#
#   IT IS UNDECIDED if the two corrections disagree - in which case E1 applies
#   and neither is reportable without the other.
#
# EXPLORATORY. Nothing here is pre-registered. The claim was formed after seeing
# E10's unadjusted numbers, which is exactly why it needs a control and a
# falsifier rather than a p-value.
#
# SCALE: linear DESeq2-normalised at gene level, and every correlation is
# rank-based, so the log-versus-linear question does not arise (see the header
# of functions/correlation_engine.R). SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))
source(here::here("functions", "gene_matrix.R"))
source(here::here("functions", "strata.R"))

message("\nE11: does the apoptotic machinery still track OXPHOS, and MYC not,",
        "\n     once proliferation is taken out?\n", strrep("=", 78))

PATH_E11     <- file.path(DIR_RESULTS, "prolif_adjusted_machinery.rds")
PATH_E11_CSV <- file.path(DIR_TABLES,  "E11_gene_rho_by_adjustment.csv")

NULL_DRAWS <- 2000L
N_BINS     <- 20L
PROLIF_REF_COV <- "PROLIF_DISJOINT"   # the clean one - see the header

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

e08  <- readRDS(file.path(DIR_RESULTS, "strata_and_death_genes.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)
frames <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames

if (!"canonical" %in% names(e08) || nrow(e08$canonical) != 44L) {
  stop("results/strata_and_death_genes.rds does not carry the 44-gene ",
       "`canonical` table. Re-run E08.", call. = FALSE)
}

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)
tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
LT <- tcga_lin$mat[, ID_T, drop = FALSE]
LS <- scanb_lin$mat[, ID_S, drop = FALSE]
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

RES_T <- .symbol_resolver(rownames(LT), NULL)
RES_S <- .symbol_resolver(rownames(LS), sc$symbol_map)

# --- the gene sets -----------------------------------------------------------
CANON  <- sort(e08$canonical$gene)
MITORIBO <- sd_$arm_sets[["Mitochondrial ribosome"]]
OXARM    <- sd_$arm_sets[["OXPHOS subunits"]]
PROLIF_COV_GENES <- sd_$cov_sets[[PROLIF_REF_COV]]

SETS <- list(`apoptotic machinery (44)` = CANON,
             `mitoribosome arm (83)`    = MITORIBO,
             `OXPHOS subunits arm (89)` = OXARM)

# The overlap audit, printed rather than assumed. Anything non-zero here is a
# cell that is partly a correlation with itself.
overlap_audit <- tibble::tibble(
  set = names(SETS),
  n = lengths(SETS),
  in_prolif_covariate = vapply(SETS, function(g)
    length(intersect(g, PROLIF_COV_GENES)), integer(1)),
  in_MYC_ref = vapply(SETS, function(g)
    length(intersect(g, sd_$myc_sets[[MYC_REF]])), integer(1)),
  in_OXPHOS_arm = vapply(SETS, function(g)
    length(intersect(g, OXARM)), integer(1)))
message("\n   overlap audit - every non-zero cell is partial self-correlation:")
overlap_audit %>% as.data.frame() %>% print(row.names = FALSE)
SELF_ADJ <- intersect(CANON, PROLIF_COV_GENES)
message("   of the 44, adjusted for a covariate containing themselves: ",
        paste(SELF_ADJ, collapse = ", "))

# --- the axes and the covariates ---------------------------------------------
# Two axes only. E10 carried six; this script is one comparison and extra axes
# would just be more cells to not report.
.axis_mat <- function(gsva_new, arms_obj, ids) {
  m <- rbind(
    MYC    = as.numeric(gsva_new[MYC_REF, ids]),
    OXPHOS = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]))
  colnames(m) <- ids; m
}
# The removal-based MYC variants, read raw. No OXPHOS counterpart exists.
.strip_mat <- function(gsva_new, ids) {
  m <- rbind(
    `MYC PROLIFSTRIP` = as.numeric(gsva_new["FELSHER__PROLIFSTRIP", ids]),
    `MYC BOTHSTRIP`   = as.numeric(gsva_new["FELSHER__BOTHSTRIP", ids]))
  colnames(m) <- ids; m
}

COH <- list(
  TCGA = list(
    L = LT, res = RES_T, ids = ID_T,
    ax = .axis_mat(nw$tcga_gsva_new, mito, ID_T),
    strip = .strip_mat(nw$tcga_gsva_new, ID_T),
    cov = t(mito$gsva_cov[, ID_T, drop = FALSE])),
  `SCAN-B` = list(
    L = LS, res = RES_S, ids = ID_S,
    ax = .axis_mat(sc$gsva_new, sc, ID_S),
    strip = .strip_mat(sc$gsva_new, ID_S),
    cov = t(sc$gsva_cov[, ID_S, drop = FALSE])))
for (coh in names(COH)) {
  stopifnot(identical(rownames(COH[[coh]]$cov), COH[[coh]]$ids),
            all(c("PROLIF_STD", "PROLIF_DISJOINT") %in%
                  colnames(COH[[coh]]$cov)))
}
message("   axes and covariates built for both cohorts")

# TCGA only: purity and leukocyte fraction, for the secondary check. SCAN-B has
# no purity estimate at all (trap 2) and it is never imputed.
PL <- frames %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::select(sample_id, purity, leuko)
PL <- PL[match(ID_T, PL$sample_id), ]
PURITY_OK <- stats::complete.cases(PL[, c("purity", "leuko")])
message("   TCGA samples with purity AND leukocyte fraction: ", sum(PURITY_OK),
        " of ", length(ID_T))

# =============================================================================
# 2. Every gene of the matrix, under every axis and every adjustment
# =============================================================================
# ONE FULL-MATRIX PASS PER (cohort, axis, adjustment), not one per gene set.
# Everything below - the three sets, the background, the matched null - is then
# a subset of the same vector, so a set and its null are guaranteed to have been
# computed the same way. Building them separately is how a set and its
# comparator quietly stop being comparable.
message("\n2. per-gene correlations over the whole matrix")

ADJUSTMENTS <- list(
  raw                  = NULL,
  `adj. PROLIF_DISJOINT` = "PROLIF_DISJOINT",
  `adj. PROLIF_STD`      = "PROLIF_STD")

per_gene <- list()
for (coh in names(COH)) {
  C <- COH[[coh]]
  for (ax in rownames(C$ax)) {
    for (adj in names(ADJUSTMENTS)) {
      cv <- ADJUSTMENTS[[adj]]
      Z  <- if (is.null(cv)) NULL else C$cov[, cv, drop = FALSE]
      per_gene[[paste(coh, ax, adj, sep = "|")]] <-
        .per_gene_rho(C$L, C$ax[ax, ], cov = Z)
      message("   ", coh, " ", ax, " ", adj, " - done")
    }
  }
  # The removal-based MYC variants, raw only.
  for (ax in rownames(C$strip)) {
    per_gene[[paste(coh, ax, "raw", sep = "|")]] <-
      .per_gene_rho(C$L, C$strip[ax, ], cov = NULL)
    message("   ", coh, " ", ax, " raw - done")
  }
}

# =============================================================================
# 3. The expression-matched null, and what a set scores above it
# =============================================================================
# THREE STATISTICS, BECAUSE "CORRELATES WITH" IS AMBIGUOUS AND THE THREE READ
# DIFFERENTLY FOR A SET LIKE THIS ONE:
#
#   mean |rho|   "do these genes track this axis at all". The blunt version of
#                the claim.
#   mean rho     the SIGNED mean. The machinery is roughly half positive and
#                half negative, so this is near zero however strong the
#                individual correlations are - it is here to make that visible,
#                not because it is the statistic of interest.
#   SD of rho    "are these genes SPREAD along the axis" - the structural
#                version, and the one that matches what E10 fig1 shows.
#
# Each is compared with the SAME statistic on 2,000 expression-matched draws of
# the same size, so a set is measured against what a set of its size and
# expression profile scores by construction. A gene with no variance comes back
# NA from the engine and is dropped from the set AND the background alike, or
# the draws are biased - that is what `keep` is.
message("\n3. expression-matched nulls (", NULL_DRAWS, " draws, ", N_BINS,
        " ventiles)")
set.seed(PROJECT_SEED)

.set_stats <- function(v) c(mean_abs = mean(abs(v)), mean_signed = mean(v),
                            sd = stats::sd(v))
.z <- function(obs, nd) (obs - mean(nd)) / stats::sd(nd)

null_tests <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))
  B <- .expression_bins(C$L, keep, n_bins = N_BINS)
  gr <- lapply(SETS, function(g) .gene_rows(g, C$L, C$res))
  idx <- lapply(gr, function(x) {
    i <- match(rownames(x$mat), rownames(C$L)); intersect(i, keep) })
  draws <- lapply(idx, function(i)
    replicate(NULL_DRAWS, .matched_draw(i, B), simplify = FALSE))

  dplyr::bind_rows(lapply(names(per_gene)[startsWith(names(per_gene),
                                                     paste0(coh, "|"))],
    function(key) {
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]
      v <- per_gene[[key]]
      dplyr::bind_rows(lapply(names(SETS), function(sn) {
        obs <- .set_stats(v[idx[[sn]]])
        nd  <- vapply(draws[[sn]], function(d) .set_stats(v[d]), numeric(3))
        tibble::tibble(
          cohort = coh, axis = parts[2], adjustment = parts[3], set = sn,
          n_genes = length(idx[[sn]]),
          mean_abs_rho = obs[["mean_abs"]],
          null_mean_abs = mean(nd["mean_abs", ]),
          z_mean_abs = .z(obs[["mean_abs"]], nd["mean_abs", ]),
          mean_signed_rho = obs[["mean_signed"]],
          z_mean_signed = .z(obs[["mean_signed"]], nd["mean_signed", ]),
          sd_rho = obs[["sd"]],
          null_sd = mean(nd["sd", ]),
          z_sd = .z(obs[["sd"]], nd["sd", ]))
      }))
    }))
}))

message("\n   MEAN |rho| AND ITS EXPRESSION-MATCHED NULL:")
null_tests %>%
  dplyr::filter(adjustment %in% c("raw", "adj. PROLIF_DISJOINT"),
                axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(cohort, axis, adjustment, set, mean_abs_rho, null_mean_abs,
                z_mean_abs, sd_rho, null_sd, z_sd) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  dplyr::arrange(axis, set, cohort, adjustment) %>% as.data.frame() %>%
  print(row.names = FALSE)
message("\n   READ THE NULL COLUMN, NOT ONLY THE OBSERVED ONE. An OXPHOS score",
        " correlates\n   with a large fraction of the transcriptome, so a set",
        " can carry a substantial\n   mean |rho| and still be doing nothing a",
        " random set of the same size and\n   expression profile would not do.",
        " That is what z is for.")
message("\n   THE CONTROL. The claim is about apoptosis only if adjusted MYC",
        " still tracks\n   the MITORIBOSOME - 83 genes sharing none with the",
        " 44, none with the OXPHOS\n   arm and none with the covariate. If that",
        " collapses too, the finding is that\n   the adjustment emptied the",
        " estimator.")

message("\n   the other correction - removing the genes from the estimator:")
null_tests %>%
  dplyr::filter(axis %in% c("MYC PROLIFSTRIP", "MYC BOTHSTRIP")) %>%
  dplyr::select(cohort, axis, set, mean_abs_rho, z_mean_abs, sd_rho, z_sd) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  dplyr::arrange(set, cohort) %>% as.data.frame() %>% print(row.names = FALSE)
message("   E1: the two corrections disagree in sign elsewhere. If they",
        " disagree here,\n   neither is reportable without the other.")

# =============================================================================
# 3.1 THE HARDER NULL: is the localisation split specific to apoptosis?
# =============================================================================
# E10 R1 reported that among the 44, whether a gene is in MitoCarta predicts the
# sign of its OXPHOS correlation at 0.453, replicated under an independent
# localisation source. THAT NUMBER HAS AN OBVIOUS COMPETING EXPLANATION AND E10
# DID NOT TEST IT: mitochondrial genes correlate with an OXPHOS score more than
# non-mitochondrial genes do, for reasons that have nothing to do with death.
# Any 20-mitochondrial-plus-24-other set would then show the same split.
#
# So the null here holds the COMPOSITION fixed and matches expression WITHIN
# each half: 20 genes drawn from MitoCarta and 24 from outside it, each half
# expression-matched to the corresponding half of the machinery. The question
# it answers is exactly the right one - given 20 mitochondrial and 24
# non-mitochondrial genes of this expression profile, how large a split arises
# by construction?
#
# IF z IS NEAR ZERO, S6 IS GENERIC MITOCHONDRIAL BIOLOGY and E10 R1 must be
# restated. If it is large, the machinery's split is steeper than composition
# alone explains. Either answer is worth having; this is written before it.
message("\n3.1 is the localisation split specific to the machinery?")

# TWO NULL POOLS, because the first one can be accused of being unfair to the
# machinery. MitoCarta's 1,136 genes are dominated by the respiratory chain and
# the mitoribosome, which correlate with an OXPHOS score by definition, whereas
# the machinery's 20 mitochondrial genes are BCL2-family members and caspases.
# The stricter pool removes every OXPHOS and mitoribosome gene, so the null is
# drawn from mitochondrial proteins that are NOT part of the measured programme.
# If the observed split sits inside THAT null too, the conclusion is firm.
MC_ALL <- sd_$strip_refs$MITOCARTA_ALL
MC_STRICT <- setdiff(MC_ALL, unique(c(sd_$arm_sets[["OXPHOS umbrella"]],
                                      sd_$arm_sets[["OXPHOS subunits"]],
                                      sd_$arm_sets[["Mitochondrial ribosome"]],
                                      sd_$arm_sets[["mtDNA-encoded OXPHOS"]])))
# Short names because these become facet strips, and a clipped strip label is
# a figure that does not say what it is showing.
POOLS <- list(`null pool: all MitoCarta` = MC_ALL,
              `null pool: no OXPHOS/mitoribo` = MC_STRICT)
message("   null pools: all MitoCarta ", length(MC_ALL),
        " genes | strict ", length(MC_STRICT))

split_null <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))
  gr <- .gene_rows(CANON, C$L, C$res)
  own <- match(rownames(gr$mat), rownames(C$L))
  own <- own[!is.na(own)]

  dplyr::bind_rows(lapply(names(POOLS), function(pn) {
    is_mc <- rownames(C$L) %in% POOLS[[pn]]
    keep_mc  <- intersect(keep, which(is_mc))
    keep_non <- intersect(keep, which(!is_mc))
    # Ten bins for the mitochondrial pool - it is far smaller than the
    # background and 20 ventiles would leave too few genes per bin to draw
    # from without replacement.
    B_mc  <- .expression_bins(C$L, keep_mc,  n_bins = 10L)
    B_non <- .expression_bins(C$L, keep_non, n_bins = N_BINS)

    # The OBSERVED split always uses MitoCarta membership as E10 defined it;
    # only the null POOL changes between the two rows.
    i_mc  <- intersect(own, intersect(keep, which(rownames(C$L) %in% MC_ALL)))
    i_non <- setdiff(intersect(own, keep), i_mc)
    # Draw the mitochondrial half from this pool and the rest from outside it,
    # each expression-matched, and never draw one of the 44 themselves.
    d_mc  <- intersect(i_mc,  keep_mc)
    d_non <- intersect(i_non, keep_non)
    if (length(d_mc) < 3L) return(NULL)
    dr <- replicate(NULL_DRAWS,
                    c(.matched_draw(d_mc,  B_mc,  exclude = own),
                      .matched_draw(d_non, B_non, exclude = own)),
                    simplify = FALSE)
    lab_draw <- c(rep(1, length(d_mc)), rep(0, length(d_non)))
    lab_obs  <- c(rep(1, length(i_mc)), rep(0, length(i_non)))

    dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
      dplyr::bind_rows(lapply(c("raw", "adj. PROLIF_DISJOINT"), function(adj) {
        v <- per_gene[[paste(coh, ax, adj, sep = "|")]]
        obs <- stats::cor(v[c(i_mc, i_non)], lab_obs, method = "spearman")
        nd  <- vapply(dr, function(d)
          stats::cor(v[d], lab_draw, method = "spearman"), numeric(1))
        tibble::tibble(cohort = coh, axis = ax, adjustment = adj, pool = pn,
                       n_mito = length(i_mc), n_other = length(i_non),
                       n_drawn_mito = length(d_mc),
                       observed_split = obs, null_mean = mean(nd),
                       null_sd = stats::sd(nd), z = .z(obs, nd),
                       pct_of_draws_below = mean(nd < obs))
      }))))
  }))
}))
message("\n   the MitoCarta split, against a COMPOSITION-matched null:")
split_null %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(cohort, axis, adjustment, pool, observed_split, null_mean,
                null_sd, z) %>%
  dplyr::arrange(axis, pool, cohort, adjustment) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   `null_mean` IS THE NUMBER TO READ FIRST. It is how large the",
        " split is for a\n   set with the machinery's mitochondrial",
        " composition and expression profile and\n   nothing else in common",
        " with it. E10 R1's 0.453 is only a statement about\n   apoptosis to",
        " the extent that it exceeds that.")

# =============================================================================
# 3.2 THE STRICTER NULL AGAIN: match the SUB-COMPARTMENT too
# =============================================================================
# Section 3.1 matched MitoCarta MEMBERSHIP. That is not the whole of the
# composition: E10 found a ladder inside the organelle - not in MitoCarta
# -0.100, outer membrane +0.085, intermembrane space +0.115, inner membrane
# +0.443 - so 20 mitochondrial genes drawn without regard to compartment are
# not the same object as the machinery's 20, which are 13 MOM, 5 IMS and 2 MIM.
#
# This null matches that distribution exactly, expression-matching within each
# compartment separately, still drawing from the strict pool with every OXPHOS
# and mitoribosome gene removed. It is the hardest test of E10 R1 available
# without a new annotation source.
#
# IF THE OBSERVED SPLIT SITS INSIDE THIS NULL, the whole of E10 R1's 0.453 is
# expression plus sub-mitochondrial location and none of it is apoptosis.
#
# THE POOLS ARE SMALL AND THE BINNING IS COARSE BECAUSE OF IT. IMS holds 39
# usable genes, so it gets 2 expression bins rather than 20; MOM 110 gets 4 and
# MIM 213 gets 8. That is weaker expression matching than section 3, and it is a
# limitation of the annotation, not a choice.
message("\n3.2 the same null, matching sub-mitochondrial compartment as well")

mitocarta_sheet <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
SUBMITO <- stats::setNames(
  mitocarta_sheet[["MitoCarta3.0_SubMitoLocalization"]], mitocarta_sheet$Symbol)

# The machinery's own composition, which is what the draws must reproduce.
#
# NOTE THE unname() AND THE names() PUT BACK. Subsetting a named vector by a
# name it does not carry returns an element named <NA>, not one named by the
# key that was asked for. Leaving that in place makes the 24 non-MitoCarta
# genes unnamed, the later lookup by gene symbol returns NA for all of them,
# split() drops NA groups silently, and the label vector comes out constant -
# which yields an NA correlation and no error anywhere.
canon_sub <- unname(SUBMITO[CANON])
canon_sub[is.na(canon_sub)] <- "(not in MitoCarta)"
names(canon_sub) <- CANON
stopifnot(!anyNA(names(canon_sub)), !anyNA(canon_sub),
          length(canon_sub) == length(CANON))
message("   the 44 by compartment: ",
        paste(names(table(canon_sub)), table(canon_sub), sep = " ",
              collapse = " | "))

# Bins scale with the pool. A 39-gene pool cannot support 20 ventiles and
# .matched_draw would stop rather than quietly draw a worse match.
.nbins <- function(n) max(2L, min(N_BINS, n %/% 25L))

# Build one composition-matched draw set for a cohort, reusable across axes and
# adjustments because the draws do not depend on which correlation is read.
.compartment_draws <- function(C, keep, own, want, pool_genes) {
  parts <- list()
  for (cmp in names(want)) {
    i_want <- want[[cmp]]
    if (!length(i_want)) next
    pool <- if (cmp == "(not in MitoCarta)") {
      setdiff(keep, which(rownames(C$L) %in% MC_ALL))
    } else {
      intersect(keep, which(rownames(C$L) %in%
                              intersect(pool_genes,
                                        names(SUBMITO)[!is.na(SUBMITO) &
                                                         SUBMITO == cmp])))
    }
    # The machinery's own genes must be IN the binning (or their bin is NA and
    # the draw silently returns nothing) and OUT of the selection. CYCS is the
    # gene that makes this necessary: it is IMS but sits in the OXPHOS arm, so
    # the strict pool excludes it.
    pool <- union(pool, i_want)
    B <- .expression_bins(C$L, pool, n_bins = .nbins(length(pool)))
    parts[[cmp]] <- list(want = i_want, B = B)
  }
  parts
}

split_null_submito <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))
  gr <- .gene_rows(CANON, C$L, C$res)
  own <- match(rownames(gr$mat), rownames(C$L)); own <- own[!is.na(own)]
  sub_of <- canon_sub[rownames(gr$mat)]
  want <- lapply(split(own, sub_of), function(v) intersect(v, keep))
  want <- want[lengths(want) > 0L]

  parts <- .compartment_draws(C, keep, own, want, MC_STRICT)
  dr <- replicate(NULL_DRAWS,
                  unlist(lapply(parts, function(pp)
                    .matched_draw(pp$want, pp$B, exclude = own)),
                    use.names = FALSE),
                  simplify = FALSE)
  # Labels follow the same compartment order in the observed set and the draws.
  obs_idx <- unlist(lapply(parts, function(pp) pp$want), use.names = FALSE)
  lab <- as.numeric(rep(names(parts) != "(not in MitoCarta)",
                        vapply(parts, function(pp) length(pp$want), integer(1))))
  # A constant label gives an NA correlation and no error. Stop instead.
  if (length(unique(lab)) < 2L || length(obs_idx) != length(lab)) {
    stop("the composition labels are degenerate in ", coh, " - ",
         length(obs_idx), " genes, ", length(unique(lab)), " label value(s). ",
         "The sub-compartment lookup has lost genes.", call. = FALSE)
  }

  dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
    dplyr::bind_rows(lapply(c("raw", "adj. PROLIF_DISJOINT"), function(adj) {
      v <- per_gene[[paste(coh, ax, adj, sep = "|")]]
      obs <- stats::cor(v[obs_idx], lab, method = "spearman")
      nd  <- vapply(dr, function(d) stats::cor(v[d], lab, method = "spearman"),
                    numeric(1))
      tibble::tibble(cohort = coh, axis = ax, adjustment = adj,
                     pool = "sub-compartment matched",
                     observed_split = obs, null_mean = mean(nd),
                     null_sd = stats::sd(nd), z = .z(obs, nd),
                     pct_of_draws_below = mean(nd < obs))
    }))))
}))
message("\n   the split against a SUB-COMPARTMENT-matched null:")
split_null_submito %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   Compare with section 3.1. A null that RISES when the compartment",
        " is matched\n   means the compartment ladder was carrying the split;",
        " a z that falls towards\n   zero means E10 R1 is composition and",
        " nothing else.")

# =============================================================================
# 4. The 44 genes one by one, on both axes, before and after
# =============================================================================
message("\n4. the 44 genes, gene by gene")

gene_tab <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  gr <- .gene_rows(CANON, C$L, C$res)
  i  <- match(rownames(gr$mat), rownames(C$L))
  dplyr::bind_rows(lapply(names(per_gene)[startsWith(names(per_gene),
                                                     paste0(coh, "|"))],
    function(key) {
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]
      tibble::tibble(cohort = coh, axis = parts[2], adjustment = parts[3],
                     gene = rownames(gr$mat), rho = per_gene[[key]][i])
    }))
})) %>%
  dplyr::left_join(
    e08$canonical %>% dplyr::select(gene, effect, cdc_module = module,
                                    cdc_acts_at_mito = acts_at_mito),
    by = "gene") %>%
  dplyr::mutate(
    mitocarta       = gene %in% sd_$strip_refs$MITOCARTA_ALL,
    in_prolif_cov   = gene %in% PROLIF_COV_GENES,
    in_oxphos_arm   = gene %in% OXARM)

wide <- gene_tab %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"),
                adjustment %in% c("raw", "adj. PROLIF_DISJOINT")) %>%
  dplyr::select(cohort, gene, axis, adjustment, rho, mitocarta, in_prolif_cov,
                effect) %>%
  tidyr::pivot_wider(names_from = c(axis, adjustment), values_from = rho)

message("\n   spread of the 44 per-gene correlations (SD), which is the whole",
        " claim in one number:")
spread <- gene_tab %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::group_by(cohort, axis, adjustment) %>%
  dplyr::summarise(sd_all = stats::sd(rho),
                   sd_excl_self = stats::sd(rho[!in_prolif_cov]),
                   frac_above_0.2 = mean(abs(rho) > 0.2),
                   .groups = "drop")
spread %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   `sd_excl_self` drops TP53 and BIRC5, which are IN the proliferation",
        "\n   covariate and are therefore partly adjusted for themselves.")

message("\n   does the LOCALISATION SPLIT (S6) survive the adjustment?")
s6_adj <- gene_tab %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"), !in_prolif_cov) %>%
  dplyr::group_by(cohort, axis, adjustment) %>%
  dplyr::summarise(
    vs_mitocarta = stats::cor(rho, as.numeric(mitocarta), method = "spearman"),
    vs_pro_death = stats::cor(rho, as.numeric(effect == "pro-death"),
                              method = "spearman"),
    median_mito = stats::median(rho[mitocarta]),
    median_nonmito = stats::median(rho[!mitocarta]), .groups = "drop")
s6_adj %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   cross-cohort replication of the 44 per-gene values:")
replication <- gene_tab %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS")) %>%
  dplyr::select(cohort, axis, adjustment, gene, rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
  dplyr::group_by(axis, adjustment) %>%
  dplyr::summarise(spearman_TCGA_vs_SCANB =
                     stats::cor(TCGA, `SCAN-B`, method = "spearman"),
                   .groups = "drop")
replication %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   A correlation structure that VANISHES should also stop replicating.",
        "\n   If an axis keeps a high replication with a collapsed spread, the",
        " signal moved\n   rather than disappeared and this reading is wrong.")

message("\n   the genes that move most when proliferation is taken out:")
gene_tab %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"),
                adjustment %in% c("raw", "adj. PROLIF_DISJOINT")) %>%
  dplyr::select(cohort, axis, gene, adjustment, rho, in_prolif_cov) %>%
  tidyr::pivot_wider(names_from = adjustment, values_from = rho) %>%
  dplyr::mutate(shift = `adj. PROLIF_DISJOINT` - raw) %>%
  dplyr::arrange(dplyr::desc(abs(shift))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 4.1 the TCGA-only purity check -----------------------------------------
# Trap 2: breast is the worst TCGA tissue for infiltrate, adipose is
# OXPHOS-high, and immune cells carry their own BCL2-family profile. SCAN-B has
# no purity estimate, so this is one cohort only and cannot replicate.
message("\n4.1 TCGA only: proliferation PLUS purity and leukocyte fraction")
C <- COH$TCGA
gr <- .gene_rows(CANON, C$L, C$res)
ids_p <- C$ids[PURITY_OK]
cov_p <- cbind(C$cov[PURITY_OK, PROLIF_REF_COV, drop = FALSE],
               as.matrix(PL[PURITY_OK, c("purity", "leuko")]))
purity_tab <- dplyr::bind_rows(lapply(rownames(C$ax), function(ax)
  .atlas_block(C$ax[ax, , drop = FALSE], gr$mat, ids_p, cov = cov_p,
               min_n = 100L) %>%
    dplyr::rename(axis = myc_estimator, gene = measure) %>%
    dplyr::mutate(adjustment = "adj. prolif + purity + leuko"))) %>%
  dplyr::mutate(mitocarta = gene %in% sd_$strip_refs$MITOCARTA_ALL)
purity_tab %>% dplyr::group_by(axis) %>%
  dplyr::summarise(n_samples = dplyr::first(n), sd_rho = stats::sd(rho),
                   vs_mitocarta = stats::cor(rho, as.numeric(mitocarta),
                                             method = "spearman"),
                   .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   TCGA only, n = ", sum(PURITY_OK), ". SCAN-B has no purity estimate",
        " and this row\n   cannot be replicated - trap 2.")

# =============================================================================
# 4.2 THE 44 INSIDE THE LUMINAL AND BASAL COMPARTMENTS
# =============================================================================
# P1 - the claim that the OXPHOS/MYC contrast survives proliferation - was
# measured on all samples. E10 fig9 showed that for the BCL2-family ratios the
# POOLED MYC value sits outside the range of both its compartments for 27 of 39
# ratios in TCGA and 26 in SCAN-B, which is the D3/S1 signature of a
# between-subtype effect read as a within-subtype one. If the same is true of
# the 44, then P1's MYC column is describing the difference between luminal and
# basal tumours rather than anything inside either.
#
# THE COMPARISON IS UNEVEN AND SAYING SO IS PART OF THE RESULT. Luminal is 696
# TCGA and 2,436 SCAN-B samples; Basal is 171 and 317. A 171-sample stratum
# gives a 95% interval about +/- 0.15 wide on a single gene's rho, so a
# Luminal-minus-Basal difference under roughly 0.2 is not separable from
# sampling in TCGA. The nulls absorb this - a null draw in Basal is exactly as
# noisy as the observed set - but the per-gene numbers do not.
message("\n4.2 the 44 by compartment")

STRATA_E11 <- c("Luminal", "Basal")
STR <- list(TCGA = .build_strata(frames, "TCGA", ID_T),
            `SCAN-B` = .build_strata(frames, "SCAN-B", ID_S))
tibble::tibble(stratum = c("all", STRATA_E11),
               TCGA = lengths(STR$TCGA[c("all", STRATA_E11)]),
               `SCAN-B` = lengths(STR$`SCAN-B`[c("all", STRATA_E11)])) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Per-gene correlations again, but within each compartment. PROLIF_STD is
# dropped here - section 4 showed it and PROLIF_DISJOINT agree to the third
# decimal - and so are the stripped estimators, which are a check on the
# correction rather than a question about subtype.
message("\n   per-gene correlations within each compartment")
per_gene_str <- list()
for (coh in names(COH)) {
  C <- COH[[coh]]
  for (st in STRATA_E11) {
    ids <- STR[[coh]][[st]]
    Lsub <- C$L[, ids, drop = FALSE]
    for (ax in c("MYC", "OXPHOS")) {
      for (adj in c("raw", "adj. PROLIF_DISJOINT")) {
        cv <- ADJUSTMENTS[[adj]]
        Z  <- if (is.null(cv)) NULL else C$cov[ids, cv, drop = FALSE]
        per_gene_str[[paste(coh, st, ax, adj, sep = "|")]] <-
          .per_gene_rho(Lsub, C$ax[ax, ids], cov = Z)
      }
    }
    rm(Lsub); invisible(gc(verbose = FALSE))
    message("   ", coh, " ", st, " (n = ", length(ids), ") - done")
  }
}

# The two summaries, and the composition-matched null for the split. Without
# the null this section would repeat exactly the error section 3.1 found.
compartment <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  gr <- .gene_rows(CANON, C$L, C$res)
  own <- match(rownames(gr$mat), rownames(C$L)); own <- own[!is.na(own)]
  sub_of <- canon_sub[rownames(gr$mat)]
  is_mc <- rownames(C$L)[own] %in% MC_ALL

  dplyr::bind_rows(lapply(STRATA_E11, function(st) {
    keep <- which(!is.na(per_gene_str[[paste(coh, st, "MYC", "raw",
                                             sep = "|")]]))
    want <- lapply(split(own, sub_of), function(v) intersect(v, keep))
    want <- want[lengths(want) > 0L]
    parts <- .compartment_draws(C, keep, own, want, MC_STRICT)
    dr <- replicate(NULL_DRAWS,
                    unlist(lapply(parts, function(pp)
                      .matched_draw(pp$want, pp$B, exclude = own)),
                      use.names = FALSE),
                    simplify = FALSE)
    obs_idx <- unlist(lapply(parts, function(pp) pp$want), use.names = FALSE)
    lab <- as.numeric(rep(names(parts) != "(not in MitoCarta)",
                          vapply(parts, function(pp) length(pp$want),
                                 integer(1))))
    if (length(unique(lab)) < 2L || length(obs_idx) != length(lab)) {
      stop("degenerate composition labels in ", coh, " ", st, call. = FALSE)
    }

    dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
      dplyr::bind_rows(lapply(c("raw", "adj. PROLIF_DISJOINT"), function(adj) {
        v <- per_gene_str[[paste(coh, st, ax, adj, sep = "|")]]
        obs <- stats::cor(v[obs_idx], lab, method = "spearman")
        nd  <- vapply(dr, function(d) stats::cor(v[d], lab, method = "spearman"),
                      numeric(1))
        tibble::tibble(
          cohort = coh, stratum = st, n_samples = length(STR[[coh]][[st]]),
          axis = ax, adjustment = adj,
          sd_rho = stats::sd(v[own]),
          median_mito = stats::median(v[own][is_mc]),
          median_nonmito = stats::median(v[own][!is_mc]),
          observed_split = obs, null_mean = mean(nd), null_sd = stats::sd(nd),
          z = .z(obs, nd))
      }))))
  }))
}))

# The pooled row, computed the same way, so the three can be compared directly.
pooled_row <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  gr <- .gene_rows(CANON, C$L, C$res)
  own <- match(rownames(gr$mat), rownames(C$L)); own <- own[!is.na(own)]
  is_mc <- rownames(C$L)[own] %in% MC_ALL
  dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
    dplyr::bind_rows(lapply(c("raw", "adj. PROLIF_DISJOINT"), function(adj) {
      v <- per_gene[[paste(coh, ax, adj, sep = "|")]]
      tibble::tibble(cohort = coh, stratum = "all",
                     n_samples = length(COH[[coh]]$ids), axis = ax,
                     adjustment = adj, sd_rho = stats::sd(v[own]),
                     median_mito = stats::median(v[own][is_mc]),
                     median_nonmito = stats::median(v[own][!is_mc]),
                     observed_split = stats::cor(v[own], as.numeric(is_mc),
                                                 method = "spearman"),
                     null_mean = NA_real_, null_sd = NA_real_, z = NA_real_)
    }))))
}))
compartment <- dplyr::bind_rows(pooled_row, compartment) %>%
  dplyr::mutate(stratum = factor(stratum, levels = c("all", STRATA_E11)))

message("\n   spread and split of the 44, by compartment (adjusted):")
compartment %>% dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::select(cohort, stratum, n_samples, axis, sd_rho, observed_split,
                null_mean, z) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(cohort, axis, stratum) %>% as.data.frame() %>%
  print(row.names = FALSE)

message("\n   IS THE POOLED VALUE OUTSIDE BOTH ITS COMPARTMENTS?")
message("   (that is the D3/S1 signature - a between-subtype effect read as a",
        " within one)")
outside <- compartment %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::select(cohort, axis, stratum, sd_rho, observed_split) %>%
  tidyr::pivot_wider(names_from = stratum,
                     values_from = c(sd_rho, observed_split)) %>%
  dplyr::mutate(
    sd_pooled_outside = sd_rho_all < pmin(sd_rho_Luminal, sd_rho_Basal) |
                        sd_rho_all > pmax(sd_rho_Luminal, sd_rho_Basal),
    split_pooled_outside =
      observed_split_all < pmin(observed_split_Luminal, observed_split_Basal) |
      observed_split_all > pmax(observed_split_Luminal, observed_split_Basal))
outside %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   and the per-gene values, Luminal minus Basal, both cohorts",
        " agreeing on sign:")
gene_compartment <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  gr <- .gene_rows(CANON, C$L, C$res)
  own <- match(rownames(gr$mat), rownames(C$L))
  dplyr::bind_rows(lapply(STRATA_E11, function(st)
    dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
      tibble::tibble(cohort = coh, stratum = st, axis = ax,
                     gene = rownames(gr$mat),
                     rho = per_gene_str[[paste(coh, st, ax,
                                               "adj. PROLIF_DISJOINT",
                                               sep = "|")]][own])))))
})) %>%
  dplyr::left_join(dplyr::distinct(dplyr::select(gene_tab, gene, mitocarta,
                                                 effect)), by = "gene")
gene_compartment %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho) %>%
  dplyr::mutate(diff = Luminal - Basal) %>%
  dplyr::select(cohort, axis, gene, mitocarta, diff) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = diff) %>%
  dplyr::filter(sign(TCGA) == sign(`SCAN-B`)) %>%
  dplyr::arrange(dplyr::desc(abs(TCGA) + abs(`SCAN-B`))) %>% utils::head(12) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   Basal is 171 TCGA samples. Nothing under about 0.2 here is",
        " separable from\n   sampling in that cohort - read the column, not",
        " the row.")

# =============================================================================
# 4.3 THE TEST THE CLAIM ACTUALLY NEEDS: condition each axis on the other
# =============================================================================
# Everything up to here compares two axes SEPARATELY. That is not the claim.
# "The machinery tracks OXPHOS rather than MYC" is a claim about which axis
# carries the association once the other is allowed for - and the two axes are
# themselves correlated (Spearman 0.39 TCGA / 0.43 SCAN-B raw, 0.32 / 0.26 after
# proliferation), so each one's apparent ordering could be entirely inherited
# from the other.
#
# The two columns are in fact near-copies: across the 44 genes the MYC and
# OXPHOS per-gene correlations rank together at 0.61 to 0.73. So the separate
# comparison CANNOT distinguish "OXPHOS orders them and MYC does not" from
# "both order them and OXPHOS is simply the better-measured of two versions of
# the same thing".
#
# FOUR MODELS, AND THE COMPARISON IS BETWEEN ROWS 2 AND 4:
#   OXPHOS | prolif           the section 4 value
#   OXPHOS | prolif + MYC     does OXPHOS survive conditioning on MYC?
#   MYC    | prolif           the section 4 value
#   MYC    | prolif + OXPHOS  does MYC survive conditioning on OXPHOS?
#
# IF BOTH SURVIVE, the two axes carry separable information and neither claim is
# available. IF ONLY OXPHOS SURVIVES, MYC's ordering was inherited and the
# author's statement holds. IF ONLY MYC SURVIVES, the statement is backwards.
# Written before the numbers were looked at.
#
# The localisation split carries its sub-compartment-matched null in every row,
# because a split that survives conditioning is still only a statement about
# apoptosis to the extent that it exceeds what composition gives (3.2).
message("\n4.3 conditioning each axis on the other")

COND_MODELS <- list(
  `OXPHOS | prolif`       = list(y = "OXPHOS", extra = character(0)),
  `OXPHOS | prolif + MYC` = list(y = "OXPHOS", extra = "MYC"),
  `MYC | prolif`          = list(y = "MYC",    extra = character(0)),
  `MYC | prolif + OXPHOS` = list(y = "MYC",    extra = "OXPHOS"))

conditioned <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))
  gr <- .gene_rows(CANON, C$L, C$res)
  own <- match(rownames(gr$mat), rownames(C$L)); own <- own[!is.na(own)]
  is_mc <- rownames(C$L)[own] %in% MC_ALL
  sub_of <- canon_sub[rownames(gr$mat)]
  want <- lapply(split(own, sub_of), function(v) intersect(v, keep))
  want <- want[lengths(want) > 0L]
  parts <- .compartment_draws(C, keep, own, want, MC_STRICT)
  dr <- replicate(NULL_DRAWS,
                  unlist(lapply(parts, function(pp)
                    .matched_draw(pp$want, pp$B, exclude = own)),
                    use.names = FALSE), simplify = FALSE)
  obs_idx <- unlist(lapply(parts, function(pp) pp$want), use.names = FALSE)
  lab <- as.numeric(rep(names(parts) != "(not in MitoCarta)",
                        vapply(parts, function(pp) length(pp$want), integer(1))))

  dplyr::bind_rows(lapply(names(COND_MODELS), function(mn) {
    m <- COND_MODELS[[mn]]
    Z <- cbind(C$cov[, PROLIF_REF_COV, drop = FALSE])
    for (e in m$extra) Z <- cbind(Z, stats::setNames(C$ax[e, ], NULL))
    colnames(Z) <- c(PROLIF_REF_COV, m$extra)
    v <- .per_gene_rho(C$L, C$ax[m$y, ], cov = Z)
    obs <- stats::cor(v[obs_idx], lab, method = "spearman")
    nd  <- vapply(dr, function(d) stats::cor(v[d], lab, method = "spearman"),
                  numeric(1))
    tibble::tibble(
      cohort = coh, model = mn, axis = m$y,
      conditioned_on = if (length(m$extra)) m$extra else "-",
      sd_rho = stats::sd(v[own]), frac_gt_0.2 = mean(abs(v[own]) > 0.2),
      median_mito = stats::median(v[own][is_mc]),
      median_nonmito = stats::median(v[own][!is_mc]),
      split = obs, null_mean = mean(nd), null_sd = stats::sd(nd),
      z_split = .z(obs, nd))
  }))
}))

message("\n   the 44 under each model - read rows 2 and 4:")
conditioned %>%
  dplyr::select(cohort, model, sd_rho, frac_gt_0.2, split, null_mean, z_split) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   A LOCALISATION SPLIT THAT GOES TO ZERO UNDER CONDITIONING MEANS",
        " THAT AXIS\n   NEVER ORDERED THE GENES - it was reading the other one.",
        " A split that is\n   UNCHANGED means the axis carries the ordering",
        " itself.")

# =============================================================================
# 5. Figures
# =============================================================================
message("\n5. figures")

COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
theme_e11 <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey92",
                                                          colour = NA),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8, colour = "grey40",
                                                      hjust = 0))
.save <- function(p, nm, w, h) {
  for (ext in c("png", "pdf"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0(nm, ".", ext)), p,
                    width = w, height = h, dpi = 300)
  message("   ", nm); invisible(p)
}

# --- FIG 1: THE PLANE. One picture, and it is the whole claim ----------------
# Each point is one of the 44 genes. x is what it does with MYC, y is what it
# does with OXPHOS, and THE AXES ARE ON THE SAME SCALE AND SQUARE, so "spreads
# vertically and not horizontally" is a fact about the picture and not about
# how it was drawn. If the claim holds, the right-hand column is a vertical
# band: the machinery is organised along OXPHOS and not along MYC.
ADJ_LEVELS <- c("raw", "adj. PROLIF_DISJOINT")
g1dat <- wide %>%
  # names_sep = "_" would split "adj. PROLIF_DISJOINT" as well and produce an
  # NA facet. The pattern anchors the split to the axis prefix only.
  tidyr::pivot_longer(dplyr::starts_with(c("MYC_", "OXPHOS_")),
                      names_to = c("axis", "adjustment"),
                      names_pattern = "^(MYC|OXPHOS)_(.*)$",
                      values_to = "rho") %>%
  tidyr::pivot_wider(names_from = axis, values_from = rho) %>%
  dplyr::mutate(adjustment = factor(adjustment, levels = ADJ_LEVELS),
                cohort = factor(cohort, levels = names(COHORT_COLS)))
LIM1 <- max(abs(c(g1dat$MYC, g1dat$OXPHOS))) * 1.05
sd_lab <- g1dat %>% dplyr::group_by(cohort, adjustment) %>%
  dplyr::summarise(lab = sprintf("SD along MYC %.3f\nSD along OXPHOS %.3f",
                                 stats::sd(MYC), stats::sd(OXPHOS)),
                   .groups = "drop")
g1 <- ggplot2::ggplot(g1dat, ggplot2::aes(MYC, OXPHOS)) +
  ggplot2::annotate("rect", xmin = -0.1, xmax = 0.1, ymin = -LIM1, ymax = LIM1,
                    fill = "grey80", alpha = 0.35) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(ggplot2::aes(colour = mitocarta,
                                   shape = in_prolif_cov), size = 1.9) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(g1dat, abs(OXPHOS) > 0.3 | abs(MYC) > 0.3),
    ggplot2::aes(label = gene), size = 2.2, max.overlaps = 18,
    show.legend = FALSE, seed = PROJECT_SEED) +
  ggplot2::geom_text(data = sd_lab, ggplot2::aes(label = lab), x = -LIM1 * 0.95,
                     y = LIM1 * 0.92, hjust = 0, vjust = 1, size = 2.6,
                     colour = "grey25", inherit.aes = FALSE) +
  ggplot2::facet_grid(cohort ~ adjustment) +
  ggplot2::coord_fixed(xlim = c(-LIM1, LIM1), ylim = c(-LIM1, LIM1)) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey55",
                                          `TRUE` = "#d7191c"),
                               name = "in MitoCarta 3.0") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                              name = "in the proliferation covariate") +
  ggplot2::labs(
    title = "The apoptotic machinery spreads along OXPHOS, and much less along MYC",
    subtitle = paste("EXPLORATORY - not pre-registered | 44 canonical genes |",
                     "partial Spearman on PROLIF_DISJOINT (318 genes)"),
    x = "per-gene correlation with MYC activity (FELSHER__MITOSTRIP)",
    y = "per-gene correlation with OXPHOS subunits",
    caption = paste0(
      "The axes are square and identically scaled, so the shape of the cloud is\n",
      "the result: it is taller than it is wide in all four panels. Partialling\n",
      "out proliferation NARROWS it horizontally and leaves the vertical spread\n",
      "alone - the SD figures in each panel are the claim. It does not flatten\n",
      "MYC to nothing; the shaded band is +/- 0.1 and several genes sit outside\n",
      "it. Crosses are TP53 and BIRC5, which are IN the covariate and so are\n",
      "partly adjusted for themselves.\n",
      "TWO FIGURES MUST BE READ WITH THIS ONE. Figure 2 shows the adjustment does\n",
      "NOT flatten MYC against the mitoribosome, so this is not simply an emptied\n",
      "score. Figure 5 shows the vertical spread is what ANY 20 mitochondrial\n",
      "plus 24 other genes of this expression profile give, so it is not by\n",
      "itself evidence of an apoptotic programme.")) +
  theme_e11
.save(g1, "E11_fig1_machinery_plane_before_after", 8.5, 8.5)

# --- FIG 2: THE CONTROL that makes figure 1 mean anything --------------------
g2dat <- null_tests %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"),
                adjustment %in% ADJ_LEVELS) %>%
  dplyr::mutate(adjustment = factor(adjustment, levels = ADJ_LEVELS),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                set = factor(set, levels = rev(names(SETS))),
                self = axis == "OXPHOS" & set == "OXPHOS subunits arm (89)")
g2 <- ggplot2::ggplot(g2dat, ggplot2::aes(z_mean_abs, set,
                                          fill = adjustment)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                    width = 0.68) +
  ggplot2::geom_point(data = dplyr::filter(g2dat, self),
                      ggplot2::aes(x = z_mean_abs), shape = 4, size = 2.6,
                      colour = "black", stroke = 0.8,
                      position = ggplot2::position_dodge(width = 0.75),
                      inherit.aes = TRUE, show.legend = FALSE) +
  ggplot2::facet_grid(cohort ~ axis, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c(raw = "grey65",
                                        `adj. PROLIF_DISJOINT` = "#1b9e77"),
                             name = NULL) +
  ggplot2::labs(
    title = "The control, and the correction it forces",
    subtitle = paste("EXPLORATORY - not pre-registered | mean |rho| above an",
                     "expression-matched null, as z |", NULL_DRAWS, "draws"),
    x = "z of mean |rho| against an expression-matched null of the same size",
    y = NULL,
    caption = paste0(
      "TWO THINGS TO READ. First, the mitoribosome bar stays tall on MYC after\n",
      "adjustment - the adjustment did not empty the estimator, so a collapse\n",
      "elsewhere is about the set and not about the score. Second, and it is a\n",
      "correction to the obvious reading of figure 1: the machinery's bar is\n",
      "SHORT ON BOTH AXES. Its mean |rho| with OXPHOS is no larger than a random\n",
      "expression-matched set of 44 genes achieves, because an OXPHOS score\n",
      "correlates with much of the transcriptome. What differs between the axes\n",
      "is the STRUCTURE of those correlations, not their size - figures 1 and 4.\n",
      "Crossed bars are the OXPHOS arm against the OXPHOS score, which is that\n",
      "score's own genes. z is a distance from a null, not a p-value.")) +
  theme_e11
.save(g2, "E11_fig2_control_vs_matched_null", 9, 5.5)

# --- FIG 3: what the adjustment does to each gene, one line each -------------
g3dat <- gene_tab %>%
  dplyr::filter(axis %in% c("MYC", "OXPHOS"), adjustment %in% ADJ_LEVELS) %>%
  dplyr::mutate(adjustment = factor(adjustment, levels = ADJ_LEVELS),
                cohort = factor(cohort, levels = names(COHORT_COLS)))
g3 <- ggplot2::ggplot(g3dat, ggplot2::aes(adjustment, rho, group = gene)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_line(ggplot2::aes(colour = mitocarta), linewidth = 0.5,
                     alpha = 0.75) +
  ggplot2::geom_point(ggplot2::aes(colour = mitocarta), size = 1.1) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey55",
                                          `TRUE` = "#d7191c"),
                               name = "in MitoCarta 3.0") +
  ggplot2::labs(
    title = "What proliferation adjustment does to each of the 44 genes",
    subtitle = paste("EXPLORATORY - not pre-registered | one line per gene;",
                     "partial Spearman on PROLIF_DISJOINT"),
    x = NULL, y = "per-gene Spearman rho with the axis",
    caption = paste0(
      "READ THE CROSSINGS, NOT ONLY THE HEIGHT. On OXPHOS the lines are almost\n",
      "horizontal: proliferation explains nearly none of the ordering, and the\n",
      "gene at the top before is the gene at the top after. On MYC the lines\n",
      "cross heavily - the adjustment does not so much shrink that column as\n",
      "REARRANGE it, which is what a correlation carried by a shared third\n",
      "factor looks like. Read beside figure 2: an axis that had been emptied\n",
      "would collapse rather than reshuffle, and the control shows it was not.")) +
  theme_e11
.save(g3, "E11_fig3_adjustment_slopegraph", 7.5, 6)

# --- FIG 4: does the localisation split survive? -----------------------------
g4dat <- dplyr::bind_rows(
  s6_adj %>% dplyr::select(cohort, axis, adjustment, vs_mitocarta,
                           vs_pro_death) %>%
    tidyr::pivot_longer(c(vs_mitocarta, vs_pro_death), names_to = "split",
                        values_to = "value")) %>%
  dplyr::filter(adjustment %in% ADJ_LEVELS) %>%
  dplyr::mutate(adjustment = factor(adjustment, levels = ADJ_LEVELS),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                split = dplyr::recode(split,
                  vs_mitocarta = "acts at the mitochondrion (MitoCarta 3.0)",
                  vs_pro_death = "annotated pro-death"))
g4 <- ggplot2::ggplot(g4dat, ggplot2::aes(value, split, colour = adjustment)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55),
                      size = 2.6) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::scale_colour_manual(values = c(raw = "grey45",
                                          `adj. PROLIF_DISJOINT` = "#1b9e77"),
                               name = NULL) +
  ggplot2::labs(
    title = "Does the localisation split survive proliferation adjustment?",
    subtitle = paste("EXPLORATORY - not pre-registered | Spearman of the 44",
                     "per-gene rho against each binary split; TP53 and BIRC5",
                     "excluded"),
    x = "Spearman of the per-gene rho with the split", y = NULL,
    caption = paste0(
      "E10 found this split at 0.453 on OXPHOS and 0.140 on MYC, unadjusted,\n",
      "with MitoCarta supplying the localisation independently of the death\n",
      "curation. A split that holds under adjustment is not proliferation.\n",
      "TP53 and BIRC5 are dropped here because they are in the covariate.")) +
  theme_e11
.save(g4, "E11_fig4_localisation_split_under_adjustment", 8, 4.5)

# --- FIG 5: is the localisation split anything more than composition? -------
g5dat <- split_null %>%
  dplyr::mutate(adjustment = factor(adjustment, levels = ADJ_LEVELS),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                lo = null_mean - 1.96 * null_sd,
                hi = null_mean + 1.96 * null_sd)
g5dat <- g5dat %>% dplyr::mutate(
  pool = factor(pool, levels = names(POOLS)))
g5 <- ggplot2::ggplot(g5dat, ggplot2::aes(y = adjustment)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi), linewidth = 3,
                          colour = "grey80") +
  ggplot2::geom_point(ggplot2::aes(x = null_mean), shape = 124, size = 4,
                      colour = "grey40") +
  ggplot2::geom_point(ggplot2::aes(x = observed_split), size = 3,
                      colour = "#d7191c") +
  ggplot2::facet_grid(cohort + pool ~ axis) +
  ggplot2::labs(
    title = "Is the localisation split about apoptosis, or about composition?",
    subtitle = paste("EXPLORATORY - not pre-registered | red = observed over",
                     "the 44 | grey = null mean +/- 2 SD, 20 MitoCarta + 24",
                     "other genes"),
    x = "Spearman of the per-gene rho with MitoCarta membership", y = NULL,
    caption = paste0(
      "The grey bar is the null's mean plus or minus two SD: how large a split\n",
      "arises from ANY 20 mitochondrial and 24 non-mitochondrial genes of the\n",
      "machinery's expression profile. A red point inside the grey bar means the\n",
      "split is generic mitochondrial biology and E10 R1 must be restated; a red\n",
      "point outside it means the machinery is steeper than composition alone.\n",
      "The strict pool removes every OXPHOS and mitoribosome gene from the null,\n",
      "because those correlate with an OXPHOS score by definition while the\n",
      "machinery's 20 mitochondrial genes are BCL2-family members and caspases.\n",
      "This null was not run in E10, which is why R1 needs this figure beside it.")) +
  theme_e11
.save(g5, "E11_fig5_split_vs_composition_null", 8, 6.5)

# --- FIG 6: the two nulls side by side --------------------------------------
g6dat <- dplyr::bind_rows(
  split_null %>% dplyr::select(cohort, axis, adjustment, pool, observed_split,
                               null_mean, null_sd, z),
  split_null_submito %>% dplyr::select(cohort, axis, adjustment, pool,
                                       observed_split, null_mean, null_sd, z)) %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::mutate(pool = factor(pool, levels = c(names(POOLS),
                                               "sub-compartment matched")),
                cohort = factor(cohort, levels = names(COHORT_COLS)),
                lo = null_mean - 1.96 * null_sd,
                hi = null_mean + 1.96 * null_sd)
g6 <- ggplot2::ggplot(g6dat, ggplot2::aes(y = pool)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi), linewidth = 3.5,
                          colour = "grey80") +
  ggplot2::geom_point(ggplot2::aes(x = null_mean), shape = 124, size = 4,
                      colour = "grey40") +
  ggplot2::geom_point(ggplot2::aes(x = observed_split), size = 3.2,
                      colour = "#d7191c") +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::labs(
    title = "Three nulls for one number, each stricter than the last",
    subtitle = paste("EXPLORATORY - not pre-registered | proliferation-adjusted",
                     "| red = observed over the 44 | grey = null mean +/- 2 SD"),
    x = "Spearman of the per-gene rho with MitoCarta membership", y = NULL,
    caption = paste0(
      "Top row matches only MitoCarta membership. The middle row removes every\n",
      "OXPHOS and mitoribosome gene from the pool, because those track an OXPHOS\n",
      "score by definition. The bottom row also matches the SUB-COMPARTMENT -\n",
      "13 outer membrane, 5 intermembrane space, 2 inner membrane, as the\n",
      "machinery is - because E10 found a ladder by depth into the organelle.\n",
      "E10 R1 read 0.453 as a fact about apoptosis. It is a statement about\n",
      "apoptosis only to the extent the red point escapes the grey bar.")) +
  theme_e11
.save(g6, "E11_fig6_three_nulls", 8.5, 5)

# --- FIG 7: does the contrast hold inside a subtype? -------------------------
g7dat <- compartment %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS))) %>%
  tidyr::pivot_longer(c(sd_rho, observed_split), names_to = "statistic",
                      values_to = "value") %>%
  dplyr::mutate(statistic = dplyr::recode(statistic,
    sd_rho = "spread of the 44 (SD of rho)",
    observed_split = "localisation split"))
g7null <- compartment %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT", !is.na(null_mean)) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                statistic = "localisation split",
                lo = null_mean - 1.96 * null_sd, hi = null_mean + 1.96 * null_sd)
g7 <- ggplot2::ggplot(g7dat, ggplot2::aes(value, stratum, colour = axis)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_linerange(data = g7null,
                          ggplot2::aes(xmin = lo, xmax = hi, y = stratum),
                          inherit.aes = FALSE, linewidth = 3,
                          colour = "grey85") +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_grid(cohort ~ statistic, scales = "free_x") +
  ggplot2::scale_colour_manual(values = c(MYC = "#7570b3", OXPHOS = "#1b9e77"),
                               name = NULL) +
  ggplot2::labs(
    title = "Does the contrast survive inside a single subtype?",
    subtitle = paste("EXPLORATORY - not pre-registered | proliferation-adjusted",
                     "| Luminal = LumA + LumB (696 / 2,436), Basal (171 / 317)"),
    x = NULL, y = NULL,
    caption = paste0(
      "The `all` row is the pooled value and is NOT an average of the two below\n",
      "it. A pooled point outside the range of its own two compartments is\n",
      "reading a difference BETWEEN subtypes, which is what E10 fig9 found for\n",
      "the BCL2-family ratios on the MYC axis - though here it runs the other\n",
      "way for MYC, where pooling UNDERSTATES both compartments. Grey bars on\n",
      "the localisation-split panel are the sub-compartment-matched null for\n",
      "that stratum: a point inside its own bar is composition, and all of them\n",
      "are. Basal is 171 TCGA samples and everything about it is correspondingly\n",
      "wide.")) +
  theme_e11
.save(g7, "E11_fig7_contrast_by_compartment", 9, 5)

# --- FIG 8: the conditioning ladder, which is the paper's test ---------------
g8dat <- conditioned %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                model = factor(model, levels = rev(names(COND_MODELS))),
                lo = null_mean - 1.96 * null_sd,
                hi = null_mean + 1.96 * null_sd,
                axis = factor(axis, levels = c("OXPHOS", "MYC")))
g8 <- ggplot2::ggplot(g8dat, ggplot2::aes(y = model)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi), linewidth = 3.5,
                          colour = "grey85") +
  ggplot2::geom_point(ggplot2::aes(x = null_mean), shape = 124, size = 4,
                      colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(x = split, colour = axis), size = 3.4) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_colour_manual(values = c(OXPHOS = "#1b9e77", MYC = "#7570b3"),
                               name = NULL) +
  ggplot2::labs(
    title = "Which axis actually orders the apoptotic machinery?",
    subtitle = paste("EXPLORATORY - not pre-registered | 44 genes | coloured",
                     "point = observed localisation split | grey = null mean",
                     "+/- 2 SD"),
    x = "Spearman of the per-gene correlation with MitoCarta membership",
    y = NULL,
    caption = paste0(
      "Every model also conditions on proliferation. The two axes are\n",
      "themselves correlated and rank these 44 genes almost identically, so\n",
      "comparing them SEPARATELY cannot say which one carries the ordering.\n",
      "Conditioning OXPHOS on MYC leaves the split intact - it rises. \n",
      "Conditioning MYC on OXPHOS takes it to zero. MYC's apparent ordering\n",
      "was inherited from OXPHOS; the reverse is not true.\n",
      "The grey band is the sub-compartment-matched null: even the surviving\n",
      "OXPHOS split is only about one SD above what any 13 outer-membrane, 5\n",
      "intermembrane, 2 inner-membrane and 24 non-mitochondrial genes give.\n",
      "The ordering is real and it is OXPHOS's, but it is not specific to\n",
      "apoptosis.")) +
  theme_e11
.save(g8, "E11_fig8_conditioning_ladder", 9, 4.5)

# =============================================================================
# 6. Save
# =============================================================================
message("\n6. save")
saveRDS(list(
  null_tests = null_tests, split_null = split_null,
  split_null_submito = split_null_submito, compartment = compartment,
  gene_compartment = gene_compartment, pooled_outside = outside,
  conditioned = conditioned,
  gene_tab = gene_tab, wide = wide,
  spread = spread, s6_adj = s6_adj, replication = replication,
  purity_tab = purity_tab, overlap_audit = overlap_audit,
  sets = SETS, self_adjusted = SELF_ADJ,
  settings = list(null_draws = NULL_DRAWS, n_bins = N_BINS,
                  prolif_covariate = PROLIF_REF_COV, myc_axis = MYC_REF,
                  adjustments = names(ADJUSTMENTS), seed = PROJECT_SEED,
                  strata = c("all", STRATA_E11),
                  gene_scale = "linear DESeq2-normalised, rank-transformed"),
  rules = list(
    claim = paste("the question was: after correcting for proliferation, does",
                  "OXPHOS still correlate with the apoptotic machinery while",
                  "MYC does not. Formed AFTER seeing E10's unadjusted numbers,",
                  "so it is exploratory and carries a control rather than a",
                  "p-value."),
    control = paste("the mitoribosome arm is the control that decides whether",
                    "a MYC collapse is about apoptosis or about the estimator.",
                    "83 genes, sharing none with the 44, none with the OXPHOS",
                    "arm and none with the proliferation covariate. Figure 1",
                    "must never be shown without figure 2."),
    covariate = paste("PROLIF_DISJOINT is the primary covariate because it is",
                      "disjoint from FELSHER__MITOSTRIP, from the 89 OXPHOS-arm",
                      "genes and from the 83 mitoribosome genes, so the same",
                      "adjustment can be applied to both axes. PROLIF_STD",
                      "shares 9 genes with the MYC reference and is secondary."),
    two_corrections = paste("E1: partialling on a proliferation score and",
                            "removing proliferation genes from the estimator",
                            "disagree in sign. Both are computed here and",
                            "neither is reportable without the other."),
    self_adjustment = paste("TP53 and BIRC5 are IN the proliferation covariate",
                            "and are partly adjusted for themselves. Every",
                            "summary is given with and without them and they",
                            "are marked on figure 1."),
    submito = paste("section 3.2 is stricter again: it matches the",
                    "sub-mitochondrial compartment of the machinery's 20 genes",
                    "(13 MOM, 5 IMS, 2 MIM) as well as membership, because E10",
                    "found a ladder by depth into the organelle. Its pools are",
                    "small - IMS has 39 usable genes - so its expression",
                    "matching is coarser than section 3's, which is a limit of",
                    "the annotation and not a choice."),
    conditioning = paste("section 4.3 is the test the claim needs. The two",
                         "axes correlate at 0.26-0.32 after proliferation and",
                         "rank these 44 genes at 0.61-0.73, so comparing them",
                         "separately cannot say which carries the ordering.",
                         "Conditioning OXPHOS on MYC leaves the localisation",
                         "split intact; conditioning MYC on OXPHOS takes it to",
                         "zero. That, not the separate comparison, is what",
                         "licenses 'OXPHOS rather than MYC'."),
    compartments = paste("section 4.2 asks whether P1 holds inside a subtype.",
                         "Basal is 171 TCGA and 317 SCAN-B samples and a",
                         "171-sample stratum gives a 95% interval about +/-",
                         "0.15 wide on one gene's rho, so a Luminal-minus-Basal",
                         "difference under roughly 0.2 is not separable from",
                         "sampling there. A POOLED value outside the range of",
                         "both compartments is a between-subtype effect - the",
                         "D3/S1 artefact, which E10 fig9 found for the",
                         "BCL2-family ratios on the MYC axis."),
    composition = paste("section 3.1 is the null E10 R1 did not have: 20",
                        "MitoCarta plus 24 other genes, expression-matched",
                        "within each half. If the observed split sits inside",
                        "that null, the localisation finding is generic",
                        "mitochondrial biology rather than a fact about",
                        "apoptosis."),
    self_correlation = paste("CYCS is in the OXPHOS subunits arm, and the whole",
                             "OXPHOS-arm set read against the OXPHOS score is",
                             "that score's own genes. Both are marked, neither",
                             "is a result."),
    purity = paste("the purity and leukocyte adjustment is TCGA-only, n =",
                   sum(PURITY_OK), ". SCAN-B has no purity estimate and it is",
                   "never imputed, so that row cannot replicate.")),
  built = Sys.time()), PATH_E11)
readr::write_csv(gene_tab, PATH_E11_CSV)

message("\nE11: done.")
message("    results/prolif_adjusted_machinery.rds")
message("    outputs/tables/E11_gene_rho_by_adjustment.csv")
message("    8 figures in outputs/figures/:")
message("      fig1 THE PICTURE - the 44 genes on the MYC-OXPHOS plane,")
message("           before and after proliferation is partialled out")
message("      fig2 THE CONTROL - the same adjustment against the mitoribosome,")
message("           which is what says fig1 is about apoptosis and not the score")
message("      fig3 the per-gene slopegraph of what the adjustment does")
message("      fig4 whether the localisation split survives the adjustment")
message("      fig5 whether that split is more than mitochondrial composition")
message("      fig6 the same number against three nulls, each stricter")
message("      fig7 whether the contrast survives inside a single subtype")
message("      fig8 WHICH AXIS ORDERS THE MACHINERY - each conditioned on the")
message("           other, which is the test the comparison actually needs")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E11)

  # THE CLAIM, as one table. `z_mean_abs` is the set's mean |rho| against an
  # expression-matched null; `z_sd` is its SPREAD against the same null, and the
  # two say different things for a set that is half positive and half negative.
  x$null_tests %>%
    dplyr::filter(axis %in% c("MYC", "OXPHOS"),
                  adjustment %in% c("raw", "adj. PROLIF_DISJOINT")) %>%
    dplyr::select(cohort, axis, adjustment, set, mean_abs_rho, null_mean_abs,
                  z_mean_abs, sd_rho, z_sd) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # THE CONTROL on its own - does adjusted MYC still track anything?
  x$null_tests %>%
    dplyr::filter(set == "mitoribosome arm (83)") %>%
    dplyr::select(cohort, axis, adjustment, mean_abs_rho, z_mean_abs) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # do the two corrections agree? (partialling vs removing the genes)
  x$null_tests %>%
    dplyr::filter(set == "apoptotic machinery (44)", grepl("^MYC", axis)) %>%
    dplyr::select(cohort, axis, adjustment, mean_abs_rho, z_mean_abs) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # the three nulls for the localisation split, stacked
  dplyr::bind_rows(
    dplyr::select(x$split_null, cohort, axis, adjustment, pool,
                  observed_split, null_mean, null_sd, z),
    dplyr::select(x$split_null_submito, cohort, axis, adjustment, pool,
                  observed_split, null_mean, null_sd, z)) %>%
    dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::arrange(axis, cohort, pool) %>% as.data.frame()

  # spread, replication, and whether S6 survives adjustment
  x$spread %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$replication %>% as.data.frame()
  x$s6_adj %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # inside a subtype - and is the pooled value outside both compartments?
  x$compartment %>% dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
    dplyr::select(cohort, stratum, n_samples, axis, sd_rho, observed_split,
                  null_mean, z) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()
  x$pooled_outside %>% as.data.frame()

  # which genes moved when proliferation came out
  x$gene_tab %>%
    dplyr::filter(axis == "MYC",
                  adjustment %in% c("raw", "adj. PROLIF_DISJOINT")) %>%
    tidyr::pivot_wider(id_cols = c(cohort, gene), names_from = adjustment,
                       values_from = rho) %>%
    dplyr::mutate(shift = `adj. PROLIF_DISJOINT` - raw) %>%
    dplyr::arrange(shift) %>% as.data.frame()

}
