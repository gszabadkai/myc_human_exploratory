# E06_estimator_anatomy.R
# =============================================================================
# PHASE 2. Why do the MYC estimators disagree the way they do?
#
# Answers two questions from the 2026-09-01 handoff section 5:
#   - MENSSEN_MYC_TARGETS is the TOP signature on OXPHOS (0.670 / 0.616) at only
#     13.2% proliferation entanglement and 53 genes, above every Hallmark set.
#     ELLWOOD (13 genes) and ALFANO are what make F1's entanglement slope
#     positive. Are they informative or merely weak?
#   - M_b, the CollecTRI regulon estimator, is the WEAKEST on the OXPHOS plane
#     (+0.240 / +0.117, and it sign-flips under proliferation adjustment in
#     SCAN-B) yet gives by far the SHARPEST BCL2-family signal (BCL2A1 +0.658,
#     BID +0.556, BCL2 -0.550). Nothing in phase 1 explains that.
#
# =============================================================================
# THE THING FOUND WHILE WRITING THIS SCRIPT, WHICH BEARS DIRECTLY ON F1
# =============================================================================
# D0 found that the death sets' correlation with OXPHOS tracks how MITOCHONDRIAL
# the gene set is, at Spearman 0.60-0.81. Nobody had asked whether the same is
# true of the MYC signatures. It is - and worse, THE PANEL IS NOT HOMOGENEOUS.
#
#   `FELSHER_61` and `M_b` were STRIPPED of every MitoCarta gene by the
#   validation study, against the full 1,136-gene inventory, precisely so the
#   MYC estimator could not overlap its own exposure. FELSHER lost 6 of 67 and
#   the CollecTRI regulon 75 of 886. Both now carry ZERO mitochondrial genes.
#
#   THE OTHER SIXTEEN SIGNATURES WERE NEVER STRIPPED. HALLMARK_MYC_TARGETS_V1
#   carries 23 MitoCarta genes, BILD 20, DANG_MYC_TARGETS_UP 19, MYC_UP.V1_UP
#   18, MENSSEN 12 of only 53.
#
# So phase 1's estimator panel mixes two confound-free estimators with sixteen
# that share genes with the arm they are correlated against, and F1 read the
# spread across all eighteen as if they were the same kind of object. That much
# is simply a fact about the inputs.
#
# WHETHER IT MATTERS IS A SEPARATE QUESTION, AND THIS SCRIPT MUST NOT PREJUDGE
# IT. D0's analogy predicts contamination: a signature correlates with OXPHOS
# because it contains OXPHOS genes. But there is a second explanation that fits
# the same association, and it is not a confound at all - a signature that
# happens to contain many mitochondrial genes may simply BE a more
# mitochondrially-inclined description of MYC, its other genes correlating with
# OXPHOS just as much.
#
# Section 3.1 separates them by deleting each signature's mitochondrial genes
# and recomputing:
#
#   if the panel's spread COLLAPSES        -> contamination; F1 must be re-read
#   if the spread and the frac_mito
#     association BOTH SURVIVE             -> the mitochondrial genes are a
#                                             MARKER of an inclination, not the
#                                             CAUSE of the correlation, and F1
#                                             stands
#
# Both outcomes are informative and only one of them is bad news. Read the
# section-3.1 output before believing either.
#
# SCALE: gene-set membership is scale-free; the score correlations are read from
# E03's atlas; the two new mean-z scores are built on the LOG matrix and stated
# where they are built. SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE06: the anatomy of the MYC estimator panel\n", strrep("=", 78))

PATH_E06 <- file.path(DIR_RESULTS, "estimator_anatomy.rds")

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

a    <- readRDS(file.path(DIR_RESULTS, "correlation_atlas.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
mito <- readRDS(PATH_TCGA_MITO)
g1   <- readRDS(PATH_G1)

ID_T <- colnames(mito$gsva_arms); ID_S <- colnames(sc$gsva_arms)
MYC_SETS  <- sd_$myc_sets
MYC_PANEL <- sd_$myc_panel

# The reference sets every membership question is asked against.
MITO_ALL  <- unique(unlist(mito$mito_paths, use.names = FALSE))
ARM_REF   <- sd_$arm_sets[c("OXPHOS subunits", "OXPHOS umbrella",
                            "Mitochondrial ribosome", "TCA cycle",
                            "mtDNA-encoded OXPHOS")]
message("   ", length(MYC_SETS), " MYC signatures | MitoCarta pathway genes: ",
        length(MITO_ALL))

# =============================================================================
# 2. What is actually in these signatures
# =============================================================================
message("\n2. composition of the panel")

# Which estimators the validation study stripped, and against what. This is the
# distinction phase 1 did not know it was making.
STRIPPED <- c("FELSHER_61", "M_b")
strip_summary <- g1$strip_summary
message("   the validation study's overlap stripping, at d3ac60e:")
strip_summary %>% as.data.frame() %>% print(row.names = FALSE)
message("   stripped against: ", g1$meta$strip_set)
message("   mitochondrial genes remaining after stripping: ",
        paste(vapply(names(g1$estimators_stripped), function(n)
          paste0(n, "=", sum(g1$estimators_stripped[[n]] %in% MITO_ALL)),
          character(1)), collapse = ", "))

composition <- tibble::tibble(signature = names(MYC_SETS)) %>%
  dplyr::mutate(
    stripped     = signature %in% STRIPPED,
    n            = vapply(MYC_SETS, length, integer(1))[signature],
    frac_prolif  = MYC_PANEL$frac_prolif[match(signature, MYC_PANEL$signature)],
    n_mito       = vapply(MYC_SETS[signature],
                          function(g) sum(g %in% MITO_ALL), integer(1)),
    frac_mito    = n_mito / n,
    n_oxphos     = vapply(MYC_SETS[signature],
                          function(g) sum(g %in% ARM_REF[["OXPHOS subunits"]]),
                          integer(1)),
    n_mitoribo   = vapply(MYC_SETS[signature],
                          function(g) sum(g %in% ARM_REF[["Mitochondrial ribosome"]]),
                          integer(1)),
    n_tca        = vapply(MYC_SETS[signature],
                          function(g) sum(g %in% ARM_REF[["TCA cycle"]]), integer(1)))

# The atlas value each signature is being explained
rho_ox <- a$atlas %>%
  dplyr::filter(arm == "OXPHOS subunits", instrument == "gsva",
                stratum == "all", adjusted == "raw",
                kind == "signature (GSVA)") %>%
  dplyr::select(signature = myc_estimator, cohort, rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
  dplyr::rename(rho_TCGA = TCGA, rho_SCANB = `SCAN-B`)
composition <- composition %>% dplyr::left_join(rho_ox, by = "signature")

composition %>% dplyr::arrange(dplyr::desc(rho_TCGA)) %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                frac_mito = round(100 * frac_mito, 1),
                dplyr::across(c(rho_TCGA, rho_SCANB), ~ round(.x, 3))) %>%
  dplyr::select(signature, stripped, n, frac_prolif, n_mito, frac_mito,
                n_mitoribo, n_oxphos, n_tca, rho_TCGA, rho_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   `stripped` = TRUE means the estimator has NO mitochondrial genes",
        " by\n   construction. Only two of the eighteen do.")

# =============================================================================
# 3. THE D0 TEST, APPLIED TO THE MYC PANEL
# =============================================================================
# Does a signature's OXPHOS correlation track how mitochondrial the signature
# is, the way the death sets' did? And does mitochondrial content explain more
# of it than proliferation entanglement does?
message("\n3. does the panel's OXPHOS correlation track MITOCHONDRIAL content?")

d0_test <- tibble::tibble(
  cohort   = c("TCGA", "SCAN-B"),
  vs_mito_fraction = c(
    stats::cor(composition$frac_mito, composition$rho_TCGA,  method = "spearman"),
    stats::cor(composition$frac_mito, composition$rho_SCANB, method = "spearman")),
  vs_prolif_entanglement = c(
    stats::cor(composition$frac_prolif, composition$rho_TCGA,  method = "spearman"),
    stats::cor(composition$frac_prolif, composition$rho_SCANB, method = "spearman")),
  vs_set_size = c(
    stats::cor(composition$n, composition$rho_TCGA,  method = "spearman"),
    stats::cor(composition$n, composition$rho_SCANB, method = "spearman")))
d0_test %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the two signatures F1 turns on, and the one it is questioned by:")
composition %>%
  dplyr::filter(signature %in% c("MYC_UP.V1_UP", "MENSSEN_MYC_TARGETS",
                                 "ELLWOOD_MYC_TARGETS_UP", "ALFANO_MYC_TARGETS",
                                 "FELSHER_61", "HALLMARK_MYC_TARGETS_V1")) %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                frac_mito = round(100 * frac_mito, 1),
                dplyr::across(c(rho_TCGA, rho_SCANB), ~ round(.x, 3))) %>%
  dplyr::select(signature, n, frac_prolif, frac_mito, n_mitoribo, n_oxphos,
                rho_TCGA, rho_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)

# The direct question: drop every mitochondrial gene from each signature and
# ask how many genes are left. A signature that loses its correlation when its
# mitochondrial genes go was never measuring MYC-and-OXPHOS separately.
# --- 3.1 the cheap version of the fix ----------------------------------------
# Properly, every signature would be re-scored by GSVA with its mitochondrial
# genes removed - one call, under the pins, 15-25 minutes. This is the gene-level
# approximation, which needs no re-scoring: the mean per-gene Spearman of a
# signature's genes against the OXPHOS arm, computed with and without its
# mitochondrial members.
#
# IT IS AN APPROXIMATION AND NOT THE SAME STATISTIC. A GSVA score is a
# relative-rank enrichment, not a mean of per-gene correlations, and E05 showed
# the two can diverge by a lot. What it CAN do is rank the signatures by how
# much they stand to lose, which is what decides whether the full re-score is
# worth running.
message("\n3.1 delete each signature's mitochondrial genes (gene-level proxy)")

tcga_linm  <- readRDS(PATH_TCGA_LINEAR)
LT <- tcga_linm$mat[, ID_T, drop = FALSE]; rm(tcga_linm); invisible(gc(verbose = FALSE))
scanb_linm <- readRDS(PATH_SCANB_LINEAR)
LS <- scanb_linm$mat[, ID_S, drop = FALSE]; rm(scanb_linm); invisible(gc(verbose = FALSE))

.inL_t <- function(g) intersect(unique(g), rownames(LT))
.inL_s <- function(g) {
  h <- sc$symbol_map[unique(g)]
  intersect(unname(ifelse(is.na(h), unique(g), h)), rownames(LS))
}
ox_t <- as.numeric(mito$gsva_arms["OXPHOS subunits", ID_T])
ox_s <- as.numeric(sc$gsva_arms["OXPHOS subunits", ID_S])
grho_t <- .per_gene_rho(LT, ox_t)
grho_s <- .per_gene_rho(LS, ox_s)

.mean_rho <- function(genes, grho, inf) {
  i <- inf(genes); i <- i[is.finite(grho[i])]
  if (length(i) < 5L) return(NA_real_)
  mean(grho[i])
}
strip_test <- tibble::tibble(signature = names(MYC_SETS)) %>%
  dplyr::mutate(
    n_mito     = composition$n_mito[match(signature, composition$signature)],
    with_TCGA  = vapply(MYC_SETS[signature], .mean_rho, numeric(1), grho_t, .inL_t),
    without_TCGA = vapply(lapply(MYC_SETS[signature], setdiff, MITO_ALL),
                          .mean_rho, numeric(1), grho_t, .inL_t),
    with_SCANB = vapply(MYC_SETS[signature], .mean_rho, numeric(1), grho_s, .inL_s),
    without_SCANB = vapply(lapply(MYC_SETS[signature], setdiff, MITO_ALL),
                           .mean_rho, numeric(1), grho_s, .inL_s)) %>%
  dplyr::mutate(loss_TCGA = with_TCGA - without_TCGA,
                loss_SCANB = with_SCANB - without_SCANB)
strip_test %>% dplyr::arrange(dplyr::desc(loss_TCGA)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   how much of the PANEL SPREAD survives the deletion:")
tibble::tibble(
  quantity = c("spread across signatures, with mitochondrial genes",
               "spread across signatures, without them",
               "rho(frac_mito, mean per-gene rho), with",
               "rho(frac_mito, mean per-gene rho), without"),
  TCGA = c(diff(range(strip_test$with_TCGA, na.rm = TRUE)),
           diff(range(strip_test$without_TCGA, na.rm = TRUE)),
           stats::cor(composition$frac_mito, strip_test$with_TCGA,
                      method = "spearman", use = "pairwise.complete.obs"),
           stats::cor(composition$frac_mito, strip_test$without_TCGA,
                      method = "spearman", use = "pairwise.complete.obs")),
  `SCAN-B` = c(diff(range(strip_test$with_SCANB, na.rm = TRUE)),
               diff(range(strip_test$without_SCANB, na.rm = TRUE)),
               stats::cor(composition$frac_mito, strip_test$with_SCANB,
                          method = "spearman", use = "pairwise.complete.obs"),
               stats::cor(composition$frac_mito, strip_test$without_SCANB,
                          method = "spearman", use = "pairwise.complete.obs"))) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message(sprintf(
  "   the proxy tracks the GSVA correlation it stands in for at rho %.3f (TCGA) / %.3f (SCAN-B)",
  stats::cor(strip_test$with_TCGA, composition$rho_TCGA, method = "spearman",
             use = "pairwise.complete.obs"),
  stats::cor(strip_test$with_SCANB, composition$rho_SCANB, method = "spearman",
             use = "pairwise.complete.obs")))
message("   If the second row is much smaller than the first, most of what F1",
        " read as\n   disagreement between MYC signatures was disagreement",
        " about how many\n   mitochondrial genes they happen to contain. If it",
        " is barely smaller, the\n   mitochondrial genes MARK an inclination",
        " rather than causing the correlation.")
message("\n   NOTE: a mean of per-gene correlations is not a GSVA score, and the",
        "\n   two can diverge. This ranks the signatures by what they stand to",
        " lose;\n   only a mitochondria-stripped GSVA rescore - one call, under",
        " the pins -\n   settles it. The proxy's agreement with the real",
        " statistic is printed above.")
rm(LT, LS, grho_t, grho_s); invisible(gc(verbose = FALSE))

message("\n   MENSSEN's mitochondrial genes, named:")
menssen_mito <- intersect(MYC_SETS[["MENSSEN_MYC_TARGETS"]], MITO_ALL)
message("   ", length(menssen_mito), " of ",
        length(MYC_SETS[["MENSSEN_MYC_TARGETS"]]), ": ",
        paste(menssen_mito, collapse = ", "))
message("   MYC_UP.V1_UP's: ", length(intersect(MYC_SETS[["MYC_UP.V1_UP"]], MITO_ALL)),
        " of ", length(MYC_SETS[["MYC_UP.V1_UP"]]), ": ",
        paste(intersect(MYC_SETS[["MYC_UP.V1_UP"]], MITO_ALL), collapse = ", "))

# =============================================================================
# 4. Are ELLWOOD and ALFANO informative, or just weak?
# =============================================================================
# Three ways a signature can be "weak": too few genes, genes that do not move
# together (low internal coherence), or genes that move together but disagree
# with every other MYC signature.
message("\n4. ELLWOOD and ALFANO: weak, or measuring something else?")

tcga_vst <- readRDS(PATH_TCGA_VST); ET <- tcga_vst$mat[, ID_T, drop = FALSE]
rm(tcga_vst); invisible(gc(verbose = FALSE))
scanb_vst <- readRDS(PATH_SCANB_VST); ES <- scanb_vst$mat[, ID_S, drop = FALSE]
rm(scanb_vst); invisible(gc(verbose = FALSE))

.in_t <- function(g) intersect(unique(g), rownames(ET))
.in_s <- function(g) {
  h <- sc$symbol_map[unique(g)]
  intersect(unname(ifelse(is.na(h), unique(g), h)), rownames(ES))
}

# Internal coherence: the median pairwise gene-gene Spearman WITHIN the
# signature. A signature whose genes do not co-vary is not measuring a
# programme, whatever it is named.
.coherence <- function(genes, E, inf) {
  g <- inf(genes)
  if (length(g) < 3L) return(NA_real_)
  C <- stats::cor(t(.rank_rows(E[g, , drop = FALSE])))
  stats::median(C[lower.tri(C)])
}
coherence <- tibble::tibble(signature = names(MYC_SETS)) %>%
  dplyr::mutate(
    coh_TCGA  = vapply(MYC_SETS[signature], .coherence, numeric(1), ET, .in_t),
    coh_SCANB = vapply(MYC_SETS[signature], .coherence, numeric(1), ES, .in_s))

# Agreement with the rest of the panel: correlation of this signature's SCORE
# with the other 17 signatures' scores, within cohort.
.agreement <- function(gsva_new, sigs) {
  S <- gsva_new[sigs, , drop = FALSE]
  C <- stats::cor(t(.rank_rows(S)))
  diag(C) <- NA
  apply(C, 1L, function(v) stats::median(v, na.rm = TRUE))
}
SIGS <- MYC_PANEL$signature
agree_t <- .agreement(nw$tcga_gsva_new, SIGS)
agree_s <- .agreement(sc$gsva_new, SIGS)

panel_quality <- composition %>%
  dplyr::left_join(coherence, by = "signature") %>%
  dplyr::mutate(agree_TCGA  = agree_t[signature],
                agree_SCANB = agree_s[signature],
                coverage_TCGA = sd_$coverage$frac[match(paste("TCGA", signature),
                    paste(sd_$coverage$cohort, sd_$coverage$set))])
panel_quality %>% dplyr::arrange(frac_prolif) %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                dplyr::across(c(coh_TCGA, coh_SCANB, agree_TCGA, agree_SCANB,
                                rho_TCGA, rho_SCANB, coverage_TCGA),
                              ~ round(.x, 3))) %>%
  dplyr::select(signature, n, frac_prolif, coverage_TCGA, coh_TCGA, coh_SCANB,
                agree_TCGA, agree_SCANB, rho_TCGA, rho_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   If ELLWOOD and ALFANO sit at the BOTTOM of coherence and",
        " agreement, they\n   are weak signatures rather than informative ones,",
        " and F1's positive\n   entanglement slope is an artefact of including",
        " them.")

# =============================================================================
# 5. What M_b is, and why it behaves differently
# =============================================================================
# TWO STRUCTURAL DIFFERENCES, AND THE FIRST REVERSES THE QUESTION.
#
# (i) M_b IS MITOCHONDRIA-FREE BY CONSTRUCTION. 75 of the regulon's 886 targets
#     were stripped as MitoCarta genes, leaving zero. Section 3 shows the
#     panel's OXPHOS correlation tracks mitochondrial content at Spearman
#     0.6-0.7. So M_b being the WEAKEST estimator on the OXPHOS plane
#     (+0.240 / +0.117) is not a defect to be explained away - it may be the
#     least contaminated number in the study, and the sixteen unstripped
#     signatures above it may simply be reading their own mitochondrial genes.
#     The handoff asked "why does M_b behave differently"; the answer may be
#     that M_b behaves correctly.
#
# (ii) IT IS SIGNED. CollecTRI carries a mode of regulation per edge and ULM
#     uses it, so a repressed target going DOWN pushes M_b UP. No GSVA signature
#     has any such notion. Section 5.2 tests whether that is what makes M_b
#     move differently, by splitting the regulon by sign and scoring each half
#     UNSIGNED.
#
# AND THE CIRCULARITY WORRY IS THE OTHER WAY ROUND FROM THE OBVIOUS ONE. None of
# the 15 BCL2-family genes is in the stripped regulon - BAX, BCL2, BCL2L1, BBC3
# and PMAIP1 were in it and were removed as mitochondrial; BAK1, BID, BCL2A1 and
# MCL1 were never in it at all. So M_b's sharp BCL2-family signal (BCL2A1 +0.658,
# BID +0.556, BCL2 -0.550) is NOT self-referential. Section 5.3 verifies that
# rather than assuming it.
message("\n5. M_b: mitochondria-free, and signed")

# THE SIGN RULE IS NOT INVENTED HERE. It is copied from the script that built
# M_b - myc_human_validation script 06 at d3ac60e, and E02 section 4.5 for
# SCAN-B - so this split describes the estimator that exists rather than a
# different one:
#
#     mor = if_else(as.logical(is_stimulation), 1, -1)
#
# The snapshot is raw OmniPath, so there is no `mor` column: the sign lives in
# `is_stimulation` / `is_inhibition` and symbols in `*_genesymbol`. 88 of the
# 811 edges are flagged BOTH stimulatory and inhibitory. The rule above sends
# every one of them to +1, silently, because `is_stimulation` is TRUE. That is
# recorded upstream as "both-flagged edges take mor = +1" and it is reproduced
# here rather than improved on - but it is also reported, because a repressed-
# target analysis in which 88 ambiguous edges are counted as activating is a
# thing a reader should be told.
collectri <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE,
                             progress = FALSE)
COLL_ALL <- g1$estimators_stripped$COLLECTRI_MYC_ALL
myc_reg <- collectri %>%
  dplyr::filter(source_genesymbol == "MYC", !is.na(target_genesymbol),
                target_genesymbol != "") %>%
  dplyr::transmute(target = target_genesymbol,
                   stim = as.logical(is_stimulation),
                   inhib = as.logical(is_inhibition),
                   mor = dplyr::if_else(as.logical(is_stimulation), 1, -1)) %>%
  dplyr::distinct(target, .keep_all = TRUE) %>%
  dplyr::filter(target %in% COLL_ALL)
stopifnot(nrow(myc_reg) > 0)
message("   CollecTRI MYC regulon: ", nrow(myc_reg), " of ", length(COLL_ALL),
        " estimator targets carry an edge")
myc_reg %>% dplyr::count(stim, inhib, mor) %>% as.data.frame() %>%
  print(row.names = FALSE)

POS <- sort(unique(myc_reg$target[myc_reg$mor > 0]))
NEG <- sort(unique(myc_reg$target[myc_reg$mor < 0]))
N_AMBIG <- sum(myc_reg$stim & myc_reg$inhib)
message("   activated ", length(POS), " (of which ", N_AMBIG,
        " are BOTH-flagged and only counted as activating by the rule)",
        " | repressed ", length(NEG))
message("   mitochondrial fraction: activated ",
        round(mean(POS %in% MITO_ALL), 3), " | repressed ",
        round(mean(NEG %in% MITO_ALL), 3),
        " | panel median ", round(stats::median(composition$frac_mito), 3))

# --- 5.1 overlap with the GSVA panel ----------------------------------------
overlap <- tibble::tibble(signature = names(MYC_SETS)) %>%
  dplyr::mutate(
    n = vapply(MYC_SETS[signature], length, integer(1)),
    in_regulon = vapply(MYC_SETS[signature],
                        function(g) sum(g %in% COLL_ALL), integer(1)),
    frac_in_regulon = in_regulon / n)
message("\n   how much of each signature is inside the CollecTRI regulon:")
overlap %>% dplyr::arrange(dplyr::desc(frac_in_regulon)) %>%
  dplyr::mutate(frac_in_regulon = round(frac_in_regulon, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.2 the signed halves, scored unsigned ---------------------------------
# LOG matrix (VST). A mean of per-gene z across samples - the `zmean`
# construction from E02 section 4.4, reused so the halves are comparable to the
# arms.
.zmean <- function(genes, E, inf) {
  g <- inf(genes)
  sub <- E[g, , drop = FALSE]
  v <- apply(sub, 1L, stats::var)
  sub <- sub[v > 0, , drop = FALSE]
  colMeans((sub - rowMeans(sub)) / apply(sub, 1L, stats::sd))
}
halves_t <- rbind(regulon_activated = .zmean(POS, ET, .in_t),
                  regulon_repressed = .zmean(NEG, ET, .in_t))
halves_s <- rbind(regulon_activated = .zmean(POS, ES, .in_s),
                  regulon_repressed = .zmean(NEG, ES, .in_s))
halves_t <- rbind(halves_t, regulon_difference = halves_t[1, ] - halves_t[2, ])
halves_s <- rbind(halves_s, regulon_difference = halves_s[1, ] - halves_s[2, ])
colnames(halves_t) <- ID_T; colnames(halves_s) <- ID_S

# What each half tracks: M_b itself, the MYC signatures, and OXPHOS.
.probe <- function(halves, gsva_new, arms_obj, mb, ids, coh) {
  TARGETS <- rbind(
    M_b            = as.numeric(mb[ids]),
    FELSHER_61     = as.numeric(gsva_new["FELSHER_61", ids]),
    MYC_UP.V1_UP   = as.numeric(gsva_new["MYC_UP.V1_UP", ids]),
    OXPHOS_gsva    = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]),
    MITORIBO_gsva  = as.numeric(arms_obj$gsva_arms["Mitochondrial ribosome", ids]),
    PROLIF         = as.numeric(arms_obj$gsva_cov["PROLIF_DISJOINT", ids]))
  colnames(TARGETS) <- ids
  .atlas_block(TARGETS, halves[, ids, drop = FALSE], ids, NULL, min_n = 30L) %>%
    dplyr::rename(target = myc_estimator, half = measure) %>%
    dplyr::mutate(cohort = coh) %>%
    dplyr::select(cohort, half, target, rho, ci_lo, ci_hi)
}
est_t  <- readRDS(PATH_TCGA_MYC)$estimators
mb_t   <- stats::setNames(est_t$M_b, est_t$patient)
mb_probe <- dplyr::bind_rows(
  .probe(halves_t, nw$tcga_gsva_new, mito, mb_t,    ID_T, "TCGA"),
  .probe(halves_s, sc$gsva_new,      sc,   sc$M_b,  ID_S, "SCAN-B"))
mb_probe %>%
  tidyr::pivot_wider(id_cols = c(half, target), names_from = cohort,
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.3 the BCL2 family, which is where M_b was sharpest --------------------
message("\n   the BCL2 family: in the RAW regulon, and in the STRIPPED one?")
raw_reg <- collectri %>%
  dplyr::filter(source_genesymbol == "MYC") %>% dplyr::pull(target_genesymbol)
bcl2_in <- tibble::tibble(gene = sd_$bcl2_family) %>%
  dplyr::mutate(
    in_raw_regulon      = gene %in% raw_reg,
    in_stripped_regulon = gene %in% COLL_ALL,
    is_mitocarta        = gene %in% MITO_ALL,
    mor                 = myc_reg$mor[match(gene, myc_reg$target)],
    both_flagged        = myc_reg$stim[match(gene, myc_reg$target)] &
                          myc_reg$inhib[match(gene, myc_reg$target)])
bcl2_in %>% as.data.frame() %>% print(row.names = FALSE)
if (any(bcl2_in$in_stripped_regulon)) {
  stop("a BCL2-family gene survives into the estimator that scores it: ",
       paste(bcl2_in$gene[bcl2_in$in_stripped_regulon], collapse = ", "),
       ". M_b's correlation with that gene is partly self-referential and must ",
       "not be reported as an observation about it.", call. = FALSE)
}
message("   NONE of the 15 survives into the estimator, so M_b's BCL2-family",
        " numbers\n   are independent observations. The ones that WERE in the",
        " raw regulon\n   (", paste(bcl2_in$gene[bcl2_in$in_raw_regulon], collapse = ", "),
        ")\n   were removed as MitoCarta genes - the strip did the work here",
        " even though\n   it was aimed at something else.")

# The same question for every GSVA signature, which was NOT stripped.
message("\n   and the unstripped signatures - do any of them contain the BCL2",
        " genes\n   whose correlation with them E05 reports?")
bcl2_in_sigs <- tibble::tibble(signature = names(MYC_SETS)) %>%
  dplyr::mutate(n_bcl2 = vapply(MYC_SETS[signature],
                                function(g) sum(g %in% sd_$bcl2_family),
                                integer(1)),
                which = vapply(MYC_SETS[signature],
                               function(g) paste(intersect(g, sd_$bcl2_family),
                                                 collapse = ","), character(1))) %>%
  dplyr::filter(n_bcl2 > 0)
if (nrow(bcl2_in_sigs)) {
  bcl2_in_sigs %>% as.data.frame() %>% print(row.names = FALSE)
  message("   Those cells of E05's overlay are NOT independent observations.")
} else {
  message("   none.")
}

# =============================================================================
# 6. Save
# =============================================================================
message("\n6. save")
saveRDS(list(composition = composition, d0_test = d0_test,
             strip_summary = strip_summary, strip_test = strip_test,
             bcl2_in_signatures = bcl2_in_sigs,
             panel_quality = panel_quality, overlap = overlap,
             regulon = list(activated = POS, repressed = NEG,
                            n_both_flagged = N_AMBIG, table = myc_reg),
             mb_probe = mb_probe, bcl2_in_regulon = bcl2_in,
             halves = list(TCGA = halves_t, `SCAN-B` = halves_s),
             rules = list(
               d0_applies_here = paste("if a MYC signature's OXPHOS correlation",
                                       "tracks its mitochondrial fraction, that",
                                       "is the same confound D0 found in the",
                                       "death sets and F1 must be re-read"),
               circularity = paste("no BCL2-family gene survives into the",
                                   "stripped regulon, so M_b's BCL2 numbers are",
                                   "independent; the unstripped GSVA signatures",
                                   "are checked separately and some are not"),
               panel_not_homogeneous = paste("FELSHER_61 and M_b were stripped",
                                             "of every MitoCarta gene; the other",
                                             "16 signatures were not. Phase 1",
                                             "read the panel as if they were the",
                                             "same kind of object")),
             built = Sys.time()), PATH_E06)
message("\nE06: done.  results/estimator_anatomy.rds")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  e <- readRDS(PATH_E06)

  # THE question: mitochondrial content vs entanglement as explanations
  e$d0_test %>% as.data.frame()

  e$composition %>% dplyr::arrange(dplyr::desc(frac_mito)) %>%
    dplyr::select(signature, n, frac_prolif, frac_mito, rho_TCGA, rho_SCANB) %>%
    as.data.frame()

  # what each signature stands to lose to a proper mitochondria-stripped rescore
  e$strip_test %>% dplyr::arrange(dplyr::desc(loss_TCGA)) %>% as.data.frame()
  e$strip_summary %>% as.data.frame()

  # are the low-entanglement outliers weak?
  e$panel_quality %>% dplyr::arrange(coh_TCGA) %>%
    dplyr::select(signature, n, coh_TCGA, agree_TCGA, rho_TCGA) %>% head(8) %>%
    as.data.frame()

  # M_b
  e$mb_probe %>% tidyr::pivot_wider(id_cols = c(half, target),
                                    names_from = cohort, values_from = rho) %>%
    as.data.frame()
  e$bcl2_in_regulon %>% as.data.frame()

}
