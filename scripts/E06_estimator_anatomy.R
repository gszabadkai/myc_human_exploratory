# E06_estimator_anatomy.R
# =============================================================================
# PHASE 2. Why do the MYC estimators disagree the way they do - and how much of
# F1 survives once they are all labelled?
#
# =============================================================================
# WHAT WENT WRONG IN PHASE 1, AND WHAT E02 NOW DOES ABOUT IT
# =============================================================================
# F1 was written on a panel of 18 sets under bare names. One of them,
# `FELSHER_61`, had been stripped of every MitoCarta gene by the validation
# study; the other 17 had not. Nothing in any name said so. On top of that the
# "proliferation-adjusted" column used `PROLIF_DISJOINT`, which is PROLIF_STD
# minus the 9 proliferation genes of FELSHER_61 and nothing else - disjoint from
# M_a ALONE, and sharing 97 genes with the 18-signature union. So for 17 of 18
# signatures F1's proliferation adjustment was partly adjusting each signature
# for itself.
#
# E02 now scores FOUR EXPLICITLY LABELLED VARIANTS of every signature -
# __FULL, __MITOSTRIP, __PROLIFSTRIP, __BOTHSTRIP - so both questions can be
# asked by construction instead of by adjustment. This script reads them.
#
# THE FOUR ROWS AND WHAT EACH ONE DECIDES
# ---------------------------------------
#   __FULL         what F1 actually reported for 17 of the 18
#   __MITOSTRIP    if the correlation collapses here, it was the signature
#                  reading its own mitochondrial genes (D0's confound)
#   __PROLIFSTRIP  if it collapses here, it was proliferation - and this is the
#                  CLEAN version of that question, unlike the PROLIF_DISJOINT
#                  adjustment
#   __BOTHSTRIP    what is left when neither explanation is available. THIS ROW
#                  IS F1'S ANSWER.
#
# Sections 4 and 5 then take the two remaining handoff questions: whether
# ELLWOOD and ALFANO are informative or merely weak, and what M_b is.
#
# SCALE: correlations are read from E03's atlas. The two mean-z scores built in
# section 5 are on the LOG (VST) matrix and say so. SPECIES: human, natively.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))

message("\nE06: the anatomy of the MYC estimator panel\n", strrep("=", 78))

PATH_E06 <- file.path(DIR_RESULTS, "estimator_anatomy.rds")
FOCUS_ARM <- "OXPHOS subunits"

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
MITO_ALL  <- sd_$strip_refs$MITOCARTA_ALL
PROLIF_REF <- sd_$strip_refs$PROLIF_REF
stopifnot(length(MITO_ALL) == EXPECT_MITOCARTA_ALL,
          length(PROLIF_REF) == EXPECT_PROLIF_REF)
message("   ", dplyr::n_distinct(MYC_PANEL$base), " signatures x 4 variants = ",
        nrow(MYC_PANEL), " estimators")

message("\n   the validation study's overlap stripping, at d3ac60e:")
g1$strip_summary %>% as.data.frame() %>% print(row.names = FALSE)
message("   stripped against: ", g1$meta$strip_set)
message("   -> that is a MITOCHONDRIAL strip. No proliferation-stripped",
        " estimator\n      existed upstream; __PROLIFSTRIP is built here for",
        " the first time.")

# PROLIF_DISJOINT, named and measured, because F1 leaned on it.
pd <- sd_$cov_sets$PROLIF_DISJOINT
panel_union <- unique(unlist(MYC_SETS[MYC_PANEL$signature[
  MYC_PANEL$strip_status == "FULL"]], use.names = FALSE))
message("\n   PROLIF_DISJOINT: ", length(pd), " genes, ",
        sum(pd %in% PROLIF_REF), " of them HALLMARK E2F/G2M, and it shares ",
        sum(pd %in% panel_union), " genes\n   with the union of the FULL",
        " signatures. It is disjoint from M_a alone.")

# =============================================================================
# 2. The panel, with every label attached
# =============================================================================
message("\n2. the panel")

rho_ox <- a$atlas %>%
  dplyr::filter(arm == FOCUS_ARM, instrument == "gsva", stratum == "all",
                adjusted == "raw", kind == "signature (GSVA)") %>%
  dplyr::select(signature = myc_estimator, cohort, rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho) %>%
  dplyr::rename(rho_TCGA = TCGA, rho_SCANB = `SCAN-B`)
panel <- MYC_PANEL %>% dplyr::left_join(rho_ox, by = "signature")

panel %>% dplyr::filter(strip_status == "FULL") %>%
  dplyr::arrange(dplyr::desc(rho_TCGA)) %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                frac_mito = round(100 * frac_mito, 1),
                dplyr::across(c(rho_TCGA, rho_SCANB), ~ round(.x, 3))) %>%
  dplyr::select(base, n, frac_prolif, frac_mito, rho_TCGA, rho_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)

# D0's test, applied to the FULL variants. Does a signature's OXPHOS
# correlation track how mitochondrial it is, or how proliferative?
full <- panel %>% dplyr::filter(strip_status == "FULL")
d0_test <- tibble::tibble(
  cohort = c("TCGA", "SCAN-B"),
  vs_mito_fraction = c(
    stats::cor(full$frac_mito, full$rho_TCGA, method = "spearman"),
    stats::cor(full$frac_mito, full$rho_SCANB, method = "spearman")),
  vs_prolif_entanglement = c(
    stats::cor(full$frac_prolif, full$rho_TCGA, method = "spearman"),
    stats::cor(full$frac_prolif, full$rho_SCANB, method = "spearman")),
  vs_set_size = c(
    stats::cor(full$n, full$rho_TCGA, method = "spearman"),
    stats::cor(full$n, full$rho_SCANB, method = "spearman")))
message("\n   what the FULL panel's spread tracks:")
d0_test %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 3. THE TEST. Four variants, one table, and F1's answer in the last row.
# =============================================================================
message("\n3. what survives each strip")

# `panel` is E02's myc_panel joined to the atlas: 18 signatures x 4 variants and
# nothing else, so there is no `kind` to filter on. M_b lives in the atlas, not
# here, and section 5 reads it from there.
variant_table <- panel %>%
  dplyr::select(base, strip_status, n, rho_TCGA, rho_SCANB) %>%
  tidyr::pivot_wider(id_cols = base, names_from = strip_status,
                     values_from = c(n, rho_TCGA, rho_SCANB)) %>%
  dplyr::mutate(
    loss_mito_TCGA    = rho_TCGA_FULL - rho_TCGA_MITOSTRIP,
    loss_prolif_TCGA  = rho_TCGA_FULL - rho_TCGA_PROLIFSTRIP,
    loss_both_TCGA    = rho_TCGA_FULL - rho_TCGA_BOTHSTRIP,
    loss_mito_SCANB   = rho_SCANB_FULL - rho_SCANB_MITOSTRIP,
    loss_prolif_SCANB = rho_SCANB_FULL - rho_SCANB_PROLIFSTRIP,
    loss_both_SCANB   = rho_SCANB_FULL - rho_SCANB_BOTHSTRIP)

message("\n   rho with ", FOCUS_ARM, ", TCGA, by variant:")
variant_table %>% dplyr::arrange(dplyr::desc(rho_TCGA_FULL)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(base, rho_TCGA_FULL, rho_TCGA_MITOSTRIP, rho_TCGA_PROLIFSTRIP,
                rho_TCGA_BOTHSTRIP, loss_both_TCGA) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   and SCAN-B:")
variant_table %>% dplyr::arrange(dplyr::desc(rho_SCANB_FULL)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(base, rho_SCANB_FULL, rho_SCANB_MITOSTRIP,
                rho_SCANB_PROLIFSTRIP, rho_SCANB_BOTHSTRIP, loss_both_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   THE PANEL-LEVEL ANSWER:")
survival <- panel %>%
  dplyr::group_by(strip_status) %>%
  dplyr::summarise(
    n_sig = dplyr::n(),
    median_TCGA = stats::median(rho_TCGA), min_TCGA = min(rho_TCGA),
    max_TCGA = max(rho_TCGA), n_ge_0.2_TCGA = sum(rho_TCGA >= 0.2),
    median_SCANB = stats::median(rho_SCANB), min_SCANB = min(rho_SCANB),
    max_SCANB = max(rho_SCANB), n_ge_0.2_SCANB = sum(rho_SCANB >= 0.2),
    .groups = "drop") %>%
  dplyr::arrange(match(strip_status, c("FULL", "MITOSTRIP", "PROLIFSTRIP",
                                       "BOTHSTRIP")))
survival %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   The BOTHSTRIP row is F1's answer. A signature there contains no",
        "\n   MitoCarta gene and no HALLMARK E2F/G2M gene, so its correlation",
        " with\n   OXPHOS is neither self-overlap nor proliferation.")

# The flawed adjustment against the clean strip, side by side.
message("\n   the PROLIF_DISJOINT adjustment against the clean PROLIFSTRIP:")
adj_vs_strip <- a$atlas %>%
  dplyr::filter(arm == FOCUS_ARM, instrument == "gsva", stratum == "all",
                kind == "signature (GSVA)",
                (strip_status == "FULL" & adjusted %in% c("raw", "prolif")) |
                (strip_status == "PROLIFSTRIP" & adjusted == "raw")) %>%
  dplyr::mutate(what = dplyr::case_when(
    strip_status == "FULL" & adjusted == "raw"    ~ "FULL_raw",
    strip_status == "FULL" & adjusted == "prolif" ~ "FULL_prolifADJUSTED",
    TRUE                                          ~ "PROLIFSTRIP_raw")) %>%
  dplyr::select(cohort, base, what, rho) %>%
  tidyr::pivot_wider(names_from = what, values_from = rho)
adj_vs_strip %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(dplyr::desc(FULL_raw)) %>% as.data.frame() %>%
  print(row.names = FALSE)
message("   Adjusting removes shared VARIANCE; stripping removes shared GENES.",
        "\n   Where the two disagree, the adjustment was doing something else.")

# =============================================================================
# 4. Are ELLWOOD and ALFANO informative, or just weak?
# =============================================================================
# Three ways a signature can be weak: too few genes, genes that do not move
# together, or genes that move together but disagree with every other MYC
# signature. Measured on the FULL variants.
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
.coherence <- function(genes, E, inf) {
  g <- inf(genes); if (length(g) < 3L) return(NA_real_)
  C <- stats::cor(t(.rank_rows(E[g, , drop = FALSE])))
  stats::median(C[lower.tri(C)])
}
FULL_SIGS <- full$signature
.agreement <- function(gsva_new) {
  S <- gsva_new[FULL_SIGS, , drop = FALSE]
  C <- stats::cor(t(.rank_rows(S))); diag(C) <- NA
  apply(C, 1L, function(v) stats::median(v, na.rm = TRUE))
}
agree_t <- .agreement(nw$tcga_gsva_new); agree_s <- .agreement(sc$gsva_new)

panel_quality <- full %>%
  dplyr::mutate(
    coh_TCGA  = vapply(MYC_SETS[signature], .coherence, numeric(1), ET, .in_t),
    coh_SCANB = vapply(MYC_SETS[signature], .coherence, numeric(1), ES, .in_s),
    agree_TCGA  = agree_t[signature],
    agree_SCANB = agree_s[signature],
    coverage_TCGA = sd_$coverage$frac[match(paste("TCGA", signature),
                        paste(sd_$coverage$cohort, sd_$coverage$set))])
panel_quality %>% dplyr::arrange(coh_TCGA) %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                dplyr::across(c(coh_TCGA, coh_SCANB, agree_TCGA, agree_SCANB,
                                rho_TCGA, rho_SCANB, coverage_TCGA),
                              ~ round(.x, 3))) %>%
  dplyr::select(base, n, frac_prolif, coverage_TCGA, coh_TCGA, coh_SCANB,
                agree_TCGA, agree_SCANB, rho_TCGA, rho_SCANB) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   If ELLWOOD and ALFANO sit at the BOTTOM of coherence AND",
        " agreement, they\n   are weak signatures rather than informative ones,",
        " and F1's positive\n   entanglement slope was an artefact of including",
        " them.")

# =============================================================================
# 5. What M_b is, and why it behaves differently
# =============================================================================
# TWO STRUCTURAL DIFFERENCES, AND THE FIRST REVERSES THE QUESTION.
#
# (i) The version phase 1 used was the MITO-STRIPPED regulon, and section 2
#     shows the panel's OXPHOS correlation tracks mitochondrial content. So M_b
#     being the weakest estimator may mean it is the least contaminated one.
#     Now that all four variants are scored, the comparison is direct.
#
# (ii) It is SIGNED. CollecTRI carries a mode of regulation per edge and ULM
#     uses it, so a repressed target going DOWN pushes M_b UP. No GSVA signature
#     has any such notion. Section 5.2 splits the regulon by sign and scores
#     each half UNSIGNED to see whether that is what makes M_b move differently.
message("\n5. M_b: mito-stripped, and signed")

mb_variants <- a$atlas %>%
  dplyr::filter(arm == FOCUS_ARM, instrument == "gsva", stratum == "all",
                adjusted == "raw", kind == "CollecTRI regulon (ULM)") %>%
  dplyr::select(myc_estimator, strip_status, cohort, n_genes, rho) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = rho)
message("\n   M_b against ", FOCUS_ARM, ", by variant:")
mb_variants %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   If M_b__FULL is much higher than M_b__MITOSTRIP, then M_b's",
        " weakness in\n   phase 1 was the strip and not the estimator.")

# --- 5.1 the sign rule, reproduced exactly -----------------------------------
# From myc_human_validation script 06 at d3ac60e and E02 section 3.1b:
# mor = if_else(as.logical(is_stimulation), 1, -1). The snapshot is raw
# OmniPath, so the sign lives in is_stimulation / is_inhibition and symbols in
# *_genesymbol. Edges flagged BOTH take +1, silently, and that is reported here
# rather than improved on.
collectri <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE,
                             progress = FALSE)
COLL_ALL <- sd_$collectri_sets[[MB_REF]]
myc_reg <- collectri %>%
  dplyr::filter(source_genesymbol == "MYC", !is.na(target_genesymbol),
                target_genesymbol != "") %>%
  dplyr::transmute(target = target_genesymbol,
                   stim  = as.logical(is_stimulation),
                   inhib = as.logical(is_inhibition),
                   mor   = dplyr::if_else(as.logical(is_stimulation), 1, -1)) %>%
  dplyr::distinct(target, .keep_all = TRUE) %>%
  dplyr::filter(target %in% COLL_ALL)
POS <- sort(unique(myc_reg$target[myc_reg$mor > 0]))
NEG <- sort(unique(myc_reg$target[myc_reg$mor < 0]))
N_AMBIG <- sum(myc_reg$stim & myc_reg$inhib)
message("\n   ", MB_REF, ": ", nrow(myc_reg), " edges | activated ", length(POS),
        " (", N_AMBIG, " of them BOTH-flagged and counted as activating)",
        " | repressed ", length(NEG))

# --- 5.2 the signed halves, scored unsigned ---------------------------------
# LOG matrix (VST), mean of per-gene z - the `zmean` construction from E02
# section 4.4, reused so the halves are comparable to the arms.
.zmean <- function(genes, E, inf) {
  g <- inf(genes); sub <- E[g, , drop = FALSE]
  v <- apply(sub, 1L, stats::var); sub <- sub[v > 0, , drop = FALSE]
  colMeans((sub - rowMeans(sub)) / apply(sub, 1L, stats::sd))
}
.halves <- function(E, inf, ids) {
  h <- rbind(regulon_activated = .zmean(POS, E, inf),
             regulon_repressed = .zmean(NEG, E, inf))
  h <- rbind(h, regulon_difference = h[1, ] - h[2, ])
  colnames(h) <- ids; h
}
halves_t <- .halves(ET, .in_t, ID_T); halves_s <- .halves(ES, .in_s, ID_S)

.probe <- function(halves, gsva_new, arms_obj, mb, ids, coh) {
  TARGETS <- rbind(
    M_b           = as.numeric(mb[ids]),
    MYC_REF_score = as.numeric(gsva_new[MYC_REF, ids]),
    MYC_LOW_ENT   = as.numeric(gsva_new[MYC_LOW_ENTANG, ids]),
    OXPHOS_gsva   = as.numeric(arms_obj$gsva_arms[FOCUS_ARM, ids]),
    MITORIBO_gsva = as.numeric(arms_obj$gsva_arms["Mitochondrial ribosome", ids]),
    PROLIF        = as.numeric(arms_obj$gsva_cov["PROLIF_DISJOINT", ids]))
  colnames(TARGETS) <- ids
  .atlas_block(TARGETS, halves[, ids, drop = FALSE], ids, NULL, min_n = 30L) %>%
    dplyr::rename(target = myc_estimator, half = measure) %>%
    dplyr::mutate(cohort = coh) %>%
    dplyr::select(cohort, half, target, rho, ci_lo, ci_hi)
}
mb_probe <- dplyr::bind_rows(
  .probe(halves_t, nw$tcga_gsva_new, mito, nw$tcga_M_b_variants[MB_REF, ], ID_T, "TCGA"),
  .probe(halves_s, sc$gsva_new,      sc,   sc$M_b_variants[MB_REF, ],      ID_S, "SCAN-B"))
message("\n   what each half of the regulon tracks:")
mb_probe %>%
  tidyr::pivot_wider(id_cols = c(half, target), names_from = cohort,
                     values_from = rho) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 5.3 is M_b's BCL2-family signal independent? ---------------------------
message("\n   the BCL2 family: in the RAW regulon, and in the one M_b uses?")
raw_reg <- collectri %>% dplyr::filter(source_genesymbol == "MYC") %>%
  dplyr::pull(target_genesymbol)
bcl2_in <- tibble::tibble(gene = sd_$bcl2_family) %>%
  dplyr::mutate(in_raw_regulon = gene %in% raw_reg,
                in_M_b_ref     = gene %in% COLL_ALL,
                is_mitocarta   = gene %in% MITO_ALL)
bcl2_in %>% as.data.frame() %>% print(row.names = FALSE)
if (any(bcl2_in$in_M_b_ref)) {
  stop("a BCL2-family gene survives into ", MB_REF, ": ",
       paste(bcl2_in$gene[bcl2_in$in_M_b_ref], collapse = ", "),
       ". Its correlation with that gene is partly self-referential.",
       call. = FALSE)
}
message("   None survives into ", MB_REF, " - all 15 are MitoCarta genes and",
        " the\n   mitochondrial strip removed every one that was in the regulon.",
        " So M_b's\n   BCL2-family numbers are independent observations. NOTE",
        " that M_b__FULL is\n   NOT: it retains them, and must not be used for",
        " a BCL2-family claim.")

message("\n   and the FULL GSVA signatures, which were never stripped:")
bcl2_in_sigs <- tibble::tibble(signature = FULL_SIGS) %>%
  dplyr::mutate(n_bcl2 = vapply(MYC_SETS[signature],
                                function(g) sum(g %in% sd_$bcl2_family),
                                integer(1)),
                which = vapply(MYC_SETS[signature],
                               function(g) paste(intersect(g, sd_$bcl2_family),
                                                 collapse = ","), character(1))) %>%
  dplyr::filter(n_bcl2 > 0)
if (nrow(bcl2_in_sigs)) {
  bcl2_in_sigs %>% as.data.frame() %>% print(row.names = FALSE)
  message("   Those cells of E05's overlay are NOT independent observations.",
          "\n   Their __MITOSTRIP variants are.")
} else message("   none.")

# =============================================================================
# 6. Save
# =============================================================================
message("\n6. save")
saveRDS(list(panel = panel, d0_test = d0_test, variant_table = variant_table,
             survival = survival, adj_vs_strip = adj_vs_strip,
             panel_quality = panel_quality, mb_variants = mb_variants,
             regulon = list(activated = POS, repressed = NEG,
                            n_both_flagged = N_AMBIG, table = myc_reg),
             mb_probe = mb_probe, bcl2_in_regulon = bcl2_in,
             bcl2_in_signatures = bcl2_in_sigs,
             halves = list(TCGA = halves_t, `SCAN-B` = halves_s),
             rules = list(
               f1_answer = paste("the BOTHSTRIP row of `survival` is F1's",
                                 "answer: those signatures contain no MitoCarta",
                                 "gene and no HALLMARK E2F/G2M gene"),
               adjust_vs_strip = paste("adjusting removes shared VARIANCE,",
                                       "stripping removes shared GENES;",
                                       "PROLIF_DISJOINT is disjoint from M_a",
                                       "alone and shares 97 genes with the",
                                       "FULL panel union"),
               mb_bcl2 = paste("no BCL2-family gene survives into M_b__MITOSTRIP",
                               "so its BCL2 numbers are independent; M_b__FULL",
                               "retains them and must not be used for a",
                               "BCL2-family claim")),
             built = Sys.time()), PATH_E06)
message("\nE06: done.  results/estimator_anatomy.rds")

# =============================================================================
# Sandbox
# =============================================================================
if (FALSE) {

  e <- readRDS(PATH_E06)

  # F1's answer, in one table
  e$survival %>% as.data.frame()

  # per signature, what each strip costs
  e$variant_table %>% dplyr::arrange(dplyr::desc(loss_both_TCGA)) %>%
    dplyr::select(base, rho_TCGA_FULL, rho_TCGA_MITOSTRIP,
                  rho_TCGA_PROLIFSTRIP, rho_TCGA_BOTHSTRIP) %>% as.data.frame()

  # adjustment vs stripping - do they agree?
  e$adj_vs_strip %>% as.data.frame()

  # is M_b weak because it was stripped?
  e$mb_variants %>% as.data.frame()
  e$mb_probe %>% tidyr::pivot_wider(id_cols = c(half, target),
                                    names_from = cohort, values_from = rho) %>%
    as.data.frame()

  # are the low-entanglement outliers weak?
  e$panel_quality %>% dplyr::arrange(coh_TCGA) %>%
    dplyr::select(base, n, coh_TCGA, agree_TCGA, rho_TCGA) %>% head(8) %>%
    as.data.frame()

  e$bcl2_in_signatures %>% as.data.frame()

}
