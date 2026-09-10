# E24_mitopps_rank_comparison.R
# =============================================================================
# WHERE DOES OXPHOS SIT WITHIN THE MITOCHONDRIAL PANEL, IN MOUSE GLAND AND IN
# HUMAN TUMOURS?
#
# Development lowers the priority of the respiratory chain within the
# mitochondrial compartment in mouse mammary epithelium between 6W and 12W. This
# asks whether human breast tumours sit at that LOWERED position or back at the
# PRE-WINDOW one.
#
# =============================================================================
# RANK, NEVER VALUE. THIS IS THE RULE THE WHOLE SCRIPT IS BUILT AROUND.
# =============================================================================
# mitoPPS scores every mitochondrial pathway WITHIN a sample, so "where does
# OXPHOS sit among the pathways in this sample" is a within-sample ORDERING.
# Comparing orderings across cohorts is permitted; comparing mitoPPS VALUES
# across species is forbidden - the baseline is composition-dependent and a
# species is a cohort (CLAUDE.md traps 3 and 6).
#
# THEREFORE: this script carries PERCENTILE RANK ONLY. No raw mitoPPS value is
# saved into any element of the output object, and no cross-species value
# appears in any table or figure. Section 0 fixes the percentile definition
# before anything is computed.
#
# =============================================================================
# EXPLORATORY, AND IT STAYS EXPLORATORY
# =============================================================================
# Post-hoc relative to the pre-registered work. BOTH outcomes are pre-declared
# in the reading rules below and NEITHER IS PREDICTED. Nothing here is a
# hypothesis test, no p-value is computed for the placement, and the note must
# say so in its own words.
#
# N3. These are transcript associations. The word "primed" is not written of a
# transcript anywhere in this script, its comments, its labels or its outputs.
#
# =============================================================================
# GATE 0 - ANSWERED BEFORE THIS FILE WAS WRITTEN. THE PANELS ARE THE SAME OBJECT
# =============================================================================
#                                  human (both cohorts)   mouse
#   pathways before filtering      150                    152   (+2, explained)
#   after filtering                142                    144   (same +2)
#   filter                         min_genes = 3L         the IDENTICAL 8 drop
#   MT-* / mt-*                    one synthetic pathway  one synthetic pathway
#                                  `mtDNA-encoded OXPHOS subunits`, n = 13, both
#   MT-* inside OXPHOS subunits    none                   none
#   MT-* inside the OXPHOS umbrella none                  none
#   `OXPHOS subunits`              89                     89    (exact)
#   input scale                    linear DESeq2-normalised, both
#
# EVERY one of the human's 150 pre-filter pathways is present in the mouse's
# 152, and THE SAME EIGHT DROP IN BOTH - Cholesterol-associated, Cytochrome C,
# Glycerol phosphate shuttle, mtDNA modifications, OXA, Vitamin B1 / B6 / C
# metabolism. That the filter lands on an identical set is stronger evidence of
# a shared rule than the threshold value.
#
# The mt- question was answered AGAINST THE SAVED OBJECT, never against a copy
# of the mouse script - the strip postdates the script's first version, so an
# old script reads wrong with no sign that it is wrong.
# `docs/2026-09-10_mouse_corrections.md` section 7.2.
#
# TWO RESIDUAL ONTOLOGY DIFFERENCES, recorded and immaterial: the OXPHOS
# umbrella is 156 against 155 and OXPHOS assembly factors 68 against 67. The
# PRIMARY set is identical at 89, so G1 is not downgraded on their account.
#
# THE +2 is `Apoptosis-PRO` (25) and `Apoptosis-ANTI` (9), a disjoint partition
# of mouse `Apoptosis` (34) - REDUNDANT panel members. The intersection drops
# them automatically, which is why PART A is primary.
#
# =============================================================================
# READING RULES - FIXED HERE, BEFORE ANY NUMBER IN THIS SCRIPT IS SEEN
# =============================================================================
# G1  PANEL COMPARABILITY. Answered at GATE 0 above: PASSED. Recorded as a
#     constant so section 8's verdict cannot drift from it.
#
# G2  THE SIZE GATE, and it runs BEFORE any placement is read. `OXPHOS subunits`
#     carries 89 members against 10-30 for most pathways, and a large set can
#     rank differently for structural reasons. If |rho(rank, log size)| >= 0.5
#     in EITHER human cohort, OXPHOS's rank is substantially structural and the
#     comparison is readable only against SIZE-MATCHED pathways, not against the
#     panel. Report it that way or report that it cannot be read.
#
#     THE REMEDY IS PRE-SPECIFIED HERE, so that firing G2 does not become an
#     invitation to improvise. If G2 fires, R1 is scored on a SIZE-MATCHED
#     panel: the pathways of the intersection whose size is within a factor f
#     of the primary's, re-ranked among themselves. The band is defined ONCE on
#     the human resolved sizes and the SAME pathway list is used for both
#     species, so the comparator set is identical across the comparison. TWO
#     widths (f = 1.5 and f = 2) are reported so that the width, which is a
#     post-hoc choice, is visible and can be seen not to matter. If the wider
#     band holds fewer than 10 pathways, the placement CANNOT BE READ and that
#     is the report.
#
# G3  THE MOUSE INTERVAL MUST EXIST. If OXPHOS's percentile does not fall from
#     6W WT to 12W WT there is no reprioritisation IN RANK TERMS and the human
#     placement has nothing to be placed against. THE ANALYSIS STOPS THERE.
#
#     G3 IS GENUINELY OPEN, NOT A FORMALITY. The mouse developmental result was
#     established on COMPOSITES against expression-matched nulls, NOT as a rank
#     within the mitochondrial panel. It has never been shown to survive that
#     translation, and this is the first time it is asked to.
#
# R1  THE PLACEMENT. BOTH READINGS PRE-DECLARED, NEITHER PREDICTED.
#       near the 12W position -> consistent with the lowered priority persisting
#         into established tumours
#       near the 6W position  -> the proportions are restored in tumours and the
#         developmental change does not persist
#       between, or outside both -> REPORTED AS SUCH, never assigned to the
#         nearer end
#     Whichever it is, it is ONE DESCRIPTIVE OBSERVATION, NOT A TEST.
#
# R2  COHORT AGREEMENT. TCGA and SCAN-B reported separately and NEVER averaged
#     or pooled. If they disagree, THAT IS THE RESULT.
#
# NO RULE IS PRE-SPECIFIED FOR THE APOPTOSIS READOUT (PART A2). It is reported
# and not scored, for the reason in the next block.
#
# =============================================================================
# PART A2, AND WHY IT IS A THIRD LEG RATHER THAN PART OF THE PRIMARY
# =============================================================================
# The human side gets `Apoptosis-PRO` / `Apoptosis-ANTI` too, so the two
# apoptosis arms can be seen to move. Three things travel with it and all three
# are stated before any number.
#
#  (1) IT REQUIRES RECOMPUTING HUMAN mitoPPS. `functions/mitopps.R` puts the
#      universe size P in the denominator of BOTH `.mitopps_universe` (P - 1)
#      and `.mitopps_query` (P). Adding two pathways changes EVERY human score.
#      The saved 142-row `$mitopps_universe` cannot be extended.
#      SO PART A2 CARRIES ITS OWN CONTROL AND IT RUNS FIRST: rebuild the human
#      142 universe from the linear matrices and assert it reproduces the SAVED
#      object bit-equal, in BOTH cohorts, before any 144 number is read. If that
#      fails the recompute machinery is wrong and nothing in A2 is readable.
#
#  (2) IT IS HUMAN-NATIVE, so no rule is broken. The split comes from
#      `data/genesets_celldeath_human/cell_death_genes_consolidated.csv`, which
#      carries first-class human `effect` calls. NO ORTHOLOG FUNCTION IS CALLED
#      and no mouse-derived gene set is loaded.
#
#  (3) IT IS NOT THE MOUSE'S SET, and this is why it is not scored.
#        curation   mouse MC_Apoptosis_Pro/_Anti   human `effect` column
#        PRO        25                            25   (agrees exactly)
#        ANTI        9                             6   (does not)
#        coverage   34 of 34, 100 pct             31 of 35, 89 pct
#      An independently curated ANALOGUE, not a translation. The two splits are
#      comparable in CONSTRUCTION RULE but not in MEMBERSHIP, so any apoptosis
#      rank read across species is DESCRIPTIVE ONLY.
#
# =============================================================================
# SCALES AND INPUTS
# =============================================================================
#   human mitoPPS      READ from the saved objects; and REBUILT in A2 from
#                      linear DESeq2-normalised counts, asserted
#   mouse mitoPPS      READ ONLY from data/from_myc_mouse/mitopps_scores.rds,
#                      md5 checked against that directory's README.
#                      NEVER RECOMPUTED HERE. The mouse repo is not read.
#   the apoptosis sets human-native, from the cell-death snapshot
#
# NO POOLING of mouse and human samples into one scoring run, ever. Cohorts are
# scored separately and TCGA and SCAN-B are never averaged. No FDR across
# pathways - this is a ranking, not a screen; estimates with intervals.
# SPECIES: the human side is human. No ortholog function is called here or
# anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "mitopps.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE24: where does OXPHOS sit in the mitochondrial panel?\n", strrep("=", 78))

PATH_E24     <- file.path(DIR_RESULTS, "mitopps_rank_comparison.rds")
PATH_E24_TAB <- file.path(DIR_TABLES,  "E24_ranks_by_group.csv")
DIR_E24_FIG  <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_MOUSE   <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")
PATH_MOUSE_README <- here::here("data", "from_myc_mouse", "README.md")

set.seed(1)
NBOOT <- 5000L

# =============================================================================
# 0. CONSTANTS AND THE PERCENTILE DEFINITION, FIXED BEFORE ANY NUMBER
# =============================================================================
PRIMARY   <- "OXPHOS subunits"          # what the mouse developmental result is on
SECONDARY <- "OXPHOS"                   # the umbrella; includes assembly factors
APO_READ  <- c("Apoptosis-PRO", "Apoptosis-ANTI")   # A2 only, descriptive
MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"
CUT_SIZE  <- 0.5    # G2: |rho(rank, log size)| at which the rank is structural

# GATE 0's answer, transcribed. A constant so section 8 cannot drift from it.
G1_PASSED <- TRUE
G1_NOTE <- paste0(
  "PASSED. Same catalogue; the IDENTICAL 8 pathways drop under min_genes = 3 in ",
  "both species; mt-/MT- stripped into one synthetic pathway of 13 in both with ",
  "none inside OXPHOS subunits or the umbrella; OXPHOS subunits 89 in both; ",
  "linear DESeq2-normalised input in both. Two residual ontology differences ",
  "(umbrella 156/155, assembly factors 68/67) sit OUTSIDE the primary set and ",
  "do NOT downgrade the placements to descriptive.")

# THE PERCENTILE, defined once. Higher mitoPPS = higher priority within the
# sample, so percentile 100 is the MOST prioritised pathway in that sample.
# Mid-rank form, so the endpoints are not 0 and 100 and the measure is symmetric.
.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}
PCT_NOTE <- "100 * (rank - 0.5) / P within each sample; 100 = most prioritised"

message("\n0. rules fixed before any number")
message("   G1  ", G1_NOTE)
message("   G2  |rho(rank, log size)| >= ", CUT_SIZE,
        " in EITHER human cohort -> the rank is structural")
message("   G3  if OXPHOS's percentile does not FALL 6W WT -> 12W WT, STOP")
message("   R1  near 12W / near 6W / between-or-outside; neither predicted")
message("   R2  TCGA and SCAN-B separate, never averaged")
message("   pct ", PCT_NOTE)

# =============================================================================
# 1. INPUTS
# =============================================================================
message("\n1. inputs")

# --- 1.1 the mouse artefact, md5-checked against its own README --------------
if (!file.exists(PATH_MOUSE)) {
  stop("data/from_myc_mouse/mitopps_scores.rds is missing. It is a snapshot; ",
       "see that directory's README for the source path.", call. = FALSE)
}
md5_now <- unname(tools::md5sum(PATH_MOUSE))
if (!identical(md5_now, MOUSE_MD5)) {
  stop("the mouse artefact's md5 is ", md5_now, ", not the ", MOUSE_MD5,
       " recorded in data/from_myc_mouse/README.md. The snapshot has changed ",
       "under this script. Re-read the README before going further.",
       call. = FALSE)
}
mo <- readRDS(PATH_MOUSE)
message("   mouse artefact md5 matches its README: ", MOUSE_MD5)

# --- 1.2 the human objects ---------------------------------------------------
m  <- readRDS(PATH_TCGA_MITO)
sc <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
ID_T <- colnames(m$gsva_arms); ID_S <- colnames(sc$gsva_arms)
stopifnot(nrow(m$mitopps_universe) == 142L, nrow(sc$mitopps_universe) == 142L,
          identical(rownames(m$mitopps_universe), rownames(sc$mitopps_universe)))
message("   human mitopps_universe: TCGA ", ncol(m$mitopps_universe),
        " samples | SCAN-B ", ncol(sc$mitopps_universe),
        " samples | 142 pathways, identical rownames")

# --- 1.3 the mouse panel and its annotation ---------------------------------
MOUSE_ANN <- c("sample", "group", "timepoint", "myc_status")
mouse_paths <- setdiff(names(mo$mitopps_scores), MOUSE_ANN)
stopifnot(length(mouse_paths) == mo$n_pathways, mo$n_pathways == 144L)
mo_ann <- mo$mitopps_scores[, MOUSE_ANN]
# WT is `neg`, Myc+ is `pos`. Named here so no downstream line has to guess.
mo_ann$arm <- dplyr::case_when(
  mo_ann$group == "6W_neg"  ~ "mouse 6W WT",
  mo_ann$group == "12W_neg" ~ "mouse 12W WT",
  mo_ann$group == "6W_pos"  ~ "mouse 6W Myc+",
  mo_ann$group == "12W_pos" ~ "mouse 12W Myc+",
  TRUE ~ NA_character_)
stopifnot(!anyNA(mo_ann$arm))
message("   mouse: ", nrow(mo_ann), " samples, 144 pathways, groups ",
        paste(unique(mo_ann$group), collapse = ", "))

# =============================================================================
# 2. PART 0 - PANEL RECONCILIATION, SAVED SO THE COMPARISON CARRIES PROVENANCE
# =============================================================================
message("\n2. PART 0: panel reconciliation")

hum_paths_all <- names(m$mito_paths)                       # 150, pre-filter
hum_paths     <- rownames(m$mitopps_universe)               # 142, scored
mou_paths_all <- unique(mo$gene_to_pathway$Pathway)         # 152, pre-filter
mou_paths     <- mouse_paths                                # 144, scored

hsz <- vapply(m$mito_paths, function(g) length(unique(g)), integer(1))
msz <- tapply(mo$gene_to_pathway$Gene, mo$gene_to_pathway$Pathway,
              function(v) length(unique(v)))

panel_reconciliation <- tibble::tibble(
  row = c("pathways before filtering", "pathways after filtering",
          "the filter", "dropped by the filter",
          "mt- / MT- handling", "mt- inside OXPHOS subunits",
          "OXPHOS subunits size", "OXPHOS umbrella size",
          "OXPHOS assembly factors size", "input scale"),
  human = c(as.character(length(hum_paths_all)), as.character(length(hum_paths)),
            "min_genes = 3L", paste(sort(setdiff(hum_paths_all, hum_paths)), collapse = ", "),
            "one synthetic pathway, n = 13", "none",
            as.character(hsz[["OXPHOS subunits"]]), as.character(hsz[["OXPHOS"]]),
            as.character(hsz[["OXPHOS assembly factors"]]),
            "linear DESeq2-normalised"),
  mouse = c(as.character(length(mou_paths_all)), as.character(length(mou_paths)),
            "the identical 8 drop", paste(sort(setdiff(mou_paths_all, mou_paths)), collapse = ", "),
            "one synthetic pathway, n = 13", "none",
            as.character(msz[["OXPHOS subunits"]]), as.character(msz[["OXPHOS"]]),
            as.character(msz[["OXPHOS assembly factors"]]),
            "linear DESeq2-normalised"))
panel_reconciliation %>% as.data.frame() %>% print(row.names = FALSE)

# The gate's own assertions, re-run here rather than trusted from the header.
stopifnot(
  identical(sort(setdiff(hum_paths_all, hum_paths)),
            sort(setdiff(mou_paths_all, mou_paths))),          # same 8 drop
  msz[["OXPHOS subunits"]] == hsz[["OXPHOS subunits"]],        # 89 == 89
  !any(grepl("^MT-", m$mito_paths[["OXPHOS subunits"]])),
  !any(grepl("^mt-", mo$gene_to_pathway$Gene[
    mo$gene_to_pathway$Pathway == "OXPHOS subunits"])))
message("   GATE 0 re-asserted in code: the same 8 drop, OXPHOS subunits 89 = 89,",
        "\n   no mt-/MT- inside the primary set in either species.")

# =============================================================================
# 3. PART A - THE INTERSECTED PANEL. THIS IS THE PRIMARY CONSTRUCTION.
# =============================================================================
# Ranking within each species' OWN full panel lets panel composition drive the
# comparison. The intersection removes that. Full-panel ranks are computed too,
# as a NAMED SENSITIVITY, and reported beside the primary - never merged with it.
message("\n3. PART A: the intersected panel")

shared <- intersect(hum_paths, mou_paths)
intersection <- list(
  n_shared = length(shared),
  shared = sort(shared),
  human_only = sort(setdiff(hum_paths, mou_paths)),
  mouse_only = sort(setdiff(mou_paths, hum_paths)))
message("   shared: ", length(shared),
        " | human-only: ", length(intersection$human_only),
        " | mouse-only: ", paste(intersection$mouse_only, collapse = ", "))
if (length(shared) < 15L) {
  stop("the intersected panel is ", length(shared), " pathways. Percentiles over ",
       "a panel this small are too coarse to read. Stop.", call. = FALSE)
}
if (length(shared) < 100L) {
  message("   *** WARNING: the intersection is below 100 pathways. That is itself",
          "\n   *** a finding - report the dropped sets and what drove it (build",
          "\n   *** version, species ontology, filtering) before reading anything.")
}
stopifnot(PRIMARY %in% shared, SECONDARY %in% shared)

# --- the rank engine ---------------------------------------------------------
# HUMAN: mitopps_universe is pathways x samples, so rank down each COLUMN.
# MOUSE: mitopps_scores is samples x pathways, so rank across each ROW.
.pct_human <- function(U, keep) {
  U <- U[keep, , drop = FALSE]
  out <- apply(U, 2L, .pct)
  rownames(out) <- rownames(U); out                     # pathways x samples
}
.pct_mouse <- function(df, keep) {
  M <- as.matrix(df[, keep, drop = FALSE])
  out <- t(apply(M, 1L, .pct))
  colnames(out) <- keep; out                            # samples x pathways
}

# --- primary: within the intersection ---------------------------------------
pct_T_int <- .pct_human(m$mitopps_universe,  shared)
pct_S_int <- .pct_human(sc$mitopps_universe, shared)
pct_M_int <- .pct_mouse(mo$mitopps_scores,   shared)
# --- sensitivity: each species' OWN full panel -------------------------------
pct_T_full <- .pct_human(m$mitopps_universe,  hum_paths)
pct_S_full <- .pct_human(sc$mitopps_universe, hum_paths)
pct_M_full <- .pct_mouse(mo$mitopps_scores,   mou_paths)
message("   ranked: intersection (", length(shared), ") and each species' own ",
        "full panel (human ", length(hum_paths), ", mouse ", length(mou_paths), ")")

# =============================================================================
# 4. PART A2 - THE APOPTOSIS-SPLIT PANEL. ITS CONTROL RUNS FIRST.
# =============================================================================
message("\n4. PART A2: the apoptosis-split panel")

# --- 4.1 THE CONTROL. Rebuild the human 142 universe and prove it bit-equal. --
# Nothing in A2 is readable unless this passes, because A2's 144 run uses the
# same machinery with two extra rows.
message("   4.1 control: rebuilding the human 142 universe from the linear matrices")
tl <- readRDS(PATH_TCGA_LINEAR)
stopifnot(identical(tl$scale, "linear_deseq2_normalised"))
Lt <- tl$mat[, ID_T, drop = FALSE]; rm(tl); invisible(gc(verbose = FALSE))
sl <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(sl$scale, "linear_deseq2_normalised"))
Ls <- sl$mat[, ID_S, drop = FALSE]; rm(sl); invisible(gc(verbose = FALSE))

.remap_s <- function(g) { h <- sc$symbol_map[g]; unname(ifelse(is.na(h), g, h)) }
paths_T <- lapply(m$mito_paths, function(g) intersect(unique(g), rownames(Lt)))
paths_S <- lapply(m$mito_paths, function(g) intersect(.remap_s(unique(g)), rownames(Ls)))

.rebuild <- function(paths, L) .mitopps_universe(.path_scores(paths, L, MIN_SET_GENES))
U_T142 <- .rebuild(paths_T, Lt)
U_S142 <- .rebuild(paths_S, Ls)
ctrl <- tibble::tibble(
  cohort = c("TCGA", "SCAN-B"),
  n_pathways = c(nrow(U_T142), nrow(U_S142)),
  max_abs_delta = c(
    max(abs(U_T142 - m$mitopps_universe[rownames(U_T142), colnames(U_T142)])),
    max(abs(U_S142 - sc$mitopps_universe[rownames(U_S142), colnames(U_S142)]))))
ctrl %>% as.data.frame() %>% print(row.names = FALSE)
if (max(ctrl$max_abs_delta) > 1e-10) {
  stop("the rebuilt human 142 universe does not reproduce the saved object ",
       "(max |delta| = ", signif(max(ctrl$max_abs_delta), 3), "). The recompute ",
       "machinery is wrong and PART A2 is unreadable. Stop.", call. = FALSE)
}
message("   4.1 CONTROL PASSED: both cohorts reproduce bit-equal (max |delta| = ",
        signif(max(ctrl$max_abs_delta), 3), ")")

# --- 4.2 the human-native apoptosis split ------------------------------------
# From the cell-death snapshot's first-class human `effect` column. NO ORTHOLOG
# FUNCTION and no mouse-derived set. It is an independently curated ANALOGUE of
# the mouse split, not a translation - see the header.
cd <- readr::read_csv(
  here::here("data", "genesets_celldeath_human", "cell_death_genes_consolidated.csv"),
  show_col_types = FALSE, progress = FALSE)
apo <- m$mito_paths[["Apoptosis"]]
apo_pro  <- intersect(apo, cd$human_symbol[cd$effect == "pro-death"])
apo_anti <- intersect(apo, cd$human_symbol[cd$effect == "pro-survival"])
apoptosis_split <- tibble::tibble(
  set = c("Apoptosis", "Apoptosis-PRO", "Apoptosis-ANTI",
          "unclassified or unannotated"),
  human_n = c(length(apo), length(apo_pro), length(apo_anti),
              length(apo) - length(apo_pro) - length(apo_anti)),
  mouse_n = c(msz[["Apoptosis"]], msz[["Apoptosis-PRO"]], msz[["Apoptosis-ANTI"]],
              msz[["Apoptosis"]] - msz[["Apoptosis-PRO"]] - msz[["Apoptosis-ANTI"]]),
  note = c("the MitoCarta pathway", "human: effect == pro-death",
           "human: effect == pro-survival",
           "human only; the mouse split is a complete partition"))
apoptosis_split %>% as.data.frame() %>% print(row.names = FALSE)
stopifnot(length(apo_pro) >= MIN_SET_GENES, length(apo_anti) >= MIN_SET_GENES,
          length(intersect(apo_pro, apo_anti)) == 0L)
message("   the human split is an ANALOGUE, not a translation: PRO agrees at ",
        length(apo_pro), ", ANTI does not (", length(apo_anti), " against ",
        msz[["Apoptosis-ANTI"]], "). Descriptive only, no rule.")

# --- 4.3 the human 144 universe ---------------------------------------------
paths_T144 <- c(paths_T, list(`Apoptosis-PRO`  = intersect(apo_pro,  rownames(Lt)),
                              `Apoptosis-ANTI` = intersect(apo_anti, rownames(Lt))))
paths_S144 <- c(paths_S, list(`Apoptosis-PRO`  = intersect(.remap_s(apo_pro),  rownames(Ls)),
                              `Apoptosis-ANTI` = intersect(.remap_s(apo_anti), rownames(Ls))))
U_T144 <- .rebuild(paths_T144, Lt)
U_S144 <- .rebuild(paths_S144, Ls)
stopifnot(nrow(U_T144) == 144L, nrow(U_S144) == 144L,
          all(APO_READ %in% rownames(U_T144)), all(APO_READ %in% rownames(U_S144)))
message("   4.3 human universes rebuilt at 144 pathways, both cohorts")
rm(Lt, Ls); invisible(gc(verbose = FALSE))

# how far the recompute moved the primary's percentile - the cost of adding two
a2_shift <- tibble::tibble(
  cohort = c("TCGA", "SCAN-B"),
  median_pct_142 = c(stats::median(.pct_human(m$mitopps_universe, hum_paths)[PRIMARY, ]),
                     stats::median(.pct_human(sc$mitopps_universe, hum_paths)[PRIMARY, ])),
  median_pct_144 = c(stats::median(.pct_human(U_T144, rownames(U_T144))[PRIMARY, ]),
                     stats::median(.pct_human(U_S144, rownames(U_S144))[PRIMARY, ]))) %>%
  dplyr::mutate(shift = median_pct_144 - median_pct_142)
message("   what adding two pathways costs the primary's percentile:")
a2_shift %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

pct_T_a2 <- .pct_human(U_T144, rownames(U_T144))
pct_S_a2 <- .pct_human(U_S144, rownames(U_S144))
pct_M_a2 <- pct_M_full     # the mouse already has its 144

# =============================================================================
# 5. PART B - THE SIZE GATE (G2). RUNS BEFORE ANY PLACEMENT IS READ.
# =============================================================================
# Human only, because G2 is written on the human cohorts. The mouse is reported
# beside it and is DESCRIPTIVE: `gene_to_pathway` is catalogue membership, not
# detected genes, so its sizes are not the same quantity as the human ones.
message("\n5. PART B: the size gate")

# TWO DIFFERENT QUANTITIES, and conflating them is easy, so both are computed
# and the gate is scored on the one the rule names.
#
#   rho_panel  -- G2's STATISTIC. The unit is the PATHWAY: each pathway's TYPICAL
#                 rank (its mean percentile across the cohort's samples) against
#                 its log size, Spearman over the whole intersected panel.
#                 Bootstrapped over PATHWAYS, which is the unit.
#   rho_within -- reported beside it, NOT the gate. The mean over samples of the
#                 within-sample rho. It answers a different question - "in any
#                 ONE sample, does size predict rank" - and it is attenuated by
#                 per-sample noise, so it is systematically the smaller number.
#
# The rule says "across the whole intersected panel", so `rho_panel` scores G2.
.size_gate <- function(P, sizes, lab) {
  s <- sizes[rownames(P)]
  stopifnot(!anyNA(s))
  mp <- rowMeans(P)
  rho_panel  <- stats::cor(mp, log(s), method = "spearman")
  rho_within <- mean(apply(P, 2L, function(col)
    stats::cor(col, log(s), method = "spearman")))
  bs <- vapply(seq_len(NBOOT), function(i) {
    k <- sample.int(length(s), replace = TRUE)
    if (length(unique(s[k])) < 3L) return(NA_real_)
    stats::cor(mp[k], log(s)[k], method = "spearman") }, 0)
  tibble::tibble(panel = lab, n_pathways = nrow(P),
                 rho_panel = rho_panel,
                 ci_lo = unname(stats::quantile(bs, 0.025, na.rm = TRUE)),
                 ci_hi = unname(stats::quantile(bs, 0.975, na.rm = TRUE)),
                 rho_within_sample_mean = rho_within,
                 structural = abs(rho_panel) >= CUT_SIZE) }

hsz_used <- vapply(paths_T, length, integer(1))   # sizes AS SCORED, not catalogue
size_gate <- dplyr::bind_rows(
  .size_gate(pct_T_int, hsz_used, "TCGA, intersected panel"),
  .size_gate(pct_S_int, hsz_used, "SCAN-B, intersected panel"))
size_gate %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
G2_STRUCTURAL <- any(size_gate$structural)
message("   G2: ", if (G2_STRUCTURAL)
  "*** |rho| >= 0.5 in at least one cohort - the rank is SUBSTANTIALLY STRUCTURAL. Read against size-matched pathways or report that it cannot be read."
  else "passed; no cohort reaches |rho| = 0.5, so the rank is not substantially structural.")
message("   where the primary sits by size: ", PRIMARY, " has ",
        hsz_used[[PRIMARY]], " genes as scored, against a panel median of ",
        stats::median(hsz_used[shared]), " -- it is the ",
        sprintf("%.0f", 100 * mean(hsz_used[shared] <= hsz_used[[PRIMARY]])),
        "th percentile of the panel BY SIZE.")
if (G2_STRUCTURAL) {
  message("   *** G2 FAILED. Size predicts a pathway's typical rank across the",
          "\n   *** panel, and the primary is one of the largest sets in it. The",
          "\n   *** placement below is READABLE ONLY AGAINST SIZE-MATCHED",
          "\n   *** PATHWAYS, not against the panel. Section 9 records this and",
          "\n   *** the note must carry it ahead of any R1 sentence.")
}

# --- PART B2: the size-matched panel, the remedy G2's rule prescribes --------
# Built unconditionally so it is in the object either way, and SCORED only if
# G2 fired. The band is defined on the HUMAN resolved sizes and the same
# pathway list is applied to both species, so the comparator set is identical.
SIZE_BANDS <- c(1.5, 2)
ref_size <- hsz_used[[PRIMARY]]
.band <- function(f) {
  keep <- shared[hsz_used[shared] >= ref_size / f & hsz_used[shared] <= ref_size * f]
  sort(keep) }
bands <- lapply(SIZE_BANDS, .band)
names(bands) <- sprintf("size-matched (f=%.1f, n=%d)", SIZE_BANDS, lengths(bands))
message("\n   size-matched bands around ", PRIMARY, " (", ref_size, " genes):")
for (i in seq_along(bands))
  message("     ", names(bands)[i], "  range ",
          min(hsz_used[bands[[i]]]), "-", max(hsz_used[bands[[i]]]), " genes")
BAND_MAIN <- names(bands)[which.max(lengths(bands))]
if (max(lengths(bands)) < 10L) {
  message("   *** the widest size-matched band holds ", max(lengths(bands)),
          " pathways. THE PLACEMENT CANNOT BE READ.")
}
BAND_READABLE <- max(lengths(bands)) >= 10L

# =============================================================================
# 6. PART C - RANKS PER SAMPLE, THEN SUMMARISED BY GROUP
# =============================================================================
# Percentile computed WITHIN EACH SAMPLE, then summarised. Group means are NEVER
# ranked: that discards the spread and the comparison has no error bar without
# it.
message("\n6. PART C: ranks per sample, then by group")

.long_human <- function(P, cohort, panel, sets) {
  dplyr::bind_rows(lapply(sets, function(s) if (!s %in% rownames(P)) NULL else
    tibble::tibble(arm = cohort, panel = panel, pathway = s,
                   sample = colnames(P), pct = as.numeric(P[s, ])))) }
.long_mouse <- function(P, panel, sets) {
  dplyr::bind_rows(lapply(sets, function(s) if (!s %in% colnames(P)) NULL else
    tibble::tibble(arm = mo_ann$arm, panel = panel, pathway = s,
                   sample = mo_ann$sample, pct = as.numeric(P[, s])))) }

SETS_MAIN <- c(PRIMARY, SECONDARY)
ranks_per_sample <- dplyr::bind_rows(
  .long_human(pct_T_int,  "TCGA",   "intersected (142)", SETS_MAIN),
  .long_human(pct_S_int,  "SCAN-B", "intersected (142)", SETS_MAIN),
  .long_mouse(pct_M_int,            "intersected (142)", SETS_MAIN),
  .long_human(pct_T_full, "TCGA",   "own full panel",    SETS_MAIN),
  .long_human(pct_S_full, "SCAN-B", "own full panel",    SETS_MAIN),
  .long_mouse(pct_M_full,           "own full panel",    SETS_MAIN),
  # A2: the apoptosis readout, descriptive, plus the primary on the same panel
  .long_human(pct_T_a2,   "TCGA",   "apoptosis-split (144)", c(SETS_MAIN, APO_READ)),
  .long_human(pct_S_a2,   "SCAN-B", "apoptosis-split (144)", c(SETS_MAIN, APO_READ)),
  .long_mouse(pct_M_a2,             "apoptosis-split (144)", c(SETS_MAIN, APO_READ)),
  # PART B2's size-matched panels, one per band
  dplyr::bind_rows(lapply(seq_along(bands), function(i) {
    bn <- names(bands)[i]; bk <- bands[[i]]
    dplyr::bind_rows(
      .long_human(.pct_human(m$mitopps_universe,  bk), "TCGA",   bn, PRIMARY),
      .long_human(.pct_human(sc$mitopps_universe, bk), "SCAN-B", bn, PRIMARY),
      .long_mouse(.pct_mouse(mo$mitopps_scores,   bk),           bn, PRIMARY)) })))

ARM_ORDER <- c("mouse 6W WT", "mouse 12W WT", "mouse 6W Myc+", "mouse 12W Myc+",
               "TCGA", "SCAN-B")
ranks_per_sample <- ranks_per_sample %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER))
stopifnot(!anyNA(ranks_per_sample$arm))

.boot_median <- function(v) {
  b <- vapply(seq_len(NBOOT), function(i) stats::median(sample(v, replace = TRUE)), 0)
  stats::quantile(b, c(0.025, 0.975), names = FALSE) }
ranks_by_group <- ranks_per_sample %>%
  dplyr::group_by(panel, pathway, arm) %>%
  dplyr::summarise(n = dplyr::n(), median_pct = stats::median(pct),
                   q25 = stats::quantile(pct, 0.25, names = FALSE),
                   q75 = stats::quantile(pct, 0.75, names = FALSE),
                   ci_lo = .boot_median(pct)[1], ci_hi = .boot_median(pct)[2],
                   .groups = "drop") %>%
  dplyr::arrange(panel, pathway, arm)

message("\n   ", PRIMARY, ", intersected panel - THE PRIMARY READOUT:")
ranks_by_group %>%
  dplyr::filter(panel == "intersected (142)", pathway == PRIMARY) %>%
  dplyr::transmute(arm, n, median_pct = round(median_pct, 1),
                   IQR = sprintf("[%.1f, %.1f]", q25, q75),
                   boot_ci = sprintf("[%.1f, %.1f]", ci_lo, ci_hi)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   ", SECONDARY, " umbrella, intersected panel - SECONDARY:")
ranks_by_group %>%
  dplyr::filter(panel == "intersected (142)", pathway == SECONDARY) %>%
  dplyr::transmute(arm, median_pct = round(median_pct, 1),
                   IQR = sprintf("[%.1f, %.1f]", q25, q75)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   ", PRIMARY, " on the SIZE-MATCHED panels",
        if (G2_STRUCTURAL) " - THIS IS WHERE R1 IS SCORED, because G2 fired:"
        else " - reported; G2 did not fire, so R1 is scored on the intersection:")
ranks_by_group %>%
  dplyr::filter(pathway == PRIMARY, grepl("^size-matched", panel)) %>%
  dplyr::transmute(panel, arm, median_pct = round(median_pct, 1),
                   IQR = sprintf("[%.1f, %.1f]", q25, q75)) %>%
  dplyr::arrange(panel, arm) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the apoptosis readout, apoptosis-split panel - DESCRIPTIVE, NO RULE:")
ranks_by_group %>%
  dplyr::filter(panel == "apoptosis-split (144)", pathway %in% APO_READ) %>%
  dplyr::transmute(pathway, arm, median_pct = round(median_pct, 1),
                   IQR = sprintf("[%.1f, %.1f]", q25, q75)) %>%
  dplyr::arrange(pathway, arm) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 7. PART D - THE MOUSE INTERNAL CONTRAST (G3)
# =============================================================================
# 6W WT against 12W WT on the OXPHOS percentile. This is the reference the human
# placement is read against and IT MUST EXIST INDEPENDENTLY before the human
# number means anything.
message("\n7. PART D: the mouse internal contrast, 6W WT vs 12W WT")

.contrast <- function(d, hi, lo, panel_lab, path_lab) {
  a <- d$pct[d$arm == hi]; b <- d$pct[d$arm == lo]
  dd <- vapply(seq_len(NBOOT), function(i)
    stats::median(sample(a, replace = TRUE)) - stats::median(sample(b, replace = TRUE)), 0)
  q <- stats::quantile(dd, c(0.025, 0.975), names = FALSE)
  tibble::tibble(panel = panel_lab, pathway = path_lab, hi = hi, lo = lo,
                 n_hi = length(a), n_lo = length(b),
                 median_hi = stats::median(a), median_lo = stats::median(b),
                 diff = stats::median(a) - stats::median(b),
                 ci_lo = q[1], ci_hi = q[2]) }

mouse_interval <- dplyr::bind_rows(lapply(
  c("intersected (142)", "own full panel", "apoptosis-split (144)", names(bands)), function(p)
    dplyr::bind_rows(lapply(SETS_MAIN, function(s) {
      d <- ranks_per_sample %>% dplyr::filter(panel == p, pathway == s)
      if (!nrow(d)) return(NULL)
      .contrast(d, "mouse 12W WT", "mouse 6W WT", p, s) }))))
mouse_interval %>%
  dplyr::transmute(panel, pathway, median_6W = round(median_lo, 1),
                   median_12W = round(median_hi, 1),
                   diff = sprintf("%+.1f", diff),
                   ci = sprintf("[%+.1f, %+.1f]", ci_lo, ci_hi)) %>%
  as.data.frame() %>% print(row.names = FALSE)

g3row <- mouse_interval %>%
  dplyr::filter(panel == "intersected (142)", pathway == PRIMARY)
# the same contrast on the panel R1 will actually be scored on, computed after
# PANEL_FOR_R1 exists; see section 9.
G3_FELL <- g3row$diff < 0
G3_EXCLUDES_ZERO <- g3row$ci_hi < 0
message("\n   G3: OXPHOS percentile ", if (G3_FELL) "FALLS" else "DOES NOT FALL",
        " from 6W WT to 12W WT (", sprintf("%+.1f", g3row$diff),
        sprintf(" [%+.1f, %+.1f])", g3row$ci_lo, g3row$ci_hi))
if (!G3_FELL) {
  message("   *** G3 FAILS. There is no reprioritisation IN RANK TERMS, so the",
          "\n   *** human placement has nothing to be placed against. The",
          "\n   *** placement below is NOT to be read; see the verdict.")
}

# =============================================================================
# 8. PART E - THE FIGURE
# =============================================================================
message("\n8. PART E: figure")
dir.create(DIR_E24_FIG, showWarnings = FALSE, recursive = TRUE)

fig_d <- ranks_per_sample %>%
  dplyr::filter(pathway == PRIMARY,
                panel %in% c("intersected (142)", "own full panel")) %>%
  dplyr::mutate(panel = factor(panel, levels = c("intersected (142)", "own full panel")))
fig_s <- ranks_by_group %>%
  dplyr::filter(pathway == PRIMARY,
                panel %in% c("intersected (142)", "own full panel")) %>%
  dplyr::mutate(panel = factor(panel, levels = c("intersected (142)", "own full panel")))

p <- ggplot2::ggplot(fig_d, ggplot2::aes(x = arm, y = pct)) +
  ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 1) +
  ggplot2::geom_linerange(data = fig_s,
                          ggplot2::aes(x = arm, ymin = q25, ymax = q75),
                          inherit.aes = FALSE, linewidth = 0.9) +
  ggplot2::geom_point(data = fig_s, ggplot2::aes(x = arm, y = median_pct),
                      inherit.aes = FALSE, size = 3) +
  ggplot2::facet_wrap(~ panel) +
  ggplot2::labs(
    title = paste0("Percentile rank of `", PRIMARY, "` within the mitochondrial panel"),
    subtitle = paste0("within-sample rank, ", PCT_NOTE,
                      ". Points are samples, bars the IQR. RANK ONLY - no ",
                      "mitoPPS value is compared across species."),
    x = NULL, y = "percentile within sample") +
  ggplot2::ylim(0, 100) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(DIR_E24_FIG, "E24_oxphos_percentile.pdf"), p,
                width = 11, height = 5)
message("   wrote ", file.path("outputs", "mitopps_rank", "E24_oxphos_percentile.pdf"))

# =============================================================================
# 9. PART F - VERDICT ON THE HEADER'S RULES
# =============================================================================
message("\n9. PART F: verdict")

# WHICH PANEL R1 IS SCORED ON is decided by G2's outcome, per the rule fixed in
# the header - not chosen after seeing the placement.
PANEL_FOR_R1 <- if (G2_STRUCTURAL) BAND_MAIN else "intersected (142)"
message("   R1 is scored on: ", PANEL_FOR_R1,
        if (G2_STRUCTURAL) "  (because G2 fired)" else "")
prim <- ranks_by_group %>%
  dplyr::filter(panel == PANEL_FOR_R1, pathway == PRIMARY)
gv <- function(a) prim$median_pct[prim$arm == a]
m6 <- gv("mouse 6W WT"); m12 <- gv("mouse 12W WT")
tc <- gv("TCGA"); sb <- gv("SCAN-B")

# R1 applied mechanically: which mouse end is each human cohort nearer, and is
# it INSIDE the 6W-12W interval at all? "Between or outside" is reported as such
# and never assigned to the nearer end.
.place <- function(h) {
  lo <- min(m6, m12); hi <- max(m6, m12)
  if (h < lo) return("OUTSIDE, below both mouse positions")
  if (h > hi) return("OUTSIDE, above both mouse positions")
  d6 <- abs(h - m6); d12 <- abs(h - m12)
  span <- abs(m12 - m6)
  if (span > 0 && min(d6, d12) / span > 0.25)
    return("BETWEEN the two mouse positions, not near either")
  if (d12 < d6) "near the 12W position" else "near the 6W position" }
R1_TCGA <- .place(tc); R1_SCANB <- .place(sb)
R2_AGREE <- identical(R1_TCGA, R1_SCANB)

verdict <- tibble::tibble(
  rule = c("G1 panel comparability", "G2 the size gate", "G3 the mouse interval",
           "R1 TCGA placement", "R1 SCAN-B placement", "R2 cohort agreement"),
  outcome = c(
    if (G1_PASSED) "PASSED - the panels are the same object" else "FAILED",
    if (G2_STRUCTURAL)
      paste0("FAILED - rho(rank, log size) >= 0.5 in both cohorts and the primary is ",
             "one of the largest sets in the panel. R1 is therefore scored on ",
             BAND_MAIN, ", the remedy the rule prescribes, NOT on the panel")
      else "passed - the rank is not substantially structural",
    if (!G3_FELL) "FAILED - no reprioritisation in rank terms; the placement is NOT to be read"
      else if (G3_EXCLUDES_ZERO) "HOLDS - OXPHOS falls 6W WT to 12W WT, interval excluding zero"
      else "FALLS but the interval covers zero - a weak reference",
    R1_TCGA, R1_SCANB,
    if (R2_AGREE) "the cohorts agree" else "THE COHORTS DISAGREE - that is the result"))
verdict %>% as.data.frame() %>% print(row.names = FALSE, right = FALSE)

# READABLE is the conjunction the header fixed: the reference must exist (G3)
# AND the rank must not be substantially structural (G2). Neither was relaxed
# after the numbers were seen.
READABLE <- G3_FELL && (!G2_STRUCTURAL || BAND_READABLE)
message("\n   IS THE PLACEMENT READABLE? ", if (READABLE) "yes" else "NO",
        if (!READABLE) " - see G2/G3 above; the note must say so before quoting R1." else "")
message("   AND WHATEVER IT SAYS, IT IS ONE DESCRIPTIVE OBSERVATION, NOT A TEST.")

NOTES <- c(
  paste0("Percentile: ", PCT_NOTE, "."),
  paste0("GATE 0: ", G1_NOTE),
  paste0("PART A2's control passed at max |delta| = ",
         signif(max(ctrl$max_abs_delta), 3),
         " - the human 142 universe rebuilds bit-equal in both cohorts, which is",
         " what licenses reading the 144 run."),
  paste0("The human apoptosis split is an independently curated ANALOGUE of the",
         " mouse's, not a translation: PRO agrees at ", length(apo_pro),
         ", ANTI does not (", length(apo_anti), " against ", msz[["Apoptosis-ANTI"]],
         "). Descriptive only; no reading rule is pre-specified for it."),
  "RANK ONLY. No mitoPPS value is saved in this object and none appears in the figure.",
  "Exploratory and post-hoc. Both R1 outcomes were pre-declared; neither was predicted.")

saveRDS(list(
  panel_reconciliation = panel_reconciliation,
  intersection = intersection, size_bands = bands,
  panel_for_r1 = PANEL_FOR_R1,
  a2_control = ctrl, a2_shift = a2_shift, apoptosis_split = apoptosis_split,
  size_gate = size_gate,
  ranks_per_sample = ranks_per_sample, ranks_by_group = ranks_by_group,
  mouse_interval = mouse_interval,
  verdict = verdict, readable = READABLE,
  g2_structural = G2_STRUCTURAL, g3_fell = G3_FELL,
  r1 = c(TCGA = R1_TCGA, `SCAN-B` = R1_SCANB), r2_agree = R2_AGREE,
  settings = list(
    question = paste("where does OXPHOS sit within the mitochondrial pathway",
                     "panel, in mouse gland and in human tumours - RANK, never value"),
    percentile = PCT_NOTE, primary = PRIMARY, secondary = SECONDARY,
    apoptosis_readout = APO_READ, cut_size = CUT_SIZE, nboot = NBOOT, seed = 1,
    mouse_artefact = list(path = "data/from_myc_mouse/mitopps_scores.rds",
                          md5 = MOUSE_MD5, recomputed = FALSE),
    arm_order = ARM_ORDER),
  rules = list(
    rank_only = paste("mitoPPS values are never compared across species; this",
                      "object carries percentile rank only"),
    no_pooling = "mouse and human are never scored together; TCGA and SCAN-B never averaged",
    exploratory = "post-hoc and descriptive; not a test of any hypothesis",
    n3 = "transcript associations; 'primed' is never written of a transcript",
    apoptosis = "reported and not scored - the two splits are different curations"),
  analysis_date = Sys.Date(), notes = NOTES), PATH_E24)

readr::write_csv(ranks_by_group %>% dplyr::mutate(arm = as.character(arm)), PATH_E24_TAB)

message("\nE24: done.")
message("    results/mitopps_rank_comparison.rds")
message("    outputs/tables/E24_ranks_by_group.csv")
message("    outputs/mitopps_rank/E24_oxphos_percentile.pdf")

# =============================================================================
# SANDBOX - run line by line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E24)

  ## 1. GATE 0 FIRST. If the panels are not the same object nothing else reads.
  x$panel_reconciliation |> as.data.frame()
  x$intersection$n_shared; x$intersection$mouse_only; x$intersection$human_only

  ## 2. PART A2's CONTROL. Must be ~0 in both cohorts, else A2 is unreadable.
  x$a2_control |> as.data.frame()
  ## and what adding two pathways cost the primary's percentile
  x$a2_shift |> as.data.frame()

  ## 3. G2, THE SIZE GATE - read BEFORE any placement. If structural, the
  ## placement is readable only against size-matched pathways.
  x$size_gate |> as.data.frame()

  ## 4. G3 - the mouse interval. If OXPHOS does not FALL 6W -> 12W, stop here.
  x$mouse_interval |> as.data.frame()

  ## 5. THE PRIMARY READOUT. Intersected panel, OXPHOS subunits, by group.
  ## Read the IQR, not just the median - the comparison has no error bar without it.
  subset(x$ranks_by_group,
         panel == "intersected (142)" & pathway == "OXPHOS subunits") |>
    as.data.frame()

  ## 6. The sensitivity: each species' own full panel. If it disagrees with the
  ## primary, panel composition is driving the comparison.
  subset(x$ranks_by_group,
         panel == "own full panel" & pathway == "OXPHOS subunits") |>
    as.data.frame()

  ## 7. The apoptosis readout. DESCRIPTIVE, NO RULE - the two splits are
  ## different curations (PRO agrees, ANTI does not). Read x$apoptosis_split first.
  x$apoptosis_split |> as.data.frame()
  subset(x$ranks_by_group, panel == "apoptosis-split (144)") |> as.data.frame()

  ## 8. The verdict, on the rules fixed in the header before any number.
  x$verdict |> as.data.frame()
  x$readable; x$r1; x$r2_agree
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")

}
