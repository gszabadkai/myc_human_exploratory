# E02_score_cohorts.R
# =============================================================================
# The one heavy script. Brings both cohorts to the same scoring, on the same
# set definitions, so their CORRELATION STRUCTURES can be compared.
#
# Built to: docs/2026-08-31_phase1_plan.md section 2 (E02), sections 1.3, 1.4
#
# =============================================================================
# TWO JOBS, BECAUSE THE SNAPSHOT IS ASYMMETRIC
# =============================================================================
#   SCAN-B needs the FULL mitochondrial scoring - it has almost none.
#   TCGA already carries 18 arms on 4 instruments, and needs only the NEW sets
#     (the death axis and the MYC compendium).
#
# THE TCGA MITOCHONDRIAL ARMS ARE REUSED, NOT RECOMPUTED. Their values are
# reported in the validation study and E01 asserts three of them. Re-scoring
# would risk a silently different answer for no gain.
#
# =============================================================================
# SCALE - opposite requirements, never the same object
# =============================================================================
#   GSVA     -> *_vst$mat      LOG (variance-stabilised), kcdf = "Gaussian"
#   mitoPPS  -> *_linear$mat   LINEAR (DESeq2-normalised)
#   content  -> LINEAR         log2(colSums + 1)
#   zmean    -> LOG            mean of per-gene z across samples
#   log2(MYC)-> LINEAR
#
# WHY EVERYTHING GOES IN ONE GSVA CALL PER COHORT. GSVA ranks each sample across
# the gene universe the call is given. Two calls carrying different set
# collections walk DIFFERENT universes and their scores are not comparable. The
# .PIN_A/.PIN_B half-matrix sets pin the universe to the full matrix in every
# call, which is what makes a second call (TCGA's new sets, added beside arms
# scored in 2026-08) legitimate. Two halves rather than one all-genes set,
# because a set containing every gene leaves the walk's miss-penalty at 0/0.
#
# COHORT-RELATIVITY. Scores are comparable WITHIN a cohort and never between.
# E03 compares correlations, never values. See CLAUDE.md.
#
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "mitopps.R"))
suppressPackageStartupMessages(library(GSVA))

message("\nE02: score both cohorts on one set of definitions\n", strrep("=", 78))

# =============================================================================
# 0. Constants
# =============================================================================
PATH_SCANB_SCORES <- file.path(DIR_RESULTS, "scanb_scores.rds")
PATH_NEW_SETS     <- file.path(DIR_RESULTS, "new_set_scores.rds")
PATH_SETDEFS      <- file.path(DIR_RESULTS, "set_definitions.rds")

MTDNA_PATHWAY <- "mtDNA-encoded OXPHOS subunits"

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

tcga_vst  <- readRDS(PATH_TCGA_VST)
tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
tcga_mito <- readRDS(PATH_TCGA_MITO)
tcga_myc  <- readRDS(PATH_TCGA_MYC)$estimators
scanb_vst <- readRDS(PATH_SCANB_VST)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
scanb_ph  <- readRDS(PATH_SCANB_PHENO)
g1        <- readRDS(PATH_G1)

stopifnot(identical(tcga_vst$scale, "log_vst"),
          identical(tcga_lin$scale, "linear_deseq2_normalised"),
          identical(scanb_vst$scale, "log_vst"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))

felsher_61      <- g1$estimators_stripped$FELSHER
collectri_full  <- g1$estimators_raw$COLLECTRI_MYC_ALL
collectri_strip <- g1$estimators_stripped$COLLECTRI_MYC_ALL
stopifnot(length(felsher_61) == EXPECT_FELSHER_STRIP,
          length(collectri_full) == EXPECT_COLLECTRI_FULL,
          length(collectri_strip) == EXPECT_COLLECTRI_STRIP)
message("   Felsher: ", length(g1$estimators_raw$FELSHER), " raw / ",
        length(felsher_61), " mito-stripped | CollecTRI: ",
        length(collectri_full), " raw / ", length(collectri_strip),
        " mito-stripped")

mitocarta_inventory  <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
mitocarta_background <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 3))
stopifnot(all(c("Symbol", "Synonyms") %in% colnames(mitocarta_background)))

# =============================================================================
# 2. The set definitions
# =============================================================================
# Built ONCE, in human symbols, and used for BOTH cohorts. Each cohort then
# harmonises them to its own vocabulary (section 3) - the sets themselves never
# change.
message("\n2. set definitions")

# --- 2.1 mitochondrial arms and pathway universe, reused from the snapshot ---
# NOT rebuilt. These are the definitions the validation study's 18 arms were
# scored on, and reusing them is what makes TCGA's existing scores usable
# alongside SCAN-B's new ones.
arm_sets     <- tcga_mito$arm_sets
cov_sets     <- tcga_mito$covariate_sets
mito_paths   <- tcga_mito$mito_paths
arm_univ_map <- tcga_mito$arm_universe_path
message("   mitochondrial arms: ", length(arm_sets),
        " | MitoPathway universe: ", length(mito_paths))

# --- 2.2 the MYC estimator panel, in FOUR EXPLICITLY LABELLED VARIANTS -------
# See the naming contract in E00. Phase 1 scored 18 sets under bare names, one
# of which (FELSHER_61) had already been stripped of every MitoCarta gene while
# the other 17 had not, and F1 was written on that mixture. From here every
# scored estimator carries a suffix and there are no bare names.
#
# THE TWO REFERENCE SETS, TAKEN FROM THE VALIDATION STUDY'S OWN AUDIT so the
# strip is reproducible rather than re-derived:
#   MITOCARTA_ALL          1,136 genes - the set it actually stripped against
#   HALLMARK E2F + G2M       327 genes - what `frac_prolif` is measured against
#
# GMX layout: row 1 set names, row 2 a description line (all NA here), rows 3+
# genes, one column per set, NA-padded.
message("\n2.2 MYC estimator panel, four variants each")

gmx <- utils::read.delim(PATH_MYC_GMX, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE, na.strings = c("NA", ""))
gmx <- gmx[-1, , drop = FALSE]
myc_base <- lapply(gmx, function(v) sort(unique(v[!is.na(v) & v != ""])))
if (length(myc_base) != EXPECT_MYC_SETS) {
  stop("the MYC GMX carries ", length(myc_base), " signatures, expected ",
       EXPECT_MYC_SETS, ". Re-snapshot and update E00.", call. = FALSE)
}

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
.hm <- function(nm) unique(hallmark$gene_symbol[hallmark$gs_name == nm])
myc_base[["HALLMARK_MYC_TARGETS_V1"]] <- .hm("HALLMARK_MYC_TARGETS_V1")

# FELSHER enters as its RAW 67 genes. The 61-gene version everyone has been
# calling FELSHER_61 is the MITOSTRIP variant and is reconstructed below, not
# imported, so the reconstruction can be checked against the original.
felsher_full <- g1$estimators_raw$FELSHER
stopifnot(length(felsher_full) == EXPECT_FELSHER_FULL)
myc_base[["FELSHER"]] <- sort(unique(felsher_full))

MITOCARTA_ALL <- g1$reference_sets$MITOCARTA_ALL
PROLIF_REF    <- unique(c(.hm("HALLMARK_E2F_TARGETS"),
                          .hm("HALLMARK_G2M_CHECKPOINT")))
if (length(MITOCARTA_ALL) != EXPECT_MITOCARTA_ALL ||
    length(PROLIF_REF) != EXPECT_PROLIF_REF) {
  stop("a strip reference set changed size: MITOCARTA_ALL ",
       length(MITOCARTA_ALL), " (expected ", EXPECT_MITOCARTA_ALL,
       "), PROLIF_REF ", length(PROLIF_REF), " (expected ", EXPECT_PROLIF_REF,
       "). Re-snapshot and update E00 rather than adjusting the expectation.",
       call. = FALSE)
}
message("   strip references: MITOCARTA_ALL ", length(MITOCARTA_ALL),
        " | HALLMARK E2F+G2M ", length(PROLIF_REF))

.variants <- function(g) list(
  `__FULL`        = g,
  `__MITOSTRIP`   = setdiff(g, MITOCARTA_ALL),
  `__PROLIFSTRIP` = setdiff(g, PROLIF_REF),
  `__BOTHSTRIP`   = setdiff(g, union(MITOCARTA_ALL, PROLIF_REF)))
myc_sets <- list()
for (b in names(myc_base)) {
  v <- .variants(myc_base[[b]])
  for (sfx in names(v)) myc_sets[[paste0(b, sfx)]] <- sort(v[[sfx]])
}

# THE RECONSTRUCTION CHECK. If setdiff(FELSHER-67, MITOCARTA_ALL) is not the
# validation study's 61-gene M_a set, then this script's understanding of what
# was stripped is wrong and every label below is wrong with it.
if (!identical(myc_sets[["FELSHER__MITOSTRIP"]], sort(felsher_61))) {
  stop("FELSHER__MITOSTRIP (", length(myc_sets[["FELSHER__MITOSTRIP"]]),
       " genes) does not reproduce the validation study's stripped FELSHER (",
       length(felsher_61), "). The strip reference or the raw set has moved; ",
       "do not trust any variant label until this is resolved.", call. = FALSE)
}
message("   FELSHER__MITOSTRIP reproduces the validation study's M_a set ",
        "exactly (", length(felsher_61), " genes)")

# Proliferation entanglement and mitochondrial content, both carried with the
# scores so no downstream script has to re-derive either (CLAUDE.md trap 3, and
# D0's confound applied to the MYC panel).
myc_panel <- tibble::tibble(signature = names(myc_sets)) %>%
  dplyr::mutate(
    base         = sub("__[A-Z]+$", "", signature),
    strip_status = sub("^.*__", "", signature),
    n            = vapply(myc_sets[signature], length, integer(1)),
    n_full       = vapply(myc_base[base], length, integer(1)),
    n_removed    = n_full - n,
    n_prolif     = vapply(myc_sets[signature],
                          function(g) sum(g %in% PROLIF_REF), integer(1)),
    n_mito       = vapply(myc_sets[signature],
                          function(g) sum(g %in% MITOCARTA_ALL), integer(1)),
    frac_prolif  = n_prolif / n,
    frac_mito    = n_mito / n,
    thin         = n < MIN_GSVA_N) %>%
  dplyr::arrange(base, match(strip_status,
                             c("FULL", "MITOSTRIP", "PROLIFSTRIP", "BOTHSTRIP")))

# The strips must actually have worked.
stopifnot(all(myc_panel$n_mito[myc_panel$strip_status %in%
                                 c("MITOSTRIP", "BOTHSTRIP")] == 0L),
          all(myc_panel$n_prolif[myc_panel$strip_status %in%
                                   c("PROLIFSTRIP", "BOTHSTRIP")] == 0L))

myc_panel %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1),
                frac_mito = round(100 * frac_mito, 1)) %>%
  dplyr::select(signature, strip_status, n, n_removed, frac_prolif, frac_mito,
                thin) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   ", nrow(myc_panel), " estimators = ", length(myc_base),
        " signatures x 4 variants")
if (any(myc_panel$thin)) {
  message("   THIN (n < ", MIN_GSVA_N, "), scored but not to be read alone: ",
          paste(myc_panel$signature[myc_panel$thin], collapse = ", "))
}
fp <- myc_panel$frac_prolif[myc_panel$strip_status == "FULL"]
message("   entanglement across the FULL variants spans ",
        round(100 * min(fp), 1), "% to ", round(100 * max(fp), 1),
        "%. NEVER report a MYC-OXPHOS correlation from one signature, and",
        " never\n   compare a FULL to a stripped one without saying so.")

# --- 2.2b the CollecTRI regulon, same four variants --------------------------
# M_b is the odd member of the panel twice over: it is a signed regulon scored
# by ULM rather than a gene set scored by GSVA, AND the version phase 1 used was
# the mito-stripped one. Both facts were invisible in the name. It now gets the
# same four labels as everything else so "is M_b weak because it is stripped?"
# is a question the atlas can answer instead of a question about the atlas.
collectri_sets <- list()
for (sfx in names(.variants(collectri_full)))
  collectri_sets[[paste0("M_b", sfx)]] <- sort(.variants(collectri_full)[[sfx]])
if (!identical(collectri_sets[["M_b__MITOSTRIP"]], sort(collectri_strip))) {
  stop("M_b__MITOSTRIP does not reproduce the validation study's stripped ",
       "CollecTRI regulon. Do not trust the variant labels.", call. = FALSE)
}
message("   CollecTRI variants: ",
        paste(sprintf("%s=%d", names(collectri_sets),
                      vapply(collectri_sets, length, integer(1))),
              collapse = " | "))

# --- 2.3 the death axis ------------------------------------------------------
# Two sources, both human-native. See data/genesets_celldeath_human/README.md
# for why taking the human column is not an ortholog round trip, and why no
# remapping is applied.
message("\n2.3 death axis")

cdc <- readr::read_csv(PATH_CDC, show_col_types = FALSE, progress = FALSE)
if (nrow(cdc) != EXPECT_CDC_ROWS) {
  stop("cell_death_genes_consolidated.csv has ", nrow(cdc), " rows, expected ",
       EXPECT_CDC_ROWS, ". The upstream file moved; re-snapshot rather than ",
       "patching, and update E00.", call. = FALSE)
}
stopifnot(!anyNA(cdc$human_symbol), !any(cdc$human_symbol == ""),
          dplyr::n_distinct(cdc$human_symbol) == nrow(cdc))
if (!all(cdc$effect %in% c("pro-death", "pro-survival", "ambiguous", "unclassified")) ||
    !all(cdc$pathway %in% c("apoptosis", "CICD", "both"))) {
  stop("unexpected value in cdc$effect or cdc$pathway. Inspect before using.",
       call. = FALSE)
}
message("   consolidated table: ", nrow(cdc), " genes, all with a human symbol")

# The 8 CDC strata. `pathway == "both"` contributes to BOTH, as the mouse arm
# did. `ambiguous` and `unclassified` are excluded from scores and reported
# separately - they are not a third category, they are missing information.
.cdc_stratum <- function(eff, pw, mito) {
  s <- cdc[cdc$effect == eff & cdc$pathway %in% c(pw, "both"), ]
  if (mito) s <- s[s$is_mitochondrial %in% c(TRUE, "TRUE"), ]
  sort(unique(s$human_symbol))
}
cdc_sets <- list()
for (eff in c("pro-death", "pro-survival")) {
  for (pw in c("apoptosis", "CICD")) {
    for (mt in c(FALSE, TRUE)) {
      nm <- sprintf("CDC_%s_%s%s",
                    ifelse(eff == "pro-death", "PRODEATH", "PROSURVIVAL"),
                    toupper(pw), ifelse(mt, "_MITO", ""))
      cdc_sets[[nm]] <- .cdc_stratum(eff, pw, mt)
    }
  }
}
cdc_tab <- tibble::tibble(set = names(cdc_sets),
                          n = vapply(cdc_sets, length, integer(1))) %>%
  dplyr::mutate(scored = n >= MIN_SCORE_N)
cdc_tab %>% as.data.frame() %>% print(row.names = FALSE)

# CLAUDE.md trap 9. CICD is the axis of most interest and the weakest measured.
below <- cdc_tab$set[!cdc_tab$scored]
if (length(below)) {
  message("   NOT SCORED (n < ", MIN_SCORE_N, "), carried as gene lists: ",
          paste(below, collapse = ", "))
  message("   This is deliberate. A 4-gene GSVA score is not a programme.")
}

# The 15 Tang modalities.
tang_files <- list.files(DIR_TANG, pattern = "\\.csv$", full.names = TRUE)
tang_sets  <- stats::setNames(
  lapply(tang_files, function(f) {
    d <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
    sort(unique(d$gene[!is.na(d$gene) & d$gene != ""]))
  }),
  paste0("TANG_", toupper(tools::file_path_sans_ext(basename(tang_files)))))

tang_n <- vapply(tang_sets, length, integer(1))
names(tang_n) <- tools::file_path_sans_ext(basename(tang_files))
if (!identical(sort(names(tang_n)), sort(names(EXPECT_TANG))) ||
    any(tang_n[names(EXPECT_TANG)] != EXPECT_TANG)) {
  stop("the Tang modality sizes do not match E00's expectations.\n",
       paste(sprintf("  %-32s got %5d, expected %5d", names(EXPECT_TANG),
                     tang_n[names(EXPECT_TANG)], EXPECT_TANG), collapse = "\n"),
       "\nRe-snapshot rather than patching, and update E00.", call. = FALSE)
}
tang_tab <- tibble::tibble(set = names(tang_sets), n = unname(tang_n)) %>%
  dplyr::mutate(scored = n >= MIN_GSVA_N,
                huge   = n >= 500L)   # ~3% of the matrix and up
tang_tab %>% dplyr::arrange(dplyr::desc(n)) %>% as.data.frame() %>%
  print(row.names = FALSE)
message("   not scored (n < ", MIN_GSVA_N, "): ",
        paste(tang_tab$set[!tang_tab$scored], collapse = ", "))
message("   NEAR-TRANSCRIPTOME-WIDE (n >= 500, CLAUDE.md trap 10): ",
        paste(tang_tab$set[tang_tab$huge], collapse = ", "))

# Gene-level overlay labels. Too small to score; carried as annotation and used
# by E05 as points on the MYC-OXPHOS plane.
family_labels <- cdc %>%
  dplyr::filter(!is.na(family_pathway), family_pathway != "NA") %>%
  dplyr::select(gene = human_symbol, family_pathway, effect, pathway,
                is_mitochondrial)
BCL2_FAMILY <- c("BCL2", "BCL2L1", "MCL1", "BBC3", "BCL2L11", "BAX", "BAK1",
                 "BID", "BCL2L2", "BCL2A1", "BAD", "BIK", "BMF", "HRK", "PMAIP1")
message("   gene-level overlay: ", nrow(family_labels), " genes with a ",
        "family_pathway label; ", sum(BCL2_FAMILY %in% cdc$human_symbol),
        " of ", length(BCL2_FAMILY), " BCL2-family genes in the table")

death_sets <- c(cdc_sets[cdc_tab$scored], tang_sets[tang_tab$scored])
message("   -> ", length(death_sets), " death sets will be SCORED (",
        sum(cdc_tab$scored), " CDC + ", sum(tang_tab$scored), " Tang)")

# The full collection every cohort gets, in ONE call.
NEW_SETS <- c(myc_sets, death_sets)
message("\n   new sets to score in both cohorts: ", length(NEW_SETS),
        " (", length(myc_sets), " MYC + ", length(death_sets), " death)")

# =============================================================================
# 3. Scoring machinery
# =============================================================================
message("\n3. machinery")

# --- 3.1 symbol harmonisation, per cohort ------------------------------------
# Script 07 section 2's map, applied to whichever matrix is being scored. It
# matters for SCAN-B (a 2014 UCSC build: 19 of 89 OXPHOS subunit genes are
# pre-2018 ATP-synthase names) and is a no-op for TCGA.
#
# Forward only: current MitoCarta symbol -> its listed alias. The reverse is
# dangerous - COX1/COX2/COX3 resolve to the PROSTAGLANDIN SYNTHASES. A candidate
# that is itself a current symbol for a different gene is refused, and so is any
# symbol with more than one surviving candidate.
#
# ONE DEPRECATED ROW BELONGS TO EXACTLY ONE GENE. MitoCarta's Synonyms column
# does not respect this: it lists QARS under BOTH `QARS1` (the true rename) and
# `EPRS1` (which carries the whole multi-synthetase complex's abbreviations,
# EARS|PARS|QARS|QPRS), and GARS under both `GARS1` and `GART`. So a candidate
# is admitted only when it is claimed by exactly ONE symbol of the gene universe
# being mapped. QARS and GARS are contested and go to NEITHER claimant - the
# conservative direction, costing QARS1 and GARS1 their SCAN-B rows.
#
# This replaces an earlier `setdiff(cand, genes)`, which refused any candidate
# appearing in the universe at all. That punished a gene for its own old name
# being present: HALLMARK_MYC_TARGETS_V1 spells it EPRS1 and MENSSEN_MYC_TARGETS
# spells it EPRS, and the old rule therefore denied EPRS1 the row EPRS and left
# it the wrong row QARS. It also silently cost H2AZ1, POLR1G and VARS1 their
# resolutions, which is why the agreement check below is load-bearing.
#
# This is NOT an ortholog mapping. It is a within-species vocabulary
# reconciliation. No ortholog function is called (CLAUDE.md trap 11).
MC_SYMBOLS <- unique(mitocarta_background$Symbol)
.syn_map <- local({
  syn <- strsplit(ifelse(is.na(mitocarta_background$Synonyms), "",
                         mitocarta_background$Synonyms), "|", fixed = TRUE)
  tibble::tibble(symbol = rep(mitocarta_background$Symbol, lengths(syn)),
                 alias  = unlist(syn)) %>%
    dplyr::filter(!is.na(alias), alias != "", alias != symbol) %>%
    dplyr::distinct()
})
# alias -> every current symbol that lists it. The contest test.
.alias_owners <- split(.syn_map$symbol, .syn_map$alias)

.build_symbol_map <- function(genes, matrix_symbols) {
  genes   <- sort(unique(genes[!is.na(genes) & genes != ""]))
  present <- genes %in% matrix_symbols
  out     <- stats::setNames(genes, genes)
  status  <- ifelse(present, "matched", "unresolved")
  for (i in which(!present)) {
    g    <- genes[[i]]
    cand <- .syn_map$alias[.syn_map$symbol == g]
    cand <- cand[cand %in% matrix_symbols]
    cand <- cand[!(cand %in% MC_SYMBOLS & cand != g)]
    n_reachable <- length(cand)
    cand <- cand[vapply(cand, function(cc)
      length(intersect(.alias_owners[[cc]], genes)) == 1L, logical(1))]
    if (length(cand) == 1L)     { out[[g]] <- cand; status[[i]] <- "resolved" }
    else if (length(cand) > 1L)   status[[i]] <- "ambiguous"
    else if (n_reachable > 0L)    status[[i]] <- "contested"
  }
  list(map = out, report = tibble::tibble(input_symbol = genes, status = status,
                                          resolved_to = unname(out[genes])))
}

# Reports the map, and stops if it ever puts two genes on one row.
#
# A row reached by BOTH a synonym and its own current name is one gene spelled
# two ways (EPRS1 -> EPRS, where EPRS is also an input) - legitimate, named, and
# harmless because every set is intersected against the matrix. A row reached by
# two DIFFERENT synonym-resolvers is two genes collapsed into one, and the
# candidate filter above is what makes it unreachable. The guard stays anyway.
.check_symbol_map <- function(bm, matrix_symbols, label) {
  print(table(bm$report$status))
  contested <- bm$report$input_symbol[bm$report$status == "contested"]
  if (length(contested)) {
    message("   contested (an alias in the matrix, but claimed by another gene ",
            "of this universe; deliberately unresolved): ",
            paste(contested, collapse = ", "))
  }
  m     <- bm$map[unname(bm$map) %in% matrix_symbols]
  claim <- split(names(m), unname(m))
  claim <- claim[lengths(claim) > 1L]
  if (!length(claim)) return(invisible(NULL))
  n_res <- vapply(names(claim), function(t) length(setdiff(claim[[t]], t)),
                  integer(1))
  bad <- names(claim)[n_res > 1L]
  if (length(bad)) {
    stop("the ", label, " symbol map sends two DIFFERENT genes to one row: ",
         paste(sprintf("%s <- %s", bad,
                       vapply(bad, function(t) paste(claim[[t]], collapse = " + "),
                              character(1))), collapse = "; "), call. = FALSE)
  }
  message("   ", length(claim), " row(s) reached by a current symbol and its own ",
          "old name (one gene, two spellings): ",
          paste(sprintf("%s <- %s", names(claim),
                        vapply(names(claim),
                               function(t) paste(setdiff(claim[[t]], t), collapse = ","),
                               character(1))), collapse = ", "))
  invisible(NULL)
}

# --- 3.1b CollecTRI ULM, the M_b machinery -----------------------------------
# SIGN RULE, identical to myc_human_validation script 06 at d3ac60e:
# mor = +1 for a stimulatory edge, -1 for an inhibitory one, and an edge flagged
# BOTH takes +1. Reproduced rather than improved on, because M_b in the snapshot
# was built with it and the variants must be comparable to it.
ct <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE, progress = FALSE)
myc_net <- ct %>%
  dplyr::filter(source_genesymbol == "MYC", !is.na(target_genesymbol),
                target_genesymbol != "") %>%
  dplyr::transmute(source = "MYC", target = target_genesymbol,
                   mor = dplyr::if_else(as.logical(is_stimulation), 1, -1),
                   likelihood = 1) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)

.ulm <- function(targets, E, label) {
  net <- myc_net %>% dplyr::filter(target %in% targets, target %in% rownames(E))
  message(sprintf("   %-24s %3d targets in the matrix", label, nrow(net)))
  res <- decoupleR::run_ulm(mat = E, network = net, .source = source,
                            .target = target, .mor = mor, minsize = 5L)
  v <- res %>% dplyr::filter(statistic == "ulm", source == "MYC") %>%
    dplyr::select(condition, score)
  stats::setNames(v$score, v$condition)[colnames(E)]
}
# All four variants at once, as a matrix with the variant names as rownames.
.ulm_variants <- function(E, remap, label) {
  m <- t(vapply(names(collectri_sets), function(nm)
    .ulm(remap(collectri_sets[[nm]]), E, paste(label, nm)), numeric(ncol(E))))
  colnames(m) <- colnames(E); m
}

# --- 3.2 GSVA, with the universe pinned --------------------------------------
.gsva_batch <- function(sets, E, label) {
  n_ok <- vapply(sets, length, integer(1))
  if (any(n_ok < MIN_SET_GENES)) {
    stop("GSVA batch '", label, "': set(s) below the size floor -> ",
         paste(names(sets)[n_ok < MIN_SET_GENES], collapse = ", "), call. = FALSE)
  }
  wanted <- names(sets)
  syms <- rownames(E)
  sets[[".PIN_A"]] <- syms[c(TRUE, FALSE)]
  sets[[".PIN_B"]] <- syms[c(FALSE, TRUE)]
  par <- GSVA::gsvaParam(exprData = E, geneSets = sets, kcdf = "Gaussian",
                         minSize = MIN_SET_GENES, maxSize = Inf)
  s <- GSVA::gsva(par, verbose = FALSE)
  dropped <- setdiff(wanted, rownames(s))
  if (length(dropped)) {
    stop("GSVA silently dropped set(s) in '", label, "': ",
         paste(utils::head(dropped, 10), collapse = ", "), call. = FALSE)
  }
  s[wanted, , drop = FALSE]
}

# --- 3.3 mitoPPS, script 07 section 5 ----------------------------------------
# Reports the SHAPE of the mitochondrial programme, blind to total content.
# .mitopps_universe, .mitopps_query and .path_scores now live in
# functions/mitopps.R, sourced at the top, because E03b tests one of the scores
# they produce and the two must not drift apart.

# =============================================================================
# 4. SCAN-B: the full scoring
# =============================================================================
message("\n", strrep("=", 78), "\n4. SCAN-B\n", strrep("=", 78))

Es <- scanb_vst$mat
Ls <- scanb_lin$mat
SYM_S <- rownames(Es)

# --- 4.1 harmonise ------------------------------------------------------------
message("\n4.1 symbol harmonisation")
bm_s <- .build_symbol_map(
  c(unlist(arm_sets, use.names = FALSE), unlist(cov_sets, use.names = FALSE),
    unlist(mito_paths, use.names = FALSE), unlist(NEW_SETS, use.names = FALSE)),
  SYM_S)
.check_symbol_map(bm_s, SYM_S, "SCAN-B")
map_s <- bm_s$map
# Agreement with the map script 16 built in the validation study, over the
# genes both cover. A wider gene universe changing a resolution would mean the
# map is context-dependent and must be understood, not absorbed.
common <- intersect(names(map_s), names(scanb_ph$symbol_map))
dis <- common[map_s[common] != scanb_ph$symbol_map[common]]
if (length(dis)) {
  stop("the SCAN-B symbol map disagrees with the validation study's on ",
       length(dis), " symbol(s): ", paste(utils::head(dis, 10), collapse = ", "),
       call. = FALSE)
}
message("   agrees with the validation study's map on all ", length(common),
        " shared symbols")

.remap_s <- function(g) { h <- map_s[g]; unname(ifelse(is.na(h), g, h)) }
.in_s    <- function(g) intersect(.remap_s(unique(g)), SYM_S)

arm_s   <- lapply(arm_sets,   .in_s)
cov_s   <- lapply(cov_sets,   .in_s)
paths_s <- lapply(mito_paths, .in_s)
new_s   <- lapply(NEW_SETS,   .in_s)

ox_frac <- length(arm_s[["OXPHOS subunits"]]) / length(unique(arm_sets[["OXPHOS subunits"]]))
message(sprintf("   OXPHOS subunits coverage after harmonisation: %.3f", ox_frac))
if (ox_frac < 0.98) {
  stop("OXPHOS subunits covers only ", round(ox_frac, 3), " of its genes in ",
       "SCAN-B. Unharmonised it is 0.775; this suggests the map did not apply. ",
       "Do not score - 70 of 89 genes is a Complex V with no F1 head.",
       call. = FALSE)
}

# --- 4.2 GSVA, one call ------------------------------------------------------
message("\n4.2 GSVA (log VST, kcdf Gaussian) - one call, universe pinned")
all_s <- c(arm_s, cov_s, new_s)
message("   ", length(all_s), " sets x ", ncol(Es), " samples; a few minutes")
gsva_s <- .gsva_batch(all_s, Es, "SCAN-B")
message("   scored ", nrow(gsva_s), " sets")

# --- 4.3 mitoPPS -------------------------------------------------------------
message("\n4.3 mitoPPS (linear DESeq2-normalised)")
S_univ_s <- .path_scores(paths_s, Ls, MIN_SET_GENES)
if (min(S_univ_s) <= 0) {
  stop("a MitoPathway score is zero or negative in some SCAN-B sample; the ",
       "pairwise ratio is undefined.", call. = FALSE)
}
message("   universe: ", nrow(S_univ_s), " MitoPathways of ", length(paths_s))
mpps_univ_s <- .mitopps_universe(S_univ_s)
message(sprintf("   global mean %.4f (should be ~1)", mean(mpps_univ_s)))

# Each arm is queried against the universe with its OWN pathway held out, so an
# arm that IS a universe pathway reproduces its canonical value exactly. That
# identity is asserted - it is the check that the query form is right.
S_arms_s <- .path_scores(arm_s, Ls, MIN_SET_GENES)
mitopps_s <- t(vapply(rownames(S_arms_s), function(a) {
  hold <- arm_univ_map[[a]]
  U <- if (is.null(hold) || is.na(hold)) S_univ_s
       else S_univ_s[setdiff(rownames(S_univ_s), hold), , drop = FALSE]
  as.numeric(.mitopps_query(S_arms_s[a, , drop = FALSE], U))
}, numeric(ncol(Ls))))
colnames(mitopps_s) <- colnames(Ls)

n_canon <- 0L
for (a in rownames(mitopps_s)) {
  hold <- arm_univ_map[[a]]
  if (is.null(hold) || is.na(hold) || !hold %in% rownames(mpps_univ_s)) next
  d <- max(abs(mitopps_s[a, ] - mpps_univ_s[hold, ]))
  if (d > 1e-8) {
    stop("arm '", a, "' does not reproduce its canonical mitoPPS (diff ",
         signif(d, 3), ")", call. = FALSE)
  }
  n_canon <- n_canon + 1L
}
message("   ", n_canon, " arms reproduce their canonical value exactly")

# --- 4.4 content and zmean ---------------------------------------------------
# The two descriptive instruments. `content` is abundance-weighted by
# construction (a sum of linear counts); `zmean` is not. Together with GSVA
# (level) and mitoPPS (composition) they are the four rulers.
message("\n4.4 content (linear) and zmean (log)")
content_s <- t(vapply(arm_s, function(g)
  log2(colSums(Ls[g, , drop = FALSE]) + 1), numeric(ncol(Ls))))
colnames(content_s) <- colnames(Ls)

.zmean <- function(g, E) {
  sub <- E[g, , drop = FALSE]
  v <- apply(sub, 1L, stats::var)
  sub <- sub[v > 0, , drop = FALSE]      # a zero-variance row makes scale() NaN
  colMeans(t(scale(t(sub))))
}
zmean_s <- t(vapply(arm_s, .zmean, numeric(ncol(Es)), E = Es))
colnames(zmean_s) <- colnames(Es)

# --- 4.5 MYC, the estimators that are not gene-set scores --------------------
message("\n4.5 MYC: CollecTRI ULM (four variants) and raw expression")
# NOTE: the CollecTRI targets are NOT in the symbol map's input union, so
# .remap_s leaves any unmatched one alone rather than resolving it. That is
# deliberate, and the reason survives the 2026-08-31 rewrite of the candidate
# filter even though the mechanism changed. A candidate is now admitted only
# when exactly ONE symbol of the universe being mapped claims it, so adding 811
# CollecTRI targets could make an alias newly CONTESTED and withdraw a
# resolution the map already makes - which would break the agreement assertion
# in 4.1 for a secondary estimator's sake. M_b is reported alongside the panel,
# never alone. See docs/2026-08-31_symbol_map_note.md.
M_b_s <- .ulm_variants(Es, .remap_s, "SCAN-B")
if (!"MYC" %in% rownames(Ls)) stop("MYC absent from SCAN-B.", call. = FALSE)
log2MYC_s <- as.numeric(log2(Ls["MYC", ]))
names(log2MYC_s) <- colnames(Ls)

message(sprintf("   rho(%s, %s) in SCAN-B = %+.3f", MYC_REF, MB_REF,
                .rho(gsva_s[MYC_REF, ], M_b_s[MB_REF, ])))
message(sprintf("   rho(%s, log2(MYC))          = %+.3f", MYC_REF,
                .rho(gsva_s[MYC_REF, ], log2MYC_s)))

# =============================================================================
# 5. TCGA: only the new sets
# =============================================================================
message("\n", strrep("=", 78), "\n5. TCGA: the new sets only\n", strrep("=", 78))

Et <- tcga_vst$mat
Lt <- tcga_lin$mat
SYM_T <- rownames(Et)

bm_t <- .build_symbol_map(unlist(NEW_SETS, use.names = FALSE), SYM_T)
.check_symbol_map(bm_t, SYM_T, "TCGA")
map_t <- bm_t$map
.in_t <- function(g) intersect(unname(ifelse(is.na(map_t[unique(g)]),
                                             unique(g), map_t[unique(g)])), SYM_T)
new_t <- lapply(NEW_SETS, .in_t)

message("\n   GSVA on ", length(new_t), " new sets x ", ncol(Et), " samples")
message("   The 18 mitochondrial arms are NOT rescored - they are reused from")
message("   the snapshot, and the pins make these scores comparable to them.")
gsva_t_new <- .gsva_batch(new_t, Et, "TCGA new sets")

log2MYC_t <- as.numeric(log2(Lt["MYC", ]))
names(log2MYC_t) <- colnames(Lt)

# --- 5.1 the M_b variants for TCGA, and the check that they are the same
#         estimator the snapshot carries ---------------------------------------
# The snapshot holds ONE M_b, built by myc_human_validation script 06 from the
# mito-stripped regulon. Recomputing it here alongside the other three variants
# is what makes the four comparable; reproducing the snapshot's value is what
# proves the reconstruction is the same estimator and not merely a similar one.
message("\n5.1 M_b, four variants, TCGA")
M_b_t <- .ulm_variants(Et, function(g) g, "TCGA")
mb_snap <- tcga_myc$M_b[match(colnames(Et), tcga_myc$patient)]
mb_rho  <- .rho(M_b_t[MB_REF, ], mb_snap)
mb_diff <- max(abs(M_b_t[MB_REF, ] - mb_snap))
message(sprintf("   %s vs the snapshot's M_b: rho %.6f, max abs diff %.3g",
                MB_REF, mb_rho, mb_diff))
if (!is.finite(mb_rho) || mb_rho < 0.999) {
  stop("the recomputed ", MB_REF, " does not reproduce the snapshot's M_b ",
       "(rho ", signif(mb_rho, 4), "). The regulon, the sign rule or the ",
       "matrix has moved; the variant labels cannot be trusted.", call. = FALSE)
}

# =============================================================================
# 6. Coverage, both cohorts
# =============================================================================
message("\n6. coverage of the new sets in both cohorts")
coverage <- dplyr::bind_rows(
  tibble::tibble(cohort = "TCGA", set = names(NEW_SETS),
                 n_defined = vapply(NEW_SETS, function(g) length(unique(g)), integer(1)),
                 n_present = vapply(new_t, length, integer(1))),
  tibble::tibble(cohort = "SCAN-B", set = names(NEW_SETS),
                 n_defined = vapply(NEW_SETS, function(g) length(unique(g)), integer(1)),
                 n_present = vapply(new_s, length, integer(1)))) %>%
  dplyr::mutate(frac = n_present / n_defined)
coverage %>%
  tidyr::pivot_wider(id_cols = c(set, n_defined), names_from = cohort,
                     values_from = frac) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

low <- coverage %>% dplyr::filter(frac < 0.80)
if (nrow(low)) {
  warning(nrow(low), " set/cohort combination(s) below 0.80 coverage: ",
          paste(unique(low$set), collapse = ", "),
          ". Reported, not fatal - E03 flags them.", call. = FALSE)
}

# =============================================================================
# 7. Save
# =============================================================================
message("\n7. save")

saveRDS(list(
  gsva_arms    = gsva_s[names(arm_s), , drop = FALSE],
  gsva_cov     = gsva_s[names(cov_s), , drop = FALSE],
  gsva_new     = gsva_s[names(new_s), , drop = FALSE],
  mitopps_arms = mitopps_s,
  content_arms = content_s,
  zmean_arms   = zmean_s,
  mitopps_universe = mpps_univ_s,
  M_b_variants = M_b_s, log2MYC = log2MYC_s,
  symbol_map = map_s, symbol_report = bm_s$report,
  scale = list(gsva = "log_vst", mitopps = "linear_deseq2_normalised",
               content = "log2 of summed linear", zmean = "log_vst"),
  built = Sys.time()), PATH_SCANB_SCORES)

saveRDS(list(
  tcga_gsva_new = gsva_t_new,
  tcga_M_b_variants = M_b_t,
  tcga_log2MYC  = log2MYC_t,
  tcga_symbol_map = map_t,
  note = paste("TCGA's 18 mitochondrial arms on 4 instruments live in the",
               "snapshot (tcga_brca_mito_scores.rds) and are NOT rescored."),
  built = Sys.time()), PATH_NEW_SETS)

saveRDS(list(
  myc_sets = myc_sets, myc_panel = myc_panel, myc_base = myc_base,
  collectri_sets = collectri_sets,
  strip_refs = list(MITOCARTA_ALL = MITOCARTA_ALL, PROLIF_REF = PROLIF_REF),
  cdc_sets = cdc_sets, cdc_table = cdc_tab,
  tang_sets = tang_sets, tang_table = tang_tab,
  family_labels = family_labels, bcl2_family = BCL2_FAMILY,
  arm_sets = arm_sets, cov_sets = cov_sets,
  coverage = coverage,
  floors = c(min_set_genes = MIN_SET_GENES, min_score_n = MIN_SCORE_N,
             min_gsva_n = MIN_GSVA_N),
  rules = list(
    death_provenance = paste("human-native columns; no ortholog mapping;",
                             "see data/genesets_celldeath_human/README.md"),
    cicd = paste("pro-survival CICD (4) and CICD_MITO (2) are NOT scored;",
                 "carried as gene lists"),
    huge = paste("TANG_AUTOPHAGY_DEPENDENT_CELL_DEATH and TANG_FERROPTOSIS are",
                 "5-6.5% of the matrix; read against a size-matched comparator"),
    myc  = paste("report the panel ordered by proliferation entanglement,",
                 "never a single signature"),
    naming = paste("every MYC estimator carries an explicit suffix:",
                   "__FULL as distributed, __MITOSTRIP minus MITOCARTA_ALL,",
                   "__PROLIFSTRIP minus HALLMARK E2F+G2M, __BOTHSTRIP minus",
                   "both. There are no bare names. FELSHER__MITOSTRIP IS the",
                   "validation study's M_a and reproduces it exactly."),
    prolif_disjoint = paste("PROLIF_DISJOINT is PROLIF_STD minus the 9",
                            "proliferation genes of FELSHER_61 and nothing",
                            "else. It is disjoint from M_a ALONE and shares 97",
                            "genes with the 18-signature union, so a",
                            "proliferation-ADJUSTED rho for any other signature",
                            "is partly adjusting it for itself. Use the",
                            "__PROLIFSTRIP variants to ask it cleanly.")),
  built = Sys.time()), PATH_SETDEFS)

message("\nE02: done.")
message("    results/scanb_scores.rds       SCAN-B, 4 instruments + new sets")
message("    results/new_set_scores.rds     TCGA, the new sets only")
message("    results/set_definitions.rds    every set, panel and floor")
message("    NEXT: E03, the correlation atlas.")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  s <- readRDS(PATH_SCANB_SCORES)
  n <- readRDS(PATH_NEW_SETS)
  d <- readRDS(PATH_SETDEFS)
  m <- readRDS(PATH_TCGA_MITO)

  # --- did SCAN-B reproduce TCGA's shape? ----------------------------------
  # NOT the values - GSVA is cohort-relative and mitoPPS is composition. The
  # question is whether the ORDERING of arms by MYC correlation agrees.
  ma_t <- readRDS(PATH_TCGA_MYC)$estimators
  ma_t <- ma_t$M_a[match(colnames(m$gsva_arms), ma_t$patient)]
  ma_s <- s$gsva_new[MYC_REF, ]
  cmp <- data.frame(
    arm   = rownames(m$gsva_arms),
    TCGA  = round(sapply(rownames(m$gsva_arms),
                         function(a) .rho(m$gsva_arms[a, ], ma_t)), 3),
    SCANB = round(sapply(rownames(m$gsva_arms),
                         function(a) .rho(s$gsva_arms[a, ], ma_s)), 3))
  cmp[order(-cmp$TCGA), ]
  with(cmp, plot(TCGA, SCANB, pch = 16, xlim = c(-.2,.7), ylim = c(-.2,.7),
                 xlab = "rho(M_a, arm) TCGA", ylab = "SCAN-B"))
  abline(0, 1, lty = 2); with(cmp, text(TCGA, SCANB, arm, pos = 4, cex = .55))

  # --- the MYC panel, ordered by entanglement ------------------------------
  # Trap 3. If the OXPHOS correlation decays as frac_prolif falls, it is
  # proliferation.
  d$myc_panel %>% as.data.frame() %>% print(row.names = FALSE)

  # --- the death sets that were and were not scored ------------------------
  d$cdc_table %>% as.data.frame() %>% print(row.names = FALSE)
  d$tang_table %>% dplyr::arrange(dplyr::desc(n)) %>% as.data.frame() %>%
    print(row.names = FALSE)

  # --- the symbol map, which is the thing most likely to be wrong ----------
  table(s$symbol_report$status)
  s$symbol_report %>% dplyr::filter(status == "resolved") %>%
    as.data.frame() %>% head(25)

  # --- coverage ------------------------------------------------------------
  d$coverage %>% dplyr::filter(frac < 0.9) %>% as.data.frame() %>%
    print(row.names = FALSE)

}
