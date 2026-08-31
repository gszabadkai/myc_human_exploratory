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
scanb_vst <- readRDS(PATH_SCANB_VST)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
scanb_ph  <- readRDS(PATH_SCANB_PHENO)
g1        <- readRDS(PATH_G1)

stopifnot(identical(tcga_vst$scale, "log_vst"),
          identical(tcga_lin$scale, "linear_deseq2_normalised"),
          identical(scanb_vst$scale, "log_vst"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))

felsher_61 <- g1$estimators_stripped$FELSHER
collectri  <- g1$estimators_stripped$COLLECTRI_MYC_ALL
stopifnot(length(felsher_61) == EXPECT_FELSHER_STRIP)
message("   Felsher-61: ", length(felsher_61), " genes | CollecTRI: ",
        length(collectri))

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

# --- 2.2 the MYC estimator panel --------------------------------------------
# The 16-signature compendium (native human GMX) plus HALLMARK_MYC_TARGETS_V1,
# which is NOT in the GMX. Felsher-61 is scored alongside so this study's panel
# contains the validation study's M_a as one member rather than as the default.
#
# GMX layout: row 1 set names, row 2 a description line (all NA here), rows 3+
# genes, one column per set, NA-padded.
message("\n2.2 MYC estimator panel")

gmx <- utils::read.delim(PATH_MYC_GMX, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE, na.strings = c("NA", ""))
gmx <- gmx[-1, , drop = FALSE]
myc_sets <- lapply(gmx, function(v) sort(unique(v[!is.na(v) & v != ""])))
if (length(myc_sets) != EXPECT_MYC_SETS) {
  stop("the MYC GMX carries ", length(myc_sets), " signatures, expected ",
       EXPECT_MYC_SETS, ". Re-snapshot and update E00.", call. = FALSE)
}

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
.hm <- function(nm) unique(hallmark$gene_symbol[hallmark$gs_name == nm])
myc_sets[["HALLMARK_MYC_TARGETS_V1"]] <- .hm("HALLMARK_MYC_TARGETS_V1")
myc_sets[["FELSHER_61"]]              <- felsher_61

# Proliferation entanglement, the ordering variable for the whole panel
# (CLAUDE.md trap 3). Computed here so it travels with the scores.
prolif_hallmark <- unique(c(.hm("HALLMARK_E2F_TARGETS"),
                            .hm("HALLMARK_G2M_CHECKPOINT")))
myc_panel <- tibble::tibble(
  signature   = names(myc_sets),
  n           = vapply(myc_sets, length, integer(1)),
  n_prolif    = vapply(myc_sets, function(g) sum(g %in% prolif_hallmark), integer(1))) %>%
  dplyr::mutate(frac_prolif = n_prolif / n) %>%
  dplyr::arrange(frac_prolif)
myc_panel %>%
  dplyr::mutate(frac_prolif = round(100 * frac_prolif, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("   entanglement spans ", round(100 * min(myc_panel$frac_prolif), 1),
        "% to ", round(100 * max(myc_panel$frac_prolif), 1),
        "%. NEVER report a MYC-OXPHOS correlation from one signature.")

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
    cand <- setdiff(cand, genes)
    if (length(cand) == 1L) { out[[g]] <- cand; status[[i]] <- "resolved" }
    else if (length(cand) > 1L) status[[i]] <- "ambiguous"
  }
  list(map = out, report = tibble::tibble(input_symbol = genes, status = status,
                                          resolved_to = unname(out[genes])))
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
.mitopps_universe <- function(S) {
  N <- ncol(S); P <- nrow(S)
  Bi <- 1 / S
  A  <- (S %*% t(Bi)) / N
  out <- S * (((1 / A) %*% Bi) - Bi) / (P - 1)
  dimnames(out) <- dimnames(S)
  out
}
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
print(table(bm_s$report$status))
map_s <- bm_s$map

mapped <- unname(map_s)[unname(map_s) %in% SYM_S]
if (anyDuplicated(mapped)) {
  stop("the SCAN-B symbol map sends two inputs to the same row: ",
       paste(unique(mapped[duplicated(mapped)]), collapse = ", "), call. = FALSE)
}
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
.path_scores <- function(sets, L) {
  keep <- sets[vapply(sets, length, integer(1)) >= MIN_SET_GENES]
  out  <- t(vapply(keep, function(g) colMeans(L[g, , drop = FALSE]),
                   numeric(ncol(L))))
  rownames(out) <- names(keep); colnames(out) <- colnames(L)
  out
}
S_univ_s <- .path_scores(paths_s, Ls)
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
S_arms_s <- .path_scores(arm_s, Ls)
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
message("\n4.5 MYC: CollecTRI ULM and raw expression")
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
  message(sprintf("   %-14s %3d targets in the matrix", label, nrow(net)))
  res <- decoupleR::run_ulm(mat = E, network = net, .source = source,
                            .target = target, .mor = mor, minsize = 5L)
  v <- res %>% dplyr::filter(statistic == "ulm", source == "MYC") %>%
    dplyr::select(condition, score)
  stats::setNames(v$score, v$condition)[colnames(E)]
}
# NOTE: the CollecTRI targets are NOT in the symbol map's input union, so
# .remap_s leaves any unmatched one alone rather than resolving it. That is
# deliberate: widening the map's input universe can change a resolution for a
# symbol it already covers (the setdiff guard sees a different "already in the
# set" universe), which would break the agreement assertion in 4.1 for a
# secondary estimator's sake. M_b is reported alongside the panel, never alone.
M_b_s <- .ulm(.remap_s(collectri), Es, "SCAN-B M_b")
if (!"MYC" %in% rownames(Ls)) stop("MYC absent from SCAN-B.", call. = FALSE)
log2MYC_s <- as.numeric(log2(Ls["MYC", ]))
names(log2MYC_s) <- colnames(Ls)

message(sprintf("   rho(M_a GSVA, M_b ULM) in SCAN-B = %+.3f",
                .rho(gsva_s["FELSHER_61", ], M_b_s)))
message(sprintf("   rho(M_a GSVA, log2(MYC))         = %+.3f",
                .rho(gsva_s["FELSHER_61", ], log2MYC_s)))

# =============================================================================
# 5. TCGA: only the new sets
# =============================================================================
message("\n", strrep("=", 78), "\n5. TCGA: the new sets only\n", strrep("=", 78))

Et <- tcga_vst$mat
Lt <- tcga_lin$mat
SYM_T <- rownames(Et)

bm_t <- .build_symbol_map(unlist(NEW_SETS, use.names = FALSE), SYM_T)
print(table(bm_t$report$status))
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
  M_b = M_b_s, log2MYC = log2MYC_s,
  symbol_map = map_s, symbol_report = bm_s$report,
  scale = list(gsva = "log_vst", mitopps = "linear_deseq2_normalised",
               content = "log2 of summed linear", zmean = "log_vst"),
  built = Sys.time()), PATH_SCANB_SCORES)

saveRDS(list(
  tcga_gsva_new = gsva_t_new,
  tcga_log2MYC  = log2MYC_t,
  tcga_symbol_map = map_t,
  note = paste("TCGA's 18 mitochondrial arms on 4 instruments live in the",
               "snapshot (tcga_brca_mito_scores.rds) and are NOT rescored."),
  built = Sys.time()), PATH_NEW_SETS)

saveRDS(list(
  myc_sets = myc_sets, myc_panel = myc_panel,
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
                 "never a single signature")),
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
  ma_s <- s$gsva_new["FELSHER_61", ]
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
