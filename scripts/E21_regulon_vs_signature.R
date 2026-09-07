# E21_regulon_vs_signature.R
# =============================================================================
# WHY DOES THE REGULON ESTIMATOR DISSENT? A NEW ANALYSIS, NOT A RE-READ.
#
# E20 found that "does MYC track compartment CONTENT or SHARE" is
# ESTIMATOR-DEPENDENT: the three GSVA signature estimators say content, and the
# CollecTRI regulon estimator `M_b` says the opposite, in both cohorts,
# materially. The statistic that carries that is
#
#     d = rho(est, ox_lvl) - rho(est, ox_rel)
#
# positive on every signature (+0.083 to +0.142) and negative on the regulon
# (-0.089 / -0.098). This script asks WHICH PROPERTY OF `M_b` PRODUCES IT.
#
# `M_b` differs from the signatures in FOUR ways at once, and nobody has ever
# separated them:
#
#   (1) INSTRUMENT   M_b is decoupleR ULM - a per-sample regression of the whole
#                    expression vector on a target indicator. The signatures are
#                    GSVA, a rank-based enrichment. Different mathematics.
#   (2) GENE SET     811 curated TF-target edges against 61-178 experimentally
#                    derived co-expression signature genes.
#   (3) SIGN         ULM is signed: mor = +1 on a stimulatory edge, -1 on an
#                    inhibitory one. No GSVA set carries a sign.
#   (4) SIZE         811 genes is 4.5% of the matrix. A large mean-z composite
#                    carries a global-expression component that a 61-gene one
#                    does not, and `ox_rel` is a WITHIN-SAMPLE CONTRAST that
#                    cancels exactly that component. CLAUDE.md trap 10.
#
# THE DESIGN. Score every gene set on ONE COMMON INSTRUMENT, so the gene set
# varies and nothing else; then apply the OTHER instrument to the signature
# sets, so the instrument varies and nothing else. The common instrument is
# `zmean` - mean of per-gene z - because it is what the four rulers are built
# from and what E06 already used to split this regulon.
#
#   gene set                 zmean      ULM        GSVA
#   regulon_all   (811)      built      READ M_b   -
#   regulon_stim  (736)      built      built      -      <- "COLLECTRI_MYC_STIM"
#   regulon_repr  (72)       built      built      -
#   regulon_diff  (stim-repr) built     -          -      <- E06's contrast
#   sig_msigdb    (177)      built      built      READ
#   sig_felsher   (61)       built      built      READ
#   sig_lowent    (178)      built      built      READ
#
# GSVA IS NOT RUN HERE. Using `zmean` in its place is licensed only if the two
# agree on the sets where both exist, and section 6 MEASURES that rather than
# assuming it. If they disagree, the instrument leg is reported as incomplete
# and the missing run is named.
#
# =============================================================================
# THE CONTROL THAT MAY MATTER MORE THAN THE DECOMPOSITION
# =============================================================================
# `ox_lvl` is a mean-z composite and `ox_rel` is `ox_lvl` minus a 1,047-gene
# mean-z. So `ox_rel` cancels whatever the two composites share - and what a
# LARGE random gene set shares with any other large composite is a global
# expression component that has nothing to do with MYC.
#
# THEREFORE: a size-matched RANDOM gene set may produce a positive `d` all on
# its own. If it does, "the signatures say content" is partly an artefact of
# composite size, and the estimator that dissents is the informative one rather
# than the odd one out. Section 7 draws random sets at each observed size from
# the non-MitoCarta universe and reports the null distribution of `d`.
#
# This control is not a formality. It can overturn E20's reading, and the rule
# for when it does is fixed in section 0 before any number exists.
#
# =============================================================================
# SCALES, STATED ONCE
# =============================================================================
#   the four rulers, and every zmean estimator
#                    BUILT from log2(linear DESeq2-normalised + 1) - E16/E20's
#                    `.comp`, copied verbatim
#   the E06 control  BUILT from the VST, because that is what E06 used and a
#                    reproduction check must use the same input
#   ULM              BUILT from the VST, because that is what E02 fed run_ulm
#   GSVA, M_b        READ from their saved objects; neither is recomputed except
#                    as the reproduction control in section 4
#
# Every correlation is rank-based. NEVER ACROSS COHORTS, and the mouse is not
# read by this script at all. SCAN-B is resolved through
# `scanb_scores.rds$symbol_map` throughout.
#
# ADJUSTMENT: partial Spearman on PROLIF_DISJOINT, as E16, E20 and E10. Raw is
# reported beside it. N3 THROUGHOUT: these are transcript correlations.
# SPECIES: human. No ortholog function is called here or anywhere in this repo.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE21: why does the regulon estimator dissent?\n", strrep("=", 78))

PATH_E21       <- file.path(DIR_RESULTS, "regulon_vs_signature.rds")
PATH_E21_PANEL <- file.path(DIR_TABLES,  "E21_regulon_vs_signature.csv")
PATH_E20       <- file.path(DIR_RESULTS, "share_vs_content.rds")
PATH_E06       <- file.path(DIR_RESULTS, "estimator_anatomy.rds")

# =============================================================================
# 0. THE RULES, FIXED BEFORE ANY NUMBER IN THIS SCRIPT EXISTS
# =============================================================================
# THE STATISTIC. One number per estimator per cohort, and the whole script is
# about its sign:
#     d = rho(est, ox_lvl) - rho(est, ox_rel)
# POSITIVE = the signature behaviour ("content"). NEGATIVE = the M_b dissent.
D_LABEL <- "d = rho(ox_lvl) - rho(ox_rel); positive = content, negative = the M_b dissent"

# Materiality, borrowed from E16's R1_MAX_DELTA exactly as E20 borrowed it.
MATERIAL_DELTA <- 0.05
BOOT_B         <- 2000L

# THE ATTRIBUTION RULE. Applied mechanically in section 8.
#   SIGN      if regulon_stim flips d positive while regulon_all stays negative.
#             This is the branch in which "score COLLECTRI_MYC_STIM" RESOLVES it.
#   GENE SET  if the regulon sets stay negative on zmean while the signature
#             sets stay positive on zmean - i.e. the dissent survives holding
#             the instrument fixed.
#   INSTRUMENT if the signature sets go negative under ULM.
#   SIZE      if the size-matched random null already produces the signature's
#             own d, so "content" is not a MYC property at all (section 7).
# More than one may fire. They are not exclusive and the script does not force
# a single winner.
ATTRIB_RULE <- paste0(
  "SIGN if regulon_stim flips d positive while regulon_all stays negative; ",
  "GENE SET if the regulon sets stay negative on zmean and the signatures stay ",
  "positive on zmean; INSTRUMENT if the signatures go negative under ULM; ",
  "SIZE if the size-matched random null reproduces the signatures' own d. ",
  "More than one may fire; no single winner is forced.")

# THE SIZE NULL, and what it would mean. Fixed here so section 7 cannot be
# read after the fact as a footnote.
#   If the random null at a signature's size has a median d of the same sign AND
#   within MATERIAL_DELTA of that signature's observed d, the signature's d is
#   NOT evidence about MYC. Reported per estimator, not pooled.
N_RANDOM   <- 200L
SIZE_RULE  <- paste0(
  "a signature's d is NOT evidence about MYC if the size-matched random null's ",
  "median d has the same sign and lies within ", MATERIAL_DELTA, " of it.")

message("\n0. rules fixed before any number")
message("   statistic       ", D_LABEL)
message("   attribution     ", ATTRIB_RULE)
message("   size null       ", N_RANDOM, " random sets per observed size. ",
        SIZE_RULE)
message("   materiality     |d| >= ", MATERIAL_DELTA,
        "; contrast intervals from a paired bootstrap, B = ", BOOT_B)

# =============================================================================
# 1. Inputs
# =============================================================================
message("\n1. inputs")

for (p in c(PATH_E20, PATH_E06)) if (!file.exists(p)) {
  stop("missing ", basename(p), ". E21 reproduces both E20's panel and E06's ",
       "regulon halves as its controls and cannot run without them.",
       call. = FALSE)
}
e20 <- readRDS(PATH_E20)
e06 <- readRDS(PATH_E06)
stopifnot("panel" %in% names(e20), "regulon" %in% names(e06),
          "halves" %in% names(e06))

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))

ID_T <- colnames(mito$gsva_arms)
ID_S <- colnames(sc$gsva_arms)

# --- 1.1 the gene sets -------------------------------------------------------
# CONSUMED, NOT REBUILT. The two regulon halves are E06's saved objects - the
# same 736 / 72 every earlier note refers to - and the signature sets are the
# snapshot's. Rebuilding either would put a second definition in the repo.
POS <- e06$regulon$activated      # 736, the "COLLECTRI_MYC_STIM" object
NEG <- e06$regulon$repressed      # 72
MB_ALL <- sd_$collectri_sets[[MB_REF]]
stopifnot(length(POS) == 736L, length(NEG) == 72L,
          length(MB_ALL) == EXPECT_COLLECTRI_STRIP,
          length(intersect(POS, NEG)) == 0L,
          all(c(POS, NEG) %in% MB_ALL))
message("   regulon halves from E06: activated ", length(POS),
        " (", e06$regulon$n_both_flagged, " BOTH-flagged, counted as activating)",
        " | repressed ", length(NEG))
message("   M_b__MITOSTRIP ", length(MB_ALL), " targets; ",
        length(MB_ALL) - length(POS) - length(NEG),
        " carry no usable edge after the distinct() and are in neither half")

SET_GENES <- list(
  regulon_all  = MB_ALL,
  regulon_stim = POS,
  regulon_repr = NEG,
  sig_msigdb   = sd_$myc_sets[["HALLMARK_MYC_TARGETS_V1__MITOSTRIP"]],
  sig_felsher  = sd_$myc_sets[[MYC_REF]],
  sig_lowent   = sd_$myc_sets[["MYC_UP.V1_UP__MITOSTRIP"]])
SET_FAMILY <- c(regulon_all = "regulon", regulon_stim = "regulon",
                regulon_repr = "regulon", sig_msigdb = "signature",
                sig_felsher = "signature", sig_lowent = "signature")

# The rulers' own genes, so nothing is measured partly against itself.
OX_SUB   <- mito$arm_sets[["OXPHOS subunits"]]
MITO_ALL <- sd_$strip_refs$MITOCARTA_ALL
REST     <- setdiff(MITO_ALL, OX_SUB)
MITORIB  <- mito$arm_sets[["Mitochondrial ribosome"]]
PROLIF_COV <- "PROLIF_DISJOINT"
PD <- mito$covariate_sets[[PROLIF_COV]]

overlap_audit <- tibble::tibble(
  set = names(SET_GENES), n = lengths(SET_GENES),
  x_numerator   = vapply(SET_GENES, function(s) length(intersect(s, OX_SUB)), integer(1)),
  x_rest_mito   = vapply(SET_GENES, function(s) length(intersect(s, REST)), integer(1)),
  x_mitoribosome = vapply(SET_GENES, function(s) length(intersect(s, MITORIB)), integer(1)),
  x_prolif      = vapply(SET_GENES, function(s) length(intersect(s, PD)), integer(1)))
message("\n   overlap of every estimator gene set with the rulers' own genes:")
overlap_audit %>% as.data.frame() %>% print(row.names = FALSE)
# All six are __MITOSTRIP or subsets of one, and MITOCARTA_ALL is the strip set,
# so the three mitochondrial columns are 0 BY CONSTRUCTION. Asserted, because a
# non-zero would mean an estimator is partly correlated with its own ruler.
stopifnot(all(overlap_audit$x_numerator == 0L),
          all(overlap_audit$x_rest_mito == 0L),
          all(overlap_audit$x_mitoribosome == 0L))
message("   every mitochondrial column is 0 by construction (MITOCARTA_ALL is",
        " the strip set).\n   the proliferation column is NOT zero and is why",
        " the adjusted panel is the one read.")

# --- 1.2 the matrices --------------------------------------------------------
tcga_lin  <- readRDS(PATH_TCGA_LINEAR)
scanb_lin <- readRDS(PATH_SCANB_LINEAR)
stopifnot(identical(tcga_lin$scale, "linear_deseq2_normalised"),
          identical(scanb_lin$scale, "linear_deseq2_normalised"))
GT <- log2(tcga_lin$mat[, ID_T, drop = FALSE] + 1)
GS <- log2(scanb_lin$mat[, ID_S, drop = FALSE] + 1)
rm(tcga_lin, scanb_lin); invisible(gc(verbose = FALSE))

RES_T <- .symbol_resolver(rownames(GT), NULL)
RES_S <- .symbol_resolver(rownames(GS), sc$symbol_map)
COH <- list(
  TCGA     = list(G = GT, res = RES_T, ids = ID_T, arms = mito,
                  gsva_new = nw$tcga_gsva_new, mb = nw$tcga_M_b_variants),
  `SCAN-B` = list(G = GS, res = RES_S, ids = ID_S, arms = sc,
                  gsva_new = sc$gsva_new, mb = sc$M_b_variants))
message("\n   TCGA ", length(ID_T), " samples | SCAN-B ", length(ID_S), " samples")

# =============================================================================
# 2. The rulers - E20's, rebuilt to E20's recipe and checked against it
# =============================================================================
message("\n2. the rulers")

.comp <- function(genes, C) {
  gr <- .gene_rows(genes, C$G, C$res)
  M  <- gr$mat
  v  <- apply(M, 1L, stats::var)
  colMeans(t(scale(t(M[v > 0, , drop = FALSE]))))
}

rulers <- lapply(COH, function(C) {
  num <- .comp(OX_SUB, C); mtr <- .comp(MITORIB, C)
  M <- rbind(ox_lvl       = num,
             ox_rel       = num - .comp(REST, C),
             ox_nuc_mtrib = num - mtr,
             mtrib_lvl    = mtr)
  colnames(M) <- C$ids
  stopifnot(!anyNA(M)); M
})
message("   four rulers rebuilt to E20's recipe; section 4 checks them against",
        " E20's saved panel")

# =============================================================================
# 3. The estimators - one gene set at a time, on each instrument
# =============================================================================
# INSTRUMENT 1, `zmean`: the common instrument. Identical construction to the
# rulers, so a gene-set comparison across it varies the gene set and nothing
# else.
#
# INSTRUMENT 2, `ulm`: decoupleR's univariate linear model, on the VST, with
# E02's network and E02's sign rule reproduced exactly. For the five UNSIGNED
# sets every mor is +1 - that is what "the same set through the other
# instrument" means. For `regulon_all` the real signed mor is used, because
# that IS `M_b` and section 4 requires it to come back bit-equal.
message("\n3. the estimators")

# --- 3.1 zmean, on log2(linear + 1) -----------------------------------------
zmean <- lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  M <- do.call(rbind, lapply(SET_GENES, function(g) .comp(g, C)))
  rownames(M) <- names(SET_GENES)
  # E06's contrast, built from the two halves rather than re-derived.
  M <- rbind(M, regulon_diff = M["regulon_stim", ] - M["regulon_repr", ])
  colnames(M) <- C$ids; M
})
names(zmean) <- names(COH)
message("   zmean: ", nrow(zmean$TCGA), " estimators per cohort (6 sets + the",
        " stim-minus-repr contrast)")

# --- 3.2 ULM, on the VST, with E02's network and sign rule ------------------
# E02 section 3.1b, transcribed. The sign rule is REPRODUCED, not improved on:
# mor = +1 where is_stimulation, -1 otherwise, and an edge flagged BOTH takes
# +1. That rule is why `M_b` is empirically its own activated half (E06: rho
# 0.955 / 0.934), and changing it here would break the section 4 control.
ct <- readr::read_tsv(PATH_COLLECTRI, show_col_types = FALSE, progress = FALSE)
myc_net <- ct %>%
  dplyr::filter(source_genesymbol == "MYC", !is.na(target_genesymbol),
                target_genesymbol != "") %>%
  dplyr::transmute(source = "MYC", target = target_genesymbol,
                   mor = dplyr::if_else(as.logical(is_stimulation), 1, -1),
                   likelihood = 1) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)
rm(ct); invisible(gc(verbose = FALSE))

.ulm_one <- function(targets, E, signed, remap) {
  tg <- remap(unique(targets))
  net <- if (signed) {
    myc_net %>% dplyr::filter(target %in% targets, target %in% rownames(E))
  } else {
    tibble::tibble(source = "MYC", target = intersect(tg, rownames(E)),
                   mor = 1, likelihood = 1)
  }
  res <- decoupleR::run_ulm(mat = E, network = net, .source = source,
                            .target = target, .mor = mor, minsize = 5L)
  v <- res %>% dplyr::filter(statistic == "ulm", source == "MYC") %>%
    dplyr::select(condition, score)
  stats::setNames(v$score, v$condition)[colnames(E)]
}

# The VST is loaded ONE COHORT AT A TIME and dropped immediately - both VSTs
# and both linear matrices together is more memory than this needs.
ulm <- list(); e06_control <- list()
.zmean_vst <- function(genes, E, remap) {
  g <- intersect(remap(unique(genes)), rownames(E))
  sub <- E[g, , drop = FALSE]
  v <- apply(sub, 1L, stats::var); sub <- sub[v > 0, , drop = FALSE]
  colMeans((sub - rowMeans(sub)) / apply(sub, 1L, stats::sd))
}
for (coh in names(COH)) {
  C <- COH[[coh]]
  vst <- readRDS(if (coh == "TCGA") PATH_TCGA_VST else PATH_SCANB_VST)
  stopifnot(identical(vst$scale, "log_vst"))
  E <- vst$mat[, C$ids, drop = FALSE]
  rm(vst); invisible(gc(verbose = FALSE))
  remap <- if (coh == "TCGA") function(g) g else
    function(g) { h <- sc$symbol_map[g]; unname(ifelse(is.na(h), g, h)) }

  message("   ULM, ", coh, ": ", length(SET_GENES), " sets")
  M <- do.call(rbind, lapply(names(SET_GENES), function(nm)
    .ulm_one(SET_GENES[[nm]], E, signed = identical(nm, "regulon_all"), remap)))
  dimnames(M) <- list(names(SET_GENES), C$ids)
  ulm[[coh]] <- M

  # E06's halves, on E06's own input, for the section 4 control.
  H <- rbind(regulon_activated = .zmean_vst(POS, E, remap),
             regulon_repressed = .zmean_vst(NEG, E, remap))
  H <- rbind(H, regulon_difference = H[1, ] - H[2, ])
  colnames(H) <- C$ids
  e06_control[[coh]] <- H
  rm(E); invisible(gc(verbose = FALSE))
}

# --- 3.3 GSVA, READ - the incumbent instrument, for the three signatures -----
gsva_est <- lapply(COH, function(C) {
  M <- rbind(
    sig_msigdb  = as.numeric(C$gsva_new["HALLMARK_MYC_TARGETS_V1__MITOSTRIP", C$ids]),
    sig_felsher = as.numeric(C$gsva_new[MYC_REF, C$ids]),
    sig_lowent  = as.numeric(C$gsva_new["MYC_UP.V1_UP__MITOSTRIP", C$ids]),
    regulon_all = as.numeric(C$mb[MB_REF, C$ids]))   # ULM, not GSVA - labelled
  colnames(M) <- C$ids; M
})
message("   read: 3 GSVA signatures + the saved M_b ULM, exactly the four",
        " estimators E20 used")

# =============================================================================
# 4. THE CONTROLS
# =============================================================================
message("\n4. controls")

# --- 4.1 does the rebuilt ULM reproduce the saved M_b? ----------------------
mb_repro <- dplyr::bind_rows(lapply(names(COH), function(coh) tibble::tibble(
  cohort = coh,
  max_abs_diff = max(abs(ulm[[coh]]["regulon_all", ] -
                           as.numeric(COH[[coh]]$mb[MB_REF, COH[[coh]]$ids]))),
  rho = .rho(ulm[[coh]]["regulon_all", ],
             as.numeric(COH[[coh]]$mb[MB_REF, COH[[coh]]$ids])))))
mb_repro %>% dplyr::mutate(rho = round(rho, 6)) %>%
  as.data.frame() %>% print(row.names = FALSE)
if (max(mb_repro$max_abs_diff) > 1e-8) {
  stop("the rebuilt signed ULM does not reproduce the saved M_b__MITOSTRIP ",
       "(max abs diff ", signif(max(mb_repro$max_abs_diff), 3), "). The network ",
       "or the sign rule has drifted from E02; stop and reconcile.",
       call. = FALSE)
}
message("   4.1 the rebuilt signed ULM IS M_b, bit-equal in both cohorts")

# --- 4.2 do the rebuilt halves reproduce E06's? -----------------------------
halves_repro <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  A <- e06_control[[coh]]; B <- e06$halves[[coh]][, COH[[coh]]$ids, drop = FALSE]
  tibble::tibble(cohort = coh, max_abs_diff = max(abs(A[rownames(B), ] - B)))
}))
halves_repro %>% as.data.frame() %>% print(row.names = FALSE)
if (max(halves_repro$max_abs_diff) > 1e-10) {
  stop("the rebuilt regulon halves do not reproduce E06's `$halves` ",
       "(max abs diff ", signif(max(halves_repro$max_abs_diff), 3), "). ",
       "The halves or the zmean recipe has drifted; stop and reconcile.",
       call. = FALSE)
}
message("   4.2 the two regulon halves reproduce E06 bit-equal")

# =============================================================================
# 5. THE PANEL, and the statistic
# =============================================================================
# `.cor_block` is E10's, copied verbatim as E16, E18, E19 and E20 copy it.
message("\n5. the panel")

.rank_rows <- function(M) {
  out <- t(apply(M, 1L, rank))
  if (nrow(M) == 1L) out <- matrix(out, nrow = 1L)
  dimnames(out) <- dimnames(M); out
}
.cor_block <- function(A, B, cov = NULL) {
  n <- ncol(A)
  stopifnot(identical(colnames(A), colnames(B)))
  RA <- .rank_rows(A); RB <- .rank_rows(B)
  k  <- 0L
  if (!is.null(cov)) {
    stopifnot(nrow(cov) == n)
    H  <- qr(cbind(`(Intercept)` = 1, apply(cov, 2L, rank)))
    k  <- ncol(cov)
    RA <- RA - t(qr.fitted(H, t(RA)))
    RB <- RB - t(qr.fitted(H, t(RB)))
  }
  R  <- suppressWarnings(stats::cor(t(RA), t(RB)))
  z  <- atanh(pmin(pmax(R, -0.999999999), 0.999999999))
  se <- sqrt((1 + R^2 / 2) / (n - 3 - k))
  tibble::tibble(
    ruler = rep(rownames(R), times = ncol(R)),
    gene  = rep(colnames(R), each  = nrow(R)),
    n = n, k_cov = k, rho = as.vector(R),
    ci_lo = as.vector(tanh(z - 1.959964 * se)),
    ci_hi = as.vector(tanh(z + 1.959964 * se)))
}

CV <- lapply(names(COH), function(coh)
  t(COH[[coh]]$arms$gsva_cov[PROLIF_COV, COH[[coh]]$ids, drop = FALSE]))
names(CV) <- names(COH)

EST <- lapply(names(COH), function(coh) {
  Z <- zmean[[coh]]; U <- ulm[[coh]]; G <- gsva_est[[coh]]
  rownames(Z) <- paste0(rownames(Z), "|zmean")
  rownames(U) <- paste0(rownames(U), "|ulm")
  rownames(G) <- paste0(rownames(G), "|gsva")
  # `regulon_all|gsva` would be a lie - that row is the saved ULM. Renamed so
  # no table can imply a GSVA score for the regulon exists.
  rownames(G)[rownames(G) == "regulon_all|gsva"] <- "regulon_all|ulm_saved"
  rbind(Z, U, G)
})
names(EST) <- names(COH)

panel <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  M <- rulers[[coh]]; E <- EST[[coh]]
  dplyr::bind_rows(
    .cor_block(M, E)                  %>% dplyr::mutate(adjustment = "raw"),
    .cor_block(M, E, cov = CV[[coh]]) %>% dplyr::mutate(adjustment = "PROLIF_DISJOINT")) %>%
    dplyr::mutate(cohort = coh)
})) %>%
  dplyr::rename(estimator = gene) %>%
  tidyr::separate_wider_delim(estimator, delim = "|",
                              names = c("set", "instrument"), cols_remove = FALSE) %>%
  dplyr::mutate(family = dplyr::if_else(grepl("^regulon", set), "regulon", "signature"))

# --- 5.1 REPRODUCTION AGAINST E20 -------------------------------------------
# The three GSVA signatures and the saved M_b, on ox_lvl / ox_rel /
# ox_nuc_mtrib / mtrib_lvl, are exactly cells E20 already carries. Same vectors,
# same covariate, same estimator: they must come back bit-equal.
E20_MAP <- c(sig_msigdb = "myc_msigdb", sig_felsher = "myc_felsher",
             sig_lowent = "myc_lowent", regulon_all = "myc_regulon")
repro20 <- panel %>%
  dplyr::filter(instrument %in% c("gsva", "ulm_saved")) %>%
  dplyr::mutate(against = unname(E20_MAP[set])) %>%
  dplyr::inner_join(
    e20$panel %>%
      dplyr::mutate(ruler = as.character(ruler), against = as.character(against)) %>%
      dplyr::select(cohort, adjustment, ruler, against, e20_rho = rho),
    by = c("cohort", "adjustment", "ruler", "against")) %>%
  dplyr::mutate(delta = rho - e20_rho)
stopifnot(nrow(repro20) == 2L * 2L * 4L * 4L)
if (max(abs(repro20$delta)) >= 1e-12) {
  stop("E21 does not reproduce E20's cells (max |delta| = ",
       signif(max(abs(repro20$delta)), 3), "). An input has moved; stop.",
       call. = FALSE)
}
message("   5.1 ", nrow(repro20), " of E20's cells reproduce exactly (max |delta| = ",
        signif(max(abs(repro20$delta)), 3), ")")

# --- 5.2 the statistic ------------------------------------------------------
.dstat <- function(P) P %>%
  dplyr::select(cohort, adjustment, set, instrument, estimator, family, ruler, rho) %>%
  tidyr::pivot_wider(names_from = ruler, values_from = rho) %>%
  dplyr::mutate(d = ox_lvl - ox_rel,
                d_mtrib = ox_lvl - ox_nuc_mtrib,
                mtrib_pref = mtrib_lvl - ox_lvl)
dstat <- .dstat(panel)

for (adj in c("raw", "PROLIF_DISJOINT")) {
  message("\n   ", adj, " - rho against each ruler, and d:")
  dstat %>%
    dplyr::filter(adjustment == adj) %>%
    dplyr::transmute(cohort, estimator,
                     says = dplyr::if_else(d > 0, "content", "SHARE"),
                     ox_lvl = sprintf("%+.3f", ox_lvl),
                     ox_rel = sprintf("%+.3f", ox_rel),
                     mtrib_lvl = sprintf("%+.3f", mtrib_lvl),
                     d = sprintf("%+.3f", d)) %>%
    dplyr::arrange(cohort, estimator) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

# --- 5.3 paired bootstrap on d ----------------------------------------------
message("\n   paired bootstrap on d, B = ", BOOT_B, " - the slow part")
.boot_d <- function(M, E, cov, B, seed) {
  n <- ncol(M); set.seed(seed)
  out <- list(raw = matrix(NA_real_, B, nrow(E)),
              adj = matrix(NA_real_, B, nrow(E)))
  for (b in seq_len(B)) {
    i  <- sample.int(n, n, replace = TRUE)
    RA <- .rank_rows(M[c("ox_lvl", "ox_rel"), i, drop = FALSE])
    RB <- .rank_rows(E[, i, drop = FALSE])
    r  <- suppressWarnings(stats::cor(t(RA), t(RB)))
    out$raw[b, ] <- r["ox_lvl", ] - r["ox_rel", ]
    H  <- qr(cbind(`(Intercept)` = 1, apply(cov[i, , drop = FALSE], 2L, rank)))
    r2 <- suppressWarnings(stats::cor(t(RA - t(qr.fitted(H, t(RA)))),
                                      t(RB - t(qr.fitted(H, t(RB))))))
    out$adj[b, ] <- r2["ox_lvl", ] - r2["ox_rel", ]
  }
  for (nm in names(out)) colnames(out[[nm]]) <- rownames(E)
  out
}
bd <- lapply(names(COH), function(coh)
  .boot_d(rulers[[coh]], EST[[coh]], CV[[coh]], BOOT_B, PROJECT_SEED))
names(bd) <- names(COH)

dstat <- dstat %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    boot_lo = stats::quantile(bd[[cohort]][[if (adjustment == "raw") "raw" else "adj"]][, estimator],
                              0.025, names = FALSE),
    boot_hi = stats::quantile(bd[[cohort]][[if (adjustment == "raw") "raw" else "adj"]][, estimator],
                              0.975, names = FALSE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(excludes_zero = (boot_lo > 0) | (boot_hi < 0),
                material = excludes_zero & abs(d) >= MATERIAL_DELTA)

message("\n   d with its paired bootstrap interval, adjusted:")
dstat %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
  dplyr::transmute(cohort, estimator, family,
                   d = sprintf("%+.3f", d),
                   boot = sprintf("[%+.3f, %+.3f]", boot_lo, boot_hi),
                   material) %>%
  dplyr::arrange(estimator, cohort) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. IS `zmean` A FAIR STAND-IN FOR GSVA?
# =============================================================================
# The whole gene-set leg rests on zmean, and GSVA is not run here. That is
# licensed only if the two instruments give the same answer where both exist -
# the three signature sets. MEASURED, not assumed.
message("\n6. is zmean a fair stand-in for GSVA?")

instr_check <- dstat %>%
  dplyr::filter(set %in% c("sig_msigdb", "sig_felsher", "sig_lowent"),
                instrument %in% c("zmean", "gsva")) %>%
  dplyr::select(cohort, adjustment, set, instrument, d) %>%
  tidyr::pivot_wider(names_from = instrument, values_from = d) %>%
  dplyr::mutate(gap = zmean - gsva, same_sign = sign(zmean) == sign(gsva))
instr_check %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
ZMEAN_OK <- all(instr_check$same_sign) &&
  max(abs(instr_check$gap)) < MATERIAL_DELTA
message("   zmean is a fair stand-in for GSVA on these sets: ", ZMEAN_OK,
        "  (all same sign: ", all(instr_check$same_sign),
        "; largest gap ", sprintf("%.3f", max(abs(instr_check$gap))),
        " against ", MATERIAL_DELTA, ")")
if (!ZMEAN_OK) {
  message("   -> the GENE-SET leg is INCOMPLETE. What is missing is a GSVA run",
          " over the\n      three regulon sets, in one batch with the pins,",
          " which this script does not do.")
}

# =============================================================================
# 7. THE SIZE NULL - and it can overturn the reading
# =============================================================================
# Random gene sets at each observed size, drawn from the SAME universe the
# estimators live in: matrix genes that are not MitoCarta members, so a random
# set cannot overlap either half of any ruler. zmean, the common instrument.
#
# If a random set of a signature's size already produces that signature's d,
# then "the signatures say content" is a property of composite SIZE against a
# within-sample contrast ruler, not a property of MYC.
message("\n7. the size-matched random null (", N_RANDOM, " sets per size)")

null_d <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  pool <- setdiff(rownames(C$G), C$res(MITO_ALL))
  R2   <- .rank_rows(rulers[[coh]][c("ox_lvl", "ox_rel"), , drop = FALSE])
  H    <- qr(cbind(`(Intercept)` = 1, rank(CV[[coh]][, 1])))
  R2a  <- R2 - t(qr.fitted(H, t(R2)))
  set.seed(PROJECT_SEED)
  dplyr::bind_rows(lapply(unique(lengths(SET_GENES)), function(k) {
    D <- vapply(seq_len(N_RANDOM), function(i) {
      g <- sample(pool, k)
      sub <- C$G[g, , drop = FALSE]
      v <- apply(sub, 1L, stats::var)
      s <- colMeans(t(scale(t(sub[v > 0, , drop = FALSE]))))
      rs <- rank(s)
      rr <- suppressWarnings(stats::cor(t(R2), rs))
      ra <- suppressWarnings(stats::cor(t(R2a), rs - qr.fitted(H, rs)))
      c(raw = rr["ox_lvl", 1] - rr["ox_rel", 1],
        adj = ra["ox_lvl", 1] - ra["ox_rel", 1])
    }, c(raw = NA_real_, adj = NA_real_))
    tibble::tibble(cohort = coh, size = k,
                   adjustment = rep(c("raw", "PROLIF_DISJOINT"), each = N_RANDOM),
                   d = c(D["raw", ], D["adj", ]))
  }))
}))

null_summary <- null_d %>%
  dplyr::group_by(cohort, adjustment, size) %>%
  dplyr::summarise(median_d = stats::median(d),
                   lo = stats::quantile(d, 0.025, names = FALSE),
                   hi = stats::quantile(d, 0.975, names = FALSE),
                   frac_positive = mean(d > 0), .groups = "drop")
message("\n   null distribution of d by set size, adjusted:")
null_summary %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(cohort, size) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Every observed zmean estimator against its own size-matched null.
vs_null <- dstat %>%
  dplyr::filter(instrument == "zmean", set != "regulon_diff") %>%
  dplyr::mutate(size = unname(lengths(SET_GENES)[set])) %>%
  dplyr::inner_join(null_summary, by = c("cohort", "adjustment", "size")) %>%
  dplyr::mutate(
    null_explains = sign(median_d) == sign(d) & abs(d - median_d) < MATERIAL_DELTA,
    outside_null  = d < lo | d > hi)
message("\n   each zmean estimator against its OWN size-matched null, adjusted:")
vs_null %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
  dplyr::transmute(cohort, set, size,
                   observed_d = sprintf("%+.3f", d),
                   null_median = sprintf("%+.3f", median_d),
                   null_95 = sprintf("[%+.3f, %+.3f]", lo, hi),
                   outside_null, null_explains) %>%
  dplyr::arrange(cohort, size) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 8. ATTRIBUTION, on the rule fixed in section 0
# =============================================================================
message("\n8. attribution")

A <- dstat %>% dplyr::filter(adjustment == "PROLIF_DISJOINT")
.d <- function(s, i, coh) A$d[A$set == s & A$instrument == i & A$cohort == coh]
.both <- function(s, i, f) all(vapply(names(COH), function(c) f(.d(s, i, c)), logical(1)))

sign_fires <- .both("regulon_stim", "zmean", function(x) x > 0) &&
  .both("regulon_all", "zmean", function(x) x < 0)
geneset_fires <- .both("regulon_all", "zmean", function(x) x < 0) &&
  all(vapply(c("sig_msigdb", "sig_felsher", "sig_lowent"), function(s)
    .both(s, "zmean", function(x) x > 0), logical(1)))
instrument_fires <- any(vapply(c("sig_msigdb", "sig_felsher", "sig_lowent"),
                               function(s) .both(s, "ulm", function(x) x < 0),
                               logical(1)))
size_fires <- vs_null %>%
  dplyr::filter(adjustment == "PROLIF_DISJOINT", family == "signature") %>%
  dplyr::pull(null_explains) %>% all()

attribution <- tibble::tibble(
  cause = c("SIGN (COLLECTRI_MYC_STIM resolves it)", "GENE SET", "INSTRUMENT",
            "SIZE (the null already gives the signatures' d)"),
  fires = c(sign_fires, geneset_fires, instrument_fires, size_fires))
attribution %>% as.data.frame() %>% print(row.names = FALSE)

message("\n   DOES SCORING COLLECTRI_MYC_STIM RESOLVE THE DISSENT? ",
        if (sign_fires) "YES" else "NO")
message("   regulon_stim, zmean, d = ",
        paste(sprintf("%s %+.3f", names(COH),
                      vapply(names(COH), function(c) .d("regulon_stim", "zmean", c),
                             numeric(1))), collapse = " | "))
message("   regulon_all,  zmean, d = ",
        paste(sprintf("%s %+.3f", names(COH),
                      vapply(names(COH), function(c) .d("regulon_all", "zmean", c),
                             numeric(1))), collapse = " | "))

verdicts <- tibble::tibble(
  item = c("does COLLECTRI_MYC_STIM resolve it", "SIGN", "GENE SET",
           "INSTRUMENT", "SIZE", "zmean a fair stand-in for GSVA"),
  value = as.character(c(sign_fires, sign_fires, geneset_fires,
                         instrument_fires, size_fires, ZMEAN_OK)))
message("")
verdicts %>% as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 9. Save
# =============================================================================
saveRDS(list(
  overlap_audit = overlap_audit, panel = panel, dstat = dstat,
  mb_repro = mb_repro, halves_repro = halves_repro, repro20 = repro20,
  instr_check = instr_check, zmean_ok = ZMEAN_OK,
  null_d = null_d, null_summary = null_summary, vs_null = vs_null,
  attribution = attribution, verdicts = verdicts,
  set_genes_n = lengths(SET_GENES), set_family = SET_FAMILY,
  settings = list(
    what_this_is = paste("a NEW analysis: which property of M_b - sign, gene",
                         "set, instrument or size - produces E20's dissent."),
    statistic = D_LABEL,
    instruments = c(zmean = "mean per-gene z of log2(linear DESeq2-normalised + 1)",
                    ulm = "decoupleR run_ulm on the VST, E02's network and sign rule",
                    gsva = "READ from E02; not recomputed",
                    ulm_saved = "the saved M_b__MITOSTRIP, read"),
    sign_rule = paste("mor = +1 where is_stimulation, -1 otherwise; an edge",
                      "flagged BOTH takes +1. E02's rule, reproduced not fixed."),
    covariate = PROLIF_COV,
    attribution_rule = ATTRIB_RULE, size_rule = SIZE_RULE,
    decision_rules = c(MATERIAL_DELTA = MATERIAL_DELTA, BOOT_B = BOOT_B,
                       N_RANDOM = N_RANDOM),
    seed = PROJECT_SEED),
  rules = list(
    scope = paste("No interaction model, no MYC stratification, no subtype",
                  "split, no new endpoint. No GSVA run - section 6 tests",
                  "whether that was allowed."),
    n3 = "these are transcript correlations; the word 'primed' is never used",
    cohorts = "never pooled and never compared as values",
    mouse = "the mouse is not read, compared or cited by this script"),
  built = Sys.time()), PATH_E21)

readr::write_csv(dstat, PATH_E21_PANEL)

message("\nE21: done.")
message("    results/regulon_vs_signature.rds")
message("    outputs/tables/E21_regulon_vs_signature.csv")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E21)

  # The answer.
  x$verdicts %>% as.data.frame()
  x$attribution %>% as.data.frame()

  # THE CONTROLS. All three must hold before anything else is read.
  x$mb_repro %>% as.data.frame()       # the rebuilt ULM IS M_b
  x$halves_repro %>% as.data.frame()   # the halves ARE E06's
  summary(abs(x$repro20$delta))        # E20's cells come back

  # d, every set on every instrument, adjusted. The whole script in one table.
  x$dstat %>%
    dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
    dplyr::transmute(cohort, set, instrument, family,
                     d = round(d, 3),
                     boot = sprintf("[%+.3f, %+.3f]", boot_lo, boot_hi),
                     material) %>%
    dplyr::arrange(set, instrument, cohort) %>%
    as.data.frame()

  # Was zmean allowed to stand in for GSVA?
  x$instr_check %>%
    dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

  # THE SIZE NULL. If the signatures sit inside their own null, E20's reading
  # is about composite size and not about MYC.
  x$vs_null %>%
    dplyr::filter(adjustment == "PROLIF_DISJOINT") %>%
    dplyr::transmute(cohort, set, size, d = round(d, 3),
                     null_median = round(median_d, 3),
                     null_95 = sprintf("[%+.3f, %+.3f]", lo, hi),
                     outside_null, null_explains) %>%
    as.data.frame()

  # The null itself, if the summary looks surprising.
  x$null_summary %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame()

}
