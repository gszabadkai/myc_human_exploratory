# E38_normal_diagnostics.R
# =============================================================================
# FOUR DIAGNOSTICS UNDER E36'S "UNINTERPRETABLE"
#
# EXPLORATORY AND POST-HOC. Nothing here is pre-registered, and TWO OF THE FOUR
# LEGS ARE SECOND ATTEMPTS AFTER SEEING A FIRST RESULT - the luminal-only score
# after the declared nine-gene EPITHELIAL failed, and the compartment halves
# after the signed composite split by adjustment. Both are reported BESIDE the
# originals, never instead of them.
#
# Data note: docs/2026-10-01_e38_data.md at 94d5fd7, amended at a6125f7, BOTH
# BEFORE THIS FILE EXISTED. It names the arbitrary sign splits gene by gene,
# which is the only thing that makes diagnostic 1(b) worth running. Read it
# first; this header restates its rules and invents nothing.
#
# =============================================================================
# THIS IS A DIAGNOSTIC. IT CANNOT OVERTURN E36.
# =============================================================================
# E38 decomposes quantities E36 already computed. It adds no hypothesis and
# tests none. E36's verdict - UNINTERPRETABLE, on a declared control and a
# declared reading rule - STANDS AS RECORDED whatever appears below. Nothing
# here licenses a claim about normal breast, about the coupling, or about MYC.
#
# =============================================================================
# THE SNAPSHOT - SAY IT ONCE, LOUDLY
# =============================================================================
# Reads the JOINT snapshot E37 built, and nothing else:
#   results/joint_tcga_linear.rds   LINEAR DESeq2-normalised
#   results/joint_tcga_vst.rds      VST, log scale
#   results/joint_tcga_scores.rds   GSVA arms + covariates, mitoPPS arms
# NO VALUE HERE IS COMPARABLE WITH ONE FROM data/from_validation/.
#
# ONE OBJECT IS READ FROM data/from_validation/ AND IT IS NOT AN EXPRESSION
# VALUE: MYC_amp, a per-patient GISTIC call, for diagnostic 4. E37 section 4.2
# kept the genomics untouched on purpose and E36/E37 already carry
# BUFFER_gistic across by the same reasoning. No expression value crosses.
#
# =============================================================================
# WHAT IS FORBIDDEN
# =============================================================================
#   - No outcome, survival or treatment variable. None is read or joined.
#   - No product term. No causal language. No killing language.
#   - N3: transcript scores. Nothing is "primed", in either tissue.
#   - No per-gene FDR, and no p-value plucked from the grid.
#   - Nothing written to myc_human_validation or myc_mouse.
#
# =============================================================================
# THE FOUR DIAGNOSTICS, AND THE GATE
# =============================================================================
# 1  SIGNING. In normals every unsigned composite sits at +0.729 to +0.902
#    against OXPHOS while the signed configuration sits at -0.011. Is that
#    cancellation or biology? (a) the same six genes unsigned; (b) arbitrarily
#    SIGNED comparators, split by the gene lists fixed in the data note; and
#    (c) the sign-permutation of the configuration's own six genes, which needs
#    no threshold. DIAGNOSTIC 1 GATES THE READING OF 2 AND 3.
# 2  COMPARTMENT. E14's claim is about halves, not a signed composite. THE
#    CONFIGURATION HAS NO CYTOSOLIC HALF - all six genes are in MitoCarta - so
#    this runs on E14's own apoptotic machinery (44), 20/24, with CYCS dropped
#    as the one member also in the OXPHOS arm. Fe-S and mitophagy alongside.
# 3  LUMINAL-ONLY. E36's declared EPITHELIAL failed its own validity check
#    because its luminal and basal halves cancel. The luminal half alone, under
#    E36's criterion UNCHANGED. Post-hoc. Reported beside the declared version.
# 4  MYC TRANSCRIPT. M_b__PROLIFSTRIP reverses the paired MYC result. MYC,
#    MYCN and MYCL transcripts paired; the 82 stripped targets printed; MYC
#    amplification and composition of the 113.
#
# =============================================================================
# SCALE DISCIPLINE - CLAUDE.md
# =============================================================================
#   GSVA / mitoPPS   read as scored by E37, on their own opposite scales
#   all composites   log2(linear + 1), then z ACROSS THE JOINT 1,208, exactly
#                    as E36 built them. Spearman is invariant to that choice.
#   NO BLOM TRANSFORM ANYWHERE. E36's Q1 path applies none - the partial
#   Spearman ranks internally - so E38 applies none either.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

message("\nE38: four diagnostics under E36's UNINTERPRETABLE\n", strrep("=", 78))
message("POST-HOC AND EXPLORATORY. This cannot overturn E36.")

PATH_JOINT_LINEAR <- file.path(DIR_RESULTS, "joint_tcga_linear.rds")
PATH_JOINT_VST    <- file.path(DIR_RESULTS, "joint_tcga_vst.rds")
PATH_JOINT_SCORES <- file.path(DIR_RESULTS, "joint_tcga_scores.rds")
PATH_E14          <- file.path(DIR_RESULTS, "curated_comparators.rds")
PATH_COVARS <- here::here("data", "from_validation", "tcga_brca_covariates.rds")

PATH_E38         <- file.path(DIR_RESULTS, "e38_normal_diagnostics.rds")
PATH_E38_CSV     <- file.path(DIR_TABLES,  "E38_normal_side_correlations.csv")
PATH_E38_FIG     <- file.path(DIR_FIGURES, "E38_normal_diagnostics.png")
PATH_E38_FIG_PDF <- file.path(DIR_FIGURES, "E38_normal_diagnostics.pdf")
DIR_DOC_FIGURES  <- here::here("docs", "figures")
PATH_E38_FIG_DOC     <- file.path(DIR_DOC_FIGURES, "E38_normal_diagnostics.png")
PATH_E38_FIG_DOC_PDF <- file.path(DIR_DOC_FIGURES, "E38_normal_diagnostics.pdf")

# =============================================================================
# 0. CONSTANTS, AND THE GATE - all from the data note, nothing invented here
# =============================================================================
N_BOOT   <- 10000L
CI_LEVEL <- 0.95
SEED     <- 20261001L

ARM_OX <- "OXPHOS subunits"
INSTR  <- c(gsva = "OX_gsva", mitopps = "OX_mitopps")
COMPARATOR_LEV <- c("all", "site", "paired")

CONFIG_SIGN <- c(BBC3 = 1, BID = 1, BIK = 1, BAD = 1, BCL2L1 = 1, MCL1 = -1)
# E10's twelve, verbatim from E34 line 144.
PRIMING_PRO  <- c("BCL2L11", "BMF", "PMAIP1", "BBC3", "BID", "BAD", "BIK")
PRIMING_ANTI <- c("BCL2", "BCL2L1", "MCL1", "BCL2L2", "BCL2A1")
TWELVE       <- c(PRIMING_PRO, PRIMING_ANTI)

EPI_GENES  <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19",
                "KRT5", "KRT14", "KRT17", "TP63")      # E36's DECLARED nine
LUM_GENES  <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19")   # POST-HOC five
ADIPOSE_GENES    <- c("ADIPOQ", "FABP4", "PLIN1", "LEP", "CIDEC")
FIBROBLAST_GENES <- c("COL1A1", "COL1A2", "FAP", "PDGFRB", "THY1")

MYC_TX <- c("MYC", "MYCN", "MYCL")

# The one gene dropped anywhere in E36 or E38, declared in the data note s4.
DROP_FROM_44 <- "CYCS"

ADJUST <- list(raw = character(0), epithelial = "EPITHELIAL",
               epi_prolif = c("EPITHELIAL", "PROLIF"))

# --- THE ARBITRARY SIGN SPLITS ----------------------------------------------
# Rule from the data note section 3: sort alphabetically, take the negative set
# BY POSITION so it spreads across the set rather than forming a block.
# BALANCED every 2nd gene; RATIO-MATCHED every 6th. The gene lists the rule
# produces are written out in the data note and are ASSERTED below, so a change
# of sort order or set membership stops the script rather than silently
# re-drawing the split after the fact.
.neg_balanced <- function(g) { g <- sort(g); g[seq(2L, length(g), by = 2L)] }
.neg_ratio51  <- function(g) { g <- sort(g); g[seq(6L, length(g), by = 6L)] }

DECLARED_NEG <- list(
  `Fe-S cluster assembly` = list(
    balanced = c("BRIP1", "CIAPIN1", "FDX1", "FXN", "HSCB", "ISCA2", "LYRM4",
                 "NDOR1", "NUBP1", "POLD1", "SLC25A28"),
    ratio51  = c("FDX1", "ISCA2", "NUBP1")),
  mitophagy = list(
    balanced = c("ATG5", "CSNK2A1", "CSNK2B", "MAP1LC3A", "MFN1", "MTERF3",
                 "PGAM5", "PRKN", "SQSTM1", "TBK1", "TOMM22", "TOMM5",
                 "TOMM70", "UBB", "UBE2D2", "UBE2L3", "UBE2V1", "VDAC1",
                 "VDAC3"),
    ratio51  = c("CSNK2B", "MTERF3", "SQSTM1", "TOMM5", "UBE2D2", "VDAC1")))

# --- THE GATE, fixed in the data note section 7 before any value existed -----
GATE_UNSIGNED_FLOOR <- 0.729   # E36's five unsigned companions ran +0.729 to
                               # +0.902 in normals. An ALREADY-COMMITTED number.
GATE_COLLAPSE_MAX   <- 0.30    # under half that floor
GATE_RULE <- paste0(
  "LEG A passes iff rho(OXPHOS, UNSIGNED configuration) in normals is >= ",
  GATE_UNSIGNED_FLOOR, ", the floor of E36's five unsigned companions. LEG B ",
  "passes iff BOTH ratio-matched 5:1 signed comparators have |rho| < ",
  GATE_COLLAPSE_MAX, " in normals. A and B both pass -> CANCELLATION: the ",
  "signed near-zero is arithmetic, the normal comparison is closed, and ",
  "diagnostics 2 and 3 are REPORTED BUT NOT READ AS EVIDENCE. Either fails -> ",
  "NOT CANCELLATION ALONE, and 2 and 3 are worth reading - still post-hoc, ",
  "still unable to change E36's verdict.")

message("\n0. the gate, fixed at 94d5fd7 before any value existed:\n   ",
        paste(strwrap(GATE_RULE, width = 74), collapse = "\n   "))
set.seed(SEED)

# =============================================================================
# 1. The joint snapshot, and E36's composite constructor verbatim
# =============================================================================
message("\n1. reading the joint snapshot (E37)")

for (p in c(PATH_JOINT_LINEAR, PATH_JOINT_VST, PATH_JOINT_SCORES)) {
  if (!file.exists(p)) {
    stop("the joint snapshot is absent: ", basename(p),
         ". Source scripts/E37_joint_snapshot_and_gate.R, then E36, first.",
         call. = FALSE)
  }
}
jl <- readRDS(PATH_JOINT_LINEAR)
jv <- readRDS(PATH_JOINT_VST)
js <- readRDS(PATH_JOINT_SCORES)
stopifnot(identical(jl$scale, "linear_deseq2_normalised"),
          identical(jv$scale, "vst_log"),
          identical(colnames(jl$mat), colnames(jv$mat)))
smap <- jv$sample_map
sd_  <- readRDS(file.path(DIR_RESULTS, "set_definitions.rds"))
e14  <- readRDS(PATH_E14)
message("   ", nrow(smap), " samples: ", sum(smap$tissue == "tumour"),
        " tumour, ", sum(smap$tissue == "normal"), " normal")

# SCALE: log2(LINEAR + 1), z over all 1,208. E36 section 3, verbatim.
GL_ALL <- log2(jl$mat + 1)
.composite <- function(genes, signs = NULL, label) {
  miss <- setdiff(genes, rownames(GL_ALL))
  if (length(miss)) stop("composite '", label, "' missing: ",
                         paste(miss, collapse = ", "), call. = FALSE)
  G  <- GL_ALL[genes, , drop = FALSE]
  gv <- apply(G, 1L, stats::var)
  if (any(!(gv > 0))) stop("composite '", label, "' zero-variance gene",
                           call. = FALSE)
  Z <- t(scale(t(G)))
  w <- if (is.null(signs)) stats::setNames(rep(1, length(genes)), genes) else
    signs[rownames(Z)]
  as.numeric(colMeans(sweep(Z, 1L, w[rownames(Z)], "*")))
}

mito_inventory <- unique(suppressWarnings(readxl::read_excel(
  PATH_MITOCARTA, sheet = "A Human MitoCarta3.0"))$Symbol)
.halves <- function(set_name, drop = character(0)) {
  g <- setdiff(intersect(e14$comparators[[set_name]], rownames(GL_ALL)), drop)
  list(mito = intersect(g, mito_inventory), cyto = setdiff(g, mito_inventory))
}

# =============================================================================
# 2. The composites - every one the four diagnostics need
# =============================================================================
message("\n2. the composites")

# --- the compartment sets, and the fact that gates diagnostic 2 --------------
# THE CONFIGURATION HAS NO CYTOSOLIC HALF. Established in the data note and
# asserted here, so a change in MitoCarta would stop the script rather than
# silently produce a one-sided "split".
cfg_mito <- intersect(names(CONFIG_SIGN), mito_inventory)
twelve_mito <- intersect(TWELVE, mito_inventory)
CONFIG_SPLITTABLE <- length(cfg_mito) < length(CONFIG_SIGN)
if (CONFIG_SPLITTABLE) {
  stop("the configuration now HAS a cytosolic half (", length(CONFIG_SIGN) -
       length(cfg_mito), " gene(s)). The data note's section 4 is out of date ",
       "and diagnostic 2 must be redesigned, not quietly run.", call. = FALSE)
}
message("   configuration: ", length(cfg_mito), "/", length(CONFIG_SIGN),
        " in MitoCarta -> NO cytosolic half. E10's twelve: ",
        length(twelve_mito), "/12. Diagnostic 2 cannot run on either.")

AP44 <- .halves("apoptotic machinery (44)", drop = DROP_FROM_44)
ov_ox <- intersect(unlist(AP44), sd_$arm_sets[[ARM_OX]])
if (length(ov_ox)) {
  stop("after dropping ", DROP_FROM_44, " the apoptotic 44 still overlaps the ",
       "OXPHOS arm: ", paste(ov_ox, collapse = ", "), call. = FALSE)
}
FES <- .halves("Fe-S cluster assembly")
MPG <- .halves("mitophagy")
MPP <- .halves("mitophagy, PINK1/PRKN only")
message("   apoptotic machinery (44) less ", DROP_FROM_44, ": ",
        length(AP44$mito), " mitochondrial / ", length(AP44$cyto),
        " cytosolic; ", sum(names(CONFIG_SIGN) %in% unlist(AP44)),
        " of the configuration's 6 sit inside it, so its mitochondrial half",
        " is NOT a control")
message("   Fe-S ", length(FES$mito), "/", length(FES$cyto), "; mitophagy ",
        length(MPG$mito), "/", length(MPG$cyto), "; PINK1/PRKN ",
        length(MPP$mito), "/", length(MPP$cyto))

# --- the arbitrary signed comparators, ASSERTED against the data note --------
.signed_set <- function(set_name, which) {
  g   <- sort(intersect(e14$comparators[[set_name]], rownames(GL_ALL)))
  neg <- if (which == "balanced") .neg_balanced(g) else .neg_ratio51(g)
  declared <- DECLARED_NEG[[set_name]][[which]]
  if (!identical(sort(neg), sort(declared))) {
    stop("the ", which, " split of '", set_name, "' does not reproduce the ",
         "gene list declared at 94d5fd7. Computed: ",
         paste(sort(neg), collapse = ","), ". Declared: ",
         paste(sort(declared), collapse = ","),
         ". STOP - a re-drawn arbitrary split is worthless.", call. = FALSE)
  }
  stats::setNames(ifelse(g %in% neg, -1, 1), g)
}
SIGNED_COMPARATORS <- list()
for (sn in names(DECLARED_NEG)) for (wh in c("balanced", "ratio51")) {
  SIGNED_COMPARATORS[[paste0(sn, " [", wh, "]")]] <- .signed_set(sn, wh)
}
message("   arbitrary sign splits reproduce the declared gene lists exactly")
for (nm in names(SIGNED_COMPARATORS)) {
  w <- SIGNED_COMPARATORS[[nm]]
  message("      ", sprintf("%-34s %2d pos : %2d neg, mean weight %+.3f", nm,
          sum(w > 0), sum(w < 0), mean(w)))
}
message("      ", sprintf("%-34s %2d pos : %2d neg, mean weight %+.3f",
        "the configuration, for reference", sum(CONFIG_SIGN > 0),
        sum(CONFIG_SIGN < 0), mean(CONFIG_SIGN)))

# --- the sign-permutation of the configuration's own six ---------------------
PERM_SIGNS <- lapply(sort(names(CONFIG_SIGN)), function(g) {
  s <- stats::setNames(rep(1, length(CONFIG_SIGN)), names(CONFIG_SIGN))
  s[g] <- -1; s })
names(PERM_SIGNS) <- paste0("config, -", sort(names(CONFIG_SIGN)))

# --- build everything --------------------------------------------------------
CM <- tibble::tibble(
  sample_id = smap$sample_id, patient = smap$patient,
  tissue    = factor(smap$tissue, levels = c("normal", "tumour")),
  OX_gsva    = as.numeric(js$gsva[ARM_OX, smap$sample_id]),
  OX_mitopps = as.numeric(js$mitopps_arms[ARM_OX, smap$sample_id]),
  PROLIF     = as.numeric(js$gsva["PROLIF_DISJOINT", smap$sample_id]),
  EPITHELIAL = .composite(EPI_GENES, NULL, "EPITHELIAL (declared 9)"),
  LUMINAL    = .composite(LUM_GENES, NULL, "LUMINAL (post-hoc 5)"),
  ADIPOSE    = .composite(ADIPOSE_GENES, NULL, "adipose"),
  FIBROBLAST = .composite(FIBROBLAST_GENES, NULL, "fibroblast"))

# The readouts, with a role and a label that travel with every number.
RO <- tibble::tribble(
  ~readout,              ~label,                                  ~family,      ~signed,
  "config_signed",       "configuration, SIGNED (E36's readout)",  "configuration", TRUE,
  "config_unsigned",     "configuration, UNSIGNED (6)",            "configuration", FALSE,
  "twelve_unsigned",     "E10's twelve, unsigned",                 "configuration", FALSE,
  "twelve_signed",       "E10's twelve, signed 7:5",               "configuration", TRUE,
  "ap44_mito",           "apoptotic 44, mitochondrial half (19)",  "compartment",   FALSE,
  "ap44_cyto",           "apoptotic 44, cytosolic half (24)",      "compartment",   FALSE,
  "fes_mito",            "Fe-S, mitochondrial half",               "compartment",   FALSE,
  "fes_cyto",            "Fe-S, cytosolic half (E36's CONTROL)",   "compartment",   FALSE,
  "mitophagy_mito",      "mitophagy, mitochondrial half",          "compartment",   FALSE,
  "mitophagy_cyto",      "mitophagy, cytosolic half",              "compartment",   FALSE,
  "mitophagy_pp_mito",   "PINK1/PRKN, mitochondrial half",         "compartment",   FALSE,
  "mitophagy_pp_cyto",   "PINK1/PRKN, cytosolic half",             "compartment",   FALSE,
  "fes_whole",           "Fe-S assembly, whole, unsigned",         "unsigned",      FALSE,
  "mitophagy_whole",     "mitophagy, whole, unsigned",             "unsigned",      FALSE)

CM$config_signed   <- .composite(names(CONFIG_SIGN), CONFIG_SIGN, "config signed")
CM$config_unsigned <- .composite(names(CONFIG_SIGN), NULL, "config unsigned")
CM$twelve_unsigned <- .composite(TWELVE, NULL, "twelve unsigned")
CM$twelve_signed   <- .composite(
  TWELVE, stats::setNames(c(rep(1, length(PRIMING_PRO)),
                            rep(-1, length(PRIMING_ANTI))), TWELVE),
  "twelve signed")
CM$ap44_mito <- .composite(AP44$mito, NULL, "ap44 mito")
CM$ap44_cyto <- .composite(AP44$cyto, NULL, "ap44 cyto")
CM$fes_mito  <- .composite(FES$mito, NULL, "fes mito")
CM$fes_cyto  <- .composite(FES$cyto, NULL, "fes cyto")
CM$mitophagy_mito    <- .composite(MPG$mito, NULL, "mpg mito")
CM$mitophagy_cyto    <- .composite(MPG$cyto, NULL, "mpg cyto")
CM$mitophagy_pp_mito <- .composite(MPP$mito, NULL, "mpp mito")
CM$mitophagy_pp_cyto <- .composite(MPP$cyto, NULL, "mpp cyto")
CM$fes_whole       <- .composite(unlist(FES, use.names = FALSE), NULL, "fes whole")
CM$mitophagy_whole <- .composite(unlist(MPG, use.names = FALSE), NULL, "mpg whole")

for (nm in names(SIGNED_COMPARATORS)) {
  key <- paste0("sgn_", gsub("[^A-Za-z0-9]+", "_", nm))
  CM[[key]] <- .composite(names(SIGNED_COMPARATORS[[nm]]),
                          SIGNED_COMPARATORS[[nm]], nm)
  RO <- dplyr::bind_rows(RO, tibble::tibble(
    readout = key, label = paste0("ARBITRARY signed: ", nm),
    family = "arbitrary_signed", signed = TRUE))
}
for (nm in names(PERM_SIGNS)) {
  key <- paste0("perm_", gsub("[^A-Za-z0-9]+", "_", nm))
  CM[[key]] <- .composite(names(PERM_SIGNS[[nm]]), PERM_SIGNS[[nm]], nm)
  RO <- dplyr::bind_rows(RO, tibble::tibble(
    readout = key, label = nm, family = "sign_permutation", signed = TRUE))
}
for (g in MYC_TX) CM[[paste0("tx_", g)]] <- as.numeric(GL_ALL[g, ])
message("   ", nrow(RO), " readouts built; ", length(MYC_TX),
        " MYC-family transcripts carried separately")

# =============================================================================
# 3. The comparators and the patient-level resamples - E36 section 4, verbatim
# =============================================================================
message("\n3. the comparators and the patient-level resamples")

NOR_PAT   <- CM$patient[CM$tissue == "normal"]
NOR_SITES <- sort(unique(substr(NOR_PAT, 6L, 7L)))
tum_rows  <- which(CM$tissue == "tumour")
nor_rows  <- which(CM$tissue == "normal")
stopifnot(!anyDuplicated(CM$patient[tum_rows]),
          !anyDuplicated(CM$patient[nor_rows]))
T_IX <- stats::setNames(tum_rows, CM$patient[tum_rows])
N_IX <- stats::setNames(nor_rows, CM$patient[nor_rows])
POOL <- list(
  all    = names(T_IX),
  site   = names(T_IX)[substr(names(T_IX), 6L, 7L) %in% NOR_SITES],
  paired = intersect(names(T_IX), NOR_PAT))
for (cm in COMPARATOR_LEV) message("   ", sprintf("%-7s n_tumour %5d",
                                                  cm, length(POOL[[cm]])))

message("   drawing ", format(N_BOOT, big.mark = ","),
        " patient-level resamples per comparator")
BOOT <- stats::setNames(lapply(COMPARATOR_LEV, function(cm) {
  P <- POOL[[cm]]; np <- length(P)
  lapply(seq_len(N_BOOT), function(b) {
    d  <- P[sample.int(np, np, replace = TRUE)]
    nn <- N_IX[d]
    list(t = unname(T_IX[d]), n = unname(nn[!is.na(nn)]))
  })
}), COMPARATOR_LEV)

.partial_spearman <- function(x, y, Z = NULL) {
  ok <- stats::complete.cases(x, y, Z)
  if (sum(ok) < 10L) return(NA_real_)
  rx <- rank(x[ok]); ry <- rank(y[ok])
  if (is.null(Z)) return(stats::cor(rx, ry))
  RZ <- apply(as.matrix(Z)[ok, , drop = FALSE], 2L, rank)
  H  <- qr(cbind(1, RZ))
  stats::cor(qr.resid(H, rx), qr.resid(H, ry))
}
.ci <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, c((1 - CI_LEVEL) / 2, 1 - (1 - CI_LEVEL) / 2)))
}
.excl0 <- function(lo, hi) isTRUE(lo > 0) || isTRUE(hi < 0)

# One readout, one adjustment, one instrument, one comparator. Returns both
# tissues and their difference, with intervals when asked for.
.block <- function(cm, inst, readout, adj_name, do_boot, zcols = NULL) {
  ox <- INSTR[[inst]]
  zc <- if (is.null(zcols)) ADJUST[[adj_name]] else zcols
  pt <- POOL[[cm]]
  rt <- unname(T_IX[pt]); rn <- unname(N_IX[intersect(pt, names(N_IX))])
  .ps <- function(rows) .partial_spearman(
    CM[[ox]][rows], CM[[readout]][rows],
    if (length(zc)) as.matrix(CM[rows, zc, drop = FALSE]) else NULL)
  et <- .ps(rt); en <- .ps(rn)
  if (do_boot) {
    bs <- BOOT[[cm]]
    bt <- vapply(bs, function(z) .ps(z$t), numeric(1))
    bn <- vapply(bs, function(z) .ps(z$n), numeric(1))
    ct <- .ci(bt); cn <- .ci(bn); cd <- .ci(bt - bn)
  } else ct <- cn <- cd <- c(NA_real_, NA_real_)
  tibble::tibble(
    comparator = cm, instrument = inst, readout = readout,
    adjustment = adj_name,
    adjusted_for = if (length(zc)) paste(zc, collapse = " + ") else "-",
    rho_normal = en, n_lo = cn[1], n_hi = cn[2],
    rho_tumour = et, t_lo = ct[1], t_hi = ct[2],
    diff = et - en, d_lo = cd[1], d_hi = cd[2],
    normal_excludes_0 = .excl0(cn[1], cn[2]),
    tumour_excludes_0 = .excl0(ct[1], ct[2]),
    diff_excludes_0   = .excl0(cd[1], cd[2]),
    bootstrapped = do_boot)
}

# =============================================================================
# 4. DIAGNOSTIC 1 - signing. This gates the reading of 2 and 3.
# =============================================================================
message("\n4. DIAGNOSTIC 1: is the signed near-zero cancellation or biology?")

G1 <- tidyr::expand_grid(comparator = COMPARATOR_LEV,
                         instrument = names(INSTR),
                         readout = RO$readout,
                         adjustment = names(ADJUST)) %>%
  # Bootstrapped: the PAIRED comparator, UNADJUSTED. That is exactly what the
  # gate, the sign-permutation and the figure read; the adjusted columns are
  # reported as point estimates. A COST DECISION with no bearing on any
  # reading, narrowed here rather than discovered in the runtime.
  dplyr::mutate(do_boot = comparator == "paired" & adjustment == "raw")
message("   ", nrow(G1), " cells, ", sum(G1$do_boot),
        " bootstrapped (paired comparator, unadjusted - what the gate reads)")

D1 <- vector("list", nrow(G1))
for (i in seq_len(nrow(G1))) {
  g <- G1[i, ]
  D1[[i]] <- .block(g$comparator, g$instrument, g$readout, g$adjustment,
                    g$do_boot)
  if (i %% 60L == 0L) message("      ", i, " / ", nrow(G1))
}
D1 <- dplyr::bind_rows(D1) %>% dplyr::left_join(RO, by = "readout")

# --- THE GATE ----------------------------------------------------------------
# Read on the UNADJUSTED normal side, which is where E36's +0.729 to +0.902
# range was measured and so the only scale on which the floor is meaningful.
.nrho <- function(key, inst = "gsva")
  D1$rho_normal[D1$readout == key & D1$instrument == inst &
                  D1$adjustment == "raw" & D1$comparator == "paired"]

leg_a_val <- .nrho("config_unsigned")
ratio_keys <- RO$readout[RO$family == "arbitrary_signed" &
                           grepl("ratio51", RO$readout)]
leg_b_vals <- vapply(ratio_keys, .nrho, numeric(1))

LEG_A <- isTRUE(leg_a_val >= GATE_UNSIGNED_FLOOR)
LEG_B <- all(abs(leg_b_vals) < GATE_COLLAPSE_MAX)
CANCELLATION <- LEG_A && LEG_B

GATE_TEXT <- if (CANCELLATION) paste0(
  "CANCELLATION. The unsigned configuration reaches ",
  sprintf("%+.3f", leg_a_val), " in normals, at or above E36's unsigned floor ",
  "of ", GATE_UNSIGNED_FLOOR, ", AND both ratio-matched arbitrary signed ",
  "comparators collapse below ", GATE_COLLAPSE_MAX, " (",
  paste(sprintf("%+.3f", leg_b_vals), collapse = ", "), "). THE SIGNED ",
  "NEAR-ZERO IS ARITHMETIC. The normal comparison is closed, and diagnostics ",
  "2 and 3 are REPORTED BUT NOT READ AS EVIDENCE.") else paste0(
  "NOT CANCELLATION ALONE. Leg A (unsigned configuration >= ",
  GATE_UNSIGNED_FLOOR, "): ", LEG_A, ", observed ",
  sprintf("%+.3f", leg_a_val), ". Leg B (both ratio-matched signed ",
  "comparators below ", GATE_COLLAPSE_MAX, "): ", LEG_B, ", observed ",
  paste(sprintf("%+.3f", leg_b_vals), collapse = ", "),
  ". Something is not explained by signing alone, so diagnostics 2 and 3 are ",
  "worth reading - STILL POST-HOC, and still unable to change E36's verdict.")

message("\n   ", paste(strwrap(GATE_TEXT, width = 74), collapse = "\n   "))

# --- the sign-permutation, which needs no threshold --------------------------
PERM <- D1 %>%
  dplyr::filter(family == "sign_permutation", adjustment == "raw",
                comparator == "paired") %>%
  dplyr::select(instrument, label, rho_normal, n_lo, n_hi, rho_tumour)
message("\n   the sign-permutation of the configuration's own six (normals,",
        " unadjusted):")
PERM %>% dplyr::filter(instrument == "gsva") %>%
  dplyr::transmute(variant = label, normal = round(rho_normal, 3),
                   ci = paste0("[", round(n_lo, 2), ", ", round(n_hi, 2), "]"),
                   tumour = round(rho_tumour, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
PERM_ALL_COLLAPSE <- all(abs(PERM$rho_normal) < GATE_COLLAPSE_MAX)
PERM_ONLY_MCL1 <- (abs(PERM$rho_normal[grepl("MCL1", PERM$label)]) <
                     GATE_COLLAPSE_MAX) |> all() &&
  !any(abs(PERM$rho_normal[!grepl("MCL1", PERM$label)]) < GATE_COLLAPSE_MAX)
message("   all six collapse below ", GATE_COLLAPSE_MAX, ": ",
        PERM_ALL_COLLAPSE, " | only the MCL1 variant does: ", PERM_ONLY_MCL1)

# =============================================================================
# 5. DIAGNOSTIC 2 - compartment halves. POST-HOC, and gated by diagnostic 1.
# =============================================================================
# E14's claim is that the mitochondrial half rises along the OXPHOS axis while
# the cytosolic half falls. Asked directly, in both tissues, each half UNSIGNED
# so it is on the same footing as every comparator.
#
# THE CONFIGURATION ITSELF CANNOT BE DECOMPOSED - section 2 asserted it. This
# runs on E14's apoptotic machinery (44), which is the object E14's claim is
# about, with Fe-S and mitophagy as the ranks it is a rank against.
#
# LIMITATION, not fixable here: the mitochondrial half IS mitochondrial, so
# part of any correlation with OXPHOS is shared mitochondrial content. E14
# built a compartment-matched null for exactly that; E38 does not rebuild it,
# so a mitochondrial half is read only as a RANK against the other sets' halves.
message("\n5. DIAGNOSTIC 2: compartment halves (POST-HOC)")

HALF_PAIRS <- tibble::tribble(
  ~set,                   ~mito,               ~cyto,
  "apoptotic 44",         "ap44_mito",         "ap44_cyto",
  "Fe-S assembly",        "fes_mito",          "fes_cyto",
  "mitophagy",            "mitophagy_mito",    "mitophagy_cyto",
  "mitophagy PINK1/PRKN", "mitophagy_pp_mito", "mitophagy_pp_cyto")

# mito-minus-cyto within each tissue, and the difference of that difference
# between tissues, on the SAME resamples so the contrast needs no refit.
.half_block <- function(set_lab, k_mito, k_cyto, inst, adj_name) {
  ox <- INSTR[[inst]]; zc <- ADJUST[[adj_name]]
  rt <- unname(T_IX[POOL$paired]); rn <- unname(N_IX[POOL$paired])
  .ps <- function(rows, key) .partial_spearman(
    CM[[ox]][rows], CM[[key]][rows],
    if (length(zc)) as.matrix(CM[rows, zc, drop = FALSE]) else NULL)
  mt <- .ps(rt, k_mito); ct <- .ps(rt, k_cyto)
  mn <- .ps(rn, k_mito); cn <- .ps(rn, k_cyto)
  bs <- BOOT$paired
  gap_t <- vapply(bs, function(z) .ps(z$t, k_mito) - .ps(z$t, k_cyto), numeric(1))
  gap_n <- vapply(bs, function(z) .ps(z$n, k_mito) - .ps(z$n, k_cyto), numeric(1))
  did   <- gap_t - gap_n
  cgt <- .ci(gap_t); cgn <- .ci(gap_n); cdd <- .ci(did)
  tibble::tibble(
    set = set_lab, instrument = inst, adjustment = adj_name,
    rho_mito_normal = mn, rho_cyto_normal = cn, gap_normal = mn - cn,
    gn_lo = cgn[1], gn_hi = cgn[2],
    rho_mito_tumour = mt, rho_cyto_tumour = ct, gap_tumour = mt - ct,
    gt_lo = cgt[1], gt_hi = cgt[2],
    did = (mt - ct) - (mn - cn), did_lo = cdd[1], did_hi = cdd[2],
    gap_normal_excludes_0 = .excl0(cgn[1], cgn[2]),
    gap_tumour_excludes_0 = .excl0(cgt[1], cgt[2]),
    did_excludes_0        = .excl0(cdd[1], cdd[2]))
}
D2 <- dplyr::bind_rows(lapply(seq_len(nrow(HALF_PAIRS)), function(i)
  dplyr::bind_rows(lapply(names(INSTR), function(inst)
    dplyr::bind_rows(lapply(names(ADJUST), function(a)
      .half_block(HALF_PAIRS$set[i], HALF_PAIRS$mito[i], HALF_PAIRS$cyto[i],
                  inst, a)))))))
message("   ", nrow(D2), " set x instrument x adjustment cells, paired",
        " comparator, bootstrapped")
D2 %>% dplyr::filter(adjustment == "raw", instrument == "gsva") %>%
  dplyr::transmute(set, normal_gap = round(gap_normal, 3),
                   tumour_gap = round(gap_tumour, 3),
                   DiD = round(did, 3),
                   DiD_ci = paste0("[", round(did_lo, 2), ", ",
                                   round(did_hi, 2), "]"),
                   DiD_excl_0 = did_excludes_0) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 6. DIAGNOSTIC 3 - the luminal-only score. POST-HOC, SECOND ATTEMPT.
# =============================================================================
# E36's DECLARED nine-gene EPITHELIAL failed its own validity check because its
# luminal and basal halves separate the tissues in opposite directions and
# cancel. This tries the luminal half alone, under E36'S CRITERION UNCHANGED.
# It NEVER replaces the declared version, which is reported beside it.
message("\n6. DIAGNOSTIC 3: the luminal-only score (POST-HOC, second attempt)")

.auc <- function(v) {
  r <- rank(v); nt <- length(tum_rows); nn <- length(nor_rows)
  (sum(r[tum_rows]) - nt * (nt + 1) / 2) / (nt * nn)
}
D3_VALID <- tibble::tibble(
  score = c("EPITHELIAL (DECLARED 9)", "LUMINAL (POST-HOC 5)"),
  key   = c("EPITHELIAL", "LUMINAL"),
  auc_tumour_higher = vapply(c("EPITHELIAL", "LUMINAL"),
                             function(k) .auc(CM[[k]]), numeric(1)),
  rho_adipose = vapply(c("EPITHELIAL", "LUMINAL"), function(k)
    stats::cor(CM[[k]], CM$ADIPOSE, method = "spearman"), numeric(1)),
  rho_fibroblast = vapply(c("EPITHELIAL", "LUMINAL"), function(k)
    stats::cor(CM[[k]], CM$FIBROBLAST, method = "spearman"), numeric(1))) %>%
  dplyr::mutate(higher_in_tumour = auc_tumour_higher > 0.5,
                anticorrelated_with_stroma = rho_adipose < 0 &
                  rho_fibroblast < 0,
                passes = higher_in_tumour & anticorrelated_with_stroma)
D3_VALID %>%
  dplyr::transmute(score, AUC = round(auc_tumour_higher, 3),
                   rho_adipose = round(rho_adipose, 3),
                   rho_fibroblast = round(rho_fibroblast, 3),
                   higher_in_tumour, anticorrelated_with_stroma, passes) %>%
  as.data.frame() %>% print(row.names = FALSE)

LUMINAL_PASSES <- D3_VALID$passes[D3_VALID$key == "LUMINAL"]
message("   LUMINAL under E36's criterion, UNCHANGED: ",
        if (LUMINAL_PASSES) "PASSES" else "FAILS")

if (LUMINAL_PASSES) {
  message("   repeating E36's Q1 and the Fe-S control check with LUMINAL as",
          " the adjustment, on all three comparators")
  G3 <- tidyr::expand_grid(comparator = COMPARATOR_LEV,
                           instrument = names(INSTR),
                           readout = c("config_signed", "fes_cyto"),
                           adjustment = c("raw", "luminal"))
  D3 <- dplyr::bind_rows(lapply(seq_len(nrow(G3)), function(i) {
    g <- G3[i, ]
    .block(g$comparator, g$instrument, g$readout, g$adjustment, TRUE,
           zcols = if (g$adjustment == "luminal") "LUMINAL" else character(0))
  }))
  .read_cell <- function(r) {
    if (is.na(r$rho_normal) || is.na(r$diff)) return(NA_character_)
    if (!r$normal_excludes_0 && r$tumour_excludes_0 && r$diff_excludes_0)
      return("ABSENT IN NORMAL")
    if (r$normal_excludes_0 && r$tumour_excludes_0 && r$diff_excludes_0)
      return("WEAKER IN NORMAL")
    if (r$normal_excludes_0 && !r$diff_excludes_0) return("PRESENT IN NORMAL")
    if (!r$normal_excludes_0 && !r$diff_excludes_0) return("INDETERMINATE")
    "COUPLING IN NORMAL ONLY"
  }
  D3_CELL <- D3 %>% dplyr::filter(readout == "config_signed") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(reading = .read_cell(dplyr::pick(dplyr::everything()))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(comparator, instrument) %>%
    dplyr::summarise(readings = paste(sort(unique(reading)), collapse = " / "),
                     agreed = dplyr::n_distinct(reading) == 1L,
                     .groups = "drop")
  D3_STILL_SPLIT <- !all(D3_CELL$agreed)
  message("   do the six cells STILL split by adjustment, with LUMINAL in",
          " place of the declared score? ", D3_STILL_SPLIT)
  D3_CELL %>% as.data.frame() %>% print(row.names = FALSE)
} else {
  D3 <- NULL; D3_CELL <- NULL; D3_STILL_SPLIT <- NA
  message("   LUMINAL fails, so it is NOT used as an adjustment. E36's",
          " reading is not repeated with it.")
}

# =============================================================================
# 7. DIAGNOSTIC 4 - the MYC transcript, the 82, amplification, composition
# =============================================================================
# M_b__PROLIFSTRIP reverses the sign of E36's paired MYC result against
# M_b__FULL, and the whole reversal lives in the 82 stripped targets. These
# four are cheap and none depends on a regulon or a scoring choice.
message("\n7. DIAGNOSTIC 4: the MYC transcript, paired")

PAIR <- tibble::tibble(patient = POOL$paired) %>%
  dplyr::mutate(row_t = unname(T_IX[patient]), row_n = unname(N_IX[patient]))
.wilson <- function(k, n) {
  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2); p <- k / n; d <- 1 + z^2 / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / d,
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / d)
}
BOOT_PAIR <- replicate(N_BOOT, sample.int(nrow(PAIR), replace = TRUE),
                       simplify = FALSE)

# --- 7.1 MYC, MYCN, MYCL transcripts, paired. One number, no regulon. --------
D4_TX <- dplyr::bind_rows(lapply(MYC_TX, function(g) {
  d <- CM[[paste0("tx_", g)]][PAIR$row_t] - CM[[paste0("tx_", g)]][PAIR$row_n]
  k <- sum(d > 0); n <- length(d); w <- .wilson(k, n)
  bm <- vapply(BOOT_PAIR, function(ix) stats::median(d[ix]), numeric(1))
  ci <- .ci(bm)
  tibble::tibble(gene = g, n_pairs = n, median_diff = stats::median(d),
                 med_lo = ci[1], med_hi = ci[2],
                 iqr_lo = unname(stats::quantile(d, 0.25)),
                 iqr_hi = unname(stats::quantile(d, 0.75)),
                 n_tumour_higher = k, pct_tumour_higher = 100 * k / n,
                 pct_lo = 100 * w[1], pct_hi = 100 * w[2])
}))
message("   7.1 log2 expression, tumour minus that patient's own normal:")
D4_TX %>% dplyr::transmute(gene, n_pairs, median = round(median_diff, 3),
    median_ci = paste0("[", round(med_lo, 2), ", ", round(med_hi, 2), "]"),
    pct_higher = round(pct_tumour_higher, 1),
    pct_ci = paste0("[", round(pct_lo, 1), ", ", round(pct_hi, 1), "]")) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 7.2 the 82 stripped targets, in full -----------------------------------
FULL82 <- sd_$collectri_sets[["M_b__FULL"]]
KEPT   <- sd_$collectri_sets[["M_b__PROLIFSTRIP"]]
S82    <- sort(setdiff(FULL82, KEPT))
D4_STRIP <- tibble::tibble(
  n_full = length(FULL82), n_kept = length(KEPT), n_stripped = length(S82),
  in_PROLIF_REF = sum(S82 %in% sd_$strip_refs$PROLIF_REF),
  in_PROLIF_DISJOINT = sum(S82 %in% sd_$cov_sets$PROLIF_DISJOINT))
message("\n   7.2 ", length(S82), " targets stripped from ", length(FULL82),
        "; ", D4_STRIP$in_PROLIF_REF, " of them in PROLIF_REF (n = ",
        length(sd_$strip_refs$PROLIF_REF), "), the recorded strip reference")
message("      ", paste(strwrap(paste(S82, collapse = ", "), width = 72),
                        collapse = "\n      "))

# Enrichment against INDEPENDENT MYC target sets - the __FULL variants only.
# Every __PROLIFSTRIP variant was stripped with the SAME reference, so a
# comparison against one would be circular and is not made.
MYC_REF_SETS <- grep("__FULL$", names(sd_$myc_sets), value = TRUE)
D4_ENRICH <- dplyr::bind_rows(lapply(MYC_REF_SETS, function(nm) {
  g <- sd_$myc_sets[[nm]]
  a <- sum(S82 %in% g); b <- sum(KEPT %in% g)
  tibble::tibble(myc_set = sub("__FULL$", "", nm), set_n = length(g),
                 n_of_82 = a, pct_of_82 = 100 * a / length(S82),
                 n_of_kept = b, pct_of_kept = 100 * b / length(KEPT),
                 fold = (a / length(S82)) / max(b / length(KEPT), 1e-9))
})) %>% dplyr::arrange(dplyr::desc(fold))
message("\n   the 82 against independent MYC target sets (__FULL only, so the",
        " comparison is not circular):")
D4_ENRICH %>% dplyr::filter(n_of_82 >= 3) %>% dplyr::slice_head(n = 8) %>%
  dplyr::transmute(myc_set, set_n, pct_of_82 = round(pct_of_82, 1),
                   pct_of_kept = round(pct_of_kept, 1),
                   fold = round(fold, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 7.3 MYC amplification of the 113 ---------------------------------------
# MYC_amp is a per-patient GISTIC call, NOT an expression value. See the data
# note section 6.3 for why joining it is not snapshot-mixing.
co <- readRDS(PATH_COVARS)$covariates
D4_AMP <- dplyr::bind_rows(
  tibble::tibble(group = "the 113 paired",
                 MYC_amp = co$MYC_amp[match(POOL$paired, co$patient)]),
  tibble::tibble(group = "all 1,095",
                 MYC_amp = co$MYC_amp[match(POOL$all, co$patient)])) %>%
  dplyr::count(group, MYC_amp) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(pct_of_called = dplyr::if_else(
    is.na(MYC_amp), NA_real_, 100 * n / sum(n[!is.na(MYC_amp)]))) %>%
  dplyr::ungroup()
message("\n   7.3 MYC amplification:")
D4_AMP %>% dplyr::mutate(pct_of_called = round(pct_of_called, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- 7.4 composition of the 113, which must REPRODUCE data note 2 -----------
fr <- readRDS(file.path(DIR_RESULTS, "frames.rds"))$frames
fr <- fr[fr$cohort == "TCGA", ]
CMP113 <- tibble::tibble(patient = POOL$all) %>%
  dplyr::mutate(paired = patient %in% POOL$paired,
                PAM50 = as.character(fr$PAM50[match(patient, fr$sample_id)]),
                TSS = substr(patient, 6L, 7L))
D4_PAM50 <- CMP113 %>% dplyr::count(PAM50, paired) %>%
  tidyr::pivot_wider(names_from = paired, values_from = n, values_fill = 0L) %>%
  dplyr::rename(other = `FALSE`, paired = `TRUE`) %>%
  dplyr::mutate(pct_other = 100 * other / sum(other),
                pct_paired = 100 * paired / sum(paired))
m5 <- as.matrix(D4_PAM50[!is.na(D4_PAM50$PAM50), c("other", "paired")])
D4_TSS <- CMP113 %>% dplyr::count(TSS, paired) %>%
  tidyr::pivot_wider(names_from = paired, values_from = n, values_fill = 0L) %>%
  dplyr::rename(other = `FALSE`, paired = `TRUE`) %>%
  dplyr::arrange(dplyr::desc(paired))
D4_COMPOSITION <- list(
  pam50 = D4_PAM50, tss = D4_TSS,
  pam50_fisher_p = stats::fisher.test(m5, simulate.p.value = TRUE,
                                      B = 1e5)$p.value,
  tss_fisher_p = stats::fisher.test(as.matrix(D4_TSS[, c("other", "paired")]),
                                    simulate.p.value = TRUE, B = 1e5)$p.value,
  n_sites_total = nrow(D4_TSS),
  n_sites_with_normal = sum(D4_TSS$paired > 0),
  pct_from_top3 = 100 * sum(utils::head(D4_TSS$paired, 3)) / length(POOL$paired))
message("\n   7.4 composition of the 113 - MUST reproduce data note 2 s3:")
message("      PAM50 Fisher p = ", signif(D4_COMPOSITION$pam50_fisher_p, 3),
        " (data note 2: 0.76)")
message("      sites ", D4_COMPOSITION$n_sites_with_normal, " of ",
        D4_COMPOSITION$n_sites_total, ", top 3 carry ",
        round(D4_COMPOSITION$pct_from_top3, 1), "% (data note 2: 6 of 40, ",
        "87.6%); TSS Fisher p = ", signif(D4_COMPOSITION$tss_fisher_p, 3))
REPRODUCES_NOTE2 <- D4_COMPOSITION$n_sites_with_normal == 6L &&
  D4_COMPOSITION$n_sites_total == 40L &&
  abs(D4_COMPOSITION$pct_from_top3 - 87.6) < 0.1
if (!REPRODUCES_NOTE2) {
  message("      !! DOES NOT REPRODUCE data note 2. That is a DEFECT, not a ",
          "finding, and must be chased before anything here is read.")
}

# =============================================================================
# 8. The figure - one panel, which is the whole diagnostic
# =============================================================================
message("\n8. the figure")

.ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
for (d in c(DIR_FIGURES, DIR_TABLES, DIR_DOC_FIGURES)) .ensure_dir(d)
.sub <- function(x, width = 124) paste(strwrap(x, width = width),
                                       collapse = "\n")

FAM_LAB <- c(unsigned = "UNSIGNED, whole sets",
             compartment = "compartment halves (unsigned)",
             configuration = "the configuration, and E10's twelve",
             arbitrary_signed = "ARBITRARILY signed comparators",
             sign_permutation = "sign-permutation of the six")
FAM_ORDER <- names(FAM_LAB)

FIG <- D1 %>%
  dplyr::filter(comparator == "paired", adjustment == "raw") %>%
  dplyr::mutate(fam = factor(FAM_LAB[family], levels = FAM_LAB)) %>%
  dplyr::arrange(match(family, FAM_ORDER), rho_normal) %>%
  dplyr::mutate(label = factor(label, levels = unique(label)))

FIG_LONG <- FIG %>%
  tidyr::pivot_longer(c(rho_normal, rho_tumour), names_to = "tissue",
                      values_to = "rho") %>%
  dplyr::mutate(
    lo = dplyr::if_else(tissue == "rho_normal", n_lo, t_lo),
    hi = dplyr::if_else(tissue == "rho_normal", n_hi, t_hi),
    tissue = factor(dplyr::if_else(tissue == "rho_normal",
                                   "normal (n = 113)", "tumour (n = 113)"),
                    levels = c("normal (n = 113)", "tumour (n = 113)")))

pFig <- ggplot2::ggplot(FIG_LONG,
    ggplot2::aes(x = rho, y = label, colour = tissue)) +
  ggplot2::annotate("rect", xmin = -GATE_COLLAPSE_MAX, xmax = GATE_COLLAPSE_MAX,
                    ymin = -Inf, ymax = Inf, fill = "grey88", alpha = 0.55) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  ggplot2::geom_vline(xintercept = GATE_UNSIGNED_FLOOR, linewidth = 0.4,
                      linetype = "22", colour = "#2166ac") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi),
                          position = ggplot2::position_dodge(width = 0.65),
                          linewidth = 0.45, na.rm = TRUE) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.65),
                      size = 1.7, na.rm = TRUE) +
  ggplot2::facet_grid(fam ~ instrument, scales = "free_y", space = "free_y",
                      labeller = ggplot2::label_wrap_gen(22)) +
  ggplot2::scale_colour_manual(values = c("normal (n = 113)" = "#1f78b4",
                                          "tumour (n = 113)" = "#c2272d"),
                               name = NULL) +
  ggplot2::labs(
    title = paste0("E38 diagnostic: every composite against OXPHOS, in both ",
                   "tissues. GATE: ",
                   if (CANCELLATION) "CANCELLATION" else "NOT CANCELLATION ALONE"),
    subtitle = .sub(paste0(
      "POST-HOC AND EXPLORATORY; this cannot overturn E36's UNINTERPRETABLE ",
      "verdict. Partial Spearman, unadjusted, on the PAIRED comparator so ",
      "both tissues are n = 113; ", format(N_BOOT, big.mark = ","),
      " patient-level resamples. The dashed blue line is E36's unsigned floor ",
      "of ", GATE_UNSIGNED_FLOOR, " (leg A); the grey band is the collapse ",
      "window of +/-", GATE_COLLAPSE_MAX, " (leg B). Both were fixed at ",
      "94d5fd7 before any value here existed.")),
    x = "partial Spearman against the OXPHOS subunits arm", y = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
    strip.text.y = ggplot2::element_text(size = 7, angle = 0),
    strip.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 6.6),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold", size = 10),
    plot.subtitle = ggplot2::element_text(size = 7.2),
    plot.title.position = "plot",
    plot.caption = ggplot2::element_text(size = 5.8, hjust = 0,
                                         lineheight = 1.15),
    plot.caption.position = "plot") +
  ggplot2::labs(caption = paste(strwrap(paste0(
    "E38. EXPLORATORY AND POST-HOC. A DIAGNOSTIC that decomposes quantities ",
    "E36 already computed; it adds no hypothesis and CANNOT OVERTURN E36'S ",
    "UNINTERPRETABLE VERDICT. The luminal-only score and the compartment ",
    "halves are SECOND ATTEMPTS after seeing a first result and are reported ",
    "beside the originals, never instead of them. The arbitrary sign splits ",
    "were named gene by gene at 94d5fd7 before this ran. THE CONFIGURATION ",
    "HAS NO CYTOSOLIC HALF - all six of its genes are in MitoCarta 3.0 - so ",
    "the compartment panel uses E14's apoptotic machinery (44) less CYCS, ",
    "whose mitochondrial half contains 5 of the 6 configuration genes and is ",
    "therefore NOT a control. A mitochondrial half correlates with OXPHOS ",
    "partly through shared mitochondrial content, and E14's compartment-",
    "matched null is NOT rebuilt here. TCGA-BRCA only, unreplicated; the ",
    "normals are FIELD-ADJACENT tissue; the composition difference between ",
    "the tissues remains enormous whatever any adjustment does. N3; no ",
    "product term; no outcome variable."), width = 182), collapse = "\n"))

for (p in c(PATH_E38_FIG, PATH_E38_FIG_DOC)) {
  ggplot2::ggsave(p, pFig, width = 10.2, height = 10.4, dpi = 300)
}
for (p in c(PATH_E38_FIG_PDF, PATH_E38_FIG_DOC_PDF)) {
  ggplot2::ggsave(p, pFig, width = 10.2, height = 10.4,
                  device = grDevices::cairo_pdf)
}
message("   ", PATH_E38_FIG, "\n   ", PATH_E38_FIG_DOC, "  (tracked)")

# =============================================================================
# 9. Save
# =============================================================================
message("\n9. saving")

readr::write_csv(D1 %>% dplyr::filter(comparator == "paired",
                                      adjustment == "raw"), PATH_E38_CSV)

saveRDS(list(
  posture = paste0(
    "EXPLORATORY AND POST-HOC. A DIAGNOSTIC that decomposes quantities E36 ",
    "already computed. It adds no hypothesis, tests none, and CANNOT OVERTURN ",
    "E36'S UNINTERPRETABLE VERDICT. Diagnostics 2 and 3 are SECOND ATTEMPTS ",
    "after seeing a first result."),
  gate_rule = GATE_RULE, gate_text = GATE_TEXT,
  cancellation = CANCELLATION, leg_a = LEG_A, leg_b = LEG_B,
  leg_a_value = leg_a_val, leg_b_values = leg_b_vals,
  thresholds = c(unsigned_floor = GATE_UNSIGNED_FLOOR,
                 collapse_max = GATE_COLLAPSE_MAX),
  d1_signing = D1, d1_permutation = PERM,
  perm_all_collapse = PERM_ALL_COLLAPSE, perm_only_mcl1 = PERM_ONLY_MCL1,
  d2_compartment = D2, half_pairs = HALF_PAIRS,
  d3_validity = D3_VALID, d3_luminal_passes = LUMINAL_PASSES,
  d3_q1 = D3, d3_cells = D3_CELL, d3_still_split = D3_STILL_SPLIT,
  d4_transcripts = D4_TX, d4_stripped_genes = S82, d4_strip_counts = D4_STRIP,
  d4_enrichment = D4_ENRICH, d4_amplification = D4_AMP,
  d4_composition = D4_COMPOSITION, d4_reproduces_note2 = REPRODUCES_NOTE2,
  readouts = RO,
  sets = list(
    config_sign = CONFIG_SIGN, twelve = TWELVE,
    config_all_mitocarta = cfg_mito, twelve_all_mitocarta = twelve_mito,
    ap44_mito = AP44$mito, ap44_cyto = AP44$cyto, dropped_from_44 = DROP_FROM_44,
    fes_mito = FES$mito, fes_cyto = FES$cyto,
    mitophagy_mito = MPG$mito, mitophagy_cyto = MPG$cyto,
    mitophagy_pp_mito = MPP$mito, mitophagy_pp_cyto = MPP$cyto,
    signed_comparators = SIGNED_COMPARATORS, sign_permutations = PERM_SIGNS,
    epithelial_declared = EPI_GENES, luminal_posthoc = LUM_GENES,
    adipose = ADIPOSE_GENES, fibroblast = FIBROBLAST_GENES),
  spec = list(
    snapshot = "JOINT tumour+normal, built by E37. NEVER data/from_validation/",
    one_exception = paste0(
      "MYC_amp is read from data/from_validation/tcga_brca_covariates.rds. It ",
      "is a per-patient GISTIC call, not an expression value; E36 and E37 ",
      "already carry BUFFER_gistic across by the same reasoning."),
    n_boot = N_BOOT, ci_level = CI_LEVEL, seed = SEED,
    boot_scope = paste0(
      "Diagnostics 1 and 2: the paired comparator only. Diagnostic 3: all ",
      "three, because its question is defined over comparator x instrument. ",
      "Declared at 94d5fd7 and widened at a6125f7, both before the script."),
    no_blom = paste0(
      "No Blom transform anywhere. E36's Q1 path applies none - the partial ",
      "Spearman ranks internally - so E38 applies none either."),
    config_has_no_cytosolic_half = paste0(
      "All 6 configuration genes and all 12 of E10's twelve are in MitoCarta ",
      "3.0, so neither can be decomposed by compartment. Diagnostic 2 runs on ",
      "E14's apoptotic machinery (44) instead, less CYCS."),
    data_note = "docs/2026-10-01_e38_data.md at 94d5fd7, amended a6125f7",
    cannot = paste0(
      "Cannot overturn E36. Cannot make the tissues comparable - adipose AUC ",
      "0.054 and fibroblast 0.757 in E36. One cohort, unreplicated.")),
  built = Sys.time()), PATH_E38)

message("   ", PATH_E38, "\n   ", PATH_E38_CSV)
message("\nE38 done. GATE: ",
        if (CANCELLATION) "CANCELLATION" else "NOT CANCELLATION ALONE",
        " | LUMINAL: ", if (LUMINAL_PASSES) "passes" else "fails",
        "\n", strrep("=", 78))

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E38)

  cat(x$posture, "\n\n"); cat(x$gate_text, "\n")

  # DIAGNOSTIC 1. The normal side is the whole question.
  x$d1_signing %>%
    dplyr::filter(comparator == "paired", adjustment == "raw",
                  instrument == "gsva") %>%
    dplyr::arrange(dplyr::desc(rho_normal)) %>%
    dplyr::transmute(label, family,
                     normal = round(rho_normal, 3),
                     n_ci = paste0("[", round(n_lo, 2), ",", round(n_hi, 2), "]"),
                     tumour = round(rho_tumour, 3)) %>%
    as.data.frame()

  # The sign-permutation, which needs no threshold.
  x$d1_permutation %>% as.data.frame()
  x$perm_all_collapse; x$perm_only_mcl1

  # DIAGNOSTIC 2. Read ONLY if the gate did not say CANCELLATION.
  x$d2_compartment %>% dplyr::filter(adjustment == "raw") %>%
    dplyr::transmute(set, instrument, gap_normal = round(gap_normal, 3),
                     gap_tumour = round(gap_tumour, 3),
                     DiD = round(did, 3), did_excludes_0) %>%
    as.data.frame()

  # DIAGNOSTIC 3. The declared score is the first row, always.
  x$d3_validity %>% as.data.frame()
  x$d3_still_split
  if (!is.null(x$d3_cells)) x$d3_cells %>% as.data.frame()

  # DIAGNOSTIC 4. Independent of everything above.
  x$d4_transcripts %>% as.data.frame()
  cat(length(x$d4_stripped_genes), "stripped:\n",
      paste(x$d4_stripped_genes, collapse = ", "), "\n")
  x$d4_enrichment %>% dplyr::slice_head(n = 10) %>% as.data.frame()
  x$d4_amplification %>% as.data.frame()
  x$d4_reproduces_note2
  utils::str(x$d4_composition[c("pam50_fisher_p", "tss_fisher_p",
                                "n_sites_with_normal", "pct_from_top3")])
}
