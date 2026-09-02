# E14_curated_comparators.R
# =============================================================================
# THE OPEN QUESTION E11 LEFT, AND THE ONLY WAY TO CLOSE IT
# =============================================================================
# E11 P4/P4a: the 44 canonical apoptosis genes split by MitoCarta membership at
# 0.453 (TCGA) and 0.489 (SCAN-B) on the OXPHOS axis, and an expression- and
# sub-compartment-matched RANDOM draw of 20 mitochondrial plus 24 cytosolic
# genes gives 0.276 and 0.312. The observed value sits about 1.3 SD above that,
# at the 90th percentile of draws, in both cohorts and under both adjustments.
#
# That residue is the whole distance between two sentences:
#
#   "OXPHOS orders these genes."                        (composition)
#   "OXPHOS orders these genes BECAUSE they are apoptotic."  (a finding)
#
# MORE RANDOM DRAWS WILL NOT SETTLE IT, AND NOT BECAUSE OF POWER. The two
# questions are different:
#
#   THE RANDOM NULL ASKS   would ANY 20 mitochondrial plus 24 cytosolic genes
#                          of this expression profile and this compartment
#                          distribution split this way? Answer, near enough:
#                          yes.
#   THE QUESTION ASKS      does any other PROGRAMME that spans the outer
#                          membrane split this way? Never asked.
#
# The difference matters because a random draw's members share nothing but
# expression level and compartment, while a real pathway's members are
# co-regulated WITH EACH OTHER. Co-regulation within the set is exactly the
# property that could produce a split, and the random null does not have it. So
# the null used so far is, in the one way that counts, too lenient. What
# replaces it is a NAMED, CURATED, NON-APOPTOTIC comparator programme, and the
# reading is the head-to-head, not a p-value.
#
# =============================================================================
# WHAT A COMPARATOR HAS TO SATISFY, AND WHAT FAILS
# =============================================================================
# Six criteria. The statistic is a split between a set's mitochondrial and
# non-mitochondrial halves, so a comparator that lives entirely inside
# MitoCarta cannot produce one at all:
#
#   1. non-apoptotic
#   2. one curated, mechanistically coherent programme - not a category
#   3. SPANS THE BOUNDARY: a MitoCarta half AND a non-MitoCarta half, both
#      >= MIN_HALF genes present in both matrices
#   4. its mitochondrial half is outer-membrane-weighted, as the machinery's
#      is (13 MOM, 5 IMS, 2 MIM)
#   5. zero overlap with the 44
#   6. no OXPHOS subunit and no mitoribosome gene, which correlate with an
#      OXPHOS score by construction
#
# MITOCHONDRIAL PROTEIN IMPORT FAILS CRITERION 3, and it was the first
# candidate. Checked against the TCGA matrix on 2026-09-02:
# REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT is 62 genes and 62 OF 62 ARE IN
# MITOCARTA - TOM, TIM, SAM, PAM, the presequence proteases and TOMM34 with
# them. There is no cytosolic half to split against. The same is true of
# cristae formation (31 of 31) and mitochondrial calcium transport (21 of 22).
# Those three are LEVEL comparators only and are used as such in section 6.
#
# MITOPHAGY PASSES ALL SIX and is the primary comparator. REACTOME_MITOPHAGY is
# 38 genes in the TCGA matrix: 16 MitoCarta (15 MOM, 1 matrix) against 22
# cytosolic - the autophagy conjugation machinery, LC3, the ubiquitin-binding
# adaptors, ULK1, TBK1, CK2, SRC - and shares NO gene with the 44. Its shape is
# the machinery's shape: an outer-membrane-weighted mitochondrial half against
# a cytosolic effector half whose members TRANSLOCATE to the outer membrane on
# activation. That last property is the one that made a static localisation
# reading of BAX, BID and BAD unsafe, so the comparator carries the same
# hazard as the target, which is what a comparator is for.
#
# ITS ONE HAZARD, DECLARED: mitophagy is death-adjacent in the literature, and
# BNIP3, BNIP3L and BCL2L13 carry BH3 domains. None of the three is in the 44
# and none is in Reactome's PINK1/PRKN arm, which is scored separately for that
# reason. If the two mitophagy variants disagree, neither is reportable alone.
#
# THE PINK1/PRKN VARIANT IS COMPUTED AND CHECKED BUT NOT DRAWN. It is
# mitophagy minus the receptor arm, so on every figure it lands on top of
# mitophagy and costs a row of ink for no information. Section 4.1 asserts the
# two variants agree in sign and stops if they do not, so the control is live
# rather than decorative; the figures then plot `PLOT_SETS`, which is
# `SPLIT_SETS` without it. Deleting the variant outright would delete a
# declared control, which is not the same thing as simplifying a figure.
#
# =============================================================================
# THE FALSIFIER, WRITTEN BEFORE THE ANSWER IS SEEN
# =============================================================================
#   THE RESIDUE IS COMPOSITION - and E11 P4 becomes final - if mitophagy's
#   split lands at 0.40-0.55 in both cohorts, i.e. where the machinery's does.
#   The paper sentence does not change; it gets a named comparator in place of
#   2,000 synthetic draws, which is a much stronger way to write the same
#   claim.
#
#   THE RESIDUE IS APOPTOSIS-SPECIFIC - and the claim upgrades - if mitophagy
#   sits at its matched null (about 0.28-0.31) while the machinery is at
#   0.45-0.49, IN BOTH COHORTS and under both adjustments.
#
#   IT IS UNDECIDED, and must be written as a bounded percentile rather than a
#   result, if the comparators straddle the machinery or the two cohorts
#   disagree. FOUR STRUCTURAL COMPARATORS CANNOT PRODUCE A P-VALUE. The
#   deliverable is a RANKING THAT REPRODUCES ACROSS COHORTS, which is what this
#   repo counts as evidence.
#
# =============================================================================
# THE FREE STEP THAT COMES FIRST - WHICH HALF CARRIES THE RESIDUE
# =============================================================================
# Section 5, and it needs no new gene set. A split can sit above its null for
# two opposite reasons:
#
#   (a) the 20 mitochondrial genes track OXPHOS MORE than matched
#       mitochondrial genes do - apoptotic proteins at the organelle are extra
#       coupled to the OXPHOS programme; or
#   (b) the 24 cytosolic genes track OXPHOS LESS than matched cytosolic genes
#       do - the split is made by the cytosolic half being unusually flat.
#
# These are different claims with different follow-ups, and the split statistic
# cannot tell them apart. Two one-sided z values can, immediately.
#
# =============================================================================
# THE CONFOUND THAT COMES WITH THE CYTOSOLIC HALF, AND ITS FALSIFIER
# =============================================================================
# Section 5.1, declared here because it must be settled before section 4's
# verdict can be read at all. THE MACHINERY'S 24 NON-MITOCARTA GENES ARE NOT A
# NEUTRAL CYTOSOLIC SET. They are FAS, FASLG, TNF, TNFRSF1A, TNFRSF10A,
# TNFRSF10B, TRADD, FADD, CASP10, CFLAR, BIRC2, BIRC3, XIAP, NFKB1 and RELA -
# the death-receptor and NF-kB arm, which is as much an inflammatory module as
# an apoptotic one. Trap 2 says rho(OXPHOS subunits, leukocyte fraction) =
# -0.158 in TCGA, so an infiltrate-high tumour is OXPHOS-low AND
# death-receptor-high BY CONSTRUCTION.
#
# That single fact would produce a negative cytosolic half, a large split, and
# an "apoptosis-specific" verdict, with no apoptosis in it anywhere. And it
# would spare the comparators: mitophagy's cytosolic arm is ATG5, ATG12, LC3,
# SQSTM1, ULK1 and the ubiquitin conjugation machinery, which is housekeeping
# and carries no immune signal at all. THE COMPARATOR CANNOT CONTROL FOR THIS.
# Only the covariate can.
#
#   IT SURVIVES if the machinery's cytosolic half stays negative and below its
#   matched null in TCGA after purity AND leukocyte fraction are partialled
#   out on top of proliferation, AND if the cytosolic genes that are NOT in
#   the death-receptor module are negative on their own.
#
#   IT FAILS, and section 4's verdict is an infiltrate artefact that must be
#   reported as one, if the negative is carried by the death-receptor and
#   NF-kB genes and does not survive the adjustment.
#
# SCAN-B HAS NO PURITY ESTIMATE (trap 2) and it is never imputed, so this test
# is TCGA-only and the cross-cohort claim rests on the module breakdown, which
# both cohorts can carry.
#
# =============================================================================
# THE LADDER, WHICH IS THE OTHER THING THAT WAS ASKED FOR
# =============================================================================
# Section 6 answers a second and easier question: IS THE MACHINERY'S
# MITOCHONDRIAL HALF UNUSUAL AMONG MITOCHONDRIAL PROGRAMMES? For that the split
# is not needed and every MitoCarta leaf pathway qualifies. Each pathway of
# n >= LADDER_MIN gets its mean rho, mean |rho| and SD on both axes, and the
# machinery's 20 are placed on the resulting ladder as a PERCENTILE AMONG REAL
# MITOCHONDRIAL PROGRAMMES rather than as a z against a synthetic null. A
# reader can check a percentile against biology they already know; they cannot
# check a z against 2,000 draws they cannot see.
#
# It also delivers what P4a asked for and section 3.2 could not give: a
# sub-compartment ladder over hundreds of genes in curated functional groups -
# outer-membrane-weighted (fission, mitophagy, contact sites), matrix-weighted
# (Fe-S, proteases, chaperones), inner-membrane-weighted (cristae, calcium,
# the carriers) - instead of over the machinery's 2 inner-membrane genes.
#
# =============================================================================
# WHAT "MITOCHONDRIAL" MEANS HERE - CARRIED FORWARD FROM E11
# =============================================================================
# MitoCarta is a PROTEOME catalogue. Membership means a protein was detected at
# mitochondria, not that it sits there constitutively, and a transcript has no
# idea where its protein ends up. At this level of measurement MitoCarta
# membership marks MEMBERSHIP OF THE NUCLEAR-ENCODED MITOCHONDRIAL REGULON -
# transcripts co-regulated with mitochondrial biogenesis - and that is the
# reading throughout. It is also why the comparator has to be a programme:
# regulon membership is a co-regulation claim, and only a co-regulated set
# tests it.
#
# EXPLORATORY. Nothing here is pre-registered. SCALE: linear DESeq2-normalised
# at gene level, every correlation rank-based. SPECIES: human, natively; no
# ortholog mapping in either direction.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
source(here::here("functions", "correlation_engine.R"))
source(here::here("functions", "gene_matrix.R"))

message("\nE14: is the OXPHOS ordering of the apoptotic machinery anything a",
        "\n     curated non-apoptotic mitochondrial programme does not do?\n",
        strrep("=", 78))

PATH_E14      <- file.path(DIR_RESULTS, "curated_comparators.rds")
PATH_E14_CSV  <- file.path(DIR_TABLES,  "E14_comparator_splits.csv")
PATH_E14_CSV2 <- file.path(DIR_TABLES,  "E14_mitocarta_pathway_ladder.csv")

NULL_DRAWS <- 2000L
N_BINS     <- 20L
MIN_HALF   <- 8L     # criterion 3: a half smaller than this cannot be split on
LADDER_MIN <- 12L    # section 6: pathways below this are too few to rank
PROLIF_REF_COV <- "PROLIF_DISJOINT"

# =============================================================================
# 1. Inputs - the same objects E11 reads, built the same way
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

CANON    <- sort(e08$canonical$gene)
# E08's Reactome module per gene, used only to ask WHICH cytosolic genes carry
# the negative. Taken as it stands; nothing here re-derives it.
CANON_MODULE <- stats::setNames(e08$canonical$module, e08$canonical$gene)
stopifnot(!anyNA(CANON_MODULE[CANON]))
MITORIBO <- sd_$arm_sets[["Mitochondrial ribosome"]]
OXARM    <- sd_$arm_sets[["OXPHOS subunits"]]

.axis_mat <- function(gsva_new, arms_obj, ids) {
  m <- rbind(
    MYC    = as.numeric(gsva_new[MYC_REF, ids]),
    OXPHOS = as.numeric(arms_obj$gsva_arms["OXPHOS subunits", ids]))
  colnames(m) <- ids; m
}
COH <- list(
  TCGA = list(L = LT, res = RES_T, ids = ID_T,
              ax = .axis_mat(nw$tcga_gsva_new, mito, ID_T),
              cov = t(mito$gsva_cov[, ID_T, drop = FALSE])),
  `SCAN-B` = list(L = LS, res = RES_S, ids = ID_S,
                  ax = .axis_mat(sc$gsva_new, sc, ID_S),
                  cov = t(sc$gsva_cov[, ID_S, drop = FALSE])))
for (coh in names(COH)) {
  stopifnot(identical(rownames(COH[[coh]]$cov), COH[[coh]]$ids),
            PROLIF_REF_COV %in% colnames(COH[[coh]]$cov))
}
message("   axes and covariates built for both cohorts")

# --- MitoCarta: the inventory, the compartments and the pathway sheet --------
mitocarta_sheet <- suppressWarnings(readxl::read_xls(PATH_MITOCARTA, sheet = 2))
stopifnot(all(c("Symbol", "MitoCarta3.0_SubMitoLocalization") %in%
                names(mitocarta_sheet)))
MC_ALL  <- unique(mitocarta_sheet$Symbol)
SUBMITO <- stats::setNames(
  mitocarta_sheet[["MitoCarta3.0_SubMitoLocalization"]], mitocarta_sheet$Symbol)
# The strict pool, as in E11 section 3.1: no OXPHOS subunit, no mitoribosome
# gene, because both correlate with an OXPHOS score by construction.
MC_STRICT <- setdiff(MC_ALL, union(OXARM, MITORIBO))
message("   MitoCarta 3.0: ", length(MC_ALL), " symbols, strict pool ",
        length(MC_STRICT), " after removing the OXPHOS arm and the ",
        "mitoribosome")

mito_paths_sheet <- suppressWarnings(
  readxl::read_xls(PATH_MITOCARTA, sheet = "C MitoPathways"))
names(mito_paths_sheet) <- make.names(names(mito_paths_sheet))
mito_paths_sheet <- mito_paths_sheet[!is.na(mito_paths_sheet$MitoPathway), ]
MITO_PATHS <- stats::setNames(
  lapply(strsplit(mito_paths_sheet$Genes, ",\\s*"),
         function(z) unique(trimws(z[z != ""]))),
  mito_paths_sheet$MitoPathway)
MITO_PATH_HIER <- stats::setNames(mito_paths_sheet$MitoPathways.Hierarchy,
                                  mito_paths_sheet$MitoPathway)
message("   MitoCarta pathway sheet: ", length(MITO_PATHS), " pathways")

# =============================================================================
# 2. The comparator sets
# =============================================================================
# Reactome supplies the programmes and MitoCarta supplies the membership label,
# so the set definition and the thing being tested come from different
# catalogues. That independence is the same one E10 used and it is the reason
# the split is not circular: nothing about Reactome's mitophagy pathway knows
# which of its genes MitoCarta lists.
message("\n2. comparator sets from Reactome, labelled by MitoCarta")

reactome <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                             subcollection = "CP:REACTOME")
MSIGDB_VERSION <- unique(reactome$db_version)
message("   MSigDB ", paste(MSIGDB_VERSION, collapse = "/"), ", ",
        dplyr::n_distinct(reactome$gs_name), " Reactome sets")
RGS <- split(reactome$gene_symbol, reactome$gs_name)

WANT_REACTOME <- c(
  "REACTOME_MITOPHAGY",
  "REACTOME_PINK1_PRKN_MEDIATED_MITOPHAGY",
  "REACTOME_MITOCHONDRIAL_BIOGENESIS",
  "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
  "REACTOME_CRISTAE_FORMATION",
  "REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT",
  "REACTOME_MITOCHONDRIAL_IRON_SULFUR_CLUSTER_BIOGENESIS",
  "REACTOME_CYTOSOLIC_IRON_SULFUR_CLUSTER_ASSEMBLY")
missing_sets <- setdiff(WANT_REACTOME, names(RGS))
if (length(missing_sets)) {
  stop("these Reactome sets are not in the installed MSigDB (",
       paste(MSIGDB_VERSION, collapse = "/"), "): ",
       paste(missing_sets, collapse = ", "),
       ". The set names moved between releases - do not substitute a ",
       "near-match silently, look them up.", call. = FALSE)
}

# Fe-S cluster assembly is ONE pathway that Reactome stores as two, split at
# exactly the boundary this script is testing: the ISC arm builds clusters in
# the matrix, the CIA arm builds them in the cytosol from a sulfur species the
# matrix exports. Joining them gives a genuinely two-compartment programme with
# no signalling and no translocation - the cleanest possible comparator, and
# matrix-weighted rather than outer-membrane-weighted, so it tests the
# compartment gradient rather than matching it.
FES <- union(RGS[["REACTOME_MITOCHONDRIAL_IRON_SULFUR_CLUSTER_BIOGENESIS"]],
             RGS[["REACTOME_CYTOSOLIC_IRON_SULFUR_CLUSTER_ASSEMBLY"]])

# THE ONE SET IN THIS SCRIPT NOT TAKEN FROM A PINNED CATALOGUE, AND IT IS
# MARKED AS SUCH EVERYWHERE. Cytosolic/mitochondrial isozyme pairs: two genes
# catalysing the SAME reaction in the two compartments. Function is held
# exactly constant and only compartment varies, so this is the strongest
# logical test of "regulon membership, not function" available - and it is a
# CEILING, not a null. Several members are TCA enzymes and the mitochondrial
# aminoacyl-tRNA synthetases co-regulate with the mitoribosome, so a large
# split here is expected and is the point of reference, not a finding.
# Criterion: both members annotated, distinct genes, same catalytic activity.
PARALOGUE_PAIRS <- rbind(
  c("SHMT1", "SHMT2"),   c("MTHFD1", "MTHFD2"), c("ALDH1L1", "ALDH1L2"),
  c("GOT1",  "GOT2"),    c("MDH1",   "MDH2"),   c("IDH1",    "IDH2"),
  c("ACO1",  "ACO2"),    c("SOD1",   "SOD2"),   c("TXN",     "TXN2"),
  c("TXNRD1", "TXNRD2"), c("PRDX1",  "PRDX3"),  c("HSPA8",   "HSPA9"),
  c("PCK1",  "PCK2"),    c("ACSS2",  "ACSS1"),  c("AK1",     "AK2"),
  c("AARS1", "AARS2"),   c("CARS1",  "CARS2"),  c("DARS1",   "DARS2"),
  c("HARS1", "HARS2"),   c("IARS1",  "IARS2"),  c("LARS1",   "LARS2"),
  c("MARS1", "MARS2"),   c("NARS1",  "NARS2"),  c("RARS1",   "RARS2"),
  c("SARS1", "SARS2"),   c("TARS1",  "TARS2"),  c("VARS1",   "VARS2"),
  c("WARS1", "WARS2"),   c("YARS1",  "YARS2"),  c("FARSB",   "FARS2"))
colnames(PARALOGUE_PAIRS) <- c("cytosolic", "mitochondrial")

COMPARATORS <- list(
  `apoptotic machinery (44)`   = CANON,
  `mitophagy`                  = RGS[["REACTOME_MITOPHAGY"]],
  `mitophagy, PINK1/PRKN only` = RGS[["REACTOME_PINK1_PRKN_MEDIATED_MITOPHAGY"]],
  `Fe-S cluster assembly`      = FES,
  `mitochondrial biogenesis`   = RGS[["REACTOME_MITOCHONDRIAL_BIOGENESIS"]],
  `isozyme pairs (CEILING)`    = as.vector(PARALOGUE_PAIRS),
  `protein import`             = RGS[["REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT"]],
  `cristae formation`          = RGS[["REACTOME_CRISTAE_FORMATION"]],
  `mitochondrial Ca transport` =
    RGS[["REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT"]])
COMPARATOR_SOURCE <- c(
  `apoptotic machinery (44)`   = "E08 canonical, curated",
  `mitophagy`                  = "MSigDB C2:CP:REACTOME",
  `mitophagy, PINK1/PRKN only` = "MSigDB C2:CP:REACTOME",
  `Fe-S cluster assembly`      = "MSigDB C2:CP:REACTOME, ISC + CIA joined",
  `mitochondrial biogenesis`   = "MSigDB C2:CP:REACTOME",
  `isozyme pairs (CEILING)`    = "AUTHOR-CURATED IN THIS SCRIPT, see header",
  `protein import`             = "MSigDB C2:CP:REACTOME",
  `cristae formation`          = "MSigDB C2:CP:REACTOME",
  `mitochondrial Ca transport` = "MSigDB C2:CP:REACTOME")

# --- the row index of a set, WITH ITS GENE NAMES STILL ATTACHED -------------
# `.gene_rows` names its matrix by the ORIGINAL symbol, not by the matrix row
# it resolved to. Taking match(rownames(mat), rownames(L)) after that returns
# NA for every gene SCAN-B's symbol map renamed, and dropping those NAs makes
# the index vector SHORTER than the compartment vector computed beside it. R
# then RECYCLES the labels in split() with a warning and no error, and the
# compartments come out attached to the wrong genes. This helper resolves and
# indexes in one step so the two vectors cannot come apart, and the callers
# assert their lengths anyway.
.row_index <- function(genes, C) {
  g <- unique(genes)
  h <- vapply(g, function(x) {
    r <- C$res(x); if (length(r) != 1L) NA_character_ else r }, character(1))
  i <- match(h, rownames(C$L))
  ok <- !is.na(i)
  list(gene = g[ok], row = unname(i[ok]), missing = g[!ok])
}

# --- the composition of a set, with the <NA>-name trap closed ----------------
# Subsetting a named vector by a name it does not carry returns an element
# named <NA>. Left in place, split() drops the group silently, the label vector
# comes out constant and cor() returns NA with a warning and no error. E11 trap
# 4; unname() then names()<- is the fix.
.composition <- function(genes) {
  v <- unname(SUBMITO[genes])
  v[is.na(v)] <- "(not in MitoCarta)"
  names(v) <- genes
  stopifnot(!anyNA(v), !anyNA(names(v)), length(v) == length(genes))
  v
}

# --- EVERY COMPARATOR IS RESTRICTED TO GENES PRESENT IN BOTH COHORTS --------
# Trap 7 with a new face. SCAN-B's symbols are a 2014 UCSC build and the
# symbol map was built over the arm sets, the MitoCarta pathways and the
# estimators - not over sets invented later in this script. Eleven cytosolic
# aminoacyl-tRNA synthetases are `AARS`, `DARS`, `LARS` there and `AARS1`,
# `DARS1`, `LARS1` here, so the isozyme set arrived as 60 genes in TCGA and 49
# in SCAN-B, with the loss falling ENTIRELY on its cytosolic half. A
# comparator whose composition depends on the cohort cannot be compared across
# cohorts, and the difference would have read as biology.
#
# So each set is cut to its intersection across both matrices before anything
# is measured, and what was cut is printed. This is not the symbol map's job -
# the map is a pinned artefact of E02 and is not extended here.
.resolvable <- function(genes, C) .row_index(genes, C)$gene
drops <- dplyr::bind_rows(lapply(names(COMPARATORS), function(nm) {
  g <- unique(COMPARATORS[[nm]])
  keepg <- Reduce(intersect, lapply(COH, function(C) .resolvable(g, C)))
  lost <- setdiff(g, keepg)
  tibble::tibble(set = nm, n_declared = length(g), n_kept = length(keepg),
                 n_dropped = length(lost),
                 dropped = paste(sort(lost), collapse = ", "))
}))
COMPARATORS <- lapply(COMPARATORS, function(g) {
  g <- unique(g)
  Reduce(intersect, lapply(COH, function(C) .resolvable(g, C)))
})
message("\n   restricted to genes present in BOTH cohorts:")
drops %>% dplyr::select(set, n_declared, n_kept, n_dropped) %>%
  as.data.frame() %>% print(row.names = FALSE)
for (i in which(drops$n_dropped > 0L)) {
  message("     ", drops$set[i], " lost ", drops$n_dropped[i], ": ",
          drops$dropped[i])
}

# --- the audit, which is what decides eligibility ---------------------------
# Every criterion is a printed column, not an assumption. A comparator that
# fails 3 or 5 is dropped from section 4 by the filter below and kept for
# section 6, where only the level is read.
comparator_audit <- dplyr::bind_rows(lapply(names(COMPARATORS), function(nm) {
  dplyr::bind_rows(lapply(names(COH), function(coh) {
    C <- COH[[coh]]
    ri <- .row_index(COMPARATORS[[nm]], C)
    g  <- ri$gene
    cmp <- .composition(g)
    inmc <- cmp != "(not in MitoCarta)"
    tibble::tibble(
      set = nm, cohort = coh, source = unname(COMPARATOR_SOURCE[[nm]]),
      n_declared = length(unique(COMPARATORS[[nm]])), n_in_matrix = length(g),
      n_mito = sum(inmc), n_cytosolic = sum(!inmc),
      MOM = sum(cmp == "MOM"), IMS = sum(cmp == "IMS"),
      MIM = sum(cmp == "MIM"), Matrix = sum(cmp == "Matrix"),
      overlap_44 = length(intersect(g, CANON)),
      in_OXPHOS_arm = length(intersect(g, OXARM)),
      in_mitoribosome = length(intersect(g, MITORIBO)))
  }))
}))
message("\n   the comparator audit - every criterion is a column:")
comparator_audit %>%
  dplyr::select(set, cohort, n_in_matrix, n_mito, n_cytosolic, MOM, IMS, MIM,
                Matrix, overlap_44, in_OXPHOS_arm, in_mitoribosome) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Criterion 3 and criterion 5, applied. The machinery itself is always kept -
# it is the target, not a comparator.
elig <- comparator_audit %>%
  dplyr::group_by(set) %>%
  dplyr::summarise(min_mito = min(n_mito), min_cyt = min(n_cytosolic),
                   max_overlap = max(overlap_44),
                   max_ceiling = max(in_OXPHOS_arm + in_mitoribosome),
                   .groups = "drop") %>%
  dplyr::mutate(
    splittable = min_mito >= MIN_HALF & min_cyt >= MIN_HALF,
    clean = max_overlap == 0L | set == "apoptotic machinery (44)",
    no_ceiling_genes = max_ceiling == 0L | set == "apoptotic machinery (44)",
    eligible = splittable & clean & no_ceiling_genes)
message("\n   eligibility for the head-to-head (criterion 3: both halves >= ",
        MIN_HALF, "; criterion 5: no overlap with the 44; criterion 6: no",
        "\n   OXPHOS-subunit or mitoribosome gene, which correlate with the",
        " axis by construction):")
elig %>% as.data.frame() %>% print(row.names = FALSE)
# The 44 are exempt from 5 and 6 because they ARE the target, not a
# comparator - and they carry exactly one OXPHOS-arm gene, CYCS, which E11
# already handles by keeping it in the binning and out of the selection.
# Ordered deliberately, not alphabetically: the target first, the comparators
# in the middle, the declared ceiling last, so every table and every figure
# reads top-to-bottom in the order the argument is made.
SET_ORDER <- names(COMPARATORS)
SPLIT_SETS <- SET_ORDER[SET_ORDER %in% elig$set[elig$eligible]]
LEVEL_ONLY <- setdiff(names(COMPARATORS), SPLIT_SETS)
# Measured in full, plotted in part. See the header: the PINK1/PRKN variant is
# a mutual control on full mitophagy and lands on top of it in every panel.
NOT_PLOTTED <- "mitophagy, PINK1/PRKN only"
PLOT_SETS   <- setdiff(SPLIT_SETS, NOT_PLOTTED)
message("   head-to-head sets: ", paste(SPLIT_SETS, collapse = ", "))
message("   measured but not plotted (a control, see 4.1): ",
        paste(intersect(NOT_PLOTTED, SPLIT_SETS), collapse = ", "))
message("   level-only (no cytosolic half to split against): ",
        paste(LEVEL_ONLY, collapse = ", "))
if (!"mitophagy" %in% SPLIT_SETS) {
  stop("mitophagy is the primary comparator and it did not clear the ",
       "eligibility filter. Read the audit above before going further - ",
       "either the Reactome set moved or the matrix did.", call. = FALSE)
}

# =============================================================================
# 3. Per-gene correlations over the whole matrix
# =============================================================================
# ONE FULL-MATRIX PASS PER (cohort, axis, adjustment). Every set, every half and
# every null draw below is then a subset of the SAME vector, so a set and its
# comparator are guaranteed to have been computed identically. Building them
# per set is how a set and its comparator quietly stop being comparable.
message("\n3. per-gene correlations over the whole matrix")

ADJUSTMENTS <- list(raw = NULL, `adj. PROLIF_DISJOINT` = PROLIF_REF_COV)

per_gene <- list()
for (coh in names(COH)) {
  C <- COH[[coh]]
  for (ax in rownames(C$ax)) {
    for (adj in names(ADJUSTMENTS)) {
      cv <- ADJUSTMENTS[[adj]]
      Z  <- if (is.null(cv)) NULL else C$cov[, cv, drop = FALSE]
      per_gene[[paste(coh, ax, adj, sep = "|")]] <-
        .per_gene_rho(C$L, C$ax[ax, ], cov = Z)
    }
  }
  message("   ", coh, " - 2 axes x ", length(ADJUSTMENTS), " adjustments done")
}

# =============================================================================
# 4. THE HEAD-TO-HEAD. Every eligible programme, the same statistic
# =============================================================================
# Two statistics, both carried, for the reason recorded in E11: `split` is the
# statistic every number written up before 2026-09-02 was computed on and it
# stays the anchor; `med_diff` is in the units of the correlations themselves
# and is what the figures plot. Choosing between them after seeing their z
# values would be statistic-shopping, so no choice is made.
#
# EACH SET GETS ITS OWN SUB-COMPARTMENT-MATCHED NULL, matched to ITS OWN
# composition, so its z is on the same footing as the machinery's. The
# head-to-head of the OBSERVED values is the primary reading; the z values say
# how much of each observed value its own composition already accounts for.
message("\n4. the head-to-head: the split, per programme, against its own null")
set.seed(PROJECT_SEED)

.z <- function(obs, nd) (obs - mean(nd)) / stats::sd(nd)
.nbins <- function(n) max(2L, min(N_BINS, n %/% 25L))

# One composition-matched draw set for one set in one cohort. It does not
# depend on which correlation is read, so it is built once and reused across
# axes and adjustments.
.compartment_draws <- function(C, keep, own, want) {
  parts <- list()
  for (cmp in names(want)) {
    i_want <- want[[cmp]]
    if (!length(i_want)) next
    pool <- if (cmp == "(not in MitoCarta)") {
      setdiff(keep, which(rownames(C$L) %in% MC_ALL))
    } else {
      intersect(keep, which(rownames(C$L) %in%
                              intersect(MC_STRICT,
                                        names(SUBMITO)[!is.na(SUBMITO) &
                                                         SUBMITO == cmp])))
    }
    # The set's own genes must be IN the binning - otherwise their bin is NA and
    # the draw silently returns nothing - and OUT of the selection, which is
    # what `exclude` does at draw time. CYCS is the gene that forces this: IMS,
    # but in the OXPHOS arm, so the strict pool excludes it.
    pool <- union(pool, i_want)
    parts[[cmp]] <- list(want = i_want,
                         B = .expression_bins(C$L, pool,
                                              n_bins = .nbins(length(pool))))
  }
  parts
}

.split_stats <- function(v, idx, lab, draws, ...) {
  .md <- function(w) stats::median(w[lab == 1]) - stats::median(w[lab == 0])
  o_s <- stats::cor(v[idx], lab, method = "spearman")
  o_m <- .md(v[idx])
  n_s <- vapply(draws, function(d) stats::cor(v[d], lab, method = "spearman"),
                numeric(1))
  n_m <- vapply(draws, function(d) .md(v[d]), numeric(1))
  tibble::tibble(..., split = o_s, null_mean = mean(n_s),
                 null_sd = stats::sd(n_s), z = .z(o_s, n_s),
                 pct_of_draws_below = mean(n_s < o_s),
                 med_diff = o_m, null_mean_md = mean(n_m),
                 null_sd_md = stats::sd(n_m), z_md = .z(o_m, n_m))
}

comparator_splits <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))

  dplyr::bind_rows(lapply(SPLIT_SETS, function(nm) {
    ri  <- .row_index(COMPARATORS[[nm]], C)
    own <- ri$row
    cmp <- .composition(ri$gene)
    stopifnot(length(own) == length(cmp))
    want <- lapply(split(own, cmp), function(v) intersect(v, keep))
    want <- want[lengths(want) > 0L]

    parts <- .compartment_draws(C, keep, own, want)
    dr <- replicate(NULL_DRAWS,
                    unlist(lapply(parts, function(pp)
                      .matched_draw(pp$want, pp$B, exclude = own)),
                      use.names = FALSE),
                    simplify = FALSE)
    obs_idx <- unlist(lapply(parts, function(pp) pp$want), use.names = FALSE)
    lab <- as.numeric(rep(names(parts) != "(not in MitoCarta)",
                          vapply(parts, function(pp) length(pp$want),
                                 integer(1))))
    # A constant label gives NA with a warning and no error. Stop instead.
    if (length(unique(lab)) < 2L || length(obs_idx) != length(lab)) {
      stop("the composition labels are degenerate for ", nm, " in ", coh,
           " - ", length(obs_idx), " genes, ", length(unique(lab)),
           " label value(s).", call. = FALSE)
    }
    dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
      dplyr::bind_rows(lapply(names(ADJUSTMENTS), function(adj)
        .split_stats(per_gene[[paste(coh, ax, adj, sep = "|")]],
                     obs_idx, lab, dr, cohort = coh, set = nm, axis = ax,
                     adjustment = adj, n_mito = sum(lab == 1),
                     n_cytosolic = sum(lab == 0))))))
  }))
}))

message("\n   THE ANSWER, adjusted, OXPHOS axis - read the `split` column ",
        "ACROSS SETS first:")
comparator_splits %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::select(cohort, set, n_mito, n_cytosolic, split, null_mean, null_sd, z,
                pct_of_draws_below, med_diff, z_md) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(cohort, dplyr::desc(split)) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the same on the MYC axis:")
comparator_splits %>%
  dplyr::filter(axis == "MYC", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::select(cohort, set, split, null_mean, z, med_diff, z_md) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(cohort, dplyr::desc(split)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- the falsifier, evaluated in code rather than by eye --------------------
# The header committed to three outcomes before the numbers existed. This block
# reads them off, so the verdict is not a matter of which row the eye lands on.
verdict_tbl <- comparator_splits %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT",
                set %in% c("apoptotic machinery (44)", "mitophagy")) %>%
  dplyr::select(cohort, set, split) %>%
  tidyr::pivot_wider(names_from = set, values_from = split) %>%
  dplyr::rename(machinery = `apoptotic machinery (44)`) %>%
  dplyr::mutate(gap = machinery - mitophagy)
VERDICT <- with(verdict_tbl, {
  if (all(mitophagy >= 0.40 & mitophagy <= 0.55)) {
    "COMPOSITION: mitophagy splits where the machinery does. E11 P4 is final."
  } else if (all(gap > 0.12)) {
    paste("APOPTOSIS-SPECIFIC in both cohorts: the machinery out-splits",
          "mitophagy by", paste(round(gap, 3), collapse = " and "))
  } else if (all(gap < -0.12)) {
    "REVERSED: mitophagy out-splits the machinery in both cohorts."
  } else {
    "UNDECIDED: the cohorts disagree or the gap is inside sampling noise."
  }
})
message("\n   VERDICT AGAINST THE FALSIFIER WRITTEN IN THE HEADER:\n   ",
        VERDICT)
verdict_tbl %>% dplyr::mutate(dplyr::across(where(is.numeric),
                                            ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

# =============================================================================
# 4.1 THE MUTUAL CONTROL ON MITOPHAGY, so that not drawing it is not hiding it
# =============================================================================
# The header made full mitophagy reportable only if its PINK1/PRKN subset
# agrees with it - the subset drops BNIP3, BNIP3L and the other BH3-domain
# receptors, which are the reason mitophagy is death-adjacent at all.
#
# THE CRITERION IS SIGN AGREEMENT ON `split`, AND THE HONEST VERSION OF WHY:
# `split` is this study's anchor statistic, fixed in E11 before E14 existed
# and for reasons that had nothing to do with mitophagy. It is NOT chosen here
# because it is the one that passes - and it needs saying, because the OTHER
# statistic does not pass. `med_diff` is +0.184 for full mitophagy and -0.030
# for the subset in TCGA: opposite signs. The subset has 14 mitochondrial
# genes against 16, so a median moves easily, but that is an explanation and
# not a defence. The check below prints both, asserts on the anchor, and says
# out loud when the second one disagrees, so that "the variants agree" is
# never read as more than it is.
message("\n4.1 do the two mitophagy variants agree? (the control for not ",
        "plotting one)")
mito_pair <- comparator_splits %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT",
                set %in% c("mitophagy", NOT_PLOTTED)) %>%
  dplyr::select(cohort, set, split, med_diff, z)
mito_pair %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
# Defined before the block, or the save below fails when the variant is not
# eligible at all - which is the case this script must survive, not assume.
MITO_PAIR_MD_AGREE <- NA
if (NOT_PLOTTED %in% SPLIT_SETS) {
  agree <- mito_pair %>%
    dplyr::select(cohort, set, split) %>%
    tidyr::pivot_wider(names_from = set, values_from = split)
  same_sign <- sign(agree[["mitophagy"]]) == sign(agree[[NOT_PLOTTED]])
  if (!all(same_sign)) {
    stop("the two mitophagy variants disagree in sign on the split. The ",
         "header made each reportable only with the other, so neither can be ",
         "read and the PINK1/PRKN row must go back into the figures.",
         call. = FALSE)
  }
  message("   they agree in sign on `split` in both cohorts, so full ",
          "mitophagy is reportable and the\n   subset is left off every ",
          "figure as duplicated ink rather than as a hidden result.")

  # The second statistic, reported whether or not it agrees. Silence here
  # would make the removal look better supported than it is.
  agree_md <- mito_pair %>%
    dplyr::select(cohort, set, med_diff) %>%
    tidyr::pivot_wider(names_from = set, values_from = med_diff)
  md_same <- sign(agree_md[["mitophagy"]]) == sign(agree_md[[NOT_PLOTTED]])
  MITO_PAIR_MD_AGREE <- all(md_same)
  if (!MITO_PAIR_MD_AGREE) {
    message("   BUT `med_diff` DISAGREES IN SIGN in ",
            paste(agree_md$cohort[!md_same], collapse = ", "), ": ",
            paste(round(agree_md[["mitophagy"]][!md_same], 3), "vs",
                  round(agree_md[[NOT_PLOTTED]][!md_same], 3),
                  collapse = "; "), ".")
    message("   The subset has 14 mitochondrial genes against 16, so a ",
            "median moves easily - but\n   that is an explanation, not a ",
            "defence. It means the mitophagy comparator is\n   SOFTER than ",
            "the figures alone suggest, and any sentence resting on it must ",
            "say so.")
  }
}

# =============================================================================
# 5. WHICH HALF CARRIES THE RESIDUE
# =============================================================================
# A split above its null can come from the mitochondrial half being high or the
# cytosolic half being low, and the split statistic cannot say which. Each half
# is therefore compared with its OWN matched null, one-sided, using mean rho -
# signed, because the direction is the question. No new gene set is needed and
# the answer decides which comparator matters if section 4 comes out undecided.
message("\n5. which half carries the residue - each half against its own null")
set.seed(PROJECT_SEED + 1L)

half_tests <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))
  dplyr::bind_rows(lapply(SPLIT_SETS, function(nm) {
    ri  <- .row_index(COMPARATORS[[nm]], C)
    own <- ri$row
    cmp <- .composition(ri$gene)
    stopifnot(length(own) == length(cmp))
    want <- lapply(split(own, cmp), function(v) intersect(v, keep))
    want <- want[lengths(want) > 0L]
    parts <- .compartment_draws(C, keep, own, want)

    halves <- list(
      mitochondrial = names(parts)[names(parts) != "(not in MitoCarta)"],
      cytosolic     = intersect("(not in MitoCarta)", names(parts)))
    dplyr::bind_rows(lapply(names(halves), function(hn) {
      pp <- parts[halves[[hn]]]
      if (!length(pp)) return(NULL)
      idx <- unlist(lapply(pp, function(q) q$want), use.names = FALSE)
      dr <- replicate(NULL_DRAWS,
                      unlist(lapply(pp, function(q)
                        .matched_draw(q$want, q$B, exclude = own)),
                        use.names = FALSE),
                      simplify = FALSE)
      dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
        dplyr::bind_rows(lapply(names(ADJUSTMENTS), function(adj) {
          v <- per_gene[[paste(coh, ax, adj, sep = "|")]]
          nd <- vapply(dr, function(d) mean(v[d]), numeric(1))
          tibble::tibble(cohort = coh, set = nm, half = hn, axis = ax,
                         adjustment = adj, n = length(idx),
                         mean_rho = mean(v[idx]), null_mean = mean(nd),
                         null_sd = stats::sd(nd),
                         z = .z(mean(v[idx]), nd))
        }))))
    }))
  }))
}))
message("\n   each half against its own matched null, OXPHOS axis, adjusted:")
half_tests %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::select(cohort, set, half, n, mean_rho, null_mean, null_sd, z) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(set, half, cohort) %>% as.data.frame() %>%
  print(row.names = FALSE)
message("\n   READ THE SIGN OF z, NOT ITS SIZE. A positive z on the",
        " mitochondrial half\n   and a z near zero on the cytosolic half says",
        " the mitochondrial genes are\n   extra-coupled. The reverse says the",
        " split is made by the cytosolic half\n   being unusually flat. They",
        " are different findings with different follow-ups.")

# =============================================================================
# 5.1 THE CONFOUND: is the cytosolic half immune infiltrate?
# =============================================================================
# The falsifier declared in the header. TCGA only - SCAN-B has no purity
# estimate and it is never imputed - so the cross-cohort half of the argument
# is the module breakdown below, which both cohorts carry.
message("\n5.1 is the cytosolic half immune infiltrate? (TCGA only, trap 2)")
set.seed(PROJECT_SEED + 2L)

PL <- frames %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::select(sample_id, purity, leuko)
PL <- PL[match(ID_T, PL$sample_id), ]
PURITY_OK <- stats::complete.cases(PL[, c("purity", "leuko")])
message("   TCGA samples with purity AND leukocyte fraction: ", sum(PURITY_OK),
        " of ", length(ID_T))

CP <- list(L = COH$TCGA$L[, PURITY_OK, drop = FALSE], res = COH$TCGA$res,
           ax = COH$TCGA$ax[, PURITY_OK, drop = FALSE])
COV_PL <- cbind(COH$TCGA$cov[PURITY_OK, PROLIF_REF_COV, drop = FALSE],
                purity = PL$purity[PURITY_OK], leuko = PL$leuko[PURITY_OK])
stopifnot(!anyNA(COV_PL), nrow(COV_PL) == ncol(CP$L))

per_gene_pl <- list()
for (ax in rownames(CP$ax)) {
  per_gene_pl[[paste(ax, "adj. PROLIF_DISJOINT", sep = "|")]] <-
    .per_gene_rho(CP$L, CP$ax[ax, ],
                  cov = COV_PL[, PROLIF_REF_COV, drop = FALSE])
  per_gene_pl[[paste(ax, "adj. PROLIF + purity + leuko", sep = "|")]] <-
    .per_gene_rho(CP$L, CP$ax[ax, ], cov = COV_PL)
}
keep_pl <- which(!is.na(per_gene_pl[["MYC|adj. PROLIF_DISJOINT"]]))

purity_halves <- dplyr::bind_rows(lapply(SPLIT_SETS, function(nm) {
  ri  <- .row_index(COMPARATORS[[nm]], CP)
  own <- ri$row
  cmp <- .composition(ri$gene)
  stopifnot(length(own) == length(cmp))
  want <- lapply(split(own, cmp), function(v) intersect(v, keep_pl))
  want <- want[lengths(want) > 0L]
  parts <- .compartment_draws(CP, keep_pl, own, want)
  halves <- list(
    mitochondrial = names(parts)[names(parts) != "(not in MitoCarta)"],
    cytosolic     = intersect("(not in MitoCarta)", names(parts)))
  dplyr::bind_rows(lapply(names(halves), function(hn) {
    pp <- parts[halves[[hn]]]
    if (!length(pp)) return(NULL)
    idx <- unlist(lapply(pp, function(q) q$want), use.names = FALSE)
    dr <- replicate(NULL_DRAWS,
                    unlist(lapply(pp, function(q)
                      .matched_draw(q$want, q$B, exclude = own)),
                      use.names = FALSE),
                    simplify = FALSE)
    dplyr::bind_rows(lapply(names(per_gene_pl), function(key) {
      parts_k <- strsplit(key, "|", fixed = TRUE)[[1]]
      v <- per_gene_pl[[key]]
      nd <- vapply(dr, function(d) mean(v[d]), numeric(1))
      tibble::tibble(cohort = "TCGA", set = nm, half = hn, axis = parts_k[1],
                     adjustment = parts_k[2], n = length(idx),
                     mean_rho = mean(v[idx]), null_mean = mean(nd),
                     null_sd = stats::sd(nd), z = .z(mean(v[idx]), nd))
    }))
  }))
}))
message("\n   each half before and after purity + leukocyte fraction, ",
        "OXPHOS axis:")
purity_halves %>%
  dplyr::filter(axis == "OXPHOS") %>%
  dplyr::select(set, half, adjustment, n, mean_rho, null_mean, z) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(set, half, adjustment) %>% as.data.frame() %>%
  print(row.names = FALSE)

# --- WHICH cytosolic genes carry it, by Reactome module ---------------------
# Both cohorts, because the module label needs no covariate. If the negative is
# the death-receptor module alone, the infiltrate reading wins whatever the
# TCGA adjustment says; if the other modules are negative too, it does not.
module_break <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  ri <- .row_index(CANON, C)
  cmp <- .composition(ri$gene)
  cyt <- cmp == "(not in MitoCarta)"
  dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
    dplyr::bind_rows(lapply(names(ADJUSTMENTS), function(adj) {
      v <- per_gene[[paste(coh, ax, adj, sep = "|")]][ri$row[cyt]]
      tibble::tibble(cohort = coh, axis = ax, adjustment = adj,
                     gene = ri$gene[cyt],
                     module = unname(CANON_MODULE[ri$gene[cyt]]), rho = v)
    }))))
}))
message("\n   the 24 cytosolic genes by Reactome module (OXPHOS, adjusted):")
module_break %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::group_by(cohort, module) %>%
  dplyr::summarise(n = dplyr::n(), mean_rho = mean(rho),
                   median_rho = stats::median(rho),
                   n_negative = sum(rho < 0), .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\n   the same with the death-receptor module removed entirely:")
module_break %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT",
                module != "death receptor / extrinsic") %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(n = dplyr::n(), mean_rho = round(mean(rho), 3),
                   n_negative = sum(rho < 0), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

CYT_VERDICT <- {
  a <- purity_halves %>%
    dplyr::filter(axis == "OXPHOS", half == "cytosolic",
                  set == "apoptotic machinery (44)")
  after <- a$mean_rho[a$adjustment == "adj. PROLIF + purity + leuko"]
  before <- a$mean_rho[a$adjustment == "adj. PROLIF_DISJOINT"]
  nodr <- module_break %>%
    dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT",
                  module != "death receptor / extrinsic") %>%
    dplyr::group_by(cohort) %>% dplyr::summarise(m = mean(rho),
                                                 .groups = "drop")
  if (after < 0 && after <= before / 2 && all(nodr$m < 0)) {
    "SURVIVES: negative after purity and leukocyte fraction, and negative without the death-receptor module."
  } else if (after >= 0 || after > before / 2) {
    paste0("INFILTRATE: the cytosolic negative goes from ", round(before, 3),
           " to ", round(after, 3),
           " once purity and leukocyte fraction are removed.")
  } else {
    paste0("PARTIAL: survives the covariate (", round(after, 3),
           ") but the non-death-receptor modules are not all negative.")
  }
}
message("\n   VERDICT ON THE CONFOUND:\n   ", CYT_VERDICT)

# =============================================================================
# 6. THE LADDER - the machinery's 20 among real mitochondrial programmes
# =============================================================================
# A percentile among named pathways, not a z against synthetic draws. Pathways
# containing OXPHOS subunits or mitoribosome genes are KEPT AND MARKED rather
# than dropped: they are the ceiling of the scale and a ladder without its
# ceiling cannot be read. `is_apoptosis` marks MitoCarta's own Apoptosis
# pathway, which is very largely the machinery's mitochondrial half seen again
# through a different catalogue and is NOT an independent comparator.
message("\n6. the ladder: every MitoCarta leaf pathway of n >= ", LADDER_MIN)

# Leaf pathways only. The hierarchy nests - "Signaling" contains "Apoptosis" -
# and ranking a parent beside its own child double-counts the child's genes.
is_leaf <- vapply(names(MITO_PATHS), function(nm) {
  h <- MITO_PATH_HIER[[nm]]
  !any(startsWith(unlist(MITO_PATH_HIER, use.names = FALSE),
                  paste0(h, " > ")))
}, logical(1))
message("   ", sum(is_leaf), " leaf pathways of ", length(MITO_PATHS))

CANON_COMP <- .composition(CANON)
CANON_MITO <- names(CANON_COMP)[CANON_COMP != "(not in MitoCarta)"]
CANON_CYT  <- setdiff(CANON, CANON_MITO)
LADDER_SETS <- c(MITO_PATHS[is_leaf],
                 list(`* machinery, mitochondrial half` = CANON_MITO,
                      `* machinery, cytosolic half` = CANON_CYT,
                      `* machinery, all 44` = CANON))

ladder <- dplyr::bind_rows(lapply(names(COH), function(coh) {
  C <- COH[[coh]]
  keep <- which(!is.na(per_gene[[paste(coh, "MYC", "raw", sep = "|")]]))
  dplyr::bind_rows(lapply(names(LADDER_SETS), function(nm) {
    ri <- .row_index(LADDER_SETS[[nm]], C)
    inkeep <- ri$row %in% keep
    idx <- ri$row[inkeep]
    if (length(idx) < LADDER_MIN) return(NULL)
    # The compartments are taken over the SAME genes the correlations are, so
    # frac_mito describes what was measured and not what was asked for.
    cmp <- .composition(ri$gene[inkeep])
    stopifnot(length(idx) == length(cmp))
    tb <- table(cmp[cmp != "(not in MitoCarta)"])
    dplyr::bind_rows(lapply(c("MYC", "OXPHOS"), function(ax)
      dplyr::bind_rows(lapply(names(ADJUSTMENTS), function(adj) {
        v <- per_gene[[paste(coh, ax, adj, sep = "|")]][idx]
        tibble::tibble(
          cohort = coh, pathway = nm, axis = ax, adjustment = adj,
          n = length(idx),
          hierarchy = if (nm %in% names(MITO_PATH_HIER))
            MITO_PATH_HIER[[nm]] else "(this study)",
          dominant_compartment = if (length(tb)) names(tb)[which.max(tb)]
            else "(none)",
          frac_mito = mean(cmp != "(not in MitoCarta)"),
          is_target = startsWith(nm, "* "),
          is_apoptosis = nm == "Apoptosis",
          has_oxphos_gene = length(intersect(ri$gene[inkeep], OXARM)) > 0L,
          has_mitoribo_gene = length(intersect(ri$gene[inkeep],
                                               MITORIBO)) > 0L,
          frac_ceiling = length(intersect(ri$gene[inkeep],
                                          union(OXARM, MITORIBO))) /
            length(idx),
          mean_rho = mean(v), mean_abs_rho = mean(abs(v)),
          sd_rho = stats::sd(v))
      }))))
  }))
}))

# The percentile the machinery's mitochondrial half occupies, computed against
# pathways that are neither the target itself nor an OXPHOS/mitoribosome
# ceiling nor MitoCarta's own re-statement of apoptosis.
# TWO REFERENCE SETS, BOTH REPORTED, BECAUSE THE STRICT ONE IS THIN. Excluding
# a pathway for holding a SINGLE OXPHOS-arm gene removes Fe-S biosynthesis for
# NDUFAB1 alone and leaves about 30 pathways to take a percentile against. The
# lenient set excludes only pathways that are SUBSTANTIALLY ceiling. Neither is
# chosen after the fact: both columns are printed for every row.
ladder <- ladder %>%
  dplyr::mutate(mostly_ceiling = frac_ceiling > 0.2)
ladder_ref <- ladder %>%
  dplyr::filter(!is_target, !is_apoptosis, !has_oxphos_gene,
                !has_mitoribo_gene)
ladder_ref_lenient <- ladder %>%
  dplyr::filter(!is_target, !is_apoptosis, !mostly_ceiling)
ladder_pct <- ladder %>%
  dplyr::filter(is_target) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n_reference_pathways = sum(ladder_ref$cohort == cohort &
                                 ladder_ref$axis == axis &
                                 ladder_ref$adjustment == adjustment),
    pct_mean_rho = mean(ladder_ref$mean_rho[
      ladder_ref$cohort == cohort & ladder_ref$axis == axis &
        ladder_ref$adjustment == adjustment] < mean_rho),
    pct_sd_rho = mean(ladder_ref$sd_rho[
      ladder_ref$cohort == cohort & ladder_ref$axis == axis &
        ladder_ref$adjustment == adjustment] < sd_rho),
    n_reference_lenient = sum(ladder_ref_lenient$cohort == cohort &
                                ladder_ref_lenient$axis == axis &
                                ladder_ref_lenient$adjustment == adjustment),
    pct_mean_rho_lenient = mean(ladder_ref_lenient$mean_rho[
      ladder_ref_lenient$cohort == cohort & ladder_ref_lenient$axis == axis &
        ladder_ref_lenient$adjustment == adjustment] < mean_rho)) %>%
  dplyr::ungroup()
message("\n   where the machinery sits among mitochondrial programmes ",
        "(adjusted):")
ladder_pct %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::select(cohort, axis, pathway, n, mean_rho, pct_mean_rho,
                pct_mean_rho_lenient, sd_rho, pct_sd_rho,
                n_reference_pathways, n_reference_lenient) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::arrange(axis, pathway, cohort) %>% as.data.frame() %>%
  print(row.names = FALSE)

message("\n   the compartment gradient over curated pathways, not over 2 ",
        "genes (OXPHOS,\n   adjusted, ceiling pathways excluded):")
ladder_ref %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::group_by(cohort, dominant_compartment) %>%
  dplyr::summarise(n_pathways = dplyr::n(),
                   median_mean_rho = stats::median(mean_rho),
                   .groups = "drop") %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\n   the ten pathways that track OXPHOS most and least strongly ",
        "(TCGA, adjusted):")
ladder %>%
  dplyr::filter(cohort == "TCGA", axis == "OXPHOS",
                adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::arrange(dplyr::desc(mean_rho)) %>%
  dplyr::select(pathway, n, dominant_compartment, mean_rho, sd_rho,
                has_oxphos_gene, has_mitoribo_gene) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
  { rbind(head(., 10), tail(., 10)) } %>% as.data.frame() %>%
  print(row.names = FALSE)

# =============================================================================
# 7. Figures
# =============================================================================
message("\n7. figures")

COHORT_COLS <- c(TCGA = "#1f4e79", `SCAN-B` = "#b8541a")
theme_e14 <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey92",
                                                          colour = NA),
                 legend.position = "bottom",
                 plot.caption = ggplot2::element_text(size = 8,
                                                      colour = "grey40",
                                                      hjust = 0))
.save <- function(p, nm, w, h) {
  for (ext in c("png", "pdf"))
    ggplot2::ggsave(file.path(DIR_FIGURES, paste0(nm, ".", ext)), p,
                    width = w, height = h, dpi = 300)
  message("   ", nm); invisible(p)
}

# --- FIG 1: THE HEAD-TO-HEAD, and it is the whole point of the script -------
# Each row is one programme. The point is its observed split; the grey bar
# behind it is its OWN composition-matched null, mean +/- 1 SD. If the
# machinery's point sits where the other programmes' points sit, the ordering
# is what any outer-membrane-spanning programme shows.
fig1_d <- comparator_splits %>%
  dplyr::filter(set %in% PLOT_SETS) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                is_target = set == "apoptotic machinery (44)",
                set = factor(set, levels = rev(PLOT_SETS)))
p1 <- ggplot2::ggplot(fig1_d, ggplot2::aes(y = set)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_hline(
    yintercept = which(rev(PLOT_SETS) == "apoptotic machinery (44)"),
    colour = "grey85", linewidth = 6) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = null_mean - null_sd, xmax = null_mean + null_sd),
    colour = "grey65", linewidth = 3, alpha = 0.55) +
  ggplot2::geom_point(ggplot2::aes(x = null_mean), colour = "grey35",
                      shape = 124, size = 2.5) +
  ggplot2::geom_point(ggplot2::aes(x = split, colour = cohort,
                                   shape = is_target), size = 2.6) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16),
                              guide = "none") +
  ggplot2::facet_grid(cohort ~ axis + adjustment) +
  ggplot2::labs(
    x = "split: Spearman of per-gene rho against the MitoCarta 0/1 label",
    y = NULL,
    title = "Does any non-apoptotic programme spanning the outer membrane split the same way?",
    caption = paste0(
      "Grey bar: that programme's OWN sub-compartment- and expression-matched null, mean +/- 1 SD,\n",
      NULL_DRAWS, " draws from MitoCarta minus the OXPHOS arm and the mitoribosome. Coloured point:\n",
      "observed, diamond and grey band for the target. The reading is ACROSS ROWS, not against the\n",
      "grey. `isozyme pairs` is a declared CEILING - the same reaction either side of the membrane -\n",
      "and is author-curated, not from a pinned catalogue. The PINK1/PRKN subset of mitophagy is\n",
      "measured and sign-checked in section 4.1 but not drawn: it lands on top of mitophagy.\n",
      "EXPLORATORY: nothing pre-registered; three comparators cannot make a p-value.")) +
  theme_e14
.save(p1, "E14_fig1_head_to_head_split", 12, 6)

# --- FIG 2: which half carries it ------------------------------------------
fig2_d <- half_tests %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT", set %in% PLOT_SETS) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                set = factor(set, levels = rev(PLOT_SETS)))
p2 <- ggplot2::ggplot(fig2_d, ggplot2::aes(y = set)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = null_mean - null_sd, xmax = null_mean + null_sd),
    colour = "grey65", linewidth = 3, alpha = 0.55) +
  ggplot2::geom_point(ggplot2::aes(x = mean_rho, colour = cohort,
                                   shape = set == "apoptotic machinery (44)"),
                      size = 2.6) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16),
                              guide = "none") +
  ggplot2::facet_grid(half ~ axis, scales = "free_x") +
  ggplot2::labs(
    x = "mean per-gene rho (signed), proliferation-adjusted", y = NULL,
    title = "A split has two halves, and they are different findings",
    caption = paste0(
      "Top: the mitochondrial half against mitochondrial genes matched on sub-compartment and\n",
      "expression. Bottom: the cytosolic half against non-MitoCarta genes matched on expression.\n",
      "A high mitochondrial half means apoptotic proteins at the organelle are extra-coupled to the\n",
      "OXPHOS programme; a low cytosolic half means the split is made by the cytosolic genes being\n",
      "unusually flat. The split statistic alone cannot distinguish these.")) +
  theme_e14
.save(p2, "E14_fig2_which_half", 10, 6)

# --- FIG 3: the ladder ------------------------------------------------------
# The percentile figure. Every MitoCarta leaf pathway is a point; the
# machinery's halves are marked. Ceiling pathways are drawn in a different
# shape rather than hidden, because a ladder without its top rung is unreadable.
fig3_d <- ladder %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                kind = dplyr::case_when(
                  is_target ~ "the machinery",
                  is_apoptosis ~ "MitoCarta Apoptosis (not independent)",
                  has_oxphos_gene | has_mitoribo_gene ~ "ceiling: contains OXPHOS/mitoribosome genes",
                  TRUE ~ "other mitochondrial programme"))
ord <- fig3_d %>% dplyr::filter(cohort == "TCGA") %>%
  dplyr::arrange(mean_rho) %>% dplyr::pull(pathway)
p3 <- ggplot2::ggplot(
  fig3_d %>% dplyr::mutate(pathway = factor(pathway, levels = ord)),
  ggplot2::aes(x = mean_rho, y = pathway, colour = kind, size = kind)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_point() +
  ggplot2::scale_colour_manual(values = c(
    "the machinery" = "#b8541a",
    "MitoCarta Apoptosis (not independent)" = "#d9a441",
    "ceiling: contains OXPHOS/mitoribosome genes" = "grey55",
    "other mitochondrial programme" = "#1f4e79"), name = NULL) +
  ggplot2::scale_size_manual(values = c(
    "the machinery" = 2.4,
    "MitoCarta Apoptosis (not independent)" = 2.0,
    "ceiling: contains OXPHOS/mitoribosome genes" = 1.1,
    "other mitochondrial programme" = 1.4), guide = "none") +
  ggplot2::facet_wrap(~ cohort, nrow = 1) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
  ggplot2::labs(
    x = "mean per-gene rho with OXPHOS, proliferation-adjusted", y = NULL,
    title = "The machinery among named mitochondrial programmes",
    caption = paste0(
      "Every MitoCarta 3.0 leaf pathway with at least ", LADDER_MIN, " genes in the matrix. Ordered by TCGA.\n",
      "A percentile among programmes a reader already knows, in place of a z against draws they\n",
      "cannot see. Ceiling pathways are marked, not dropped: a ladder without its top rung is not\n",
      "readable. MitoCarta's own Apoptosis pathway is the machinery seen through a second catalogue\n",
      "and is NOT an independent comparator.")) +
  theme_e14 +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.5))
.save(p3, "E14_fig3_mitocarta_ladder", 11, 11)

# --- FIG 4: the compartment gradient, over pathways rather than over 2 genes -
fig4_d <- ladder_ref %>%
  dplyr::filter(adjustment == "adj. PROLIF_DISJOINT",
                dominant_compartment %in% c("MOM", "IMS", "MIM", "Matrix")) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                dominant_compartment = factor(dominant_compartment,
                                              levels = c("MOM", "IMS", "MIM",
                                                         "Matrix")))
mach_line <- ladder %>%
  dplyr::filter(is_target, adjustment == "adj. PROLIF_DISJOINT",
                pathway == "* machinery, mitochondrial half") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
p4 <- ggplot2::ggplot(fig4_d, ggplot2::aes(x = dominant_compartment,
                                           y = mean_rho)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey93", colour = "grey45",
                        width = 0.55) +
  ggplot2::geom_jitter(ggplot2::aes(colour = cohort), width = 0.13,
                       height = 0, size = 1.3, alpha = 0.85) +
  ggplot2::geom_hline(data = mach_line,
                      ggplot2::aes(yintercept = mean_rho), colour = "#b8541a",
                      linetype = "dashed", linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::facet_grid(cohort ~ axis) +
  ggplot2::labs(
    x = "the pathway's dominant sub-mitochondrial compartment",
    y = "mean per-gene rho, proliferation-adjusted",
    title = "The compartment gradient, measured over pathways instead of over two genes",
    caption = paste0(
      "One point per MitoCarta leaf pathway, grouped by where most of its members sit. Dashed line:\n",
      "the machinery's own mitochondrial half. E11 section 3.2 had to read this gradient off the\n",
      "machinery's 2 inner-membrane genes; here it is read off ", nrow(fig4_d) / 4,
      " pathways per cell.\n",
      "Ceiling pathways (OXPHOS subunits, mitoribosome) are excluded from this panel.")) +
  theme_e14
.save(p4, "E14_fig4_compartment_gradient", 9, 7)

# --- FIG 5: the falsifier for figure 1, and it is not optional --------------
# Figure 1 says the machinery out-splits every comparator. Section 5 says the
# reason is entirely its cytosolic half. This panel is what stops that being an
# immune-infiltrate result: the same cytosolic halves before and after purity
# and leukocyte fraction, beside the 24 genes broken down by which apoptotic
# module they belong to. Figure 1 must not be shown without it.
fig5a <- purity_halves %>%
  dplyr::filter(axis == "OXPHOS", half == "cytosolic", set %in% PLOT_SETS) %>%
  dplyr::mutate(set = factor(set, levels = rev(PLOT_SETS)),
                adjustment = factor(adjustment,
                                    levels = c("adj. PROLIF_DISJOINT",
                                               "adj. PROLIF + purity + leuko")))
p5a <- ggplot2::ggplot(fig5a, ggplot2::aes(x = mean_rho, y = set)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = null_mean - null_sd, xmax = null_mean + null_sd),
    colour = "grey65", linewidth = 3, alpha = 0.55) +
  ggplot2::geom_point(ggplot2::aes(shape = adjustment), size = 2.4,
                      colour = COHORT_COLS[["TCGA"]]) +
  ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
  ggplot2::labs(x = "mean rho with OXPHOS, cytosolic half only", y = NULL,
                subtitle = "TCGA, n = 1007 with both estimates") +
  theme_e14
fig5b <- module_break %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT") %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)))
p5b <- ggplot2::ggplot(fig5b, ggplot2::aes(x = rho,
                                           y = stats::reorder(module, rho))) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_point(ggplot2::aes(colour = cohort), size = 1.8, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.4)) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::labs(x = "per-gene rho with OXPHOS, proliferation-adjusted",
                y = NULL,
                subtitle = "the 24 cytosolic genes, one point per gene") +
  theme_e14
p5 <- patchwork::wrap_plots(p5a, p5b, nrow = 2, heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    title = "Is the cytosolic half immune infiltrate rather than apoptosis?",
    caption = paste0(
      "Top: filled = proliferation only, open = proliferation plus purity plus leukocyte fraction.\n",
      "Grey bar is each half's own expression-matched null of non-MitoCarta genes. Bottom: the same\n",
      "24 genes grouped by apoptotic module. The death-receptor and NF-kB genes are the infiltrate\n",
      "hazard (trap 2: rho(OXPHOS, leukocyte fraction) = -0.158); the intrinsic-pathway cytosolic\n",
      "genes are not. SCAN-B has no purity estimate and it is never imputed."),
    theme = theme_e14)
.save(p5, "E14_fig5_infiltrate_falsifier", 9, 8)

# --- FIG 6: THE SUPPLEMENTARY PANEL. One picture, the specificity claim ------
# Figures 1, 2 and 5 make the argument in three steps. This makes it in one,
# for a supplementary slot: each programme is a point whose x is what its
# MITOCHONDRIAL members do with OXPHOS and whose y is what its CYTOSOLIC
# members do.
#
# THE AXES ARE ON THE SAME SCALE AND THE PANEL IS SQUARE, so "they agree
# horizontally and disagree vertically" is a fact about the picture rather
# than about how it was drawn. Both zeros are in view for the same reason.
#
# FOUR PROGRAMMES, NOT FIVE - and by section 4.1 that is a decision about ink,
# not about evidence. The PINK1/PRKN subset is mitophagy minus the receptor
# arm, agrees with it in sign in both cohorts, and lands on top of it in every
# panel. EVERY figure in this script now omits it, so PANEL_SETS is PLOT_SETS.
#
# THE LEADER LINES ANCHOR TO THE NEAREST COHORT POINT, not to the midpoint of
# the pair, because midpoint anchors made three leaders cross the Fe-S
# segment. Label positions are FIXED COORDINATES, not repelled: a seed-driven
# layout moves between runs and a figure that goes in a paper must not.
message("   composing the one-panel supplementary figure")

PANEL_SETS <- PLOT_SETS
PANEL_SHORT <- c(`apoptotic machinery (44)` = "apoptotic machinery",
                 mitophagy = "mitophagy",
                 `Fe-S cluster assembly` = "Fe-S cluster assembly",
                 `isozyme pairs (CEILING)` = "isozyme pairs")
LABPOS <- tibble::tribble(
  ~set,                       ~lx,    ~ly,   ~hj,
  "isozyme pairs (CEILING)",  0.055,  0.222, 0.0,
  "Fe-S cluster assembly",    0.385,  0.170, 1.0,
  "mitophagy",                0.245,  0.052, 0.0,
  "apoptotic machinery (44)", 0.235, -0.120, 0.0)
stopifnot(all(PANEL_SETS %in% SPLIT_SETS), all(LABPOS$set %in% PANEL_SETS))

panel_data <- half_tests %>%
  dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT",
                set %in% PANEL_SETS) %>%
  dplyr::select(cohort, set, half, n, mean_rho, null_mean, null_sd, z) %>%
  tidyr::pivot_wider(names_from = half,
                     values_from = c(mean_rho, null_mean, null_sd, n, z)) %>%
  dplyr::mutate(cohort = factor(cohort, levels = names(COHORT_COLS)),
                is_target = set == "apoptotic machinery (44)")

# The cytosolic matched null is nearly the same for all four, so it is drawn
# once as a band. The MITOCHONDRIAL null is not (0.105 to 0.184, because the
# sets differ in sub-compartment) and is reported in the caption instead of
# drawn as a band that would be wrong for three of the four rows.
CYT_NULL <- panel_data %>%
  dplyr::summarise(m = mean(null_mean_cytosolic),
                   s = mean(null_sd_cytosolic))

panel_labs <- LABPOS %>%
  dplyr::rowwise() %>%
  dplyr::mutate({
    d <- panel_data[panel_data$set == set, ]
    near <- which.min((d$mean_rho_mitochondrial - lx)^2 +
                        (d$mean_rho_cytosolic - ly)^2)
    tibble::tibble(px = d$mean_rho_mitochondrial[near],
                   py = d$mean_rho_cytosolic[near],
                   nm = d$n_mitochondrial[near], nc = d$n_cytosolic[near],
                   is_target = d$is_target[near])
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lab = paste0(unname(PANEL_SHORT[set]), " (", nm, "/", nc, ")"))

p6 <- ggplot2::ggplot(panel_data, ggplot2::aes(x = mean_rho_mitochondrial,
                                               y = mean_rho_cytosolic)) +
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                    ymin = CYT_NULL$m - CYT_NULL$s,
                    ymax = CYT_NULL$m + CYT_NULL$s,
                    fill = "grey55", alpha = 0.16) +
  ggplot2::annotate("text", x = 0.393, y = CYT_NULL$m + CYT_NULL$s - 0.006,
                    label = "matched null for a cytosolic gene set, mean +/- 1 SD",
                    hjust = 1, vjust = 1, size = 2.5, colour = "grey35") +
  ggplot2::geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  ggplot2::geom_segment(data = panel_labs,
                        ggplot2::aes(x = lx, y = ly, xend = px, yend = py),
                        colour = "grey72", linewidth = 0.3,
                        inherit.aes = FALSE) +
  ggplot2::geom_line(ggplot2::aes(group = set), colour = "grey62",
                     linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(colour = cohort, shape = is_target),
                      size = 3.3) +
  ggplot2::geom_text(data = panel_labs,
                     ggplot2::aes(x = lx, y = ly, label = lab, hjust = hj),
                     size = 3.0, inherit.aes = FALSE,
                     colour = ifelse(panel_labs$is_target, "#8a3d12",
                                     "grey15"),
                     fontface = ifelse(panel_labs$is_target, "bold",
                                       "plain")) +
  ggplot2::scale_colour_manual(values = COHORT_COLS, name = NULL) +
  ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16),
                              guide = "none") +
  ggplot2::coord_fixed(xlim = c(-0.02, 0.40), ylim = c(-0.17, 0.25),
                       expand = FALSE) +
  ggplot2::labs(
    x = "mitochondrial (MitoCarta) members: mean rho with OXPHOS",
    y = "cytosolic members: mean rho with OXPHOS",
    title = "Apoptosis is the only programme whose cytosolic members run against OXPHOS",
    caption = paste0(
"Four curated programmes that span the outer mitochondrial membrane, each split into its\n",
"MitoCarta and non-MitoCarta members (n mito / n cytosolic in the label). Proliferation-adjusted\n",
"partial Spearman; the grey line joins the two cohorts of one programme.\n",
"HORIZONTALLY the programmes agree: every mitochondrial half lies between +0.14 and +0.31, at\n",
"z -0.9 to +2.1 against its own expression- and sub-compartment-matched null.\n",
"VERTICALLY they do not: every cytosolic half is positive and above its null (z +1.7 to +3.4)\n",
"except the apoptotic machinery's, which is negative in both cohorts (-0.106 TCGA, -0.094\n",
"SCAN-B). That negative survives purity and leukocyte fraction (-0.091, TCGA n = 1007) and\n",
"survives deleting the death-receptor module entirely (9 of the 14 remaining genes still\n",
"negative in both cohorts), so it is not immune infiltrate.\n",
"Fe-S cluster assembly is the one cohort-inconsistent comparator - its long grey line. The\n",
"PINK1/PRKN-only subset of mitophagy is mitophagy minus its BH3-domain receptors; it agrees with\n",
"mitophagy in sign in both cohorts (+0.155, +0.104 cytosolic) and is on no figure here. Its\n",
"values are in the E14 table.\n",
"THE CLAIM IS THE CONTRAST BETWEEN PROGRAMMES, not the distance from the null band.\n",
"EXPLORATORY: nothing pre-registered; four programmes cannot make a p-value.")) +
  theme_e14 +
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 7.2,
                                                      colour = "grey40",
                                                      hjust = 0))
.save(p6, "E14_fig6_specificity_one_panel", 8.2, 8.6)

# =============================================================================
# 8. Save
# =============================================================================
message("\n8. save")

readr::write_csv(comparator_splits, PATH_E14_CSV)
readr::write_csv(ladder, PATH_E14_CSV2)

saveRDS(list(
  comparator_audit = comparator_audit, eligibility = elig,
  comparator_splits = comparator_splits, verdict_tbl = verdict_tbl,
  verdict = VERDICT, half_tests = half_tests,
  ladder = ladder, ladder_pct = ladder_pct, drops = drops,
  purity_halves = purity_halves, module_break = module_break,
  panel_data = panel_data, panel_labs = panel_labs,
  cytosolic_verdict = CYT_VERDICT,
  comparators = COMPARATORS, comparator_source = COMPARATOR_SOURCE,
  paralogue_pairs = PARALOGUE_PAIRS,
  split_sets = SPLIT_SETS, plot_sets = PLOT_SETS, not_plotted = NOT_PLOTTED,
  mito_pair = mito_pair, mito_pair_md_agree = MITO_PAIR_MD_AGREE,
  level_only = LEVEL_ONLY,
  settings = list(null_draws = NULL_DRAWS, n_bins = N_BINS,
                  min_half = MIN_HALF, ladder_min = LADDER_MIN,
                  prolif_covariate = PROLIF_REF_COV, myc_axis = MYC_REF,
                  adjustments = names(ADJUSTMENTS), seed = PROJECT_SEED,
                  msigdb_version = MSIGDB_VERSION,
                  gene_scale = "linear DESeq2-normalised, rank-transformed"),
  rules = list(
    question = paste("E11 P4a left a 1.3-1.6 SD residue above a",
                     "sub-compartment-matched RANDOM null. A random draw's",
                     "members share only expression and compartment; a real",
                     "pathway's members are co-regulated with each other, and",
                     "that is the property that could make a split. So the",
                     "residue is closed by a NAMED non-apoptotic programme,",
                     "not by more draws."),
    primary_comparator = paste("REACTOME_MITOPHAGY. Outer-membrane-weighted",
                               "mitochondrial half against a cytosolic",
                               "effector half whose members translocate to the",
                               "outer membrane on activation - the same shape",
                               "and the same hazard as the 44, and zero genes",
                               "in common with them."),
    import_fails = paste("mitochondrial protein import was the first",
                         "candidate and it cannot be used: 62 of 62 genes are",
                         "in MitoCarta, so there is no cytosolic half to split",
                         "against. Same for cristae formation (31 of 31) and",
                         "mitochondrial calcium transport (21 of 22). They are",
                         "level comparators only, in section 6."),
    ceiling = paste("the isozyme pairs are AUTHOR-CURATED IN THIS SCRIPT and",
                    "are the only set here not from a pinned catalogue. Same",
                    "reaction either side of the membrane, so function is held",
                    "constant and only compartment varies. Several are TCA",
                    "enzymes and the mitochondrial tRNA synthetases",
                    "co-regulate with the mitoribosome: a large split there is",
                    "expected and is a scale marker, not a finding."),
    two_halves = paste("section 5 exists because a split above its null can",
                       "come from the mitochondrial half being high OR the",
                       "cytosolic half being low. Those are different findings",
                       "and the split statistic cannot tell them apart."),
    statistic = paste("split stays the anchor because it is what every number",
                      "written up before 2026-09-02 was computed on;",
                      "med_diff is carried beside it and is what the figures",
                      "plot. Choosing between them after seeing their z values",
                      "would be statistic-shopping."),
    not_plotted = paste("the PINK1/PRKN mitophagy variant is measured, saved",
                        "and sign-checked in section 4.1 but is on no figure.",
                        "It is mitophagy minus the BH3-domain receptors, so it",
                        "lands on top of mitophagy in every panel. Section 4.1",
                        "stops the script if the two disagree in sign on the",
                        "anchor statistic, which is what keeps 'not drawn'",
                        "different from 'dropped'. THEY DO AGREE ON `split`",
                        "AND DISAGREE IN SIGN ON `med_diff` IN TCGA (+0.184",
                        "against -0.030), so the mitophagy comparator is",
                        "softer than the figures alone suggest."),
    supplementary_panel = paste("E14_fig6 is the one-panel version of the",
                                "whole argument, for a supplementary slot. Its",
                                "axes are on the same scale and the panel is",
                                "square, so 'they agree horizontally and",
                                "disagree vertically' is a fact about the",
                                "picture. Label positions are fixed",
                                "coordinates, not repelled - a seed-driven",
                                "layout moves between runs."),
    no_p_value = paste("four structural comparators cannot produce a p-value.",
                       "The deliverable is a ranking that reproduces across",
                       "two cohorts, which is what this repo counts as",
                       "evidence."),
    infiltrate = paste("section 5.1 is the falsifier for section 4. The 44's",
                       "cytosolic half is the death-receptor and NF-kB arm -",
                       "FAS, FASLG, TNF, TNFRSF10A/B, TRADD, NFKB1, RELA -",
                       "which is an inflammatory module as much as an",
                       "apoptotic one, and trap 2 puts rho(OXPHOS, leukocyte",
                       "fraction) at -0.158. Infiltrate alone would give a",
                       "negative cytosolic half and an apoptosis-specific",
                       "verdict with no apoptosis in it. The comparator cannot",
                       "control for this; only the covariate can, and SCAN-B",
                       "has no purity estimate."),
    both_cohorts = paste("every comparator is cut to the genes resolvable in",
                         "BOTH matrices before anything is measured. The",
                         "isozyme set arrived as 60 genes in TCGA and 49 in",
                         "SCAN-B - eleven cytosolic tRNA synthetases whose",
                         "2014 names the pinned symbol map does not carry -",
                         "and the loss fell entirely on one half, which would",
                         "have read as biology."),
    mitocarta_meaning = paste("MitoCarta is a proteome catalogue and a",
                              "transcript has no idea where its protein ends",
                              "up. Membership here marks membership of the",
                              "nuclear-encoded mitochondrial regulon, not",
                              "localisation."))), PATH_E14)
message("   ", PATH_E14)
message("   ", PATH_E14_CSV)
message("   ", PATH_E14_CSV2)
message("\nE14 done.\n", strrep("=", 78))

# =============================================================================
# SANDBOX - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  x <- readRDS(file.path(DIR_RESULTS, "curated_comparators.rds"))

  # THE ANSWER. Read across the `set` column, not down it.
  x$comparator_splits %>%
    dplyr::filter(axis == "OXPHOS",
                  adjustment == "adj. PROLIF_DISJOINT") %>%
    dplyr::select(cohort, set, split, null_mean, z, med_diff) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::arrange(cohort, dplyr::desc(split)) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # The verdict against the falsifier that was written before the run.
  cat(x$verdict, "\n")
  x$verdict_tbl %>% as.data.frame() %>% print(row.names = FALSE)

  # Why protein import could not be used - the eligibility table says so.
  x$comparator_audit %>%
    dplyr::filter(cohort == "TCGA") %>%
    dplyr::select(set, n_in_matrix, n_mito, n_cytosolic, MOM, IMS, MIM,
                  Matrix, overlap_44) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # Which half carries the residue.
  x$half_tests %>%
    dplyr::filter(axis == "OXPHOS",
                  adjustment == "adj. PROLIF_DISJOINT",
                  set %in% c("apoptotic machinery (44)", "mitophagy")) %>%
    dplyr::select(cohort, set, half, n, mean_rho, null_mean, z) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # Where the machinery sits among named mitochondrial programmes.
  x$ladder_pct %>%
    dplyr::filter(adjustment == "adj. PROLIF_DISJOINT",
                  axis == "OXPHOS") %>%
    dplyr::select(cohort, pathway, n, mean_rho, pct_mean_rho, sd_rho,
                  pct_sd_rho) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # The compartment gradient over pathways rather than over 2 genes.
  x$ladder %>%
    dplyr::filter(axis == "OXPHOS", adjustment == "adj. PROLIF_DISJOINT",
                  !is_target, !has_oxphos_gene, !has_mitoribo_gene) %>%
    dplyr::group_by(cohort, dominant_compartment) %>%
    dplyr::summarise(n_pathways = dplyr::n(),
                     median_rho = round(stats::median(mean_rho), 3),
                     .groups = "drop") %>%
    as.data.frame() %>% print(row.names = FALSE)

  # THE FALSIFIER FOR THE HEAD-TO-HEAD: is the cytosolic half infiltrate?
  cat(x$cytosolic_verdict, "\n")
  x$purity_halves %>%
    dplyr::filter(axis == "OXPHOS", half == "cytosolic") %>%
    dplyr::select(set, adjustment, n, mean_rho, null_mean, z) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # Which cytosolic genes carry the negative.
  x$module_break %>%
    dplyr::filter(axis == "OXPHOS",
                  adjustment == "adj. PROLIF_DISJOINT") %>%
    dplyr::group_by(cohort, module) %>%
    dplyr::summarise(n = dplyr::n(), mean_rho = round(mean(rho), 3),
                     n_negative = sum(rho < 0), .groups = "drop") %>%
    as.data.frame() %>% print(row.names = FALSE)

  # The control that licenses leaving PINK1/PRKN off every figure.
  x$mito_pair %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    as.data.frame() %>% print(row.names = FALSE)

  # The supplementary panel's own numbers, four programmes, two cohorts.
  x$panel_data %>%
    dplyr::select(cohort, set, mean_rho_mitochondrial, z_mitochondrial,
                  mean_rho_cytosolic, z_cytosolic) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::arrange(set, cohort) %>% as.data.frame() %>%
    print(row.names = FALSE)

  # The author-curated set, named in full so it can be argued with.
  x$paralogue_pairs %>% as.data.frame() %>% print(row.names = FALSE)
}
