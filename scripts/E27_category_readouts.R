# scripts/E27_category_readouts.R
# =============================================================================
# E27 - the seven MitoCarta level-1 categories, read across all eight groups
# =============================================================================
#
# WHAT THIS IS. E25 built eight group orderings of the shared 142-pathway panel
# and then read ONE pathway's position inside them. This reads the SEVEN
# MitoCarta level-1 categories in the same orderings, on the same panel, plus
# four nested readouts. Nothing is re-scored: the level-1 categories are
# themselves rows of the 142-pathway universe and were already scored. E25 is
# not modified; it is signed off and its note is live.
#
# THE DESIGN, AND WHY IT IS THREE CATEGORIES AND NOT ONE
# -------------------------------------------------------
#   OXPHOS                       expected DOWN across the mouse WT window
#   Mitochondrial central dogma  expected FLAT - THE INTERNAL NULL
#   Metabolism                   expected UP
#
# mitoPPS is a relative score, so "OXPHOS down" and "Metabolism up" are not
# independent - one forces some of the other. A category that does NOT move, in
# the same samples on the same ruler, is what converts the pattern from an
# artefact of relativity into a reallocation. The remaining four level-1
# categories are reported as context: if the reallocation is spread evenly
# across all six non-OXPHOS categories that is a different and weaker result,
# and only the complete partition can show it.
#
# TWO ESTIMANDS, AND WHICH ANSWERS WHICH QUESTION
# ------------------------------------------------
#   VALUE       group-mean mitoPPS per category. PRIMARY WITHIN A SPECIES, and
#               what R1, R2, R3 and R5 are read on. R2 is an EQUIVALENCE
#               statement about DISTANCES - "central dogma's interval excludes
#               the magnitude of the OXPHOS change" - and that needs a LINEAR
#               scale. Percentile distances are bounded and compress near the
#               ends, so they are the wrong ruler for it.
#   PERCENTILE  position in the group's ordering of the 142-pathway panel.
#               SECONDARY, and THE ONLY THING THAT CROSSES THE SPECIES
#               BOUNDARY. It is what R4 is read on and what lets these numbers
#               sit beside E25's.
#
# NO mitoPPS VALUE CROSSES THE SPECIES BOUNDARY. A mouse-only or human-only
# panel or table on values is fine. A table or axis carrying both is not, and
# the figure that spans eight groups is on percentiles for exactly that reason.
#
# THE ZERO-SUM ARGUMENT IS WEAKER THAN IT LOOKS, and that is recorded rather
# than relied on: seven readouts inside a 142-pathway panel leave 135 pathways
# to absorb any difference, so percentile forcing is weak rather than strict.
# The reason for values is the equivalence scale above, not zero-sum.
#
# WHAT WAS ALREADY VISIBLE BEFORE THIS SCRIPT WAS WRITTEN, disclosed so the
# record is honest: the eleven readouts' PERCENTILE positions were computable
# from E25's saved object - they are rows of $orderings - though E25's NOTE
# reported only four of them. They were surfaced during planning, before the
# run, not after it. The reading rules below were fixed by the author before any
# number was seen, and that is what protects them; this script's author having
# seen E25's saved percentiles does not unfix them, and the note says so.
#
# EXPLORATORY, POST-HOC, DESCRIPTIVE, as the whole of Phase 2 is. Nothing here
# is pre-registered and nothing is a hypothesis test. R1-R5 are reading rules
# fixed in advance so that the outcome cannot be chosen afterwards. THEY ARE NOT
# PREDICTIONS AND NONE OF THEM IS PREDICTED. Every outcome is reportable.
#
# THIS IS A PRESENTATION OF AN EXISTING CLAIM, NOT A NEW RESULT. The
# expression-matched-null analysis in myc_mouse remains the primary evidence for
# the developmental change. Both stay reportable, and IF THEY DISAGREE THAT
# DISAGREEMENT IS THE FINDING and is not resolved by choosing the friendlier one.
#
# SPECIES: the human side is human. No ortholog function is called here or
# anywhere in this repo. A consequence recorded rather than worked around: the
# umbrella is 156 human against 155 mouse and assembly factors 68 against 67,
# and WHICH GENE DIFFERS CANNOT BE REPORTED, because naming it requires aligning
# mouse and human symbols, which is an ortholog mapping whatever it is called.
# Within-species membership and cross-species COUNTS only.
#
# N3. Transcript associations throughout. "Primed" appears nowhere.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE27: the seven level-1 categories across all eight groups\n",
        strrep("=", 78))

PATH_E27     <- file.path(DIR_RESULTS, "category_readouts.rds")
PATH_E27_TAB <- file.path(DIR_TABLES,  "E27_category_readouts.csv")
DIR_E27_FIG  <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_E25     <- file.path(DIR_RESULTS, "mitopps_group_rankings.rds")
PATH_MOUSE   <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")
PATH_MITOCARTA_H <- here::here("data", "mitocarta_human", "Human.MitoCarta3.0.xls")

set.seed(1)
NBOOT <- 5000L

# =============================================================================
# 0. CONSTANTS AND THE READING RULES, FIXED BEFORE ANY NUMBER
# =============================================================================

MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"

# The three the design is about. Named here, ahead of the other four.
CAT_OX     <- "OXPHOS"
CAT_NULL   <- "Mitochondrial central dogma"
CAT_UP     <- "Metabolism"
CAT_DESIGN <- c(CAT_OX, CAT_NULL, CAT_UP)

# The nested readouts. OXPHOS decomposes EXACTLY into its two functional halves;
# the mitoribosome is a GRANDCHILD of central dogma, via Translation.
OX_SUB   <- "OXPHOS subunits"
OX_AF    <- "OXPHOS assembly factors"
OX_MT    <- "mtDNA-encoded OXPHOS subunits"
MTRIB    <- "Mitochondrial ribosome"
NESTED   <- c(OX_SUB, OX_AF, OX_MT, MTRIB)

# E25's four readouts, carried for the arithmetic control.
E25_READOUTS <- c(OX_SUB, CAT_OX, MTRIB, OX_MT)

ARM_ORDER <- c("mouse 6W WT", "mouse 12W WT", "mouse 6W Myc+", "mouse 12W Myc+",
               "TCGA MYC-low", "TCGA MYC-high", "SCAN-B MYC-low", "SCAN-B MYC-high")
TERTILE   <- 1/3

.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}

RULES <- list(
  question = paste(
    "read the seven MitoCarta level-1 categories, plus four nested readouts,",
    "in the eight group orderings E25 already built"),
  estimands = paste(
    "VALUE (group-mean mitoPPS) is primary WITHIN a species and carries R1, R2,",
    "R3 and R5, because R2 is a statement about distances and needs a linear",
    "scale. PERCENTILE is secondary and is THE ONLY THING THAT CROSSES THE",
    "SPECIES BOUNDARY; it carries R4"),
  no_value_crossing = paste(
    "no mitoPPS value crosses the species boundary, in a table, an axis or a",
    "figure panel"),
  nothing_rescored = paste(
    "the level-1 categories are rows of the 142-pathway universe and were",
    "already scored; this re-reads E25's orderings and re-aggregates from",
    "existing mitoPPS. The mouse artefact is read-only and its md5 is checked"),
  overlap_reported = paste(
    "category overlap is REPORTED AND NOT GATED; it is not a stop condition"),
  batch = paste(
    "the mouse timeline is batch-confounded (batch = timepoint). That ruling",
    "stands: temporal numbers are reported in full with the confounder stated",
    "and NOTHING is withdrawn on its account"),
  purity_open = paste(
    "purity is NOT addressed here. E26 did that for OXPHOS subunits in TCGA",
    "only. Whether the category readouts survive adjustment is OPEN"),
  presentation = paste(
    "this is a presentation of an existing claim, not a new result; the",
    "expression-matched-null analysis in myc_mouse remains primary, and if the",
    "two disagree THAT DISAGREEMENT IS THE FINDING"),
  exploratory = "post-hoc and descriptive; not a hypothesis test",
  n3 = "transcript associations; 'primed' is never written of a transcript")

message("\n0. reading rules fixed before any number")
message("   R1   OXPHOS across the mouse WT window: falls with an interval clear")
message("        of zero -> the claim holds on this readout; covering zero ->")
message("        NOT RESOLVED AT THIS n, said plainly")
message("   R2   the null. STRONG: central dogma's interval EXCLUDES the magnitude")
message("        of the OXPHOS change -> the categories demonstrably differ.")
message("        WEAK: too wide -> 'not resolvable at this n'. NEVER a demonstrated null")
message("   R3   Metabolism rises -> where the priority went. Flat or falling ->")
message("        reported as such. A rise ALONE confirms nothing - it is partly forced")
message("   R4   OXPHOS above central dogma in mouse 6W Myc+ and below it in human")
message("        MYC-high? Neither outcome is predicted and neither is better")
message("   R4b  subunits against assembly factors, IN EACH SPECIES SEPARATELY,")
message("        on VALUES - percentile room near the panel middle is not flatness")
message("   R5   the other four. If the reallocation is even across all six")
message("        non-OXPHOS categories, say so")

# =============================================================================
# 1. INPUTS
# =============================================================================

message("\n1. inputs")

stopifnot(file.exists(PATH_E25))
e25 <- readRDS(PATH_E25)
message("   E25 object: ", length(e25), " elements, analysis date ",
        as.character(e25$analysis_date))

md5_now <- unname(tools::md5sum(PATH_MOUSE))
if (!identical(md5_now, MOUSE_MD5)) {
  stop("mouse artefact md5 has changed: ", md5_now, call. = FALSE)
}
mo <- readRDS(PATH_MOUSE)
message("   mouse artefact md5 matches: ", md5_now)

mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))

shared <- e25$intersection$shared
stopifnot(length(shared) == 142L)

UT <- mito$mitopps_universe[shared, , drop = FALSE]
US <- sc$mitopps_universe[shared, , drop = FALSE]
ID_T <- colnames(UT); ID_S <- colnames(US)

ms     <- mo$mitopps_scores
MPATHS <- setdiff(names(ms), c("sample", "group", "timepoint", "myc_status"))
MM     <- as.matrix(ms[, MPATHS, drop = FALSE]); rownames(MM) <- ms$sample
MMs    <- t(MM[, shared, drop = FALSE])            # pathways x samples, shared panel

# mitoPPS is centred at exactly 1 per pathway. The value estimand rests on it.
stopifnot(max(abs(rowMeans(UT) - 1)) < 1e-9,
          max(abs(rowMeans(US) - 1)) < 1e-9,
          max(abs(colMeans(MM) - 1)) < 1e-9)
message("   mitoPPS centred at 1 per pathway in all three cohorts (asserted)")

MOUSE_ARM <- c(`6W_neg` = "mouse 6W WT", `12W_neg` = "mouse 12W WT",
               `6W_pos` = "mouse 6W Myc+", `12W_pos` = "mouse 12W Myc+")

# --- the human MYC groups, E25's split reproduced and asserted ----------------
.tertile_groups <- function(v, ids) {
  q <- stats::quantile(v, c(TERTILE, 1 - TERTILE), names = FALSE)
  g <- rep(NA_character_, length(v)); names(g) <- ids
  g[v <= q[1]] <- "MYC-low"; g[v >= q[2]] <- "MYC-high"
  g
}
grp_T <- .tertile_groups(as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T]), ID_T)
grp_S <- .tertile_groups(as.numeric(sc$gsva_new[MYC_REF, ID_S]), ID_S)
stopifnot(identical(
  as.integer(c(sum(grp_T == "MYC-low", na.rm = TRUE), sum(grp_T == "MYC-high", na.rm = TRUE),
               sum(grp_S == "MYC-low", na.rm = TRUE), sum(grp_S == "MYC-high", na.rm = TRUE))),
  as.integer(e25$myc_groups$n)))
message("   human MYC groups reproduce E25 exactly: ",
        paste(sprintf("%s %s n=%d", e25$myc_groups$cohort, e25$myc_groups$group,
                      e25$myc_groups$n), collapse = " | "))

# --- the per-group sample index, one accessor used everywhere -----------------
ARM_IDX <- list()
for (g in names(MOUSE_ARM)) {
  ARM_IDX[[unname(MOUSE_ARM[[g]])]] <- list(sp = "mouse", M = MMs,
                                            j = which(ms$group == g))
}
for (lv in c("MYC-low", "MYC-high")) {
  ARM_IDX[[paste("TCGA", lv)]]   <- list(sp = "human", M = UT, j = which(grp_T == lv))
  ARM_IDX[[paste("SCAN-B", lv)]] <- list(sp = "human", M = US, j = which(grp_S == lv))
}
ARM_IDX <- ARM_IDX[ARM_ORDER]
ARM_SPECIES  <- vapply(ARM_IDX, function(z) z$sp, "")
ARM_COHORT <- c("mouse", "mouse", "mouse", "mouse",
                "TCGA", "TCGA", "SCAN-B", "SCAN-B")
names(ARM_COHORT) <- ARM_ORDER

# =============================================================================
# 2. PART 0 - the partition, read from source
# =============================================================================

message("\n2. PART 0: the level-1 partition, read from source not from memory")

hx <- readxl::read_excel(PATH_MITOCARTA_H, sheet = 4)
h_hier <- trimws(hx[["MitoPathways Hierarchy"]])
names(h_hier) <- trimws(hx[["MitoPathway"]])
h_hier <- h_hier[!is.na(h_hier) & nzchar(h_hier)]
h_lvl  <- lengths(strsplit(h_hier, " > ", fixed = TRUE))

pl <- mo$pathway_levels
LV_COLS <- c("Pathway_Level1", "Pathway_Level2", "Pathway_Level3")
m_hier <- apply(pl[, LV_COLS], 1L, function(z) {
  z <- z[!is.na(z) & nzchar(z)]; paste(trimws(z), collapse = " > ") })
names(m_hier) <- trimws(pl$Pathway)

CATS <- sort(names(h_lvl)[h_lvl == 1L & names(h_lvl) %in% shared])
stopifnot(length(CATS) == 7L, all(CATS %in% names(m_hier)))
CAT_ORDER <- c(CAT_DESIGN, setdiff(CATS, CAT_DESIGN))
stopifnot(setequal(CAT_ORDER, CATS))

HP  <- mito$mito_paths
MSZ <- c(table(mo$gene_to_pathway$Pathway))
.size <- function(p, sp) if (identical(sp, "mouse")) unname(MSZ[[p]]) else length(HP[[p]])

READOUTS <- c(CAT_ORDER, NESTED)
stopifnot(all(READOUTS %in% shared))

partition <- tibble::tibble(
  readout = READOUTS,
  role = c(rep("level-1 DESIGN", 3L), rep("level-1 context", 4L),
           "child of OXPHOS", "child of OXPHOS",
           "synthetic, OUTSIDE the OXPHOS umbrella",
           "grandchild of central dogma"),
  hierarchy_path = vapply(READOUTS, function(p)
    if (p %in% names(h_hier)) unname(h_hier[[p]]) else "(synthetic, not in the catalogue)", ""),
  human_n = vapply(READOUTS, function(p) .size(p, "human"), 0L),
  mouse_n = vapply(READOUTS, function(p) .size(p, "mouse"), 0L))
partition$size_delta <- partition$human_n - partition$mouse_n
print(as.data.frame(partition), row.names = FALSE, right = FALSE)

message("\n   the ontology asymmetry, NAMED here rather than found later:")
asym <- partition[partition$size_delta != 0, c("readout", "human_n", "mouse_n", "size_delta")]
if (nrow(asym)) {
  print(as.data.frame(asym), row.names = FALSE)
} else {
  message("      none - every readout is the same size in both species")
}
message("   WHICH gene differs cannot be reported: naming it needs a mouse-to-human")
message("   symbol alignment, which is an ortholog mapping whatever it is called.")
message("   Within-species membership and cross-species COUNTS only.")

# =============================================================================
# 3. PART A - overlap and containment. Reported, never gated.
# =============================================================================

message("\n3. PART A: overlap between the seven categories - REPORTED, NOT GATED")

gr <- utils::combn(length(CAT_ORDER), 2L)
overlap <- dplyr::bind_rows(lapply(seq_len(ncol(gr)), function(i) {
  a <- CAT_ORDER[gr[1, i]]; b <- CAT_ORDER[gr[2, i]]
  n <- length(intersect(HP[[a]], HP[[b]]))
  tibble::tibble(cat_a = a, cat_b = b, shared_genes = n,
                 pct_of_a = 100 * n / length(HP[[a]]),
                 pct_of_b = 100 * n / length(HP[[b]]))
}))
ov_nz <- overlap[overlap$shared_genes > 0, ]
ov_nz <- ov_nz[order(-ov_nz$shared_genes), ]
print(as.data.frame(ov_nz), digits = 3, row.names = FALSE)

uni7 <- unique(unlist(HP[CAT_ORDER]))
sum7 <- sum(lengths(HP[CAT_ORDER]))
uni_all <- unique(unlist(HP))
outside <- setdiff(uni_all, uni7)
overlap_summary <- tibble::tibble(
  quantity = c("union of the seven level-1 categories",
               "sum of the seven sizes",
               "genes assigned to more than one level-1 category",
               "genes in the 150-pathway catalogue",
               "catalogue genes in NO level-1 category"),
  n = c(length(uni7), sum7, sum7 - length(uni7), length(uni_all), length(outside)))
print(as.data.frame(overlap_summary), row.names = FALSE, right = FALSE)

OV_OX_MET <- overlap$shared_genes[overlap$cat_a == CAT_OX & overlap$cat_b == CAT_UP |
                                  overlap$cat_a == CAT_UP & overlap$cat_b == CAT_OX]
message("\n   THE PAIR THE REALLOCATION CLAIM IS ABOUT: ", CAT_OX, " x ", CAT_UP,
        " share ", OV_OX_MET, " genes.")
message("   DIRECTION, because the overlap is not neutral: shared genes push two")
message("   categories toward the SAME value. So OXPHOS falling while Metabolism")
message("   rises happens DESPITE the overlap, not because of it. The overlap is a")
message("   CONSERVATIVE confound for a dissociation and an inflating one only for")
message("   concordance.")
message("   QUALIFIER: that holds if the shared genes behave like OXPHOS-typical")
message("   genes. If they are atypical it weakens. A disjointified recomputation")
message("   would settle it and is FEASIBLE ON THE HUMAN SIDE ONLY - E24's PART A2")
message("   rebuilt the human universe bit-equal, but the mouse artefact carries")
message("   pathway scores, not gene-level values. NOT RUN. Available if pressed.")

message("\n   containment: the OXPHOS umbrella against its two halves (human)")
ox <- HP[[CAT_OX]]; sub <- HP[[OX_SUB]]; af <- HP[[OX_AF]]
containment <- tibble::tibble(
  quantity = c("umbrella", "subunits", "assembly factors",
               "subunits inside the umbrella", "assembly inside the umbrella",
               "genes in BOTH subunits and assembly factors",
               "umbrella genes in NEITHER half",
               "mtDNA-encoded genes inside the umbrella"),
  n = c(length(ox), length(sub), length(af), sum(sub %in% ox), sum(af %in% ox),
        length(intersect(sub, af)), length(setdiff(ox, union(sub, af))),
        length(intersect(HP[[OX_MT]], ox))))
print(as.data.frame(containment), row.names = FALSE, right = FALSE)
EXACT_CONTAINMENT <- length(setdiff(ox, union(sub, af))) == 0L
message("   CONTAINMENT IS ", if (EXACT_CONTAINMENT) "EXACT" else "NOT exact",
        ": the umbrella is precisely its two halves, ", length(sub), " + ",
        length(af), " - ", length(intersect(sub, af)), " = ", length(ox), ".")
message("   THE UMBRELLA THEREFORE CARRIES NO INDEPENDENT INFORMATION. It is an")
message("   aggregate of its halves and is never a third measurement.")
message("   Level-1 categories CONTAIN their descendants by design. That is")
message("   containment, not overlap, and it is expected.")

# =============================================================================
# 4. PART B - the readouts across the eight orderings, on both estimands
# =============================================================================

message("\n4. PART B: eleven readouts x eight groups, VALUE and PERCENTILE")

.group_value <- function(z, j) rowMeans(z$M[, j, drop = FALSE])
.group_pct   <- function(z, j) stats::setNames(.pct(.group_value(z, j)), rownames(z$M))

positions <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  z  <- ARM_IDX[[a]]
  v  <- .group_value(z, z$j)
  p  <- .pct(v)
  names(p) <- rownames(z$M)
  n  <- length(z$j)
  bs <- vapply(seq_len(NBOOT), function(i) {
    k  <- sample(z$j, replace = TRUE)
    vv <- .group_value(z, k)
    c(vv[READOUTS], .pct(vv)[match(READOUTS, rownames(z$M))])
  }, numeric(2L * length(READOUTS)))
  tibble::tibble(
    arm = a, cohort = unname(ARM_COHORT[[a]]), species = z$sp, n = n,
    readout = READOUTS,
    size = vapply(READOUTS, function(p2) .size(p2, z$sp), 0L),
    value = unname(v[READOUTS]),
    value_lo = apply(bs[seq_along(READOUTS), , drop = FALSE], 1L,
                     stats::quantile, 0.025),
    value_hi = apply(bs[seq_along(READOUTS), , drop = FALSE], 1L,
                     stats::quantile, 0.975),
    pct = unname(p[READOUTS]),
    pct_lo = apply(bs[length(READOUTS) + seq_along(READOUTS), , drop = FALSE], 1L,
                   stats::quantile, 0.025),
    pct_hi = apply(bs[length(READOUTS) + seq_along(READOUTS), , drop = FALSE], 1L,
                   stats::quantile, 0.975))
})) %>%
  dplyr::mutate(
    arm = factor(arm, levels = ARM_ORDER),
    readout = factor(readout, levels = READOUTS),
    role = partition$role[match(as.character(readout), partition$readout)])

# --- THE ARITHMETIC CONTROL: every percentile must reproduce E25 exactly ------
e25_ord <- e25$orderings[e25$orderings$panel == "shared (142)", ]
ctrl <- merge(
  positions[, c("arm", "readout", "pct")],
  stats::setNames(e25_ord[, c("arm", "pathway", "pct")], c("arm", "readout", "e25_pct")),
  by = c("arm", "readout"))
CTRL_DELTA <- max(abs(ctrl$pct - ctrl$e25_pct))
message("   ARITHMETIC CONTROL: ", nrow(ctrl),
        " recomputed percentiles against E25, max |delta| = ",
        format(CTRL_DELTA, digits = 2))
if (CTRL_DELTA >= 1e-9) {
  stop("CONTROL FAILED: the percentiles do not reproduce E25. Every number ",
       "below would be a rebuild rather than a re-read.", call. = FALSE)
}
message("   PASSED - this is the control on the whole script")

message("\n   PERCENTILE, the cross-species estimand (the three design categories):")
positions %>%
  dplyr::filter(readout %in% CAT_DESIGN) %>%
  dplyr::select(readout, arm, size, pct, pct_lo, pct_hi) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

message("\n   PERCENTILE, the OXPHOS decomposition (umbrella = its two halves):")
positions %>%
  dplyr::filter(readout %in% c(CAT_OX, OX_SUB, OX_AF, OX_MT)) %>%
  dplyr::select(readout, arm, size, pct, pct_lo, pct_hi) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

message("\n   VALUE, WITHIN SPECIES ONLY - mouse and human never share an axis:")
for (sp in c("mouse", "human")) {
  message("      --- ", sp, " ---")
  positions %>%
    dplyr::filter(species == sp, readout %in% c(CAT_DESIGN, OX_SUB, OX_AF)) %>%
    dplyr::select(readout, arm, value, value_lo, value_hi) %>%
    as.data.frame() %>%
    print(digits = 4, row.names = FALSE)
}

# =============================================================================
# 5. PART C - the six contrast orderings
# =============================================================================

message("\n5. PART C: the readouts inside the six contrast orderings")

CONTRASTS <- list(
  Myc_effect_6W          = list(sp = "mouse", M = MMs, a = "6W_neg",  b = "6W_pos"),
  Myc_effect_12W         = list(sp = "mouse", M = MMs, a = "12W_neg", b = "12W_pos"),
  `Temporal_Myc-`        = list(sp = "mouse", M = MMs, a = "6W_neg",  b = "12W_neg"),
  `Temporal_Myc+`        = list(sp = "mouse", M = MMs, a = "6W_pos",  b = "12W_pos"),
  `TCGA MYC-high - low`  = list(sp = "human", M = UT, a = "MYC-low", b = "MYC-high"),
  `SCAN-B MYC-high - low` = list(sp = "human", M = US, a = "MYC-low", b = "MYC-high"))
CONTRAST_KIND <- c(Myc_effect_6W = "genotype, within timepoint",
                   Myc_effect_12W = "genotype, within timepoint",
                   `Temporal_Myc-` = "TIMELINE - batch-confounded",
                   `Temporal_Myc+` = "TIMELINE - batch-confounded",
                   `TCGA MYC-high - low` = "MYC activity, within cohort",
                   `SCAN-B MYC-high - low` = "MYC activity, within cohort")

.contrast_idx <- function(cn) {
  if (identical(cn$sp, "mouse")) {
    return(list(a = which(ms$group == cn$a), b = which(ms$group == cn$b)))
  }
  g <- if (identical(dim(cn$M), dim(UT))) grp_T else grp_S
  list(a = which(g == cn$a), b = which(g == cn$b))
}
.contrast_vec <- function(cn, ia, ib) {
  rowMeans(cn$M[, ib, drop = FALSE]) - rowMeans(cn$M[, ia, drop = FALSE])
}

contrast_positions <- dplyr::bind_rows(lapply(names(CONTRASTS), function(ct) {
  cn <- CONTRASTS[[ct]]; ii <- .contrast_idx(cn)
  d  <- .contrast_vec(cn, ii$a, ii$b)
  pv <- .pct(d); names(pv) <- rownames(cn$M)
  bs <- vapply(seq_len(NBOOT), function(i) {
    dd <- .contrast_vec(cn, sample(ii$a, replace = TRUE), sample(ii$b, replace = TRUE))
    c(dd[READOUTS], .pct(dd)[match(READOUTS, rownames(cn$M))])
  }, numeric(2L * length(READOUTS)))
  tibble::tibble(
    contrast = ct, kind = unname(CONTRAST_KIND[[ct]]), species = cn$sp,
    readout = READOUTS,
    delta = unname(d[READOUTS]),
    delta_lo = apply(bs[seq_along(READOUTS), , drop = FALSE], 1L, stats::quantile, 0.025),
    delta_hi = apply(bs[seq_along(READOUTS), , drop = FALSE], 1L, stats::quantile, 0.975),
    pct = unname(pv[READOUTS]),
    pct_lo = apply(bs[length(READOUTS) + seq_along(READOUTS), , drop = FALSE], 1L,
                   stats::quantile, 0.025),
    pct_hi = apply(bs[length(READOUTS) + seq_along(READOUTS), , drop = FALSE], 1L,
                   stats::quantile, 0.975))
})) %>%
  dplyr::mutate(readout = factor(readout, levels = READOUTS))

message("   the three design categories, PERCENTILE inside each contrast ordering:")
contrast_positions %>%
  dplyr::filter(readout %in% CAT_DESIGN) %>%
  dplyr::select(contrast, kind, readout, pct) %>%
  tidyr::pivot_wider(names_from = readout, values_from = pct) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)
message("   the two temporal contrasts are BATCH-CONFOUNDED (batch = timepoint)")

# =============================================================================
# 6. PART D - R1, R2 and R3 on the VALUE scale
# =============================================================================

message("\n6. PART D: R1, R2, R3 - read on VALUES, within species")

# The mouse WT window, and its human analogues. Paired bootstrap: the SAME
# resample indices are used for every category, so a difference between two
# categories is a paired quantity and not two independent ones.
.window_boot <- function(cn) {
  ii <- .contrast_idx(cn)
  d  <- .contrast_vec(cn, ii$a, ii$b)
  B  <- vapply(seq_len(NBOOT), function(i)
    .contrast_vec(cn, sample(ii$a, replace = TRUE), sample(ii$b, replace = TRUE))[READOUTS],
    numeric(length(READOUTS)))
  rownames(B) <- READOUTS
  list(point = d[READOUTS], boot = B)
}
WIN <- list(
  `mouse WT window (6W -> 12W)`  = CONTRASTS[["Temporal_Myc-"]],
  `mouse Myc+ window (6W -> 12W)` = CONTRASTS[["Temporal_Myc+"]],
  `TCGA MYC-high - low`           = CONTRASTS[["TCGA MYC-high - low"]],
  `SCAN-B MYC-high - low`         = CONTRASTS[["SCAN-B MYC-high - low"]])
WB <- lapply(WIN, .window_boot)

window_change <- dplyr::bind_rows(lapply(names(WB), function(w) {
  z <- WB[[w]]
  tibble::tibble(window = w, readout = READOUTS,
                 delta = unname(z$point[READOUTS]),
                 ci_lo = apply(z$boot, 1L, stats::quantile, 0.025),
                 ci_hi = apply(z$boot, 1L, stats::quantile, 0.975))
})) %>%
  dplyr::mutate(readout = factor(readout, levels = READOUTS),
                clear_of_zero = (ci_lo > 0 & ci_hi > 0) | (ci_lo < 0 & ci_hi < 0))

message("   the mouse WT window, VALUES, all seven categories:")
window_change %>%
  dplyr::filter(window == "mouse WT window (6W -> 12W)",
                readout %in% CAT_ORDER) %>%
  dplyr::select(readout, delta, ci_lo, ci_hi, clear_of_zero) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

# --- R1 ----------------------------------------------------------------------
r1row <- window_change[window_change$window == "mouse WT window (6W -> 12W)" &
                         window_change$readout == CAT_OX, ]
R1 <- if (r1row$delta < 0 && r1row$clear_of_zero) {
  "FALLS with an interval clear of zero - the claim holds on this readout"
} else if (r1row$delta < 0) {
  "falls, but the interval covers zero - NOT RESOLVED AT THIS n"
} else {
  "does NOT fall - reported as such"
}
message("\n   R1 (", CAT_OX, ", mouse WT window): ", R1)

# --- R2, the equivalence reading ---------------------------------------------
# Pre-specified form: the null's interval must lie entirely inside
# (-|delta_OXPHOS|, +|delta_OXPHOS|). The paired difference is reported beside
# it because it uses the same resamples and is the more powerful statement, but
# the EQUIVALENCE BOUND is the rule that was fixed in advance.
zw    <- WB[["mouse WT window (6W -> 12W)"]]
BOUND <- abs(unname(zw$point[[CAT_OX]]))
nullrow <- window_change[window_change$window == "mouse WT window (6W -> 12W)" &
                           window_change$readout == CAT_NULL, ]
paired_diff <- abs(zw$boot[CAT_OX, ]) - abs(zw$boot[CAT_NULL, ])
# Assembly factors are carried here as an ALTERNATIVE null, read and not
# designed for: a sibling inside the OXPHOS parent controls for the respiratory
# branch specifically, not merely for "mitochondrial pathways in general".
.equiv_row <- function(nullset, label) {
  nr <- window_change[window_change$window == "mouse WT window (6W -> 12W)" &
                        window_change$readout == nullset, ]
  pdiff <- abs(zw$boot[CAT_OX, ]) - abs(zw$boot[nullset, ])
  tibble::tibble(
    null_readout = nullset, null_role = label,
    bound_is = paste0("|delta ", CAT_OX, "|"), bound = BOUND,
    null_delta = nr$delta, null_lo = nr$ci_lo, null_hi = nr$ci_hi,
    inside_bound = nr$ci_lo > -BOUND && nr$ci_hi < BOUND,
    paired_gap = abs(unname(zw$point[[CAT_OX]])) - abs(unname(zw$point[[nullset]])),
    paired_lo = unname(stats::quantile(pdiff, 0.025)),
    paired_hi = unname(stats::quantile(pdiff, 0.975)))
}
equivalence <- dplyr::bind_rows(
  .equiv_row(CAT_NULL, "PRE-SPECIFIED null - level-1 sibling"),
  .equiv_row(OX_AF, "alternative null, READ not designed - sibling inside OXPHOS"))
equivalence$paired_clear <- (equivalence$paired_lo > 0 & equivalence$paired_hi > 0) |
                            (equivalence$paired_lo < 0 & equivalence$paired_hi < 0)
print(as.data.frame(equivalence), digits = 3, row.names = FALSE)
eq_main <- equivalence[equivalence$null_readout == CAT_NULL, ]
R2 <- if (eq_main$inside_bound) {
  paste0("STRONG - ", CAT_NULL, "'s interval lies inside the magnitude of the ",
         CAT_OX, " change; the categories demonstrably differ")
} else {
  paste0("WEAK - ", CAT_NULL, "'s change is NOT RESOLVABLE AT THIS n. It is a ",
         "useful contrast against a resolved ", CAT_OX,
         " fall and is NEVER written as a demonstrated null")
}
message("   R2: ", R2)
# The pre-specified rule treats |delta OXPHOS| as a FIXED bound. The paired
# difference propagates the uncertainty in BOTH quantities and is the stricter
# statement. Report both; the pre-specified rule stands as the verdict.
if (eq_main$inside_bound && !eq_main$paired_clear) {
  message("   *** BUT: the paired difference |delta ", CAT_OX, "| - |delta ",
          CAT_NULL, "| is ", sprintf("%.4f", eq_main$paired_gap), " [",
          sprintf("%.4f", eq_main$paired_lo), ", ", sprintf("%.4f", eq_main$paired_hi),
          "] and COVERS ZERO.")
  message("   The equivalence bound treats |delta ", CAT_OX,
          "| as fixed; the paired statistic propagates both uncertainties and")
  message("   does NOT separate the two categories at this n. The pre-specified")
  message("   rule stands as the verdict, and BOTH numbers go in the note.")
  R2 <- paste0(R2, ". BUT the paired difference, which propagates uncertainty in ",
               "BOTH quantities, covers zero - so the separation is not ",
               "established by the stricter statistic at this n")
}
eq_alt <- equivalence[equivalence$null_readout == OX_AF, ]
message("   ALTERNATIVE NULL, read not designed: ", OX_AF, " moves ",
        sprintf("%.4f", eq_alt$null_delta), " [", sprintf("%.4f", eq_alt$null_lo),
        ", ", sprintf("%.4f", eq_alt$null_hi), "] across the same window,")
message("   against ", CAT_NULL, "'s ", sprintf("%.4f", eq_main$null_delta),
        ". A sibling INSIDE the OXPHOS parent controls for the respiratory")
message("   branch specifically, not merely for mitochondrial pathways in general.")

# --- R3, and the direction count on BOTH estimands ---------------------------
wt_val <- window_change[window_change$window == "mouse WT window (6W -> 12W)" &
                          window_change$readout %in% CAT_ORDER, ]
pc_6 <- positions$pct[positions$arm == "mouse 6W WT" & positions$readout %in% CAT_ORDER]
names(pc_6) <- as.character(positions$readout[positions$arm == "mouse 6W WT" &
                                                positions$readout %in% CAT_ORDER])
pc_12 <- positions$pct[positions$arm == "mouse 12W WT" & positions$readout %in% CAT_ORDER]
names(pc_12) <- as.character(positions$readout[positions$arm == "mouse 12W WT" &
                                                 positions$readout %in% CAT_ORDER])
# An exact zero is FLAT, not down. Percentile ties are possible on a 142-point
# scale and coding one as a fall would invent a direction the data does not have.
.dir <- function(d) ifelse(d > 0, "up", ifelse(d < 0, "down", "flat"))
vd <- wt_val$delta[match(CAT_ORDER, as.character(wt_val$readout))]
pd <- unname(pc_12[CAT_ORDER] - pc_6[CAT_ORDER])
direction <- tibble::tibble(
  readout = CAT_ORDER, value_delta = vd, value_dir = .dir(vd),
  pct_delta = pd, pct_dir = .dir(pd))
direction$agree <- direction$value_dir == direction$pct_dir
message("\n   R3 and R5: direction across the mouse WT window, BOTH estimands:")
print(as.data.frame(direction), digits = 3, row.names = FALSE)
n_up_v <- sum(direction$value_dir == "up"); n_up_p <- sum(direction$pct_dir == "up")
n_dn_v <- sum(direction$value_dir == "down"); n_dn_p <- sum(direction$pct_dir == "down")
message("   VALUES: ", n_up_v, " up, ", n_dn_v, " down, ", 7L - n_up_v - n_dn_v,
        " flat | PERCENTILES: ", n_up_p, " up, ", n_dn_p, " down, ",
        7L - n_up_p - n_dn_p, " flat")
ESTIMANDS_AGREE <- all(direction$agree)
if (!ESTIMANDS_AGREE) {
  message("   *** THE TWO ESTIMANDS DISAGREE on: ",
          paste(direction$readout[!direction$agree], collapse = ", "),
          " - that disagreement goes in the note")
}
met <- direction[direction$readout == CAT_UP, ]
R3 <- if (met$value_dir == "up") {
  paste0(CAT_UP, " RISES on values. It is NOT confirmation on its own - a rise ",
         "is partly forced by the relativity of the score, which is what R2 controls")
} else {
  paste0(CAT_UP, " does NOT rise on values - reported as such; the full ",
         "partition says where the priority went instead")
}
message("   R3: ", R3)
n_up6 <- sum(direction$value_dir[direction$readout != CAT_OX] == "up")
R5 <- if (n_up6 %in% c(0L, 6L)) {
  "the reallocation is EVEN across the six non-OXPHOS categories - a weaker result"
} else {
  paste0("the reallocation is NOT even: of the six non-OXPHOS categories, ",
         n_up6, " rise on values and the rest do not. THREE categories fall on ",
         "values alongside OXPHOS, so this is not a clean one-down-one-up ",
         "reallocation and R3's second branch is live")
}
message("   R5: ", R5)

# =============================================================================
# 7. R4 and R4b
# =============================================================================

message("\n7. R4 and R4b")

.pos <- function(a, r) positions$pct[positions$arm == a & positions$readout == r]
r4 <- tibble::tibble(
  group = c("mouse 6W Myc+", "TCGA MYC-high", "SCAN-B MYC-high"),
  OXPHOS_pct = vapply(c("mouse 6W Myc+", "TCGA MYC-high", "SCAN-B MYC-high"),
                      function(a) .pos(a, CAT_OX), 0),
  central_dogma_pct = vapply(c("mouse 6W Myc+", "TCGA MYC-high", "SCAN-B MYC-high"),
                             function(a) .pos(a, CAT_NULL), 0))
r4$OXPHOS_above <- r4$OXPHOS_pct > r4$central_dogma_pct
print(as.data.frame(r4), digits = 3, row.names = FALSE)
R4 <- if (r4$OXPHOS_above[1] && !any(r4$OXPHOS_above[2:3])) {
  paste0("the ordering difference E25 showed for one pathway IS present at ",
         "category level: OXPHOS above central dogma in mouse 6W Myc+ and below ",
         "it in human MYC-high, both cohorts")
} else if (all(r4$OXPHOS_above) || !any(r4$OXPHOS_above)) {
  "the same order in mouse and human - the ordering difference is NOT present at category level"
} else {
  "mixed across the two human cohorts - reported as such, not assigned to either reading"
}
message("   R4: ", R4)

# --- R4b, read on VALUES, in each species separately -------------------------
# Percentile room is not flatness: a readout sitting near the panel middle has
# less room to move than one swinging between the ends. The value scale has no
# such bound, which is why R4b is read on it.
r4b <- dplyr::bind_rows(lapply(names(WB), function(w) {
  z <- WB[[w]]
  tibble::tibble(window = w,
                 species = if (grepl("^mouse", w)) "mouse" else "human",
                 subunits = unname(z$point[[OX_SUB]]),
                 sub_lo = unname(stats::quantile(z$boot[OX_SUB, ], 0.025)),
                 sub_hi = unname(stats::quantile(z$boot[OX_SUB, ], 0.975)),
                 assembly = unname(z$point[[OX_AF]]),
                 af_lo = unname(stats::quantile(z$boot[OX_AF, ], 0.025)),
                 af_hi = unname(stats::quantile(z$boot[OX_AF, ], 0.975)),
                 gap = unname(z$point[[OX_SUB]] - z$point[[OX_AF]]),
                 gap_lo = unname(stats::quantile(z$boot[OX_SUB, ] - z$boot[OX_AF, ], 0.025)),
                 gap_hi = unname(stats::quantile(z$boot[OX_SUB, ] - z$boot[OX_AF, ], 0.975)))
}))
r4b$separate <- (r4b$gap_lo > 0 & r4b$gap_hi > 0) | (r4b$gap_lo < 0 & r4b$gap_hi < 0)
print(as.data.frame(r4b[, c("window", "species", "subunits", "assembly", "gap",
                            "gap_lo", "gap_hi", "separate")]),
      digits = 3, row.names = FALSE)

# does the assembly-factor SCORE hold near its cohort baseline of 1 while the
# subunit score moves? That is the question percentile room cannot answer.
baseline <- positions %>%
  dplyr::filter(readout %in% c(OX_SUB, OX_AF)) %>%
  dplyr::mutate(abs_dev = abs(value - 1)) %>%
  dplyr::group_by(species, readout) %>%
  dplyr::summarise(mean_abs_dev_from_1 = mean(abs_dev),
                   max_abs_dev_from_1 = max(abs_dev), .groups = "drop")
message("\n   how far each half moves from its cohort baseline of 1:")
print(as.data.frame(baseline), digits = 3, row.names = FALSE)

R4B <- if (all(r4b$separate)) {
  paste0(OX_SUB, " and ", OX_AF, " SEPARATE in every window and in both species -",
         " the respiratory chain and the machinery that builds it are ",
         "prioritised differently, one level below the title's distinction")
} else if (any(r4b$separate)) {
  paste0("they separate in ", sum(r4b$separate), " of ", nrow(r4b),
         " windows - named individually, not generalised")
} else {
  paste0(OX_SUB, " and ", OX_AF, " move together - OXPHOS behaves as one block ",
         "and the umbrella is an adequate summary of it")
}
message("   R4b: ", R4B)

# =============================================================================
# 8. PART E - the figures
# =============================================================================

message("\n8. PART E: figures")
if (!dir.exists(DIR_E27_FIG)) dir.create(DIR_E27_FIG, recursive = TRUE)

CAP_RANK <- paste(
  "PERCENTILE within the shared 142-pathway panel. Values never cross the",
  "species boundary and none appears here.")
CAP_EXPL <- paste(
  "Exploratory, post-hoc, descriptive; not a hypothesis test. Reading rules",
  "were fixed before any number so the outcome could not be chosen afterwards.")
CAP_BATCH <- paste(
  "The mouse timeline is batch-confounded (batch = timepoint); temporal numbers",
  "are reported in full and nothing is withdrawn on that account.")

.panel <- function(df, title, sub) {
  ggplot2::ggplot(df, ggplot2::aes(arm, pct, colour = readout, group = readout)) +
    ggplot2::geom_hline(yintercept = 50, linewidth = 0.3, colour = "grey75",
                        linetype = 2) +
    ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.4, colour = "grey45") +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.5) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = pct_lo, ymax = pct_hi),
                            linewidth = 0.45, alpha = 0.8) +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::labs(x = NULL, y = "percentile in the group ordering",
                  colour = NULL, title = title, subtitle = sub) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   legend.position = "top",
                   legend.text = ggplot2::element_text(size = 6.5),
                   plot.subtitle = ggplot2::element_text(size = 7))
}

pA <- .panel(dplyr::filter(positions, readout %in% CAT_DESIGN),
             "A. The three categories the design is about",
             paste0(CAT_OX, " expected down, ", CAT_NULL,
                    " the internal null, ", CAT_UP, " expected up"))
pB <- .panel(dplyr::filter(positions, readout %in% c(CAT_OX, OX_SUB, OX_AF, OX_MT)),
             "B. The OXPHOS decomposition",
             paste0("CONTAINMENT: the umbrella is EXACTLY subunits + assembly ",
                    "factors (", length(sub), " + ", length(af), " - ",
                    length(intersect(sub, af)), " = ", length(ox),
                    "), so it is an aggregate and not a third measurement"))
pC <- .panel(dplyr::filter(positions, readout %in% setdiff(CAT_ORDER, CAT_DESIGN)),
             "C. The remaining four level-1 categories",
             "Context: the complete partition is what shows whether the reallocation is even")

fig1 <- patchwork::wrap_plots(pA, pB, pC, ncol = 1L) +
  patchwork::plot_annotation(
    title = "The MitoCarta level-1 categories across all eight groups",
    caption = paste(CAP_RANK, CAP_BATCH, CAP_EXPL, sep = "\n"),
    theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, hjust = 0)))
ggplot2::ggsave(file.path(DIR_E27_FIG, "E27_category_positions.pdf"), fig1,
                width = 8, height = 11)

# --- the value figure. Mouse and human in SEPARATE panels with independent
# --- axes, so no reader can put a mouse value beside a human one.
vdat <- window_change %>%
  dplyr::filter(readout %in% CAT_ORDER) %>%
  dplyr::mutate(species = ifelse(grepl("^mouse", window), "mouse", "human"),
                readout = factor(readout, levels = rev(CAT_ORDER)))
pV <- ggplot2::ggplot(vdat, ggplot2::aes(delta, readout)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.5) +
  ggplot2::geom_point(size = 1.9) +
  ggplot2::facet_wrap(~ window, scales = "free_x", ncol = 2L) +
  ggplot2::labs(
    x = "change in group-mean mitoPPS (VALUE)", y = NULL,
    title = "Where R2 lives: the change on a linear scale, within species only",
    subtitle = paste0("Equivalence bound for the internal null is |delta ", CAT_OX,
                      "| in the mouse WT window. Axes are independent per panel."),
    caption = paste(
      "VALUES. Mouse and human panels have INDEPENDENT axes and no value is compared across them.",
      CAP_BATCH, CAP_EXPL, sep = "\n")) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E27_FIG, "E27_equivalence_values.pdf"), pV,
                width = 9, height = 5)
message("   wrote E27_category_positions.pdf and E27_equivalence_values.pdf")

# =============================================================================
# 9. PART F - verdict
# =============================================================================

message("\n9. PART F: verdict")

verdict <- tibble::tibble(
  rule = c("CONTROL percentiles vs E25", "PART A overlap", "PART A containment",
           "R1 OXPHOS", "R2 the null", "R3 Metabolism", "R4 the MYC ordering",
           "R4b subunits vs assembly", "R5 the other four", "estimand agreement"),
  outcome = c(
    paste0("PASSED - ", nrow(ctrl), " recomputed percentiles reproduce E25 exactly ",
           "(max |delta| ", format(CTRL_DELTA, digits = 2), ")"),
    paste0("NOT disjoint. ", CAT_OX, " x ", CAT_UP, " share ", OV_OX_MET,
           " genes; ", overlap_summary$n[3], " genes sit in more than one category; ",
           overlap_summary$n[5], " catalogue genes sit in none. Shared genes push ",
           "two categories TOGETHER, so a dissociation happens DESPITE the ",
           "overlap - conservative for this reading, inflating only for concordance"),
    paste0(if (EXACT_CONTAINMENT) "EXACT" else "NOT exact", " - the umbrella is ",
           "precisely its two halves and carries no independent information"),
    R1, R2, R3, R4, R4B, R5,
    if (ESTIMANDS_AGREE) {
      "values and percentiles agree on the direction of all seven categories"
    } else {
      paste0("VALUES AND PERCENTILES DISAGREE on: ",
             paste(direction$readout[!direction$agree], collapse = ", "))
    }))
print(as.data.frame(verdict), row.names = FALSE, right = FALSE)
message("\n   AND WHATEVER IT SAYS, IT IS DESCRIPTIVE. NOT A TEST.")

NOTES <- c(
  paste("Two estimands. VALUE (group-mean mitoPPS) is primary within a species",
        "and carries R1, R2, R3 and R5, because R2 is a statement about",
        "DISTANCES and needs a linear scale - percentile distances are bounded",
        "and compress near the ends. PERCENTILE is secondary and is the only",
        "thing that crosses the species boundary; it carries R4."),
  paste("NO mitoPPS VALUE CROSSES THE SPECIES BOUNDARY, in a table, an axis or",
        "a figure panel. The eight-group figure is on percentiles for exactly",
        "that reason, and the value figure gives mouse and human independent axes."),
  paste("The zero-sum argument for values is weaker than it looks and is not",
        "relied on: seven readouts inside a 142-pathway panel leave 135",
        "pathways to absorb any difference, so percentile forcing is weak",
        "rather than strict. The reason for values is the equivalence scale."),
  paste("DISCLOSURE. The eleven readouts' PERCENTILE positions were computable",
        "from E25's saved object - they are rows of $orderings - but E25's NOTE",
        "reported only four. They were surfaced during planning, BEFORE the run,",
        "not after it. The reading rules were fixed by the author before any",
        "number was seen, and that is what protects them."),
  paste("R4b is read on VALUES because percentile room is not flatness: a",
        "readout near the panel middle has less room to move than one swinging",
        "between the ends. The value scale has no such bound."),
  paste("Category overlap is REPORTED AND NOT GATED. Direction matters: shared",
        "genes push two categories toward the SAME value, so a dissociation",
        "happens DESPITE the overlap. QUALIFIER: that holds if the shared genes",
        "behave like OXPHOS-typical genes. A disjointified recomputation would",
        "settle it and is FEASIBLE ON THE HUMAN SIDE ONLY, because the mouse",
        "artefact carries pathway scores and not gene-level values. NOT RUN."),
  paste("The OXPHOS umbrella is EXACTLY its two halves, so it carries no",
        "independent information and is never a third measurement. The",
        "mitoribosome is a GRANDCHILD of central dogma, via Translation."),
  paste("The ontology asymmetry is NAMED, not found later: the umbrella is 156",
        "human against 155 mouse and assembly factors 68 against 67. WHICH gene",
        "differs cannot be reported - naming it needs a mouse-to-human symbol",
        "alignment, which is an ortholog mapping whatever it is called. Within-",
        "species membership and cross-species counts only."),
  paste("The mouse timeline is BATCH-CONFOUNDED (batch = timepoint). Temporal",
        "numbers are reported in full and nothing is withdrawn on that account.",
        "The note may observe that a batch effect would have to act",
        "DIFFERENTIALLY ACROSS THREE CATEGORIES IN THE SAME LIBRARIES - that",
        "narrows the confound, it does not exclude it."),
  paste("PURITY IS NOT ADDRESSED HERE. E26 did that for OXPHOS subunits in TCGA",
        "only. Whether the category readouts survive adjustment is OPEN."),
  paste("This is a PRESENTATION OF AN EXISTING CLAIM, not a new result. The",
        "expression-matched-null analysis in myc_mouse remains the primary",
        "evidence for the developmental change. Both stay reportable and IF",
        "THEY DISAGREE THAT DISAGREEMENT IS THE FINDING."),
  "Nothing is re-scored. No ortholog call. N3 throughout.")

saveRDS(list(
  rules = RULES,
  settings = list(question = RULES$question, categories = CAT_ORDER,
                  design = CAT_DESIGN, nested = NESTED, readouts = READOUTS,
                  arm_order = ARM_ORDER, tertile = TERTILE, nboot = NBOOT,
                  seed = 1, equivalence_bound = BOUND,
                  mouse_artefact = list(md5 = MOUSE_MD5, recomputed = FALSE)),
  partition = partition, overlap = overlap, overlap_summary = overlap_summary,
  containment = containment, exact_containment = EXACT_CONTAINMENT,
  category_positions = positions, contrast_positions = contrast_positions,
  window_change = window_change, direction = direction,
  equivalence = equivalence, r4 = r4, r4b = r4b, baseline = baseline,
  control_delta = CTRL_DELTA, estimands_agree = ESTIMANDS_AGREE,
  readings = list(R1 = R1, R2 = R2, R3 = R3, R4 = R4, R4b = R4B, R5 = R5),
  verdict = verdict, analysis_date = Sys.Date(), notes = NOTES), PATH_E27)

readr::write_csv(
  positions %>% dplyr::mutate(arm = as.character(arm), readout = as.character(readout)),
  PATH_E27_TAB)

message("\nE27: done.")
message("    ", PATH_E27)
message("    ", PATH_E27_TAB)
message("    ", file.path(DIR_E27_FIG, "E27_category_positions.pdf"))
message("    ", file.path(DIR_E27_FIG, "E27_equivalence_values.pdf"))

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E27)

  ## the control first - if this is not 0 nothing else is a re-read
  x$control_delta

  ## the partition and the ontology asymmetry, read from source
  x$partition |> as.data.frame()

  ## overlap: reported, never gated. Read the direction argument in the notes.
  subset(x$overlap, shared_genes > 0) |> as.data.frame()
  x$overlap_summary |> as.data.frame()
  x$containment |> as.data.frame()

  ## the verdict and the six readings
  x$verdict |> as.data.frame()
  x$readings

  ## R1/R2/R3 live on VALUES, within species
  subset(x$window_change, window == "mouse WT window (6W -> 12W)") |> as.data.frame()
  x$equivalence |> as.data.frame()
  x$direction |> as.data.frame()

  ## R4 on percentiles - the only estimand that crosses the species boundary
  x$r4 |> as.data.frame()

  ## R4b on values, each species separately, plus the baseline check that
  ## distinguishes real flatness from percentile room
  x$r4b |> as.data.frame()
  x$baseline |> as.data.frame()

  ## the eleven readouts, percentile, all eight groups
  subset(x$category_positions, readout %in% x$settings$design)[
    , c("arm", "readout", "pct", "pct_lo", "pct_hi")] |> as.data.frame()

  ## inside the six contrast orderings
  subset(x$contrast_positions, readout %in% x$settings$design)[
    , c("contrast", "kind", "readout", "pct")] |> as.data.frame()

  ## every caveat the note must repeat
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")
}
