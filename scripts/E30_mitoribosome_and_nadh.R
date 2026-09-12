# scripts/E30_mitoribosome_and_nadh.R
# =============================================================================
# E30 - the mitoribosome gap, and the NADH candidate panel
# =============================================================================
#
# TWO PARTS WITH DIFFERENT STANDING. THEY DO NOT MERGE.
#
#   PART A is a DECLARED CHECK of an existing finding on a second ruler. E20 and
#   E21 established on a LEVEL ruler that MYC activity tracks the MITORIBOSOME
#   more than it tracks OXPHOS. PART A asks whether the same ordering appears on
#   the RANK ruler. The comparator is named in advance and is not chosen from
#   the panel.
#
#   PART B is a DESCRIPTIVE PANEL AND NOTHING MORE. The annotation is the
#   author's, it is ambiguous and contested, and IT IS EXPLICITLY NOT A
#   PRE-DECLARED RULE. No hypothesis is tested, no set membership is defended,
#   and nothing from PART B may be promoted to a readout without a separately
#   declared validation.
#
# EXPLORATORY, POST-HOC, DESCRIPTIVE throughout, as all of Phase 2 is.
#
# FIVE THINGS ESTABLISHED BEFORE THIS SCRIPT WAS WRITTEN, ALL OF WHICH CHANGE
# HOW IT MUST BE READ
# -----------------------------------------------------------------------------
# 1. PART A IS A REFINEMENT OF THE STATISTIC THAT FAILED, NOT AN INDEPENDENT ONE.
#    The mitoribosome is 100 per cent contained in Mitochondrial central dogma -
#    all 83 of its genes - so the new denominator sits INSIDE the denominator
#    E28's gap readout failed with. Across the eight groups the two gaps
#    correlate at rho +0.905. The PRIOR is independent (E20/E21 named the
#    mitoribosome in advance, on a different ruler); the STATISTIC is nested.
#    CONSEQUENCE FOR R1: if this gap also fails the ladder that is close to
#    expected and adds little. THE INFORMATIVE CASE IS IF IT PASSES WHERE THE
#    PARENT FAILED.
#
# 2. THE NUMERATOR IS BETTER HERE THAN IN E28, and that is not a variation.
#    E28's gap used the OXPHOS UMBRELLA. This uses OXPHOS SUBUNITS. E27 showed
#    the umbrella is exactly subunits plus assembly factors with assembly
#    factors barely moving, and E28's own aggregate measured the umbrella
#    DILUTING by 8.5 percentile points. So both terms change and one change is
#    an improvement.
#
# 3. E20's FINDING IS ALREADY ESTIMATOR-DEPENDENT, AND R1 IS BUILT KNOWING IT.
#    E20: the THREE SIGNATURE estimators track the mitoribosome more than OXPHOS
#    (6 of 6). THE REGULON IS THE ONE ESTIMATOR THAT DOES NOT (+0.004 / +0.003,
#    both intervals covering zero). E21 explained that dissent as a SIZE
#    artefact - the regulon sits inside its own size-matched null. So a regulon
#    dissent on the rank ruler is CONSISTENT WITH THE PRIOR and R1 exempts it,
#    with the reason stated rather than the exemption hidden.
#    E20 also already built this comparator on the level ruler as
#    `ox_nuc_mtrib = comp(OXPHOS subunits) - comp(mitoribosome)`, where the
#    adjudicating comparison came out 8 OF 8 INCLUDING THE REGULON. That is the
#    concrete number PART A is checked against.
#
# 4. THE "TWO TIERS" ARE EXACT TO THE RANK STEP, AND THAT IS A COINCIDENCE.
#    Both mouse windows land on -54 steps from DIFFERENT components (-76 and -22
#    against -80 and -26). Both within-timepoint MYC contrasts land on -24, from
#    +40 and +64 against +36 and +60. Four contrasts hitting two integers makes
#    the tier story look structural when the components say it is not. The
#    components are reported IN RANK STEPS so the exactness cannot be read as
#    meaning. R2's tier comparison is run WITHIN THE MOUSE, which avoids a
#    cross-species magnitude comparison in percentile space.
#
# 5. ONE PART B ENTRY IS FULLY TAUTOLOGICAL AND TWO ARE PARTLY SO.
#    `CI subunits` is 100 per cent contained in OXPHOS subunits, so a gap
#    between them is ARITHMETIC. `Complex I` shares 63 per cent, `TCA cycle` 20,
#    `Electron carriers` 14; everything else shares zero, INCLUDING
#    `CI assembly factors`. Gene sharing is marked on every PART B row and
#    panel, and `CI subunits` is kept but greyed with its containment stated.
#
# WHAT E28 ESTABLISHED ABOUT GAPS, AND WHAT IS DONE ABOUT IT HERE
# ----------------------------------------------------------------
#   A GAP DISCARDS THE THING IT IS BUILT FROM. E28: mouse 12W WT and TCGA
#   MYC-high gave an identical central-dogma gap from OXPHOS positions of 14 and
#   62. THEREFORE EVERY GAP IN THIS SCRIPT IS REPORTED WITH BOTH COMPONENT
#   POSITIONS BESIDE IT, in every table and every figure. NEVER A GAP ALONE.
#
#   THE GAP IS NOT PROTECTED AGAINST COMPOSITION. E28 measured it moving 2.1
#   percentile points against 1.4 and 0.7 for its components - differencing ADDS
#   rather than cancels. That check is rerun here for the mitoribosome gap.
#   TCGA ONLY: SCAN-B has no purity estimate and nothing is imputed.
#
#   E28's GAP READOUT FAILED AND IS NOT REHABILITATED BY THIS SCRIPT. If PART A
#   fails its ladder, it joins E28 as a negative.
#
# NO mitoPPS VALUE CROSSES THE SPECIES BOUNDARY - percentiles only. TCGA and
# SCAN-B are never averaged. The mouse timeline is BATCH-CONFOUNDED (batch =
# timepoint): reported in full, confounder stated, nothing withdrawn. Nothing is
# re-scored. No ortholog call. No FDR. N3 throughout.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))

message("\nE30: the mitoribosome gap, and the NADH candidate panel\n", strrep("=", 78))

PATH_E30      <- file.path(DIR_RESULTS, "mitoribosome_and_nadh.rds")
PATH_E30_A    <- file.path(DIR_TABLES,  "E30_mitoribosome_gap.csv")
PATH_E30_B    <- file.path(DIR_TABLES,  "E30_nadh_panel.csv")
DIR_E30_FIG   <- file.path(DIR_OUTPUTS, "mitopps_rank")
PATH_E25 <- file.path(DIR_RESULTS, "mitopps_group_rankings.rds")
PATH_E27 <- file.path(DIR_RESULTS, "category_readouts.rds")
PATH_E28 <- file.path(DIR_RESULTS, "gap_readout.rds")
PATH_E29 <- file.path(DIR_RESULTS, "oxphos_neighbourhood.rds")
PATH_E26 <- file.path(DIR_RESULTS, "purity_and_composition.rds")
PATH_MOUSE <- here::here("data", "from_myc_mouse", "mitopps_scores.rds")

set.seed(1)
NBOOT <- 5000L

# =============================================================================
# 0. CONSTANTS AND THE READING RULES, FIXED BEFORE ANY NUMBER
# =============================================================================

MOUSE_MD5 <- "b8a125af4bc0d5909ad03f6e126e2890"
REF   <- "OXPHOS subunits"              # the numerator, in BOTH parts
MTRIB <- "Mitochondrial ribosome"       # PART A's declared comparator
CD    <- "Mitochondrial central dogma"  # E28's comparator, for comparison only
TERTILE <- 1/3

ARM_ORDER <- c("mouse 6W WT", "mouse 12W WT", "mouse 6W Myc+", "mouse 12W Myc+",
               "TCGA MYC-low", "TCGA MYC-high", "SCAN-B MYC-low", "SCAN-B MYC-high")

# The author's PART B candidate set. PROVISIONAL ANNOTATION, NOT A CLASSIFICATION.
NADH <- list(
  `plausible matrix NADH producers` = c(
    "TCA cycle", "TCA-associated", "Pyruvate metabolism", "Fatty acid oxidation",
    "Glycine cleavage system", "Branched-chain amino acid dehydrogenase complex",
    "Glutamate metabolism", "Ketone metabolism", "Propanoate metabolism"),
  `shuttle / pool` = c(
    "Malate-aspartate shuttle", "NAD biosynthesis and metabolism", "Electron carriers"),
  `the consumer side` = c("Complex I", "CI subunits", "CI assembly factors"),
  `ambiguous, named as such` = c(
    "Folate and 1-C metabolism", "Sulfur metabolism", "Coenzyme Q metabolism"))
NADH_NOTE <- c(
  `Folate and 1-C metabolism` = "runs to NADPH as often as NADH",
  `Sulfur metabolism` = "ambiguous",
  `Coenzyme Q metabolism` = "ambiguous")

# The prompt's supplied point estimates, asserted rather than trusted.
SUPPLIED <- c(`WT window (6W -> 12W)` = -38.0, `Myc+ window (6W -> 12W)` = -38.0,
              `MYC at 6W` = -16.9, `MYC at 12W` = -16.9,
              `TCGA MYC-high - low` = -17.6, `SCAN-B MYC-high - low` = -21.1)
SUPPLIED_TOL <- 0.05

# E20's level-ruler number, quoted so R1's agreement branch is concrete.
E20_LEVEL <- paste(
  "E20, LEVEL ruler: rho(ox_lvl) - rho(ox_nuc_mtrib) positive in 8 of 8",
  "adjudicating cells, both cohorts, raw and adjusted, INCLUDING the regulon.",
  "Separately, mtrib_lvl - ox_lvl positive 6 of 6 for the three SIGNATURE",
  "estimators and NOT for the regulon (+0.004 / +0.003, intervals covering 0),",
  "which E21 explained as a SIZE artefact.")

RULES <- list(
  standing = paste(
    "PART A is a DECLARED CHECK of an existing finding on a second ruler.",
    "PART B is a DESCRIPTIVE PANEL on the author's PROVISIONAL annotation and",
    "is not a classification, not a hypothesis and not promotable to a readout"),
  gap_components = paste(
    "EVERY gap is reported with BOTH component positions, in every table and",
    "every figure. A gap discards the thing it is built from (E28)"),
  nested = paste(
    "PART A's comparator is 100 per cent CONTAINED in E28's failed comparator",
    "and the two gaps correlate at rho +0.905. The prior is independent; the",
    "statistic is not. The informative case is PASSING where the parent FAILED"),
  numerator = paste(
    "the numerator is OXPHOS SUBUNITS here and was the UMBRELLA in E28, which",
    "E28 measured as diluting by 8.5 percentile points"),
  regulon = paste(
    "R1 EXEMPTS myc_regulon: E20 already recorded it as the one estimator that",
    "does not show the level-ruler effect, and E21 explained that as size"),
  tiers = paste(
    "the two tiers are exact to the RANK STEP from different components - a",
    "coincidence, not structure. Components are reported in steps"),
  tautology = paste(
    "CI subunits is 100 per cent contained in the numerator so its gap is",
    "ARITHMETIC; gene sharing is marked on every PART B row"),
  no_value_crossing = "no mitoPPS value crosses the species boundary; percentiles only",
  batch = "the mouse timeline is batch-confounded; nothing withdrawn on that account",
  not_rehabilitated = paste(
    "E28's gap readout FAILED and is not rehabilitated here. If PART A fails",
    "its ladder it joins E28 as a negative"),
  exploratory = "post-hoc and descriptive; not a hypothesis test",
  n3 = "transcript associations; 'primed' is never written of a transcript")

message("\n0. reading rules fixed before any number")
message("   R1  the ladder. Sign held across the estimators AND agreeing with the")
message("       COPY-NUMBER split -> the rank ruler corroborates E20/E21. Sign flips")
message("       or copy number disagrees -> it FAILS AS E28's DID, which is itself")
message("       informative: the failure is then a property of the RANK-GAP")
message("       STATISTIC rather than of the comparator. myc_regulon is EXEMPT")
message("       because E20/E21 already recorded and explained its dissent")
message("   R2  do the tiers separate WITH INTERVALS, within the mouse? If not, say")
message("       the separation is not established rather than describing two tiers")
message("   R3  the mouse contrasts. Covering zero does not refute the point estimate;")
message("       it means the mouse side CANNOT STAND ALONE and must be read beside")
message("       the level-ruler evidence in E20/E21")
message("   R4  composition. Gap moves LESS than its components -> differencing")
message("       protects here where it did not for central dogma. As much or more ->")
message("       it does not, and every gap number inherits that")
message("   R5  the NADH panel is DESCRIPTIVE. No set-level statistic, no counting")
message("       members by direction as though the count were a result, no programme")
message("       named. E29 already found the structure diffuse")

.pct <- function(v) {
  stopifnot(!anyNA(v))
  100 * (rank(v, ties.method = "average") - 0.5) / length(v)
}

# =============================================================================
# 1. INPUTS
# =============================================================================

message("\n1. inputs")

e25 <- readRDS(PATH_E25); e26 <- readRDS(PATH_E26); e27 <- readRDS(PATH_E27)
e28 <- readRDS(PATH_E28); e29 <- readRDS(PATH_E29)
stopifnot(identical(unname(tools::md5sum(PATH_MOUSE)), MOUSE_MD5))
mo   <- readRDS(PATH_MOUSE)
mito <- readRDS(PATH_TCGA_MITO)
sc   <- readRDS(file.path(DIR_RESULTS, "scanb_scores.rds"))
nw   <- readRDS(file.path(DIR_RESULTS, "new_set_scores.rds"))
cvd  <- readRDS(here::here("data", "from_validation",
                          "tcga_brca_covariates.rds"))$covariates
message("   E25 ", length(e25), " | E26 ", length(e26), " | E27 ", length(e27),
        " | E28 ", length(e28), " | E29 ", length(e29), " elements")

shared <- e25$intersection$shared
P_ALL  <- length(shared)
STEP   <- 100 / P_ALL
HP     <- mito$mito_paths
stopifnot(all(c(REF, MTRIB, CD) %in% shared))

UT <- mito$mitopps_universe[shared, , drop = FALSE]
US <- sc$mitopps_universe[shared, , drop = FALSE]
ID_T <- colnames(UT); ID_S <- colnames(US)
ms <- mo$mitopps_scores
MPATHS <- setdiff(names(ms), c("sample", "group", "timepoint", "myc_status"))
MM <- as.matrix(ms[, MPATHS, drop = FALSE]); rownames(MM) <- ms$sample
MMs <- t(MM[, shared, drop = FALSE])

MOUSE_ARM <- c(`6W_neg` = "mouse 6W WT", `12W_neg` = "mouse 12W WT",
               `6W_pos` = "mouse 6W Myc+", `12W_pos` = "mouse 12W Myc+")
.tertiles <- function(v, ids) {
  q <- stats::quantile(v, c(TERTILE, 1 - TERTILE), names = FALSE)
  g <- rep(NA_character_, length(v)); names(g) <- ids
  g[v <= q[1]] <- "MYC-low"; g[v >= q[2]] <- "MYC-high"; g
}
est_T <- as.numeric(nw$tcga_gsva_new[MYC_REF, ID_T])
est_S <- as.numeric(sc$gsva_new[MYC_REF, ID_S])
grp_T <- .tertiles(est_T, ID_T); grp_S <- .tertiles(est_S, ID_S)
stopifnot(identical(as.integer(c(sum(grp_T == "MYC-low", na.rm = TRUE),
                                 sum(grp_T == "MYC-high", na.rm = TRUE),
                                 sum(grp_S == "MYC-low", na.rm = TRUE),
                                 sum(grp_S == "MYC-high", na.rm = TRUE))),
                    as.integer(e25$myc_groups$n)))

ARM <- list()
for (g in names(MOUSE_ARM)) {
  ARM[[unname(MOUSE_ARM[[g]])]] <- list(sp = "mouse", M = MMs, j = which(ms$group == g))
}
for (lv in c("MYC-low", "MYC-high")) {
  ARM[[paste("TCGA", lv)]]   <- list(sp = "human", M = UT, j = which(grp_T == lv))
  ARM[[paste("SCAN-B", lv)]] <- list(sp = "human", M = US, j = which(grp_S == lv))
}
ARM <- ARM[ARM_ORDER]
.gp <- function(z, j) {
  v <- rowMeans(z$M[, j, drop = FALSE])
  stats::setNames(.pct(v), rownames(z$M))
}

CONTRASTS <- list(
  `WT window (6W -> 12W)`   = c(a = "mouse 6W WT",    b = "mouse 12W WT"),
  `Myc+ window (6W -> 12W)` = c(a = "mouse 6W Myc+",  b = "mouse 12W Myc+"),
  `MYC at 6W`               = c(a = "mouse 6W WT",    b = "mouse 6W Myc+"),
  `MYC at 12W`              = c(a = "mouse 12W WT",   b = "mouse 12W Myc+"),
  `TCGA MYC-high - low`     = c(a = "TCGA MYC-low",   b = "TCGA MYC-high"),
  `SCAN-B MYC-high - low`   = c(a = "SCAN-B MYC-low", b = "SCAN-B MYC-high"))
CONTRAST_KIND <- c(`WT window (6W -> 12W)` = "TIMELINE - batch-confounded",
                   `Myc+ window (6W -> 12W)` = "TIMELINE - batch-confounded",
                   `MYC at 6W` = "genotype, within timepoint",
                   `MYC at 12W` = "genotype, within timepoint",
                   `TCGA MYC-high - low` = "MYC activity, within cohort",
                   `SCAN-B MYC-high - low` = "MYC activity, within cohort")

# =============================================================================
# 2. PART 0 - reproduce, and the nesting this script sits inside
# =============================================================================

message("\n2. PART 0: reproduce")

pos <- vapply(ARM_ORDER, function(a) .gp(ARM[[a]], ARM[[a]]$j), numeric(P_ALL))
rownames(pos) <- shared
e25o <- e25$orderings[e25$orderings$panel == "shared (142)", ]
c25 <- max(abs(pos[cbind(match(as.character(e25o$pathway), shared),
                         match(as.character(e25o$arm), ARM_ORDER))] - e25o$pct))
e27p <- e27$category_positions
c27 <- max(abs(pos[cbind(match(as.character(e27p$readout), shared),
                         match(as.character(e27p$arm), ARM_ORDER))] - e27p$pct))
message("   ARITHMETIC CONTROL: ", nrow(e25o), " against E25 max |delta| ",
        format(c25, digits = 2), "; ", nrow(e27p), " against E27 max |delta| ",
        format(c27, digits = 2))
if (max(c25, c27) >= 1e-9) stop("CONTROL FAILED - stop.", call. = FALSE)

message("\n   THE NESTING, measured rather than asserted:")
nest <- tibble::tibble(
  quantity = c(paste(MTRIB, "genes"), paste(CD, "genes"),
               "shared", "pct of the mitoribosome inside central dogma"),
  value = c(length(HP[[MTRIB]]), length(HP[[CD]]),
            length(intersect(HP[[MTRIB]], HP[[CD]])),
            round(100 * length(intersect(HP[[MTRIB]], HP[[CD]])) / length(HP[[MTRIB]]))))
print(as.data.frame(nest), row.names = FALSE, right = FALSE)
gap_mr <- pos[REF, ] - pos[MTRIB, ]
gap_cd_sub <- pos[REF, ] - pos[CD, ]
gap_cd_umb <- pos["OXPHOS", ] - pos[CD, ]           # E28 AS RUN
NEST_RHO <- stats::cor(gap_mr, gap_cd_umb, method = "spearman")
message("   this gap against E28's AS RUN (umbrella - central dogma): rho ",
        sprintf("%+.3f", NEST_RHO))
message("   against the like-for-like (subunits - central dogma):      rho ",
        sprintf("%+.3f", stats::cor(gap_mr, gap_cd_sub, method = "spearman")))
message("   PART A is a REFINEMENT of the statistic E28 failed with, not an")
message("   independent one. The informative case is PASSING where it FAILED.")

# =============================================================================
# 3. PART A - the mitoribosome gap
# =============================================================================

message("\n3. PART A: the mitoribosome gap. BOTH COMPONENTS ON EVERY ROW.")

.boot_group <- function(a, paths) {
  z <- ARM[[a]]
  vapply(seq_len(NBOOT), function(i) .gp(z, sample(z$j, replace = TRUE))[paths],
         numeric(length(paths)))
}
BG <- lapply(ARM_ORDER, function(a) .boot_group(a, c(REF, MTRIB)))
names(BG) <- ARM_ORDER

mitorib_by_group <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  b <- BG[[a]]; g <- b[REF, ] - b[MTRIB, ]
  tibble::tibble(
    arm = a, species = ARM[[a]]$sp, n = length(ARM[[a]]$j),
    oxphos_pct = pos[REF, a],
    ox_lo = unname(stats::quantile(b[REF, ], 0.025)),
    ox_hi = unname(stats::quantile(b[REF, ], 0.975)),
    mitorib_pct = pos[MTRIB, a],
    mt_lo = unname(stats::quantile(b[MTRIB, ], 0.025)),
    mt_hi = unname(stats::quantile(b[MTRIB, ], 0.975)),
    gap = pos[REF, a] - pos[MTRIB, a],
    gap_lo = unname(stats::quantile(g, 0.025)),
    gap_hi = unname(stats::quantile(g, 0.975)),
    gap_steps = round((pos[REF, a] - pos[MTRIB, a]) / STEP))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER),
                gap_clear = (gap_lo > 0 & gap_hi > 0) | (gap_lo < 0 & gap_hi < 0))
message("   A1: positions and gap per group")
mitorib_by_group %>%
  dplyr::select(arm, n, oxphos_pct, mitorib_pct, gap, gap_lo, gap_hi,
                gap_steps, gap_clear) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

.boot_contrast <- function(ct, paths) {
  cn <- CONTRASTS[[ct]]; za <- ARM[[cn[["a"]]]]; zb <- ARM[[cn[["b"]]]]
  vapply(seq_len(NBOOT), function(i) {
    pa <- .gp(za, sample(za$j, replace = TRUE))[paths]
    pb <- .gp(zb, sample(zb$j, replace = TRUE))[paths]
    pb - pa
  }, numeric(length(paths)))
}
BC <- lapply(names(CONTRASTS), function(ct) .boot_contrast(ct, c(REF, MTRIB)))
names(BC) <- names(CONTRASTS)

mitorib_by_contrast <- dplyr::bind_rows(lapply(names(CONTRASTS), function(ct) {
  cn <- CONTRASTS[[ct]]; b <- BC[[ct]]; g <- b[REF, ] - b[MTRIB, ]
  d_ox <- pos[REF, cn[["b"]]] - pos[REF, cn[["a"]]]
  d_mt <- pos[MTRIB, cn[["b"]]] - pos[MTRIB, cn[["a"]]]
  tibble::tibble(
    contrast = ct, kind = unname(CONTRAST_KIND[[ct]]),
    species = ARM[[cn[["a"]]]]$sp,
    d_oxphos = d_ox, d_oxphos_steps = round(d_ox / STEP),
    d_mitorib = d_mt, d_mitorib_steps = round(d_mt / STEP),
    gap_change = d_ox - d_mt, gap_change_steps = round((d_ox - d_mt) / STEP),
    ci_lo = unname(stats::quantile(g, 0.025)),
    ci_hi = unname(stats::quantile(g, 0.975)))
})) %>%
  dplyr::mutate(clear_of_zero = (ci_lo > 0 & ci_hi > 0) | (ci_lo < 0 & ci_hi < 0))

message("\n   A2: the gap change per contrast, with BOTH components in RANK STEPS")
mitorib_by_contrast %>%
  dplyr::select(contrast, d_oxphos_steps, d_mitorib_steps, gap_change,
                gap_change_steps, ci_lo, ci_hi, clear_of_zero) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

assert_supplied <- tibble::tibble(
  contrast = names(SUPPLIED), supplied = unname(SUPPLIED),
  computed = mitorib_by_contrast$gap_change[match(names(SUPPLIED),
                                                  mitorib_by_contrast$contrast)])
assert_supplied$delta <- assert_supplied$computed - assert_supplied$supplied
SUPPLIED_OK <- max(abs(assert_supplied$delta)) < SUPPLIED_TOL
print(as.data.frame(assert_supplied), digits = 4, row.names = FALSE)
message("   the supplied estimates reproduce: ", SUPPLIED_OK,
        " (max |delta| ", sprintf("%.3f", max(abs(assert_supplied$delta))),
        ", tolerance ", SUPPLIED_TOL, ")")
if (!SUPPLIED_OK) message("   *** DISCREPANCY - reported, not accepted silently")

message("\n   THE TIER COINCIDENCE, named so the exactness is not read as meaning:")
message("   the two mouse windows both land on ",
        mitorib_by_contrast$gap_change_steps[1], " steps from components (",
        mitorib_by_contrast$d_oxphos_steps[1], ", ",
        mitorib_by_contrast$d_mitorib_steps[1], ") and (",
        mitorib_by_contrast$d_oxphos_steps[2], ", ",
        mitorib_by_contrast$d_mitorib_steps[2], ");")
message("   the two within-timepoint MYC contrasts both land on ",
        mitorib_by_contrast$gap_change_steps[3], " from (",
        mitorib_by_contrast$d_oxphos_steps[3], ", ",
        mitorib_by_contrast$d_mitorib_steps[3], ") and (",
        mitorib_by_contrast$d_oxphos_steps[4], ", ",
        mitorib_by_contrast$d_mitorib_steps[4], ").")
message("   DIFFERENT COMPONENTS, IDENTICAL DIFFERENCE. A coincidence, not structure.")

# --- R2, within the mouse, on one paired resample of the four mouse groups ----
MOUSE_ARMS <- ARM_ORDER[1:4]
MC <- names(CONTRASTS)[1:4]
tier_boot <- vapply(seq_len(NBOOT), function(i) {
  p <- vapply(MOUSE_ARMS, function(a) {
    z <- ARM[[a]]; .gp(z, sample(z$j, replace = TRUE))[c(REF, MTRIB)] }, numeric(2))
  gg <- p[1, ] - p[2, ]; names(gg) <- MOUSE_ARMS
  ch <- vapply(MC, function(ct) {
    cn <- CONTRASTS[[ct]]; gg[[cn[["b"]]]] - gg[[cn[["a"]]]] }, 0)
  mean(abs(ch[1:2])) - mean(abs(ch[3:4]))
}, 0)
tier <- tibble::tibble(
  comparison = "mouse developmental tier minus mouse MYC tier, absolute magnitudes",
  dev_mean = mean(abs(mitorib_by_contrast$gap_change[1:2])),
  myc_mean = mean(abs(mitorib_by_contrast$gap_change[3:4])),
  difference = mean(abs(mitorib_by_contrast$gap_change[1:2])) -
               mean(abs(mitorib_by_contrast$gap_change[3:4])),
  ci_lo = unname(stats::quantile(tier_boot, 0.025)),
  ci_hi = unname(stats::quantile(tier_boot, 0.975)))
tier$separates <- tier$ci_lo > 0
message("\n   R2: the tiers, WITHIN THE MOUSE (one paired resample of all four groups)")
print(as.data.frame(tier), digits = 3, row.names = FALSE)
R2 <- if (tier$separates) {
  paste0("the developmental tier SEPARATES from the MYC tier within the mouse (",
         sprintf("%.1f", tier$difference), " [", sprintf("%.1f", tier$ci_lo),
         ", ", sprintf("%.1f", tier$ci_hi), "] percentile points)")
} else {
  paste0("the separation is NOT ESTABLISHED at this n (", sprintf("%.1f", tier$difference),
         " [", sprintf("%.1f", tier$ci_lo), ", ", sprintf("%.1f", tier$ci_hi),
         "]). The point estimates are reported as two tiers ONLY as point ",
         "estimates, and the exactness to the rank step is a coincidence")
}
message("   R2: ", R2)

mouse_rows <- mitorib_by_contrast[mitorib_by_contrast$species == "mouse", ]
R3 <- paste0(sum(mouse_rows$clear_of_zero), " of 4 mouse contrasts have intervals ",
             "clear of zero at n = 6. ",
             if (all(mouse_rows$clear_of_zero)) {
               "The mouse side stands on its own here"
             } else {
               paste0("Covering zero does not refute the point estimate; it means ",
                      "THE MOUSE SIDE CANNOT STAND ALONE and must be read beside ",
                      "the level-ruler evidence in E20/E21")
             })
message("   R3: ", R3)

# --- A3 the estimator ladder --------------------------------------------------
message("\n   A3: the estimator ladder - what the central-dogma gap failed")
EST_PANEL <- e25$settings$estimator_panel
.est_vec <- function(coh, w) {
  if (identical(coh, "TCGA")) {
    if (identical(w, "log2MYC")) return(as.numeric(nw$tcga_log2MYC[ID_T]))
    if (w %in% rownames(nw$tcga_M_b_variants)) return(as.numeric(nw$tcga_M_b_variants[w, ID_T]))
    return(as.numeric(nw$tcga_gsva_new[w, ID_T]))
  }
  if (identical(w, "log2MYC")) return(as.numeric(sc$log2MYC[ID_S]))
  if (w %in% rownames(sc$M_b_variants)) return(as.numeric(sc$M_b_variants[w, ID_S]))
  as.numeric(sc$gsva_new[w, ID_S])
}
.gap_diff <- function(U, g, comparator) {
  pl <- .pct(rowMeans(U[, which(g == "MYC-low"), drop = FALSE]))
  ph <- .pct(rowMeans(U[, which(g == "MYC-high"), drop = FALSE]))
  names(pl) <- names(ph) <- shared
  unname((ph[REF] - ph[comparator]) - (pl[REF] - pl[comparator]))
}
amp <- cvd$MYC_amp[match(ID_T, cvd$patient)]
amp_g <- ifelse(is.na(amp), NA_character_, ifelse(amp, "MYC-high", "MYC-low"))

estimator_ladder <- dplyr::bind_rows(lapply(names(EST_PANEL), function(en) {
  gT <- .tertiles(.est_vec("TCGA", EST_PANEL[[en]]), ID_T)
  gS <- .tertiles(.est_vec("SCAN-B", EST_PANEL[[en]]), ID_S)
  tibble::tibble(instrument = en, set = unname(EST_PANEL[[en]]),
                 exempt = en == "myc_regulon",
                 mitorib_TCGA = .gap_diff(UT, gT, MTRIB),
                 `mitorib_SCAN-B` = .gap_diff(US, gS, MTRIB),
                 cd_TCGA = .gap_diff(UT, gT, CD),
                 `cd_SCAN-B` = .gap_diff(US, gS, CD))
}))
estimator_ladder <- dplyr::bind_rows(
  estimator_ladder,
  tibble::tibble(instrument = "TCGA copy number (MYC_amp)",
                 set = "DNA-level, NOT an expression programme", exempt = FALSE,
                 mitorib_TCGA = .gap_diff(UT, amp_g, MTRIB),
                 `mitorib_SCAN-B` = NA_real_,
                 cd_TCGA = .gap_diff(UT, amp_g, CD), `cd_SCAN-B` = NA_real_))
print(as.data.frame(estimator_ladder[, -2]), digits = 3, row.names = FALSE)
message("   the cd_ columns use the SAME numerator (", REF,
        ") so the two comparators are like for like.")
message("   E28's AS-RUN central-dogma ladder used the UMBRELLA and gave ",
        "FELSHER -20.4 / -14.1, HALLMARK +23.2 / +21.8, copy number +25.4.")

lad <- estimator_ladder[!estimator_ladder$exempt &
                          estimator_ladder$instrument != "TCGA copy number (MYC_amp)", ]
PRIMARY_SIGN <- sign(lad$mitorib_TCGA[lad$instrument == "myc_felsher"])
agree_both <- lad$instrument[sign(lad$mitorib_TCGA) == PRIMARY_SIGN &
                               sign(lad$`mitorib_SCAN-B`) == PRIMARY_SIGN]
AMP_ROW <- estimator_ladder[estimator_ladder$instrument == "TCGA copy number (MYC_amp)", ]
AMP_AGREES <- sign(AMP_ROW$mitorib_TCGA) == PRIMARY_SIGN
REG <- estimator_ladder[estimator_ladder$exempt, ]
LADDER_HOLDS <- length(agree_both) == nrow(lad) && AMP_AGREES
R1 <- if (LADDER_HOLDS) {
  paste0("HOLDS. All ", nrow(lad), " non-exempt estimators keep the primary's ",
         "sign in BOTH cohorts and the COPY-NUMBER split agrees (",
         sprintf("%+.1f", AMP_ROW$mitorib_TCGA),
         "). The rank ruler CORROBORATES E20/E21, and it passes where the ",
         "parent comparator failed - which is the informative case, given the ",
         "two statistics are nested at rho ", sprintf("%+.3f", NEST_RHO))
} else {
  paste0("FAILS. ", length(agree_both), " of ", nrow(lad),
         " non-exempt estimators keep the primary's sign in both cohorts and ",
         "the copy-number split ", if (AMP_AGREES) "agrees" else "DISAGREES",
         " (", sprintf("%+.1f", AMP_ROW$mitorib_TCGA),
         "). It fails as E28's did, and given the two statistics are nested at ",
         "rho ", sprintf("%+.3f", NEST_RHO),
         " that is close to expected: THE FAILURE IS A PROPERTY OF THE RANK-GAP ",
         "STATISTIC, not of the comparator")
}
message("\n   R1: ", R1)
message("   myc_regulon is EXEMPT and reported separately: ",
        sprintf("%+.1f / %+.1f", REG$mitorib_TCGA, REG$`mitorib_SCAN-B`),
        ". E20 already recorded it as the one estimator not showing the")
message("   level-ruler effect and E21 explained that as a size artefact.")

# --- A4 composition, TCGA only ------------------------------------------------
message("\n   A4: composition, TCGA only (SCAN-B has no purity estimate)")
ADIPO <- e26$settings$adipo_markers
tl <- readRDS(PATH_TCGA_LINEAR)
stopifnot(identical(tl$scale, "linear_deseq2_normalised"))
adipose <- colMeans(t(scale(t(log2(tl$mat[ADIPO, ID_T, drop = FALSE] + 1)))))
rm(tl); invisible(gc(verbose = FALSE))
kk <- match(ID_T, cvd$patient)
COVT <- data.frame(id = ID_T, purity = cvd$purity[kk],
                   leukocyte = cvd$leukocyte_fraction[kk], adipose = adipose,
                   grp = grp_T, stringsAsFactors = FALSE)
keep <- !is.na(COVT$grp) & !is.na(COVT$purity) & !is.na(COVT$leukocyte)
FT <- COVT[keep, , drop = FALSE]; UF <- UT[, FT$id, drop = FALSE]
.residualise <- function(M, D) {
  X <- if (length(D)) cbind(1, as.matrix(D)) else matrix(1, ncol(M), 1L)
  R <- t(qr.resid(qr(X), t(M))); dimnames(R) <- dimnames(M); R
}
MODELS <- e26$settings$models
composition <- dplyr::bind_rows(lapply(names(MODELS), function(m) {
  R <- .residualise(UF, if (length(MODELS[[m]])) FT[, MODELS[[m]], drop = FALSE] else NULL)
  dplyr::bind_rows(lapply(c("MYC-low", "MYC-high"), function(g) {
    p <- .pct(rowMeans(R[, FT$grp == g, drop = FALSE])); names(p) <- shared
    tibble::tibble(model = m, group = g, oxphos = unname(p[REF]),
                   mitorib = unname(p[MTRIB]), gap = unname(p[REF] - p[MTRIB]))
  }))
}))
print(as.data.frame(composition), digits = 3, row.names = FALSE)
.maxmove <- function(col) {
  m0 <- composition[[col]][composition$model == "M0 unadjusted"]
  max(abs(vapply(setdiff(names(MODELS), "M0 unadjusted"), function(m)
    max(abs(composition[[col]][composition$model == m] - m0)), 0)))
}
comp_summary <- tibble::tibble(
  quantity = c(paste(REF, "position"), paste(MTRIB, "position"), "THE GAP"),
  max_move_from_M0 = c(.maxmove("oxphos"), .maxmove("mitorib"), .maxmove("gap")))
print(as.data.frame(comp_summary), digits = 3, row.names = FALSE)
PROTECTS <- comp_summary$max_move_from_M0[3] < min(comp_summary$max_move_from_M0[1:2])
R4 <- if (PROTECTS) {
  paste0("differencing PROTECTS here where it did not for central dogma: the gap ",
         "moves ", sprintf("%.1f", comp_summary$max_move_from_M0[3]),
         " percentile points against ",
         sprintf("%.1f and %.1f", comp_summary$max_move_from_M0[1],
                 comp_summary$max_move_from_M0[2]), " for its components")
} else {
  paste0("differencing does NOT protect: the gap moves ",
         sprintf("%.1f", comp_summary$max_move_from_M0[3]), " against ",
         sprintf("%.1f and %.1f", comp_summary$max_move_from_M0[1],
                 comp_summary$max_move_from_M0[2]),
         " for its components, as E28's did (2.1 against 1.4 and 0.7). Every gap ",
         "number inherits that")
}
message("   R4: ", R4)

# --- A5 proliferation strata, both cohorts, NO adjustment ---------------------
message("\n   A5: proliferation STRATIFICATION only. E28 showed the adjustment")
message("   uninformative here - its positive control GREW 2.2x under it.")
PROLIF <- "PROLIF_DISJOINT"
pro_T <- as.numeric(mito$gsva_cov[PROLIF, ID_T])
pro_S <- as.numeric(sc$gsva_cov[PROLIF, ID_S])
.strata <- function(U, pro, est, ids) {
  q <- stats::quantile(pro, c(TERTILE, 1 - TERTILE), names = FALSE)
  s <- ifelse(pro <= q[1], "low prolif", ifelse(pro >= q[2], "high prolif", "mid prolif"))
  dplyr::bind_rows(lapply(c("low prolif", "mid prolif", "high prolif"), function(lv) {
    j <- which(s == lv); g <- .tertiles(est[j], ids[j])
    tibble::tibble(stratum = lv, n_stratum = length(j),
                   n_low = sum(g == "MYC-low", na.rm = TRUE),
                   n_high = sum(g == "MYC-high", na.rm = TRUE),
                   gap_diff = .gap_diff(U[, j, drop = FALSE], g, MTRIB))
  }))
}
strata <- dplyr::bind_rows(
  dplyr::mutate(.strata(UT, pro_T, est_T, ID_T), cohort = "TCGA"),
  dplyr::mutate(.strata(US, pro_S, est_S, ID_S), cohort = "SCAN-B"))
strata <- strata[, c("cohort", "stratum", "n_stratum", "n_low", "n_high", "gap_diff")]
print(as.data.frame(strata), digits = 3, row.names = FALSE)
STRATA_HOLD <- sum(sign(strata$gap_diff) == PRIMARY_SIGN)
message("   ", STRATA_HOLD, " of ", nrow(strata),
        " strata keep the primary's sign")

# =============================================================================
# 4. PART B - the NADH candidate panel. DESCRIPTIVE ONLY.
# =============================================================================

message("\n4. PART B: the NADH candidate panel")
message("   THE ANNOTATION IS THE AUTHOR'S, PROVISIONAL, AMBIGUOUS AND CONTESTED.")
message("   It is NOT a classification and NOT a declared rule. No set-level")
message("   statistic is computed and no programme is named.")

NADH_FLAT <- unlist(NADH, use.names = FALSE)
NADH_GRP  <- rep(names(NADH), lengths(NADH))
names(NADH_GRP) <- NADH_FLAT
ABSENT <- setdiff(NADH_FLAT, shared)
message("   members: ", length(NADH_FLAT), "; absent from the shared panel: ",
        if (length(ABSENT)) paste(ABSENT, collapse = ", ") else "NONE")
NADH_IN <- intersect(NADH_FLAT, shared)

refg <- HP[[REF]]
nadh_meta <- tibble::tibble(
  pathway = NADH_IN,
  annotation = unname(NADH_GRP[NADH_IN]),
  caveat = unname(ifelse(is.na(NADH_NOTE[NADH_IN]), "", NADH_NOTE[NADH_IN])),
  n_genes = vapply(NADH_IN, function(p) length(HP[[p]]), 0L),
  shared_with_numerator = vapply(NADH_IN, function(p)
    length(intersect(HP[[p]], refg)), 0L))
nadh_meta$pct_shared <- 100 * nadh_meta$shared_with_numerator / nadh_meta$n_genes
nadh_meta$tautology <- dplyr::case_when(
  nadh_meta$pct_shared >= 99 ~ "FULLY CONTAINED - the gap is ARITHMETIC",
  nadh_meta$pct_shared >= 10 ~ "partly shared - the gap is partly arithmetic",
  TRUE ~ "")
print(as.data.frame(nadh_meta[, c("pathway", "annotation", "n_genes",
                                  "shared_with_numerator", "pct_shared",
                                  "tautology")]), digits = 3, row.names = FALSE)

BGN <- lapply(ARM_ORDER, function(a) .boot_group(a, c(REF, NADH_IN)))
names(BGN) <- ARM_ORDER
nadh_panel <- dplyr::bind_rows(lapply(ARM_ORDER, function(a) {
  b <- BGN[[a]]
  dplyr::bind_rows(lapply(NADH_IN, function(p) {
    g <- b[REF, ] - b[p, ]
    tibble::tibble(arm = a, species = ARM[[a]]$sp, pathway = p,
                   annotation = unname(NADH_GRP[[p]]),
                   oxphos_pct = pos[REF, a], pathway_pct = pos[p, a],
                   gap = pos[REF, a] - pos[p, a],
                   gap_lo = unname(stats::quantile(g, 0.025)),
                   gap_hi = unname(stats::quantile(g, 0.975)))
  }))
})) %>%
  dplyr::mutate(arm = factor(arm, levels = ARM_ORDER),
                pct_shared = nadh_meta$pct_shared[match(pathway, nadh_meta$pathway)],
                tautology = nadh_meta$tautology[match(pathway, nadh_meta$pathway)])

BCN <- lapply(names(CONTRASTS), function(ct) .boot_contrast(ct, c(REF, NADH_IN)))
names(BCN) <- names(CONTRASTS)
nadh_contrasts <- dplyr::bind_rows(lapply(names(CONTRASTS), function(ct) {
  cn <- CONTRASTS[[ct]]; b <- BCN[[ct]]
  dplyr::bind_rows(lapply(NADH_IN, function(p) {
    g <- b[REF, ] - b[p, ]
    d_ox <- pos[REF, cn[["b"]]] - pos[REF, cn[["a"]]]
    d_p  <- pos[p, cn[["b"]]] - pos[p, cn[["a"]]]
    tibble::tibble(contrast = ct, kind = unname(CONTRAST_KIND[[ct]]),
                   pathway = p, annotation = unname(NADH_GRP[[p]]),
                   d_oxphos = d_ox, d_pathway = d_p, gap_change = d_ox - d_p,
                   ci_lo = unname(stats::quantile(g, 0.025)),
                   ci_hi = unname(stats::quantile(g, 0.975)))
  }))
}))

message("\n   B1/B2: every row carries BOTH components. Human MYC contrasts:")
nadh_contrasts %>%
  dplyr::filter(grepl("MYC-high - low", contrast)) %>%
  dplyr::select(contrast, annotation, pathway, d_oxphos, d_pathway, gap_change) %>%
  as.data.frame() %>%
  print(digits = 3, row.names = FALSE)

aligned_carryover <- e29$aligned_by_species %>%
  dplyr::filter(pathway %in% NADH_IN) %>%
  dplyr::mutate(annotation = unname(NADH_GRP[pathway])) %>%
  dplyr::select(annotation, pathway, mouse_mean, human_mean, quadrant)
message("\n   B3: E29's aligned-change means carried forward")
print(as.data.frame(aligned_carryover), digits = 3, row.names = FALSE)
message("   THIS PANEL SITS BESIDE E29'S MAP, IT DOES NOT REPLACE IT. No set-level")
message("   statistic, no direction count, no programme named (R5).")

# =============================================================================
# 5. PART C - figures
# =============================================================================

message("\n5. PART C: figures")
if (!dir.exists(DIR_E30_FIG)) dir.create(DIR_E30_FIG, recursive = TRUE)

CAP <- paste(
  "EVERY GAP IS SHOWN WITH BOTH COMPONENT POSITIONS - a gap discards the thing it is built from (E28).",
  "Percentiles only; no mitoPPS value crosses the species boundary.",
  "Mouse contrasts come from n = 6; the mouse timeline is batch-confounded (batch = timepoint).",
  "Exploratory, post-hoc, descriptive; not a hypothesis test.", sep = "\n")

f1 <- dplyr::bind_rows(
  mitorib_by_group %>% dplyr::transmute(arm, species, series = REF,
                                        y = oxphos_pct, lo = ox_lo, hi = ox_hi),
  mitorib_by_group %>% dplyr::transmute(arm, species, series = MTRIB,
                                        y = mitorib_pct, lo = mt_lo, hi = mt_hi),
  mitorib_by_group %>% dplyr::transmute(arm, species, series = "the GAP",
                                        y = gap, lo = gap_lo, hi = gap_hi)) %>%
  dplyr::mutate(series = factor(series, levels = c(REF, MTRIB, "the GAP")))
p1 <- ggplot2::ggplot(f1, ggplot2::aes(arm, y, colour = series, group = series)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey45") +
  ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.4, colour = "grey70",
                      linetype = 3) +
  ggplot2::geom_line(linewidth = 0.4, alpha = 0.6) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi),
                          position = ggplot2::position_dodge(width = 0.45),
                          linewidth = 0.6) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45), size = 2.2) +
  ggplot2::scale_colour_manual(values = stats::setNames(
    c("#0072B2", "#009E73", "#D55E00"), c(REF, MTRIB, "the GAP"))) +
  ggplot2::labs(x = NULL, y = "percentile, and the gap in percentile points",
                colour = NULL,
                title = "PART A: the mitoribosome gap, with both components",
                subtitle = paste0("Three series, not one: the two positions are ",
                                  "visible because a gap alone hides them (E28)."),
                caption = CAP) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                 legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E30_FIG, "E30_C1_mitoribosome_gap.pdf"), p1,
                width = 8, height = 5.5)

f2 <- dplyr::bind_rows(
  nadh_panel %>% dplyr::transmute(arm, pathway, annotation, series = REF, y = oxphos_pct),
  nadh_panel %>% dplyr::transmute(arm, pathway, annotation, series = "the pathway",
                                  y = pathway_pct),
  nadh_panel %>% dplyr::transmute(arm, pathway, annotation, series = "the GAP", y = gap)) %>%
  dplyr::mutate(series = factor(series, levels = c(REF, "the pathway", "the GAP")),
                facet = paste0(pathway, "\n[", annotation, "]",
                               ifelse(pathway %in% nadh_meta$pathway[nadh_meta$tautology != ""],
                                      "  *shares genes with the numerator*", "")))
p2 <- ggplot2::ggplot(f2, ggplot2::aes(arm, y, colour = series, group = series)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  ggplot2::geom_vline(xintercept = 4.5, linewidth = 0.35, colour = "grey75",
                      linetype = 3) +
  ggplot2::geom_line(linewidth = 0.35) +
  ggplot2::geom_point(size = 0.9) +
  ggplot2::scale_colour_manual(values = stats::setNames(
    c("#0072B2", "#009E73", "#D55E00"), c(REF, "the pathway", "the GAP"))) +
  ggplot2::facet_wrap(~ facet, ncol = 4L) +
  ggplot2::labs(x = NULL, y = "percentile, and the gap in percentile points",
                colour = NULL,
                title = "PART B: the NADH candidate panel - DESCRIPTIVE ONLY",
                subtitle = paste0("The grouping is the AUTHOR'S PROVISIONAL ",
                                  "ANNOTATION, not a result and not a classification. ",
                                  "No set-level statistic is computed."),
                caption = paste(CAP,
                  "Panels marked *shares genes with the numerator* have a gap that is partly or wholly ARITHMETIC; CI subunits is 100 pct contained.",
                  sep = "\n")) +
  ggplot2::theme_bw(base_size = 7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5),
                 strip.text = ggplot2::element_text(size = 5.5),
                 legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 5.5, hjust = 0))
ggplot2::ggsave(file.path(DIR_E30_FIG, "E30_C2_nadh_panel.pdf"), p2,
                width = 11, height = 12)

f3 <- estimator_ladder %>%
  tidyr::pivot_longer(c("mitorib_TCGA", "mitorib_SCAN-B", "cd_TCGA", "cd_SCAN-B"),
                      names_to = "key", values_to = "gap_diff") %>%
  tidyr::separate(key, into = c("comparator", "cohort"), sep = "_") %>%
  dplyr::filter(!is.na(gap_diff)) %>%
  dplyr::mutate(
    comparator = ifelse(comparator == "mitorib", MTRIB, CD),
    instrument = factor(instrument, levels = rev(estimator_ladder$instrument)),
    kind = ifelse(grepl("copy number", instrument), "DNA-level",
                  ifelse(exempt, "exempt (E20/E21)", "expression")))
p3 <- ggplot2::ggplot(f3, ggplot2::aes(gap_diff, instrument, colour = kind)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey45") +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_colour_manual(values = c(expression = "#0072B2",
                                          `DNA-level` = "#D55E00",
                                          `exempt (E20/E21)` = "grey60")) +
  ggplot2::facet_grid(comparator ~ cohort) +
  ggplot2::labs(x = "MYC-high minus MYC-low change in the gap (percentile points)",
                y = NULL, colour = NULL,
                title = "PART A3: the estimator ladder, both comparators, same numerator",
                subtitle = paste0("The copy-number split is the only instrument not ",
                                  "defined on the expression matrix. myc_regulon is ",
                                  "EXEMPT: E20 recorded its dissent and E21 explained it as size."),
                caption = CAP) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(legend.position = "top",
                 plot.caption = ggplot2::element_text(size = 6, hjust = 0))
ggplot2::ggsave(file.path(DIR_E30_FIG, "E30_C3_estimator_ladder.pdf"), p3,
                width = 9, height = 6)
message("   wrote three figures under ", DIR_E30_FIG)

# =============================================================================
# 6. PART D - verdict
# =============================================================================

message("\n6. PART D: verdict")

verdict <- tibble::tibble(
  rule = c("STANDING", "CONTROL", "the nesting", "the supplied estimates",
           "the tier coincidence", "R1 the ladder", "R2 the tiers",
           "R3 the mouse contrasts", "R4 composition", "A5 strata",
           "R5 the NADH panel"),
  outcome = c(
    paste("PART A is a declared check of E20/E21 on a second ruler.",
          "PART B is a descriptive panel on the author's PROVISIONAL annotation",
          "and is not a classification, a hypothesis or a readout."),
    paste0("PASSED - ", nrow(e25o), " positions reproduce E25 and ", nrow(e27p),
           " reproduce E27, max |delta| ", format(max(c25, c27), digits = 2)),
    paste0("PART A's comparator is 100 pct contained in E28's failed one (all ",
           length(HP[[MTRIB]]), " genes) and the two gaps correlate at rho ",
           sprintf("%+.3f", NEST_RHO),
           ". The prior is independent; the statistic is not"),
    paste0(if (SUPPLIED_OK) "reproduce" else "DO NOT reproduce",
           " within ", SUPPLIED_TOL, " (max |delta| ",
           sprintf("%.3f", max(abs(assert_supplied$delta))), ")"),
    paste0("the two mouse windows both land on ",
           mitorib_by_contrast$gap_change_steps[1],
           " rank steps and the two within-timepoint MYC contrasts both on ",
           mitorib_by_contrast$gap_change_steps[3],
           ", from DIFFERENT components. A COINCIDENCE, not structure"),
    R1, R2, R3, R4,
    paste0(STRATA_HOLD, " of ", nrow(strata),
           " proliferation strata keep the primary's sign; the adjustment was",
           " NOT run because E28 showed it uninformative here"),
    paste("DESCRIPTIVE ONLY. No set-level statistic, no direction count, no",
          "programme named. One member is FULLY CONTAINED in the numerator and",
          "two are partly, so those gaps are arithmetic and are marked")))
print(as.data.frame(verdict), row.names = FALSE, right = FALSE)
message("\n   AND WHATEVER IT SAYS, IT IS DESCRIPTIVE. NOT A TEST.")

NOTES <- c(
  paste("TWO PARTS, DIFFERENT STANDING. PART A is a declared check of E20/E21 on",
        "a second ruler with the comparator named in advance. PART B is a",
        "descriptive panel on the author's PROVISIONAL, ambiguous and contested",
        "annotation; it is not a classification, nothing in it is defended as",
        "set membership, and nothing may be promoted to a readout."),
  paste("EVERY GAP IS REPORTED WITH BOTH COMPONENT POSITIONS, in every table and",
        "every figure. E28's failure mode was that a gap discards the thing it",
        "is built from - two groups gave an identical gap from positions of 14",
        "and 62."),
  paste("PART A IS A REFINEMENT OF THE STATISTIC E28 FAILED WITH, NOT AN",
        "INDEPENDENT ONE. The mitoribosome is 100 per cent contained in",
        "Mitochondrial central dogma and the two gaps correlate at rho",
        sprintf("%+.3f", NEST_RHO),
        "across the eight groups. The PRIOR is independent - E20/E21 named the",
        "mitoribosome in advance on a level ruler - but the STATISTIC is nested.",
        "The informative case is PASSING where the parent FAILED."),
  paste("The numerator here is", REF, "where E28's was the OXPHOS UMBRELLA,",
        "which E27 showed is exactly subunits plus assembly factors with",
        "assembly factors barely moving, and which E28 measured as diluting by",
        "8.5 percentile points. Both terms changed and one change is an",
        "improvement."),
  paste("R1 EXEMPTS myc_regulon, with the reason stated rather than the",
        "exemption hidden.", E20_LEVEL),
  paste("THE TWO TIERS ARE EXACT TO THE RANK STEP AND THAT IS A COINCIDENCE.",
        "Both mouse windows land on the same integer from different components,",
        "and so do both within-timepoint MYC contrasts. The components are",
        "reported in rank steps so the exactness cannot be read as meaning.",
        "R2's tier comparison is run WITHIN THE MOUSE, on one paired resample of",
        "all four mouse groups, which avoids a cross-species magnitude",
        "comparison in percentile space."),
  paste("ONE PART B ENTRY IS FULLY TAUTOLOGICAL AND TWO ARE PARTLY SO.",
        "CI subunits is 100 per cent contained in the numerator so its gap is",
        "ARITHMETIC; Complex I shares 63 per cent, TCA cycle 20 and Electron",
        "carriers 14. Everything else shares zero, INCLUDING CI assembly",
        "factors. Gene sharing is on every row and marked on every panel."),
  paste("Composition is TCGA ONLY: SCAN-B has no purity estimate and nothing is",
        "imputed. The proliferation ADJUSTMENT was NOT run - E28 showed it",
        "uninformative here because its positive control grew 2.2 times under",
        "it - and stratification is reported instead, in both cohorts."),
  paste("E28's GAP READOUT FAILED AND IS NOT REHABILITATED BY THIS SCRIPT.",
        "If PART A fails its ladder it joins E28 as a negative."),
  paste("R5: the NADH panel is descriptive. No set-level statistic, no counting",
        "members by direction as though the count were a result, no programme",
        "named. E29 already found the off-diagonal structure diffuse and the",
        "sign counts incoherent; this is finer-grained description of the same",
        "object, not a second attempt at the same claim."),
  paste("The mouse timeline is BATCH-CONFOUNDED (batch = timepoint); reported in",
        "full and nothing withdrawn on that account. n = 6 per mouse group."),
  "Nothing re-scored. No mitoPPS value crosses the species boundary. No ortholog call. No FDR. N3 throughout.")

saveRDS(list(
  rules = RULES, e20_level_ruler = E20_LEVEL,
  settings = list(reference = REF, comparator = MTRIB, e28_comparator = CD,
                  arm_order = ARM_ORDER, contrasts = CONTRASTS,
                  nadh = NADH, nboot = NBOOT, seed = 1, tertile = TERTILE,
                  mouse_artefact = list(md5 = MOUSE_MD5, recomputed = FALSE)),
  control = list(vs_e25 = c25, vs_e27 = c27,
                 supplied = assert_supplied, supplied_ok = SUPPLIED_OK,
                 nesting = nest, nesting_rho = NEST_RHO),
  mitorib_by_group = mitorib_by_group, mitorib_by_contrast = mitorib_by_contrast,
  tier = tier, estimator_ladder = estimator_ladder,
  ladder_holds = LADDER_HOLDS, amp_agrees = AMP_AGREES,
  composition = composition, composition_summary = comp_summary,
  protects = PROTECTS, strata = strata,
  nadh_meta = nadh_meta, nadh_panel = nadh_panel, nadh_contrasts = nadh_contrasts,
  aligned_carryover = aligned_carryover,
  readings = list(R1 = R1, R2 = R2, R3 = R3, R4 = R4),
  verdict = verdict, analysis_date = Sys.Date(), notes = NOTES), PATH_E30)

readr::write_csv(mitorib_by_group %>% dplyr::mutate(arm = as.character(arm)), PATH_E30_A)
readr::write_csv(nadh_panel %>% dplyr::mutate(arm = as.character(arm)), PATH_E30_B)

message("\nE30: done.")
message("    ", PATH_E30); message("    ", PATH_E30_A); message("    ", PATH_E30_B)
message("    three figures under ", DIR_E30_FIG)

# =============================================================================
# SANDBOX -- run line-by-line in Positron; skipped by source()
# =============================================================================
if (FALSE) {

  x <- readRDS(PATH_E30)

  ## controls, and the nesting this whole script sits inside
  x$control$vs_e25; x$control$vs_e27
  x$control$nesting |> as.data.frame(); x$control$nesting_rho
  x$control$supplied |> as.data.frame()

  ## the verdict and the four readings
  x$verdict |> as.data.frame()
  x$readings

  ## PART A - BOTH COMPONENTS ARE ON EVERY ROW. Never read the gap alone.
  x$mitorib_by_group |> as.data.frame()
  x$mitorib_by_contrast |> as.data.frame()
  x$tier |> as.data.frame()

  ## A3 - the ladder that killed the central-dogma gap. cd_ columns use the
  ## SAME numerator, so the two comparators are like for like.
  x$estimator_ladder |> as.data.frame()
  x$ladder_holds; x$amp_agrees
  x$e20_level_ruler

  ## A4/A5
  x$composition |> as.data.frame()
  x$composition_summary |> as.data.frame(); x$protects
  x$strata |> as.data.frame()

  ## PART B - DESCRIPTIVE. Read nadh_meta FIRST: some gaps are arithmetic.
  x$nadh_meta |> as.data.frame()
  subset(x$nadh_contrasts, grepl("MYC-high - low", contrast)) |> as.data.frame()
  x$aligned_carryover |> as.data.frame()

  ## every caveat the note must repeat
  cat(paste0("- ", x$notes, collapse = "\n"), "\n")
}
