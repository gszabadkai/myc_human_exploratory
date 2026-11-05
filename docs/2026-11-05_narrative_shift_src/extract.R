# Read-only extraction of verified, saved numbers for the narrative figure.
suppressMessages({library(dplyr); library(jsonlite)})
out <- commandArgs(TRUE)[1]
x49 <- readRDS("results/e49_signature_and_fig4g.rds")
x50 <- readRDS("results/e50_fig4g_in_fig1g_style.rds")
x47 <- readRDS("results/e47_mutual_adjustment.rds")
x40 <- readRDS("results/e40_ladder.rds")
x43 <- readRDS("results/e43_proliferation_coefficient.rds")
x44 <- readRDS("results/e44_contrast_test.rds")
raw <- x49$raw_correlations %>%
  filter(readout %in% c("S6", "BCL2L1", "BBC3"),
         exposure %in% c("OX_gsva", "MYC_M_a", "MYC_clean", "PROLIF"))
sig <- x49$signature_coefficients %>% filter(readout == "S6", instrument == "gsva") %>%
  select(cohort, estimator, void, b_MYC_m1, myc_lo_m1, myc_hi_m1, b_MYC_m3, myc_lo_m3, myc_hi_m3, b_PR, pr_lo, pr_hi)
g3 <- x47$gates$g3 %>% filter(instrument == "gsva") %>% select(cohort, estimator, rho_MYC_PROLIF, rho_OX_PROLIF, rho_OX_MYC)
cells <- x50$values %>% filter(panel == "fig4g")
cat("E50 fig4g readouts:", unique(cells$readout), "\n")
ar <- x49$all_readouts_m3; cat("all_readouts_m3 cols:", names(ar), "\n")
pcr <- list(ox = x40$meta %>% select(rung, fe_estimate, fe_ci_lo, fe_ci_hi, I2, signs),
            prolif = x43$meta %>% filter(term == "PROLIF") %>% select(re_estimate, re_ci_lo, re_ci_hi, I2, signs),
            delta = x44$meta %>% select(fe_estimate, fe_ci_lo, fe_ci_hi, I2, signs))
write_json(list(raw = raw, sig = sig, g3 = g3, cells = cells, pcr = pcr), out, digits = NA, auto_unbox = TRUE, pretty = TRUE)
cat("written\n")
