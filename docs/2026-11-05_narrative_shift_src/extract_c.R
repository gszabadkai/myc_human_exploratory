suppressMessages({library(dplyr); library(jsonlite)})
x49 <- readRDS("results/e49_signature_and_fig4g.rds"); x50 <- readRDS("results/e50_fig4g_in_fig1g_style.rds")
keep <- c("BCL2L1", "BBC3", "S6")
cells <- x50$values %>% filter(panel == "fig4g", readout %in% keep) %>% select(readout, cohort, term, value, printed, fill, ink)
ci <- x49$all_readouts_m3 %>% filter(estimator == "clean", instrument == "gsva", readout %in% keep) %>%
  select(readout, cohort, b_OX_m3, ox_lo_m3, ox_hi_m3, b_MYC_m3, myc_lo_m3, myc_hi_m3, b_PR, pr_lo, pr_hi)
# Check: E50's drawn values equal E49's stored coefficients.
chk <- cells %>% left_join(ci, by = c("readout", "cohort")) %>%
  mutate(stored = dplyr::case_when(term == "OX" ~ b_OX_m3, term == "MYC" ~ b_MYC_m3, term == "PROLIF" ~ b_PR))
stopifnot(nrow(chk) == 18, max(abs(chk$value - chk$stored)) < 1e-12)
cat("panel c: 18 cells, E50 values equal E49 coefficients\n")
write_json(list(cells = cells, ci = ci, lim = 0.6, pal = list(neg = "#4A3525", zero = "#FAFAFA", pos = "#2A8A6D")),
           commandArgs(TRUE)[1], digits = NA, auto_unbox = TRUE, pretty = TRUE)
