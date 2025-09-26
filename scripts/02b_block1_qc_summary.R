packages <- c("tidyverse","readr")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))
stopifnot(exists("config"))

metab_feat <- read_csv(file.path(config$paths$results, "metabolome_feature_missingness.csv"), show_col_types = FALSE)
prot_feat  <- read_csv(file.path(config$paths$results, "proteome_feature_missingness.csv"),  show_col_types = FALSE)
prot_row   <- read_csv(file.path(config$paths$results, "proteome_sample_missingness.csv"),   show_col_types = FALSE)
metab_miss_vs <- suppressWarnings(try(read_csv(file.path(config$paths$results, "metabolome_missingness_vs_covariates.csv"), show_col_types = FALSE), silent = TRUE))
prot_miss_vs  <- suppressWarnings(try(read_csv(file.path(config$paths$results, "proteome_missingness_vs_covariates.csv"),  show_col_types = FALSE), silent = TRUE))
retention <- read_csv(file.path(config$paths$results, "block1_retention_counts.csv"), show_col_types = FALSE)

qnums <- function(x) tibble(min = min(x), q25 = quantile(x, .25), median = median(x), q75 = quantile(x, .75), q90 = quantile(x, .90), q95 = quantile(x, .95), max = max(x))
mf <- qnums(metab_feat$pct_missing)
pf <- qnums(prot_feat$pct_missing)
pr <- qnums(prot_row$pct_missing)

thr_m <- 100*config$thresholds$metab_feat
thr_p <- 100*config$thresholds$prot_feat
thr_s <- 100*config$thresholds$prot_row

mf_over <- mean(metab_feat$pct_missing > thr_m)
pf_over <- mean(prot_feat$pct_missing  > thr_p)
pr_over <- mean(prot_row$pct_missing   > thr_s)

sex_age_counts <- function(df){
  if (inherits(df, "try-error") || nrow(df)==0) return(tibble(n = 0, n_sex_q = 0, n_age_q = 0, frac_sex_q = 0, frac_age_q = 0))
  df <- df %>% mutate(q_sex = ifelse(is.finite(q_sex), q_sex, NA_real_), q_age = ifelse(is.finite(q_age), q_age, NA_real_))
  n <- nrow(df); ns <- sum(df$q_sex < 0.05, na.rm = TRUE); na <- sum(df$q_age < 0.05, na.rm = TRUE)
  tibble(n = n, n_sex_q = ns, n_age_q = na, frac_sex_q = ns/n, frac_age_q = na/n)
}
metab_sa <- sex_age_counts(metab_miss_vs)
prot_sa  <- sex_age_counts(prot_miss_vs)

txt <- c(
  "BLOCK 1 QC SUMMARY",
  paste0("Metabolome features: n=", nrow(metab_feat)),
  paste0("  % missing (pct): min=", round(mf$min,2), " med=", round(mf$median,2), " p95=", round(mf$q95,2), " max=", round(mf$max,2)),
  paste0("  Over ", thr_m, "% missing: ", scales::percent(mf_over, accuracy = 0.1)),
  paste0("Proteome features: n=", nrow(prot_feat)),
  paste0("  % missing (pct): min=", round(pf$min,2), " med=", round(pf$median,2), " p95=", round(pf$q95,2), " max=", round(pf$max,2)),
  paste0("  Over ", thr_p, "% missing: ", scales::percent(pf_over, accuracy = 0.1)),
  paste0("Proteome samples: n=", nrow(prot_row)),
  paste0("  % missing (pct): min=", round(pr$min,2), " med=", round(pr$median,2), " p95=", round(pr$q95,2), " max=", round(pr$max,2)),
  paste0("  Over ", thr_s, "% missing: ", scales::percent(pr_over, accuracy = 0.1)),
  paste0("Missingness ~ covariates (metabolome): features=", metab_sa$n, " with q_sex<0.05=", metab_sa$n_sex_q, " (", scales::percent(metab_sa$frac_sex_q), "); q_age<0.05=", metab_sa$n_age_q, " (", scales::percent(metab_sa$frac_age_q), ")"),
  paste0("Missingness ~ covariates (proteome): features=", prot_sa$n, " with q_sex<0.05=", prot_sa$n_sex_q, " (", scales::percent(prot_sa$frac_sex_q), "); q_age<0.05=", prot_sa$n_age_q, " (", scales::percent(prot_sa$frac_age_q), ")"),
  "Retention counts:",
  paste(capture.output(print(retention, n=Inf)), collapse = "\n")
)
writeLines(txt, file.path(config$paths$results, "block1_qc_summary.txt"))
cat(paste(txt, collapse = "\n"), "\n")





