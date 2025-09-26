library(readr); library(dplyr); library(here)

res <- here("results")

read_csv(file.path(res,"proteome_sample_missingness.csv"), show_col_types = FALSE) %>%
  summarise(n = n(), median = median(pct_missing),
            p95 = quantile(pct_missing, .95), max = max(pct_missing),
            over30 = sum(pct_missing > 30), over40 = sum(pct_missing > 40)) %>% print()

read_csv(file.path(res,"proteome_feature_missingness.csv"), show_col_types = FALSE) %>%
  summarise(n = n(), median = median(pct_missing),
            p95 = quantile(pct_missing, .95), max = max(pct_missing),
            over30 = sum(pct_missing > 30)) %>% print()

read_csv(file.path(res,"metabolome_feature_missingness.csv"), show_col_types = FALSE) %>%
  summarise(n = n(), median = median(pct_missing),
            p95 = quantile(pct_missing, .95), max = max(pct_missing),
            over5 = sum(pct_missing > 5)) %>% print()

read_csv(file.path(res,"proteome_missingness_vs_covariates.csv"), show_col_types = FALSE) %>%
  summarise(n = n(), sex_q_lt05 = sum(q_sex < 0.05, na.rm = TRUE),
            age_q_lt05 = sum(q_age < 0.05, na.rm = TRUE)) %>% print()

read_csv(file.path(res,"metabolome_missingness_vs_covariates.csv"), show_col_types = FALSE) %>%
  summarise(n = n(), sex_q_lt05 = sum(q_sex < 0.05, na.rm = TRUE),
            age_q_lt05 = sum(q_age < 0.05, na.rm = TRUE)) %>% print()
