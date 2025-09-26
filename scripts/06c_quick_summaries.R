packages <- c("readr","dplyr","here")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

ed <- readr::read_csv(here::here("results","block5c_edges_top20k_by_pdiff.csv"), show_col_types = FALSE)

overall <- ed %>% summarize(
  n = n(),
  sign_flip = sum(sign(r_F) != sign(r_M)),
  pct_flip = round(100*sign_flip/n, 1),
  stronger_in_male = sum(abs(r_M) > abs(r_F)),
  stronger_in_female = sum(abs(r_F) > abs(r_M))
)

top_proteins <- ed %>%
  mutate(sign_flip = sign(r_F) != sign(r_M)) %>%
  count(protein, wt = (q_diff <= 0.05 & abs(r_diff) >= 0.10)) %>%
  arrange(desc(n)) %>% head(10)

print(overall); print(top_proteins)

readr::write_csv(overall,     here::here("results","block5c_quick_summary_overall.csv"))
readr::write_csv(top_proteins,here::here("results","block5c_quick_summary_top_proteins.csv"))

cat("Saved to results/: block5c_quick_summary_overall.csv and ..._top_proteins.csv\n")
