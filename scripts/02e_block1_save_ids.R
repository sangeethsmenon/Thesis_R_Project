library(readr); library(here); library(dplyr)

summ <- readr::read_csv(here("results","block1_sample_missingness_proteome.csv"), show_col_types = FALSE)
raw  <- readr::read_tsv(here("data","raw","NMR_proteome_n2923_SN13411_20250829.txt"), show_col_types = FALSE)

keep <- summ |> mutate(keep = pct_missing <= 50) |> select(sid, keep)
readr::write_csv(
  keep |> filter(keep) |> select(sid),
  here("results","block1_ids_after_sample_filter.csv")
)
