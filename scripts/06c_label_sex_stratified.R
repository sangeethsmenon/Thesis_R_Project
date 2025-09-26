# scripts/06c_label_sex_stratified.R
packages <- c("tidyverse","here","readr")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

info <- read_tsv(here("data","raw","metabolite_info.txt"),
                 col_types = "cc", show_col_types = FALSE) %>%
  setNames(c("metabolite","metabolite_label"))

resdir <- here("results")
edges_diff <- read_csv(file.path(resdir,"block5c_edges_top20k_by_pdiff.csv"), show_col_types = FALSE)
edges_F    <- read_csv(file.path(resdir,"block5c_edges_F_q_le_alpha_absr.csv"), show_col_types = FALSE)
edges_M    <- read_csv(file.path(resdir,"block5c_edges_M_q_le_alpha_absr.csv"), show_col_types = FALSE)

label <- function(df){
  df %>% left_join(info, by = "metabolite") %>%
    mutate(metabolite_label = if_else(is.na(metabolite_label), metabolite, metabolite_label))
}

readr::write_csv(label(edges_diff),
                 file.path(resdir,"block5c_edges_top20k_by_pdiff_labeled.csv"))
readr::write_csv(label(edges_F),
                 file.path(resdir,"block5c_edges_F_q_le_alpha_absr_labeled.csv"))
readr::write_csv(label(edges_M),
                 file.path(resdir,"block5c_edges_M_q_le_alpha_absr_labeled.csv"))

cat("Sex-stratified tables labeled.\n")
