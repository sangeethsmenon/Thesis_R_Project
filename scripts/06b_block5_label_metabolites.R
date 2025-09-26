packages <- c("tidyverse","here","readr")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))
suppressPackageStartupMessages({ library(jsonlite) })

if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found. Run Block 0.")
}

raw_dir <- tryCatch(config$paths$data_raw, error = function(e) here::here("data","raw"))
info_path <- file.path(raw_dir, "metabolite_info.txt")
stopifnot(file.exists(info_path))

meta_info <- read_tsv(info_path, show_col_types = FALSE, col_types = "cc")
names(meta_info) <- c("metabolite","metabolite_name")

res_dir <- config$paths$results
edges_es_path <- file.path(res_dir, "block5_edges_pearson_q_le_alpha_absr.csv")
edges_top_path<- file.path(res_dir, "block5_edges_top10k_by_q.csv")
stopifnot(file.exists(edges_es_path), file.exists(edges_top_path))

edges_es  <- read_csv(edges_es_path,  show_col_types = FALSE)
edges_top <- read_csv(edges_top_path, show_col_types = FALSE)

label_edges <- function(df){
  df %>% left_join(meta_info, by = "metabolite") %>%
    mutate(metabolite_label = if_else(is.na(metabolite_name), metabolite, metabolite_name)) %>%
    select(protein, metabolite, metabolite_label, r, q, r_s, q_s, sign, abs_r)
}

edges_es_lab  <- label_edges(edges_es)
edges_top_lab <- label_edges(edges_top)

readr::write_csv(edges_es_lab,  file.path(res_dir, "block5_edges_labeled_q_le_alpha_absr.csv"))
readr::write_csv(edges_top_lab, file.path(res_dir, "block5_edges_labeled_top10k_by_q.csv"))

nodes_met <- read_csv(file.path(res_dir, "block5_nodes_metabolites.csv"), show_col_types = FALSE)
nodes_met_lab <- nodes_met %>% rename(metabolite = name) %>%
  left_join(meta_info, by = "metabolite") %>%
  mutate(name = if_else(is.na(metabolite_name), metabolite, metabolite_name)) %>%
  select(name, type)
readr::write_csv(nodes_met_lab, file.path(res_dir, "block5_nodes_metabolites_labeled.csv"))

top_by_protein <- edges_es_lab %>%
  arrange(protein, desc(abs_r)) %>%
  group_by(protein) %>%
  slice_head(n = 5) %>%
  ungroup()
readr::write_csv(top_by_protein, file.path(res_dir, "block5_top5_metabolites_per_protein_labeled.csv"))

cat("Labeling complete.\n")
