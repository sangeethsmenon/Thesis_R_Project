packages <- c("tidyverse","here","readr")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))
suppressPackageStartupMessages({ library(jsonlite) })

# ---- Load config ----
if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found. Run Block 0.")
}

# ---- Inputs ----
raw_dir  <- tryCatch(config$paths$data_raw, error = function(e) here::here("data","raw"))
res_dir  <- config$paths$results

info_path   <- file.path(raw_dir, "metabolite_info.txt")
edges_path  <- file.path(res_dir, "block5_edges_filtered_pruned.csv")   # <- use PRUNED edges
nodes_met_p <- file.path(res_dir, "block5_nodes_metabolites.csv")

stopifnot(file.exists(info_path), file.exists(edges_path), file.exists(nodes_met_p))

# Expected columns in edges: protein, metabolite, n_pair, r, p, q, r_s, p_s, q_s, sign, abs_r
edges <- readr::read_csv(edges_path, show_col_types = FALSE)

# ---- Metabolite dictionary ----
# Two columns: fid and human-readable name
meta_info <- readr::read_tsv(info_path, show_col_types = FALSE, col_types = "cc")
names(meta_info) <- c("metabolite","metabolite_name")

# ---- Apply labels to edges ----
edges_lab <- edges %>%
  left_join(meta_info, by = "metabolite") %>%
  mutate(metabolite_label = if_else(is.na(metabolite_name), metabolite, metabolite_name)) %>%
  select(protein, metabolite, metabolite_label,
         n_pair, r, p, q, r_s, p_s, q_s, sign, abs_r)

# ---- Labeled metabolite nodes (for plots/Cytoscape) ----
nodes_met <- readr::read_csv(nodes_met_p, show_col_types = FALSE) %>%
  rename(metabolite = name)

nodes_met_lab <- nodes_met %>%
  left_join(meta_info, by = "metabolite") %>%
  transmute(name = if_else(is.na(metabolite_name), metabolite, metabolite_name),
            type)

# ---- Exports ----
readr::write_csv(edges_lab,   file.path(res_dir, "block5_edges_filtered_pruned_labeled.csv"))
readr::write_csv(nodes_met_lab, file.path(res_dir, "block5_nodes_metabolites_labeled.csv"))

# ---- Quick sanity checks ----
cat("Labeling complete.\n")
cat("Unique metabolites in edges (raw IDs): ", dplyr::n_distinct(edges$metabolite), "\n")
cat("Unique metabolite labels written:      ", dplyr::n_distinct(edges_lab$metabolite_label), "\n")
cat("Unmapped fids (rows):                  ", sum(is.na(edges_lab$metabolite_label)), "\n")
