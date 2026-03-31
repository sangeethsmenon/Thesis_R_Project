# ===============================================================
# 09_block8_robustness_protein_first.R
#
# Robustness & Sensitivity Analyses (Protein-first pipeline)
#
# A) PM association robustness
# B) PP network robustness
# C) Cross-omic dominance stability
#
# Output:
#  - CSV summary tables
#  - One figure: dominant module fraction across PM thresholds
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(jsonlite)
  library(igraph)
})

# -----------------------------
# Load config
# -----------------------------
if (!exists("config")) {
  cfg_rds  <- here("results", "config.rds")
  cfg_json <- here("results", "config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found.")
}

resdir <- config$paths$results
figdir <- config$paths$figures
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load MAIN results (baseline)
# -----------------------------
pm_main_path <- file.path(resdir, "block5_edges_filtered_metab_quantile.csv")
pp_main_path <- file.path(resdir, "block5b_pp_edges.csv")
prot_mod_path <- file.path(resdir, "block5b_protein_modules.csv")

stopifnot(
  file.exists(pm_main_path),
  file.exists(pp_main_path),
  file.exists(prot_mod_path)
)

pm_main   <- readr::read_csv(pm_main_path, show_col_types = FALSE)
pp_main   <- readr::read_csv(pp_main_path, show_col_types = FALSE)
prot_mod  <- readr::read_csv(prot_mod_path, show_col_types = FALSE)

# -----------------------------
# Helper: dominant module fraction
# -----------------------------
compute_dominant_fraction <- function(pm_edges, prot_mod_tbl) {
  
  pm_mapped <- pm_edges %>%
    inner_join(prot_mod_tbl, by = c("protein" = "protein"))
  
  if (nrow(pm_mapped) == 0) return(NA_real_)
  
  dominant_tbl <- pm_mapped %>%
    count(module_protein, name = "pm_edges") %>%
    arrange(desc(pm_edges)) %>%
    slice(1)
  
  dominant_edges <- dominant_tbl$pm_edges
  dominant_edges / nrow(pm_mapped)
}

# ===============================================================
# A) PM ASSOCIATION ROBUSTNESS
# ===============================================================

message("Running PM robustness checks...")

pm_threshold_grid <- expand_grid(
  min_abs_r = c(0.15, 0.20, 0.25),
  metab_quantile = c(0.90, 0.95, 0.98),
  min_edges_per_metab = c(15, 25, 35)
)

pm_robust_results <- pm_threshold_grid %>%
  mutate(
    result = pmap(
      list(min_abs_r, metab_quantile, min_edges_per_metab),
      function(r_cut, q_cut, min_deg) {
        
        pm_tmp <- pm_main %>%
          dplyr::filter(abs_r >= r_cut)
        
        # metabolite-wise quantile filter
        pm_tmp <- pm_tmp %>%
          dplyr::group_by(metabolite_label) %>%
          dplyr::filter(abs_r >= quantile(abs_r, q_cut)) %>%
          dplyr::ungroup()
        
        
        
        # min edges per metabolite
        pm_tmp <- pm_tmp %>%
          dplyr::group_by(metabolite_label) %>%
          dplyr::filter(dplyr::n() >= min_deg) %>%
          dplyr::ungroup()
        
        if (nrow(pm_tmp) == 0) {
          return(tibble(
            n_edges = 0,
            n_proteins = 0,
            n_metabolites = 0,
            dominant_frac = NA_real_
          ))
        }
        
        tibble(
          n_edges = nrow(pm_tmp),
          n_proteins = n_distinct(pm_tmp$protein),
          n_metabolites = n_distinct(pm_tmp$metabolite_label),
          dominant_frac = compute_dominant_fraction(pm_tmp, prot_mod)
        )
      }
    )
  ) %>%
  unnest(result)

write_csv(
  pm_robust_results,
  file.path(resdir, "robustness_pm_thresholds_summary.csv")
)

# ===============================================================
# B) PP NETWORK ROBUSTNESS
# ===============================================================

message("Running PP robustness checks...")

pp_thresholds <- expand_grid(
  r_pp = c(0.25, 0.30),
  q_pp = c(0.01, 0.05)
)

pp_robust_results <- pp_thresholds %>%
  mutate(
    result = pmap(
      list(r_pp, q_pp),
      function(r_cut, q_cut) {
        
        pp_tmp <- pp_main %>%
          dplyr::filter(abs_r >= r_cut, q <= q_cut)
        
        if (nrow(pp_tmp) == 0) {
          return(tibble(
            n_edges = 0,
            n_nodes = 0,
            n_modules = 0
          ))
        }
        
        g <- graph_from_data_frame(
          pp_tmp %>% transmute(from = protein1, to = protein2, weight = abs_r),
          directed = FALSE
        )
        
        cl <- cluster_louvain(g, weights = E(g)$weight)
        
        tibble(
          n_edges = gsize(g),
          n_nodes = gorder(g),
          n_modules = length(unique(membership(cl)))
        )
      }
    )
  ) %>%
  unnest(result)

write_csv(
  pp_robust_results,
  file.path(resdir, "robustness_pp_thresholds_summary.csv")
)



message("Robustness analysis complete.")
