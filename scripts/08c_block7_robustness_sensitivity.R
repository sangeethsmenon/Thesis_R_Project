# --- 08c_block7_robustness_sensitivity.R --------------------------------------
# Sweeps (alpha, |r|) thresholds, rebuilds networks, compares to baseline.
# Outputs per-grid summaries + overlap metrics (edges/node overlap, NMI of modules).
# ------------------------------------------------------------------------------

packages <- c("tidyverse","here","igraph","jsonlite")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# Load config
if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found. Run Block 0.")
}
resdir <- config$paths$results; figdir <- config$paths$figures
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)

# Load a generic Block 5 edge list (prefer labeled)
cand <- c("block5_edges_labeled_q_le_alpha_absr.csv",
          "block5_edges_filtered_pruned.csv",
          "block5_edges_filtered.csv",
          "block5_edges_pearson_q_le_alpha_absr.csv")
edges_path <- file.path(resdir, cand)
edges_path <- edges_path[file.exists(edges_path)]
stopifnot(length(edges_path) >= 1)
edges_all <- readr::read_csv(edges_path[1], show_col_types = FALSE)

# Harmonize columns
nm <- names(edges_all)
if (!"protein"    %in% nm) stop("Edges need 'protein'.")
if (!"metabolite" %in% nm) stop("Edges need 'metabolite'.")
if (!"r" %in% nm) stop("Edges need 'r' (Pearson).")
if (!"q" %in% nm) {
  if ("q_value" %in% nm) edges_all <- dplyr::rename(edges_all, q = q_value)
  else if ("q_s" %in% nm) edges_all <- dplyr::rename(edges_all, q = q_s)
  else stop("Need q/q_value/q_s.")
}
if (!"abs_r" %in% nm) edges_all <- dplyr::mutate(edges_all, abs_r = abs(.data$r))

# Parameter grid
alphas <- c(0.01, 0.005, 0.001)
rmins  <- c(0.20, 0.25, 0.30)
grid   <- tidyr::crossing(alpha = alphas, min_abs_r = rmins)

# Helper: build graph + Louvain, return list
build_graph <- function(ed) {
  nodes <- dplyr::bind_rows(
    tibble::tibble(name = unique(ed$protein),   type = "protein"),
    tibble::tibble(name = unique(ed$metabolite),type = "metabolite")
  ) %>% dplyr::distinct(name, .keep_all = TRUE)
  
  g <- igraph::graph_from_data_frame(
    d = ed %>% dplyr::transmute(from = protein, to = metabolite, r = r, q = q, abs_r = abs_r),
    directed = FALSE, vertices = nodes
  )
  louv <- igraph::cluster_louvain(g)
  list(
    g = g,
    louv = louv,
    nodes_df = tibble::tibble(
      name   = V(g)$name,
      type   = ifelse(V(g)$type, "metabolite","protein"),
      module = as.integer(igraph::membership(louv))
    ),
    edges_df = igraph::as_data_frame(g, what = "edges") %>% tibble::as_tibble()
  )
}

# Baseline (from config, or default)
alpha_base <- tryCatch(config$block7$q_alpha, error=function(e) 0.01); if (is.null(alpha_base)) alpha_base <- 0.01
rmin_base  <- tryCatch(config$block7$min_abs_r, error=function(e) 0.20); if (is.null(rmin_base))  rmin_base  <- 0.20

baseline <- edges_all %>% dplyr::filter(q <= alpha_base, abs_r >= rmin_base)
stopifnot(nrow(baseline) > 0)
fit_base <- build_graph(baseline)

# Storage
sum_rows <- vector("list", nrow(grid))
memb_files <- character(0)

for (i in seq_len(nrow(grid))) {
  a <- grid$alpha[i]; rmin <- grid$min_abs_r[i]
  ed <- edges_all %>% dplyr::filter(q <= a, abs_r >= rmin)
  if (!nrow(ed)) next
  fit <- build_graph(ed)
  
  # Edge Jaccard (undirected)
  ep <- fit$edges_df %>% dplyr::transmute(e = paste(pmin(from,to), pmax(from,to), sep="||"))
  eb <- fit_base$edges_df %>% dplyr::transmute(e = paste(pmin(from,to), pmax(from,to), sep="||"))
  inter <- length(intersect(ep$e, eb$e))
  union  <- length(union(ep$e, eb$e))
  jacc   <- if (union == 0) NA_real_ else inter/union
  
  # Node overlap
  node_inter <- length(intersect(fit$nodes_df$name, fit_base$nodes_df$name))
  node_union <- length(union(fit$nodes_df$name, fit_base$nodes_df$name))
  node_jacc  <- if (node_union == 0) NA_real_ else node_inter / node_union
  
  # NMI between Louvain partitions on shared nodes (numeric via igraph::compare)
  shared <- intersect(fit$nodes_df$name, fit_base$nodes_df$name)
  nmi <- NA_real_
  if (length(shared) >= 5) {
    lab_cur  <- fit$nodes_df$module[match(shared, fit$nodes_df$name)]
    lab_base <- fit_base$nodes_df$module[match(shared, fit_base$nodes_df$name)]
    nmi <- igraph::compare(lab_cur, lab_base, method = "nmi")
    nmi <- as.numeric(nmi)
  }
  
  sum_rows[[i]] <- tibble::tibble(
    alpha = as.numeric(a),
    min_abs_r = as.numeric(rmin),
    n_edges = as.integer(igraph::gsize(fit$g)),
    n_nodes = as.integer(igraph::gorder(fit$g)),
    n_modules = as.integer(dplyr::n_distinct(fit$nodes_df$module)),
    modularity = as.numeric(igraph::modularity(fit$louv)),
    edge_jaccard_vs_baseline = as.numeric(jacc),
    node_jaccard_vs_baseline = as.numeric(node_jacc),
    nmi_vs_baseline = as.numeric(nmi)
  )
  
  # Save membership for this grid
  key <- sprintf("a%g_r%.2f", a, rmin)
  out_mem <- file.path(resdir, sprintf("block7_membership_%s.csv", gsub("[.]", "", key)))
  readr::write_csv(fit$nodes_df, out_mem)
  memb_files <- c(memb_files, out_mem)
}

summary_tbl <- dplyr::bind_rows(sum_rows) %>% dplyr::arrange(alpha, min_abs_r)
readr::write_csv(summary_tbl, file.path(resdir, "block7_sensitivity_summary.csv"))

message("Sensitivity sweep written: results/block7_sensitivity_summary.csv")
if (length(memb_files)) message("Per-grid membership CSVs:\n - ", paste(memb_files, collapse="\n - "))
