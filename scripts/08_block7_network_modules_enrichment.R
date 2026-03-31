# --- 08_block7_network_modules_enrichment.R -----------------------------------
# Builds bipartite graph from Block 5 edges, runs community detection (Louvain),
# exports node/edge tables and module summaries, and (optionally) runs GO/KEGG.
# ------------------------------------------------------------------------------

packages <- c("tidyverse","here","igraph","jsonlite")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# --- Load config
if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found. Run Block 0.")
}

dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)
dir.create(config$paths$figures, recursive = TRUE, showWarnings = FALSE)

# --- Thresholds (fallbacks if missing in config$block7)
alpha        <- tryCatch(config$block7$q_alpha, error=function(e) 0.01); if (is.null(alpha)) alpha <- 0.01
min_abs_r    <- tryCatch(config$block7$min_abs_r, error=function(e) 0.20); if (is.null(min_abs_r)) min_abs_r <- 0.20
topn_feat    <- tryCatch(config$block7$topn_features, error=function(e) 20); if (is.null(topn_feat)) topn_feat <- 20
require_pair <- tryCatch(config$block7$require_pearson_spearman_agree, error=function(e) FALSE)
if (is.null(require_pair)) require_pair <- FALSE

# --- Load a labeled Block 5 edge file (auto-detect among common names)
cand <- c("block5_edges_labeled_q_le_alpha_absr.csv",
          "block5_edges_filtered_pruned.csv",
          "block5_edges_filtered.csv",
          "block5_edges_pearson_q_le_alpha_absr.csv")
edges0_path <- file.path(config$paths$results, cand)
edges0_path <- edges0_path[file.exists(edges0_path)]
stopifnot(length(edges0_path) >= 1)
edges0 <- readr::read_csv(edges0_path[1], show_col_types = FALSE)

# Harmonize column names
nm <- names(edges0)
if (!"protein"    %in% nm) stop("Edges file must contain 'protein'.")
if (!"metabolite" %in% nm) stop("Edges file must contain 'metabolite'.")
if (!"r"          %in% nm) stop("Edges file must contain 'r' (Pearson).")
if (!"q"          %in% nm) {
  if ("q_value" %in% nm) edges0 <- dplyr::rename(edges0, q = q_value)
  else if ("q_s" %in% nm) edges0 <- dplyr::rename(edges0, q = q_s)
  else stop("Edges file must contain q/q_value/q_s.")
}
if (!"abs_r" %in% nm) edges0 <- dplyr::mutate(edges0, abs_r = abs(.data$r))

# Optional metabolite labels (if present)
lab_col <- if ("metabolite_label" %in% names(edges0)) "metabolite_label" else "metabolite"

# Optional Pearson–Spearman agreement (if your file has both)
if (require_pair && all(c("r_s","q_s") %in% names(edges0))) {
  # keep only edges that are significant/consistent in both metrics
  edges0 <- edges0 %>%
    dplyr::filter(.data$q <= alpha, abs(.data$r) >= min_abs_r,
                  .data$q_s <= alpha, abs(.data$r_s) >= min_abs_r,
                  sign(.data$r) == sign(.data$r_s))
} else {
  edges0 <- edges0 %>%
    dplyr::filter(.data$q <= alpha, .data$abs_r >= min_abs_r)
}

stopifnot(nrow(edges0) > 0)

# --- Build nodes and igraph
edges <- edges0 %>%
  dplyr::transmute(
    from = protein, to = metabolite,
    r, q, abs_r, sign = ifelse(r >= 0, "+", "-"),
    to_label = .data[[lab_col]]
  )

prot_nodes <- tibble::tibble(name = unique(edges$from), type = "protein",  display_label = unique(edges$from))
met_nodes  <- tibble::tibble(name = unique(edges$to),   type = "metabolite",
                             display_label = edges %>%
                               dplyr::distinct(to, to_label) %>%
                               { setNames(.$to_label, .$to) } %>% { .[unique(edges$to)] } %>% unname())
nodes <- dplyr::bind_rows(prot_nodes, met_nodes)

g <- igraph::graph_from_data_frame(
  d = edges %>% dplyr::select(from, to, r, q, abs_r, sign),
  directed = FALSE, vertices = nodes
)

# igraph 'type' is logical; TRUE often used for one side
V(g)$type  <- ifelse(nodes$type[match(V(g)$name, nodes$name)] == "metabolite", TRUE, FALSE)
V(g)$label <- nodes$display_label[match(V(g)$name, nodes$name)]

# --- Community detection
set.seed(42)
louv <- igraph::cluster_louvain(g)
V(g)$module_louvain <- as.integer(igraph::membership(louv))
mod_val <- igraph::modularity(louv)

# An alternative (walktrap) if you want a second view
wt <- igraph::cluster_walktrap(g)
wt_mem <- as.integer(igraph::membership(wt))
V(g)$module_walktrap <- wt_mem

# --- Export nodes/edges and summaries
nodes_df <- tibble::tibble(
  name    = V(g)$name,
  type    = ifelse(V(g)$type, "metabolite", "protein"),
  display_label = V(g)$label,
  degree  = igraph::degree(g),
  strength = igraph::strength(g, weights = abs(E(g)$r)),
  module_louvain  = V(g)$module_louvain,
  module_walktrap = V(g)$module_walktrap
)
edges_df <- igraph::as_data_frame(g, what = "edges") %>%
  tibble::as_tibble() %>%
  dplyr::rename(protein_metabolite_from = from,
                protein_metabolite_to   = to)

readr::write_csv(nodes_df, file.path(config$paths$results, "block7_network_nodes.csv"))
readr::write_csv(edges_df, file.path(config$paths$results, "block7_network_edges.csv"))

mod_summary <- nodes_df %>%
  dplyr::count(module_louvain, type, name = "n_type") %>%
  tidyr::pivot_wider(names_from = type, values_from = n_type, values_fill = 0) %>%
  dplyr::rename(n_metabolites = metabolite, n_proteins = protein) %>%
  dplyr::mutate(n_total = n_metabolites + n_proteins) %>%
  dplyr::arrange(dplyr::desc(n_total))
readr::write_csv(mod_summary, file.path(config$paths$results, "block7_module_summary.csv"))

net_summary <- tibble::tibble(
  alpha = alpha, min_abs_r = min_abs_r,
  n_edges = igraph::gsize(g), n_nodes = igraph::gorder(g),
  n_modules_louvain = dplyr::n_distinct(nodes_df$module_louvain),
  n_modules_walktrap = dplyr::n_distinct(nodes_df$module_walktrap),
  modularity_louvain = mod_val
)
readr::write_csv(net_summary, file.path(config$paths$results, "block7_network_summary.csv"))

# --- Cytoscape-friendly exports
cyto_nodes <- nodes_df %>% dplyr::transmute(
  id = name, label = display_label, type,
  degree, strength, module_louvain, module_walktrap
)
cyto_edges <- edges_df %>%
  dplyr::transmute(source = protein_metabolite_from, target = protein_metabolite_to,
                   weight = r, abs_weight = abs_r, sign = sign, q = q)
readr::write_csv(cyto_nodes, file.path(config$paths$results, "block7_cytoscape_nodes.csv"))
readr::write_csv(cyto_edges, file.path(config$paths$results, "block7_cytoscape_edges.csv"))

# --- GO / KEGG enrichment (optional; runs only if Bioconductor pkgs exist)
has_bioc <- FALSE
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  has_bioc <- TRUE
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
}
min_proteins_for_enrich <- tryCatch(config$block7$min_proteins_for_enrich, error=function(e) 10); if (is.null(min_proteins_for_enrich)) min_proteins_for_enrich <- 10

if (has_bioc) {
  dir.create(file.path(config$paths$results,"block7_enrichment_go"),   recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(config$paths$results,"block7_enrichment_kegg"), recursive = TRUE, showWarnings = FALSE)
  
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = nodes_df$name[nodes_df$type=="protein"],
    columns = c("ENTREZID","SYMBOL"),
    keytype = "SYMBOL"
  ) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
    dplyr::rename(name = SYMBOL)
  prot2entrez <- gene_map %>% dplyr::select(name, ENTREZID)
  
  modules <- sort(unique(nodes_df$module_louvain))
  enrich_go_all <- list(); enrich_kegg_all <- list()
  
  for (m in modules) {
    prots <- nodes_df %>%
      dplyr::filter(module_louvain == m, type == "protein") %>%
      dplyr::left_join(prot2entrez, by = "name") %>%
      dplyr::pull(ENTREZID) %>% unique() %>% { .[!is.na(.)] }
    if (length(prots) >= min_proteins_for_enrich) {
      ego <- tryCatch(
        clusterProfiler::enrichGO(gene = prots, OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05, readable = TRUE),
        error=function(e) NULL
      )
      ekg <- tryCatch(
        clusterProfiler::enrichKEGG(gene = prots, organism = "hsa",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05),
        error=function(e) NULL
      )
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        df_go <- tibble::as_tibble(as.data.frame(ego)) %>% dplyr::mutate(module = m, .before = 1)
        readr::write_csv(df_go, file.path(config$paths$results, "block7_enrichment_go", sprintf("go_module_%03d.csv", m)))
        enrich_go_all[[as.character(m)]] <- df_go
      }
      if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0) {
        df_kegg <- tibble::as_tibble(as.data.frame(ekg)) %>% dplyr::mutate(module = m, .before = 1)
        readr::write_csv(df_kegg, file.path(config$paths$results, "block7_enrichment_kegg", sprintf("kegg_module_%03d.csv", m)))
        enrich_kegg_all[[as.character(m)]] <- df_kegg
      }
    }
  }
  if (length(enrich_go_all))
    readr::write_csv(dplyr::bind_rows(enrich_go_all),   file.path(config$paths$results, "block7_enrichment_go_all.csv"))
  if (length(enrich_kegg_all))
    readr::write_csv(dplyr::bind_rows(enrich_kegg_all), file.path(config$paths$results, "block7_enrichment_kegg_all.csv"))
} else {
  readr::write_lines("Enrichment skipped: Bioconductor packages not available.",
                     file.path(config$paths$results, "block7_enrichment_SKIPPED.txt"))
}

message("Block 7 complete.")
message("Vertices: ", igraph::gorder(g), "  Edges: ", igraph::gsize(g))
message("Louvain modules: ", length(unique(V(g)$module_louvain)))
message("Walktrap modules: ", length(unique(V(g)$module_walktrap)))
