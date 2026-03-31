# ---- Block 8: Robustness, summaries, and exportables (auto-detect columns) --

packages <- c("tidyverse","here","igraph","ggplot2")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

res_dir <- here::here("results")
fig_dir <- here::here("figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Block 7 outputs
edges_raw <- readr::read_csv(file.path(res_dir,"block7_network_edges.csv"), show_col_types = FALSE)
nodes_raw <- readr::read_csv(file.path(res_dir,"block7_network_nodes.csv"), show_col_types = FALSE)

go_all   <- suppressWarnings(tryCatch(readr::read_csv(file.path(res_dir,"block7_enrichment_go_all.csv"),   show_col_types = FALSE), error = function(e) NULL))
kegg_all <- suppressWarnings(tryCatch(readr::read_csv(file.path(res_dir,"block7_enrichment_kegg_all.csv"), show_col_types = FALSE), error = function(e) NULL))

# --- Helper: auto-detect edge columns no matter how they were named
std_edges <- function(df) {
  nm <- names(df); low <- tolower(nm)
  
  # 1) from/to candidates by name
  from_idx <- which(low %in% c("from","protein","source"))
  to_idx   <- which(low %in% c("to","metabolite","target"))
  
  # If still missing, pick first two character columns
  if (length(from_idx)==0 || length(to_idx)==0) {
    chr_cols <- which(sapply(df, function(x) is.character(x) || is.factor(x)))
    if (length(chr_cols) >= 2) {
      if (length(from_idx)==0) from_idx <- chr_cols[1]
      if (length(to_idx)==0)   to_idx   <- chr_cols[2]
    }
  }
  
  # 2) weight candidates by name, else by numeric heuristics
  weight_idx <- which(low %in% c("weight","r","pearson_r","rho","corr"))
  if (length(weight_idx)==0) {
    num_cols <- which(sapply(df, is.numeric))
    if (length(num_cols)) {
      # Prefer a column with values in [-1,1] and some negatives (looks like a correlation)
      cand <- num_cols[
        sapply(num_cols, function(i){
          x <- suppressWarnings(as.numeric(df[[i]]))
          rng <- range(x, na.rm=TRUE)
          is.finite(rng[1]) && is.finite(rng[2]) && rng[1] >= -1.05 && rng[2] <= 1.05 && any(x < 0, na.rm=TRUE)
        })
      ]
      if (length(cand)==0) {
        # fallback: any numeric with range ≤ 1
        cand <- num_cols[sapply(num_cols, function(i){
          x <- suppressWarnings(as.numeric(df[[i]])); rng <- range(x, na.rm=TRUE)
          is.finite(rng[1]) && is.finite(rng[2]) && rng[1] >= 0 && rng[2] <= 1.05
        })]
      }
      if (length(cand)) weight_idx <- cand[1]
    }
  }
  
  # 3) q-value & abs_weight (optional)
  q_idx  <- which(low %in% c("q","q_s","fdr","padj","qvalue"))
  abs_idx<- which(low %in% c("abs_weight","abs_r","abs_corr"))
  
  if (length(from_idx)==0 || length(to_idx)==0 || length(weight_idx)==0)
    stop("Edges file lacks recognizable columns for nodes and a correlation/weight. ",
         "Found columns: ", paste(nm, collapse=", "))
  
  tibble(
    from = as.character(df[[from_idx[1]]]),
    to   = as.character(df[[to_idx[1]]]),
    weight     = as.numeric(df[[weight_idx[1]]]),
    q          = if (length(q_idx))   as.numeric(df[[q_idx[1]]]) else NA_real_,
    abs_weight = if (length(abs_idx)) as.numeric(df[[abs_idx[1]]]) else abs(as.numeric(df[[weight_idx[1]]])),
    sign       = ifelse(as.numeric(df[[weight_idx[1]]]) >= 0, "+", "-")
  ) %>%
    filter(!is.na(from), !is.na(to), !is.na(weight))
}

edges0 <- std_edges(edges_raw)

# --- Standardize node table (fill missing label/type if needed)
nodes0 <- nodes_raw %>%
  mutate(
    label = dplyr::coalesce(label, name),
    type  = dplyr::coalesce(type, if_else(grepl("^fid_", name), "metabolite", "protein"))
  ) %>% distinct(name, .keep_all = TRUE)
stopifnot(all(c("name","label","type") %in% names(nodes0)))

# --- Builder for a network at (alpha, rmin)
build_mods <- function(alpha, rmin) {
  ef <- edges0 %>% filter((is.na(q) | q <= alpha), abs_weight >= rmin)
  keep <- unique(c(ef$from, ef$to))
  vf <- nodes0 %>% filter(name %in% keep)
  g <- igraph::graph_from_data_frame(d = ef %>% select(from,to,weight,abs_weight,q,sign),
                                     directed = FALSE, vertices = vf)
  
  clu <- igraph::cluster_louvain(g)
  V(g)$module <- igraph::membership(clu)
  
  nodes_tbl <- tibble::as_tibble(igraph::as_data_frame(g, what = "vertices"))
  edges_tbl <- tibble::as_tibble(igraph::as_data_frame(g, what = "edges"))
  m_from <- nodes_tbl$module[match(edges_tbl$from, nodes_tbl$name)]
  m_to   <- nodes_tbl$module[match(edges_tbl$to,   nodes_tbl$name)]
  edges_tbl <- edges_tbl %>% mutate(module_from = m_from, module_to = m_to, is_intra = module_from == module_to)
  
  list(
    g = g, nodes = nodes_tbl, edges = edges_tbl,
    summary = tibble(alpha=alpha, r_min=rmin,
                     n_edges=nrow(edges_tbl), n_nodes=nrow(nodes_tbl),
                     n_modules=dplyr::n_distinct(nodes_tbl$module),
                     modularity=igraph::modularity(clu))
  )
}

# --- Two thresholds → robustness snapshot
m1 <- build_mods(alpha = 0.01, rmin = 0.25)  # main
m2 <- build_mods(alpha = 0.05, rmin = 0.30)  # sensitivity

bind_rows(m1$summary, m2$summary) %>% readr::write_csv(file.path(res_dir,"block8_robustness_network_sizes.csv"))

stab_tab <- full_join(m1$nodes %>% select(name, module1 = module),
                      m2$nodes %>% select(name, module2 = module), by="name")
stab <- stab_tab %>% summarise(n_overlap=sum(!is.na(module1)&!is.na(module2)),
                               same_module=sum(module1==module2, na.rm=TRUE),
                               pct_same=100*same_module/n_overlap)
readr::write_csv(stab, file.path(res_dir,"block8_robustness_module_stability.csv"))

# --- Per-module summary (main threshold)
module_summary <- m1$nodes %>%
  group_by(module) %>%
  summarise(n_total=n(),
            n_metabolites=sum(type=="metabolite"),
            n_proteins=sum(type=="protein"),
            .groups="drop") %>%
  left_join(
    m1$edges %>% filter(is_intra) %>% group_by(module = module_from) %>%
      summarise(n_intra_edges=n(), mean_abs_r=mean(abs_weight, na.rm=TRUE), .groups="drop"),
    by="module"
  ) %>% arrange(module)

readr::write_csv(module_summary, file.path(res_dir,"block8_module_summary_main.csv"))

# --- Top hubs and edges
deg <- igraph::degree(m1$g)
top_hubs <- m1$nodes %>%
  mutate(degree = deg[name]) %>%
  group_by(module) %>% slice_max(order_by = degree, n = 15, with_ties = FALSE) %>% ungroup()
readr::write_csv(top_hubs, file.path(res_dir,"block8_top15_hubs_per_module.csv"))

top_edges <- m1$edges %>%
  filter(is_intra) %>% group_by(module = module_from) %>%
  slice_max(order_by = abs_weight, n = 50, with_ties = FALSE) %>% ungroup()
readr::write_csv(top_edges, file.path(res_dir,"block8_top50_edges_per_module.csv"))

# --- Compact GO/KEGG “top terms per module” (if available)
pick_top_terms <- function(df, k=10) df %>% group_by(module) %>% arrange(p.adjust, .by_group=TRUE) %>% slice_head(n=k) %>% ungroup()
if (!is.null(go_all))   readr::write_csv(pick_top_terms(go_all,   10), file.path(res_dir,"block8_go_top10_per_module.csv"))
if (!is.null(kegg_all)) readr::write_csv(pick_top_terms(kegg_all, 10), file.path(res_dir,"block8_kegg_top10_per_module.csv"))

# --- Quick figure: module sizes
p_sizes <- m1$nodes %>%
  count(module, type) %>%
  ggplot(aes(x=factor(module), y=n, fill=type)) +
  geom_col() + labs(x="Module", y="Nodes", fill="Type",
                    title="Module sizes (q ≤ 0.01, |r| ≥ 0.25)") +
  theme_minimal(base_size = 12)
ggsave(file.path(fig_dir,"block8_module_sizes.png"), p_sizes, width=7, height=4, dpi=300)

# --- Optional: sex-specific Cytoscape files from Block 5c (if present)
sx_file <- file.path(res_dir,"block5c_edges_top20k_by_pdiff.csv")
if (file.exists(sx_file)) {
  sx <- readr::read_csv(sx_file, show_col_types = FALSE)
  mk <- function(col_r, col_q, tag) {
    e <- sx %>% transmute(from = protein, to = metabolite,
                          weight = .data[[col_r]], q = .data[[col_q]],
                          abs_weight = abs(.data[[col_r]]),
                          sign = if_else(.data[[col_r]] >= 0, "+", "-")) %>%
      filter(!is.na(weight), !is.na(q), q <= 0.01, abs_weight >= 0.25)
    keep <- unique(c(e$from, e$to)); v <- nodes0 %>% filter(name %in% keep)
    readr::write_csv(v, file.path(res_dir, paste0("block8_",tag,"_cytoscape_nodes.csv")))
    readr::write_csv(e, file.path(res_dir, paste0("block8_",tag,"_cytoscape_edges.csv")))
  }
  mk("r_F","q_F","female"); mk("r_M","q_M","male")
}

cat("Block 8 complete: robustness tables, module summaries, hubs/edges, top GO/KEGG, and figure.\n")
