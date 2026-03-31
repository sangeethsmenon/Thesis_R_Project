# --- 08_block7_network.R ------------------------------------------------------
# Block 7: Build bipartite network from Block 5 edges and detect communities
# Exports CSVs + figures into config$paths$results / config$paths$figures
# ------------------------------------------------------------------------------

# helper first (used early)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Packages ----
pkgs <- c("tidyverse","here","igraph","jsonlite")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Load config ----
if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds))       config <- readRDS(cfg_rds)
  else if (file.exists(cfg_json)) config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  else stop("Config not found. Run Block 0.")
}
paths <- config$paths
stopifnot(!is.null(paths$results), !is.null(paths$figures))

dir.create(paths$results,  recursive = TRUE, showWarnings = FALSE)
dir.create(paths$figures,  recursive = TRUE, showWarnings = FALSE)

# ---- Locate Block 5 output (prefer pruned) ----
cand <- c(file.path(paths$results, "block5_edges_filtered_pruned.csv"),
          file.path(paths$results, "block5_edges_filtered.csv"))
edge_fp <- cand[file.exists(cand)]
if (!length(edge_fp)) stop("Block 5 edges not found in results/.")
edge_fp <- edge_fp[1]
message("Using edges file: ", edge_fp)

edges_raw <- suppressMessages(readr::read_csv(edge_fp, show_col_types = FALSE))

# ---- Heuristic column mapping (pick ONE per role) ----
cols <- names(edges_raw)
pick_col <- function(candidates, allow_prefix = TRUE) {
  hit <- intersect(candidates, cols)
  if (length(hit)) return(hit[1])
  if (allow_prefix) {
    for (p in candidates) {
      m <- cols[startsWith(cols, p)]
      if (length(m)) return(m[1])
    }
  }
  NA_character_
}

prot_col <- pick_col(c("protein","prot","x_name","a_name","p_protein","p"))
met_col  <- pick_col(c("metabolite","metab","y_name","b_name","m_metabolite","m"))
r_col    <- pick_col(c("r","pearson_r","corr","cor","r_abs","abs_r","abs.r"))
q_col    <- pick_col(c("q","q_s","fdr","bh_q","q_value","qval"))
n_col    <- pick_col(c("Npair","n_pair","npair","N_pair","n"))

if (is.na(prot_col) || is.na(met_col)) {
  stop("Could not identify protein/metabolite columns. Found: ", paste(cols, collapse = ", "))
}
if (is.na(r_col)) stop("Could not identify correlation column (r).")

edges <- edges_raw %>%
  transmute(
    protein    = as.character(.data[[prot_col]]),
    metabolite = as.character(.data[[met_col]]),
    r          = as.numeric(.data[[r_col]]),
    q          = if (!is.na(q_col)) as.numeric(.data[[q_col]]) else NA_real_,
    Npair      = if (!is.na(n_col)) as.integer(.data[[n_col]]) else NA_integer_
  ) %>%
  mutate(weight = abs(r)) %>%
  filter(!is.na(protein), !is.na(metabolite), !is.na(weight)) %>%
  distinct(protein, metabolite, .keep_all = TRUE)

stopifnot(nrow(edges) > 0)

# ---- Optional metabolite label mapping ----
label_fp <- file.path(paths$data_raw %||% "data/raw", "metabolite_info.txt")
if (file.exists(label_fp)) {
  lab <- suppressMessages(readr::read_tsv(label_fp, col_names = FALSE, show_col_types = FALSE))
  if (ncol(lab) >= 2) {
    names(lab)[1:2] <- c("metabolite","metabolite_label")
    lab <- lab %>% transmute(metabolite = as.character(metabolite),
                             metabolite_label = as.character(metabolite_label))
    edges <- edges %>%
      left_join(lab, by = "metabolite") %>%
      mutate(metabolite_disp = if_else(!is.na(metabolite_label), metabolite_label, metabolite))
  } else {
    edges <- edges %>% mutate(metabolite_disp = metabolite)
  }
} else {
  edges <- edges %>% mutate(metabolite_disp = metabolite)
}

# ---- Build bipartite graph ----
prot_nodes <- edges %>% distinct(name = protein)     %>% mutate(type = TRUE)
met_nodes  <- edges %>% distinct(name = metabolite)  %>% mutate(type = FALSE)
vertices <- bind_rows(prot_nodes, met_nodes) %>% distinct(name, .keep_all = TRUE)

g <- graph_from_data_frame(
  d = edges %>% transmute(from = protein, to = metabolite, weight = weight),
  directed = FALSE,
  vertices = vertices
)
stopifnot(!is.null(V(g)$type))

# degrees/strengths
V(g)$degree   <- degree(g)
V(g)$strength <- strength(g, weights = E(g)$weight)

# ---- Community detection ----
set.seed(42)
louv <- cluster_louvain(g, weights = E(g)$weight)
V(g)$module_louvain <- membership(louv)

wt <- cluster_walktrap(g, weights = E(g)$weight)
V(g)$module_walktrap <- membership(wt)

# ---- Export node & module tables ----
nodes_df <- tibble(
  name     = V(g)$name,
  type     = ifelse(V(g)$type, "protein", "metabolite"),
  degree   = V(g)$degree,
  strength = V(g)$strength,
  module_louvain = V(g)$module_louvain,
  module_walktrap = V(g)$module_walktrap
)

readr::write_csv(nodes_df, file.path(paths$results, "block7_nodes_louvain.csv"))
readr::write_csv(nodes_df, file.path(paths$results, "block7_nodes_walktrap.csv"))

summarise_modules <- function(nodes, label) {
  nodes %>%
    group_by(.data[[label]]) %>%
    summarise(
      module  = cur_group()[[1]],
      size    = n(),
      n_prot  = sum(type == "protein"),
      n_metab = sum(type == "metabolite"),
      .groups = "drop"
    ) %>%
    arrange(desc(size))
}
mod_louv <- summarise_modules(nodes_df, "module_louvain")
mod_wt   <- summarise_modules(nodes_df, "module_walktrap")

readr::write_csv(mod_louv, file.path(paths$results, "block7_modules_louvain.csv"))
readr::write_csv(mod_wt,   file.path(paths$results, "block7_modules_walktrap.csv"))

# ---- Annotate edges with community labels ----
memb_louv <- nodes_df %>% select(name, module_louvain)
memb_wt   <- nodes_df %>% select(name, module_walktrap)

edges_annot <- edges %>%
  select(protein, metabolite, metabolite_disp, r, q, Npair, weight) %>%
  left_join(memb_louv %>% rename(protein = name, modP_louvain = module_louvain), by = "protein") %>%
  left_join(memb_louv %>% rename(metabolite = name, modM_louvain = module_louvain), by = "metabolite") %>%
  left_join(memb_wt   %>% rename(protein = name, modP_walktrap = module_walktrap), by = "protein") %>%
  left_join(memb_wt   %>% rename(metabolite = name, modM_walktrap = module_walktrap), by = "metabolite") %>%
  mutate(
    intra_louvain  = modP_louvain  == modM_louvain,
    intra_walktrap = modP_walktrap == modM_walktrap
  )
readr::write_csv(edges_annot, file.path(paths$results, "block7_edges_annotated.csv"))

# ---- Plots ----
# Make top-N constants BEFORE slice_* calls
topN <- 20L
nP <- sum(nodes_df$type == "protein")
nM <- sum(nodes_df$type == "metabolite")
takeP <- min(topN, nP)
takeM <- min(topN, nM)

# 1) Module sizes (Louvain)
p1 <- mod_louv %>%
  ggplot(aes(x = reorder(factor(module), -size), y = size)) +
  geom_col() +
  labs(x = "Module (Louvain)", y = "Size", title = "Block 7: Module sizes (Louvain)") +
  theme_minimal(base_size = 12)
ggsave(file.path(paths$figures, "block7_module_sizes_louvain.png"), p1, width = 8, height = 4, dpi = 200)

# 2) Module sizes (Walktrap)
p2 <- mod_wt %>%
  ggplot(aes(x = reorder(factor(module), -size), y = size)) +
  geom_col() +
  labs(x = "Module (Walktrap)", y = "Size", title = "Block 7: Module sizes (Walktrap)") +
  theme_minimal(base_size = 12)
ggsave(file.path(paths$figures, "block7_module_sizes_walktrap.png"), p2, width = 8, height = 4, dpi = 200)

# 3) Top-20 degree by type (use constants, not n())
topP <- nodes_df %>%
  filter(type == "protein") %>%
  slice_max(order_by = degree, n = takeP, with_ties = FALSE) %>%
  arrange(degree)

p3 <- topP %>%
  ggplot(aes(x = degree, y = reorder(name, degree))) +
  geom_col() +
  labs(y = "Protein", x = "Degree", title = "Block 7: Top-20 protein degrees") +
  theme_minimal(base_size = 11)
ggsave(file.path(paths$figures, "block7_degree_top20_proteins.png"), p3, width = 7, height = 6, dpi = 200)

topM <- nodes_df %>%
  filter(type == "metabolite") %>%
  slice_max(order_by = degree, n = takeM, with_ties = FALSE) %>%
  arrange(degree)

p4 <- topM %>%
  ggplot(aes(x = degree, y = reorder(name, degree))) +
  geom_col() +
  labs(y = "Metabolite", x = "Degree", title = "Block 7: Top-20 metabolite degrees") +
  theme_minimal(base_size = 11)
ggsave(file.path(paths$figures, "block7_degree_top20_metabolites.png"), p4, width = 7, height = 6, dpi = 200)

# ---- Console summary ----
message("Block 7 complete.")
message("Vertices: ", gorder(g), "  Edges: ", gsize(g))
message("Louvain modules: ", length(unique(V(g)$module_louvain)))
message("Walktrap modules: ", length(unique(V(g)$module_walktrap)))
# ------------------------------------------------------------------------------
