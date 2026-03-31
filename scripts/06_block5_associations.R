# scripts/06_block5_associations.R  (UPDATED with per-metabolite threshold option)

packages <- c("tidyverse","here")
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
  } else {
    stop("Config not found. Run Block 0.")
  }
}

# ---- Load Block 4 residualized data ----
b4_path <- file.path(config$paths$data_processed, "block4_residuals_objects.rds")
if (!file.exists(b4_path)) stop("Block 4 outputs not found: ", b4_path)
b4 <- readRDS(b4_path)

M <- as.matrix(b4$metabolome_resid_std)  # n x m  (residualized, with NAs)
P <- as.matrix(b4$proteome_resid_std)    # n x p  (residualized, with NAs)
metab_names <- b4$metabolome_features
prot_names  <- b4$proteome_features

stopifnot(nrow(M) == nrow(P))
stopifnot(ncol(M) == length(metab_names), ncol(P) == length(prot_names))

n <- nrow(M); m <- ncol(M); p <- ncol(P)

# ---- Thresholds (strict, from config) ----
thr <- config$thresholds
fget <- function(lst, nm, default) { if (is.null(lst[[nm]])) default else lst[[nm]] }

alpha            <- fget(thr, "assoc_alpha",           0.01)
min_abs_r        <- fget(thr, "assoc_min_abs_r",       0.20)
min_pair_n       <- fget(thr, "assoc_min_pair_n",      9500L)
require_overlap  <- fget(thr, "assoc_require_overlap", TRUE)
per_node_top_k   <- fget(thr, "assoc_per_node_top_k",  25L)

# NEW: per-metabolite distribution filter controls (Xingyue suggestion)
use_metab_quantile_filter <- fget(thr, "assoc_use_metab_quantile_filter", TRUE)
metab_quantile            <- fget(thr, "assoc_metab_quantile",            0.95)  # keep top 5% |r| per metabolite
min_edges_per_metab       <- fget(thr, "assoc_min_edges_per_metab",       25L)  # fallback keep at least this many
metab_quantile_use_floor  <- fget(thr, "assoc_metab_quantile_use_floor",  TRUE) # enforce global min_abs_r too

cat("Using thresholds:\n",
    "  alpha (FDR)        =", alpha, "\n",
    "  min |r|            =", min_abs_r, "\n",
    "  min pairwise N     =", min_pair_n, "\n",
    "  require overlap    =", require_overlap, "\n",
    "  per-node top-k     =", per_node_top_k, "\n",
    "  use metab-quantile =", use_metab_quantile_filter, "\n",
    "  metab quantile     =", metab_quantile, "\n",
    "  min edges/metab    =", min_edges_per_metab, "\n")

# ---- Helpers ----
zscore_mat <- function(X) {
  mu  <- colMeans(X, na.rm = TRUE)
  sdv <- apply(X, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
}

rank_zscore_mat <- function(X) {
  R <- apply(X, 2, function(x) rank(x, ties.method = "average", na.last = "keep"))
  zscore_mat(R)
}

# ---- Pairwise non-missing counts (proteins x metabolites) ----
okP <- !is.na(P)
okM <- !is.na(M)
N_pair <- crossprod(okP, okM)   # p x m
stopifnot(all(dim(N_pair) == c(p, m)))

# ---- Pearson correlations ----
ZP <- zscore_mat(P)
ZM <- zscore_mat(M)
ZP[is.na(ZP)] <- 0
ZM[is.na(ZM)] <- 0

R_num <- crossprod(ZP, ZM)         # p x m
R_den <- pmax(1, N_pair - 1)
R_pear <- R_num / R_den
R_pear[R_pear >  1] <-  1
R_pear[R_pear < -1] <- -1

df_pear <- pmax(1, N_pair - 2)
T_pear  <- R_pear * sqrt(df_pear / pmax(1e-12, 1 - R_pear^2))
P_pear  <- 2 * pt(abs(T_pear), df = df_pear, lower.tail = FALSE)
Q_pear  <- matrix(p.adjust(as.vector(P_pear), method = "BH"), nrow = p, ncol = m)

# ---- Spearman correlations ----
ZP_s <- rank_zscore_mat(P); ZP_s[is.na(ZP_s)] <- 0
ZM_s <- rank_zscore_mat(M); ZM_s[is.na(ZM_s)] <- 0

R_num_s <- crossprod(ZP_s, ZM_s)
R_spear <- R_num_s / R_den
R_spear[R_spear >  1] <-  1
R_spear[R_spear < -1] <- -1

T_spear <- R_spear * sqrt(df_pear / pmax(1e-12, 1 - R_spear^2))
P_spear <- 2 * pt(abs(T_spear), df = df_pear, lower.tail = FALSE)
Q_spear <- matrix(p.adjust(as.vector(P_spear), method = "BH"), nrow = p, ncol = m)

# ---- Build long table of edges ----
edges <- tibble(
  protein     = rep(prot_names, times = m),
  metabolite  = rep(metab_names, each  = p),
  n_pair      = as.vector(N_pair),
  r           = as.vector(R_pear),
  p           = as.vector(P_pear),
  q           = as.vector(Q_pear),
  r_s         = as.vector(R_spear),
  p_s         = as.vector(P_spear),
  q_s         = as.vector(Q_spear)
) %>%
  mutate(sign = ifelse(r >= 0, "+", "-"),
         abs_r = abs(r))

# ---- Base filtering by thresholds (your original logic) ----
edges_base <- edges %>%
  filter(n_pair >= min_pair_n,
         q <= alpha,
         abs_r >= min_abs_r,
         if (require_overlap) q_s <= alpha else TRUE)


# ============================================================
# Metabolite-wise filtering (Xingyue suggestion)
# Keep only strongest edges PER metabolite based on |r| quantile
# ============================================================

# Add these in config$thresholds if you want, otherwise defaults below
metab_q <- fget(thr, "assoc_metab_absr_quantile", 0.95)   # e.g., 0.95 = top 5%
metab_min_keep <- fget(thr, "assoc_metab_min_keep", 10L)  # always keep at least this many per metabolite

cat("Metabolite-wise filter:\n",
    "  quantile =", metab_q, "\n",
    "  min_keep =", metab_min_keep, "\n")

# For each metabolite, compute its |r| quantile cutoff (within edges_base only),
# then keep edges with |r| >= that cutoff.
edges_metab_quant <- edges_base %>%
  group_by(metabolite) %>%
  mutate(
    metab_absr_cut = quantile(abs_r, probs = metab_q, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(abs_r >= metab_absr_cut)

# Guarantee at least metab_min_keep edges per metabolite (if available)
# This prevents a metabolite from becoming "empty" due to quantile ties or small counts.
edges_metab_topN <- edges_base %>%
  group_by(metabolite) %>%
  arrange(desc(abs_r), .by_group = TRUE) %>%
  slice_head(n = metab_min_keep) %>%
  ungroup()

edges_metab_filtered <- bind_rows(edges_metab_quant, edges_metab_topN) %>%
  distinct()

cat("Edges after metabolite-wise filtering:", nrow(edges_metab_filtered), "\n")

# Write this as a NEW primary PM-edge file (for Block 5b to use)
readr::write_csv(edges_metab_filtered %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_filtered_metab_quantile.csv"))

# ------------------------------------------------------------
# NEW: Per-metabolite distribution filter (Xingyue suggestion)
# ------------------------------------------------------------
edges_metab_balanced <- edges_base
if (isTRUE(use_metab_quantile_filter)) {
  
  # Quantile cutoff per metabolite
  edges_metab_quant <- edges_base %>%
    group_by(metabolite) %>%
    mutate(metab_cut = as.numeric(stats::quantile(abs_r, probs = metab_quantile, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(metab_cut = if (isTRUE(metab_quantile_use_floor)) pmax(metab_cut, min_abs_r) else metab_cut) %>%
    filter(abs_r >= metab_cut) %>%
    select(-metab_cut)
  
  # Fallback: ensure at least N edges per metabolite by top-|r|
  edges_metab_topN <- edges_base %>%
    group_by(metabolite) %>%
    arrange(desc(abs_r), .by_group = TRUE) %>%
    slice_head(n = min_edges_per_metab) %>%
    ungroup()
  
  # Union
  edges_metab_balanced <- bind_rows(edges_metab_quant, edges_metab_topN) %>%
    distinct(protein, metabolite, .keep_all = TRUE)
  
  cat("Metab-balanced PM edges:\n",
      "  base edges      =", nrow(edges_base), "\n",
      "  balanced edges  =", nrow(edges_metab_balanced), "\n")
  
  cat("Top metabolites by edge count (balanced):\n")
  print(edges_metab_balanced %>% count(metabolite, sort = TRUE) %>% slice_head(n = 15))
}

# ---- Optional pruning: top-k per node by |r|, then union ----
# Apply to BOTH the original base and the balanced set (so you can choose later).
prune_topk <- function(df) {
  if (!is.null(per_node_top_k) && per_node_top_k > 0) {
    topP <- df %>%
      group_by(protein) %>% arrange(desc(abs_r), .by_group = TRUE) %>%
      slice_head(n = per_node_top_k) %>% ungroup()
    topM <- df %>%
      group_by(metabolite) %>% arrange(desc(abs_r), .by_group = TRUE) %>%
      slice_head(n = per_node_top_k) %>% ungroup()
    bind_rows(topP, topM) %>% distinct()
  } else {
    df
  }
}

edges_pruned          <- prune_topk(edges_base)
edges_metab_bal_pruned <- prune_topk(edges_metab_balanced)

# ---- Summary ----
summ <- tibble(
  n_samples       = n,
  n_proteins      = p,
  n_metabolites   = m,
  n_tests         = p * m,
  alpha           = alpha,
  min_abs_r       = min_abs_r,
  min_pair_n      = min_pair_n,
  require_overlap = require_overlap,
  per_node_top_k  = per_node_top_k > 0,
  metab_quantile_filter = isTRUE(use_metab_quantile_filter),
  metab_quantile  = metab_quantile,
  min_edges_per_metab = min_edges_per_metab,
  pearson_sig     = sum(edges$q <= alpha, na.rm = TRUE),
  spearman_sig    = sum(edges$q_s <= alpha, na.rm = TRUE),
  edges_base_n    = nrow(edges_base),
  edges_pruned_n  = nrow(edges_pruned),
  edges_metab_balanced_n = nrow(edges_metab_balanced),
  edges_metab_balanced_pruned_n = nrow(edges_metab_bal_pruned)
)

# ------------------------------------------------------------
# Attach human-readable metabolite names to PM edge tables
# ------------------------------------------------------------

data_raw <- config$paths$data_raw
if (is.null(data_raw)) data_raw <- here::here("data", "raw")

map_fp <- file.path(data_raw, "metabolite_info.txt")
if (!file.exists(map_fp)) {
  stop("metabolite_info.txt not found at: ", map_fp)
}

metab_map <- readr::read_tsv(map_fp, show_col_types = FALSE) %>%
  dplyr::transmute(
    metabolite = as.character(fid),
    metabolite_label = as.character(metabolite_name)
  ) %>%
  dplyr::distinct(metabolite, .keep_all = TRUE)

cat("Metabolite mapping loaded:", nrow(metab_map), "rows\n")

attach_metab_label <- function(df) {
  df %>%
    dplyr::left_join(metab_map, by = "metabolite") %>%
    dplyr::mutate(
      metabolite_label = dplyr::coalesce(metabolite_label, metabolite)
    )
}

# Apply mapping to all PM edge tables that will be written out
edges_base               <- attach_metab_label(edges_base)
edges_pruned             <- attach_metab_label(edges_pruned)
edges_metab_balanced     <- attach_metab_label(edges_metab_balanced)
edges_metab_bal_pruned   <- attach_metab_label(edges_metab_bal_pruned)





# ---- Write outputs ----
dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(summ, file.path(config$paths$results, "block5_assoc_summary.csv"))

# original outputs (unchanged)
readr::write_csv(edges_base %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_filtered.csv"))
readr::write_csv(edges_pruned %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_filtered_pruned.csv"))

# NEW outputs (Xingyue)
readr::write_csv(edges_metab_balanced %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_filtered_metab_quantile.csv"))
readr::write_csv(edges_metab_bal_pruned %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_filtered_metab_quantile_pruned.csv"))

readr::write_csv(tibble(name = prot_names,  type = "protein"),
                 file.path(config$paths$results, "block5_nodes_proteins.csv"))
readr::write_csv(tibble(name = metab_names, type = "metabolite"),
                 file.path(config$paths$results, "block5_nodes_metabolites.csv"))

saveRDS(list(
  N_pair      = N_pair,
  R_pearson   = R_pear,
  R_spearman  = R_spear,
  edges_all   = edges,
  edges_base  = edges_base,
  edges_pruned= edges_pruned,
  edges_metab_balanced = edges_metab_balanced,
  edges_metab_balanced_pruned = edges_metab_bal_pruned,
  summary     = summ,
  thresholds  = list(
    alpha=alpha, min_abs_r=min_abs_r, min_pair_n=min_pair_n,
    require_overlap=require_overlap, per_node_top_k=per_node_top_k,
    use_metab_quantile_filter=use_metab_quantile_filter,
    metab_quantile=metab_quantile,
    min_edges_per_metab=min_edges_per_metab,
    metab_quantile_use_floor=metab_quantile_use_floor
  )
), file = file.path(config$paths$data_processed, "block5_associations.rds"))

cat("Block 5 complete (with metab-quantile outputs).\n")
