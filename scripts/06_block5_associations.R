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

b4_path <- file.path(config$paths$data_processed, "block4_residuals_objects.rds")
if (!file.exists(b4_path)) stop("Block 4 outputs not found: ", b4_path)
b4 <- readRDS(b4_path)

dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)

id_cols <- c("sid","sex","age_in0")
M <- b4$metabolome_resid_std
P <- b4$proteome_resid_std
stopifnot(nrow(M) == nrow(P))
n <- nrow(M)
metab_names <- b4$metabolome_features
prot_names  <- b4$proteome_features
stopifnot(ncol(M) == length(metab_names), ncol(P) == length(prot_names))

alpha <- tryCatch(config$thresholds$assoc_alpha, error = function(e) 0.05); if (is.null(alpha)) alpha <- 0.05
min_abs_r <- tryCatch(config$thresholds$assoc_min_abs_r, error = function(e) 0.10); if (is.null(min_abs_r)) min_abs_r <- 0.10

R_pear <- crossprod(P, M) / (n - 1)
R_pear[R_pear >  1] <-  1
R_pear[R_pear < -1] <- -1
df <- n - 2
Tmat <- R_pear * sqrt(df / pmax(1e-12, 1 - R_pear^2))
Pval <- 2 * pt(abs(Tmat), df = df, lower.tail = FALSE)
Qval <- matrix(p.adjust(as.vector(Pval), method = "BH"), nrow = nrow(Pval), ncol = ncol(Pval))

edges_pear <- tibble(
  protein    = rep(prot_names, times = length(metab_names)),
  metabolite = rep(metab_names, each  = length(prot_names)),
  r          = as.vector(R_pear),
  p          = as.vector(Pval),
  q          = as.vector(Qval)
)

rank_std <- function(X){
  R <- apply(X, 2, rank, ties.method = "average")
  R <- scale(R, center = TRUE, scale = TRUE)
  R
}
M_rank <- rank_std(M)
P_rank <- rank_std(P)
R_spear <- crossprod(P_rank, M_rank) / (n - 1)
R_spear[R_spear >  1] <-  1
R_spear[R_spear < -1] <- -1
Tmat_s <- R_spear * sqrt(df / pmax(1e-12, 1 - R_spear^2))
Pval_s <- 2 * pt(abs(Tmat_s), df = df, lower.tail = FALSE)
Qval_s <- matrix(p.adjust(as.vector(Pval_s), method = "BH"), nrow = nrow(Pval_s), ncol = ncol(Pval_s))

edges_spear <- tibble(
  protein    = rep(prot_names, times = length(metab_names)),
  metabolite = rep(metab_names, each  = length(prot_names)),
  r_s        = as.vector(R_spear),
  p_s        = as.vector(Pval_s),
  q_s        = as.vector(Qval_s)
)

edges <- edges_pear %>%
  inner_join(edges_spear, by = c("protein","metabolite")) %>%
  mutate(sign = ifelse(r >= 0, "+", "-"),
         abs_r = abs(r))

edges_sig    <- edges %>% filter(q <= alpha)
edges_sig_es <- edges %>% filter(q <= alpha, abs_r >= min_abs_r)

summ <- tibble(
  n_samples     = n,
  n_proteins    = length(prot_names),
  n_metabolites = length(metab_names),
  n_tests       = nrow(edges),
  alpha         = alpha,
  min_abs_r     = min_abs_r,
  pearson_sig   = nrow(edges_sig),
  pearson_sig_es= nrow(edges_sig_es),
  spearman_sig  = sum(edges$q_s <= alpha),
  overlap_sig   = sum(edges$q <= alpha & edges$q_s <= alpha)
)

readr::write_csv(summ, file.path(config$paths$results, "block5_assoc_summary.csv"))
readr::write_csv(edges_sig %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_pearson_q_le_alpha.csv"))
readr::write_csv(edges_sig_es %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5_edges_pearson_q_le_alpha_absr.csv"))
readr::write_csv(edges %>% arrange(q) %>% slice_head(n = 10000),
                 file.path(config$paths$results, "block5_edges_top10k_by_q.csv"))
readr::write_csv(tibble(name = prot_names, type = "protein"),
                 file.path(config$paths$results, "block5_nodes_proteins.csv"))
readr::write_csv(tibble(name = metab_names, type = "metabolite"),
                 file.path(config$paths$results, "block5_nodes_metabolites.csv"))

saveRDS(list(
  R_pearson   = R_pear,
  R_spearman  = R_spear,
  edges       = edges,
  edges_sig   = edges_sig,
  edges_sig_es= edges_sig_es,
  summary     = summ
), file = file.path(config$paths$data_processed, "block5_associations.rds"))

cat("Block 5 complete.\n")
