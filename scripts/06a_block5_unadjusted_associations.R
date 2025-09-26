packages <- c("tidyverse","here")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos="https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))
suppressPackageStartupMessages({ library(jsonlite) })

if (!exists("config")) {
  cfg_rds <- here::here("results","config.rds"); cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) config <- readRDS(cfg_rds) else if (file.exists(cfg_json)) config <- jsonlite::read_json(cfg_json, simplifyVector=TRUE) else stop("Config not found.")
}

b3_path <- file.path(config$paths$data_processed, "block3_scaled_pca_outliers.rds")
stopifnot(file.exists(b3_path))
b3 <- readRDS(b3_path)

dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)

id_cols <- c("sid","sex","age_in0")
M_raw <- b3$metabolome_clean %>% select(-all_of(id_cols))
P_raw <- b3$proteome_clean   %>% select(-all_of(id_cols))
stopifnot(nrow(M_raw)==nrow(P_raw))
n <- nrow(M_raw)
metab_names <- colnames(M_raw)
prot_names  <- colnames(P_raw)

scale_cols <- function(M){M <- as.matrix(M); mu <- colMeans(M); sdv <- apply(M,2,sd); sdv[sdv==0] <- 1; sweep(sweep(M,2,mu,"-"),2,sdv,"/")}
M <- scale_cols(M_raw); P <- scale_cols(P_raw)

alpha <- tryCatch(config$thresholds$assoc_alpha, error=function(e) 0.05); if (is.null(alpha)) alpha <- 0.05
min_abs_r <- tryCatch(config$thresholds$assoc_min_abs_r, error=function(e) 0.10); if (is.null(min_abs_r)) min_abs_r <- 0.10

R_pear <- crossprod(P, M) / (n - 1); R_pear[R_pear>1] <- 1; R_pear[R_pear< -1] <- -1
df <- n - 2
Tmat <- R_pear * sqrt(df / pmax(1e-12, 1 - R_pear^2))
Pval <- 2 * pt(abs(Tmat), df = df, lower.tail = FALSE)
Qval <- matrix(p.adjust(as.vector(Pval), "BH"), nrow = nrow(Pval), ncol = ncol(Pval))

edges_pear <- tibble(
  protein    = rep(prot_names, times = length(metab_names)),
  metabolite = rep(metab_names, each  = length(prot_names)),
  r          = as.vector(R_pear),
  p          = as.vector(Pval),
  q          = as.vector(Qval)
)

rank_std <- function(X){R <- apply(X,2,rank,ties.method="average"); scale(R,center=TRUE,scale=TRUE)}
M_rank <- rank_std(M); P_rank <- rank_std(P)
R_spear <- crossprod(P_rank, M_rank) / (n - 1); R_spear[R_spear>1] <- 1; R_spear[R_spear< -1] <- -1
Tmat_s <- R_spear * sqrt(df / pmax(1e-12, 1 - R_spear^2))
Pval_s <- 2 * pt(abs(Tmat_s), df = df, lower.tail = FALSE)
Qval_s <- matrix(p.adjust(as.vector(Pval_s), "BH"), nrow = nrow(Pval_s), ncol = ncol(Pval_s))

edges <- edges_pear %>% bind_cols(
  tibble(r_s = as.vector(R_spear), p_s = as.vector(Pval_s), q_s = as.vector(Qval_s))
) %>% mutate(sign = ifelse(r>=0,"+","-"), abs_r = abs(r))

edges_sig    <- edges %>% filter(q <= alpha)
edges_sig_es <- edges %>% filter(q <= alpha, abs_r >= min_abs_r)

summ <- tibble(
  view = "unadjusted",
  n_samples = n,
  n_proteins = length(prot_names),
  n_metabolites = length(metab_names),
  n_tests = nrow(edges),
  alpha = alpha,
  min_abs_r = min_abs_r,
  pearson_sig = nrow(edges_sig),
  pearson_sig_es = nrow(edges_sig_es),
  spearman_sig = sum(edges$q_s <= alpha),
  overlap_sig = sum(edges$q <= alpha & edges$q_s <= alpha)
)

readr::write_csv(summ, file.path(config$paths$results, "block5a_unadjusted_assoc_summary.csv"))
readr::write_csv(edges_sig %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5a_edges_pearson_q_le_alpha.csv"))
readr::write_csv(edges_sig_es %>% arrange(q, desc(abs_r)),
                 file.path(config$paths$results, "block5a_edges_pearson_q_le_alpha_absr.csv"))
readr::write_csv(edges %>% arrange(q) %>% slice_head(n=10000),
                 file.path(config$paths$results, "block5a_edges_top10k_by_q.csv"))

saveRDS(list(edges=edges, summary=summ),
        file = file.path(config$paths$data_processed, "block5a_unadjusted_associations.rds"))

cat("Block 5A (unadjusted) complete.\n")
