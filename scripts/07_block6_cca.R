# --- 07_block6_cca.R ---------------------------------------------------------
# Block 6: CCA layer
#
# This script ALWAYS runs:
#   1. Ridge-regularized CCA on complete cases  (primary, stable)
#   2. Classic CCA on same matrix               (for visualization / comparison)
#   3. Sparse CCA via PMA::CCA (if available)   (for interpretability)
#
# Design:
#   - Ridge CCA: main evidence of global cross-omic coupling.
#   - Classic CCA: check that ridge structure is consistent with standard CCA.
#   - Sparse CCA: highlight a subset of variables; compare overlaps with ridge/classic.
#
# Controlled via config$block6:
#   Common:
#     k_max
#     n_perms
#     topn
#     complete_cases
#     max_missing_frac_proteome
#     max_missing_frac_metabolome
#     min_complete_cases
#
#   Ridge:
#     ridge_lambda
#
#   Sparse (if PMA available):
#     scca_penaltyx
#     scca_penaltyy
#     scca_K  (optional; defaults to k_max)
#
# Outputs:
#   Ridge:
#     results/block6_cca_summary_ridge.csv
#     results/block6_cca_loadings_proteins_ridge.csv
#     results/block6_cca_loadings_metabolites_ridge.csv
#     results/block6_toploadings/proteins_ridge_kXX.csv
#     results/block6_toploadings/metabolites_ridge_kXX.csv
#     data/processed/block6_cca_scores_ridge.rds
#
#   Classic:
#     results/block6_cca_summary_classic.csv
#     results/block6_cca_loadings_proteins_classic.csv
#     results/block6_cca_loadings_metabolites_classic.csv
#     results/block6_toploadings/proteins_classic_kXX.csv
#     results/block6_toploadings/metabolites_classic_kXX.csv
#     data/processed/block6_cca_scores_classic.rds
#
#   Sparse (if PMA installed):
#     results/block6_scca_summary.csv
#     results/block6_scca_loadings_proteins.csv
#     results/block6_scca_loadings_metabolites.csv
#     results/block6_toploadings/proteins_scca_kXX.csv
#     results/block6_toploadings/metabolites_scca_kXX.csv
#     data/processed/block6_scca_scores.rds
#
#   Plus (if metabolite_info.txt exists):
#     *_labeled.csv versions of metabolite loadings/toploadings with human-readable labels.
# ------------------------------------------------------------------------------

packages <- c("tidyverse", "here", "matrixStats", "jsonlite")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

has_RS   <- requireNamespace("RSpectra", quietly = TRUE)
has_PMA  <- requireNamespace("PMA", quietly = TRUE)

fast_svd_vals <- function(M, k) {
  if (has_RS) RSpectra::svds(M, k = k, nu = 0, nv = 0)$d
  else svd(M, nu = 0, nv = 0)$d
}

safe_svd_uv <- function(M, k) {
  if (has_RS) {
    out <- RSpectra::svds(M, k = k)
    list(u = out$u, v = out$v, d = out$d)
  } else {
    out <- svd(M, nu = k, nv = k)
    list(
      u = out$u[, seq_len(k), drop = FALSE],
      v = out$v[, seq_len(k), drop = FALSE],
      d = out$d[seq_len(k)]
    )
  }
}

safe_names <- function(nm, prefix) {
  nm <- as.character(nm)
  bad <- which(is.na(nm) | nm == "")
  if (length(bad)) nm[bad] <- paste0(prefix, "_unnamed_", seq_along(bad))
  nm
}

top_tbl <- function(W, names_vec, side, kk, topn, method_tag) {
  tibble(name = names_vec, weight = W[, kk]) %>%
    mutate(
      abs_w = abs(weight),
      rank  = rank(-abs_w, ties.method = "first")
    ) %>%
    arrange(desc(abs_w)) %>%
    slice_head(n = topn) %>%
    mutate(
      component = kk,
      side      = side,
      method    = method_tag,
      .before   = 1
    ) %>%
    select(component, side, method, name, weight, abs_w, rank)
}

# ---- Load config ----
if (!exists("config")) {
  cfg_rds  <- here::here("results", "config.rds")
  cfg_json <- here::here("results", "config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else {
    stop("Config not found. Run Block 0 before Block 6.")
  }
}

paths <- config$paths
stopifnot(!is.null(paths$data_processed), !is.null(paths$data_raw), !is.null(paths$results))

b4_path <- file.path(paths$data_processed, "block4_residuals_objects.rds")
stopifnot(file.exists(b4_path))
b4 <- readRDS(b4_path)

X0 <- as.matrix(b4$proteome_resid_std)    # samples x proteins
Y0 <- as.matrix(b4$metabolome_resid_std)  # samples x metabolites
stopifnot(nrow(X0) == nrow(Y0))

sid <- if (!is.null(b4$covars_clean$sid)) b4$covars_clean$sid else paste0("S", seq_len(nrow(X0)))
prot_names0  <- b4$proteome_features
metab_names0 <- b4$metabolome_features

blk6 <- tryCatch(config$block6, error = function(e) NULL)
if (is.null(blk6)) blk6 <- list()

ridge_lambda <- if (!is.null(blk6$ridge_lambda)) blk6$ridge_lambda else 1e-4
k_max        <- if (!is.null(blk6$k_max))        as.integer(blk6$k_max) else 5
n_perms      <- if (!is.null(blk6$n_perms))      as.integer(blk6$n_perms) else 500
topn         <- if (!is.null(blk6$topn))         as.integer(blk6$topn) else 25
use_cc       <- if (!is.null(blk6$complete_cases)) isTRUE(blk6$complete_cases) else TRUE
max_miss_P   <- if (!is.null(blk6$max_missing_frac_proteome))   blk6$max_missing_frac_proteome else 0.05
max_miss_M   <- if (!is.null(blk6$max_missing_frac_metabolome)) blk6$max_missing_frac_metabolome else 0.05
min_cc       <- if (!is.null(blk6$min_complete_cases)) as.integer(blk6$min_complete_cases) else 2000

scca_penaltyx <- if (!is.null(blk6$scca_penaltyx)) blk6$scca_penaltyx else 0.3
scca_penaltyy <- if (!is.null(blk6$scca_penaltyy)) blk6$scca_penaltyy else 0.3
scca_K        <- if (!is.null(blk6$scca_K))        as.integer(blk6$scca_K) else k_max

message(sprintf(
  "Block 6: running ridge + classic + sparse (if PMA available). RSpectra=%s, PMA=%s",
  has_RS, has_PMA
))
message(sprintf(
  "Missingness thresholds: proteome=%.2f, metabolome=%.2f; min_complete_cases=%d",
  max_miss_P, max_miss_M, min_cc
))

# ---- Helper: filter features by missingness threshold ----
filter_by_missing <- function(X, names_vec, max_frac) {
  miss_frac <- colMeans(is.na(X))
  keep <- miss_frac <= max_frac
  list(
    X         = X[, keep, drop = FALSE],
    names     = names_vec[keep],
    keep      = keep,
    miss_frac = miss_frac
  )
}

# Try thresholds until enough complete cases
tries <- tibble(
  p = c(max_miss_P, 0.10, 0.20),
  m = c(max_miss_M, 0.10, 0.20)
)
success <- FALSE

for (i in seq_len(nrow(tries))) {
  thrP <- tries$p[i]
  thrM <- tries$m[i]
  
  fp <- filter_by_missing(X0, prot_names0, thrP)
  fm <- filter_by_missing(Y0, metab_names0, thrM)
  
  X1 <- fp$X; prot_names1 <- fp$names
  Y1 <- fm$X; metab_names1 <- fm$names
  
  nzv <- function(M) matrixStats::colSds(M, na.rm = TRUE) > 0
  if (ncol(X1)) {
    keepX <- nzv(X1); X1 <- X1[, keepX, drop = FALSE]; prot_names1 <- prot_names1[keepX]
  }
  if (ncol(Y1)) {
    keepY <- nzv(Y1); Y1 <- Y1[, keepY, drop = FALSE]; metab_names1 <- metab_names1[keepY]
  }
  
  if (ncol(X1) < 2 || ncol(Y1) < 2) {
    message(sprintf("Try %d: too few features after filtering (p1=%d, p2=%d).", i, ncol(X1), ncol(Y1)))
    next
  }
  
  if (!use_cc) stop("Block 6 requires complete_cases=TRUE (no imputation here).")
  
  ok_rows <- rowSums(is.na(cbind(X1, Y1))) == 0
  n_cc <- sum(ok_rows)
  
  message(sprintf(
    "Try %d: thresholds (P=%.2f, M=%.2f) => p1=%d, p2=%d, complete-cases n=%d",
    i, thrP, thrM, ncol(X1), ncol(Y1), n_cc
  ))
  
  if (n_cc >= min_cc) {
    X <- X1[ok_rows, , drop = FALSE]
    Y <- Y1[ok_rows, , drop = FALSE]
    sid_cc <- sid[ok_rows]
    prot_names <- prot_names1
    metab_names <- metab_names1
    chosen_thrP <- thrP
    chosen_thrM <- thrM
    success <- TRUE
    break
  }
}

if (!success) {
  stop(
    "Too few complete-case samples after feature-missingness filtering. ",
    "Tried thresholds up to 20%. Adjust thresholds or min_complete_cases."
  )
}

# Center (residuals already standardized; here we just recenter after subsetting)
X <- scale(X, center = TRUE, scale = FALSE)
Y <- scale(Y, center = TRUE, scale = FALSE)

n  <- nrow(X)
p1 <- ncol(X)
p2 <- ncol(Y)
k  <- min(k_max, p1, p2, n - 1)
if (k < 1) stop("k < 1 after dimension checks.")

dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$data_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$results, "block6_toploadings"), showWarnings = FALSE)

prot_names2  <- safe_names(prot_names,  "protein")
metab_names2 <- safe_names(metab_names, "metabolite")

# ---------------- 1) Ridge CCA (primary) ----------------
message("Block 6: running ridge CCA")

Sxx <- crossprod(X) / (n - 1)
Syy <- crossprod(Y) / (n - 1)
diag(Sxx) <- diag(Sxx) + ridge_lambda
diag(Syy) <- diag(Syy) + ridge_lambda

Rx <- chol(Sxx)
Ry <- chol(Syy)

Xw <- X %*% solve(Rx)
Yw <- Y %*% solve(Ry)

Kmat <- crossprod(Xw, Yw) / (n - 1)
sv_obs <- safe_svd_uv(Kmat, k)
cancor_ridge <- pmin(pmax(sv_obs$d, 0), 1)

A <- solve(Rx, sv_obs$u)
B <- solve(Ry, sv_obs$v)

Xs <- X %*% A
Ys <- Y %*% B
sdx <- sqrt(colSums(Xs^2) / (n - 1))
sdy <- sqrt(colSums(Ys^2) / (n - 1))
A <- sweep(A, 2, sdx, "/")
B <- sweep(B, 2, sdy, "/")
Xs <- X %*% A
Ys <- Y %*% B

set.seed(123)
perm_vals <- matrix(NA_real_, nrow = n_perms, ncol = k)
for (b in seq_len(n_perms)) {
  if (b %% 10 == 0) message("ridge perm ", b, "/", n_perms)
  idx <- sample.int(n)
  Kb <- crossprod(Xw, Yw[idx, , drop = FALSE]) / (n - 1)
  db <- fast_svd_vals(Kb, k)
  len <- min(length(db), k)
  perm_vals[b, seq_len(len)] <- db[seq_len(len)]
}
p_perm <- sapply(seq_len(k), function(j) {
  obs <- cancor_ridge[j]
  vals <- perm_vals[, j]
  vals <- vals[!is.na(vals)]
  (1 + sum(vals >= obs)) / (length(vals) + 1)
})
q_perm <- p.adjust(p_perm, "BH")

summ_ridge <- tibble(
  method                 = "ridge",
  k                      = seq_len(k),
  cancor                 = cancor_ridge,
  p_perm                 = p_perm,
  q_perm                 = q_perm,
  n                      = n,
  p_proteins             = p1,
  p_metabolites          = p2,
  ridge_lambda           = ridge_lambda,
  n_perms                = n_perms,
  complete_cases         = TRUE,
  thr_missing_proteome   = chosen_thrP,
  thr_missing_metabolome = chosen_thrM
)
readr::write_csv(summ_ridge,
                 file.path(paths$results, "block6_cca_summary_ridge.csv"))

loadings_prot_ridge <- as_tibble(A) %>%
  `colnames<-`(paste0("wX_", sprintf("%02d", seq_len(ncol(A))))) %>%
  mutate(name = prot_names2, .before = 1)
readr::write_csv(loadings_prot_ridge,
                 file.path(paths$results, "block6_cca_loadings_proteins_ridge.csv"))

loadings_met_ridge <- as_tibble(B) %>%
  `colnames<-`(paste0("wY_", sprintf("%02d", seq_len(ncol(B))))) %>%
  mutate(name = metab_names2, .before = 1)
readr::write_csv(loadings_met_ridge,
                 file.path(paths$results, "block6_cca_loadings_metabolites_ridge.csv"))

for (kk in seq_len(k)) {
  readr::write_csv(
    top_tbl(A, prot_names2,  "proteins",   kk, topn, "ridge"),
    file.path(paths$results, "block6_toploadings",
              sprintf("proteins_ridge_k%02d.csv", kk))
  )
  readr::write_csv(
    top_tbl(B, metab_names2, "metabolites",kk, topn, "ridge"),
    file.path(paths$results, "block6_toploadings",
              sprintf("metabolites_ridge_k%02d.csv", kk))
  )
}

saveRDS(
  list(X_scores = Xs, Y_scores = Ys, sample_ids = sid_cc),
  file.path(paths$data_processed, "block6_cca_scores_ridge.rds")
)

message(sprintf(
  "Ridge CCA done. n=%d, p1=%d, p2=%d, k=%d, perms=%d, RSpectra=%s",
  n, p1, p2, k, n_perms, has_RS
))

# ---------------- 2) Classic CCA ----------------
message("Block 6: running classic CCA")

kk <- min(k, p1, p2, n - 1)
if (kk >= 1) {
  Sxx_c <- crossprod(X) / (n - 1)
  Syy_c <- crossprod(Y) / (n - 1)
  diag(Sxx_c) <- diag(Sxx_c) + 1e-8
  diag(Syy_c) <- diag(Syy_c) + 1e-8
  Rx_c <- chol(Sxx_c)
  Ry_c <- chol(Syy_c)
  
  Xw_c <- X %*% solve(Rx_c)
  Yw_c <- Y %*% solve(Ry_c)
  K_c  <- crossprod(Xw_c, Yw_c) / (n - 1)
  sv_c <- safe_svd_uv(K_c, kk)
  cancor_classic <- pmin(pmax(sv_c$d, 0), 1)
  
  A_c <- solve(Rx_c, sv_c$u)
  B_c <- solve(Ry_c, sv_c$v)
  
  Xs_c <- X %*% A_c
  Ys_c <- Y %*% B_c
  sdx_c <- sqrt(colSums(Xs_c^2) / (n - 1))
  sdy_c <- sqrt(colSums(Ys_c^2) / (n - 1))
  A_c <- sweep(A_c, 2, sdx_c, "/")
  B_c <- sweep(B_c, 2, sdy_c, "/")
  Xs_c <- X %*% A_c
  Ys_c <- Y %*% B_c
  
  summ_classic <- tibble(
    method                 = "classic",
    k                      = seq_len(kk),
    cancor                 = cancor_classic,
    n                      = n,
    p_proteins             = p1,
    p_metabolites          = p2,
    complete_cases         = TRUE,
    thr_missing_proteome   = chosen_thrP,
    thr_missing_metabolome = chosen_thrM
  )
  readr::write_csv(summ_classic,
                   file.path(paths$results, "block6_cca_summary_classic.csv"))
  
  loadings_prot_classic <- as_tibble(A_c) %>%
    `colnames<-`(paste0("wX_", sprintf("%02d", seq_len(ncol(A_c))))) %>%
    mutate(name = prot_names2, .before = 1)
  readr::write_csv(loadings_prot_classic,
                   file.path(paths$results, "block6_cca_loadings_proteins_classic.csv"))
  
  loadings_met_classic <- as_tibble(B_c) %>%
    `colnames<-`(paste0("wY_", sprintf("%02d", seq_len(ncol(B_c))))) %>%
    mutate(name = metab_names2, .before = 1)
  readr::write_csv(loadings_met_classic,
                   file.path(paths$results, "block6_cca_loadings_metabolites_classic.csv"))
  
  for (kk2 in seq_len(kk)) {
    readr::write_csv(
      top_tbl(A_c, prot_names2,  "proteins",   kk2, topn, "classic"),
      file.path(paths$results, "block6_toploadings",
                sprintf("proteins_classic_k%02d.csv", kk2))
    )
    readr::write_csv(
      top_tbl(B_c, metab_names2, "metabolites",kk2, topn, "classic"),
      file.path(paths$results, "block6_toploadings",
                sprintf("metabolites_classic_k%02d.csv", kk2))
    )
  }
  
  saveRDS(
    list(X_scores = Xs_c, Y_scores = Ys_c, sample_ids = sid_cc),
    file.path(paths$data_processed, "block6_cca_scores_classic.rds")
  )
  
  message(sprintf(
    "Classic CCA done. n=%d, p1=%d, p2=%d, k=%d",
    n, p1, p2, kk
  ))
} else {
  warning("Classic CCA skipped: insufficient dimension (kk < 1).")
}

# ---------------- 3) Sparse CCA (if PMA available) ----------------
if (!has_PMA) {
  warning("PMA package not installed: skipping sparse CCA. Install PMA to enable it.")
} else {
  message("Block 6: running sparse CCA via PMA::CCA")
  
  k_scca <- max(1L, min(as.integer(scca_K), p1, p2, n - 1))
  
  scca_fit <- PMA::CCA(
    x = X,
    z = Y,
    typex = "standard",
    typez = "standard",
    penaltyx = scca_penaltyx,
    penaltyz = scca_penaltyy,
    K = k_scca
  )
  
  U <- as.matrix(scca_fit$u)[, seq_len(k_scca), drop = FALSE]
  V <- as.matrix(scca_fit$v)[, seq_len(k_scca), drop = FALSE]
  
  Xs_s <- X %*% U
  Ys_s <- Y %*% V
  cancor_scca <- vapply(seq_len(k_scca), function(j) cor(Xs_s[, j], Ys_s[, j]), numeric(1))
  
  summ_scca <- tibble(
    method                 = "sparse",
    k                      = seq_len(k_scca),
    cancor                 = cancor_scca,
    n                      = n,
    p_proteins             = p1,
    p_metabolites          = p2,
    penaltyx               = scca_penaltyx,
    penaltyy               = scca_penaltyy,
    complete_cases         = TRUE,
    thr_missing_proteome   = chosen_thrP,
    thr_missing_metabolome = chosen_thrM
  )
  readr::write_csv(summ_scca,
                   file.path(paths$results, "block6_scca_summary.csv"))
  
  loadings_prot_scca <- as_tibble(U) %>%
    `colnames<-`(paste0("u_", sprintf("%02d", seq_len(ncol(U))))) %>%
    mutate(name = prot_names2, .before = 1)
  readr::write_csv(loadings_prot_scca,
                   file.path(paths$results, "block6_scca_loadings_proteins.csv"))
  
  loadings_met_scca <- as_tibble(V) %>%
    `colnames<-`(paste0("v_", sprintf("%02d", seq_len(ncol(V))))) %>%
    mutate(name = metab_names2, .before = 1)
  readr::write_csv(loadings_met_scca,
                   file.path(paths$results, "block6_scca_loadings_metabolites.csv"))
  
  for (kk3 in seq_len(k_scca)) {
    readr::write_csv(
      top_tbl(U, prot_names2,  "proteins",   kk3, topn, "scca"),
      file.path(paths$results, "block6_toploadings",
                sprintf("proteins_scca_k%02d.csv", kk3))
    )
    readr::write_csv(
      top_tbl(V, metab_names2, "metabolites",kk3, topn, "scca"),
      file.path(paths$results, "block6_toploadings",
                sprintf("metabolites_scca_k%02d.csv", kk3))
    )
  }
  
  saveRDS(
    list(X_scores = Xs_s, Y_scores = Ys_s, sample_ids = sid_cc),
    file.path(paths$data_processed, "block6_scca_scores.rds")
  )
  
  message(sprintf(
    "Sparse CCA done. n=%d, p1=%d, p2=%d, k=%d, penalties=(%.3f, %.3f)",
    n, p1, p2, k_scca, scca_penaltyx, scca_penaltyy
  ))
}

# ---------------- Metabolite label mapping (for human-readable outputs) ------
info_path <- file.path(paths$data_raw, "metabolite_info.txt")
if (file.exists(info_path)) {
  message("Block 6: attaching metabolite labels from metabolite_info.txt")
  meta_info <- readr::read_tsv(info_path, show_col_types = FALSE)
  if (ncol(meta_info) >= 2) {
    names(meta_info)[1:2] <- c("name", "metabolite_label")
    
    # Label full metabolite loading tables
    met_load_files <- list.files(
      paths$results,
      pattern = "block6_cca_loadings_metabolites_.*\\.csv$|block6_scca_loadings_metabolites\\.csv$",
      full.names = TRUE
    )
    for (fp in met_load_files) {
      df <- readr::read_csv(fp, show_col_types = FALSE) %>%
        left_join(meta_info, by = "name") %>%
        mutate(metabolite_label = if_else(is.na(metabolite_label), name, metabolite_label))
      readr::write_csv(df, sub("\\.csv$", "_labeled.csv", fp))
    }
    
    # Label top-loading metabolite tables (ridge / classic / scca)
    met_top_files <- list.files(
      file.path(paths$results, "block6_toploadings"),
      pattern = "^metabolites_.*\\.csv$",
      full.names = TRUE
    )
    for (fp in met_top_files) {
      df <- readr::read_csv(fp, show_col_types = FALSE) %>%
        left_join(meta_info, by = "name") %>%
        mutate(metabolite_label = if_else(is.na(metabolite_label), name, metabolite_label))
      readr::write_csv(df, sub("\\.csv$", "_labeled.csv", fp))
    }
  } else {
    warning("metabolite_info.txt found but has <2 columns; skip label mapping.")
  }
} else {
  message("Block 6: metabolite_info.txt not found; keeping fid_ IDs only.")
}

message("Block 6 complete: ridge + classic + (if available) sparse CCA run.")
# ---------------------------------------------------------------------------
