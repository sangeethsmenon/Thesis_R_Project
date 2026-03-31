# --- 07b_block6_cca_completecases.R -----------------------------------------
packages <- c("tidyverse","here","matrixStats")
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

# ---- Parameters (with safe defaults) ----
cc <- config$block6
prot_min_cov   <- if (!is.null(cc$prot_min_cov))   cc$prot_min_cov   else 0.99  # ≥99% non-missing
metab_min_cov  <- if (!is.null(cc$metab_min_cov))  cc$metab_min_cov  else 0.99
max_proteins   <- if (!is.null(cc$max_proteins))   cc$max_proteins   else 1500  # cap to keep CCA light
max_metabolites<- if (!is.null(cc$max_metabolites))cc$max_metabolites else 250
ridge_lambda   <- if (!is.null(cc$ridge_lambda))   cc$ridge_lambda   else 1e-4
k_max          <- if (!is.null(cc$k_max))          cc$k_max          else 10
n_perms        <- if (!is.null(cc$n_perms))        cc$n_perms        else 200
topn           <- if (!is.null(cc$topn))           cc$topn           else 25
set.seed(if (!is.null(config$seed)) config$seed else 42)

# ---- Load Block 4 residuals (still contain NA by design) ----
b4_path <- file.path(config$paths$data_processed, "block4_residuals_objects.rds")
stopifnot(file.exists(b4_path))
b4 <- readRDS(b4_path)

X_full <- as.matrix(b4$proteome_resid_std)     # n x p_prot  (may contain NAs)
Y_full <- as.matrix(b4$metabolome_resid_std)   # n x p_meta  (may contain NAs)
stopifnot(nrow(X_full) == nrow(Y_full))
n_all <- nrow(X_full)

prot_names  <- b4$proteome_features
metab_names <- b4$metabolome_features

# ---- 1) Choose high-coverage features ----
cov_prot  <- 1 - (colSums(is.na(X_full)) / n_all)
cov_metab <- 1 - (colSums(is.na(Y_full)) / n_all)

prot_keep  <- which(cov_prot  >= prot_min_cov)
metab_keep <- which(cov_metab >= metab_min_cov)

# optional: within the high-coverage set, keep the best up to a cap
if (length(prot_keep)  > max_proteins)
  prot_keep  <- prot_keep[order(cov_prot[prot_keep], decreasing = TRUE)][seq_len(max_proteins)]
if (length(metab_keep) > max_metabolites)
  metab_keep <- metab_keep[order(cov_metab[metab_keep], decreasing = TRUE)][seq_len(max_metabolites)]

X <- X_full[, prot_keep, drop = FALSE]
Y <- Y_full[, metab_keep, drop = FALSE]
p1 <- ncol(X); p2 <- ncol(Y)

# ---- 2) Keep only samples with no NA across the chosen features (complete cases) ----
ok_samples <- rowSums(is.na(X)) == 0 & rowSums(is.na(Y)) == 0
X <- X[ok_samples, , drop = FALSE]
Y <- Y[ok_samples, , drop = FALSE]
n <- nrow(X)
stopifnot(n == nrow(Y))

message("Complete-case CCA: kept ", n, " samples, ",
        p1, " proteins, ", p2, " metabolites.")

if (n < 200) warning("Very few complete-case samples left (", n, "). Consider relaxing coverage thresholds.")

# ---- 3) Standardize columns (z-score) just in case ----
zscale <- function(M){
  mu <- colMeans(M); sdv <- matrixStats::colSds(M)
  sdv[sdv == 0] <- 1
  sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
}
X <- zscale(X); Y <- zscale(Y)

# ---- 4) Ridge-regularized CCA on complete cases (no NA now) ----
k <- min(k_max, p1, p2)

Sxx <- crossprod(X) / (n - 1)
Syy <- crossprod(Y) / (n - 1)
diag(Sxx) <- diag(Sxx) + ridge_lambda
diag(Syy) <- diag(Syy) + ridge_lambda
Rx <- chol(Sxx)
Ry <- chol(Syy)

Xw <- X %*% solve(Rx)
Yw <- Y %*% solve(Ry)

Kmat <- crossprod(Xw, Yw) / (n - 1)
sv  <- svd(Kmat, nu = k, nv = k)
cancor <- pmin(pmax(sv$d[seq_len(k)], 0), 1)
Ux <- sv$u[, seq_len(k), drop = FALSE]
Vy <- sv$v[, seq_len(k), drop = FALSE]

A <- solve(Rx, Ux)
B <- solve(Ry, Vy)

# scale canonical variates to unit variance
Xs <- X %*% A
Ys <- Y %*% B
sdx <- sqrt(colSums(Xs^2) / (n - 1))
sdy <- sqrt(colSums(Ys^2) / (n - 1))
A <- sweep(A, 2, sdx, "/")
B <- sweep(B, 2, sdy, "/")
Xs <- X %*% A
Ys <- Y %*% B

# ---- 5) Permutation p values (optional but useful) ----
perm_vals <- matrix(NA_real_, nrow = n_perms, ncol = k)
for (b in seq_len(n_perms)) {
  idx <- sample.int(n)
  Kb <- crossprod(Xw, Yw[idx, , drop = FALSE]) / (n - 1)
  svb <- svd(Kb, nu = 0, nv = 0)
  db  <- svb$d
  len <- min(length(db), k)
  perm_vals[b, seq_len(len)] <- db[seq_len(len)]
}
p_perm <- sapply(seq_len(k), function(j) {
  obs  <- cancor[j]
  vals <- perm_vals[, j]; vals <- vals[!is.na(vals)]
  (1 + sum(vals >= obs)) / (length(vals) + 1)
})
q_perm <- p.adjust(p_perm, method = "BH")

# ---- 6) Save outputs ----
res_dir <- config$paths$results
dp_dir  <- config$paths$data_processed
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dp_dir,  recursive = TRUE, showWarnings = FALSE)

summ <- tibble(
  k = seq_len(k),
  cancor = cancor,
  p_perm = p_perm,
  q_perm = q_perm,
  n_complete = n,
  p_proteins = p1,
  p_metabolites = p2,
  prot_min_cov = prot_min_cov,
  metab_min_cov = metab_min_cov,
  ridge_lambda = ridge_lambda,
  n_perms = n_perms
)
readr::write_csv(summ, file.path(res_dir, "block6_cca_summary_completecases.csv"))

safe_names <- function(nm, prefix){
  nm <- as.character(nm)
  bad <- which(is.na(nm) | nm == "")
  if (length(bad)) nm[bad] <- paste0(prefix, "_unnamed_", seq_along(bad))
  nm
}
prot_names2  <- safe_names(colnames(X_full)[prot_keep],  "protein")
metab_names2 <- safe_names(colnames(Y_full)[metab_keep], "metabolite")

loadings_prot <- as_tibble(A)
colnames(loadings_prot) <- paste0("wX_", sprintf("%02d", seq_len(ncol(A))))
loadings_prot <- loadings_prot %>% mutate(name = prot_names2, .before = 1)
readr::write_csv(loadings_prot, file.path(res_dir, "block6_cca_loadings_proteins_completecases.csv"))

loadings_met <- as_tibble(B)
colnames(loadings_met) <- paste0("wY_", sprintf("%02d", seq_len(ncol(B))))
loadings_met <- loadings_met %>% mutate(name = metab_names2, .before = 1)
readr::write_csv(loadings_met, file.path(res_dir, "block6_cca_loadings_metabolites_completecases.csv"))

dir.create(file.path(res_dir, "block6_toploadings_completecases"), showWarnings = FALSE)
top_tbl <- function(W, names_vec, side, kk, topn) {
  tibble(name = names_vec, weight = W[, kk]) %>%
    mutate(abs_w = abs(weight), rank = rank(-abs_w, ties.method = "first")) %>%
    arrange(desc(abs_w)) %>% slice_head(n = topn) %>%
    mutate(component = kk, side = side, .before = 1) %>%
    select(component, side, name, weight, abs_w, rank)
}
for (kk in seq_len(k)) {
  readr::write_csv(top_tbl(A, prot_names2,  "proteins",   kk, topn),
                   file.path(res_dir, "block6_toploadings_completecases", sprintf("proteins_k%02d.csv", kk)))
  readr::write_csv(top_tbl(B, metab_names2, "metabolites",kk, topn),
                   file.path(res_dir, "block6_toploadings_completecases", sprintf("metabolites_k%02d.csv", kk)))
}

# matrices you can reuse for UMAP complete-case runs
saveRDS(list(X_cc = X, Y_cc = Y, sample_mask = ok_samples),
        file = file.path(dp_dir, "block6_completecase_matrices.rds"))

cat("Complete-case CCA finished. Kept ", n, " samples.\n")
# ---------------------------------------------------------------------------
