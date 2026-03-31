# --- 04b_block3_scaling_pca_outliers_after_samplefilter.R --------------------
# Auto-k PCA (NIPALS) + Mahalanobis outlier flagging after sample-missingness filter
# - Reads thresholds from results/config.json (with sensible defaults)
# - Saves cleaned tables and summary compatible with downstream blocks

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(jsonlite)
})

# ---- Packages for NIPALS PCA (handles NA) ----
if (!requireNamespace("pcaMethods", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("pcaMethods", ask = FALSE, update = FALSE)
}
library(pcaMethods)

# ---- Load config ----
cfg_rds  <- here("results","config.rds")
cfg_json <- here("results","config.json")
if (file.exists(cfg_rds)) {
  config <- readRDS(cfg_rds)
} else if (file.exists(cfg_json)) {
  config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
} else {
  stop("Config not found. Run 00_env_setup.R first.")
}

thr <- config$thresholds %||% list()

get_thr <- function(key, default) {
  if (!is.null(thr[[key]])) thr[[key]] else default
}

# Block 3 thresholds (with safe defaults)
pca_target_metab <- as.numeric(get_thr("block3_pca_cumvar_target_metabolome", 0.80))
pca_target_prot  <- as.numeric(get_thr("block3_pca_cumvar_target_proteome",   0.70))
pca_min_k        <- as.integer(get_thr("block3_pca_min_k", 3))
pca_max_k        <- as.integer(get_thr("block3_pca_max_k", 10))
chi_q            <- as.numeric(get_thr("block3_chi2_quantile", 0.999))

# ---- Paths ----
raw_dir  <- config$paths$data_raw       %||% here("data","raw")
proc_dir <- config$paths$data_processed %||% here("data","processed")
res_dir  <- config$paths$results        %||% here("results")
fig_dir  <- config$paths$figures        %||% here("figures")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# These files must exist (created by Block 1/2 filtering)
proteome_path   <- file.path(proc_dir, "proteome_qc.tsv")
metabolome_path <- file.path(proc_dir, "metabolome_qc.tsv")
stopifnot(file.exists(proteome_path), file.exists(metabolome_path))

proteome   <- readr::read_tsv(proteome_path,   show_col_types = FALSE)
metabolome <- readr::read_tsv(metabolome_path, show_col_types = FALSE)
stopifnot(all(metabolome$sid == proteome$sid))

id_cols <- c("sid","sex","age_in0")
covars  <- metabolome %>% select(all_of(id_cols))

# Matrices (keep columns with at least one finite value)
drop_all_na_cols <- function(M) {
  keep <- colSums(is.finite(M)) > 0
  M[, keep, drop = FALSE]
}
metab_M <- metabolome %>% select(-all_of(id_cols)) %>% as.matrix() %>% drop_all_na_cols()
prot_M  <- proteome   %>% select(-all_of(id_cols)) %>% as.matrix() %>% drop_all_na_cols()

n_before <- nrow(metabolome)

# ---- Helper: run NIPALS PCA + choose k by cumulative variance target ----
pcak_choose_k <- function(M, target, min_k, max_k) {
  # NIPALS with UV scaling + centering (scores tolerate NAs)
  fit <- pcaMethods::pca(M, method = "nipals",
                         nPcs = max_k, center = TRUE, scale = "uv")
  
  # Get per-component explained variance fraction
  get_var <- function(f) {
    v <- suppressWarnings(try(as.numeric(f@R2), silent = TRUE))
    if (inherits(v, "try-error") || anyNA(v) || length(v) == 0) {
      # fallback: approximate from scores variance
      S  <- as.matrix(f@scores)
      s2 <- apply(S, 2, stats::var, na.rm = TRUE)
      s2 / sum(s2)
    } else v
  }
  var_frac <- get_var(fit)
  var_cum  <- cumsum(var_frac)
  
  # auto k: first PC count reaching target, clamped to [min_k, max_k] and available PCs
  k_auto <- which(var_cum >= target)[1]
  if (is.na(k_auto)) k_auto <- length(var_cum)
  k_auto <- max(min_k, min(k_auto, max_k, ncol(fit@scores)))
  
  list(fit = fit, var_frac = var_frac, var_cum = var_cum, k = k_auto)
}

# ---- PCA + auto-k for each omics ----
metab_pca <- pcak_choose_k(metab_M, pca_target_metab, pca_min_k, pca_max_k)
prot_pca  <- pcak_choose_k(prot_M,  pca_target_prot,  pca_min_k, pca_max_k)

scores_m <- as.matrix(metab_pca$fit@scores)[, 1:metab_pca$k, drop = FALSE]
scores_p <- as.matrix(prot_pca$fit@scores )[, 1:prot_pca$k,  drop = FALSE]

cumvar_m <- metab_pca$var_cum[metab_pca$k]
cumvar_p <- prot_pca$var_cum[prot_pca$k]

# ---- Mahalanobis distances on chosen scores ----
md_fun <- function(S) {
  m <- colMeans(S)
  V <- stats::cov(S)
  stats::mahalanobis(S, m, V)
}
md_m <- md_fun(scores_m)
md_p <- md_fun(scores_p)

cut_m <- stats::qchisq(chi_q, df = ncol(scores_m))
cut_p <- stats::qchisq(chi_q, df = ncol(scores_p))

flag_m <- md_m > cut_m
flag_p <- md_p > cut_p
out_idx <- which(flag_m | flag_p)
keep_idx <- setdiff(seq_len(nrow(metabolome)), out_idx)

# ---- Filter data ----
metabolome_clean <- metabolome[keep_idx, , drop = FALSE]
proteome_clean   <- proteome  [keep_idx, , drop = FALSE]
covars_clean     <- covars    [keep_idx, , drop = FALSE]
stopifnot(all(metabolome_clean$sid == proteome_clean$sid))

# ---- Write per-sample MD table ----
md_tbl <- tibble(
  sid             = covars$sid,
  md_metabolome   = as.numeric(md_m),
  md_proteome     = as.numeric(md_p),
  cut_metabolome  = cut_m,
  cut_proteome    = cut_p,
  flag_metabolome = as.logical(flag_m),
  flag_proteome   = as.logical(flag_p),
  outlier         = as.logical(flag_m | flag_p)
)
readr::write_csv(md_tbl, file.path(res_dir, "block3_mahalanobis_per_sample.csv"))

# ---- Summary CSV ----
summ <- tibble(
  n_before            = n_before,
  n_removed           = length(out_idx),
  n_after             = nrow(metabolome_clean),
  k_pcs_metabolome    = metab_pca$k,
  k_pcs_proteome      = prot_pca$k,
  cumvar_metabolome   = round(cumvar_m, 3),
  cumvar_proteome     = round(cumvar_p, 3),
  chi2_quantile       = chi_q,
  chi2_cut_metabolome = cut_m,
  chi2_cut_proteome   = cut_p
)
readr::write_csv(summ, file.path(res_dir, "block3_outlier_summary.csv"))

# ---- Scree plots (bars = per-PC %, line = cumulative), with vertical line at k ----
plot_scree <- function(var_frac, var_cum, k, title, out_png) {
  df <- tibble(PC = seq_along(var_frac),
               frac = var_frac,
               cum  = var_cum)
  g <- ggplot(df %>% slice(1:pca_max_k), aes(PC, frac)) +
    geom_col() +
    geom_point(aes(y = cum)) +
    geom_line(aes(y = cum)) +
    geom_vline(xintercept = k, linetype = "dashed") +
    labs(title = title, x = "PC", y = "% variance / cumulative") +
    theme_minimal()
  ggsave(out_png, g, width = 7, height = 5, dpi = 300)
}

plot_scree(metab_pca$var_frac, metab_pca$var_cum, metab_pca$k,
           "Metabolome PCA (NIPALS): % variance (bars) & cumulative (line)",
           file.path(fig_dir, "block3_scree_metabolome.png"))

plot_scree(prot_pca$var_frac, prot_pca$var_cum, prot_pca$k,
           "Proteome PCA (NIPALS): % variance (bars) & cumulative (line)",
           file.path(fig_dir, "block3_scree_proteome.png"))

# ---- Save RDS used by Block 4 ----
saveRDS(list(
  metabolome_clean      = metabolome_clean,
  proteome_clean        = proteome_clean,
  covars_clean          = covars_clean,
  pca_k_metabolome      = metab_pca$k,
  pca_k_proteome        = prot_pca$k,
  pca_varfrac_metabolome= metab_pca$var_frac,
  pca_varfrac_proteome  = prot_pca$var_frac,
  mahalanobis           = md_tbl,
  removed_idx           = out_idx,
  kept_idx              = keep_idx
), file = file.path(proc_dir, "block3_scaled_pca_outliers.rds"))

cat(sprintf(
  "Block 3 complete. Removed outliers: %d | k_m=%d (cum=%.1f%%), k_p=%d (cum=%.1f%%)\n",
  length(out_idx), metab_pca$k, 100*cumvar_m, prot_pca$k, 100*cumvar_p))
# -------------------------------------------------------------------------------
