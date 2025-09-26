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

b2_path <- file.path(config$paths$data_processed, "block2_imputed_objects.rds")
if (!file.exists(b2_path)) stop("Block 2 outputs not found: ", b2_path)
b2 <- readRDS(b2_path)

id_cols <- c("sid","sex","age_in0")
metabolome_imp <- b2$metabolome_imp
proteome_imp   <- b2$proteome_imp
stopifnot(all(metabolome_imp$sid == proteome_imp$sid))

covars <- metabolome_imp %>% select(all_of(id_cols))

scale_mat <- function(m){
  m <- as.matrix(m)
  mu <- colMeans(m)
  sdv <- apply(m, 2, sd)
  sdv[sdv == 0] <- 1
  sweep(sweep(m, 2, mu, "-"), 2, sdv, "/")
}

metab_scaled <- scale_mat(metabolome_imp %>% select(-all_of(id_cols)))
prot_scaled  <- scale_mat(proteome_imp   %>% select(-all_of(id_cols)))

pca_metab <- prcomp(metab_scaled, center = FALSE, scale. = FALSE)
pca_prot  <- prcomp(prot_scaled,  center = FALSE, scale. = FALSE)

k_m <- min(5, ncol(pca_metab$x))
k_p <- min(5, ncol(pca_prot$x))
pcm <- pca_metab$x[, 1:k_m, drop = FALSE]
pcp <- pca_prot$x[,  1:k_p, drop = FALSE]

mah <- function(X){m <- colMeans(X); S <- cov(X); mahalanobis(X, m, S)}
md_m <- mah(pcm)
md_p <- mah(pcp)

cut_m <- qchisq(0.999, df = k_m)
cut_p <- qchisq(0.999, df = k_p)
flag_m <- md_m > cut_m
flag_p <- md_p > cut_p
out_idx <- which(flag_m | flag_p)
keep_idx <- setdiff(seq_len(nrow(metabolome_imp)), out_idx)

metabolome_clean <- metabolome_imp[keep_idx, , drop = FALSE]
proteome_clean   <- proteome_imp  [keep_idx, , drop = FALSE]
covars_clean     <- covars        [keep_idx, , drop = FALSE]
stopifnot(all(metabolome_clean$sid == proteome_clean$sid))

metab_scaled_clean <- scale_mat(metabolome_clean %>% select(-all_of(id_cols)))
prot_scaled_clean  <- scale_mat(proteome_clean   %>% select(-all_of(id_cols)))

var_m <- (pca_metab$sdev^2) / sum(pca_metab$sdev^2)
var_p <- (pca_prot$sdev^2)  / sum(pca_prot$sdev^2)
var_tbl_m <- tibble(PC = seq_along(var_m), pct = 100*var_m)
var_tbl_p <- tibble(PC = seq_along(var_p), pct = 100*var_p)

out_tbl <- tibble(
  sid = covars$sid,
  md_metabolome = as.numeric(md_m),
  md_proteome   = as.numeric(md_p),
  cut_metabolome = cut_m,
  cut_proteome   = cut_p,
  flag_metabolome = as.logical(flag_m),
  flag_proteome   = as.logical(flag_p),
  outlier = as.logical(flag_m | flag_p)
)

readr::write_csv(out_tbl, file.path(config$paths$results, "block3_mahalanobis_per_sample.csv"))

summ <- tibble(
  n_before = nrow(metabolome_imp),
  n_removed = length(out_idx),
  n_after = nrow(metabolome_clean),
  k_pcs_metabolome = k_m,
  k_pcs_proteome   = k_p,
  chi2_cut_metabolome = cut_m,
  chi2_cut_proteome   = cut_p
)
readr::write_csv(summ, file.path(config$paths$results, "block3_outlier_summary.csv"))

dir.create(config$paths$figures, showWarnings = FALSE, recursive = TRUE)
p1 <- ggplot(var_tbl_m %>% slice(1:30), aes(PC, pct)) + geom_col() + labs(title = "Metabolome PCA: % variance by PC (top 30)", x = "PC", y = "% variance") + theme_minimal()
p2 <- ggplot(var_tbl_p %>% slice(1:30), aes(PC, pct)) + geom_col() + labs(title = "Proteome PCA: % variance by PC (top 30)", x = "PC", y = "% variance") + theme_minimal()
ggsave(file.path(config$paths$figures, "block3_scree_metabolome.png"), p1, width = 7, height = 5, dpi = 300)
ggsave(file.path(config$paths$figures, "block3_scree_proteome.png"),  p2, width = 7, height = 5, dpi = 300)

saveRDS(list(
  metabolome_clean = metabolome_clean,
  proteome_clean   = proteome_clean,
  covars_clean     = covars_clean,
  metab_scaled     = metab_scaled_clean,
  prot_scaled      = prot_scaled_clean,
  pca_metabolome   = list(sdev = pca_metab$sdev, rotation = pca_metab$rotation),
  pca_proteome     = list(sdev = pca_prot$sdev, rotation = pca_prot$rotation),
  mahalanobis      = out_tbl,
  removed_idx      = out_idx,
  kept_idx         = keep_idx
), file = file.path(config$paths$data_processed, "block3_scaled_pca_outliers.rds"))

cat("Block 3 complete. Removed outliers:", length(out_idx), "\n")
