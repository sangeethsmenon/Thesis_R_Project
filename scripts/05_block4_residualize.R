# --- 05_block4_residualize.R (NA-aware residualization, no imputation) -------

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(jsonlite); library(splines)
})

# --------- Loader (config + QC TSVs + kept sample index from Block 3) --------
if (!exists("config")) {
  cfg_rds  <- here("results","config.rds")
  cfg_json <- here("results","config.json")
  if (file.exists(cfg_rds)) config <- readRDS(cfg_rds)
  else if (file.exists(cfg_json)) config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  else stop("Config not found. Run Block 0.")
}

prot_qc_path  <- file.path(config$paths$data_processed, "proteome_qc.tsv")
metab_qc_path <- file.path(config$paths$data_processed, "metabolome_qc.tsv")
b3_rds        <- file.path(config$paths$data_processed, "block3_scaled_pca_outliers.rds")
stopifnot(file.exists(prot_qc_path), file.exists(metab_qc_path), file.exists(b3_rds))

proteome   <- readr::read_tsv(prot_qc_path,  show_col_types = FALSE)
metabolome <- readr::read_tsv(metab_qc_path, show_col_types = FALSE)
b3 <- readRDS(b3_rds)

keep_idx <- b3$kept_idx
proteome   <- proteome  [keep_idx, , drop = FALSE]
metabolome <- metabolome[keep_idx, , drop = FALSE]
stopifnot(nrow(proteome) == nrow(metabolome), all(proteome$sid == metabolome$sid))

# ------------------------- Design and helpers --------------------------------
id_cols <- c("sid","sex","age_in0")
covars <- metabolome %>% select(all_of(id_cols))
covars$sex <- factor(covars$sex)

X_full  <- model.matrix(~ sex + splines::ns(age_in0, df = 3), data = covars)

# Residualize one column y against X_full on rows where y is observed
residualize_one <- function(y, X_full) {
  ok <- !is.na(y)
  if (sum(ok) <= ncol(X_full) + 1L) return(rep(NA_real_, length(y)))
  yj <- y[ok]
  Xf <- X_full[ok, , drop = TRUE]
  betaf <- solve(crossprod(Xf), crossprod(Xf, yj))
  ef <- yj - drop(Xf %*% betaf)
  out <- rep(NA_real_, length(y))
  out[ok] <- ef
  out
}

# Standardize a numeric vector while keeping NAs
z_std <- function(v) {
  mu <- mean(v, na.rm = TRUE); sdv <- sd(v, na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
  (v - mu) / sdv
}

# ------------------------- Residualize (NA-aware) -----------------------------
metab_num <- metabolome %>% select(-all_of(id_cols))
prot_num  <- proteome   %>% select(-all_of(id_cols))

cat("Residualizing metabolome (NA-aware)…\n")
metab_resid <- as_tibble(lapply(metab_num, residualize_one, X_full = X_full))
cat("Residualizing proteome (NA-aware)…\n")
prot_resid  <- as_tibble(lapply(prot_num,  residualize_one, X_full = X_full))

# Standardize columns (z-score) without touching NAs
metab_resid_std <- metab_resid %>% mutate(across(everything(), z_std))
prot_resid_std  <- prot_resid  %>% mutate(across(everything(), z_std))

# ------------------------------- Save ----------------------------------------
saveRDS(list(
  covars_clean = covars,
  X_full = X_full,
  metabolome_resid_std = as.matrix(metab_resid_std),
  proteome_resid_std   = as.matrix(prot_resid_std),
  metabolome_features  = colnames(metab_num),
  proteome_features    = colnames(prot_num)
), file = file.path(config$paths$data_processed, "block4_residuals_objects.rds"))

# Light summary (how many non-NA residuals per feature)
summ <- bind_rows(
  tibble(dataset="metabolome", feature=names(metab_resid_std),
         n_obs = colSums(!is.na(metab_resid_std))),
  tibble(dataset="proteome",   feature=names(prot_resid_std),
         n_obs = colSums(!is.na(prot_resid_std)))
)
readr::write_csv(summ, file.path(config$paths$results, "block4_residuals_nobs.csv"))

cat("Block 4 residualization complete (NA-aware).\n")
