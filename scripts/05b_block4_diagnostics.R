# --- 05b_block4_diagnostics.R (NA-aware, no imputation) ----------------------

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(jsonlite); library(splines)
})

# ---------------- Standard loader (config + QC TSVs + keep_idx) --------------
if (!exists("config")) {
  cfg_rds  <- here("results","config.rds")
  cfg_json <- here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else {
    stop("Config not found. Run Block 0.")
  }
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

# ----------------------------- Setup -----------------------------------------
id_cols <- c("sid","sex","age_in0")
covars <- metabolome %>% select(all_of(id_cols))
covars$sex <- factor(covars$sex)
sex_levels <- levels(covars$sex)

# Model matrices (on the full covariate frame; we will subset rows per feature)
X_full  <- model.matrix(~ sex + splines::ns(age_in0, df = 3), data = covars)
X_nosex <- model.matrix(~ splines::ns(age_in0, df = 3), data = covars)
X_noage <- model.matrix(~ sex, data = covars)

k_full  <- ncol(X_full)
k_nosex <- ncol(X_nosex)
k_noage <- ncol(X_noage)

# --------------- Per-feature NA-aware diagnostics helper ---------------------
diag_one <- function(y, X_full, X_nosex, X_noage) {
  ok <- !is.na(y)
  n_j <- sum(ok)
  
  if (n_j <= max(k_full, k_nosex, k_noage) + 1L) {
    return(tibble(
      n = n_j, R2_all = NA_real_, partR2_sex = NA_real_, partR2_age = NA_real_,
      p_sex = NA_real_, p_age = NA_real_, beta_sex = NA_real_
    ))
  }
  
  yj      <- y[ok]
  Xf      <- X_full [ok, , drop = FALSE]
  Xns     <- X_nosex[ok, , drop = FALSE]
  Xna     <- X_noage[ok, , drop = FALSE]
  
  # Full model
  betaf   <- solve(crossprod(Xf), crossprod(Xf, yj))
  ef      <- yj - drop(Xf %*% betaf)
  SSEf    <- sum(ef^2)
  
  # Reduced models
  betans  <- solve(crossprod(Xns), crossprod(Xns, yj))
  ens     <- yj - drop(Xns %*% betans)
  SSEns   <- sum(ens^2)
  
  betana  <- solve(crossprod(Xna), crossprod(Xna, yj))
  ena     <- yj - drop(Xna %*% betana)
  SSEna   <- sum(ena^2)
  
  # R2 and partial R2 (all computed on the same subset of rows)
  ybar    <- mean(yj)
  TSS     <- sum((yj - ybar)^2)
  R2_all  <- max(0, 1 - SSEf / max(1e-12, TSS))
  partR2_sex <- max(0, 1 - SSEf / max(1e-12, SSEns))
  partR2_age <- max(0, 1 - SSEf / max(1e-12, SSEna))
  
  # F-tests
  df1_sex <- k_full  - k_nosex
  df1_age <- k_full  - k_noage
  df2     <- n_j - k_full
  if (df2 <= 0) {
    p_sex <- NA_real_; p_age <- NA_real_
  } else {
    F_sex <- ((SSEns - SSEf) / df1_sex) / (SSEf / df2)
    F_age <- ((SSEna - SSEf) / df1_age) / (SSEf / df2)
    p_sex <- pf(F_sex, df1_sex, df2, lower.tail = FALSE)
    p_age <- pf(F_age, df1_age, df2, lower.tail = FALSE)
  }
  
  # Beta for sex (if a single coefficient — i.e., two-level sex with reference)
  sex_cols <- grepl("^sex", colnames(Xf))
  beta_sex <- if (sum(sex_cols) == 1) betaf[sex_cols] else NA_real_
  
  tibble(
    n = n_j,
    R2_all = R2_all,
    partR2_sex = partR2_sex,
    partR2_age = partR2_age,
    p_sex = p_sex,
    p_age = p_age,
    beta_sex = as.numeric(beta_sex)
  )
}

analyze_df <- function(df_num, dataset_label) {
  out <- lapply(seq_len(ncol(df_num)), function(j) {
    y <- df_num[[j]]
    res <- diag_one(y, X_full, X_nosex, X_noage)
    res$feature <- colnames(df_num)[j]
    res
  })
  bind_rows(out) %>%
    mutate(dataset = dataset_label,
           q_sex = p.adjust(p_sex, "BH"),
           q_age = p.adjust(p_age, "BH")) %>%
    relocate(dataset, feature, n)
}

# ------------------------- Run diagnostics -----------------------------------
metab_num <- metabolome %>% select(-all_of(id_cols))
prot_num  <- proteome   %>% select(-all_of(id_cols))

cat("Running diagnostics (NA-aware)…\n")
stats_m <- analyze_df(metab_num, "metabolome")
stats_p <- analyze_df(prot_num,  "proteome")
stats_all <- bind_rows(stats_m, stats_p)

# ----------------------------- Save outputs ----------------------------------
readr::write_csv(stats_all, file.path(config$paths$results, "block4_diagnostics_partialR2.csv"))

top20_sex <- stats_all %>%
  group_by(dataset) %>% arrange(desc(partR2_sex), .by_group = TRUE) %>%
  slice_head(n = 20) %>% ungroup()
top20_age <- stats_all %>%
  group_by(dataset) %>% arrange(desc(partR2_age), .by_group = TRUE) %>%
  slice_head(n = 20) %>% ungroup()

readr::write_csv(top20_sex, file.path(config$paths$results, "block4_top20_by_partialR2_sex.csv"))
readr::write_csv(top20_age, file.path(config$paths$results, "block4_top20_by_partialR2_age.csv"))

dir.create(config$paths$figures, showWarnings = FALSE, recursive = TRUE)
g_m <- ggplot(stats_m, aes(R2_all)) + geom_histogram(bins = 30) +
  labs(title = "Metabolome: R² (age+sex) distribution", x = "R²", y = "Features") +
  theme_minimal()
g_p <- ggplot(stats_p, aes(R2_all)) + geom_histogram(bins = 30) +
  labs(title = "Proteome: R² (age+sex) distribution", x = "R²", y = "Features") +
  theme_minimal()
ggsave(file.path(config$paths$figures, "block4_R2_hist_metabolome.png"), g_m, width = 7, height = 5, dpi = 300)
ggsave(file.path(config$paths$figures, "block4_R2_hist_proteome.png"),  g_p, width = 7, height = 5, dpi = 300)

# Small README to interpret tables
txt <- c(
  paste0("Sex factor levels (reference = first level): ", paste(sex_levels, collapse = ", ")),
  "Interpretation: beta_sex > 0 means mean(level!=ref) > mean(ref) after adjusting for age spline.",
  "Files:",
  "  - results/block4_top20_by_partialR2_sex.csv",
  "  - results/block4_top20_by_partialR2_age.csv",
  "  - results/block4_diagnostics_partialR2.csv",
  "  - figures/block4_R2_hist_[metabolome|proteome].png"
)
writeLines(txt, file.path(config$paths$results, "block4_top20_tables_readme.txt"))

cat("Diagnostics complete.\n")
cat("Sex levels: ", paste(sex_levels, collapse = ", "), "\n")
