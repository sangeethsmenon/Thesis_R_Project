packages <- c("tidyverse","here","splines")
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

b3_path <- file.path(config$paths$data_processed, "block3_scaled_pca_outliers.rds")
if (!file.exists(b3_path)) stop("Block 3 outputs not found: ", b3_path)
b3 <- readRDS(b3_path)

id_cols <- c("sid","sex","age_in0")
covars <- b3$covars_clean
covars$sex <- factor(covars$sex)
sex_levels <- levels(covars$sex)

X_full   <- model.matrix(~ sex + splines::ns(age_in0, df = 3), data = covars)
X_nosex  <- model.matrix(~ splines::ns(age_in0, df = 3), data = covars)
X_noage  <- model.matrix(~ sex, data = covars)

XtX_f  <- crossprod(X_full);  XtX_f_inv  <- solve(XtX_f)
XtX_ns <- crossprod(X_nosex); XtX_ns_inv <- solve(XtX_ns)
XtX_na <- crossprod(X_noage); XtX_na_inv <- solve(XtX_na)

n <- nrow(X_full)
k_full  <- ncol(X_full)
df_sex  <- ncol(X_full) - ncol(X_nosex)
df_age  <- ncol(X_full) - ncol(X_noage)

analyze <- function(df_num, dataset_label){
  Y <- as.matrix(df_num)
  beta_f  <- XtX_f_inv  %*% crossprod(X_full,  Y)
  e_full  <- Y - X_full  %*% beta_f
  SSE_f   <- colSums(e_full^2)
  tss     <- colSums((Y - matrix(colMeans(Y), n, ncol(Y), byrow = TRUE))^2)
  R2_all  <- pmax(0, 1 - SSE_f / pmax(1e-12, tss))
  
  beta_ns <- XtX_ns_inv %*% crossprod(X_nosex, Y)
  e_ns    <- Y - X_nosex %*% beta_ns
  SSE_ns  <- colSums(e_ns^2)
  
  beta_na <- XtX_na_inv %*% crossprod(X_noage, Y)
  e_na    <- Y - X_noage %*% beta_na
  SSE_na  <- colSums(e_na^2)
  
  F_sex <- ((SSE_ns - SSE_f) / df_sex) / (SSE_f / pmax(1, (n - k_full)))
  F_age <- ((SSE_na - SSE_f) / df_age) / (SSE_f / pmax(1, (n - k_full)))
  p_sex <- pf(F_sex, df_sex, n - k_full, lower.tail = FALSE)
  p_age <- pf(F_age, df_age, n - k_full, lower.tail = FALSE)
  
  partR2_sex <- pmax(0, 1 - SSE_f / pmax(1e-12, SSE_ns))
  partR2_age <- pmax(0, 1 - SSE_f / pmax(1e-12, SSE_na))
  
  sex_cols <- grepl("^sex", colnames(X_full))
  beta_sex <- rep(NA_real_, ncol(Y))
  if (sum(sex_cols) == 1) beta_sex <- as.numeric(beta_f[sex_cols, , drop = FALSE])
  
  tibble(
    dataset = dataset_label,
    feature = colnames(df_num),
    R2_all = R2_all,
    partR2_sex = partR2_sex,
    partR2_age = partR2_age,
    p_sex = p_sex,
    p_age = p_age,
    q_sex = p.adjust(p_sex, "BH"),
    q_age = p.adjust(p_age, "BH"),
    beta_sex = beta_sex
  )
}

metab_num <- b3$metabolome_clean %>% select(-all_of(id_cols))
prot_num  <- b3$proteome_clean   %>% select(-all_of(id_cols))

stats_m <- analyze(metab_num, "metabolome")
stats_p <- analyze(prot_num,  "proteome")
stats_all <- bind_rows(stats_m, stats_p)

readr::write_csv(stats_all, file.path(config$paths$results, "block4_diagnostics_partialR2.csv"))

top20_sex <- stats_all %>% group_by(dataset) %>% arrange(desc(partR2_sex), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()
top20_age <- stats_all %>% group_by(dataset) %>% arrange(desc(partR2_age), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()

readr::write_csv(top20_sex, file.path(config$paths$results, "block4_top20_by_partialR2_sex.csv"))
readr::write_csv(top20_age, file.path(config$paths$results, "block4_top20_by_partialR2_age.csv"))

dir.create(config$paths$figures, showWarnings = FALSE, recursive = TRUE)
g_m <- ggplot(stats_m, aes(R2_all)) + geom_histogram(bins = 30) + labs(title = "Metabolome: R² (age+sex) distribution", x = "R²", y = "Features") + theme_minimal()
g_p <- ggplot(stats_p, aes(R2_all)) + geom_histogram(bins = 30) + labs(title = "Proteome: R² (age+sex) distribution", x = "R²", y = "Features") + theme_minimal()
ggsave(file.path(config$paths$figures, "block4_R2_hist_metabolome.png"), g_m, width = 7, height = 5, dpi = 300)
ggsave(file.path(config$paths$figures, "block4_R2_hist_proteome.png"),  g_p, width = 7, height = 5, dpi = 300)

top_dir <- file.path(config$paths$results, "block4_top20_tables_readme.txt")
txt <- c(
  paste0("Sex factor levels (reference = first level): ", paste(sex_levels, collapse = ", ")),
  "Interpretation: beta_sex > 0 means mean(Male) > mean(reference level).",
  "Files:",
  "  - block4_top20_by_partialR2_sex.csv",
  "  - block4_top20_by_partialR2_age.csv",
  "  - block4_diagnostics_partialR2.csv",
  "  - figures/block4_R2_hist_[metabolome|proteome].png"
)
writeLines(txt, top_dir)

cat("Diagnostics complete.\n")
cat("Sex levels: ", paste(sex_levels, collapse = ", "), "\n")
