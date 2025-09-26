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
metab_df <- b3$metabolome_clean
prot_df  <- b3$proteome_clean
covars   <- b3$covars_clean
covars$sex <- factor(covars$sex)

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

scale_cols <- function(M){M <- as.matrix(M); mu <- colMeans(M); sdv <- apply(M,2,sd); sdv[sdv==0] <- 1; sweep(sweep(M,2,mu,"-"),2,sdv,"/")}
std_cols   <- scale_cols

resid_and_stats <- function(Y){
  Y <- as.matrix(Y)
  beta_f  <- XtX_f_inv  %*% crossprod(X_full,  Y)
  e_full  <- Y - X_full  %*% beta_f
  SSE_f   <- colSums(e_full^2)
  sigma2f <- SSE_f / pmax(1, (n - k_full))
  tss     <- colSums((Y - matrix(colMeans(Y), n, ncol(Y), byrow = TRUE))^2)
  R2      <- pmax(0, 1 - SSE_f / pmax(1e-12, tss))
  
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
  
  sex_cols <- grepl("^sex", colnames(X_full))
  beta_sex <- rep(NA_real_, ncol(Y))
  if (sum(sex_cols) == 1) {
    beta_sex <- as.numeric(beta_f[sex_cols, , drop = FALSE])
  }
  
  list(resid = e_full,
       stats = tibble(beta_sex = beta_sex, p_sex = p_sex, p_age = p_age, R2 = R2))
}

metab_num <- metab_df %>% select(-all_of(id_cols))
prot_num  <- prot_df  %>% select(-all_of(id_cols))

cat("Residualizing metabolome...\n")
m_res <- resid_and_stats(metab_num)
cat("Residualizing proteome...\n")
p_res <- resid_and_stats(prot_num)

metab_resid_std <- std_cols(m_res$resid)
prot_resid_std  <- std_cols(p_res$resid)

metab_stats <- m_res$stats %>% mutate(feature = colnames(metab_num), dataset = "metabolome",
                                      q_sex = p.adjust(p_sex, "BH"), q_age = p.adjust(p_age, "BH")) %>%
  select(dataset, feature, beta_sex, p_sex, q_sex, p_age, q_age, R2)
prot_stats  <- p_res$stats %>% mutate(feature = colnames(prot_num),  dataset = "proteome",
                                      q_sex = p.adjust(p_sex, "BH"), q_age = p.adjust(p_age, "BH")) %>%
  select(dataset, feature, beta_sex, p_sex, q_sex, p_age, q_age, R2)

stats_all <- bind_rows(metab_stats, prot_stats)

readr::write_csv(stats_all, file.path(config$paths$results, "block4_covariate_effects_per_feature.csv"))

summ <- stats_all %>%
  group_by(dataset) %>%
  summarise(
    n = n(),
    n_q_sex = sum(q_sex < 0.05, na.rm = TRUE),
    n_q_age = sum(q_age < 0.05, na.rm = TRUE),
    med_R2 = median(R2, na.rm = TRUE),
    p90_R2 = quantile(R2, .90, na.rm = TRUE),
    p95_R2 = quantile(R2, .95, na.rm = TRUE),
    max_R2 = max(R2, na.rm = TRUE)
  )
readr::write_csv(summ, file.path(config$paths$results, "block4_residualization_summary.csv"))

saveRDS(list(
  covars_clean = covars,
  X_full = X_full,
  metabolome_resid_std = metab_resid_std,
  proteome_resid_std   = prot_resid_std,
  metabolome_features = colnames(metab_num),
  proteome_features   = colnames(prot_num),
  covariate_stats = stats_all
), file = file.path(config$paths$data_processed, "block4_residuals_objects.rds"))

cat("Block 4 complete.\n")
