packages <- c("tidyverse","here","softImpute","Matrix")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

stopifnot(exists("config"))
b1 <- readRDS(file.path(config$paths$data_processed, "block1_filtered_objects.rds"))

id_cols <- c("sid","sex","age_in0")
keep_sids <- intersect(b1$metabolome_filt$sid, b1$proteome_filt$sid)
metabolome_filt <- b1$metabolome_filt[match(keep_sids, b1$metabolome_filt$sid), ]
proteome_filt   <- b1$proteome_filt  [match(keep_sids, b1$proteome_filt$sid), ]
metabolome_filt$sex <- factor(metabolome_filt$sex)
proteome_filt$sex   <- factor(proteome_filt$sex)

metab_num <- metabolome_filt %>% select(-all_of(id_cols))
prot_num  <- proteome_filt   %>% select(-all_of(id_cols))

metab_na_before <- sum(is.na(metab_num))
prot_na_before  <- sum(is.na(prot_num))

set.seed(config$seed)
metab_imputed <- as_tibble(apply(metab_num, 2, function(x){x[is.na(x)] <- median(x, na.rm = TRUE); x}))

X <- as.matrix(prot_num)
mu <- colMeans(X, na.rm = TRUE)
sdv <- apply(X, 2, sd, na.rm = TRUE); sdv[sdv == 0] <- 1
Z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
lam <- softImpute::lambda0(Z) * 0.1
t0 <- proc.time()
fit <- softImpute::softImpute(Z, rank.max = 50, lambda = lam, type = "als", maxit = 200, trace.it = TRUE)
elapsed <- (proc.time() - t0)["elapsed"]
Zhat <- softImpute::complete(Z, fit)
Xhat <- sweep(sweep(Zhat, 2, sdv, "*"), 2, mu, "+")
prot_imputed <- as_tibble(Xhat)

metabolome_imp <- bind_cols(metabolome_filt %>% select(all_of(id_cols)), metab_imputed)
proteome_imp   <- bind_cols(proteome_filt   %>% select(all_of(id_cols)), prot_imputed)

metab_na_after <- sum(is.na(metabolome_imp %>% select(-all_of(id_cols))))
prot_na_after  <- sum(is.na(proteome_imp   %>% select(-all_of(id_cols))))
stopifnot(metab_na_after == 0, prot_na_after == 0)

n_m <- nrow(metab_imputed); p_m <- ncol(metab_imputed)
n_p <- nrow(prot_imputed);  p_p <- ncol(prot_imputed)

sum_tbl <- tibble(
  dataset = c("metabolome","proteome"),
  n_samples = c(n_m, n_p),
  n_features = c(p_m, p_p),
  n_cells = c(n_m*p_m, n_p*p_p),
  n_imputed_cells = c(metab_na_before, prot_na_before),
  pct_imputed_cells = 100 * c(metab_na_before/(n_m*p_m), prot_na_before/(n_p*p_p)),
  method = c("median-per-feature", paste0("softImpute ALS; rank.max=50; lambda=", signif(lam,3))),
  softImpute_elapsed_sec = c(NA_real_, as.numeric(elapsed))
)

write.csv(sum_tbl, file.path(config$paths$results, "block2_imputation_summary.csv"), row.names = FALSE)
saveRDS(list(
  metabolome_imp = metabolome_imp,
  proteome_imp   = proteome_imp,
  keep_sids      = keep_sids,
  imputation_summary = sum_tbl,
  softImpute_fit = fit
), file = file.path(config$paths$data_processed, "block2_imputed_objects.rds"))

print(sum_tbl)
cat("Imputation complete.\n")
