packages <- c("tidyverse","pheatmap","splines")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

stopifnot(exists("config"), exists("metabolome"), exists("proteome"))

id_cols <- c("sid","sex","age_in0")
metab_mat  <- metabolome %>% select(-all_of(id_cols))
proteo_mat <- proteome   %>% select(-all_of(id_cols))

miss_feature <- function(df) {
  tibble(feature = colnames(df),
         n_missing = colSums(is.na(df)),
         n = nrow(df),
         pct_missing = 100 * n_missing / n)
}
miss_sample <- function(df) {
  tibble(row = seq_len(nrow(df)),
         n_missing = rowSums(is.na(df)),
         p = ncol(df),
         pct_missing = 100 * n_missing / p)
}

metab_feat_miss  <- miss_feature(metab_mat)  %>% mutate(dataset = "metabolome")
proteo_feat_miss <- miss_feature(proteo_mat) %>% mutate(dataset = "proteome")
proteo_row_miss  <- miss_sample(proteo_mat)  %>% mutate(dataset = "proteome")

write.csv(metab_feat_miss,  file.path(config$paths$results, "metabolome_feature_missingness.csv"), row.names = FALSE)
write.csv(proteo_feat_miss, file.path(config$paths$results, "proteome_feature_missingness.csv"),  row.names = FALSE)
write.csv(proteo_row_miss,  file.path(config$paths$results, "proteome_sample_missingness.csv"),   row.names = FALSE)

p_metab_hist <- ggplot(metab_feat_miss, aes(pct_missing)) +
  geom_histogram(bins = 30) + labs(title = "Metabolome: % missing per feature", x = "% missing", y = "Features") +
  theme_minimal()
p_prot_hist <- ggplot(proteo_feat_miss, aes(pct_missing)) +
  geom_histogram(bins = 30) + labs(title = "Proteome: % missing per feature", x = "% missing", y = "Features") +
  theme_minimal()

ggsave(file.path(config$paths$figures, "metabolome_missing_feature_hist.png"), p_metab_hist, width = 7, height = 5, dpi = 300)
ggsave(file.path(config$paths$figures, "proteome_missing_feature_hist.png"),  p_prot_hist,  width = 7, height = 5, dpi = 300)

set.seed(config$seed)
n_rows <- min(400, nrow(metabolome))
n_metab_cols <- min(100, ncol(metab_mat))
n_prot_cols  <- min(200, ncol(proteo_mat))
rows_idx <- sample(seq_len(nrow(metabolome)), n_rows)
metab_cols <- sample(colnames(metab_mat), n_metab_cols)
prot_cols  <- sample(colnames(proteo_mat), n_prot_cols)

metab_miss_pct  <- is.na(metabolome[rows_idx, metab_cols, drop = FALSE]) * 100
proteo_miss_pct <- is.na(proteome  [rows_idx, prot_cols,  drop = FALSE]) * 100

bk <- seq(0, 100, by = 10)
pheatmap(metab_miss_pct,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         legend = TRUE, legend_breaks = c(0, 50, 100),
         legend_labels = c("0%","50%","100%"),
         breaks = bk, main = "Metabolome missingness (%), sampled",
         filename = file.path(config$paths$figures, "metabolome_missing_heatmap.png"),
         width = 9, height = 6)
pheatmap(proteo_miss_pct,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         legend = TRUE, legend_breaks = c(0, 50, 100),
         legend_labels = c("0%","50%","100%"),
         breaks = bk, main = "Proteome missingness (%), sampled",
         filename = file.path(config$paths$figures, "proteome_missing_heatmap.png"),
         width = 10, height = 6)

assess_missing_assoc <- function(df_full, feat_mat) {
  if (!is.factor(df_full$sex)) df_full$sex <- factor(df_full$sex)
  age_spline <- splines::ns(df_full$age_in0, df = 3)
  res <- lapply(seq_len(ncol(feat_mat)), function(j){
    miss <- as.integer(is.na(feat_mat[, j]))
    if (sum(miss) == 0L || sum(miss) == length(miss)) return(NULL)
    fit <- try(glm(miss ~ df_full$sex + age_spline, family = binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    s <- summary(fit)$coefficients
    p_sex <- if ("df_full$sexMale" %in% rownames(s)) s["df_full$sexMale","Pr(>|z|)"] else NA_real_
    p_age <- max(NA_real_, suppressWarnings(anova(fit, test = "Chisq"))$`Pr(>Chi)`[2], na.rm = TRUE)
    tibble(feature = colnames(feat_mat)[j], p_sex = p_sex, p_age = p_age)
  })
  res <- bind_rows(res)
  if (nrow(res)) res <- res %>% mutate(q_sex = p.adjust(p_sex, "BH"),
                                       q_age = p.adjust(p_age, "BH"))
  res
}

metab_miss_assoc  <- assess_missing_assoc(metabolome, metab_mat)  %>% mutate(dataset = "metabolome")
proteo_miss_assoc <- assess_missing_assoc(proteome,   proteo_mat) %>% mutate(dataset = "proteome")

write.csv(metab_miss_assoc,  file.path(config$paths$results, "metabolome_missingness_vs_covariates.csv"), row.names = FALSE)
write.csv(proteo_miss_assoc, file.path(config$paths$results, "proteome_missingness_vs_covariates.csv"),  row.names = FALSE)

thr_metab_feat <- config$thresholds$metab_feat
thr_prot_feat  <- config$thresholds$prot_feat
thr_prot_row   <- config$thresholds$prot_row

metab_keep_feat_names  <- metab_feat_miss %>% filter(pct_missing <= 100*thr_metab_feat) %>% pull(feature)
proteo_keep_feat_names <- proteo_feat_miss %>% filter(pct_missing <= 100*thr_prot_feat)  %>% pull(feature)
proteo_keep_rows_idx   <- proteo_row_miss %>% filter(pct_missing <= 100*thr_prot_row)    %>% pull(row)

metabolome_filt <- bind_cols(metabolome[id_cols], metab_mat[, metab_keep_feat_names, drop = FALSE])
proteome_filt   <- bind_cols(proteome[proteo_keep_rows_idx, id_cols],
                             proteo_mat[proteo_keep_rows_idx, proteo_keep_feat_names, drop = FALSE])

retention <- tibble(
  dataset = c("metabolome","proteome_features","proteome_samples"),
  kept = c(ncol(metabolome_filt)-length(id_cols),
           length(proteo_keep_feat_names),
           length(proteo_keep_rows_idx)),
  total = c(ncol(metabolome)-length(id_cols),
            ncol(proteome)-length(id_cols),
            nrow(proteome))
)
write.csv(retention, file.path(config$paths$results, "block1_retention_counts.csv"), row.names = FALSE)

saveRDS(list(
  metabolome_filt = metabolome_filt,
  proteome_filt   = proteome_filt,
  metab_feat_miss = metab_feat_miss,
  proteo_feat_miss = proteo_feat_miss,
  proteo_row_miss  = proteo_row_miss,
  metab_miss_assoc = metab_miss_assoc,
  proteo_miss_assoc = proteo_miss_assoc
), file = file.path(config$paths$data_processed, "block1_filtered_objects.rds"))

cat("Kept metabolite features:", ncol(metabolome_filt)-length(id_cols), "\n")
cat("Kept protein features:", length(proteo_keep_feat_names), "\n")
cat("Kept proteome samples:", length(proteo_keep_rows_idx), "\n")
