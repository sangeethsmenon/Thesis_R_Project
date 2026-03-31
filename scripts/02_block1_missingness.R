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


# ---- Heatmaps: missingness indicator (NOT percent) ----
# Each tile represents one measurement (sample × feature):
# 0 = observed value present, 1 = missing (NA)

set.seed(config$seed)
n_rows <- min(400, nrow(metabolome))
n_metab_cols <- min(100, ncol(metab_mat))
n_prot_cols  <- min(200, ncol(proteo_mat))

rows_idx  <- sample(seq_len(nrow(metabolome)), n_rows)
metab_cols <- sample(colnames(metab_mat), n_metab_cols)
prot_cols  <- sample(colnames(proteo_mat), n_prot_cols)

# Binary indicator matrices
metab_miss_ind  <- is.na(metabolome[rows_idx, metab_cols, drop = FALSE]) * 1L
proteo_miss_ind <- is.na(proteome  [rows_idx, prot_cols,  drop = FALSE]) * 1L

# Legend + colors for binary matrix
# breaks must bracket 0 and 1 correctly
bk_bin <- c(-0.5, 0.5, 1.5)
col_bin <- c("white", "#2C7FB8")  # observed, missing

# Metabolome heatmap
pheatmap::pheatmap(
  metab_miss_ind,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = FALSE,
  legend = TRUE,
  breaks = bk_bin,
  color = col_bin,
  legend_breaks = c(0, 1),
  legend_labels = c("Observed", "Missing"),
  main = "Metabolome missingness (sampled)",
  fontsize = 10,
  filename = file.path(config$paths$figures, "metabolome_missing_heatmap.png"),
  width = 9, height = 6
)

# Proteome heatmap
pheatmap::pheatmap(
  proteo_miss_ind,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = FALSE,
  legend = TRUE,
  breaks = bk_bin,
  color = col_bin,
  legend_breaks = c(0, 1),
  legend_labels = c("Observed", "Missing"),
  main = "Proteome missingness (sampled)",
  fontsize = 10,
  filename = file.path(config$paths$figures, "proteome_missing_heatmap.png"),
  width = 10, height = 6
)

cat("Heatmaps saved with binary legend (Observed/Missing).\n")
cat("Rows = sampled individuals (samples), columns = sampled features.\n")
