## ===============================
## RAW DATA MISSINGNESS ASSESSMENT
## ===============================

packages <- c("tidyverse", "pheatmap")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

stopifnot(exists("metabolome"), exists("proteome"))

## -------------------------------
## Identify ID columns explicitly
## -------------------------------
id_cols <- c("sid", "sex", "age_in0")

# Ensure base data.frame compatibility
metabolome <- as.data.frame(metabolome)
proteome   <- as.data.frame(proteome)

# Remove ID columns safely
metab_mat  <- metabolome[, !(colnames(metabolome) %in% id_cols), drop = FALSE]
proteo_mat <- proteome[,   !(colnames(proteome)   %in% id_cols), drop = FALSE]

## -------------------------------
## Missingness summaries
## -------------------------------
miss_feature <- function(mat) {
  data.frame(
    feature      = colnames(mat),
    n_missing    = colSums(is.na(mat)),
    n_samples    = nrow(mat),
    pct_missing  = 100 * colSums(is.na(mat)) / nrow(mat)
  )
}

miss_sample <- function(mat) {
  data.frame(
    sample       = seq_len(nrow(mat)),
    n_missing    = rowSums(is.na(mat)),
    n_features   = ncol(mat),
    pct_missing  = 100 * rowSums(is.na(mat)) / ncol(mat)
  )
}

metab_feat_miss  <- miss_feature(metab_mat)
proteo_feat_miss <- miss_feature(proteo_mat)
proteo_row_miss  <- miss_sample(proteo_mat)

## -------------------------------
## Feature-level histograms
## -------------------------------
ggsave(
  "metabolome_missing_feature_hist.png",
  ggplot(metab_feat_miss, aes(pct_missing)) +
    geom_histogram(bins = 30, fill = "steelblue") +
    labs(
      title = "Metabolome: feature-level missingness (raw data)",
      x = "% missing across samples",
      y = "Number of features"
    ) +
    theme_minimal(),
  width = 7, height = 5, dpi = 300
)

ggsave(
  "proteome_missing_feature_hist.png",
  ggplot(proteo_feat_miss, aes(pct_missing)) +
    geom_histogram(bins = 30, fill = "steelblue") +
    labs(
      title = "Proteome: feature-level missingness (raw data)",
      x = "% missing across samples",
      y = "Number of features"
    ) +
    theme_minimal(),
  width = 7, height = 5, dpi = 300
)

## -------------------------------
## Heatmaps (BINARY missingness)
## -------------------------------
# Each tile = one sample × feature
# 0 = observed, 1 = missing

set.seed(123)

n_rows <- min(400, nrow(metabolome))
n_metab_cols <- min(100, ncol(metab_mat))
n_prot_cols  <- min(200, ncol(proteo_mat))

rows_idx  <- sample(seq_len(nrow(metabolome)), n_rows)
metab_cols <- sample(colnames(metab_mat), n_metab_cols)
prot_cols  <- sample(colnames(proteo_mat), n_prot_cols)

metab_miss_bin  <- is.na(metab_mat [rows_idx, metab_cols, drop = FALSE]) * 1L
proteo_miss_bin <- is.na(proteo_mat[rows_idx, prot_cols,  drop = FALSE]) * 1L

# Color mapping (binary)
col_bin <- c("#2166AC", "white")
bk_bin  <- c(-0.5, 0.5, 1.5)


pheatmap(
  metab_miss_bin,
  color = col_bin,
  breaks = bk_bin,
  legend = TRUE,
  legend_breaks = c(0, 1),
  legend_labels = c("Observed", "Missing"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Metabolome missingness (raw data, sampled)\nRows = samples, columns = metabolites",
  fontsize = 10,
  border_color = NA,
  filename = "metabolomes_missing_heatmap.png",
  width = 9,
  height = 6
)

pheatmap(
  proteo_miss_bin,
  color = col_bin,
  breaks = bk_bin,
  legend = TRUE,
  legend_breaks = c(0, 1),
  legend_labels = c("Observed", "Missing"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Proteome missingness (raw data, sampled)\nRows = samples, columns = proteins",
  fontsize = 10,
  border_color = NA,
  filename = "proteome_missing_heatmap.png",
  width = 10,
  height = 6
)
cat("✔ Raw-data missingness heatmaps generated.\n")
cat("Rows = samples, columns = features.\n")
cat("White = observed, blue = missing.\n")
