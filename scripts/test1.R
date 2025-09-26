# For metabolome
metabolome <- read.table("~/732A64 THESIS/Data_Sangeeth_20250829/NMR_metabolome_n249_SN13411_20250829.txt", header = TRUE, sep = "\t")
# For proteome
proteome <- read.table("~/732A64 THESIS/Data_Sangeeth_20250829/NMR_proteome_n2923_SN13411_20250829.txt", header = TRUE, sep = "\t")

# --- Step 1: Basic structure ---
# Dimensions
dim(metabolome)   # Rows x Cols
dim(proteome)

# Column names (first 10 for quick check)
head(colnames(metabolome), 10)
head(colnames(proteome), 10)

# Data types
str(metabolome[, 1:10])  
str(proteome[, 1:10])  

# --- Step 2: Sample alignment ---
# Check if sample IDs match in both datasets
all(metabolome$sid == proteome$sid)   # TRUE if perfect match
length(intersect(metabolome$sid, proteome$sid))  # Number of common samples

# --- Step 3: Summary statistics ---
# Age distribution
summary(metabolome$age_in0)

# Sex distribution
table(metabolome$sex)

# Missing values
colSums(is.na(metabolome))[1:10]  # first 10 columns
colSums(is.na(proteome))[1:10]

# Proportion missing per dataset
mean(is.na(metabolome))
mean(is.na(proteome))


# --- Setup ---
packages <- c("tidyverse", "naniar", "pheatmap")
need <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
library(tidyverse)
library(naniar)
library(pheatmap)

# If you haven't already:
# metabolome <- read.table(".../NMR_metabolome_n249_SN13411_20250829.txt", header=TRUE, sep="\t")
# proteome   <- read.table(".../NMR_proteome_n2923_SN13411_20250829.txt",  header=TRUE, sep="\t")

# --- 1) Age distribution (metabolome; same ages as proteome) ---
ggplot(metabolome, aes(x = age_in0)) +
  geom_histogram(bins = 30) +
  labs(title = "Age distribution", x = "Age at baseline", y = "Count")

# --- 2) Sex distribution ---
ggplot(metabolome, aes(x = sex)) +
  geom_bar() +
  labs(title = "Sex distribution", x = "Sex", y = "Count")

# --- 3) Missingness overview: % missing per feature (quick scan) ---
miss_pct <- function(df) {
  tibble(
    feature = names(df),
    missing = colSums(is.na(df)),
    total   = nrow(df),
    pct     = 100 * missing / total
  ) %>% arrange(desc(pct))
}

metab_miss <- miss_pct(metabolome) %>% filter(!feature %in% c("sid","sex","age_in0"))
proteo_miss <- miss_pct(proteome)   %>% filter(!feature %in% c("sid","sex","age_in0"))

# Distributions of % missing across features
ggplot(metab_miss, aes(x = pct)) +
  geom_histogram(bins = 30) +
  labs(title = "Metabolome: % missing per feature", x = "% missing", y = "Number of features")

ggplot(proteo_miss, aes(x = pct)) +
  geom_histogram(bins = 30) +
  labs(title = "Proteome: % missing per feature", x = "% missing", y = "Number of features")

# --- 4) Missingness HEATMAPS (sampled; % scale; no row/col names) ---
set.seed(42)

n_rows <- min(400, nrow(metabolome))
n_metab_cols <- min(100, ncol(metabolome)-3)
n_prot_cols  <- min(200, ncol(proteome)-3)

metab_feat <- metabolome %>% select(-sid, -sex, -age_in0)
proteo_feat <- proteome   %>% select(-sid, -sex, -age_in0)

rows_idx <- sample(seq_len(nrow(metabolome)), n_rows)
metab_cols <- sample(names(metab_feat), n_metab_cols)
prot_cols  <- sample(names(proteo_feat), n_prot_cols)

# Convert to percentage (0 or 100 per cell) so legend is in %
metab_miss_mat_pct  <- is.na(metabolome[rows_idx, metab_cols, drop = FALSE]) * 100
proteo_miss_mat_pct <- is.na(proteome[rows_idx,  prot_cols,  drop = FALSE]) * 100

bk <- seq(0, 100, by = 10)

pheatmap(metab_miss_mat_pct,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         legend = TRUE,
         legend_breaks = c(0, 50, 100),
         legend_labels = c("0%", "50%", "100%"),
         breaks = bk,
         main = "Metabolome missingness (%), sampled")

pheatmap(proteo_miss_mat_pct,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         legend = TRUE,
         legend_breaks = c(0, 50, 100),
         legend_labels = c("0%", "50%", "100%"),
         breaks = bk,
         main = "Proteome missingness (%), sampled")