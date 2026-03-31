packages <- c("tidyverse","here","readr")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

raw_dir  <- here::here("data","raw")
proc_dir <- here::here("data","processed")
res_dir  <- here::here("results")
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir,  showWarnings = FALSE, recursive = TRUE)

proteome <- readr::read_tsv(
  file.path(raw_dir, "NMR_proteome_n2923_SN13411_20250829.txt"),
  show_col_types = FALSE
)
metabolome <- readr::read_tsv(
  file.path(raw_dir, "NMR_metabolome_n249_SN13411_20250829.txt"),
  show_col_types = FALSE
)

stopifnot(all(metabolome$sid == proteome$sid))

prot_mat <- proteome %>% select(-sid, -sex, -age_in0)
samp_miss_frac <- rowMeans(is.na(prot_mat))
samp_miss_tbl <- tibble(
  sid = proteome$sid,
  pct_missing = 100 * samp_miss_frac
)
readr::write_csv(samp_miss_tbl, file.path(res_dir, "block1_sample_missingness_proteome.csv"))

thr <- 0.50   # keep samples with ≤50% NPX missing (change here for sensitivity)
keep <- samp_miss_frac <= thr

summary_tbl <- tibble(
  rule     = "proteome_missingness_per_sample <= 50%",
  n_before = nrow(proteome),
  n_drop   = sum(!keep),
  n_after  = sum(keep)
)
readr::write_csv(summary_tbl, file.path(res_dir, "block1_sample_filter_summary.csv"))

proteome_qc   <- proteome[keep, , drop = FALSE]
metabolome_qc <- metabolome[keep, , drop = FALSE]
stopifnot(all(metabolome_qc$sid == proteome_qc$sid))

readr::write_tsv(proteome_qc,   file.path(proc_dir, "proteome_qc.tsv"))
readr::write_tsv(metabolome_qc, file.path(proc_dir, "metabolome_qc.tsv"))

cat("Block 1 (sample filter) complete.\n")
print(summary_tbl)
