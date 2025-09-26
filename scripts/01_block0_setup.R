
# BLOCK 0 Packages, thresholds, filtering (features & samples)

packages <- c("tidyverse","data.table","VIM","splines","matrixStats")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")
library(here)

proj_root <- here::here()

fs::dir_create(here::here("R"))
fs::dir_create(here::here("scripts"))
fs::dir_create(here::here("data","raw"))
fs::dir_create(here::here("data","processed"))
fs::dir_create(here::here("data","metadata"))
fs::dir_create(here::here("results"))
fs::dir_create(here::here("figures"))
fs::dir_create(here::here("reports"))
fs::dir_create(here::here("cache"))
fs::dir_create(here::here("logs"))

data_root <- Sys.getenv("THESIS_DATA_DIR", unset = here::here("data","raw"))

config <- list(
  project_name = basename(proj_root),
  seed = 42,
  paths = list(
    data_raw = here::here("data","raw"),
    data_processed = here::here("data","processed"),
    results = here::here("results"),
    figures = here::here("figures"),
    logs = here::here("logs")
  ),
  files = list(
    metabolome = "NMR_metabolome_n249_SN13411_20250829.txt",
    proteome   = "NMR_proteome_n2923_SN13411_20250829.txt"
  ),
  thresholds = list(metab_feat = 0.05, prot_feat = 0.30, prot_row = 0.40)
)

set.seed(config$seed)

metabolome <- data.table::fread(file.path(data_root, config$files$metabolome), sep = "\t", header = TRUE) |> as.data.frame()
proteome   <- data.table::fread(file.path(data_root, config$files$proteome),   sep = "\t", header = TRUE) |> as.data.frame()

stopifnot(all(c("sid","sex","age_in0") %in% names(metabolome)))
stopifnot(all(c("sid","sex","age_in0") %in% names(proteome)))
stopifnot(nrow(metabolome) == nrow(proteome))
stopifnot(all(metabolome$sid == proteome$sid))
metabolome$sex <- factor(metabolome$sex)
proteome$sex   <- factor(proteome$sex)

jsonlite::write_json(config, file.path(config$paths$results, "config.json"), auto_unbox = TRUE, pretty = TRUE)
saveRDS(config, file.path(config$paths$results, "config.rds"))
capture.output(sessionInfo(), file = file.path(config$paths$logs, "sessionInfo.txt"))

cat("Loaded metabolome:", nrow(metabolome), "rows ×", ncol(metabolome), "cols\n")
cat("Loaded proteome:",   nrow(proteome),   "rows ×", ncol(proteome),   "cols\n")
cat("Data root:", data_root, "\n")







