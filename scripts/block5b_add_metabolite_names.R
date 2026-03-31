# --- 06b_block5b_add_metabolite_names.R ---------------------------------------
# Adds human-readable metabolite names (from metabolite_info.txt)
# to Block 5B outputs:
#   - block5b_pm_edges_with_protein_module.csv
#   - block5b_metabolite_module_enrichment.csv
# Writes:
#   - block5b_pm_edges_with_protein_module_named.csv
#   - block5b_metabolite_module_enrichment_named.csv
# -----------------------------------------------------------------------------

packages <- c("tidyverse", "here", "jsonlite")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# ---- Load config (same pattern as your other scripts) ----
if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found. Run Block 0.")
}

stopifnot(!is.null(config$paths$results))
dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)

# ---- Locate metabolite mapping file ----
# Prefer config$paths$data_raw if it exists; otherwise fallback to data/raw
data_raw <- config$paths$data_raw
if (is.null(data_raw)) data_raw <- here::here("data","raw")

map_fp <- file.path(data_raw, "metabolite_info.txt")
if (!file.exists(map_fp)) {
  stop("Could not find metabolite_info.txt at: ", map_fp,
       "\nCheck config$paths$data_raw or place it in data/raw/")
}

# ---- Load metabolite ID -> name mapping ----
met_map <- readr::read_tsv(
  map_fp,
  col_names = c("metabolite", "metabolite_name"),
  show_col_types = FALSE
) %>%
  mutate(
    metabolite = as.character(metabolite),
    metabolite_name = as.character(metabolite_name)
  ) %>%
  distinct(metabolite, .keep_all = TRUE)

cat("Loaded metabolite mapping rows:", nrow(met_map), "\n")

# ---- Helper: annotate a CSV with metabolite_name ----
annotate_metabolites <- function(infile, outfile) {
  in_fp  <- file.path(config$paths$results, infile)
  out_fp <- file.path(config$paths$results, outfile)
  
  if (!file.exists(in_fp)) stop("Input file not found: ", in_fp)
  
  df <- readr::read_csv(in_fp, show_col_types = FALSE)
  
  if (!"metabolite" %in% names(df)) {
    stop("File does not have a 'metabolite' column: ", infile)
  }
  
  df2 <- df %>%
    mutate(metabolite = as.character(metabolite)) %>%
    left_join(met_map, by = "metabolite") %>%
    relocate(metabolite_name, .after = metabolite)
  
  # quick sanity check: unmapped IDs
  n_miss <- sum(is.na(df2$metabolite_name))
  if (n_miss > 0) {
    cat("WARNING:", n_miss, "metabolites did not map in", infile, "\n")
    cat("Example unmapped IDs:\n")
    print(head(df2$metabolite[is.na(df2$metabolite_name)], 10))
  } else {
    cat("All metabolites mapped in", infile, "\n")
  }
  
  readr::write_csv(df2, out_fp)
  cat("Wrote:", out_fp, "\n")
}

# ---- Apply to Block 5B outputs ----
annotate_metabolites(
  infile  = "block5b_pm_edges_with_protein_module.csv",
  outfile = "block5b_pm_edges_with_protein_module_named.csv"
)

annotate_metabolites(
  infile  = "block5b_metabolite_module_enrichment.csv",
  outfile = "block5b_metabolite_module_enrichment_named.csv"
)

cat("Done.\n")
# -----------------------------------------------------------------------------
