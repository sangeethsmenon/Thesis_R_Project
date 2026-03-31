# ------------------------------------------------------------
# Script 3: Block 5b — Per-module integrative report
# Combines:
#   - protein modules
#   - GO / KEGG enrichment
#   - metabolite associations
# ------------------------------------------------------------

# ---------------------------
# 0) Packages
# ---------------------------
packages <- c("dplyr", "readr", "tibble", "stringr", "purrr")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

# ---------------------------
# 1) Paths
# ---------------------------
resdir <- "results"

fp_modules   <- file.path(resdir, "block5b_protein_modules.csv")
fp_pm        <- file.path(resdir, "block5b_pm_edges_with_protein_module.csv")
fp_metab_enr <- file.path(resdir, "block5b_metabolite_module_enrichment.csv")
fp_go        <- file.path(resdir, "block5b_enrichment_go_all.csv")
fp_kegg      <- file.path(resdir, "block5b_enrichment_kegg_all.csv")

out_dir <- file.path(resdir, "block5b_module_reports")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 2) Load data
# ---------------------------
prot_modules <- read_csv(fp_modules, show_col_types = FALSE)
pm_edges     <- read_csv(fp_pm, show_col_types = FALSE)
metab_enr    <- read_csv(fp_metab_enr, show_col_types = FALSE)

go_all   <- if (file.exists(fp_go))   read_csv(fp_go, show_col_types = FALSE)   else tibble()
kegg_all <- if (file.exists(fp_kegg)) read_csv(fp_kegg, show_col_types = FALSE) else tibble()

modules <- sort(unique(prot_modules$module_protein))

cat("Loaded",
    length(modules), "modules\n")

# ---------------------------
# 3) Helper functions
# ---------------------------
collapse_top <- function(x, n = 5) {
  if (length(x) == 0) return(NA_character_)
  paste(head(unique(x), n), collapse = "; ")
}

# ---------------------------
# 4) Build per-module summaries
# ---------------------------
module_summaries <- map_dfr(modules, function(m) {
  
  # Proteins
  prots <- prot_modules %>%
    filter(module_protein == m) %>%
    pull(protein)
  
  # GO terms
  go_terms <- go_all %>%
    filter(module_protein == m) %>%
    arrange(p.adjust) %>%
    pull(Description)
  
  # KEGG pathways
  kegg_terms <- kegg_all %>%
    filter(module_protein == m) %>%
    arrange(p.adjust) %>%
    pull(Description)
  
  # Metabolites
  mets <- metab_enr %>%
    filter(module_protein == m, q_fisher <= 0.05) %>%
    arrange(q_fisher) %>%
    pull(metabolite_label)
  
  # Build one-row summary
  tibble(
    module_protein       = m,
    n_proteins           = length(prots),
    proteins             = paste(prots, collapse = "; "),
    top_GO_terms         = collapse_top(go_terms, 6),
    top_KEGG_pathways    = collapse_top(kegg_terms, 6),
    associated_metabolites = collapse_top(mets, 8)
  )
})

# ---------------------------
# 5) Write master table
# ---------------------------
write_csv(
  module_summaries,
  file.path(resdir, "block5b_module_integrative_summary.csv")
)

# ---------------------------
# 6) Write per-module files
# ---------------------------
walk(modules, function(m) {
  
  prots <- prot_modules %>%
    filter(module_protein == m)
  
  go_m <- go_all %>%
    filter(module_protein == m)
  
  kegg_m <- kegg_all %>%
    filter(module_protein == m)
  
  mets <- metab_enr %>%
    filter(module_protein == m, q_fisher <= 0.05) %>%
    arrange(q_fisher)
  
  write_csv(
    tibble(
      section = c(
        "Proteins",
        "Top GO terms",
        "Top KEGG pathways",
        "Associated metabolites"
      ),
      content = c(
        paste(prots$protein, collapse = "; "),
        collapse_top(go_m$Description, 15),
        collapse_top(kegg_m$Description, 15),
        collapse_top(mets$metabolite_label, 15)
      )
    ),
    file.path(out_dir, sprintf("module_%02d_summary.csv", m))
  )
})

cat(
  "Script 3 complete.\n",
  "Outputs:\n",
  " - results/block5b_module_integrative_summary.csv\n",
  " - results/block5b_module_reports/\n"
)
