# ------------------------------------------------------------
# Script 4 (OPTIONAL): Pathway → Protein modules → Metabolites
# Purpose:
#   For each GO / KEGG pathway, list:
#     - enriched protein modules
#     - proteins in the pathway
#     - metabolites associated with those proteins
#
# This is a pathway-centric integrative view.
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

resdir <- "results"

# ---------------------------
# 1) Load required inputs
# ---------------------------
go_fp   <- file.path(resdir, "block5b_enrichment_go_all.csv")
kegg_fp <- file.path(resdir, "block5b_enrichment_kegg_all.csv")
pm_fp   <- file.path(resdir, "block5b_pm_edges_with_protein_module.csv")

stopifnot(file.exists(go_fp), file.exists(kegg_fp), file.exists(pm_fp))

go_all   <- readr::read_csv(go_fp, show_col_types = FALSE)
kegg_all <- readr::read_csv(kegg_fp, show_col_types = FALSE)
pm       <- readr::read_csv(pm_fp, show_col_types = FALSE)

# Safety checks
req_pm <- c("protein", "module_protein", "metabolite_label")
if (!all(req_pm %in% names(pm))) {
  stop("PM file must contain: ", paste(req_pm, collapse = ", "))
}

# ---------------------------
# 2) Helper: expand gene list
# ---------------------------
# clusterProfiler stores genes as "A/B/C"
expand_genes <- function(df) {
  df %>%
    mutate(gene = str_split(geneID, "/")) %>%
    unnest(gene)
}

# ---------------------------
# 3) GO → proteins → metabolites
# ---------------------------
go_long <- expand_genes(go_all)

go_pathway_map <- go_long %>%
  inner_join(pm, by = c("gene" = "protein", "module_protein")) %>%
  distinct(
    ontology = "GO",
    pathway_id = ID,
    pathway_description = Description,
    module_protein,
    protein = gene,
    metabolite_label
  )

# ---------------------------
# 4) KEGG → proteins → metabolites
# ---------------------------
kegg_long <- expand_genes(kegg_all)

kegg_pathway_map <- kegg_long %>%
  inner_join(pm, by = c("gene" = "protein", "module_protein")) %>%
  distinct(
    ontology = "KEGG",
    pathway_id = ID,
    pathway_description = Description,
    module_protein,
    protein = gene,
    metabolite_label
  )

# ---------------------------
# 5) Combine + summarize
# ---------------------------
pathway_map <- bind_rows(go_pathway_map, kegg_pathway_map)

# Compact summary per pathway + module
pathway_summary <- pathway_map %>%
  group_by(
    ontology,
    pathway_id,
    pathway_description,
    module_protein
  ) %>%
  summarise(
    n_proteins    = n_distinct(protein),
    n_metabolites = n_distinct(metabolite_label),
    proteins      = paste(sort(unique(protein)), collapse = ", "),
    metabolites   = paste(sort(unique(metabolite_label)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(ontology, module_protein, desc(n_metabolites))

# ---------------------------
# 6) Write outputs
# ---------------------------
write_csv(
  pathway_map,
  file.path(resdir, "block5b_pathway_protein_metabolite_map.csv")
)

write_csv(
  pathway_summary,
  file.path(resdir, "block5b_pathway_module_summary.csv")
)

cat(
  "Script 4 complete.\n",
  "Outputs:\n",
  " - results/block5b_pathway_protein_metabolite_map.csv\n",
  " - results/block5b_pathway_module_summary.csv\n"
)
