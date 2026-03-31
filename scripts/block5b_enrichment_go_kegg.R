# ------------------------------------------------------------
# Script 2: Block 5b (Part 2) — GO / KEGG enrichment on protein modules
# Inputs:
#   results/block5b_protein_modules.csv   (from Script 1)
# Outputs:
#   results/block5b_enrichment_go_all.csv
#   results/block5b_enrichment_kegg_all.csv
#   plus per-module files in:
#     results/block5b_enrichment_go/
#     results/block5b_enrichment_kegg/
# Notes:
# - Assumes `protein` column contains human gene SYMBOLS (e.g., LDLR, APOA1, AAMDC).
# - If your identifiers are NOT gene symbols, you must provide a mapping file first.
# ------------------------------------------------------------

# ---------------------------
# 0) Packages
# ---------------------------
cran_pkgs <- c("dplyr", "readr", "tibble", "stringr")
need_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need_cran)) install.packages(need_cran, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# Bioconductor packages (should already be installed in your renv; we check + fail clearly)
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc)) {
  stop(
    "Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "), "\n",
    "Install with:\n",
    "  install.packages('BiocManager')\n",
    "  BiocManager::install(c('clusterProfiler','org.Hs.eg.db','AnnotationDbi'))\n"
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---------------------------
# 1) Paths
# ---------------------------
resdir <- "results"
in_fp  <- file.path(resdir, "block5b_protein_modules.csv")

if (!file.exists(in_fp)) {
  stop("Cannot find: ", in_fp, "\nRun Script 1 first to generate block5b_protein_modules.csv")
}

out_go_dir   <- file.path(resdir, "block5b_enrichment_go")
out_kegg_dir <- file.path(resdir, "block5b_enrichment_kegg")
dir.create(out_go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_kegg_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 2) Load protein modules
# ---------------------------
prot_modules <- readr::read_csv(in_fp, show_col_types = FALSE)

req_cols <- c("protein", "module_protein")
if (!all(req_cols %in% names(prot_modules))) {
  stop("Input file must contain columns: ", paste(req_cols, collapse = ", "),
       "\nFound: ", paste(names(prot_modules), collapse = ", "))
}

prot_modules <- prot_modules %>%
  mutate(
    protein = as.character(protein),
    module_protein = as.integer(module_protein)
  ) %>%
  filter(!is.na(protein), protein != "", !is.na(module_protein)) %>%
  distinct(protein, module_protein, .keep_all = TRUE)

cat("Loaded protein modules:", nrow(prot_modules), "rows; modules:",
    n_distinct(prot_modules$module_protein), "\n")

# ---------------------------
# 3) Map SYMBOL -> ENTREZID
# ---------------------------
# This is the correct mapping for your case (protein names are gene symbols).
keys <- sort(unique(prot_modules$protein))

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = keys,
  keytype = "SYMBOL",
  columns = c("ENTREZID")
)

# AnnotationDbi::select returns a data.frame with columns: SYMBOL, ENTREZID
# Keep 1 ENTREZID per SYMBOL (deduplicate)
gene_map <- gene_map %>%
  as_tibble() %>%
  rename(protein = SYMBOL) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(protein, .keep_all = TRUE)

cat("Mapped proteins to ENTREZID:",
    n_distinct(gene_map$protein), "out of", length(keys), "\n")

# Join mapping into module table
prot_modules_entrez <- prot_modules %>%
  left_join(gene_map, by = "protein")

# Optional: save mapping for transparency
write_csv(
  prot_modules_entrez %>% select(protein, module_protein, ENTREZID) %>% distinct(),
  file.path(resdir, "block5b_symbol_to_entrez_map.csv")
)

# ---------------------------
# 4) Enrichment parameters
# ---------------------------
min_proteins_for_enrich <- 10     # change to 15/20 if you want stricter modules only
q_cutoff <- 0.05                  # BH q-value cutoff used by clusterProfiler
ont_go   <- "BP"                  # BP is most interpretable for thesis; can also try "MF","CC"

# ---------------------------
# 5) Run GO and KEGG per module
# ---------------------------
modules <- sort(unique(prot_modules_entrez$module_protein))

go_all   <- list()
kegg_all <- list()

for (m in modules) {
  ents <- prot_modules_entrez %>%
    filter(module_protein == m) %>%
    pull(ENTREZID) %>%
    unique() %>%
    { .[!is.na(.)] }
  
  if (length(ents) < min_proteins_for_enrich) next
  
  # ---- GO (Biological Process)
  ego <- tryCatch(
    enrichGO(
      gene          = ents,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont_go,
      pAdjustMethod = "BH",
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    df_go <- as_tibble(as.data.frame(ego)) %>%
      mutate(module_protein = m, .before = 1) %>%
      arrange(p.adjust, pvalue)
    
    write_csv(df_go, file.path(out_go_dir, sprintf("go_module_%03d.csv", m)))
    go_all[[as.character(m)]] <- df_go
  }
  
  # ---- KEGG
  ekg <- tryCatch(
    enrichKEGG(
      gene          = ents,
      organism      = "hsa",
      pAdjustMethod = "BH",
      qvalueCutoff  = q_cutoff
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0) {
    df_kegg <- as_tibble(as.data.frame(ekg)) %>%
      mutate(module_protein = m, .before = 1) %>%
      arrange(p.adjust, pvalue)
    
    write_csv(df_kegg, file.path(out_kegg_dir, sprintf("kegg_module_%03d.csv", m)))
    kegg_all[[as.character(m)]] <- df_kegg
  }
}

# ---------------------------
# 6) Write combined outputs + simple summary
# ---------------------------
if (length(go_all)) {
  go_all_df <- bind_rows(go_all)
  write_csv(go_all_df, file.path(resdir, "block5b_enrichment_go_all.csv"))
  cat("GO enrichment: modules with results =", n_distinct(go_all_df$module_protein),
      "; total terms =", nrow(go_all_df), "\n")
} else {
  write_lines("No GO results produced (maybe modules too small or mapping too low).",
              file.path(resdir, "block5b_enrichment_go_NONE.txt"))
  cat("GO enrichment: no results.\n")
}

if (length(kegg_all)) {
  kegg_all_df <- bind_rows(kegg_all)
  write_csv(kegg_all_df, file.path(resdir, "block5b_enrichment_kegg_all.csv"))
  cat("KEGG enrichment: modules with results =", n_distinct(kegg_all_df$module_protein),
      "; total pathways =", nrow(kegg_all_df), "\n")
} else {
  write_lines("No KEGG results produced (maybe modules too small or mapping too low).",
              file.path(resdir, "block5b_enrichment_kegg_NONE.txt"))
  cat("KEGG enrichment: no results.\n")
}

cat(
  "Script 2 complete.\n",
  "Outputs:\n",
  " - results/block5b_symbol_to_entrez_map.csv\n",
  " - results/block5b_enrichment_go_all.csv (and per-module files in results/block5b_enrichment_go/)\n",
  " - results/block5b_enrichment_kegg_all.csv (and per-module files in results/block5b_enrichment_kegg/)\n"
)
