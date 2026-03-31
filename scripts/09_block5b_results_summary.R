#!/usr/bin/env Rscript

# 09_block5b_results_summary.R
# Purpose: Generate a compact, thesis-ready “Results snapshot” (Markdown + JSON)
# from the Block5b protein-first + enrichment outputs, without sharing large CSVs.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(glue)
  library(jsonlite)
})

# ========== USER SETTINGS (edit if needed) ==========
Q_PP   <- 0.01
R_PP   <- 0.30
Q_FISH <- 0.05
ONT_GO <- "BP"

TOP_MODULES_N <- 8
TOP_ENRICH_N  <- 10
TOP_GO_PER_MOD   <- 5
TOP_KEGG_PER_MOD <- 5

# ========== INPUT PATHS ==========
in_paths <- list(
  pp_edges   = here("results", "block5b_pp_edges.csv"),
  modules    = here("results", "block5b_protein_modules.csv"),
  pm_mapped  = here("results", "block5b_pm_edges_with_protein_module.csv"),
  fish       = here("results", "block5b_metabolite_module_enrichment.csv"),
  map_entrez = here("results", "block5b_symbol_to_entrez_map.csv"),
  go_all     = here("results", "block5b_enrichment_go_all.csv"),
  kegg_all   = here("results", "block5b_enrichment_kegg_all.csv"),
  integr     = here("results", "block5b_module_integrative_summary.csv")
)

out_md   <- here("results", "block5b_results_summary.md")
out_json <- here("results", "block5b_results_summary.json")
out_dir  <- here("results", "block5b_results_summary_tables")

# ========== HELPERS ==========
read_csv_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
           error = function(e) NULL)
}

pick_col <- function(df, candidates) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

as_int <- function(x) suppressWarnings(as.integer(x))
as_num <- function(x) suppressWarnings(as.numeric(x))

fmt_int <- function(x) if (is.na(x)) "NA" else format(as.integer(x), big.mark = ",")
fmt_num <- function(x, digits = 3) if (is.na(x)) "NA" else format(round(as.numeric(x), digits), trim = TRUE)

ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

topn_table <- function(df, n = 10) {
  if (is.null(df) || nrow(df) == 0) return(df)
  head(df, n)
}

# ========== LOAD ==========
pp_edges  <- read_csv_safely(in_paths$pp_edges)
modules   <- read_csv_safely(in_paths$modules)
pm_mapped <- read_csv_safely(in_paths$pm_mapped)
fish      <- read_csv_safely(in_paths$fish)
map_entrez<- read_csv_safely(in_paths$map_entrez)
go_all    <- read_csv_safely(in_paths$go_all)
kegg_all  <- read_csv_safely(in_paths$kegg_all)
integr    <- read_csv_safely(in_paths$integr)

missing_inputs <- names(in_paths)[!file.exists(unlist(in_paths))]

# ========== SUMMARIES ==========

# 1) PP network size
pp_summary <- list(
  pp_edges_n = if (!is.null(pp_edges)) nrow(pp_edges) else NA_integer_,
  pp_nodes_n = {
    if (is.null(pp_edges)) NA_integer_ else {
      from_col <- pick_col(pp_edges, c("from", "protein1", "protein_a", "node1", "src"))
      to_col   <- pick_col(pp_edges, c("to", "protein2", "protein_b", "node2", "dst"))
      if (is.na(from_col) || is.na(to_col)) NA_integer_
      else length(unique(c(pp_edges[[from_col]], pp_edges[[to_col]])))
    }
  }
)

# 2) Module counts and sizes
module_cols <- list(
  protein = if (!is.null(modules)) pick_col(modules, c("protein", "Protein", "symbol", "SYMBOL", "node", "id")) else NA_character_,
  module  = if (!is.null(modules)) pick_col(modules, c("module_protein", "module", "community", "cluster", "louvain", "group")) else NA_character_
)

module_size_tbl <- NULL
modules_summary <- list(
  modules_n = NA_integer_,
  proteins_in_modules_n = NA_integer_,
  module_size_min = NA_integer_,
  module_size_median = NA_integer_,
  module_size_max = NA_integer_
)

if (!is.null(modules) && !is.na(module_cols$protein) && !is.na(module_cols$module)) {
  module_size_tbl <- modules %>%
    filter(!is.na(.data[[module_cols$module]]), !is.na(.data[[module_cols$protein]])) %>%
    group_by(module = .data[[module_cols$module]]) %>%
    summarise(n_proteins = n_distinct(.data[[module_cols$protein]]), .groups = "drop") %>%
    arrange(desc(n_proteins))
  
  msizes <- module_size_tbl$n_proteins
  modules_summary <- list(
    modules_n = nrow(module_size_tbl),
    proteins_in_modules_n = sum(module_size_tbl$n_proteins),
    module_size_min = min(msizes),
    module_size_median = as.integer(median(msizes)),
    module_size_max = max(msizes)
  )
}

top_modules_tbl <- if (!is.null(module_size_tbl)) topn_table(module_size_tbl, TOP_MODULES_N) else NULL

# 3) PM edges mapped to modules
pm_mapped_summary <- list(
  pm_edges_mapped_n = if (!is.null(pm_mapped)) nrow(pm_mapped) else NA_integer_,
  pm_unique_metabolites_n = {
    if (is.null(pm_mapped)) NA_integer_ else {
      mcol <- pick_col(pm_mapped, c("metabolite_label", "metabolite", "met", "FID", "met_id"))
      if (is.na(mcol)) NA_integer_ else n_distinct(pm_mapped[[mcol]])
    }
  },
  pm_unique_proteins_n = {
    if (is.null(pm_mapped)) NA_integer_ else {
      pcol <- pick_col(pm_mapped, c("protein", "Protein", "symbol", "SYMBOL", "prot"))
      if (is.na(pcol)) NA_integer_ else n_distinct(pm_mapped[[pcol]])
    }
  },
  pm_unique_modules_n = {
    if (is.null(pm_mapped)) NA_integer_ else {
      mcol <- pick_col(pm_mapped, c("module_protein", "module", "community", "cluster"))
      if (is.na(mcol)) NA_integer_ else n_distinct(pm_mapped[[mcol]])
    }
  }
)

# 4) Fisher metabolite->module enrichment (significant + top pairs)
fish_summary <- list(
  fish_sig_pairs_n = NA_integer_,
  fish_sig_unique_metabolites_n = NA_integer_,
  fish_sig_unique_modules_n = NA_integer_
)

fish_top_tbl <- NULL

if (!is.null(fish)) {
  qcol <- pick_col(fish, c("q_fisher", "q", "qvalue", "padj", "p.adjust"))
  mcol <- pick_col(fish, c("metabolite_label", "metabolite", "met", "FID", "met_id"))
  modc <- pick_col(fish, c("module_protein", "module", "community", "cluster"))
  
  if (!is.na(qcol)) {
    fish2 <- fish %>% mutate(qv = as_num(.data[[qcol]]))
    sig  <- fish2 %>% filter(!is.na(qv), qv <= Q_FISH)
    
    fish_summary <- list(
      fish_sig_pairs_n = nrow(sig),
      fish_sig_unique_metabolites_n = if (!is.na(mcol) && nrow(sig)>0) n_distinct(sig[[mcol]]) else NA_integer_,
      fish_sig_unique_modules_n = if (!is.na(modc) && nrow(sig)>0) n_distinct(sig[[modc]]) else NA_integer_
    )
    
    # keep interpretable columns if present
    keep_cols <- c(mcol, modc, qcol, "k", "K", "n", "N", "p_fisher", "p")
    keep_cols <- keep_cols[!is.na(keep_cols) & keep_cols %in% names(fish)]
    fish_top_tbl <- sig %>%
      arrange(qv) %>%
      select(any_of(keep_cols)) %>%
      head(TOP_ENRICH_N)
  }
}

# 5) Mapping success: SYMBOL -> ENTREZ
map_summary <- list(
  mapping_total_symbols = NA_integer_,
  mapping_mapped_entrez = NA_integer_,
  mapping_rate = NA_real_
)

if (!is.null(map_entrez)) {
  sym <- pick_col(map_entrez, c("SYMBOL", "symbol", "protein", "Protein"))
  ent <- pick_col(map_entrez, c("ENTREZID", "entrez", "entrez_id", "ENTREZ"))
  if (!is.na(sym) && !is.na(ent)) {
    total <- n_distinct(map_entrez[[sym]])
    mapped <- map_entrez %>% filter(!is.na(.data[[ent]]), .data[[ent]] != "") %>% summarise(n = n_distinct(.data[[sym]])) %>% pull(n)
    map_summary <- list(
      mapping_total_symbols = total,
      mapping_mapped_entrez = mapped,
      mapping_rate = if (total > 0) mapped / total else NA_real_
    )
  }
}

# 6) GO + KEGG summaries (modules with >=1 term, totals, top terms per module)
enrich_go_summary <- list(modules_with_terms = NA_integer_, total_terms = NA_integer_)
go_top_by_module_tbl <- NULL

if (!is.null(go_all) && nrow(go_all) > 0) {
  modc <- pick_col(go_all, c("module_protein", "module", "community", "cluster"))
  pcol <- pick_col(go_all, c("p.adjust", "padj", "qvalue", "q", "p_adj"))
  dcol <- pick_col(go_all, c("Description", "description", "term"))
  if (!is.na(modc) && !is.na(pcol) && !is.na(dcol)) {
    go2 <- go_all %>% mutate(padj = as_num(.data[[pcol]]))
    enrich_go_summary <- list(
      modules_with_terms = n_distinct(go2[[modc]]),
      total_terms = nrow(go2)
    )
    go_top_by_module_tbl <- go2 %>%
      filter(!is.na(padj)) %>%
      arrange(.data[[modc]], padj) %>%
      group_by(module = .data[[modc]]) %>%
      slice_head(n = TOP_GO_PER_MOD) %>%
      ungroup() %>%
      transmute(module,
                term = .data[[dcol]],
                p_adjust = padj)
  }
}

enrich_kegg_summary <- list(modules_with_terms = NA_integer_, total_terms = NA_integer_)
kegg_top_by_module_tbl <- NULL

if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
  modc <- pick_col(kegg_all, c("module_protein", "module", "community", "cluster"))
  pcol <- pick_col(kegg_all, c("p.adjust", "padj", "qvalue", "q", "p_adj"))
  dcol <- pick_col(kegg_all, c("Description", "description", "term", "Pathway", "pathway"))
  if (!is.na(modc) && !is.na(pcol) && !is.na(dcol)) {
    k2 <- kegg_all %>% mutate(padj = as_num(.data[[pcol]]))
    enrich_kegg_summary <- list(
      modules_with_terms = n_distinct(k2[[modc]]),
      total_terms = nrow(k2)
    )
    kegg_top_by_module_tbl <- k2 %>%
      filter(!is.na(padj)) %>%
      arrange(.data[[modc]], padj) %>%
      group_by(module = .data[[modc]]) %>%
      slice_head(n = TOP_KEGG_PER_MOD) %>%
      ungroup() %>%
      transmute(module,
                pathway = .data[[dcol]],
                p_adjust = padj)
  }
}

# 7) Integrative summary (already “thesis-ready” but we extract top lines)
integr_summary <- list(
  integr_rows = if (!is.null(integr)) nrow(integr) else NA_integer_,
  integr_modules = {
    if (is.null(integr)) NA_integer_ else {
      modc <- pick_col(integr, c("module_protein", "module", "community", "cluster"))
      if (is.na(modc)) NA_integer_ else n_distinct(integr[[modc]])
    }
  }
)

# ========== WRITE SMALL TABLES ==========
ensure_dir(out_dir)

if (!is.null(top_modules_tbl)) {
  readr::write_csv(top_modules_tbl, file.path(out_dir, "top_modules_by_size.csv"))
}
if (!is.null(fish_top_tbl)) {
  readr::write_csv(fish_top_tbl, file.path(out_dir, "top_metabolite_module_enrichment.csv"))
}
if (!is.null(go_top_by_module_tbl)) {
  readr::write_csv(go_top_by_module_tbl, file.path(out_dir, "top_GO_terms_by_module.csv"))
}
if (!is.null(kegg_top_by_module_tbl)) {
  readr::write_csv(kegg_top_by_module_tbl, file.path(out_dir, "top_KEGG_pathways_by_module.csv"))
}

# ========== BUILD MARKDOWN ==========
md_lines <- c(
  "# Block5b Results Summary (protein-first modules + enrichment)",
  "",
  glue("Generated: {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"),
  "",
  "## Inputs checked",
  "",
  if (length(missing_inputs) == 0) "All expected inputs were found." else
    glue("Missing inputs: {paste(missing_inputs, collapse = ', ')}"),
  "",
  "## Protein–protein network (PP)",
  "",
  glue("Thresholds (reported): q ≤ {Q_PP}, |r| ≥ {R_PP}."),
  glue("PP edges retained: {fmt_int(pp_summary$pp_edges_n)}."),
  glue("PP nodes (unique proteins in PP edges): {fmt_int(pp_summary$pp_nodes_n)}."),
  "",
  "## Protein modules (Louvain on PP graph)",
  "",
  glue("Number of modules: {fmt_int(modules_summary$modules_n)}."),
  glue("Proteins assigned to modules: {fmt_int(modules_summary$proteins_in_modules_n)}."),
  glue("Module size (min / median / max): {fmt_int(modules_summary$module_size_min)} / {fmt_int(modules_summary$module_size_median)} / {fmt_int(modules_summary$module_size_max)}."),
  "",
  glue("Top modules by size (see `{file.path('results','block5b_results_summary_tables','top_modules_by_size.csv')}` for the table)."),
  "",
  "## PM edges mapped onto protein modules",
  "",
  glue("PM edges with module annotation: {fmt_int(pm_mapped_summary$pm_edges_mapped_n)}."),
  glue("Unique metabolites represented: {fmt_int(pm_mapped_summary$pm_unique_metabolites_n)}."),
  glue("Unique proteins represented: {fmt_int(pm_mapped_summary$pm_unique_proteins_n)}."),
  glue("Unique modules represented: {fmt_int(pm_mapped_summary$pm_unique_modules_n)}."),
  "",
  "## Metabolite → module enrichment (Fisher)",
  "",
  glue("Significance cutoff (reported): q ≤ {Q_FISH}."),
  glue("Significant metabolite–module pairs: {fmt_int(fish_summary$fish_sig_pairs_n)}."),
  glue("Unique significant metabolites: {fmt_int(fish_summary$fish_sig_unique_metabolites_n)}."),
  glue("Unique significant modules: {fmt_int(fish_summary$fish_sig_unique_modules_n)}."),
  "",
  glue("Top signals table: `{file.path('results','block5b_results_summary_tables','top_metabolite_module_enrichment.csv')}`."),
  "",
  "## ID mapping (SYMBOL → ENTREZ)",
  "",
  glue("Symbols in mapping table: {fmt_int(map_summary$mapping_total_symbols)}."),
  glue("Symbols mapped to an ENTREZ ID: {fmt_int(map_summary$mapping_mapped_entrez)}."),
  glue("Mapping rate: {fmt_num(100 * map_summary$mapping_rate, 1)}%."),
  "",
  "## GO enrichment (clusterProfiler)",
  "",
  glue("Ontology (reported): {ONT_GO}."),
  glue("Modules with ≥1 GO term (in results table): {fmt_int(enrich_go_summary$modules_with_terms)}."),
  glue("Total GO terms reported (all modules combined): {fmt_int(enrich_go_summary$total_terms)}."),
  glue("Top GO terms per module table: `{file.path('results','block5b_results_summary_tables','top_GO_terms_by_module.csv')}`."),
  "",
  "## KEGG enrichment (clusterProfiler)",
  "",
  glue("Modules with ≥1 KEGG pathway (in results table): {fmt_int(enrich_kegg_summary$modules_with_terms)}."),
  glue("Total KEGG pathways reported (all modules combined): {fmt_int(enrich_kegg_summary$total_terms)}."),
  glue("Top KEGG pathways per module table: `{file.path('results','block5b_results_summary_tables','top_KEGG_pathways_by_module.csv')}`."),
  "",
  "## Integrative summary table",
  "",
  glue("Rows in `block5b_module_integrative_summary.csv`: {fmt_int(integr_summary$integr_rows)}."),
  glue("Distinct modules in integrative summary: {fmt_int(integr_summary$integr_modules)}."),
  "",
  "## Notes for thesis writing",
  "",
  "Use this file for the Results chapter narrative (counts + headline findings), and move full enrichment tables and full module tables to the Appendix. The small tables in `results/block5b_results_summary_tables/` are already in ‘top-only’ form for main-text reporting."
)

writeLines(md_lines, out_md)

# ========== BUILD JSON (machine-readable snapshot) ==========
snapshot <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  thresholds = list(q_pp = Q_PP, r_pp = R_PP, q_fisher = Q_FISH, ont_go = ONT_GO),
  inputs = list(paths = in_paths, missing = missing_inputs),
  pp = pp_summary,
  modules = modules_summary,
  pm_mapped = pm_mapped_summary,
  fisher = fish_summary,
  mapping = map_summary,
  go = enrich_go_summary,
  kegg = enrich_kegg_summary,
  integrative = integr_summary,
  outputs = list(
    markdown = out_md,
    json = out_json,
    tables_dir = out_dir
  )
)

write_json(snapshot, out_json, pretty = TRUE, auto_unbox = TRUE)

cat(glue("Wrote:\n- {out_md}\n- {out_json}\n- {out_dir}/(top tables)\n"))
