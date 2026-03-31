suppressPackageStartupMessages({
  library(tidyverse); library(here); library(jsonlite)
})

# --- load config (same logic as your pipeline) ---
if (!exists("config")) {
  cfg_rds  <- here("results","config.rds")
  cfg_json <- here("results","config.json")
  if (file.exists(cfg_rds)) {
    config <- readRDS(cfg_rds)
  } else if (file.exists(cfg_json)) {
    config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  } else stop("Config not found.")
}

resdir <- config$paths$results
figdir <- config$paths$figures
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# ---- expected Block5b outputs (module projection + enrichment) ----
fp_pm_map   <- file.path(resdir, "block5b_pm_edges_with_protein_module.csv")
fp_met_enr  <- file.path(resdir, "block5b_metabolite_module_enrichment.csv")
fp_go_all   <- file.path(resdir, "block5b_enrichment_go_all.csv")
fp_kegg_all <- file.path(resdir, "block5b_enrichment_kegg_all.csv")

stopifnot(file.exists(fp_pm_map), file.exists(fp_met_enr))

pm_map  <- readr::read_csv(fp_pm_map,  show_col_types = FALSE)
met_enr <- readr::read_csv(fp_met_enr, show_col_types = FALSE)

go_all   <- if (file.exists(fp_go_all))   readr::read_csv(fp_go_all,   show_col_types = FALSE) else tibble()
kegg_all <- if (file.exists(fp_kegg_all)) readr::read_csv(fp_kegg_all, show_col_types = FALSE) else tibble()

# -----------------------------
# 1) Core counts for Section 4.4
# -----------------------------
# pm_map expected columns include: protein, metabolite (or metabolite_label), module_protein, etc.
# met_enr expected columns include: module_protein, metabolite_label, p_fisher, q_fisher (names may vary).
col_metlab <- if ("metabolite_label" %in% names(pm_map)) "metabolite_label" else "metabolite"
col_qfish  <- if ("q_fisher" %in% names(met_enr)) "q_fisher" else if ("q" %in% names(met_enr)) "q" else NA_character_
if (is.na(col_qfish)) stop("Cannot find q_fisher (or q) column in block5b_metabolite_module_enrichment.csv")

n_edges_mapped <- nrow(pm_map)
n_modules_with_pm <- pm_map %>% distinct(module_protein) %>% nrow()
n_pm_proteins <- pm_map %>% distinct(protein) %>% nrow()
n_pm_mets <- pm_map %>% distinct(.data[[col_metlab]]) %>% nrow()

n_tests_met_module <- nrow(met_enr)
n_sig_met_module <- met_enr %>% filter(.data[[col_qfish]] <= 0.05) %>% nrow()

cat("\n--- Section 4.4 core counts ---\n")
cat("Mapped PM edges:", n_edges_mapped, "\n")
cat("Proteins with >=1 mapped PM edge:", n_pm_proteins, "\n")
cat("Metabolites with >=1 mapped PM edge:", n_pm_mets, "\n")
cat("Modules receiving >=1 PM edge:", n_modules_with_pm, "\n")
cat("Metabolite–module enrichment tests:", n_tests_met_module, "\n")
cat("Significant metabolite–module pairs (q<=0.05):", n_sig_met_module, "\n")

# Save a LaTeX-friendly summary table (you can include this via \input or copy)
tbl_44_overview <- tibble(
  metric = c(
    "Mapped protein--metabolite edges",
    "Proteins with ≥1 mapped edge",
    "Metabolites with ≥1 mapped edge",
    "Protein modules receiving ≥1 mapped edge",
    "Metabolite--module tests (Fisher)",
    "Significant metabolite--module pairs (q ≤ 0.05)"
  ),
  value = c(
    n_edges_mapped,
    n_pm_proteins,
    n_pm_mets,
    n_modules_with_pm,
    n_tests_met_module,
    n_sig_met_module
  )
)
readr::write_csv(tbl_44_overview, file.path(resdir, "results_sec44_overview_counts.csv"))
cat("Saved:", file.path(resdir, "results_sec44_overview_counts.csv"), "\n")

# -----------------------------------------
# 2) Top enriched metabolite–module pairs
# -----------------------------------------
top_pairs <- met_enr %>%
  mutate(q_fisher = .data[[col_qfish]]) %>%
  arrange(q_fisher) %>%
  slice_head(n = 20)

readr::write_csv(top_pairs, file.path(resdir, "results_sec44_top20_metabolite_module_pairs.csv"))
cat("Saved:", file.path(resdir, "results_sec44_top20_metabolite_module_pairs.csv"), "\n")

# -----------------------------------------
# 3) GO/KEGG compact summaries (if present)
# -----------------------------------------
if (nrow(go_all) > 0) {
  go_mod_summary <- go_all %>%
    group_by(module_protein) %>%
    summarise(n_terms = n(), best_q = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_q)
  readr::write_csv(go_mod_summary, file.path(resdir, "results_sec44_go_module_summary.csv"))
  cat("Saved:", file.path(resdir, "results_sec44_go_module_summary.csv"), "\n")
} else {
  cat("GO file not found (block5b_enrichment_go_all.csv). If you want it, run your GO script.\n")
}

if (nrow(kegg_all) > 0) {
  kegg_mod_summary <- kegg_all %>%
    group_by(module_protein) %>%
    summarise(n_pathways = n(), best_q = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_q)
  readr::write_csv(kegg_mod_summary, file.path(resdir, "results_sec44_kegg_module_summary.csv"))
  cat("Saved:", file.path(resdir, "results_sec44_kegg_module_summary.csv"), "\n")
} else {
  cat("KEGG file not found (block5b_enrichment_kegg_all.csv). If you want it, run your KEGG script.\n")
}



# -----------------------------------------
# 4) Clean heatmap: module × metabolite (sig only)
#    FINAL VERSION: module 3 only + metabolites ranked by score
# -----------------------------------------
hm <- met_enr %>%
  mutate(q_fisher = .data[[col_qfish]]) %>%
  filter(q_fisher <= 0.05) %>%
  transmute(
    module_protein = module_protein,
    metabolite_label = if ("metabolite_label" %in% names(met_enr)) metabolite_label else metabolite,
    score = -log10(pmax(q_fisher, 1e-300))
  )

if (nrow(hm) == 0) {
  cat("No significant metabolite–module pairs at q<=0.05, so skipping heatmap.\n")
} else {
  
  # --- keep only module 3 (dominant enrichment module) ---
  hm3 <- hm %>% filter(module_protein == 3)
  
  if (nrow(hm3) == 0) {
    cat("No significant pairs in module 3, so skipping module-3 heatmap.\n")
  } else {
    
    # --- choose top metabolites by score (readability + ranking) ---
    top_mets <- hm3 %>%
      arrange(desc(score)) %>%
      distinct(metabolite_label, .keep_all = TRUE) %>%  # if duplicates exist, keep best score
      slice_head(n = 30) %>%                            # change 30 -> 20 if you want tighter
      pull(metabolite_label)
    
    hm3_plot <- hm3 %>%
      filter(metabolite_label %in% top_mets)
    
    # --- rank x-axis by score (highest score first) ---
    met_order <- hm3_plot %>%
      group_by(metabolite_label) %>%
      summarise(score_max = max(score), .groups = "drop") %>%
      arrange(desc(score_max)) %>%
      pull(metabolite_label)
    
    hm3_plot <- hm3_plot %>%
      mutate(
        metabolite_label = factor(metabolite_label, levels = met_order),
        module_protein = factor(module_protein)   # y-axis single row "3"
      )
    
    p_hm <- ggplot(hm3_plot, aes(x = metabolite_label, y = module_protein, fill = score)) +
      geom_tile() +
      labs(
        title = "Metabolite–module enrichment (significant pairs only)",
        subtitle = "Module 3 shown; tile value = -log10(q) for Fisher enrichment",
        x = "Metabolite (ranked by score)",
        y = "Protein module",
        fill = "score"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
    out_png <- file.path(figdir, "results_sec44_metabolite_module_enrichment_heatmap_ranked_module3.png")
    ggsave(out_png, p_hm, width = 20, height = 8, dpi = 300)
    cat("Saved heatmap:", out_png, "\n")
  }
}
