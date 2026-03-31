# ------------------------------------------------------------
# scripts/06b_block5b_protein_first_modules.R
# Protein-first pipeline (Block 5b, Part 1):
#  1) Load Block 4 residualized proteome
#  2) Build protein–protein correlation network + Louvain modules
#  3) Load Block 5 protein–metabolite edges (PM)
#  4) Map metabolite fid -> metabolite_name (metabolite_info.txt)
#  5) Attach protein modules to PM edges
#  6) Metabolite -> protein-module enrichment (Fisher)
#
# Outputs (in results/):
#  - block5b_pp_edges.csv
#  - block5b_protein_modules.csv
#  - block5b_pm_edges_with_protein_module.csv
#  - block5b_metabolite_module_enrichment.csv
# ------------------------------------------------------------

packages <- c("tidyverse","here","igraph","jsonlite")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# ---------------------------
# 1) Load config
# ---------------------------
if (!exists("config")) {
  cfg_rds  <- here::here("results","config.rds")
  cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) config <- readRDS(cfg_rds)
  else if (file.exists(cfg_json)) config <- jsonlite::read_json(cfg_json, simplifyVector = TRUE)
  else stop("Config not found. Run Block 0.")
}

resdir <- config$paths$results
if (is.null(resdir)) resdir <- here::here("results")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)

datadir <- config$paths$data_processed
if (is.null(datadir)) datadir <- here::here("data","processed")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 2) Load Block 4 residualized data
# ---------------------------
b4_path <- file.path(datadir, "block4_residuals_objects.rds")
if (!file.exists(b4_path)) stop("Block 4 outputs not found: ", b4_path)
b4 <- readRDS(b4_path)

P <- as.matrix(b4$proteome_resid_std)     # n x p (residualized; may contain NAs)
prot_names <- b4$proteome_features
stopifnot(ncol(P) == length(prot_names))
colnames(P) <- prot_names

cat("Loaded residualized proteome:", nrow(P), "samples x", ncol(P), "proteins\n")

# ---------------------------
# 3) Protein–protein correlations -> edge list -> filtered PP network
# ---------------------------
Rpp <- cor(P, use = "pairwise.complete.obs", method = "pearson")
diag(Rpp) <- NA_real_

# IMPORTANT: use base::as.table to avoid weird masks
pp_edges <- as.data.frame(base::as.table(Rpp), stringsAsFactors = FALSE)
colnames(pp_edges) <- c("protein1","protein2","r")

pp_edges <- pp_edges %>%
  filter(!is.na(r), protein1 < protein2) %>%
  mutate(abs_r = abs(r))

# Filtering thresholds (tune as needed)
q_pp <- 0.01
r_pp <- 0.30

# Approx p-values using n_eff; (screening for modules; okay as first pass)
n_eff <- nrow(P)
tval <- pp_edges$r * sqrt((n_eff - 2) / pmax(1e-12, (1 - pp_edges$r^2)))
pp_edges$p <- 2 * pt(-abs(tval), df = n_eff - 2)
pp_edges$q <- p.adjust(pp_edges$p, method = "BH")

pp_edges_f <- pp_edges %>% filter(q <= q_pp, abs_r >= r_pp)
if (nrow(pp_edges_f) == 0) stop("Protein–protein filtering gave 0 edges. Try r_pp=0.25 or q_pp=0.05.")

readr::write_csv(pp_edges_f, file.path(resdir, "block5b_pp_edges.csv"))
cat("PP edges kept:", nrow(pp_edges_f), "\n")

# ---------------------------
# 4) Build protein-only network + Louvain modules
# ---------------------------
g_pp <- igraph::graph_from_data_frame(
  d = pp_edges_f %>% transmute(from = protein1, to = protein2, weight = abs_r),
  directed = FALSE
)

set.seed(42)
louv_pp <- igraph::cluster_louvain(g_pp, weights = E(g_pp)$weight)

prot_modules <- tibble(
  protein = V(g_pp)$name,
  module_protein = as.integer(igraph::membership(louv_pp)),
  degree_pp = igraph::degree(g_pp)
)

readr::write_csv(prot_modules, file.path(resdir, "block5b_protein_modules.csv"))
cat("Protein modules:", dplyr::n_distinct(prot_modules$module_protein), "\n")

# ---------------------------
# 5) Load PM edges from your original Block 5 outputs
# ---------------------------
#pm_path <- file.path(resdir, "block5_edges_filtered_pruned.csv")
pm_path <- file.path(resdir, "block5_edges_filtered_metab_quantile.csv")
if (!file.exists(pm_path)) pm_path <- file.path(resdir, "block5_edges_filtered.csv")
if (!file.exists(pm_path)) stop("Could not find block5_edges_filtered(_pruned).csv in results/")

pm <- readr::read_csv(pm_path, show_col_types = FALSE)
stopifnot(all(c("protein","metabolite","r","q") %in% names(pm)))
if (!"abs_r" %in% names(pm)) pm <- pm %>% mutate(abs_r = abs(r))

# ---------------------------
# 6) Map metabolite fid -> metabolite_name using metabolite_info.txt
#    (your format: fid <tab> metabolite_name)
# ---------------------------
data_raw <- config$paths$data_raw
if (is.null(data_raw)) data_raw <- here::here("data","raw")

map_fp <- file.path(data_raw, "metabolite_info.txt")
if (!file.exists(map_fp)) stop("metabolite_info.txt not found at: ", map_fp)

metab_map <- readr::read_tsv(map_fp, show_col_types = FALSE) %>%
  transmute(
    metabolite = as.character(fid),
    metabolite_label = as.character(metabolite_name)
  ) %>%
  distinct(metabolite, .keep_all = TRUE)

pm <- pm %>%
  left_join(metab_map, by = "metabolite") %>%
  mutate(metabolite_label = coalesce(metabolite_label, metabolite))

cat("Metabolite mapping loaded:", nrow(metab_map), "rows\n")

# ---------------------------
# 7) Attach protein-module to each PM edge
# ---------------------------
pm_annot <- pm %>%
  inner_join(prot_modules %>% select(protein, module_protein), by = "protein")

readr::write_csv(pm_annot, file.path(resdir, "block5b_pm_edges_with_protein_module.csv"))
cat("PM edges that map to PP-module proteins:", nrow(pm_annot), "\n")

# ---------------------------
# 8) Metabolite -> module enrichment (Fisher)
#    Universe = proteins present in PP network (have a module)
# ---------------------------
universe <- unique(prot_modules$protein)
N <- length(universe)

module_sizes <- prot_modules %>% count(module_protein, name = "K")

metab_totals <- pm_annot %>%
  distinct(metabolite, metabolite_label, protein) %>%
  count(metabolite, metabolite_label, name = "n")

metab_module_counts <- pm_annot %>%
  distinct(metabolite, metabolite_label, protein, module_protein) %>%
  count(metabolite, metabolite_label, module_protein, name = "k") %>%
  left_join(module_sizes, by = "module_protein") %>%
  left_join(metab_totals, by = c("metabolite","metabolite_label")) %>%
  mutate(
    p_fisher = purrr::pmap_dbl(list(k, K, n), function(k, K, n) {
      a <- k
      b <- K - k
      c <- n - k
      d <- (N - K) - (n - k)
      if (any(c(a,b,c,d) < 0)) return(NA_real_)
      fisher.test(matrix(c(a,b,c,d), nrow = 2), alternative = "greater")$p.value
    }),
    q_fisher = p.adjust(p_fisher, method = "BH")
  ) %>%
  arrange(q_fisher, desc(k))

readr::write_csv(metab_module_counts, file.path(resdir, "block5b_metabolite_module_enrichment.csv"))

cat(
  "Block 5b (Part 1) complete.\n",
  "PP nodes:", igraph::gorder(g_pp), " PP edges:", igraph::gsize(g_pp), "\n",
  "Saved:\n",
  " - results/block5b_pp_edges.csv\n",
  " - results/block5b_protein_modules.csv\n",
  " - results/block5b_pm_edges_with_protein_module.csv\n",
  " - results/block5b_metabolite_module_enrichment.csv\n"
)
