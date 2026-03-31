packages <- c("tidyverse","here","fs")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

res <- here::here("results")
fig <- here::here("figures")
out <- here::here("deliverables","v1")
dir_create(out, recurse = TRUE)

read_if <- function(path) if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else NULL
copy_if <- function(path, dest = out) if (file.exists(path)) file_copy(path, dest, overwrite = TRUE)

# --- Collect summaries (exists-checks baked in)
s5   <- read_if(file.path(res,"block5_assoc_summary.csv"))
s5a  <- read_if(file.path(res,"block5a_unadjusted_assoc_summary.csv"))
s5c  <- read_if(file.path(res,"block5c_sex_stratified_summary.csv"))
s6   <- read_if(file.path(res,"block6_cca_summary.csv"))
s7n  <- read_if(file.path(res,"block7_network_summary.csv"))
s7m  <- read_if(file.path(res,"block7_module_summary.csv"))

# Headline numbers (NA-safe)
get1 <- function(df, col) if (!is.null(df) && col %in% names(df)) df[[col]][1] else NA

n_samp     <- get1(s5,  "n_samples")
n_prot     <- get1(s5,  "n_proteins")
n_metab    <- get1(s5,  "n_metabolites")
n_tests    <- get1(s5,  "n_tests")
n_edges_q  <- get1(s5,  "pearson_sig")
n_edges_qes<- get1(s5,  "pearson_sig_es")
n_edges_s  <- get1(s5,  "spearman_sig")
overlap    <- get1(s5,  "overlap_sig")
cca_k1     <- get1(s6,  "cancor")
cca_p1     <- get1(s6,  "p_perm")
net_edges  <- get1(s7n, "n_edges")
net_nodes  <- get1(s7n, "n_nodes")
net_mods   <- get1(s7n, "n_modules")
net_modq   <- get1(s7n, "modularity")
n_F_sig    <- get1(s5c, "F_sig")
n_M_sig    <- get1(s5c, "M_sig")
n_diff_sig <- get1(s5c, "diff_sig")

# Write a small key-numbers CSV
tibble::tibble(
  metric = c("samples","proteins","metabolites","pairwise_tests",
             "sig_edges_q<=0.05","sig_edges_q<=0.05_|r|>=0.1","sig_edges_spearman",
             "pearson_spearman_overlap","CCA_k1_correlation","CCA_k1_perm_p",
             "network_edges","network_nodes","network_modules","network_modularity",
             "F_sig_edges","M_sig_edges","sex_diff_edges"),
  value  = c(n_samp, n_prot, n_metab, n_tests,
             n_edges_q, n_edges_qes, n_edges_s, overlap, cca_k1, cca_p1,
             net_edges, net_nodes, net_mods, net_modq, n_F_sig, n_M_sig, n_diff_sig)
) %>% readr::write_csv(file.path(out,"key_numbers.csv"))

# Compose a short README-style summary
md <- c(
  "# Results snapshot (Blocks 0–7)\n",
  "## Data\n",
  sprintf("- Samples used: **%s**; Proteins: **%s**; Metabolites: **%s**.\n", n_samp, n_prot, n_metab),
  "## Pairwise protein–metabolite associations\n",
  sprintf("- Tests: **%s**. Significant (Pearson, FDR≤0.05): **%s** (of which |r|≥0.1: **%s**). Spearman sig: **%s**. Overlap: **%s**.\n",
          n_tests, n_edges_q, n_edges_qes, n_edges_s, overlap),
  "## Sex-stratified associations (age-adjusted)\n",
  sprintf("- Female sig edges: **%s**; Male sig edges: **%s**; Edges with significant F–M difference: **%s**.\n",
          n_F_sig, n_M_sig, n_diff_sig),
  "## Canonical correlation analysis (age/sex residualized)\n",
  sprintf("- First canonical correlation: **%.3f** (perm p = **%.4g**). Subsequent components also significant.\n",
          cca_k1 %||% NA_real_, cca_p1 %||% NA_real_),
  "## Network (FDR≤0.01, |r|≥0.20)\n",
  sprintf("- Edges: **%s**; Nodes: **%s**; Modules: **%s**; Modularity: **%.3f**.\n",
          net_edges, net_nodes, net_mods, net_modq),
  "## Files included in this bundle\n",
  "- `key_numbers.csv` (this table)\n",
  "- `block5_assoc_summary.csv`, `block5a_unadjusted_assoc_summary.csv`\n",
  "- `block5c_sex_stratified_summary.csv`, `block5c_edges_top20k_by_pdiff.csv`\n",
  "- `block6_cca_summary.csv`, top loadings under `results/block6_toploadings/`\n",
  "- `block7_network_summary.csv`, `block7_module_summary.csv`\n",
  "- Optional figures: main network PNG/SVG if present\n"
)
writeLines(md, con = file.path(out, "RESULTS_SNAPSHOT.md"))

# Copy key CSVs (when present)
to_copy <- c(
  "block5_assoc_summary.csv",
  "block5a_unadjusted_assoc_summary.csv",
  "block5c_sex_stratified_summary.csv",
  "block5c_edges_top20k_by_pdiff.csv",
  "block6_cca_summary.csv",
  "block7_network_summary.csv",
  "block7_module_summary.csv",
  "block7_enrichment_go_all.csv",
  "block7_enrichment_kegg_all.csv"
)
invisible(lapply(file.path(res, to_copy), copy_if))

# Try copying a couple of figures if they exist
invisible(lapply(c("block7_network_all.png","block7_network_by_module.png","block8_module_sizes.png"),
                 function(f) copy_if(file.path(fig, f))))

cat("Packaged deliverables to ", out, "\n", sep = "")
