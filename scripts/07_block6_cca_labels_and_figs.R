# --- 07_block6_cca_labels_and_figs.R -----------------------------------------
# Label metabolites (replace fid_* by names) for ridge/classic/scca top-loadings
# + Export figures for canonical correlations and k=1 top-loadings
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(jsonlite)
  library(scales)
})

# -- Load config / paths
cfg_rds  <- here("results","config.rds")
cfg_json <- here("results","config.json")
config <- if (file.exists(cfg_rds)) readRDS(cfg_rds) else
  jsonlite::read_json(cfg_json, simplifyVector = TRUE)

paths <- config$paths
stopifnot(!is.null(paths$results), !is.null(paths$figures), !is.null(paths$data_raw))
dir.create(paths$figures, showWarnings = FALSE, recursive = TRUE)

# -- Load metabolite info
info_path <- file.path(paths$data_raw, "metabolite_info.txt")
stopifnot(file.exists(info_path))
meta_info <- readr::read_tsv(info_path, col_names = FALSE, show_col_types = FALSE)
stopifnot(ncol(meta_info) >= 2)
names(meta_info)[1:2] <- c("name","metabolite_label")

# Helper: safe label join
label_metabolites <- function(df){
  df %>%
    left_join(meta_info, by = "name") %>%
    mutate(label = if_else(is.na(metabolite_label) | metabolite_label=="", name, metabolite_label)) %>%
    select(-metabolite_label)
}

# --------------------------
# 1) Label ALL top-loadings
# --------------------------
top_dir <- file.path(paths$results, "block6_toploadings")
stopifnot(dir.exists(top_dir))

# Find all metabolite top-loading files across methods
met_files <- list.files(top_dir, pattern = "^metabolites_(ridge|classic|scca)_k\\d+\\.csv$", full.names = TRUE)

for (fp in met_files){
  df <- readr::read_csv(fp, show_col_types = FALSE)
  # expect columns: component, side, method, name, weight, abs_w, rank
  if (!all(c("name","weight") %in% names(df))) next
  df_lab <- label_metabolites(df)
  out <- sub("\\.csv$", "_labeled.csv", fp)
  readr::write_csv(df_lab, out)
}

message(sprintf("Labeled %d metabolite top-loading files.", length(met_files)))

# --------------------------------------
# 2) Canonical correlation barplot (k<=5)
# --------------------------------------
summ_paths <- c(
  ridge   = file.path(paths$results, "block6_cca_summary_ridge.csv"),
  classic = file.path(paths$results, "block6_cca_summary_classic.csv"),
  scca    = file.path(paths$results, "block6_scca_summary.csv")
)

load_summ <- function(path, method_name){
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE) %>%
    transmute(method = method_name, k = as.integer(k), cancor = as.numeric(cancor))
}

cca_summ <- bind_rows(
  load_summ(summ_paths["classic"], "Classic CCA"),
  load_summ(summ_paths["ridge"],   "Ridge CCA"),
  load_summ(summ_paths["scca"],    "Sparse CCA")
) %>% filter(!is.na(k), k <= 5)

p_cc <- ggplot(cca_summ, aes(x = factor(k), y = cancor, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ method, ncol = 1) +
  scale_y_continuous("Canonical correlation", limits = c(0, 1), breaks = pretty_breaks()) +
  scale_x_discrete("Component (k)") +
  ggtitle("Canonical correlations across methods (complete cases)") +
  theme_minimal(base_size = 12)

ggsave(file.path(paths$figures, "block6_cancor_by_method.pdf"), p_cc, width = 6, height = 8, device = cairo_pdf)
ggsave(file.path(paths$figures, "block6_cancor_by_method.png"), p_cc, width = 6, height = 8, dpi = 300)

# ---------------------------------------------------------
# 3) k=1 top-loadings lollipop plots (proteins & metabolites)
#    for each method; metabolites use labeled *_labeled.csv
# ---------------------------------------------------------
make_lolli <- function(df, title, n = 20, flip = TRUE){
  d <- df %>%
    mutate(sign = if_else(weight >= 0, "positive", "negative")) %>%
    arrange(desc(abs_w)) %>% slice_head(n = n) %>%
    mutate(name_plot = if ("label" %in% names(.)) label else name) %>%
    mutate(name_plot = factor(name_plot, levels = rev(name_plot)))
  
  g <- ggplot(d, aes(x = name_plot, y = weight)) +
    geom_segment(aes(x = name_plot, xend = name_plot, y = 0, yend = weight)) +
    geom_point(size = 2) +
    labs(x = NULL, y = "Loading weight", title = title) +
    theme_minimal(base_size = 11)
  if (flip) g <- g + coord_flip()
  g
}

plot_top_for <- function(side = c("proteins","metabolites"), method_tag = c("ridge","classic","scca")){
  side <- match.arg(side); method_tag <- match.arg(method_tag)
  base <- file.path(paths$results, "block6_toploadings")
  fn <- if (side == "metabolites") {
    # Prefer the labeled file if present
    cand <- file.path(base, sprintf("metabolites_%s_k01_labeled.csv", method_tag))
    if (file.exists(cand)) cand else file.path(base, sprintf("metabolites_%s_k01.csv", method_tag))
  } else {
    file.path(base, sprintf("proteins_%s_k01.csv", method_tag))
  }
  if (!file.exists(fn)) return(NULL)
  df <- readr::read_csv(fn, show_col_types = FALSE)
  title <- sprintf("CCA k=1 top loadings — %s (%s)", if (side=="proteins") "Proteins" else "Metabolites", toupper(method_tag))
  make_lolli(df, title, n = 20, flip = TRUE)
}

plots <- list(
  prot_ridge = plot_top_for("proteins","ridge"),
  prot_class = plot_top_for("proteins","classic"),
  prot_scca  = plot_top_for("proteins","scca"),
  met_ridge  = plot_top_for("metabolites","ridge"),
  met_class  = plot_top_for("metabolites","classic"),
  met_scca   = plot_top_for("metabolites","scca")
)

# Save each plot if it exists
for (nm in names(plots)){
  gp <- plots[[nm]]
  if (inherits(gp, "gg")) {
    ggsave(file.path(paths$figures, paste0("block6_", nm, ".pdf")), gp, width = 7, height = 6, device = cairo_pdf)
    ggsave(file.path(paths$figures, paste0("block6_", nm, ".png")), gp, width = 7, height = 6, dpi = 300)
  }
}

# ---------------------------------------------------
# 4) Export a LaTeX-ready table snippet for summaries
# ---------------------------------------------------
tex_dir <- file.path(paths$results, "tex")
dir.create(tex_dir, showWarnings = FALSE, recursive = TRUE)

tab <- cca_summ %>%
  mutate(method = factor(method, levels = c("Classic CCA","Ridge CCA","Sparse CCA"))) %>%
  pivot_wider(names_from = k, values_from = cancor) %>%
  arrange(method)

# Write as a tiny LaTeX table (booktabs)
tex <- c(
  "\\begin{table}[h]",
  "\\centering",
  "\\caption{Canonical correlations (k=1--5) across methods on complete cases ($n=3642$). Ridge adds mild regularization; sparse CCA enforces sparsity for interpretability.}",
  "\\label{tab:cca-summary}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Method & k=1 & k=2 & k=3 & k=4 & k=5 \\\\",
  "\\midrule"
)
row_to_line <- function(m, v){
  sprintf("%s & %s \\\\", m, paste(sprintf("%.3f", v), collapse = " & "))
}
for (i in seq_len(nrow(tab))){
  vals <- as.numeric(tab[i, -1, drop = TRUE])
  tex <- c(tex, row_to_line(as.character(tab$method[i]), vals))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(tex_dir, "cca_summary_table.tex"))

message("Done: labels + figures + LaTeX table snippet exported.")
# ----------------------------------------------------------------------------- 
