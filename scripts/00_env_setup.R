pkgs <- c(
  "tidyverse","data.table","here","fs","jsonlite",
  "VIM","splines","matrixStats",
  "igraph","tidygraph","ggraph","patchwork",
  "uwot","cluster","mixOmics","enrichplot"
)
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","org.Hs.eg.db","ReactomePA"), ask = FALSE, update = FALSE)
renv::snapshot(prompt = FALSE)
