# minimal prerequisites
pkgs <- c("readr","dplyr","here")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

# coverage check
lab <- readr::read_csv(here::here("results","block5_edges_labeled_q_le_alpha_absr.csv"), show_col_types = FALSE)
cat("Unique metabolites in edges:", dplyr::n_distinct(lab$metabolite), "\n")
cat("Labeled (non-NA) names:", dplyr::n_distinct(lab$metabolite_label[!is.na(lab$metabolite_label)]), "\n")
cat("Unmapped fids (rows):", sum(is.na(lab$metabolite_label)), "\n")

lab %>%
  dplyr::filter(protein == "LDLR") %>%
  dplyr::arrange(dplyr::desc(abs_r)) %>%
  dplyr::select(protein, metabolite_label, r, q) %>%
  head(10) %>% print()

nodes_met <- readr::read_csv(here::here("results","block5_nodes_metabolites_labeled.csv"), show_col_types = FALSE)
cat("Metabolite nodes still showing as fid_*:", sum(grepl("^fid_", nodes_met$name)), "\n")
