library(tidyverse)

# Load enrichment results
met_enr <- read_csv(
  "results/block5b_metabolite_module_enrichment.csv",
  show_col_types = FALSE
)

# Parameters
MODULE_ID <- 3
N_TOP <- 200   # change to 10 / 20 if you want

# Create compact table for Module 3
module3_table <- met_enr %>%
  filter(
    module_protein == MODULE_ID,
    q_fisher <= 0.05
  ) %>%
  mutate(
    score = -log10(q_fisher)
  ) %>%
  arrange(q_fisher) %>%          # strongest enrichment first
  select(
    metabolite = metabolite_label,
    q_value = q_fisher,
    score
  ) %>%
  slice_head(n = N_TOP)

# View table
print(module3_table)

# Save for LaTeX / thesis
write_csv(
  module3_table,
  "results/results_sec44_module3_top_metabolites.csv"
)
