packages <- c("tidyverse","here","ggrepel","jsonlite")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# config + paths
cfg_rds  <- here::here("results","config.rds")
cfg_json <- here::here("results","config.json")
config <- if (file.exists(cfg_rds)) readRDS(cfg_rds) else jsonlite::read_json(cfg_json, simplifyVector = TRUE)

res_dir <- config$paths$results
fig_dir <- config$paths$figures
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Inputs produced by 06_block5_associations.R + 06b_label script
edges <- readr::read_csv(file.path(res_dir,"block5_edges_filtered_pruned_labeled.csv"), show_col_types = FALSE)
summ  <- readr::read_csv(file.path(res_dir,"block5_assoc_summary.csv"), show_col_types = FALSE)

#############################################################


p1 <- edges %>%
  mutate(minus_log10q = -log10(q),
         sign = if_else(r >= 0, "positive","negative")) %>%
  ggplot(aes(x = abs(r), y = minus_log10q, color = sign)) +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(title = "Block 5: Effect size vs significance",
       x = "|Pearson r|", y = "-log10(FDR q)") +
  theme_minimal()
ggsave(file.path(fig_dir,"block5_volcano_absr_vs_logq.png"), p1, width=7, height=5, dpi=300)


##########################################################
topP <- edges %>% count(protein, sort = TRUE) %>% slice_head(n=20)
topM <- edges %>% count(metabolite_label, sort = TRUE) %>% slice_head(n=20)

gp <- ggplot(topP, aes(x = reorder(protein, n), y = n)) +
  geom_col() + coord_flip() +
  labs(title="Top proteins by degree", x="", y="Edges kept") +
  theme_minimal()

gm <- ggplot(topM, aes(x = reorder(metabolite_label, n), y = n)) +
  geom_col() + coord_flip() +
  labs(title="Top metabolites by degree", x="", y="Edges kept") +
  theme_minimal()

ggsave(file.path(fig_dir,"block5_degree_proteins_top20.png"), gp, width=6, height=6, dpi=300)
ggsave(file.path(fig_dir,"block5_degree_metabolites_top20.png"), gm, width=6, height=6, dpi=300)

###############################################################


topP30 <- edges %>% group_by(protein) %>% summarise(max_abs = max(abs(r))) %>%
  slice_max(max_abs, n=30) %>% pull(protein)
topM20 <- edges %>% group_by(metabolite_label) %>% summarise(max_abs = max(abs(r))) %>%
  slice_max(max_abs, n=20) %>% pull(metabolite_label)

hm_df <- edges %>%
  filter(protein %in% topP30, metabolite_label %in% topM20)

p_hm <- ggplot(hm_df, aes(x = metabolite_label, y = protein, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(limits=c(-1,1)) +
  labs(title="Top edges heatmap (r)", x="", y="") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(file.path(fig_dir,"block5_heatmap_top30x20.png"), p_hm, width=9, height=8, dpi=300)
##############################################################


# pick the single strongest |r|
best <- edges %>% arrange(desc(abs(r))) %>% slice(1)

b4 <- readRDS(file.path(config$paths$data_processed, "block4_residuals_objects.rds"))
Mx <- b4$metabolome_resid_std; Px <- b4$proteome_resid_std
met_names <- b4$metabolome_features; prot_names <- b4$proteome_features

ixP <- match(best$protein, prot_names)
ixM <- match(best$metabolite, met_names)

df_sc <- tibble(
  prot = Px[, ixP],
  metab = Mx[, ixM]
) %>% drop_na()

p_sc <- ggplot(df_sc, aes(x = prot, y = metab)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = paste0("Example pair: ", best$protein, " vs ", best$metabolite_label,
                      " (r=", round(best$r,3), ", q=", signif(best$q,3), ")"),
       x = paste0(best$protein," residual"), y = paste0(best$metabolite_label," residual")) +
  theme_minimal()

ggsave(file.path(fig_dir,"block5_scatter_strongest_pair.png"), p_sc, width=6.5, height=5, dpi=300)
############################################################

p_agree <- edges %>%
  ggplot(aes(x = r, y = r_s)) +
  geom_abline(linetype="dashed") +
  geom_point(alpha = 0.4, size = 0.7) +
  coord_equal(xlim=c(-1,1), ylim=c(-1,1)) +
  labs(title = "Pearson vs Spearman r (kept edges)", x="Pearson r", y="Spearman r") +
  theme_minimal()
ggsave(file.path(fig_dir,"block5_pearson_vs_spearman.png"), p_agree, width=6, height=6, dpi=300)



