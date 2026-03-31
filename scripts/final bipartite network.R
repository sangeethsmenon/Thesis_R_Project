# ============================================================
# Bipartite two-layer subgraph (balanced pos + neg)
# Proteins (top) = circles, Metabolites (bottom) = squares
# Edge colour by sign: Negative=red, Positive=green
# Straight lines, white background, no duplicate edges
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# ---------------------------
# INPUT
# ---------------------------
in_file  <- "results/block5_edges_filtered_metab_quantile.csv"
out_file <- "figures/block5_pm_bipartite_subgraph_pos_neg.png"

# display controls
n_proteins <- 3
n_mets_total <- 10              # total metabolites shown
top_edge_pool <- 5000           # bigger pool => more chance to find negatives

# try to balance metabolites across signs
# (if there are fewer negatives available, it will fall back to what exists)
prop_neg <- 0.35                # 35% metabolites chosen from negative edges
n_mets_neg <- max(1, round(n_mets_total * prop_neg))
n_mets_pos <- n_mets_total - n_mets_neg

# label controls
wrap_width <- 18
met_angle  <- 35

# ---------------------------
# READ + CLEAN
# ---------------------------
edges <- read_csv(in_file, show_col_types = FALSE) %>%
  mutate(
    protein = as.character(protein),
    metabolite_label = as.character(metabolite_label),
    r = as.numeric(r),
    abs_r = abs(r),
    sign = case_when(
      sign %in% c("pos","neg") ~ sign,
      r >= 0 ~ "pos",
      TRUE ~ "neg"
    )
  ) %>%
  filter(!is.na(protein), !is.na(metabolite_label), !is.na(r))

# Remove duplicate protein–metabolite pairs (keep strongest |r|)
edges <- edges %>%
  group_by(protein, metabolite_label) %>%
  slice_max(order_by = abs_r, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---------------------------
# PICK TOP PROTEINS (from large pool of strongest edges)
# ---------------------------
pool <- edges %>%
  arrange(desc(abs_r)) %>%
  slice_head(n = min(top_edge_pool, nrow(edges)))

prot_rank <- pool %>%
  count(protein, name = "deg") %>%
  arrange(desc(deg), protein)

top_prots <- prot_rank %>%
  slice_head(n = n_proteins) %>%
  pull(protein)

if (length(top_prots) < n_proteins) {
  message("Pool did not contain enough unique proteins; using global degree ranking.")
  top_prots <- edges %>%
    count(protein, name = "deg") %>%
    arrange(desc(deg), protein) %>%
    slice_head(n = n_proteins) %>%
    pull(protein)
}

# ---------------------------
# FILTER TO SELECTED PROTEINS
# ---------------------------
edges_focus <- edges %>%
  filter(protein %in% top_prots) %>%
  arrange(desc(abs_r))

# ---------------------------
# CHOOSE METABOLITES BALANCED BY SIGN
# ---------------------------
mets_pos <- edges_focus %>%
  filter(sign == "pos") %>%
  distinct(metabolite_label, .keep_all = TRUE) %>%
  slice_head(n = n_mets_pos) %>%
  pull(metabolite_label)

mets_neg <- edges_focus %>%
  filter(sign == "neg") %>%
  distinct(metabolite_label, .keep_all = TRUE) %>%
  slice_head(n = n_mets_neg) %>%
  pull(metabolite_label)

# If not enough negatives (common), fill remaining from positives (or vice versa)
mets <- unique(c(mets_pos, mets_neg))

if (length(mets) < n_mets_total) {
  need <- n_mets_total - length(mets)
  fill <- edges_focus %>%
    filter(!metabolite_label %in% mets) %>%
    distinct(metabolite_label, .keep_all = TRUE) %>%
    slice_head(n = need) %>%
    pull(metabolite_label)
  mets <- c(mets, fill)
}

# keep exact requested count
mets <- mets[seq_len(min(length(mets), n_mets_total))]

# final edges
sub_edges <- edges_focus %>%
  filter(metabolite_label %in% mets) %>%
  mutate(sign = factor(sign, levels = c("neg","pos")))

# ---------------------------
# LAYOUT (even spacing)
# ---------------------------
x_mets <- seq(1, length(mets) * 1.35, length.out = length(mets))
x_prot <- seq(min(x_mets), max(x_mets), length.out = length(top_prots))

nodes_mets <- tibble(name = mets, type = "Metabolite", x = x_mets, y = 0)
nodes_prot <- tibble(name = top_prots, type = "Protein", x = x_prot, y = 1)

nodes <- bind_rows(nodes_prot, nodes_mets) %>%
  mutate(label = ifelse(type == "Metabolite",
                        str_wrap(name, width = wrap_width),
                        name))

sub_edges_xy <- sub_edges %>%
  left_join(nodes %>% select(name, x, y), by = c("protein" = "name")) %>%
  rename(x_prot = x, y_prot = y) %>%
  left_join(nodes %>% select(name, x, y), by = c("metabolite_label" = "name")) %>%
  rename(x_met = x, y_met = y)

# ---------------------------
# PLOT
# ---------------------------
p <- ggplot() +
  geom_segment(
    data = sub_edges_xy,
    aes(x = x_prot, y = y_prot, xend = x_met, yend = y_met, colour = sign),
    linewidth = 0.9,
    alpha = 0.95
  ) +
  # metabolites squares
  geom_point(
    data = nodes %>% filter(type == "Metabolite"),
    aes(x = x, y = y),
    shape = 15, size = 3
  ) +
  # proteins circles
  geom_point(
    data = nodes %>% filter(type == "Protein"),
    aes(x = x, y = y),
    shape = 16, size = 4.2
  ) +
  geom_text(
    data = nodes %>% filter(type == "Protein"),
    aes(x = x, y = y, label = label),
    vjust = -0.9,
    fontface = "bold",
    size = 5
  ) +
  geom_text(
    data = nodes %>% filter(type == "Metabolite"),
    aes(x = x, y = y, label = label),
    angle = met_angle,
    hjust = 1,
    vjust = 1.1,
    size = 3.6
  ) +
  scale_colour_manual(
    values = c("neg" = "red", "pos" = "green"),
    labels = c("neg" = "Negative", "pos" = "Positive"),
    name = "Association sign"
  ) +
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("Metabolites", "Proteins"),
    limits = c(-0.70, 1.25)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.margin      = margin(20, 40, 100, 40)
  )

# save
ggsave(out_file, p, width = 18, height = 8, dpi = 300, bg = "white")
message("Saved: ", out_file)

p
