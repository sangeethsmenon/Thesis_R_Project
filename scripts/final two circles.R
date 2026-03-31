# ============================================================
# Two-ring network + separate legend panel (NO OVERLAP)
# Output:
#   results/block5_two_ring_top100_FINAL_no_overlap.png
# ============================================================

library(tidyverse)
library(scales)
library(cowplot)

# ----------------------------
# USER SETTINGS
# ----------------------------
topN <- 100
n_proteins_inner <- 30
inner_radius <- 0.62
outer_radius <- 1.05
label_push <- 1.20

# Main figure size (network panel)
main_w <- 22
main_h <- 16

# Legend panel size
leg_h  <- 4

# Total output size
save_w <- main_w
save_h <- main_h + leg_h

# Big margins for network only (labels need room)
mar_top <- 320
mar_right <- 240
mar_bottom <- 280   # smaller now, because legends are not inside this panel
mar_left <- 170

set.seed(1)

# ----------------------------
# LOAD DATA
# ----------------------------
edges <- read_csv("results/block5b_pm_edges_with_protein_module.csv",
                  show_col_types = FALSE) %>%
  mutate(
    sign  = ifelse(r >= 0, "Positive", "Negative"),
    abs_r = abs(r)
  )

# ----------------------------
# SELECT EDGES (topN by |r|)
# ----------------------------
edges_top <- edges %>%
  arrange(desc(abs_r)) %>%
  slice(1:topN)

# ensure at least n_proteins_inner proteins appear
prot_in_top <- edges_top %>% distinct(protein) %>% nrow()
if (prot_in_top < n_proteins_inner) {
  current_prots <- unique(edges_top$protein)
  
  extra <- edges %>%
    arrange(desc(abs_r)) %>%
    anti_join(edges_top %>% select(protein, metabolite_label, r),
              by = c("protein","metabolite_label","r")) %>%
    filter(!(protein %in% current_prots)) %>%
    slice_head(n = 300)
  
  edges_top <- bind_rows(edges_top, extra) %>%
    distinct(protein, metabolite_label, r, .keep_all = TRUE)
}

message("Edges plotted: ", nrow(edges_top))
message("Proteins plotted: ", edges_top %>% distinct(protein) %>% nrow())
message("Metabolites plotted: ", edges_top %>% distinct(metabolite_label) %>% nrow())

# ----------------------------
# NODE LISTS
# ----------------------------
proteins <- edges_top %>% distinct(protein) %>%
  rename(name = protein) %>%
  mutate(type = "Protein")

metabolites <- edges_top %>% distinct(metabolite_label) %>%
  rename(name = metabolite_label) %>%
  mutate(type = "Metabolite")

nodes <- bind_rows(proteins, metabolites)

# ----------------------------
# DEGREE
# ----------------------------
deg <- edges_top %>%
  select(protein, metabolite_label) %>%
  pivot_longer(cols = everything(),
               names_to = "node_type",
               values_to = "node_name") %>%
  count(node_name, name = "degree")

nodes <- nodes %>%
  left_join(deg, by = c("name" = "node_name")) %>%
  mutate(degree = replace_na(degree, 0))

# ----------------------------
# CIRCULAR LAYOUT
# ----------------------------
nodes <- nodes %>%
  group_by(type) %>%
  arrange(desc(degree), name, .by_group = TRUE) %>%
  mutate(
    angle = 2 * pi * (row_number() - 1) / n(),
    radius = ifelse(type == "Protein", inner_radius, outer_radius),
    x = radius * cos(angle),
    y = radius * sin(angle),
    angle_deg = angle * 180 / pi,
    hjust = ifelse(angle_deg > 90 & angle_deg < 270, 1, 0),
    angle_text = ifelse(angle_deg > 90 & angle_deg < 270, angle_deg + 180, angle_deg)
  ) %>%
  ungroup()

# ----------------------------
# EDGE COORDINATES
# ----------------------------
edges_xy <- edges_top %>%
  left_join(nodes %>% filter(type == "Protein") %>%
              select(protein = name, x1 = x, y1 = y),
            by = "protein") %>%
  left_join(nodes %>% filter(type == "Metabolite") %>%
              select(metabolite_label = name, x2 = x, y2 = y),
            by = "metabolite_label")

# ----------------------------
# BASE PLOT (with legends)
# We'll reuse this to extract a clean legend panel.
# ----------------------------
base_plot <- ggplot() +
  geom_curve(
    data = edges_xy,
    aes(x = x1, y = y1, xend = x2, yend = y2, color = sign, size = abs_r),
    curvature = 0.25,
    alpha = 0.25
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y, shape = type, fill = degree),
    size = 3.0,
    stroke = 0.35,
    color = "black"
  ) +
  geom_text(
    data = nodes,
    aes(
      x = x * label_push,
      y = y * label_push,
      label = name,
      angle = angle_text,
      hjust = hjust
    ),
    size = 4.2
  ) +
  coord_equal(clip = "off") +
  scale_color_manual(values = c(Positive = "#2ca25f", Negative = "#de2d26")) +
  scale_size_continuous(range = c(0.25, 1.4), breaks = pretty_breaks(3)) +
  scale_shape_manual(values = c(Protein = 21, Metabolite = 22)) +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    color = "Association sign",
    size  = "|r|",
    shape = "Node type",
    fill  = "Node degree"
  ) +
  theme_void(base_size = 18) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(mar_top, mar_right, mar_bottom, mar_left),
    
    # Legends in the base plot (we will extract them)
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    
    legend.key.height = unit(0.9, "cm"),
    legend.key.width  = unit(1.4, "cm")
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    size  = guide_legend(order = 2, nrow = 1),
    shape = guide_legend(order = 3, nrow = 1),
    fill  = guide_colorbar(order = 4,
                           barwidth = unit(16, "cm"),
                           barheight = unit(0.9, "cm"))
  )

# ----------------------------
# 1) MAIN NETWORK PANEL (NO LEGENDS)
# ----------------------------
network_panel <- base_plot + theme(legend.position = "none")

# ----------------------------
# 2) LEGEND-ONLY PANEL
# (extract legend, then draw it on a blank canvas with padding)
# ----------------------------
leg <- cowplot::get_legend(base_plot)

legend_panel <- ggdraw() +
  draw_plot(leg, x = 0.02, y = 0.10, width = 0.96, height = 0.85) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# ----------------------------
# 3) COMBINE (VERTICAL)
# ----------------------------
final_plot <- plot_grid(
  network_panel,
  legend_panel,
  ncol = 1,
  rel_heights = c(main_h, leg_h)
)

# ----------------------------
# SAVE
# ----------------------------
dir.create("figures", showWarnings = FALSE)
out_fp <- "figures/block5_two_ring_top100_FINAL_no_overlap.png"

ggsave(
  filename = out_fp,
  plot = final_plot,
  width = save_w,
  height = save_h,
  dpi = 300,
  bg = "white"
)

message("Saved final: ", out_fp)
final_plot
