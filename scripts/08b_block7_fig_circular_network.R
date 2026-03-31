# ============================================================
# 08b_block7_fig_circular_network.R  (UPDATED)
# Thesis figures for Protein–Metabolite (PM) bipartite network:
#   (A) Edge-weight histogram (choose cutoff)
#   (B) Ring / inner–outer overview (proteins inner, metabolites outer)
#   (C1) Two-layer bipartite subgraph: top-|r| edges overall (readable)
#   (C2) Two-layer bipartite subgraph: diversified top edges (spreads across proteins)
#
# Key user-requested fixes:
# - Ring plots: reduce gap between inner and outer circles (tracks much closer)
# - Make labels larger & more legible (bigger canvas, smaller track heights, smaller offsets)
# - Bipartite: reduce edges (e.g., 10) and ensure labels don’t cross border (more margin + wrap + angle)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(scales)
  library(rlang)
})

# ---------------------------
# 0) Resolve paths (config + here fallback)
# ---------------------------
if (!exists("config")) {
  config <- list(paths = list(results = NULL, figures = NULL, data_raw = NULL))
}

resdir <- config$paths$results; if (is.null(resdir)) resdir <- here::here("results")
figdir <- config$paths$figures; if (is.null(figdir)) figdir <- here::here("figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

message("resdir = ", resdir)
message("figdir = ", figdir)

# ---------------------------
# 1) Input: Block 5 pruned PM edge table
# ---------------------------
pm_fp <- file.path(resdir, "block5_edges_filtered_metab_quantile_pruned.csv")
if (!file.exists(pm_fp)) {
  alt_fp <- file.path(resdir, "block5_edges_filtered_metab_quantile.csv")
  if (!file.exists(alt_fp)) {
    stop("PM edge file not found.\nTried:\n  ", pm_fp, "\n  ", alt_fp,
         "\n\nCheck your Block 5 outputs in: ", resdir)
  }
  pm_fp <- alt_fp
}

pm <- readr::read_csv(pm_fp, show_col_types = FALSE)
message("Loaded PM edges: ", nrow(pm), " rows from ", basename(pm_fp))

# ---------------------------
# 2) Standardize column names robustly
# ---------------------------
rename_if_present <- function(df, old, new) {
  if (old %in% names(df) && !(new %in% names(df))) dplyr::rename(df, !!new := !!rlang::sym(old)) else df
}

pm <- pm %>%
  rename_if_present("prot", "protein") %>%
  rename_if_present("met", "metabolite") %>%
  rename_if_present("met_name", "metabolite_label") %>%
  rename_if_present("r_pearson", "r") %>%
  rename_if_present("cor", "r") %>%
  mutate(
    protein = as.character(protein),
    metabolite = as.character(metabolite),
    metabolite_label = if ("metabolite_label" %in% names(.)) as.character(metabolite_label) else metabolite,
    r = as.numeric(r),
    abs_r = if ("abs_r" %in% names(.)) as.numeric(abs_r) else abs(r),
    sign = ifelse(is.na(r), NA_character_, ifelse(r >= 0, "pos", "neg")),
    met_display = dplyr::coalesce(metabolite_label, metabolite)
  ) %>%
  filter(!is.na(r), !is.na(abs_r), !is.na(sign))

needed <- c("protein", "metabolite", "metabolite_label", "met_display", "r", "abs_r", "sign")
missing_cols <- setdiff(needed, names(pm))
if (length(missing_cols) > 0) {
  stop("PM table is missing required columns: ", paste(missing_cols, collapse = ", "),
       "\nColumns available: ", paste(names(pm), collapse = ", "))
}

# ---------------------------
# 3) Figure A: Histogram of |r|
# ---------------------------
hist_fp <- file.path(figdir, "block5_pm_edgeweight_hist_absr.png")

p_hist <- ggplot(pm, aes(x = abs_r)) +
  geom_histogram(bins = 50) +
  labs(
    x = expression("|r| (absolute Pearson correlation)"),
    y = "Number of protein–metabolite edges",
    title = "Distribution of protein–metabolite association strengths"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

ggsave(hist_fp, p_hist, width = 9.0, height = 6.5, dpi = 300)
message("Wrote: ", hist_fp)

# ---------------------------
# 4) Choose a cutoff for a “thresholded” view (optional)
# ---------------------------
cut_strategy <- "quantile"  # "quantile" or "hard"
cut_quantile <- 0.95        # top 5% strongest edges by |r|
cut_hard <- 0.30            # used only if cut_strategy == "hard"

cut_absr <- if (cut_strategy == "hard") {
  cut_hard
} else {
  as.numeric(stats::quantile(pm$abs_r, probs = cut_quantile, na.rm = TRUE))
}
message("Cutoff chosen: |r| >= ", round(cut_absr, 3),
        ifelse(cut_strategy == "hard", " (hard)", paste0(" (quantile=", cut_quantile, ")")))

pm_cut <- pm %>% filter(abs_r >= cut_absr)

# ============================================================
# 5) Bipartite two-layer plots (readable labels, no clipping)
# ============================================================

# ---- USER TUNABLE SETTINGS (BIPARTITE) ----
top_edges_for_twolayer_main <- 10     # Plot 1: top edges overall (reduce crowding)
top_edges_for_twolayer_div  <- 10     # Plot 2: diversified across proteins
wrap_width_bipartite        <- 22     # wrap metabolite labels
met_label_angle             <- 35     # rotate metabolite labels
bipartite_width_inches      <- 16
bipartite_height_inches     <- 11

# Helper: build a two-layer plot from an edge table
make_two_layer_plot <- function(pm_edges, out_fp, subtitle_text,
                                wrap_width = 22,
                                angle = 35,
                                width = 16,
                                height = 11) {
  pm_edges <- pm_edges %>%
    mutate(
      met_display_wrapped = stringr::str_wrap(met_display, width = wrap_width)
    )

  deg_prot <- pm_edges %>% count(protein, name = "deg_protein")
  deg_met  <- pm_edges %>% count(met_display_wrapped, name = "deg_metabolite")

  prot_order <- deg_prot %>% arrange(desc(deg_protein), protein) %>% pull(protein)
  met_order  <- deg_met  %>% arrange(desc(deg_metabolite), met_display_wrapped) %>% pull(met_display_wrapped)

  prot_pos <- tibble(protein = prot_order) %>%
    mutate(x = row_number(), y = 1, label = protein)

  met_pos <- tibble(met_display_wrapped = met_order) %>%
    mutate(x = row_number(), y = 0, label = met_display_wrapped)

  edges_xy <- pm_edges %>%
    inner_join(prot_pos, by = "protein") %>%
    rename(x1 = x, y1 = y, prot_label = label) %>%
    inner_join(met_pos, by = "met_display_wrapped") %>%
    rename(x2 = x, y2 = y, met_label = label)

  p_two <- ggplot() +
    geom_segment(
      data = edges_xy,
      aes(x = x1, y = y1, xend = x2, yend = y2, colour = sign, alpha = abs_r),
      linewidth = 1.1
    ) +
    geom_point(data = prot_pos, aes(x = x, y = y), size = 4.2, colour = "red") +
    geom_point(data = met_pos,  aes(x = x, y = y), size = 4.2, shape = 15, colour = "blue") +
    geom_text(
      data = prot_pos,
      aes(x = x, y = y + 0.10, label = label),
      size = 5.0
    ) +
    geom_text(
      data = met_pos,
      aes(x = x, y = y - 0.12, label = label),
      size = 4.3,
      angle = angle,
      hjust = 1,
      vjust = 1.5
    ) +
    scale_colour_manual(values = c(pos = "darkgreen", neg = "red3")) +
    scale_alpha(range = c(0.35, 0.95), guide = "none") +
    labs(
      x = NULL, y = NULL,
      title = "Protein–metabolite bipartite subgraph (two-layer layout)",
      subtitle = subtitle_text,
      colour = "sign"
    ) +
    coord_cartesian(ylim = c(-0.60, 1.30), clip = "off") +
    theme_bw(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 24),
      plot.subtitle = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      # Large bottom margin ensures rotated labels never cross the border
      plot.margin = margin(t = 14, r = 24, b = 170, l = 24)
    )

  ggsave(out_fp, p_two, width = width, height = height, dpi = 300)
  message("Wrote: ", out_fp)
}

# ---- Plot 1: strongest edges overall ----
pm_main <- pm %>%
  arrange(desc(abs_r)) %>%
  slice_head(n = top_edges_for_twolayer_main)

two_layer_fp_main <- file.path(figdir, sprintf("block5_pm_bipartite_twolayer_top%d.png", top_edges_for_twolayer_main))
make_two_layer_plot(
  pm_edges = pm_main,
  out_fp = two_layer_fp_main,
  subtitle_text = paste0("Top ", top_edges_for_twolayer_main, " edges by |r| (overall); proteins (red) vs metabolites (blue); edges coloured by sign"),
  wrap_width = wrap_width_bipartite,
  angle = met_label_angle,
  width = bipartite_width_inches,
  height = bipartite_height_inches
)

# ---- Plot 2: diversified top edges across proteins ----
# Logic:
# 1) take the top 1 edge per protein (by |r|)
# 2) rank those edges; take as many as needed
# 3) if still short, fill with remaining strongest edges overall (excluding already picked)
get_diversified_edges <- function(pm_df, n_total) {
  per_prot <- pm_df %>%
    group_by(protein) %>%
    slice_max(order_by = abs_r, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(abs_r))

  picked <- per_prot %>% slice_head(n = min(n_total, nrow(per_prot)))

  if (nrow(picked) >= n_total) return(picked)

  remaining <- pm_df %>%
    anti_join(picked %>% select(protein, metabolite), by = c("protein", "metabolite")) %>%
    arrange(desc(abs_r)) %>%
    slice_head(n = n_total - nrow(picked))

  bind_rows(picked, remaining)
}

pm_div <- get_diversified_edges(pm, top_edges_for_twolayer_div)

two_layer_fp_div <- file.path(figdir, sprintf("block5_pm_bipartite_twolayer_diverse_top%d.png", top_edges_for_twolayer_div))
make_two_layer_plot(
  pm_edges = pm_div,
  out_fp = two_layer_fp_div,
  subtitle_text = paste0("Diversified set (n=", top_edges_for_twolayer_div, "): spreads top edges across proteins, then fills by |r|; proteins (red) vs metabolites (blue); edges coloured by sign"),
  wrap_width = wrap_width_bipartite,
  angle = met_label_angle,
  width = bipartite_width_inches,
  height = bipartite_height_inches
)

# ============================================================
# 6) Ring plots (circlize) — metabolites distributed around 360°
# ============================================================

have_circlize <- requireNamespace("circlize", quietly = TRUE)
if (!have_circlize) {
  message("Package 'circlize' not installed -> skipping ring plots.\nInstall once with: install.packages('circlize')")
} else {
  suppressPackageStartupMessages(library(circlize))
  
  # ---- USER TUNABLE SETTINGS (RING) ----
  ring_png_wh        <- 5200
  ring_png_res       <- 300
  
  # Make OUTER (metabolite) label track bigger and INNER (protein) label track smaller
  ring_track_h_outer <- 0.10   # bigger outer track = more room for metabolite labels
  ring_track_h_inner <- 0.04   # smaller inner track = effectively "smaller inner radius"
  
  ring_track_margin  <- c(0.001, 0.001)
  ring_outer_label_offset <- 0.08
  ring_inner_label_offset <- 0.05
  
  # Sector spacing
  start_degree <- 90
  gap_small    <- 0.5   # small uniform gaps between all sectors (no big gaps)
  
  # Interleave ordering so metabolites spread around the full circle
  interleave <- function(a, b) {
    na <- length(a); nb <- length(b)
    n  <- max(na, nb)
    a2 <- c(a, rep(NA_character_, n - na))
    b2 <- c(b, rep(NA_character_, n - nb))
    out <- as.vector(rbind(a2, b2))
    out[!is.na(out)]
  }
  
  make_ring_plot_inner_outer <- function(df_edges, out_fp, title_text,
                                         cex_prot = 0.65,
                                         cex_met  = 0.55,
                                         wrap_met = 18,
                                         wrap_prot = 999) {
    
    df_edges <- df_edges %>%
      mutate(
        from  = as.character(from),
        to    = as.character(to),
        abs_r = as.numeric(abs_r),
        sign  = as.character(sign)
      ) %>%
      filter(!is.na(from), !is.na(to), !is.na(abs_r), !is.na(sign))
    
    if (nrow(df_edges) == 0) {
      message("No edges for: ", out_fp)
      return(invisible(NULL))
    }
    
    proteins    <- unique(df_edges$from)
    metabolites <- unique(df_edges$to)
    
    prot_lab_map <- setNames(stringr::str_wrap(proteins, width = wrap_prot), proteins)
    met_lab_map  <- setNames(stringr::str_wrap(metabolites, width = wrap_met), metabolites)
    
    # KEY CHANGE: interleave to distribute metabolites around 360°
    order_levels <- interleave(metabolites, proteins)
    
    # small uniform gaps only (removes "metabolites only on one side" effect)
    gap_after <- rep(gap_small, length(order_levels))
    
    link_cols <- ifelse(df_edges$sign == "pos", "darkgreen", "red3")
    link_lwd  <- scales::rescale(df_edges$abs_r, to = c(0.7, 2.8))
    
    circos.clear()
    circos.par(
      start.degree = start_degree,
      gap.after = gap_after,
      points.overflow.warning = FALSE,
      track.margin = ring_track_margin,
      cell.padding = c(0, 0, 0, 0)
    )
    
    png(out_fp, width = ring_png_wh, height = ring_png_wh, res = ring_png_res)
    
    chordDiagram(
      x = df_edges %>% select(from, to, abs_r),
      order = order_levels,
      annotationTrack = c("grid"),
      preAllocateTracks = list(
        list(track.height = ring_track_h_outer),  # track 1 = outer labels (metabolites)
        list(track.height = ring_track_h_inner)   # track 2 = inner labels (proteins)
      ),
      transparency = 0.70,
      col = link_cols,
      link.lwd = link_lwd
    )
    
    # OUTER: metabolites
    circos.trackPlotRegion(
      track.index = 1, bg.border = NA,
      panel.fun = function(x, y) {
        sector <- get.cell.meta.data("sector.index")
        if (!(sector %in% metabolites)) return(NULL)
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(
          x = mean(xlim),
          y = ylim[2] + ring_outer_label_offset,
          labels = met_lab_map[[sector]],
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = cex_met
        )
      }
    )
    
    # INNER: proteins
    circos.trackPlotRegion(
      track.index = 2, bg.border = NA,
      panel.fun = function(x, y) {
        sector <- get.cell.meta.data("sector.index")
        if (!(sector %in% proteins)) return(NULL)
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(
          x = mean(xlim),
          y = ylim[1] - ring_inner_label_offset,
          labels = prot_lab_map[[sector]],
          facing = "reverse.clockwise",
          niceFacing = TRUE,
          adj = c(1, 0.5),
          cex = cex_prot
        )
      }
    )
    
    title(main = title_text, cex.main = 1.35)
    dev.off()
    circos.clear()
    message("Wrote: ", out_fp)
  }
  
  # (B1) Overview ring (top N edges by |r|)
  top_edges_for_ring <- 100
  
  pm_ring <- pm %>%
    arrange(desc(abs_r)) %>%
    slice_head(n = top_edges_for_ring) %>%
    transmute(from = protein, to = met_display, abs_r = abs_r, sign = sign)
  
  ring_fp <- file.path(figdir, sprintf("block5_pm_ring_top%d.png", top_edges_for_ring))
  make_ring_plot_inner_outer(
    df_edges = pm_ring,
    out_fp = ring_fp,
    title_text = paste0("Protein–metabolite network overview (top ", top_edges_for_ring, " edges by |r|)"),
    cex_prot = 0.70,
    cex_met  = 0.55,
    wrap_met = 18
  )
}


# ---------------------------
# 7) Save cutoff used (reproducibility)
# ---------------------------
cut_fp <- file.path(figdir, "block5_pm_plot_cutoff.txt")
readr::write_lines(
  c(
    paste0("PM edge table: ", basename(pm_fp)),
    paste0("Cut strategy: ", cut_strategy),
    paste0("Quantile (if used): ", cut_quantile),
    paste0("Hard cutoff (if used): ", cut_hard),
    paste0("Chosen cutoff abs_r: ", cut_absr),
    paste0("Two-layer Plot 1: top_edges_for_twolayer_main = ", top_edges_for_twolayer_main),
    paste0("Two-layer Plot 2: top_edges_for_twolayer_div  = ", top_edges_for_twolayer_div),
    paste0("Two-layer metabolite wrap_width_bipartite = ", wrap_width_bipartite),
    "Ring plot: proteins labeled on INNER track; metabolites labeled on OUTER track",
    paste0("Ring tracks: outer height=", ring_track_h_outer, ", inner height=", ring_track_h_inner,
           ", track.margin=", paste(ring_track_margin, collapse = ", "),
           ", outer_label_offset=", ring_outer_label_offset,
           ", inner_label_offset=", ring_inner_label_offset,
           ", png_wh=", ring_png_wh)
  ),
  cut_fp
)
message("Wrote: ", cut_fp)

message("Done: PM figure generation complete.")
