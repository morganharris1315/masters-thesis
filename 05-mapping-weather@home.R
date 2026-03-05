# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Mapping weather@home probability ratios (future/current) for:
# - >= 4 exceedance days
# - >= 5 exceedance days
# - top 10% RX1day years
# - joint event (top 10% RX1day and >= 4 exceedance days)
# Saves each map separately and also saves a combined 3-panel figure
# (>=4 exceedance days, top 10% RX1day, and joint event).
# -------------------------------------------------------------------------

# Packages -----------------------------------------------------------------
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(dplyr)
library(patchwork)

# Input/output paths --------------------------------------------------------
weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"

input_file <- file.path(weatherathome_dir,"weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv")

# Read outputs from the processing script (04) ----------------------------------------------
if (!file.exists(input_file)) {
  stop("Input file not found. Run 04-processing-weather@home.R first.")
}

grid_results <- read.csv(input_file)

required_cols <- c(
  "global_longitude0",
  "global_latitude0",
  "probability_ratio_ge4_future_over_current",
  "probability_ratio_ge5_future_over_current",
  "probability_ratio_rx1day_top10_future_over_current",
  "probability_ratio_joint_top10_ge4_future_over_current"
)

get_ratio_bin_spec <- function(x, n_bins = 6) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r))) {
    stop("Non-finite values found while building discrete ratio bins.")
  }
  
  if (diff(r) == 0) {
    brks <- c(r[1] - 1e-8, r[2] + 1e-8)
  } else {
    brks <- pretty(r, n = n_bins)
    if (min(brks) > r[1]) brks <- c(r[1], brks)
    if (max(brks) < r[2]) brks <- c(brks, r[2])
    brks <- sort(unique(brks))
    if (length(brks) < 2) brks <- c(r[1], r[2])
  }
  
  labels <- paste0(
    format(head(brks, -1), trim = TRUE, scientific = FALSE),
    " to ",
    format(tail(brks, -1), trim = TRUE, scientific = FALSE)
  )
  
  list(
    breaks = brks,
    labels = labels
  )
}

get_fixed_width_bin_spec <- function(x, bin_width = 0.5, min_value = NULL, max_value = NULL) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r))) {
    stop("Non-finite values found while building fixed-width ratio bins.")
  }
  
  if (is.null(min_value)) {
    min_break <- floor(r[1] / bin_width) * bin_width
  } else {
    min_break <- min_value
  }
  
  if (is.null(max_value)) {
    max_break <- ceiling(r[2] / bin_width) * bin_width
  } else {
    max_break <- max_value
  }
  
  if (min_break == max_break) {
    max_break <- min_break + bin_width
  }
  
  brks <- seq(min_break, max_break, by = bin_width)
  
  labels <- paste0(
    format(head(brks, -1), trim = TRUE, scientific = FALSE, nsmall = 1),
    " to ",
    format(tail(brks, -1), trim = TRUE, scientific = FALSE, nsmall = 1)
  )
  
  list(
    breaks = brks,
    labels = labels
  )
}

build_discrete_ratio_bins <- function(x, bin_spec = NULL, n_bins = 6) {
  if (is.null(bin_spec)) {
    bin_spec <- get_ratio_bin_spec(x, n_bins = n_bins)
  }
  
  cut(
    x,
    breaks = bin_spec$breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = bin_spec$labels,
    ordered_result = TRUE
  )
}

# Keep the highest ratio class at the top of the legend for quick scanning.
ratio_legend_guide <- guide_legend(reverse = TRUE)

# Build touching polygons from local grid vectors so the map
build_cell_polygons <- function(df, ratio_col, ratio_bin_col) {
  centres <- df |>
    select(lon_index, lat_index, lon, lat, !!sym(ratio_col), !!sym(ratio_bin_col)) |>
    distinct()
  
  names(centres)[5:6] <- c("ratio_value", "ratio_bin")
  
  key <- paste(centres$lon_index, centres$lat_index, sep = "_")
  key_lookup <- setNames(seq_len(nrow(centres)), key)
  
  get_xy <- function(i, j) {
    k <- paste(i, j, sep = "_")
    idx <- unname(key_lookup[k])
    if (length(idx) == 0L || is.na(idx)) {
      return(c(NA_real_, NA_real_))
    }
    c(centres$lon[idx], centres$lat[idx])
  }
  
  polygon_parts <- vector("list", nrow(centres))
  part_i <- 0L
  
  for (r in seq_len(nrow(centres))) {
    i <- centres$lon_index[r]
    j <- centres$lat_index[r]
    c0 <- c(centres$lon[r], centres$lat[r])
    
    c_w <- get_xy(i - 1, j)
    c_e <- get_xy(i + 1, j)
    c_s <- get_xy(i, j - 1)
    c_n <- get_xy(i, j + 1)
    
    v_i <- if (all(is.finite(c_w)) && all(is.finite(c_e))) {
      (c_e - c_w) / 2
    } else if (all(is.finite(c_e))) {
      c_e - c0
    } else if (all(is.finite(c_w))) {
      c0 - c_w
    } else {
      c(NA_real_, NA_real_)
    }
    
    v_j <- if (all(is.finite(c_s)) && all(is.finite(c_n))) {
      (c_n - c_s) / 2
    } else if (all(is.finite(c_n))) {
      c_n - c0
    } else if (all(is.finite(c_s))) {
      c0 - c_s
    } else {
      c(NA_real_, NA_real_)
    }
    
    if (!all(is.finite(v_i)) || !all(is.finite(v_j))) {
      next
    }
    
    corners <- rbind(
      c0 - 0.5 * v_i - 0.5 * v_j,
      c0 + 0.5 * v_i - 0.5 * v_j,
      c0 + 0.5 * v_i + 0.5 * v_j,
      c0 - 0.5 * v_i + 0.5 * v_j,
      c0 - 0.5 * v_i - 0.5 * v_j
    )
    
    part_i <- part_i + 1L
    polygon_parts[[part_i]] <- data.frame(
      cell_id = paste(i, j, sep = "_"),
      vertex_id = seq_len(5),
      lon = corners[, 1],
      lat = corners[, 2],
      ratio_value = centres$ratio_value[r],
      ratio_bin = centres$ratio_bin[r]
    )
  }
  
  if (part_i == 0L) {
    return(data.frame())
  }
  
  do.call(rbind, polygon_parts[seq_len(part_i)])
}

make_ratio_plot <- function(df_finite, cell_polygons, plot_mode, title_text) {
  nz_outline <- map_data("nz")
  ratio_levels <- levels(df_finite$ratio_bin)
  ratio_palette <- setNames(
    viridisLite::viridis(length(ratio_levels), option = "magma", direction = -1),
    ratio_levels
  )
  
  if (plot_mode == "rotated_polygon") {
    ggplot(cell_polygons, aes(x = lon, y = lat, fill = ratio_bin)) +
      geom_polygon(aes(group = cell_id), colour = NA, linewidth = 0) +
      geom_path(
        data = nz_outline,
        aes(x = long, y = lat, group = group),
        inherit.aes = FALSE,
        colour = "white",
        linewidth = 0.45,
        alpha = 0.9
      ) +
      coord_fixed() +
      scale_fill_manual(
        values = ratio_palette,
        limits = ratio_levels,
        name = "Probability\nratio",
        drop = FALSE
      ) +
      guides(fill = ratio_legend_guide) +
      labs(
        title = title_text,
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 11, lineheight = 0.95),
        axis.title = element_text(face = "bold")
      )
  } else {
    ggplot(df_finite, aes(x = lon, y = lat, fill = ratio_bin)) +
      geom_point(shape = 22, size = 2.2, stroke = 0, colour = NA) +
      geom_path(
        data = nz_outline,
        aes(x = long, y = lat, group = group),
        inherit.aes = FALSE,
        colour = "white",
        linewidth = 0.45,
        alpha = 0.9
      ) +
      coord_fixed() +
      scale_fill_manual(
        values = ratio_palette,
        limits = ratio_levels,
        name = "Probability\nratio",
        drop = FALSE
      ) +
      guides(fill = ratio_legend_guide) +
      labs(
        title = title_text,
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 11, lineheight = 0.95),
        axis.title = element_text(face = "bold")
      )
  }
}

build_metric_layers <- function(ratio_col, bin_spec = NULL) {
  base_data <- data.frame(
    lon_index = grid_results$lon_index,
    lat_index = grid_results$lat_index,
    lon = grid_results$global_longitude0,
    lat = grid_results$global_latitude0,
    ratio_value = grid_results[[ratio_col]]
  )
  
  finite_data <- base_data[is.finite(base_data$ratio_value), ]
  if (nrow(finite_data) == 0) {
    stop(sprintf("No finite probability-ratio values available to plot for %s.", ratio_col))
  }
  
  finite_data$ratio_bin <- build_discrete_ratio_bins(finite_data$ratio_value, bin_spec = bin_spec)
  
  plot_mode <- "point"
  cell_polygons <- data.frame()
  
  has_rotated_metadata <- all(c("lon_index", "lat_index", "longitude0", "latitude0") %in% names(grid_results))
  if (has_rotated_metadata) {
    cell_polygons <- build_cell_polygons(finite_data, "ratio_value", "ratio_bin")
    if (nrow(cell_polygons) > 0) {
      plot_mode <- "rotated_polygon"
    }
  }
  
  list(
    base_data = base_data,
    finite_data = finite_data,
    cell_polygons = cell_polygons,
    plot_mode = plot_mode
  )
}

# Build plot layers for all requested metrics ------------------------------
shared_ratio_values <- c(
  grid_results$probability_ratio_ge4_future_over_current,
  grid_results$probability_ratio_rx1day_top10_future_over_current,
  grid_results$probability_ratio_joint_top10_ge4_future_over_current
)
shared_ratio_values <- shared_ratio_values[is.finite(shared_ratio_values)]

if (length(shared_ratio_values) == 0) {
  stop("No finite values found for shared binning across >=4, top 10%, and joint ratios.")
}

shared_ratio_bin_spec <- get_fixed_width_bin_spec(
  shared_ratio_values,
  bin_width = 0.5,
  min_value = 0,
  max_value = 6.5
)

layers_ge4 <- build_metric_layers(
  "probability_ratio_ge4_future_over_current",
  bin_spec = shared_ratio_bin_spec
)
layers_ge5 <- build_metric_layers("probability_ratio_ge5_future_over_current")
layers_top10 <- build_metric_layers(
  "probability_ratio_rx1day_top10_future_over_current",
  bin_spec = shared_ratio_bin_spec
)
layers_joint <- build_metric_layers(
  "probability_ratio_joint_top10_ge4_future_over_current",
  bin_spec = shared_ratio_bin_spec
)

# Save plotting data for reproducibility
output_map_data <- grid_results |>
  select(
    lon_index,
    lat_index,
    longitude0,
    latitude0,
    global_longitude0,
    global_latitude0,
    probability_ratio_ge4_future_over_current,
    probability_ratio_ge5_future_over_current,
    probability_ratio_rx1day_top10_future_over_current,
    probability_ratio_joint_top10_ge4_future_over_current
  )
# Creating maps -------------------------------------------
p_ge4_ratio <- make_ratio_plot(
  df_finite = layers_ge4$finite_data,
  cell_polygons = layers_ge4$cell_polygons,
  plot_mode = layers_ge4$plot_mode,
  title_text = "≥ 4 Exceedance Days")

p_ge5_ratio <- make_ratio_plot(
  df_finite = layers_ge5$finite_data,
  cell_polygons = layers_ge5$cell_polygons,
  plot_mode = layers_ge5$plot_mode,
  title_text = "≥ 5 Exceedance Days")

p_top10_ratio <- make_ratio_plot(
  df_finite = layers_top10$finite_data,
  cell_polygons = layers_top10$cell_polygons,
  plot_mode = layers_top10$plot_mode,
  title_text = "Top 10% RX1day")

p_joint_ratio <- make_ratio_plot(
  df_finite = layers_joint$finite_data,
  cell_polygons = layers_joint$cell_polygons,
  plot_mode = layers_joint$plot_mode,
  title_text = "Top 10% RX1day and ≥ 4 Exceedance Days")

p_ge4_ratio
p_ge5_ratio
p_top10_ratio
p_joint_ratio

# Combined 3-panel figure ≥4, top 10%, and joint event --------
p_combined <- (p_ge4_ratio + p_top10_ratio + p_joint_ratio) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "right")

p_combined

# Saving outputs ----------------------------------------------------------
output_csv <- file.path(weatherathome_dir, "weather@home_probability_ratio_map_data.csv")
output_png_ge4 <- file.path(weatherathome_dir, "weather@home_probability_ratio_ge4_map.png")
output_png_ge5 <- file.path(weatherathome_dir, "weather@home_probability_ratio_ge5_map.png")
output_png_top10 <- file.path(weatherathome_dir, "weather@home_probability_ratio_rx1day_top10_map.png")
output_png_joint <- file.path(weatherathome_dir, "weather@home_probability_ratio_joint_top10_ge4_map.png")
output_png_combined <- file.path(weatherathome_dir, "weather@home_probability_ratio_ge4_top10_joint_combined_map.png")

write.csv(output_map_data, output_csv, row.names = FALSE)

ggsave(filename = output_png_ge4, plot = p_ge4_ratio, width = 8, height = 7, dpi = 300)
ggsave(filename = output_png_ge5, plot = p_ge5_ratio, width = 8, height = 7, dpi = 300)
ggsave(filename = output_png_top10, plot = p_top10_ratio, width = 8, height = 7, dpi = 300)
ggsave(filename = output_png_joint, plot = p_joint_ratio, width = 8, height = 7, dpi = 300)
ggsave(filename = output_png_combined, plot = p_combined, width = 18, height = 6.5, dpi = 300)

# Diagnostics ---------------------------------------------------------
cat("Finite cell count (>=4):", nrow(layers_ge4$finite_data), "\n")
cat("Finite cell count (>=5):", nrow(layers_ge5$finite_data), "\n")
cat("Finite cell count (top 10%):", nrow(layers_top10$finite_data), "\n")
cat("Finite cell count (joint top10 + >=4):", nrow(layers_joint$finite_data), "\n")

cat("Plot mode (>=4):", layers_ge4$plot_mode, "\n")
cat("Plot mode (>=5):", layers_ge5$plot_mode, "\n")
cat("Plot mode (top 10%):", layers_top10$plot_mode, "\n")
cat("Plot mode (joint top10 + >=4):", layers_joint$plot_mode, "\n")

cat("Saved map to:", output_png_ge4, "\n")
cat("Saved map to:", output_png_ge5, "\n")
cat("Saved map to:", output_png_top10, "\n")
cat("Saved map to:", output_png_joint, "\n")
cat("Saved combined map to:", output_png_combined, "\n")
