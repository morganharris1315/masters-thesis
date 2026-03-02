# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Mapping the probability ratio (future/current) for years with
# greater than or equal to 4 exceedance days at each weather@home grid cell.
# -------------------------------------------------------------------------

# Packages -----------------------------------------------------------------
# install.packages(c("ggplot2", "viridis", "maps", "mapdata", "tidyr"))
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(tidyr)

# Input/output paths --------------------------------------------------------
weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"

input_file <- file.path(
  weatherathome_dir,
  "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv"
)

output_png <- file.path(
  weatherathome_dir,
  "weather@home_probability_ratio_ge4_map.png"
)

output_csv <- file.path(
  weatherathome_dir,
  "weather@home_probability_ratio_ge4_map_data.csv"
)

# Read outputs from Script 04 ----------------------------------------------
if (!file.exists(input_file)) {
  stop("Input file not found. Run 04-processing-weather@home.R first.")
}

grid_results <- read.csv(input_file)

required_cols <- c(
  "global_longitude0",
  "global_latitude0",
  "probability_ratio_ge4_future_over_current"
)

missing_cols <- setdiff(required_cols, names(grid_results))
if (length(missing_cols) > 0) {
  stop(
    sprintf(
      "Input file is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Keep finite values for colour scaling; keep Inf separately for diagnostics.
map_data <- data.frame(
  lon = grid_results$global_longitude0,
  lat = grid_results$global_latitude0,
  probability_ratio_ge4 = grid_results$probability_ratio_ge4_future_over_current
)

map_data$ratio_class <- ifelse(
  is.infinite(map_data$probability_ratio_ge4),
  "Inf (current=0, future>0)",
  "Finite"
)

finite_map_data <- map_data[is.finite(map_data$probability_ratio_ge4), ]

if (nrow(finite_map_data) == 0) {
  stop("No finite probability-ratio values available to plot.")
}

# For aligned model-grid plotting we need rotated-grid coordinates and matrix
# indices exported by 04-processing-weather@home.R.
plot_mode <- "point"
if (all(c("lon_index", "lat_index", "longitude0", "latitude0") %in% names(grid_results))) {
  finite_map_data <- grid_results[is.finite(grid_results$probability_ratio_ge4_future_over_current), c(
    "lon_index", "lat_index", "longitude0", "latitude0",
    "global_longitude0", "global_latitude0",
    "probability_ratio_ge4_future_over_current"
  )]

  finite_map_data <- finite_map_data |>
    rename(
      lon = global_longitude0,
      lat = global_latitude0,
      probability_ratio_ge4 = probability_ratio_ge4_future_over_current
    )

  # Reconstruct ordered matrix for the georeferenced raster and keep one value
  # per cell even if duplicate lon/lat values exist.
  ratio_wide <- finite_map_data |>
    select(lon_index, lat_index, probability_ratio_ge4) |>
    distinct() |>
    pivot_wider(
      names_from = lat_index,
      values_from = probability_ratio_ge4
    ) |>
    arrange(lon_index)

  lon_lookup <- finite_map_data |>
    select(lon_index, longitude0) |>
    distinct() |>
    arrange(lon_index)

  lat_lookup <- finite_map_data |>
    select(lat_index, latitude0) |>
    distinct() |>
    arrange(lat_index)

  ratio_matrix <- as.matrix(ratio_wide[, -1, drop = FALSE])
  rownames(ratio_matrix) <- lon_lookup$longitude0
  colnames(ratio_matrix) <- lat_lookup$latitude0

  # geom_raster requires equal spacing, so use geom_tile with explicit width /
  # height from rotated-grid spacing.
  tile_width <- median(diff(sort(unique(finite_map_data$longitude0))), na.rm = TRUE)
  tile_height <- median(diff(sort(unique(finite_map_data$latitude0))), na.rm = TRUE)

  # Keep rotated-grid coordinates for plotting so model cells remain aligned
  # and touching in the native weather@home grid.
  rotated_grid_tiles <- finite_map_data |>
    transmute(
      x = lon,
      y = lat,
      probability_ratio_ge4
    )

  # Build a land-mask outline directly on the rotated grid by testing whether
  # each cell centre falls on NZ land in geographic lon/lat space. Drawing the
  # contour in rotated x/y keeps the coastline aligned with the rotated tiles.
  rotated_grid_tiles$nz_land <- ifelse(
    is.na(map.where("nz", rotated_grid_tiles$lon, rotated_grid_tiles$lat)),
    0,
    1
  )

  plot_mode <- "rotated_tile"
}

# Save plotting data for reproducibility
write.csv(finite_map_data, output_csv, row.names = FALSE)

# Outline of New Zealand in lon/lat for overlay (fallback mode).
nz_outline <- map_data("nz")

# Plot ---------------------------------------------------------------------
if (plot_mode == "rotated_tile") {
  p_ge4_ratio <- ggplot(rotated_grid_tiles, aes(x = x, y = y, fill = probability_ratio_ge4)) +
    geom_tile(width = tile_width, height = tile_height) +
    geom_contour(
      aes(z = nz_land),
      breaks = 0.5,
      colour = "white",
      linewidth = 0.45,
      alpha = 0.9
    ) +
    coord_fixed() +
    scale_fill_viridis(
      option = "magma",
      name = "Probability\nratio",
      direction = 1
    ) +
    labs(
      title = "weather@home: Probability ratio for >=4 exceedance days",
      subtitle = "Future (3k warmer) / Current decade, rotated grid with NZ land-mask outline",
      x = "Rotated longitude0",
      y = "Rotated latitude0"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
} else {
  # Fallback for older CSVs that do not contain rotated-grid metadata.
  p_ge4_ratio <- ggplot(finite_map_data, aes(x = lon, y = lat, color = probability_ratio_ge4)) +
    geom_point(shape = 15, size = 2.2, stroke = 0) +
    geom_path(
      data = nz_outline,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      colour = "white",
      linewidth = 0.45,
      alpha = 0.9
    ) +
    coord_fixed() +
    scale_color_viridis(
      option = "magma",
      name = "Probability\nratio",
      direction = 1
    ) +
    labs(
      title = "weather@home: Probability ratio for >=4 exceedance days",
      subtitle = "Future (3k warmer) / Current decade, by model grid cell",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

if (interactive()) {
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  print(p_ge4_ratio)
}

ggsave(
  filename = output_png,
  plot = p_ge4_ratio,
  width = 8,
  height = 7,
  dpi = 300
)

# Quick diagnostics ---------------------------------------------------------
cat("Finite cell count:", nrow(finite_map_data), "\n")
cat("Infinite cell count:", sum(is.infinite(map_data$probability_ratio_ge4)), "\n")
cat("Unique longitudes:", length(unique(map_data$lon)), "\n")
cat("Unique latitudes:", length(unique(map_data$lat)), "\n")
cat("Plot mode:", plot_mode, "\n")
cat("Saved map to:", output_png, "\n")
