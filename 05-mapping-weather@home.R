# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Mapping the probability ratio (future/current) for years with
# greater than or equal to 4 exceedance days at each weather@home grid cell.
# -------------------------------------------------------------------------

# Packages -----------------------------------------------------------------
# install.packages(c("ggplot2", "viridis", "maps", "mapdata", "tidyr", "dplyr"))
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(tidyr)
library(dplyr)

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
cell_polygons <- NULL
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

  # Build touching quadrilateral polygons from local grid vectors so the map
  # honours the rotated weather@home grid even after conversion to lon/lat.
  build_cell_polygons <- function(df) {
    centres <- df |>
      select(lon_index, lat_index, lon, lat, probability_ratio_ge4) |>
      distinct()

    key <- paste(centres$lon_index, centres$lat_index, sep = "_")
    key_lookup <- setNames(seq_len(nrow(centres)), key)

    get_xy <- function(i, j) {
      k <- paste(i, j, sep = "_")
      idx <- key_lookup[[k]]
      if (is.null(idx) || is.na(idx)) {
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
        probability_ratio_ge4 = centres$probability_ratio_ge4[r]
      )
    }

    if (part_i == 0L) {
      return(data.frame())
    }

    do.call(rbind, polygon_parts[seq_len(part_i)])
  }

  cell_polygons <- build_cell_polygons(finite_map_data)

  if (nrow(cell_polygons) > 0) {
    plot_mode <- "rotated_polygon"
  }
}

# Save plotting data for reproducibility
write.csv(finite_map_data, output_csv, row.names = FALSE)

# Outline of New Zealand in lon/lat for overlay (fallback mode).
nz_outline <- map_data("nz")

# Plot ---------------------------------------------------------------------
if (plot_mode == "rotated_polygon") {
  p_ge4_ratio <- ggplot(cell_polygons, aes(x = lon, y = lat, fill = probability_ratio_ge4)) +
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
    scale_fill_viridis(
      option = "magma",
      name = "Probability\nratio",
      direction = 1
    ) +
    labs(
      title = "weather@home: Probability ratio for >=4 exceedance days",
      subtitle = "Future (3k warmer) / Current decade, by touching rotated-grid polygons",
      x = "Longitude",
      y = "Latitude"
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
