# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Build weather@home probability-ratio maps using NCIP5 grid + LSE land mask.
# Outputs:
# - main combined map (>=4, top 10% Rx1day, and joint)
# - >=5 exceedance map (for supplementary material)
# - masked CSV + masked NetCDF on NCIP5 grid
# -------------------------------------------------------------------------

library(ggplot2)
library(viridis)
library(maps)
library(dplyr)
library(patchwork)
library(grid)
library(colorspace)
library(RNetCDF)

# Colouring  --------------------------------------------------------------
#blue_palette <- choose_palette()
blues <- blue_palette(14)

# Input/output paths -------------------------------------------------------
weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"
nc_template_file <- file.path(weatherathome_dir, "current_decade", "item5216_daily_mean_a000_2006-07_2007-06-NZtrim-mm.nc")
lse_mask_file <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/Land-Sea Mask for Weather@home Data.xlsx"
lse_mask_sheet <- "capro"
lse_mask_transpose <- TRUE

# Reading processed data ---------------------------------------------------
grid_results <- read.csv(
  file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv")
)

ratio_columns <- c(
  "probability_ratio_ge4_future_over_current",
  "probability_ratio_ge5_future_over_current",
  "probability_ratio_rx1day_top10_future_over_current",
  "probability_ratio_joint_top10_ge4_future_over_current"
)

load_nc_grid_template <- function(nc_file) {
  nc <- open.nc(nc_file)
  on.exit(close.nc(nc), add = TRUE)

  longitude0 <- var.get.nc(nc, "longitude0")
  latitude0 <- var.get.nc(nc, "latitude0")
  global_longitude0 <- var.get.nc(nc, "global_longitude0")
  global_latitude0 <- var.get.nc(nc, "global_latitude0")

  grid_template <- expand.grid(
    lon_index = seq_along(longitude0),
    lat_index = seq_along(latitude0)
  )

  data.frame(
    lon_index = grid_template$lon_index,
    lat_index = grid_template$lat_index,
    longitude0 = longitude0[grid_template$lon_index],
    latitude0 = latitude0[grid_template$lat_index],
    global_longitude0 = global_longitude0[cbind(grid_template$lon_index, grid_template$lat_index)],
    global_latitude0 = global_latitude0[cbind(grid_template$lon_index, grid_template$lat_index)]
  )
}

load_nc_cell_polygon_template <- function(nc_file, grid_template) {
  nc <- open.nc(nc_file)
  on.exit(close.nc(nc), add = TRUE)

  nlon <- length(unique(grid_template$lon_index))
  nlat <- length(unique(grid_template$lat_index))

  read_var_if_exists <- function(var_name) {
    tryCatch(var.get.nc(nc, var_name), error = function(e) NULL)
  }

  build_from_1d_bounds <- function(lon_bounds, lat_bounds) {
    if (length(dim(lon_bounds)) != 2 || length(dim(lat_bounds)) != 2) {
      return(NULL)
    }

    if (nrow(lon_bounds) != nlon || nrow(lat_bounds) != nlat || ncol(lon_bounds) < 2 || ncol(lat_bounds) < 2) {
      return(NULL)
    }

    polygon_parts <- vector("list", nlon * nlat)
    part_i <- 0L

    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        lon_lo <- lon_bounds[i, 1]
        lon_hi <- lon_bounds[i, 2]
        lat_lo <- lat_bounds[j, 1]
        lat_hi <- lat_bounds[j, 2]

        if (!all(is.finite(c(lon_lo, lon_hi, lat_lo, lat_hi)))) next

        corners <- rbind(
          c(lon_lo, lat_lo),
          c(lon_hi, lat_lo),
          c(lon_hi, lat_hi),
          c(lon_lo, lat_hi),
          c(lon_lo, lat_lo)
        )

        part_i <- part_i + 1L
        polygon_parts[[part_i]] <- data.frame(
          lon_index = i,
          lat_index = j,
          cell_id = paste(i, j, sep = "_"),
          vertex_id = seq_len(5),
          lon = corners[, 1],
          lat = corners[, 2]
        )
      }
    }

    if (part_i == 0L) return(NULL)
    do.call(rbind, polygon_parts[seq_len(part_i)])
  }

  build_from_2d_bounds <- function(lon_bounds, lat_bounds) {
    d_lon <- dim(lon_bounds)
    d_lat <- dim(lat_bounds)

    if (length(d_lon) != 3 || length(d_lat) != 3 || d_lon[3] < 3 || d_lat[3] < 3) {
      return(NULL)
    }

    if (!all(d_lon[1:2] == c(nlon, nlat)) || !all(d_lat[1:2] == c(nlon, nlat))) {
      return(NULL)
    }

    n_vertices <- min(d_lon[3], d_lat[3])
    polygon_parts <- vector("list", nlon * nlat)
    part_i <- 0L

    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        lon_seq <- as.numeric(lon_bounds[i, j, seq_len(n_vertices)])
        lat_seq <- as.numeric(lat_bounds[i, j, seq_len(n_vertices)])

        if (!all(is.finite(lon_seq)) || !all(is.finite(lat_seq))) next

        if (lon_seq[1] != lon_seq[length(lon_seq)] || lat_seq[1] != lat_seq[length(lat_seq)]) {
          lon_seq <- c(lon_seq, lon_seq[1])
          lat_seq <- c(lat_seq, lat_seq[1])
        }

        part_i <- part_i + 1L
        polygon_parts[[part_i]] <- data.frame(
          lon_index = i,
          lat_index = j,
          cell_id = paste(i, j, sep = "_"),
          vertex_id = seq_along(lon_seq),
          lon = lon_seq,
          lat = lat_seq
        )
      }
    }

    if (part_i == 0L) return(NULL)
    do.call(rbind, polygon_parts[seq_len(part_i)])
  }

  bounds_candidates <- list(
    c("global_longitude0_bounds", "global_latitude0_bounds"),
    c("global_longitude0_bnds", "global_latitude0_bnds"),
    c("longitude0_bounds", "latitude0_bounds"),
    c("longitude0_bnds", "latitude0_bnds")
  )

  for (pair in bounds_candidates) {
    lon_bounds <- read_var_if_exists(pair[1])
    lat_bounds <- read_var_if_exists(pair[2])

    if (is.null(lon_bounds) || is.null(lat_bounds)) next

    polygons <- build_from_2d_bounds(lon_bounds, lat_bounds)
    if (is.null(polygons)) {
      polygons <- build_from_1d_bounds(lon_bounds, lat_bounds)
    }

    if (!is.null(polygons)) {
      message(sprintf("Using NC cell bounds from variables '%s' and '%s'.", pair[1], pair[2]))
      return(polygons)
    }
  }

  message("No NC cell-bound variables found with supported shapes; using neighbour-based polygon reconstruction.")
  NULL
}

read_lse_mask <- function(mask_file, nlon, nlat) {
  ext <- tolower(tools::file_ext(mask_file))

  if (ext == "csv") {
    raw <- read.csv(mask_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("readxl package is required to read Excel LSE mask files.")
    }

    available_sheets <- readxl::excel_sheets(mask_file)
    sheet_match <- available_sheets[tolower(available_sheets) == tolower(lse_mask_sheet)]
    sheet_to_read <- if (length(sheet_match) > 0) sheet_match[1] else available_sheets[1]

    raw <- readxl::read_excel(mask_file, sheet = sheet_to_read, col_names = TRUE)
    raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported LSE mask format. Use .csv, .xlsx or .xls")
  }

  numeric_raw <- suppressWarnings(as.data.frame(lapply(raw, as.numeric)))

  to_matrix <- function(df) as.matrix(df)
  candidates <- list(
    to_matrix(numeric_raw),
    to_matrix(numeric_raw[, -1, drop = FALSE]),
    to_matrix(numeric_raw[-1, , drop = FALSE]),
    to_matrix(numeric_raw[-1, -1, drop = FALSE])
  )

  dims_ok <- vapply(candidates, function(m) {
    is.matrix(m) && nrow(m) == nlon && ncol(m) == nlat
  }, logical(1))

  if (!any(dims_ok)) {
    stop(sprintf("Could not coerce LSE mask to %d x %d matrix. Check the file layout.", nlon, nlat))
  }

  mask_matrix <- candidates[[which(dims_ok)[1]]]

  if (isTRUE(lse_mask_transpose)) {
    mask_matrix <- t(mask_matrix)
  }

  data.frame(
    lon_index = rep(seq_len(nlon), times = nlat),
    lat_index = rep(seq_len(nlat), each = nlon),
    lse_mask_value = as.vector(mask_matrix)
  )
}

apply_lse_mask_to_metrics <- function(grid_df, ratio_cols) {
  masked <- grid_df
  mask_is_land <- !is.na(masked$lse_mask_value)

  for (col_name in ratio_cols) {
    masked[[col_name]][!mask_is_land] <- NA_real_
  }

  masked$lse_is_land <- mask_is_land
  masked
}

# Build NC-anchored table --------------------------------------------------
nc_grid_template <- load_nc_grid_template(nc_template_file)
nc_cell_polygon_template <- load_nc_cell_polygon_template(nc_template_file, nc_grid_template)

nc_grid_metrics <- nc_grid_template |>
  left_join(grid_results, by = c("lon_index", "lat_index"))

missing_metrics <- rowSums(is.na(nc_grid_metrics[, ratio_columns]))
if (any(missing_metrics == length(ratio_columns))) {
  stop("Some NC grid cells did not match CSV probability-ratio rows. Check lon_index/lat_index alignment.")
}

lse_mask_df <- read_lse_mask(
  mask_file = lse_mask_file,
  nlon = length(unique(nc_grid_template$lon_index)),
  nlat = length(unique(nc_grid_template$lat_index))
)

nc_grid_metrics <- nc_grid_metrics |>
  left_join(lse_mask_df, by = c("lon_index", "lat_index"))

if (all(is.na(nc_grid_metrics$lse_mask_value))) {
  stop("LSE mask contains only NA values after alignment.")
}

nc_grid_metrics_masked <- apply_lse_mask_to_metrics(nc_grid_metrics, ratio_columns)

extract_masked_ratio_grid <- function(grid_df, ratio_col) {
  grid_df |>
    select(
      lon_index, lat_index,
      longitude0, latitude0,
      global_longitude0, global_latitude0,
      lse_mask_value, lse_is_land,
      all_of(ratio_col)
    )
}

# Shared helpers -----------------------------------------------------------
get_fixed_width_bin_spec <- function(x, bin_width = 0.5, min_value = 0, max_value = NULL) {
  r <- range(x, na.rm = TRUE)
  min_break <- min_value
  max_break <- if (is.null(max_value)) ceiling(r[2] / bin_width) * bin_width else max_value

  if (min_break == max_break) {
    max_break <- min_break + bin_width
  }

  brks <- seq(min_break, max_break, by = bin_width)

  list(
    breaks = brks,
    labels = format(head(brks, -1), trim = TRUE, scientific = FALSE, nsmall = 1)
  )
}

build_discrete_ratio_bins <- function(x, bin_spec) {
  cut(
    x,
    breaks = bin_spec$breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = bin_spec$labels,
    ordered_result = TRUE
  )
}

build_cell_polygons <- function(df) {
  has_ratio_bin <- "ratio_bin" %in% names(df)
  centres <- df |>
    select(lon_index, lat_index, lon, lat, ratio_value) |>
    distinct()

  if (has_ratio_bin) {
    ratio_bins <- df |>
      select(lon_index, lat_index, ratio_bin) |>
      distinct()

    centres <- centres |>
      left_join(ratio_bins, by = c("lon_index", "lat_index"))
  } else {
    centres$ratio_bin <- NA_character_
  }

  key <- paste(centres$lon_index, centres$lat_index, sep = "_")
  key_lookup <- setNames(seq_len(nrow(centres)), key)

  get_xy <- function(i, j) {
    k <- paste(i, j, sep = "_")
    idx <- unname(key_lookup[k])
    if (length(idx) == 0L || is.na(idx)) return(c(NA_real_, NA_real_))
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

    if (!all(is.finite(v_i)) || !all(is.finite(v_j))) next

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

  if (part_i == 0L) return(data.frame())
  do.call(rbind, polygon_parts[seq_len(part_i)])
}

build_metric_layer <- function(ratio_col) {
  base_data <- data.frame(
    lon_index = nc_grid_metrics_masked$lon_index,
    lat_index = nc_grid_metrics_masked$lat_index,
    lon = nc_grid_metrics_masked$global_longitude0,
    lat = nc_grid_metrics_masked$global_latitude0,
    ratio_value = nc_grid_metrics_masked[[ratio_col]]
  )

  finite_data <- base_data[is.finite(base_data$ratio_value), ]
  if (nrow(finite_data) == 0) {
    stop(sprintf("No finite values to plot for %s after LSE masking.", ratio_col))
  }

  cell_polygons <- if (is.null(nc_cell_polygon_template)) {
    build_cell_polygons(finite_data)
  } else {
    nc_cell_polygon_template |>
      inner_join(
        finite_data |>
          select(lon_index, lat_index, ratio_value),
        by = c("lon_index", "lat_index")
      )
  }

  list(
    finite_data = finite_data,
    cell_polygons = cell_polygons,
    keep_cell_ids = unique(cell_polygons$cell_id)
  )
}

make_nz_ratio_plot <- function(
    layer_obj,
    title_text,
    ratio_breaks,
    blues,
    show_legend = TRUE,
    legend_height_cm = 18,
    legend_width_cm = 0.65) {
  cell_polygons_nz <- layer_obj$cell_polygons[layer_obj$cell_polygons$cell_id %in% layer_obj$keep_cell_ids, ]

  if (nrow(cell_polygons_nz) == 0) {
    stop(sprintf("No cells found for '%s'.", title_text))
  }

  ratio_limits <- range(ratio_breaks)
  ratio_labels <- format(ratio_breaks, trim = TRUE, scientific = FALSE, nsmall = 1)
  nz_outline <- map_data("nz")

  legend_position <- if (isTRUE(show_legend)) "right" else "none"
  legend_height <- unit(legend_height_cm, "cm")
  legend_width <- unit(legend_width_cm, "cm")

  ggplot(cell_polygons_nz, aes(x = lon, y = lat, fill = ratio_value)) +
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
    scale_fill_stepsn(
      colours = blues,
      breaks = ratio_breaks,
      labels = ratio_labels,
      limits = ratio_limits,
      oob = scales::squish,
      show.limits = TRUE,
      name = "Probability Ratio",
      guide = guide_coloursteps(
        reverse = FALSE,
        even.steps = TRUE,
        show.limits = TRUE,
        barheight = legend_height,
        barwidth = legend_width,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 9),
      legend.key.height = unit(1.05, "cm"),
      legend.key.width = legend_width,
      legend.margin = margin(4, 8, 4, 2),
      legend.position = legend_position
    )
}

# Build layers --------------------------------------------------------------
layers_ge4 <- build_metric_layer("probability_ratio_ge4_future_over_current")
layers_ge5 <- build_metric_layer("probability_ratio_ge5_future_over_current")
layers_top10 <- build_metric_layer("probability_ratio_rx1day_top10_future_over_current")
layers_joint <- build_metric_layer("probability_ratio_joint_top10_ge4_future_over_current")

# Shared bins using LSE-masked cells only ----------------------------------
intersection_values <- c(
  layers_ge4$finite_data$ratio_value,
  layers_ge5$finite_data$ratio_value,
  layers_top10$finite_data$ratio_value,
  layers_joint$finite_data$ratio_value
)

intersection_values <- intersection_values[is.finite(intersection_values)]
if (length(intersection_values) == 0) {
  stop("No finite probability-ratio values found after LSE masking.")
}

ratio_bin_spec <- get_fixed_width_bin_spec(intersection_values, bin_width = 0.5, min_value = 0)
ratio_breaks <- ratio_bin_spec$breaks

add_ratio_bins <- function(layer_obj) {
  layer_obj$finite_data$ratio_bin <- build_discrete_ratio_bins(layer_obj$finite_data$ratio_value, ratio_bin_spec)
  layer_obj$cell_polygons$ratio_bin <- build_discrete_ratio_bins(layer_obj$cell_polygons$ratio_value, ratio_bin_spec)
  layer_obj
}

layers_ge4 <- add_ratio_bins(layers_ge4)
layers_ge5 <- add_ratio_bins(layers_ge5)
layers_top10 <- add_ratio_bins(layers_top10)
layers_joint <- add_ratio_bins(layers_joint)

# Build requested plots -----------------------------------------------------
p_ge4 <- make_nz_ratio_plot(
  layers_ge4, "(b) Years with ≥4 exceedances",
  ratio_breaks, blues
)

p_top10 <- make_nz_ratio_plot(
  layers_top10, "(a) Years with extreme Rx1day",
  ratio_breaks, blues, show_legend = FALSE
)

p_joint <- make_nz_ratio_plot(
  layers_joint, "(c) Years with extreme Rx1day AND ≥4 exceedances",
  ratio_breaks, blues, show_legend = FALSE
) +
  theme(plot.title = element_text(hjust = 0.5))

p_ge5 <- make_nz_ratio_plot(
  layers_ge5, "Years with ≥5 exceedances",
  ratio_breaks, blues, legend_height_cm = 13
)

combined_design <- c(
  area(t = 1, l = 1, b = 1, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 2, l = 1, b = 2, r = 2),
  area(t = 1, l = 3, b = 2, r = 3)
)

p_combined <- (p_top10 + p_ge4 + p_joint + guide_area()) +
  plot_layout(
    design = combined_design,
    guides = "collect",
    widths = c(1, 1, 0.18),
    heights = c(1, 1)
  ) +
  plot_annotation(
    theme = theme(
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(5, 5, 5, 5)
    )
  )

p_combined
p_ge5

# Save outputs --------------------------------------------------------------
output_map_data <- nc_grid_metrics_masked |>
  select(
    lon_index, lat_index,
    longitude0, latitude0,
    global_longitude0, global_latitude0,
    lse_mask_value, lse_is_land,
    all_of(ratio_columns)
  )

write.csv(
  output_map_data,
  file.path(weatherathome_dir, "weather@home_probability_ratio_map_data_lse_masked.csv"),
  row.names = FALSE
)

write.csv(
  extract_masked_ratio_grid(output_map_data, "probability_ratio_ge4_future_over_current"),
  file.path(weatherathome_dir, "weather@home_probability_ratio_ge4_lse_masked.csv"),
  row.names = FALSE
)

write.csv(
  extract_masked_ratio_grid(output_map_data, "probability_ratio_ge5_future_over_current"),
  file.path(weatherathome_dir, "weather@home_probability_ratio_ge5_lse_masked.csv"),
  row.names = FALSE
)

write.csv(
  extract_masked_ratio_grid(output_map_data, "probability_ratio_rx1day_top10_future_over_current"),
  file.path(weatherathome_dir, "weather@home_probability_ratio_rx1day_top10_lse_masked.csv"),
  row.names = FALSE
)

write.csv(
  extract_masked_ratio_grid(output_map_data, "probability_ratio_joint_top10_ge4_future_over_current"),
  file.path(weatherathome_dir, "weather@home_probability_ratio_joint_top10_ge4_lse_masked.csv"),
  row.names = FALSE
)

write_probability_ratio_netcdf <- function(template_nc_file, grid_df, output_nc_file, ratio_cols) {
  template_nc <- open.nc(template_nc_file)
  on.exit(close.nc(template_nc), add = TRUE)

  longitude0 <- var.get.nc(template_nc, "longitude0")
  latitude0 <- var.get.nc(template_nc, "latitude0")
  global_longitude0 <- var.get.nc(template_nc, "global_longitude0")
  global_latitude0 <- var.get.nc(template_nc, "global_latitude0")

  nlon <- length(longitude0)
  nlat <- length(latitude0)

  build_grid_matrix <- function(col_name) {
    mat <- matrix(NA_real_, nrow = nlon, ncol = nlat)
    mat[cbind(grid_df$lon_index, grid_df$lat_index)] <- grid_df[[col_name]]
    mat
  }

  nc_out <- create.nc(output_nc_file)
  on.exit(close.nc(nc_out), add = TRUE)

  dim.def.nc(nc_out, "longitude0", nlon)
  dim.def.nc(nc_out, "latitude0", nlat)

  var.def.nc(nc_out, "longitude0", "NC_DOUBLE", "longitude0")
  var.def.nc(nc_out, "latitude0", "NC_DOUBLE", "latitude0")
  var.def.nc(nc_out, "global_longitude0", "NC_DOUBLE", c("longitude0", "latitude0"))
  var.def.nc(nc_out, "global_latitude0", "NC_DOUBLE", c("longitude0", "latitude0"))
  var.def.nc(nc_out, "lse_mask_value", "NC_DOUBLE", c("longitude0", "latitude0"))

  for (col_name in ratio_cols) {
    var.def.nc(nc_out, col_name, "NC_DOUBLE", c("longitude0", "latitude0"))
  }

  att.put.nc(nc_out, "NC_GLOBAL", "title", "NC_CHAR", "Weather@home probability-ratio grids on NCIP5 template grid with LSE masking")
  att.put.nc(nc_out, "NC_GLOBAL", "source_template", "NC_CHAR", basename(template_nc_file))
  att.put.nc(nc_out, "NC_GLOBAL", "mask_file", "NC_CHAR", basename(lse_mask_file))
  att.put.nc(nc_out, "NC_GLOBAL", "mask_sheet", "NC_CHAR", lse_mask_sheet)
  att.put.nc(nc_out, "NC_GLOBAL", "mask_transposed", "NC_CHAR", if (isTRUE(lse_mask_transpose)) "TRUE" else "FALSE")

  var.put.nc(nc_out, "longitude0", longitude0)
  var.put.nc(nc_out, "latitude0", latitude0)
  var.put.nc(nc_out, "global_longitude0", global_longitude0)
  var.put.nc(nc_out, "global_latitude0", global_latitude0)
  var.put.nc(nc_out, "lse_mask_value", build_grid_matrix("lse_mask_value"))

  for (col_name in ratio_cols) {
    var.put.nc(nc_out, col_name, build_grid_matrix(col_name))
  }

  sync.nc(nc_out)
}

write_probability_ratio_netcdf(
  template_nc_file = nc_template_file,
  grid_df = output_map_data,
  output_nc_file = file.path(weatherathome_dir, "weather@home_probability_ratio_grid_from_NCIP5_lse_masked.nc"),
  ratio_cols = ratio_columns
)

ggsave(
  filename = file.path(weatherathome_dir, "weather@home_probability_ratio_ge4_top10_joint_lse_masked_combined_map.png"),
  plot = p_combined, width = 12, height = 10, dpi = 300
)

ggsave(
  filename = file.path(weatherathome_dir, "weather@home_probability_ratio_ge5_lse_masked_map.png"),
  plot = p_ge5, width = 8, height = 7, dpi = 300
)

# Quick checks --------------------------------------------------------------
cat("LSE land cell count:", sum(output_map_data$lse_is_land, na.rm = TRUE), "\n")
cat("Finite LSE-masked cells (>=4):", sum(is.finite(output_map_data$probability_ratio_ge4_future_over_current)), "\n")
cat("Finite LSE-masked cells (>=5):", sum(is.finite(output_map_data$probability_ratio_ge5_future_over_current)), "\n")
cat("Finite LSE-masked cells (top 10%):", sum(is.finite(output_map_data$probability_ratio_rx1day_top10_future_over_current)), "\n")
cat("Finite LSE-masked cells (joint):", sum(is.finite(output_map_data$probability_ratio_joint_top10_ge4_future_over_current)), "\n")
