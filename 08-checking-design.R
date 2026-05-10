# -------------------------------------------------------------------------
# 08-checking-design.R
# -------------------------------------------------------------------------
# April 2026
# Aligning Extreme threshold with % of at least 4 excedance day in current day. 
# -------------------------------------------------------------------------

# Matching Rx1day Extreme threshold with 1.1% -----------------------------

rx1day_threshold_98.9_current <- apply(
  current_rx_array,
  c(1, 2),
  quantile,
  probs = 0.989,
  na.rm = TRUE,
  type = 7
)

rx1day_threshold_98.9_future <- apply(
  future_rx_array,
  c(1, 2),
  quantile,
  probs = 0.989,
  na.rm = TRUE,
  type = 7
)

current_prop_rx1day_98.9 <- calc_rx1day_top10_proportion(
  current_rx_array,
  rx1day_threshold_98.9_current
)
future_prop_rx1day_98.9 <- calc_rx1day_top10_proportion(
  future_rx_array,
  rx1day_threshold_98.9_current
)

probability_ratio_rx1day_98.9 <- calc_probability_ratio(
  current_prop_rx1day_98.9,
  future_prop_rx1day_98.9
)

dim(current_prop_rx1day_98.9)
dim(future_prop_rx1day_98.9)
# both should be 44 X 44


# Matching Rx1day Extreme threshold with 3.4% -----------------------------

rx1day_threshold_96.6_current <- apply(
  current_rx_array,
  c(1, 2),
  quantile,
  probs = 0.966,
  na.rm = TRUE,
  type = 7
)

rx1day_threshold_96.6_future <- apply(
  future_rx_array,
  c(1, 2),
  quantile,
  probs = 0.966,
  na.rm = TRUE,
  type = 7
)

current_prop_rx1day_96.6 <- calc_rx1day_top10_proportion(
  current_rx_array,
  rx1day_threshold_96.6_current
)
future_prop_rx1day_96.6 <- calc_rx1day_top10_proportion(
  future_rx_array,
  rx1day_threshold_96.6_current
)

probability_ratio_rx1day_96.6 <- calc_probability_ratio(
  current_prop_rx1day_96.6,
  future_prop_rx1day_96.6
)

dim(current_prop_rx1day_96.6)
dim(future_prop_rx1day_96.6)
# both should be 44 X 44

# Joint proportions with >=4 heavy days -----------------------------------

current_prop_joint_98.9_ge4 <- calc_joint_top10_ge4_proportion(
  rx_array = current_rx_array,
  exceedance_array = current_exceedance_array,
  threshold_matrix = rx1day_threshold_98.9_current,
  min_days = 4
)
future_prop_joint_98.9_ge4 <- calc_joint_top10_ge4_proportion(
  rx_array = future_rx_array,
  exceedance_array = future_exceedance_array,
  threshold_matrix = rx1day_threshold_98.9_current,
  min_days = 4
)
probability_ratio_joint_98.9_ge4 <- calc_probability_ratio(
  current_prop_joint_98.9_ge4,
  future_prop_joint_98.9_ge4
)

current_prop_joint_96.6_ge4 <- calc_joint_top10_ge4_proportion(
  rx_array = current_rx_array,
  exceedance_array = current_exceedance_array,
  threshold_matrix = rx1day_threshold_96.6_current,
  min_days = 4
)
future_prop_joint_96.6_ge4 <- calc_joint_top10_ge4_proportion(
  rx_array = future_rx_array,
  exceedance_array = future_exceedance_array,
  threshold_matrix = rx1day_threshold_96.6_current,
  min_days = 4
)
probability_ratio_joint_96.6_ge4 <- calc_probability_ratio(
  current_prop_joint_96.6_ge4,
  future_prop_joint_96.6_ge4
)

# Saving alternative-threshold CSV --------------------------------

nc_grid <- open.nc(current_day_files[1])
longitude0 <- var.get.nc(nc_grid, "longitude0")
latitude0 <- var.get.nc(nc_grid, "latitude0")
global_longitude0 <- var.get.nc(nc_grid, "global_longitude0")
global_latitude0 <- var.get.nc(nc_grid, "global_latitude0")
close.nc(nc_grid)

grid_template <- expand.grid(
  lon_index = seq_along(longitude0),
  lat_index = seq_along(latitude0)
)

grid_results_design <- data.frame(
  lon_index = grid_template$lon_index,
  lat_index = grid_template$lat_index,
  longitude0 = longitude0[grid_template$lon_index],
  latitude0 = latitude0[grid_template$lat_index],
  global_longitude0 = global_longitude0[cbind(grid_template$lon_index, grid_template$lat_index)],
  global_latitude0 = global_latitude0[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_33_current = rx1day_threshold_33_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_98.9_current = rx1day_threshold_98.9_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_96.6_current = rx1day_threshold_96.6_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge4_current = current_prop_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge4_future = future_prop_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_rx1day_98.9_current = current_prop_rx1day_98.9[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_rx1day_98.9_future = future_prop_rx1day_98.9[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_rx1day_98.9_future_over_current = probability_ratio_rx1day_98.9[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_98.9_ge4_current = current_prop_joint_98.9_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_98.9_ge4_future = future_prop_joint_98.9_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_joint_98.9_ge4_future_over_current = probability_ratio_joint_98.9_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_rx1day_96.6_current = current_prop_rx1day_96.6[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_rx1day_96.6_future = future_prop_rx1day_96.6[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_rx1day_96.6_future_over_current = probability_ratio_rx1day_96.6[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_96.6_ge4_current = current_prop_joint_96.6_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_96.6_ge4_future = future_prop_joint_96.6_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_joint_96.6_ge4_future_over_current = probability_ratio_joint_96.6_ge4[cbind(grid_template$lon_index, grid_template$lat_index)]
)

weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"

design_csv_file <- file.path(
  weatherathome_dir,
  "weather@home_design_check_threshold_98.9_96.6_grid.csv"
)

write.csv(grid_results_design, design_csv_file, row.names = FALSE)

# NZ-only mapping using same style as 05 ----------------------------

model_data_dir <- weatherathome_dir
nc_file <- file.path(model_data_dir, "current_decade",
                     "item5216_daily_mean_a000_2006-07_2007-06-NZtrim-mm.nc")
lse_mask_file <- file.path(model_data_dir, "Land-Sea Mask for Weather@home Data.csv")
mask_transpose <- TRUE

load_lse_mask_matrix <- function(mask_file, nlon, nlat, transpose_mask = FALSE) {
  raw <- utils::read.csv(
    file = mask_file,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NaN", "NA", "")
  )
  
  numeric_raw <- suppressWarnings(as.data.frame(lapply(raw, as.numeric)))
  
  as_m <- function(x) as.matrix(x)
  base_candidates <- list(
    as_m(numeric_raw),
    as_m(numeric_raw[-1, , drop = FALSE]),
    as_m(numeric_raw[, -1, drop = FALSE]),
    as_m(numeric_raw[-1, -1, drop = FALSE])
  )
  
  dims_ok <- vapply(base_candidates, function(m) {
    is.matrix(m) && nrow(m) == nlon && ncol(m) == nlat
  }, logical(1))
  
  if (!any(dims_ok)) {
    stop(sprintf("Could not coerce LSE mask to %d x %d matrix.", nlon, nlat))
  }
  
  mask <- base_candidates[[which(dims_ok)[1]]]
  if (isTRUE(transpose_mask)) {
    mask <- t(mask)
  }
  
  if (nrow(mask) != nlon || ncol(mask) != nlat) {
    stop("Mask dimensions do not match rainfall grid after transpose setting.")
  }
  mask
}

build_cell_polygons <- function(lon_mat, lat_mat, value_mat, value_name = "value") {
  centres <- expand.grid(lon_index = seq_len(nrow(lon_mat)), lat_index = seq_len(ncol(lon_mat)))
  centres$lon <- as.vector(lon_mat)
  centres$lat <- as.vector(lat_mat)
  centres$value <- as.vector(value_mat)
  centres <- centres[is.finite(centres$lon) & is.finite(centres$lat), ]
  
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
    
    v_i <- if (all(is.finite(c_w)) && all(is.finite(c_e))) (c_e - c_w) / 2 else if (all(is.finite(c_e))) c_e - c0 else if (all(is.finite(c_w))) c0 - c_w else c(NA_real_, NA_real_)
    v_j <- if (all(is.finite(c_s)) && all(is.finite(c_n))) (c_n - c_s) / 2 else if (all(is.finite(c_n))) c_n - c0 else if (all(is.finite(c_s))) c0 - c_s else c(NA_real_, NA_real_)
    
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
      id = paste(i, j, sep = "_"),
      lon_index = i,
      lat_index = j,
      lon = corners[, 1],
      lat = corners[, 2],
      vertex_id = seq_len(nrow(corners)),
      value = ifelse(is.nan(centres$value[r]), NA_real_, centres$value[r])
    )
  }
  
  if (part_i == 0L) return(data.frame())
  polygons <- do.call(rbind, polygon_parts[seq_len(part_i)])
  names(polygons)[names(polygons) == "value"] <- value_name
  polygons
}

sanitize_geometry <- function(x) {
  if (nrow(x) == 0) return(x)
  is_bad <- !st_is_valid(x)
  if (any(is_bad)) x[is_bad, ] <- st_make_valid(x[is_bad, ])
  x
}

cell_polygons_to_sf <- function(cell_polygons_df) {
  if (nrow(cell_polygons_df) == 0) return(st_sf(id = character(0), geometry = st_sfc(crs = 4326)))
  split_polys <- split(cell_polygons_df, cell_polygons_df$id)
  sf_list <- lapply(names(split_polys), function(id) {
    piece <- split_polys[[id]]
    coords <- as.matrix(piece[order(piece$vertex_id), c("lon", "lat")])
    coords <- coords[stats::complete.cases(coords), , drop = FALSE]
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1, ] == coords[nrow(coords), ])) coords <- rbind(coords, coords[1, ])
    st_sf(id = id, geometry = st_sfc(st_polygon(list(coords)), crs = 4326))
  })
  sf_list <- Filter(Negate(is.null), sf_list)
  if (length(sf_list) == 0) return(st_sf(id = character(0), geometry = st_sfc(crs = 4326)))
  do.call(rbind, sf_list)
}

get_nz_intersecting_cell_ids <- function(cell_polygons_df) {
  if (nrow(cell_polygons_df) == 0) return(character(0))
  old_s2 <- sf_use_s2(); on.exit(sf_use_s2(old_s2), add = TRUE); sf_use_s2(FALSE)
  nz_map <- maps::map("nz", fill = TRUE, plot = FALSE)
  nz_sf <- st_as_sf(nz_map) |> st_set_crs(4326) |> sanitize_geometry()
  nz_union <- st_union(nz_sf)
  nz_union <- st_sf(geometry = nz_union) |> sanitize_geometry()
  cells_sf <- cell_polygons_to_sf(cell_polygons_df) |> sanitize_geometry()
  if (nrow(cells_sf) == 0 || nrow(nz_union) == 0) return(character(0))
  intersects <- st_intersects(cells_sf, nz_union, sparse = FALSE)[, 1]
  cells_sf$id[intersects]
}

make_nz_ratio_plot <- function(poly_df, keep_ids, title_text, ratio_breaks, ratio_palette, output_file) {
  poly_nz <- poly_df[poly_df$id %in% keep_ids, ]
  poly_nz <- poly_nz[is.finite(poly_nz$ratio_value), ]
  plot_breaks <- c(-Inf, ratio_breaks, Inf)
  bin_levels <- paste0("bin_", seq_len(length(plot_breaks) - 1))
  poly_nz$ratio_bin <- cut(poly_nz$ratio_value, breaks = plot_breaks, include.lowest = TRUE, right = FALSE, labels = bin_levels)
  poly_nz$ratio_bin <- factor(poly_nz$ratio_bin, levels = bin_levels)
  
  nz_outline <- map_data("nz")
  p <- ggplot(poly_nz, aes(x = lon, y = lat, fill = ratio_bin)) +
    geom_polygon(aes(group = id), colour = NA, linewidth = 0) +
    geom_path(data = nz_outline, aes(x = long, y = lat, group = group), inherit.aes = FALSE, colour = "black", linewidth = 0.45, alpha = 0.9) +
    coord_fixed() +
    scale_fill_manual(values = setNames(ratio_palette, bin_levels), drop = FALSE, guide = "none") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14), axis.title = element_blank())
  
  ggsave(output_file, p, width = 6.7, height = 8.0, dpi = 600)
}

nc <- open.nc(nc_file)
on.exit(close.nc(nc), add = TRUE)
lon <- var.get.nc(nc, "global_longitude0")
lat <- var.get.nc(nc, "global_latitude0")
rain <- var.get.nc(nc, "item5216_daily_mean")[, , 1]
if (length(dim(lon)) == 3) lon <- lon[, , 1]
if (length(dim(lat)) == 3) lat <- lat[, , 1]
nlon <- dim(rain)[1]
nlat <- dim(rain)[2]

mask_matrix <- load_lse_mask_matrix(lse_mask_file, nlon, nlat, mask_transpose)
mask_is_land <- !is.na(mask_matrix)

ratio_breaks <- c(1, 1.5, 2, 2.5, 3, 4, 5)
ratio_palette <- c("#D0D4DA", "#D7E8FF", "#BFD9FF", "#7FB3FF", "#3F8BE6", "#0B4FAF", "#08306B", "#041F4A")

plot_ratio_surface <- function(df, ratio_col, title_text, output_file) {
  ratio_mat <- matrix(NA_real_, nrow = nlon, ncol = nlat)
  idx <- cbind(as.integer(df$lon_index), as.integer(df$lat_index))
  ratio_mat[idx] <- as.numeric(df[[ratio_col]])
  ratio_mat[!mask_is_land] <- NA_real_
  
  ratio_poly <- build_cell_polygons(lon, lat, ratio_mat, "ratio_value")
  keep_ids <- get_nz_intersecting_cell_ids(ratio_poly)
  make_nz_ratio_plot(ratio_poly, keep_ids, title_text, ratio_breaks, ratio_palette, output_file)
  list(poly = ratio_poly, keep_ids = keep_ids)
}

ratio_layers_design <- list()

ratio_layers_design[["probability_ratio_rx1day_98.9_future_over_current"]] <- plot_ratio_surface(
  df = grid_results_design,
  ratio_col = "probability_ratio_rx1day_98.9_future_over_current",
  title_text = "(a) Years with Extreme Rx1day (98.9th)",
  output_file = file.path(weatherathome_dir, "weather@home_design_check_ratio_map_98.9.png")
)

ratio_layers_design[["probability_ratio_rx1day_96.6_future_over_current"]] <- plot_ratio_surface(
  df = grid_results_design,
  ratio_col = "probability_ratio_rx1day_96.6_future_over_current",
  title_text = "(a) Years with Extreme Rx1day (96.6th)",
  output_file = file.path(weatherathome_dir, "weather@home_design_check_ratio_map_96.6.png")
)

calc_finite_stats <- function(x) {
  x_finite <- x[is.finite(x)]
  if (length(x_finite) == 0) {
    return(c(min = NA_real_, mean = NA_real_, max = NA_real_))
  }
  c(min = min(x_finite), mean = mean(x_finite), max = max(x_finite))
}

get_mapped_values <- function(layer_entry) {
  keep_rows <- layer_entry$poly$id %in% layer_entry$keep_ids
  layer_entry$poly$ratio_value[keep_rows & is.finite(layer_entry$poly$ratio_value)]
}

nz_probability_ratio_rx1day_98.9_stats <- calc_finite_stats(
  get_mapped_values(ratio_layers_design[["probability_ratio_rx1day_98.9_future_over_current"]]))
nz_probability_ratio_rx1day_96.6_stats <- calc_finite_stats(
  get_mapped_values(ratio_layers_design[["probability_ratio_rx1day_96.6_future_over_current"]]))

summary_stats_design <- data.frame(
  metric = c(
    "probability_ratio_rx1day_98.9_future_over_current",
    "probability_ratio_rx1day_96.6_future_over_current"
  ),
  min = c(nz_probability_ratio_rx1day_98.9_stats[["min"]], nz_probability_ratio_rx1day_96.6_stats[["min"]]),
  mean = c(nz_probability_ratio_rx1day_98.9_stats[["mean"]], nz_probability_ratio_rx1day_96.6_stats[["mean"]]),
  max = c(nz_probability_ratio_rx1day_98.9_stats[["max"]], nz_probability_ratio_rx1day_96.6_stats[["max"]])
)

write.csv(
  summary_stats_design,
  file = file.path(weatherathome_dir, "weather@home_design_check_summary_stats.csv"),
  row.names = FALSE
)

nz_probability_ratio_rx1day_98.9_stats
nz_probability_ratio_rx1day_96.6_stats
summary_stats_design


# Heatmap Style -----------------------------------------------------------


if (!exists("matched_cell")) {
  matched_cell <- data.frame(lon_index = 30L, lat_index = 16L)
}

if (!exists("box_colour_light")) box_colour_light <- "lavender"
if (!exists("box_colour_dark")) box_colour_dark <- "lightslateblue"
if (!exists("theme_thesis")) {
  theme_thesis <- ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 7),
      panel.grid.major = ggplot2::element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(colour = "grey92", linewidth = 0.2)
    )
}

if (!exists("build_quadrant_df")) {
  build_quadrant_df <- function(df_period, heavy_cutoff = 4L) {
    exceed_max <- max(df_period$heavy_days, na.rm = TRUE)
    
    tile_df <- df_period %>%
      dplyr::mutate(
        exceed_group = dplyr::if_else(
          heavy_days < heavy_cutoff,
          paste0("<", heavy_cutoff, " heavy days"),
          paste0(">=", heavy_cutoff, " heavy days")
        ),
        rx_group = dplyr::if_else(RX1day < extreme_threshold, "below Extreme", "at/above Extreme")
      ) %>%
      dplyr::count(exceed_group, rx_group, name = "n_years") %>%
      tidyr::complete(
        exceed_group = c(paste0("<", heavy_cutoff, " heavy days"), paste0(">=", heavy_cutoff, " heavy days")),
        rx_group = c("below Extreme", "at/above Extreme"),
        fill = list(n_years = 0)
      ) %>%
      dplyr::mutate(
        pct_years = 100 * n_years / sum(n_years),
        pct_label = paste0(round(pct_years, 1), "%"),
        xmin = dplyr::if_else(exceed_group == paste0("<", heavy_cutoff, " heavy days"), -0.5, heavy_cutoff),
        xmax = dplyr::if_else(exceed_group == paste0("<", heavy_cutoff, " heavy days"), heavy_cutoff, exceed_max + 0.5),
        ymin = dplyr::if_else(rx_group == "below Extreme", 0, unique(df_period$extreme_threshold)),
        ymax = dplyr::if_else(rx_group == "below Extreme", unique(df_period$extreme_threshold), max(df_period$RX1day, na.rm = TRUE) * 1.02),
        xmid = (xmin + xmax) / 2,
        ymid = (ymin + ymax) / 2
      )
    
    list(tile_df = tile_df, exceed_max = exceed_max)
  }
}

if (!exists("plot_quadrant_heatmap")) {
  plot_quadrant_heatmap <- function(df_period, panel_tag, heavy_cutoff = 4L, x_max = NULL) {
    hm <- build_quadrant_df(df_period, heavy_cutoff = heavy_cutoff)
    if (is.null(x_max)) {
      x_max <- hm$exceed_max
    }
    
    tile_df_plot <- hm$tile_df %>%
      dplyr::mutate(
        xmax = dplyr::if_else(exceed_group == paste0(">=", heavy_cutoff, " heavy days"), x_max + 0.5, xmax),
        xmid = (xmin + xmax) / 2
      )
    
    ggplot2::ggplot(tile_df_plot) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pct_years),
        colour = "white",
        linewidth = 0.3
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = xmid, y = ymid, label = pct_label),
        colour = "black",
        size = 3.0,
        fontface = "bold"
      ) +
      ggplot2::geom_vline(xintercept = heavy_cutoff, colour = "black", linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = unique(df_period$extreme_threshold), colour = "black", linewidth = 0.35) +
      ggplot2::scale_x_continuous(breaks = 0:x_max, limits = c(-0.5, x_max + 0.5), expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_y_continuous(limits = c(0, max(df_period$RX1day, na.rm = TRUE) * 1.02), expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_fill_gradient(low = box_colour_light, high = box_colour_dark, limits = c(0, 100), guide = "none") +
      ggplot2::labs(title = panel_tag, x = "Number of Heavy days", y = "RX1day (mm)") +
      theme_thesis +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "grey45", fill = NA, linewidth = 0.35),
        axis.line = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0)
      )
  }
}

if (!exists("make_column_header")) {
  make_column_header <- function(label_text) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = label_text, fontface = "bold", size = 3.5) +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
  }
}

build_design_cell_df <- function(rx_array, exceedance_array, extreme_threshold, period_label) {
  lon_index <- matched_cell$lon_index
  lat_index <- matched_cell$lat_index
  n_years <- dim(rx_array)[3]
  data.frame(
    Year = seq_len(n_years),
    RX1day = as.numeric(rx_array[lon_index, lat_index, ]),
    Period = period_label,
    heavy_threshold = as.numeric(rx1day_threshold_33_current[lon_index, lat_index]),
    extreme_threshold = as.numeric(extreme_threshold),
    heavy_days = as.numeric(exceedance_array[lon_index, lat_index, ]),
    stringsAsFactors = FALSE
  )
}

build_alternative_threshold_heatmap_pair <- function(extreme_threshold, threshold_label, output_file) {
  current_df <- build_design_cell_df(
    rx_array = current_rx_array,
    exceedance_array = current_exceedance_array,
    extreme_threshold = extreme_threshold,
    period_label = "Current Day"
  )
  future_df <- build_design_cell_df(
    rx_array = future_rx_array,
    exceedance_array = future_exceedance_array,
    extreme_threshold = extreme_threshold,
    period_label = "Future Projection"
  )
  heavy_days_max <- max(c(current_df$heavy_days, future_df$heavy_days), na.rm = TRUE)
  
  current_panel <- plot_quadrant_heatmap(current_df, panel_tag = "(e)", heavy_cutoff = 4L, x_max = heavy_days_max)
  future_panel <- plot_quadrant_heatmap(future_df, panel_tag = "(f)", heavy_cutoff = 4L, x_max = heavy_days_max)
  
  current_column <- patchwork::wrap_plots(
    make_column_header("Current Day"),
    current_panel,
    ncol = 1,
    heights = c(0.09, 1)
  )
  future_column <- patchwork::wrap_plots(
    make_column_header("Future Projection"),
    future_panel,
    ncol = 1,
    heights = c(0.09, 1)
  )
  heatmap_plot <- patchwork::wrap_plots(current_column, future_column, ncol = 2)
  heatmap_plot <- heatmap_plot + patchwork::plot_annotation(
    title = paste0(threshold_label)
  )
  
  ggplot2::ggsave(filename = output_file, plot = heatmap_plot, width = 7, height = 2.8, dpi = 2000)
  heatmap_plot
}

matched_cell_heatmap_98.9_file <- file.path(
  weatherathome_dir,
  "weather@home_design_check_matched_cell_heatmap_98.9.png"
)
matched_cell_heatmap_96.6_file <- file.path(
  weatherathome_dir,
  "weather@home_design_check_matched_cell_heatmap_96.6.png"
)

matched_cell_heatmap_98.9 <- build_alternative_threshold_heatmap_pair(
  extreme_threshold = rx1day_threshold_98.9_current[matched_cell$lon_index, matched_cell$lat_index],
  threshold_label = "(98.9th)",
  output_file = matched_cell_heatmap_98.9_file
)

matched_cell_heatmap_96.6 <- build_alternative_threshold_heatmap_pair(
  extreme_threshold = rx1day_threshold_96.6_current[matched_cell$lon_index, matched_cell$lat_index],
  threshold_label = "(96.6th)",
  output_file = matched_cell_heatmap_96.6_file
)

#Matched grid-cell RX1day threshold comparison ---------------------------
  
  if (!exists("rx1day_threshold_90_current")) {
    rx1day_threshold_90_current <- apply(
      current_rx_array,
      c(1, 2),
      quantile,
      probs = 0.90,
      na.rm = TRUE,
      type = 7
    )
  }

if (!exists("rx1day_threshold_90_future")) {
  rx1day_threshold_90_future <- apply(
    future_rx_array,
    c(1, 2),
    quantile,
    probs = 0.90,
    na.rm = TRUE,
    type = 7
  )
}

matched_cell_rx1day_threshold_comparison <- data.frame(
  lon_index = matched_cell$lon_index,
  lat_index = matched_cell$lat_index,
  percentile = c("90th", "96.6th", "98.9th"),
  probability = c(0.900, 0.966, 0.989),
  current_threshold_mm = c(
    rx1day_threshold_90_current[matched_cell$lon_index, matched_cell$lat_index],
    rx1day_threshold_96.6_current[matched_cell$lon_index, matched_cell$lat_index],
    rx1day_threshold_98.9_current[matched_cell$lon_index, matched_cell$lat_index]
  ),
  future_threshold_mm = c(
    rx1day_threshold_90_future[matched_cell$lon_index, matched_cell$lat_index],
    rx1day_threshold_96.6_future[matched_cell$lon_index, matched_cell$lat_index],
    rx1day_threshold_98.9_future[matched_cell$lon_index, matched_cell$lat_index]
  )
)

matched_cell_rx1day_threshold_comparison$current_threshold_mm <- round(
  matched_cell_rx1day_threshold_comparison$current_threshold_mm,
  2
)
matched_cell_rx1day_threshold_comparison$future_threshold_mm <- round(
  matched_cell_rx1day_threshold_comparison$future_threshold_mm,
  2
)

matched_cell_rx1day_threshold_comparison_file <- file.path(
  weatherathome_dir,
  "weather@home_design_check_matched_cell_rx1day_threshold_comparison.csv"
)

write.csv(
  matched_cell_rx1day_threshold_comparison,
  matched_cell_rx1day_threshold_comparison_file,
  row.names = FALSE
)

matched_cell_rx1day_threshold_comparison


