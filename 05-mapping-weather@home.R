# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Build NZ-intersection weather@home probability-ratio maps.
# Outputs:
# - main combined map (>=4, top 10% Rx1day, and joint)
# - >=5 exceedance map (for supplementary material)
# -------------------------------------------------------------------------

library(ggplot2)
library(viridis)
library(maps)
library(dplyr)
library(patchwork)
library(sf)
library(grid)
library(colorspace)

# Colouring  --------------------------------------------------------------
#blue_palette <- choose_palette()
blues <- blue_palette(14)

#alternative colouring replace "blues" with "ratio_palette"
#ratio_palette <- viridisLite::viridis(length(ratio_breaks) - 1, option = "viridis", direction = -1)

# Input/output paths -------------------------------------------------------
weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"
# Reading processed data ------------------------------------------------------
grid_results <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv")

# Shared helpers -----------------------------------------------------------
get_fixed_width_bin_spec <- function(x, bin_width = 0.5, min_value = 0, max_value = NULL) {
  r <- range(x, na.rm = TRUE)
  min_break <- min_value
  max_break <- if (is.null(max_value)) ceiling(r[2] / bin_width) * bin_width else max_value
  
  if (min_break == max_break) {
    max_break <- min_break + bin_width}
  
  brks <- seq(min_break, max_break, by = bin_width)
  
  list(
    breaks = brks,
    labels = format(head(brks, -1), trim = TRUE, scientific = FALSE, nsmall = 1))}

build_discrete_ratio_bins <- function(x, bin_spec) {
  cut(x, breaks = bin_spec$breaks, include.lowest = TRUE,
    right = TRUE, labels = bin_spec$labels, ordered_result = TRUE)}

build_cell_polygons <- function(df) {
  has_ratio_bin <- "ratio_bin" %in% names(df)
  centres <- df |>select(lon_index, lat_index, lon, lat, ratio_value) |>
    distinct()
  
  if (has_ratio_bin) {
    ratio_bins <- df |>
      select(lon_index, lat_index, ratio_bin) |>
      distinct()
    
    centres <- centres |>
      left_join(ratio_bins, by = c("lon_index", "lat_index"))} else {
    centres$ratio_bin <- NA_character_}
  
  key <- paste(centres$lon_index, centres$lat_index, sep = "_")
  key_lookup <- setNames(seq_len(nrow(centres)), key)
  
  get_xy <- function(i, j) {
    k <- paste(i, j, sep = "_")
    idx <- unname(key_lookup[k])
    if (length(idx) == 0L || is.na(idx)) return(c(NA_real_, NA_real_))
    c(centres$lon[idx], centres$lat[idx])}
  
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
      c(NA_real_, NA_real_)}
    
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
      c0 - 0.5 * v_i - 0.5 * v_j)
    
    part_i <- part_i + 1L
    polygon_parts[[part_i]] <- data.frame(
      cell_id = paste(i, j, sep = "_"),
      vertex_id = seq_len(5),
      lon = corners[, 1],
      lat = corners[, 2],
      ratio_value = centres$ratio_value[r],
      ratio_bin = centres$ratio_bin[r])}
  
  if (part_i == 0L) return(data.frame())
  do.call(rbind, polygon_parts[seq_len(part_i)])}

cell_polygons_to_sf <- function(cell_polygons_df) {
  if (nrow(cell_polygons_df) == 0) {
    return(st_sf(cell_id = character(0), geometry = st_sfc(crs = 4326)))}
  
  split_polys <- split(cell_polygons_df, cell_polygons_df$cell_id)
  
  sf_list <- lapply(names(split_polys), function(id) {
    piece <- split_polys[[id]]
    coords <- as.matrix(piece[order(piece$vertex_id), c("lon", "lat")])
    coords <- coords[stats::complete.cases(coords), , drop = FALSE]
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1, ] == coords[nrow(coords), ])) coords <- rbind(coords, coords[1, ])
    
    st_sf(cell_id = id, geometry = st_sfc(st_polygon(list(coords)), crs = 4326))})
  
  sf_list <- Filter(Negate(is.null), sf_list)
  if (length(sf_list) == 0) {
    return(st_sf(cell_id = character(0), geometry = st_sfc(crs = 4326)))}
  
  do.call(rbind, sf_list)}

sanitize_geometry <- function(x) {
  x <- suppressWarnings(st_make_valid(x))
  x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  non_empty <- !st_is_empty(x)
  non_empty[is.na(non_empty)] <- FALSE
  x[non_empty, ]}

get_nz_intersecting_cell_ids <- function(cell_polygons_df) {
  if (nrow(cell_polygons_df) == 0) return(character(0))
  
  old_s2 <- sf_use_s2()
  on.exit(sf_use_s2(old_s2), add = TRUE)
  sf_use_s2(FALSE)
  
  nz_map <- maps::map("nz", fill = TRUE, plot = FALSE)
  nz_sf <- st_as_sf(nz_map) |> st_set_crs(4326) |> sanitize_geometry()
  nz_union <- st_union(nz_sf)
  nz_union <- st_sf(geometry = nz_union) |> sanitize_geometry()
  
  cells_sf <- cell_polygons_to_sf(cell_polygons_df) |> sanitize_geometry()
  if (nrow(cells_sf) == 0 || nrow(nz_union) == 0) return(character(0))
  
  intersects <- st_intersects(cells_sf, nz_union, sparse = FALSE)[, 1]
  cells_sf$cell_id[intersects]}

build_metric_layer <- function(ratio_col) {
  base_data <- data.frame(
    lon_index = grid_results$lon_index,
    lat_index = grid_results$lat_index,
    lon = grid_results$global_longitude0,
    lat = grid_results$global_latitude0,
    ratio_value = grid_results[[ratio_col]])
  
  finite_data <- base_data[is.finite(base_data$ratio_value), ]
  if (nrow(finite_data) == 0) {
    stop(sprintf("No finite values to plot for %s.", ratio_col))}
  
  cell_polygons <- build_cell_polygons(finite_data)
  keep_cell_ids <- get_nz_intersecting_cell_ids(cell_polygons)
  
  list(
    finite_data = finite_data,
    cell_polygons = cell_polygons,
    keep_cell_ids = keep_cell_ids)}

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
    stop(sprintf("No NZ-intersecting cells found for '%s'.", title_text))}
  
  ratio_limits <- range(ratio_breaks)
  ratio_labels <- format(ratio_breaks, trim = TRUE, scientific = FALSE, nsmall = 1)
  nz_outline <- map_data("nz")
  
  legend_position <- if (isTRUE(show_legend)) "right" else "none"
  legend_height <- unit(legend_height_cm, "cm")
  legend_width <- unit(legend_width_cm, "cm")
  
  ggplot(cell_polygons_nz, aes(x = lon, y = lat, fill = ratio_value)) +
    geom_polygon(aes(group = cell_id), colour = NA, linewidth = 0) +
    geom_path(data = nz_outline, aes(x = long, y = lat, group = group),
    inherit.aes = FALSE, colour = "white", linewidth = 0.45, alpha = 0.9) +
    coord_fixed() +
    scale_fill_stepsn(
      colours = blues, breaks = ratio_breaks, 
      labels = ratio_labels, limits = ratio_limits,
      oob = scales::squish, show.limits = TRUE,
      name = "Probability Ratio",
      guide = guide_coloursteps(reverse = FALSE,
      even.steps = TRUE, show.limits = TRUE,
      barheight = legend_height,
      barwidth = legend_width,
      title.position = "top", title.hjust = 0.5)) +
    guides(
      fill = guide_coloursteps(reverse = FALSE, even.steps = TRUE,
      show.limits = TRUE, barheight = legend_height,
       barwidth = legend_width, title.position = "top",
       title.hjust = 0.5)) +
    labs(
      title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 9),
      legend.key.height = unit(1.05, "cm"),
      legend.key.width = legend_width,
      legend.margin = margin(4, 8, 4, 2),
      legend.position = legend_position)}

# Build layers --------------------------------------------------------------
layers_ge4 <- build_metric_layer("probability_ratio_ge4_future_over_current")
layers_ge5 <- build_metric_layer("probability_ratio_ge5_future_over_current")
layers_top10 <- build_metric_layer("probability_ratio_rx1day_top10_future_over_current")
layers_joint <- build_metric_layer("probability_ratio_joint_top10_ge4_future_over_current")

# Shared bins using NZ-intersecting cells only -----------------------------
get_intersection_values <- function(layer_obj) {
  keep_rows <- paste(layer_obj$finite_data$lon_index, layer_obj$finite_data$lat_index, sep = "_") %in% layer_obj$keep_cell_ids
  layer_obj$finite_data$ratio_value[keep_rows]}

intersection_values <- c(
  get_intersection_values(layers_ge4),
  get_intersection_values(layers_ge5),
  get_intersection_values(layers_top10),
  get_intersection_values(layers_joint))

intersection_values <- intersection_values[is.finite(intersection_values)]

if (length(intersection_values) == 0) {
  stop("No finite probability-ratio values found for NZ-intersecting cells.")}

ratio_bin_spec <- get_fixed_width_bin_spec(intersection_values, bin_width = 0.5, min_value = 0)
ratio_breaks <- ratio_bin_spec$breaks

# Attach ratio bins ---------------------------------------------------------
add_ratio_bins <- function(layer_obj) {
  layer_obj$finite_data$ratio_bin <- build_discrete_ratio_bins(layer_obj$finite_data$ratio_value, ratio_bin_spec)
  layer_obj$cell_polygons$ratio_bin <- build_discrete_ratio_bins(layer_obj$cell_polygons$ratio_value, ratio_bin_spec)
  layer_obj}

layers_ge4 <- add_ratio_bins(layers_ge4)
layers_ge5 <- add_ratio_bins(layers_ge5)
layers_top10 <- add_ratio_bins(layers_top10)
layers_joint <- add_ratio_bins(layers_joint)

# Build requested plots -----------------------------------------------------
p_ge4 <- make_nz_ratio_plot(
  layers_ge4, "(b) Years with ≥4 exceedances",
  ratio_breaks, blues)

p_top10 <- make_nz_ratio_plot(
  layers_top10, "(a) Years with extreme Rx1day",
  ratio_breaks, blues, show_legend = FALSE)

p_joint <- make_nz_ratio_plot(
  layers_joint, "(c) Years with extreme Rx1day AND ≥4 exceedances",
  ratio_breaks,blues, show_legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

p_ge5 <- make_nz_ratio_plot(
  layers_ge5,"Years with ≥5 exceedances",
  ratio_breaks, blues, legend_height_cm = 13)

# Top row with ≥4 exceedances and extreme Rx1day, bottom row the joint.
combined_design <- c(
  area(t = 1, l = 1, b = 1, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 2, l = 1, b = 2, r = 2),
  area(t = 1, l = 3, b = 2, r = 3))

p_combined <- (p_top10 + p_ge4 + p_joint + guide_area()) +
  plot_layout(
    design = combined_design,
    guides = "collect",
    widths = c(1, 1, 0.18),
    heights = c(1, 1)) +
  plot_annotation(
    theme = theme(
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(5, 5, 5, 5)))

p_combined
p_ge5

# Save outputs --------------------------------------------------------------
output_map_data <- grid_results |>
  select(
    lon_index, lat_index,
    global_longitude0, global_latitude0,
    probability_ratio_ge4_future_over_current,
    probability_ratio_ge5_future_over_current,
    probability_ratio_rx1day_top10_future_over_current,
    probability_ratio_joint_top10_ge4_future_over_current)

write.csv(
  output_map_data,
  file.path(weatherathome_dir, "weather@home_probability_ratio_map_data.csv"),
  row.names = FALSE)

ggsave(
  filename = file.path(weatherathome_dir, "weather@home_probability_ratio_ge4_top10_joint_nz_intersection_combined_map.png"),
  plot = p_combined, width = 12, height = 10, dpi = 300)

ggsave(
  filename = file.path(weatherathome_dir, "weather@home_probability_ratio_ge5_nz_intersection_map.png"),
  plot = p_ge5, width = 8, height = 7, dpi = 300)

# Quick checks --------------------------------------------------------------
cat("NZ-intersecting cell count (>=4):", length(layers_ge4$keep_cell_ids), "\n")
cat("NZ-intersecting cell count (>=5):", length(layers_ge5$keep_cell_ids), "\n")
cat("NZ-intersecting cell count (top 10%):", length(layers_top10$keep_cell_ids), "\n")
cat("NZ-intersecting cell count (joint):", length(layers_joint$keep_cell_ids), "\n")
