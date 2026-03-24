# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Mapping weather@home data and saving NZ probability-ratio maps.
# -------------------------------------------------------------------------

# Loading packages ---------------------------------------------------------
library(RNetCDF)
library(ncdf4)
library(ggplot2)
library(maps)
library(dplyr)
library(patchwork)
library(sf)
library(grid)

# Setting input and output paths -------------------------------------------
model_data_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"
nc_file <- file.path(model_data_dir,"current_decade",
                     "item5216_daily_mean_a000_2006-07_2007-06-NZtrim-mm.nc")
lse_mask_file <- file.path(model_data_dir, "Land-Sea Mask for Weather@home Data.csv")
ratio_grid_file <- file.path(model_data_dir, "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv")

# Transposing the land mask
# If it looks weird and the wrong way around can can from True to False 
mask_transpose <- TRUE

masked_nc_file <- file.path(model_data_dir, "current_decade",
                            "item5216_daily_mean_a000_2006-07_2007-06-NZtrim-mm-masked.nc")

# Setting ratio map output names -------------------------------------------
combined_ratio_output_png <- file.path(model_data_dir,
                                       "weather@home_probability_ratio_ge4_top10_joint_combined_map.png")

ge5_ratio_output_png <- file.path(model_data_dir,
                                  "weather@home_probability_ratio_ge5_map.png")

# Loading the mask into an nlon x nlat matrix ------------------------------
load_lse_mask_matrix <- function(mask_file, nlon, nlat, transpose_mask = FALSE) {
  raw <- utils::read.csv(
    file = mask_file,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NaN", "NA", ""))
  
  numeric_raw <- suppressWarnings(as.data.frame(lapply(raw, as.numeric)))
  
  as_m <- function(x) as.matrix(x)
  base_candidates <- list(
    as_m(numeric_raw),
    as_m(numeric_raw[-1, , drop = FALSE]),
    as_m(numeric_raw[, -1, drop = FALSE]),
    as_m(numeric_raw[-1, -1, drop = FALSE]))
  
  dims_ok <- vapply(base_candidates, function(m) {
    is.matrix(m) && nrow(m) == nlon && ncol(m) == nlat
  }, logical(1))
  
  if (!any(dims_ok)) {
    stop(sprintf("Could not coerce LSE mask to %d x %d matrix.", nlon, nlat))}
  
  mask <- base_candidates[[which(dims_ok)[1]]]
  if (isTRUE(transpose_mask)) {
    mask <- t(mask)}
  
  if (nrow(mask) != nlon || ncol(mask) != nlat) {
    stop("Mask dimensions do not match rainfall grid after transpose setting.")}
  
  message(sprintf(
    "Loaded LSE mask as %d x %d matrix (transpose=%s).",
    nrow(mask),
    ncol(mask),
    ifelse(isTRUE(transpose_mask), "TRUE", "FALSE")))
  mask}

# Building touching cell polygons from the grid ----------------------------
build_cell_polygons <- function(lon_mat, lat_mat, value_mat, value_name = "value") {
  centres <- expand.grid(
    lon_index = seq_len(nrow(lon_mat)),
    lat_index = seq_len(ncol(lon_mat)))
  
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
      c0 - 0.5 * v_i - 0.5 * v_j)
    
    part_i <- part_i + 1L
    polygon_parts[[part_i]] <- data.frame(
      id = paste(i, j, sep = "_"),
      lon_index = i,
      lat_index = j,
      lon = corners[, 1],
      lat = corners[, 2],
      vertex_id = seq_len(nrow(corners)),
      value = ifelse(is.nan(centres$value[r]), NA_real_, centres$value[r]))}
  
  if (part_i == 0L) return(data.frame())
  polygons <- do.call(rbind, polygon_parts[seq_len(part_i)])
  names(polygons)[names(polygons) == "value"] <- value_name
  polygons}

# Placing indexed values into an nlon x nlat matrix ------------------------
matrix_from_indexed_values <- function(df, value_col, nlon, nlat) {
  out <- matrix(NA_real_, nrow = nlon, ncol = nlat)
  
  needed <- c("lon_index", "lat_index", value_col)
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns in ratio grid CSV: %s", paste(missing, collapse = ", ")))}
  
  idx_ok <- is.finite(df$lon_index) & is.finite(df$lat_index)
  idx_ok <- idx_ok & df$lon_index >= 1 & df$lon_index <= nlon
  idx_ok <- idx_ok & df$lat_index >= 1 & df$lat_index <= nlat
  
  if (!any(idx_ok)) {
    stop(sprintf("No valid lon_index/lat_index rows found for %s.", value_col))}
  
  mat_idx <- cbind(as.integer(df$lon_index[idx_ok]), as.integer(df$lat_index[idx_ok]))
  out[mat_idx] <- as.numeric(df[[value_col]][idx_ok])
  out}

# Keeping NZ-intersection helpers -------------------------------------------
sanitize_geometry <- function(x) {
  if (nrow(x) == 0) return(x)
  is_bad <- !st_is_valid(x)
  if (any(is_bad)) {
    x[is_bad, ] <- st_make_valid(x[is_bad, ])}
  x}

cell_polygons_to_sf <- function(cell_polygons_df) {
  if (nrow(cell_polygons_df) == 0) {
    return(st_sf(id = character(0), geometry = st_sfc(crs = 4326)))}
  
  split_polys <- split(cell_polygons_df, cell_polygons_df$id)
  
  sf_list <- lapply(names(split_polys), function(id) {
    piece <- split_polys[[id]]
    coords <- as.matrix(piece[order(piece$vertex_id), c("lon", "lat")])
    coords <- coords[stats::complete.cases(coords), , drop = FALSE]
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])}
    
    st_sf(id = id, geometry = st_sfc(st_polygon(list(coords)), crs = 4326))})
  
  sf_list <- Filter(Negate(is.null), sf_list)
  if (length(sf_list) == 0) {
    return(st_sf(id = character(0), geometry = st_sfc(crs = 4326)))}
  
  do.call(rbind, sf_list)}

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
  cells_sf$id[intersects]}

get_fixed_width_bin_spec <- function(x, bin_width = 1, min_value = 1, max_value = 6.5) {
  r <- range(x, na.rm = TRUE)
  min_break <- min_value
  max_break <- if (is.null(max_value)) ceiling(r[2] / bin_width) * bin_width else max_value
  
  if (min_break == max_break) {
    max_break <- min_break + bin_width}
  
  brks <- seq(min_break, max_break, by = bin_width)
  list(breaks = brks)}

# Making NZ probability-ratio plots ----------------------------------------
make_nz_ratio_plot <- function(poly_df, keep_ids, title_text, ratio_breaks, ratio_palette) {
  poly_nz <- poly_df[poly_df$id %in% keep_ids, ]
  poly_nz <- poly_nz[is.finite(poly_nz$ratio_value), ]
  
  if (nrow(poly_nz) == 0) {
    stop(sprintf("No NZ-intersecting cells found for '%s'.", title_text))}
  
  core_breaks <- ratio_breaks[ratio_breaks >= 1 & ratio_breaks <= 6]
  if (length(core_breaks) < 2) {
    stop("`ratio_breaks` must include at least two values between 1 and 6.")}
  
  plot_breaks <- c(-Inf, core_breaks, Inf)
  bin_levels <- paste0("bin_", seq_len(length(plot_breaks) - 1))
  
  poly_nz$ratio_bin <- cut(
    poly_nz$ratio_value,
    breaks = plot_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = bin_levels)
  poly_nz$ratio_bin <- factor(poly_nz$ratio_bin, levels = bin_levels)
  
  palette_for_bins <- setNames(ratio_palette, bin_levels)
  nz_outline <- map_data("nz")
  
  ggplot(poly_nz, aes(x = lon, y = lat, fill = ratio_bin)) +
    geom_polygon(aes(group = id), colour = NA, linewidth = 0) +
    geom_path(
      data = nz_outline,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.45,
      alpha = 0.9) +
    coord_fixed() +
    scale_fill_manual(
      values = palette_for_bins,
      drop = FALSE,
      guide = "none") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_blank(),
      legend.position = "none")}

make_triangle_colorbar_plot <- function(ratio_breaks, ratio_palette, legend_title = "Probability Ratio") {
  if (length(ratio_breaks) < 2) {
    stop("`ratio_breaks` must contain at least two values.")}
  
  core_breaks <- sort(unique(ratio_breaks))
  if (length(core_breaks) < 2) {
    stop("`ratio_breaks` must include at least two finite values.")}
  
  lower_cap <- 0
  upper_step <- diff(tail(core_breaks, 2))
  upper_cap <- max(core_breaks)
  
  interval_min <- c(lower_cap, head(core_breaks, -1))
  interval_max <- core_breaks
  
  if (length(ratio_palette) < (length(interval_min) + 1)) {
    stop("`ratio_palette` must include one colour per bar segment plus one for the top arrow head.")}
  
  bar_df <- data.frame(
    ymin = interval_min,
    ymax = interval_max,
    fill_col = ratio_palette[seq_len(length(interval_min))])
  
  ratio_min <- lower_cap
  ratio_max <- upper_cap
  ratio_span <- ratio_max - ratio_min
  top_triangle_height <- upper_step
  bar_xmin <- 0.31
  bar_xmax <- 0.57
  
  tri_df <- data.frame(
    x = c(bar_xmin, (bar_xmin + bar_xmax) / 2, bar_xmax),
    y = c(ratio_max, ratio_max + top_triangle_height, ratio_max),
    group = "top",
    fill_col = ratio_palette[length(ratio_palette)])
  
  tick_df <- data.frame(
    y = c(0, ratio_breaks),
    label = scales::label_number(accuracy = 0.1)(c(0, ratio_breaks)))
  
  ggplot() +
    geom_rect(
      data = bar_df,
      aes(xmin = bar_xmin, xmax = bar_xmax, ymin = ymin, ymax = ymax, fill = fill_col),
      colour = NA) +
    geom_polygon(
      data = tri_df,
      aes(x = x, y = y, group = group, fill = fill_col),
      colour = NA) +
    geom_path(
      data = data.frame(
        x = c(bar_xmin, bar_xmin, (bar_xmin + bar_xmax) / 2, bar_xmax, bar_xmax),
        y = c(ratio_min, ratio_max, ratio_max + top_triangle_height, ratio_max, ratio_min)),
      aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black") +
    geom_segment(
      aes(x = bar_xmin, xend = bar_xmax, y = ratio_min, yend = ratio_min),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black") +
    geom_segment(
      data = tick_df,
      aes(x = bar_xmax, xend = bar_xmax + 0.11, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3,
      colour = "black") +
    geom_text(
      data = tick_df,
      aes(x = bar_xmax + 0.17, y = y, label = label),
      hjust = 0,
      size = 3.8) +
    annotate(
      "text",
      x = bar_xmin,
      y = ratio_max + top_triangle_height + (0.05 * ratio_span),
      label = legend_title,
      hjust = 0,
      vjust = 0,
      fontface = "bold",
      size = 5.5) +
    scale_fill_identity() +
    coord_cartesian(
      xlim = c(0, 1.95),
      ylim = c(ratio_min - (0.04 * ratio_span), ratio_max + top_triangle_height + (0.16 * ratio_span)),
      clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(14, 16, 10, 8))}

# Reading one .nc file -----------------------------------------------------
nc <- open.nc(nc_file)
on.exit(close.nc(nc), add = TRUE)

lon <- var.get.nc(nc, "global_longitude0")
lat <- var.get.nc(nc, "global_latitude0")
rain <- var.get.nc(nc, "item5216_daily_mean")[, , 1]

if (length(dim(lon)) == 3) lon <- lon[, , 1]
if (length(dim(lat)) == 3) lat <- lat[, , 1]

nlon <- dim(rain)[1]
nlat <- dim(rain)[2]

# Applying the LSE mask ----------------------------------------------------
mask_matrix <- load_lse_mask_matrix(
  mask_file = lse_mask_file,
  nlon = nlon,
  nlat = nlat,
  transpose_mask = mask_transpose)

mask_is_land <- !is.na(mask_matrix)
rain[!mask_is_land] <- NaN

non_missing_after_mask <- sum(is.finite(rain))
if (non_missing_after_mask == 0) {
  stop("Masking removed all rainfall cells. Toggle mask_transpose and rerun.")}

message(sprintf(
  "Mask applied successfully. Finite rainfall cells after mask: %d of %d. Mask TRUE cells: %d.",
  non_missing_after_mask,
  length(rain),
  sum(mask_is_land, na.rm = TRUE)))

# Saving a masked NetCDF copy ----------------------------------------------
file.copy(nc_file, masked_nc_file, overwrite = TRUE)

nc_masked <- ncdf4::nc_open(masked_nc_file, write = TRUE)
precip_all <- ncdf4::ncvar_get(nc_masked, "item5216_daily_mean")

if (length(dim(precip_all)) == 3) {
  for (t in seq_len(dim(precip_all)[3])) {
    layer <- precip_all[, , t]
    layer[!mask_is_land] <- NaN
    precip_all[, , t] <- layer}
} else if (length(dim(precip_all)) == 4) {
  for (z in seq_len(dim(precip_all)[3])) {
    for (t in seq_len(dim(precip_all)[4])) {
      layer <- precip_all[, , z, t]
      layer[!mask_is_land] <- NaN
      precip_all[, , z, t] <- layer}}
} else {
  ncdf4::nc_close(nc_masked)
  stop("Unexpected dimensions for item5216_daily_mean. Expected 3D or 4D variable.")}

ncdf4::ncvar_put(nc_masked, "item5216_daily_mean", precip_all)
ncdf4::nc_close(nc_masked)

# Reading the probability-ratio grid ---------------------------------------
ratio_grid <- utils::read.csv(ratio_grid_file, stringsAsFactors = FALSE)
ratio_vars <- c(
  "probability_ratio_ge4_future_over_current",
  "probability_ratio_ge5_future_over_current",
  "probability_ratio_joint_top10_ge5_future_over_current",
  "probability_ratio_rx1day_top10_future_over_current",
  "probability_ratio_joint_top10_ge4_future_over_current")

missing_ratio_vars <- setdiff(c("lon_index", "lat_index", ratio_vars), names(ratio_grid))
if (length(missing_ratio_vars) > 0) {
  stop(sprintf("Missing required columns in ratio grid file: %s", paste(missing_ratio_vars, collapse = ", ")))}

# Building ratio layers and NZ intersections -------------------------------
ratio_layers <- list()
for (ratio_var in ratio_vars) {
  ratio_mat <- matrix_from_indexed_values(
    df = ratio_grid,
    value_col = ratio_var,
    nlon = nlon,
    nlat = nlat)
  
  ratio_mat[!mask_is_land] <- NA_real_
  
  ratio_poly <- build_cell_polygons(
    lon_mat = lon,
    lat_mat = lat,
    value_mat = ratio_mat,
    value_name = "ratio_value")
  
  keep_ids <- get_nz_intersecting_cell_ids(ratio_poly)
  
  ratio_layers[[ratio_var]] <- list(
    poly = ratio_poly,
    keep_ids = keep_ids)}

# Setting shared breaks across NZ cells ------------------------------------
intersection_values <- c()
for (ratio_var in ratio_vars) {
  layer <- ratio_layers[[ratio_var]]
  keep_rows <- layer$poly$id %in% layer$keep_ids
  intersection_values <- c(intersection_values, layer$poly$ratio_value[keep_rows])}
intersection_values <- intersection_values[is.finite(intersection_values)]

if (length(intersection_values) == 0) {
  stop("No finite probability-ratio values found for NZ-intersecting cells.")}

ratio_breaks <- c(1, 1.5, 2, 2.5, 3, 4, 5)

# Setting probability-ratio colours
ratio_palette <- c(
  "#D0D4DA", # 0-1
  "#D7E8FF", # 1-1.5
  "#BFD9FF", # 1.5-2
  "#7FB3FF", # 2-2.5
  "#3F8BE6", # 2.5-3
  "#0B4FAF", # 3-4
  "#08306B", # 4-5
  "#041F4A"  # >5
)

# Building the plots --------------------------------------------------------
p_top10 <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_rx1day_top10_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_rx1day_top10_future_over_current"]]$keep_ids,
  "(a) Years with extreme Rx1day",
  ratio_breaks,
  ratio_palette)

p_ge4 <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_ge4_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_ge4_future_over_current"]]$keep_ids,
  "(b) Years with ≥4 exceedances",
  ratio_breaks,
  ratio_palette)

p_joint <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_joint_top10_ge4_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_joint_top10_ge4_future_over_current"]]$keep_ids,
  "(c) Years with extreme Rx1day AND ≥4 exceedances",
  ratio_breaks,
  ratio_palette) +
  theme(
    plot.title = element_text(hjust = 0.5))

p_ge5 <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_ge5_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_ge5_future_over_current"]]$keep_ids,
  "(a) Years with ≥5 exceedances",
  ratio_breaks,
  ratio_palette) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = -18)))

p_ge5_joint <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_joint_top10_ge5_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_joint_top10_ge5_future_over_current"]]$keep_ids,
  "(b) Years with extreme Rx1day AND ≥5 exceedances",
  ratio_breaks,
  ratio_palette) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = -18)))

combined_design <- c(
  area(t = 1, l = 1, b = 1, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 2, l = 1, b = 2, r = 2),
  area(t = 1, l = 3, b = 2, r = 3))

p_ratio_legend <- make_triangle_colorbar_plot(ratio_breaks, ratio_palette)

p_combined <- (p_top10 + p_ge4 + p_joint + p_ratio_legend) +
  plot_layout(design = combined_design,
              widths = c(1, 1, 0.40), heights = c(1, 1))

p_ge5_with_legend <- p_ge5 + p_ge5_joint + p_ratio_legend +
  plot_layout(widths = c(1, 1, 0.40), heights = c(11,11,0.8))

print(p_combined)
print(p_ge5_with_legend)

# Saving outputs ------------------------------------------------------------
ggsave(filename = combined_ratio_output_png, plot = p_combined,
       width = 16, height = 11, dpi = 300)

ggsave(filename = ge5_ratio_output_png, plot = p_ge5_with_legend,
       width = 16, height = 11, dpi = 300)

# Running quick checks ------------------------------------------------------
cat("NZ-intersecting cell count (>=4):", length(ratio_layers[["probability_ratio_ge4_future_over_current"]]$keep_ids),"\n")
# 185 cells 

cat("NZ-intersecting cell count (>=5):", length(ratio_layers[["probability_ratio_ge5_future_over_current"]]$keep_ids),"\n")
# 185 cells

cat("NZ-intersecting cell count (>=5 + top 10%):", length(ratio_layers[["probability_ratio_joint_top10_ge5_future_over_current"]]$keep_ids),"\n")
# 185 cells

cat("NZ-intersecting cell count (top 10%):", length(ratio_layers[["probability_ratio_rx1day_top10_future_over_current"]]$keep_ids),"\n")
# 185 cells

cat("NZ-intersecting cell count (joint):", length(ratio_layers[["probability_ratio_joint_top10_ge4_future_over_current"]]$keep_ids),"\n")
# 185 cells

