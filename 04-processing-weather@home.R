# -------------------------------------------------------------------------
# 04-processing-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Loading in .nc files and understanding file structure
# Calculating RX1day, 33rd Percentile and exceedance days 
# Getting the proportion of greater than or equal to 4 exceedance days
# Creating probability ratio (Future projection/ Current Day) 
# -------------------------------------------------------------------------

# Loading in RNetCDF ------------------------------------------------------
#install.packages("RNetCDF")
library(RNetCDF)
library(sf)
library(maps)

# Exploring the structure of the data -------------------------------------
#Looking at a single file in the list of files in the folder current_day
current_day_files <- list.files ("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/current_decade",
                                 pattern = "\\.nc$",
                                 full.names = TRUE)

nc <- open.nc(current_day_files[1])

file.inq.nc(nc)
file.inq.nc(nc)['format']
print.nc(nc)

# 4 dimiensions: longitude0 (44), latitude0 (44), z0 (1), time1 (360)
# 8 vaiables, with rain in `item5216_daily_mean'

var.get.nc(nc,"longitude0")
var.get.nc(nc,"latitude0")

# Have rotated-model grid coordinates (not plain lat/lon)
# longitude0 runs ~194.52 to 213.44, latitude0 runs ~-4.84 to -23.76.

var.inq.nc(nc, "item5216_daily_mean")
precip <- var.get.nc(nc, "item5216_daily_mean")
dim(precip)

# `item5216_daily_mean` is precipitation with dimensions
# Dimensions are (longitude0, latitude0, z0, time1). But as z0 is length 1, R drops it when read and dimesions are...
# 44 x 44 x 360

#Random single grid cell
single_gird_cell <- precip[10, 20, ]

plot(single_gird_cell, type = "l")
max(single_gird_cell, na.rm = TRUE)

#Getting RX1day from precip varaible 
rx1day <- apply(precip, c(1,2), max, na.rm = TRUE)
dim(rx1day)

time <- var.get.nc(nc, "time1")
head(time)
tail(time)
# Time uses a 360_day calendar and units "days since 1840-12-01 00:00:00".
# Values ending in .5 are expected from daily means centred on each day.

close.nc(nc)

# Calculating RX1day -------------------------------------------------------
# For all years (files), all grid cells and both current day and future projections

# All current day files
current_day_files <- list.files ("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/current_decade",
                                 pattern = "\\.nc$",
                                 full.names = TRUE)

# Function to caculate RX1day
compute_rx1day_NetCDF <- function(file) {
  nc <- open.nc(file)
  pr <- var.get.nc(nc, "item5216_daily_mean")
  rx <- apply(pr, c(1,2), max, na.rm = TRUE)
  close.nc(nc)
  return(rx)
}

#Applying function to current day files
current_rx_list <- lapply(current_day_files, compute_rx1day_NetCDF)

length(current_rx_list)
# the length is 3226 e.g there is an RX1day value for every one of the 3226 years. 

current_rx_array <- simplify2array(current_rx_list)
dim(current_rx_array)
# 44 X 44 X 3226 

#Applying function to current day files
future_projection_files <- list.files("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/3k_warmer",
                                      pattern = "\\.nc$",
                                      full.names = TRUE)

future_rx_list <- lapply(future_projection_files, compute_rx1day_NetCDF)
length(future_rx_list)
# the length is 2535 e.g there is an RX1day value for every one of the 2535 years. 

future_rx_array <- simplify2array(future_rx_list)
dim(future_rx_array)
# 44 X 44 X 2535

# Calculating 33rd percentile threshold from current-day RX1day ------------
# Threshold is calculated for each grid cell using only current-day RX1day.
# Result is a 44 x 44 matrix where each cell has its own threshold value.

rx1day_threshold_33_current <- apply(
  current_rx_array,
  c(1, 2),
  quantile,
  probs = 0.33,
  na.rm = TRUE,
  type = 7
)

dim(rx1day_threshold_33_current)
# 44 X 44

# Calculating exceedance days above threshold for each year and grid cell --
# For each year file, count daily precipitation values above the cell-specific
# current-day RX1day 33rd percentile threshold.

count_exceedance_days <- function(file, threshold_matrix) {
  nc <- open.nc(file)
  pr <- var.get.nc(nc, "item5216_daily_mean")
  close.nc(nc)
  
  # Ensure precipitation is [lon, lat, time]. Some NetCDF reads retain a
  # singleton z-dimension (e.g., [lon, lat, 1, time]).
  pr_dim <- dim(pr)
  if (length(pr_dim) == 4 && pr_dim[3] == 1) {
    pr <- pr[, , 1, ]
    pr_dim <- dim(pr)
  }
  
  if (length(pr_dim) != 3) {
    stop("Expected precipitation array with dimensions [lon, lat, time].")
  }
  
  if (!all(pr_dim[1:2] == dim(threshold_matrix))) {
    stop("threshold_matrix dimensions do not match precipitation lon/lat grid.")
  }
  
  # Compare each day to the threshold matrix and count exceedance days.
  # sweep() applies the 2D threshold across the time dimension safely.
  exceedance_binary <- sweep(pr, c(1, 2), threshold_matrix, FUN = ">")
  exceedance_days <- apply(exceedance_binary, c(1, 2), sum, na.rm = TRUE)
  return(exceedance_days)
}

# Exceedance day counts for all current-day years (44 x 44 x n_years)
current_exceedance_list <- lapply(current_day_files,count_exceedance_days, threshold_matrix = rx1day_threshold_33_current)

current_exceedance_array <- simplify2array(current_exceedance_list)
dim(current_exceedance_array)
# 44 X 44 X 3226

# Exceedance day counts for all future years using the same current-day
# threshold baseline (44 x 44 x n_years)
future_exceedance_list <- lapply(future_projection_files, count_exceedance_days, threshold_matrix = rx1day_threshold_33_current)

#############################################
future_exceedance_array <- simplify2array(future_exceedance_list)
dim(future_exceedance_array)
# expect 44 X 44 X 2535

# Proportion of years with threshold exceedance days -----------------------
# Reusable helper for each threshold block below.

calc_exceedance_proportion <- function(exceedance_array, min_days) {
  apply(exceedance_array, c(1, 2), function(x) {
    mean(x >= min_days, na.rm = TRUE)
  })
}

calc_probability_ratio <- function(current_prop, future_prop) {
  probability_ratio <- future_prop / current_prop
  
  # Handle divide-by-zero cases explicitly:
  # - if current proportion is 0 and future > 0, ratio is Inf
  # - if both are 0, set ratio to NA (undefined)
  probability_ratio[current_prop == 0 & future_prop > 0] <- Inf
  probability_ratio[current_prop == 0 & future_prop == 0] <- NA_real_
  probability_ratio
}

# Probability-ratio block for >= 4 exceedance days -------------------------
# Ratio > 1 means years with >=4 exceedance days are more frequent in the
# future simulation than in current climate.

current_prop_ge4 <- calc_exceedance_proportion(current_exceedance_array, min_days = 4)
future_prop_ge4 <- calc_exceedance_proportion(future_exceedance_array, min_days = 4)
probability_ratio_ge4 <- calc_probability_ratio(current_prop_ge4, future_prop_ge4)

dim(current_prop_ge4)
dim(future_prop_ge4)
# both should be 44 X 44

# Probability-ratio block for >= 5 exceedance days -------------------------
# This is intentionally a separate run block from >=4 to execute independently.

current_prop_ge5 <- calc_exceedance_proportion(current_exceedance_array, min_days = 5)
future_prop_ge5 <- calc_exceedance_proportion(future_exceedance_array, min_days = 5)
probability_ratio_ge5 <- calc_probability_ratio(current_prop_ge5, future_prop_ge5)

dim(current_prop_ge5)
dim(future_prop_ge5)
# both should be 44 X 44

# Probability-ratio block for top 10% RX1day years -------------------------
# Threshold is the cell-specific 90th percentile from current-day RX1day,
# and this same threshold is used for both current-day and future projections.

rx1day_threshold_90_current <- apply(
  current_rx_array,
  c(1, 2),
  quantile,
  probs = 0.90,
  na.rm = TRUE,
  type = 7
)

rx1day_threshold_90_future <- apply(
  future_rx_array,
  c(1, 2),
  quantile,
  probs = 0.90,
  na.rm = TRUE,
  type = 7
)

rx1day_threshold_90_pct_change_future_over_current <-
  100 * (rx1day_threshold_90_future - rx1day_threshold_90_current) / rx1day_threshold_90_current
rx1day_threshold_90_pct_change_future_over_current[!is.finite(rx1day_threshold_90_pct_change_future_over_current)] <- NA_real_

calc_rx1day_top10_proportion <- function(rx_array, threshold_matrix) {
  exceed_top10 <- sweep(rx_array, c(1, 2), threshold_matrix, FUN = ">=")
  apply(exceed_top10, c(1, 2), mean, na.rm = TRUE)
}

current_prop_rx1day_top10 <- calc_rx1day_top10_proportion(
  current_rx_array,
  rx1day_threshold_90_current
)
future_prop_rx1day_top10 <- calc_rx1day_top10_proportion(
  future_rx_array,
  rx1day_threshold_90_current
)
probability_ratio_rx1day_top10 <- calc_probability_ratio(
  current_prop_rx1day_top10,
  future_prop_rx1day_top10
)

dim(current_prop_rx1day_top10)
dim(future_prop_rx1day_top10)
# both should be 44 X 44

# Probability-ratio block for joint event: top 10% RX1day and >=4 exceedances
# Uses:
# - top 10% threshold from current-day RX1day
# - >=4 exceedance days from the exceedance arrays
# Joint proportion is the share of years where both conditions occur together.

calc_joint_top10_ge4_proportion <- function(rx_array, exceedance_array, threshold_matrix, min_days = 4) {
  rx_top10 <- sweep(rx_array, c(1, 2), threshold_matrix, FUN = ">=")
  ge_min_days <- exceedance_array >= min_days
  joint_condition <- rx_top10 & ge_min_days
  apply(joint_condition, c(1, 2), mean, na.rm = TRUE)
}

current_prop_joint_top10_ge4 <- calc_joint_top10_ge4_proportion(
  rx_array = current_rx_array,
  exceedance_array = current_exceedance_array,
  threshold_matrix = rx1day_threshold_90_current,
  min_days = 4
)

future_prop_joint_top10_ge4 <- calc_joint_top10_ge4_proportion(
  rx_array = future_rx_array,
  exceedance_array = future_exceedance_array,
  threshold_matrix = rx1day_threshold_90_current,
  min_days = 4
)

probability_ratio_joint_top10_ge4 <- calc_probability_ratio(
  current_prop_joint_top10_ge4,
  future_prop_joint_top10_ge4
)

dim(current_prop_joint_top10_ge4)
dim(future_prop_joint_top10_ge4)
# both should be 44 X 44

# Probability-ratio block for joint event: top 10% RX1day and >=5 exceedances
# Uses the same top 10% threshold and joint-event logic, but with min_days = 5.

current_prop_joint_top10_ge5 <- calc_joint_top10_ge4_proportion(
  rx_array = current_rx_array,
  exceedance_array = current_exceedance_array,
  threshold_matrix = rx1day_threshold_90_current,
  min_days = 5
)

future_prop_joint_top10_ge5 <- calc_joint_top10_ge4_proportion(
  rx_array = future_rx_array,
  exceedance_array = future_exceedance_array,
  threshold_matrix = rx1day_threshold_90_current,
  min_days = 5
)

probability_ratio_joint_top10_ge5 <- calc_probability_ratio(
  current_prop_joint_top10_ge5,
  future_prop_joint_top10_ge5
)

dim(current_prop_joint_top10_ge5)
dim(future_prop_joint_top10_ge5)
# both should be 44 X 44

# Build a per-grid-cell output table with lon/lat and all required metrics --
# Includes:
# - rotated-model longitude0 and latitude0
# - current-day 33rd and 90th percentile RX1day thresholds
# - proportion of years with >=4 and >=5 exceedance days (current, future)
# - proportion of years with top 10% RX1day and joint top10%+>=4 event
# - probability ratios (future/current) for all metrics

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

grid_results <- data.frame(
  lon_index = grid_template$lon_index,
  lat_index = grid_template$lat_index,
  longitude0 = longitude0[grid_template$lon_index],
  latitude0 = latitude0[grid_template$lat_index],
  global_longitude0 = global_longitude0[cbind(grid_template$lon_index, grid_template$lat_index)],
  global_latitude0 = global_latitude0[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_33_current = rx1day_threshold_33_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_90_current = rx1day_threshold_90_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_90_future = rx1day_threshold_90_future[cbind(grid_template$lon_index, grid_template$lat_index)],
  rx1day_threshold_90_pct_change_future_over_current = rx1day_threshold_90_pct_change_future_over_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge4_current = current_prop_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge4_future = future_prop_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_ge4_future_over_current = probability_ratio_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge5_current = current_prop_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge5_future = future_prop_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_ge5_future_over_current = probability_ratio_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_rx1day_top10_current = current_prop_rx1day_top10[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_rx1day_top10_future = future_prop_rx1day_top10[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_rx1day_top10_future_over_current = probability_ratio_rx1day_top10[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_top10_ge4_current = current_prop_joint_top10_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_top10_ge4_future = future_prop_joint_top10_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_joint_top10_ge4_future_over_current = probability_ratio_joint_top10_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_top10_ge5_current = current_prop_joint_top10_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_joint_top10_ge5_future = future_prop_joint_top10_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_joint_top10_ge5_future_over_current = probability_ratio_joint_top10_ge5[cbind(grid_template$lon_index, grid_template$lat_index)]
)

head(grid_results)

# Saving outputs -----------------------------------------------------------

weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"

write.csv(
  grid_results,
  file = file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv"),
  row.names = FALSE
)

# NZ land-cell filtering for summary statistics ----------------------------
lse_mask_file <- file.path(weatherathome_dir, "Land-Sea Mask for Weather@home Data.csv")
mask_transpose <- TRUE

load_lse_mask_matrix <- function(mask_file, nlon, nlat, transpose_mask = FALSE) {
  mask_tbl <- utils::read.csv(mask_file, stringsAsFactors = FALSE)
  vals <- suppressWarnings(as.numeric(unlist(mask_tbl, use.names = FALSE)))
  
  base_candidates <- list(
    matrix(vals, nrow = nlon, ncol = nlat, byrow = TRUE),
    matrix(vals, nrow = nlon, ncol = nlat, byrow = FALSE),
    matrix(vals, nrow = nlat, ncol = nlon, byrow = TRUE),
    matrix(vals, nrow = nlat, ncol = nlon, byrow = FALSE)
  )
  
  dims_ok <- vapply(base_candidates, function(m) {
    nr <- nrow(m)
    nc <- ncol(m)
    (nr == nlon && nc == nlat) || (nr == nlat && nc == nlon)
  }, logical(1))
  
  if (!any(dims_ok)) {
    stop(sprintf("Could not coerce LSE mask to %d x %d matrix.", nlon, nlat))}
  
  mask <- base_candidates[[which(dims_ok)[1]]]
  if (isTRUE(transpose_mask)) {
    mask <- t(mask)
  }
  
  if (nrow(mask) != nlon || ncol(mask) != nlat) {
    stop("Mask dimensions do not match rainfall grid after transpose setting.")
  } 
  mask
}
  
calc_finite_stats <- function(x) {
  x_finite <- x[is.finite(x)]
  if (length(x_finite) == 0) {
    return(c(min = NA_real_, mean = NA_real_, max = NA_real_))
  }
  c(min = min(x_finite), mean = mean(x_finite), max = max(x_finite))
}

nlon <- nrow(rx1day_threshold_33_current)
nlat <- ncol(rx1day_threshold_33_current)
mask_matrix <- load_lse_mask_matrix(
  mask_file = lse_mask_file,
  nlon = nlon,
  nlat = nlat,
  transpose_mask = mask_transpose
)
mask_is_land <- !is.na(mask_matrix)

grid_results$nz_land_cell <- mask_is_land[cbind(grid_results$lon_index, grid_results$lat_index)]
grid_results_nz <- grid_results[which(grid_results$nz_land_cell), , drop = FALSE]

write.csv(
  grid_results_nz,
  file = file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid_nz_land.csv"),
  row.names = FALSE
)

nz_probability_ratio_ge4_stats <- calc_finite_stats(grid_results_nz$probability_ratio_ge4_future_over_current)
nz_probability_ratio_ge5_stats <- calc_finite_stats(grid_results_nz$probability_ratio_ge5_future_over_current)
nz_probability_ratio_rx1day_top10_stats <- calc_finite_stats(grid_results_nz$probability_ratio_rx1day_top10_future_over_current)
nz_probability_ratio_joint_top10_ge4_stats <- calc_finite_stats(grid_results_nz$probability_ratio_joint_top10_ge4_future_over_current)
nz_probability_ratio_joint_top10_ge5_stats <- calc_finite_stats(grid_results_nz$probability_ratio_joint_top10_ge5_future_over_current)

mean_rx1day_threshold_90_pct_change_future_over_current_nz <- mean(
  grid_results_nz$rx1day_threshold_90_pct_change_future_over_current,
  na.rm = TRUE
)

nz_probability_ratio_ge4_stats
nz_probability_ratio_ge5_stats
nz_probability_ratio_rx1day_top10_stats
nz_probability_ratio_joint_top10_ge4_stats
nz_probability_ratio_joint_top10_ge5_stats
mean_rx1day_threshold_90_pct_change_future_over_current_nz

saveRDS(current_rx_array, "current_rx_array.rds")
saveRDS(future_rx_array, "future_rx_array.rds")

saveRDS(current_exceedance_array, "current_exceedance_array.rds")
saveRDS(future_exceedance_array, "future_exceedance_array.rds")
  
