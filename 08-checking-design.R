# -------------------------------------------------------------------------
# 04b-processing-weather@home-extreme-percentiles.R
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

calc_extreme_prob_blocks <- function(current_rx_array,
                                     future_rx_array,
                                     current_exceedance_array,
                                     future_exceedance_array,
                                     extreme_prob) {
  extreme_threshold_current <- apply(
    current_rx_array,
    c(1, 2),
    quantile,
    probs = extreme_prob,
    na.rm = TRUE,
    type = 7
  )
  
  current_prop_extreme <- calc_rx1day_top10_proportion(
    current_rx_array,
    extreme_threshold_current
  )
  future_prop_extreme <- calc_rx1day_top10_proportion(
    future_rx_array,
    extreme_threshold_current
  )
  
  current_prop_joint_ge4 <- calc_joint_top10_ge4_proportion(
    rx_array = current_rx_array,
    exceedance_array = current_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 4
  )
  
  future_prop_joint_ge4 <- calc_joint_top10_ge4_proportion(
    rx_array = future_rx_array,
    exceedance_array = future_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 4
  )
  
  current_prop_joint_ge5 <- calc_joint_top10_ge4_proportion(
    rx_array = current_rx_array,
    exceedance_array = current_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 5
  )
  
  future_prop_joint_ge5 <- calc_joint_top10_ge4_proportion(
    rx_array = future_rx_array,
    exceedance_array = future_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 5
  )
  
  list(
    threshold_current = extreme_threshold_current,
    current_prop_extreme = current_prop_extreme,
    future_prop_extreme = future_prop_extreme,
    probability_ratio_extreme = calc_probability_ratio(current_prop_extreme, future_prop_extreme),
    current_prop_joint_ge4 = current_prop_joint_ge4,
    future_prop_joint_ge4 = future_prop_joint_ge4,
    probability_ratio_joint_ge4 = calc_probability_ratio(current_prop_joint_ge4, future_prop_joint_ge4),
    current_prop_joint_ge5 = current_prop_joint_ge5,
    future_prop_joint_ge5 = future_prop_joint_ge5,
    probability_ratio_joint_ge5 = calc_probability_ratio(current_prop_joint_ge5, future_prop_joint_ge5)
  )
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

# Additional percentile-specific outputs for mapping scripts ----------------
build_extreme_grid_results <- function(grid_template, base_grid_results, extreme_blocks, suffix) {
  data.frame(
    lon_index = grid_template$lon_index,
    lat_index = grid_template$lat_index,
    longitude0 = base_grid_results$longitude0,
    latitude0 = base_grid_results$latitude0,
    global_longitude0 = base_grid_results$global_longitude0,
    global_latitude0 = base_grid_results$global_latitude0,
    rx1day_threshold_33_current = base_grid_results$rx1day_threshold_33_current,
    prop_years_ge4_current = base_grid_results$prop_years_ge4_current,
    prop_years_ge4_future = base_grid_results$prop_years_ge4_future,
    probability_ratio_ge4_future_over_current = base_grid_results$probability_ratio_ge4_future_over_current,
    prop_years_ge5_current = base_grid_results$prop_years_ge5_current,
    prop_years_ge5_future = base_grid_results$prop_years_ge5_future,
    probability_ratio_ge5_future_over_current = base_grid_results$probability_ratio_ge5_future_over_current,
    extreme_threshold_current = extreme_blocks$threshold_current[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_rx1day_extreme_current = extreme_blocks$current_prop_extreme[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_rx1day_extreme_future = extreme_blocks$future_prop_extreme[cbind(grid_template$lon_index, grid_template$lat_index)],
    probability_ratio_rx1day_extreme_future_over_current = extreme_blocks$probability_ratio_extreme[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge4_current = extreme_blocks$current_prop_joint_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge4_future = extreme_blocks$future_prop_joint_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
    probability_ratio_joint_extreme_ge4_future_over_current = extreme_blocks$probability_ratio_joint_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge5_current = extreme_blocks$current_prop_joint_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge5_future = extreme_blocks$future_prop_joint_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
    probability_ratio_joint_extreme_ge5_future_over_current = extreme_blocks$probability_ratio_joint_ge5[cbind(grid_template$lon_index, grid_template$lat_index)]
  ) |>
    dplyr::rename(
      !!paste0("rx1day_threshold_", suffix, "_current") := extreme_threshold_current,
      !!paste0("prop_years_rx1day_", suffix, "_current") := prop_years_rx1day_extreme_current,
      !!paste0("prop_years_rx1day_", suffix, "_future") := prop_years_rx1day_extreme_future,
      !!paste0("probability_ratio_rx1day_", suffix, "_future_over_current") := probability_ratio_rx1day_extreme_future_over_current,
      !!paste0("prop_years_joint_", suffix, "_ge4_current") := prop_years_joint_extreme_ge4_current,
      !!paste0("prop_years_joint_", suffix, "_ge4_future") := prop_years_joint_extreme_ge4_future,
      !!paste0("probability_ratio_joint_", suffix, "_ge4_future_over_current") := probability_ratio_joint_extreme_ge4_future_over_current,
      !!paste0("prop_years_joint_", suffix, "_ge5_current") := prop_years_joint_extreme_ge5_current,
      !!paste0("prop_years_joint_", suffix, "_ge5_future") := prop_years_joint_extreme_ge5_future,
      !!paste0("probability_ratio_joint_", suffix, "_ge5_future_over_current") := probability_ratio_joint_extreme_ge5_future_over_current
    )
}

extreme_blocks_p98_9 <- calc_extreme_prob_blocks(
  current_rx_array = current_rx_array,
  future_rx_array = future_rx_array,
  current_exceedance_array = current_exceedance_array,
  future_exceedance_array = future_exceedance_array,
  extreme_prob = 0.989
)

grid_results_p98_9 <- build_extreme_grid_results(
  grid_template = grid_template,
  base_grid_results = grid_results,
  extreme_blocks = extreme_blocks_p98_9,
  suffix = "p98_9"
)

write.csv(
  grid_results_p98_9,
  file = file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_p98_9_joint_probability_ratio_grid.csv"),
  row.names = FALSE
)

extreme_blocks_p96_6 <- calc_extreme_prob_blocks(
  current_rx_array = current_rx_array,
  future_rx_array = future_rx_array,
  current_exceedance_array = current_exceedance_array,
  future_exceedance_array = future_exceedance_array,
  extreme_prob = 0.966
)

grid_results_p96_6 <- build_extreme_grid_results(
  grid_template = grid_template,
  base_grid_results = grid_results,
  extreme_blocks = extreme_blocks_p96_6,
  suffix = "p96_6"
)

write.csv(
  grid_results_p96_6,
  file = file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_p96_6_joint_probability_ratio_grid.csv"),
  row.names = FALSE
)

max(probability_ratio_ge4, na.rm = TRUE) # Max ge4
mean(probability_ratio_ge4[is.finite(probability_ratio_ge4)], na.rm = TRUE) # Mean ge4
min(probability_ratio_ge4, na.rm = TRUE) # Min ge4

max(probability_ratio_ge5, na.rm = TRUE) # Max ge5
mean(probability_ratio_ge5[is.finite(probability_ratio_ge5)], na.rm = TRUE) # Mean ge5
min(probability_ratio_ge5, na.rm = TRUE) # Min ge5


max(probability_ratio_rx1day_top10, na.rm = TRUE) # Max top10 RX1day
mean(probability_ratio_rx1day_top10[is.finite(probability_ratio_rx1day_top10)], na.rm = TRUE) # Mean top10 RX1day
min(probability_ratio_rx1day_top10, na.rm = TRUE) # Min top10 RX1day

max(probability_ratio_joint_top10_ge4, na.rm = TRUE) # Max joint top10 + ge4
mean(probability_ratio_joint_top10_ge4[is.finite(probability_ratio_joint_top10_ge4)], na.rm = TRUE) # Mean joint top10 + ge4
min(probability_ratio_joint_top10_ge4, na.rm = TRUE) # Min joint top10 + ge4

max(probability_ratio_joint_top10_ge5, na.rm = TRUE) # Max joint top10 + ge5
mean(probability_ratio_joint_top10_ge5[is.finite(probability_ratio_joint_top10_ge5)], na.rm = TRUE) # Mean joint top10 + ge5
min(probability_ratio_joint_top10_ge5, na.rm = TRUE) # Min joint top10 + ge5


# -------------------------------------------------------------------------
# 04b-processing-weather@home-extreme-percentiles.R
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

calc_extreme_prob_blocks <- function(current_rx_array,
                                     future_rx_array,
                                     current_exceedance_array,
                                     future_exceedance_array,
                                     extreme_prob) {
  extreme_threshold_current <- apply(
    current_rx_array,
    c(1, 2),
    quantile,
    probs = extreme_prob,
    na.rm = TRUE,
    type = 7
  )
  
  current_prop_extreme <- calc_rx1day_top10_proportion(
    current_rx_array,
    extreme_threshold_current
  )
  future_prop_extreme <- calc_rx1day_top10_proportion(
    future_rx_array,
    extreme_threshold_current
  )
  
  current_prop_joint_ge4 <- calc_joint_top10_ge4_proportion(
    rx_array = current_rx_array,
    exceedance_array = current_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 4
  )
  
  future_prop_joint_ge4 <- calc_joint_top10_ge4_proportion(
    rx_array = future_rx_array,
    exceedance_array = future_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 4
  )
  
  current_prop_joint_ge5 <- calc_joint_top10_ge4_proportion(
    rx_array = current_rx_array,
    exceedance_array = current_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 5
  )
  
  future_prop_joint_ge5 <- calc_joint_top10_ge4_proportion(
    rx_array = future_rx_array,
    exceedance_array = future_exceedance_array,
    threshold_matrix = extreme_threshold_current,
    min_days = 5
  )
  
  list(
    threshold_current = extreme_threshold_current,
    current_prop_extreme = current_prop_extreme,
    future_prop_extreme = future_prop_extreme,
    probability_ratio_extreme = calc_probability_ratio(current_prop_extreme, future_prop_extreme),
    current_prop_joint_ge4 = current_prop_joint_ge4,
    future_prop_joint_ge4 = future_prop_joint_ge4,
    probability_ratio_joint_ge4 = calc_probability_ratio(current_prop_joint_ge4, future_prop_joint_ge4),
    current_prop_joint_ge5 = current_prop_joint_ge5,
    future_prop_joint_ge5 = future_prop_joint_ge5,
    probability_ratio_joint_ge5 = calc_probability_ratio(current_prop_joint_ge5, future_prop_joint_ge5)
  )
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

# Additional percentile-specific outputs for mapping scripts ----------------
build_extreme_grid_results <- function(grid_template, base_grid_results, extreme_blocks, suffix) {
  data.frame(
    lon_index = grid_template$lon_index,
    lat_index = grid_template$lat_index,
    longitude0 = base_grid_results$longitude0,
    latitude0 = base_grid_results$latitude0,
    global_longitude0 = base_grid_results$global_longitude0,
    global_latitude0 = base_grid_results$global_latitude0,
    rx1day_threshold_33_current = base_grid_results$rx1day_threshold_33_current,
    prop_years_ge4_current = base_grid_results$prop_years_ge4_current,
    prop_years_ge4_future = base_grid_results$prop_years_ge4_future,
    probability_ratio_ge4_future_over_current = base_grid_results$probability_ratio_ge4_future_over_current,
    prop_years_ge5_current = base_grid_results$prop_years_ge5_current,
    prop_years_ge5_future = base_grid_results$prop_years_ge5_future,
    probability_ratio_ge5_future_over_current = base_grid_results$probability_ratio_ge5_future_over_current,
    extreme_threshold_current = extreme_blocks$threshold_current[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_rx1day_extreme_current = extreme_blocks$current_prop_extreme[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_rx1day_extreme_future = extreme_blocks$future_prop_extreme[cbind(grid_template$lon_index, grid_template$lat_index)],
    probability_ratio_rx1day_extreme_future_over_current = extreme_blocks$probability_ratio_extreme[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge4_current = extreme_blocks$current_prop_joint_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge4_future = extreme_blocks$future_prop_joint_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
    probability_ratio_joint_extreme_ge4_future_over_current = extreme_blocks$probability_ratio_joint_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge5_current = extreme_blocks$current_prop_joint_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
    prop_years_joint_extreme_ge5_future = extreme_blocks$future_prop_joint_ge5[cbind(grid_template$lon_index, grid_template$lat_index)],
    probability_ratio_joint_extreme_ge5_future_over_current = extreme_blocks$probability_ratio_joint_ge5[cbind(grid_template$lon_index, grid_template$lat_index)]
  ) |>
    dplyr::rename(
      !!paste0("rx1day_threshold_", suffix, "_current") := extreme_threshold_current,
      !!paste0("prop_years_rx1day_", suffix, "_current") := prop_years_rx1day_extreme_current,
      !!paste0("prop_years_rx1day_", suffix, "_future") := prop_years_rx1day_extreme_future,
      !!paste0("probability_ratio_rx1day_", suffix, "_future_over_current") := probability_ratio_rx1day_extreme_future_over_current,
      !!paste0("prop_years_joint_", suffix, "_ge4_current") := prop_years_joint_extreme_ge4_current,
      !!paste0("prop_years_joint_", suffix, "_ge4_future") := prop_years_joint_extreme_ge4_future,
      !!paste0("probability_ratio_joint_", suffix, "_ge4_future_over_current") := probability_ratio_joint_extreme_ge4_future_over_current,
      !!paste0("prop_years_joint_", suffix, "_ge5_current") := prop_years_joint_extreme_ge5_current,
      !!paste0("prop_years_joint_", suffix, "_ge5_future") := prop_years_joint_extreme_ge5_future,
      !!paste0("probability_ratio_joint_", suffix, "_ge5_future_over_current") := probability_ratio_joint_extreme_ge5_future_over_current
    )
}

extreme_blocks_p98_9 <- calc_extreme_prob_blocks(
  current_rx_array = current_rx_array,
  future_rx_array = future_rx_array,
  current_exceedance_array = current_exceedance_array,
  future_exceedance_array = future_exceedance_array,
  extreme_prob = 0.989
)

grid_results_p98_9 <- build_extreme_grid_results(
  grid_template = grid_template,
  base_grid_results = grid_results,
  extreme_blocks = extreme_blocks_p98_9,
  suffix = "p98_9"
)

write.csv(
  grid_results_p98_9,
  file = file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_p98_9_joint_probability_ratio_grid.csv"),
  row.names = FALSE
)

extreme_blocks_p96_6 <- calc_extreme_prob_blocks(
  current_rx_array = current_rx_array,
  future_rx_array = future_rx_array,
  current_exceedance_array = current_exceedance_array,
  future_exceedance_array = future_exceedance_array,
  extreme_prob = 0.966
)

grid_results_p96_6 <- build_extreme_grid_results(
  grid_template = grid_template,
  base_grid_results = grid_results,
  extreme_blocks = extreme_blocks_p96_6,
  suffix = "p96_6"
)

write.csv(
  grid_results_p96_6,
  file = file.path(weatherathome_dir, "weather@home_exceedance_ge4_ge5_p96_6_joint_probability_ratio_grid.csv"),
  row.names = FALSE
)

max(probability_ratio_ge4, na.rm = TRUE) # Max ge4
mean(probability_ratio_ge4[is.finite(probability_ratio_ge4)], na.rm = TRUE) # Mean ge4
min(probability_ratio_ge4, na.rm = TRUE) # Min ge4

max(probability_ratio_ge5, na.rm = TRUE) # Max ge5
mean(probability_ratio_ge5[is.finite(probability_ratio_ge5)], na.rm = TRUE) # Mean ge5
min(probability_ratio_ge5, na.rm = TRUE) # Min ge5


max(probability_ratio_rx1day_top10, na.rm = TRUE) # Max top10 RX1day
mean(probability_ratio_rx1day_top10[is.finite(probability_ratio_rx1day_top10)], na.rm = TRUE) # Mean top10 RX1day
min(probability_ratio_rx1day_top10, na.rm = TRUE) # Min top10 RX1day

max(probability_ratio_joint_top10_ge4, na.rm = TRUE) # Max joint top10 + ge4
mean(probability_ratio_joint_top10_ge4[is.finite(probability_ratio_joint_top10_ge4)], na.rm = TRUE) # Mean joint top10 + ge4
min(probability_ratio_joint_top10_ge4, na.rm = TRUE) # Min joint top10 + ge4

max(probability_ratio_joint_top10_ge5, na.rm = TRUE) # Max joint top10 + ge5
mean(probability_ratio_joint_top10_ge5[is.finite(probability_ratio_joint_top10_ge5)], na.rm = TRUE) # Mean joint top10 + ge5
min(probability_ratio_joint_top10_ge5, na.rm = TRUE) # Min joint top10 + ge5

# -------------------------------------------------------------------------
# 05-mapping-weather@home-extreme-p98_9.R
# -------------------------------------------------------------------------
# Feb 2026
# Mapping weather@home data and saving NZ probability-ratio maps
# using an extreme Rx1day threshold at the 98.9th percentile.
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
ratio_grid_file_primary <- file.path(model_data_dir, "weather@home_exceedance_ge4_ge5_p98_9_joint_probability_ratio_grid.csv")
ratio_grid_file_fallback <- file.path(model_data_dir, "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv")

# Transposing the land mask
# If it looks weird and the wrong way around can can from True to False 
mask_transpose <- TRUE

masked_nc_file <- file.path(model_data_dir, "current_decade",
                            "item5216_daily_mean_a000_2006-07_2007-06-NZtrim-mm-masked.nc")

# Setting ratio map output names -------------------------------------------
combined_ratio_output_png <- file.path(model_data_dir,
                                       "weather@home_probability_ratio_ge4_p98_9_joint_combined_map.png")

ge5_ratio_output_png <- file.path(model_data_dir,
                                  "weather@home_probability_ratio_ge5_map.png")

# Cell highlight ---------------------------------------------------
matched_cell <- data.frame(lon_index = 30L, lat_index = 16L)
matched_land_mask_value <- 294.9544983

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
  
  p <- ggplot(poly_nz, aes(x = lon, y = lat, fill = ratio_bin)) +
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
      legend.position = "none")
  
  selected_id <- paste(matched_cell$lon_index, matched_cell$lat_index, sep = "_")
  selected_poly <- poly_df[poly_df$id == selected_id, ]
  if (nrow(selected_poly) > 0) {
    p <- p + geom_polygon(
      data = selected_poly,
      aes(x = lon, y = lat, group = id),
      inherit.aes = FALSE,
      fill = NA,
      colour = "#f03b20",
      linewidth = 0.7
    )
  } else {
    warning(sprintf("Selected cell id '%s' was not found in polygon data.", selected_id))
  }
  
  p}

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
      inherit.aes = FALSE,
      colour = NA) +
    geom_polygon(
      data = tri_df,
      aes(x = x, y = y, group = group, fill = fill_col),
      inherit.aes = FALSE,
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

selected_cell_id <- paste(matched_cell$lon_index, matched_cell$lat_index, sep = "_")
if (matched_cell$lon_index >= 1 && matched_cell$lon_index <= nlon &&
    matched_cell$lat_index >= 1 && matched_cell$lat_index <= nlat) {
  selected_mask_value <- mask_matrix[matched_cell$lon_index, matched_cell$lat_index]
  message(sprintf(
    "Selected map cell %s has land-mask value: %.10f",
    selected_cell_id,
    selected_mask_value
  ))
  
  matched_by_value <- which(
    mask_matrix == matched_land_mask_value,
    arr.ind = TRUE
  )
  if (nrow(matched_by_value) > 0) {
    message(sprintf(
      "Land-mask value %.7f found at %d cell(s): %s",
      matched_land_mask_value,
      nrow(matched_by_value),
      paste(apply(matched_by_value, 1, function(idx) sprintf("(%d,%d)", idx[1], idx[2])), collapse = ", ")
    ))
  } else {
    warning(sprintf(
      "Land-mask value %.7f was not found in the loaded mask.",
      matched_land_mask_value
    ))
  }
} else {
  warning(sprintf(
    "Selected map cell indices are out of range: lon_index=%d (max=%d), lat_index=%d (max=%d).",
    matched_cell$lon_index, nlon, matched_cell$lat_index, nlat
  ))
}

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
ratio_grid_file <- if (file.exists(ratio_grid_file_primary)) {
  ratio_grid_file_primary
} else {
  warning(sprintf(
    "Primary percentile-specific ratio grid not found (%s). Falling back to top10 ratio grid (%s).",
    ratio_grid_file_primary,
    ratio_grid_file_fallback
  ))
  ratio_grid_file_fallback
}

ratio_grid <- utils::read.csv(ratio_grid_file, stringsAsFactors = FALSE)
ratio_vars_preferred <- c(
  "probability_ratio_ge4_future_over_current",
  "probability_ratio_ge5_future_over_current",
  "probability_ratio_joint_p98_9_ge5_future_over_current",
  "probability_ratio_rx1day_p98_9_future_over_current",
  "probability_ratio_joint_p98_9_ge4_future_over_current")
ratio_vars_fallback <- c(
  "probability_ratio_ge4_future_over_current",
  "probability_ratio_ge5_future_over_current",
  "probability_ratio_joint_top10_ge5_future_over_current",
  "probability_ratio_rx1day_top10_future_over_current",
  "probability_ratio_joint_top10_ge4_future_over_current")

has_preferred <- all(c("lon_index", "lat_index", ratio_vars_preferred) %in% names(ratio_grid))
has_fallback <- all(c("lon_index", "lat_index", ratio_vars_fallback) %in% names(ratio_grid))

if (has_preferred) {
  ratio_vars <- ratio_vars_preferred
  extreme_rx_ratio_var <- "probability_ratio_rx1day_p98_9_future_over_current"
  extreme_joint_ge4_ratio_var <- "probability_ratio_joint_p98_9_ge4_future_over_current"
  extreme_joint_ge5_ratio_var <- "probability_ratio_joint_p98_9_ge5_future_over_current"
} else if (has_fallback) {
  ratio_vars <- ratio_vars_fallback
  extreme_rx_ratio_var <- "probability_ratio_rx1day_top10_future_over_current"
  extreme_joint_ge4_ratio_var <- "probability_ratio_joint_top10_ge4_future_over_current"
  extreme_joint_ge5_ratio_var <- "probability_ratio_joint_top10_ge5_future_over_current"
  warning(sprintf(
    "Percentile-specific ratio columns for %sth percentile not found. Using top10 ratio columns instead.",
    "98.9"
  ))
} else {
  missing_ratio_vars <- setdiff(
    c("lon_index", "lat_index", ratio_vars_preferred, ratio_vars_fallback),
    names(ratio_grid)
  )
  stop(sprintf("Missing required columns in ratio grid file: %s", paste(unique(missing_ratio_vars), collapse = ", ")))
}

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
  ratio_layers[[extreme_rx_ratio_var]]$poly,
  ratio_layers[[extreme_rx_ratio_var]]$keep_ids,
  "(a) Years with Extreme Rx1day (98.9th percentile)",
  ratio_breaks,
  ratio_palette)

p_ge4 <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_ge4_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_ge4_future_over_current"]]$keep_ids,
  "(b) Years with ≥4 Heavy days",
  ratio_breaks,
  ratio_palette)

p_joint <- make_nz_ratio_plot(
  ratio_layers[[extreme_joint_ge4_ratio_var]]$poly,
  ratio_layers[[extreme_joint_ge4_ratio_var]]$keep_ids,
  "(c) Years with Extreme Rx1day (98.9th percentile) AND ≥4 Heavy days",
  ratio_breaks,
  ratio_palette) +
  theme(
    plot.title = element_text(hjust = 0.5))

p_ge5 <- make_nz_ratio_plot(
  ratio_layers[["probability_ratio_ge5_future_over_current"]]$poly,
  ratio_layers[["probability_ratio_ge5_future_over_current"]]$keep_ids,
  "(a) Years with ≥5 Heavy days",
  ratio_breaks,
  ratio_palette) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = -18)))

p_ge5_joint <- make_nz_ratio_plot(
  ratio_layers[[extreme_joint_ge5_ratio_var]]$poly,
  ratio_layers[[extreme_joint_ge5_ratio_var]]$keep_ids,
  "(b) Years with Extreme Rx1day (98.9th percentile) AND ≥5 Heavy days",
  ratio_breaks,
  ratio_palette) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = -18)))

combined_design <- c(
  patchwork::area(t = 1, l = 1, b = 1, r = 1),
  patchwork::area(t = 1, l = 2, b = 1, r = 2),
  patchwork::area(t = 2, l = 1, b = 2, r = 2),
  patchwork::area(t = 1, l = 3, b = 2, r = 3))

p_ratio_legend <- make_triangle_colorbar_plot(ratio_breaks, ratio_palette)

p_combined <- (p_top10 + p_ge4 + p_joint + p_ratio_legend) +
  patchwork::plot_layout(
    design = combined_design,
    widths = c(1, 1, 0.40),
    heights = c(1, 1))

p_ge5_with_legend <- p_ge5 + p_ge5_joint + p_ratio_legend +
  patchwork::plot_layout(widths = c(1, 1, 0.40))

print(p_combined)
print(p_ge5_with_legend)

# Saving outputs ------------------------------------------------------------
ggsave(filename = combined_ratio_output_png, plot = p_combined,
       width = 16, height = 11, dpi = 300)

ggsave(filename = ge5_ratio_output_png, plot = p_ge5_with_legend,
       width = 16, height = 8.5, dpi = 300)

# Running quick checks ------------------------------------------------------
cat("NZ-intersecting cell count (>=4):", length(ratio_layers[["probability_ratio_ge4_future_over_current"]]$keep_ids),"\n")
# 185 cells 

cat("NZ-intersecting cell count (>=5):", length(ratio_layers[["probability_ratio_ge5_future_over_current"]]$keep_ids),"\n")
# 185 cells

cat("NZ-intersecting cell count (>=5 + extreme p98.9):", length(ratio_layers[[extreme_joint_ge5_ratio_var]]$keep_ids),"\n")
# 185 cells

cat("NZ-intersecting cell count (extreme p98.9):", length(ratio_layers[[extreme_rx_ratio_var]]$keep_ids),"\n")
# 185 cells

cat("NZ-intersecting cell count (joint):", length(ratio_layers[[extreme_joint_ge4_ratio_var]]$keep_ids),"\n")
# 185 cells
