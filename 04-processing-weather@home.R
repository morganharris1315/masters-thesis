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
current_exceedance_list <- lapply(
  current_day_files,
  count_exceedance_days,
  threshold_matrix = rx1day_threshold_33_current
)

current_exceedance_array <- simplify2array(current_exceedance_list)
dim(current_exceedance_array)
# 44 X 44 X 3226

# Exceedance day counts for all future years using the same current-day
# threshold baseline (44 x 44 x n_years)
future_exceedance_list <- lapply(
  future_projection_files,
  count_exceedance_days,
  threshold_matrix = rx1day_threshold_33_current
)

future_exceedance_array <- simplify2array(future_exceedance_list)
dim(future_exceedance_array)
# 44 X 44 X 2535

# Proportion of years with >= 4 exceedance days ----------------------------
# For each grid cell, calculate the proportion of years where exceedance
# day count is greater than or equal to 4.

calc_exceedance_proportion <- function(exceedance_array, min_days = 4) {
  apply(exceedance_array, c(1, 2), function(x) {
    mean(x >= min_days, na.rm = TRUE)
  })
}

current_prop_ge4 <- calc_exceedance_proportion(current_exceedance_array, min_days = 4)
future_prop_ge4 <- calc_exceedance_proportion(future_exceedance_array, min_days = 4)

dim(current_prop_ge4)
dim(future_prop_ge4)
# both are 44 X 44

# Probability ratio (future/current) ---------------------------------------
# Ratio > 1 means years with >=4 exceedance days are more frequent in the
# future simulation than in current climate.

probability_ratio_ge4 <- future_prop_ge4 / current_prop_ge4

# Handle divide-by-zero cases explicitly:
# - if current proportion is 0 and future > 0, ratio is Inf
# - if both are 0, set ratio to NA (undefined)
probability_ratio_ge4[current_prop_ge4 == 0 & future_prop_ge4 > 0] <- Inf
probability_ratio_ge4[current_prop_ge4 == 0 & future_prop_ge4 == 0] <- NA_real_

# Build a per-grid-cell output table with lon/lat and all required metrics --
# Includes:
# - rotated-model longitude0 and latitude0
# - current-day 33rd percentile RX1day threshold
# - proportion of years with >=4 exceedance days (current, future)
# - probability ratio (future/current)

nc_grid <- open.nc(current_day_files[1])
longitude0 <- var.get.nc(nc_grid, "longitude0")
latitude0 <- var.get.nc(nc_grid, "latitude0")
close.nc(nc_grid)

grid_template <- expand.grid(
  lon_index = seq_along(longitude0),
  lat_index = seq_along(latitude0)
)

grid_results <- data.frame(
  longitude0 = longitude0[grid_template$lon_index],
  latitude0 = latitude0[grid_template$lat_index],
  rx1day_threshold_33_current = rx1day_threshold_33_current[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge4_current = current_prop_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  prop_years_ge4_future = future_prop_ge4[cbind(grid_template$lon_index, grid_template$lat_index)],
  probability_ratio_future_over_current = probability_ratio_ge4[cbind(grid_template$lon_index, grid_template$lat_index)]
)

head(grid_results)

# Saving outputs -----------------------------------------------------------
# Save as both CSV (easy to inspect/share) and RDS (preserves types exactly).

output_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  grid_results,
  file = file.path(output_dir, "weatherathome_exceedance_ge4_probability_ratio_grid.csv"),
  row.names = FALSE
)

saveRDS(
  grid_results,
  file = file.path(output_dir, "weatherathome_exceedance_ge4_probability_ratio_grid.rds")
)
